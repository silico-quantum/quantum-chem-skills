#!/usr/bin/env python3
"""Run a version-piloted Multiwfn response file with fail-closed gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import select
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path


REJECT_PATTERNS = {
    "fatal": re.compile(r"\bfatal(?:\s+error)?\b", re.IGNORECASE),
    "segmentation_fault": re.compile(r"segmentation fault", re.IGNORECASE),
    "fortran_runtime": re.compile(r"fortran runtime error", re.IGNORECASE),
    "unexpected_eof": re.compile(
        r"unexpected (?:end of file|eof)|read past end", re.IGNORECASE
    ),
    "invalid_selection": re.compile(
        r"invalid (?:menu )?(?:selection|choice|input)", re.IGNORECASE
    ),
    "read_error": re.compile(r"(?:input|read) error", re.IGNORECASE),
}
SHA256_PATTERN = re.compile(r"[0-9a-fA-F]{64}\Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_sha256(value: str) -> str:
    if not SHA256_PATTERN.fullmatch(value):
        raise argparse.ArgumentTypeError("expected exactly 64 hexadecimal characters")
    return value.lower()


def identity_record(metadata: os.stat_result) -> dict[str, int]:
    return {
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "size_bytes": metadata.st_size,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
        "mode": stat.S_IMODE(metadata.st_mode),
    }


def read_stable_regular_file(path: Path) -> tuple[bytes, os.stat_result]:
    """Read one unchanged regular file without following a symbolic link."""
    if path.is_symlink():
        raise ValueError(f"output must not be a symbolic link: {path}")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"output is not a regular file: {path}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)

    current = path.lstat()
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    current_identity = (
        current.st_dev,
        current.st_ino,
        current.st_size,
        current.st_mtime_ns,
        current.st_ctime_ns,
    )
    if identity_before != identity_after or identity_after != current_identity:
        raise ValueError(f"output changed while it was being validated: {path}")
    return b"".join(chunks), after


def resolve_executable(value: str) -> str:
    candidate = Path(value).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        resolved = candidate.resolve()
        if not resolved.is_file() or not os.access(resolved, os.X_OK):
            raise FileNotFoundError(f"Multiwfn executable is not runnable: {resolved}")
        return str(resolved)
    located = shutil.which(value)
    if located is None:
        raise FileNotFoundError(f"Multiwfn executable was not found: {value}")
    return str(Path(located).resolve())


def sentinel_is_closed(descriptor: int) -> bool:
    readable, _, _ = select.select([descriptor], [], [], 0)
    if not readable:
        return False
    if os.read(descriptor, 1):
        raise RuntimeError("process-tree sentinel received unexpected data")
    return True


def probe_process_group(process_group_id: int) -> str:
    """Return present, absent, or permission_denied for a POSIX process group."""
    try:
        os.killpg(process_group_id, 0)
    except ProcessLookupError:
        return "absent"
    except PermissionError:
        return "permission_denied"
    return "present"


def process_tree_has_survivors(process_group_id: int, sentinel: int) -> bool:
    # The sentinel detects ordinary descendants, while killpg(0) also detects a
    # same-group descendant that deliberately closes inherited descriptors.
    return (
        not sentinel_is_closed(sentinel)
        or probe_process_group(process_group_id) != "absent"
    )


def wait_for_process_tree_exit(
    process_group_id: int, sentinel: int, timeout: float
) -> bool:
    deadline = time.monotonic() + timeout
    while process_tree_has_survivors(process_group_id, sentinel):
        if time.monotonic() >= deadline:
            return False
        time.sleep(0.02)
    return True


def wait_for_inherited_sentinel_close(sentinel: int, timeout: float) -> bool:
    """Confirm every descendant retaining the launch sentinel has exited."""
    deadline = time.monotonic() + timeout
    while not sentinel_is_closed(sentinel):
        if time.monotonic() >= deadline:
            return False
        time.sleep(0.02)
    return True


def signal_process_group(process_group_id: int, signal_number: int) -> str:
    try:
        os.killpg(process_group_id, signal_number)
    except ProcessLookupError:
        return "absent"
    except PermissionError:
        return "permission_denied"
    return "sent"


def terminate_process_group(process: subprocess.Popen[bytes], sentinel: int) -> bool:
    """Terminate the whole group even after its session leader has exited."""
    process_group_id = process.pid
    term_result = signal_process_group(process_group_id, signal.SIGTERM)
    if term_result == "permission_denied":
        return False
    if process.poll() is None:
        try:
            process.wait(timeout=1)
        except subprocess.TimeoutExpired:
            pass
    if wait_for_process_tree_exit(process_group_id, sentinel, 1):
        return True
    kill_result = signal_process_group(process_group_id, signal.SIGKILL)
    if kill_result == "permission_denied":
        return False
    if process.poll() is None:
        try:
            process.wait(timeout=1)
        except subprocess.TimeoutExpired:
            pass
    if kill_result == "absent":
        return sentinel_is_closed(sentinel)
    # A successful group-wide SIGKILL prevents further userspace writes. The
    # sentinel closing proves inherited descendants have exited; any remaining
    # process-group entry is therefore an unreaped zombie.
    return wait_for_inherited_sentinel_close(sentinel, 4)


def path_in_workdir(value: Path, workdir: Path) -> Path:
    candidate = value.expanduser() if value.is_absolute() else workdir / value
    return candidate.resolve()


def require_path_in_workdir(value: Path, workdir: Path, label: str) -> Path:
    path = path_in_workdir(value, workdir)
    if not path.is_relative_to(workdir):
        raise ValueError(f"{label} path escapes the workdir: {path}")
    return path


def write_atomic_json(path: Path, payload: dict[str, object]) -> None:
    """Publish a new JSON file atomically without replacing an existing path."""
    if path.exists():
        raise FileExistsError(f"refusing to overwrite report: {path}")
    if not path.parent.is_dir():
        raise FileNotFoundError(f"report parent does not exist: {path.parent}")
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def write_execution_snapshots(
    workdir: Path,
    input_file: Path,
    input_bytes: bytes,
    response_bytes: bytes,
) -> tuple[Path, Path]:
    snapshot_dir = Path(tempfile.mkdtemp(prefix=".multiwfn-snapshot-", dir=workdir))
    input_snapshot = snapshot_dir / f"input{input_file.suffix}"
    responses_snapshot = snapshot_dir / "responses.in"
    input_snapshot.write_bytes(input_bytes)
    responses_snapshot.write_bytes(response_bytes)
    input_snapshot.chmod(0o444)
    responses_snapshot.chmod(0o444)
    return input_snapshot, responses_snapshot


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--responses", type=Path, required=True)
    parser.add_argument("--executable", default="Multiwfn")
    parser.add_argument(
        "--expected-executable-sha256",
        type=parse_sha256,
        required=True,
        help="exact SHA-256 recorded for the executable used by the accepted pilot",
    )
    parser.add_argument("--timeout-seconds", type=float, required=True)
    parser.add_argument("--load-marker", required=True)
    parser.add_argument(
        "--version-marker",
        required=True,
        help="unique exact full banner/version line captured by the accepted pilot",
    )
    parser.add_argument("--task-marker", required=True)
    parser.add_argument("--exit-marker", required=True)
    parser.add_argument(
        "--final-response",
        required=True,
        help="last non-empty response captured for clean program exit",
    )
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        action="append",
        required=True,
        help="expected fresh output; repeat for every primary or auxiliary file",
    )
    return parser


def execute(args: argparse.Namespace) -> tuple[dict[str, object], int]:
    workdir = args.workdir.expanduser().resolve()
    input_file = args.input_file.expanduser().resolve()
    responses = args.responses.expanduser().resolve()
    log = path_in_workdir(args.log, workdir)
    report = path_in_workdir(args.report, workdir)
    outputs = [path_in_workdir(path, workdir) for path in args.output]

    base_report: dict[str, object] = {
        "status": "rejected",
        "process_status": "rejected",
        "scientific_status": "not_evaluated",
        "workdir": str(workdir),
        "input": str(input_file),
        "responses": str(responses),
        "log": str(log),
        "expected_outputs": [str(path) for path in outputs],
        "timeout_seconds": args.timeout_seconds,
        "expected_executable_sha256": args.expected_executable_sha256,
        "markers": {
            "version": args.version_marker,
            "load": args.load_marker,
            "task": args.task_marker,
            "exit": args.exit_marker,
        },
        "checks": {},
        "errors": [],
    }
    errors = base_report["errors"]
    checks = base_report["checks"]
    assert isinstance(errors, list)
    assert isinstance(checks, dict)

    input_bytes: bytes | None = None
    response_bytes: bytes | None = None
    if not workdir.is_dir():
        errors.append("workdir does not exist")
    if input_file.is_file():
        input_bytes = input_file.read_bytes()
    if not input_bytes:
        errors.append("input file is absent or empty")
    if responses.is_file():
        response_bytes = responses.read_bytes()
    if not response_bytes:
        errors.append("response file is absent or empty")
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        errors.append("timeout must be finite and positive")
    contained_paths = [
        ("log", log),
        ("report", report),
        *(("output", path) for path in outputs),
    ]
    if len({path for _, path in contained_paths}) != len(contained_paths):
        errors.append("log, report, and every analysis output must be distinct paths")
    for name, path in contained_paths:
        if not path.is_relative_to(workdir):
            errors.append(f"{name} path escapes the workdir: {path}")
    for marker_name in (
        "version_marker",
        "load_marker",
        "task_marker",
        "exit_marker",
    ):
        if not getattr(args, marker_name):
            errors.append(f"{marker_name} is empty")
    markers = [
        args.version_marker,
        args.load_marker,
        args.task_marker,
        args.exit_marker,
    ]
    if any(
        left in right or right in left
        for index, left in enumerate(markers)
        for right in markers[index + 1 :]
    ):
        errors.append(
            "version, load, task, and exit markers must be distinct and nonoverlapping"
        )
    process_group_supported = os.name == "posix" and hasattr(os, "killpg")
    checks["posix_process_group_isolation_available"] = process_group_supported
    if not process_group_supported:
        errors.append("POSIX process-group isolation is required for batch execution")
    if report.exists() or report.is_symlink():
        raise FileExistsError(f"refusing to overwrite report: {report}")
    if not report.parent.is_dir():
        raise FileNotFoundError(f"report parent does not exist: {report.parent}")
    for path in [log, *outputs]:
        if path.exists() or path.is_symlink():
            errors.append(f"refusing stale or existing path: {path}")
        if not path.parent.is_dir():
            errors.append(f"output parent does not exist: {path.parent}")

    executable: str | None = None
    executable_before_identity: dict[str, int] | None = None
    try:
        executable = resolve_executable(args.executable)
        executable_before_bytes, executable_before_stat = read_stable_regular_file(
            Path(executable)
        )
        executable_before_sha256 = sha256_bytes(executable_before_bytes)
        executable_before_identity = identity_record(executable_before_stat)
        executable_matches_pilot = (
            executable_before_sha256 == args.expected_executable_sha256
        )
        base_report["executable"] = {
            "requested": args.executable,
            "resolved_path": executable,
            "expected_sha256": args.expected_executable_sha256,
            "sha256_before": executable_before_sha256,
            "identity_before": executable_before_identity,
            "sha256_after": None,
            "identity_after": None,
        }
    except (OSError, ValueError) as exc:
        executable_matches_pilot = False
        base_report["executable"] = {
            "requested": args.executable,
            "resolved_path": executable,
            "expected_sha256": args.expected_executable_sha256,
            "sha256_before": None,
            "identity_before": None,
            "sha256_after": None,
            "identity_after": None,
            "error": f"{type(exc).__name__}: {exc}",
        }
    checks["executable_matches_pilot"] = executable_matches_pilot
    if not executable_matches_pilot:
        errors.append("Multiwfn executable does not match the accepted pilot SHA-256")

    if response_bytes:
        try:
            response_text = response_bytes.decode("utf-8")
        except UnicodeDecodeError:
            checks["final_response_matches_pilot"] = False
            errors.append("response file is not valid UTF-8")
        else:
            response_lines = [
                line.strip() for line in response_text.splitlines() if line.strip()
            ]
            checks["final_response_matches_pilot"] = bool(
                response_lines and response_lines[-1] == args.final_response
            )
            if not checks["final_response_matches_pilot"]:
                errors.append(
                    "response file does not end with the declared clean-exit response"
                )

    if errors:
        return base_report, 1

    assert input_bytes is not None
    assert response_bytes is not None
    input_snapshot, responses_snapshot = write_execution_snapshots(
        workdir, input_file, input_bytes, response_bytes
    )
    input_digest = hashlib.sha256(input_bytes).hexdigest()
    responses_digest = hashlib.sha256(response_bytes).hexdigest()

    assert executable is not None
    command = [executable, str(input_snapshot)]
    base_report["command"] = command
    base_report["input_snapshot"] = str(input_snapshot)
    base_report["responses_snapshot"] = str(responses_snapshot)
    base_report["input_sha256"] = input_digest
    base_report["responses_sha256"] = responses_digest
    base_report["input_snapshot_sha256"] = sha256_file(input_snapshot)
    base_report["responses_snapshot_sha256"] = sha256_file(responses_snapshot)

    timed_out = False
    sentinel_read, sentinel_write = os.pipe()
    os.set_inheritable(sentinel_write, True)
    try:
        with responses_snapshot.open("rb") as response_handle, log.open("xb") as log_handle:
            try:
                process = subprocess.Popen(
                    command,
                    cwd=workdir,
                    stdin=response_handle,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    pass_fds=(sentinel_write,),
                )
            finally:
                os.close(sentinel_write)
            try:
                return_code = process.wait(timeout=args.timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                process_group_cleanup_confirmed = terminate_process_group(
                    process, sentinel_read
                )
                return_code = process.returncode

        process_group_had_no_survivors = not process_tree_has_survivors(
            process.pid, sentinel_read
        )
        if process_group_had_no_survivors:
            process_group_cleanup_confirmed = True
        else:
            process_group_cleanup_confirmed = terminate_process_group(
                process, sentinel_read
            )
            errors.append(
                "Multiwfn left a live descendant in its isolated process group"
            )
    finally:
        os.close(sentinel_read)
    checks["process_group_had_no_survivors"] = process_group_had_no_survivors
    checks["process_group_cleanup_confirmed"] = process_group_cleanup_confirmed
    base_report["process_group"] = {
        "id": process.pid,
        "had_no_survivors_after_leader_exit": process_group_had_no_survivors,
        "cleanup_confirmed_before_artifact_capture": process_group_cleanup_confirmed,
    }
    if not process_group_cleanup_confirmed:
        errors.append("could not confirm that the Multiwfn process group exited")
        base_report["return_code"] = return_code
        checks["completed_before_timeout"] = not timed_out
        checks["return_code_zero"] = return_code == 0
        return base_report, 1

    base_report["return_code"] = return_code
    checks["completed_before_timeout"] = not timed_out
    checks["return_code_zero"] = return_code == 0
    if timed_out:
        errors.append("Multiwfn exceeded the wall-clock timeout and was terminated")
    if return_code != 0:
        errors.append(f"Multiwfn exited with status {return_code}")

    executable_record = base_report["executable"]
    assert isinstance(executable_record, dict)
    try:
        executable_after_bytes, executable_after_stat = read_stable_regular_file(
            Path(executable)
        )
        executable_after_sha256 = sha256_bytes(executable_after_bytes)
        executable_after_identity = identity_record(executable_after_stat)
        executable_record["sha256_after"] = executable_after_sha256
        executable_record["identity_after"] = executable_after_identity
        checks["executable_unchanged_after_process"] = (
            executable_after_sha256 == args.expected_executable_sha256
            and executable_after_identity == executable_before_identity
        )
    except (OSError, ValueError) as exc:
        executable_record["post_run_error"] = f"{type(exc).__name__}: {exc}"
        checks["executable_unchanged_after_process"] = False
    if not checks["executable_unchanged_after_process"]:
        errors.append("Multiwfn executable changed during execution")

    input_snapshot_after = (
        sha256_file(input_snapshot) if input_snapshot.is_file() else None
    )
    responses_snapshot_after = (
        sha256_file(responses_snapshot) if responses_snapshot.is_file() else None
    )
    base_report["input_snapshot_sha256_after_process"] = input_snapshot_after
    base_report["responses_snapshot_sha256_after_process"] = (
        responses_snapshot_after
    )
    checks["input_snapshot_unchanged"] = input_snapshot_after == input_digest
    checks["responses_snapshot_unchanged"] = (
        responses_snapshot_after == responses_digest
    )
    if not checks["input_snapshot_unchanged"]:
        errors.append("execution input snapshot changed during the child process")
    if not checks["responses_snapshot_unchanged"]:
        errors.append("execution response snapshot changed during the child process")

    log_text = log.read_text(encoding="utf-8", errors="replace")
    log_lines = [line.strip() for line in log_text.splitlines()]
    marker_positions: list[int] = []
    for name, marker in (
        ("version_marker", args.version_marker),
        ("load_marker", args.load_marker),
        ("task_marker", args.task_marker),
        ("clean_exit_marker", args.exit_marker),
    ):
        positions = [index for index, line in enumerate(log_lines) if line == marker]
        unique = len(positions) == 1
        checks[name] = unique
        if not unique:
            errors.append(
                f"expected exactly one full-line {name} {marker!r}; found {len(positions)}"
            )
        marker_positions.append(positions[0] if unique else -1)
    markers_ordered = (
        all(position >= 0 for position in marker_positions)
        and marker_positions == sorted(marker_positions)
        and len(set(marker_positions)) == 4
    )
    checks["markers_unique_and_ordered"] = markers_ordered
    if not markers_ordered:
        errors.append(
            "version, load, task, and clean-exit markers are not unique and ordered"
        )

    rejected = [name for name, pattern in REJECT_PATTERNS.items() if pattern.search(log_text)]
    checks["rejection_diagnostics_absent"] = not rejected
    if rejected:
        errors.append("rejected transcript diagnostics: " + ", ".join(rejected))

    output_records = []
    for path in outputs:
        try:
            output_bytes, output_stat = read_stable_regular_file(path)
            valid = bool(output_bytes)
            output_error = None if valid else "output is empty"
        except (OSError, ValueError) as exc:
            output_bytes = b""
            output_stat = None
            valid = False
            output_error = str(exc)
        output_records.append(
            {
                "path": str(path),
                "non_empty": valid,
                "regular_file": output_stat is not None,
                "size_bytes": len(output_bytes) if output_stat is not None else None,
                "sha256": hashlib.sha256(output_bytes).hexdigest() if valid else None,
            }
        )
        if not valid:
            errors.append(
                f"expected output is absent, empty, or invalid: {path}: {output_error}"
            )
    base_report["outputs"] = output_records
    checks["all_expected_outputs_non_empty"] = all(
        record["non_empty"] for record in output_records
    )
    base_report["log_sha256"] = sha256_file(log)

    accepted = not errors and all(bool(value) for value in checks.values())
    base_report["status"] = "process_validated" if accepted else "rejected"
    base_report["process_status"] = "accepted" if accepted else "rejected"
    base_report["scientific_status"] = (
        "pending_independent_validation" if accepted else "not_evaluated"
    )
    return base_report, 0 if accepted else 1


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workdir = args.workdir.expanduser().resolve()
    try:
        report_path = require_path_in_workdir(args.report, workdir, "report")
    except Exception as exc:
        print(f"REJECTED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    try:
        report, exit_code = execute(args)
    except Exception as exc:
        report = {
            "status": "rejected",
            "process_status": "rejected",
            "scientific_status": "not_evaluated",
            "error": f"{type(exc).__name__}: {exc}",
        }
        exit_code = 2
    write_atomic_json(report_path, report)
    print(f"{str(report.get('status', 'rejected')).upper()}: {report_path}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
