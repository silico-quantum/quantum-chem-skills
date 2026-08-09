#!/usr/bin/env python3
"""Run one xTB MD command without a shell and write an immutable receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import select
import shutil
import signal
import stat
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def resolve_inside(run_dir: Path, value: Path) -> Path:
    path = value if value.is_absolute() else run_dir / value
    resolved = path.resolve()
    if resolved != run_dir and run_dir not in resolved.parents:
        raise ValueError(f"path escapes run directory: {value}")
    return resolved


def resolve_input_source(run_dir: Path, value: Path) -> Path:
    candidate = value if value.is_absolute() else run_dir / value
    parent = candidate.parent.resolve()
    if parent != run_dir and run_dir not in parent.parents:
        raise ValueError(f"input path escapes run directory: {value}")
    return parent / candidate.name


def read_stable_regular_file(path: Path) -> tuple[bytes, os.stat_result]:
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
            raise ValueError(f"input file is absent, empty, or not regular: {path}")
        blocks: list[bytes] = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, name) != getattr(after, name) for name in stable_fields):
        raise ValueError(f"input changed while it was being captured: {path}")
    data = b"".join(blocks)
    if len(data) != after.st_size:
        raise ValueError(f"input size changed while it was being captured: {path}")
    return data, after


def identity_record(path: Path, data: bytes, metadata: os.stat_result) -> dict[str, object]:
    return {
        "path": str(path),
        "sha256": sha256_bytes(data),
        "size_bytes": len(data),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
    }


def write_private_snapshot(path: Path, data: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o400)
    finally:
        os.close(descriptor)


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


def terminate_process_group(process: subprocess.Popen[str], sentinel: int) -> bool:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--log", type=Path, default=Path("md.log"))
    parser.add_argument(
        "--receipt", type=Path, default=Path("run_receipt.json")
    )
    parser.add_argument(
        "--trajectory", type=Path, default=Path("xtb.trj")
    )
    parser.add_argument(
        "--success-marker", type=Path, default=Path("xtbmdok")
    )
    parser.add_argument(
        "--input-file",
        type=Path,
        action="append",
        required=True,
        help="run-dir-relative input to hash before launch; repeat for every input",
    )
    parser.add_argument("--timeout-seconds", type=float)
    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="exact xTB argv after --, for example: -- xtb input.xyz --gfnff --md ...",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    command = list(args.command)
    if command and command[0] == "--":
        command.pop(0)
    if not command:
        raise ValueError("an exact command is required after --")
    if "--md" not in command:
        raise ValueError("the command must explicitly contain --md")
    if args.timeout_seconds is not None and (
        not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0
    ):
        raise ValueError("timeout-seconds must be positive")

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"run directory does not exist: {run_dir}")
    log_path = resolve_inside(run_dir, args.log)
    receipt_path = resolve_inside(run_dir, args.receipt)
    trajectory_path = resolve_inside(run_dir, args.trajectory)
    marker_path = resolve_inside(run_dir, args.success_marker)
    publication_paths = (log_path, receipt_path, trajectory_path, marker_path)
    if len(set(publication_paths)) != len(publication_paths):
        raise ValueError("log, receipt, trajectory, and success marker must be distinct")
    for output in publication_paths:
        if output.exists():
            raise FileExistsError(f"refusing to overwrite output: {output}")
        if not output.parent.is_dir():
            raise FileNotFoundError(f"output parent does not exist: {output.parent}")
    if os.name != "posix" or not hasattr(os, "killpg"):
        raise RuntimeError(
            "POSIX process-group isolation is required for xTB execution"
        )

    input_paths = [resolve_input_source(run_dir, value) for value in args.input_file]
    if len(set(input_paths)) != len(input_paths):
        raise ValueError("input-file values must be unique")
    if set(input_paths) & set(publication_paths):
        raise ValueError("input files must be distinct from generated artifacts")
    captured_inputs: list[tuple[Path, bytes, os.stat_result]] = []
    for path in input_paths:
        try:
            captured_inputs.append((path, *read_stable_regular_file(path)))
        except OSError as exc:
            raise ValueError(f"cannot securely capture input file {path}: {exc}") from exc

    snapshot_dir = run_dir / ".xtb-input-snapshots"
    snapshot_dir.mkdir(mode=0o700)
    os.chmod(snapshot_dir, 0o700)
    input_records: list[dict[str, object]] = []
    for index, (source_path, source_bytes, source_stat) in enumerate(captured_inputs):
        snapshot_path = snapshot_dir / f"{index:03d}-{source_path.name}"
        write_private_snapshot(snapshot_path, source_bytes)
        snapshot_token = str(snapshot_path.relative_to(run_dir))
        input_records.append(
            {
                "source": identity_record(source_path, source_bytes, source_stat),
                "snapshot": {
                    "path": str(snapshot_path),
                    "argv_token": snapshot_token,
                    "sha256": sha256_bytes(source_bytes),
                    "size_bytes": len(source_bytes),
                    "mode": "0o400",
                },
            }
        )

    rewritten_command = list(command)
    for source_path, record in zip(input_paths, input_records, strict=True):
        matching_indices: list[int] = []
        for index, token in enumerate(rewritten_command[1:], start=1):
            try:
                if resolve_input_source(run_dir, Path(token)) == source_path:
                    matching_indices.append(index)
            except ValueError:
                continue
        if len(matching_indices) != 1:
            raise ValueError(
                f"command must reference input file exactly once: {source_path}; "
                f"found {len(matching_indices)}"
            )
        rewritten_command[matching_indices[0]] = record["snapshot"]["argv_token"]

    requested_executable = Path(rewritten_command[0]).expanduser()
    if requested_executable.is_absolute() or requested_executable.parent != Path("."):
        candidate = (
            requested_executable
            if requested_executable.is_absolute()
            else run_dir / requested_executable
        ).resolve()
        resolved_executable = (
            str(candidate)
            if candidate.is_file() and os.access(candidate, os.X_OK)
            else None
        )
    else:
        located = shutil.which(rewritten_command[0])
        resolved_executable = str(Path(located).resolve()) if located else None
    exact_argv = [resolved_executable or rewritten_command[0], *rewritten_command[1:]]
    executable_path = Path(resolved_executable) if resolved_executable else None
    executable_sha256 = (
        sha256_file(executable_path)
        if executable_path is not None and executable_path.is_file()
        else None
    )
    executable_size = (
        executable_path.stat().st_size
        if executable_path is not None and executable_path.is_file()
        else None
    )
    executable_mtime_ns = (
        executable_path.stat().st_mtime_ns
        if executable_path is not None and executable_path.is_file()
        else None
    )
    started = datetime.now(timezone.utc)
    started_time_ns = time.time_ns()
    start_clock = time.monotonic()
    status = "launch_failed"
    returncode: int | None = None
    caught_signal: int | None = None
    launch_error: str | None = None
    process_group_id: int | None = None
    process_group_had_no_survivors: bool | None = None
    process_group_cleanup_confirmed: bool | None = None

    sentinel_read, sentinel_write = os.pipe()
    os.set_inheritable(sentinel_write, True)
    try:
        with log_path.open("x", encoding="utf-8") as log_handle:
            try:
                try:
                    process = subprocess.Popen(
                        exact_argv,
                        cwd=run_dir,
                        stdin=subprocess.DEVNULL,
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                        start_new_session=True,
                        text=True,
                        pass_fds=(sentinel_write,),
                    )
                finally:
                    os.close(sentinel_write)
                process_group_id = process.pid
                try:
                    returncode = process.wait(timeout=args.timeout_seconds)
                    process_group_had_no_survivors = not process_tree_has_survivors(
                        process.pid, sentinel_read
                    )
                    if process_group_had_no_survivors:
                        process_group_cleanup_confirmed = True
                        status = "completed"
                    else:
                        process_group_cleanup_confirmed = terminate_process_group(
                            process, sentinel_read
                        )
                        status = (
                            "rejected_descendants"
                            if process_group_cleanup_confirmed
                            else "cleanup_failed"
                        )
                except subprocess.TimeoutExpired:
                    process_group_had_no_survivors = False
                    process_group_cleanup_confirmed = terminate_process_group(
                        process, sentinel_read
                    )
                    returncode = process.returncode
                    status = (
                        "timed_out"
                        if process_group_cleanup_confirmed
                        else "cleanup_failed"
                    )
                except KeyboardInterrupt:
                    process_group_had_no_survivors = False
                    process_group_cleanup_confirmed = terminate_process_group(
                        process, sentinel_read
                    )
                    returncode = process.returncode
                    status = (
                        "interrupted"
                        if process_group_cleanup_confirmed
                        else "cleanup_failed"
                    )
            except OSError as exc:
                launch_error = f"{type(exc).__name__}: {exc}"
    finally:
        os.close(sentinel_read)

    if process_group_cleanup_confirmed is False:
        raise RuntimeError(
            "could not confirm process-group exit; artifact capture is unsafe"
        )

    if returncode is not None and returncode < 0:
        caught_signal = -returncode
    finished = datetime.now(timezone.utc)

    for record in input_records:
        source_path = Path(record["source"]["path"])
        source_after: dict[str, object]
        unchanged = False
        try:
            current_bytes, current_stat = read_stable_regular_file(source_path)
            source_after = identity_record(source_path, current_bytes, current_stat)
            unchanged = (
                all(
                    source_after.get(name) == record["source"].get(name)
                    for name in (
                        "path",
                        "sha256",
                        "size_bytes",
                        "device",
                        "inode",
                        "mtime_ns",
                        "ctime_ns",
                    )
                )
            )
        except (OSError, ValueError) as exc:
            source_after = {"path": str(source_path), "error": f"{type(exc).__name__}: {exc}"}
        record["source_after_run"] = source_after
        record["source_unchanged_after_run"] = unchanged

    def artifact_record(path: Path) -> dict[str, object]:
        regular = path.is_file() and not path.is_symlink()
        size = path.stat().st_size if regular else None
        non_empty = regular and isinstance(size, int) and size > 0
        return {
            "path": str(path),
            "regular_file": regular,
            "size_bytes": size,
            "mtime_ns": path.stat().st_mtime_ns if regular else None,
            "sha256": sha256_file(path) if non_empty else None,
        }

    receipt = {
        "schema_version": 3,
        "status": status,
        "working_directory": str(run_dir),
        "argv": exact_argv,
        "requested_argv": command,
        "returncode": returncode,
        "signal": caught_signal,
        "started_at_utc": started.isoformat(),
        "started_time_ns": started_time_ns,
        "finished_at_utc": finished.isoformat(),
        "elapsed_seconds": time.monotonic() - start_clock,
        "timeout_seconds": args.timeout_seconds,
        "launch_error": launch_error,
        "process_group": {
            "id": process_group_id,
            "had_no_survivors": process_group_had_no_survivors,
            "cleanup_confirmed": process_group_cleanup_confirmed,
            "cleanup_timing": "before input and artifact capture",
        },
        "input_files": input_records,
        "executable": {
            "requested": command[0],
            "resolved_path": resolved_executable,
            "sha256": executable_sha256,
            "size_bytes": executable_size,
            "mtime_ns": executable_mtime_ns,
            "hash_timing": "captured before process launch",
        },
        "log": {
            "path": str(log_path),
            "sha256": sha256_file(log_path),
            "size_bytes": log_path.stat().st_size,
            "mtime_ns": log_path.stat().st_mtime_ns,
        },
        "generated_artifacts": {
            "trajectory": artifact_record(trajectory_path),
            "success_marker": artifact_record(marker_path),
        },
    }
    with receipt_path.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"{status.upper()}: returncode={returncode!r}; receipt={receipt_path}")
    if status == "completed" and returncode == 0:
        return 0
    if isinstance(returncode, int) and 0 < returncode <= 125:
        return returncode
    if caught_signal is not None:
        return min(255, 128 + caught_signal)
    return 124 if status == "timed_out" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
