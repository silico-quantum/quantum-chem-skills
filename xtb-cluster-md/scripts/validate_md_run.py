#!/usr/bin/env python3
"""Validate an xTB MD receipt, log, marker, input, and full trajectory."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import stat
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO


ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$|^X$")
REJECT_PATTERNS = {
    "md_unstable": re.compile(r"\bMD is unstable\b", re.IGNORECASE),
    "thermostat_problem": re.compile(r"\bthermostating problem\b", re.IGNORECASE),
    "fatal_error": re.compile(r"\bfatal error\b", re.IGNORECASE),
    "runtime_exception": re.compile(r"\bruntime exception occurred\b", re.IGNORECASE),
    "scc_not_converged": re.compile(r"\bSCC\b.*\bnot converg", re.IGNORECASE),
}


@dataclass
class TrajectoryStats:
    frame_count: int = 0
    labels: list[str] | None = None
    last_comment: str = ""
    global_minimum_inter_distance_angstrom: float | None = None
    global_minimum_pair: dict[str, object] | None = None
    maximum_extent_angstrom: float = 0.0
    distinct_coordinate_frames: int = 0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_stable_regular_file(path: Path) -> tuple[bytes, os.stat_result]:
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
            raise ValueError(f"bound file is empty or not regular: {path}")
        blocks: list[bytes] = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    identity_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, field) != getattr(after, field) for field in identity_fields):
        raise ValueError(f"bound file changed while it was being validated: {path}")
    data = b"".join(blocks)
    if len(data) != after.st_size:
        raise ValueError(f"bound file size changed while it was being validated: {path}")
    return data, after


def resolve_run_path(run_dir: Path, value: str | Path) -> Path:
    candidate = Path(value)
    path = candidate if candidate.is_absolute() else run_dir / candidate
    resolved = path.resolve()
    if resolved != run_dir and run_dir not in resolved.parents:
        raise ValueError(f"artifact path escapes run directory: {value}")
    return resolved


def resolve_run_leaf_no_follow(run_dir: Path, value: str | Path) -> Path:
    candidate = Path(value)
    path = candidate if candidate.is_absolute() else run_dir / candidate
    parent = path.parent.resolve()
    if parent != run_dir and run_dir not in parent.parents:
        raise ValueError(f"artifact path escapes run directory: {value}")
    return parent / path.name


def unique_matches(text: str, pattern: re.Pattern[str], label: str) -> list[str]:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise ValueError(f"log must contain exactly one {label}; found {len(matches)}")
    first = matches[0]
    return list(first) if isinstance(first, tuple) else [first]


def parse_unique_log_value(text: str, label: str) -> float:
    pattern = re.compile(
        rf"^\s*{re.escape(label)}\s*:\s*([-+0-9.eEdD]+)", re.MULTILINE
    )
    value = unique_matches(text, pattern, label)[0]
    parsed = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(parsed):
        raise ValueError(f"logged {label} is not finite")
    return parsed


def parse_unique_log_integer(text: str, label: str) -> int:
    value = parse_unique_log_value(text, label)
    if not value.is_integer():
        raise ValueError(f"logged {label} is not an integer")
    return int(value)


def exact_option_integer(argv: list[str], option: str) -> int | None:
    indices = [index for index, token in enumerate(argv) if token == option]
    if len(indices) != 1 or indices[0] + 1 >= len(argv):
        return None
    try:
        return int(argv[indices[0] + 1])
    except ValueError:
        return None


def read_required_line(handle: TextIO, description: str) -> str:
    line = handle.readline()
    if line == "":
        raise ValueError(f"trajectory ended while reading {description}")
    return line


def read_xyz_frame(
    handle: TextIO, frame_number: int, expected_atoms: int
) -> tuple[str, list[str], list[tuple[float, float, float]]] | None:
    atom_count_line = handle.readline()
    if atom_count_line == "":
        return None
    if not atom_count_line.strip():
        raise ValueError("blank lines between XYZ frames are not accepted")
    try:
        atom_count = int(atom_count_line.strip())
    except ValueError as exc:
        raise ValueError(f"frame {frame_number} atom count is not an integer") from exc
    if atom_count != expected_atoms:
        raise ValueError(
            f"frame {frame_number} has {atom_count} atoms; expected {expected_atoms}"
        )
    comment = read_required_line(handle, f"frame {frame_number} comment").rstrip("\n")
    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for atom_index in range(1, atom_count + 1):
        line = read_required_line(handle, f"frame {frame_number} atom {atom_index}")
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"frame {frame_number} atom {atom_index} is malformed")
        if not ELEMENT_RE.fullmatch(fields[0]):
            raise ValueError(
                f"frame {frame_number} atom {atom_index} has invalid element {fields[0]!r}"
            )
        try:
            xyz = tuple(float(value) for value in fields[1:4])
        except ValueError as exc:
            raise ValueError(
                f"frame {frame_number} atom {atom_index} has non-numeric coordinates"
            ) from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(
                f"frame {frame_number} atom {atom_index} has non-finite coordinates"
            )
        labels.append(fields[0])
        coordinates.append(xyz)
    return comment, labels, coordinates


def minimum_inter_distance(
    coordinates: list[tuple[float, float, float]],
    molecules: int,
    atoms_per_molecule: int,
) -> tuple[float, dict[str, object]]:
    minimum = math.inf
    closest: dict[str, object] = {}
    for left_molecule in range(molecules):
        left_start = left_molecule * atoms_per_molecule
        for right_molecule in range(left_molecule + 1, molecules):
            right_start = right_molecule * atoms_per_molecule
            for left_atom in range(atoms_per_molecule):
                for right_atom in range(atoms_per_molecule):
                    distance = math.dist(
                        coordinates[left_start + left_atom],
                        coordinates[right_start + right_atom],
                    )
                    if distance < minimum:
                        minimum = distance
                        closest = {
                            "molecules_zero_based": [left_molecule, right_molecule],
                            "atoms_within_molecule_zero_based": [left_atom, right_atom],
                        }
    return minimum, closest


def system_extent(coordinates: list[tuple[float, float, float]]) -> float:
    spans = [
        max(point[axis] for point in coordinates)
        - min(point[axis] for point in coordinates)
        for axis in range(3)
    ]
    return math.sqrt(sum(span * span for span in spans))


def validate_single_xyz(path: Path, expected_atoms: int) -> tuple[list[str], str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("input XYZ is absent or empty")
    with path.open(encoding="utf-8") as handle:
        frame = read_xyz_frame(handle, 1, expected_atoms)
        if frame is None:
            raise ValueError("input XYZ contains no frame")
        if handle.read().strip():
            raise ValueError("input XYZ must contain exactly one frame")
    comment, labels, _ = frame
    return labels, comment


def validate_trajectory(
    path: Path,
    molecules: int,
    atoms_per_molecule: int,
    minimum_distance: float,
    maximum_extent: float,
) -> TrajectoryStats:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("trajectory is absent or empty")
    expected_atoms = molecules * atoms_per_molecule
    stats = TrajectoryStats()
    reference_labels: list[str] | None = None
    coordinate_hashes: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        while True:
            frame = read_xyz_frame(handle, stats.frame_count + 1, expected_atoms)
            if frame is None:
                break
            comment, labels, coordinates = frame
            if reference_labels is None:
                reference_labels = labels
                reference_block = labels[:atoms_per_molecule]
                for molecule_index in range(molecules):
                    start = molecule_index * atoms_per_molecule
                    if labels[start : start + atoms_per_molecule] != reference_block:
                        raise ValueError(
                            "trajectory molecule blocks do not share the monomer atom labels/order"
                        )
            elif labels != reference_labels:
                raise ValueError(
                    f"atom labels or ordering changed in frame {stats.frame_count + 1}"
                )

            minimum, pair = minimum_inter_distance(
                coordinates, molecules, atoms_per_molecule
            )
            if minimum < minimum_distance:
                raise ValueError(
                    f"frame {stats.frame_count + 1} minimum intermolecular distance "
                    f"{minimum:.8g} is below {minimum_distance:.8g} Angstrom"
                )
            if (
                stats.global_minimum_inter_distance_angstrom is None
                or minimum < stats.global_minimum_inter_distance_angstrom
            ):
                stats.global_minimum_inter_distance_angstrom = minimum
                stats.global_minimum_pair = {
                    "frame_one_based": stats.frame_count + 1,
                    **pair,
                    "elements": [
                        labels[
                            pair["molecules_zero_based"][0] * atoms_per_molecule
                            + pair["atoms_within_molecule_zero_based"][0]
                        ],
                        labels[
                            pair["molecules_zero_based"][1] * atoms_per_molecule
                            + pair["atoms_within_molecule_zero_based"][1]
                        ],
                    ],
                }
            extent = system_extent(coordinates)
            if extent > maximum_extent:
                raise ValueError(
                    f"frame {stats.frame_count + 1} extent {extent:.8g} exceeds "
                    f"{maximum_extent:.8g} Angstrom"
                )
            stats.maximum_extent_angstrom = max(stats.maximum_extent_angstrom, extent)
            digest = hashlib.sha256()
            for point in coordinates:
                for value in point:
                    digest.update(value.hex().encode("ascii"))
                    digest.update(b"\0")
            coordinate_hashes.add(digest.hexdigest())
            stats.last_comment = comment
            stats.frame_count += 1

    if stats.frame_count == 0 or reference_labels is None:
        raise ValueError("trajectory contains no complete frame")
    if stats.frame_count < 2:
        raise ValueError("positive-time MD trajectory must contain at least two frames")
    if len(coordinate_hashes) < 2:
        raise ValueError("all trajectory coordinate frames are identical")
    stats.labels = reference_labels
    stats.distinct_coordinate_frames = len(coordinate_hashes)
    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--receipt", default="run_receipt.json", help="path relative to run-dir")
    parser.add_argument("--log", default="md.log", help="path relative to run-dir")
    parser.add_argument("--trajectory", default="xtb.trj", help="path relative to run-dir")
    parser.add_argument("--success-marker", default="xtbmdok", help="path relative to run-dir")
    parser.add_argument("--input-xyz", required=True, help="path relative to run-dir")
    parser.add_argument("--md-input", required=True, help="path relative to run-dir")
    parser.add_argument(
        "--expected-method", choices=("gfnff", "gfn0", "gfn1", "gfn2"), required=True
    )
    parser.add_argument("--expected-charge", type=int, required=True)
    parser.add_argument("--expected-uhf", type=int, required=True)
    parser.add_argument("--expected-xtb-version", required=True)
    parser.add_argument("--expected-molecules", type=int, required=True)
    parser.add_argument("--atoms-per-molecule", type=int, required=True)
    parser.add_argument("--requested-time-ps", type=float, required=True)
    parser.add_argument("--requested-step-fs", type=float, required=True)
    parser.add_argument("--requested-dump-fs", type=float, required=True)
    parser.add_argument("--requested-temperature-k", type=float, required=True)
    parser.add_argument("--requested-scc-accuracy", type=float, required=True)
    parser.add_argument("--requested-hydrogen-mass", type=float, required=True)
    parser.add_argument("--requested-shake-bonds", type=int, required=True)
    parser.add_argument("--expect-thermostat", choices=("on", "off"), required=True)
    parser.add_argument(
        "--minimum-intermolecular-distance", type=float, required=True
    )
    parser.add_argument("--maximum-system-extent", type=float, required=True)
    parser.add_argument("--report", type=Path, required=True, help="path relative to run-dir")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    positive_values = (
        args.expected_molecules,
        args.atoms_per_molecule,
        args.requested_time_ps,
        args.requested_step_fs,
        args.requested_dump_fs,
        args.requested_temperature_k,
        args.requested_scc_accuracy,
        args.requested_hydrogen_mass,
        args.minimum_intermolecular_distance,
        args.maximum_system_extent,
    )
    if (
        any(
            not math.isfinite(float(value)) or value <= 0
            for value in positive_values
        )
        or args.requested_shake_bonds < 0
        or args.expected_uhf < 0
    ):
        raise ValueError(
            "counts and physical thresholds must be positive; shake bonds and UHF cannot be negative"
        )
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+){1,3}(?:[-+][A-Za-z0-9._-]+)?", args.expected_xtb_version):
        raise ValueError("expected xTB version must be an explicit version identifier")

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"run directory does not exist: {run_dir}")
    report_path = resolve_run_path(run_dir, args.report)
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite report: {report_path}")
    if not report_path.parent.is_dir():
        raise FileNotFoundError(f"report parent does not exist: {report_path.parent}")

    receipt_path = resolve_run_path(run_dir, args.receipt)
    log_path = resolve_run_path(run_dir, args.log)
    trajectory_path = resolve_run_path(run_dir, args.trajectory)
    input_path = resolve_run_leaf_no_follow(run_dir, args.input_xyz)
    md_input_path = resolve_run_leaf_no_follow(run_dir, args.md_input)
    marker_path = resolve_run_path(run_dir, args.success_marker)
    bound_paths = (
        receipt_path,
        log_path,
        trajectory_path,
        input_path,
        md_input_path,
        marker_path,
        report_path,
    )
    if len(set(bound_paths)) != len(bound_paths):
        raise ValueError("receipt, log, trajectory, inputs, marker, and report must be distinct")
    checks: dict[str, bool] = {}
    errors: list[str] = []

    if not log_path.is_file() or log_path.stat().st_size == 0:
        log_text = ""
        errors.append("log is absent or empty")
    else:
        log_text = log_path.read_text(encoding="utf-8", errors="replace")

    receipt: dict[str, object] = {}
    if not receipt_path.is_file() or receipt_path.stat().st_size == 0:
        errors.append("process receipt is absent or empty")
    else:
        try:
            loaded = json.loads(receipt_path.read_text(encoding="utf-8"))
            if not isinstance(loaded, dict):
                raise ValueError("receipt root is not an object")
            receipt = loaded
        except (json.JSONDecodeError, ValueError) as exc:
            errors.append(f"process receipt is invalid: {exc}")

    checks["receipt_completed"] = receipt.get("status") == "completed"
    if not checks["receipt_completed"]:
        errors.append("process receipt status is not completed")
    checks["process_exit_zero"] = receipt.get("returncode") == 0
    if not checks["process_exit_zero"]:
        errors.append(f"xTB process return code is {receipt.get('returncode')!r}, not zero")
    checks["receipt_working_directory"] = receipt.get("working_directory") == str(run_dir)
    if not checks["receipt_working_directory"]:
        errors.append("process receipt working directory does not match run-dir")
    checks["receipt_schema"] = receipt.get("schema_version") == 3
    if not checks["receipt_schema"]:
        errors.append("process receipt schema_version is not 3")
    started_time_ns = receipt.get("started_time_ns")
    checks["receipt_start_time"] = (
        isinstance(started_time_ns, int) and started_time_ns > 0
    )
    if not checks["receipt_start_time"]:
        errors.append("process receipt lacks a valid started_time_ns")
    receipt_executable = receipt.get("executable")
    executable_ok = False
    if isinstance(receipt_executable, dict):
        resolved_executable = receipt_executable.get("resolved_path")
        executable_digest = receipt_executable.get("sha256")
        executable_ok = (
            isinstance(resolved_executable, str)
            and Path(resolved_executable).is_absolute()
            and isinstance(executable_digest, str)
            and bool(re.fullmatch(r"[0-9a-f]{64}", executable_digest))
            and isinstance(receipt_executable.get("size_bytes"), int)
            and receipt_executable.get("size_bytes") > 0
            and isinstance(receipt_executable.get("mtime_ns"), int)
            and receipt_executable.get("hash_timing")
            == "captured before process launch"
        )
    checks["receipt_executable_identity"] = executable_ok
    if not executable_ok:
        errors.append("process receipt does not identify the resolved executable and hash")
    argv = receipt.get("argv")
    checks["receipt_md_command"] = (
        isinstance(argv, list)
        and all(isinstance(value, str) for value in argv)
        and argv.count("--md") == 1
    )
    if not checks["receipt_md_command"]:
        errors.append("process receipt does not contain an exact argv list with --md")
    checks["receipt_executable_argv_binding"] = (
        executable_ok
        and isinstance(argv, list)
        and bool(argv)
        and argv[0] == receipt_executable.get("resolved_path")
    )
    if not checks["receipt_executable_argv_binding"]:
        errors.append("receipt executable identity does not match argv[0]")
    argv_strings = (
        argv
        if isinstance(argv, list) and all(isinstance(value, str) for value in argv)
        else []
    )

    expected_method_tokens = {
        "gfnff": ("--gfnff", None),
        "gfn0": ("--gfn", "0"),
        "gfn1": ("--gfn", "1"),
        "gfn2": ("--gfn", "2"),
    }
    method_option, method_value = expected_method_tokens[args.expected_method]
    method_indices = [
        index for index, token in enumerate(argv_strings) if token == method_option
    ]
    method_matches = len(method_indices) == 1
    if method_matches and method_value is not None:
        index = method_indices[0]
        method_matches = index + 1 < len(argv_strings) and argv_strings[index + 1] == method_value
    competing_method = "--gfn" if method_option == "--gfnff" else "--gfnff"
    method_matches = method_matches and competing_method not in argv_strings
    checks["argv_method_matches_expected"] = method_matches
    if not method_matches:
        errors.append("receipt argv does not uniquely match the expected xTB method")

    logged_charge = exact_option_integer(argv_strings, "--chrg")
    checks["argv_charge_matches_expected"] = logged_charge == args.expected_charge
    if not checks["argv_charge_matches_expected"]:
        errors.append("receipt argv does not uniquely match the expected total charge")
    logged_uhf = exact_option_integer(argv_strings, "--uhf")
    checks["argv_uhf_matches_expected"] = logged_uhf == args.expected_uhf
    if not checks["argv_uhf_matches_expected"]:
        errors.append("receipt argv does not uniquely match the expected UHF count")
    receipt_log = receipt.get("log")
    checks["receipt_log_binding"] = (
        isinstance(receipt_log, dict)
        and receipt_log.get("path") == str(log_path)
        and log_path.is_file()
        and receipt_log.get("sha256") == sha256_file(log_path)
    )
    if not checks["receipt_log_binding"]:
        errors.append("receipt log path/hash does not bind the validated log")

    generated_artifacts = receipt.get("generated_artifacts")

    def artifact_is_bound(name: str, path: Path) -> bool:
        if not isinstance(generated_artifacts, dict):
            return False
        record = generated_artifacts.get(name)
        if not isinstance(record, dict):
            return False
        if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
            return False
        recorded_mtime = record.get("mtime_ns")
        freshness_floor = (
            started_time_ns - 2_000_000_000
            if isinstance(started_time_ns, int)
            else None
        )
        return (
            record.get("path") == str(path)
            and record.get("regular_file") is True
            and record.get("size_bytes") == path.stat().st_size
            and recorded_mtime == path.stat().st_mtime_ns
            and isinstance(recorded_mtime, int)
            and isinstance(freshness_floor, int)
            and recorded_mtime >= freshness_floor
            and record.get("sha256") == sha256_file(path)
        )

    checks["receipt_trajectory_binding"] = artifact_is_bound(
        "trajectory", trajectory_path
    )
    if not checks["receipt_trajectory_binding"]:
        errors.append("receipt does not bind a fresh current trajectory")
    checks["receipt_success_marker_binding"] = artifact_is_bound(
        "success_marker", marker_path
    )
    if not checks["receipt_success_marker_binding"]:
        errors.append("receipt does not bind a fresh current success marker")
    receipt_inputs = receipt.get("input_files")
    expected_sources = {str(input_path), str(md_input_path)}
    observed_sources: set[str] = set()
    input_binding_ok = isinstance(receipt_inputs, list) and len(receipt_inputs) == 2
    if isinstance(receipt_inputs, list):
        for record in receipt_inputs:
            if not isinstance(record, dict):
                input_binding_ok = False
                continue
            source = record.get("source")
            source_after = record.get("source_after_run")
            snapshot = record.get("snapshot")
            if not all(isinstance(value, dict) for value in (source, source_after, snapshot)):
                input_binding_ok = False
                continue
            source_path_value = source.get("path")
            snapshot_path_value = snapshot.get("path")
            snapshot_token = snapshot.get("argv_token")
            if not all(
                isinstance(value, str)
                for value in (source_path_value, snapshot_path_value, snapshot_token)
            ):
                input_binding_ok = False
                continue
            try:
                current_source = resolve_run_leaf_no_follow(run_dir, source_path_value)
                current_snapshot = resolve_run_leaf_no_follow(run_dir, snapshot_path_value)
                token_snapshot = resolve_run_leaf_no_follow(run_dir, snapshot_token)
            except ValueError:
                input_binding_ok = False
                continue
            observed_sources.add(str(current_source))
            try:
                source_bytes, source_stat = read_stable_regular_file(current_source)
                snapshot_bytes, snapshot_stat = read_stable_regular_file(current_snapshot)
                source_regular = not current_source.is_symlink()
                snapshot_regular = not current_snapshot.is_symlink()
            except (OSError, ValueError):
                source_bytes = b""
                snapshot_bytes = b""
                source_stat = None
                snapshot_stat = None
                source_regular = False
                snapshot_regular = False
            snapshot_directory_private = (
                current_snapshot.parent == run_dir / ".xtb-input-snapshots"
                and current_snapshot.parent.is_dir()
                and not current_snapshot.parent.is_symlink()
                and current_snapshot.parent.stat().st_mode & 0o777 == 0o700
            )
            source_digest = hashlib.sha256(source_bytes).hexdigest() if source_regular else None
            snapshot_digest = hashlib.sha256(snapshot_bytes).hexdigest() if snapshot_regular else None
            source_argv_matches = 0
            snapshot_argv_matches = 0
            for token in argv_strings[1:]:
                try:
                    resolved_token = resolve_run_leaf_no_follow(run_dir, token)
                except ValueError:
                    continue
                source_argv_matches += resolved_token == current_source
                snapshot_argv_matches += resolved_token == current_snapshot
            record_ok = (
                str(current_source) in expected_sources
                and source_regular
                and snapshot_regular
                and snapshot_directory_private
                and token_snapshot == current_snapshot
                and source.get("sha256") == source_digest
                and source.get("size_bytes") == source_stat.st_size
                and source.get("device") == source_stat.st_dev
                and source.get("inode") == source_stat.st_ino
                and source.get("mtime_ns") == source_stat.st_mtime_ns
                and source.get("ctime_ns") == source_stat.st_ctime_ns
                and source_after.get("path") == str(current_source)
                and source_after.get("sha256") == source_digest
                and all(
                    source_after.get(field) == source.get(field)
                    for field in (
                        "path",
                        "sha256",
                        "size_bytes",
                        "device",
                        "inode",
                        "mtime_ns",
                        "ctime_ns",
                    )
                )
                and record.get("source_unchanged_after_run") is True
                and snapshot.get("sha256") == snapshot_digest == source_digest
                and snapshot.get("size_bytes") == snapshot_stat.st_size
                and snapshot.get("mode") == "0o400"
                and current_snapshot.stat().st_mode & 0o777 == 0o400
                and source_argv_matches == 0
                and snapshot_argv_matches == 1
                and argv_strings.count(snapshot_token) == 1
            )
            input_binding_ok = input_binding_ok and record_ok
    checks["receipt_input_binding"] = input_binding_ok and observed_sources == expected_sources
    if not checks["receipt_input_binding"]:
        errors.append(
            "receipt does not bind unchanged current sources to private read-only argv snapshots"
        )

    marker_counts = {
        "calculation_setup": len(re.findall(r"\bCalculation Setup\b", log_text)),
        "md_setup": len(re.findall(r"\bMolecular Dynamics\b", log_text)),
        "md_normal_exit": len(re.findall(r"\bnormal exit of md\(\)", log_text, re.IGNORECASE)),
        "normal_termination": len(
            re.findall(r"\bnormal termination of xtb\b", log_text, re.IGNORECASE)
        ),
    }
    for name, count in marker_counts.items():
        checks[name] = count == 1
        if count != 1:
            errors.append(f"log must contain exactly one {name.replace('_', ' ')} marker; found {count}")
    marker_order = [
        log_text.find("Calculation Setup"),
        log_text.find("Molecular Dynamics"),
        log_text.lower().find("normal exit of md()"),
        log_text.lower().find("normal termination of xtb"),
    ]
    checks["log_markers_ordered"] = (
        all(index >= 0 for index in marker_order)
        and marker_order == sorted(marker_order)
    )
    if not checks["log_markers_ordered"]:
        errors.append("calculation, MD, MD-exit, and program-exit markers are not ordered")
    checks["xtbmdok_marker"] = marker_path.is_file() and marker_path.stat().st_size > 0
    if not checks["xtbmdok_marker"]:
        errors.append("xtbmdok success marker is absent or empty")

    rejected_messages: list[str] = []
    for name, pattern in REJECT_PATTERNS.items():
        present = bool(pattern.search(log_text))
        checks[f"reject_pattern_absent:{name}"] = not present
        if present:
            rejected_messages.append(name)
    if rejected_messages:
        errors.append("rejected log diagnostics: " + ", ".join(rejected_messages))

    logged_values: dict[str, float | int | str | None] = {}
    value_specs = (
        ("time_ps", "MD time /ps", args.requested_time_ps),
        ("step_fs", "dt /fs", args.requested_step_fs),
        ("dump_fs", "dumpstep(trj) /fs", args.requested_dump_fs),
        ("temperature_k", "temperature /K", args.requested_temperature_k),
        ("scc_accuracy", "SCC accuracy", args.requested_scc_accuracy),
        ("hydrogen_mass", "H atoms mass (amu)", args.requested_hydrogen_mass),
    )
    for name, label, expected in value_specs:
        try:
            actual = parse_unique_log_value(log_text, label)
            matches = math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-8)
        except ValueError as exc:
            actual = None
            matches = False
            errors.append(str(exc))
        logged_values[name] = actual
        checks[f"logged_{name}_matches_request"] = matches
        if actual is not None and not matches:
            errors.append(f"logged {name} {actual!r} does not match requested {expected}")

    expected_steps_float = args.requested_time_ps * 1000.0 / args.requested_step_fs
    expected_steps = round(expected_steps_float)
    if not math.isclose(expected_steps_float, expected_steps, abs_tol=1.0e-8):
        raise ValueError("requested time divided by step must be an integer number of MD steps")
    try:
        logged_steps = parse_unique_log_integer(log_text, "max steps")
        checks["logged_max_steps_matches_request"] = logged_steps == expected_steps
    except ValueError as exc:
        logged_steps = None
        checks["logged_max_steps_matches_request"] = False
        errors.append(str(exc))
    if logged_steps is not None and logged_steps != expected_steps:
        errors.append(
            f"logged max steps {logged_steps} does not match derived request {expected_steps}"
        )
    logged_values["max_steps"] = logged_steps

    shake_matches = re.findall(
        r"^\s*SHAKE on\. # bonds\s*:\s*(\d+)\b", log_text, re.MULTILINE
    )
    checks["logged_shake_bonds_matches_request"] = (
        len(shake_matches) == 1 and int(shake_matches[0]) == args.requested_shake_bonds
    )
    if not checks["logged_shake_bonds_matches_request"]:
        errors.append(
            "log must contain exactly one SHAKE bond count matching the request"
        )
    logged_values["shake_bonds"] = (
        int(shake_matches[0]) if len(shake_matches) == 1 else None
    )

    thermostat_on_lines = re.findall(
        r"^\s*(.+THERMOSTAT on)\s*$", log_text, re.MULTILINE | re.IGNORECASE
    )
    thermostat_off_lines = re.findall(
        r"^\s*(.+THERMOSTAT off)\s*$", log_text, re.MULTILINE | re.IGNORECASE
    )
    if args.expect_thermostat == "on":
        thermostat_ok = len(thermostat_on_lines) == 1 and not thermostat_off_lines
        thermostat_value = (
            thermostat_on_lines[0].strip() if len(thermostat_on_lines) == 1 else None
        )
    else:
        thermostat_ok = len(thermostat_off_lines) == 1 and not thermostat_on_lines
        thermostat_value = (
            thermostat_off_lines[0].strip()
            if len(thermostat_off_lines) == 1
            else None
        )
    checks["logged_thermostat_matches_request"] = thermostat_ok
    if not thermostat_ok:
        errors.append("logged thermostat state does not match the request")
    logged_values["thermostat"] = thermostat_value

    version_matches = re.findall(
        r"^\s*\*?\s*xTB version\s+([^\s]+)\s*$", log_text, re.MULTILINE | re.IGNORECASE
    )
    checks["log_xtb_version_matches_expected"] = (
        len(version_matches) == 1 and version_matches[0] == args.expected_xtb_version
    )
    if not checks["log_xtb_version_matches_expected"]:
        errors.append("log must contain exactly one xTB version matching the expectation")

    program_calls = re.findall(r"^\s*program call\s*:\s*(.+?)\s*$", log_text, re.MULTILINE)
    checks["unique_program_call"] = len(program_calls) == 1 and "--md" in program_calls[0].split()
    if not checks["unique_program_call"]:
        errors.append("log must contain exactly one MD program call")
    try:
        program_tokens = shlex.split(program_calls[0]) if len(program_calls) == 1 else []
    except ValueError as exc:
        program_tokens = []
        errors.append(f"logged program call cannot be parsed: {exc}")
    checks["program_call_matches_receipt_argv"] = (
        bool(program_tokens)
        and bool(argv_strings)
        and Path(program_tokens[0]).name == Path(argv_strings[0]).name
        and program_tokens[1:] == argv_strings[1:]
    )
    if not checks["program_call_matches_receipt_argv"]:
        errors.append("logged program call does not exactly match the receipt argv")

    input_labels: list[str] = []
    input_comment = ""
    expected_atoms = args.expected_molecules * args.atoms_per_molecule
    try:
        input_labels, input_comment = validate_single_xyz(input_path, expected_atoms)
        checks["input_xyz_structure"] = True
    except ValueError as exc:
        checks["input_xyz_structure"] = False
        errors.append(str(exc))

    stats = TrajectoryStats()
    try:
        stats = validate_trajectory(
            trajectory_path,
            args.expected_molecules,
            args.atoms_per_molecule,
            args.minimum_intermolecular_distance,
            args.maximum_system_extent,
        )
        checks["trajectory_structure_and_physics"] = True
    except ValueError as exc:
        checks["trajectory_structure_and_physics"] = False
        errors.append(str(exc))
    checks["trajectory_matches_input_labels"] = bool(
        input_labels and stats.labels and input_labels == stats.labels
    )
    if not checks["trajectory_matches_input_labels"]:
        errors.append("trajectory atom labels/order do not match the input XYZ")

    effective_dump = (
        float(logged_values["dump_fs"])
        if logged_values.get("dump_fs") is not None
        else args.requested_dump_fs
    )
    implied_last_time_ps = max(0, stats.frame_count - 1) * effective_dump / 1000.0
    frame_tolerance_ps = effective_dump / 1000.0
    time_covered = implied_last_time_ps + frame_tolerance_ps >= args.requested_time_ps
    checks["trajectory_frame_coverage"] = stats.frame_count > 1 and time_covered
    if not checks["trajectory_frame_coverage"]:
        errors.append(
            "trajectory frame count does not cover the requested time within one dump interval"
        )

    accepted = not errors and all(checks.values())
    report = {
        "schema_version": 2,
        "status": "accepted" if accepted else "rejected",
        "run_directory": str(run_dir),
        "inputs": {
            "expected_molecules": args.expected_molecules,
            "atoms_per_molecule": args.atoms_per_molecule,
            "expected_atoms": expected_atoms,
            "expected_method": args.expected_method,
            "expected_total_charge": args.expected_charge,
            "expected_uhf": args.expected_uhf,
            "expected_xtb_version": args.expected_xtb_version,
            "requested_time_ps": args.requested_time_ps,
            "requested_step_fs": args.requested_step_fs,
            "requested_dump_fs": args.requested_dump_fs,
            "requested_temperature_k": args.requested_temperature_k,
            "requested_scc_accuracy": args.requested_scc_accuracy,
            "requested_hydrogen_mass": args.requested_hydrogen_mass,
            "requested_shake_bonds": args.requested_shake_bonds,
            "expected_thermostat": args.expect_thermostat,
            "minimum_inter_distance_angstrom": args.minimum_intermolecular_distance,
            "maximum_system_extent_angstrom": args.maximum_system_extent,
        },
        "receipt": {
            "path": str(receipt_path),
            "sha256": sha256_file(receipt_path) if receipt_path.is_file() else None,
            "started_time_ns": started_time_ns,
            "returncode": receipt.get("returncode"),
            "signal": receipt.get("signal"),
            "argv": argv,
            "executable": receipt_executable,
            "generated_artifacts": generated_artifacts,
        },
        "log": {
            "path": str(log_path),
            "sha256": sha256_file(log_path) if log_path.is_file() else None,
            "program_call": program_calls[0] if len(program_calls) == 1 else None,
            "logged_settings": logged_values,
        },
        "input_xyz": {
            "path": str(input_path),
            "sha256": sha256_file(input_path) if input_path.is_file() else None,
            "comment": input_comment,
        },
        "md_input": {
            "path": str(md_input_path),
            "sha256": sha256_file(md_input_path) if md_input_path.is_file() else None,
        },
        "trajectory": {
            "path": str(trajectory_path),
            "sha256": sha256_file(trajectory_path) if trajectory_path.is_file() else None,
            **asdict(stats),
            "implied_last_time_ps": implied_last_time_ps,
            "time_basis": "derived from frame count and the unique logged dump interval; tolerance is one dump interval",
        },
        "checks": checks,
        "errors": errors,
    }
    with report_path.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"{report['status'].upper()}: {report_path}")
    return 0 if accepted else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
