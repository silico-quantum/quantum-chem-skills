"""Strict accepted-trajectory loading and atomic GIF publication helpers."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import TextIO


ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$|^X$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_required_line(handle: TextIO, description: str) -> str:
    line = handle.readline()
    if line == "":
        raise ValueError(f"trajectory ended while reading {description}")
    return line


def iter_strict_xyz_frames(path: Path, expected_atoms: int):
    """Yield every complete frame while rejecting malformed data and label drift."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("trajectory is absent or empty")
    reference_labels: list[str] | None = None
    frame_number = 0
    with path.open(encoding="utf-8") as handle:
        while True:
            count_line = handle.readline()
            if count_line == "":
                break
            frame_number += 1
            if not count_line.strip():
                raise ValueError("blank lines between XYZ frames are not accepted")
            try:
                atom_count = int(count_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"frame {frame_number} atom count is not an integer"
                ) from exc
            if atom_count != expected_atoms:
                raise ValueError(
                    f"frame {frame_number} has {atom_count} atoms; expected {expected_atoms}"
                )
            comment = read_required_line(handle, f"frame {frame_number} comment").rstrip(
                "\n"
            )
            labels: list[str] = []
            coordinates: list[tuple[float, float, float]] = []
            for atom_index in range(1, atom_count + 1):
                fields = read_required_line(
                    handle, f"frame {frame_number} atom {atom_index}"
                ).split()
                if len(fields) != 4 or not ELEMENT_RE.fullmatch(fields[0]):
                    raise ValueError(
                        f"frame {frame_number} atom {atom_index} is malformed"
                    )
                try:
                    point = tuple(float(value) for value in fields[1:])
                except ValueError as exc:
                    raise ValueError(
                        f"frame {frame_number} atom {atom_index} is non-numeric"
                    ) from exc
                if not all(math.isfinite(value) for value in point):
                    raise ValueError(
                        f"frame {frame_number} atom {atom_index} is non-finite"
                    )
                labels.append(fields[0])
                coordinates.append(point)
            if reference_labels is None:
                reference_labels = labels
            elif labels != reference_labels:
                raise ValueError(
                    f"atom labels or ordering changed in frame {frame_number}"
                )
            yield frame_number - 1, comment, labels, coordinates
    if frame_number == 0:
        raise ValueError("trajectory contains no complete frame")


def selected_indices(frame_count: int, stride: int, max_frames: int) -> list[int]:
    candidates = list(range(0, frame_count, stride))
    if candidates[-1] != frame_count - 1:
        candidates.append(frame_count - 1)
    if len(candidates) <= max_frames:
        return candidates
    if max_frames == 1:
        return [frame_count - 1]
    return [
        candidates[round(index * (len(candidates) - 1) / (max_frames - 1))]
        for index in range(max_frames)
    ]


def load_accepted_sampled_trajectory(
    trajectory: Path,
    validation: Path,
    molecules: int,
    atoms_per_molecule: int,
    stride: int,
    max_frames: int,
) -> tuple[list[tuple[str, list[str], list[tuple[float, float, float]]]], dict[str, object]]:
    """Verify the accepted report/hash, scan all frames, then load a bounded sample."""
    if not validation.is_file() or validation.stat().st_size == 0:
        raise ValueError("validation report is absent or empty")
    report = json.loads(validation.read_text(encoding="utf-8"))
    if not isinstance(report, dict) or report.get("status") != "accepted":
        raise ValueError("validation report does not mark the trajectory accepted")
    checks = report.get("checks")
    if not isinstance(checks, dict) or not checks or not all(
        value is True for value in checks.values()
    ):
        raise ValueError("validation report contains a failed or non-boolean check")
    report_inputs = report.get("inputs")
    report_trajectory = report.get("trajectory")
    if not isinstance(report_inputs, dict) or not isinstance(report_trajectory, dict):
        raise ValueError("validation report is missing inputs or trajectory metadata")
    if (
        report_inputs.get("expected_molecules") != molecules
        or report_inputs.get("atoms_per_molecule") != atoms_per_molecule
    ):
        raise ValueError("animation molecule counts do not match the validation report")
    expected_hash = report_trajectory.get("sha256")
    actual_hash = sha256_file(trajectory)
    reported_path = report_trajectory.get("path")
    if not isinstance(reported_path, str) or Path(reported_path).resolve() != trajectory.resolve():
        raise ValueError("trajectory path does not match the accepted validation report")
    if not isinstance(expected_hash, str) or expected_hash != actual_hash:
        raise ValueError("trajectory hash does not match the accepted validation report")

    expected_atoms = molecules * atoms_per_molecule
    frame_count = 0
    reference_labels: list[str] | None = None
    for frame_index, _, labels, _ in iter_strict_xyz_frames(
        trajectory, expected_atoms
    ):
        frame_count = frame_index + 1
        if reference_labels is None:
            reference_labels = labels
            block = labels[:atoms_per_molecule]
            for molecule_index in range(molecules):
                start = molecule_index * atoms_per_molecule
                if labels[start : start + atoms_per_molecule] != block:
                    raise ValueError(
                        "trajectory molecule blocks do not share one atom-label ordering"
                    )
    if report_trajectory.get("frame_count") != frame_count:
        raise ValueError("trajectory frame count does not match the validation report")

    indices = selected_indices(frame_count, stride, max_frames)
    selected = set(indices)
    frames = [
        (comment, labels, coordinates)
        for frame_index, comment, labels, coordinates in iter_strict_xyz_frames(
            trajectory, expected_atoms
        )
        if frame_index in selected
    ]
    if len(frames) != len(indices):
        raise RuntimeError("failed to reload every selected trajectory frame")
    provenance = {
        "trajectory_path": str(trajectory.resolve()),
        "trajectory_sha256": actual_hash,
        "validation_path": str(validation.resolve()),
        "validation_sha256": sha256_file(validation),
        "total_frame_count": frame_count,
        "selected_frame_indices_zero_based": indices,
    }
    return frames, provenance


def save_verified_gif(
    frames,
    output: Path,
    manifest: Path,
    duration_ms: int,
    provenance: dict[str, object],
    rendering: dict[str, object],
) -> None:
    """Write a GIF through a temporary file, reopen it, and publish a manifest."""
    from PIL import Image

    if not frames:
        raise ValueError("cannot write a GIF without frames")
    for path in (output, manifest):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"refusing to overwrite output: {path}")
        if not path.parent.is_dir():
            raise FileNotFoundError(f"output parent does not exist: {path.parent}")
    temporary = output.with_name(f".{output.name}.partial")
    temporary_manifest = manifest.with_name(f".{manifest.name}.partial")
    for path in (temporary, temporary_manifest):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"stale partial output exists: {path}")
    output_published = False
    try:
        frames[0].save(
            temporary,
            save_all=True,
            append_images=frames[1:],
            duration=duration_ms,
            loop=0,
            optimize=False,
            format="GIF",
        )
        with Image.open(temporary) as image:
            actual_frames = getattr(image, "n_frames", 1)
            dimensions = list(image.size)
            if image.format != "GIF" or actual_frames != len(frames):
                raise RuntimeError("reopened GIF format/frame count failed validation")
            if dimensions[0] <= 0 or dimensions[1] <= 0:
                raise RuntimeError("reopened GIF has invalid dimensions")
        receipt = {
            "schema_version": 1,
            "status": "accepted",
            "source": provenance,
            "rendering": rendering,
            "output": {
                "path": str(output.resolve()),
                "sha256": sha256_file(temporary),
                "size_bytes": temporary.stat().st_size,
                "frame_count": actual_frames,
                "dimensions_pixels": dimensions,
                "duration_ms_per_frame": duration_ms,
            },
        }
        with temporary_manifest.open("x", encoding="utf-8") as handle:
            json.dump(receipt, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, output)
        output_published = True
        try:
            os.link(temporary_manifest, manifest)
        except Exception:
            if output_published and output.exists() and not output.is_symlink():
                temporary_stat = temporary.stat()
                output_stat = output.stat()
                if (
                    output_stat.st_dev == temporary_stat.st_dev
                    and output_stat.st_ino == temporary_stat.st_ino
                ):
                    output.unlink()
            raise
    finally:
        for path in (temporary, temporary_manifest):
            if path.exists():
                path.unlink()
