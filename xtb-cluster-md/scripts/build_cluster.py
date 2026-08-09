#!/usr/bin/env python3
"""Build a deterministic molecular cluster with an atom-pair distance cutoff."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path


ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$|^X$")


@dataclass(frozen=True)
class SdfRecord:
    title: str
    program_line: str
    comment: str
    atoms: list[str]
    coordinates: list[tuple[float, float, float]]
    bonds: list[tuple[int, int, int]]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def parse_sdf_text(text: str) -> SdfRecord:
    """Parse exactly one complete V2000 SDF record from captured text."""
    lines = text.splitlines()
    delimiters = [index for index, line in enumerate(lines) if line.strip() == "$$$$"]
    if len(delimiters) != 1:
        raise ValueError("SDF must contain exactly one record and one $$$$ delimiter")
    delimiter = delimiters[0]
    if any(line.strip() for line in lines[delimiter + 1 :]):
        raise ValueError("SDF contains content after its single record delimiter")
    lines = lines[:delimiter]
    if len(lines) < 5:
        raise ValueError("SDF is too short")
    if "V2000" not in lines[3]:
        raise ValueError("only a V2000 SDF atom block is supported")
    try:
        atom_count = int(lines[3][0:3])
        bond_count = int(lines[3][3:6])
    except ValueError as exc:
        raise ValueError("invalid V2000 atom or bond count") from exc
    if atom_count <= 0 or bond_count < 0:
        raise ValueError("SDF atom count must be positive and bond count non-negative")
    if len(lines) < 4 + atom_count + bond_count + 1:
        raise ValueError("SDF atom/bond block is truncated")

    atoms: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for offset, line in enumerate(lines[4 : 4 + atom_count], start=5):
        if len(line) < 34:
            raise ValueError(f"truncated SDF atom line {offset}")
        try:
            xyz = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
        except ValueError as exc:
            raise ValueError(f"invalid coordinate on SDF atom line {offset}") from exc
        symbol = line[31:34].strip()
        if not ELEMENT_RE.fullmatch(symbol):
            raise ValueError(f"invalid element on SDF atom line {offset}: {symbol}")
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"non-finite coordinate on SDF atom line {offset}")
        atoms.append(symbol)
        coords.append(xyz)

    bonds: list[tuple[int, int, int]] = []
    bond_start = 4 + atom_count
    for offset, line in enumerate(
        lines[bond_start : bond_start + bond_count], start=bond_start + 1
    ):
        if len(line) < 9:
            raise ValueError(f"truncated SDF bond line {offset}")
        try:
            first = int(line[0:3])
            second = int(line[3:6])
            order = int(line[6:9])
        except ValueError as exc:
            raise ValueError(f"invalid SDF bond line {offset}") from exc
        if not (1 <= first <= atom_count and 1 <= second <= atom_count):
            raise ValueError(f"bond atom index out of range on SDF line {offset}")
        if first == second or order <= 0:
            raise ValueError(f"invalid bond on SDF line {offset}")
        bonds.append((first - 1, second - 1, order))

    property_lines = lines[bond_start + bond_count :]
    end_markers = [line for line in property_lines if line.strip() == "M  END"]
    if len(end_markers) != 1:
        raise ValueError("SDF record must contain exactly one M  END marker")

    return SdfRecord(
        title=lines[0].strip(),
        program_line=lines[1].rstrip(),
        comment=lines[2].rstrip(),
        atoms=atoms,
        coordinates=coords,
        bonds=bonds,
    )


def read_sdf_record(path: Path) -> SdfRecord:
    """Read and parse one V2000 SDF record."""
    return parse_sdf_text(path.read_text(encoding="utf-8"))


def read_sdf_coords(path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    """Compatibility wrapper returning coordinates from one strict SDF record."""
    record = read_sdf_record(path)
    return record.atoms, record.coordinates


def coordinate_centroid(
    coords: list[tuple[float, float, float]],
) -> tuple[float, float, float]:
    count = len(coords)
    return tuple(sum(point[axis] for point in coords) / count for axis in range(3))


def random_rotation_matrix(rng: random.Random) -> tuple[tuple[float, ...], ...]:
    """Return a uniformly sampled rotation using a random unit quaternion."""
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    x = math.sqrt(1.0 - u1) * math.sin(2.0 * math.pi * u2)
    y = math.sqrt(1.0 - u1) * math.cos(2.0 * math.pi * u2)
    z = math.sqrt(u1) * math.sin(2.0 * math.pi * u3)
    w = math.sqrt(u1) * math.cos(2.0 * math.pi * u3)
    return (
        (1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
        (2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
        (2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)),
    )


def rotate_translate(
    coords: list[tuple[float, float, float]],
    rotation: tuple[tuple[float, ...], ...],
    translation: tuple[float, float, float],
) -> list[tuple[float, float, float]]:
    return [
        tuple(
            sum(point[column] * rotation[row][column] for column in range(3))
            + translation[row]
            for row in range(3)
        )
        for point in coords
    ]


def minimum_pair(
    left: list[tuple[float, float, float]],
    right: list[tuple[float, float, float]],
) -> tuple[float, int, int]:
    minimum_squared = math.inf
    pair = (-1, -1)
    for first_index, first in enumerate(left):
        for second_index, second in enumerate(right):
            distance_squared = sum((a - b) ** 2 for a, b in zip(first, second))
            if distance_squared < minimum_squared:
                minimum_squared = distance_squared
                pair = (first_index, second_index)
    return math.sqrt(minimum_squared), pair[0], pair[1]


def minimum_pair_distance(
    left: list[tuple[float, float, float]],
    right: list[tuple[float, float, float]],
) -> float:
    return minimum_pair(left, right)[0]


def clears_cutoff(
    candidate: list[tuple[float, float, float]],
    placements: list[list[tuple[float, float, float]]],
    cutoff: float,
) -> bool:
    cutoff_squared = cutoff * cutoff
    for placed in placements:
        for first in candidate:
            for second in placed:
                if sum((a - b) ** 2 for a, b in zip(first, second)) < cutoff_squared:
                    return False
    return True


def nearest_placement_distance(
    candidate: list[tuple[float, float, float]],
    placements: list[list[tuple[float, float, float]]],
) -> float:
    return min(minimum_pair_distance(candidate, placed) for placed in placements)


def write_xyz(
    path: Path,
    atoms: list[str],
    placements: list[list[tuple[float, float, float]]],
    seed: int,
    box: float,
    requested_cutoff: float,
    actual_minimum: float | None,
) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(f"{len(atoms) * len(placements)}\n")
        handle.write(
            f"connected_placement_ensemble molecules={len(placements)} seed={seed} "
            f"centroid_translation_box_A={box} "
            f"requested_min_inter_distance_A={requested_cutoff} "
            f"actual_min_inter_distance_A={actual_minimum}\n"
        )
        for coords in placements:
            for symbol, point in zip(atoms, coords):
                handle.write(
                    f"{symbol:2s} {point[0]:14.8f} {point[1]:14.8f} {point[2]:14.8f}\n"
                )


def validate_xyz_output(
    path: Path, expected_atoms: list[str], molecule_count: int
) -> list[tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    expected_count = len(expected_atoms) * molecule_count
    if len(lines) != expected_count + 2:
        raise RuntimeError("written XYZ has an unexpected line count")
    try:
        declared_count = int(lines[0].strip())
    except ValueError as exc:
        raise RuntimeError("written XYZ atom count is invalid") from exc
    if declared_count != expected_count:
        raise RuntimeError("written XYZ atom count does not match the build request")
    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for atom_index, line in enumerate(lines[2:], start=1):
        fields = line.split()
        if len(fields) != 4:
            raise RuntimeError(f"written XYZ atom {atom_index} is malformed")
        try:
            xyz = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise RuntimeError(f"written XYZ atom {atom_index} is non-numeric") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise RuntimeError(f"written XYZ atom {atom_index} is non-finite")
        labels.append(fields[0])
        coordinates.append(xyz)
    if labels != expected_atoms * molecule_count:
        raise RuntimeError("written XYZ atom labels/order do not match the monomer blocks")
    return coordinates


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdf", type=Path, required=True, help="verified 3D V2000 SDF")
    parser.add_argument("-n", "--molecules", type=int, required=True)
    parser.add_argument(
        "--box",
        type=float,
        required=True,
        help="cubic centroid-translation sampling range in Angstrom",
    )
    parser.add_argument(
        "--min-distance",
        type=float,
        required=True,
        help="minimum intermolecular atom-atom distance in Angstrom",
    )
    parser.add_argument(
        "--max-neighbor-distance",
        type=float,
        required=True,
        help="maximum atom-pair distance connecting each new molecule to the ensemble",
    )
    parser.add_argument("--monomer-id", required=True)
    parser.add_argument(
        "--coordinate-provenance",
        required=True,
        help="source or method establishing that the SDF coordinates are 3D",
    )
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--max-tries", type=int, default=20000)
    parser.add_argument("-o", "--out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.molecules <= 0
        or not math.isfinite(args.box)
        or args.box <= 0
        or not math.isfinite(args.min_distance)
        or args.min_distance <= 0
        or not math.isfinite(args.max_neighbor_distance)
        or args.max_neighbor_distance <= 0
    ):
        raise ValueError(
            "molecules, box, min-distance, and max-neighbor-distance must be positive"
        )
    if args.max_neighbor_distance < args.min_distance:
        raise ValueError("max-neighbor-distance cannot be smaller than min-distance")
    if not args.monomer_id.strip() or not args.coordinate_provenance.strip():
        raise ValueError("monomer-id and coordinate-provenance cannot be blank")
    if args.max_tries <= 0:
        raise ValueError("max-tries must be positive")
    output_partial = args.out.with_name(f".{args.out.name}.partial")
    manifest_partial = args.manifest.with_name(f".{args.manifest.name}.partial")
    incomplete_marker = args.manifest.with_name(
        f".{args.manifest.name}.RUN_INCOMPLETE"
    )
    for output in (
        args.out,
        args.manifest,
        output_partial,
        manifest_partial,
        incomplete_marker,
    ):
        if output.exists():
            raise FileExistsError(f"refusing to overwrite output: {output}")
        if not output.parent.is_dir():
            raise FileNotFoundError(f"output parent does not exist: {output.parent}")
    incomplete_marker.touch(exist_ok=False)

    if not args.sdf.is_file() or args.sdf.stat().st_size == 0:
        raise ValueError("SDF input is absent or empty")
    source_bytes = args.sdf.read_bytes()
    try:
        source_text = source_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("SDF input is not valid UTF-8") from exc
    source_sha256 = sha256_bytes(source_bytes)
    record = parse_sdf_text(source_text)
    atoms, source_coords = record.atoms, record.coordinates
    unique_source_pairs = [
        (
            math.dist(source_coords[first], source_coords[second]),
            first,
            second,
        )
        for first in range(len(source_coords))
        for second in range(first + 1, len(source_coords))
    ]
    if unique_source_pairs:
        source_minimum, source_first, source_second = min(unique_source_pairs)
        if source_minimum <= 1.0e-6:
            raise ValueError("SDF contains duplicate intramolecular atom coordinates")
    else:
        source_minimum, source_first, source_second = None, None, None
    centroid = coordinate_centroid(source_coords)
    centered = [
        tuple(value - centroid[axis] for axis, value in enumerate(point))
        for point in source_coords
    ]
    rng = random.Random(args.seed)
    placements: list[list[tuple[float, float, float]]] = []
    attempts_per_molecule: list[int] = []

    for molecule_index in range(args.molecules):
        for attempt in range(1, args.max_tries + 1):
            rotation = random_rotation_matrix(rng)
            translation = tuple((rng.random() - 0.5) * args.box for _ in range(3))
            candidate = rotate_translate(centered, rotation, translation)
            connected = (
                not placements
                or nearest_placement_distance(candidate, placements)
                <= args.max_neighbor_distance
            )
            if connected and clears_cutoff(candidate, placements, args.min_distance):
                placements.append(candidate)
                attempts_per_molecule.append(attempt)
                break
        else:
            raise RuntimeError(
                f"failed to place molecule {molecule_index + 1}/{args.molecules}; "
                "increase --box/--max-neighbor-distance, reduce --min-distance, "
                "or increase --max-tries"
            )

    intermolecular_pairs = [
        (minimum_pair(placements[i], placements[j]), i, j)
        for i in range(len(placements))
        for j in range(i + 1, len(placements))
    ]
    if intermolecular_pairs:
        (actual_minimum, first_atom, second_atom), first_molecule, second_molecule = min(
            intermolecular_pairs, key=lambda item: item[0][0]
        )
        nearest_by_molecule = []
        for molecule_index in range(len(placements)):
            nearest_by_molecule.append(
                min(
                    minimum_pair_distance(placements[molecule_index], placements[other])
                    for other in range(len(placements))
                    if other != molecule_index
                )
            )
        largest_nearest_neighbor = max(nearest_by_molecule)
    else:
        actual_minimum = None
        first_atom = second_atom = first_molecule = second_molecule = None
        largest_nearest_neighbor = None
    if actual_minimum is not None and actual_minimum + 1.0e-10 < args.min_distance:
        raise RuntimeError("internal validation found an intermolecular cutoff violation")
    if (
        largest_nearest_neighbor is not None
        and largest_nearest_neighbor > args.max_neighbor_distance + 1.0e-10
    ):
        raise RuntimeError("internal validation found a disconnected molecule placement")

    write_xyz(
        output_partial,
        atoms,
        placements,
        args.seed,
        args.box,
        args.min_distance,
        actual_minimum,
    )
    written_coordinates = validate_xyz_output(output_partial, atoms, args.molecules)
    written_placements = [
        written_coordinates[index : index + len(atoms)]
        for index in range(0, len(written_coordinates), len(atoms))
    ]
    written_pairs = [
        (minimum_pair(written_placements[i], written_placements[j]), i, j)
        for i in range(len(written_placements))
        for j in range(i + 1, len(written_placements))
    ]
    if written_pairs:
        (
            written_minimum,
            written_first_atom,
            written_second_atom,
        ), written_first_molecule, written_second_molecule = min(
            written_pairs, key=lambda item: item[0][0]
        )
        if written_minimum + 1.0e-10 < args.min_distance:
            raise RuntimeError("rounded XYZ violates the intermolecular cutoff")
    else:
        written_minimum = None
        written_first_atom = written_second_atom = None
        written_first_molecule = written_second_molecule = None
    manifest = {
        "status": "accepted",
        "runtime": {
            "python_version": platform.python_version(),
            "builder_sha256": sha256_file(Path(__file__)),
            "random_generator": "python random.Random with explicit seed",
        },
        "source": {
            "path": str(args.sdf),
            "sha256": source_sha256,
            "hash_basis": "the exact captured bytes decoded and parsed for this build",
            "format": "SDF V2000",
            "coordinate_unit": "Angstrom",
            "monomer_id": args.monomer_id,
            "coordinate_provenance": args.coordinate_provenance,
            "title": record.title,
            "program_line": record.program_line,
            "comment": record.comment,
            "atoms_per_molecule": len(atoms),
            "atom_labels": atoms,
            "bond_count": len(record.bonds),
            "bonds_zero_based": [list(bond) for bond in record.bonds],
            "minimum_intramolecular_distance_angstrom": source_minimum,
            "minimum_intramolecular_pair_zero_based": (
                [source_first, source_second]
                if source_first is not None and source_second is not None
                else None
            ),
        },
        "sampling": {
            "molecule_count": args.molecules,
            "centroid_translation_box_angstrom": args.box,
            "requested_minimum_interatomic_distance_angstrom": args.min_distance,
            "actual_minimum_interatomic_distance_angstrom": written_minimum,
            "pre_rounding_minimum_interatomic_distance_angstrom": actual_minimum,
            "actual_minimum_pair": (
                {
                    "molecules_zero_based": [
                        written_first_molecule,
                        written_second_molecule,
                    ],
                    "atoms_within_molecule_zero_based": [
                        written_first_atom,
                        written_second_atom,
                    ],
                    "elements": [
                        atoms[written_first_atom],
                        atoms[written_second_atom],
                    ],
                }
                if written_first_atom is not None and written_second_atom is not None
                else None
            ),
            "requested_maximum_neighbor_distance_angstrom": args.max_neighbor_distance,
            "largest_nearest_neighbor_distance_angstrom": largest_nearest_neighbor,
            "connectivity": "each molecule after the first was placed within the declared atom-pair neighbor distance of an earlier molecule",
            "seed": args.seed,
            "max_tries_per_molecule": args.max_tries,
            "attempts_per_molecule": attempts_per_molecule,
            "centering": "arithmetic coordinate centroid; not center of mass",
        },
        "output": {
            "path": str(args.out),
            "sha256": sha256_file(output_partial),
            "atom_count": len(atoms) * args.molecules,
        },
        "publication": {
            "incomplete_marker": str(incomplete_marker),
            "acceptance_requires_marker_absent": True,
        },
    }
    with manifest_partial.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.link(output_partial, args.out)
    os.link(manifest_partial, args.manifest)
    output_partial.unlink()
    manifest_partial.unlink()
    incomplete_marker.unlink()
    print(f"Wrote {args.out} with {args.molecules} molecules")
    if actual_minimum is not None:
        print(f"Minimum intermolecular atom distance: {actual_minimum:.6f} Angstrom")
    print(f"Manifest: {args.manifest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
