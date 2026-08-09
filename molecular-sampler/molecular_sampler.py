#!/usr/bin/env python3
"""Extract molecular fragments and deterministic nearest-neighbor samples."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


COVALENT_RADII = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
    "Se": 1.20,
}
CHARGE_MULTIPLICITY_PATTERN = re.compile(r"^[+-]?\d+(?:\s+[+-]?\d+)+$")
ELEMENT_PATTERN = re.compile(r"^([A-Z][a-z]?)")
NON_ANGSTROM_GAUSSIAN_UNITS_PATTERN = re.compile(
    r"\bunits\s*=\s*(?:bohr|au)\b", re.IGNORECASE
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def normalize_element(token: str, *, line_number: int) -> str:
    match = ELEMENT_PATTERN.match(token)
    if not match:
        raise ValueError(f"Line {line_number}: invalid atom token {token!r}")
    element = match.group(1)
    if element not in COVALENT_RADII:
        raise ValueError(
            f"Line {line_number}: no bundled covalent radius for element {element}"
        )
    return element


def parse_xyz_text(text: str) -> list[dict[str, object]]:
    """Parse exactly one strict XYZ frame in angstrom from captured text."""
    lines = text.splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ input must contain an atom count and comment line")

    try:
        expected_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("XYZ first line must be an integer atom count") from exc
    if expected_count < 1:
        raise ValueError("XYZ atom count must be positive")

    coordinate_lines = [
        (line_number, line)
        for line_number, line in enumerate(lines[2:], start=3)
        if line.strip()
    ]
    if len(coordinate_lines) != expected_count:
        raise ValueError(
            "XYZ atom-count mismatch: "
            f"header declares {expected_count}, found {len(coordinate_lines)} records"
        )

    atoms: list[dict[str, object]] = []
    for source_index, (line_number, line) in enumerate(coordinate_lines):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Line {line_number}: incomplete XYZ coordinate")
        element = normalize_element(parts[0], line_number=line_number)
        try:
            x, y, z = (float(parts[index]) for index in range(1, 4))
        except ValueError as exc:
            raise ValueError(f"Line {line_number}: invalid XYZ coordinate") from exc
        if not all(math.isfinite(value) for value in (x, y, z)):
            raise ValueError(f"Line {line_number}: non-finite XYZ coordinate")
        atoms.append(
            {
                "element": element,
                "layer": "L",
                "x": x,
                "y": y,
                "z": z,
                "source_index": source_index,
            }
        )
    return atoms


def parse_xyz_file(filepath: str | Path) -> list[dict[str, object]]:
    """Parse one file; callers needing provenance should capture bytes first."""
    return parse_xyz_text(Path(filepath).read_text(encoding="utf-8"))


def _is_charge_multiplicity_line(line: str) -> bool:
    if not CHARGE_MULTIPLICITY_PATTERN.fullmatch(line):
        return False
    values = [int(value) for value in line.split()]
    return len(values) >= 2 and len(values) % 2 == 0 and all(
        values[index] > 0 for index in range(1, len(values), 2)
    )


def parse_gjf_text(text: str) -> list[dict[str, object]]:
    """Parse one standard or ONIOM Gaussian Cartesian geometry block."""
    lines = text.splitlines()

    if NON_ANGSTROM_GAUSSIAN_UNITS_PATTERN.search(text):
        raise ValueError(
            "Gaussian inputs using bohr/atomic units are unsupported; "
            "convert the Cartesian geometry to angstrom first"
        )

    if any(line.strip().lower() == "--link1--" for line in lines):
        raise ValueError(
            "Input must contain exactly one Gaussian job; --Link1-- is unsupported"
        )
    charge_line_indices = [
        index
        for index, line in enumerate(lines)
        if _is_charge_multiplicity_line(line.strip())
    ]
    if not charge_line_indices:
        raise ValueError("Gaussian input has no recognizable charge/multiplicity line")
    if len(charge_line_indices) != 1:
        raise ValueError(
            "Input must contain exactly one Gaussian job and geometry block"
        )
    charge_line_index = charge_line_indices[0]

    atoms: list[dict[str, object]] = []
    for index in range(charge_line_index + 1, len(lines)):
        line_number = index + 1
        line = lines[index].strip()
        if not line:
            if atoms:
                break
            continue

        parts = line.split()
        atom_token = parts[0]
        element = normalize_element(atom_token, line_number=line_number)

        is_oniom = (
            len(parts) >= 6
            and re.fullmatch(r"[+-]?\d+", parts[1]) is not None
            and parts[5].upper() in {"H", "M", "L"}
        )
        coordinate_offset = 2 if is_oniom else 1
        if len(parts) < coordinate_offset + 3:
            raise ValueError(f"Line {line_number}: incomplete Gaussian coordinate")
        try:
            x, y, z = (
                float(parts[coordinate_offset + axis]) for axis in range(3)
            )
        except ValueError as exc:
            raise ValueError(f"Line {line_number}: invalid Gaussian coordinate") from exc
        if not all(math.isfinite(value) for value in (x, y, z)):
            raise ValueError(f"Line {line_number}: non-finite Gaussian coordinate")

        atoms.append(
            {
                "element": element,
                "layer": parts[5].upper() if is_oniom else "L",
                "x": x,
                "y": y,
                "z": z,
                "source_index": len(atoms),
            }
        )

    if not atoms:
        raise ValueError("Gaussian geometry block is empty")
    return atoms


def parse_gjf_file(filepath: str | Path) -> list[dict[str, object]]:
    """Parse one file; callers needing provenance should capture bytes first."""
    return parse_gjf_text(Path(filepath).read_text(encoding="utf-8"))


def bond_distance(element_a: str, element_b: str, scale: float = 1.3) -> float:
    """Return the covalent-radius distance cutoff in angstrom."""
    if scale <= 0:
        raise ValueError("Bond-distance scale must be positive")
    try:
        return (COVALENT_RADII[element_a] + COVALENT_RADII[element_b]) * scale
    except KeyError as exc:
        raise ValueError(f"No bundled covalent radius for element {exc.args[0]}") from exc


def atom_distance(atom_a: dict[str, object], atom_b: dict[str, object]) -> float:
    return math.sqrt(
        (float(atom_a["x"]) - float(atom_b["x"])) ** 2
        + (float(atom_a["y"]) - float(atom_b["y"])) ** 2
        + (float(atom_a["z"]) - float(atom_b["z"])) ** 2
    )


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def identify_molecules(
    atoms: list[dict[str, object]],
    layer_filter: str = "L",
    *,
    bond_scale: float = 1.3,
    min_fragment_atoms: int = 1,
) -> list[dict[str, object]]:
    """Identify connected components without silently losing atoms."""
    if layer_filter not in {"H", "M", "L", "all"}:
        raise ValueError(f"Unsupported layer filter: {layer_filter}")
    if min_fragment_atoms < 1:
        raise ValueError("Minimum fragment size must be at least one atom")
    if not math.isfinite(bond_scale) or bond_scale <= 0:
        raise ValueError("Bond-distance scale must be finite and positive")

    filtered_atoms = [
        atom
        for atom in atoms
        if layer_filter == "all" or atom["layer"] == layer_filter
    ]
    if not filtered_atoms:
        return []

    union_find = UnionFind(len(filtered_atoms))
    for left in range(len(filtered_atoms)):
        for right in range(left + 1, len(filtered_atoms)):
            cutoff = bond_distance(
                str(filtered_atoms[left]["element"]),
                str(filtered_atoms[right]["element"]),
                scale=bond_scale,
            )
            if atom_distance(filtered_atoms[left], filtered_atoms[right]) < cutoff:
                union_find.union(left, right)

    components: dict[int, list[dict[str, object]]] = defaultdict(list)
    for index, atom in enumerate(filtered_atoms):
        components[union_find.find(index)].append(atom)

    component_atoms = sorted(
        components.values(),
        key=lambda component: min(int(atom["source_index"]) for atom in component),
    )
    molecules: list[dict[str, object]] = []
    for component in component_atoms:
        if len(component) < min_fragment_atoms:
            continue
        center = tuple(
            sum(float(atom[axis]) for atom in component) / len(component)
            for axis in ("x", "y", "z")
        )
        molecules.append(
            {
                "id": len(molecules),
                "atoms": component,
                "elements": dict(Counter(str(atom["element"]) for atom in component)),
                "count": len(component),
                "center": center,
                "source_indices": [int(atom["source_index"]) for atom in component],
            }
        )
    return molecules


def center_distance(left: dict[str, object], right: dict[str, object]) -> float:
    left_center = tuple(float(value) for value in left["center"])
    right_center = tuple(float(value) for value in right["center"])
    return math.dist(left_center, right_center)


def sample_molecules(
    molecules: list[dict[str, object]], n_samples: int = 20
) -> dict[str, list[list[dict[str, object]]]]:
    """Build unique, deterministic centroid-nearest samples."""
    if n_samples < 1:
        raise ValueError("Number of samples must be positive")

    sample_names = {
        1: "monomers",
        2: "dimers",
        3: "trimers",
        4: "tetramers",
        5: "pentamers",
    }
    samples: dict[str, list[list[dict[str, object]]]] = {
        name: [] for name in sample_names.values()
    }
    samples["monomers"] = [[molecule] for molecule in molecules]

    neighbors: dict[int, list[tuple[int, float]]] = {
        index: [] for index in range(len(molecules))
    }
    for left in range(len(molecules)):
        for right in range(left + 1, len(molecules)):
            distance = center_distance(molecules[left], molecules[right])
            neighbors[left].append((right, distance))
            neighbors[right].append((left, distance))
    for entries in neighbors.values():
        entries.sort(key=lambda item: (item[1], item[0]))

    for target_size in range(2, 6):
        name = sample_names[target_size]
        if len(molecules) < target_size:
            continue
        seen: set[tuple[int, ...]] = set()
        for start_index in range(len(molecules)):
            member_indices = [start_index]
            member_indices.extend(
                neighbor_index
                for neighbor_index, _ in neighbors[start_index][
                    : target_size - 1
                ]
            )
            key = tuple(sorted(member_indices))
            if len(key) != target_size or key in seen:
                continue
            seen.add(key)
            samples[name].append([molecules[index] for index in key])
            if len(samples[name]) >= n_samples:
                break
    return samples


def write_xyz(
    molecules: list[dict[str, object]], filename: str | Path, comment: str
) -> None:
    atoms = [atom for molecule in molecules for atom in molecule["atoms"]]
    lines = [str(len(atoms)), comment]
    lines.extend(
        f"{atom['element']:2s} {float(atom['x']):16.10f} "
        f"{float(atom['y']):16.10f} {float(atom['z']):16.10f}"
        for atom in atoms
    )
    Path(filename).write_text("\n".join(lines) + "\n", encoding="utf-8")


def create_output_directory(path: str | Path) -> Path:
    output = Path(path)
    output.mkdir(parents=True, exist_ok=False)
    return output


def atomic_write_text(path: Path, content: str) -> None:
    """Write a new text file through a same-directory partial file."""
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    partial.write_text(content, encoding="utf-8")
    partial.replace(path)


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract fragments and deterministic nearest-neighbor samples"
    )
    parser.add_argument("input_file", type=Path, help="Gaussian GJF/COM or XYZ input")
    parser.add_argument(
        "--output-dir", type=Path, default=Path("molecular_samples")
    )
    parser.add_argument("--samples", type=positive_int, default=20)
    parser.add_argument(
        "--layer", choices=["H", "M", "L", "all"], default="L"
    )
    parser.add_argument("--bond-scale", type=float, default=1.3)
    parser.add_argument("--min-fragment-atoms", type=positive_int, default=1)
    parser.add_argument(
        "--expected-fragments",
        type=positive_int,
        required=True,
        help="fail unless this independently expected fragment count remains",
    )
    parser.add_argument(
        "--xyz-units",
        choices=["angstrom"],
        help="required declaration for XYZ inputs because XYZ has no unit field",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    effective_argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(effective_argv)
    input_path = args.input_file.resolve()
    if not input_path.is_file():
        raise ValueError(f"Input file does not exist: {input_path}")
    input_bytes = input_path.read_bytes()
    try:
        input_text = input_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("Input file is not valid UTF-8 text") from exc
    input_sha256 = sha256_bytes(input_bytes)

    suffix = input_path.suffix.lower()
    if suffix in {".gjf", ".com"}:
        atoms = parse_gjf_text(input_text)
        input_format = "gaussian"
    elif suffix == ".xyz":
        if args.xyz_units != "angstrom":
            raise ValueError(
                "XYZ has no unit metadata; declare an angstrom input with "
                "--xyz-units angstrom"
            )
        atoms = parse_xyz_text(input_text)
        input_format = "xyz"
    else:
        raise ValueError("Unsupported input format; use .gjf, .com, or .xyz")

    all_fragments = identify_molecules(
        atoms,
        layer_filter=args.layer,
        bond_scale=args.bond_scale,
        min_fragment_atoms=1,
    )
    expected_source_indices = sorted(
        int(atom["source_index"])
        for atom in atoms
        if args.layer == "all" or atom["layer"] == args.layer
    )
    detected_source_indices = sorted(
        int(source_index)
        for molecule in all_fragments
        for source_index in molecule["source_indices"]
    )
    if detected_source_indices != expected_source_indices:
        raise RuntimeError("Detected fragments do not exactly cover selected atoms")
    molecules = [
        molecule
        for molecule in all_fragments
        if int(molecule["count"]) >= args.min_fragment_atoms
    ]
    for new_id, molecule in enumerate(molecules):
        molecule["id"] = new_id
    dropped = [
        molecule
        for molecule in all_fragments
        if int(molecule["count"]) < args.min_fragment_atoms
    ]
    if not molecules:
        raise ValueError("No fragments remain after applying the selected filters")
    if len(molecules) != args.expected_fragments:
        raise ValueError(
            "Fragment-count mismatch: "
            f"expected {args.expected_fragments}, detected {len(molecules)}"
        )

    samples = sample_molecules(molecules, n_samples=args.samples)
    output_dir = create_output_directory(args.output_dir)
    incomplete_marker = output_dir / "RUN_INCOMPLETE"
    incomplete_marker.write_text(
        "Artifact publication is incomplete until this marker is removed.\n",
        encoding="utf-8",
    )

    warnings: list[str] = []
    if dropped:
        warnings.append(
            f"{len(dropped)} fragment(s) were excluded by min_fragment_atoms"
        )

    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "internal_validation_passed",
        "scientific_status": "pending_independent_fragment_review",
        "runtime": {
            "python_version": platform.python_version(),
            "argv": effective_argv,
            "sampler_sha256": sha256_file(Path(__file__)),
        },
        "warnings": warnings,
        "input": {
            "name": input_path.name,
            "sha256": input_sha256,
            "hash_basis": "the exact bytes decoded and parsed for this run",
            "format": input_format,
            "units": "angstrom",
            "units_source": (
                "declared_by_cli" if input_format == "xyz" else "gaussian_default"
            ),
            "parsed_atom_count": len(atoms),
            "selected_atom_count": len(expected_source_indices),
            "layer": args.layer,
        },
        "parameters": {
            "bond_scale": args.bond_scale,
            "min_fragment_atoms": args.min_fragment_atoms,
            "maximum_unique_samples_per_size": args.samples,
            "expected_fragments": args.expected_fragments,
            "selection": "centroid-nearest-neighbors",
        },
        "fragment_summary": {
            "detected_before_size_filter": len(all_fragments),
            "retained_after_size_filter": len(molecules),
            "dropped_by_size_filter": len(dropped),
        },
        "fragments": [
            {
                "id": molecule["id"],
                "atom_count": molecule["count"],
                "elements": molecule["elements"],
                "source_indices": molecule["source_indices"],
            }
            for molecule in molecules
        ],
        "dropped_fragments": [
            {
                "atom_count": molecule["count"],
                "elements": molecule["elements"],
                "source_indices": molecule["source_indices"],
            }
            for molecule in dropped
        ],
        "samples": [],
        "validation": {
            "selected_source_index_partition": "passed",
            "generated_xyz_atom_counts": "passed",
            "member_sets_unique_within_each_size": "passed",
        },
    }

    sample_records: list[dict[str, object]] = []
    for sample_type, molecule_lists in samples.items():
        sample_dir = output_dir / sample_type
        sample_dir.mkdir()
        singular = sample_type[:-1]
        for index, molecule_list in enumerate(molecule_lists, start=1):
            member_ids = [int(molecule["id"]) for molecule in molecule_list]
            atom_count = sum(int(molecule["count"]) for molecule in molecule_list)
            filename = sample_dir / f"{singular}_{index:02d}.xyz"
            comment = (
                f"{singular} {index}; members={','.join(map(str, member_ids))}; "
                f"molecules={len(molecule_list)}; atoms={atom_count}"
            )
            write_xyz(molecule_list, filename, comment)
            reparsed = parse_xyz_file(filename)
            if len(reparsed) != atom_count:
                raise RuntimeError(f"Written XYZ failed atom-count check: {filename}")
            sample_records.append(
                {
                    "type": singular,
                    "member_ids": member_ids,
                    "atom_count": atom_count,
                    "path": filename.relative_to(output_dir).as_posix(),
                    "sha256": sha256_file(filename),
                }
            )

    manifest["samples"] = sample_records
    manifest_path = output_dir / "sampling_manifest.json"
    atomic_write_text(
        manifest_path,
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )

    summary_lines = [
        "=== Molecular Sampling Summary ===",
        "",
        f"Input name: {input_path.name}",
        f"Input SHA-256: {manifest['input']['sha256']}",
        f"Layer: {args.layer}",
        f"Parsed atoms: {len(atoms)}",
        f"Detected fragments before size filter: {len(all_fragments)}",
        f"Retained fragments: {len(molecules)}",
        f"Dropped fragments: {len(dropped)}",
        f"Status: {manifest['status']}",
        f"Scientific status: {manifest['scientific_status']}",
        "",
        "Sampling results:",
    ]
    summary_lines.extend(
        f"  - {name}: {len(values)} files" for name, values in samples.items()
    )
    if warnings:
        summary_lines.extend(["", "Warnings:"])
        summary_lines.extend(f"  - {warning}" for warning in warnings)
    atomic_write_text(
        output_dir / "sampling_summary.txt",
        "\n".join(summary_lines) + "\n",
    )
    incomplete_marker.unlink()

    print(
        f"Wrote {len(sample_records)} internally validated XYZ files to {output_dir}"
    )
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileExistsError, OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
