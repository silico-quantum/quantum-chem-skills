#!/usr/bin/env python3
"""Run a PySCF SCF calculation and generate explicitly indexed MO cubes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
from pathlib import Path
from typing import Any


ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$|^X$")
OCCUPATION_SUM_ABS_TOL = 1.0e-8


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    """Hash the exact in-memory source bytes captured at run start."""
    return hashlib.sha256(payload).hexdigest()


def publish_bytes_new(path: Path, payload: bytes) -> None:
    """Durably stage bytes and publish a new file without overwriting."""
    partial_path = path.with_name(f"{path.name}.partial")
    if path.exists() or partial_path.exists():
        raise FileExistsError(f"refusing to overwrite input snapshot: {path}")
    with partial_path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.link(partial_path, path)
    partial_path.unlink()


def flush_runtime_streams(*objects: Any) -> None:
    """Flush each distinct PySCF output stream before hashing the log."""
    flushed: set[int] = set()
    for obj in objects:
        stream = getattr(obj, "stdout", None)
        if stream is None or id(stream) in flushed:
            continue
        try:
            stream.flush()
        except Exception as exc:
            raise RuntimeError("failed to flush the PySCF log stream") from exc
        flushed.add(id(stream))


def validate_runtime_artifact(
    path: Path, description: str, run_started_ns: int
) -> dict[str, Any]:
    """Require a fresh, regular, non-empty artifact and record its hash."""
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"{description} is missing or is not a regular file: {path}")
    stat = path.stat()
    if stat.st_size <= 0:
        raise RuntimeError(f"{description} must be non-empty: {path}")
    if stat.st_mtime_ns < run_started_ns:
        raise RuntimeError(f"{description} is not fresh for this run: {path}")
    return {
        "status": "valid",
        "path": path.name,
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(path),
    }


def audit_runtime_artifacts(
    checkpoint_path: Path, log_path: Path, run_started_ns: int
) -> dict[str, dict[str, Any]]:
    """Audit both required artifacts so failure manifests remain complete."""
    artifacts: dict[str, dict[str, Any]] = {}
    for key, path, description in (
        ("checkpoint", checkpoint_path, "SCF checkpoint"),
        ("pyscf_log", log_path, "PySCF log"),
    ):
        try:
            artifacts[key] = validate_runtime_artifact(
                path, description, run_started_ns
            )
        except Exception as exc:
            artifacts[key] = {
                "status": "invalid",
                "path": path.name,
                "error": f"{type(exc).__name__}: {exc}",
            }
    return artifacts


def require_valid_runtime_artifacts(
    artifacts: dict[str, dict[str, Any]],
) -> None:
    """Fail acceptance when either audited runtime artifact is invalid."""
    errors = [
        record["error"]
        for record in artifacts.values()
        if record["status"] != "valid"
    ]
    if errors:
        raise RuntimeError("; ".join(errors))


def read_xyz(path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    """Read exactly one strict XYZ frame."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ must contain an atom-count line and a comment line")
    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("XYZ atom count is not an integer") from exc
    if atom_count <= 0:
        raise ValueError("XYZ atom count must be positive")
    if len(lines) != atom_count + 2:
        raise ValueError(
            f"XYZ declares {atom_count} atoms but contains {len(lines) - 2} records"
        )

    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line_number, line in enumerate(lines[2:], start=3):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"XYZ line {line_number} must have exactly four fields")
        symbol = fields[0]
        if not ELEMENT_RE.fullmatch(symbol):
            raise ValueError(f"invalid element symbol on XYZ line {line_number}: {symbol}")
        try:
            xyz = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise ValueError(f"non-numeric coordinate on XYZ line {line_number}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"non-finite coordinate on XYZ line {line_number}")
        atoms.append((symbol, xyz))
    return atoms


def parse_orbital_request(token: str) -> tuple[str, int]:
    """Parse CHANNEL:INDEX, where INDEX is deliberately one-based."""
    try:
        channel, raw_index = token.lower().split(":", maxsplit=1)
        index = int(raw_index)
    except (ValueError, AttributeError) as exc:
        raise argparse.ArgumentTypeError(
            "orbital requests must use CHANNEL:INDEX, for example spatial:5"
        ) from exc
    if channel not in {"spatial", "alpha", "beta"}:
        raise argparse.ArgumentTypeError("channel must be spatial, alpha, or beta")
    if index <= 0:
        raise argparse.ArgumentTypeError("orbital INDEX is one-based and must be positive")
    return channel, index


def parse_finite_float(token: str, description: str) -> float:
    """Parse a finite decimal, including Fortran D exponents."""
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(f"{description} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{description} is non-finite: {token!r}")
    return value


def validate_cube_file(path: Path, expected_atoms: int) -> dict[str, Any]:
    """Validate the complete scalar payload of a single-dataset cube file."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"cube is absent or empty: {path}")
    with path.open(encoding="utf-8") as handle:
        comments = [handle.readline(), handle.readline()]
        if any(line == "" for line in comments):
            raise ValueError(f"cube comments are truncated: {path}")

        origin_line = handle.readline()
        if origin_line == "":
            raise ValueError(f"cube origin line is missing: {path}")
        origin_fields = origin_line.split()
        if len(origin_fields) not in {4, 5}:
            raise ValueError(f"cube origin line is malformed: {path}")
        try:
            atom_count = int(origin_fields[0])
        except ValueError as exc:
            raise ValueError(f"cube atom count is malformed: {path}") from exc
        if atom_count <= 0:
            raise ValueError(
                "cube must use a positive atom count and one scalar dataset"
            )
        for axis, token in enumerate(origin_fields[1:4], start=1):
            parse_finite_float(token, f"cube origin coordinate {axis}")
        data_set_count = 1
        if len(origin_fields) == 5:
            try:
                data_set_count = int(origin_fields[4])
            except ValueError as exc:
                raise ValueError(f"cube NVAL is malformed: {path}") from exc
            if data_set_count != 1:
                raise ValueError("orbital cube must contain exactly one scalar dataset")

        grid: list[int] = []
        for axis in range(1, 4):
            axis_line = handle.readline()
            if axis_line == "":
                raise ValueError(f"cube grid axis {axis} is missing: {path}")
            fields = axis_line.split()
            if len(fields) != 4:
                raise ValueError(f"cube grid axis {axis} is malformed: {path}")
            try:
                point_count = int(fields[0])
            except ValueError as exc:
                raise ValueError(
                    f"cube grid axis {axis} point count is malformed: {path}"
                ) from exc
            if point_count <= 0:
                raise ValueError(
                    f"cube grid axis {axis} point count must be positive"
                )
            vector = [
                parse_finite_float(token, f"cube grid axis {axis} vector")
                for token in fields[1:]
            ]
            if math.sqrt(sum(value * value for value in vector)) == 0.0:
                raise ValueError(f"cube grid axis {axis} vector must be non-zero")
            grid.append(point_count)

        if atom_count != expected_atoms:
            raise ValueError(
                f"cube atom count {atom_count} does not match molecule {expected_atoms}"
            )

        for atom_index in range(1, atom_count + 1):
            atom_line = handle.readline()
            if atom_line == "":
                raise ValueError(f"cube atom record {atom_index} is missing: {path}")
            fields = atom_line.split()
            if len(fields) != 5:
                raise ValueError(
                    f"cube atom record {atom_index} must contain five fields"
                )
            try:
                atomic_number = int(fields[0])
            except ValueError as exc:
                raise ValueError(
                    f"cube atom record {atom_index} has an invalid atomic number"
                ) from exc
            if atomic_number < 0:
                raise ValueError(
                    f"cube atom record {atom_index} has a negative atomic number"
                )
            for field_index, token in enumerate(fields[1:], start=2):
                parse_finite_float(
                    token, f"cube atom record {atom_index} field {field_index}"
                )

        payload_tokens = handle.read().split()

    expected_values = math.prod(grid) * data_set_count
    if len(payload_tokens) != expected_values:
        raise ValueError(
            f"cube payload has {len(payload_tokens)} values; expected {expected_values}"
        )
    payload = [
        parse_finite_float(token, f"cube payload value {index}")
        for index, token in enumerate(payload_tokens, start=1)
    ]
    return {
        "bytes": path.stat().st_size,
        "atom_count": atom_count,
        "grid": grid,
        "data_set_count": data_set_count,
        "value_count": len(payload),
        "value_min": min(payload),
        "value_max": max(payload),
    }


def validate_cube_header(path: Path, expected_atoms: int) -> dict[str, Any]:
    """Backward-compatible name for complete cube validation."""
    return validate_cube_file(path, expected_atoms)


def array_shape(array: Any, description: str) -> tuple[int, ...]:
    """Return a positive integer shape without importing NumPy in this module."""
    raw_shape = getattr(array, "shape", None)
    if raw_shape is None:
        raise ValueError(f"{description} has no shape")
    try:
        shape = tuple(int(value) for value in raw_shape)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{description} has an invalid shape") from exc
    if not shape or any(value <= 0 for value in shape):
        raise ValueError(f"{description} has a non-positive shape: {shape}")
    return shape


def validate_finite_array(
    array: Any, description: str, expected_count: int
) -> list[float]:
    """Return exactly ``expected_count`` real finite array entries."""
    flat = getattr(array, "flat", None)
    if flat is None:
        raise ValueError(f"{description} cannot be flattened")
    values: list[float] = []
    for count, raw_value in enumerate(flat, start=1):
        try:
            value = float(raw_value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{description} value {count} is not real numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"{description} value {count} is non-finite")
        values.append(value)
    if len(values) != expected_count:
        raise ValueError(
            f"{description} stores {len(values)} values but its shape requires "
            f"{expected_count}"
        )
    return values


def validate_electronic_result(
    total_energy: float,
    channels: dict[str, tuple[Any, Any, Any]],
    electron_count: Any,
) -> dict[str, dict[str, int]]:
    """Validate total energy and all MO coefficient/energy/occupation arrays."""
    try:
        energy = float(total_energy)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("SCF total energy is not real numeric") from exc
    if not math.isfinite(energy):
        raise ValueError("SCF total energy is non-finite")

    if isinstance(electron_count, bool) or type(electron_count).__name__ == "bool_":
        raise ValueError("molecular electron count must be a non-negative integer")
    try:
        expected_electrons = int(electron_count)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "molecular electron count must be a non-negative integer"
        ) from exc
    if expected_electrons < 0 or electron_count != expected_electrons:
        raise ValueError("molecular electron count must be a non-negative integer")

    channel_names = set(channels)
    if channel_names not in ({"spatial"}, {"alpha", "beta"}):
        raise ValueError(
            "orbital channels must be spatial or the complete alpha/beta pair"
        )

    summary: dict[str, dict[str, int]] = {}
    ao_counts: set[int] = set()
    orbital_counts: set[int] = set()
    occupation_values: dict[str, list[float]] = {}
    for channel, (coefficients, energies, occupations) in channels.items():
        coefficient_shape = array_shape(coefficients, f"{channel} coefficients")
        energy_shape = array_shape(energies, f"{channel} energies")
        occupation_shape = array_shape(occupations, f"{channel} occupations")
        if len(coefficient_shape) != 2:
            raise ValueError(
                f"{channel} coefficients shape must be two-dimensional"
            )
        ao_count, orbital_count = coefficient_shape
        if energy_shape != (orbital_count,) or occupation_shape != (orbital_count,):
            raise ValueError(
                f"{channel} orbital array shape mismatch: coefficients "
                f"{coefficient_shape}, energies {energy_shape}, occupations "
                f"{occupation_shape}"
            )
        validate_finite_array(
            coefficients, f"{channel} coefficients", ao_count * orbital_count
        )
        validate_finite_array(energies, f"{channel} energies", orbital_count)
        values = validate_finite_array(
            occupations, f"{channel} occupations", orbital_count
        )
        upper_bound = 2.0 if channel == "spatial" else 1.0
        if any(value < 0.0 or value > upper_bound for value in values):
            raise ValueError(
                f"{channel} occupations must lie in the physical range "
                f"[0, {upper_bound:g}]"
            )
        occupation_values[channel] = values
        ao_counts.add(ao_count)
        orbital_counts.add(orbital_count)
        summary[channel] = {
            "ao_count": ao_count,
            "orbital_count": orbital_count,
        }
    if len(ao_counts) != 1 or len(orbital_counts) != 1:
        raise ValueError("alpha and beta orbital array shapes do not match")
    observed_electrons = sum(
        value for values in occupation_values.values() for value in values
    )
    if not math.isclose(
        observed_electrons,
        float(expected_electrons),
        rel_tol=0.0,
        abs_tol=OCCUPATION_SUM_ABS_TOL,
    ):
        raise ValueError(
            "orbital occupation sum does not match the molecular electron count "
            f"within {OCCUPATION_SUM_ABS_TOL:.1e} electrons"
        )
    return summary


def build_mean_field(mol: Any, method: str, xc: str | None) -> Any:
    from pyscf import dft, scf

    if method == "rhf":
        if mol.spin != 0:
            raise ValueError("RHF requires spin=0; use ROHF or UHF for an open shell")
        return scf.RHF(mol)
    if method == "rohf":
        if mol.spin == 0:
            raise ValueError("ROHF is reserved here for spin>0; use RHF for spin=0")
        return scf.ROHF(mol)
    if method == "uhf":
        return scf.UHF(mol)
    if method == "rks":
        if mol.spin != 0:
            raise ValueError("RKS requires spin=0; use UKS for an open shell")
        mean_field = dft.RKS(mol)
        mean_field.xc = xc
        return mean_field
    if method == "uks":
        mean_field = dft.UKS(mol)
        mean_field.xc = xc
        return mean_field
    raise ValueError(f"unsupported method: {method}")


def orbital_channels(mean_field: Any) -> dict[str, tuple[Any, Any, Any]]:
    """Return channel -> (coefficients, energies, occupations)."""
    coeff = mean_field.mo_coeff
    if getattr(coeff, "ndim", 0) == 2:
        return {"spatial": (coeff, mean_field.mo_energy, mean_field.mo_occ)}
    if len(coeff) != 2:
        raise ValueError("unexpected unrestricted MO coefficient layout")
    return {
        "alpha": (coeff[0], mean_field.mo_energy[0], mean_field.mo_occ[0]),
        "beta": (coeff[1], mean_field.mo_energy[1], mean_field.mo_occ[1]),
    }


def write_orbital_table(
    path: Path, channels: dict[str, tuple[Any, Any, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["channel", "orbital_index_1based", "energy_hartree", "occupation"])
        for channel, (_, energies, occupations) in channels.items():
            for index, (energy, occupation) in enumerate(
                zip(energies, occupations), start=1
            ):
                writer.writerow(
                    [channel, index, f"{float(energy):.12g}", f"{float(occupation):.8g}"]
                )


def preflight_orbital_requests(
    requests: list[tuple[str, int]],
    channels: dict[str, tuple[Any, Any, Any]],
) -> list[dict[str, Any]]:
    """Resolve every request before any cube generator can create a file."""
    plan: list[dict[str, Any]] = []
    seen: set[tuple[str, int]] = set()
    for channel, one_based_index in requests:
        request = (channel, one_based_index)
        if request in seen:
            raise ValueError(f"duplicate orbital request: {channel}:{one_based_index}")
        seen.add(request)
        if channel not in channels:
            available = ", ".join(sorted(channels))
            raise ValueError(
                f"channel {channel!r} is invalid for this calculation; use {available}"
            )
        if one_based_index <= 0:
            raise IndexError("orbital indices are one-based and must be positive")
        coefficients, energies, occupations = channels[channel]
        coefficient_shape = array_shape(coefficients, f"{channel} coefficients")
        if len(coefficient_shape) != 2:
            raise ValueError(
                f"{channel} coefficients shape must be two-dimensional"
            )
        zero_based_index = one_based_index - 1
        if zero_based_index >= coefficient_shape[1]:
            raise IndexError(
                f"orbital {channel}:{one_based_index} exceeds "
                f"{coefficient_shape[1]} MOs"
            )
        try:
            coefficient_column = coefficients[:, zero_based_index]
            orbital_energy = float(energies[zero_based_index])
            occupation = float(occupations[zero_based_index])
        except (IndexError, TypeError, ValueError) as exc:
            raise ValueError(
                f"failed to resolve orbital {channel}:{one_based_index}"
            ) from exc
        if not math.isfinite(orbital_energy) or not math.isfinite(occupation):
            raise ValueError(
                f"orbital {channel}:{one_based_index} metadata is non-finite"
            )
        plan.append(
            {
                "channel": channel,
                "orbital_index_1based": one_based_index,
                "coefficients": coefficient_column,
                "energy_hartree": orbital_energy,
                "occupation": occupation,
            }
        )
    return plan


def generate_cube_artifacts(
    mol: Any,
    cubegen: Any,
    channels: dict[str, tuple[Any, Any, Any]],
    requests: list[tuple[str, int]],
    output_dir: Path,
    *,
    grid_points: int,
    margin_bohr: float,
) -> list[dict[str, Any]]:
    """Stage, completely validate, and atomically publish requested cubes."""
    plan = preflight_orbital_requests(requests, channels)
    staged: list[tuple[Path, Path, dict[str, Any]]] = []
    for item in plan:
        channel = item["channel"]
        one_based_index = item["orbital_index_1based"]
        final_path = output_dir / f"mo_{channel}_{one_based_index:04d}.cube"
        partial_path = output_dir / f"{final_path.name}.partial"
        if final_path.exists() or partial_path.exists():
            raise FileExistsError(
                f"refusing to overwrite cube artifact for {channel}:{one_based_index}"
            )
        cubegen.orbital(
            mol,
            str(partial_path),
            item["coefficients"],
            nx=grid_points,
            ny=grid_points,
            nz=grid_points,
            margin=margin_bohr,
        )
        cube_record = validate_cube_file(partial_path, mol.natm)
        cube_record.update(
            {
                "path": final_path.name,
                "sha256": sha256_file(partial_path),
                "channel": channel,
                "orbital_index_1based": one_based_index,
                "energy_hartree": item["energy_hartree"],
                "occupation": item["occupation"],
            }
        )
        staged.append((partial_path, final_path, cube_record))

    published: list[Path] = []
    try:
        for partial_path, final_path, _ in staged:
            os.link(partial_path, final_path)
            published.append(final_path)
            partial_path.unlink()
    except Exception:
        for published_path in published:
            try:
                published_path.unlink()
            except FileNotFoundError:
                pass
        raise
    return [record for _, _, record in staged]


def inventory_artifacts(output_dir: Path) -> list[dict[str, Any]]:
    """Best-effort inventory of diagnostics and rejected files after failure."""
    inventory: list[dict[str, Any]] = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file() or path.name == "run.json":
            continue
        if path.name.endswith(".partial"):
            status = "rejected_partial"
        elif path.suffix == ".cube":
            status = "rejected_unaccepted"
        else:
            status = "diagnostic"
        record: dict[str, Any] = {
            "path": str(path.relative_to(output_dir)),
            "status": status,
        }
        try:
            record["bytes"] = path.stat().st_size
            record["sha256"] = sha256_file(path)
        except OSError as exc:
            record["inventory_error"] = f"{type(exc).__name__}: {exc}"
            record["bytes"] = None
            record["sha256"] = None
        inventory.append(record)
    return inventory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xyz", type=Path, help="single-frame XYZ geometry")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--unit", choices=["Angstrom", "Bohr"], required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument(
        "--spin",
        type=int,
        required=True,
        help="PySCF spin = Nalpha - Nbeta = 2S, not multiplicity",
    )
    parser.add_argument("--basis", required=True)
    parser.add_argument(
        "--method", choices=["rhf", "rohf", "uhf", "rks", "uks"], required=True
    )
    parser.add_argument(
        "--xc",
        default=None,
        help="explicit XC functional required for RKS/UKS",
    )
    parser.add_argument("--max-cycle", type=int, default=100)
    parser.add_argument("--grid-points", type=int, default=80)
    parser.add_argument("--margin-bohr", type=float, default=3.0)
    parser.add_argument(
        "--orbital",
        action="append",
        default=[],
        type=parse_orbital_request,
        metavar="CHANNEL:INDEX",
        help="one-based explicit index; repeat to request multiple cubes",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.spin < 0:
        raise ValueError("spin must be a non-negative Nalpha-Nbeta value")
    if args.method in {"rks", "uks"} and (
        not isinstance(args.xc, str) or not args.xc.strip()
    ):
        raise ValueError(f"{args.method.upper()} requires an explicit --xc value")
    if args.max_cycle <= 0 or args.grid_points <= 0 or args.margin_bohr <= 0:
        raise ValueError("max-cycle, grid-points, and margin-bohr must be positive")
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to reuse output directory: {args.output_dir}")
    if not args.output_dir.parent.is_dir():
        raise FileNotFoundError(f"output parent does not exist: {args.output_dir.parent}")

    source_bytes = args.xyz.read_bytes()
    source_sha256 = sha256_bytes(source_bytes)
    run_started_ns = time.time_ns()
    args.output_dir.mkdir()
    manifest_path = args.output_dir / "run.json"
    snapshot_path = args.output_dir / "input.xyz"
    log_path = args.output_dir / "pyscf.log"
    checkpoint_path = args.output_dir / "scf.chk"
    phase = "publish_input_snapshot"
    snapshot_sha256: str | None = None
    pyscf_version: str | None = None
    energy: float | None = None
    converged: bool | None = None
    channel_summary: dict[str, dict[str, int]] | None = None
    runtime_artifacts: dict[str, dict[str, Any]] | None = None
    mol: Any = None
    mean_field: Any = None

    try:
        publish_bytes_new(snapshot_path, source_bytes)
        snapshot_sha256 = sha256_file(snapshot_path)
        if snapshot_sha256 != source_sha256:
            raise RuntimeError("published input snapshot hash does not match source bytes")

        phase = "parse_input_snapshot"
        atoms = read_xyz(snapshot_path)

        phase = "import_dependencies"
        import pyscf
        from pyscf import gto
        from pyscf.tools import cubegen

        pyscf_version = pyscf.__version__
        phase = "build_molecule"
        mol = gto.M(
            atom=atoms,
            basis=args.basis,
            charge=args.charge,
            spin=args.spin,
            unit=args.unit,
            verbose=4,
            output=str(log_path),
        )
        if mol.nelectron < args.spin or (mol.nelectron - args.spin) % 2:
            raise ValueError(
                "charge and spin are inconsistent with the molecular electron count"
            )

        phase = "run_scf"
        mean_field = build_mean_field(mol, args.method, args.xc)
        mean_field.max_cycle = args.max_cycle
        mean_field.chkfile = str(checkpoint_path)
        energy = float(mean_field.kernel())
        converged = bool(mean_field.converged)
        if not converged:
            raise RuntimeError("SCF did not converge; no orbital cube was accepted")

        flush_runtime_streams(mol, mean_field)
        phase = "validate_runtime_artifacts"
        runtime_artifacts = audit_runtime_artifacts(
            checkpoint_path, log_path, run_started_ns
        )
        require_valid_runtime_artifacts(runtime_artifacts)

        phase = "validate_electronic_result"
        channels = orbital_channels(mean_field)
        channel_summary = validate_electronic_result(
            energy, channels, electron_count=mol.nelectron
        )
        phase = "write_orbital_table"
        table_path = args.output_dir / "orbitals.csv"
        write_orbital_table(table_path, channels)

        phase = "generate_cube_artifacts"
        cubes = generate_cube_artifacts(
            mol,
            cubegen,
            channels,
            args.orbital,
            args.output_dir,
            grid_points=args.grid_points,
            margin_bohr=args.margin_bohr,
        )

        flush_runtime_streams(mol, mean_field)
        phase = "finalize_runtime_artifacts"
        runtime_artifacts = audit_runtime_artifacts(
            checkpoint_path, log_path, run_started_ns
        )
        require_valid_runtime_artifacts(runtime_artifacts)

        phase = "write_accepted_manifest"
        manifest = {
            "status": "accepted",
            "run_started_ns": run_started_ns,
            "input": {
                "path": snapshot_path.name,
                "source_path": str(args.xyz),
                "sha256": snapshot_sha256,
                "bytes": len(source_bytes),
                "atom_count": mol.natm,
                "unit": args.unit,
                "charge": args.charge,
                "spin_nalpha_minus_nbeta": args.spin,
                "electron_count": mol.nelectron,
            },
            "calculation": {
                "python_version": platform.python_version(),
                "python_executable": sys.executable,
                "pyscf_version": pyscf_version,
                "method": args.method,
                "xc": args.xc if args.method in {"rks", "uks"} else None,
                "basis": args.basis,
                "converged": True,
                "energy_hartree": energy,
                "max_cycle": args.max_cycle,
                "orbital_array_shapes": channel_summary,
            },
            "runtime_artifacts": runtime_artifacts,
            "cube_grid": {
                "points_per_axis": args.grid_points,
                "margin_bohr": args.margin_bohr,
            },
            "orbital_table": {
                "path": table_path.name,
                "sha256": sha256_file(table_path),
            },
            "cubes": cubes,
        }
        manifest_path.write_text(
            json.dumps(manifest, allow_nan=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except Exception as exc:
        failure_flush_error: str | None = None
        try:
            flush_runtime_streams(mol, mean_field)
        except Exception as flush_exc:
            failure_flush_error = f"{type(flush_exc).__name__}: {flush_exc}"
        runtime_artifacts = audit_runtime_artifacts(
            checkpoint_path, log_path, run_started_ns
        )
        finite_energy = (
            energy if energy is not None and math.isfinite(energy) else None
        )
        failure_snapshot_sha256 = snapshot_sha256
        if failure_snapshot_sha256 is None and snapshot_path.is_file():
            try:
                failure_snapshot_sha256 = sha256_file(snapshot_path)
            except OSError:
                pass
        failure = {
            "status": "failed",
            "run_started_ns": run_started_ns,
            "failure_phase": phase,
            "error_type": type(exc).__name__,
            "error": f"{type(exc).__name__}: {exc}",
            "failure_flush_error": failure_flush_error,
            "input": {
                "path": snapshot_path.name if snapshot_path.is_file() else None,
                "source_path": str(args.xyz),
                "source_sha256_at_start": source_sha256,
                "sha256": failure_snapshot_sha256,
                "bytes": len(source_bytes),
                "unit": args.unit,
                "charge": args.charge,
                "spin_nalpha_minus_nbeta": args.spin,
            },
            "calculation": {
                "python_version": platform.python_version(),
                "python_executable": sys.executable,
                "pyscf_version": pyscf_version,
                "method": args.method,
                "xc": args.xc if args.method in {"rks", "uks"} else None,
                "basis": args.basis,
                "max_cycle": args.max_cycle,
                "converged": converged,
                "energy_hartree": finite_energy,
                "energy_observed": repr(energy) if energy is not None else None,
                "orbital_array_shapes": channel_summary,
            },
            "runtime_artifacts": runtime_artifacts,
            "cube_grid": {
                "points_per_axis": args.grid_points,
                "margin_bohr": args.margin_bohr,
            },
            "requested_orbitals": [
                {"channel": channel, "orbital_index_1based": index}
                for channel, index in args.orbital
            ],
            "artifacts": inventory_artifacts(args.output_dir),
        }
        manifest_path.write_text(
            json.dumps(failure, allow_nan=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        raise

    print(f"Accepted output: {args.output_dir}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
