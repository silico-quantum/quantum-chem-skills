#!/usr/bin/env python3
"""Run one fail-closed closed-shell RKS/TDA calculation.

PySCF and NumPy are imported only after the input contract is validated and a
fresh run directory contains a ``running`` manifest.  This lets ``--help`` and
static contract tests run in environments where PySCF is not installed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CONFIG_KEYS = (
    "atom",
    "unit",
    "basis",
    "charge",
    "spin",
    "xc",
    "grid_level",
    "conv_tol",
    "max_cycle",
    "nstates",
)

SCF_ENERGY_ABS_TOL_HARTREE = 1.0e-10
OCCUPATION_SUM_ABS_TOL = 1.0e-8


def utc_now() -> str:
    """Return a second-resolution UTC timestamp."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_bytes(data: bytes) -> str:
    """Return the lowercase SHA-256 digest of *data*."""

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    """Hash a regular, non-empty artifact."""

    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"required artifact is missing or empty: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def flush_pyscf_stdout(mol: Any, mf: Any) -> None:
    """Flush the molecule and mean-field log streams, rejecting any failure."""

    flushed: set[int] = set()
    for name, obj in (("mol", mol), ("mf", mf)):
        stream = getattr(obj, "stdout", None)
        flush = getattr(stream, "flush", None)
        if not callable(flush):
            raise RuntimeError(f"PySCF {name}.stdout is not flushable")
        if id(stream) in flushed:
            continue
        try:
            flush()
        except Exception as exc:
            raise RuntimeError(f"failed to flush PySCF {name}.stdout") from exc
        flushed.add(id(stream))


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    """Write JSON through a sibling temporary file and atomically replace."""

    temporary = path.with_name(f".{path.name}.tmp")
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=True,
        indent=2,
        sort_keys=True,
    )
    with temporary.open("x", encoding="utf-8") as handle:
        handle.write(encoded)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _require_int(value: Any, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _require_finite_float(
    value: Any, name: str, *, positive: bool = False
) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be finite")
    if positive and converted <= 0:
        raise ValueError(f"{name} must be positive")
    return converted


def validate_config(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the intentionally narrow RKS/TDA input contract."""

    if not isinstance(raw, Mapping):
        raise ValueError("the configuration root must be a JSON object")
    missing = sorted(set(CONFIG_KEYS) - set(raw))
    unknown = sorted(set(raw) - set(CONFIG_KEYS))
    if missing:
        raise ValueError(f"missing configuration fields: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"unknown configuration fields: {', '.join(unknown)}")

    atom = raw["atom"]
    if not isinstance(atom, str) or not atom.strip():
        raise ValueError("atom must be a non-empty PySCF atom string")
    unit = raw["unit"]
    if unit not in {"Angstrom", "Bohr"}:
        raise ValueError("unit must be exactly 'Angstrom' or 'Bohr'")
    basis = raw["basis"]
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError("basis must be a non-empty string")
    xc = raw["xc"]
    if not isinstance(xc, str) or not xc.strip():
        raise ValueError("xc must be a non-empty string")

    charge = _require_int(raw["charge"], "charge")
    spin = _require_int(raw["spin"], "spin", minimum=0)
    if spin != 0:
        raise ValueError(
            "this supported runner is closed-shell RKS/TDA only; use spin=0"
        )
    grid_level = _require_int(raw["grid_level"], "grid_level", minimum=0)
    if grid_level > 9:
        raise ValueError("grid_level must be between 0 and 9")
    max_cycle = _require_int(raw["max_cycle"], "max_cycle", minimum=1)
    nstates = _require_int(raw["nstates"], "nstates", minimum=1)
    conv_tol = _require_finite_float(raw["conv_tol"], "conv_tol", positive=True)

    return {
        "atom": atom,
        "unit": unit,
        "basis": basis,
        "charge": charge,
        "spin": spin,
        "xc": xc,
        "grid_level": grid_level,
        "conv_tol": conv_tol,
        "max_cycle": max_cycle,
        "nstates": nstates,
    }


def _iter_numeric_scalars(value: Any, name: str) -> Iterable[complex]:
    """Yield numeric scalars from nested Python or NumPy sequences."""

    if isinstance(value, (str, bytes, bytearray)):
        raise ValueError(f"{name} contains non-numeric data")
    try:
        converted = complex(value)
    except (TypeError, ValueError, OverflowError):
        try:
            iterator = iter(value)
        except TypeError as exc:
            raise ValueError(f"{name} contains non-numeric data") from exc
        found = False
        for child in iterator:
            found = True
            yield from _iter_numeric_scalars(child, name)
        if not found:
            raise ValueError(f"{name} is empty")
    else:
        yield converted


def _require_finite_array(value: Any, name: str) -> None:
    found = False
    for scalar in _iter_numeric_scalars(value, name):
        found = True
        if not (math.isfinite(scalar.real) and math.isfinite(scalar.imag)):
            raise ValueError(f"{name} contains a non-finite value")
    if not found:
        raise ValueError(f"{name} is empty")


def _array_shape(value: Any, name: str) -> tuple[int, ...]:
    """Return a positive rectangular shape without importing NumPy."""

    raw_shape = getattr(value, "shape", None)
    if raw_shape is not None:
        try:
            shape = tuple(int(dimension) for dimension in raw_shape)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} has an invalid shape") from exc
    else:
        def infer(sequence: Any) -> tuple[int, ...]:
            if not isinstance(sequence, (list, tuple)):
                return ()
            if not sequence:
                return (0,)
            child_shapes = [infer(child) for child in sequence]
            if any(child != child_shapes[0] for child in child_shapes[1:]):
                raise ValueError(f"{name} has a ragged shape")
            return (len(sequence), *child_shapes[0])

        shape = infer(value)
    if not shape or any(dimension <= 0 for dimension in shape):
        raise ValueError(f"{name} has a non-positive shape: {shape}")
    return shape


def _require_finite_real_vector(value: Any, name: str) -> list[float]:
    """Return a non-empty flat vector of finite real scalar values."""

    result: list[float] = []
    for scalar in _iter_numeric_scalars(value, name):
        if not (math.isfinite(scalar.real) and math.isfinite(scalar.imag)):
            raise ValueError(f"{name} contains a non-finite value")
        if scalar.imag != 0.0:
            raise ValueError(f"{name} contains a non-real value")
        result.append(float(scalar.real))
    if not result:
        raise ValueError(f"{name} is empty")
    return result


def validate_scf_result(mf: Any, energy: Any, checkpoint: Path) -> None:
    """Reject incomplete, unconverged, or non-finite SCF results."""

    if not bool(getattr(mf, "converged", False)):
        raise ValueError("SCF did not converge")
    kernel_energy = _require_finite_float(energy, "SCF kernel energy")
    final_energy = _require_finite_float(getattr(mf, "e_tot", None), "mf.e_tot")
    if abs(kernel_energy - final_energy) > SCF_ENERGY_ABS_TOL_HARTREE:
        raise ValueError(
            "SCF kernel energy does not match mf.e_tot within "
            f"{SCF_ENERGY_ABS_TOL_HARTREE:.1e} Hartree"
        )
    mo_energy = getattr(mf, "mo_energy", None)
    mo_coeff = getattr(mf, "mo_coeff", None)
    mo_occ = getattr(mf, "mo_occ", None)
    energy_shape = _array_shape(mo_energy, "mf.mo_energy")
    occupation_shape = _array_shape(mo_occ, "mf.mo_occ")
    coefficient_shape = _array_shape(mo_coeff, "mf.mo_coeff")
    if len(energy_shape) != 1:
        raise ValueError(f"mf.mo_energy shape must be one-dimensional: {energy_shape}")
    orbital_count = energy_shape[0]
    if occupation_shape != (orbital_count,):
        raise ValueError(
            "mf.mo_occ shape must match mf.mo_energy: "
            f"{occupation_shape} != {(orbital_count,)}"
        )
    if len(coefficient_shape) != 2 or coefficient_shape[1] != orbital_count:
        raise ValueError(
            "mf.mo_coeff shape must be (positive AO count, orbital count): "
            f"{coefficient_shape}, expected second dimension {orbital_count}"
        )
    _require_finite_array(mo_energy, "mf.mo_energy")
    _require_finite_array(mo_coeff, "mf.mo_coeff")
    occupations = _require_finite_real_vector(
        mo_occ, "mf.mo_occ"
    )
    if any(value < 0.0 or value > 2.0 for value in occupations):
        raise ValueError("mf.mo_occ values must be in the physical range [0, 2]")
    molecule = getattr(mf, "mol", None)
    electron_count = _require_int(
        getattr(molecule, "nelectron", None), "mf.mol.nelectron", minimum=1
    )
    if not math.isclose(
        sum(occupations),
        float(electron_count),
        rel_tol=0.0,
        abs_tol=OCCUPATION_SUM_ABS_TOL,
    ):
        raise ValueError(
            "mf.mo_occ sum does not match the molecular electron count within "
            f"{OCCUPATION_SUM_ABS_TOL:.1e} electrons"
        )
    sha256_file(checkpoint)


def validate_scf_stability(mf: Any) -> dict[str, bool]:
    """Require successful internal and external SCF stability analyses."""

    stability = getattr(mf, "stability", None)
    if not callable(stability):
        raise ValueError("SCF stability analysis is unavailable")
    result = stability(internal=True, external=True, return_status=True)
    if isinstance(result, (str, bytes, bytearray)):
        raise ValueError("SCF stability analysis returned an invalid result")
    try:
        _, _, internal, external = result
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "SCF stability analysis did not return internal and external status"
        ) from exc
    internal_status = _as_bool(internal, "internal SCF stability")
    external_status = _as_bool(external, "external SCF stability")
    if not internal_status or not external_status:
        failed = []
        if not internal_status:
            failed.append("internal")
        if not external_status:
            failed.append("external")
        raise ValueError(f"SCF stability gate failed: {', '.join(failed)}")
    return {"internal": internal_status, "external": external_status}


def _as_scalar_vector(values: Any, name: str, expected: int) -> list[Any]:
    if isinstance(values, (str, bytes, bytearray)):
        raise ValueError(f"{name} must be a one-dimensional vector")
    try:
        vector = list(values)
    except TypeError as exc:
        raise ValueError(f"{name} must be a one-dimensional vector") from exc
    if len(vector) != expected:
        raise ValueError(f"{name} has {len(vector)} roots; expected {expected}")
    for value in vector:
        if isinstance(value, (str, bytes, bytearray)):
            raise ValueError(f"{name} must contain scalar values")
        try:
            complex(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} must contain scalar values") from exc
    return vector


def _as_bool(value: Any, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if type(value).__name__ == "bool_":
        return bool(value)
    raise ValueError(f"{name} must contain Boolean convergence flags")


def validate_tda_result(
    converged: Any,
    energies: Any,
    oscillator_strengths: Any,
    nstates: int,
) -> None:
    """Require exactly ``nstates`` converged, physical finite TDA roots."""

    requested = _require_int(nstates, "nstates", minimum=1)
    convergence_vector = _as_scalar_vector(
        converged, "td.converged", requested
    )
    energy_vector = _as_scalar_vector(energies, "td.e", requested)
    strength_vector = _as_scalar_vector(
        oscillator_strengths, "oscillator strengths", requested
    )
    if not all(_as_bool(value, "td.converged") for value in convergence_vector):
        raise ValueError("at least one requested TDA root did not converge")
    energy_values = [
        _require_finite_float(value, "td.e") for value in energy_vector
    ]
    strength_values = [
        _require_finite_float(value, "oscillator strength")
        for value in strength_vector
    ]
    if any(value <= 0.0 for value in energy_values):
        raise ValueError("every TDA excitation energy must be positive in Hartree")
    if any(value < 0.0 for value in strength_values):
        raise ValueError("every oscillator strength must be non-negative")


def initialize_run(
    output_dir: Path, config_bytes: bytes, config_sha256: str
) -> dict[str, Any]:
    """Create a fresh output directory and its running manifest."""

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=False)
    input_copy = output_dir / "input_config.json"
    with input_copy.open("xb") as handle:
        handle.write(config_bytes)
        handle.flush()
        os.fsync(handle.fileno())

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "state": "running",
        "started_at_utc": utc_now(),
        "python_version": platform.python_version(),
        "pyscf_version": "not_loaded",
        "runner": "scripts/run_safe_dft_tda.py",
        "runner_sha256": sha256_file(Path(__file__)),
        "input_config_sha256": config_sha256,
    }
    _write_json_atomic(output_dir / "run_manifest.partial.json", manifest)
    return manifest


def _artifact_record(path: Path, output_dir: Path) -> dict[str, Any]:
    resolved_output = output_dir.resolve()
    resolved_path = path.resolve()
    try:
        relative = resolved_path.relative_to(resolved_output)
    except ValueError as exc:
        raise ValueError(f"artifact is outside the run directory: {path}") from exc
    return {
        "path": relative.as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def finalize_accepted_run(
    *,
    output_dir: Path,
    manifest: Mapping[str, Any],
    results: Mapping[str, Any],
    artifacts: Iterable[Path],
    log_path: Path,
    mol: Any,
    mf: Any,
    pyscf_version: str,
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Flush logs, publish acceptance, then prove the log hash stayed stable."""

    output_dir = Path(output_dir)
    log_path = Path(log_path)
    artifact_paths = tuple(Path(artifact) for artifact in artifacts)
    if log_path.resolve() not in {
        artifact.resolve() for artifact in artifact_paths
    }:
        raise ValueError("the PySCF log must be included in accepted artifacts")

    flush_pyscf_stdout(mol, mf)
    results_partial = output_dir / "results.partial.json"
    results_final = output_dir / "results.json"
    _write_json_atomic(results_partial, results)
    os.replace(results_partial, results_final)

    records: dict[str, dict[str, Any]] = {}
    for artifact in (*artifact_paths, results_final):
        record = _artifact_record(Path(artifact), output_dir)
        records[record["path"]] = record

    log_record = _artifact_record(log_path, output_dir)
    recorded_log_hash = records[log_record["path"]]["sha256"]

    accepted = dict(manifest)
    accepted.update(
        {
            "state": "accepted",
            "completed_at_utc": utc_now(),
            "elapsed_seconds": _require_finite_float(
                elapsed_seconds, "elapsed_seconds"
            ),
            "pyscf_version": str(pyscf_version),
            "artifacts": records,
        }
    )
    accepted_path = output_dir / "run_manifest.json"
    try:
        _write_json_atomic(accepted_path, accepted)
        flush_pyscf_stdout(mol, mf)
        post_publish_log_hash = sha256_file(log_path)
        if post_publish_log_hash != recorded_log_hash:
            raise RuntimeError(
                "PySCF log changed after the accepted manifest was published"
            )
    except Exception:
        accepted_path.unlink(missing_ok=True)
        results_final.unlink(missing_ok=True)
        raise
    (output_dir / "run_manifest.partial.json").unlink()
    return accepted


def finalize_failed_run(
    *,
    output_dir: Path,
    manifest: Mapping[str, Any],
    error: BaseException,
    pyscf_version: str,
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Keep a partial manifest that explicitly records a failed state."""

    failed = dict(manifest)
    failed.update(
        {
            "state": "failed",
            "completed_at_utc": utc_now(),
            "elapsed_seconds": _require_finite_float(
                elapsed_seconds, "elapsed_seconds"
            ),
            "pyscf_version": str(pyscf_version),
            "failure": {
                "type": type(error).__name__,
                "message": str(error),
            },
        }
    )
    _write_json_atomic(Path(output_dir) / "run_manifest.partial.json", failed)
    return failed


def _load_json_without_duplicate_keys(data: bytes) -> Mapping[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON field: {key}")
            result[key] = value
        return result

    try:
        parsed = json.loads(data.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("configuration must be valid UTF-8 JSON") from exc
    if not isinstance(parsed, Mapping):
        raise ValueError("the configuration root must be a JSON object")
    return parsed


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    """Execute the supported closed-shell RKS/TDA workflow."""

    config_path = Path(config_path)
    config_bytes = config_path.read_bytes()
    config = validate_config(_load_json_without_duplicate_keys(config_bytes))
    config_hash = sha256_bytes(config_bytes)
    started = time.monotonic()
    manifest = initialize_run(output_dir, config_bytes, config_hash)
    pyscf_version = "not_loaded"

    try:
        import pyscf
        from pyscf import dft, gto, lib

        pyscf_version = str(pyscf.__version__)
        manifest["pyscf_version"] = pyscf_version
        _write_json_atomic(
            Path(output_dir) / "run_manifest.partial.json", manifest
        )

        checkpoint = Path(output_dir) / "scf.chk"
        log_file = Path(output_dir) / "pyscf.log"
        mol = gto.M(
            atom=config["atom"],
            unit=config["unit"],
            basis=config["basis"],
            charge=config["charge"],
            spin=config["spin"],
            symmetry=False,
            verbose=4,
            output=str(log_file),
        )
        if mol.nelectron <= 0 or mol.nelectron < abs(mol.spin):
            raise ValueError("charge and spin imply an invalid electron count")
        if (mol.nelectron - mol.spin) % 2:
            raise ValueError("charge and spin imply non-integral spin counts")

        mf = dft.RKS(mol)
        mf.xc = config["xc"]
        mf.grids.level = config["grid_level"]
        mf.conv_tol = config["conv_tol"]
        mf.max_cycle = config["max_cycle"]
        mf.chkfile = str(checkpoint)
        energy_hartree = mf.kernel()
        validate_scf_result(mf, energy_hartree, checkpoint)
        scf_stability = validate_scf_stability(mf)

        td = mf.TDA()
        td.nstates = config["nstates"]
        td.kernel()
        oscillator_strength = td.oscillator_strength()
        validate_tda_result(
            td.converged,
            td.e,
            oscillator_strength,
            td.nstates,
        )

        excitation_hartree = [float(value) for value in td.e]
        strengths = [float(value) for value in oscillator_strength]
        results = {
            "method": "RKS/TDA",
            "xc": config["xc"],
            "basis": config["basis"],
            "charge": config["charge"],
            "spin": config["spin"],
            "coordinate_units": config["unit"],
            "atom_count": int(mol.natm),
            "electron_count": int(mol.nelectron),
            "scf_converged": True,
            "scf_energy_hartree": float(mf.e_tot),
            "scf_stability": scf_stability,
            "tda_roots_requested": td.nstates,
            "tda_roots_accepted": len(excitation_hartree),
            "excitation_energy_hartree": excitation_hartree,
            "excitation_energy_ev": [
                value * float(lib.param.HARTREE2EV)
                for value in excitation_hartree
            ],
            "oscillator_strength_length_gauge": strengths,
        }
        return finalize_accepted_run(
            output_dir=Path(output_dir),
            manifest=manifest,
            results=results,
            artifacts=(checkpoint, log_file),
            log_path=log_file,
            mol=mol,
            mf=mf,
            pyscf_version=pyscf_version,
            elapsed_seconds=time.monotonic() - started,
        )
    except Exception as error:
        try:
            finalize_failed_run(
                output_dir=Path(output_dir),
                manifest=manifest,
                error=error,
                pyscf_version=pyscf_version,
                elapsed_seconds=time.monotonic() - started,
            )
        except Exception as manifest_error:
            raise RuntimeError(
                f"run failed and the failed manifest could not be written: "
                f"{manifest_error}"
            ) from error
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the supported fail-closed, closed-shell PySCF RKS/TDA "
            "workflow in a fresh output directory."
        )
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        accepted = run(args.config, args.output)
    except Exception as error:
        print(f"PySCF run failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(accepted, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
