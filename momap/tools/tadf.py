#!/usr/bin/env python3
"""
momap-tadf — Full TADF photophysics pipeline using MOMAP.

Given a directory of Gaussian S0/T1/S1 outputs, runs:
  1. EVC (S1→S0) → Duschinsky + HR for fluorescence
  2. spec_tvcf (S1→S0) → fluorescence spectrum
  3. Optional EVC + ISC rate (S1→T1) when an explicit Hso is supplied
  4. Summary → ΔE_ST and spectral peak in the blue window
"""
import sys
import os
import re
import json
import math
import argparse
import hashlib
import shutil
import tempfile
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from extract import (
    extract_scf_energy,
    extract_last_excitation_ev,
    extract_s1_total_energy,
    extract_state1_transition_endpoints,
    validate_gaussian_log_contract,
    generate_spec_tvcf_input,
)
from oled import generate_isc_input, parse_isc_output
from runner import (
    stage_gaussian_inputs,
    run_local,
    validate_execution_evidence,
    validate_expected_identity,
)

AU2DEBYE = 2.541746
HA2EV = 27.2114
_MOL_ID_PATTERN = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$')
_EVC_NORMAL_MARKER = 'Normal finish of evc calculation'
_SPEC_NORMAL_MARKER = 'Normal finish of spec_tvcf calculation'
_ISC_NORMAL_MARKER = 'Normal finish of isc_tvcf calculation'
_HC_EV_NM = 1239.8419843320026
_EV_TO_WAVENUMBER_CM1 = 8065.544005
_SPECTRUM_COLUMN_COUNT = 7
_SPECTRUM_MIN_POINTS = 3
_FATAL_DIAGNOSTIC_PATTERN = re.compile(
    r'\bfatal\b|segmentation fault|error termination|\baborted\b|'
    r'traceback \(most recent call last\)|forrtl:',
    flags=re.IGNORECASE,
)

def run_momap_in_dir(
    inputfile,
    workdir,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Run hash-bound MOMAP in workdir and return execution evidence."""
    workdir = Path(workdir).expanduser().resolve()
    inputfile = Path(inputfile).expanduser().resolve(strict=True)
    return run_local(
        inputfile,
        workdir,
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
        expected_version_banner=expected_version_banner,
    )

def parse_evc_out(evc_out):
    """Extract reorganization energy and mode count from evc.out."""
    with open(evc_out) as f:
        content = f.read()
    
    m = re.search(r'(\d+)\s+# num of atoms', content)
    natoms = int(m.group(1)) if m else 0
    m = re.search(r'(\d+)\s+# num of modes', content)
    nmodes = int(m.group(1)) if m else 0
    
    return {'natoms': natoms, 'nmodes': nmodes}


def sha256_file(path):
    """Return the SHA-256 digest of one regular file."""
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _prepare_json_output(value, protected_paths):
    """Resolve a fresh JSON target without creating or replacing anything."""
    lexical_target = Path(value).expanduser()
    if lexical_target.exists() or lexical_target.is_symlink():
        raise FileExistsError(
            f'Refusing to overwrite or follow JSON output: {lexical_target}'
        )
    parent = lexical_target.parent.resolve(strict=True)
    if not parent.is_dir():
        raise FileNotFoundError(
            f'JSON output parent directory does not exist: {parent}'
        )
    target = parent / lexical_target.name
    protected = {
        Path(path).expanduser().resolve()
        for path in protected_paths
        if path is not None
    }
    if target in protected:
        raise ValueError('JSON output must be distinct from every Gaussian input')
    return target


def _write_json_exclusive(path, payload):
    """Atomically publish JSON without replacing an existing path."""
    path = Path(path)
    fd, partial_name = tempfile.mkstemp(
        prefix=f'.{path.name}.', suffix='.partial', dir=path.parent
    )
    partial = Path(partial_name)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, indent=2, default=str)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.link(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def _fresh_nonempty_file(path, started_ns):
    """Validate that ``path`` was produced after the current stage began."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f'Required output is missing: {path}')
    stat_result = path.stat()
    if stat_result.st_size <= 0:
        raise ValueError(f'Required output must be non-empty: {path}')
    if stat_result.st_mtime_ns < started_ns:
        raise ValueError(
            f'Required output is not fresh (older than stage start): {path}'
        )
    return path


def _require_marker(path, marker):
    text = Path(path).read_text(encoding='utf-8', errors='replace')
    _reject_fatal_diagnostics(text, path)
    matches = [
        line
        for line in text.splitlines()
        if line.strip().casefold() == marker.casefold()
    ]
    if len(matches) != 1:
        raise ValueError(
            f'{path} must contain exactly one full-line normal completion '
            f'marker {marker!r}; found {len(matches)}'
        )
    return text


def _reject_fatal_diagnostics(text, path):
    match = _FATAL_DIAGNOSTIC_PATTERN.search(text)
    if match:
        raise ValueError(
            f'Fatal diagnostic in {path}: {match.group(0)!r}'
        )


def validate_evc_stage(
    workdir,
    started_ns,
    execution,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Accept EVC only from fresh artifacts with positive atom/mode counts."""
    workdir = Path(workdir)
    execution = validate_execution_evidence(
        execution,
        workdir,
        started_ns,
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
        expected_version_banner=expected_version_banner,
    )
    if execution['process_exit_code'] != 0:
        raise ValueError('MOMAP EVC process returned nonzero')
    evc_out = _fresh_nonempty_file(workdir / 'evc.out', started_ns)
    evc_cart = _fresh_nonempty_file(workdir / 'evc.cart.dat', started_ns)
    _require_marker(evc_out, _EVC_NORMAL_MARKER)
    info = parse_evc_out(evc_out)
    if info['natoms'] <= 0 or info['nmodes'] <= 0:
        raise ValueError(
            'EVC output must report a positive atom count and mode count'
        )
    info['evc_out'] = str(evc_out)
    info['evc_cart_dat'] = str(evc_cart)
    info['execution'] = execution
    return info


def parse_isc_rate_output(rate_file, *, expected_build):
    """Use the standalone OLED parser as the single ISC rate contract."""
    return parse_isc_output(rate_file, expected_build=expected_build)


def validate_isc_stage(
    workdir,
    started_ns,
    execution,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Accept ISC only after validating fresh log and finite rate artifacts."""
    workdir = Path(workdir)
    execution = validate_execution_evidence(
        execution,
        workdir,
        started_ns,
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
        expected_version_banner=expected_version_banner,
    )
    if execution['process_exit_code'] != 0:
        raise ValueError('MOMAP ISC process returned nonzero')
    log_file = _fresh_nonempty_file(workdir / 'isc.tvcf.log', started_ns)
    rate_file = _fresh_nonempty_file(workdir / 'isc.tvcf.fo.dat', started_ns)
    _require_marker(log_file, _ISC_NORMAL_MARKER)
    rates = parse_isc_rate_output(rate_file, expected_build=expected_build)
    rates['execution'] = execution
    return rates


def validate_spectrum_stage(
    workdir,
    started_ns,
    execution,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Accept a spectrum only from fresh artifacts and a normal marker."""
    workdir = Path(workdir)
    execution = validate_execution_evidence(
        execution,
        workdir,
        started_ns,
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
        expected_version_banner=expected_version_banner,
    )
    if execution['process_exit_code'] != 0:
        raise ValueError('MOMAP spectrum process returned nonzero')
    log_file = _fresh_nonempty_file(workdir / 'spec.tvcf.log', started_ns)
    spec_file = _fresh_nonempty_file(
        workdir / 'spec.tvcf.spec.dat', started_ns
    )
    _require_marker(log_file, _SPEC_NORMAL_MARKER)
    spectrum = parse_spec_output(spec_file, expected_build=expected_build)
    spectrum['execution'] = execution
    return spectrum


def _new_stage_directory(mol_dir, name):
    stage_dir = Path(mol_dir) / name
    stage_dir.mkdir()
    return stage_dir


def _stage_files(stage_dir, paths):
    staged = []
    for source in paths:
        source = Path(source).resolve(strict=True)
        target = Path(stage_dir) / source.name
        if target.exists():
            raise FileExistsError(f'Stage input would overwrite {target}')
        shutil.copy2(source, target)
        staged.append(target)
    return staged


def write_stage_receipt(
    stage_dir,
    stage,
    started_ns,
    inputs,
    outputs,
    metadata=None,
    *,
    executable_identity,
):
    """Write a hash-bound receipt after a stage has passed validation."""
    stage_dir = Path(stage_dir)
    receipt_path = stage_dir / 'stage_receipt.json'
    if receipt_path.exists():
        raise FileExistsError(f'Receipt already exists: {receipt_path}')

    def describe(path):
        path = Path(path).resolve(strict=True)
        stat_result = path.stat()
        return {
            'path': str(path.relative_to(stage_dir.resolve())),
            'sha256': sha256_file(path),
            'size_bytes': stat_result.st_size,
            'mtime_ns': stat_result.st_mtime_ns,
        }

    payload = {
        'stage': stage,
        'status': 'accepted',
        'started_ns': started_ns,
        'recorded_ns': time.time_ns(),
        'inputs': [describe(path) for path in inputs],
        'outputs': [describe(path) for path in outputs],
        'metadata': metadata or {},
        'executable_identity': executable_identity,
    }
    receipt_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8'
    )
    return receipt_path

def parse_spec_output(spec_file, *, expected_build):
    """Parse the strict seven-column local MOMAP 2024A spectrum fixture."""
    if expected_build != '2024A':
        raise ValueError(
            'The strict seven-column spectrum parser requires explicit '
            'expected_build="2024A"'
        )
    spec_file = Path(spec_file)
    if not spec_file.is_file():
        raise ValueError(f'Spectrum output is missing: {spec_file}')
    text = spec_file.read_text(encoding='utf-8', errors='replace')
    _reject_fatal_diagnostics(text, spec_file)

    data = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = stripped.split()
        if len(parts) != _SPECTRUM_COLUMN_COUNT:
            raise ValueError(
                f'Spectrum line {line_number} must contain exactly '
                f'{_SPECTRUM_COLUMN_COUNT} columns, found {len(parts)}'
            )
        try:
            row = [float(value) for value in parts]
        except ValueError as exc:
            raise ValueError(
                f'Spectrum line {line_number} contains a non-numeric value'
            ) from exc
        if not all(math.isfinite(value) for value in row):
            raise ValueError(
                f'Spectrum line {line_number} contains a non-finite value'
            )
        data.append(row)

    if len(data) < _SPECTRUM_MIN_POINTS:
        raise ValueError(
            f'Spectrum requires at least {_SPECTRUM_MIN_POINTS} data points; '
            f'found {len(data)}'
        )

    energy_ev = [row[1] for row in data]
    wavelengths = [row[3] for row in data]
    energy_hartree = [row[0] for row in data]
    wavenumbers_cm1 = [row[2] for row in data]
    intensities = [row[index] for row in data for index in (4, 5, 6)]
    if any(value <= 0 for value in energy_hartree):
        raise ValueError('Spectrum Hartree energies must be positive')
    if any(value <= 0 for value in energy_ev):
        raise ValueError('Spectrum eV energies must be positive')
    if any(value <= 0 for value in wavenumbers_cm1):
        raise ValueError('Spectrum wavenumbers must be positive')
    if any(value <= 0 for value in wavelengths):
        raise ValueError('Spectrum wavelengths must be positive')
    if any(value < 0 for value in intensities):
        raise ValueError('Spectrum intensity columns must be non-negative')

    energy_differences = [
        right - left for left, right in zip(energy_ev, energy_ev[1:])
    ]
    increasing = all(difference > 0 for difference in energy_differences)
    decreasing = all(difference < 0 for difference in energy_differences)
    if not (increasing or decreasing):
        raise ValueError('Spectrum eV energy axis must be strictly monotonic')

    wavelength_differences = [
        right - left for left, right in zip(wavelengths, wavelengths[1:])
    ]
    if increasing:
        wavelength_consistent_order = all(
            difference < 0 for difference in wavelength_differences
        )
    else:
        wavelength_consistent_order = all(
            difference > 0 for difference in wavelength_differences
        )
    if not wavelength_consistent_order:
        raise ValueError(
            'Spectrum wavelength axis must be strictly monotonic opposite to eV'
        )

    for row_number, (energy, wavelength) in enumerate(
        zip(energy_ev, wavelengths), start=1
    ):
        relative_error = abs(energy * wavelength - _HC_EV_NM) / _HC_EV_NM
        if relative_error > 0.01:
            raise ValueError(
                f'Spectrum eV/nm values are inconsistent at data row '
                f'{row_number}: {energy:g} eV and {wavelength:g} nm'
            )
    for row_number, row in enumerate(data, start=1):
        hartree_relative_error = abs(row[0] * HA2EV - row[1]) / row[1]
        if hartree_relative_error > 0.01:
            raise ValueError(
                f'Spectrum Hartree/eV values are inconsistent at data row '
                f'{row_number}: {row[0]:g} Hartree and {row[1]:g} eV'
            )
        wavenumber_relative_error = (
            abs(row[1] * _EV_TO_WAVENUMBER_CM1 - row[2]) / row[2]
        )
        if wavenumber_relative_error > 0.01:
            raise ValueError(
                f'Spectrum eV/wavenumber values are inconsistent at data row '
                f'{row_number}: {row[1]:g} eV and {row[2]:g} cm^-1'
            )

    emi = [row[5] for row in data]
    peaks = []
    for index in range(1, len(emi) - 1):
        if emi[index] > emi[index - 1] and emi[index] > emi[index + 1] and emi[index] > 0:
            peaks.append({
                'wavelength': wavelengths[index],
                'intensity': emi[index],
            })
    peaks.sort(key=lambda item: item['intensity'], reverse=True)
    max_index = max(range(len(emi)), key=lambda index: emi[index])

    return {
        'peak_wavelength': wavelengths[max_index],
        'peak_intensity': emi[max_index],
        'top_peaks': peaks[:5],
        'data_points': len(data),
        'spectrum_contract': 'MOMAP_2024A_exact_7_columns',
    }

def process_molecule(
    mol_id,
    s0_log,
    s1_log,
    t1_log,
    output_dir,
    temperature=300,
    hso_cm1=None,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Run the TADF MOMAP pipeline for one molecule.

    ISC is intentionally opt-in: ``hso_cm1`` must be supplied from a measured
    or computed spin-orbit coupling. A missing value is reported as
    ``not_computed`` and is never replaced by a numerical placeholder.
    """
    results = {'mol_id': mol_id, 'success': False}
    if not isinstance(mol_id, str) or not _MOL_ID_PATTERN.fullmatch(mol_id):
        results['error'] = (
            'mol_id must contain only letters, digits, dots, underscores, or '
            'hyphens and must not contain a path'
        )
        return results

    try:
        expected_launcher_sha256 = validate_expected_identity(
            expected_build,
            expected_launcher_sha256,
            expected_version_banner,
        )
    except ValueError as exc:
        results['error'] = f'MOMAP executable expectation failed: {exc}'
        return results
    identity_expectations = {
        'expected_build': expected_build,
        'expected_launcher_sha256': expected_launcher_sha256,
        'expected_version_banner': expected_version_banner,
    }

    try:
        temperature = float(temperature)
    except (TypeError, ValueError):
        results['error'] = 'temperature must be a finite positive value in K'
        return results
    if not math.isfinite(temperature) or temperature <= 0:
        results['error'] = 'temperature must be a finite positive value in K'
        return results

    if hso_cm1 is not None:
        try:
            hso_cm1 = float(hso_cm1)
        except (TypeError, ValueError):
            results['error'] = 'hso_cm1 must be a finite positive number'
            return results
        if not math.isfinite(hso_cm1) or hso_cm1 <= 0:
            results['error'] = 'hso_cm1 must be a finite positive number'
            return results
        if t1_log is None:
            results['error'] = 'A T1 Gaussian log is required when hso_cm1 is set'
            return results

    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    mol_dir = output_root / mol_id
    try:
        mol_dir.mkdir()
    except FileExistsError:
        results['error'] = f'Output directory already exists: {mol_dir}'
        return results

    print(f"\n{'='*60}")
    print(f"🧬 Processing {mol_id}")
    print(f"{'='*60}")

    source_logs = {'S0': s0_log, 'S1': s1_log, 'T1': t1_log}
    staged_logs = {}
    gaussian_contracts = {}
    try:
        for label, path in source_logs.items():
            if path is None:
                print(f"  ⚠️  {label}: no file provided")
                continue
            source = Path(path).expanduser().resolve(strict=True)
            staged_log, staged_fchk = stage_gaussian_inputs(
                source, mol_dir, target_stem=label.lower()
            )
            contract = validate_gaussian_log_contract(
                staged_log, staged_fchk, state_label=label
            )
            print(
                f"  📄 {label}: {source} "
                f"({contract['normal_job_count']} validated opt/freq jobs)"
            )
            staged_logs[label] = staged_log
            gaussian_contracts[label] = contract

        reference_contract = gaussian_contracts.get('S0')
        if reference_contract is not None:
            for label in ('S1', 'T1'):
                contract = gaussian_contracts.get(label)
                if contract is None:
                    continue
                if contract['charge'] != reference_contract['charge']:
                    raise ValueError(
                        f'{label} and S0 charges differ; all states must belong '
                        'to the same molecular system'
                    )
                if contract['atomic_numbers'] != reference_contract['atomic_numbers']:
                    raise ValueError(
                        f'{label} and S0 atom identities/order differ; all states '
                        'must belong to the same molecular system'
                    )
    except (OSError, ValueError) as exc:
        results['error'] = f"Gaussian input staging failed: {exc}"
        print(f"  ❌ {results['error']}")
        return results

    if 'S0' not in staged_logs or 'S1' not in staged_logs:
        results['error'] = 'Both S0 and S1 Gaussian logs are required'
        print(f"  ❌ {results['error']}")
        return results

    s0_staged = staged_logs['S0']
    s1_staged = staged_logs['S1']
    t1_staged = staged_logs.get('T1')

    try:
        E_S0 = extract_scf_energy(s0_staged)
        E_S1_scf = extract_scf_energy(s1_staged)
        E_exc_ev = extract_last_excitation_ev(s1_staged)
        E_S1 = extract_s1_total_energy(s1_staged)
        E_T1 = extract_scf_energy(t1_staged) if t1_staged else None
        absorption, emission = extract_state1_transition_endpoints(s1_staged)
    except (OSError, ValueError) as exc:
        results['error'] = f"Gaussian parameter extraction failed: {exc}"
        print(f"  ❌ {results['error']}")
        return results

    delta_EST = (E_S1 - E_T1) * HA2EV if E_T1 is not None else None
    state_order_status = None
    if delta_EST is not None:
        state_order_status = (
            'expected_S1_above_T1'
            if delta_EST > 0
            else 'unexpected_S1_not_above_T1'
        )
    results.update({
        'E_S0': E_S0,
        'E_S1': E_S1,
        'E_S1_scf': E_S1_scf,
        'E_exc_ev': E_exc_ev,
        'E_T1': E_T1,
        'E_T1_energy_type': (
            'Gaussian SCF total energy' if E_T1 is not None else None
        ),
        'Ead_S1_S0': (E_S1 - E_S0) * HA2EV,
        'delta_EST_eV': delta_EST,
        'delta_EST_signed_eV': delta_EST,
        'state_order_status': state_order_status,
        'EDMA': absorption['magnitude_debye'],
        'EDME': emission['magnitude_debye'],
        'f_abs': absorption['Osc'],
        'f_emi': emission['Osc'],
        'gaussian_contracts': gaussian_contracts,
    })

    print(f"\n  📊 Extracted parameters:")
    print(f"     E(S0)   = {E_S0:.8f} au")
    print(f"     E(S1)   = {E_S1:.8f} au")
    if E_T1 is not None:
        print(f"     E(T1)   = {E_T1:.8f} au")
        print(f"     Signed ΔE_ST = {delta_EST:.4f} eV")
    print(f"     EDMA    = {results['EDMA']:.4f} debye")
    print(f"     EDME    = {results['EDME']:.4f} debye")

    if hso_cm1 is not None and delta_EST is not None and delta_EST <= 0:
        results['isc'] = {
            'status': 'failed',
            'reason': (
                'Signed S1-T1 gap is non-positive; verify state identity and '
                'energy conventions before ISC'
            ),
            'hso_cm1': hso_cm1,
        }
        results['error'] = results['isc']['reason']
        return results

    print(f"\n  📐 Step 1: EVC (S1→S0)")
    try:
        evc_s1_dir = _new_stage_directory(mol_dir, 'evc_s1_s0')
        evc_s1_files = _stage_files(
            evc_s1_dir,
            (
                s0_staged,
                s0_staged.with_suffix('.fchk'),
                s1_staged,
                s1_staged.with_suffix('.fchk'),
            ),
        )
        evc_s1_input = evc_s1_dir / 'momap_evc_s1.inp'
        evc_s1_input.write_text(
            f"""do_evc = 1

&evc
  ffreq(1) = "{s0_staged.name}"
  ffreq(2) = "{s1_staged.name}"
  sort_mode = 1
/
""",
            encoding='utf-8',
        )
        evc_s1_started = time.time_ns()
        evc_s1_execution = run_momap_in_dir(
            evc_s1_input, evc_s1_dir, **identity_expectations
        )
        evc_s1_info = validate_evc_stage(
            evc_s1_dir,
            evc_s1_started,
            evc_s1_execution,
            **identity_expectations,
        )
        evc_s1_execution = evc_s1_info.pop('execution')
        evc_s1_receipt = write_stage_receipt(
            evc_s1_dir,
            'evc_s1_s0',
            evc_s1_started,
            [evc_s1_input, *evc_s1_files],
            [
                evc_s1_dir / 'evc.out',
                evc_s1_dir / 'evc.cart.dat',
                Path(evc_s1_execution['run_log']),
            ],
            {
                'process_exit_code': evc_s1_execution['process_exit_code'],
                'natoms': evc_s1_info['natoms'],
                'nmodes': evc_s1_info['nmodes'],
            },
            executable_identity=evc_s1_execution['executable_identity'],
        )
    except (OSError, RuntimeError, ValueError) as exc:
        results['error'] = f'EVC (S1 to S0) failed acceptance: {exc}'
        print(f"  ❌ {results['error']}")
        return results
    print(f"  ✅ {evc_s1_info['natoms']} atoms, {evc_s1_info['nmodes']} modes")

    print(f"\n  🌈 Step 2: Spectrum (S1→S0)")
    try:
        spectrum_dir = _new_stage_directory(mol_dir, 'spectrum')
        spectrum_evc = _stage_files(
            spectrum_dir, [evc_s1_dir / 'evc.cart.dat']
        )[0]
        spec_params = {
            'temperature': temperature,
            'Ead': E_S1 - E_S0,
            'EDMA': results['EDMA'],
            'EDME': results['EDME'],
            'dsfile': spectrum_evc.name,
        }
        spec_input = spectrum_dir / 'momap_spec.inp'
        generate_spec_tvcf_input(spec_params, spec_input)
        spectrum_started = time.time_ns()
        spectrum_execution = run_momap_in_dir(
            spec_input, spectrum_dir, **identity_expectations
        )
        spec_info = validate_spectrum_stage(
            spectrum_dir,
            spectrum_started,
            spectrum_execution,
            **identity_expectations,
        )
        spectrum_execution = spec_info.pop('execution')
        peak_nm = spec_info.get('peak_wavelength')
        peak_intensity = spec_info.get('peak_intensity')
        valid_peak = (
            isinstance(peak_nm, (int, float))
            and not isinstance(peak_nm, bool)
            and math.isfinite(peak_nm)
            and peak_nm > 0
            and isinstance(peak_intensity, (int, float))
            and not isinstance(peak_intensity, bool)
            and math.isfinite(peak_intensity)
            and peak_intensity > 0
        )
        if not valid_peak:
            raise ValueError(
                'spectrum output contains no valid positive-intensity peak'
            )
        spectrum_receipt = write_stage_receipt(
            spectrum_dir,
            'spectrum_s1_s0',
            spectrum_started,
            [spec_input, spectrum_evc],
            [
                spectrum_dir / 'spec.tvcf.log',
                spectrum_dir / 'spec.tvcf.spec.dat',
                Path(spectrum_execution['run_log']),
            ],
            {
                'process_exit_code': spectrum_execution['process_exit_code'],
                'peak_wavelength_nm': peak_nm,
                'peak_intensity': peak_intensity,
                'data_points': spec_info.get('data_points'),
            },
            executable_identity=spectrum_execution['executable_identity'],
        )
    except (OSError, RuntimeError, ValueError) as exc:
        results['spectrum'] = {'status': 'failed'}
        results['error'] = f'Spectrum calculation failed acceptance: {exc}'
        print(f"  ❌ {results['error']}")
        return results

    results['spectrum'] = spec_info
    print(f"  ✅ Peak: {peak_nm:.1f} nm")
    if 450 <= peak_nm <= 490:
        print("  🔵 WITHIN BLUE WINDOW (450–490 nm)!")
        results['blue_window'] = True
    elif 400 <= peak_nm <= 500:
        print(f"  🔹 Near blue window ({peak_nm:.0f} nm)")
        results['blue_window'] = 'near'
    else:
        print("  ⚪ Outside blue window")
        results['blue_window'] = False
    if spec_info.get('top_peaks'):
        results['top_spectral_peaks'] = spec_info['top_peaks'][:3]

    if t1_staged is None:
        results['isc'] = {
            'status': 'not_computed',
            'reason': 'T1 Gaussian log not provided',
        }
    elif hso_cm1 is None:
        results['isc'] = {
            'status': 'not_computed',
            'reason': 'hso_cm1 not provided',
        }
        print("\n  ⏭️  ISC not computed: hso_cm1 was not provided")
    else:
        print(f"\n  🔀 Step 3: ISC TVCF (S1→T1 phosphorescence)")
        try:
            evc_isc_dir = _new_stage_directory(mol_dir, 'evc_s1_t1')
            evc_isc_files = _stage_files(
                evc_isc_dir,
                (
                    s1_staged,
                    s1_staged.with_suffix('.fchk'),
                    t1_staged,
                    t1_staged.with_suffix('.fchk'),
                ),
            )
            evc_isc_input = evc_isc_dir / 'momap_evc_isc.inp'
            evc_isc_input.write_text(
                f"""do_evc = 1

&evc
  ffreq(1) = "{s1_staged.name}"
  ffreq(2) = "{t1_staged.name}"
  sort_mode = 1
/
""",
                encoding='utf-8',
            )
            evc_isc_started = time.time_ns()
            evc_isc_execution = run_momap_in_dir(
                evc_isc_input, evc_isc_dir, **identity_expectations
            )
            evc_isc_info = validate_evc_stage(
                evc_isc_dir,
                evc_isc_started,
                evc_isc_execution,
                **identity_expectations,
            )
            evc_isc_execution = evc_isc_info.pop('execution')
            evc_isc_receipt = write_stage_receipt(
                evc_isc_dir,
                'evc_s1_t1',
                evc_isc_started,
                [evc_isc_input, *evc_isc_files],
                [
                    evc_isc_dir / 'evc.out',
                    evc_isc_dir / 'evc.cart.dat',
                    Path(evc_isc_execution['run_log']),
                ],
                {
                    'process_exit_code': evc_isc_execution['process_exit_code'],
                    'natoms': evc_isc_info['natoms'],
                    'nmodes': evc_isc_info['nmodes'],
                    'delta_EST_signed_eV': delta_EST,
                },
                executable_identity=evc_isc_execution['executable_identity'],
            )

            isc_dir = _new_stage_directory(mol_dir, 'isc')
            isc_evc = _stage_files(
                isc_dir, [evc_isc_dir / 'evc.cart.dat']
            )[0]
            isc_input = isc_dir / 'momap_isc.inp'
            generate_isc_input(
                isc_evc,
                E_S1 - E_T1,
                hso_cm1,
                isc_input,
                temp=temperature,
                tmax=5000,
            )
            isc_started = time.time_ns()
            isc_execution = run_momap_in_dir(
                isc_input, isc_dir, **identity_expectations
            )
            isc_rates = validate_isc_stage(
                isc_dir,
                isc_started,
                isc_execution,
                **identity_expectations,
            )
            isc_execution = isc_rates.pop('execution')
            isc_receipt = write_stage_receipt(
                isc_dir,
                'isc_s1_t1',
                isc_started,
                [isc_input, isc_evc],
                [
                    isc_dir / 'isc.tvcf.log',
                    isc_dir / 'isc.tvcf.fo.dat',
                    Path(isc_execution['run_log']),
                ],
                {
                    'process_exit_code': isc_execution['process_exit_code'],
                    'hso_cm-1': hso_cm1,
                    'Ead_signed_au': E_S1 - E_T1,
                    **isc_rates,
                },
                executable_identity=isc_execution['executable_identity'],
            )
        except (OSError, RuntimeError, ValueError) as exc:
            results['isc'] = {
                'status': 'failed',
                'reason': str(exc),
                'hso_cm1': hso_cm1,
            }
            results['error'] = f'ISC calculation failed acceptance: {exc}'
            print(f"  ❌ {results['error']}")
            return results

        results['isc'] = {
            'status': 'computed',
            'hso_cm1': hso_cm1,
            'Ead_signed_au': E_S1 - E_T1,
            **isc_rates,
        }
        print("  ✅ ISC TVCF accepted with explicit s^-1 rates")

    results['stage_receipts'] = {
        'evc_s1_s0': str(evc_s1_receipt),
        'spectrum': str(spectrum_receipt),
    }
    if hso_cm1 is not None:
        results['stage_receipts'].update({
            'evc_s1_t1': str(evc_isc_receipt),
            'isc': str(isc_receipt),
        })

    results['success'] = True
    results['output_dir'] = str(mol_dir)

    print(f"\n  {'='*50}")
    print(f"  📋 {mol_id} Summary")
    print(f"  {'='*50}")
    print(
        f"  ΔE_ST     = {delta_EST:.4f} eV"
        if delta_EST is not None
        else "  ΔE_ST     = N/A"
    )
    if results.get('spectrum', {}).get('peak_wavelength'):
        spectrum = results['spectrum']
        print(f"  λ_emi      = {spectrum['peak_wavelength']:.1f} nm")
        blue = results.get('blue_window')
        if blue is True:
            print("  Blue Window = ✅ YES")
        elif blue == 'near':
            print("  Blue Window = 🟡 NEAR")
        else:
            print("  Blue Window = ❌ NO")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='MOMAP TADF pipeline — auto EVC → spectrum → ISC for TADF candidates'
    )
    parser.add_argument('mol_id', help='Molecule identifier')
    parser.add_argument('--s0', required=True, help='S0 ground state Gaussian log')
    parser.add_argument('--s1', required=True, help='S1 excited state Gaussian log')
    parser.add_argument('--t1', help='T1 triplet state Gaussian log')
    parser.add_argument('--output', '-o', default='./momap_tadf_output', help='Output directory')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument(
        '--hso-cm1',
        type=float,
        help='Explicit S1-T1 spin-orbit coupling in cm^-1; omit to skip ISC',
    )
    parser.add_argument(
        '--json-output',
        help='Write machine-readable results to this JSON file',
    )
    parser.add_argument(
        '--expected-build',
        required=True,
        choices=['2024A'],
        help='Verified MOMAP build contract (currently 2024A only)',
    )
    parser.add_argument(
        '--expected-launcher-sha256',
        required=True,
        help='Expected SHA-256 of the original licensed momap launcher',
    )
    parser.add_argument(
        '--expected-version-banner',
        required=True,
        help='Exact full version-banner line expected once per stage run log',
    )
    args = parser.parse_args()

    json_output = None
    if args.json_output:
        json_output = _prepare_json_output(
            args.json_output, (args.s0, args.s1, args.t1)
        )

    try:
        results = process_molecule(
            mol_id=args.mol_id,
            s0_log=args.s0,
            s1_log=args.s1,
            t1_log=args.t1,
            output_dir=args.output,
            temperature=args.temperature,
            hso_cm1=args.hso_cm1,
            expected_build=args.expected_build,
            expected_launcher_sha256=args.expected_launcher_sha256,
            expected_version_banner=args.expected_version_banner,
        )
    except Exception as exc:
        results = {
            'mol_id': args.mol_id,
            'success': False,
            'error': f'Unexpected pipeline failure: {exc}',
        }

    if json_output is not None:
        _write_json_exclusive(json_output, results)
    
    if results.get('success'):
        return 0
    else:
        print(f"\n❌ Processing failed for {args.mol_id}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
