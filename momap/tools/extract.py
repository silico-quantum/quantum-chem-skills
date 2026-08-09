#!/usr/bin/env python3
"""
momap-extract — Parse Gaussian output files, extract MOMAP parameters.
Generates a complete spec_tvcf input file ready to run.
"""
import sys
import re
import os
import argparse
import math
from pathlib import Path

AU2DEBYE = 2.541746
HA2EV = 27.2114
BOHR_TO_ANGSTROM = 0.529177210903
FCHK_GEOMETRY_ATOL_ANGSTROM = 1.0e-6

_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_GAUSSIAN_NORMAL_LINE = re.compile(
    r"^\s*Normal termination of Gaussian(?:\s+\S+)?(?:\s+at\s+.+)?\s*$",
    re.IGNORECASE,
)
_GAUSSIAN_FATAL = re.compile(
    r"Error termination|l9999\.exe|segmentation fault|\bfatal\b|"
    r"traceback \(most recent call last\)",
    re.IGNORECASE,
)
_CHARGE_MULTIPLICITY = re.compile(
    r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)",
    re.IGNORECASE,
)
_ORIENTATION_ROW = re.compile(
    rf"^\s*\d+\s+(\d+)\s+\d+\s+({_FLOAT_PATTERN})\s+"
    rf"({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s*$"
)
_ORIENTATION_UNITS = re.compile(
    r"Coordinates\s*\(\s*(Angstroms|Bohr)\s*\)",
    re.IGNORECASE,
)


def _as_float(value):
    """Parse Gaussian-style floating-point text, including Fortran D exponents."""
    return float(value.replace('D', 'E').replace('d', 'e'))


def _last_standard_orientation(lines):
    """Return atom order and coordinates from the last Standard orientation table."""
    blocks = []
    for start, line in enumerate(lines):
        if "Standard orientation:" not in line:
            continue
        atomic_numbers = []
        coordinates = []
        unit = None
        for candidate in lines[start + 1 :]:
            if "Standard orientation:" in candidate:
                break
            unit_match = _ORIENTATION_UNITS.search(candidate)
            if unit_match:
                unit = unit_match.group(1).lower()
                continue
            match = _ORIENTATION_ROW.fullmatch(candidate)
            if match:
                if unit is None:
                    raise ValueError(
                        "Gaussian Standard orientation lacks an explicit coordinate unit"
                    )
                xyz = tuple(_as_float(match.group(index)) for index in range(2, 5))
                if not all(math.isfinite(value) for value in xyz):
                    raise ValueError(
                        "Gaussian Standard orientation contains non-finite coordinates"
                    )
                scale = BOHR_TO_ANGSTROM if unit == "bohr" else 1.0
                atomic_numbers.append(int(match.group(1)))
                coordinates.append(tuple(value * scale for value in xyz))
                continue
            if atomic_numbers:
                break
        if atomic_numbers:
            blocks.append((atomic_numbers, coordinates))
    if not blocks:
        raise ValueError("Gaussian job lacks a complete Standard orientation table")
    return blocks[-1]


def _fchk_scalar(text, label):
    match = re.search(
        rf"^\s*{re.escape(label)}\s+I\s+(-?\d+)\s*$",
        text,
        re.MULTILINE,
    )
    if not match:
        raise ValueError(f"formatted checkpoint lacks {label}")
    return int(match.group(1))


def _fchk_atomic_numbers(text):
    match = re.search(
        r"^\s*Atomic numbers\s+I\s+N=\s*(\d+)\s*$",
        text,
        re.MULTILINE,
    )
    if not match:
        raise ValueError("formatted checkpoint lacks Atomic numbers")
    expected = int(match.group(1))
    values = []
    tail = text[match.end() :]
    for line in tail.splitlines():
        if not line.strip():
            if values:
                break
            continue
        if not re.fullmatch(r"\s*(?:\d+\s*)+", line):
            if values:
                break
            raise ValueError("formatted checkpoint Atomic numbers payload is malformed")
        values.extend(int(value) for value in line.split())
        if len(values) >= expected:
            break
    if len(values) != expected or any(value <= 0 for value in values):
        raise ValueError(
            "formatted checkpoint Atomic numbers count/content is inconsistent"
        )
    return values


def _fchk_cartesian_coordinates(text, atom_count):
    """Return fchk Current cartesian coordinates converted from Bohr to Angstrom."""
    match = re.search(
        r"^\s*Current cartesian coordinates\s+R\s+N=\s*(\d+)\s*$",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if not match:
        raise ValueError("formatted checkpoint lacks Current cartesian coordinates")
    declared = int(match.group(1))
    expected = 3 * atom_count
    if declared != expected:
        raise ValueError(
            "formatted checkpoint Current cartesian coordinates count is inconsistent"
        )

    values = []
    for line in text[match.end() :].splitlines():
        tokens = line.split()
        if not tokens:
            if values:
                break
            continue
        if not all(re.fullmatch(_FLOAT_PATTERN, token) for token in tokens):
            raise ValueError(
                "formatted checkpoint Current cartesian coordinates payload is malformed"
            )
        values.extend(_as_float(token) for token in tokens)
        if len(values) >= declared:
            break
    if len(values) != declared or not all(math.isfinite(value) for value in values):
        raise ValueError(
            "formatted checkpoint Current cartesian coordinates are incomplete or non-finite"
        )
    converted = [value * BOHR_TO_ANGSTROM for value in values]
    return [tuple(converted[index : index + 3]) for index in range(0, expected, 3)]


def validate_gaussian_log_contract(logpath, fchk_path, state_label):
    """Validate the supported two-job Gaussian opt/freq state chain.

    The MOMAP adapter deliberately supports one conservative layout: a normally
    terminated optimization job followed by a normally terminated frequency
    job for S0, S1, or T1. It binds charge, multiplicity, atom order, and the
    matching formatted checkpoint before any MOMAP stage can be accepted.
    """
    if state_label not in {"S0", "S1", "T1"}:
        raise ValueError(f"unsupported Gaussian state label: {state_label!r}")
    logpath = Path(logpath)
    fchk_path = Path(fchk_path)
    text = logpath.read_text(encoding="utf-8", errors="replace")
    if _GAUSSIAN_FATAL.search(text):
        raise ValueError(f"{state_label} Gaussian log contains a fatal termination")
    lines = text.splitlines()
    normal_indices = [
        index for index, line in enumerate(lines) if _GAUSSIAN_NORMAL_LINE.fullmatch(line)
    ]
    if len(normal_indices) != 2:
        raise ValueError(
            f"{state_label} Gaussian chain must contain exactly two exact normal "
            f"termination lines (optimization and frequency); found {len(normal_indices)}"
        )
    last_nonempty = max(
        (index for index, line in enumerate(lines) if line.strip()), default=-1
    )
    if normal_indices[-1] != last_nonempty:
        raise ValueError(
            f"{state_label} Gaussian final job does not end at the final normal termination"
        )
    jobs = [
        lines[: normal_indices[0] + 1],
        lines[normal_indices[0] + 1 : normal_indices[1] + 1],
    ]
    if not any("Optimization completed." in line for line in jobs[0]):
        raise ValueError(f"{state_label} Gaussian optimization is not completed")
    if not any("Harmonic frequencies (cm**-1)" in line for line in jobs[1]):
        raise ValueError(f"{state_label} Gaussian frequency job lacks its analysis header")
    if not any(re.search(r"^\s*Frequencies\s+--\s+", line) for line in jobs[1]):
        raise ValueError(f"{state_label} Gaussian frequency job has no frequencies")

    charge_multiplicity = []
    atom_orders = []
    geometries = []
    for job_number, job_lines in enumerate(jobs, start=1):
        pairs = {
            (int(charge), int(multiplicity))
            for charge, multiplicity in _CHARGE_MULTIPLICITY.findall("\n".join(job_lines))
        }
        if len(pairs) != 1:
            raise ValueError(
                f"{state_label} Gaussian job {job_number} lacks one unambiguous charge/multiplicity"
            )
        charge_multiplicity.append(next(iter(pairs)))
        atom_order, geometry = _last_standard_orientation(job_lines)
        atom_orders.append(atom_order)
        geometries.append(geometry)
    if charge_multiplicity[0] != charge_multiplicity[1]:
        raise ValueError(f"{state_label} charge/multiplicity changes between jobs")
    if atom_orders[0] != atom_orders[1]:
        raise ValueError(f"{state_label} atom order changes between optimization and frequency")

    charge, multiplicity = charge_multiplicity[0]
    expected_multiplicity = {"S0": 1, "S1": 1, "T1": 3}[state_label]
    if multiplicity != expected_multiplicity:
        raise ValueError(
            f"{state_label} requires multiplicity {expected_multiplicity}; found {multiplicity}"
        )
    if not any("SCF Done:" in line for line in jobs[1]):
        raise ValueError(f"{state_label} final frequency job lacks an SCF energy")

    target_root = None
    if state_label == "S1":
        for job_number, job_lines in enumerate(jobs, start=1):
            job_text = "\n".join(job_lines)
            state_one_records = re.findall(
                r"Excited\s+State\s+1\s*:\s*([^\n]+)",
                job_text,
                re.IGNORECASE,
            )
            if not state_one_records:
                raise ValueError(f"S1 Gaussian job {job_number} lacks state-1 excitation data")
            if any(
                not re.search(r"\bSinglet(?:[-\s]|$)", record, re.IGNORECASE)
                for record in state_one_records
            ):
                raise ValueError(
                    f"S1 Gaussian job {job_number} contains a non-singlet state-1 record"
                )
            if "Ground to excited state transition electric dipole moments" not in job_text:
                raise ValueError(f"S1 Gaussian job {job_number} lacks a transition-dipole table")
        targeted_roots = []
        for index, line in enumerate(jobs[0]):
            if "This state for optimization" not in line:
                continue
            preceding = "\n".join(jobs[0][max(0, index - 3) : index])
            matches = re.findall(
                r"Excited\s+State\s+(\d+)\s*:", preceding, re.IGNORECASE
            )
            if not matches:
                raise ValueError("S1 optimization marker is not bound to an excited-state root")
            targeted_roots.append(int(matches[-1]))
        if not targeted_roots or any(root != 1 for root in targeted_roots):
            raise ValueError("S1 optimization must remain on target root 1")
        target_root = 1

    fchk_text = fchk_path.read_text(encoding="utf-8", errors="replace")
    fchk_atoms = _fchk_scalar(fchk_text, "Number of atoms")
    fchk_charge = _fchk_scalar(fchk_text, "Charge")
    fchk_multiplicity = _fchk_scalar(fchk_text, "Multiplicity")
    fchk_atomic_numbers = _fchk_atomic_numbers(fchk_text)
    fchk_geometry = _fchk_cartesian_coordinates(fchk_text, fchk_atoms)
    if fchk_atoms != len(atom_orders[-1]) or fchk_atoms != len(fchk_atomic_numbers):
        raise ValueError("Gaussian log and formatted checkpoint atom counts do not match")
    if fchk_atomic_numbers != atom_orders[-1]:
        raise ValueError("Gaussian log and formatted checkpoint Atomic numbers/order do not match")
    if (fchk_charge, fchk_multiplicity) != (charge, multiplicity):
        raise ValueError("Gaussian log and formatted checkpoint charge/multiplicity do not match")
    final_frequency_geometry = geometries[-1]
    for atom_index, (log_xyz, fchk_xyz) in enumerate(
        zip(final_frequency_geometry, fchk_geometry), start=1
    ):
        for axis, (log_value, fchk_value) in zip("xyz", zip(log_xyz, fchk_xyz)):
            if not math.isclose(
                log_value,
                fchk_value,
                rel_tol=0.0,
                abs_tol=FCHK_GEOMETRY_ATOL_ANGSTROM,
            ):
                raise ValueError(
                    "Gaussian final frequency geometry and formatted checkpoint "
                    f"coordinates differ at atom {atom_index} {axis} by more than "
                    f"{FCHK_GEOMETRY_ATOL_ANGSTROM:.1e} Angstrom"
                )

    return {
        "state": state_label,
        "normal_job_count": 2,
        "charge": charge,
        "multiplicity": multiplicity,
        "atomic_numbers": atom_orders[-1],
        "geometry_angstrom": [list(xyz) for xyz in final_frequency_geometry],
        "geometry_tolerance_angstrom": FCHK_GEOMETRY_ATOL_ANGSTROM,
        "target_root": target_root,
    }

def extract_scf_energy(logpath):
    """Extract last SCF Done energy from Gaussian log."""
    energies = []
    with open(logpath) as f:
        for line in f:
            m = re.search(r'SCF Done:\s+E\([^)]+\)\s*=\s*(-?\d+\.\d+)', line)
            if m:
                energies.append(float(m.group(1)))
    if not energies:
        raise ValueError(f"No SCF Done found in {logpath}")
    return energies[-1]

def extract_last_excitation_ev(logpath):
    """Extract the last S1 (state 1) excitation energy (eV) from TDDFT log.
    Returns the adiabatic excitation energy = the very last Excited State 1 value."""
    last_exc_ev = None
    with open(logpath) as f:
        for line in f:
            # Example: "Excited State   1: Singlet-A  1.6143 eV ..."
            match = re.search(
                rf'Excited\s+State\s+(\d+)\s*:.*?({_FLOAT_PATTERN})\s+eV\b',
                line,
                flags=re.IGNORECASE,
            )
            if match and int(match.group(1)) == 1:
                last_exc_ev = _as_float(match.group(2))
    return last_exc_ev


def extract_s1_total_energy(logpath):
    """Return the final S1 total energy in Hartree.

    Gaussian TDDFT logs report the reference-state SCF energy and the
    excitation energy separately. The S1 total energy is therefore the last
    SCF energy plus the last state-1 excitation energy. Missing excitation
    data is an error; silently substituting the SCF energy would label an S0
    reference energy as S1.
    """
    scf_energy = extract_scf_energy(logpath)
    excitation_ev = extract_last_excitation_ev(logpath)
    if excitation_ev is None:
        raise ValueError(f"No S1 state 1 excitation energy found in {logpath}")
    return scf_energy + excitation_ev / HA2EV


def extract_transition_dipole_blocks(logpath):
    """Extract separate Gaussian transition-dipole tables in file order."""
    blocks = []
    current = None
    row_pattern = re.compile(
        rf'^\s*(\d+)\s+({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s+'
        rf'({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})(?:\s|$)'
    )

    with open(logpath) as f:
        for line in f:
            if 'Ground to excited state transition electric dipole moments' in line:
                if current:
                    blocks.append(current)
                current = []
                continue

            if current is None:
                continue

            match = row_pattern.match(line)
            if match:
                current.append({
                    'state': int(match.group(1)),
                    'X': _as_float(match.group(2)),
                    'Y': _as_float(match.group(3)),
                    'Z': _as_float(match.group(4)),
                    'DipS': _as_float(match.group(5)),
                    'Osc': _as_float(match.group(6)),
                })
            elif current and (
                not line.strip()
                or 'Ground to excited state transition velocity' in line
            ):
                blocks.append(current)
                current = None

    if current:
        blocks.append(current)
    return blocks

def extract_transition_dipoles(logpath):
    """Extract transition electric dipole moments for all excited states.
    Returns list of dicts with state, X, Y, Z, Dip.S., Osc."""
    return [
        transition
        for block in extract_transition_dipole_blocks(logpath)
        for transition in block
    ]


def _with_transition_magnitude(transition):
    """Add transition-dipole magnitudes in atomic units and Debye."""
    coordinates = [transition.get(axis) for axis in ('X', 'Y', 'Z')]
    if all(value is not None for value in coordinates):
        magnitude_au = math.sqrt(sum(value * value for value in coordinates))
    elif transition.get('DipS') is not None and transition['DipS'] >= 0:
        magnitude_au = math.sqrt(transition['DipS'])
    else:
        raise ValueError("Transition dipole has neither XYZ components nor valid DipS")

    enriched = dict(transition)
    enriched['magnitude_au'] = magnitude_au
    enriched['magnitude_debye'] = magnitude_au * AU2DEBYE
    return enriched


def extract_state1_transition_endpoints(logpath):
    """Return state-1 transitions from the first and last TD dipole blocks."""
    blocks = extract_transition_dipole_blocks(logpath)
    if not blocks:
        raise ValueError(f"No TD transition-dipole blocks found in {logpath}")

    endpoints = []
    for label, block in (('first', blocks[0]), ('last', blocks[-1])):
        transition = next((item for item in block if item['state'] == 1), None)
        if transition is None:
            raise ValueError(
                f"No state 1 transition in {label} TD dipole block of {logpath}"
            )
        endpoints.append(_with_transition_magnitude(transition))
    return tuple(endpoints)

def extract_dipole_moment(logpath):
    """Extract last dipole moment from Gaussian log."""
    dipoles = []
    with open(logpath) as f:
        for line in f:
            m = re.match(r'\s+X=\s*(-?\d+\.\d+)\s+Y=\s*(-?\d+\.\d+)\s+Z=\s*(-?\d+\.\d+)\s+Tot=\s*(-?\d+\.\d+)', line)
            if m:
                dipoles.append({
                    'X': float(m.group(1)), 'Y': float(m.group(2)),
                    'Z': float(m.group(3)), 'Tot': float(m.group(4))
                })
    return dipoles if dipoles else None

def count_normal_terminations(logpath):
    """Count exact Gaussian normal-termination lines."""
    count = 0
    with open(logpath) as f:
        for line in f:
            if _GAUSSIAN_NORMAL_LINE.fullmatch(line.rstrip("\n")):
                count += 1
    return count

def generate_spec_tvcf_input(params, output_path='momap_spec.inp'):
    """Generate MOMAP spec_tvcf input file."""
    missing = [key for key in ('Ead', 'EDMA', 'EDME') if key not in params]
    if missing:
        raise ValueError(
            'Missing molecule-specific spectrum parameters: ' + ', '.join(missing)
        )
    for key in ('Ead', 'EDMA', 'EDME'):
        value = params[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f'{key} must be a finite numeric value')
        if not math.isfinite(value) or value < 0:
            raise ValueError(f'{key} must be finite and non-negative')

    output_path = Path(output_path).expanduser().resolve()
    if output_path.exists() or output_path.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"Output parent does not exist: {output_path.parent}")

    template = f"""do_spec_tvcf_ft   = 1
do_spec_tvcf_spec = 1

&spec_tvcf
  DUSHIN   = .t.
  HERZ     = .t.
  Temp     = {params.get('temperature', 300)} K
  tmax     = {params.get('tmax', 5000)} fs
  dt       = {params.get('dt', 0.001)} fs
  Ead      = {params['Ead']:.8f} au
  EDMA     = {params['EDMA']:.6f} debye
  EDME     = {params['EDME']:.6f} debye
  FreqScale  = 1.0
  DSFile     = "{params.get('dsfile', 'evc.cart.dat')}"
  Emax       = {params.get('emax', 0.3)} au
  dE         = {params.get('de', 0.00001)} au
  logFile    = "spec.tvcf.log"
  FtFile     = "spec.tvcf.ft.dat"
  FoFile     = "spec.tvcf.fo.dat"
  FoSFile    = "spec.tvcf.spec.dat"
/
"""
    with open(output_path, 'x', encoding='utf-8') as f:
        f.write(template)
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Extract MOMAP parameters from Gaussian logs')
    parser.add_argument('--s0', required=True, help='S0 ground state Gaussian log')
    parser.add_argument('--s1', required=True, help='S1 excited state Gaussian log')
    parser.add_argument('--t1', help='T1 triplet state Gaussian log (optional)')
    parser.add_argument('--s0-fchk', help='S0 formatted checkpoint (default: paired .fchk)')
    parser.add_argument('--s1-fchk', help='S1 formatted checkpoint (default: paired .fchk)')
    parser.add_argument('--t1-fchk', help='T1 formatted checkpoint (default: paired .fchk)')
    parser.add_argument('--output', '-o', default='momap_spec.inp', help='Output MOMAP input file')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--json', action='store_true', help='Output JSON instead of input file')
    import json as json_mod

    args = parser.parse_args()

    def paired_fchk(log_value, explicit_value):
        return (
            Path(explicit_value).expanduser()
            if explicit_value
            else Path(log_value).expanduser().with_suffix('.fchk')
        )

    validate_gaussian_log_contract(
        args.s0, paired_fchk(args.s0, args.s0_fchk), 'S0'
    )
    validate_gaussian_log_contract(
        args.s1, paired_fchk(args.s1, args.s1_fchk), 'S1'
    )
    if args.t1:
        validate_gaussian_log_contract(
            args.t1, paired_fchk(args.t1, args.t1_fchk), 'T1'
        )

    # Extract SCF energies + adiabatic excitation
    E_S0 = extract_scf_energy(args.s0)
    E_S1_scf = extract_scf_energy(args.s1)  # S0 reference at S1 geometry
    exc_ev = extract_last_excitation_ev(args.s1)
    E_S1_total = extract_s1_total_energy(args.s1)
    Ead = E_S1_total - E_S0

    # Extract transition dipoles from S1 log
    # EDMA: from first TDM block (S0 geometry → vertical absorption)
    # EDME: from last TDM block (S1 minimum → emission)
    absorption, emission = extract_state1_transition_endpoints(args.s1)
    EDMA = absorption['magnitude_debye']
    EDME = emission['magnitude_debye']

    # Extract oscillator strengths
    f_abs = absorption['Osc']
    f_emi = emission['Osc']

    # Count normal terminations
    nts_s0 = count_normal_terminations(args.s0)
    nts_s1 = count_normal_terminations(args.s1)
    nts_t1 = count_normal_terminations(args.t1) if args.t1 else None

    params = {
        'Ead': Ead,
        'EDMA': EDMA,
        'EDME': EDME,
        'f_abs': f_abs,
        'f_emi': f_emi,
        'E_S0': E_S0,
        'E_S1': E_S1_total,
        'E_S1_scf': E_S1_scf,
        'E_exc_ev': exc_ev,
        'temperature': args.temperature,
        'dsfile': 'evc.cart.dat',
        'tmax': 5000,
        'dt': 0.001,
        'emax': 0.3,
        'de': 0.00001,
        'nts_s0': nts_s0,
        'nts_s1': nts_s1,
        'nts_t1': nts_t1,
    }

    if args.json:
        print(json_mod.dumps(params, indent=2))
    else:
        path = generate_spec_tvcf_input(params, args.output)
        print(f"✅ MOMAP spec_tvcf input → {path}")
        print(f"   Ead  = {Ead:.8f} au ({Ead*27.2114:.4f} eV)")
        print(f"   EDMA = {EDMA:.4f} debye (f={f_abs:.4f})")
        print(f"   EDME = {EDME:.4f} debye (f={f_emi:.4f})")
        print(f"   S0 NTs: {nts_s0}, S1 NTs: {nts_s1}" + (f", T1 NTs: {nts_t1}" if nts_t1 else ""))

if __name__ == '__main__':
    main()
