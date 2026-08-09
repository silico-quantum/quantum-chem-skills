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

_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _as_float(value):
    """Parse Gaussian-style floating-point text, including Fortran D exponents."""
    return float(value.replace('D', 'E').replace('d', 'e'))

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
    """Count Normal termination occurrences."""
    count = 0
    with open(logpath) as f:
        for line in f:
            if 'Normal termination' in line:
                count += 1
    return count

def generate_spec_tvcf_input(params, output_path='momap_spec.inp'):
    """Generate MOMAP spec_tvcf input file."""
    template = f"""do_spec_tvcf_ft   = 1
do_spec_tvcf_spec = 1

&spec_tvcf
  DUSHIN   = .t.
  HERZ     = .t.
  Temp     = {params.get('temperature', 300)} K
  tmax     = {params.get('tmax', 5000)} fs
  dt       = {params.get('dt', 0.001)} fs
  Ead      = {params.get('Ead', 0.07509):.8f} au
  EDMA     = {params.get('EDMA', 0.92694):.6f} debye
  EDME     = {params.get('EDME', 0.64751):.6f} debye
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
    with open(output_path, 'w') as f:
        f.write(template)
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Extract MOMAP parameters from Gaussian logs')
    parser.add_argument('--s0', required=True, help='S0 ground state Gaussian log')
    parser.add_argument('--s1', required=True, help='S1 excited state Gaussian log')
    parser.add_argument('--t1', help='T1 triplet state Gaussian log (optional)')
    parser.add_argument('--output', '-o', default='momap_spec.inp', help='Output MOMAP input file')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--json', action='store_true', help='Output JSON instead of input file')
    import json as json_mod

    args = parser.parse_args()

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
