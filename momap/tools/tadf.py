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
import subprocess
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from extract import (
    extract_scf_energy,
    extract_last_excitation_ev,
    extract_s1_total_energy,
    extract_state1_transition_endpoints,
    count_normal_terminations,
    generate_spec_tvcf_input,
)
from runner import (
    patch_momap_for_mpi3,
    stage_gaussian_inputs,
    create_nodefile,
    MOMAP_ROOT,
)

AU2DEBYE = 2.541746
HA2EV = 27.2114
_MOL_ID_PATTERN = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$')

def run_momap_in_dir(inputfile, workdir):
    """Run patched MOMAP in workdir."""
    workdir = Path(workdir).expanduser().resolve()
    inputfile = Path(inputfile).expanduser().resolve(strict=True)
    patched = patch_momap_for_mpi3()
    create_nodefile(workdir)
    
    env = os.environ.copy()
    env['MOMAP_ROOT'] = MOMAP_ROOT
    env['MOMAP_LICENSE'] = os.path.join(MOMAP_ROOT, 'license', 'hzwtech.lic')
    
    cmd = [sys.executable, patched, '-i', str(inputfile)]
    print(f"  🚀 {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(workdir), env=env,
                          capture_output=False)
    return result.returncode == 0

def parse_evc_out(evc_out):
    """Extract reorganization energy and mode count from evc.out."""
    with open(evc_out) as f:
        content = f.read()
    
    m = re.search(r'(\d+)\s+# num of atoms', content)
    natoms = int(m.group(1)) if m else 0
    m = re.search(r'(\d+)\s+# num of modes', content)
    nmodes = int(m.group(1)) if m else 0
    
    return {'natoms': natoms, 'nmodes': nmodes}

def parse_spec_output(spec_file):
    """Parse spec.tvcf.spec.dat for peak wavelength."""
    if not Path(spec_file).exists():
        return {}
    
    try:
        data = []
        with open(spec_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    data.append([float(x) for x in parts])
        
        if not data:
            return {}
        
        emi = [row[5] for row in data]  # FC_emi
        wavelengths = [row[3] for row in data]
        
        # Find peaks
        peaks = []
        for i in range(1, len(emi)-1):
            if emi[i] > emi[i-1] and emi[i] > emi[i+1] and emi[i] > 0:
                peaks.append({'wavelength': wavelengths[i], 'intensity': emi[i]})
        
        peaks.sort(key=lambda x: x['intensity'], reverse=True)
        
        max_idx = max(range(len(emi)), key=lambda i: emi[i])
        
        return {
            'peak_wavelength': wavelengths[max_idx],
            'peak_intensity': emi[max_idx],
            'top_peaks': peaks[:5],
            'data_points': len(data),
        }
    except Exception as e:
        return {'error': str(e)}

def process_molecule(
    mol_id,
    s0_log,
    s1_log,
    t1_log,
    output_dir,
    temperature=300,
    hso_cm1=None,
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

    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    mol_dir = output_root / mol_id
    mol_dir.mkdir(parents=True, exist_ok=True)

    if hso_cm1 is not None:
        try:
            hso_cm1 = float(hso_cm1)
        except (TypeError, ValueError):
            results['error'] = 'hso_cm1 must be a finite positive number'
            return results
        if not math.isfinite(hso_cm1) or hso_cm1 <= 0:
            results['error'] = 'hso_cm1 must be a finite positive number'
            return results

    print(f"\n{'='*60}")
    print(f"🧬 Processing {mol_id}")
    print(f"{'='*60}")

    source_logs = {'S0': s0_log, 'S1': s1_log, 'T1': t1_log}
    staged_logs = {}
    try:
        for label, path in source_logs.items():
            if path is None:
                print(f"  ⚠️  {label}: no file provided")
                continue
            source = Path(path).expanduser().resolve(strict=True)
            nts = count_normal_terminations(source)
            print(f"  📄 {label}: {source} ({nts} Normal terminations)")
            if nts < 2:
                print(f"     ⚠️  Warning: expected 2 NTs (opt+freq), got {nts}")
            staged_log, _ = stage_gaussian_inputs(
                source, mol_dir, target_stem=label.lower()
            )
            staged_logs[label] = staged_log
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
    results.update({
        'E_S0': E_S0,
        'E_S1': E_S1,
        'E_S1_scf': E_S1_scf,
        'E_exc_ev': E_exc_ev,
        'E_T1': E_T1,
        'Ead_S1_S0': (E_S1 - E_S0) * HA2EV,
        'delta_EST_eV': delta_EST,
        'EDMA': absorption['magnitude_debye'],
        'EDME': emission['magnitude_debye'],
        'f_abs': absorption['Osc'],
        'f_emi': emission['Osc'],
    })

    print(f"\n  📊 Extracted parameters:")
    print(f"     E(S0)   = {E_S0:.8f} au")
    print(f"     E(S1)   = {E_S1:.8f} au")
    if E_T1 is not None:
        print(f"     E(T1)   = {E_T1:.8f} au")
        print(f"     ΔE_ST   = {delta_EST:.4f} eV")
    print(f"     EDMA    = {results['EDMA']:.4f} debye")
    print(f"     EDME    = {results['EDME']:.4f} debye")

    print(f"\n  📐 Step 1: EVC (S1→S0)")
    evc_s1_input = mol_dir / 'momap_evc_s1.inp'
    with open(evc_s1_input, 'w') as f:
        f.write(f"""do_evc = 1

&evc
  ffreq(1) = "{s0_staged.name}"
  ffreq(2) = "{s1_staged.name}"
  sort_mode = 1
/
""")
    if not run_momap_in_dir(evc_s1_input, mol_dir):
        results['error'] = 'EVC (S1 to S0) failed'
        print(f"  ❌ {results['error']}")
        return results
    evc_s1_info = parse_evc_out(mol_dir / 'evc.out')
    print(f"  ✅ {evc_s1_info['natoms']} atoms, {evc_s1_info['nmodes']} modes")

    print(f"\n  🌈 Step 2: Spectrum (S1→S0)")
    spec_params = {
        'temperature': temperature,
        'Ead': E_S1 - E_S0,
        'EDMA': results['EDMA'],
        'EDME': results['EDME'],
        'dsfile': 'evc.cart.dat',
    }
    spec_input = mol_dir / 'momap_spec.inp'
    generate_spec_tvcf_input(spec_params, spec_input)
    if run_momap_in_dir(spec_input, mol_dir):
        spec_info = parse_spec_output(mol_dir / 'spec.tvcf.spec.dat')
        results['spectrum'] = spec_info
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
            results['error'] = 'Spectrum output contains no valid positive-intensity peak'
            print(f"  ❌ {results['error']}")
            return results

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
    else:
        results['spectrum'] = {'status': 'failed'}
        results['error'] = 'Spectrum calculation failed'
        print(f"  ❌ {results['error']}")
        return results

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
        evc_isc_input = mol_dir / 'momap_evc_isc.inp'
        with open(evc_isc_input, 'w') as f:
            f.write(f"""do_evc = 1

&evc
  ffreq(1) = "{s1_staged.name}"
  ffreq(2) = "{t1_staged.name}"
  sort_mode = 1
/
""")
        if not run_momap_in_dir(evc_isc_input, mol_dir):
            results['isc'] = {'status': 'failed', 'reason': 'EVC (S1 to T1) failed'}
            print("  ❌ EVC (S1→T1) failed")
        else:
            isc_input = mol_dir / 'momap_isc.inp'
            with open(isc_input, 'w') as f:
                f.write(f"""do_isc_tvcf_ft   = 1
do_isc_tvcf_spec = 1

&isc_tvcf
  DUSHIN  = .t.
  Temp    = {temperature} K
  tmax    = 5000 fs
  dt      = 0.001 fs
  Ead     = {abs(E_S1 - E_T1):.8f} au
  Hso     = {hso_cm1:.8g} cm-1
  DSFile  = "evc.cart.dat"
  Emax    = 0.3 au
  logFile = "isc.tvcf.log"
  FtFile  = "isc.tvcf.ft.dat"
  FoFile  = "isc.tvcf.fo.dat"
/
""")
            ok_isc = run_momap_in_dir(isc_input, mol_dir)
            results['isc'] = {
                'status': 'computed' if ok_isc else 'failed',
                'hso_cm1': hso_cm1,
            }
            print(
                f"  {'✅' if ok_isc else '⚠️ '} ISC TVCF "
                f"{'complete' if ok_isc else 'had issues'}"
            )

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
    args = parser.parse_args()

    try:
        results = process_molecule(
            mol_id=args.mol_id,
            s0_log=args.s0,
            s1_log=args.s1,
            t1_log=args.t1,
            output_dir=args.output,
            temperature=args.temperature,
            hso_cm1=args.hso_cm1,
        )
    except Exception as exc:
        results = {
            'mol_id': args.mol_id,
            'success': False,
            'error': f'Unexpected pipeline failure: {exc}',
        }

    if args.json_output:
        json_output = Path(args.json_output).expanduser()
        json_output.parent.mkdir(parents=True, exist_ok=True)
        with open(json_output, 'w') as f:
            json.dump(results, f, indent=2, default=str)
            f.write('\n')
    
    if results.get('success'):
        return 0
    else:
        print(f"\n❌ Processing failed for {args.mol_id}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
