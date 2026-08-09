#!/usr/bin/env python3
"""
momap-oled — Extended OLED tool wrappers for isc_tvcf, ic_tvcf, pysoc, sumstat, transport.
"""
import os
import sys
import subprocess
import argparse
import math
import re
from pathlib import Path

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
_RATE_NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
_RATE_UNIT = r's(?:\^-1|-1|⁻¹)'
_FATAL_DIAGNOSTIC_PATTERN = re.compile(
    r'\bfatal\b|segmentation fault|error termination|\baborted\b|'
    r'traceback \(most recent call last\)|forrtl:',
    flags=re.IGNORECASE,
)
_SUPPORTED_MOMAP_BUILD = '2024A'

# ─── isc_tvcf ───────────────────────────────────────────────────────────
def _finite_positive(name, value):
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a finite positive number')
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be a finite positive number') from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a finite positive number')
    return value


def generate_isc_input(
    evc_dat,
    ead_au,
    hso_cm1,
    output='momap_isc.inp',
    temp=298,
    tmax=5000,
):
    """Generate isc_tvcf input for phosphorescence spectrum + k_ISC."""
    evc_path = Path(evc_dat).expanduser().resolve(strict=True)
    if not evc_path.is_file() or evc_path.stat().st_size <= 0:
        raise ValueError(f'evc_dat must be a non-empty regular file: {evc_path}')
    ead_au = _finite_positive('ead_au', ead_au)
    hso_cm1 = _finite_positive('hso_cm1', hso_cm1)
    temp = _finite_positive('temp', temp)
    tmax = _finite_positive('tmax', tmax)
    output_path = Path(output).expanduser().resolve()
    if output_path.exists():
        raise FileExistsError(f'Refusing to overwrite existing output: {output_path}')
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f'Output parent does not exist: {output_path.parent}')
    content = f"""do_isc_tvcf_ft   = 1
do_isc_tvcf_spec = 1

&isc_tvcf
  DUSHIN  = .t.
  Temp    = {temp} K
  tmax    = {tmax} fs
  dt      = 0.001 fs
  Ead     = {ead_au:.8f} au
  Hso     = {hso_cm1:.6f} cm-1
  DSFile  = "{evc_path}"
  Emax    = 0.3 au
  logFile = "isc.tvcf.log"
  FtFile  = "isc.tvcf.ft.dat"
  FoFile  = "isc.tvcf.fo.dat"
/
"""
    with open(output_path, 'x', encoding='utf-8') as f:
        f.write(content)
    return output_path


def parse_isc_output(fo_file='isc.tvcf.fo.dat', *, expected_build=None):
    """Parse finite ISC/RISC rates using the local 2024A output contract."""
    path = Path(fo_file).expanduser().resolve(strict=True)
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f'ISC rate output must be non-empty: {path}')
    text = path.read_text(encoding='utf-8', errors='replace')
    fatal_match = _FATAL_DIAGNOSTIC_PATTERN.search(text)
    if fatal_match:
        raise ValueError(
            f'Fatal diagnostic in ISC rate output: {fatal_match.group(0)!r}'
        )

    explicit_pattern = re.compile(
        rf'^\s*(ISC|RISC)\b.*?\brate\s+is\s+({_RATE_NUMBER})\s*'
        rf'{_RATE_UNIT}\s*$',
        flags=re.IGNORECASE,
    )
    unlabeled_pattern = re.compile(
        rf'^\s*rate\s+is\s+({_RATE_NUMBER})\s*{_RATE_UNIT}\s*$',
        flags=re.IGNORECASE,
    )
    rate_trigger = re.compile(r'\brate\s+is\b', flags=re.IGNORECASE)
    explicit_records = {'ISC': [], 'RISC': []}
    unlabeled_records = []

    def parse_value(raw_value, label):
        value = float(raw_value.replace('D', 'E').replace('d', 'e'))
        if not math.isfinite(value) or value < 0:
            raise ValueError(f'{label} rate must be finite and non-negative')
        return value

    for line_number, line in enumerate(text.splitlines(), start=1):
        if not rate_trigger.search(line):
            continue
        explicit_match = explicit_pattern.fullmatch(line)
        if explicit_match:
            label = explicit_match.group(1).upper()
            explicit_records[label].append(
                parse_value(explicit_match.group(2), label)
            )
            continue
        unlabeled_match = unlabeled_pattern.fullmatch(line)
        if unlabeled_match:
            unlabeled_records.append(
                parse_value(unlabeled_match.group(1), 'unlabeled')
            )
            continue
        raise ValueError(
            f'Malformed or unitless ISC rate record at line {line_number}'
        )

    explicit_count = sum(len(records) for records in explicit_records.values())
    if explicit_count:
        if unlabeled_records:
            raise ValueError(
                'ISC rate output mixes explicit labels with unlabeled fallback records'
            )
        for label in ('ISC', 'RISC'):
            count = len(explicit_records[label])
            if count != 1:
                raise ValueError(
                    f'Explicit {label} rate must appear exactly once; found {count}'
                )
        return {
            'k_ISC_s-1': explicit_records['ISC'][0],
            'k_RISC_s-1': explicit_records['RISC'][0],
            'rate_parse_contract': 'explicit_ISC_RISC_labels_exactly_once',
        }

    if len(unlabeled_records) != 2:
        raise ValueError(
            'MOMAP 2024A unlabeled fallback requires exactly two rate records '
            f'with s^-1 units; found {len(unlabeled_records)}'
        )
    if expected_build != _SUPPORTED_MOMAP_BUILD:
        raise ValueError(
            'The unlabeled ISC/RISC ordering is supported only when '
            'expected_build is explicitly MOMAP 2024A'
        )
    return {
        'k_ISC_s-1': unlabeled_records[0],
        'k_RISC_s-1': unlabeled_records[1],
        'rate_parse_contract': 'MOMAP_2024A_ordered_first_ISC_second_RISC',
    }


# ─── ic_tvcf ────────────────────────────────────────────────────────────
def generate_ic_input(evc_dat, nac_coul, ead_au, output='momap_ic.inp', temp=300, tmax=655):
    """Generate ic_tvcf input for internal conversion rate."""
    content = f"""do_ic_tvcf_ft   = 1
do_ic_tvcf_spec = 1

&ic_tvcf
  DUSHIN   = .t.
  Temp     = {temp} K
  tmax     = {tmax} fs
  dt       = 0.01 fs
  Ead      = {ead_au:.8f} au
  DSFile   = "{evc_dat}"
  CoulFile = "{nac_coul}"
  logFile  = "ic.tvcf.log"
  FtFile   = "ic.tvcf.ft.dat"
  FoFile   = "ic.tvcf.fo.dat"
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── pysoc ──────────────────────────────────────────────────────────────
def generate_pysoc_input(qm_input_com, output='momap_pysoc.inp',
                         n_singlets=4, n_triplets=4,
                         qc_exe='g16', qc_ppn=8):
    """Generate PySOC input for spin-orbit coupling calculation."""
    content = f"""do_pysoc = 1

&pysoc
  sched_type = local
  qc_exe     = {qc_exe}
  qc_ppn     = {qc_ppn}
  pysoc_QM_code = 'gauss_tddft'
  pysoc_QM_input_file = {qm_input_com}
  n_excited_singlets = {n_singlets}
  n_excited_triplets = {n_triplets}
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── spec_sums ──────────────────────────────────────────────────────────
def generate_sums_input(evc_dat, ead_au, dipole_abs, dipole_emi,
                        output='momap_sums.inp', fwhm=500, maxvib=10):
    """Generate spec_sums input — all rates + quantum yield in one pass."""
    content = f"""do_spec_sums = 1

&spec_sums
  DSFile     = "{evc_dat}"
  Ead        = {ead_au:.8f} au
  dipole_abs = {dipole_abs:.6f} debye
  dipole_emi = {dipole_emi:.6f} debye
  maxvib     = {maxvib}
  if_cal_ic  = .t.
  FWHM       = {fwhm} cm-1
  flog       = "spec.sums.log"
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


def parse_sums_output(logfile='spec.sums.log'):
    """Parse spec_sums output for all rates."""
    with open(logfile) as f:
        text = f.read()
    import re
    k_emi = re.search(r'Emission rate is\s+([\d.E+-]+)\s+s-1', text)
    return {
        'k_r_s1': float(k_emi.group(1)) if k_emi else None,
    }


# ─── transport ──────────────────────────────────────────────────────────
def generate_transport_input(cif_file, output='momap_transport.inp',
                             basis='b3lyp STO-3g', temp=300,
                             qc_exe='g16', ratetype='marcus',
                             nsimu=2000, tsimu=1000):
    """Generate transport input for charge carrier mobility."""
    content = f"""&transport
  do_transport_prepare              = 1
  do_transport_submit_HL_job        = 1
  do_transport_get_transferintegral = 1
  do_transport_submit_RE_job        = 1
  do_transport_get_re_evc           = 1
  do_transport_run_MC               = 1
  do_transport_get_mob_MC           = 1
  do_transport_gather_momap_data    = 1

  sched_type     = local
  queue_name     = localhost
  compute_engine = 1
  qc_exe         = {qc_exe}
  basis_name     = {basis}
  basis_name_re  = {basis}
  qc_memory      = 4096
  qc_nodes       = 1
  qc_ppn         = 2
  temp           = {temp}
  ratetype       = {ratetype}
  lat_cutoff     = 4
  nsimu          = {nsimu}
  tsimu          = {tsimu}
  crystal        = {cif_file}
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output


# ─── CLI ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='MOMAP OLED tool wrappers')
    sub = parser.add_subparsers(dest='cmd')

    p = sub.add_parser('isc', help='Generate isc_tvcf input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument(
        '--hso', type=float, required=True,
        help='Computed or measured spin-orbit coupling in cm^-1',
    )
    p.add_argument('--temp', type=float, default=298, help='Temperature in K')
    p.add_argument('--tmax', type=float, default=5000, help='Correlation time in fs')
    p.add_argument('-o', default='momap_isc.inp')

    p = sub.add_parser('ic', help='Generate ic_tvcf input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--nac', default='evc.cart.nac')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument('-o', default='momap_ic.inp')

    p = sub.add_parser('pysoc', help='Generate PySOC input')
    p.add_argument('--com', required=True, help='Gaussian TDDFT input file')
    p.add_argument('--n-singlets', type=int, default=4)
    p.add_argument('--n-triplets', type=int, default=4)
    p.add_argument('-o', default='momap_pysoc.inp')

    p = sub.add_parser('sums', help='Generate spec_sums input')
    p.add_argument('--evc-dat', default='evc.cart.dat')
    p.add_argument('--ead', type=float, required=True)
    p.add_argument('--dipole-abs', type=float, required=True)
    p.add_argument('--dipole-emi', type=float, required=True)
    p.add_argument('-o', default='momap_sums.inp')

    p = sub.add_parser('transport', help='Generate transport input')
    p.add_argument('--cif', required=True, help='Crystal structure .cif file')
    p.add_argument('-o', default='momap_transport.inp')

    args = parser.parse_args()
    if not args.cmd:
        parser.print_help()
        return 1

    if args.cmd == 'isc':
        path = generate_isc_input(
            args.evc_dat,
            args.ead,
            args.hso,
            args.o,
            temp=args.temp,
            tmax=args.tmax,
        )
        print(f"✅ isc_tvcf input → {path}")
    elif args.cmd == 'ic':
        path = generate_ic_input(args.evc_dat, args.nac, args.ead, args.o)
        print(f"✅ ic_tvcf input → {path}")
    elif args.cmd == 'pysoc':
        path = generate_pysoc_input(args.com, args.o, args.n_singlets, args.n_triplets)
        print(f"✅ pysoc input → {path}")
    elif args.cmd == 'sums':
        path = generate_sums_input(args.evc_dat, args.ead, args.dipole_abs, args.dipole_emi, args.o)
        print(f"✅ spec_sums input → {path}")
    elif args.cmd == 'transport':
        path = generate_transport_input(args.cif, args.o)
        print(f"✅ transport input → {path}")

    return 0

if __name__ == '__main__':
    sys.exit(main())
