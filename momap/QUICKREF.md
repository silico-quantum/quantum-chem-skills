# MOMAP 2024A Quick Reference

## Load

```bash
: "${MOMAP_ENV:?Set MOMAP_ENV to the licensed installation env.sh}"
source "$MOMAP_ENV"
# Or load the site-specific module documented by the local administrator.
MOMAP_BUILD=2024A
MOMAP_LAUNCHER=$(command -v momap)
test -f "$MOMAP_LAUNCHER" && test ! -L "$MOMAP_LAUNCHER" || exit 1
MOMAP_LAUNCHER_SHA256=$(python3 -c \
  'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' \
  "$MOMAP_LAUNCHER")
: "${MOMAP_VERSION_BANNER:?Set the exact verified 2024A full-line runtime banner}"
```

## Input Format Summary

| Calculation | Top Flag | Namelist Block | Key Files | OLED Use |
|------------|----------|---------------|-----------|----------|
| EVC | `do_evc = 1` | `&evc` | `.log` + `.fchk` × 2 | Duschinsky + HR |
| Fluorescence | `do_spec_tvcf_ft = 1` + `do_spec_tvcf_spec = 1` | `&spec_tvcf` | `evc.cart.dat` | k_r, emission λ |
| Internal Conv. | `do_ic_tvcf_ft = 1` + `do_ic_tvcf_spec = 1` | `&ic_tvcf` | EVC + NACME | k_IC |
| ISC / Phosphor. | `do_isc_tvcf_ft = 1` + `do_isc_tvcf_spec = 1` | `&isc_tvcf` | EVC + Hso | k_ISC, Φ_phos |
| Sum-over-states | `do_spec_sums = 1` | `&spec_sums` | EVC + all dipoles | Φ, all rates |
| Spin-orbit coupling | `do_pysoc = 1` | `&pysoc` | Gaussian TDDFT com | Hso values |
| Charge transport | `do_transport_* = 1` | `&transport` | Dimer logs + .cif | μ_h, μ_e |

## EVC Quick Template

```
do_evc = 1
&evc
  ffreq(1) = "s0.log"
  ffreq(2) = "s1.log"
/
```

Output: `evc.cart.dat`, `evc.dint.dat`, `evc.out`

Automated acceptance requires fresh, non-empty `evc.out` and `evc.cart.dat`,
positive atom/mode counts, and `Normal finish of evc calculation` in `evc.out`.

## spec_tvcf Quick Template

```
do_spec_tvcf_ft   = 1
do_spec_tvcf_spec = 1
&spec_tvcf
  DUSHIN  = .t.
  Temp    = 300 K
  tmax    = 5000 fs
  dt      = 0.001 fs
  Ead     = 0.075 au      ← E_S1,total - E_S0
  EDMA    = 0.93 debye    ← from S0→S1 vertical TDM
  EDME    = 0.65 debye    ← from S1→S0 emission TDM
  DSFile  = "evc.cart.dat"
  Emax    = 0.3 au
  dE      = 0.00001 au
  FoSFile = "spec.tvcf.spec.dat"
/
```

## Extract from Gaussian

```bash
# Ead: S1 total energy is the final SCF reference plus the final S1 excitation
E_S0=$(grep "SCF Done" s0.log | tail -1 | awk '{print $5}')
E_S1_SCF=$(grep "SCF Done" s1.log | tail -1 | awk '{print $5}')
E_S1_EXC_EV=$(grep "Excited State[[:space:]]*1:" s1.log | tail -1 | awk '{for (i=1;i<=NF;i++) if ($i=="eV") print $(i-1)}')
# E_S1_total = E_S1_SCF + E_S1_EXC_EV / 27.2114
# Ead = E_S1_total - E_S0; fail if either S1 term is unavailable

# State-1 transition dipole: norm(X,Y,Z) × 2.541746 au→debye.
# Use the first table for absorption and the last table for emission.
grep -A5 "transition electric dipole" s1.log

# Generate fchk
test ! -e job.fchk && test ! -L job.fchk || exit 1
formchk job.chk job.fchk
```

## Output Files

| File | Content | Use |
|------|---------|-----|
| `evc.cart.dat` | Duschinsky + HR factors | Input to spec_tvcf / ISC |
| `evc.dint.dat` | Internal-coord Duschinsky | Verification |
| `evc.out` | Frequencies, reorg energies, mode pairing | Analysis |
| `spec.tvcf.spec.dat` | Spectrum: eV, nm, FC_abs, FC_emi | Plot |
| `spec.tvcf.ft.dat` | Time correlation function FT | Debug |
| `evc.vib1.xyz` | Vibration mode XYZ | Visualization |

## Fail-closed TADF stage layout

`tools/tadf.py` refuses an existing molecule output directory and isolates each
producer of fixed MOMAP filenames:

```text
<output>/<molecule>/
  evc_s1_s0/   stage_receipt.json
  spectrum/    stage_receipt.json
  evc_s1_t1/   stage_receipt.json
  isc/         stage_receipt.json
```

The receipts bind accepted inputs and outputs to SHA-256 hashes, sizes, and
timestamps. They also record build, exact version banner, original launcher
SHA-256, patched launcher SHA-256, and the MPI replacement contract.
`evc.cart.dat` is copied forward only after its producing EVC stage passes
acceptance.

Completion markers for the repository's MOMAP 2024A-style fixture are:

| Stage | File | Required marker |
|---|---|---|
| EVC | `evc.out` | `Normal finish of evc calculation` |
| spectrum | `spec.tvcf.log` | `Normal finish of spec_tvcf calculation` |
| ISC | `isc.tvcf.log` | `Normal finish of isc_tvcf calculation` |

Other releases fail closed until their local output is reviewed and tested.

## MPI Launcher Compatibility

Do not modify the licensed launcher or keep a permanent hand-edited copy. The
bundled runner creates a private mode-0700 temporary launcher, applies the
`-machinefile` to `--hostfile` compatibility substitution there, and leaves the
original untouched. Give it a fresh work directory because MOMAP uses fixed
output names:

```bash
SOURCE_INPUT="$PWD/momap.inp"
RUN_DIR="$PWD/momap_mpi_run_001"
test -s "$SOURCE_INPUT"
test ! -e "$RUN_DIR" && test ! -L "$RUN_DIR" || exit 1
mkdir "$RUN_DIR"
python3 "${MOMAP_SKILL_DIR:?Set MOMAP_SKILL_DIR}/tools/runner.py" \
  "$SOURCE_INPUT" --workdir "$RUN_DIR" \
  --expected-build "$MOMAP_BUILD" \
  --expected-launcher-sha256 "$MOMAP_LAUNCHER_SHA256" \
  --expected-version-banner "$MOMAP_VERSION_BANNER"
```

## Spectrum Data Columns

```
#1 Energy(Ha)  2 Energy(eV)  3 Wavenumber(cm-1)  4 Wavelength(nm)  5 FC_abs  6 FC_emi  7 FC_emi_intensity
```

Plot column 4 vs column 6 (emission) or 4 vs 5 (absorption).

The automated parser requires explicit `expected_build=2024A` plus verified
launcher/run evidence before accepting exactly seven numeric columns on
every non-comment data row, at least three finite rows, a strictly monotonic eV
axis, an oppositely ordered wavelength axis, and `E(eV) * wavelength(nm)` within
1% of `hc`. Tail diagnostics, partial rows, extra columns, non-finite values,
and fatal text fail the spectrum stage.

## Bundled Tests

```bash
ls "${MOMAP_ROOT:?Set MOMAP_ROOT}/tests/"
# azulene/  Irppy3/  porphine/  transport/  numfreq/  pysoc/
# Each: gaussian/  kr/  kic/  kisc/  evc/  sumstat/
```

## OLED Quick Templates

### isc_tvcf (phosphorescence + k_ISC)
```
do_isc_tvcf_ft   = 1
do_isc_tvcf_spec = 1
&isc_tvcf
  DUSHIN  = .t.           Temp  = 298 K
  tmax    = 5000 fs       dt    = 0.001 fs
  Ead     = 0.094 au      Hso   = 116.9 cm-1  ← named computed/measured source
  DSFile  = "evc.cart.dat"
  Emax    = 0.3 au
/
```

For the bundled pipeline, `Ead = E_S1,total - E_T1,SCF` is signed and expected
to be positive. `E_T1,SCF` is the final Gaussian `SCF Done` total electronic
energy at the T1 geometry. Never use `abs(E_S1 - E_T1)`. A requested ISC stage
is accepted only when fresh non-empty files contain the completion marker and
two finite rates with explicit `s^-1` units. The local unlabeled two-line 2024A
format is enabled only with an explicit verified 2024A build and is recorded as
`first=ISC, second=RISC`. Labelled output must contain
exactly one ISC and one RISC record; duplicates, conflicts, mixed labelled and
unlabelled records, or extra fallback rate lines are rejected.

### ic_tvcf (internal conversion)
```
do_ic_tvcf_ft   = 1
do_ic_tvcf_spec = 1
&ic_tvcf
  DUSHIN  = .t.           Temp  = 300 K
  tmax    = 655 fs        dt    = 0.01 fs
  Ead     = 0.094 au      DSFile  = "evc.cart.dat"
  CoulFile = "evc.cart.nac"  ← NACME required
/
```

### spec_sums (all rates + Φ)
```
do_spec_sums = 1
&spec_sums
  DSFile     = "evc.cart.dat"
  Ead        = 0.094 au
  dipole_abs = 0.09 debye   dipole_emi = 0.44 debye
  maxvib     = 10           if_cal_ic  = .t.
  FWHM       = 500 cm-1
/
```

### pysoc (spin-orbit coupling)
```
do_pysoc = 1
&pysoc
  sched_type = local        qc_exe = g16
  pysoc_QM_code = 'gauss_tddft'
  pysoc_QM_input_file = mol.com
  n_excited_singlets = 4    n_excited_triplets = 4
/
```

### transport (carrier mobility)
```
&transport
  do_transport_* = 1       # all stages enabled
  compute_engine = 1       qc_exe = g16
  basis_name = b3lyp STO-3g
  temp = 300               ratetype = marcus
  crystal = molecule.cif
/
```
