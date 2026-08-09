# MOMAP 2024A Quick Reference

## Load

```bash
source "${MOMAP_ENV:-/opt/MOMAP-2024A/env.sh}"  # Installation environment
module load "${MOMAP_MODULE:-momap/2024A-openmpi}"  # Site-specific alternative
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

## MPI Patch (for OpenMPI 3.x)

```bash
cp $MOMAP_BIN/momap momap_patched
sed -i 's/-machinefile/--hostfile/g' momap_patched
echo "localhost slots=4" > nodefile
python momap_patched -i momap.inp
```

## Spectrum Data Columns

```
#1 Energy(Ha)  2 Energy(eV)  3 Wavenumber(cm-1)  4 Wavelength(nm)  5 FC_abs  6 FC_emi  7 FC_emi_intensity
```

Plot column 4 vs column 6 (emission) or 4 vs 5 (absorption).

## Bundled Tests

```bash
ls "$MOMAP_ROOT/tests/"
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
  Ead     = 0.094 au      Hso   = 116.9 cm-1  ← from PySOC
  DSFile  = "evc.cart.dat"
  Emax    = 0.3 au
/
```

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
