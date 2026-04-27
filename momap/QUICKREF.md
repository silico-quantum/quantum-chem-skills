# MOMAP Quick Reference (v2 — 2024A Verified)

## Load

```bash
source /opt/MOMAP-2024A/env.sh           # full env
module load momap/2024A-openmpi          # or via module
```

## Input Format Summary

| Calculation | Top Flag | Namelist Block | Key Files |
|------------|----------|---------------|-----------|
| EVC (electron-vibration coupling) | `do_evc = 1` | `&evc` | `.log` + `.fchk` × 2 |
| Spectrum (TVCF) | `do_spec_tvcf_ft = 1` + `do_spec_tvcf_spec = 1` | `&spec_tvcf` | `evc.cart.dat` |
| ISC rate | `do_isc = 1` | `&isc` | S1 + T1 `.log`/`.fchk` |
| IC rate | `do_ic = 1` | `&ic` | S0 + S1 `.log`/`.fchk` + NACME |
| Transport | `do_transport = 1` | `&transport` | Dimer logs |

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
  Ead     = 0.075 au      ← from S0/S1 SCF difference
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
# Ead
E_S0=$(grep "SCF Done" s0.log | tail -1 | awk '{print $5}')
E_S1=$(grep "SCF Done" s1.log | tail -1 | awk '{print $5}')
# Ead = E_S1 - E_S0

# Transition dipole (au → debye: ×2.5417)
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
ls /opt/MOMAP-2024A/tests/
# azulene/  Irppy3/  porphine/  transport/  numfreq/
# Each: gaussian/  kr/  kic/  kisc/  evc/  ...
```
