# MOMAP Skill (v2 — 2024A Verified)

**MOMAP (Molecular Materials Property Prediction Package)** — 分子材料光物理/传输性质预测
Version 2024A (2.3.7), Hongzhiwei Technology / Z.G. Shuai Group

## Server Setup

```bash
# On marcus / marcus2
source /opt/MOMAP-2024A/env.sh
# or
module load momap/2024A-openmpi
```

**Key env vars set by env.sh:**
- `MOMAP_ROOT=/opt/MOMAP-2024A`
- `MOMAP_BIN=/opt/MOMAP-2024A/bin`
- `MOMAP_MPI_BIN=/opt/MOMAP-2024A/bin/openmpi/bin`
- `MOMAP_LICENSE=/opt/MOMAP-2024A/license/hzwtech.lic`
- `LD_LIBRARY_PATH` includes `.../openmpi/lib` and `.../lib`

**MPI compatibility note:** MOMAP ships OpenMPI 1.x but the system OpenMPI is 3.1.4.
The `-machinefile` flag is not supported in OpenMPI 3.x (use `--hostfile`).
→ Copy and patch `momap` script: `sed 's/-machinefile/--hostfile/g'`
→ Run `python momap_patched -i momap.inp`
→ Nodefile format for 3.x: `localhost slots=4`

## Real Input Format (MOMAP 2024A)

**No `&control` block, no `jobtype`.** Each calculation type is a top-level flag + a named block.

### Step 1: Electron-Vibration Coupling (EVC)

**Prerequisites:** Gaussian `.log` **AND** `.fchk` files for both states.
Generate `.fchk` with: `formchk job.chk job.fchk`

```
do_evc = 1

&evc
  ffreq(1) = "mol-s0.log"    # Ground state freq (S0)
  ffreq(2) = "mol-s1.log"    # Excited state freq (S1/T1)
  ftdipd   = "mol-s1.log"    # Optional: transition dipole from this file
  sort_mode = 1              # Mode sorting (1=default)
/
```

```bash
momap -i momap.inp
```

**Output files:**
| File | Content |
|------|---------|
| `evc.cart.dat` | Duschinsky matrix + Huang-Rhys factors (Cartesian) |
| `evc.dint.dat` | Duschinsky matrix (internal coords) |
| `evc.cart.abs` | Absorption data |
| `evc.dint.abs` | Absorption (internal coords) |
| `evc.out` | Full log with frequencies, mode pairing, reorganization energies |
| `evc.vib*.xyz` | Vibration mode visualizations |
| `evc.dx.x.xyz` | Displacement vectors |
| `evc.dx.x.com` | Displaced structure (Gaussian input) |

### Step 2: Spectrum (spec_tvcf)

```
do_spec_tvcf_ft   = 1       # Fourier transform of time correlation function
do_spec_tvcf_spec = 1       # Generate spectrum

&spec_tvcf
  DUSHIN   = .t.            # Use Duschinsky rotation
  HERZ     = .t.            # Herzberg-Teller effect (optional)
  Temp     = 300 K
  tmax     = 5000 fs        # Correlation time
  dt       = 0.001 fs        # Time step
  Ead      = 0.07509 au     # ⚠️ Adiabatic excitation energy (S1_min - S0_min)
  EDMA     = 0.92694 debye  # ⚠️ Absorption transition dipole (vertical)
  EDME     = 0.64751 debye  # ⚠️ Emission transition dipole (from S1 min)
  FreqScale = 1.0
  DSFile   = "evc.cart.dat" # From EVC step
  Emax     = 0.3 au         # Max energy for output
  dE       = 0.00001 au     # Energy resolution
  logFile  = "spec.tvcf.log"
  FtFile   = "spec.tvcf.ft.dat"
  FoFile   = "spec.tvcf.fo.dat"
  FoSFile  = "spec.tvcf.spec.dat"
/
```

```bash
momap -i momap.inp
```

**Output (spec.tvcf.spec.dat):**
```
#1Energy(Hartree) 2Energy(eV) 3Wavenumber(cm-1) 4Wavelength(nm) 5FC_abs 6FC_emi 7FC_emi_intensity
```

### Step 3: ISC / IC Rates

```
do_isc = 1

&isc
  ffreq(1) = "mol-s1.log"   # S1 frequency
  ffreq(2) = "mol-t1.log"   # T1 frequency
  Temp     = 300 K
/
```

## Extracting Parameters from Gaussian Output

### Ead (adiabatic excitation energy)
```
# From S0 log:
grep "SCF Done" s0.log | tail -1    →  E_SCF_S0

# From S1 log (last SCF Done = S1 min):
grep "SCF Done" s1.log | tail -1    →  E_SCF_S1

# Ead = E_SCF_S1 - E_SCF_S0  (in Hartree)
```

### EDMA (absorption transition dipole moment)
From S0 geometry TDDFT output:
```
grep -A5 "transition electric dipole" s1.log | head -10
#  state    X        Y        Z       Dip.S.    Osc.
#    1   0.2169   0.2932   0.0000   0.1330    0.0079
# Dip.S. (au) → multiply by 2.5417 to convert to debye
```

### EDME (emission transition dipole moment)
From S1 minimum geometry TDDFT (last Excited State block):
```
# Same as EDMA but from the last TDDFT block in s1.log
# Use Dip.S. of state 1 at S1-optimized geometry
```

## Example: Azulene (bundled test)

Test directory: `/opt/MOMAP-2024A/tests/azulene/gaussian/`

```bash
mkdir ~/momap_azulene && cd ~/momap_azulene
# Copy ref log + fchk files
cp /opt/MOMAP-2024A/tests/azulene/gaussian/ref/azulene-s0.* .
cp /opt/MOMAP-2024A/tests/azulene/gaussian/ref/azulene-s1.* .

# Step 1: EVC
cat > momap_evc.inp << 'EOF'
do_evc = 1
&evc
  ffreq(1) = "azulene-s0.log"
  ffreq(2) = "azulene-s1.log"
/
EOF
source /opt/MOMAP-2024A/env.sh
momap -i momap_evc.inp

# Step 2: Spectrum (copy ref params)
cp /opt/MOMAP-2024A/tests/azulene/kr/momap.inp momap_spec.inp
momap -i momap_spec.inp
```

## TADF Workflow Integration

For TADF molecules, the MOMAP pipeline is:

```
Gaussian S0 opt+freq  →  s0.log + s0.fchk
Gaussian T1 opt+freq  →  t1.log + t1.fchk
Gaussian S1 TDDFT     →  s1.log + s1.fchk (vertical + adiabatic)
       ↓
   EVC (S0 vs S1)     →  evc_s1.cart.dat
   EVC (S1 vs T1)     →  evc_t1.cart.dat
       ↓
   spec_tvcf (S1→S0)  →  fluorescence spectrum
   ISC (S1→T1)        →  k_ISC rate
       ↓
   Quantum yield: Φ = k_r / (k_r + k_IC + k_ISC)
```

## Troubleshooting

### `mpiexec: Error: unknown option "-machinefile"`
MOMAP was compiled with OpenMPI 1.x, system has 3.x.
**Fix:** Copy + patch `momap` script:
```bash
cp $MOMAP_BIN/momap momap_patched
sed -i 's/-machinefile/--hostfile/g' momap_patched
echo "localhost slots=4" > nodefile
python momap_patched -i momap.inp
```

### `Can not find Gaussian fchk file!`
MOMAP needs **both** `.log` and `.fchk` files in the same directory.
```bash
formchk job.chk job.fchk
```

### spec_tvcf killed (OOM)
Use MPI on a compute node with more memory:
```bash
#SBATCH --mem=64G
mpirun -np 4 momap -i momap.inp
```

### Reference files
All test examples with pre-computed Gaussian logs:
- `/opt/MOMAP-2024A/tests/azulene/gaussian/ref/` — azulene S0/S1
- `/opt/MOMAP-2024A/tests/porphine/gaussian_g16/` — porphine
- `/opt/MOMAP-2024A/tests/Irppy3/` — Ir(ppy)₃

## Citations

- Y. Niu et al., *Molecular Physics*, **2018**, doi: 10.1080/00268976.2017.1402966
- Q. Peng et al., *J. Am. Chem. Soc.*, **2007**, 129, 9333-9339
- Z. Shuai, Q. Peng, *Phys. Rep.*, **2014**, 537, 123
- Z. Shuai, Q. Peng, *Nat. Sci. Rev.*, **2017**, 4, 224

Official: http://www.momap.cn
