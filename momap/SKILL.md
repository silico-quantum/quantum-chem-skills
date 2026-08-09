---
name: momap
description: Use when preparing, running, or interpreting MOMAP workflows for molecular photophysics, nonradiative rates, spectra, spin-orbit coupling, or charge transport.
license: MIT
compatibility: Requires a licensed MOMAP installation; selected workflows also require Gaussian outputs, MPI, PySOC, or a Slurm environment.
---

# MOMAP Workflow Skill

**MOMAP (Molecular Materials Property Prediction Package)** predicts molecular-material photophysical and transport properties.
Version 2024A (2.3.7), Hongzhiwei Technology / Z.G. Shuai Group

## Environment Setup

```bash
# Source the licensed installation's environment script.
source "${MOMAP_ENV:-/opt/MOMAP-2024A/env.sh}"

# Alternatively, use the site-specific environment module when available.
module load "${MOMAP_MODULE:-momap/2024A-openmpi}"
```

**Key env vars set by env.sh:**
- `MOMAP_ROOT` points to the installation root.
- `MOMAP_BIN=$MOMAP_ROOT/bin`
- `MOMAP_MPI_BIN=$MOMAP_ROOT/bin/openmpi/bin`
- `MOMAP_LICENSE` points to the licensed installation's license file.
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

### Step 3: ISC Rate (isc_tvcf)

```
do_isc_tvcf_ft   = 1       # Fourier transform
do_isc_tvcf_spec = 1       # Phosphorescence spectrum

&isc_tvcf
  DUSHIN   = .t.
  Temp     = 298 K
  tmax     = 5000 fs
  dt       = 0.001 fs
  Ead      = 0.094 au     # S1-T1 adiabatic gap
  Hso      = 116.9 cm-1   # ⚠️ Spin-orbit coupling (from PySOC or exp.)
  DSFile   = "evc.cart.dat"
  Emax     = 0.3 au
  logFile  = "isc.tvcf.log"
  FtFile   = "isc.tvcf.ft.dat"
  FoFile   = "isc.tvcf.fo.dat"
/
```

**Output:** ISC rate (k_ISC), phosphorescence spectrum.

### Step 4: Internal Conversion (ic_tvcf)

```
do_ic_tvcf_ft   = 1
do_ic_tvcf_spec = 1

&ic_tvcf
  DUSHIN   = .t.
  Temp     = 300 K
  tmax     = 655 fs
  dt       = 0.01 fs
  Ead      = 0.094 au
  DSFile   = "evc.cart.dat"
  CoulFile = "evc.cart.nac"  # ⚠️ Non-adiabatic coupling
  logFile  = "ic.tvcf.log"
  FtFile   = "ic.tvcf.ft.dat"
/
```

**Output:** IC rate (k_IC) — non-radiative decay S₁→S₀.

### Step 5: Spin-Orbit Coupling (PySOC)

Auto-generates SOC matrix elements from TDDFT:

```
do_pysoc = 1

&pysoc
  sched_type = local
  qc_exe     = g16
  qc_ppn     = 8
  pysoc_QM_code = 'gauss_tddft'
  n_excited_singlets = 4
  n_excited_triplets = 4
/
```

**Output:** Hso values (cm⁻¹) for S₁→T₁, S₁→T₂, ... → feeds into isc_tvcf.

### Step 6: Full Rate Summary (spec_sums)

Sum-over-states method — computes k_r, k_IC, k_ISC simultaneously:

```
do_spec_sums = 1

&spec_sums
  DSFile     = "evc.cart.dat"
  Ead        = 0.094 au
  dipole_abs = 0.092 debye
  dipole_emi = 0.441 debye
  maxvib     = 10
  if_cal_ic  = .t.
  FWHM       = 500 cm-1
/
```

**Output:** All rates + quantum yield Φ in one pass.

### Step 7: Charge Transport (transport)

Carrier mobility for OLED charge balance:

```
&transport
  do_transport_prepare              = 1
  do_transport_get_transferintegral = 1
  do_transport_get_re_evc           = 1
  do_transport_run_MC               = 1
  do_transport_get_mob_MC           = 1
  compute_engine  = 1            # Gaussian
  qc_exe          = g16
  basis_name      = b3lyp STO-3g
  temp            = 300
  ratetype        = marcus       # marcus or quanta
  lat_cutoff      = 4
  nsimu           = 2000
  tsimu           = 1000         # ns
  crystal         = molecule.cif
/
```

**Output:** Hole/electron mobility (cm²/V·s), transfer integrals, reorganization energies.

## OLED TADF Design — Complete Workflow

```
       Gaussian                          MOMAP
  ┌──────────────┐              ┌─────────────────────┐
  │ S0 opt+freq  │──────┐       │ EVC (S₁→S₀)         │
  │ S1 opt+freq  │      ├──────→│ EVC (S₁→T₁)         │
  │ T1 opt+freq  │──────┘       │                     │
  │ NACT calc    │──→ NACME ──→│ ic_tvcf   → k_IC    │
  │ PySOC calc   │──→ Hso  ───→│ isc_tvcf  → k_ISC   │
  │ TDDFT (vert) │──→ EDMA/EDME│ spec_tvcf → k_r     │
  └──────────────┘              │ spec_sums → Φ (QY)  │
                                │ transport → μ_h, μ_e│
                                └─────────────────────┘
```

**Key OLED figures of merit:**
| Quantity | MOMAP Module | Required Input |
|----------|-------------|---------------|
| k_r (radiative) | spec_tvcf | EVC + EDMA/EDME |
| k_IC (internal conv.) | ic_tvcf | EVC + NACME |
| k_ISC (S₁→T₁) | isc_tvcf | EVC + Hso (SOC) |
| k_RISC (T₁→S₁) | isc_tvcf | k_ISC × exp(-ΔE_ST/kT) |
| Φ (quantum yield) | spec_sums | all above |
| μ_h (hole mobility) | transport | Transfer integral + λ |
| μ_e (electron mob.) | transport | Transfer integral + λ |

## Extracting Parameters from Gaussian Output

### Ead (adiabatic excitation energy)
```
# From S0 log:
grep "SCF Done" s0.log | tail -1    →  E_SCF_S0

# From S1 log at the S1 minimum:
grep "SCF Done" s1.log | tail -1    →  E_SCF_reference
grep "Excited State   1:" s1.log | tail -1  →  E_exc,S1 in eV

# E_S1,total = E_SCF_reference + E_exc,S1 / 27.2114
# Ead = E_S1,total - E_SCF_S0  (in Hartree)
# Do not substitute E_SCF_reference when the state-1 excitation is missing.
```

### EDMA (absorption transition dipole moment)
From S0 geometry TDDFT output:
```
grep -A5 "transition electric dipole" s1.log | head -10
#  state    X        Y        Z       Dip.S.    Osc.
#    1   0.2169   0.2932   0.0000   0.1330    0.0079
# sqrt(X² + Y² + Z²) × 2.541746 → debye
# Equivalently, use sqrt(Dip.S.) × 2.541746 if XYZ is unavailable.
```

### EDME (emission transition dipole moment)
From S1 minimum geometry TDDFT (last Excited State block):
```
# Same as EDMA but from the last TDDFT block in s1.log
# Use the vector norm for state 1 at the S1-optimized geometry
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
- `$MOMAP_ROOT/tests/azulene/gaussian/ref/` — azulene S0/S1
- `$MOMAP_ROOT/tests/porphine/gaussian_g16/` — porphine
- `$MOMAP_ROOT/tests/Irppy3/` — Ir(ppy)₃

Repository guidance:
- [Worked examples](EXAMPLES.md)
- [Command quick reference](QUICKREF.md)

## Citations

- Y. Niu et al., *Molecular Physics*, **2018**, doi: 10.1080/00268976.2017.1402966
- Q. Peng et al., *J. Am. Chem. Soc.*, **2007**, 129, 9333-9339
- Z. Shuai, Q. Peng, *Phys. Rep.*, **2014**, 537, 123
- Z. Shuai, Q. Peng, *Nat. Sci. Rev.*, **2017**, 4, 224

Official: http://www.momap.cn

## Toolkit (momap/tools/)

Five Python CLI tools for automated MOMAP workflows:

### momap-extract — Gaussian → MOMAP params
```bash
python3 tools/extract.py --s0 mol-s0.log --s1 mol-s1.log [--t1 mol-t1.log]
# Output: spec_tvcf input + JSON params (Ead, EDMA, EDME)
# Auto-handles TDDFT: Ead = SCF_S1 + E_exc_adiabatic - SCF_S0
```

### momap-run — One-command MOMAP wrapper  
```bash
python3 tools/runner.py momap.inp [--slurm] [--nprocs 4]
# Auto: MPI --hostfile patch, formchk (.chk→.fchk), Slurm submit
```

### momap-tadf — Full TADF pipeline
```bash
python3 tools/tadf.py mol_id --s0 s0.log --s1 s1.log --t1 t1.log \
  --json-output result.json
# Runs: EVC(S1→S0) → spec_tvcf → summary
# Output: spectrum, peak λ, blue window check, ΔE_ST, ISC=not_computed

# Add ISC only when a computed or measured SOC value is available:
python3 tools/tadf.py mol_id --s0 s0.log --s1 s1.log --t1 t1.log \
  --hso-cm1 12.5 --json-output result.json
```

The TADF tool stages each Gaussian log with its matching `.fchk` in the MOMAP
work directory. It computes the S1 total energy as the last SCF reference
energy plus the last state-1 excitation energy and fails if either quantity is
missing. It uses state 1 from the first and last TD transition-dipole tables
for absorption and emission, respectively, and converts the XYZ vector norm
from atomic units to Debye. ISC is not run without an explicit `--hso-cm1`;
no placeholder spin-orbit coupling is assumed.

### momap-plot — Spectrum visualization
```bash
python3 tools/plot.py spec.tvcf.spec.dat -D --blue 450 490
# Pillow fallback if no matplotlib. -D = save to Desktop.
# Prints peak analysis with blue window highlighting.
```

### momap-oled — OLED extended modules (v3)
```bash
python3 tools/oled.py isc  --ead 0.094 --hso 116.9          # isc_tvcf → k_ISC
python3 tools/oled.py ic   --ead 0.094 --nac evc.cart.nac   # ic_tvcf → k_IC
python3 tools/oled.py pysoc --com mol.com                    # PySOC → Hso values
python3 tools/oled.py sums --ead 0.094 --dipole-abs 0.09 --dipole-emi 0.44  # Φ, all rates
python3 tools/oled.py transport --cif crystal.cif            # μ_h, μ_e
# Each generates a ready-to-run MOMAP input file.
```
