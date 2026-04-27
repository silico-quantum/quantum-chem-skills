---
name: gaussian
version: 2.0.0
description: Gaussian 16 quantum chemistry software. Electronic structure calculations, geometry optimization, frequency analysis, TDDFT excited states, and more. Based on official Gaussian.com documentation.
homepage: https://gaussian.com
metadata:
  category: computational_chemistry
  tags: [quantum_chemistry, DFT, TDDFT, frequency, optimization, excited_states, ab_initio]
---

# Gaussian 16 — Quantum Chemistry Software

## Overview

Gaussian is the industry-standard quantum chemistry software for electronic structure calculations. It supports DFT (all major functionals), HF, MPn, CCSD(T), CASSCF, TDDFT, geometry optimization, frequency analysis, reaction pathways, and much more.

**Installed on:** marcus HPC cluster (remote server access via marcus2)

---

## Connection

```bash
ssh -p 8722 openclaw@124.16.75.110   # Login to marcus2
ssh marcus                               # Jump to compute node
module load gaussian/g16.c01-avx2       # Load Gaussian 16 Rev. C.01
```

---

## Input File Format (.gjf)

```
%chk=filename.chk          # Checkpoint file
%mem=60GB                  # Memory allocation
%nproc=64                  # Number of CPU cores
# route_section             # Method, basis set, job type, options

Title
Charge SpinMultiplicity     # Molecular charge and spin (0 1 = neutral singlet)
Cartesian or Z-Matrix coordinates
```

---

## Link0 Commands

| Command | Description |
|---------|-------------|
| `%chk=file.chk` | Checkpoint file path |
| `%mem=60GB` | Memory (any unit: MB, GB) |
| `%nproc=64` | CPU cores |
| `%nprocshared=32` | Shared memory parallelism |
| `%oldchk=file.chk` | Continue from existing checkpoint |
| `%save` | Force checkpoint saving |

---

## Job Types

| Keyword | Description | Common Options |
|---------|-------------|----------------|
| `SP` | Single-point energy | Default calculation |
| `Opt` | Geometry optimization | `Opt=CalcFC`, `Opt=TS`, `Opt=Tight`, `Opt=MaxStep=N` |
| `Freq` | Frequency analysis | `Freq=Raman`, `Freq=HPC` |
| `Opt Freq` | Optimize then compute frequencies | Verify minima/transition states |
| `Scan` | Potential energy surface scan | `Scan=ModRedundant`, `Scan=Bond`, `Scan=Angle` |
| `IRC` | Intrinsic reaction coordinate | Follow reaction path from TS |
| `Polar` | Polarizability | `Polar=EnOnly` |
| `ADMP` | Car-Parrinello molecular dynamics | |
| `BOMD` | Born-Oppenheimer MD | |
| `stable` | Wavefunction stability test | |

---

## Methods

### Hartree-Fock
```
# HF/6-31G*
```
- `HF` — Restricted (RHF), `UHF` — Unrestricted, `ROHF` — Restricted open-shell

### DFT Functionals

| Functional | Type | HF Exchange | Notes |
|------------|------|-------------|-------|
| `B3LYP` | Hybrid GGA | 20% | Most popular, good all-round |
| `PBE0` | Hybrid GGA | 25% | Based on PBE, reliable |
| `M06-2X` | Meta-GGA hybrid | 54% | Good for non-covalent, thermochemistry |
| `M06` | Meta-GGA hybrid | 27% | Minnesota functional |
| `CAM-B3LYP` | Range-separated hybrid | Variable | Long-range corrected |
| `WB97X-D` | Range-separated + D3 | Variable | Dispersion-corrected |
| `WB97X-D3BJ` | Range-separated + D3(BJ) | Variable | Better dispersion |
| `LC-WPBE` | Long-range corrected | 100% at long range | Charge-transfer states |
| `wPBE` | Long-range corrected | 100% at long range | Same as LC-WPBE |
| `TPSSH` | Hybrid meta-GGA | 10% | |
| `BP86` | GGA | 0% | |
| `PBE` | GGA | 0% | |
| `BLYP` | GGA | 0% | |
| `SVWN5` | LDA | 0% | |
| `HSEHJS` | Range-separated | 50% at short range | Screened Coulomb |

### Wave Function Methods

| Method | Description |
|--------|-------------|
| `MP4` | 4th-order Møller-Plesset perturbation |
| `MP3` | 3rd-order MP |
| `MP2` | 2nd-order MP (also `MP2=Full`) |
| `CCSD` | Coupled-cluster singles + doubles |
| `CCSD(T)` | CCSD with perturbative triples |
| `CCSD=T` | CCSD with full triples (non-perturbative) |
| `CCD`, `CCSD`, `CID`, `CISD` | Other coupled-cluster variants |
| `CASSCF(n,m)` | Complete active space SCF (n electrons, m orbitals) |
| `CASPT2` | Multi-state second-order perturbation |
| `EOMCCSD` | Equation-of-motion CCSD (excited states) |
| `SAC-CI` | Symmetry-adapted cluster CI |
| `GVB` | Generalized valence bond |
| `CIS` | Configuration interaction singles |
| `CIS(D)` | CIS with perturbative doubles |

### Semi-Empirical
```
# CNDO, INDO, MINDO3, MNDO, AM1, PM3, ZINDO
```

---

## Basis Sets

### Pople Split-Valence
| Basis | Description |
|-------|-------------|
| `STO-3G` | Minimal basis (fast) |
| `3-21G` | Split-valence, lightest |
| `6-21G` | Pople split-valence |
| `6-31G` | Double-zeta, standard for organics |
| `6-31G*` | With polarization on heavy atoms |
| `6-31G**` | With polarization on all atoms |
| `6-31+G*` | With diffusion on heavy atoms |
| `6-31++G**` | With diffusion on all |
| `6-311G` | Triple-zeta |
| `6-311G*` | Triple-zeta + polarization |
| `6-311++G(2d,2p)` | Highly extended |

### Correlation-Consistent (cc-pVnZ)
```
cc-pVDZ, cc-pVTZ, cc-pVQZ, cc-pV5Z, cc-pV6Z
aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ    # Augmented with diffuse functions
```

### Def2 / Karlsruhe
```
def2-SVP, def2-TZVP, def2-QZVP
def2-TZVPP — with polarization functions
```

### Effective Core Potentials (ECPs)
| ECP | Description |
|-----|-------------|
| `LANL2DZ` | Los Alamos ECP + DZ basis |
| `LANL2DZdp` | LANL2DZ with polarization |
| `Stuttgart` | Stuttgart-Cologne ECPs |
| `CRENBL` | Relativistic ECP |

### Other
| Basis | Description |
|-------|-------------|
| `cc-pCVnZ` | Core-valence (for core excitations) |
| `pc-n` | Polarization consistent basis |
| `pcseg-n` | Segmented polarization consistent |
| `N07D` | Density-fitting basis sets |

---

## TDDFT — Excited States

```
# td=(nstates=10,root=1,singlets) b3lyp/6-31G*
# td=(nstates=10,triplets) b3lyp/6-31G*   # Triplet states only
```

**Options:**
- `nstates=N` — Number of excited states (default 10)
- `root=N` — Optimize/persist on state N
- `singlets` — Singlet excitations (default)
- `triplets` — Triplet excitations
- `alltrans` — Include all triplet-singlet transitions
- `50-50` — 50% CI coefficients for each determinant
- `50-50TDDFT` — Special for degenerate states

**For vertical emission:**
```
# td=(nstates=5,singlets) b3lyp/6-31G*     # Excitation
# b3lyp/6-31G* td=(root=1) geom=check     # Emission from S1
```

---

## Solvent Effects (SCRF)

```
# b3lyp/6-31G* scrf=(smd,solvent=acetonitrile)
# b3lyp/6-31G* scrf=(pcm,solvent=water)
```

**SCRF Models:**
| Model | Keyword | Notes |
|-------|---------|-------|
| PCM | `SCRF=PCM` | Polarizable continuum |
| SMD | `SCRF=SMD` | Solvent model based on density (G16 default) |
| CPMF | `SCRF=CM5` | |
| HCTH | `SCRF=HCTH` | |
| SWIG | `SCRF=SWIG` | |

**Common Solvents:** `water`, `ethanol`, `methanol`, `acetonitrile`, `dmso`, `dmf`, `chloroform`, `toluene`, `benzene`, `diethylether`, `thf`, `acetone`

---

## Geometry Optimization Options

```
# opt=calcfc           # Compute force constants initially (default for TS)
# opt=ts               # Transition state search
# opt=tight            # Tight convergence criteria
# opt=gdiis            # Use GDIIS for optimization
# opt=modredundant     # Modredundant coordinates
# opt=maxstep=n        # Maximum step size (default 0.3 bohr)
# opt=readvariance     # Read variance data
# freq=hpc             # HPC-mode frequency scaling
```

**ModRedundant examples (in input):**
```
C 1 2 R    # Define bond between atoms 1-2
D 1 2 3 4 180.0 F    # Scan dihedral angle 1-2-3-4 to 180°
```

---

## Frequency Options

```
# freq=hindered rotor   # Include hindered rotor corrections
# freq=raman             # Raman activity
# freq=check            # Read frequency data from checkpoint
# freq=(hinderedrot,temperature=298.15,pressure=1.0)
```

---

## Population Analysis

```
# pop=full              # Full AO/MO analysis
# pop=full,nbomo        # Include NBOs
# pop=nbo               # Natural Bond Orbital analysis
# pop=npa               # Natural Population Analysis
# pop=mk                # Merz-Kollman ESP charges
# pop=chelpg            # CHELPG ESP charges
# pop=(bondorder)       # Wiberg/Mayer bond order
# pop=density          # Output density matrices
```

---

## NMR Properties

```
# nmr=cpa              # Continuous set of gauge transformations
# nmr=spin-spin        # Spin-spin coupling constants
# spinnuclear=         # Specify magnetic nuclei
# mag=                  # Magnetic properties
```

---

## SCF Convergence

```
# scf=qc               # Quadratic convergence (most robust)
# scf=xqc              # Extra quadratic convergence (try QC if this fails)
# scf=diis             # DIIS (default)
# scf=maxcyc=n         # Max SCF cycles (default 128)
# scf=conver=n         # Convergence criterion (10^-n)
# scf=tight            # Tight SCF convergence
# scf=novariational    # Non-variational MOs
# scf=direct           # Direct SCF (on-the-fly integral recomputation)
# scf=fermi            # Fermi-VWN5 for metals
```

---

## CBS Extrapolation

```
# cbs=QB3              # CBS-QB3 composite method
# cbs=APF              # CBS-APF method
# cbs=4s               # 4-parameter CBS extrapolation
# CBSExtrapolate=2/3   # Two-point extrapolation n=2,3
```

**Composite Methods:**
| Method | Description |
|--------|-------------|
| `CBS-QB3` | CBS-QB3 composite |
| `CBS-APF` | CBS-APF composite |
| `W1U` | W1 ultra-fine |
| `W1` | W1 method |
| `W2` | W2 method |
| `W3` | W3 method |
| `W4` | W4 method |
| `CBS-4s` | CBS with 4-parameter extrapolation |
| `G1` | G1 composite |
| `G2` | G2 composite |
| `G3` | G3 composite |
| `G4` | G4 composite |

---

## IRC (Intrinsic Reaction Coordinate)

```
# irc=kalida            # Use Kupidonov-Khailenko algorithm
# irc=gs2               # Gonzalez-Schlegel second-order algorithm (default)
# irc=nproc=n           # Points per process in IRC
# ircmax=n             # Maximum points in IRC path
```

---

## External Potentials / Basis

```
# ExtraBasis=file        # Additional basis functions
# Gen                   # User-defined basis in input
# GenECP               # User-defined ECP in input
# Pseudo=read          # Read ECP from checkpoint
# densityfit           # Use density fitting (DF) for efficiency
# nodensityfit         # Disable density fitting
```

---

## Wavefunction Analysis

```
# Density=current       # Use current density
# guess=modify,alter    # Modify/alter initial MOs
# guess=hop            # Hückel initial guess
# guess=read           # Read MOs from checkpoint
# stable=opt           # Optimize unstable wavefunction
# prop=flex           # Flexible property output
```

---

## Stabilty Analysis

```
# stable                # Test wavefunction stability
# stable=opt            # Optimize to stable wavefunction
```

---

## Common Task Templates

### Ground-State Optimization + Frequency
```
%chk=opt_freq.chk
%mem=60GB
%nproc=64
# opt freq b3lyp/6-31G**

Molecule optimization
0 1
C  0.0  0.0  0.0
...
```

### TDDFT Excited States
```
%chk=tddft.chk
%mem=60GB
%nproc=64
# td=(nstates=10,singlets) b3lyp/6-31G*

S1 vertical excitations
0 1
C  0.0  0.0  0.0
...
```

### S1 Optimization + Vertical Emission
```
%chk=s1_opt.chk
%mem=60GB
%nproc=64
# td=(nstates=5,root=1) b3lyp/6-31G* opt

S1 optimization
0 1
... (ground-state optimized geometry)
```

```
%chk=emission.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* td=(root=1) geom=check

S1 vertical emission
0 1
... (S1 optimized geometry)
```

### Solvent Effect (SMD)
```
%chk=smd.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* scrf=(smd,solvent=acetonitrile)

SMD solvent model
0 1
...
```

### CBS-QB3 Composite
```
%chk=cbs-qb3.chk
%mem=60GB
%nproc=64
# cbs-qb3

CBS-QB3 calculation
0 1
...
```

### CASSCF
```
%chk=cas.chk
%mem=60GB
%nproc=64
# casSCF(10,8)/6-31G*

CAS(10,8) — 10 electrons, 8 orbitals in active space
0 1
...
```

### ModRedundant Scan
```
%chk=scan.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* opt=modredundant

Dihedral scan
0 1
C 0 0 0
C 1 2 1.4
...
1 2 3 4 90.0 S 20 5.0   # Scan dihedral 1-2-3-4, 20 steps, 5° each
```

---

## Useful IOp Values

| IOp | Description |
|-----|-------------|
| `IOp(3/76=055000)` | Scale HF exchange in hybrid functionals |
| `IOp(3/77=055000)` | Scale correlation |
| `IOp(3/78=...)` | DFT integration grid |
| `IOp(5/33=...)` | Frozen core approximation |
| `IOp(5/34=...)` | Use spherical harmonics |
| `IOp(5/41=...)` | Population analysis options |

---

## Converting Checkpoint Files

```bash
formchk filename.chk filename.fchk      # CHK → FCHK
# or
unfchk filename.chk filename.fchk       # Alternative
```

---

## Output Files

| Extension | Description |
|-----------|-------------|
| `.log` | Main output file |
| `.chk` | Binary checkpoint file |
| `.fchk` | Formatted ASCII checkpoint |
| `.out` | (sometimes used) |
| `.cub` | Cube file for visualization |
| `.int` | Two-electron integrals |
| `.sbn` | Scratch binary file |

---

## Official Keyword Index

Full keyword list from gaussian.com/keywords/:

ADMP, BD, BOMD, CacheSize, CASSCF, CBS Methods, CBSExtrapolate, CCD and CCSD, Charge, ChkBasis, CID and CISD, CIS, CNDO, Complex, Constants, Counterpoise, CPHF, Density, DensityFit and NoDensityFit, DFT Methods, DFTB and DFTBA, EET, EOMCCSD, EPT, External, ExtraBasis & ExtraDensityBasis, Field, FMM, Force, Freq, Gen and GenECP, GenChk, Geom, GFInput, GFPrint, Gn Methods, Guess, GVB, HF, Huckel, INDO, Integral, IOp, IRC, IRCmax, Link0 Commands, LSDA, MaxDisk, MINDO3, MNDO, MM Methods, MP Methods, Name, NMR, ONIOM, Opt, Output, PBC, Polar, Population, Pressure, Prop, Pseudo, Punch, QCI, Restart, SAC-CI, Scale, Scan, SCF, SCRF, Semi-Empirical Methods, SP, Sparse, Stable, Symmetry, TD, Temperature, Test, TestMO, TrackIO, Transformation, Units, Volume, W1 Methods, Window Keyword and Frozen Core Options, ZIndo

---

## Advanced Topics

###ONIOM (Embedded Fragment)
```
# oniom(b3lyp/6-31G*:pm3)=(...)

Layer definition in molecule section:
H  -1 0.0 0.0 0.0 q=0   # High layer
L  0  0.0 0.0 0.0 q=+1  # Low layer
...
```

### EOM-CCSD (Excited States / Ionization)
```
# eomccsd                # EOM-CCSD for excited states
# eom(ip)                # Ionization potentials
# eom(ea)                # Electron affinities
# eom(ee)                # Two-electron excitations
```

### Density Fitting
```
# densityfit             # Enable DF for all calculations
# scf=(nodensityfit)     # Disable DF for SCF only
# mp2=fulldensityfit     # Full density fitting for MP2
```

### PBC (Periodic Boundary Conditions)
```
# pbc载体                 # Enable PBC calculations
# pbe0/pcs Solvent=...   # With SCRF
```
