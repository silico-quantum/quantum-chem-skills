---
name: pyscf
version: 2.0.0
description: PySCF (Python-based Simulations of Chemistry Framework) — modular pure-Python quantum chemistry library. Supports HF, DFT, MP2, CCSD, CCSD(T), CASSCF, TDDFT, geometry optimization, PES scanning, and spectroscopy. Based on pyscf.org documentation.
homepage: https://pyscf.org
metadata:
  category: computational_chemistry
  tags: [quantum_chemistry, DFT, TDDFT, MP2, CCSD, CASSCF, excited_states, python]
---

# PySCF — Python-based Simulations of Chemistry Framework

## Overview

PySCF is an open-source pure-Python quantum chemistry package featuring:
- **Wave function methods:** HF, MP2, MP3, CCSD, CCSD(T), CASSCF, CASCI, NEVPT2, DMRG
- **DFT:** All standard functionals (LDA, GGA, meta-GGA, hybrid, range-separated)
- **Excited states:** TDDFT, TDA, CIS, SF-DFT, EOM-CCSD
- **Geometry optimization:** Energy minimization, transition states, numerical Hessians
- **PES scanning:** Rigid scans, relaxed 1D/2D scans
- **Spectroscopy:** UV-Vis absorption, emission, CD, Raman
- **Periodic systems:** PBC-DFT, PBC-TDDFT

**Installation:** `pip install pyscf`

---

## Core Concepts

### Molecule Definition

```python
from pyscf import gto

mol = gto.M(
    atom='H 0 0 0; F 0 0 1.0',   # Inline coordinates
    # atom=open('h2o.xyz').read(),  # Or read from file
    basis='cc-pvdz',              # Basis set
    charge=0,                    # Molecular charge
    spin=0,                      # Spin multiplicity (0=closed-shell singlet)
    verbose=4                    # Output level
)
```

### Basis Set Selection

| Task | Recommended Basis |
|------|-----------------|
| Rapid screening | `sto-3g`, `3-21G` |
| Balanced accuracy | `6-31G*`, `cc-pVDZ` |
| High accuracy | `cc-pVTZ`, `cc-pVQZ`, `def2-TZVP` |
| Heavy atoms | `LANL2DZ`, `Stuttgart RSC` |
| Core excitations | `cc-pCVTZ` |

### Common Functionals

```python
# GGA functionals
'pbe', 'blyp', 'bp86', 'pbesol'

# Meta-GGA functionals
'scan', 'tpss', 'm06-l'

# Hybrid functionals
'b3lyp', 'pbe0', 'm06', 'm06-2x', 'wb97x', 'wb97x-d', 'wb97x-d3bj'

# Range-separated functionals (best for CT/excited states)
'wb97x', 'wb97x-d', 'wb97x-d3bj', 'camb3lyp', 'lc-wpbe'

# Dispersion correction
'b3lyp-d3', 'pbe-d3', 'wb97x-d3bj'
```

### Running Modes

```python
mf = scf.RHF(mol)
mf.kernel()              # Explicit run
mf.run()                 # Same as above
mf.converged            # Check convergence
mf.e_tot                # Total energy

# Progressive options
mf = scf.RHF(mol).newton()          # Newton-Raphson acceleration
mf = scf.RHF(mol).density_fit()     # Density fitting (faster)
mf = scf.RHF(mol).x2c()             # Scalar relativistic (2-component)
```

---

## 1. Wave Function Methods

### 1.1 Hartree-Fock (HF)

**Restricted closed-shell RHF:**
```python
from pyscf import scf, gto

mol = gto.M(atom='H 0 0 0; F 0 0 1.0', basis='cc-pvdz')
mf = scf.RHF(mol).run()
print('E(RHF) =', mf.e_tot)
```

**Open-shell UHF/ROHF:**
```python
mol = gto.M(atom='O 0 0 0; O 0 0 1.2', spin=2, basis='cc-pvdz')  # O2 triplet
mf = scf.UHF(mol).run()     # Unrestricted
mf = scf.ROHF(mol).run()    # Restricted open-shell
```

**Dirac-Hartree-Fock (relativistic):**
```python
from pyscf import scf
mol = gto.M(atom='Hg 0 0 0', basis='cc-pvdz')
mf = scf.DHF(mol).run()
```

### 1.2 MP2 (Second-Order Møller-Plesset)

```python
from pyscf import mp, scf, gto

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# Closed-shell MP2
mmp = mp.MP2(mf).run()
print('E(MP2) =', mmp.e_tot)
print('E(corr) =', mmp.e_corr)     # Correlation energy
print('T2 shape:', mmp.t2.shape)

# Open-shell UMP2
mmp = mp.UMP2(mf).run()
```

### 1.3 CCSD / CCSD(T)

**CCSD (Coupled-Cluster Singles and Doubles):**
```python
from pyscf import cc, scf, gto

mol = gto.M(atom='N 0 0 0; N 0 0 1.1', basis='cc-pvdz')
mf = scf.RHF(mol).run()

mycc = cc.CCSD(mf).run()
print('E(CCSD) =', mycc.e_tot)
print('T1 shape:', mycc.t1.shape)
print('T2 shape:', mycc.t2.shape)

# One-particle density matrix
dm = mycc.make_rdm1()
```

**CCSD(T) (perturbative triples):**
```python
et = cc.CCSD(mf).ccsd_t()    # Returns ΔE
e_total = mycc.e_tot + et
```

**Unrestricted UCCSD:**
```python
mol = gto.M(atom='N 0 0 0', spin=1, basis='cc-pvdz')
mf = scf.UHF(mol).run()
mycc = cc.UCCSD(mf).run()
```

### 1.4 CASSCF / CASCI

**CASSCF (Complete Active Space SCF):**
```python
from pyscf import mcscf, scf, gto

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# CAS(nelec, ncas) — 10 electrons, 8 orbital active space
mc = mcscf.CASSCF(mf, 8, 10).run()
print('E(CASSCF) =', mc.e_tot)
print('CAS orbitals:', mc.mo_coeff.shape)

# Custom active space (manual orbital selection)
mc = mcscf.CASSCF(mf, 8, 10)
mo = mcscf.sort_mo(mc, mf.mo_coeff, [15, 16, 17, 20, 21, 22, 25, 26])
mc.run(mo)
```

**CASCI:**
```python
mc = mcscf.CASCI(mf, 8, 10).run()
```

**CASPT2 (multi-state second-order perturbation):**
```python
mc = mcscf.CASSCF(mf, 8, 10).run()
mcpt2 = mcscf.CASPT2(mc).run()
print('E(CASPT2) =', mcpt2.e_tot)
```

**NEVPT2 (N-Electron Valence Second-Order PT):**
```python
mcpt2 = mcscf.NEVPT2(mc).run()
```

---

## 2. DFT Methods

### 2.1 Basic DFT
```python
from pyscf import dft

mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')

# RKS = Restricted Kohn-Sham (closed-shell)
mf = dft.RKS(mol, xc='pbe').run()
print('E(PBE) =', mf.e_tot)

mf = dft.RKS(mol, xc='b3lyp').run()
print('E(B3LYP) =', mf.e_tot)

# UKS = Unrestricted Kohn-Sham (open-shell)
mol = gto.M(atom='O 0 0 0', spin=2, basis='cc-pvdz')
mf = dft.UKS(mol, xc='pbe').run()
```

### 2.2 Functional Selection Guide

```python
xc_options = {
    # LDA
    'lda': 'svwn5',           # VWN5 Slater-VWN (RPA)

    # GGA
    'pbe': 'pbe',             # Perdew-Burke-Ernzerhof
    'blyp': 'blyp',           # Becke-Lee-Yang-Parr
    'bp86': 'bp86',           # B88 + P86

    # Meta-GGA
    'scan': 'scan',           # Strongly Constrained

    # Hybrid
    'b3lyp': 'b3lyp',         # Most popular
    'pbe0': 'pbe0',           # 25% HF exchange
    'm06': 'm06',             # Minnesota 2006
    'm06-2x': 'm06-2x',      # 54% HF exchange
    'wb97x': 'wb97x',        # Range-separated
    'wb97x-d3': 'wb97x-d3bj', # with dispersion

    # Dispersion correction
    'pbe-d3': 'pbe,pbe-d3bj',
    'b3lyp-d3': 'b3lyp,b3lyp-d3bj',
}
```

### 2.3 Range-Separated Functionals
```python
# ωB97X-D: Long-range corrected, ideal for CT and charge-transfer states
mf = dft.RKS(mol, xc='wb97x-d').run()

# CAM-B3LYP: Variable range separation
mf = dft.RKS(mol, xc='camb3lyp').run()
```

### 2.4 Orbital Visualization / DOS
```python
# MO coefficients and energies
mo_coeff = mf.mo_coeff
mo_energy = mf.mo_energy  # in Hartree

# Mulliken population
from pyscf import lo
occ = lo.orth_ao(mol, 'minao')
dm = mf.make_rdm1()
mulliken = lo.mulliken(mol, dm, occ)
```

---

## 3. Excited States: TDDFT / TDA

### 3.1 TDDFT (Linear Response)
```python
from pyscf import tddft, dft, gto

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# TDDFT calculation
td = tddft.TDDFT(mf).run()
print('Excitation energies (eV):', td.e * 27.2114)
print('Oscillator strengths:', td.f)
print('S1 energy (eV):', td.e[0] * 27.2114)

# Natural transition orbitals (NTO)
td.analyze()
```

### 3.2 TDA (Tamm-Dancoff Approximation)
```python
# TDA is faster than TDDFT but less accurate for double excitations
td = tddft.TDA(mf).run(nstates=5)   # default: 10 states
td.analyze()
```

### 3.3 UV-Vis Spectrum
```python
import numpy as np
import matplotlib.pyplot as plt

# Compute 50 excited states
td = tddft.TDDFT(mf).run(nstates=50)

# Broadened spectrum
e = td.e * 27.2114   # eV
f = td.f
sigma = 0.1           # eV

energies = np.linspace(0, 10, 500)
spectrum = np.zeros_like(energies)
for ei, fi in zip(e, f):
    spectrum += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))

plt.plot(energies, spectrum)
plt.xlabel('Energy (eV)')
plt.ylabel('Oscillator strength')
plt.show()
```

### 3.4 Open-Shell TDDFT
```python
mol = gto.M(atom='O 0 0 0', spin=2, basis='cc-pvdz')
mf = dft.UKS(mol, xc='pbe').run()
td = tddft.TDDFT(mf).run(nstates=5)
```

### 3.5 SF-DFT (Spin-Flip DFT, for diradicals)
```python
from pyscf import tddft, dft
mf = dft.RKS(mol, xc='b3lyp').run()
td = tddft.TDDFT(mf).run(nstates=3)

from pyscf import df
td_sf = df.DFTStateTransfer(mf, ['4.0', '5.0'])  # energy gap range
```

---

## 4. Geometry Optimization

### 4.1 Energy Minimization
```python
from pyscf import dft, geomopt

mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# Standard optimization
opt = geomopt.GeometryOptimizer(mf)
mol_opt = opt.run()

print('Optimized coordinates:')
for i in range(mol_opt.natm):
    print(mol_opt.atom_coord(i))

# Direct method
mol_opt = dft.RKS(mol, xc='b3lyp').geomopt().run()
```

### 4.2 Transition State Optimization
```python
from pyscf import dft, geomopt

mol = gto.M(atom='H 0 0 0; Cl 0 0 1.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# Dimer method for transition states
opt = geomopt.DimerOptimizer(mf).run()
```

### 4.3 Frequency Calculation (Verify Minima / TS)
```python
from pyscf import hessian

mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()
hess = hessian.RHF(mf).run()

# Frequencies
freq = hess.kernel()
print('Frequencies (cm^-1):', freq)
```

---

## 5. Potential Energy Surface (PES) Scanning

### 5.1 Relaxed Scan (Stepwise Optimization)
```python
from pyscf import dft
import numpy as np

mol = gto.M(atom='C 0 0 0; H 0 0 1.1', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# Scan C-H bond length
scan_coords = []
scan_energies = []

for r in np.linspace(0.9, 1.3, 10):
    mol_temp = gto.M(
        atom='C 0 0 0; H 0 0 %.4f' % r,
        basis='cc-pvdz'
    )
    mf_temp = dft.RKS(mol_temp, xc='pbe').run(chkfile=False)
    scan_coords.append(r)
    scan_energies.append(mf_temp.e_tot)

print('Scan complete')
for c, e in zip(scan_coords, scan_energies):
    print(f'r={c:.3f} E={e:.6f}')
```

### 5.2 2D Relaxed Scan
```python
# Combine two internal coordinates in nested loops
```

---

## 6. Spectroscopy

### 6.1 UV-Vis Absorption Spectrum
```python
import numpy as np
import matplotlib.pyplot as plt

def uv_vis_spectrum(td, nstates=50, sigma=0.05, emin=0, emax=10):
    """Generate broadened UV-Vis spectrum"""
    e = td.e[:nstates] * 27.2114   # eV
    f = td.f[:nstates]

    energies = np.linspace(emin, emax, 1000)
    absorption = np.zeros_like(energies)

    for ei, fi in zip(e, f):
        if fi > 0:
            absorption += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))

    return energies, absorption

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='b3lyp').run()
td = tddft.TDDFT(mf).run(nstates=50)

e, spec = uv_vis_spectrum(td, sigma=0.05)
plt.plot(e, spec)
plt.xlabel('Photon energy (eV)')
plt.ylabel('Absorption (arb. units)')
plt.show()
```

### 6.2 Emission Spectrum (from S1)
```python
from pyscf import tddft, dft

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()
td = tddft.TDDFT(mf).run(nstates=5)

# S1 emission = S1 → S0 vertical transition
td.emit()   # Compute emission (experimental feature)
```

### 6.3 CD Spectrum (Circular Dichroism)
```python
from pyscf import tddft, gdft

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = gdft.RKS(mol, xc='pbe').run()   # Use GDFT for CD
td = tddft.TDDFT(mf).run(nstates=20)

# Rotatory strength
td.analyze()
```

---

## 7. Parallelization & Performance

### 7.1 Density Fitting (DF)
```python
from pyscf import df, scf

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='def2-tzvp')
mf = scf.RHF(mol).density_fit(auxbasis='def2-tzvp-ri').run()
```

### 7.2 Multi-threaded BLAS
```python
import os
os.environ['OMP_NUM_THREADS'] = '16'

# Or at runtime
from pyscf.lib import num_threads
num_threads(16)
```

### 7.3 CHK Files (Checkpoint / Restart)
```python
mf = scf.RHF(mol)
mf.chkfile = 'hf.chk'
mf.run()

# Reload
mf = scf.RHF(mol)
mf.init_guess = 'chk'
mf.run(chkfile='hf.chk')
```

---

## 8. Common Basis Sets

| Basis Set | Description | Use Case |
|-----------|-------------|----------|
| `sto-3g` | Minimal basis, fastest | Quick testing |
| `3-21G` | Split-valence | Small molecule screening |
| `6-31G*` | Double-zeta + polarization | Standard organic molecules |
| `cc-pVDZ` | Correlation-consistent DZ | MP2/CCSD |
| `cc-pVTZ` | Correlation-consistent TZ | High precision |
| `def2-TZVP` | Pseudopotential basis | Heavy atoms |
| `LANL2DZ` | Effective core potential | Transition metals |

---

## 9. Troubleshooting

**Q: SCF not converging?**
```python
# Method 1: Change initial guess
mf.init_guess = 'atom'    # Atomic density
mf.init_guess = 'hcore'   # Core Hamiltonian
mf.diis = True            # DIIS acceleration (default on)

# Method 2: Newton-Raphson
mf = scf.RHF(mol).newton().run()

# Method 3: Change SCF method
mf = scf.DHF(mol).run()   # Dirac-HF (heavy atoms)
```

**Q: Out of memory?**
```python
# Use density fitting to reduce memory
mf = scf.RHF(mol).density_fit().run()

# Use direct SCF
mf.direct_scf = True
```

**Q: How to print all orbital energies?**
```python
mf = scf.RHF(mol).run()
print('HOMO:', mf.mo_energy[mf.mo_occ > 0].max())
print('LUMO:', mf.mo_energy[mf.mo_occ == 0].min())
```

---

## Tools Script Index

| Script | Description |
|--------|-------------|
| `tools/scf.py` | HF/DFT basic framework |
| `tools/tddft.py` | TDDFT/TDA calculation |
| `tools/dft.py` | DFT functional selection & comparison |
| `tools/mp2.py` | MP2 (RMP2/UMP2) |
| `tools/ccsd.py` | CCSD/CCSD(T) framework |
| `tools/cascf.py` | CASSCF/CASPT2/NEVPT2 framework |
| `tools/geometry.py` | Geometry optimization, TS, frequencies |
| `tools/spectrum.py` | UV-Vis/CD spectrum generation |
| `tools/pes.py` | 1D/2D potential energy surface scanning |
| `tools/analysis.py` | Wavefunction analysis (Mulliken, Löwdin, DOS, NBO) |
