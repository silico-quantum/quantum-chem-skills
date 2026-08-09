# PySCF Emission-Spectrum Calculation Guide

## Overview

An emission spectrum (fluorescence spectrum) describes light emitted when a molecule returns from an excited state to the ground state. Unlike an absorption spectrum, an emission spectrum must be calculated at the **geometry optimized for the excited state**.

## Physical Process

```text
Ground state (S₀) --[absorption]--> Excited state (S_n)
       ↓                                  ↓
Ground-state geometry              Excited-state geometry (relaxed)
                                          ↓
                                      [emission]
                                          ↓
                                   Ground state (S₀)
```

## Key Concepts

### 1. Stokes Shift

- **Definition**: absorption energy - emission energy
- **Cause**: excited-state geometry relaxation
- **Formula**: ΔE_Stokes = E_abs - E_em

### 2. Vertical Transitions

- **Vertical absorption**: calculate the excitation energy at the ground-state geometry
- **Vertical emission**: calculate the emission energy at the excited-state geometry
- Franck–Condon principle: an electronic transition is much faster than nuclear motion

### 3. Geometric Relaxation

An excited-state molecule undergoes:

1. Electronic excitation (~femtoseconds)
2. Vibrational relaxation (~picoseconds)
3. Geometric rearrangement (~nanoseconds)
4. Fluorescence emission (~nanoseconds)

## PySCF Calculation Workflow

### Standard Workflow

```python
from pyscf import gto, dft, tdscf

# 1. Define the molecule at its ground-state geometry
mol_gs = gto.M(atom='...', basis='cc-pvdz')

# 2. Ground-state DFT calculation
mf_gs = dft.RKS(mol_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# 3. Calculate the absorption spectrum at the ground-state geometry
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()
absorption_energies = td_abs.e * 27.2114  # eV

# 4. Optimize the first excited-state geometry
# Method 1: use TDA excited-state optimization
from pyscf import tdscf
td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
td_s1.kernel()
s1_gradient = td_s1.nuc_grad_method().as_scanner(state=1)
mol_ex = s1_gradient.optimizer().kernel()

# Method 2: adjust the geometry manually (approximation)
# Bond lengths are usually slightly longer in an excited state

# 5. Perform a ground-state calculation at the excited-state geometry
mf_ex = dft.RKS(mol_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# 6. Calculate the emission spectrum at the excited-state geometry
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()
emission_energies = td_em.e * 27.2114  # eV

# 7. Calculate the Stokes shift
stokes_shift = absorption_energies[0] - emission_energies[0]
```

## Complete Example: Benzene

```python
#!/usr/bin/env python3
"""
Example emission-spectrum calculation for benzene.
"""

import numpy as np
from pyscf import gto, dft, tdscf

# 1. Ground-state geometry
benzene_gs = gto.M(
    atom='''
    C  0.000  1.396  0.000
    C  1.210  0.698  0.000
    C  1.210 -0.698  0.000
    C  0.000 -1.396  0.000
    C -1.210 -0.698  0.000
    C -1.210  0.698  0.000
    H  0.000  2.476  0.000
    H  2.147  1.239  0.000
    H  2.147 -1.239  0.000
    H  0.000 -2.476  0.000
    H -2.147 -1.239  0.000
    H -2.147  1.239  0.000
    ''',
    basis='cc-pvdz'
)

# 2. Ground-state DFT
mf_gs = dft.RKS(benzene_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# 3. Absorption spectrum
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()
abs_e = td_abs.e * 27.2114

print(f"Absorption S₁: {abs_e[0]:.2f} eV ({1240/abs_e[0]:.0f} nm)")

# 4. Approximate excited-state geometry: increase bond lengths by 2-3%
benzene_ex = gto.M(
    atom='''
    C  0.000  1.430  0.000
    C  1.239  0.715  0.000
    C  1.239 -0.715  0.000
    C  0.000 -1.430  0.000
    C -1.239 -0.715  0.000
    C -1.239  0.715  0.000
    H  0.000  2.510  0.000
    H  2.190  1.255  0.000
    H  2.190 -1.255  0.000
    H  0.000 -2.510  0.000
    H -2.190 -1.255  0.000
    H -2.190  1.255  0.000
    ''',
    basis='cc-pvdz'
)

# 5. Ground-state calculation at the excited-state geometry
mf_ex = dft.RKS(benzene_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# 6. Emission spectrum
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()
em_e = td_em.e * 27.2114

print(f"Emission S₁→S₀: {em_e[0]:.2f} eV ({1240/em_e[0]:.0f} nm)")
print(f"Stokes shift: {abs_e[0] - em_e[0]:.3f} eV")
```

## Excited-State Optimization Methods

### Method 1: Built-In PySCF Optimization (Recommended)

```python
# TDA excited-state optimization (more stable)
td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
td_s1.kernel()
s1_gradient = td_s1.nuc_grad_method().as_scanner(state=1)
mol_ex = s1_gradient.optimizer().kernel()
```

### Method 2: Use an External Optimizer

```python
# Use the PyBerny optimizer
from pyscf.geomopt import berny_solver

s1_gradient = td_s1.nuc_grad_method().as_scanner(state=1)
mol_ex = berny_solver.optimize(s1_gradient)
```

### Method 3: Adjust the Geometry Manually (Quick Approximation)

Empirical rules for an excited-state geometry:

- **C-C bond**: increase by 2-4%
- **C-H bond**: increase by 1-2%
- **C=O bond**: increase by 3-5%

## Common Questions

### Q1: What if the excited-state optimization does not converge?

**Possible remedies:**

1. Use TDA instead of full TDDFT; TDA is more stable.
2. Improve the initial guess by starting from the ground-state geometry.
3. Loosen the convergence threshold (`conv_tol=1e-4`).

```python
td_s1 = tdscf.TDA(mf_gs)
td_s1.conv_tol = 1e-4
```

### Q2: Why is an oscillator strength zero (a forbidden transition)?

**Cause:** symmetry-forbidden transition

**Checks:**

```python
# Inspect the molecular symmetry
print(f"Molecular symmetry: {mol.symmetry}")

# Inspect excited-state symmetries
for i, e in enumerate(td.e):
    print(f"S{i+1}: {e*27.2114:.2f} eV, f={td.oscillator_strength()[i]:.4f}")
```

### Q3: What if the Stokes shift is too large or too small?

**Typical range:** 0.1-1.0 eV

**Potential anomalies:**

- >1.0 eV: check whether the geometry optimization is correct
- <0.05 eV: the excited-state geometry may have changed only slightly

## Performance Optimization

### 1. Calculate Fewer States

```python
# Calculate only the first three states
td.nstates = 3  # Instead of 10
```

### 2. Use Density Fitting

```python
mf = dft.RKS(mol).density_fit(auxbasis='def2-universal-jkfit')
```

### 3. Use a Coarse Grid for Testing

```python
mf.grids.atom_grid = (50, 194)  # Default: (75, 302)
```

## Validate the Results

### Checklist

- [ ] Ground-state SCF calculation converged
- [ ] TDDFT calculation converged
- [ ] Excited-state optimization converged
- [ ] Stokes shift is reasonable (0.1-1.0 eV)
- [ ] Oscillator strengths are nonnegative
- [ ] Wavelengths are in a reasonable ultraviolet/visible range

### Compare with Experimental Data

Typical Stokes shifts for aromatic molecules:

- Benzene: ~0.1-0.3 eV
- Naphthalene: ~0.2-0.4 eV
- Anthracene: ~0.3-0.5 eV

## Visualization

```python
import matplotlib.pyplot as plt

# Gaussian broadening
def gaussian(x, mu, sigma=0.15):
    return np.exp(-((x - mu)**2) / (2 * sigma**2))

energy_range = np.linspace(3, 8, 1000)

# Absorption spectrum
abs_spectrum = sum(gaussian(energy_range, e) * f 
                   for e, f in zip(abs_e, td_abs.oscillator_strength()))

# Emission spectrum
em_spectrum = sum(gaussian(energy_range, e) * f 
                  for e, f in zip(em_e, td_em.oscillator_strength()))

plt.plot(energy_range, abs_spectrum, 'b-', label='Absorption')
plt.plot(energy_range, em_spectrum, 'r-', label='Emission')
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity')
plt.legend()
plt.show()
```

## References

1. **TDDFT theory**
   - Casida, M. E. (1995). "Time-dependent density functional response theory for molecules"

2. **Excited-state optimization**
   - Furche, F. & Ahlrichs, R. (2002). "Adiabatic time-dependent density functional methods"

3. **Stokes shift**
   - Lakowicz, J. R. (2006). "Principles of Fluorescence Spectroscopy"

## Related Resources

- PySCF documentation: http://www.pyscf.org
- TDDFT examples: `pyscf/examples/tdscf/`
- Geometry-optimization examples: `pyscf/examples/geomopt/`
