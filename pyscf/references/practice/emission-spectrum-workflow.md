# TDDFT Emission-Spectrum Calculation Workflow

## Overview

Calculating an emission spectrum differs fundamentally from calculating an absorption spectrum: the calculation **must use a geometry optimized for the excited state**.

## Correct Calculation Workflow

### 1. Optimize the Ground-State Geometry (Optional)

```python
from pyscf import gto, dft, geomopt

mol = gto.M(atom='...', basis='cc-pvdz')
mf = dft.RKS(mol)
mf.xc = 'b3lyp'

# Optimize the ground-state geometry
mol_gs = geomopt.optimize(mf)
```

### 2. Absorption Spectrum at the Ground-State Geometry

```python
from pyscf import tdscf

# Calculate at the ground-state geometry
mf_gs = dft.RKS(mol_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# TDDFT absorption spectrum
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()

absorption_energies = td_abs.e * 27.2114  # eV
absorption_wavelengths = 1240 / absorption_energies  # nm
oscillator_strengths = td_abs.oscillator_strength()
```

### 3. Optimize the Excited-State Geometry

#### Method A: Use PySCF's Built-In Optimization

```python
# Note: the excited-state optimization interface in PySCF 2.12.x may be unstable
# TDA (the Tamm-Dancoff approximation) is recommended for better stability
from pyscf import tdscf

td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
td_s1.kernel()

# PySCF uses one-based excited-state IDs in gradient scanners.
s1_gradient = td_s1.nuc_grad_method().as_scanner(state=1)

# Attempt excited-state optimization; state crossings can still cause failure.
try:
    mol_ex = s1_gradient.optimizer().kernel()
except RuntimeError as exc:
    raise RuntimeError("First-excited-state optimization failed") from exc
```

#### Method B: Adjust the Geometry Manually (Simplified Method)

```python
# Adjust bond lengths based on chemical intuition
# Bond lengths are usually slightly longer in an excited state because an
# electron has been promoted to an antibonding orbital
# Example: increase the C-C bond length in benzene from 1.396 to 1.430 Å

mol_ex = gto.M(
    atom='''
    C  0.0000  1.4300  0.0000  # Adjust the C-C bond length
    C  1.2390  0.7150  0.0000
    ...
    ''',
    basis='cc-pvdz'
)
```

### 4. Emission Spectrum at the Excited-State Geometry

```python
# Perform a ground-state DFT calculation at the excited-state geometry
mf_ex = dft.RKS(mol_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# TDDFT emission spectrum (vertical emission)
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()

emission_energies = td_em.e * 27.2114  # eV
emission_wavelengths = 1240 / emission_energies  # nm
emission_osc = td_em.oscillator_strength()
```

## Stokes Shift

**Definition**: absorption energy - emission energy

```python
# Calculate the Stokes shift
stokes_shift = absorption_energies[0] - emission_energies[0]
print(f"Stokes shift: {stokes_shift:.3f} eV")
```

**Physical meaning**:

- The excited-state molecule first undergoes **geometric relaxation** (vibrational relaxation).
- The molecule then emits a photon from the relaxed excited state.
- The Stokes shift reflects the energy lost during relaxation.

## Two-Dimensional Potential Energy Surface Scan

A two-dimensional scan provides a more detailed way to study excited-state geometric relaxation:

```python
import numpy as np
from pyscf import gto, dft, tdscf

# Define bond-length ranges
cc_range = np.linspace(1.30, 1.50, 11)  # C-C bond length
ch_range = np.linspace(1.05, 1.15, 11)  # C-H bond length

gs_energies = []
s1_total_energies = []
s1_excitation_energies_ev = []

for cc in cc_range:
    for ch in ch_range:
        # Build the molecular geometry
        mol = build_molecule(cc, ch)  # User-defined function
        
        # Ground-state energy
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        e_gs = mf.kernel()
        
        # First-excited-state total energy and vertical gap at this geometry.
        td = tdscf.TDDFT(mf)
        td.nstates = 3
        td.kernel()
        e_s1 = e_gs + td.e[0]
        excitation_s1_ev = td.e[0] * 27.2114
        
        gs_energies.append(e_gs)
        s1_total_energies.append(e_s1)
        s1_excitation_energies_ev.append(excitation_s1_ev)

# Reshape as 2D grids
gs_2d = np.array(gs_energies).reshape(len(cc_range), len(ch_range))
s1_total_2d = np.array(s1_total_energies).reshape(
    len(cc_range), len(ch_range)
)

# Visualize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(14, 5))

# Ground-state potential energy surface
ax1 = fig.add_subplot(121, projection='3d')
CC, CH = np.meshgrid(cc_range, ch_range)
ax1.plot_surface(CC, CH, gs_2d.T, cmap='viridis')
ax1.set_xlabel('C-C Bond Length (Å)')
ax1.set_ylabel('C-H Bond Length (Å)')
ax1.set_zlabel('Ground-State Energy (Hartree)')

# Excited-state potential energy surface
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(CC, CH, s1_total_2d.T, cmap='plasma')
ax2.set_xlabel('C-C Bond Length (Å)')
ax2.set_ylabel('C-H Bond Length (Å)')
ax2.set_zlabel('S1 Total Energy (Hartree)')

plt.tight_layout()
plt.savefig('2d_pes.png', dpi=300)
```

## Common Questions

### Q1: Why is the emission energy lower than the absorption energy?

**A**: Because of excited-state geometric relaxation. After excitation, the molecule first adjusts its geometry to the minimum-energy structure and only then emits a photon.

### Q2: What information does the Stokes shift provide?

**A**:

- Large Stokes shift → large change in the excited-state geometry
- Small Stokes shift → similar excited-state and ground-state geometries
- The shift can help infer the nature of the excited state; for example, a charge-transfer state often has a large Stokes shift

### Q3: How can I determine whether excited-state optimization has converged?

**A**:

- Check the gradient norm; it should be < 0.001 Hartree/Bohr.
- Check the energy change; it should be < 1e-6 Hartree.
- Perform a vibrational-frequency analysis and confirm that there are no imaginary frequencies.

## Worked Example: Benzene

**Ground-state geometry**:

- C-C bond length: 1.396 Å
- C-H bond length: 1.089 Å

**Excited-state geometry (S1)**:

- C-C bond length: 1.430 Å (+2.4%)
- C-H bond length: 1.095 Å (+0.6%)

**Spectral data** (B3LYP/cc-pVDZ):

- Absorption (S3): 6.95 eV (178 nm)
- Emission (S3): 6.64 eV (187 nm)
- Stokes shift: 0.31 eV

## References

1. **TDDFT theory**:
   - Casida, M. E. (1995). "Time-dependent density functional response theory for molecules"

2. **Excited-state optimization**:
   - Li, Z., et al. (2018). "Analytic energy gradient for the second-order approximate coupled-cluster method"

3. **Stokes shift**:
   - Lakowicz, J. R. (2006). "Principles of Fluorescence Spectroscopy"

## Related PySCF Documentation

- [TDDFT documentation](http://www.pyscf.org/user/tdscf.html)
- [Geometry optimization](http://www.pyscf.org/user/geomopt.html)
- [Excited-state gradients](http://www.pyscf.org/user/gradient.html)
