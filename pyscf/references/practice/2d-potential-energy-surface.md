# Two-Dimensional Potential Energy Surface Scan

## Overview

A two-dimensional potential-energy-surface (2D PES) scan samples the
relationship between selected molecular coordinates and electronic energy. A
constrained scan can help with:

- Comparing ground- and excited-state energy trends over selected coordinates
- Analyzing geometric relaxation effects
- Investigating the physical origin of a Stokes shift
- Exploring chemical reaction pathways

## Theoretical Background

### Potential Energy Surface (PES)

A potential energy surface describes the molecular energy, E, as a function of the nuclear coordinates, {R}:

```text
E = E(R₁, R₂, ..., R₃N)
```

For a molecule containing N atoms, the PES is a 3N-6-dimensional hypersurface for a nonlinear molecule or a 3N-5-dimensional hypersurface for a linear molecule.

### Two-Dimensional Scan

Fix two key coordinates, such as bond lengths r₁ and r₂, and calculate the energy:

```text
E(r₁, r₂) = E(r₁, r₂, all other coordinates fixed)
```

## Grid-Design Strategy

### 1. Selecting the Range

**Principle**: Start from physically plausible structures and choose a range
that brackets the region of interest. A pilot scan can identify whether the
range needs to be extended.

**Example (benzene)**:

```python
# Illustrative ranges only; justify them for the selected method and state.
cc_range = np.linspace(1.25, 1.80, 15)  # C-C bond length
ch_range = np.linspace(0.85, 1.20, 15)  # C-H bond length
```

### 2. Grid Density

| Grid | Number of points | Typical role |
|---|---:|---|
| 5×5 | 25 | Pilot scan |
| 7×7 | 49 | Refine a broad region |
| 10×10 | 100 | Inspect a narrower region |
| 15×15 | 225 | Higher-resolution sampling after convergence checks |

Runtime depends on the molecule, method, basis, state count, convergence, and
hardware. Grid density alone does not make a result publication-ready.

### 3. Basis-Set Selection

| Basis set | Possible role in this example |
|---|---|
| STO-3G | Interface or debugging test |
| 3-21G | Low-cost exploratory scan |
| cc-pVDZ | Larger comparison calculation |

Basis-set suitability is property- and system-dependent. Converge the result
with respect to basis and method rather than treating this ordering as a
universal accuracy ranking.

## Illustrative Code Template

### Basic Version (7×7 Grid)

```python
import numpy as np
from pyscf import gto, dft, tdscf

# Define the grid
cc_range = np.linspace(1.25, 1.80, 7)
ch_range = np.linspace(0.85, 1.20, 7)

gs_energies = []
s1_total_energies = []
excitation_energies_ev = []

for cc in cc_range:
    for ch in ch_range:
        # Build the molecule
        atoms = []
        for i in range(6):
            angle = i * np.pi / 3
            atoms.append(f"C {cc*np.cos(angle):.6f} {cc*np.sin(angle):.6f} 0.0")
            atoms.append(f"H {(cc+ch)*np.cos(angle):.6f} {(cc+ch)*np.sin(angle):.6f} 0.0")
        
        mol = gto.M(atom='\n'.join(atoms), basis='3-21g', verbose=0)
        
        # Ground-state calculation
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        e_gs = mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"SCF failed at C-C={cc:.4f}, C-H={ch:.4f} Å")
        
        # Excited-state calculation
        td = tdscf.TDDFT(mf)
        td.nstates = 3
        td.kernel()
        # td.e[0] is an excitation energy in Hartree, not a total energy.
        e_s1 = e_gs + td.e[0]
        
        gs_energies.append(e_gs)
        s1_total_energies.append(e_s1)
        excitation_energies_ev.append(td.e[0] * 27.211386245988)

# Reshape as 2D arrays
shape = (len(cc_range), len(ch_range))
gs_2d = np.array(gs_energies).reshape(shape)
s1_2d = np.array(s1_total_energies).reshape(shape)
excitation_2d = np.array(excitation_energies_ev).reshape(shape)
```

This template follows the lowest returned TDDFT root at each point. Near state
crossings, root order can change; production work needs state-character checks
or an overlap-based state-tracking procedure. A constrained two-coordinate
scan also does not optimize the remaining internal coordinates.

## Visualization

### Contour Plot (Recommended)

```python
import matplotlib.pyplot as plt

CC, CH = np.meshgrid(cc_range, ch_range)
gs_relative_kcal = (gs_2d - gs_2d.min()) * 627.509

fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contourf(CC, CH, gs_relative_kcal.T, levels=25, cmap='viridis')
plt.colorbar(contour, ax=ax, label='Relative energy (kcal mol⁻¹)')
ax.contour(CC, CH, gs_relative_kcal.T, levels=12, colors='white', alpha=0.4, linewidths=0.5)

# Mark the minimum
min_idx = np.unravel_index(gs_2d.argmin(), gs_2d.shape)
ax.plot(cc_range[min_idx[0]], ch_range[min_idx[1]], 'r*', markersize=20)

ax.set_xlabel('C-C Bond Length (Å)', fontsize=14)
ax.set_ylabel('C-H Bond Length (Å)', fontsize=14)
ax.set_title('Ground State Potential Energy Surface', fontsize=16, fontweight='bold')
ax.set_aspect('equal', adjustable='box')
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('2d_pes.png', dpi=300)
```

## Analysis

### Locate the Minima

```python
# Ground-state minimum
min_gs = np.unravel_index(gs_2d.argmin(), gs_2d.shape)
gs_min_cc = cc_range[min_gs[0]]
gs_min_ch = ch_range[min_gs[1]]

# Excited-state minimum
min_s1 = np.unravel_index(s1_2d.argmin(), s1_2d.shape)
s1_min_cc = cc_range[min_s1[0]]
s1_min_ch = ch_range[min_s1[1]]

print(f"Ground state S₀: C-C = {gs_min_cc:.4f} Å, C-H = {gs_min_ch:.4f} Å")
print(f"Excited state S₁: C-C = {s1_min_cc:.4f} Å, C-H = {s1_min_ch:.4f} Å")

def lies_on_boundary(index, shape):
    return any(i == 0 or i == n - 1 for i, n in zip(index, shape))

if lies_on_boundary(min_gs, gs_2d.shape) or lies_on_boundary(min_s1, s1_2d.shape):
    print("Warning: a grid minimum lies on the boundary; extend the scan before interpreting it.")
```

### Stokes Shift

```python
# Vertical absorption energy at the ground-state grid minimum
vertical_abs = excitation_2d[min_gs]

# Vertical S1-S0 gap at the excited-state grid minimum
vertical_em = excitation_2d[min_s1]

# Vertical Stokes-shift estimate
stokes = vertical_abs - vertical_em

print(f"Vertical absorption: {vertical_abs:.2f} eV ({1240/vertical_abs:.1f} nm)")
print(f"Vertical emission: {vertical_em:.2f} eV ({1240/vertical_em:.1f} nm)")
print(f"Vertical Stokes-shift estimate: {stokes:.3f} eV ({stokes*8065.5:.0f} cm⁻¹)")
```

## Interpreting a Benzene Scan

The previous version of this guide reported a boundary point as an S₁ minimum
after minimizing `td.e[0]` alone. That quantity is the vertical excitation
energy, not the excited-state total energy, so those numerical conclusions
have been removed.

Before reporting a minimum or a Stokes shift from a new scan:

1. minimize `E_S0(R)` and `E_S0(R) + omega_S1(R)` on their respective grids;
2. reject a claimed minimum that lies on a grid boundary and extend the range;
3. check SCF and response convergence at every point;
4. track the character of the target state rather than assuming root 1 remains
   the same state across the grid;
5. refine the grid and relax unscanned coordinates where the scientific
   question requires an equilibrium geometry;
6. distinguish a vertical energy-gap estimate from a vibronic or 0-0 Stokes
   shift.

---

**Created**: 2026-03-09

**Original draft environment**: PySCF v2.12.1. Revalidate API behavior and
numerical settings for the installed version.
