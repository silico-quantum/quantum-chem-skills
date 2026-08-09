---
name: molecular-orbital-analysis
description: Use when analyzing and visualizing molecular orbitals with PySCF, Multiwfn, and PyMOL-compatible outputs.
license: MIT
compatibility: Requires Python, PySCF, Multiwfn, and optionally PyMOL; installation details are platform-specific.
---

# Molecular Orbital Analysis Workflow Guide

This skill describes a workflow for analyzing molecular orbitals with quantum
chemistry calculations and visualization. It is a guide: this directory does
not contain a complete automation program.

## Overview

The workflow consists of three main steps:
1. **PySCF** - Quantum chemistry calculation (HF/DFT)
2. **Multiwfn** - Generate 3D orbital data (cube files)
3. **PyMOL** - Professional visualization with molecular structure

## Requirements

### Software Installation

```bash
# Python packages used by the examples
python -m pip install pyscf numpy matplotlib

# Verify separately installed external programs.
command -v Multiwfn
command -v pymol
```

Install Multiwfn and PyMOL from their upstream distribution channels or a
trusted platform package. Confirm the executable names and versions locally.

### Environment Setup

```bash
# Add to ~/.zshrc or ~/.bashrc
export OMP_STACKSIZE=64000000
```

## Complete Workflow

### Step 1: Prepare Molecular Structure

Create an XYZ file with molecular coordinates:

```xyz
<n_atoms>
<molecule_name>
<element> <x> <y> <z>
...
```

**Example - Water:**
```xyz
3
Water molecule
O    0.000000    0.000000    0.117489
H    0.000000    0.757210   -0.469957
H    0.000000   -0.757210   -0.469957
```

### Step 2: Quantum Chemistry Calculation (PySCF)

Use the template script:

```python
#!/usr/bin/env python3
from pyscf import dft, gto, scf
from pyscf.tools import molden

# Create molecule
mol = gto.M(
    atom='<xyz_coordinates_or_file>',
    basis='6-31G*',  # or 'cc-pVDZ', 'def2-TZVP', etc.
    charge=0,
    spin=0,  # 0 for closed-shell
    verbose=3
)

# Perform calculation
# Hartree-Fock example:
mf = scf.RHF(mol)

# For a DFT calculation instead, use for example:
# mf = dft.RKS(mol)
# mf.xc = 'PBE0'
mf.kernel()

# Output results
print(f"Total energy: {mf.e_tot:.8f} Hartree")
nocc = mol.nelectron // 2
print(f"HOMO energy: {mf.mo_energy[nocc-1]*27.2114:.3f} eV")
print(f"LUMO energy: {mf.mo_energy[nocc]*27.2114:.3f} eV")

# Generate molden file
molden.from_mo(mol, 'molecule.molden', mf.mo_coeff)
```

### Step 3: Generate Orbital Data (Multiwfn)

Create input file for Multiwfn:

```bash
cat > multiwfn_input.txt << 'EOF'
molecule.molden
200
3
<orbital_index>
2
1
0
q
EOF

OMP_STACKSIZE=64000000 Multiwfn < multiwfn_input.txt
```

**Finding orbital indices:**
- HOMO index = `nelectron // 2`
- LUMO index = HOMO + 1

**Output:** `orb<index>.cub` file containing 3D orbital data

### Step 4: Visualization (PyMOL)

Create PyMOL script:

```python
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])

# Load molecule
cmd.load('molecule.xyz', 'mol')
cmd.show('spheres', 'mol')
cmd.show('sticks', 'mol')
cmd.color('red', 'elem O')
cmd.color('white', 'elem H')
cmd.set('sphere_scale', 0.3)

# Load orbital
cmd.load('orbital.cub', 'orbital')
cmd.isosurface('surf_pos', 'orbital', level=0.05)
cmd.isosurface('surf_neg', 'orbital', level=-0.05)
cmd.color('blue', 'surf_pos')
cmd.color('red', 'surf_neg')
cmd.set('transparency', 0.4, 'surf_pos')
cmd.set('transparency', 0.4, 'surf_neg')

# Render
cmd.bg_color('white')
cmd.zoom(complete=1)
cmd.png('orbital.png', width=1400, height=1050, dpi=150, ray=1)

cmd.quit()
```

## Automation Skeleton

The following skeleton illustrates orchestration only. The Multiwfn and PyMOL
sections are placeholders and must be implemented and validated before use.

```python
#!/usr/bin/env python3
"""
Molecular orbital analysis pipeline skeleton
Usage: python analyze_molecule.py molecule.xyz
"""
import sys
import subprocess
from pathlib import Path
from pyscf import gto, scf
from pyscf.tools import molden
import pymol
from pymol import cmd

def analyze_molecule(xyz_file, basis='6-31G*', homo_only=False):
    """Run the implemented portions of the orbital-analysis skeleton."""

    # Step 1: Quantum chemistry calculation
    print("=" * 60)
    print("Step 1: Quantum Chemistry Calculation")
    print("=" * 60)

    mol = gto.M(atom=xyz_file, basis=basis, charge=0, spin=0, verbose=3)
    mf = scf.RHF(mol)
    mf.kernel()

    # Generate molden file
    molden_file = Path(xyz_file).stem + '.molden'
    molden.from_mo(mol, molden_file, mf.mo_coeff)

    # Step 2: Generate cube files with Multiwfn
    print("\n" + "=" * 60)
    print("Step 2: Generating Orbital Data")
    print("=" * 60)

    nocc = mol.nelectron // 2
    orbitals = [nocc, nocc+1] if not homo_only else [nocc]

    for orb_idx in orbitals:
        print(f"Generating cube file for orbital {orb_idx}...")
        # Multiwfn input
        # ... (automation code)

    # Step 3: Visualization with PyMOL
    print("\n" + "=" * 60)
    print("Step 3: Visualization")
    print("=" * 60)

    # PyMOL visualization code
    # ... (automation code)

    print("\nImplemented steps complete; validate generated files before use.")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python analyze_molecule.py molecule.xyz")
        sys.exit(1)

    analyze_molecule(sys.argv[1])
```

## Common Parameters

### Basis Sets
- **6-31G***: Pople split-valence basis with polarization
- **cc-pVDZ**: Correlation-consistent double-zeta basis
- **def2-TZVP**: Karlsruhe triple-zeta valence basis with polarization
- **STO-3G**: Minimal basis, generally suitable only for demonstrations

Choose and justify a basis for the target property and chemical system; this
list is not an accuracy ranking.

### Isosurface Levels
- **0.05**: example starting value
- **0.03**: lower-value surface that may show more diffuse regions
- **0.08**: higher-value surface that emphasizes larger amplitudes

Report the chosen value and use the same absolute positive and negative levels
when comparing orbitals.

### PyMOL Customization

**Colors:**
- Carbon: gray, cyan, green
- Oxygen: red
- Nitrogen: blue
- Hydrogen: white

**Zoom level:**
- `cmd.zoom(complete=1)`: Fit entire molecule
- `cmd.move('z', -5)`: Zoom in
- `cmd.move('z', 5)`: Zoom out

## Output Files

Expected output structure after the placeholder steps have been implemented:
```
<molecule>/
├── molecule.xyz           # Input structure
├── molecule.molden        # Wavefunction file
├── molecule_HOMO.cub     # HOMO 3D data
├── molecule_LUMO.cub     # LUMO 3D data
├── molecule_homo.png     # HOMO visualization
└── molecule_lumo.png     # LUMO visualization
```

## Examples

### Water (H₂O)
- 10 electrons, 5 occupied orbitals
- HOMO: Non-bonding orbital (lone pair on O)
- LUMO: Antibonding orbital

### Benzene (C₆H₆)
- 42 electrons, 21 occupied orbitals
- HOMO: π bonding orbital
- LUMO: π* antibonding orbital
- Characteristic delocalized π system

## Troubleshooting

### Multiwfn Issues
- **"settings.ini not found"**: Warning only, uses defaults
- **Slow calculation**: Use lower quality grid (option 1 or 2)

### PyMOL Issues
- **No output image**: Check if running in GUI mode
- **Black images**: OpenGL issue, try `ray=0` instead of `ray=1`

### PySCF Issues
- **SCF not converging**: Try different initial guess or use DIIS
- **Memory error**: Reduce basis set size or use density fitting

## References

- PySCF: https://pyscf.org/
- Multiwfn: http://sobereva.com/multiwfn/
- PyMOL: https://pymol.org/2/

## Citation

If using this workflow for publications, cite:

1. **PySCF**: Q. Sun et al., PySCF: the Python‐based simulations of chemistry framework, Wiley Interdiscip. Rev. Comput. Mol. Sci. 8, e1340 (2018)

2. **Multiwfn**: T. Lu, F. Chen, J. Comput. Chem. 33, 580 (2012) and T. Lu, J. Chem. Phys. 161, 082503 (2024)

3. **PyMOL**: The PyMOL Molecular Graphics System, Version 3.1 Schrödinger, LLC.
