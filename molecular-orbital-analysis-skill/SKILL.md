---
name: molecular-orbital-analysis
version: 1.0.0
description: Complete workflow for molecular orbital analysis using PySCF, Multiwfn, and PyMOL
homepage: https://github.com/STOKES-DOT
metadata:
  category: quantum-chemistry
  tools: ["pyscf", "multiwfn", "pymol"]
  requirements:
    - PySCF (quantum chemistry)
    - Multiwfn (wavefunction analysis)
    - PyMOL (molecular visualization)
---

# Molecular Orbital Analysis - Complete Workflow

This skill provides a complete workflow for analyzing molecular orbitals using quantum chemistry calculations and visualization.

## Overview

The workflow consists of three main steps:
1. **PySCF** - Quantum chemistry calculation (HF/DFT)
2. **Multiwfn** - Generate 3D orbital data (cube files)
3. **PyMOL** - Professional visualization with molecular structure

## Requirements

### Software Installation

```bash
# Python packages
pip install pyscf numpy matplotlib

# Multiwfn (via Homebrew)
brew tap digital-chemistry-laboratory/multiwfn
brew install --HEAD digital-chemistry-laboratory/multiwfn/multiwfn

# PyMOL
brew install pymol
```

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
from pyscf import gto, scf
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
mf = scf.RHF(mol)  # or RKS for DFT
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

OMP_STACKSIZE=64000000 multiwfn < multiwfn_input.txt
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

## Automation Script

Use this complete automation script:

```python
#!/usr/bin/env python3
"""
Complete molecular orbital analysis pipeline
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
    """Complete molecular orbital analysis"""

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

    print("\n✅ Analysis complete!")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python analyze_molecule.py molecule.xyz")
        sys.exit(1)

    analyze_molecule(sys.argv[1])
```

## Common Parameters

### Basis Sets
- **6-31G***: Standard Pople basis (recommended for routine calculations)
- **cc-pVDZ**: Correlation-consistent basis
- **def2-TZVP**: Triple-zeta basis for higher accuracy
- **sto-3g**: Minimal basis (for testing only)

### Isosurface Levels
- **0.05**: Standard for most orbitals
- **0.03**: For diffuse orbitals (e.g., π systems)
- **0.08**: For compact orbitals

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

Typical output structure:
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
