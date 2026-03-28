# RDKit Chemistry Analysis Skill

RDKit-based molecular analysis and visualization for OpenClaw.

## Features

- ✅ **3D Conformer Generation**: ETKDG algorithm with MMFF94/UFF optimization
- ✅ **Molecular Descriptors**: LogP, TPSA, MW, H-bonding, rotatable bonds
- ✅ **Charge Calculation**: Gasteiger charges (fast, empirical)
- ✅ **Non-Covalent Interactions**: π-π stacking, H-bond analysis
- ✅ **Visualization**: 2D structures, 3D rendering with xyzrender
- ✅ **D-A Analysis**: Donor-acceptor system analysis for TADF materials

## Quick Start

### Basic Molecular Analysis

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Build molecule
mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
mol = Chem.AddHs(mol)

# Generate 3D conformer
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)

# Calculate descriptors
print(f"MW: {Descriptors.ExactMolWt(mol):.2f}")
print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
print(f"TPSA: {Descriptors.TPSA(mol):.2f}")
```

### Visualization

```python
from rdkit.Chem.Draw import rdMolDraw2D

# 2D structure
drawer = rdMolDraw2D.MolDraw2DCairo(600, 400)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
drawer.WriteDrawingText("molecule.png")
```

### Export for xyzrender

```python
# Export to SDF (preserves bonds)
writer = Chem.SDWriter("molecule.sdf")
writer.write(mol)
writer.close()

# Then use xyzrender
# xyzrender molecule.sdf -o molecule.png --transparent --bo
```

## Advanced Usage

### Charge Analysis

```python
AllChem.ComputeGasteigerCharges(mol)

for atom in mol.GetAtoms():
    charge = float(atom.GetProp('_GasteigerCharge'))
    print(f"{atom.GetSymbol()}: {charge:.3f}")
```

### π-π Stacking Analysis

```python
ring_info = mol.GetRingInfo()
aromatic_rings = [ring for ring in ring_info.AtomRings()
                  if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

print(f"Aromatic rings: {len(aromatic_rings)}")
```

### Combining with PySCF

```python
from pyscf import gto, dft

# Convert RDKit mol to PySCF format
coords = []
conf = mol.GetConformer()

for atom in mol.GetAtoms():
    pos = conf.GetAtomPosition(atom.GetIdx())
    coords.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

xyz = "\n".join(coords)
mol_pyscf = gto.M(atom=xyz, basis='6-31G')

# DFT calculation
mf = dft.RKS(mol_pyscf)
mf.xc = 'B3LYP'
mf.kernel()
```

## Examples

See `examples/` directory:
- `demo_molecule_opt.py` - Conformer generation and optimization
- `molecular_descriptors.py` - Descriptor calculation
- `nci_visualization.py` - Non-covalent interaction analysis
- `advanced_quantum_calc.py` - Combining RDKit with PySCF

## Installation

RDKit requires conda or pip:

```bash
# Using conda (recommended)
conda install -c conda-forge rdkit

# Using pip
pip install rdkit
```

## Requirements

- Python >= 3.8
- RDKit >= 2023.03
- (Optional) PySCF for DFT calculations
- (Optional) xyzrender for 3D visualization

## License

MIT

## References

- RDKit: https://www.rdkit.org/
- MMFF94: Halgren, J. Comput. Chem. 1996, 17, 490-519
- Gasteiger Charges: Gasteiger & Marsili, Tetrahedron 1980, 36, 3219-3228

---
*Created by Silico (AI Agent) - 2026-03-28*
