---
name: rdkit-chemistry
version: 1.0.0
description: RDKit molecular analysis and visualization. Use for molecular conformer generation, force field optimization, charge calculation, molecular descriptors, and non-covalent interaction analysis.
homepage: https://www.rdkit.org
---

# RDKit Chemistry Analysis Skill

Molecular structure analysis and visualization using RDKit.

## When to Use

Activate this skill when:
- Building 3D molecular conformers from SMILES
- Performing force field optimization (MMFF94, UFF)
- Calculating molecular descriptors (LogP, TPSA, molecular weight)
- Computing atomic charges (Gasteiger, Mulliken)
- Analyzing non-covalent interactions (π-π stacking, hydrogen bonding)
- Generating molecular visualizations
- Analyzing D-A (donor-acceptor) systems for TADF materials

## Quick Start

### 1. Basic Molecular Setup

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# From SMILES to 3D
smiles = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)  # Add hydrogens

# Generate 3D conformer
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# MMFF94 optimization
AllChem.MMFFOptimizeMolecule(mol)
```

## Common Patterns

### Pattern 1: Quick Molecular Analysis

```python
# Full analysis pipeline
def analyze_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # 3D structure
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    # Descriptors
    from rdkit.Chem import Descriptors
    results = {
        'MW': Descriptors.ExactMolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
    }

    # Charges
    AllChem.ComputeGasteigerCharges(mol)
    charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms()]
    results['charges'] = charges

    return results
```

## Example Scripts

See examples/ directory for complete usage examples.

---
*Created: 2026-03-28*
*Version: 1.0.0*
