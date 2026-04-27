#!/usr/bin/env python3
"""
Benzene (C₆H₆) - Comprehensive RDKit Visualization
A showcase of molecular analysis and visualization capabilities
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np

print("=" * 60)
print("Benzene (C₆H₆) - RDKit Visualization Showcase")
print("=" * 60)

# Step 1: Build molecule
smiles = "c1ccccc1"
mol = Chem.MolFromSmiles(smiles)
mol_3d = Chem.AddHs(mol)

# Generate 3D conformer
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)

print(f"\n✓ Molecular formula: {rdMolDescriptors.CalcMolFormula(mol)}")
print(f"✓ MW: {Descriptors.ExactMolWt(mol):.2f} Da")
print(f"✓ Atoms: {mol_3d.GetNumAtoms()}")

# Step 2: 2D Structure with atom indices
print("\n" + "=" * 60)
print("Generating visualizations...")
print("=" * 60)

# 2D structure with indices
drawer = rdMolDraw2D.MolDraw2DCairo(800, 800)
drawer.SetFontSize(0.9)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
drawer.WriteDrawingText("benzene_showcase_2d.png")
print("✓ 2D structure: benzene_showcase_2d.png")

# 2D structure with atom indices (highlighted)
drawer2 = rdMolDraw2D.MolDraw2DCairo(800, 800)
drawer2.SetFontSize(0.9)
drawer2.DrawMolecule(mol_3d)
drawer2.FinishDrawing()
drawer2.WriteDrawingText("benzene_showcase_2d_h.png")
print("✓ 2D structure (with H): benzene_showcase_2d_h.png")

# Step 3: Calculate Gasteiger charges
AllChem.ComputeGasteigerCharges(mol_3d)

charges = []
for atom in mol_3d.GetAtoms():
    charge = float(atom.GetProp('_GasteigerCharge'))
    charges.append(charge)

charge_min, charge_max = min(charges), max(charges)

# 2D structure with charge coloring
atom_colors = {}
for i, atom in enumerate(mol_3d.GetAtoms()):
    charge = charges[i]
    if charge_max != charge_min:
        norm_charge = (charge - charge_min) / (charge_max - charge_min)
    else:
        norm_charge = 0.5
    
    # Red (negative) → White → Blue (positive)
    if norm_charge < 0.5:
        r, g, b = 1.0, norm_charge * 2, norm_charge * 2
    else:
        r, g, b = (1 - norm_charge) * 2, (1 - norm_charge) * 2, 1.0
    
    atom_colors[i] = (r, g, b, 1.0)

drawer3 = rdMolDraw2D.MolDraw2DCairo(1000, 800)
drawer3.DrawMolecule(mol_3d, highlightAtoms=list(atom_colors.keys()), 
                     highlightAtomColors=atom_colors)
drawer3.FinishDrawing()
drawer3.WriteDrawingText("benzene_showcase_charges.png")
print("✓ Charge distribution: benzene_showcase_charges.png")

# Step 4: Aromatic ring highlighting
ring_info = mol.GetRingInfo()
aromatic_rings = [ring for ring in ring_info.AtomRings() 
                  if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

# Highlight aromatic ring
ring_colors = {}
for ring in aromatic_rings:
    for idx in ring:
        ring_colors[idx] = (0.3, 0.5, 1.0, 0.8)  # Blue

drawer4 = rdMolDraw2D.MolDraw2DCairo(800, 800)
drawer4.DrawMolecule(mol, highlightAtoms=list(ring_colors.keys()),
                     highlightAtomColors=ring_colors)
drawer4.FinishDrawing()
drawer4.WriteDrawingText("benzene_showcase_aromatic.png")
print("✓ Aromatic system: benzene_showcase_aromatic.png")

# Step 5: Export to SDF for xyzrender
writer = Chem.SDWriter("benzene_showcase.sdf")
writer.write(mol_3d)
writer.close()
print("✓ SDF file: benzene_showcase.sdf")

# Step 6: Molecular descriptors
print("\n" + "=" * 60)
print("Molecular Descriptors")
print("=" * 60)
print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
print(f"TPSA: {Descriptors.TPSA(mol):.2f} Å²")
print(f"HBD: {Descriptors.NumHDonors(mol)}")
print(f"HBA: {Descriptors.NumHAcceptors(mol)}")
print(f"Rotatable bonds: {Descriptors.NumRotatableBonds(mol)}")
print(f"Aromatic rings: {len(aromatic_rings)}")

# Step 7: Charge analysis
print("\n" + "=" * 60)
print("Gasteiger Charges")
print("=" * 60)
print("C atoms (negative):")
for i, atom in enumerate(mol_3d.GetAtoms()):
    if atom.GetSymbol() == 'C':
        print(f"  C[{i}]: {charges[i]:+.3f}")

print("\nH atoms (positive):")
for i, atom in enumerate(mol_3d.GetAtoms()):
    if atom.GetSymbol() == 'H':
        print(f"  H[{i}]: {charges[i]:+.3f}")

print("\n" + "=" * 60)
print("✅ All visualizations generated!")
print("=" * 60)
print("\nFiles created:")
print("  - benzene_showcase_2d.png")
print("  - benzene_showcase_2d_h.png")
print("  - benzene_showcase_charges.png")
print("  - benzene_showcase_aromatic.png")
print("  - benzene_showcase.sdf")
print("\nNext: Use xyzrender for 3D visualization")
