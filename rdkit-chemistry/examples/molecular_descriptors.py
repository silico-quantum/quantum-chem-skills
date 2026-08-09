#!/usr/bin/env python3
"""
RDKit molecular feature calculation and visualization.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import json

# Use DMAC-TRZ as the example
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== Molecular Feature Calculation: {NAME} ===\n")

# Build the molecule
mol = Chem.MolFromSmiles(SMILES)
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)

print("=" * 60)
print("1. Basic Molecular Information")
print("=" * 60)
print(f"Molecular formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Molecular weight: {Descriptors.ExactMolWt(mol):.2f} Da")
print(f"Total atom count: {mol_3d.GetNumAtoms()} (including H)")
print(f"Heavy atom count: {mol.GetNumAtoms()}")
print(f"Total bond count: {mol.GetNumBonds()}")

print("\n" + "=" * 60)
print("2. Physicochemical Properties")
print("=" * 60)
print(f"LogP (lipophilicity): {Descriptors.MolLogP(mol):.2f}")
print(f"TPSA (topological polar surface area): {Descriptors.TPSA(mol):.2f} Å²")
print(f"Rotatable bond count: {Descriptors.NumRotatableBonds(mol)}")
print(f"Hydrogen-bond donors: {Descriptors.NumHDonors(mol)}")
print(f"Hydrogen-bond acceptors: {Descriptors.NumHAcceptors(mol)}")

print("\n" + "=" * 60)
print("3. Electronic Properties")
print("=" * 60)

# Gasteiger charges
AllChem.ComputeGasteigerCharges(mol_3d)
charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol_3d.GetAtoms()]
print(f"Atomic charge range: [{min(charges):.3f}, {max(charges):.3f}]")
print(f"Mean charge: {np.mean(charges):.3f}")

# Identify the positive and negative charge centers
most_positive_idx = int(np.argmax(charges))
most_negative_idx = int(np.argmin(charges))
most_positive_atom = mol_3d.GetAtomWithIdx(most_positive_idx)
most_negative_atom = mol_3d.GetAtomWithIdx(most_negative_idx)

print(f"\nMost positive atom: {most_positive_atom.GetSymbol()} (idx {most_positive_idx}, charge={charges[most_positive_idx]:.3f})")
print(f"Most negative atom: {most_negative_atom.GetSymbol()} (idx {most_negative_idx}, charge={charges[most_negative_idx]:.3f})")

print("\n" + "=" * 60)
print("4. 3D Descriptors")
print("=" * 60)

# Molecular volume and surface area
volume = AllChem.ComputeMolVolume(mol_3d)
print(f"Molecular volume: {volume:.2f} Å³")

# Calculate the inertia tensor and molecular shape
conf = mol_3d.GetConformer()
coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
center = coords.mean(axis=0)

# Calculate the principal axes
coords_centered = coords - center
cov_matrix = np.cov(coords_centered.T)
eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
eigenvalues = np.real(eigenvalues)
eigenvalues = np.sort(eigenvalues)[::-1]

print(f"Principal-axis lengths (σ): {np.sqrt(eigenvalues[0]):.2f}, {np.sqrt(eigenvalues[1]):.2f}, {np.sqrt(eigenvalues[2]):.2f} Å")
print(f"Flatness (c/a): {np.sqrt(eigenvalues[2]/eigenvalues[0]):.3f}")

print("\n" + "=" * 60)
print("5. Conjugated-System Analysis")
print("=" * 60)

# Aromatic rings
num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
print(f"Aromatic ring count: {num_aromatic_rings}")

# Ring information
ring_info = mol.GetRingInfo()
print(f"Total ring count: {ring_info.NumRings()}")

# Approximate π-electron count
aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
print(f"Aromatic atom count: {aromatic_atoms}")

print("\n" + "=" * 60)
print("6. Donor-Acceptor Features (TADF-Relevant)")
print("=" * 60)

# Identify donor and acceptor groups
donor_pattern = Chem.MolFromSmarts('[N,n,O,o,S,s]')  # N, O, S
acceptor_pattern = Chem.MolFromSmarts('[N+!$(*[O-]),n+!$(*[O-]),O,o]')  # N+ or O

donor_matches = mol.GetSubstructMatches(donor_pattern)
acceptor_matches = mol.GetSubstructMatches(acceptor_pattern)

print(f"Potential donor sites: {len(donor_matches)}")
print(f"Potential acceptor sites: {len(acceptor_matches)}")

# Identify characteristic functional groups
functional_groups = {
    'tertiary amine': Chem.MolFromSmarts('[N;R0](C)(C)C'),
    'triazine': Chem.MolFromSmarts('c1ncncn1'),
    'carbazole': Chem.MolFromSmarts('c1ccc2c(c1)[nH]c1ccccc12'),
}

print("\nFunctional groups:")
for name, pattern in functional_groups.items():
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        print(f"  {name}: {len(matches)}")

print("\n" + "=" * 60)
print("7. Generate Visualizations")
print("=" * 60)

# 7.1 2D structure with atom indices
print("Generating 2D structure...")
drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
drawer.SetFontSize(0.8)
opts = drawer.drawOptions()
opts.addAtomIndices = True
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
drawer.WriteDrawingText(f"{NAME}_2d.png")
print(f"✓ Saved: {NAME}_2d.png")

# 7.2 2D structure with Gasteiger charges
print("Generating charge-distribution plot...")
for atom in mol_3d.GetAtoms():
    charge = float(atom.GetProp('_GasteigerCharge'))
    # Normalize to [-1, 1]
    norm_charge = np.clip(charge / 0.5, -1, 1)
    atom.SetProp('Note', f"{charge:.2f}")

drawer2 = rdMolDraw2D.MolDraw2DCairo(800, 600)
drawer2.DrawMolecule(mol_3d)
drawer2.FinishDrawing()
drawer2.WriteDrawingText(f"{NAME}_charges.png")
print(f"✓ Saved: {NAME}_charges.png")

# 7.3 3D conformation in SDF format
sdf_file = f"{NAME}_3d.sdf"
writer = Chem.SDWriter(sdf_file)
writer.write(mol_3d)
writer.close()
print(f"✓ Saved: {sdf_file}")

# 7.4 Export features to JSON
features = {
    "Molecular information": {
        "Molecular formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight": round(Descriptors.ExactMolWt(mol), 2),
        "Heavy atom count": mol.GetNumAtoms(),
    },
    "Physicochemical properties": {
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "Rotatable bond count": Descriptors.NumRotatableBonds(mol),
        "Hydrogen-bond donors": Descriptors.NumHDonors(mol),
        "Hydrogen-bond acceptors": Descriptors.NumHAcceptors(mol),
    },
    "Electronic properties": {
        "Charge range": [round(min(charges), 3), round(max(charges), 3)],
        "Most positive atom": f"{most_positive_atom.GetSymbol()} (idx {most_positive_idx})",
        "Most negative atom": f"{most_negative_atom.GetSymbol()} (idx {most_negative_idx})",
    },
    "3D descriptors": {
        "Molecular volume_Å3": round(volume, 2),
        "Principal-axis lengths_Å": [round(np.sqrt(e), 2) for e in eigenvalues],
        "Flatness": round(np.sqrt(eigenvalues[2]/eigenvalues[0]), 3),
    },
    "Conjugated system": {
        "Aromatic ring count": num_aromatic_rings,
        "Aromatic atom count": aromatic_atoms,
    },
    "TADF features": {
        "Donor site count": len(donor_matches),
        "Acceptor site count": len(acceptor_matches),
    }
}

json_file = f"{NAME}_features.json"
with open(json_file, 'w', encoding='utf-8') as f:
    json.dump(features, f, ensure_ascii=False, indent=2)
print(f"✓ Saved: {json_file}")

print("\n" + "=" * 60)
print("✅ Complete!")
print("=" * 60)
print("Generated files:")
print(f"  - {NAME}_2d.png (2D structure)")
print(f"  - {NAME}_charges.png (charge distribution)")
print(f"  - {NAME}_3d.sdf (3D conformation)")
print(f"  - {NAME}_features.json (feature data)")
print("=" * 60)
