#!/usr/bin/env python3
"""
Simplified charge and noncovalent-interaction visualization using RDKit only.

This qualitative analysis does not replace rigorous ESP or RDG calculations.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np

E_ANGSTROM_TO_DEBYE = 4.8032047

# Use DMAC-TRZ
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== Charge and noncovalent-interaction analysis: {NAME} ===\n")

# Build the molecule
mol = Chem.MolFromSmiles(SMILES)
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)

# Compute atomic charges
AllChem.ComputeGasteigerCharges(mol_3d)

def getCharge(atom):
    """Return the atomic charge."""
    try:
        return float(atom.GetProp('_GasteigerCharge'))
    except:
        return 0.0

print("=" * 60)
print("1. Identification of potential noncovalent-interaction sites")
print("=" * 60)

# Identify potential noncovalent-interaction sites
nci_sites = {
    "H-bond acceptors": [],
    "H-bond donors": [],
    "pi systems": [],
    "lone pairs": []
}

# Analyze the explicit-hydrogen molecule on which the charges were computed.
for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    symbol = atom.GetSymbol()
    charge = getCharge(atom)
    
    # H-bond acceptors (electronegative atoms with negative charge)
    if symbol in ['N', 'O', 'F'] and charge < -0.05:
        nci_sites["H-bond acceptors"].append((idx, symbol, charge))
    
    # H-bond donors (N or O bonded to H)
    if symbol in ['N', 'O'] and any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors()):
        nci_sites["H-bond donors"].append((idx, symbol, charge))
    
    # Lone pairs (nonbonding electrons on N, O, or S)
    if symbol in ['N', 'O', 'S']:
        nci_sites["lone pairs"].append((idx, symbol, charge))
    
    # Pi systems (aromatic atoms)
    if atom.GetIsAromatic():
        nci_sites["pi systems"].append((idx, symbol, charge))

# Print results
for site_type, sites in nci_sites.items():
    if sites:
        print(f"\n{site_type}:")
        for idx, symbol, charge in sites[:8]:  # Show the first eight
            print(f"  {symbol} (idx {idx}): charge {charge:.3f}")
        if len(sites) > 8:
            print(f"  ... {len(sites)} total")

print("\n" + "=" * 60)
print("2. Pi-pi stacking-propensity analysis")
print("=" * 60)

# Compute geometric features of aromatic rings
ring_info = mol_3d.GetRingInfo()
aromatic_rings = [ring for ring in ring_info.AtomRings() if all(mol_3d.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

print(f"Number of aromatic rings: {len(aromatic_rings)}")

# Compute each ring's center and normal vector
conf = mol_3d.GetConformer()
ring_centers = []
ring_normals = []

for ring in aromatic_rings:
    coords = np.array([list(conf.GetAtomPosition(i)) for i in ring])
    center = coords.mean(axis=0)
    ring_centers.append(center)
    
    # Compute the normal vector from the first three atoms
    if len(ring) >= 3:
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)
        ring_normals.append(normal)

# Analyze inter-ring angles as a proxy for pi-pi stacking propensity
print("\nRelative orientations of aromatic rings:")
for i in range(len(ring_normals)):
    for j in range(i+1, len(ring_normals)):
        cosine = np.clip(np.abs(np.dot(ring_normals[i], ring_normals[j])), 0.0, 1.0)
        angle = np.degrees(np.arccos(cosine))
        distance = np.linalg.norm(ring_centers[i] - ring_centers[j])
        
        if distance < 6.0:  # Rings within 6 Å
            print(f"  Ring {i} - ring {j}: angle {angle:.1f}°, distance {distance:.2f} Å")
            if angle < 20:
                print(f"    → Possible pi-pi stacking")

print("\n" + "=" * 60)
print("3. Electrostatic-potential proxy (atomic-charge coloring)")
print("=" * 60)

# Generate an atom-color map based on charge
charges = [getCharge(atom) for atom in mol_3d.GetAtoms()]
charge_min, charge_max = min(charges), max(charges)

print(f"Charge range: [{charge_min:.3f}, {charge_max:.3f}]")
print("Color map: red (negative) → white (neutral) → blue (positive)")

# Create atom colors
atom_colors = {}
for i, atom in enumerate(mol_3d.GetAtoms()):
    charge = charges[i]
    # Normalize to [0, 1]
    if charge_max != charge_min:
        norm_charge = (charge - charge_min) / (charge_max - charge_min)
    else:
        norm_charge = 0.5
    
    # RGB: red (negative) -> white (neutral) -> blue (positive)
    if norm_charge < 0.5:
        # Negative charge: red
        r = 1.0
        g = norm_charge * 2
        b = norm_charge * 2
    else:
        # Positive charge: blue
        r = (1 - norm_charge) * 2
        g = (1 - norm_charge) * 2
        b = 1.0
    
    atom_colors[i] = (r, g, b, 1.0)

# Generate the visualization
drawer = rdMolDraw2D.MolDraw2DCairo(1000, 800)
drawer.SetFontSize(0.9)
drawer.DrawMolecule(mol_3d, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors)
drawer.FinishDrawing()
drawer.WriteDrawingText(f"{NAME}_charge_surface.png")
print(f"✓ Saved: {NAME}_charge_surface.png")

print("\n" + "=" * 60)
print("4. Noncovalent-interaction-site visualization")
print("=" * 60)

# Mark key sites
hb_acceptors = [site[0] for site in nci_sites["H-bond acceptors"]]
pi_centers = []

# Select a representative atom from each aromatic ring
for i, ring in enumerate(aromatic_rings):
    if len(ring) > 0:
        # Select the ring's representative atom
        center_atom = ring[len(ring)//2]
        pi_centers.append(center_atom)

# Mark H-bond acceptors in red and pi-system representatives in blue
atom_colors2 = {}
for idx in hb_acceptors:
    atom_colors2[idx] = (1.0, 0.2, 0.2, 0.8)  # Red

for idx in pi_centers:
    if idx not in atom_colors2:
        atom_colors2[idx] = (0.2, 0.4, 1.0, 0.8)  # Blue

# Generate the visualization
drawer2 = rdMolDraw2D.MolDraw2DCairo(1000, 800)
drawer2.SetFontSize(0.9)
drawer2.DrawMolecule(mol_3d, highlightAtoms=list(atom_colors2.keys()), highlightAtomColors=atom_colors2)
drawer2.FinishDrawing()
drawer2.WriteDrawingText(f"{NAME}_nci_sites.png")
print(f"✓ Saved: {NAME}_nci_sites.png")
print(f"  Red: H-bond acceptors ({len(hb_acceptors)})")
print(f"  Blue: pi-system representatives ({len(pi_centers)})")

print("\n" + "=" * 60)
print("5. Molecular packing-propensity estimate")
print("=" * 60)

# Compute molecular shape factors
coords_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
cov_matrix = np.cov((coords_3d - coords_3d.mean(axis=0)).T)
eigenvalues = np.real(np.sort(np.linalg.eigvals(cov_matrix))[::-1])

L1, L2, L3 = np.sqrt(eigenvalues)
aspect_ratios = (L1/L3, L2/L3)

print(f"Molecular shape factors:")
print(f"  Length-to-width ratio (L1/L3): {aspect_ratios[0]:.2f}")
print(f"  Width-to-thickness ratio (L2/L3): {aspect_ratios[1]:.2f}")

if aspect_ratios[0] > 2 and aspect_ratios[1] > 1.5:
    print(f"  → Flat shape; may favor face-to-face packing")
elif aspect_ratios[0] > 3:
    print(f"  → Elongated shape; may favor edge-to-face packing")
else:
    print(f"  → Relatively spherical shape; weaker directional packing propensity")

# Estimate the dipole moment
print(f"\nDipole-moment estimate:")
dipole = np.zeros(3)
for i, atom in enumerate(mol_3d.GetAtoms()):
    charge = getCharge(atom)
    pos = np.array(list(conf.GetAtomPosition(i)))
    dipole += charge * pos

dipole_magnitude = np.linalg.norm(dipole) * E_ANGSTROM_TO_DEBYE
print(f"  |μ| ≈ {dipole_magnitude:.2f} Debye (estimated from Gasteiger charges)")
if dipole_magnitude > 3:
    print(f"  → A relatively large dipole moment may affect solid-state packing")

print("\n" + "=" * 60)
print("✅ Complete")
print("=" * 60)
print(f"Generated visualizations:")
print(f"  - {NAME}_charge_surface.png (atomic-charge coloring)")
print(f"  - {NAME}_nci_sites.png (potential noncovalent-interaction sites)")
print("\nLimitation: this simplified RDKit analysis does not replace rigorous ESP or RDG calculations; use Gaussian and Multiwfn for quantitative analysis.")
print("  Gaussian: #P B3LYP/6-31G* SCF=Tight")
print("  Multiwfn: function 5 (ESP), function 20 (RDG)")
print("=" * 60)
