#!/usr/bin/env python3
"""
RDKit molecular optimization and visualization demonstration.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import os

# Use the representative TADF molecule DMAC-TRZ (DMAC donor and TRZ acceptor).
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== Molecular optimization demonstration: {NAME} ===\n")

# Step 1: Build a molecule from SMILES.
print("Step 1: Build a 3D conformer from SMILES...")
mol = Chem.MolFromSmiles(SMILES)
if mol is None:
    print("❌ Failed to parse SMILES")
    exit(1)

mol = Chem.AddHs(mol)
print(f"✓ Molecule built: {mol.GetNumAtoms()} atoms")

# Step 2: Generate 3D coordinates.
print("\nStep 2: Generate 3D coordinates...")
res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
if res == -1:
    print("⚠️  Default embedding failed; trying random coordinates...")
    res = AllChem.EmbedMolecule(mol, AllChem.ETKDG(), useRandomCoords=True)

if res == -1:
    print("❌ Failed to generate 3D coordinates")
    exit(1)
print("✓ Generated 3D coordinates")

# Step 3: Optimize with the MMFF94 force field.
print("\nStep 3: MMFF94 force-field optimization...")
try:
    # Construct the force field.
    mmff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    if mmff is None:
        raise Exception("Failed to initialize the MMFF force field")

    # Optimize the geometry.
    initial_energy = mmff.CalcEnergy()
    print(f"  Initial energy: {initial_energy:.2f} kcal/mol")

    converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    if converged == 0:
        final_energy = mmff.CalcEnergy()
        print(f"  Final energy: {final_energy:.2f} kcal/mol")
        print(f"  Energy decrease: {initial_energy - final_energy:.2f} kcal/mol")
        print("✓ MMFF optimization converged")
    else:
        print("⚠️  MMFF optimization did not fully converge")

except Exception as e:
    print(f"⚠️  MMFF failed: {e}")
    print("  Trying the UFF force field...")
    try:
        AllChem.UFFOptimizeMolecule(mol)
        print("✓ UFF optimization completed")
    except Exception as e2:
        print(f"❌ UFF also failed: {e2}")
        exit(1)

# Step 4: Export XYZ and SDF files.
print("\nStep 4: Export files...")
xyz_file = f"{NAME}.xyz"
sdf_file = f"{NAME}.sdf"

# XYZ file.
Chem.MolToXYZFile(mol, xyz_file)
print(f"✓ Saved: {xyz_file}")

# SDF file, including bond information.
writer = Chem.SDWriter(sdf_file)
writer.write(mol)
writer.close()
print(f"✓ Saved: {sdf_file}")

# Read the XYZ file for the preview used before rendering.
with open(xyz_file, 'r') as f:
    xyz_content = f.read()
    print("\nXYZ file preview:")
    print("```")
    lines = xyz_content.split('\n')
    print(lines[0])  # Atom count.
    print(lines[1])  # Comment line.
    for i in range(2, min(7, len(lines))):
        print(lines[i])
    if len(lines) > 7:
        print("...")
    print("```")

# Step 5: Render with xyzrender.
print("\nStep 5: Render with xyzrender...")
png_file = f"{NAME}.png"

try:
    # Render the SDF file to retain bond information.
    result = subprocess.run(
        ['xyzrender', sdf_file, '-o', png_file, '--transparent', '--bo'],
        capture_output=True,
        text=True,
        timeout=30
    )

    if result.returncode == 0 and os.path.exists(png_file):
        file_size = os.path.getsize(png_file)
        print(f"✓ Rendering succeeded: {png_file} ({file_size} bytes)")
    else:
        print("❌ Rendering failed:")
        print(result.stderr)
        exit(1)

except FileNotFoundError:
    print("❌ xyzrender is not installed")
    print("  Install with: pip install xyzrender")
    exit(1)
except Exception as e:
    print(f"❌ Rendering error: {e}")
    exit(1)

print("\n" + "="*50)
print("✅ Complete!")
print(f"   Molecule: {NAME}")
print(f"   Atoms: {mol.GetNumAtoms()}")
print(f"   XYZ: {xyz_file}")
print(f"   PNG: {png_file}")
print("="*50)
