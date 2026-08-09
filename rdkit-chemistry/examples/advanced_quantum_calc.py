#!/usr/bin/env python3
"""Illustrative benzene DFT, orbital, and vertical delta-SCF workflow."""

from _legacy_guard import refuse_legacy_example

refuse_legacy_example(__file__)

import json
import time

import numpy as np
from pyscf import dft, gto
from pyscf.tools import cubegen

print("=" * 60)
print("Illustrative Benzene Quantum Chemistry Calculation")
print("=" * 60)
print("Method: DFT (B3LYP/6-31G)")
print("Purpose: molecular orbitals and vertical delta-SCF descriptors")
print("=" * 60)

# Use a simple test molecule: benzene (12 atoms and an even electron count).
xyz = """
C     0.000000   1.400000   0.000000
C     1.212400   0.700000   0.000000
C     1.212400  -0.700000   0.000000
C     0.000000  -1.400000   0.000000
C    -1.212400  -0.700000   0.000000
C    -1.212400   0.700000   0.000000
H     0.000000   2.490000   0.000000
H     2.156000   1.245000   0.000000
H     2.156000  -1.245000   0.000000
H     0.000000  -2.490000   0.000000
H    -2.156000  -1.245000   0.000000
H    -2.156000   1.245000   0.000000
"""

print("\nUsing benzene as the example molecule (C6H6)")
print("Number of atoms: 12 (6 C + 6 H)")

# Build the PySCF molecule object.
# Build it once first to obtain the electron count.
mol_temp = gto.M(atom=xyz.strip(), basis='6-31G', charge=0, verbose=0)
nelectron = mol_temp.nelectron
spin = nelectron % 2  # An odd electron count requires an open-shell calculation.

# Rebuild with the correct spin.
mol = gto.M(
    atom=xyz.strip(),
    basis='6-31G',  # Use 6-31G to balance accuracy and speed.
    charge=0,
    spin=spin,
    verbose=4
)

# Check the electron count and spin.
nelectron = mol.nelectron
spin = nelectron % 2  # An odd electron count requires an open-shell calculation.
print(f"\nMolecular information:")
print(f"  Number of atoms: {mol.natm}")
print(f"  Number of electrons: {nelectron}")
print(f"  Spin: {spin}")
print(f"  Basis set: {mol.basis}")

# Step 1: DFT single-point calculation
print("\n" + "=" * 60)
print("Step 1: DFT Single-Point Calculation (B3LYP)")
print("=" * 60)

start_time = time.time()

# Use the B3LYP functional, selecting RKS or UKS according to the spin.
if spin == 0:
    mf = dft.RKS(mol)
else:
    mf = dft.UKS(mol)
mf.xc = 'B3LYP'
mf.conv_tol = 1e-8
mf.kernel()
if not mf.converged:
    raise RuntimeError("Neutral B3LYP calculation did not converge")

elapsed = time.time() - start_time
print(f"✓ SCF converged: {elapsed:.2f} seconds")
print(f"  Total energy: {mf.e_tot:.6f} Hartree")
print(f"  ({mf.e_tot * 627.509:.2f} kcal/mol)")

# Step 2: Molecular-orbital analysis
print("\n" + "=" * 60)
print("Step 2: Molecular-Orbital Analysis")
print("=" * 60)

# Obtain the molecular-orbital coefficients.
mo_coeff = mf.mo_coeff
mo_occ = mf.mo_occ
mo_energy = mf.mo_energy

# Locate the HOMO and LUMO.
if spin == 0:
    # Closed shell
    nocc = int(mol.nelectron / 2)  # Number of occupied orbitals
    homo_idx = nocc - 1
    lumo_idx = nocc
else:
    # Open shell: locate the last occupied orbital.
    if isinstance(mo_occ, tuple):
        # For UKS, mo_occ is (alpha, beta).
        occ_alpha = mo_occ[0]
        occ_beta = mo_occ[1]
        nocc_alpha = np.sum(occ_alpha > 0)
        nocc_beta = np.sum(occ_beta > 0)
        # Use the alpha orbitals.
        homo_idx = nocc_alpha - 1
        lumo_idx = nocc_alpha
    else:
        nocc = np.sum(mo_occ > 0)
        homo_idx = nocc - 1
        lumo_idx = nocc

homo_energy = mo_energy[homo_idx] if not isinstance(mo_energy, tuple) else mo_energy[0][homo_idx]
lumo_energy = mo_energy[lumo_idx] if not isinstance(mo_energy, tuple) else mo_energy[0][lumo_idx]
gap = lumo_energy - homo_energy

print(f"HOMO energy: {homo_energy:.4f} Hartree ({homo_energy * 27.2114:.2f} eV)")
print(f"LUMO energy: {lumo_energy:.4f} Hartree ({lumo_energy * 27.2114:.2f} eV)")
print(f"HOMO-LUMO Gap: {gap:.4f} Hartree ({gap * 27.2114:.2f} eV)")

# Orbital-composition analysis
print("\nOrbital-composition analysis:")
print("  HOMO (largest contributions):")
homo_coeff = mo_coeff[:, homo_idx] if not isinstance(mo_coeff, tuple) else mo_coeff[0][:, homo_idx]
# Locate the atoms with the largest contributions.
contributions = []
for i in range(mol.natm):
    ao_slice = mol.aoslice_by_atom()[i]
    coeff_sq = (homo_coeff[ao_slice[2]:ao_slice[3]] ** 2).sum()
    contributions.append((i, mol.atom_pure_symbol(i), coeff_sq))

contributions.sort(key=lambda x: x[2], reverse=True)
for i, symbol, contrib in contributions[:3]:
    print(f"    {symbol} (idx {i}): {contrib * 100:.1f}%")

print("  LUMO (largest contributions):")
lumo_coeff = mo_coeff[:, lumo_idx] if not isinstance(mo_coeff, tuple) else mo_coeff[0][:, lumo_idx]
contributions = []
for i in range(mol.natm):
    ao_slice = mol.aoslice_by_atom()[i]
    coeff_sq = (lumo_coeff[ao_slice[2]:ao_slice[3]] ** 2).sum()
    contributions.append((i, mol.atom_pure_symbol(i), coeff_sq))

contributions.sort(key=lambda x: x[2], reverse=True)
for i, symbol, contrib in contributions[:3]:
    print(f"    {symbol} (idx {i}): {contrib * 100:.1f}%")

# Step 3: Electron-density analysis
print("\n" + "=" * 60)
print("Step 3: Electron Density and Mulliken Charges")
print("=" * 60)

# Calculate the electron density.
dm = mf.make_rdm1()

# Calculate Mulliken charges.
S = mol.intor('int1e_ovlp')
C = mf.mo_coeff
occ = mf.mo_occ

# Mulliken population analysis
dm_mulliken = np.dot(C * occ, C.T)
P_mul = np.dot(dm_mulliken, S)
mulliken_charges = []
for i in range(mol.natm):
    ao_slice = mol.aoslice_by_atom()[i]
    nuc_charge = mol.atom_charge(i)
    elec_charge = np.trace(P_mul[ao_slice[2]:ao_slice[3], ao_slice[2]:ao_slice[3]])
    mulliken_charges.append(nuc_charge - elec_charge)

print("Mulliken atomic charges:")
for i in range(mol.natm):
    symbol = mol.atom_pure_symbol(i)
    print(f"  {symbol} (idx {i}): {mulliken_charges[i]:+.3f}")

# Calculate the dipole moment.
dipole = mf.dip_moment()
dipole_magnitude = np.linalg.norm(dipole)
print(f"\nDipole moment: {dipole_magnitude:.3f} Debye")
print(f"  Components: [{dipole[0]:.3f}, {dipole[1]:.3f}, {dipole[2]:.3f}]")

# Step 4: vertical delta-SCF reactivity descriptors
print("\n" + "=" * 60)
print("Step 4: Vertical Delta-SCF Reactivity Descriptors")
print("=" * 60)

print("Calculating N-1, N, and N+1 single-point energies...")

# N-1 electron system (cation)
mol_plus = gto.M(
    atom=xyz.strip(),
    basis='6-31G',
    charge=+1,
    spin=1,  # Open shell
    verbose=0
)
mf_plus = dft.UKS(mol_plus)
mf_plus.xc = 'B3LYP'
mf_plus.kernel()
if not mf_plus.converged:
    raise RuntimeError("Cation B3LYP calculation did not converge")

# N+1 electron system (anion)
mol_minus = gto.M(
    atom=xyz.strip(),
    basis='6-31G',
    charge=-1,
    spin=1,  # Open shell
    verbose=0
)
mf_minus = dft.UKS(mol_minus)
mf_minus.xc = 'B3LYP'
mf_minus.kernel()
if not mf_minus.converged:
    raise RuntimeError("Anion B3LYP calculation did not converge")

print(f"✓ Neutral-molecule energy: {mf.e_tot:.6f} Hartree")
print(f"✓ Cation energy: {mf_plus.e_tot:.6f} Hartree")
print(f"✓ Anion energy: {mf_minus.e_tot:.6f} Hartree")

# Calculate the ionization potential and electron affinity.
IP = mf_plus.e_tot - mf.e_tot  # Ionization potential
EA = mf.e_tot - mf_minus.e_tot  # Electron affinity

print(f"\nVertical ionization potential (IP): {IP:.4f} Hartree ({IP * 27.2114:.2f} eV)")
print(f"Vertical electron affinity (EA): {EA:.4f} Hartree ({EA * 27.2114:.2f} eV)")
print(f"Electronegativity (χ): {(IP + EA) / 2 * 27.2114:.2f} eV")
print(f"Chemical hardness (η): {(IP - EA) / 2 * 27.2114:.2f} eV")

print(
    "Atom-condensed Fukui functions are not reported: they require a "
    "consistent population analysis of the converged N-1, N, and N+1 densities."
)

# Step 5: Generate cube files for visualization.
print("\n" + "=" * 60)
print("Step 5: Generate Cube Files")
print("=" * 60)

print("Generating molecular-orbital cube files...")
# HOMO
cubegen.orbital(mol, "benzene_HOMO.cube", mo_coeff[:, homo_idx])
print("✓ benzene_HOMO.cube")

# LUMO
cubegen.orbital(mol, "benzene_LUMO.cube", mo_coeff[:, lumo_idx])
print("✓ benzene_LUMO.cube")

# Electron density
cubegen.density(mol, "benzene_density.cube", dm)
print("✓ benzene_density.cube")

# Summary
print("\n" + "=" * 60)
print("✅ Calculation complete!")
print("=" * 60)
print(f"Total elapsed time: {time.time() - start_time:.2f} seconds")
print("\nGenerated files:")
print("  - benzene_HOMO.cube")
print("  - benzene_LUMO.cube")
print("  - benzene_density.cube")
print("\nNext steps:")
print("  1. Analyze the cube files with Multiwfn")
print("  2. Generate an ESP isosurface plot")
print("  3. Calculate excited states with TDDFT")
print("=" * 60)

# Save the result summary.
summary = {
    "method": "B3LYP/6-31G",
    "total_energy_Hartree": mf.e_tot,
    "HOMO_eV": homo_energy * 27.2114,
    "LUMO_eV": lumo_energy * 27.2114,
    "Gap_eV": gap * 27.2114,
    "dipole_moment_Debye": dipole_magnitude,
    "IP_eV": IP * 27.2114,
    "EA_eV": EA * 27.2114,
}

with open("benzene_quantum_calc_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nResult summary saved to benzene_quantum_calc_summary.json")
