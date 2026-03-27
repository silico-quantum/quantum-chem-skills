#!/usr/bin/env python3
"""PySCF Benzene: DFT + TDDFT (PySCF 2.12.1)
Tested: 2026-03-27 on macOS arm64, PySCF 2.12.1, Python 3.14
"""
from pyscf import gto, dft, tdscf
from pyscf.tools import cubegen
from pyscf import solvent
import numpy as np

benzene_xyz = """
C   0.000000   1.395000   0.000000
C   1.208543   0.697500   0.000000
C   1.208543  -0.697500   0.000000
C   0.000000  -1.395000   0.000000
C  -1.208543  -0.697500   0.000000
C  -1.208543   0.697500   0.000000
H   0.000000   2.479000   0.000000
H   2.150000   1.239500   0.000000
H   2.150000  -1.239500   0.000000
H   0.000000  -2.479000   0.000000
H  -2.150000  -1.239500   0.000000
H  -2.150000   1.239500   0.000000
"""

mol = gto.M(atom=benzene_xyz, basis='cc-pvdz', charge=0, spin=0,
            symmetry=True, verbose=4)
print(f"\n{'='*60}")
print(f"Benzene: {mol.natm} atoms, {mol.nao_nr()} AOs, {mol.groupname}")
print(f"{'='*60}")

# === 1. Ground State DFT ===
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.grids.atom_grid = (99, 590)
mf.conv_tol = 1e-10
mf.kernel()

occ = np.where(mf.mo_occ > 0)[0]
virt = np.where(mf.mo_occ == 0)[0]
homo_e, lumo_e = mf.mo_energy[occ[-1]], mf.mo_energy[virt[0]]

print(f"\n--- Ground State ---")
print(f"Energy: {mf.e_tot:.8f} Ha ({mf.e_tot*27.2114:.4f} eV)")
print(f"HOMO:  {homo_e:.6f} Ha ({homo_e*27.2114:.4f} eV)")
print(f"LUMO:  {lumo_e:.6f} Ha ({lumo_e*27.2114:.4f} eV)")
print(f"Gap:   {lumo_e-homo_e:.6f} Ha ({(lumo_e-homo_e)*27.2114:.4f} eV)")

# === 2. LR-TDDFT ===
td = tdscf.TDDFT(mf)
td.nstates = 10
td.kernel()
osc = td.oscillator_strength()

print(f"\n--- LR-TDDFT (B3LYP/cc-pVDZ) ---")
print(f"{'State':>6} {'E(eV)':>10} {'nm':>8} {'f':>10}")
print("-" * 38)
for i in range(td.nstates):
    e = td.e[i] * 27.2114
    print(f"  S{i+1:>2} {e:>8.4f} {1240.0/e:>7.1f} {osc[i]:>10.4f}")

# === 3. TDA ===
td2 = tdscf.TDA(mf)
td2.nstates = 6
td2.kernel()
osc2 = td2.oscillator_strength()

print(f"\n--- TDA (B3LYP/cc-pVDZ) ---")
print(f"{'State':>6} {'E(eV)':>10} {'nm':>8} {'f':>10}")
print("-" * 38)
for i in range(td2.nstates):
    e = td2.e[i] * 27.2114
    print(f"  S{i+1:>2} {e:>8.4f} {1240.0/e:>7.1f} {osc2[i]:>10.4f}")

# === 4. NTO ===
b = int(np.argmax(osc))
w, _ = td.get_nto(state=b+1)
print(f"\n--- NTO (S{b+1}, f={osc[b]:.4f}) ---")
print(f"Weights: {w[:5]}")

# === 5. Mulliken ===
pop = mf.mulliken_pop()
print(f"\n--- Mulliken ---")
for i in range(mol.natm):
    print(f"  {mol.atom_symbol(i)}{i+1:>2}: q={pop[0][i]:+.4f}")

# === 6. Density Fitting ===
mf_df = dft.RKS(mol)
mf_df.xc = 'b3lyp'
mf_df = mf_df.density_fit(auxbasis='cc-pvdz-jkfit')
mf_df.kernel()
print(f"\n--- DF: {mf_df.e_tot:.8f} Ha, diff={(mf_df.e_tot-mf.e_tot)*27.2114:.4f} eV ---")

# === 7. PCM (ddPCM) ===
try:
    mf_pcm = solvent.ddPCM(dft.RKS(mol))
    mf_pcm.xc = 'b3lyp'
    mf_pcm.with_solvent.eps = 2.0
    mf_pcm.kernel()
    print(f"--- PCM: {mf_pcm.e_tot:.8f} Ha, shift={(mf_pcm.e_tot-mf.e_tot)*27.2114:.4f} eV ---")
except Exception as e:
    print(f"--- PCM skipped: {e} ---")

# === 8. Cube files ===
import os
# Cube files need non-symmetrized mol for cubegen compatibility
mol_nosym = gto.M(atom=benzene_xyz, basis='cc-pvdz', charge=0, spin=0, symmetry=False, verbose=0)
mf_nosym = dft.RKS(mol_nosym)
mf_nosym.xc = 'b3lyp'
mf_nosym.kernel()
occ_n = np.where(mf_nosym.mo_occ > 0)[0]
virt_n = np.where(mf_nosym.mo_occ == 0)[0]
outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'examples', 'pyscf')
os.makedirs(outdir, exist_ok=True)
homo_path = os.path.join(outdir, 'benzene_homo.cube')
lumo_path = os.path.join(outdir, 'benzene_lumo.cube')
cubegen.orbital(mol_nosym, homo_path, mf_nosym.mo_coeff[:, occ_n[-1]])
cubegen.orbital(mol_nosym, lumo_path, mf_nosym.mo_coeff[:, virt_n[0]])
print(f"\nCube files: {homo_path}, {lumo_path}")
