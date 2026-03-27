#!/usr/bin/env python3
"""
PySCF Benzene Example: DFT Ground State + LR-TDDFT Excited States

Demonstrates:
  - Molecular definition with symmetry (D6h → D2h)
  - RKS-DFT with B3LYP/cc-pVDZ
  - Orbital analysis (HOMO/LUMO, gap)
  - LR-TDDFT excited states with oscillator strengths
  - Tamm-Dancoff Approximation (TDA)
  - Natural Transition Orbitals (NTO) analysis
  - Density fitting for efficiency
  - Solvent effects (PCM)

Usage:
    python benzene-dft-tddft.py
"""

from pyscf import gto, dft, tdscf, lo
import numpy as np

# ── 1. Molecular Definition ──────────────────────────────────────
# Benzene geometry (D6h symmetry, Å)
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

mol = gto.M(
    atom=benzene_xyz,
    basis='cc-pvdz',
    charge=0,
    spin=0,
    symmetry=True,     # PySCF auto-detects D2h (max Abelian subgroup)
    verbose=4,
)

print(f"\n{'='*60}")
print(f"Benzene: {mol.natm} atoms, {mol.nao_nr()} basis functions")
print(f"Point group: {mol.groupname}")
print(f"{'='*60}")

# ── 2. Ground State DFT ─────────────────────────────────────────
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.grids.atom_grid = (99, 590)   # (radial, angular) integration grid
mf.conv_tol = 1e-10              # SCF convergence threshold (Hartree)
mf.kernel()

print(f"\n--- Ground State ---")
print(f"B3LYP/cc-pVDZ energy: {mf.e_tot:.8f} Hartree ({mf.e_tot*27.2114:.4f} eV)")

# Orbital analysis
occ_idx = np.where(mf.mo_occ > 0)[0]
virt_idx = np.where(mf.mo_occ == 0)[0]
homo_energy = mf.mo_energy[occ_idx[-1]]
lumo_energy = mf.mo_energy[virt_idx[0]]
gap = lumo_energy - homo_energy

print(f"HOMO energy:  {homo_energy:.6f} Hartree ({homo_energy*27.2114:.4f} eV)")
print(f"LUMO energy:  {lumo_energy:.6f} Hartree ({lumo_energy*27.2114:.4f} eV)")
print(f"HOMO-LUMO gap: {gap:.6f} Hartree ({gap*27.2114:.4f} eV)")

# ── 3. LR-TDDFT Excited States ──────────────────────────────────
td = tdscf.TDDFT(mf)
td.nstates = 10
td.kernel()

print(f"\n--- LR-TDDFT Excited States (B3LYP/cc-pVDZ) ---")
print(f"{'State':>6} {'Energy (eV)':>12} {'λ (nm)':>10} {'f':>10} {'Sym':>8}")
print("-" * 50)
for i in range(td.nstates):
    e_ev = td.e[i] * 27.2114
    wl = 1240.0 / e_ev if e_ev > 0 else float('inf')
    print(f"  S{i+1:>2}  {e_ev:>10.4f}  {wl:>8.1f}  {td.oscillator_strength[i]:>10.4f}  "
          f"{td.get_symmetry(i)}")

# ── 4. TDA Comparison ───────────────────────────────────────────
td_tda = tdscf.TDA(mf)
td_tda.nstates = 6
td_tda.kernel()

print(f"\n--- TDA Excited States (B3LYP/cc-pVDZ) ---")
print(f"{'State':>6} {'Energy (eV)':>12} {'λ (nm)':>10} {'f':>10}")
print("-" * 42)
for i in range(td_tda.nstates):
    e_ev = td_tda.e[i] * 27.2114
    wl = 1240.0 / e_ev if e_ev > 0 else float('inf')
    print(f"  S{i+1:>2}  {e_ev:>10.4f}  {wl:>8.1f}  {td_tda.oscillator_strength[i]:>10.4f}")

# ── 5. NTO Analysis (first bright state) ────────────────────────
bright_idx = np.argmax(td.oscillator_strength)
weights, nto = td.get_nto(state=bright_idx)
print(f"\n--- NTO Analysis (brightest state S{bright_idx+1}) ---")
print(f"Oscillator strength: {td.oscillator_strength[bright_idx]:.4f}")
print(f"Excitation energy: {td.e[bright_idx]*27.2114:.4f} eV")
print(f"NTO weights: {weights[:5]}")  # top 5

# ── 6. Mulliken Population Analysis ─────────────────────────────
pop = mf.mulliken_pop(mol, mf.make_rdm1())
print(f"\n--- Mulliken Population ---")
for i, (sym, (q, s)) in enumerate(zip(mol.atom_symbol_list(), pop)):
    print(f"  {sym}{i+1:>2}: charge = {q:+.4f}, spin = {s:+.4f}")

# ── 7. Density Fitting (faster for larger systems) ──────────────
print(f"\n--- Density Fitting ---")
mf_df = dft.RKS(mol)
mf_df.xc = 'b3lyp'
mf_df = mf_df.density_fit(auxbasis='cc-pvdz-jkfit')
mf_df.kernel()
print(f"DF-B3LYP/cc-pVDZ energy: {mf_df.e_tot:.8f} Hartree")
print(f"Energy difference: {(mf_df.e_tot - mf.e_tot)*27.2114:.4f} eV")

# ── 8. Solvent Effect (PCM, cyclohexane) ────────────────────────
print(f"\n--- PCM Solvent (cyclohexane, ε=2.0) ---")
mf_pcm = dft.RKS(mol)
mf_pcm.xc = 'b3lyp'
mf_pcm = mf_pcm.ddPCM()
mf_pcm.with_solvent.eps = 2.0
mf_pcm.kernel()
print(f"B3LYP/cc-pVDZ(PCM) energy: {mf_pcm.e_tot:.8f} Hartree")
print(f"Solvation shift: {(mf_pcm.e_tot - mf.e_tot)*27.2114:.4f} eV")

# ── 9. Generate Cube Files for HOMO/LUMO ───────────────────────
from pyscf.tools import cubegen
cubegen.orbital(mol, 'benzene_homo.cube', mf.mo_coeff[:, occ_idx[-1]], mf.make_rdm1())
cubegen.orbital(mol, 'benzene_lumo.cube', mf.mo_coeff[:, virt_idx[0]], mf.make_rdm1())
print(f"\nCube files written: benzene_homo.cube, benzene_lumo.cube")
