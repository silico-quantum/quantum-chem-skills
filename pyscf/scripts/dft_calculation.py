"""
Complete PySCF DFT example.

Demonstrates a water-molecule DFT calculation, excited-state analysis, and
selected molecular properties.
"""

from pyscf import gto, dft, tdscf, lib
import numpy as np

def water_dft():
    """Run the water-molecule DFT example."""

    print("="*60)
    print("1. Molecule definition")
    print("="*60)

    # Define the water molecule.
    mol = gto.M(
        atom='''
        O  0.000000  0.000000  0.117790
        H  0.000000  0.755453 -0.471161
        H  0.000000 -0.755453 -0.471161
        ''',
        basis='cc-pvdz',
        charge=0,
        spin=0,
        unit='Ang',
        symmetry='c2v'
    )

    print(f"Molecule: {mol.atom[0][0]}{'H'*2}")
    print(f"Electrons: {mol.nelectron}")
    print(f"Atomic orbitals: {mol.nao}")
    print(f"Symmetry: {mol.symmetry}")
    print()

    print("="*60)
    print("2. DFT calculation (B3LYP)")
    print("="*60)

    # RKS-DFT
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.grids.atom_grid = (75, 302)

    # Run the SCF calculation.
    e_tot = mf.kernel()

    print(f"\nGround-state energy: {e_tot:.10f} Hartree")
    print(f"Energy: {e_tot*27.2114:.6f} eV")
    print(f"Converged: {mf.converged}")
    print()

    # Molecular-orbital analysis.
    print("="*60)
    print("3. Molecular-orbital analysis")
    print("="*60)

    # Orbital energies and occupations.
    mo_energy = mf.mo_energy
    mo_occ = mf.mo_occ

    # HOMO
    homo_idx = np.where(mo_occ > 0)[0][-1]
    homo_energy = mo_energy[homo_idx]

    # LUMO
    lumo_idx = np.where(mo_occ == 0)[0][0]
    lumo_energy = mo_energy[lumo_idx]

    # Orbital-energy gap.
    gap = lumo_energy - homo_energy

    print(f"HOMO energy: {homo_energy:.6f} Hartree ({homo_energy*27.2114:.3f} eV)")
    print(f"LUMO energy: {lumo_energy:.6f} Hartree ({lumo_energy*27.2114:.3f} eV)")
    print(f"Orbital-energy gap: {gap:.6f} Hartree ({gap*27.2114:.3f} eV)")
    print()

    # Occupied orbitals.
    print("Occupied orbitals:")
    for i in range(mol.nao):
        if mo_occ[i] > 0:
            print(f"  MO {i+1:2d}: E = {mo_energy[i]:10.6f} Ha,  occ = {mo_occ[i]:.2f}")
    print()

    # Virtual orbitals.
    print("Virtual orbitals (first five):")
    for i in range(lumo_idx, min(lumo_idx+5, mol.nao)):
        print(f"  MO {i+1:2d}: E = {mo_energy[i]:10.6f} Ha,  occ = {mo_occ[i]:.2f}")
    print()

    print("="*60)
    print("4. LR-TDDFT excited states")
    print("="*60)

    # TDDFT calculation.
    td = tdscf.TDDFT(mf)
    td.nstates = 6
    td.kernel()
    oscillator_strengths = td.oscillator_strength()

    # Analyze the excited states.
    print("\nExcited-state analysis:")
    print(f"{'State':<6} {'Energy (eV)':<12} {'Wavelength (nm)':<16} {'Osc. strength':<14}")
    print("-"*60)

    for i, excitation_energy in enumerate(td.e):
        e_ev = excitation_energy * 27.2114
        wavelength = 1240.0 / e_ev
        f = oscillator_strengths[i]
        print(f"{i+1:<4} {e_ev:<10.3f} {wavelength:<12.2f} {f:<12.3f}")
    print()

    # Brightest computed excited state.
    strongest = np.argmax(oscillator_strengths)
    print(f"Brightest excited state: state {strongest+1}, oscillator strength = {oscillator_strengths[strongest]:.3f}")
    print()

    print("="*60)
    print("5. NTO analysis (first excited state)")
    print("="*60)

    # Natural transition orbitals.
    weights, nto_coeff = td.get_nto(state=1)
    nocc = np.count_nonzero(mf.mo_occ > 0)
    occupied_ntos = nto_coeff[:, :nocc]
    virtual_ntos = nto_coeff[:, nocc:]

    print(f"Leading transition weight: {weights.max():.6f}")
    print(f"Occupied NTOs: {occupied_ntos.shape[1]}")
    print(f"Virtual NTOs: {virtual_ntos.shape[1]}")
    print()

    print("="*60)
    print("6. Dipole moment")
    print("="*60)

    # Electric dipole moment.
    dip = mf.dip_moment(unit='Debye')
    print("Electric dipole moment (Debye):")
    print(f"  μ_x = {dip[0]:.6f}")
    print(f"  μ_y = {dip[1]:.6f}")
    print(f"  μ_z = {dip[2]:.6f}")
    print(f"  |μ| = {np.linalg.norm(dip):.6f} Debye")
    print()

    print("="*60)
    print("7. Charge-population analysis")
    print("="*60)

    # Mulliken population analysis. The second return value contains one
    # atomic charge per atom; the first contains AO populations.
    _, atomic_charges = mf.mulliken_pop(mol, mf.make_rdm1())

    print("Mulliken atomic charges:")
    for i in range(mol.natm):
        atom = mol.atom_pure_symbol(i)
        charge = atomic_charges[i]
        print(f"  {atom}: {charge:.6f} e")
    print()

    print("="*60)
    print("8. Functional comparison")
    print("="*60)

    functionals = ['pbe', 'b3lyp', 'wb97x-d', 'cam-b3lyp']
    results = {}

    for xc in functionals:
        mf_xc = dft.RKS(mol)
        mf_xc.xc = xc
        e_xc = mf_xc.kernel()
        results[xc] = e_xc
        print(f"{xc:12s}: {e_xc:.10f} Hartree")

    # Energy differences.
    print("\nRelative energies (with respect to PBE):")
    e_ref = results['pbe']
    for xc, e in results.items():
        delta = (e - e_ref) * 627.509  # Hartree to kcal/mol.
        print(f"{xc:12s}: {delta:10.3f} kcal/mol")
    print()

    print("="*60)
    print("9. Save results")
    print("="*60)

    # Save a checkpoint.
    mf.chkfile = 'water_dft.chk'
    mf.dump_chk(mf.chkfile)

    print(f"Checkpoint file: {mf.chkfile}")
    print("Contents include MO coefficients, orbital energies, and the density matrix.")
    print()

    return mf, td

if __name__ == '__main__':
    # Set the output level.
    lib.logger.TIMER_LEVEL = 3
    lib.logger.INFO = True

    # Run the calculation.
    mf, td = water_dft()

    print("="*60)
    print("Calculation complete!")
    print("="*60)
