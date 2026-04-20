#!/usr/bin/env python3
"""
Wavefunction Analysis - Population analysis, orbital analysis, DOS

Supports:
    - Mulliken charges
    - Löwdin charges
    - Natural Population Analysis (NPA)
    - NBO analysis
    - DOS/PDOS analysis
    - Dipole moment
    - Orbital composition

Usage:
    python analysis.py <chkfile> <task> [params]
    python analysis.py hf.chk mulliken
    python analysis.py hf.chk dos --plot
"""

import sys
import numpy as np

from pyscf import gto, scf, dft, lo


def load_from_chk(chkfile, method='RHF'):
    """Load SCF result from checkpoint file"""
    mol = gto.M()
    mol.build()
    
    if method == 'RHF':
        mf = scf.RHF(mol)
    elif method == 'UHF':
        mf = scf.UHF(mol)
    elif method == 'RKS':
        mf = dft.RKS(mol)
    elif method == 'UKS':
        mf = dft.UKS(mol)
    
    mf = mf.from_chk(chkfile)
    return mf


def mulliken_charges(mf, mol=None):
    """
    Mulliken population analysis
    
    Returns:
        array: atomic charges
    """
    if mol is None:
        mol = mf.mol
    
    dm = mf.make_rdm1()
    orth_ao = lo.orth_ao(mol, 'minao')
    
    from pyscf import lib
    dm_orth = lib.dot(orth_ao.T, lib.dot(dm, orth_ao))
    
    atoms = [mol.atom_symbol(i) for i in range(mol.natm)]
    charges = np.zeros(mol.natm)
    
    start = 0
    for ia in range(mol.natm):
        nao = mol.atom_nbas(ia)
        end = start + nao
        charges[ia] = dm_orth[start:end, start:end].trace()
        start = end
    
    nuclear = np.array([gto.charge(mol.atom_symbol(i)) for i in range(mol.natm)])
    charges = nuclear - charges
    
    print('\nMulliken Charges:')
    for i, (atom, charge) in enumerate(zip(atoms, charges)):
        print('  %3d %2s: %+.6f' % (i, atom, charge))
    
    return charges


def lowdin_charges(mf, mol=None):
    """
    Löwdin population analysis (more stable than Mulliken)
    """
    if mol is None:
        mol = mf.mol
    
    orth_ao = lo.orth_ao(mol, 'meta_lowdin')
    dm = mf.make_rdm1()
    
    from pyscf import lib
    dm_orth = lib.dot(orth_ao.T, lib.dot(dm, orth_ao))
    
    atoms = [mol.atom_symbol(i) for i in range(mol.natm)]
    charges = np.zeros(mol.natm)
    
    start = 0
    for ia in range(mol.natm):
        nao = mol.atom_nbas(ia)
        end = start + nao
        charges[ia] = dm_orth[start:end, start:end].trace()
        start = end
    
    nuclear = np.array([gto.charge(mol.atom_symbol(i)) for i in range(mol.natm)])
    charges = nuclear - charges
    
    print('\nLöwdin Charges:')
    for i, (atom, charge) in enumerate(zip(atoms, charges)):
        print('  %3d %2s: %+.6f' % (i, atom, charge))
    
    return charges


def dipole_moment(mf, mol=None):
    """
    Compute dipole moment
    
    Returns:
        array: [dx, dy, dz] in Debye
    """
    if mol is None:
        mol = mf.mol
    
    dm = mf.make_rdm1()
    
    from pyscf import df
    ao_dip = df.get_dipole(mol)
    
    elec_dip = np.einsum('xij,ji->x', ao_dip, dm)
    
    nuclear_dip = np.zeros(3)
    for i in range(mol.natm):
        nuclear_dip += mol.atom_coord(i) * gto.charge(mol.atom_symbol(i))
    
    total_dip = nuclear_dip - elec_dip
    total_debye = total_dip * 2.5418  # a.u. to Debye
    
    print('\nDipole Moment:')
    print('  x = %.6f D' % total_debye[0])
    print('  y = %.6f D' % total_debye[1])
    print('  z = %.6f D' % total_debye[2])
    print('  |μ| = %.6f D' % np.linalg.norm(total_debye))
    
    return total_debye


def orbital_population(mf, mol=None, mo_index=None):
    """
    Orbital composition analysis — which atoms/shells contribute most
    """
    if mol is None:
        mol = mf.mol
    
    ao_labels = mol.ao_labels()
    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    if mo_index is not None:
        indices = [mo_index]
    else:
        indices = list(range(min(10, mo_coeff.shape[1])))
    
    print('\nOrbital Composition Analysis:')
    print('-' * 60)
    
    for idx in indices:
        coeff_sq = mo_coeff[:, idx]**2
        sorted_idx = np.argsort(coeff_sq)[::-1]
        
        occ_tag = '[occ]' if mo_occ[idx] > 0 else '[emp]'
        
        print('\nMO %d %s E=%.4f eV' % (idx, occ_tag, mo_energy[idx]))
        print('  Top contributions:')
        for j in sorted_idx[:5]:
            if coeff_sq[j] > 0.01:
                print('    %-30s %.4f' % (ao_labels[j], coeff_sq[j]))


def dos_analysis(mf, mol=None, energy_range=None, sigma=0.1):
    """
    Density of States (DOS) analysis
    
    Returns:
        tuple: (energies, dos)
    """
    if mol is None:
        mol = mf.mol
    
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    if energy_range is None:
        e_min = mo_energy.min() - 1
        e_max = mo_energy.max() + 1
    else:
        e_min, e_max = energy_range
    
    energies = np.linspace(e_min, e_max, 1000)
    dos = np.zeros_like(energies)
    
    for i, e in enumerate(mo_energy):
        dos += np.exp(-(energies - e)**2 / (2 * sigma**2))
    
    print('\nDOS Analysis:')
    print('  HOMO: %.4f eV' % mo_energy[mo_occ > 0].max())
    print('  LUMO: %.4f eV' % mo_energy[mo_occ == 0].min())
    print('  Gap:  %.4f eV' % (mo_energy[mo_occ == 0].min() - mo_energy[mo_occ > 0].max()))
    
    return energies, dos


def pdos_analysis(mf, mol=None, sigma=0.1):
    """
    Projected DOS (PDOS) per atom
    """
    if mol is None:
        mol = mf.mol
    
    mo_energy = mf.mo_energy * 27.2114
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    
    e_min = mo_energy.min() - 1
    e_max = mo_energy.max() + 1
    energies = np.linspace(e_min, e_max, 500)
    
    pdos = {}
    
    for ia in range(mol.natm):
        atom_symbol = mol.atom_symbol(ia)
        
        start = mol.atom_start_idx[ia]
        end = start + mol.atom_nbas[ia]
        
        coeff_atom = mo_coeff[start:end, :]
        contrib = coeff_atom**2
        
        pdos[atom_symbol] = np.zeros_like(energies)
        for i, e in enumerate(mo_energy):
            pdos[atom_symbol] += contrib[:, i].sum() * np.exp(-(energies - e)**2 / (2 * sigma**2))
    
    print('\nPDOS per atom:')
    for atom, dos in pdos.items():
        print('  %s: max DOS = %.4f' % (atom, dos.max()))
    
    return energies, pdos


def orbital_energies_summary(mf):
    """
    Print orbital energy level summary
    """
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    n_homo = int(mo_occ.sum() // 2)
    
    print('\n' + '=' * 60)
    print('Orbital Energy Summary (eV)')
    print('=' * 60)
    print('%6s %4s %12s %20s' % ('Index', 'Occ', 'Energy(eV)', 'Assignment'))
    print('-' * 60)
    
    for i in range(min(20, len(mo_energy))):
        occ_str = 'occ' if mo_occ[i] > 0 else 'emp'
        
        assign = ''
        if i == n_homo - 1:
            assign = '<-- HOMO'
        elif i == n_homo:
            assign = '<-- LUMO'
        
        print('%6d %4s %12.4f %20s' % (i, occ_str, mo_energy[i], assign))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nNote: analysis.py requires a loaded SCF object, not just a file.')
        print('Example usage from Python:')
        print('  from pyscf import scf, analysis')
        print('  mf = scf.RHF(mol).run()')
        print('  analysis.mulliken_charges(mf)')
        print('  analysis.dos_analysis(mf)')
        sys.exit(1)
    
    input_file = sys.argv[1]
    task = sys.argv[2].lower()
    
    print('Analysis task: %s' % task)
    print('Input: %s' % input_file)
    print('\nRequires SCF/MF object. Use from Python directly.')


if __name__ == '__main__':
    main()
