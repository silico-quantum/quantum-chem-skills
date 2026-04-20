#!/usr/bin/env python3
"""
SCF - Self-Consistent Field calculations
Hartree-Fock (RHF, UHF, ROHF) framework

Usage:
    python scf.py <xyz_file> <method> [basis]
    python scf.py water.xyz RHF cc-pvdz
    python scf.py molecule.xyz UHF sto-3g
"""

import sys
import numpy as np

from pyscf import gto, scf


def read_xyz(xyz_file):
    """Read XYZ file and return (symbols, coords)"""
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    return symbols, np.array(coords)


def build_mol(symbols, coords, basis='cc-pvdz', charge=0, spin=0):
    """Build PySCF molecule object"""
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                         for s, c in zip(symbols, coords)])
    mol = gto.M(
        atom=atom_str,
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=4
    )
    return mol


def run_scf(xyz_file, method='RHF', basis='cc-pvdz', 
            charge=0, spin=0, chkfile=None):
    """
    Run SCF calculation
    
    Args:
        xyz_file: XYZ file path
        method: RHF, UHF, ROHF
        basis: basis set
        charge: molecular charge
        spin: spin multiplicity (0=singlet, 1=doublet, 2=triplet)
        chkfile: checkpoint file (optional)
    
    Returns:
        dict with 'mf', 'mol', 'symbols'
    """
    symbols, coords = read_xyz(xyz_file)
    mol = build_mol(symbols, coords, basis, charge, spin)
    
    method = method.upper()
    if method == 'RHF':
        mf = scf.RHF(mol)
    elif method == 'UHF':
        mf = scf.UHF(mol)
    elif method == 'ROHF':
        mf = scf.ROHF(mol)
    else:
        raise ValueError('Unknown method: %s' % method)
    
    if chkfile:
        mf.chkfile = chkfile
    
    mf.run()
    
    # Summary
    n electrons = int(mol.nelectron)
    n_bf = mol.nao
    mo_e = mf.mo_energy * 27.2114  # eV
    
    print('\n' + '=' * 50)
    print('SCF Results: %s/%s' % (method, basis))
    print('=' * 50)
    print('Total energy: %.10f Eh' % mf.e_tot)
    print('Nuclear repulsion: %.10f Eh' % mol.energy_nuc())
    print('Electrons: %d' % n_electrons)
    print('Basis functions: %d' % n_bf)
    
    occ = mf.mo_occ
    homo_idx = np.where(occ > 0)[0].max()
    lumo_idx = np.where(occ == 0)[0].min()
    
    print('\nOrbital energies (eV):')
    print('  HOMO (%d): %.4f eV' % (homo_idx, mo_e[homo_idx]))
    print('  LUMO (%d): %.4f eV' % (lumo_idx, mo_e[lumo_idx]))
    print('  Gap: %.4f eV' % (mo_e[lumo_idx] - mo_e[homo_idx]))
    
    return {'mf': mf, 'mol': mol, 'symbols': symbols}


def compare_basis(xyz_file, basis_sets, method='RHF'):
    """
    Compare energies across basis sets
    
    Args:
        basis_sets: list of basis set names
    """
    symbols, coords = read_xyz(xyz_file)
    
    print('\nBasis set comparison:')
    print('%-15s %20s %10s' % ('Basis', 'E(tot)', 'Gap(eV)'))
    print('-' * 50)
    
    for basis in basis_sets:
        try:
            mol = build_mol(symbols, coords, basis)
            if method == 'RHF':
                mf = scf.RHF(mol).run()
            else:
                mf = scf.UHF(mol).run()
            
            mo_e = mf.mo_energy * 27.2114
            occ = mf.mo_occ
            gap = mo_e[occ == 0].min() - mo_e[occ > 0].max()
            
            print('%-15s %20.10f %10.4f' % (basis, mf.e_tot, gap))
        except Exception as e:
            print('%-15s ERROR: %s' % (basis, str(e)[:40]))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nExamples:')
        print('  python scf.py water.xyz RHF cc-pvdz')
        print('  python scf.py mol.xyz UHF sto-3g --charge -1')
        print('  python scf.py test.xyz RHF 6-31G* --spin 1')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    method = sys.argv[2].upper()
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    charge, spin = 0, 0
    if '--charge' in sys.argv:
        idx = sys.argv.index('--charge')
        charge = int(sys.argv[idx + 1])
    if '--spin' in sys.argv:
        idx = sys.argv.index('--spin')
        spin = int(sys.argv[idx + 1])
    
    run_scf(xyz_file, method, basis, charge, spin)


if __name__ == '__main__':
    main()
