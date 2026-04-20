#!/usr/bin/env python3
"""
MP2 - Second-Order Møller-Plesset Perturbation Theory

Supports:
    - RMP2 (restricted closed-shell)
    - UMP2 (unrestricted open-shell)
    - RI-MP2 (density fitting)

Usage:
    python mp2.py <xyz_file> [basis] [method]
    python mp2.py water.csv cc-pvdz rmp2
    python mp2.py o2_triplet.csv sto-3g ump2
"""

import sys
import numpy as np

from pyscf import gto, scf, mp


def read_xyz(xyz_file):
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    return symbols, np.array(coords)


def run_rmp2(xyz_file, basis='cc-pvdz', density_fit=False, auxbasis=None):
    """
    Closed-shell RMP2 calculation
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    mf = scf.RHF(mol).run()
    
    if density_fit:
        mmp = mp.MP2(mf).density_fit(auxbasis=auxbasis).run()
    else:
        mmp = mp.MP2(mf).run()
    
    print('\n' + '=' * 50)
    print('RMP2 Results: %s' % basis)
    print('=' * 50)
    print('E(SCF)   = %.10f Eh' % mf.e_tot)
    print('E(MP2)   = %.10f Eh' % mmp.e_tot)
    print('E(corr)  = %.10f Eh' % mmp.e_corr)
    print('T2 shape:', mmp.t2.shape)
    
    return {'mp': mmp, 'mf': mf, 'mol': mol}


def run_ump2(xyz_file, basis='cc-pvdz'):
    """
    Open-shell UMP2 calculation
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, spin=1, verbose=4)
    mf = scf.UHF(mol).run()
    
    mmp = mp.UMP2(mf).run()
    
    print('\n' + '=' * 50)
    print('UMP2 Results: %s' % basis)
    print('=' * 50)
    print('E(SCF)   = %.10f Eh' % mf.e_tot)
    print('E(MP2)   = %.10f Eh' % mmp.e_tot)
    print('E(corr)  = %.10f Eh' % mmp.e_corr)
    
    return {'mp': mmp, 'mf': mf, 'mol': mol}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print('\nExamples:')
        print('  python mp2.py water.csv cc-pvdz')
        print('  python mp2.py n2.csv aug-cc-pvdz rmp2')
        print('  python mp2.py radical.csv sto-3g ump2')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    basis = sys.argv[2] if len(sys.argv) > 2 else 'cc-pvdz'
    method = sys.argv[3].lower() if len(sys.argv) > 3 else 'rmp2'
    
    if method == 'ump2':
        run_ump2(xyz_file, basis)
    else:
        run_rmp2(xyz_file, basis)


if __name__ == '__main__':
    main()
