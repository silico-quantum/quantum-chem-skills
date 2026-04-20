#!/usr/bin/env python3
"""
CCSD - Coupled-Cluster Singles and Doubles

Supports:
    - RCCSD (restricted closed-shell)
    - UCCSD (unrestricted open-shell)
    - CCSD(T) perturbative triples

Usage:
    python ccsd.py <xyz_file> [basis]
    python ccsd.py n2.csv cc-pvdz
"""

import sys
import numpy as np

from pyscf import gto, scf, cc


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


def run_ccsd(xyz_file, basis='cc-pvdz', with_t=True, charge=0, spin=0):
    """
    Run CCSD calculation
    
    Args:
        xyz_file: XYZ file
        basis: basis set
        with_t: if True, compute CCSD(T) correction
        charge: molecular charge
        spin: spin multiplicity (0=closed-shell singlet)
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=4)
    
    if spin == 0:
        mf = scf.RHF(mol).run()
        mycc = cc.CCSD(mf).run()
    else:
        mf = scf.UHF(mol).run()
        mycc = cc.UCCSD(mf).run()
    
    print('\n' + '=' * 50)
    print('CCSD Results: %s' % basis)
    print('=' * 50)
    print('E(SCF)   = %.10f Eh' % mf.e_tot)
    print('E(CCSD)  = %.10f Eh' % mycc.e_tot)
    
    if hasattr(mycc, 't1'):
        print('T1 shape:', mycc.t1.shape)
    if hasattr(mycc, 't2'):
        print('T2 shape:', mycc.t2.shape)
    
    # CCSD(T) correction
    if with_t and spin == 0:
        et = mycc.ccsd_t()
        e_total = mycc.e_tot + et
        print('\nE(CCSD(T)) = %.10f Eh' % e_total)
        print('ΔE(T)     = %.10f Eh' % et)
    
    # One-particle density matrix
    dm = mycc.make_rdm1()
    
    return {'cc': mycc, 'mf': mf, 'mol': mol}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print('\nExamples:')
        print('  python ccsd.py n2.csv cc-pvdz')
        print('  python ccsd.py h2o.csv 6-31G* --no-t')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    basis = sys.argv[2] if len(sys.argv) > 2 else 'cc-pvdz'
    with_t = '--no-t' not in sys.argv
    
    charge, spin = 0, 0
    if '--charge' in sys.argv:
        idx = sys.argv.index('--charge')
        charge = int(sys.argv[idx + 1])
    if '--spin' in sys.argv:
        idx = sys.argv.index('--spin')
        spin = int(sys.argv[idx + 1])
    
    run_ccsd(xyz_file, basis, with_t, charge, spin)


if __name__ == '__main__':
    main()
