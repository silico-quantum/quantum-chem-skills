#!/usr/bin/env python3
"""
Geometry Optimization - Molecular geometry optimization and transition state search

Supports:
    - Energy minimization (BFGS, CG, Newton-Raphson)
    - Transition state (dimer method)
    - Frequency calculation (verify minima/TS)
    - Constrained optimization

Usage:
    python geometry.py <xyz_file> <functional> [basis]
    python geometry.py water.csv b3lyp cc-pvdz
    python geometry.py ts_guess.xyz pbe 6-31G* --ts
"""

import sys
import numpy as np

from pyscf import gto, dft, geomopt


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


def run_optimization(xyz_file, functional='b3lyp', basis='cc-pvdz',
                     charge=0, spin=0, method='bfgs'):
    """
    Geometry optimization
    
    Args:
        xyz_file: XYZ file
        functional: DFT functional
        basis: basis set
        charge: molecular charge
        spin: spin multiplicity
        method: 'bfgs', 'cg', 'newton'
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=4)
    
    mf = dft.RKS(mol, xc=functional).run()
    
    # Choose optimizer
    if method == 'newton':
        opt = geomopt.NewtonOptimizer(mf)
    elif method == 'cg':
        opt = geomopt.CGOptimizer(mf)
    else:
        opt = geomopt.GeometryOptimizer(mf)
    
    mol_opt = opt.run()
    
    print('\n' + '=' * 50)
    print('Optimization Results: %s/%s' % (functional, basis))
    print('=' * 50)
    print('Initial E = %.10f Eh' % mf.e_tot)
    print('Final E   = %.10f Eh' % opt.e_final)
    print('\nOptimized geometry (Angstrom):')
    for i in range(mol_opt.natm):
        coord = mol_opt.atom_coord(i) * 0.529177
        print('  %2s %12.6f %12.6f %12.6f' % (
            mol_opt.atom_symbol(i), coord[0], coord[1], coord[2]))
    
    return {'mol_opt': mol_opt, 'mf': mf}


def run_ts_search(xyz_file, functional='b3lyp', basis='cc-pvdz'):
    """
    Transition state search using dimer method
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    mf = dft.RKS(mol, xc=functional).run()
    
    opt = geomopt.DimerOptimizer(mf).run()
    
    print('\n' + '=' * 50)
    print('Transition State Search')
    print('=' * 50)
    print('E(TS) = %.10f Eh' % opt.e_final)
    print('\nTS geometry (Angstrom):')
    for i in range(mol.natm):
        coord = mol.atom_coord(i) * 0.529177
        print('  %2s %12.6f %12.6f %12.6f' % (
            mol.atom_symbol(i), coord[0], coord[1], coord[2]))
    
    return {'mol': mol, 'mf': mf}


def run_frequency(mf_or_mol, functional='b3lyp', basis='cc-pvdz'):
    """
    Frequency calculation to verify stationary point
    
    Returns:
        list: vibrational frequencies (cm^-1)
    """
    from pyscf import hessian
    
    if hasattr(mf_or_mol, 'mol'):
        mol = mf_or_mol.mol
        mf = mf_or_mol
    else:
        mol = mf_or_mol
        mf = dft.RKS(mol, xc=functional).run()
    
    h = hessian.RHF(mf).run()
    freq = h.hessian()
    
    # Convert to vibrational frequencies
    from pyscf.hessian import thermo
    freq_cm = thermo.vibration(mol, freq)
    
    print('\n' + '=' * 50)
    print('Frequency Analysis')
    print('=' * 50)
    print('Frequencies (cm^-1):')
    for i, f in enumerate(freq_cm):
        tag = ''
        if f < 0:
            tag = ' <-- imaginary (TS)'
        elif i < 6:
            tag = ' <-- translation/rotation'
        print('  %3d %10.4f%s' % (i+1, f, tag))
    
    # Thermodynamic analysis
    print('\nThermodynamic properties (298.15 K, 1 atm):')
    thermo_info = thermo.thermo(mf, freq_cm, mol.natm)
    print('  ZPE: %.4f Eh' % thermo_info['ZPE'])
    print('  H:   %.4f Eh' % thermo_info['H'])
    print('  S:   %.4f Eh/K' % thermo_info['S'])
    
    return freq_cm


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nExamples:')
        print('  python geometry.py water.csv b3lyp cc-pvdz')
        print('  python geometry.py ts_guess.xyz pbe 6-31G* --ts')
        print('  python geometry.py optmol.xyz b3lyp def2-tzvp --freq')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    ts_mode = '--ts' in sys.argv
    freq_mode = '--freq' in sys.argv
    
    if ts_mode:
        result = run_ts_search(xyz_file, functional, basis)
    else:
        result = run_optimization(xyz_file, functional, basis)
    
    if freq_mode:
        run_frequency(result['mf'])


if __name__ == '__main__':
    main()
