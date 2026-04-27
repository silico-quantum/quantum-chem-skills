#!/usr/bin/env python3
"""
DFT - Density Functional Theory functional selection and comparison

Functional categories:
    LDA: SVWN5
    GGA: PBE, BLYP, BP86
    meta-GGA: TPSS, SCAN
    hybrid: B3LYP, PBE0, M06, M06-2X
    range-separated: WB97X, WB97X-D, CAM-B3LYP
    dispersion: WB97X-D3BJ

Usage:
    python dft.py <xyz_file> <functional> [basis]
    python dft.py water.csv pbe cc-pvdz
    python dft.py mol.csv wb97x-d3bj def2-tzvp
    python dft.py compare.csv "b3lyp,pbe0,wb97x-d3bj" 6-31G*
"""

import sys
import numpy as np

from pyscf import gto, dft


FUNCTIONAL_CATEGORIES = {
    'lda': ['svwn5', 'lda'],
    'gga': ['pbe', 'blyp', 'bp86', 'pbesol'],
    'meta_gga': ['tpss', 'scan', 'm06-l', 'tpss0'],
    'hybrid': ['b3lyp', 'pbe0', 'm06', 'm06-2x', 'b3pw91', 'wb97x', 'scan0'],
    'range_separated': ['wb97x', 'wb97x-d', 'wb97x-d3bj', 'camb3lyp', 'lc-wpbe'],
    'dispersion': ['pbe-d3', 'b3lyp-d3', 'wb97x-d3bj'],
}

RECOMMENDED = {
    'Organic molecules (balanced)': 'b3lyp',
    'Organic molecules (fast)': 'pbe',
    'Charge transfer / CT states': 'wb97x-d3bj',
    'Weak interactions / vdW': 'wb97x-d3bj',
    'Excited states (TDDFT)': 'b3lyp',
    'Solvation': 'pbe0',
    'Transition metals': 'b3lyp-d3',
    'Heavy atoms (relativistic)': 'wpbesol',
}


def list_functionals(category=None):
    """List all available functionals"""
    if category:
        return FUNCTIONAL_CATEGORIES.get(category, [])
    else:
        all_funcs = set()
        for funcs in FUNCTIONAL_CATEGORIES.values():
            all_funcs.update(funcs)
        return sorted(all_funcs)


def recommend_functional(use_case):
    """Get recommended functional for a use case"""
    return RECOMMENDED.get(use_case, 'b3lyp')


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


def run_dft(xyz_file, functional='pbe', basis='cc-pvdz',
            restricted=True, spin=0, charge=0):
    """
    Run DFT calculation
    
    Args:
        xyz_file: XYZ file
        functional: functional name
        basis: basis set
        restricted: True = RKS, False = UKS
        spin: spin multiplicity
        charge: molecular charge
    
    Returns:
        dict: {'mf': SCF object, 'mol': molecule object}
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=4)
    
    if restricted:
        mf = dft.RKS(mol, xc=functional).run()
    else:
        mf = dft.UKS(mol, xc=functional).run()
    
    print('\n' + '=' * 50)
    print('DFT Results: %s/%s' % (functional, basis))
    print('=' * 50)
    print('Total energy: %.10f Eh' % mf.e_tot)
    print('Total electrons: %d' % mol.nelectron)
    
    mo_e = mf.mo_energy * 27.2114
    occ = mf.mo_occ
    
    homo_idx = np.where(occ > 0)[0].max()
    lumo_idx = np.where(occ == 0)[0].min()
    
    print('\nHOMO: %.4f eV' % mo_e[homo_idx])
    print('LUMO: %.4f eV' % mo_e[lumo_idx])
    print('Gap:  %.4f eV' % (mo_e[lumo_idx] - mo_e[homo_idx]))
    
    return {'mf': mf, 'mol': mol, 'symbols': symbols}


def compare_functionals(xyz_file, functionals, basis='cc-pvdz'):
    """
    Compare energy and HOMO-LUMO gap across functionals
    
    Args:
        functionals: comma-separated list or list of functional names
    """
    if isinstance(functionals, str):
        functionals = [f.strip() for f in functionals.split(',')]
    
    print('\nFunctional Comparison:')
    print('%-20s %15s %10s %10s %10s' % ('Functional', 'E(tot)', 'HOMO', 'LUMO', 'Gap'))
    print('-' * 70)
    
    results = []
    for xc in functionals:
        try:
            result = run_dft(xyz_file, xc, basis, output_file=None)
            mf = result['mf']
            mo_e = mf.mo_energy * 27.2114
            occ = mf.mo_occ
            
            homo = mo_e[occ > 0].max()
            lumo = mo_e[occ == 0].min()
            gap = lumo - homo
            
            print('%-20s %15.8f %10.4f %10.4f %10.4f' % (xc, mf.e_tot, homo, lumo, gap))
            results.append({'xc': xc, 'e_tot': mf.e_tot, 'homo': homo,
                          'lumo': lumo, 'gap': gap})
        except Exception as e:
            print('%-20s ERROR: %s' % (xc, str(e)[:40]))
    
    return results


def grid_convergence_test(xyz_file, functional='pbe',
                         grids=[(50, 110), (75, 195), (100, 290)]):
    """
    Test SCF convergence vs integration grid size
    
    grids: list of (atomic_radii, lebedev_n) tuples
    """
    print('\nGrid Convergence Test:')
    
    for rad, leb in grids:
        try:
            mol = gto.M(atom=open(xyz_file).read(), basis='cc-pvdz')
            mf = dft.RKS(mol, xc=functional)
            mf.grids.atom_radius = rad
            mf.grids.lebedev_grids = {'': (leb, 194)}
            mf.run()
            
            print('  Grid (%d, %d): E = %.10f' % (rad, leb, mf.e_tot))
        except Exception as e:
            print('  Grid (%d, %d): FAILED - %s' % (rad, leb, str(e)[:40]))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nRecommended functionals:')
        for use, func in RECOMMENDED.items():
            print('  %s: %s' % (use, func))
        print('\nFunctional categories:')
        for cat in FUNCTIONAL_CATEGORIES:
            funcs = FUNCTIONAL_CATEGORIES[cat]
            print('  %s: %s' % (cat, ', '.join(funcs[:5])))
        print('\nExamples:')
        print('  python dft.py water.csv b3lyp cc-pvdz')
        print('  python dft.py mol.csv wb97x-d3bj def2-tzvp')
        print('  python dft.py compare.csv "b3lyp,pbe0,m06-2x" 6-31G*')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    if ',' in functional:
        compare_functionals(xyz_file, functional, basis)
    else:
        run_dft(xyz_file, functional, basis)


if __name__ == '__main__':
    main()
