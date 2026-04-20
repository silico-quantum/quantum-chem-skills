#!/usr/bin/env python3
"""
TDDFT - Time-Dependent DFT excited state calculations

Supports:
    - TDDFT (linear response)
    - TDA (Tamm-Dancoff approximation)
    - CIS (Configuration Interaction Singles)
    - UV-Vis spectrum generation

Usage:
    python tddft.py <xyz_file> <functional> [basis] [nstates]
    python tddft.py water.csv b3lyp cc-pvdz 10
"""

import sys
import numpy as np

from pyscf import gto, dft, tddft


def read_xyz(xyz_file):
    """Read XYZ file"""
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    return symbols, np.array(coords)


def run_tddft(xyz_file, functional='b3lyp', basis='cc-pvdz',
               nstates=10, restricted=True, charge=0, spin=0):
    """
    Run TDDFT calculation
    
    Args:
        xyz_file: XYZ file
        functional: DFT functional
        basis: basis set
        nstates: number of excited states
        restricted: RKS (True) or UKS (False)
        charge: molecular charge
        spin: spin multiplicity
    
    Returns:
        dict with 'td', 'mf', 'mol'
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=4)
    
    if restricted:
        mf = dft.RKS(mol, xc=functional).run()
    else:
        mf = dft.UKS(mol, xc=functional).run()
    
    td = tddft.TDDFT(mf).run(nstates=nstates)
    
    e = td.e * 27.2114  # eV
    
    print('\n' + '=' * 60)
    print('TDDFT Results: %s/%s' % (functional, basis))
    print('=' * 60)
    print('%-6s %12s %10s %10s' % ('State', 'Energy(eV)', 'f', 'Character'))
    print('-' * 60)
    
    for i in range(min(nstates, len(e))):
        occ = mf.mo_occ
        # Identify dominant excitations
        print('%-6d %12.4f %10.4f' % (i+1, e[i], td.f[i]))
    
    # Analyze dominant excitations
    print('\nExcitation analysis:')
    td.analyze()
    
    return {'td': td, 'mf': mf, 'mol': mol}


def run_tda(xyz_file, functional='b3lyp', basis='cc-pvdz',
            nstates=10, restricted=True):
    """Run TDA calculation"""
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    
    if restricted:
        mf = dft.RKS(mol, xc=functional).run()
    else:
        mf = dft.UKS(mol, xc=functional).run()
    
    td = tddft.TDA(mf).run(nstates=nstates)
    
    e = td.e * 27.2114
    print('\nTDA Results:')
    print('%-6s %12s %10s' % ('State', 'Energy(eV)', 'f'))
    print('-' * 35)
    for i in range(min(nstates, len(e))):
        print('%-6d %12.4f %10.4f' % (i+1, e[i], td.f[i]))
    
    return {'td': td, 'mf': mf, 'mol': mol}


def generate_spectrum(td, sigma=0.1, emin=0, emax=15, npoints=1000):
    """
    Generate broadened UV-Vis spectrum
    
    Returns:
        tuple: (energies, spectrum)
    """
    e = td.e * 27.2114
    f = td.f
    
    energies = np.linspace(emin, emax, npoints)
    spectrum = np.zeros_like(energies)
    
    for ei, fi in zip(e, f):
        if fi > 0:
            spectrum += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))
    
    return energies, spectrum


def compare_functionals(xyz_file, functionals, basis='cc-pvdz', nstates=5):
    """Compare TDDFT results across functionals"""
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=0)
    mf = dft.RKS(mol, xc='b3lyp').run()
    
    print('\n%-10s %8s %10s %10s' % ('Functional', 'S1(eV)', 'f(S1)', 'T1(eV)'))
    print('-' * 45)
    
    for xc in functionals:
        try:
            mf.xc = xc
            mf.run()
            td = tddft.TDDFT(mf).run(nstates=nstates)
            e = td.e * 27.2114
            f = td.f
            
            s1 = e[0] if len(e) > 0 else 0
            t1 = e[nstates] if len(e) > nstates else 0
            f1 = f[0] if len(f) > 0 else 0
            
            print('%-10s %8.3f %10.4f %10.3f' % (xc, s1, f1, t1))
        except Exception as ex:
            print('%-10s ERROR: %s' % (xc, str(ex)[:30]))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nExamples:')
        print('  python tddft.py water.csv b3lyp cc-pvdz 10')
        print('  python tddft.py mol.csv wb97x-d3bj def2-tzvp 20')
        print('  python tddft.py test.xyz pbe 6-31G* 5 --tda')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    nstates = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    
    use_tda = '--tda' in sys.argv
    
    if use_tda:
        run_tda(xyz_file, functional, basis, nstates)
    else:
        run_tddft(xyz_file, functional, basis, nstates)


if __name__ == '__main__':
    main()
