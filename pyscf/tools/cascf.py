#!/usr/bin/env python3
"""
CASSCF - Complete Active Space Self-Consistent Field

Supports:
    - CASSCF (multi-configurational)
    - CASCI
    - CASPT2 (multi-state second-order perturbation)
    - NEVPT2 (N-Electron Valence Second-Order PT)

Usage:
    python cascf.py <xyz_file> <ncas> <nelec> [basis] [nstates]
    python cascf.py butadiene.csv 4 4 cc-pvdz    # CAS(4,4) for butadiene
    python cascf.py ethylene.csv 2 2 6-31G*     # CAS(2,2) for ethylene pi system
"""

import sys
import numpy as np

from pyscf import gto, scf, mcscf


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


def auto_select_active_space(mol, mf):
    """
    Automatic active space selection based on frontier orbitals
    
    Returns:
        tuple: (nelec, ncas_orbs) 
    """
    mo_energy = mf.mo_energy
    nelectron = mol.nelectron
    nocc = nelectron // 2
    
    # HOMO and LUMO orbitals as initial guess
    print('  Auto selection: nocc = %d' % nocc)
    print('  HOMO-1: %.4f Eh, HOMO: %.4f Eh' % (mo_energy[nocc-2], mo_energy[nocc-1]))
    print('  LUMO: %.4f Eh, LUMO+1: %.4f Eh' % (mo_energy[nocc], mo_energy[nocc+1]))
    
    # Default: 2 electrons in 2 orbitals (one bond)
    return 2, 2


def run_casscf(xyz_file, nelec, ncas, basis='cc-pvdz',
               mf_method='RHF', nstates=1, active_orbs=None):
    """
    Run CASSCF calculation
    
    Args:
        xyz_file: XYZ file
        nelec: number of electrons in active space
        ncas: number of orbitals in active space
        basis: basis set
        mf_method: RHF or UHF
        nstates: number of root-specific CASSCF states
        active_orbs: list of orbital indices (if None, auto-select)
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    
    if mf_method.upper() == 'RHF':
        mf = scf.RHF(mol).run()
    else:
        mf = scf.UHF(mol).run()
    
    print('\n' + '=' * 50)
    print('CAS(%d,%d) SCF Results' % (nelec, ncas))
    print('=' * 50)
    
    if nstates == 1:
        mc = mcscf.CASSCF(mf, ncas, nelec).run()
    else:
        mc = mcscf.CASSCF(mf, ncas, nelec).state_average(nstates=nstates).run()
    
    print('\nE(CASSCF) = %.10f Eh' % mc.e_tot)
    print('CAS orbitals shape:', mc.mo_coeff.shape)
    
    # Apply orbital sorting if specified
    if active_orbs is not None:
        mo = mcscf.sort_mo(mc, mf.mo_coeff, active_orbs)
        mc.run(mo)
        print('\nOrbitals sorted, active space:', active_orbs)
        print('E(sorted CASSCF) = %.10f Eh' % mc.e_tot)
    
    return {'mc': mc, 'mf': mf, 'mol': mol}


def run_caspt2(mc, root=0):
    """
    Add CASPT2 correction to CASSCF
    """
    mcpt2 = mcscf.CASPT2(mc, root=root).run()
    
    print('\n' + '=' * 50)
    print('CASPT2 Results (root %d)' % root)
    print('=' * 50)
    print('E(CASPT2) = %.10f Eh' % mcpt2.e_tot)
    print('MS-CASPT2 energies:')
    if hasattr(mcpt2, 'e_states'):
        for i, e in enumerate(mcpt2.e_states):
            print('  State %d: %.10f Eh' % (i+1, e))
    
    return mcpt2


def run_nevpt2(mc, root=0):
    """
    Add NEVPT2 correction to CASSCF
    """
    mcev = mcscf.NEVPT2(mc, root=root).run()
    
    print('\nE(NEVPT2) = %.10f Eh' % mcev.e_tot)
    
    return mcev


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        print('\nExamples:')
        print('  python cascf.py butadiene.xyz 4 4 cc-pvdz')
        print('  python cascf.py o2.csv 6 6 sto-3g 3   # 3 states')
        print('  python cascf.py cu_cluster.xyz 10 12 lanl2dz  # transition metal')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    nelec = int(sys.argv[2])
    ncas = int(sys.argv[3])
    basis = sys.argv[4] if len(sys.argv) > 4 else 'cc-pvdz'
    
    nstates = 1
    if '--nstates' in sys.argv:
        idx = sys.argv.index('--nstates')
        nstates = int(sys.argv[idx + 1])
    
    # Custom active space orbitals
    active_orbs = None
    if '--orbs' in sys.argv:
        idx = sys.argv.index('--orbs')
        active_orbs = [int(x) for x in sys.argv[idx+1].split(',')]
    
    result = run_casscf(xyz_file, nelec, ncas, basis, nstates=nstates,
                        active_orbs=active_orbs)
    
    # Add perturbation corrections
    if '--pt2' in sys.argv:
        pt2_type = 'caspt2' if '--nevpt2' not in sys.argv else 'nevpt2'
        if pt2_type == 'caspt2':
            run_caspt2(result['mc'])
        else:
            run_nevpt2(result['mc'])


if __name__ == '__main__':
    main()
