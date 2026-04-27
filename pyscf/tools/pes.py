#!/usr/bin/env python3
"""
PES - Potential Energy Surface scanning

Supports:
    - 1D rigid scan (fixed geometry)
    - 1D relaxed scan (optimize at each point)
    - 2D grid scan

Usage:
    python pes.py <xyz_file> <functional> [basis]
    python pes.py h2_scan.csv b3lyp cc-pvdz --scan-type bond --atom1 1 --atom2 2 --rmin 0.8 --rmax 1.5 --npoints 20
"""

import sys
import numpy as np

from pyscf import gto, dft


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


def rigid_scan_1d(xyz_file, functional='b3lyp', basis='cc-pvdz',
                  atom1=None, atom2=None, coord_type='bond',
                  rmin=0.8, rmax=1.5, npoints=20):
    """
    1D rigid scan (no optimization at each step)
    
    Args:
        atom1, atom2: atom indices (0-based) for bond/distance
        coord_type: 'bond', 'angle', 'dihedral'
        rmin, rmax: range (Angstrom or degree)
        npoints: number of scan points
    """
    symbols, coords = read_xyz(xyz_file)
    
    if atom1 is None:
        atom1 = 0
    if atom2 is None:
        atom2 = 1
    
    if coord_type == 'bond' or coord_type == 'distance':
        values = np.linspace(rmin, rmax, npoints)
        results = []
        
        for r in values:
            coords_new = coords.copy()
            vec = coords_new[atom2] - coords_new[atom1]
            d = np.linalg.norm(vec)
            if d > 0:
                coords_new[atom2] = coords_new[atom1] + vec / d * r
            
            atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                                  for s, c in zip(symbols, coords_new)])
            mol = gto.M(atom=atom_str, basis=basis, verbose=0)
            mf = dft.RKS(mol, xc=functional).run()
            results.append((r, mf.e_tot))
            print('r=%.3f E=%.10f' % (r, mf.e_tot))
        
        print('\n1D Rigid Scan Results:')
        print('%-12s %20s' % ('R (Ang)', 'E (Eh)'))
        print('-' * 35)
        for r, e in results:
            print('%.6f  %.15f' % (r, e))
        
        return np.array(results)
    
    return None


def relaxed_scan_1d(xyz_file, functional='b3lyp', basis='cc-pvdz',
                    scan_range=(0.9, 1.5), npoints=10):
    """
    1D relaxed scan (optimize at each point)
    
    Much slower than rigid scan but gives minimum energy path
    """
    from pyscf import geomopt
    
    symbols, coords = read_xyz(xyz_file)
    r_values = np.linspace(scan_range[0], scan_range[1], npoints)
    energies = []
    
    print('\n1D Relaxed Scan:')
    print('%-12s %20s' % ('R (Ang)', 'E (Eh)'))
    print('-' * 35)
    
    for i, r in enumerate(r_values):
        # Modify geometry (example: scale first bond)
        coords_mod = coords.copy()
        # For a diatomic: scale the bond
        if len(coords) == 2:
            vec = coords[1] - coords[0]
            bond_len = np.linalg.norm(vec)
            coords_mod[1] = coords[0] + vec / bond_len * r
        
        atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                              for s, c in zip(symbols, coords_mod)])
        mol = gto.M(atom=atom_str, basis=basis, verbose=0)
        
        try:
            mf = dft.RKS(mol, xc=functional).run()
            opt = geomopt.GeometryOptimizer(mf)
            mol_opt = opt.run()
            energies.append((r, opt.e_final))
            print('%.6f  %.15f (converged)' % (r, opt.e_final))
        except Exception as e:
            print('%.6f  FAILED: %s' % (r, str(e)[:40]))
            energies.append((r, None))
    
    return energies


def relaxed_scan_2d(xyz_file, functional='b3lyp', basis='cc-pvdz',
                    param1=(0, 1, 0.9, 1.3, 5),
                    param2=(0, 2, 90.0, 120.0, 5)):
    """
    2D relaxed grid scan
    
    param format: (atom_a, atom_b, atom_c, min_val, max_val, npoints)
    
    For dihedral: (atom1, atom2, atom3, atom4, min_deg, max_deg, npoints)
    """
    from pyscf import geomopt
    
    symbols, coords = read_xyz(xyz_file)
    
    ia, ib, v1_min, v1_max, n1 = param1
    ja, jb, jc, v2_min, v2_max, n2 = param2
    
    values1 = np.linspace(v1_min, v1_max, n1)
    values2 = np.linspace(v2_min, v2_max, n2)
    
    print('\n2D Relaxed Scan: %dx%d grid' % (n1, n2))
    print('Parameter 1: atoms %d-%d, range [%.2f, %.2f]' % (ia, ib, v1_min, v1_max))
    print('Parameter 2: atoms %d-%d-%d, range [%.2f, %.2f]' % (ja, jb, jc, v2_min, v2_max))
    
    results = np.full((n1, n2), np.nan)
    
    for i, v1 in enumerate(values1):
        for j, v2 in enumerate(values2):
            coords_mod = apply_internal_coords(coords, ia, ib, ja, jb, jc, v1, v2)
            atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                                  for s, c in zip(symbols, coords_mod)])
            
            try:
                mol = gto.M(atom=atom_str, basis=basis, verbose=0)
                mf = dft.RKS(mol, xc=functional).run()
                opt = geomopt.GeometryOptimizer(mf).run()
                results[i, j] = opt.e_final
            except:
                pass
    
    print('\nEnergy grid (Eh):')
    header = '%8s' % 'v1\\v2'
    for v2 in values2:
        header += '%12.4f' % v2
    print(header)
    for i, v1 in enumerate(values1):
        row = '%8.4f' % v1
        for j in range(n2):
            if np.isnan(results[i, j]):
                row += '%12s' % 'N/A'
            else:
                row += '%12.6f' % results[i, j]
        print(row)
    
    return values1, values2, results


def apply_internal_coords(coords, ia, ib, ja, jb, jc, r_ab, angle_b):
    """Apply bond length and angle constraints"""
    coords_new = coords.copy()
    
    # Fix atom ia at origin
    # Place atom ib at distance r_ab along x-axis
    coords_new[ib] = coords_new[ia] + np.array([r_ab, 0.0, 0.0])
    
    # Place atom jc at angle from ib-ja axis
    angle_rad = np.deg2rad(angle_b)
    coords_new[jc] = coords_new[ib] + np.array([
        np.cos(angle_rad), np.sin(angle_rad), 0.0
    ])
    
    return coords_new


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nExamples:')
        print('  python pes.py h2.csv b3lyp cc-pvdz --rigid --a1 0 --a2 1 --rmin 0.6 --rmax 2.0 -n 20')
        print('  python pes.py h2o.csv pbe 6-31G* --relaxed')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    scan_type = 'rigid'
    if '--relaxed' in sys.argv:
        scan_type = 'relaxed'
    
    if scan_type == 'rigid':
        rigid_scan_1d(xyz_file, functional, basis,
                      atom1=0, atom2=1,
                      rmin=0.8, rmax=2.0, npoints=20)
    else:
        relaxed_scan_1d(xyz_file, functional, basis,
                       scan_range=(0.8, 2.0), npoints=10)


if __name__ == '__main__':
    main()
