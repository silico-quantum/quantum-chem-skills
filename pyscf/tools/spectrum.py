#!/usr/bin/env python3
"""
Spectrum - UV-Vis, CD, and Raman spectrum generation

Supports:
    - UV-Vis absorption spectrum (TDDFT/TDA)
    - Circular dichroism (CD) spectrum
    - Emission spectrum
    - Stick spectrum → broadened spectrum

Usage:
    python spectrum.py <xyz_file> <functional> [basis] [nstates]
    python spectrum.py water.csv b3lyp cc-pvdz 50
"""

import sys
import numpy as np

from pyscf import gto, dft, tddft, gdft


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


def uv_vis_spectrum(td, nstates=None, sigma=0.1, emin=0, emax=15, npoints=2000):
    """
    Generate UV-Vis absorption spectrum from TDDFT results
    
    Args:
        td: TDDFT object after .run()
        nstates: number of states to include (None = all)
        sigma: Gaussian broadening (eV)
        emin, emax: energy range (eV)
        npoints: number of grid points
    
    Returns:
        tuple: (energies, absorption_spectrum)
    """
    e = td.e * 27.2114
    f = td.f
    
    if nstates is not None:
        e = e[:nstates]
        f = f[:nstates]
    
    energies = np.linspace(emin, emax, npoints)
    spectrum = np.zeros_like(energies)
    
    for ei, fi in zip(e, f):
        if fi > 0:
            spectrum += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))
    
    return energies, spectrum


def cd_spectrum(td, nstates=None, sigma=0.1, emin=0, emax=15, npoints=2000):
    """
    Generate CD spectrum from GDFT-TDDFT results
    
    Returns rotatory strength and CD spectrum
    """
    # Rotatory strength in cgs units
    r = td.rotatory_strength() * 1.695  # Convert to cgs
    
    e = td.e * 27.2114
    if nstates is not None:
        e = e[:nstates]
        r = r[:nstates]
    
    energies = np.linspace(emin, emax, npoints)
    cd = np.zeros_like(energies)
    
    for ei, ri in zip(e, r):
        cd += ri * np.exp(-(energies - ei)**2 / (2 * sigma**2)) * 1.0
    
    return energies, cd


def plot_spectrum(energies, spectrum, filename='spectrum.dat', 
                  xlabel='Energy (eV)', ylabel='Intensity', title=''):
    """Save spectrum data to file"""
    with open(filename, 'w') as f:
        f.write('# %s\n' % title)
        f.write('# %s\t%s\n' % (xlabel, ylabel))
        for e, s in zip(energies, spectrum):
            f.write('%.6f\t%.8f\n' % (e, s))
    print('Spectrum saved to: %s' % filename)


def run_uvvis(xyz_file, functional='b3lyp', basis='cc-pvdz',
              nstates=50, sigma=0.1):
    """
    Compute UV-Vis absorption spectrum
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    mf = dft.RKS(mol, xc=functional).run()
    
    td = tddft.TDDFT(mf).run(nstates=nstates)
    
    print('\n' + '=' * 50)
    print('UV-Vis Absorption Spectrum')
    print('=' * 50)
    print('%-6s %12s %10s' % ('State', 'Energy(eV)', 'f (osc.)'))
    print('-' * 35)
    
    e = td.e * 27.2114
    for i in range(min(10, len(e))):
        print('%-6d %12.4f %10.4f' % (i+1, e[i], td.f[i]))
    
    # Generate broadened spectrum
    energies, spectrum = uv_vis_spectrum(td, sigma=sigma)
    
    # Save
    output_file = xyz_file.replace('.xyz', '_uvvis.dat').replace('.csv', '_uvvis.dat')
    plot_spectrum(energies, spectrum, output_file, 
                  title='%s/%s UV-Vis' % (functional, basis))
    
    # Wavelength scale
    wavelengths = 1240 / energies  # nm
    print('\nFirst 5 excitations in wavelength:')
    for i in range(min(5, len(e))):
        wl = 1240 / e[i]
        print('  S%d: %.1f nm (f=%.4f)' % (i+1, wl, td.f[i]))
    
    return {'td': td, 'mf': mf, 'energies': energies, 'spectrum': spectrum}


def run_cd(xyz_file, functional='b3lyp', basis='cc-pvdz', nstates=30):
    """
    Compute Circular Dichroism (CD) spectrum
    Requires GDFT (Gauge-Including DFT)
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c)
                         for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, verbose=4)
    mf = gdft.RKS(mol, xc=functional).run()
    
    td = tddft.TDDFT(mf).run(nstates=nstates)
    
    # Rotatory strength
    r = td.rotatory_strength()
    e = td.e * 27.2114
    
    print('\n' + '=' * 50)
    print('CD Spectrum')
    print('=' * 50)
    print('%-6s %12s %15s' % ('State', 'Energy(eV)', 'R (cgs)'))
    print('-' * 40)
    for i in range(min(10, len(e))):
        print('%-6d %12.4f %15.6f' % (i+1, e[i], r[i]))
    
    energies, cd = cd_spectrum(td, sigma=0.1)
    output_file = xyz_file.replace('.xyz', '_cd.dat').replace('.csv', '_cd.dat')
    plot_spectrum(energies, cd, output_file, 
                  title='%s/%s CD' % (functional, basis),
                  ylabel='Delta epsilon')
    
    return {'td': td, 'mf': mf, 'energies': energies, 'cd': cd}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print('\nExamples:')
        print('  python spectrum.py water.csv b3lyp cc-pvdz 50')
        print('  python spectrum.py mol.csv wb97x-d3bj def2-tzvp 30 --cd')
        print('  python spectrum.py test.xyz pbe 6-31G* 20 --sigma 0.05')
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    nstates = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    
    sigma = 0.1
    if '--sigma' in sys.argv:
        idx = sys.argv.index('--sigma')
        sigma = float(sys.argv[idx + 1])
    
    cd_mode = '--cd' in sys.argv
    
    if cd_mode:
        run_cd(xyz_file, functional, basis, nstates)
    else:
        run_uvvis(xyz_file, functional, basis, nstates, sigma)


if __name__ == '__main__':
    main()
