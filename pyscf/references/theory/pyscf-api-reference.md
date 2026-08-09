# Detailed PySCF API Reference

> **Unvalidated historical draft.** Some names and signatures in this document
> may not exist in the installed PySCF version. Use the official versioned API
> documentation and verify each call before execution. See
> [`../../tools/README.md`](../../tools/README.md).

## Module Imports

```python
from pyscf import gto      # Basis-set and molecular definitions
from pyscf import scf      # SCF methods
from pyscf import dft      # DFT methods
from pyscf import mp       # Perturbation theory (MP2, etc.)
from pyscf import cc       # Coupled cluster (CCSD, etc.)
from pyscf import mcscf    # Multiconfigurational SCF
from pyscf import fci      # Full CI
from pyscf import tdscf    # Time-dependent SCF (TDDFT)
from pyscf import ao2mo    # AO-to-MO integral transformation
from pyscf import df       # Density fitting
from pyscf import grad     # Energy gradients
from pyscf import geomopt  # Geometry optimization
from pyscf import solvent  # Solvent effects
from pyscf.pbc import gto as pbcgto  # GTOs for periodic systems
from pyscf.pbc import scf as pbcscf  # SCF for periodic systems
from pyscf.pbc import dft as pbcdft  # DFT for periodic systems
```

## `gto` Module: Molecules and Basis Sets

### `M` Class: Molecule/Cell Objects

#### Constructor
```python
mol = gto.M(
    atom=None,              # Atomic coordinates
    basis=None,             # Basis set
    spin=0,                 # Spin, 2S
    charge=0,               # Total charge
    symmetry=False,         # Symmetry
    unit='Ang',             # Coordinate unit ('Ang' or 'Bohr')
    parse_arg=False,        # Whether to parse arguments
    verbose=4,              # Output verbosity
)
```

#### Common Attributes
```python
mol.atom          # Atom list
mol.natm          # Number of atoms
mol.nel           # Number of electrons
mol.nao           # Number of AOs
mol.nelectron     # Number of electrons
mol.spin          # Spin, 2S
mol.charge        # Charge
mol.basis         # Basis-set name
mol.symmetry      # Symmetry
mol.irrep_id      # Irreducible-representation IDs
mol.symm_orb      # Symmetry-adapted orbitals
mol.cart          # Whether Cartesian Gaussian functions are used
```

#### Common Methods

**Molecular Geometry**
```python
# Cartesian-coordinate format
mol = gto.M(
    atom='''
    O  0.0  0.0  0.0
    H  0.0  0.76 0.59
    H  0.0 -0.76 0.59
    ''',
    basis='cc-pvdz'
)

# Z-matrix format
mol = gto.M(
    atom='''
    O
    H 1 0.96
    H 1 0.96 2 104.5
    ''',
    basis='cc-pvdz'
)

# Read from an XYZ file
mol = gto.M(atom='h2o.xyz', basis='cc-pvdz')
```

**Basis-Set Configuration**
```python
# Use a built-in basis set
mol = gto.M(atom='...', basis='cc-pvdz')

# Define a custom basis set
mol = gto.M(
    atom='...',
    basis={
        'H': gto.parse('''
    H    S
      3.42525091              0.15432897
      0.62391373              0.53532814
      0.16885540              0.44463454
    '''),
        'O': gto.parse('...')
    }
)

# Mix basis sets
mol = gto.M(
    atom='O 0 0 0; H 0 1 0; H 0 0 1',
    basis={'O': 'cc-pvdz', 'H': 'cc-pvsz'}
)
```

**Symmetry**
```python
# Enable symmetry
mol = gto.M(atom='...', symmetry=True)

# Specify a point group
mol = gto.M(atom='...', symmetry='d2h')

# Available point groups: 'd2h', 'c2v', 'c2h', 'd2', 'ci', 'cs', 'c1'
```

**Integral Evaluation**
```python
# Overlap matrix
S = mol.intor('int1e_ovlp')

# Kinetic-energy matrix
T = mol.intor('int1e_kin')

# Nuclear-attraction matrix
V = mol.intor('int1e_nuc')

# Nuclear-repulsion energy
Vnn = mol.energy_nuc()

# Two-electron integrals (eightfold symmetry)
eri = mol.intor('int2e')

# Two-electron integrals (fourfold symmetry)
eri_4 = mol.intor('int2e_sph', aosym=4)
```

**Periodic Cell (PBC)**
```python
from pyscf.pbc import gto as pbcgto

cell = pbcgto.M(
    atom='...',
    basis='gth-szv',
    a=[[3.5, 0, 0], [0, 3.5, 0], [0, 0, 3.5]],  # Lattice vectors
    pseudo='gth-pade',  # Pseudopotential
    ke_cutoff=100,      # Kinetic-energy cutoff
)

# k-point sampling
kpts = cell.make_kpts([2, 2, 2])
```

## `scf` Module: SCF Methods

### SCF Base Classes

#### Common Classes
```python
from pyscf import dft, scf

# Hartree-Fock methods
mf = scf.RHF(mol)  # Restricted HF
mf = scf.UHF(mol)  # Unrestricted HF
mf = scf.ROHF(mol) # Restricted open-shell HF

# Kohn-Sham DFT methods
mf = dft.RKS(mol)  # Restricted KS-DFT
mf = dft.UKS(mol)  # Unrestricted KS-DFT
```

#### Core Attributes
```python
mf.mol           # Molecule object
mf.mo_coeff      # MO coefficients (nao, nmo)
mf.mo_energy     # MO energies
mf.mo_occ        # MO occupations
mf.e_tot         # Total energy
mf.converged     # Convergence status
mf.scf_summary   # SCF summary dictionary
mf.chkfile       # Checkpoint file
```

#### Core Methods

**Run SCF**
```python
# Standard run
mf.kernel()
mf.kernel(dm0=None)  # Specify an initial density matrix

# Return the energy and inspect the convergence flag
e_conv = mf.kernel()
print(f'Energy: {e_conv:.6f}')
print(f'Converged: {mf.converged}')
```

**Density Matrices**
```python
# Total density matrix
dm = mf.make_rdm1()

# Alpha and beta density matrices (open-shell)
dma, dmb = mf.make_rdm1()
```

**Fock Matrices**
```python
# Effective Fock matrix
fock = mf.get_fock()

# Fock matrix in the AO basis
fock_ao = mf.get_fock(dm=mf.make_rdm1())
```

**Orbital Analysis**
```python
# HOMO energy
homo_idx = np.where(mf.mo_occ > 0)[0][-1]
homo_energy = mf.mo_energy[homo_idx]

# LUMO energy
lumo_idx = np.where(mf.mo_occ == 0)[0][0]
lumo_energy = mf.mo_energy[lumo_idx]

# Energy gap
gap = lumo_energy - homo_energy
```

#### SCF Control

**Initial Guess**
```python
# Core Hamiltonian
mf.init_guess = '1e'

# Hückel method
mf.init_guess = 'huckel'

# Read from a checkpoint
mf.init_guess = 'chkfile'
mf.chkfile = 'previous.chk'

# Start from a density matrix
dm_guess = mf.get_init_guess()
mf.kernel(dm0=dm_guess)
```

**Convergence Acceleration**
```python
# DIIS
mf.diis_start_cycle = 3  # DIIS starting cycle
mf.diis = True

# Damping
mf.damp_factor = 0.2  # Damping factor

# Level shift
mf.level_shift = 0.5  # Virtual-orbital shift

# Maximum number of iterations
mf.max_cycle = 100

# Convergence thresholds
mf.conv_tol = 1e-8  # Energy threshold
mf.conv_tol_grad = 1e-5  # Gradient threshold
```

**Save/Load**
```python
# Save to a checkpoint
mf.chkfile = 'calc.chk'
mf.dump_chk(mf.chkfile)  # Save the current SCF object manually

# Read from a checkpoint
mf = scf.RHF(mol)
mf.chkfile = 'calc.chk'
mf.kernel()  # Read automatically

# Read only the density
mf_chk = scf.RHF(mol)
dm = mf_chk.from_chk('calc.chk')
```

## `dft` Module: DFT Methods

### RKS/UKS Classes

#### Basic Usage
```python
from pyscf import dft

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()
```

#### Functional Configuration

**Built-In Functionals**
```python
# LDA
mf.xc = 'svwn'
mf.xc = 'vwn'

# GGA
mf.xc = 'pbe'
mf.xc = 'blyp'
mf.xc = 'bop'
mf.xc = 'bp86'

# meta-GGA
mf.xc = 'tpss'
mf.xc = 'm06-l'
mf.xc = 'scan'

# Hybrid functionals
mf.xc = 'b3lyp'
mf.xc = 'pbe0'
mf.xc = 'wb97x-d'
mf.xc = 'm06-2x'

# Range-separated functionals
mf.xc = 'cam-b3lyp'
mf.xc = 'wb97x'
mf.xc = 'lc-wpbe'
```

**Custom Functionals**
```python
# Expression-based definition
mf.xc = '.2*HF + .08*LDA + .72*B88, .81*LYP + .19*VWN'

# Fully custom function
def my_xc(rho, spin=0, deriv=1, **kwargs):
    # rho: (N, 4) array [density, grad_x, grad_y, grad_z]
    # Return: (ex+vc, vrho, vgamma)
    e_xc = ...  # Exchange-correlation energy density
    vrho = ...  # Density potential
    vgamma = ...  # Gradient potential
    return e_xc, vrho, vgamma

mf._numint._xc = my_xc
```

#### Grid Configuration

**Atomic Grid**
```python
# (number of angular grid points, number of radial grid points)
mf.grids.atom_grid = (75, 302)  # Standard
mf.grids.atom_grid = (99, 590)  # Fine
mf.grids.atom_grid = (250, 974) # Ultrafine

# Grid pruning
mf.grids.prune = None  # No pruning
mf.grids.prune = dft.gen_grid.sg1_prune  # SG1 pruning
```

**Nonlocal-Correlation Grid**
```python
# Used for nonlocal dispersion
mf.nlc = 'vv10'  # Or 'rvv10'
mf.nlcgrids.atom_grid = (50, 194)
```

#### Molecular Integrals

**Effective Potential and Coulomb/Exchange Matrices**
```python
# Effective potential used by the Kohn-Sham build. It includes Coulomb and
# exchange-correlation contributions and, for hybrids, exact exchange.
veff = mf.get_veff(mol, mf.make_rdm1())

# Coulomb and exact-exchange matrix builders; these are not an
# exchange-correlation decomposition of veff.
vj = mf.get_j(mol, mf.make_rdm1())  # Coulomb matrix
vk = mf.get_k(mol, mf.make_rdm1())  # Exchange matrix
```

## `tdscf` Module: Time-Dependent SCF

### TDDFT Class

#### Basic Usage
```python
from pyscf import tdscf

# TDDFT
td = tdscf.TDDFT(mf)
td.nstates = 6
td.kernel()

# TDA approximation
td_tda = tdscf.TDA(mf)
td_tda.nstates = 4
td_tda.kernel()
```

#### Core Attributes
```python
td.nstates         # Number of excited states
td.converged       # Convergence status for the requested states
td.e               # Excitation energies (Hartree)
td.xy              # Excitation/de-excitation amplitudes
```

#### Extract Excited-State Information

```python
# Excitation energies (eV)
for i, e in enumerate(td.e):
    print(f'State {i+1}: {e*27.2114:.2f} eV')

# Wavelengths (nm)
for i, e in enumerate(td.e):
    wavelength = 1240/(e*27.2114)
    print(f'State {i+1}: {wavelength:.1f} nm')

# Oscillator strengths
oscillator_strengths = td.oscillator_strength()
for i, f in enumerate(oscillator_strengths):
    print(f'State {i+1}: f = {f:.3f}')

# Transition dipole moments
transition_dipoles = td.transition_dipole()
for i, transition_dipole in enumerate(transition_dipoles):
    print(f'State {i+1}: μ = {transition_dipole}')
```

#### NTO Analysis

```python
# Natural transition orbitals
weights, nto_coeff = td.get_nto(state=1)

# For a restricted reference, occupied NTOs precede virtual NTOs.
nocc = np.count_nonzero(mf.mo_occ > 0)
occupied_ntos = nto_coeff[:, :nocc]
virtual_ntos = nto_coeff[:, nocc:]

# Dominant occupied-virtual NTO-pair weight
dominant_weight = weights.max()
print(f'Dominant weight: {dominant_weight:.3f}')
```

#### State Averaging

```python
# State-averaged TDDFT
td = tdscf.TDDFT(mf)
td.state_average = True
td.nstates = 5
td.kernel()
```

## `mcscf` Module: Multiconfigurational SCF

### CASSCF Class

#### Basic Usage
```python
from pyscf import mcscf

# CAS(nele, norb): number of electrons, number of orbitals
cas = mcscf.CASSCF(mf, 6, 8)
e_cas, ci_vec, mo, mo_occ = cas.kernel()
```

#### Core Attributes
```python
cas.ncas           # Number of active orbitals
cas.nelecas        # Number of active electrons
cas.mo_coeff       # Optimized MO coefficients
cas.ci             # CI vector
cas.e_cas          # CASSCF energy
cas.frozen         # Number of frozen orbitals
```

#### Initial-Orbital Selection

```python
# Automatic selection (near the HOMO/LUMO)
cas = mcscf.CASSCF(mf, 6, 8)

# Manual specification
cas = mcscf.CASSCF(mf, 6, 8)
cas.mo_coeff = mf.mo_coeff  # Start from SCF orbitals

# Natural orbitals (recommended)
no_coeff, no_occ = mcscf.cas_natorb(cas)
cas.mo_coeff = no_coeff
```

#### State Averaging

```python
# State-averaged CASSCF
cas = mcscf.CASSCF(mf, 6, 8)
cas.state_average([0.5, 0.3, 0.2])  # Three states with weights 0.5/0.3/0.2
cas.kernel()
```

#### Excited States

```python
# State-specific CASSCF
cas = mcscf.CASSCF(mf, 6, 8)
cas.state_specific(2)  # Third excited state
cas.kernel()
```

#### Density-Fitted CASSCF

```python
# Density-fitting acceleration
cas_df = mcscf.DFCASSCF(mf, 6, 8, auxbasis='cc-pvtz-jkfit')
cas_df.kernel()
```

## `grad` Module: Energy Gradients

### Gradient Methods

#### HF/DFT Gradients

```python
from pyscf import grad

# Get the gradient method
g = mf.nuc_grad_method()
grad = g.kernel()  # (natm, 3)

# Print the gradient
print('Energy gradient:')
for i in range(mol.natm):
    atom = mol.atom_symbol(i)
    print(f'{atom}: {grad[i]}')
```

#### CASSCF Gradients

```python
# CASSCF gradient
g = cas.nuc_grad_method()
grad = g.kernel()
```

#### Gradient-Based Optimization

```python
# Use gradients for geometry optimization
scanner = mf.as_scanner()
optimizer = scanner.nuc_grad_method().optimizer()
optimized_mol = optimizer.kernel()
```

## `ao2mo` Module: Integral Transformation

### AO-to-MO Transformation

```python
from pyscf import ao2mo

# Full transformation
eri_mo = ao2mo.incore.full(eri_ao, mo_coeff)

# Partial transformation (saves memory)
eri_mo = ao2mo.incore.general(eri_ao, [mo_occ, mo_occ, mo_occ, mo_occ])

# Out-of-core transformation (large molecules)
ao2mo.outcore.full(mol, mo_coeff, 'eri_mo.h5')
```

## `df` Module: Density Fitting

### Density-Fitted SCF

```python
from pyscf import df

# Automatic density fitting
mf_df = mf.density_fit(auxbasis='def2-universal-jfit')

# Manual density fitting
mf_df = df.density_fit(scf.RHF(mol), auxbasis='cc-pvtz-jkfit')

# Common auxiliary basis sets
auxbasis = 'def2-universal-jfit'  # General-purpose auxiliary basis
auxbasis = 'cc-pvtz-jkfit'       # JK-fitting basis
auxbasis = 'cc-pvtz-ri'          # RI-fitting basis
```

## Utility Functions

### Coordinate and File Operations

```python
from pyscf.tools import cubegen, molden

# Cube files (for visualization)
cubegen.orbital(mol, 'h2o_homo.cube',
                mo_coeff[:, mo_occ>0][:, -1])
cubegen.density(mol, 'h2o_density.cube',
                mf.make_rdm1())

# Molden file (orbital visualization)
molden.from_mo(mol, 'h2o.molden', mf.mo_coeff)
```

### Population Analysis

```python
from pyscf import lo

# Mulliken population
pop = mf.mulliken_pop(mol, mf.make_rdm1())

# Löwdin population
pop = lo.vec_lowdin(mf.mo_coeff, mf.get_ovlp())

# Natural population analysis (NPA)
# Requires an external tool or a manual implementation
```

## Parallel Computing

### OpenMP Parallelism

```python
import os

# Set the number of threads
os.environ['OMP_NUM_THREADS'] = '8'

# Run with automatic parallelism
mf.kernel()
```

### MPI Parallelism

```python
# Run through mpi4py
# mpirun -np 4 python script.py

from pyscf import lib

# Inspect MPI settings
print(lib.num_threads())  # Number of threads
print(lib.nproc)          # Number of processes (when using MPI)
```

## Common Troubleshooting

### SCF Does Not Converge

```python
# Inspect convergence information
print(mf.converged)
print(mf.scf_summary)

# Common remedies
mf.init_guess = 'huckel'   # Use a better initial guess
mf.diis_start_cycle = 5    # Delay DIIS
mf.level_shift = 0.5       # Apply a level shift
mf.damp_factor = 0.2       # Apply damping

# Newton-Raphson method
mf_newton = mf.newton()
mf_newton.kernel()
```

### Insufficient Memory

```python
# Density fitting
mf = mf.density_fit()

# Out-of-core calculation
mf.direct_scf = False
mf.chkfile = 'calc.chk'

# Reduce the grid size
mf.grids.atom_grid = (75, 302)
```

### Excited-State Problems

```python
# Try TDA first
td_tda = tdscf.TDA(mf)

# Check the quality of the ground state
print(f'Gap: {gap*27.2114:.2f} eV')

# Increase grid accuracy
mf.grids.atom_grid = (99, 590)
```

## Checkpoint Files

### Save All Data

```python
# Save automatically
mf.chkfile = 'calculation.chk'
mf.kernel()

# Save manually
mf.dump_chk(mf.chkfile)
```

### Read from a Checkpoint

```python
# Read the molecule and saved SCF record.
mol_chk, scf_data = scf.chkfile.load_scf('calculation.chk')
mo_coeff = scf_data['mo_coeff']
mo_energy = scf_data['mo_energy']
mo_occ = scf_data['mo_occ']

# Reconstruct the RHF/RKS density matrix from the saved orbitals.
mf_chk = scf.RHF(mol_chk)
dm = mf_chk.make_rdm1(mo_coeff, mo_occ)
```

## Performance-Optimization Suggestions

### Memory-versus-Speed Tradeoff

```python
# More memory: store all integrals
mf.direct_scf = False

# Less memory: evaluate integrals on demand
mf.direct_scf = True
```

### Basis-Set Selection

```python
# Fast preliminary calculation
basis = 'sto-3g'

# Standard calculation
basis = 'def2-svp'

# High accuracy
basis = 'def2-tzvp'

# Highest accuracy
basis = 'def2-qzvp'
```

### Parallel Efficiency

```python
# SCF parallelism: good efficiency
mf.kernel()

# Density fitting: good parallel efficiency
mf_df = mf.density_fit()

# MP2/CCSD: memory-limited, moderate parallel efficiency
```

## Extending PySCF

### Custom Hamiltonian

```python
import numpy as np

# 1D Hubbard model
h1 = np.zeros((N, N))
for i in range(N-1):
    h1[i, i+1] = h1[i+1, i] = -1.0

eri = np.zeros((N, N, N, N))
for i in range(N):
    eri[i, i, i, i] = U

# Construct an SCF calculation
mol = gto.M()
mol.nelectron = N
mf = scf.RHF(mol)
mf.get_hcore = lambda *args: h1
mf.get_ovlp = lambda *args: np.eye(N)
mf._eri = ao2mo.restore(8, eri, N)
mf.kernel()
```

### Integration with JAX

```python
import jax
import jax.numpy as jnp
from pyscf.jax import scf as jax_scf

# Convert an SCF object
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol)
mf_jax = jax_scf.RHF(mol)

# Automatic-differentiation gradient
grad_func = jax.grad(lambda coords: mf_jax.kernel()[0])
grad = grad_func(coords)
```

## References

1. PySCF documentation: https://pyscf.org/
2. PySCF GitHub: https://github.com/pyscf/pyscf
3. PySCF examples: https://github.com/pyscf/pyscf/tree/master/examples
4. Sun, Q. et al. "Recent developments in the PySCF program package." WIREs Comput. Mol. Sci. 2020.
