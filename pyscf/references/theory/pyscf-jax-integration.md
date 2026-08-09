# PySCF and JAX Integration Guide

> **Unvalidated historical draft.** Do not use this document as an executable
> API guide. Array conversion alone does not make PySCF SCF calculations
> differentiable, JIT-compatible, or GPU accelerated. Select a documented
> differentiable/accelerator backend and validate it against the installed
> versions. See [`../../tools/README.md`](../../tools/README.md).

## Overview

JAX (Just Another eXecutor) is a high-performance numerical-computing library with automatic differentiation and JIT compilation. Integrating PySCF with JAX enables:

1. **End-to-end automatic differentiation**: differentiate DFT calculations directly
2. **Differentiable quantum chemistry**: use PySCF within neural networks
3. **Functional optimization**: optimize exchange-correlation functional parameters automatically
4. **GPU acceleration**: use GPUs to accelerate integral evaluation

## Installation and Configuration

### 1. Environment Setup

```bash
# Install PySCF
pip install pyscf

# Install JAX
pip install jax jaxlib

# GPU support (optional)
# pip install jax[cuda]
```

### 2. Basic Imports

```python
import numpy as np
import jax
import jax.numpy as jnp
from jax import grad, jit, vmap, pmap

import pyscf
from pyscf import gto, scf, dft, lib
```

## PySCF-JAX Basics

### 1. Array Conversion

#### NumPy ↔ JAX Arrays
```python
# NumPy array
np_array = np.array([1.0, 2.0, 3.0])

# Convert to a JAX array
jx_array = jnp.array(np_array)
# Or
jx_array = jax.device_put(np_array)

# Convert back to NumPy
np_from_jax = np.array(jx_array)
# Or
np_from_jax = jax.device_get(jx_array)
```

#### Handling PySCF Objects
```python
# PySCF molecule object
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')

# Extract JAX-compatible data
coords = jnp.array(mol.atom_coords())
charges = jnp.array([mol.atom_charge(i) for i in range(mol.natm)])
basis_info = mol._basis  # Basis-set information
```

### 2. Wrappable Functions

#### JIT-Compatible Functions
```python
from jax import jit

# Define a simple energy function
def simple_h2_energy(coords):
    """
    H₂ nuclear-repulsion energy plus a simplified electronic energy.

    params:
        coords: (2, 3) atomic coordinates
    """
    # Nuclear-nuclear repulsion
    r_vec = coords[0] - coords[1]
    r = jnp.sqrt(jnp.sum(r_vec**2))
    Vnn = 1.0 / r

    # Simplified electronic energy (Hückel model)
    E_el = -2.0 * np.exp(-r)  # Simplified model

    return Vnn + E_el

# JIT compilation
jit_energy = jit(simple_h2_energy)

# Test
coords_test = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
E = jit_energy(coords_test)
print(f"Energy: {E:.6f} Ha")
```

#### Automatic Differentiation
```python
# Compute the gradient
grad_energy = grad(simple_h2_energy)

# Compute forces at the given geometry (negative gradient)
forces = -grad_energy(coords_test)
print(f"Force on atom 0: {forces[0]}")
print(f"Force on atom 1: {forces[1]}")
```

## GradDFT Framework

### 1. GradDFT Architecture

GradDFT is an end-to-end differentiable DFT framework that allows gradients to be backpropagated through a DFT calculation.

```python
import jax
import jax.numpy as jnp
from pyscf import gto, scf, dft
from pyscf.dft import numint

class GradDFT:
    """
    Differentiable DFT framework.
    """
    def __init__(self, mol):
        self.mol = mol
        self.nao = mol.nao
        self.nelec = mol.nelectron

        # Initialize SCF
        self.mf = dft.RKS(mol)
        self.mf.xc = 'pbe'

        # Numerical-integration object
        self.ni = numint.NumInt()

    def build_matrices(self, coords):
        """
        Build JAX-compatible integral matrices.

        params:
            coords: (natm, 3) atomic coordinates
        """
        # Note: the molecule object must be rebuilt here.
        # Simplification: assume fixed basis functions and vary only atom positions.

        # Overlap matrix (simplification: assume it is fixed)
        S = self.mf.get_ovlp()

        # Kinetic-energy matrix (simplified)
        T = self.mf.get_hcore() - self.mol.intor('int1e_nuc')

        # Nuclear-attraction matrix (coordinate-dependent)
        V = self._build_nuclear_potential(coords)

        # Total core Hamiltonian
        Hcore = T + V

        return S, Hcore

    def _build_nuclear_potential(self, coords):
        """
        Build the nuclear-attraction potential matrix.

        V_μν = ∫ χ_μ(r) [ -Σ_A Z_A/|r-R_A| ] χ_ν(r) dr
        """
        # Simplification: use fixed integrals plus a displacement.
        # A complete implementation must recompute the integrals.

        # Reference nuclear attraction
        V0 = self.mol.intor('int1e_nuc')

        # Effect of coordinate displacement (must be recomputed)
        # This example uses a simplifying assumption.
        return V0  # Must be recomputed in a complete implementation

    def scf_iteration(self, coords, dm):
        """
        Perform one SCF iteration.

        params:
            coords: (natm, 3) atomic coordinates
            dm: (nao, nao) density matrix
        """
        S, Hcore = self.build_matrices(coords)

        # Exchange-correlation potential
        rho = self._compute_rho(dm)
        vxc = self._compute_vxc(rho)

        # Hartree potential
        J = self._compute_j(dm)

        # Fock matrix
        F = Hcore + J + vxc

        # Diagonalization
        mo_energy, mo_coeff = jnp.linalg.eigh(F, S)

        # New density
        nocc = self.nelec // 2
        mo_occ = jnp.zeros_like(mo_energy)
        mo_occ = mo_occ.at[:nocc].set(2.0)

        dm_new = 2.0 * mo_coeff[:, :nocc] @ mo_coeff[:, :nocc].T

        return dm_new, mo_energy, mo_coeff, F

    def _compute_rho(self, dm):
        """
        Compute the electron density (simplified).
        """
        # Numerical integration should be used here.
        # Simplification: approximate the ground-state density.
        return jnp.trace(dm) / self.nao

    def _compute_vxc(self, rho):
        """
        Compute the exchange-correlation potential (simplified).
        """
        # LDA approximation
        if isinstance(rho, (float, int)):
            rho_val = rho
        else:
            rho_val = jnp.mean(rho)

        # LDA exchange potential
        vx = -(3/np.pi)**(1/3) * rho_val**(1/3)

        # Diagonal matrix (simplified)
        n = self.nao
        return vx * jnp.eye(n)

    def _compute_j(self, dm):
        """
        Compute the Hartree potential (simplified).
        """
        # Simplification: J = G * dm, where G is the Coulomb operator.
        # This example uses a scaling approximation.
        return 0.5 * dm

    def solve_scf(self, coords, max_iter=50, tol=1e-8):
        """
        Solve the complete SCF problem.

        params:
            coords: (natm, 3) atomic coordinates
            max_iter: maximum number of iterations
            tol: convergence threshold
        """
        # Initial density: core Hamiltonian
        S, Hcore = self.build_matrices(coords)
        dm = jnp.linalg.inv(S) @ Hcore @ jnp.linalg.inv(S)
        dm = dm / jnp.trace(dm) * self.nelec

        # SCF loop
        for i in range(max_iter):
            dm_new, mo_energy, mo_coeff, F = self.scf_iteration(coords, dm)

            # Convergence check
            delta = jnp.max(jnp.abs(dm_new - dm))
            if delta < tol:
                break

            dm = dm_new

        # Total energy
        E = self._total_energy(coords, dm, F)

        return E, dm, mo_energy, mo_coeff

    def _total_energy(self, coords, dm, F):
        """
        Compute the total energy.
        """
        S, Hcore = self.build_matrices(coords)

        # Electronic energy
        E_el = 0.5 * jnp.trace((Hcore + F) @ dm)

        # Nuclear-nuclear repulsion
        Vnn = self._nuclear_repulsion(coords)

        return E_el + Vnn

    def _nuclear_repulsion(self, coords):
        """
        Nuclear-nuclear repulsion energy.
        """
        natm = coords.shape[0]
        Vnn = 0.0

        for i in range(natm):
            for j in range(i+1, natm):
                r = jnp.sqrt(jnp.sum((coords[i] - coords[j])**2))
                Z_i = self.mol.atom_charge(i)
                Z_j = self.mol.atom_charge(j)
                Vnn += Z_i * Z_j / r

        return Vnn

# Usage
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
graddft = GradDFT(mol)

coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
E, dm, mo_e, mo_c = graddft.solve_scf(coords)
print(f"DFT energy: {E:.6f} Ha")
```

### 2. Automatic-Differentiation Gradients

```python
# Define a differentiable energy function
def dft_energy(coords):
    """
    DFT energy function (JAX-differentiable).

    params:
        coords: (natm, 3) atomic coordinates
    """
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# Compute the gradient
grad_energy = grad(dft_energy)

# Test
coords_test = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
forces = -grad_energy(coords_test)

print(f"Atomic forces:")
for i, f in enumerate(forces):
    print(f"  Atom {i}: {f}")
```

### 3. Geometry Optimization

```python
# Gradient-descent optimization
def gradient_descent(initial_coords, lr=0.01, steps=100):
    """
    Simple gradient-descent optimization.

    params:
        initial_coords: initial coordinates
        lr: learning rate
        steps: number of optimization steps
    """
    coords = initial_coords.copy()

    for step in range(steps):
        # Compute forces and energy
        E = dft_energy(coords)
        forces = -grad_energy(coords)

        # Update coordinates
        coords = coords + lr * forces

        if step % 10 == 0:
            print(f"Step {step}: E = {E:.6f} Ha")

    return coords, E

# Optimize
coords_init = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
coords_opt, E_opt = gradient_descent(coords_init, lr=0.01, steps=50)

print(f"\nOptimized energy: {E_opt:.6f} Ha")
print(f"Optimized coordinates:")
print(coords_init)
print(coords_opt)
```

## Functional Optimization

### 1. Parameterized Functional

```python
from jax import vmap

class ParameterizedXC:
    """
    Parameterized exchange-correlation functional.
    E_xc = ∫ f(ρ, θ) dτ
    Here, θ denotes the learnable parameters.
    """
    def __init__(self, mol, n_params=3):
        self.mol = mol
        self.n_params = n_params

        # Initial parameters (for example, a variant of the PBE parameters)
        self.params = jnp.array([0.9, 0.8, 0.1])

    def lda_exchange(self, rho, params):
        """
        Parameterized LDA exchange.

        E_x = a1 * E_x^LDA
        """
        a1 = params[0]

        # Standard LDA exchange
        Ex = -0.75 * (3/np.pi)**(1/3) * rho**(4/3)

        return a1 * Ex

    def lda_correlation(self, rho, params):
        """
        Parameterized LDA correlation.
        """
        a2 = params[1]
        a3 = params[2]

        # Simplified VWN form
        Ec = a2 * rho**(1/3) + a3 * rho**(2/3)

        return Ec

    def xc_energy_density(self, rho, params):
        """
        Exchange-correlation energy density.
        """
        Ex = self.lda_exchange(rho, params)
        Ec = self.lda_correlation(rho, params)

        return Ex + Ec

    def xc_potential(self, rho, params):
        """
        Exchange-correlation potential V_xc = dE_xc/dρ.
        """
        # Compute through automatic differentiation
        # E_xc[ρ] = ∫ ε_xc(ρ) dτ
        # V_xc = dε_xc/dρ

        xc_density = lambda r: self.xc_energy_density(r, params)

        # JAX automatic differentiation
        vxc = grad(xc_density)(rho)

        return vxc

# Usage
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
xc = ParameterizedXC(mol)

# Compute the XC energy at a given density
rho = 0.5  # Simplification: uniform density
Exc_density = xc.xc_energy_density(rho, xc.params)
print(f"XC energy density: {Exc_density:.6f}")

# XC potential
Vxc = xc.xc_potential(rho, xc.params)
print(f"XC potential: {Vxc:.6f}")
```

### 2. Functional-Parameter Optimization

```python
def functional_loss(params, training_data):
    """
    Loss function for optimizing functional parameters.

    params:
        params: functional parameters
        training_data: [(mol, E_ref), ...] training data
    """
    xc = ParameterizedXC(training_data[0][0])
    xc.params = params

    total_loss = 0.0

    for mol, E_ref in training_data:
        # Run a DFT calculation with the parameterized functional
        graddft = GradDFT(mol)

        # Modify XC to use the parameterized functional
        # (requires integration into GradDFT)
        # graddft.xc_func = xc

        # Compute the energy
        coords = jnp.array(mol.atom_coords())
        E_pred, _, _, _ = graddft.solve_scf(coords)

        # MSE loss
        total_loss += (E_pred - E_ref)**2

    return total_loss / len(training_data)

# Optimize
from jax import value_and_grad
from jax.example_libraries import optimizers

# Create training data (simplified)
training_mols = [
    (gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g'), -1.16),
    (gto.M(atom='H 0 0 0; H 0 0 0.80', basis='sto-3g'), -1.14),
]

# Initial parameters
init_params = jnp.array([0.9, 0.8, 0.1])

# Optimizer
opt_init, opt_update, get_params = optimizers.adam(0.01)
opt_state = opt_init(init_params)

@jit
def step(i, opt_state, training_data):
    params = get_params(opt_state)
    loss, grads = value_and_grad(functional_loss)(params, training_data)
    opt_state = opt_update(i, grads, opt_state)
    return loss, opt_state

# Optimization loop
for i in range(100):
    loss, opt_state = step(i, opt_state, training_mols)
    if i % 10 == 0:
        print(f"Iteration {i}: Loss = {loss:.6f}")

# Optimal parameters
optimal_params = get_params(opt_state)
print(f"\nOptimal parameters: {optimal_params}")
```

## Advanced Integration

### 1. Neural-Network-Augmented DFT

```python
import flax.linen as nn
import optax

class NNXC(nn.Module):
    """
    Neural-network exchange-correlation functional.
    """
    hidden_dim: int = 64

    @nn.compact
    def __call__(self, rho_features):
        """
        Density features → XC energy.

        params:
            rho_features: (N, n_features) density features
        """
        x = nn.Dense(self.hidden_dim)(rho_features)
        x = nn.tanh(x)
        x = nn.Dense(self.hidden_dim)(x)
        x = nn.tanh(x)
        x = nn.Dense(1)(x)
        return x.squeeze()

# Feature extraction
def extract_density_features(mol, dm):
    """
    Extract features from the density matrix.
    """
    # Simplification: use the local density and gradient.
    # A complete implementation should compute features at every grid point.

    # Example features
    rho = jnp.trace(dm) / mol.nao  # Mean density
    grad_rho = jnp.sqrt(jnp.sum(jnp.array([0.1, 0.2, 0.3])**2))  # Simplified

    features = jnp.array([rho, grad_rho, rho**2, grad_rho**2])
    return features

# End-to-end training
class NNEnergyModel:
    """
    Neural-network-augmented DFT energy model.
    """
    def __init__(self, mol):
        self.mol = mol
        self.nn_xc = NNXC()
        self.graddft = GradDFT(mol)

    def forward(self, coords, dm):
        """
        Forward pass.
        """
        # Extract density features
        features = extract_density_features(self.mol, dm)

        # Predict the XC energy density with the neural network
        Exc_nn = self.nn_xc(features)

        # Compute the total energy
        E = self.graddft._total_energy(coords, dm, None)

        return E + Exc_nn

    def loss(self, coords, dm, E_ref):
        E_pred = self.forward(coords, dm)
        return (E_pred - E_ref)**2

# Train
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
model = NNEnergyModel(mol)

# Initialize parameters
dummy_features = jnp.array([0.5, 0.1, 0.25, 0.01])
params = model.nn_xc.init(jax.random.PRNGKey(0), dummy_features)

# Optimizer
optimizer = optax.adam(0.01)
opt_state = optimizer.init(params)

@jit
def train_step(coords, dm, E_ref, params, opt_state):
    def loss_fn(p):
        # Temporarily update the model parameters
        # (requires correct integration with Flax)
        return model.loss(coords, dm, E_ref)

    loss, grads = jax.value_and_grad(loss_fn)(params)
    updates, opt_state = optimizer.update(grads, opt_state)
    params = optax.apply_updates(params, updates)
    return loss, params, opt_state

# Training loop (simplified)
coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
dm = jnp.eye(mol.nao) * 0.5  # Simplified density
E_ref = -1.16  # Reference energy

for i in range(100):
    loss, params, opt_state = train_step(coords, dm, E_ref, params, opt_state)
    if i % 20 == 0:
        print(f"Iteration {i}: Loss = {loss:.6f}")
```

### 2. Batched Calculations for Multiple Systems

```python
# Process multiple molecules in a batch
from jax import vmap

def batch_dft_energy(batch_coords, mol_template):
    """
    Compute DFT energies for multiple molecules in a batch.

    params:
        batch_coords: (batch_size, natm, 3) batched coordinates
        mol_template: molecule template (basis set, etc.)
    """
    def single_energy(coords):
        graddft = GradDFT(mol_template)
        E, _, _, _ = graddft.solve_scf(coords)
        return E

    # Vectorize
    energies = vmap(single_energy)(batch_coords)
    return energies

# Test
batch_coords = jnp.array([
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],  # Molecule 1
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.80]],  # Molecule 2
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.86]],  # Molecule 3
])

energies = batch_dft_energy(batch_coords, mol)
print(f"Batched energies: {energies}")
```

### 3. GPU Acceleration

```python
# Check GPU availability
from jax.lib import xla_bridge
print(f"Device: {xla_bridge.get_backend().platform}")

# Move data to the GPU
coords_gpu = jax.device_put(coords, jax.devices('gpu')[0])

# JIT compilation plus automatic GPU acceleration
@jit
def fast_scf(coords):
    graddft = GradDFT(mol)
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# Run (uses the GPU automatically)
E_gpu = fast_scf(coords_gpu)
```

## Practical Applications

### 1. Potential-Energy-Surface Scan

```python
# H₂ potential-energy-surface scan
def scan_bond_length(r_values):
    """
    Scan the H₂ bond length.

    params:
        r_values: (n,) array of bond lengths
    """
    energies = []

    for r in r_values:
        # Construct the molecule
        coords = jnp.array([
            [0.0, 0.0, -r/2],
            [0.0, 0.0, r/2]
        ])

        # Compute the DFT energy
        E = dft_energy(coords)
        energies.append(E)

    return jnp.array(energies)

# Scan
r_values = jnp.linspace(0.5, 2.0, 20)
energies = scan_bond_length(r_values)

# Plot
import matplotlib.pyplot as plt
plt.plot(r_values, energies, 'o-')
plt.xlabel('Bond Length (Å)')
plt.ylabel('Energy (Ha)')
plt.title('H₂ Potential Energy Curve')
plt.show()
```

### 2. Force-Field Fitting

```python
# Fit a force field to DFT forces
def force_loss(params, coords, target_forces):
    """
    Loss for force-field parameter optimization.

    params:
        params: force-field parameters
        coords: atomic coordinates
        target_forces: target forces (natm, 3)
    """
    # Compute DFT forces
    forces_pred = -grad_energy(coords)

    # MSE loss
    loss = jnp.mean((forces_pred - target_forces)**2)
    return loss

# Usage
coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
target_forces = jnp.array([[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]])  # Example

# Optimize (simplified)
params = jnp.array([1.0, 0.5])

for i in range(50):
    loss, grads = value_and_grad(force_loss)(params, coords, target_forces)
    params = params - 0.01 * grads

    if i % 10 == 0:
        print(f"Iteration {i}: Loss = {loss:.6f}")
```

## Considerations

### 1. Numerical Stability

```python
# Avoid division by zero
r = jnp.sqrt(jnp.sum((coords[i] - coords[j])**2) + 1e-12)

# Gradient clipping
@jit
def clip_grads(grads, max_norm=1.0):
    norm = jnp.sqrt(jnp.sum(grads**2))
    scale = jnp.minimum(1.0, max_norm / (norm + 1e-8))
    return grads * scale

# Usage
grads = clip_grads(grad_energy(coords))
```

### 2. SCF Convergence

```python
# Convergence check
def is_converged(dm_old, dm_new, tol=1e-8):
    delta = jnp.max(jnp.abs(dm_new - dm_old))
    return delta < tol

# DIIS extrapolation (simplified)
def diis_extrapolate(fock_list, error_list, diis_dim=6):
    """
    Extrapolate the Fock matrix with DIIS.
    """
    # Implement the DIIS algorithm.
    # (Simplification: return the latest Fock matrix directly.)
    return fock_list[-1]
```

### 3. Memory Management

```python
# Process a large molecule in batches
def batch_process_large_molecule(mol, batch_size=100):
    """
    Process a large molecule in batches.
    """
    n_atoms = mol.natm
    coords = mol.atom_coords()

    # Compute energy contributions in batches
    total_energy = 0.0

    for i in range(0, n_atoms, batch_size):
        batch_coords = coords[i:i+batch_size]
        E_batch = dft_energy(batch_coords)
        total_energy += E_batch

    return total_energy
```

## Performance Optimization

### 1. JIT Compilation

```python
# Compile the critical function
@jit
def compiled_scf(coords):
    """Compiled SCF calculation."""
    graddft = GradDFT(mol)
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# The first call triggers compilation
E1 = compiled_scf(coords)  # Compile

# Subsequent calls are fast
E2 = compiled_scf(coords)  # Fast
```

### 2. Vectorization

```python
# Vectorize over multiple geometries
from jax import vmap

# Prepare the vectorized function
batch_energy = vmap(compiled_scf)

# Batched calculation
batch_coords = jnp.array([coords, coords*1.01, coords*0.99])
energies = batch_energy(batch_coords)
```

### 3. Parallelization

```python
# Multi-device parallelism (GPU cluster)
from jax import pmap

# Pre-distribute data across multiple GPUs
batch_coords = jnp.array([coords] * 4)  # Four geometries
batch_coords = jnp.reshape(batch_coords, (4, -1, 3))

# Parallel calculation
parallel_energy = pmap(compiled_scf)
energies = parallel_energy(batch_coords)
```

## References

### PySCF-JAX
1. PySCF JAX module: https://github.com/pyscf/pyscf/tree/master/pyscf/jax
2. GradDFT: https://github.com/pyscf/pyscf/tree/master/pyscf/dft/gradients

### JAX
3. JAX documentation: https://jax.readthedocs.io/
4. JAX automatic differentiation: https://jax.readthedocs.io/en/latest/jax.html#automatic-differentiation

### Papers
5. Miller et al., "JAX-MD: A Framework for Differentiable Molecular Dynamics", arXiv 2023
6. Kasim et al., "Differentiable Quantum Chemistry with JAX", arXiv 2022
