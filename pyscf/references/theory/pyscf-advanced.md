# PySCF Advanced Features Reference

## Custom Exchange-Correlation Functionals

### 1. Basic Functional Structure

#### Functional Interface
```python
from pyscf.dft import libxc, numint

# Functional definition
def my_xc(rho, spin=0, deriv=1, **kwargs):
    """
    Custom exchange-correlation functional.

    Parameters:
        rho: (N, 4) or (N, 2) array
              spin=0: [rho, grad_x, grad_y, grad_z]
              spin=1: [rho_a, rho_b, grad_ax, grad_ay, grad_az,
                       grad_bx, grad_by, grad_bz]
        spin: Spin polarization (0: restricted, 1: unrestricted)
        deriv: Derivative order (1: energy and potential;
               2: derivatives of the potential, and so on)

    Returns:
        deriv=1: (ex+vc, vrho, vgamma)
            ex+vc: (N,) exchange-correlation energy density
            vrho: (N, 2) density potential (dE/drho_a, dE/drho_b)
            vgamma: (N, 3) gradient potential
                    (dE/dgamma_aa, dE/dgamma_ab, dE/dgamma_bb)

            where gamma = grad_rho·grad_rho
    """
    # Extract the density
    if spin == 0:
        r = rho[:, 0]
        grho2 = rho[:, 1]**2 + rho[:, 2]**2 + rho[:, 3]**2
    else:
        ra = rho[:, 0]
        rb = rho[:, 1]
        grho2_a = rho[:, 2]**2 + rho[:, 3]**2 + rho[:, 4]**2
        grho2_b = rho[:, 5]**2 + rho[:, 6]**2 + rho[:, 7]**2

    # Implement the functional here
    exc = ...  # Exchange-correlation energy density
    vrho = ...  # Density potential
    vgamma = ...  # Gradient potential

    return exc, vrho, vgamma
```

### 2. Implementing LDA Functionals

#### Slater Exchange (LDA)
```python
def lda_exchange(rho, spin=0, deriv=1):
    """
    Slater exchange functional.
    Ex = -3/4 * (3/π)^(1/3) * ∫ ρ^(4/3) dr
    """
    if spin == 0:
        r = rho[:, 0]
        # Exchange energy density
        exc = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        if deriv == 1:
            # Vx = dEx/drho = 4/3 * Ex/rho
            vrho = exc * (4/3) / (r + 1e-12)
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    else:
        ra = rho[:, 0]
        rb = rho[:, 1]
        # Spin-polarized LDA exchange
        exc = -0.75 * (3/np.pi)**(1/3) * (
            ra**(4/3) + rb**(4/3)
        )

        if deriv == 1:
            vrho_a = exc * (4/3) / (ra + 1e-12)
            vrho_b = exc * (4/3) / (rb + 1e-12)
            vrho = np.stack([vrho_a, vrho_b], axis=1)
            vgamma = np.zeros((len(rho), 3))
            return exc, vrho, vgamma

    return exc
```

#### VWN Correlation (LDA)
```python
def vwn_correlation(rho, spin=0, deriv=1):
    """
    Vosko-Wilk-Nusair (VWN) correlation functional.
    """
    # VWN parameters from the original paper
    A = 0.0310907
    x0 = -0.10498
    b = 3.72744
    c = 12.9352

    def vwn_func(rs):
        """
        rs = (3/(4πρ))^(1/3), the Wigner-Seitz radius.
        """
        x = np.sqrt(rs)
        X = x + b*x0 + c*x0**2
        Q = np.sqrt(4*c - b**2)

        fx = X**2 + Q**2
        fx0 = x0**2 + b*x0 + c*x0**2

        ec = A * (np.log(x**2/fx) +
                  2*b/Q * np.arctan((2*x + b)/Q) -
                  (b*x0)/Q * np.arctan((2*x0 + b)/Q) -
                  (x0*x)/fx)

        # Derivative
        dec_drs = -A/3 * (
            1/x - (2*x + b)/X +
            2*b/Q * (1/(2*x + b + Q*x) + 1/(2*x + b - Q*x))
        )

        return ec, dec_drs

    if spin == 0:
        r = rho[:, 0]
        rs = (3/(4*np.pi*r))**(1/3)
        ec, dec_drs = vwn_func(rs)

        exc = ec

        if deriv == 1:
            # Vc = dEc/drho = dEc/drs * drs/drho
            # drs/drho = -rs/(3*rho)
            drs_drho = -rs/(3*r)
            vrho = dec_drs * drs_drho
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    else:
        # Spin-polarized VWN (more involved)
        ra = rho[:, 0]
        rb = rho[:, 1]
        r = ra + rb
        zeta = (ra - rb) / r  # Spin-polarization parameter

        rs = (3/(4*np.pi*r))**(1/3)
        ec0, dec_drs = vwn_func(rs)  # Unpolarized value

        # Simplified spin-polarization correction
        alpha = 0.001  # Parameter
        exc = ec0 * (1 + alpha * f(zeta))

        if deriv == 1:
            # The complete spin-polarized derivative is more involved
            vrho_a = ...  # dE/drho_a
            vrho_b = ...  # dE/drho_b
            vrho = np.stack([vrho_a, vrho_b], axis=1)
            vgamma = np.zeros((len(rho), 3))
            return exc, vrho, vgamma

    return exc
```

### 3. Implementing GGA Functionals

#### Becke88 Exchange (B88)
```python
def becke88_exchange(rho, spin=0, deriv=1):
    """
    Becke 1988 GGA exchange functional.
    Ex = ∫ εx^LDA * Fx(s) dr
    where s = |∇ρ| / (2kF ρ) is the reduced density gradient.
    """
    beta = 0.0042  # Becke parameter

    if spin == 0:
        r = rho[:, 0]
        grx = rho[:, 1]
        gry = rho[:, 2]
        grz = rho[:, 3]

        grho2 = grx**2 + gry**2 + grz**2

        # Reduced density gradient
        kF = (3*np.pi**2*r)**(1/3)
        s = np.sqrt(grho2) / (2*kF*r + 1e-12)

        # LDA exchange
        ex_lda = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        # Becke enhancement factor
        s2 = s**2
        asinh_s = np.arcsinh(s)
        Fx = 1 + beta * s2 / (1 + 6*beta*s*asinh_s)

        # Exchange energy density
        exc = ex_lda * Fx

        if deriv == 1:
            # Computing the potential requires dFx/ds
            # Vx = dEx/drho - ∇·(dEx/d∇ρ)
            # This is a nontrivial functional-derivative calculation

            # Simplified version: return only the LDA potential
            vrho = ex_lda * (4/3) / (r + 1e-12) * Fx
            vgamma = np.zeros_like(rho[:, 0])

            # The complete GGA potential requires the tau term
            return exc, vrho, vgamma

    return exc
```

#### PBE Exchange
```python
def pbe_exchange(rho, spin=0, deriv=1):
    """
    Perdew-Burke-Ernzerhof exchange functional.
    """
    kappa = 0.804
    mu = 0.21951

    if spin == 0:
        r = rho[:, 0]
        grho2 = rho[:, 1]**2 + rho[:, 2]**2 + rho[:, 3]**2

        kF = (3*np.pi**2*r)**(1/3)
        s = np.sqrt(grho2) / (2*kF*r + 1e-12)

        # PBE enhancement factor
        Fx = 1 + kappa - kappa / (1 + mu*s**2/kappa)

        # LDA exchange
        ex_lda = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        exc = ex_lda * Fx

        if deriv == 1:
            # Simplified PBE potential
            vrho = ex_lda * (4/3) / (r + 1e-12) * Fx
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    return exc
```

### 4. Hybrid Functionals

#### Manual B3LYP Implementation
```python
def b3lyp_functional(mol):
    """
    Manual implementation of the B3LYP functional.
    Ex = 0.2*Ex^HF + 0.08*Ex^LDA + 0.72*Ex^B88
    Ec = 0.19*Ec^VWN + 0.81*Ec^LYP
    """
    from pyscf import dft

    mf = dft.RKS(mol)

    # Set the hybrid parameters
    mf.xc = {
        'hybrid': {
            'hf_fraction': 0.2,  # 20% HF exchange
        },
        'RXC': {
            'type': 'GGA',
            'exch': {
                'LDA': 0.08,
                'B88': 0.72,
            },
            'corr': {
                'VWN': 0.19,
                'LYP': 0.81,
            }
        }
    }

    # Or use the simpler form
    mf.xc = 'b3lyp'

    return mf
```

#### Custom Hybrid Functional
```python
def custom_hybrid(mol, hf_fraction=0.25, xc_dft='pbe'):
    """
    Custom hybrid functional.
    Ex = α*Ex^HF + (1-α)*Ex^DFT
    """
    mf = dft.RKS(mol)

    # Specify it with a string
    mf.xc = f'{hf_fraction}*HF + {1-hf_fraction}*{xc_dft}'

    # Or use the comma-separated form
    mf.xc = f'{hf_fraction}*HF,{xc_dft}'

    return mf
```

### 5. Applying a Custom Functional

```python
from pyscf import gto, dft, lib

# Define the molecule
mol = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# Method 1: assign the functional directly
mf = dft.RKS(mol)
mf._numint._xc = my_xc  # Replace the functional-evaluation function
e = mf.kernel()

# Method 2: register it through libxc
# Register the new functional
libxc.define_xc_ = libxc.define_xc_.copy()
libxc.define_xc_('MY_XC', my_xc)

# Use it
mf = dft.RKS(mol)
mf.xc = 'MY_XC'
e = mf.kernel()

# Method 3: use a NumInt object
from pyscf.dft import numint

ni = numint.NumInt()
ni._xc = my_xc

mf = dft.RKS(mol)
mf._numint = ni
e = mf.kernel()
```

## Integral Operations

### 1. Two-Electron Integrals

#### Raw Two-Electron Integrals (8-Fold Symmetry)
```python
from pyscf import gto

mol = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# Full two-electron integrals (8-fold symmetry)
# (ij|kl) ~ (ji|kl) ~ (ij|lk) ~ ...
eri = mol.intor('int2e')
print(f"Integral shape: {eri.shape}")  # (nao, nao, nao, nao)

# A specific integral
i, j, k, l = 0, 0, 0, 0
value = eri[i, j, k, l]
print(f"({i}{j}|{k}{l}) = {value:.6f}")
```

#### 4-Fold-Symmetric Integrals
```python
# 4-fold symmetry: (ij|kl) = (kl|ij)
eri_4 = mol.intor('int2e_sph', aosym=4)
print(f"4-fold-symmetric shape: {eri_4.shape}")  # (nao*(nao+1)/2, nao*(nao+1)/2)
```

#### Density-Fitting Integrals
```python
from pyscf import df

# Create the DF object
dfobj = df.DF(mol)
dfobj.auxbasis = 'cc-pvdz-ri'

# Compute 3-center integrals (ij|A)
ijA = dfobj.get_2c2e()  # (naux, nao, nao)

# Compute 2-center integrals (A|B)
AB = dfobj.get_2c2e_aux()  # (naux, naux)

# Fit the two-electron integrals
eri_df = dfobj.get_j()  # Approximate two-electron integrals
```

### 2. AO-to-MO Integral Transformation

#### Full Transformation (Memory Intensive)
```python
from pyscf import ao2mo

# Two-electron integrals in the AO basis
eri_ao = mol.intor('int2e')

# MO coefficients
C = mf.mo_coeff  # (nao, nmo)

# Full transformation to the MO basis
eri_mo = ao2mo.incore.full(eri_ao, C)
print(f"MO integral shape: {eri_mo.shape}")  # (nmo, nmo, nmo, nmo)

# A specific MO integral
i, j, k, l = 0, 1, 2, 3  # Orbital indices
value = eri_mo[i, j, k, l]
print(f"({i}{j}|{k}{l})_MO = {value:.6f}")
```

#### Partial Transformation (Occupied-Occupied)
```python
# Occupied orbitals
nocc = mol.nelectron // 2
C_occ = C[:, :nocc]  # (nao, nocc)

# Compute only occupied-occupied integrals
eri_oooo = ao2mo.incore.general(eri_ao,
                                [C_occ, C_occ, C_occ, C_occ])
print(f"(oo|oo) shape: {eri_oooo.shape}")  # (nocc*nocc, nocc*nocc)
```

#### Out-of-Core Transformation (Large Molecules)
```python
# Save to disk
ao2mo.outcore.full(mol, C, 'eri_mo.h5')

# Read a specific block
from pyscf import lib
with lib.H5File('eri_mo.h5', 'r') as f:
    eri_block = f['eri_mo'][:10, :10, :10, :10]

# Or use memory mapping
eri_mo = ao2mo.load('eri_mo.h5')
```

### 3. Integral Derivatives

#### Basis-Function Derivative Integrals
```python
# Basis-function gradient: d/dx χ_a(r)
# dS_ij/dx_A = ∫ [dχ_i/dx_A χ_j + χ_i dχ_j/dx_A] dr

# Overlap-matrix derivative
dS_dx = mol.intor('int1e_ovlp_ip1', comp=3)  # (3, nao, nao)
# dS_dx[0] = d/dx, dS_dx[1] = d/dy, dS_dx[2] = d/dz

# Kinetic-energy matrix derivative
dT_dx = mol.intor('int1e_kin_ip1', comp=3)

# Nuclear-attraction matrix derivative
dV_dx = mol.intor('int1e_nuc_ip1', comp=3)
```

#### Two-Electron Integral Derivatives
```python
# Two-electron integral derivative (3-center)
# d(ij|kl)/dA = d/dx_A (ij|kl)
eri_ip1 = mol.intor('int2e_ip1', comp=3)
print(f"Integral-derivative shape: {eri_ip1.shape}")  # (3, nao, nao, nao, nao)
```

### 4. Relativistic Integrals

#### DKH (Douglas-Kroll-Hess)
```python
# Douglas-Kroll-Hess Hamiltonian
from pyscf.dkh import dkh

# DKH2
dkh2 = dkh.make_dkh2_hamiltonian(mol)

# DKH-corrected kinetic and potential energies
T_dkh, V_dkh = dkh2

# Use in SCF
mf = scf.RHF(mol)
mf.get_hcore = lambda *args: T_dkh + V_dkh
mf.kernel()
```

#### ZORA (Zeroth-Order Regular Approximation)
```python
from pyscf.sfx2c1e import sfx2c1e

# Spin-free X2C1e
mf = scf.RHF(mol)
mf = sfx2c1e.get_x2c(mf)  # Convert to relativistic SCF
mf.kernel()
```

## Advanced SCF Control

### 1. DIIS Acceleration

#### Custom DIIS
```python
from pyscf.scf import diis

# DIIS object
my_diis = diis.DIIS()

# SCF cycle
for cycle in range(max_cycle):
    # Build the Fock matrix
    fock = build_fock(dm)

    # DIIS extrapolation
    if cycle >= diis_start:
        fock = my_diis.update(fock, x=dm)

    # Diagonalize
    mo_energy, mo_coeff = eig(fock, S)

    # New density
    dm_new = build_density(mo_coeff)

    # Check convergence
    if converged(dm, dm_new):
        break

    dm = dm_new
```

#### DIIS Error Vector
```python
# Custom error vector
def error_vector(fock, dm, S):
    """
    DIIS error: [F,D]S = FDS - SDF
    """
    return fock @ dm @ S - S @ dm @ fock

# Use it in DIIS
my_diis = diis.DIIS()
error = error_vector(fock, dm, S)
fock = my_diis.update(fock, x=dm, errvec=error)
```

### 2. Direct SCF

#### Controlling Direct Evaluation
```python
mf = scf.RHF(mol)

# Fully direct: recompute integrals in every cycle
mf.direct_scf = True

# Semidirect: cache selected integrals
mf.direct_scf = False
mf._eri = mol.intor('int2e')  # Precompute all integrals

# Hybrid mode
mf.direct_scf = True
mf.direct_scf_tol = 1e-12  # Integral threshold
```

#### Memory Management
```python
# Set the thread count
lib.num_threads(4)  # 4 threads
lib.chkfile.save_temp(mol.chkfile, mf.get_fock())  # Save intermediate results

# Out-of-core SCF
mf.chkfile = 'calc.chk'
mf.kernel()
# Intermediate results are saved to the checkpoint file
```

### 3. SCF Stabilization

#### Stability Analysis
```python
# Check the stability of the SCF solution
mf = scf.RHF(mol)
mf.kernel()

# Stability analysis
from pyscf.scf import stability

# Internal stability (possible UHF instability)
stable = stability.stability_internal(mf)
print(f"Internally stable: {stable}")

# External stability (is there a lower-energy solution?)
stable2, e_new, mf_new = stability.stability_external(mf)
print(f"Externally stable: {stable2}")

# If unstable, use the new SCF solution
if not stable2:
    mf_new.kernel()
```

#### Damping and Level Shifting
```python
# Damping for oscillatory convergence
mf = scf.RHF(mol)
mf.damp = 0.2  # 20% damping
mf.kernel()

# Level shifting for slow convergence
mf = scf.RHF(mol)
mf.level_shift = 0.5  # Shift virtual orbitals by 0.5 Ha
mf.kernel()

# Combine both methods
mf.damp = 0.2
mf.level_shift = 0.5
mf.diis_start_cycle = 5  # Delay DIIS
```

### 4. Newton-Raphson SCF

#### Second-Order SCF
```python
# Newton-Raphson method (fast convergence, but potentially unstable)
mf = scf.RHF(mol)
mf_nr = mf.newton()  # Create the NR-SCF object
mf_nr.kernel()

# Constrained NR-SCF
mf_nr.constrain_occ = True  # Preserve occupations
mf_nr.kernel()
```

#### Odd-Electron Systems
```python
# ROHF stability
mol = gto.M(atom='C 0 0 0', charge=0, spin=2, basis='sto-3g')
mf = scf.ROHF(mol)
mf.kernel()

# Convert to UHF and stabilize
mf_uhf = scf.UHF(mol)
mf_uhf.kernel()

# Stability check
mf_uhf_stable = stability.stability_internal(mf_uhf)
```

## Advanced TDDFT Features

### 1. Natural Transition Orbitals (NTOs)

#### Detailed NTO Analysis
```python
from pyscf import tdscf

td = tdscf.TDDFT(mf)
td.nstates = 5
td.kernel()

# NTOs for all states
nocc = np.count_nonzero(mf.mo_occ > 0)
for n in range(td.nstates):
    weights, nto_coeff = td.get_nto(state=n + 1)

    print(f"\nState {n+1}:")
    print(f"  Excitation energy: {td.e[n]*27.2114:.2f} eV")
    print(f"  Leading weight: {weights.max():.4f}")

    # PySCF returns occupied NTOs followed by virtual NTOs.
    occupied_ntos = nto_coeff[:, :nocc]
    virtual_ntos = nto_coeff[:, nocc:]

    # Dominant NTO pair
    dominant_pair = np.argmax(weights)
    print(f"  Dominant-pair weight: {weights[dominant_pair]:.4f}")

    # Coefficient vectors for a downstream, separately validated population
    # or fragment-projection analysis.
    hole_orb = occupied_ntos[:, dominant_pair]  # (nao,)
    electron_orb = virtual_ntos[:, dominant_pair]
```

#### NTO Visualization
```python
# Save NTOs to Molden files
from pyscf.tools import molden

nocc = np.count_nonzero(mf.mo_occ > 0)
for n in range(td.nstates):
    weights, nto_coeff = td.get_nto(state=n + 1)
    occupied_ntos = nto_coeff[:, :nocc]
    virtual_ntos = nto_coeff[:, nocc:]

    # Dominant hole orbital
    hole = occupied_ntos[:, np.argmax(weights)]
    # Dominant electron orbital
    electron = virtual_ntos[:, np.argmax(weights)]

    # Save
    molden.from_mo(mol, f'hole_state{n+1}.molden', hole)
    molden.from_mo(mol, f'electron_state{n+1}.molden', electron)
```

### 2. TDDFT Density Analysis

The core TD-SCF interface exposes excitation energies, response amplitudes,
transition moments, oscillator strengths, gradients, and NTOs. It does not
provide one method-independent excited-state density interface suitable for
the executable example previously shown here. Use a density implementation
documented for the selected PySCF method and version, then validate the density
trace, spin convention, AO ordering, coordinate units, and state identity
before deriving atomic charge differences or a charge-transfer distance.

Do not substitute AO populations for atomic charges. A charge-transfer metric
also needs an explicitly defined scalar reduction of the charge-displacement
vector and an explicit Bohr-to-angstrom conversion when coordinates are in
atomic units. See the [PySCF TD-SCF documentation](https://pyscf.org/user/tddft.html)
for the currently supported public properties.

### 3. Spin-Orbit Coupling (SOC)

The repository does not provide a validated, general TDDFT spin-orbit coupling
implementation. Do not infer a coupling matrix from an undocumented convenience
module or insert a placeholder coupling into a rate calculation. Select and
document a supported relativistic/SOC workflow, state convention, and unit
conversion, and retain its source output as provenance before using the result.

## MP2 and Post-HF Methods

### 1. MP2

#### RMP2 Calculation
```python
from pyscf import mp

# MP2 object
mp2 = mp.MP2(mf)
e_mp2, t2 = mp2.kernel()

print(f"MP2 correlation energy: {e_mp2:.6f} Ha")
print(f"Total energy (HF+MP2): {mf.e_tot + e_mp2:.6f} Ha")

# MP2 energy decomposition
# E(MP2) = E(HF) + E_2
# E_2 = Σ_{ijab} |<ij||ab>|^2 / (ε_i+ε_j-ε_a-ε_b)

# Orbital energies
eps_i = mf.mo_energy[:nocc]  # Occupied
eps_a = mf.mo_energy[nocc:]  # Virtual

# Vertex: <ij||ab> = <ij|ab> - <ij|ba>
# Stored in t2

# Diagnostic parameters
# D1 = Σ_{ijab} |t2_{ijab}|^2
# T1 = Σ_{ia} |t1_{ia}|  (for MP3+)

d1 = mp2.get_d1()
t1 = mp2.get_t1()

print(f"D1 diagnostic: {d1:.3f} (ideal < 0.02)")
print(f"T1 diagnostic: {t1:.3f} (ideal < 0.02)")

# Multireference warning
if d1 > 0.02 or t1 > 0.02:
    print("Warning: Significant multireference character; consider CASSCF")
```

#### Local MP2 (LMP2)
```python
# Density-fitted MP2
mp2_df = mp.DFMP2(mf)
mp2_df.auxbasis = 'cc-pvtz-ri'
e_mp2_df = mp2_df.kernel()[0]

print(f"DF-MP2 correlation energy: {e_mp2_df:.6f} Ha")
print(f"Difference from conventional MP2: {abs(e_mp2 - e_mp2_df):.6f} Ha")
```

#### RI-MP2 (Density Fitting)
```python
# Use the resolution of the identity
mp2_ri = mp.MP2(mf)
mp2_ri.with_df = df.DF(mol)
mp2_ri.with_df.auxbasis = 'cc-pvdz-ri'
e_mp2_ri = mp2_ri.kernel()[0]

print(f"RI-MP2 correlation energy: {e_mp2_ri:.6f} Ha")
```

### 2. Coupled Cluster

#### CCSD
```python
from pyscf import cc

# CCSD object
ccsd = cc.CCSD(mf)
e_ccsd, t1, t2 = ccsd.kernel()

print(f"CCSD correlation energy: {e_ccsd:.6f} Ha")
print(f"Total energy: {mf.e_tot + e_ccsd:.6f} Ha")

# T1 diagnostic
t1_norm = np.linalg.norm(t1)
print(f"T1 diagnostic: {t1_norm:.3f}")

# CCSD(T), the gold standard
e_ccsdt = ccsd.ccsd_t()
print(f"(T) correction: {e_ccsdt:.6f} Ha")
print(f"CCSD(T) total energy: {mf.e_tot + e_ccsd + e_ccsdt:.6f} Ha")
```

#### Decomposition
```python
# E(CCSD) = Σ_{ia} t_{ia} <ia||ia> + Σ_{ijab} t_{ijab} <ij||ab>

# Component energies
e_1 = ccsd.energy1()  # Singles
e_2 = ccsd.energy2()  # Doubles

print(f"E(singles): {e_1:.6f} Ha")
print(f"E(doubles): {e_2:.6f} Ha")
```

#### Lambda Equations (Excited States)
```python
# CCSD Lambda equations for excited states
ccsd_lambda = cc.ccsd_lambda.CCSD_Lambda(ccsd)

# Ground-state Lambda amplitudes
l1, l2 = ccsd_lambda.kernel()

# EOM-CCSD (equation-of-motion coupled cluster)
from pyscf import cc

eom = cc.eom_rccsd.EOMIP(ccsd)  # Ionization
# eom = cc.eom_rccsd.EOMEA(ccsd)  # Electron attachment
# eom = cc.eom_rccsd.EOMEES(ccsd)  # Excited states

ip_e, ip_vec = eom.ipccsd(nroots=3)  # 3 ionized states
print("Vertical ionization energies (eV):")
for i, e in enumerate(ip_e):
    print(f"  State {i+1}: {e*27.2114:.2f} eV")
```

### 3. CI Methods

#### CIS (Configuration Interaction Singles)
```python
from pyscf import ci

# CIS object
cis = ci.CIS(mf)
cis.nstates = 5
e_cis, civec = cis.kernel()

print("CIS excited states:")
for i in range(cis.nstates):
    print(f"  State {i+1}: {e_cis[i]*27.2114:.2f} eV")
```

#### CISD
```python
# CISD (CIS plus double excitations)
cisd = ci.CISD(mf)
e_cisd, ci_vec = cisd.kernel()

print(f"CISD correlation energy: {e_cisd - mf.e_tot:.6f} Ha")

# Davidson algorithm
cisd.max_space = 100  # Davidson subspace size
e_cisd, ci_vec = cisd.kernel()
```

### 4. FCI (Full CI)

#### FCI for a Small System
```python
from pyscf import fci

# Small molecule
mol_small = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)
mf_small = scf.RHF(mol_small)
mf_small.kernel()

# FCI calculation
cisolver = fci.FCI(mf_small)
e_fci, ci_vec = cisolver.kernel()

print(f"FCI energy: {e_fci:.6f} Ha")
print(f"FCI vs HF: {(e_fci - mf_small.e_tot)*27.2114:.2f} eV")

# Exact wavefunction analysis
print(f"Number of CI coefficients: {len(ci_vec)}")
print(f"Leading-configuration weight: {abs(ci_vec).max():.4f}")
```

#### Selected CI
```python
# Selected CI in a finite space
fcisolver = fci.FCI(mf_small)
fcisolver.nroots = 3  # Multiple states

# Selected space
norb = mol_small.nao
nelec = mol_small.nelectron

# Select the orbital space manually
# This example uses full CI
e_fci, civec = fcisolver.kernel()
```

## Density Fitting (DF)

### 1. Basic DF

#### DF-SCF
```python
from pyscf import df

# Automatic DF
mf = scf.RHF(mol)
mf_df = mf.density_fit()

# Or specify the fitting basis manually
mf_df = df.density_fit(mf, auxbasis='def2-universal-jfit')

e_df = mf_df.kernel()
print(f"DF-RHF energy: {e_df:.6f} Ha")

# Compare with conventional RHF
print(f"DF vs RHF energy difference: {abs(e_df - mf.kernel()):.8f} Ha")
```

#### DF Object
```python
# DF object
dfobj = df.DF(mol)
dfobj.auxbasis = 'cc-pvdz-ri'

# 3-center integrals (ij|P)
ijP = dfobj.get_2c2e()
print(f"(ij|P) shape: {ijP.shape}")  # (naux, nao, nao)

# 2-center integrals (P|Q)
PQ = dfobj.get_2c2e_aux()
print(f"(P|Q) shape: {PQ.shape}")  # (naux, naux)

# Solve (P|Q) * J^P = (ij|P)
# J = Σ_PQ (ij|P) * (P|Q)^(-1)
```

### 2. MP2-DF

#### DF-MP2
```python
from pyscf import mp

# MP2 plus DF
mp2 = mp.MP2(mf)
mp2_df = mp.density_fit(mp2, auxbasis='cc-pvdz-ri')

e_mp2_df = mp2_df.kernel()[0]
print(f"DF-MP2 correlation energy: {e_mp2_df:.6f} Ha")

# Check the difference from conventional MP2
e_mp2_std = mp.MP2(mf).kernel()[0]
print(f"DF vs conventional MP2: {abs(e_mp2_df - e_mp2_std):.6f} Ha")
```

### 3. CCSD-DF

#### DF-CCSD
```python
# CCSD plus DF
ccsd = cc.CCSD(mf)
ccsd_df = ccsd.density_fit()

e_ccsd_df = ccsd_df.kernel()[0]
print(f"DF-CCSD correlation energy: {e_ccsd_df:.6f} Ha")
```

## Periodic Systems

### 1. Defining Periodic Systems

#### 1D Periodicity
```python
from pyscf.pbc import gto as pbcgto

# 1D cell (linear chain)
cell = pbcgto.M(
    atom='H 0 0 0',
    basis='sto-3g',
    a=[[2.0, 0, 0], [0, 20, 0], [0, 0, 20]],  # Lattice vectors
    unit='Ang',
)

# k-point sampling
kpts = cell.make_kpts([2, 1, 1])  # 2 k-points
```

#### 2D Periodicity
```python
# 2D cell (graphene layer)
cell = pbcgto.M(
    atom='''
    C  0.0000  0.0000  0.0000
    C  0.0000  1.4200  0.0000
    C  1.2308  0.7100  0.0000
    C  1.2308  2.1300  0.0000
    ''',
    basis='gth-szv',
    a=[[2.4616, 0, 0], [0, 4.2600, 0], [0, 0, 20]],  # In-plane periodicity
    pseudo='gth-pade',
    ke_cutoff=100,  # Kinetic-energy cutoff (eV)
)

# k-points
kpts = cell.make_kpts([4, 4, 1])
```

#### 3D Periodicity
```python
# 3D cell (silicon)
cell = pbcgto.M(
    atom='Si 0 0 0; Si 0.25 0.25 0.25',
    basis='gth-szv',
    a=[[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
    pseudo='gth-pade',
    ke_cutoff=100,
)

# k-point mesh
kpts = cell.make_kpts([4, 4, 4])
```

### 2. Periodic SCF

#### Gamma-only
```python
from pyscf.pbc import scf as pbcscf

# Gamma point only
mf = pbcscf.RHF(cell)
e_gam = mf.kernel()
print(f"Gamma-only energy: {e_gam:.6f} Ha")
```

#### Multiple-k-Point SCF
```python
# Multiple k-points
mf = pbcscf.RHF(cell, kpts)
e_kpts = mf.kernel()
print(f"Multiple-k-point energy: {e_kpts:.6f} Ha")

# k-point weights
print(f"Number of k-points: {len(kpts)}")
print(f"Total weight: {mf.kpts.sum()}")
```

### 3. Periodic DFT

#### PBC-DFT
```python
from pyscf.pbc import dft as pbcdft

# Periodic DFT
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'pbe'
mf.grids.level = 3  # Grid accuracy
e_pbe = mf.kernel()
print(f"PBE energy: {e_pbe:.6f} Ha")
```

#### Hybrid Functionals
```python
# Periodic hybrid functional (computationally expensive)
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'hse06'  # Screened hybrid
e_hse = mf.kernel()
print(f"HSE06 energy: {e_hse:.6f} Ha")

# Or use density fitting for exact exchange
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'pbe0'
mf = mf.density_fit()  # Accelerate HF exchange
e_pbe0 = mf.kernel()
```

### 4. Band Structure

#### Band-Structure Calculation
```python
# Compute the band path
from pyscf.pbc import tools

# High-symmetry-point path
Gamma = [0, 0, 0]
X = [0.5, 0, 0]
M = [0.5, 0.5, 0]

# Generate the k-path
kpath = tools.get_bandpath([Gamma, X, M], cell, 20)

# Compute bands at every k-point
band_energies = []
for k in kpath:
    mf_k = pbcdft.RKS(cell, k.reshape(1, 3))
    mf_k.xc = 'pbe'
    mf_k.kernel()
    band_energies.append(mf_k.mo_energy[0])

band_energies = np.array(band_energies)  # (nk, nmo)

# Plot the bands
import matplotlib.pyplot as plt
plt.plot(band_energies)
plt.xlabel('k-path')
plt.ylabel('Energy (Ha)')
plt.title('Band Structure')
plt.show()
```

## Solvent Effects

### 1. PCM (Polarizable Continuum Model)

#### PCM-SCF
```python
from pyscf import solvent

# PCM model
pcm = solvent.PCM(mol)
pcm.eps = 78.4  # Dielectric constant of water
pcm.method = 'IEFPCM'  # Integral-equation-formalism PCM

# SCF with PCM
mf = scf.RHF(mol)
mf = pcm.run(mf)
e_pcm = mf.e_tot
print(f"PCM energy: {e_pcm:.6f} Ha")

# Solvation energy
e_gas = scf.RHF(mol).kernel()
G_sol = e_pcm - e_gas
print(f"Solvation energy: {G_sol*27.2114:.2f} eV")
```

#### Different Solvents
```python
# Different solvents
solvents = {
    'water': 78.4,
    'ethanol': 24.3,
    'methanol': 32.6,
    'acetone': 20.7,
    'dichloromethane': 8.93,
}

for name, eps in solvents.items():
    pcm = solvent.PCM(mol)
    pcm.eps = eps
    mf_pcm = scf.RHF(mol)
    mf_pcm = pcm.run(mf_pcm)
    print(f"{name:15s}: {mf_pcm.e_tot:.6f} Ha")
```

### 2. SMD Model

#### SMD-SCF
```python
# SMD (Solvation Model based on Density)
pcm = solvent.PCM(mol)
pcm.method = 'SMD'
pcm.solvent = 'water'  # Solvent name

mf = scf.RHF(mol)
mf = pcm.run(mf)
print(f"SMD energy: {mf.e_tot:.6f} Ha")
```

### 3. Explicit Solvent

#### Microcluster Model
```python
# Microcluster consisting of water plus 6 water molecules
cluster = gto.M(
    atom='''
    O  0.000000  0.000000  0.000000
    H  0.000000  0.758602 -0.504284
    H  0.000000 -0.758602 -0.504284
    O  2.800000  0.000000  1.500000
    H  3.558602  0.000000  1.995716
    H  2.041398  0.000000  1.995716
    ... (additional water molecules)
    ''',
    basis='cc-pvdz'
)

mf_cluster = scf.RHF(cluster)
e_cluster = mf_cluster.kernel()
print(f"Cluster energy: {e_cluster:.6f} Ha")
```

## Parallel Computation

### 1. OpenMP Parallelism

#### Setting the Thread Count
```python
import os

# Set the OpenMP thread count
os.environ['OMP_NUM_THREADS'] = '8'

# Or configure it through PySCF
lib.num_threads(8)

# SCF parallelizes automatically
mf = scf.RHF(mol)
e = mf.kernel()
print(f"Number of threads used: {lib.num_threads()}")
```

#### Performance Analysis
```python
# Enable timing
lib.logger.TIMER_LEVEL = 3

mf = scf.RHF(mol)
e = mf.kernel()

# Inspect timing information
# Timing information is written to the log
```

### 2. MPI Parallelism

#### MPI-SCF
```python
# Launch with mpi4py
# mpirun -np 4 python script.py

from pyscf import lib

# Check MPI
if hasattr(lib, 'mpi'):
    print(f"MPI process: {lib.mpi.rank} of {lib.mpi.size}")
else:
    print("MPI is not enabled")

# SCF parallelizes automatically
mf = scf.RHF(mol)
e = mf.kernel()
```

## References and Resources

### PySCF Documentation
1. Official PySCF documentation: https://pyscf.org/
2. PySCF examples: https://github.com/pyscf/pyscf/tree/master/examples
3. PySCF API documentation: https://pyscf.org/api.html

### Papers
4. Sun et al., "Recent developments in the PySCF program package", WIREs Comput Mol Sci, 2020
5. Sun et al., "PySCF: the Python-based simulations of chemistry framework", WIREs Comput Mol Sci, 2018

### Tutorials and Courses
6. PySCF workshop: https://pyscf.org/workshop/
7. Quantum chemistry with PySCF: https://pyscf.org/tutorial.html
