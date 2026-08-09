---
name: pyscf
description: Use when building, reviewing, or troubleshooting PySCF calculations for HF, DFT, correlated wavefunction methods, excited states, geometry optimization, PES scans, or spectroscopy.
license: MIT
compatibility: Requires Python and PySCF; optional workflows require additional optimizer, plotting, or accelerator packages.
---

# PySCF

## Overview

Use this skill to turn a quantum-chemistry question into a traceable PySCF
calculation. The installed PySCF version and its documentation are authoritative;
examples here are starting points, not evidence that a method is appropriate for
a particular system.

**Core rule:** define the physical quantity first, run the smallest defensible
calculation, and reject unconverged or ambiguously defined results.

Install the selected environment explicitly:

```bash
python -m pip install pyscf
python -c "import pyscf; print(pyscf.__version__)"
```

## When to Use

Use this skill for:

- molecular HF or Kohn-Sham DFT;
- MP2, coupled-cluster, CASCI, CASSCF, or related post-SCF work;
- TDA, TDHF, or TDDFT excited states;
- gradients, geometry optimization, Hessians, or PES scans;
- periodic PySCF workflows;
- review of PySCF inputs, outputs, units, convergence, and provenance.

Do not use a molecular PySCF template unchanged for periodic systems, unsupported
relativistic Hamiltonians, multireference states, or production benchmarks.
Route those cases to the relevant specialized API and validate the method first.

## Non-Negotiable Scientific Boundaries

Before reporting a result:

1. Record geometry source and coordinate unit, charge, and `spin = Nalpha - Nbeta = 2S`.
2. Record state/root, method, functional, basis/ECP, numerical grid, solvent model,
   thresholds, PySCF version, and relevant optional backends.
3. Check SCF and every downstream convergence flag that the selected method exposes.
4. Keep total energies, excitation energies, gaps, and correlation energies distinct.
5. Track excited-state character during scans or optimization; a root number alone
   does not prove state identity.
6. Treat density fitting, finite grids, frozen cores, ECPs, and truncated active
   spaces as declared approximations.
7. Never label an optimization as a minimum or transition state without an
   appropriate stationary-point and vibrational analysis.

If any required definition or convergence evidence is missing, return the result
as incomplete rather than filling in a plausible value.

## Core Workflow

### 1. Define the target

Specify the quantity and reference state before choosing code:

- electronic total energy or relative energy;
- vertical or adiabatic excitation/emission energy;
- optimized geometry or stationary-point class;
- response property, spectrum, population, or orbital diagnostic;
- rigid or relaxed PES, including which coordinates are constrained.

### 2. Build and validate the molecule

```python
from pyscf import gto

mol = gto.M(
    atom="H 0.0 0.0 0.0; F 0.0 0.0 0.92",
    unit="Angstrom",
    basis="def2-svp",
    charge=0,
    spin=0,
    symmetry=False,
    verbose=4,
)

assert mol.nelectron > 0
assert mol.nelectron % 2 == mol.spin % 2
```

Use an explicit coordinate unit. Confirm the electron count, charge, spin, atom
order, basis availability, and ECP assignment before launching an expensive job.

### 3. Establish a converged reference

```python
from pyscf import dft

mf = dft.RKS(mol)
mf.xc = "pbe0"
mf.grids.level = 4
mf.conv_tol = 1e-9
mf.max_cycle = 100
energy_hartree = mf.kernel()

if not mf.converged:
    raise RuntimeError("SCF did not converge; do not use downstream properties")

print(f"E = {energy_hartree:.12f} Eh")
```

Use RHF/RKS only for an appropriate closed-shell reference. For an open-shell
system, choose ROHF/UHF or ROKS/UKS deliberately and inspect spin contamination,
orbital occupations, and possible broken-symmetry solutions.

### 4. Add one downstream layer at a time

Choose the least complex method that answers the stated question:

| Question | Starting family | Required checks |
|---|---|---|
| Closed-shell mean-field baseline | RHF or RKS | SCF convergence, stability, basis/grid convergence |
| Open-shell baseline | ROHF/UHF or ROKS/UKS | `<S^2>`, occupations, stability |
| Dynamic correlation | MP2 or CC methods | reference quality, frozen-core choice, basis convergence |
| Strong static correlation | CASCI/CASSCF | active-space definition, roots, occupations |
| Single-excitation-dominated states | TDA/TDHF/TDDFT | ground-state quality, root convergence and character |
| Periodic material | `pyscf.pbc` | cell, k mesh, Coulomb treatment, finite-size convergence |

Do not infer that a functional, basis, or method is accurate because it is listed.
Benchmark or justify it for the target property and chemical regime.

### 5. Verify, then report

Perform a cheaper smoke test first, followed by systematic convergence or method
checks that can change the conclusion. Save inputs, output, checkpoints, software
versions, and a compact table of assumptions and results.

See [workflow and safety details](references/theory/workflow-and-safety.md) for the
full decision checklist.

## Excited States

Start from a converged ground-state reference. TDA is often a useful diagnostic,
but it is not universally more accurate than full TDDFT.

```python
from pyscf import lib

td = mf.TDA()
td.nstates = 5
td.kernel()

if not all(td.converged):
    raise RuntimeError("One or more requested excited states did not converge")

excitation_ev = td.e * lib.param.HARTREE2EV
oscillator_strength = td.oscillator_strength()
nto_weights, nto_coeff = td.get_nto(state=1)  # one-based state index

s1_total_energy_hartree = mf.e_tot + td.e[0]
```

Keep these quantities separate:

- `mf.e_tot`: ground-state total energy in Hartree;
- `td.e[i]`: excitation energy relative to that reference in Hartree;
- `mf.e_tot + td.e[i]`: corresponding excited-state total-energy estimate;
- `td.oscillator_strength()`: dimensionless oscillator strength.

For spectra, document line shape, width, energy grid, temperature or environment,
and whether the sticks are absorption or emission transitions. Do not call a
ground-geometry excitation an emission calculation.

## Geometry, Frequencies, and PES

- Confirm that analytic gradients/Hessians exist for the exact method and state.
- Geometry optimizers are optional integrations; install and record the selected
  optimizer separately.
- Check both electronic convergence and optimizer termination at every geometry.
- Convert a mass-weighted Hessian to frequencies with the correct units; raw
  Hessian eigenvalues are not frequencies in `cm^-1`.
- For an excited-state optimization, use a state-specific gradient scanner and
  monitor state character for root flips or crossings.
- For a PES, state whether each point is rigid or relaxed, preserve constraints,
  and compute the full energy at every geometry.
- An S1 PES uses `E_S1(R) = E_S0(R) + omega_S1(R)`, not `omega_S1(R)` alone.

Use the dedicated [2D PES guide](references/practice/2d-potential-energy-surface.md)
and [emission workflow](references/practice/emission-spectrum-workflow.md) before
adapting those calculations.

## Units and Numerical Meaning

| Quantity | Typical PySCF value | Reporting requirement |
|---|---|---|
| Electronic energy | Hartree (`Eh`) | State the conversion used for eV or kJ/mol |
| TD excitation energy | Hartree | Convert once; do not mix with total energy |
| Nuclear gradient | `Eh/Bohr` | Match optimizer and coordinate conventions |
| Input coordinates | User-selected | Set `unit` explicitly |
| Electric dipole | Atomic units | State conversion if reporting Debye |
| Oscillator strength | Dimensionless | Keep state index and gauge/protocol |

Useful constants are exposed by `pyscf.lib.param`; avoid copying rounded constants
through multiple post-processing layers.

## Convergence and Failure Handling

| Symptom | Safe response |
|---|---|
| SCF oscillates or stalls | Recheck charge/spin/geometry; inspect occupations; then test damping, level shift, a new guess, or Newton SCF |
| Different guesses give different energies | Run stability analysis and characterize competing solutions |
| Open-shell result has unexpected `<S^2>` | Inspect spin contamination and compare appropriate restricted/open-shell references |
| TD roots reorder | Track transition character/NTOs, not only root number |
| Optimization stops unexpectedly | Inspect electronic convergence, gradients, constraints, and optimizer status |
| Memory is exhausted | Estimate tensor sizes; consider density fitting, direct algorithms, frozen cores, or smaller tests |
| Result changes with grid/basis | Report the convergence study; do not select the preferred value silently |

Changing the Hamiltonian or electronic-structure method is a scientific decision,
not merely a convergence fix. Keep failed attempts and their diagnostics in the
calculation record.

## Performance and Restart

- Set thread counts before expensive work and record them with hardware/backend.
- Use `max_memory` and estimate AO/MO integral storage before post-HF jobs.
- Density fitting can reduce cost, but validate its error for the target observable.
- Save checkpoints for restart; verify that the loaded molecule, orbitals, method,
  and basis match the intended calculation.
- Benchmark with a fixed geometry, method, basis, thread count, and warm-up policy.

## Included Tools

| Script | Intended starting point |
|---|---|
| `tools/scf.py` | HF/DFT reference calculations |
| `tools/dft.py` | DFT setup and comparison |
| `tools/tddft.py` | TDA/TDDFT calculations |
| `tools/mp2.py` | MP2 calculations |
| `tools/ccsd.py` | CCSD and perturbative triples |
| `tools/cascf.py` | Active-space workflows |
| `tools/geometry.py` | Optimization, transition states, and frequencies |
| `tools/pes.py` | One- and two-dimensional scans |
| `tools/spectrum.py` | Spectrum post-processing |
| `tools/analysis.py` | Population, orbital, and wavefunction analysis |

Inspect a script before use. Its presence does not establish method support,
scientific validity, or compatibility with the installed PySCF version.

## References

- [Workflow, method selection, convergence, and reporting safety](references/theory/workflow-and-safety.md)
- [Detailed PySCF API reference](references/theory/pyscf-api-reference.md)
- [Advanced methods, integrals, periodic systems, and performance](references/theory/pyscf-advanced.md)
- [PySCF and JAX integration](references/theory/pyscf-jax-integration.md)
- [Two-dimensional potential-energy surfaces](references/practice/2d-potential-energy-surface.md)
- [Emission-spectrum guide](references/practice/emission-spectrum-guide.md)
- [Emission-spectrum workflow](references/practice/emission-spectrum-workflow.md)
- [Benzene DFT and TDDFT example](references/benzene-dft-tddft.py)
- [Version notes](VERSION_UPDATE.md)
