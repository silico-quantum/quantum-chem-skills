# PySCF Workflow and Scientific-Safety Reference

## Purpose

This reference expands the compact workflow in `pyscf/SKILL.md`. Use it when
selecting a method, diagnosing convergence, designing a scan, or deciding whether
a result is ready to report.

## Calculation Definition

Write a calculation specification before writing code:

| Field | Minimum information |
|---|---|
| System | Structure identifier, geometry source, atom order, coordinate unit |
| Electronic state | Charge, `2S`, multiplicity label, restricted/open-shell reference |
| Target | Total energy, relative energy, gap, excitation, gradient, geometry, or property |
| Hamiltonian | Nonrelativistic, scalar relativistic, spin-orbit, ECP, external field |
| Method | HF/DFT/post-HF/multireference/excited-state method and approximations |
| Basis | Orbital basis, auxiliary basis, ECP, frozen-core convention |
| Environment | Vacuum, PCM/SMD parameters, explicit environment, periodic cell |
| Numerical controls | SCF threshold, grid, integral threshold, optimizer criteria |
| Software | PySCF and plugin versions, backend, dtype, threads, hardware |

If two results differ in any field that can affect the comparison, describe the
difference rather than presenting a single unlabeled energy gap.

## Method-Selection Ladder

### Mean-field reference

1. Determine whether a closed-shell, open-shell, broken-symmetry, or periodic
   reference is physically intended.
2. For DFT, choose a functional using property- and regime-specific evidence.
   Charge transfer, diffuse states, barriers, transition metals, dispersion, and
   near-degeneracy require different validation.
3. Converge the orbital basis, ECP/relativistic treatment, and DFT integration grid
   to the precision required by the scientific conclusion.
4. Inspect occupations, orbital character, symmetry, `<S^2>` where applicable,
   and SCF stability.

### Single-reference correlation

Use MP2 or coupled-cluster methods only after checking that the reference is
qualitatively adequate. Record frozen-core and density-fitting approximations.
Assess basis-set convergence and, when relevant, diagnostics for strong
correlation. Perturbative triples do not repair an invalid single-reference model.

### Multireference calculations

Define the active space by orbitals and electrons, not only by two integers. Save
the selected orbitals and report state averaging, weights, roots, natural
occupations, and any orbital localization or reordering. Test whether the result
is stable to reasonable active-space changes.

### Excited states

State the response theory, reference, spin channel, number of roots, and root
indexing convention. Inspect transition vectors, dominant configurations, NTOs,
oscillator strengths, and state symmetry. Linear-response TDDFT/TDA is not a
reliable description of every double-excitation or strong multireference state.

### Periodic systems

Use `pyscf.pbc` APIs and document lattice vectors, dimensionality, vacuum,
k-point mesh, pseudopotentials, Coulomb treatment, smearing, and finite-size
convergence. A molecular input with repeated coordinates is not a periodic model.

## Basis, Grid, and Model Convergence

Run a convergence sequence that targets the reported observable:

1. Start with a small basis and moderate grid to validate input and workflow.
2. Increase basis quality and add diffuse or core-valence functions when the
   property requires them.
3. Increase the DFT grid until the energy difference, gradient, or response
   property is stable at the intended precision.
4. Compare conventional and density-fitted results for at least representative
   points if the approximation could affect the conclusion.
5. For relative quantities, use a consistent method and basis across all systems.

Do not call a basis “high accuracy” without a convergence result for the selected
property. Heavy elements may require an ECP or relativistic Hamiltonian and a
compatible basis; choose them as one protocol.

## Convergence Decision Procedure

### SCF

1. Validate geometry, units, charge, spin, electron parity, basis, and ECP.
2. Inspect the iteration log, orbital occupations, energy changes, and density
   residual rather than only the final printed energy.
3. Try a physically motivated initial guess or continuation from a nearby geometry.
4. If oscillatory, test damping or level shifting; if close to a solution, test
   Newton SCF.
5. Run stability analysis and compare plausible restricted, unrestricted, or
   broken-symmetry solutions.
6. Preserve all competing solutions and select among them using energy and state
   character, not convergence speed.

Do not switch to an unrelated Hamiltonian, lower the convergence threshold, or
discard a lower-energy solution merely to obtain a successful status.

### Post-SCF and response calculations

Downstream methods inherit defects in the reference. Confirm reference convergence
before each post-SCF calculation, then inspect the method-specific convergence
status, residual, root count, amplitudes, and diagnostics. A finite number in an
output file is not proof of convergence.

### Geometry optimization

At every step, distinguish electronic failure from optimizer failure. Record the
maximum gradient, displacement, energy change, constraints, and termination reason.
After convergence, recompute the electronic state and use a compatible vibrational
analysis to classify the stationary point.

### Excited-state optimization

Track the physical state across geometries using transition character, NTO overlap,
symmetry, or another declared criterion. Near crossings, a fixed root index can
follow a different state. Report crossings and failed tracking explicitly.

## Energy and Unit Ledger

Maintain a machine-readable or tabular ledger with separate columns for:

- ground-state total energy in Hartree;
- correlation correction in Hartree;
- excitation energy in Hartree and its single converted value in eV;
- excited-state total-energy estimate in Hartree;
- relative energy with declared zero and conversion;
- geometry unit and gradient unit;
- state index, state character, and convergence status.

For a vertical excitation at geometry `R`:

```text
E_excitation(R) = omega_i(R)
E_excited,total(R) = E_ground,total(R) + omega_i(R)
```

Do not plot `omega_i(R)` as an excited-state total PES. Do not apply the same unit
conversion in both the calculation layer and plotting layer.

## PES and Spectroscopy Design

For every scan, define coordinate ranges, grid spacing, constraints, rigid versus
relaxed treatment, state-tracking rule, and failed-point policy before execution.
Never replace failed points with zero or interpolate them silently.

For absorption and emission:

- absorption is normally evaluated from the initial-state geometry;
- vertical fluorescence emission requires an excited-state geometry and the
  corresponding downward energy difference;
- adiabatic gaps require separately optimized states and a clearly stated
  zero-point-energy convention;
- broadened spectra require a line shape, width, grid, and normalization rule;
- solvent and vibronic effects are separate modeling choices.

## Reproducibility Record

Archive:

- exact input structure and script;
- stdout/stderr and termination status;
- checkpoints and hashes for important artifacts;
- PySCF, Python, BLAS, optimizer, and plugin versions;
- CPU/GPU model, thread count, memory limit, and elapsed time for benchmarks;
- random seeds or stochastic settings, if any;
- every approximation and convergence threshold;
- failed calculations that influence point selection or interpretation.

Before publishing a value, rerun a representative calculation from the archived
input and confirm that the parsed result, unit conversion, and provenance match the
reported table or figure.

## Reporting Gate

A result is ready to report only when all applicable items are true:

- [ ] The requested physical quantity is unambiguous.
- [ ] Geometry, units, charge, spin, state, method, and basis are recorded.
- [ ] The SCF reference and downstream method converged.
- [ ] State identity and stationary-point character were checked where relevant.
- [ ] Basis, grid, finite-size, or model sensitivity supports the claimed precision.
- [ ] Unit conversions occur exactly once and retain sufficient precision.
- [ ] Logs, checkpoints, versions, and failures are archived.
- [ ] The text distinguishes computed evidence from interpretation.

If a box cannot be checked, report the missing validation and keep the conclusion
provisional.
