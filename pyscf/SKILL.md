---
name: pyscf
description: Use when defining, running, validating, or reporting molecular PySCF HF, DFT, or linear-response calculations with explicit state, units, convergence, and software provenance.
license: MIT
compatibility: Requires Python 3 and a separately installed PySCF environment; optional methods, optimizers, accelerators, and plugins require their own version-specific validation.
---

# PySCF

Use this skill to build a small, traceable calculation from documented PySCF
core APIs. Choose a method only after defining the physical quantity. A finite
energy is not an acceptance criterion, and a bundled historical prototype is
not evidence that an API is supported.

## Prerequisites

Create an isolated environment, install the intended PySCF build, and record
the exact interpreter and PySCF version:

```bash
python3 -m pip install pyscf
python3 -c 'import platform, pyscf; print(platform.python_version(), pyscf.__version__)'
```

Before production work, also record BLAS/backend information, thread settings,
hardware, and optional plugin versions. Do not assume that an example written
for one PySCF version is compatible with another.

The only supported bundled executable is
[`scripts/run_safe_dft_tda.py`](scripts/run_safe_dft_tda.py), an intentionally
narrow closed-shell RKS/TDA runner. The files under `tools/`, the old
`scripts/dft_calculation.py`, and executable examples under `references/` are
historical prototypes. Their direct CLIs fail closed; see
[legacy status](tools/README.md).

## Input contract

Define these fields before constructing `gto.M`:

- target quantity and comparison being made;
- structure identifier, geometry source, atom order, and coordinate units;
- total charge and spin, where PySCF `spin = Nalpha - Nbeta = 2S`, not the
  multiplicity;
- electronic state and restricted, restricted-open-shell, or unrestricted
  reference choice;
- Hamiltonian, method or XC functional, orbital basis, auxiliary basis, ECP,
  relativistic treatment, solvent/environment, and symmetry choice;
- SCF threshold, maximum cycles, DFT grid, memory, thread count, checkpoint,
  and all downstream method controls;
- output directory, overwrite policy, and intended unit conversions.

Reject an input when the electron count is non-positive, `Nalpha` or `Nbeta`
would be negative/non-integral, the charge and spin parity is inconsistent, a
basis/ECP assignment is missing, or the geometry units are unknown.

The supported runner accepts exactly these JSON fields: `atom`, `unit`,
`basis`, `charge`, `spin`, `xc`, `grid_level`, `conv_tol`, `max_cycle`, and
`nstates`. It rejects unknown or missing fields, non-finite controls, ambiguous
units, and nonzero `spin`. That last restriction is a runner boundary, not a
claim that PySCF lacks open-shell methods.

## Workflow

### 1. Prepare one explicit configuration

For a small interface test, save a JSON object such as:

```json
{
  "atom": "H 0.0 0.0 0.0; H 0.0 0.0 0.74",
  "unit": "Angstrom",
  "basis": "sto-3g",
  "charge": 0,
  "spin": 0,
  "xc": "pbe0",
  "grid_level": 3,
  "conv_tol": 1e-9,
  "max_cycle": 80,
  "nstates": 2
}
```

This is an interface example, not a recommended scientific protocol. Select and
converge the functional, basis, grid, state count, and thresholds for the
actual quantity.

### 2. Run only into a fresh directory

```bash
python3 scripts/run_safe_dft_tda.py \
  --config h2-pbe0.json \
  --output h2-pbe0-run
```

The runner creates the output directory with `exist_ok=False`. Before importing
PySCF it writes the exact input copy and a `run_manifest.partial.json` whose
state is `running`. A failure changes that partial manifest to `failed` and
never publishes accepted results. Acceptance atomically publishes
`results.json` and `run_manifest.json`, removes the partial manifest, and stores
the elapsed time, Python/PySCF version, and SHA-256 hashes for the result,
checkpoint, and log.

### 3. Apply every acceptance gate

The supported runner sets the checkpoint before SCF, rejects unconverged SCF,
and requires finite SCF energy, orbital energies, orbital coefficients, and
occupations. The energy returned by `mf.kernel()` must match `mf.e_tot` within
`1e-10` Hartree. Every RKS occupation must lie in `[0, 2]`, and their sum must
match `mol.nelectron` within `1e-8` electrons. The checkpoint must exist and be
non-empty before TDA begins. The runner then requires PySCF's internal and
external restricted-SCF stability statuses to both be true; an unavailable,
malformed, or unstable result remains a failed partial run and never reaches
TDA or accepted publication.

For TDA, all three returned vectors must be one-dimensional, finite, exact, and
inside the accepted physical domain:

```python
len(td.converged) == td.nstates
len(td.e) == td.nstates
len(oscillator_strength) == td.nstates
all(td.converged)
all(value > 0.0 for value in td.e)
all(value >= 0.0 for value in oscillator_strength)
```

Returning fewer roots than requested is a failed run, even if every returned
root converged. Zero or negative excitation energies and negative oscillator
strengths also fail the run and cannot reach accepted result publication. Use
`td.oscillator_strength()`; do not substitute an undocumented `td.f` attribute.
`td.e[i]` is an excitation energy in Hartree, not a total energy or an emission
result.

### 4. Extend only behind a new tested boundary

For open-shell references, MP2, CC, multireference, gradients, geometry
optimization, solvent, periodic, relativistic, GPU, or other response
workflows, consult the official module documentation and installed examples.
Add a separate version-pinned runner and RED-to-GREEN contract tests before
advertising support. Do not re-enable a historical CLI as a shortcut.

## Validation and acceptance

Accept a molecular result only when all applicable checks pass:

- atom order, geometry source, units, charge and spin reproduce the manifest;
- electron parity and basis/ECP coverage are valid;
- the intended restricted/open-shell reference was actually constructed;
- `mf.converged` is true and the final energy and orbitals are finite;
- occupations, symmetry, `<S^2>`, and stability are appropriate for the state;
- every downstream method/root reports its own convergence status;
- every accepted TDA root has a strictly positive excitation energy in Hartree
  and a non-negative oscillator strength;
- basis, DFT grid, reference, and numerical sensitivity support the reported
  precision;
- checkpoint and output hashes resolve to the accepted run;
- the report contains the Python and PySCF version and all units.

For spectra, record whether transitions are absorption or emission, state/root,
oscillator-strength convention, line shape, width, energy grid, and
normalization. A ground-geometry vertical excitation is not an emission
calculation.

## Failure handling

- SCF failure: recheck geometry, units, charge and spin first; inspect guesses,
  occupations, and competing solutions before trying damping, level shifting,
  or Newton SCF.
- Different solutions from different guesses: run stability analysis and
  characterize each solution; do not select only by convergence speed.
- Unexpected spin contamination: retain the failed result and compare justified
  restricted-open-shell/unrestricted references.
- Response root failure or reordering: increase diagnostics and track physical
  character; never substitute another root silently.
- Unsupported attribute/import/API: classify the script as incompatible with
  the installed PySCF version and consult official documentation. Do not guess a
  replacement API.
- Optimization or frequency failure: separate electronic convergence from the
  optimizer/Hessian status. Do not label a stationary point without compatible
  vibrational analysis.
- Memory failure: estimate tensor sizes and consider a smaller basis, direct or
  density-fitted algorithm, or frozen-core approximation only as a declared
  protocol change.

Return `failed` or `not_computed` for missing validation. Preserve logs and
partial artifacts, but exclude them from accepted result tables.

## Output and reporting

Report, at minimum:

- structure ID/hash, geometry source, atom count, and coordinate units;
- charge and spin (`Nalpha - Nbeta`), multiplicity label, state, and reference;
- method/XC, basis/ECP, grid, solvent/environment, thresholds, and
  approximations;
- Python and PySCF version, plugins, backend, threads, memory, hardware, command,
  and elapsed time;
- convergence flags, stability/spin diagnostics, warnings, and failed attempts;
- total energies in Hartree; each converted value with one named conversion;
- excitation energies and oscillator strengths separately from total energies;
- checkpoint/log/result paths and hashes;
- sensitivity checks, limitations, and any `not_computed` properties.

Keep theoretical interpretation, implementation status, and measured output
separate. Never describe an unexecuted example as a computed result.

## References

- [Official PySCF quickstart](https://pyscf.org/quickstart.html)
- [Official SCF user guide](https://pyscf.org/user/scf.html)
- [Official DFT user guide](https://pyscf.org/user/dft.html)
- [Official TD-SCF API](https://pyscf.org/pyscf_api_docs/pyscf.tdscf.html)
- [Official density-fitting guide](https://pyscf.org/user/df.html)
- [Repository scientific-safety checklist](references/theory/workflow-and-safety.md)
- [Status of bundled legacy tools and references](tools/README.md)

Quarantined historical material is linked here for repository completeness,
not as an executable or authoritative API guide:

- [Historical version note](VERSION_UPDATE.md)
- [Historical 2D PES draft](references/practice/2d-potential-energy-surface.md)
- [Historical emission guide](references/practice/emission-spectrum-guide.md)
- [Historical emission workflow](references/practice/emission-spectrum-workflow.md)
- [Historical advanced-API draft](references/theory/pyscf-advanced.md)
- [Historical API-reference draft](references/theory/pyscf-api-reference.md)
- [Historical JAX-integration draft](references/theory/pyscf-jax-integration.md)

Do not execute those historical documents until each call, unit convention, and
failure path has been verified against the installed PySCF version.
