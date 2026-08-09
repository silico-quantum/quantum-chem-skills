---
name: gaussian
description: Use when preparing, reviewing, or troubleshooting licensed Gaussian electronic-structure calculations, especially jobs involving optimization, frequencies, excited states, checkpoint reuse, or Gaussian input and output files.
license: MIT
compatibility: Requires an authorized Gaussian installation; supported keywords and utilities depend on the installed product revision and site configuration.
---

# Gaussian electronic-structure calculations

Treat every calculation as a provenance-bearing scientific run. Never infer
success from an output file or checkpoint file merely existing, and never add
keywords whose behavior has not been checked against the installed revision.

## Prerequisites

- Confirm access to a separately licensed, authorized Gaussian installation.
- Record the product, revision, executable, scheduler, and site module used.
- Locate the matching `formchk` and `unfchk` utilities when checkpoint
  conversion is required. Do not mix utilities from an unrelated installation.
- Obtain a writable scratch directory with sufficient quota and a scheduler
  allocation consistent with `%Mem` and `%NProcShared`.
- Keep hosts, user names, ports, credentials, module names, and scratch paths in
  private site configuration, never in this reusable skill.

Do not claim a Gaussian calculation was executed when the licensed executable
was unavailable. Input review and static validation are distinct from a run.

## Input contract

Before writing the input, require and record:

1. the target property or decision and its acceptance criteria;
2. geometry source, coordinate units, molecular identity, and conformer;
3. molecular charge and multiplicity, including an electron-count and spin
   consistency check;
4. electronic state and, for excited states, the requested manifold, number of
   roots, target root, and how state identity will be tracked;
5. method, basis set or ECP, solvent/environment model, and the source of that
   protocol;
6. job type, numerical settings, resource limits, and output units;
7. new input, log, and checkpoint names plus any explicitly approved old
   checkpoint.

Fail closed when any item needed to interpret the result is unknown. Do not
silently choose a functional, basis, charge, multiplicity, solvent, state, or
unit. Refuse an existing output target unless overwrite was explicitly approved
and its provenance is preserved.

A Gaussian input has Link 0 commands, a route section, a title, the charge and
multiplicity line, and a molecule specification unless `Geom` explicitly reads
some or all of those data from a compatible checkpoint. Use
[`WORKFLOWS.md`](WORKFLOWS.md) for minimal templates.

## Workflow

1. **Freeze the scientific protocol.** Document the target, molecular state,
   geometry, method, basis/ECP, environment, and acceptance thresholds before
   editing Gaussian input.
2. **Create a fresh run identity.** Use unique input, log, and `%Chk` names.
   When reading a prior job, use `%OldChk` for the immutable source and `%Chk`
   for a new destination checkpoint. For S0-to-emission work, use the explicit
   five-stage chain in [`WORKFLOWS.md`](WORKFLOWS.md); do not re-paste geometry
   between accepted stages.
3. **Build and review the input.** Check section order, blank lines, atom count,
   geometry units, charge and multiplicity, route syntax, and resource values.
   Resolve every placeholder in a template.
4. **Check revision-specific syntax.** Consult the installed Gaussian manual or
   the official keyword pages. [`KEYWORDS.md`](KEYWORDS.md) lists only a small,
   conservative orientation set; it is not a production keyword catalogue.
5. **Submit through the authorized local or scheduler procedure.** Capture the
   exact command, exit status, scheduler job identifier, environment, and
   standard error without embedding private connection details in the skill.
6. **Parse the complete log.** Check termination, SCF and geometry convergence,
   state tracking, warnings, and task-specific acceptance criteria. For every
   excited-state stage, record excitation energy, oscillator strength, leading
   configuration/NTO character, and a state-continuity decision rather than
   following the root number alone.
7. **Post-process only an accepted source job.** Convert a binary checkpoint to
   a formatted checkpoint with:

   ```bash
   test ! -e molecule.fchk && test ! -L molecule.fchk || exit 1
   formchk molecule.chk molecule.fchk
   ```

   The reverse operation is explicitly:

   ```bash
   test ! -e molecule-restored.chk && \
     test ! -L molecule-restored.chk || exit 1
   unfchk molecule.fchk molecule-restored.chk
   ```

   `formchk` is binary checkpoint to formatted checkpoint; `unfchk` is
   formatted checkpoint to a machine-local binary checkpoint. Never describe
   them as interchangeable alternatives.

## Validation and acceptance

Apply all relevant gates:

- The process exit status is successful and the final job step contains
  `Normal termination of Gaussian`; no later error or truncated section exists.
- Every required SCF step converged. An `SCF Done` line by itself is not enough
  to accept an optimization, frequency, or excited-state workflow.
- An optimization reports convergence to a stationary point and the accepted
  geometry is the final converged geometry, not merely the last printed one.
- A claimed minimum has no chemically meaningful imaginary frequencies. A
  claimed first-order transition state has exactly one intended negative mode,
  and its displacement is inspected; use IRC when reaction-path assignment is
  required. Report the numerical threshold used to distinguish noise.
- Excited-state jobs retain the intended state character, not just the same
  root number, across the geometry path.
- A claimed S1 minimum has a compatible accepted excited-state frequency
  analysis. Otherwise label it `optimized S1 stationary geometry` and record
  `minimum_not_verified`.
- A reported vertical emission energy is the target-state TD excitation energy
  at the accepted excited-state geometry, with oscillator strength, units, and
  state-character evidence. It is not a difference of raw SCF reference
  energies from separate jobs.
- Charge, multiplicity, atom count, method, basis/ECP, solvent, and units in the
  output agree with the input contract.
- A required checkpoint is non-empty and belongs to this accepted job. Its
  existence does not replace log validation.
- A generated formatted checkpoint is non-empty and `formchk` exits
  successfully before it is supplied to downstream software.

Record warnings and deviations even when all acceptance gates pass.

## Failure handling

- **No normal termination:** mark the run failed or incomplete; inspect the end
  of the log and scheduler diagnostics before changing anything.
- **Input or keyword error:** correct only the diagnosed syntax using the
  installed revision's documentation. Do not pile on speculative keywords.
- **SCF failure:** first check geometry, charge/spin, method suitability, linear
  dependence, and numerical diagnostics. Change one justified setting at a
  time and create a new run identity.
- **Optimization failure:** inspect forces, steps, constraints, and the actual
  structure. Do not accept the final printed geometry as optimized.
- **Unexpected frequencies:** inspect normal modes and geometry; do not delete
  modes, take absolute values, or relabel a saddle point as a minimum.
- **Excited-state root change:** stop state-specific interpretation, identify
  state character, and restart only with a documented tracking strategy.
- **Checkpoint read failure:** verify revision compatibility, source-job
  acceptance, and direction of conversion. Preserve the original checkpoint.
- **Resource or license failure:** resolve the site allocation or authorization;
  it is not a scientific convergence problem.

Never overwrite the failed log. Keep it linked to the corrected rerun.

## Output and reporting

Return a compact run record containing:

- status: `accepted`, `failed`, `incomplete`, or `not_run`;
- molecular identity, geometry source, charge, multiplicity, and target state;
- Gaussian product/revision and executable provenance;
- method, basis/ECP, environment, numerical settings, and units;
- scheduler/job identifier, resources, command, exit status, and elapsed time;
- paths and hashes for input, complete log, accepted geometry, checkpoint, and
  formatted checkpoint when created;
- the exact termination, SCF, optimization, frequency, and state-tracking gates
  applied;
- requested observables with units, warnings, deviations, and unresolved
  limitations.

Separate measured output from interpretation. Leave unavailable values
explicitly `not_computed`; never fill them with estimates or placeholders.

## References

- [Gaussian official keyword index](https://gaussian.com/keywords/)
- [Gaussian Link 0 commands](https://gaussian.com/link0/)
- [Gaussian geometry keyword](https://gaussian.com/geom/)
- [Gaussian optimization keyword](https://gaussian.com/opt/)
- [Gaussian formatted-checkpoint guidance](https://gaussian.com/wp-content/uploads/dl/remote.pdf)
- [Local conservative keyword orientation](KEYWORDS.md)
- [Local reviewed input templates](WORKFLOWS.md)
