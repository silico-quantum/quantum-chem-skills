---
name: multiwfn
description: Use when inspecting or analyzing molecular wavefunctions with Multiwfn, especially for orbital, population, bond-order, real-space, spectrum, excited-state, or grid analyses from fch, wfn, wfx, mwfn, Molden, cube, or coordinate files.
license: MIT
compatibility: Requires POSIX process-group support, Python 3.10+, and a local Multiwfn installation; menu numbers, prompts, supported fields, and output names may change with the version and development-build update date.
---

# Multiwfn wavefunction analysis

Choose an analysis only after proving that the input contains the required
information. Treat menu sequences as versioned executable interfaces, not as
stable prose shortcuts.

## Prerequisites

- Locate the executable with `command -v Multiwfn` or an explicit installation
  path obtained from the user or site configuration.
- Start the program once and record the complete version and update date from
  its banner. A `3.8(dev)` build from one date may not share prompts with
  another `3.8(dev)` build.
- Keep the official manual shipped with that build available. Use its section
  for the requested analysis and input format.
- Create a fresh working directory; preserve the source file unchanged.
- Confirm sufficient memory, thread, and disk limits for grids and topology
  calculations.

Do not claim Multiwfn was run when the executable or required input was absent.

## Input contract

Require a scientific target, requested quantities and units, source file,
expected molecule, charge, multiplicity, electron count, electronic state,
method/basis or density provenance, and acceptance criteria. Record complete
wavefunction provenance: generating program/version, calculation type,
geometry/state, conversion chain, and source-file hash.

### Capability matrix

| Input class | Information normally available | Safe scope | Not justified from this file alone |
|---|---|---|---|
| `.mwfn`, `.wfn`, `.wfx` | Wavefunction representation, atoms, and basis-dependent data | Supported wavefunction and real-space analyses after load validation | Properties absent from the source or generation method |
| `.fch`, `.fchk` | Gaussian formatted-checkpoint fields present in that file | Orbital, density, population, grid, and related analyses when required arrays load | Assuming every checkpoint contains excited-state, frequency, or response data |
| Molden input | Geometry, basis, and orbitals represented by the exporter | Analyses supported after checking basis conventions, occupations, ECP treatment, and electron count | Blind trust in exporter-specific normalization or missing metadata |
| `.cub`, `.cube` | Geometry plus one sampled scalar grid | Grid processing and visualization of that scalar field | Full wavefunction analysis, orbital reconstruction, population analysis, or a different field |
| `.xyz` | Element labels and coordinates only | Geometry-only functions | Electron density, orbitals, charges, bond order, topology, or spectrum |

XYZ and cube inputs cannot support complete wavefunction analysis. Request a
wavefunction-bearing source instead of fabricating missing data. Even a capable
format must be rejected when its actual fields are incomplete.

For a Gaussian source, preserve this chain:

1. include `%Chk=<name>.chk` in the Gaussian input;
2. accept the Gaussian job only after its complete log has normal termination
   and passes its scientific convergence gates;
3. run `formchk <name>.chk <name>.fchk` with the matching Gaussian utility;
4. require a successful conversion and non-empty `.fchk` before Multiwfn.

`unfchk` performs the reverse conversion and is not part of this analysis chain.

## Workflow

1. **Freeze provenance and scope.** Hash the source file and write down the
   exact requested analysis, units, expected molecule, charge, multiplicity,
   electron count, state, and source calculation.
2. **Pilot interactively.** Run `Multiwfn <input-file>` in a fresh directory.
   Capture the banner, load summary, displayed menus, prompts, selected values,
   and terminal transcript.
3. **Validate the load before choosing a menu.** Compare formula, atom count,
   total/alpha/beta electrons, net charge, expected multiplicity, restricted or
   unrestricted character, basis functions, orbital count, and occupations to
   provenance. Stop on unexplained differences.
4. **Use the installed manual.** Select the task from the menu actually printed
   by this build. Confirm grid definition, orbital/spin selection, numerical
   settings, units, output name, and overwrite behavior at every prompt.
5. **Write to a new output file.** Avoid default names shared with previous
   runs. Preserve the interactive transcript with the result.
6. **Validate the scientific result.** Apply format-specific and task-specific
   checks from the next section.
7. **Automate only after the pilot passes.** Convert the exact accepted prompt
   sequence to a response file for the same build and input class. Menu paths in
   examples or another version are not sufficient evidence.

### Fail-closed batch pattern

Run a batch only in a fresh directory, with an absolute input path and a
response file copied from an accepted interactive pilot. The response file
must end with the exact response that cleanly exits the program in that pilot.
Record the executable's resolved-path SHA-256 plus exact full-line banner,
successful-load, task-completion, and clean-exit markers. Export the resolved
path, SHA-256, and banner line from the accepted pilot; do not recompute an
expected hash from an unpiloted replacement binary.

```bash
set -e
test ! -e batch_run_001 && test ! -L batch_run_001 || exit 1
mkdir batch_run_001
: "${MULTIWFN_BIN:?Set the resolved executable path saved by the pilot}"
: "${MULTIWFN_SHA256:?Set the 64-hex executable SHA-256 saved by the pilot}"
: "${MULTIWFN_VERSION_MARKER:?Set the exact full banner line saved by the pilot}"
python3 <skill-dir>/scripts/run_batch.py /absolute/path/to/source.fchk \
  --workdir batch_run_001 \
  --responses /absolute/path/to/accepted-commands.in \
  --executable "$MULTIWFN_BIN" \
  --expected-executable-sha256 "$MULTIWFN_SHA256" \
  --timeout-seconds 1800 \
  --version-marker "$MULTIWFN_VERSION_MARKER" \
  --load-marker "<exact successful-load line from pilot>" \
  --task-marker "<exact analysis-complete line from pilot>" \
  --exit-marker "<exact clean-exit line from pilot>" \
  --final-response "<exact final quit response from pilot>" \
  --log multiwfn-transcript.log \
  --report batch-validation.json \
  --output task-specific-output.cub \
  --output task-specific-summary.txt
```

Repeat `--output` for every primary and auxiliary artifact. The helper refuses
stale paths, applies a fixed wall-clock deadline, starts a separate process
session, rejects and terminates a process group that retains any descendant
after the leader exits, and confirms the group is quiescent before reading any
artifact. It nofollow-reads the resolved executable before and after execution,
requires the pilot SHA-256 and file identity to remain bound, checks all four
exact full-line markers once and in order, rejects common invalid-input/EOF/
runtime diagnostics, hashes every non-empty output, and writes a
machine-readable report. A successful helper report uses `status:
process_validated` and
`scientific_status: pending_independent_validation`; it does not certify the
scientific contents. The accepted pilot must establish the clean-exit response
and marker; the helper cannot guess either one.

After the helper accepts the process-level run, parse every output with a
task-appropriate independent reader and apply the scientific checks below.

## Validation and acceptance

- The banner records the executable version and update date used.
- The loaded molecule, atom count, electron count, charge, multiplicity,
  alpha/beta populations, orbital occupations, and basis/ECP treatment match
  the wavefunction provenance.
- The selected function is present in the displayed menu and every response
  matches the prompt captured for this exact build.
- Orbital indices come from the loaded occupation list. For open-shell data,
  report spin channel explicitly rather than using a single ambiguous HOMO or
  LUMO label.
- Grid analyses report field identity, units, origin, axes, spacing, dimensions,
  and boundary coverage. Confirm nonzero dimensions and a complete data payload.
- Integrated or partitioned quantities include the scheme, numerical settings,
  units, convergence/integration diagnostics, and a physically motivated
  conservation or sum-rule check when available.
- A spectrum requires source transitions or modes appropriate to that spectrum;
  report broadening model and parameter, x/y units, range, and point count.
- Every requested output file is newly created, non-empty, parseable by an
  independent reader when practical, and linked to the transcript and source
  hash.
- Batch runs finish before a declared wall-clock deadline, contain exact
  successful-load, task-completion, and clean-exit markers from the interactive
  pilot, contain no invalid-selection/EOF/runtime diagnostic, and produce every
  declared fresh output.

For neutral closed-shell benzene, `6*6 + 6*1 = 42` electrons gives 21 doubly
occupied spatial orbitals; therefore HOMO is orbital 21 and LUMO is orbital 22
when Multiwfn uses one-based spatial-orbital numbering. Still verify the loaded
occupation list rather than transferring those indices to another molecule.

## Failure handling

- **Coordinate or scalar-grid input for a wavefunction task:** stop and request
  `.mwfn`, `.wfn`, `.wfx`, `.fch/.fchk`, or another demonstrably sufficient
  source. Do not synthesize orbitals or density.
- **Load-summary mismatch:** reject the analysis until electron count, charge,
  spin, ECP/core treatment, basis, and orbital occupations are reconciled.
- **Menu or prompt mismatch:** stop the batch, delete no evidence, and repeat an
  interactive pilot using the installed manual and build.
- **Missing task data:** mark the quantity `not_computed`; do not substitute a
  value from a different file or method.
- **Nonzero exit, absent marker, fatal diagnostic, or missing output file:** mark
  the batch failed even if another artifact exists.
- **Unexpected overwrite prompt or stale default output:** abort and rerun in a
  fresh directory with a unique name.
- **Numerical or physical validation failure:** retain raw outputs, report the
  failed gate, and change one justified setting at a time in a new run.

## Output and reporting

Return a manifest or equivalent record with:

- process status: `process_validated`, `rejected`, `incomplete`, or `not_run`;
- scientific status: `pending_independent_validation`, `accepted` only after
  task-specific checks, or `not_evaluated`;
- Multiwfn executable, version/update date, working directory, and thread count;
- source path/hash and complete wavefunction provenance;
- molecule, charge, multiplicity, electron count, state, and loaded summary;
- analysis name, manual section/build consulted, interactive menu sequence or
  response-file hash, and numerical settings;
- transcript path/hash and every output file path/hash;
- units, definitions, validation checks, warnings, and unresolved limitations;
- results separated from scientific interpretation.

Do not report undocumented menu paths or method-dependent reference numbers as
universal expected values.

## References

- [Multiwfn official website and current manual](https://sobereva.com/multiwfn/)
- [Official Multiwfn quick start](https://sobereva.com/multiwfn/misc/Multiwfn%20quick%20start.pdf)
- [Multiwfn English forum](https://sobereva.com/wfnbbs/)
- [Lu and Chen, Journal of Computational Chemistry 33, 580-592 (2012)](https://doi.org/10.1002/jcc.22885)
- [Lu, Journal of Chemical Physics 161, 082503 (2024)](https://doi.org/10.1063/5.0216272)
- [Local benzene provenance example](references/benzene-analysis.md)
