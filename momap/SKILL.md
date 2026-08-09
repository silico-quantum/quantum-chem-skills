---
name: momap
description: Use when preparing, running, validating, or reporting licensed MOMAP photophysics workflows from Gaussian outputs, including EVC, vibronic spectra, and explicitly parameterized ISC calculations.
license: MIT
compatibility: Requires Python 3, a licensed MOMAP installation, and Gaussian log plus formatted-checkpoint files; selected workflows also require formchk, MPI, PySOC, or Slurm.
---

# MOMAP

Use MOMAP only through the syntax and executables shipped with the licensed
local installation. This skill documents the repository's MOMAP 2024A-style
interface; confirm it against the local examples before adapting another
release. Never convert a missing parameter or failed stage into a numerical
result.

## Prerequisites

- Source the licensed installation environment or set `MOMAP_ROOT` so that the
  `momap` launcher is discoverable. Preserve the administrator-provided
  `MOMAP_LICENSE`; the bundled runner does not replace it.
- Use the Python interpreter for the active environment. `formchk` is required
  only when a Gaussian `.chk` must be converted to `.fchk`.
- For Slurm, use the site's partition and MPI policy. The wrapper has no
  site-specific partition default.
- Record the MOMAP release, Gaussian release, Python version, MPI
  implementation, host, and command before production work.

Bind every run to the reviewed licensed launcher and its exact runtime banner:

```bash
MOMAP_BUILD=2024A
MOMAP_LAUNCHER=$(command -v momap)
test -f "$MOMAP_LAUNCHER" && test ! -L "$MOMAP_LAUNCHER" || exit 1
MOMAP_LAUNCHER_SHA256=$(python3 -c \
  'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' \
  "$MOMAP_LAUNCHER")
: "${MOMAP_VERSION_BANNER:?Set the exact full-line banner verified for this licensed 2024A build}"
```

Do not infer `2024A` from filenames or example output. The runner re-reads the
original launcher without following symbolic links, matches its SHA-256,
captures `momap_runner.log`, and requires the exact banner once.

Preflight without running a calculation:

```bash
command -v momap
python3 tools/extract.py --help
python3 tools/runner.py --help
python3 tools/tadf.py --help
```

Run these relative paths from the `momap/` skill directory, or replace them
with the absolute path to that directory when operating in a calculation
worktree.

## Input contract

Create a calculation manifest before execution. It must identify:

- molecule and geometry provenance for every state;
- charge, multiplicity, electronic state, method, basis, solvent model, and
  Gaussian route for S0, S1, and T1 where applicable;
- the state/root mapping used throughout optimization and property extraction;
- temperature, time step, correlation time, broadening, and all MOMAP units;
- the source and uncertainty of every transition dipole, nonadiabatic
  coupling, and spin-orbit coupling;
- expected output paths and whether each output may be created.

EVC requires matched Gaussian `.log` and `.fchk` files for both states. Confirm
successful Gaussian termination, compatible atom order, atom count, isotope
assignment, method, and frequency treatment. A file name is not evidence that
an optimization, frequency calculation, or target state succeeded.

For a Gaussian linear-response TDDFT S1 log, compute the state total energy on
a common reference:

```text
E_S1,total(R_S1) = E_reference,SCF(R_S1) + omega_S1(R_S1)
Ead = E_S1,total(R_S1) - E_S0,total(R_S0)
```

`omega_S1` is the final state-1 vertical excitation at the S1 geometry. Both
total-energy terms must use the same Hamiltonian and energy convention. Never
use `E_reference,SCF(R_S1) - E_S0,total(R_S0)` as Ead. Apply zero-point or
thermal corrections only when they are defined consistently for both states.

ISC additionally requires a finite positive `Hso` in `cm^-1` from a named
computed or measured source and a signed S1-T1 adiabatic gap in Hartree. In the
bundled pipeline, `E_T1` is the final Gaussian `SCF Done` total electronic
energy at the T1 geometry; it is not a TD excitation, zero-point-corrected
energy, or free energy. `E_S1` is the reconstructed linear-response total
energy defined above. Both must use compatible electronic-structure
conventions. The expected ordering is `E_S1 - E_T1 > 0`. Never apply `abs()` to
hide the opposite ordering: flag it, verify state identity, and fail a requested
ISC calculation. Missing SOC means `not_computed`, not a default coupling.

## Workflow

1. Copy source data into a new work directory without modifying the originals.
2. Verify Gaussian termination, state identity, atom order, and matching
   `.fchk` files. Hash the source files.
3. Extract Ead and transition dipoles:

   ```bash
   set -e
   test ! -e parameters.json && test ! -L parameters.json || exit 1
   test ! -e parameters.json.partial && \
     test ! -L parameters.json.partial || exit 1
   python3 tools/extract.py --s0 s0.log --s1 s1.log --json \
     > parameters.json.partial
   test -s parameters.json.partial
   mv parameters.json.partial parameters.json
   ```

   The extractor uses the first and last state-1 transition-dipole tables for
   absorption and emission. Review that mapping in the actual Gaussian log;
   reject root changes or ambiguous tables.

4. Prefer the fail-closed pipeline for a complete S1-to-S0 spectrum and optional
   S1-to-T1 ISC calculation:

   ```bash
   test ! -e ./momap-results && test ! -L ./momap-results || exit 1
   test ! -e ./molecule-001-result.json && \
     test ! -L ./molecule-001-result.json || exit 1
   python3 tools/tadf.py molecule-001 \
     --s0 s0.log --s1 s1.log --t1 t1.log \
     --hso-cm1 "$HSO_CM1" --output ./momap-results \
     --json-output ./molecule-001-result.json \
     --expected-build "$MOMAP_BUILD" \
     --expected-launcher-sha256 "$MOMAP_LAUNCHER_SHA256" \
     --expected-version-banner "$MOMAP_VERSION_BANNER"
   ```

   Use a new molecule ID and output directory. The pipeline refuses an existing
   molecule directory. It isolates generated files under `evc_s1_s0/`,
   `spectrum/`, `evc_s1_t1/`, and `isc/`; each accepted stage receives a
   `stage_receipt.json` containing input/output SHA-256 hashes, sizes, and
   timestamps. This prevents the two EVC calculations from overwriting each
   other's `evc.out` or `evc.cart.dat`.
5. For a manual EVC stage, use a dedicated new directory and the local interface:

   ```text
   do_evc = 1
   &evc
     ffreq(1) = "s0.log"
     ffreq(2) = "s1.log"
     sort_mode = 1
   /
   ```

   ```bash
   for target in nodefile evc.out evc.cart.dat evc.dint.dat; do
     test ! -e "$target" && test ! -L "$target" || exit 1
   done
   python3 tools/runner.py momap_evc.inp --workdir "$PWD" \
     --expected-build "$MOMAP_BUILD" \
     --expected-launcher-sha256 "$MOMAP_LAUNCHER_SHA256" \
     --expected-version-banner "$MOMAP_VERSION_BANNER"
   ```

6. Generate and review the fluorescence input from the extracted values. Run
   `spec_tvcf` only after `evc.cart.dat` is accepted. Do not reuse another
   molecule's Ead, EDMA, or EDME.
7. For standalone ISC input generation, first build and accept the S1-T1 EVC
   file in its own directory, then use the repository's local
   `isc_tvcf` interface. The generator requires explicit SOC:

   ```bash
   test ! -e momap_isc.inp && test ! -L momap_isc.inp || exit 1
   python3 tools/oled.py isc --evc-dat evc.cart.dat \
     --ead "$ISC_EAD_AU" --hso "$HSO_CM1" -o momap_isc.inp
   for target in nodefile isc.tvcf.log isc.tvcf.spec.dat; do
     test ! -e "$target" && test ! -L "$target" || exit 1
   done
   python3 tools/runner.py momap_isc.inp --workdir "$PWD" \
     --expected-build "$MOMAP_BUILD" \
     --expected-launcher-sha256 "$MOMAP_LAUNCHER_SHA256" \
     --expected-version-banner "$MOMAP_VERSION_BANNER"
   ```

   The generator rejects a missing or empty EVC file, non-positive/non-finite
   Ead, temperature, correlation time, or SOC, and any existing output path.
   The generated input uses `do_isc_tvcf_ft`, `do_isc_tvcf_spec`, and the
   `&isc_tvcf` block used by the bundled local examples. Confirm those names
   against the licensed installation before execution.
8. Plot only a validated, non-empty spectrum to a new path:

   ```bash
   test ! -e spectrum.png && test ! -L spectrum.png || exit 1
   python3 tools/plot.py spec.tvcf.spec.dat \
     --energy-window 1 5 --output spectrum.png
   ```

   The plotter refuses overwrite. Its title is user-supplied; it does not invent
   a MOMAP version or temperature. Put the actual version and temperature in
   the manifest or `--title` only after verifying them.
9. Archive input, stdout/stderr, exit status, outputs, source hashes, stage
   receipts, and the
   final manifest before interpreting rates or spectral peaks.

## Validation and acceptance

Accept a stage only when all applicable checks pass. The automated validators
currently implement this MOMAP 2024A-style local fixture contract:

- each Gaussian state uses the supported two-job chain: one completed
  optimization followed by one harmonic-frequency job, both ending in unique
  exact Gaussian normal-termination lines, with no fatal tail;
- S0/S1/T1 multiplicities are 1/1/3, respectively; charge and atom order remain
  unchanged across the two jobs and match `Number of atoms`, `Atomic numbers`,
  `Charge`, and `Multiplicity` in the paired formatted checkpoint;
- the S1 optimization remains on target root 1, and both jobs contain state-1
  excitation and transition-dipole evidence;

- EVC: `evc.out` contains `Normal finish of evc calculation`;
- spectrum: `spec.tvcf.log` contains `Normal finish of spec_tvcf calculation`;
- ISC: `isc.tvcf.log` contains `Normal finish of isc_tvcf calculation`.

These strings are case-insensitive but otherwise exact. They are repository
compatibility contracts, not universal MOMAP promises. A different release
must fail closed until its local output is reviewed and a tested marker is
added. In addition:

- the process exit code is zero; `momap_runner.log` contains the caller-supplied
  exact full-line 2024A version banner exactly once; and the module log shows
  its documented normal completion without fatal diagnostics;
- every required output exists, is non-empty, and was created after that stage
  began;
- EVC reports the expected atom and mode counts and no unexplained mode loss;
- extracted state energies reproduce the declared Ead formula and units;
- every non-empty, non-comment spectrum row has exactly seven numeric columns
  for the local 2024A fixture; at least three rows are present, all values are
  finite, the Hartree/eV and eV/wavenumber conversions agree within 1%, all
  three intensity columns are non-negative, the eV axis is strictly monotonic,
  the wavelength axis runs in the opposite direction, and each eV/nm pair
  agrees with `hc` within 1%. Reject tail text, partial rows, and fatal
  diagnostics even if a completion marker also appears;
- an ISC result names the SOC source and value and reports finite ISC and RISC
  rates with explicit `s^-1` units. When labels are present, `ISC` and `RISC`
  must each occur exactly once; duplicates, conflicts, or mixed labelled and
  unlabelled records fail. The versioned local 2024A fallback accepts exactly
  two unlabelled rate records and records `first=ISC, second=RISC`;
- rerunning the parser on archived outputs reproduces the reported values.

Do not infer a rate unit from magnitude. Do not call a spectrum, rate, quantum
yield, or mobility successful merely because a file was created.

## Failure handling

- Missing `.fchk`: choose a new target, reject both an existing path and a
  dangling symbolic link, run `formchk` on the matching `.chk`, then verify the
  target is non-empty and belongs to the same Gaussian job.
- Missing final state-1 excitation: stop; an S1 total energy cannot be recovered
  from the SCF reference alone.
- Root change or inconsistent atom order: stop and repair the electronic-state
  or geometry provenance before rerunning EVC.
- Missing SOC or T1 data: report ISC as `not_computed`. Never insert a test value.
- Non-positive signed S1-T1 gap: record
  `unexpected_S1_not_above_T1`; fail requested ISC and review state/energy
  provenance. Never replace the signed gap with its absolute value.
- Missing, empty, stale, unitless, non-finite, or unparsable ISC output: report
  ISC as `failed`, keep global `success: false`, and retain the stage directory.
- Empty or all-zero spectrum: retain the raw file and mark the spectrum failed;
  do not normalize or select a peak.
- MOMAP, EVC, spectrum, or requested ISC nonzero exit: keep `success: false` and
  record the failing stage, command, exit code, and log path.
- MPI launcher incompatibility: use `tools/runner.py`, which reads the original
  launcher without following links, verifies its expected SHA-256, and creates
  a private mode-0700 temporary patched launcher. Receipts record original and
  patched hashes plus the exact replacement count; do not maintain a permanent
  hand-edited copy.
- License or scheduler error: stop and ask the site administrator; do not alter
  license paths or choose a guessed partition.

## Output and reporting

Write a machine-readable result and a short human summary. At minimum include:

- molecule ID, source hashes, geometry/state provenance, and atom count;
- software versions, command, working directory, and exit codes;
- method, basis, charge, multiplicity, state/root mapping, and units;
- `E_S0`, reference SCF energy at S1, `omega_S1`, reconstructed `E_S1,total`,
  Ead, EDMA, and EDME with their sources;
- EVC output paths, validation status, receipt paths, and SHA-256 hashes;
- spectrum path, energy window, peak rule, and peak only if validated;
- every rate with the exact output label, rate units, and parser contract;
- ISC status as `computed`, `failed`, or `not_computed`, plus the SOC provenance;
- warnings, failed stages, approximations, and unresolved questions.

Keep computed evidence separate from interpretation. A blue-window flag, TADF
assessment, or device conclusion is downstream interpretation, not a MOMAP
termination criterion.

## References

- [Repository worked examples](EXAMPLES.md)
- [Repository MOMAP 2024A-style quick reference](QUICKREF.md)
- [MOMAP product page](https://www.hzw.ai/momap/index)
- Y. Niu et al., *Molecular Physics* 116 (2018), DOI:
  `10.1080/00268976.2017.1402966`

The repository examples are compatibility aids, not evidence that a numerical
result has been reproduced on the current machine.
