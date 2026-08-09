# Benzene provenance and orbital-index example

This example demonstrates a Gaussian checkpoint chain and Multiwfn load
validation. The scientific method and basis are illustrative, not a recommended
production protocol. The example has not been executed in this repository.

## 1. Create a Gaussian checkpoint

Save the following as `benzene.com` after confirming the method, basis, geometry,
resources, and Gaussian revision for the intended use:

```bash
for target in benzene.com benzene.log benzene.chk; do
  test ! -e "$target" && test ! -L "$target" || exit 1
done
```

Only after these guards succeed, create the input at `benzene.com` and submit
the job with the site's authorized procedure, directing its log to
`benzene.log`.

```text
%Chk=benzene.chk
%NProcShared=8
%Mem=8GB
#P B3LYP/cc-pVDZ SP

Benzene single point for wavefunction analysis

0 1
C   0.000000   1.395000   0.000000
C   1.208543   0.697500   0.000000
C   1.208543  -0.697500   0.000000
C   0.000000  -1.395000   0.000000
C  -1.208543  -0.697500   0.000000
C  -1.208543   0.697500   0.000000
H   0.000000   2.479000   0.000000
H   2.150000   1.239500   0.000000
H   2.150000  -1.239500   0.000000
H   0.000000  -2.479000   0.000000
H  -2.150000  -1.239500   0.000000
H  -2.150000   1.239500   0.000000

```

Run it only on an authorized Gaussian installation. Accept the source job only
after checking the complete log for normal termination, SCF convergence,
molecular identity, charge, multiplicity, method, basis, and warnings. Then use
the matching Gaussian utility:

```bash
test ! -e benzene.fchk && test ! -L benzene.fchk || exit 1
formchk benzene.chk benzene.fchk
test -s benzene.fchk
shasum -a 256 benzene.com benzene.log benzene.chk benzene.fchk
```

`formchk` converts binary `.chk` to formatted `.fchk`. It cannot repair a failed
source calculation, and `unfchk` is the reverse conversion rather than an
alternative spelling.

## 2. Derive the expected load invariants

Neutral benzene has:

- nuclear charge: `6 carbon * 6 + 6 hydrogen * 1 = 42`;
- molecular charge: 0;
- total electron count: 42;
- singlet restricted closed shell: 21 doubly occupied spatial orbitals;
- one-based HOMO index: 21;
- one-based LUMO index: 22.

These indices follow from this input and occupation model. They are not generic
benzene constants to reuse after changing charge, multiplicity, wavefunction
type, or orbital representation.

## 3. Pilot Multiwfn interactively

Create a fresh directory, copy the immutable `.fchk` into it without replacing
an existing path or dangling symbolic link, and start:

```bash
set -e
SOURCE_FCHK="$PWD/benzene.fchk"
RUN_DIR="$PWD/multiwfn_benzene_run_001"
MULTIWFN_COMMAND="$(command -v Multiwfn)"
MULTIWFN_BIN="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$MULTIWFN_COMMAND")"
test -s "$SOURCE_FCHK"
test -x "$MULTIWFN_BIN"
MULTIWFN_SHA256="$(python3 -c 'import hashlib, pathlib, sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$MULTIWFN_BIN")"
test ! -e "$RUN_DIR" && test ! -L "$RUN_DIR" || exit 1
mkdir "$RUN_DIR"
test ! -e "$RUN_DIR/benzene.fchk" && \
  test ! -L "$RUN_DIR/benzene.fchk" || exit 1
cp "$SOURCE_FCHK" "$RUN_DIR/benzene.fchk"
cd "$RUN_DIR"
"$MULTIWFN_BIN" benzene.fchk
```

Save `MULTIWFN_BIN` and the 64-hex `MULTIWFN_SHA256` with the accepted pilot.
Before selecting an analysis, save one unique exact full banner/version line as
`MULTIWFN_VERSION_MARKER` and save the load summary. Require:

- formula and atom count consistent with six carbon and six hydrogen atoms;
- total electron count 42, net charge 0, and expected multiplicity 1;
- alpha and beta electron counts consistent with 21 and 21. If the Multiwfn
  banner does not print multiplicity directly, cross-check these counts against
  the accepted Gaussian charge/multiplicity provenance rather than deriving a
  singlet label from total electron count alone;
- a restricted single-determinant wavefunction unless the source calculation
  intentionally differs;
- orbitals 1 through 21 occupied and orbital 22 unoccupied;
- basis and orbital counts consistent with the formatted checkpoint.

Use the menu displayed by the installed build to inspect orbital energies or
perform a population/bond-order analysis. Record every prompt and response in a
transcript. Do not copy a numeric menu path from another version.

## 4. Validate a selected analysis

- For frontier-orbital output, record orbital 21 and 22 energies with the units
  printed by Multiwfn. Under a declared energy tolerance, identify every
  degenerate partner belonging to the HOMO and LUMO frontier subspaces; orbital
  21 and 22 are occupancy boundaries, not necessarily complete subspaces.
- For atomic partitions or bond orders, name the exact scheme and numerical
  settings. Symmetry-equivalent atoms should agree within a declared tolerance;
  do not use hard-coded charge or bond-order values as an acceptance oracle.
- For a cube output, record the orbital/spin identity, field units, origin,
  axes, spacing, dimensions, atom count, and isosurface convention used by the
  downstream viewer. Require a new, non-empty, parseable output file.
- A UV-Vis or excited-state task requires an input that actually contains the
  relevant transition data. Do not infer it from this ground-state `.fchk`.

Only after this interactive pilot passes may its exact response sequence be
used for a fail-closed batch on the same Multiwfn build and input class.
The batch command must bind both human-readable and byte-exact build evidence:

```bash
: "${MULTIWFN_BIN:?Use the resolved executable path saved by the pilot}"
: "${MULTIWFN_SHA256:?Use the 64-hex executable SHA-256 saved by the pilot}"
: "${MULTIWFN_VERSION_MARKER:?Use the unique exact banner line saved by the pilot}"
BATCH_DIR="${RUN_DIR}_batch_001"
test ! -e "$BATCH_DIR" && test ! -L "$BATCH_DIR" || exit 1
mkdir "$BATCH_DIR"
python3 <skill-dir>/scripts/run_batch.py "$RUN_DIR/benzene.fchk" \
  --workdir "$BATCH_DIR" \
  --responses /absolute/path/to/accepted-benzene-commands.in \
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
  --output benzene-analysis-output.txt
```

## 5. Report provenance

Report Gaussian and Multiwfn versions, all source and output hashes, charge,
multiplicity, electron count, method/basis, exact analysis, menu transcript,
units, validation results, and any unresolved warnings. Mark unexecuted values
`not_computed` rather than presenting the illustrative setup as measured data.
