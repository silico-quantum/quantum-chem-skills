---
name: xtb-cluster-md
description: Use when building non-overlapping molecular clusters, running explicit-charge xTB molecular dynamics, validating trajectories, or rendering cluster-motion animations.
license: MIT
compatibility: Requires POSIX process-group support, xTB, and Python 3.10+; cluster building and validation use the standard library, while animation helpers require NumPy, Matplotlib, and Pillow.
---

# xTB cluster molecular dynamics

Build a deterministic molecular cluster, run xTB MD in a clean directory, and
accept the result only after log, success-marker, trajectory, and requested-time
checks pass. A visually plausible animation is not evidence of a successful MD
run.

## Prerequisites

1. Install xTB from an upstream release or a trusted package and capture the
   actual binary and version:

   ```bash
   command -v xtb
   xtb --version
   python3 --version
   ```

2. Locate this skill directory. Its cluster builder and validator require only
   Python's standard library. Animation helpers additionally require NumPy,
   Matplotlib, and Pillow.
3. Provide a verified three-dimensional monomer structure. Network lookup is
   optional; if used, save the source URL, retrieval date, and file hash.
4. Estimate storage and compute cost before production MD. Start with a short
   smoke run using the same method, charge, spin state, temperature, and input
   preparation as production.

## Input contract

Declare and record:

- exactly one complete 3D V2000 SDF record in Angstrom, including its bond
  block, `M  END`, `$$$$`, explicit molecular identity, and coordinate
  provenance;
- molecule count, cubic placement-box length in Angstrom, random seed, and a
  physically justified minimum intermolecular atom-atom distance cutoff plus
  a maximum atom-pair neighbor distance that keeps the placement graph
  connected;
- total-system **charge and unpaired electrons** as integers; xTB `--uhf` is
  `Nalpha-Nbeta`, which equals the number of unpaired electrons in the usual
  high-spin convention;
- xTB method and version, temperature in K, total time in ps, propagation
  **time step** in fs, and trajectory dump interval in fs;
- atom count per molecule and total atom count;
- optimization/equilibration decisions, ensemble, SHAKE setting, hydrogen-mass
  setting, SCC accuracy, thread count, and solvent/environment model; and
- new `build`, `opt`, and `md` directories with no stale xTB files.

The bundled builder centers the monomer on its arithmetic coordinate centroid,
not its center of mass. Placement acceptance uses the true minimum atom-atom
distance between different molecules; centroid separation is not an overlap
test.

GFN-FF is non-electronic. Record the intended unpaired-electron count even when
the selected force-field model cannot represent that electronic state. If spin
physics matters to the conclusion, use an appropriate electronic GFN-xTB method
and justify it rather than treating GFN-FF as spin-resolved.

## Workflow

### 1. Create a fresh run layout

Choose a new project path and fail if it already exists:

```bash
set -e
test ! -e cluster_run_001 && test ! -L cluster_run_001 || exit 1
mkdir cluster_run_001
mkdir cluster_run_001/build cluster_run_001/opt cluster_run_001/md
```

Copy the verified monomer into `build/`. Do not run in a directory containing
`.CHRG`, `.UHF`, `gfnff_topo`, `xtbrestart`, `mdrestart`, `xtb.trj`, or earlier
logs because xTB may read or overwrite conventional filenames.

```bash
VERIFIED_MONOMER=/absolute/path/to/verified-monomer.sdf
test -s "$VERIFIED_MONOMER"
test ! -e cluster_run_001/build/monomer.sdf && \
  test ! -L cluster_run_001/build/monomer.sdf || exit 1
cp "$VERIFIED_MONOMER" cluster_run_001/build/monomer.sdf
```

### 2. Build and validate the initial cluster

Replace `<skill-dir>` with the resolved skill path:

```bash
python3 <skill-dir>/scripts/build_cluster.py \
  --sdf cluster_run_001/build/monomer.sdf \
  --monomer-id "benzene; verified neutral structure" \
  --coordinate-provenance "saved 3D source record and preparation method" \
  --molecules 8 --box 24.0 \
  --min-distance 2.0 --max-neighbor-distance 8.0 --seed 42 \
  --out cluster_run_001/build/cluster_N8.xyz \
  --manifest cluster_run_001/build/cluster_build.json
```

Both distance thresholds are examples and must be chosen for the chemistry.
The script rejects multiple/truncated SDF records, refuses existing outputs,
preserves a deterministic seed, requires every molecule after the first to
connect to an earlier molecule within the declared neighbor threshold, and
records the closest intramolecular/intermolecular atom pairs. The result is a
connected placement ensemble, not an equilibrated condensed-phase structure.
The XYZ and manifest are published from sibling partial files; acceptance also
requires the manifest's `RUN_INCOMPLETE` marker to be absent. Inspect the full
structure and manifest before calculating.

### 3. Optimize in a separate directory

Copy the accepted cluster into the new `opt/` directory. Run the chosen method
with explicit total charge and `--uhf`, even when both are zero:

```bash
(
cd cluster_run_001/opt
for target in input.xyz opt.log; do
  test ! -e "$target" && test ! -L "$target" || exit 1
done
cp ../build/cluster_N8.xyz input.xyz
xtb input.xyz --gfnff --opt --chrg 0 --uhf 0 --parallel 4 > opt.log 2>&1
)
```

Before continuing, require a zero shell exit status, a non-empty optimized
geometry with the same atom count/order, and `normal termination of xtb` in
`opt.log`. Reject an optimization with fatal errors, non-convergence, or
unphysical fragmentation/contacts even if an output geometry exists.

Copy only the accepted optimized geometry into the fresh `md/` directory. Do
not copy `gfnff_topo` after a material geometry or method change; upstream xTB
documentation warns that stale topology data must be removed after significant
changes.

### 4. Define the MD settings

For a GFN-FF run, an explicit example is:

```bash
test ! -e cluster_run_001/md/md.inp && \
  test ! -L cluster_run_001/md/md.inp || exit 1
```

Only after that guard succeeds, save this content as
`cluster_run_001/md/md.inp`:

```text
$md
   temp=300.0
   time=5.0
   step=2.0
   dump=50.0
   nvt=true
   hmass=4
   shake=0
   sccacc=2.0
$end
```

Units are K, ps, fs, and fs, respectively. The official GFN-FF documentation
recommends `step=2.0`, `hmass=4.0`, and `shake=0` instead of the generic MD
defaults. Do not copy those choices to another method or system without review.
Ensure that `dump / step` gives an intentional integer-like output cadence.

### 5. Run MD with explicit state

Run through the bundled no-shell wrapper so the exact argv, working directory,
binary hash, immutable private input snapshots, source hashes before and after
execution, log hash, return code, signal, elapsed time, and post-run
trajectory/success-marker hashes and timestamps are written to a fresh receipt.
The child process receives the read-only snapshots, not the mutable source
paths. The runner rejects a leader that exits while any descendant remains,
terminates the whole process group, and confirms descendant exit before hashing
the log or generated artifacts:

```bash
for target in input.xyz md.log run_receipt.json xtb.trj xtbmdok; do
  test ! -e "cluster_run_001/md/$target" && \
    test ! -L "cluster_run_001/md/$target" || exit 1
done
cp cluster_run_001/opt/xtbopt.xyz cluster_run_001/md/input.xyz
python3 <skill-dir>/scripts/run_xtb_md.py \
  --run-dir cluster_run_001/md \
  --input-file input.xyz --input-file md.inp \
  --trajectory xtb.trj --success-marker xtbmdok \
  --log md.log --receipt run_receipt.json -- \
  xtb input.xyz --gfnff --md -I md.inp \
    --chrg 0 --uhf 0 --parallel 4
```

The wrapper refuses pre-existing `xtb.trj` and `xtbmdok`; custom artifact names
must be declared explicitly. Do not accept the run merely because a trajectory
exists or because the log contains a normal-looking final energy. A
nonzero/unknown receipt return code, timeout, signal, or receipt/log/artifact
hash mismatch is a hard rejection.

### 6. Validate before visualization

Use the total atom count and the exact requested MD values:

```bash
: "${XTB_VERSION:?Set XTB_VERSION to the exact version recorded in md.log}"
test ! -e cluster_run_001/md/validation.json && \
  test ! -L cluster_run_001/md/validation.json || exit 1
python3 <skill-dir>/scripts/validate_md_run.py \
  --run-dir cluster_run_001/md \
  --receipt run_receipt.json --log md.log --trajectory xtb.trj \
  --success-marker xtbmdok \
  --input-xyz input.xyz --md-input md.inp \
  --expected-molecules 8 --atoms-per-molecule 12 \
  --expected-method gfnff \
  --expected-charge 0 --expected-uhf 0 \
  --expected-xtb-version "$XTB_VERSION" \
  --requested-time-ps 5.0 \
  --requested-step-fs 2.0 \
  --requested-dump-fs 50.0 \
  --requested-temperature-k 300.0 \
  --requested-scc-accuracy 2.0 \
  --requested-hydrogen-mass 4.0 \
  --requested-shake-bonds 0 \
  --expect-thermostat on \
  --minimum-intermolecular-distance 0.7 \
  --maximum-system-extent 100.0 \
  --report validation.json
```

Set `XTB_VERSION` to the exact version token reported by the accepted MD log,
for example `6.7.1`; do not infer it from a different installation. The two
geometry thresholds are explicit run-specific rejection limits, not universal
defaults. The validator requires a schema-3 receipt with a resolved executable
hash, zero return code, unchanged current input hashes, matching private
read-only argv snapshots, a current log hash, and fresh bound trajectory and
success-marker records. It also binds the method, charge, unpaired-electron
count, and xTB version and requires one ordered
calculation/MD setup, one `normal exit of md()`, one
`normal termination of xtb`, a non-empty `xtbmdok`, unique matching logged
settings, and a complete trajectory. Every frame must preserve input labels and
monomer blocks, remain finite, exceed the declared intermolecular collision
threshold, and stay inside the declared extent limit. It rejects instability,
thermostat problems, fatal errors, runtime exceptions, and SCC non-convergence.
Coverage is derived from frame count and the unique logged dump interval with
one dump interval of tolerance; no absent timestamp is invented.

### 7. Render only an accepted trajectory

After `validation.json` reports `accepted`, choose one visualization:

```bash
# All atoms, colored by molecule index.
for target in cluster_atoms.gif cluster_atoms.json; do
  test ! -e "cluster_run_001/md/$target" && \
    test ! -L "cluster_run_001/md/$target" || exit 1
done
python3 <skill-dir>/scripts/make_atom_animation.py \
  --traj cluster_run_001/md/xtb.trj \
  --validation cluster_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --stride 5 \
  --title "Cluster MD, accepted run 001" \
  --out cluster_run_001/md/cluster_atoms.gif \
  --manifest cluster_run_001/md/cluster_atoms.json

# Arithmetic coordinate centroids, not centers of mass.
for target in cluster_centroids.gif cluster_centroids.json; do
  test ! -e "cluster_run_001/md/$target" && \
    test ! -L "cluster_run_001/md/$target" || exit 1
done
python3 <skill-dir>/scripts/make_animation.py \
  --traj cluster_run_001/md/xtb.trj \
  --validation cluster_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --stride 5 \
  --title "Molecular coordinate centroids, run 001" \
  --out cluster_run_001/md/cluster_centroids.gif \
  --manifest cluster_run_001/md/cluster_centroids.json

# A last-frame local subset selected by minimum atom-atom distance.
for target in cluster_local.gif cluster_local.json; do
  test ! -e "cluster_run_001/md/$target" && \
    test ! -L "cluster_run_001/md/$target" || exit 1
done
python3 <skill-dir>/scripts/make_local_cluster_animation.py \
  --traj cluster_run_001/md/xtb.trj \
  --validation cluster_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --cluster-size 4 \
  --distance minatom --stride 2 --max-frames 100 \
  --title "Local cluster, accepted run 001" \
  --out cluster_run_001/md/cluster_local.gif \
  --manifest cluster_run_001/md/cluster_local.json
```

Each helper requires an accepted validation report whose trajectory hash and
counts match, scans the entire trajectory before sampling, writes the GIF via a
temporary file, reopens it, and emits a rendering manifest. When drawing bonds
in the local helper, also pass `--bonds --build-manifest <build-manifest>`; only
the explicit source-SDF topology is allowed. Record downsampling and selection.
A selected local cluster is a visualization, not an unbiased aggregation
statistic.

## Validation and acceptance

Accept the workflow only if all applicable checks pass:

1. The one-record SDF, monomer identity, 3D coordinate provenance, source hash,
   atom/bond counts, explicit topology, and Angstrom units are documented.
2. The cluster XYZ contains exactly `molecule_count * atoms_per_molecule`
   finite records, its manifest hash matches, and no incomplete marker remains.
3. The measured minimum intermolecular atom-atom distance meets the declared
   cutoff, the placement graph meets the maximum neighbor threshold, and
   closest atom/molecule pairs are recorded; centroid distance is not
   substituted.
4. Total charge, unpaired electrons, method, xTB version, temperature, time,
   time step, dump interval, and environment are explicit.
5. Optimization and MD use fresh directories and have preserved command/log
   provenance.
6. The immutable schema-3 receipt records successful process completion and
   binds the exact MD argv, private read-only input snapshots, unchanged source
   inputs, run directory, return code, executable, and validated log hash.
7. The non-empty `xtbmdok` exists; the log has unique calculation/MD sections,
   MD-specific and program termination markers, and no rejection diagnostic.
8. Logged time, step, dump, temperature, SCC, hydrogen mass, SHAKE, thermostat,
   and derived maximum steps match the requested values.
9. Every trajectory frame is complete, matches input/monomer atom ordering,
   contains finite coordinates, and passes declared collision/extent limits.
10. The final frame count covers the requested time within the documented xTB
   dump-frequency tolerance.
11. Any claimed aggregation metric is computed over an explicitly defined
    trajectory window and is not inferred from a single chosen animation.

## Failure handling

| Failure | Required response |
|---|---|
| SDF is 2D, truncated, or lacks provenance | Stop and obtain a verified 3D monomer; do not let xTB conversion silently redefine the input. |
| SDF has multiple records, an invalid bond block, or no exact `M  END`/`$$$$` boundary | Reject it; select and save one authoritative monomer record explicitly. |
| Cluster placement exhausts attempts | Increase the box, review the cutoff, or reduce molecule count in a new run; never bypass atom-pair checks. |
| Placement cannot satisfy the maximum neighbor distance | Review whether the requested connected ensemble and box are physically sensible; never relabel an isolated gas ensemble as a connected cluster. |
| Intermolecular cutoff fails after writing | Reject the XYZ and investigate the builder/settings. |
| Charge or unpaired-electron state is unknown | Stop and obtain the intended total state; do not assume neutral singlet. |
| Optimization fails or creates implausible contacts | Reject it, revise preparation/method, and preserve the failed log. |
| Receipt is missing, nonzero, interrupted, timed out, or has a mismatched log hash | Reject the run regardless of trajectory presence. |
| MD lacks non-empty `xtbmdok`, `normal exit of md()`, or program termination | Reject the run regardless of trajectory presence. |
| Log reports instability, thermostat problem, runtime exception, fatal error, or SCC failure | Reject the run even if xTB later prints a termination line. Reduce the time step or repair the model only in a new documented run. |
| Trajectory is truncated, atom ordering changes, collides, or exceeds the declared extent | Reject it; do not animate a partial or unphysical frame set as a complete simulation. |
| Frame coverage is shorter than requested | Report the achieved coverage and mark production MD incomplete. |
| Animation helper receives no valid frames or would overwrite output | Stop, fix inputs or choose a new output path, and rerun after validation. |

## Output and reporting

Report:

- monomer identity/source/hash, coordinate provenance, format, units, atom/bond
  counts, topology, and closest intramolecular pair;
- cluster molecule/atom counts, seed, centroid-translation box, requested lower
  cutoff/upper neighbor threshold, measured nearest distances and atom pairs;
- total charge and unpaired electrons, xTB version/method, exact command,
  environment, thread count, and input hashes;
- optimization acceptance and optimized-geometry hash;
- MD temperature, ensemble, total time, time step, dump interval, SHAKE,
  hydrogen mass, SCC settings, exact argv, executable/receipt/log hashes, process
  status, exit status, signal, and elapsed time;
- success-marker and normal-termination checks plus every rejected warning;
- trajectory hash, frame count, atoms per frame, distinct-frame count, minimum
  intermolecular distance/pair, maximum extent, derived coverage, and final
  frame/comment; and
- visualization type, selected molecules/window, stride, frame count, title,
  dimensions, dependency versions, output path/hash, and limitations.

Keep requested settings, observed outputs, validation results, and scientific
interpretation separate. Do not report a nominal total time as achieved when
the accepted trajectory does not cover it.

## References

- [xTB command-line documentation](https://xtb-docs.readthedocs.io/en/latest/commandline.html)
  — `--chrg`, `--uhf`, method, input, and parallel options.
- [xTB molecular-dynamics documentation](https://xtb-docs.readthedocs.io/en/latest/md.html)
  — MD units, trajectory cadence, `xtb.trj`, `mdrestart`, and `xtbmdok`.
- [xTB GFN-FF documentation](https://xtb-docs.readthedocs.io/en/latest/gfnff.html)
  — GFN-FF scope, MD-specific settings, and topology-file warning.
- [xTB upstream repository](https://github.com/grimme-lab/xtb) — source,
  releases, and citation information.
- [Benzene reference audit](references/benzene-cluster-example.md) — historical
  bundled artifacts and the current acceptance boundary.
