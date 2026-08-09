# Benzene cluster reference audit

This page is a reproducible template and an audit of the historical benzene
artifacts under `examples/xtb-cluster-md/`. It is not a claim that a current
xTB installation reproduced those artifacts.

## Historical artifact status

The bundled example includes an initial 96-atom XYZ, two copies of an xTB log,
and three GIF files. It does not include the original `xtb.trj` or the
`xtbmdok` success marker required by the current acceptance contract. The log
contains both `thermostating problem` and `Runtime exception occurred`, so the
current validator rejects it even though it also contains a normal-termination
line. The GIFs are therefore historical visual artifacts, not accepted
trajectory evidence.

The old wall-time observation is also not a portable benchmark because it lacks
enough retained hardware and environment metadata. Measure a new run locally.

## Reproducible input template

### 1. Obtain and identify the monomer

Use an authoritative 3D benzene record and preserve its URL, retrieval date,
identifier, and SHA-256. For example, PubChem PUG REST can provide an SDF, but
network results must still be validated before use:

```bash
test ! -e benzene_3d.sdf && test ! -L benzene_3d.sdf || exit 1
curl --fail --location --remove-on-error --output benzene_3d.sdf \
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/241/SDF?record_type=3d"
shasum -a 256 benzene_3d.sdf
```

Confirm CID 241, formula C6H6, one V2000 record, 12 atoms, and a source record
identified as three-dimensional. A planar molecule may legitimately have one
constant coordinate, so dimensional provenance cannot be inferred from a
nonzero-z heuristic alone.

### 2. Build eight molecules

From the repository root, replace `<skill-dir>` with `xtb-cluster-md` or its
resolved installed path:

```bash
set -e
test ! -e benzene_run_001 && test ! -L benzene_run_001 || exit 1
mkdir benzene_run_001
mkdir benzene_run_001/build benzene_run_001/opt benzene_run_001/md
test ! -e benzene_run_001/build/monomer.sdf && \
  test ! -L benzene_run_001/build/monomer.sdf || exit 1
cp benzene_3d.sdf benzene_run_001/build/monomer.sdf

python3 <skill-dir>/scripts/build_cluster.py \
  --sdf benzene_run_001/build/monomer.sdf \
  --monomer-id "benzene; PubChem CID 241" \
  --coordinate-provenance "PubChem 3D SDF record; saved URL/date/hash" \
  --molecules 8 --box 24.0 \
  --min-distance 2.0 --max-neighbor-distance 8.0 --seed 42 \
  --out benzene_run_001/build/benzene_N8.xyz \
  --manifest benzene_run_001/build/cluster_build.json
```

The 2.0 Angstrom lower cutoff, 8.0 Angstrom neighbor threshold, and 24 Angstrom
centroid-translation box are template values, not validated universal settings.
Inspect the measured atom pairs, connected placement, and full geometry.

### 3. Optimize with explicit total state

Neutral benzene N8 has total charge 0 and zero unpaired electrons in this
template. Preserve those explicit values:

```bash
(
cd benzene_run_001/opt
for target in input.xyz opt.log; do
  test ! -e "$target" && test ! -L "$target" || exit 1
done
cp ../build/benzene_N8.xyz input.xyz
xtb input.xyz --gfnff --opt --chrg 0 --uhf 0 --parallel 4 > opt.log 2>&1
)
```

Continue only after the optimized 96-atom geometry and log pass the acceptance
criteria in the parent skill.

### 4. Run a short MD smoke test

First refuse an existing path or dangling symbolic link:

```bash
test ! -e benzene_run_001/md/md.inp && \
  test ! -L benzene_run_001/md/md.inp || exit 1
```

Only after that guard succeeds, create `benzene_run_001/md/md.inp` with:

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

Then run through the receipt-producing wrapper from the repository root:

```bash
for target in input.xyz md.log run_receipt.json xtb.trj xtbmdok; do
  test ! -e "benzene_run_001/md/$target" && \
    test ! -L "benzene_run_001/md/$target" || exit 1
done
cp benzene_run_001/opt/xtbopt.xyz benzene_run_001/md/input.xyz
python3 <skill-dir>/scripts/run_xtb_md.py \
  --run-dir benzene_run_001/md \
  --input-file input.xyz --input-file md.inp \
  --trajectory xtb.trj --success-marker xtbmdok \
  --log md.log --receipt run_receipt.json -- \
  xtb input.xyz --gfnff --md -I md.inp \
    --chrg 0 --uhf 0 --parallel 4
```

The requested 5 ps time, 2 fs step, and 50 fs dump interval are not evidence
that 5 ps completed. The wrapper must first reject and terminate any surviving
process-group descendant and confirm exit before its artifact hashes are
stable. Set `XTB_VERSION` to the exact version token recorded in the MD log,
then validate the actual command state, log, and trajectory:

```bash
: "${XTB_VERSION:?Set XTB_VERSION to the exact version recorded in md.log}"
test ! -e benzene_run_001/md/validation.json && \
  test ! -L benzene_run_001/md/validation.json || exit 1
python3 <skill-dir>/scripts/validate_md_run.py \
  --run-dir benzene_run_001/md \
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

The 0.7 and 100.0 Angstrom limits are explicit example rejection bounds. Review
and justify them for the intended system before use. The validator also binds
the schema-3 receipt to private read-only input snapshots and rejects changed
source inputs, a method/state mismatch, or a different xTB version.

### 5. Render accepted output

Only after `validation.json` reports `accepted`:

```bash
for target in benzene_atoms.gif benzene_atoms.json \
  benzene_centroids.gif benzene_centroids.json \
  benzene_local.gif benzene_local.json; do
  test ! -e "benzene_run_001/md/$target" && \
    test ! -L "benzene_run_001/md/$target" || exit 1
done

python3 <skill-dir>/scripts/make_atom_animation.py \
  --traj benzene_run_001/md/xtb.trj \
  --validation benzene_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --stride 5 \
  --title "Benzene N8 atoms, accepted run 001" \
  --out benzene_run_001/md/benzene_atoms.gif \
  --manifest benzene_run_001/md/benzene_atoms.json

python3 <skill-dir>/scripts/make_animation.py \
  --traj benzene_run_001/md/xtb.trj \
  --validation benzene_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --stride 5 \
  --title "Benzene N8 coordinate centroids, accepted run 001" \
  --out benzene_run_001/md/benzene_centroids.gif \
  --manifest benzene_run_001/md/benzene_centroids.json

python3 <skill-dir>/scripts/make_local_cluster_animation.py \
  --traj benzene_run_001/md/xtb.trj \
  --validation benzene_run_001/md/validation.json \
  --molecules 8 --nat-per-mol 12 --cluster-size 4 \
  --distance minatom --stride 2 --max-frames 100 \
  --title "Benzene local subset, accepted run 001" \
  --out benzene_run_001/md/benzene_local.gif \
  --manifest benzene_run_001/md/benzene_local.json
```

The local subset is selected from last-frame intermolecular atom distances. It
must not be presented as an unbiased cluster-size distribution or kinetic
analysis.
