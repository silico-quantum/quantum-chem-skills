# Usage

Each installable skill is a self-contained directory with a `SKILL.md`. After placing a skill in a location supported by your agent client, start a new session or follow that client's reload procedure. See [`INSTALL.md`](INSTALL.md) for OpenAI Codex, Claude Code, OpenClaw, and GitHub Copilot paths.

## Invoke a Skill Through an Agent

Describe both the scientific task and its constraints. For example:

```text
Use the PySCF skill to draft a neutral-singlet B3LYP/def2-SVP calculation
for molecule.xyz. Keep the charge, spin, units, convergence thresholds, and
software assumptions explicit. Do not run the calculation yet.
```

```text
Use the Multiwfn skill to propose a reproducible Hirshfeld charge analysis
for sample.molden. Separate the exact menu inputs from the interpretation.
```

Treat generated commands and scientific interpretations as proposals until their inputs, units, methods, and local software versions have been checked.

## Run Included Helpers Directly

Agent installation is not required to run repository scripts. Install only the dependencies for the helper you use.

### Molecular sampler

Inspect the interface before processing a structure:

```bash
python3 molecular-sampler/molecular_sampler.py --help
```

Example sampling command:

```bash
python3 molecular-sampler/molecular_sampler.py \
  cluster.xyz \
  --xyz-units angstrom \
  --layer all \
  --expected-fragments 12 \
  --samples 20 \
  --output-dir ./samples_run01
```

Replace `12` with an independently justified fragment count. Confirm the input
units, inferred molecular segmentation, manifest hashes, and actual sample
counts before interpreting the generated structures.

### xTB cluster construction and animation

Build a cluster from an SDF input:

```bash
set -e
test ! -e cluster_run01
mkdir cluster_run01
mkdir cluster_run01/build cluster_run01/opt cluster_run01/md
python3 xtb-cluster-md/scripts/build_cluster.py \
  --sdf molecule.sdf \
  --monomer-id "verified monomer identifier" \
  --coordinate-provenance "saved 3D source and preparation record" \
  --molecules 24 --box 40.0 \
  --min-distance 2.0 --max-neighbor-distance 8.0 --seed 42 \
  --out cluster_run01/build/cluster.xyz \
  --manifest cluster_run01/build/cluster_build.json
```

The box and both distance values are starting assumptions, not universal
chemistry settings. Inspect the measured atom pairs, placement connectivity,
and full geometry before optimization.

After a separately accepted optimization, copy only its accepted geometry to
`cluster_run01/md/input.xyz`. Then use the wrapper so the exact argv, exit code,
log hash, unchanged source hashes, and private read-only input snapshots are
retained:

```bash
python3 xtb-cluster-md/scripts/run_xtb_md.py \
  --run-dir cluster_run01/md \
  --input-file input.xyz --input-file md.inp \
  --trajectory xtb.trj --success-marker xtbmdok \
  --log md.log --receipt run_receipt.json -- \
  xtb input.xyz --gfnff --md -I md.inp --chrg 0 --uhf 0
```

Validate the schema-3 receipt, expected method, charge, unpaired-electron
count, exact xTB version, log, marker, settings, input/snapshot identity, every
trajectory frame, and run-specific geometry bounds before visualization; see
[`xtb-cluster-md/SKILL.md`](xtb-cluster-md/SKILL.md) for the complete command.
Then create an animation bound to the accepted report:

```bash
python3 xtb-cluster-md/scripts/make_animation.py \
  --traj cluster_run01/md/xtb.trj \
  --validation cluster_run01/md/validation.json \
  -n 24 --nat-per-mol 12 \
  -o cluster_run01/md/animation.gif \
  --manifest cluster_run01/md/animation.json
```

Use each script's `--help` output as the authority for its current option names.

### Molecular rendering

After installing `xyzrender`:

```bash
set -e
test ! -e molecule_render_run01
mkdir molecule_render_run01
xyzrender molecule.xyz -o molecule_render_run01/molecule.png -t --bo
```

Inspect inferred bonds and bond orders before using a rendering as structural evidence.

### Supported PySCF runner

Inspect the fail-closed closed-shell RKS/TDA runner without starting a
calculation:

```bash
python3 pyscf/scripts/run_safe_dft_tda.py --help
```

Prepare the exact JSON contract described in [`pyscf/SKILL.md`](pyscf/SKILL.md),
then run only into a fresh output directory. Files under `pyscf/tools/`, the old
DFT script, and executable reference examples are quarantined historical
prototypes and intentionally fail closed; they are not supported calculation
entry points.

### Multiwfn

Launch Multiwfn with a supported wavefunction file:

```bash
Multiwfn molecule.molden
```

Menu numbers and available analyses can change between releases. Follow [`multiwfn/SKILL.md`](multiwfn/SKILL.md) together with the documentation for the installed Multiwfn version.

## TADF Stage 4 Adapter

The local [`tadf-screening/stage4_momap.py`](tadf-screening/stage4_momap.py) file is a Stage 4 MOMAP adapter/example. It does not implement the complete candidate-generation, xTB prescreening, and Gaussian workflow.

Inspect its current command-line interface before use:

```bash
python3 tadf-screening/stage4_momap.py --help
```

For the complete pipeline, use the external [`silico-quantum/tadf-screening`](https://github.com/silico-quantum/tadf-screening) repository and verify that its inputs and software assumptions match your environment.

## Saved Examples

[`examples/`](examples/README.md) contains saved benzene logs, cube files, structures, renderings, and animations from a previously recorded environment. Use them as inspection and reproducibility fixtures; they are not automatically regenerated by the repository validation checks.
