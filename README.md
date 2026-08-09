# Quantum Chemistry Skills

Reusable agent skills, reference notes, and small helper scripts for computational chemistry workflows.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Agent Skills](https://img.shields.io/badge/format-Agent%20Skills-blue.svg)](https://agentskills.io/)

This repository is a toolbox, not a single Python package. Each installable skill is a directory containing a `SKILL.md`; some skills also include reference material or executable helpers. Install only the skills you need, then install their scientific software dependencies separately.

## Skills

| Skill | Primary use | External software commonly required |
|---|---|---|
| [`gaussian`](gaussian/) | Gaussian input, job, and analysis guidance | Licensed Gaussian installation |
| [`molecular-orbital-analysis`](molecular-orbital-analysis/) | Indexed orbital cubes and visualization | PySCF; cube viewer; optional Multiwfn |
| [`molecular-sampler`](molecular-sampler/) | Extract and sample molecular assemblies | Python 3.10+ standard library |
| [`momap`](momap/) | MOMAP EVC, spectrum, and ISC workflows | Licensed MOMAP and Gaussian installations |
| [`multiwfn`](multiwfn/) | Wavefunction and real-space analysis | Multiwfn |
| [`pyscf`](pyscf/) | Electronic-structure calculations with PySCF | PySCF, NumPy, SciPy |
| [`rdkit-chemistry`](rdkit-chemistry/) | Molecular construction, conformers, and descriptors | RDKit |
| [`xtb-cluster-md`](xtb-cluster-md/) | Cluster construction and xTB/GFN-FF dynamics | xTB; optional plotting packages |
| [`xyzrender`](xyzrender/) | Molecular and orbital visualization | xyzrender |

The repository also contains [`tadf-screening/`](tadf-screening/README.md), a local Stage 4 MOMAP adapter and example assets. It is not the complete four-stage screening pipeline. See the external [`silico-quantum/tadf-screening`](https://github.com/silico-quantum/tadf-screening) repository for the end-to-end workflow.

## Install a Skill

Clone the repository and choose a skill directory. The following Codex example installs `pyscf` by symbolic link so repository updates are reflected without copying files:

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
cd quantum-chem-skills
mkdir -p /path/to/project/.agents/skills
test ! -e /path/to/project/.agents/skills/pyscf && \
  test ! -L /path/to/project/.agents/skills/pyscf || exit 1
ln -s "$PWD/pyscf" /path/to/project/.agents/skills/pyscf
```

Use a different destination for each agent client:

| Agent client | Project or workspace scope | Personal scope |
|---|---|---|
| [OpenAI Codex](https://developers.openai.com/codex/skills/) | `.agents/skills` | `~/.agents/skills` |
| [Claude Code](https://code.claude.com/docs/en/skills) | `.claude/skills` | `~/.claude/skills` |
| [OpenClaw](https://docs.openclaw.ai/tools/skills) | `<workspace>/skills` or `<workspace>/.agents/skills` | `~/.agents/skills` for the default state, or managed `<state-dir>/skills` (default: `~/.openclaw/skills`) |
| [GitHub Copilot](https://docs.github.com/en/copilot/reference/customization-cheat-sheet) | `.github/skills`, `.agents/skills`, or `.claude/skills` | `~/.copilot/skills` or `~/.agents/skills` |

OpenClaw applies realpath-containment checks to workspace and project skill
roots. Use its local installer or copy the skill into those roots unless the
external source directory is explicitly trusted in OpenClaw configuration.
For cloud-hosted GitHub Copilot agents and code review, copy and commit the
skill directory; a local external symbolic-link target is not included in the
repository checkout.

See [`INSTALL.md`](INSTALL.md) for safe client-specific copy, installer, and
symbolic-link commands, plus optional scientific software setup.

## Use the Repository Directly

The helpers can also be run without installing an agent skill:

```bash
# Inspect the molecular sampler interface; an actual XYZ run must declare units.
python3 molecular-sampler/molecular_sampler.py --help

# Build an initial molecular cluster.
set -e
test ! -e cluster_build_run01
mkdir cluster_build_run01
python3 xtb-cluster-md/scripts/build_cluster.py \
  --sdf molecule.sdf \
  --monomer-id "verified monomer identifier" \
  --coordinate-provenance "saved 3D source and preparation record" \
  --molecules 24 --box 40.0 \
  --min-distance 2.0 --max-neighbor-distance 8.0 --seed 42 \
  --out cluster_build_run01/cluster.xyz \
  --manifest cluster_build_run01/cluster_build.json

# Render an XYZ structure after installing xyzrender.
test ! -e molecule_render_run01
mkdir molecule_render_run01
xyzrender molecule.xyz -o molecule_render_run01/molecule.png -t --bo
```

Commands, input assumptions, and scientific caveats are collected in [`USAGE.md`](USAGE.md) and in each skill's own documentation.

## Repository Layout

```text
quantum-chem-skills/
├── gaussian/                       Gaussian guidance
├── molecular-orbital-analysis/     Orbital analysis workflow
├── molecular-sampler/              Molecular assembly sampler
├── momap/                          MOMAP workflows
├── multiwfn/                       Multiwfn workflows
├── pyscf/                          PySCF workflows and references
├── rdkit-chemistry/                RDKit workflows
├── xtb-cluster-md/                 xTB cluster-MD helpers
├── xyzrender/                      Visualization guidance
├── tadf-screening/                 Stage 4 adapter and examples
├── examples/                       Saved example artifacts
├── scripts/                        Repository validation tools
└── tests/                          Contract tests
```

## Example Artifacts

[`examples/`](examples/README.md) contains saved benzene calculation outputs, renderings, sampled structures, and MD animations. These artifacts document a prior environment and should be treated as reproducibility references, not as claims that every result is regenerated in continuous integration.

## Validation

Run the repository's structural checks without installing scientific packages:

```bash
python3 -m unittest discover -s tests -v
python3 scripts/validate_repository.py .
```

These checks cover repository contracts such as skill metadata, entrypoint
context limits, Markdown links, English-only documentation, Python syntax, and
license completeness. They do not validate scientific accuracy or reproduce
quantum-chemistry results.

## License

Released under the [MIT License](LICENSE).
