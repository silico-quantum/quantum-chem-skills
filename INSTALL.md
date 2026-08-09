# Installation

There are two independent installation layers:

1. Make one or more skill directories discoverable by your agent client.
2. Install the scientific programs required by those skills.

Installing a skill does not install Gaussian, MOMAP, Multiwfn, PySCF, RDKit, xTB, or other scientific software.

## 1. Clone the Repository

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
cd quantum-chem-skills
```

The installable directories are:

```text
gaussian
molecular-orbital-analysis
molecular-sampler
momap
multiwfn
pyscf
rdkit-chemistry
xtb-cluster-md
xyzrender
```

Every directory in this list contains a `SKILL.md`. The `tadf-screening` directory is an adapter/example directory rather than an installable agent skill.

## 2. Choose Copy or Symbolic Link

A symbolic link is convenient for a development checkout because `git pull` updates the installed skill. A copy is independent of the checkout but must be refreshed manually.

```bash
# Preferred for a development checkout. The destination must not already exist.
mkdir -p /path/to/project/.agents/skills
ln -s "$PWD/pyscf" /path/to/project/.agents/skills/pyscf

# Copy instead. Run this only when the destination does not already contain pyscf.
mkdir -p /path/to/project/.agents/skills
cp -R pyscf /path/to/project/.agents/skills/
```

Replace `pyscf` with another directory from the list above. Check an existing destination before copying or linking so local changes are not merged or overwritten unintentionally.

## 3. Agent-Specific Locations

### OpenAI Codex

[Codex skill documentation](https://developers.openai.com/codex/skills/) defines `.agents/skills` for project skills and `~/.agents/skills` for personal skills.

```bash
# Project-level installation
mkdir -p /path/to/project/.agents/skills
ln -s "$PWD/pyscf" /path/to/project/.agents/skills/pyscf

# Personal installation
mkdir -p ~/.agents/skills
ln -s "$PWD/pyscf" ~/.agents/skills/pyscf
```

### Claude Code

[Claude Code skill documentation](https://code.claude.com/docs/en/skills) defines `.claude/skills` for project skills and `~/.claude/skills` for personal skills.

```bash
# Project-level installation
mkdir -p /path/to/project/.claude/skills
ln -s "$PWD/pyscf" /path/to/project/.claude/skills/pyscf

# Personal installation
mkdir -p ~/.claude/skills
ln -s "$PWD/pyscf" ~/.claude/skills/pyscf
```

### OpenClaw

[OpenClaw skill documentation](https://docs.openclaw.ai/tools/skills) supports a native workspace directory at `<workspace>/skills`, a project-compatible directory at `<workspace>/.agents/skills`, personal skills at `~/.agents/skills` for the default state, and managed skills at `<state-dir>/skills`. The default state directory is `~/.openclaw`.

```bash
# Install the local skill into the active OpenClaw workspace.
openclaw skills install ./pyscf --as pyscf

# Or copy it into an explicit native workspace directory.
mkdir -p /path/to/workspace/skills
cp -R pyscf /path/to/workspace/skills/

# Project-compatible installation by copy
mkdir -p /path/to/workspace/.agents/skills
cp -R pyscf /path/to/workspace/.agents/skills/

# Personal installation for the default state; this root may contain symlinks.
mkdir -p ~/.agents/skills
ln -s "$PWD/pyscf" ~/.agents/skills/pyscf

# Managed installation in the default state directory
mkdir -p ~/.openclaw/skills
ln -s "$PWD/pyscf" ~/.openclaw/skills/pyscf
```

OpenClaw requires workspace, project-agent, and extra-directory skill targets
to resolve inside the configured root. An external symbolic link works only
when its target root is explicitly listed in `skills.load.allowSymlinkTargets`.
The local installer copies into the active workspace by default. If
`OPENCLAW_STATE_DIR` selects a custom state, use that `<state-dir>/skills`
directory instead of the default `~/.openclaw/skills`; home-scoped personal
skills may not be visible in a custom state.

### GitHub Copilot

[GitHub Copilot skill documentation](https://docs.github.com/en/copilot/customizing-copilot/extending-copilot-chat-with-skills) supports `.github/skills`, `.agents/skills`, and `.claude/skills` at project scope. It supports `~/.copilot/skills` and `~/.agents/skills` at personal scope.

```bash
# Repository-contained installation for cloud-hosted Copilot agents and review
mkdir -p /path/to/project/.github/skills
cp -R pyscf /path/to/project/.github/skills/

# A project symlink is suitable only for local IDE or CLI sessions.
mkdir -p /path/to/project/.agents/skills
ln -s "$PWD/pyscf" /path/to/project/.agents/skills/pyscf

# Personal installation using the Copilot-specific location
mkdir -p ~/.copilot/skills
ln -s "$PWD/pyscf" ~/.copilot/skills/pyscf

# Personal installation using the shared Agent Skills location
mkdir -p ~/.agents/skills
ln -s "$PWD/pyscf" ~/.agents/skills/pyscf
```

Commit the copied skill directory when a cloud-hosted Copilot agent or code
review must load it. External symbolic-link targets exist only on the local
machine and are not part of a remote repository checkout. The same distinction
applies when choosing `.agents/skills` or `.claude/skills` at project scope.

Start a new agent session, or use the client's documented reload behavior, after adding or updating a skill.

## 4. Install Scientific Dependencies Selectively

Follow each skill's documentation and the upstream software documentation. A single environment is not required for every skill.

Examples for the open-source Python helpers are:

```bash
conda create -n quantum-chem-skills -y python=3.11
conda activate quantum-chem-skills

# PySCF references and the molecular sampler
python -m pip install pyscf numpy scipy

# RDKit and xTB, when needed
conda install -c conda-forge rdkit xtb

# Animation helpers, when needed
python -m pip install matplotlib pillow

# Rendering skill, when needed
python -m pip install xyzrender
```

Gaussian and MOMAP are separately distributed scientific programs and may require a license or site-specific module setup. Multiwfn and PyMOL should likewise be installed from their upstream distribution channels. Do not copy machine-specific login commands or module names from another environment without checking your local installation.

## 5. Verify the Installation

First verify that the selected skill is visible at the intended path:

```bash
test -f /path/to/project/.agents/skills/pyscf/SKILL.md
```

Then run only the relevant software checks:

```bash
python3 molecular-sampler/molecular_sampler.py --help
python3 xtb-cluster-md/scripts/build_cluster.py --help
python -c "import pyscf; print(pyscf.__version__)"
xtb --version
```

Licensed or site-managed programs should be checked with the command provided by the local administrator. A version command or import check confirms availability only; it does not validate a scientific workflow.

To validate the repository itself:

```bash
python3 -m unittest discover -s tests -v
python3 scripts/validate_repository.py .
```
