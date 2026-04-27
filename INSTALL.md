# Installation

## Prerequisites

- Python ≥ 3.10
- OpenClaw (for agent skill integration)

## Software Setup

```bash
# Create environment
conda create -n qc -y python=3.11
conda activate qc

# Core tools
pip install pyscf numpy scipy
pip install xyzrender
conda install -c conda-forge xtb

# Multiwfn (macOS)
brew install multiwfn
# Or download from http://sobereva.com/multiwfn/

# MOMAP (on HPC with module system)
module load momap/2024A-openmpi
```

## Skill Installation

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
cd quantum-chem-skills
cp -r pyscf multiwfn momap molecular-sampler xyzrender xtb-cluster-md ~/.openclaw/skills/
```

## Verify Installation

```bash
python3 pyscf/references/benzene-dft-tddft.py   # Should complete in ~30s
xyzrender -h                                       # Should show help
xtb --version                                      # Should show version
Multiwfn --version                                 # Should show version
```
