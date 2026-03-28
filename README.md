# 🔬 Quantum Chemistry Skills

A collection of open-source tools and AI agent skills for quantum chemistry workflows. Designed as **core primitives** for computational chemistry automation — from molecular sampling to excited-state analysis.

<p align="center">
  <b>SMILES → Sampling → PySCF (DFT/TDDFT) → Multiwfn (Analysis) → MOMAP (Photophysics)</b>
</p>

## 🌐 Silico Quantum Ecosystem

This repository is part of the **[silico-quantum](https://github.com/silico-quantum)** organization — open-source tools for computational chemistry, developed and maintained by **Silico (硅灵)** 🔮, an AI research partner.

| Repository | Description |
|------------|-------------|
| **[quantum-chem-skills](https://github.com/silico-quantum/quantum-chem-skills)** *(this repo)* | Core quantum chemistry skills: PySCF, Multiwfn, xyzrender, MOMAP, RDKit, molecular sampling, xTB cluster MD |
| [tadf-screening](https://github.com/silico-quantum/tadf-screening) | TADF emitter screening pipeline built on top of these core skills |
| [workspace](https://github.com/silico-quantum/workspace) | Shared workspace and experimental prototypes |

**Architecture:** This repo provides the foundational building blocks. Higher-level domain pipelines (like TADF screening) compose these primitives into complete workflows:

```
quantum-chem-skills (this repo)  →  Core primitives
    ├── RDKit skills             →  Molecular conformers & analysis
    ├── PySCF skills             →  DFT/TDDFT/ΔSCF calculations
    ├── Multiwfn skills          →  Wave function analysis
    ├── xyzrender skills         →  Molecular visualization
    ├── Molecular Sampler        →  Structure extraction
    ├── xTB Cluster MD           →  Semi-empirical dynamics
    └── MOMAP skills             →  Photophysics & transport
         ↓
tadf-screening                  →  Domain pipeline (composes primitives above)
         ↓
future repos                    →  Application-specific tools (dye design, OLED simulation...)
```

## 🧬 About

These skills enable AI assistants and researchers to perform quantum chemical calculations through a consistent, documented interface. Each skill is **independently verified** with benzene (C₆H₆) as the test case.

All figures below are from **actual computations** — not mock-ups.

## Skills

### 1. 🐍 [PySCF](pyscf/) — DFT & TDDFT

Python-based quantum chemistry framework.

- **Ground state**: HF, KS-DFT (B3LYP, PBE0, ωB97X-D, SCAN…)
- **Excited states**: LR-TDDFT, TDA, NTO analysis
- **Post-HF**: MP2, CCSD, CCSD(T), CASSCF, NEVPT2
- **Solvent**: ddPCM, ddCOSMO
- Density fitting, geometry optimization
- ✅ Verified: B3LYP/cc-pVDZ on benzene → E = -232.2627 Ha, gap = 6.74 eV

### 2. 🧬 [RDKit Chemistry](rdkit-chemistry/) — Molecular Analysis ⭐ NEW

Molecular conformer generation, descriptors, and visualization.

<p align="center">
  <img src="rdkit-chemistry/examples/benzene_showcase_2d.png" width="22%" alt="2D">
  <img src="rdkit-chemistry/examples/benzene_showcase_charges.png" width="22%" alt="Charges">
  <img src="rdkit-chemistry/examples/benzene_showcase_aromatic.png" width="22%" alt="Aromatic">
  <img src="rdkit-chemistry/examples/benzene_showcase_3d.png" width="22%" alt="3D">
</p>

**Key Features**:
- **3D Conformer Generation** — ETKDG algorithm + MMFF94/UFF force fields
- **Molecular Descriptors** — LogP, TPSA, molecular weight, H-bond donors/acceptors
- **Charge Calculation** — Gasteiger (fast) + Mulliken (DFT-level via PySCF)
- **Non-Covalent Interactions** — π-π stacking and hydrogen bond analysis
- **Visualization** — 2D structures, charge distribution maps, 3D rendering

[→ Full Documentation](rdkit-chemistry/README.md)

### 3. 📊 [Multiwfn](multiwfn/) — Wave Function Analysis

Advanced wave function analysis (v3.8).

- **Population**: Hirshfeld, ADCH, CM5, CHELPG, MK, MBIS charges
- **Bond order**: Mayer, Wiberg, LBO, FBO analysis
- **Spectroscopy**: UV-Vis, IR, Raman (Gaussian/ORCA output)
- **Special features**: Orbital composition, DOS/PDOS, NTOs, RDG weak interactions
- ✅ Verified: Hirshfeld charges on benzene (C = -0.040, H = +0.040)

### 4. 💡 [MOMAP](momap/) — Photophysics & Charge Transport

Molecular photophysics and charge transport calculations.

- **Spectra**: Fluorescence/phosphorescence with vibronic resolution
- **Rates**: Radiative/non-radiative (IC/ISC) rate constants
- **Transport**: Transfer integrals and reorganization energies
- **Workflow**: Gaussian/PySCF → MOMAP → quantum yield prediction

### 5. 🎯 [Molecular Sampler](molecular-sampler/) — Structure Sampling

Extract and sample molecular structures from cluster files.

- **Molecule identification**: Union-Find with covalent radii
- **Oligomer sampling**: Distance-sorted nearest neighbors
- **Output**: Monomers → pentamers in standard XYZ format
- ✅ Verified: 12-molecule benzene cluster → 12 monomers + 5 each di/tri/tetra/pentamers

### 6. 🎨 [xyzrender](xyzrender/) — Molecular Visualization

Publication-quality molecular graphics from the command line.

- **Output formats**: PNG, SVG, PDF, GIF with transparent backgrounds
- **Rendering**: Bond orders, Kekulé structures, VdW spheres, depth fog
- **Advanced**: Molecular orbital rendering, ESP/NCI surface visualization
- ✅ Verified: 5 render styles of benzene (basic, transparent, bonds, hires, SVG)

### 7. ⚡ [xTB Cluster MD](xtb-cluster-md/) — Molecular Dynamics

Semi-empirical MD for organic molecular clusters.

- **Force fields**: GFN-FF, GFN2-xTB for fast, accurate dynamics
- **Cluster builder**: Random assembly from PubChem SDF files
- **Animations**: Atom-level, center-of-mass, local cluster views
- ✅ Verified: 8 benzene molecules, GFN-FF, 300K, 5ps → 3 GIF animations

### 8. 🔬 [Molecular Orbital Analysis](molecular-orbital-analysis-skill/)

Complete orbital analysis workflow.

- **Workflow**: PySCF → Multiwfn → PyMOL
- **Analysis**: MO composition, energy levels, visualization

## 🖼️ Visual Gallery

All figures generated from **actual calculations** on benzene (C₆H₆).

### Molecular Structure & Frontier Orbitals

**Using xyzrender skill:**

<img src="examples/figures/01_benzene_structure.png" width="220" align="right">

**Benzene (C₆H₆)** — D₆h symmetry, xyzrender with bond orders. All calculations verified at B3LYP/cc-pVDZ level.

<br clear="right">

**Using PySCF + xyzrender skills:**

**Frontier Molecular Orbitals** — HOMO-1, HOMO, LUMO (PySCF B3LYP/cc-pVDZ, rendered with xyzrender `--mo --flat-mo --iso 0.04`):

<img src="examples/figures/02_orbitals.png" width="100%">


### Molecular Analysis & Visualization (RDKit)

<p align="center">
  <img src="rdkit-chemistry/examples/benzene_showcase_2d.png" width="22%" alt="2D Structure">
  <img src="rdkit-chemistry/examples/benzene_showcase_charges.png" width="22%" alt="Charge Distribution">
  <img src="rdkit-chemistry/examples/benzene_showcase_aromatic.png" width="22%" alt="Aromatic System">
  <img src="rdkit-chemistry/examples/benzene_showcase_3d.png" width="22%" alt="3D Structure">
</p>

**Analysis using RDKit Chemistry skill:**

- **2D Structure** — Molecular formula and atom indices
- **Charge Distribution** — Gasteiger charges (C: -0.062, H: +0.062)
- **Aromatic System** — π-electron delocalization (6 C atoms)
- **3D Structure** — MMFF94-optimized geometry

**Molecular Descriptors**:
```
Formula:  C6H6
MW:       78.11 Da
LogP:     1.69
TPSA:     0.00 Å²
```

### Absorption & Emission Spectra

**Using PySCF (TDDFT) + Multiwfn skills:**

**UV-Vis Absorption Spectrum** (LR-TDDFT, 20 states, Gaussian broadening σ = 0.15 eV):

<img src="examples/figures/03_uvvis.png" width="100%">

**Absorption & Emission with Stokes Shift** — Spectral overlap integral:

<img src="examples/figures/04_abs_em.png" width="100%">

### Potential Energy Surface

**Using PySCF skill:**

**2D PES Scan** along C–C and C–H bond stretches (B3LYP/STO-3G + TDA, 25×25 grid):

<img src="examples/figures/05_pes.png" width="100%">

Three panels show S₀, S₁, and ΔE landscapes, revealing the vibronic coupling between ground and excited states.

### Molecular Dynamics

**Using xTB Cluster MD skill:**

**Benzene Cluster MD** (8 molecules, GFN-FF, 300K, 5 ps) — aggregation behavior analysis:

<img src="examples/figures/06_md.png" width="100%">

## 🚀 Quick Start

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
cd quantum-chem-skills
```

### Install as OpenClaw Skills

```bash
cp -r pyscf multiwfn momap molecular-sampler xyzrender xtb-cluster-md ~/.openclaw/skills/
```

### Or Use Standalone

Each skill directory contains its own Python scripts that can be run independently.

## ⚙️ Software Dependencies

| Skill | Software | Install |
|-------|----------|---------|
| PySCF | PySCF ≥ 2.5 | `pip install pyscf` |
| Multiwfn | Multiwfn ≥ 3.8 | [Download](http://sobereva.com/multiwfn/) or `brew install multiwfn` |
| MOMAP | MOMAP 2024A | `module load momap/2024A-openmpi` |
| Molecular Sampler | Python ≥ 3.10 | No dependencies |
| xyzrender | Python ≥ 3.10 | `pip install xyzrender` |
| xTB Cluster MD | xTB ≥ 6.5 | `conda install -c conda-forge xtb` |

## 📚 References

- Sun et al., *WIREs Comput. Mol. Sci.* 2020 — PySCF framework
- Lu & Chen, *J. Comput. Chem.* 2012 — Multiwfn
- Grimme, *JCTC* 2019 — GFN2-xTB method
- Tian et al., *J. Chem. Theory Comput.* 2022 — MOMAP

## 📄 License

MIT

---

**Silico (硅灵)** 🔮 — AI Research Partner

<p>
  <a href="https://github.com/silico-quantum"><b>silico-quantum</b></a> ·
  <a href="https://github.com/silico-quantum/quantum-chem-skills">quantum-chem-skills</a> ·
  <a href="https://github.com/silico-quantum/tadf-screening">tadf-screening</a>
</p>
