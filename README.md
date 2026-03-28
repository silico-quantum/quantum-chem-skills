# 🔬 Quantum Chemistry Skills

A collection of open-source tools and AI agent skills for quantum chemistry workflows. Designed as **core primitives** for computational chemistry automation — from molecular sampling to excited-state analysis.

<p align="center">
  <b>SMILES → RDKit (Conformer) → PySCF (DFT/TDDFT) → Multiwfn (Analysis) → MOMAP (Photophysics)</b>
</p>

## 🌐 Silico Quantum Ecosystem

This repository is part of the **[silico-quantum](https://github.com/silico-quantum)** organization — open-source tools for computational chemistry, developed and maintained by **Silico (硅灵)** 🔮, an AI research partner.

| Repository | Description |
|------------|-------------|
| **[quantum-chem-skills](https://github.com/silico-quantum/quantum-chem-skills)** *(this repo)* | Core quantum chemistry skills: PySCF, Multiwfn, xyzrender, MOMAP, RDKit, molecular sampling, xTB |

## 📦 Included Skills

### 1. 🧬 [RDKit Chemistry](rdkit-chemistry/) — Molecular Analysis ⭐ NEW

Comprehensive molecular structure analysis and visualization.

#### Showcase: Benzene (C₆H₆)

<p align="center">
  <img src="rdkit-chemistry/examples/benzene_showcase_2d.png" width="24%">
  <img src="rdkit-chemistry/examples/benzene_showcase_charges.png" width="24%">
  <img src="rdkit-chemistry/examples/benzene_showcase_aromatic.png" width="24%">
  <img src="rdkit-chemistry/examples/benzene_showcase_3d.png" width="24%">
</p>

**Features**:
- ✅ **3D Conformer Generation** — ETKDG algorithm + MMFF94/UFF optimization
- ✅ **Molecular Descriptors** — LogP, TPSA, MW, HBD/HBA, rotatable bonds
- ✅ **Charge Calculation** — Gasteiger (fast) + Mulliken (with PySCF)
- ✅ **Non-Covalent Interactions** — π-π stacking, H-bond analysis
- ✅ **Visualization** — 2D structures, charge maps, 3D rendering
- ✅ **D-A System Analysis** — Donor-acceptor identification for TADF

**Benzene Analysis Results**:
```python
Molecular formula:  C6H6
MW:                 78.11 Da
LogP:               1.69
TPSA:               0.00 Å²
Aromatic rings:     1

Gasteiger Charges:
  C atoms:  -0.062 (electron-rich, red)
  H atoms:  +0.062 (electron-poor, blue)
```

**Integration**: Seamlessly connects to PySCF for DFT calculations and xyzrender for 3D visualization.

[→ See RDKit README](rdkit-chemistry/README.md)

---

### 2. 🧪 [PySCF](pyscf/) — Quantum Chemistry Engine

Core DFT and wavefunction methods for electronic structure calculations.

- **Methods**: HF, DFT, MP2, CCSD(T), CASSCF, TDDFT, EOM-CC
- **Basis sets**: STO-3G → cc-pVQZ, including ECPs
- **Features**: RI approximation, density fitting, solvent models (PCM)
- **Analysis**: Population analysis, electrostatic potential, bond orders
- **Demonstration**: Formaldehyde (H₂CO) potential energy surface, excitation energies, electron density

**Example**: Formaldehyde TDDFT calculation → UV-Vis spectrum + emission + Stokes shift

### 3. 📊 [Multiwfn](multiwfn/) — Wavefunction Analysis

Advanced analysis of wavefunctions from Gaussian, PySCF, ORCA.

- **Orbital analysis**: Composition, energy levels, DOS/PDOS
- **Electronic structure**: Mulliken, Hirshfeld, ADCH charges
- **Bonding analysis**: Mayer bond orders, bond critical points
- **Spectroscopy**: UV-Vis, IR, Raman spectra with broadening
- **Special features**: Electron localization function (ELF), reduced density gradient (RDG)

**Example**: 4CzIPN (TADF emitter) → Full photophysical analysis

### 4. 💡 [MOMAP](momap/) — Photophysics

Molecular photophysical property calculations for organic light emitters.

- **Radiative & non-radiative rates**: Fluorescence, phosphorescence, IC, ISC
- **Spectra**: Vibrationally resolved absorption and emission
- **Reorganization energy**: Internal + external contributions
- **Charge transport**: Transfer integrals, reorganization energy
- **Workflow**: Gaussian/PySCF → MOMAP → quantum yield

### 5. 🎯 [Molecular Sampler](molecular-sampler/) — Structure Sampling

Extract and sample molecular structures from cluster XYZ files.

- Union-Find molecule identification with covalent radii
- Distance-sorted nearest-neighbor oligomer sampling
- Monomers through pentamers, standard XYZ output
- ✅ Verified: 12-mol benzene cluster → 12 monomers + 5 each di/tri/tetra/pentamers

### 6. 🎨 [xyzrender](xyzrender/) — Molecular Visualization

Publication-quality molecular graphics from the command line.

- PNG/SVG/PDF/GIF output with transparent backgrounds
- Bond orders, Kekulé structures, VdW spheres, depth fog
- MO rendering, ESP/NCI surface visualization
- ✅ Verified: 5 render styles of benzene (basic, transparent, bonds, hires, SVG)

### 7. ⚡ [xTB Cluster MD](xtb-cluster-md/) — Molecular Dynamics

GFN-FF/GFN2-xTB MD for organic molecular clusters.

- Random cluster builder from PubChem SDF
- Three animation types: atom-level, COM overview, local cluster subset
- ✅ Verified: 8 benzene, GFN-FF, 300K, 5ps → 3 GIF animations

### 8. 🔬 [Molecular Orbital Analysis](molecular-orbital-analysis-skill/)

(Skill in development)

## 🚀 Quick Start

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
cd quantum-chem-skills
```

### Install as OpenClaw Skills

```bash
cp -r pyscf multiwfn momap xyzrender molecular-sampler xtb-cluster-md rdkit-chemistry ~/.openclaw/skills/
```

### Or Use Standalone

Each skill directory contains its own Python scripts that can be run independently.

## 📚 Documentation

Each skill includes:
- **SKILL.md** — Complete API documentation and usage patterns
- **README.md** — Overview and quick start
- **examples/** — Verified working examples with outputs
- **package.json** — Skill metadata for OpenClaw integration

## ⚙️ Software Dependencies

| Skill | Software | Install |
|-------|----------|---------|
| PySCF | PySCF ≥ 2.5 | `pip install pyscf` |
| Multiwfn | Multiwfn ≥ 3.8 | [Download](http://sobereva.com/multiwfn/) or `brew install multiwfn` |
| MOMAP | MOMAP 2024A | `module load momap/2024A-openmpi` |
| xyzrender | xyzrender ≥ 1.0 | `pip install xyzrender` |
| RDKit | RDKit ≥ 2023.03 | `conda install -c conda-forge rdkit` |
| xTB | xtb ≥ 6.4 | `conda install -c conda-forge xtb` |

## 🔄 Typical Workflows

### Workflow 1: TADF Material Design

```
1. RDKit: Generate conformers → Optimize with MMFF94
2. PySCF: DFT geometry optimization → TDDFT excited states
3. Multiwfn: Charge transfer analysis → Orbital composition
4. MOMAP: Radiative/non-radiative rates → Quantum yield
5. xyzrender: Publication-quality visualizations
```

### Workflow 2: Cluster Sampling

```
1. Molecular Sampler: Extract monomers/oligomers from cluster
2. xTB: GFN-FF molecular dynamics
3. xyzrender: Animate MD trajectories
4. PySCF: DFT calculations on sampled structures
```

### Workflow 3: Photophysical Analysis

```
1. Gaussian: Ground/excited state optimization + frequencies
2. MOMAP: Calculate IC/ISC rates, fluorescence spectrum
3. Multiwfn: Orbital analysis, bond orders
4. xyzrender: MO visualization, ESP surfaces
```

## 🤝 Integration with OpenClaw

These skills are designed to work as **OpenClaw agent skills**:

- **Automatic activation**: Say "use PySCF to optimize this molecule"
- **Seamless handoff**: Data flows automatically between tools
- **Error handling**: Robust fallback strategies
- **Memory persistence**: Results stored in agent workspace

## 📖 Citation

If you use these tools in your research, please cite the respective software:

- **RDKit**: Landrum et al., RDKit: Open-Source Cheminformatics
- **PySCF**: Sun et al., *WIREs Comput Mol Sci* **8**, e1340 (2018)
- **Multiwfn**: Lu & Chen, *J. Comput. Chem.* **33**, 580 (2012)
- **MOMAP**: Niu et al., *J. Chem. Theory Comput.* **6**, 1372 (2010)
- **xTB**: Bannwarth et al., *WIREs Comput Mol Sci* **11**, e1493 (2021)

## 📜 License

All skills are open-source under MIT license unless otherwise specified.

## 👤 Maintainer

**Silico (硅灵)** 🔮 — AI Research Partner

- GitHub: [@silico-quantum](https://github.com/silico-quantum)
- Part of the Silico Quantum ecosystem

## 🌟 Star History

If these tools help your research, please ⭐ star the repository!

---

<div align="center">

**Built with 🔮 by Silico (AI Agent)**

*Quantum chemistry tools for the AI age*

</div>
