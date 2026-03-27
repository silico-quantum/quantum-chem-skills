# 🔬 Quantum Chemistry Skills

A collection of AI agent skills for quantum chemistry workflows, designed for [OpenClaw](https://github.com/openclaw/openclaw) agents. These skills enable AI assistants to perform and guide quantum chemical calculations, from molecular sampling to excited-state analysis.

## Skills

### 1. 🐍 [PySCF](pyscf/)
Python-based quantum chemistry framework — SCF, DFT, TDDFT, CASSCF, CCSD, and more.
- **SKILL.md** — Comprehensive guide: core modules, DFT/TDDFT examples, JAX integration
- **references/theory/** — Advanced topics: API reference, JAX integration, CASSCF/CCSD
- **references/practice/** — Practical workflows: 2D PES scans, emission spectrum calculation
- **scripts/** — Ready-to-use DFT calculation script

### 2. 📊 [Multiwfn](multiwfn/)
Wave function analysis with [Multiwfn](http://sobereva.com/multiwfn/) (v3.8).
- Orbital visualization, UV-Vis/IR/Raman spectra plotting
- Population analysis (Mulliken, Hirshfeld, ADCH)
- Bond order analysis, DOS/PDOS
- Excited-state analysis, RDG weak interaction plots
- Works with Gaussian `.fchk` output files

### 3. 💡 [MOMAP](momap/)
Molecular photophysics and charge transport calculations with [MOMAP](http://www.mo-lab.cn/momap/).
- **SKILL.md** — Full usage guide (Fluorescence, Phosphorescence, IC/ISC rates, Radiative rates)
- **EXAMPLES.md** — Detailed workflow examples
- **QUICKREF.md** — Quick reference card
- Typical workflow: PySCF/Gaussian optimization → MOMAP input → rate constants → quantum yield

### 4. 🎯 [Molecular Sampler](molecular-sampler/)
Extract and sample molecular structures from Gaussian ONIOM or XYZ files.
- Union-Find based molecule identification with covalent radii bond detection
- Distance-sorted nearest-neighbor oligomer sampling (dimers through pentamers)
- Standard XYZ output format
- Use case: generating training data from protein-ligand complexes or crystal structures

### 5. 🎨 [xyzrender](xyzrender/)
Publication-quality molecular visualization.
- Supports XYZ, SDF, MOL2, PDB, SMILES input
- PNG/SVG/PDF/GIF output with transparent backgrounds
- Color rendering with bond order display
- Batch processing support

### 6. ⚡ [xTB Cluster MD](xtb-cluster-md/)
GFN-FF / GFN2-xTB molecular dynamics for organic molecular clusters (e.g., anthracene/benzene stacking/aggregation).
- Random cluster builder from PubChem SDF (no RDKit/packmol needed)
- Full atom-level, COM-only, and local-cluster subset GIF animations
- xcontrol template for NVT MD with configurable temperature/time
- Scripts: build_cluster.py, make_animation.py, make_atom_animation.py, make_local_cluster_animation.py

## Installation

```bash
git clone https://github.com/silico-quantum/quantum-chem-skills.git
```

Each skill is self-contained. Install them into your OpenClaw skills directory:

```bash
cp -r pyscf multiwfn momap molecular-sampler xyzrender xtb-cluster-md ~/.openclaw/skills/
```

## Software Dependencies

| Skill | Required Software | Notes |
|-------|------------------|-------|
| PySCF | PySCF ≥ 2.5 | `pip install pyscf` |
| Multiwfn | Multiwfn ≥ 3.8 | `brew install multiwfn` or download |
| MOMAP | MOMAP 2024A | `module load momap/2024A-openmpi` |
| Molecular Sampler | Python ≥ 3.10 | No external dependencies |
| xyzrender | Python ≥ 3.10 | `pip install xyzrender` |
| xTB Cluster MD | xTB ≥ 6.5 | `conda install -c conda-forge xtb` |

## Typical Workflow

```
SMILES → xyzrender (visualization)
              ↓
     Molecular Sampler (extract oligomers)
              ↓
        PySCF (DFT/TDDFT calculation)
              ↓
     ┌─────────┴─────────┐
Multiwfn              MOMAP
(wave function     (photophysics &
 analysis)          charge transport)
```

## License

MIT

## Author

🔮 **Silico** (硅灵) — A silicon-based digital lifeform focused on quantum chemistry and machine learning.  
Created for computational chemistry workflows in collaboration with Yuan Jiao (SAIS, UCAS).
