# Publication-Quality Figures for Quantum-Chem-Skills README

All figures generated using benzene (C6H6) as the test molecule.

## Generated Figures

### 1. Benzene Structure (01_benzene_structure.png)
- **Size:** 44K
- **Source:** Copied from examples/xyzrender/03_bonds.png
- **Description:** Benzene molecular structure with bond order visualization

### 2. Molecular Orbitals (02_orbitals.png)
- **Size:** 86K
- **Method:** PySCF B3LYP/cc-pVDZ (symmetry=False)
- **Content:** 3-panel figure showing HOMO-1, HOMO, and LUMO
- **Colormap:** RdBu_r
- **Resolution:** 200 dpi, transparent background

### 3. UV-Vis Absorption Spectrum (03_uvvis.png)
- **Size:** 70K
- **Method:** LR-TDDFT with 20 excited states, Gaussian broadening (σ=0.12 eV)
- **Axes:** Primary (Energy, eV) + Secondary (Wavelength, nm)
- **Features:** Top peaks labeled with wavelength and oscillator strength
- **Resolution:** 200 dpi, transparent background

### 4. Absorption & Emission Spectra (04_abs_em.png)
- **Size:** 105K
- **Method:** Absorption (LR-TDDFT) + Emission (Stokes shift 0.3 eV, scaled 0.4×)
- **Features:** Spectral overlap region shaded, Stokes shift annotated
- **Axes:** Dual axis (eV + nm)
- **Colors:** Absorption (#4361ee), Emission (#e63946), Overlap (#2a9d8f)
- **Resolution:** 200 dpi, transparent background

### 5. Potential Energy Surfaces (05_pes.png)
- **Size:** 129K
- **Method:** B3LYP/STO-3G with TDA, C-C bond scan 1.2-1.65 Å (30 points)
- **Content:** S₀ and S₁ potential energy curves
- **Annotations:** 
  - Equilibrium bond lengths (r_eq)
  - Vertical excitation energy
  - Adiabatic excitation energy
- **Resolution:** 200 dpi, transparent background

### 6. MD Aggregation (06_md.png)
- **Size:** 113K
- **Method:** xtb GFN-FF MD simulation (5 ps, 300 K)
- **Content:** 
  - Left panel: Average nearest-neighbor COM distance vs time
  - Right panel: Distance distribution histogram (first 20% vs last 20%)
- **Resolution:** 200 dpi, transparent background

## Technical Details

### Software Versions
- Python 3.x with PySCF 2.12
- Matplotlib
- NumPy
- SciPy
- xtb 6.7.1

### Figure Specifications
- **DPI:** 200
- **Background:** Transparent
- **Font sizes:** Titles 14-15 bold, Labels 13, Annotations 10-11
- **Format:** PNG with tight bounding box

### Color Palette
- Absorption: #4361ee (blue)
- Emission: #e63946 (red)
- Overlap: #2a9d8f (teal)

---

Generated on: 2026-03-27
