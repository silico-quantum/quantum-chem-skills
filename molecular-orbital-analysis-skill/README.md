# Molecular Orbital Analysis Skill

This skill contains tools and scripts for analyzing molecular orbitals and electronic structure outputs. Typical features:

- Parse quantum chemistry output (e.g., from PySCF/ Gaussian) to extract orbital energies and coefficients
- Visualize orbitals (isosurfaces/slices) using common viewers or convert to standard formats
- Compute Mulliken / Löwdin populations and generate summary tables

Usage

1. Place calculation output files into `inputs/` directory.
2. Run: `python analyze_orbitals.py inputs/ output_summary.md`

Requirements

- Python 3.11+
- numpy, scipy, matplotlib

See scripts in this directory for examples.
