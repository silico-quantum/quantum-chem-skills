---
name: xyzrender
version: 1.0.0
description: Publication-quality molecular graphics. Render molecular structures as SVG, PNG, PDF, and animated GIF from XYZ, mol/SDF, MOL2, PDB, SMILES, CIF, cube files, or quantum chemistry output.
homepage: https://github.com/aligfellow/xyzrender
metadata:
  category: visualization
  tags: [chemistry, molecular-graphics, svg, gif, xyz, quantum-chemistry, 分子可视化, molecule visualization, 渲染分子, render molecule, xyz渲染, xyz render, 分子出图, molecular graphics]
---

# xyzrender - Publication-Quality Molecular Graphics

Render molecular structures as **publication-quality SVG, PNG, PDF, and animated GIF** from XYZ files and quantum chemistry input/output.

## Installation

```bash
# Using uv (recommended)
uv tool install xyzrender

# Using pip
pip install xyzrender
```

## Quick Start

### Command Line

```bash
# Basic rendering (outputs SVG by default)
xyzrender molecule.xyz

# Specify output format
xyzrender molecule.xyz -o molecule.png

# Rotation GIF
xyzrender molecule.xyz --gif-rot -o movie.gif
```

### Python API

```python
from xyzrender import load, render, render_gif

mol = load("molecule.xyz")
render(mol, output="output.svg")
render_gif("molecule.xyz", gif_rot="y", output="movie.gif")
```

## Key Features

- **Automatic detection**: Bond orders, aromaticity, TS bonds, NCI interactions
- **Multiple formats**: SVG, PNG, PDF, GIF
- **Styling presets**: default, flat, paton (PyMOL-like)
- **Annotations**: atom indices, bond measurements, custom labels
- **Surfaces**: MO, density, ESP, NCI from cube files

## Common Use Cases

```bash
# Publication figure
xyzrender molecule.xyz --config paton -o paper.svg

# Rotating molecule
xyzrender molecule.xyz --gif-rot -o rotating.gif

# Transition state
xyzrender ts.out --ts --gif-ts -o ts_animation.gif

# With vdW surface
xyzrender molecule.xyz --vdw --gif-rot -o vdw.gif
```

## Citation

> A.S. Goodfellow* and B.N. Nguyen, *J. Chem. Theory Comput.*, 2026, DOI: 10.1021/acs.jctc.5c02073
