# Examples

Tested examples using benzene (C6H6) as the reference molecule.
All examples were verified on macOS arm64, 2026-03-27.

## Environment

| Tool | Version |
|------|---------|
| Python | 3.14.2 |
| PySCF | 2.12.1 |
| Multiwfn | 3.8 |
| xyzrender | latest |
| xTB | 6.7.1 |

## pyscf/

B3LYP/cc-pVDZ calculation on benzene (D2h symmetry).

### Results

| Property | Value |
|----------|-------|
| SCF Energy | -232.2627 Ha |
| HOMO (B3g) | -6.8956 eV |
| LUMO (B1u) | -0.1571 eV |
| HOMO-LUMO gap | 6.7386 eV |

### LR-TDDFT Excited States (B3LYP/cc-pVDZ)

| State | E (eV) | λ (nm) | f |
|-------|--------|--------|---|
| S1 | 5.4951 | 225.6 | 0.0000 |
| S2 | 6.2307 | 199.0 | 0.0000 |
| S3 | 7.2775 | 170.4 | 0.5737 |
| S4 | 7.2781 | 170.4 | 0.5737 |
| S5 | 7.4489 | 166.5 | 0.0000 |

### Other Results
- **Density Fitting**: -232.2627 Ha, diff = -0.0001 eV from canonical
- **PCM (cyclohexane, ε=2.0)**: -232.2630 Ha, shift = -0.0083 eV
- **NTO** (S3, brightest): weight = 1.000 (single-pair transition)

### Files
- `benzene_output.log` — Full calculation output
- `benzene_homo.cube` — HOMO orbital cube file
- `benzene_lumo.cube` — LUMO orbital cube file

## multiwfn/

Wave function analysis from PySCF-generated molden file.

### Orbital Info
- HOMO: orbital 21, -6.8956 eV
- LUMO: orbital 22, -0.1570 eV
- Gap: 6.7385 eV

### Hirshfeld Charges
| Atom | Charge |
|------|--------|
| C1-6 | -0.040 (all equivalent) |
| H1-6 | +0.040 (all equivalent) |

### Files
- `benzene.molden` — Input wavefunction file
- `benzene_hirshfeld.log` — Hirshfeld population analysis
- `benzene_orbitals.log` — Orbital energies

## xyzrender/

Molecular visualization from XYZ file.

### Rendered Files
| File | Options |
|------|---------|
| `01_basic.png` | Default |
| `02_transparent.png` | `-t` (transparent background) |
| `03_bonds.png` | `-t --bo` (bond orders) |
| `04_hires.png` | `-t --bo -S 1000` (high-res) |
| `05_basic.svg` | SVG vector format |

## molecular-sampler/

Molecular sampling from benzene cluster (12 molecules, 144 atoms).

### Results
- 12 monomers extracted
- 5 dimers, trimers, tetramers, pentamers sampled (distance-sorted nearest neighbors)

### Files
- `benzene_cluster_N12.xyz` — Input cluster
- `monomers/` — 12 individual benzene molecules
- `dimers/` — 5 nearest-neighbor pairs
- `trimers/` — 5 trimers
- `tetramers/` — 5 tetramers
- `pentamers/` — 5 pentamers
- `sampling_summary.txt` — Sampling statistics

## xtb-cluster-md/

GFN-FF molecular dynamics of benzene cluster (8 molecules, 96 atoms, 300 K, 5 ps).

### xTB MD Settings
- Method: GFN-FF
- Temperature: 300 K
- Time: 5 ps, step: 2 fs
- NVT ensemble

### Generated Animations
| File | Description | Frames |
|------|-------------|--------|
| `benzene_com.gif` | COM motion overview | 20 |
| `benzene_full.gif` | Full atom-level | 20 |
| `benzene_local.gif` | Local cluster subset with bonds | 50 |

### Files
- `benzene_N8_init.xyz` — Initial random cluster
- `xtb_md.log` — xTB MD output log
- `md.inp` — xTB MD input parameters
