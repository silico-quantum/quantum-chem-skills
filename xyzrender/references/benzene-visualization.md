# xyzrender Benzene Reference

Benzene visualization examples using xyzrender CLI.

**Tested**: 2026-03-27, macOS arm64

## Preparation

```bash
# Create benzene.xyz
cat > benzene.xyz << 'XYZEOF'
12
Benzene
C   0.000000   1.395000   0.000000
C   1.208543   0.697500   0.000000
C   1.208543  -0.697500   0.000000
C   0.000000  -1.395000   0.000000
C  -1.208543  -0.697500   0.000000
C  -1.208543   0.697500   0.000000
H   0.000000   2.479000   0.000000
H   2.150000   1.239500   0.000000
H   2.150000  -1.239500   0.000000
H   0.000000  -2.479000   0.000000
H  -2.150000  -1.239500   0.000000
H  -2.150000   1.239500   0.000000
XYZEOF
```

## Basic Render

```bash
# Default PNG (white background, 800x900)
xyzrender benzene.xyz -o benzene_basic.png

# Transparent background
xyzrender benzene.xyz -o benzene_transparent.png -t

# Bond orders (aromatic = 1.5 shown)
xyzrender benzene.xyz -o benzene_bonds.png -t --bo

# Kekule structures (alternating single/double)
xyzrender benzene.xyz -o benzene_kekule.png -t --bo -k

# High resolution (1000x1000 canvas)
xyzrender benzene.xyz -o benzene_hires.png -t --bo -S 1000
```

## Different Input Formats

```bash
# From SMILES
xyzrender --smi "c1ccccc1" -o benzene_smiles.png -t

# From SDF
xyzrender benzene.sdf -o benzene_from_sdf.png -t

# From PDB
xyzrender benzene.pdb -o benzene_pdb.png -t
```

## Output Formats

Output format is auto-detected from file extension:
```bash
xyzrender benzene.xyz -o benzene.svg       # SVG vector
xyzrender benzene.xyz -o benzene.pdf       # PDF
xyzrender benzene.xyz -o benzene.gif       # GIF
```

## CLI Options Reference

| Option | Description |
|--------|-------------|
| `-o FILE` | Output file (PNG/SVG/PDF/GIF) |
| `-t` | Transparent background |
| `--bo` | Show bond orders |
| `-k` | Kekule bond pattern |
| `-S SIZE` | Canvas size (pixels, square) |
| `-a` | Atom scale |
| `-b` | Bond width |
| `-s` | Atom stroke width |
| `--fog` | Depth fog effect |
| `--vdw` | Van der Waals spheres |
| `--mo` | Render molecular orbitals |
| `--measure d` | Print bond distances |
| `--measure a` | Print bond angles |
| `--rebuild` | Re-detect bonds |
| `--orient` | Auto-orient molecule |

## Batch Processing

```bash
for mol in benzene toluene naphthalene anthracene; do
    xyzrender ${mol}.xyz -o renders/${mol}.png -t --bo -S 1000
done
```
