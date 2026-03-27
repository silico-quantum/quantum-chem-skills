# xyzrender Benzene Reference

Benzene visualization examples using xyzrender CLI.

## Preparation

```bash
# Create benzene.xyz
cat > benzene.xyz << 'XYZEOF'
12
Benzene D6h
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
xyzrender benzene.xyz --png benzene_basic.png
xyzrender benzene.xyz --png benzene_color.png --color
xyzrender benzene.xyz --png benzene_transparent.png --transparent --color
xyzrender benzene.xyz --png benzene_bonds.png --transparent --color --bond-order
```

## High Resolution

```bash
xyzrender benzene.xyz --png benzene_hires.png --transparent --color --bond-order --scale 2.0
xyzrender benzene.xyz --png benzene_ultra.png --transparent --color --bond-order --scale 3.0
```

## Different Input Formats

```bash
xyzrender benzene.sdf --png benzene_from_sdf.png --transparent --color --bond-order
xyzrender "SMILES:c1ccccc1" --png benzene_smiles.png --transparent --color
xyzrender benzene.pdb --png benzene_pdb.png --transparent --color --bond-order
```

## Output Formats

```bash
xyzrender benzene.xyz --svg benzene.svg --transparent --color --bond-order
xyzrender benzene.xyz --pdf benzene.pdf --transparent --color --bond-order
xyzrender benzene.xyz --gif benzene.gif --transparent --color --bond-order
```

## Option Reference

| Option | Effect | Recommended |
|--------|--------|-------------|
| `--png` | PNG output | Default |
| `--transparent` | Transparent background | Papers |
| `--color` | CPK element coloring | Clarity |
| `--bond-order` | Show single/double bonds | Organic molecules |
| `--scale 2.0` | 2x resolution | Publication |
| `--svg` | Vector format | Scalable figures |

## Batch Processing

```bash
mkdir -p renders
for mol in benzene toluene naphthalene anthracene; do
    xyzrender ${mol}.xyz --png renders/${mol}.png --transparent --color --bond-order --scale 2.0
done
```
