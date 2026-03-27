# xTB Cluster MD Benzene Reference

Complete example: building a benzene cluster, running GFN-FF MD, and generating animations.

## Step 1: Get Benzene Monomer

```bash
curl -L -o benzene.sdf "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/benzene/SDF?record_type=3d"
```

## Step 2: Build Random Cluster

```bash
python3 scripts/build_cluster.py \
    --sdf benzene.sdf \
    -n 36 \
    --box 40 \
    --min-com 6.0 \
    --seed 42 \
    -o benzene_N36_init.xyz
```

Parameters explained:
- `-n 36`: 36 benzene molecules (432 atoms)
- `--box 40`: 40 Å cubic box
- `--min-com 6.0`: minimum 6 Å between molecule centers
- `--seed 42`: reproducible random placement

## Step 3: Run xTB MD

Create `md.inp`:
```
$md
   temp=300
   time=100.0
   step=1.0
   dump=50.0
   sdump=500.0
   nvt=1
   shake=1
   sccacc=2.0
   forcewrrestart=true
$end
```

Run GFN-FF MD (recommended for 36 molecules):
```bash
xtb benzene_N36_init.xyz --gfnff --md -I md.inp --parallel 10 > md.log 2>&1
```

Key parameters:
- `temp=300`: 300 K (room temperature)
- `time=100.0`: 100 ps total simulation
- `step=1.0`: 1 fs timestep
- `dump=50.0`: save frame every 50 fs (0.05 ps) → ~2000 frames
- `nvt=1`: NVT ensemble (thermostat on)
- `shake=1`: SHAKE constraint for C-H bonds (faster)
- `--gfnff`: GFN-FF force field (much faster than GFN2-xTB for large N)

For higher accuracy (slower):
```bash
xtb benzene_N36_init.xyz --gfn2 --md -I md.inp --parallel 10 > md.log 2>&1
```

## Step 4: Generate Animations

### Full Atom-Level GIF
```bash
python3 scripts/make_atom_animation.py \
    --traj xtb.trj \
    -n 36 --nat-per-mol 12 \
    --stride 10 --zoom 1.2 \
    -o benzene_full.gif
```

### COM Motion Overview
```bash
python3 scripts/make_animation.py \
    --traj xtb.trj \
    -n 36 --nat-per-mol 12 \
    --stride 10 \
    -o benzene_com.gif
```

### Local Cluster Subset (aggregation detail)
```bash
python3 scripts/make_local_cluster_animation.py \
    --traj xtb.trj \
    -n 36 --nat-per-mol 12 \
    --k 6 --dist minatom --drop-outlier --bonds \
    --pad 2.2 --stride 2 --max-frames 300 \
    -o benzene_local.gif
```

## Step 5: Analysis

### Check convergence from log:
```bash
grep -E "energy|temperature" md.log | tail -20
```

### Extract geometry snapshots:
```bash
# Frame at 50 ps (frame ~1000)
head -2 xtb.trj | head -1  # get atom count
# Use trjconv or custom script to extract specific frames
```

## Temperature Variations

| System | temp (K) | Behavior |
|--------|----------|----------|
| 10 | Low-T freezing | Molecules lock into stacking |
| 300 | Room temperature | Dynamic equilibrium |
| 500 | High-T | More disordered, less aggregation |

## Full Pipeline Script

```bash
#!/bin/bash
set -e

MOL="benzene"
N=36
BOX=40
TEMP=300
TIME=100.0

echo "=== Step 1: Get monomer ==="
curl -sL -o ${MOL}.sdf \
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${MOL}/SDF?record_type=3d"

echo "=== Step 2: Build cluster ==="
python3 scripts/build_cluster.py \
    --sdf ${MOL}.sdf -n ${N} --box ${BOX} --min-com 6.0 --seed 42 \
    -o ${MOL}_N${N}_init.xyz

echo "=== Step 3: Run xTB MD ==="
xtb ${MOL}_N${N}_init.xyz --gfnff --md -I md.inp --parallel 10 > md.log 2>&1

echo "=== Step 4: Animations ==="
python3 scripts/make_atom_animation.py \
    --traj xtb.trj -n ${N} --nat-per-mol 12 \
    --stride 10 --zoom 1.2 -o ${MOL}_full.gif

python3 scripts/make_local_cluster_animation.py \
    --traj xtb.trj -n ${N} --nat-per-mol 12 \
    --k 6 --dist minatom --drop-outlier --bonds \
    --pad 2.2 --stride 2 --max-frames 300 -o ${MOL}_local.gif

echo "=== Done! ==="
ls -lh *.gif
```

## xTB Methods Comparison

| Method | Speed | Accuracy | Use Case |
|--------|-------|----------|----------|
| GFN-FF | Very fast | Semi-quantitative | Large clusters (N>20), screening |
| GFN2-xTB | Medium | Good | Medium clusters (N<20), benchmarking |
| GFN1-xTB | Fast | Moderate | Quick tests |
