---
name: xtb-cluster-md
description: Run xTB (GFN-FF / GFN2-xTB) molecular dynamics for organic molecular clusters (e.g., anthracene/benzene), starting from random packed geometries; then generate animated GIFs (full trajectory, COM view, local clustered subset, bond-emphasized). Use when asked to simulate many-molecule stacking/aggregation dynamics, low-temperature “freezing/condensation” behavior, or to visualize xTB MD trajectories as animations.
---

# xTB cluster MD + visualization

## Quick workflow

### 1) Get a 3D monomer geometry (PubChem SDF)
Use PubChem PUG REST and save as `molecule.sdf`:

```bash
curl -L -o benzene.sdf "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/benzene/SDF?record_type=3d"
curl -L -o anthracene.sdf "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/anthracene/SDF?record_type=3d"
```

### 2) Build a random cluster XYZ (no RDKit/packmol needed)

```bash
python3 scripts/build_cluster.py --sdf benzene.sdf -n 96 --box 90 --min-com 6.0 --seed 1 -o benzene_N96_init.xyz
```

Tips:
- If placement fails: increase `--box` or decrease `--min-com`.
- For aromatics: `--min-com ~ 6–8 Å` is usually reasonable.

### 3) Run xTB MD (GFN-FF recommended for large N)
Create an xcontrol file, e.g. `md.inp`:

```text
$md
   temp=279
   time=50.0
   step=1.0
   dump=50.0
   sdump=500.0
   nvt=1
   shake=1
   sccacc=2.0
   forcewrrestart=true
$end
```

Run:

```bash
xtb benzene_N96_init.xyz --gfnff --md -I md.inp --parallel 10 > md.log 2>&1
```

Notes:
- `dump=50 fs` ⇒ each XYZ frame corresponds to **0.05 ps**.
- `time=200 ps` ⇒ about **4000 frames**.

### 4) Make animations

#### (A) Full-system atom-level GIF (downsampled)
```bash
python3 scripts/make_atom_animation.py --traj xtb.trj -n 96 --nat-per-mol 12 --stride 20 --zoom 1.2 -o full.gif
```

#### (B) COM-only “motion overview” GIF
```bash
python3 scripts/make_animation.py --traj xtb.trj -n 96 --nat-per-mol 12 --stride 20 -o com.gif
```

#### (C) Local cluster (picked from last frame) with clearer bonds
Best for “final aggregation / hugging together”:

```bash
python3 scripts/make_local_cluster_animation.py \
  --traj xtb.trj \
  -n 96 --nat-per-mol 12 \
  --k 8 --dist minatom --drop-outlier --bonds \
  --pad 2.2 --stride 2 --max-frames 300 \
  -o local.gif
```

Local-cluster selection options:
- `--dist minatom` (recommended): molecule–molecule distance = **minimum atom–atom distance** in the last frame.
- `--drop-outlier`: after selecting the k-cluster, drop the farthest molecule from the cluster centroid (helps remove stragglers).
- `--drop-leftmost`: drop leftmost by COM-x (useful for camera framing).

## Included scripts
- `scripts/build_cluster.py`: build N-molecule random cluster XYZ from PubChem SDF.
- `scripts/make_atom_animation.py`: full atom-level GIF (colored by molecule index).
- `scripts/make_animation.py`: COM-only 3D motion GIF.
- `scripts/make_local_cluster_animation.py`: pick a “final cluster” from last frame and render only that subset; optional bond drawing.
