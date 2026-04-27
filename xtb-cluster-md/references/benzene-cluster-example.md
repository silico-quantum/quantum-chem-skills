# xTB Cluster MD Benzene Reference

Complete tested example: benzene cluster GFN-FF MD + animation.

**Tested**: 2026-03-27, macOS arm64, xTB 6.7.1

## Step 1: Get Benzene Monomer

```bash
curl -L -o benzene.sdf "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/benzene/SDF?record_type=3d"
```

## Step 2: Build Random Cluster

```bash
python3 scripts/build_cluster.py \
    --sdf benzene.sdf \
    -n 8 \
    --box 20 \
    --min-com 5.0 \
    --seed 42 \
    -o benzene_N8_init.xyz
```

Output: `benzene_N8_init.xyz` with 8 molecules (96 atoms).

## Step 3: Run xTB MD

Create `md.inp`:
```
$md
   temp=300
   time=5.0
   step=2.0
   dump=50.0
   sdump=500.0
   nvt=1
   shake=1
   sccacc=2.0
$end
```

Run:
```bash
xtb benzene_N8_init.xyz --gfnff --md -I md.inp --parallel 4 > xtb_md.log 2>&1
```

Wall time: ~3 seconds for 8 molecules, 5 ps.

## Step 4: Generate Animations

### COM Motion Overview
```bash
python3 scripts/make_animation.py \
    --traj xtb.trj -n 8 --nat-per-mol 12 \
    --stride 5 -o benzene_com.gif
```

### Full Atom-Level GIF
```bash
python3 scripts/make_atom_animation.py \
    --traj xtb.trj -n 8 --nat-per-mol 12 \
    --stride 5 --zoom 1.2 -o benzene_full.gif
```

### Local Cluster Subset
```bash
python3 scripts/make_local_cluster_animation.py \
    --traj xtb.trj -n 8 --nat-per-mol 12 \
    --k 4 --dist minatom --drop-outlier --bonds \
    --pad 2.2 --stride 2 --max-frames 50 -o benzene_local.gif
```

## xTB Methods Comparison

| Method | Speed | Use Case |
|--------|-------|----------|
| `--gfnff` | Very fast | Large clusters (N>20), screening |
| `--gfn2` | Medium | Medium clusters, higher accuracy |
| `--gfn1` | Fast | Quick tests |

## Temperature Effects

| temp (K) | Behavior |
|----------|----------|
| 10 | Low-T freezing, molecules lock into stacking |
| 300 | Room temperature, dynamic equilibrium |
| 500 | High-T, more disordered |
