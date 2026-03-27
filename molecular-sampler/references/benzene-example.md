# Molecular Sampler Benzene Reference

Example: extracting benzene molecules from a cluster and sampling oligomers.

## Preparation

Create a benzene dimer cluster (benzene_dimer.xyz):

```bash
cat > benzene_dimer.xyz << 'EOF'
24
Benzene Dimer (T-shaped)
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
C   0.000000   5.300000   0.000000
C   1.208543   6.100000   0.000000
C   1.208543   7.500000   0.000000
C   0.000000   8.300000   0.000000
C  -1.208543   7.500000   0.000000
C  -1.208543   6.100000   0.000000
H   0.000000   4.200000   0.000000
H   2.150000   5.600000   0.000000
H   2.150000   7.900000   0.000000
H   0.000000   9.400000   0.000000
H  -2.150000   7.900000   0.000000
H  -2.150000   5.600000   0.000000
EOF
```

## Basic Usage

```bash
# Extract all monomers (12 atoms each)
python3 molecular_sampler.py benzene_dimer.xyz --output-dir ./samples

# Extract 20 dimers (nearest neighbor pairs)
python3 molecular_sampler.py benzene_dimer.xyz --samples 20 --output-dir ./samples
```

## Create a Benzene Cluster (for more interesting sampling)

Build a larger cluster with multiple benzene molecules:

```python
#!/usr/bin/env python3
"""build_benzene_cluster.py — Create a benzene cluster XYZ file"""

import numpy as np
import random

# Benzene monomer template (centered at origin)
benzene = np.array([
    [0.000000, 1.395000, 0.000000],   # C
    [1.208543, 0.697500, 0.000000],   # C
    [1.208543,-0.697500, 0.000000],   # C
    [0.000000,-1.395000, 0.000000],   # C
    [-1.208543,-0.697500, 0.000000],  # C
    [-1.208543, 0.697500, 0.000000],  # C
    [0.000000, 2.479000, 0.000000],   # H
    [2.150000, 1.239500, 0.000000],   # H
    [2.150000,-1.239500, 0.000000],   # H
    [0.000000,-2.479000, 0.000000],   # H
    [-2.150000,-1.239500, 0.000000],  # H
    [-2.150000, 1.239500, 0.000000],  # H
])

random.seed(42)
n_mol = 24
nat_per_mol = 12
min_dist = 5.0  # minimum COM-COM distance (Angstrom)

# Place molecules randomly
positions = []
atoms = []
for i in range(n_mol):
    while True:
        offset = np.array([
            random.uniform(-15, 15),
            random.uniform(-15, 15),
            random.uniform(-15, 15)
        ])
        # Check minimum distance
        ok = True
        for prev in positions:
            if np.linalg.norm(offset - prev) < min_dist:
                ok = False
                break
        if ok:
            positions.append(offset)
            for atom in benzene:
                shifted = atom + offset
                atoms.append(f"C {shifted[0]:.6f} {shifted[1]:.6f} {shifted[2]:.6f}"
                             if len(atoms) % nat_per_mol < 6 else
                             f"H {shifted[0]:.6f} {shifted[1]:.6f} {shifted[2]:.6f}")
            break

# Write XYZ
with open('benzene_cluster_N24.xyz', 'w') as f:
    f.write(f"{len(atoms)}\n")
    f.write(f"Benzene cluster, {n_mol} molecules\n")
    for line in atoms:
        f.write(line + '\n')

print(f"Created benzene_cluster_N24.xyz with {n_mol} benzene molecules")
```

## Run Sampling on Benzene Cluster

```bash
# Extract all 24 monomers
python3 molecular_sampler.py benzene_cluster_N24.xyz \
    --output-dir ./benzene_samples \
    --layer all

# Sample dimers (20 nearest-neighbor pairs)
python3 molecular_sampler.py benzene_cluster_N24.xyz \
    --samples 20 \
    --output-dir ./benzene_samples

# Check results
ls ./benzene_samples/
# monomers/  dimers/  trimers/  tetramers/  pentamers/
ls ./benzene_samples/dimers/ | wc -l
# 20
```

## Typical Output

Each sampled file is in standard XYZ format:
```
24
Dimer: mol_0 + mol_1, dist=5.23 A
C   0.000000   1.395000   0.000000
...
C   0.000000   5.300000   0.000000
...
```

## Pipeline Integration

Combine with other skills:

```bash
# 1. Sample dimers from cluster
python3 molecular_sampler.py benzene_cluster_N24.xyz --samples 20

# 2. Visualize sampled structures
for f in ./benzene_samples/dimers/*.xyz; do
    xyzrender "$f" --png "renders/$(basename $f .xyz).png" --transparent --color --bond-order
done

# 3. Run PySCF on each dimer
for f in ./benzene_samples/dimers/*.xyz; do
    python3 pyscf_dimer_calc.py "$f"
done
```
