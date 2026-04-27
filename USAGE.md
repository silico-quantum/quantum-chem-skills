# Usage

Each skill directory is self-contained. Install into your OpenClaw skills directory:

```bash
cp -r pyscf multiwfn momap molecular-sampler xyzrender xtb-cluster-md ~/.openclaw/skills/
```

## Quick Start

### PySCF — DFT calculation
```bash
python3 pyscf/references/benzene-dft-tddft.py
```

### xyzrender — Molecular visualization
```bash
xyzrender molecule.xyz -o output.png -t --bo
```

### Molecular Sampler — Extract oligomers
```bash
python3 molecular-sampler/molecular_sampler.py cluster.xyz --samples 20 --output-dir ./samples
```

### xTB Cluster MD
```bash
python3 xtb-cluster-md/scripts/build_cluster.py --sdf molecule.sdf -n 24 -o cluster.xyz
xtb cluster.xyz --gfnff --md -I md.inp
python3 xtb-cluster-md/scripts/make_animation.py --traj xtb.trj -n 24 --nat-per-mol 12 -o anim.gif
```

### Multiwfn — Wave function analysis
```
Multiwfn molecule.molden
# Then select functions from the menu (7=population, 9=bond order, etc.)
```

## Verified Examples

See [`examples/`](examples/README.md) for tested benzene examples with actual output.
