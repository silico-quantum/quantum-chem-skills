# Molecular Sampler

This directory contains a Python 3.10+ command-line tool for deterministic,
distance-inferred fragment sampling from one-frame XYZ or Cartesian Gaussian
GJF/COM inputs. It has no third-party Python dependency.

Start with the operational contract in [SKILL.md](SKILL.md). A safe XYZ run
declares angstrom units, asserts the expected fragment count, and writes to a
new directory:

```bash
python3 molecular_sampler.py aggregate.xyz \
  --xyz-units angstrom \
  --layer all \
  --expected-fragments 6 \
  --output-dir aggregate_samples_run01
```

The command writes monomer through pentamer directories plus
`sampling_manifest.json` and `sampling_summary.txt`. It never treats XYZ bond
inference as authoritative, never silently drops small fragments by default,
and never reuses an existing output directory.

Limitations: no periodic boundaries, explicit-bond input, Z-matrix support, or
charge/spin propagation. Review every inferred fragment before downstream
quantum-chemistry calculations.

See [the benzene dimer example](references/benzene-example.md) for a bounded
walkthrough.
