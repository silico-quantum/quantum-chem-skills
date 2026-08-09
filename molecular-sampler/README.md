# Molecular Sampler

Molecular Sampler extracts connected molecular fragments from Gaussian ONIOM GJF or XYZ geometries. It writes every detected monomer and generates deterministic nearest-neighbor dimers, trimers, tetramers, and pentamers.

## Quick Start

From this directory, run:

```bash
python3 molecular_sampler.py <input_file> [options]
```

For example:

```bash
python3 molecular_sampler.py guest_monomer.gjf \
  --layer L \
  --samples 20 \
  --output-dir ./molecular_samples
```

## Options

- `input_file`: Gaussian GJF or standard XYZ input.
- `--output-dir`: output directory; default: `./molecular_samples`.
- `--samples`: maximum number of starting molecules sampled for each complex size; default: `20`.
- `--layer`: `H`, `L`, or `all`; default: `L`. XYZ atoms are assigned to the `L` layer.

## Method

### Fragment Detection

- Infers bonds from tabulated covalent radii with a `1.3 × (r1 + r2)` cutoff.
- Uses Union-Find connected components to identify molecules.
- Excludes connected components with fewer than five atoms.
- Uses each molecule's coordinate centroid for neighbor-distance calculations.

### Complex Sampling

- Sorts molecules deterministically by centroid coordinates.
- Sorts each molecule's neighbors by centroid distance.
- Combines a starting molecule with its nearest neighbors to form each requested complex size.
- Uses up to the first `--samples` starting molecules for dimers through pentamers.

The sampler is deterministic and distance-based. It does not perform random sampling, enforce structural diversity, or account for periodic boundary conditions.

## Output

```text
molecular_samples/
├── monomers/
│   ├── monomer_01.xyz
│   ├── monomer_02.xyz
│   └── ...
├── dimers/
├── trimers/
├── tetramers/
├── pentamers/
└── sampling_summary.txt
```

Every generated structure uses standard XYZ format. The summary records the input path, selected layer, detected-molecule count, and number of files generated for each sample type.

## Dependencies

- Python 3.7 or later
- NumPy

Install NumPy in your chosen environment, then run the script directly:

```bash
python3 -m pip install numpy
python3 molecular_sampler.py structure.xyz --layer all
```

## Limitations

- The built-in covalent-radius table covers `H`, `C`, `N`, `O`, `S`, `P`, `F`, `Cl`, `Br`, and `I`; other elements use a fallback radius.
- Fragments with fewer than five atoms are omitted.
- Bond detection evaluates every atom pair and can be slow for large systems.
- The GJF parser supports only the geometry layouts recognized by the bundled implementation.
- Charge, multiplicity, bonding, and chemical validity must be checked before downstream quantum-chemistry calculations.

## Example

See [Benzene cluster sampling](references/benzene-example.md) for a complete walkthrough.

## Version History

- **v1.0.0** (2026-03-03): initial release with monomer extraction, nearest-neighbor oligomer sampling, and standard XYZ output.
