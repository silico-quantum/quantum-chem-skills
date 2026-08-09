---
name: molecular-sampler
description: Use when extracting molecular fragments from Gaussian GJF or XYZ geometries and generating deterministic nearest-neighbor monomer, dimer, trimer, tetramer, and pentamer samples.
license: MIT
compatibility: Requires Python 3.7 or later, NumPy, and local access to the input geometry file.
---

# Molecular Sampler

Extract connected molecular fragments from a Gaussian GJF or XYZ geometry, write every detected monomer, and construct nearby multi-molecule complexes.

## Workflow

1. Confirm that the input is a `.gjf` or `.xyz` file.
2. For an ONIOM GJF file, select the `H`, `L`, or `all` layer scope. XYZ atoms are assigned to the `L` layer by default.
3. Choose an output directory and the maximum number of starting molecules to sample for each complex size.
4. Run the sampler.
5. Inspect `sampling_summary.txt` and spot-check the generated XYZ files before using them in downstream calculations.

## Basic Usage

Run the script from this skill directory:

```bash
python3 molecular_sampler.py <input_file> [options]
```

### Arguments

- `input_file`: path to a Gaussian GJF or standard XYZ file.
- `--output-dir`: output directory; default: `./molecular_samples`.
- `--samples`: maximum number of starting molecules used for each multi-molecule sample type; default: `20`.
- `--layer`: ONIOM layer selection: `H`, `L`, or `all`; default: `L`.

### Examples

```bash
# Sample the low-level layer from an ONIOM input.
python3 molecular_sampler.py guest_monomer.gjf --layer L --samples 20

# Process all atoms in an XYZ file and choose an output directory.
python3 molecular_sampler.py structure.xyz --layer all --output-dir ./my_samples
```

## Output Layout

```text
molecular_samples/
├── monomers/           # Every detected monomer
│   ├── monomer_01.xyz
│   ├── monomer_02.xyz
│   └── ...
├── dimers/             # Up to --samples structures
├── trimers/            # Up to --samples structures
├── tetramers/          # Up to --samples structures
├── pentamers/          # Up to --samples structures
└── sampling_summary.txt
```

Each generated structure is a standard XYZ file. Its comment line records the sample index, molecule count, and atom count.

## Implemented Sampling Method

### Molecule Identification

- Filter atoms by the selected ONIOM layer.
- Infer bonds when the interatomic distance is less than `1.3 × (r1 + r2)`, where `r1` and `r2` are tabulated covalent radii.
- Use a Union-Find connected-component algorithm to identify fragments.
- Discard connected components containing fewer than five atoms.
- Represent each retained molecule by the arithmetic mean of its atomic coordinates.

### Multi-Molecule Sampling

1. Calculate all pairwise distances between molecular coordinate centroids.
2. Sort the neighbor list of each molecule by distance.
3. For an N-molecule complex, combine each starting molecule with its nearest `N - 1` neighbors.
4. Iterate over at most the first `--samples` molecules in deterministic coordinate-sorted order.

This is deterministic nearest-neighbor sampling, not random or diversity-optimized sampling. The requested oligomer size also requires enough detected molecules; inspect the output rather than assuming every generated sample contains the target count.

## Input Format

A standard XYZ file has this form:

```text
<atom_count>
<comment>
<element> <x> <y> <z>
...
```

Gaussian inputs must contain a geometry block that the bundled parser recognizes. ONIOM inputs must include valid `H` or `L` layer markers.

## Dependencies

- Python 3.7 or later
- NumPy
- Python standard-library modules: `argparse`, `collections`, `os`, `re`, and `sys`

## Appropriate Uses

- Extracting molecular clusters from crystal or aggregate geometries
- Preparing nearby fragment complexes for quantum-chemistry calculations
- Separating sufficiently disconnected molecular fragments in ONIOM or XYZ geometries

## Validation Notes

- Bond perception depends on the built-in covalent-radius table and a fixed `1.3` tolerance factor.
- Components with fewer than five atoms are intentionally omitted, so the tool is unsuitable for retaining small molecules or ions without code changes.
- The implementation evaluates all atom pairs during bond detection; large systems may require substantial runtime.
- Review charge, multiplicity, chemical identity, and periodic-boundary effects independently before downstream calculations.

## Extended Example

See [Benzene cluster sampling](references/benzene-example.md) for a complete input-generation and sampling walkthrough.

## Version History

- **v1.0.0** (2026-03-03): initial release with monomer extraction, distance-sorted oligomer sampling, and standard XYZ output.
