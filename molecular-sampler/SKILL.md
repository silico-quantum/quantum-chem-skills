---
name: molecular-sampler
description: Use when separating nonperiodic molecular aggregates into distance-inferred fragments and preparing deterministic nearest-neighbor monomer through pentamer XYZ samples.
license: MIT
compatibility: Requires Python 3.10 or newer; accepts one-frame XYZ and Cartesian Gaussian GJF or COM files in angstrom for a limited element set.
metadata:
  version: "2.0.0"
---

# Molecular Sampler

Use the bundled standard-library CLI to split a nonperiodic aggregate into
distance-inferred connected components and write deterministic local
neighborhood samples. Treat every inferred fragment as a hypothesis until the
component inventory has been checked against the source structure.

## Prerequisites

- Python 3.10 or newer; no third-party Python package is required.
- A readable `.xyz`, `.gjf`, or `.com` file containing Cartesian coordinates.
- An independently justified expected fragment count; the CLI requires it and
  still requires review of the detailed inferred partition.

Inspect the exact local interface before running it:

```bash
python3 molecular_sampler.py --help
```

Run from this skill directory or address `molecular_sampler.py` by an absolute
path. Use a new output directory for every attempt.

## Input contract

### XYZ

- Exactly one XYZ frame is accepted.
- The first-line atom count must equal the number of non-empty coordinate
  records after the comment line.
- XYZ has no reliable unit or connectivity field. Confirm that coordinates are
  in angstrom and pass `--xyz-units angstrom`; this is a declaration, not an
  automatic unit check.
- Every XYZ atom is assigned to the logical `L` layer. Use `--layer all` or the
  default `--layer L`.

### Gaussian GJF or COM

- A recognizable charge and multiplicity line must precede a Cartesian
  geometry block.
- Exactly one Gaussian job section is accepted. Inputs containing `--Link1--`
  or multiple recognizable geometry sections are rejected rather than
  silently selecting one.
- Standard Cartesian and ONIOM-style atom records with trailing `H`, `M`, or
  `L` layer markers are supported. Z-matrices are not supported.
- Gaussian input coordinates are treated as the default angstrom convention.
  Files declaring `Units=Bohr` or `Units=AU` are rejected; convert them first.
- The parser extracts geometry only. Charge, multiplicity, ONIOM link atoms,
  periodic boundaries, and fragment charge/spin assignments are not propagated
  into output XYZ files.

The bundled covalent-radius table covers H, B, C, N, O, F, Si, P, S, Cl, Br,
I, and Se. An unsupported element is an error. Connectivity is inferred using
the distance cutoff
`distance < bond_scale * (covalent_radius_i + covalent_radius_j)`; the default
`bond_scale` is `1.3`. This rule is unsuitable as authoritative connectivity
for salts, coordination complexes, reactive structures, periodic systems, or
compressed snapshots.

## Workflow

### 1. Preserve and characterize the source

Record the source filename, SHA-256, coordinate provenance, expected molecule
count, selected ONIOM layer, and any known ions or coordination bonds. Do not
edit the only copy of the source to make distance inference pass.

### 2. Run a bounded sample

For a six-molecule XYZ aggregate:

```bash
python3 molecular_sampler.py aggregate.xyz \
  --xyz-units angstrom \
  --layer all \
  --expected-fragments 6 \
  --bond-scale 1.3 \
  --samples 20 \
  --output-dir aggregate_samples_run01
```

For an ONIOM low layer:

```bash
python3 molecular_sampler.py model.com \
  --layer L \
  --expected-fragments 12 \
  --output-dir model_L_samples_run01
```

`--min-fragment-atoms` defaults to `1`, so small molecules and ions are
retained. Setting a larger value is an explicit exclusion policy; every
excluded component remains listed under `dropped_fragments` in the manifest.

### 3. Understand the deterministic selection

The tool:

1. filters atoms by the selected layer;
2. infers a bond graph using the declared distance cutoff;
3. orders connected components by their smallest source atom index;
4. computes each component's arithmetic coordinate centroid;
5. for each starting component, selects its nearest `N - 1` centroid
   neighbors for target sizes two through five;
6. removes duplicate member-ID sets and stops at `--samples` unique sets for
   each size.

This is deterministic nearest-neighbor sampling, not random sampling, exhaustive
combination enumeration, RMSD diversity selection, or thermodynamic sampling.
If fewer than `N` fragments exist, the `N`-mer directory is intentionally empty;
the tool never writes an undersized structure under an `N`-mer name.

### 4. Review before downstream computation

Inspect `sampling_manifest.json` first, then representative monomers and every
complex that will be submitted. Assign charge and spin separately for each
downstream quantum-chemistry job; XYZ output does not encode them.

## Validation and acceptance

Treat the generated files as internally validated only when all of the following
hold:

1. The process exits zero and the new output directory contains both
   `sampling_manifest.json` and `sampling_summary.txt`.
2. The manifest input SHA-256 matches the source, and `parsed_atom_count` and
   `selected_atom_count` match the declared layer scope.
3. The fragment count equals `--expected-fragments`; per-fragment atom counts,
   formulas, and source indices match an independent chemical review.
4. Retained and explicitly dropped source indices form a disjoint cover of all
   selected atoms. No atom disappears silently.
5. Every generated XYZ reparses, has a finite coordinate for every atom, and
   its header atom count equals the sum for its member fragments.
6. Every dimer through pentamer has exactly the molecule count indicated by its
   directory, and member-ID tuples are unique within that size.
7. Every generated XYZ is non-empty and its SHA-256 matches its manifest
   record. Hash the final manifest and summary externally when archiving the
   run because a file cannot contain its own stable checksum.

The CLI performs strict input parsing, source-index coverage checks, output XYZ
atom-count checks, and output hashing. It creates `RUN_INCOMPLETE` before
writing artifacts and removes that marker only after the manifest and summary
are completely published. A completed manifest uses
`status: internal_validation_passed` and
`scientific_status: pending_independent_fragment_review`. The chemical
partition and inferred connectivity still require independent review; only that
separate review may advance the scientific status to accepted.

## Failure handling

| Failure | Required response |
|---|---|
| Atom count, coordinate, unit, or element validation fails | Stop and repair or convert the source; do not weaken parsing. |
| Detected fragment count differs from expectation | Preserve the report, inspect close contacts and source connectivity, and do not tune the cutoff solely to force the desired count. |
| Salt, metal, periodic, or reactive connectivity is ambiguous | Use a tool that accepts authoritative molecule IDs or explicit bonds; this sampler has no explicit-connectivity input. |
| A requested oligomer directory is empty | Confirm that enough fragments exist; never reinterpret a smaller complex as the requested size. |
| Output directory already exists | Choose a new run directory; stale files are never accepted as current output. |
| A write fails after directory creation | Treat any directory retaining `RUN_INCOMPLETE` as failed even if other files exist. Preserve it for diagnosis or move it aside before a new run. |

## Output and reporting

The output layout is:

```text
run_directory/
├── monomers/
├── dimers/
├── trimers/
├── tetramers/
├── pentamers/
├── sampling_manifest.json
└── sampling_summary.txt
```

Report the Python version, exact argv, sampler SHA-256, input name and SHA-256
of the exact captured bytes that were parsed, declared
units, parsed and selected atom counts, layer, bond-scale distance cutoff,
minimum fragment size, expected and detected fragment counts, per-fragment
formula/source indices, selection policy, member IDs, output atom counts and
SHA-256 values, dropped fragments, warnings, internal validation status, and
scientific review status. Report
charge, multiplicity, and periodicity as `not_propagated` rather than guessing.

## References

- [Benzene dimer example](references/benzene-example.md)
- Covalent radii embedded in [`molecular_sampler.py`](molecular_sampler.py)
