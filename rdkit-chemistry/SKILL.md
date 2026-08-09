---
name: rdkit-chemistry
description: Use when parsing chemical structures, generating reproducible conformers, applying RDKit force fields, calculating descriptors or charges, or preparing structures for downstream simulation.
license: MIT
compatibility: Requires Python and RDKit; optional rendering and quantum-chemistry steps require their own separately installed packages.
---

# RDKit chemistry

Use RDKit for cheminformatics, conformer generation, force-field relaxation,
descriptors, chemical-feature perception, and structure interchange. Treat
force-field geometries and heuristic features as preparation or screening
results, not as quantum-chemical or solid-state evidence.

## Prerequisites

Install RDKit in an isolated environment and record the exact version:

```bash
python -m pip install rdkit
python -c "from rdkit import rdBase; print(rdBase.rdkitVersion)"
```

Confirm optional packages separately before using the bundled PySCF or
xyzrender demonstrations. Do not assume that installing this skill installs
any scientific software.

## Input contract

Before changing a structure, record:

- the source string or file and, for files, a checksum;
- whether salts, solvents, isotopes, atom maps, and explicit hydrogens must be
  retained;
- formal charge, radical electrons, and stereochemical intent;
- whether unassigned stereocenters represent a racemate, an unknown, or an
  input error;
- the required output format and downstream atom-order constraints.

`Chem.MolFromSmiles()` performs sanitization by default and returns `None` on
failure. Never continue after a parse or sanitization failure. SDF preserves
explicit connectivity and properties; XYZ contains elements and coordinates
but no bond orders or formal charges. RDKit conformer coordinates are in
angstrom.

Charge and spin for a later quantum calculation remain explicit user inputs.
Formal charge and radical perception can inform that decision but do not by
themselves define a defensible electronic state.

## Workflow

### 1. Parse, sanitize, and inspect stereochemistry

```python
from rdkit import Chem, rdBase

smiles = "C[C@H](O)c1ccccc1"
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    raise ValueError("SMILES parsing or sanitization failed")

canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
stereo = Chem.FindMolChiralCenters(
    mol, includeUnassigned=True, useLegacyImplementation=False
)
if any(label == "?" for _, label in stereo):
    raise ValueError("Resolve or explicitly accept unassigned stereochemistry")
```

Do not repair a chemically invalid input by silently disabling sanitization.
Partial sanitization is an expert diagnostic path and must be reported as such.

### 2. Generate a reproducible conformer ensemble

```python
from rdkit.Chem import AllChem

mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 0xF00D
params.numThreads = 1
params.enforceChirality = True
params.pruneRmsThresh = 0.5

conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
if not conf_ids:
    raise RuntimeError("ETKDG did not generate a conformer")
```

The random seed makes the request reproducible within the recorded RDKit
environment; do not promise bitwise-identical coordinates across RDKit
versions or platforms. For `EmbedMolecule`, inspect the embedding return code:
`-1` is failure and a non-negative value is a conformer ID.

### 3. Check force-field parameters and optimize

```python
if not AllChem.MMFFHasAllMoleculeParams(mol):
    raise RuntimeError("MMFF94 does not cover every atom in this molecule")

results = AllChem.MMFFOptimizeMoleculeConfs(
    mol,
    numThreads=1,
    maxIters=1000,
    mmffVariant="MMFF94",
)

converged = [
    (conf_id, energy)
    for conf_id, (status, energy) in zip(conf_ids, results)
    if status == 0
]
if not converged:
    raise RuntimeError("No MMFF94 conformer converged")

best_conf_id, best_energy = min(converged, key=lambda item: item[1])
```

MMFF energies are reported in kcal/mol. Do not compare MMFF and UFF energies.
Use UFF only when explicitly selected and only after
`AllChem.UFFHasAllMoleculeParams(mol)` succeeds. A nonzero optimization status
must not be relabeled as convergence.

### 4. Select one conformer and write a self-describing SDF

```python
import math
import os
from pathlib import Path

output = Path("selected_conformer.sdf")
temporary = output.with_suffix(".sdf.partial")
if output.exists() or temporary.exists():
    raise FileExistsError("Refusing to overwrite an existing output")

selected = Chem.Mol(mol)
selected_conf = Chem.Conformer(mol.GetConformer(best_conf_id))
selected.RemoveAllConformers()
selected.AddConformer(selected_conf, assignId=True)

coordinates = selected.GetConformer(0).GetPositions()
if not all(math.isfinite(value) for row in coordinates for value in row):
    raise RuntimeError("Selected conformer contains a non-finite coordinate")

selected.SetProp("CanonicalSMILES", canonical_smiles)
selected.SetProp("RDKitVersion", rdBase.rdkitVersion)
selected.SetIntProp("SourceConformerId", int(best_conf_id))
selected.SetDoubleProp("MMFF94Energy_kcal_mol", float(best_energy))

writer = Chem.SDWriter(str(temporary))
writer.write(selected, confId=0)
writer.close()

reloaded = Chem.SDMolSupplier(
    str(temporary), removeHs=False, sanitize=True
)[0]
if reloaded is None or reloaded.GetNumAtoms() != selected.GetNumAtoms():
    raise RuntimeError("SDF round-trip validation failed")
if not reloaded.GetConformer(0).Is3D():
    raise RuntimeError("Reloaded SDF is not marked as 3D")
if temporary.stat().st_size == 0:
    raise RuntimeError("SDF writer produced an empty file")

os.link(temporary, output)  # Fails if another process created the target.
temporary.unlink()
```

Retain explicit hydrogens for downstream 3D work unless the receiving program
has a documented reason to regenerate them. Preserve atom order or provide an
explicit atom-index map when crossing tool boundaries.

### 5. Calculate only defined properties

- Use RDKit descriptors for their documented definitions and record the RDKit
  version.
- Treat Gasteiger charges as empirical partial charges. Reject missing, `NaN`,
  or infinite values instead of replacing them with zero.
- Use RDKit chemical-feature definitions for donor/acceptor perception rather
  than broad element-only SMARTS.
- Do not infer intermolecular pi stacking, crystal packing, ESP, RDG, or a TADF
  mechanism from an isolated RDKit conformer.

## Validation and acceptance

Accept a generated structure only when all applicable checks pass:

1. Parsing and sanitization succeeded, and canonical isomeric SMILES matches
   the intended connectivity and stereochemistry.
2. Explicit-H atom count, formal charge, radical count, and atom mapping are
   unchanged from the documented transformation plan.
3. At least one conformer was embedded; the selected conformer has finite 3D
   coordinates, no implausibly short nonbonded distances, and no unintended
   stereochemical inversion. To test coordinate-derived chirality, assign
   chiral tags from the selected 3D structure on a copy and compare its assigned
   centers with the sanitized input; `enforceChirality=True` alone is not the
   acceptance check.
4. The chosen force field has all parameters, optimization status is zero, and
   the selected energy is finite.
5. Reloading the SDF with sanitization enabled succeeds, gives the expected atom
   count and connectivity, and contains a 3D conformer.
6. The SDF and report are non-empty and identify the same selected conformer.
   Treat `SourceConformerId` as provenance because SDF round-trip conformer IDs
   can be renumbered. Record an output SHA-256 after the validated rename.

For an ensemble, report the number requested, number returned after pruning,
number converged, energy range, pruning threshold, and selection rule.

## Failure handling

| Failure | Required response |
|---|---|
| Parse or sanitization failure | Stop and return the offending input plus the RDKit error; do not guess a repair. |
| Unassigned stereochemistry | Ask whether it is intentional; enumerate only under an explicit stereoisomer policy. |
| Embedding failure | Record the embedding return code and parameters. Retry first with `useRandomCoords=True` and the same seed; only then try a second declared seed. Record each attempt separately. |
| Missing force-field parameters | Stop, or switch to a user-approved method and label it separately. |
| Optimization did not converge | Increase iterations or regenerate the ensemble; retain the failed status. |
| Non-finite charge, coordinate, or energy | Reject the affected conformer or property and report it as not computed. |
| Output path exists | Refuse to overwrite unless the user explicitly requested replacement. |

## Output and reporting

Report at minimum:

- source value or path and input checksum;
- canonical isomeric SMILES, formal charge, radical electrons, and stereo
  policy;
- RDKit version, ETKDG variant, random seed, thread count, requested conformer
  count, and pruning threshold;
- force field, parameter-coverage result, maximum iterations, optimization
  status, selected conformer ID, energy, and energy units;
- atom count with and without explicit hydrogens, coordinate units, output path,
  and output checksum;
- every warning, fallback, or rejected conformer.

Separate facts computed by RDKit from downstream assumptions such as the
quantum-chemistry charge/spin state or a claim about solid-state packing.
For a literal SMILES input, define the input checksum as SHA-256 over its exact
UTF-8 bytes and report the canonical isomeric SMILES separately.

## References

- [RDKit getting started guide](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [RDKit Python API](https://www.rdkit.org/docs/api-docs.html)
- [RDKit book](https://www.rdkit.org/docs/RDKit_Book.html)
- [Bundled examples](examples/) are quarantined historical demonstrations for
  static review only. Direct execution stops before scientific imports or
  output writes; they are not supported interfaces or acceptance fixtures.
