# RDKit Chemistry Skill

This skill provides a fail-closed workflow for molecular parsing, conformer
generation, force-field relaxation, descriptors, empirical charges, and
structure interchange with RDKit.

Read [`SKILL.md`](SKILL.md) before using a generated structure in a scientific
workflow. Installing this skill does not install RDKit, PySCF, xyzrender, or
other optional software.

## Scope

The supported guidance covers:

- sanitized SMILES or file parsing with explicit failure checks;
- stereochemistry inspection and policy for unassigned centers;
- reproducible ETKDGv3 conformer generation with a recorded seed;
- MMFF94 or UFF parameter-coverage and convergence checks;
- finite-coordinate validation and sanitized SDF round trips;
- documented RDKit descriptors and finite Gasteiger charges; and
- provenance for source structures, transformations, versions, and outputs.

RDKit conformers and force-field energies are preparation or screening results.
They are not quantum-chemical geometries, crystal-packing evidence, or a basis
for inferring pi stacking, noncovalent interaction energies, ESP/RDG fields, or
TADF mechanisms.

## Installation

Use an isolated environment and record the installed version:

```bash
conda install -c conda-forge rdkit
python -c "from rdkit import rdBase; print(rdBase.rdkitVersion)"
```

The upstream project also publishes a pip package:

```bash
python -m pip install rdkit
```

Optional PySCF or rendering demonstrations require their own independently
validated environments.

## Safe Interface Example

```python
import math

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

source_smiles = "C[C@H](O)c1ccccc1"
mol = Chem.MolFromSmiles(source_smiles)
if mol is None:
    raise ValueError("SMILES parsing or sanitization failed")

stereo = Chem.FindMolChiralCenters(
    mol, includeUnassigned=True, useLegacyImplementation=False
)
if any(label == "?" for _, label in stereo):
    raise ValueError("Resolve or explicitly accept unassigned stereochemistry")

canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
mol = Chem.AddHs(mol)

parameters = AllChem.ETKDGv3()
parameters.randomSeed = 0xF00D
parameters.numThreads = 1
parameters.enforceChirality = True
conformer_id = AllChem.EmbedMolecule(mol, parameters)
if conformer_id < 0:
    raise RuntimeError("ETKDG embedding failed")

if not AllChem.MMFFHasAllMoleculeParams(mol):
    raise RuntimeError("MMFF94 does not cover every atom")
status = AllChem.MMFFOptimizeMolecule(
    mol, confId=conformer_id, maxIters=1000, mmffVariant="MMFF94"
)
if status != 0:
    raise RuntimeError(f"MMFF94 optimization did not converge: status={status}")

coordinates = mol.GetConformer(conformer_id).GetPositions()
if not all(math.isfinite(value) for row in coordinates for value in row):
    raise RuntimeError("Conformer contains a non-finite coordinate")

print(rdBase.rdkitVersion, canonical_smiles, conformer_id)
```

This snippet demonstrates interface checks only. The complete skill also
requires nonbonded-contact review, coordinate-derived stereochemistry checks,
fresh output paths, SDF round-trip validation, checksums, and an explicit
report.

## Bundled Examples

[`examples/`](examples/) contains historical demonstration molecules and saved
figures from prior environments. Their Python files are quarantined and exit
before importing scientific packages or writing outputs because they do not
satisfy this skill's current fail-closed contract. Retain them for static review
only; implement a new workflow from `SKILL.md` in a fresh run directory.

In particular, any replacement workflow that invokes PySCF or xyzrender crosses
into separate software boundaries. Its outputs must satisfy those tools' own
state, units, convergence, version, and provenance contracts.

## Acceptance Boundary

Do not accept an RDKit-generated structure unless parsing/sanitization,
stereochemistry, embedding, parameter coverage, optimization status, finite
coordinates, nonbonded-contact checks, and sanitized output round trip all pass.
Record the exact source, canonical isomeric SMILES, formal charge, radicals,
RDKit version, seed, conformer selection rule, force field, energy units,
warnings, output path, and SHA-256.

## References

- [RDKit getting started guide](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [RDKit Python API](https://www.rdkit.org/docs/api-docs.html)
- [RDKit book](https://www.rdkit.org/docs/RDKit_Book.html)

Released under the repository [MIT License](../LICENSE).
