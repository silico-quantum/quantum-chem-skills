---
name: molecular-orbital-analysis
description: Use when generating, validating, and visualizing explicitly indexed molecular orbitals from converged closed-shell or open-shell PySCF calculations.
license: MIT
compatibility: Requires Python 3.10+ and PySCF; rendering requires a cube-capable viewer, while Multiwfn is optional and version-dependent.
---

# Molecular orbital analysis

Generate molecular-orbital cube files only from a converged, fully specified
electronic-structure calculation. Keep the calculation, orbital identity, cube
grid, and visualization settings traceable. Do not infer an open-shell frontier
orbital from a closed-shell electron-count formula.

## Prerequisites

1. Use a Python environment with PySCF installed. Record the interpreter and
   package versions before calculating:

   ```bash
   python3 --version
   python3 -c 'import pyscf; print(pyscf.__version__)'
   ```

2. Locate this skill directory and use
   `scripts/generate_orbital_cubes.py`; the repository does not install a
   global `molecular-orbital-analysis` executable.
3. For images, install and version a viewer that reads Gaussian cube files,
   such as VMD, ChimeraX, or PyMOL. Verify cube support in that exact build.
4. Treat Multiwfn as optional. If it is used, record its version and follow the
   matching upstream manual because interactive menu paths can change.
5. Work in a new output directory. The bundled generator refuses to reuse one.
   It reads the source XYZ bytes once, atomically publishes them as
   `input.xyz` inside that directory, and parses only this immutable run
   snapshot.

## Input contract

Supply all of the following before running:

- a single-frame XYZ file with an exact atom count, one comment line, finite
  coordinates, and an explicit coordinate unit (`Angstrom` or `Bohr`);
- integer molecular **charge and spin**, where PySCF `spin` means
  `Nalpha - Nbeta = 2S`, not spin multiplicity;
- method, basis, and, for RKS or UKS, the exchange-correlation functional;
- an SCF cycle limit and a new output directory;
- every requested **orbital index** as a one-based `CHANNEL:INDEX` token;
- cube grid points and margin; and
- the intended positive and negative **isosurface** values for rendering.

Restricted calculations use the `spatial` channel. Unrestricted calculations
use separate `alpha` and `beta` channels. ROHF has spatial orbitals with
open-shell occupations. Never label an orbital only as “HOMO” or “LUMO” when
the spin channel and indexing convention are ambiguous.

The XYZ format has no charge, spin, method, or provenance fields. Obtain those
from an authoritative source or request them from the user; do not guess them
from the geometry.

RKS and UKS have no implicit functional in the supported runner. Supply
`--xc` explicitly; omission fails before the output directory is created.

The accepted manifest binds `input.sha256` to the run-local `input.xyz`, not to
a later reread of the source path. A concurrent source edit therefore cannot
change the geometry/hash pair used by the calculation.

## Workflow

### 1. Validate the electronic state

Check that the charge and spin are chemically intended and compatible with the
electron count:

```text
electron_count >= spin
(electron_count - spin) is even
```

Use RHF/RKS only for `spin=0`. Use ROHF or UHF for an open-shell HF state and
UKS for an open-shell Kohn-Sham state. Method selection is a scientific choice,
not an automatic fallback.

### 2. Run PySCF and inspect the orbital table

First perform a calculation without requesting cubes. Replace `<skill-dir>`
with the resolved path to this skill:

```bash
python3 <skill-dir>/scripts/generate_orbital_cubes.py molecule.xyz \
  --output-dir inspect_run \
  --unit Angstrom --charge 0 --spin 0 \
  --method rks --xc pbe0 --basis def2-svp
```

The script writes `orbitals.csv` with channel, one-based index, energy in
Hartree, and occupation. It writes cubes only after `mf.converged` is true.
Select orbitals from this table using the scientific question and occupation;
do not use `electron_count // 2` for unrestricted frontier selection.
An accepted calculation also requires the explicitly configured `scf.chk` and
`pyscf.log` to be fresh regular files, non-empty, and SHA-256-bound in
`run.json`.

### 3. Generate explicitly indexed cubes

Run a new calculation directory and repeat each desired orbital:

```bash
python3 <skill-dir>/scripts/generate_orbital_cubes.py molecule.xyz \
  --output-dir cube_run \
  --unit Angstrom --charge 0 --spin 0 \
  --method rks --xc pbe0 --basis def2-svp \
  --orbital spatial:5 --orbital spatial:6 \
  --grid-points 100 --margin-bohr 4.0
```

For UHF or UKS, use explicit spin channels, for example `alpha:8` and
`beta:7`. Each output name contains the channel and one-based index. The
manifest records the input hash, PySCF version, convergence, energy,
occupation, cube grid, exact finite payload count/range, and output hashes.
All requests are checked before cube generation. Cubes are generated under
`.partial` names and receive their final names only after every requested cube
passes complete validation.

### 4. Render without changing scientific identity

Load the original geometry and its matching cube in a versioned viewer. Render
paired surfaces at equal magnitudes, for example `+0.05` and `-0.05` in the
cube field's value units. The value is only an example starting point, not a
universal standard. Record:

- cube filename and SHA-256;
- channel and one-based orbital index;
- positive and negative isosurface values;
- colors, opacity, camera, crop, and viewer version; and
- image dimensions and background.

Use identical grid and isosurface settings for comparative figures. A visually
larger lobe after changing the isovalue is not evidence of a larger orbital
coefficient.

### 5. Optional Multiwfn analysis

If analysis beyond direct PySCF cube generation is needed, export a Molden file
from the same accepted PySCF state with the official PySCF Molden API, then use
that file with Multiwfn. Before automation, save an
interactive transcript for the exact Multiwfn version, confirm the selected
orbital/channel in its output, and test the batch input on a small fixture.
Do not reuse an undocumented numeric menu sequence from another release.

## Validation and acceptance

Accept a calculation only when all checks pass:

1. The strict XYZ parser reads only run-local `input.xyz`; its declared atom
   count and values are valid, and its hash matches `input.sha256`.
2. Charge, spin, method, basis, functional, units, and PySCF version are present
   in `run.json`.
3. `calculation.converged` is `true`; a finite printed energy alone is not SCF
   convergence.
4. Every cube request matches the calculation type: `spatial` for restricted,
   and `alpha` or `beta` for unrestricted orbitals.
5. Each one-based orbital index exists and its energy and occupation agree
   with `orbitals.csv`.
6. Every cube has complete finite origin, axis, atom-record, and scalar-payload
   fields; its atom count matches the input, its grid dimensions are positive,
   its payload contains exactly `nx * ny * nz` values, and its SHA-256 matches
   the manifest.
7. The viewer loads the molecule and both signed isosurfaces without parser
   errors; the rendered structure aligns with the cube atom records.
8. Comparative images use the declared common grid and absolute isosurface
   magnitude, unless a difference is explicitly justified.
9. `scf.chk` and `pyscf.log` are fresh for this run, non-empty, regular files,
   and their byte sizes, modification times, and SHA-256 values appear under
   `runtime_artifacts`.

For publication, also preserve the input, `pyscf.log`, `scf.chk`,
`orbitals.csv`, cubes, `run.json`, render settings, and final images.

## Failure handling

| Failure | Required response |
|---|---|
| XYZ atom count, symbol, coordinate, or unit is invalid | Stop before PySCF; repair the source input and retain its provenance. |
| Charge and spin fail the electron-parity check | Stop and obtain the intended electronic state; do not change spin automatically. |
| Restricted method requested with nonzero spin | Stop and select an appropriate open-shell method explicitly. |
| SCF does not converge | Mark the run failed and generate no accepted cube. Diagnose the state, initial guess, method, basis, and convergence settings in a new run. |
| Run-local input snapshot publication or hash check fails | Reject the run and preserve the partial/snapshot inventory; never fall back to rereading a mutable source path. |
| `scf.chk` or `pyscf.log` is missing, empty, stale, or not a regular file | Reject the run even when `mf.converged` is true; retain the artifact inventory and diagnose the PySCF output configuration. |
| Requested channel or orbital index does not exist | Use `orbitals.csv` to select a valid explicit channel and index; do not renumber silently. |
| Cube header, atom records, payload, grid, or hash fails | Reject every unpublished `.partial` cube, preserve the failed-run artifact inventory, and rerun in a new output directory after resolving the cause. The bundled CLI does not reload a checkpoint. |
| Positive or negative isosurface is missing | Check sign, level, and viewer parsing; do not present a one-sign image as the complete orbital. |
| Multiwfn batch path differs from the installed version | Stop automation and re-establish the menu path interactively against the matching manual. |

## Output and reporting

Report, at minimum:

- original source path, run-local input snapshot path/hash/byte size, atom
  count, coordinate units, charge, and PySCF spin;
- PySCF/Python versions, method, functional, basis, convergence, cycle limit,
  and total energy in Hartree;
- orbital channel, one-based orbital index, occupation, and energy in Hartree;
- cube path/hash, atom count, grid dimensions, grid margin, exact payload value
  count/range, and byte size;
- fresh `scf.chk` and `pyscf.log` paths, byte sizes, modification times, and
  SHA-256 values;
- viewer/version, signed isosurface values, colors, opacity, camera, and output
  image dimensions; and
- acceptance status, failure phase, rejected/partial artifact inventory, failed
  checks, scientific assumptions, and any optional Multiwfn
  version/transcript.

Separate computed quantities from interpretation. Orbital shape and energy do
not by themselves establish a reaction mechanism, excitation assignment, or
charge-transfer rate.

## References

- [PySCF quickstart](https://pyscf.org/quickstart.html) — restricted and
  open-shell SCF examples.
- [PySCF tools API](https://pyscf.org/pyscf_api_docs/pyscf.tools.html) —
  `pyscf.tools.cubegen.orbital` and cube-grid parameters.
- [Multiwfn official site](http://sobereva.com/multiwfn/) — obtain the manual
  matching the installed release.
- [Gaussian Cube specification](https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html)
  — file structure and interoperability considerations.
