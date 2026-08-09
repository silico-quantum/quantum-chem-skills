---
name: xyzrender
description: Use when rendering molecular structures, trajectories, transition states, or volumetric cube data as reproducible SVG, PNG, PDF, or GIF artifacts with xyzrender.
license: MIT
compatibility: Requires Python 3.10 or newer and xyzrender; SMILES, crystal, SHELXL, and interactive-orientation features require optional extras.
---

# xyzrender

Render molecular structures and scalar fields while preserving the distinction
between input evidence and graph or bond-order inference. A polished image is a
visual artifact, not proof that inferred connectivity is chemically correct.

## Prerequisites

Install the required feature set in an isolated environment:

```bash
python -m pip install xyzrender
# Optional format support:
python -m pip install "xyzrender[smi,cif,shelxl]"

xyzrender --version
xyzrender --help
```

Record the xyzrender version before relying on an option. The optional `v`
orientation interface is Linux-specific; automatic noninteractive rendering
does not require it.

## Input contract

Record the input path, SHA-256, input format, atom count, coordinate units, and
whether it contains one structure, multiple SDF records, a trajectory, a
vibrational mode, or a scalar field.

Connectivity rules are format-dependent:

- XYZ provides elements and coordinates but no bonds, bond orders, formal
  charges, or aromaticity. xyzrender/xyzgraph inference must be reviewed and
  reported as inferred.
- Mol/SDF and MOL2 can provide explicit connectivity and bond orders. Preserve
  those records by default; `--rebuild` intentionally discards them and
  re-detects the graph.
- PDB connectivity may be incomplete or convention-dependent.
- QM output and transition-state rendering require the relevant geometry and
  frequency data, not just a file extension.
- Cube rendering requires a valid volumetric grid. State whether it contains an
  orbital, density, electrostatic potential, or interaction field and retain
  its sign and units.

Do not infer coordinate units from an extension alone. Confirm them from the
generating program or source record. If an XYZ has no trustworthy provenance,
record `units=assumed_angstrom` and keep inferred connectivity unverified. Use
`--bohr` only when the source coordinates are actually in bohr. For charged or
open-shell graph reconstruction, provide `--charge` and `--multiplicity` when
required by the installed version.

## Workflow

### 1. Preflight the input and output

1. Reject missing, empty, truncated, or mixed-format input.
2. For XYZ, verify that the first-line atom count equals the number of coordinate
   records in every frame and that all coordinates are finite.
3. For a multi-record SDF, require the user to identify the intended record and
   document why it was selected before using `--mol-frame`. For a multi-frame
   XYZ static figure, extract and hash one declared frame; use the original only
   for a trajectory or ensemble render.
4. Choose a new output path. Do not accept an older artifact left by a previous
   run.
5. Save `xyzrender --version`, the exact command, and any custom JSON config.

### 2. Render a static structure

```bash
set -e
test ! -e renders_run_001
mkdir renders_run_001

# XYZ: use explicit filenames that disclose inferred connectivity.
xyzrender molecule.xyz --config paton --bo -S 1600 \
  -o renders_run_001/molecule_xyz_inferred.svg
xyzrender molecule.xyz --config paton --bo -S 2400 \
  -o renders_run_001/molecule_xyz_inferred.png

# SDF record 0: preserve file connectivity unless --rebuild is requested.
xyzrender molecule.sdf --mol-frame 0 --config paton --bo -S 1600 \
  -o renders_run_001/molecule_sdf_explicit.svg
xyzrender molecule.sdf --mol-frame 0 --config paton --bo -S 2400 \
  -o renders_run_001/molecule_sdf_explicit.png
```

Use a new run directory and fail if it already exists. This four-command
matrix prevents format/input combinations from overwriting one another.

Use `--bond` and `--unbond` only as documented visualization overrides. Record
them; they do not alter or validate the source structure.

For inferred XYZ connectivity, make the status visible in the filename and in
an adjacent manifest or figure caption. Add a visible in-figure note only when
the publication design requires it; never imply that inference is explicit
source evidence.

### 3. Render animations

```bash
set -e
test ! -e animation_run_001
mkdir animation_run_001

# Rotation GIF; GIF output uses -go/--gif-output.
xyzrender molecule.sdf --gif-rot y -go animation_run_001/molecule_rotation.gif

# Multi-frame trajectory.
xyzrender trajectory.xyz --gif-trj -go animation_run_001/trajectory.gif
```

Use `--trj-bonds` only when connectivity is expected to change between frames.
For a TS animation, require a parsed imaginary mode before using `--gif-ts`; a
geometry alone cannot define the reaction coordinate.

### 4. Render volumetric data

```bash
set -e
test ! -e cube_run_001
mkdir cube_run_001
xyzrender orbital.cube --mo --iso 0.05 -o cube_run_001/orbital.svg
xyzrender density.cube --dens --iso 0.001 -o cube_run_001/density.png
xyzrender density.cube --esp potential.cube --iso 0.001 \
  -o cube_run_001/esp.png
```

Confirm that density and ESP grids share atoms, origin, axes, dimensions, and
units before combining them. Positive and negative MO lobes depend on an
arbitrary global orbital phase; preserve the color convention across a figure
set.

### 5. Use the Python API when orchestration is needed

```python
from pathlib import Path

from xyzrender import load, render

output = Path("molecule.svg")
if output.exists():
    raise FileExistsError(f"Refusing to overwrite: {output}")
mol = load("molecule.sdf")
render(mol, output=str(output), config="paton", bo=True)
```

Use the installed API documentation as authoritative because keyword support
can change between versions.

## Validation and acceptance

Accept an artifact only when:

1. xyzrender exits successfully and the requested output was newly created and
   is non-empty.
2. The input atom count, selected record/frame count, and confirmed or assumed
   units match the run report.
3. PNG/PDF/GIF/SVG signatures match their extensions; SVG parses as XML and a
   structure-only SVG contains vector primitives rather than a raster wrapper.
4. Requested pixel size, transparency, and GIF frame count are verified from
   the artifact rather than inferred from the command. `--canvas-size` controls
   the nominal canvas setting but does not guarantee a square final raster.
5. A visual review at fit-to-page and 100% scale finds no missing atoms,
   clipping, implausible bonds, mislabeled indices, unreadable text, or hidden
   positive/negative surface lobes.
6. For XYZ, inferred bonds and bond orders are compared with an authoritative
   structure or an explicitly reviewed graph. Without that comparison, record
   `connectivity_status=unverified` even if the image is visually clean.
7. For SDF, extract an atom-pair/bond-type table from the selected source record
   and review the displayed edges and bond-order styling against it. Absence of
   `--rebuild` is necessary but not proof of fidelity. If xyzrender exposes no
   machine-readable rendered graph, report `connectivity_status=visually-reviewed`
   rather than `machine-verified`, including aromatic and stereobond limitations.

Run concrete artifact checks for every output:

```bash
test -s renders_run_001/molecule_xyz_inferred.svg
test -s renders_run_001/molecule_xyz_inferred.png
test -s renders_run_001/molecule_sdf_explicit.svg
test -s renders_run_001/molecule_sdf_explicit.png
file renders_run_001/*
xmllint --noout renders_run_001/*.svg
rg -q '<(path|circle|line|polygon|ellipse|use)\b' renders_run_001/*.svg
shasum -a 256 molecule.xyz molecule.sdf renders_run_001/*
```

For structure-only SVG, also reject unexpected embedded raster `<image>` nodes.
Use Pillow, ImageMagick, or another recorded parser to check actual PNG width,
height, alpha channel, and color mode. Generate one manifest record per output;
do not validate four artifacts with one aggregate boolean.

For publication output, obtain the target physical width, aspect ratio,
background, palette, label size, and journal requirements. For raster figures,
derive required pixels from physical inches times target DPI and verify the
actual file dimensions. Prefer SVG when a scalable vector artifact is accepted.

Two byte-identical outputs do not demonstrate that two rendering options work.
Regenerate and inspect each artifact before using it as comparative evidence.

## Failure handling

| Failure | Required response |
|---|---|
| Input parser rejects or truncates a structure | Stop and report the format, record/frame, and parser message. |
| XYZ inference gives implausible connectivity | Require reviewed explicit connectivity, or document precise `--bond`/`--unbond` overrides. |
| SDF bond orders disappear | Confirm that `--rebuild` was not used and inspect the source bond block. |
| Output is missing, empty, stale, or unparsable | Reject it even if the process returned zero; retain stdout and stderr. |
| PNG/PDF backend is unavailable | Keep a valid SVG or install the documented backend; do not relabel a raster or another format. |
| Cube grids are incompatible | Do not overlay them; regenerate compatible fields from the same geometry/grid. |
| TS or NCI auto-detection is unsupported by the input | Report not computed and use manual annotations only when independently justified. |

## Output and reporting

Report:

- xyzrender version, Python version, exact command or API call, preset, and
  custom-config checksum;
- input path, input format, input checksum, atom count, units, selected record
  or frame, and charge/multiplicity when supplied;
- connectivity source (`explicit`, `inferred`, or `manual override`) and every
  `--rebuild`, `--bond`, or `--unbond` decision;
- surface type, cube-grid provenance, isovalue, opacity, and color convention;
- output path, format, byte size, dimensions or frame count, background,
  output checksum, and visual-review status;
- connectivity status (`unverified`, `visually-reviewed`, or
  `machine-verified`) and the evidence supporting that status;
- warnings, fallbacks, and any feature reported as not computed.

## References

- [xyzrender documentation](https://xyzrender.readthedocs.io/en/latest/)
- [CLI reference](https://xyzrender.readthedocs.io/en/latest/cli_reference.html)
- [Python API guide](https://xyzrender.readthedocs.io/en/stable/python_api.html)
- [Package metadata and optional extras](https://pypi.org/project/xyzrender/)
- [Local benzene example and known limitations](references/benzene-visualization.md)

When publishing graph, TS, or NCI results derived by xyzgraph/graphRC, follow
the citation guidance in the installed xyzrender documentation.
