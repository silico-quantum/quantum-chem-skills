# xyzrender benzene reference

This page records a small benzene rendering workflow and the limitations of
the repository's saved artifacts. It is not a version-pinned acceptance test.

## Input

Use a validated 12-atom benzene XYZ in angstrom. An XYZ file does not contain
bond or bond-order records, so connectivity and aromaticity in these commands
are inferred by xyzrender/xyzgraph.

Before rendering, check that the XYZ has exactly 12 coordinate records, finite
coordinates, and the expected `C6H6` composition. Record its SHA-256.

## Commands

Confirm the installed interface first:

```bash
xyzrender --version
xyzrender --help
```

Render separate, newly named outputs:

```bash
set -e
test ! -e benzene_render_run_001
mkdir benzene_render_run_001
xyzrender benzene.xyz -o benzene_render_run_001/benzene_basic.png
xyzrender benzene.xyz -t --bo -o benzene_render_run_001/benzene_bonds.png
xyzrender benzene.xyz -t --bo -S 1600 \
  -o benzene_render_run_001/benzene_hires.png
xyzrender benzene.xyz --bo -o benzene_render_run_001/benzene.svg
xyzrender benzene.xyz --gif-rot y \
  -go benzene_render_run_001/benzene_rotation.gif
```

`-S`/`--canvas-size` controls the nominal canvas size; it does not guarantee a
square final raster. Measure the generated artifact instead of reporting the
requested setting as the actual dimensions.

For explicit connectivity, prepare and validate an SDF and render it without
`--rebuild`:

```bash
set -e
test ! -e benzene_sdf_run_001
mkdir benzene_sdf_run_001
xyzrender benzene.sdf --bo -o benzene_sdf_run_001/benzene_from_sdf.svg
```

## Batch pattern

```bash
set -e
test ! -e renders
mkdir renders
for mol in benzene toluene naphthalene anthracene; do
    input="${mol}.xyz"
    output="renders/${mol}.svg"
    test -s "$input" || exit 1
    test ! -e "$output" || exit 1
    xyzrender "$input" --bo -S 1600 -o "$output" || exit 1
    test -s "$output" || exit 1
done
```

## Artifact checks

```bash
file benzene_render_run_001/benzene_basic.png \
  benzene_render_run_001/benzene.svg \
  benzene_render_run_001/benzene_rotation.gif
xmllint --noout benzene_render_run_001/benzene.svg
shasum -a 256 benzene.xyz benzene_render_run_001/*
```

Also inspect the images at 100% scale and check atom count, inferred bonds,
bond-order styling, clipping, background, and labels.

## Repository snapshot limitations

The files under [`examples/xyzrender/`](../../examples/xyzrender/) are a
historical snapshot without a recorded xyzrender version, source-input hash,
or generation log. In that snapshot:

- `02_transparent.png` and `03_bonds.png` are byte-identical, so they do not
  demonstrate an effect from `--bo`;
- the raster called high resolution is 1000 by 1126 pixels, not 1000 by 1000;
- transparent images can appear to have a black background in viewers that
  composite alpha against black;
- the saved images must be regenerated before they are used as current
  validation evidence.

Treat these files as display examples only. A reproducible result needs the
input, exact command, xyzrender version, output dimensions, visual review, and
input/output checksums.
