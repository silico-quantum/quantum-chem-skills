# Example Artifact Index

Preview files from the saved benzene examples. These artifacts document a prior run; this index does not independently reproduce or validate the calculations.

## xyzrender Renders

| Style | Description |
|-------|-------------|
| Basic | Default white background |
| Transparent | No background (for papers) |
| Bond orders | Aromatic bonds shown |
| Hi-res | 1000×1000 canvas |

The saved `02_transparent.png` and `03_bonds.png` files are byte-identical in
the current snapshot. The style labels describe the intended commands; they do
not establish that the two checked-in images differ.

See [`../../examples/xyzrender/`](../../examples/xyzrender/) for the complete files.

## Rejected Historical xTB Cluster MD Animations

- **Legacy `benzene_com.gif`** — Arithmetic coordinate-centroid motion, not
  center-of-mass motion despite the filename
- **Full atom** — All atoms colored by molecule index
- **Local cluster** — Subset with bond drawing

The retained logs contain a thermostat problem and runtime exception, while
the original trajectory and `xtbmdok` marker are absent. These GIFs are visual
artifacts from a rejected historical run, not accepted MD evidence.

See [`../../examples/xtb-cluster-md/`](../../examples/xtb-cluster-md/) for the GIF files.

## Molecular Sampling

12 benzene molecules → monomers, dimers, trimers, tetramers, pentamers.

See [`../../examples/molecular-sampler/`](../../examples/molecular-sampler/) for the XYZ files.
