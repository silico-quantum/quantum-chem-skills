#!/usr/bin/env python3
"""Atom-level animated GIF for xtb.trj (multi-frame XYZ).

This renders each molecule as a colored point cloud (atoms), using a fixed
camera and fixed bounds across frames so the stacking process is visible.
No bond perception is needed.

Notes:
- Assumes constant ordering: N molecules, nat_per_mol atoms each.
- Colors by molecule index.
"""

from __future__ import annotations

import argparse
from io import BytesIO
from pathlib import Path

from trajectory_io import load_accepted_sampled_trajectory, save_verified_gif


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", type=Path, default=Path("xtb.trj"))
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("-n", "--molecules", dest="molecules", type=int, required=True)
    ap.add_argument("--nat-per-mol", type=int, default=24)
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--max-frames", type=int, default=100)
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--size", type=int, default=800)
    ap.add_argument("--view", default="18,35", help="elev,azim")
    ap.add_argument("--zoom", type=float, default=1.0, help=">1 zooms in (reduces plot bounds)")
    ap.add_argument("--title", default="Cluster molecular dynamics")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()

    if args.molecules <= 0 or args.nat_per_mol <= 0 or args.stride <= 0 or args.max_frames <= 0:
        raise ValueError("molecules, nat-per-mol, stride, and max-frames must be positive")
    if args.zoom <= 0:
        raise ValueError("zoom must be positive")
    try:
        import numpy as np
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from PIL import Image
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "animation requires numpy, matplotlib, and Pillow"
        ) from exc
    elev, azim = [float(x) for x in args.view.split(",")]

    series = []
    atom_types = None
    comments = []
    source_frames, provenance = load_accepted_sampled_trajectory(
        args.traj,
        args.validation,
        args.molecules,
        args.nat_per_mol,
        args.stride,
        args.max_frames,
    )
    for comment, labels, coordinates in source_frames:
        atoms = np.array(labels, dtype=object)
        xyz = np.array(coordinates, float)
        if atom_types is None:
            atom_types = atoms.reshape(args.molecules, args.nat_per_mol)
        series.append(xyz.reshape(args.molecules, args.nat_per_mol, 3))
        comments.append(comment)

    if not series or atom_types is None:
        raise ValueError("trajectory yielded no selected frames")
    series = np.array(series)  # (T,n,nat,3)

    mins = series.min(axis=(0, 1, 2))
    maxs = series.max(axis=(0, 1, 2))
    center = 0.5 * (mins + maxs)
    span = (maxs - mins).max() * 0.60
    if span < 10:
        span = 10
    # zoom in by shrinking bounds
    if args.zoom and args.zoom > 0:
        span = span / args.zoom

    def lim(c):
        return (c - span, c + span)

    cmap = plt.get_cmap("tab20")
    frames = []
    for t in range(series.shape[0]):
        xyz_m = series[t]

        fig = plt.figure(figsize=(args.size / args.dpi, args.size / args.dpi), dpi=args.dpi)
        ax = fig.add_subplot(111, projection="3d")

        for m in range(args.molecules):
            atoms_m = atom_types[m]
            pts = xyz_m[m]
            is_h = (atoms_m == 'H')
            # heavy atoms
            ax.scatter(pts[~is_h, 0], pts[~is_h, 1], pts[~is_h, 2],
                       s=10, alpha=0.95, color=cmap(m % 20), linewidths=0)
            # hydrogens faint
            ax.scatter(pts[is_h, 0], pts[is_h, 1], pts[is_h, 2],
                       s=2, alpha=0.25, color=cmap(m % 20), linewidths=0)

        ax.set_xlim(*lim(center[0])); ax.set_ylim(*lim(center[1])); ax.set_zlim(*lim(center[2]))
        ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
        ax.set_title(f"{args.title} — frame {t+1}/{series.shape[0]}")
        ax.view_init(elev=elev, azim=azim)
        ax.grid(False)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

        fig.tight_layout()
        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        plt.close(fig)

        buffer.seek(0)
        with Image.open(buffer) as rendered:
            frames.append(rendered.convert("P", palette=Image.Palette.ADAPTIVE).copy())

    save_verified_gif(
        frames,
        args.out,
        args.manifest,
        70,
        provenance,
        {
            "type": "atom_point_cloud",
            "molecules": args.molecules,
            "atoms_per_molecule": args.nat_per_mol,
            "stride": args.stride,
            "maximum_rendered_frames": args.max_frames,
            "size_pixels": args.size,
            "dpi": args.dpi,
            "view": args.view,
            "zoom": args.zoom,
            "title": args.title,
            "limitation": "points are colored by molecule block; no bonds are drawn",
        },
    )
    print(f"Wrote {args.out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
