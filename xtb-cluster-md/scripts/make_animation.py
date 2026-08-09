#!/usr/bin/env python3
"""Make an animated GIF of molecular coordinate centroids in an XYZ trajectory.

Reads xtb.trj (multi-frame XYZ). Assumes N molecules, each with nat_per_mol atoms
in the same ordering across frames.

Visualization: one arithmetic coordinate centroid per molecule plus faint lines
from the previous frame. These points are not mass-weighted centers of mass.

Output: GIF via Pillow.
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
    ap.add_argument("--stride", type=int, default=5, help="use every k-th frame")
    ap.add_argument("--max-frames", type=int, default=240)
    ap.add_argument("--size", type=int, default=700)
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--title", default="Molecular coordinate centroids")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()

    if args.molecules <= 0 or args.nat_per_mol <= 0 or args.stride <= 0 or args.max_frames <= 0:
        raise ValueError("molecules, nat-per-mol, stride, and max-frames must be positive")
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
    frames = []
    centroid_series = []
    comments = []
    source_frames, provenance = load_accepted_sampled_trajectory(
        args.traj,
        args.validation,
        args.molecules,
        args.nat_per_mol,
        args.stride,
        args.max_frames,
    )
    for comment, _, coordinates in source_frames:
        xyz = np.array(coordinates, float)
        xyz_m = xyz.reshape(args.molecules, args.nat_per_mol, 3)
        centroids = xyz_m.mean(axis=1)
        centroid_series.append(centroids)
        comments.append(comment)

    if not centroid_series:
        raise ValueError("trajectory yielded no selected frames")
    centroid_series = np.array(centroid_series)  # (T, n, 3)

    # determine plot bounds
    mins = centroid_series.min(axis=(0, 1))
    maxs = centroid_series.max(axis=(0, 1))
    center = 0.5 * (mins + maxs)
    span = (maxs - mins).max() * 0.65
    if span <= 1e-6:
        span = 10.0

    def lim(c):
        return (c - span, c + span)

    for t in range(centroid_series.shape[0]):
        centroids = centroid_series[t]

        fig = plt.figure(figsize=(args.size / args.dpi, args.size / args.dpi), dpi=args.dpi)
        ax = fig.add_subplot(111, projection="3d")

        # plot previous positions faintly
        if t > 0:
            prev = centroid_series[t - 1]
            for k in range(args.molecules):
                ax.plot([prev[k, 0], centroids[k, 0]], [prev[k, 1], centroids[k, 1]], [prev[k, 2], centroids[k, 2]],
                        color="#999999", linewidth=0.8, alpha=0.6)

        ax.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2], s=35, c=np.linspace(0, 1, args.molecules), cmap="viridis")

        ax.set_xlim(*lim(center[0])); ax.set_ylim(*lim(center[1])); ax.set_zlim(*lim(center[2]))
        ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
        ax.set_title(f"{args.title} — frame {t+1}/{centroid_series.shape[0]}")
        ax.view_init(elev=18, azim=35)
        ax.grid(False)

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
        60,
        provenance,
        {
            "type": "arithmetic_coordinate_centroids",
            "molecules": args.molecules,
            "atoms_per_molecule": args.nat_per_mol,
            "stride": args.stride,
            "maximum_rendered_frames": args.max_frames,
            "size_pixels": args.size,
            "dpi": args.dpi,
            "title": args.title,
            "limitation": "centroids are arithmetic coordinate centroids, not centers of mass",
        },
    )
    print(f"Wrote {args.out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
