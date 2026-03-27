#!/usr/bin/env python3
"""Atom-level animated GIF for xtb.trj (multi-frame XYZ).

This renders each molecule as a colored point cloud (atoms), using a fixed
camera and fixed bounds across frames so the stacking process is visible.
No bond perception is needed.

Notes:
- Assumes constant ordering: N molecules, nat_per_mol atoms each.
- Colors by molecule index.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


def iter_xyz_frames(path: Path):
    with path.open() as f:
        while True:
            line = f.readline()
            if not line:
                return
            line = line.strip()
            if not line:
                continue
            nat = int(line)
            comment = f.readline().rstrip("\n")
            atoms = []
            coords = []
            for _ in range(nat):
                parts = f.readline().split()
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            yield comment, np.array(atoms, dtype=object), np.array(coords, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="xtb.trj")
    ap.add_argument("-n", type=int, default=12)
    ap.add_argument("--nat-per-mol", type=int, default=24)
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--max-frames", type=int, default=100)
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--size", type=int, default=800)
    ap.add_argument("--view", default="18,35", help="elev,azim")
    ap.add_argument("--zoom", type=float, default=1.0, help=">1 zooms in (reduces plot bounds)")
    ap.add_argument("-o", "--out", default="anthracene_atoms.gif")
    args = ap.parse_args()

    elev, azim = [float(x) for x in args.view.split(",")]

    traj = Path(args.traj)
    series = []
    atom_types = None
    comments = []

    for i, (comment, atoms, xyz) in enumerate(iter_xyz_frames(traj)):
        if i % args.stride != 0:
            continue
        if xyz.shape[0] != args.n * args.nat_per_mol:
            raise ValueError(f"Unexpected nat={xyz.shape[0]} expected {args.n*args.nat_per_mol}")
        if atom_types is None:
            atom_types = atoms.reshape(args.n, args.nat_per_mol)
        series.append(xyz.reshape(args.n, args.nat_per_mol, 3))
        comments.append(comment)
        if len(series) >= args.max_frames:
            break

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
    tmp_dir = Path("_atom_frames")
    tmp_dir.mkdir(exist_ok=True)

    for t in range(series.shape[0]):
        xyz_m = series[t]

        fig = plt.figure(figsize=(args.size / args.dpi, args.size / args.dpi), dpi=args.dpi)
        ax = fig.add_subplot(111, projection="3d")

        for m in range(args.n):
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
        ax.set_title(f"Anthracene MD (xTB GFN-FF, 300K)  frame {t+1}/{series.shape[0]}")
        ax.view_init(elev=elev, azim=azim)
        ax.grid(False)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

        out_png = tmp_dir / f"atom_{t:04d}.png"
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)

        frames.append(Image.open(out_png).convert("P", palette=Image.Palette.ADAPTIVE))

    out = Path(args.out)
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=70, loop=0, optimize=False)
    print(f"Wrote {out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
