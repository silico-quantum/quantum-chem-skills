#!/usr/bin/env python3
"""Make an animated GIF showing the stacking dynamics of an anthracene cluster.

Reads xtb.trj (multi-frame XYZ). Assumes N molecules, each with nat_per_mol atoms
in the same ordering across frames.

Visualization: 3D COM trajectory (one point per molecule) + faint lines from
previous frame to current to emphasize motion.

Output: GIF via Pillow.
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
                if len(parts) < 4:
                    raise ValueError("Bad XYZ line")
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            yield comment, atoms, np.array(coords, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="xtb.trj")
    ap.add_argument("-n", type=int, default=12, help="number of molecules")
    ap.add_argument("--nat-per-mol", type=int, default=24)
    ap.add_argument("--stride", type=int, default=5, help="use every k-th frame")
    ap.add_argument("--max-frames", type=int, default=240)
    ap.add_argument("--size", type=int, default=700)
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("-o", "--out", default="anthracene_N12_md.gif")
    args = ap.parse_args()

    traj = Path(args.traj)
    frames = []
    com_prev = None

    # collect a bounded set of COM frames
    com_series = []
    comments = []
    for i, (comment, atoms, xyz) in enumerate(iter_xyz_frames(traj)):
        if i % args.stride != 0:
            continue
        if xyz.shape[0] != args.n * args.nat_per_mol:
            raise ValueError(f"Unexpected nat={xyz.shape[0]} (expected {args.n*args.nat_per_mol})")
        xyz_m = xyz.reshape(args.n, args.nat_per_mol, 3)
        com = xyz_m.mean(axis=1)
        com_series.append(com)
        comments.append(comment)
        if len(com_series) >= args.max_frames:
            break

    com_series = np.array(com_series)  # (T, n, 3)

    # determine plot bounds
    mins = com_series.min(axis=(0, 1))
    maxs = com_series.max(axis=(0, 1))
    center = 0.5 * (mins + maxs)
    span = (maxs - mins).max() * 0.65
    if span <= 1e-6:
        span = 10.0

    def lim(c):
        return (c - span, c + span)

    tmp_dir = Path("_frames")
    tmp_dir.mkdir(exist_ok=True)

    for t in range(com_series.shape[0]):
        com = com_series[t]

        fig = plt.figure(figsize=(args.size / args.dpi, args.size / args.dpi), dpi=args.dpi)
        ax = fig.add_subplot(111, projection="3d")

        # plot previous positions faintly
        if t > 0:
            prev = com_series[t - 1]
            for k in range(args.n):
                ax.plot([prev[k, 0], com[k, 0]], [prev[k, 1], com[k, 1]], [prev[k, 2], com[k, 2]],
                        color="#999999", linewidth=0.8, alpha=0.6)

        ax.scatter(com[:, 0], com[:, 1], com[:, 2], s=35, c=np.linspace(0, 1, args.n), cmap="viridis")

        ax.set_xlim(*lim(center[0])); ax.set_ylim(*lim(center[1])); ax.set_zlim(*lim(center[2]))
        ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
        ax.set_title(f"Anthracene cluster MD (xTB GFN-FF)  frame {t+1}/{com_series.shape[0]}")
        ax.view_init(elev=18, azim=35)
        ax.grid(False)

        out_png = tmp_dir / f"frame_{t:04d}.png"
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)

        frames.append(Image.open(out_png).convert("P", palette=Image.Palette.ADAPTIVE))

    out = Path(args.out)
    # duration in ms per frame: adjust for smoothness
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=60, loop=0, optimize=False)
    print(f"Wrote {out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
