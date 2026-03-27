#!/usr/bin/env python3
"""Render a local 'cluster' (subset of molecules) from an xTB XYZ trajectory.

Strategy:
- Read a window of frames from a traj file (recommend a tail-truncated file).
- Determine molecule COMs in the LAST frame.
- Pick a seed molecule with the smallest average distance to its k-1 nearest neighbors.
- Select that molecule + its (k-1) nearest neighbors as the cluster.
- Render only atoms from those molecules across frames, with tight bounds.

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
                atoms.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            yield comment, np.array(atoms, dtype=object), np.array(coords, float)


def pairwise_min_atom_dist(last_xyz: np.ndarray):
    """Compute molecule-molecule distance matrix using minimum atom-atom distance.

    last_xyz: (n, nat, 3)
    returns: (n,n) symmetric matrix with inf on diagonal.
    """
    n, nat, _ = last_xyz.shape
    dmat = np.full((n, n), np.inf, float)
    for i in range(n):
        xi = last_xyz[i]  # (nat,3)
        for j in range(i + 1, n):
            xj = last_xyz[j]
            # brute force min over atom pairs
            # (nat,nat,3) -> (nat,nat)
            diff = xi[:, None, :] - xj[None, :, :]
            dij = np.sqrt((diff * diff).sum(axis=-1)).min()
            dmat[i, j] = dmat[j, i] = float(dij)
    return dmat


def pick_cluster_by_dmat(dmat: np.ndarray, k: int):
    """Pick a fixed-size cluster using a precomputed molecule distance matrix."""
    n = dmat.shape[0]
    scores = []
    neighs = []
    for i in range(n):
        idx = np.argsort(dmat[i])
        idx = idx[idx != i][: max(0, k - 1)]
        neighs.append(idx)
        scores.append(dmat[i, idx].mean() if len(idx) else 1e9)
    seed = int(np.argmin(scores))
    cluster = [seed] + [int(x) for x in neighs[seed]]
    return cluster


def infer_bonds(atoms_m, xyz_m):
    # simple distance-based bonds, good enough for anthracene
    # cutoffs in Å
    cut = {
        ('C','C'): 1.75,
        ('C','H'): 1.25,
        ('H','C'): 1.25,
    }
    bonds=[]
    n=len(atoms_m)
    for i in range(n):
        for j in range(i+1,n):
            a,b=atoms_m[i], atoms_m[j]
            if (a,b) not in cut:
                continue
            d=np.linalg.norm(xyz_m[i]-xyz_m[j])
            if d <= cut[(a,b)]:
                bonds.append((i,j))
    return bonds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", default="xtb_tail_2000plus.trj")
    ap.add_argument("-n", type=int, default=24)
    ap.add_argument("--nat-per-mol", type=int, default=24)
    ap.add_argument("--k", type=int, default=6, help="number of molecules in local cluster")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--max-frames", type=int, default=160)
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--size", type=int, default=1000)
    ap.add_argument("--view", default="20,35", help="elev,azim")
    ap.add_argument("--pad", type=float, default=3.5, help="extra Å padding around cluster bounds")
    ap.add_argument("--drop-leftmost", action="store_true", help="drop the leftmost molecule in the selected cluster (by COM x in last frame)")
    ap.add_argument("--drop-outlier", action="store_true", help="drop the most distant molecule from cluster centroid (by COM distance in last frame)")
    ap.add_argument("--bonds", action="store_true", help="draw intra-molecular bonds (distance-based) to make structure clearer")
    ap.add_argument("--dist", choices=["com", "minatom"], default="minatom", help="molecule-molecule distance metric used to pick the cluster from the last frame")
    ap.add_argument("-o", "--out", default="anthracene_local_cluster.gif")
    args = ap.parse_args()

    elev, azim = [float(x) for x in args.view.split(",")]

    series = []
    atom_types = None

    # load frames (already should be tail window)
    for i, (comment, atoms, xyz) in enumerate(iter_xyz_frames(Path(args.traj))):
        if i % args.stride != 0:
            continue
        if xyz.shape[0] != args.n * args.nat_per_mol:
            raise ValueError(f"Unexpected nat={xyz.shape[0]} expected {args.n*args.nat_per_mol}")
        if atom_types is None:
            atom_types = atoms.reshape(args.n, args.nat_per_mol)
        series.append(xyz.reshape(args.n, args.nat_per_mol, 3))
        if len(series) >= args.max_frames:
            break

    series = np.array(series)  # (T,n,nat,3)
    com_last = series[-1].mean(axis=1)  # (n,3)

    if args.dist == "com":
        dmat = np.linalg.norm(com_last[:, None, :] - com_last[None, :, :], axis=-1)
        np.fill_diagonal(dmat, np.inf)
    else:
        dmat = pairwise_min_atom_dist(series[-1])

    cluster = pick_cluster_by_dmat(dmat, args.k)
    cluster = sorted(cluster)

    if args.drop_leftmost and len(cluster) > 1:
        xs = com_last[cluster, 0]
        drop = int(cluster[int(np.argmin(xs))])
        cluster = [m for m in cluster if m != drop]
        print(f"dropped leftmost molecule (0-based): {drop}")

    if args.drop_outlier and len(cluster) > 1:
        cen = com_last[cluster].mean(axis=0)
        ds = np.linalg.norm(com_last[cluster] - cen[None, :], axis=1)
        drop = int(cluster[int(np.argmax(ds))])
        cluster = [m for m in cluster if m != drop]
        print(f"dropped outlier molecule (0-based): {drop}")

    print(f"distance metric: {args.dist}")
    print("cluster mol indices (0-based):", cluster)

    # infer bonds once from first frame of first molecule in cluster
    bond_list = None
    if args.bonds:
        m0 = cluster[0]
        bond_list = infer_bonds(atom_types[m0], series[0, m0])
        print(f"inferred {len(bond_list)} bonds per molecule")

    # compute bounds from all frames for selected molecules
    pts = series[:, cluster, :, :].reshape(series.shape[0], -1, 3)
    mins = pts.min(axis=(0, 1))
    maxs = pts.max(axis=(0, 1))
    mins = mins - args.pad
    maxs = maxs + args.pad

    cmap = plt.get_cmap("tab10")
    frames = []
    tmp_dir = Path("_local_frames")
    tmp_dir.mkdir(exist_ok=True)

    for t in range(series.shape[0]):
        fig = plt.figure(figsize=(args.size / args.dpi, args.size / args.dpi), dpi=args.dpi)
        ax = fig.add_subplot(111, projection="3d")

        for j, m in enumerate(cluster):
            atoms_m = atom_types[m]
            pts_m = series[t, m]
            is_h = (atoms_m == 'H')
            col = cmap(j % 10)

            # bonds (within each molecule)
            if bond_list is not None:
                for (a, b) in bond_list:
                    ax.plot([pts_m[a, 0], pts_m[b, 0]],
                            [pts_m[a, 1], pts_m[b, 1]],
                            [pts_m[a, 2], pts_m[b, 2]],
                            color=col, alpha=0.75, linewidth=1.6)

            ax.scatter(pts_m[~is_h, 0], pts_m[~is_h, 1], pts_m[~is_h, 2],
                       s=30, alpha=0.98, color=col, linewidths=0)
            ax.scatter(pts_m[is_h, 0], pts_m[is_h, 1], pts_m[is_h, 2],
                       s=6, alpha=0.18, color=col, linewidths=0)

        ax.set_xlim(mins[0], maxs[0]); ax.set_ylim(mins[1], maxs[1]); ax.set_zlim(mins[2], maxs[2])
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"Local cluster (k={args.k}) — frame {t+1}/{series.shape[0]}")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.grid(False)

        out_png = tmp_dir / f"local_{t:04d}.png"
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)

        frames.append(Image.open(out_png).convert("P", palette=Image.Palette.ADAPTIVE))

    out = Path(args.out)
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=60, loop=0, optimize=False)
    print(f"Wrote {out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
