#!/usr/bin/env python3
"""Render a local 'cluster' (subset of molecules) from an xTB XYZ trajectory.

Strategy:
- Read a window of frames from a traj file (recommend a tail-truncated file).
- Determine arithmetic coordinate centroids in the LAST frame.
- Pick a seed molecule with the smallest average distance to its k-1 nearest neighbors.
- Select that molecule + its (k-1) nearest neighbors as the cluster.
- Render only atoms from those molecules across frames, with tight bounds.

Output: GIF via Pillow.
"""

from __future__ import annotations

import argparse
import json
from io import BytesIO
from pathlib import Path

from trajectory_io import (
    load_accepted_sampled_trajectory,
    save_verified_gif,
    sha256_file,
)


def pairwise_min_atom_dist(last_xyz, np_module):
    """Compute molecule-molecule distance matrix using minimum atom-atom distance.

    last_xyz: (n, nat, 3)
    returns: (n,n) symmetric matrix with inf on diagonal.
    """
    n, nat, _ = last_xyz.shape
    dmat = np_module.full((n, n), np_module.inf, float)
    for i in range(n):
        xi = last_xyz[i]  # (nat,3)
        for j in range(i + 1, n):
            xj = last_xyz[j]
            # brute force min over atom pairs
            # (nat,nat,3) -> (nat,nat)
            diff = xi[:, None, :] - xj[None, :, :]
            dij = np_module.sqrt((diff * diff).sum(axis=-1)).min()
            dmat[i, j] = dmat[j, i] = float(dij)
    return dmat


def pick_cluster_by_dmat(dmat, k: int, np_module):
    """Pick a fixed-size cluster using a precomputed molecule distance matrix."""
    n = dmat.shape[0]
    scores = []
    neighs = []
    for i in range(n):
        idx = np_module.argsort(dmat[i])
        idx = idx[idx != i][: max(0, k - 1)]
        neighs.append(idx)
        scores.append(dmat[i, idx].mean() if len(idx) else 1e9)
    seed = int(np_module.argmin(scores))
    cluster = [seed] + [int(x) for x in neighs[seed]]
    return cluster


def load_explicit_bonds(
    manifest_path: Path, atoms_per_molecule: int, expected_labels: list[str]
) -> tuple[list[tuple[int, int]], dict[str, str]]:
    """Load display topology from the original accepted SDF build manifest."""
    if not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        raise ValueError("build manifest is absent or empty")
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("status") != "accepted":
        raise ValueError("build manifest does not mark the source structure accepted")
    source = data.get("source")
    if not isinstance(source, dict):
        raise ValueError("build manifest has no source metadata")
    if source.get("atoms_per_molecule") != atoms_per_molecule:
        raise ValueError("build manifest atom count does not match the animation")
    if source.get("atom_labels") != expected_labels:
        raise ValueError("build manifest atom labels/order do not match the trajectory")
    raw_bonds = source.get("bonds_zero_based")
    if not isinstance(raw_bonds, list):
        raise ValueError("build manifest has no explicit SDF bond topology")
    bonds: list[tuple[int, int]] = []
    for entry in raw_bonds:
        if (
            not isinstance(entry, list)
            or len(entry) != 3
            or not all(isinstance(value, int) for value in entry)
        ):
            raise ValueError("build manifest contains a malformed bond record")
        first, second, order = entry
        if not (
            0 <= first < atoms_per_molecule
            and 0 <= second < atoms_per_molecule
            and first != second
            and order > 0
        ):
            raise ValueError("build manifest contains an invalid bond record")
        bonds.append((first, second))
    return bonds, {
        "build_manifest_path": str(manifest_path.resolve()),
        "build_manifest_sha256": sha256_file(manifest_path),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("-n", "--molecules", dest="molecules", type=int, required=True)
    ap.add_argument("--nat-per-mol", type=int, default=24)
    ap.add_argument("--k", "--cluster-size", dest="cluster_size", type=int, required=True, help="number of molecules in local cluster")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--max-frames", type=int, default=160)
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--size", type=int, default=1000)
    ap.add_argument("--view", default="20,35", help="elev,azim")
    ap.add_argument("--pad", type=float, default=3.5, help="extra Å padding around cluster bounds")
    ap.add_argument("--drop-leftmost", action="store_true", help="drop the leftmost molecule by coordinate-centroid x in the last frame")
    ap.add_argument("--drop-outlier", action="store_true", help="drop the molecule farthest from the cluster's coordinate centroid in the last frame")
    ap.add_argument("--bonds", action="store_true", help="draw SDF bonds from --build-manifest")
    ap.add_argument("--build-manifest", type=Path, help="accepted cluster builder manifest required with --bonds")
    ap.add_argument("--dist", "--distance", dest="distance", choices=["centroid", "minatom"], default="minatom", help="molecule-distance metric used for last-frame subset selection")
    ap.add_argument("--title", default="Local molecular cluster")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()

    if args.molecules <= 0 or args.nat_per_mol <= 0 or args.cluster_size <= 0:
        raise ValueError("molecules, nat-per-mol, and cluster-size must be positive")
    if args.cluster_size > args.molecules:
        raise ValueError("cluster-size cannot exceed molecule count")
    if args.stride <= 0 or args.max_frames <= 0 or args.pad < 0:
        raise ValueError("stride and max-frames must be positive; pad cannot be negative")
    if args.bonds and args.build_manifest is None:
        raise ValueError("--bonds requires --build-manifest")
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

    source_frames, provenance = load_accepted_sampled_trajectory(
        args.traj,
        args.validation,
        args.molecules,
        args.nat_per_mol,
        args.stride,
        args.max_frames,
    )
    for _, labels, coordinates in source_frames:
        atoms = np.array(labels, dtype=object)
        xyz = np.array(coordinates, float)
        if atom_types is None:
            atom_types = atoms.reshape(args.molecules, args.nat_per_mol)
        series.append(xyz.reshape(args.molecules, args.nat_per_mol, 3))

    if not series or atom_types is None:
        raise ValueError("trajectory yielded no selected frames")
    series = np.array(series)  # (T,n,nat,3)
    centroid_last = series[-1].mean(axis=1)  # (n,3)

    if args.distance == "centroid":
        dmat = np.linalg.norm(centroid_last[:, None, :] - centroid_last[None, :, :], axis=-1)
        np.fill_diagonal(dmat, np.inf)
    else:
        dmat = pairwise_min_atom_dist(series[-1], np)

    cluster = pick_cluster_by_dmat(dmat, args.cluster_size, np)
    cluster = sorted(cluster)

    if args.drop_leftmost and len(cluster) > 1:
        xs = centroid_last[cluster, 0]
        drop = int(cluster[int(np.argmin(xs))])
        cluster = [m for m in cluster if m != drop]
        print(f"dropped leftmost molecule (0-based): {drop}")

    if args.drop_outlier and len(cluster) > 1:
        cen = centroid_last[cluster].mean(axis=0)
        ds = np.linalg.norm(centroid_last[cluster] - cen[None, :], axis=1)
        drop = int(cluster[int(np.argmax(ds))])
        cluster = [m for m in cluster if m != drop]
        print(f"dropped outlier molecule (0-based): {drop}")

    print(f"distance metric: {args.distance}")
    print("cluster mol indices (0-based):", cluster)

    bond_list = None
    build_provenance: dict[str, str] = {}
    if args.bonds:
        m0 = cluster[0]
        bond_list, build_provenance = load_explicit_bonds(
            args.build_manifest,
            args.nat_per_mol,
            [str(value) for value in atom_types[m0]],
        )
        print(f"loaded {len(bond_list)} explicit SDF bonds per molecule")

    # compute bounds from all frames for selected molecules
    pts = series[:, cluster, :, :].reshape(series.shape[0], -1, 3)
    mins = pts.min(axis=(0, 1))
    maxs = pts.max(axis=(0, 1))
    mins = mins - args.pad
    maxs = maxs + args.pad

    cmap = plt.get_cmap("tab10")
    frames = []
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
        ax.set_title(f"{args.title} (selected={len(cluster)}) — frame {t+1}/{series.shape[0]}")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.grid(False)

        fig.tight_layout()
        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        plt.close(fig)

        buffer.seek(0)
        with Image.open(buffer) as rendered:
            frames.append(rendered.convert("P", palette=Image.Palette.ADAPTIVE).copy())

    provenance.update(build_provenance)
    save_verified_gif(
        frames,
        args.out,
        args.manifest,
        60,
        provenance,
        {
            "type": "local_molecule_subset",
            "molecules": args.molecules,
            "atoms_per_molecule": args.nat_per_mol,
            "selected_molecules_zero_based": cluster,
            "distance_metric": args.distance,
            "stride": args.stride,
            "maximum_rendered_frames": args.max_frames,
            "size_pixels": args.size,
            "dpi": args.dpi,
            "view": args.view,
            "padding_angstrom": args.pad,
            "title": args.title,
            "bonds": "explicit SDF topology" if args.bonds else "not drawn",
            "limitation": "the selected local subset is a visualization, not an unbiased aggregation statistic",
        },
    )
    print(f"Wrote {args.out} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
