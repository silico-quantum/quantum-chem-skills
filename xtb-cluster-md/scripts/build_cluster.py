#!/usr/bin/env python3
"""Build a random cluster of N anthracene molecules into an XYZ file.

No RDKit/packmol dependency: parses 3D coordinates from PubChem SDF.
Overlap avoidance is coarse: enforces COM-COM minimum distance.
"""

import argparse
import math
import random
from pathlib import Path

import numpy as np


def read_sdf_coords(path: Path):
    lines = path.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError("SDF too short")
    counts = lines[3]
    nat = int(counts[0:3])
    atoms = []
    coords = []
    for i in range(4, 4 + nat):
        x = float(lines[i][0:10])
        y = float(lines[i][10:20])
        z = float(lines[i][20:30])
        el = lines[i][31:34].strip()
        atoms.append(el)
        coords.append([x, y, z])
    return atoms, np.array(coords, dtype=float)


def random_rotation_matrix(rng=random):
    # Uniform random rotation (Shoemake)
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q1 = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    q2 = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    x, y, z, w = q1, q2, q3, q4
    R = np.array([
        [1 - 2*(y*y + z*z),     2*(x*y - z*w),     2*(x*z + y*w)],
        [    2*(x*y + z*w), 1 - 2*(x*x + z*z),     2*(y*z - x*w)],
        [    2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ], dtype=float)
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdf", default="anthracene.sdf")
    ap.add_argument("-n", type=int, default=12)
    ap.add_argument("--box", type=float, default=35.0, help="cubic box side (Å) for initial placement")
    ap.add_argument("--min-com", type=float, default=7.0, help="minimum COM-COM distance (Å)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out", default="anthracene_N12_init.xyz")
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    atoms, base = read_sdf_coords(Path(args.sdf))
    base = base - base.mean(axis=0)  # center

    placements = []
    coms = []

    max_tries = 20000
    for k in range(args.n):
        placed = False
        for _ in range(max_tries):
            R = random_rotation_matrix()
            t = (np.random.rand(3) - 0.5) * args.box
            coords = base @ R.T + t
            com = coords.mean(axis=0)
            if all(np.linalg.norm(com - c) >= args.min_com for c in coms):
                placements.append(coords)
                coms.append(com)
                placed = True
                break
        if not placed:
            raise RuntimeError(f"Failed to place molecule {k+1}/{args.n}; try larger --box or smaller --min-com")

    all_atoms = []
    all_coords = []
    for m, coords in enumerate(placements, start=1):
        for el, (x, y, z) in zip(atoms, coords):
            all_atoms.append(el)
            all_coords.append((x, y, z))

    out = Path(args.out)
    with out.open("w") as f:
        f.write(f"{len(all_atoms)}\n")
        f.write(f"anthracene cluster N={args.n} seed={args.seed} box={args.box} minCOM={args.min_com}\n")
        for el, (x, y, z) in zip(all_atoms, all_coords):
            f.write(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

    print(f"Wrote {out} with {args.n} molecules, {len(all_atoms)} atoms")


if __name__ == "__main__":
    main()
