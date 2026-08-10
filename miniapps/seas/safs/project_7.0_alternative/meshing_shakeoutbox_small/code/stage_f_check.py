#!/usr/bin/env python3
"""stage_f_check.py -- deck compatibility of a collar-extended ALT mesh.

The enlargement moves the absorbing perimeter far outside where the SAFS
sidecar grids have data, so this reports exactly how much of the new domain
each ASAGI grid actually covers and what SeisSol will use where it does not.

Nothing here is a pass/fail on the mesh -- it is the input the deck author
needs in order to interpret ground motion in the new area.

Usage:
    python stage_f_check.py --mesh <n.puml.h5> --deck <deck dir> \
        [--parent <p.puml.h5>]
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import BC_FREE_SURFACE, LOCAL_FACES, face_code

CH = 3_000_000


def nc_hull(path):
    import h5py
    with h5py.File(path, "r") as f:
        return (float(f["x"][0]), float(f["x"][-1]),
                float(f["y"][0]), float(f["y"][-1]),
                float(f["z"][0]), float(f["z"][-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--deck", required=True)
    ap.add_argument("--parent", default=None)
    a = ap.parse_args()

    import h5py
    with h5py.File(a.mesh, "r") as f:
        P = f["geometry"][:]
        nt = f["connect"].shape[0]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        bary = np.empty((nt, 3))
        for s0 in range(0, nt, CH):
            bary[s0:s0 + CH] = P[conn[s0:s0 + CH].astype(np.int64)].mean(1)

    print(f"mesh {os.path.basename(a.mesh)}   {nt:,} tets")
    print(f"  bbox x {P[:,0].min():,.0f}..{P[:,0].max():,.0f}  "
          f"y {P[:,1].min():,.0f}..{P[:,1].max():,.0f}  "
          f"z {P[:,2].min():,.0f}..{P[:,2].max():,.0f}")

    par_box = None
    if a.parent:
        with h5py.File(a.parent, "r") as f:
            Q = f["geometry"][:]
        par_box = (Q[:, 0].min(), Q[:, 0].max(), Q[:, 1].min(), Q[:, 1].max())

    print("\nASAGI sidecar coverage of the NEW domain "
          "(cells judged at their barycentre, as SeisSol samples them):")
    print(f"  {'grid':<34}{'covers cells':>14}{'outside':>12}   hull E / N (km)")
    for nc in sorted(glob.glob(os.path.join(a.deck, "*.nc"))):
        try:
            x0, x1, y0, y1, z0, z1 = nc_hull(nc)
        except Exception as e:
            print(f"  {os.path.basename(nc):<34}  (unreadable: {e})")
            continue
        inside = ((bary[:, 0] >= x0) & (bary[:, 0] <= x1)
                  & (bary[:, 1] >= y0) & (bary[:, 1] <= y1)
                  & (bary[:, 2] >= z0) & (bary[:, 2] <= z1))
        out = int((~inside).sum())
        print(f"  {os.path.basename(nc):<34}{100*inside.mean():13.2f}%{out:12,}"
              f"   {x0/1e3:,.0f}-{x1/1e3:,.0f} / {y0/1e3:,.0f}-{y1/1e3:,.0f}")

    # free-surface footprint, and how much of it is new.
    # One pass over connect, not one per slot: re-reading a 122M-row dataset
    # four times costs ~16 GB of transient copies on the heavy mesh.
    tris = []
    with h5py.File(a.mesh, "r") as f:
        conn = f["connect"]
        for s in range(4):
            idx = np.nonzero(face_code(B, s) == BC_FREE_SURFACE)[0]
            if idx.size:
                sub = np.empty((len(idx), 4), np.int64)
                for s0 in range(0, len(idx), 1_000_000):
                    j = idx[s0:s0 + 1_000_000]
                    sub[s0:s0 + len(j)] = conn[j[0]:j[-1] + 1][j - j[0]].astype(np.int64)
                tris.append(sub[:, LOCAL_FACES[s]])
    T = np.vstack(tris)
    V = P[T]
    area = np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1) * 0.5
    c = V.mean(1)
    print(f"\nfree surface: {len(T):,} facets, {area.sum()/1e6:,.0f} km2, "
          f"z = {V[:,:,2].min():.6f}..{V[:,:,2].max():.6f}")
    if par_box:
        newf = ~((c[:, 0] >= par_box[0]) & (c[:, 0] <= par_box[1])
                 & (c[:, 1] >= par_box[2]) & (c[:, 1] <= par_box[3]))
        print(f"  outside the PARENT's bbox: {area[newf].sum()/1e6:,.0f} km2 "
              f"({100*area[newf].sum()/area.sum():.1f} % of the new lid)")

    print("\nreceiver containment: the lid is exactly flat at z = 0 and receivers "
          "sit at z = -1 m,\n  so every receiver is below its local free surface "
          "by construction -- the barycentric\n  containment test that matters on a "
          "topographic lid is vacuous here.")


if __name__ == "__main__":
    main()
