#!/usr/bin/env python3
"""check_receivers_local_top.py -- Stage F receiver containment on a NON-UNIFORM lid.

SeisSol v1.1.3 silently DROPS receivers that sit above their local free surface, so
"inside the horizontal footprint" is not sufficient on a mesh whose top is not flat
(the PREFERRED lid carries fault-trace vertices from -49.9 to +25.8 m while
receivers are placed at z = -1 m).

For each receiver this locates the free-surface triangle whose (x,y) projection
contains it, interpolates the lid height there barycentrically, and checks the
receiver is BELOW it.  Vertex-proximity is NOT good enough -- the far-field lid is
kilometres coarse, so the nearest top VERTEX can belong to a triangle that does not
contain the point at all.

    python check_receivers_local_top.py --mesh X.puml.h5 --receivers r.dat [--margin 0.0]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import h5py
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from puml_io import LOCAL_FACES  # noqa: E402

CH = 4_000_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--receivers", required=True)
    ap.add_argument("--margin", type=float, default=0.0,
                    help="required clearance below the local lid [m]")
    a = ap.parse_args()

    rec = np.loadtxt(a.receivers)
    if rec.ndim == 1:
        rec = rec[None, :]
    print(f"[recv] {len(rec):,} receivers, z in "
          f"[{rec[:,2].min():.3f}, {rec[:,2].max():.3f}]", flush=True)

    f = h5py.File(a.mesh, "r")
    P = f["geometry"][:]
    ds = f["connect"]
    nt = ds.shape[0]
    bnd = f["boundary"][:]
    bu = np.ascontiguousarray(bnd, np.int32).view(np.uint32)
    codes = np.stack([((bu >> np.uint32(8 * s)) & np.uint32(0xFF)).astype(np.uint8)
                      for s in range(4)], 1)
    del bu, bnd
    tris = []
    for s0 in range(0, nt, CH):
        T = ds[s0:s0 + CH].astype(np.int64)
        cc = codes[s0:s0 + CH]
        for s in range(4):
            sel = np.nonzero(cc[:, s] == 1)[0]
            if sel.size:
                tris.append(P[T[sel][:, LOCAL_FACES[s]]])
    f.close()
    A = np.concatenate(tris, 0)                       # (m,3,3)
    del tris
    print(f"[lid] {len(A):,} free-surface triangles, z in "
          f"[{A[:,:,2].min():.3f}, {A[:,:,2].max():.3f}]", flush=True)

    cen = A[:, :, :2].mean(1)
    rad = np.linalg.norm(A[:, :, :2] - cen[:, None, :], axis=2).max(1)
    tree = cKDTree(cen)
    rmax = float(rad.max())
    xy = np.ascontiguousarray(rec[:, :2])
    cand = tree.query_ball_point(xy, r=rmax, workers=-1)

    unlocated, above = [], []
    zsurf = np.full(len(rec), np.nan)
    for i, cl in enumerate(cand):
        if not cl:
            unlocated.append(i)
            continue
        c = np.asarray(cl)
        Q = A[c]
        v0 = Q[:, 1, :2] - Q[:, 0, :2]
        v1 = Q[:, 2, :2] - Q[:, 0, :2]
        v2 = xy[i] - Q[:, 0, :2]
        den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
        den = np.where(np.abs(den) < 1e-12, 1e-12, den)
        u = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
        v = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
        tol = 1e-9
        hit = np.nonzero((u >= -tol) & (v >= -tol) & (u + v <= 1 + tol))[0]
        if not hit.size:
            unlocated.append(i)
            continue
        # highest containing triangle -> the most conservative lid
        zs = (Q[hit, 0, 2] + u[hit] * (Q[hit, 1, 2] - Q[hit, 0, 2])
              + v[hit] * (Q[hit, 2, 2] - Q[hit, 0, 2]))
        zsurf[i] = float(zs.max())
        if rec[i, 2] > zsurf[i] - a.margin:
            above.append(i)

    ok = np.isfinite(zsurf)
    clear = zsurf[ok] - rec[ok, 2]
    print()
    print(f"[PASS] located          {int(ok.sum()):,} / {len(rec):,}"
          if not unlocated else f"[FAIL] UNLOCATED        {len(unlocated):,}")
    print(f"[PASS] below local lid  {len(rec)-len(above):,} / {len(rec):,}"
          if not above else f"[FAIL] ABOVE local lid  {len(above):,}")
    if ok.any():
        print(f"       clearance below lid: min {clear.min():.3f} m, "
              f"p01 {np.percentile(clear,1):.3f}, med {np.median(clear):.3f} m")
    bad = len(unlocated) + len(above)
    print()
    print("RECEIVER CONTAINMENT PASS" if bad == 0 else f"RECEIVER CONTAINMENT FAILED ({bad})")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
