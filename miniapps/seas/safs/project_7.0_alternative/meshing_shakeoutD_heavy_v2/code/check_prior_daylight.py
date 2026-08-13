#!/usr/bin/env python3
"""check_prior_daylight.py -- what did the DEPLOYED ALT mesh actually do at the top?

"Snap the fault up, or leave it buried?" should be answered from the record, not
from taste. The record is the deployed small-domain heavy mesh: for each vertex
on the fault's boundary, is that vertex ALSO a free-surface (BC 1) vertex?

  shared with the free surface  -> the fault DAYLIGHTS there (conformal trace)
  fault-only, near z = 0        -> the fault is BURIED there (blind tip)

The original builder had an explicit per-strand policy (`--daylight-faults`:
"daylight only the STEEP strands (clean wedge) and keep a SHALLOW thrust blind
... avoids the flat wedge tets that force a size-map remesh + mesh bloat"), so
the question is which side of that policy the SE stretch fell on.

Reports the split, and where the buried part is, so the ShakeOut-D build can
reproduce the same choice rather than silently change the physics.
"""
import argparse

import h5py
import numpy as np

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
CH = 8_000_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    a = ap.parse_args()

    ftri, top_v = [], []
    with h5py.File(a.mesh) as f:
        V = f["geometry"][:].astype(np.float64)
        C = f["connect"]
        B = f["boundary"][:].astype(np.int64)
        for s0 in range(0, C.shape[0], CH):
            c = C[s0:s0 + CH][:].astype(np.int64)
            b = B[s0:s0 + CH]
            for k in range(4):
                code = (b >> (8 * k)) & 0xFF
                m3 = code == 3
                if m3.any():
                    ftri.append(c[m3][:, list(FACE[k])])
                m1 = code == 1
                if m1.any():
                    top_v.append(np.unique(c[m1][:, list(FACE[k])]))
    F = np.vstack(ftri)
    key = np.sort(F, axis=1)
    o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    T = F[o][0::2]
    topv = np.unique(np.concatenate(top_v))
    print(f"[mesh] fault {len(T):,} unique triangles; free surface {len(topv):,} vertices")

    e = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [0, 2]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    bnd = u[c == 1]
    bv = np.unique(bnd)
    zb = V[bv, 2]
    print(f"[fault boundary] {len(bv):,} vertices, z {zb.min():,.2f} .. {zb.max():,.2f}")

    isday = np.isin(bv, topv)
    shallow = zb > -200.0
    print(f"\n[daylight test] fault-boundary vertices that are ALSO free-surface vertices:")
    print(f"  daylighting        {int(isday.sum()):>8,}")
    print(f"  not daylighting    {int((~isday).sum()):>8,}")
    print(f"  ... of those, shallow (z > -200 m): {int((~isday & shallow).sum()):>8,}")

    bur = bv[~isday & shallow]
    if len(bur):
        Pb = V[bur]
        print(f"\n[BURIED SHALLOW TIP] {len(bur):,} vertices")
        print(f"  z    {Pb[:,2].min():,.2f} .. {Pb[:,2].max():,.2f} m")
        print(f"  E    {Pb[:,0].min():,.0f} .. {Pb[:,0].max():,.0f}")
        print(f"  N    {Pb[:,1].min():,.0f} .. {Pb[:,1].max():,.0f}")
        h, edg = np.histogram(Pb[:, 2], bins=[-200, -100, -60, -50, -25, -10, -5, -1, 0])
        for i in range(len(h)):
            if h[i]:
                print(f"    z {edg[i]:>7,.0f} .. {edg[i+1]:<7,.0f} : {h[i]:>6,}")
    dv = bv[isday]
    if len(dv):
        print(f"\n[DAYLIGHTING TRACE] {len(dv):,} vertices, "
              f"z {V[dv,2].min():,.3f} .. {V[dv,2].max():,.3f} m")


if __name__ == "__main__":
    main()
