#!/usr/bin/env python3
"""check_fault_in_domain.py -- does the ALT fault fit inside the trimmed box?

The domain is a ROTATED rectangle (azimuth 130 deg) that was trimmed on the
southeast to stay inside MUSCAL's footprint, so an axis-aligned bbox test would
both over- and under-state the fit. Test in the rotated (s,t) frame the build
actually uses:

    u = (S-W)/|S-W|   L1 = 520.000 km
    v = (N-W)/|N-W|   L2 = 300.009 km
    s = (p-W).u   must lie in [0, L1]
    t = (p-W).v   must lie in [0, L2]

Reports the fault's extent in that frame and its clearance to each of the four
walls. A fault that touches or crosses a wall would make the PLC unbuildable --
the fault has to be an interior crack, not a boundary.

Memory-light on purpose: only fault vertices are read, never the volume.
"""
import argparse

import h5py
import numpy as np

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
CH = 8_000_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--corners", required=True)
    ap.add_argument("--depth", type=float, default=80000.0)
    a = ap.parse_args()

    W, N, E, S = np.load(a.corners)
    u = (S - W) / np.linalg.norm(S - W)
    v = (N - W) / np.linalg.norm(N - W)
    L1 = float(np.linalg.norm(S - W))
    L2 = float(np.linalg.norm(N - W))
    print(f"[domain] L1 {L1/1000:,.3f} km along u, L2 {L2/1000:,.3f} km along v, "
          f"depth {a.depth/1000:,.0f} km")

    with h5py.File(a.mesh) as f:
        V = f["geometry"][:].astype(np.float64)
        C = f["connect"]
        B = f["boundary"][:].astype(np.int64)
        vid = []
        for s0 in range(0, C.shape[0], CH):
            c = C[s0:s0 + CH][:].astype(np.int64)
            b = B[s0:s0 + CH]
            for k in range(4):
                m = ((b >> (8 * k)) & 0xFF) == 3
                if m.any():
                    vid.append(np.unique(c[m][:, list(FACE[k])]))
    fv = np.unique(np.concatenate(vid))
    P = V[fv]
    print(f"[fault] {len(P):,} distinct fault vertices")
    print(f"  UTM  E {P[:,0].min():,.1f} .. {P[:,0].max():,.1f}   "
          f"N {P[:,1].min():,.1f} .. {P[:,1].max():,.1f}   "
          f"z {P[:,2].min():,.1f} .. {P[:,2].max():,.1f}")

    d = P[:, :2] - W
    s = d @ u
    t = d @ v
    print(f"\n[rotated frame]  s {s.min():,.1f} .. {s.max():,.1f} m   "
          f"(domain 0 .. {L1:,.1f})")
    print(f"                 t {t.min():,.1f} .. {t.max():,.1f} m   "
          f"(domain 0 .. {L2:,.1f})")
    print(f"\n[clearance to walls]")
    print(f"  s = 0    (NW end)  {s.min():>12,.1f} m")
    print(f"  s = L1   (SE end)  {L1 - s.max():>12,.1f} m")
    print(f"  t = 0    (SW side) {t.min():>12,.1f} m")
    print(f"  t = L2   (NE side) {L2 - t.max():>12,.1f} m")
    print(f"  bottom z = -{a.depth/1000:,.0f} km  {a.depth + P[:,2].min():>12,.1f} m")
    print(f"  top    z = 0       {-P[:,2].max():>12,.1f} m")

    inside = (s >= 0) & (s <= L1) & (t >= 0) & (t <= L2) & \
             (P[:, 2] >= -a.depth) & (P[:, 2] <= 0)
    n_out = int((~inside).sum())
    print(f"\n[verdict] fault vertices inside the domain: "
          f"{int(inside.sum()):,} / {len(P):,}   OUTSIDE: {n_out:,}")
    if n_out:
        o = P[~inside]
        print(f"  !! offenders bbox E {o[:,0].min():,.1f}..{o[:,0].max():,.1f} "
              f"N {o[:,1].min():,.1f}..{o[:,1].max():,.1f} z {o[:,2].min():,.1f}..{o[:,2].max():,.1f}")
    else:
        print("  the fault is a strictly INTERIOR crack -- PLC is buildable")


if __name__ == "__main__":
    main()
