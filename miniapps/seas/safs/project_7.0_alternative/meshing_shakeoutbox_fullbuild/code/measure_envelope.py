#!/usr/bin/env python3
"""measure_envelope.py -- recover the DEPLOYED size design from the shipped mesh.

A from-scratch rebuild has to reproduce the fault-normal refinement that the
current heavy mesh actually has, and the honest way to get that spec is to
MEASURE it rather than to re-derive it from a remembered "200 m within 1 km".
The deployed mesh is the result of a base fb200 build plus two red refinements
plus a gate campaign, so its real envelope is not any single documented rule.

Reports element size (max edge) as a function of distance from the fault sheet,
which is the h(d) envelope a new build must match, plus the same split by depth
so the free-surface behaviour is visible separately.

Distance is to the fault TRIANGLE CENTROIDS via a KD tree, always with
`distance_upper_bound` and `workers=-1` -- an unbounded single-threaded query on
a 2.5M-triangle sheet costs ~4 ms each, i.e. days over 10^8 cells.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_gate" / "code"))
LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EI = np.array([p[0] for p in PAIRS])
EJ = np.array([p[1] for p in PAIRS])
BC_DR = 3
CH = 1_000_000


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--band-max", type=float, default=12000.0,
                    help="report out to this distance; beyond it the far field "
                         "is set by the frequency gate, not by the fault")
    ap.add_argument("--sample", type=int, default=6,
                    help="take every Nth cell (the envelope is smooth; this is a "
                         "shape measurement, not a census)")
    a = ap.parse_args()

    with h5py.File(a.mesh, "r") as f:
        G = f["geometry"][:]
        nt = f["connect"].shape[0]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        # fault triangle centroids
        cen = []
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                hit |= face_code(Bc, s) == BC_DR
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            for s in range(4):
                m = face_code(Bs, s) == BC_DR
                if m.any():
                    cen.append(G[Tc[m][:, list(LOCAL_FACES[s])]].mean(1))
        C = np.vstack(cen)
    print(f"[fault] {len(C):,} DR facet centroids", flush=True)
    tree = cKDTree(C)

    d_all, h_all, z_all = [], [], []
    with h5py.File(a.mesh, "r") as f:
        conn = f["connect"]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)[::a.sample]
            if not len(T):
                continue
            P = G[T]
            bary = P.mean(1)
            hmax = np.linalg.norm(P[:, EI] - P[:, EJ], axis=2).max(1)
            d, _ = tree.query(bary, distance_upper_bound=a.band_max,
                              workers=-1)
            d_all.append(d.astype(np.float32))
            h_all.append(hmax.astype(np.float32))
            z_all.append(bary[:, 2].astype(np.float32))
            del T, P
            if (s0 // CH) % 30 == 0:
                print(f"    ... {s1:,}/{nt:,}", flush=True)
    D = np.concatenate(d_all)
    H = np.concatenate(h_all)
    Z = np.concatenate(z_all)
    fin = np.isfinite(D)
    print(f"\nsampled {len(D):,} cells; {int(fin.sum()):,} within "
          f"{a.band_max:,.0f} m of the fault\n")

    edges = [0, 100, 200, 400, 700, 1000, 1500, 2000, 3000, 5000, 8000, 12000]
    print("  ELEMENT SIZE vs DISTANCE FROM FAULT (max edge, m)")
    print(f"    {'d band [m]':<16}{'n':>12}{'p5':>9}{'p50':>9}{'p95':>9}{'max':>9}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = fin & (D >= lo) & (D < hi)
        if not m.any():
            continue
        h = H[m]
        print(f"    {lo:>6,}-{hi:<9,}{int(m.sum()):>12,}{np.percentile(h,5):>9,.0f}"
              f"{np.median(h):>9,.0f}{np.percentile(h,95):>9,.0f}{h.max():>9,.0f}")
    far = ~fin
    if far.any():
        h = H[far]
        print(f"    {'>'+format(a.band_max,',.0f'):<16}{int(far.sum()):>12,}"
              f"{np.percentile(h,5):>9,.0f}{np.median(h):>9,.0f}"
              f"{np.percentile(h,95):>9,.0f}{h.max():>9,.0f}")

    print("\n  NEAR-FAULT (d < 1 km) size vs DEPTH")
    print(f"    {'z band [m]':<20}{'n':>12}{'p50 h':>9}{'p95 h':>9}")
    zed = [0, -1000, -3000, -6000, -12000, -20000, -40001]
    nf = fin & (D < 1000.0)
    for lo, hi in zip(zed[:-1], zed[1:]):
        m = nf & (Z <= lo) & (Z > hi)
        if m.any():
            print(f"    {lo:>7,} .. {hi:<10,}{int(m.sum()):>12,}"
                  f"{np.median(H[m]):>9,.0f}{np.percentile(H[m],95):>9,.0f}")

    print("\n  FAR FIELD (d > 5 km) size vs DEPTH -- this is the gate's territory")
    ff = (~fin) | (D > 5000.0)
    for lo, hi in zip(zed[:-1], zed[1:]):
        m = ff & (Z <= lo) & (Z > hi)
        if m.any():
            print(f"    {lo:>7,} .. {hi:<10,}{int(m.sum()):>12,}"
                  f"{np.median(H[m]):>9,.0f}{np.percentile(H[m],95):>9,.0f}")


if __name__ == "__main__":
    main()
