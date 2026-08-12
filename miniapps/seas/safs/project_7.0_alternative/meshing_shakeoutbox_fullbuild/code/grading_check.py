#!/usr/bin/env python3
"""grading_check.py -- how abrupt are the element-size changes, and WHERE?

The stated reason to rebuild from scratch is to get a gradual size transition
everywhere. Whether that is needed is measurable on the deployed mesh: for every
INTERIOR face, take the ratio of the two neighbouring cells' max edges. A
conforming tet mesh built from a Lipschitz-limited field sits near 1; a weld
between two independently-built blocks, or a red-refinement boundary, shows a
tail at 2 and beyond.

Reports the ratio distribution globally, and separately for faces near the
frozen-parent WELD (the collar/parent interface, which is the seam a rebuild
would remove) so the two causes are not confused.

Face-neighbour pairing is done by sorting canonical vertex triples -- O(n log n)
and chunk-free on the key array, which at 10^8 cells is the only part that must
be watched (4 x nt keys).
"""
import argparse

import h5py
import numpy as np

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EI = np.array([p[0] for p in PAIRS])
EJ = np.array([p[1] for p in PAIRS])
CH = 2_000_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--parent-tets", type=int, default=0,
                    help="tets [0,N) came from the frozen parent; faces straddling "
                         "the boundary are the WELD")
    ap.add_argument("--sample", type=int, default=1)
    a = ap.parse_args()

    with h5py.File(a.mesh, "r") as f:
        G = f["geometry"][:]
        nt = f["connect"].shape[0]
        conn = f["connect"]
        hmax = np.empty(nt, np.float32)
        keys = np.empty((nt * 4, 3), np.int64)
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            P = G[T]
            hmax[s0:s1] = np.linalg.norm(P[:, EI] - P[:, EJ], axis=2).max(1)
            for s, fv in enumerate(LOCAL_FACES):
                tri = np.sort(T[:, list(fv)], axis=1)
                keys[(s0 * 4 + s * (s1 - s0)):(s0 * 4 + (s + 1) * (s1 - s0))] = tri
            del T, P
            if (s0 // CH) % 20 == 0:
                print(f"    ... {s1:,}/{nt:,}", flush=True)
    print(f"[mesh] {nt:,} tets", flush=True)

    owner = np.repeat(np.arange(nt, dtype=np.int64), 4)
    # rebuild owner ordering to match how keys were filled (tet-major per chunk)
    owner = np.empty(nt * 4, np.int64)
    for s0 in range(0, nt, CH):
        s1 = min(s0 + CH, nt)
        n = s1 - s0
        for s in range(4):
            owner[(s0 * 4 + s * n):(s0 * 4 + (s + 1) * n)] = np.arange(s0, s1)

    print("[faces] sorting…", flush=True)
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
    keys = keys[order]
    owner = owner[order]
    same = np.all(keys[1:] == keys[:-1], axis=1)
    i = np.nonzero(same)[0]
    A, B = owner[i], owner[i + 1]
    del keys, owner, order, same
    print(f"[faces] {len(A):,} interior faces", flush=True)

    hi = np.maximum(hmax[A], hmax[B])
    lo = np.minimum(hmax[A], hmax[B])
    r = hi / np.maximum(lo, 1e-9)
    print("\n  ADJACENT-CELL SIZE RATIO (max edge, larger/smaller), ALL interior faces")
    for q in (50, 90, 99, 99.9):
        print(f"    p{q:<6} {np.percentile(r, q):.3f}")
    print(f"    max     {r.max():.3f}")
    for t in (1.5, 2.0, 3.0, 4.0):
        print(f"    ratio > {t}: {int((r > t).sum()):,} "
              f"({100*(r>t).mean():.4f} %)")

    if a.parent_tets:
        weld = ((A < a.parent_tets) != (B < a.parent_tets))
        print(f"\n  WELD faces (parent block <-> collar block): {int(weld.sum()):,}")
        if weld.any():
            rw = r[weld]
            for q in (50, 90, 99):
                print(f"    p{q:<6} {np.percentile(rw, q):.3f}")
            print(f"    max     {rw.max():.3f}")
            for t in (1.5, 2.0, 3.0):
                print(f"    ratio > {t}: {int((rw > t).sum()):,} "
                      f"({100*(rw>t).mean():.4f} %)")
        interior = ~weld
        print(f"\n  NON-weld interior faces: {int(interior.sum()):,}")
        for q in (99, 99.9):
            print(f"    p{q:<6} {np.percentile(r[interior], q):.3f}")
        print(f"    max     {r[interior].max():.3f}")


if __name__ == "__main__":
    main()
