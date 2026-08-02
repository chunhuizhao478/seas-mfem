#!/usr/bin/env python3
"""check_fault_identity.py -- prove the refined mesh's fault surface is the SAME
SURFACE as the parent's, and that the PUML BC contract still round-trips.

Every fault-referenced deck input (stress nc, friction nc, the [rs_muw] LuaMap,
nucleation, pickpoints) is evaluated at fault facet positions, so a refinement
that perturbed the fault -- even slightly -- would silently change the physics.

Checks:
  F1  fault triangle MULTISET is identical (same triangles, same coordinates)
  F2  total fault area identical
  F3  every fault face is INTERIOR (appears exactly twice); boundary faces once
  F4  free surface is exactly flat
  F5  0 inverted tets
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from puml_io import LOCAL_FACES, BC_DYNAMIC_RUPTURE  # noqa: E402

CH = 4_000_000


def load(path):
    f = h5py.File(path, "r")
    P = f["geometry"][:]
    ds = f["connect"]
    nt = ds.shape[0]
    C = np.empty((nt, 4), np.int32)
    for s0 in range(0, nt, 8_000_000):
        C[s0:s0 + 8_000_000] = ds[s0:s0 + 8_000_000].astype(np.int32)
    B = f["boundary"][:].astype(np.int32)
    f.close()
    return P, C, B


def codes_of(B):
    bu = np.ascontiguousarray(B, np.int32).view(np.uint32)
    return np.stack([((bu >> np.uint32(8 * s)) & np.uint32(0xFF)).astype(np.uint8)
                     for s in range(4)], 1)


def fault_tris(P, C, B):
    """Sorted-by-coordinate fault triangle vertex coordinates, (n,3,3)."""
    cd = codes_of(B)
    out = []
    for s in range(4):
        sel = np.nonzero(cd[:, s] == BC_DYNAMIC_RUPTURE)[0]
        if sel.size:
            out.append(P[C[sel][:, LOCAL_FACES[s]]])
    T = np.concatenate(out, 0)
    T = np.sort(T.reshape(len(T), 3, 3).view([('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
                .reshape(len(T), 3), axis=1)
    key = np.lexsort((T['z'][:, 2], T['y'][:, 2], T['x'][:, 2],
                      T['z'][:, 1], T['y'][:, 1], T['x'][:, 1],
                      T['z'][:, 0], T['y'][:, 0], T['x'][:, 0]))
    return T[key]


def area(T):
    a = np.stack([T['x'][:, 0], T['y'][:, 0], T['z'][:, 0]], 1)
    b = np.stack([T['x'][:, 1], T['y'][:, 1], T['z'][:, 1]], 1)
    c = np.stack([T['x'][:, 2], T['y'][:, 2], T['z'][:, 2]], 1)
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--new", required=True)
    a = ap.parse_args()
    fails = []

    def chk(n, ok, d=""):
        print(f"[{'PASS' if ok else 'FAIL'}] {n}" + (f"  {d}" if d else ""))
        if not ok:
            fails.append(n)

    Pp, Cp, Bp = load(a.parent)
    Tp = fault_tris(Pp, Cp, Bp)
    ap_area = area(Tp)
    del Pp, Cp, Bp
    Pn, Cn, Bn = load(a.new)
    Tn = fault_tris(Pn, Cn, Bn)
    an_area = area(Tn)

    chk("F1 fault triangle multiset identical",
        Tp.shape == Tn.shape and np.array_equal(Tp.view('f8'), Tn.view('f8')),
        f"{len(Tp):,} vs {len(Tn):,} triangles")
    chk("F2 fault area identical", abs(ap_area - an_area) < 1e-6,
        f"{ap_area/1e6:.9f} vs {an_area/1e6:.9f} km^2  (delta {abs(ap_area-an_area):.3e} m^2)")

    # F3 BC round-trip: fault faces interior (x2), boundary faces (x1)
    cd = codes_of(Bn)
    nt = len(Cn)
    keys, cods = [], []
    for s in range(4):
        tri = np.sort(Cn[:, LOCAL_FACES[s]].astype(np.int64), axis=1)
        keys.append(tri)
        cods.append(cd[:, s])
    K = np.concatenate(keys, 0)
    CO = np.concatenate(cods, 0)
    ordr = np.lexsort((K[:, 2], K[:, 1], K[:, 0]))
    K, CO = K[ordr], CO[ordr]
    same = np.all(K[1:] == K[:-1], axis=1)
    start = np.concatenate([[0], np.nonzero(~same)[0] + 1])
    cnt = np.diff(np.concatenate([start, [len(K)]]))
    first = CO[start]
    n_fault_faces = int((first == BC_DYNAMIC_RUPTURE).sum())
    ok_fault = bool(np.all(cnt[first == BC_DYNAMIC_RUPTURE] == 2))
    ok_bnd = bool(np.all(cnt[(first != 0) & (first != BC_DYNAMIC_RUPTURE)] == 1))
    chk("F3 fault faces interior (x2)", ok_fault, f"{n_fault_faces:,} unique fault faces")
    chk("F3 boundary faces appear once", ok_bnd)

    topv = np.unique(np.concatenate(
        [Cn[cd[:, s] == 1][:, LOCAL_FACES[s]].ravel() for s in range(4)]))
    zr = (float(Pn[topv, 2].min()), float(Pn[topv, 2].max()))
    chk("F4 free surface exactly flat", max(abs(zr[0]), abs(zr[1])) < 1e-6,
        f"z in [{zr[0]:.9f}, {zr[1]:.9f}], {len(topv):,} verts")

    ninv = 0
    for s0 in range(0, nt, CH):
        T = Cn[s0:s0 + CH].astype(np.int64)
        P4 = Pn[T]
        v = np.einsum('ij,ij->i', P4[:, 1] - P4[:, 0],
                      np.cross(P4[:, 2] - P4[:, 0], P4[:, 3] - P4[:, 0])) / 6.0
        ninv += int((v < 0).sum())
    chk("F5 no inverted tets", ninv == 0, f"{ninv} inverted")

    print()
    print("ALL FAULT-IDENTITY CHECKS PASS" if not fails else f"FAILED: {fails}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
