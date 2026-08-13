#!/usr/bin/env python3
"""rv_fixability.py -- how much of the rv>target population is even MOVABLE?

rv = h_A/h_B over the shared fault triangle, so the only lever is the position of
the two APEX vertices.  `fix_rv_apex_*.py` pins every vertex carried by a boundary
face (BC 1 free surface, 3 fault, 5 absorbing) -- moving those would change the
free surface or the dynamic-rupture surface itself.  A face whose BOTH apexes are
pinned therefore cannot be repaired by vertex motion at any effort.

Reports the ceiling before any optimisation is spent.
"""
import argparse
import numpy as np, h5py

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
BC_PINNED = (1, 3, 5)
ap = argparse.ArgumentParser()
ap.add_argument("mesh"); ap.add_argument("--target", type=float, default=2.0)
a = ap.parse_args()

with h5py.File(a.mesh) as f:
    V = f["geometry"][:].astype(np.float64); T = f["connect"][:].astype(np.int64)
    B = f["boundary"][:].astype(np.int64)
nV = len(V)

pinned = np.zeros(nV, bool)
for k in range(4):
    fc = (B >> (8 * k)) & 0xFF
    for bc in BC_PINNED:
        m = fc == bc
        if m.any():
            pinned[np.unique(T[m][:, list(FACE[k])])] = True
print(f"[pinned] {int(pinned.sum()):,} of {nV:,} vertices ({100*pinned.mean():.2f} %) "
      f"lie on a free-surface / fault / absorbing face")

tri, own, apex = [], [], []
for k in range(4):
    m = ((B >> (8 * k)) & 0xFF) == 3
    if m.any():
        idx = np.flatnonzero(m)
        tri.append(T[idx][:, list(FACE[k])]); own.append(idx)
        ap_local = [v for v in range(4) if v not in FACE[k]][0]
        apex.append(T[idx][:, ap_local])
tri = np.vstack(tri); own = np.concatenate(own); apex = np.concatenate(apex)

key = np.sort(tri, axis=1)
o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
key, own, apex, tri = key[o], own[o], apex[o], tri[o]
tA, tB, aA, aB = own[0::2], own[1::2], apex[0::2], apex[1::2]

def vol(idx):
    p = V[T[idx]]; d = p[:, 1:] - p[:, 0:1]
    return np.abs(np.einsum('ij,ij->i', np.cross(d[:, 0], d[:, 1]), d[:, 2])) / 6.0
vA, vB = vol(tA), vol(tB)
rv = np.maximum(vA / vB, vB / vA)
cen = V[tri[0::2]].mean(1)

bad = rv > a.target
pA, pB = pinned[aA], pinned[aB]
both = bad & pA & pB
one = bad & (pA ^ pB)
free = bad & ~pA & ~pB
print(f"\nrv > {a.target}: {int(bad.sum()):,} faces")
print(f"  BOTH apexes pinned  : {int(both.sum()):8,} ({100*both.sum()/bad.sum():5.2f} %)  -- UNFIXABLE by vertex motion")
print(f"  ONE apex free       : {int(one.sum()):8,} ({100*one.sum()/bad.sum():5.2f} %)")
print(f"  BOTH apexes free    : {int(free.sum()):8,} ({100*free.sum()/bad.sum():5.2f} %)")
print(f"  => movable ceiling  : {int((one|free).sum()):8,} ({100*(one|free).sum()/bad.sum():5.2f} %)")
print("\n  by depth (unfixable share):")
for lo, hi, lab in ((0, -1500, "  0 .. -1.5 km"), (-1500, -12000, "-1.5 .. -12 km"),
                    (-12000, -1e9, " -12 .. -40 km")):
    m = bad & (cen[:, 2] <= lo) & (cen[:, 2] > hi)
    if m.any():
        print(f"   {lab}: {int(m.sum()):8,} bad, {int((m&both).sum()):7,} unfixable "
              f"({100*(m&both).sum()/m.sum():5.2f} %), max rv {rv[m].max():8.3f}")
w = np.argsort(-rv)[:8]
print("\n  worst 8:")
for i in w:
    s = "BOTH PINNED" if (pA[i] and pB[i]) else ("one free" if (pA[i] ^ pB[i]) else "both free")
    print(f"   rv {rv[i]:9.3f}  z {cen[i,2]:8.1f} m  {s}")
