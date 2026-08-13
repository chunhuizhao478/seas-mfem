#!/usr/bin/env python3
"""rv_census.py -- across-fault tet volume ratio (Zhang et al. 2023 JGR Eq.18).

rv = max(V_A/V_B, V_B/V_A) for the tet pair sharing a dynamic-rupture face.
Because both tets share that triangle as their base, rv is identically the ratio
of the two APEX HEIGHTS -- so repairing rv moves apex VERTICES and never touches
the fault triangulation.

Zhang's thresholds: >1.5 SSO onset, >3 severe.  Measured on this project
(PRE_k_2_5_CASE2): deep P(tensile flip) vs the rv<1.5 baseline is 3.0x at rv 2-2.5,
48.8x at rv 3-5, 198.8x above 5 -- a tensile flip clamps strength to zero and the
GP slides forever.

Usage: rv_census.py <mesh.puml.h5> [--csv out.csv] [--top N]
"""
import argparse
import numpy as np, h5py

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
CH = 4_000_000

ap = argparse.ArgumentParser()
ap.add_argument("mesh"); ap.add_argument("--csv"); ap.add_argument("--top", type=int, default=0)
a = ap.parse_args()

with h5py.File(a.mesh) as f:
    V = f["geometry"][:].astype(np.float64)
    T = f["connect"][:].astype(np.int64)
    B = f["boundary"][:].astype(np.int64)
nt = len(T)

tri, own = [], []
for k in range(4):
    m = ((B >> (8 * k)) & 0xFF) == 3
    if m.any():
        idx = np.flatnonzero(m)
        tri.append(T[idx][:, list(FACE[k])]); own.append(idx)
tri = np.vstack(tri); own = np.concatenate(own)
print(f"[mesh] {nt:,} tets   {len(tri):,} DR face instances -> {len(tri)//2:,} faces")

# pair the two instances of each face
key = np.sort(tri, axis=1)
order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
key, own_s, tri_s = key[order], own[order], tri[order]
same = np.all(key[0::2] == key[1::2], axis=1)
if not same.all():
    raise RuntimeError(f"{int((~same).sum())} DR faces are not shared by exactly 2 tets")
tA, tB = own_s[0::2], own_s[1::2]

def vol(idx):
    out = np.empty(len(idx))
    for s in range(0, len(idx), CH):
        p = V[T[idx[s:s + CH]]]
        d = p[:, 1:] - p[:, 0:1]
        out[s:s + CH] = np.abs(np.einsum('ij,ij->i', np.cross(d[:, 0], d[:, 1]), d[:, 2])) / 6.0
    return out

vA, vB = vol(tA), vol(tB)
rv = np.maximum(vA / vB, vB / vA)
cen = V[tri_s[0::2]].mean(1)

n = len(rv)
print(f"\nrv  min {rv.min():.4f}  median {np.median(rv):.4f}  mean {rv.mean():.4f}  max {rv.max():.3f}")
print(f"    p90 {np.percentile(rv,90):.3f}  p99 {np.percentile(rv,99):.3f}  p99.9 {np.percentile(rv,99.9):.3f}")
print("\n threshold |     faces |    share | Zhang")
for t, lab in ((1.5, "SSO onset"), (2.0, "TARGET"), (3.0, "severe"), (5.0, ""), (10.0, "")):
    m = rv > t
    print(f"   rv > {t:<4} | {int(m.sum()):9,} | {100*m.mean():7.4f} % | {lab}")
for lo, hi in ((0, -1500), (-1500, -12000), (-12000, -1e9)):
    m = (cen[:, 2] <= lo) & (cen[:, 2] > hi)
    if m.any():
        print(f"  depth {lo/1000:6.1f}..{max(hi,-40000)/1000:7.1f} km: {int(m.sum()):8,} faces, "
              f"rv>2 {int((rv[m]>2).sum()):7,} ({100*(rv[m]>2).mean():6.3f} %), max {rv[m].max():7.3f}")
if a.csv:
    k = np.flatnonzero(rv > 2.0)
    o = k[np.argsort(-rv[k])]
    np.savetxt(a.csv, np.column_stack([o, rv[o], cen[o]]), delimiter=",",
               header="face_idx,rv,x,y,z", comments="", fmt="%d,%.6f,%.2f,%.2f,%.2f")
    print(f"\n[csv] {len(o):,} faces with rv > 2 -> {a.csv}")
if a.top:
    o = np.argsort(-rv)[:a.top]
    print(f"\n worst {a.top}:")
    for i in o:
        print(f"   rv {rv[i]:9.3f}  at ({cen[i,0]:.0f}, {cen[i,1]:.0f}, {cen[i,2]:.0f})")
