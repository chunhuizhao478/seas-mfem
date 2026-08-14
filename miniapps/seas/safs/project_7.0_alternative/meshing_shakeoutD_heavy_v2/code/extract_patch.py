#!/usr/bin/env python3
"""extract_patch.py -- cut a box submesh out of a PUML as a MEDIT testbed.

Purpose: A/B repair candidates on a REAL piece of the production mesh (its
actual sliver population, its actual frozen fault) at ~1M-tet scale, minutes per
run, before committing to an 87M-tet pass. The repo's own lesson: smoke-test the
repairer on a small mesh of the same family first.

The patch's SKIN is frozen so a repaired patch can later be stitched back:
  * every face of a selected tet whose neighbour is NOT selected and which has
    no BC code is an INTERFACE face -> ref 999, RequiredTriangles
  * real surfaces keep their SAFS refs (fault 101 / top 102 / absorbing 104),
    also Required
  * fault faces interior to the patch (both sides selected) are emitted ONCE
    (an interior crack appears on both tets) with ref 101, Required
With -nosurf + Required, mmg may not move any of these vertices, so the patch
boundary stays bit-compatible with the surrounding mesh.
"""
import argparse
import numpy as np, h5py, pandas as pd

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
CH = 4_000_000

ap = argparse.ArgumentParser()
ap.add_argument('--mesh', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--meta', required=True)
ap.add_argument('--cx', type=float, required=True)
ap.add_argument('--cy', type=float, required=True)
ap.add_argument('--hx', type=float, default=3000.0)
ap.add_argument('--hy', type=float, default=3000.0)
ap.add_argument('--z0', type=float, default=-19500.0)
ap.add_argument('--z1', type=float, default=-14500.0)
a = ap.parse_args()

f = h5py.File(a.mesh)
V = f['geometry'][:]
C = f['connect']
B = f['boundary']
NT = C.shape[0]

selT, selB = [], []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    b = B[s:s + CH][:].astype(np.int64)
    cen = V[t].mean(1)
    m = ((np.abs(cen[:, 0] - a.cx) < a.hx) & (np.abs(cen[:, 1] - a.cy) < a.hy)
         & (cen[:, 2] > a.z0) & (cen[:, 2] < a.z1))
    if m.any():
        selT.append(t[m]); selB.append(b[m])
    del t, b, cen
T = np.vstack(selT)
Bb = np.concatenate(selB)
del selT, selB
print(f'[patch] {len(T):,} tets selected', flush=True)

# faces + their BC codes (0 = plain interior face of the global mesh)
fc = np.concatenate([T[:, list(FACE[k])] for k in range(4)])
code = np.concatenate([((Bb >> (8 * k)) & 0xFF) for k in range(4)]).astype(np.int32)
key = np.sort(fc, axis=1)
u, first, inv, cnt = np.unique(key, axis=0, return_index=True,
                               return_inverse=True, return_counts=True)
agg = np.zeros(len(u), np.int32)
np.maximum.at(agg, inv, code)
ref = np.where(agg == 3, 101, np.where(agg == 1, 102, np.where(agg == 5, 104, 0)))
ref = np.where((cnt == 1) & (ref == 0), 999, ref)     # skin
keep = ref > 0
tri = fc[first][keep]
rf = ref[keep]
for r in (101, 102, 104, 999):
    n = int((rf == r).sum())
    if n:
        print(f'  ref {r}: {n:,} triangles')

# renumber
vid = np.unique(np.concatenate([T.ravel(), tri.ravel()]))
rm = -np.ones(len(V), np.int64)
rm[vid] = np.arange(len(vid))
P = V[vid]
Tl = rm[T]
tril = rm[tri]
skinv = np.unique(tril)
print(f'[patch] {len(P):,} verts, {len(skinv):,} on the frozen skin')

with open(a.out, 'w') as fo:
    fo.write('MeshVersionFormatted 2\nDimension 3\n\nVertices\n%d\n' % len(P))
    pd.DataFrame({0: P[:, 0], 1: P[:, 1], 2: P[:, 2], 3: 0}).to_csv(
        fo, sep=' ', header=False, index=False, float_format='%.10g')
    fo.write('\nTriangles\n%d\n' % len(tril))
    pd.DataFrame({0: tril[:, 0] + 1, 1: tril[:, 1] + 1, 2: tril[:, 2] + 1,
                  3: rf}).to_csv(fo, sep=' ', header=False, index=False)
    fo.write('\nRequiredTriangles\n%d\n' % len(tril))
    pd.DataFrame({0: np.arange(1, len(tril) + 1)}).to_csv(
        fo, sep=' ', header=False, index=False)
    fo.write('\nTetrahedra\n%d\n' % len(Tl))
    pd.DataFrame({0: Tl[:, 0] + 1, 1: Tl[:, 1] + 1, 2: Tl[:, 2] + 1,
                  3: Tl[:, 3] + 1, 4: np.ones(len(Tl), np.int32)}).to_csv(
        fo, sep=' ', header=False, index=False)
    fo.write('\nEnd\n')
np.savez(a.meta, skinP=P[skinv], nfault=int((rf == 101).sum()),
         gids=vid, box=[a.cx, a.cy, a.hx, a.hy, a.z0, a.z1])
print(f'[out] {a.out}  +  {a.meta}')
