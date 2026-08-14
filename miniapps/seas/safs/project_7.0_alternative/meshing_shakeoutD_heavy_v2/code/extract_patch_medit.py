#!/usr/bin/env python3
"""extract_patch_medit.py -- box submesh from a MEDIT mesh, written as gmsh22 .msh.

Same job as extract_patch.py but for the mid-pipeline MEDIT state, and emitting
the .msh that remove_sliver_tets.py consumes. Surface refs are carried over from
the global Triangles section; skin faces (patch boundary that is interior to the
global mesh) get ref 999. Pass --fault-tags 101,999 downstream so BOTH are
protected -- the fault because it is the fault, the skin so the patch stitches
back without touching its surroundings.
"""
import argparse
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
REC = np.dtype([('a', np.int32), ('b', np.int32), ('c', np.int32)])
CH = 8_000_000


def rec(A):
    return np.ascontiguousarray(np.sort(A, axis=1).astype(np.int32)).view(REC).ravel()


ap = argparse.ArgumentParser()
ap.add_argument('--mesh', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--meta', required=True)
ap.add_argument('--cx', type=float, required=True)
ap.add_argument('--cy', type=float, required=True)
ap.add_argument('--hx', type=float, default=2500.0)
ap.add_argument('--hy', type=float, default=2500.0)
ap.add_argument('--z0', type=float, required=True)
ap.add_argument('--z1', type=float, required=True)
a = ap.parse_args()

sec = medit_sections(a.mesh)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
LS, NS = sec['Triangles']
P = pd.read_csv(a.mesh, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
S = pd.read_csv(a.mesh, sep=r'\s+', header=None, skiprows=LS + 1, nrows=NS,
                usecols=range(4), dtype=np.int64, engine='c').to_numpy()
gtri, gref = S[:, :3] - 1, S[:, 3].astype(np.int32)
gk = rec(gtri)
o = np.argsort(gk, kind='stable')
gk, gref_s = gk[o], gref[o]

sel = []
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    T = pd.read_csv(a.mesh, sep=r'\s+', header=None, skiprows=LT + 1 + s, nrows=n,
                    usecols=range(4), dtype=np.int32, engine='c').to_numpy() - 1
    cen = P[T].mean(1)
    m = ((np.abs(cen[:, 0] - a.cx) < a.hx) & (np.abs(cen[:, 1] - a.cy) < a.hy)
         & (cen[:, 2] > a.z0) & (cen[:, 2] < a.z1))
    if m.any():
        sel.append(T[m])
    del T, cen
T = np.vstack(sel)
del sel
print(f'[patch] {len(T):,} tets', flush=True)

fc = np.concatenate([T[:, list(FACE[k])] for k in range(4)])
key = rec(fc)
u, first, inv, cnt = np.unique(key, return_index=True, return_inverse=True,
                               return_counts=True)
pos = np.searchsorted(gk, u)
np.clip(pos, 0, len(gk) - 1, out=pos)
hit = gk[pos] == u
ref = np.where(hit, gref_s[pos], 0).astype(np.int32)
ref = np.where((cnt == 1) & (ref == 0), 999, ref)
keep = ref > 0
tri = fc[first][keep]
rf = ref[keep]
for r in (101, 102, 104, 999):
    n = int((rf == r).sum())
    if n:
        print(f'  ref {r}: {n:,}')

vid = np.unique(np.concatenate([T.ravel(), tri.ravel()]))
rm = -np.ones(NV, np.int64)
rm[vid] = np.arange(len(vid))
Pl, Tl, tril = P[vid], rm[T], rm[tri]
skinv = np.unique(tril)
print(f'[patch] {len(Pl):,} verts, {len(skinv):,} frozen-skin verts')

import meshio
mesh = meshio.Mesh(Pl, [('triangle', tril), ('tetra', Tl)],
                   cell_data={'gmsh:physical': [rf, np.ones(len(Tl), np.int32)],
                              'gmsh:geometrical': [rf, np.ones(len(Tl), np.int32)]})
meshio.write(a.out, mesh, file_format='gmsh22', binary=False)
np.savez(a.meta, skinP=Pl[skinv], nfault=int((rf == 101).sum()), gids=vid,
         box=[a.cx, a.cy, a.hx, a.hy, a.z0, a.z1])
print(f'[out] {a.out}  +  {a.meta}')
