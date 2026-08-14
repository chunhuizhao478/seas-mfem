#!/usr/bin/env python3
"""extract_multi.py -- extract ALL surgery tiles from one pass over the mesh.

The per-box extractor scans the whole 90M-tet mesh per box; at ~2,000 tiles that
is a day of redundant IO. This makes ONE chunked pass, assigns every tet to at
most one tile (tiles are disjoint by construction), and then emits each tile's
.msh + meta exactly like extract_patch_medit would have.
"""
import argparse
import json
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections
import meshio

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
REC = np.dtype([('a', np.int32), ('b', np.int32), ('c', np.int32)])
CH = 6_000_000


def rec(A):
    return np.ascontiguousarray(np.sort(A, axis=1).astype(np.int32)).view(REC).ravel()


ap = argparse.ArgumentParser()
ap.add_argument('--mesh', required=True)
ap.add_argument('--boxes', required=True)
ap.add_argument('--outdir', default='build_tmp/tiles')
a = ap.parse_args()

import os
os.makedirs(a.outdir, exist_ok=True)
boxes = json.load(open(a.boxes))
NB = len(boxes)
bx = np.array([[b[0], b[1], b[2], b[3], b[4], b[5]] for b in boxes])
print(f'[tiles] {NB}', flush=True)

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
del S, gtri, gref

# tile assignment key: tiles are an axis-aligned grid in x,y (all same half-size)
H = bx[0, 2]
T = 2 * H
# vectorized assignment: 64-bit tile key + searchsorted (a per-tet Python loop
# at 90M tets would take the better part of an hour by itself)
KMUL = 1_000_000
tk64 = (np.floor(bx[:, 0] / T).astype(np.int64) + 500_000) * KMUL \
     + (np.floor(bx[:, 1] / T).astype(np.int64) + 500_000)
tord = np.argsort(tk64)
tk64s = tk64[tord]
tets = [[] for _ in range(NB)]
tids = [[] for _ in range(NB)]
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    Tt = pd.read_csv(a.mesh, sep=r'\s+', header=None, skiprows=LT + 1 + s, nrows=n,
                     usecols=range(4), dtype=np.int32, engine='c').to_numpy() - 1
    cen = P[Tt].mean(1)
    k64 = (np.floor(cen[:, 0] / T).astype(np.int64) + 500_000) * KMUL \
        + (np.floor(cen[:, 1] / T).astype(np.int64) + 500_000)
    pos = np.searchsorted(tk64s, k64)
    np.clip(pos, 0, NB - 1, out=pos)
    hit = tk64s[pos] == k64
    b = tord[pos]
    zin = hit & (cen[:, 2] > bx[b, 4]) & (cen[:, 2] < bx[b, 5])
    sel = np.flatnonzero(zin)
    bs = b[sel]
    o2 = np.argsort(bs, kind='stable')
    sel, bs = sel[o2], bs[o2]
    cuts = np.searchsorted(bs, np.arange(NB + 1))
    for bb in np.unique(bs):
        sl = sel[cuts[bb]:cuts[bb + 1]]
        tets[bb].append(Tt[sl])
        tids[bb].append(sl.astype(np.int64) + s)
    del Tt, cen, k64, pos, hit, b, zin, sel, bs
    print(f'  ..{s+n:,}/{NT:,}', flush=True)
for bb in range(NB):
    if tets[bb]:
        tets[bb] = [np.vstack(tets[bb])]
        tids[bb] = [np.concatenate(tids[bb])]

nsk = 0
for b in range(NB):
    if not tets[b]:
        print(f'[tile {b}] EMPTY -- skipped')
        continue
    Tl = tets[b][0].astype(np.int64)
    tids[b] = tids[b][0]
    fc = np.concatenate([Tl[:, list(FACE[k])] for k in range(4)])
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
    vid = np.unique(np.concatenate([Tl.ravel(), tri.ravel()]))
    rm = -np.ones(NV, np.int64)
    rm[vid] = np.arange(len(vid))
    Pl, Tll, tril = P[vid], rm[Tl], rm[tri]
    mesh = meshio.Mesh(Pl, [('triangle', tril), ('tetra', Tll)],
                       cell_data={'gmsh:physical': [rf, np.ones(len(Tll), np.int32)],
                                  'gmsh:geometrical': [rf, np.ones(len(Tll), np.int32)]})
    meshio.write(f'{a.outdir}/tile_{b}.msh', mesh, file_format='gmsh22', binary=False)
    np.savez(f'{a.outdir}/tile_{b}_meta.npz', skinP=Pl[np.unique(tril)],
             nfault=int((rf == 101).sum()), gids=vid,
             tids=tids[b], box=bx[b], et=boxes[b][7])
    nsk += 1
print(f'[out] {nsk} tiles written to {a.outdir}/')
