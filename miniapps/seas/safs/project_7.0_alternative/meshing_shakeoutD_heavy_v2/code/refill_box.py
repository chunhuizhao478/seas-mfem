#!/usr/bin/env python3
"""refill_box.py -- the surgery operator of last resort: cavity RE-FILL.

Both splice-safe operators refused the fault-locked residue (9 mmg -optim passes:
zero progress at 435; remove_sliver_tets: "NO splice-safe remesh found" on every
worst cluster). But the small-domain mesh carries the IDENTICAL fault facets with
zero sub-0.05 tets, so a valid connectivity exists -- the residue is a local
minimum, not geometry. A local operation that must keep most of the old
connectivity cannot escape it; throwing the box's interior away entirely can.

Method: take the extracted box (.msh from extract_patch_medit), keep ONLY its
surface triangles -- the frozen skin (999) as the closed outer boundary and the
fault (101) as the internal crack -- and hand tetgen just those facets and their
vertices. `-pq1.4/10Y`: -Y freezes every input facet (fault and skin survive
verbatim, so the stitch contract holds), -q fills the interior with fresh
well-shaped Steiner points. Old interior vertices are NOT passed in -- the whole
point is that tetgen owns the interior.

Verified afterwards, not assumed: every input facet must survive as a tet face,
and the refilled volume must match (checked again by stitch_patch at splice time).
"""
import argparse
import sys
import time

import meshio
import numpy as np
import tetgen
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument('--box', required=True, help='.msh from extract_patch_medit')
ap.add_argument('--out', required=True, help='refilled patch .msh (same tags)')
ap.add_argument('--shuffle', type=int, default=0)
ap.add_argument('--minratio', type=float, default=1.4)
ap.add_argument('--mindihedral', type=float, default=10.0)
a = ap.parse_args()
t0 = time.time()

m = meshio.read(a.box)
P = m.points
tri = np.vstack([c.data for c in m.cells if c.type == 'triangle']).astype(np.int64)
ref = np.concatenate([d for c, d in zip(m.cells, m.cell_data['gmsh:physical'])
                      if c.type == 'triangle']).astype(np.int32)
oldT = np.vstack([c.data for c in m.cells if c.type == 'tetra']).astype(np.int64)
print(f'[box] {len(oldT):,} old tets, {len(tri):,} facets '
      f'(101:{int((ref==101).sum()):,} 999:{int((ref==999).sum()):,})')

# facet-only point set
vid = np.unique(tri)
rm = -np.ones(len(P), np.int64)
rm[vid] = np.arange(len(vid))
FP = P[vid]
FT = rm[tri]

if a.shuffle:
    rng = np.random.default_rng(a.shuffle)
    vp = rng.permutation(len(FP))
    inv = np.empty_like(vp); inv[vp] = np.arange(len(vp))
    FP = FP[vp]; FT = inv[FT]
    fp = rng.permutation(len(FT)); FT = FT[fp]; ref = ref[fp]

shift = FP.mean(0)
sw = f'pq{a.minratio}/{a.mindihedral}Y'
print(f'[tetgen] -{sw}', flush=True)
tg = tetgen.TetGen(np.ascontiguousarray(FP - shift), np.ascontiguousarray(FT))
tg.tetrahedralize(switches=sw)
TP = np.asarray(tg.node) + shift
TT = np.asarray(tg.elem)
print(f'[tetgen] {len(TT):,} tets, {len(TP):,} verts ({time.time()-t0:.1f} s)')

# every input facet must survive as a tet face (that is what -Y promises; verify)
kd = cKDTree(TP)
d, mp = kd.query(FP, k=1)
if d.max() > 1e-6:
    raise SystemExit(f'input vertices moved: max {d.max():.3e} m')
faces = np.sort(np.concatenate([TT[:, [0, 1, 2]], TT[:, [0, 1, 3]],
                                TT[:, [0, 2, 3]], TT[:, [1, 2, 3]]]), axis=1)
fset = set(map(tuple, np.unique(faces, axis=0).tolist()))
want = np.sort(mp[FT], axis=1)
lost = sum(1 for t in map(tuple, want.tolist()) if t not in fset)
print(f'[check] facets lost: {lost} of {len(FT):,}')
if lost:
    raise SystemExit('facet loss -- retry with a different --shuffle seed')

# eta of the refill
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
p = TP[TT]
e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                      p[:, 3] - p[:, 0])) / 6.0
eta = np.where((e ** 2).sum(1) > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / (e ** 2).sum(1), 0.0)
print(f'[eta] min {eta.min():.4f}   <0.05 {int((eta<0.05).sum()):,}   '
      f'<0.1 {int((eta<0.1).sum()):,}   med {np.median(eta):.3f}   edge_min {e.min():.2f}')

out = meshio.Mesh(TP, [('triangle', mp[FT]), ('tetra', TT)],
                  cell_data={'gmsh:physical': [ref, np.ones(len(TT), np.int32)],
                             'gmsh:geometrical': [ref, np.ones(len(TT), np.int32)]})
meshio.write(a.out, out, file_format='gmsh22', binary=False)
print(f'[out] {a.out}')
