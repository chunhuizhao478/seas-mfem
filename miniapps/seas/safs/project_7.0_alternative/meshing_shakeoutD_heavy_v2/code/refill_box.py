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
ap.add_argument('--polish', type=int, default=2,
                help='box-local mmg -optim passes after the refill. Raw tetgen -q '
                     'bounds radius-edge, NOT Joe-Liu eta: measured, unpolished '
                     'refills minted 51 new sub-0.05 cells (one at 0.0001) across '
                     '17 boxes. The minted cells are INTERIOR, which optim eats '
                     'in seconds at box scale. 0 disables.')
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
lostm = np.array([t not in fset for t in map(tuple, want.tolist())], bool)
lost = int(lostm.sum())
print(f'[check] facets lost: {lost} of {len(FT):,}')
if lost:
    # say WHICH facets, so a grow-the-box retry can be judged: a lost 999 (skin)
    # facet often lands interior to a larger box; a lost 101 (fault) facet is the
    # same recovery-residue class as the base fill's pinholes.
    for r in np.unique(ref[lostm]):
        print(f'  lost ref {r}: {int((ref[lostm] == r).sum())}')
    A = FP[FT[lostm]]
    ar = 0.5 * np.linalg.norm(np.cross(A[:, 1] - A[:, 0], A[:, 2] - A[:, 0]), axis=1)
    for i in range(min(4, lost)):
        print(f'  lost facet area {ar[i]:.1f} m2  centroid {A[i].mean(0)}')
    raise SystemExit('facet loss -- retry with a different --shuffle seed or a grown box')

# eta of the refill
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def eta_of(Pts, Tets):
    p = Pts[Tets]
    e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
    d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0])) / 6.0
    ss = (e ** 2).sum(1)
    return np.where(ss > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ss, 0.0), e


eta, e = eta_of(TP, TT)
print(f'[eta] refill min {eta.min():.4f}   <0.05 {int((eta<0.05).sum()):,}   '
      f'<0.1 {int((eta<0.1).sum()):,}   med {np.median(eta):.3f}   edge_min {e.min():.2f}')

# ---- box-local optim polish -------------------------------------------------
if a.polish and int((eta < 0.1).sum()):
    import os
    import subprocess
    import pandas as pd
    sys.path.insert(0, 'code')
    from medit_hdr import medit_sections
    MMG = '/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3'
    tmp = a.out + '.polish.mesh'
    tmpo = a.out + '.polish_o.mesh'
    for it in range(a.polish):
        if int((eta < 0.1).sum()) == 0:
            break
        for f_ in (tmp, tmpo, tmp[:-5] + '.sol', tmpo[:-5] + '.sol'):
            try:
                os.unlink(f_)
            except OSError:
                pass
        with open(tmp, 'w') as fo:
            fo.write('MeshVersionFormatted 2\nDimension 3\n\nVertices\n%d\n' % len(TP))
            pd.DataFrame({0: TP[:, 0], 1: TP[:, 1], 2: TP[:, 2], 3: 0}).to_csv(
                fo, sep=' ', header=False, index=False, float_format='%.17g')
            tl = mp[FT]
            fo.write('\nTriangles\n%d\n' % len(tl))
            pd.DataFrame({0: tl[:, 0] + 1, 1: tl[:, 1] + 1, 2: tl[:, 2] + 1,
                          3: ref}).to_csv(fo, sep=' ', header=False, index=False)
            fo.write('\nRequiredTriangles\n%d\n' % len(tl))
            pd.DataFrame({0: np.arange(1, len(tl) + 1)}).to_csv(
                fo, sep=' ', header=False, index=False)
            fo.write('\nTetrahedra\n%d\n' % len(TT))
            pd.DataFrame({0: TT[:, 0] + 1, 1: TT[:, 1] + 1, 2: TT[:, 2] + 1,
                          3: TT[:, 3] + 1, 4: np.ones(len(TT), np.int32)}).to_csv(
                fo, sep=' ', header=False, index=False)
            fo.write('\nEnd\n')
        r = subprocess.run([MMG, '-in', tmp, '-out', tmpo, '-opnbdy', '-optim',
                            '-nosurf', '-hmin', '20', '-hmax', '200000',
                            '-hgrad', '3', '-m', '3000', '-v', '0'],
                           capture_output=True)
        if r.returncode != 0 or not os.path.exists(tmpo):
            print(f'[polish {it}] mmg failed (rc {r.returncode}) -- keeping pre-polish state')
            break
        sec = medit_sections(tmpo)
        pLV, pNV = sec['Vertices']
        pLT, pNT = sec['Tetrahedra']
        TPn = pd.read_csv(tmpo, sep=r'\s+', header=None, skiprows=pLV + 1, nrows=pNV,
                          usecols=range(3), dtype=np.float64, engine='c').to_numpy()
        TTn = pd.read_csv(tmpo, sep=r'\s+', header=None, skiprows=pLT + 1, nrows=pNT,
                          usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
        kd2 = cKDTree(TPn)
        d2, mp2 = kd2.query(FP, k=1)
        if d2.max() > 1e-6:
            print(f'[polish {it}] moved input verts ({d2.max():.2e}) -- keeping pre-polish state')
            break
        faces2 = np.sort(np.concatenate([TTn[:, [0, 1, 2]], TTn[:, [0, 1, 3]],
                                         TTn[:, [0, 2, 3]], TTn[:, [1, 2, 3]]]), axis=1)
        fset2 = set(map(tuple, np.unique(faces2, axis=0).tolist()))
        want2 = np.sort(mp2[FT], axis=1)
        lost2 = sum(1 for t in map(tuple, want2.tolist()) if t not in fset2)
        if lost2:
            print(f'[polish {it}] lost {lost2} facets -- keeping pre-polish state')
            break
        TP, TT, mp = TPn, TTn, mp2
        eta, e = eta_of(TP, TT)
        print(f'[polish {it}] min {eta.min():.4f}   <0.05 {int((eta<0.05).sum()):,}   '
              f'<0.1 {int((eta<0.1).sum()):,}   tets {len(TT):,}')
    for f_ in (tmp, tmpo, tmp[:-5] + '.sol', tmpo[:-5] + '.sol'):
        try:
            os.unlink(f_)
        except OSError:
            pass

out = meshio.Mesh(TP, [('triangle', mp[FT]), ('tetra', TT)],
                  cell_data={'gmsh:physical': [ref, np.ones(len(TT), np.int32)],
                             'gmsh:geometrical': [ref, np.ones(len(TT), np.int32)]})
meshio.write(a.out, out, file_format='gmsh22', binary=False)
print(f'[out] {a.out}')
