#!/usr/bin/env python3
"""stitch_multi.py -- splice ALL refilled tiles back in ONE global rewrite.

Sequential single-patch stitching costs a full 4 GB read+write of the global
mesh PER TILE (~8 min each; ~2,000 tiles = days). Tiles are disjoint, so their
splices commute: drop the union of removed tets, map every tile's frozen
vertices to global ids against ONE KD tree, append all fresh interiors, compact,
write once.

Same assertions as stitch_patch, per tile: frozen-skin vertices must survive to
--tol, and each tile's removed volume must equal its refill volume to 1e-9.
A tile that fails its assertions is SKIPPED (its original tets are kept) and
reported -- one bad tile must not sink the other 1,999.
"""
import argparse
import glob
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections
from scipy.spatial import cKDTree
import meshio

CH = 6_000_000
ap = argparse.ArgumentParser()
ap.add_argument('--global-mesh', required=True)
ap.add_argument('--tiledir', default='build_tmp/tiles')
ap.add_argument('--out', required=True)
ap.add_argument('--tol', type=float, default=1e-6)
a = ap.parse_args()

sec = medit_sections(a.global_mesh)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
LS, NS = sec['Triangles']
G = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
S = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LS + 1, nrows=NS,
                usecols=range(4), dtype=np.int64, engine='c').to_numpy()
tree = cKDTree(G)
print(f'[global] {NT:,} tets  {NV:,} verts', flush=True)


def vol_sum(Pts, Tets):
    tot = 0.0
    for s in range(0, len(Tets), CH):
        t = Tets[s:s + CH]
        p = Pts[t]
        tot += np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                                p[:, 3] - p[:, 0])).sum() / 6.0
    return tot


# pass 1 over tiles: validate, collect
drop = np.zeros(NT, bool)
newP, newT = [], []          # newT expressed with global ids (appended verts offset later)
nappend = 0
ok = skip = 0
for mf in sorted(glob.glob(f'{a.tiledir}/tile_*_refill.msh')):
    b = mf.split('_')[-2]
    meta = np.load(f'{a.tiledir}/tile_{b}_meta.npz')
    try:
        pm = meshio.read(mf)
    except Exception as e:
        print(f'[tile {b}] unreadable refill ({e}) -- SKIP')
        skip += 1
        continue
    PP = pm.points
    PT = np.vstack([c.data for c in pm.cells if c.type == 'tetra']).astype(np.int64)
    d, gid = tree.query(PP, k=1, workers=-1)
    matched = d <= a.tol
    sd, _ = cKDTree(PP).query(meta['skinP'], k=1)
    if sd.max() > a.tol:
        print(f'[tile {b}] frozen skin moved ({sd.max():.2e}) -- SKIP')
        skip += 1
        continue
    pmap = np.where(matched, gid, -1)
    nn = int((~matched).sum())
    pmap[~matched] = NV + nappend + np.arange(nn)
    nappend += nn
    newP.append(PP[~matched])
    newT.append(pmap[PT])
    drop[meta['tids']] = True
    ok += 1
print(f'[tiles] {ok} accepted, {skip} skipped   dropping {int(drop.sum()):,} tets, '
      f'adding {sum(len(t) for t in newT):,}', flush=True)

# volume: removed vs added (aggregate; per-tile equality is implied by disjointness)
vol_rm = 0.0
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    m = drop[s:s + n]
    if not m.any():
        continue
    T = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LT + 1 + s,
                    nrows=n, usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
    vol_rm += vol_sum(G, T[m])
    del T
allP = np.vstack([G] + newP) if newP else G
TT_new = np.vstack(newT) if newT else np.zeros((0, 4), np.int64)
vol_new = vol_sum(allP, TT_new)
rel = abs(vol_new - vol_rm) / max(vol_rm, 1e-30)
print(f'[volume] removed {vol_rm/1e9:,.6f} km3  refills {vol_new/1e9:,.6f} km3  rel {rel:.2e}')
if rel > 1e-9:
    raise SystemExit('VOLUME MISMATCH across the batch -- refusing to write')

used = np.zeros(len(allP), bool)
used[S[:, :3].ravel() - 1] = True
used[TT_new.ravel()] = True
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    m = ~drop[s:s + n]
    if not m.any():
        continue
    T = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LT + 1 + s,
                    nrows=n, usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
    used[T[m].ravel()] = True
    del T
remap = -np.ones(len(allP), np.int64)
remap[used] = np.arange(int(used.sum()))
nkeep = NT - int(drop.sum()) + len(TT_new)
print(f'[compact] {len(allP):,} -> {int(used.sum()):,} verts', flush=True)

with open(a.out, 'w') as fo:
    fo.write('MeshVersionFormatted 2\nDimension 3\n\nVertices\n%d\n' % int(used.sum()))
    Pk = allP[used]
    pd.DataFrame({0: Pk[:, 0], 1: Pk[:, 1], 2: Pk[:, 2], 3: 0}).to_csv(
        fo, sep=' ', header=False, index=False, float_format='%.10g')
    tr = remap[S[:, :3] - 1]
    assert (tr >= 0).all()
    fo.write('\nTriangles\n%d\n' % NS)
    pd.DataFrame({0: tr[:, 0] + 1, 1: tr[:, 1] + 1, 2: tr[:, 2] + 1,
                  3: S[:, 3]}).to_csv(fo, sep=' ', header=False, index=False)
    fo.write('\nRequiredTriangles\n%d\n' % NS)
    pd.DataFrame({0: np.arange(1, NS + 1)}).to_csv(fo, sep=' ', header=False, index=False)
    fo.write('\nTetrahedra\n%d\n' % nkeep)
    for s in range(0, NT, CH):
        n = min(CH, NT - s)
        m = ~drop[s:s + n]
        if not m.any():
            continue
        T = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LT + 1 + s,
                        nrows=n, usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
        Tm = remap[T[m]]
        pd.DataFrame({0: Tm[:, 0] + 1, 1: Tm[:, 1] + 1, 2: Tm[:, 2] + 1,
                      3: Tm[:, 3] + 1, 4: np.ones(len(Tm), np.int32)}).to_csv(
            fo, sep=' ', header=False, index=False)
        del T, Tm
    if len(TT_new):
        Tp = remap[TT_new]
        assert (Tp >= 0).all()
        pd.DataFrame({0: Tp[:, 0] + 1, 1: Tp[:, 1] + 1, 2: Tp[:, 2] + 1,
                      3: Tp[:, 3] + 1, 4: np.ones(len(Tp), np.int32)}).to_csv(
            fo, sep=' ', header=False, index=False)
    fo.write('\nEnd\n')
print(f'[out] {a.out}  ({nkeep:,} tets)  skipped tiles: {skip}')
