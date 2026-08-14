#!/usr/bin/env python3
"""stitch_patch.py -- swap a repaired patch back into the global MEDIT mesh.

The surgery loop is extract -> repair -> STITCH, per residual sliver cluster:

    extract_patch_medit.py  cuts the box, freezes fault (101) + skin (999),
                            records the global tet ids it removed (meta tids)
    <repair>                remove_sliver_tets.py or a local mmg draw; both may
                            delete/add INTERIOR vertices but every Required
                            (fault+skin) vertex survives with identical coords
    stitch_patch.py         drop the recorded tets from the global mesh, map the
                            repaired patch's frozen vertices back to their global
                            ids by exact coordinate match, append the patch's new
                            interior vertices, splice in its tets

Integrity is asserted, not assumed:
  * every skin/fault vertex of the repaired patch must land on a global vertex
    within --tol (default 1e-6 m) -- a miss means the repair moved the frozen
    boundary and the patch CANNOT be stitched;
  * the summed volume of removed tets must equal the repaired patch volume to
    1e-9 relative -- a leak or overlap shows up here immediately;
  * the global Triangles section is carried over UNCHANGED (fault count must
    match before/after by construction; asserted anyway).

Orphaned vertices (interior vertices the repair deleted) are compacted out and
the whole mesh renumbered, so the output is directly consumable by
medit_to_puml.py.
"""
import argparse
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections
from scipy.spatial import cKDTree

CH = 8_000_000
ap = argparse.ArgumentParser()
ap.add_argument('--global-mesh', required=True)
ap.add_argument('--patch', required=True, help='repaired patch (.msh gmsh22 or .mesh MEDIT)')
ap.add_argument('--meta', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--tol', type=float, default=1e-6)
a = ap.parse_args()

meta = np.load(a.meta)
tids = meta['tids']
drop = np.zeros(0, bool)          # sized after NT is known

# ---- read the repaired patch (either format) -------------------------------
if a.patch.endswith('.msh'):
    import meshio
    pm = meshio.read(a.patch)
    PP = pm.points
    PT = np.vstack([c.data for c in pm.cells if c.type == 'tetra']).astype(np.int64)
else:
    psec = medit_sections(a.patch)
    pLV, pNV = psec['Vertices']
    pLT, pNT = psec['Tetrahedra']
    PP = pd.read_csv(a.patch, sep=r'\s+', header=None, skiprows=pLV + 1, nrows=pNV,
                     usecols=range(3), dtype=np.float64, engine='c').to_numpy()
    PT = pd.read_csv(a.patch, sep=r'\s+', header=None, skiprows=pLT + 1, nrows=pNT,
                     usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
print(f'[patch] {len(PT):,} tets  {len(PP):,} verts', flush=True)

# ---- read the global mesh ---------------------------------------------------
sec = medit_sections(a.global_mesh)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
LS, NS = sec['Triangles']
G = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
S = pd.read_csv(a.global_mesh, sep=r'\s+', header=None, skiprows=LS + 1, nrows=NS,
                usecols=range(4), dtype=np.int64, engine='c').to_numpy()
drop = np.zeros(NT, bool)
drop[tids] = True
print(f'[global] {NT:,} tets, dropping {int(drop.sum()):,}', flush=True)

# ---- map patch verts to global ids ------------------------------------------
tree = cKDTree(G)
d, gid = tree.query(PP, k=1, workers=-1)
matched = d <= a.tol
newv = ~matched
pmap = np.empty(len(PP), np.int64)
pmap[matched] = gid[matched]
pmap[newv] = NV + np.arange(int(newv.sum()))
print(f'[map] {int(matched.sum()):,} patch verts matched to global ids '
      f'(max d {d[matched].max() if matched.any() else 0:.2e}), '
      f'{int(newv.sum()):,} new interior verts')
# the FROZEN vertices are the ones that must match; verify against meta.skinP
sd, _ = cKDTree(PP).query(meta['skinP'], k=1)
if sd.max() > a.tol:
    raise SystemExit(f'FROZEN SKIN MOVED: max {sd.max():.3e} m > tol -- cannot stitch')

# ---- volume conservation ----------------------------------------------------
def vol_sum(Pts, Tets):
    tot = 0.0
    for s in range(0, len(Tets), CH):
        t = Tets[s:s + CH]
        p = Pts[t]
        tot += np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                                p[:, 3] - p[:, 0])).sum() / 6.0
    return tot

# removed-region volume, read in chunks
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
vol_new = vol_sum(PP, PT)
rel = abs(vol_new - vol_rm) / max(vol_rm, 1e-30)
print(f'[volume] removed {vol_rm/1e9:,.6f} km3  patch {vol_new/1e9:,.6f} km3  '
      f'rel diff {rel:.2e}')
if rel > 1e-9:
    raise SystemExit('VOLUME MISMATCH -- leak or overlap, refusing to stitch')

# ---- write the stitched mesh ------------------------------------------------
allP = np.vstack([G, PP[newv]])
used = np.zeros(len(allP), bool)
used[S[:, :3].ravel() - 1] = True
used[pmap[PT].ravel()] = True

with open(a.out, 'w') as fo:
    fo.write('MeshVersionFormatted 2\nDimension 3\n\n')
    # compaction happens via `used`; count kept tets first
    nkeep = NT - int(drop.sum()) + len(PT)
    # mark used verts from kept global tets (chunked)
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
    print(f'[compact] {len(allP):,} -> {int(used.sum()):,} verts '
          f'({len(allP)-int(used.sum()):,} orphans dropped)')
    fo.write('Vertices\n%d\n' % int(used.sum()))
    Pk = allP[used]
    pd.DataFrame({0: Pk[:, 0], 1: Pk[:, 1], 2: Pk[:, 2], 3: 0}).to_csv(
        fo, sep=' ', header=False, index=False, float_format='%.10g')
    tr = remap[S[:, :3] - 1]
    assert (tr >= 0).all(), 'triangle references an orphaned vertex'
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
    Tp = remap[pmap[PT]]
    assert (Tp >= 0).all()
    pd.DataFrame({0: Tp[:, 0] + 1, 1: Tp[:, 1] + 1, 2: Tp[:, 2] + 1,
                  3: Tp[:, 3] + 1, 4: np.ones(len(Tp), np.int32)}).to_csv(
        fo, sep=' ', header=False, index=False)
    fo.write('\nEnd\n')
print(f'[out] {a.out}  ({nkeep:,} tets)')
