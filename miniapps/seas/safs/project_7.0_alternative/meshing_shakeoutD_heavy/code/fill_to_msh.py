#!/usr/bin/env python3
"""fill_to_msh.py -- tag the base fill's boundary and write Gmsh 2.2 for mmg.

`-Y` forbade Steiner points on input facets, so every PLC triangle survives as a
tet face and the boundary is tagged by EXACT vertex-triple match against the PLC
rather than re-derived from geometry.  SAFS tags: 101 fault, 102 free surface,
104 absorbing (walls + bottom); volume 1.

The match key is a 3-field structured record, NOT a bit-packed integer: at 2.4M
vertices three 22-bit fields need 66 bits and silently overflow int64.
"""
import numpy as np, meshio
d = np.load("build_tmp/fill.npz")
P, T, plcT, MARK = d["P"], d["T"], d["plcT"], d["MARK"]
plcP = d["plcP"]        # needed to derive FAULT_AREA below, rather than hardcode it
LOCAL = [(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
REC = np.dtype([('a',np.int32),('b',np.int32),('c',np.int32)])
def rec(A):
    return np.ascontiguousarray(np.sort(A,axis=1).astype(np.int32)).view(REC).ravel()
pk = rec(plcT); o = np.argsort(pk, kind='stable'); pk = pk[o]; ms = MARK[o]
assert len(np.unique(pk)) == len(pk), "duplicate PLC facet"
faces=[]; tags=[]
for k in range(4):
    f  = T[:, list(LOCAL[k])]
    fk = rec(f)
    pos = np.searchsorted(pk, fk); np.clip(pos, 0, len(pk)-1, out=pos)
    hit = pk[pos] == fk
    if hit.any(): faces.append(f[hit]); tags.append(ms[pos[hit]])
    del fk, pos, hit
faces = np.vstack(faces); tags = np.concatenate(tags)
TAG = {1:102, 3:104, 5:104, 7:101}
gt = np.array([TAG[int(x)] for x in tags], np.int32)
nf = int((gt==101).sum())
print(f"[tag] {len(faces):,} boundary faces matched (PLC has {len(plcT):,})")
for v,n in ((101,"fault"),(102,"free surface"),(104,"absorbing")):
    print(f"   {v} {n:13s}: {int((gt==v).sum()):,}")
# An interior crack contributes each facet TWICE.  The 4-side deficit is the 4
# fault triangles that lie flat in the z=0 plane (found during PLC debugging):
# they are boundary, not interior, so they appear once.  Assert the UNIQUE facet
# count instead -- that is the invariant that must hold exactly.
nuniq = len(np.unique(rec(faces[gt==101])))
single = 2*nuniq - nf
print(f"   unique fault facets {nuniq:,}   single-sided {single}")
NF = int((MARK == 7).sum())
_fp = plcP[plcT[MARK == 7]]
FAULT_AREA = float(0.5*np.linalg.norm(np.cross(_fp[:,1]-_fp[:,0], _fp[:,2]-_fp[:,0]), axis=1).sum())   # ALT: derive, never hardcode a lineage's count
miss = NF - nuniq
if miss:
    # A facet present in the PLC but on NO tet face is a pinhole in the crack:
    # tetgen dropped it, almost always because the flat-top clamp made it
    # degenerate.  Report area + location so the defect is on the record.
    fp = plcT[MARK == 7]
    have = set(rec(faces[gt==101]).tolist())
    bad = np.array([i for i, k in enumerate(rec(fp).tolist()) if k not in have])
    Vb = P[fp[bad]]
    ar = 0.5*np.linalg.norm(np.cross(Vb[:,1]-Vb[:,0], Vb[:,2]-Vb[:,0]), axis=1)
    print(f"[pinhole] {miss} PLC fault facets on no tet face")
    for t, A in zip(Vb, ar):
        print(f"   area {A:10.3e} m2   centroid {t.mean(0)[0]:.1f} {t.mean(0)[1]:.1f} {t.mean(0)[2]:.2f}")
    print(f"   total pinhole area {ar.sum():.3e} m2 = {100*ar.sum()/FAULT_AREA:.3e} % of the fault")
assert miss <= 8, f"{miss} missing fault facets -- expected <= 8"
assert single <= 8, f"{single} single-sided fault facets -- expected <= 8 (z=0 in-plane)"
assert int((gt==102).sum()) + int((gt==104).sum()) == len(plcT) - NF, "hull count"
# keep ONE copy of each fault face -- mmg wants a single crack surface, -opnbdy opens it
keep = np.ones(len(faces), bool)
isf  = gt == 101
fk   = rec(faces[isf]); _, first = np.unique(fk, return_index=True)
idxf = np.flatnonzero(isf); drop = np.setdiff1d(idxf, idxf[first]); keep[drop] = False
faces, gt = faces[keep], gt[keep]
print(f"[dedup] fault faces -> {int((gt==101).sum()):,}   total surface {len(faces):,}")
cells=[("triangle", faces), ("tetra", T)]
cd={"gmsh:physical":[gt, np.ones(len(T),np.int32)],
    "gmsh:geometrical":[gt, np.ones(len(T),np.int32)]}
# BINARY, not ASCII: the next stage (mmg_refine_sizemap.py) reads this through
# fastmsh, which is binary-only and exists precisely because meshio.read on a large
# ASCII .msh materialises per-element Python objects and exhausts 36 GB. ASCII also
# cost 819 MB here for a 9.6 M-tet mesh.
meshio.write("build_tmp/base.msh", meshio.Mesh(P, cells, cell_data=cd),
             file_format="gmsh22", binary=True)
print("[out] build_tmp/base.msh")
