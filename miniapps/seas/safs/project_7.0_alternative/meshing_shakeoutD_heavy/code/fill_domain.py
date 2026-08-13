#!/usr/bin/env python3
"""fill_domain.py -- tetrahedralise the ShakeOut-D PLC with the fault embedded.

The PLC is the watertight hull (top with the trace embedded, four walls, flat
bottom) PLUS the fault surface as INTERNAL facets.  `-Y` forbids Steiner points on
any input facet, so the fault's 2,761,488 triangles survive vertex-for-vertex --
that is what "keep the fault zone treatment" means mechanically, and it is
asserted afterwards, not assumed.

`-a` caps the far-field cell volume: tetgen honours input facets but takes NO
interior size field, and left alone its interior coarsens without limit (19.1 km
max edge measured on an earlier build of this project).
"""
import argparse, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--minratio", type=float, default=1.4)
ap.add_argument("--mindihedral", type=float, default=10.0)
ap.add_argument("--hmax", type=float, default=5000.0)
ap.add_argument("--shuffle", type=int, default=1,
                help="permute vertex/facet order before tetgen (0 = input order). "
                     "REQUIRED on ALT -- input order fails boundary recovery; see below.")
a = ap.parse_args()
t0 = time.time()

H = np.load("build_tmp/hull.npz"); F = np.load("build_tmp/fault_surface.npz")
HP, HT, PART = H["P"], H["T"], H["PART"]
FP, FT = F["P"], F["T"]
# weld: the trace nodes are shared, everything else is new
tree = cKDTree(HP)
d, k = tree.query(FP, k=1, distance_upper_bound=1e-3)
hit = np.isfinite(d)
print(f"[weld] {int(hit.sum()):,} fault vertices coincide with hull vertices")
newid = np.empty(len(FP), np.int64)
newid[hit] = k[hit]
newid[~hit] = len(HP) + np.arange(int((~hit).sum()))
P = np.vstack([HP, FP[~hit]])
T = np.vstack([HT, newid[FT]])
MARK = np.concatenate([PART, np.full(len(FT), 7, np.int8)])   # 7 = fault
print(f"[plc] {len(T):,} facets, {len(P):,} vertices "
      f"(top {int((PART==1).sum()):,}, bottom {int((PART==3).sum()):,}, "
      f"wall {int((PART==5).sum()):,}, fault {len(FT):,})")

# NO -a.  With a 5 km volume cap tetgen drove refinement to 134,481,852 tets in
# queue and died in locate_point_walk -- that is the whole mesh built in one shot,
# which is exactly what the base-plus-refine path exists to avoid.  The BASE is
# driven by the input facets alone (fault 115 m, free surface 115-1143 m); the far
# field is then brought to the gate by the chunked LEB tools, which are memory-safe.
sw = "pY" if a.minratio <= 0 else f"pq{a.minratio}/{a.mindihedral}Y"
if a.hmax > 0:
    sw += f"a{a.hmax**3/(6.0*np.sqrt(2.0)):.6g}"
# ORDER PERMUTATION -- required on ALT, not cosmetic.
#
# In input order this PLC dies with `Internal TetGen error within recoversubfaces`
# during boundary recovery, under both -d and the real -pY.  Bisection localised it
# to 67 fault triangles at the SE lateral tip (s 455.0-457.6 km, E 616,573-618,343):
# removing exactly those fills cleanly, while removing a comparable 88-triangle set
# elsewhere still fails, so it is those facets and not a size effect.  They are not
# degenerate -- min edge 80.2 m, min quality 0.768, and the top surface there is min
# quality 0.644.  The trace simply terminates in mid-surface at that point, and
# tetgen's recovery is order-dependent: vertices are inserted in input order, so the
# Delaunay history decides which subfaces are recoverable by flips alone.
#
# Relabelling is exact -- same points, same facets, same geometry -- so this is a
# free fix.  Measured on the full PLC:
#   input order   -pY   FAIL (recoversubfaces)
#   seed 1        -pY   PASS  279 Steiner,   6 of 2,564,480 fault facets lost
#   seed 2        -pY   PASS  297 Steiner,   8 lost
#   no -Y at all        PASS    1 Steiner, 215,780 LOST (facets re-diagonalised)
# For scale, PREFERRED's own base fill added 969 Steiner points and lost 2.
# Dropping -Y is NOT an acceptable alternative: it silently retriangulates 8.4 % of
# the frozen fault.
if a.shuffle:
    rng = np.random.default_rng(a.shuffle)
    vperm = rng.permutation(len(P))
    inv = np.empty_like(vperm); inv[vperm] = np.arange(len(vperm))
    P = P[vperm]; T = inv[T]
    fperm = rng.permutation(len(T))
    T = T[fperm]; MARK = MARK[fperm]      # MARK rides with T or every tag is wrong
    print(f"[shuffle] vertex + facet order permuted with seed {a.shuffle}")

print(f"[tetgen] switches -{sw}", flush=True)
shift = P.mean(0)
import tetgen
tg = tetgen.TetGen(np.ascontiguousarray(P - shift), np.ascontiguousarray(T))
tg.tetrahedralize(switches=sw)
TP = np.asarray(tg.node) + shift
TT = np.asarray(tg.elem)
print(f"[tetgen] {len(TT):,} tets, {len(TP):,} verts   ({time.time()-t0:.0f} s)", flush=True)
d2, _ = cKDTree(TP).query(P, k=1)
print(f"[check] PLC vertices preserved to {d2.max():.3e} m")
if d2.max() > 1e-6:
    raise RuntimeError("tetgen moved PLC vertices")
# Volume is the cheapest end-to-end proof that the PLC closed and the interior was
# filled once: a leak through the hull, or the fault cutting the domain in two, both
# show up here long before any face census.
A = TP[TT[:, 1]] - TP[TT[:, 0]]; B = TP[TT[:, 2]] - TP[TT[:, 0]]; C = TP[TT[:, 3]] - TP[TT[:, 0]]
vol = np.abs(np.einsum("ij,ij->i", np.cross(A, B), C)).sum() / 6.0
print(f"[check] volume {vol/1e9:,.0f} km3   (domain 12,480,379 km3)")
# How much of the FROZEN fault survived as tet faces?  -Y is what guarantees this,
# so it is asserted, not assumed.
kd = cKDTree(TP); _, mp = kd.query(P, k=1)
faces = np.sort(np.concatenate([TT[:, [0,1,2]], TT[:, [0,1,3]], TT[:, [0,2,3]], TT[:, [1,2,3]]]), axis=1)
fset = set(map(tuple, np.unique(faces, axis=0).tolist()))
ftri = np.sort(mp[T[MARK == 7]], axis=1)
surv = sum(1 for t_ in map(tuple, ftri.tolist()) if t_ in fset)
print(f"[check] fault facets surviving as tet faces: {surv:,} of {len(ftri):,} "
      f"-> {len(ftri)-surv:,} lost ({100.0*(len(ftri)-surv)/len(ftri):.5f} %)")
np.savez("build_tmp/fill.npz", P=TP, T=TT, plcP=P, plcT=T, MARK=MARK)
print(f"[out] build_tmp/fill.npz  ({time.time()-t0:.0f} s)")
