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
np.savez("build_tmp/fill.npz", P=TP, T=TT, plcP=P, plcT=T, MARK=MARK)
print(f"[out] build_tmp/fill.npz  ({time.time()-t0:.0f} s)")
