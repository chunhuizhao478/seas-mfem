#!/usr/bin/env python3
"""fill_bisect.py -- WHERE does `-pY` boundary recovery fail on the ALT PLC?

The full fill dies with `Internal TetGen error within recoversubfaces` during
"Recovering boundaries...". tetgen writes no skipped-face dump for an internal
error, so the offending region has to be found by bisection.

Method: keep the hull WHOLE (it fills on its own) and admit only a SUBSET of the
fault, selected by depth and/or along-strike position. If a run succeeds, the
excluded part contains the defect.

Dropping shallow fault triangles is safe for this purpose: the trace edges stay
embedded in the top surface as constrained edges with no fault attached, which
tetgen handles, so a pass isolates the fault side rather than changing the hull.

Coordinates: `s` is along-strike distance from the W corner in the rotated frame
the domain is built in (u = (S-W)/|S-W|), matching check_fault_in_domain.py.

  --zmax Z    keep fault triangles whose SHALLOWEST vertex is below Z (Z<=0)
  --smin/--smax  keep only triangles whose centroid s lies in the window
  --invert    keep the COMPLEMENT of the s-window (to confirm a hit)
"""
import argparse
import time

import numpy as np
import tetgen
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--zmax", type=float, default=0.0)
ap.add_argument("--smin", type=float, default=-1e18)
ap.add_argument("--smax", type=float, default=1e18)
ap.add_argument("--invert", action="store_true")
ap.add_argument("--shuffle", type=int, default=0,
                help="permute vertex/facet order with this seed (0 = keep input order)")
ap.add_argument("--switches", default="pY",
                help="tetgen switches; 'p' (no Y) lets it add Steiner points on facets, "
                     "which recovers boundaries it otherwise cannot -- at the cost of "
                     "altering the frozen fault triangulation, so the damage is counted")
a = ap.parse_args()

H = np.load("build_tmp/hull.npz")
F = np.load("build_tmp/fault_surface.npz")
HP, HT, PART = H["P"], H["T"], H["PART"]
FP, FT = F["P"], F["T"]

W, N, E, S = np.load("build_tmp/domain_corners_utm.npy")
u = (S - W) / np.linalg.norm(S - W)

zmaxtri = FP[FT][:, :, 2].max(1)
cen = FP[FT].mean(1)
s = (cen[:, :2] - W) @ u
inwin = (s >= a.smin) & (s <= a.smax)
if a.invert:
    inwin = ~inwin
# The s-window is applied ONLY to the shallow band. Measured: admitting every
# triangle below -1 m passes, so the deep fault is innocent and must stay WHOLE --
# otherwise each window would also carve the deep sheet and the bisection would be
# testing two changes at once.
deep = zmaxtri <= a.zmax
keep = deep | (~deep & inwin)
FTk = FT[keep]
print(f"[select] zmax<={a.zmax:,.0f}  s in [{a.smin:,.0f},{a.smax:,.0f}]"
      f"{' INVERTED' if a.invert else ''}  ->  {len(FTk):,} of {len(FT):,} fault triangles")

tree = cKDTree(HP)
d, k = tree.query(FP, k=1, distance_upper_bound=1e-3)
hit = np.isfinite(d)
newid = np.empty(len(FP), np.int64)
newid[hit] = k[hit]
newid[~hit] = len(HP) + np.arange(int((~hit).sum()))
P = np.vstack([HP, FP[~hit]])
T = np.vstack([HT, newid[FTk]])
print(f"[plc] {len(T):,} facets, {len(P):,} vertices")

# tetgen's boundary recovery is ORDER-DEPENDENT: vertices are inserted in input
# order, so the Delaunay history -- and therefore which subfaces are recoverable by
# flips alone -- changes with the numbering.  Relabelling is an exact operation
# (identical point set, identical facets), so if a permutation recovers under -Y it
# is a free fix: nothing about the geometry or the frozen fault has changed.
if a.shuffle:
    rng = np.random.default_rng(a.shuffle)
    perm = rng.permutation(len(P))
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    P = P[perm]
    T = inv[T]
    newid = inv[newid]
    T = T[rng.permutation(len(T))]
    print(f"[shuffle] vertex + facet order permuted with seed {a.shuffle}")

t0 = time.time()
shift = P.mean(0)
tg = tetgen.TetGen(np.ascontiguousarray(P - shift), np.ascontiguousarray(T))
try:
    tg.tetrahedralize(switches=a.switches)
    TP = np.asarray(tg.node) + shift
    TT = np.asarray(tg.elem)
    print(f"[RESULT] PASS -- {len(TT):,} tets, {len(TP):,} verts ({time.time() - t0:.0f} s)")
    print(f"   Steiner points added: {len(TP) - len(P):,}")
    # How much of the FROZEN fault survived?  A surviving facet is one whose vertex
    # triple is still a face of some tet.  Counted, never assumed -- -Y is what
    # normally guarantees this, so without it the number is the whole question.
    kd = cKDTree(TP)
    dd, mp = kd.query(P, k=1)
    print(f"   original PLC vertices preserved to {dd.max():.3e} m")
    faces = np.sort(np.concatenate([TT[:, [0, 1, 2]], TT[:, [0, 1, 3]],
                                    TT[:, [0, 2, 3]], TT[:, [1, 2, 3]]]), axis=1)
    fset = set(map(tuple, np.unique(faces, axis=0).tolist()))
    ftri = np.sort(mp[newid[FTk]], axis=1)
    surv = sum(1 for t_ in map(tuple, ftri.tolist()) if t_ in fset)
    print(f"   fault facets surviving as tet faces: {surv:,} of {len(FTk):,} "
          f"({100.0 * surv / len(FTk):.4f} %)  -> {len(FTk) - surv:,} LOST")
except RuntimeError as e:
    print(f"[RESULT] FAIL -- {str(e).strip()[:80]} ({time.time() - t0:.0f} s)")
