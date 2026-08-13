#!/usr/bin/env python3
"""localbox_fill.py -- can tetgen recover THIS fault geometry under the real `-pY`?

The full ShakeOut-D fill is the one stage that must never overlap another heavy job
on this box, so it is a bad way to answer "is the fault triangulation recoverable?".
This runs the *identical* switches on a small closed box around a chosen piece of
the fault -- same code path, ~1 % of the memory, under a second.

Why it is worth having: `isolate_selfint.py` (bare `-d`) died with an *internal*
TetGen error inside `recoversubfaces` on the assembled ALT PLC, while every
structural audit of that PLC came back clean (exact weld, no coincident vertices,
no T-junctions, all 8,137 trace edges embedded).  `-d` is a detection-only path and
is not what the build actually runs.  Boxing the worst geometry and running `-pY`
separates "the geometry is unbuildable" from "the diagnostic crashed", and the
answer was the latter -- see docs/PLC_GATE_DIAGNOSIS.md.

The fault patch inside the box is an INTERNAL OPEN CRACK (it terminates in mid-air
at the box walls), which is strictly harder for facet recovery than the real case
where the crack reaches the free surface, so a pass here is a conservative result.

  usage: localbox_fill.py HALF_WIDTH_M Z_LO Z_HI [CX CY]

Interior Steiner points are expected and fine: `-Y` forbids them on input facets
only, so a pass means the frozen fault triangulation survived vertex-for-vertex.
"""
import sys
import time

import numpy as np
import tetgen

F = np.load("build_tmp/fault_surface.npz")
P, T = F["P"], F["T"]

H = float(sys.argv[1])
Z0, Z1 = float(sys.argv[2]), float(sys.argv[3])
cx = float(sys.argv[4]) if len(sys.argv) > 4 else 573455.0
cy = float(sys.argv[5]) if len(sys.argv) > 5 else 3735271.0

lo = np.array([cx - H, cy - H, Z0])
hi = np.array([cx + H, cy + H, Z1])

# 1 m inset: a fault vertex sitting exactly ON a box wall would be a degeneracy of
# this harness, not of the mesh under test.
ins = np.all((P > lo + 1.0) & (P < hi - 1.0), axis=1)
Tk = T[ins[T].all(1)]
if len(Tk) == 0:
    sys.exit("no fault triangles in box")

vid = np.unique(Tk)
rm = -np.ones(len(P), np.int64)
rm[vid] = np.arange(len(vid))
FP, FT = P[vid], rm[Tk]

a = np.linalg.norm(FP[FT[:, 1]] - FP[FT[:, 0]], axis=1)
b = np.linalg.norm(FP[FT[:, 2]] - FP[FT[:, 1]], axis=1)
c = np.linalg.norm(FP[FT[:, 0]] - FP[FT[:, 2]], axis=1)
ar = 0.5 * np.linalg.norm(np.cross(FP[FT[:, 1]] - FP[FT[:, 0]],
                                   FP[FT[:, 2]] - FP[FT[:, 0]]), axis=1)
q = 4 * np.sqrt(3) * ar / (a * a + b * b + c * c)
e = np.sort(np.concatenate([FT[:, [0, 1]], FT[:, [1, 2]], FT[:, [0, 2]]]), axis=1)
_, cnt = np.unique(e, axis=0, return_counts=True)
print(f"[patch] {len(FT):,} fault tris, {len(FP):,} verts, "
      f"min edge {min(a.min(), b.min(), c.min()):.3f} m, min q {q.min():.5f}, "
      f"needles(q<0.1) {int((q < 0.1).sum())}, non-manifold edges {int((cnt > 2).sum())}")

X, Y, Zb = [lo[0], hi[0]], [lo[1], hi[1]], [lo[2], hi[2]]
BP = np.array([[X[i], Y[j], Zb[k]] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float)


def idx(i, j, k):
    return i * 4 + j * 2 + k


BT = []
for k in (0, 1):
    BT += [[idx(0, 0, k), idx(1, 0, k), idx(1, 1, k)],
           [idx(0, 0, k), idx(1, 1, k), idx(0, 1, k)]]
for j in (0, 1):
    BT += [[idx(0, j, 0), idx(1, j, 0), idx(1, j, 1)],
           [idx(0, j, 0), idx(1, j, 1), idx(0, j, 1)]]
for i in (0, 1):
    BT += [[idx(i, 0, 0), idx(i, 1, 0), idx(i, 1, 1)],
           [idx(i, 0, 0), idx(i, 1, 1), idx(i, 0, 1)]]

AP = np.vstack([BP, np.asarray(FP)])
AT = np.vstack([np.array(BT, np.int64), FT + len(BP)])
print(f"[plc] {len(AT):,} facets, {len(AP):,} vertices")

t0 = time.time()
shift = AP.mean(0)
tg = tetgen.TetGen(np.ascontiguousarray(AP - shift), np.ascontiguousarray(AT))
tg.tetrahedralize(switches="pY")
TT, TP = np.asarray(tg.elem), np.asarray(tg.node) + shift
print(f"[RESULT] -pY SUCCEEDED: {len(TT):,} tets, {len(TP):,} verts ({time.time() - t0:.1f} s)")
print(f"   interior Steiner points added: {len(TP) - len(AP)}")
