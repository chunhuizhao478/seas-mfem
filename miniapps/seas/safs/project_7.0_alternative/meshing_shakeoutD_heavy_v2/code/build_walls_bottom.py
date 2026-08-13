#!/usr/bin/env python3
"""build_walls_bottom.py -- the four absorbing walls and the flat bottom.

Every curve shared by two parts is generated ONCE and handed to both, with each
segment pinned `setTransfiniteCurve(...,2)`, so gmsh reproduces the node set
vertex-for-vertex rather than remeshing it.  That is what makes the assembled PLC
watertight, and it is asserted (every edge in exactly two triangles, Euler
V-E+F = 2 for a sphere-topology shell) rather than hoped for.

Order matters: TOP is already meshed, so the walls take their top edge from it;
the bottom then takes its boundary from the walls.
"""
import numpy as np, gmsh, os

Q = np.load("build_tmp/domain_corners_utm.npy")          # W,N,E,S
top = np.load("build_tmp/top_surface.npz"); TP, TT = top["P"], top["T"]
ZB = -80000.0
HMAX = 5000.0
G = 0.15

# ---- top boundary loop, split per side ------------------------------------
e = np.sort(np.concatenate([TT[:, [0, 1]], TT[:, [1, 2]], TT[:, [0, 2]]]), axis=1)
u_, c_ = np.unique(e, axis=0, return_counts=True)
bnd = u_[c_ == 1]
bv = np.unique(bnd)
print(f"[top] boundary {len(bnd):,} edges, {len(bv):,} nodes")

sides = []
for i in range(4):
    A, B = Q[i], Q[(i + 1) % 4]
    d = B - A; L = np.linalg.norm(d); t = d / L
    p = TP[bv, :2] - A
    s = p @ t; perp = np.abs(p @ np.array([-t[1], t[0]]))
    m = (perp < 1.0) & (s > -1.0) & (s < L + 1.0)
    idx = bv[m][np.argsort(s[m])]
    sides.append(idx)
    print(f"  side {i}: {len(idx):,} nodes over {L/1000:.1f} km "
          f"(spacing med {np.median(np.diff(np.sort(s[m]))):.0f} m)")
assert sum(len(s) for s in sides) - 4 == len(bv), \
    f"side split lost nodes: {sum(len(s) for s in sides)-4} vs {len(bv)}"

def zlevels():
    z = [0.0]
    while z[-1] > ZB:
        h = min(HMAX, 400.0 + G * abs(z[-1]))
        z.append(max(ZB, z[-1] - h))
    return np.array(z)
ZL = zlevels()
print(f"[walls] {len(ZL)} z levels, spacing {abs(np.diff(ZL)).min():.0f}..{abs(np.diff(ZL)).max():.0f} m")

allP = [TP]; allT = [TT]; nP = len(TP)
part = [np.full(len(TT), 1, np.int8)]          # 1 = top

# Vertical columns are created ONCE PER BOUNDARY NODE and shared.  Building them
# per-wall duplicates the corner columns that two walls have in common, which
# leaves the hull open along four vertical seams (measured: 232 edges used once,
# Euler -2 instead of +2).
col = {}
for v in bv:
    ids = [int(v)]
    pts = []
    for k in range(1, len(ZL)):
        ids.append(nP + len(pts)); pts.append((TP[v, 0], TP[v, 1], ZL[k]))
    allP.append(np.array(pts, float)); nP += len(pts)
    col[int(v)] = np.array(ids, np.int64)
print(f"[walls] {len(col):,} shared vertical columns x {len(ZL)} levels")

bottom_edges = []
for i in range(4):
    idx = sides[i]
    ns = len(idx)
    rows = [np.array([col[int(v)][k] for v in idx], np.int64) for k in range(len(ZL))]
    tri = []
    for k in range(len(ZL) - 1):
        a, b = rows[k], rows[k + 1]
        for j in range(ns - 1):
            tri.append((a[j], a[j + 1], b[j]))
            tri.append((a[j + 1], b[j + 1], b[j]))
    allT.append(np.array(tri, np.int64)); part.append(np.full(len(tri), 5, np.int8))
    bottom_edges.append(rows[-1])
    print(f"  wall {i}: {len(tri):,} triangles")

# ---- bottom: flat, boundary pinned to the walls' bottom rows ---------------
loop = []
for i in range(4):
    r = bottom_edges[i]
    loop.extend(list(r[:-1]))
loop = np.array(loop, np.int64)
P_all = np.vstack(allP)
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("bot"); geo = gmsh.model.geo
gp = [geo.addPoint(float(P_all[v, 0]), float(P_all[v, 1]), ZB, HMAX) for v in loop]
gl = [geo.addLine(gp[k], gp[(k + 1) % len(gp)]) for k in range(len(gp))]
sb = geo.addPlaneSurface([geo.addCurveLoop(gl)])
geo.synchronize()
for l in gl: gmsh.model.mesh.setTransfiniteCurve(l, 2)
for o in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature",
          "Mesh.MeshSizeExtendFromBoundary"): gmsh.option.setNumber(o, 0)
gmsh.option.setNumber("Mesh.MeshSizeMin", HMAX*0.8)
gmsh.option.setNumber("Mesh.MeshSizeMax", HMAX)
gmsh.option.setNumber("Mesh.Algorithm", 5)
gmsh.model.mesh.generate(2)
nt_, nc, _ = gmsh.model.mesh.getNodes()
BP = np.asarray(nc, float).reshape(-1, 3)
order = {int(t): i for i, t in enumerate(nt_)}
et, _, en = gmsh.model.mesh.getElements(2, sb)
BT = None
for et_, en_ in zip(et, en):
    if et_ == 2: BT = np.array([order[int(q)] for q in en_], np.int64).reshape(-1, 3)
gmsh.finalize()
# weld the bottom's boundary nodes onto the existing wall-bottom nodes
from scipy.spatial import cKDTree
tree = cKDTree(P_all[loop, :2])
dd, kk = tree.query(BP[:, :2], k=1, distance_upper_bound=1.0)
hit = np.isfinite(dd)
newid = np.empty(len(BP), np.int64)
newid[hit] = loop[kk[hit]]
nnew = int((~hit).sum())
newid[~hit] = np.arange(nP, nP + nnew)
allP.append(BP[~hit]); nP += nnew
allT.append(newid[BT]); part.append(np.full(len(BT), 3, np.int8))
print(f"[bottom] {len(BT):,} triangles, {nnew:,} interior nodes, "
      f"{int(hit.sum()):,} welded to walls")

P = np.vstack(allP); T = np.vstack(allT); PART = np.concatenate(part)
print(f"\n[PLC hull] {len(T):,} triangles, {len(P):,} vertices")
ee = np.sort(np.concatenate([T[:, [0,1]], T[:, [1,2]], T[:, [0,2]]]), axis=1)
uu, cc = np.unique(ee, axis=0, return_counts=True)
V_, E_, F_ = len(np.unique(T)), len(uu), len(T)
print(f"  edges used once {int((cc==1).sum()):,}   twice {int((cc==2).sum()):,}   >2 {int((cc>2).sum()):,}")
print(f"  Euler V-E+F = {V_-E_+F_}  (2 = closed sphere-topology shell)")
np.savez_compressed("build_tmp/hull.npz", P=P, T=T, PART=PART)
print("[out] build_tmp/hull.npz")
