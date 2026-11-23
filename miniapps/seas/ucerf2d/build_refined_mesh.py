#!/usr/bin/env python3
# build_refined_mesh.py — the UCERF 2-D R 70 km mesh with a 250 m fault band
# (heterogeneity plan rev5, D5; Chunhui Zhao, 2026-11-23: "I would choose
# refine the fault band to 250 m").  Why: with the so2 stress (sigma_n' about
# 180 MPa) the dynamic process zone of H4 was 3.4-4 cells of 500 m (median)
# and 1.4-3 cells at the 5th percentile, below the 3 cells of D5.
#
# The fault traces are the committed 500 m mesh's, vertex for vertex (tag 3),
# with one new vertex at the middle of every edge: the same polylines, edges
# of 250 m.  The box is the 500 m mesh's (the four corners in its json), so
# walls, plate loading and PML band are unchanged.  Size h(d) at distance d
# from the faults: 250 m out to 5 km, then linear to 500 m at 10 km, then the
# 500 m mesh's ramp (500 m + 0.1 (d - 10 km), capped at 5 km); gmsh's Max of
# two Threshold fields.  Tags 3 fault, 5 the two long walls, 1 the two ends.
#
#   python3 build_refined_mesh.py FULL500.msh FULL500.json OUT.msh [--chain-start X Y]
# --chain-start keeps only the chain with a tip within 1 m of (X, Y) (sjf1:
# 160906.78568 -330070.69813), as build_sjf1_mesh.py.
import sys, json, math, argparse, collections
import numpy as np, gmsh
ap = argparse.ArgumentParser()
ap.add_argument('full'); ap.add_argument('fjson'); ap.add_argument('out')
ap.add_argument('--chain-start', nargs=2, type=float, default=None)
a = ap.parse_args()
LCF, LC1, D1, D2, GRADE, HMAX, LCFAR = 250.0, 500.0, 5e3, 10e3, 0.1, 5e3, 40000.0

L = open(a.full).read().split('\n'); i = L.index('$Nodes'); n = int(L[i+1]); xy = {}
for l in L[i+2:i+2+n]:
    p = l.split(); xy[int(p[0])] = (float(p[1]), float(p[2]))
j = L.index('$Elements'); m = int(L[j+1]); adj = collections.defaultdict(list)
for l in L[j+2:j+2+m]:
    p = l.split()
    if p[1] == '1' and int(p[3]) == 3:
        u, v = int(p[-2]), int(p[-1]); adj[u].append(v); adj[v].append(u)
assert all(len(v) <= 2 for v in adj.values()), "a fault vertex of degree 3"
tips = sorted(k for k, v in adj.items() if len(v) == 1)
chains, seen = [], set()
for t0 in tips:
    if t0 in seen: continue
    ch = [t0]; seen.add(t0)
    while True:
        nx = [v for v in adj[ch[-1]] if v not in seen]
        if not nx: break
        ch.append(nx[0]); seen.add(nx[0])
    chains.append([xy[k] for k in ch])
assert len(seen) == len(adj), "a closed fault loop"
if a.chain_start:
    S = tuple(a.chain_start)
    chains = [c for c in chains if min(math.dist(c[0], S), math.dist(c[-1], S)) < 1.0]
    assert len(chains) == 1, "chain start not found"
# one new vertex at the middle of every edge
fine = []
for c in chains:
    f = [c[0]]
    for p, q in zip(c[:-1], c[1:]):
        f.append(((p[0]+q[0])/2, (p[1]+q[1])/2)); f.append(q)
    fine.append(f)
flen = sum(math.dist(p, q) for c in fine for p, q in zip(c[:-1], c[1:]))

T = json.load(open(a.fjson))['tilt']; e1 = tuple(T['e1']); e2 = tuple(T['e2'])
corners = [tuple(c) for c in T['corners']]
A = [c[0]*e1[0]+c[1]*e1[1] for c in corners]; B = [c[0]*e2[0]+c[1]*e2[1] for c in corners]
a0, a1, b0, b1 = min(A), max(A), min(B), max(B)

gmsh.initialize(); gmsh.option.setNumber('General.Verbosity', 2)
gmsh.model.add('ucerf2d_250m'); occ = gmsh.model.occ
cp = [occ.addPoint(c[0], c[1], 0) for c in corners]
cl = [occ.addLine(cp[k], cp[(k+1) % 4]) for k in range(4)]
rect = occ.addPlaneSurface([occ.addCurveLoop(cl)])
cur = []
for c in fine:
    pt = [occ.addPoint(p[0], p[1], 0, LCF) for p in c]
    cur += [occ.addLine(pt[k], pt[k+1]) for k in range(len(pt)-1)]
o, omap = occ.fragment([(2, rect)], [(1, c) for c in cur]); occ.synchronize()
surf = [e[1] for e in gmsh.model.getEntities(2)]
fault = sorted({e[1] for mm in omap[1:] for e in mm if e[0] == 1})
bnd = [e[1] for e in gmsh.model.getBoundary([(2, s) for s in surf], oriented=False) if e[1] not in fault]
fset = {(round(p[0], 3), round(p[1], 3)) for c in fine for p in c}
for e in gmsh.model.getEntities(0):
    x, y, _ = gmsh.model.getValue(0, e[1], [])
    gmsh.model.mesh.setSize([e], LCF if (round(x, 3), round(y, 3)) in fset else LCFAR)
fld = gmsh.model.mesh.field
fd = fld.add('Distance'); fld.setNumbers(fd, 'CurvesList', fault)
fld.setNumber(fd, 'Sampling', int(max(occ.getMass(1, c) for c in fault)/(0.5*LCF))+2)
# g1: 250 m to 5 km, linear to 500 m at 10 km, 500 m beyond
g1 = fld.add('Threshold'); fld.setNumber(g1, 'InField', fd)
fld.setNumber(g1, 'SizeMin', LCF); fld.setNumber(g1, 'SizeMax', LC1)
fld.setNumber(g1, 'DistMin', D1); fld.setNumber(g1, 'DistMax', D2)
# g2: the 500 m mesh's ramp, 500 m + 0.1 (d - 10 km) to 5 km (below 500 m inside 10 km)
g2 = fld.add('Threshold'); fld.setNumber(g2, 'InField', fd)
d2min = D2 - (LC1 - LCF)/GRADE
fld.setNumber(g2, 'SizeMin', LCF); fld.setNumber(g2, 'SizeMax', HMAX)
fld.setNumber(g2, 'DistMin', d2min); fld.setNumber(g2, 'DistMax', d2min + (HMAX-LCF)/GRADE)
gm = fld.add('Max'); fld.setNumbers(gm, 'FieldsList', [g1, g2])
fld.setAsBackgroundMesh(gm)
for k in ('MeshSizeExtendFromBoundary', 'MeshSizeFromPoints', 'MeshSizeFromCurvature'):
    gmsh.option.setNumber('Mesh.'+k, 0)
gmsh.option.setNumber('General.NumThreads', 1)     # one thread: the same mesh on a rebuild
def is_end(c):
    bb = gmsh.model.getBoundingBox(1, c); mx = ((bb[0]+bb[3])/2, (bb[1]+bb[4])/2)
    Ax = mx[0]*e1[0]+mx[1]*e1[1]; return abs(Ax-a0) < 1 or abs(Ax-a1) < 1
ends = [c for c in bnd if is_end(c)]; walls = [c for c in bnd if c not in ends]
gmsh.model.addPhysicalGroup(1, fault, 3); gmsh.model.addPhysicalGroup(1, walls, 5)
gmsh.model.addPhysicalGroup(1, ends, 1); gmsh.model.addPhysicalGroup(2, surf, 1)
gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
gmsh.model.mesh.generate(2); gmsh.write(a.out)
nn = len(gmsh.model.mesh.getNodes()[0])
et, tags, nodes = gmsh.model.mesh.getElements(2)
nodes = np.array(nodes[0]).reshape(-1, 3)
ntag, coord, _ = gmsh.model.mesh.getNodes(); C = np.zeros((int(ntag.max())+1, 3)); C[ntag.astype(int)] = coord.reshape(-1, 3)
P0 = C[nodes[:, 0].astype(int)]; P1 = C[nodes[:, 1].astype(int)]; P2 = C[nodes[:, 2].astype(int)]
area = 0.5*np.abs((P1[:, 0]-P0[:, 0])*(P2[:, 1]-P0[:, 1])-(P1[:, 1]-P0[:, 1])*(P2[:, 0]-P0[:, 0]))
per = np.linalg.norm(P1-P0, axis=1)+np.linalg.norm(P2-P1, axis=1)+np.linalg.norm(P0-P2, axis=1)
he = 4*area/per
q = gmsh.model.mesh.getElementQualities(tags[0].tolist(), 'gamma')
fe = sum(len(t) for c in fault for t in gmsh.model.mesh.getElements(1, c)[1])
res = dict(source=a.full, chains=len(fine), fault_length_km=flen/1e3, fault_edges=fe, nodes=nn, triangles=len(tags[0]),
           he_min=float(he.min()), he_p01=float(np.percentile(he, 1)), he_med=float(np.median(he)), he_max=float(he.max()),
           gamma_min=float(min(q)), dof_dg_p1=6*len(tags[0]), corners=corners, e1=e1, e2=e2,
           tips=[[c[0], c[-1]] for c in fine], lcf=LCF, d_fine=D1, d_500=D2, grade=GRADE, hmax=HMAX,
           tilt=dict(T))
print(json.dumps({k: v for k, v in res.items() if k not in ('corners', 'tips', 'tilt')}))
json.dump(res, open(a.out.replace('.msh', '.json'), 'w'))
gmsh.finalize()
