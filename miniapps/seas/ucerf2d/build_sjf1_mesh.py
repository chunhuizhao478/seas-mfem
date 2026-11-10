#!/usr/bin/env python3
# build_sjf1_mesh.py — the San Jacinto trace sjf1 (Anza) of the UCERF 2-D R 70 km
# mesh, alone, in a box of the same tilt and wall distance (Chunhui Zhao,
# 2026-11-09: "can we just build this San Jacino fault (this trace alone) and
# load it? can we reproduce this issue").
#
# The fault vertices are read from the committed full mesh
# (ucerf2d/mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh, tag 3): the chain
# that starts at (160.907, -330.071) km, 347 vertices, 173.06 km.  So the
# trace is that mesh's sjf1, vertex for vertex (no new smoothing).  The box
# is the full mesh's box (the four corners in its json), so the walls, the
# plate loading and the PML band are those of the full model; sizing and
# tags as build_mesh_smooth.py made the full mesh
# (/mnt/project-files/ucerf-2d/smoothed-mesh-2026-11-01/): 500 m on the
# fault out to 10 km, then growing by 0.1 to 5 km; tags 3 fault, 5 the two
# long walls, 1 the two ends; gmsh 2.2.
#
#   python3 build_sjf1_mesh.py FULL.msh FULL.json OUT.msh
import sys, json, math, collections
import numpy as np, gmsh
full, fjson, out = sys.argv[1:4]
LCF, LCFAR, GRADE, D0, HMAX = 500.0, 40000.0, 0.1, 10e3, 5e3
START = (160906.78568, -330070.69813)

L = open(full).read().split('\n'); i = L.index('$Nodes'); n = int(L[i+1]); xy = {}
for l in L[i+2:i+2+n]:
    p = l.split(); xy[int(p[0])] = (float(p[1]), float(p[2]))
j = L.index('$Elements'); m = int(L[j+1]); adj = collections.defaultdict(list)
for l in L[j+2:j+2+m]:
    p = l.split()
    if p[1] == '1' and int(p[3]) == 3:
        a, b = int(p[-2]), int(p[-1]); adj[a].append(b); adj[b].append(a)
tips = [k for k, v in adj.items() if len(v) == 1]
s0 = min(tips, key=lambda k: math.dist(xy[k], START))
assert math.dist(xy[s0], START) < 1.0, "sjf1's first tip not found"
chain = [s0]; seen = {s0}
while True:
    nx = [v for v in adj[chain[-1]] if v not in seen]
    if not nx: break
    chain.append(nx[0]); seen.add(nx[0])
pl = [xy[k] for k in chain]
flen = sum(math.dist(pl[k], pl[k+1]) for k in range(len(pl)-1))
assert len(pl) == 347, len(pl)

T = json.load(open(fjson))['tilt']; e1 = tuple(T['e1']); e2 = tuple(T['e2'])
corners = [tuple(c) for c in T['corners']]          # the full model's box, unchanged
A = [c[0]*e1[0]+c[1]*e1[1] for c in corners]; B = [c[0]*e2[0]+c[1]*e2[1] for c in corners]
a0, a1, b0, b1 = min(A), max(A), min(B), max(B)
xi = [p[0]*e1[0]+p[1]*e1[1] for p in pl]; et = [p[0]*e2[0]+p[1]*e2[1] for p in pl]
wall = min(min(xi)-a0, a1-max(xi), min(et)-b0, b1-max(et))

gmsh.initialize(); gmsh.option.setNumber('General.Verbosity', 2)
gmsh.model.add('sjf1'); occ = gmsh.model.occ
cp = [occ.addPoint(c[0], c[1], 0) for c in corners]
cl = [occ.addLine(cp[k], cp[(k+1) % 4]) for k in range(4)]
rect = occ.addPlaneSurface([occ.addCurveLoop(cl)])
pt = [occ.addPoint(p[0], p[1], 0, LCF) for p in pl]
cur = [occ.addLine(pt[k], pt[k+1]) for k in range(len(pt)-1)]
o, omap = occ.fragment([(2, rect)], [(1, c) for c in cur]); occ.synchronize()
surf = [e[1] for e in gmsh.model.getEntities(2)]
fault = sorted({e[1] for mm in omap[1:] for e in mm if e[0] == 1})
bnd = [e[1] for e in gmsh.model.getBoundary([(2, s) for s in surf], oriented=False) if e[1] not in fault]
fset = {(round(p[0], 3), round(p[1], 3)) for p in pl}
for e in gmsh.model.getEntities(0):
    x, y, _ = gmsh.model.getValue(0, e[1], [])
    gmsh.model.mesh.setSize([e], LCF if (round(x, 3), round(y, 3)) in fset else LCFAR)
fld = gmsh.model.mesh.field
fd = fld.add('Distance'); fld.setNumbers(fd, 'CurvesList', fault)
fld.setNumber(fd, 'Sampling', int(max(occ.getMass(1, c) for c in fault)/(0.5*LCF))+2)
ft = fld.add('Threshold'); fld.setNumber(ft, 'InField', fd)
fld.setNumber(ft, 'SizeMin', LCF); fld.setNumber(ft, 'SizeMax', HMAX)
fld.setNumber(ft, 'DistMin', D0); fld.setNumber(ft, 'DistMax', D0+(HMAX-LCF)/GRADE)
fld.setAsBackgroundMesh(ft)
for k in ('MeshSizeExtendFromBoundary', 'MeshSizeFromPoints', 'MeshSizeFromCurvature'):
    gmsh.option.setNumber('Mesh.'+k, 0)
gmsh.option.setNumber('General.NumThreads', 1)     # one thread: the same mesh on a rebuild
def is_end(c):
    bb = gmsh.model.getBoundingBox(1, c); mx = ((bb[0]+bb[3])/2, (bb[1]+bb[4])/2)
    A = mx[0]*e1[0]+mx[1]*e1[1]; return abs(A-a0) < 1 or abs(A-a1) < 1
ends = [c for c in bnd if is_end(c)]; walls = [c for c in bnd if c not in ends]
gmsh.model.addPhysicalGroup(1, fault, 3); gmsh.model.addPhysicalGroup(1, walls, 5)
gmsh.model.addPhysicalGroup(1, ends, 1); gmsh.model.addPhysicalGroup(2, surf, 1)
gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
gmsh.model.mesh.generate(2); gmsh.write(out)
nn = len(gmsh.model.mesh.getNodes()[0]); nt = len(gmsh.model.mesh.getElementsByType(2)[0])
fe = sum(len(t) for c in fault for t in gmsh.model.mesh.getElements(1, c)[1])
res = dict(source=full, chain_vertices=len(pl), fault_length_km=flen/1e3, fault_edges=fe, nodes=nn, triangles=nt,
           domain_km=[(a1-a0)/1e3, (b1-b0)/1e3], corners=corners, e1=e1, e2=e2, min_wall_distance_km=wall/1e3,
           tips=[pl[0], pl[-1]], lcf=LCF, grade=GRADE, d0=D0, hmax=HMAX)
print(json.dumps(res)); json.dump(res, open(out.replace('.msh', '.json'), 'w'))
gmsh.finalize()
