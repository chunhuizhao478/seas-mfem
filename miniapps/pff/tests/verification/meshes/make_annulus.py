#!/usr/bin/env python3
"""Annulus r_b <= r <= R (gmsh OCC, triangles, uniform h) for the C1 cavity gate.

Physical groups: curves "borehole" (r = r_b), "outer" (r = R); surface "solid".
  python3 make_annulus.py --rb 4.15e-3 --R 16.6e-3 --h 0.3e-3 --out annulus_rb4p15_R16p6_h0p3
"""
import argparse
import gmsh

p = argparse.ArgumentParser()
p.add_argument("--rb", type=float, default=4.15e-3)
p.add_argument("--R", type=float, default=16.6e-3)
p.add_argument("--h", type=float, default=0.3e-3)
p.add_argument("--out", required=True)
a = p.parse_args()

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.NumThreads", 1)
gmsh.model.add("annulus")
occ = gmsh.model.occ
outer = occ.addDisk(0, 0, 0, a.R, a.R)
inner = occ.addDisk(0, 0, 0, a.rb, a.rb)
occ.cut([(2, outer)], [(2, inner)])
occ.synchronize()
surf = gmsh.model.getEntities(2)[0][1]
bore, out = [], []
for _, c in gmsh.model.getEntities(1):
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, c)
    (bore if (x1 - x0) < 1.5 * a.rb * 2 else out).append(c)
gmsh.model.addPhysicalGroup(1, bore, name="borehole")
gmsh.model.addPhysicalGroup(1, out, name="outer")
gmsh.model.addPhysicalGroup(2, [surf], name="solid")
for k, v in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
             "Mesh.MeshSizeFromCurvature": 0, "Mesh.MeshSizeMin": a.h, "Mesh.MeshSizeMax": a.h,
             "Mesh.Algorithm": 6, "Mesh.MshFileVersion": 2.2, "Mesh.Binary": 0}.items():
    gmsh.option.setNumber(k, v)
gmsh.model.mesh.generate(2)
gmsh.write(a.out + ".msh")
print(a.out + ".msh", len(gmsh.model.mesh.getElementsByType(2)[0]), "triangles")
gmsh.finalize()
