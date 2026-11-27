#!/usr/bin/env python3
"""Structured triangle mesh of a column [0, W] x [0, L] in Gmsh 2.2 ASCII, for the
poroelastic verification cases (tests/verification/poro/*.toml).

  python3 make_column_mesh.py --L 0.01 --W 0.0005 --ny 40 --nx 2 -o column_L10_h0p25.msh

Physical curves: bottom (y = 0), right (x = W), top (y = L), left (x = 0);
physical surface: column. Each cell is split along alternating diagonals.
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--L", type=float, required=True)
ap.add_argument("--W", type=float, required=True)
ap.add_argument("--ny", type=int, required=True)
ap.add_argument("--nx", type=int, required=True)
ap.add_argument("-o", required=True)
a = ap.parse_args()

nid = lambda i, j: j * (a.nx + 1) + i + 1
nodes = [(a.W * i / a.nx, a.L * j / a.ny) for j in range(a.ny + 1) for i in range(a.nx + 1)]
els = []
for i in range(a.nx):
    els.append((1, 1, (nid(i, 0), nid(i + 1, 0))))
for j in range(a.ny):
    els.append((1, 2, (nid(a.nx, j), nid(a.nx, j + 1))))
for i in range(a.nx):
    els.append((1, 3, (nid(i + 1, a.ny), nid(i, a.ny))))
for j in range(a.ny):
    els.append((1, 4, (nid(0, j + 1), nid(0, j))))
for j in range(a.ny):
    for i in range(a.nx):
        n00, n10, n11, n01 = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
        if (i + j) % 2 == 0:
            els += [(2, 5, (n00, n10, n11)), (2, 5, (n00, n11, n01))]
        else:
            els += [(2, 5, (n00, n10, n01)), (2, 5, (n10, n11, n01))]
with open(a.o, "w") as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n5\n")
    f.write('1 1 "bottom"\n1 2 "right"\n1 3 "top"\n1 4 "left"\n2 5 "column"\n$EndPhysicalNames\n')
    f.write("$Nodes\n%d\n" % len(nodes))
    for k, (x, y) in enumerate(nodes):
        f.write("%d %.17g %.17g 0\n" % (k + 1, x, y))
    f.write("$EndNodes\n$Elements\n%d\n" % len(els))
    for k, (t, phys, c) in enumerate(els):
        f.write("%d %d 2 %d %d %s\n" % (k + 1, t, phys, phys, " ".join(map(str, c))))
    f.write("$EndElements\n")
print("wrote", a.o, len(nodes), "nodes", len(els), "elements")
