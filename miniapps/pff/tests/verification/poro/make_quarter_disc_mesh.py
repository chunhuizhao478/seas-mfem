#!/usr/bin/env python3
"""Triangle mesh of a quarter disc (x, y >= 0, r <= a) in Gmsh 2.2 ASCII, for the
de Leeuw cylinder case (tests/verification/poro/deleeuw_cylinder.toml).

  python3 make_quarter_disc_mesh.py --a 0.01 --nr 24 -o quarter_disc_a10_nr24.msh

Ring i (radius i a / nr) carries 2 i + 1 nodes. Physical curves: bottom (y = 0),
outer (r = a), left (x = 0); physical surface: disc.
"""
import argparse, math

ap = argparse.ArgumentParser()
ap.add_argument("--a", type=float, required=True)
ap.add_argument("--nr", type=int, required=True)
ap.add_argument("-o", required=True)
A = ap.parse_args()
q = 2
nodes = [(0.0, 0.0)]
ring = [[1]]
for i in range(1, A.nr + 1):
    ids = []
    for j in range(q * i + 1):
        th = 0.5 * math.pi * j / (q * i)
        x, y = A.a * i / A.nr * math.cos(th), A.a * i / A.nr * math.sin(th)
        if j == 0: y = 0.0
        if j == q * i: x = 0.0
        nodes.append((x, y))
        ids.append(len(nodes))
    ring.append(ids)
tris = []
for i in range(1, A.nr + 1):
    inner, outer = ring[i - 1], ring[i]
    a, b = 0, 0
    while a < len(inner) - 1 or b < len(outer) - 1:
        # advance on the ring whose next node has the smaller angle
        ta = (a + 1) / max(len(inner) - 1, 1) if a < len(inner) - 1 else 2.0
        tb = (b + 1) / (len(outer) - 1) if b < len(outer) - 1 else 2.0
        if tb <= ta:
            tris.append((inner[a], outer[b], outer[b + 1])); b += 1
        else:
            tris.append((inner[a], outer[b], inner[a + 1])); a += 1
els = []
for i in range(A.nr):
    els.append((1, 1, (ring[i][0], ring[i + 1][0])))
o = ring[A.nr]
for j in range(len(o) - 1):
    els.append((1, 2, (o[j], o[j + 1])))
for i in range(A.nr, 0, -1):
    els.append((1, 3, (ring[i][-1], ring[i - 1][-1])))
els += [(2, 4, t) for t in tris]
with open(A.o, "w") as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n")
    f.write('1 1 "bottom"\n1 2 "outer"\n1 3 "left"\n2 4 "disc"\n$EndPhysicalNames\n')
    f.write("$Nodes\n%d\n" % len(nodes))
    for k, (x, y) in enumerate(nodes):
        f.write("%d %.17g %.17g 0\n" % (k + 1, x, y))
    f.write("$EndNodes\n$Elements\n%d\n" % len(els))
    for k, (t, phys, c) in enumerate(els):
        f.write("%d %d 2 %d %d %s\n" % (k + 1, t, phys, phys, " ".join(map(str, c))))
    f.write("$EndElements\n")
area = sum(0.5 * abs((nodes[b-1][0]-nodes[a-1][0])*(nodes[c-1][1]-nodes[a-1][1]) - (nodes[c-1][0]-nodes[a-1][0])*(nodes[b-1][1]-nodes[a-1][1])) for a, b, c in tris)
print("wrote", A.o, len(nodes), "nodes", len(tris), "triangles, area / (pi a^2 / 4) =", area / (math.pi * A.a**2 / 4))
