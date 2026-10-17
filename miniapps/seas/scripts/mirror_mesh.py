#!/usr/bin/env python3
"""Reflect a half tetrahedral mesh across the plane Y = 0 into a mirror-symmetric mesh.

Chunhui Zhao, 2026-10-17.  For a planar vertical fault on Y = 0 in a setup that is symmetric
about the fault plane (uniform or Y-symmetric material, prestress and boundaries), the exact
fault normal-stress change is zero.  On a mesh whose tets differ on the two sides of a fault
face the discrete one is not: it grows with the slip gradient and peaks at the fault edges
(TPV104 100 m, Vista 1076418, went tensile at the fault end).  Meshing the Y >= 0 half with
gmsh and reflecting it here gives a mesh in which every tet on the -Y side is the exact mirror
image of one on the +Y side, so the discrete normal-stress change is zero as well.
The same half + mirror construction as
safs/project_7.0_alternative/experimental_mesh_refinement/make_symmetric_unstructured_mesh.py,
which builds the TPV102 geometry in the script; this one takes the half mesh of any .geo
(tpv104/mesh/tpv104_half.geo, tpv205/mesh/tpv2053d_half.geo) and keeps its physical groups.

    python3 scripts/mirror_mesh.py half.msh full.msh [--plane-attr 99]

Input: gmsh msh 2.2 ASCII (`gmsh -3 ... -format msh22`) of the Y >= 0 half, with physical
groups: tets (type 4), boundary triangles (type 2); the fault triangles and the rest of the
Y = 0 plane (--plane-attr, default 99) lie on Y = 0, every other node has Y > 0.
Output: msh 2.2 with each node off the plane duplicated at -Y, each tet and each boundary
triangle off the plane duplicated (reflected, re-oriented to positive volume / the same
outward sense), the fault triangles kept once (they are interior faces now) and the
--plane-attr triangles dropped (welded interior faces).  Physical names are copied.
Elements of other types (points, lines) are dropped.  Checks: plane nodes |Y| <= tol (then
set to exactly 0), positive tet volumes, and that every fault triangle has a tet on each side.
"""

import argparse
import sys


def read_msh22(path):
    with open(path) as f:
        lines = f.read().split("\n")
    def section(name):
        i = lines.index("$" + name)
        j = lines.index("$End" + name)
        return lines[i + 1:j]
    fmt = section("MeshFormat")[0].split()
    if not fmt[0].startswith("2.2") or fmt[1] != "0":
        sys.exit(f"mirror_mesh: {path}: need msh 2.2 ASCII, got MeshFormat '{' '.join(fmt)}'")
    names = section("PhysicalNames")[1:] if "$PhysicalNames" in lines else []
    nl = section("Nodes")
    nodes = {}
    for l in nl[1:int(nl[0]) + 1]:
        s = l.split()
        nodes[int(s[0])] = [float(s[1]), float(s[2]), float(s[3])]
    el = section("Elements")
    elems = []
    for l in el[1:int(el[0]) + 1]:
        s = list(map(int, l.split()))
        typ, ntag = s[1], s[2]
        elems.append((typ, s[3], s[4] if ntag > 1 else s[3], s[3 + ntag:]))
    return names, nodes, elems


def vol6(X, a, b, c, d):
    pa, pb, pc, pd = X[a], X[b], X[c], X[d]
    e1 = [pb[k] - pa[k] for k in range(3)]
    e2 = [pc[k] - pa[k] for k in range(3)]
    e3 = [pd[k] - pa[k] for k in range(3)]
    return (e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
            + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("half")
    p.add_argument("out")
    p.add_argument("--plane-attr", type=int, default=99,
                   help="physical tag of the non-fault part of Y = 0 (dropped; default 99)")
    p.add_argument("--fault-attr", type=int, default=None,
                   help="physical tag of the fault triangles (default: every Y = 0 triangle "
                        "not tagged --plane-attr)")
    p.add_argument("--tol", type=float, default=1e-6,
                   help="nodes with |Y| <= tol x (mesh extent) are plane nodes (default 1e-6)")
    a = p.parse_args()

    names, X, elems = read_msh22(a.half)
    ext = max(max(abs(c) for c in x) for x in X.values())
    tol = a.tol * ext
    mirror = {}
    for v, x in X.items():
        if abs(x[1]) <= tol:
            x[1] = 0.0
            mirror[v] = v
        elif x[1] < 0:
            sys.exit(f"mirror_mesh: node {v} at Y = {x[1]} < 0; the input must be the Y >= 0 half")
    nid = max(X) + 1
    for v in sorted(X):
        if v not in mirror:
            x = X[v]
            mirror[v] = nid
            X[nid] = [x[0], -x[1], x[2]]
            nid += 1

    tris, tets = [], []
    n_plane = n_fault = n_flip = 0
    fault_keys = []
    for typ, phys, geom, vs in elems:
        if typ == 2:
            on_plane = all(mirror[v] == v for v in vs)
            if on_plane:
                if phys == a.plane_attr:
                    n_plane += 1
                    continue
                if a.fault_attr is not None and phys != a.fault_attr:
                    sys.exit(f"mirror_mesh: triangle {vs} on Y = 0 has tag {phys}, neither "
                             f"--fault-attr {a.fault_attr} nor --plane-attr {a.plane_attr}")
                tris.append((phys, geom, vs))
                fault_keys.append(tuple(sorted(vs)))
                n_fault += 1
            else:
                tris.append((phys, geom, vs))
                tris.append((phys, geom, [mirror[vs[0]], mirror[vs[2]], mirror[vs[1]]]))
        elif typ == 4:
            if vol6(X, *vs) < 0:
                vs = [vs[0], vs[1], vs[3], vs[2]]
                n_flip += 1
            mv = [mirror[vs[0]], mirror[vs[1]], mirror[vs[3]], mirror[vs[2]]]
            v0, v1 = vol6(X, *vs), vol6(X, *mv)
            if not (v0 > 0 and v1 > 0):
                sys.exit(f"mirror_mesh: degenerate tet {vs} (6V = {v0}, mirror {v1})")
            tets.append((phys, geom, vs))
            tets.append((phys, geom, mv))

    # every fault triangle must be a face of exactly one tet on each side
    want = set(fault_keys)
    sides = {k: 0 for k in want}
    for _, _, t in tets:
        a0, b, c, d = t
        for f in ((b, c, d), (a0, d, c), (a0, b, d), (a0, c, b)):
            k = tuple(sorted(f))
            if k in want:
                sides[k] += 1
    bad = [k for k, n in sides.items() if n != 2]
    if bad:
        sys.exit(f"mirror_mesh: {len(bad)} fault triangles do not have one tet on each side "
                 f"(first {bad[0]})")

    with open(a.out, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        if names:
            keep = [l for l in names if l.strip() and int(l.split()[1]) != a.plane_attr]
            f.write("$PhysicalNames\n%d\n%s\n$EndPhysicalNames\n" % (len(keep), "\n".join(keep)))
        f.write("$Nodes\n%d\n" % len(X))
        for v in sorted(X):
            x = X[v]
            f.write("%d %.17g %.17g %.17g\n" % (v, x[0], x[1], x[2]))
        f.write("$EndNodes\n$Elements\n%d\n" % (len(tris) + len(tets)))
        k = 1
        for typ, lst in ((2, tris), (4, tets)):
            for phys, geom, vs in lst:
                f.write("%d %d 2 %d %d %s\n" % (k, typ, phys, geom, " ".join(map(str, vs))))
                k += 1
        f.write("$EndElements\n")
    n_bdr = len(tris) - n_fault
    print(f"mirror_mesh: {a.half} -> {a.out}: nodes {len(X)}, tets {len(tets)}, "
          f"fault triangles {n_fault} (each with a tet on both sides), other boundary "
          f"triangles {n_bdr}, dropped {n_plane} plane triangles (tag {a.plane_attr}), "
          f"re-oriented {n_flip} input tets")


if __name__ == "__main__":
    main()
