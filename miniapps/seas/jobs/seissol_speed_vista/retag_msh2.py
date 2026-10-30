#!/usr/bin/env python3
"""Copy a Gmsh 2.2 ASCII mesh with its boundary triangles' physical tags renamed for SeisSol.

SeisSol and MFEM speed test on Vista (plan rev1, Chunhui Zhao 2026-10-30).  The TPV205 mesh
tpv2053d_200m.msh tags its boundary triangles 101 (free surface), 103 (fault) and 105 (absorbing).
MFEM reads these tags; SeisSol at a8e7f04f reads a face tag 3 or >= 65 as a rupture face
(src/Geometry/FaceMap.cpp:98-117), so 101 and 105 would become rupture faces.  The deck option
that renames tags (facemapfile) does not work at that commit (src/Initializer/InitMesh.cpp:279-287
passes the file name to YAML::Load).  This script writes a copy whose triangles carry 1, 3, 5
instead, for pumgen; nothing else in the file changes (nodes, tetrahedra, element order).

Usage: python3 retag_msh2.py in.msh out.msh 101:1 103:3 105:5
Prints how many triangles each rename touched; exits 1 if the mesh is not Gmsh 2.2 ASCII or a
triangle carries a tag not in the map.
"""
import sys


def main(argv):
    src, dst, pairs = argv[1], argv[2], argv[3:]
    rename = {int(a): int(b) for a, b in (p.split(":") for p in pairs)}
    count = {k: 0 for k in rename}
    with open(src) as f, open(dst, "w") as g:
        section = None
        for line in f:
            s = line.strip()
            if s.startswith("$"):
                section = None if s.startswith("$End") else s
                if s == "$Elements":
                    g.write(line)
                    g.write(next(f))            # the element count
                    section = "elements"
                    continue
            elif section == "$MeshFormat" and not s.startswith("2.2 0"):
                sys.exit("ERROR: %s is not Gmsh 2.2 ASCII (format line '%s')" % (src, s))
            elif section == "elements":
                w = line.split()
                if w[1] == "2":                 # triangle: number, type, ntags, physical, elementary, ...
                    tag = int(w[3])
                    if tag not in rename:
                        sys.exit("ERROR: triangle %s has tag %d, not in the map" % (w[0], tag))
                    count[tag] += 1
                    w[3] = str(rename[tag])
                    line = " ".join(w) + "\n"
            g.write(line)
    for k, n in count.items():
        print("%s: %d triangles retagged %d -> %d" % (src, n, k, rename[k]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
