#!/usr/bin/env python3
"""Compute the fault mesh-quality parameter r_v from Eq. (18) of

    Zhang et al. (2023), "A Mixed-Flux-Based Nodal Discontinuous Galerkin
    Method for 3D Dynamic Rupture Modeling", JGR Solid Earth, 128, e2022JB025817.

Eq. (18):   r_v = max{ V_A/V_B , V_B/V_A }

where A and B are the pair of neighbouring tetrahedra adjacent to the fault
surface and V_A, V_B are their volumes.  r_v = 1 means the two elements
straddling the fault have equal volume; the paper reports SSO instabilities
appear when r_v > 1.5.

This script reads a Gmsh 2.2 ASCII .msh file, identifies the fault triangles
(physical surface group "fault"), finds the two tetrahedra sharing each fault
triangle, and reports the distribution of r_v over all fault faces.
"""
import sys
from collections import defaultdict

MSH = "/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/" \
      "project_7.0_alternative/experimental_mesh_refinement/" \
      "safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"

FAULT_PHYS = 101   # physical tag of the "fault" surface group
ROCK_PHYS  = 1     # physical tag of the "rock" volume group
TRI  = 2           # Gmsh element type: 3-node triangle
TET  = 4           # Gmsh element type: 4-node tetrahedron


def tet_volume(p0, p1, p2, p3):
    """Signed-magnitude volume of a tetrahedron from its 4 vertices."""
    ax, ay, az = p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
    bx, by, bz = p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]
    cx, cy, cz = p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]
    # scalar triple product (a . (b x c)) / 6
    det = (ax*(by*cz - bz*cy)
         - ay*(bx*cz - bz*cx)
         + az*(bx*cy - by*cx))
    return abs(det) / 6.0


def main(path, fault_phys=FAULT_PHYS, rock_phys=ROCK_PHYS):
    FAULT_PHYS = fault_phys
    ROCK_PHYS = rock_phys
    nodes = {}
    fault_faces = set()       # frozenset of 3 node ids on the fault
    # face (sorted 3-node key) -> list of (tet_volume) for tets touching it
    face_to_tets = defaultdict(list)
    n_tri_fault = 0
    n_tet = 0

    with open(path) as f:
        line = f.readline()
        while line:
            tok = line.strip()
            if tok == "$Nodes":
                count = int(f.readline())
                for _ in range(count):
                    parts = f.readline().split()
                    nid = int(parts[0])
                    nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                f.readline()  # $EndNodes
            elif tok == "$Elements":
                count = int(f.readline())
                # First pass over elements: collect fault triangle faces.
                # We stream once, storing fault faces and tet connectivity,
                # then resolve volumes afterward (single pass, low memory).
                tets = []  # (v0,v1,v2,v3) node ids for rock tets
                for _ in range(count):
                    parts = f.readline().split()
                    etype = int(parts[1])
                    ntags = int(parts[2])
                    phys  = int(parts[3])
                    conn  = list(map(int, parts[3+ntags:]))
                    if etype == TRI and phys == FAULT_PHYS:
                        fault_faces.add(frozenset(conn))
                        n_tri_fault += 1
                    elif etype == TET and phys == ROCK_PHYS:
                        tets.append(conn)
                        n_tet += 1
                f.readline()  # $EndElements
                # Build adjacency only for faces that lie on the fault.
                for conn in tets:
                    v0, v1, v2, v3 = conn
                    p = (nodes[v0], nodes[v1], nodes[v2], nodes[v3])
                    vol = tet_volume(*p)
                    # 4 triangular faces of the tet
                    for face in ((v0, v1, v2), (v0, v1, v3),
                                 (v0, v2, v3), (v1, v2, v3)):
                        key = frozenset(face)
                        if key in fault_faces:
                            face_to_tets[key].append(vol)
                break
            line = f.readline()

    # Compute r_v per fault face.
    ratios = []
    n_shared2 = n_other = 0
    for key in fault_faces:
        vols = face_to_tets.get(key, [])
        if len(vols) == 2:
            vA, vB = vols
            r = max(vA/vB, vB/vA)
            ratios.append(r)
            n_shared2 += 1
        else:
            n_other += 1

    ratios.sort()
    n = len(ratios)
    print(f"Mesh: {path.split('/')[-1]}")
    print(f"  nodes                 : {len(nodes):,}")
    print(f"  rock tetrahedra       : {n_tet:,}")
    print(f"  fault triangles       : {n_tri_fault:,}")
    print(f"  fault faces w/ 2 tets : {n_shared2:,}")
    print(f"  fault faces w/ != 2   : {n_other:,}")
    if n == 0:
        print("No interior fault faces found (fault may be split/duplicated).")
        return
    import statistics
    mean = statistics.fmean(ratios)
    median = ratios[n//2]
    p95 = ratios[int(0.95*(n-1))]
    p99 = ratios[int(0.99*(n-1))]
    rmax = ratios[-1]
    n_gt15 = sum(1 for r in ratios if r > 1.5)
    n_gt20 = sum(1 for r in ratios if r > 2.0)
    print()
    print("  Eq.(18)  r_v = max{V_A/V_B, V_B/V_A}  over fault faces:")
    print(f"    min      : {ratios[0]:.6f}")
    print(f"    mean     : {mean:.6f}")
    print(f"    median   : {median:.6f}")
    print(f"    95th pct : {p95:.6f}")
    print(f"    99th pct : {p99:.6f}")
    print(f"    MAX r_v  : {rmax:.6f}")
    print()
    print(f"    fault faces with r_v > 1.5 : {n_gt15:,} ({100.0*n_gt15/n:.2f}%)")
    print(f"    fault faces with r_v > 2.0 : {n_gt20:,} ({100.0*n_gt20/n:.2f}%)")


if __name__ == "__main__":
    # usage: compute_eq18_volume_ratio.py [mesh.msh] [fault_phys] [rock_phys]
    p = sys.argv[1] if len(sys.argv) > 1 else MSH
    fp = int(sys.argv[2]) if len(sys.argv) > 2 else FAULT_PHYS
    rp = int(sys.argv[3]) if len(sys.argv) > 3 else ROCK_PHYS
    main(p, fp, rp)
