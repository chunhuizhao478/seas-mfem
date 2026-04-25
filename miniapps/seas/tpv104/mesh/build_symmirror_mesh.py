#!/usr/bin/env python3
"""Build a Y-mirror-symmetric tetrahedral mesh for TPV104 in Gmsh v2.2 ASCII.

Why a Python script and not a Gmsh .geo file?
---------------------------------------------
Gmsh's tetrahedral mesher (Delaunay / Frontal) is not deterministic on
mirror-symmetric input geometry — the mesher visits elements in an order
that depends on internal data structures, so even with a perfectly
mirrored .geo input, the output tet topology is not guaranteed to be
mirror-symmetric.  See debug doc 2026-04-25_pm + the local 2-step test
on tpv104_symmirror_1000m.msh which proved that bulk Q+ vs Q− is bit-
exact symmetric only when the tet topology itself is constructed to be
mirror-symmetric (k+_SXY = k−_SXY exactly, k+_VX = −k−_VX exactly,
diff_SXX = 0 bit-exact at step 1).

This script constructs the mesh directly:
  * Regular hex grid with a vertex plane at Y=0 (so y=0 is shared).
  * Each hex split into 6 tets with diagonal v0–v7.
  * +Y hexes use Pattern A (diagonal from low-corner v0 to high-corner v7).
  * −Y hexes use Pattern B = y-mirror image of Pattern A (swap y bit).
  * Boundary triangles derived from tet topology so they always match
    the per-side face splits.
  * Tet vertices are reoriented to positive volume so MFEM doesn't
    abort on STable3D inconsistencies.

SCEC ↔ MFEM coordinates:
  SCEC x (strike)        ↔ code X
  SCEC y (depth ≥0)      ↔ code -Z   (free surface at Z=0; bulk Z<0)
  SCEC z (fault-normal)  ↔ code Y    (fault plane Y=0)

Usage
-----
  python3 build_symmirror_mesh.py [--dx DX] [--lx LX] [--ly LY] [--lz LZ]
                                  [--out OUT_MSH]

  Default: dx=200m, lx=16km, ly=4km, lz=16km, out=tpv104_symmirror_200m.msh

  Element count:
     N_tets = 6 × (LX/DX) × (LY/DY) × (LZ/DZ)
     with DY=DX and DZ=DX:
        16x4x16 km @ 200m  →  768,000 tets
        12x4x12 km @ 200m  →  432,000 tets
         4x4x8  km @ 1000m →     768 tets   (local 2-step test)

Physical groups (matches TPV104 driver defaults):
  Surface 1 = free surface (Z=0)
  Surface 3 = fault (Y=0 interior)
  Surface 5 = absorbing outer boundaries (4 sides + bottom)
  Volume  1 = bulk
"""

import argparse
import sys
from collections import defaultdict


def build_mesh(dx, lx, ly, lz):
    """Build mesh data: vertices, tets, boundary triangles."""
    nx, ny, nz = int(round(lx / dx)), int(round(ly / dx)), int(round(lz / dx))
    if ny % 2 != 0:
        raise ValueError(f"ly/dx={ny} must be even so Y=0 is a vertex plane")

    def vid(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k

    vertices = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                vertices.append((-lx / 2 + i * dx, -ly / 2 + j * dx, -k * dx))

    def hex_corner(i, j, k, lv):
        di, dj, dk = lv & 1, (lv >> 1) & 1, (lv >> 2) & 1
        return vid(i + di, j + dj, k + dk)

    # 6-tet split with diagonal v0–v7 (true opposite corners).
    # Cube vertex labels with bit encoding (x=bit0, y=bit1, z=bit2):
    #   0=(0,0,0)  1=(1,0,0)  2=(0,1,0)  3=(1,1,0)
    #   4=(0,0,1)  5=(1,0,1)  6=(0,1,1)  7=(1,1,1)
    pattern_a = [
        (0, 1, 3, 7),
        (0, 3, 2, 7),
        (0, 2, 6, 7),
        (0, 6, 4, 7),
        (0, 4, 5, 7),
        (0, 5, 1, 7),
    ]
    swap_y = {0: 2, 1: 3, 2: 0, 3: 1, 4: 6, 5: 7, 6: 4, 7: 5}
    pattern_b = [tuple(swap_y[v] for v in tet) for tet in pattern_a]

    tets = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                yc = -ly / 2 + (j + 0.5) * dx
                pat = pattern_a if yc > 0 else pattern_b
                for tet_local in pat:
                    tets.append(tuple(hex_corner(i, j, k, lv) for lv in tet_local))

    def reorient(t):
        a, b, c, d = (vertices[v - 1] for v in t)
        e1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        e2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        e3 = (d[0] - a[0], d[1] - a[1], d[2] - a[2])
        cross = (e2[1] * e3[2] - e2[2] * e3[1],
                 e2[2] * e3[0] - e2[0] * e3[2],
                 e2[0] * e3[1] - e2[1] * e3[0])
        det = e1[0] * cross[0] + e1[1] * cross[1] + e1[2] * cross[2]
        return (t[0], t[1], t[3], t[2]) if det < 0 else t

    tets = [reorient(t) for t in tets]

    # Derive boundary triangles from tet topology.
    def tet_faces(t):
        a, b, c, d = t
        return [(b, c, d), (a, d, c), (a, b, d), (a, c, b)]

    face_map = defaultdict(list)
    for ti, t in enumerate(tets):
        for f in tet_faces(t):
            face_map[tuple(sorted(f))].append((ti, f))

    bdr_triangles = []
    eps = 1e-6
    for key, lst in face_map.items():
        cx = sum(vertices[v - 1][0] for v in key) / 3
        cy = sum(vertices[v - 1][1] for v in key) / 3
        cz = sum(vertices[v - 1][2] for v in key) / 3
        if len(lst) == 1:
            ti, oriented = lst[0]
            if abs(cz) < eps:                attr = 1   # free surface
            elif abs(cz + lz) < eps:         attr = 5   # bottom
            elif abs(cx - lx / 2) < eps:     attr = 5   # +X
            elif abs(cx + lx / 2) < eps:     attr = 5   # -X
            elif abs(cy - ly / 2) < eps:     attr = 5   # +Y
            elif abs(cy + ly / 2) < eps:     attr = 5   # -Y
            else:                            continue
            bdr_triangles.append((oriented[0], oriented[1], oriented[2], attr))
        elif len(lst) == 2:
            if abs(cy) < eps:
                ti, oriented = lst[0]
                bdr_triangles.append((oriented[0], oriented[1], oriented[2], 3))

    return vertices, tets, bdr_triangles, (nx, ny, nz)


def write_msh(path, vertices, tets, bdr_triangles):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n4\n")
        f.write('2 1 "free_surface"\n')
        f.write('2 3 "fault"\n')
        f.write('2 5 "absorbing"\n')
        f.write('3 1 "bulk"\n')
        f.write("$EndPhysicalNames\n")
        f.write("$Nodes\n")
        f.write(f"{len(vertices)}\n")
        for n, (x, y, z) in enumerate(vertices, start=1):
            f.write(f"{n} {x} {y} {z}\n")
        f.write("$EndNodes\n")
        f.write("$Elements\n")
        f.write(f"{len(bdr_triangles) + len(tets)}\n")
        eid = 1
        for (a, b, c, attr) in bdr_triangles:
            f.write(f"{eid} 2 2 {attr} {attr} {a} {b} {c}\n")
            eid += 1
        for (a, b, c, d) in tets:
            f.write(f"{eid} 4 2 1 1 {a} {b} {c} {d}\n")
            eid += 1
        f.write("$EndElements\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dx", type=float, default=200.0,
                   help="Cell size in meters (default: 200)")
    p.add_argument("--lx", type=float, default=16000.0,
                   help="Domain length along strike (X) in meters (default: 16000)")
    p.add_argument("--ly", type=float, default=4000.0,
                   help="Domain length fault-normal (Y) in meters (default: 4000)")
    p.add_argument("--lz", type=float, default=16000.0,
                   help="Domain depth (Z, downward) in meters (default: 16000)")
    p.add_argument("--out", type=str, default="tpv104_symmirror_200m.msh",
                   help="Output .msh path")
    args = p.parse_args()

    vertices, tets, bdr_tris, (nx, ny, nz) = build_mesh(
        args.dx, args.lx, args.ly, args.lz)
    write_msh(args.out, vertices, tets, bdr_tris)

    n_fault = sum(1 for t in bdr_tris if t[3] == 3)
    n_free = sum(1 for t in bdr_tris if t[3] == 1)
    n_abs = sum(1 for t in bdr_tris if t[3] == 5)
    print(f"Wrote {args.out}")
    print(f"  domain  : [-{args.lx/2:.0f}, {args.lx/2:.0f}] x "
          f"[-{args.ly/2:.0f}, {args.ly/2:.0f}] x "
          f"[-{args.lz:.0f}, 0] m")
    print(f"  cells   : {nx} x {ny} x {nz}  (dx = {args.dx:.0f} m)")
    print(f"  nodes   : {len(vertices)}")
    print(f"  tets    : {len(tets)}")
    print(f"  bdr     : free={n_free}  fault={n_fault}  absorb={n_abs}")


if __name__ == "__main__":
    main()
