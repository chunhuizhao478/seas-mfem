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


def _make_symmetric_widths(d_fine, inner_half, outer_half, ratio):
    """Build a y-mirror-symmetric cell-width sequence about the origin.

    Inner zone: |coord| <= inner_half is uniform with width d_fine.
    Outer zone: cells geometrically grow by `ratio` each step until the
    cumulative outer extent reaches (outer_half - inner_half).  The final
    outer cell is clipped if needed so the total cumulative half-extent
    equals exactly `outer_half`.

    Returns a list of widths in order from -outer_half to +outer_half.
    The list is symmetric: widths[j] == widths[len-1-j], guaranteeing
    that y-mirror partner cells (j, len-1-j) have identical geometry.

    Edge cases:
      - outer_half <= inner_half: no grading; pure uniform inner.
      - ratio == 1.0: outer cells are uniform too at width d_fine.
    """
    if outer_half <= 0 or inner_half < 0:
        raise ValueError(
            f"_make_symmetric_widths: outer_half={outer_half} must be > 0 and "
            f"inner_half={inner_half} must be >= 0")
    if inner_half > outer_half:
        raise ValueError(
            f"_make_symmetric_widths: inner_half={inner_half} > "
            f"outer_half={outer_half}")
    if d_fine <= 0:
        raise ValueError(f"_make_symmetric_widths: d_fine={d_fine} must be > 0")
    if ratio <= 0:
        raise ValueError(f"_make_symmetric_widths: ratio={ratio} must be > 0")

    # Inner zone: 2*N_inner cells of width d_fine, where 2*N_inner*d_fine
    # = 2*inner_half exactly (we require inner_half to be a multiple of
    # d_fine for clean placement of the y=0 vertex).
    n_inner_per_side = int(round(inner_half / d_fine))
    if abs(n_inner_per_side * d_fine - inner_half) > 1e-9 * max(d_fine, 1.0):
        raise ValueError(
            f"_make_symmetric_widths: inner_half={inner_half} must be an "
            f"integer multiple of d_fine={d_fine} (got {inner_half/d_fine} cells)")
    inner = [d_fine] * (2 * n_inner_per_side)

    # Outer zone: cells of width d_fine*ratio, d_fine*ratio^2, ..., up to
    # (outer_half - inner_half).  Final cell clipped to fit exactly.
    outer_span = outer_half - inner_half
    outer_per_side = []
    span = 0.0
    w = d_fine
    while span < outer_span - 1e-9 * max(d_fine, 1.0):
        w_next = w * ratio
        remaining = outer_span - span
        if w_next >= remaining:
            outer_per_side.append(remaining)
            span = outer_span
            break
        outer_per_side.append(w_next)
        span += w_next
        w = w_next

    # Assemble: reverse(outer) + inner + outer (mirror about 0).
    return list(reversed(outer_per_side)) + inner + outer_per_side


def _make_asymmetric_widths_z(d_fine, inner_depth, outer_depth, ratio):
    """Build a cell-width sequence for the z-axis (depth direction).

    Free surface at z=0 → no grading on top.  Fine 200m for z in
    [-inner_depth, 0] (the rupture-active zone).  Below z=-inner_depth,
    cells grow geometrically by `ratio` until the cumulative depth
    reaches `outer_depth` (final cell clipped if needed).

    Returns widths in order from z=0 (top) to z=-outer_depth (bottom),
    so z-coordinate of vertex k is -sum(widths[:k]).
    """
    if outer_depth <= 0 or inner_depth < 0:
        raise ValueError(
            f"_make_asymmetric_widths_z: outer_depth={outer_depth} must be > 0 "
            f"and inner_depth={inner_depth} must be >= 0")
    if inner_depth > outer_depth:
        raise ValueError(
            f"_make_asymmetric_widths_z: inner_depth={inner_depth} > "
            f"outer_depth={outer_depth}")
    if d_fine <= 0:
        raise ValueError(f"_make_asymmetric_widths_z: d_fine={d_fine} must be > 0")

    n_inner = int(round(inner_depth / d_fine))
    if abs(n_inner * d_fine - inner_depth) > 1e-9 * max(d_fine, 1.0):
        raise ValueError(
            f"_make_asymmetric_widths_z: inner_depth={inner_depth} must be an "
            f"integer multiple of d_fine={d_fine}")
    inner = [d_fine] * n_inner

    outer_span = outer_depth - inner_depth
    outer = []
    span = 0.0
    w = d_fine
    while span < outer_span - 1e-9 * max(d_fine, 1.0):
        w_next = w * ratio
        remaining = outer_span - span
        if w_next >= remaining:
            outer.append(remaining)
            span = outer_span
            break
        outer.append(w_next)
        span += w_next
        w = w_next

    return inner + outer


def build_mesh(dx, lx, ly, lz,
               x_inner_half=None, y_inner_half=None, z_inner_depth=None,
               x_ratio=1.0, y_ratio=1.0, z_ratio=1.0):
    """Build mesh data: vertices, tets, boundary triangles.

    If `*_inner_half` / `z_inner_depth` are None, defaults to lx/2, ly/2,
    lz, respectively (the entire domain is uniform — backward compatible
    with pre-grading callers).

    `*_ratio` is the geometric grading ratio outside the inner zone.
    Ratio == 1.0 reproduces uniform spacing in the outer zone.
    """
    if x_inner_half is None: x_inner_half = lx / 2
    if y_inner_half is None: y_inner_half = ly / 2
    if z_inner_depth is None: z_inner_depth = lz

    x_widths = _make_symmetric_widths(dx, x_inner_half, lx / 2, x_ratio)
    y_widths = _make_symmetric_widths(dx, y_inner_half, ly / 2, y_ratio)
    z_widths = _make_asymmetric_widths_z(dx, z_inner_depth, lz, z_ratio)

    nx = len(x_widths)
    ny = len(y_widths)
    nz = len(z_widths)
    if ny % 2 != 0:
        raise ValueError(
            f"ly/dx-equivalent ny={ny} must be even so Y=0 is a vertex plane "
            f"(check y_inner_half={y_inner_half} is multiple of dx={dx})")

    # Cumulative vertex coordinates.
    xs = [-lx / 2]
    for w in x_widths:
        xs.append(xs[-1] + w)
    ys = [-ly / 2]
    for w in y_widths:
        ys.append(ys[-1] + w)
    zs = [0.0]                         # free surface at z = 0
    for w in z_widths:
        zs.append(zs[-1] - w)          # depth increases downward

    def vid(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k

    vertices = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                vertices.append((xs[i], ys[j], zs[k]))

    def hex_corner(i, j, k, lv):
        di, dj, dk = lv & 1, (lv >> 1) & 1, (lv >> 2) & 1
        return vid(i + di, j + dj, k + dk)

    # 6-tet split with diagonal v0–v7 (true opposite corners).
    # Cube vertex labels with bit encoding (x=bit0, y=bit1, z=bit2):
    #   0=(0,0,0)  1=(1,0,0)  2=(0,1,0)  3=(1,1,0)
    #   4=(0,0,1)  5=(1,0,1)  6=(0,1,1)  7=(1,1,1)
    #
    # Note (R-104 audit, 2026-04-25): pattern_b = swap_y(pattern_a) is
    # *not* tuple-mirror at the +y / -y halves after reorient.  This is
    # mathematically unavoidable: y-mirror reverses orientation, so
    # exactly one of {pattern_a, swap_y(pattern_a)} produces det<0 tets
    # that reorient must permute (v3 <-> v4), breaking pure tuple
    # mirror.  The vertex SETS are still y-mirror; only the per-tet
    # tuple ORDER differs.  Confirming this defect is the cause of the
    # np=1 slip_dip baseline requires either (i) a self-mirror-symmetric
    # tet split (e.g. 12-tet Kuhn) or (ii) probing CalcOrtho normals
    # on non-fault interior faces for mirror consistency.
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
                # Cell-center y from cumulative vertex coordinates (handles
                # graded spacing where (j+0.5)*dx is wrong).
                yc = 0.5 * (ys[j] + ys[j + 1])
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


def factor_2d(np):
    """Find (nx_ranks, nz_ranks) closest to sqrt(np) with nx_ranks * nz_ranks == np."""
    best = (1, np)
    for a in range(1, int(np**0.5) + 1):
        if np % a == 0:
            best = (a, np // a)
    return best


def write_partition(path, nx, ny, nz, np_target):
    """Emit a partitioning sidecar file: pure Cartesian X+Z tiling, full
    Y per rank.

    For a y-mirror-symmetric problem (TPV104 fault at y=0), each rank
    owns the full y-extent of its (X, Z) tile, so the y=0 plane is
    never cut by a partition boundary.
    """
    nx_ranks, nz_ranks = factor_2d(np_target)
    parts = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                rx = (i * nx_ranks) // nx
                rz = (k * nz_ranks) // nz
                rank = rx + rz * nx_ranks
                for _ in range(6):
                    parts.append(rank)

    ne = len(parts)
    with open(path, "w") as f:
        f.write(f"np {np_target}\n")
        f.write(f"ne {ne}\n")
        for r in parts:
            f.write(f"{r}\n")
    return nx_ranks, nz_ranks


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
    p.add_argument("--emit-partition", type=int, nargs="+", metavar="NP",
                   help="Also emit <out>.np<NP>.partition for each NP")
    # Per-axis grading: defaults reproduce the uniform mesh.  When set,
    # the inner zone (|x| <= x_inner_half, |y| <= y_inner_half, depth <=
    # z_inner_depth) is uniform at `dx`; the outer zone uses geometric
    # grading with the corresponding ratio.  inner_half / inner_depth
    # MUST be integer multiples of dx.
    p.add_argument("--x-inner-half", type=float, default=None,
                   help="Inner uniform-dx half-extent along strike (X) "
                        "[m].  Default: lx/2 (entire domain uniform).")
    p.add_argument("--y-inner-half", type=float, default=None,
                   help="Inner uniform-dx half-extent fault-normal (Y) "
                        "[m].  Default: ly/2 (entire domain uniform).")
    p.add_argument("--z-inner-depth", type=float, default=None,
                   help="Inner uniform-dx depth from free surface [m]. "
                        "Default: lz (entire depth uniform).")
    p.add_argument("--x-ratio", type=float, default=1.0,
                   help="Geometric grading ratio outside x-inner zone. "
                        "1.0 = uniform.  (default: 1.0)")
    p.add_argument("--y-ratio", type=float, default=1.0,
                   help="Geometric grading ratio outside y-inner zone. "
                        "(default: 1.0)")
    p.add_argument("--z-ratio", type=float, default=1.0,
                   help="Geometric grading ratio below z-inner zone. "
                        "(default: 1.0)")
    args = p.parse_args()

    vertices, tets, bdr_tris, (nx, ny, nz) = build_mesh(
        args.dx, args.lx, args.ly, args.lz,
        x_inner_half=args.x_inner_half,
        y_inner_half=args.y_inner_half,
        z_inner_depth=args.z_inner_depth,
        x_ratio=args.x_ratio,
        y_ratio=args.y_ratio,
        z_ratio=args.z_ratio)
    write_msh(args.out, vertices, tets, bdr_tris)

    if args.emit_partition:
        for np_target in args.emit_partition:
            part_path = f"{args.out}.np{np_target}.partition"
            nx_ranks, nz_ranks = write_partition(part_path, nx, ny, nz,
                                                 np_target)
            print(f"  partition np={np_target}: "
                  f"nx_ranks={nx_ranks} nz_ranks={nz_ranks} "
                  f"(pure Cartesian X+Z, full y-extent per rank) "
                  f"-> {part_path}")

    n_fault = sum(1 for t in bdr_tris if t[3] == 3)
    n_free = sum(1 for t in bdr_tris if t[3] == 1)
    n_abs = sum(1 for t in bdr_tris if t[3] == 5)
    print(f"Wrote {args.out}")
    print(f"  domain  : [-{args.lx/2:.0f}, {args.lx/2:.0f}] x "
          f"[-{args.ly/2:.0f}, {args.ly/2:.0f}] x "
          f"[-{args.lz:.0f}, 0] m")
    print(f"  cells   : {nx} x {ny} x {nz}  (dx_fine = {args.dx:.0f} m)")
    if args.x_inner_half is not None or args.y_inner_half is not None \
       or args.z_inner_depth is not None:
        print(f"  inner   : x_half={args.x_inner_half} m  "
              f"y_half={args.y_inner_half} m  z_depth={args.z_inner_depth} m")
        print(f"  ratios  : x={args.x_ratio}  y={args.y_ratio}  "
              f"z={args.z_ratio}")
    print(f"  nodes   : {len(vertices)}")
    print(f"  tets    : {len(tets)}")
    print(f"  bdr     : free={n_free}  fault={n_fault}  absorb={n_abs}")


if __name__ == "__main__":
    main()
