#!/usr/bin/env python3
"""
ts_to_stl.py — Convert a GoCAD TSurf (.ts) fault surface to STL,
clipping all material above z = 0 (free surface).

Clipping is done per-triangle with linear interpolation at z = 0:
    - all three vertices below z = 0   -> keep triangle as-is
    - all three vertices above z = 0   -> discard
    - one vertex below, two above      -> emit one triangle (b, a', c')
                                          where a', c' are the z=0
                                          intersections on edges b-a, b-c
    - two vertices below, one above    -> emit two triangles tiling the
                                          quadrilateral {b, c, c', a'}

Triangle orientation is preserved (the resulting normal points the same
way as the original triangle's normal).

Vertices already at z = 0 are treated as on-or-below.

Optional top-edge snap (`--top-offset OFFSET`):
    After clipping, snap every output vertex with z > -OFFSET down to
    z = -OFFSET.  This produces a clean horizontal top edge at
    z = -OFFSET and keeps the fault strictly interior to a box whose
    top face sits at z = 0, which gmsh meshing requires to enclose a
    closed volume.  The original fault surface near z = 0 is jagged
    (clipped vertices at z = 0 mixed with un-clipped originals at
    z ~ -30 m), so the snap also flattens that zigzag.

Pass-through mode (`--no-clip`):
    Skip clipping entirely and emit every input triangle as-is, including
    triangles with z > 0.  Useful for inspecting the raw surface (or
    re-triangulating the cap externally, e.g. in MeshLab) before any
    z = 0 manipulation.

Usage:
    python ts_to_stl.py INPUT.ts OUTPUT.stl
        [--name solidname] [--binary] [--top-offset OFFSET] [--no-clip]
"""

import argparse
import struct
import sys
from pathlib import Path

EPS = 1.0e-9


def parse_ts(path):
    """Parse a GoCAD TSurf file. Returns (verts, tris).

    verts is a dict: vertex_id -> (x, y, z).
    tris is a list of (id_a, id_b, id_c) tuples.
    """
    verts = {}
    tris = []
    with open(path, "r") as f:
        for line in f:
            tok = line.split()
            if not tok:
                continue
            head = tok[0]
            if head == "VRTX" or head == "PVRTX":
                # VRTX <id> <x> <y> <z>           (PVRTX has extra fields after z)
                vid = int(tok[1])
                x = float(tok[2])
                y = float(tok[3])
                z = float(tok[4])
                verts[vid] = (x, y, z)
            elif head == "ATOM":
                # ATOM <new_id> <ref_id>: alias of an existing vertex
                new_id = int(tok[1])
                ref_id = int(tok[2])
                if ref_id in verts:
                    verts[new_id] = verts[ref_id]
            elif head == "TRGL":
                a = int(tok[1])
                b = int(tok[2])
                c = int(tok[3])
                tris.append((a, b, c))
    return verts, tris


def lerp_to_z0(p_below, p_above):
    """Linearly interpolate between p_below (z<=0) and p_above (z>0) to find
    the point where the segment crosses z = 0."""
    zb = p_below[2]
    za = p_above[2]
    # Segment: P(t) = p_below + t * (p_above - p_below), 0 <= t <= 1
    # P(t).z = 0  =>  t = -zb / (za - zb)
    denom = za - zb
    if abs(denom) < EPS:
        # Degenerate: both at ~0; just return one of them at z = 0.
        return (p_below[0], p_below[1], 0.0)
    t = -zb / denom
    x = p_below[0] + t * (p_above[0] - p_below[0])
    y = p_below[1] + t * (p_above[1] - p_below[1])
    return (x, y, 0.0)


def clip_triangle_at_z0(p0, p1, p2):
    """Clip a triangle to the z <= 0 half-space.

    Returns a list of triangles (each = 3 (x,y,z) tuples) lying at or below
    z = 0. The original triangle's vertex order is preserved, so the
    returned triangle(s) carry the same outward normal direction as the
    input.
    """
    pts = [p0, p1, p2]
    z = [p[2] for p in pts]
    # "Below" = at or below the clipping plane.
    below = [zi <= EPS for zi in z]
    n_below = sum(below)

    if n_below == 3:
        return [(p0, p1, p2)]
    if n_below == 0:
        return []

    if n_below == 1:
        # One below, two above: clip to a single triangle.
        # Find the one below and label the other two in CCW order.
        bi = below.index(True)
        ai = (bi + 1) % 3
        ci = (bi + 2) % 3
        b = pts[bi]
        a = pts[ai]
        c = pts[ci]
        a_clip = lerp_to_z0(b, a)  # on edge b->a
        c_clip = lerp_to_z0(b, c)  # on edge b->c -- but we want b->...->c order
        # Triangle (b, a_clip, c_clip) preserves the original orientation
        # (b, a, c) is the same cyclic order as (bi, ai, ci) of (p0, p1, p2).
        return [(b, a_clip, c_clip)]

    # n_below == 2: two below, one above; emit two triangles.
    ai = below.index(False)        # the above vertex
    bi = (ai + 1) % 3
    ci = (ai + 2) % 3
    a = pts[ai]
    b = pts[bi]
    c = pts[ci]
    ab = lerp_to_z0(b, a)  # on edge b->a, but we want a->b order; lerp is symmetric for the point
    ac = lerp_to_z0(c, a)
    # Original ordering (a, b, c) cyclic; we drop a and replace with ab, ac.
    # Resulting quad is (ab, b, c, ac); split into (ab, b, c) and (ab, c, ac).
    return [(ab, b, c), (ab, c, ac)]


def write_stl_ascii(triangles, path, name="fault"):
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for p0, p1, p2 in triangles:
            ux = p1[0] - p0[0]
            uy = p1[1] - p0[1]
            uz = p1[2] - p0[2]
            vx = p2[0] - p0[0]
            vy = p2[1] - p0[1]
            vz = p2[2] - p0[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            mag = (nx * nx + ny * ny + nz * nz) ** 0.5
            if mag > 0.0:
                nx /= mag
                ny /= mag
                nz /= mag
            f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n")
            f.write("    outer loop\n")
            for p in (p0, p1, p2):
                f.write(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")


def write_stl_binary(triangles, path, name="fault"):
    with open(path, "wb") as f:
        header = name.encode("ascii", errors="replace")[:80].ljust(80, b"\x00")
        f.write(header)
        f.write(struct.pack("<I", len(triangles)))
        for p0, p1, p2 in triangles:
            ux = p1[0] - p0[0]
            uy = p1[1] - p0[1]
            uz = p1[2] - p0[2]
            vx = p2[0] - p0[0]
            vy = p2[1] - p0[1]
            vz = p2[2] - p0[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            mag = (nx * nx + ny * ny + nz * nz) ** 0.5
            if mag > 0.0:
                nx /= mag
                ny /= mag
                nz /= mag
            f.write(struct.pack("<3f", nx, ny, nz))
            for p in (p0, p1, p2):
                f.write(struct.pack("<3f", p[0], p[1], p[2]))
            f.write(struct.pack("<H", 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="input GoCAD .ts file")
    ap.add_argument("output", type=Path, help="output .stl file")
    ap.add_argument("--name", default="fault", help="solid name (default: fault)")
    ap.add_argument("--binary", action="store_true", help="write binary STL (default: ASCII)")
    ap.add_argument("--top-offset", type=float, default=0.0,
                    help="after clipping, snap every output vertex with "
                         "z > -OFFSET down to z = -OFFSET. Default 0 = no "
                         "snap (top edge stays at z = 0).")
    ap.add_argument("--no-clip", action="store_true",
                    help="skip the z = 0 clipping pass and emit every input "
                         "triangle as-is. Use this to inspect the raw "
                         "surface (including z > 0) in MeshLab/ParaView.")
    args = ap.parse_args()

    if not args.input.is_file():
        print(f"error: input not found: {args.input}", file=sys.stderr)
        return 1

    verts, tris = parse_ts(args.input)
    print(f"Parsed {args.input.name}: {len(verts)} vertices, {len(tris)} triangles")

    n_above = sum(1 for v in verts.values() if v[2] > EPS)
    print(f"Vertices with z > 0: {n_above}")

    if args.no_clip:
        out_tris = []
        n_missing = 0
        for ia, ib, ic in tris:
            if ia not in verts or ib not in verts or ic not in verts:
                n_missing += 1
                continue
            out_tris.append((verts[ia], verts[ib], verts[ic]))
        print(f"--no-clip: emitting all input triangles as-is "
              f"(missing-vertex skips: {n_missing})")
        print(f"Output triangles     : {len(out_tris)}")
        if args.top_offset > 0.0:
            print("warning: --top-offset is ignored when --no-clip is set",
                  file=sys.stderr)
            args.top_offset = 0.0
    else:
        out_tris = []
        n_kept_full = 0
        n_clipped = 0
        n_dropped = 0
        for ia, ib, ic in tris:
            if ia not in verts or ib not in verts or ic not in verts:
                continue
            result = clip_triangle_at_z0(verts[ia], verts[ib], verts[ic])
            if not result:
                n_dropped += 1
            elif len(result) == 1 and all(abs(verts[i][2]) <= EPS or verts[i][2] < 0.0 for i in (ia, ib, ic)):
                n_kept_full += 1
            else:
                n_clipped += 1
            out_tris.extend(result)

        print(f"Triangles kept whole : {n_kept_full}")
        print(f"Triangles clipped    : {n_clipped}")
        print(f"Triangles dropped    : {n_dropped}")
        print(f"Output triangles     : {len(out_tris)}")

    if not out_tris:
        print("error: no triangles to write", file=sys.stderr)
        return 1

    if args.top_offset > 0.0:
        z_top = -float(args.top_offset)
        n_snapped = 0
        snapped_tris = []
        for tri in out_tris:
            new_tri = []
            for p in tri:
                if p[2] > z_top:
                    new_tri.append((p[0], p[1], z_top))
                    n_snapped += 1
                else:
                    new_tri.append(p)
            snapped_tris.append(tuple(new_tri))
        # Drop triangles that became degenerate (all three vertices at z_top
        # AND collinear, or exactly co-located).  We keep all triangles whose
        # area is strictly positive.
        kept = []
        for p0, p1, p2 in snapped_tris:
            ux = p1[0] - p0[0]; uy = p1[1] - p0[1]; uz = p1[2] - p0[2]
            vx = p2[0] - p0[0]; vy = p2[1] - p0[1]; vz = p2[2] - p0[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            if (nx * nx + ny * ny + nz * nz) > 1.0e-6:
                kept.append((p0, p1, p2))
        n_dropped_post = len(snapped_tris) - len(kept)
        out_tris = kept
        print(f"Top-edge snap to z = {z_top:.2f}: snapped {n_snapped} vertices, "
              f"dropped {n_dropped_post} degenerate triangles after snap, "
              f"output triangles now {len(out_tris)}")

    xs = []
    ys = []
    zs = []
    for tri in out_tris:
        for p in tri:
            xs.append(p[0])
            ys.append(p[1])
            zs.append(p[2])
    print(f"Output bbox: X=[{min(xs):.2f}, {max(xs):.2f}]")
    print(f"             Y=[{min(ys):.2f}, {max(ys):.2f}]")
    print(f"             Z=[{min(zs):.2f}, {max(zs):.2f}]")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.binary:
        write_stl_binary(out_tris, args.output, name=args.name)
    else:
        write_stl_ascii(out_tris, args.output, name=args.name)
    size_mb = args.output.stat().st_size / (1024 * 1024)
    print(f"Wrote {args.output} ({size_mb:.2f} MB, {'binary' if args.binary else 'ascii'} STL)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
