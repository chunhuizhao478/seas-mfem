#!/usr/bin/env python3
"""retriangulate_fault_rv.py — force every fault face below an r_v ceiling by a
conforming local 1->4 tetra split, for faces that node smoothing cannot fix.

Eq.(18):  r_v = max(h_A/h_B, h_B/h_A)  per fault triangle, h = apex perpendicular
distance to the (fixed) fault plane.  When a face exceeds the ceiling and is
inversion-blocked (no apex can move), we instead split the TALLER of its two
tets:

    tall tet A = (v0,v1,v2, a)   with the fault triangle f=(v0,v1,v2) and the
    tall apex a (|h_a| > |h_b|).  Insert a node p on the segment centroid(f)->a
    at height |h_b| (so it is strictly interior to A), and replace A by the four
    tets that connect p to each face of A:
        (v0,v1,v2, p)   <- new fault-adjacent tet, height |h_b|  => r_v = 1
        (v0,v1, a, p) , (v1,v2, a, p) , (v0,v2, a, p)

All four OUTER faces of A are preserved, so the mesh stays conforming (neighbours
untouched); the fault triangle f and every fault/boundary node are unchanged.
For a PLANAR fault each tet borders at most one fault face, so the tall tets of
the over-ceiling faces are disjoint and the splits are independent.

Intended input: the smoothed mesh (few over-ceiling faces left).  +1 node and
+3 tets per fixed face.

Usage:
  python retriangulate_fault_rv.py --in IN.msh --out OUT.msh \
      [--fault-phys 101] [--rock-phys 1] [--rv-ceiling 1.5] [--target-rv 1.0]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import smooth_fault_volume_ratio as S   # read_gmsh22, signed_volumes, rv_*, tet_etas


def _face_to_tets(tets, fault_tris):
    """Per fault triangle, the two (tet_index, apex_node) it bounds.

    Returns a list aligned with fault_tris rows; raises if any face is not
    bounded by exactly two tets."""
    fault_keys = {frozenset(int(v) for v in t) for t in fault_tris}
    fmap = defaultdict(list)
    faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for ti in range(tets.shape[0]):
        t = tets[ti]
        nodes = (int(t[0]), int(t[1]), int(t[2]), int(t[3]))
        for fa in faces:
            key = frozenset((nodes[fa[0]], nodes[fa[1]], nodes[fa[2]]))
            if key in fault_keys:
                apex = (set(nodes) - key).pop()
                fmap[key].append((ti, apex))
    out = []
    for tri in fault_tris:
        key = frozenset(int(v) for v in tri)
        pair = fmap.get(key, [])
        if len(pair) != 2:
            raise ValueError("a fault face is not bounded by exactly 2 tets; "
                             "fault not cleanly embedded")
        out.append(pair)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, required=True)
    ap.add_argument("--out", dest="out", type=Path, required=True)
    ap.add_argument("--fault-phys", type=int, default=101)
    ap.add_argument("--rock-phys", type=int, default=1)
    ap.add_argument("--boundary-phys", type=str, default="102,103,104")
    ap.add_argument("--rv-ceiling", type=float, default=1.5,
                    help="split every fault face with r_v above this")
    ap.add_argument("--target-rv", type=float, default=1.0,
                    help="r_v the split should achieve on each fixed face "
                         "(1.0 = match the short side; must be < rv-ceiling)")
    args = ap.parse_args(argv)
    bphys = tuple(int(x) for x in args.boundary_phys.split(","))
    if not (1.0 <= args.target_rv < args.rv_ceiling):
        print("ERROR: need 1.0 <= --target-rv < --rv-ceiling", file=sys.stderr)
        return 1

    print(f"reading {args.inp} ...")
    g = S.read_gmsh22(args.inp, fault_phys=args.fault_phys,
                      rock_phys=args.rock_phys, boundary_phys=bphys)
    topo = S.build_fault_apex_map(g.coords, g.tets, g.fault_tris)
    rv0 = S.rv_values(g.coords, topo)
    over = np.where(rv0 > args.rv_ceiling)[0]
    print(f"  tets={g.tets.shape[0]:,}  fault_tris={g.fault_tris.shape[0]:,}  "
          f"faces > {args.rv_ceiling}: {over.size} (max r_v={rv0.max():.4f})")
    if over.size == 0:
        print("  nothing to do; copying input -> output")

    print("  building fault-face -> tet map ...")
    f2t = _face_to_tets(g.tets, g.fault_tris)

    coords = g.coords
    tets = g.tets
    sgn0 = np.sign(S.signed_volumes(coords, tets))
    new_nodes = []                  # list of (x,y,z) for appended nodes
    split_tet = {}                  # tet_idx -> list of 4 new tet rows (idx)
    n_collide = 0
    for fi in over:
        tri = g.fault_tris[fi]
        nrm, cen = topo.n_f[fi], topo.c_f[fi]
        (ta, aa), (tb, ab) = f2t[fi]
        ha = abs(float((coords[aa] - cen) @ nrm))
        hb = abs(float((coords[ab] - cen) @ nrm))
        # split the TALLER tet; target height = (short side height)*target_rv
        if ha >= hb:
            tall_t, apex, h_tall, h_short = ta, aa, ha, hb
        else:
            tall_t, apex, h_tall, h_short = tb, ab, hb, ha
        if tall_t in split_tet:
            n_collide += 1
            continue
        if h_tall <= 0:
            continue
        target_h = min(h_tall, args.target_rv * h_short)
        s = target_h / h_tall                      # in (0,1] -> p interior
        p = (1.0 - s) * cen + s * coords[apex]
        pid = len(coords) + len(new_nodes)         # new node index
        new_nodes.append(p)
        v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
        a = int(apex)
        # 4 sub-tets (orientation fixed below)
        sub = [[v0, v1, v2, pid], [v0, v1, a, pid],
               [v1, v2, a, pid], [v0, v2, a, pid]]
        split_tet[tall_t] = (sub, sgn0[tall_t])

    if n_collide:
        print(f"  WARNING: {n_collide} tall-tet collisions skipped (a tet was "
              f"the taller tet of two over-ceiling faces — unexpected on a "
              f"planar fault). Those faces are NOT fixed.")

    # assemble new coords + tets, fixing sub-tet orientation to match the parent
    extra = np.asarray(new_nodes, dtype=np.float64).reshape(-1, 3)
    coords_new = np.vstack([coords, extra]) if extra.size else coords.copy()

    new_tets = []
    for ti in range(tets.shape[0]):
        if ti in split_tet:
            sub, parent_sgn = split_tet[ti]
            for row in sub:
                r = list(row)
                v = S.signed_volumes(coords_new, np.array([r]))[0]
                if np.sign(v) != parent_sgn:
                    r[1], r[2] = r[2], r[1]
                new_tets.append(r)
        else:
            new_tets.append([int(x) for x in tets[ti]])
    new_tets = np.asarray(new_tets, dtype=np.int64)
    print(f"  split {len(split_tet)} tets -> +{extra.shape[0]} nodes, "
          f"+{3*len(split_tet)} tets (now {new_tets.shape[0]:,})")

    # --- write Gmsh 2.2: triangle element lines verbatim, tets regenerated ---
    tri_lines = [ln for ln in g.elem_lines if ln.split()[1] != "4"]
    with open(args.out, "w") as o:
        o.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        if g.phys_lines:
            o.write("$PhysicalNames\n")
            o.write(f"{len(g.phys_lines)}\n")
            o.write("\n".join(g.phys_lines) + "\n")
            o.write("$EndPhysicalNames\n")
        o.write("$Nodes\n")
        n_old = len(g.node_ids)
        o.write(f"{n_old + extra.shape[0]}\n")
        for k in range(n_old):
            o.write(g.node_lines[k] + "\n")
        for j in range(extra.shape[0]):
            x, y, z = extra[j]
            o.write(f"{n_old + j + 1} {x:.12g} {y:.12g} {z:.12g}\n")
        o.write("$EndNodes\n")
        o.write("$Elements\n")
        o.write(f"{len(tri_lines) + new_tets.shape[0]}\n")
        eid = 0
        id1 = {i: nid for i, nid in enumerate(g.node_ids)}   # idx(0)->gmsh id
        next_id = int(g.node_ids.max()) + 1
        new_id = {n_old + j: next_id + j for j in range(extra.shape[0])}

        def gid(idx):
            return id1[idx] if idx < n_old else new_id[idx]

        for ln in tri_lines:                  # triangles unchanged (verbatim)
            eid += 1
            p = ln.split()
            o.write(" ".join([str(eid)] + p[1:]) + "\n")
        for r in new_tets:
            eid += 1
            o.write(f"{eid} 4 2 {args.rock_phys} {args.rock_phys} "
                    f"{gid(r[0])} {gid(r[1])} {gid(r[2])} {gid(r[3])}\n")
        o.write("$EndElements\n")
    print(f"wrote {args.out}")

    # --- verify on the written mesh ---
    print("verifying ...")
    gv = S.read_gmsh22(args.out, fault_phys=args.fault_phys,
                       rock_phys=args.rock_phys, boundary_phys=bphys)
    tv = S.build_fault_apex_map(gv.coords, gv.tets, gv.fault_tris)
    s = S.rv_stats(gv.coords, tv)
    eta = S.tet_etas(gv.coords, gv.tets)
    inv = int((S.signed_volumes(gv.coords, gv.tets) <= 0).sum())
    ff = np.array_equal(gv.coords[g.fault_node_idx], g.coords[g.fault_node_idx])
    print(f"  r_v max={s['max']:.4f}  > {args.rv_ceiling}: {s['n_gt15'] if args.rv_ceiling==1.5 else int((S.rv_values(gv.coords,tv)>args.rv_ceiling).sum())}"
          f"  mean={s['mean']:.4f}")
    print(f"  eta_min={eta.min():.4f}  eta_med={np.median(eta):.4f}  "
          f"non-positive tets={inv}  fault_nodes_frozen={ff}")
    n_bad = int((S.rv_values(gv.coords, tv) > args.rv_ceiling).sum())
    if n_bad == 0 and inv == 0 and ff:
        print(f"  SUCCESS: every fault face r_v <= {args.rv_ceiling}")
        return 0
    print(f"  INCOMPLETE: {n_bad} face(s) still over, inverted={inv}",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
