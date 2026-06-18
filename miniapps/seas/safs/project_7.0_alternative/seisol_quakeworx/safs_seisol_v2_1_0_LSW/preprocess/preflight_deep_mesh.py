#!/usr/bin/env python3
"""Phase 0 pre-flight for porting the LSW SeisSol case to the new multi-strand
mesh ``safv4_deep_500m.msh`` (PLAN_safv4_deep_mesh_port_2026-06-15.md).

Read-only.  Parses the gmsh v2.2 mesh directly (the file has NO $PhysicalNames
block) and runs the compatibility gates required before any retag / pumgen work:

  1. Tag inventory      -- triangles == {101,102,103,201,202,203}, tets == {1}
  2. Bounding box       -- of all referenced nodes
  3. CVM coverage       -- fraction of nodes outside the ASAGI grid / above z=0
                           (informational; clamping is accepted, not a gate)
  4. Hypocenter on fault-- min point-to-triangle distance from the nucleation
                           centre to any fault facet must be < 600 m   [GATE]
  5. Fault not in xy    -- fault-triangle normals must be near-horizontal
                           (faults near-vertical); large horizontal patches warn

Exit code 0 on PASS, non-zero on any FAIL.

Usage:
  preflight_deep_mesh.py [path/to/safv4_deep_500m.msh]
Default mesh path is the project_7.0_preferred result relative to this folder.
"""
import os
import sys
from collections import Counter

import numpy as np

# --- Fixed problem constants (do NOT edit; they define the reused physics) ----

# Nucleation Gaussian centre = hypocenter, from safs_fault.yaml [Tnuc_s] and
# README_port.md.  (x East, y North, z Up; UTM 11N.)
HYPOCENTER = np.array([606971.0, 3707270.0, -4965.62])

# Fault strand physical tags in the new preferred-project mesh
# (RUNLOG_safv4_remesh_phases0-5_2026-06-12.md line 264).
FAULT_TAGS = {101, 102, 103}
STRAND_NAME = {101: "San Andreas", 102: "Banning", 103: "Garnet Hill"}
EXPECTED_TRI_TAGS = {101, 102, 103, 201, 202, 203}
EXPECTED_TET_TAGS = {1}

# CVM ASAGI grid extent (read from safs_material_cvm.nc on 2026-06-15).
GRID_X = (303000.0, 697500.0)
GRID_Y = (3612000.0, 3903000.0)
GRID_Z = (-45000.0, 0.0)

# Gate: hypocenter must sit within ~one 500 m element of a fault facet.
HYPO_GATE_M = 600.0

# this script lives in safs_seisol_v2_1_0_LSW/preprocess/ -> 4 levels up to safs/
DEFAULT_MESH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "..", "project_7.0_preferred", "meshing", "results",
    "safv4_deep_500m.msh",
)

GMSH_TRI = 2   # 3-node triangle
GMSH_TET = 4   # 4-node tetrahedron


def parse_msh(path):
    """Stream a gmsh v2.2 ASCII mesh.

    Returns
      coords     : (Nn, 3) float64 node coordinates, row i = node index i
      tri_tags   : (Nt,)   int   first physical tag of each triangle
      tri_nodes  : (Nt, 3) int   node indices (into coords) of each triangle
      tet_tags   : Counter       {physical tag: count} for tetrahedra
    """
    coords_list = []
    id_to_idx = {}
    tri_tags = []
    tri_nodes = []
    tet_tags = Counter()

    with open(path) as fh:
        line = fh.readline()
        while line:
            s = line.strip()
            if s == "$Nodes":
                nn = int(fh.readline())
                for _ in range(nn):
                    p = fh.readline().split()
                    nid = int(p[0])
                    id_to_idx[nid] = len(coords_list)
                    coords_list.append((float(p[1]), float(p[2]), float(p[3])))
                end = fh.readline().strip()
                if end != "$EndNodes":
                    raise ValueError(f"expected $EndNodes, got {end!r}")
            elif s == "$Elements":
                ne = int(fh.readline())
                for _ in range(ne):
                    t = fh.readline().split()
                    etype = int(t[1])
                    ntags = int(t[2])
                    if ntags < 1:
                        continue
                    phys = int(t[3])
                    nodes = t[3 + ntags:]
                    if etype == GMSH_TRI:
                        tri_tags.append(phys)
                        tri_nodes.append((int(nodes[0]), int(nodes[1]),
                                          int(nodes[2])))
                    elif etype == GMSH_TET:
                        tet_tags[phys] += 1
                end = fh.readline().strip()
                if end != "$EndElements":
                    raise ValueError(f"expected $EndElements, got {end!r}")
            line = fh.readline()

    coords = np.asarray(coords_list, dtype=np.float64)
    tri_tags = np.asarray(tri_tags, dtype=np.int64)
    if len(tri_nodes):
        tri_idx = np.array([(id_to_idx[a], id_to_idx[b], id_to_idx[c])
                            for (a, b, c) in tri_nodes], dtype=np.int64)
    else:
        tri_idx = np.empty((0, 3), dtype=np.int64)
    return coords, tri_tags, tri_idx, tet_tags


def point_triangle_min_dist(p, A, B, C):
    """Vectorized closest-point-on-triangle distance (Ericson, RTCD).

    p : (3,) query point.  A,B,C : (N,3) triangle vertices.
    Returns (dist (N,), closest (N,3)).
    """
    N = A.shape[0]
    ab = B - A
    ac = C - A
    ap = p - A
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - B
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - C
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4

    closest = np.empty((N, 3), dtype=np.float64)
    assigned = np.zeros(N, dtype=bool)

    def safe_div(num, den):
        return num / np.where(den == 0.0, 1.0, den)

    def assign(mask, pts):
        m = mask & ~assigned
        closest[m] = pts[m]
        assigned[m] = True

    # vertex region A (1,0,0)
    assign((d1 <= 0) & (d2 <= 0), np.broadcast_to(A, (N, 3)))
    # vertex region B (0,1,0)
    assign((d3 >= 0) & (d4 <= d3), np.broadcast_to(B, (N, 3)))
    # edge region AB
    v = safe_div(d1, d1 - d3)
    assign((vc <= 0) & (d1 >= 0) & (d3 <= 0), A + v[:, None] * ab)
    # vertex region C (0,0,1)
    assign((d6 >= 0) & (d5 <= d6), np.broadcast_to(C, (N, 3)))
    # edge region AC
    w = safe_div(d2, d2 - d6)
    assign((vb <= 0) & (d2 >= 0) & (d6 <= 0), A + w[:, None] * ac)
    # edge region BC
    w_bc = safe_div(d4 - d3, (d4 - d3) + (d5 - d6))
    assign((va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0),
           B + w_bc[:, None] * (C - B))
    # interior face region (everything not yet assigned)
    denom = safe_div(np.ones(N), va + vb + vc)
    vv = vb * denom
    ww = vc * denom
    face_pts = A + ab * vv[:, None] + ac * ww[:, None]
    assign(np.ones(N, dtype=bool), face_pts)

    dist = np.linalg.norm(p - closest, axis=1)
    return dist, closest


def main(argv):
    mesh = argv[1] if len(argv) > 1 else DEFAULT_MESH
    mesh = os.path.abspath(mesh)
    print(f"mesh: {mesh}")
    if not os.path.isfile(mesh):
        print(f"FAIL: mesh not found: {mesh}")
        return 2

    coords, tri_tags, tri_nodes, tet_tags = parse_msh(mesh)
    n_nodes = coords.shape[0]
    n_tri = tri_tags.shape[0]
    print(f"nodes: {n_nodes}   triangles: {n_tri}   "
          f"tets: {sum(tet_tags.values())}")

    fails = []

    # --- Gate 1: tag inventory ------------------------------------------------
    tri_hist = Counter(int(t) for t in tri_tags.tolist())
    print("\n[1] tag inventory")
    print("  triangle tags:", dict(sorted(tri_hist.items())))
    print("  tet tags:     ", dict(sorted(tet_tags.items())))
    if set(tri_hist) != EXPECTED_TRI_TAGS:
        fails.append(f"triangle tags {set(tri_hist)} != {EXPECTED_TRI_TAGS}")
    if set(tet_tags) != EXPECTED_TET_TAGS:
        fails.append(f"tet tags {set(tet_tags)} != {EXPECTED_TET_TAGS}")

    # --- 2: bounding box ------------------------------------------------------
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    print("\n[2] bounding box")
    print(f"  x: {lo[0]:.1f} .. {hi[0]:.1f}")
    print(f"  y: {lo[1]:.1f} .. {hi[1]:.1f}")
    print(f"  z: {lo[2]:.1f} .. {hi[2]:.1f}")

    # --- 3: CVM coverage (informational) -------------------------------------
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    out_xy = ((x < GRID_X[0]) | (x > GRID_X[1]) |
              (y < GRID_Y[0]) | (y > GRID_Y[1]))
    above = z > GRID_Z[1]
    below = z < GRID_Z[0]
    print("\n[3] CVM coverage (clamping accepted; informational)")
    print(f"  grid x {GRID_X}  y {GRID_Y}  z {GRID_Z}")
    print(f"  nodes outside xy footprint: {out_xy.sum()} "
          f"({100.0 * out_xy.mean():.1f}%)")
    print(f"  nodes above z=0 (topography): {above.sum()} "
          f"({100.0 * above.mean():.1f}%)")
    print(f"  nodes below z={GRID_Z[0]:.0f}: {below.sum()} "
          f"({100.0 * below.mean():.1f}%)")

    # --- 4: hypocenter-on-fault GATE -----------------------------------------
    fmask = np.isin(tri_tags, list(FAULT_TAGS))
    fnodes = tri_nodes[fmask]
    ftags = tri_tags[fmask]
    print("\n[4] hypocenter-on-fault gate")
    if fnodes.shape[0] == 0:
        fails.append("no fault triangles found")
    else:
        A = coords[fnodes[:, 0]]
        B = coords[fnodes[:, 1]]
        C = coords[fnodes[:, 2]]
        # cheap proxy: nearest fault vertex (exact dist <= this)
        verts = coords[np.unique(fnodes)]
        vmin = np.linalg.norm(verts - HYPOCENTER, axis=1).min()
        # exact point-to-triangle
        dist, _ = point_triangle_min_dist(HYPOCENTER, A, B, C)
        j = int(np.argmin(dist))
        dmin = float(dist[j])
        strand = int(ftags[j])
        print(f"  hypocenter: {HYPOCENTER.tolist()}")
        print(f"  nearest fault vertex distance: {vmin:.2f} m")
        print(f"  exact min point-to-triangle distance: {dmin:.2f} m")
        print(f"  nearest facet strand: {strand} "
              f"({STRAND_NAME.get(strand, '?')})")
        if dmin >= HYPO_GATE_M:
            fails.append(f"hypocenter {dmin:.1f} m from nearest fault facet "
                         f">= gate {HYPO_GATE_M:.0f} m")

    # --- 5: fault not in xy-plane --------------------------------------------
    print("\n[5] fault orientation (must be near-vertical)")
    if fnodes.shape[0] > 0:
        nrm = np.cross(B - A, C - A)
        ln = np.linalg.norm(nrm, axis=1)
        good = ln > 0
        nz = np.abs(nrm[good, 2]) / ln[good]
        pct = np.percentile(nz, [5, 25, 50, 75, 95])
        print(f"  |n_z| percentiles (p5/25/50/75/95): "
              f"{pct[0]:.3f} {pct[1]:.3f} {pct[2]:.3f} {pct[3]:.3f} {pct[4]:.3f}")
        horiz = float((nz > 0.7).mean())
        print(f"  fraction near-horizontal (|n_z|>0.7): {100.0 * horiz:.2f}%")
        if np.median(nz) > 0.5:
            fails.append(f"median |n_z| {np.median(nz):.3f} > 0.5 "
                         "-- fault appears horizontal (SeisSol right-handed risk)")
        elif horiz > 0.05:
            print("  WARN: >5% of fault facets are near-horizontal "
                  "(check free-surface/DR interplay)")

    # --- verdict --------------------------------------------------------------
    print("\n" + "=" * 60)
    if fails:
        print("PREFLIGHT: FAIL")
        for f in fails:
            print("  - " + f)
        return 1
    print("PREFLIGHT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
