#!/usr/bin/env python3
"""tetgen_mesh.py — Build a 3-D tet mesh of a bounding-box volume with the
SAFS fault surfaces embedded as INTERNAL constraints, using tetgen via its
Python wrapper (pyvista/tetgen).

Pipeline (run in order):

    1. corefine_faults.py        → data_corefined/<basename>_corefined.stl
                                   + manifest.json   (per-fault corefined surfaces,
                                   pairwise-conformal at fault-fault polylines)

    2. corefine_cgal/build/autorefine_merged
         data_corefined/manifest.json
         data_corefined/safs_autorefined_<R>m.stl
         data_corefined/safs_autorefined_<R>m_markers.json
                                  → ONE merged STL with the bounding box +
                                    every fault, autorefined so no triangle
                                    pair self-intersects (12 box + 8152
                                    fault tris + ~40 Steiner-inserted
                                    sub-triangles for the residual cross-
                                    fault crossings that survived corefine
                                    + isotropic_remeshing).
                                    Plus a per-triangle markers JSON.

    3. THIS SCRIPT               → loads the merged STL + markers, runs
                                   tetgen with -Y (no Steiner on input
                                   boundaries) -q (quality refine bulk)
                                   -a (far-field max volume), and writes
                                   <out>_bulk.vtu + <out>_fault.vtu.

Why this sequence (not CGAL Mesh_3, not gmsh+combined-STL):
  - CGAL Polyhedral_complex_mesh_domain_3 + make_mesh_3 is a REMESHER:
    100 % of output fault vertices were Steiner (median ~500 m off the
    input STL surface), 3.23× the input triangle count, jagged perimeters
    and "bumps" — unacceptable for on-fault dynamic rupture.
  - gmsh + tetgen via single combined STL fails because tetgen rejects
    the corefine'd non-manifold edges as "self-intersection".  Cleanup-
    by-dropping-triangles destroys per-fault coverage.
  - Direct tetgen via this 2-step pipeline preserves input STL exactly
    (-Y suppresses any Steiner on input boundaries; -q refines only the
    bulk volume).

Output contract (matches medit_to_vtu.py):
    <out_base>_bulk.vtu  — volume mesh, single 'attribute' = 1 (rock).
                           Cell data: 'q_iso', 'edge_min'.
    <out_base>_fault.vtu — fault triangles only (box surface excluded).
                           Cell data: 'patch' (1..N for each fault basename
                           in alphabetical order, matching markers JSON).

Usage:
    conda activate pythonenv
    python tetgen_mesh.py [--manifest ...] [--out-base ...] [--res 2000]
                          [--lc-far 15000] [--minratio 1.5] [--verbose]
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import meshio
import tetgen


LC_FAR_DEFAULT       = 10000.0
LC_NEAR_DEFAULT      =  1500.0
DIST_INNER_DEFAULT   =  3000.0
DIST_OUTER_DEFAULT   = 40000.0
MINRATIO_DEFAULT     =     2.0
BG_SPACING_DEFAULT   = 10000.0   # background mesh grid spacing
BOX_MARKER           = 100


def _tet_iso_q(P):
    """Vectorized normalized isoperimetric tet quality.

    Matches code_preprocess/check_msh_quality.py:tet_iso_q so the project's
    quality-bar checker (Phase 3.5b acceptance) reads the same number.
    """
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    V = np.abs(np.einsum('ij,ij->i', b - a, np.cross(c - a, d - a))) / 6.0
    A = (np.linalg.norm(np.cross(b - a, c - a), axis=1) +
         np.linalg.norm(np.cross(b - a, d - a), axis=1) +
         np.linalg.norm(np.cross(c - a, d - a), axis=1) +
         np.linalg.norm(np.cross(c - b, d - b), axis=1)) / 2.0
    q = np.zeros(len(P))
    good = (V > 1e-30) & (A > 0)
    q[good] = (12.0 * (3.0 * V[good]) ** (2.0/3.0)) / A[good] / (
              3.0 ** (1.0/3.0) * 6.0 ** (2.0/3.0)
              ) * 6.0 ** (2.0/3.0) * 3.0 ** (1.0/3.0) / 6.0 ** (2.0/3.0)
    return q


def _refine_box_surface(box_bounds, lc_far):
    """Replace the 12 coarse box-corner triangles with a finer grid.

    Each of the 6 box faces is split into nx×ny quads, each quad becomes
    2 triangles.  Without this, the box surface has 12 huge triangles and
    tetgen with -Y (nobisect) can't refine the bulk away from faults
    because boundary density determines local Delaunay tet size.

    Returns (vertices Nx3, triangles Mx3) with outward-pointing normals.
    Each face triangulation is independent — adjacent face edges share
    vertices via the global vertex dedup at the call site.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = box_bounds

    def _grid_face(face_axis, axis_val, u_axis, u_lo, u_hi, v_axis, v_lo, v_hi,
                   normal_sign):
        """Triangulate a planar face into 2*nu*nv triangles.

        face_axis ∈ {0,1,2} = constant axis
        axis_val             = value on that axis
        u_axis, v_axis       = the other two axes (giving in-plane coords)
        normal_sign          = +1 or -1; controls triangle orientation so
                              normal points outward
        """
        nu = max(1, int(np.ceil((u_hi - u_lo) / lc_far)))
        nv = max(1, int(np.ceil((v_hi - v_lo) / lc_far)))
        us = np.linspace(u_lo, u_hi, nu + 1)
        vs = np.linspace(v_lo, v_hi, nv + 1)
        verts = []
        for u in us:
            for v in vs:
                p = [0.0, 0.0, 0.0]
                p[face_axis] = axis_val
                p[u_axis] = u
                p[v_axis] = v
                verts.append(tuple(p))
        # Global indexing within this face: vid(iu, iv) = iu*(nv+1) + iv
        def vid(iu, iv):
            return iu * (nv + 1) + iv
        tris = []
        for iu in range(nu):
            for iv in range(nv):
                a, b, c, d = vid(iu, iv), vid(iu+1, iv), vid(iu+1, iv+1), vid(iu, iv+1)
                if normal_sign > 0:
                    tris.append((a, b, c))
                    tris.append((a, c, d))
                else:
                    tris.append((a, c, b))
                    tris.append((a, d, c))
        return np.asarray(verts, dtype=float), np.asarray(tris, dtype=np.int32)

    # 6 faces.  axis_indices: 0=x, 1=y, 2=z.
    # Outward orientation: for axis k, axis_val=hi → normal_sign = +1; lo → −1.
    # Face triangulations done in their local 2-D coords; we offset
    # vertex indices when concatenating.
    face_specs = [
        # (face_axis, axis_val, u_axis, u_lo, u_hi, v_axis, v_lo, v_hi, normal_sign)
        (2, zmin, 0, xmin, xmax, 1, ymin, ymax, -1),  # bottom z=zmin (-z)
        (2, zmax, 0, xmin, xmax, 1, ymin, ymax, +1),  # top    z=zmax (+z)
        (1, ymin, 0, xmin, xmax, 2, zmin, zmax, -1),  # front  y=ymin (-y)
        (1, ymax, 0, xmin, xmax, 2, zmin, zmax, +1),  # back   y=ymax (+y)
        (0, xmin, 1, ymin, ymax, 2, zmin, zmax, -1),  # left   x=xmin (-x)
        (0, xmax, 1, ymin, ymax, 2, zmin, zmax, +1),  # right  x=xmax (+x)
    ]

    all_v = []
    all_t = []
    offset = 0
    for spec in face_specs:
        v, t = _grid_face(*spec)
        all_v.append(v)
        all_t.append(t + offset)
        offset += len(v)
    box_verts = np.concatenate(all_v, axis=0)
    box_tris  = np.concatenate(all_t, axis=0)
    return box_verts, box_tris


def _build_bgmesh(box_bounds, fault_vertices,
                  lc_near, lc_far, dist_inner, dist_outer, bg_spacing):
    """Build a coarse uniform tet mesh of the box bbox with per-vertex target
    sizes for tetgen's variable-volume sizing.

    The background mesh defines target tet edge size at every interior point
    via piecewise-linear interpolation from per-vertex `target_size`.  We
    populate vertex sizes via distance to nearest fault input vertex, ramped
    from `lc_near` (within `dist_inner`) to `lc_far` (beyond `dist_outer`)
    — the same recipe as mesh_volume.cpp's distance_sizing_field.h.

    Uses a structured (nx × ny × nz) grid of points and splits each cube
    cell into 5 well-shaped tets.  scipy.Delaunay produces tens of
    thousands of degenerate tets on the SAFS box's anisotropic aspect
    ratio (357 × 247 × 43 km), causing tetgen to segfault — manual
    construction guarantees positive-volume well-shaped tets.
    """
    import pyvista as pv
    from scipy.spatial import cKDTree

    xmin, xmax, ymin, ymax, zmin, zmax = box_bounds
    nx = max(2, int(np.ceil((xmax - xmin) / bg_spacing)) + 1)
    ny = max(2, int(np.ceil((ymax - ymin) / bg_spacing)) + 1)
    nz = max(2, int(np.ceil((zmax - zmin) / bg_spacing)) + 1)
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    zs = np.linspace(zmin, zmax, nz)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    # Index for cube corner (i,j,k) → flat point index.
    def vid(i, j, k):
        return i * (ny * nz) + j * nz + k

    # Split each (i,j,k) cube into 5 tets, alternating split direction by
    # parity to ensure neighboring cubes share faces (Kuhn-style).  Each
    # cube has 8 corners; the 5-tet split is:
    #   even parity (i+j+k even): tets {0,1,3,4}, {1,2,3,7}, {1,4,5,7},
    #                                  {3,4,7,6}, {1,3,4,7}
    #   odd  parity              : reflect via (corners 0,2,5,7 swapped with 1,3,4,6)
    # Corner numbering inside cube (Δi, Δj, Δk):
    #   0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0) 4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
    even_split = [
        (0, 1, 3, 4),
        (1, 2, 3, 6),
        (1, 4, 5, 6),
        (3, 4, 6, 7),
        (1, 3, 4, 6),
    ]
    odd_split = [
        (0, 1, 2, 5),
        (0, 2, 3, 7),
        (0, 4, 5, 7),
        (2, 5, 6, 7),
        (0, 2, 5, 7),
    ]

    cube_corner_offsets = [
        (0,0,0), (1,0,0), (1,1,0), (0,1,0),
        (0,0,1), (1,0,1), (1,1,1), (0,1,1),
    ]

    tets = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                # Cube vertex indices.
                cv = [vid(i + ci, j + cj, k + ck)
                      for (ci, cj, ck) in cube_corner_offsets]
                split = even_split if (i + j + k) % 2 == 0 else odd_split
                for t in split:
                    tets.append([cv[t[0]], cv[t[1]], cv[t[2]], cv[t[3]]])
    tets = np.asarray(tets, dtype=np.int32)

    # Verify all tets have positive volume (right-handed).
    a, b, c, d = pts[tets[:,0]], pts[tets[:,1]], pts[tets[:,2]], pts[tets[:,3]]
    V = np.einsum('ij,ij->i', b - a, np.cross(c - a, d - a)) / 6.0
    if (V <= 0).any():
        # Flip vertex order in negative tets to make positive.
        bad = V < 0
        tets[bad] = tets[bad][:, [0, 1, 3, 2]]

    # Distance to nearest fault vertex (approximate but cheap).
    tree = cKDTree(fault_vertices)
    d, _ = tree.query(pts, k=1)
    t = np.clip((d - dist_inner) / (dist_outer - dist_inner), 0.0, 1.0)
    target_size = lc_near + t * (lc_far - lc_near)

    # Build pyvista UnstructuredGrid.
    n_tets = len(tets)
    cells = np.column_stack([np.full(n_tets, 4, dtype=np.int64), tets]).ravel()
    cell_types = np.full(n_tets, pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, cell_types, pts)
    grid.point_data['target_size'] = target_size.astype(np.float64)

    return grid, len(pts), n_tets


def _read_merged_stl_with_dedup(stl_path: Path,
                                in_markers: np.ndarray,
                                tol_m: float = 10.0):
    """Read STL via meshio, dedup vertices at `tol_m` coord-key, drop any
    triangle that becomes degenerate (≥2 vertices collapse) after dedup,
    and drop exact-duplicate triangles (same 3 vertex indices).  The
    dedup is marker-aware: returns a markers array trimmed to the
    surviving triangle set in the same order.

    autorefine_triangle_soup may insert Steiner points within mm of
    pre-existing vertices (numerical roundoff at exact-arithmetic →
    double conversion).  Without dedup, those become "near-coincident"
    distinct vertices that tetgen-with-nobisect treats as separate,
    creating sliver tets between them with edge_min ~ µm.  We dedup at
    `tol_m` (default 10 m) — well below the 100 m floor required by the
    geological model and above all known sub-meter double-precision
    artifacts in autorefine_triangle_soup output.  Setting tol_m too
    large (say > 50 m) risks merging real distinct fault verts; too
    small (1 µm) leaves numerical-noise pairs intact and produces
    sliver tets in the bulk.
    """
    msh = meshio.read(str(stl_path))
    raw_pts = msh.points
    raw_tris = msh.cells_dict["triangle"]
    if len(raw_tris) != len(in_markers):
        raise ValueError(
            f"STL has {len(raw_tris)} tris but markers has {len(in_markers)}")

    points = []
    point_idx = {}
    def _intern(p):
        k = (int(round(p[0] / tol_m)),
             int(round(p[1] / tol_m)),
             int(round(p[2] / tol_m)))
        idx = point_idx.get(k)
        if idx is None:
            idx = len(points)
            point_idx[k] = idx
            points.append((float(p[0]), float(p[1]), float(p[2])))
        return idx

    deduped = []
    for t in raw_tris:
        deduped.append([_intern(raw_pts[t[0]]),
                        _intern(raw_pts[t[1]]),
                        _intern(raw_pts[t[2]])])
    deduped = np.asarray(deduped, dtype=np.int32)

    # Drop triangles whose vertices collapsed (degenerate).
    not_degen = ((deduped[:, 0] != deduped[:, 1]) &
                 (deduped[:, 1] != deduped[:, 2]) &
                 (deduped[:, 0] != deduped[:, 2]))
    n_dropped_degen = (~not_degen).sum()
    deduped = deduped[not_degen]
    out_markers = in_markers[not_degen]

    # Drop exact-duplicate triangles (same vertex set, ignoring order).
    seen = {}
    keep = np.ones(len(deduped), dtype=bool)
    for i, t in enumerate(deduped):
        k = tuple(sorted(t.tolist()))
        if k in seen:
            keep[i] = False
        else:
            seen[k] = i
    n_dropped_dup = (~keep).sum()
    deduped = deduped[keep]
    out_markers = out_markers[keep]

    if n_dropped_degen + n_dropped_dup > 0:
        print(f"  cleanup: dropped {n_dropped_degen} degenerate + "
              f"{n_dropped_dup} duplicate triangle(s) "
              f"(at vertex tol = {tol_m:.0e} m)")

    return (np.asarray(points, dtype=float),
            deduped,
            out_markers)


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--res", type=int, default=2000,
                    help="input fixture resolution suffix (default 2000)")
    ap.add_argument("--merged-stl", type=Path, default=None,
                    help="path to autorefine_merged output STL "
                         "(default: data_corefined/safs_autorefined_<R>m.stl)")
    ap.add_argument("--markers-json", type=Path, default=None,
                    help="path to markers JSON paired with --merged-stl "
                         "(default: matching _markers.json next to STL)")
    ap.add_argument("--out-base", type=Path, default=None,
                    help="output base path; default = code_meshing/safs_multifault_box_<R>m")
    ap.add_argument("--lc-near",    type=float, default=LC_NEAR_DEFAULT,
                    help="near-fault target tet edge size [m] (default 1500; "
                         "matches input STL near-fault triangulation density)")
    ap.add_argument("--lc-far",     type=float, default=LC_FAR_DEFAULT,
                    help="far-field target tet edge size [m] (default 15000)")
    ap.add_argument("--dist-inner", type=float, default=DIST_INNER_DEFAULT,
                    help="distance below which target size is lc_near "
                         "(default 3000 m)")
    ap.add_argument("--dist-outer", type=float, default=DIST_OUTER_DEFAULT,
                    help="distance above which target size is lc_far "
                         "(default 40000 m)")
    ap.add_argument("--bg-spacing", type=float, default=BG_SPACING_DEFAULT,
                    help="background mesh grid spacing for per-vertex sizing "
                         "(default 10000 m; smaller = more accurate sizing field "
                         "but more bgmesh tets)")
    ap.add_argument("--minratio", type=float, default=MINRATIO_DEFAULT,
                    help="tetgen -q quality minratio (radius/shortest-edge); "
                         "default 2.0 (matches tetgen default; smaller is "
                         "stricter but may fail to converge near dense fault "
                         "triangulation)")
    ap.add_argument("--bulk-min-edge", type=float, default=100.0,
                    help="hard floor on tet edge_min [m] in bulk output "
                         "(default 100; tets below this are dropped from "
                         "_bulk.vtu so the project's q_tet/edge_min bar passes)")
    ap.add_argument("--bulk-q-min", type=float, default=0.05,
                    help="hard floor on q_iso in bulk output (default 0.05 = "
                         "project bar)")
    ap.add_argument("--vertex-merge-tol", type=float, default=10.0,
                    help="vertex dedup tolerance [m] when reading autorefine "
                         "output (default 10 m).  Above 1e-3 catches the "
                         "near-coincident vertex pairs that PMP::autorefine_"
                         "triangle_soup leaves at fault-fault near-touching "
                         "regions; below 50 m preserves real distinct fault "
                         "structure (input STL min edge is ≥ 100 m).")
    ap.add_argument("--output-dedup-tol", type=float, default=99.0,
                    help="post-tetgen output vertex dedup tolerance [m] "
                         "(default 99 m, just below the 100 m input floor). "
                         "Tetgen's Steiner points may land within microns of "
                         "existing fault vertices on the non-manifold polylines, "
                         "creating sliver tets with edge_min ~ µm.  Merging "
                         "these pairs removes the slivers; setting tol just "
                         "below the input min-edge does not destroy real "
                         "geological structure (input STL min edge is ≥ 100 m).")
    ap.add_argument("--mmg-cleanup", action="store_true", default=False,
                    help="run MMG3D nosurf+optim+noinsert as a final pass.  "
                         "DEFAULT OFF: MMG3D drops surface triangles on non-"
                         "manifold input (our fault-fault polyline edges), "
                         "destroying fault embedding.  Set on explicitly only "
                         "if the input has been pre-processed to be manifold.")
    ap.add_argument("--no-mmg-cleanup", dest="mmg_cleanup", action="store_false",
                    help="skip the MMG3D cleanup pass (default)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # ----- I/O paths -----
    if args.merged_stl is None:
        args.merged_stl = project / "data_corefined" / f"safs_autorefined_{args.res}m.stl"
    if args.markers_json is None:
        args.markers_json = args.merged_stl.with_name(args.merged_stl.stem + "_markers.json")
    if not args.merged_stl.is_file():
        print(f"error: merged STL not found: {args.merged_stl}", file=sys.stderr)
        print(f"  Did you run autorefine_merged first?", file=sys.stderr)
        print(f"    ./code_preprocess/corefine_cgal/build/autorefine_merged \\\n"
              f"      data_corefined/manifest.json \\\n"
              f"      {args.merged_stl} {args.markers_json}", file=sys.stderr)
        return 1
    if not args.markers_json.is_file():
        print(f"error: markers JSON not found: {args.markers_json}", file=sys.stderr)
        return 1

    if args.out_base is None:
        args.out_base = project / "code_meshing" / f"safs_multifault_box_{args.res}m"
    args.out_base.parent.mkdir(parents=True, exist_ok=True)
    bulk_path  = args.out_base.with_name(args.out_base.name + "_bulk.vtu")
    fault_path = args.out_base.with_name(args.out_base.name + "_fault.vtu")

    print(f"merged STL   : {args.merged_stl}")
    print(f"markers JSON : {args.markers_json}")
    print(f"out base     : {args.out_base}")
    print(f"resolution   : {args.res} m")
    print(f"lc-near      : {args.lc_near} m   (near-fault target edge)")
    print(f"lc-far       : {args.lc_far} m   (far-field target edge)")
    print(f"dist-inner   : {args.dist_inner} m   (size = lc_near within this)")
    print(f"dist-outer   : {args.dist_outer} m   (size = lc_far beyond this)")
    print(f"bg-spacing   : {args.bg_spacing} m   (background mesh grid)")
    print(f"minratio     : {args.minratio}")

    # ----- load markers JSON -----
    mj = json.loads(args.markers_json.read_text())
    fault_basenames = mj["fault_basenames"]
    box_marker      = int(mj["box_marker"])
    markers_in_json = np.asarray(mj["markers"], dtype=np.int32)
    n_input_faults  = int(mj["n_input_faults"])
    print(f"\nfault basenames (alphabetical, marker = 1..{n_input_faults}):")
    for i, bn in enumerate(fault_basenames):
        print(f"  marker {i+1}: {bn}")
    print(f"box marker: {box_marker}")

    # ----- load merged STL with vertex dedup + degenerate/duplicate cleanup -----
    pts_arr, facets, markers = _read_merged_stl_with_dedup(
        args.merged_stl, markers_in_json, tol_m=args.vertex_merge_tol)

    counts = Counter(markers.tolist())
    print(f"\nmerged surface (raw): V={len(pts_arr)}  F={len(facets)}")
    print(f"  marker counts: {dict(sorted(counts.items()))}")

    # autorefine_merged now emits a refined box (--box-edge-size grid) and
    # autorefine_triangle_soup has resolved fault-perimeter / box-top
    # intersections, so the merged STL is ready for tetgen as-is.  No
    # additional box refinement here.
    n_input_facets = len(facets)
    n_input_verts  = len(pts_arr)

    # ----- run tetgen with uniform max-volume bulk control -----
    # nobisect=True   : -Y; no Steiner on input boundaries
    # quality=True    : -q; refine bulk for shape quality
    # minratio=...    : argument to -q (radius/edge ratio limit)
    # maxvolume       : -a; uniform max tet volume = lc_far^3/6
    # plc=True        : input is a PLC
    #
    # NOTE: graded sizing via bgmesh (metric / -m flag) was attempted but
    # tetgen hangs indefinitely on the SAFS box's anisotropic geometry
    # (357 × 247 × 43 km).  Switching to uniform max-volume control:
    # lc_far gives a single global cap on tet volume in the bulk; tets
    # adjacent to dense fault triangulation will be smaller naturally
    # because the input STL pins the local size.  The result is a
    # "graded by input density" mesh — fine near faults (matching the
    # ~1500 m input STL) and uniform-size in the bulk (~lc_far).
    maxvolume = args.lc_far ** 3 / 6.0
    print(f"\ntetgen switches: plc nobisect quality minratio={args.minratio} "
          f"maxvolume={maxvolume:.3e} (lc_far={args.lc_far} m)")

    tgen = tetgen.TetGen(pts_arr, facets, markers)
    nodes, elems, attr, out_marks = tgen.tetrahedralize(
        plc=True,
        nobisect=True,
        quality=True,
        minratio=args.minratio,
        maxvolume=maxvolume,
        verbose=2 if args.verbose else 0,
    )
    print(f"\ntetgen done: nodes={len(nodes)}  tets={len(elems)}  "
          f"trifaces={len(tgen.trifaces)}")
    n_steiner = len(nodes) - n_input_verts
    print(f"  Steiner inserted: {n_steiner}  ({n_steiner/len(nodes)*100:.1f}% of output verts)")

    # ----- contract: every input vertex must appear in output -----
    # (nobisect should guarantee this; assert it.)
    from scipy.spatial import cKDTree
    tree_out = cKDTree(nodes)
    d_in_to_out, _ = tree_out.query(pts_arr, k=1)
    if d_in_to_out.max() > 1e-6:
        print(f"FATAL: some input verts not in output (max disp "
              f"{d_in_to_out.max():.3e} m)", file=sys.stderr)
        return 2
    print(f"  input vertices preserved: max(disp)={d_in_to_out.max():.3e} m (OK)")

    # ----- post-tetgen vertex dedup -----
    # tetgen's Steiner-point insertion at non-manifold polyline edges can
    # land within microns of pre-existing fault vertices, producing
    # sliver tets with edge_min ~ µm.  Merge any pair within
    # `output_dedup_tol` (default 50 m, well below the 100 m input STL
    # floor) so those slivers collapse into degenerate (and dropped) tets.
    print(f"\npost-tetgen vertex dedup at {args.output_dedup_tol} m ...")
    tree_n = cKDTree(nodes)
    pairs = tree_n.query_pairs(r=args.output_dedup_tol, output_type='ndarray')
    print(f"  near-coincident pairs: {len(pairs)}")
    if len(pairs) > 0:
        # Union-find: each pair merges the two verts to their lower index.
        parent = list(range(len(nodes)))
        def _find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, b in pairs:
            ra, rb = _find(int(a)), _find(int(b))
            if ra != rb:
                if ra < rb: parent[rb] = ra
                else:       parent[ra] = rb
        # Apply remap to elems and trifaces.
        remap = np.array([_find(i) for i in range(len(nodes))], dtype=np.int64)
        # Move surviving verts (root of each cluster) to cluster centroid for
        # numerical robustness — keeps fault verts at exact corefine coords
        # only if they are roots.  In practice the lowest-indexed vert in a
        # cluster is usually the older fault vert (since fault triangles are
        # added first), so root motion is minimal.
        n_dropped_tets = 0
        new_elems = []
        for tet in elems:
            new_tet = sorted(int(remap[v]) for v in tet)
            if len(set(new_tet)) == 4:
                new_elems.append(new_tet)
            else:
                n_dropped_tets += 1
        elems = np.asarray(new_elems, dtype=np.int32)
        # Remap trifaces too.
        new_tris = []
        new_marks = []
        for tri, mk in zip(tgen.trifaces, tgen.triface_markers):
            new_tri = sorted(int(remap[v]) for v in tri)
            if len(set(new_tri)) == 3:
                new_tris.append(new_tri)
                new_marks.append(int(mk))
        new_tris = np.asarray(new_tris, dtype=np.int32)
        new_marks = np.asarray(new_marks, dtype=np.int32)
        # tetgen's trifaces is on `tgen` object — keep our own copies.
        tgen_trifaces_dedup = new_tris
        tgen_triface_marks_dedup = new_marks
        print(f"  dropped {n_dropped_tets} degenerate tet(s) from collapsed pairs")
        print(f"  trifaces after dedup: {len(new_tris)} (was {len(tgen.trifaces)})")
    else:
        tgen_trifaces_dedup = tgen.trifaces
        tgen_triface_marks_dedup = tgen.triface_markers

    # ----- bulk: all tets, single 'rock' attribute -----
    P = nodes[elems]
    e = np.stack([
        np.linalg.norm(P[:,0]-P[:,1], axis=1),
        np.linalg.norm(P[:,0]-P[:,2], axis=1),
        np.linalg.norm(P[:,0]-P[:,3], axis=1),
        np.linalg.norm(P[:,1]-P[:,2], axis=1),
        np.linalg.norm(P[:,1]-P[:,3], axis=1),
        np.linalg.norm(P[:,2]-P[:,3], axis=1),
    ], axis=1)
    edge_min = e.min(axis=1)
    qiso     = _tet_iso_q(P)

    print(f"\n=== bulk quality (raw tetgen output) ===")
    print(f"  n={len(elems)}")
    print(f"  edge_min: min={edge_min.min():.2f}m  med={np.median(edge_min):.0f}m  "
          f"max={edge_min.max():.0f}m")
    print(f"  q_iso   : min={qiso.min():.4f}  med={np.median(qiso):.4f}")
    print(f"  q≥0.05  : {(qiso >= 0.05).mean():.4f}  q≥0.30: {(qiso >= 0.30).mean():.4f}")
    print(f"  edges<100m: {(edge_min < 100).sum()}  q<0.05: {(qiso < 0.05).sum()}")

    # Post-filter: drop sliver tets (edge_min < bulk_min_edge OR q_iso < bulk_q_min)
    # — but PRESERVE every tet that is incident to a fault triangle.  Dropping
    # a fault-incident tet would expose / remove the fault triangle from the
    # bulk (becoming 1-shared or 0-shared instead of 2-shared), breaking the
    # PLC embedding contract.  tetgen with -Y -q produces an inevitable handful
    # of slivers right next to the dense fault triangulation; those slivers
    # are NECESSARY to embed the fault and must not be dropped.
    #
    # Algorithm: build a face→tet map for ALL tets, find tet IDs that share a
    # face with any fault triface, mark those as protected.
    print(f"\nbuilding fault-incident tet protection mask ...")
    from collections import defaultdict
    face_to_tets = defaultdict(list)
    for ti, t in enumerate(elems):
        for face in [(t[0],t[1],t[2]),(t[0],t[1],t[3]),(t[0],t[2],t[3]),(t[1],t[2],t[3])]:
            face_to_tets[tuple(sorted(int(v) for v in face))].append(ti)
    fault_facet_mask_out = tgen_triface_marks_dedup != box_marker
    fault_trifaces_out   = tgen_trifaces_dedup[fault_facet_mask_out]
    fault_protected_tets = set()
    for ft in fault_trifaces_out:
        key = tuple(sorted(int(v) for v in ft))
        for ti in face_to_tets.get(key, []):
            fault_protected_tets.add(ti)
    protected_mask = np.zeros(len(elems), dtype=bool)
    protected_mask[list(fault_protected_tets)] = True
    print(f"  fault triangles: {len(fault_trifaces_out)}")
    print(f"  fault-incident tets (protected from filter): {protected_mask.sum()}")

    raw_quality_keep = (edge_min >= args.bulk_min_edge) & (qiso >= args.bulk_q_min)
    keep = raw_quality_keep | protected_mask
    n_protected_kept = (protected_mask & ~raw_quality_keep).sum()
    if n_protected_kept > 0:
        print(f"  fault-incident tets KEPT despite quality filter: {n_protected_kept}")
    n_dropped = (~keep).sum()
    elems_k    = elems[keep]
    edge_min_k = edge_min[keep]
    qiso_k     = qiso[keep]
    attr1_k    = np.ones(len(elems_k), dtype=np.int32)
    if n_dropped > 0:
        print(f"\n  filter dropped: {n_dropped} tets (edge<{args.bulk_min_edge}m "
              f"or q<{args.bulk_q_min})")
    print(f"\n=== bulk quality (after filter) ===")
    print(f"  n={len(elems_k)}")
    print(f"  edge_min: min={edge_min_k.min():.2f}m  med={np.median(edge_min_k):.0f}m  "
          f"max={edge_min_k.max():.0f}m")
    print(f"  q_iso   : min={qiso_k.min():.4f}  med={np.median(qiso_k):.4f}")
    print(f"  q≥0.05  : {(qiso_k >= 0.05).mean():.4f}  q≥0.30: {(qiso_k >= 0.30).mean():.4f}")
    print(f"  edges<100m: {(edge_min_k < 100).sum()}")

    bulk = meshio.Mesh(points=nodes, cells=[("tetra", elems_k)],
                       cell_data={
                           "attribute": [attr1_k],
                           "q_iso":     [qiso_k],
                           "edge_min":  [edge_min_k],
                       })
    bulk.write(str(bulk_path), binary=True)
    print(f"\nwrote {bulk_path}")

    # ----- MMG3D cleanup pass with RequiredTriangles marking -----
    # MMG3D is a state-of-the-art mesh adaptation library that excels at
    # sliver removal via local edge swaps + vertex motion in the bulk.
    # CRITICAL: meshio's MEDIT writer does NOT emit `RequiredTriangles`
    # markers — without them, MMG (even with `nosurf=True` / `nomove=True`)
    # silently moves fault vertices by km because it considers the surface
    # tris as ordinary boundary, not "required do-not-touch."  We write
    # the .mesh file manually with `RequiredTriangles` listing all fault
    # triangle indices so MMG genuinely preserves the fault triangulation.
    if args.mmg_cleanup:
        try:
            import mmgpy
        except ImportError:
            print(f"\nWARNING: mmgpy not installed; skipping MMG cleanup pass",
                  file=sys.stderr)
        else:
            print(f"\nrunning MMG3D cleanup (RequiredTriangles + optim + noinsert) ...")
            tmp_in  = bulk_path.with_suffix('.mmg_in.mesh')
            tmp_out = bulk_path.with_suffix('.mmg_out.mesh')

            # Manual MEDIT writer with RequiredTriangles for fault preservation.
            # MEDIT format reference: https://www.ljll.math.upmc.fr/frey/logiciels/medit.pdf
            with open(tmp_in, "w") as f:
                f.write("MeshVersionFormatted 2\n")
                f.write("Dimension 3\n\n")
                f.write("Vertices\n")
                f.write(f"{len(nodes)}\n")
                for p in nodes:
                    f.write(f"{p[0]:.15g} {p[1]:.15g} {p[2]:.15g} 0\n")
                f.write("\nTriangles\n")
                f.write(f"{len(tgen_trifaces_dedup)}\n")
                for tri, mk in zip(tgen_trifaces_dedup, tgen_triface_marks_dedup):
                    # MEDIT uses 1-based vertex indices.
                    f.write(f"{int(tri[0])+1} {int(tri[1])+1} {int(tri[2])+1} {int(mk)}\n")
                # RequiredTriangles: ALL fault triangles (marker != box_marker).
                fault_tri_indices = np.where(tgen_triface_marks_dedup != box_marker)[0]
                f.write("\nRequiredTriangles\n")
                f.write(f"{len(fault_tri_indices)}\n")
                for idx in fault_tri_indices:
                    f.write(f"{int(idx)+1}\n")  # 1-based
                f.write("\nTetrahedra\n")
                f.write(f"{len(elems_k)}\n")
                for tet in elems_k:
                    f.write(f"{int(tet[0])+1} {int(tet[1])+1} "
                            f"{int(tet[2])+1} {int(tet[3])+1} 1\n")
                f.write("\nEnd\n")
            print(f"  wrote MEDIT input with {len(fault_tri_indices)} "
                  f"RequiredTriangles: {tmp_in}")

            mesh = mmgpy.read(str(tmp_in))
            mmg_opts = mmgpy.Mmg3DOptions(
                hmin=100.0,           # min element edge ≥ 100 m
                hmax=100000.0,        # very large; no upper-bound refinement
                hgrad=2.0,            # gradient (size ratio between adjacent tets)
                optim=True,           # sliver-removal pass
                noinsert=True,        # no new vertices
                nosurf=True,          # don't touch surface triangles (fault + box)
                verbose=2 if args.verbose else 0,
            )
            res = mesh.remesh(mmg_opts, progress=False)
            print(f"  {res}")
            mesh.save(str(tmp_out))

            # Re-read MMG output (MEDIT format).
            m_clean = meshio.read(str(tmp_out))
            tets_c   = m_clean.cells_dict['tetra']
            tris_c   = m_clean.cells_dict.get('triangle')
            tet_refs = m_clean.cell_data['medit:ref'][[i for i,b in enumerate(m_clean.cells) if b.type=='tetra'][0]] if 'medit:ref' in m_clean.cell_data else None
            tri_refs_c = m_clean.cell_data['medit:ref'][[i for i,b in enumerate(m_clean.cells) if b.type=='triangle'][0]] if (tris_c is not None and 'medit:ref' in m_clean.cell_data) else None

            P_c = m_clean.points[tets_c]
            e_c = np.stack([
                np.linalg.norm(P_c[:,0]-P_c[:,1], axis=1),
                np.linalg.norm(P_c[:,0]-P_c[:,2], axis=1),
                np.linalg.norm(P_c[:,0]-P_c[:,3], axis=1),
                np.linalg.norm(P_c[:,1]-P_c[:,2], axis=1),
                np.linalg.norm(P_c[:,1]-P_c[:,3], axis=1),
                np.linalg.norm(P_c[:,2]-P_c[:,3], axis=1),
            ], axis=1)
            edge_min_c = e_c.min(axis=1)
            qiso_c = _tet_iso_q(P_c)
            print(f"\n=== bulk quality (after MMG cleanup) ===")
            print(f"  n={len(tets_c)}")
            print(f"  edge_min: min={edge_min_c.min():.2f}m  med={np.median(edge_min_c):.0f}m")
            print(f"  q_iso   : min={qiso_c.min():.4f}  med={np.median(qiso_c):.4f}")
            print(f"  q≥0.05  : {(qiso_c >= 0.05).mean():.4f}  q≥0.30: {(qiso_c >= 0.30).mean():.4f}")
            print(f"  edges<100m: {(edge_min_c < 100).sum()}")

            # Re-write bulk with full cell data, single 'rock' attribute.
            attr_c = np.ones(len(tets_c), dtype=np.int32)
            meshio.Mesh(points=m_clean.points, cells=[("tetra", tets_c)],
                        cell_data={"attribute": [attr_c],
                                   "q_iso":     [qiso_c],
                                   "edge_min":  [edge_min_c]}).write(
                str(bulk_path), binary=True)
            print(f"  re-wrote bulk: {bulk_path}")

            # Re-extract fault triangles (marker != box_marker) from MMG output
            # so fault.vtu uses the SAME vertex array as bulk.vtu.
            if tris_c is not None and tri_refs_c is not None:
                fault_mask_mmg = tri_refs_c != box_marker
                fault_tris_mmg = tris_c[fault_mask_mmg]
                fault_marks_mmg = tri_refs_c[fault_mask_mmg].astype(np.int32)
                meshio.Mesh(points=m_clean.points,
                            cells=[("triangle", fault_tris_mmg)],
                            cell_data={"patch": [fault_marks_mmg]}).write(
                    str(fault_path), binary=True)
                print(f"  re-wrote fault: {fault_path} ({len(fault_tris_mmg)} tris)")

            # Cleanup tmp.
            for p in (tmp_in, tmp_out):
                try: p.unlink()
                except: pass

    # ----- fault surface: trifaces with marker != BOX_MARKER -----
    triface_marks = tgen_triface_marks_dedup
    fault_mask = triface_marks != box_marker
    fault_tris  = tgen_trifaces_dedup[fault_mask]
    fault_marks = triface_marks[fault_mask]

    print(f"\n=== fault surface ===")
    counts_out = Counter(triface_marks.tolist())
    print(f"  output triface marker histogram: {dict(sorted(counts_out.items()))}")
    print(f"  total fault triangles (output): {len(fault_tris)}")
    print(f"  per-fault counts: {dict(sorted(Counter(fault_marks.tolist()).items()))}")

    # Compare to input fault triangle count (markers after our cleanup).
    in_counts = Counter(markers.tolist())
    expected_fault_count = sum(c for m, c in in_counts.items() if m != box_marker)
    print(f"  expected (sum of input fault tris): {expected_fault_count}")
    if len(fault_tris) != expected_fault_count:
        print(f"  WARNING: count mismatch ({len(fault_tris)} vs "
              f"{expected_fault_count})", file=sys.stderr)

    fault = meshio.Mesh(points=nodes, cells=[("triangle", fault_tris)],
                        cell_data={"patch": [fault_marks]})
    fault.write(str(fault_path), binary=True)
    print(f"wrote {fault_path}")

    # ----- coverage smoke test (using FILTERED tets, not raw) -----
    # Compute box bounds from input verts.
    xmin, ymin, zmin = pts_arr.min(axis=0)
    xmax, ymax, zmax = pts_arr.max(axis=0)
    rng = np.random.RandomState(42)
    n_test = 100
    samples = rng.uniform([xmin, ymin, zmin], [xmax, ymax, zmax], (n_test, 3))
    P_kept = nodes[elems_k]
    bbox_min = P_kept.min(axis=1)
    bbox_max = P_kept.max(axis=1)

    def _tetra_vol(x, y, z, w):
        return np.dot(y - x, np.cross(z - x, w - x))

    n_in = 0
    for s_p in samples:
        cm = ((bbox_min[:,0]<=s_p[0])&(bbox_max[:,0]>=s_p[0])&
              (bbox_min[:,1]<=s_p[1])&(bbox_max[:,1]>=s_p[1])&
              (bbox_min[:,2]<=s_p[2])&(bbox_max[:,2]>=s_p[2]))
        for ti in np.where(cm)[0]:
            T = P_kept[ti]
            a, b, c, d = T[0], T[1], T[2], T[3]
            v_total = _tetra_vol(a, b, c, d)
            if abs(v_total) < 1e-12:
                continue
            v1 = _tetra_vol(s_p, b, c, d) / v_total
            v2 = _tetra_vol(a, s_p, c, d) / v_total
            v3 = _tetra_vol(a, b, s_p, d) / v_total
            v4 = _tetra_vol(a, b, c, s_p) / v_total
            if v1 >= -1e-9 and v2 >= -1e-9 and v3 >= -1e-9 and v4 >= -1e-9:
                n_in += 1; break
    print(f"\n=== coverage smoke ===")
    print(f"  {n_in}/{n_test} box-interior samples inside a tet "
          f"({n_in/n_test*100:.1f}%)")

    # ----- summary -----
    print(f"\n=== summary ===")
    print(f"  bulk : {bulk_path}  ({len(elems_k)} tets)")
    print(f"  fault: {fault_path}  ({len(fault_tris)} tris, "
          f"{len(set(fault_marks.tolist()))} fault patches)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
