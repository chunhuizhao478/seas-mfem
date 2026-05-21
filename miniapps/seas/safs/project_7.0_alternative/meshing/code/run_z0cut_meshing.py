#!/usr/bin/env python3
"""run_z0cut_meshing.py — mesh the SAFS NW-cut fault cut EXACTLY at z = 0,
with z = 0 as the domain top boundary (Option A: trace-embed, gmsh Python API).

Background (2026-05-20 session; see
docs/EXPLORE_tmop_mesh_optimization.md and PLAN_mesh_quality.md):

  The production pipeline (run_nwcut_meshing.py + safs_fault_box_nwcut.geo)
  pulls the fault top DOWN to z = -fault_top_clamp (default -100 m) and lifts
  the box top by pad_top so the box top lands at z = 0.  That 100 m "headroom
  wedge" between the buried fault top and the free surface is the source of the
  near-surface sliver tets (bulk Joe-Liu eta ~ 0.09) and the fault-top "strip"
  triangles tied to the 2026-05-20 slip-rate blow-up.

  The user requires instead a surface-rupturing fault: cut EXACTLY at z = 0,
  nothing above z = 0, z = 0 IS the domain top boundary.  The naive
  `--fault-top-clamp 0 --pad-top 0` path through the CLI .geo PLC-errors
  ("A segment and a facet intersect ... on surface 6") because the OCC box top
  face is triangulated WITHOUT knowledge of the fault trace, so fault edge
  segments slice through top-face triangles -> invalid PLC -> empty volume.

  This script fixes that by reproducing the .geo's box + embedded-fault +
  size-field + physical-group setup in the gmsh Python API, and ADDITIONALLY
  building a discrete curve from the fault's z = 0 boundary edges (reusing the
  fault's own merged node tags) and embedding it into the box top face via
  gmsh.model.mesh.embed.  The top face then triangulates conformingly to the
  fault trace, eliminating the segment-facet intersection.

Hard invariants preserved (PLAN_mesh_quality.md):
  - z = 0 is the domain top boundary; mesh_zmax = 0 exactly; nothing above 0.
  - Fault stays embedded (every fault triangle is a face of two bulk tets).
  - Physical tags: rock=1, fault=101, top=102, bottom=103, sides=104.

Usage (defaults reproduce the lcfar3000 1000 m recipe):
    conda activate pythonenv
    python run_z0cut_meshing.py \
        --stl ../results/stl_nwcut/SAFS-...-ALT6_1000m_clean_clip_nwcut.stl \
        --out ../results/msh/safs_fault_box_nwcut_1000m_lcfar3000_z0embed.msh
    # fast feasibility pass (skip the production field):
    python run_z0cut_meshing.py ... --uniform-lc 3000
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import meshio
import numpy as np

import gmsh

EPS_Z = 1.0e-6          # tolerance [m] for "a node lies on the z = 0 plane"
GMSH_TRI = 2            # gmsh element type: 3-node triangle
GMSH_LINE = 1          # gmsh element type: 2-node line


def stl_bbox(stl_path: Path) -> np.ndarray:
    """Return [[xmin,ymin,zmin],[xmax,ymax,zmax]] for the STL (matches the
    convention used by run_nwcut_meshing.py:stl_bbox)."""
    m = meshio.read(str(stl_path))
    pts = m.points
    return np.array([pts.min(axis=0), pts.max(axis=0)])


def identify_box_faces(box_face_tags: list[int], zmin: float, zmax: float):
    """Classify the six OCC box faces into (top, bottom, [sides]) by their
    z bounding box.  Raises if the box is not a clean axis-aligned box with one
    face at z = zmax (top) and one at z = zmin (bottom)."""
    top, bottom, sides = None, None, []
    tol = 1.0e-6 * max(1.0, abs(zmax - zmin))
    for tag in box_face_tags:
        bb = gmsh.model.getBoundingBox(2, tag)  # (xmin,ymin,zmin,xmax,ymax,zmax)
        f_zmin, f_zmax = bb[2], bb[5]
        if abs(f_zmin - zmax) < tol and abs(f_zmax - zmax) < tol:
            if top is not None:
                raise RuntimeError("two faces found at z = zmax; box is not "
                                   "a clean axis-aligned box")
            top = tag
        elif abs(f_zmin - zmin) < tol and abs(f_zmax - zmin) < tol:
            if bottom is not None:
                raise RuntimeError("two faces found at z = zmin")
            bottom = tag
        else:
            sides.append(tag)
    if top is None or bottom is None or len(sides) != 4:
        raise RuntimeError(
            f"box-face classification failed: top={top}, bottom={bottom}, "
            f"n_sides={len(sides)} (expected top, bottom, 4 sides)")
    return top, bottom, sides


def extract_z0_trace(fault_surf: int) -> list[tuple[int, int]]:
    """Return the fault's z = 0 boundary edges as (nodeTagA, nodeTagB) pairs,
    reusing the node tags gmsh assigned when the STL was merged.

    A trace edge is a fault-triangle edge that (a) is used by exactly one
    triangle (a boundary edge of the fault patch) AND (b) has both endpoints on
    the z = 0 plane.  These are exactly the segments where the fault meets the
    free surface."""
    node_tags, coords, _ = gmsh.model.mesh.getNodes(2, fault_surf, includeBoundary=True)
    coords = np.asarray(coords).reshape(-1, 3)
    zmap = {int(t): float(coords[i, 2]) for i, t in enumerate(node_tags)}

    etypes, _, enodes = gmsh.model.mesh.getElements(2, fault_surf)
    tris = None
    for et, en in zip(etypes, enodes):
        if et == GMSH_TRI:
            tris = np.asarray(en, dtype=np.int64).reshape(-1, 3)
            break
    if tris is None or tris.shape[0] == 0:
        raise RuntimeError(f"fault surface {fault_surf} has no triangles")

    edge_count: dict[frozenset, int] = {}
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            k = frozenset((int(a), int(b)))
            edge_count[k] = edge_count.get(k, 0) + 1

    trace: list[tuple[int, int]] = []
    for k, c in edge_count.items():
        if c != 1:
            continue
        a, b = tuple(k)
        if abs(zmap.get(a, 1.0)) < EPS_Z and abs(zmap.get(b, 1.0)) < EPS_Z:
            trace.append((a, b))
    return trace


def build_size_field(fault_surf: int, args) -> int:
    """Create the background size field and return its tag.

    Reproduces safs_fault_box_nwcut.geo for the lcfar3000 (no z-graded, no
    sidecar) variant: Distance(fault) -> Threshold -> Max(., lc_min).
    With --uniform-lc > 0, use a constant field instead (feasibility runs)."""
    field = gmsh.model.mesh.field
    if args.uniform_lc > 0.0:
        f_const = field.add("MathEval")
        field.setString(f_const, "F", repr(float(args.uniform_lc)))
        return f_const

    f_dist = field.add("Distance")
    field.setNumbers(f_dist, "SurfacesList", [fault_surf])
    field.setNumber(f_dist, "Sampling", 100)

    f_thr = field.add("Threshold")
    field.setNumber(f_thr, "InField", f_dist)
    field.setNumber(f_thr, "SizeMin", args.lc_near)
    field.setNumber(f_thr, "SizeMax", args.lc_far)
    field.setNumber(f_thr, "DistMin", args.dist_inner)
    field.setNumber(f_thr, "DistMax", args.dist_outer)

    base = f_thr
    if args.lc_min > 0.0:
        f_floor = field.add("MathEval")
        field.setString(f_floor, "F", repr(float(args.lc_min)))
        f_max = field.add("Max")
        field.setNumbers(f_max, "FieldsList", [f_thr, f_floor])
        base = f_max

    # Near-surface coarsening buffer (Box field — fast, no transcendental).
    # Cutting EXACTLY at z = 0 leaves the shallow-dipping fault facets nearly
    # coplanar with the z = 0 top face; a FINE top-face mesh there triggers a
    # TetGen PLC degeneracy ("a segment and a facet intersect" / "a vertex lies
    # in a segment") — see docs/DEBUG_z0cut_fine_resolution.md.  Coarsening a
    # thin layer at z = 0 keeps the top face coarse (no degeneracy) while the
    # seismogenic bulk below stays at lc_near.  Off by default (size 0); needed
    # for lc_near << 1 km.  Box: VIn (coarse) for z > -depth, VOut=0 below, with
    # a linear transition of width `depth`; Max keeps the deeper bulk at base.
    if args.surface_buffer_size > 0.0:
        f_box = field.add("Box")
        field.setNumber(f_box, "VIn", float(args.surface_buffer_size))
        field.setNumber(f_box, "VOut", 0.0)
        field.setNumber(f_box, "XMin", -1.0e12)
        field.setNumber(f_box, "XMax", 1.0e12)
        field.setNumber(f_box, "YMin", -1.0e12)
        field.setNumber(f_box, "YMax", 1.0e12)
        field.setNumber(f_box, "ZMin", -float(args.surface_buffer_depth))
        field.setNumber(f_box, "ZMax", 1.0e12)
        field.setNumber(f_box, "Thickness", float(args.surface_buffer_depth))
        f_sbuf = field.add("Max")
        field.setNumbers(f_sbuf, "FieldsList", [base, f_box])
        base = f_sbuf
    return base


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stl", type=Path, required=True,
                    help="cut fault STL (must already reach z = 0)")
    ap.add_argument("--out", type=Path, required=True,
                    help="output .msh path (written as Gmsh v2.2)")
    ap.add_argument("--pad-x", type=float, default=50_000.0)
    ap.add_argument("--pad-y", type=float, default=50_000.0)
    ap.add_argument("--pad-bottom", type=float, default=25_000.0)
    ap.add_argument("--lc-near", type=float, default=1500.0)
    ap.add_argument("--lc-far", type=float, default=3000.0)
    ap.add_argument("--dist-inner", type=float, default=3000.0)
    ap.add_argument("--dist-outer", type=float, default=40_000.0)
    ap.add_argument("--lc-min", type=float, default=100.0,
                    help="hard lower bound on local size [m] (Field[Max]); "
                         "set 0 to disable")
    ap.add_argument("--uniform-lc", type=float, default=0.0,
                    help="if > 0, use a constant size field of this value "
                         "(skips the distance/threshold field); for fast "
                         "trace-embed feasibility checks")
    ap.add_argument("--surface-buffer-size", type=float, default=0.0,
                    help="if > 0, coarsen a near-surface layer to this size "
                         "[m] (Box field Max'd into the background) so the "
                         "exact-z=0 top face stays coarse and avoids the "
                         "TetGen PLC degeneracy at fine lc_near; 0 disables. "
                         "Use ~lc_near of a mesh that meshes cleanly (e.g. "
                         "1500).")
    ap.add_argument("--surface-buffer-depth", type=float, default=1000.0,
                    help="depth [m] of the --surface-buffer-size coarsening "
                         "layer at z=0 (linear transition over the same depth)")
    ap.add_argument("--threads", type=int, default=7)
    ap.add_argument("--no-optimize", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    if not args.stl.is_file():
        print(f"ERROR: STL not found: {args.stl}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # --- bounding box (fault bbox + pads), top pad fixed at 0 (z = 0 top) ---
    bbox = stl_bbox(args.stl)
    fxmin, fymin, fzmin = (float(v) for v in bbox[0])
    fxmax, fymax, fzmax = (float(v) for v in bbox[1])
    if abs(fzmax) > EPS_Z:
        print(f"ERROR: STL fault_zmax = {fzmax:.6f} m != 0; this mesher "
              f"requires a fault already clipped to z = 0 (use the "
              f"_clean_clip_nwcut STL, not a _zclamp variant).",
              file=sys.stderr)
        return 1
    xmin = fxmin - args.pad_x
    xmax = fxmax + args.pad_x
    ymin = fymin - args.pad_y
    ymax = fymax + args.pad_y
    zmin = fzmin - args.pad_bottom
    zmax = 0.0                      # PAD_TOP = 0: box top = z = 0 exactly
    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    print(f"fault bbox z=[{fzmin:.2f}, {fzmax:.2f}]; "
          f"mesh bbox x=[{xmin:.0f},{xmax:.0f}] y=[{ymin:.0f},{ymax:.0f}] "
          f"z=[{zmin:.0f},{zmax:.0f}]  ({dx/1e3:.1f} x {dy/1e3:.1f} x "
          f"{dz/1e3:.1f} km)")

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0 if args.quiet else 1)
        gmsh.option.setNumber("General.NumThreads", args.threads)
        gmsh.model.add("safs_z0cut")

        # --- OCC box first (faces 1..6, volume 1), matching the .geo order ---
        vol = gmsh.model.occ.addBox(xmin, ymin, zmin, dx, dy, dz)
        gmsh.model.occ.synchronize()
        box_faces = [t for (d, t) in gmsh.model.getEntities(2)]
        if len(box_faces) != 6:
            raise RuntimeError(f"expected 6 box faces, got {len(box_faces)}")
        top, bottom, sides = identify_box_faces(box_faces, zmin, zmax)

        # --- merge the cut fault (discrete surface = the new dim-2 entity) ---
        gmsh.merge(str(args.stl))
        all_surfs = [t for (d, t) in gmsh.model.getEntities(2)]
        new_surfs = [t for t in all_surfs if t not in box_faces]
        if len(new_surfs) != 1:
            raise RuntimeError(
                f"expected exactly 1 merged fault surface, got "
                f"{len(new_surfs)} (tags {new_surfs}); the STL may contain "
                f"multiple solids")
        fault_surf = new_surfs[0]

        # --- trace = fault's z = 0 boundary edges, reusing merged node tags --
        trace = extract_z0_trace(fault_surf)
        if not trace:
            raise RuntimeError(
                "no z = 0 trace edges found on the fault surface; the STL "
                "does not reach z = 0 as a boundary, so there is nothing to "
                "embed (cannot cut at z = 0)")
        trace_nodes = sorted({n for seg in trace for n in seg})
        print(f"fault surface tag={fault_surf}; box top={top}, bottom={bottom}, "
              f"sides={sides}")
        print(f"z=0 trace: {len(trace)} boundary segments, "
              f"{len(trace_nodes)} nodes")

        # GEOMETRIC trace embedded in the top face (NOT a discrete pre-meshed
        # curve).  A discrete trace curve does not constrain the OCC top-face
        # 2-D mesher (it is ignored: 0/N segments conform) and its node tags
        # are unstable across generate(), so the z = 0 seam stays open and the
        # 3-D mesher encloses nothing (empty volume).  A CAD (OCC) embed IS
        # honored, so we build geometric points+lines at the trace coords,
        # embed them in the top face, keep the fault discrete, and weld the
        # coincident top-face/fault z = 0 nodes after the 2-D mesh (below).
        # Full root-cause: docs/DEBUG_z0cut_option_a_failure.md.
        _ntags, _ncoords, _ = gmsh.model.mesh.getNodes(
            2, fault_surf, includeBoundary=True)
        _ncoords = _ncoords.reshape(-1, 3)
        tag2xy = {int(t): (float(_ncoords[i, 0]), float(_ncoords[i, 1]))
                  for i, t in enumerate(_ntags)}
        gp = {tg: gmsh.model.occ.addPoint(tag2xy[tg][0], tag2xy[tg][1], 0.0)
              for tg in trace_nodes}            # snap each trace point to z = 0
        trace_lines = [gmsh.model.occ.addLine(gp[a], gp[b]) for (a, b) in trace]
        gmsh.model.occ.synchronize()

        # --- conforming constraints ---
        # fault surface lives inside the bulk volume (== Surface{} In Volume{})
        gmsh.model.mesh.embed(2, [fault_surf], 3, vol)
        # top face must triangulate through the geometric trace (the new piece)
        gmsh.model.mesh.embed(1, trace_lines, 2, top)

        # --- size field ---
        bg = build_size_field(fault_surf, args)
        gmsh.model.mesh.field.setAsBackgroundMesh(bg)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.lc_near)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.lc_far)

        # --- physical groups (tags match the .geo) ---
        gmsh.model.addPhysicalGroup(3, [vol], 1, "rock")
        gmsh.model.addPhysicalGroup(2, [fault_surf], 101, "fault")
        gmsh.model.addPhysicalGroup(2, [top], 102, "top")
        gmsh.model.addPhysicalGroup(2, [bottom], 103, "bottom")
        gmsh.model.addPhysicalGroup(2, sides, 104, "sides")

        # --- algorithms / optimization (match the .geo) ---
        gmsh.option.setNumber("Mesh.Algorithm", 6)     # 2D Frontal-Delaunay
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # 3D Delaunay
        gmsh.option.setNumber("Mesh.Optimize", 0 if args.no_optimize else 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
        gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.3)
        gmsh.option.setNumber("Mesh.Smoothing", 10)
        gmsh.option.setNumber("Mesh.QualityType", 2)

        # 2-D mesh, then WELD the z = 0 seam before 3-D.  The top-face trace
        # nodes (on the geometric trace) and the fault's z = 0 nodes are
        # coincident but carry distinct tags; removeDuplicateNodes() merges
        # them so the volume boundary is watertight.  Without the weld the
        # 3-D mesher leaves the volume empty (Option-A failure mode).
        gmsh.model.mesh.generate(2)
        _n_pre = gmsh.model.mesh.getNodes()[0].size
        gmsh.model.mesh.removeDuplicateNodes()
        _n_welded = _n_pre - gmsh.model.mesh.getNodes()[0].size
        print(f"welded {_n_welded} coincident z=0 seam nodes "
              f"(expected ~{len(trace_nodes)})")
        gmsh.model.mesh.generate(3)

        # --- report volume fill (the PLC failure mode leaves volume empty) ---
        et, etags, _ = gmsh.model.mesh.getElements(3, vol)
        n_tet = sum(len(t) for t in etags) if etags else 0
        nb = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
        mesh_zmax = float(nb[:, 2].max()) if nb.size else float("nan")
        print(f"RESULT: volume tets = {n_tet}; mesh_zmax = {mesh_zmax:.6f}")
        if n_tet == 0:
            print("ERROR: volume 1 is EMPTY -- the trace-embed did not produce "
                  "a valid PLC (Option A failed; fall back to Option B).",
                  file=sys.stderr)
            return 2

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(args.out))
        print(f"wrote {args.out}")
    finally:
        gmsh.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main())
