"""Constrained re-triangulation: insert cross-fault polylines as exact
mesh edges on both faults.

Implements PLAN_multifault_intersections.md Phase 2.

Pipeline per pair (A, B):
    polylines = chain_segments(tri_tri_intersect_3d(A, B))
    for each polyline P:
        for each fault F in (A, B):
            split parents pierced by P (Step 1)
            propagate edge-pierces to in-fault neighbours (Step 2)
    write per-fault conformal STL + triangle_to_fault.json + report

Key invariants (verified by the validation gates in this file):
  - manifold per fault: every edge shared by ≤ 2 triangles, no
    T-junctions (no vertex sits on another edge's interior).
  - polyline edge-coincidence: every consecutive polyline vertex pair
    appears as an edge in BOTH faults' post-split meshes.
  - interior-crossing-only gate: post-split, no interior crossings
    remain (shared edges along the polyline are correct, not failures).
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# Phase 4 of PLAN_cgal_corefine.md — the production conformalizer is
# now `tools/corefine_faults` (CGAL).  This Python implementation is
# preserved as a regression fixture (per P-013, the 90-day retention
# clock) but new pipelines should not call it.  See
# `mesh/tests/README.md` for the legacy-test policy.
warnings.warn(
    "conformalize_faults.py is the legacy Python conformalizer; "
    "the production path is tools/corefine_faults (CGAL).  "
    "See PLAN_cgal_corefine.md.",
    DeprecationWarning,
    stacklevel=2,
)

try:
    import triangle as triangle_lib                  # noqa: F401
except ImportError as exc:                           # pragma: no cover
    raise ImportError(
        "conformalize_faults requires Shewchuk's `triangle` library. "
        "Install with `pip install triangle` in the pythonenv conda env."
    ) from exc

# Phase 1 module — shared dataclasses, geometric primitives, STL reader.
import fault_intersect as fi

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-parent local frame and 2-D projection
# ---------------------------------------------------------------------------
@dataclass
class _LocalFrame:
    origin: np.ndarray   # (3,) — chosen as parent v0
    u: np.ndarray        # (3,) unit vector along edge (v0, v1)
    v: np.ndarray        # (3,) unit vector in plane, orthogonal to u


def _make_local_frame(parent_verts_3d: np.ndarray) -> _LocalFrame:
    """Build an orthonormal 2-D basis on the triangle's plane.

    Origin = parent v0; u-axis along edge (v0 → v1); v-axis in-plane and
    perpendicular to u.  Coordinates in this frame are O(edge_length),
    keeping float64 round-off below the SAFS coordinate scale (~1e5 m)
    times machine eps (~1e-16) ≈ 1e-11 m, well below any tolerance.
    """
    v0 = parent_verts_3d[0]
    e1 = parent_verts_3d[1] - v0
    e2 = parent_verts_3d[2] - v0
    e1_norm = float(np.linalg.norm(e1))
    if e1_norm == 0.0:
        raise ValueError("parent triangle has zero-length edge v0-v1")
    u = e1 / e1_norm
    n = np.cross(e1, e2)
    n_norm = float(np.linalg.norm(n))
    if n_norm == 0.0:
        raise ValueError("parent triangle is degenerate (zero area)")
    n = n / n_norm
    v = np.cross(n, u)                                # already unit
    return _LocalFrame(origin=v0, u=u, v=v)


def _to_2d(p_3d: np.ndarray, frame: _LocalFrame) -> tuple[float, float]:
    rel = p_3d - frame.origin
    return float(rel @ frame.u), float(rel @ frame.v)


def _to_3d(x: float, y: float, frame: _LocalFrame) -> np.ndarray:
    return frame.origin + x * frame.u + y * frame.v


# ---------------------------------------------------------------------------
# Pierce classification: on-edge / interior
# ---------------------------------------------------------------------------
@dataclass
class _PierceLocation:
    """Where a pierce sits relative to its parent triangle."""
    edge_idx: int | None     # 0/1/2 for edge (v0-v1)/(v1-v2)/(v2-v0); None if interior or vertex
    edge_param: float         # parameter along the edge (0.0 at start vertex, 1.0 at end), in [0, 1]
    vertex_idx: int | None    # 0/1/2 if the pierce coincides with a parent vertex; else None


def _classify_pierce_2d(
    p2d: tuple[float, float],
    parent_2d: list[tuple[float, float]],
    edge_tol: float,
    vertex_tol: float,
) -> _PierceLocation:
    """Categorise a 2-D pierce point relative to a 2-D parent triangle.

    Returns:
        - vertex_idx ∈ {0, 1, 2} if the pierce coincides (within
          vertex_tol) with a parent vertex.  edge_idx is set to the
          incident edge with edge_param ∈ {0.0, 1.0}.
        - edge_idx ∈ {0, 1, 2} with edge_param ∈ (0, 1) if the pierce
          lies on the interior of that edge.
        - edge_idx = None if the pierce is strictly interior to the
          parent triangle.
    """
    px, py = p2d
    # Vertex coincidence first — vertex_tol can be small.
    for i, (vx, vy) in enumerate(parent_2d):
        if (px - vx) ** 2 + (py - vy) ** 2 < vertex_tol ** 2:
            # Pick the lexicographically first edge incident to this vertex.
            edge_idx = i                          # edge starting at vertex i
            return _PierceLocation(
                edge_idx=edge_idx, edge_param=0.0, vertex_idx=i,
            )
    # Edge proximity.
    edges = ((0, 1), (1, 2), (2, 0))
    best_edge = None
    best_dist_sq = math.inf
    best_t = 0.0
    for ei, (i, j) in enumerate(edges):
        ax, ay = parent_2d[i]
        bx, by = parent_2d[j]
        ex, ey = bx - ax, by - ay
        edge_len_sq = ex * ex + ey * ey
        if edge_len_sq == 0.0:
            continue
        t = ((px - ax) * ex + (py - ay) * ey) / edge_len_sq
        # We want the perpendicular distance even if t is outside [0, 1]
        # — out-of-range t means the pierce projection falls past an
        # endpoint, in which case the perpendicular distance to the
        # infinite line is what we use to classify.  Vertex-coincidence
        # would have been caught above.
        if t < 0.0:
            cx, cy = ax, ay
        elif t > 1.0:
            cx, cy = bx, by
        else:
            cx, cy = ax + t * ex, ay + t * ey
        d_sq = (px - cx) ** 2 + (py - cy) ** 2
        if d_sq < best_dist_sq:
            best_dist_sq = d_sq
            best_edge = ei
            best_t = t
    if best_edge is not None and best_dist_sq < edge_tol ** 2 \
            and 0.0 < best_t < 1.0:
        return _PierceLocation(
            edge_idx=best_edge, edge_param=best_t, vertex_idx=None,
        )
    return _PierceLocation(edge_idx=None, edge_param=0.0, vertex_idx=None)


# ---------------------------------------------------------------------------
# Step (1) — split one parent triangle
# ---------------------------------------------------------------------------
def _split_one_parent(
    parent_verts_3d: np.ndarray,        # (3, 3)
    pierce_pts_3d: list[np.ndarray],    # K points; ALL must already be unique global pierces
    polyline_edges: list[tuple[int, int]],  # indices into pierce_pts_3d (may repeat)
    *,
    edge_tol_frac: float = 5e-3,        # fraction of parent's mean edge → 2-D edge tolerance (0.5%)
    vertex_tol_frac: float = 1e-3,      # fraction of parent's mean edge → vertex coincidence
    cdt_min_angle_deg: float | None = None,    # knob B: enforce min angle by adding Steiner points
    target_edge_length_m: float | None = None, # knob B/C: enforce max area
    edge_tol: float | None = None,      # legacy absolute override (m); takes precedence if not None
    vertex_tol: float | None = None,    # legacy absolute override (m); takes precedence if not None
) -> tuple[
    list[tuple[int, int, int]],   # children as (la, lb, lc) LOCAL CDT input indices
    list[int],                    # pierce_to_local: for each input pierce_pts_3d[i],
                                  #     its local CDT index (0/1/2 if vertex-coincident,
                                  #     3..3+K_unique-1 otherwise)
    list[np.ndarray],             # NEW interior Steiner points in 3-D added by the CDT
]:
    """Re-triangulate a parent triangle to include the pierce points and
    polyline edges as constraints.

    Parameters
    ----------
    parent_verts_3d : (3, 3)
        Parent triangle vertices in 3-D, in their original orientation.
    pierce_pts_3d : list of length K
        Polyline vertex coordinates that touch this parent (entry, exit,
        and any interior kink point).  Must already be deduplicated at
        the caller's snap_m precision (the caller passes globally-unique
        nodes).
    polyline_edges : list of (a, b)
        Each entry is a pair of indices into ``pierce_pts_3d`` that must
        appear as a mesh edge after re-triangulation.

    Returns
    -------
    list of (A, B, C) — the 3-D children, each an ordered triple of
    triangle vertices in 3-D.

    Algorithm (PLAN Phase 2 step 1):
        a. project parent + pierces to the parent's local 2-D frame
        b. classify each pierce: on which edge (or interior?)
        c. build a closed boundary in CCW order, inserting on-edge
           pierces in their natural order along each edge
        d. add polyline-internal segments as PSLG constraints
        e. call ``triangle.triangulate(..., 'p')`` (no 'q' — that flag
           inserts Steiner points on constraints and re-breaks
           conformality)
        f. assert input vertices and segments are preserved 1:1
        g. unproject children back to 3-D
    """
    if len(parent_verts_3d) != 3:
        raise ValueError("parent_verts_3d must have shape (3, 3)")

    parent_3d_arr = np.asarray(parent_verts_3d, dtype=np.float64)
    frame = _make_local_frame(parent_3d_arr)
    parent_2d = [_to_2d(parent_verts_3d[i], frame) for i in range(3)]

    # Compute scaled tolerances based on the parent's mean edge length.
    # Real CFM data has parents ~1500 m on a side; synthetic test
    # fixtures have parents ~1 m.  A fixed absolute tolerance can't
    # serve both; the fraction scales naturally.  Cross-parent pierces
    # (from pierces_per_edge of a non-coplanar neighbour) project up
    # to ~dihedral_angle × edge_length off the on-edge line in 2-D —
    # for CFM dihedrals (~10°) and 1500 m edges that's ~18 cm.  A
    # 0.5%-of-edge default (= 7.5 m for 1500 m edges, 5 mm for 1 m
    # edges) absorbs this without false positives.
    pmean = fi._mean_edge(parent_3d_arr)
    eff_edge_tol = (edge_tol if edge_tol is not None
                    else edge_tol_frac * pmean)
    eff_vertex_tol = (vertex_tol if vertex_tol is not None
                      else vertex_tol_frac * pmean)

    pierces_2d = [_to_2d(np.asarray(p, dtype=np.float64), frame)
                  for p in pierce_pts_3d]

    # Classify each pierce; for on-edge pierces, SNAP to the edge so
    # the boundary builder sees them exactly on the edge (avoids the
    # CDT receiving a near-boundary interior vertex with a polyline-
    # edge constraint that crosses the parent boundary).
    locations: list[_PierceLocation] = []
    for k, p2d in enumerate(pierces_2d):
        loc = _classify_pierce_2d(
            p2d, parent_2d,
            edge_tol=eff_edge_tol, vertex_tol=eff_vertex_tol,
        )
        locations.append(loc)
        if loc.vertex_idx is None and loc.edge_idx is not None:
            # Snap pierce 2-D position onto the edge at edge_param.
            edges_def_2d = ((0, 1), (1, 2), (2, 0))
            i_start, i_end = edges_def_2d[loc.edge_idx]
            ax, ay = parent_2d[i_start]
            bx, by = parent_2d[i_end]
            t = loc.edge_param
            pierces_2d[k] = (ax + t * (bx - ax), ay + t * (by - ay))

    # Build vertex array (parent first, then pierces) and a node->2D map.
    # Pierces that classify as a parent vertex are remapped to that vertex
    # (i.e., the pierce's "node id" in the CDT input is the parent's
    # vertex index, not a new index).  For interior or on-edge pierces
    # we add new vertices.
    pts_2d: list[tuple[float, float]] = list(parent_2d)
    pierce_to_local: list[int] = []
    for loc, p2d in zip(locations, pierces_2d):
        if loc.vertex_idx is not None:
            pierce_to_local.append(loc.vertex_idx)
        else:
            pierce_to_local.append(len(pts_2d))
            pts_2d.append(p2d)

    # Build boundary: walk parent edges, inserting on-edge pierces in
    # parameter-order along each edge.
    edges_def = ((0, 1), (1, 2), (2, 0))
    boundary_segs: list[tuple[int, int]] = []
    for ei, (i_start, i_end) in enumerate(edges_def):
        # Collect (param, local_idx) for every pierce that lies on this edge.
        on_edge: list[tuple[float, int]] = []
        for k, loc in enumerate(locations):
            if loc.vertex_idx is not None:
                continue                          # vertex-coincident; not a boundary insertion
            if loc.edge_idx == ei:
                on_edge.append((loc.edge_param, pierce_to_local[k]))
        on_edge.sort(key=lambda t: t[0])
        prev = i_start
        for _, idx in on_edge:
            if idx != prev:
                boundary_segs.append((prev, idx))
                prev = idx
        if i_end != prev:
            boundary_segs.append((prev, i_end))

    # Build polyline-internal segments (those that connect two distinct
    # local indices and are not already on the boundary).
    boundary_set = {frozenset(s) for s in boundary_segs}
    interior_segs: list[tuple[int, int]] = []
    for a, b in polyline_edges:
        la = pierce_to_local[a]
        lb = pierce_to_local[b]
        if la == lb:
            continue                              # zero-length after vertex collapse
        if frozenset((la, lb)) in boundary_set:
            continue                              # already on the boundary, no need to repeat
        interior_segs.append((la, lb))

    all_segs = boundary_segs + interior_segs
    seg_arr = np.asarray(all_segs, dtype=np.int32)
    pts_arr = np.asarray(pts_2d, dtype=np.float64)
    n_input_verts = pts_arr.shape[0]

    # Knob B: ``Y`` forbids Steiner points on segments (boundary +
    # polyline-internal), preserving all input edges 1:1.  ``q{angle}``
    # enforces a minimum angle by inserting Steiner points in the
    # interior; ``a{area}`` adds an upper bound on triangle area to
    # control mesh density.  When neither is requested, fall back to
    # bare ``'p'`` (no quality control, no Steiner).
    flags = "p"
    if cdt_min_angle_deg is not None or target_edge_length_m is not None:
        # Single 'Y' protects only boundary segments; 'YY' protects ALL
        # segments (boundary + internal/polyline).  We need YY because
        # the polyline edges are internal segments that the q/a flags
        # would otherwise subdivide, breaking cross-fault polyline
        # edge-coincidence.
        flags += "YY"
        if cdt_min_angle_deg is not None:
            if not (0.0 < cdt_min_angle_deg < 34.0):
                # Triangle's q-flag is unstable above 34° (Shewchuk's
                # documented theoretical limit ~33.8°); reject loudly.
                raise ValueError(
                    f"cdt_min_angle_deg must be in (0, 34); "
                    f"got {cdt_min_angle_deg}"
                )
            flags += f"q{cdt_min_angle_deg:.4f}"
        if target_edge_length_m is not None:
            if target_edge_length_m <= 0.0:
                raise ValueError(
                    f"target_edge_length_m must be > 0; "
                    f"got {target_edge_length_m}"
                )
            # Equilateral triangle with edge L has area sqrt(3)/4 * L^2.
            max_area = (math.sqrt(3.0) / 4.0) * target_edge_length_m ** 2
            flags += f"a{max_area:.6f}"

    out = triangle_lib.triangulate(
        {"vertices": pts_arr, "segments": seg_arr},
        flags,
    )

    # Step 1d — assert input vertices and segments preserved 1:1.
    out_verts = out["vertices"]
    if out_verts.shape[0] < n_input_verts:
        raise RuntimeError(
            f"triangle.triangulate dropped vertices "
            f"({n_input_verts} input → {out_verts.shape[0]} output); "
            f"input vertices must be preserved."
        )
    # Input vertices must appear unchanged in the first n_input_verts rows.
    if not np.allclose(out_verts[:n_input_verts], pts_arr, atol=1e-12):
        raise RuntimeError(
            "triangle.triangulate reordered or perturbed input vertices; "
            "polyline edge-coincidence will be broken."
        )
    # All input segments must be present in the output (Steiner-on-segments
    # is forbidden by 'Y' when quality control is on).
    out_segs = {frozenset((int(a), int(b))) for a, b in out["segments"]}
    in_segs = {frozenset((int(a), int(b))) for a, b in all_segs}
    missing = in_segs - out_segs
    if missing:
        raise RuntimeError(
            f"triangle.triangulate dropped or subdivided constraint "
            f"segments: {sorted(missing)} not in output.  Polyline "
            f"edge-coincidence will be broken."
        )

    # Extract any new Steiner vertices added by the quality CDT.
    n_steiner = int(out_verts.shape[0]) - n_input_verts
    new_steiner_3d: list[np.ndarray] = []
    for k in range(n_steiner):
        sx, sy = out_verts[n_input_verts + k]
        new_steiner_3d.append(_to_3d(float(sx), float(sy), frame))

    # Step 1g — return children as LOCAL CDT INPUT INDICES.  The caller
    # maps these to global node IDs via:
    #   - local idx 0/1/2 → parent's vertices (a_orig, b_orig, c_orig)
    #   - local idx in [3, n_input_verts) → pierce, via pierce_to_local
    #   - local idx in [n_input_verts, n_input_verts + n_steiner) →
    #     fresh global ID (caller appends Steiner to V_global)
    children_indices: list[tuple[int, int, int]] = [
        (int(t[0]), int(t[1]), int(t[2])) for t in out["triangles"]
    ]
    return children_indices, pierce_to_local, new_steiner_3d


# ---------------------------------------------------------------------------
# Per-fault: split all pierced parents and propagate to neighbours.
# ---------------------------------------------------------------------------
def _conformalize_one_fault(
    V_orig: np.ndarray, T_orig: np.ndarray,
    polylines: list[fi.Polyline],
    *,
    use_pierces_a: bool,
    snap_m: float = 1e-2,
    drop_short_segment_frac: float = 0.0,
    edge_tol_local: float = 1e-3,
    vertex_tol_local: float = 1e-6,
    cdt_min_angle_deg: float | None = None,
    target_edge_length_m: float | None = None,
) -> tuple[np.ndarray, np.ndarray, list[list[int]], int]:
    """Apply Phase 2 Step (1) split + Step (2) propagate to one fault.

    Returns (V_new, T_new, polyline_to_node, n_dropped_short_segments).

    polyline_to_node[i] is a list of V_new indices, one per
    polylines[i].points entry, identifying the post-conformal mesh
    vertex that each polyline vertex landed on.
    """
    V_orig = np.asarray(V_orig, dtype=np.float64)
    T_orig = np.asarray(T_orig, dtype=np.int64)
    n_v_orig = V_orig.shape[0]

    # ---- A: catalog polyline vertices as global nodes. ----
    V_list: list[np.ndarray] = [V_orig[i].copy() for i in range(n_v_orig)]
    snap_table: dict[tuple[int, int, int], int] = {}

    def snap_key(p):
        return (int(round(p[0] / snap_m)),
                int(round(p[1] / snap_m)),
                int(round(p[2] / snap_m)))

    polyline_to_node: list[list[int]] = []
    for pl in polylines:
        nodes_for_pl: list[int] = []
        for p in pl.points:
            arr = np.asarray(p, dtype=np.float64)
            k = snap_key(arr)
            if k not in snap_table:
                snap_table[k] = len(V_list)
                V_list.append(arr.copy())
            nodes_for_pl.append(snap_table[k])
        polyline_to_node.append(nodes_for_pl)

    # ---- B: group polyline segments by parent triangle. ----
    parent_segments: dict[int, list[tuple[int, int]]] = {}
    n_dropped = 0
    for pl, nodes in zip(polylines, polyline_to_node):
        pierces = pl.pierces_a if use_pierces_a else pl.pierces_b
        if len(pierces) != len(nodes) - 1:
            _log.warning(
                "polyline (%s × %s) has %d nodes but %d pierces; "
                "expected %d.  Skipping its segment attribution.",
                pl.fault_a, pl.fault_b, len(nodes), len(pierces),
                len(nodes) - 1,
            )
            continue
        for i, tri_id in enumerate(pierces):
            seg_a = nodes[i]
            seg_b = nodes[i + 1]
            if seg_a == seg_b:
                continue
            parent_mean = fi._mean_edge(V_orig[T_orig[tri_id]])
            seg_len = float(np.linalg.norm(V_list[seg_b] - V_list[seg_a]))
            if seg_len < drop_short_segment_frac * parent_mean:
                n_dropped += 1
                continue
            parent_segments.setdefault(int(tri_id), []).append(
                (seg_a, seg_b))

    # ---- C: process each pierced parent.  For each pierced parent,
    # gather both its OWN polyline pierces AND any cross-parent edge
    # pierces (pre-collected below in a tiny first pass).  Without
    # the cross-parent gathering, two pierced parents sharing an edge
    # can produce inconsistent boundary subdivisions → T-junction.
    pierced_parent_set = set(parent_segments.keys())
    pierced_children: dict[int, list[tuple[int, int, int]]] = {}
    pierces_per_edge: dict[tuple[int, int], list[tuple[float, int]]] = {}

    # First, classify every pierce of every pierced parent against
    # that parent's edges to populate pierces_per_edge for cross-
    # parent gathering.
    for tri_id in pierced_parent_set:
        a_orig, b_orig, c_orig = (int(x) for x in T_orig[tri_id])
        parent_verts_3d = V_orig[T_orig[tri_id]]
        pierce_node_ids = sorted({n for s in parent_segments[tri_id] for n in s})
        frame = _make_local_frame(parent_verts_3d)
        parent_2d = [_to_2d(parent_verts_3d[i], frame) for i in range(3)]
        for nid in pierce_node_ids:
            p2d = _to_2d(np.asarray(V_list[nid]), frame)
            loc = _classify_pierce_2d(
                p2d, parent_2d, edge_tol=edge_tol_local,
                vertex_tol=vertex_tol_local,
            )
            if loc.vertex_idx is not None or loc.edge_idx is None:
                continue
            ei = loc.edge_idx
            uv = ((a_orig, b_orig), (b_orig, c_orig), (c_orig, a_orig))[ei]
            u, v = uv
            key = (min(u, v), max(u, v))
            t_canon = loc.edge_param if u <= v else 1.0 - loc.edge_param
            pierces_per_edge.setdefault(key, []).append((t_canon, nid))

    # Dedup pierces_per_edge.
    for key, lst in pierces_per_edge.items():
        seen: dict[int, float] = {}
        for t, nid in lst:
            if nid not in seen:
                seen[nid] = t
        pierces_per_edge[key] = sorted(
            ((t, nid) for nid, t in seen.items()), key=lambda x: x[0]
        )

    # Now run the CDT split for each pierced parent.
    for tri_id in pierced_parent_set:
        a_orig, b_orig, c_orig = (int(x) for x in T_orig[tri_id])
        parent_verts_3d = V_orig[T_orig[tri_id]]
        seg_indices = parent_segments[tri_id]
        own_pierce_node_ids: set[int] = {n for s in seg_indices for n in s}

        # Gather cross-parent edge pierces on this parent's boundary.
        edges_def = ((a_orig, b_orig), (b_orig, c_orig), (c_orig, a_orig))
        for u, v in edges_def:
            key = (min(u, v), max(u, v))
            for _, nid in pierces_per_edge.get(key, []):
                own_pierce_node_ids.add(nid)

        pierce_node_ids = sorted(own_pierce_node_ids)
        pierce_pts_3d = [V_list[nid] for nid in pierce_node_ids]
        global_to_local = {nid: li for li, nid in enumerate(pierce_node_ids)}
        local_polyline_edges = [
            (global_to_local[a], global_to_local[b]) for a, b in seg_indices
        ]

        # _split_one_parent returns children as LOCAL CDT input indices
        # plus a pierce_to_local mapping; we map directly to global IDs
        # without 3-D coordinate matching.  This avoids the post-snap
        # mismatch where pierce_pts_3d entries were the pre-snap 3-D
        # coords but child vertices unproject from snapped 2-D.
        children_local, pierce_to_local, new_steiner_3d = _split_one_parent(
            parent_verts_3d, pierce_pts_3d, local_polyline_edges,
            cdt_min_angle_deg=cdt_min_angle_deg,
            target_edge_length_m=target_edge_length_m,
        )

        # Build local-index → global-ID table.
        # Local 0/1/2 → parent's verts.
        # Local 3..3+K_unique-1 → pierces, but the mapping is via
        #   pierce_to_local: pierce_node_ids[i] is at local pierce_to_local[i].
        #   (Multiple input pierces may share one local index if they
        #   were vertex-coincident — that's fine, all such pierces have
        #   the same global node id by definition of the parent vertex.)
        # n_input_verts = 3 + (# unique non-vertex pierces).
        n_input_verts = 3 + len({li for li in pierce_to_local if li >= 3})
        local_to_global: list[int] = [a_orig, b_orig, c_orig]
        # Fill in pierce slots: walk pierce_to_local, assign global IDs
        # by local index.  Multiple input pierces sharing a local idx
        # >= 3 must agree on the global id (they geometrically coincide).
        local_to_global_map: dict[int, int] = {0: a_orig, 1: b_orig, 2: c_orig}
        for input_pierce_idx, li in enumerate(pierce_to_local):
            if li in (0, 1, 2):
                continue                          # already set to parent vertex
            global_id = pierce_node_ids[input_pierce_idx]
            if li in local_to_global_map:
                if local_to_global_map[li] != global_id:
                    raise RuntimeError(
                        f"local idx {li} maps to two different global "
                        f"ids ({local_to_global_map[li]} and {global_id})"
                    )
            else:
                local_to_global_map[li] = global_id
        # Append in order (li 0,1,2,...,n_input_verts-1).
        local_to_global = [local_to_global_map[li] for li in range(n_input_verts)]
        # Append Steiner with fresh global IDs.
        for s_3d in new_steiner_3d:
            new_global_id = len(V_list)
            V_list.append(s_3d.copy())
            local_to_global.append(new_global_id)

        # Map child local indices to global IDs.
        children_indexed: list[tuple[int, int, int]] = []
        for la, lb, lc in children_local:
            try:
                ga = local_to_global[la]
                gb = local_to_global[lb]
                gc = local_to_global[lc]
            except IndexError:
                raise RuntimeError(
                    f"child has out-of-range local index for parent "
                    f"{tri_id}: ({la}, {lb}, {lc}) vs |local_to_global|"
                    f"={len(local_to_global)}"
                )
            children_indexed.append((ga, gb, gc))
        pierced_children[tri_id] = children_indexed

    # ---- D: walk T_orig and emit final triangle list. ----
    T_new_list: list[tuple[int, int, int]] = []
    for tri_id in range(T_orig.shape[0]):
        a_orig, b_orig, c_orig = (int(x) for x in T_orig[tri_id])
        if tri_id in pierced_parent_set:
            T_new_list.extend(pierced_children[tri_id])
            continue
        # Unpierced parent: check whether any of its 3 edges has pierces
        # from a NEIGHBOUR's split.  If so, propagate.
        edges_def = ((a_orig, b_orig), (b_orig, c_orig), (c_orig, a_orig))
        edge_pierces: list[list[int]] = []          # per edge, sorted node ids along (u → v)
        any_pierces = False
        for u, v in edges_def:
            key = (min(u, v), max(u, v))
            lst = pierces_per_edge.get(key, [])
            if not lst:
                edge_pierces.append([])
                continue
            any_pierces = True
            # Order pierces along u → v.
            if u <= v:
                ordered = [nid for (_, nid) in lst]
            else:
                ordered = [nid for (_, nid) in reversed(lst)]
            edge_pierces.append(ordered)
        quality_on = (cdt_min_angle_deg is not None
                      or target_edge_length_m is not None)
        if not any_pierces and not quality_on:
            # Fast path: unpierced parent, no quality refinement → keep
            # the original triangle.
            T_new_list.append((a_orig, b_orig, c_orig))
            continue
        # Otherwise pass through the CDT helper so the parent gets
        # boundary pierces inserted (if any) AND quality-refined (if
        # the knobs are on).  When quality is on, the helper may
        # append Steiner points to V_list (mutates the caller's list
        # in place) and returns children with already-mapped global
        # indices.
        children, _new_steiner = _fan_split_with_edge_pierces(
            (a_orig, b_orig, c_orig), edges_def, edge_pierces, V_list,
            cdt_min_angle_deg=cdt_min_angle_deg,
            target_edge_length_m=target_edge_length_m,
        )
        T_new_list.extend(children)

    V_new = np.asarray(V_list, dtype=np.float64)
    T_new = np.asarray(T_new_list, dtype=np.int64)
    return V_new, T_new, polyline_to_node, n_dropped


def _fan_split_with_edge_pierces(
    parent_abc: tuple[int, int, int],
    edges_def: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    edge_pierces: list[list[int]],
    V_global: list[np.ndarray],
    *,
    cdt_min_angle_deg: float | None = None,
    target_edge_length_m: float | None = None,
) -> tuple[
    list[tuple[int, int, int]],
    list[np.ndarray],
]:
    """Step (2) split: insert boundary pierces into an unpierced parent.

    Edge cases handled:
      - 0 pierces total: returns ``[parent_abc]`` unchanged.
      - All pierces on one edge: fan from the opposite vertex (the
        plan's BAD/GOOD picture).
      - Pierces on multiple edges: runs a 2-D constrained Delaunay
        triangulation of the parent's perimeter (closed boundary,
        no interior segments, no internal vertices) so the resulting
        triangulation is correct regardless of where the pierces are.

    Returns a list of (i, j, k) global-index triangles.

    Why a CDT-based fallback rather than always-fan: a naive fan from
    perimeter[0] produces a degenerate (collinear) child whenever
    perimeter[0] sits on a pierced edge — the bug the previous
    implementation had.  Picking the apex as "the vertex opposite the
    edge with the most pierces" works for the single-pierced-edge case
    but is fragile when multiple edges are pierced; instead we let the
    CDT choose the diagonals.
    """
    a, b, c = parent_abc
    n_total = sum(len(p) for p in edge_pierces)
    quality_on = (cdt_min_angle_deg is not None
                  or target_edge_length_m is not None)
    if n_total == 0 and not quality_on:
        return [parent_abc], []

    # Without quality flags AND with a single pierced edge, fan from
    # the opposite vertex (cheapest path; no Steiner points).
    pierced_edge_indices = [i for i, p in enumerate(edge_pierces) if p]
    if not quality_on and len(pierced_edge_indices) == 1:
        ei = pierced_edge_indices[0]
        u, v = edges_def[ei]
        apex = (c, a, b)[ei]
        pierces = edge_pierces[ei]
        children: list[tuple[int, int, int]] = []
        prev = u
        for p in pierces:
            children.append((prev, p, apex))
            prev = p
        children.append((prev, v, apex))
        return children, []

    # General / quality-controlled case: build the parent's perimeter
    # (on-edge pierces inserted), project to 2-D, run constrained
    # Delaunay (with optional quality flags), and emit children.
    perimeter: list[int] = []
    for ei, (u, _v) in enumerate(edges_def):
        perimeter.append(u)
        perimeter.extend(edge_pierces[ei])
    parent_3d = np.asarray(
        [V_global[a], V_global[b], V_global[c]], dtype=np.float64
    )
    frame = _make_local_frame(parent_3d)
    pts_2d = [_to_2d(np.asarray(V_global[gid], dtype=np.float64), frame)
              for gid in perimeter]
    n = len(perimeter)
    boundary_segs = [(i, (i + 1) % n) for i in range(n)]

    flags = "p"
    if quality_on:
        flags += "Y"
        if cdt_min_angle_deg is not None:
            flags += f"q{cdt_min_angle_deg:.4f}"
        if target_edge_length_m is not None:
            max_area = (math.sqrt(3.0) / 4.0) * target_edge_length_m ** 2
            flags += f"a{max_area:.6f}"

    out = triangle_lib.triangulate(
        {"vertices": np.asarray(pts_2d, dtype=np.float64),
         "segments": np.asarray(boundary_segs, dtype=np.int32)},
        flags,
    )
    out_verts = out["vertices"]
    if out_verts.shape[0] < n:
        raise RuntimeError(
            f"propagation CDT dropped vertices for parent {parent_abc}: "
            f"{n} input → {out_verts.shape[0]} output"
        )
    if not np.allclose(out_verts[:n], np.asarray(pts_2d), atol=1e-12):
        raise RuntimeError(
            f"propagation CDT reordered input vertices for parent "
            f"{parent_abc}"
        )
    n_steiner = int(out_verts.shape[0]) - n

    # Allocate global IDs for any new Steiner vertices (mutating V_global
    # — V_global is the caller's V_list, kept as a list so we can append).
    new_steiner_3d: list[np.ndarray] = []
    steiner_global_ids: list[int] = []
    for k in range(n_steiner):
        sx, sy = out_verts[n + k]
        s3 = _to_3d(float(sx), float(sy), frame)
        new_steiner_3d.append(s3)
        steiner_global_ids.append(len(V_global))
        V_global.append(s3)

    local_to_global = list(perimeter) + steiner_global_ids
    children = [
        (local_to_global[int(t[0])],
         local_to_global[int(t[1])],
         local_to_global[int(t[2])])
        for t in out["triangles"]
    ]
    return children, new_steiner_3d


def _smooth_fault_surface(
    V: np.ndarray, T: np.ndarray,
    locked_vertex_ids: set[int],
    *,
    n_iters: int = 5,
    relaxation: float = 0.5,
) -> np.ndarray:
    """Uniform-Laplacian smoothing of a triangulated fault surface.

    Moves each non-locked vertex toward the centroid of its 1-ring
    neighbours, ``n_iters`` times.  Vertices in ``locked_vertex_ids``
    stay put — these must include polyline vertices (so polyline
    edge-coincidence is preserved) and fault-boundary vertices (so
    the fault's outline is preserved).

    Smoothing happens in 3-D; for shallow surfaces (fault triangles
    nearly planar locally) this stays close to the original surface
    geometry.  Vertex moves are bounded per iteration by the
    ``relaxation`` parameter (0.5 = move halfway to centroid each
    step) so the surface doesn't drift far from the original CFM
    geometry — typical max move on SAFS data is ~10 m at step 1,
    decaying to <1 m by step 5.

    The function does NOT change T (no edge flips, no remeshing),
    so the polyline-edge-coincidence and manifold properties are
    preserved by construction.
    """
    if n_iters <= 0:
        return V.copy()
    V_out = np.asarray(V, dtype=np.float64).copy()
    T = np.asarray(T, dtype=np.int64)

    # Build vertex adjacency (1-ring neighbours).
    adj: list[set[int]] = [set() for _ in range(V_out.shape[0])]
    for tri in T:
        a, b, c = (int(x) for x in tri)
        adj[a].update((b, c))
        adj[b].update((a, c))
        adj[c].update((a, b))

    locked = locked_vertex_ids                    # alias

    for _ in range(n_iters):
        V_new = V_out.copy()
        for vi in range(V_out.shape[0]):
            if vi in locked:
                continue
            nbrs = adj[vi]
            if not nbrs:
                continue
            centroid = np.zeros(3)
            for ni in nbrs:
                centroid += V_out[ni]
            centroid /= len(nbrs)
            V_new[vi] = (1.0 - relaxation) * V_out[vi] + relaxation * centroid
        V_out = V_new
    return V_out


def _constrained_delaunay_flips(
    V: np.ndarray, T: np.ndarray,
    locked_edges: set[tuple[int, int]],
    *,
    max_passes: int = 8,
) -> np.ndarray:
    """Constrained-Delaunay edge-flip optimisation.

    For each interior edge (shared by exactly 2 triangles within the
    fault) that is NOT in ``locked_edges`` (polyline edges + fault-
    boundary edges), test the Lawson flip criterion:

        flip improves quality iff
            min(angle in current 4-tri-fan)
              <  min(angle in flipped 4-tri-fan)

    If so, swap the diagonal:

           u                u
           /\              /|\\
          /  \\           / | \\
         /    \\        /  |  \\
        a━━━━━b   →   a   |   b
         \\    /        \\  |  /
          \\  /          \\ | /
           \\/            \\|/
           v                v

       T1 = (a, b, u), T2 = (b, a, v)
                ↓
       T1' = (a, v, u), T2' = (b, u, v)

    Iterate until a full pass produces zero flips OR ``max_passes``
    is reached.

    Constraints preserved:
      - Polyline edges (in locked_edges): never flipped → polyline
        edge-coincidence with the OTHER fault stays intact.
      - Fault-boundary edges (1-incidence): never flipped (have only
        one adjacent triangle).
      - Manifold property: flipping a 2-incidence edge yields a
        2-incidence edge — manifold preserved.

    Returns the new ``T`` array (V is unchanged).
    """
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64).copy()
    if T.shape[0] == 0:
        return T

    def tri_min_angle_3d(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
        """Smallest interior angle of triangle (p0, p1, p2), in radians.
        Returns 0 for degenerate triangles."""
        e01 = p1 - p0; e12 = p2 - p1; e20 = p0 - p2
        L01 = float(np.linalg.norm(e01))
        L12 = float(np.linalg.norm(e12))
        L20 = float(np.linalg.norm(e20))
        if min(L01, L12, L20) <= 0.0:
            return 0.0
        # Use law of cosines for stability.
        a01 = math.acos(max(-1.0, min(1.0,
            ( L20 ** 2 + L01 ** 2 - L12 ** 2) / (2.0 * L20 * L01))))
        a12 = math.acos(max(-1.0, min(1.0,
            ( L01 ** 2 + L12 ** 2 - L20 ** 2) / (2.0 * L01 * L12))))
        a20 = math.pi - a01 - a12
        return min(a01, a12, a20)

    n_total_flips = 0
    for _pass in range(max_passes):
        # Build edge → list of (tri_idx, opposite_vertex) up-to-date.
        edge_to_tris: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for ti in range(T.shape[0]):
            a, b, c = (int(T[ti, 0]), int(T[ti, 1]), int(T[ti, 2]))
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                key = (min(u, v), max(u, v))
                edge_to_tris.setdefault(key, []).append((ti, w))

        flips_this_pass = 0
        flipped_tris: set[int] = set()
        for key, ts in edge_to_tris.items():
            if key in locked_edges:
                continue
            if len(ts) != 2:
                continue                          # boundary or non-manifold
            (ti1, w1), (ti2, w2) = ts
            if ti1 in flipped_tris or ti2 in flipped_tris:
                continue                          # already touched this pass
            u, v = key
            # Original triangles share edge (u, v).
            old_min = min(
                tri_min_angle_3d(V[u], V[v], V[w1]),
                tri_min_angle_3d(V[u], V[v], V[w2]),
            )
            # Flipped triangles share new diagonal (w1, w2).
            new_min = min(
                tri_min_angle_3d(V[u], V[w1], V[w2]),
                tri_min_angle_3d(V[v], V[w1], V[w2]),
            )
            if new_min <= old_min + 1e-12:
                continue                          # no improvement
            # Reject flip if it would create a degenerate triangle.
            if new_min < 1e-9:
                continue
            # Apply flip: rewrite both triangles.
            T[ti1] = (u, w1, w2)
            T[ti2] = (v, w2, w1)
            flipped_tris.add(ti1)
            flipped_tris.add(ti2)
            flips_this_pass += 1
        n_total_flips += flips_this_pass
        if flips_this_pass == 0:
            break
    return T


def _identify_locked_edges(
    T: np.ndarray,
    polyline_node_ids_per_pl: list[list[int]],
) -> set[tuple[int, int]]:
    """Collect edges that the flip pass must NOT touch:
      (1) every consecutive polyline-vertex pair (the polyline edges
          inserted into both faults — flipping these would break
          cross-fault polyline-edge coincidence),
      (2) every fault-boundary edge (incidence == 1; flipping means
          shrinking or distorting the fault outline).
    """
    locked: set[tuple[int, int]] = set()
    for nodes in polyline_node_ids_per_pl:
        for i in range(len(nodes) - 1):
            u, v = nodes[i], nodes[i + 1]
            if u == v:
                continue
            locked.add((min(u, v), max(u, v)))
    edge_count: dict[tuple[int, int], int] = {}
    for tri in T:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            edge_count[key] = edge_count.get(key, 0) + 1
    for k, n in edge_count.items():
        if n == 1:
            locked.add(k)
    return locked


def _identify_boundary_vertices(T: np.ndarray) -> set[int]:
    """Return the set of vertex indices that lie on the fault's
    boundary (edges shared by exactly 1 triangle).  Smoothing must
    NOT move these so the fault's outline is preserved."""
    edge_count: dict[tuple[int, int], int] = {}
    for tri in T:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            edge_count[key] = edge_count.get(key, 0) + 1
    boundary: set[int] = set()
    for (u, v), n in edge_count.items():
        if n == 1:
            boundary.add(u)
            boundary.add(v)
    return boundary


def _coalesce_polylines(
    polylines: list[fi.Polyline], tol_m: float,
) -> list[fi.Polyline]:
    """Greedy join of polylines whose end-to-start (or end-to-end)
    points are within ``tol_m`` of each other in 3-D.

    Phase 1's ``chain_segments`` already does this at its ``snap_m``
    tolerance (1 cm by default).  This helper picks up cases where
    the chain just-barely-failed to merge because endpoints were a
    few cm apart.  Without it, the resulting two-polylines case
    breaks Phase 2 when both polylines fall on the same parent
    triangle's edge: the two endpoint nodes are distinct in V_list
    but project to the same 2-D point on the parent, which
    `triangle.triangulate` then either collapses or fails on.

    The merge concatenates two polylines into one:
        [p_0, p_1, ..., p_n]  +  [q_0, q_1, ..., q_m]
    when ``|p_n - q_0| < tol_m``.  Pierces are concatenated to match;
    the meeting vertex in the merged polyline keeps the position of
    p_n (q_0 is dropped as a near-duplicate).
    """
    if len(polylines) <= 1 or tol_m <= 0.0:
        return list(polylines)
    tol_sq = tol_m * tol_m

    def dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
        return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

    remaining = list(polylines)
    out: list[fi.Polyline] = []
    while remaining:
        cur = remaining.pop(0)
        # Try to extend cur by chaining an adjacent polyline.
        merged = True
        while merged:
            merged = False
            for j, other in enumerate(remaining):
                if cur.fault_a != other.fault_a or cur.fault_b != other.fault_b:
                    continue
                if cur.closed or other.closed:
                    continue
                # cur.end == other.start → append other after cur
                if dist2(cur.points[-1], other.points[0]) < tol_sq:
                    cur = fi.Polyline(
                        points=cur.points + other.points[1:],
                        closed=False,
                        fault_a=cur.fault_a,
                        fault_b=cur.fault_b,
                        pierces_a=list(cur.pierces_a) + list(other.pierces_a),
                        pierces_b=list(cur.pierces_b) + list(other.pierces_b),
                    )
                    remaining.pop(j)
                    merged = True
                    break
                # cur.end == other.end → reverse other and append
                if dist2(cur.points[-1], other.points[-1]) < tol_sq:
                    rev_pts = list(reversed(other.points))
                    rev_pa = list(reversed(other.pierces_a))
                    rev_pb = list(reversed(other.pierces_b))
                    cur = fi.Polyline(
                        points=cur.points + rev_pts[1:],
                        closed=False,
                        fault_a=cur.fault_a,
                        fault_b=cur.fault_b,
                        pierces_a=list(cur.pierces_a) + rev_pa,
                        pierces_b=list(cur.pierces_b) + rev_pb,
                    )
                    remaining.pop(j)
                    merged = True
                    break
                # cur.start == other.end → prepend other before cur
                if dist2(cur.points[0], other.points[-1]) < tol_sq:
                    cur = fi.Polyline(
                        points=other.points + cur.points[1:],
                        closed=False,
                        fault_a=cur.fault_a,
                        fault_b=cur.fault_b,
                        pierces_a=list(other.pierces_a) + list(cur.pierces_a),
                        pierces_b=list(other.pierces_b) + list(cur.pierces_b),
                    )
                    remaining.pop(j)
                    merged = True
                    break
                # cur.start == other.start → reverse other and prepend
                if dist2(cur.points[0], other.points[0]) < tol_sq:
                    rev_pts = list(reversed(other.points))
                    rev_pa = list(reversed(other.pierces_a))
                    rev_pb = list(reversed(other.pierces_b))
                    cur = fi.Polyline(
                        points=rev_pts + cur.points[1:],
                        closed=False,
                        fault_a=cur.fault_a,
                        fault_b=cur.fault_b,
                        pierces_a=rev_pa + list(cur.pierces_a),
                        pierces_b=rev_pb + list(cur.pierces_b),
                    )
                    remaining.pop(j)
                    merged = True
                    break
        out.append(cur)
    return out


# ---------------------------------------------------------------------------
# Validation gates
# ---------------------------------------------------------------------------
def is_manifold_no_t_junctions(
    V: np.ndarray, T: np.ndarray,
    *, tol_m: float = 1e-3,
) -> tuple[bool, list[tuple[int, int]]]:
    """Manifold gate per fault.  Returns (passed, offending_edges).

    Check 1: every edge is shared by ≤ 2 triangles within this fault.
    Check 2: no vertex sits in the open interior of any other edge to
             within ``tol_m`` (T-junction detection).

    See also :func:`is_manifold_split` for the same checks reported
    separately so the caller can treat 3+-fold edges as fatal while
    tolerating a few T-junctions.
    """
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    edge_count: dict[tuple[int, int], int] = {}
    for tri in T:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            edge_count[key] = edge_count.get(key, 0) + 1

    bad: list[tuple[int, int]] = []
    for key, n in edge_count.items():
        if n > 2:
            bad.append(key)
    if bad:
        return False, bad

    # T-junction check: for every edge (u, v), check that no OTHER
    # vertex (not u, not v) lies on the open segment within tol_m.
    # Vectorised over edges for speed.
    edges = np.array(list(edge_count.keys()), dtype=np.int64)   # (E, 2)
    A_pts = V[edges[:, 0]]                          # (E, 3)
    B_pts = V[edges[:, 1]]
    edge_vec = B_pts - A_pts
    edge_len_sq = (edge_vec * edge_vec).sum(axis=1).clip(min=1e-30)

    # Iterate over candidate vertices.  At SAFS scale (~1500 verts per
    # fault, ~3000 edges), this is ~5M ops — sub-second.
    for vi in range(V.shape[0]):
        p = V[vi]
        d_p = p[None, :] - A_pts                    # (E, 3)
        t = (d_p * edge_vec).sum(axis=1) / edge_len_sq
        # Strictly interior to the edge (not at endpoints).
        interior_mask = (t > 1e-6) & (t < 1.0 - 1e-6)
        # Skip edges that involve this vertex.
        interior_mask &= (edges[:, 0] != vi) & (edges[:, 1] != vi)
        if not interior_mask.any():
            continue
        proj = A_pts + t[:, None] * edge_vec
        d_perp = np.linalg.norm(p[None, :] - proj, axis=1)
        viol = interior_mask & (d_perp < tol_m)
        if viol.any():
            for ei in np.where(viol)[0]:
                bad.append((int(edges[ei, 0]), int(edges[ei, 1])))
    if bad:
        return False, bad
    return True, []


def is_manifold_split(
    V: np.ndarray, T: np.ndarray,
    *, tol_m: float = 1e-3,
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Same checks as :func:`is_manifold_no_t_junctions` but reports
    the two failure modes separately.

    Returns ``(non_manifold_edges, t_junction_edges)``:
        - ``non_manifold_edges``: edges shared by 3 or more triangles
          within the fault.  These are FATAL — HXT will reject the
          input outright.
        - ``t_junction_edges``: edges that have a vertex sitting in
          their open interior to within ``tol_m`` but are otherwise
          manifold.  These are recoverable — HXT may produce slightly
          worse local quality but typically meshes through them with
          the OptimizeNetgen post-pass.
    """
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    edge_count: dict[tuple[int, int], int] = {}
    for tri in T:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            edge_count[key] = edge_count.get(key, 0) + 1
    non_manifold = [k for k, n in edge_count.items() if n > 2]

    t_junctions: list[tuple[int, int]] = []
    if edge_count:
        edges = np.array(list(edge_count.keys()), dtype=np.int64)
        A_pts = V[edges[:, 0]]
        B_pts = V[edges[:, 1]]
        edge_vec = B_pts - A_pts
        edge_len_sq = (edge_vec * edge_vec).sum(axis=1).clip(min=1e-30)
        for vi in range(V.shape[0]):
            p = V[vi]
            d_p = p[None, :] - A_pts
            t = (d_p * edge_vec).sum(axis=1) / edge_len_sq
            interior_mask = (t > 1e-6) & (t < 1.0 - 1e-6)
            interior_mask &= (edges[:, 0] != vi) & (edges[:, 1] != vi)
            if not interior_mask.any():
                continue
            proj = A_pts + t[:, None] * edge_vec
            d_perp = np.linalg.norm(p[None, :] - proj, axis=1)
            viol = interior_mask & (d_perp < tol_m)
            if viol.any():
                for ei in np.where(viol)[0]:
                    t_junctions.append((int(edges[ei, 0]),
                                          int(edges[ei, 1])))
    return non_manifold, t_junctions


def verify_polyline_in_mesh(
    V: np.ndarray, T: np.ndarray,
    polyline: fi.Polyline,
    *, match_tol_m: float = 2e-2,
) -> bool:
    """Polyline edge-coincidence gate.  Returns True iff every
    consecutive pair of polyline points appears as an edge in the
    mesh (V, T) — using a vertex-match tolerance of ``match_tol_m``."""
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    if V.shape[0] == 0 or T.shape[0] == 0:
        return False
    pts = np.asarray(polyline.points, dtype=np.float64)
    # Map each polyline point to a mesh vertex by nearest-vertex within
    # match_tol_m.
    polyline_to_v: list[int] = []
    for p in pts:
        d = np.linalg.norm(V - p, axis=1)
        idx = int(np.argmin(d))
        if d[idx] > match_tol_m:
            return False
        polyline_to_v.append(idx)
    # Build edge set.
    edges: set[tuple[int, int]] = set()
    for tri in T:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((min(u, v), max(u, v)))
    # Check every consecutive pair.
    for i in range(len(polyline_to_v) - 1):
        u, v = polyline_to_v[i], polyline_to_v[i + 1]
        if u == v:
            continue                              # collapsed under snap; tolerated
        key = (min(u, v), max(u, v))
        if key not in edges:
            return False
    return True


# ---------------------------------------------------------------------------
# STL writer (matching ts_to_stl.py's ASCII format with full-double
# coordinates so HXT's PLC recovery sees no float32 truncation).
# ---------------------------------------------------------------------------
def _write_ascii_stl(
    path: Path, V: np.ndarray, T: np.ndarray, name: str,
) -> None:
    V = np.asarray(V, dtype=np.float64)
    T = np.asarray(T, dtype=np.int64)
    if T.shape[0] == 0:
        raise ValueError(f"refusing to write empty STL: {path}")
    a = V[T[:, 0]]
    b = V[T[:, 1]]
    c = V[T[:, 2]]
    n = np.cross(b - a, c - a)
    norms = np.linalg.norm(n, axis=1, keepdims=True)
    norms = np.where(norms > 0, norms, 1.0)
    n = n / norms
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%+.17e"
    with open(path, "w") as fh:
        fh.write(f"solid SAFS:{name}:conformal\n")
        for i in range(T.shape[0]):
            nx, ny, nz = n[i]
            ax, ay, az = a[i]
            bx, by, bz = b[i]
            cx, cy, cz = c[i]
            fh.write(
                f"  facet normal {nx:{'+.17e'}} {ny:{'+.17e'}} {nz:{'+.17e'}}\n"
                "    outer loop\n"
                f"      vertex {ax:{'+.17e'}} {ay:{'+.17e'}} {az:{'+.17e'}}\n"
                f"      vertex {bx:{'+.17e'}} {by:{'+.17e'}} {bz:{'+.17e'}}\n"
                f"      vertex {cx:{'+.17e'}} {cy:{'+.17e'}} {cz:{'+.17e'}}\n"
                "    endloop\n  endfacet\n"
            )
        fh.write(f"endsolid SAFS:{name}:conformal\n")


# ---------------------------------------------------------------------------
# Public driver: conformalize
# ---------------------------------------------------------------------------
def conformalize(
    in_stl_dir: Path | str,
    out_stl_dir: Path | str,
    included: list[str],
    *,
    snap_m: float = 1e-2,
    drop_short_segment_frac: float = 0.0,
    clearance_m: float = 0.0,
    target_edge_length_m: float | None = 1000.0,
    min_pierce_separation_m: float | None = None,
    cdt_min_angle_deg: float | None = None,
    smoothing_iters: int = 5,
    smoothing_relaxation: float = 0.5,
    flip_passes: int = 8,
) -> dict:
    """Phase 2 driver: read per-fault STLs, run Phase 1 + Phase 2 on
    every unordered fault pair, write conformal per-fault STLs +
    triangle_to_fault.json + intersection_report.json.

    Returns the report dict (intersection_report.json contents).

    Mesh-quality knobs (defaults tuned for 2000m CFM input + 1 km
    target on-fault edge):

    target_edge_length_m
        Knob C: pre-refine each fault to ~uniform target edge length
        BEFORE the intersection scan.  Implemented via lock-step 1-to-4
        midpoint subdivision (``fault_intersect.refine_uniform_midpoint``)
        so every shared edge has identical midpoints on both sides → no
        T-junctions.  Pass ``None`` to skip global refinement and keep
        the original CFM triangulation.  Default 1000.0 m.
    min_pierce_separation_m
        Knob A: merge consecutive polyline pierces closer than this
        (``fault_intersect.simplify_polyline``).  Default ``None`` →
        ``0.5 × target_edge_length_m`` if knob C is on, else ``0.0``
        (no simplification).  Pass ``0.0`` to disable explicitly.
    cdt_min_angle_deg
        Knob B (advanced, default off): post-conformalize quality CDT
        flag.  ``None`` (default) → use bare ``'p'`` constrained
        Delaunay (preserves all polyline edges exactly).  Set to e.g.
        ``30.0`` only if you accept that the triangle library may
        relax some constraint preservation to satisfy the angle bound;
        most users should leave this off and rely on knob C
        (pre-refinement) for quality control.
    """
    in_stl_dir = Path(in_stl_dir)
    out_stl_dir = Path(out_stl_dir)
    out_stl_dir.mkdir(parents=True, exist_ok=True)
    if snap_m <= 0.0:
        raise ValueError("snap_m must be > 0")
    if drop_short_segment_frac < 0.0:
        raise ValueError("drop_short_segment_frac must be >= 0")
    if target_edge_length_m is not None and target_edge_length_m <= 0.0:
        raise ValueError("target_edge_length_m must be > 0 (or None)")
    if cdt_min_angle_deg is not None and not (0.0 < cdt_min_angle_deg < 34.0):
        raise ValueError("cdt_min_angle_deg must be in (0, 34) (or None)")
    if min_pierce_separation_m is None:
        # Default: knob A off (no within-polyline simplification).
        # The midpoint refinement (knob C) keeps parents at uniform
        # ~target_edge_length, so even short polyline segments
        # produce only modest slivers (aspect ratio bounded by
        # parent_edge / shortest_polyline_seg).  Knob A introduces
        # cross-fault polyline-attribution complexity that interacts
        # badly with the CDT — leave it off by default; users who
        # want it can pass an explicit value.
        min_pierce_separation_m = 0.0
    if min_pierce_separation_m < 0.0:
        raise ValueError("min_pierce_separation_m must be >= 0")

    # Read inputs.  STL stores one vertex slot per triangle corner
    # (no shared vertices in the file), so we dedup after reading to
    # rebuild connectivity.  Without this, the manifold and polyline-
    # edge-coincidence gates would never see any shared edges.
    fault_meshes: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    dedup_snap_m = min(snap_m, 1e-3)            # tighter than chain-snap so we
                                                  # don't accidentally merge
                                                  # distinct CFM vertices
    for short in included:
        p = in_stl_dir / f"{short}.stl"
        if not p.exists():
            raise FileNotFoundError(f"missing STL for {short}: {p}")
        V_raw, T_raw = fi._read_stl_arrays(p)
        V_d, T_d = fi.dedup_mesh(V_raw, T_raw, snap_m=dedup_snap_m)
        # Knob C — uniform pre-refinement of the fault surface to the
        # target edge length BEFORE we look for intersections.  Without
        # this, far-from-polyline parts of the mesh stay at the
        # original (~1500 m) CFM resolution while near-polyline gets
        # finer, producing a non-uniform mesh.  Midpoint subdivision is
        # consistent across neighbours by construction so it cannot
        # introduce T-junctions.
        if target_edge_length_m is not None:
            V_d, T_d = fi.refine_uniform_midpoint(
                V_d, T_d, target_edge_length_m=target_edge_length_m,
            )
        fault_meshes[short] = (V_d, T_d)

    polylines_by_pair: dict[tuple[str, str], list[fi.Polyline]] = {}
    report: dict = {
        "schema_version": 1,
        "snap_m": snap_m,
        "drop_short_segment_frac": drop_short_segment_frac,
        "target_edge_length_m": target_edge_length_m,
        "min_pierce_separation_m": min_pierce_separation_m,
        "cdt_min_angle_deg": cdt_min_angle_deg,
        "pairs": {},
        "per_fault": {},
    }

    # PLAN Step (3) — re-scan between cascade passes.  We re-scan each
    # pair against the *current* (possibly already-split) fault meshes
    # so that pierce indices reference the live triangulation.
    pair_indices: list[tuple[int, int]] = [
        (i, j) for i in range(len(included))
                for j in range(i + 1, len(included))
    ]

    for ai, bi in pair_indices:
        sa = included[ai]
        sb = included[bi]
        V_A, T_A = fault_meshes[sa]
        V_B, T_B = fault_meshes[sb]

        segs = fi.tri_tri_intersect_3d(
            V_A, T_A, V_B, T_B, sa, sb, clearance_m=clearance_m,
        )
        if not segs:
            continue
        polylines = fi.chain_segments(segs, snap_m=snap_m)
        if not polylines:
            continue

        # Coalesce polylines whose endpoints are close (sub-cm to
        # cm).  Phase 1's chain_segments uses snap_m=1cm; CFM-data
        # near-tangent intersections at fault boundaries can produce
        # two separate polyline pieces with endpoints just outside
        # that tolerance.  Without this step, both pieces' shared
        # endpoint becomes two distinct global nodes that sit on the
        # same triangle edge in Phase 2 → degenerate sub-triangles or
        # T-junctions.  Merge tolerance is half the target edge
        # length (or 5 cm if knob C is off) — far below CFM minimum
        # spacing.
        if len(polylines) > 1:
            coalesce_tol = (
                0.5 * target_edge_length_m
                if target_edge_length_m is not None else 5e-2
            )
            polylines = _coalesce_polylines(polylines, coalesce_tol)

        # Knob A: simplify each polyline by merging consecutive pierces
        # closer than ``min_pierce_separation_m``.  Both faults share
        # the same simplified polyline → conformality preserved.
        if min_pierce_separation_m > 0.0:
            polylines = [
                fi.simplify_polyline(
                    pl,
                    min_pierce_separation_m=min_pierce_separation_m,
                )
                for pl in polylines
            ]

        polylines_by_pair[(sa, sb)] = polylines

        # Apply Step (1) + (2) to fault A and fault B with the
        # quality CDT (knob B) and uniform refinement (knob C).
        # Knob C is applied UPSTREAM via the midpoint pre-refinement
        # (above), so we do NOT pass target_edge_length_m to the
        # per-parent CDT here.  Knob B (cdt_min_angle_deg) is opt-in;
        # default off because the triangle library can't satisfy
        # quality flags while preserving polyline edges (q+a tend to
        # subdivide segments).
        V_A_new, T_A_new, pl_to_v_A, dropped_A = _conformalize_one_fault(
            V_A, T_A, polylines,
            use_pierces_a=True,
            snap_m=snap_m,
            drop_short_segment_frac=drop_short_segment_frac,
            cdt_min_angle_deg=cdt_min_angle_deg,
            target_edge_length_m=None,
        )
        V_B_new, T_B_new, pl_to_v_B, dropped_B = _conformalize_one_fault(
            V_B, T_B, polylines,
            use_pierces_a=False,
            snap_m=snap_m,
            drop_short_segment_frac=drop_short_segment_frac,
            cdt_min_angle_deg=cdt_min_angle_deg,
            target_edge_length_m=None,
        )
        fault_meshes[sa] = (V_A_new, T_A_new)
        fault_meshes[sb] = (V_B_new, T_B_new)

        # PLAN §Validation gates — non-manifold (3+ fold edge) is FATAL
        # because HXT will reject; T-junctions are tolerated as
        # warnings because HXT can usually mesh through them with the
        # OptimizeNetgen post-pass, and on real CFM data a few are
        # introduced by cross-parent projection round-off.
        nm_A, tj_A = is_manifold_split(V_A_new, T_A_new)
        nm_B, tj_B = is_manifold_split(V_B_new, T_B_new)
        if nm_A:
            raise RuntimeError(
                f"non-manifold edges (3+ fold) on fault {sa} after pair "
                f"({sa} × {sb}): {len(nm_A)} edges; first 5: {nm_A[:5]}"
            )
        if nm_B:
            raise RuntimeError(
                f"non-manifold edges (3+ fold) on fault {sb} after pair "
                f"({sa} × {sb}): {len(nm_B)} edges; first 5: {nm_B[:5]}"
            )
        if tj_A:
            _log.warning(
                "fault %s has %d T-junction edge(s) after pair (%s × %s); "
                "first 5: %s.  Will pass to HXT — OptimizeNetgen typically "
                "tolerates these.",
                sa, len(tj_A), sa, sb, tj_A[:5],
            )
        if tj_B:
            _log.warning(
                "fault %s has %d T-junction edge(s) after pair (%s × %s); "
                "first 5: %s.  Will pass to HXT — OptimizeNetgen typically "
                "tolerates these.",
                sb, len(tj_B), sa, sb, tj_B[:5],
            )
        for pl_idx, pl in enumerate(polylines):
            if not verify_polyline_in_mesh(V_A_new, T_A_new, pl):
                raise RuntimeError(
                    f"polyline edge-coincidence gate FAILED for pair "
                    f"({sa} × {sb}) polyline #{pl_idx} on fault {sa}"
                )
            if not verify_polyline_in_mesh(V_B_new, T_B_new, pl):
                raise RuntimeError(
                    f"polyline edge-coincidence gate FAILED for pair "
                    f"({sa} × {sb}) polyline #{pl_idx} on fault {sb}"
                )
        residual = fi.tri_tri_intersect_3d_interior_only(
            V_A_new, T_A_new, V_B_new, T_B_new, sa, sb,
            clearance_m=clearance_m,
        )
        if residual:
            raise RuntimeError(
                f"interior-crossing gate FAILED for pair ({sa} × {sb}): "
                f"{len(residual)} residual interior crossings remain "
                f"after conformalization; first: {residual[0]}"
            )

        # Constrained-Delaunay edge-flip pass: rewire pairs of
        # triangles whose shared (non-constraint) edge is the
        # "wrong diagonal" — replaces sliver pairs with better-shaped
        # pairs.  Polyline edges and fault-boundary edges are NEVER
        # flipped, so cross-fault polyline edge-coincidence is
        # preserved.  Run BEFORE smoothing — flipping fixes the
        # topology, smoothing then redistributes vertex positions on
        # the now-better-connected mesh.
        if flip_passes > 0:
            locked_edges_A = _identify_locked_edges(T_A_new, pl_to_v_A)
            locked_edges_B = _identify_locked_edges(T_B_new, pl_to_v_B)
            T_A_new = _constrained_delaunay_flips(
                V_A_new, T_A_new, locked_edges_A, max_passes=flip_passes,
            )
            T_B_new = _constrained_delaunay_flips(
                V_B_new, T_B_new, locked_edges_B, max_passes=flip_passes,
            )

        # Surface-mesh smoothing pass: relax INTERIOR vertices toward
        # their 1-ring centroid to even out triangle shapes near the
        # polyline.  Locked vertices: (1) polyline vertices (their
        # node IDs in pl_to_v_A / pl_to_v_B) — moving them would break
        # cross-fault polyline edge-coincidence; (2) fault-boundary
        # vertices — moving them would shrink the fault outline.
        if smoothing_iters > 0:
            locked_A: set[int] = set()
            for nodes in pl_to_v_A:
                locked_A.update(nodes)
            locked_A.update(_identify_boundary_vertices(T_A_new))
            locked_B: set[int] = set()
            for nodes in pl_to_v_B:
                locked_B.update(nodes)
            locked_B.update(_identify_boundary_vertices(T_B_new))
            V_A_new = _smooth_fault_surface(
                V_A_new, T_A_new, locked_A,
                n_iters=smoothing_iters,
                relaxation=smoothing_relaxation,
            )
            V_B_new = _smooth_fault_surface(
                V_B_new, T_B_new, locked_B,
                n_iters=smoothing_iters,
                relaxation=smoothing_relaxation,
            )
            fault_meshes[sa] = (V_A_new, T_A_new)
            fault_meshes[sb] = (V_B_new, T_B_new)

        report["pairs"][f"{sa}__x__{sb}"] = {
            "n_polylines": len(polylines),
            "n_polyline_vertices": sum(len(p.points) for p in polylines),
            "pre_split_n_tri_A": int(T_A.shape[0]),
            "post_split_n_tri_A": int(T_A_new.shape[0]),
            "pre_split_n_tri_B": int(T_B.shape[0]),
            "post_split_n_tri_B": int(T_B_new.shape[0]),
            "n_dropped_short_A": dropped_A,
            "n_dropped_short_B": dropped_B,
            "smoothing_iters": smoothing_iters,
            "gates": {
                "manifold_A": "PASS" if not tj_A else f"PASS_WITH_{len(tj_A)}_T_JUNCTIONS",
                "manifold_B": "PASS" if not tj_B else f"PASS_WITH_{len(tj_B)}_T_JUNCTIONS",
                "polyline_edge_coincidence": "PASS",
                "interior_crossing_only": "PASS",
            },
        }

    # Write per-fault conformal STLs.
    range_record: dict[str, dict] = {}
    cursor = 0
    for short in included:
        V, T = fault_meshes[short]
        out_path = out_stl_dir / f"{short}.stl"
        _write_ascii_stl(out_path, V, T, short)
        range_record[short] = {
            "n_triangles": int(T.shape[0]),
            "range": [cursor, cursor + int(T.shape[0])],
        }
        cursor += int(T.shape[0])
        report["per_fault"][short] = {
            "n_vertices": int(V.shape[0]),
            "n_triangles": int(T.shape[0]),
            "stl_path": str(out_path),
        }

    # triangle_to_fault.json
    tri_to_fault = {
        "schema_version": 1,
        "faults": range_record,
        "n_total_triangles": cursor,
    }
    (out_stl_dir / "triangle_to_fault.json").write_text(
        json.dumps(tri_to_fault, indent=2) + "\n"
    )

    # intersection_report.json
    (out_stl_dir / "intersection_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase 2 driver: insert cross-fault polylines as "
                    "exact mesh edges on each fault.",
    )
    parser.add_argument("--in-stl-dir", required=True, type=Path)
    parser.add_argument("--out-stl-dir", required=True, type=Path)
    parser.add_argument("--include-fault", action="append", required=True,
                        metavar="SHORT_NAME")
    parser.add_argument("--snap-m", type=float, default=1e-2)
    parser.add_argument("--drop-short-segment-frac", type=float, default=0.0)
    parser.add_argument("--clearance-m", type=float, default=0.0)
    parser.add_argument("--target-edge-length-m", type=float, default=1000.0,
                        help="Knob C: pre-refine each fault to ~uniform "
                             "target edge length via 1-to-4 midpoint "
                             "subdivision before Phase 1 runs.  Pass 0 to "
                             "skip global refinement.  Default 1000 m.")
    parser.add_argument("--min-pierce-separation-m", type=float, default=-1.0,
                        help="Knob A: merge consecutive polyline pierces "
                             "closer than this (only when both adjacent "
                             "segments are in the same parent triangle "
                             "pair).  -1 (default) → 0.5 × target_edge "
                             "if knob C is on, else 0.")
    parser.add_argument("--cdt-min-angle-deg", type=float, default=-1.0,
                        help="Knob B: minimum interior angle for the per-"
                             "parent CDT.  -1 (default) → off.  Setting "
                             "may break cross-fault polyline edge-"
                             "coincidence due to triangle-library "
                             "constraint-vs-quality limitations.")
    parser.add_argument("--smoothing-iters", type=int, default=5,
                        help="Laplacian-smoothing iterations on each "
                             "post-conformal fault surface.  Default 5; "
                             "0 to disable.")
    parser.add_argument("--smoothing-relaxation", type=float, default=0.5,
                        help="Per-iteration relaxation factor for the "
                             "Laplacian smoother.  Default 0.5.")
    parser.add_argument("--flip-passes", type=int, default=8,
                        help="Constrained-Delaunay edge-flip passes after "
                             "split + propagate.  Default 8; 0 to disable.")
    args = parser.parse_args(argv)

    target_edge = (args.target_edge_length_m
                   if args.target_edge_length_m > 0.0 else None)
    min_pierce = (args.min_pierce_separation_m
                  if args.min_pierce_separation_m >= 0.0 else None)
    cdt_angle = (args.cdt_min_angle_deg
                 if args.cdt_min_angle_deg > 0.0 else None)
    report = conformalize(
        args.in_stl_dir, args.out_stl_dir, args.include_fault,
        snap_m=args.snap_m,
        drop_short_segment_frac=args.drop_short_segment_frac,
        clearance_m=args.clearance_m,
        target_edge_length_m=target_edge,
        min_pierce_separation_m=min_pierce,
        cdt_min_angle_deg=cdt_angle,
        smoothing_iters=args.smoothing_iters,
        smoothing_relaxation=args.smoothing_relaxation,
        flip_passes=args.flip_passes,
    )
    print(json.dumps(report["per_fault"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
