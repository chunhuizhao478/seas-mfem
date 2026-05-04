#!/usr/bin/env python
"""Refine the per-fault SURFACE triangulation in a band around every
cross-fault polyline vertex, by subdividing fault triangles 1→4
(Loop-style midpoint subdivision).

This operates on the per-fault STL output of `corefine_faults` and
produces refined STLs ready for `generate_safs_mesh.py`.  It does NOT
touch the bulk volume mesh — that's HXT's job downstream.

SCOPE LIMIT (read this before assuming the script will fix slivers):
This script changes per-fault triangle DENSITY in the polyline-band
region.  It does NOT change the cross-fault polyline TOPOLOGY (the
sequence and length of shared edges between two faults).  Slivers in
the SAFS dataset are 4-fault-vertex wedge tets driven by SHALLOW
DIHEDRAL ANGLES at polyline edges (REVIEW_sliver_classification.md
R-001, REVIEW_intersection_refinement_investigation.md R-007).  The
wedge configuration is determined by the polyline itself —
adjacent-triangle density does not change it.  To kill those slivers,
use `break_fault_wedges.py`, which inserts a midpoint vertex ON the
polyline edge, splitting each wedge tet into two 5-vertex
non-coplanar tets.

This script is useful for sliver kinds that depend on
adjacent-triangle density (e.g., needle tets where one vertex is far
from the polyline) — those are NOT the dominant kind in the SAFS
dataset.

Why subdivide here, not via gmsh size fields:
  - The polyline-adjacent bulk volume is ALREADY at fault-triangulation-
    driven density (HXT honors fault edges as hard constraints).
  - The slivers come from 4-fault-vertex wedge tets where the fault
    triangulation imposes 4 near-coplanar vertices.  More vertices on
    the fault surface in the wedge zone means HXT has MORE choices for
    apex vertices (c, d) — some of which are off the near-coplanar
    line and produce well-shaped tets.
  - Subdividing 1→4 with midpoint vertices keeps the surface
    piecewise-flat on each original triangle (the new vertices are
    coplanar with the original) so the SURFACE GEOMETRY is preserved
    bit-identically.

Cross-fault conformity invariant:
  - Polyline edges (= edges shared between fault A and fault B's
    triangulations at coord-precision `snap_m`) MUST be split in BOTH
    faults consistently, with the SAME midpoint coord, otherwise HXT
    sees mismatched polyline geometry and rejects the PLC.
  - Implementation: collect ALL edges that need to be split across all
    faults FIRST (union over per-fault subdivision plans), then apply
    the splits using a SHARED edge-midpoint table.

Algorithm:
  1. Read per-fault STLs.
  2. Identify polyline-endpoint vertex coords (≥2 fault appearances).
  3. For each fault triangle whose centroid is within `band_radius`
     of any polyline-endpoint vertex, MARK FOR SUBDIVISION.
  4. Collect all unique edges that need to be split (across all
     marked triangles in all faults).
  5. For each shared edge, compute midpoint ONCE.  When fault A and
     fault B share the edge (it's a polyline edge), both faults use
     the same midpoint coord.
  6. Apply 1→4 Loop subdivision on every marked triangle, using the
     shared edge midpoints.  Adjacent UNMARKED triangles whose edge
     was split need a 1→2 split to maintain manifoldness (T-junction
     fix).
  7. Write modified STLs.

Usage:
  python refine_fault_near_intersections.py \\
      --in-stl-dir output/.../stl_conformal \\
      --out-stl-dir output/.../stl_refined \\
      --include-fault A --include-fault B [...] \\
      [--band-radius-m 1500] [--snap-m 0.1]
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np


# ---------------------------------------------------------------------------
# STL I/O — same precision as the rest of the pipeline.
# ---------------------------------------------------------------------------
def _read_ascii_stl(p: Path) -> tuple[np.ndarray, np.ndarray, str]:
    verts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    cur: list[tuple[float, float, float]] = []
    solid_name = ""
    with p.open() as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            if s[0] == "solid" and len(s) > 1:
                solid_name = " ".join(s[1:])
            elif s[0] == "vertex":
                cur.append((float(s[1]), float(s[2]), float(s[3])))
            elif s[0] == "endloop":
                if len(cur) != 3:
                    raise ValueError(f"{p}: outer loop with {len(cur)} vertices")
                a = len(verts); b = a + 1; c = a + 2
                verts.extend(cur)
                tris.append((a, b, c))
                cur = []
    return (np.asarray(verts, dtype=np.float64),
            np.asarray(tris, dtype=np.int64),
            solid_name)


def _dedup_vertices(V: np.ndarray, T: np.ndarray, tol: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    if V.shape[0] == 0 or tol <= 0:
        return V, T
    rounded = np.round(V / tol).astype(np.int64)
    _, inverse = np.unique(rounded, axis=0, return_inverse=True)
    new_V = np.zeros((inverse.max() + 1, 3), dtype=np.float64)
    seen = np.zeros(inverse.max() + 1, dtype=bool)
    for i, idx in enumerate(inverse):
        if not seen[idx]:
            new_V[idx] = V[i]
            seen[idx] = True
    new_T = inverse[T]
    keep = (new_T[:, 0] != new_T[:, 1]) & \
           (new_T[:, 1] != new_T[:, 2]) & \
           (new_T[:, 0] != new_T[:, 2])
    return new_V, new_T[keep]


def _write_ascii_stl(p: Path, V: np.ndarray, T: np.ndarray,
                     solid_name: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        f.write(f"solid {solid_name}\n")
        for i in range(T.shape[0]):
            v0 = V[T[i, 0]]; v1 = V[T[i, 1]]; v2 = V[T[i, 2]]
            n = np.cross(v1 - v0, v2 - v0)
            nlen = np.linalg.norm(n)
            n = n / nlen if nlen > 0 else np.zeros(3)
            f.write(f"  facet normal {n[0]:+.17e} {n[1]:+.17e} {n[2]:+.17e}\n")
            f.write("    outer loop\n")
            for v in (v0, v1, v2):
                f.write(f"      vertex {v[0]:+.17e} {v[1]:+.17e} {v[2]:+.17e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")


# ---------------------------------------------------------------------------
# Snap-key helpers for cross-fault edge matching.
# ---------------------------------------------------------------------------
def _snap_key(p: np.ndarray, snap_m: float) -> tuple[int, int, int]:
    return (int(round(p[0] / snap_m)),
            int(round(p[1] / snap_m)),
            int(round(p[2] / snap_m)))


def _edge_snap_key(p1: np.ndarray, p2: np.ndarray, snap_m: float
                   ) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    k1 = _snap_key(p1, snap_m)
    k2 = _snap_key(p2, snap_m)
    return (k1, k2) if k1 < k2 else (k2, k1)


# ---------------------------------------------------------------------------
# Step 1: identify cross-fault polyline vertices (vertex coords appearing in
# ≥2 faults).  This includes EVERY vertex along EVERY cross-fault
# intersection polyline — both interior and curve-endpoints — not just
# curve endpoints.  The legacy name `collect_polyline_endpoints` is kept
# as an alias for backward compatibility.
# ---------------------------------------------------------------------------
def collect_polyline_vertices(meshes: dict[str, tuple[np.ndarray, np.ndarray]],
                               snap_m: float
                               ) -> set[tuple[int, int, int]]:
    key_to_faults: dict[tuple[int, int, int], set[str]] = {}
    for name, (V, _) in meshes.items():
        for v in V:
            key_to_faults.setdefault(_snap_key(v, snap_m), set()).add(name)
    return {k for k, fs in key_to_faults.items() if len(fs) >= 2}


# Backward-compat alias.  Prefer `collect_polyline_vertices` in new code.
collect_polyline_endpoints = collect_polyline_vertices


# ---------------------------------------------------------------------------
# Step 2: per-fault, mark triangles whose centroid is within band_radius of
# any polyline endpoint vertex coord.
# ---------------------------------------------------------------------------
def mark_triangles_in_band(V: np.ndarray, T: np.ndarray,
                            polyline_pts_xyz: np.ndarray,
                            polyline_keys: set[tuple[int, int, int]],
                            snap_m: float,
                            band_radius: float,
                            include_polyline_triangles: bool = True
                            ) -> np.ndarray:
    """Mark triangles whose centroid is within `band_radius` of any
    cross-fault polyline vertex.

    R-001 fix (REVIEW_intersection_refinement_investigation.md):
    polyline-bordering triangles ARE the triangles that bound the
    wedge configurations driving slivers.  Earlier revisions excluded
    them defensively against a cross-fault asymmetry concern that is
    already addressed by the shared `edge_midpoint_coord` map keyed on
    snap_key (subdivide_fault below uses this map so both faults that
    share a polyline edge get bit-identical midpoint coords).  Removing
    the exclusion is required for the refinement to reach the slivers.

    `include_polyline_triangles=False` restores the legacy behaviour
    (mark only triangles whose 3 vertices are interior to this fault)
    in case a future caller needs it for debugging.
    """
    if V.shape[0] == 0 or T.shape[0] == 0 or polyline_pts_xyz.shape[0] == 0:
        return np.zeros(T.shape[0], dtype=bool)
    centroids = V[T].mean(axis=1)
    marked = np.zeros(T.shape[0], dtype=bool)
    for p in polyline_pts_xyz:
        d = np.linalg.norm(centroids - p, axis=1)
        marked |= (d <= band_radius)
    if not include_polyline_triangles:
        for ti in range(T.shape[0]):
            if not marked[ti]:
                continue
            for vidx in T[ti]:
                if _snap_key(V[vidx], snap_m) in polyline_keys:
                    marked[ti] = False
                    break
    return marked


# ---------------------------------------------------------------------------
# Step 3 & 4: Loop-style 1→4 subdivision with cross-fault edge sharing.
#
# For each triangle (a, b, c) marked for subdivision:
#   - Compute midpoints m_ab, m_bc, m_ca (each is the average of two
#     existing vertex coords, hence ON the fault's piecewise plane).
#   - Replace triangle with 4 sub-triangles: (a, m_ab, m_ca),
#     (m_ab, b, m_bc), (m_ca, m_bc, c), (m_ab, m_bc, m_ca).
#
# Adjacent UNMARKED triangles whose edge was split need a 1→2 split
# to maintain manifoldness (T-junction fix).  For an unmarked triangle
# (p, q, r) where edge (p, q) was split at midpoint m_pq, replace with
# (p, m_pq, r) and (m_pq, q, r).
#
# Cross-fault edge sharing: when an edge (a, b) is a polyline edge
# (snap-key matches across two faults), both faults use the SAME
# midpoint coord — computed ONCE in the global edge_midpoint_coord map.
# ---------------------------------------------------------------------------
def subdivide_fault(V: np.ndarray, T: np.ndarray,
                     marked: np.ndarray,
                     edge_midpoint_coord: dict[tuple[tuple[int, int, int],
                                                      tuple[int, int, int]],
                                                np.ndarray],
                     edges_to_split: set[tuple[tuple[int, int, int],
                                                 tuple[int, int, int]]],
                     snap_m: float
                     ) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Apply 1→4 subdivision to every `marked` triangle and 1→2 T-junction
    fix to every unmarked triangle that has an edge in `edges_to_split`.

    Returns (new_V, new_T, n_loop_splits, n_t_junction_fixes).
    """
    # Vertex map: snap-key (of midpoint coord) -> vertex index in new_V.
    new_V: list[np.ndarray] = [v for v in V]
    midkey_to_vidx: dict[tuple[int, int, int], int] = {}
    new_T: list[list[int]] = []
    n_loop = 0
    n_tjunc = 0

    def _get_or_add_midpoint(p1: np.ndarray, p2: np.ndarray) -> int:
        ek = _edge_snap_key(p1, p2, snap_m)
        # Use the canonical midpoint coord from edge_midpoint_coord if
        # this edge is in the global split-set (cross-fault sharing).
        if ek in edge_midpoint_coord:
            mid = edge_midpoint_coord[ek]
        else:
            mid = (p1 + p2) / 2.0
        midkey = _snap_key(mid, snap_m)
        if midkey in midkey_to_vidx:
            return midkey_to_vidx[midkey]
        idx = len(new_V)
        new_V.append(mid)
        midkey_to_vidx[midkey] = idx
        return idx

    for ti in range(T.shape[0]):
        a, b, c = int(T[ti, 0]), int(T[ti, 1]), int(T[ti, 2])
        pa, pb, pc = V[a], V[b], V[c]

        if marked[ti]:
            # 1->4 Loop subdivision.
            m_ab = _get_or_add_midpoint(pa, pb)
            m_bc = _get_or_add_midpoint(pb, pc)
            m_ca = _get_or_add_midpoint(pc, pa)
            new_T.append([a,    m_ab, m_ca])
            new_T.append([m_ab, b,    m_bc])
            new_T.append([m_ca, m_bc, c])
            new_T.append([m_ab, m_bc, m_ca])
            n_loop += 1
            continue

        # Unmarked triangle: check if any edge was split (T-junction fix).
        ek_ab = _edge_snap_key(pa, pb, snap_m)
        ek_bc = _edge_snap_key(pb, pc, snap_m)
        ek_ca = _edge_snap_key(pc, pa, snap_m)
        split_ab = ek_ab in edges_to_split
        split_bc = ek_bc in edges_to_split
        split_ca = ek_ca in edges_to_split
        n_split = int(split_ab) + int(split_bc) + int(split_ca)

        if n_split == 0:
            new_T.append([a, b, c])
            continue

        if n_split == 1:
            # 1→2 split.
            if split_ab:
                m = _get_or_add_midpoint(pa, pb)
                new_T.append([a, m, c])
                new_T.append([m, b, c])
            elif split_bc:
                m = _get_or_add_midpoint(pb, pc)
                new_T.append([a, b, m])
                new_T.append([a, m, c])
            else:  # split_ca
                m = _get_or_add_midpoint(pc, pa)
                new_T.append([a, b, m])
                new_T.append([m, b, c])
            n_tjunc += 1
            continue

        if n_split == 2:
            # 1→3 split.  Two edges share an apex; the third vertex of the
            # triangle is the apex.  Identify the apex and form three
            # sub-tris like in break_fault_wedges.
            if split_ab and split_bc:
                # apex = b
                m1 = _get_or_add_midpoint(pa, pb)
                m2 = _get_or_add_midpoint(pb, pc)
                new_T.append([a,  m1, c])
                new_T.append([m1, b,  m2])
                new_T.append([m1, m2, c])
            elif split_bc and split_ca:
                # apex = c
                m1 = _get_or_add_midpoint(pb, pc)
                m2 = _get_or_add_midpoint(pc, pa)
                new_T.append([a, b,  m1])
                new_T.append([a, m1, m2])
                new_T.append([m2, m1, c])
            else:  # split_ab and split_ca
                # apex = a
                m1 = _get_or_add_midpoint(pa, pb)
                m2 = _get_or_add_midpoint(pc, pa)
                new_T.append([a,  m1, m2])
                new_T.append([m1, b,  c])
                new_T.append([m2, m1, c])
            n_tjunc += 1
            continue

        # n_split == 3: full 1→4 subdivision (treat as marked).
        m_ab = _get_or_add_midpoint(pa, pb)
        m_bc = _get_or_add_midpoint(pb, pc)
        m_ca = _get_or_add_midpoint(pc, pa)
        new_T.append([a,    m_ab, m_ca])
        new_T.append([m_ab, b,    m_bc])
        new_T.append([m_ca, m_bc, c])
        new_T.append([m_ab, m_bc, m_ca])
        n_tjunc += 1  # counts as T-junction fix even if topology is 1→4

    new_V_arr = np.asarray(new_V, dtype=np.float64)
    new_T_arr = np.asarray(new_T, dtype=np.int64)
    return new_V_arr, new_T_arr, n_loop, n_tjunc


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-stl-dir", required=True, type=Path)
    p.add_argument("--out-stl-dir", required=True, type=Path)
    p.add_argument("--include-fault", action="append", default=[],
                   metavar="SHORT_NAME")
    p.add_argument("--band-radius-m", type=float, default=1500.0,
                   help="Distance from polyline endpoints within which "
                        "fault triangles are flagged for 1->4 subdivision. "
                        "Default 1500 m.")
    p.add_argument("--snap-m", type=float, default=0.1,
                   help="Coord-precision for cross-fault edge matching. "
                        "Default 0.1 m.")
    p.add_argument("--legacy-exclude-polyline-triangles",
                   dest="include_polyline_triangles",
                   action="store_false", default=True,
                   help="Restore the pre-R-001 behaviour: triangles "
                        "with any polyline-vertex are excluded from "
                        "refinement.  Default (R-001 fix) is to "
                        "include them; cross-fault midpoint conformity "
                        "is preserved by the shared edge_midpoint_coord "
                        "map keyed on snap_key.")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    if not args.in_stl_dir.is_dir():
        raise SystemExit(f"--in-stl-dir not found: {args.in_stl_dir}")
    if not args.include_fault:
        raise SystemExit("--include-fault must appear at least once")
    args.out_stl_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: load and dedup.
    raw: dict[str, tuple[np.ndarray, np.ndarray, str]] = {}
    for short in args.include_fault:
        path = args.in_stl_dir / f"{short}.stl"
        if not path.exists():
            raise SystemExit(f"missing input STL: {path}")
        V, T, solid = _read_ascii_stl(path)
        V, T = _dedup_vertices(V, T, tol=args.snap_m)
        raw[short] = (V, T, solid)
        print(f"  load {short}: V={V.shape[0]}, T={T.shape[0]}",
              file=sys.stderr)

    # Step 2: collect cross-fault polyline vertex coords.
    pl_keys = collect_polyline_vertices(
        {n: (V, T) for n, (V, T, _) in raw.items()}, args.snap_m)
    pl_pts_xyz = np.asarray(
        [[k[0] * args.snap_m, k[1] * args.snap_m, k[2] * args.snap_m]
         for k in pl_keys], dtype=np.float64).reshape(-1, 3)
    print(f"  cross-fault polyline vertices (coords in ≥2 faults): "
          f"{pl_pts_xyz.shape[0]}", file=sys.stderr)

    # Step 3: per-fault, mark triangles in band.
    marks: dict[str, np.ndarray] = {}
    for short, (V, T, _) in raw.items():
        m = mark_triangles_in_band(V, T, pl_pts_xyz, pl_keys,
                                     args.snap_m, args.band_radius_m,
                                     include_polyline_triangles=
                                         args.include_polyline_triangles)
        marks[short] = m
        legacy_note = ("" if args.include_polyline_triangles
                       else " AND no polyline-vertex vertex")
        print(f"  mark {short}: {int(m.sum())} of {T.shape[0]} "
              f"triangles within {args.band_radius_m:.0f} m of any "
              f"cross-fault polyline vertex{legacy_note}",
              file=sys.stderr)

    # Step 4: enumerate edges to split.  An edge is split if EITHER:
    #   (a) it belongs to a marked triangle (its 3 edges all get split),
    # OR
    #   (b) it's the OTHER side of a polyline edge that's being split in
    #       the partner fault — this is automatic if both faults' edge
    #       maps overlap at the same snap-key.
    # We collect all edges-to-split GLOBALLY across faults, then each
    # fault's `subdivide_fault` will check its triangles' edges against
    # this global set.
    edges_to_split: set[tuple[tuple[int, int, int],
                                tuple[int, int, int]]] = set()
    for short, (V, T, _) in raw.items():
        m = marks[short]
        for ti in range(T.shape[0]):
            if not m[ti]:
                continue
            a, b, c = int(T[ti, 0]), int(T[ti, 1]), int(T[ti, 2])
            edges_to_split.add(_edge_snap_key(V[a], V[b], args.snap_m))
            edges_to_split.add(_edge_snap_key(V[b], V[c], args.snap_m))
            edges_to_split.add(_edge_snap_key(V[c], V[a], args.snap_m))
    print(f"  global edges-to-split: {len(edges_to_split)}",
          file=sys.stderr)

    # Step 5: compute canonical midpoint coords for every edge in
    # `edges_to_split`.  When an edge appears in multiple faults
    # (polyline edge), we compute the midpoint from FAULT A's coords
    # (the first fault that has the edge) — same canonicalization as
    # in break_fault_wedges.  Cross-fault sharing is preserved because
    # the snap-key is the same and we reuse the same midpoint coord.
    edge_midpoint_coord: dict[tuple[tuple[int, int, int],
                                     tuple[int, int, int]],
                               np.ndarray] = {}
    for short, (V, T, _) in raw.items():
        for ti in range(T.shape[0]):
            a, b, c = int(T[ti, 0]), int(T[ti, 1]), int(T[ti, 2])
            for p1, p2 in ((V[a], V[b]), (V[b], V[c]), (V[c], V[a])):
                ek = _edge_snap_key(p1, p2, args.snap_m)
                if ek in edges_to_split and ek not in edge_midpoint_coord:
                    edge_midpoint_coord[ek] = (p1 + p2) / 2.0

    # Step 6: subdivide each fault.
    n_loop_total = 0
    n_tjunc_total = 0
    out_meshes: dict[str, tuple[np.ndarray, np.ndarray, str]] = {}
    for short, (V, T, solid) in raw.items():
        new_V, new_T, n_loop, n_tjunc = subdivide_fault(
            V, T, marks[short], edge_midpoint_coord, edges_to_split,
            args.snap_m)
        n_loop_total += n_loop
        n_tjunc_total += n_tjunc
        out_meshes[short] = (new_V, new_T, solid)
        print(f"  refine {short}: V={V.shape[0]}->{new_V.shape[0]}, "
              f"T={T.shape[0]}->{new_T.shape[0]} "
              f"({n_loop} loop, {n_tjunc} t-junction)",
              file=sys.stderr)

    # Step 7: write modified STLs.
    for short, (V, T, solid) in out_meshes.items():
        out_path = args.out_stl_dir / f"{short}.stl"
        _write_ascii_stl(out_path, V, T,
                          solid or f"SAFS:{short}:refined")

    # Report.
    report_path = args.report_json or (
        args.out_stl_dir / "refine_report.json")
    with report_path.open("w") as fh:
        json.dump({
            "snap_m": args.snap_m,
            "band_radius_m": args.band_radius_m,
            "include_polyline_triangles": bool(
                args.include_polyline_triangles),
            "n_polyline_vertices": int(pl_pts_xyz.shape[0]),
            "n_polyline_endpoints": int(pl_pts_xyz.shape[0]),  # legacy alias
            "n_global_edges_to_split": len(edges_to_split),
            "n_loop_subdivisions": n_loop_total,
            "n_t_junction_fixes": n_tjunc_total,
            "per_fault": {
                short: {
                    "V_in":  int(raw[short][0].shape[0]),
                    "V_out": int(out_meshes[short][0].shape[0]),
                    "T_in":  int(raw[short][1].shape[0]),
                    "T_out": int(out_meshes[short][1].shape[0]),
                }
                for short in raw
            },
        }, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
