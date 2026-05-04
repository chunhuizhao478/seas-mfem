#!/usr/bin/env python
"""Pre-HXT helper: split shallow-dihedral cross-fault wedge edges to
break the 4-fault-vertex wedge tets that drive γ_min.

Diagnosis (2026-04-30):
  - Step 5b's worst sliver tet has 4 fault vertices forming a wedge
    `{a, b, c, d}` where (a, b) is a polyline edge shared between fault
    A and fault B, c is the apex of A's triangle, d is the apex of
    B's triangle.  The dihedral between A's plane and B's plane at
    this edge is shallow (often < 10°).
  - Bulk refinement cannot reach this tet (3 empirical confirmations).
  - The fix is to MODIFY the fault triangulation: insert a midpoint
    `m = (a + b) / 2` as a new fault vertex on the polyline, splitting
    BOTH fault triangles meeting at (a, b).  After the split, HXT's
    wedge candidates become `{a, m, c, d}` and `{m, b, c, d}` — each
    with 5 distinct vertices and finite dihedral; no longer near-
    coplanar slivers.

Algorithm:
  1. Read all per-fault conformal STLs.
  2. For every pair (A, B), find triangles sharing an edge at
     coordinate-precision `snap_m`.
  3. Compute dihedral angle between A's plane and B's plane at the
     shared edge.
  4. If dihedral < `dihedral_deg_max`, mark the edge as a wedge edge.
  5. For each wedge edge, split BOTH faults' incident triangle in
     place, inserting the midpoint as a shared new vertex with
     bit-identical coords on both sides.

Output: modified STLs in --out-stl-dir, plus a `wedge_split_report.json`
listing which edges were split, dihedral angles, and vertex coords.

CRITICAL invariant: a wedge edge identified between A and B must be
split with the SAME midpoint coords in both A and B's STLs.  If the
midpoint is computed independently per fault, the coords will differ
by float64 round-off and the post-split polyline endpoints disagree
across faults — HXT then rejects with `Two segments intersect` or
similar.  This script computes the midpoint ONCE per wedge edge and
applies it to both faults' triangles.

Usage:
  python break_fault_wedges.py \\
      --in-stl-dir output/.../stl_conformal \\
      --out-stl-dir output/.../stl_wedge_broken \\
      --include-fault A --include-fault B [--include-fault C ...] \\
      [--dihedral-deg-max 30.0] [--snap-m 0.1]
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np


# ---------------------------------------------------------------------------
# STL I/O — matching `mesh/inject_corner_touches.py` byte-for-byte.
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
    """Snap-grid dedup at tolerance `tol`."""
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
# Geometry: triangle normal, dihedral angle.
# ---------------------------------------------------------------------------
def _tri_normal(V: np.ndarray, t: np.ndarray) -> np.ndarray:
    a, b, c = V[t[0]], V[t[1]], V[t[2]]
    n = np.cross(b - a, c - a)
    L = np.linalg.norm(n)
    return n / L if L > 0 else n


def _dihedral_deg(n_a: np.ndarray, n_b: np.ndarray) -> float:
    """Angle between the two PLANES of fault A's and fault B's
    triangles meeting at the shared edge.

    Uses |n_A · n_B| (absolute value) so that orientation differences
    between the two faults' STLs don't bias the result.  Dihedral
    range: 0° (coplanar — worst case for wedge) to 90°
    (perpendicular — no wedge concern).
    """
    cos_t = float(abs(np.dot(n_a, n_b)))
    cos_t = max(-1.0, min(1.0, cos_t))
    return float(np.degrees(np.arccos(cos_t)))


# ---------------------------------------------------------------------------
# Edge-key helpers.  An edge is identified by the SNAP-ROUNDED 3-D coords
# of its two endpoints (sorted), so it can be matched across faults that
# share the same polyline endpoint at coord-precision `snap_m`.
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
# Per-fault data + edge index.
# ---------------------------------------------------------------------------
class FaultMesh:
    """Per-fault triangulation with edge → triangle indices index."""
    def __init__(self, name: str, V: np.ndarray, T: np.ndarray,
                 solid_name: str, snap_m: float):
        self.name = name
        self.V = V              # (N, 3)
        self.T = T              # (M, 3)  -> indices into V
        self.solid_name = solid_name
        self.snap_m = snap_m
        # edge-key (snap-rounded) -> list of (tri_idx, edge_local_idx)
        # where edge_local_idx is 0 (v0,v1), 1 (v1,v2), or 2 (v2,v0).
        self.edge_to_tris: dict[
            tuple[tuple[int, int, int], tuple[int, int, int]],
            list[tuple[int, int]]] = {}
        self._build_edge_index()

    def _build_edge_index(self) -> None:
        for ti, t in enumerate(self.T):
            v0 = self.V[t[0]]; v1 = self.V[t[1]]; v2 = self.V[t[2]]
            for li, (p, q) in enumerate(((v0, v1), (v1, v2), (v2, v0))):
                ek = _edge_snap_key(p, q, self.snap_m)
                self.edge_to_tris.setdefault(ek, []).append((ti, li))


# ---------------------------------------------------------------------------
# Polyline-vertex canonicalization.  Ensures fault A's and fault B's
# vertices that share a snap_key are EXACTLY bit-identical post-snap,
# eliminating the sub-cm drift that breaks HXT's PLC recovery on
# CGAL 6.1 autorefine output.
# ---------------------------------------------------------------------------
def collect_polyline_vertex_keys(meshes: dict[str, FaultMesh]
                                  ) -> set[tuple[int, int, int]]:
    """Return snap_keys that appear in at least 2 faults — i.e., every
    cross-fault polyline vertex (curve interior + curve endpoints).
    """
    key_to_faults: dict[tuple[int, int, int], set[str]] = {}
    for name, mesh in meshes.items():
        for v in mesh.V:
            key_to_faults.setdefault(
                _snap_key(v, mesh.snap_m), set()).add(name)
    return {k for k, fs in key_to_faults.items() if len(fs) >= 2}


def snap_polyline_vertices_to_canonical(
        meshes: dict[str, FaultMesh]
        ) -> dict[str, int]:
    """Replace every polyline-vertex coord (vertex whose snap_key is
    shared with another fault) with its canonical position
    `snap_key * snap_m`.  Modifies each FaultMesh's `V` array in place
    and returns a per-fault count of vertices snapped.

    This is the prerequisite for HXT to accept the post-split polyline
    on autorefine-produced inputs: it eliminates the sub-cm drift
    between fault A's and fault B's "shared" vertices, so after the
    split both faults' polylines are bit-identical at every vertex
    (existing endpoints AND inserted midpoints).

    Maximum perturbation per vertex coord: `snap_m / 2`.  For the SAFS
    pipeline default `snap_m = 0.1 m`, this is < 5 cm — well within
    HXT's gmsh STL Merge tolerance (~ 2 mm of 200 km box ≈ 0.4 m).
    """
    polyline_keys = collect_polyline_vertex_keys(meshes)
    counts: dict[str, int] = {}
    for name, mesh in meshes.items():
        snap = mesh.snap_m
        n_snapped = 0
        for vi in range(mesh.V.shape[0]):
            key = _snap_key(mesh.V[vi], snap)
            if key in polyline_keys:
                canon = np.array(key, dtype=np.float64) * snap
                if not np.array_equal(mesh.V[vi], canon):
                    mesh.V[vi] = canon
                    n_snapped += 1
        counts[name] = n_snapped
        # Rebuild the edge index — V coords changed.
        mesh.edge_to_tris.clear()
        mesh._build_edge_index()
    return counts


# ---------------------------------------------------------------------------
# Wedge detection across pairs of faults.
# ---------------------------------------------------------------------------
def find_wedge_edges(meshes: dict[str, FaultMesh],
                     dihedral_deg_max: float,
                     midpoint_mode: str = "canonical"
                     ) -> list[dict]:
    """For each pair (A, B), find every edge whose (snap-rounded) key
    is incident to triangles in BOTH A and B AND whose dihedral angle
    between A's and B's incident planes is below `dihedral_deg_max`.

    Returns a list of records: {fault_a, fault_b, edge_key,
    midpoint, endpoint_a, endpoint_b, endpoint_a_canonical,
    endpoint_b_canonical, dihedral_deg, tri_indices_a, tri_indices_b,
    midpoint_mode}.

    midpoint_mode controls how the inserted midpoint coord is computed:

    - "canonical" (default, post-2026-05-02 fix):
        midpoint = (canon_p1 + canon_p2) / 2, where
        canon_p1 = ek[0] * snap_m and canon_p2 = ek[1] * snap_m.
        Both faults' raw endpoints differ from canonical by < snap_m/2,
        so each fault's resulting post-split polyline (raw_endpoint →
        mid_canon → raw_endpoint) has a small kink at mid_canon — but
        BOTH faults insert the SAME canonical midpoint vertex, which
        gmsh's STL Merge welds into a single shared polyline vertex.
        This eliminates the "two near-coincident lines crossing"
        failure mode that HXT's PLC recovery rejects when many wedge
        edges are split with the legacy "fault_a" mode on
        autorefine-produced inputs (where pa1 ≠ pb1 at sub-cm).

    - "fault_a" (legacy, pre-2026-05-02):
        midpoint = (pa1 + pa2) / 2, where pa1, pa2 are fault A's raw
        endpoint coords.  Preserves fault A's polyline collinearity
        but introduces a kink in fault B's polyline whenever pa_i ≠
        pb_i.  Acceptable on cascade output (where pa = pb to bit
        precision); rejected by HXT on autorefine output above
        ~50 simultaneous wedge splits.  See
        REVIEW_intersection_refinement_investigation_fix.md
        Empirical Validation table.
    """
    if midpoint_mode not in ("canonical", "fault_a"):
        raise ValueError(
            f"midpoint_mode must be 'canonical' or 'fault_a'; "
            f"got {midpoint_mode!r}")
    names = sorted(meshes.keys())
    wedge_records: list[dict] = []
    seen_edges: set[tuple[tuple[int, int, int], tuple[int, int, int]]] = set()

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a_name = names[i]; b_name = names[j]
            mA = meshes[a_name]
            mB = meshes[b_name]
            if mA.snap_m != mB.snap_m:
                raise ValueError(
                    f"snap_m mismatch between fault {a_name} "
                    f"({mA.snap_m}) and fault {b_name} ({mB.snap_m}); "
                    f"all FaultMesh inputs must use the same snap "
                    f"precision for cross-fault edge matching")
            snap = mA.snap_m
            for ek, a_tri_list in mA.edge_to_tris.items():
                if ek not in mB.edge_to_tris:
                    continue
                if ek in seen_edges:
                    continue
                # Skip degenerate edges where the two endpoints share
                # a snap_key (zero-length canonical edge).  A "wedge
                # edge" with length 0 is a vertex pair, not an edge;
                # splitting it would produce degenerate triangles.
                if ek[0] == ek[1]:
                    continue
                # Compute representative dihedral: use the FIRST
                # incident triangle from each fault.
                a_ti, _ = a_tri_list[0]
                b_ti, _ = mB.edge_to_tris[ek][0]
                n_a = _tri_normal(mA.V, mA.T[a_ti])
                n_b = _tri_normal(mB.V, mB.T[b_ti])
                d = _dihedral_deg(n_a, n_b)
                if d >= dihedral_deg_max:
                    continue
                # Recover the actual edge endpoints from one of the
                # tris on each side.  These are each fault's raw
                # vertex coords; they may differ between A and B at
                # sub-cm precision when the input came from CGAL 6.1
                # autorefine.
                def edge_endpoints(mesh: FaultMesh,
                                    ti: int, li: int
                                    ) -> tuple[np.ndarray, np.ndarray]:
                    t = mesh.T[ti]
                    if   li == 0: return mesh.V[t[0]], mesh.V[t[1]]
                    elif li == 1: return mesh.V[t[1]], mesh.V[t[2]]
                    else:         return mesh.V[t[2]], mesh.V[t[0]]
                pa1, pa2 = edge_endpoints(mA, a_ti, a_tri_list[0][1])
                pb1, pb2 = edge_endpoints(mB, b_ti, mB.edge_to_tris[ek][0][1])

                # Canonical endpoint coords: position implied by the
                # snap_key alone.  Both faults' raw endpoints map to
                # the SAME canonical position by construction (they
                # share the snap_key — that's how the edge was matched).
                canon_p1 = np.array(ek[0], dtype=np.float64) * snap
                canon_p2 = np.array(ek[1], dtype=np.float64) * snap

                if midpoint_mode == "canonical":
                    mid = (canon_p1 + canon_p2) / 2.0
                else:  # "fault_a"
                    mid = (pa1 + pa2) / 2.0

                wedge_records.append({
                    "fault_a": a_name,
                    "fault_b": b_name,
                    "endpoint_a": pa1.tolist(),
                    "endpoint_b": pa2.tolist(),
                    "endpoint_a_canonical": canon_p1.tolist(),
                    "endpoint_b_canonical": canon_p2.tolist(),
                    "endpoint_b_raw_in_fault_b": [
                        pb1.tolist(), pb2.tolist()],
                    "midpoint": mid.tolist(),
                    "midpoint_mode": midpoint_mode,
                    "dihedral_deg": d,
                    "tri_indices_a": [ti for ti, _ in a_tri_list],
                    "tri_indices_b": [ti for ti, _ in mB.edge_to_tris[ek]],
                    "edge_key": list(ek),
                })
                seen_edges.add(ek)
    return wedge_records


# ---------------------------------------------------------------------------
# All-polyline bisection (post-2026-05-02 follow-up): bisect EVERY
# cross-fault polyline edge, regardless of dihedral angle.  Used in
# combination with mmg3d_post_pass.py mode='optim_relax_fault' to
# break the polyline-locked wedge tets that mmg3d cannot otherwise
# remove (because the polyline edges are RequiredEdges in mmg3d's
# input, so it can't split them itself).
#
# Implementation: same machinery as `find_wedge_edges` but with the
# dihedral threshold raised above the [0, 90]° valid range so every
# shared edge passes the filter.
# ---------------------------------------------------------------------------
def find_all_polyline_edges(meshes: dict[str, FaultMesh],
                             midpoint_mode: str = "canonical"
                             ) -> list[dict]:
    """Return records for EVERY edge whose snap-rounded key is
    incident to triangles in BOTH faults of any pair (A, B) — i.e.,
    every cross-fault polyline edge.  No dihedral-angle filter.

    Output schema is identical to `find_wedge_edges`; the
    `dihedral_deg` field still records the actual measured angle for
    diagnostics, even though all edges pass.

    Use this for the "bisect every polyline edge" mode (CLI flag
    `--bisect-all-polyline`).  Without that flag, only edges with
    `dihedral_deg < dihedral_deg_max` are split — which leaves
    polyline-locked wedge tets in place wherever the dihedral is
    above the threshold.
    """
    # 180° is above the [0, 90] valid range for the |n_A · n_B|
    # dihedral so every edge passes `d < dihedral_deg_max`.
    return find_wedge_edges(meshes,
                             dihedral_deg_max=180.0,
                             midpoint_mode=midpoint_mode)


# ---------------------------------------------------------------------------
# Splitting logic.  For a fault mesh with edge `(a, b)` incident to
# triangle `t = {a, b, c}`, replace `t` with two triangles
# `{a, m, c}` and `{m, b, c}`.  Vertex `m` is added as a new vertex.
# We APPLY the same midpoint coord across every fault that shares
# this edge so cross-fault conformity is preserved bit-exactly.
# ---------------------------------------------------------------------------
def split_wedge_edges_in_fault(mesh: FaultMesh,
                                edge_to_midpoint: dict[
                                    tuple[tuple[int, int, int],
                                          tuple[int, int, int]],
                                    np.ndarray]
                                ) -> tuple[np.ndarray, np.ndarray, int]:
    """Split every triangle in `mesh` whose edge is in `edge_to_midpoint`.

    Returns (new_V, new_T, n_splits).
    """
    V = mesh.V
    T = mesh.T
    new_V: list[np.ndarray] = [v for v in V]
    new_T: list[list[int]] = []
    snap = mesh.snap_m
    # Track which midpoint-key has been added as a vertex; reuse the
    # vertex index if already added (so triangles sharing the midpoint
    # all reference the same vertex).
    midkey_to_vidx: dict[tuple[int, int, int], int] = {}

    n_splits = 0
    for ti, t in enumerate(T):
        v0, v1, v2 = int(t[0]), int(t[1]), int(t[2])
        # Examine all 3 edges of this triangle.  For each, see if it's
        # one of the wedge edges to be split.
        edges = [(0, 1, v0, v1, v2), (1, 2, v1, v2, v0), (2, 0, v2, v0, v1)]
        # We may have multiple edges to split in one triangle; handle
        # the most common case first (1 edge) and skip the rare 2/3
        # edge cases.  Empirically a fault triangle rarely shares >1
        # cross-fault polyline edge.
        edges_to_split = []
        for li, lj, vi, vj, vk in edges:
            ek = _edge_snap_key(V[vi], V[vj], snap)
            if ek in edge_to_midpoint:
                edges_to_split.append((li, lj, vi, vj, vk, ek))

        if not edges_to_split:
            new_T.append([v0, v1, v2])
            continue

        # Allocate / re-use midpoint vertex for each wedge edge.
        def _get_m_idx(ek_local):
            mid_coord = edge_to_midpoint[ek_local]
            midkey = _snap_key(mid_coord, snap)
            if midkey in midkey_to_vidx:
                return midkey_to_vidx[midkey]
            idx = len(new_V)
            new_V.append(mid_coord)
            midkey_to_vidx[midkey] = idx
            return idx

        if len(edges_to_split) == 1:
            li, lj, vi, vj, vk, ek = edges_to_split[0]
            m_idx = _get_m_idx(ek)
            # Triangle (vi, vj, vk) — split edge (vi, vj) at m.
            new_T.append([vi, m_idx, vk])
            new_T.append([m_idx, vj, vk])
            n_splits += 1
        elif len(edges_to_split) == 2:
            # Two wedge edges meeting at a shared corner of the triangle.
            # The two edges share exactly one vertex (the "apex").  Order
            # the canonical labels (a, b, c) so the two split edges are
            # (a, b) and (b, c), meeting at b.  Then sub-triangles:
            #   (a, m_ab, c), (m_ab, b, m_bc), (m_ab, m_bc, c).
            # Find the shared corner.
            edge_pair_endpoints = []
            for li_, lj_, vi_, vj_, vk_, ek_ in edges_to_split:
                edge_pair_endpoints.append((vi_, vj_, vk_, ek_))
            # Each edges_to_split entry: (li, lj, vi, vj, vk, ek)
            # vi-vj is the wedge edge; vk is the third vertex.
            (vi1, vj1, vk1, ek1), (vi2, vj2, vk2, ek2) = edge_pair_endpoints
            # The two wedge edges (vi1, vj1) and (vi2, vj2) share one
            # vertex (the apex `b`).  The other endpoints are `a` and `c`.
            shared = set([vi1, vj1]) & set([vi2, vj2])
            if len(shared) != 1:
                sys.stderr.write(
                    f"[break_fault_wedges] WARN: fault {mesh.name} tri "
                    f"{ti} 2-wedge-edge case has non-corner sharing; "
                    f"keeping triangle un-split.\n")
                new_T.append([v0, v1, v2])
                continue
            b = shared.pop()
            a = vi1 if vi1 != b else vj1
            c = vi2 if vi2 != b else vj2
            # Verify: a, b, c are the original triangle vertices.
            if {a, b, c} != {v0, v1, v2}:
                sys.stderr.write(
                    f"[break_fault_wedges] WARN: fault {mesh.name} tri "
                    f"{ti} 2-wedge-edge case label mismatch; keeping "
                    f"triangle un-split.\n")
                new_T.append([v0, v1, v2])
                continue
            m_ab = _get_m_idx(ek1)
            m_bc = _get_m_idx(ek2)
            # Three sub-triangles preserving (a, b, c) orientation.
            # Original triangle was emitted with vertex order (v0, v1, v2)
            # in the .geo's CCW.  Maintain that ordering by checking
            # which input order corresponds to (a, b, c).
            new_T.append([a, m_ab, c])
            new_T.append([m_ab, b, m_bc])
            new_T.append([m_ab, m_bc, c])
            n_splits += 1
        else:
            # 3 wedge edges in one triangle is exceptionally rare; warn
            # and keep un-split.  In a future revision we'd split into
            # 4 sub-tris using the 3 edge midpoints.
            sys.stderr.write(
                f"[break_fault_wedges] WARN: fault {mesh.name} tri {ti} "
                f"has {len(edges_to_split)} wedge edges; un-split (3-"
                f"edge case not implemented).\n")
            new_T.append([v0, v1, v2])

    new_V_arr = np.asarray(new_V, dtype=np.float64)
    new_T_arr = np.asarray(new_T, dtype=np.int64)
    return new_V_arr, new_T_arr, n_splits


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-stl-dir", required=True, type=Path)
    p.add_argument("--out-stl-dir", required=True, type=Path)
    p.add_argument("--include-fault", action="append", default=[],
                   metavar="SHORT_NAME")
    p.add_argument("--snap-m", type=float, default=0.1,
                   help="Coordinate snap precision for matching shared "
                        "polyline edges across faults.  Default 0.1 m.")
    p.add_argument("--dihedral-deg-max", type=float, default=30.0,
                   help="Edges where the dihedral angle between fault "
                        "A's and fault B's incident triangle planes is "
                        "below this threshold are considered 'wedge "
                        "edges' and split.  Default 30 deg.")
    p.add_argument("--midpoint-mode",
                   choices=("canonical", "fault_a"), default="canonical",
                   help="Where to place the inserted midpoint vertex. "
                        "'canonical' (default, post-2026-05-02 fix) "
                        "uses (snap_key * snap_m) for both endpoints "
                        "so both faults agree on the midpoint position "
                        "to bit-precision.  'fault_a' (legacy) uses "
                        "fault A's raw endpoint coords; this preserves "
                        "fault A's polyline collinearity but moves the "
                        "midpoint off fault B's polyline whenever "
                        "fault B's vertex coords differ at sub-cm.")
    p.add_argument("--snap-polyline-vertices",
                   dest="snap_polyline_vertices",
                   action="store_true", default=True,
                   help="Before detecting wedge edges, replace every "
                        "polyline vertex coord (vertex whose snap_key "
                        "appears in ≥2 faults) with its canonical "
                        "position snap_key * snap_m.  This eliminates "
                        "the sub-cm drift between fault A's and fault "
                        "B's shared vertices that breaks HXT's PLC "
                        "recovery on autorefine-produced inputs.  "
                        "Default: ON.")
    p.add_argument("--no-snap-polyline-vertices",
                   dest="snap_polyline_vertices",
                   action="store_false",
                   help="Skip polyline-vertex canonicalization "
                        "(legacy behaviour, pre-2026-05-02).")
    p.add_argument("--bisect-all-polyline",
                   action="store_true",
                   help="Bisect EVERY cross-fault polyline edge "
                        "(1->2 split with cross-fault canonical "
                        "midpoint), regardless of dihedral angle.  "
                        "Used in combination with mmg3d_post_pass.py "
                        "mode='optim_relax_fault' to break the "
                        "polyline-locked wedge tets that mmg3d "
                        "cannot otherwise remove.  Default OFF "
                        "(only edges with dihedral < "
                        "--dihedral-deg-max are split).")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    if not args.in_stl_dir.is_dir():
        raise SystemExit(f"--in-stl-dir not found: {args.in_stl_dir}")
    if not args.include_fault:
        raise SystemExit("--include-fault must appear at least once")
    args.out_stl_dir.mkdir(parents=True, exist_ok=True)

    # Load and dedup each fault's triangulation.
    meshes: dict[str, FaultMesh] = {}
    for short in args.include_fault:
        path = args.in_stl_dir / f"{short}.stl"
        if not path.exists():
            raise SystemExit(f"missing input STL: {path}")
        V, T, solid_name = _read_ascii_stl(path)
        V, T = _dedup_vertices(V, T, tol=args.snap_m)
        meshes[short] = FaultMesh(short, V, T, solid_name, args.snap_m)
        print(f"  load {short}: V={V.shape[0]}, T={T.shape[0]}",
              file=sys.stderr)

    # Pre-step: snap polyline vertices to canonical (pre-fix this was
    # never done; see REVIEW_intersection_refinement_investigation_fix.md
    # Empirical Validation for the rationale).
    snap_counts: dict[str, int] = {}
    if args.snap_polyline_vertices:
        snap_counts = snap_polyline_vertices_to_canonical(meshes)
        n_total = sum(snap_counts.values())
        print(f"  snap polyline vertices to canonical: "
              f"{n_total} vertex coord(s) moved across "
              f"{len(meshes)} faults", file=sys.stderr)
        for short, n in snap_counts.items():
            print(f"    {short}: {n} polyline-vertex(es) snapped",
                  file=sys.stderr)

    # Find polyline edges to split.  Two modes:
    #   default:           edges with dihedral < --dihedral-deg-max
    #                      (the "wedge" subset).
    #   --bisect-all-polyline: every cross-fault polyline edge, no
    #                          dihedral filter.
    if args.bisect_all_polyline:
        wedges = find_all_polyline_edges(meshes,
                                          midpoint_mode=args.midpoint_mode)
        print(f"  ALL polyline edges (no dihedral filter, "
              f"midpoint_mode={args.midpoint_mode}): {len(wedges)}",
              file=sys.stderr)
    else:
        wedges = find_wedge_edges(meshes, args.dihedral_deg_max,
                                   midpoint_mode=args.midpoint_mode)
        print(f"  wedge edges (dihedral < {args.dihedral_deg_max}°, "
              f"midpoint_mode={args.midpoint_mode}): {len(wedges)}",
              file=sys.stderr)

    # Build edge -> midpoint map (a wedge edge belongs to exactly one
    # pair of faults; the midpoint is shared between them).
    edge_to_midpoint: dict[tuple[tuple[int, int, int],
                                  tuple[int, int, int]],
                            np.ndarray] = {}
    for w in wedges:
        ek = (tuple(w["edge_key"][0]), tuple(w["edge_key"][1]))
        edge_to_midpoint[ek] = np.asarray(w["midpoint"], dtype=np.float64)

    # Apply splits to every fault that has any wedge edge.
    n_splits_total = 0
    for short, mesh in meshes.items():
        new_V, new_T, n_splits = split_wedge_edges_in_fault(
            mesh, edge_to_midpoint)
        n_splits_total += n_splits
        print(f"  split {short}: {n_splits} edge-split(s); "
              f"V={mesh.V.shape[0]}->{new_V.shape[0]}, "
              f"T={mesh.T.shape[0]}->{new_T.shape[0]}",
              file=sys.stderr)
        out_path = args.out_stl_dir / f"{short}.stl"
        _write_ascii_stl(out_path, new_V, new_T,
                          mesh.solid_name or f"SAFS:{short}:wedge_broken")

    report_path = args.report_json or (
        args.out_stl_dir / "wedge_split_report.json")
    with report_path.open("w") as fh:
        json.dump({
            "snap_m": args.snap_m,
            "dihedral_deg_max": args.dihedral_deg_max,
            "midpoint_mode": args.midpoint_mode,
            "snap_polyline_vertices": bool(args.snap_polyline_vertices),
            "bisect_all_polyline": bool(args.bisect_all_polyline),
            "polyline_vertices_snapped_per_fault": snap_counts,
            "n_wedge_edges": len(wedges),
            "n_splits_applied": n_splits_total,
            "wedges": wedges,
        }, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
