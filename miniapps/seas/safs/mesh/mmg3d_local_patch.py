#!/usr/bin/env python
"""Local sliver-targeted mmg3d patch pass.

Implements R-001 from `REVIEW_local_sliver_targeting.md`: the slivers
in the Step-6 production mesh are 0.027 % of the mesh and 58.5 % are
in 5 spatial buckets at the cross-fault polyline triple-junction
(`coav_missioncreek × sbmt_garnethill`).  A global mmg3d pass spends
99.97 % of its work on tets that don't need help.

This script:

  1. Identifies tets with γ < `gamma_thresh` (default 0.01).
  2. Builds a "subdomain" containing every tet whose centroid is
     within `halo_radius_m` of any sliver tet's centroid.
  3. Identifies halo-boundary triangles (shared between an
     in-subdomain tet and an outside-subdomain tet) and tags them
     with a special tag (default 200) marked `RequiredTriangles` in
     mmg3d's input.  This LOCKS the halo boundary so mmg3d's output
     can be merged back without geometric drift.
  4. Runs mmg3d on the subdomain medit file with tight tolerances
     (`-hmin 50`, `-hgradreq 1.1`, `-hausd 5`).  Polyline / fault
     protection is preserved (re-uses `optim_relax_fault` mode
     plumbing from `mmg3d_post_pass.py`).
  5. Stitches the mmg3d-modified subdomain back into the input
     mesh: outside-subdomain tets / vertices keep their original
     coords; in-subdomain tets are replaced by the mmg3d output;
     halo-boundary vertices match by snap-key precision.

Pipeline position (post `mmg3d_post_pass.py`):
    ... → mmg3d_post_pass (global) → THIS → validate_msh

Usage:
    python mmg3d_local_patch.py \\
        --in-msh safs_mmg.msh --out-msh safs_local_patch.msh \\
        --provenance-json fault_provenance.json \\
        --transform-json transform.json \\
        [--gamma-thresh 0.01] [--halo-radius-m 2000] \\
        [--hmin 50] [--hgrad 1.1] [--hausd 5]

Requires `mmg3d_O3` on PATH.
"""
from __future__ import annotations
import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import mmg3d_post_pass as m3p  # noqa: E402


# ---------------------------------------------------------------------------
# Mesh reading helpers.
# ---------------------------------------------------------------------------
def _read_msh(msh_path: Path) -> dict:
    """Read a gmsh .msh file via meshio and return a flat dict."""
    import meshio
    m = meshio.read(msh_path)
    points = np.asarray(m.points, dtype=np.float64)
    tris: list[tuple[int, int, int]] = []
    tri_tags: list[int] = []
    tets: list[tuple[int, int, int, int]] = []
    tet_tags: list[int] = []
    for cb, tags in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        if cb.type == "triangle":
            for tri, tag in zip(cb.data, tags):
                tris.append((int(tri[0]), int(tri[1]), int(tri[2])))
                tri_tags.append(int(tag))
        elif cb.type == "tetra":
            for tet, tag in zip(cb.data, tags):
                tets.append((int(tet[0]), int(tet[1]),
                              int(tet[2]), int(tet[3])))
                tet_tags.append(int(tag))
    return {
        "points":   points,
        "tris":     np.asarray(tris, dtype=np.int64),
        "tri_tags": np.asarray(tri_tags, dtype=np.int32),
        "tets":     np.asarray(tets, dtype=np.int64),
        "tet_tags": np.asarray(tet_tags, dtype=np.int32),
    }


def _min_edge_of_tets(points: np.ndarray, tets: np.ndarray) -> float:
    """Global minimum tet edge length (m).  Returns 0.0 for empty input."""
    if tets.shape[0] == 0:
        return 0.0
    p0 = points[tets[:, 0]]; p1 = points[tets[:, 1]]
    p2 = points[tets[:, 2]]; p3 = points[tets[:, 3]]
    e_sq = np.stack([
        np.sum((p1 - p0) ** 2, axis=1),
        np.sum((p2 - p0) ** 2, axis=1),
        np.sum((p3 - p0) ** 2, axis=1),
        np.sum((p2 - p1) ** 2, axis=1),
        np.sum((p3 - p1) ** 2, axis=1),
        np.sum((p3 - p2) ** 2, axis=1),
    ], axis=1)
    return float(np.sqrt(e_sq.min()))


def _gamma_per_tet(points: np.ndarray, tets: np.ndarray
                    ) -> np.ndarray:
    """Vectorized γ per tet (same metric as `_gamma_min_of_msh`)."""
    if tets.shape[0] == 0:
        return np.zeros(0, dtype=np.float64)
    p0 = points[tets[:, 0]]
    p1 = points[tets[:, 1]]
    p2 = points[tets[:, 2]]
    p3 = points[tets[:, 3]]
    e2 = (np.sum((p1 - p0) ** 2, axis=1)
           + np.sum((p2 - p0) ** 2, axis=1)
           + np.sum((p3 - p0) ** 2, axis=1)
           + np.sum((p2 - p1) ** 2, axis=1)
           + np.sum((p3 - p1) ** 2, axis=1)
           + np.sum((p3 - p2) ** 2, axis=1))
    vol = np.abs(
        np.einsum('ij,ij->i', p1 - p0,
                   np.cross(p2 - p0, p3 - p0))) / 6.0
    g = np.zeros(tets.shape[0], dtype=np.float64)
    mask = (vol > 0) & (e2 > 0)
    g[mask] = (12.0 * (3.0 * vol[mask]) ** (2.0 / 3.0)) / e2[mask]
    return g


# ---------------------------------------------------------------------------
# Subdomain extraction.
# ---------------------------------------------------------------------------
def _identify_subdomain_tets(
        points: np.ndarray,
        tets: np.ndarray,
        bad_tet_idx: np.ndarray,
        halo_radius_m: float
        ) -> np.ndarray:
    """Return a boolean mask over `tets`: True iff the tet's centroid
    is within `halo_radius_m` of any bad tet's centroid.
    """
    if bad_tet_idx.size == 0:
        return np.zeros(tets.shape[0], dtype=bool)
    centroids = points[tets].mean(axis=1)        # (N_t, 3)
    bad_centroids = centroids[bad_tet_idx]       # (N_bad, 3)
    in_subdomain = np.zeros(tets.shape[0], dtype=bool)
    r2 = halo_radius_m * halo_radius_m
    # Per bad-centroid expansion.  N_bad is small (≤ few hundred);
    # N_t is large (~1M).  Vectorize over tets per bad centroid.
    for bc in bad_centroids:
        d2 = np.sum((centroids - bc) ** 2, axis=1)
        in_subdomain |= (d2 < r2)
    return in_subdomain


def _build_tet_face_adjacency(tets: np.ndarray
                                ) -> dict[tuple[int, int, int], list[int]]:
    """Map each (sorted-vertex) face triple to the list of tet indices
    that contain it.  An interior face has 2 tets; a boundary face has 1.
    """
    face_to_tets: dict[tuple[int, int, int], list[int]] = {}
    for ti in range(tets.shape[0]):
        a, b, c, d = (int(tets[ti, 0]), int(tets[ti, 1]),
                       int(tets[ti, 2]), int(tets[ti, 3]))
        for face in (
            tuple(sorted((b, c, d))),
            tuple(sorted((a, c, d))),
            tuple(sorted((a, b, d))),
            tuple(sorted((a, b, c))),
        ):
            face_to_tets.setdefault(face, []).append(ti)
    return face_to_tets


def _identify_halo_boundary_faces(
        tets: np.ndarray,
        in_subdomain: np.ndarray,
        face_to_tets: dict[tuple[int, int, int], list[int]]
        ) -> list[tuple[int, int, int]]:
    """Return the list of (sorted-vertex) faces that are shared
    between exactly one in-subdomain tet and one outside-subdomain
    tet.  These are the halo boundary triangles to lock.
    """
    halo: list[tuple[int, int, int]] = []
    for face, ti_list in face_to_tets.items():
        if len(ti_list) != 2:
            continue
        t0, t1 = ti_list
        in0 = in_subdomain[t0]
        in1 = in_subdomain[t1]
        if in0 != in1:
            halo.append(face)
    return halo


# ---------------------------------------------------------------------------
# Subdomain medit writer.
# ---------------------------------------------------------------------------
def _write_subdomain_medit(
        out_path: Path,
        points: np.ndarray,
        in_subdomain: np.ndarray,
        tets: np.ndarray,
        tet_tags: np.ndarray,
        tris: np.ndarray,
        tri_tags: np.ndarray,
        halo_faces: list[tuple[int, int, int]],
        halo_tag: int,
        polyline_keys: set[tuple[int, int, int]] | None,
        snap_m: float,
        free_surface_clearance_m: float,
        z_top: float,
        fault_tri_tag: int,
        polyline_relax: bool = False,
        ) -> dict:
    """Write the subdomain to medit format.

    Triangles emitted:
      - All halo-boundary triangles, tagged `halo_tag`,
        marked RequiredTriangles (locks the halo).
      - All input fault triangles whose 3 vertices are all in
        subdomain (tag 100).  Marked RequiredTriangles ONLY in
        the perimeter / free-surface sense via RequiredEdges,
        following the optim_relax_fault convention.
      - All input box-face triangles whose 3 vertices are all in
        subdomain (tags 1-6).
    Edges emitted as RequiredEdges:
      - Polyline edges (both endpoints have polyline snap_key)
        — UNLESS `polyline_relax=True`, in which case polyline
        edges are NOT marked Required (mmg3d may collapse /
        subdivide them; the post-stitch canonical-snap restores
        cross-fault conformity).
      - Free-surface trace edges (both endpoints near z=z_top).
      - Fault-perimeter edges (incident to exactly one fault tri).
    Vertices emitted as RequiredVertices:
      - Polyline vertices (always — even with polyline_relax=True,
        polyline VERTICES stay pinned at canonical position so
        cross-fault conformity is preserved at the vertex level
        even when polyline edges can collapse).

    Returns a dict mapping local→global vertex index plus counters.
    """
    # Map subdomain tet → local 0..N_sub_t-1.
    sub_tet_idx = np.flatnonzero(in_subdomain)
    sub_tets = tets[sub_tet_idx]  # (N_sub_t, 4) in global vertex IDs
    sub_tet_tags = tet_tags[sub_tet_idx]

    # Subdomain vertex set: union of all sub_tets vertices.
    sub_v_global = np.unique(sub_tets.flatten())
    # Map global → local.
    global_to_local = -np.ones(points.shape[0], dtype=np.int64)
    global_to_local[sub_v_global] = np.arange(sub_v_global.shape[0])

    # Re-index sub_tets to local IDs.
    sub_tets_local = global_to_local[sub_tets]
    # Sanity check.
    if (sub_tets_local < 0).any():
        raise RuntimeError(
            "subdomain tets reference vertices outside subdomain "
            "vertex set — internal error")

    sub_points = points[sub_v_global]  # (N_sub_v, 3)

    # Filter triangles: keep tris whose all 3 vertices are in subdomain
    # AND whose vertex-set matches at least one subdomain tet's face.
    #
    # The vertex-only filter (the previous version) was insufficient:
    # a triangle whose 3 vertices all lie on the subdomain boundary
    # (e.g., on the halo) can have its 2 incident tets BOTH outside
    # the subdomain.  In that case the triangle's vertex-set does not
    # match any subdomain tet's face, so mmg3d rejects it via
    # `MMG5_chkBdryTria_deleteExtraTriangles`.  After mmg3d returns,
    # the dropped triangle is gone, the corresponding tet-boundary
    # face has no tagged entry → SURFACE HOLE.
    #
    # Build the set of subdomain tet faces (frozensets of vertex IDs),
    # and only emit a triangle if its 3 vertices match such a face.
    sub_tet_faces: set[frozenset[int]] = set()
    for ti_local in range(sub_tets.shape[0]):
        a, b, c, d = (int(sub_tets[ti_local, 0]),
                       int(sub_tets[ti_local, 1]),
                       int(sub_tets[ti_local, 2]),
                       int(sub_tets[ti_local, 3]))
        sub_tet_faces.add(frozenset((b, c, d)))
        sub_tet_faces.add(frozenset((a, c, d)))
        sub_tet_faces.add(frozenset((a, b, d)))
        sub_tet_faces.add(frozenset((a, b, c)))
    sub_tris_global: list[tuple[int, int, int]] = []
    sub_tri_tags_out: list[int] = []
    n_skipped_extra_bdry_tri = 0
    for ti in range(tris.shape[0]):
        a, b, c = int(tris[ti, 0]), int(tris[ti, 1]), int(tris[ti, 2])
        if not (global_to_local[a] >= 0
                 and global_to_local[b] >= 0
                 and global_to_local[c] >= 0):
            continue
        if frozenset((a, b, c)) not in sub_tet_faces:
            # Triangle has all vertices in subdomain but no subdomain
            # tet uses it as a face.  Skip — mmg3d would drop it as
            # an "extra boundary" anyway, and including it in the
            # medit causes silent triangle loss in the output.
            n_skipped_extra_bdry_tri += 1
            continue
        sub_tris_global.append((a, b, c))
        sub_tri_tags_out.append(int(tri_tags[ti]))

    # Add halo-boundary triangles with halo_tag.  Halo faces have
    # vertices that are all in subdomain (because both incident tets
    # share the face, and at least one of those tets is in the
    # subdomain so its vertices are in the subdomain set).
    #
    # BUG FIX (2026-05-02 surface-holes regression): if a halo face
    # COINCIDES with an existing tagged surface triangle (box face
    # tag 1-6, fault tag 100), do NOT add a halo_tag duplicate —
    # instead, mark that existing triangle as a RequiredTriangle so
    # mmg3d preserves it under its ORIGINAL tag.  Adding a duplicate
    # caused the medit input to contain the same face twice with
    # different tags; mmg3d kept the Required (halo_tag) version
    # and could drop/modify the non-Required (original tag) version.
    # After stitch-back drops halo_tag triangles, the original
    # tagged triangle was lost → surface hole.
    sub_tri_face_to_local_idx: dict[
        frozenset[int], int] = {}
    for li, (a, b, c) in enumerate(sub_tris_global):
        sub_tri_face_to_local_idx[frozenset((a, b, c))] = li

    halo_tris_global: list[tuple[int, int, int]] = []
    halo_collisions_with_existing_tri: list[int] = []
    for face in halo_faces:
        a, b, c = face
        if (global_to_local[a] < 0
                or global_to_local[b] < 0
                or global_to_local[c] < 0):
            # Should not happen if subdomain is well-formed.
            continue
        existing_li = sub_tri_face_to_local_idx.get(
            frozenset((a, b, c)))
        if existing_li is not None:
            # The halo face IS an existing tagged surface triangle
            # (box / fault).  Mark it Required under its original
            # tag rather than adding a halo_tag duplicate.
            halo_collisions_with_existing_tri.append(existing_li)
        else:
            halo_tris_global.append((a, b, c))

    # Combine with halo last so we know their indices.
    all_tris_global = sub_tris_global + halo_tris_global
    all_tri_tags = sub_tri_tags_out + [halo_tag] * len(halo_tris_global)
    halo_tri_local_indices = list(range(
        len(sub_tris_global),
        len(sub_tris_global) + len(halo_tris_global)))

    # Required triangle set:
    #   1. ALL halo triangles (locks the halo).
    #   2. Halo faces that coincided with existing tagged surface
    #      triangles — Required under their original tag (prior fix).
    #   3. ALL box-face triangles (tags 1-6) in the subdomain —
    #      mmg3d with polyline-relax (no `-nosurf`) and `-hmin 50`
    #      can otherwise collapse short edges in box-face triangles,
    #      orphaning the original tag-5 entry in `tris[]` and
    #      producing surface holes at z=0 (empirical regression from
    #      Step 6 production runs).
    required_tri_local: list[int] = list(halo_tri_local_indices)
    required_tri_local.extend(halo_collisions_with_existing_tri)
    # Add every box-face (tag != fault_tri_tag) sub-triangle.
    box_face_local_v: set[int] = set()
    for li, tag in enumerate(sub_tri_tags_out):
        if int(tag) != int(fault_tri_tag):
            # Avoid duplicates: only add if not already in
            # halo_collisions_with_existing_tri (set membership).
            if li not in halo_collisions_with_existing_tri:
                required_tri_local.append(li)
            # Also collect every vertex of every box-face triangle
            # — these become RequiredVertices below so mmg3d cannot
            # insert NEW Steiner vertices that land on box-face
            # planes (which would create unlabeled boundary faces
            # in the output).  Required for the box-face hole fix
            # in polyline-relax mode.
            (a, b, c) = sub_tris_global[li]
            box_face_local_v.add(int(global_to_local[a]))
            box_face_local_v.add(int(global_to_local[b]))
            box_face_local_v.add(int(global_to_local[c]))

    # Polyline / free-surface / fault-perimeter protections, in
    # local vertex indices.
    required_vertex_local: set[int] = set()
    required_edge_pairs: set[tuple[int, int]] = set()

    if polyline_keys:
        def _snap_key(v: np.ndarray) -> tuple[int, int, int]:
            return (int(round(v[0] / snap_m)),
                    int(round(v[1] / snap_m)),
                    int(round(v[2] / snap_m)))
        # Local-vertex snap_keys.
        is_polyline_v = np.array(
            [_snap_key(sub_points[li]) in polyline_keys
             for li in range(sub_points.shape[0])], dtype=bool)
        # Polyline VERTICES are ALWAYS marked Required.  Even with
        # polyline_relax=True, vertex positions stay pinned at
        # canonical positions — only polyline EDGES become
        # collapsible / splittable.  This preserves cross-fault
        # conformity at the vertex level: shared-snap-key vertices
        # remain bit-identical between faults across the patch.
        for li in np.flatnonzero(is_polyline_v):
            required_vertex_local.add(int(li))
        # Polyline edges in subdomain fault triangles.
        # If polyline_relax=True, SKIP this loop — polyline edges
        # are left non-Required so mmg3d can collapse / split
        # them (this is what kills the 0.20 m needle edges).
        if not polyline_relax:
            for (a, b, c), tag in zip(sub_tris_global,
                                       sub_tri_tags_out):
                if tag != fault_tri_tag:
                    continue
                la = int(global_to_local[a])
                lb = int(global_to_local[b])
                lc = int(global_to_local[c])
                for u, v in ((la, lb), (lb, lc), (lc, la)):
                    if is_polyline_v[u] and is_polyline_v[v]:
                        required_edge_pairs.add(
                            (u, v) if u < v else (v, u))
    else:
        is_polyline_v = np.zeros(sub_points.shape[0], dtype=bool)

    # Free-surface vertices (fault-touching only, per R-005).
    fault_vertex_set: set[int] = set()
    for (a, b, c), tag in zip(sub_tris_global, sub_tri_tags_out):
        if tag == fault_tri_tag:
            fault_vertex_set.update(
                (int(global_to_local[a]),
                 int(global_to_local[b]),
                 int(global_to_local[c])))
    is_freesurface_v = np.array(
        [(li in fault_vertex_set
          and sub_points[li, 2] >= z_top - free_surface_clearance_m)
         for li in range(sub_points.shape[0])], dtype=bool)
    # Free-surface trace edges.
    for (a, b, c), tag in zip(sub_tris_global, sub_tri_tags_out):
        if tag != fault_tri_tag:
            continue
        la = int(global_to_local[a])
        lb = int(global_to_local[b])
        lc = int(global_to_local[c])
        for u, v in ((la, lb), (lb, lc), (lc, la)):
            if is_freesurface_v[u] and is_freesurface_v[v]:
                required_edge_pairs.add(
                    (u, v) if u < v else (v, u))

    # Fault perimeter (within subdomain): edges incident to exactly
    # one subdomain fault triangle.  Caveat: this is a SUBDOMAIN-
    # local count.  A polyline edge that's incident to fault A and
    # fault B in the full mesh may appear as count=1 here if the
    # subdomain extraction kept fault A's triangle but not fault B's
    # triangle on the other side of the polyline.  Such edges are
    # NOT true fault-patch perimeter — they're just where the
    # subdomain boundary cuts through the fault triangulation.
    # When polyline_relax=True we explicitly want polyline edges to
    # be modifiable; we therefore EXCLUDE polyline edges from the
    # perimeter set in relax mode.
    fault_edge_count_local: dict[tuple[int, int], int] = {}
    for (a, b, c), tag in zip(sub_tris_global, sub_tri_tags_out):
        if tag != fault_tri_tag:
            continue
        la = int(global_to_local[a])
        lb = int(global_to_local[b])
        lc = int(global_to_local[c])
        for u, v in ((la, lb), (lb, lc), (lc, la)):
            ek = (u, v) if u < v else (v, u)
            fault_edge_count_local[ek] = (
                fault_edge_count_local.get(ek, 0) + 1)
    for ek, n in fault_edge_count_local.items():
        if n == 1:
            u, v = ek
            # Skip polyline edges in relax mode (see comment above).
            if (polyline_relax
                    and is_polyline_v[u]
                    and is_polyline_v[v]):
                continue
            required_edge_pairs.add(ek)

    # Augment RequiredVertices with all box-face vertices so mmg3d
    # cannot insert a new Steiner point that lands on a box-face
    # plane — preventing unlabeled boundary faces at z = z_top in
    # polyline-relax mode (empirical fix for Step 6 box-top holes).
    required_vertex_local.update(box_face_local_v)

    # Write medit.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_v = sub_points.shape[0]
    n_t_tri = len(all_tris_global)
    n_tetra = sub_tets_local.shape[0]
    required_edges_list = sorted(required_edge_pairs)
    required_vertex_idx_1based = sorted(
        vi + 1 for vi in required_vertex_local)
    with out_path.open("w") as f:
        f.write("MeshVersionFormatted 2\n\nDimension 3\n\n")

        f.write(f"Vertices\n{n_v}\n")
        for v in sub_points:
            f.write(f"{v[0]:+.17e} {v[1]:+.17e} "
                    f"{v[2]:+.17e} 0\n")
        f.write("\n")

        f.write(f"Triangles\n{n_t_tri}\n")
        for (a, b, c), tag in zip(all_tris_global, all_tri_tags):
            la = int(global_to_local[a]) + 1
            lb = int(global_to_local[b]) + 1
            lc = int(global_to_local[c]) + 1
            f.write(f"{la} {lb} {lc} {tag}\n")
        f.write("\n")

        if required_tri_local:
            f.write(f"RequiredTriangles\n{len(required_tri_local)}\n")
            for li in required_tri_local:
                f.write(f"{li + 1}\n")  # 1-based
            f.write("\n")

        if required_edges_list:
            f.write(f"Edges\n{len(required_edges_list)}\n")
            for (u, v) in required_edges_list:
                f.write(f"{u + 1} {v + 1} 0\n")
            f.write("\n")
            f.write(f"RequiredEdges\n{len(required_edges_list)}\n")
            for i in range(len(required_edges_list)):
                f.write(f"{i + 1}\n")
            f.write("\n")

        if required_vertex_idx_1based:
            f.write(f"RequiredVertices\n"
                    f"{len(required_vertex_idx_1based)}\n")
            for vi in required_vertex_idx_1based:
                f.write(f"{vi}\n")
            f.write("\n")

        f.write(f"Tetrahedra\n{n_tetra}\n")
        for (a, b, c, d), tag in zip(sub_tets_local, sub_tet_tags):
            f.write(f"{a + 1} {b + 1} {c + 1} {d + 1} {tag}\n")
        f.write("\n")

        f.write("End\n")

    return {
        "n_subdomain_vertices":  int(n_v),
        "n_subdomain_tets":      int(n_tetra),
        "n_subdomain_triangles": int(n_t_tri),
        "n_halo_boundary_triangles": int(len(halo_tris_global)),
        "n_required_triangles":  int(len(required_tri_local)),
        "n_required_edges":      int(len(required_edges_list)),
        "n_required_vertices":   int(len(required_vertex_local)),
        "polyline_relax":        bool(polyline_relax),
        "n_skipped_extra_bdry_tri": int(n_skipped_extra_bdry_tri),
        "global_to_local":       global_to_local,
        "sub_v_global":          sub_v_global,
        "sub_tet_idx":           sub_tet_idx,
    }


# ---------------------------------------------------------------------------
# Stitch-back.
# ---------------------------------------------------------------------------
def _read_subdomain_medit(mesh_path: Path
                            ) -> tuple[np.ndarray, np.ndarray,
                                        np.ndarray, np.ndarray,
                                        np.ndarray]:
    """Parse medit `.mesh` output: returns (points, tris, tri_tags,
    tets, tet_tags), all 0-based vertex indices.  Reuses
    `_convert_medit_to_msh`'s parser logic but returns raw arrays.
    """
    points: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    tri_tags: list[int] = []
    tets: list[tuple[int, int, int, int]] = []
    tet_tags: list[int] = []
    tokens: list[str] = []
    with mesh_path.open() as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if s:
                tokens.extend(s.split())
    sections = {"Vertices", "Edges", "Triangles", "Quadrilaterals",
                "Tetrahedra", "Hexahedra", "RequiredVertices",
                "RequiredEdges", "RequiredTriangles", "Ridges",
                "Normals", "Tangents", "NormalAtVertices",
                "TangentAtVertices", "End"}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == "Vertices":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                x = float(tokens[i]); y = float(tokens[i + 1])
                z = float(tokens[i + 2])
                i += 4
                points.append((x, y, z))
        elif tok == "Triangles":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                a = int(tokens[i]) - 1
                b = int(tokens[i + 1]) - 1
                c = int(tokens[i + 2]) - 1
                tag = int(tokens[i + 3])
                i += 4
                tris.append((a, b, c))
                tri_tags.append(tag)
        elif tok == "Tetrahedra":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                a = int(tokens[i]) - 1
                b = int(tokens[i + 1]) - 1
                c = int(tokens[i + 2]) - 1
                d = int(tokens[i + 3]) - 1
                tag = int(tokens[i + 4])
                i += 5
                tets.append((a, b, c, d))
                tet_tags.append(tag)
        elif tok == "End":
            break
        else:
            i += 1
            if i >= len(tokens):
                break
            try:
                int(tokens[i]); i += 1
            except ValueError:
                continue
            while i < len(tokens) and tokens[i] not in sections:
                i += 1
    return (np.asarray(points, dtype=np.float64),
            np.asarray(tris, dtype=np.int64),
            np.asarray(tri_tags, dtype=np.int32),
            np.asarray(tets, dtype=np.int64),
            np.asarray(tet_tags, dtype=np.int32))


def _stitch_back(
        in_data: dict,
        in_subdomain: np.ndarray,
        sub_out_points: np.ndarray,
        sub_out_tets: np.ndarray,
        sub_out_tet_tags: np.ndarray,
        sub_out_tris: np.ndarray,
        sub_out_tri_tags: np.ndarray,
        halo_tag: int,
        snap_m: float,
        ) -> dict:
    """Merge mmg3d-modified subdomain back into the input mesh.

    Strategy:
      - Outside-subdomain tets keep their original vertex coords +
        connectivity.
      - In-subdomain tets are REPLACED by mmg3d's output tets.
      - Halo-boundary vertices match by snap_key (the halo lock
        ensures mmg3d preserved them bit-exactly... up to medit
        write/read round-off, which we resolve via snap_key).

    Returns a dict with `points, tris, tri_tags, tets, tet_tags`
    suitable for `meshio.write`.
    """
    points_in = in_data["points"]
    tris_in = in_data["tris"]
    tri_tags_in = in_data["tri_tags"]
    tets_in = in_data["tets"]
    tet_tags_in = in_data["tet_tags"]

    def _snap_key(v: np.ndarray) -> tuple[int, int, int]:
        return (int(round(v[0] / snap_m)),
                int(round(v[1] / snap_m)),
                int(round(v[2] / snap_m)))

    # 1. Collect vertices that are referenced by any OUTSIDE tet.
    outside_tet_idx = np.flatnonzero(~in_subdomain)
    outside_tets = tets_in[outside_tet_idx]
    outside_v_global = np.unique(outside_tets.flatten())
    outside_v_set = set(int(v) for v in outside_v_global)

    # Also keep vertices referenced by any outside-subdomain triangle
    # (box / fault triangles whose 3 vertices are all NOT in
    # subdomain).  We need them in the merged vertex array even if
    # no outside tet references them (rare but possible).
    in_subdomain_v = np.zeros(points_in.shape[0], dtype=bool)
    for ti in np.flatnonzero(in_subdomain):
        in_subdomain_v[tets_in[ti]] = True
    for ti in range(tris_in.shape[0]):
        a, b, c = (int(tris_in[ti, 0]), int(tris_in[ti, 1]),
                    int(tris_in[ti, 2]))
        if (not in_subdomain_v[a] and not in_subdomain_v[b]
                and not in_subdomain_v[c]):
            outside_v_set.update((a, b, c))

    # 2. Build the merged vertex list.  Start with outside vertices
    # (in their original order so outside tet indexing only needs a
    # remap on outside_v_global).
    outside_v_sorted = np.asarray(sorted(outside_v_set), dtype=np.int64)
    merged_points: list[tuple[float, float, float]] = []
    snap_to_merged: dict[tuple[int, int, int], int] = {}
    old_global_to_merged = -np.ones(points_in.shape[0], dtype=np.int64)
    for old_gi in outside_v_sorted:
        coord = tuple(points_in[old_gi])
        merged_points.append(coord)
        old_global_to_merged[old_gi] = len(merged_points) - 1
        snap_to_merged[_snap_key(points_in[old_gi])] = (
            len(merged_points) - 1)

    # 3. Insert mmg3d output vertices.  Halo-boundary vertices match
    # by snap_key with outside vertices; reuse those merged indices.
    # Interior subdomain vertices are new.
    sub_old_to_merged = -np.ones(sub_out_points.shape[0],
                                   dtype=np.int64)
    for li in range(sub_out_points.shape[0]):
        k = _snap_key(sub_out_points[li])
        if k in snap_to_merged:
            sub_old_to_merged[li] = snap_to_merged[k]
        else:
            merged_points.append(tuple(sub_out_points[li]))
            mi = len(merged_points) - 1
            sub_old_to_merged[li] = mi
            snap_to_merged[k] = mi

    # 4. Build merged tets.
    merged_tets: list[tuple[int, int, int, int]] = []
    merged_tet_tags: list[int] = []
    # Outside tets (re-indexed via old_global_to_merged).
    for ti in outside_tet_idx:
        a, b, c, d = (int(tets_in[ti, 0]), int(tets_in[ti, 1]),
                       int(tets_in[ti, 2]), int(tets_in[ti, 3]))
        merged_tets.append((
            int(old_global_to_merged[a]),
            int(old_global_to_merged[b]),
            int(old_global_to_merged[c]),
            int(old_global_to_merged[d])))
        merged_tet_tags.append(int(tet_tags_in[ti]))
    # In-subdomain tets from mmg3d.
    for ti in range(sub_out_tets.shape[0]):
        a, b, c, d = (int(sub_out_tets[ti, 0]),
                       int(sub_out_tets[ti, 1]),
                       int(sub_out_tets[ti, 2]),
                       int(sub_out_tets[ti, 3]))
        merged_tets.append((
            int(sub_old_to_merged[a]),
            int(sub_old_to_merged[b]),
            int(sub_old_to_merged[c]),
            int(sub_old_to_merged[d])))
        merged_tet_tags.append(int(sub_out_tet_tags[ti]))

    # 5. Build merged triangles.
    # Outside triangles: original tri's 3 vertices ALL outside
    # subdomain → keep.  If ALL 3 are in subdomain → drop (mmg3d
    # owns it).  If MIXED → keep (it's a halo / boundary-area tri
    # not modified by mmg3d).
    merged_tris: list[tuple[int, int, int]] = []
    merged_tri_tags: list[int] = []
    for ti in range(tris_in.shape[0]):
        a, b, c = (int(tris_in[ti, 0]), int(tris_in[ti, 1]),
                    int(tris_in[ti, 2]))
        all_in = (in_subdomain_v[a] and in_subdomain_v[b]
                   and in_subdomain_v[c])
        if all_in:
            continue
        merged_tris.append((
            int(old_global_to_merged[a]),
            int(old_global_to_merged[b]),
            int(old_global_to_merged[c])))
        merged_tri_tags.append(int(tri_tags_in[ti]))
    # Subdomain triangles from mmg3d.  Drop the halo-boundary
    # triangles (they only existed to lock the mmg3d boundary; they
    # are not part of the production tag set).
    for ti in range(sub_out_tris.shape[0]):
        tag = int(sub_out_tri_tags[ti])
        if tag == halo_tag:
            continue
        a, b, c = (int(sub_out_tris[ti, 0]),
                    int(sub_out_tris[ti, 1]),
                    int(sub_out_tris[ti, 2]))
        merged_tris.append((
            int(sub_old_to_merged[a]),
            int(sub_old_to_merged[b]),
            int(sub_old_to_merged[c])))
        merged_tri_tags.append(tag)

    # R-006: producer-side topology assertion.  Every 1-tet boundary
    # face of the stitched mesh MUST have a matching tagged triangle.
    # Catches surface-hole regressions locally (e.g., the 2026-05-02
    # all-in-subdomain skip dropping the only tag for a model-bdry
    # face that an outside-cavity tet still uses) before they escape
    # to validate_msh.py downstream.
    face_count: dict[frozenset, int] = {}
    for (a, b, c, d) in merged_tets:
        for face in (frozenset((a, b, c)), frozenset((a, b, d)),
                      frozenset((a, c, d)), frozenset((b, c, d))):
            face_count[face] = face_count.get(face, 0) + 1
    bdry = {f for f, n in face_count.items() if n == 1}
    tagged = {frozenset((a, b, c)) for (a, b, c) in merged_tris}
    holes = bdry - tagged
    print(f"  [R-006 DIAG] merged_tets={len(merged_tets)}, "
          f"merged_tris={len(merged_tris)}, "
          f"bdry={len(bdry)}, tagged={len(tagged)}, "
          f"holes={len(holes)}", file=sys.stderr)
    if holes:
        merged_pts = np.asarray(merged_points)
        ex = []
        for f in list(holes)[:5]:
            vs = list(f)
            cc = merged_pts[vs].mean(axis=0)
            ex.append(f"verts={sorted(vs)} centroid="
                       f"({cc[0]:.0f},{cc[1]:.0f},{cc[2]:.0f})")
        raise RuntimeError(
            f"mmg3d_local_patch._stitch_back: produced "
            f"{len(holes)} unlabeled 1-tet bdry face(s); "
            f"first 5: {'; '.join(ex)}.  This is a stitch bug — "
            f"outside-cavity tet has a face whose 3 vertices are "
            f"all in subdomain but mmg3d's output didn't produce "
            f"a matching tagged tri.")

    return {
        "points":   np.asarray(merged_points, dtype=np.float64),
        "tris":     np.asarray(merged_tris, dtype=np.int64),
        "tri_tags": np.asarray(merged_tri_tags, dtype=np.int32),
        "tets":     np.asarray(merged_tets, dtype=np.int64),
        "tet_tags": np.asarray(merged_tet_tags, dtype=np.int32),
    }


def _canonical_snap_polyline_in_mesh(
        merged: dict,
        polyline_keys: set[tuple[int, int, int]],
        snap_m: float
        ) -> int:
    """Snap every polyline-vertex coord in `merged['points']` to its
    canonical position `snap_key * snap_m`.

    Used after the local-patch stitch when `polyline_relax=True`:
    mmg3d may have introduced machine-epsilon perturbation on
    RequiredVertex coords (and could even re-locate them within
    `-hausd` tolerance).  Snapping back to canonical guarantees
    that EVERY vertex with a polyline snap_key has bit-identical
    coords across faults — restoring the cross-fault conformity
    invariant we deliberately relaxed at the edge level.

    The same mechanism is used by
    `break_fault_wedges.snap_polyline_vertices_to_canonical`
    (the pattern this implementation mirrors).

    Modifies `merged['points']` in place; returns the count of
    vertices actually snapped (== count of polyline vertices that
    had non-canonical coords).
    """
    if not polyline_keys:
        return 0
    P = merged["points"]
    n_snapped = 0
    for vi in range(P.shape[0]):
        key = (int(round(P[vi, 0] / snap_m)),
               int(round(P[vi, 1] / snap_m)),
               int(round(P[vi, 2] / snap_m)))
        if key in polyline_keys:
            canon = np.array(key, dtype=np.float64) * snap_m
            if not np.array_equal(P[vi], canon):
                P[vi] = canon
                n_snapped += 1
    return n_snapped


def _write_msh(out_path: Path, data: dict) -> None:
    """Write a merged mesh dict to gmsh .msh v2.2 ASCII via meshio."""
    import meshio
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mesh = meshio.Mesh(
        points=data["points"],
        cells=[("triangle", data["tris"]),
               ("tetra",    data["tets"])],
        cell_data={"gmsh:physical":
                   [data["tri_tags"], data["tet_tags"]],
                   "gmsh:geometrical":
                   [data["tri_tags"], data["tet_tags"]]},
    )
    meshio.write(out_path, mesh, file_format="gmsh22", binary=False)


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------
def local_patch(in_msh: Path, out_msh: Path,
                gamma_thresh: float = 0.01,
                halo_radius_m: float = 2000.0,
                hmin: float = 50.0,
                hmax: float = 1000.0,
                hgrad: float = 1.1,
                hausd: float = 5.0,
                hgradreq: float = 1.1,
                fault_tri_tag: int = 100,
                halo_tag: int = 200,
                snap_m: float = 0.1,
                provenance_path: Path | None = None,
                transform_path: Path | None = None,
                binary: str = "mmg3d_O3",
                keep_intermediate: bool = False,
                polyline_relax: bool = False,
                ) -> dict:
    """Run a local mmg3d patch on the slivers of `in_msh`, writing
    the stitched result to `out_msh`.  Returns a JSON-serializable
    report.
    """
    if not in_msh.exists():
        raise SystemExit(f"--in-msh not found: {in_msh}")
    if gamma_thresh <= 0 or gamma_thresh >= 1:
        raise ValueError(
            f"gamma_thresh must be in (0, 1); got {gamma_thresh}")
    if halo_radius_m <= 0:
        raise ValueError(
            f"halo_radius_m must be > 0; got {halo_radius_m}")
    if hmin <= 0 or hmax <= 0:
        raise ValueError(f"hmin/hmax must be > 0")
    if hmin > hmax:
        raise ValueError(f"hmin must be <= hmax")
    if hgrad < 1.0:
        raise ValueError(f"hgrad must be >= 1.0")
    if hausd <= 0:
        raise ValueError(f"hausd must be > 0")
    if hgradreq < 1.0:
        raise ValueError(f"hgradreq must be >= 1.0")

    out_msh.parent.mkdir(parents=True, exist_ok=True)
    work_dir = out_msh.parent / "_mmg3d_patch_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Read input.
    in_data = _read_msh(in_msh)
    points, tets = in_data["points"], in_data["tets"]
    print(f"  read {in_msh}: V={points.shape[0]}, T={tets.shape[0]}, "
          f"tris={in_data['tris'].shape[0]}", file=sys.stderr)

    # 2. Identify bad tets.
    g = _gamma_per_tet(points, tets)
    bad_tet_idx = np.flatnonzero(g < gamma_thresh)
    print(f"  bad tets (γ < {gamma_thresh}): {bad_tet_idx.size}",
          file=sys.stderr)
    if bad_tet_idx.size == 0:
        # Nothing to do — copy input to output.
        shutil.copy(in_msh, out_msh)
        if not keep_intermediate:
            shutil.rmtree(work_dir, ignore_errors=True)
        me_noop = _min_edge_of_tets(points, tets)
        return {
            "in_msh":  str(in_msh),
            "out_msh": str(out_msh),
            "n_bad_tets":          0,
            "n_subdomain_tets":    0,
            "gamma_min_in":        float(g.min()) if g.size else 0.0,
            "gamma_min_out":       float(g.min()) if g.size else 0.0,
            "min_edge_in_m":       me_noop,
            "min_edge_out_m":      me_noop,
            "no_op":               True,
        }

    # 3. Subdomain.
    in_subdomain = _identify_subdomain_tets(
        points, tets, bad_tet_idx, halo_radius_m)
    n_sub = int(in_subdomain.sum())
    print(f"  subdomain (within {halo_radius_m} m of any bad tet): "
          f"{n_sub} tets ({100*n_sub/tets.shape[0]:.3f} %)",
          file=sys.stderr)

    # 4. Halo boundary faces.
    face_to_tets = _build_tet_face_adjacency(tets)
    halo_faces = _identify_halo_boundary_faces(
        tets, in_subdomain, face_to_tets)
    print(f"  halo boundary triangles: {len(halo_faces)}",
          file=sys.stderr)

    # 5. Polyline keys (re-use `optim_relax_fault` logic).
    polyline_keys: set[tuple[int, int, int]] | None = None
    if provenance_path is not None and provenance_path.exists():
        # Build the (a, b, c, tag) list that
        # _identify_polyline_keys_from_provenance expects.
        triangles_with_tags = [
            (int(in_data["tris"][ti, 0]),
             int(in_data["tris"][ti, 1]),
             int(in_data["tris"][ti, 2]),
             int(in_data["tri_tags"][ti]))
            for ti in range(in_data["tris"].shape[0])]
        polyline_keys = m3p._identify_polyline_keys_from_provenance(
            points, triangles_with_tags,
            provenance_path, snap_m, fault_tri_tag)
        print(f"  polyline keys: {len(polyline_keys)}",
              file=sys.stderr)

    # 6. Free-surface clearance.
    free_surface_clearance_m = 0.0
    z_top = 0.0
    if transform_path is not None and transform_path.exists():
        tx = json.loads(transform_path.read_text())
        free_surface_clearance_m = float(
            tx.get("free_surface_clearance_m", 0.0))

    # 7. Write subdomain medit.
    sub_in_medit = work_dir / "subdomain.mesh"
    sub_meta = _write_subdomain_medit(
        sub_in_medit,
        points=points,
        in_subdomain=in_subdomain,
        tets=tets,
        tet_tags=in_data["tet_tags"],
        tris=in_data["tris"],
        tri_tags=in_data["tri_tags"],
        halo_faces=halo_faces,
        halo_tag=halo_tag,
        polyline_keys=polyline_keys,
        snap_m=snap_m,
        free_surface_clearance_m=free_surface_clearance_m,
        z_top=z_top,
        fault_tri_tag=fault_tri_tag,
        polyline_relax=polyline_relax,
    )
    print(f"  wrote subdomain medit: {sub_meta['n_subdomain_vertices']} V, "
          f"{sub_meta['n_subdomain_tets']} T, "
          f"{sub_meta['n_subdomain_triangles']} tris "
          f"({sub_meta['n_halo_boundary_triangles']} halo, "
          f"{sub_meta['n_required_triangles']} required, "
          f"polyline_relax={polyline_relax})",
          file=sys.stderr)

    # 8. Run mmg3d on subdomain.
    #
    # `-nosurf` is REQUIRED here.  Without it, mmg3d creates NEW
    # Steiner vertices that can land on box-face planes (e.g.,
    # z = z_top), generating new triangular faces on the model
    # boundary that we cannot tag (mmg3d doesn't know they should
    # be tag-5).  Stitch-back then leaves these as unlabeled
    # boundary faces — surface holes.  Empirical: `-nosurf` brings
    # the post-patch hole count from 3 to 0 on Step 6.
    #
    # `-nosurf` does prevent some quality optimizations on the
    # fault surface (mmg3d can no longer collapse short fault
    # edges).  For polyline-relax=True users who explicitly want
    # mmg3d to modify the fault surface, the polyline edges are
    # still freed via the RequiredEdges relaxation; the surface
    # itself stays bit-exact.  The trade-off is acceptable: surface
    # integrity > polyline-edge-collapse opportunity.
    sub_out_medit = work_dir / "subdomain_post.mesh"
    log_path = work_dir / "mmg3d.log"
    mode_flags = ["-optim", "-opnbdy", "-nosurf",
                   "-hgradreq", repr(float(hgradreq))]
    m3p._run_mmg3d(sub_in_medit, sub_out_medit,
                    hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
                    extra_args=mode_flags,
                    binary=binary, log_path=log_path)

    # 9. Read mmg3d output.
    (sub_out_points, sub_out_tris, sub_out_tri_tags,
     sub_out_tets, sub_out_tet_tags) = _read_subdomain_medit(sub_out_medit)
    print(f"  mmg3d output: V={sub_out_points.shape[0]}, "
          f"T={sub_out_tets.shape[0]}, "
          f"tris={sub_out_tris.shape[0]}", file=sys.stderr)

    # 10. Stitch back.
    merged = _stitch_back(
        in_data, in_subdomain,
        sub_out_points, sub_out_tets, sub_out_tet_tags,
        sub_out_tris, sub_out_tri_tags,
        halo_tag=halo_tag, snap_m=snap_m)

    # 10b. Canonical-snap of polyline vertices.  In polyline_relax
    # mode mmg3d may have collapsed a polyline edge OR perturbed a
    # RequiredVertex coord at machine epsilon.  Even outside relax
    # mode, mmg3d's medit write/read round-trip can introduce
    # sub-femtometer drift on vertices.  Snapping every polyline-
    # snap-key vertex to canonical position guarantees cross-fault
    # conformity is bit-identical in the output.
    n_snapped = 0
    if polyline_keys:
        n_snapped = _canonical_snap_polyline_in_mesh(
            merged, polyline_keys, snap_m)
        print(f"  canonical-snap polyline vertices: "
              f"{n_snapped} of {merged['points'].shape[0]} "
              f"vertices snapped to canonical positions",
              file=sys.stderr)

    # 11. Write merged mesh.
    _write_msh(out_msh, merged)
    print(f"  wrote merged mesh: V={merged['points'].shape[0]}, "
          f"T={merged['tets'].shape[0]}, "
          f"tris={merged['tris'].shape[0]}", file=sys.stderr)

    # 12. Report.
    g_in_min = float(g.min())
    me_in = _min_edge_of_tets(points, tets)
    g_out = _gamma_per_tet(merged["points"], merged["tets"])
    g_out_min = float(g_out.min()) if g_out.size else 0.0
    me_out = _min_edge_of_tets(merged["points"], merged["tets"])
    n_below_thresh_out = int((g_out < gamma_thresh).sum())
    print(f"  γ_min: {g_in_min:.3e} → {g_out_min:.3e}  "
          f"min_edge: {me_in:.3f} m → {me_out:.3f} m",
          file=sys.stderr)

    if not keep_intermediate:
        shutil.rmtree(work_dir, ignore_errors=True)

    return {
        "in_msh":  str(in_msh),
        "out_msh": str(out_msh),
        "gamma_thresh":          gamma_thresh,
        "halo_radius_m":         halo_radius_m,
        "hmin": hmin, "hmax": hmax, "hgrad": hgrad,
        "hausd": hausd, "hgradreq": hgradreq,
        "fault_tri_tag":         fault_tri_tag,
        "halo_tag":              halo_tag,
        "polyline_relax":        bool(polyline_relax),
        "n_bad_tets":            int(bad_tet_idx.size),
        "n_subdomain_tets":      int(n_sub),
        "n_halo_boundary_tris":  int(len(halo_faces)),
        "n_polyline_canonical_snapped": int(n_snapped),
        "subdomain_meta":        {
            "n_vertices": sub_meta["n_subdomain_vertices"],
            "n_tets":     sub_meta["n_subdomain_tets"],
            "n_required_triangles": sub_meta["n_required_triangles"],
            "n_required_edges":     sub_meta["n_required_edges"],
            "n_required_vertices":  sub_meta["n_required_vertices"],
        },
        "gamma_min_in":          g_in_min,
        "gamma_min_out":         g_out_min,
        "min_edge_in_m":         me_in,
        "min_edge_out_m":        me_out,
        "n_below_thresh_out":    n_below_thresh_out,
        "tets_in":               int(tets.shape[0]),
        "tets_out":              int(merged["tets"].shape[0]),
        "no_op":                 False,
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Local sliver-targeted mmg3d patch pass.")
    p.add_argument("--in-msh", required=True, type=Path)
    p.add_argument("--out-msh", required=True, type=Path)
    p.add_argument("--gamma-thresh", type=float, default=0.01,
                   help="Tets with γ below this threshold are "
                        "candidates for local patch (default 0.01).")
    p.add_argument("--halo-radius-m", type=float, default=2000.0,
                   help="Halo radius around each bad tet's centroid "
                        "(meters).  Tets with centroid within this "
                        "distance of any bad tet are included in "
                        "the subdomain.  Default 2000 m.")
    p.add_argument("--hmin", type=float, default=50.0,
                   help="mmg3d -hmin (m).  Default 50 — allows mmg3d "
                        "to collapse polyline-locked needle edges.")
    p.add_argument("--hmax", type=float, default=1000.0,
                   help="mmg3d -hmax (m).  Default 1000 — match the "
                        "in-tube target.")
    p.add_argument("--hgrad", type=float, default=1.1,
                   help="mmg3d -hgrad.  Default 1.1 — tight gradient.")
    p.add_argument("--hausd", type=float, default=5.0,
                   help="mmg3d -hausd (m).  Default 5 — tight "
                        "geometry-deviation bound.")
    p.add_argument("--hgradreq", type=float, default=1.1,
                   help="mmg3d -hgradreq (gradient from required "
                        "entities).  Default 1.1.")
    p.add_argument("--fault-tri-tag", type=int, default=100)
    p.add_argument("--halo-tag", type=int, default=200,
                   help="Triangle tag assigned to halo-boundary "
                        "triangles (RequiredTriangles in mmg3d "
                        "input).  Must NOT collide with any "
                        "existing physical tag in the input.  "
                        "Default 200.")
    p.add_argument("--snap-m", type=float, default=0.1)
    p.add_argument("--provenance-json", type=Path, default=None)
    p.add_argument("--transform-json", type=Path, default=None)
    p.add_argument("--polyline-relax", action="store_true",
                   help="Drop polyline-edge protection in the "
                        "subdomain so mmg3d can collapse / split "
                        "polyline edges (essential to eliminate "
                        "sub-meter polyline-locked needle edges).  "
                        "Polyline VERTICES remain RequiredVertices "
                        "(pinned at canonical positions).  After "
                        "stitch-back, every polyline-snap-key "
                        "vertex is canonicalized to "
                        "`snap_key * snap_m` to restore "
                        "cross-fault conformity bit-exactly.")
    p.add_argument("--mmg3d-binary", default="mmg3d_O3")
    p.add_argument("--keep-intermediate", action="store_true")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = local_patch(
        in_msh=args.in_msh, out_msh=args.out_msh,
        gamma_thresh=args.gamma_thresh,
        halo_radius_m=args.halo_radius_m,
        hmin=args.hmin, hmax=args.hmax,
        hgrad=args.hgrad, hausd=args.hausd, hgradreq=args.hgradreq,
        fault_tri_tag=args.fault_tri_tag,
        halo_tag=args.halo_tag,
        snap_m=args.snap_m,
        provenance_path=args.provenance_json,
        transform_path=args.transform_json,
        binary=args.mmg3d_binary,
        keep_intermediate=args.keep_intermediate,
        polyline_relax=args.polyline_relax,
    )

    report_path = args.report_json or (
        args.out_msh.parent / "mmg3d_local_patch_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    print(f"  γ_min: {report['gamma_min_in']:.3e} → "
          f"{report['gamma_min_out']:.3e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
