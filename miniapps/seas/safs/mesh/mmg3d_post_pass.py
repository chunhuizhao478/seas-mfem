#!/usr/bin/env python
"""Post-HXT volume-mesh sliver fixer via the MMG Platform's `mmg3d_O3`.

Tier 1 of `REVIEW_interior_subdivision_and_collapse.md`, fallback path:
when surface-level remeshing is blocked (per-fault mmgs cannot see
other faults' surfaces, so its inserted vertices intersect neighbour
faults), run mmg3d on the full HXT output instead.  mmg3d sees the
entire volume mesh — including ALL fault surfaces simultaneously —
and can flip / collapse / Steiner-insert to remove sliver tets while
preserving fault triangles as RequiredTriangles.

Pipeline position:
    ... → generate_safs_mesh → THIS → validate_msh

What mmg3d does (with fault triangles marked Required):
    - Insert Steiner points inside bad-shape tets (sliver elimination).
    - Collapse short interior edges.
    - Flip interior faces / edges to improve dihedral angles.
    - Smooth interior vertex locations (geometry-bounded by `-hausd`).
    - DOES NOT modify fault triangles (RequiredTriangles), free
      surface (z = z_top), or domain box faces.

Cross-fault conformity: mmg3d operates on the post-HXT mesh where
fault triangles are already shared between two tets (validator
check_5).  Marking them Required preserves this two-tet adjacency
through the post-pass.

Usage:
    python mmg3d_post_pass.py --in-msh safs.msh --out-msh safs_mmg.msh \\
        [--hmin 100] [--hmax 25000] [--hgrad 1.3] [--hausd 50]

Requires `mmg3d_O3` on PATH (install via
`conda install -c conda-forge mmgsuite`).
"""
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np


def _run_mmg3d(in_mesh: Path, out_mesh: Path,
                hmin: float, hmax: float, hgrad: float, hausd: float,
                extra_args: list[str] | None = None,
                binary: str = "mmg3d_O3",
                log_path: Path | None = None) -> None:
    """Invoke `mmg3d_O3` on a medit input.

    Always passes `-nr` so that surface ridges are determined by the
    explicit RequiredEdges block (or absent — mmg3d treats all
    boundary edges with two distinct triangle references as ridges
    by default if -nr is not passed; with -nr we rely entirely on
    explicit Required* markers).
    """
    cmd = [binary,
           "-in", str(in_mesh),
           "-out", str(out_mesh),
           "-hmin", repr(float(hmin)),
           "-hmax", repr(float(hmax)),
           "-hgrad", repr(float(hgrad)),
           "-hausd", repr(float(hausd))]
    if extra_args:
        cmd.extend(extra_args)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w") as fh:
            fh.write("CMD: " + " ".join(cmd) + "\n\n")
            fh.flush()
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        msg = (f"mmg3d_O3 failed (rc={r.returncode}) on {in_mesh}; "
               f"see log at {log_path}" if log_path else
               f"mmg3d_O3 failed (rc={r.returncode}) on {in_mesh}")
        raise SystemExit(msg)


def _identify_polyline_keys_from_provenance(
        points: np.ndarray,
        triangles: list[tuple[int, int, int, int]],
        provenance_path: Path,
        snap_m: float,
        fault_tri_tag: int = 100
        ) -> set[tuple[int, int, int]]:
    """Return the set of vertex snap_keys that appear in ≥ 2 faults
    (i.e., cross-fault polyline vertices).

    Reads `fault_provenance.json` (the schema written by
    write_fault_provenance.py) which maps each source fault to a
    contiguous range of triangle indices in the .msh's tag-100
    triangle list.  A vertex is a polyline vertex if it is referenced
    by triangles from ≥ 2 distinct source faults.
    """
    if not provenance_path.exists():
        return set()
    prov = json.loads(provenance_path.read_text())
    faults_dict = prov.get("faults", {})
    if not faults_dict:
        return set()

    # Filter to fault-tag triangles only, in the same order as the .msh.
    fault_tris: list[tuple[int, int, int]] = [
        (a, b, c) for (a, b, c, tag) in triangles if tag == fault_tri_tag]

    if not fault_tris:
        return set()

    # R-002 fix: build an explicit global-index → local-tag100-index
    # map by walking `triangles` in order.  This handles non-contiguous
    # tag-100 placement (e.g., chained mmg3d output where triangles
    # may interleave tags) without relying on a single contiguous
    # offset.  The previous implementation assumed all tag-100
    # triangles formed one contiguous block starting at `offset`,
    # which silently dropped polyline vertices when the assumption
    # broke.
    global_to_local: dict[int, int] = {}
    local_i = 0
    for gi, (_a, _b, _c, tag) in enumerate(triangles):
        if tag == fault_tri_tag:
            global_to_local[gi] = local_i
            local_i += 1

    def _snap_key(v: np.ndarray) -> tuple[int, int, int]:
        return (int(round(v[0] / snap_m)),
                int(round(v[1] / snap_m)),
                int(round(v[2] / snap_m)))

    # vertex snap-key -> set of fault names that touch it.
    key_to_faults: dict[tuple[int, int, int], set[str]] = {}
    for fault_name, entry in faults_dict.items():
        msh_indices = entry.get("triangle_indices_in_msh", [])
        for gi in msh_indices:
            li = global_to_local.get(gi)
            if li is None or li >= len(fault_tris):
                continue
            tri = fault_tris[li]
            for vi in tri:
                k = _snap_key(points[vi])
                key_to_faults.setdefault(k, set()).add(fault_name)

    # Hard-fail sanity check: if provenance lists faults but we found
    # ZERO assigned vertices, the indexing is mismatched.  Raise
    # rather than silently returning an empty polyline_keys set
    # (which would cause mmg3d to under-protect the polyline).
    n_assigned = sum(len(s) for s in key_to_faults.values())
    if n_assigned == 0 and faults_dict:
        raise RuntimeError(
            f"polyline detection produced ZERO assigned vertices "
            f"despite {len(faults_dict)} faults in provenance.  "
            f"Likely cause: `triangle_indices_in_msh` does not match "
            f"the layout of `triangles` argument.  Re-emit "
            f"fault_provenance.json against the same .msh that "
            f"mmg3d_post_pass is consuming.")

    return {k for k, fs in key_to_faults.items() if len(fs) >= 2}


def _convert_msh_to_medit(msh_path: Path, mesh_path: Path,
                            fault_tri_tag: int = 100,
                            mode: str = "optim",
                            provenance_path: Path | None = None,
                            transform_path: Path | None = None,
                            snap_m: float = 0.1,
                            free_surface_clearance_m: float | None = None
                            ) -> dict:
    """Convert a gmsh `.msh` to mmg's medit `.mesh` format using meshio.

    For `mode="optim"`/`"adapt"` (default): all triangles tagged
    `fault_tri_tag` (100) are marked RequiredTriangles, so mmg3d
    preserves the fault interface bit-exactly through its
    sliver-removal pass.

    For `mode="optim_relax_fault"`: fault triangles are NOT marked
    Required.  Instead, the cross-fault polyline (edges shared
    between faults), the free-surface trace (edges where both
    endpoints have z within `free_surface_clearance_m` of the top of
    the model), and the fault-perimeter (edges incident to only one
    fault triangle = the manifold boundary of each fault patch) are
    marked RequiredEdges.  Polyline + free-surface vertices are
    marked RequiredVertices.  This relaxation lets mmg3d collapse
    needle fault triangles and Steiner-insert in the fault interior
    while preserving the SEAS-relevant fault geometry (perimeter +
    polyline + free-surface trace).

    Box-face triangles (tags 1-6) are kept and NOT marked Required —
    mmg3d preserves them as the domain boundary by default with
    `-opnbdy`.
    """
    import meshio  # local import: optional dep at module level
    m = meshio.read(msh_path)

    points = np.asarray(m.points, dtype=np.float64)
    if points.shape[1] != 3:
        raise ValueError(
            f"expected 3-D points; got shape {points.shape}")

    # Collect triangles + tetrahedra and their physical tags.
    tris: list[tuple[int, int, int, int]] = []
    tets: list[tuple[int, int, int, int, int]] = []

    for cb, tags in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        if cb.type == "triangle":
            for tri, tag in zip(cb.data, tags):
                tris.append((int(tri[0]), int(tri[1]),
                              int(tri[2]), int(tag)))
        elif cb.type == "tetra":
            for tet, tag in zip(cb.data, tags):
                tets.append((int(tet[0]), int(tet[1]),
                              int(tet[2]), int(tet[3]), int(tag)))

    if not tets:
        raise ValueError(f"no tetrahedra in {msh_path}; nothing to fix")

    # Identify required triangle indices (1-based for medit).  For
    # the relaxed-fault mode, we do NOT mark fault triangles as
    # Required — mmg3d is allowed to modify them (collapse needles,
    # Steiner-insert, edge-flip), constrained only by per-edge /
    # per-vertex Required markers.
    required_tri_idx: list[int] = []
    if mode in ("optim", "adapt"):
        for i, (_, _, _, tag) in enumerate(tris):
            if tag == fault_tri_tag:
                required_tri_idx.append(i + 1)  # 1-based

    # Collect protected (Required) edges + vertices for
    # `optim_relax_fault` mode.  An edge is identified by its sorted
    # vertex-index pair (min, max).
    required_edge_pairs: set[tuple[int, int]] = set()
    required_vertex_set: set[int] = set()
    polyline_keys: set[tuple[int, int, int]] = set()
    n_polyline_edges = 0
    n_freesurface_edges = 0
    n_perimeter_edges = 0

    if mode == "optim_relax_fault":
        # 1. Cross-fault polyline edges: extracted via
        #    fault_provenance.json (vertices touching ≥ 2 faults).
        if provenance_path is None:
            raise ValueError(
                "mode='optim_relax_fault' requires provenance_path "
                "(read from fault_provenance.json) to identify "
                "cross-fault polyline vertices")
        polyline_keys = _identify_polyline_keys_from_provenance(
            points, tris, provenance_path, snap_m, fault_tri_tag)

        def _snap_key(v: np.ndarray) -> tuple[int, int, int]:
            return (int(round(v[0] / snap_m)),
                    int(round(v[1] / snap_m)),
                    int(round(v[2] / snap_m)))

        # vertex_idx -> True iff snap_key is in polyline_keys.
        is_polyline_vertex = np.array(
            [_snap_key(points[vi]) in polyline_keys
             for vi in range(points.shape[0])], dtype=bool)
        for vi in np.flatnonzero(is_polyline_vertex):
            required_vertex_set.add(int(vi))

        # 2. Free-surface trace: vertices within clearance of the
        #    free surface (top of the model, z = z_top).  z_top is
        #    typically 0 for SAFS; clearance from transform.json.
        z_top = 0.0
        clearance = float(free_surface_clearance_m or 0.0)
        if transform_path is not None and transform_path.exists():
            tx = json.loads(transform_path.read_text())
            if free_surface_clearance_m is None:
                clearance = float(
                    tx.get("free_surface_clearance_m", clearance))
            # `z_top` is the model's top elevation in local frame;
            # for SAFS it's 0.0 (free surface).  No transform.json
            # field exposes it directly, so we keep z_top = 0.0
            # (the SAFS convention documented in safs.geo line 16).
        # R-005 fix: build the set of vertices touched by FAULT
        # triangles only.  The "free-surface trace" is the curve
        # where the fault meets z=0; box-top vertices (tag 5) at
        # z=0 are a different thing and must NOT be flagged as
        # free-surface for the polyline-protection logic that
        # consumes this array.  Without the fault-only filter,
        # box-top vertices that happen to belong to a fault triangle
        # (because the fault triangle reaches z≈0) get flagged
        # alongside legitimate free-surface trace vertices,
        # producing extra RequiredEdges that block mmg3d
        # optimization in the bulk near the box top.
        fault_vertex_set: set[int] = set()
        for (a, b, c, tag) in tris:
            if tag == fault_tri_tag:
                fault_vertex_set.update((a, b, c))
        is_freesurface_vertex = np.array(
            [(vi in fault_vertex_set
              and points[vi, 2] >= z_top - clearance)
             for vi in range(points.shape[0])], dtype=bool)

        # 3. Per-fault perimeter: an edge is a fault-perimeter edge
        #    if it is incident to exactly ONE tag=100 triangle in the
        #    .msh (i.e., a boundary edge of the fault manifold).
        #    Walk the fault triangles and count occurrences per edge.
        fault_edge_count: dict[tuple[int, int], int] = {}
        for (a, b, c, tag) in tris:
            if tag != fault_tri_tag:
                continue
            for (u, v) in ((a, b), (b, c), (c, a)):
                key = (u, v) if u < v else (v, u)
                fault_edge_count[key] = fault_edge_count.get(key, 0) + 1
        perimeter_edges: set[tuple[int, int]] = {
            ek for ek, n in fault_edge_count.items() if n == 1}

        # 4. Fault-fault polyline edges: edges of fault triangles
        #    where BOTH endpoints are polyline vertices.
        polyline_edges: set[tuple[int, int]] = set()
        for (a, b, c, tag) in tris:
            if tag != fault_tri_tag:
                continue
            for (u, v) in ((a, b), (b, c), (c, a)):
                if (is_polyline_vertex[u]
                        and is_polyline_vertex[v]):
                    polyline_edges.add(
                        (u, v) if u < v else (v, u))

        # 5. Free-surface trace edges: edges of fault triangles
        #    where BOTH endpoints are free-surface vertices.
        freesurface_edges: set[tuple[int, int]] = set()
        for (a, b, c, tag) in tris:
            if tag != fault_tri_tag:
                continue
            for (u, v) in ((a, b), (b, c), (c, a)):
                if (is_freesurface_vertex[u]
                        and is_freesurface_vertex[v]):
                    freesurface_edges.add(
                        (u, v) if u < v else (v, u))

        # Union all three protected edge classes.  Required-vertex
        # set: polyline vertices are required (must remain at
        # canonical position; free-surface vertices may slide
        # along the surface, so we DON'T mark them Required).
        n_polyline_edges    = len(polyline_edges)
        n_freesurface_edges = len(freesurface_edges)
        n_perimeter_edges   = len(perimeter_edges)
        required_edge_pairs = (polyline_edges
                                | freesurface_edges
                                | perimeter_edges)

    # Convert RequiredVertices to 1-based, sorted list.
    required_vertex_idx = sorted(vi + 1 for vi in required_vertex_set)
    # Convert RequiredEdges to a 1-based ordered list of (a, b)
    # plus a parallel list of edge-block 1-based indices.
    required_edges_list = sorted(required_edge_pairs)

    mesh_path.parent.mkdir(parents=True, exist_ok=True)
    with mesh_path.open("w") as f:
        f.write("MeshVersionFormatted 2\n\nDimension 3\n\n")

        f.write(f"Vertices\n{points.shape[0]}\n")
        for v in points:
            f.write(f"{v[0]:+.17e} {v[1]:+.17e} {v[2]:+.17e} 0\n")
        f.write("\n")

        f.write(f"Triangles\n{len(tris)}\n")
        for (a, b, c, tag) in tris:
            f.write(f"{a+1} {b+1} {c+1} {tag}\n")
        f.write("\n")

        if required_tri_idx:
            f.write(f"RequiredTriangles\n{len(required_tri_idx)}\n")
            for ti in required_tri_idx:
                f.write(f"{ti}\n")
            f.write("\n")

        if required_edges_list:
            f.write(f"Edges\n{len(required_edges_list)}\n")
            for (u, v) in required_edges_list:
                f.write(f"{u+1} {v+1} 0\n")
            f.write("\n")
            f.write(f"RequiredEdges\n{len(required_edges_list)}\n")
            for i in range(len(required_edges_list)):
                f.write(f"{i+1}\n")
            f.write("\n")

        if required_vertex_idx:
            f.write(f"RequiredVertices\n{len(required_vertex_idx)}\n")
            for vi in required_vertex_idx:
                f.write(f"{vi}\n")
            f.write("\n")

        f.write(f"Tetrahedra\n{len(tets)}\n")
        for (a, b, c, d, tag) in tets:
            f.write(f"{a+1} {b+1} {c+1} {d+1} {tag}\n")
        f.write("\n")

        f.write("End\n")

    return {
        "n_vertices":   int(points.shape[0]),
        "n_triangles":  len(tris),
        "n_tets":       len(tets),
        "n_required_triangles": len(required_tri_idx),
        "n_required_edges":     len(required_edges_list),
        "n_required_vertices":  len(required_vertex_idx),
        "n_polyline_edges":     int(n_polyline_edges),
        "n_freesurface_edges":  int(n_freesurface_edges),
        "n_perimeter_edges":    int(n_perimeter_edges),
        "n_polyline_keys":      int(len(polyline_keys)),
        "fault_tri_tag":        int(fault_tri_tag),
        "mode":                 mode,
    }


def _convert_medit_to_msh(mesh_path: Path, msh_path: Path,
                            expected_tri_tags: set[int] | None = None,
                            expected_tet_tags: set[int] | None = None,
                            ) -> dict:
    """Convert mmg3d's medit `.mesh` output back to gmsh `.msh` v2.2,
    preserving the physical tags via `gmsh:physical` cell_data.

    Triangle tags and tet tags survive in mmg3d's output medit "ref"
    field, so we recover them and re-emit a tagged .msh consumable
    by validate_msh.py / convert_msh.py.

    R-008 defense: if `expected_tri_tags` and / or `expected_tet_tags`
    are provided, the output's tags must be a SUBSET of the expected
    set.  mmg3d should not introduce new physical tags; if it does,
    raise rather than silently propagating an unexpected tag that
    would cause `validate_msh.py`'s check_1 to fail with a
    confusing message.
    """
    import meshio
    points: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    tri_tags: list[int] = []
    tets: list[tuple[int, int, int, int]] = []
    tet_tags: list[int] = []

    tokens: list[str] = []
    with mesh_path.open() as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
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
                x = float(tokens[i]); y = float(tokens[i+1])
                z = float(tokens[i+2])
                i += 4
                points.append((x, y, z))
        elif tok == "Triangles":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                a = int(tokens[i]) - 1
                b = int(tokens[i+1]) - 1
                c = int(tokens[i+2]) - 1
                tag = int(tokens[i+3])
                i += 4
                tris.append((a, b, c))
                tri_tags.append(tag)
        elif tok == "Tetrahedra":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                a = int(tokens[i]) - 1
                b = int(tokens[i+1]) - 1
                c = int(tokens[i+2]) - 1
                d = int(tokens[i+3]) - 1
                tag = int(tokens[i+4])
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

    # R-008 defense: validate that mmg3d did not introduce any
    # physical tag not present in the input.  Such tags would
    # silently propagate into the output .msh and cause
    # `validate_msh.py` check_1 to fail with a confusing message.
    if expected_tri_tags is not None:
        out_tri_tags = set(tri_tags)
        new_tri_tags = out_tri_tags - expected_tri_tags
        if new_tri_tags:
            raise RuntimeError(
                f"mmg3d output introduced unknown triangle tags "
                f"{sorted(new_tri_tags)}; expected subset of "
                f"{sorted(expected_tri_tags)}.  Investigate "
                f"{mesh_path}")
    if expected_tet_tags is not None:
        out_tet_tags = set(tet_tags)
        new_tet_tags = out_tet_tags - expected_tet_tags
        if new_tet_tags:
            raise RuntimeError(
                f"mmg3d output introduced unknown tet tags "
                f"{sorted(new_tet_tags)}; expected subset of "
                f"{sorted(expected_tet_tags)}.  Investigate "
                f"{mesh_path}")

    mesh = meshio.Mesh(
        points=np.asarray(points, dtype=np.float64),
        cells=[("triangle", np.asarray(tris, dtype=np.int64)),
               ("tetra",    np.asarray(tets, dtype=np.int64))],
        cell_data={"gmsh:physical":
                   [np.asarray(tri_tags, dtype=np.int32),
                    np.asarray(tet_tags, dtype=np.int32)],
                   "gmsh:geometrical":
                   [np.asarray(tri_tags, dtype=np.int32),
                    np.asarray(tet_tags, dtype=np.int32)]},
    )
    msh_path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(msh_path, mesh, file_format="gmsh22", binary=False)
    return {
        "n_points":     len(points),
        "n_triangles":  len(tris),
        "n_tets":       len(tets),
    }


def _resolve_mode_flags(mode: str,
                          extra_mmg_args: list[str] | None
                          ) -> list[str]:
    """Compute the mmg3d invocation flags for `mode`."""
    if mode == "optim":
        mode_flags = ["-optim", "-nosurf", "-opnbdy"]
    elif mode == "adapt":
        mode_flags = ["-opnbdy"]
    elif mode == "optim_relax_fault":
        # See post_pass() docstring for the comment block on
        # `-optim -opnbdy -hgradreq 1.3` rationale.
        mode_flags = ["-optim", "-opnbdy", "-hgradreq", "1.3"]
    else:
        raise ValueError(
            f"unknown mode {mode!r}")
    if extra_mmg_args:
        mode_flags = mode_flags + list(extra_mmg_args)
    return mode_flags


def _min_edge_of_msh(msh_path: Path) -> float:
    """Global minimum tet edge length (m) for a .msh file.  Returns
    -1.0 for empty / unreadable input.
    """
    import meshio
    m = meshio.read(msh_path)
    tets = None
    for cb in m.cells:
        if cb.type == "tetra":
            tets = cb.data
            break
    if tets is None or tets.shape[0] == 0:
        return -1.0
    P = np.asarray(m.points, dtype=np.float64)
    p0 = P[tets[:, 0]]; p1 = P[tets[:, 1]]
    p2 = P[tets[:, 2]]; p3 = P[tets[:, 3]]
    e_sq = np.stack([
        np.sum((p1 - p0) ** 2, axis=1),
        np.sum((p2 - p0) ** 2, axis=1),
        np.sum((p3 - p0) ** 2, axis=1),
        np.sum((p2 - p1) ** 2, axis=1),
        np.sum((p3 - p1) ** 2, axis=1),
        np.sum((p3 - p2) ** 2, axis=1),
    ], axis=1)
    return float(np.sqrt(e_sq.min()))


def _gamma_min_of_msh(msh_path: Path,
                       fault_tri_tag: int = 100) -> float:
    """Compute the worst-tet γ over a .msh file's tetrahedra.

    Uses the same metric as `validate_msh.py` check_10
    (Mesh.QualityType=2, isotropic gamma normalized to 1).  Returns
    the minimum over all tets; for an empty mesh, returns -1.0.

    Used by the best-of-N trial selector to pick the best output
    across mmg3d's non-deterministic runs.
    """
    import meshio
    m = meshio.read(msh_path)
    tets = None
    for cb in m.cells:
        if cb.type == "tetra":
            tets = cb.data
            break
    if tets is None or tets.shape[0] == 0:
        return -1.0
    P = np.asarray(m.points, dtype=np.float64)
    g_min = np.inf
    for i in range(tets.shape[0]):
        p0, p1, p2, p3 = (P[tets[i, 0]], P[tets[i, 1]],
                           P[tets[i, 2]], P[tets[i, 3]])
        e_lens_sq = (
            np.sum((p1 - p0) ** 2) + np.sum((p2 - p0) ** 2)
            + np.sum((p3 - p0) ** 2) + np.sum((p2 - p1) ** 2)
            + np.sum((p3 - p1) ** 2) + np.sum((p3 - p2) ** 2)
        )
        vol = abs(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
        if vol == 0.0 or e_lens_sq == 0.0:
            return 0.0
        g = 12.0 * (3.0 * vol) ** (2.0 / 3.0) / e_lens_sq
        if g < g_min:
            g_min = g
    return float(g_min)


def _run_one_mmg3d_pass(
        in_msh: Path, out_msh: Path,
        work_dir: Path, log_suffix: str,
        hmin: float, hmax: float, hgrad: float, hausd: float,
        fault_tri_tag: int, mode: str,
        provenance_path: Path | None,
        transform_path: Path | None,
        snap_m: float,
        binary: str,
        mode_flags: list[str]
        ) -> tuple[dict, dict]:
    """Run a single mmg3d pass: msh → medit → mmg3d → medit → msh.

    Returns (input_counts, output_counts) dicts for the report.
    Captures the input's tag set and validates the output via R-008.
    """
    in_medit = work_dir / (in_msh.stem + f"_{log_suffix}_in.mesh")
    out_medit = work_dir / (in_msh.stem + f"_{log_suffix}_post.mesh")
    log_path = work_dir / (in_msh.stem + f"_{log_suffix}_mmg3d.log")

    convert_in = _convert_msh_to_medit(
        in_msh, in_medit,
        fault_tri_tag=fault_tri_tag,
        mode=mode,
        provenance_path=provenance_path,
        transform_path=transform_path,
        snap_m=snap_m)

    # Capture input tag sets for R-008 output validation.
    import meshio
    m_in = meshio.read(in_msh)
    in_tri_tags: set[int] = set()
    in_tet_tags: set[int] = set()
    for cb, tags in zip(m_in.cells,
                         m_in.cell_data.get("gmsh:physical", [])):
        if cb.type == "triangle":
            in_tri_tags.update(int(t) for t in tags)
        elif cb.type == "tetra":
            in_tet_tags.update(int(t) for t in tags)

    _run_mmg3d(in_medit, out_medit,
                hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
                extra_args=mode_flags,
                binary=binary, log_path=log_path)

    convert_out = _convert_medit_to_msh(
        out_medit, out_msh,
        expected_tri_tags=in_tri_tags,
        expected_tet_tags=in_tet_tags)
    return convert_in, convert_out


def post_pass(in_msh: Path, out_msh: Path,
              hmin: float = 100.0, hmax: float = 25000.0,
              hgrad: float = 1.3, hausd: float = 50.0,
              fault_tri_tag: int = 100,
              binary: str = "mmg3d_O3",
              keep_intermediate: bool = False,
              extra_mmg_args: list[str] | None = None,
              mode: str = "optim",
              provenance_path: Path | None = None,
              transform_path: Path | None = None,
              snap_m: float = 0.1,
              n_passes: int = 1,
              n_trials: int = 1
              ) -> dict:
    """Run mmg3d on `in_msh`, preserving fault triangles (tag
    `fault_tri_tag`) as RequiredTriangles, and write the result to
    `out_msh` in gmsh .msh v2.2 format.

    Multi-pass / multi-trial behaviour:

    - `n_passes >= 1`: run mmg3d sequentially, feeding pass i's
      output as pass (i+1)'s input.  Each pass's medit conversion is
      fresh, so RequiredEdges/RequiredVertices reflect the current
      polyline geometry.  R-003 fix.
    - `n_trials >= 1`: repeat the multi-pass sequence N times
      independently and pick the trial with the highest γ_min.
      Exploits mmg3d's run-to-run non-determinism.  R-006 fix.

    The total number of mmg3d invocations is `n_passes * n_trials`.

    Returns a JSON-serializable report.
    """
    if mode not in ("optim", "adapt", "optim_relax_fault"):
        raise ValueError(
            f"mode must be 'optim', 'adapt', or 'optim_relax_fault'; "
            f"got {mode!r}")
    if mode == "optim_relax_fault" and provenance_path is None:
        raise ValueError(
            "mode='optim_relax_fault' requires --provenance-json "
            "(reads fault_provenance.json to identify the cross-fault "
            "polyline vertices that must be protected as "
            "RequiredVertices)")
    if not in_msh.exists():
        raise SystemExit(f"--in-msh not found: {in_msh}")
    if hmin <= 0 or hmax <= 0:
        raise ValueError(f"hmin/hmax must be > 0; got hmin={hmin} hmax={hmax}")
    if hmin > hmax:
        raise ValueError(f"hmin must be <= hmax; got hmin={hmin} hmax={hmax}")
    if hgrad < 1.0:
        raise ValueError(f"hgrad must be >= 1.0; got {hgrad}")
    if hausd <= 0:
        raise ValueError(f"hausd must be > 0; got {hausd}")
    if n_passes < 1:
        raise ValueError(f"n_passes must be >= 1; got {n_passes}")
    if n_trials < 1:
        raise ValueError(f"n_trials must be >= 1; got {n_trials}")
    out_msh.parent.mkdir(parents=True, exist_ok=True)

    work_dir = out_msh.parent / "_mmg3d_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    mode_flags = _resolve_mode_flags(mode, extra_mmg_args)

    # Multi-trial loop.  Each trial is an independent sequence of
    # `n_passes` mmg3d invocations.  We pick the trial whose final
    # γ_min is highest.  R-006 fix.
    best: dict | None = None
    trials_log: list[dict] = []
    for trial_idx in range(n_trials):
        # Multi-pass loop: mmg3d_i feeds mmg3d_{i+1}.  R-003 fix.
        current_in = in_msh
        passes_log: list[dict] = []
        first_input_counts: dict | None = None
        last_output_counts: dict | None = None
        for pass_idx in range(n_passes):
            is_last = (pass_idx == n_passes - 1)
            # Final pass writes to a per-trial candidate location;
            # earlier passes write to intermediate medit-derived msh.
            if is_last:
                trial_out = (work_dir
                             / f"trial{trial_idx}_final.msh")
            else:
                trial_out = (work_dir
                             / f"trial{trial_idx}_p{pass_idx}.msh")
            log_suffix = f"t{trial_idx}_p{pass_idx}"
            ci, co = _run_one_mmg3d_pass(
                current_in, trial_out, work_dir, log_suffix,
                hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
                fault_tri_tag=fault_tri_tag, mode=mode,
                provenance_path=provenance_path,
                transform_path=transform_path,
                snap_m=snap_m, binary=binary,
                mode_flags=mode_flags)
            if first_input_counts is None:
                first_input_counts = ci
            last_output_counts = co
            passes_log.append({
                "pass_idx": pass_idx,
                "input_counts": ci,
                "output_counts": co,
            })
            current_in = trial_out
            print(f"  trial {trial_idx} pass {pass_idx}: "
                  f"{co['n_tets']} tets, "
                  f"{co['n_triangles']} tris",
                  file=sys.stderr)

        # Final γ_min and min element edge for this trial.
        gamma = _gamma_min_of_msh(current_in,
                                    fault_tri_tag=fault_tri_tag)
        min_edge = _min_edge_of_msh(current_in)
        trials_log.append({
            "trial_idx":   trial_idx,
            "final_msh":   str(current_in),
            "gamma_min":   gamma,
            "min_edge_m":  min_edge,
            "passes":      passes_log,
        })
        print(f"  trial {trial_idx} γ_min: {gamma:.6e}  "
              f"min_edge: {min_edge:.3f} m", file=sys.stderr)

        if best is None or gamma > best["gamma_min"]:
            best = {
                "trial_idx": trial_idx,
                "gamma_min": gamma,
                "min_edge_m": min_edge,
                "final_msh": current_in,
                "first_input_counts": first_input_counts,
                "last_output_counts": last_output_counts,
            }

    # Copy the best trial's final mesh to the user-requested out_msh.
    assert best is not None
    shutil.copy(best["final_msh"], out_msh)
    print(f"  best trial: {best['trial_idx']} "
          f"(γ_min = {best['gamma_min']:.6e}, "
          f"min_edge = {best['min_edge_m']:.3f} m)",
          file=sys.stderr)

    if not keep_intermediate:
        shutil.rmtree(work_dir, ignore_errors=True)

    return {
        "in_msh":       str(in_msh),
        "out_msh":      str(out_msh),
        "hmin": hmin, "hmax": hmax, "hgrad": hgrad, "hausd": hausd,
        "fault_tri_tag": fault_tri_tag,
        "mode":         mode,
        "binary":       binary,
        "n_passes":     n_passes,
        "n_trials":     n_trials,
        "best_trial_idx":   best["trial_idx"],
        "best_gamma_min":   best["gamma_min"],
        "best_min_edge_m":  best["min_edge_m"],
        "trials":           trials_log,
        "input_counts":  best["first_input_counts"],
        "output_counts": best["last_output_counts"],
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Post-HXT mmg3d sliver-removal pass.")
    p.add_argument("--in-msh", required=True, type=Path)
    p.add_argument("--out-msh", required=True, type=Path)
    p.add_argument("--hmin", type=float, default=100.0,
                   help="mmg3d -hmin (m).  Default 100.")
    p.add_argument("--hmax", type=float, default=25000.0,
                   help="mmg3d -hmax (m).  Default 25000 (= res_ff).")
    p.add_argument("--hgrad", type=float, default=1.3,
                   help="mmg3d -hgrad.  Default 1.3.")
    p.add_argument("--hausd", type=float, default=50.0,
                   help="mmg3d -hausd (m).  Default 50.")
    p.add_argument("--fault-tri-tag", type=int, default=100,
                   help="Triangle physical tag identifying fault "
                        "triangles to preserve as RequiredTriangles.  "
                        "Default 100 (matches safs.geo).")
    p.add_argument("--mode",
                   choices=("optim", "adapt", "optim_relax_fault"),
                   default="optim",
                   help="'optim' (default): mmg3d -optim -nosurf -opnbdy "
                        "— improve tet quality WITHOUT surface "
                        "modification or uniform re-sizing.  "
                        "'adapt': full mesh adaptation enforcing "
                        "hmin/hmax/hgrad uniformly (use only when the "
                        "input mesh sizing is wrong).  "
                        "'optim_relax_fault': allow mmg3d to modify "
                        "fault surface triangles for sliver removal, "
                        "while protecting cross-fault polyline + "
                        "free-surface trace + fault perimeter via "
                        "explicit RequiredEdges (requires "
                        "--provenance-json).")
    p.add_argument("--provenance-json", type=Path, default=None,
                   help="Path to fault_provenance.json (required for "
                        "mode='optim_relax_fault').  Used to identify "
                        "cross-fault polyline vertices.")
    p.add_argument("--transform-json", type=Path, default=None,
                   help="Path to transform.json (used by "
                        "mode='optim_relax_fault' to read "
                        "free_surface_clearance_m).")
    p.add_argument("--snap-m", type=float, default=0.1,
                   help="Snap precision for cross-fault vertex "
                        "matching (default 0.1 m).")
    p.add_argument("--n-passes", type=int, default=1,
                   help="Number of sequential mmg3d passes "
                        "(R-003 fix).  Each pass feeds the next; γ_min "
                        "typically improves 5-50× per pass for the "
                        "first 2-3 iterations.  Default 1.")
    p.add_argument("--n-trials", type=int, default=1,
                   help="Number of independent multi-pass trials "
                        "(R-006 fix).  mmg3d is non-deterministic "
                        "across runs; trials with different seeds "
                        "yield different γ_min by 5-55×.  The trial "
                        "with the highest γ_min is selected.  "
                        "Total mmg3d invocations = n_passes * "
                        "n_trials.  Default 1.")
    p.add_argument("--mmg3d-binary", default="mmg3d_O3")
    p.add_argument("--keep-intermediate", action="store_true")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = post_pass(
        in_msh=args.in_msh, out_msh=args.out_msh,
        hmin=args.hmin, hmax=args.hmax,
        hgrad=args.hgrad, hausd=args.hausd,
        fault_tri_tag=args.fault_tri_tag,
        binary=args.mmg3d_binary,
        keep_intermediate=args.keep_intermediate,
        mode=args.mode,
        provenance_path=args.provenance_json,
        transform_path=args.transform_json,
        snap_m=args.snap_m,
        n_passes=args.n_passes,
        n_trials=args.n_trials,
    )

    report_path = args.report_json or (
        args.out_msh.parent / "mmg3d_post_pass_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
