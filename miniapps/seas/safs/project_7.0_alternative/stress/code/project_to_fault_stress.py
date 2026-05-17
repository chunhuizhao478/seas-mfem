#!/usr/bin/env python3
"""
project_to_fault_stress.py — Project a regional bulk Cauchy stress
tensor σ⁰ onto the SAFS fault and bulk meshes.

Phase 1 of the pipeline (this file's current scope) supplies the mesh
I/O and per-cell geometry foundations. Phases 2–5 add the fault basis,
the stress rotation, and the VTU / JSON / schema-v1 HDF5 writers.

Conventions
-----------
- Coordinate frame: UTM Zone 11 N metres; (x = east, y = north, z = up).
- Sign convention split between the H&Z and SEAS conventions follows
  the contract in PLAN_onfaultstress.md §"Constraints":
    * The bulk Cauchy field stays in H&Z's compression-NEGATIVE
      convention throughout `bulk_stress_tensor_field` and the internal
      computation of `resolve_traction_per_cell`.
    * Compression-POSITIVE SEAS convention enters at two write
      boundaries only: (i) on the resolved σ_n scalar in
      `resolve_traction_per_cell`, and (ii) at the bulk VTU / schema-v1
      sidecar write sites.
- Fault-local frame uses the Tandem / SEAS convention `(s, d, n)` with
  `s = up × n`, `d = s × n` (down-dip), matching CLAUDE.md's
  "Fault-local tangent frame" rule. This differs from H&Z's
  `d_HZ = n × s` (up-dip).

Usage
-----
CLI (Phase 4 entry point — runs on the SAFS fault mesh):

    python project_to_fault_stress.py <INPUT_FAULT_VTU> <MESH_BASE>
        [--write-bulk] [--SHmax 113.0] [--Shmin 49.0] [--Sv 45.0]
        [--P_p 16.0] [--SHmax-az 23.0]
        [--depth-model {constant,lithostatic_sv}]
        [--SHmax-grad 0.0] [--Shmin-grad 0.0] [--Sv-grad 0.0]
        [--P_p-grad 0.0]
        [--hz-dump-file PATH]
        [--strike-hint-az 314.0] [--no-strike-hint]
        [--rake-sense right-lateral]
        [--out-dir DATA_PROJECTED_ONFAULTSTRESS]

    python project_to_fault_stress.py --batch [same flags]
    python project_to_fault_stress.py --print-info INPUT_FAULT_VTU

ParaView verification recipe (Phase 4 acceptance criterion):
1. Open the produced `<base>_fault_stress.vtu` in ParaView.
2. Apply the "Surface" representation; colour by `tau_magnitude_MPa`.
3. Add an arrow glyph filter coloured by `tau_strike_MPa` to inspect
   the right-lateral shear distribution.

Library use:

    from project_to_fault_stress import (
        load_fault_mesh, load_bulk_mesh,
        triangle_geometry, tet_geometry,
        extract_cells_by_physical,
        build_fault_basis, cell_to_node_average, basis_to_node,
        resolve_traction_per_cell, bulk_stress_tensor_field,
        pore_pressure_field, project_stress_onto_fault,
        write_fault_vtu, write_bulk_vtu, write_summary_json,
    )

Tests
-----
    cd project_7.0_alternative/stress/code && pytest -q test_project_to_fault_stress.py
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import meshio
import numpy as np


# ----------------------------------------------------------------------
# Sibling import of the H&Z regional-stress module
# ----------------------------------------------------------------------

# H&Z sibling lives in the same stress/code/ directory after the
# project_7.0_alternative reorg.
_HZ_DIR = Path(__file__).resolve().parent
if str(_HZ_DIR) not in sys.path:
    sys.path.insert(0, str(_HZ_DIR))
from hickman_and_zoback_regional_stress_projection import (  # noqa: E402
    ResolvedTraction,
    build_bulk_stress_tensor,
    fault_basis_vectors,
    resolve_traction,
)


# ----------------------------------------------------------------------
# Module-level constants
# ----------------------------------------------------------------------

EPS: float = 1.0e-9
UP_VECTOR: np.ndarray = np.array([0.0, 0.0, 1.0], dtype=np.float64)

# |n_z| above this means dip ≈ 0 and strike is degenerate (Phase 2).
NEAR_HORIZONTAL_NZ_TOL: float = 1.0 - 1.0e-6

DEFAULT_BULK_NAME: str = "rock"
DEFAULT_FAULT_NAME: str = "fault"

DEFAULT_OUT_DIR: Path = (
    Path(__file__).resolve().parent.parent / "results"
)

DEFAULT_MESH_GLOB: str = "safs_fault_box_nwcut_*m.msh"

DEFAULT_INPUT_BASES: tuple[str, ...] = (
    "safs_fault_box_nwcut_500m",
    "safs_fault_box_nwcut_1000m",
    "safs_fault_box_nwcut_2000m",
)

# CRS sanity-check threshold (X-coordinate magnitude in metres). UTM
# coordinates are typically ≥ 100 km; values below this likely mean a
# local or non-UTM frame and the caller is warned (not failed) per the
# Phase 1 §"Edge Cases" requirement.
_NON_UTM_X_MAX_THRESHOLD_M: float = 1.0e4


# ----------------------------------------------------------------------
# Dataclasses
# ----------------------------------------------------------------------


@dataclass
class FaultCellGeometry:
    """Per-triangle geometry for the fault surface.

    Notes
    -----
    ``n_degenerate`` conflates two distinct degeneracy modes that
    the plan (``PLAN_onfaultstress.md`` Phase 1 §"Edge Cases" and
    Phase 2 §"Edge Cases") describes as separate concerns:

    1. **Zero-area triangles** (mesh pathology): triangle_geometry
       emits NaN normal/centroid; these propagate to NaN basis rows
       in per_triangle_basis_raw and increment ``n_degenerate``.
    2. **Near-horizontal triangles** (geometric pathology): valid
       normal but |up × n| < tolerance, strike undefined; also
       counted in ``n_degenerate``.

    The dataclass provides only a single integer field. For the
    SAFS 2 km mesh both counts are zero, so the conflation has no
    operational impact. A future mesh exhibiting both kinds of
    pathology will need a richer diagnostic structure.
    """

    centroids: np.ndarray   # (N_tri, 3)
    normals: np.ndarray     # (N_tri, 3), unit, orientation-harmonised
    strikes: np.ndarray     # (N_tri, 3), unit; NaN for degenerate rows
    dips: np.ndarray        # (N_tri, 3), unit; NaN for degenerate rows
    areas: np.ndarray       # (N_tri,)
    n_degenerate: int       # see class docstring


@dataclass
class BulkCellGeometry:
    """Per-tet geometry for the bulk volume."""

    centroids: np.ndarray   # (N_tet, 3)
    volumes: np.ndarray     # (N_tet,)


# ----------------------------------------------------------------------
# Mesh I/O
# ----------------------------------------------------------------------


def load_fault_mesh(path: Path) -> meshio.Mesh:
    """Load a `*_fault.vtu` (written by `msh_to_vtu.py`).

    Returns the full `meshio.Mesh` with its point cloud and triangle
    cell blocks. Caller is responsible for selecting which physical
    group to use; see `extract_cells_by_physical`.

    Raises
    ------
    FileNotFoundError
        If the path does not exist (with the absolute path in the
        message — mirrors `nw_cut_strip.py:load_vertices` pattern).
    ValueError
        If the file contains no triangle cells.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"fault mesh not found: {p.resolve()}")

    mesh = meshio.read(str(p))
    if not any(cb.type == "triangle" for cb in mesh.cells):
        raise ValueError(
            f"fault mesh has no triangle cells: {p.resolve()} "
            f"(cell types present: {[cb.type for cb in mesh.cells]})"
        )
    _maybe_warn_non_utm(mesh.points, p)
    return mesh


def load_bulk_mesh(path: Path) -> meshio.Mesh:
    """Load a `*_bulk.vtu` (written by `msh_to_vtu.py`).

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the file contains no tetra cells.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"bulk mesh not found: {p.resolve()}")

    mesh = meshio.read(str(p))
    if not any(cb.type == "tetra" for cb in mesh.cells):
        raise ValueError(
            f"bulk mesh has no tetra cells: {p.resolve()} "
            f"(cell types present: {[cb.type for cb in mesh.cells]})"
        )
    _maybe_warn_non_utm(mesh.points, p)
    return mesh


def _maybe_warn_non_utm(
    points: np.ndarray,
    path: Path,
    threshold_m: float = _NON_UTM_X_MAX_THRESHOLD_M,
) -> None:
    """Emit a stderr warning if the mesh's X-coordinate magnitude is
    too small to be UTM metres. Does not raise.

    The ``threshold_m`` kwarg is exposed for callers who load meshes
    in non-default UTM zones / clipped local frames; the default is
    permissive enough to accept any UTM Zone 1–60 mesh
    (typical UTM X-coords are 150–850 km)."""
    if points.size == 0:
        return
    xmax = float(np.max(np.abs(points[:, 0])))
    if xmax < threshold_m:
        print(
            f"warning: mesh X-coordinate max magnitude {xmax:.2e} m is "
            f"below the UTM sanity threshold "
            f"({threshold_m:.0e} m); the file "
            f"{path.resolve()} may be in a non-UTM (local or rotated) "
            f"frame. The geometry math is linear-algebraic and still "
            f"works, but the H&Z `(east, north, up)` convention will "
            f"be mis-applied if the frame differs.",
            file=sys.stderr,
        )


# ----------------------------------------------------------------------
# Per-cell geometry
# ----------------------------------------------------------------------


def triangle_geometry(
    points: np.ndarray,
    tri_conn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute centroids, unit normals, and areas for triangles.

    Mirrors the algebra in `msh_to_vtu.py:43-58` (vol/area/normal via
    cross product) and `ts_to_stl.py:154-160` (per-triangle normal).

    Parameters
    ----------
    points : (N_pt, 3) float
        Mesh point cloud.
    tri_conn : (N_tri, 3) int
        Triangle vertex indices into `points`.

    Returns
    -------
    centroids : (N_tri, 3) float
        Mean of the three vertices per triangle.
    normals : (N_tri, 3) float
        Unit normal, computed via `(p1 - p0) × (p2 - p0)` then
        normalised. Degenerate (zero-area) triangles get a NaN normal
        and a NaN centroid (the centroid is technically well-defined,
        but we mark the whole row to make NaN propagation in Phase 2
        explicit).
    areas : (N_tri,) float
        `0.5 * |(p1 - p0) × (p2 - p0)|`. Zero for degenerate
        triangles (not NaN — caller can sum areas without filtering).
    """
    points = np.asarray(points, dtype=np.float64)
    tri_conn = np.asarray(tri_conn, dtype=np.int64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(
            f"points must be (N_pt, 3); got shape {points.shape}"
        )
    if tri_conn.ndim != 2 or tri_conn.shape[1] != 3:
        raise ValueError(
            f"tri_conn must be (N_tri, 3); got shape {tri_conn.shape}"
        )
    n_pt = points.shape[0]
    if tri_conn.size > 0:
        idx_min = int(tri_conn.min())
        idx_max = int(tri_conn.max())
        if idx_min < 0 or idx_max >= n_pt:
            raise ValueError(
                f"tri_conn indices out of bounds: min={idx_min}, "
                f"max={idx_max}, but points has only {n_pt} rows"
            )

    p = points[tri_conn]                          # (N_tri, 3, 3)
    centroids = p.mean(axis=1)                    # (N_tri, 3)
    e1 = p[:, 1] - p[:, 0]
    e2 = p[:, 2] - p[:, 0]
    cross = np.cross(e1, e2)                      # (N_tri, 3)
    twice_area = np.linalg.norm(cross, axis=1)
    areas = 0.5 * twice_area
    # Degenerate-triangle mask: zero area.
    degen = twice_area < EPS
    # Unit normal; NaN where degenerate.
    safe = np.where(degen[:, None], 1.0, twice_area[:, None])
    normals = cross / safe
    if degen.any():
        normals[degen] = np.nan
        centroids[degen] = np.nan
    return centroids, normals, areas


def tet_geometry(
    points: np.ndarray,
    tet_conn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute centroids and absolute volumes for tetrahedra.

    Mirrors `msh_to_vtu.py:tet_quality:40-45` for the volume formula.

    Parameters
    ----------
    points : (N_pt, 3) float
        Mesh point cloud.
    tet_conn : (N_tet, 4) int
        Tetrahedron vertex indices.

    Returns
    -------
    centroids : (N_tet, 3) float
        Mean of the four vertices per tetrahedron.
    volumes : (N_tet,) float
        `|det([p1-p0, p2-p0, p3-p0])| / 6`. Always non-negative
        (caller-friendly); use the signed scalar-triple-product
        directly if a chirality check is needed.

    Notes
    -----
    The Phase 1 plan's signature comment says "signed volumes (N,)"
    but its body says "use |V_signed| as the cell weight". This
    implementation follows the body rule (always-positive output).
    A caller that needs the signed scalar-triple-product (e.g. for
    a mesh chirality check) should compute it directly via
    ``np.einsum("ij,ij->i", a, np.cross(b, c))`` and divide by 6.
    """
    points = np.asarray(points, dtype=np.float64)
    tet_conn = np.asarray(tet_conn, dtype=np.int64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(
            f"points must be (N_pt, 3); got shape {points.shape}"
        )
    if tet_conn.ndim != 2 or tet_conn.shape[1] != 4:
        raise ValueError(
            f"tet_conn must be (N_tet, 4); got shape {tet_conn.shape}"
        )
    n_pt = points.shape[0]
    if tet_conn.size > 0:
        idx_min = int(tet_conn.min())
        idx_max = int(tet_conn.max())
        if idx_min < 0 or idx_max >= n_pt:
            raise ValueError(
                f"tet_conn indices out of bounds: min={idx_min}, "
                f"max={idx_max}, but points has only {n_pt} rows"
            )

    p = points[tet_conn]                          # (N_tet, 4, 3)
    centroids = p.mean(axis=1)                    # (N_tet, 3)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    triple = np.einsum("ij,ij->i", a, np.cross(b, c))
    volumes = np.abs(triple) / 6.0
    return centroids, volumes


# ----------------------------------------------------------------------
# Physical-group extraction
# ----------------------------------------------------------------------


def extract_cells_by_physical(
    mesh: meshio.Mesh,
    cell_type: str,
    phys_name: str,
) -> np.ndarray:
    """Return the connectivity (rows) of cells belonging to a named
    gmsh physical group.

    Factored out of `msh_to_vtu.py:120-145`. Reads `mesh.field_data`
    to map `phys_name → tag`; falls back to the unique-value heuristic
    if the VTU lost `field_data` on write (meshio does not always carry
    it through).

    Parameters
    ----------
    mesh : meshio.Mesh
        Mesh loaded from a `_fault.vtu` or `_bulk.vtu`.
    cell_type : str
        Either ``"triangle"`` or ``"tetra"``.
    phys_name : str
        Physical-group name (e.g. ``"fault"`` or ``"rock"``).

    Returns
    -------
    conn : (N, D) int
        Concatenated connectivity rows from every matching cell block
        (D = 3 for triangles, 4 for tetra).

    Raises
    ------
    ValueError
        If no matching cell block is found.
    """
    if cell_type not in ("triangle", "tetra"):
        raise ValueError(
            f"cell_type must be 'triangle' or 'tetra'; got {cell_type!r}"
        )

    field_data = getattr(mesh, "field_data", None) or {}
    name_to_tag = {
        name: tag for name, (tag, _dim) in field_data.items()
    }
    phys_tag: Optional[int] = name_to_tag.get(phys_name)
    # R-705: distinguish "field_data missing entirely" (silent
    # unique-value fallback is OK; documented in Phase 1
    # §"Implementation note") from "field_data present but
    # missing the requested name" (probably a typo or
    # multi-fault setup — surface the error instead of guessing).
    if field_data and phys_tag is None:
        raise ValueError(
            f"physical group {phys_name!r} not found in mesh "
            f"field_data; known groups: {sorted(name_to_tag)!r}"
        )

    phys_arrays = mesh.cell_data.get("gmsh:physical") if mesh.cell_data else None

    collected: list[np.ndarray] = []
    for i, cb in enumerate(mesh.cells):
        if cb.type != cell_type:
            continue
        if phys_arrays is None or i >= len(phys_arrays):
            # No physical tag available — fall back to taking the
            # whole block (the field_data block may have been dropped
            # on VTU write). Only safe when there is exactly one
            # matching block; if there are several we cannot
            # disambiguate without phys tags.
            collected.append(np.asarray(cb.data, dtype=np.int64))
            continue
        phys_arr = np.asarray(phys_arrays[i])
        if phys_tag is not None:
            mask = phys_arr == phys_tag
            if mask.any():
                collected.append(
                    np.asarray(cb.data[mask], dtype=np.int64)
                )
        else:
            # field_data didn't include phys_name — fall back to the
            # unique-value heuristic per the Phase 1 §"Implementation
            # note".
            unique = np.unique(phys_arr)
            if len(unique) == 1:
                collected.append(np.asarray(cb.data, dtype=np.int64))

    if not collected:
        available = {
            name: tag for name, (tag, _) in field_data.items()
        }
        raise ValueError(
            f"no {cell_type!r} cells found for physical group "
            f"{phys_name!r}; available physical groups in mesh: "
            f"{available!r}"
        )

    return np.concatenate(collected, axis=0)


# ----------------------------------------------------------------------
# Phase 2 — Per-triangle fault basis
# ----------------------------------------------------------------------

# Tolerance on |up × n| below which the strike direction is considered
# degenerate (triangle near-horizontal). Derived from
# NEAR_HORIZONTAL_NZ_TOL via |up × n| = sqrt(1 - n_z²) for unit
# vectors and up = ẑ. NEAR_HORIZONTAL_NZ_TOL is the *lower bound*
# on |n_z| at which we flag a triangle as horizontal — by default
# NEAR_HORIZONTAL_NZ_TOL = 1 - 1e-6 = 0.999999. The corresponding
# strike-norm threshold is sqrt(1 - NEAR_HORIZONTAL_NZ_TOL²) ≈
# sqrt(2 · 1e-6) ≈ 1.4e-3. Per the Phase 2 §1 spec.
_DEGEN_STRIKE_NORM_TOL: float = float(
    np.sqrt(1.0 - NEAR_HORIZONTAL_NZ_TOL ** 2)
)

# Maximum allowed ratio (smallest / largest PCA eigenvalue) for the
# PCA-fallback fault-normal estimation. Above this the mesh's
# triangulated centroids don't approximate a plane and the smallest
# eigenvector is unreliable; we warn but continue.
_PCA_PLANARITY_RATIO_MAX: float = 1.0e-2


def per_triangle_basis_raw(
    normals: np.ndarray,
    up: np.ndarray = UP_VECTOR,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a per-triangle right-handed orthonormal basis (s, d, n).

    Follows the Tandem ``facetBasis`` algorithm
    (`tandem/src/geometry/Curvilinear.cpp:281-292`) with the
    **down-dip** convention `d = s × n` (CLAUDE.md "Fault-local tangent
    frame" rule `can_t1 = (0, 0, -1)`).

    For each triangle ``i``:

    - ``s_i = up × n_i``; if ``|s_i| < _DEGEN_STRIKE_NORM_TOL`` then
      the triangle is near-horizontal and the strike direction is
      undefined → mark `degen[i] = True` and set strike / dip rows to
      NaN.
    - Otherwise ``strikes[i] = s_i / |s_i|`` (unit, along strike) and
      ``dips[i] = strikes[i] × n_i`` (unit, down-dip; this differs
      from H&Z's `n × s` up-dip convention by a sign).

    NaN normals (from degenerate zero-area triangles in Phase 1) flow
    through to NaN basis rows; they are counted as degenerate.

    Parameters
    ----------
    normals : (N, 3) float
        Unit normals per triangle.
    up : (3,) float, default UP_VECTOR
        Global up vector (typically `ẑ`).

    Returns
    -------
    strikes : (N, 3) float
        Unit along-strike vectors; NaN rows where degenerate.
    dips : (N, 3) float
        Unit down-dip vectors; NaN rows where degenerate.
    degen : (N,) bool
        True where the strike direction was undefined.
    """
    normals = np.asarray(normals, dtype=np.float64)
    up = np.asarray(up, dtype=np.float64)
    if normals.ndim != 2 or normals.shape[1] != 3:
        raise ValueError(
            f"normals must be (N, 3); got shape {normals.shape}"
        )
    if up.shape != (3,):
        raise ValueError(f"up must be (3,); got shape {up.shape}")

    N = normals.shape[0]
    # Strike = up × n
    s_raw = np.cross(np.broadcast_to(up, (N, 3)), normals)
    s_norm = np.linalg.norm(s_raw, axis=1)
    # Triangles with near-zero |up × n| are horizontal (strike
    # undefined); triangles whose normal itself is NaN propagate to
    # NaN s_norm and fall into the same degenerate bucket.
    degen = ~(s_norm > _DEGEN_STRIKE_NORM_TOL)  # NaN-safe: NaN > x is False
    strikes = np.full_like(normals, np.nan)
    dips = np.full_like(normals, np.nan)
    good = ~degen
    if good.any():
        strikes[good] = s_raw[good] / s_norm[good, None]
        # Down-dip: d = s × n (Tandem / SEAS convention).
        dips[good] = np.cross(strikes[good], normals[good])
    return strikes, dips, degen


def harmonise_normal_orientation(
    normals: np.ndarray,
    rake_sense: str,
    trace_centroids: np.ndarray,
    fault_strike_azimuth_hint_deg: Optional[float],
) -> np.ndarray:
    """Flip per-triangle normals to a consistent global orientation.

    For a right-lateral SAF with NW strike (azimuth 314°), the H&Z
    convention places ``n`` in the SW half-space everywhere (n_hat
    "to the left of strike in map view"). A raw mesh whose triangles
    were CCW-oriented gives a mix of SW- and NE-pointing normals; this
    function detects the global outward direction and flips
    individual normals to agree.

    Parameters
    ----------
    normals : (N, 3) float
        Unit normals (may contain NaN rows; passed through unchanged).
    rake_sense : str
        ``"right-lateral"`` or ``"left-lateral"``.
    trace_centroids : (N, 3) float
        Triangle centroids; used for the PCA fallback. The current
        implementation uses only the rake-hint path; the centroids
        are accepted for API symmetry with the PCA fallback.
    fault_strike_azimuth_hint_deg : float or None
        - If a float: use H&Z's ``fault_basis_vectors`` at this strike
          (with ``dip_deg=90``) to derive the reference outward
          direction.
        - If None: fall back to PCA on ``trace_centroids`` to derive
          the plane normal, then sign-match against the majority of
          ``normals``.

    Returns
    -------
    harmonised : (N, 3) float
        Normals flipped so that ``dot(n_i, n_global) > 0`` for every
        non-NaN row. Idempotent: calling twice yields the same result.

    Raises
    ------
    ValueError
        For an unrecognised ``rake_sense``.
    """
    normals = np.asarray(normals, dtype=np.float64)
    trace_centroids = np.asarray(trace_centroids, dtype=np.float64)
    # Up-front shape validation: callers passing 1-D arrays should
    # see a clear ValueError, not numpy.AxisError later (R-403).
    if normals.ndim != 2 or normals.shape[1] != 3:
        raise ValueError(
            f"normals must be (N, 3); got shape {normals.shape}"
        )
    if trace_centroids.ndim != 2 or trace_centroids.shape[1] != 3:
        raise ValueError(
            f"trace_centroids must be (N, 3); got shape "
            f"{trace_centroids.shape}"
        )
    if rake_sense not in ("right-lateral", "left-lateral"):
        raise ValueError(
            f"rake_sense must be 'right-lateral' or 'left-lateral'; "
            f"got {rake_sense!r}"
        )
    # trace_centroids is only consulted when the hint is None
    # (PCA fallback); validate shape MATCH only in that branch so
    # the hint path is robust to caller-supplied dummies (R-305).
    if (fault_strike_azimuth_hint_deg is None
            and normals.shape != trace_centroids.shape):
        raise ValueError(
            f"PCA fallback requires normals shape {normals.shape} "
            f"== trace_centroids shape {trace_centroids.shape}"
        )

    n_global = _global_outward_normal(
        normals, trace_centroids, rake_sense,
        fault_strike_azimuth_hint_deg,
    )

    # Flip each non-NaN row so dot(n_i, n_global) > 0.
    # Rows already aligned (dot > 0) keep their sign; rows opposite
    # (dot < 0) flip; ties (dot == 0) keep their sign (the
    # PCA-anchored reference makes exact zero vanishingly unlikely).
    out = normals.copy()
    finite = ~np.isnan(normals).any(axis=1)
    dot = np.einsum("ij,j->i", out, n_global)
    flip = finite & (dot < 0)
    if flip.any():
        out[flip] = -out[flip]
    return out


def _global_outward_normal(
    normals: np.ndarray,
    trace_centroids: np.ndarray,
    rake_sense: str,
    fault_strike_azimuth_hint_deg: Optional[float],
) -> np.ndarray:
    """Compute the global outward fault normal for orientation harmonisation."""
    if fault_strike_azimuth_hint_deg is not None:
        _s_hat, _d_hat, n_hat_global = fault_basis_vectors(
            float(fault_strike_azimuth_hint_deg),
            dip_deg=90.0,
            rake_sense=rake_sense,
        )
        return np.asarray(n_hat_global, dtype=np.float64)

    # PCA fallback on triangle centroids.
    finite = ~np.isnan(trace_centroids).any(axis=1)
    if not finite.any():
        raise ValueError(
            "PCA fallback for global fault normal requires at least one "
            "finite centroid row; all rows are NaN"
        )
    pts = trace_centroids[finite]
    centred = pts - pts.mean(axis=0)
    # eigvalsh / eigh because the covariance is symmetric.
    cov = centred.T @ centred / max(1, pts.shape[0] - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # eigvalsh returns ascending order; smallest eigenvector is column 0.
    n_pca = eigvecs[:, 0]
    # Planarity check: smallest / largest eigenvalue ratio.
    if eigvals[-1] > 0:
        ratio = eigvals[0] / eigvals[-1]
        if ratio > _PCA_PLANARITY_RATIO_MAX:
            print(
                f"warning: PCA-fallback fault-normal estimation found "
                f"poor planarity (eigenvalue ratio {ratio:.3e} > "
                f"threshold {_PCA_PLANARITY_RATIO_MAX:.0e}); the "
                f"smallest eigenvector may not be a meaningful global "
                f"fault normal. Consider passing "
                f"`fault_strike_azimuth_hint_deg` explicitly.",
                file=sys.stderr,
            )
    # Sign-match against the majority of (finite) input normals.
    n_in_finite = normals[finite]
    finite_in = ~np.isnan(n_in_finite).any(axis=1)
    if finite_in.any():
        dots = n_in_finite[finite_in] @ n_pca
        if (dots < 0).sum() > (dots > 0).sum():
            n_pca = -n_pca
    return n_pca.astype(np.float64)


def build_fault_basis(
    fault_mesh: meshio.Mesh,
    fault_phys_name: str = DEFAULT_FAULT_NAME,
    rake_sense: str = "right-lateral",
    fault_strike_azimuth_hint_deg: Optional[float] = 314.0,
    up: np.ndarray = UP_VECTOR,
) -> FaultCellGeometry:
    """Build the full per-triangle fault basis from a fault VTU mesh.

    Wires Phase 1 geometry (`triangle_geometry`,
    `extract_cells_by_physical`) + Phase 2 basis helpers
    (`per_triangle_basis_raw`, `harmonise_normal_orientation`).

    Parameters
    ----------
    fault_mesh : meshio.Mesh
        Loaded fault mesh (use ``load_fault_mesh``).
    fault_phys_name : str
        Physical-group name; defaults to ``DEFAULT_FAULT_NAME``.
    rake_sense : str
        ``"right-lateral"`` or ``"left-lateral"``.
    fault_strike_azimuth_hint_deg : float or None
        Optional global strike hint for orientation harmonisation;
        default 314° (SAF N46W). Pass ``None`` to use the PCA
        fallback.
    up : (3,) float
        Up vector for the basis algebra (default ẑ).

    Returns
    -------
    FaultCellGeometry
        Centroids, harmonised normals, strikes, dips, areas, and
        ``n_degenerate`` count.

    Notes
    -----
    Multi-fault future case: this function operates on a single
    physical group at a time. A caller with multiple fault groups
    should loop over them. Do not silently merge groups here.
    """
    tri_conn = extract_cells_by_physical(
        fault_mesh, "triangle", fault_phys_name
    )
    centroids, normals_raw, areas = triangle_geometry(
        fault_mesh.points, tri_conn
    )
    normals = harmonise_normal_orientation(
        normals_raw, rake_sense, centroids,
        fault_strike_azimuth_hint_deg,
    )
    strikes, dips, degen = per_triangle_basis_raw(normals, up)
    return FaultCellGeometry(
        centroids=centroids,
        normals=normals,
        strikes=strikes,
        dips=dips,
        areas=areas,
        n_degenerate=int(degen.sum()),
    )


# ----------------------------------------------------------------------
# Phase 2 — Cell→node averaging helpers
# ----------------------------------------------------------------------


def cell_to_node_average(
    points: np.ndarray,
    tri_conn: np.ndarray,
    cell_values: np.ndarray,
    areas: np.ndarray,
) -> np.ndarray:
    """Area-weighted cell→vertex average for a per-triangle field.

    For each vertex ``v`` incident on triangles ``T_v``,

        f_node[v] = (Σ_{t ∈ T_v} A_t · f_cell[t]) / (Σ_{t ∈ T_v} A_t)

    Works uniformly for scalar (``cell_values.shape == (N_tri,)``)
    and vector (``cell_values.shape == (N_tri, K)``) inputs.

    Parameters
    ----------
    points : (N_pt, 3) float
        Mesh point cloud.
    tri_conn : (N_tri, 3) int
        Triangle vertex indices.
    cell_values : (N_tri, ...) float
        Per-cell values. ``cell_values.shape[0]`` must equal N_tri;
        ``cell_values.shape[1:]`` is the "trailing" broadcast partner
        and is preserved on output.
    areas : (N_tri,) float
        Per-triangle areas (used as weights).

    Returns
    -------
    node_values : (N_pt, ...) float
        Area-weighted averaged values at each vertex; ``np.nan`` rows
        for vertices not referenced by any triangle.

    Raises
    ------
    ValueError
        If shapes are inconsistent.
    """
    points = np.asarray(points, dtype=np.float64)
    tri_conn = np.asarray(tri_conn, dtype=np.int64)
    cell_values = np.asarray(cell_values, dtype=np.float64)
    areas = np.asarray(areas, dtype=np.float64)

    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(
            f"points must be (N_pt, 3); got shape {points.shape}"
        )
    if tri_conn.ndim != 2 or tri_conn.shape[1] != 3:
        raise ValueError(
            f"tri_conn must be (N_tri, 3); got shape {tri_conn.shape}"
        )
    N_tri = tri_conn.shape[0]
    if cell_values.shape[0] != N_tri:
        raise ValueError(
            f"cell_values.shape[0]={cell_values.shape[0]} != "
            f"N_tri={N_tri}"
        )
    if areas.shape != (N_tri,):
        raise ValueError(
            f"areas shape {areas.shape} != ({N_tri},)"
        )

    N_pt = points.shape[0]
    K = cell_values.shape[1:]   # trailing shape

    # Filter out cells whose values contain NaN: such cells must NOT
    # pollute their neighbouring vertices' averages.  A "good" cell
    # is one whose entire trailing-shape slice is finite AND whose
    # area is finite.  (R-302: previously, a single NaN cell would
    # poison every vertex it touched.)
    if len(K) == 0:
        cell_finite = np.isfinite(cell_values)
    else:
        cell_finite = np.isfinite(cell_values).reshape(
            N_tri, -1
        ).all(axis=1)
    cell_finite &= np.isfinite(areas)

    # Zero-out the area weight on bad cells so they contribute neither
    # to the area sum nor to the accumulator.
    good_areas = np.where(cell_finite, areas, 0.0)

    # Per-vertex area sum: w[v] = Σ_{t ∈ T_v, good} A_t
    w = np.zeros(N_pt, dtype=np.float64)
    np.add.at(w, tri_conn.flatten(), np.repeat(good_areas, 3))

    # Zero-out the cell values on bad cells too — combined with
    # good_areas = 0 this is doubly defensive (the area weight alone
    # already excludes them, but zeroing values guards against NaN
    # leaking through `0 * NaN = NaN`).
    cv_safe = np.where(
        cell_finite.reshape((N_tri,) + (1,) * len(K)),
        cell_values,
        0.0,
    )

    # Weighted accumulation:
    # weighted[3*t + j] = good_areas[t] · cv_safe[t]  for j = 0, 1, 2.
    # The broadcast handles K = () (scalar) and any higher trailing K
    # uniformly.
    areas_b = good_areas.reshape(
        (-1,) + (1,) * (1 + len(K))
    )                                                       # (N_tri, 1, *K)
    cv_b = cv_safe.reshape((N_tri, 1) + K)                  # (N_tri, 1, *K)
    weighted = (areas_b * cv_b).repeat(3, axis=1)           # (N_tri, 3, *K)
    weighted = weighted.reshape((3 * N_tri,) + K)           # (3*N_tri, *K)

    acc = np.zeros((N_pt,) + K, dtype=np.float64)
    np.add.at(acc, tri_conn.flatten(), weighted)

    # Divide by per-vertex area; orphan vertices (w == 0) → NaN.
    w_b = w.reshape((-1,) + (1,) * len(K))
    out = np.divide(
        acc, w_b,
        out=np.full_like(acc, np.nan),
        where=w_b > 0,
    )
    return out


def basis_to_node(
    points: np.ndarray,
    tri_conn: np.ndarray,
    strikes_cell: np.ndarray,
    dips_cell: np.ndarray,
    normals_cell: np.ndarray,
    areas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project the per-triangle (s, d, n) basis onto per-vertex slots.

    Used by the fault-VTU writer to emit a continuous (point-data)
    basis field. Re-orthonormalises at each vertex via Gram–Schmidt:

    1. ``n_node = normalize(Σ A_t · n_cell)``
    2. ``s_avg  = Σ A_t · s_cell`` (raw vector sum)
    3. ``s_proj = s_avg − (s_avg · n_node) · n_node``  (Gram–Schmidt)
    4. ``s_node = s_proj / |s_proj|``
    5. ``d_node = s_node × n_node``  (Tandem down-dip; already unit
       by construction; matches `per_triangle_basis_raw`'s
       ``d = s × n`` per CLAUDE.md "Fault-local tangent frame" rule
       ``can_t1 = (0, 0, -1)``. NOT ``n × s`` — that is H&Z up-dip
       and would invert the dip sign relative to the cell-data dip
       (round-2 R-301).)

    Vertices with degenerate Gram–Schmidt (`|s_proj| < EPS`) get NaN
    rows in all three returned arrays.

    Parameters
    ----------
    points : (N_pt, 3) float
    tri_conn : (N_tri, 3) int
    strikes_cell, dips_cell, normals_cell : (N_tri, 3) float
        Per-cell basis vectors (output of `per_triangle_basis_raw`
        after orientation harmonisation).
    areas : (N_tri,) float

    Returns
    -------
    strikes_node, dips_node, normals_node : (N_pt, 3) float
        Per-vertex orthonormal triple. NaN rows where degenerate.
    """
    strikes_cell = np.asarray(strikes_cell, dtype=np.float64)
    dips_cell = np.asarray(dips_cell, dtype=np.float64)
    normals_cell = np.asarray(normals_cell, dtype=np.float64)
    # Argument-shape validation deferred to cell_to_node_average,
    # which will raise on inconsistent inputs.
    _ = dips_cell  # not used directly; passing through cell-data
                   # averaging is unnecessary because dip is
                   # recomputed via s × n after re-orthonormalisation
                   # (Tandem down-dip convention; round-2 R-301).

    # Step 1: average normals, then renormalise.
    n_avg = cell_to_node_average(points, tri_conn, normals_cell, areas)
    n_norm = np.linalg.norm(n_avg, axis=1)
    n_node = np.full_like(n_avg, np.nan)
    safe_n = n_norm > EPS
    n_node[safe_n] = n_avg[safe_n] / n_norm[safe_n, None]

    # Step 2: average strikes.
    s_avg = cell_to_node_average(points, tri_conn, strikes_cell, areas)

    # Step 3-4: Gram–Schmidt onto plane ⊥ n_node, then normalise.
    s_node = np.full_like(s_avg, np.nan)
    d_node = np.full_like(s_avg, np.nan)

    # Project s_avg onto plane perpendicular to n_node (only where n
    # is safe). Use einsum to compute dot per row.
    if safe_n.any():
        s_avg_sn = s_avg[safe_n]
        n_node_sn = n_node[safe_n]
        dot_sn = np.einsum("ij,ij->i", s_avg_sn, n_node_sn)
        s_proj = s_avg_sn - dot_sn[:, None] * n_node_sn
        s_proj_norm = np.linalg.norm(s_proj, axis=1)
        safe_s_local = s_proj_norm > EPS

        # Materialise per-row outputs back into the full arrays.
        sn_indices = np.nonzero(safe_n)[0]
        good_indices = sn_indices[safe_s_local]
        s_node[good_indices] = (
            s_proj[safe_s_local] / s_proj_norm[safe_s_local, None]
        )
        # Step 5: down-dip d = s × n (Tandem / SEAS convention,
        # matches `per_triangle_basis_raw`'s `d = s × n` so that
        # cell-data and point-data dip vectors carry the same sign).
        # NB the plan is internally inconsistent: Phase 2 §5
        # line 675 says `cross(strikes_node, normals_node)` =
        # `s × n` (Tandem down-dip — correct, matches this code),
        # while the Phase 4 §1 continuous-field math summary at
        # plan line 1057 writes `d_node = n_node × s_node` =
        # `n × s` (H&Z up-dip — wrong).  This code follows
        # Phase 2 §5; plan line 1057 has been updated to match
        # (R-704).  Tandem is canonical per CLAUDE.md
        # "Fault-local tangent frame" rule (R-301).
        d_node[good_indices] = np.cross(
            s_node[good_indices], n_node[good_indices]
        )

    # Vertices where s_proj is degenerate must have ALL three rows
    # NaN — n_node row alone is meaningless without a usable strike.
    bad = np.isnan(s_node).any(axis=1)
    n_node[bad] = np.nan
    d_node[bad] = np.nan

    return s_node, d_node, n_node


# ----------------------------------------------------------------------
# Phase 3 — Stress projection (constant + depth-dependent)
# ----------------------------------------------------------------------

# Stress-tensor symmetry tolerance (MPa). The H&Z rotation should
# produce machine-precision-symmetric tensors; we re-check after the
# one-shot sign flip as a single source-of-truth defence.
_SIGMA_SYMMETRY_TOL_MPA: float = 1.0e-9


def _assert_symmetric(sigma: np.ndarray, tol: float) -> None:
    """Raise ValueError if a (3, 3) tensor is not symmetric within tol."""
    asym = float(np.max(np.abs(sigma - sigma.T)))
    if asym > tol:
        raise ValueError(
            f"sigma is not symmetric; max off-diagonal asymmetry "
            f"{asym:.3e} > tol {tol:.3e}"
        )


def resolve_traction_per_cell(
    sigma_global: np.ndarray,
    strikes: np.ndarray,
    dips: np.ndarray,
    normals: np.ndarray,
    P_p_per_cell=0.0,
) -> dict:
    """Rotate a Cauchy stress tensor onto the per-triangle Tandem
    fault basis and emit the SEAS-convention pre-stress components.

    **`sigma_global` is in compression-POSITIVE SEAS convention here**
    (σ_seas = −σ_HZ). The single-site sign flip from H&Z
    continuum-mechanics to SEAS happens at the source in
    `bulk_stress_tensor_field` (Phase 3 §2); this resolver is a pure
    rotation onto Phase 2's Tandem basis with **no further sign
    flips**.

    Under the unified σ_seas convention with Tandem basis
    (``s = up × n``, ``d = s × n``; locked in by Phase 2 R-301/R-102),
    all three emitted scalar components carry the SEAS-internal sign:

    - ``sigma_n_total > 0`` ↔ compression
    - ``tau_strike > 0`` ↔ right-lateral
    - ``tau_dip > 0`` ↔ reverse

    Parameters
    ----------
    sigma_global : (3, 3) or (N_tri, 3, 3) ndarray
        σ_seas in compression-positive SEAS convention (MPa).
        Constant-field case: pass a single (3, 3) tensor; the
        function broadcasts to every triangle. Depth-dependent case:
        pass a (N_tri, 3, 3) stack.
    strikes, dips, normals : (N_tri, 3) ndarray
        Per-triangle Tandem-basis unit vectors (output of
        ``build_fault_basis``).
    P_p_per_cell : float or (N_tri,) ndarray, default 0.0
        Pore pressure in MPa (positive value). Subtracted from
        ``sigma_n_total`` to give ``sigma_n_eff``.

    Returns
    -------
    dict[str, np.ndarray]
        Keys (each value is (N_tri,) except ``traction_vec`` which is
        (N_tri, 3)):

        - ``sigma_n_total`` — n · σ_seas · n (MPa, compression positive)
        - ``sigma_n_eff``   — sigma_n_total − P_p (MPa, compression positive)
        - ``tau_strike``    — s · σ_seas · n (MPa, right-lateral positive)
        - ``tau_dip``       — d · σ_seas · n (MPa, reverse positive)
        - ``tau_magnitude`` — sqrt(τ_s² + τ_d²) (MPa, ≥ 0)
        - ``rake_deg``      — atan2(τ_d, τ_s) (degrees)
        - ``mu_apparent``   — tau_magnitude / sigma_n_eff; NaN where
          sigma_n_eff ≤ 0 (effective tension — frictional Mohr-Coulomb
          model not applicable; logged in the JSON summary's
          ``stats.mu_apparent.nan_count``)
        - ``traction_vec``  — t = σ_seas · n (N_tri, 3); SEAS
          convention so n·t = sigma_n_total, s·t = tau_strike,
          d·t = tau_dip.

    Raises
    ------
    ValueError
        If shapes are inconsistent or σ_global is not symmetric
        within ``_SIGMA_SYMMETRY_TOL_MPA``.
    """
    sigma_global = np.asarray(sigma_global, dtype=np.float64)
    strikes = np.asarray(strikes, dtype=np.float64)
    dips = np.asarray(dips, dtype=np.float64)
    normals = np.asarray(normals, dtype=np.float64)

    if strikes.ndim != 2 or strikes.shape[1] != 3:
        raise ValueError(
            f"strikes must be (N_tri, 3); got shape {strikes.shape}"
        )
    if dips.shape != strikes.shape:
        raise ValueError(
            f"dips shape {dips.shape} != strikes shape {strikes.shape}"
        )
    if normals.shape != strikes.shape:
        raise ValueError(
            f"normals shape {normals.shape} != strikes shape "
            f"{strikes.shape}"
        )

    N_tri = strikes.shape[0]

    if sigma_global.shape == (3, 3):
        _assert_symmetric(sigma_global, _SIGMA_SYMMETRY_TOL_MPA)
        # t[k, i] = Σ_j σ[i, j] n[k, j]
        t = np.einsum("ij,kj->ki", sigma_global, normals)   # (N, 3)
    elif sigma_global.shape == (N_tri, 3, 3):
        # Per-cell symmetry check (vectorised). Skipped when N_tri=0
        # because np.max on an empty reduction raises an opaque
        # numpy error (R-601).
        if N_tri > 0:
            asym = float(np.max(
                np.abs(sigma_global - np.swapaxes(sigma_global, -1, -2))
            ))
            if asym > _SIGMA_SYMMETRY_TOL_MPA:
                raise ValueError(
                    f"sigma_global has non-symmetric rows; max "
                    f"off-diagonal asymmetry {asym:.3e} > "
                    f"{_SIGMA_SYMMETRY_TOL_MPA:.3e}"
                )
        # t[k, i] = Σ_j σ[k, i, j] n[k, j]
        t = np.einsum("kij,kj->ki", sigma_global, normals)
    else:
        raise ValueError(
            f"sigma_global must be (3, 3) or ({N_tri}, 3, 3); "
            f"got shape {sigma_global.shape}"
        )

    # Pure rotation onto Tandem basis — σ_seas already compression
    # positive. No further sign manipulation.
    sigma_n_total = np.einsum("ki,ki->k", normals, t)
    tau_strike = np.einsum("ki,ki->k", strikes, t)
    tau_dip = np.einsum("ki,ki->k", dips, t)
    tau_magnitude = np.sqrt(tau_strike ** 2 + tau_dip ** 2)
    rake_deg = np.degrees(np.arctan2(tau_dip, tau_strike))

    # Pore pressure: scalar or per-cell broadcast.
    P_p_arr = np.asarray(P_p_per_cell, dtype=np.float64)
    if P_p_arr.ndim == 0:
        pass  # scalar broadcast
    elif P_p_arr.shape != (N_tri,):
        raise ValueError(
            f"P_p_per_cell must be scalar or (N_tri={N_tri},); "
            f"got shape {P_p_arr.shape}"
        )
    sigma_n_eff = sigma_n_total - P_p_arr

    # μ_apparent: NaN where σ_n_eff ≤ 0 (effective tension,
    # Mohr-Coulomb model not applicable per Phase 3 §"Edge Cases").
    mu_apparent = np.full(N_tri, np.nan, dtype=np.float64)
    safe = sigma_n_eff > 0.0
    if safe.any():
        mu_apparent[safe] = tau_magnitude[safe] / sigma_n_eff[safe]

    return {
        "sigma_n_total": sigma_n_total,
        "sigma_n_eff": sigma_n_eff,
        "tau_strike": tau_strike,
        "tau_dip": tau_dip,
        "tau_magnitude": tau_magnitude,
        "rake_deg": rake_deg,
        "mu_apparent": mu_apparent,
        "traction_vec": t,
    }


def bulk_stress_tensor_field(
    centroids_z: np.ndarray,
    *,
    SHmax_az_deg: float,
    depth_model: str = "constant",
    SHmax_top: float,
    Shmin_top: float,
    Sv_top: float,
    SHmax_grad: float = 0.0,
    Shmin_grad: float = 0.0,
    Sv_grad: float = 0.0,
) -> np.ndarray:
    """Build the bulk σ⁰ field in compression-POSITIVE SEAS convention.

    **This is the SINGLE sign-flip site for the entire pipeline**
    (R-501/R-502). The function calls H&Z's
    ``build_bulk_stress_tensor`` (which returns compression-negative
    tensors per continuum-mechanics convention) and multiplies by
    −1 once at the source, returning σ_seas. Every downstream
    consumer (resolver, fault VTU, bulk VTU, sidecar, C++ projector,
    verifier) sees compression-positive MPa with no further sign
    flips.

    Parameters
    ----------
    centroids_z : (N,) ndarray
        z-coordinates of the cells where σ⁰ should be evaluated
        (typically ≤ 0 in the subsurface; positive values are
        clamped to depth 0 in the lithostatic branch).
    SHmax_az_deg : float
        Geological azimuth of SHmax (clockwise from north, degrees).
    depth_model : {"constant", "lithostatic_sv"}
        - ``"constant"`` (default): every row equals
          ``-build_bulk_stress_tensor(SHmax_top, Shmin_top, Sv_top,
          SHmax_az_deg)``.
        - ``"lithostatic_sv"``: principal magnitudes scale linearly
          with depth ``d = max(0, -z)``:
          ``SHmax(d) = SHmax_top + SHmax_grad * d`` (and analogously
          for Shmin / Sv).
    SHmax_top, Shmin_top, Sv_top : float
        Principal-stress magnitudes (POSITIVE values, MPa) at z = 0.
    SHmax_grad, Shmin_grad, Sv_grad : float, default 0.0
        Linear depth gradients (MPa/m). Default 0.0 → constant.
        **No hardcoded ρ̄·g gradient**; depth-dependence is opt-in.

    Returns
    -------
    sigma_seas : (N, 3, 3) ndarray
        σ_seas in compression-POSITIVE SEAS convention, MPa.

    Raises
    ------
    ValueError
        For an unrecognised ``depth_model``.
    """
    centroids_z = np.asarray(centroids_z, dtype=np.float64)
    if centroids_z.ndim != 1:
        raise ValueError(
            f"centroids_z must be 1-D; got shape {centroids_z.shape}"
        )
    if depth_model not in ("constant", "lithostatic_sv"):
        raise ValueError(
            f"depth_model must be 'constant' or 'lithostatic_sv'; "
            f"got {depth_model!r}"
        )

    N = centroids_z.shape[0]
    sigma_seas = np.empty((N, 3, 3), dtype=np.float64)

    if depth_model == "constant":
        sigma_HZ = build_bulk_stress_tensor(
            SHmax_top, Shmin_top, Sv_top, SHmax_az_deg
        )
        sigma_seas[:] = -np.asarray(sigma_HZ, dtype=np.float64)
    else:
        # lithostatic_sv: linear depth dependence on the principal
        # magnitudes. Depth d = max(0, -z) (above-surface cells use
        # the surface value).
        depths = np.maximum(0.0, -centroids_z)
        for i in range(N):
            d = depths[i]
            sh = SHmax_top + SHmax_grad * d
            sh_min = Shmin_top + Shmin_grad * d
            sv = Sv_top + Sv_grad * d
            row = build_bulk_stress_tensor(sh, sh_min, sv, SHmax_az_deg)
            sigma_seas[i] = -np.asarray(row, dtype=np.float64)

    # Defensive symmetry check on the full stack after the flip.
    # Skipped when empty (np.max on an empty reduction raises an
    # opaque numpy error — R-601). Also explicitly surfaces NaN
    # entries (which would silently bypass `> tol` via NaN > x =
    # False — R-605).
    #
    # R-702: NaN rows that correspond 1-to-1 with NaN centroids
    # are per-spec (Phase 3 §"Edge Cases": a degenerate triangle
    # in Phase 1 emits a NaN centroid; the rotation propagates
    # NaN through that row). Only raise when the NaN pattern is
    # unexpected (e.g., NaN-corrupted SHmax_top / Shmin_top /
    # Sv_top / *_grad).
    if N > 0:
        nan_row = np.isnan(sigma_seas).any(axis=(1, 2))
        nan_centroid = np.isnan(centroids_z)
        unexpected_nan = nan_row & ~nan_centroid
        if np.any(unexpected_nan) or np.any(np.isinf(sigma_seas)):
            raise ValueError(
                "bulk_stress_tensor_field produced non-finite (NaN or "
                "inf) entries on rows where centroids_z is finite; "
                "check SHmax_top / Shmin_top / Sv_top / *_grad for "
                "non-finite inputs"
            )
        # Symmetry check on the finite rows only (NaN rows are
        # per-spec and have NaN - NaN = NaN in the diff).
        finite_rows = ~nan_row
        if np.any(finite_rows):
            diff = (
                sigma_seas[finite_rows]
                - np.swapaxes(sigma_seas[finite_rows], -1, -2)
            )
            asym = float(np.max(np.abs(diff)))
        else:
            asym = 0.0
        if asym > _SIGMA_SYMMETRY_TOL_MPA:
            raise ValueError(
                f"bulk_stress_tensor_field produced non-symmetric "
                f"tensors; max off-diagonal asymmetry {asym:.3e}"
            )
    return sigma_seas


def pore_pressure_field(
    centroids_z: np.ndarray,
    *,
    P_p_top: float = 0.0,
    P_p_grad: float = 0.0,
) -> np.ndarray:
    """Build a per-cell pore-pressure field with optional linear depth scaling.

    Returns ``P_p_top + P_p_grad * max(0, -z)`` per cell. Default
    (P_p_top = 0, P_p_grad = 0) reproduces "no pore pressure"; the
    H&Z demo case (P_p = 16 MPa applied uniformly) is recovered by
    passing ``P_p_top=16.0, P_p_grad=0.0``.

    Parameters
    ----------
    centroids_z : (N,) ndarray
        z-coordinates of the cells (≤ 0 in the subsurface).
    P_p_top : float, default 0.0
        Pore pressure at z = 0 (POSITIVE value, MPa).
    P_p_grad : float, default 0.0
        Linear depth gradient (MPa/m).

    Returns
    -------
    P_p : (N,) ndarray
        Per-cell pore pressure (MPa).
    """
    centroids_z = np.asarray(centroids_z, dtype=np.float64)
    if centroids_z.ndim != 1:
        raise ValueError(
            f"centroids_z must be 1-D; got shape {centroids_z.shape}"
        )
    depths = np.maximum(0.0, -centroids_z)
    return P_p_top + P_p_grad * depths


def project_stress_onto_fault(
    fault_geom: FaultCellGeometry,
    *,
    SHmax: float,
    Shmin: float,
    Sv: float,
    P_p: float,
    SHmax_az_deg: float,
    depth_model: str = "constant",
    SHmax_grad: float = 0.0,
    Shmin_grad: float = 0.0,
    Sv_grad: float = 0.0,
    P_p_grad: float = 0.0,
) -> dict:
    """Wire-up: build σ_seas, build P_p, rotate onto the fault basis.

    Convenience entry point that combines `bulk_stress_tensor_field`
    (Phase 3 §2) + `pore_pressure_field` (Phase 3 §3) +
    `resolve_traction_per_cell` (Phase 3 §1) into a single call.

    Returns the same dict as `resolve_traction_per_cell` plus a
    ``"sigma_field"`` key (shape ``(N_tri, 3, 3)``) so the bulk path
    (Phase 4) can reuse the rotation results without re-evaluating
    the depth model.

    Parameters
    ----------
    fault_geom : FaultCellGeometry
        Per-triangle Tandem basis + centroids from
        ``build_fault_basis``.
    SHmax, Shmin, Sv, P_p : float
        Principal-stress magnitudes and pore pressure at z = 0
        (POSITIVE values, MPa).
    SHmax_az_deg : float
        Geological azimuth of SHmax (clockwise from north, degrees).
    depth_model : {"constant", "lithostatic_sv"}
    SHmax_grad, Shmin_grad, Sv_grad, P_p_grad : float, default 0.0
        Linear depth gradients (MPa/m).

    Returns
    -------
    dict
        ``resolve_traction_per_cell`` output + ``"sigma_field"``
        ``(N_tri, 3, 3)``.
    """
    # R-602: warn when caller passes depth gradients but did not
    # opt into the depth-dependent model. Symmetric with Phase 5's
    # warning for `lithostatic_sv` + all-zero gradients.
    if depth_model == "constant":
        nonzero_grads = [
            (name, g) for name, g in (
                ("SHmax_grad", SHmax_grad),
                ("Shmin_grad", Shmin_grad),
                ("Sv_grad", Sv_grad),
                ("P_p_grad", P_p_grad),
            ) if g != 0.0
        ]
        if nonzero_grads:
            print(
                f"warning: project_stress_onto_fault received "
                f"non-zero depth gradients ({nonzero_grads}) but "
                f"depth_model='constant'; gradients are silently "
                f"ignored. Pass depth_model='lithostatic_sv' to "
                f"apply them.",
                file=sys.stderr,
            )
    centroids_z = fault_geom.centroids[:, 2]
    sigma_field = bulk_stress_tensor_field(
        centroids_z,
        SHmax_az_deg=SHmax_az_deg,
        depth_model=depth_model,
        SHmax_top=SHmax,
        Shmin_top=Shmin,
        Sv_top=Sv,
        SHmax_grad=SHmax_grad,
        Shmin_grad=Shmin_grad,
        Sv_grad=Sv_grad,
    )
    P_p_arr = pore_pressure_field(
        centroids_z, P_p_top=P_p, P_p_grad=P_p_grad,
    )
    resolved = resolve_traction_per_cell(
        sigma_field,
        fault_geom.strikes,
        fault_geom.dips,
        fault_geom.normals,
        P_p_per_cell=P_p_arr,
    )
    resolved["sigma_field"] = sigma_field
    return resolved


# ----------------------------------------------------------------------
# Phase 4 — VTU + JSON writers + CLI
# ----------------------------------------------------------------------

# H&Z `demo_safod` SAFOD reference parameters are introspected from
# the `dump_safod_sigma0` function signature so the H&Z module
# remains the single source of truth (no hardcoded numerics here —
# user-feedback memory `feedback_no_hardcoded_numbers.md`).
import inspect as _inspect

from hickman_and_zoback_regional_stress_projection import (  # noqa: E402
    dump_safod_sigma0 as _hz_dump_safod_sigma0,
)


def _hz_demo_safod_defaults() -> dict:
    """Introspect H&Z's `dump_safod_sigma0` for the SAFOD demo defaults."""
    sig = _inspect.signature(_hz_dump_safod_sigma0)
    return {
        "SHmax": sig.parameters["SHmax"].default,
        "Shmin": sig.parameters["Shmin"].default,
        "Sv": sig.parameters["Sv"].default,
        "P_p": sig.parameters["P_p"].default,
        "SHmax_az_deg": sig.parameters["SHmax_azimuth_deg"].default,
        "fault_strike_az_deg": (
            sig.parameters["fault_strike_azimuth_deg"].default
        ),
        "rake_sense": sig.parameters["rake_sense"].default,
    }


DEFAULTS = _hz_demo_safod_defaults()

# Suffix conventions for the output files. Mirror `nw_cut_strip.py`'s
# DEFAULT_SUFFIX pattern.
_FAULT_VTU_SUFFIX: str = "_fault_stress.vtu"
_BULK_VTU_SUFFIX: str = "_bulk_stress.vtu"
_SUMMARY_JSON_SUFFIX: str = "_summary.json"

# JSON summary "convention" string — pinned by Phase 4's plan §3.
_CONVENTION_STRING: str = (
    "compression POSITIVE (SEAS internal); P_p positive; "
    "sigma_n_eff = sigma_n_total - P_p; "
    "tau_strike + = right-lateral (under Tandem fault basis "
    "s = up x n, d = s x n); frame (east, north, up) UTM Zone 11N; "
    "units MPa"
)


def write_fault_vtu(
    fault_mesh,                          # meshio.Mesh
    tri_conn: np.ndarray,
    fault_geom: FaultCellGeometry,
    resolved: dict,
    out_path: Path,
    convention_str: str = _CONVENTION_STRING,
) -> int:
    """Write the fault VTU with continuous (point-data) primary
    output and raw (cell-data) `_cell` secondary channel.

    Continuous field: per-vertex values are produced by
    area-weighted node-averaging (`cell_to_node_average`) of the
    cell-centred resolved values. The basis at each vertex is
    Gram-Schmidt re-orthonormalised via `basis_to_node` (Tandem
    convention `d = s × n`, R-301).

    **`rake_deg`, `tau_magnitude_MPa`, `mu_apparent` are RECOMPUTED
    at the node from the node-averaged shear / normal components,
    NOT averaged directly.** Averaging an angle wraps the periodic
    ±180° discontinuity; averaging a magnitude double-counts the
    cell-level cancellation; averaging a ratio inflates the result.

    Parameters
    ----------
    fault_mesh : meshio.Mesh
        Input fault VTU (point cloud reused for the output).
    tri_conn : (N_tri, 3) ndarray
        Fault-phys triangle connectivity into ``fault_mesh.points``.
    fault_geom : FaultCellGeometry
        Per-cell basis + centroids + areas (Phase 2 output).
    resolved : dict
        Output of ``project_stress_onto_fault`` — cell-centred
        compression-positive SEAS values.
    out_path : Path
        Destination VTU path. Parent directory is created if missing.
    convention_str : str
        Metadata string attached as ``mesh.info["stress_convention"]``.

    Returns
    -------
    int
        Number of NaN-tainted vertices in the point-data (vertices
        where the basis is degenerate after node averaging, or all
        incident cells had NaN values).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    points = np.asarray(fault_mesh.points, dtype=np.float64)
    tri_conn = np.asarray(tri_conn, dtype=np.int64)
    N_tri = tri_conn.shape[0]
    N_pt = points.shape[0]

    # --- Cell→Node: scalar fields ---
    sigma_n_total_node = cell_to_node_average(
        points, tri_conn, resolved["sigma_n_total"], fault_geom.areas,
    )
    sigma_n_eff_node = cell_to_node_average(
        points, tri_conn, resolved["sigma_n_eff"], fault_geom.areas,
    )
    tau_strike_node = cell_to_node_average(
        points, tri_conn, resolved["tau_strike"], fault_geom.areas,
    )
    tau_dip_node = cell_to_node_average(
        points, tri_conn, resolved["tau_dip"], fault_geom.areas,
    )

    # --- Cell→Node: vector fields ---
    traction_vec_node = cell_to_node_average(
        points, tri_conn, resolved["traction_vec"], fault_geom.areas,
    )

    # --- Basis at vertex (Gram–Schmidt) ---
    s_node, d_node, n_node = basis_to_node(
        points, tri_conn,
        fault_geom.strikes, fault_geom.dips, fault_geom.normals,
        fault_geom.areas,
    )

    # --- Recompute periodic / non-linear quantities at the node ---
    # tau_magnitude: |τ| = sqrt(τ_s² + τ_d²); avoids the
    # averaging-magnitude-double-counting bug.
    tau_magnitude_node = np.sqrt(
        tau_strike_node ** 2 + tau_dip_node ** 2
    )
    # rake_deg: atan2(τ_d, τ_s); avoids periodic-average wrapping.
    rake_deg_node = np.degrees(
        np.arctan2(tau_dip_node, tau_strike_node)
    )
    # mu_apparent: |τ| / σ_n_eff with the same > 0 guard as Phase 3.
    mu_apparent_node = np.full(N_pt, np.nan, dtype=np.float64)
    finite = np.isfinite(sigma_n_eff_node) & (sigma_n_eff_node > 0.0)
    if finite.any():
        mu_apparent_node[finite] = (
            tau_magnitude_node[finite] / sigma_n_eff_node[finite]
        )

    # --- Cell-data: raw scalars + carried-through gmsh tags ---
    cell_data: dict = {
        "sigma_n_total_MPa_cell": [resolved["sigma_n_total"]],
        "sigma_n_eff_MPa_cell": [resolved["sigma_n_eff"]],
        "tau_strike_MPa_cell": [resolved["tau_strike"]],
        "tau_dip_MPa_cell": [resolved["tau_dip"]],
    }
    # Carry through gmsh:physical / gmsh:geometrical if present in
    # the input mesh. Note: the fault VTU we write has a single
    # triangle block (the fault group), so we pick the corresponding
    # cell-data slice that matches the tri_conn rows. The simplest
    # robust approach is to look up the per-cell tags by mask.
    input_cell_data = fault_mesh.cell_data or {}
    if "gmsh:physical" in input_cell_data:
        _carry = _slice_cell_data_for_phys(
            fault_mesh, input_cell_data["gmsh:physical"], tri_conn,
        )
        if _carry is not None:
            cell_data["gmsh:physical"] = [_carry]
    if "gmsh:geometrical" in input_cell_data:
        _carry = _slice_cell_data_for_phys(
            fault_mesh, input_cell_data["gmsh:geometrical"], tri_conn,
        )
        if _carry is not None:
            cell_data["gmsh:geometrical"] = [_carry]

    # --- Point-data: continuous primary field ---
    point_data: dict = {
        "sigma_n_total_MPa": sigma_n_total_node,
        "sigma_n_eff_MPa": sigma_n_eff_node,
        "tau_strike_MPa": tau_strike_node,
        "tau_dip_MPa": tau_dip_node,
        "tau_magnitude_MPa": tau_magnitude_node,
        "rake_deg": rake_deg_node,
        "mu_apparent": mu_apparent_node,
        "strike_vec": s_node,
        "dip_vec": d_node,
        "normal_vec": n_node,
        "traction_vec_MPa": traction_vec_node,
    }

    # Count NaN-tainted vertices that are actually REFERENCED by a
    # fault triangle (true degeneracy). Orphan vertices — those that
    # do not appear in `tri_conn` because msh_to_vtu.py emits the
    # fault VTU with the FULL bulk point cloud — are correctly NaN
    # by spec and excluded from this count.
    referenced = np.zeros(N_pt, dtype=bool)
    referenced[np.unique(tri_conn.flatten())] = True
    nan_vertex_count = int(
        (np.isnan(s_node).any(axis=1) & referenced).sum()
    )

    # Build and write the mesh.
    fault_out = meshio.Mesh(
        points=points,
        cells=[("triangle", tri_conn)],
        point_data=point_data,
        cell_data=cell_data,
    )
    # Metadata.
    fault_out.info = {
        "stress_convention": convention_str,
        "units": "MPa",
        "frame": "(east, north, up) / UTM Zone 11N",
    }
    meshio.write(str(out_path), fault_out)
    return nan_vertex_count


def _slice_cell_data_for_phys(
    fault_mesh,
    phys_arrays_list,
    tri_conn: np.ndarray,
) -> Optional[np.ndarray]:
    """Recover the per-row gmsh:physical / gmsh:geometrical tag values
    for the (N_tri, 3) triangles we're writing out.

    The input meshio.Mesh has one cell-data array per cell block; the
    triangle block's array gives the per-triangle tag. Our `tri_conn`
    rows were extracted from that block via
    `extract_cells_by_physical`, which may or may not have filtered.
    For the SAFS use case (single triangle block, single fault group,
    fault_mesh has empty field_data), every triangle in the block
    was kept; the array maps 1-to-1 onto `tri_conn` rows.

    Returns the matching slice (shape (N_tri,)) or None if no match.
    """
    # Find the triangle block index.
    for i, cb in enumerate(fault_mesh.cells):
        if cb.type == "triangle" and len(cb.data) == tri_conn.shape[0]:
            return np.asarray(phys_arrays_list[i])
    return None


def write_bulk_vtu(
    bulk_mesh,                           # meshio.Mesh
    tet_conn: np.ndarray,
    bulk_geom: BulkCellGeometry,
    sigma_field: np.ndarray,
    out_path: Path,
    convention_str: str = _CONVENTION_STRING,
) -> None:
    """Write the bulk VTU as a pass-through cell-data emission.

    Pass-through writer (R-501/R-502): ``sigma_field`` arrives in
    compression-positive SEAS convention from Phase 3's
    `bulk_stress_tensor_field`. **No sign flip applied** — split the
    (N_tet, 3, 3) tensor into six scalar `sigma_*_MPa` cell-data
    fields plus a `sigma_tensor_MPa` (N_tet, 9) flat tensor for
    ParaView's TensorGlyph filter.

    Important note re: dynamic-rupture use — the bulk stress field
    is NOT consumed by the dynamic-rupture run because that run
    assumes ∇·σ⁰ ≈ 0 (equilibrium); this artefact is for the
    downstream static-solve step.

    Parameters
    ----------
    bulk_mesh : meshio.Mesh
        Input bulk VTU (point cloud reused for the output).
    tet_conn : (N_tet, 4) ndarray
        Bulk-phys tetrahedral connectivity.
    bulk_geom : BulkCellGeometry
        Per-cell centroids and volumes (Phase 1 output). Not used in
        the writer body but reserved for downstream attachment.
    sigma_field : (N_tet, 3, 3) ndarray
        Cell-centred Cauchy tensors, compression-POSITIVE SEAS.
    out_path : Path
        Destination VTU path. Parent directory created if missing.
    convention_str : str
        Metadata string.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    points = np.asarray(bulk_mesh.points, dtype=np.float64)
    tet_conn = np.asarray(tet_conn, dtype=np.int64)
    sigma_field = np.asarray(sigma_field, dtype=np.float64)
    N_tet = tet_conn.shape[0]
    if sigma_field.shape != (N_tet, 3, 3):
        raise ValueError(
            f"sigma_field must be (N_tet={N_tet}, 3, 3); got "
            f"{sigma_field.shape}"
        )
    # Suppress the unused-argument warning for bulk_geom (reserved
    # for future attachment of centroid / volume cell-data).
    _ = bulk_geom

    # Cell-data: six scalar components + 9-component flat tensor.
    cell_data: dict = {
        "sigma_xx_MPa": [sigma_field[:, 0, 0]],
        "sigma_yy_MPa": [sigma_field[:, 1, 1]],
        "sigma_zz_MPa": [sigma_field[:, 2, 2]],
        "sigma_xy_MPa": [sigma_field[:, 0, 1]],
        "sigma_xz_MPa": [sigma_field[:, 0, 2]],
        "sigma_yz_MPa": [sigma_field[:, 1, 2]],
        "sigma_tensor_MPa": [sigma_field.reshape(N_tet, 9)],
    }
    input_cell_data = bulk_mesh.cell_data or {}
    if "gmsh:physical" in input_cell_data:
        _carry = _slice_cell_data_for_tet(
            bulk_mesh, input_cell_data["gmsh:physical"], tet_conn,
        )
        if _carry is not None:
            cell_data["gmsh:physical"] = [_carry]
    if "gmsh:geometrical" in input_cell_data:
        _carry = _slice_cell_data_for_tet(
            bulk_mesh, input_cell_data["gmsh:geometrical"], tet_conn,
        )
        if _carry is not None:
            cell_data["gmsh:geometrical"] = [_carry]

    bulk_out = meshio.Mesh(
        points=points,
        cells=[("tetra", tet_conn)],
        cell_data=cell_data,
    )
    bulk_out.info = {
        "stress_convention": convention_str,
        "units": "MPa",
        "frame": "(east, north, up) / UTM Zone 11N",
    }
    meshio.write(str(out_path), bulk_out)


def _slice_cell_data_for_tet(
    bulk_mesh,
    phys_arrays_list,
    tet_conn: np.ndarray,
) -> Optional[np.ndarray]:
    """Analog of `_slice_cell_data_for_phys` for tetrahedra."""
    for i, cb in enumerate(bulk_mesh.cells):
        if cb.type == "tetra" and len(cb.data) == tet_conn.shape[0]:
            return np.asarray(phys_arrays_list[i])
    return None


def _field_stats(arr: np.ndarray) -> dict:
    """Compute min / median / max / nan_count for a 1-D numeric
    field. Used by `write_summary_json` per the schema."""
    arr = np.asarray(arr, dtype=np.float64).ravel()
    finite_mask = np.isfinite(arr)
    nan_count = int(arr.size - finite_mask.sum())
    if finite_mask.any():
        finite_vals = arr[finite_mask]
        return {
            "min": float(np.min(finite_vals)),
            "median": float(np.median(finite_vals)),
            "max": float(np.max(finite_vals)),
            "nan_count": nan_count,
        }
    return {
        "min": float("nan"),
        "median": float("nan"),
        "max": float("nan"),
        "nan_count": nan_count,
    }


def write_summary_json(
    out_path: Path,
    *,
    input_mesh_path: Path,
    fault_geom: FaultCellGeometry,
    bulk_geom: Optional[BulkCellGeometry],
    sigma_global_at_z0: np.ndarray,
    params: dict,
    resolved: dict,
) -> None:
    """Write the per-run JSON summary.

    Schema follows Phase 4 §3 — top-level keys: ``input_mesh``,
    ``fault_n_cells``, ``fault_n_degenerate_basis``, ``bulk_n_cells``,
    ``params``, ``sigma_global_at_z0_MPa``, ``stats``, ``convention``.
    Uses ``json.dump(..., indent=2, allow_nan=True)`` so NaN
    ``mu_apparent`` entries serialise correctly.

    Parameters
    ----------
    out_path : Path
        Destination JSON path. Parent directory created if missing.
    input_mesh_path : Path
        Path of the input fault VTU (recorded as ``input_mesh``).
    fault_geom : FaultCellGeometry
        Phase 2 output (provides cell counts).
    bulk_geom : BulkCellGeometry or None
        Phase 1 output for the bulk mesh; ``None`` if bulk was not
        loaded (``--write-bulk`` not set).
    sigma_global_at_z0 : (3, 3) ndarray
        σ_seas at the free surface, MPa, compression-positive.
    params : dict
        Full CLI input echo (parameter values).
    resolved : dict
        Output of ``project_stress_onto_fault``; ``stats`` are
        computed over its cell-centred arrays.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sigma_3x3 = np.asarray(sigma_global_at_z0, dtype=np.float64)
    if sigma_3x3.shape != (3, 3):
        raise ValueError(
            f"sigma_global_at_z0 must be (3, 3); got {sigma_3x3.shape}"
        )

    payload = {
        "input_mesh": str(input_mesh_path),
        "fault_n_cells": int(fault_geom.strikes.shape[0]),
        "fault_n_degenerate_basis": int(fault_geom.n_degenerate),
        "bulk_n_cells": (
            int(bulk_geom.centroids.shape[0]) if bulk_geom is not None
            else None
        ),
        "params": {k: _json_safe(v) for k, v in params.items()},
        "sigma_global_at_z0_MPa": [
            [float(sigma_3x3[i, j]) for j in range(3)] for i in range(3)
        ],
        "stats": {
            "sigma_n_total_MPa": _field_stats(resolved["sigma_n_total"]),
            "sigma_n_eff_MPa": _field_stats(resolved["sigma_n_eff"]),
            "tau_strike_MPa": _field_stats(resolved["tau_strike"]),
            "tau_dip_MPa": _field_stats(resolved["tau_dip"]),
            "tau_magnitude_MPa": _field_stats(resolved["tau_magnitude"]),
            "rake_deg": _field_stats(resolved["rake_deg"]),
            "mu_apparent": _field_stats(resolved["mu_apparent"]),
        },
        "convention": _CONVENTION_STRING,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, allow_nan=True)


def _json_safe(value):
    """Convert numpy scalars / Path objects to JSON-serialisable types."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    return value


# ----------------------------------------------------------------------
# Phase 4 — CLI
# ----------------------------------------------------------------------


def _load_hz_dump_file(path: Path) -> dict:
    """Load and validate a Phase 0 H&Z dump JSON.

    Returns the ``params`` block; raises ``ValueError`` on schema
    mismatch.
    """
    with open(path) as f:
        data = json.load(f)
    if data.get("schema") != "hickman_zoback_sigma0_v1":
        raise ValueError(
            f"H&Z dump file {path} has unexpected schema "
            f"{data.get('schema')!r}; expected "
            f"'hickman_zoback_sigma0_v1'"
        )
    return data.get("params", {})


def _build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Project a regional bulk Cauchy stress tensor onto the "
            "SAFS fault and bulk meshes. Outputs are ParaView VTU "
            "files plus a JSON summary in compression-positive SEAS "
            "convention. See PLAN_onfaultstress.md."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "input_mesh", nargs="?", type=Path,
        help="path to a *_fault.vtu file (single-file mode)",
    )
    ap.add_argument(
        "mesh_base", nargs="?", type=str,
        help="output filename stem for single-file mode "
             "(e.g. 'safs_fault_box_nwcut_2000m')",
    )
    ap.add_argument(
        "--batch", action="store_true",
        help="batch-process every mesh in DEFAULT_INPUT_BASES",
    )
    ap.add_argument(
        "--print-info", type=Path, default=None, metavar="PATH",
        help="print σ⁰ and fault cell statistics, then exit "
             "without writing any artefacts",
    )
    ap.add_argument(
        "--bulk-name", type=str, default=DEFAULT_BULK_NAME,
        help="physical-group name of the bulk volume in the VTU "
             f"(default {DEFAULT_BULK_NAME!r})",
    )
    ap.add_argument(
        "--fault-name", type=str, default=DEFAULT_FAULT_NAME,
        help="physical-group name of the fault surface in the VTU "
             f"(default {DEFAULT_FAULT_NAME!r})",
    )
    ap.add_argument(
        "--hz-dump-file", type=Path, default=None,
        help="optional path to a Phase 0 H&Z dump JSON "
             "(hickman_zoback_sigma0_v1 schema); its params block "
             "populates SHmax/Shmin/Sv/P_p/SHmax_az unless "
             "explicit --flag overrides are also given.",
    )
    # Stress parameters with H&Z demo_safod defaults (via introspection).
    ap.add_argument("--SHmax", type=float, default=None,
                    help=f"SHmax (MPa, positive). Default "
                         f"{DEFAULTS['SHmax']} (H&Z SAFOD).")
    ap.add_argument("--Shmin", type=float, default=None,
                    help=f"Shmin (MPa). Default {DEFAULTS['Shmin']}.")
    ap.add_argument("--Sv", type=float, default=None,
                    help=f"Sv (MPa). Default {DEFAULTS['Sv']}.")
    ap.add_argument("--P_p", type=float, default=None,
                    help=f"Pore pressure (MPa). Default "
                         f"{DEFAULTS['P_p']}.")
    ap.add_argument("--SHmax-az", type=float, default=None,
                    dest="SHmax_az",
                    help=f"SHmax azimuth (deg cw from N). Default "
                         f"{DEFAULTS['SHmax_az_deg']}.")
    ap.add_argument(
        "--rake-sense", type=str, default=DEFAULTS["rake_sense"],
        choices=("right-lateral", "left-lateral"),
        help=f"rake convention. Default {DEFAULTS['rake_sense']!r}.",
    )
    ap.add_argument(
        "--strike-hint-az", type=float,
        default=DEFAULTS["fault_strike_az_deg"],
        help=f"strike-hint azimuth for orientation harmonisation. "
             f"Default {DEFAULTS['fault_strike_az_deg']} (SAF N46W).",
    )
    ap.add_argument(
        "--no-strike-hint", action="store_true",
        help="disable the strike-hint path and use PCA fallback "
             "for orientation harmonisation",
    )
    ap.add_argument(
        "--depth-model", type=str, default="constant",
        choices=("constant", "lithostatic_sv"),
        help="bulk σ⁰ depth model (default 'constant')",
    )
    ap.add_argument("--SHmax-grad", type=float, default=0.0,
                    dest="SHmax_grad",
                    help="SHmax depth gradient (MPa/m), used with "
                         "depth-model=lithostatic_sv (default 0)")
    ap.add_argument("--Shmin-grad", type=float, default=0.0,
                    dest="Shmin_grad",
                    help="Shmin depth gradient (MPa/m) (default 0)")
    ap.add_argument("--Sv-grad", type=float, default=0.0,
                    dest="Sv_grad",
                    help="Sv depth gradient (MPa/m) (default 0)")
    ap.add_argument("--P_p-grad", type=float, default=0.0,
                    dest="P_p_grad",
                    help="pore-pressure depth gradient (MPa/m) "
                         "(default 0)")
    ap.add_argument(
        "--write-bulk", action="store_true",
        help="also load the *_bulk.vtu sibling and write a "
             "<base>_bulk_stress.vtu artefact",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=DEFAULT_OUT_DIR,
        help=f"output directory (default {DEFAULT_OUT_DIR})",
    )
    return ap


def _resolve_stress_params(args) -> dict:
    """Resolve the SHmax / Shmin / Sv / P_p / SHmax_az parameters
    from (1) CLI overrides, (2) H&Z dump file, (3) DEFAULTS.

    Precedence: CLI flag > dump-file > DEFAULTS (H&Z demo_safod).
    """
    dump_params: dict = {}
    if args.hz_dump_file is not None:
        dump_params = _load_hz_dump_file(args.hz_dump_file)

    def _pick(cli_value, dump_key, default_key):
        if cli_value is not None:
            return cli_value
        if dump_key in dump_params:
            return dump_params[dump_key]
        return DEFAULTS[default_key]

    return {
        "SHmax": _pick(args.SHmax, "SHmax_MPa", "SHmax"),
        "Shmin": _pick(args.Shmin, "Shmin_MPa", "Shmin"),
        "Sv": _pick(args.Sv, "Sv_MPa", "Sv"),
        "P_p": _pick(args.P_p, "P_p_MPa", "P_p"),
        "SHmax_az_deg": _pick(
            args.SHmax_az, "SHmax_azimuth_deg", "SHmax_az_deg"
        ),
    }


def _print_info(input_mesh_path: Path, args) -> int:
    """Implement --print-info: load the fault mesh, print σ⁰ and
    fault-cell statistics, exit without writing."""
    print(f"Reading {input_mesh_path} ...")
    mesh = load_fault_mesh(input_mesh_path)
    p = _resolve_stress_params(args)
    print(f"  SHmax = {p['SHmax']} MPa, Shmin = {p['Shmin']} MPa, "
          f"Sv = {p['Sv']} MPa, P_p = {p['P_p']} MPa")
    print(f"  SHmax_az = {p['SHmax_az_deg']}° (cw from north)")
    sigma_seas_at_z0 = bulk_stress_tensor_field(
        np.array([0.0]),
        SHmax_az_deg=p["SHmax_az_deg"], depth_model="constant",
        SHmax_top=p["SHmax"], Shmin_top=p["Shmin"], Sv_top=p["Sv"],
    )[0]
    print("Bulk σ⁰ at z = 0 (compression positive, MPa):")
    np.set_printoptions(precision=3, suppress=True, sign="+")
    print(sigma_seas_at_z0)
    np.set_printoptions()  # restore
    tri_conn = extract_cells_by_physical(
        mesh, "triangle", args.fault_name,
    )
    print(f"Fault triangles: {tri_conn.shape[0]}")
    # R-703: strike-azimuth statistics per plan §Phase 4 §4 line
    # 1253-1255.  Built from the harmonised per-triangle strikes
    # (Tandem basis); azimuth is measured clockwise from north,
    # i.e. north = +y, east = +x.
    fault_geom = build_fault_basis(
        mesh, fault_phys_name=args.fault_name,
        rake_sense=args.rake_sense,
        fault_strike_azimuth_hint_deg=(
            None if args.no_strike_hint else args.strike_hint_az
        ),
    )
    s = fault_geom.strikes
    finite = ~np.isnan(s).any(axis=1)
    if finite.any():
        s_f = s[finite]
        az = (
            90.0 - np.degrees(np.arctan2(s_f[:, 1], s_f[:, 0]))
        ) % 360.0
        print(
            f"Strike-azimuth statistics: "
            f"min={az.min():.1f}°, "
            f"median={float(np.median(az)):.1f}°, "
            f"max={az.max():.1f}°, "
            f"stddev={az.std():.1f}°"
        )
    return 0


def _run_single(
    input_mesh_path: Path,
    mesh_base: str,
    args,
) -> int:
    """Implement single-file mode: load fault VTU, build basis,
    project stress, write the three artefacts under out_dir/<lc>m/.

    Returns 0 on success, non-zero on error.
    """
    print(f"Reading {input_mesh_path} ...")
    fault_mesh = load_fault_mesh(input_mesh_path)

    print("Building per-cell fault geometry + Tandem basis ...")
    strike_hint = (
        None if args.no_strike_hint else args.strike_hint_az
    )
    fault_geom = build_fault_basis(
        fault_mesh,
        fault_phys_name=args.fault_name,
        rake_sense=args.rake_sense,
        fault_strike_azimuth_hint_deg=strike_hint,
    )
    tri_conn = extract_cells_by_physical(
        fault_mesh, "triangle", args.fault_name,
    )
    print(f"  {tri_conn.shape[0]} triangles, "
          f"{fault_geom.n_degenerate} degenerate basis row(s)")

    # Resolve stress parameters (CLI > dump > defaults).
    p = _resolve_stress_params(args)
    if args.depth_model == "lithostatic_sv" and all(
        g == 0.0 for g in (args.SHmax_grad, args.Shmin_grad,
                            args.Sv_grad, args.P_p_grad)
    ):
        print(
            "warning: depth-model=lithostatic_sv with all gradients "
            "= 0 is equivalent to 'constant'",
            file=sys.stderr,
        )

    print(f"Rotating stress onto fault basis "
          f"({args.depth_model}) ...")
    resolved = project_stress_onto_fault(
        fault_geom,
        SHmax=p["SHmax"], Shmin=p["Shmin"], Sv=p["Sv"], P_p=p["P_p"],
        SHmax_az_deg=p["SHmax_az_deg"],
        depth_model=args.depth_model,
        SHmax_grad=args.SHmax_grad, Shmin_grad=args.Shmin_grad,
        Sv_grad=args.Sv_grad, P_p_grad=args.P_p_grad,
    )

    # Output directory: out_dir/<lc>m/  where <lc> is the trailing
    # "<N>m" tag in the mesh_base (e.g. "..._2000m" → "2000m/").
    lc_tag = _extract_lc_tag(mesh_base)
    out_dir = args.out_dir / lc_tag
    out_dir.mkdir(parents=True, exist_ok=True)
    fault_out = out_dir / f"{mesh_base}{_FAULT_VTU_SUFFIX}"
    bulk_out = out_dir / f"{mesh_base}{_BULK_VTU_SUFFIX}"
    summary_out = out_dir / f"{mesh_base}{_SUMMARY_JSON_SUFFIX}"

    print(f"Writing fault VTU {fault_out} ...")
    nan_vertex_count = write_fault_vtu(
        fault_mesh, tri_conn, fault_geom, resolved, fault_out,
    )
    print(f"  fault VTU NaN-tainted vertices: {nan_vertex_count}")

    bulk_geom: Optional[BulkCellGeometry] = None
    if args.write_bulk:
        bulk_path = _derive_bulk_path(input_mesh_path)
        if not bulk_path.is_file():
            raise FileNotFoundError(
                f"--write-bulk requested but bulk VTU not found: "
                f"{bulk_path}. Run msh_to_vtu.py first."
            )
        print(f"Reading bulk mesh {bulk_path} ...")
        bulk_mesh = load_bulk_mesh(bulk_path)
        tet_conn = extract_cells_by_physical(
            bulk_mesh, "tetra", args.bulk_name,
        )
        bulk_centroids, bulk_volumes = tet_geometry(
            bulk_mesh.points, tet_conn,
        )
        bulk_geom = BulkCellGeometry(
            centroids=bulk_centroids, volumes=bulk_volumes,
        )
        print(f"  {tet_conn.shape[0]} tetrahedra")
        print(f"Building bulk σ_seas field ({args.depth_model}) ...")
        bulk_sigma = bulk_stress_tensor_field(
            bulk_centroids[:, 2],
            SHmax_az_deg=p["SHmax_az_deg"],
            depth_model=args.depth_model,
            SHmax_top=p["SHmax"], Shmin_top=p["Shmin"],
            Sv_top=p["Sv"],
            SHmax_grad=args.SHmax_grad, Shmin_grad=args.Shmin_grad,
            Sv_grad=args.Sv_grad,
        )
        print(f"Writing bulk VTU {bulk_out} ...")
        write_bulk_vtu(bulk_mesh, tet_conn, bulk_geom, bulk_sigma,
                       bulk_out)

    # σ⁰ at z=0 for the summary (for the deprecated "constant" model
    # this matches every cell; for "lithostatic_sv" it's the surface
    # reference value).
    sigma_at_z0 = bulk_stress_tensor_field(
        np.array([0.0]),
        SHmax_az_deg=p["SHmax_az_deg"], depth_model="constant",
        SHmax_top=p["SHmax"], Shmin_top=p["Shmin"], Sv_top=p["Sv"],
    )[0]
    params_echo = {
        "SHmax_MPa": p["SHmax"],
        "Shmin_MPa": p["Shmin"],
        "Sv_MPa": p["Sv"],
        "P_p_MPa": p["P_p"],
        "SHmax_azimuth_deg": p["SHmax_az_deg"],
        "rake_sense": args.rake_sense,
        "depth_model": args.depth_model,
        "SHmax_grad_MPa_per_m": args.SHmax_grad,
        "Shmin_grad_MPa_per_m": args.Shmin_grad,
        "Sv_grad_MPa_per_m": args.Sv_grad,
        "P_p_grad_MPa_per_m": args.P_p_grad,
        "fault_strike_azimuth_hint_deg": (
            None if args.no_strike_hint else args.strike_hint_az
        ),
        "fault_vtu_nan_vertex_count": nan_vertex_count,
    }
    print(f"Writing summary JSON {summary_out} ...")
    write_summary_json(
        summary_out,
        input_mesh_path=input_mesh_path,
        fault_geom=fault_geom,
        bulk_geom=bulk_geom,
        sigma_global_at_z0=sigma_at_z0,
        params=params_echo,
        resolved=resolved,
    )
    print("Done.")
    return 0


def _derive_bulk_path(fault_path: Path) -> Path:
    """Given a `*_fault.vtu`, return the sibling `*_bulk.vtu`.

    Only the canonical `_fault.vtu` stem suffix is accepted; any
    other use of `_fault` raises (R-706: the prior substring-replace
    fallback could mangle paths like ``safs_fault_box.vtu`` →
    ``safs_bulk_box.vtu``, producing a sibling path that almost
    certainly does not exist).
    """
    name = fault_path.name
    if not name.endswith("_fault.vtu"):
        raise ValueError(
            f"cannot derive bulk path from {fault_path}: filename "
            f"does not end with '_fault.vtu' (plan §Phase 4 §5 "
            f"requires the `_fault.vtu` suffix)"
        )
    new_name = name[: -len("_fault.vtu")] + "_bulk.vtu"
    return fault_path.with_name(new_name)


def _extract_lc_tag(mesh_base: str) -> str:
    """Extract the trailing `<N>m` tag (optionally followed by a single
    variant suffix) from a mesh-base name.

    Examples
    --------
    ``safs_fault_box_nwcut_2000m`` → ``2000m``
    ``safs_fault_box_nwcut_500m_lcfar3000`` → ``500m_lcfar3000``
    ``safs_fault_box_nwcut_1000m_zgraded`` → ``1000m_zgraded``

    Falls back to the full ``mesh_base`` if no trailing ``<N>m`` is
    found.
    """
    import re
    m = re.search(r"_(\d+m(?:_[A-Za-z0-9]+)?)$", mesh_base)
    if m:
        return m.group(1)
    return mesh_base


def _run_batch(args) -> int:
    """Implement --batch: iterate over DEFAULT_INPUT_BASES, locating
    each ``meshing/results/vtu/<base>_fault.vtu`` and running the
    single-file path. Returns 0 on success."""
    code_meshing_dir = (
        Path(__file__).resolve().parent.parent.parent
        / "meshing" / "results" / "vtu"
    )
    for base in DEFAULT_INPUT_BASES:
        input_path = code_meshing_dir / f"{base}_fault.vtu"
        if not input_path.is_file():
            print(
                f"warning: skipping {base} — fault VTU not on disk: "
                f"{input_path}",
                file=sys.stderr,
            )
            continue
        print(f"=== batch mode: processing {base} ===")
        _run_single(input_path, base, args)
    return 0


def main(argv: Optional[list] = None) -> int:
    args = _build_argparser().parse_args(argv)

    has_print = args.print_info is not None
    has_batch = bool(args.batch)
    has_pos = (args.input_mesh is not None) or (args.mesh_base is not None)
    n_modes = sum((has_print, has_batch, has_pos))
    if n_modes > 1:
        print(
            "--print-info / --batch / (input_mesh, mesh_base) are "
            "mutually exclusive",
            file=sys.stderr,
        )
        return 2
    if n_modes == 0:
        print(
            "no mode given: use --print-info PATH, --batch, or "
            "input_mesh + mesh_base",
            file=sys.stderr,
        )
        return 2
    if has_pos and (args.input_mesh is None or args.mesh_base is None):
        print(
            "single-file mode requires BOTH input_mesh AND mesh_base",
            file=sys.stderr,
        )
        return 2

    if has_print:
        return _print_info(args.print_info, args)
    if has_batch:
        return _run_batch(args)
    return _run_single(args.input_mesh, args.mesh_base, args)


if __name__ == "__main__":
    sys.exit(main())
