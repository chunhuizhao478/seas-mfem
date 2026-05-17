"""fault_zone_metric.py — per-tet integrated L² of (proj − sidecar)
stratified by distance to the cut fault STL.

Ports the cubature + STL-distance stratification path out of
``plot_comparison_with_dg0.py:226-279`` into a matplotlib-free,
JSON-emitting CLI tool. Same numerics, same band definitions, same
4-point Gauss tet cubature.

CLI::

    python fault_zone_metric.py
        --vtu        PATH         (required)
        --sidecar    PATH         (required)
        --stl        PATH         (required)
        --field      Vs           (default Vs)
        --z-target   -1000.0      (default -1000 m)
        --half-thick 200.0        (default 200 m, matches v3 §0)
        --bands      0,500,1500,3000,inf
        --json       PATH         (optional)
        --parse-mfem-vtu PATH     (optional override of /tmp helper)

Importable surface (used by the Phase-2/3 acceptance tests in
test_fault_zone_metric.py)::

    per_tet_l2_by_band(vtu_path, sidecar_path, stl_path, *,
                       field='Vs', z_target=-1000.0, half_thick=200.0,
                       bands=DEFAULT_BANDS) -> list[BandRow]

Exit codes:
    0   success
    1   non-fatal error (bad CLI, missing field, …)
    2   sidecar bbox does not contain the mesh footprint at the slab z

This script intentionally avoids matplotlib so it imports on thin
remote nodes (e.g. Frontera login).  PNGs stay in
``plot_comparison_with_dg0.py``.
"""
from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from plot_comparison_with_dg0 import (                              # noqa: E402
    TET_CUB_BC,
    TET_CUB_W,
    trilinear,
    load_stl_vertices,
)


DEFAULT_BANDS: tuple[tuple[float, float], ...] = (
    (0.0, 500.0),
    (500.0, 1500.0),
    (1500.0, 3000.0),
    (3000.0, math.inf),
)


@dataclasses.dataclass(frozen=True)
class BandRow:
    """One distance-band aggregate.

    ``l2_proj`` is the RMS over the band of per-tet integrated L²
    residuals, in the same units as the underlying field (m/s for Vs).
    """
    lo: float
    hi: float
    n_tets: int
    l2_proj: float


# ---------------------------------------------------------------------
# Optional /tmp helper loader
# ---------------------------------------------------------------------

def _load_parse_vtu(parse_mfem_vtu: Path | None):
    """Load ``parse_mfem_vtu.parse_vtu`` exactly the way
    ``plot_comparison_with_dg0.py`` does (``sys.path.insert(0, '/tmp')``)
    so behaviour is reproducible across environments.  Allow an explicit
    override path for CI / sandboxed runs.
    """
    if parse_mfem_vtu is not None:
        path = Path(parse_mfem_vtu)
        if not path.is_file():
            raise FileNotFoundError(
                f"--parse-mfem-vtu: file not found: {path}")
        spec = importlib.util.spec_from_file_location(
            "parse_mfem_vtu_local", str(path))
        mod = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(mod)
        return mod.parse_vtu

    # Default: same as plot_comparison_with_dg0.py.
    if "/tmp" not in sys.path:
        sys.path.insert(0, "/tmp")
    try:
        from parse_mfem_vtu import parse_vtu                        # noqa: E402
    except ImportError as exc:
        raise RuntimeError(
            "missing parse_mfem_vtu helper; copy parse_mfem_vtu.py to "
            "/tmp/ or pass --parse-mfem-vtu PATH"
        ) from exc
    return parse_vtu


# ---------------------------------------------------------------------
# Per-tet evaluation kernels
# ---------------------------------------------------------------------

# MFEM Lagrange-tet (P2) edge ordering on the reference simplex
# vertices (v0, v1, v2, v3): edges enumerated as
#     (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).
# Closed-form basis on barycentric (l0, l1, l2, l3):
#     phi_corner_i = l_i (2 l_i - 1)               i in {0,1,2,3}
#     phi_edge_(i,j) = 4 l_i l_j                   for the 6 edges above
# Result is the 10-DOF Lagrange tet of order 2.
_P2_EDGES: tuple[tuple[int, int], ...] = (
    (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
)


def _p2_basis(bc: np.ndarray) -> np.ndarray:
    """Evaluate the 10 P2 tet basis functions at the supplied
    barycentric coordinates.

    Parameters
    ----------
    bc : (Q, 4) float
        Per-quadrature-point barycentric coords (l0, l1, l2, l3),
        ordered as (corner0, corner1, corner2, corner3).

    Returns
    -------
    phi : (Q, 10) float
        First 4 columns: corner basis.  Last 6 columns: edge basis in
        ``_P2_EDGES`` order.
    """
    if bc.ndim != 2 or bc.shape[1] != 4:
        raise ValueError(f"_p2_basis: bc shape must be (Q, 4); got "
                         f"{bc.shape}")
    Q = bc.shape[0]
    phi = np.empty((Q, 10), dtype=np.float64)
    # corners
    for i in range(4):
        phi[:, i] = bc[:, i] * (2.0 * bc[:, i] - 1.0)
    # edges
    for k, (i, j) in enumerate(_P2_EDGES):
        phi[:, 4 + k] = 4.0 * bc[:, i] * bc[:, j]
    return phi


def _p1_proj_at_cubature(vert_vals: np.ndarray, conn: np.ndarray
                         ) -> np.ndarray:
    """H1-P1 evaluation at the 4 cubature points of every tet.

    Identical to ``plot_comparison_with_dg0.py:248-249``::

        h1_at_q = einsum('qi,ti->tq', TET_CUB_BC, vert_vals[conn])

    Returns array of shape (n_tets, 4).
    """
    if conn.shape[1] != 4:
        raise ValueError(
            f"_p1_proj_at_cubature: P1 conn must be (n_tets, 4); got "
            f"{conn.shape}")
    Vverts = vert_vals[conn]                                # (n, 4)
    return np.einsum("qi,ti->tq", TET_CUB_BC, Vverts)       # (n, 4)


def _p2_proj_at_cubature(node_vals: np.ndarray,
                         conn_p2: np.ndarray) -> np.ndarray:
    """H1-P2 evaluation at the 4 cubature points of every tet.

    Parameters
    ----------
    node_vals : (Nn,) float
        Field value at every P2 node (corner + edge-mid).  Indexed by
        the same DOF id that appears in ``conn_p2``.
    conn_p2 : (n_tets, 10) int
        MFEM P2 connectivity, ordering [c0, c1, c2, c3, e01, e02,
        e03, e12, e13, e23].

    Returns
    -------
    p2_at_q : (n_tets, 4) float
    """
    if conn_p2.shape[1] != 10:
        raise ValueError(
            f"_p2_proj_at_cubature: P2 conn must be (n_tets, 10); got "
            f"{conn_p2.shape}")
    phi = _p2_basis(TET_CUB_BC)                             # (4, 10)
    # node values per element: (n_tets, 10)
    Vnodes = node_vals[conn_p2]
    # P2 value at every cubature point: (n_tets, 4)
    return np.einsum("qd,td->tq", phi, Vnodes)


# ---------------------------------------------------------------------
# .vtu helpers
# ---------------------------------------------------------------------

def _split_conn_by_offset(conn_flat: np.ndarray,
                          offsets: np.ndarray) -> tuple[int, np.ndarray]:
    """Detect P1 (4 DOFs) vs P2 (10 DOFs) tet connectivity from an
    MFEM-emitted vtu.  Returns (dofs_per_elem, conn_2d).

    Raises ValueError on a heterogeneous (mixed cell-type) mesh — the
    SAFS preview pipeline never emits one.
    """
    if offsets.size == 0:
        raise ValueError(
            "_split_conn_by_offset: empty cell-offset array")
    cell_lengths = np.diff(np.concatenate([[0], offsets]))
    unique_lengths = np.unique(cell_lengths)
    if unique_lengths.size != 1:
        raise ValueError(
            "_split_conn_by_offset: non-homogeneous cell connectivity "
            f"lengths {unique_lengths.tolist()}; this CLI requires a "
            "single tet order across the mesh")
    dofs = int(unique_lengths[0])
    if dofs not in (4, 10):
        raise ValueError(
            f"_split_conn_by_offset: unsupported DOFs/elem={dofs}; "
            f"only 4 (P1) and 10 (P2) tets are supported")
    n_tets = offsets.size
    conn_2d = conn_flat.reshape(n_tets, dofs)
    return dofs, conn_2d


def _mesh_xyz_bbox(pts: np.ndarray) -> tuple[float, float, float,
                                              float, float, float]:
    """Return (xmin, xmax, ymin, ymax, zmin, zmax) over a (N, 3) array."""
    return (float(pts[:, 0].min()), float(pts[:, 0].max()),
            float(pts[:, 1].min()), float(pts[:, 1].max()),
            float(pts[:, 2].min()), float(pts[:, 2].max()))


def _print_bbox_violation(field_name: str,
                          mxyz: tuple[float, float, float, float,
                                      float, float],
                          dxyz: tuple[float, float, float, float,
                                      float, float]) -> None:
    """Print a violation-message in the same shape as
    ``io/field_coefficient.cpp::AbortContainmentFailure``."""
    print(f"ERROR: sidecar field '{field_name}' bbox does NOT "
          f"contain the mesh footprint at the slab z.")
    print(f"  mesh bbox  (UTM 11 N, m): "
          f"x=[{mxyz[0]}, {mxyz[1]}] "
          f"y=[{mxyz[2]}, {mxyz[3]}] "
          f"z=[{mxyz[4]}, {mxyz[5]}]")
    print(f"  data bbox  (UTM 11 N, m): "
          f"x=[{dxyz[0]}, {dxyz[1]}] "
          f"y=[{dxyz[2]}, {dxyz[3]}] "
          f"z=[{dxyz[4]}, {dxyz[5]}]")
    print("v1 schema enforces interpolation-only; rebuild the sidecar "
          "or reduce the mesh padding.")


# ---------------------------------------------------------------------
# Public function (used by tests + CLI)
# ---------------------------------------------------------------------

def per_tet_l2_by_band(
    vtu_path: Path,
    sidecar_path: Path,
    stl_path: Path,
    *,
    field: str = "Vs",
    z_target: float = -1000.0,
    half_thick: float = 200.0,
    bands: Iterable[tuple[float, float]] = DEFAULT_BANDS,
    parse_mfem_vtu: Path | None = None,
) -> list[BandRow]:
    """Compute the per-tet integrated L² of (proj − sidecar) on the
    slab ``z = z_target ± half_thick``, stratified by distance to the
    fault STL into ``bands``.

    Per-band aggregate is RMS over the tets in the band of the per-tet
    integrated L² residuals.

    Returns one :class:`BandRow` per supplied band.
    """
    parse_vtu = _load_parse_vtu(parse_mfem_vtu)

    bands_t = tuple(bands)
    if not bands_t:
        raise ValueError("per_tet_l2_by_band: bands must be non-empty")

    vtu_path = Path(vtu_path)
    sidecar_path = Path(sidecar_path)
    stl_path = Path(stl_path)
    if not vtu_path.is_file():
        raise FileNotFoundError(f"--vtu: file not found: {vtu_path}")
    if not sidecar_path.is_file():
        raise FileNotFoundError(
            f"--sidecar: file not found: {sidecar_path}")
    if not stl_path.is_file():
        raise FileNotFoundError(f"--stl: file not found: {stl_path}")

    # Sidecar.
    with h5py.File(sidecar_path, "r") as h5:
        if "fields" not in h5 or field not in h5["fields"]:
            available = (list(h5["fields"].keys())
                         if "fields" in h5 else [])
            raise KeyError(
                f"sidecar '{sidecar_path}' has no field '{field}'; "
                f"available fields: {available}")
        gx = np.asarray(h5["grid/x"][...], dtype=np.float64)
        gy = np.asarray(h5["grid/y"][...], dtype=np.float64)
        gz = np.asarray(h5["grid/z"][...], dtype=np.float64)
        F = np.asarray(h5[f"fields/{field}"][...], dtype=np.float64)

    # vtu.
    m = parse_vtu(str(vtu_path))
    pts = np.ascontiguousarray(m["points"], dtype=np.float64)
    if pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError(f"vtu has unexpected points shape {pts.shape}")
    conn_flat = np.asarray(m["conn"], dtype=np.int64)
    offsets = np.asarray(m["offsets"], dtype=np.int64)
    dofs, conn_2d = _split_conn_by_offset(conn_flat, offsets)
    pd = m["point_data"]
    if field not in pd:
        raise KeyError(
            f"vtu has no PointData '{field}'; available: "
            f"{list(pd.keys())}")
    vert_vals = np.asarray(pd[field], dtype=np.float64)

    # Slab mask via tet centroids (pts at corners 0..3 always exist on
    # both P1 and P2 — the first 4 columns of conn_2d are the corner
    # vertices for any MFEM Lagrange tet).
    corner_conn = conn_2d[:, :4]
    V_xyz = pts[corner_conn]                               # (n_t_all, 4, 3)
    centroids = V_xyz.mean(axis=1)                         # (n_t_all, 3)
    mask = ((centroids[:, 2] >= z_target - half_thick)
            & (centroids[:, 2] <= z_target + half_thick))

    if not mask.any():
        return [BandRow(lo, hi, 0, float("nan")) for (lo, hi) in bands_t]

    V_xyz = V_xyz[mask]
    conn_slab = conn_2d[mask]
    centroids = centroids[mask]
    n_t = V_xyz.shape[0]

    # --- Strict-interpolation guard against the sidecar bbox ----------
    mxyz = _mesh_xyz_bbox(V_xyz.reshape(-1, 3))
    dxyz = (float(gx[0]), float(gx[-1]),
            float(gy[0]), float(gy[-1]),
            float(gz[0]), float(gz[-1]))
    if not (dxyz[0] <= mxyz[0] and mxyz[1] <= dxyz[1] and
            dxyz[2] <= mxyz[2] and mxyz[3] <= dxyz[3] and
            dxyz[4] <= mxyz[4] and mxyz[5] <= dxyz[5]):
        _print_bbox_violation(field, mxyz, dxyz)
        raise BBoxOutOfRange(
            f"sidecar bbox does not contain the slab mesh footprint")

    # --- Cubature points -------------------------------------------------
    qpts = np.einsum("qi,tij->tqj", TET_CUB_BC, V_xyz)     # (n, 4, 3)
    f_truth = trilinear(gx, gy, gz, F,
                        qpts.reshape(-1, 3)).reshape(n_t, 4)

    # --- Projection at the same cubature points -------------------------
    if dofs == 4:
        proj_at_q = _p1_proj_at_cubature(vert_vals, conn_slab)
    else:  # dofs == 10
        proj_at_q = _p2_proj_at_cubature(vert_vals, conn_slab)

    res = proj_at_q - f_truth
    l2_per_tet = np.sqrt((res ** 2 * TET_CUB_W).sum(axis=1))

    # --- STL distance ---------------------------------------------------
    stl_v = load_stl_vertices(str(stl_path))
    if stl_v is None or len(stl_v) == 0:
        raise ValueError(f"STL '{stl_path}' parsed empty")
    stl_xy = stl_v.mean(axis=1)[:, :2]
    tree = cKDTree(stl_xy)
    d, _ = tree.query(centroids[:, :2], k=1)

    rows: list[BandRow] = []
    for lo, hi in bands_t:
        sel = (d >= lo) & (d < hi)
        n_sel = int(sel.sum())
        if n_sel == 0:
            rows.append(BandRow(lo=float(lo), hi=float(hi),
                                n_tets=0, l2_proj=float("nan")))
            continue
        rms = float(np.sqrt(np.mean(l2_per_tet[sel] ** 2)))
        rows.append(BandRow(lo=float(lo), hi=float(hi),
                            n_tets=n_sel, l2_proj=rms))
    return rows


class BBoxOutOfRange(RuntimeError):
    """Raised when the slab mesh footprint exceeds the sidecar bbox."""


# ---------------------------------------------------------------------
# Sidecar-grid-point sampling (no interpolation in the truth)
# ---------------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class GridPointSamples:
    """Per-sidecar-grid-point evaluation of the projection on a slab.

    Truth values are bit-identical to the raw sidecar voxel values
    ``Vs[i,j,k]`` -- no trilinear interpolation.  Projection values are
    the H1 finite-element evaluation at the sidecar voxel's 3-D
    coordinate, computed on the enclosing mesh tet via barycentric
    coords (P1) or the 10-DOF Lagrange basis (P2).
    """
    xyz:    np.ndarray   # (N, 3) UTM-11N coords of every sampled voxel
    proj:   np.ndarray   # (N,)   projection eval at xyz
    truth:  np.ndarray   # (N,)   raw sidecar value (== Vs[i,j,k])
    res:    np.ndarray   # (N,)   proj - truth
    band:   np.ndarray   # (N,)   index into the supplied bands tuple


def _locate_tets(query_xyz: np.ndarray,
                 pts: np.ndarray,
                 corner_conn: np.ndarray,
                 *,
                 k_candidates: int = 24,
                 tol: float = 1.0e-6
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Locate the enclosing tet for every query point.

    Returns ``(tet_idx, bary)`` where ``tet_idx[n] >= 0`` is the index
    of the enclosing tet (or -1 if the point is outside the mesh) and
    ``bary[n]`` is the (4,) barycentric coordinate of the query in
    that tet.  ``corner_conn`` is the (n_tets, 4) corner-vertex
    connectivity (the first 4 columns of any P1/P2 tet conn).

    Algorithm: cKDTree on tet centroids; for each query, evaluate the
    k_candidates nearest tets and pick the first whose barycentric
    coords are all >= -tol and sum to 1 +- tol.  Geometry primitive
    matches plot_comparison_with_dg0.py:render_h1_at_z_via_tet_slice.
    """
    n_q = query_xyz.shape[0]
    n_tets = corner_conn.shape[0]
    if n_q == 0 or n_tets == 0:
        return (np.full(n_q, -1, dtype=np.int64),
                np.full((n_q, 4), np.nan, dtype=np.float64))

    V_xyz = pts[corner_conn]                      # (n_tets, 4, 3)
    centroids = V_xyz.mean(axis=1)                # (n_tets, 3)
    tree = cKDTree(centroids)
    k = min(k_candidates, n_tets)
    _, candidates = tree.query(query_xyz, k=k)    # (n_q, k)
    if candidates.ndim == 1:
        candidates = candidates[:, None]

    # Vectorised barycentric for each query against each candidate tet.
    # We loop over k since each iteration is a clean (n_q, ...) op,
    # and k is small (24).
    tet_idx_out = np.full(n_q, -1, dtype=np.int64)
    bary_out = np.full((n_q, 4), np.nan, dtype=np.float64)
    found = np.zeros(n_q, dtype=bool)

    for kk in range(k):
        unresolved = ~found
        if not unresolved.any():
            break
        cand = candidates[unresolved, kk]                       # (n_u,)
        Vc = V_xyz[cand]                                        # (n_u,4,3)
        A = np.stack([Vc[:, 1] - Vc[:, 0],
                      Vc[:, 2] - Vc[:, 0],
                      Vc[:, 3] - Vc[:, 0]], axis=2)             # (n_u,3,3)
        rhs = (query_xyz[unresolved] - Vc[:, 0])[:, :, None]    # (n_u,3,1)
        try:
            bc = np.linalg.solve(A, rhs).squeeze(-1)             # (n_u,3)
        except np.linalg.LinAlgError:
            # Singular tets in candidate list -- skip this candidate.
            continue
        b0 = 1.0 - bc.sum(axis=1)
        bary = np.column_stack([b0, bc])                         # (n_u,4)
        ok = ((bary >= -tol).all(axis=1) &
              (np.abs(bary.sum(axis=1) - 1.0) < tol))
        # Map results back to global query indices.
        unresolved_idx = np.flatnonzero(unresolved)
        hits = unresolved_idx[ok]
        tet_idx_out[hits] = cand[ok]
        bary_out[hits] = bary[ok]
        found[hits] = True
    return tet_idx_out, bary_out


def per_grid_point_residuals(
    vtu_path: Path,
    sidecar_path: Path,
    stl_path: Path,
    *,
    field: str = "Vs",
    z_target: float = -1000.0,
    half_thick: float = 200.0,
    bands: Iterable[tuple[float, float]] = DEFAULT_BANDS,
    parse_mfem_vtu: Path | None = None,
) -> GridPointSamples:
    """Sample the projected GF at every sidecar grid voxel that lies
    inside the mesh on the slab ``z = z_target +- half_thick``.

    The "truth" at each sample is the **raw sidecar voxel value**
    ``F[i,j,k]`` -- no trilinear interpolation.  This is the
    bit-identical comparison the user asked for: every truth datum is
    a literal entry in the .h5.

    Voxels outside the mesh (or inside it but with ill-conditioned
    enclosing tet locator) are dropped from the returned arrays.
    """
    parse_vtu = _load_parse_vtu(parse_mfem_vtu)
    bands_t = tuple(bands)

    vtu_path = Path(vtu_path)
    sidecar_path = Path(sidecar_path)
    stl_path = Path(stl_path)

    with h5py.File(sidecar_path, "r") as h5:
        if "fields" not in h5 or field not in h5["fields"]:
            raise KeyError(
                f"sidecar has no field '{field}'")
        gx = np.asarray(h5["grid/x"][...], dtype=np.float64)
        gy = np.asarray(h5["grid/y"][...], dtype=np.float64)
        gz = np.asarray(h5["grid/z"][...], dtype=np.float64)
        F = np.asarray(h5[f"fields/{field}"][...], dtype=np.float64)

    # Z-slab voxel selection: every sidecar k whose gz[k] is in band.
    k_in_slab = np.flatnonzero(
        (gz >= z_target - half_thick) & (gz <= z_target + half_thick))
    if k_in_slab.size == 0:
        return GridPointSamples(
            xyz=np.empty((0, 3)), proj=np.empty(0),
            truth=np.empty(0), res=np.empty(0),
            band=np.empty(0, dtype=np.int64))

    m = parse_vtu(str(vtu_path))
    pts = np.ascontiguousarray(m["points"], dtype=np.float64)
    conn_flat = np.asarray(m["conn"], dtype=np.int64)
    offsets = np.asarray(m["offsets"], dtype=np.int64)
    dofs, conn_2d = _split_conn_by_offset(conn_flat, offsets)
    pd = m["point_data"]
    if field not in pd:
        raise KeyError(f"vtu has no PointData '{field}'")
    vert_vals = np.asarray(pd[field], dtype=np.float64)

    corner_conn = conn_2d[:, :4]

    # Build (X, Y, Z, truth) records for every (i, j, k) with k in slab.
    II, JJ = np.meshgrid(np.arange(gx.size), np.arange(gy.size),
                          indexing="ij")
    II = II.ravel(); JJ = JJ.ravel()
    Nxy = II.size

    all_xyz: list[np.ndarray] = []
    all_truth: list[np.ndarray] = []
    for k in k_in_slab:
        x = gx[II]
        y = gy[JJ]
        z = np.full(Nxy, gz[k])
        truth = F[II, JJ, k]
        all_xyz.append(np.column_stack([x, y, z]))
        all_truth.append(truth)
    query_xyz = np.vstack(all_xyz)
    truth_all = np.concatenate(all_truth)

    # Locate enclosing tets.
    tet_idx, bary = _locate_tets(query_xyz, pts, corner_conn)
    inside = tet_idx >= 0
    if not inside.any():
        return GridPointSamples(
            xyz=np.empty((0, 3)), proj=np.empty(0),
            truth=np.empty(0), res=np.empty(0),
            band=np.empty(0, dtype=np.int64))

    query_xyz = query_xyz[inside]
    truth_all = truth_all[inside]
    tet_idx = tet_idx[inside]
    bary = bary[inside]

    # Evaluate the projection at the query barycentric coords.
    if dofs == 4:
        F_corner = vert_vals[corner_conn[tet_idx]]            # (n,4)
        proj = (bary * F_corner).sum(axis=1)
    elif dofs == 10:
        # P2 closed-form: phi_corner_i = l_i (2 l_i - 1);
        # phi_edge_(i,j) = 4 l_i l_j.
        n = tet_idx.size
        phi = np.empty((n, 10), dtype=np.float64)
        for i in range(4):
            phi[:, i] = bary[:, i] * (2.0 * bary[:, i] - 1.0)
        for k_e, (i, j) in enumerate(_P2_EDGES):
            phi[:, 4 + k_e] = 4.0 * bary[:, i] * bary[:, j]
        node_vals = vert_vals[conn_2d[tet_idx]]               # (n,10)
        proj = (phi * node_vals).sum(axis=1)
    else:
        raise ValueError(f"unsupported dofs/elem={dofs}")

    res = proj - truth_all

    # Distance to fault for band binning.
    stl_v = load_stl_vertices(str(stl_path))
    if stl_v is None or len(stl_v) == 0:
        raise ValueError(f"STL '{stl_path}' parsed empty")
    stl_xy = stl_v.mean(axis=1)[:, :2]
    tree = cKDTree(stl_xy)
    d, _ = tree.query(query_xyz[:, :2], k=1)
    band_idx = np.full(query_xyz.shape[0], -1, dtype=np.int64)
    for bi, (lo, hi) in enumerate(bands_t):
        mask = (d >= lo) & (d < hi)
        band_idx[mask] = bi

    return GridPointSamples(
        xyz=query_xyz, proj=proj, truth=truth_all,
        res=res, band=band_idx)


def grid_point_band_rms(samples: GridPointSamples,
                        bands: Iterable[tuple[float, float]]
                        ) -> list[BandRow]:
    """Per-band RMS of (proj - raw_sidecar) at sidecar grid points."""
    rows: list[BandRow] = []
    for bi, (lo, hi) in enumerate(tuple(bands)):
        sel = samples.band == bi
        n = int(sel.sum())
        if n == 0:
            rows.append(BandRow(lo=float(lo), hi=float(hi),
                                n_tets=0, l2_proj=float("nan")))
            continue
        rms = float(np.sqrt(np.mean(samples.res[sel] ** 2)))
        rows.append(BandRow(lo=float(lo), hi=float(hi),
                            n_tets=n, l2_proj=rms))
    return rows


# ---------------------------------------------------------------------
# Metric A: sidecar voxel midpoints (no vertex-grid alignment bonus)
# ---------------------------------------------------------------------

def per_voxel_midpoint_residuals(
    vtu_path: Path,
    sidecar_path: Path,
    stl_path: Path,
    *,
    field: str = "Vs",
    z_target: float = -1000.0,
    half_thick: float = 200.0,
    bands: Iterable[tuple[float, float]] = DEFAULT_BANDS,
    parse_mfem_vtu: Path | None = None,
) -> GridPointSamples:
    """Sample the projected GF at every sidecar voxel **midpoint**
    (instead of corner) on the slab z = z_target +/- half_thick.

    For an axis with grid x_0 < x_1 < ... < x_{Nx-1}, the midpoints
    are ((x_i + x_{i+1}) / 2) for i = 0 .. Nx-2.  At each midpoint,
    truth = trilinear evaluation of the sidecar at that midpoint
    (this is interpolation **within the sidecar**, which the user has
    permitted), and proj = H1 in-tet evaluation at the same midpoint.

    Why this metric: ``per_grid_point_residuals`` rewards meshes whose
    vertices coincidentally align with the sidecar grid (e.g.,
    BASELINE's lc_near=1500m matches sidecar dx=1500m, so many BASELINE
    vertices land on sidecar voxel corners and the metric reads near
    zero there).  Voxel midpoints break this alignment so the metric
    measures the projection's off-grid accuracy — what the wave
    operator actually sees at quadrature points.
    """
    parse_vtu = _load_parse_vtu(parse_mfem_vtu)
    bands_t = tuple(bands)

    vtu_path = Path(vtu_path)
    sidecar_path = Path(sidecar_path)
    stl_path = Path(stl_path)

    with h5py.File(sidecar_path, "r") as h5:
        if "fields" not in h5 or field not in h5["fields"]:
            raise KeyError(f"sidecar has no field '{field}'")
        gx = np.asarray(h5["grid/x"][...], dtype=np.float64)
        gy = np.asarray(h5["grid/y"][...], dtype=np.float64)
        gz = np.asarray(h5["grid/z"][...], dtype=np.float64)
        F = np.asarray(h5[f"fields/{field}"][...], dtype=np.float64)

    # Voxel midpoints along each axis: midpoint i = (axis[i] +
    # axis[i+1]) / 2 for i = 0 .. N-2.  These are NEVER on a sidecar
    # grid corner so the projection has to actually interpolate
    # in-tet to evaluate them.
    mx = 0.5 * (gx[:-1] + gx[1:])
    my = 0.5 * (gy[:-1] + gy[1:])
    mz = 0.5 * (gz[:-1] + gz[1:])

    # Z-slab voxel-midpoint selection.
    k_in_slab = np.flatnonzero(
        (mz >= z_target - half_thick) & (mz <= z_target + half_thick))
    if k_in_slab.size == 0:
        return GridPointSamples(
            xyz=np.empty((0, 3)), proj=np.empty(0),
            truth=np.empty(0), res=np.empty(0),
            band=np.empty(0, dtype=np.int64))

    m = parse_vtu(str(vtu_path))
    pts = np.ascontiguousarray(m["points"], dtype=np.float64)
    conn_flat = np.asarray(m["conn"], dtype=np.int64)
    offsets = np.asarray(m["offsets"], dtype=np.int64)
    dofs, conn_2d = _split_conn_by_offset(conn_flat, offsets)
    pd = m["point_data"]
    if field not in pd:
        raise KeyError(f"vtu has no PointData '{field}'")
    vert_vals = np.asarray(pd[field], dtype=np.float64)

    corner_conn = conn_2d[:, :4]

    # Build (X, Y, Z) records at every (i, j, k) midpoint with
    # k in slab.  The (i, j) midpoints span all of mx x my.
    II, JJ = np.meshgrid(np.arange(mx.size), np.arange(my.size),
                          indexing="ij")
    II = II.ravel(); JJ = JJ.ravel()
    Nxy = II.size

    all_xyz: list[np.ndarray] = []
    for k in k_in_slab:
        x = mx[II]
        y = my[JJ]
        z = np.full(Nxy, mz[k])
        all_xyz.append(np.column_stack([x, y, z]))
    query_xyz = np.vstack(all_xyz)

    # Truth at midpoint: trilinear interpolation of F at the midpoint.
    # The trilinear evaluator from plot_comparison_with_dg0 needs the
    # full F array + axis coords; reuse it.
    from plot_comparison_with_dg0 import trilinear                  # noqa: E402
    truth_all = trilinear(gx, gy, gz, F, query_xyz)

    # Locate enclosing tets.
    tet_idx, bary = _locate_tets(query_xyz, pts, corner_conn)
    inside = tet_idx >= 0
    if not inside.any():
        return GridPointSamples(
            xyz=np.empty((0, 3)), proj=np.empty(0),
            truth=np.empty(0), res=np.empty(0),
            band=np.empty(0, dtype=np.int64))

    query_xyz = query_xyz[inside]
    truth_all = truth_all[inside]
    tet_idx = tet_idx[inside]
    bary = bary[inside]

    # Evaluate the projection at the query barycentric coords.
    if dofs == 4:
        F_corner = vert_vals[corner_conn[tet_idx]]
        proj = (bary * F_corner).sum(axis=1)
    elif dofs == 10:
        n = tet_idx.size
        phi = np.empty((n, 10), dtype=np.float64)
        for i in range(4):
            phi[:, i] = bary[:, i] * (2.0 * bary[:, i] - 1.0)
        for k_e, (i, j) in enumerate(_P2_EDGES):
            phi[:, 4 + k_e] = 4.0 * bary[:, i] * bary[:, j]
        node_vals = vert_vals[conn_2d[tet_idx]]
        proj = (phi * node_vals).sum(axis=1)
    else:
        raise ValueError(f"unsupported dofs/elem={dofs}")

    res = proj - truth_all

    # Distance-to-fault binning, same as per_grid_point_residuals.
    stl_v = load_stl_vertices(str(stl_path))
    if stl_v is None or len(stl_v) == 0:
        raise ValueError(f"STL '{stl_path}' parsed empty")
    stl_xy = stl_v.mean(axis=1)[:, :2]
    tree = cKDTree(stl_xy)
    d, _ = tree.query(query_xyz[:, :2], k=1)
    band_idx = np.full(query_xyz.shape[0], -1, dtype=np.int64)
    for bi, (lo, hi) in enumerate(bands_t):
        mask = (d >= lo) & (d < hi)
        band_idx[mask] = bi

    return GridPointSamples(
        xyz=query_xyz, proj=proj, truth=truth_all,
        res=res, band=band_idx)


# ---------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------

def _format_band_label(lo: float, hi: float) -> str:
    hi_str = "inf" if math.isinf(hi) else str(int(hi))
    return f"[{int(lo)}-{hi_str}]"


def _print_table(rows: list[BandRow]) -> None:
    """Match the layout of ``plot_comparison_with_dg0.py:469-486`` on
    the H1 column (which is the only projection this CLI evaluates)."""
    print("per-tet integrated L2 |proj - source|  (m/s):")
    print(f"    {'band (m)':<14} {'n_tets':>6} {'proj':>9}")
    n_all = 0
    sumsq_all = 0.0
    for r in rows:
        print(f"    {_format_band_label(r.lo, r.hi):<14} "
              f"{r.n_tets:>6} "
              f"{(f'{r.l2_proj:>9.1f}' if not math.isnan(r.l2_proj) else '      nan')}")
        if r.n_tets > 0 and not math.isnan(r.l2_proj):
            n_all += r.n_tets
            sumsq_all += r.n_tets * r.l2_proj ** 2
    if n_all > 0:
        global_l2 = math.sqrt(sumsq_all / n_all)
        print(f"    {'GLOBAL':<14} {n_all:>6} {global_l2:>9.1f}")
    else:
        print(f"    {'GLOBAL':<14} {0:>6}       nan")


def _write_json(rows: list[BandRow], json_path: Path,
                *, vtu_path: Path, sidecar_path: Path,
                stl_path: Path, field: str,
                z_target: float, half_thick: float) -> None:
    bands_records = [
        {"lo_m": r.lo,
         "hi_m": (None if math.isinf(r.hi) else r.hi),
         "n_tets": r.n_tets,
         "l2": (None if math.isnan(r.l2_proj) else r.l2_proj)}
        for r in rows
    ]
    n_all = sum(r.n_tets for r in rows
                if r.n_tets > 0 and not math.isnan(r.l2_proj))
    if n_all > 0:
        sumsq_all = sum(r.n_tets * r.l2_proj ** 2
                        for r in rows
                        if r.n_tets > 0 and not math.isnan(r.l2_proj))
        global_l2 = math.sqrt(sumsq_all / n_all)
    else:
        global_l2 = None
    record = {
        "vtu":          str(vtu_path),
        "sidecar":      str(sidecar_path),
        "stl":          str(stl_path),
        "field":        field,
        "z_target_m":   z_target,
        "half_thick_m": half_thick,
        "bands":        bands_records,
        "global":       {"n_tets": n_all, "l2": global_l2},
    }
    Path(json_path).write_text(json.dumps(record, indent=2))


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def _parse_bands(spec: str) -> tuple[tuple[float, float], ...]:
    """Parse '0,500,1500,3000,inf' → ((0,500),(500,1500),...)."""
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    if len(parts) < 2:
        raise argparse.ArgumentTypeError(
            f"--bands needs at least two edges; got '{spec}'")
    edges: list[float] = []
    for p in parts:
        if p.lower() in ("inf", "+inf", "infinity"):
            edges.append(math.inf)
        else:
            edges.append(float(p))
    out = []
    for i in range(len(edges) - 1):
        if edges[i] >= edges[i + 1]:
            raise argparse.ArgumentTypeError(
                f"--bands edges must be strictly increasing; got "
                f"{edges[i]} >= {edges[i+1]} at index {i}")
        out.append((edges[i], edges[i + 1]))
    return tuple(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vtu", type=Path, required=True)
    ap.add_argument("--sidecar", type=Path, required=True)
    ap.add_argument("--stl", type=Path, required=True)
    ap.add_argument("--field", default="Vs")
    ap.add_argument("--z-target", type=float, default=-1000.0,
                    dest="z_target")
    ap.add_argument("--half-thick", type=float, default=200.0,
                    dest="half_thick")
    ap.add_argument("--bands", type=_parse_bands,
                    default=DEFAULT_BANDS,
                    help="comma-separated list of band edges in m; "
                         "trailing 'inf' allowed (default: "
                         "0,500,1500,3000,inf)")
    ap.add_argument("--json", type=Path, default=None,
                    dest="json_path")
    ap.add_argument("--parse-mfem-vtu", type=Path, default=None,
                    dest="parse_mfem_vtu")
    args = ap.parse_args(argv)

    try:
        rows = per_tet_l2_by_band(
            args.vtu, args.sidecar, args.stl,
            field=args.field,
            z_target=args.z_target,
            half_thick=args.half_thick,
            bands=args.bands,
            parse_mfem_vtu=args.parse_mfem_vtu)
    except BBoxOutOfRange:
        return 2
    except (FileNotFoundError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    _print_table(rows)
    if args.json_path is not None:
        _write_json(rows, args.json_path,
                    vtu_path=args.vtu, sidecar_path=args.sidecar,
                    stl_path=args.stl, field=args.field,
                    z_target=args.z_target,
                    half_thick=args.half_thick)
    return 0


if __name__ == "__main__":
    sys.exit(main())
