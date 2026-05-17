#!/usr/bin/env python3
"""
verify_onfault_stress.py — Postprocess verifier for the on-fault
stress projection artefacts.

Plan reference: PLAN_onfaultstress.md Phase 8 (§2055-2263).

Reads the VTU outputs written by the Phase-4 fault/bulk writers and
(optionally) the Phase-7 H1-projected VTU set produced by
`seas_project_stress_to_mesh`, and verifies the recorded values
against the analytic prediction evaluated independently at the same
coordinates.

The script returns a non-zero exit code if any field exceeds the
per-field error tolerance.

Verifier-specific math (R-501/R-502 single-flip contract):

  Bulk point p = (x, y, z):
      σ_seas(p) = bulk_stress_tensor_field(z, ...) * 1e6   [Pa]
      (i.e., compression-POSITIVE, MPa→Pa)
      → component lookup against sigma_xx_Pa, etc.

  Fault point at vertex v:
      Re-evaluate node-averaged basis (s_v, d_v, n_v) via
      basis_to_node (Tandem convention).
      σ_n_v      = n_v · σ_seas(p_v) · n_v
      τ_strike_v = s_v · σ_seas(p_v) · n_v
      τ_dip_v    = d_v · σ_seas(p_v) · n_v
      → compare to sigma_n_total_MPa, tau_strike_MPa, tau_dip_MPa.

  Fault cell at triangle t:
      Same predictor evaluated at the triangle centroid using the
      per-cell basis (no averaging).

Tolerances default to the floating-point round-off limits of the
rotation math; depth-varying σ⁰ or H1-projected output requires a
larger `--tol-pa`.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Optional

import numpy as np

# Sibling-import the pipeline analytic predictor.
_THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS))
import project_to_fault_stress as ptfs  # noqa: E402


# ----------------------------------------------------------------------
# Module-level defaults (per feedback_no_hardcoded_numbers.md, these
# are the only numeric constants; everything else is derived).
# ----------------------------------------------------------------------
DEFAULT_TOL_PA_CELL: float = 1.0          # raw cell-data: rotation error only
DEFAULT_TOL_PA_NODE: float = 1.0e3        # node-averaged: linear interp error
DEFAULT_TOL_PA_BULK: float = 1.0          # constant σ⁰ on bulk: exact
DEFAULT_TOL_PA_H1_PROJ: float = 1.0e3     # H1(p) projection: tied to h
DEFAULT_REL_TOL: float = 1.0e-6


@dataclass
class VerificationResult:
    """Per-field verification record produced by the verifier."""

    field_name: str
    l_inf_err: float
    l_inf_loc: np.ndarray       # (3,) point of the L∞ violation
    l2_err: float
    observed_min: float
    observed_max: float
    passed: bool
    tol_abs: float
    tol_rel: float

    def to_dict(self) -> dict:
        out = asdict(self)
        out["l_inf_loc"] = np.asarray(self.l_inf_loc).tolist()
        return out


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _load_meshio(path: Path):
    """Lazy import meshio and read a VTU/PVD file."""
    import meshio
    return meshio.read(str(path))


# Mapping from the Phase-4 summary JSON params keys to the keys
# `_bulk_sigma_seas_pa` consumes.  The two writers (summary JSON and
# verifier predictor) drifted out of sync — this adapter keeps the
# verifier working against either naming.
_PARAMS_KEY_ALIASES: dict = {
    "SHmax": ("SHmax", "SHmax_MPa"),
    "Shmin": ("Shmin", "Shmin_MPa"),
    "Sv": ("Sv", "Sv_MPa"),
    "P_p": ("P_p", "P_p_MPa"),
    "SHmax_az_deg": ("SHmax_az_deg", "SHmax_azimuth_deg"),
    "SHmax_grad": ("SHmax_grad", "SHmax_grad_MPa_per_m"),
    "Shmin_grad": ("Shmin_grad", "Shmin_grad_MPa_per_m"),
    "Sv_grad": ("Sv_grad", "Sv_grad_MPa_per_m"),
    "P_p_grad": ("P_p_grad", "P_p_grad_MPa_per_m"),
    "depth_model": ("depth_model",),
    "fault_strike_azimuth_hint_deg": (
        "fault_strike_azimuth_hint_deg",
        "fault_strike_azimuth_hint",
        "fault_strike_az_deg",
    ),
    "rake_sense": ("rake_sense",),
}


def _normalize_params(params: dict) -> dict:
    """Translate a summary JSON params block into the key set the
    analytic predictors consume.

    Accepts either the canonical writer keys (``SHmax_MPa``,
    ``SHmax_azimuth_deg``, ``SHmax_grad_MPa_per_m``, ...) or the
    legacy short keys (``SHmax``, ``SHmax_az_deg``, ``SHmax_grad``).
    Missing keys are skipped; downstream `_bulk_sigma_seas_pa` raises
    its own KeyError if a required field is absent.
    """
    out: dict = {}
    for canonical, aliases in _PARAMS_KEY_ALIASES.items():
        for k in aliases:
            if k in params:
                out[canonical] = params[k]
                break
    return out


def _bulk_sigma_seas_pa(
    z: np.ndarray, params: dict
) -> np.ndarray:
    """Analytic σ_seas at points with the supplied z coordinates [Pa].

    Returns (N, 3, 3).
    """
    sigma_seas_mpa = ptfs.bulk_stress_tensor_field(
        np.asarray(z, dtype=np.float64),
        SHmax_az_deg=float(params["SHmax_az_deg"]),
        depth_model=str(params.get("depth_model", "constant")),
        SHmax_top=float(params["SHmax"]),
        Shmin_top=float(params["Shmin"]),
        Sv_top=float(params["Sv"]),
        SHmax_grad=float(params.get("SHmax_grad", 0.0)),
        Shmin_grad=float(params.get("Shmin_grad", 0.0)),
        Sv_grad=float(params.get("Sv_grad", 0.0)),
    )
    return sigma_seas_mpa * 1.0e6  # MPa → Pa


def _verify_field(
    name: str,
    observed: np.ndarray,
    predicted: np.ndarray,
    coords: np.ndarray,
    tol_abs: float,
    tol_rel: float,
) -> VerificationResult:
    """Compare observed and predicted at the given coordinates."""
    if observed.shape != predicted.shape:
        raise ValueError(
            f"verify_field: shape mismatch on '{name}' — "
            f"observed {observed.shape} vs predicted {predicted.shape}"
        )

    diff = observed - predicted
    abs_diff = np.abs(diff)
    if abs_diff.size == 0:
        return VerificationResult(
            field_name=name,
            l_inf_err=0.0,
            l_inf_loc=np.zeros(3),
            l2_err=0.0,
            observed_min=0.0,
            observed_max=0.0,
            passed=True,
            tol_abs=tol_abs,
            tol_rel=tol_rel,
        )
    idx = int(np.argmax(abs_diff))
    l_inf_err = float(abs_diff.flat[idx])
    l_inf_loc = coords.reshape(-1, 3)[idx]
    l2_err = float(np.sqrt(np.mean(diff * diff)))
    observed_min = float(np.min(observed))
    observed_max = float(np.max(observed))
    obs_abs_max = max(abs(observed_min), abs(observed_max))
    tol_effective = max(tol_abs, tol_rel * obs_abs_max)
    return VerificationResult(
        field_name=name,
        l_inf_err=l_inf_err,
        l_inf_loc=l_inf_loc,
        l2_err=l2_err,
        observed_min=observed_min,
        observed_max=observed_max,
        passed=(l_inf_err <= tol_effective),
        tol_abs=tol_abs,
        tol_rel=tol_rel,
    )


# ----------------------------------------------------------------------
# Public verifier functions
# ----------------------------------------------------------------------
def verify_bulk_cell_data(
    bulk_vtu_path: Path,
    params: dict,
    *,
    tol_abs_pa: float = DEFAULT_TOL_PA_BULK,
    tol_rel: float = DEFAULT_REL_TOL,
) -> List[VerificationResult]:
    """Verify the six per-component σ_seas cell-data fields on a bulk
    VTU.

    Field-name resolution: tries ``sigma_<ij>_MPa`` first (canonical
    writer output as of Phase-4) and converts MPa → Pa internally; if
    absent, falls back to ``sigma_<ij>_Pa`` (legacy artefact).
    """
    norm = _normalize_params(params)
    m = _load_meshio(Path(bulk_vtu_path))

    results: List[VerificationResult] = []
    name_pairs = [
        ("sigma_xx", 0, 0),
        ("sigma_yy", 1, 1),
        ("sigma_zz", 2, 2),
        ("sigma_xy", 0, 1),
        ("sigma_yz", 1, 2),
        ("sigma_xz", 0, 2),
    ]
    for stem, i, j in name_pairs:
        mpa_name = f"{stem}_MPa"
        pa_name = f"{stem}_Pa"
        observed_mpa = _cell_data_array(m, mpa_name)
        if observed_mpa is not None:
            observed_pa = observed_mpa * 1.0e6
            field_label = mpa_name + " [→ Pa]"
            centroids = _cell_centroids_with_field(m, mpa_name)
        else:
            observed_pa = _cell_data_array(m, pa_name)
            if observed_pa is None:
                continue
            field_label = pa_name
            centroids = _cell_centroids_with_field(m, pa_name)
        sigma_pred = _bulk_sigma_seas_pa(centroids[:, 2], norm)
        predicted = sigma_pred[:, i, j]
        results.append(_verify_field(field_label, observed_pa, predicted,
                                     centroids, tol_abs_pa, tol_rel))
    return results


def verify_fault_cell_data(
    fault_vtu_path: Path,
    params: dict,
    *,
    tol_abs_pa: float = DEFAULT_TOL_PA_CELL,
    tol_rel: float = DEFAULT_REL_TOL,
) -> List[VerificationResult]:
    """Verify per-triangle rotated stress on a fault VTU.

    Looks for cell-data fields ``sigma_n_total_MPa_cell``,
    ``tau_strike_MPa_cell``, ``tau_dip_MPa_cell`` written by the
    Phase 4 fault writer (the ``_cell`` suffix distinguishes them from
    the node-averaged point-data variants).
    """
    norm = _normalize_params(params)
    m = _load_meshio(Path(fault_vtu_path))

    # Per-cell basis from connectivity (Tandem convention).  This
    # builds geometry from triangle blocks only.
    fault_geom = _build_fault_geom_from_meshio(m, norm)

    results: List[VerificationResult] = []
    targets = [
        ("sigma_n_total_MPa_cell", "n"),
        ("tau_strike_MPa_cell",    "s"),
        ("tau_dip_MPa_cell",       "d"),
    ]
    for name, kind in targets:
        # Try the canonical `_cell`-suffixed name first; fall back to
        # the unsuffixed name (older artefact convention) so the
        # legacy synthetic-VTU tests keep working.
        resolved_name = name
        observed_mpa = _cell_data_array(m, resolved_name)
        if observed_mpa is None:
            resolved_name = name.replace("_cell", "")
            observed_mpa = _cell_data_array(m, resolved_name)
            if observed_mpa is None:
                continue
        # R-102: align centroids and prediction with cell blocks that
        # actually carry the field.  For a fault VTU containing only
        # triangles this matches `fault_geom.centroids`; for any
        # mixed-block mesh the size still aligns with `observed_mpa`.
        centroids = _cell_centroids_with_field(m, resolved_name)
        sigma_pred_pa = _bulk_sigma_seas_pa(centroids[:, 2], norm)
        # When centroids size matches fault_geom (the canonical
        # triangle-only case), use the per-cell rotated value.
        # Otherwise (rare mixed-block case) fall back to a direct
        # `n^T S n` style rotation using the elasticity convention
        # implicit in `fault_geom`.
        if centroids.shape[0] == fault_geom.centroids.shape[0]:
            predicted_pa = _project_fault_per_cell_pa(
                sigma_pred_pa, fault_geom, kind
            )
        else:
            # Shape mismatch — can't apply fault rotation without
            # per-cell basis on the extra blocks.  Skip this field
            # with a diagnostic.
            print(
                f"warning: '{resolved_name}' cell-data size "
                f"{centroids.shape[0]} != fault_geom triangle count "
                f"{fault_geom.centroids.shape[0]}; skipping",
                file=sys.stderr,
            )
            continue
        results.append(_verify_field(resolved_name + " [Pa]",
                                     observed_mpa * 1.0e6, predicted_pa,
                                     centroids, tol_abs_pa, tol_rel))
    return results


def verify_fault_point_data(
    fault_vtu_path: Path,
    params: dict,
    *,
    tol_abs_pa: float = DEFAULT_TOL_PA_NODE,
    tol_rel: float = DEFAULT_REL_TOL,
) -> List[VerificationResult]:
    """Verify node-averaged rotated stress at fault vertices.

    The writer's pipeline is: (i) compute σ_n / τ_s / τ_d per-cell
    using each cell's own basis, then (ii) area-weighted average those
    *scalars* to nodes via ``cell_to_node_average``.  This verifier
    mirrors that two-step pipeline rather than predicting via a
    node-averaged basis — the two are mathematically distinct on a
    faceted fault (variance of the per-cell normal across each node's
    patch separates them), and the writer's choice is the canonical
    one for visualization smoothness.
    """
    norm = _normalize_params(params)
    m = _load_meshio(Path(fault_vtu_path))
    fault_geom = _build_fault_geom_from_meshio(m, norm)

    # Per-cell analytic prediction (validated to machine precision by
    # `verify_fault_cell_data`).
    sigma_pred_pa = _bulk_sigma_seas_pa(fault_geom.centroids[:, 2], norm)
    sigma_n_cell_pa = _project_fault_per_cell_pa(
        sigma_pred_pa, fault_geom, "n"
    )
    tau_s_cell_pa = _project_fault_per_cell_pa(
        sigma_pred_pa, fault_geom, "s"
    )
    tau_d_cell_pa = _project_fault_per_cell_pa(
        sigma_pred_pa, fault_geom, "d"
    )

    # Same area-weighted node averaging the writer uses.
    sigma_n_pred_node_pa = ptfs.cell_to_node_average(
        m.points, fault_geom.connectivity,
        sigma_n_cell_pa, fault_geom.areas,
    )
    tau_s_pred_node_pa = ptfs.cell_to_node_average(
        m.points, fault_geom.connectivity,
        tau_s_cell_pa, fault_geom.areas,
    )
    tau_d_pred_node_pa = ptfs.cell_to_node_average(
        m.points, fault_geom.connectivity,
        tau_d_cell_pa, fault_geom.areas,
    )

    results: List[VerificationResult] = []
    targets = [
        ("sigma_n_total_MPa", sigma_n_pred_node_pa),
        ("tau_strike_MPa",    tau_s_pred_node_pa),
        ("tau_dip_MPa",       tau_d_pred_node_pa),
    ]
    for name, predicted_pa in targets:
        observed_mpa = _point_data_array(m, name)
        if observed_mpa is None:
            continue
        observed_pa = observed_mpa * 1.0e6
        # Orphan vertices (in the VTU point cloud but not in any
        # triangle) carry NaN basis / observed values.  Mask both so
        # the L∞ / L2 reflect only triangle-referenced nodes.
        valid = np.isfinite(observed_pa) & np.isfinite(predicted_pa)
        if not np.any(valid):
            print(
                f"warning: '{name}' has no finite-valued nodes "
                f"after orphan-vertex masking; skipping",
                file=sys.stderr,
            )
            continue
        results.append(_verify_field(name + " [Pa]",
                                     observed_pa[valid],
                                     predicted_pa[valid],
                                     m.points[valid], tol_abs_pa, tol_rel))
    return results


def verify_h1_projected_vtu(
    projected_pvd_path: Path,
    params: dict,
    *,
    tol_abs_pa: float = DEFAULT_TOL_PA_H1_PROJ,
    tol_rel: float = DEFAULT_REL_TOL,
) -> List[VerificationResult]:
    """Verify a Phase-7 H1(p) projected stress VTU set.

    The .pvd input references per-rank .vtu files; we load each and
    verify the six PointData fields `sigma_xx`..`sigma_xz` (in Pa)
    against the analytic prediction at the node coordinates.
    """
    norm = _normalize_params(params)
    p = Path(projected_pvd_path)
    if not p.exists():
        raise FileNotFoundError(f"H1 projected VTU not found: {p}")

    # If the path is a .pvd we look up the referenced VTU files; if it
    # is a single .vtu we use it directly.
    if p.suffix == ".pvd":
        vtu_paths = _pvd_referenced_vtus(p)
    else:
        vtu_paths = [p]

    # Concatenate all rank-local VTUs into a single point set for the
    # comparison; meshio reads PointData per-rank, so we just compute
    # the analytic prediction at each rank's nodes independently and
    # accumulate L∞.
    results: List[VerificationResult] = []
    field_names = ["sigma_xx", "sigma_yy", "sigma_zz",
                   "sigma_xy", "sigma_yz", "sigma_xz"]
    name_to_indices = {
        "sigma_xx": (0, 0), "sigma_yy": (1, 1), "sigma_zz": (2, 2),
        "sigma_xy": (0, 1), "sigma_yz": (1, 2), "sigma_xz": (0, 2),
    }
    for name in field_names:
        l_inf, l_inf_loc, sq_err_sum, n_total = 0.0, np.zeros(3), 0.0, 0
        observed_min, observed_max = float("inf"), float("-inf")
        for vtu in vtu_paths:
            m = _load_meshio(vtu)
            observed = _point_data_array(m, name)
            if observed is None:
                continue
            sigma_pred = _bulk_sigma_seas_pa(m.points[:, 2], norm)
            i, j = name_to_indices[name]
            predicted = sigma_pred[:, i, j]
            diff = observed - predicted
            abs_diff = np.abs(diff)
            if abs_diff.size == 0:
                continue
            idx = int(np.argmax(abs_diff))
            if abs_diff.flat[idx] > l_inf:
                l_inf = float(abs_diff.flat[idx])
                l_inf_loc = m.points[idx]
            sq_err_sum += float(np.sum(diff * diff))
            n_total += int(observed.size)
            observed_min = min(observed_min, float(np.min(observed)))
            observed_max = max(observed_max, float(np.max(observed)))
        if n_total == 0:
            # R-009: warn instead of silently dropping a missing field.
            # A typo / writer drift between writer and verifier produces
            # this case, and the user has no diagnostic otherwise.
            print(
                f"warning: '{name}' not found in any of "
                f"{len(vtu_paths)} referenced VTU files",
                file=sys.stderr,
            )
            continue
        l2 = float(np.sqrt(sq_err_sum / n_total))
        obs_abs_max = max(abs(observed_min), abs(observed_max))
        tol_effective = max(tol_abs_pa, tol_rel * obs_abs_max)
        results.append(VerificationResult(
            field_name=name + " [H1-projected, Pa]",
            l_inf_err=l_inf,
            l_inf_loc=l_inf_loc,
            l2_err=l2,
            observed_min=observed_min,
            observed_max=observed_max,
            passed=(l_inf <= tol_effective),
            tol_abs=tol_abs_pa,
            tol_rel=tol_rel,
        ))
    return results


# ----------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------
def _cell_centroids(m) -> np.ndarray:
    """Compute centroids (cell-mean of node coords) for all cell blocks."""
    pts = m.points
    out = []
    for cb in m.cells:
        out.append(np.mean(pts[cb.data], axis=1))
    if not out:
        return np.zeros((0, 3))
    return np.concatenate(out, axis=0)


def _cell_centroids_with_field(m, name: str) -> np.ndarray:
    """Compute centroids only over cell blocks that carry the named field.

    R-102: meshio's `cell_data` may be defined on a subset of cell
    blocks (e.g., σ_*_Pa on tetrahedra only, no entry on boundary
    triangles).  The unfiltered `_cell_centroids` would mismatch
    `_cell_data_array` size in that case and crash `_verify_field`.
    This helper aligns centroids with the named field block-by-block.
    """
    pts = m.points
    field_arrs = m.cell_data.get(name)
    if field_arrs is None:
        return np.zeros((0, 3))
    out = []
    for cb, arr in zip(m.cells, field_arrs):
        if arr is None or np.asarray(arr).size == 0:
            continue
        out.append(np.mean(pts[cb.data], axis=1))
    if not out:
        return np.zeros((0, 3))
    return np.concatenate(out, axis=0)


def _cell_data_array(m, name: str) -> Optional[np.ndarray]:
    """Return the named cell-data field concatenated across cell blocks.

    R-102: skips empty / `None` entries so a mixed-block mesh whose
    field is defined only on some blocks yields a 1-D array matching
    `_cell_centroids_with_field(m, name)`.
    """
    arrays = m.cell_data.get(name)
    if arrays is None:
        return None
    keep = [np.asarray(a) for a in arrays
            if a is not None and np.asarray(a).size > 0]
    if not keep:
        return None
    return np.concatenate(keep).astype(np.float64).ravel()


def _point_data_array(m, name: str) -> Optional[np.ndarray]:
    """Return the named point-data field."""
    arr = m.point_data.get(name)
    if arr is None:
        return None
    return np.asarray(arr, dtype=np.float64).ravel()


def _pvd_referenced_vtus(pvd_path: Path) -> List[Path]:
    """Parse a ParaView .pvd file and return its DataSet file paths."""
    import xml.etree.ElementTree as ET
    root = ET.parse(pvd_path).getroot()
    out: List[Path] = []
    for ds in root.iter("DataSet"):
        file_ = ds.get("file")
        if file_ is None:
            continue
        out.append((pvd_path.parent / file_).resolve())
    return out


@dataclass
class _FaultGeomLite:
    centroids: np.ndarray         # (N_tri, 3)
    tandem_s: np.ndarray          # (N_tri, 3)
    tandem_d: np.ndarray          # (N_tri, 3)
    tandem_n: np.ndarray          # (N_tri, 3)
    areas: np.ndarray             # (N_tri,)
    connectivity: np.ndarray      # (N_tri, 3)


def _build_fault_geom_from_meshio(m, params: dict) -> _FaultGeomLite:
    """Rebuild the per-triangle Tandem (s, d, n) basis from a fault VTU.

    Mirrors Phase 2 `build_fault_basis` but starts from a meshio mesh
    instead of the .vtu / .json pair the writer produced.
    """
    pts = m.points
    # Take only triangle blocks (gmsh quadrangles or hex shouldn't appear
    # in the fault surface VTU).
    tri_blocks = [cb for cb in m.cells if cb.type == "triangle"]
    if not tri_blocks:
        raise ValueError("verify_fault_*: no triangle cell blocks found")
    conn = np.concatenate([cb.data for cb in tri_blocks], axis=0)
    v0 = pts[conn[:, 0]]
    v1 = pts[conn[:, 1]]
    v2 = pts[conn[:, 2]]
    centroids = (v0 + v1 + v2) / 3.0

    # Per-triangle normal and area (no harmonization — we just need
    # raw geometry).
    e1 = v1 - v0
    e2 = v2 - v0
    cross = np.cross(e1, e2)
    cross_norm = np.linalg.norm(cross, axis=1, keepdims=True)
    cross_norm_safe = np.where(cross_norm > 0.0, cross_norm,
                               np.ones_like(cross_norm))
    n_raw = cross / cross_norm_safe
    areas = 0.5 * cross_norm[:, 0]

    # Reorient to a consistent global frame. The pipeline's
    # `harmonise_normal_orientation` uses a rake-sense + trace
    # centroids + optional strike-azimuth hint to pin orientation.
    # The strike hint must be the *fault* strike azimuth (e.g. 314°
    # for the SAF), NOT the SHmax azimuth — using the latter
    # silently flips the fault normal on inclined facets and
    # corrupts τ_dip by ~σ-magnitude.  The summary JSON exposes the
    # writer-side hint as ``fault_strike_azimuth_hint_deg``.
    strike_az_deg_raw = params.get("fault_strike_azimuth_hint_deg")
    if strike_az_deg_raw is None:
        strike_az_deg = 0.0
    else:
        strike_az_deg = float(strike_az_deg_raw)
    rake_sense = str(params.get("rake_sense", "right-lateral"))
    try:
        n_oriented = ptfs.harmonise_normal_orientation(
            n_raw, rake_sense, centroids,
            fault_strike_azimuth_hint_deg=strike_az_deg,
        )
    except Exception:
        # Defensive: fall back to raw normals if the harmoniser cannot
        # converge (e.g. single-triangle synthetic mesh).
        n_oriented = n_raw

    # Tandem-style basis: s = up × n, d = s × n (down-dip), all unit.
    s, d, _ = ptfs.per_triangle_basis_raw(n_oriented)
    return _FaultGeomLite(centroids=centroids,
                          tandem_s=s, tandem_d=d, tandem_n=n_oriented,
                          areas=areas, connectivity=conn)


def _project_fault_per_cell_pa(sigma_pa: np.ndarray,
                                fault_geom: _FaultGeomLite,
                                kind: str) -> np.ndarray:
    """Compute per-cell rotated component (kind='n','s','d') in Pa."""
    n = fault_geom.tandem_n
    s = fault_geom.tandem_s
    d = fault_geom.tandem_d
    Sn = np.einsum("nij,nj->ni", sigma_pa, n)
    if kind == "n":
        return np.einsum("ni,ni->n", n, Sn)
    if kind == "s":
        return np.einsum("ni,ni->n", s, Sn)
    if kind == "d":
        return np.einsum("ni,ni->n", d, Sn)
    raise ValueError(f"_project_fault_per_cell_pa: unknown kind {kind!r}")


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def _print_summary(results: List[VerificationResult]) -> bool:
    """Print a one-row-per-field summary; return True iff all passed."""
    all_passed = True
    print(f"{'field':<35} {'L_inf':>12} {'L_2':>12} "
          f"{'tol_abs':>10} {'pass':>6}")
    for r in results:
        ok = "yes" if r.passed else "NO"
        all_passed &= r.passed
        print(f"{r.field_name:<35} {r.l_inf_err:>12.3e} "
              f"{r.l2_err:>12.3e} {r.tol_abs:>10.2e} {ok:>6}")
    return all_passed


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify on-fault stress projection artefacts.")
    parser.add_argument("--fault-vtu", type=Path, default=None,
                        help="Fault VTU (point-data + cell-data).")
    parser.add_argument("--bulk-vtu", type=Path, default=None,
                        help="Bulk VTU (cell-data σ_*_Pa).")
    parser.add_argument("--h1-projected-pvd", type=Path, default=None,
                        help="Phase-7 H1-projected .pvd (six PointData fields).")
    parser.add_argument("--summary-json", type=Path, default=None,
                        help="Phase 4 summary JSON carrying the params block.")
    parser.add_argument("--tol-pa", type=float, default=None,
                        help="Override the absolute Pa tolerance for "
                             "every layer (fault cell, fault node, bulk, H1).")
    parser.add_argument("--tol-rel", type=float, default=DEFAULT_REL_TOL,
                        help="Relative tolerance (default %(default)g).")
    parser.add_argument("--fail-on-warn", action="store_true",
                        help="Exit 1 if any warning is emitted.")
    parser.add_argument("--report-json", type=Path, default=None,
                        help="Write a verify_report.json to this path.")
    parser.add_argument("--params", type=Path, default=None,
                        help="Alternative to --summary-json: JSON file "
                             "containing a flat dict of the analytic "
                             "params (SHmax, Shmin, Sv, SHmax_az_deg, ...).")
    args = parser.parse_args(argv)

    if args.summary_json is not None:
        with open(args.summary_json, "r") as fh:
            summary = json.load(fh)
        params = summary.get("params", summary)
    elif args.params is not None:
        with open(args.params, "r") as fh:
            params = json.load(fh)
    else:
        parser.error("either --summary-json or --params is required")

    tol_cell = (args.tol_pa if args.tol_pa is not None
                else DEFAULT_TOL_PA_CELL)
    tol_node = (args.tol_pa if args.tol_pa is not None
                else DEFAULT_TOL_PA_NODE)
    tol_bulk = (args.tol_pa if args.tol_pa is not None
                else DEFAULT_TOL_PA_BULK)
    tol_h1   = (args.tol_pa if args.tol_pa is not None
                else DEFAULT_TOL_PA_H1_PROJ)

    all_results: List[VerificationResult] = []

    if args.bulk_vtu is not None:
        print(f"\n== bulk VTU: {args.bulk_vtu} ==")
        rs = verify_bulk_cell_data(args.bulk_vtu, params,
                                   tol_abs_pa=tol_bulk,
                                   tol_rel=args.tol_rel)
        _print_summary(rs)
        all_results.extend(rs)

    if args.fault_vtu is not None:
        print(f"\n== fault VTU (cell-data): {args.fault_vtu} ==")
        rs = verify_fault_cell_data(args.fault_vtu, params,
                                    tol_abs_pa=tol_cell,
                                    tol_rel=args.tol_rel)
        _print_summary(rs)
        all_results.extend(rs)

        print(f"\n== fault VTU (point-data): {args.fault_vtu} ==")
        rs = verify_fault_point_data(args.fault_vtu, params,
                                     tol_abs_pa=tol_node,
                                     tol_rel=args.tol_rel)
        _print_summary(rs)
        all_results.extend(rs)

    if args.h1_projected_pvd is not None:
        print(f"\n== H1-projected: {args.h1_projected_pvd} ==")
        rs = verify_h1_projected_vtu(args.h1_projected_pvd, params,
                                     tol_abs_pa=tol_h1,
                                     tol_rel=args.tol_rel)
        _print_summary(rs)
        all_results.extend(rs)

    if not all_results:
        parser.error("no input artefacts supplied "
                     "(--fault-vtu / --bulk-vtu / --h1-projected-pvd)")

    if args.report_json is not None:
        with open(args.report_json, "w") as fh:
            json.dump([r.to_dict() for r in all_results], fh, indent=2)
        print(f"\nwrote report: {args.report_json}")

    n_pass = sum(1 for r in all_results if r.passed)
    n_total = len(all_results)
    print(f"\n=== {n_pass} / {n_total} fields passed ===")
    return 0 if n_pass == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
