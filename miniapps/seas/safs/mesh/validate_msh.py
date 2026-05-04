"""Post-mesh validation for SAFS .msh files.

Implements PLAN.md Phase 6 + PLAN_smoke_millcreek.md Phase 3.

Runs a numbered suite of checks (1..10).  In single-fault mode (when
--expected-faults has length 1) the suite is the smoke-test invariants
defined in PLAN_smoke_millcreek.md.  Each check pass/fails independently;
the script exits 0 only if every check passes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from safs_origin import (
    FAULT_SHORT_NAMES,
    ORIGIN_UTM,
    local_to_utm,
    read_transform_json,
)

EXPECTED_TAGS_TRIANGLE = {1, 2, 3, 4, 5, 6, 100}
EXPECTED_TAGS_TETRA = {10}
FAULT_TAG = 100
ZTOP_TAG = 5
ZBOT_TAG = 6


@dataclass
class CheckResult:
    name: str
    passed: bool
    message: str = ""
    metrics: dict = field(default_factory=dict)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _gmsh_version() -> str:
    try:
        out = subprocess.run(["gmsh", "-version"], check=True,
                              capture_output=True, text=True)
        return (out.stdout + out.stderr).strip().splitlines()[0]
    except Exception:
        return "unknown"


def _git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True,
            capture_output=True, text=True,
            cwd=Path(__file__).resolve().parent,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def _mesh_data(msh_path: Path):
    """Return a flat description of the mesh for the checks below."""
    import meshio
    m = meshio.read(msh_path)
    points = np.asarray(m.points, dtype=np.float64)
    gp = m.cell_data_dict.get("gmsh:physical", {})

    # Triangles
    tri_blocks = [(i, cb) for i, cb in enumerate(m.cells)
                  if cb.type == "triangle"]
    tri_data = []
    tri_tags = []
    if tri_blocks:
        cursor = 0
        all_tags = gp.get("triangle")
        for _, cb in tri_blocks:
            n = len(cb.data)
            tri_data.append(np.asarray(cb.data, dtype=np.int64))
            if all_tags is not None:
                tri_tags.append(np.asarray(all_tags[cursor:cursor + n],
                                            dtype=np.int64))
                cursor += n
            else:
                tri_tags.append(np.full(n, -1, dtype=np.int64))
        tris = np.concatenate(tri_data, axis=0)
        ttags = np.concatenate(tri_tags, axis=0)
    else:
        tris = np.empty((0, 3), dtype=np.int64)
        ttags = np.empty(0, dtype=np.int64)

    # Tetrahedra
    tet_blocks = [(i, cb) for i, cb in enumerate(m.cells)
                  if cb.type == "tetra"]
    tet_data = []
    tet_tags = []
    if tet_blocks:
        cursor = 0
        all_tags = gp.get("tetra")
        for _, cb in tet_blocks:
            n = len(cb.data)
            tet_data.append(np.asarray(cb.data, dtype=np.int64))
            if all_tags is not None:
                tet_tags.append(np.asarray(all_tags[cursor:cursor + n],
                                            dtype=np.int64))
                cursor += n
            else:
                tet_tags.append(np.full(n, -1, dtype=np.int64))
        tets = np.concatenate(tet_data, axis=0)
        tetags = np.concatenate(tet_tags, axis=0)
    else:
        tets = np.empty((0, 4), dtype=np.int64)
        tetags = np.empty(0, dtype=np.int64)

    return points, tris, ttags, tets, tetags


def _build_face_tet_adjacency(tets: np.ndarray) -> dict:
    """faces (frozensets of 3 vertex indices) → list of tet indices."""
    adj = defaultdict(list)
    for i, tet in enumerate(tets):
        a, b, c, d = (int(x) for x in tet)
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            adj[frozenset(face)].append(i)
    return adj


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------
def check_1_tag_inventory(ttags, tetags) -> CheckResult:
    tri_present = set(int(t) for t in np.unique(ttags))
    tet_present = set(int(t) for t in np.unique(tetags))
    missing_tri = EXPECTED_TAGS_TRIANGLE - tri_present
    missing_tet = EXPECTED_TAGS_TETRA - tet_present
    extra_tri = tri_present - EXPECTED_TAGS_TRIANGLE
    extra_tet = tet_present - EXPECTED_TAGS_TETRA
    metrics = {
        "tri_tags": sorted(tri_present),
        "tet_tags": sorted(tet_present),
        "tri_counts": {int(t): int(c) for t, c in
                       zip(*np.unique(ttags, return_counts=True))},
        "tet_counts": {int(t): int(c) for t, c in
                       zip(*np.unique(tetags, return_counts=True))},
    }
    if missing_tri or missing_tet or extra_tri or extra_tet:
        return CheckResult(
            "1_tag_inventory", False,
            f"missing tri tags={sorted(missing_tri)}, "
            f"extra tri tags={sorted(extra_tri)}, "
            f"missing tet tags={sorted(missing_tet)}, "
            f"extra tet tags={sorted(extra_tet)}",
            metrics,
        )
    # Also: every expected tag is non-empty
    for t in EXPECTED_TAGS_TRIANGLE:
        if metrics["tri_counts"].get(int(t), 0) == 0:
            return CheckResult(
                "1_tag_inventory", False,
                f"triangle tag {t} is empty", metrics,
            )
    for t in EXPECTED_TAGS_TETRA:
        if metrics["tet_counts"].get(int(t), 0) == 0:
            return CheckResult(
                "1_tag_inventory", False,
                f"tetra tag {t} is empty", metrics,
            )
    return CheckResult("1_tag_inventory", True,
                       "all expected tags present and non-empty", metrics)


def check_2_box_arithmetic(domain_box: dict, expected: dict
                            ) -> CheckResult:
    """Compare domain_box.json to caller's expected dict.  Equality to ±1 m."""
    metrics = {"actual": domain_box, "expected": expected}
    for k, v in expected.items():
        if abs(float(domain_box.get(k, math.nan)) - float(v)) > 1.0:
            return CheckResult(
                "2_box_arithmetic", False,
                f"key {k}: expected {v}, got {domain_box.get(k)}",
                metrics,
            )
    return CheckResult("2_box_arithmetic", True,
                       "domain box matches expected snap-to-grid output",
                       metrics)


def check_3_fault_clearance(domain_box: dict, fault_bbox: dict,
                              buf_x: float, buf_y: float, depth: float,
                              min_clearance: float = 5_000.0) -> CheckResult:
    metrics = {"min_clearance_required_m": min_clearance}
    issues = []
    for short, b in fault_bbox.items():
        c_xmin = b["xmin"] - domain_box["x_min_m"]
        c_xmax = domain_box["x_max_m"] - b["xmax"]
        c_ymin = b["ymin"] - domain_box["y_min_m"]
        c_ymax = domain_box["y_max_m"] - b["ymax"]
        c_zbot = b["zmin"] - domain_box["z_bot_m"]
        for label, c in [("xmin", c_xmin), ("xmax", c_xmax),
                          ("ymin", c_ymin), ("ymax", c_ymax),
                          ("zbot", c_zbot)]:
            if c < min_clearance:
                issues.append(f"{short}.{label}: {c:.0f} m")
    if issues:
        return CheckResult("3_fault_clearance", False,
                           "; ".join(issues), metrics)
    return CheckResult("3_fault_clearance", True,
                       "every fault is >= 5 km clear of every box face",
                       metrics)


def check_4_freesurface_trace(points, tris, ttags,
                               clearance: float
                               ) -> CheckResult:
    """At least one tag-100 triangle reaches z >= -(clearance + slack) AND
    no tag-100 vertex sits above z = 0.

    Threshold derives from the ts_to_stl.py free-surface clearance: a
    triangle within (clearance + 10 m) of z=0 is treated as part of the
    surface trace.  10 m absorbs HXT-introduced node-snap noise.
    """
    slack = 10.0
    trace_z_threshold = -(float(clearance) + slack)
    mask = ttags == FAULT_TAG
    if not np.any(mask):
        return CheckResult("4_freesurface_trace", False,
                            "no tag-100 triangles", {})
    fault_tris = tris[mask]
    fault_v_z = points[fault_tris, 2]
    near_surface = (fault_v_z > trace_z_threshold).all(axis=1)
    above_zero = (points[fault_tris, 2] > 0).any(axis=1)
    n_near = int(np.sum(near_surface))
    n_above = int(np.sum(above_zero))
    metrics = {"n_near_surface": n_near, "n_above_zero": n_above,
               "trace_z_threshold_m": trace_z_threshold,
               "clearance_m": float(clearance)}
    if n_near == 0:
        return CheckResult("4_freesurface_trace", False,
                            f"no fault triangle has all 3 vertices with "
                            f"z > {trace_z_threshold} m — fault may not "
                            f"reach the free surface (clearance="
                            f"{clearance:.0f} m)",
                            metrics)
    if n_above > 0:
        return CheckResult("4_freesurface_trace", False,
                            f"{n_above} fault triangle(s) have at least "
                            f"one vertex above z=0 — z>0 clamp failed",
                            metrics)
    return CheckResult("4_freesurface_trace", True,
                       f"{n_near} fault triangle(s) reach within "
                       f"{abs(trace_z_threshold):.0f} m of free surface; "
                       f"no z>0 leakage",
                       metrics)


def check_5_internal_interface(points, tris, ttags, tets,
                                 clearance: float) -> CheckResult:
    """Every tag-100 triangle has exactly 2 adjacent tets, except trace
    triangles whose all 3 vertices lie within `clearance + 1 m` of z = 0
    (which legitimately have 1)."""
    if tets.shape[0] == 0:
        return CheckResult("5_internal_interface", False,
                            "no tetrahedra", {})
    adj = _build_face_tet_adjacency(tets)
    mask = ttags == FAULT_TAG
    trace_zmax = -(float(clearance) - 1.0)  # near-trace boundary
    bad: list[int] = []
    n_trace = 0
    n_ok = 0
    for tri_idx in np.where(mask)[0]:
        a, b, c = (int(x) for x in tris[tri_idx])
        n_tets = len(adj.get(frozenset((a, b, c)), []))
        if n_tets == 2:
            n_ok += 1
        elif n_tets == 1:
            zs = points[[a, b, c], 2]
            if np.all(zs >= trace_zmax):
                n_trace += 1
            else:
                bad.append(tri_idx)
        else:
            bad.append(tri_idx)
    metrics = {"n_internal_with_2_tets": n_ok,
               "n_trace_with_1_tet": n_trace,
               "n_offending": len(bad),
               "clearance_m": float(clearance)}
    if bad:
        return CheckResult("5_internal_interface", False,
                           f"{len(bad)} fault triangle(s) have wrong tet "
                           f"adjacency (first 5 indices: {bad[:5]}); "
                           f"fault may have collapsed to a boundary",
                           metrics)
    return CheckResult("5_internal_interface", True,
                       f"{n_ok} internal + {n_trace} trace fault triangles "
                       f"have correct tet adjacency",
                       metrics)


def _edge_lengths_of_triangles(points, tris):
    a = points[tris[:, 0]]
    b = points[tris[:, 1]]
    c = points[tris[:, 2]]
    e1 = np.linalg.norm(b - a, axis=1)
    e2 = np.linalg.norm(c - b, axis=1)
    e3 = np.linalg.norm(a - c, axis=1)
    return np.concatenate([e1, e2, e3])


def check_6_fault_edge_length(points, tris, ttags, res_f: float
                                ) -> CheckResult:
    mask = ttags == FAULT_TAG
    fault_tris = tris[mask]
    if fault_tris.shape[0] == 0:
        return CheckResult("6_fault_edge_length", False,
                            "no fault triangles", {})
    edges = _edge_lengths_of_triangles(points, fault_tris)
    mean_e = float(edges.mean())
    max_e = float(edges.max())
    # Tolerance: mean ∈ [0.5*res_f, 2*res_f], max < 4*res_f.
    metrics = {"mean_m": mean_e, "max_m": max_e,
               "res_f_m": res_f,
               "p50_m": float(np.percentile(edges, 50)),
               "p95_m": float(np.percentile(edges, 95))}
    if not (0.5 * res_f <= mean_e <= 2.0 * res_f):
        return CheckResult("6_fault_edge_length", False,
                           f"mean fault edge {mean_e:.0f} m outside "
                           f"[{0.5*res_f:.0f}, {2*res_f:.0f}]",
                           metrics)
    if max_e > 4.0 * res_f:
        return CheckResult("6_fault_edge_length", False,
                           f"max fault edge {max_e:.0f} m exceeds "
                           f"{4*res_f:.0f}",
                           metrics)
    return CheckResult("6_fault_edge_length", True,
                       f"fault edge mean={mean_e:.0f} m, max={max_e:.0f} m",
                       metrics)


def check_7_far_field_edge(points, tets, tris, ttags,
                             res_ff: float,
                             distance_threshold: float = 30_000.0
                             ) -> CheckResult:
    """Tets with all 4 vertices > threshold from any tag-100 face have
    mean edge length > 0.6*res_ff."""
    mask = ttags == FAULT_TAG
    fault_tris = tris[mask]
    if fault_tris.shape[0] == 0:
        return CheckResult("7_far_field_edge", False,
                            "no fault triangles", {})
    fault_centroids = points[fault_tris].mean(axis=1)
    from scipy.spatial import cKDTree
    t = cKDTree(fault_centroids)
    tet_v_dists = []
    for tet in tets:
        d, _ = t.query(points[tet], k=1)
        tet_v_dists.append(d)
    tet_v_dists = np.asarray(tet_v_dists)
    far_mask = (tet_v_dists.min(axis=1) > distance_threshold)
    far_tets = tets[far_mask]
    if far_tets.shape[0] == 0:
        return CheckResult("7_far_field_edge", False,
                           f"no tets are further than "
                           f"{distance_threshold} m from a fault — "
                           f"box may be too small or threshold too high",
                           {})
    a, b, c, d = (points[far_tets[:, k]] for k in range(4))
    edges = np.concatenate([
        np.linalg.norm(b - a, axis=1),
        np.linalg.norm(c - a, axis=1),
        np.linalg.norm(d - a, axis=1),
        np.linalg.norm(c - b, axis=1),
        np.linalg.norm(d - b, axis=1),
        np.linalg.norm(d - c, axis=1),
    ])
    mean_e = float(edges.mean())
    metrics = {"n_far_tets": int(far_tets.shape[0]),
               "mean_far_edge_m": mean_e,
               "res_ff_m": res_ff,
               "distance_threshold_m": distance_threshold}
    if mean_e < 0.6 * res_ff:
        return CheckResult("7_far_field_edge", False,
                           f"far-field mean edge {mean_e:.0f} m < "
                           f"{0.6*res_ff:.0f} — threshold field not "
                           f"reaching res_ff",
                           metrics)
    return CheckResult("7_far_field_edge", True,
                       f"{far_tets.shape[0]} far tets, mean edge {mean_e:.0f} m",
                       metrics)


def check_8_provenance_partition(provenance: dict, ttags) -> CheckResult:
    n_tag100 = int(np.sum(ttags == FAULT_TAG))
    n_in_prov = int(provenance["n_tag100_triangles_total"])
    sum_per_fault = sum(e["n_triangles_in_msh"]
                         for e in provenance["faults"].values())
    metrics = {"n_tag100_in_msh": n_tag100,
               "n_in_provenance": n_in_prov,
               "sum_per_fault": sum_per_fault}
    if not (n_tag100 == n_in_prov == sum_per_fault):
        return CheckResult("8_provenance_partition", False,
                           f"counts disagree: msh={n_tag100}, "
                           f"prov.total={n_in_prov}, "
                           f"sum-per-fault={sum_per_fault}",
                           metrics)
    # Also check no duplicate indices and indices are within range
    seen: set[int] = set()
    for short, e in provenance["faults"].items():
        for idx in e["triangle_indices_in_msh"]:
            if idx in seen:
                return CheckResult("8_provenance_partition", False,
                                   f"index {idx} appears in multiple "
                                   f"faults (most recent: {short})",
                                   metrics)
            seen.add(int(idx))
    return CheckResult("8_provenance_partition", True,
                       f"{n_tag100} tag-100 triangles partitioned with "
                       f"no duplicates and no missing",
                       metrics)


def check_9_utm_round_trip(points, tris, ttags, transform: dict,
                             cfm_bbox: dict | None,
                             rng: np.random.Generator) -> CheckResult:
    mask = ttags == FAULT_TAG
    fault_v = np.unique(tris[mask].ravel())
    if fault_v.size == 0:
        return CheckResult("9_utm_round_trip", False,
                            "no fault vertices", {})
    pick = rng.choice(fault_v, size=min(5, fault_v.size), replace=False)
    sample_local = points[pick]
    sample_utm = local_to_utm(sample_local,
                               origin=tuple(transform["origin_utm_m"]))
    metrics = {"sample_utm": sample_utm.tolist()}
    if cfm_bbox is None:
        return CheckResult("9_utm_round_trip", True,
                            f"recovered UTM coords for {len(pick)} sampled "
                            f"vertices (no CFM bbox to compare)",
                            metrics)
    # Combine all included CFM bboxes into a union, allow 0 m margin.
    ux_min = min(b["xmin"] for b in cfm_bbox.values())
    ux_max = max(b["xmax"] for b in cfm_bbox.values())
    uy_min = min(b["ymin"] for b in cfm_bbox.values())
    uy_max = max(b["ymax"] for b in cfm_bbox.values())
    uz_min = min(b["zmin"] for b in cfm_bbox.values())
    uz_max = max(b["zmax"] for b in cfm_bbox.values())
    # Note: the local fault bbox is shrunk vs CFM (clearance + clamp), so
    # local_to_utm should land inside CFM bbox after re-adding origin.
    # Allow slop = 1 m for floating-point.
    slop = 1.0
    for v in sample_utm:
        x, y, z = float(v[0]), float(v[1]), float(v[2])
        # x, y must lie inside CFM bbox; z is clamped so local zmax = -1 m
        # ⇒ utm zmax = -1 m, which may be below cfm zmin if cfm zmax > 0.
        # Hence we only require utm.x ∈ [ux_min, ux_max] and utm.y ∈ [uy_min, uy_max].
        if not (ux_min - slop <= x <= ux_max + slop):
            return CheckResult("9_utm_round_trip", False,
                               f"recovered x={x:.0f} outside CFM "
                               f"[{ux_min:.0f}, {ux_max:.0f}]",
                               metrics)
        if not (uy_min - slop <= y <= uy_max + slop):
            return CheckResult("9_utm_round_trip", False,
                               f"recovered y={y:.0f} outside CFM "
                               f"[{uy_min:.0f}, {uy_max:.0f}]",
                               metrics)
    return CheckResult("9_utm_round_trip", True,
                       f"{len(pick)} sampled vertices land inside CFM bbox",
                       metrics)


def check_10_tet_quality(points, tets, tris=None, ttags=None,
                          tube_radius: float = 0.0) -> CheckResult:
    """Tet quality (gamma).  R-402: when `tube_radius > 0` and
    fault triangles are available, ALSO assert that no tet within
    `tube_radius` of any fault triangle has gamma < 0.05 — slivers
    in the SEAS friction-active zone are fatal even if they are a
    tiny fraction of the global mesh."""
    if tets.shape[0] == 0:
        return CheckResult("10_tet_quality", False, "no tets", {})
    # gamma = 12 * (3*V)^(2/3) / sum(edge^2) -- a standard tet quality measure
    # Use the "minimum dihedral angle" bound via volume / (sum of edge sqs).
    V = points
    a = V[tets[:, 0]]; b = V[tets[:, 1]]
    c = V[tets[:, 2]]; d = V[tets[:, 3]]
    signed_vol = np.einsum("ij,ij->i",
                              np.cross(b - a, c - a), d - a) / 6.0
    vol = np.abs(signed_vol)
    edges_sq = (
        np.sum((b - a) ** 2, axis=1) +
        np.sum((c - a) ** 2, axis=1) +
        np.sum((d - a) ** 2, axis=1) +
        np.sum((c - b) ** 2, axis=1) +
        np.sum((d - b) ** 2, axis=1) +
        np.sum((d - c) ** 2, axis=1)
    )
    # gamma_tet ≈ 12 * (3*V)^(2/3) / sum(edge^2); equilateral -> 1
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = 12.0 * np.cbrt((3.0 * vol) ** 2) / edges_sq
        gamma = np.where(np.isfinite(gamma), gamma, 0.0)
    g_min = float(gamma.min())
    g_mean = float(gamma.mean())
    n_bad = int(np.sum(gamma < 0.05))
    n_total = int(tets.shape[0])
    bad_frac = n_bad / max(n_total, 1)
    # Global minimum tet edge — useful for CFL diagnostics and as a
    # second quality metric independent of γ.
    edge_lengths_sq = np.stack([
        np.sum((b - a) ** 2, axis=1),
        np.sum((c - a) ** 2, axis=1),
        np.sum((d - a) ** 2, axis=1),
        np.sum((c - b) ** 2, axis=1),
        np.sum((d - b) ** 2, axis=1),
        np.sum((d - c) ** 2, axis=1),
    ], axis=1)
    min_edge_m = float(np.sqrt(edge_lengths_sq.min())) \
        if edge_lengths_sq.size else 0.0
    # Inverted (negative-volume) tet count.  γ uses |vol| so a fully
    # inverted tet can still pass the sliver gate; an inverted tet is
    # a HARD FAIL for any DG solver because the Jacobian determinant
    # is negative on a finite-measure set.  Threshold is strict
    # (signed_vol <= 0); a "zero-volume" tet is geometrically
    # degenerate and equally fatal.
    n_inverted = int(np.sum(signed_vol <= 0.0))
    worst_inverted_idx = -1
    worst_inverted_signed_vol = 0.0
    if n_inverted > 0:
        idx = int(np.argmin(signed_vol))
        worst_inverted_idx = idx
        worst_inverted_signed_vol = float(signed_vol[idx])
    metrics = {"gamma_min": g_min, "gamma_mean": g_mean,
                "n_tets": n_total,
                "n_gamma_below_0p05": n_bad,
                "frac_gamma_below_0p05": bad_frac,
                "min_edge_m": min_edge_m,
                "n_inverted_tets": n_inverted,
                "worst_inverted_tet_idx": worst_inverted_idx,
                "worst_inverted_tet_signed_vol_m3":
                    worst_inverted_signed_vol}
    # Inverted tets are categorical: even one is fatal.
    if n_inverted > 0:
        # Worst-case centroid for downstream debugging.
        ctr = (a[worst_inverted_idx] + b[worst_inverted_idx]
               + c[worst_inverted_idx] + d[worst_inverted_idx]) / 4.0
        return CheckResult(
            "10_tet_quality", False,
            f"{n_inverted} inverted tet(s) (signed_vol <= 0); worst "
            f"signed_vol={worst_inverted_signed_vol:.3e} m^3 at "
            f"centroid ({float(ctr[0]):.0f},{float(ctr[1]):.0f},"
            f"{float(ctr[2]):.0f}).  This is fatal for any DG solver.",
            metrics)
    # Aggregate quality is what matters for solver stability.  HXT
    # constraint recovery on realistic CFM geometry is allowed to leave
    # a handful of slivers as long as they are not dominant.
    if g_mean < 0.30:
        return CheckResult("10_tet_quality", False,
                           f"mean gamma {g_mean:.3f} below 0.30 — overall "
                           f"poor tet quality",
                           metrics)
    if bad_frac > 0.01:
        return CheckResult("10_tet_quality", False,
                           f"{bad_frac*100:.2f}% of tets have gamma<0.05 "
                           f"({n_bad}/{n_total}); above 1% threshold",
                           metrics)

    # R-402: near-fault-sliver gate.  Sliver tets in the far field
    # are tolerable; sliver tets within `tube_radius` of the fault
    # live in the SEAS friction-active zone where DG flux assembly
    # produces large numerical noise on poorly-shaped tets.  This
    # check is opt-in (only runs when fault tris + tube_radius are
    # supplied by the caller).
    if (tube_radius > 0.0 and tris is not None and ttags is not None
            and tris.shape[0] > 0):
        fault_mask = ttags == FAULT_TAG
        fault_tris = tris[fault_mask]
        if fault_tris.shape[0] > 0:
            from scipy.spatial import cKDTree
            fault_centroids = points[fault_tris].mean(axis=1)
            tree = cKDTree(fault_centroids)
            tet_centroids = points[tets].mean(axis=1)
            d_to_fault, _ = tree.query(tet_centroids, k=1)
            in_tube = d_to_fault <= tube_radius
            near_fault_slivers = int(np.sum(in_tube & (gamma < 0.05)))
            metrics["near_fault_slivers"] = near_fault_slivers
            metrics["near_fault_sliver_tube_radius_m"] = float(tube_radius)
            if near_fault_slivers > 0:
                near_idx = np.where(in_tube & (gamma < 0.05))[0]
                worst = int(near_idx[np.argmin(gamma[near_idx])])
                worst_gamma = float(gamma[worst])
                worst_d = float(d_to_fault[worst])
                metrics["worst_near_fault_sliver_gamma"] = worst_gamma
                metrics["worst_near_fault_sliver_distance_m"] = worst_d
                return CheckResult(
                    "10_tet_quality", False,
                    f"{near_fault_slivers} sliver tet(s) within "
                    f"tube_radius={tube_radius} m of the fault have "
                    f"gamma<0.05 (worst: gamma={worst_gamma:.4f} at "
                    f"d={worst_d:.0f} m); SEAS DG flux assembly will "
                    f"produce numerical noise on these — re-mesh with "
                    f"tighter size-field bias or run safs.geo's "
                    f"OptimizeNetgen at a higher threshold.",
                    metrics)

    msg = f"min gamma={g_min:.4f}, mean gamma={g_mean:.3f}"
    if g_min < 0.10:
        msg += (f" — WARNING: {n_bad} sliver tet(s) below 0.05 "
                f"(<{bad_frac*100:.2f}% of mesh)")
    return CheckResult("10_tet_quality", True, msg, metrics)


def check_11_tube_uniformity(points, tets, tris, ttags,
                              res_f: float, tube_radius: float,
                              ) -> CheckResult:
    """PLAN Phase 3 — within `tube_radius` of any fault triangle, the
    bulk tet edge length should be approximately uniform at `res_f`.

    Skipped (PASS trivially) when ``tube_radius <= 0``: the legacy
    "ramp from 0" sizing field doesn't have a uniform shell, so the
    test does not apply.
    """
    if tube_radius <= 0.0:
        return CheckResult(
            "11_tube_uniformity", True,
            "tube_radius=0 → legacy ramp-from-0 size field; check skipped",
            {"tube_radius_m": tube_radius},
        )
    mask = ttags == FAULT_TAG
    fault_tris = tris[mask]
    if fault_tris.shape[0] == 0:
        return CheckResult("11_tube_uniformity", False,
                            "no fault triangles", {})
    fault_centroids = points[fault_tris].mean(axis=1)
    from scipy.spatial import cKDTree
    t = cKDTree(fault_centroids)

    # Tet centroid → distance to nearest fault triangle.
    tet_centroids = points[tets].mean(axis=1)
    d_centroid, _ = t.query(tet_centroids, k=1)
    in_tube = d_centroid <= tube_radius
    in_tube_tets = tets[in_tube]
    metrics: dict = {
        "tube_radius_m": tube_radius,
        "res_f_m": res_f,
        "n_in_tube_tets": int(in_tube_tets.shape[0]),
    }
    if in_tube_tets.shape[0] == 0:
        return CheckResult(
            "11_tube_uniformity", False,
            f"no tets inside tube_radius={tube_radius} — sizing field "
            f"may not have produced any near-fault refinement",
            metrics,
        )
    a, b, c, d = (points[in_tube_tets[:, k]] for k in range(4))
    edges = np.concatenate([
        np.linalg.norm(b - a, axis=1),
        np.linalg.norm(c - a, axis=1),
        np.linalg.norm(d - a, axis=1),
        np.linalg.norm(c - b, axis=1),
        np.linalg.norm(d - b, axis=1),
        np.linalg.norm(d - c, axis=1),
    ])
    mean_e = float(edges.mean())
    p95_e = float(np.percentile(edges, 95))
    metrics["mean_edge_in_tube_m"] = mean_e
    metrics["p95_edge_in_tube_m"] = p95_e
    # R-401: tightened acceptance.  Was [0.5, 1.5] × res_f and
    # p95 ≤ 2.5 × res_f — too loose to detect a 38 % uniform-shift
    # (the original Phase-3 run had mean=1383 m for res_f=1000 m,
    # which the old bound passed silently).  The new bound asserts
    # the mean is within ±25 % of res_f AND that at least 70 % of
    # in-tube tets are individually within ±25 % of res_f (the
    # tube is uniform, not just centred-on-target on average).
    if not (0.75 * res_f <= mean_e <= 1.25 * res_f):
        return CheckResult(
            "11_tube_uniformity", False,
            f"mean tube tet edge {mean_e:.0f} m outside "
            f"[{0.75*res_f:.0f}, {1.25*res_f:.0f}] (res_f={res_f}); "
            f"sizing field is shifting the mean by more than ±25 % — "
            f"the tube exists but is not at res_f.  Adjust "
            f"`Mesh.MeshSizeFactor` and `Field[2].SizeMin` in "
            f"safs.geo (see R-401).",
            metrics,
        )
    if p95_e > 1.5 * res_f:
        return CheckResult(
            "11_tube_uniformity", False,
            f"p95 tube tet edge {p95_e:.0f} m exceeds {1.5*res_f:.0f} "
            f"— too much variation inside the uniform shell.",
            metrics,
        )
    n_uniform = int(np.sum((edges >= 0.75 * res_f) & (edges <= 1.25 * res_f)))
    frac_uniform = n_uniform / max(int(edges.size), 1)
    metrics["frac_in_tube_within_25pct"] = frac_uniform
    if frac_uniform < 0.70:
        return CheckResult(
            "11_tube_uniformity", False,
            f"only {frac_uniform*100:.1f}% of in-tube tet edges are "
            f"within ±25 % of res_f (need >= 70 %); the tube has the "
            f"right mean but is internally non-uniform.",
            metrics,
        )
    return CheckResult(
        "11_tube_uniformity", True,
        f"{in_tube_tets.shape[0]} in-tube tets, mean edge "
        f"{mean_e:.0f} m, p95 {p95_e:.0f} m, "
        f"{frac_uniform*100:.0f}% within ±25%",
        metrics,
    )


def check_12_surface_closure(points, tris, ttags, tets) -> CheckResult:
    """Every 1-tet boundary face MUST have a matching tagged triangle.

    Catches:
      - Surface holes (mmg3d_local_patch regression of 2026-05-02)
      - Lost box-face triangles after re-meshing
      - Cavity-retet stitch errors
    """
    if tets.shape[0] == 0:
        return CheckResult("12_surface_closure", False, "no tetrahedra", {})
    face_count: dict[frozenset, int] = {}
    for ti in range(tets.shape[0]):
        a, b, c, d = (int(x) for x in tets[ti])
        for face in (frozenset((a, b, c)), frozenset((a, b, d)),
                      frozenset((a, c, d)), frozenset((b, c, d))):
            face_count[face] = face_count.get(face, 0) + 1
    bdry = {f for f, n in face_count.items() if n == 1}
    tagged = {frozenset((int(tris[ti, 0]), int(tris[ti, 1]),
                          int(tris[ti, 2])))
              for ti in range(tris.shape[0])}
    holes = bdry - tagged
    metrics = {"n_bdry_faces": len(bdry),
               "n_tagged_tris": len(tagged),
               "n_unlabeled_holes": len(holes)}
    if holes:
        ex = []
        for f in list(holes)[:3]:
            vs = list(f)
            c = points[vs].mean(axis=0)
            ex.append(f"({c[0]:.0f},{c[1]:.0f},{c[2]:.0f})")
        return CheckResult(
            "12_surface_closure", False,
            f"{len(holes)} unlabeled 1-tet bdry face(s); "
            f"first centroids: {', '.join(ex)}",
            metrics)
    return CheckResult("12_surface_closure", True,
                       f"all {len(bdry)} bdry faces have tagged tris",
                       metrics)


def check_13_fault_orientation(points, tris, ttags) -> CheckResult:
    """Fault tris must be consistently wound: every interior fault edge
    is traversed in OPPOSITE directions by the two incident fault tris.
    Counts the number of winding-flip mismatches via BFS propagation
    from a seed tri per connected component.  Also reports non-manifold
    fault edges (>2 fault tris share the edge) which cannot be cleanly
    oriented without splitting the fault into separate tags.
    """
    fault_idx = np.where(ttags == FAULT_TAG)[0]
    if fault_idx.size == 0:
        return CheckResult("13_fault_orientation", False,
                            "no fault triangles", {})
    fault_tris = [(int(tris[ti, 0]), int(tris[ti, 1]),
                    int(tris[ti, 2])) for ti in fault_idx]
    n_ft = len(fault_tris)
    edge2tri: dict[frozenset, list[int]] = defaultdict(list)
    for ti, (a, b, c) in enumerate(fault_tris):
        for u, v in ((a, b), (b, c), (c, a)):
            edge2tri[frozenset((u, v))].append(ti)
    n_nonmanifold = sum(1 for vs in edge2tri.values() if len(vs) > 2)

    # Dihedral-pair table for non-manifold (count >= 3) edges.
    # Mirrors orient_fault_surface._orient_fault: pair tris by
    # anti-parallel third-vertex bearings (smooth manifold
    # continuation across the edge).  Through-going pairs at a
    # T-junction (count==3) and the two pairs at an X-junction
    # (count==4) are propagated; terminator tris are skipped.
    partner: dict[tuple[int, frozenset], int] = {}
    for e, tilist in edge2tri.items():
        if len(tilist) <= 2:
            continue
        uv = list(e)
        u_, v_ = uv[0], uv[1]
        edge_vec = points[v_] - points[u_]
        edge_len = float(np.linalg.norm(edge_vec))
        if edge_len < 1e-12:
            continue
        edge_unit = edge_vec / edge_len
        bearings: list = []
        for ti in tilist:
            a, b, c = fault_tris[ti]
            third = (a if a not in (u_, v_)
                     else (b if b not in (u_, v_) else c))
            r = points[third] - points[u_]
            rp = r - float(np.dot(r, edge_unit)) * edge_unit
            n_ = float(np.linalg.norm(rp))
            bearings.append(rp / n_ if n_ > 1e-12 else None)
        used = [False] * len(tilist)
        for i in range(len(tilist)):
            if used[i] or bearings[i] is None:
                continue
            best_j, best_d = -1, 1.0
            for j in range(i + 1, len(tilist)):
                if used[j] or bearings[j] is None:
                    continue
                d = float(np.dot(bearings[i], bearings[j]))
                if d < best_d:
                    best_d = d
                    best_j = j
            if best_j < 0 or best_d > -0.5:
                used[i] = True
                continue
            partner[(tilist[i], e)] = tilist[best_j]
            partner[(tilist[best_j], e)] = tilist[i]
            used[i] = True
            used[best_j] = True

    visited = [False] * n_ft
    n_flips = 0
    n_components = 0
    for start in range(n_ft):
        if visited[start]:
            continue
        n_components += 1
        visited[start] = True
        queue = [start]
        while queue:
            ti = queue.pop()
            a, b, c = fault_tris[ti]
            for u, v in ((a, b), (b, c), (c, a)):
                e = frozenset((u, v))
                e_tris = edge2tri[e]
                cnt = len(e_tris)
                if cnt == 1:
                    continue
                if cnt == 2:
                    nti_candidates = [n for n in e_tris if n != ti]
                else:
                    p = partner.get((ti, e))
                    if p is None:
                        continue
                    nti_candidates = [p]
                for nti in nti_candidates:
                    if nti == ti or visited[nti]:
                        continue
                    na, nb, nc = fault_tris[nti]
                    same_dir = ((na == u and nb == v)
                                 or (nb == u and nc == v)
                                 or (nc == u and na == v))
                    if same_dir:
                        n_flips += 1
                    visited[nti] = True
                    queue.append(nti)
    metrics = {"n_fault_tris": n_ft,
               "n_components": n_components,
               "n_winding_flips_needed": int(n_flips),
               "n_nonmanifold_edges": int(n_nonmanifold)}
    # R-002 spec: HARD-fail only on winding flips (orientation defect that
    # orient_fault_surface.py must repair upstream).  Non-manifold edges
    # are unavoidable for branching multi-fault SAFS surfaces — report
    # the count as a metric but do NOT fail validation.  Per REVIEW.md
    # R-002: "Hard fail if n_flips_needed > 0 ...; soft warn if non-
    # manifold edges > 0 (cannot always be repaired without changing
    # the fault geometry)."
    if n_flips > 0:
        return CheckResult(
            "13_fault_orientation", False,
            f"{n_flips} winding flips + {n_nonmanifold} non-manifold "
            f"edges across {n_components} fault component(s)",
            metrics)
    msg = (f"{n_components} fault component(s) consistently wound; "
           f"{n_nonmanifold} non-manifold edge(s) (branching, ok)")
    return CheckResult("13_fault_orientation", True, msg, metrics)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate a SAFS .msh against the smoke-test suite.",
    )
    parser.add_argument("--msh", required=True, type=Path)
    parser.add_argument("--transform-json", default=None, type=Path)
    parser.add_argument("--provenance-json", default=None, type=Path)
    parser.add_argument("--bbox-json", default=None, type=Path,
                        help="Local-frame per-fault bbox.json from "
                             "ts_to_stl.py.  Used by check 3 (clearance) "
                             "and check 9 (UTM round-trip).")
    parser.add_argument("--domain-box-json", default=None, type=Path)
    parser.add_argument("--sizing-json", default=None, type=Path)
    parser.add_argument("--expected-faults", action="append", default=None,
                        metavar="SHORT_NAME",
                        help="Short names that must appear in the "
                             "provenance JSON.  May be repeated.")
    parser.add_argument("--expected-box", default=None, type=Path,
                        help="JSON file with expected domain_box keys "
                             "(used by check 2).  If omitted, check 2 is "
                             "passed trivially with the actual values.")
    parser.add_argument("--report", default=None, type=Path)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args(argv)

    # Defaults
    out_dir = args.msh.parent
    transform_json = args.transform_json or (out_dir.parent /
                                              "transform.json")
    provenance_json = args.provenance_json or (out_dir / "fault_provenance.json")
    bbox_json = args.bbox_json or (out_dir.parent / "bbox.json")
    domain_box_json = args.domain_box_json or (out_dir / "domain_box.json")
    sizing_json = args.sizing_json or (out_dir / "sizing.json")

    # Mesh data
    points, tris, ttags, tets, tetags = _mesh_data(args.msh)

    # Side-cars
    transform = read_transform_json(transform_json)
    fault_bbox = json.loads(bbox_json.read_text()) if bbox_json.exists() else {}
    if args.expected_faults:
        keep = set(args.expected_faults)
        fault_bbox = {k: v for k, v in fault_bbox.items() if k in keep}
    domain_box = json.loads(domain_box_json.read_text()) \
        if domain_box_json.exists() else {}
    sizing = json.loads(sizing_json.read_text()) \
        if sizing_json.exists() else {}
    res_f = float(sizing.get("res_f_m", 1000.0))
    res_ff = float(sizing.get("res_ff_m", 20000.0))

    provenance = {}
    if provenance_json.exists():
        provenance = json.loads(provenance_json.read_text())

    # CFM-UTM bbox for check 9 (recover via local + origin)
    cfm_bbox = None
    if fault_bbox:
        origin = tuple(transform["origin_utm_m"])
        cfm_bbox = {}
        for short, b in fault_bbox.items():
            cfm_bbox[short] = {
                "xmin": b["xmin"] + origin[0], "xmax": b["xmax"] + origin[0],
                "ymin": b["ymin"] + origin[1], "ymax": b["ymax"] + origin[1],
                "zmin": b["zmin"] + origin[2], "zmax": b["zmax"] + origin[2],
            }

    rng = np.random.default_rng(args.seed)

    # Run checks
    results: list[CheckResult] = []

    results.append(check_1_tag_inventory(ttags, tetags))

    if args.expected_box and args.expected_box.exists():
        expected = json.loads(args.expected_box.read_text())
        results.append(check_2_box_arithmetic(domain_box, expected))
    elif domain_box:
        results.append(CheckResult(
            "2_box_arithmetic", True,
            "no --expected-box provided; trivially passed with actual",
            {"actual": domain_box}))
    else:
        results.append(CheckResult(
            "2_box_arithmetic", False,
            "domain_box.json not found", {}))

    if domain_box and fault_bbox:
        results.append(check_3_fault_clearance(
            domain_box, fault_bbox,
            float(domain_box.get("buf_x_m", 50000.0)),
            float(domain_box.get("buf_y_m", 50000.0)),
            float(domain_box.get("depth_m", 50000.0)),
        ))
    else:
        results.append(CheckResult(
            "3_fault_clearance", False,
            "missing domain_box.json or bbox.json", {}))

    clearance = float(transform.get("free_surface_clearance_m", 1.0))
    results.append(check_4_freesurface_trace(points, tris, ttags,
                                              clearance=clearance))
    results.append(check_5_internal_interface(points, tris, ttags, tets,
                                                 clearance=clearance))
    results.append(check_6_fault_edge_length(points, tris, ttags, res_f))
    results.append(check_7_far_field_edge(points, tets, tris, ttags,
                                            res_ff))

    if provenance:
        results.append(check_8_provenance_partition(provenance, ttags))
    else:
        results.append(CheckResult(
            "8_provenance_partition", False,
            "fault_provenance.json not found at "
            f"{provenance_json}", {}))

    results.append(check_9_utm_round_trip(points, tris, ttags,
                                            transform, cfm_bbox, rng))
    tube_radius = float(sizing.get("tube_radius_m", 0.0))
    # R-402: pass fault tris + tube_radius into check_10 so the
    # near-fault-sliver gate fires.  Sliver tets in the friction-active
    # zone (within tube_radius of any fault triangle) are fatal even
    # when their global fraction is small.
    results.append(check_10_tet_quality(
        points, tets, tris=tris, ttags=ttags, tube_radius=tube_radius))
    results.append(check_11_tube_uniformity(
        points, tets, tris, ttags, res_f, tube_radius,
    ))
    results.append(check_12_surface_closure(points, tris, ttags, tets))
    results.append(check_13_fault_orientation(points, tris, ttags))

    # Report
    report_lines = []
    report_lines.append(f"validate_msh.py report — {args.msh}")
    report_lines.append("=" * 64)
    n_pass = 0
    for r in results:
        flag = "PASS" if r.passed else "FAIL"
        if r.passed:
            n_pass += 1
        report_lines.append(f"[{flag}] {r.name}: {r.message}")
        for k, v in r.metrics.items():
            report_lines.append(f"        {k}: {v}")
    report_lines.append("=" * 64)
    report_lines.append(f"{n_pass}/{len(results)} checks passed")
    report = "\n".join(report_lines) + "\n"
    print(report)

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(report)

    # Provenance seal
    seal = {
        "git_commit": _git_commit(),
        "gmsh_version": _gmsh_version(),
        "msh_sha256": _sha256(args.msh),
        "n_points": int(points.shape[0]),
        "n_triangles": int(tris.shape[0]),
        "n_tetrahedra": int(tets.shape[0]),
    }
    seal_path = args.msh.with_suffix(".provenance.json")
    seal_path.write_text(json.dumps(seal, indent=2) + "\n")

    return 0 if n_pass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
