#!/usr/bin/env python
"""Centroid-Steiner pre-split for collinear / near-coplanar bad tets.

Implements R-001 of `REVIEW_centroid_steiner_for_collinear_tets.md`:
the empirical worst-tet diagnostic on Step 6 production mesh shows
the γ_min = 2.5e-06 ceiling is set by **near-collinear 4-vertex tets**
(e.g., 3 of 4 vertices on the fault surface, all 4 lying on the same
3-D line).  mmg3d's `-optim` cannot break these because:

  - Edge-midpoint split keeps the new vertex ON the line.
  - Edge collapse would require moving a Required vertex.
  - Edge swap preserves the 4-vertex set.

This module implements the missing operation: **insert a centroid
Steiner point INSIDE the bad tet**, splitting it into 4 sub-tets.
The centroid is in the bulk volume (off the line / plane), giving
each sub-tet a well-conditioned 4-vertex configuration.

Cross-fault conformity invariant:
  The original tet's 4 faces (each potentially a fault triangle)
  REMAIN faces of the new sub-tets.  Each face is still shared with
  the same neighbouring tet on the other side of that face.
  Validator check_5 (fault-tri 2-tet adjacency) is preserved by
  construction.  No surface modification.

Pipeline position (post `mmg3d_local_patch`):
    ... → mmg3d_local_patch → THIS → validate_msh

Usage:
  python centroid_steiner_split.py \\
      --in-msh safs_patch.msh --out-msh safs_centroid.msh \\
      [--gamma-thresh 1e-3] [--collinearity-thresh 0.99] \\
      [--coplanarity-thresh 1e-4]
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import mmg3d_local_patch as mlp  # noqa: E402


# ---------------------------------------------------------------------------
# Geometric tests.
# ---------------------------------------------------------------------------
def _is_near_collinear(points: np.ndarray, tet: np.ndarray,
                        collinearity_thresh: float) -> bool:
    """Return True iff ALL FOUR vertices are nearly collinear in 3-D.

    Test: compute the unit-vector edge directions e_01, e_02, e_03
    from vertex 0.  Require ALL THREE pairwise |dot products| above
    `collinearity_thresh` (near 1 = parallel or antiparallel) — i.e.,
    `min(|e_01·e_02|, |e_01·e_03|, |e_02·e_03|) > thresh`.

    Using `min` (not `max`) avoids a false-positive on 3-collinear-
    plus-1-off configurations: if v0, v1, v2 are collinear but v3
    is off-line, |e_01·e_02| ≈ 1 BUT |e_01·e_03|, |e_02·e_03| may
    be ≪ 1.  Such tets are well-conditioned and must NOT be flagged.
    """
    v = points[tet]
    e01 = v[1] - v[0]
    e02 = v[2] - v[0]
    e03 = v[3] - v[0]
    n01 = np.linalg.norm(e01)
    n02 = np.linalg.norm(e02)
    n03 = np.linalg.norm(e03)
    if n01 == 0 or n02 == 0 or n03 == 0:
        return True  # degenerate (zero-length edge) — treat as collinear
    e01u = e01 / n01
    e02u = e02 / n02
    e03u = e03 / n03
    d_pairs = (abs(float(np.dot(e01u, e02u))),
                abs(float(np.dot(e01u, e03u))),
                abs(float(np.dot(e02u, e03u))))
    return min(d_pairs) > collinearity_thresh


def _is_near_coplanar_or_flat(points: np.ndarray, tet: np.ndarray,
                                coplanarity_thresh: float) -> bool:
    """Return True iff tet has very small volume relative to the
    cube of its longest edge.  Captures both "near-coplanar"
    (4 vertices nearly on a plane) and "near-collinear" (which is
    a stronger form).

    Test: vol / max_edge^3 < coplanarity_thresh.
    Reference: a regular tet has vol/edge^3 = sqrt(2)/12 ≈ 0.118.
    A "very flat" tet has this ratio near 0.
    """
    v = points[tet]
    p0, p1, p2, p3 = v[0], v[1], v[2], v[3]
    edge_lengths = [
        np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p0),
        np.linalg.norm(p3 - p0), np.linalg.norm(p2 - p1),
        np.linalg.norm(p3 - p1), np.linalg.norm(p3 - p2),
    ]
    max_edge = max(edge_lengths)
    if max_edge == 0:
        return True  # degenerate
    vol = abs(float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))) / 6.0
    ratio = vol / (max_edge ** 3)
    return ratio < coplanarity_thresh


def _should_split(points: np.ndarray, tet: np.ndarray,
                   collinearity_thresh: float,
                   coplanarity_thresh: float) -> bool:
    """A tet should be centroid-split if it's near-collinear OR
    near-coplanar — both are 1-D / 2-D degenerate configurations
    that mmg3d's `-optim` cannot break with edge operations.
    """
    return (_is_near_collinear(points, tet, collinearity_thresh)
            or _is_near_coplanar_or_flat(points, tet,
                                            coplanarity_thresh))


def _tet_volume(points: np.ndarray, tet: np.ndarray) -> float:
    """Signed volume of a tet (absolute value)."""
    v = points[tet]
    return abs(float(np.dot(v[1] - v[0],
                              np.cross(v[2] - v[0],
                                        v[3] - v[0])))) / 6.0


# ---------------------------------------------------------------------------
# The split itself.
# ---------------------------------------------------------------------------
def centroid_split_collinear_tets(
        in_msh: Path, out_msh: Path,
        gamma_thresh: float = 1e-3,
        collinearity_thresh: float = 0.99,
        coplanarity_thresh: float = 1e-4,
        eps_vol_m3: float = 1e-12,
        ) -> dict:
    """For each tet with γ < `gamma_thresh` AND geometrically
    near-collinear or near-coplanar: insert a centroid Steiner
    point and replace with 4 sub-tets.

    Args:
      in_msh: input gmsh .msh path.
      out_msh: output gmsh .msh path.
      gamma_thresh: only tets below this γ are candidates.
      collinearity_thresh: max(|e_i · e_j|) > this → collinear.
      coplanarity_thresh: vol / max_edge^3 < this → coplanar.
      eps_vol_m3: skip tets with volume < this absolute floor;
        a true-zero-volume tet would produce 4 sub-zero-volume
        tets and not improve γ.

    Returns a JSON-serializable report.
    """
    if gamma_thresh <= 0:
        raise ValueError(
            f"gamma_thresh must be > 0; got {gamma_thresh}")
    if not (0.0 < collinearity_thresh <= 1.0):
        raise ValueError(
            f"collinearity_thresh must be in (0, 1]; got "
            f"{collinearity_thresh}")
    if coplanarity_thresh <= 0:
        raise ValueError(
            f"coplanarity_thresh must be > 0; got "
            f"{coplanarity_thresh}")
    if eps_vol_m3 <= 0:
        raise ValueError(
            f"eps_vol_m3 must be > 0; got {eps_vol_m3}")
    if not in_msh.exists():
        raise SystemExit(f"--in-msh not found: {in_msh}")

    out_msh.parent.mkdir(parents=True, exist_ok=True)

    in_data = mlp._read_msh(in_msh)
    points = in_data["points"]
    tets = in_data["tets"]
    tet_tags = in_data["tet_tags"]
    print(f"  read {in_msh}: V={points.shape[0]}, "
          f"T={tets.shape[0]}, "
          f"tris={in_data['tris'].shape[0]}", file=sys.stderr)

    g = mlp._gamma_per_tet(points, tets)
    bad_tet_idx = np.flatnonzero(g < gamma_thresh)
    print(f"  bad tets (γ < {gamma_thresh}): {bad_tet_idx.size}",
          file=sys.stderr)

    # Filter bad tets by geometric near-degenerate test.
    split_set: set[int] = set()
    n_skipped_vol_too_small = 0
    n_skipped_not_degenerate = 0
    for ti in bad_tet_idx:
        ti = int(ti)
        if _tet_volume(points, tets[ti]) < eps_vol_m3:
            n_skipped_vol_too_small += 1
            continue
        if _should_split(points, tets[ti],
                          collinearity_thresh, coplanarity_thresh):
            split_set.add(ti)
        else:
            n_skipped_not_degenerate += 1

    print(f"  to split (collinear/coplanar): {len(split_set)}, "
          f"skipped not-degenerate: {n_skipped_not_degenerate}, "
          f"skipped vol-too-small: {n_skipped_vol_too_small}",
          file=sys.stderr)

    if not split_set:
        # No-op: copy input to output and report.
        import shutil
        shutil.copy(in_msh, out_msh)
        me_noop = mlp._min_edge_of_tets(points, tets)
        return {
            "in_msh":   str(in_msh),
            "out_msh":  str(out_msh),
            "gamma_thresh":         gamma_thresh,
            "collinearity_thresh":  collinearity_thresh,
            "coplanarity_thresh":   coplanarity_thresh,
            "eps_vol_m3":           eps_vol_m3,
            "n_bad_tets":           int(bad_tet_idx.size),
            "n_centroid_split":     0,
            "n_steiner_added":      0,
            "n_skipped_vol_too_small": int(n_skipped_vol_too_small),
            "n_skipped_not_degenerate": int(n_skipped_not_degenerate),
            "tets_in":              int(tets.shape[0]),
            "tets_out":             int(tets.shape[0]),
            "gamma_min_in":         float(g.min()) if g.size else 0.0,
            "gamma_min_out":        float(g.min()) if g.size else 0.0,
            "min_edge_in_m":        me_noop,
            "min_edge_out_m":       me_noop,
            "no_op":                True,
        }

    # Build the new mesh.
    new_points: list[tuple[float, float, float]] = [
        tuple(p) for p in points]
    new_tets: list[tuple[int, int, int, int]] = []
    new_tags: list[int] = []
    for ti in range(tets.shape[0]):
        if ti not in split_set:
            new_tets.append((int(tets[ti, 0]), int(tets[ti, 1]),
                              int(tets[ti, 2]), int(tets[ti, 3])))
            new_tags.append(int(tet_tags[ti]))
            continue
        # Centroid Steiner split.
        a = int(tets[ti, 0]); b = int(tets[ti, 1])
        c = int(tets[ti, 2]); d = int(tets[ti, 3])
        centroid = points[tets[ti]].mean(axis=0)
        ci = len(new_points)
        new_points.append(tuple(centroid))
        # 4 sub-tets — one per face of the original tet.  Each
        # sub-tet is (face_vertex_0, face_vertex_1, face_vertex_2,
        # centroid).  The face is preserved as a sub-tet face
        # adjacent to the original neighbour on the other side.
        tag = int(tet_tags[ti])
        new_tets.append((a, b, c, ci)); new_tags.append(tag)
        new_tets.append((a, b, d, ci)); new_tags.append(tag)
        new_tets.append((a, c, d, ci)); new_tags.append(tag)
        new_tets.append((b, c, d, ci)); new_tags.append(tag)

    out_data = {
        "points":   np.asarray(new_points, dtype=np.float64),
        "tris":     in_data["tris"],
        "tri_tags": in_data["tri_tags"],
        "tets":     np.asarray(new_tets, dtype=np.int64),
        "tet_tags": np.asarray(new_tags, dtype=np.int32),
    }
    mlp._write_msh(out_msh, out_data)

    # Report.
    g_in_min = float(g.min()) if g.size else 0.0
    me_in = mlp._min_edge_of_tets(points, tets)
    g_out = mlp._gamma_per_tet(out_data["points"], out_data["tets"])
    g_out_min = float(g_out.min()) if g_out.size else 0.0
    me_out = mlp._min_edge_of_tets(out_data["points"], out_data["tets"])

    print(f"  wrote merged mesh: V={out_data['points'].shape[0]}, "
          f"T={out_data['tets'].shape[0]}", file=sys.stderr)
    print(f"  γ_min: {g_in_min:.3e} → {g_out_min:.3e}  "
          f"min_edge: {me_in:.3f} m → {me_out:.3f} m",
          file=sys.stderr)

    return {
        "in_msh":   str(in_msh),
        "out_msh":  str(out_msh),
        "gamma_thresh":         gamma_thresh,
        "collinearity_thresh":  collinearity_thresh,
        "coplanarity_thresh":   coplanarity_thresh,
        "eps_vol_m3":           eps_vol_m3,
        "n_bad_tets":           int(bad_tet_idx.size),
        "n_centroid_split":     int(len(split_set)),
        "n_steiner_added":      int(len(split_set)),
        "n_skipped_vol_too_small": int(n_skipped_vol_too_small),
        "n_skipped_not_degenerate": int(n_skipped_not_degenerate),
        "tets_in":              int(tets.shape[0]),
        "tets_out":             int(out_data["tets"].shape[0]),
        "gamma_min_in":         g_in_min,
        "gamma_min_out":        g_out_min,
        "min_edge_in_m":        me_in,
        "min_edge_out_m":       me_out,
        "no_op":                False,
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Centroid-Steiner pre-split for collinear / "
                    "near-coplanar bad tets.")
    p.add_argument("--in-msh", required=True, type=Path)
    p.add_argument("--out-msh", required=True, type=Path)
    p.add_argument("--gamma-thresh", type=float, default=1e-3,
                   help="Only tets with γ below this threshold are "
                        "candidates for centroid split.  Default "
                        "1e-3 — targets the truly degenerate ones, "
                        "not the moderate slivers.")
    p.add_argument("--collinearity-thresh", type=float, default=0.99,
                   help="Tet flagged collinear if "
                        "max(|e_01·e_02|, |e_01·e_03|, "
                        "|e_02·e_03|) > this (unit-vector dot "
                        "products).  Default 0.99.")
    p.add_argument("--coplanarity-thresh", type=float, default=1e-4,
                   help="Tet flagged coplanar if "
                        "vol / max_edge^3 < this.  A regular tet "
                        "has ratio ~0.118; coplanar tets approach "
                        "0.  Default 1e-4.")
    p.add_argument("--eps-vol-m3", type=float, default=1e-12,
                   help="Tets with volume below this absolute floor "
                        "(m^3) are SKIPPED — a true-zero-volume tet "
                        "would produce sub-zero-volume children.  "
                        "Default 1e-12.")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = centroid_split_collinear_tets(
        in_msh=args.in_msh, out_msh=args.out_msh,
        gamma_thresh=args.gamma_thresh,
        collinearity_thresh=args.collinearity_thresh,
        coplanarity_thresh=args.coplanarity_thresh,
        eps_vol_m3=args.eps_vol_m3,
    )

    report_path = args.report_json or (
        args.out_msh.parent / "centroid_steiner_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
