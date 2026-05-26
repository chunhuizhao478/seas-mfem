#!/usr/bin/env python3
"""
Plot along-fault strike- vs dip-shear dominance from a fault-stress VTU.

The script reads cell-centred `tau_strike_MPa_cell` / `tau_dip_MPa_cell`,
projects triangle centroids onto a horizontal along-strike axis, bins the
fault along that coordinate, writes an augmented VTU, and writes:

- `<out_prefix>.vtu`   : fault VTU with added cell-data fields
- `<out_prefix>.png`   : two-panel diagnostic figure
- `<out_prefix>.csv`   : per-bin summary table
- `<out_prefix>.json`  : overall counts / extents / dominant intervals

Default along-fault axis uses the strike-hint azimuth recorded in the
summary JSON when available; otherwise it falls back to a horizontal PCA
fit to the triangle centroids.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np

_RATIO_EPS = 1.0e-12


def _load_fault_triangles(mesh: meshio.Mesh) -> np.ndarray:
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            return np.asarray(cell_block.data, dtype=np.int64)
    raise ValueError("fault VTU contains no triangle cells")


def _load_cell_array(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.cell_data:
        raise KeyError(f"fault VTU missing cell-data field {name!r}")
    for cell_block, values in zip(mesh.cells, mesh.cell_data[name]):
        if cell_block.type == "triangle":
            arr = np.asarray(values, dtype=np.float64).reshape(-1)
            if arr.shape[0] != cell_block.data.shape[0]:
                raise ValueError(
                    f"cell-data field {name!r} length {arr.shape[0]} does not "
                    f"match triangle count {cell_block.data.shape[0]}"
                )
            return arr
    raise ValueError(f"fault VTU has field {name!r} but no triangle block")


def _triangle_centroids(points: np.ndarray, tri_conn: np.ndarray) -> np.ndarray:
    tri_pts = np.asarray(points, dtype=np.float64)[tri_conn]
    return tri_pts.mean(axis=1)


def _strike_hat_xy(summary_json: Path | None, centroids_xy: np.ndarray) -> tuple[np.ndarray, str]:
    if summary_json is not None:
        with open(summary_json) as f:
            summary = json.load(f)
        az = summary.get("params", {}).get("fault_strike_azimuth_hint_deg")
        if az is not None and np.isfinite(az):
            az_rad = np.deg2rad(float(az))
            return np.array([np.sin(az_rad), np.cos(az_rad)], dtype=np.float64), "summary_strike_hint"

    xy = np.asarray(centroids_xy, dtype=np.float64)
    xy0 = xy - xy.mean(axis=0)
    _, _, vh = np.linalg.svd(xy0, full_matrices=False)
    axis = vh[0]
    if axis[1] < 0.0:
        axis = -axis
    return axis / np.linalg.norm(axis), "pca"


def _dominant_label(abs_tau_strike: float, abs_tau_dip: float) -> str:
    return "strike" if abs_tau_strike >= abs_tau_dip else "dip"


def _intervals_from_labels(starts: np.ndarray, ends: np.ndarray, labels: list[str]) -> list[dict]:
    intervals: list[dict] = []
    if len(labels) == 0:
        return intervals
    cur_label = labels[0]
    cur_start = float(starts[0])
    cur_end = float(ends[0])
    for start, end, label in zip(starts[1:], ends[1:], labels[1:]):
        if label == cur_label:
            cur_end = float(end)
            continue
        intervals.append(
            {
                "dominant": cur_label,
                "start_km": cur_start,
                "end_km": cur_end,
                "length_km": cur_end - cur_start,
            }
        )
        cur_label = label
        cur_start = float(start)
        cur_end = float(end)
    intervals.append(
        {
            "dominant": cur_label,
            "start_km": cur_start,
            "end_km": cur_end,
            "length_km": cur_end - cur_start,
        }
    )
    return intervals


def _build_bin_rows(
    along_km: np.ndarray,
    tau_strike: np.ndarray,
    tau_dip: np.ndarray,
    n_bins: int,
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    edges = np.linspace(float(along_km.min()), float(along_km.max()), n_bins + 1)
    rows: list[dict] = []
    for i in range(n_bins):
        left = edges[i]
        right = edges[i + 1]
        if i == n_bins - 1:
            mask = (along_km >= left) & (along_km <= right)
        else:
            mask = (along_km >= left) & (along_km < right)
        if not np.any(mask):
            continue
        ats = np.abs(tau_strike[mask])
        atd = np.abs(tau_dip[mask])
        strike_count = int(np.count_nonzero(ats >= atd))
        dip_count = int(mask.sum() - strike_count)
        rows.append(
            {
                "bin_start_km": float(left),
                "bin_end_km": float(right),
                "bin_center_km": float(0.5 * (left + right)),
                "n_cells": int(mask.sum()),
                "median_abs_tau_strike_MPa": float(np.median(ats)),
                "median_abs_tau_dip_MPa": float(np.median(atd)),
                "mean_abs_tau_strike_MPa": float(np.mean(ats)),
                "mean_abs_tau_dip_MPa": float(np.mean(atd)),
                "strike_dominant_fraction": float(strike_count / mask.sum()),
                "strike_dominant_count": strike_count,
                "dip_dominant_count": dip_count,
                "dominant": _dominant_label(float(np.median(ats)), float(np.median(atd))),
            }
        )
    return rows, edges[:-1], edges[1:]


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("no populated bins available to write")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict) -> None:
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _write_augmented_vtu(
    out_vtu: Path,
    mesh: meshio.Mesh,
    along_km: np.ndarray,
    tau_strike: np.ndarray,
    tau_dip: np.ndarray,
) -> None:
    abs_tau_strike = np.abs(tau_strike)
    abs_tau_dip = np.abs(tau_dip)
    margin = abs_tau_strike - abs_tau_dip
    ratio = abs_tau_strike / np.maximum(abs_tau_dip, _RATIO_EPS)
    flag = np.where(abs_tau_strike >= abs_tau_dip, 1, -1).astype(np.int8)

    cell_data = {
        key: [np.array(block, copy=True) for block in values]
        for key, values in (mesh.cell_data or {}).items()
    }

    float_fields = {
        "along_fault_km_cell": along_km.astype(np.float64),
        "abs_tau_strike_MPa_cell": abs_tau_strike.astype(np.float64),
        "abs_tau_dip_MPa_cell": abs_tau_dip.astype(np.float64),
        "shear_dominance_margin_MPa_cell": margin.astype(np.float64),
        "shear_dominance_ratio_cell": ratio.astype(np.float64),
    }
    int_fields = {
        "shear_dominance_flag_cell": flag,
    }

    for name, tri_values in float_fields.items():
        values = []
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                values.append(np.array(tri_values, copy=True))
            else:
                values.append(np.full(len(cell_block.data), np.nan, dtype=np.float64))
        cell_data[name] = values
    for name, tri_values in int_fields.items():
        values = []
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                values.append(np.array(tri_values, copy=True))
            else:
                values.append(np.zeros(len(cell_block.data), dtype=np.int8))
        cell_data[name] = values

    out_mesh = meshio.Mesh(
        points=np.array(mesh.points, copy=True),
        cells=[(cell_block.type, np.array(cell_block.data, copy=True)) for cell_block in mesh.cells],
        point_data={
            key: np.array(values, copy=True)
            for key, values in (mesh.point_data or {}).items()
        },
        cell_data=cell_data,
        field_data={
            key: np.array(values, copy=True)
            for key, values in (mesh.field_data or {}).items()
        },
    )
    out_mesh.info = dict(getattr(mesh, "info", {}) or {})
    out_mesh.info["shear_dominance_flag_cell_mapping"] = "-1=dip, +1=strike"
    out_mesh.info["along_fault_axis"] = "horizontal projection onto chosen strike axis"
    meshio.write(str(out_vtu), out_mesh)


def _plot(
    out_png: Path,
    along_km: np.ndarray,
    tau_strike: np.ndarray,
    tau_dip: np.ndarray,
    rows: list[dict],
    title: str,
) -> None:
    centers = np.array([row["bin_center_km"] for row in rows], dtype=np.float64)
    med_abs_tau_strike = np.array(
        [row["median_abs_tau_strike_MPa"] for row in rows], dtype=np.float64
    )
    med_abs_tau_dip = np.array(
        [row["median_abs_tau_dip_MPa"] for row in rows], dtype=np.float64
    )
    strike_frac = np.array(
        [row["strike_dominant_fraction"] for row in rows], dtype=np.float64
    )
    dominant = [row["dominant"] for row in rows]

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1.4]}
    )

    ax0.scatter(along_km, np.abs(tau_strike), s=4, alpha=0.18, c="#1f77b4", label="|tau_strike| cells")
    ax0.scatter(along_km, np.abs(tau_dip), s=4, alpha=0.18, c="#d95f02", label="|tau_dip| cells")
    ax0.plot(centers, med_abs_tau_strike, lw=2.2, c="#08519c", label="median |tau_strike| by bin")
    ax0.plot(centers, med_abs_tau_dip, lw=2.2, c="#a63603", label="median |tau_dip| by bin")
    ax0.set_ylabel("Shear magnitude (MPa)")
    ax0.grid(True, alpha=0.25)
    ax0.legend(loc="upper right", ncol=2, fontsize=9)
    ax0.set_title(title)

    for row in rows:
        color = "#9ecae1" if row["dominant"] == "strike" else "#fdd0a2"
        ax1.axvspan(row["bin_start_km"], row["bin_end_km"], color=color, alpha=0.45, lw=0)
    ax1.plot(centers, strike_frac, c="black", lw=1.5)
    ax1.axhline(0.5, c="0.35", ls="--", lw=1.0)
    ax1.set_ylim(-0.02, 1.02)
    ax1.set_ylabel("Strike-dom. frac")
    ax1.set_xlabel("Along-fault coordinate (km)")
    ax1.grid(True, alpha=0.25)

    ax1.text(0.01, 0.92, "blue: strike-dominant bins", transform=ax1.transAxes, fontsize=9)
    ax1.text(0.01, 0.80, "orange: dip-dominant bins", transform=ax1.transAxes, fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot along-fault strike- vs dip-shear dominance from a fault-stress VTU."
    )
    ap.add_argument("--fault-vtu", type=Path, required=True, help="path to *_fault_stress.vtu")
    ap.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="optional summary JSON; used to recover the strike-hint azimuth",
    )
    ap.add_argument("--out-prefix", type=Path, required=True, help="output prefix for .png/.csv/.json")
    ap.add_argument(
        "--out-vtu",
        type=Path,
        default=None,
        help="optional explicit path for the augmented VTU; defaults to <out_prefix>.vtu",
    )
    ap.add_argument("--n-bins", type=int, default=80, help="number of along-fault bins")
    args = ap.parse_args()

    if args.n_bins < 2:
        raise ValueError("--n-bins must be at least 2")

    mesh = meshio.read(str(args.fault_vtu))
    tri_conn = _load_fault_triangles(mesh)
    tau_strike = _load_cell_array(mesh, "tau_strike_MPa_cell")
    tau_dip = _load_cell_array(mesh, "tau_dip_MPa_cell")
    centroids = _triangle_centroids(np.asarray(mesh.points, dtype=np.float64), tri_conn)
    strike_hat_xy, axis_source = _strike_hat_xy(args.summary_json, centroids[:, :2])

    proj_m = (centroids[:, :2] - centroids[:, :2].mean(axis=0)) @ strike_hat_xy
    along_km = (proj_m - proj_m.min()) / 1000.0
    rows, _, _ = _build_bin_rows(along_km, tau_strike, tau_dip, args.n_bins)
    if not rows:
        raise ValueError("no populated bins produced; fault mesh may be empty")

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_png = args.out_prefix.with_suffix(".png")
    out_csv = args.out_prefix.with_suffix(".csv")
    out_json = args.out_prefix.with_suffix(".json")
    out_vtu = args.out_vtu if args.out_vtu is not None else args.out_prefix.with_suffix(".vtu")

    title = (
        f"{args.fault_vtu.name}\n"
        f"axis={axis_source}, cells={tri_conn.shape[0]}, "
        f"strike-dominant={(np.abs(tau_strike) >= np.abs(tau_dip)).sum()} / {tri_conn.shape[0]}"
    )
    _write_augmented_vtu(out_vtu, mesh, along_km, tau_strike, tau_dip)
    _plot(out_png, along_km, tau_strike, tau_dip, rows, title)
    _write_csv(out_csv, rows)

    dominant = np.where(np.abs(tau_strike) >= np.abs(tau_dip), "strike", "dip")
    starts = np.array([row["bin_start_km"] for row in rows], dtype=np.float64)
    ends = np.array([row["bin_end_km"] for row in rows], dtype=np.float64)
    labels = [row["dominant"] for row in rows]
    payload = {
        "fault_vtu": str(args.fault_vtu),
        "augmented_vtu": str(out_vtu),
        "summary_json": None if args.summary_json is None else str(args.summary_json),
        "n_cells": int(tri_conn.shape[0]),
        "n_bins": int(len(rows)),
        "along_fault_range_km": [float(along_km.min()), float(along_km.max())],
        "axis_source": axis_source,
        "strike_hat_xy": [float(strike_hat_xy[0]), float(strike_hat_xy[1])],
        "cell_counts": {
            "strike_dominant": int(np.count_nonzero(dominant == "strike")),
            "dip_dominant": int(np.count_nonzero(dominant == "dip")),
        },
        "bin_counts": {
            "strike_dominant": int(sum(label == "strike" for label in labels)),
            "dip_dominant": int(sum(label == "dip" for label in labels)),
        },
        "vtu_fields": {
            "along_fault_km_cell": "triangle-centroid along-fault coordinate",
            "abs_tau_strike_MPa_cell": "|tau_strike| at triangle",
            "abs_tau_dip_MPa_cell": "|tau_dip| at triangle",
            "shear_dominance_margin_MPa_cell": "|tau_strike| - |tau_dip|",
            "shear_dominance_ratio_cell": "|tau_strike| / max(|tau_dip|, eps)",
            "shear_dominance_flag_cell": "-1=dip-dominant, +1=strike-dominant",
        },
        "dominant_intervals_km": _intervals_from_labels(starts, ends, labels),
    }
    _write_json(out_json, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
