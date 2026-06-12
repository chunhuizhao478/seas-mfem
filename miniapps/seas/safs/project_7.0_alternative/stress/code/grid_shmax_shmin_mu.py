#!/usr/bin/env python3
"""grid_shmax_shmin_mu.py — grid-search (SHmax, Shmin) within the
Hickman & Zoback (2004) Fig. 4a error bars at the ~2 km anchor depth and
report the mu_apparent field statistics per pair, at fixed psi / Sv / P_p.

Companion to `sweep_psi_max_mu.py` (which varies psi at fixed magnitudes);
this varies the horizontal principal magnitudes at fixed orientation.
The figure of merit follows the same convention: CELL-level field maximum
of mu_apparent = |tau| / sigma_n_eff, plus area-weighted quantiles and
critically-stressed area fractions.

Analytic context (used as a cross-check, printed per pair): for a
constant tensor the field max is bounded by the orientation bound

    mu* = (s1' - s3') / (2 sqrt(s1' s3'))

with s1' = SHmax - P_p and s3' = min(Shmin, Sv) - P_p; the bound is
attained whenever the meshed fault contains a critically oriented patch.

Outputs (in --out-dir):
    grid_shmax_shmin_mu.csv    one row per (SHmax, Shmin)
    grid_shmax_shmin_mu.json   parameter echo + argmin rows per metric
    grid_shmax_shmin_mu.png    heatmaps of mu_max and mu_median

Usage:
    conda activate pythonenv
    python grid_shmax_shmin_mu.py \
        --fault-vtu ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
        --out-dir   ../results/grid_shmax_shmin_triq_psi69_pp20 \
        --SHmax-range 100 130 2.5 --Shmin-range 40 60 2
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from project_to_fault_stress import (
    DEFAULT_FAULT_NAME,
    build_fault_basis,
    load_fault_mesh,
    project_stress_onto_fault,
)
from sweep_psi_max_mu import _area_weighted_quantiles


def orientation_bound(SHmax, Shmin, Sv, P_p):
    """Analytic max of |tau|/sigma_n_eff over ALL plane orientations."""
    s1 = SHmax - P_p
    s3 = min(Shmin, Sv) - P_p
    if s3 <= 0:
        return float("inf")
    return (s1 - s3) / (2.0 * np.sqrt(s1 * s3))


def grid_search(
    fault_geom,
    shmax_values,
    shmin_values,
    *,
    psi: float,
    strike_hint_az: float,
    Sv: float,
    P_p: float,
) -> list:
    """Evaluate the projection for every (SHmax, Shmin) pair."""
    az = (strike_hint_az + psi) % 360.0
    areas = fault_geom.areas
    total_area = float(areas.sum())
    rows = []
    for SHmax in shmax_values:
        for Shmin in shmin_values:
            if Shmin > SHmax:
                raise ValueError(
                    f"Shmin ({Shmin}) > SHmax ({SHmax}) in grid")
            resolved = project_stress_onto_fault(
                fault_geom,
                SHmax=float(SHmax), Shmin=float(Shmin),
                Sv=Sv, P_p=P_p, SHmax_az_deg=az,
            )
            mu = resolved["mu_apparent"]
            n_nan = int(np.isnan(mu).sum())
            if n_nan == mu.shape[0]:
                raise ValueError(
                    f"SHmax={SHmax} Shmin={Shmin}: sigma_n_eff <= 0 "
                    f"over the whole field")
            med, p90, p99 = _area_weighted_quantiles(
                mu, areas, (0.5, 0.90, 0.99))
            rows.append({
                "SHmax_MPa": float(SHmax),
                "Shmin_MPa": float(Shmin),
                "mu_max": float(np.nanmax(mu)),
                "mu_bound_analytic": orientation_bound(
                    SHmax, Shmin, Sv, P_p),
                "mu_p99": p99,
                "mu_p90": p90,
                "mu_median": med,
                "area_pct_mu_gt_060": float(
                    areas[mu > 0.60].sum() / total_area * 100.0),
                "area_pct_mu_gt_055": float(
                    areas[mu > 0.55].sum() / total_area * 100.0),
                "area_pct_mu_gt_050": float(
                    areas[mu > 0.50].sum() / total_area * 100.0),
                "tau_mag_max_MPa": float(np.max(resolved["tau_magnitude"])),
                "sigma_n_eff_min_MPa": float(
                    np.min(resolved["sigma_n_eff"])),
                "mu_nan_cells": n_nan,
            })
    return rows


def _write_outputs(rows, out_dir: Path, params: dict,
                   shmax_values, shmin_values) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "grid_shmax_shmin_mu.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {csv_path}")

    argmins = {}
    for key in ("mu_max", "mu_p99", "mu_p90", "mu_median",
                "area_pct_mu_gt_055", "area_pct_mu_gt_050"):
        r = min(rows, key=lambda x: x[key])
        argmins[key] = {"SHmax_MPa": r["SHmax_MPa"],
                        "Shmin_MPa": r["Shmin_MPa"],
                        "value": r[key]}
    best = min(rows, key=lambda r: r["mu_max"])

    json_path = out_dir / "grid_shmax_shmin_mu.json"
    with open(json_path, "w") as f:
        json.dump({"params": params, "n_pairs": len(rows),
                   "argmin_per_metric": argmins,
                   "best_by_mu_max": best}, f, indent=2)
    print(f"wrote {json_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        nx, ny = len(shmin_values), len(shmax_values)
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
        for ax, key, title in (
            (axes[0], "mu_max", "field-max mu_apparent"),
            (axes[1], "mu_median", "median (area-wt) mu_apparent"),
        ):
            grid = np.array([r[key] for r in rows]).reshape(ny, nx)
            im = ax.imshow(
                grid, origin="lower", aspect="auto", cmap="viridis",
                extent=[shmin_values[0], shmin_values[-1],
                        shmax_values[0], shmax_values[-1]])
            kbest = int(np.argmin([r[key] for r in rows]))
            ax.plot(rows[kbest]["Shmin_MPa"], rows[kbest]["SHmax_MPa"],
                    "r*", ms=16)
            ax.set_xlabel("Shmin [MPa]")
            ax.set_ylabel("SHmax [MPa]")
            ax.set_title(f"{title} (red star = argmin)")
            fig.colorbar(im, ax=ax)
        fig.suptitle(
            f"psi={params['psi_deg']:g} deg, Sv={params['Sv_MPa']:g}, "
            f"P_p={params['P_p_MPa']:g} MPa (constant with depth)",
            fontsize=10)
        fig.tight_layout()
        png_path = out_dir / "grid_shmax_shmin_mu.png"
        fig.savefig(png_path, dpi=160)
        print(f"wrote {png_path}")
    except ImportError as exc:
        print(f"matplotlib unavailable ({exc}); skipping PNG",
              file=sys.stderr)

    return best


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Grid-search (SHmax, Shmin) for minimum mu_apparent "
                    "field at fixed psi / Sv / P_p.")
    ap.add_argument("--fault-vtu", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--SHmax-range", type=float, nargs=3,
                    default=[100.0, 130.0, 2.5],
                    metavar=("MIN", "MAX", "STEP"))
    ap.add_argument("--Shmin-range", type=float, nargs=3,
                    default=[40.0, 60.0, 2.0],
                    metavar=("MIN", "MAX", "STEP"))
    ap.add_argument("--psi", type=float, default=69.0,
                    help="angle SHmax-to-strike-hint (deg, default 69)")
    ap.add_argument("--Sv", type=float, default=50.0)
    ap.add_argument("--P_p", type=float, default=20.0)
    ap.add_argument("--strike-hint-az", type=float, default=314.0)
    ap.add_argument("--fault-name", type=str, default=DEFAULT_FAULT_NAME)
    args = ap.parse_args(argv)

    def _expand(rng, name):
        lo, hi, step = rng
        if step <= 0 or hi < lo:
            ap.error(f"--{name}-range must be MIN MAX STEP with "
                     f"STEP > 0 and MAX >= MIN; got {rng}")
        n = int(round((hi - lo) / step)) + 1
        vals = lo + step * np.arange(n)
        return vals[vals <= hi + 1e-9]

    shmax_values = _expand(args.SHmax_range, "SHmax")
    shmin_values = _expand(args.Shmin_range, "Shmin")

    print(f"Reading {args.fault_vtu} ...")
    mesh = load_fault_mesh(args.fault_vtu)
    geom = build_fault_basis(
        mesh,
        fault_phys_name=args.fault_name,
        fault_strike_azimuth_hint_deg=args.strike_hint_az,
    )
    print(f"  {geom.centroids.shape[0]} triangles, "
          f"{geom.n_degenerate} degenerate")
    print(f"  grid: {len(shmax_values)} SHmax x {len(shmin_values)} "
          f"Shmin = {len(shmax_values) * len(shmin_values)} pairs")

    rows = grid_search(
        geom, shmax_values, shmin_values,
        psi=args.psi, strike_hint_az=args.strike_hint_az,
        Sv=args.Sv, P_p=args.P_p,
    )

    params = {
        "fault_vtu": str(args.fault_vtu),
        "psi_deg": args.psi,
        "SHmax_az_deg": (args.strike_hint_az + args.psi) % 360.0,
        "strike_hint_az_deg": args.strike_hint_az,
        "Sv_MPa": args.Sv,
        "P_p_MPa": args.P_p,
        "depth_model": "constant",
        "SHmax_range": list(args.SHmax_range),
        "Shmin_range": list(args.Shmin_range),
    }
    best = _write_outputs(rows, args.out_dir, params,
                          shmax_values, shmin_values)

    print(f"\n{'SHmax':>6} {'Shmin':>6} {'mu_max':>7} {'bound':>7} "
          f"{'p99':>7} {'p90':>7} {'median':>7} {'A%>.55':>7} "
          f"{'A%>.50':>7}")
    for r in rows:
        mark = "  <-- min(mu_max)" if r is best else ""
        print(f"{r['SHmax_MPa']:6.1f} {r['Shmin_MPa']:6.1f} "
              f"{r['mu_max']:7.4f} {r['mu_bound_analytic']:7.4f} "
              f"{r['mu_p99']:7.4f} {r['mu_p90']:7.4f} "
              f"{r['mu_median']:7.4f} "
              f"{r['area_pct_mu_gt_055']:7.3f} "
              f"{r['area_pct_mu_gt_050']:7.3f}{mark}")

    print(f"\nBEST by mu_max: SHmax = {best['SHmax_MPa']:g}, "
          f"Shmin = {best['Shmin_MPa']:g} "
          f"(field-max mu = {best['mu_max']:.4f}, "
          f"median = {best['mu_median']:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
