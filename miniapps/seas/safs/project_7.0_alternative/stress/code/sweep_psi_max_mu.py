#!/usr/bin/env python3
"""sweep_psi_max_mu.py — sweep Ψ (angle between SHmax and the SAF strike
hint) and report the field-maximum apparent friction per angle.

For each Ψ in the sweep the SHmax azimuth is placed at

    SHmax_az = (strike_hint_az + Ψ) mod 360

(the same branch used by the psi37/az351 and psi69/az23 runs: 314 + 37 =
351, 314 + 69 = 23 mod 360), the constant-depth regional tensor is rebuilt,
and the per-cell `resolve_traction` projection is re-evaluated on the SAME
per-triangle Tandem basis (the basis depends only on the mesh + strike
hint, not on Ψ, so it is built once).

The figure of merit is the CELL-level field maximum of

    mu_apparent = tau_magnitude / sigma_n_eff ,  sigma_n_eff = sigma_n - P_p

which is the same quantity `write_summary_json` reports as
`stats.mu_apparent.max` — so any sweep row is directly comparable to a
full projection run at that angle.

Outputs (in --out-dir):
    sweep_psi_max_mu.csv   one row per Ψ
    sweep_psi_max_mu.json  parameter echo + the argmin row
    sweep_psi_max_mu.png   max/median mu vs Ψ (skipped if matplotlib
                           is unavailable)

Usage:
    conda activate pythonenv
    python sweep_psi_max_mu.py \
        --fault-vtu ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
        --out-dir   ../results/psi_sweep_triq_shmax120_pp20 \
        --psi-min 52 --psi-max 86 --psi-step 0.5
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


def _area_weighted_quantiles(values, areas, qs):
    """Area-weighted quantiles of a per-cell field (NaN cells dropped)."""
    finite = np.isfinite(values)
    v = values[finite]
    w = areas[finite]
    order = np.argsort(v)
    v = v[order]
    cw = np.cumsum(w[order]) / w[order].sum()
    return [float(np.interp(q, cw, v)) for q in qs]


def sweep_psi(
    fault_geom,
    psi_values: np.ndarray,
    *,
    strike_hint_az: float,
    SHmax: float,
    Shmin: float,
    Sv: float,
    P_p: float,
) -> list:
    """Evaluate the on-fault projection for every Ψ and collect stats.

    Returns a list of per-Ψ dicts (CSV-ready, plain floats). Quantiles
    and tail fractions are AREA-weighted so they measure fault surface,
    not triangle count.
    """
    psi_values = np.asarray(psi_values, dtype=np.float64)
    if psi_values.ndim != 1 or psi_values.size == 0:
        raise ValueError(
            f"psi_values must be a non-empty 1-D array; got shape "
            f"{psi_values.shape}"
        )
    areas = fault_geom.areas
    total_area = float(areas.sum())
    rows = []
    for psi in psi_values:
        az = (strike_hint_az + psi) % 360.0
        resolved = project_stress_onto_fault(
            fault_geom,
            SHmax=SHmax,
            Shmin=Shmin,
            Sv=Sv,
            P_p=P_p,
            SHmax_az_deg=az,
        )
        mu = resolved["mu_apparent"]
        n_nan = int(np.isnan(mu).sum())
        if n_nan == mu.shape[0]:
            raise ValueError(
                f"psi={psi}: every cell has sigma_n_eff <= 0; "
                f"mu_apparent undefined over the whole field"
            )
        med, p90, p99, p999 = _area_weighted_quantiles(
            mu, areas, (0.5, 0.90, 0.99, 0.999))
        i_max = int(np.nanargmax(mu))
        c = fault_geom.centroids[i_max]
        rows.append({
            "psi_deg": float(psi),
            "SHmax_az_deg": float(az),
            "mu_max": float(np.nanmax(mu)),
            "mu_p999": p999,
            "mu_p99": p99,
            "mu_p90": p90,
            "mu_median": med,
            "area_pct_mu_gt_060": float(
                areas[mu > 0.60].sum() / total_area * 100.0),
            "area_pct_mu_gt_055": float(
                areas[mu > 0.55].sum() / total_area * 100.0),
            "area_pct_mu_gt_050": float(
                areas[mu > 0.50].sum() / total_area * 100.0),
            "mu_nan_cells": n_nan,
            "sigma_n_eff_min_MPa": float(np.min(resolved["sigma_n_eff"])),
            "tau_mag_max_MPa": float(np.max(resolved["tau_magnitude"])),
            "tau_strike_median_MPa": float(
                np.median(resolved["tau_strike"])
            ),
            "mu_max_cell_x_m": float(c[0]),
            "mu_max_cell_y_m": float(c[1]),
            "mu_max_cell_z_m": float(c[2]),
        })
    return rows


def _write_outputs(rows: list, out_dir: Path, params: dict) -> dict:
    """Write CSV + JSON (+ optional PNG); return the argmin row."""
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "sweep_psi_max_mu.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {csv_path}")

    best = min(rows, key=lambda r: r["mu_max"])
    json_path = out_dir / "sweep_psi_max_mu.json"
    with open(json_path, "w") as f:
        json.dump({"params": params, "n_psi": len(rows), "best": best},
                  f, indent=2)
    print(f"wrote {json_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        psi = [r["psi_deg"] for r in rows]
        fig, (ax, ax2) = plt.subplots(
            1, 2, figsize=(12.5, 4.8), sharex=True)

        for key, style, label in (
            ("mu_max", "k-", "max (cell)"),
            ("mu_p999", "v-", "p99.9 (area-wt)"),
            ("mu_p99", "^-", "p99 (area-wt)"),
            ("mu_p90", "o-", "p90 (area-wt)"),
            ("mu_median", "s--", "median (area-wt)"),
        ):
            ax.plot(psi, [r[key] for r in rows], style, ms=3, label=label)
        ax.set_xlabel("psi = angle(SHmax, SAF strike hint)  [deg]")
        ax.set_ylabel("mu_apparent = |tau| / sigma_n_eff")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_title("mu_apparent quantiles vs psi")

        for key, style, label in (
            ("area_pct_mu_gt_060", "o-", "mu > 0.60"),
            ("area_pct_mu_gt_055", "^-", "mu > 0.55"),
            ("area_pct_mu_gt_050", "s--", "mu > 0.50"),
        ):
            vals = [r[key] for r in rows]
            ax2.plot(psi, vals, style, ms=3, label=label)
            k = int(np.argmin(vals))
            ax2.plot(psi[k], vals[k], "r*", ms=12)
        ax2.set_xlabel("psi = angle(SHmax, SAF strike hint)  [deg]")
        ax2.set_ylabel("% of fault area above threshold")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8)
        ax2.set_title("critically-stressed area fraction vs psi "
                      "(red star = argmin)")

        fig.suptitle(
            f"SHmax={params['SHmax_MPa']:g} Shmin={params['Shmin_MPa']:g} "
            f"Sv={params['Sv_MPa']:g} P_p={params['P_p_MPa']:g} MPa "
            f"(constant with depth)", fontsize=10)
        fig.tight_layout()
        png_path = out_dir / "sweep_psi_max_mu.png"
        fig.savefig(png_path, dpi=160)
        print(f"wrote {png_path}")
    except ImportError as exc:
        print(f"matplotlib unavailable ({exc}); skipping PNG",
              file=sys.stderr)

    return best


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Sweep psi = angle(SHmax, SAF strike) and find the "
                    "angle minimising the field-max apparent friction.")
    ap.add_argument("--fault-vtu", type=Path, required=True,
                    help="path to a *_fault.vtu (msh_to_vtu.py output)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--psi-min", type=float, default=52.0)
    ap.add_argument("--psi-max", type=float, default=86.0)
    ap.add_argument("--psi-step", type=float, default=0.5)
    ap.add_argument("--SHmax", type=float, default=120.0)
    ap.add_argument("--Shmin", type=float, default=55.0)
    ap.add_argument("--Sv", type=float, default=50.0)
    ap.add_argument("--P_p", type=float, default=20.0)
    ap.add_argument("--strike-hint-az", type=float, default=314.0)
    ap.add_argument("--fault-name", type=str, default=DEFAULT_FAULT_NAME)
    args = ap.parse_args(argv)

    if args.psi_step <= 0:
        ap.error(f"--psi-step must be > 0; got {args.psi_step}")
    if args.psi_max < args.psi_min:
        ap.error(f"--psi-max ({args.psi_max}) < --psi-min ({args.psi_min})")

    n = int(round((args.psi_max - args.psi_min) / args.psi_step)) + 1
    psi_values = args.psi_min + args.psi_step * np.arange(n)
    psi_values = psi_values[psi_values <= args.psi_max + 1e-9]

    print(f"Reading {args.fault_vtu} ...")
    mesh = load_fault_mesh(args.fault_vtu)
    geom = build_fault_basis(
        mesh,
        fault_phys_name=args.fault_name,
        fault_strike_azimuth_hint_deg=args.strike_hint_az,
    )
    print(f"  {geom.centroids.shape[0]} triangles, "
          f"{geom.n_degenerate} degenerate")

    rows = sweep_psi(
        geom, psi_values,
        strike_hint_az=args.strike_hint_az,
        SHmax=args.SHmax, Shmin=args.Shmin, Sv=args.Sv, P_p=args.P_p,
    )

    params = {
        "fault_vtu": str(args.fault_vtu),
        "strike_hint_az_deg": args.strike_hint_az,
        "SHmax_MPa": args.SHmax,
        "Shmin_MPa": args.Shmin,
        "Sv_MPa": args.Sv,
        "P_p_MPa": args.P_p,
        "depth_model": "constant",
        "psi_min_deg": args.psi_min,
        "psi_max_deg": args.psi_max,
        "psi_step_deg": args.psi_step,
    }
    best = _write_outputs(rows, args.out_dir, params)

    print(f"\n{'psi':>6} {'az':>6} {'mu_max':>7} {'p99.9':>7} {'p99':>7} "
          f"{'p90':>7} {'median':>7} {'A%>.60':>7} {'A%>.55':>7}")
    for r in rows:
        mark = "  <-- min(mu_max)" if r is best else ""
        print(f"{r['psi_deg']:6.1f} {r['SHmax_az_deg']:6.1f} "
              f"{r['mu_max']:7.4f} {r['mu_p999']:7.4f} "
              f"{r['mu_p99']:7.4f} {r['mu_p90']:7.4f} "
              f"{r['mu_median']:7.4f} "
              f"{r['area_pct_mu_gt_060']:7.3f} "
              f"{r['area_pct_mu_gt_055']:7.3f}{mark}")

    print(f"\nBEST by mu_max: psi = {best['psi_deg']:g} deg "
          f"(SHmax_az = {best['SHmax_az_deg']:g} deg cw-from-N), "
          f"field-max mu_apparent = {best['mu_max']:.4f}")
    for key in ("mu_p999", "mu_p99", "mu_p90", "mu_median",
                "area_pct_mu_gt_060", "area_pct_mu_gt_055"):
        r = min(rows, key=lambda x: x[key])
        print(f"  argmin {key:22s}: psi = {r['psi_deg']:5.1f}  "
              f"(value {r[key]:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
