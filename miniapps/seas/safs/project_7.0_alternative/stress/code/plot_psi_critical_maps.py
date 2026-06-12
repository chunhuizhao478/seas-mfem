#!/usr/bin/env python3
"""plot_psi_critical_maps.py — plan-view maps of the critically-stressed
fault area (mu_apparent thresholds) at selected Ψ angles.

Companion to `sweep_psi_max_mu.py`: that script shows HOW MUCH fault area
is critically stressed per Ψ; this one shows WHERE it sits. One panel per
requested Ψ; cells are coloured by mu_apparent band:

    grey   mu <= 0.55
    orange 0.55 < mu <= 0.60
    red    mu > 0.60   (near the orientation bound)

Usage:
    conda activate pythonenv
    python plot_psi_critical_maps.py \
        --fault-vtu ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
        --out ../results/psi_sweep_triq_shmax120_pp20/psi_critical_maps.png \
        --psi 52 69 72.5 83.5
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from project_to_fault_stress import (
    build_fault_basis,
    load_fault_mesh,
    project_stress_onto_fault,
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Plan-view maps of critically-stressed fault area "
                    "at selected psi angles.")
    ap.add_argument("--fault-vtu", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--psi", type=float, nargs="+",
                    default=[52.0, 69.0, 72.5, 83.5])
    ap.add_argument("--SHmax", type=float, default=120.0)
    ap.add_argument("--Shmin", type=float, default=55.0)
    ap.add_argument("--Sv", type=float, default=50.0)
    ap.add_argument("--P_p", type=float, default=20.0)
    ap.add_argument("--strike-hint-az", type=float, default=314.0)
    args = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print(f"Reading {args.fault_vtu} ...")
    mesh = load_fault_mesh(args.fault_vtu)
    geom = build_fault_basis(
        mesh, fault_strike_azimuth_hint_deg=args.strike_hint_az)
    c = geom.centroids
    areas = geom.areas
    total_area = float(areas.sum())

    n = len(args.psi)
    ncol = min(n, 2)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(7.0 * ncol, 5.6 * nrow),
        sharex=True, sharey=True, squeeze=False)

    for k, psi in enumerate(args.psi):
        ax = axes[k // ncol][k % ncol]
        az = (args.strike_hint_az + psi) % 360.0
        resolved = project_stress_onto_fault(
            geom, SHmax=args.SHmax, Shmin=args.Shmin, Sv=args.Sv,
            P_p=args.P_p, SHmax_az_deg=az)
        mu = resolved["mu_apparent"]

        lo = mu <= 0.55
        mid = (mu > 0.55) & (mu <= 0.60)
        hi = mu > 0.60
        pct_hi = areas[hi].sum() / total_area * 100.0
        pct_mid = areas[mid].sum() / total_area * 100.0

        ax.scatter(c[lo, 0] / 1e3, c[lo, 1] / 1e3, s=1.0,
                   c="0.8", rasterized=True)
        ax.scatter(c[mid, 0] / 1e3, c[mid, 1] / 1e3, s=1.2,
                   c="orange", rasterized=True)
        ax.scatter(c[hi, 0] / 1e3, c[hi, 1] / 1e3, s=1.4,
                   c="red", rasterized=True)
        ax.set_title(
            f"psi = {psi:g} deg (SHmax az {az:g})   "
            f"red mu>0.60: {pct_hi:.2f}%   orange 0.55-0.60: "
            f"{pct_mid:.1f}%", fontsize=10)
        ax.set_xlabel("UTM easting [km]")
        ax.set_ylabel("UTM northing [km]")
        ax.set_aspect("equal")
        ax.grid(alpha=0.25)

    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)

    fig.suptitle(
        f"Critically-stressed fault area, SHmax={args.SHmax:g} "
        f"Shmin={args.Shmin:g} Sv={args.Sv:g} P_p={args.P_p:g} MPa",
        fontsize=11)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=170)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
