#!/usr/bin/env python3
"""
BP5-QD 2x2 strike-component comparison: MFEM-SEAS vs Tandem (benchmark).

Panels (clockwise from top-left):
    [0,0] Slip rate V along strike (m/s, log scale)
    [0,1] Slip along strike (m)
    [1,0] Shear stress along strike (MPa)
    [1,1] State variable log10(theta) (s)

MFEM data:    /Users/chunhuizhao/Downloads/seas-mfem/results_phase7_production_job7650438
Tandem data:  ./benchmark_data/scec_tandem  (Tandem p4 reference)

Usage:
    python plot_bp5_strike_2x2.py                # all 10 stations, save PNGs
    python plot_bp5_strike_2x2.py --show         # display instead of save
    python plot_bp5_strike_2x2.py --stations 1 3 # subset
"""

import argparse
import os
import sys

import numpy as np

SCEC_STATIONS = [
    ("fltst_strk-36dp+00", -36, 0),
    ("fltst_strk-16dp+00", -16, 0),
    ("fltst_strk+00dp+00", 0, 0),
    ("fltst_strk+16dp+00", 16, 0),
    ("fltst_strk+36dp+00", 36, 0),
    ("fltst_strk-24dp+10", -24, 10),
    ("fltst_strk-16dp+10", -16, 10),
    ("fltst_strk+00dp+10", 0, 10),
    ("fltst_strk+16dp+10", 16, 10),
    ("fltst_strk+00dp+22", 0, 22),
]

SEC_PER_YR = 3.15576e7

DEFAULT_MFEM_DIR = (
    "/Users/chunhuizhao/Downloads/seas-mfem/results_phase7_production_job7650438"
)
DEFAULT_MFEM_PREFIX = "bp5_phase7_prod_job7650438"
DEFAULT_TANDEM_SUBDIR = os.path.join("benchmark_data", "golden_tandem_p1_1000m")
DEFAULT_TANDEM_ORDER = 1


def load_bp5_file(filepath):
    """Load 8-column BP5 SCEC text file: t, slip_s, slip_d, log10V_s, log10V_d,
    tau_s, tau_d, log10(theta)."""
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    rows.append([float(x) for x in parts[:8]])
                except ValueError:
                    continue
    if not rows:
        return None
    arr = np.array(rows)
    return {
        "time_yr": arr[:, 0] / SEC_PER_YR,
        "slip_strike": arr[:, 1],
        "V_strike": 10.0 ** arr[:, 3],
        "tau_strike": arr[:, 5],
        "log10_theta": arr[:, 7],
    }


def coord_str(v):
    return str(int(v))


def mfem_path(mfem_dir, prefix, station):
    return os.path.join(mfem_dir, f"{prefix}_{station}.txt")


def tandem_path(tandem_dir, x2_km, x3_km, order):
    return os.path.join(
        tandem_dir,
        f"bp5qd_tandem_p{order}_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt",
    )


def plot_station_2x2(mfem, tandem, station, x2_km, x3_km, save_path=None):
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
            "axes.linewidth": 1.1,
            "axes.grid": True,
            "grid.alpha": 0.3,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"BP5-QD strike-component comparison — {station}  "
        f"(x2={x2_km} km, x3={x3_km} km)",
        fontsize=15,
        fontweight="bold",
    )

    panels = [
        (axes[0, 0], "V_strike", "Slip rate $V_{\\mathrm{strike}}$ (m/s)", True),
        (axes[0, 1], "slip_strike", "Slip $s_{\\mathrm{strike}}$ (m)", False),
        (axes[1, 0], "tau_strike", "Shear stress $\\tau_{\\mathrm{strike}}$ (MPa)", False),
        (axes[1, 1], "log10_theta", "State $\\log_{10}\\theta$ (s)", False),
    ]

    series = []
    if tandem is not None:
        series.append(("Tandem (benchmark)", tandem, "#000000", "--", 1.6))
    if mfem is not None:
        series.append(("MFEM-SEAS", mfem, "#d62728", "-", 1.4))

    legend_panel_idx = 1  # only top-right panel (slip strike) carries the legend
    for idx, (ax, key, ylabel, use_log) in enumerate(panels):
        handles, labels = [], []
        for label, data, color, ls, lw in series:
            (line,) = ax.plot(
                data["time_yr"],
                data[key],
                ls,
                color=color,
                label=label,
                linewidth=lw,
                alpha=0.9,
            )
            handles.append(line)
            labels.append(label)
        ax.set_xlabel("Time (years)")
        ax.set_ylabel(ylabel)
        if use_log:
            ax.set_yscale("log")
        if idx == legend_panel_idx:
            ax.legend(
                handles,
                labels,
                loc="lower right",
                frameon=True,
                framealpha=0.95,
                fontsize=14,
                handlelength=3.0,
                borderpad=0.8,
            )
        ax.tick_params(direction="in", length=4, top=True, right=True)

    plt.tight_layout(rect=(0, 0, 1, 0.96))
    if save_path:
        plt.savefig(save_path, dpi=180, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="2x2 BP5 strike-component comparison: MFEM-SEAS vs Tandem (benchmark)."
    )
    parser.add_argument("--mfem-dir", default=DEFAULT_MFEM_DIR)
    parser.add_argument("--mfem-prefix", default=DEFAULT_MFEM_PREFIX)
    parser.add_argument(
        "--tandem-dir",
        default=None,
        help=f"Default: <script_dir>/{DEFAULT_TANDEM_SUBDIR}",
    )
    parser.add_argument("--tandem-order", type=int, default=DEFAULT_TANDEM_ORDER)
    parser.add_argument("--stations", nargs="+", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "plots_strike_2x2_phase7_production_job7650438",
        ),
    )
    parser.add_argument("--show", action="store_true", help="Display instead of save")
    args = parser.parse_args()

    import matplotlib

    if not args.show:
        matplotlib.use("Agg")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    tandem_dir = args.tandem_dir or os.path.join(script_dir, DEFAULT_TANDEM_SUBDIR)

    stations = (
        [SCEC_STATIONS[i - 1] for i in args.stations if 1 <= i <= len(SCEC_STATIONS)]
        if args.stations
        else SCEC_STATIONS
    )

    if not args.show:
        os.makedirs(args.output_dir, exist_ok=True)

    print("BP5-QD strike-component 2x2 comparison")
    print(f"  MFEM-SEAS:        {args.mfem_dir} (prefix: {args.mfem_prefix})")
    print(f"  Tandem benchmark: {tandem_dir} (p{args.tandem_order})")
    print(f"  Output dir:       {args.output_dir if not args.show else '(display)'}")
    print(f"  Stations:         {len(stations)}")
    print()

    n_plotted = 0
    for station, x2_km, x3_km in stations:
        mfem_p = mfem_path(args.mfem_dir, args.mfem_prefix, station)
        tand_p = tandem_path(tandem_dir, x2_km, x3_km, args.tandem_order)

        mfem = load_bp5_file(mfem_p) if os.path.exists(mfem_p) else None
        tandem = load_bp5_file(tand_p) if os.path.exists(tand_p) else None

        if mfem is None and tandem is None:
            print(f"  {station}: no data (missing both MFEM and Tandem files)")
            continue
        if mfem is None:
            print(f"  {station}: MFEM file missing -> {mfem_p}")
        if tandem is None:
            print(f"  {station}: Tandem file missing -> {tand_p}")

        save_path = (
            None
            if args.show
            else os.path.join(args.output_dir, f"bp5_strike_2x2_{station}.png")
        )
        plot_station_2x2(mfem, tandem, station, x2_km, x3_km, save_path=save_path)
        n_plotted += 1

    print(f"\nPlotted {n_plotted} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
