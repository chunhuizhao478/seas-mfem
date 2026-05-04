#!/usr/bin/env python3
"""
TPV104 2x2 strike-component comparison: MFEM-SEAS vs DRDG3D (benchmark).

Panels (clockwise from top-left):
    [0,0] Slip rate V_strike (m/s, log scale)
    [0,1] Slip s_strike (m)
    [1,0] Shear stress tau_strike (MPa)
    [1,1] State variable psi (dimensionless)

MFEM data:    <gold>/results/tpv104_<prefix>_station_x2_<x2>_x3_<x3>.dat   (stresses in Pa)
DRDG3D data:  ./benchmark_data/DRDG3D/tpv104_drdg3d_x2_<x2>_x3_<x3>.txt    (stresses in MPa)

Both files share the SCEC TPV104 9-column layout:
    t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  psi

Usage:
    python plot_tpv104_strike_2x2.py                       # all 9 stations, save PNGs
    python plot_tpv104_strike_2x2.py --show                # display instead of save
    python plot_tpv104_strike_2x2.py --stations 1 4 6      # subset (1-based indexing)
"""

import argparse
import os
import sys

import numpy as np

# (label, x2_km, x3_km) -- TPV104 SCEC fault stations.
TPV104_STATIONS = [
    ("x2_-12_x3_3",   -12, 3),
    ("x2_-12_x3_12",  -12, 12),
    ("x2_-9_x3_7.5",   -9, 7.5),
    ("x2_0_x3_3",       0, 3),
    ("x2_0_x3_7.5",     0, 7.5),
    ("x2_0_x3_12",      0, 12),
    ("x2_9_x3_7.5",     9, 7.5),
    ("x2_12_x3_3",     12, 3),
    ("x2_12_x3_12",    12, 12),
]

DEFAULT_MFEM_DIR = (
    "/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tpv104/gold/"
    "results_nonsymmetric_mesh_mixedflux_none_O2_job7680536/results"
)
DEFAULT_MFEM_PREFIX = "tpv104_mfnone_p1_O2"
DEFAULT_BENCHMARK_SUBDIR = os.path.join("benchmark_data", "DRDG3D")
DEFAULT_BENCHMARK_PREFIX = "tpv104_drdg3d"


def coord_str(v):
    """Format station coordinate: integers as '12', halves as '7.5'."""
    return str(int(v)) if float(v).is_integer() else str(v)


def load_tpv104_file(filepath, stress_in_pa):
    """Load 9-column SCEC TPV104 text file.

    Columns: t, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate,
             v-shear-stress, n-stress, psi.

    `stress_in_pa=True` for MFEM-SEAS output (Pa); converts to MPa.
    `stress_in_pa=False` for DRDG3D / SCEC benchmark output (already MPa).
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("t "):
                continue
            parts = line.split()
            if len(parts) >= 9:
                try:
                    rows.append([float(x) for x in parts[:9]])
                except ValueError:
                    continue
    if not rows:
        return None
    arr = np.array(rows)
    stress_scale = 1e-6 if stress_in_pa else 1.0
    return {
        "time_s": arr[:, 0],
        "slip_strike": arr[:, 1],
        "V_strike": np.abs(arr[:, 2]),
        "tau_strike": arr[:, 3] * stress_scale,
        "psi": arr[:, 8],
    }


def mfem_path(mfem_dir, prefix, station_label):
    return os.path.join(mfem_dir, f"{prefix}_station_{station_label}.dat")


def benchmark_path(bench_dir, prefix, x2_km, x3_km):
    return os.path.join(
        bench_dir, f"{prefix}_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt"
    )


def plot_station_2x2(mfem, bench, station_label, x2_km, x3_km, save_path=None):
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
        f"TPV104 strike-component comparison — {station_label}  "
        f"(x2={x2_km} km, x3={x3_km} km)",
        fontsize=15,
        fontweight="bold",
    )

    panels = [
        (axes[0, 0], "V_strike", "Slip Rate $V_{\\mathrm{strike}}$ (m/s)", False),
        (axes[0, 1], "slip_strike", "Slip Strike (m)", False),
        (axes[1, 0], "tau_strike", "Shear Stress $\\tau_{\\mathrm{strike}}$ (MPa)", False),
        (axes[1, 1], "psi", "State Variable $\\psi$", False),
    ]

    series = []
    if bench is not None:
        series.append(("DRDG3D (benchmark)", bench, "#000000", "--", 1.6))
    if mfem is not None:
        series.append(("MFEM-SEAS", mfem, "#d62728", "-", 1.4))

    legend_panel_idx = 1  # top-right panel carries the legend
    for idx, (ax, key, ylabel, use_log) in enumerate(panels):
        handles, labels = [], []
        for label, data, color, ls, lw in series:
            (line,) = ax.plot(
                data["time_s"],
                data[key],
                ls,
                color=color,
                label=label,
                linewidth=lw,
                alpha=0.9,
            )
            handles.append(line)
            labels.append(label)
        ax.set_xlabel("Time (s)")
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
        description="2x2 TPV104 strike-component comparison: MFEM-SEAS vs DRDG3D."
    )
    parser.add_argument("--mfem-dir", default=DEFAULT_MFEM_DIR)
    parser.add_argument("--mfem-prefix", default=DEFAULT_MFEM_PREFIX)
    parser.add_argument(
        "--benchmark-dir",
        default=None,
        help=f"Default: <script_dir>/{DEFAULT_BENCHMARK_SUBDIR}",
    )
    parser.add_argument("--benchmark-prefix", default=DEFAULT_BENCHMARK_PREFIX)
    parser.add_argument("--stations", nargs="+", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "plots_strike_2x2_drdg3d_mixedflux_none_O2_job7680536",
        ),
    )
    parser.add_argument("--show", action="store_true", help="Display instead of save")
    args = parser.parse_args()

    import matplotlib

    if not args.show:
        matplotlib.use("Agg")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    bench_dir = args.benchmark_dir or os.path.join(script_dir, DEFAULT_BENCHMARK_SUBDIR)

    stations = (
        [TPV104_STATIONS[i - 1] for i in args.stations if 1 <= i <= len(TPV104_STATIONS)]
        if args.stations
        else TPV104_STATIONS
    )

    if not args.show:
        os.makedirs(args.output_dir, exist_ok=True)

    print("TPV104 strike-component 2x2 comparison")
    print(f"  MFEM-SEAS:        {args.mfem_dir} (prefix: {args.mfem_prefix})")
    print(f"  DRDG3D benchmark: {bench_dir} (prefix: {args.benchmark_prefix})")
    print(f"  Output dir:       {args.output_dir if not args.show else '(display)'}")
    print(f"  Stations:         {len(stations)}")
    print()

    n_plotted = 0
    for station_label, x2_km, x3_km in stations:
        mfem_p = mfem_path(args.mfem_dir, args.mfem_prefix, station_label)
        bench_p = benchmark_path(bench_dir, args.benchmark_prefix, x2_km, x3_km)

        mfem = load_tpv104_file(mfem_p, stress_in_pa=True) if os.path.exists(mfem_p) else None
        bench = load_tpv104_file(bench_p, stress_in_pa=False) if os.path.exists(bench_p) else None

        if mfem is None and bench is None:
            print(f"  {station_label}: no data (missing both MFEM and benchmark files)")
            continue
        if mfem is None:
            print(f"  {station_label}: MFEM file missing -> {mfem_p}")
        if bench is None:
            print(f"  {station_label}: benchmark file missing -> {bench_p}")

        save_path = (
            None
            if args.show
            else os.path.join(args.output_dir, f"tpv104_strike_2x2_{station_label}.png")
        )
        plot_station_2x2(mfem, bench, station_label, x2_km, x3_km, save_path=save_path)
        n_plotted += 1

    print(f"\nPlotted {n_plotted} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
