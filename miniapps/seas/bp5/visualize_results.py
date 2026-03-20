#!/usr/bin/env python3
"""
BP5 Benchmark Visualization Script

Plots MFEM SEAS miniapp SCEC-format output alongside Tandem reference data.
BP5 uses 8-column vector format (3D with strike + dip components):
  time(s), slip_strike(m), slip_dip(m), log10(V_strike)(m/s),
  log10(V_dip)(m/s), tau_strike(MPa), tau_dip(MPa), log10(state)(s)

Usage:
    python visualize_results.py <mfem_prefix> [options]

Examples:
    # MFEM vs Tandem benchmark (default)
    python visualize_results.py results_1000m/bp5_full --tandem --save

    # MFEM only, no benchmarks
    python visualize_results.py results_1000m/bp5_full --no-benchmark --save

    # Compare multiple resolutions against Tandem
    python visualize_results.py results_250m/bp5_full --tandem \\
        --compare 500m:results_500m/bp5_full --save --output-dir plots_comparison
"""

import argparse
import os
import sys

import numpy as np

# SCEC BP5 standard 10 on-fault stations: (name, x2_km, x3_km)
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


def load_bp5_file(filepath):
    """Load BP5-format 8-column text file.

    Returns dict with physical quantities (converting log10 columns).
    """
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    data.append([float(x) for x in parts[:8]])
                except ValueError:
                    continue

    if len(data) == 0:
        return None

    arr = np.array(data)

    return {
        "time_s": arr[:, 0],
        "time_yr": arr[:, 0] / SEC_PER_YR,
        "slip_strike": arr[:, 1],
        "slip_dip": arr[:, 2],
        "V_strike": 10.0 ** arr[:, 3],
        "V_dip": 10.0 ** arr[:, 4],
        "log10_V_strike": arr[:, 3],
        "log10_V_dip": arr[:, 4],
        "tau_strike": arr[:, 5],
        "tau_dip": arr[:, 6],
        "theta_s": 10.0 ** arr[:, 7],
        "log10_theta": arr[:, 7],
    }


def coord_str(val_km):
    """Format coordinate for Tandem filenames: -36 -> '-36', 0 -> '0', 10 -> '10'."""
    return str(int(val_km))


def tandem_filename(directory, x2_km, x3_km, order=4):
    """Generate Tandem benchmark filename for a given station and polynomial order."""
    fname = f"bp5qd_tandem_p{order}_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt"
    return os.path.join(directory, fname)


def eqsim_filename(directory, x2_km, x3_km):
    """Generate EQSim benchmark filename for a given station."""
    fname = f"bp5qd_eqsim_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt"
    return os.path.join(directory, fname)


def tribie_filename(directory, x2_km, x3_km):
    """Generate TriBIE benchmark filename for a given station."""
    fname = f"bp5qd_tribie_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt"
    return os.path.join(directory, fname)


def mfem_filename(prefix, station_name):
    """Generate MFEM output filename for a given station."""
    return f"{prefix}_{station_name}.txt"


def plot_station(datasets, station_name, x2_km, x3_km, save_path=None):
    """Plot 8-panel comparison for one BP5 station.

    datasets: list of (label, data_dict, color, linestyle) tuples.
    Shows: slip_strike, slip_dip, V_strike, V_dip, tau_strike, tau_dip, log10(state).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    fig.suptitle(
        f"BP5-QD: {station_name}  (x2={x2_km} km, x3={x3_km} km)",
        fontsize=14,
        fontweight="bold",
    )

    panels = [
        ("slip_strike", "Slip Strike (m)", False),
        ("slip_dip", "Slip Dip (m)", False),
        ("V_strike", "Slip Rate V_strike (m/s)", True),
        ("V_dip", "Slip Rate V_dip (m/s)", True),
        ("tau_strike", "Shear Stress \u03c4_strike (MPa)", False),
        ("tau_dip", "Shear Stress \u03c4_dip (MPa)", False),
        ("log10_theta", "log\u2081\u2080(State) (s)", False),
    ]

    for ax, (key, ylabel, use_log) in zip(axes.flat, panels):
        for label, data, color, ls in datasets:
            if data is not None:
                ax.plot(
                    data["time_yr"],
                    data[key],
                    ls,
                    color=color,
                    label=label,
                    linewidth=0.8,
                    alpha=0.85,
                )
        ax.set_xlabel("Time (years)")
        ax.set_ylabel(ylabel)
        if use_log:
            ax.set_yscale("log")
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)

    # Hide unused subplot
    axes[3, 1].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close()


def plot_closeup(datasets, station_name, x2_km, x3_km, t_max_yr=0.1, save_path=None):
    """Plot 8-panel comparison zoomed to the first t_max_yr years."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    fig.suptitle(
        f"BP5-QD Close-up (0\u2013{t_max_yr:.0f} yr): {station_name}  "
        f"(x2={x2_km} km, x3={x3_km} km)",
        fontsize=14,
        fontweight="bold",
    )

    panels = [
        ("slip_strike", "Slip Strike (m)", False),
        ("slip_dip", "Slip Dip (m)", False),
        ("V_strike", "Slip Rate V_strike (m/s)", True),
        ("V_dip", "Slip Rate V_dip (m/s)", True),
        ("tau_strike", "Shear Stress \u03c4_strike (MPa)", False),
        ("tau_dip", "Shear Stress \u03c4_dip (MPa)", False),
        ("log10_theta", "log\u2081\u2080(State) (s)", False),
    ]

    for ax, (key, ylabel, use_log) in zip(axes.flat, panels):
        for label, data, color, ls in datasets:
            if data is not None:
                mask = data["time_yr"] <= t_max_yr
                ax.plot(
                    data["time_yr"][mask],
                    data[key][mask],
                    ls,
                    color=color,
                    label=label,
                    linewidth=0.8,
                    alpha=0.85,
                )
        ax.set_xlabel("Time (years)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, t_max_yr)
        if use_log:
            ax.set_yscale("log")
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)

    # Hide unused subplot
    axes[3, 1].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close()


def plot_overview(all_results, save_path=None):
    """Plot slip rate (strike component) across all stations."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    fig.suptitle(
        "BP5-QD: Slip Rate (strike) at All Stations",
        fontsize=14,
        fontweight="bold",
    )

    cmap = plt.cm.tab10
    n = len(all_results)

    for i, res in enumerate(all_results):
        name = res["station_name"]
        color = cmap(i / max(n - 1, 1))
        first = True
        for label, data, _, ls in res["datasets"]:
            if data is not None:
                lbl = name if first else None
                ax.plot(
                    data["time_yr"],
                    data["V_strike"],
                    ls,
                    color=color,
                    linewidth=0.6,
                    label=lbl,
                )
                first = False

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Slip Rate V_strike (m/s)")
    ax.set_yscale("log")
    ax.legend(fontsize=7, ncol=2, loc="best")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close()


COLORS = [
    "#000000",  # black
    "#d62728",  # red
    "#1f77b4",  # blue
    "#2ca02c",  # green
    "#9467bd",  # purple
    "#ff7f0e",  # orange
    "#8c564b",  # brown
    "#e377c2",  # pink
]


def main():
    parser = argparse.ArgumentParser(
        description="Visualize MFEM SEAS BP5 output vs Tandem benchmark data"
    )
    parser.add_argument(
        "mfem_prefix", help="MFEM output file prefix (e.g., results_1000m/bp5_full)"
    )
    parser.add_argument(
        "--tandem-p4",
        action="store_true",
        help="Include Tandem p4 benchmark data",
    )
    parser.add_argument(
        "--tandem-p6",
        action="store_true",
        help="Include Tandem p6 benchmark data",
    )
    parser.add_argument(
        "--tandem",
        action="store_true",
        help="Include both Tandem p4 and p6 benchmark data",
    )
    parser.add_argument(
        "--eqsim",
        action="store_true",
        help="Include EQSim benchmark data",
    )
    parser.add_argument(
        "--tribie",
        action="store_true",
        help="Include TriBIE benchmark data",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing benchmark files (Tandem, EQSim)",
    )
    parser.add_argument(
        "--compare",
        action="append",
        metavar="LABEL:PREFIX",
        help="Additional MFEM dataset to overlay (e.g., 500m:results_500m/bp5_full). "
        "Can be repeated for multiple comparisons.",
    )
    parser.add_argument(
        "--stations",
        nargs="+",
        type=int,
        default=None,
        help="Specific station indices (1-10) to plot. Default: all available",
    )
    parser.add_argument(
        "--save", action="store_true", help="Save plots as PNG (default: display)"
    )
    parser.add_argument("--output-dir", default=".", help="Directory for output plots")
    parser.add_argument(
        "--no-benchmark",
        action="store_true",
        help="Plot MFEM data only (no benchmark overlay)",
    )
    parser.add_argument(
        "--flip-dip",
        action="store_true",
        help="Flip sign of slip_dip and tau_dip for MFEM datasets (debug convention mismatch)",
    )
    args = parser.parse_args()

    # --tandem enables both p4 and p6
    if args.tandem:
        args.tandem_p4 = True
        args.tandem_p6 = True

    # Default: Tandem p4 if no benchmark flags specified and not --no-benchmark
    if not args.no_benchmark and not args.tandem_p4 and not args.tandem_p6 and not args.eqsim and not args.tribie:
        args.tandem_p4 = True

    try:
        import matplotlib

        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        return 1

    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), args.benchmark_dir
    )
    if not os.path.isabs(args.benchmark_dir) and not os.path.isdir(data_dir):
        data_dir = args.benchmark_dir

    # Select stations
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build list of data sources
    sources = []
    color_idx = 0

    if args.tandem_p4 and not args.no_benchmark:
        sources.append(("Tandem p4", "tandem_p4", None, "#000000", "--"))  # black dashed

    if args.tandem_p6 and not args.no_benchmark:
        sources.append(("Tandem p6", "tandem_p6", None, "#2ca02c", "--"))  # green dashed

    if args.eqsim and not args.no_benchmark:
        sources.append(("EQSim", "eqsim", None, "#1f77b4", "--"))  # blue dashed

    if args.tribie and not args.no_benchmark:
        sources.append(("TriBIE", "tribie", None, "#9467bd", "--"))  # purple dashed

    # Primary MFEM dataset
    primary_label = f"MFEM {os.path.basename(args.mfem_prefix)}"
    sources.append(
        (primary_label, "mfem", args.mfem_prefix, "#d62728", "-")  # red solid
    )
    color_idx = 3

    # Additional --compare MFEM datasets
    compare_colors = ["#000000", "#ff7f0e"]  # 1st: black, 2nd: orange
    if args.compare:
        for i, spec in enumerate(args.compare):
            if ":" in spec:
                label, prefix = spec.split(":", 1)
            else:
                label = os.path.basename(spec)
                prefix = spec
            color = compare_colors[i] if i < len(compare_colors) else COLORS[(i + 3) % len(COLORS)]
            sources.append(
                (f"MFEM {label}", "mfem", prefix, color, "-")
            )

    print("=" * 60)
    print("BP5-QD Visualization")
    print("=" * 60)
    for label, stype, info, color, ls in sources:
        if stype in ("tandem_p4", "tandem_p6", "eqsim", "tribie"):
            print(f"  {label}: {data_dir}/")
        else:
            print(f"  {label}: {info}")
    print(f"  Stations: {len(stations)}")
    print()

    all_results = []

    for station_name, x2_km, x3_km in stations:
        datasets = []
        for label, stype, info, color, ls in sources:
            data = None
            if stype == "tandem_p4":
                path = tandem_filename(data_dir, x2_km, x3_km, order=4)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "tandem_p6":
                path = tandem_filename(data_dir, x2_km, x3_km, order=6)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "tribie":
                path = tribie_filename(data_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "eqsim":
                path = eqsim_filename(data_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            else:
                path = mfem_filename(info, station_name)
                if os.path.exists(path):
                    data = load_bp5_file(path)
                    if data is not None and args.flip_dip:
                        data["slip_dip"] = -data["slip_dip"]
                        data["tau_dip"] = -data["tau_dip"]
            datasets.append((label, data, color, ls))

        # Skip if no data at this station
        if all(d is None for _, d, _, _ in datasets):
            print(f"  {station_name}: no data found")
            continue

        # Print info
        pts_info = []
        for label, data, _, _ in datasets:
            if data is not None:
                pts_info.append(f"{label}: {len(data['time_s'])} pts")
        print(f"  {station_name} (x2={x2_km}, x3={x3_km}): {', '.join(pts_info)}")

        result = {
            "station_name": station_name,
            "x2_km": x2_km,
            "x3_km": x3_km,
            "datasets": datasets,
        }
        all_results.append(result)

        # Plot full time range
        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"bp5_{station_name}.png")
            plot_station(datasets, station_name, x2_km, x3_km, save_path=fname)
        else:
            plot_station(datasets, station_name, x2_km, x3_km)

        # Plot close-up (first 10 years)
        if args.save:
            fname_close = os.path.join(
                args.output_dir, f"bp5_{station_name}_closeup.png"
            )
            plot_closeup(
                datasets,
                station_name,
                x2_km,
                x3_km,
                t_max_yr=0.1,
                save_path=fname_close,
            )
        else:
            plot_closeup(datasets, station_name, x2_km, x3_km, t_max_yr=0.1)

    # Overview plot
    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "bp5_overview.png")
            plot_overview(all_results, save_path=fname)
        else:
            plot_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
