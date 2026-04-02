#!/usr/bin/env python3
"""
BP5 Benchmark Visualization Script

Plots MFEM and/or Tandem SCEC-format output alongside benchmark reference data.
BP5 uses 8-column vector format (3D with strike + dip components):
  time(s), slip_strike(m), slip_dip(m), log10(V_strike)(m/s),
  log10(V_dip)(m/s), tau_strike(MPa), tau_dip(MPa), log10(state)(s)

Usage:
    # MFEM results vs Tandem p4 benchmark
    python visualize_results.py --mfem results_1000m/bp5_full --tandem --save

    # Tandem cluster results vs benchmark
    python visualize_results.py --tandem-results /path/to/bp5qd_tandem_p1 --tandem --save

    # Compare MFEM vs Tandem cluster results vs benchmarks
    python visualize_results.py --mfem results_1000m/bp5_full \\
        --tandem-results /path/to/bp5qd_tandem_p1 --tandem --save

    # Multiple datasets
    python visualize_results.py --mfem "1000m:results_1000m/bp5_full" \\
        --mfem "500m:results_500m/bp5_full" --tandem --save

    # Legacy mode (positional argument)
    python visualize_results.py results_1000m/bp5_full --tandem --save
"""

import argparse
import os
import re
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
    """Format coordinate for filenames: -36 -> '-36', 0 -> '0', 10 -> '10'."""
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


def tandem_results_filename(prefix, x2_km, x3_km):
    """Generate Tandem results filename from prefix.

    Prefix format: /path/to/bp5qd_tandem_p1
    -> /path/to/bp5qd_tandem_p1_x2_{x2}_x3_{x3}.txt
    """
    return f"{prefix}_x2_{coord_str(x2_km)}_x3_{coord_str(x3_km)}.txt"


def plot_station(datasets, station_name, x2_km, x3_km, save_path=None):
    """Plot 8-panel comparison for one BP5 station.

    datasets: list of (label, data_dict, color, linestyle) tuples.
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


def plot_closeup(datasets, station_name, x2_km, x3_km, t_max_yr=1.0, save_path=None):
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
                if np.any(mask):
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
    "#ff7f0e",  # orange
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # cyan
    "#bcbd22",  # olive
    "#7f7f7f",  # gray
    "#aec7e8",  # light blue
]


def parse_labeled_arg(spec):
    """Parse 'label:value' or just 'value' -> (label, value)."""
    if ":" in spec:
        label, value = spec.split(":", 1)
        return label, value
    return None, spec


def main():
    parser = argparse.ArgumentParser(
        description="Visualize BP5 output: MFEM, Tandem cluster results, and benchmarks"
    )
    # Primary datasets
    parser.add_argument(
        "mfem_prefix_positional",
        nargs="?",
        default=None,
        help="(Legacy) MFEM output prefix (e.g., results_1000m/bp5_full)",
    )
    parser.add_argument(
        "--mfem",
        action="append",
        metavar="[LABEL:]PREFIX",
        help="MFEM output prefix. Use 'label:prefix' for custom label. Repeatable.",
    )
    parser.add_argument(
        "--tandem-results",
        action="append",
        metavar="[LABEL:]PREFIX",
        help="Tandem results prefix (e.g., 'p1 1000m:/path/to/bp5qd_tandem_p1'). "
        "Files: {prefix}_x2_{x2}_x3_{x3}.txt. Repeatable.",
    )

    # Benchmark references
    parser.add_argument(
        "--tandem-p4", action="store_true", help="Include Tandem p4 benchmark data"
    )
    parser.add_argument(
        "--tandem-p6", action="store_true", help="Include Tandem p6 benchmark data"
    )
    parser.add_argument(
        "--tandem",
        action="store_true",
        help="Include both Tandem p4 and p6 benchmark data",
    )
    parser.add_argument(
        "--eqsim", action="store_true", help="Include EQSim benchmark data"
    )
    parser.add_argument(
        "--tribie", action="store_true", help="Include TriBIE benchmark data"
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing benchmark files",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true", help="Skip all benchmark references"
    )

    # Output options
    parser.add_argument(
        "--stations",
        nargs="+",
        type=int,
        default=None,
        help="Specific station indices (1-10). Default: all",
    )
    parser.add_argument(
        "--save", action="store_true", help="Save plots as PNG (default: display)"
    )
    parser.add_argument("--output-dir", default=".", help="Directory for output plots")
    parser.add_argument(
        "--flip-dip",
        action="store_true",
        help="Flip sign of slip_dip and tau_dip for MFEM datasets",
    )

    # Legacy compat
    parser.add_argument(
        "--compare",
        action="append",
        metavar="LABEL:PREFIX",
        help="(Legacy) Additional MFEM dataset to overlay",
    )

    args = parser.parse_args()

    # Handle --tandem shorthand
    if args.tandem:
        args.tandem_p4 = True
        args.tandem_p6 = True

    # Default: Tandem p4 if no benchmark flags and not --no-benchmark
    if (
        not args.no_benchmark
        and not args.tandem_p4
        and not args.tandem_p6
        and not args.eqsim
        and not args.tribie
    ):
        args.tandem_p4 = True

    # Build ordered source list following command-line order.
    # Scan sys.argv to determine the order of data source flags.
    ordered_sources = []  # list of (source_type, spec_or_none)
    if args.mfem_prefix_positional:
        ordered_sources.insert(0, ("mfem", args.mfem_prefix_positional))

    mfem_iter = iter(args.mfem or [])
    tandem_iter = iter(args.tandem_results or [])
    compare_iter = iter(args.compare or [])

    for i, arg in enumerate(sys.argv[1:]):
        if arg == "--tandem":
            ordered_sources.append(("tandem_p4", None))
            ordered_sources.append(("tandem_p6", None))
        elif arg == "--tandem-p4":
            ordered_sources.append(("tandem_p4", None))
        elif arg == "--tandem-p6":
            ordered_sources.append(("tandem_p6", None))
        elif arg == "--eqsim":
            ordered_sources.append(("eqsim", None))
        elif arg == "--tribie":
            ordered_sources.append(("tribie", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))
        elif arg == "--tandem-results":
            ordered_sources.append(("tandem_results", next(tandem_iter)))
        elif arg == "--compare":
            ordered_sources.append(("mfem", next(compare_iter)))

    # Check we have at least one data source
    if not ordered_sources:
        parser.error(
            "No data sources specified. Use --mfem, --tandem-results, "
            "or benchmark flags (--tandem, --eqsim, etc.)"
        )

    try:
        import matplotlib

        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib required. Install with: pip install matplotlib")
        return 1

    # Resolve benchmark directory
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

    # Build sources list in command-line order
    # Each dataset gets a unique color from the shared COLORS pool
    # Each entry: (label, source_type, extra_info, color, linestyle)
    sources = []
    ci = 0

    benchmark_types = {"tandem_p4", "tandem_p6", "eqsim", "tribie"}
    benchmark_labels = {
        "tandem_p4": "Tandem p4 (ref)",
        "tandem_p6": "Tandem p6 (ref)",
        "eqsim": "EQSim (ref)",
        "tribie": "TriBIE (ref)",
    }

    for stype, spec in ordered_sources:
        if args.no_benchmark and stype in benchmark_types:
            continue

        color = COLORS[ci % len(COLORS)]
        ci += 1

        if stype in benchmark_types:
            sources.append((benchmark_labels[stype], stype, None, color, "--"))
        elif stype == "tandem_results":
            label, prefix = parse_labeled_arg(spec)
            if label is None:
                label = os.path.basename(prefix)
            sources.append((f"Tandem {label}", "tandem_results", prefix, color, "-"))
        elif stype == "mfem":
            label, prefix = parse_labeled_arg(spec)
            if label is None:
                label = os.path.basename(prefix)
            sources.append((f"MFEM {label}", "mfem", prefix, color, "-"))

    # Print summary
    print("=" * 60)
    print("BP5-QD Visualization")
    print("=" * 60)
    for label, stype, info, color, ls in sources:
        style = "dashed" if ls == "--" else "solid"
        if stype in ("tandem_p4", "tandem_p6", "eqsim", "tribie"):
            print(f"  [{style}] {label}: {data_dir}/")
        elif stype == "tandem_results":
            print(f"  [{style}] {label}: {info}_x2_*_x3_*.txt")
        else:
            print(f"  [{style}] {label}: {info}_*.txt")
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
            elif stype == "eqsim":
                path = eqsim_filename(data_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "tribie":
                path = tribie_filename(data_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "tandem_results":
                path = tandem_results_filename(info, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_bp5_file(path)
            elif stype == "mfem":
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
                t_max = data["time_yr"][-1]
                pts_info.append(f"{label}: {len(data['time_s'])} pts ({t_max:.1f} yr)")
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

        # Plot close-up (first year)
        if args.save:
            fname_close = os.path.join(
                args.output_dir, f"bp5_{station_name}_closeup.png"
            )
            plot_closeup(
                datasets,
                station_name,
                x2_km,
                x3_km,
                t_max_yr=1e-5,
                save_path=fname_close,
            )
        else:
            plot_closeup(datasets, station_name, x2_km, x3_km, t_max_yr=1.0)

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
