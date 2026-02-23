#!/usr/bin/env python3
"""
BP2 Benchmark Visualization Script

Plots MFEM SEAS miniapp SCEC-format output alongside reference benchmark data.
Displays: slip rate (V, m/s), slip (m), shear stress (MPa), state variable (theta, s).

Both MFEM and benchmark files use the SCEC 5-column format:
  time(s)  slip(m)  log10(V)(m/s)  tau(MPa)  log10(theta)(s)

Usage:
    python visualize_results.py <mfem_prefix> [options]

Examples:
    # MFEM vs Erickson benchmark (default if no benchmark flags given)
    python visualize_results.py ../mfem_bp2qd --erickson --save

    # MFEM vs both benchmarks
    python visualize_results.py ../mfem_bp2qd --erickson --unicycle --save

    # MFEM only, no benchmarks
    python visualize_results.py ../mfem_bp2qd --no-benchmark --save

    # Compare multiple resolutions against Erickson
    python visualize_results.py results_50m/bp2_full --erickson \
        --compare 200m:results_200m/bp2_full --compare 100m:results_100m/bp2_full \
        --save --output-dir plots_comparison

    # Specific depths only
    python visualize_results.py ../mfem_bp2qd --erickson --depths 0 7.2 12 --save
"""

import argparse
import os
import sys

import numpy as np


def load_scec_file(filepath):
    """Load SCEC-format 5-column text file.

    Returns dict with physical quantities (converting log10 columns).
    """
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    data.append([float(x) for x in parts[:5]])
                except ValueError:
                    continue

    if len(data) == 0:
        return None

    arr = np.array(data)
    sec_per_yr = 3.15576e7

    return {
        "time_s": arr[:, 0],
        "time_yr": arr[:, 0] / sec_per_yr,
        "slip_m": arr[:, 1],
        "slip_rate": 10.0 ** arr[:, 2],  # Convert log10(V) -> V
        "log10_V": arr[:, 2],
        "tau_MPa": arr[:, 3],
        "theta_s": 10.0 ** arr[:, 4],  # Convert log10(theta) -> theta
        "log10_theta": arr[:, 4],
    }


# SCEC benchmark depth stations (km)
ALL_DEPTHS_KM = [0, 2.4, 4.8, 7.2, 9.6, 12, 14.4, 16.8, 19.2, 24, 28.8, 36]


def mfem_filename(prefix, depth_km):
    """Generate MFEM SCEC filename for a given depth."""
    if abs(depth_km - round(depth_km)) < 1e-6:
        return f"{prefix}_z{int(round(depth_km))}km.txt"
    else:
        return f"{prefix}_z{depth_km:.1f}km.txt"


def find_benchmark_file(directory, code, depth_km):
    """Find benchmark file, handling naming inconsistencies (z24km vs z24.0km).

    code: 'erickson' or 'unicycle'
    """
    if abs(depth_km - round(depth_km)) < 1e-6:
        candidates = [
            f"bp2-qd-{code}-z{int(round(depth_km))}km-res.txt",
            f"bp2-qd-{code}-z{depth_km:.1f}km-res.txt",
        ]
    else:
        candidates = [
            f"bp2-qd-{code}-z{depth_km:.1f}km-res.txt",
            f"bp2-qd-{code}-z{int(round(depth_km))}km-res.txt",
        ]
    for c in candidates:
        path = os.path.join(directory, c)
        if os.path.exists(path):
            return path
    return None


def plot_station(datasets, depth_km, save_path=None):
    """Plot 4-panel comparison for one depth station.

    datasets: list of (label, data_dict, color, linestyle) tuples.
    Shows: slip rate (V), slip, shear stress, state variable (theta).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"BP2-QD: z = {depth_km} km", fontsize=14, fontweight="bold")

    panels = [
        ("slip_rate", "Slip Rate V (m/s)", True),
        ("slip_m", "Slip (m)", False),
        ("tau_MPa", "Shear Stress (MPa)", False),
        ("theta_s", "State Variable \u03b8 (s)", True),
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

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close()


def plot_all_depths_overview(all_results, save_path=None):
    """Plot slip rate time series for all available depths in one figure."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    fig.suptitle("BP2-QD: Slip Rate at All Depths", fontsize=14, fontweight="bold")

    cmap = plt.cm.viridis
    n = len(all_results)

    for i, res in enumerate(all_results):
        depth_km = res["depth_km"]
        color = cmap(i / max(n - 1, 1))
        first = True
        for label, data, _, ls in res["datasets"]:
            if data is not None:
                lbl = f"z={depth_km} km" if first else None
                ax.plot(
                    data["time_yr"],
                    data["slip_rate"],
                    ls,
                    color=color,
                    linewidth=0.6,
                    label=lbl,
                )
                first = False

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Slip Rate V (m/s)")
    ax.set_yscale("log")
    ax.legend(fontsize=7, ncol=3, loc="best")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    else:
        plt.show()
    plt.close()


COLORS = [
    "#000000",  # black  (first benchmark)
    "#d62728",  # red    (second benchmark)
    "#1f77b4",  # blue
    "#2ca02c",  # green
    "#9467bd",  # purple
    "#ff7f0e",  # orange
    "#8c564b",  # brown
    "#e377c2",  # pink
]

LINE_STYLES = ["-"]


def main():
    parser = argparse.ArgumentParser(
        description="Visualize MFEM SEAS BP2 output vs benchmark data"
    )
    parser.add_argument(
        "mfem_prefix", help="MFEM output file prefix (e.g., results_50m/bp2_full)"
    )
    parser.add_argument(
        "--erickson",
        action="store_true",
        help="Include Erickson (FD) benchmark",
    )
    parser.add_argument(
        "--unicycle",
        action="store_true",
        help="Include Unicycle (BEM) benchmark",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data_200m",
        help="Directory containing benchmark reference files",
    )
    parser.add_argument(
        "--compare",
        action="append",
        metavar="LABEL:PREFIX",
        help="Additional MFEM dataset to overlay (e.g., 200m:../bp2_full). "
        "Can be repeated for multiple comparisons.",
    )
    parser.add_argument(
        "--depths",
        nargs="+",
        type=float,
        default=None,
        help="Specific depths to plot (km). Default: all available",
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
    args = parser.parse_args()

    # Default: Erickson if no benchmark flags specified and not --no-benchmark
    if not args.no_benchmark and not args.erickson and not args.unicycle:
        args.erickson = True

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

    depths = args.depths if args.depths else ALL_DEPTHS_KM

    # Build list of data sources: (label, type, extra_info, color, linestyle)
    sources = []
    color_idx = 0

    if args.erickson and not args.no_benchmark:
        sources.append(
            (
                "Erickson (FD)",
                "benchmark",
                "erickson",
                COLORS[color_idx],
                LINE_STYLES[0],
            )
        )
        color_idx += 1

    if args.unicycle and not args.no_benchmark:
        sources.append(
            (
                "Unicycle (BEM)",
                "benchmark",
                "unicycle",
                COLORS[color_idx],
                LINE_STYLES[1 if color_idx > 0 else 0],
            )
        )
        color_idx += 1

    # Primary MFEM dataset
    sources.append(
        (
            "MFEM",
            "mfem",
            args.mfem_prefix,
            COLORS[color_idx % len(COLORS)],
            LINE_STYLES[min(color_idx, len(LINE_STYLES) - 1)],
        )
    )
    color_idx += 1

    # Additional --compare MFEM datasets
    if args.compare:
        for spec in args.compare:
            if ":" in spec:
                label, prefix = spec.split(":", 1)
            else:
                label = os.path.basename(spec)
                prefix = spec
            ls = LINE_STYLES[min(color_idx, len(LINE_STYLES) - 1)]
            sources.append(
                (f"MFEM {label}", "mfem", prefix, COLORS[color_idx % len(COLORS)], ls)
            )
            color_idx += 1

    print("=" * 60)
    print("BP2-QD Visualization")
    print("=" * 60)
    for label, stype, info, color, ls in sources:
        if stype == "benchmark":
            print(f"  {label}: {data_dir}/{info}")
        else:
            print(f"  {label}: {info}")
    print(f"  Depths: {depths}")
    print()

    all_results = []

    for depth_km in depths:
        # Load data for each source
        datasets = []  # (label, data, color, linestyle)
        for label, stype, info, color, ls in sources:
            data = None
            if stype == "benchmark":
                path = find_benchmark_file(data_dir, info, depth_km)
                if path:
                    data = load_scec_file(path)
            else:
                path = mfem_filename(info, depth_km)
                if os.path.exists(path):
                    data = load_scec_file(path)
            datasets.append((label, data, color, ls))

        # Skip if no data at this depth
        if all(d is None for _, d, _, _ in datasets):
            continue

        # Print info
        pts_info = []
        for label, data, _, _ in datasets:
            if data is not None:
                pts_info.append(f"{label}: {len(data['time_s'])} pts")
        print(f"  z={depth_km:5.1f} km: {', '.join(pts_info)}")

        result = {"depth_km": depth_km, "datasets": datasets}
        all_results.append(result)

        # Plot
        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"bp2_z{depth_km}km.png")
            plot_station(datasets, depth_km, save_path=fname)
        else:
            plot_station(datasets, depth_km)

    # Overview plot
    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "bp2_overview.png")
            plot_all_depths_overview(all_results, save_path=fname)
        else:
            plot_all_depths_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
