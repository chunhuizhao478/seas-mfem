#!/usr/bin/env python3
"""
Compare BP2-QD benchmark and simulation results.

Plots 4-panel comparison (slip rate, slip, shear stress, state variable)
at each of the 12 SCEC probe depths.

Supports any combination of:
  - Erickson (FD) benchmark
  - Unicycle (BEM) benchmark
  - MFEM simulation results at various resolutions

Usage:
    # Compare two benchmarks (default)
    python3 compare_benchmarks.py --erickson --unicycle

    # Compare MFEM 50m results against Erickson
    python3 compare_benchmarks.py --erickson --mfem 50m:../results_50m/bp2_full

    # Compare all resolutions against Erickson
    python3 compare_benchmarks.py --erickson \
        --mfem 200m:../../bp2_full \
        --mfem 50m:../results_50m/bp2_full \
        --mfem 25m:../results_25m/bp2_full

    # Just MFEM results, no benchmarks
    python3 compare_benchmarks.py \
        --mfem 50m:../results_50m/bp2_full \
        --mfem 25m:../results_25m/bp2_full

    # Specific depths only
    python3 compare_benchmarks.py --erickson --unicycle --depths 0 7.2 12 --save
"""

import argparse
import os
import sys
import numpy as np

DEPTHS_KM = [0, 2.4, 4.8, 7.2, 9.6, 12, 14.4, 16.8, 19.2, 24, 28.8, 36]
SECONDS_PER_YEAR = 365.25 * 24 * 3600

# Distinct colors for datasets
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

LINE_STYLES = ["-", "--", "-.", ":"]


def load_scec_file(filepath):
    """Load SCEC-format 5-column file: time, slip, log10(V), tau, log10(theta)."""
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("t ") or line.startswith("t\t"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    data.append([float(x) for x in parts[:5]])
                except ValueError:
                    continue
    if not data:
        return None
    arr = np.array(data)
    return {
        "time_yr": arr[:, 0] / SECONDS_PER_YEAR,
        "time_s": arr[:, 0],
        "slip_m": arr[:, 1],
        "slip_rate": 10.0 ** arr[:, 2],
        "log10_V": arr[:, 2],
        "tau_MPa": arr[:, 3],
        "theta_s": 10.0 ** arr[:, 4],
        "log10_theta": arr[:, 4],
    }


def find_benchmark_file(directory, code, depth_km):
    """Find benchmark file, handling naming inconsistencies (z24km vs z24.0km)."""
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


def mfem_filename(prefix, depth_km):
    """Generate MFEM SCEC filename for a given depth."""
    if abs(depth_km - round(depth_km)) < 1e-6:
        return f"{prefix}_z{int(round(depth_km))}km.txt"
    else:
        return f"{prefix}_z{depth_km:.1f}km.txt"


def plot_comparison(datasets, depth_km, save_path=None):
    """4-panel comparison for one depth.

    datasets: list of (label, data_dict, color, linestyle)
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
                ax.plot(data["time_yr"], data[key], ls,
                        color=color, label=label, linewidth=0.8, alpha=0.85)
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


def plot_overview(all_results, dataset_names, save_path=None):
    """Slip rate time series for all depths, all datasets overlaid."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    fig.suptitle("BP2-QD: Slip Rate at All Depths", fontsize=14, fontweight="bold")

    cmap = plt.cm.viridis
    n_depths = len(all_results)

    for i, res in enumerate(all_results):
        color = cmap(i / max(n_depths - 1, 1))
        depth = res["depth_km"]
        first = True
        for label, data, _, ls in res["datasets"]:
            if data is not None:
                lbl = f"z={depth} km" if first else None
                ax.plot(data["time_yr"], data["slip_rate"], ls,
                        color=color, linewidth=0.6, label=lbl)
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


def relative_l2_error(a, b, time_a, time_b):
    """Interpolate b onto a's time grid and compute relative L2 error."""
    b_interp = np.interp(time_a, time_b, b)
    num = np.sqrt(np.sum((a - b_interp) ** 2))
    den = np.sqrt(np.sum(b_interp ** 2))
    if den < 1e-30:
        return 0.0 if num < 1e-30 else 1e30
    return num / den


def nucleation_time(data, threshold_log10V=-3.0):
    """First time log10(V) exceeds threshold (in years)."""
    if data is None:
        return -1.0
    mask = data["log10_V"] > threshold_log10V
    if np.any(mask):
        return data["time_yr"][np.argmax(mask)]
    return -1.0


def main():
    parser = argparse.ArgumentParser(
        description="Compare BP2-QD benchmark and simulation results")
    parser.add_argument("--erickson", action="store_true",
                        help="Include Erickson (FD) benchmark")
    parser.add_argument("--unicycle", action="store_true",
                        help="Include Unicycle (BEM) benchmark")
    parser.add_argument("--mfem", action="append", metavar="LABEL:PREFIX",
                        help="MFEM result to include (e.g., 50m:../results_50m/bp2_full). "
                             "Can be repeated.")
    parser.add_argument("--depths", nargs="+", type=float, default=None,
                        help="Depths to plot (km). Default: all 12")
    parser.add_argument("--save", action="store_true",
                        help="Save plots as PNG (default: display)")
    parser.add_argument("--output-dir", default="benchmark_comparison_plots",
                        help="Directory for output plots")
    args = parser.parse_args()

    # Default: both benchmarks if nothing specified
    if not args.erickson and not args.unicycle and not args.mfem:
        args.erickson = True
        args.unicycle = True

    data_dir = os.path.dirname(os.path.abspath(__file__))
    depths = args.depths if args.depths else DEPTHS_KM

    try:
        import matplotlib
        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib required. Install with: pip install matplotlib")
        return 1

    # Build list of data sources: (label, type, extra_info)
    sources = []
    color_idx = 0

    if args.erickson:
        sources.append(("Erickson (FD)", "benchmark", "erickson",
                         COLORS[color_idx], LINE_STYLES[0]))
        color_idx += 1

    if args.unicycle:
        sources.append(("Unicycle (BEM)", "benchmark", "unicycle",
                         COLORS[color_idx], LINE_STYLES[1 if color_idx > 0 else 0]))
        color_idx += 1

    if args.mfem:
        for spec in args.mfem:
            if ":" in spec:
                label, prefix = spec.split(":", 1)
            else:
                label = os.path.basename(spec)
                prefix = spec
            ls = LINE_STYLES[min(color_idx, len(LINE_STYLES) - 1)]
            sources.append((f"MFEM {label}", "mfem", prefix,
                            COLORS[color_idx % len(COLORS)], ls))
            color_idx += 1

    print("=" * 60)
    print("BP2-QD Comparison")
    print("=" * 60)
    for label, stype, info, color, ls in sources:
        if stype == "benchmark":
            print(f"  {label}: {data_dir}/{info}")
        else:
            print(f"  {label}: {info}")
    print(f"  Depths: {depths}")
    print()

    # Print error table header (pairwise against first source)
    ref_label = sources[0][0] if sources else ""
    if len(sources) >= 2:
        print(f"  Relative L2 errors vs {ref_label}:")
        print(f"  {'Dataset':>20s} {'Depth(km)':>10s} {'Slip':>12s} "
              f"{'Rate':>12s} {'Stress':>12s} {'State':>12s}")
        print("  " + "-" * 78)

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

        # Compute L2 errors against first source
        ref_data = datasets[0][1] if datasets else None
        if ref_data is not None and len(sources) >= 2:
            for label, data, _, _ in datasets[1:]:
                if data is not None:
                    err_slip = relative_l2_error(
                        ref_data["slip_m"], data["slip_m"],
                        ref_data["time_s"], data["time_s"])
                    err_rate = relative_l2_error(
                        ref_data["log10_V"], data["log10_V"],
                        ref_data["time_s"], data["time_s"])
                    err_stress = relative_l2_error(
                        ref_data["tau_MPa"], data["tau_MPa"],
                        ref_data["time_s"], data["time_s"])
                    err_state = relative_l2_error(
                        ref_data["log10_theta"], data["log10_theta"],
                        ref_data["time_s"], data["time_s"])
                    print(f"  {label:>20s} {depth_km:10.1f} {err_slip:12.4e} "
                          f"{err_rate:12.4e} {err_stress:12.4e} {err_state:12.4e}")

        result = {"depth_km": depth_km, "datasets": datasets}
        all_results.append(result)

        # Plot
        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"bp2_z{depth_km}km.png")
            plot_comparison(datasets, depth_km, save_path=fname)
        else:
            plot_comparison(datasets, depth_km)

    # Overview plot
    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "bp2_overview.png")
            plot_overview(all_results, [s[0] for s in sources], save_path=fname)
        else:
            plot_overview(all_results, [s[0] for s in sources])

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
