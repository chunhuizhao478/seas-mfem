#!/usr/bin/env python3
"""
BP1 Mesh Convergence Check

Compares MFEM SEAS BP1 results at 100m, 50m, and 25m resolutions
alongside Binhao Wang benchmark to verify mesh convergence.

Plots 4-panel comparison (slip rate, slip, shear stress, state variable)
at each probe depth, with all three resolutions overlaid.

Usage:
    # Default: plot all depths, display interactively
    python mesh_convergence.py

    # Save all plots as PNG
    python mesh_convergence.py --save

    # Specific depths only
    python mesh_convergence.py --depths 0 7.5 12.5 --save

    # Without benchmark
    python mesh_convergence.py --no-benchmark

    # Custom output directory
    python mesh_convergence.py --save --output-dir convergence_plots
"""

import argparse
import os
import sys

import numpy as np

SECONDS_PER_YEAR = 365.25 * 24 * 3600
ALL_DEPTHS_KM = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35]

# Data paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.join(SCRIPT_DIR, "benchmark_data")

# RESULTS_200M = "/Users/chunhuizhao/Downloads/seas-mfem/results_200m"
# RESULTS_100M = "/Users/chunhuizhao/Downloads/seas-mfem/results_100m_ss_bdrload"
RESULTS_75M = "/Users/chunhuizhao/Downloads/seas-mfem/results_75m_ss_bdrload"
RESULTS_50M = "/Users/chunhuizhao/Downloads/seas-mfem/results_50m_ss_bdrload"
RESULTS_25M = "/Users/chunhuizhao/Downloads/seas-mfem/results_25m_ss_bdrload"
RESULTS_12d5M = "/Users/chunhuizhao/Downloads/seas-mfem/results_12.5m_ss_bdrload"

RESULTS_50Mip = "/Users/chunhuizhao/Downloads/seas-mfem/results_50m_ip"
RESULTS_25Mip = "/Users/chunhuizhao/Downloads/seas-mfem/results_25m_ip"
RESULTS_12d5Mip = "/Users/chunhuizhao/Downloads/seas-mfem/results_12.5m_ip"


def load_scec_file(filepath):
    """Load SCEC-format 5-column text file."""
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


def depth_str(depth_km):
    """Format depth for filenames."""
    if abs(depth_km - round(depth_km)) < 1e-6:
        return str(int(round(depth_km)))
    else:
        return f"{depth_km:.1f}"


def mfem_filepath(results_dir, depth_km):
    """Build MFEM output filepath for a given depth."""
    return os.path.join(results_dir, f"bp1_bdrload_z{depth_str(depth_km)}km.txt")


BENCHMARKS = {
    "tandem": ("Tandem", "bp1-qd-tandem"),
    "binhaowang": ("Binhao Wang", "bp1-qd-binhaowang"),
    "ozawa": ("Ozawa", "bp1-qd-ozawa"),
    "junlejiang": ("Junle Jiang", "bp1-qd-junlejiang"),
}


def benchmark_filepath(depth_km, prefix):
    """Build benchmark filepath for a given depth and file prefix."""
    return os.path.join(BENCHMARK_DIR, f"{prefix}-z{depth_str(depth_km)}km-res.txt")


def relative_l2_error(a, b, time_a, time_b):
    """Interpolate b onto a's time grid and compute relative L2 error."""
    t_max = min(time_a[-1], time_b[-1])
    mask = time_a <= t_max
    t = time_a[mask]
    a_vals = a[mask]
    b_interp = np.interp(t, time_b, b)
    num = np.sqrt(np.sum((a_vals - b_interp) ** 2))
    den = np.sqrt(np.sum(b_interp**2))
    if den < 1e-30:
        return 0.0 if num < 1e-30 else 1e30
    return num / den


def plot_convergence(datasets, depth_km, save_path=None):
    """4-panel comparison for one depth with all resolutions."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"BP1-QD Mesh Convergence: z = {depth_km} km", fontsize=14, fontweight="bold"
    )

    panels = [
        ("slip_rate", "Slip Rate V (m/s)", True),
        ("slip_m", "Slip (m)", False),
        ("tau_MPa", "Shear Stress (MPa)", False),
        ("theta_s", "State Variable \u03b8 (s)", True),
    ]

    for ax, (key, ylabel, use_log) in zip(axes.flat, panels):
        for label, data, color, ls, lw in datasets:
            if data is not None:
                ax.plot(
                    data["time_yr"],
                    data[key],
                    ls,
                    color=color,
                    label=label,
                    linewidth=lw,
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


def plot_overview(all_results, save_path=None):
    """Slip rate overview for all depths, colored by depth."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    fig.suptitle(
        "BP1-QD Mesh Convergence: Slip Rate at All Depths",
        fontsize=14,
        fontweight="bold",
    )

    cmap = plt.cm.viridis
    n = len(all_results)

    for i, res in enumerate(all_results):
        depth_km = res["depth_km"]
        color = cmap(i / max(n - 1, 1))
        first = True
        for label, data, _, ls, lw in res["datasets"]:
            if data is not None:
                lbl = f"z={depth_km} km" if first else None
                ax.plot(
                    data["time_yr"],
                    data["slip_rate"],
                    ls,
                    color=color,
                    linewidth=0.5,
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


def main():
    parser = argparse.ArgumentParser(
        description="BP1 mesh convergence check: 100m vs 50m vs 25m"
    )
    parser.add_argument(
        "--depths",
        nargs="+",
        type=float,
        default=None,
        help="Specific depths to plot (km). Default: all 15",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true", help="Exclude all benchmark data"
    )
    parser.add_argument(
        "--benchmark",
        nargs="+",
        default=["tandem"],
        choices=list(BENCHMARKS.keys()),
        help="Benchmark datasets to include (default: tandem)",
    )
    parser.add_argument(
        "--output-dir", default="plots_convergence", help="Directory for output plots"
    )
    args = parser.parse_args()

    depths = args.depths if args.depths else ALL_DEPTHS_KM

    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib required. Install with: pip install matplotlib")
        return 1

    # Define MFEM sources: (label, results_dir, color, linestyle, linewidth)
    resolutions = [
        # ("MFEM 200m", RESULTS_200M, "k", "-", 0.8),
        # ("MFEM 75m", RESULTS_75M, "#d62728", "-", 0.8),
        # ("MFEM 50m", RESULTS_50M, "#1f77b4", "-", 0.8),
        # ("MFEM 25m", RESULTS_25M, "#2ca02c", "-", 0.8),
        # ("MFEM 12.5m", RESULTS_12d5M, "k", "-", 0.8),
        ("MFEM 50m", RESULTS_50Mip, "k", "-", 0.8),
        ("MFEM 25m", RESULTS_25Mip, "#2ca02c", "-", 0.8),
        ("MFEM 12.5m", RESULTS_12d5Mip, "#1f77b4", "-", 0.8),
    ]

    # Resolve benchmark list
    bench_list = []
    if not args.no_benchmark:
        for bkey in args.benchmark:
            label, prefix = BENCHMARKS[bkey]
            bench_list.append((label, prefix))

    print("=" * 60)
    print("BP1-QD Mesh Convergence")
    print("=" * 60)
    for label, _ in bench_list:
        print(f"  {label}: {BENCHMARK_DIR}")
    for label, rdir, _, _, _ in resolutions:
        print(f"  {label}: {rdir}")
    print(f"  Depths: {depths}")
    print()

    # Error table
    print(
        f"  {'':>12s} {'Depth':>8s} {'100m-50m':>12s} {'50m-25m':>12s} {'100m-25m':>12s}"
    )
    print(f"  {'':>12s} {'(km)':>8s} {'rel L2':>12s} {'rel L2':>12s} {'rel L2':>12s}")
    print("  " + "-" * 56)

    all_results = []

    for depth_km in depths:
        datasets = []  # (label, data, color, linestyle, linewidth)

        # Load benchmark data first
        bench_colors = ["#d62728", "#e377c2", "#8c564b", "#9467bd"]
        for i, (blabel, bprefix) in enumerate(bench_list):
            bpath = benchmark_filepath(depth_km, bprefix)
            bdata = None
            if os.path.exists(bpath):
                bdata = load_scec_file(bpath)
            datasets.append(
                (blabel, bdata, bench_colors[i % len(bench_colors)], "-", 1.2)
            )

        # Load each MFEM resolution
        res_data = {}
        for label, rdir, color, ls, lw in resolutions:
            fpath = mfem_filepath(rdir, depth_km)
            data = None
            if os.path.exists(fpath):
                data = load_scec_file(fpath)
            datasets.append((label, data, color, ls, lw))
            res_data[label] = data

        # Skip if no data
        if all(d is None for _, d, _, _, _ in datasets):
            continue

        # Print point counts
        pts_info = []
        for label, data, _, _, _ in datasets:
            if data is not None:
                pts_info.append(f"{label}: {len(data['time_s'])} pts")
            else:
                pts_info.append(f"{label}: not found")
        print(f"  z={depth_km:5.1f} km: {', '.join(pts_info)}")

        # Compute pairwise L2 errors for slip rate
        d100 = res_data.get("MFEM 100m")
        d50 = res_data.get("MFEM 50m")
        d25 = res_data.get("MFEM 25m")

        err_100_50 = err_50_25 = err_100_25 = "    ---     "
        if d100 is not None and d50 is not None:
            err_100_50 = f"{relative_l2_error(d100['log10_V'], d50['log10_V'], d100['time_s'], d50['time_s']):12.4e}"
        if d50 is not None and d25 is not None:
            err_50_25 = f"{relative_l2_error(d50['log10_V'], d25['log10_V'], d50['time_s'], d25['time_s']):12.4e}"
        if d100 is not None and d25 is not None:
            err_100_25 = f"{relative_l2_error(d100['log10_V'], d25['log10_V'], d100['time_s'], d25['time_s']):12.4e}"

        print(
            f"  {'log10(V)':>12s} {depth_km:8.1f} {err_100_50} {err_50_25} {err_100_25}"
        )

        result = {"depth_km": depth_km, "datasets": datasets}
        all_results.append(result)

        # Save plot
        os.makedirs(args.output_dir, exist_ok=True)
        fname = os.path.join(
            args.output_dir, f"bp1_convergence_z{depth_str(depth_km)}km.png"
        )
        plot_convergence(datasets, depth_km, save_path=fname)

    # Overview plot
    if len(all_results) > 1:
        fname = os.path.join(args.output_dir, "bp1_convergence_overview.png")
        plot_overview(all_results, save_path=fname)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
