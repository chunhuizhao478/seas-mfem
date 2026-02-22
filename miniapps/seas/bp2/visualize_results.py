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
    python visualize_results.py ../mfem_bp2qd
    python visualize_results.py ../mfem_bp2qd --depths 0 12 --save
    python visualize_results.py ../mfem_bp2qd --benchmark-dir benchmark_data --save
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
ALL_DEPTHS_KM = [0, 2.4, 4.8, 7.2, 9.6, 12, 14.4, 16.8, 19.2, 21.6, 24, 36]


def mfem_filename(prefix, depth_km):
    """Generate MFEM SCEC filename for a given depth."""
    if abs(depth_km - round(depth_km)) < 1e-6:
        return f"{prefix}_z{int(round(depth_km))}km.txt"
    else:
        return f"{prefix}_z{depth_km:.1f}km.txt"


def benchmark_filename(depth_km, bench_dir=None):
    """Generate benchmark reference filename for a given depth.

    Tries two naming conventions:
      1. bp2-qd-erickson-z{depth}km-res.txt  (Erickson reference data)
      2. bp2-qd-z{depth}km-res.txt           (generic benchmark data)

    If bench_dir is provided, returns the first file that exists.
    """
    if abs(depth_km - round(depth_km)) < 1e-6:
        d = f"{int(round(depth_km))}"
    else:
        d = f"{depth_km:.1f}"

    candidates = [
        f"bp2-qd-erickson-z{d}km-res.txt",
        f"bp2-qd-z{d}km-res.txt",
    ]

    if bench_dir is not None:
        for c in candidates:
            if os.path.exists(os.path.join(bench_dir, c)):
                return c

    return candidates[0]


def plot_station(mfem_data, bench_data, depth_km, save_path=None):
    """Plot comparison for one depth station.

    Shows: slip rate (V), slip, shear stress, state variable (theta)
    Uses log scale for slip rate and state variable (matching benchmark convention).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"BP2-QD: z = {depth_km} km", fontsize=14, fontweight="bold")

    # Define the 4 panels matching benchmark format
    panels = [
        ("slip_rate", "Slip Rate V (m/s)", True),  # log y-axis
        ("slip_m", "Slip (m)", False),  # linear y-axis
        ("tau_MPa", "Shear Stress (MPa)", False),  # linear y-axis
        ("theta_s", "State Variable theta (s)", True),  # log y-axis
    ]

    for ax, (key, ylabel, use_log) in zip(axes.flat, panels):
        if bench_data is not None:
            ax.plot(
                bench_data["time_yr"],
                bench_data[key],
                "k-",
                label="Benchmark",
                linewidth=0.8,
            )
        if mfem_data is not None:
            ax.plot(
                mfem_data["time_yr"],
                mfem_data[key],
                "r-",
                label="MFEM",
                linewidth=0.8,
                alpha=0.8,
            )

        ax.set_xlabel("Time (years)")
        ax.set_ylabel(ylabel)
        if use_log:
            ax.set_yscale("log")
        ax.legend(fontsize=9, loc="best")
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

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    fig.suptitle("BP2-QD: Slip Rate at All Depths", fontsize=14, fontweight="bold")

    cmap = plt.cm.viridis
    n = len(all_results)

    for i, res in enumerate(all_results):
        depth_km = res["depth_km"]
        color = cmap(i / max(n - 1, 1))

        if res.get("bench_data") is not None:
            ax.plot(
                res["bench_data"]["time_yr"],
                res["bench_data"]["slip_rate"],
                "-",
                color=color,
                linewidth=0.5,
                alpha=0.5,
            )
        if res.get("mfem_data") is not None:
            ax.plot(
                res["mfem_data"]["time_yr"],
                res["mfem_data"]["slip_rate"],
                "--",
                color=color,
                linewidth=0.8,
                label=f"z={depth_km} km",
            )

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
        description="Visualize MFEM SEAS BP2 output vs benchmark data"
    )
    parser.add_argument(
        "mfem_prefix", help="MFEM output file prefix (e.g., ../mfem_bp2qd)"
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing benchmark reference files",
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

    try:
        import matplotlib
        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        return 1

    depths = args.depths if args.depths else ALL_DEPTHS_KM

    print("=" * 60)
    print("BP2-QD Visualization")
    print("=" * 60)
    print(f"MFEM prefix:   {args.mfem_prefix}")
    if not args.no_benchmark:
        print(f"Benchmark dir: {args.benchmark_dir}")
    print(f"Depths (km):   {depths}")
    print()

    all_results = []

    for depth_km in depths:
        mfem_path = mfem_filename(args.mfem_prefix, depth_km)
        bench_path = os.path.join(
            args.benchmark_dir,
            benchmark_filename(depth_km, bench_dir=args.benchmark_dir))

        mfem_data = None
        bench_data = None

        if os.path.exists(mfem_path):
            mfem_data = load_scec_file(mfem_path)
            if mfem_data is not None:
                print(
                    f"z={depth_km:5.1f} km: MFEM {len(mfem_data['time_s'])} pts, "
                    f"t=[{mfem_data['time_yr'][0]:.1f}, "
                    f"{mfem_data['time_yr'][-1]:.1f}] yr"
                )
        else:
            print(f"z={depth_km:5.1f} km: MFEM file not found ({mfem_path})")

        if not args.no_benchmark and os.path.exists(bench_path):
            bench_data = load_scec_file(bench_path)

        if mfem_data is None and bench_data is None:
            continue

        result = {
            "depth_km": depth_km,
            "mfem_data": mfem_data,
            "bench_data": bench_data,
        }
        all_results.append(result)

        # Plot individual station
        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"bp2_z{depth_km}km.png")
            plot_station(mfem_data, bench_data, depth_km, save_path=fname)
        else:
            plot_station(mfem_data, bench_data, depth_km)

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
