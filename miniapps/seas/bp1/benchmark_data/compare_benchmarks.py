#!/usr/bin/env python3
"""
Compare BP1-QD benchmark and simulation results.

Plots 4-panel comparison (slip rate, slip, shear stress, state variable)
at each of the BP1 SCEC probe depths.

Supports any combination of:
  - Binhao Wang (SBIEM) benchmark
  - Junle Jiang (SBIM) benchmark
  - Dal Zilio (GARNET) benchmark
  - Harvey (FDCycle) benchmark
  - SCycle benchmark
  - MFEM simulation results at various resolutions

Usage:
    # Default: plot all benchmarks at z=0km
    python3 compare_benchmarks.py

    # Include MFEM results
    python3 compare_benchmarks.py \
        --mfem 100m:../results_100m/bp1_full \
        --mfem 50m:../results_50m/bp1_full

    # Only specific benchmarks
    python3 compare_benchmarks.py --binhaowang --junlejiang

    # All benchmarks at multiple depths (only binhaowang has all depths)
    python3 compare_benchmarks.py --depths 0 7.5 12.5 --save

    # Save plots as PNG
    python3 compare_benchmarks.py --save
"""

import argparse
import glob
import os
import re
import sys
import numpy as np

DEPTHS_KM = [0, 2.5, 5.0, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35]
SECONDS_PER_YEAR = 365.25 * 24 * 3600

# Known benchmark codes and their display labels
BENCHMARK_CODES = {
    "binhaowang": "Binhao Wang (SBIEM)",
    "junlejiang": "Junle Jiang (SBIM)",
    "dalzilio":   "Dal Zilio (GARNET)",
    "harvey":     "Harvey (FDCycle)",
    "scycle":     "SCycle",
    "ozawa":      "So Ozawa (HBI)",
    "tandem":     "Tandem (Uphoff)",
}

# Distinct colors for datasets
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

LINE_STYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":"]


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
    """Find benchmark file, handling naming inconsistencies (z25km vs z25.0km)."""
    if abs(depth_km - round(depth_km)) < 1e-6:
        candidates = [
            f"bp1-qd-{code}-z{int(round(depth_km))}km-res.txt",
            f"bp1-qd-{code}-z{depth_km:.1f}km-res.txt",
        ]
    else:
        candidates = [
            f"bp1-qd-{code}-z{depth_km:.1f}km-res.txt",
            f"bp1-qd-{code}-z{int(round(depth_km))}km-res.txt",
        ]
    for c in candidates:
        path = os.path.join(directory, c)
        if os.path.exists(path):
            return path
    return None


def detect_benchmark_codes(directory, depth_km):
    """Auto-detect all benchmark codes that have data at a given depth."""
    codes = []
    for code in BENCHMARK_CODES:
        if find_benchmark_file(directory, code, depth_km) is not None:
            codes.append(code)
    # Also scan for unknown codes
    pattern = os.path.join(directory, "bp1-qd-*-res.txt")
    for path in glob.glob(pattern):
        fname = os.path.basename(path)
        m = re.match(r"bp1-qd-(.+)-z[\d.]+km-res\.txt", fname)
        if m:
            code = m.group(1)
            if code not in codes and code not in BENCHMARK_CODES:
                codes.append(code)
    return codes


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
    fig.suptitle(f"BP1-QD: z = {depth_km} km", fontsize=14, fontweight="bold")

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
    fig.suptitle("BP1-QD: Slip Rate at All Depths", fontsize=14, fontweight="bold")

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


def main():
    parser = argparse.ArgumentParser(
        description="Compare BP1-QD benchmark and simulation results")
    parser.add_argument("--binhaowang", action="store_true",
                        help="Include Binhao Wang (SBIEM) benchmark")
    parser.add_argument("--junlejiang", action="store_true",
                        help="Include Junle Jiang (SBIM) benchmark")
    parser.add_argument("--dalzilio", action="store_true",
                        help="Include Dal Zilio (GARNET) benchmark")
    parser.add_argument("--harvey", action="store_true",
                        help="Include Harvey (FDCycle) benchmark")
    parser.add_argument("--scycle", action="store_true",
                        help="Include SCycle benchmark")
    parser.add_argument("--ozawa", action="store_true",
                        help="Include So Ozawa (HBI) benchmark")
    parser.add_argument("--tandem", action="store_true",
                        help="Include Tandem (Uphoff) benchmark")
    parser.add_argument("--all-benchmarks", action="store_true",
                        help="Include all available benchmarks")
    parser.add_argument("--mfem", action="append", metavar="LABEL:PREFIX",
                        help="MFEM result to include (e.g., 100m:../results_100m/bp1_full). "
                             "Can be repeated.")
    parser.add_argument("--depths", nargs="+", type=float, default=None,
                        help="Depths to plot (km). Default: 0")
    parser.add_argument("--save", action="store_true",
                        help="Save plots as PNG (default: display)")
    parser.add_argument("--output-dir", default="benchmark_comparison_plots",
                        help="Directory for output plots")
    args = parser.parse_args()

    data_dir = os.path.dirname(os.path.abspath(__file__))

    # Determine which benchmarks to include
    any_benchmark = (args.binhaowang or args.junlejiang or args.dalzilio
                     or args.harvey or args.scycle or args.ozawa
                     or args.tandem or args.all_benchmarks)

    # Default: include all benchmarks at z=0km
    if not any_benchmark:
        args.all_benchmarks = True
        if args.depths is None:
            args.depths = [0]

    depths = args.depths if args.depths else [0]

    try:
        import matplotlib
        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib required. Install with: pip install matplotlib")
        return 1

    # Build list of data sources: (label, type, extra_info, color, linestyle)
    sources = []
    color_idx = 0

    # Order of benchmarks to add
    benchmark_flags = [
        ("binhaowang", args.binhaowang or args.all_benchmarks),
        ("junlejiang", args.junlejiang or args.all_benchmarks),
        ("dalzilio", args.dalzilio or args.all_benchmarks),
        ("harvey", args.harvey or args.all_benchmarks),
        ("scycle", args.scycle or args.all_benchmarks),
        ("ozawa", args.ozawa or args.all_benchmarks),
        ("tandem", args.tandem or args.all_benchmarks),
    ]

    for code, enabled in benchmark_flags:
        if enabled:
            label = BENCHMARK_CODES.get(code, code)
            ls = LINE_STYLES[color_idx % len(LINE_STYLES)]
            sources.append((label, "benchmark", code,
                            COLORS[color_idx % len(COLORS)], ls))
            color_idx += 1

    if args.mfem:
        for spec in args.mfem:
            if ":" in spec:
                label, prefix = spec.split(":", 1)
            else:
                label = os.path.basename(spec)
                prefix = spec
            sources.append((f"MFEM {label}", "mfem", prefix,
                            COLORS[color_idx % len(COLORS)], "-."))
            color_idx += 1

    print("=" * 60)
    print("BP1-QD Comparison")
    print("=" * 60)
    for label, stype, info, color, ls in sources:
        if stype == "benchmark":
            print(f"  {label}")
        else:
            print(f"  {label}: {info}")
    print(f"  Depths: {depths}")
    print()

    # Print error table header (pairwise against first source)
    ref_label = sources[0][0] if sources else ""
    if len(sources) >= 2:
        print(f"  Relative L2 errors vs {ref_label}:")
        print(f"  {'Dataset':>25s} {'Depth(km)':>10s} {'Slip':>12s} "
              f"{'Rate':>12s} {'Stress':>12s} {'State':>12s}")
        print("  " + "-" * 83)

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
            print(f"  z={depth_km:5.1f} km: no data found, skipping")
            continue

        # Print info
        pts_info = []
        for label, data, _, _ in datasets:
            if data is not None:
                pts_info.append(f"{label}: {len(data['time_s'])} pts")
            else:
                pts_info.append(f"{label}: not found")
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
                    print(f"  {label:>25s} {depth_km:10.1f} {err_slip:12.4e} "
                          f"{err_rate:12.4e} {err_stress:12.4e} {err_state:12.4e}")

        result = {"depth_km": depth_km, "datasets": datasets}
        all_results.append(result)

        # Plot
        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"bp1_z{depth_km}km.png")
            plot_comparison(datasets, depth_km, save_path=fname)
        else:
            plot_comparison(datasets, depth_km)

    # Overview plot
    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "bp1_overview.png")
            plot_overview(all_results, [s[0] for s in sources], save_path=fname)
        else:
            plot_overview(all_results, [s[0] for s in sources])

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
