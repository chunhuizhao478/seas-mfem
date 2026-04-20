#!/usr/bin/env python3
"""
TPV102 Benchmark Visualization Script

Plots MFEM dynamic-rupture station output alongside SCEC reference data
(DRDG3D, PyLith). Nine on-fault stations are compared.

MFEM station format (tpv102_station_flt_{x}_{z}.dat), columns:
  time(s) slip1(m) slip2(m) V1(m/s) V2(m/s) tau1(Pa) tau2(Pa) sigma_n(Pa) log10_theta
  where 1 = dip (tangent1) and 2 = strike (tangent2). TPV102 is pure
  strike-slip so the interesting physics is in column 2 (strike).

SCEC benchmark format (tpv102_{code}_x2_{x2}_x3_{x3}.txt), columns:
  time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(MPa)
         v-slip(m) v-slip-rate(m/s) v-shear-stress(MPa)
         n-stress(MPa) log-theta
  where h = horizontal (strike) and v = vertical (dip).

Usage:
    # MFEM results vs both DRDG3D and PyLith benchmarks
    python visualize_results.py \\
        --mfem /path/to/results_200m_p1_1.5s_400r_v2 --benchmarks --save

    # Only DRDG3D benchmark
    python visualize_results.py \\
        --mfem /path/to/results_dir --drdg3d --save

    # Compare two MFEM runs plus benchmarks
    python visualize_results.py \\
        --mfem "p1 200m:/path/to/results_200m_p1" \\
        --mfem "p2 100m:/path/to/results_100m_p2" \\
        --benchmarks --save

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --benchmarks --save
"""

import argparse
import os
import sys

import numpy as np

# SCEC TPV102 on-fault stations: (mfem_name, bench_x2_km, bench_x3_km)
# The MFEM filename encodes negative along-strike as 'n' (e.g., n12 == -12).
# Benchmark files use "-12". Both use decimal for 7.5.
SCEC_STATIONS = [
    ("flt_0_3", 0, 3),
    ("flt_0_7.5", 0, 7.5),
    ("flt_0_12", 0, 12),
    ("flt_9_7.5", 9, 7.5),
    ("flt_12_3", 12, 3),
    ("flt_12_12", 12, 12),
    ("flt_n9_7.5", -9, 7.5),
    ("flt_n12_3", -12, 3),
    ("flt_n12_12", -12, 12),
]

PA_TO_MPA = 1.0e-6


def _parse_numeric_table(filepath, min_cols=9):
    """Load a whitespace-delimited numeric table, skipping '#' comments.

    Also tolerates SCEC-style header lines that begin with a letter
    (e.g., 't h-slip h-slip-rate ...').
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < min_cols:
                continue
            try:
                rows.append([float(x) for x in parts[:min_cols]])
            except ValueError:
                # Skip header/text lines like 't h-slip ...'
                continue
    if not rows:
        return None
    return np.array(rows)


def load_mfem_file(filepath):
    """Load MFEM TPV102 fault-station file.

    Columns in file:
      0: time (s)
      1: slip1 (m)      == dip slip
      2: slip2 (m)      == strike slip
      3: V1 (m/s)       == dip slip rate
      4: V2 (m/s)       == strike slip rate
      5: tau1 (Pa)      == dip shear stress
      6: tau2 (Pa)      == strike shear stress
      7: sigma_n (Pa)   (MFEM convention: positive = compression)
      8: log10(theta)

    Returns dict in a unified (strike/dip, stresses in MPa) schema.
    """
    arr = _parse_numeric_table(filepath, min_cols=9)
    if arr is None:
        return None
    return {
        "time_s": arr[:, 0],
        "slip_strike": arr[:, 2],
        "slip_dip": arr[:, 1],
        "V_strike": arr[:, 4],
        "V_dip": arr[:, 3],
        "tau_strike": arr[:, 6] * PA_TO_MPA,
        "tau_dip": arr[:, 5] * PA_TO_MPA,
        "sigma_n": np.abs(arr[:, 7]) * PA_TO_MPA,
        "log10_theta": arr[:, 8],
    }


def load_scec_file(filepath):
    """Load SCEC benchmark file (DRDG3D or PyLith).

    SCEC standard columns (text; header lines prefixed with '#'):
      0: time (s)
      1: h-slip (m)            == strike slip (horizontal)
      2: h-slip-rate (m/s)     == strike slip rate
      3: h-shear-stress (MPa)  == strike shear stress
      4: v-slip (m)            == dip slip (vertical)
      5: v-slip-rate (m/s)     == dip slip rate
      6: v-shear-stress (MPa)  == dip shear stress
      7: n-stress (MPa)        (sign varies by code; plotted as magnitude)
      8: log10(theta)

    Returns dict in the same unified schema as load_mfem_file.
    """
    arr = _parse_numeric_table(filepath, min_cols=9)
    if arr is None:
        return None
    return {
        "time_s": arr[:, 0],
        "slip_strike": arr[:, 1],
        "slip_dip": arr[:, 4],
        "V_strike": arr[:, 2],
        "V_dip": arr[:, 5],
        "tau_strike": arr[:, 3],
        "tau_dip": arr[:, 6],
        "sigma_n": np.abs(arr[:, 7]),
        "log10_theta": arr[:, 8],
    }


def mfem_coord_str(val_km):
    """MFEM filename coord: -12 -> 'n12', 0 -> '0', 7.5 -> '7.5'."""
    s = f"{val_km:g}"
    if s.startswith("-"):
        s = "n" + s[1:]
    return s


def bench_coord_str(val_km):
    """Benchmark filename coord: -12 -> '-12', 0 -> '0', 7.5 -> '7.5'."""
    return f"{val_km:g}"


def mfem_filename(results_dir, x2_km, x3_km):
    fname = f"tpv102_station_flt_{mfem_coord_str(x2_km)}_{mfem_coord_str(x3_km)}.dat"
    return os.path.join(results_dir, fname)


def drdg3d_filename(bench_dir, x2_km, x3_km):
    fname = f"tpv102_drdg3d_x2_{bench_coord_str(x2_km)}_x3_{bench_coord_str(x3_km)}.txt"
    return os.path.join(bench_dir, fname)


def pylith_filename(bench_dir, x2_km, x3_km):
    fname = f"tpv102_pylith_x2_{bench_coord_str(x2_km)}_x3_{bench_coord_str(x3_km)}.txt"
    return os.path.join(bench_dir, fname)


PANELS = [
    ("V_strike", "Slip Rate V_strike (m/s)", False),
    ("slip_strike", "Slip Strike (m)", False),
    ("tau_strike", "Shear Stress \u03c4_strike (MPa)", False),
    ("V_dip", "Slip Rate V_dip (m/s)", False),
    ("slip_dip", "Slip Dip (m)", False),
    ("tau_dip", "Shear Stress \u03c4_dip (MPa)", False),
    ("sigma_n", "|Normal Stress| (MPa)", False),
    ("log10_theta", "log\u2081\u2080(State) (s)", False),
]


def plot_station(datasets, station_name, x2_km, x3_km, save_path=None, t_max=None):
    """Plot 8-panel station comparison. Time axis is seconds (dynamic rupture)."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    title = f"TPV102: {station_name}  (x2={x2_km} km, x3={x3_km} km)"
    if t_max is not None:
        title += f"   [0\u2013{t_max:g} s close-up]"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    for ax, (key, ylabel, _use_log) in zip(axes.flat, PANELS):
        for label, data, color, ls in datasets:
            if data is None:
                continue
            t = data["time_s"]
            y = data[key]
            if t_max is not None:
                mask = t <= t_max
                if not np.any(mask):
                    continue
                t = t[mask]
                y = y[mask]
            ax.plot(t, y, ls, color=color, label=label, linewidth=0.9, alpha=0.85)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel)
        if t_max is not None:
            ax.set_xlim(0, t_max)
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
    """Plot V_strike at all stations on one axis."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    fig.suptitle(
        "TPV102: Slip Rate V_strike at All Stations",
        fontsize=14,
        fontweight="bold",
    )

    cmap = plt.cm.tab10
    n = len(all_results)
    for i, res in enumerate(all_results):
        name = res["station_name"]
        color = cmap(i / max(n - 1, 1))
        first = True
        for _label, data, _c, ls in res["datasets"]:
            if data is None:
                continue
            lbl = name if first else None
            ax.plot(
                data["time_s"], data["V_strike"], ls,
                color=color, linewidth=0.6, label=lbl,
            )
            first = False

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Slip Rate V_strike (m/s)")
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
]


def parse_labeled_arg(spec):
    """Parse 'label:value' or just 'value' -> (label, value)."""
    if ":" in spec:
        label, value = spec.split(":", 1)
        return label, value
    return None, spec


def main():
    parser = argparse.ArgumentParser(
        description="Visualize TPV102 output: MFEM vs SCEC reference data"
    )
    # Primary datasets
    parser.add_argument(
        "mfem_dir_positional",
        nargs="?",
        default=None,
        help="(Legacy) MFEM results directory containing tpv102_station_flt_*.dat",
    )
    parser.add_argument(
        "--mfem",
        action="append",
        metavar="[LABEL:]DIR",
        help="MFEM results directory. Use 'label:dir' for a custom label. Repeatable.",
    )

    # Benchmark references
    parser.add_argument(
        "--drdg3d", action="store_true",
        help="Include DRDG3D (Wenqiang Zhang) benchmark data",
    )
    parser.add_argument(
        "--pylith", action="store_true",
        help="Include PyLith (Brad Aagaard) benchmark data",
    )
    parser.add_argument(
        "--benchmarks", action="store_true",
        help="Include both DRDG3D and PyLith benchmark data",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing scec_drdg3d/ and scec_pylith/ subfolders",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip all benchmark references",
    )

    # Output options
    parser.add_argument(
        "--stations", nargs="+", type=int, default=None,
        help="Specific station indices (1-9). Default: all",
    )
    parser.add_argument(
        "--save", action="store_true",
        help="Save plots as PNG (default: display interactively)",
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory for output plots",
    )
    parser.add_argument(
        "--closeup-t", type=float, default=None,
        help="If set, additionally produce a close-up plot truncated to this time [s]",
    )

    args = parser.parse_args()

    if args.benchmarks:
        args.drdg3d = True
        args.pylith = True

    # Default: both benchmarks if none specified and not --no-benchmark
    if (
        not args.no_benchmark
        and not args.drdg3d
        and not args.pylith
    ):
        args.drdg3d = True
        args.pylith = True

    # Build ordered sources following command-line order so colors/legend match
    ordered_sources = []
    if args.mfem_dir_positional:
        ordered_sources.append(("mfem", args.mfem_dir_positional))

    mfem_iter = iter(args.mfem or [])
    for arg in sys.argv[1:]:
        if arg == "--benchmarks":
            ordered_sources.append(("drdg3d", None))
            ordered_sources.append(("pylith", None))
        elif arg == "--drdg3d":
            ordered_sources.append(("drdg3d", None))
        elif arg == "--pylith":
            ordered_sources.append(("pylith", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))

    if not ordered_sources:
        parser.error(
            "No data sources specified. Provide --mfem DIR and/or benchmark flags."
        )

    # If user only passed benchmark flags with no MFEM, still allow it.
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
    drdg3d_dir = os.path.join(data_dir, "scec_drdg3d")
    pylith_dir = os.path.join(data_dir, "scec_pylith")

    # Select stations
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations
            if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build sources list in command-line order; drop benchmarks if --no-benchmark
    # Each tuple: (label, source_type, dir_or_None, color, linestyle)
    benchmark_types = {"drdg3d", "pylith"}
    benchmark_labels = {
        "drdg3d": "DRDG3D (ref)",
        "pylith": "PyLith (ref)",
    }
    sources = []
    ci = 0
    for stype, spec in ordered_sources:
        if args.no_benchmark and stype in benchmark_types:
            continue
        color = COLORS[ci % len(COLORS)]
        ci += 1
        if stype in benchmark_types:
            sources.append((benchmark_labels[stype], stype, None, color, "--"))
        elif stype == "mfem":
            label, directory = parse_labeled_arg(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            sources.append((f"MFEM {label}", "mfem", directory, color, "-"))

    # Print summary
    print("=" * 60)
    print("TPV102 Visualization")
    print("=" * 60)
    for label, stype, info, _color, ls in sources:
        style = "dashed" if ls == "--" else "solid"
        if stype == "drdg3d":
            print(f"  [{style}] {label}: {drdg3d_dir}/")
        elif stype == "pylith":
            print(f"  [{style}] {label}: {pylith_dir}/")
        else:
            print(f"  [{style}] {label}: {info}/tpv102_station_flt_*.dat")
    print(f"  Stations: {len(stations)}")
    print()

    all_results = []
    for station_name, x2_km, x3_km in stations:
        datasets = []
        for label, stype, info, color, ls in sources:
            data = None
            if stype == "drdg3d":
                path = drdg3d_filename(drdg3d_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_scec_file(path)
            elif stype == "pylith":
                path = pylith_filename(pylith_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_scec_file(path)
            elif stype == "mfem":
                path = mfem_filename(info, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_mfem_file(path)
            datasets.append((label, data, color, ls))

        if all(d is None for _, d, _, _ in datasets):
            print(f"  {station_name}: no data found")
            continue

        pts_info = []
        for label, data, _c, _ls in datasets:
            if data is not None:
                t_last = data["time_s"][-1]
                pts_info.append(
                    f"{label}: {len(data['time_s'])} pts ({t_last:.2f} s)"
                )
        print(f"  {station_name} (x2={x2_km}, x3={x3_km}): {', '.join(pts_info)}")

        all_results.append({
            "station_name": station_name,
            "x2_km": x2_km,
            "x3_km": x3_km,
            "datasets": datasets,
        })

        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(args.output_dir, f"tpv102_{station_name}.png")
            plot_station(datasets, station_name, x2_km, x3_km, save_path=fname)
            if args.closeup_t is not None:
                fname_close = os.path.join(
                    args.output_dir, f"tpv102_{station_name}_closeup.png"
                )
                plot_station(
                    datasets, station_name, x2_km, x3_km,
                    save_path=fname_close, t_max=args.closeup_t,
                )
        else:
            plot_station(datasets, station_name, x2_km, x3_km)
            if args.closeup_t is not None:
                plot_station(
                    datasets, station_name, x2_km, x3_km, t_max=args.closeup_t,
                )

    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "tpv102_overview.png")
            plot_overview(all_results, save_path=fname)
        else:
            plot_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
