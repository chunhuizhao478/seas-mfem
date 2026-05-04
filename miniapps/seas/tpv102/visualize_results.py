#!/usr/bin/env python3
"""
TPV102 Benchmark Visualization Script

Plots MFEM dynamic-rupture station output alongside SCEC reference data
(DRDG3D, PyLith) for the SCEC TPV101/102 ageing-law benchmark.  Nine
on-fault stations are compared.

MFEM station format (set in dynamic/tpv102_setup.hpp::TPV102StationWriter):
    <output_prefix>_station_flt_<X>_<Y>.dat
columns:
    time(s) slip1(m) slip2(m) V1(m/s) V2(m/s) tau1(Pa) tau2(Pa) sigma_n(Pa) log10_theta
where 1 = dip (tangent1) and 2 = strike (tangent2).  TPV102 is pure
strike-slip so the interesting physics is in column 2 (strike); column 1
remains ~0.  Filename uses 'n' as the minus-sign prefix
(e.g. flt_n12_3.dat for x2 = -12 km, x3 = 3 km).

DRDG3D reference format (benchmark_data/scec_drdg3d/tpv102_drdg3d_x2_<X>_x3_<Y>.txt):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(MPa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(MPa)
            n-stress(MPa) log10_theta
SCEC standard column order; stresses already in MPa.

PyLith reference format (benchmark_data/scec_pylith/tpv102_pylith_x2_<X>_x3_<Y>.txt):
    Identical 9-column SCEC layout to DRDG3D.

Usage:
    # Single MFEM run vs both references — easiest form, prefix is
    # auto-detected from station files in the directory:
    python visualize_results.py \\
        --mfem /path/to/results_dir --benchmarks --save

    # Compare TWO MFEM runs (each prefix auto-detected; directory basename
    # becomes the legend label):
    python visualize_results.py \\
        --mfem /path/to/run_a \\
        --mfem /path/to/run_b \\
        --benchmarks --save

    # Same, but with custom legend labels:
    python visualize_results.py \\
        --mfem "MF=adjacent:/path/to/run_a" \\
        --mfem "MF=none:/path/to/run_b" \\
        --benchmarks --save

    # Force a specific prefix when auto-detection can't disambiguate
    # (multiple distinct prefixes in the same dir): append @PREFIX.
    python visualize_results.py \\
        --mfem "/path/to/run_a@tpv102_mfadj_p1_O2" \\
        --drdg3d --save

    # Default benchmark directory is tpv102/benchmark_data/{scec_drdg3d,scec_pylith}
    # relative to this script.  Override with --benchmark-dir.

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --benchmarks --save
"""

import argparse
import os
import sys

import numpy as np

# SCEC TPV102 on-fault stations: (mfem_name, x2_km, x3_km).
# The MFEM filename suffix encodes negative along-strike as 'n' (e.g.,
# flt_n12_3.dat for x2 = -12 km).  Benchmark filenames use the
# conventional "-12" minus sign.  Both use a decimal point for 7.5.
SCEC_STATIONS = [
    ("flt_0_3",     0,    3),
    ("flt_0_7.5",   0,    7.5),
    ("flt_0_12",    0,   12),
    ("flt_9_7.5",   9,    7.5),
    ("flt_12_3",   12,    3),
    ("flt_12_12",  12,   12),
    ("flt_n9_7.5", -9,    7.5),
    ("flt_n12_3", -12,    3),
    ("flt_n12_12",-12,   12),
]

PA_TO_MPA = 1.0e-6


def _parse_numeric_table(filepath, min_cols=9):
    """Load a whitespace-delimited numeric table, skipping '#' comments
    and SCEC-style header lines that begin with a letter.
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
                continue
    if not rows:
        return None
    return np.array(rows)


def load_mfem_file(filepath):
    """Load MFEM TPV102 fault-station file.

    Columns in file (TPV102StationWriter convention; BP5 frame where
    1 = dip and 2 = strike):
      0: time (s)
      1: slip1 (m)      == dip slip
      2: slip2 (m)      == strike slip
      3: V1 (m/s)       == dip slip rate
      4: V2 (m/s)       == strike slip rate
      5: tau1 (Pa)      == dip shear stress
      6: tau2 (Pa)      == strike shear stress
      7: sigma_n (Pa)   (positive = compression in MFEM's convention)
      8: log10(theta)

    Returns dict in a unified (strike/dip, stresses in MPa) schema so
    plot_station / plot_overview can compare against the SCEC references
    without further remapping.
    """
    arr = _parse_numeric_table(filepath, min_cols=9)
    if arr is None:
        return None
    return {
        "time_s":      arr[:, 0],
        "slip_strike": arr[:, 2],
        "slip_dip":    arr[:, 1],
        "V_strike":    arr[:, 4],
        "V_dip":       arr[:, 3],
        "tau_strike":  arr[:, 6] * PA_TO_MPA,
        "tau_dip":     arr[:, 5] * PA_TO_MPA,
        "sigma_n":     np.abs(arr[:, 7]) * PA_TO_MPA,
        "log10_theta": arr[:, 8],
    }


def load_scec_file(filepath):
    """Load SCEC TPV102 reference file (DRDG3D or PyLith).

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
        "time_s":      arr[:, 0],
        "slip_strike": arr[:, 1],
        "slip_dip":    arr[:, 4],
        "V_strike":    arr[:, 2],
        "V_dip":       arr[:, 5],
        "tau_strike":  arr[:, 3],
        "tau_dip":     arr[:, 6],
        "sigma_n":     np.abs(arr[:, 7]),
        "log10_theta": arr[:, 8],
    }


def load_drdg3d_file(filepath):
    """Load DRDG3D TPV102 reference file.

    DRDG3D's column layout matches the SCEC standard exactly; wrapper
    exists for symmetry with load_mfem_file / load_pylith_file in case
    the format ever diverges.
    """
    return load_scec_file(filepath)


def load_pylith_file(filepath):
    """Load PyLith TPV102 reference file (same SCEC standard layout)."""
    return load_scec_file(filepath)


def mfem_coord_str(val_km):
    """MFEM filename coord: -12 -> 'n12', 0 -> '0', 7.5 -> '7.5'."""
    s = f"{val_km:g}"
    if s.startswith("-"):
        s = "n" + s[1:]
    return s


def bench_coord_str(val_km):
    """Benchmark filename coord: -12 -> '-12', 0 -> '0', 7.5 -> '7.5'."""
    return f"{val_km:g}"


def mfem_filename(results_dir, prefix, x2_km, x3_km):
    """Build the MFEM station filename for this prefix and station coords.

    Pattern: <prefix>_station_flt_<mfem_coord(x2)>_<mfem_coord(x3)>.dat
    e.g. tpv102_station_flt_n12_3.dat or tpv102_mfadj_p1_O2_station_flt_0_7.5.dat.
    """
    suffix = f"flt_{mfem_coord_str(x2_km)}_{mfem_coord_str(x3_km)}"
    fname = f"{prefix}_station_{suffix}.dat"
    return os.path.join(results_dir, fname)


def drdg3d_filename(bench_dir, x2_km, x3_km):
    fname = (f"tpv102_drdg3d_x2_{bench_coord_str(x2_km)}_"
             f"x3_{bench_coord_str(x3_km)}.txt")
    return os.path.join(bench_dir, fname)


def pylith_filename(bench_dir, x2_km, x3_km):
    fname = (f"tpv102_pylith_x2_{bench_coord_str(x2_km)}_"
             f"x3_{bench_coord_str(x3_km)}.txt")
    return os.path.join(bench_dir, fname)


PANELS = [
    ("V_strike",    "Slip Rate V_strike (m/s)",       False),
    ("slip_strike", "Slip Strike (m)",                 False),
    ("tau_strike",  "Shear Stress τ_strike (MPa)", False),
    ("V_dip",       "Slip Rate V_dip (m/s)",           False),
    ("slip_dip",    "Slip Dip (m)",                    False),
    ("tau_dip",     "Shear Stress τ_dip (MPa)",   False),
    ("sigma_n",     "|Normal Stress| (MPa)",           False),
    ("log10_theta", "log₁₀(State) (s)",      False),
]


def plot_station(datasets, station_name, x2_km, x3_km, save_path=None,
                 t_max=None):
    """Plot 8-panel station comparison."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    title = f"TPV102: {station_name}  (x2={x2_km} km, x3={x3_km} km)"
    if t_max is not None:
        title += f"   [0–{t_max:g} s close-up]"
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
            ax.plot(t, y, ls, color=color, label=label, linewidth=0.9,
                    alpha=0.85)
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


def parse_mfem_spec(spec):
    """Parse '[LABEL:]DIR[@PREFIX]' -> (label_or_None, dir, prefix_or_None).

    Splits the optional '@PREFIX' suffix from the DIR before applying the
    label parser, so a per-source prefix can be carried alongside the
    directory.  Filesystem paths are not expected to contain '@' on the
    target platforms; if they do, pass DIR via --mfem-prefix and skip the
    suffix.
    """
    label, value = parse_labeled_arg(spec)
    if "@" in value:
        directory, prefix = value.rsplit("@", 1)
    else:
        directory, prefix = value, None
    return label, directory, prefix


def detect_mfem_prefix(directory):
    """Auto-detect the station-file prefix used in `directory`.

    Scans for files matching '*_station_flt_*.dat' and extracts the
    common stem (everything before '_station_').  Returns the prefix
    string when exactly one prefix is found; returns None when zero
    matches OR when multiple distinct prefixes coexist (the caller
    should fall back to --mfem-prefix in that case).
    """
    import glob
    if not os.path.isdir(directory):
        return None
    files = glob.glob(os.path.join(directory, "*_station_flt_*.dat"))
    prefixes = set()
    for f in files:
        name = os.path.basename(f)
        if "_station_" not in name:
            continue
        prefixes.add(name.split("_station_", 1)[0])
    if len(prefixes) == 1:
        return next(iter(prefixes))
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Visualize TPV102 output: MFEM vs SCEC reference data"
    )
    parser.add_argument(
        "mfem_dir_positional",
        nargs="?",
        default=None,
        help="(Legacy) MFEM results directory containing "
             "<prefix>_station_flt_*.dat",
    )
    parser.add_argument(
        "--mfem",
        action="append",
        metavar="[LABEL:]DIR[@PREFIX]",
        help="MFEM results directory.  Use 'label:dir' for a custom legend "
             "label, '@prefix' to force a specific filename prefix.  "
             "Repeatable.",
    )
    parser.add_argument(
        "--mfem-prefix",
        default="tpv102",
        help="MFEM file prefix matching --output-prefix passed to "
             "seas_tpv102_driver (default: 'tpv102').  Used as the "
             "fallback when the per-source @PREFIX is absent and "
             "auto-detection finds zero or multiple candidate prefixes.",
    )

    parser.add_argument(
        "--drdg3d", action="store_true",
        help="Include DRDG3D (Wenqiang Zhang) reference data.",
    )
    parser.add_argument(
        "--pylith", action="store_true",
        help="Include PyLith (Brad Aagaard) reference data.",
    )
    parser.add_argument(
        "--benchmarks", action="store_true",
        help="Shortcut: enable both --drdg3d and --pylith.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing scec_drdg3d/ and scec_pylith/ "
             "subfolders with reference traces.  Resolved relative to "
             "this script if not absolute.",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip ALL benchmark references (DRDG3D, PyLith) even if "
             "--drdg3d / --pylith / --benchmarks is given.",
    )

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
        help="Additionally produce a close-up plot truncated to this time [s]",
    )

    args = parser.parse_args()

    if args.benchmarks:
        args.drdg3d = True
        args.pylith = True

    # Default to including BOTH benchmarks when no benchmark flag was
    # given and --no-benchmark wasn't requested.  If the user opted in
    # to ANY single benchmark (e.g. --drdg3d alone), the default-on rule
    # does NOT fire — they get only what they asked for.
    if (
        not args.no_benchmark
        and not args.drdg3d
        and not args.pylith
    ):
        args.drdg3d = True
        args.pylith = True

    # Build ordered sources from command-line argv so legend colors match
    # the order the user typed.
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

    # If a benchmark was implicit (default), append it at the end so it
    # still appears in the typed-source list when programmatic callers
    # never put the flag in sys.argv.
    if (
        args.drdg3d
        and not args.no_benchmark
        and not any(s == "drdg3d" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("drdg3d", None))
    if (
        args.pylith
        and not args.no_benchmark
        and not any(s == "pylith" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("pylith", None))

    if not ordered_sources:
        parser.error(
            "No data sources specified.  Provide --mfem DIR and/or "
            "--drdg3d / --pylith / --benchmarks."
        )

    try:
        import matplotlib
        if args.save:
            matplotlib.use("Agg")
    except ImportError:
        print("Error: matplotlib required. Install with: pip install matplotlib")
        return 1

    # Resolve benchmark directory.
    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), args.benchmark_dir
    )
    if not os.path.isabs(args.benchmark_dir) and not os.path.isdir(data_dir):
        data_dir = args.benchmark_dir
    drdg3d_dir = os.path.join(data_dir, "scec_drdg3d")
    pylith_dir = os.path.join(data_dir, "scec_pylith")

    # Select stations.
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations
            if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build the typed-source list with colors and line styles.
    benchmark_types = {"drdg3d", "pylith"}
    benchmark_labels = {
        "drdg3d": "DRDG3D (ref)",
        "pylith": "PyLith (ref)",
    }
    benchmark_linestyle = {"drdg3d": "--", "pylith": ":"}
    sources = []
    ci = 0
    for stype, spec in ordered_sources:
        if args.no_benchmark and stype in benchmark_types:
            continue
        color = COLORS[ci % len(COLORS)]
        ci += 1
        if stype in benchmark_types:
            sources.append(
                (benchmark_labels[stype], stype, None, color,
                 benchmark_linestyle[stype])
            )
        elif stype == "mfem":
            label, directory, prefix = parse_mfem_spec(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            # Prefix-resolution priority:
            #   1. explicit '@PREFIX' in the spec
            #   2. auto-detected unique prefix in <directory>/*_station_flt_*.dat
            #   3. --mfem-prefix CLI flag (default 'tpv102')
            if prefix is None:
                prefix = detect_mfem_prefix(directory)
            if prefix is None:
                prefix = args.mfem_prefix
            sources.append(
                (f"MFEM {label}", "mfem",
                 (directory, prefix), color, "-")
            )

    print("=" * 60)
    print("TPV102 Visualization")
    print("=" * 60)
    style_name = {"-": "solid", "--": "dashed", ":": "dotted"}
    for label, stype, info, _color, ls in sources:
        style = style_name.get(ls, ls)
        if stype == "drdg3d":
            print(f"  [{style}] {label}: {drdg3d_dir}/")
        elif stype == "pylith":
            print(f"  [{style}] {label}: {pylith_dir}/")
        else:
            directory, prefix = info
            print(f"  [{style}] {label}: {directory}/{prefix}_station_flt_*.dat")
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
                    data = load_drdg3d_file(path)
            elif stype == "pylith":
                path = pylith_filename(pylith_dir, x2_km, x3_km)
                if os.path.exists(path):
                    data = load_pylith_file(path)
            elif stype == "mfem":
                directory, prefix = info
                path = mfem_filename(directory, prefix, x2_km, x3_km)
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
        print(f"  {station_name} (x2={x2_km}, x3={x3_km}): "
              f"{', '.join(pts_info)}")

        all_results.append({
            "station_name": station_name,
            "x2_km": x2_km,
            "x3_km": x3_km,
            "datasets": datasets,
        })

        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(
                args.output_dir, f"tpv102_{station_name}.png"
            )
            plot_station(
                datasets, station_name, x2_km, x3_km, save_path=fname
            )
            if args.closeup_t is not None:
                fname_close = os.path.join(
                    args.output_dir,
                    f"tpv102_{station_name}_closeup.png"
                )
                plot_station(
                    datasets, station_name, x2_km, x3_km,
                    save_path=fname_close, t_max=args.closeup_t,
                )
        else:
            plot_station(datasets, station_name, x2_km, x3_km)
            if args.closeup_t is not None:
                plot_station(
                    datasets, station_name, x2_km, x3_km,
                    t_max=args.closeup_t,
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
