#!/usr/bin/env python3
"""
TPV205 (SCEC TPV5) Benchmark Visualization Script

Plots MFEM dynamic-rupture station output alongside the DRDG3D
reference traces (Wenqiang Zhang, on-fault 200 m / 100 m at order 4)
for the SCEC TPV5 linear-slip-weakening benchmark.  16 on-fault
stations are compared.

MFEM station format (set in dynamic/tpv205_setup.hpp::TPV205StationWriter):
    <output_prefix>_station_<x2_x3_label>.dat
columns (LSW; final column is μ_eff(δ), not ψ as in TPV104):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(Pa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(Pa)
            n-stress(Pa) mu_eff
where h = horizontal = strike (BP5 component 2 in MFEM's convention),
v = vertical = dip (BP5 component 1).  TPV205 is right-lateral
strike-slip so v-* channels remain near 0.  The MFEM writer emits the
SCEC column order — no in-Python remapping needed.

DRDG3D reference format (benchmark_data/DRDG3D_{200m,100m}_O4/
tpv205_drdg3d_<x2_x3_label>.txt):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(MPa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(MPa)
The DRDG3D file has only 7 columns; n-stress and mu_eff are NOT
provided.  Those panels are blank for the reference dataset; the MFEM
trace fills them.

Usage:
    # Single MFEM run vs the 200 m DRDG3D reference (default).  Prefix
    # is auto-detected from station files in the directory:
    python visualize_results.py \\
        --mfem /path/to/results_dir --save

    # Compare TWO MFEM runs (each prefix auto-detected; directory
    # basename becomes the legend label):
    python visualize_results.py \\
        --mfem /path/to/run_a \\
        --mfem /path/to/run_b \\
        --save

    # Compare against BOTH reference resolutions (200 m and 100 m):
    python visualize_results.py \\
        --mfem /path/to/run --drdg3d-200m --drdg3d-100m --save

    # Per-source label + forced filename prefix:
    python visualize_results.py \\
        --mfem "MF=adjacent:/path/to/run_a@tpv205_mfadj_p1_O2" \\
        --save

    # Default benchmark directory is tpv205/benchmark_data/ relative to
    # this script.  Override with --benchmark-dir.

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --save
"""

import argparse
import os
import sys

import numpy as np

# SCEC TPV5 / TPV205 on-fault stations: (mfem_label, x2_km, x3_km).
# Label matches both the MFEM TPV205StationWriter (via kStationsTPV205
# in config/tpv205_params.hpp) and the DRDG3D reference filenames in
# benchmark_data/DRDG3D_{200m,100m}_O4/.
SCEC_STATIONS = [
    ("x2_-12_x3_0",   -12.0,  0.0),
    ("x2_-12_x3_7.5", -12.0,  7.5),
    ("x2_-7.5_x3_0",   -7.5,  0.0),
    ("x2_-7.5_x3_7.5", -7.5,  7.5),
    ("x2_-4.5_x3_0",   -4.5,  0.0),
    ("x2_-4.5_x3_7.5", -4.5,  7.5),
    ("x2_0_x3_0",       0.0,  0.0),
    ("x2_0_x3_3",       0.0,  3.0),
    ("x2_0_x3_7.5",     0.0,  7.5),
    ("x2_0_x3_12",      0.0, 12.0),
    ("x2_4.5_x3_0",     4.5,  0.0),
    ("x2_4.5_x3_7.5",   4.5,  7.5),
    ("x2_7.5_x3_0",     7.5,  0.0),
    ("x2_7.5_x3_7.5",   7.5,  7.5),
    ("x2_12_x3_0",     12.0,  0.0),
    ("x2_12_x3_7.5",   12.0,  7.5),
]

PA_TO_MPA = 1.0e-6


def _parse_numeric_table(filepath, min_cols):
    """Load a whitespace-delimited numeric table, skipping '#' comments
    and DRDG3D-style header lines that begin with a letter.
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
                # DRDG3D names its columns on a header line that begins
                # with `t h-slip ...` — skip it.
                continue
    if not rows:
        return None
    return np.array(rows)


def load_mfem_file(filepath):
    """Load MFEM TPV205 fault-station file.

    MFEM writes columns in SCEC TPV5 order (per TPV205StationWriter)
    with stresses in **Pa** instead of MPa, and a final μ_eff column
    in place of TPV104's ψ:
      0: time (s)
      1: h-slip (m)        == strike slip
      2: h-slip-rate (m/s) == strike slip rate
      3: h-shear (Pa)      == strike shear stress
      4: v-slip (m)        == dip slip
      5: v-slip-rate (m/s) == dip slip rate
      6: v-shear (Pa)      == dip shear stress
      7: n-stress (Pa)     (positive = compression in MFEM's convention)
      8: mu_eff            (LSW μ(δ); 0.677 → 0.525 over d_c = 0.4 m)

    Returns dict with stresses converted to MPa for plotting.
    """
    arr = _parse_numeric_table(filepath, min_cols=9)
    if arr is None:
        return None
    return {
        "time_s":      arr[:, 0],
        "slip_strike": arr[:, 1],
        "V_strike":    arr[:, 2],
        "tau_strike":  arr[:, 3] * PA_TO_MPA,
        "slip_dip":    arr[:, 4],
        "V_dip":       arr[:, 5],
        "tau_dip":     arr[:, 6] * PA_TO_MPA,
        "sigma_n":     np.abs(arr[:, 7]) * PA_TO_MPA,
        "mu_eff":      arr[:, 8],
    }


def load_drdg3d_file(filepath):
    """Load DRDG3D TPV205 reference station file.

    DRDG3D's TPV5 trace has 7 columns (no n-stress, no μ_eff):
      0: time (s)
      1: h-slip (m)            == strike slip
      2: h-slip-rate (m/s)     == strike slip rate
      3: h-shear-stress (MPa)  == strike shear stress
      4: v-slip (m)            == dip slip
      5: v-slip-rate (m/s)     == dip slip rate
      6: v-shear-stress (MPa)  == dip shear stress
    Stresses are already in MPa.
    """
    arr = _parse_numeric_table(filepath, min_cols=7)
    if arr is None:
        return None
    n = arr.shape[0]
    return {
        "time_s":      arr[:, 0],
        "slip_strike": arr[:, 1],
        "V_strike":    arr[:, 2],
        "tau_strike":  arr[:, 3],
        "slip_dip":    arr[:, 4],
        "V_dip":       arr[:, 5],
        "tau_dip":     arr[:, 6],
        # DRDG3D does not report σ_n or μ_eff for TPV5; leave blank so
        # the corresponding panels show only the MFEM trace.
        "sigma_n":     np.full(n, np.nan),
        "mu_eff":      np.full(n, np.nan),
    }


def mfem_filename(results_dir, prefix, label):
    fname = f"{prefix}_station_{label}.dat"
    return os.path.join(results_dir, fname)


def drdg3d_filename(bench_dir, label):
    fname = f"tpv205_drdg3d_{label}.txt"
    return os.path.join(bench_dir, fname)


PANELS = [
    ("V_strike",    "Slip Rate V_strike (m/s)",      False),
    ("slip_strike", "Slip Strike (m)",                False),
    ("tau_strike",  "Shear Stress τ_strike (MPa)", False),
    ("V_dip",       "Slip Rate V_dip (m/s)",          False),
    ("slip_dip",    "Slip Dip (m)",                   False),
    ("tau_dip",     "Shear Stress τ_dip (MPa)",  False),
    # Normal stress: MFEM-only channel (DRDG3D's 7-column TPV5 trace does
    # not report σ_n, so its sigma_n is all-NaN and is skipped by the
    # all-NaN guard in plot_station — this panel shows only the MFEM run).
    ("sigma_n",     "Normal Stress σ_n (MPa)",       False),
]


def plot_station(datasets, station_label, x2_km, x3_km, save_path=None,
                 t_max=None):
    """Plot the per-station panel comparison (one panel per PANELS entry).

    Panels: V/slip/τ for the strike and dip components, plus normal
    stress σ_n.  σ_n is an MFEM-only channel — the DRDG3D TPV5 reference
    is a 7-column file with no σ_n, so its trace is all-NaN and skipped
    by the all-NaN guard below (the σ_n panel shows only the MFEM run).
    The LSW μ_eff channel remains omitted for the same reason (DRDG3D
    does not report it and it adds no cross-code comparison).

    The grid is sized to hold len(PANELS) panels in 2 columns; any unused
    trailing cell is hidden.
    """
    import matplotlib.pyplot as plt

    ncols = 2
    nrows = (len(PANELS) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows))
    title = (f"TPV205: {station_label}  "
             f"(x2={x2_km} km, x3={x3_km} km)")
    if t_max is not None:
        title += f"   [0–{t_max:g} s close-up]"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    for ax, (key, ylabel, _use_log) in zip(axes.flat, PANELS):
        for label, data, color, ls in datasets:
            if data is None:
                continue
            t = data["time_s"]
            y = data[key]
            # Skip channels that are all-NaN (DRDG3D σ_n / μ_eff).
            if np.all(np.isnan(y)):
                continue
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

    # Hide any trailing cells not backed by a PANELS entry (e.g. the 8th
    # cell when there are 7 panels in a 4x2 grid).
    for ax in axes.flat[len(PANELS):]:
        ax.axis("off")

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
        "TPV205: Slip Rate V_strike at All Stations",
        fontsize=14,
        fontweight="bold",
    )

    cmap = plt.cm.tab20
    n = len(all_results)
    for i, res in enumerate(all_results):
        name = res["station_label"]
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
    target platforms; if they do, pass DIR via --mfem-prefix and skip
    the suffix.
    """
    label, value = parse_labeled_arg(spec)
    if "@" in value:
        directory, prefix = value.rsplit("@", 1)
    else:
        directory, prefix = value, None
    return label, directory, prefix


def detect_mfem_prefix(directory):
    """Auto-detect the station-file prefix used in `directory`.

    Scans for files matching '*_station_x2_*_x3_*.dat' and extracts the
    common stem (everything before '_station_').  Returns the prefix
    string when exactly one prefix is found; returns None when zero
    matches OR when multiple distinct prefixes coexist (the caller
    should fall back to --mfem-prefix in that case).
    """
    import glob
    if not os.path.isdir(directory):
        return None
    files = glob.glob(os.path.join(directory, "*_station_x2_*_x3_*.dat"))
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
        description="Visualize TPV205 output: MFEM vs DRDG3D reference"
    )
    parser.add_argument(
        "mfem_dir_positional",
        nargs="?",
        default=None,
        help="(Legacy) MFEM results directory containing "
             "<prefix>_station_x2_*_x3_*.dat",
    )
    parser.add_argument(
        "--mfem",
        action="append",
        metavar="[LABEL:]DIR[@PREFIX]",
        help="MFEM results directory.  Use 'label:dir' for a custom "
             "legend label, '@prefix' to force a specific filename "
             "prefix.  Repeatable.",
    )
    parser.add_argument(
        "--mfem-prefix",
        default="tpv205",
        help="MFEM file prefix matching --output-prefix passed to "
             "seas_tpv205_driver (default: 'tpv205').  Used as the "
             "fallback when the per-source @PREFIX is absent and "
             "auto-detection finds zero or multiple candidate prefixes.",
    )

    parser.add_argument(
        "--drdg3d-200m", action="store_true",
        help="Include DRDG3D reference at 200 m / O4 resolution.",
    )
    parser.add_argument(
        "--drdg3d-100m", action="store_true",
        help="Include DRDG3D reference at 100 m / O4 resolution.",
    )
    parser.add_argument(
        "--drdg3d", action="store_true",
        help="Shortcut: enable both --drdg3d-200m and --drdg3d-100m.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing DRDG3D_{200m,100m}_O4/ subfolders.  "
             "Resolved relative to this script if not absolute.",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip ALL benchmark references (DRDG3D 200 m, DRDG3D 100 m)"
             " even if --drdg3d-* / --drdg3d is given.",
    )

    parser.add_argument(
        "--stations", nargs="+", type=int, default=None,
        help="Specific station indices (1-16). Default: all",
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

    if args.drdg3d:
        args.drdg3d_200m = True
        args.drdg3d_100m = True

    # Default to including the 200 m DRDG3D reference when no benchmark
    # flag is given and --no-benchmark wasn't requested.  If the user
    # opted in to ANY single benchmark, the default-on rule does NOT
    # fire — they get only what they asked for.
    if (
        not args.no_benchmark
        and not args.drdg3d_200m
        and not args.drdg3d_100m
    ):
        args.drdg3d_200m = True

    # Build ordered sources from command-line argv so legend colors
    # match the order the user typed.
    ordered_sources = []
    if args.mfem_dir_positional:
        ordered_sources.append(("mfem", args.mfem_dir_positional))

    mfem_iter = iter(args.mfem or [])
    for arg in sys.argv[1:]:
        if arg == "--drdg3d":
            ordered_sources.append(("drdg3d-200m", None))
            ordered_sources.append(("drdg3d-100m", None))
        elif arg == "--drdg3d-200m":
            ordered_sources.append(("drdg3d-200m", None))
        elif arg == "--drdg3d-100m":
            ordered_sources.append(("drdg3d-100m", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))

    # If a benchmark was implicit (default), append it at the end so it
    # still appears in the typed-source list when programmatic callers
    # never put the flag in sys.argv.
    if (
        args.drdg3d_200m
        and not args.no_benchmark
        and not any(s == "drdg3d-200m" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("drdg3d-200m", None))
    if (
        args.drdg3d_100m
        and not args.no_benchmark
        and not any(s == "drdg3d-100m" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("drdg3d-100m", None))

    if not ordered_sources:
        parser.error(
            "No data sources specified.  Provide --mfem DIR and/or "
            "--drdg3d-200m / --drdg3d-100m / --drdg3d."
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
    drdg3d_200m_dir = os.path.join(data_dir, "DRDG3D_200m_O4")
    drdg3d_100m_dir = os.path.join(data_dir, "DRDG3D_100m_O4")

    # Select stations.
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations
            if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build the typed-source list with colors and line styles.
    benchmark_types = {"drdg3d-200m", "drdg3d-100m"}
    benchmark_labels = {
        "drdg3d-200m": "DRDG3D 200 m O4",
        "drdg3d-100m": "DRDG3D 100 m O4",
    }
    benchmark_dirs = {
        "drdg3d-200m": drdg3d_200m_dir,
        "drdg3d-100m": drdg3d_100m_dir,
    }
    benchmark_linestyle = {"drdg3d-200m": "--", "drdg3d-100m": ":"}
    sources = []
    ci = 0
    for stype, spec in ordered_sources:
        if args.no_benchmark and stype in benchmark_types:
            continue
        color = COLORS[ci % len(COLORS)]
        ci += 1
        if stype in benchmark_types:
            sources.append(
                (benchmark_labels[stype], stype, benchmark_dirs[stype],
                 color, benchmark_linestyle[stype])
            )
        elif stype == "mfem":
            label, directory, prefix = parse_mfem_spec(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            # Prefix-resolution priority:
            #   1. explicit '@PREFIX' in the spec
            #   2. auto-detected unique prefix in
            #      <directory>/*_station_x2_*_x3_*.dat
            #   3. --mfem-prefix CLI flag (default 'tpv205')
            if prefix is None:
                prefix = detect_mfem_prefix(directory)
            if prefix is None:
                prefix = args.mfem_prefix
            sources.append(
                (f"MFEM {label}", "mfem",
                 (directory, prefix), color, "-")
            )

    print("=" * 60)
    print("TPV205 Visualization")
    print("=" * 60)
    style_name = {"-": "solid", "--": "dashed", ":": "dotted"}
    for label, stype, info, _color, ls in sources:
        style = style_name.get(ls, ls)
        if stype in benchmark_types:
            print(f"  [{style}] {label}: {info}/")
        else:
            directory, prefix = info
            print(f"  [{style}] {label}: "
                  f"{directory}/{prefix}_station_x2_*_x3_*.dat")
    print(f"  Stations: {len(stations)}")
    print()

    all_results = []
    for station_label, x2_km, x3_km in stations:
        datasets = []
        for label, stype, info, color, ls in sources:
            data = None
            if stype in benchmark_types:
                bench_dir = info
                path = drdg3d_filename(bench_dir, station_label)
                if os.path.exists(path):
                    data = load_drdg3d_file(path)
            elif stype == "mfem":
                directory, prefix = info
                path = mfem_filename(directory, prefix, station_label)
                if os.path.exists(path):
                    data = load_mfem_file(path)
            datasets.append((label, data, color, ls))

        if all(d is None for _, d, _, _ in datasets):
            print(f"  {station_label}: no data found")
            continue

        pts_info = []
        for label, data, _c, _ls in datasets:
            if data is not None:
                t_last = data["time_s"][-1]
                pts_info.append(
                    f"{label}: {len(data['time_s'])} pts ({t_last:.2f} s)"
                )
        print(f"  {station_label} (x2={x2_km}, x3={x3_km}): "
              f"{', '.join(pts_info)}")

        all_results.append({
            "station_label": station_label,
            "x2_km": x2_km,
            "x3_km": x3_km,
            "datasets": datasets,
        })

        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(
                args.output_dir, f"tpv205_{station_label}.png"
            )
            plot_station(
                datasets, station_label, x2_km, x3_km, save_path=fname
            )
            if args.closeup_t is not None:
                fname_close = os.path.join(
                    args.output_dir,
                    f"tpv205_{station_label}_closeup.png"
                )
                plot_station(
                    datasets, station_label, x2_km, x3_km,
                    save_path=fname_close, t_max=args.closeup_t,
                )
        else:
            plot_station(datasets, station_label, x2_km, x3_km)
            if args.closeup_t is not None:
                plot_station(
                    datasets, station_label, x2_km, x3_km,
                    t_max=args.closeup_t,
                )

    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "tpv205_overview.png")
            plot_overview(all_results, save_path=fname)
        else:
            plot_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
