#!/usr/bin/env python3
"""
TPV104 Benchmark Visualization Script

Plots MFEM dynamic-rupture station output alongside the SeisSol reference
(SCEC TPV104 / FL=103, slip law with strong rate weakening). Nine on-fault
stations are compared.

MFEM station format (set in dynamic/tpv104_setup.hpp::TPV104StationWriter):
    <output_prefix>_station_x2_<X>_x3_<Y>.dat
columns:
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(Pa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(Pa)
            n-stress(Pa) psi
where h = horizontal = strike (BP5 component 2 in MFEM's convention),
v = vertical = dip (BP5 component 1).  TPV104 is pure strike-slip so the
v-* channels remain ~0.  The MFEM writer already emits the SCEC column
order — no in-Python remapping needed.

SeisSol reference format (benchmark_data/seisol/tpv104_seisol_x2_<X>_x3_<Y>.txt):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(MPa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(MPa)
            n-stress(MPa) psi
identical column order; stresses are in MPa instead of Pa.

DRDG3D reference format (benchmark_data/DRDG3D/tpv104_drdg3d_x2_<X>_x3_<Y>.txt):
    SCEC TPV104 / FL=103 reference produced by Wenqiang Zhang's DRDG3D
    (DG-on-fault, 200 m, O4).  Same 9-column layout as SeisSol with
    stresses already in MPa.  Opt in with --drdg3d.

Usage:
    # MFEM results vs SeisSol + DRDG3D references, save plots
    python visualize_results.py \\
        --mfem /path/to/results_dir --seisol --drdg3d --save

    # Compare two MFEM runs against both references
    python visualize_results.py \\
        --mfem "p1 1000m:/path/to/run_a" \\
        --mfem "p1 500m:/path/to/run_b" \\
        --seisol --drdg3d --save

    # Default benchmark directory is tpv104/benchmark_data/{seisol,DRDG3D}
    # relative to this script.  Override with --benchmark-dir.

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --seisol --save
"""

import argparse
import os
import sys

import numpy as np

# SCEC TPV104 on-fault stations: (mfem_label, x2_km, x3_km).
# The label is the suffix used both by MFEM's TPV104StationWriter (via
# kStationsTPV104 in config/tpv104_params.hpp) and by the SeisSol reference
# filenames in benchmark_data/seisol/.
SCEC_STATIONS = [
    ("x2_0_x3_3",     0,    3),
    ("x2_0_x3_7.5",   0,    7.5),
    ("x2_0_x3_12",    0,   12),
    ("x2_9_x3_7.5",   9,    7.5),
    ("x2_12_x3_3",   12,    3),
    ("x2_12_x3_12",  12,   12),
    ("x2_-9_x3_7.5", -9,    7.5),
    ("x2_-12_x3_3", -12,    3),
    ("x2_-12_x3_12",-12,   12),
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
    """Load MFEM TPV104 station file.

    MFEM writes columns in SCEC order (per TPV104StationWriter), but with
    stresses in **Pa** instead of MPa:
      0: time (s)
      1: h-slip (m)        == strike slip
      2: h-slip-rate (m/s) == strike slip rate
      3: h-shear (Pa)      == strike shear stress
      4: v-slip (m)        == dip slip
      5: v-slip-rate (m/s) == dip slip rate
      6: v-shear (Pa)      == dip shear stress
      7: n-stress (Pa)     (positive = compression in MFEM's convention)
      8: psi (linear, not log10)

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
        "psi":         arr[:, 8],
    }


def load_seisol_file(filepath):
    """Load SeisSol TPV104 reference station file.

    Columns:
      0: time (s)
      1: h-slip (m)            == strike slip
      2: h-slip-rate (m/s)     == strike slip rate
      3: h-shear-stress (MPa)  == strike shear stress
      4: v-slip (m)            == dip slip
      5: v-slip-rate (m/s)     == dip slip rate
      6: v-shear-stress (MPa)  == dip shear stress
      7: n-stress (MPa)        (sign convention varies; plotted as |.|)
      8: psi (linear)
    """
    arr = _parse_numeric_table(filepath, min_cols=9)
    if arr is None:
        return None
    return {
        "time_s":      arr[:, 0],
        "slip_strike": arr[:, 1],
        "V_strike":    arr[:, 2],
        "tau_strike":  arr[:, 3],
        "slip_dip":    arr[:, 4],
        "V_dip":       arr[:, 5],
        "tau_dip":     arr[:, 6],
        "sigma_n":     np.abs(arr[:, 7]),
        "psi":         arr[:, 8],
    }


def load_drdg3d_file(filepath):
    """Load DRDG3D TPV104 reference station file.

    DRDG3D's column layout matches SeisSol's exactly (9 columns, stresses
    in MPa, same SCEC ordering).  Wrapper exists for symmetry with
    load_mfem_file / load_seisol_file in case the format ever diverges.
    """
    return load_seisol_file(filepath)


def mfem_filename(results_dir, prefix, label):
    fname = f"{prefix}_station_{label}.dat"
    return os.path.join(results_dir, fname)


def seisol_filename(bench_dir, label):
    fname = f"tpv104_seisol_{label}.txt"
    return os.path.join(bench_dir, fname)


def drdg3d_filename(bench_dir, label):
    fname = f"tpv104_drdg3d_{label}.txt"
    return os.path.join(bench_dir, fname)


PANELS = [
    ("V_strike",    "Slip Rate V_strike (m/s)",      False),
    ("slip_strike", "Slip Strike (m)",                False),
    ("tau_strike",  "Shear Stress τ_strike (MPa)", False),
    ("V_dip",       "Slip Rate V_dip (m/s)",          False),
    ("slip_dip",    "Slip Dip (m)",                   False),
    ("tau_dip",     "Shear Stress τ_dip (MPa)",   False),
    ("sigma_n",     "|Normal Stress| (MPa)",           False),
    ("psi",         "State Variable ψ",           False),
]


def plot_station(datasets, station_label, x2_km, x3_km, save_path=None,
                 t_max=None):
    """Plot 8-panel station comparison."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    title = f"TPV104: {station_label}  (x2={x2_km} km, x3={x3_km} km)"
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
        "TPV104: Slip Rate V_strike at All Stations",
        fontsize=14,
        fontweight="bold",
    )

    cmap = plt.cm.tab10
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


def main():
    parser = argparse.ArgumentParser(
        description="Visualize TPV104 output: MFEM vs SeisSol reference"
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
        metavar="[LABEL:]DIR",
        help="MFEM results directory. Use 'label:dir' for a custom legend "
             "label. Repeatable.",
    )
    parser.add_argument(
        "--mfem-prefix",
        default="tpv104",
        help="MFEM file prefix matching --output-prefix passed to "
             "seas_tpv104_driver (default: 'tpv104').",
    )

    parser.add_argument(
        "--seisol", action="store_true",
        help="Include SeisSol (SCEC TPV104 FL=103) reference data.",
    )
    parser.add_argument(
        "--drdg3d", action="store_true",
        help="Include DRDG3D (Wenqiang Zhang, 200 m, O4) reference data.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing seisol/ and DRDG3D/ subfolders with "
             "reference traces.  Resolved relative to this script if not "
             "absolute.",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip ALL benchmark references (SeisSol, DRDG3D) even if "
             "--seisol or --drdg3d is given.",
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

    # Default to including SeisSol when no benchmark was specified at all,
    # unless --no-benchmark was requested explicitly.  If the user opted in
    # to ANY benchmark (e.g. --drdg3d), the default-on rule does NOT fire —
    # they get only what they asked for.
    if not args.no_benchmark and not args.seisol and not args.drdg3d:
        args.seisol = True

    # Build ordered sources from command-line argv so legend colors match
    # the order the user typed.
    ordered_sources = []
    if args.mfem_dir_positional:
        ordered_sources.append(("mfem", args.mfem_dir_positional))

    mfem_iter = iter(args.mfem or [])
    for arg in sys.argv[1:]:
        if arg == "--seisol":
            ordered_sources.append(("seisol", None))
        elif arg == "--drdg3d":
            ordered_sources.append(("drdg3d", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))

    # If --seisol was implicit (default), append it at the end.
    if (
        args.seisol
        and not args.no_benchmark
        and not any(s == "seisol" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("seisol", None))
    # --drdg3d is opt-in only (no default-on); honor explicit flag here in
    # case argparse parsed it but it never appeared in sys.argv (e.g. when
    # invoked programmatically via main(['--drdg3d', ...])).
    if (
        args.drdg3d
        and not args.no_benchmark
        and not any(s == "drdg3d" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("drdg3d", None))

    if not ordered_sources:
        parser.error(
            "No data sources specified. Provide --mfem DIR and/or --seisol."
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
    seisol_dir = os.path.join(data_dir, "seisol")
    drdg3d_dir = os.path.join(data_dir, "DRDG3D")

    # Select stations.
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations
            if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build the typed-source list with colors and line styles.
    benchmark_types = {"seisol", "drdg3d"}
    benchmark_labels = {
        "seisol": "SeisSol (ref)",
        "drdg3d": "DRDG3D (ref)",
    }
    benchmark_linestyle = {"seisol": "--", "drdg3d": ":"}
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
            label, directory = parse_labeled_arg(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            sources.append((f"MFEM {label}", "mfem", directory, color, "-"))

    print("=" * 60)
    print("TPV104 Visualization")
    print("=" * 60)
    style_name = {"-": "solid", "--": "dashed", ":": "dotted"}
    for label, stype, info, _color, ls in sources:
        style = style_name.get(ls, ls)
        if stype == "seisol":
            print(f"  [{style}] {label}: {seisol_dir}/")
        elif stype == "drdg3d":
            print(f"  [{style}] {label}: {drdg3d_dir}/")
        else:
            print(f"  [{style}] {label}: {info}/{args.mfem_prefix}_station_*.dat")
    print(f"  Stations: {len(stations)}")
    print()

    all_results = []
    for station_label, x2_km, x3_km in stations:
        datasets = []
        for label, stype, info, color, ls in sources:
            data = None
            if stype == "seisol":
                path = seisol_filename(seisol_dir, station_label)
                if os.path.exists(path):
                    data = load_seisol_file(path)
            elif stype == "drdg3d":
                path = drdg3d_filename(drdg3d_dir, station_label)
                if os.path.exists(path):
                    data = load_drdg3d_file(path)
            elif stype == "mfem":
                path = mfem_filename(info, args.mfem_prefix, station_label)
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
                args.output_dir, f"tpv104_{station_label}.png"
            )
            plot_station(
                datasets, station_label, x2_km, x3_km, save_path=fname
            )
            if args.closeup_t is not None:
                fname_close = os.path.join(
                    args.output_dir,
                    f"tpv104_{station_label}_closeup.png"
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
            fname = os.path.join(args.output_dir, "tpv104_overview.png")
            plot_overview(all_results, save_path=fname)
        else:
            plot_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
