#!/usr/bin/env python3
"""
TPV31 (SCEC TPV31) Benchmark Visualization Script

Plots MFEM dynamic-rupture station output alongside the SCEC community
reference traces (EQdyna, Benchun Duan; and SeisSol ADER-DG, Thomas
Ulrich) for the SCEC TPV31 bi-material linear-slip-weakening benchmark.
30 on-fault stations are compared (3 along-strike × 10 down-dip).

MFEM station format (dynamic/tpv31_stations.hpp::TPV31StationWriter):
    <output_prefix>_station_<faultst###dp###>.dat
columns (LSW; final column is μ_eff(δ)):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(Pa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(Pa)
            n-stress(Pa) mu_eff
where h = horizontal = strike (canonical component 2), v = vertical =
down-dip (canonical component 1).  TPV31 is right-lateral strike-slip
so v-* channels remain near 0.  Stresses are in **Pa**; n-stress is
**compression-positive** (MFEM internal / SEAS convention).

SCEC reference format (benchmark_data/scec_{eqdyna,seisol}/
tpv31_{eqdyna,seisol}_x2_<strike_km>_x3_<depth_km>.txt):
    time(s) h-slip(m) h-slip-rate(m/s) h-shear-stress(MPa)
            v-slip(m) v-slip-rate(m/s) v-shear-stress(MPa)
            n-stress(MPa)
Stresses are in **MPa**; n-stress is **compression-negative** (SCEC
extension-positive convention).  The reference does NOT provide μ_eff.

Convention reconciliation for overlay:
  * MFEM stresses Pa → MPa (×1e-6).
  * MFEM n-stress is kept compression-positive; the SCEC reference
    n-stress is NEGATED so both overlay as compression-positive.
  * Shear stresses share a sign convention already (right-lateral
    σ_xz = +30 MPa is positive in both frames — see tpv31.toml).

Usage:
    # Single MFEM run vs both SCEC references (default).  Prefix is
    # auto-detected from station files in the directory:
    python visualize_results.py --mfem /path/to/results_dir --save

    # Compare TWO MFEM runs (each prefix auto-detected; directory
    # basename becomes the legend label):
    python visualize_results.py \\
        --mfem /path/to/run_a \\
        --mfem /path/to/run_b \\
        --save

    # Compare against a single reference code:
    python visualize_results.py --mfem /path/to/run --seisol --save

    # Per-source label + forced filename prefix:
    python visualize_results.py \\
        --mfem "p1O2:/path/to/run_a@tpv31" --save

    # Default benchmark directory is tpv31/benchmark_data/ relative to
    # this script.  Override with --benchmark-dir.

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --save
"""

import argparse
import os
import sys

import numpy as np

PA_TO_MPA = 1.0e-6


def _build_stations():
    """Return the 30 SCEC TPV31 on-fault stations.

    Each entry is (mfem_name, ref_label, strike_km, depth_km) where
      * mfem_name  = 'faultst###dp###' (hundreds-of-metres codes), the
        suffix written by TPV31StationWriter, e.g. 'faultst000dp075'.
      * ref_label  = 'x2_<strike_km>_x3_<depth_km>' (km, trailing zeros
        stripped), the suffix used by the EQdyna/SeisSol reference
        files, e.g. 'x2_0_x3_7.5'.

    Grid: along-strike {0, 6, 12} km × down-dip {0, 0.2, 0.5, 1.0, 2.4,
    3.0, 5.0, 7.5, 10.0, 12.0} km (spec TPV31_32_Description_v03 Part 5).
    """
    strikes_km = [0.0, 6.0, 12.0]
    depths_km = [0.0, 0.2, 0.5, 1.0, 2.4, 3.0, 5.0, 7.5, 10.0, 12.0]
    stations = []
    for s_km in strikes_km:
        for d_km in depths_km:
            st_code = int(round(s_km * 10.0))   # round(metres/100)
            dp_code = int(round(d_km * 10.0))
            mfem_name = f"faultst{st_code:03d}dp{dp_code:03d}"
            ref_label = f"x2_{s_km:g}_x3_{d_km:g}"
            stations.append((mfem_name, ref_label, s_km, d_km))
    return stations


SCEC_STATIONS = _build_stations()


def _parse_numeric_table(filepath, min_cols):
    """Load a whitespace-delimited numeric table, skipping comments and
    SCEC-style header lines (leading '#', possibly indented) and the
    `t h-slip ...` field-name line that begins with a letter.
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
                # SCEC names its columns on a header line that begins
                # with `t h-slip ...` — skip it.
                continue
    if not rows:
        return None
    return np.array(rows)


def load_mfem_file(filepath):
    """Load MFEM TPV31 fault-station file.

    Columns (per TPV31StationWriter), stresses in **Pa**:
      0: time (s)
      1: h-slip (m)        == strike slip
      2: h-slip-rate (m/s) == strike slip rate
      3: h-shear (Pa)      == strike shear stress
      4: v-slip (m)        == dip slip
      5: v-slip-rate (m/s) == dip slip rate
      6: v-shear (Pa)      == dip shear stress
      7: n-stress (Pa)     (compression POSITIVE in MFEM's convention)
      8: mu_eff            (LSW μ(δ); 0.580 → 0.450 over d_c = 0.18 m)

    Returns dict with stresses converted to MPa; n-stress kept
    compression-positive.
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
        "sigma_n":     arr[:, 7] * PA_TO_MPA,
        "mu_eff":      arr[:, 8],
    }


def load_reference_file(filepath):
    """Load a SCEC TPV31 reference station file (EQdyna or SeisSol).

    Columns (8), stresses already in **MPa**:
      0: time (s)
      1: h-slip (m)            == strike slip
      2: h-slip-rate (m/s)     == strike slip rate
      3: h-shear-stress (MPa)  == strike shear stress
      4: v-slip (m)            == dip slip
      5: v-slip-rate (m/s)     == dip slip rate
      6: v-shear-stress (MPa)  == dip shear stress
      7: n-stress (MPa)        (compression NEGATIVE — SCEC convention)

    n-stress is NEGATED so it overlays MFEM's compression-positive
    σ_n.  μ_eff is not provided by the reference; left as NaN so the
    corresponding panel shows only the MFEM trace.
    """
    arr = _parse_numeric_table(filepath, min_cols=8)
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
        # Flip SCEC compression-negative σ_n to compression-positive.
        "sigma_n":     -arr[:, 7],
        # SCEC reference does not report μ_eff; leave blank so the panel
        # shows only the MFEM trace.
        "mu_eff":      np.full(n, np.nan),
    }


def mfem_filename(results_dir, prefix, mfem_name):
    fname = f"{prefix}_station_{mfem_name}.dat"
    return os.path.join(results_dir, fname)


def reference_filename(bench_dir, code, ref_label):
    fname = f"tpv31_{code}_{ref_label}.txt"
    return os.path.join(bench_dir, fname)


PANELS = [
    ("V_strike",    "Slip Rate V_strike (m/s)"),
    ("slip_strike", "Slip Strike (m)"),
    ("tau_strike",  "Shear Stress τ_strike (MPa)"),
    ("V_dip",       "Slip Rate V_dip (m/s)"),
    ("slip_dip",    "Slip Dip (m)"),
    ("tau_dip",     "Shear Stress τ_dip (MPa)"),
    ("sigma_n",     "Normal Stress σ_n (MPa, compression +)"),
    ("mu_eff",      "μ_eff (LSW; MFEM only)"),
]


def plot_station(datasets, station_label, strike_km, depth_km,
                 save_path=None, t_max=None):
    """Plot 8-panel station comparison.

    σ_n is a genuine cross-code panel (both EQdyna and SeisSol provide
    it).  μ_eff is MFEM-only (the SCEC traces omit it), so that panel
    shows the MFEM curve alone — the all-NaN reference channels are
    skipped automatically.
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    title = (f"TPV31: {station_label}  "
             f"(strike={strike_km:g} km, depth={depth_km:g} km)")
    if t_max is not None:
        title += f"   [0–{t_max:g} s close-up]"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    for ax, (key, ylabel) in zip(axes.flat, PANELS):
        for label, data, color, ls in datasets:
            if data is None:
                continue
            t = data["time_s"]
            y = data[key]
            # Skip channels that are all-NaN (reference μ_eff).
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
        "TPV31: Slip Rate V_strike at All Stations",
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
    directory.
    """
    label, value = parse_labeled_arg(spec)
    if "@" in value:
        directory, prefix = value.rsplit("@", 1)
    else:
        directory, prefix = value, None
    return label, directory, prefix


def detect_mfem_prefix(directory):
    """Auto-detect the station-file prefix used in `directory`.

    Scans for files matching '*_station_faultst*dp*.dat' and extracts
    the common stem (everything before '_station_').  Returns the prefix
    string when exactly one prefix is found; returns None when zero
    matches OR when multiple distinct prefixes coexist (the caller
    should fall back to --mfem-prefix in that case).
    """
    import glob
    if not os.path.isdir(directory):
        return None
    files = glob.glob(os.path.join(directory, "*_station_faultst*dp*.dat"))
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
        description="Visualize TPV31 output: MFEM vs SCEC references"
    )
    parser.add_argument(
        "mfem_dir_positional",
        nargs="?",
        default=None,
        help="(Legacy) MFEM results directory containing "
             "<prefix>_station_faultst*dp*.dat",
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
        default="tpv31",
        help="MFEM file prefix matching --output-prefix passed to "
             "seas_spatial_dyn_driver (default: 'tpv31').  Used as the "
             "fallback when the per-source @PREFIX is absent and "
             "auto-detection finds zero or multiple candidate prefixes.",
    )

    parser.add_argument(
        "--eqdyna", action="store_true",
        help="Include the SCEC EQdyna reference (Benchun Duan).",
    )
    parser.add_argument(
        "--seisol", action="store_true",
        help="Include the SCEC SeisSol ADER-DG reference (Thomas Ulrich).",
    )
    parser.add_argument(
        "--both", action="store_true",
        help="Shortcut: enable both --eqdyna and --seisol.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="benchmark_data",
        help="Directory containing scec_eqdyna/ and scec_seisol/ "
             "subfolders.  Resolved relative to this script if not "
             "absolute.",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip ALL reference traces (EQdyna, SeisSol) even if "
             "--eqdyna / --seisol / --both is given.",
    )

    parser.add_argument(
        "--stations", nargs="+", type=int, default=None,
        help="Specific station indices (1-30). Default: all",
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

    if args.both:
        args.eqdyna = True
        args.seisol = True

    # Default to including BOTH references when no benchmark flag is
    # given and --no-benchmark wasn't requested.  If the user opted in
    # to ANY single reference, the default-on rule does NOT fire — they
    # get only what they asked for.
    if (
        not args.no_benchmark
        and not args.eqdyna
        and not args.seisol
    ):
        args.eqdyna = True
        args.seisol = True

    # Build ordered sources from command-line argv so legend colors
    # match the order the user typed.
    ordered_sources = []
    if args.mfem_dir_positional:
        ordered_sources.append(("mfem", args.mfem_dir_positional))

    mfem_iter = iter(args.mfem or [])
    for arg in sys.argv[1:]:
        if arg == "--both":
            ordered_sources.append(("eqdyna", None))
            ordered_sources.append(("seisol", None))
        elif arg == "--eqdyna":
            ordered_sources.append(("eqdyna", None))
        elif arg == "--seisol":
            ordered_sources.append(("seisol", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))

    # If a benchmark was implicit (default), append it at the end so it
    # still appears in the typed-source list when programmatic callers
    # never put the flag in sys.argv.
    if (
        args.eqdyna
        and not args.no_benchmark
        and not any(s == "eqdyna" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("eqdyna", None))
    if (
        args.seisol
        and not args.no_benchmark
        and not any(s == "seisol" for s, _ in ordered_sources)
    ):
        ordered_sources.append(("seisol", None))

    if not ordered_sources:
        parser.error(
            "No data sources specified.  Provide --mfem DIR and/or "
            "--eqdyna / --seisol / --both."
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
    eqdyna_dir = os.path.join(data_dir, "scec_eqdyna")
    seisol_dir = os.path.join(data_dir, "scec_seisol")

    # Select stations.
    if args.stations:
        stations = [
            SCEC_STATIONS[i - 1] for i in args.stations
            if 1 <= i <= len(SCEC_STATIONS)
        ]
    else:
        stations = SCEC_STATIONS

    # Build the typed-source list with colors and line styles.
    benchmark_types = {"eqdyna", "seisol"}
    benchmark_labels = {
        "eqdyna": "SCEC EQdyna",
        "seisol": "SCEC SeisSol",
    }
    benchmark_dirs = {
        "eqdyna": eqdyna_dir,
        "seisol": seisol_dir,
    }
    benchmark_codes = {"eqdyna": "eqdyna", "seisol": "seisol"}
    benchmark_linestyle = {"eqdyna": "--", "seisol": ":"}
    sources = []
    ci = 0
    for stype, spec in ordered_sources:
        if args.no_benchmark and stype in benchmark_types:
            continue
        color = COLORS[ci % len(COLORS)]
        ci += 1
        if stype in benchmark_types:
            sources.append(
                (benchmark_labels[stype], stype,
                 (benchmark_dirs[stype], benchmark_codes[stype]),
                 color, benchmark_linestyle[stype])
            )
        elif stype == "mfem":
            label, directory, prefix = parse_mfem_spec(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            # Prefix-resolution priority:
            #   1. explicit '@PREFIX' in the spec
            #   2. auto-detected unique prefix in
            #      <directory>/*_station_faultst*dp*.dat
            #   3. --mfem-prefix CLI flag (default 'tpv31')
            if prefix is None:
                prefix = detect_mfem_prefix(directory)
            if prefix is None:
                prefix = args.mfem_prefix
            sources.append(
                (f"MFEM {label}", "mfem",
                 (directory, prefix), color, "-")
            )

    print("=" * 60)
    print("TPV31 Visualization")
    print("=" * 60)
    style_name = {"-": "solid", "--": "dashed", ":": "dotted"}
    for label, stype, info, _color, ls in sources:
        style = style_name.get(ls, ls)
        if stype in benchmark_types:
            bench_dir, code = info
            print(f"  [{style}] {label}: "
                  f"{bench_dir}/tpv31_{code}_x2_*_x3_*.txt")
        else:
            directory, prefix = info
            print(f"  [{style}] {label}: "
                  f"{directory}/{prefix}_station_faultst*dp*.dat")
    print(f"  Stations: {len(stations)}")
    print()

    all_results = []
    for mfem_name, ref_label, strike_km, depth_km in stations:
        datasets = []
        for label, stype, info, color, ls in sources:
            data = None
            if stype in benchmark_types:
                bench_dir, code = info
                path = reference_filename(bench_dir, code, ref_label)
                if os.path.exists(path):
                    data = load_reference_file(path)
            elif stype == "mfem":
                directory, prefix = info
                path = mfem_filename(directory, prefix, mfem_name)
                if os.path.exists(path):
                    data = load_mfem_file(path)
            datasets.append((label, data, color, ls))

        if all(d is None for _, d, _, _ in datasets):
            print(f"  {mfem_name}: no data found")
            continue

        pts_info = []
        for label, data, _c, _ls in datasets:
            if data is not None:
                t_last = data["time_s"][-1]
                pts_info.append(
                    f"{label}: {len(data['time_s'])} pts ({t_last:.2f} s)"
                )
        print(f"  {mfem_name} (strike={strike_km:g}, depth={depth_km:g}): "
              f"{', '.join(pts_info)}")

        all_results.append({
            "station_label": mfem_name,
            "strike_km": strike_km,
            "depth_km": depth_km,
            "datasets": datasets,
        })

        if args.save:
            os.makedirs(args.output_dir, exist_ok=True)
            fname = os.path.join(
                args.output_dir, f"tpv31_{mfem_name}.png"
            )
            plot_station(
                datasets, mfem_name, strike_km, depth_km, save_path=fname
            )
            if args.closeup_t is not None:
                fname_close = os.path.join(
                    args.output_dir,
                    f"tpv31_{mfem_name}_closeup.png"
                )
                plot_station(
                    datasets, mfem_name, strike_km, depth_km,
                    save_path=fname_close, t_max=args.closeup_t,
                )
        else:
            plot_station(datasets, mfem_name, strike_km, depth_km)
            if args.closeup_t is not None:
                plot_station(
                    datasets, mfem_name, strike_km, depth_km,
                    t_max=args.closeup_t,
                )

    if len(all_results) > 1:
        if args.save:
            fname = os.path.join(args.output_dir, "tpv31_overview.png")
            plot_overview(all_results, save_path=fname)
        else:
            plot_overview(all_results)

    print(f"\nPlotted {len(all_results)} stations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
