#!/usr/bin/env python3
"""
TPV6 / TPV7 (SCEC) Benchmark Visualization Script

Plots MFEM dynamic-rupture PER-SIDE on-fault station output alongside the SCEC
DRDG3D community reference (Wenqiang Zhang, 2023) for the bi-material TPV6/TPV7
benchmarks.  Per the TPV6/7 spec, the on-fault split-node stations report, on
EACH side of the bi-material fault, the absolute particle DISPLACEMENT and
VELOCITY (NOT slip / slip-rate), plus the fault traction.

This script mirrors the CLI / structure of ``tpv31/visualize_results.py`` (so the
same invocation works across benchmarks): ``--mfem DIR`` for the run directory,
``--drdg3d`` for the reference, ``--save`` + ``--output-dir`` for PNGs, and the
``--tol-peak`` / ``--tol-rms`` quantitative gate.

near vs far material (CONFIRMED — spec p.3 + the reference data):
    nearside = STRONG/FAST (vp=6000)  ;  farside = WEAK/SLOW (vp=3750 TPV6 / 5000 TPV7)
    -> the WEAK (far) side moves ~3-4x faster (the bi-material signature).

MFEM station format (dynamic/tpv6_stations.hpp::TPV6StationWriter):
    <prefix>_{nearside,farside}_<x2_*_x3_*>.dat
DRDG3D reference format (benchmark_data/scec_drdg3d/):
    {tpv6,tpv7}_drdg3d_{nearside,farside}_<x2_*_x3_*>.txt
Both share columns (MKS-on-fault; stresses in **MPa**):
    t  h-disp  h-vel  h-stress  v-disp  v-vel  v-stress  n-disp  n-vel  n-stress
    h = along-strike, v = along-dip (down-dip), n = fault-normal.

Convention reconciliation for overlay:
  * MFEM and DRDG3D share units (m, m/s, MPa) — NO scale factor.
  * MFEM n-stress is COMPRESSION-POSITIVE; the DRDG3D reference is
    COMPRESSION-NEGATIVE, so the reference n-stress is NEGATED to overlay.

Usage:
    # Single MFEM run vs the DRDG3D reference (default), save PNGs:
    python visualize_results.py --mfem /path/to/results_dir --save \\
        --output-dir /path/to/plots

    # Compare TWO MFEM runs (each prefix auto-detected; directory basename
    # becomes the legend label):
    python visualize_results.py --mfem /path/run_a --mfem /path/run_b --save

    # Per-source label + forced filename prefix, TPV7 problem:
    python visualize_results.py --mfem "p1:/path/run@tpv7" --problem tpv7 --save

    # Quantitative pass/fail gate (headless, no plot) vs DRDG3D:
    python visualize_results.py --mfem /path/run --tol-peak 0.10 --tol-rms 0.10

    # Legacy mode (positional argument == MFEM results dir)
    python visualize_results.py /path/to/results_dir --save
"""

import argparse
import glob
import os
import sys

import numpy as np


def _build_stations():
    """Return the 5 SCEC TPV6/7 on-fault stations (DRDG3D reference set).

    Each entry is ``(name, strike_km, depth_km)`` where ``name`` is the
    filename id shared by the MFEM ``_<side>_<name>.dat`` and the reference
    ``{problem}_drdg3d_<side>_<name>.txt`` files, e.g. ``x2_0_x3_0``.

    Grid (strike, depth) km = (0,0), (-12,0), (12,0), (-12,7.5), (12,7.5).
    """
    return [
        ("x2_0_x3_0",      0.0,   0.0),
        ("x2_-12_x3_0",   -12.0,  0.0),
        ("x2_12_x3_0",     12.0,  0.0),
        ("x2_-12_x3_7.5", -12.0,  7.5),
        ("x2_12_x3_7.5",   12.0,  7.5),
    ]


DRDG3D_STATIONS = _build_stations()

#: Per-side keys + human labels.  nearside = STRONG (larger Zp), far = WEAK.
#: nearside and farside are plotted in SEPARATE figures (never overlapped) — one
#: figure per (station, side) = 5 x 2 = 10 figures.  Within each figure, color +
#: linestyle distinguish the SOURCE (MFEM solid, DRDG3D dashed).
SIDES = ("nearside", "farside")
SIDE_LABEL = {"nearside": "near (STRONG)", "farside": "far (WEAK)"}


def _parse_numeric_table(filepath, min_cols):
    """Load a whitespace-delimited numeric table, skipping comments and the
    ``t h-disp ...`` field-name header line that begins with a letter.
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
                # The DRDG3D files name their columns on a header line that
                # begins with `t h-disp ...` — skip it.
                continue
    if not rows:
        return None
    return np.array(rows)


def _load_station_file(filepath, flip_nstress):
    """Load a TPV6/7 per-side station file (MFEM or DRDG3D — same 10 columns).

    Columns (MKS-on-fault; stresses in MPa):
      0:t 1:h-disp 2:h-vel 3:h-stress 4:v-disp 5:v-vel 6:v-stress
      7:n-disp 8:n-vel 9:n-stress
    ``flip_nstress`` negates n-stress (DRDG3D is compression-NEGATIVE; MFEM is
    compression-POSITIVE) so the two overlay.
    """
    arr = _parse_numeric_table(filepath, min_cols=10)
    if arr is None:
        return None
    sigma_n = -arr[:, 9] if flip_nstress else arr[:, 9]
    return {
        "time_s":   arr[:, 0],
        "h_disp":   arr[:, 1],
        "h_vel":    arr[:, 2],
        "h_stress": arr[:, 3],
        "v_disp":   arr[:, 4],
        "v_vel":    arr[:, 5],
        "v_stress": arr[:, 6],
        "n_disp":   arr[:, 7],
        "n_vel":    arr[:, 8],
        "sigma_n":  sigma_n,
    }


def load_mfem_file(filepath):
    """MFEM per-side station file (n-stress already compression-positive)."""
    return _load_station_file(filepath, flip_nstress=False)


def load_reference_file(filepath):
    """DRDG3D per-side reference file (flip compression-negative n-stress)."""
    return _load_station_file(filepath, flip_nstress=True)


def mfem_filename(results_dir, prefix, side, station):
    return os.path.join(results_dir, f"{prefix}_{side}_{station}.dat")


def reference_filename(bench_dir, problem, side, station):
    return os.path.join(bench_dir, f"{problem}_drdg3d_{side}_{station}.txt")


# ---------------------------------------------------------------------------
# Quantitative pass/fail gate (mirrors tpv31): MFEM vs DRDG3D reference,
# per-side, per-channel.  Invoked by `--tol-peak` / `--tol-rms`; runs headless.
# ---------------------------------------------------------------------------

#: Minimum fraction of the reference time span the MFEM/reference overlap must
#: cover — a wall-truncated run must not masquerade as agreement.
MIN_COVERAGE_FRAC = 0.90

#: Small absolute denominator floor for the peak-normalized metric so
#: pre-rupture zero-amplitude windows do not blow the relative error up.
METRIC_FLOOR = 1.0e-6

#: Floor on peak |reference dip motion| below which the `v_*` (dip) channels
#: are skipped as trivial (TPV6/7 motion is strike-dominated).
DIP_MOTION_FLOOR = 1.0e-3

#: Channels ALWAYS gated (strike displacement/velocity/stress + normal stress).
STRIKE_GATE_CHANNELS = ("h_vel", "h_disp", "h_stress", "sigma_n")

#: Dip channels gated ONLY when the reference shows non-trivial dip motion.
DIP_GATE_CHANNELS = ("v_vel", "v_disp")


def _finite_pair(t, y):
    """Return (t, y) restricted to samples where BOTH are finite."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(t) & np.isfinite(y)
    return t[mask], y[mask]


def gate_channel(mfem_t, mfem_y, ref_t, ref_y, floor=METRIC_FLOOR,
                 min_coverage_frac=MIN_COVERAGE_FRAC):
    """Peak/RMS-normalized gate metric for ONE channel over the [t0,t1] overlap.

    Returns ``(peak_rel, rms_rel, covered)``; ``covered`` is True iff the
    overlap span is >= ``min_coverage_frac`` of the reference span.
    """
    mfem_t, mfem_y = _finite_pair(mfem_t, mfem_y)
    ref_t, ref_y = _finite_pair(ref_t, ref_y)
    if mfem_t.size == 0 or ref_t.size == 0:
        return float("nan"), float("nan"), False

    ref_span = ref_t[-1] - ref_t[0]
    t0 = max(mfem_t[0], ref_t[0])
    t1 = min(mfem_t[-1], ref_t[-1])
    overlap = t1 - t0
    if ref_span > 0.0:
        covered = bool(overlap >= (min_coverage_frac * ref_span))
    else:
        covered = bool(overlap >= 0.0)

    in_window = (ref_t >= t0) & (ref_t <= t1)
    ref_t_w = ref_t[in_window]
    ref_y_w = ref_y[in_window]
    if ref_t_w.size == 0:
        return float("nan"), float("nan"), covered

    order = np.argsort(mfem_t)
    mfem_interp = np.interp(ref_t_w, mfem_t[order], mfem_y[order])
    diff = mfem_interp - ref_y_w
    peak_abs = float(np.max(np.abs(diff)))
    rms_abs = float(np.sqrt(np.mean(diff * diff)))
    peak_ref = max(float(np.max(np.abs(ref_y_w))), floor)
    rms_ref = max(float(np.sqrt(np.mean(ref_y_w * ref_y_w))), floor)
    return peak_abs / peak_ref, rms_abs / rms_ref, covered


def gated_channels_for_station(ref_data, dip_motion_floor=DIP_MOTION_FLOOR):
    """Channel keys to gate: strike+normal always; dip only if non-trivial."""
    channels = list(STRIKE_GATE_CHANNELS)
    v_vel = np.asarray(ref_data.get("v_vel"), dtype=float)
    v_vel = v_vel[np.isfinite(v_vel)]
    if v_vel.size and float(np.max(np.abs(v_vel))) > dip_motion_floor:
        channels.extend(DIP_GATE_CHANNELS)
    return channels


def _is_all_nan_or_empty(y):
    y = np.asarray(y, dtype=float)
    return y.size == 0 or np.all(np.isnan(y))


def run_tolerance_gate(station_pairs, tol_peak, tol_rms,
                       floor=METRIC_FLOOR,
                       min_coverage_frac=MIN_COVERAGE_FRAC,
                       dip_motion_floor=DIP_MOTION_FLOOR,
                       printer=print):
    """Run the pass/fail gate over a list of ``(label, mfem_data, ref_data)``.

    Returns the process exit code: ``0`` if every gated channel of every
    (station, side) passes BOTH active tolerances and clears the coverage
    guard; ``1`` otherwise.
    """
    peak_cap = float("inf") if tol_peak is None else tol_peak
    rms_cap = float("inf") if tol_rms is None else tol_rms

    any_fail = False
    n_pairs_checked = 0
    n_channels_checked = 0

    header = (f"{'station:side':<26} {'channel':<10} "
              f"{'peak_rel':>10} {'rms_rel':>10}  result")
    printer("=" * 62)
    printer("TPV6/7 quantitative gate (MFEM vs DRDG3D reference)")
    tol_desc = []
    if tol_peak is not None:
        tol_desc.append(f"--tol-peak={tol_peak:g}")
    if tol_rms is not None:
        tol_desc.append(f"--tol-rms={tol_rms:g}")
    printer(f"  tolerances: {', '.join(tol_desc) if tol_desc else '(none)'}")
    printer(f"  min coverage fraction: {min_coverage_frac:g}")
    printer("=" * 62)
    printer(header)
    printer("-" * len(header))

    for label, mfem_data, ref_data in station_pairs:
        if mfem_data is None or ref_data is None:
            printer(f"{label:<26} {'(missing)':<10} "
                    f"{'--':>10} {'--':>10}  SKIP (no data)")
            continue

        n_pairs_checked += 1
        channels = gated_channels_for_station(
            ref_data, dip_motion_floor=dip_motion_floor)
        pair_fail = False
        ref_t = ref_data["time_s"]
        mfem_t = mfem_data["time_s"]

        for key in channels:
            ref_y = ref_data.get(key)
            mfem_y = mfem_data.get(key)
            if ref_y is None or _is_all_nan_or_empty(ref_y):
                printer(f"{label:<26} {key:<10} "
                        f"{'--':>10} {'--':>10}  SKIP (ref all-NaN)")
                continue
            peak_rel, rms_rel, covered = gate_channel(
                mfem_t, mfem_y, ref_t, ref_y,
                floor=floor, min_coverage_frac=min_coverage_frac)
            n_channels_checked += 1
            if not covered:
                pair_fail = True
                printer(f"{label:<26} {key:<10} "
                        f"{'--':>10} {'--':>10}  FAIL (insufficient coverage)")
                continue
            channel_fail = (peak_rel > peak_cap) or (rms_rel > rms_cap)
            if channel_fail:
                pair_fail = True
            printer(f"{label:<26} {key:<10} "
                    f"{peak_rel:>10.4f} {rms_rel:>10.4f}  "
                    f"{'PASS' if not channel_fail else 'FAIL'}")

        if pair_fail:
            any_fail = True

    printer("-" * len(header))
    if n_pairs_checked == 0:
        printer("SUMMARY: GATE FAILED — no station/side had both MFEM and "
                "reference data to compare.")
        return 1
    if n_channels_checked == 0:
        printer("SUMMARY: GATE FAILED — stations were found but no gatable "
                "channel had comparable reference data.")
        return 1
    if any_fail:
        printer(f"SUMMARY: GATE FAILED — {n_pairs_checked} station/side(s), "
                f"{n_channels_checked} channel(s) checked; at least one "
                f"exceeded the tolerance band or failed coverage.")
        return 1
    printer(f"SUMMARY: GATE PASSED — {n_pairs_checked} station/side(s), "
            f"{n_channels_checked} channel(s) within band "
            f"(peak<= {peak_cap:g}, rms<= {rms_cap:g}).")
    return 0


PANELS = [
    ("h_vel",    "h-vel (strike) [m/s]"),
    ("h_disp",   "h-disp (strike) [m]"),
    ("h_stress", "h-stress (strike) [MPa]"),
    ("v_vel",    "v-vel (dip) [m/s]"),
    ("v_disp",   "v-disp (dip) [m]"),
    ("sigma_n",  "n-stress [MPa, compression +]"),
]


def plot_station(datasets, station_name, side_label, strike_km, depth_km,
                 problem, save_path=None, t_max=None):
    """6-panel single-side station comparison (MFEM vs DRDG3D).

    One figure per (station, side): nearside and farside are NEVER overlapped.
    Within the figure, color + linestyle distinguish the source (MFEM solid,
    DRDG3D dashed).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    title = (f"{problem.upper()} on-fault {station_name} — {side_label}  "
             f"(strike={strike_km:g} km, depth={depth_km:g} km)")
    if t_max is not None:
        title += f"   [0-{t_max:g} s close-up]"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    for ax, (key, ylabel) in zip(axes.flat, PANELS):
        for label, data, color, ls in datasets:
            if data is None:
                continue
            t = data["time_s"]
            y = data[key]
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
    """Parse '[LABEL:]DIR[@PREFIX]' -> (label_or_None, dir, prefix_or_None)."""
    label, value = parse_labeled_arg(spec)
    if "@" in value:
        directory, prefix = value.rsplit("@", 1)
    else:
        directory, prefix = value, None
    return label, directory, prefix


def detect_mfem_prefix(directory):
    """Auto-detect the per-side station-file prefix in `directory`.

    Scans for '*_{nearside,farside}_x2_*.dat' and extracts the common stem
    (everything before '_nearside_' / '_farside_').  Returns the prefix when
    exactly one is found; None on zero or multiple (caller falls back to
    --mfem-prefix).
    """
    if not os.path.isdir(directory):
        return None
    prefixes = set()
    for side in SIDES:
        for f in glob.glob(os.path.join(directory, f"*_{side}_x2_*.dat")):
            name = os.path.basename(f)
            token = f"_{side}_"
            if token in name:
                prefixes.add(name.split(token, 1)[0])
    if len(prefixes) == 1:
        return next(iter(prefixes))
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Visualize TPV6/7 output: MFEM per-side vs DRDG3D reference"
    )
    parser.add_argument(
        "mfem_dir_positional", nargs="?", default=None,
        help="(Legacy) MFEM results directory containing "
             "<prefix>_{nearside,farside}_x2_*.dat",
    )
    parser.add_argument(
        "--mfem", action="append", metavar="[LABEL:]DIR[@PREFIX]",
        help="MFEM results directory.  Use 'label:dir' for a custom legend "
             "label, '@prefix' to force a filename prefix.  Repeatable.",
    )
    parser.add_argument(
        "--mfem-prefix", default=None,
        help="MFEM file prefix (the run's --output-prefix).  Default: the "
             "--problem value ('tpv6'/'tpv7').  Used as the fallback when the "
             "per-source @PREFIX is absent and auto-detection is ambiguous.",
    )
    parser.add_argument(
        "--problem", default="tpv6", choices=["tpv6", "tpv7"],
        help="Benchmark problem (selects the DRDG3D reference filenames and "
             "the default MFEM prefix).  Default: tpv6.",
    )
    parser.add_argument(
        "--drdg3d", action="store_true",
        help="Include the SCEC DRDG3D reference (Wenqiang Zhang).  Default-on "
             "when no reference flag is given (unless --no-benchmark).",
    )
    parser.add_argument(
        "--benchmark-dir", default="benchmark_data",
        help="Directory containing scec_drdg3d/.  Resolved relative to this "
             "script if not absolute.",
    )
    parser.add_argument(
        "--no-benchmark", action="store_true",
        help="Skip the DRDG3D reference (plot the MFEM run only).",
    )
    parser.add_argument(
        "--stations", nargs="+", type=int, default=None,
        help="Specific station indices (1-5).  Default: all.",
    )
    parser.add_argument(
        "--save", action="store_true",
        help="Save plots as PNG (default: display interactively).",
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory for output plots.",
    )
    parser.add_argument(
        "--closeup-t", type=float, default=None,
        help="Additionally produce a close-up plot truncated to this time [s].",
    )
    parser.add_argument(
        "--tol-peak", type=float, default=None,
        help="Run the quantitative pass/fail gate instead of plotting: fail "
             "(exit non-zero) if ANY gated channel of ANY station/side has "
             "peak-normalized error above this fraction (e.g. 0.10 = 10%%).  "
             "May be combined with --tol-rms.  Runs headless (no matplotlib).",
    )
    parser.add_argument(
        "--tol-rms", type=float, default=None,
        help="Run the quantitative gate (see --tol-peak): fail if any gated "
             "channel's RMS-normalized error exceeds this fraction.",
    )

    args = parser.parse_args()
    mfem_prefix_default = args.mfem_prefix or args.problem

    # Default to including the DRDG3D reference when no reference flag was
    # given and --no-benchmark wasn't requested.
    if not args.no_benchmark and not args.drdg3d:
        args.drdg3d = True

    # Build ordered sources from argv so legend colors match the typed order.
    ordered_sources = []
    if args.mfem_dir_positional:
        ordered_sources.append(("mfem", args.mfem_dir_positional))
    mfem_iter = iter(args.mfem or [])
    for arg in sys.argv[1:]:
        if arg == "--drdg3d":
            ordered_sources.append(("drdg3d", None))
        elif arg == "--mfem":
            ordered_sources.append(("mfem", next(mfem_iter)))
    if (args.drdg3d and not args.no_benchmark
            and not any(s == "drdg3d" for s, _ in ordered_sources)):
        ordered_sources.append(("drdg3d", None))

    if not ordered_sources:
        parser.error("No data sources specified.  Provide --mfem DIR "
                     "and/or --drdg3d.")

    gate_mode = (args.tol_peak is not None) or (args.tol_rms is not None)
    if not gate_mode:
        try:
            import matplotlib
            if args.save:
                matplotlib.use("Agg")
        except ImportError:
            print("Error: matplotlib required. "
                  "Install with: pip install matplotlib")
            return 1

    # Resolve benchmark directory (relative to this script unless absolute).
    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), args.benchmark_dir)
    if not os.path.isabs(args.benchmark_dir) and not os.path.isdir(data_dir):
        data_dir = args.benchmark_dir
    drdg3d_dir = os.path.join(data_dir, "scec_drdg3d")

    # Select stations (1-based indices).
    if args.stations:
        stations = [DRDG3D_STATIONS[i - 1] for i in args.stations
                    if 1 <= i <= len(DRDG3D_STATIONS)]
    else:
        stations = DRDG3D_STATIONS

    # ---- Quantitative gate ------------------------------------------------
    if gate_mode:
        mfem_sources = []
        for stype, spec in ordered_sources:
            if stype != "mfem":
                continue
            label, directory, prefix = parse_mfem_spec(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            if prefix is None:
                prefix = detect_mfem_prefix(directory)
            if prefix is None:
                prefix = mfem_prefix_default
            mfem_sources.append((label, directory, prefix))

        if not mfem_sources:
            print("Error: --tol-peak/--tol-rms requires at least one --mfem "
                  "results directory to gate against the DRDG3D reference.")
            return 1

        gate_pairs = []
        for label, directory, prefix in mfem_sources:
            for station, _s_km, _d_km in stations:
                for side in SIDES:
                    mfem_path = mfem_filename(directory, prefix, side, station)
                    ref_path = reference_filename(
                        drdg3d_dir, args.problem, side, station)
                    mfem_data = (load_mfem_file(mfem_path)
                                 if os.path.exists(mfem_path) else None)
                    ref_data = (load_reference_file(ref_path)
                                if os.path.exists(ref_path) else None)
                    if mfem_data is None and ref_data is None:
                        continue
                    base = f"{station}:{side[:4]}"
                    pair_label = (base if len(mfem_sources) == 1
                                  else f"{label}:{base}")
                    gate_pairs.append((pair_label, mfem_data, ref_data))

        return run_tolerance_gate(
            gate_pairs, tol_peak=args.tol_peak, tol_rms=args.tol_rms)

    # ---- Plot mode --------------------------------------------------------
    # One color per typed source (each --mfem run + the reference); near/far
    # are distinguished by linestyle (solid/dashed).
    sources = []
    ci = 0
    for stype, spec in ordered_sources:
        if args.no_benchmark and stype == "drdg3d":
            continue
        color = COLORS[ci % len(COLORS)]
        ci += 1
        if stype == "drdg3d":
            sources.append(("DRDG3D", "drdg3d", (drdg3d_dir, args.problem),
                            color))
        elif stype == "mfem":
            label, directory, prefix = parse_mfem_spec(spec)
            if label is None:
                label = os.path.basename(os.path.normpath(directory))
            if prefix is None:
                prefix = detect_mfem_prefix(directory)
            if prefix is None:
                prefix = mfem_prefix_default
            sources.append((f"MFEM {label}", "mfem", (directory, prefix),
                            color))

    print("=" * 62)
    print(f"{args.problem.upper()} Visualization (per-side: near=solid, "
          f"far=dashed)")
    print("=" * 62)
    for label, stype, info, _color in sources:
        if stype == "drdg3d":
            bench_dir, problem = info
            print(f"  {label}: {bench_dir}/{problem}_drdg3d_<side>_x2_*.txt")
        else:
            directory, prefix = info
            print(f"  {label}: {directory}/{prefix}_<side>_x2_*.dat")
    print(f"  Stations: {len(stations)}")
    print()

    # One figure PER (station, side) — nearside and farside are NEVER
    # overlapped.  5 stations x 2 sides = 10 figures.  Each figure overlays
    # only the SOURCES (MFEM vs DRDG3D) for that single side.
    n_figs = 0
    for station, strike_km, depth_km in stations:
        for side in SIDES:
            datasets = []
            for label, stype, info, color in sources:
                data = None
                if stype == "drdg3d":
                    bench_dir, problem = info
                    path = reference_filename(bench_dir, problem, side, station)
                    if os.path.exists(path):
                        data = load_reference_file(path)
                else:
                    directory, prefix = info
                    path = mfem_filename(directory, prefix, side, station)
                    if os.path.exists(path):
                        data = load_mfem_file(path)
                ls = "-" if stype == "mfem" else "--"
                datasets.append((label, data, color, ls))

            if all(d is None for _, d, _, _ in datasets):
                print(f"  {station} {SIDE_LABEL[side]}: no data found")
                continue

            pts = [f"{lbl}: {len(d['time_s'])} pts ({d['time_s'][-1]:.2f} s)"
                   for lbl, d, _c, _ls in datasets if d is not None]
            print(f"  {station} {SIDE_LABEL[side]} "
                  f"(strike={strike_km:g}, depth={depth_km:g}): "
                  f"{', '.join(pts)}")

            if args.save:
                os.makedirs(args.output_dir, exist_ok=True)
                fname = os.path.join(
                    args.output_dir, f"{args.problem}_{station}_{side}.png")
                plot_station(datasets, station, SIDE_LABEL[side], strike_km,
                             depth_km, args.problem, save_path=fname)
                if args.closeup_t is not None:
                    fname_c = os.path.join(
                        args.output_dir,
                        f"{args.problem}_{station}_{side}_closeup.png")
                    plot_station(datasets, station, SIDE_LABEL[side], strike_km,
                                 depth_km, args.problem, save_path=fname_c,
                                 t_max=args.closeup_t)
            else:
                plot_station(datasets, station, SIDE_LABEL[side], strike_km,
                             depth_km, args.problem)
                if args.closeup_t is not None:
                    plot_station(datasets, station, SIDE_LABEL[side], strike_km,
                                 depth_km, args.problem, t_max=args.closeup_t)
            n_figs += 1

    print(f"\nPlotted {n_figs} figures "
          f"({len(stations)} station(s) x {len(SIDES)} sides).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
