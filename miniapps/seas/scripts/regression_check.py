#!/usr/bin/env python3
"""
BP5 Regression Check Script

Compares simulation output against golden reference data for all 7 non-time
SCEC columns per station (slip_strike, slip_dip, log10_V_strike, log10_V_dip,
tau_strike, tau_dip, log10_state) plus global output (log10(Vmax)).

Interpolates to common time points and reports per-station, per-field relative
L2 errors. Returns nonzero exit code on any tolerance exceedance.

Usage:
    python3 regression_check.py \\
        --sim-dir OUTPUT_DIR --sim-prefix bp5_full \\
        --ref-dir bp5/benchmark_data/golden_serial_1r_50step \\
        [--ref-prefix bp5_verify_serial_job7643780] \\
        [--tolerance 1e-12] \\
        [--tfinal 100.0]
"""

import argparse
import os
import sys
import glob
import numpy as np


# SCEC BP5 station file columns (after time):
STATION_FIELDS = [
    "slip_strike",
    "slip_dip",
    "log10_V_strike",
    "log10_V_dip",
    "tau_strike",
    "tau_dip",
    "log10_state",
]

# BP5 default 10 on-fault stations
STATION_NAMES = [
    "fltst_strk+00dp+00",
    "fltst_strk+00dp+10",
    "fltst_strk+00dp+22",
    "fltst_strk+16dp+00",
    "fltst_strk+16dp+10",
    "fltst_strk+36dp+00",
    "fltst_strk-16dp+00",
    "fltst_strk-16dp+10",
    "fltst_strk-24dp+10",
    "fltst_strk-36dp+00",
]


def load_station_file(filepath, t_max=None):
    """Load BP5 SCEC-format station file (8 columns, # comment lines).

    Returns:
        numpy array of shape (N, 8): [time, slip_s, slip_d, logV_s, logV_d,
                                       tau_s, tau_d, log_state]
        or None on failure.
    """
    rows = []
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "slip" in line:  # skip header
                    continue
                vals = line.split()
                if len(vals) < 8:
                    continue
                try:
                    row = [float(v) for v in vals[:8]]
                except ValueError:
                    continue
                if t_max is not None and row[0] > t_max:
                    break
                rows.append(row)
    except (IOError, OSError):
        return None
    if not rows:
        return None
    return np.array(rows)


def load_global_file(filepath, t_max=None):
    """Load BP5 global output file (2 columns: time, log10(Vmax)).

    Returns:
        numpy array of shape (N, 2) or None.
    """
    rows = []
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                vals = line.split()
                if len(vals) < 2:
                    continue
                try:
                    row = [float(v) for v in vals[:2]]
                except ValueError:
                    continue
                if t_max is not None and row[0] > t_max:
                    break
                rows.append(row)
    except (IOError, OSError):
        return None
    if not rows:
        return None
    return np.array(rows)


def interpolate_onto(x_ref, y_ref, x_target):
    """Linear interpolation of (x_ref, y_ref) onto x_target grid."""
    return np.interp(x_target, x_ref, y_ref)


def relative_l2_error(sim, ref):
    """Compute relative L2 error: ||sim - ref||_2 / ||ref||_2."""
    diff = sim - ref
    num = np.sqrt(np.sum(diff ** 2))
    den = np.sqrt(np.sum(ref ** 2))
    if den < 1e-30:
        return 0.0 if num < 1e-30 else 1e30
    return num / den


def count_events(log10_vmax, threshold_log10=-3.0):
    """Count seismic events: upward crossings of threshold in log10(Vmax)."""
    count = 0
    above = False
    for v in log10_vmax:
        if v > threshold_log10 and not above:
            count += 1
            above = True
        elif v <= threshold_log10:
            above = False
    return count


def detect_ref_prefix(ref_dir):
    """Auto-detect reference file prefix from directory contents.

    Looks for *_global.txt or *_fltst_strk+00dp+00.txt.
    """
    for suffix in ["_global.txt", "_fltst_strk+00dp+00.txt"]:
        pattern = os.path.join(ref_dir, f"*{suffix}")
        matches = sorted(glob.glob(pattern))
        if len(matches) == 1:
            fname = os.path.basename(matches[0])
            return fname[: -len(suffix)]
        elif len(matches) > 1:
            fname = os.path.basename(matches[0])
            prefix = fname[: -len(suffix)]
            print(f"WARNING: Multiple prefixes in {ref_dir}, using '{prefix}'")
            return prefix
    return None


def main():
    parser = argparse.ArgumentParser(
        description="BP5 regression check: compare simulation vs golden reference"
    )
    parser.add_argument(
        "--sim-dir", required=True, help="Directory containing simulation output"
    )
    parser.add_argument(
        "--sim-prefix", required=True, help="Simulation output file prefix"
    )
    parser.add_argument(
        "--ref-dir", required=True, help="Directory containing golden reference data"
    )
    parser.add_argument(
        "--ref-prefix",
        default=None,
        help="Reference file prefix (auto-detected if omitted)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="Tolerance for pass/fail (None = informational only)",
    )
    parser.add_argument(
        "--tfinal",
        type=float,
        default=None,
        help="Maximum time to compare (seconds)",
    )
    parser.add_argument(
        "--stations",
        nargs="*",
        default=None,
        help="Station names to check (default: all 10 BP5 stations)",
    )
    args = parser.parse_args()

    # Resolve reference prefix
    ref_prefix = args.ref_prefix
    if ref_prefix is None:
        ref_prefix = detect_ref_prefix(args.ref_dir)
        if ref_prefix is None:
            print(f"ERROR: Could not auto-detect ref prefix in {args.ref_dir}")
            return 2
        print(f"Auto-detected ref prefix: {ref_prefix}")

    station_names = args.stations if args.stations else STATION_NAMES

    print("=" * 80)
    print("BP5 Regression Check")
    print("=" * 80)
    print(f"  Sim dir:    {args.sim_dir}")
    print(f"  Sim prefix: {args.sim_prefix}")
    print(f"  Ref dir:    {args.ref_dir}")
    print(f"  Ref prefix: {ref_prefix}")
    if args.tolerance is not None:
        print(f"  Tolerance:  {args.tolerance:.2e}")
    else:
        print("  Mode:       informational (no tolerance check)")
    print()

    any_fail = False
    max_error = 0.0
    worst_station = ""
    worst_field = ""
    stations_compared = 0

    # --- Station comparisons ---
    for st_name in station_names:
        sim_file = os.path.join(args.sim_dir, f"{args.sim_prefix}_{st_name}.txt")
        ref_file = os.path.join(args.ref_dir, f"{ref_prefix}_{st_name}.txt")

        sim_data = load_station_file(sim_file, args.tfinal)
        ref_data = load_station_file(ref_file, args.tfinal)

        if sim_data is None:
            print(f"  Station {st_name:>24s}: sim MISSING")
            if args.tolerance is not None:
                any_fail = True
            continue
        if ref_data is None:
            print(f"  Station {st_name:>24s}: ref MISSING ({ref_file})")
            continue

        stations_compared += 1
        print(
            f"  Station {st_name:>24s} "
            f"({sim_data.shape[0]} vs {ref_data.shape[0]} pts):"
        )

        for f_idx, f_name in enumerate(STATION_FIELDS):
            col = f_idx + 1  # skip time column
            ref_interp = interpolate_onto(
                ref_data[:, 0], ref_data[:, col], sim_data[:, 0]
            )
            err = relative_l2_error(sim_data[:, col], ref_interp)

            field_fail = args.tolerance is not None and err > args.tolerance
            status = "FAIL" if field_fail else "ok"

            print(f"    {f_name:>16s}: L2_rel = {err:.6e}  [{status}]")

            if field_fail:
                any_fail = True
            if err > max_error:
                max_error = err
                worst_station = st_name
                worst_field = f_name

    # --- Global output comparison ---
    sim_global = os.path.join(args.sim_dir, f"{args.sim_prefix}_global.txt")
    ref_global = os.path.join(args.ref_dir, f"{ref_prefix}_global.txt")

    sim_g = load_global_file(sim_global, args.tfinal)
    ref_g = load_global_file(ref_global, args.tfinal)

    if sim_g is not None and ref_g is not None:
        ref_gv_interp = interpolate_onto(ref_g[:, 0], ref_g[:, 1], sim_g[:, 0])
        err = relative_l2_error(sim_g[:, 1], ref_gv_interp)
        g_fail = args.tolerance is not None and err > args.tolerance
        status = "FAIL" if g_fail else "ok"
        print(f"  Global log10(Vmax): L2_rel = {err:.6e}  [{status}]")
        if g_fail:
            any_fail = True
        if err > max_error:
            max_error = err
            worst_station = "global"
            worst_field = "log10_Vmax"

        # Event count comparison
        sim_events = count_events(sim_g[:, 1])
        ref_events = count_events(ref_g[:, 1])
        if sim_events != ref_events:
            print(f"  Event count: {sim_events} vs {ref_events} [MISMATCH]")
            if args.tolerance is not None:
                any_fail = True
        else:
            print(f"  Event count: {sim_events} [ok]")
    else:
        if sim_g is None:
            print(f"  Global: sim file missing ({sim_global})")
        if ref_g is None:
            print(f"  Global: ref file missing ({ref_global})")

    # --- Summary ---
    print()
    print("-" * 80)
    print(f"  Stations compared: {stations_compared} / {len(station_names)}")
    print(f"  Max relative L2 error: {max_error:.6e}", end="")
    if worst_station:
        print(f" ({worst_station} / {worst_field})", end="")
    print()

    if args.tolerance is not None:
        if stations_compared == 0:
            print(
                "  REGRESSION CHECK: FAIL "
                "(no stations compared — check --ref-dir and --ref-prefix)"
            )
            any_fail = True
        elif any_fail:
            print(f"  REGRESSION CHECK: FAIL (tolerance {args.tolerance:.2e} exceeded)")
        else:
            print("  REGRESSION CHECK: PASS")
    print("=" * 80)

    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
