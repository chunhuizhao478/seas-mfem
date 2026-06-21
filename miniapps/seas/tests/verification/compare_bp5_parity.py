#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# compare_bp5_parity.py — Phase 7 of
# document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
#
# Compares BP5 SCEC station traces between two runs (column-by-column) and
# reports PASS/FAIL against a tolerance.  Two intended uses:
#   1. spatial_seas (MUMPS) vs seas_driver golden  -> golden tolerance.
#   2. spatial_seas (cg_amg) vs spatial_seas (MUMPS) -> rel < 1e-6.
#
# Station files are SCEC BP5 format (8 columns, one header line):
#   time(s) slip_strike(m) slip_dip(m) log10(V_strike)(m/s) log10(V_dip)(m/s)
#   tau_strike(MPa) tau_dip(MPa) log10(state)(s)
#
# The RUNS that produce these files are Frontera (the SAFS coupling needs the
# BP5 mesh — project memory: no local production-mesh runs).  This script is the
# comparison harness; point it at the two output directories afterwards:
#
#   python3 compare_bp5_parity.py --a out_mumps --b out_golden \
#       --prefix-a spatial_seas --prefix-b bp5 --rtol 1e-6 --atol 1e-9
#
# Standard library only (no numpy) — matches scripts/estimate_output_size.py.

import argparse
import glob
import math
import os
import sys

COLS = ["time", "slip_strike", "slip_dip", "log10V_strike",
        "log10V_dip", "tau_strike", "tau_dip", "log10state"]
NCOL = len(COLS)
# R-704: log10 columns (a log diff is already a ratio) are compared with an
# ABSOLUTE tolerance (--log-atol), not a relative one.
LOG_COLS = {3, 4, 7}   # log10V_strike, log10V_dip, log10state


def parse_station(path):
    """Return a list of rows (each a list of NCOL floats); skip non-numeric
    (header/comment) lines.  Raises ValueError on a malformed numeric row."""
    rows = []
    started = False
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            toks = s.split()
            try:
                vals = [float(t) for t in toks]
            except ValueError:
                # R-706: a non-numeric line is the (single) header ONLY before
                # numeric data begins; one appearing later is a corrupt row.
                if started:
                    raise ValueError(
                        "{}:{}: non-numeric row after data began: {!r}".format(
                            path, ln, s))
                continue
            if len(vals) != NCOL:
                raise ValueError(
                    "{}:{}: expected {} columns, got {}".format(
                        path, ln, NCOL, len(vals)))
            rows.append(vals)
            started = True
    return rows


def station_name(path, prefix):
    """{prefix}_{name}.txt -> name."""
    base = os.path.basename(path)
    base = base[:-4] if base.endswith(".txt") else base
    pre = prefix + "_"
    return base[len(pre):] if base.startswith(pre) else base


def compare_pair(rows_a, rows_b, rtol, atol, log_atol=1e-6):
    """Per-column (max_abs, max_rel) over the common rows + an aligned flag.

    Fails (ok=False) on: a row-count mismatch; a non-finite value (NaN/Inf) in
    EITHER trace (R-701 — a gate must not silently pass a blown-up run); a
    time-column (col 0) misalignment (R-705); a log10 column exceeding log_atol
    absolutely (R-704); any other column exceeding BOTH atol and rtol."""
    na, nb = len(rows_a), len(rows_b)
    n = min(na, nb)
    aligned = (na == nb)
    max_abs = [0.0] * NCOL
    max_rel = [0.0] * NCOL
    ok = aligned
    for i in range(n):
        # R-705: rows i must be at the same time (col 0) to be comparable.
        ta, tb = rows_a[i][0], rows_b[i][0]
        if not (math.isfinite(ta) and math.isfinite(tb)) \
           or abs(ta - tb) > 1e-6 * (abs(tb) + 1.0):
            ok = False
        for c in range(NCOL):
            a, b = rows_a[i][c], rows_b[i][c]
            # R-701: any non-finite value (NaN/Inf) is a hard failure — Python's
            # `NaN > atol` is False, so the plain comparison below would PASS it.
            if not (math.isfinite(a) and math.isfinite(b)):
                ok = False
                max_abs[c] = float("inf")
                max_rel[c] = float("inf")
                continue
            d = abs(a - b)
            r = d / (abs(b) + atol)
            if d > max_abs[c]:
                max_abs[c] = d
            if r > max_rel[c]:
                max_rel[c] = r
            if c in LOG_COLS:
                if d > log_atol:           # R-704: absolute tol on log columns
                    ok = False
            elif d > atol and r > rtol:
                ok = False
    return ok, max_abs, max_rel, na, nb


def main(argv=None):
    ap = argparse.ArgumentParser(description="BP5 station-trace parity check.")
    ap.add_argument("--a", required=True, help="run A output directory")
    ap.add_argument("--b", required=True, help="run B (reference) output directory")
    ap.add_argument("--prefix-a", default="spatial_seas")
    ap.add_argument("--prefix-b", default="bp5")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--log-atol", type=float, default=1e-6,
                    help="absolute tolerance on log10 columns (R-704)")
    args = ap.parse_args(argv)

    files_a = sorted(glob.glob(os.path.join(args.a, args.prefix_a + "_*.txt")))
    if not files_a:
        print("ERROR: no '{}_*.txt' station files in {}".format(
            args.prefix_a, args.a), file=sys.stderr)
        return 2

    overall_ok = True
    n_compared = 0
    for fa in files_a:
        name = station_name(fa, args.prefix_a)
        fb = os.path.join(args.b, "{}_{}.txt".format(args.prefix_b, name))
        if not os.path.exists(fb):
            print("MISSING in B: {} (station '{}')".format(fb, name))
            overall_ok = False
            continue
        try:
            ra, rb = parse_station(fa), parse_station(fb)
        except (OSError, ValueError) as e:
            print("ERROR reading {}: {}".format(name, e))
            overall_ok = False
            continue
        ok, mabs, mrel, na, nb = compare_pair(
            ra, rb, args.rtol, args.atol, args.log_atol)
        n_compared += 1
        overall_ok = overall_ok and ok
        worst_c = max(range(NCOL), key=lambda c: mrel[c]) if mrel else 0
        tag = "PASS" if ok else "FAIL"
        rowinfo = "" if na == nb else "  [ROW MISMATCH {} vs {}]".format(na, nb)
        print("[{}] {:<16} worst col '{}': abs={:.3e} rel={:.3e}{}".format(
            tag, name, COLS[worst_c], mabs[worst_c], mrel[worst_c], rowinfo))

    print("\n{} station(s) compared; rtol={:.1e} atol={:.1e} -> {}".format(
        n_compared, args.rtol, args.atol, "PASS" if overall_ok else "FAIL"))
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
