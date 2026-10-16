#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# compare_runs.py — fully dynamic BP5 cycle (plan rev 3): compare the time
# series of two runs, file by file and row by row.
#
#   python3 tests/cycle/compare_runs.py DIR_A PREFIX_A DIR_B PREFIX_B
#           [--exact] [--tol T] [--l2] [--station-columns LIST]
#           [--checkpoints] [--dedup]
#
# Files are paired by what follows the prefix: A/PREFIX_A_fltst_strk+00dp+10.txt
# with B/PREFIX_B_fltst_strk+00dp+10.txt, A/PREFIX_A_global.txt with
# B/PREFIX_B_global.txt.  Header lines ('#') are skipped, so a seas_bp5_full
# run compares with a seas_fd_cycle_driver run.  For each pair it prints the
# row counts, how many rows are identical as text, and the largest relative
# difference of each column.
#   --exact        pass only if every data row is identical as text
#   --tol T        pass if every column's max relative difference <= T
#                  (default 0: identical values)
#   --l2           with --tol: measure each column by its relative L2
#                  difference ||a - b|| / ||b|| instead (the measure of
#                  scripts/regression_check.py)
#   --station-columns LIST
#                  compare only these columns of the station files
#                  (0-based, comma separated; e.g. 0,1,2,7 = time, slip,
#                  state, leaving out log10 V and tau, which the station
#                  sigma_n fix ef70517 changes); other files: all columns
#   --checkpoints  also compare PREFIX_checkpoint_r*.txt byte for byte
#                  (implies the same commit; used by test_cycle_restart)
#   --dedup        drop a row equal to the row before it, in both files:
#                  seas_bp5_full writes the t = 0 global row twice (once
#                  before the time loop, once from the step-0 monitor) and
#                  the cycle driver once
# Exit code 0 = PASSED.

import argparse
import glob
import math
import os
import sys


def rows(path, dedup=False):
    out = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith('#'):
                if dedup and out and out[-1] == s:
                    continue
                out.append(s)
    return out


def rel(a, b):
    if a == b:
        return 0.0
    if not (math.isfinite(a) and math.isfinite(b)):
        return math.inf
    s = max(abs(a), abs(b))
    return abs(a - b) / s if s > 0 else abs(a - b)


def compare(pa, pb, tol, exact, dedup, l2=False, columns=None):
    ra, rb = rows(pa, dedup), rows(pb, dedup)
    name = os.path.basename(pa)
    if len(ra) != len(rb):
        print(f"  FAIL {name}: {len(ra)} rows vs {len(rb)}")
        return False
    ncol = len(ra[0].split()) if ra else 0
    use = [c for c in columns if c < ncol] if columns is not None else list(range(ncol))
    same = 0
    cmax = [0.0] * len(use)
    dsq = [0.0] * len(use)
    bsq = [0.0] * len(use)
    for x, y in zip(ra, rb):
        vx, vy = x.split(), y.split()
        if len(vx) != ncol or len(vy) != ncol:
            print(f"  FAIL {name}: column count differs")
            return False
        same += all(vx[c] == vy[c] for c in use)
        for k, c in enumerate(use):
            try:
                a, b = float(vx[c]), float(vy[c])
            except ValueError:
                # A text column (the switch log's regimes and reasons):
                # equal or not.
                if vx[c] != vy[c]:
                    cmax[k] = math.inf
                continue
            cmax[k] = max(cmax[k], rel(a, b))
            dsq[k] += (a - b) ** 2
            bsq[k] += b * b
    if l2:
        meas = [math.sqrt(d / q) if q > 0 else math.sqrt(d) for d, q in zip(dsq, bsq)]
    else:
        meas = cmax
    ok = (same == len(ra)) if exact else all(m <= tol for m in meas)
    cols = " ".join(f"{m:.1e}" for m in meas)
    what = "rel L2 diff" if l2 else "max rel diff"
    sel = "" if columns is None else f" (columns {','.join(map(str, use))})"
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {len(ra)} rows, {same} identical{sel}; "
          f"{what} per column: {cols}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir_a')
    ap.add_argument('prefix_a')
    ap.add_argument('dir_b')
    ap.add_argument('prefix_b')
    ap.add_argument('--exact', action='store_true')
    ap.add_argument('--tol', type=float, default=0.0)
    ap.add_argument('--checkpoints', action='store_true')
    ap.add_argument('--dedup', action='store_true')
    ap.add_argument('--l2', action='store_true')
    ap.add_argument('--station-columns', default=None)
    a = ap.parse_args()
    st_cols = ([int(c) for c in a.station_columns.split(',')]
               if a.station_columns else None)

    files_a = sorted(glob.glob(os.path.join(a.dir_a, a.prefix_a + '_*.txt')))
    files_a = [f for f in files_a if '_checkpoint_r' not in f]
    if not files_a:
        print(f"FAILED: no {a.prefix_a}_*.txt files in {a.dir_a}")
        return 1
    ok = True
    n = 0
    for fa in files_a:
        suffix = os.path.basename(fa)[len(a.prefix_a):]
        fb = os.path.join(a.dir_b, a.prefix_b + suffix)
        if not os.path.exists(fb):
            if suffix == '_catalog.txt':
                print(f"  skip {os.path.basename(fa)} (no counterpart)")
                continue
            print(f"  FAIL {os.path.basename(fa)}: {fb} missing")
            ok = False
            continue
        cols = st_cols if '_fltst_' in suffix else None
        ok = compare(fa, fb, a.tol, a.exact, a.dedup, a.l2, cols) and ok
        n += 1
    if a.checkpoints:
        cka = sorted(glob.glob(os.path.join(a.dir_a, a.prefix_a + '_checkpoint_r*.txt')))
        if not cka:
            print("  FAIL no checkpoint files in " + a.dir_a)
            ok = False
        for fa in cka:
            fb = os.path.join(a.dir_b, a.prefix_b + os.path.basename(fa)[len(a.prefix_a):])
            with open(fa, 'rb') as x:
                da = x.read()
            same = os.path.exists(fb) and open(fb, 'rb').read() == da
            print(f"  {'ok  ' if same else 'FAIL'} {os.path.basename(fa)}: "
                  f"{'identical' if same else 'differs'}")
            ok = ok and same
            n += 1
    print(f"{'PASSED' if ok else 'FAILED'}: {n} file pairs compared")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
