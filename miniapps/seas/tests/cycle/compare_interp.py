#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# compare_interp.py — fully dynamic BP5 cycle: compare the time series of two
# runs whose rows fall at different times (local time stepping against
# global stepping), file by file.  Standard library only.
#
#   python3 tests/cycle/compare_interp.py DIR_A PREFIX_A DIR_B PREFIX_B
#           [--tol T] [--files GLOB] [--t-min T0] [--t-max T1]
#           [--scale column|quantity]
#
# Files are paired by what follows the prefix (as compare_runs.py), among
# PREFIX_A_<GLOB>.txt (default: the station, off-fault and source files).
# B's columns are interpolated linearly to the times of A's rows inside
# both files' common time range (and [T0, T1]); a column's measure is
#   max_t |a(t) - b(t)| / scale
# with --scale column (the default): scale = max_t |a(t)|, the column's
# own; with --scale quantity: the largest |a| of that quantity over all
# its components and all compared files of its kind (fltst, blkst,
# source, ...), in each pair's compared time range, with log10 columns
# compared as 10^a and 10^b.  The
# quantity is the column name of the "# Columns:" header without unit,
# log10() and component (_strike, _dip, _1, _2, _3): slip, V, tau, state,
# disp, vel, Vmax_A, moment_rate_A.  So a component that stays at round-off
# (the dip slip of BP5, an off-fault station the waves have not reached)
# is measured against the size of the quantity, not against its own noise.
# Repeated times keep their first row.  Exit code 0 = every column of
# every pair <= T (default 1e-2).

import argparse
import bisect
import glob
import os
import re
import sys


def columns(path):
    """(quantity, log10?) per column from the '# Columns:' header, or None."""
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# Columns:'):
                out = []
                for name in line[len('# Columns:'):].split(','):
                    name = re.sub(r'\([^()]*\)\s*$', '', name.strip())   # unit
                    lg = name.startswith('log10(') and name.endswith(')')
                    if lg:
                        name = name[len('log10('):-1]
                    name = re.sub(r'_(strike|dip|1|2|3)$', '', name)
                    out.append((name, lg))
                return out
    return None


def table(path):
    t, rows = [], []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            v = [float(x) for x in s]
            if t and v[0] <= t[-1]:
                continue
            t.append(v[0])
            rows.append(v)
    return t, rows


def kind(path, prefix):
    rest = os.path.basename(path)[len(prefix) + 1:]
    return rest.split('_')[0].split('.')[0]


def quantity_scales(pairs, prefix, t0, t1):
    """{(kind, quantity): the largest |value| over the A rows of the pairs,
    each in its compared time range (as compare())}."""
    sc = {}
    for fa, fb in pairs:
        cols = columns(fa)
        if cols is None:
            continue
        ta, ra = table(fa)
        tb, _ = table(fb)
        if not ta or len(tb) < 2:
            continue
        lo, hi = max(ta[0], tb[0], t0), min(ta[-1], tb[-1], t1)
        k = kind(fa, prefix)
        for c in range(1, min(len(cols), len(ra[0]))):
            q, lg = cols[c]
            v = [(10.0 ** r[c] if lg else abs(r[c])) for t, r in zip(ta, ra) if lo <= t <= hi]
            if v:
                sc[(k, q)] = max(sc.get((k, q), 0.0), max(v))
    return sc


def compare(pa, pb, tol, t0, t1, scales=None, prefix=''):
    ta, ra = table(pa)
    tb, rb = table(pb)
    name = os.path.basename(pa)
    cols = columns(pa) if scales is not None else None
    if scales is not None and (cols is None or len(cols) != len(ra[0] if ra else [])):
        print(f"  FAIL {name}: no '# Columns:' header matching the rows (--scale quantity)")
        return False
    if len(ta) < 1 or len(tb) < 2:
        print(f"  FAIL {name}: {len(ta)} and {len(tb)} rows")
        return False
    lo, hi = max(ta[0], tb[0], t0), min(ta[-1], tb[-1], t1)
    ncol = len(ra[0])
    if len(rb[0]) != ncol:
        print(f"  FAIL {name}: {ncol} columns vs {len(rb[0])}")
        return False
    dmax = [0.0] * ncol
    amax = [0.0] * ncol
    n = 0
    for t, r in zip(ta, ra):
        if t < lo or t > hi:
            continue
        j = min(max(bisect.bisect_right(tb, t) - 1, 0), len(tb) - 2)
        w = (t - tb[j]) / (tb[j + 1] - tb[j])
        for c in range(1, ncol):
            b = (1.0 - w) * rb[j][c] + w * rb[j + 1][c]
            if cols is not None and cols[c][1]:
                # a log10 column: 10^b of the interpolated log, as the rows
                dmax[c] = max(dmax[c], abs(10.0 ** r[c] - 10.0 ** b))
            else:
                dmax[c] = max(dmax[c], abs(r[c] - b))
            amax[c] = max(amax[c], abs(r[c]))
        n += 1
    if cols is not None:
        k = kind(pa, prefix)
        amax = [0.0] + [scales[(k, cols[c][0])] for c in range(1, ncol)]
    meas = [d / a if a > 0 else d for d, a in zip(dmax[1:], amax[1:])]
    ok = n > 0 and all(m <= tol for m in meas)
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {n} rows in [{lo:.6g}, {hi:.6g}] s; "
          f"max |a-b| / max |a| per column: {' '.join(f'{m:.1e}' for m in meas)}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir_a')
    ap.add_argument('prefix_a')
    ap.add_argument('dir_b')
    ap.add_argument('prefix_b')
    ap.add_argument('--tol', type=float, default=1e-2)
    ap.add_argument('--files', default='{fltst,blkst,source}')
    ap.add_argument('--t-min', type=float, default=-1e300)
    ap.add_argument('--t-max', type=float, default=1e300)
    ap.add_argument('--scale', choices=('column', 'quantity'), default='column')
    a = ap.parse_args()
    pats = a.files.strip('{}').split(',')
    files = sorted(set(f for p in pats for f in
                       glob.glob(os.path.join(a.dir_a, a.prefix_a + '_' + p + '*.txt'))))
    if not files:
        print(f"FAILED: no {a.prefix_a}_{a.files}*.txt files in {a.dir_a}")
        return 1
    ok = True
    pairs = []
    for fa in files:
        fb = os.path.join(a.dir_b, a.prefix_b + os.path.basename(fa)[len(a.prefix_a):])
        if not os.path.exists(fb):
            print(f"  FAIL {os.path.basename(fa)}: {fb} missing")
            ok = False
            continue
        pairs.append((fa, fb))
    scales = (quantity_scales(pairs, a.prefix_a, a.t_min, a.t_max)
              if a.scale == 'quantity' else None)
    for fa, fb in pairs:
        ok = compare(fa, fb, a.tol, a.t_min, a.t_max, scales, a.prefix_a) and ok
    print(f"{'PASSED' if ok else 'FAILED'}: {len(files)} file pairs compared"
          + (" (scale: per quantity and file kind)" if scales is not None else ""))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
