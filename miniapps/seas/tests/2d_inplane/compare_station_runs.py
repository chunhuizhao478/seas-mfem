#!/usr/bin/env python3
"""Compare the fault-station files of two Nitsche driver runs.

    python3 tests/2d_inplane/compare_station_runs.py DIR_A DIR_B [--rtol R] [--common]

Every <prefix>_station_<name>.dat in DIR_A (not the _copies.dat envelopes)
must have a file of the same name in DIR_B with the same "# Columns:"
header and the same time column, and in every other column (named by that
header, e.g. slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path) the largest
|a - b| over the station's rows must be at most R times the largest |a| of
that column over all stations (R = 0: bit for bit).  Prints one line per
column and exits 1 when a file is missing, a header or the times differ or
a column is over the tolerance.  With --common the headers may differ: the
columns both runs have are compared by name and the others are listed
(e.g. a run with thermal pressurization, columns T and p, against one
without).  Plain Python (no numpy), so it runs wherever the make targets
run.

Used by test-2d-driver (tests/2d_inplane/2d_inplane.mk): the TPV205-2D deck
on 1 rank against 4 ranks; by test-2d-driver-tp: TPV105-2D on 1 rank
against 4, and Lambda = 0 against no thermal pressurization (--common, R =
0).  Scope: Chunhui Zhao, 2026-10-26 (2D in-plane design rev 3, phase C
part 2, phase D part 2).
"""
import argparse
import glob
import os
import sys

def read(path):
    """The column names after "time" (from the "# Columns:" header) and the rows."""
    cols, out = None, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("# Columns:"):
                cols = s.split()[2:]
            elif s and not s.startswith("#"):
                out.append([float(x) for x in s.split()])
    if cols is None or cols[0] != "time":
        sys.exit(f"{path}: no \"# Columns: time ...\" header")
    return cols[1:], out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--rtol", type=float, default=0.0)
    ap.add_argument("--common", action="store_true",
                    help="compare the columns both runs have, by name")
    args = ap.parse_args()
    files = sorted(f for f in glob.glob(os.path.join(args.a, "*_station_*.dat"))
                   if not f.endswith("_copies.dat"))
    if not files:
        sys.exit(f"no station files in {args.a}")
    bad = []
    COLS = read(files[0])[0]
    COLS_B = COLS
    if args.common:
        fb0 = os.path.join(args.b, os.path.basename(files[0]))
        if not os.path.exists(fb0):
            sys.exit(f"{os.path.basename(files[0])} missing in {args.b}")
        COLS_B = read(fb0)[0]
        COLS = [c for c in COLS if c in COLS_B]
    diff = [0.0] * len(COLS)
    scale = [0.0] * len(COLS)
    for fa in files:
        fb = os.path.join(args.b, os.path.basename(fa))
        if not os.path.exists(fb):
            bad.append(f"{os.path.basename(fa)} missing in {args.b}")
            continue
        (ca, ra), (cb, rb) = read(fa), read(fb)
        if (not args.common and (ca != COLS or cb != COLS)) or \
           (args.common and (not set(COLS) <= set(ca) or cb != COLS_B)):
            bad.append(f"{os.path.basename(fa)}: the columns differ ({' '.join(ca)} / {' '.join(cb)})")
            continue
        if len(ra) != len(rb) or any(x[0] != y[0] for x, y in zip(ra, rb)):
            bad.append(f"{os.path.basename(fa)}: the time columns differ ({len(ra)} and {len(rb)} rows)")
            continue
        ia = [ca.index(c) + 1 for c in COLS]
        ib = [cb.index(c) + 1 for c in COLS]
        for x, y in zip(ra, rb):
            for k in range(len(COLS)):
                diff[k] = max(diff[k], abs(x[ia[k]] - y[ib[k]]))
                scale[k] = max(scale[k], abs(x[ia[k]]))
    print(f"  {len(files)} station files, {args.a} against {args.b}, rtol {args.rtol:g}")
    if args.common:
        only = [c for c in read(files[0])[0] if c not in COLS] + [c for c in COLS_B if c not in COLS]
        print(f"    columns in one run only (not compared): {' '.join(only) if only else 'none'}")
    for k, c in enumerate(COLS):
        rel = diff[k] / scale[k] if scale[k] > 0.0 else diff[k]
        ok = rel <= args.rtol
        print(f"    {c:9s} max |difference| {diff[k]:.3e}, relative {rel:.3e}{'' if ok else '  OVER'}")
        if not ok:
            bad.append(f"{c}: relative difference {rel:.3e} > {args.rtol:g}")
    for b in bad:
        print("  " + b)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
