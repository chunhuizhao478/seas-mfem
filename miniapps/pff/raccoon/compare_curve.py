#!/usr/bin/env python3
"""Compare a pff reaction curve with a MOOSE CSV (time, <name>).

usage: compare_curve.py PFF_CSV MOOSE_CSV NAME [--rel TOL] [--floor F]

Steps are matched by time. The relative difference |pff - moose| / |moose| is
taken where |moose| > floor (default 1). Exit 1 when it exceeds TOL anywhere,
or when the times do not match. Prints the max difference and the peaks.
"""
import argparse
import csv
import sys


def load(path, name):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return [(float(r["time"]), float(r[name])) for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pff")
    ap.add_argument("moose")
    ap.add_argument("name")
    ap.add_argument("--rel", type=float, default=1e-4)
    ap.add_argument("--floor", type=float, default=1.0)
    a = ap.parse_args()
    p, m = load(a.pff, a.name), load(a.moose, a.name)
    if len(p) != len(m):
        print(f"FAIL: {len(p)} pff rows, {len(m)} MOOSE rows")
        return 1
    worst = (0.0, None)
    for (tp, fp), (tm, fm) in zip(p, m):
        if abs(tp - tm) > 1e-12 * max(1.0, abs(tm)):
            print(f"FAIL: time mismatch {tp} vs {tm}")
            return 1
        if abs(fm) > a.floor:
            r = abs(fp - fm) / abs(fm)
            if r > worst[0]:
                worst = (r, (tm, fp, fm))
    pp, pm = max(p, key=lambda r: r[1]), max(m, key=lambda r: r[1])
    print(f"peak pff   {pp[1]:.10g} at t = {pp[0]:g}")
    print(f"peak MOOSE {pm[1]:.10g} at t = {pm[0]:g}")
    print(f"max rel diff {worst[0]:.3e} at {worst[1]}")
    ok = worst[0] <= a.rel and abs(pp[0] - pm[0]) <= 1e-12 * max(1.0, abs(pm[0]))
    print("PASS" if ok else f"FAIL (tol {a.rel:g})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
