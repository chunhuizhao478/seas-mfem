#!/usr/bin/env python3
"""Compare pff qs_driver nodal output with a RACCOON exodus gold file.

    python3 raccoon/compare_gold.py GOLD.e NODAL.txt STEP [--rel 1e-4 --abs-zero 1e-9]

NODAL.txt is a qs_driver nodal file (x y disp_x disp_y d per vertex).  Nodes
are matched by coordinates.  A value fails as in exodiff: |diff| > abs_zero
and |diff|/|gold| > rel (the RACCOON tests use rel_err 1e-4, abs_zero 1e-9).
Exit status 1 when any value fails or a gold node has no match.
Needs numpy and netCDF4.
"""
import argparse
import sys

import netCDF4 as nc
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gold")
    ap.add_argument("nodal")
    ap.add_argument("step", type=int)
    ap.add_argument("--rel", type=float, default=1e-4)
    ap.add_argument("--abs-zero", type=float, default=1e-9)
    a = ap.parse_args()

    d = nc.Dataset(a.gold)
    names = ["".join(ch.decode() for ch in r if ch) for r in d["name_nod_var"][:]]
    g = np.c_[d["coordx"][:], d["coordy"][:]]
    m = np.loadtxt(a.nodal)
    key = {tuple(np.round(p, 10)): i for i, p in enumerate(m[:, :2])}
    try:
        idx = np.array([key[tuple(np.round(p, 10))] for p in g])
    except KeyError as e:
        print(f"gold node {e} has no match in {a.nodal}")
        return 1
    if len(g) != len(m):
        print(f"node count differs: gold {len(g)}, pff {len(m)}")
        return 1
    fails = 0
    for col, nm in [(2, "disp_x"), (3, "disp_y"), (4, "d")]:
        gv = d["vals_nod_var%d" % (names.index(nm) + 1)][a.step, :]
        mv = m[idx, col]
        diff = np.abs(mv - gv)
        rel = diff / np.maximum(np.abs(gv), 1e-300)
        bad = (diff > a.abs_zero) & (rel > a.rel)
        big = np.abs(gv) > a.abs_zero
        print(f"step {a.step} {nm:7s} max|gold| {np.abs(gv).max():.4e} max|diff| {diff.max():.3e} "
              f"max rel (|gold| > abs_zero) {rel[big].max() if big.any() else 0.0:.3e}  failing {bad.sum()}")
        fails += int(bad.sum())
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
