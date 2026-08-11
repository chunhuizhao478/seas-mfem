#!/usr/bin/env python3
"""inspect_dump.py -- slice a census failure dump by block and by gate.

The census dumps every cell failing the STRICTEST gate it was given; this pulls
out a sub-population (e.g. only the parent block, only cells below 0.6667) and
prints where they are, so a two-cell residual can be located exactly rather
than waved at.
"""
import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--block", choices=("all", "parent", "collar"), default="all")
    ap.add_argument("--parent-tets", type=int, default=0)
    ap.add_argument("--show", type=int, default=12)
    a = ap.parse_args()

    S = np.load(a.dump)
    B, IDX, DX, VS, R = S["bary"], S["idx"], S["dx"], S["vs"], S["ratio"]
    SURF = S["on_surface"]
    m = R < a.gate
    if a.block == "parent":
        m &= IDX < a.parent_tets
    elif a.block == "collar":
        m &= IDX >= a.parent_tets
    print(f"{a.dump}\n  dump holds {len(R):,} cells (gate {float(S['gate'])})")
    print(f"  selected: {a.block}, ratio < {a.gate} -> {int(m.sum()):,} cells")
    if not m.any():
        return
    b, i, d, v, r, s = B[m], IDX[m], DX[m], VS[m], R[m], SURF[m]
    o = np.argsort(r)
    print(f"\n  worst {min(a.show, len(o))}:")
    print(f"    {'tet':>12}{'x':>11}{'y':>12}{'z':>10}{'dx [m]':>10}"
          f"{'Vs':>8}{'Vs/dx':>9}  surf")
    for k in o[:a.show]:
        print(f"    {i[k]:>12,}{b[k,0]:>11,.0f}{b[k,1]:>12,.0f}{b[k,2]:>10,.0f}"
              f"{d[k]:>10,.0f}{v[k]:>8,.0f}{r[k]:>9.4f}  {'Y' if s[k] else 'n'}")
    print(f"\n  ratio  min {r.min():.4f}  p10 {np.percentile(r,10):.4f}  "
          f"med {np.median(r):.4f}")
    print(f"  depth  min {b[:,2].min():,.0f}  med {np.median(b[:,2]):,.0f}  "
          f"max {b[:,2].max():,.0f}")
    print(f"  on the free surface: {int(s.sum()):,} ({100*s.mean():.1f} %)")
    print(f"  refine factor needed (dx*gate/Vs): med {np.median(d*a.gate/v):.2f}x  "
          f"p99 {np.percentile(d*a.gate/v,99):.2f}x  max {(d*a.gate/v).max():.2f}x")
    # bisection bill: cubes of the linear factor, summed
    bill = np.sum(np.maximum(d * a.gate / v, 1.0) ** 3)
    print(f"  ideal bisection bill (sum of factor^3): {bill:,.0f} cells")


if __name__ == "__main__":
    main()
