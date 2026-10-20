#!/usr/bin/env python3
"""Merge the Nitsche driver's per-rank rupture-time files into a SCEC contour-plot file.

The driver writes <prefix>_rupture_time_r<rank>.csv per rank (x y z t_rupture per owned fault
DOF; t_rupture = the first time |V| > 1e-3 m/s, -1 = never).  A fault node carries one DOF per
adjacent fault face (DG copies), so the copies at one position are merged: the node's rupture
time is the earliest of its copies (the first time any copy's slip rate passes 1 mm/s).  The
output is the SCEC cplot of TPV26/27 description v13, Part 10:

    j  k  t        j along strike (m), k down-dip depth (m, 0 = surface), t (s; 1.0E+09 = never)

with j = x . strike and k = surface_elevation - x . up (general coordinates; TPV26/27: strike
(1, 0, 0), up (0, 0, 1), so j = X and k = -Z).

Usage (from miniapps/seas):
    python3 -I scripts/make_cplot.py --run-dir <output dir> --prefix tpv26 --problem TPV26 \\
        --element-size "200 m, p = 2" --out <output dir>/cplot

Scope: Chunhui Zhao, 2026-10-20 (rheology plan rev 1, phase R1).
"""
import argparse
import glob
import math
import os
import sys

NEVER = 1.0e9


def vec(s):
    v = [float(c) for c in s.split(",")]
    if len(v) != 3:
        raise argparse.ArgumentTypeError(f"expected three comma-separated numbers, got {s!r}")
    n = math.sqrt(sum(c * c for c in v))
    if n == 0.0:
        raise argparse.ArgumentTypeError("zero vector")
    return [c / n for c in v]


def read_rank_files(run_dir, prefix):
    paths = sorted(glob.glob(os.path.join(run_dir, f"{prefix}_rupture_time_r*.csv")))
    if not paths:
        sys.exit(f"make_cplot: no {prefix}_rupture_time_r*.csv in {run_dir}")
    rows = []
    for p in paths:
        with open(p) as f:
            for ln, line in enumerate(f, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                s = line.split()
                if len(s) != 4:
                    sys.exit(f"make_cplot: {p}:{ln}: expected 'x y z t_rupture'")
                rows.append([float(c) for c in s])
    return paths, rows


def merge(rows, tol):
    """One entry per position (positions within tol share a key); earliest positive time."""
    nodes = {}
    for x, y, z, t in rows:
        key = (round(x / tol), round(y / tol), round(z / tol))
        t = t if t >= 0.0 else NEVER
        if key not in nodes:
            nodes[key] = [x, y, z, t]
        else:
            nodes[key][3] = min(nodes[key][3], t)
    return list(nodes.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--problem", default="TPV26")
    ap.add_argument("--author", default="Chunhui Zhao")
    ap.add_argument("--code", default="seas-mfem Nitsche dynamic rupture")
    ap.add_argument("--code-version", default="", help="e.g. the commit hash of the run")
    ap.add_argument("--element-size", default="", help="free text for the header")
    ap.add_argument("--strike", type=vec, default=vec("1,0,0"))
    ap.add_argument("--up", type=vec, default=vec("0,0,1"))
    ap.add_argument("--surface-elevation", type=float, default=0.0)
    ap.add_argument("--tol", type=float, default=1e-6, help="position tolerance for DG copies (m)")
    a = ap.parse_args()

    paths, rows = read_rank_files(a.run_dir, a.prefix)
    nodes = merge(rows, a.tol)
    out = []
    for x, y, z, t in nodes:
        j = x * a.strike[0] + y * a.strike[1] + z * a.strike[2]
        k = a.surface_elevation - (x * a.up[0] + y * a.up[1] + z * a.up[2])
        out.append((j, k + 0.0, t))
    out.sort(key=lambda r: (r[1], r[0]))
    ruptured = sum(1 for r in out if r[2] < NEVER)
    with open(a.out, "w") as f:
        f.write("# Contour-plot file (SCEC TPV26/27 description v13, Part 10)\n")
        f.write(f"# problem={a.problem}\n# author={a.author}\n")
        f.write(f"# code={a.code}\n")
        if a.code_version:
            f.write(f"# code_version={a.code_version}\n")
        if a.element_size:
            f.write(f"# element_size={a.element_size}\n")
        f.write(f"# {len(out)} fault nodes from {len(rows)} DOFs in {len(paths)} rank files; "
                f"a node's time is the earliest of its DG copies; {ruptured} ruptured\n")
        f.write("# Column #1 = Distance along strike (m)\n")
        f.write("# Column #2 = Distance down-dip (m)\n")
        f.write("# Column #3 = Rupture time (s)\n")
        f.write("j k t\n")
        for j, k, t in out:
            f.write(f"{j:.6e} {k:.6e} {t:.6e}\n")
    print(f"wrote {a.out}: {len(out)} nodes ({ruptured} ruptured) from {len(rows)} DOFs, "
          f"{len(paths)} rank files")


if __name__ == "__main__":
    main()
