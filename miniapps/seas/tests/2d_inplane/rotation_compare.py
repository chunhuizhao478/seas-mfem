#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/rotation_compare.py — 2D in-plane design rev 3, phase A:
# the comparison of the BP5-mini rotation run (tests/2d_inplane/rotation_run.sh).
#
#   rotation_compare.py A B [--rotated] [--tol T] [--elem-dofs N] [--skip-checkpoint]
#
# A and B are two cycle-driver output folders.  With --rotated, B is the run
# rotated rigidly by R (rotate_case.py's R) and its global vectors are
# compared with R times A's; without it (a control, e.g. 1 rank against 2)
# R = I.  A file's columns are compared in groups, each group's largest
# difference divided by the group's largest |value| in A:
#   a vector's components together: (slip_strike, slip_dip), (V_strike,
#   V_dip) (from the log10 columns), (tau_strike, tau_dip); the off-fault
#   stations' (disp_1, disp_2, disp_3) and (vel_1, vel_2, vel_3) (along the
#   deck's volume_axes, rotated with the deck); a profile's slip_2 and
#   slip_3 (and stress_2 and stress_3) files together; the checkpoint's
#   per-DOF (s_dip, s_strike), traction and slip-rate pairs, and its
#   displacement (byNODES) as a global vector: B = R A
#   every other column alone (time, log10(Vmax), log10(state), psi,
#   t_rupture, moment rate, ...); the fault table's coordinates (x, y, z)
#   as a global vector, B = R A
# A component of a vector that is zero in the physics (slip_dip in BP5 is
# roundoff, about 1e-23 m against 1e-11 m) is thereby measured against the
# vector's size, not its own roundoff.  Text columns (the switch log's
# regimes and reasons) must be equal.
#
# The fault table's rows are matched by position (B's rows whose R^T x is
# within 1e-6 of the table's extent of A's x; of several rows at one
# position, the one nearest in the other columns), since a run on another
# rank count writes them in another order.  The checkpoint's displacement
# (N = --elem-dofs DOFs per element, the elements in mesh order) is
# compared element by element, A's DOFs matched to the element's DOFs in B
# by the permutation that fits best: MFEM's mesh loader orders each tet's
# vertices by its edges' lengths (Mesh::MarkTetMeshForRefinement), and on
# the rotated mesh edges of equal length differ in the last bit, so some
# tets list their vertices, and with them their DOFs, in another order.
# The number of such elements is printed.  --skip-checkpoint leaves the
# checkpoint out (a control on another rank count writes one file per
# rank).
#
# Exit 1 if a group exceeds --tol (default 1e-6) or a file is missing; the
# largest value of each file is printed, and the largest of all on the last
# line ("largest X").  The checkpoint's header values are compared alone
# (time, dt, the rejection count, ...), apart from the commit and the other
# files' byte sizes.
#
# Scope (Chunhui Zhao, 2026-10-25).
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rotate_case import AXIS, DEG, rotation  # noqa: E402

PAIRS = [("slip_strike", "slip_dip"), ("log10(V_strike)", "log10(V_dip)"),
         ("tau_strike", "tau_dip")]
TRIPLES = [("disp_1", "disp_2", "disp_3"), ("vel_1", "vel_2", "vel_3"), ("x", "y", "z")]
GLOBAL = {("x", "y", "z")}   # global coordinates: rotate A's


def columns(path):
    names, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            m = re.match(r"#\s*Columns:\s*(.*)$", line)
            if m:
                names = [re.sub(r"\((m|s|m/s|MPa|m\^2|N\*m/s)\)$", "", c.strip())
                         for c in m.group(1).split(",")]
            continue
        if line.strip():
            rows.append(line.split())
    return names, rows


def num(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def group_diff(a_cols, b_cols, R=None, log=False):
    """Largest |B - A'| over rows, divided by the largest |A'|, A' = R A for a
    global vector, the components taken together."""
    n = len(a_cols[0])
    worst, scale = 0.0, 0.0
    for i in range(n):
        a = [c[i] for c in a_cols]
        b = [c[i] for c in b_cols]
        if log:
            a = [10.0 ** v for v in a]
            b = [10.0 ** v for v in b]
        if R is not None:
            a = [sum(R[r][k] * a[k] for k in range(3)) for r in range(3)]
        scale = max(scale, math.sqrt(sum(v * v for v in a)))
        worst = max(worst, math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))))
    return worst / scale if scale > 0 else (0.0 if worst == 0 else math.inf)


def compare_table(fa, fb, R, extra_scale=None):
    """Grouped comparison of two column files; returns (worst, detail, text_ok)."""
    na, ra = columns(fa)
    nb, rb = columns(fb)
    if na != nb or len(ra) != len(rb) or any(len(x) != len(y) for x, y in zip(ra, rb)):
        return math.inf, "shape or columns differ", False
    if na is None:
        na = [f"c{j}" for j in range(len(ra[0]) if ra else 0)]
    ncol = len(na)
    if all(c in na for c in ("x", "y", "z")):
        ra, rb = match_rows(ra, rb, [na.index(c) for c in ("x", "y", "z")], R)
        if ra is None:
            return math.inf, "rows do not match by position", False
    text_ok = True
    numeric = []
    for j in range(ncol):
        ca = [num(r[j]) for r in ra]
        cb = [num(r[j]) for r in rb]
        if any(v is None for v in ca) or any(v is None for v in cb):
            text_ok &= all(r1[j] == r2[j] for r1, r2 in zip(ra, rb))
            numeric.append(None)
        else:
            numeric.append((ca, cb))
    done, results = set(), []
    for grp in PAIRS + TRIPLES:
        if all(g in na for g in grp):
            idx = [na.index(g) for g in grp]
            if any(numeric[j] is None for j in idx):
                continue
            Rg = R if grp in GLOBAL else None
            d = group_diff([numeric[j][0] for j in idx], [numeric[j][1] for j in idx], Rg,
                           log=grp[0].startswith("log10(V_"))
            results.append(("+".join(grp), d))
            done.update(idx)
    for j in range(ncol):
        if j in done or numeric[j] is None:
            continue
        ca, cb = numeric[j]
        s = max((abs(v) for v in ca), default=0.0)
        if extra_scale is not None and re.match(r"(slip|stress)_[23]_\d+$", na[j]):
            s = extra_scale
        w = max((abs(x - y) for x, y in zip(ca, cb)), default=0.0)
        results.append((na[j], w / s if s > 0 else (0.0 if w == 0 else math.inf)))
    worst = max(results, key=lambda r: r[1]) if results else ("-", 0.0)
    return worst[1], f"{worst[0]}", text_ok


def match_rows(ra, rb, ix, R):
    """B's rows reordered to A's by position (R^T x_B against x_A, within 1e-6
    of A's extent).  A position can hold several rows (the DOFs of adjacent
    faces at a shared node): among B's unused rows there, the one nearest in
    the other numeric columns (each scaled by its largest |value|)."""
    xa = [[float(r[j]) for j in ix] for r in ra]
    xb = [[float(r[j]) for j in ix] for r in rb]
    xb = [[sum(R[k][i] * x[k] for k in range(3)) for i in range(3)] for x in xb]
    ext = max(max(abs(v) for v in x) for x in xa) if xa else 1.0
    h = 1e-6 * ext
    vcols = [j for j in range(len(ra[0])) if j not in ix] if ra else []
    vsc = [max(abs(float(r[j])) for r in ra) or 1.0 for j in vcols]
    val = lambda r: [float(r[j]) / c for j, c in zip(vcols, vsc)]
    key = lambda x: tuple(int(math.floor(v / h)) for v in x)
    where = {}
    for i, x in enumerate(xb):
        where.setdefault(key(x), []).append(i)
    used, order = set(), []
    near = [(a, b, c) for a in (0, -1, 1) for b in (0, -1, 1) for c in (0, -1, 1)]
    for r, x in zip(ra, xa):
        k0 = key(x)
        cand = [i for d in near for i in where.get(tuple(u + e for u, e in zip(k0, d)), [])
                if i not in used and max(abs(u - v) for u, v in zip(x, xb[i])) <= h]
        if not cand:
            return None, None
        va = val(r)
        best = min(cand, key=lambda i: max(abs(u - v) for u, v in zip(va, val(rb[i]))))
        used.add(best)
        order.append(best)
    return ra, [rb[i] for i in order]


def displacement_diff(a, b, R, nd):
    """byNODES displacement: per element of nd DOFs, A's DOFs (rotated) matched
    to B's DOFs of the same element by the permutation with the smallest
    largest distance (all nd! for nd <= 6, otherwise nearest-first); (scaled
    difference, number of elements whose DOFs come in another order, i.e.
    whose identity matching is off by more than 1e-6 of the field's size)."""
    import itertools
    n = len(a) // 3
    if nd <= 0 or n % nd:
        return math.inf, -1
    va = [[sum(R[r][k] * (a[k * n + i]) for k in range(3)) for r in range(3)] for i in range(n)]
    vb = [[b[i], b[n + i], b[2 * n + i]] for i in range(n)]
    scale = max(math.sqrt(sum(v * v for v in x)) for x in va)
    if scale == 0.0:
        return (0.0, 0) if all(v == 0.0 for x in vb for v in x) else (math.inf, 0)
    dist = lambda x, y: math.sqrt(sum((u - v) ** 2 for u, v in zip(x, y)))
    perms = list(itertools.permutations(range(nd))) if nd <= 6 else None
    worst, permuted = 0.0, 0
    for e in range(n // nd):
        A = va[e * nd:(e + 1) * nd]
        B = vb[e * nd:(e + 1) * nd]
        ident = max(dist(A[j], B[j]) for j in range(nd))
        if perms is not None:
            best = min(max(dist(A[j], B[q[j]]) for j in range(nd)) for q in perms)
        else:
            free, best = set(range(nd)), 0.0
            for j in range(nd):
                k = min(free, key=lambda k: dist(A[j], B[k]))
                free.discard(k)
                best = max(best, dist(A[j], B[k]))
        if ident > 1e-6 * scale:
            permuted += 1
        worst = max(worst, best)
    return worst / scale, permuted


def profile_scale(folder, name):
    """The largest |value| of a profile's _2 and _3 files together (data rows,
    value columns)."""
    m = re.match(r"(b_profile_\w+_(slip|stress))_[23]\.txt$", name)
    s = 0.0
    for k in ("2", "3"):
        f = os.path.join(folder, f"{m.group(1)}_{k}.txt")
        names, rows = columns(f)
        for r in rows[1:]:   # the first row: the points' coordinates
            s = max(s, max(abs(float(v)) for v in r[2:]))
    return s


def read_checkpoint(path):
    lines = open(path).read().split("\n")
    blocks, scalars, i = {}, [], 0
    while i < len(lines):
        f = lines[i].split()
        if len(f) == 2 and f[0].endswith("_size"):
            n = int(f[1])
            blocks[f[0][:-5]] = [float(v) for v in lines[i + 1:i + 1 + n]]
            i += 1 + n
            continue
        if f:
            scalars.append(f)
        i += 1
    return blocks, scalars


def compare_checkpoint(fa, fb, R, nd):
    ba, sa = read_checkpoint(fa)
    bb, sb = read_checkpoint(fb)
    res = []
    for k in ba:
        a, b = ba[k], bb.get(k)
        if b is None or len(a) != len(b):
            return math.inf, f"block {k} differs in size"
        if not a:
            continue
        if k == "state":   # [s_dip, s_strike, psi] per DOF
            res.append(("state slip", group_diff([a[0::3], a[1::3]], [b[0::3], b[1::3]])))
            res.append(("state psi", group_diff([a[2::3]], [b[2::3]])))
        elif k == "displacement":   # byNODES: x block, y block, z block
            d, perm = displacement_diff(a, b, R, nd)
            res.append((f"displacement, {perm} element(s) in another DOF order", d))
        elif k in ("traction", "slip_rate"):   # (dip, strike) per DOF
            res.append((k, group_diff([a[0::2], a[1::2]], [b[0::2], b[1::2]])))
        else:
            res.append((k, group_diff([a], [b])))
    if len(sa) != len(sb):
        return math.inf, "header lines differ in number"
    for x, y in zip(sa, sb):
        if x[0] != y[0]:
            return math.inf, f"header {x[0]} vs {y[0]}"
        # the commit, and the byte sizes of the other output files (their
        # headers print the station coordinates, which the rotation changes)
        if x[0] == "commit" or x[0].endswith(".txt"):
            continue
        for u, v in zip(x[1:], y[1:]):
            nu, nv = num(u), num(v)
            if nu is None or nv is None:
                if u != v:
                    res.append((x[0], math.inf))
            elif nu != nv:
                res.append((x[0], abs(nu - nv) / max(abs(nu), 1e-300)))
    worst = max(res, key=lambda r: r[1])
    perm = [w for w, _ in res if w.startswith("displacement")]
    return worst[1], worst[0] + ("" if not perm or perm[0] == worst[0] else "; " + perm[0])


def main(argv):
    args, opts, i = [], {}, 1
    while i < len(argv):
        if argv[i] in ("--tol", "--elem-dofs"):
            opts[argv[i]] = argv[i + 1]
            i += 2
        elif argv[i].startswith("--"):
            opts[argv[i]] = True
            i += 1
        else:
            args.append(argv[i])
            i += 1
    rotated = "--rotated" in opts
    tol = float(opts.get("--tol", 1e-6))
    nd = int(opts.get("--elem-dofs", 0))
    if len(args) != 2:
        sys.exit("usage: rotation_compare.py A B [--rotated] [--tol T] [--elem-dofs N] "
                 "[--skip-checkpoint]")
    A, B = args
    R = rotation(AXIS, DEG) if rotated else [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    names = sorted(f for f in os.listdir(A) if f.startswith("b_") and f.endswith(".txt"))
    worst_all, bad = 0.0, 0
    for name in names:
        fa, fb = os.path.join(A, name), os.path.join(B, name)
        if not os.path.exists(fb):
            print(f"   MISSING  {name}")
            bad += 1
            continue
        if name.startswith("b_checkpoint"):
            if "--skip-checkpoint" in opts:
                continue
            w, what = compare_checkpoint(fa, fb, R, nd)
            text_ok = True
        else:
            es = profile_scale(A, name) if name.startswith("b_profile_") else None
            w, what, text_ok = compare_table(fa, fb, R, es)
        flag = "ok  " if (w <= tol and text_ok) else "FAIL"
        if flag == "FAIL":
            bad += 1
        worst_all = max(worst_all, w)
        print(f"   {flag} {name:34s} {w:9.2e}  ({what}{'' if text_ok else '; text differs'})")
    print(f"   largest {worst_all:.2e} (tolerance {tol:g}), {bad} file(s) over it or missing")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
