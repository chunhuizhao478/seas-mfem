#!/usr/bin/env python3
"""Check that a short Tandem BP5 QD run and an MFEM-SEAS run solve the same
problem the same way, before any speed comparison.

    python3 jobs/vista_qd_speed/tandem/tandem_check_compare.py \
        --tandem <tandem run dir> --mfem <mfem run dir> [--steps 2000]

The Tandem run is tandem_bp5.sbatch with MAXSTEPS and PROBE_TMAX=0 (a probe
row every step); the MFEM run is any production-deck run of the same
degree, mesh and TS options file that went past that many steps.  For
p=1 that is the committed gold bp5/benchmark_data/golden_vista_144r_p1_1000m
(job 1074099: same mesh, petsc_ts_rk45_tandem.cfg, 144 ranks).  The first
2000 steps cover the first 47 s: the start-up slip burst of the nucleation
patch (initial V = 0.01 m/s), the fast phase of the problem.

Prints:
1. Steps: time after each of the first STEPS accepted steps in both codes
   (Tandem from -ts_monitor in run.log, MFEM from *_global.txt, one row
   per accepted step), the largest relative difference and where it is.
2. Stations: for each of the 10 BP5 stations, Tandem's probe rows up to
   its last step, converted as bp5/convert_tandem_to_scec.py does
   (component 1 = strike; theta = L/V0 exp((psi - f0)/b)), against the
   MFEM station file interpolated linearly at the same times.  Largest
   difference per field (strike slip, log10 V strike, strike shear stress,
   log10 state) and, for slip and stress, the field's largest |value|.

Pass rule: step times agree to 1e-2 relative, and per station slip and
stress differ by at most 1e-2 of their largest |value| and the two log10
fields by at most 1e-2.  It only catches a different problem or setup
(mesh, scenario, signs, options); it is not an accuracy test.

Standard library only.
"""

import argparse
import bisect
import glob
import math
import os
import re
import sys

V0, F0, B = 1.0e-6, 0.6, 0.03
STEP_RTOL = 1e-2
FIELD_TOL = 1e-2
MONITOR_RE = re.compile(r"^\s*(\d+)\s+TS dt\s+(\S+)\s+time\s+(\S+)")


def in_nucleation(x2, x3, eps=1e-3):
    # bp5.lua BP5:in_nucleation (bp5_outside, eps = 1e-3), depth d = x3.
    return (4.0 <= x3 + eps and x3 - eps <= 16.0
            and -30.0 <= x2 + eps and x2 - eps <= -18.0)


def tandem_steps(d):
    out = {}
    with open(os.path.join(d, "run.log")) as f:
        for line in f:
            m = MONITOR_RE.match(line)
            if m:
                out[int(m.group(1))] = float(m.group(3))
    return [out[k] for k in sorted(out) if k > 0]


def mfem_prefix(d):
    pf = [f[:-len("_global.txt")] for f in os.listdir(d) if f.endswith("_global.txt")]
    if len(pf) != 1:
        sys.exit(f"ERROR: {d}: expected one *_global.txt, found {len(pf)}")
    return pf[0]


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip().strip('"') for p in s.split(",")] if "," in s else s.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue  # column header
    return rows


def mfem_steps(d, pfx):
    t = [r[0] for r in load(os.path.join(d, pfx + "_global.txt"))]
    return [x for x in t if x > 0.0]


def tandem_station(path, name):
    m = re.match(r"strk([+-]\d+)dp([+-]\d+)", name)
    x2, x3 = int(m.group(1)), int(m.group(2))
    L = 0.13 if in_nucleation(x2, x3) else 0.14
    out = []
    for r in load(path):
        if len(r) < 8:
            sys.exit(f"ERROR: {path}: {len(r)} columns, expected 8 or 9")
        t, psi, slip1, tau1, v1 = r[0], r[1], r[3], r[5], r[7]
        logv = math.log10(max(abs(v1), 1e-300))
        logth = math.log10(L / V0) + (psi - F0) / (B * math.log(10.0))
        out.append([t, slip1, logv, tau1, logth])
    return out


def interp(rows, times, col, t):
    k = bisect.bisect_left(times, t)
    if k <= 0:
        return rows[0][col]
    if k >= len(rows):
        return rows[-1][col]
    t0, t1 = times[k - 1], times[k]
    w = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
    return rows[k - 1][col] * (1 - w) + rows[k][col] * w


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tandem", required=True)
    ap.add_argument("--mfem", required=True)
    ap.add_argument("--steps", type=int, default=2000)
    a = ap.parse_args()
    ok = True

    ts = tandem_steps(a.tandem)[:a.steps]
    pfx = mfem_prefix(a.mfem)
    ms = mfem_steps(a.mfem, pfx)[:a.steps]
    n = min(len(ts), len(ms))
    print(f"1. Steps: Tandem {len(ts)}, MFEM {len(ms)} (compared {n})")
    if n < a.steps:
        print(f"   FAIL: fewer than {a.steps} steps in one run")
        ok = False
    worst, at = 0.0, 0
    for k in range(n):
        rel = abs(ts[k] - ms[k]) / max(abs(ms[k]), 1e-300)
        if rel > worst:
            worst, at = rel, k + 1
    if n:
        print(f"   t after step {n}: Tandem {ts[n-1]:.9e} s, MFEM {ms[n-1]:.9e} s")
        print(f"   largest relative difference {worst:.3e} at step {at}")
        if worst > STEP_RTOL:
            ok = False
    t_end = ts[n - 1] if n else 0.0

    print("2. Stations (strike components), Tandem rows up to its last compared step:")
    print(f"   {'station':14s} {'rows':>5s} {'d slip':>10s} {'|slip|':>10s} {'d logV':>9s}"
          f" {'d tau':>10s} {'|tau|':>9s} {'d logth':>9s}")
    fpfx = "fltst_"
    tfiles = sorted(glob.glob(os.path.join(a.tandem, fpfx + "strk*")))
    if len(tfiles) != 10:
        print(f"   FAIL: {len(tfiles)} Tandem probe files, expected 10")
        ok = False
    for tf in tfiles:
        name = re.search(r"(strk[+-]\d+dp[+-]\d+)", os.path.basename(tf)).group(1)
        mf = os.path.join(a.mfem, f"{pfx}_fltst_{name}.txt")
        if not os.path.exists(mf):
            print(f"   {name}: FAIL, no {mf}")
            ok = False
            continue
        trows = [r for r in tandem_station(tf, name) if r[0] <= t_end * (1 + 1e-12)]
        mrows = load(mf)
        mt = [r[0] for r in mrows]
        # MFEM columns: t slip_s slip_d logV_s logV_d tau_s tau_d log10state
        cols = {1: 1, 2: 3, 3: 5, 4: 7}
        d = {c: 0.0 for c in cols}
        mag = {1: 0.0, 3: 0.0}
        for r in trows:
            for c, mc in cols.items():
                d[c] = max(d[c], abs(r[c] - interp(mrows, mt, mc, r[0])))
            mag[1] = max(mag[1], abs(r[1]))
            mag[3] = max(mag[3], abs(r[3]))
        bad = (d[1] > FIELD_TOL * mag[1] + 1e-12 or d[3] > FIELD_TOL * mag[3]
               or d[2] > FIELD_TOL or d[4] > FIELD_TOL or not trows)
        ok = ok and not bad
        print(f"   {name:14s} {len(trows):5d} {d[1]:10.3e} {mag[1]:10.3e} {d[2]:9.2e}"
              f" {d[3]:10.3e} {mag[3]:9.3e} {d[4]:9.2e}{'  FAIL' if bad else ''}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
