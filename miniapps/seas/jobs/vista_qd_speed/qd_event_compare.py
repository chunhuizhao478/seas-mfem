#!/usr/bin/env python3
"""Compare a BP5 QD run with a reference run through its earthquakes.

QD speed accuracy gate: a run of seas_bp5_full (loops, default or green)
against the gold full-cycle run 1074099 (bp5/benchmark_data/
golden_vista_144r_p1_1000m), over the time the test run covers.

    python3 jobs/vista_qd_speed/qd_event_compare.py \
        --ref  bp5/benchmark_data/golden_vista_144r_p1_1000m bp5_cycle_p1 \
        --test <run dir> <output prefix>

Prints, with no judgement (the pass rule compares these numbers with the
same table for a repeat of the reference path):

1. Events, from the global files (one row per accepted step): an event is
   a run of rows with log10(Vmax) > -3 (max V above 1 mm/s, the BP5
   event rule).  Per event of the reference inside the test run's time
   span: start, end and peak log10(Vmax), test minus reference.
2. Per station, the rows both runs wrote at the same times (relative
   difference of the times below 1e-9), up to the first row where the
   times differ: how many rows that is, the time it reaches, and per
   field the largest difference over those rows relative to the largest
   |value| of the reference field at that station.
3. Per station and event: the strike slip of the event (slip at the
   event's end minus slip at its start, both interpolated linearly in
   time) and the peak log10(V_strike) inside the event, test minus
   reference.  Interpolation stays inside each run's own rows, so the
   station rows' different times after the first difference do not enter.

Standard library only (bp5_gate_check.py does the same): no numpy on the
cluster's python3 is assumed.

Dip-0 stations (the free surface) are listed after the others: their
reference output already depends on the number of ranks.
"""

import argparse
import bisect
import glob
import os
import sys

FIELDS = ["slip_s", "slip_d", "log10V_s", "log10V_d", "tau_s", "tau_d",
          "log10state"]
EVENT_LOG10V = -3.0
TIME_RTOL = 1e-9


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(x) for x in line.split()])
    if not rows:
        sys.exit(f"ERROR: no rows in {path}")
    return rows


def events(glob_rows):
    """(start, end, peak log10 Vmax, open) per run of rows above EVENT_LOG10V."""
    out = []
    i, n = 0, len(glob_rows)
    while i < n:
        if glob_rows[i][1] <= EVENT_LOG10V:
            i += 1
            continue
        j = i
        while j + 1 < n and glob_rows[j + 1][1] > EVENT_LOG10V:
            j += 1
        peak = max(r[1] for r in glob_rows[i:j + 1])
        out.append((glob_rows[i][0], glob_rows[j][0], peak, j + 1 == n))
        i = j + 1
    return out


def station_names(d, prefix):
    names = []
    for p in sorted(glob.glob(os.path.join(d, prefix + "_fltst_*.txt"))):
        names.append(os.path.basename(p)[len(prefix) + 1:-4])
    return names


def aligned_rows(r, s):
    """Number of leading rows whose times agree to TIME_RTOL."""
    n = min(len(r), len(s))
    for i in range(n):
        tr, ts = r[i][0], s[i][0]
        if abs(tr - ts) > TIME_RTOL * max(1.0, abs(tr)):
            return i
    return n


def interp(rows, times, col, t):
    """Linear interpolation of column @a col at time @a t (clamped)."""
    k = bisect.bisect_left(times, t)
    if k <= 0:
        return rows[0][col]
    if k >= len(rows):
        return rows[-1][col]
    t0, t1 = times[k - 1], times[k]
    v0, v1 = rows[k - 1][col], rows[k][col]
    if t1 == t0:
        return v1
    return v0 + (v1 - v0) * (t - t0) / (t1 - t0)


def event_slip_peak(rows, times, a, b):
    """Strike slip from @a a to @a b and peak log10 V_strike in [a, b]."""
    slip = interp(rows, times, 1, b) - interp(rows, times, 1, a)
    lo = bisect.bisect_left(times, a)
    hi = bisect.bisect_right(times, b)
    peak = max((rows[i][3] for i in range(lo, hi)), default=float("nan"))
    return slip, peak


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", nargs=2, required=True, metavar=("DIR", "PREFIX"))
    ap.add_argument("--test", nargs=2, required=True, metavar=("DIR", "PREFIX"))
    args = ap.parse_args()
    (rdir, rpfx), (tdir, tpfx) = args.ref, args.test

    rg = load(os.path.join(rdir, rpfx + "_global.txt"))
    tg = load(os.path.join(tdir, tpfx + "_global.txt"))
    t_end = tg[-1][0]
    print(f"reference {rdir}/{rpfx}: {len(rg)} global rows, to t = {rg[-1][0]:.6e} s")
    print(f"test      {tdir}/{tpfx}: {len(tg)} global rows, to t = {t_end:.6e} s")

    # ---- 1. Events ----
    rev = [e for e in events(rg) if e[0] <= t_end]
    tev = events(tg)
    print("\n== Events (max V > 1 mm/s), test minus reference ==")
    print(f"{'event':>5} {'ref start s':>16} {'d start s':>11} {'ref end s':>16}"
          f" {'d end s':>11} {'ref peak':>9} {'d peak':>10}")
    pairs = []
    for k, e in enumerate(rev):
        if k >= len(tev):
            print(f"{k + 1:>5} {e[0]:16.9e}  (no such event in the test run)")
            continue
        f = tev[k]
        open_note = ""
        if f[3] or e[1] > t_end:
            open_note = "  (not ended by the test run's last row)"
        pairs.append((e, f))
        print(f"{k + 1:>5} {e[0]:16.9e} {f[0] - e[0]:11.3e} {e[1]:16.9e}"
              f" {f[1] - e[1]:11.3e} {e[2]:9.5f} {f[2] - e[2]:10.2e}{open_note}")
    if len(tev) > len(rev):
        print(f"      the test run has {len(tev) - len(rev)} more event(s)")

    # ---- 2. and 3. Stations ----
    names = station_names(rdir, rpfx)
    tnames = set(station_names(tdir, tpfx))
    missing = [n for n in names if n not in tnames]
    if missing:
        sys.exit(f"ERROR: test run has no station file for {missing}")
    names.sort(key=lambda n: (n.endswith("dp+00"), n))

    print("\n== Stations: rows at identical times, max |diff| / max |ref| ==")
    print(f"{'station':<22} {'rows':>6} {'to t s':>12} "
          + " ".join(f"{f:>10}" for f in FIELDS))
    rows = {}
    for n in names:
        r = load(os.path.join(rdir, f"{rpfx}_{n}.txt"))
        s = load(os.path.join(tdir, f"{tpfx}_{n}.txt"))
        # The reference may run much longer: keep its rows inside the
        # test run's span (the test run's last row is an end-of-run row).
        r = [x for x in r if x[0] <= t_end * (1 + TIME_RTOL)]
        rows[n] = (r, s)
        k = aligned_rows(r, s)
        rel = []
        for c in range(1, 8):
            scale = max((abs(x[c]) for x in r), default=0.0)
            d = max((abs(s[i][c] - r[i][c]) for i in range(k)), default=0.0)
            rel.append(d / scale if scale > 0 else d)
        to_t = r[k - 1][0] if k else 0.0
        print(f"{n:<22} {k:>6} {to_t:12.5e} "
              + " ".join(f"{x:10.2e}" for x in rel))

    if pairs:
        print("\n== Stations per event: event strike slip and peak log10 V_strike,"
              " test minus reference ==")
        hdr = f"{'station':<22}"
        for k in range(len(pairs)):
            hdr += f" {'ev' + str(k + 1) + ' slip m':>14} {'d slip':>10} {'d peak':>10}"
        print(hdr)
        for n in names:
            r, s = rows[n]
            rt = [x[0] for x in r]
            st = [x[0] for x in s]
            line = f"{n:<22}"
            for e, f in pairs:
                rs, rp = event_slip_peak(r, rt, e[0], e[1])
                ss, sp = event_slip_peak(s, st, f[0], f[1])
                line += f" {rs:14.6e} {ss - rs:10.2e} {sp - rp:10.2e}"
            print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
