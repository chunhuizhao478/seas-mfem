#!/usr/bin/env python3
"""Compare the integrator-study runs with the reference run.

QD speed-up plan v2, item A: every run is the production BP5 QD deck from
t = 0 to the same end time with another PETSc TS options file
(jobs/vista_qd_speed/qd_integrator.sbatch).  The reference is today's
tableau at a much tighter tolerance; each run's error is its distance to
it.

    python3 jobs/vista_qd_speed/qd_integrator_compare.py \
        --ref ref <ref run dir> --run a0 <dir> --run a1 <dir> ... \
        --line a0 a3

Prints:

1. Cost per run: tableau, controller and atol (the driver's [TS] line),
   accepted and rejected steps, stage evaluations (the [STAGE-TIMERS]
   end-of-run count, rejected steps included), evaluations per accepted
   step, loop wall time and s/step.
2. Error per run against the reference: per event, start and end time
   differences, and the largest event strike slip difference over the
   stations; at the end time (every run ends exactly there), the largest
   strike slip and log10(state) differences over the stations.
3. Equal accuracy: the two --line runs (today's tableau and controller at
   two tolerances) give stage evaluations as a power of the error.  For
   every other run, the evaluations today's setting would need for the
   same error, divided by the run's own: above 1 the run is cheaper at
   equal accuracy.  Read for the event-2 start and the end-time slip.
   Also the strict rule of the plan: error no larger than the first --line
   run's and fewer evaluations.

Errors near the run-to-run floor (event-2 start 0.2 to 0.9 s between
round-off-level changes, accuracy gate of 2026-10-20) do not separate the
runs; the table marks them.  Standard library only.
"""

import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qd_event_compare import events, event_slip_peak, load, station_names  # noqa: E402

TIME_RTOL = 1e-9
FLOOR_EV_START = 0.9  # s, largest event-2 shift seen from round-off alone

STAGE_RE = re.compile(r"end of run: (\d+) accepted steps, (\d+) stages, "
                      r"loop wall ([0-9.]+) s \(([0-9.eE+-]+) s/step\)")


def prefix_of(d):
    pf = [f[:-len("_global.txt")] for f in os.listdir(d)
          if f.endswith("_global.txt")]
    if len(pf) != 1:
        sys.exit(f"ERROR: {d}: expected one *_global.txt, found {len(pf)}")
    return pf[0]


def read_run(tag, d):
    run = {"tag": tag, "dir": d, "pfx": prefix_of(d)}
    run["glob"] = load(os.path.join(d, run["pfx"] + "_global.txt"))
    ts = os.path.join(d, "ts.txt")
    run["ts"] = open(ts).read().strip().replace("[TS] ", "") \
        if os.path.exists(ts) else "?"
    run["steps"] = run["stages"] = run["wall"] = run["sps"] = None
    st = os.path.join(d, "stage_timers_final.txt")
    if os.path.exists(st):
        m = STAGE_RE.search(open(st).read())
        if m:
            run["steps"], run["stages"] = int(m.group(1)), int(m.group(2))
            run["wall"], run["sps"] = float(m.group(3)), float(m.group(4))
    run["rejected"] = None
    sm = os.path.join(d, "summary.txt")
    if os.path.exists(sm):
        m = re.search(r"Step rejections: (\d+)", open(sm).read())
        if m:
            run["rejected"] = int(m.group(1))
    run["events"] = events(run["glob"])
    run["stations"] = {}
    for n in station_names(d, run["pfx"]):
        rows = load(os.path.join(d, f"{run['pfx']}_{n}.txt"))
        run["stations"][n] = (rows, [x[0] for x in rows])
    return run


def errors(ref, run):
    """Error measures of @a run against @a ref (None where undefined)."""
    e = {"ev": []}
    for k, re_ in enumerate(ref["events"]):
        if k >= len(run["events"]):
            e["ev"].append(None)
            continue
        te = run["events"][k]
        dslip = 0.0
        for n, (r, rt) in ref["stations"].items():
            if n not in run["stations"]:
                sys.exit(f"ERROR: {run['tag']} has no station {n}")
            s, st = run["stations"][n]
            rs, _ = event_slip_peak(r, rt, re_[0], re_[1])
            ss, _ = event_slip_peak(s, st, te[0], te[1])
            dslip = max(dslip, abs(ss - rs))
        e["ev"].append((te[0] - re_[0], te[1] - re_[1], dslip))
    d_slip = d_state = 0.0
    for n, (r, _) in ref["stations"].items():
        s, _ = run["stations"][n]
        tr, tt = r[-1][0], s[-1][0]
        if abs(tr - tt) > TIME_RTOL * max(1.0, abs(tr)):
            sys.exit(f"ERROR: {run['tag']} station {n} ends at t = {tt:.9e},"
                     f" the reference at {tr:.9e}")
        d_slip = max(d_slip, abs(s[-1][1] - r[-1][1]))
        d_state = max(d_state, abs(s[-1][7] - r[-1][7]))
    e["slip_end"], e["state_end"] = d_slip, d_state
    return e


def ev_start(e, k=1):
    """|start difference| of event k (0-based), or None."""
    return abs(e["ev"][k][0]) if len(e["ev"]) > k and e["ev"][k] else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", nargs=2, required=True, metavar=("TAG", "DIR"))
    ap.add_argument("--run", nargs=2, action="append", required=True,
                    metavar=("TAG", "DIR"))
    ap.add_argument("--line", nargs=2, required=True, metavar=("TAG1", "TAG2"),
                    help="today's setting at two tolerances")
    args = ap.parse_args()

    ref = read_run(*args.ref)
    runs = [read_run(t, d) for t, d in args.run]
    tags = [r["tag"] for r in runs]
    for t in args.line:
        if t not in tags:
            sys.exit(f"ERROR: --line run {t} is not among --run")
    t_end = ref["glob"][-1][0]
    print(f"reference {ref['tag']}: {ref['dir']} ({ref['ts']}), to t = {t_end:.6e} s,"
          f" {len(ref['events'])} events, {len(ref['stations'])} stations")

    # ---- 1. Cost ----
    print("\n== 1. Cost ==")
    print(f"{'run':<5} {'steps':>8} {'rejected':>8} {'stages':>9} {'st/step':>7}"
          f" {'wall h':>7} {'s/step':>8}  TS")
    for r in [ref] + runs:
        def f(x, fmt):
            return format(x, fmt) if x is not None else "?"
        sps = r["stages"] / r["steps"] if r["steps"] and r["stages"] else None
        wall = r["wall"] / 3600 if r["wall"] else None
        print(f"{r['tag']:<5} {f(r['steps'], '8d')} {f(r['rejected'], '8d')}"
              f" {f(r['stages'], '9d')} {f(sps, '7.3f')} {f(wall, '7.2f')}"
              f" {f(r['sps'], '8.4f')}  {r['ts']}")

    # ---- 2. Error against the reference ----
    errs = {r["tag"]: errors(ref, r) for r in runs}
    nev = len(ref["events"])
    print("\n== 2. Error against the reference (run minus reference) ==")
    hdr = f"{'run':<5}"
    for k in range(nev):
        hdr += f" {'ev' + str(k + 1) + ' start s':>12} {'end s':>10} {'slip m':>9}"
    hdr += f" {'end slip m':>11} {'end log10 state':>15}"
    print(hdr)
    for r in runs:
        e = errs[r["tag"]]
        line = f"{r['tag']:<5}"
        for v in e["ev"]:
            line += (f" {v[0]:12.3e} {v[1]:10.3e} {v[2]:9.2e}" if v
                     else f" {'missing':>12} {'':>10} {'':>9}")
        line += f" {e['slip_end']:11.3e} {e['state_end']:15.3e}"
        print(line)

    # ---- 3. Equal accuracy ----
    print("\n== 3. Equal accuracy against today's setting"
          f" ({args.line[0]}, {args.line[1]}) ==")
    by = {r["tag"]: r for r in runs}
    l0, l1 = by[args.line[0]], by[args.line[1]]
    metrics = [("event-2 start", lambda e: ev_start(e)),
               ("end-time slip", lambda e: e["slip_end"])]
    print(f"{'run':<5} {'metric':<14} {'error':>10} {'today needs':>12}"
          f" {'run uses':>9} {'gain':>6}  strict rule (vs {args.line[0]})")
    for r in runs:
        if r["tag"] in args.line or not r["stages"]:
            continue
        for name, get in metrics:
            e_r = get(errs[r["tag"]])
            e0, e1 = get(errs[l0["tag"]]), get(errs[l1["tag"]])
            s0, s1 = l0["stages"], l1["stages"]
            note = ""
            if None in (e_r, e0, e1) or not s0 or not s1:
                print(f"{r['tag']:<5} {name:<14} {'?':>10}")
                continue
            if name == "event-2 start" and max(e_r, e0) <= FLOOR_EV_START:
                note = "  (within the run-to-run floor)"
            if e_r > 0 and e0 > 0 and e1 > 0 and e0 != e1:
                p = math.log(s1 / s0) / math.log(e1 / e0)
                need = s0 * (e_r / e0) ** p
                gain = f"{need / r['stages']:6.2f}"
                need = f"{need:12.0f}"
            else:
                need, gain = f"{'?':>12}", f"{'?':>6}"
            strict = "PASS" if (e_r <= e0 and r["stages"] < s0) else "FAIL"
            print(f"{r['tag']:<5} {name:<14} {e_r:10.3e} {need} {r['stages']:9d}"
                  f" {gain}  {strict}{note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
