#!/usr/bin/env python3
"""d2_trace_check.py — GPU dynamic-rupture plan D2: the host-device copies
in an Nsight Systems trace of seas_nitsche_dyn_driver's device step, placed
in the driver's NVTX ranges.

    python3 jobs/vista/d2_trace_check.py REPORT.sqlite --log driver.log [--gate]

REPORT.sqlite is `nsys export --type sqlite` of a CUDA build of
seas_nitsche_dyn_driver run on one rank with --device cuda --matrix-free 1
--stage-timers N (which turns its NVTX ranges on); driver.log is that run's
output.  From the log: the gather copy GATHER_BYTES = 8 x the values the
summary says the device step gathers per step ("fault gather M values per
step", which counts Udot* too with Kelvin-Voigt), and, for scale, the
volume size VOLUME_BYTES = 8 x 3 x nd x the tets of "mesh ...: NE tets;
p = P" (nd = (P+1)(P+2)(P+3)/6; the rank's own on one rank).

Each copy is placed by the CPU time of the CUDA call that issued it (the
loaders of g4_trace_check.py), in the first of these that holds it:

    step1     the first seas:dyn_step range, which also assembles the
              operators and copies them to the device
    step      a later seas:dyn_step range (one ode.Step)
    output    seas:dyn_output (stations, fault output, log)
    outside   set-up and the end of the run

Gate (--gate; plan rev2 §4.3 D2 acceptance, "an Nsight trace shows no
volume-sized copy inside a step, only the planned fault gather"): no
host-to-device or device-to-host copy in step larger than GATHER_BYTES,
which is far below VOLUME_BYTES.  The trace must also show the planned
gather, a device-to-host copy of GATHER_BYTES in step, or it cannot answer.
The copies of the slip to the device (DeviceSlipLoad) are smaller and
reported.  Without --gate the script only reports (traces of several
ranks, where the halo may stage through the host).

Exit: 0 pass (or report only), 1 gate failed, 2 the trace or the log cannot
answer (no copies table, fewer than two seas:dyn_step ranges, the device
step off, no gather seen).
"""

import argparse
import bisect
import collections
import re
import sqlite3
import sys

import g4_trace_check as g4

CONTEXT_ORDER = ["step1", "step", "output", "outside"]
HOST_DEVICE = ("HtoD", "DtoH")


def read_log(path):
    """(device step on, gathered values per step, tets, order) from the
    driver's log; None where the log does not say."""
    text = open(path, errors="replace").read()
    on = "device step on" in text
    m = re.search(r"fault gather (\d+) values per step", text)
    gather = int(m.group(1)) if m else None
    m = re.search(r": (\d+) tets; p = (\d+)", text)
    ne, p = (int(m.group(1)), int(m.group(2))) if m else (None, None)
    return on, gather, ne, p


def step_index(ranges, tid, pid_key, t):
    """Index of the seas:dyn_step range that holds time t (0 = step 1), or
    None.  A copy without a thread looks in every thread of its process."""
    by_tid = ranges.get("seas:dyn_step", {})
    keys = [tid] if tid is not None else \
        [k for k in by_tid if (k >> 24) == pid_key]
    for k in keys:
        spans = by_tid.get(k, [])
        i = bisect.bisect_right(spans, (t, float("inf"))) - 1
        if i >= 0 and spans[i][0] <= t <= spans[i][1]:
            return i
    return None


def context_of(ranges, tid, pid_key, t):
    i = step_index(ranges, tid, pid_key, t)
    if i is not None:
        return "step1" if i == 0 else "step"
    if g4.in_range(ranges, "seas:dyn_output", tid, pid_key, t):
        return "output"
    return "outside"


def mb(nbytes):
    return nbytes / 1.0e6


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sqlite")
    ap.add_argument("--log", required=True, help="the traced run's output")
    ap.add_argument("--gate", action="store_true",
                    help="exit 1 on a copy in a step larger than the gather")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()

    on, gather, ne, p = read_log(args.log)
    broken = []  # the trace or log cannot answer: exit 2 either way
    if not on:
        broken.append("the log does not say 'device step on'")
    if not gather:
        broken.append("the log has no 'fault gather M values per step'")
    gather_bytes = 8 * (gather or 0)
    volume_bytes = 8 * 3 * (p + 1) * (p + 2) * (p + 3) // 6 * ne \
        if ne else None

    db = sqlite3.connect(f"file:{args.sqlite}?mode=ro", uri=True)
    have = g4.tables(db)
    kinds = g4.copy_kinds(db, have)
    ranges = g4.load_ranges(db, have)
    copies, unmatched = g4.load_copies(db, have)

    print(f"Trace: {args.sqlite}")
    print(f"Log: {args.log}")
    print(f"Gathered values per step {gather}: GATHER_BYTES = {gather_bytes} "
          f"({mb(gather_bytes):.3f} MB)")
    if volume_bytes:
        print(f"Volume {ne} tets, p = {p}: VOLUME_BYTES = {volume_bytes} "
              f"({mb(volume_bytes):.3f} MB, {volume_bytes / max(gather_bytes, 1):.0f}"
              f" x GATHER_BYTES)")
    n_step = g4.count_ranges(ranges, "seas:dyn_step")
    n_out = g4.count_ranges(ranges, "seas:dyn_output")
    n_traction = g4.count_ranges(ranges, "seas:dyn_fault_traction")
    print(f"NVTX ranges: dyn_step {n_step}, dyn_output {n_out}, "
          f"dyn_fault_traction {n_traction}")
    if copies is None:
        broken.append("no CUDA memcpy table in the trace (was it run with "
                      "--trace=cuda?)")
        copies = []
    if n_step < 2:
        broken.append("fewer than two seas:dyn_step NVTX ranges (CUDA build "
                      "with nvtx3 headers, run with --stage-timers and nsys "
                      "--trace=nvtx?)")
    if unmatched:
        print(f"Note: {unmatched} of {len(copies)} copies had no matching "
              f"CUDA call; placed by their GPU start time")

    stats = collections.defaultdict(lambda: [0, 0, 0])
    over = []
    gathers = 0
    for t, tid, pid_key, nbytes, kind in copies:
        ctx = context_of(ranges, tid, pid_key, t)
        kname = kinds.get(kind, str(kind))
        s = stats[(ctx, kname)]
        s[0] += 1
        s[1] += nbytes
        s[2] = max(s[2], nbytes)
        if ctx == "step" and kname in HOST_DEVICE:
            if nbytes > gather_bytes:
                over.append((nbytes, kname, t))
            if kname == "DtoH" and nbytes == gather_bytes:
                gathers += 1

    steps = max(n_step - 1, 1)
    print()
    print(f"{'context':<8} {'kind':<6} {'copies':>8} {'total MB':>11} "
          f"{'max bytes':>11} {'per step':>9}")
    for ctx in CONTEXT_ORDER:
        for (c, kname), (n, tot, mx) in sorted(stats.items()):
            if c != ctx:
                continue
            per = f"{n / steps:.2f}" if ctx == "step" else ""
            print(f"{ctx:<8} {kname:<6} {n:>8} {mb(tot):>11.3f} {mx:>11} "
                  f"{per:>9}")

    print()
    print(f"Steps after step 1: {n_step - 1}; device-to-host copies of "
          f"GATHER_BYTES in them (the planned gather): {gathers}")
    print(f"Host-device copies in them larger than GATHER_BYTES: {len(over)}")
    if over:
        print(f"Largest (top {args.top}):")
        for nbytes, kname, t in sorted(over, key=lambda x: -x[0])[:args.top]:
            vol = f" ({nbytes / volume_bytes:.2f} x VOLUME_BYTES)" \
                if volume_bytes else ""
            print(f"  {kname:<5} {nbytes:>11} bytes{vol} at t = {t} ns")
    if not broken and gathers == 0:
        broken.append("no device-to-host copy of GATHER_BYTES in a step: the "
                      "trace does not show the planned gather")

    print()
    if broken:
        print("TRACE CHECK: NO RESULT (" + "; ".join(broken) + ")")
        return 2
    found = (f"{len(over)} host-device copies larger than the gather in a "
             f"step" if over else
             "no host-device copy larger than the gather in a step after "
             "step 1")
    if args.gate:
        print(f"TRACE CHECK: {'FAIL' if over else 'PASS'} ({found})")
        return 1 if over else 0
    print(f"TRACE REPORT: {found} (report only, not a gate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
