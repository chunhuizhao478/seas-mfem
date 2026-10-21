#!/usr/bin/env python3
"""g4_trace_check.py — GPU plan G4: the host-device copies in an Nsight
Systems trace of seas_bp5_full, placed in the driver's NVTX ranges.

    python3 jobs/vista/g4_trace_check.py REPORT.sqlite --fault-dofs N [--gate]

REPORT.sqlite is `nsys export --type sqlite` of a CUDA build of
seas_bp5_full run with --stage-timers (which turns its NVTX ranges on).
N is the rank's fault DOF count ("Local fault DOFs: N" in the run log).

Each copy is placed by the CPU time of the CUDA call that issued it (joined
on the correlation id; copies without one fall back to their GPU start
time), in the first of these that holds it:

    g_copy    seas:g_copy, the planned copy of g(t) to the device (G3a)
    stage     seas:stage_total, one SEASQuasiDynamicOperator::Mult
    output    seas:monitor_output or seas:monitor_checkpoint
    monitor   seas:monitor_total, outside output
    step      seas:ts_solve, outside the above (PETSc between stages)
    outside   set-up and the end of the run

Gate (--gate; plan G4, "nothing state-sized crosses between host and device
inside a stage"): no host-to-device or device-to-host copy of FAULT_BYTES =
8 * N bytes or more in stage, monitor or step.  That is one value per fault
DOF, a third of the state on a 3-D fault, so it is stricter than
state-sized.  Output events bring the state to the host by design; they are
counted, not gated.  Without --gate the script only reports (multi-GPU
traces, where halo exchanges may stage through the host).

Exit: 0 pass (or report only), 1 gate failed, 2 the trace has no copies
table or no seas:* ranges, so it cannot answer.
"""

import argparse
import bisect
import collections
import sqlite3
import sys

STAGE_NAMES = {
    "g_copy": ["seas:g_copy"],
    "stage": ["seas:stage_total"],
    "output": ["seas:monitor_output", "seas:monitor_checkpoint"],
    "monitor": ["seas:monitor_total"],
    "step": ["seas:ts_solve"],
}
CONTEXT_ORDER = ["g_copy", "stage", "output", "monitor", "step", "outside"]
GATED = ("stage", "monitor", "step")
HOST_DEVICE = ("HtoD", "DtoH")
FALLBACK_KINDS = {0: "unknown", 1: "HtoD", 2: "DtoH", 3: "HtoA", 4: "AtoH",
                  5: "AtoA", 6: "AtoD", 7: "DtoA", 8: "DtoD", 9: "HtoH",
                  10: "PtoP"}


def tables(db):
    return {r[0] for r in db.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}


def columns(db, table):
    return {r[1] for r in db.execute(f"PRAGMA table_info({table})")}


def short_kind(texts):
    """HtoD, DtoH, ... from an enum row's name or label, or None."""
    t = "".join(str(x) for x in texts if x).lower()
    for ch in "-_ ":
        t = t.replace(ch, "")
    for short, keys in (("HtoD", ("htod", "hosttodevice")),
                        ("DtoH", ("dtoh", "devicetohost")),
                        ("DtoD", ("dtod", "devicetodevice")),
                        ("HtoH", ("htoh", "hosttohost")),
                        ("PtoP", ("ptop", "peertopeer"))):
        if any(k in t for k in keys):
            return short
    return None


def copy_kinds(db, have):
    """copyKind -> HtoD, DtoH, ...: CUPTI's numbering, checked against the
    trace's own enum table when it has one."""
    kinds = dict(FALLBACK_KINDS)
    if "ENUM_CUDA_MEMCPY_OPER" in have:
        for row in db.execute("SELECT * FROM ENUM_CUDA_MEMCPY_OPER"):
            short = short_kind(row[1:])
            if short:
                kinds[row[0]] = short
    return kinds


def load_ranges(db, have):
    """{range name: {thread key: sorted [(start, end)]}} for seas:* ranges."""
    if "NVTX_EVENTS" not in have:
        return {}
    cols = columns(db, "NVTX_EVENTS")
    if "textId" in cols and "StringIds" in have:
        sql = ("SELECT e.start, e.end, COALESCE(e.text, s.value), e.globalTid "
               "FROM NVTX_EVENTS e LEFT JOIN StringIds s ON e.textId = s.id "
               "WHERE e.end IS NOT NULL")
    else:
        sql = ("SELECT start, end, text, globalTid FROM NVTX_EVENTS "
               "WHERE end IS NOT NULL")
    ranges = collections.defaultdict(lambda: collections.defaultdict(list))
    for start, end, name, tid in db.execute(sql):
        if name and name.startswith("seas:"):
            ranges[name][tid].append((start, end))
    for by_tid in ranges.values():
        for spans in by_tid.values():
            spans.sort()
    return ranges


def inside(spans, t):
    """Whether time t falls in one of the sorted, disjoint spans."""
    i = bisect.bisect_right(spans, (t, float("inf"))) - 1
    return i >= 0 and spans[i][0] <= t <= spans[i][1]


def in_range(ranges, name, tid, pid_key, t):
    by_tid = ranges.get(name, {})
    if tid is not None:
        return inside(by_tid.get(tid, []), t)
    return any(inside(spans, t) for k, spans in by_tid.items()
               if (k >> 24) == pid_key)


def load_copies(db, have):
    """[(cpu time, thread or None, pid key, bytes, copyKind)] and the count
    of copies placed by their GPU time."""
    if "CUPTI_ACTIVITY_KIND_MEMCPY" not in have:
        return None, 0
    if "CUPTI_ACTIVITY_KIND_RUNTIME" in have:
        sql = ("SELECT m.start, m.globalPid, m.bytes, m.copyKind, r.start, "
               "r.globalTid FROM CUPTI_ACTIVITY_KIND_MEMCPY m "
               "LEFT JOIN CUPTI_ACTIVITY_KIND_RUNTIME r "
               "ON r.correlationId = m.correlationId "
               "AND (r.globalTid >> 24) = (m.globalPid >> 24)")
    else:
        sql = ("SELECT start, globalPid, bytes, copyKind, NULL, NULL "
               "FROM CUPTI_ACTIVITY_KIND_MEMCPY")
    copies, unmatched = [], 0
    for gpu_t, gpid, nbytes, kind, cpu_t, tid in db.execute(sql):
        if cpu_t is None:
            unmatched += 1
            copies.append((gpu_t, None, gpid >> 24, nbytes, kind))
        else:
            copies.append((cpu_t, tid, gpid >> 24, nbytes, kind))
    return copies, unmatched


def context_of(ranges, tid, pid_key, t):
    for ctx in CONTEXT_ORDER[:-1]:
        if any(in_range(ranges, n, tid, pid_key, t)
               for n in STAGE_NAMES[ctx]):
            return ctx
    return "outside"


def count_ranges(ranges, name):
    return sum(len(spans) for spans in ranges.get(name, {}).values())


def mb(nbytes):
    return nbytes / 1.0e6


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sqlite")
    ap.add_argument("--fault-dofs", type=int, required=True,
                    help="the rank's fault DOF count from the run log")
    ap.add_argument("--gate", action="store_true",
                    help="exit 1 on a fault-sized copy in a gated context")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()

    fault_bytes = 8 * args.fault_dofs
    db = sqlite3.connect(f"file:{args.sqlite}?mode=ro", uri=True)
    have = tables(db)
    kinds = copy_kinds(db, have)
    ranges = load_ranges(db, have)
    copies, unmatched = load_copies(db, have)

    print(f"Trace: {args.sqlite}")
    print(f"Fault DOFs {args.fault_dofs}: FAULT_BYTES = {fault_bytes} "
          f"({mb(fault_bytes):.3f} MB)")
    n_stage = count_ranges(ranges, "seas:stage_total")
    n_monitor = count_ranges(ranges, "seas:monitor_total")
    n_output = count_ranges(ranges, "seas:monitor_output")
    n_ts = count_ranges(ranges, "seas:ts_solve")
    print(f"NVTX ranges: ts_solve {n_ts}, stage_total {n_stage}, "
          f"monitor_total {n_monitor}, monitor_output {n_output}, "
          f"g_copy {count_ranges(ranges, 'seas:g_copy')}")

    broken = []  # the trace cannot answer the question: exit 2 either way
    if copies is None:
        broken.append("no CUDA memcpy table in the trace (was it run "
                        "with --trace=cuda?)")
        copies = []
    if n_stage == 0 or n_ts == 0:
        broken.append("no seas:stage_total or seas:ts_solve NVTX ranges "
                        "(CUDA build with nvtx3 headers, run with "
                        "--stage-timers and nsys --trace=nvtx?)")
    if unmatched:
        print(f"Note: {unmatched} of {len(copies)} copies had no matching "
              f"CUDA call; placed by their GPU start time")

    stats = collections.defaultdict(lambda: [0, 0, 0])
    big = collections.defaultdict(list)
    output_spans_with_big = set()
    for t, tid, pid_key, nbytes, kind in copies:
        ctx = context_of(ranges, tid, pid_key, t)
        kname = kinds.get(kind, str(kind))
        s = stats[(ctx, kname)]
        s[0] += 1
        s[1] += nbytes
        s[2] = max(s[2], nbytes)
        if kname in HOST_DEVICE and nbytes >= fault_bytes:
            big[ctx].append((nbytes, kname, t))
            if ctx == "output":
                spans = ranges["seas:monitor_output"].get(tid, []) \
                    if tid is not None else []
                i = bisect.bisect_right(spans, (t, float("inf"))) - 1
                output_spans_with_big.add(i)

    print()
    print(f"{'context':<9} {'kind':<6} {'copies':>8} {'total MB':>11} "
          f"{'max bytes':>11} {'per stage':>10}")
    for ctx in CONTEXT_ORDER:
        for (c, kname), (n, tot, mx) in sorted(stats.items()):
            if c != ctx:
                continue
            per = f"{n / n_stage:.2f}" if ctx in ("stage", "g_copy") \
                and n_stage else ""
            print(f"{ctx:<9} {kname:<6} {n:>8} {mb(tot):>11.3f} {mx:>11} "
                  f"{per:>10}")

    print()
    print(f"Host-device copies of FAULT_BYTES or more, by context:")
    for ctx in CONTEXT_ORDER:
        print(f"  {ctx:<9} {len(big[ctx])}")
    if n_output:
        print(f"Output ranges that brought data to the host: "
              f"{len(output_spans_with_big)} of {n_output}")

    gated = [(ctx, b) for ctx in GATED for b in big[ctx]]
    if gated:
        print()
        print(f"Largest fault-sized copies in gated contexts "
              f"(top {args.top}):")
        for ctx, (nbytes, kname, t) in sorted(
                gated, key=lambda x: -x[1][0])[:args.top]:
            print(f"  {ctx:<8} {kname:<5} {nbytes:>11} bytes at "
                  f"t = {t} ns")

    print()
    if broken:
        print("TRACE CHECK: NO RESULT (" + "; ".join(broken) + ")")
        return 2
    found = (f"{len(gated)} fault-sized host-device copies in stage, "
             f"monitor or step" if gated else
             "no fault-sized host-device copy in a stage, in the monitor "
             "outside output, or between stages")
    if args.gate:
        print(f"TRACE CHECK: {'FAIL' if gated else 'PASS'} ({found})")
        return 1 if gated else 0
    print(f"TRACE REPORT: {found} (report only, not a gate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
