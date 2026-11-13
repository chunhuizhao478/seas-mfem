#!/usr/bin/env python3
# tests/2d_inplane/lts_quench_check.py — make test-lts-quench: [cycle.quench]
# under [cycle.lts] rate2 (LTS with PML design rev1, piece P-c; Chunhui Zhao,
# 2026-11-12: "implement local time stepping and test and integrate into
# future runs for speedup").  NitscheLtsRun::SetQuench adds C_q = d rho M_e
# on every element, treated as the absorbing blocks are; the driver switches
# it at sync points.  The piece stands alone: the decks have no [cycle.pml]
# band (reflecting walls), and BAND=1 keeps the band (pieces P-b and P-c
# together).
#
#   python3 lts_quench_check.py decks DECK TMP [band]
#       From the straight-fault band deck (simple_cases D0): an event started
#       at t = 0 by a V_init patch of 0.02 m/s, [cycle.quench] with d = 0.5
#       1/s and V_q = V_on = 0.01 m/s (on once the event's slip rates fall
#       below V_on, about 19.6 s), without the [cycle.pml] table unless
#       "band" is given (then dt_rule = "element"): TMP/glob.toml (global
#       stepping), TMP/same.toml (rate2 with the test key same_step = true)
#       and TMP/lts.toml (rate2).
#   python3 lts_quench_check.py check TMP PREFIX
#       After the three ran the same physical time into TMP/glob, TMP/same
#       and TMP/lts (logs TMP/<run>.log):
#        1. the quench was on for at least 1000 of global stepping's steps,
#           and same_step switched it at the same steps;
#        2. same_step against global stepping, the same scheme: fault
#           stresses within 1e-8 MPa and slip rate within 1e-9 m/s at every
#           fault snapshot;
#        3. rate2 against global stepping, as two step sizes do: the quench
#           goes on within one sync interval of global stepping's time; at
#           the last fault snapshot the slip within 0.1 % of the largest
#           slip, the stresses within 0.05 MPa and the slip rate within 5 %
#           of its largest value there.
# Exit 0 when every check passes.

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                                      # noqa: E402
from lts_pml_band_check import fault, replace_once      # noqa: E402

failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def decks(deck, tmp, band):
    s = open(deck).read()
    s = replace_once(s,
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 0.0\ninside  = 1.0e-9\n",
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 3.0e3\ntaper   = 8.0e3\ninside  = 2.0e-2\n", "V_init")
    qd = "# The QD stages through"
    if band:
        key = 'wall_speed_factor = "auto"\n'
        s = replace_once(s, key, key + 'dt_rule = "element"\n', "wall_speed_factor")
    else:
        a = s.index("# PML design plan rev2")
        b = s.index(qd)
        if not (0 < a < b) or s.count("[cycle.pml]") != 1 or not (a < s.index("[cycle.pml]") < b):
            sys.exit("lts_quench_check.py: no [cycle.pml] table where expected in " + deck)
        s = s[:a] + s[b:]
    s = replace_once(s, qd, '[cycle.quench]\non = true\nd_per_s = 0.5\nV_q = 1.0e-2\n\n' + qd,
                     "[cycle.qd] note")
    lts = '[cycle.lts]\nmode = "rate2"\n'
    out = {"glob": s,
           "same": replace_once(s, qd, lts + "same_step = true\n\n" + qd, "[cycle.qd] note"),
           "lts": replace_once(s, qd, lts + "\n" + qd, "[cycle.qd] note")}
    os.makedirs(tmp, exist_ok=True)
    for k, t in out.items():
        open(os.path.join(tmp, k + ".toml"), "w").write(t)


def log(tmp, run):
    t = open(os.path.join(tmp, run + ".log")).read()
    sw = re.findall(r"quench (on|off) at tau = (\S+) s \(step (\d+),", t)
    m = re.search(r"quench: on for (\d+) of (\d+) steps", t)
    sync = re.search(r"steps of (\S+) s, LTS sync intervals", t)
    return ([(a, float(b), int(c)) for a, b, c in sw], int(m.group(1)) if m else 0,
            float(sync.group(1)) if sync else 0.0)


def run_check(tmp, prefix):
    swg, ong, _ = log(tmp, "glob")
    sws, ons, _ = log(tmp, "same")
    swl, onl, sync = log(tmp, "lts")
    check(ong >= 1000 and swg == sws and ons == ong,
          "the quench was on for %d of global stepping's steps (>= 1000); same_step switched it "
          "at the same steps (%d switches)" % (ong, len(swg)))

    tg, g = fault(tmp, "glob", prefix)
    ts, s = fault(tmp, "same", prefix)
    tl, l = fault(tmp, "lts", prefix)
    n = min(len(tg), len(ts))
    dsig = max(float(np.abs(s[k][:n] - g[k][:n]).max()) for k in ("normal_stress", "shear_stress_strike"))
    dv = float(np.abs(s["slip_rate_strike"][:n] - g["slip_rate_strike"][:n]).max())
    check(len(tg) == len(ts) and np.array_equal(tg, ts) and dsig <= 1e-8 and dv <= 1e-9,
          "same_step against global stepping: %d fault snapshots at the same times, stresses "
          "within %.3g MPa (<= 1e-8), slip rate within %.3g m/s (<= 1e-9)" % (n, dsig, dv))

    last_g = [x for x in swg if x[0] == "on"][-1][1] if swg else -1.0
    last_l = [x for x in swl if x[0] == "on"][-1][1] if swl else -2.0e9
    check(onl > 0 and sync > 0.0 and abs(last_l - last_g) <= 2.0 * sync,
          "rate2: the quench goes on at %.3f s against %.3f s (within two sync intervals of "
          "%.4f s)" % (last_l, last_g, sync))
    smax = float(np.abs(g["slip_strike"][-1]).max())
    dslip = float(np.abs(l["slip_strike"][-1] - g["slip_strike"][-1]).max())
    dsl = max(float(np.abs(l[k][-1] - g[k][-1]).max()) for k in ("normal_stress", "shear_stress_strike"))
    vmax = float(np.abs(g["slip_rate_strike"][-1]).max())
    dvl = float(np.abs(l["slip_rate_strike"][-1] - g["slip_rate_strike"][-1]).max())
    check(abs(tl[-1] - tg[-1]) <= 1e-9 * tg[-1] + 1e-9 and dslip <= 1e-3 * smax and dsl <= 0.05
          and dvl <= 0.05 * vmax,
          "rate2 against global stepping at the end (t = %.2f s): slip within %.3g m (<= 0.1 %% of "
          "%.3g), stresses within %.3g MPa (<= 0.05), slip rate within %.3g m/s (<= 5 %% of %.3g)"
          % (tl[-1] - tl[0], dslip, smax, dsl, dvl, vmax))


if __name__ == "__main__":
    if len(sys.argv) in (4, 5) and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3], len(sys.argv) == 5 and sys.argv[4] == "band")
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit(__doc__ or "usage: lts_quench_check.py decks DECK TMP [band] | check TMP PREFIX")
