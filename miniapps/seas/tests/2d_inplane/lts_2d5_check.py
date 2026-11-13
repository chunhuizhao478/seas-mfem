#!/usr/bin/env python3
# tests/2d_inplane/lts_2d5_check.py — make test-lts-2d5: [cycle.finite_width]
# under [cycle.lts] rate2 (2.5D finite width with PML, quench and LTS, plan
# rev1, piece F-c; Chunhui Zhao, 2026-11-12: "we should work on 2.5D with PML,
# quenching and everything + LTS").  NitscheLtsRun adds M_beta to each
# element's own block of K and the substrate load to the element's own
# cluster; BuildNitscheLtsAssignment applies rule C6 to each element's step.
# The piece stands alone: without "band" the decks have no [cycle.pml] band
# (reflecting walls); with it they keep the band (pieces F-a and F-c
# together).
#
#   python3 lts_2d5_check.py decks DECK TMP [band]
#       From the straight-fault band deck (simple_cases D0): an event started
#       at t = 0 by a V_init patch of 0.02 m/s, [cycle.finite_width] with W =
#       15 km (cutoff 0.087 Hz; the event still ruptures) on the loading plane,
#       without the [cycle.pml] table unless "band" is given (then dt_rule =
#       "element"): TMP/glob.toml (global stepping), TMP/same.toml (rate2 with
#       the test key same_step = true), TMP/lts.toml (rate2) and TMP/off.toml
#       (global stepping without the 2.5D term, the control).
#   python3 lts_2d5_check.py check TMP PREFIX
#       After the four ran the same physical time into TMP/glob, TMP/same,
#       TMP/lts and TMP/off (logs TMP/<run>.log):
#        1. the control: the 2.5D term changes the final slip by more than
#           1 % of the largest slip (so the checks below see the term);
#        2. same_step against global stepping, the same scheme: fault
#           stresses within 1e-8 MPa and slip rate within 1e-9 m/s at every
#           fault snapshot;
#        3. rate2 against global stepping, as two step sizes do: at the last
#           fault snapshot the slip within 0.1 % of the largest slip, the
#           stresses within 0.05 MPa and the slip rate within 5 % of its
#           largest value there.
# Exit 0 when every check passes.

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py                                             # noqa: E402
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
            sys.exit("lts_2d5_check.py: no [cycle.pml] table where expected in " + deck)
        s = s[:a] + s[b:]
    fw = '[cycle.finite_width]\nW_m = 15.0e3\nsubstrate = "loading_plane"\n\n'
    lts = '[cycle.lts]\nmode = "rate2"\n'
    out = {"off": s,
           "glob": replace_once(s, qd, fw + qd, "[cycle.qd] note"),
           "same": replace_once(s, qd, fw + lts + "same_step = true\n\n" + qd, "[cycle.qd] note"),
           "lts": replace_once(s, qd, fw + lts + "\n" + qd, "[cycle.qd] note")}
    os.makedirs(tmp, exist_ok=True)
    for k, t in out.items():
        open(os.path.join(tmp, k + ".toml"), "w").write(t)


def sorted_points(tmp, run, prefix):
    f = h5py.File(os.path.join(tmp, run, prefix + "_fault_snapshots.vtkhdf"), "r")
    pts = f["VTKHDF/Points"][...]
    return pts[np.lexsort((pts[:, 1], pts[:, 0]))]


def match_ties(pts, d, ref):
    """Two fault nodes of neighbouring faces sit at the same point, and
    fault()'s sort leaves their order to the file's, which follows the mesh
    cut (rate2 has its own cut).  Order each such pair as ref has it: the
    order with the smaller shear stress misfit over all snapshots."""
    k = "shear_stress_strike"
    n = min(d[k].shape[0], ref[k].shape[0])
    swaps = 0
    for i in range(pts.shape[0] - 1):
        if not np.array_equal(pts[i], pts[i + 1]):
            continue
        keep = np.abs(d[k][:n, i] - ref[k][:n, i]).sum() + np.abs(d[k][:n, i + 1] - ref[k][:n, i + 1]).sum()
        swap = np.abs(d[k][:n, i + 1] - ref[k][:n, i]).sum() + np.abs(d[k][:n, i] - ref[k][:n, i + 1]).sum()
        if swap < keep:
            for a in d.values():
                a[:, [i, i + 1]] = a[:, [i + 1, i]]
            swaps += 1
    return swaps


def run_check(tmp, prefix):
    tg, g = fault(tmp, "glob", prefix)
    ts, s = fault(tmp, "same", prefix)
    tl, l = fault(tmp, "lts", prefix)
    to, o = fault(tmp, "off", prefix)
    pts = sorted_points(tmp, "glob", prefix)
    for run, d in (("same", s), ("lts", l), ("off", o)):
        if not np.array_equal(sorted_points(tmp, run, prefix), pts):
            sys.exit("lts_2d5_check.py: the fault points of " + run + " are not global stepping's")
        match_ties(pts, d, g)
    smax = float(np.abs(g["slip_strike"][-1]).max())
    doff = float(np.abs(o["slip_strike"][-1] - g["slip_strike"][-1]).max())
    check(abs(to[-1] - tg[-1]) <= 1e-3 * (tg[-1] - tg[0]) and doff >= 1e-2 * smax,
          "the control: the 2.5D term changes the final slip by %.3g m (>= 1 %% of %.3g)"
          % (doff, smax))

    n = min(len(tg), len(ts))
    dsig = max(float(np.abs(s[k][:n] - g[k][:n]).max()) for k in ("normal_stress", "shear_stress_strike"))
    dv = float(np.abs(s["slip_rate_strike"][:n] - g["slip_rate_strike"][:n]).max())
    check(len(tg) == len(ts) and np.array_equal(tg, ts) and dsig <= 1e-8 and dv <= 1e-9,
          "same_step against global stepping: %d fault snapshots at the same times, stresses "
          "within %.3g MPa (<= 1e-8), slip rate within %.3g m/s (<= 1e-9)" % (n, dsig, dv))

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
        sys.exit(__doc__ or "usage: lts_2d5_check.py decks DECK TMP [band] | check TMP PREFIX")
