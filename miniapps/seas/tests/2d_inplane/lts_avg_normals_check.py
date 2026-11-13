#!/usr/bin/env python3
# tests/2d_inplane/lts_avg_normals_check.py — make test-lts-avg-normals:
# [fault_geometry] vertex_normals = "averaged" under [cycle.lts] rate2 (LTS
# with PML design rev1, piece P-d; Chunhui Zhao, 2026-11-12: "implement local
# time stepping and test and integrate into future runs for speedup").  The
# per-DOF frames are an input of the fault closure and every fault tet sits
# in the fault cluster, so the piece lifts the refusal and adds this test on
# a bend.  It stands alone: no band, no quench.
#
#   python3 lts_avg_normals_check.py decks DECK TMP
#       From the D9f bend deck of the simple cases (arc, kink 8 degrees over
#       2.5 km, h = 250 m, averaged vertex normals), without its
#       [cycle.absorbing] table (reflecting walls: rate2 refuses the
#       absorbing faces): an event started at t = 0 by a V_init patch of 0.02
#       m/s (|x| <= 3 km, linear to 1e-9 m/s at 11 km): TMP/glob.toml (global
#       stepping), TMP/same.toml (rate2 with the test key same_step = true)
#       and TMP/lts.toml (rate2).
#   python3 lts_avg_normals_check.py check TMP PREFIX
#       After the three ran the same physical time into TMP/glob, TMP/same
#       and TMP/lts:
#        1. the deck has averaged vertex normals and the bend changes the
#           normal stress (its range over the fault at the end is above 1
#           MPa: the event crossed the bend);
#        2. same_step against global stepping, the same scheme: fault
#           stresses within 1e-8 MPa and slip rate within 1e-9 m/s at every
#           fault snapshot;
#        3. rate2 against global stepping, as two step sizes do: at the last
#           fault snapshot the slip within 0.2 % of the largest slip and the
#           stresses within 0.1 MPa; the peak slip rate of the phase within
#           0.5 %.
# Exit 0 when every check passes.

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                                      # noqa: E402
import h5py                                             # noqa: E402
from lts_pml_band_check import replace_once             # noqa: E402

failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def fault(tmp, run, prefix):
    """Dynamic fault snapshots with the DOFs in an order that does not depend
    on the mesh cut: by position, and the two DOFs at a vertex (one per fault
    face; they differ at the bend) by the centre of their face."""
    f = h5py.File(os.path.join(tmp, run, prefix + "_fault_snapshots.vtkhdf"), "r")
    t = f["VTKHDF/Steps/Values"][...]
    pts = f["VTKHDF/Points"][...]
    npts = pts.shape[0]
    conn = f["VTKHDF/Connectivity"][...].reshape(-1, 2)
    if conn.size != npts or not np.array_equal(np.sort(conn.ravel()), np.arange(npts)):
        sys.exit("lts_avg_normals_check.py: expected two own points per fault face")
    mid = np.empty((npts, 3))
    mid[conn[:, 0]] = mid[conn[:, 1]] = 0.5 * (pts[conn[:, 0]] + pts[conn[:, 1]])
    order = np.lexsort((mid[:, 1], mid[:, 0], pts[:, 1], pts[:, 0]))
    dyn = f["fd_cycle/regime"][...] == 1
    d = {}
    for k in ("normal_stress", "shear_stress_strike", "slip_rate_strike", "slip_strike"):
        d[k] = f["VTKHDF/PointData/" + k][...].reshape(len(t), npts)[:, order][dyn]
    return t[dyn], d


def decks(deck, tmp):
    s = open(deck).read()
    if not re.search(r'^vertex_normals\s*=\s*"averaged"', s, re.M):
        sys.exit("lts_avg_normals_check.py: " + deck + " has no averaged vertex normals")
    s = replace_once(s,
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 0.0\ninside  = 1.0e-9\n",
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 3.0e3\ntaper   = 8.0e3\ninside  = 2.0e-2\n", "V_init")
    qd = "# The QD stages through"
    a = s.index("# Every outer face absorbs")
    b = s.index(qd)
    if not (0 < a < b) or s.count("[cycle.absorbing]") != 1 or not (a < s.index("[cycle.absorbing]") < b):
        sys.exit("lts_avg_normals_check.py: no [cycle.absorbing] table where expected in " + deck)
    s = s[:a] + s[b:]
    lts = '[cycle.lts]\nmode = "rate2"\n'
    out = {"glob": s,
           "same": replace_once(s, qd, lts + "same_step = true\n\n" + qd, "[cycle.qd] note"),
           "lts": replace_once(s, qd, lts + "\n" + qd, "[cycle.qd] note")}
    os.makedirs(tmp, exist_ok=True)
    for k, t in out.items():
        open(os.path.join(tmp, k + ".toml"), "w").write(t)


def peak(tmp, run):
    t = open(os.path.join(tmp, run + ".log")).read()
    return float(re.search(r"max V: .* over the phase (\S+) m/s", t).group(1))


def run_check(tmp, prefix):
    deck = open(os.path.join(os.path.dirname(tmp.rstrip("/")), "lts.toml")).read()
    tg, g = fault(tmp, "glob", prefix)
    ts, s = fault(tmp, "same", prefix)
    tl, l = fault(tmp, "lts", prefix)
    rng = float(g["normal_stress"][-1].max() - g["normal_stress"][-1].min())
    check(re.search(r'^vertex_normals\s*=\s*"averaged"', deck, re.M) is not None and rng > 1.0,
          "averaged vertex normals, and the event crossed the bend: normal stress range %.3g MPa "
          "over the fault at the end (> 1)" % rng)

    n = min(len(tg), len(ts))
    dsig = max(float(np.abs(s[k][:n] - g[k][:n]).max()) for k in ("normal_stress", "shear_stress_strike"))
    dv = float(np.abs(s["slip_rate_strike"][:n] - g["slip_rate_strike"][:n]).max())
    check(len(tg) == len(ts) and np.array_equal(tg, ts) and dsig <= 1e-8 and dv <= 1e-9,
          "same_step against global stepping: %d fault snapshots at the same times, stresses "
          "within %.3g MPa (<= 1e-8), slip rate within %.3g m/s (<= 1e-9)" % (n, dsig, dv))

    smax = float(np.abs(g["slip_strike"][-1]).max())
    dslip = float(np.abs(l["slip_strike"][-1] - g["slip_strike"][-1]).max())
    dsl = max(float(np.abs(l[k][-1] - g[k][-1]).max()) for k in ("normal_stress", "shear_stress_strike"))
    check(abs(tl[-1] - tg[-1]) <= 1e-9 * tg[-1] + 1e-9 and dslip <= 2e-3 * smax and dsl <= 0.1,
          "rate2 against global stepping at the end (t = %.2f s): slip within %.3g m (<= 0.2 %% of "
          "%.3g), stresses within %.3g MPa (<= 0.1)" % (tl[-1] - tl[0], dslip, smax, dsl))
    pg, pl = peak(tmp, "glob"), peak(tmp, "lts")
    check(abs(pl - pg) <= 5e-3 * pg, "rate2: peak slip rate %.6g m/s against %.6g (within 0.5 %%)" % (pl, pg))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit(__doc__ or "usage: lts_avg_normals_check.py decks DECK TMP | check TMP PREFIX")
