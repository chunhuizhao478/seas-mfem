#!/usr/bin/env python3
# tests/2d_inplane/lts_static_split_check.py — make test-lts-static-split: the
# rest test R1 of local time stepping in a cycle (LTS with PML design rev1,
# Section 4.1, piece P-a; Chunhui Zhao, 2026-11-11: "let's also start working
# on LTS").  [cycle.lts] static_split = true (NitscheLtsRun::SetStaticSplit)
# steps the change from the static state handed over at the switch, so a
# loaded fault at rest stays at rest: without it the rows on a cluster seam
# ripple in proportion to the static stress.
#
#   python3 lts_static_split_check.py decks DECK TMP
#       From the straight-fault band deck (simple_cases D0), without its
#       [cycle.pml] table (reflecting walls; the LTS run refuses the band
#       today): TMP/prep_glob.toml and TMP/prep_lts.toml (the QD loading run,
#       global stepping and [cycle.lts] rate2: the mesh cut differs), and the
#       restart decks with V_on = 1e-12 m/s, below the loaded fault's max V,
#       so the dynamic regime starts at once: TMP/glob.toml (global stepping),
#       TMP/lts.toml (rate2) and TMP/split.toml (rate2, static_split = true).
#   python3 lts_static_split_check.py check TMP PREFIX
#       After the three restarts ran into TMP/glob, TMP/lts and TMP/split for
#       the same physical time.  Over every fault snapshot of the dynamic
#       phase and every fault point:
#        1. the snapshots of the two rate2 runs are at the same times, those
#           of global stepping within the snapshot interval (a run writes at
#           its own step ends);
#        2. split against global stepping: normal and shear stress within
#           1e-7 MPa, slip rate within 1e-6 of its largest value;
#        3. control: rate2 without the split differs from global stepping by
#           more than 1e-5 MPa in the shear or the normal stress (the ripple
#           the split removes; printed).
# Exit 0 when every check passes.

import os
import sys

import h5py
import numpy as np

failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def replace_once(text, old, new, what):
    n = text.count(old)
    if n != 1:
        sys.exit("lts_static_split_check.py: the deck's %s text occurs %d times (want 1)" % (what, n))
    return text.replace(old, new)


def decks(deck, tmp):
    s = open(deck).read()
    a = s.index("# PML design plan rev2")
    b = s.index("# The QD stages through")
    if not (0 < a < b) or s.count("[cycle.pml]") != 1 or not (a < s.index("[cycle.pml]") < b):
        sys.exit("lts_static_split_check.py: no [cycle.pml] table where expected in " + deck)
    lts = '[cycle.lts]\nmode = "rate2"\n'
    prep_glob = s[:a] + s[b:]
    prep_lts = s[:a] + lts + "\n" + s[b:]

    def now(t):
        t = replace_once(t, "V_on               = 1.0e-2", "V_on               = 1.0e-12", "V_on")
        return replace_once(t, "V_off              = 1.0e-3", "V_off              = 1.0e-13", "V_off")

    out = {"prep_glob": prep_glob, "prep_lts": prep_lts, "glob": now(prep_glob),
           "lts": now(prep_lts),
           "split": now(replace_once(prep_lts, lts, lts + "static_split = true\n", "[cycle.lts]"))}
    os.makedirs(tmp, exist_ok=True)
    for k, t in out.items():
        open(os.path.join(tmp, k + ".toml"), "w").write(t)


def load(tmp, run, prefix):
    f = h5py.File(os.path.join(tmp, run, prefix + "_fault_snapshots.vtkhdf"), "r")
    regime = f["fd_cycle/regime"][...]
    t = f["VTKHDF/Steps/Values"][...]
    npts = f["VTKHDF/Points"].shape[0]
    pts = f["VTKHDF/Points"][...]
    # the points sorted by position: the DOF order depends on the mesh cut
    order = np.lexsort((pts[:, 1], pts[:, 0]))
    d = {}
    for k in ("normal_stress", "shear_stress_strike", "slip_rate_strike"):
        x = f["VTKHDF/PointData/" + k][...].reshape(len(t), npts)
        d[k] = x[:, order]
    dyn = regime == regime.max()   # 0 QD, 1 dynamic
    return t[dyn], {k: x[dyn] for k, x in d.items()}, int(dyn.sum())


def run_check(tmp, prefix):
    tg, g, ng = load(tmp, "glob", prefix)
    tl, l, nl = load(tmp, "lts", prefix)
    ts, s, ns = load(tmp, "split", prefix)
    n = min(ng, nl, ns)
    check(n >= 20 and abs(ng - nl) <= 1 and abs(ng - ns) <= 1,
          "%d, %d and %d dynamic fault snapshots (global, rate2, rate2 + split)" % (ng, nl, ns))
    # a run writes a snapshot at its first step end past the due time: global
    # stepping at a step, rate2 at a sync point (8 steps)
    same_t = np.array_equal(tl[:n], ts[:n]) and np.allclose(tg[:n], tl[:n], rtol=0, atol=0.05)
    check(same_t, "the snapshots are at the same times (rate2 runs equal, global within the 0.05 s "
          "snapshot interval; %d compared, %.3f s of the phase)" % (n, tg[n - 1] - tg[0]))

    def diff(a, k):
        return float(np.abs(a[k][:n] - g[k][:n]).max())

    ds = max(diff(s, "normal_stress"), diff(s, "shear_stress_strike"))
    dl = max(diff(l, "normal_stress"), diff(l, "shear_stress_strike"))
    vmax = float(np.abs(g["slip_rate_strike"][:n]).max())
    dv = diff(s, "slip_rate_strike")
    check(ds <= 1e-7, "static split against global stepping: fault stresses within %.3g MPa (<= 1e-7)" % ds)
    check(dv <= 1e-6 * vmax, "static split against global stepping: slip rate within %.3g m/s "
          "(<= 1e-6 of %.3g)" % (dv, vmax))
    check(dl > 1e-5, "control, rate2 without the split: fault stresses off by %.3g MPa (> 1e-5), "
          "%.0f times the split's" % (dl, dl / max(ds, 1e-300)))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit(__doc__ or "usage: lts_static_split_check.py decks DECK TMP | check TMP PREFIX")
