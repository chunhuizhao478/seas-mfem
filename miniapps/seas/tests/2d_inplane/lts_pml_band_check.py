#!/usr/bin/env python3
# tests/2d_inplane/lts_pml_band_check.py — make test-lts-pml-band: the
# [cycle.pml] band under [cycle.lts] rate2 (LTS with PML design rev1, piece
# P-b; Chunhui Zhao, 2026-11-12: "implement local time stepping and test and
# integrate into future runs for speedup").  The band's elements sit in one
# cluster c_P, the band's memory steps once per step of c_P, and the band
# elements' steps carry 1 / sqrt(kappa_e) ([cycle.pml] dt_rule = "element").
#
#   python3 lts_pml_band_check.py decks DECK TMP
#       From the straight-fault band deck (simple_cases D0): an event started
#       at t = 0 by a V_init patch of 0.02 m/s (|x| <= 3 km, linear to 1e-9
#       m/s at 11 km), dt_rule = "element", volume snapshots of the whole
#       mesh every 5 s: TMP/glob.toml (global stepping), TMP/same.toml
#       ([cycle.lts] rate2 with the test key same_step = true: the LTS run's
#       clusters and split, every cluster at dt0) and TMP/lts.toml (rate2).
#   python3 lts_pml_band_check.py check TMP PREFIX
#       After the three ran the same physical time into TMP/glob, TMP/same
#       and TMP/lts (logs TMP/<run>.log):
#        1. the logs name the band's cluster, and the three runs have the
#           same dt0 (the element rule: the band does not shrink the step);
#        2. the waves are inside the band at the end (|v| there above 1e-4
#           m/s in the last volume snapshot);
#        3. same_step against global stepping, the same scheme: fault
#           stresses within 1e-8 MPa and slip rate within 1e-9 m/s at every
#           fault snapshot, the velocity of the whole mesh within 1e-10 m/s
#           at every volume snapshot;
#        4. rate2 against global stepping, as two step sizes do: at the last
#           fault snapshot the slip within 0.1 % of the largest slip and the
#           stresses within 0.05 MPa; the peak slip rate of the phase within
#           0.1 %; the band's velocity at the end within 5 % of its largest.
# Exit 0 when every check passes.

import os
import re
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
        sys.exit("lts_pml_band_check.py: the deck's %s text occurs %d times (want 1)" % (what, n))
    return text.replace(old, new)


def decks(deck, tmp):
    s = open(deck).read()
    s = replace_once(s,
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 0.0\ninside  = 1.0e-9\n",
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 3.0e3\ntaper   = 8.0e3\ninside  = 2.0e-2\n", "V_init")
    key = 'wall_speed_factor = "auto"\n'
    s = replace_once(s, key, key + 'dt_rule = "element"\n', "wall_speed_factor")
    s = replace_once(s, "box_min = [-30.0e3, -30.0e3]", "box_min = [-241.0e3, -241.0e3]", "box_min")
    s = replace_once(s, "box_max = [30.0e3, 30.0e3]", "box_max = [241.0e3, 241.0e3]", "box_max")
    s = replace_once(s, "dt_dyn  = 1.0 ", "dt_dyn  = 5.0 ", "volume dt_dyn")
    qd = "# The QD stages through"
    lts = '[cycle.lts]\nmode = "rate2"\n'
    out = {"glob": s,
           "same": replace_once(s, qd, lts + "same_step = true\n\n" + qd, "[cycle.qd] note"),
           "lts": replace_once(s, qd, lts + "\n" + qd, "[cycle.qd] note")}
    os.makedirs(tmp, exist_ok=True)
    for k, t in out.items():
        open(os.path.join(tmp, k + ".toml"), "w").write(t)


def log(tmp, run):
    t = open(os.path.join(tmp, run + ".log")).read()
    dt0 = float(re.search(r"^\s+dt = (\S+) s \(cfl", t, re.M).group(1))
    peak = float(re.search(r"max V: .* over the phase (\S+) m/s", t).group(1))
    m = re.search(r"LTS: rate2, (\d+) clusters .* dt0 = (\S+) s;.*PML band \((\d+) elements with a "
                  r"vertex in a band, kappa up to (\S+)\) in cluster (\d+)", t)
    return dt0, peak, m


def fault(tmp, run, prefix):
    f = h5py.File(os.path.join(tmp, run, prefix + "_fault_snapshots.vtkhdf"), "r")
    t = f["VTKHDF/Steps/Values"][...]
    pts = f["VTKHDF/Points"][...]
    npts = pts.shape[0]
    order = np.lexsort((pts[:, 1], pts[:, 0]))   # the DOF order depends on the mesh cut
    dyn = f["fd_cycle/regime"][...] == 1
    d = {}
    for k in ("normal_stress", "shear_stress_strike", "slip_rate_strike", "slip_strike"):
        d[k] = f["VTKHDF/PointData/" + k][...].reshape(len(t), npts)[:, order][dyn]
    return t[dyn], d


def volume(tmp, run, prefix):
    f = h5py.File(os.path.join(tmp, run, prefix + "_volume_snapshots.vtkhdf"), "r")
    t = f["VTKHDF/Steps/Values"][...]
    pts = f["VTKHDF/Points"][...]
    n = pts.shape[0] // len(t) if pts.shape[0] % len(t) == 0 and \
        f["VTKHDF/PointData/velocity"].shape[0] == pts.shape[0] else pts.shape[0]
    v = f["VTKHDF/PointData/velocity"][...]
    nt = v.shape[0] // n
    v = v.reshape(nt, n, 3)
    p = pts[:n]
    # the points sorted by position (the order depends on the mesh cut); a DG
    # field has several values at a vertex, so the values are sorted within
    # equal positions too
    dyn = f["fd_cycle/regime"][...][:nt] == 1
    key = np.lexsort((np.round(p[:, 1], 3), np.round(p[:, 0], 3)))
    band = np.abs(p[key, :2]).max(axis=1) > 9.0e4   # the band starts 90 km from the centre
    return t[:nt][dyn], v[dyn][:, key, :], band, p[key]


def vol_diff(a, b, pa):
    # largest difference of |v| per vertex position (several DG values share
    # a position; their order within it follows the mesh cut)
    na = np.linalg.norm(a, axis=2)
    nb = np.linalg.norm(b, axis=2)
    xy = np.round(pa[:, :2], 3)
    _, inv = np.unique(xy, axis=0, return_inverse=True)
    inv = inv.ravel()
    m = inv.max() + 1
    worst = 0.0
    for i in range(na.shape[0]):
        sa = np.zeros(m)
        sb = np.zeros(m)
        np.maximum.at(sa, inv, na[i])
        np.maximum.at(sb, inv, nb[i])
        worst = max(worst, float(np.abs(sa - sb).max()))
    return worst


def run_check(tmp, prefix):
    dg, pg, _ = log(tmp, "glob")
    ds, ps, ms = log(tmp, "same")
    dl, pl, ml = log(tmp, "lts")
    ok = ms is not None and ml is not None
    check(ok, "the LTS logs name the band's cluster" + (
        " (%s clusters, %s band elements, kappa up to %s, cluster %s)" % (
            ml.group(1), ml.group(3), ml.group(4), ml.group(5)) if ok else ""))
    if not ok:
        return
    check(dg == ds == dl == float(ml.group(2)),
          "the same dt0 in the three runs (%.6g s: the band does not shrink the step)" % dg)

    tg, g = fault(tmp, "glob", prefix)
    ts, s = fault(tmp, "same", prefix)
    tl, l = fault(tmp, "lts", prefix)
    vtg, vg, band, pg_ = volume(tmp, "glob", prefix)
    vts, vs, _, ps_ = volume(tmp, "same", prefix)
    vtl, vl, _, pl_ = volume(tmp, "lts", prefix)
    vband = float(np.linalg.norm(vg[-1][band], axis=1).max())
    check(vband > 1e-4, "the waves are inside the band at the end: |v| up to %.3g m/s there "
          "(t = %.2f s)" % (vband, tg[-1] - tg[0]))

    n = min(len(tg), len(ts))
    same_t = len(tg) == len(ts) and np.array_equal(tg, ts)
    dsig = max(float(np.abs(s[k][:n] - g[k][:n]).max()) for k in ("normal_stress", "shear_stress_strike"))
    dv = float(np.abs(s["slip_rate_strike"][:n] - g["slip_rate_strike"][:n]).max())
    check(same_t and dsig <= 1e-8 and dv <= 1e-9,
          "same_step against global stepping: %d fault snapshots at the same times, stresses "
          "within %.3g MPa (<= 1e-8), slip rate within %.3g m/s (<= 1e-9)" % (n, dsig, dv))
    nv = min(vg.shape[0], vs.shape[0])
    dvol = vol_diff(vg[:nv], vs[:nv], pg_)
    check(vg.shape == vs.shape and dvol <= 1e-10,
          "same_step against global stepping: velocity of the whole mesh within %.3g m/s "
          "(<= 1e-10) at %d volume snapshots" % (dvol, nv))

    smax = float(np.abs(g["slip_strike"][-1]).max())
    dslip = float(np.abs(l["slip_strike"][-1] - g["slip_strike"][-1]).max())
    dsl = max(float(np.abs(l[k][-1] - g[k][-1]).max()) for k in ("normal_stress", "shear_stress_strike"))
    check(abs(tl[-1] - tg[-1]) <= 1e-9 * tg[-1] + 1e-9 and dslip <= 1e-3 * smax and dsl <= 0.05,
          "rate2 against global stepping at the end (t = %.2f s): slip within %.3g m (<= 0.1 %% of "
          "%.3g), stresses within %.3g MPa (<= 0.05)" % (tl[-1] - tl[0], dslip, smax, dsl))
    check(abs(pl - pg) <= 1e-3 * pg, "rate2: peak slip rate %.6g m/s against %.6g (within 0.1 %%)" % (pl, pg))
    db = vol_diff(vg[-1:][:, band], vl[-1:][:, band], pg_[band])
    check(abs(vtl[-1] - vtg[-1]) <= 1e-9 * vtg[-1] + 1e-9 and db <= 0.05 * vband,
          "rate2: the band's velocity at the end within %.3g m/s (<= 5 %% of %.3g)" % (db, vband))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit(__doc__ or "usage: lts_pml_band_check.py decks DECK TMP | check TMP PREFIX")
