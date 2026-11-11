#!/usr/bin/env python3
# tests/2d_inplane/pml_dt_rule_check.py — make test-pml-dt-rule: [cycle.pml]
# dt_rule = "element" (cycle/dyn_phase.hpp; Chunhui Zhao, 2026-11-11: "let's
# focus on dynamic part speedup now").  The band's SIPG penalty is kappa
# times the plain one.  The default rule divides the whole mesh's CFL step by
# sqrt(kappa_max); the element rule takes min_e dt_e / sqrt(kappa_e), so the
# coarse band elements do not shrink the step of the fine fault elements.
#
#   python3 pml_dt_rule_check.py decks DECK TMP
#       From the straight-fault band deck (simple_cases D0): an event started
#       at t = 0 by a V_init patch of 0.02 m/s (|x| <= 3 km, linear to 1e-9
#       m/s at 11 km), TMP/max.toml (default rule) and TMP/element.toml
#       (dt_rule = "element").
#   python3 pml_dt_rule_check.py check TMP PREFIX
#       After both ran the same physical time into TMP/max and TMP/element
#       (logs TMP/max.log, TMP/element.log):
#        1. the logs name the rule; the default step is the element step /
#           sqrt(kappa_max) on this mesh (band elements 10 times the fault's),
#           so the element rule takes sqrt(kappa_max) times fewer steps;
#        2. the element run is stable: its peak slip rate over the phase is
#           the default run's to 1 %;
#        3. the fault snapshots agree as two runs with different steps do:
#           slip within 1 % of the largest slip, shear and normal stress
#           within 0.05 MPa, at every point and snapshot.
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
        sys.exit("pml_dt_rule_check.py: the deck's %s text occurs %d times (want 1)" % (what, n))
    return text.replace(old, new)


def decks(deck, tmp):
    s = open(deck).read()
    s = replace_once(s,
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 0.0\ninside  = 1.0e-9\n",
                     "[cycle.fields.V_init]\ncenter  = [0.0, 0.0]\naxis1   = [1.0, 0.0]\n"
                     "half1   = 3.0e3\ntaper   = 8.0e3\ninside  = 2.0e-2\n", "V_init")
    os.makedirs(tmp, exist_ok=True)
    open(os.path.join(tmp, "max.toml"), "w").write(s)
    key = 'wall_speed_factor = "auto"\n'
    open(os.path.join(tmp, "element.toml"), "w").write(
        replace_once(s, key, key + 'dt_rule = "element"\n', "wall_speed_factor"))


def log(tmp, run):
    t = open(os.path.join(tmp, run + ".log")).read()
    dt = float(re.search(r"^\s+dt = (\S+) s \(cfl", t, re.M).group(1))
    kmax = float(re.search(r"kappa_max = ([0-9.eE+-]+)", t).group(1))
    peak = float(re.search(r"max V: .* over the phase (\S+) m/s", t).group(1))
    return t, dt, kmax, peak


def snaps(tmp, run, prefix):
    f = h5py.File(os.path.join(tmp, run, prefix + "_fault_snapshots.vtkhdf"), "r")
    t = f["VTKHDF/Steps/Values"][...]
    n = f["VTKHDF/Points"].shape[0]
    return t, {k: f["VTKHDF/PointData/" + k][...].reshape(len(t), n)
               for k in ("slip_strike", "shear_stress_strike", "normal_stress")}


def run_check(tmp, prefix):
    tm, dtm, kmax, pm = log(tmp, "max")
    te, dte, _, pe = log(tmp, "element")
    check("the CFL step divided by sqrt(kappa_max" in tm and "dt_rule element" in te,
          "the logs name the rules (default: kappa_max; element)")
    r = dte / dtm
    check(abs(r - np.sqrt(kmax)) <= 1e-4 * r,   # the log prints six digits
          "step %.6g s (default) and %.6g s (element): ratio %.4f = sqrt(kappa_max = %.4f)"
          % (dtm, dte, r, kmax))
    check(abs(pe - pm) <= 0.01 * pm, "peak slip rate %.6f m/s (element) against %.6f m/s "
          "(default): %.2f %%" % (pe, pm, 100.0 * abs(pe - pm) / pm))
    t1, a = snaps(tmp, "max", prefix)
    t2, b = snaps(tmp, "element", prefix)
    n = min(len(t1), len(t2))
    check(n >= 20 and np.abs(t1[:n] - t2[:n]).max() <= 0.05,
          "%d fault snapshots at the same times (within the 0.05 s interval)" % n)
    smax = float(np.abs(a["slip_strike"][:n]).max())
    ds = float(np.abs(a["slip_strike"][:n] - b["slip_strike"][:n]).max())
    dtau = float(np.abs(a["shear_stress_strike"][:n] - b["shear_stress_strike"][:n]).max())
    dsig = float(np.abs(a["normal_stress"][:n] - b["normal_stress"][:n]).max())
    check(ds <= 0.01 * smax, "slip within %.3g m (<= 1 %% of %.3g m)" % (ds, smax))
    check(dtau <= 0.05 and dsig <= 0.05,
          "shear stress within %.3g MPa, normal stress within %.3g MPa (<= 0.05)" % (dtau, dsig))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit("usage: pml_dt_rule_check.py decks DECK TMP | check TMP PREFIX")
