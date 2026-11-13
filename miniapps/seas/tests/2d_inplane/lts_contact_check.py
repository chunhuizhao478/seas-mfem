#!/usr/bin/env python3
# tests/2d_inplane/lts_contact_check.py — make test-lts-contact:
# [cycle.opening] mode = "contact" under [cycle.lts] rate2 (LTS with PML design
# rev1; Chunhui Zhao, 2026-11-12: "implement local time stepping and test and
# integrate into future runs for speedup").  The contact rule is part of the
# fault closure and every fault tet sits in the fault cluster, so no solver
# code changes: this is the test that the fault opens and stays the same under
# LTS.  It stands alone: no band, no quench.
#
#   python3 lts_contact_check.py decks DECK TMP
#       From the D9f bend deck of the simple cases (arc, kink 8 degrees over
#       2.5 km, h = 250 m, averaged vertex normals), without its
#       [cycle.absorbing] table, with [cycle.opening] mode = "contact" and
#       numbers chosen only so that ONE event opens the releasing side of the
#       bend (they are not a physical model): sigma_n' 10 MPa (the floor),
#       f0 0.95, b 0.06, event started at t = 0 by a V_init patch of 0.1 m/s:
#       TMP/glob.toml (global stepping), TMP/same.toml (rate2 with the test
#       key same_step = true) and TMP/lts.toml (rate2).
#   python3 lts_contact_check.py check TMP PREFIX
#       After the three ran the same physical time into TMP/glob, TMP/same
#       and TMP/lts:
#        1. global stepping ends with open DOFs and a gap above 1 mm;
#        2. same_step against global stepping, the same scheme: the same open
#           DOFs, fault stresses within 1e-8 MPa and slip rate within 1e-9
#           m/s at every fault snapshot;
#        3. rate2 against global stepping, as two step sizes do: the same
#           number of open DOFs at the end, the largest gap within 2 %, the
#           slip within 0.2 % of the largest slip and the stresses within
#           0.1 MPa at the last fault snapshot.
# Exit 0 when every check passes.

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                                      # noqa: E402
from lts_pml_band_check import replace_once             # noqa: E402
import lts_avg_normals_check as avg                     # noqa: E402

failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def decks(deck, tmp):
    avg.decks(deck, tmp)
    for k in ("glob", "same", "lts"):
        p = os.path.join(tmp, k + ".toml")
        s = open(p).read()
        s = replace_once(s, "sigma_n_pa    = 25.0e6", "sigma_n_pa    = 10.0e6", "[stress]")
        s = replace_once(s, "sigma_n_default = 25.0e6", "sigma_n_default = 10.0e6", "sigma_n_default")
        s = replace_once(s, "f_0_default     = 0.6", "f_0_default     = 0.95", "f_0_default")
        s = replace_once(s, "b_default       = 0.03", "b_default       = 0.06", "b_default")
        s = replace_once(s, "inside  = 2.0e-2", "inside  = 1.0e-1", "V_init")
        s = replace_once(s, "check_every = 10000", "check_every = 0", "check_every")
        m = re.search(r"^\[cycle\.switch\]", s, re.M)
        if not m:
            sys.exit("lts_contact_check.py: no [cycle.switch] table in " + deck)
        s = s[:m.start()] + '[cycle.opening]\nmode        = "contact"\nfloor_taper = true\n\n' + s[m.start():]
        open(p, "w").write(s)


def opening(tmp, run):
    t = open(os.path.join(tmp, run + ".log")).read()
    n = int(re.search(r"open DOFs at the last step: (\d+)", t).group(1))
    g = float(re.search(r"largest gap after the settle (\S+) m", t).group(1))
    return n, g


def run_check(tmp, prefix):
    tg, g = avg.fault(tmp, "glob", prefix)
    ts, s = avg.fault(tmp, "same", prefix)
    tl, l = avg.fault(tmp, "lts", prefix)
    ng, gg = opening(tmp, "glob")
    ns, gs = opening(tmp, "same")
    nl, gl = opening(tmp, "lts")
    check(ng > 0 and gg > 1e-3, "global stepping ends with %d open DOFs, largest gap %.4g m (> 1e-3)" % (ng, gg))

    n = min(len(tg), len(ts))
    dsig = max(float(np.abs(s[k][:n] - g[k][:n]).max()) for k in ("normal_stress", "shear_stress_strike"))
    dv = float(np.abs(s["slip_rate_strike"][:n] - g["slip_rate_strike"][:n]).max())
    check(len(tg) == len(ts) and np.array_equal(tg, ts) and ns == ng and dsig <= 1e-8 and dv <= 1e-9,
          "same_step against global stepping: %d open DOFs (%d), %d fault snapshots at the same times, "
          "stresses within %.3g MPa (<= 1e-8), slip rate within %.3g m/s (<= 1e-9)" % (ns, ng, n, dsig, dv))

    smax = float(np.abs(g["slip_strike"][-1]).max())
    dslip = float(np.abs(l["slip_strike"][-1] - g["slip_strike"][-1]).max())
    dsl = max(float(np.abs(l[k][-1] - g[k][-1]).max()) for k in ("normal_stress", "shear_stress_strike"))
    check(nl == ng and abs(gl - gg) <= 0.02 * gg,
          "rate2: %d open DOFs at the end (%d), largest gap %.4g m against %.4g (within 2 %%)" % (nl, ng, gl, gg))
    check(abs(tl[-1] - tg[-1]) <= 1e-9 * tg[-1] + 1e-9 and dslip <= 2e-3 * smax and dsl <= 0.1,
          "rate2 against global stepping at the end (t = %.2f s): slip within %.3g m (<= 0.2 %% of "
          "%.3g), stresses within %.3g MPa (<= 0.1)" % (tl[-1] - tl[0], dslip, smax, dsl))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        run_check(sys.argv[2], sys.argv[3])
        sys.exit(1 if failed else 0)
    else:
        sys.exit("usage: lts_contact_check.py decks DECK TMP | check TMP PREFIX")
