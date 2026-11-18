#!/usr/bin/env python3
# tests/viscoelastic/visco_lts_check.py — make test-visco-lts: phase V6 of the
# viscoelastic attenuation plan rev6 (attenuation under [cycle.lts] rate2,
# section 5.8), on the 2.5D straight-fault deck of make test-lts-2d5.
#
#   python3 visco_lts_check.py decks DECK TMP
#       tests/2d_inplane/lts_2d5_check.py's decks (an event at t = 0, W = 15
#       km on the loading plane, no PML band) with [rheology.viscoelastic]
#       (the 2.5D width memory on) added: TMP/glob.toml (global stepping),
#       TMP/same.toml (rate2, same_step = true), TMP/lts.toml (rate2); both
#       rate2 decks with the test key check_ghosts = true.  TMP/off.toml is
#       global stepping without attenuation (the control).  TMP/band.toml:
#       rate2 with the deck's [cycle.pml] band kept (dt_rule "element"),
#       a start-up check of the band and attenuation together.
#   python3 visco_lts_check.py check TMP PREFIX
#       After the runs (TMP/<run>/, logs TMP/<run>.log):
#        1. the control: attenuation changes the slip rate history by >= 1 %
#           of its peak;
#        2. same_step against global stepping (the same scheme, the owner
#           split summed in another order): fault stresses within 1e-8 MPa
#           and slip rate within 1e-9 m/s at every fault snapshot;
#        3. rate2 against global stepping (two step sizes): at the last fault
#           snapshot the slip within 0.1 % of the largest slip, the stresses
#           within 0.05 MPa and the slip rate within 5 % of its largest value;
#        4. the ghosts: every rate2 log reports its check_ghosts sync points
#           with max |ghost - owner| = 0 (the run stops at any difference);
#        5. the band run ends at the switch back (exit 0 or 7) with band
#           elements excluded and the step report.
# Exit 0 when every check passes.

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "2d_inplane"))
import numpy as np                                                 # noqa: E402
from lts_pml_band_check import fault                               # noqa: E402
from lts_2d5_check import sorted_points, match_ties                # noqa: E402
from visco_cycle_check import VE, replace_once                     # noqa: E402

failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


def run_decks(script, deck, tmp, *opts):
    import subprocess
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "2d_inplane", script),
                        "decks", deck, tmp, *opts])
    if r.returncode != 0:
        sys.exit(script + " decks failed")


def decks(deck, tmp):
    os.makedirs(tmp, exist_ok=True)
    plain = os.path.join(tmp, "plain")
    band = os.path.join(tmp, "plainband")
    run_decks("lts_2d5_check.py", deck, plain)
    run_decks("lts_2d5_check.py", deck, band, "band")
    rd = lambda d, k: open(os.path.join(d, k + ".toml")).read()
    ghosts = ('mode = "rate2"\n', 'mode = "rate2"\ncheck_ghosts = true\n')
    out = {
        "off": rd(plain, "glob"),
        "glob": rd(plain, "glob") + VE,
        "same": replace_once(rd(plain, "same"), *ghosts) + VE,
        "lts": replace_once(rd(plain, "lts"), *ghosts) + VE,
        "band": replace_once(rd(band, "lts"), *ghosts) + VE,
    }
    for k, s in out.items():
        with open(os.path.join(tmp, k + ".toml"), "w") as f:
            f.write(s)
    print("  decks: " + " ".join(sorted(out)) + " in " + tmp)


def log(tmp, name):
    p = os.path.join(tmp, name + ".log")
    return open(p).read() if os.path.exists(p) else ""


def run_check(tmp, prefix):
    tg, g = fault(tmp, "glob", prefix)
    ts, s = fault(tmp, "same", prefix)
    tl, l = fault(tmp, "lts", prefix)
    to, o = fault(tmp, "off", prefix)
    pts = sorted_points(tmp, "glob", prefix)
    for run, d in (("same", s), ("lts", l), ("off", o)):
        if not np.array_equal(sorted_points(tmp, run, prefix), pts):
            sys.exit("visco_lts_check.py: the fault points of " + run + " are not global stepping's")
        match_ties(pts, d, g)

    # the control runs with the elastic step (no s_Y cut), so its snapshots
    # fall at other times: its slip rate is interpolated to global
    # stepping's snapshot times over the common interval
    sr_g, sr_o = g["slip_rate_strike"], o["slip_rate_strike"]
    keep = tg <= to[-1]
    oi = np.stack([np.interp(tg[keep], to, sr_o[:, j]) for j in range(sr_o.shape[1])], axis=1)
    vpk = float(np.abs(sr_g[keep]).max())
    dvo = float(np.abs(oi - sr_g[keep]).max())
    smax = float(np.abs(g["slip_strike"][-1]).max())
    so = float(np.abs(o["slip_strike"][-1]).max())
    check(keep.sum() > 0.9 * len(tg) and dvo >= 1e-2 * vpk,
          "the control: attenuation changes the slip rate history by %.3g m/s (>= 1 %% of its "
          "peak %.3g, %d common snapshots); largest final slip %.4g m with it, %.4g m without"
          % (dvo, vpk, int(keep.sum()), smax, so))

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

    pat = r"check_ghosts: (\d+) sync points, max \|ghost - owner\| ([0-9.eE+-]+) Pa"
    for run in ("same", "lts", "band"):
        mm = re.search(pat, log(tmp, run))
        check(mm is not None and int(mm.group(1)) > 0 and float(mm.group(2)) == 0.0,
              "%s: ghost memories equal their owners bit for bit (%s)"
              % (run, mm.group(0) if mm else "no check_ghosts report"))

    b = log(tmp, "band")
    mb = re.search(r"(\d+) band elements without memory", b)
    rc = open(os.path.join(tmp, "band.rc")).read().strip() if os.path.exists(
        os.path.join(tmp, "band.rc")) else "?"
    check(rc in ("0", "7") and mb is not None and int(mb.group(1)) > 0
          and "attenuation: steps outside the band times s_Y" in b,
          "band: rate2 with the PML band and attenuation runs (exit %s, %s band elements "
          "without memory)" % (rc, mb.group(1) if mb else "?"))
    return failed


if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 4 and sys.argv[1] == "check":
        sys.exit(1 if run_check(sys.argv[2], sys.argv[3]) else 0)
    else:
        sys.exit("usage: visco_lts_check.py decks DECK TMP | check TMP PREFIX")
