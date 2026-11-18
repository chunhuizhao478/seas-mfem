#!/usr/bin/env python3
# tests/viscoelastic/visco_cycle_check.py — make test-visco-cycle: phase V5 of
# the viscoelastic attenuation plan rev6 (the cycle wiring, section 5.9).
#
#   python3 visco_cycle_check.py decks DECK TMP
#       From the straight-fault band deck (simple_cases D0) through
#       tests/2d_inplane/lts_2d5_check.py (an event at t = 0, W = 15 km on the
#       loading plane, no PML band): TMP/off.toml (no [rheology], the
#       control), TMP/ve.toml (+ [rheology.viscoelastic] with the 2.5D width
#       memory), TMP/venow.toml (strain memory only), TMP/vetol.toml
#       (handback_tol_pa = 1, must stop), TMP/vepml.toml (ve + the deck's
#       [cycle.pml] band: band elements carry no memory, decision D4) and the
#       refusal decks TMP/r_*.toml.
#   python3 visco_cycle_check.py check TMP
#       After the runs (logs TMP/<deck>.log, exit codes TMP/<deck>.rc):
#        1. off, ve, venow, vepml end at the guard (exit 7, as the control);
#           vetol exits 9 with "attenuation_handback" in its switch log;
#        2. the dynamic step: ve's dt = off's dt x s_Y (to the printed 5
#           digits), s_Y < 1;
#        3. the hand-back report is printed (with max |eta| when the width
#           memory is on); ve on 1 and 2 ranks print the same numbers;
#        4. vepml excludes band elements (count > 0), ve none;
#        5. each refusal deck aborts with its message.
# Exit 0 when every check passes.

import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
failed = 0


def check(ok, what):
    global failed
    print(("  PASSED: " if ok else "  FAILED: ") + what)
    if not ok:
        failed += 1


VE = """
[rheology]
[rheology.viscoelastic]
enabled = true
mechanisms = 3
f_central_hz = 0.1
f_ratio = 100.0
Q_s_over_vs = 0.05
Q_p_over_Q_s = 1.5
moduli = "relaxed"
width_memory = true
in_pml_band = false
handback_tol_pa = 0.0
"""

# refusal decks: name -> (base, edit, expected message)
REFUSALS = {
    "r_fw": ("off", None, "width_memory = true needs [cycle.finite_width]"),
    "r_sw": ("ve", ("switching            = true", "switching            = false"),
             "needs [cycle].switching = true"),
    "r_probe": ("ve", ("[cycle.finite_width]",
                       "[cycle.probe]\nV_off = [1.0e-3]\nT_dwell = [10.0]\n\n[cycle.finite_width]"),
                "with [cycle.probe]"),
    "r_dp": ("ve", ("[rheology.viscoelastic]",
                    "[rheology.drucker_prager]\nenabled = true\ncohesion_pa = 1.36e6\n"
                    "bulk_friction = 0.1934\nrelaxation_time_s = 0.03\n\n"
                    "[rheology.viscoelastic]"),
             "plan section 2.14"),
    "r_odt": ("ve", ("[rheology]\n", "[rheology]\noutput_dt_s = 1.0\n"),
              "output_dt_s > 0"),
}


def replace_once(s, a, b):
    assert s.count(a) == 1, f"'{a[:40]}' occurs {s.count(a)} times"
    return s.replace(a, b)


def decks(deck, tmp):
    os.makedirs(tmp, exist_ok=True)
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "2d_inplane", "lts_2d5_check.py"),
                        "decks", deck, tmp])
    if r.returncode != 0:
        sys.exit("lts_2d5_check.py decks failed")
    src = open(deck).read()
    glob_ = open(os.path.join(tmp, "glob.toml")).read()
    off = open(os.path.join(tmp, "off.toml")).read()
    m = re.search(r"\n\[cycle\.pml\]\n(?:[^\[\n][^\n]*\n)*", src)
    assert m, "no [cycle.pml] table in " + deck
    pml = m.group(0).lstrip("\n")
    out = {
        "off": glob_,          # the control: the 2.5D deck without [rheology]
        "ve": glob_ + VE,
        "venow": glob_ + VE.replace("width_memory = true", "width_memory = false"),
        "vetol": glob_ + VE.replace("handback_tol_pa = 0.0", "handback_tol_pa = 1.0"),
        "vepml": replace_once(glob_, "[cycle.finite_width]", pml + "\n[cycle.finite_width]") + VE,
    }
    base = {"off": off + VE, "ve": out["ve"]}
    for name, (b, edit, _) in REFUSALS.items():
        s = base[b]
        if edit:
            s = replace_once(s, *edit)
        out[name] = s
    for name, s in out.items():
        with open(os.path.join(tmp, name + ".toml"), "w") as f:
            f.write(s)
    print("  decks: " + " ".join(sorted(out)) + " in " + tmp)


def log(tmp, name):
    p = os.path.join(tmp, name + ".log")
    return open(p).read() if os.path.exists(p) else ""


def rc(tmp, name):
    p = os.path.join(tmp, name + ".rc")
    return int(open(p).read().strip()) if os.path.exists(p) else None


def dt(text):
    m = re.search(r"\n\s+dt = ([0-9.eE+-]+) s \(cfl", text)
    return float(m.group(1)) if m else None


def run_check(tmp):
    for n in ("off", "ve", "ve_np2", "venow", "vepml"):
        check(rc(tmp, n) == 7, f"{n}: exit {rc(tmp, n)} (7, the guard, as the control)")
    check(rc(tmp, "vetol") == 9, f"vetol: exit {rc(tmp, 'vetol')} (9, handback_tol_pa)")
    sw = glob.glob(os.path.join(tmp, "vetol", "*switch.txt"))
    check(bool(sw) and "QD stop attenuation_handback" in open(sw[0]).read(),
          "vetol: switch log ends 'QD stop attenuation_handback'")
    check("STOP: the anelastic fault traction the hand-back would drop" in log(tmp, "vetol"),
          "vetol: STOP message")

    ve, off = log(tmp, "ve"), log(tmp, "off")
    m = re.search(r"s_Y = ([0-9.eE+-]+), the CFL step times ([0-9.eE+-]+); (\d+) band elements "
                  r"without memory", ve)
    check(m is not None, "ve: attenuation setup line")
    if m:
        sY = float(m.group(1))
        d1, d0 = dt(ve), dt(off)
        check(d1 is not None and d0 is not None and sY < 1.0
              and abs(d1 - d0 * sY) <= 5e-5 * d1,
              f"ve: dt {d1} = off dt {d0} x s_Y {sY}")
        check(abs(float(m.group(2)) - sY) <= 1e-5 * sY,
              f"ve: without a band the step factor {m.group(2)} is s_Y (to the 2.5D cut's change)")
        check(int(m.group(3)) == 0, f"ve: {m.group(3)} band elements without memory (no band)")
    mp = re.search(r"s_Y = ([0-9.eE+-]+), the CFL step times ([0-9.eE+-]+); (\d+) band elements "
                   r"without memory", log(tmp, "vepml"))
    check(mp is not None and int(mp.group(3)) > 0,
          f"vepml: {mp.group(3) if mp else '?'} band elements without memory (> 0)")
    if mp:
        # the band's own step (1 / sqrt(kappa_max)) is below the interior's
        # times s_Y, so the band sets the step and attenuation cuts nothing
        check(float(mp.group(2)) == 1.0,
              f"vepml: the band sets the step (factor {mp.group(2)}, s_Y {mp.group(1)})")
    check("No such" not in off and "[rheology" not in off, "off: no attenuation output")

    pat = r"hand-back drops max \|Q dsig\| ([0-9.eE+-]+) MPa on the fault(?:, max \|eta\| ([0-9.eE+-]+) m)?"
    h1 = re.search(pat, ve)
    h2 = re.search(pat, log(tmp, "ve_np2"))
    hn = re.search(pat, log(tmp, "venow"))
    check(h1 is not None and h1.group(2) is not None and float(h1.group(1)) > 0,
          f"ve: hand-back report {h1.group(0) if h1 else '?'}")
    check(hn is not None and hn.group(2) is None, "venow: hand-back report without eta")
    check(h1 is not None and h2 is not None and h1.groups() == h2.groups(),
          "ve: the hand-back report on 2 ranks equals 1 rank")

    for name, (_, _, msg) in REFUSALS.items():
        r = rc(tmp, name)
        check(r not in (None, 0) and msg in log(tmp, name).replace("\n", " "),
              f"{name}: aborts with '{msg}' (exit {r})")
    print(f"\n  {'FAILED' if failed else 'PASSED'}: test-visco-cycle ({failed} failed)")
    return failed


if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] == "decks":
        decks(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 3 and sys.argv[1] == "check":
        sys.exit(1 if run_check(sys.argv[2]) else 0)
    else:
        sys.exit(__doc__ if __doc__ else "usage: visco_cycle_check.py decks DECK TMP | check TMP")
