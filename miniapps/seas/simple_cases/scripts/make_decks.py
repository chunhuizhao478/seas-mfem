#!/usr/bin/env python3
# The simple cases' decks (Chunhui Zhao, 2026-10-31: "go ahead with the simple
# cases, run them on local Mac"; plan "2D plane strain simple tests plan"):
# twelve decks made from config/fd_cycle/mapview_rs_fd.toml, written to
# simple_cases/configs/<case>.toml; a thirteenth, D5f (2026-10-31: "let's
# prepare for this case, a smooth design around this bend, and rerun"), is
# D1f with its 1.5 deg corner rounded into a 10 km arc.  Run from miniapps/seas:
#   python3 simple_cases/scripts/make_decks.py [--check]
# --check: write nothing, exit 1 if a committed deck differs from what this
# script makes (the decks are committed as made).
#
# Edits to mapview_rs_fd.toml, every case:
#   [mesh] path        simple_cases/mesh/<mesh>.msh (+-240 km box, 500 m)
#   [time] tfinal      12,000 yr of 365 days (about 3 events, see README)
#   [cycle.output]     prefix = <case>; [[cycle.output.station]] on the
#                      fault by arc length (the trace's own points, so a
#                      bent limb's stations sit on it), plus the bend set
#   [cycle.benchmark]  region_x3 = [-1e4, 1e4] (the bent limb reaches
#                      y = 1.31 km), as the bend convergence decks
#   [cycle.snapshots.fault] dt_dyn = 0.05 s (was 0.1 s; output only)
#   [cycle.snapshots.volume] dt_dyn = 1.0 s (was 0.5 s; 1 MB per snapshot,
#                      so about 0.27 GB per fully dynamic run of 3 events)
# QD only (Q*):        [cycle] switching = false
# fully dynamic (D*):  [cycle.output] handback_table = true
# floor (D1f, D2f, D5f): [cycle] sigma_n_floor_pa = 1.0e6 (UCERF run8's)
import argparse
import math
import os
import sys

BASE = "config/fd_cycle/mapview_rs_fd.toml"
YR = 365.0 * 86400.0
T_FINAL_YR = 12000.0
FAULT_HALF = 15.0e3
LC_FAULT = 500.0     # the fault points' spacing in arc length [m]

# mesh: (theta_deg, arc [m], s_bend [m]), as mesh/README and bend.geo
MESHES = {
    "straight_h500": (0.0, 0.0, 0.0),
    "sharp015_h500": (1.5, 0.0, 0.0),
    "sharp05_h500": (5.0, 0.0, 0.0),
    "arc05L10_h500": (5.0, 10000.0, 0.0),
    "sharp015vs_h500": (1.5, 0.0, 12500.0),
    "arc015L10_h500": (1.5, 10000.0, 0.0),
}

# case: (mesh, regime, floor, plan row), in the run order
CASES = [
    ("D1f_sharp015_floor", "sharp015_h500", "D", True, "D1 with the 1 MPa floor: a 1.5 deg sharp bend at x = 0"),
    ("D2f_sharp05_floor", "sharp05_h500", "D", True, "D2 with the 1 MPa floor: a 5 deg sharp bend at x = 0"),
    ("Q1_sharp015", "sharp015_h500", "Q", False, "Andy's [2]: UCERF's angle, a 1.5 deg sharp bend at x = 0"),
    ("Q0_straight", "straight_h500", "Q", False, "Andy's [1]: the straight control"),
    ("Q2_sharp05", "sharp05_h500", "Q", False, "Andy's [2]: a 5 deg sharp bend at x = 0"),
    ("Q3_arc05L10", "arc05L10_h500", "Q", False, "Andy's [2]: 5 deg over a 10 km arc centred at x = 0"),
    ("Q4_sharp015vs", "sharp015vs_h500", "Q", False, "addition: a 1.5 deg sharp bend at x = 12.5 km, in the creeping part"),
    ("D0_straight", "straight_h500", "D", False, "Andy's [3]: the straight control"),
    ("D1_sharp015", "sharp015_h500", "D", False, "Andy's [4]: UCERF's angle, a 1.5 deg sharp bend at x = 0"),
    ("D2_sharp05", "sharp05_h500", "D", False, "Andy's [4]: a 5 deg sharp bend at x = 0"),
    ("D3_arc05L10", "arc05L10_h500", "D", False, "Andy's [4]: 5 deg over a 10 km arc centred at x = 0"),
    ("D4_sharp015vs", "sharp015vs_h500", "D", False, "addition: a 1.5 deg sharp bend at x = 12.5 km, in the creeping part"),
    ("D5f_arc015L10_floor", "arc015L10_h500", "D", True, "D1f with the bend rounded: 1.5 deg over a 10 km arc centred at x = 0"),
]


def trace(u, theta_deg, arc, s_bend):
    """Point (x, y) [m] of the fault at u, its arc length from the centre of
    the straight west part (u = x there), as mesh/bend.geo builds it."""
    th = math.radians(theta_deg)
    s = u - s_bend
    if arc == 0.0:
        x, y = (s, 0.0) if s <= 0.0 else (s * math.cos(th), s * math.sin(th))
    else:
        R, s0 = arc / th, -arc / 2.0
        if s <= s0:
            x, y = s, 0.0
        elif s <= -s0:
            p = (s - s0) / R
            x, y = s0 + R * math.sin(p), R * (1.0 - math.cos(p))
        else:
            x = s0 + R * math.sin(th) + (s + s0) * math.cos(th)
            y = R * (1.0 - math.cos(th)) + (s + s0) * math.sin(th)
    return x + s_bend, y


def stations(mesh):
    """Fault stations on the trace: fltst_s<km> every 3 km of arc length,
    and the bend set: bend_w and bend_e 1 m from the vertex (or the arc's
    mid point) on the west and the east face, so each face keeps its own
    sigma_n'; bend_s<+-km> along the fault from the bend; the arc's ends."""
    th, arc, sb = MESHES[mesh]
    out = []
    for k in range(-12, 13, 3):
        out.append((f"fltst_s{k:+03d}", trace(k * 1e3, th, arc, sb)))
    # 1 m along the face's chord (the fault is straight segments of
    # LC_FAULT between points on the curve), so the point lies on that face
    v = trace(sb, th, arc, sb)
    for name, u in (("bend_w", sb - LC_FAULT), ("bend_e", sb + LC_FAULT)):
        n = trace(u, th, arc, sb)
        L = math.hypot(n[0] - v[0], n[1] - v[1])
        out.append((name, (v[0] + (n[0] - v[0]) / L, v[1] + (n[1] - v[1]) / L)))
    offs = [-3.0, -1.0, 1.0, 3.0]
    if sb + 3e3 > FAULT_HALF:
        offs = [-3.0, -1.0, 1.0, 2.0]       # the tip is 2.5 km past the bend
    if arc > 0.0:
        offs = [-arc / 2e3] + offs + [arc / 2e3]
    for d in offs:
        out.append((f"bend_s{d:+.1f}", trace(sb + d * 1e3, th, arc, sb)))
    return out


# the request each deck answers (Chunhui Zhao, 2026-10-31)
GO = "# 2026-10-31: \"go ahead with the simple cases, run them on local Mac\"):"
GO_CASE = {"D5f_arc015L10_floor":
           "# 2026-10-31: \"a smooth design around this bend, and rerun\"):"}


def header(case, mesh, regime, floor, row):
    lines = [
        f"# simple_cases/configs/{case}.toml — the simple cases (Chunhui Zhao,",
        GO_CASE.get(case, GO),
        f"# {row}.",
        f"# {'QD only' if regime == 'Q' else 'Fully dynamic (QD and dynamic regimes)'}, mesh {mesh}"
        + (", sigma_n' floor 1 MPa." if floor else ", no sigma_n' floor."),
        "# Made by simple_cases/scripts/make_decks.py from",
        "# config/fd_cycle/mapview_rs_fd.toml; the edits are listed there and in",
        "# simple_cases/README.md.  Plane strain, no 2.5D term, +-240 km box, all",
        "# walls absorbing in the dynamic regime, 12,000 yr.",
        "#",
        "# Run (Mac, from miniapps/seas): simple_cases/scripts/run_simple_mac.sh",
        "#",
        "# mapview_rs_fd.toml's notes follow (its mesh: the straight fault, +-60 km).",
        "#",
    ]
    return "\n".join(lines) + "\n"


def make(case, mesh, regime, floor, row, base):
    t = base

    def sub(old, new, count=1):
        nonlocal t
        assert t.count(old) == count, (case, old, t.count(old))
        t = t.replace(old, new)

    sub('path  = "mapview_rs/mesh/mapview_rs_500m.msh"', f'path  = "simple_cases/mesh/{mesh}.msh"')
    sub("# 1500 yr of 365 days.\n[time]\ntfinal     = 4.7304e10",
        f"# {T_FINAL_YR:,.0f} yr of 365 days (mapview_rs_fd.toml: 1500 yr).\n[time]\n"
        f"tfinal     = {T_FINAL_YR * YR:.6e}")
    if regime == "Q":
        sub("switching            = true", "switching            = false   # QD only (mapview_rs_fd.toml: true)")
    if floor:
        sub("penalty_factor       = 1.0       # both regimes' K\n",
            "penalty_factor       = 1.0       # both regimes' K\n"
            "sigma_n_floor_pa     = 1.0e6     # the 1 MPa floor of UCERF run8 (not in mapview_rs_fd.toml)\n")
    sub('prefix          = "mapview_rs_fd"', f'prefix          = "{case}"'
        + ("\nhandback_table  = true" if regime == "D" else ""))
    # stations: replace the nine fltst_x entries
    i0 = t.index("# On-fault stations every 3 km from -12 to 12 km.")
    i1 = t.index("# Off-fault stations")
    st = ["# On-fault stations on the trace (make_decks.py): fltst_s<km> every 3 km",
          "# of arc length from -12 to 12 km; bend_w / bend_e 1 m from the bend",
          "# vertex (arc: its mid point) on the west / east face; bend_s<km> along",
          "# the fault from the bend (arc: also its two ends)."]
    for name, (x, y) in stations(mesh):
        st += ["[[cycle.output.station]]", f'name = "{name}"', f"x    = [{x!r}, {y!r}]"]
    t = t[:i0] + "\n".join(st) + "\n\n" + t[i1:]
    sub("region_x3       = [-1.0, 1.0]", "region_x3       = [-1.0e4, 1.0e4]   # the bent limb reaches y = 1.31 km")
    sub("dt_dyn = 0.1\n", "dt_dyn = 0.05           # mapview_rs_fd.toml: 0.1 (sigma_n' against slip at the bend)\n")
    sub("dt_dyn  = 0.5\n", "dt_dyn  = 1.0           # mapview_rs_fd.toml: 0.5 (1 MB each; disk)\n")
    return header(case, mesh, regime, floor, row) + t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    base = open(BASE).read()
    bad = 0
    for case, mesh, regime, floor, row in CASES:
        text = make(case, mesh, regime, floor, row, base)
        path = f"simple_cases/configs/{case}.toml"
        if a.check:
            same = os.path.exists(path) and open(path).read() == text
            print(f"{'same' if same else 'DIFFERS'} {path}")
            bad += not same
        else:
            with open(path, "w") as f:
                f.write(text)
            print("wrote", path)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
