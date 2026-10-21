#!/usr/bin/env python3
"""Write the on-fault station files read by the Nitsche dynamic-rupture driver.

One CSV per problem, <problem>/stations/<problem>_fault_stations.csv:

    # comment lines start with '#'
    name,x,y,z
    x2_0_x3_7.5,0,0,-7500
    ...

Coordinates are metres in the global frame of the problem's mesh (X along
strike, Y fault normal, Z up).  TPV102, TPV104 and TPV205 lie on Y = 0, so
depth = -Z.  The complex-geometry problems (TPV12, TPV15, TPV22) give each
station by (along-strike km, down-dip km) on a fault with its own frame:
an origin on the free surface, a unit strike vector and a unit down-dip
vector (complex-geometries plan rev 2, Section 5.8).
The driver snaps each station to the nearest fault node in 3-D, averages
that node's DG copies, and aborts when the node is farther than the deck's
[output.stations] snap_tolerance_m.  Station positions live here and in the
files this writes, never in the C++.

The positions are the on-fault station lists of the benchmark descriptions:
  tpv102  SCEC TPV102 (same nine points as TPV103/TPV104); names keep the
          driver's historical flt_* names so the output file names do not
          change
  tpv104  SCEC_validation_slip_law.pdf, "Time Histories of Fields at Fault
          Stations" (nine points)
  tpv205  TPV5_forwebsite.pdf, "Diagram showing locations of on-fault
          stations" (sixteen points; the same positions as the DRDG3D files
          in tpv205/benchmark_data/)
  tpv12   TPV12_13_Description_v6.pdf, Part 7 (ten points; 60 deg dip
          toward +Y, x3 = down-dip distance, not depth)
  tpv15   TPV14_15_Description_v08.pdf (fourteen points; main fault on
          Y = 0, branch at 30 deg toward -Y from the junction at the
          origin; x2 = distance from the junction along each fault)
  tpv22   TPV22_23_Description_v08.pdf (eighteen points; fault 1 on Y = 0,
          fault 2 on Y = -1.6 km; x2 = global X)
TPV104 and TPV205 stations are named like the SCEC reference files
(x2_<strike km>_x3_<depth km>); the complex-geometry ones add the fault
(main_, branch_, fault1_, fault2_) where a problem has more than one.

Usage (from miniapps/seas):
    python3 scripts/make_fault_stations.py                # all problems
    python3 scripts/make_fault_stations.py tpv205         # one problem
    python3 scripts/make_fault_stations.py --check        # files up to date?

Scope: Chunhui Zhao, 2026-10-16 (TPV205/TPV104 plan rev 1.3, decision Q4).
"""
import argparse
import math
import os
import sys

# (strike km, depth km) in the order the files list them
_TPV102_POINTS = [(0, 3), (0, 7.5), (0, 12), (9, 7.5), (12, 3), (12, 12),
                  (-9, 7.5), (-12, 3), (-12, 12)]
_TPV205_POINTS = [(-12, 0), (-12, 7.5), (-7.5, 0), (-7.5, 7.5), (-4.5, 0), (-4.5, 7.5),
                  (0, 0), (0, 3), (0, 7.5), (0, 12), (4.5, 0), (4.5, 7.5),
                  (7.5, 0), (7.5, 7.5), (12, 0), (12, 7.5)]


def _km(v):
    return f"{v:g}"


def flt_name(s, d):
    """The driver's historical TPV102 names: flt_<strike>_<depth>, n = minus."""
    return "flt_" + ("n" if s < 0 else "") + _km(abs(s)) + "_" + _km(d)


def scec_name(s, d):
    """x2_<strike km>_x3_<depth km>, as in the SCEC reference file names."""
    return f"x2_{_km(s)}_x3_{_km(d)}"


# Fault frames of the complex-geometry problems: origin (m), unit strike,
# unit down-dip; a station at (s, d) km sits at origin + 1000 (s strike + d dip).
_DIP60 = math.radians(60.0)
_BR30 = math.radians(30.0)
_FRAMES = {
    "tpv12": {"fault": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                        (0.0, math.cos(_DIP60), -math.sin(_DIP60)))},
    "tpv15": {"main": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0)),
              "branch": ((0.0, 0.0, 0.0), (math.cos(_BR30), -math.sin(_BR30), 0.0),
                         (0.0, 0.0, -1.0))},
    "tpv22": {"fault1": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0)),
              "fault2": ((0.0, -1600.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0))},
}
_FRAME_NOTES = {
    "tpv12": "# fault: trace on the X axis, 60 deg dip toward +Y; x3 = down-dip km (depth = x3 sin 60).",
    "tpv15": "# main: Y = 0, x2 from the junction (origin); branch: strike (cos 30, -sin 30, 0) from the junction.",
    "tpv22": "# fault1: Y = 0; fault2: Y = -1600 m; x2 = global X km, x3 = depth km.",
}
# (fault, strike km, down-dip km)
_TPV12_POINTS = [("fault", s, d) for s, d in
                 [(0, 0), (4.5, 0), (12, 0), (0, 1.5), (0, 3), (0, 4.5), (0, 7.5), (4.5, 7.5),
                  (12, 7.5), (0, 12)]]
_TPV15_POINTS = ([("main", s, d) for d in (0, 7.5) for s in (-2, 2, 5, 9)]
                 + [("branch", s, d) for d in (0, 7.5) for s in (2, 5, 9)])
_TPV22_POINTS = ([("fault1", s, d) for s, d in
                  [(-10, 0), (0, 0), (-10, 5), (-10, 10), (-5, 10), (0, 10), (-10, 15)]]
                 + [("fault2", s, d) for s, d in
                    [(0, 0), (5, 0), (20, 0), (4, 5), (5, 5), (0, 10), (5, 10), (6.5, 10),
                     (10, 10), (20, 10), (5, 15)]])


def framed(problem, points, tagged):
    """-> [(name, x, y, z)] in m for (fault, strike km, down-dip km) points."""
    rows = []
    for f, s, d in points:
        o, st, dp = _FRAMES[problem][f]
        x = [round(o[c] + 1000.0 * (s * st[c] + d * dp[c]), 6) for c in range(3)]   # um
        rows.append(((f + "_" if tagged else "") + scec_name(s, d), x[0], x[1], x[2]))
    return rows


def plane(points, name):
    """-> [(name, x, y, z)] in m for (strike km, depth km) points on Y = 0."""
    return [(name(s, d), 1000.0 * s, 0.0, 0.0 - 1000.0 * d) for s, d in points]


PROBLEMS = {
    "tpv102": ("SCEC TPV102 on-fault stations", plane(_TPV102_POINTS, flt_name)),
    "tpv104": ("SCEC TPV104 on-fault stations", plane(_TPV102_POINTS, scec_name)),
    "tpv205": ("SCEC TPV205 (TPV5) on-fault stations", plane(_TPV205_POINTS, scec_name)),
    "tpv12": ("SCEC TPV12 on-fault stations", framed("tpv12", _TPV12_POINTS, False)),
    "tpv15": ("SCEC TPV15 on-fault stations", framed("tpv15", _TPV15_POINTS, True)),
    "tpv22": ("SCEC TPV22 on-fault stations", framed("tpv22", _TPV22_POINTS, True)),
}


def _num(v):
    v = float(v) + 0.0          # no "-0"
    # whole metres as before; other values (the framed stations, rounded to
    # 1 um) in the shortest form that reads back to the same double
    return f"{v:.17g}" if v.is_integer() else repr(v)


def render(problem):
    title, rows = PROBLEMS[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand."]
    if problem in _FRAME_NOTES:
        out += ["# Global frame of the mesh: X, Y, Z up (m).", _FRAME_NOTES[problem]]
    else:
        out.append("# Global frame of the mesh: X along strike, Y fault normal, Z up (m); fault on Y = 0.")
    out.append("name,x,y,z")
    for name, x, y, z in rows:
        out.append(f"{name},{_num(x)},{_num(y)},{_num(z)}")
    return "\n".join(out) + "\n"


def target(root, problem):
    return os.path.join(root, problem, "stations", f"{problem}_fault_stations.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("problems", nargs="*", default=sorted(PROBLEMS))
    ap.add_argument("--out-root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    ap.add_argument("--check", action="store_true", help="exit 1 if a file differs from what would be written")
    a = ap.parse_args()
    bad = 0
    for p in a.problems:
        if p not in PROBLEMS:
            sys.exit(f"unknown problem {p!r}; known: {', '.join(sorted(PROBLEMS))}")
        path = os.path.normpath(target(a.out_root, p))
        text = render(p)
        if a.check:
            same = os.path.exists(path) and open(path).read() == text
            print(f"{'ok     ' if same else 'DIFFERS'} {path}")
            bad += 0 if same else 1
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(text)
        print(f"wrote {path} ({len(PROBLEMS[p][1])} stations)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
