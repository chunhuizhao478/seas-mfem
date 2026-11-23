#!/usr/bin/env python3
"""Write the on-fault station files read by the Nitsche dynamic-rupture driver.

One CSV per problem, <problem>/stations/<problem>_fault_stations.csv (and
for TPV26/27 also <problem>/stations/<problem>_off_fault_stations.csv, the
surface stations of [output.off_fault_stations], same format):

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
  tpv205_2d  TPV205_210_Description_v4.pdf, the 2D benchmark: the seven
          on-fault stations at 7.5 km depth (x = 0, +-4.5, +-7.5, +-12 km),
          written into tpv205/stations/ with z = 0 (the 2-D mesh is the
          horizontal slice at that depth, in the (x, y) plane)
  tpv105_2d  SCEC_validation_TPV105.pdf, "Time Histories of Fields at Fault
          Stations", the seven 2-D stations (x = 0, +-9, +-15, +-18 km on
          the y = 7.5 km line of the 3-D fault), written into tpv105/stations/
          with z = 0 (the 2-D mesh tpv105/mesh/tpv105_2d.geo, fault on y = 0)
  tpv6    2007RuthRalphletter2.pdf, on-fault stations (five points; the
          same positions as the DRDG3D files in tpv6/benchmark_data/)
  tpv7    the same five points as TPV6
  tpv31   TPV31_32_Description_v03.pdf, on-fault stations (thirty points:
          strike 0, 6, 12 km at ten depths; the same positions as the
          SeisSol and EQdyna files in tpv31/benchmark_data/)
  tpv26   TPV26_27_Description_v13.pdf, Part 8 (twelve points; the same
          positions as the DRDG3D files in tpv26/ and tpv27/benchmark_data/;
          TPV27 uses this file too)
  tpv12   TPV12_13_Description_v6.pdf, Part 7 (ten points; 60 deg dip
          toward +Y, x3 = down-dip distance, not depth)
  tpv15   TPV14_15_Description_v08.pdf (fourteen points; main fault on
          Y = 0, branch at 30 deg toward -Y from the junction at the
          origin; x2 = distance from the junction along each fault)
  tpv22   TPV22_23_Description_v08.pdf (eighteen points; fault 1 on Y = 0,
          fault 2 on Y = -1.6 km; x2 = global X)
TPV104, TPV205, TPV6, TPV7, TPV31, TPV26, TPV12, TPV15 and TPV22 stations are
named like the SCEC reference files (x2_<strike km>_x3_<depth km>); the
complex-geometry ones add the fault (main_, branch_, fault1_, fault2_) where
a problem has more than one. The 2-D ones (tpv205_2d, tpv105_2d) keep
the 3-D names of the same points.

Usage (from miniapps/seas):
    python3 scripts/make_fault_stations.py                # all problems
    python3 scripts/make_fault_stations.py tpv205         # one problem
    python3 scripts/make_fault_stations.py --check        # files up to date?

  tpv34   TPV34_Description_v10.pdf, Part 5 (thirty-five points: strike -12,
          -6, 0, 6, 12 km x depth 0, 1, 2.4, 5, 7.5, 10, 12 km), named as the
          description names them (faultst<strike>dp<depth>, units of 100 m)

Off-fault stations ([output.body_stations]) go to
<problem>/stations/<problem>_body_stations.csv, same columns:
  tpv34   TPV34_Description_v10.pdf, Part 6 (fifty-six points: 28 boreholes,
          each at the surface and at 2.4 km depth), named
          body<normal>st<strike>dp<depth>, units of 100 m.  The description's
          z (fault normal, positive on the far side) is the mesh Y.
  tpv205  uploadTPV5_v2.pdf (four: 3 km on the far side, strike -12 and
          12 km, depth 0 and 7.5 km)
  tpv102, tpv104  SCEC_validation_ageing_law.pdf / _slip_law.pdf (six
          free-surface points; names from uploadTPV101.pdf / uploadTPV103.pdf,
          whose "body-060s-120dp000" typos are written with "st")
  tpv12   TPV12_13_Description_v6.pdf, Part 8 (twelve; horizontal distance
          from the fault trace, or from the fault at 0.3 km depth for the
          four 0.3 km stations)
  tpv22   TPV22_23_Description_v08.pdf, Part 7 (six surface points;
          distance from fault 1)
  tpv31   TPV31_32_Description_v03.pdf, Part 6 (eighteen: six boreholes at
          0, 0.5 and 2.4 km depth, all on the near side)
  All of them: names body<normal>st<strike>dp<depth> (units of 100 m), the
  description's far side = +Y.

Off-fault surface stations ([output.off_fault_stations]) go to
<problem>/stations/<problem>_off_fault_stations.csv, same columns:
  tpv26   TPV26_27_Description_v13.pdf, Part 9 (six surface points),
          named body<normal>st<strike>dp000 (units of 100 m); TPV27 uses it too.
  tpv15   TPV14_15_Description_v08.pdf, Part 6 (eleven surface points;
          distance from the main fault, strike from the junction; this
          writer because TPV15 runs with LTS)

Scope: Chunhui Zhao, 2026-10-16 (TPV205/TPV104 plan rev 1.3, decision Q4); TPV6, TPV7
and TPV31 added 2026-10-19 (TPV6/TPV7/TPV31 plan rev 1); TPV34 and the off-fault
files added 2026-11-12 (TPV34 plan rev 2); the other off-fault files
2026-11-20 (SCEC submission runs); TPV26 and its
off-fault file 2026-10-20 (rheology plan rev 1, phase R0).
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
_TPV205_2D_STRIKE = [-12, -7.5, -4.5, 0, 4.5, 7.5, 12]   # all at 7.5 km depth
_TPV105_2D_STRIKE = [-18, -15, -9, 0, 9, 15, 18]         # all at y = 7.5 km
_TPV6_POINTS = [(0, 0), (-12, 0), (12, 0), (-12, 7.5), (12, 7.5)]
_TPV31_POINTS = [(s, d) for s in (0, 6, 12) for d in (0, 0.2, 0.5, 1, 2.4, 3, 5, 7.5, 10, 12)]
_TPV34_POINTS = [(s, d) for s in (-12, -6, 0, 6, 12) for d in (0, 1, 2.4, 5, 7.5, 10, 12)]
# off-fault: (fault normal km, strike km, depth km), description Part 6
_TPV34_BODY = [(n, s, d)
               for n, strikes in ((-15, (-15, 0, 15)), (-9, (-20, -10, 0, 10, 20)),
                                  (-3, (-20, -10, 0, 10, 20)), (0, (-20, 20)),
                                  (3, (-20, -10, 0, 10, 20)), (9, (-20, -10, 0, 10, 20)),
                                  (15, (-15, 0, 15)))
               for s in strikes for d in (0, 2.4)]
# Off-fault stations of the other descriptions: (fault normal km, strike km, depth km);
# + normal = the far side = +Y (SCEC submission runs, Chunhui Zhao, 2026-11-20).
# TPV5: uploadTPV5_v2.pdf (four, all on the far side)
_TPV205_BODY = [(3, -12, 0), (3, 12, 0), (3, -12, 7.5), (3, 12, 7.5)]
# TPV102 / TPV104: SCEC_validation_ageing_law.pdf / _slip_law.pdf free-surface (x, z)
# points, names from uploadTPV101.pdf / uploadTPV103.pdf
_TPV102_BODY = [(-6, -12, 0), (-9, 0, 0), (-6, 12, 0), (6, -12, 0), (9, 0, 0), (6, 12, 0)]
# TPV12: TPV12_13_Description_v6.pdf, Part 8 (twelve)
_TPV12_BODY = ([(n, 0, 0) for n in (-3, -2, -1, 1, 2, 3)]
               + [(n, 0, 0.3) for n in (-1, -0.5, 0.5, 1)] + [(-3, 12, 0), (3, 12, 0)])
# TPV15: TPV14_15_Description_v08.pdf, Part 6 (eleven, surface; normal from the main fault;
# written to the [output.off_fault_stations] file, see OFF_FAULT)
_TPV15_BODY = [(3, -2, 0), (-3, -2, 0), (3, 2, 0), (-0.6, 2, 0), (-4.2, 2, 0), (3, 5, 0),
               (-1.4, 5, 0), (-5.9, 5, 0), (3, 8, 0), (-2.3, 8, 0), (-7.6, 8, 0)]
# TPV22: TPV22_23_Description_v08.pdf, Part 7 (six, surface; normal from fault 1)
_TPV22_BODY = [(3, -5, 0), (-3, -5, 0), (3, 5, 0), (-3, 5, 0), (3, 15, 0), (-3, 15, 0)]
# TPV31: TPV31_32_Description_v03.pdf, Part 6 (eighteen: six boreholes at 0, 0.5, 2.4 km)
_TPV31_BODY = [(n, s, d) for n, s in ((-3, 0), (-3, 12), (-9, 0), (-9, 15), (-15, 0), (-15, 15))
               for d in (0, 0.5, 2.4)]
_TPV26_POINTS = [(-15, 10), (-5, 0), (-5, 5), (-5, 10), (-5, 15), (0, 10), (5, 10),
                 (10, 10), (15, 0), (15, 5), (15, 10), (15, 15)]
# TPV26/27 Part 9: (fault-normal km, strike km) on the surface; + = far side = +Y
_TPV26_OFF_POINTS = [(3, -5), (-3, -5), (3, 5), (-3, 5), (3, 15), (-3, 15)]


def body_name_surface(n, s):
    """TPV26/27 off-fault names: body<normal x10, 3 digits>st<strike x10, 3 digits>dp000."""
    def f(v):
        return ("-" if v < 0 else "") + f"{round(abs(v) * 10):03d}"
    return f"body{f(n)}st{f(s)}dp000"


def _km(v):
    return f"{v:g}"


def flt_name(s, d):
    """The driver's historical TPV102 names: flt_<strike>_<depth>, n = minus."""
    return "flt_" + ("n" if s < 0 else "") + _km(abs(s)) + "_" + _km(d)


def scec_name(s, d):
    """x2_<strike km>_x3_<depth km>, as in the SCEC reference file names."""
    return f"x2_{_km(s)}_x3_{_km(d)}"


def _h(v, signed=True):
    """SCEC station-name field: units of 100 m, three digits, minus sign when negative."""
    n = int(round(10 * abs(v)))
    return ("-" if (signed and v < 0) else "") + f"{n:03d}"


def faultst_name(s, d):
    """faultst<strike>dp<depth>, as the TPV34 description names the on-fault stations."""
    return f"faultst{_h(s)}dp{_h(d, False)}"


def body_name(n, s, d):
    """body<normal>st<strike>dp<depth>, as the TPV34 description names the off-fault stations."""
    return f"body{_h(n)}st{_h(s)}dp{_h(d, False)}"


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
    "tpv6": ("SCEC TPV6 on-fault stations", plane(_TPV6_POINTS, scec_name)),
    "tpv7": ("SCEC TPV7 on-fault stations", plane(_TPV6_POINTS, scec_name)),
    "tpv31": ("SCEC TPV31 on-fault stations", plane(_TPV31_POINTS, scec_name)),
    "tpv34": ("SCEC TPV34 on-fault stations", plane(_TPV34_POINTS, faultst_name)),
    "tpv26": ("SCEC TPV26/TPV27 on-fault stations", plane(_TPV26_POINTS, scec_name)),
    "tpv12": ("SCEC TPV12 on-fault stations", framed("tpv12", _TPV12_POINTS, False)),
    "tpv15": ("SCEC TPV15 on-fault stations", framed("tpv15", _TPV15_POINTS, True)),
    "tpv22": ("SCEC TPV22 on-fault stations", framed("tpv22", _TPV22_POINTS, True)),
    "tpv205_2d": ("SCEC TPV205-2D on-fault stations (the 7.5 km depth slice)",
                  [(scec_name(s, 7.5), 1000.0 * s, 0.0, 0.0) for s in _TPV205_2D_STRIKE]),
    "tpv105_2d": ("SCEC TPV105-2D on-fault stations (the y = 7.5 km line)",
                  [(scec_name(s, 7.5), 1000.0 * s, 0.0, 0.0) for s in _TPV105_2D_STRIKE]),
}

# Off-fault stations: (name, strike km, fault normal km, depth km).  The
# normal coordinate is the mesh Y in km (+ = the description's far side).
BODY = {
    "tpv34": ("SCEC TPV34 off-fault stations",
              [(body_name(n, s, d), s, n, d) for n, s, d in _TPV34_BODY]),
    "tpv205": ("SCEC TPV205 (TPV5) off-fault stations",
               [(body_name(n, s, d), s, n, d) for n, s, d in _TPV205_BODY]),
    "tpv102": ("SCEC TPV102 off-fault (free-surface) stations",
               [(body_name(n, s, d), s, n, d) for n, s, d in _TPV102_BODY]),
    "tpv104": ("SCEC TPV104 off-fault (free-surface) stations",
               [(body_name(n, s, d), s, n, d) for n, s, d in _TPV102_BODY]),
    "tpv31": ("SCEC TPV31 off-fault stations",
              [(body_name(n, s, d), s, n, d) for n, s, d in _TPV31_BODY]),
    # TPV12: the distance is horizontal, from the fault trace for the surface
    # stations and from the fault at the station's depth for the 0.3 km ones
    # (Part 8 diagram text); the fault dips toward +Y, so at depth d it is at
    # Y = d cot 60 and the station sits at Y = n + d cot 60.
    "tpv12": ("SCEC TPV12 off-fault stations",
              [(body_name(n, s, d), s, round(n + d / math.tan(_DIP60), 9), d)
               for n, s, d in _TPV12_BODY]),
    "tpv22": ("SCEC TPV22 off-fault stations",
              [(body_name(n, s, d), s, n, d) for n, s, d in _TPV22_BODY]),
}
_BODY_NOTES = {
    "tpv12": "# Fault trace on the X axis, 60 deg dip toward +Y (hanging wall = far side = +Y); "
             "a station n km off the fault at depth d km sits at Y = n + d cot 60.",
    "tpv22": "# Fault 1 on Y = 0, fault 2 on Y = -1600 m; distances from fault 1, x = global X.",
}
# 2-D problems: (folder, the mesh plane's z); the depth is in the name only
PLANE_2D = {"tpv205_2d": ("tpv205", 0.0), "tpv105_2d": ("tpv105", 0.0)}


def _num(v):
    v = float(v) + 0.0          # no "-0"
    # whole metres as before; other values (the framed stations, rounded to
    # 1 um) in the shortest form that reads back to the same double
    return f"{v:.17g}" if v.is_integer() else repr(v)


def render(problem):
    title, rows = PROBLEMS[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand."]
    if problem in PLANE_2D:
        out.append("# 2-D mesh in the (x, y) plane: X along strike, Y fault normal (m), z = 0; "
                   "fault on Y = 0.")
    elif problem in _FRAME_NOTES:
        out += ["# Global frame of the mesh: X, Y, Z up (m).", _FRAME_NOTES[problem]]
    else:
        out.append("# Global frame of the mesh: X along strike, Y fault normal, Z up (m); fault on Y = 0.")
    out.append("name,x,y,z")
    for name, x, y, z in rows:
        out.append(f"{name},{_num(x)},{_num(y)},{_num(z)}")
    return "\n".join(out) + "\n"


def render_body(problem):
    title, rows = BODY[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand."]
    if problem in _BODY_NOTES:
        out += ["# Global frame of the mesh: X along strike, Y fault normal (+ = far side), Z up (m).",
                _BODY_NOTES[problem]]
    else:
        out.append("# Global frame of the mesh: X along strike, Y fault normal, Z up (m); fault on Y = 0.")
    out.append("name,x,y,z")
    for name, s, n, d in rows:
        out.append(f"{name},{_num(1000.0 * s)},{_num(round(1000.0 * n, 6))},{_num(0.0 - 1000.0 * d)}")
    return "\n".join(out) + "\n"


def target(root, problem, kind="fault"):
    folder = PLANE_2D[problem][0] if problem in PLANE_2D else problem
    return os.path.join(root, folder, "stations", f"{problem}_{kind}_stations.csv")


# off-fault surface stations ([output.off_fault_stations], TPV26/27 writer):
# problem -> (title, [(name, strike km, normal km)])
OFF_FAULT = {
    "tpv26": ("SCEC TPV26/TPV27 off-fault surface stations",
              [(body_name_surface(n, s), s, n) for n, s in _TPV26_OFF_POINTS]),
    # TPV15 runs with LTS, whose volume is synchronised only at sync points; this
    # writer (unlike [output.body_stations]) writes there only.  Distances from
    # the main fault (Y = 0), x from the junction.
    "tpv15": ("SCEC TPV15 off-fault surface stations (normal distance from the main fault)",
              [(body_name(n, s, d), s, n) for n, s, d in _TPV15_BODY]),
}


def render_off(problem):
    title, rows = OFF_FAULT[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand.",
           "# Global frame of the mesh: X along strike, Y fault normal (+ = far side), Z up (m).",
           "name,x,y,z"]
    for name, s, n in rows:
        out.append(f"{name},{_num(1000.0 * s)},{_num(1000.0 * n)},{_num(0.0)}")
    return "\n".join(out) + "\n"


def target_off(root, problem):
    return os.path.join(root, problem, "stations", f"{problem}_off_fault_stations.csv")


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
        files = [(os.path.normpath(target(a.out_root, p)), render(p), len(PROBLEMS[p][1]))]
        if p in BODY:
            files.append((os.path.normpath(target(a.out_root, p, "body")), render_body(p),
                          len(BODY[p][1])))
        if p in OFF_FAULT:
            files.append((os.path.normpath(target_off(a.out_root, p)), render_off(p),
                          len(OFF_FAULT[p][1])))
        for path, text, count in files:
            if a.check:
                same = os.path.exists(path) and open(path).read() == text
                print(f"{'ok     ' if same else 'DIFFERS'} {path}")
                bad += 0 if same else 1
                continue
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                f.write(text)
            print(f"wrote {path} ({count} stations)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
