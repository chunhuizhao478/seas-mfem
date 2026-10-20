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
strike, Y fault normal, Z up; these faults lie on Y = 0, so depth = -Z).
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
  tpv26   TPV26_27_Description_v13.pdf, Part 8 (twelve points; the same
          positions as the DRDG3D files in tpv26/ and tpv27/benchmark_data/;
          TPV27 uses this file too)
TPV104, TPV205 and TPV26 stations are named like the SCEC reference files
(x2_<strike km>_x3_<depth km>).

Usage (from miniapps/seas):
    python3 scripts/make_fault_stations.py                # all problems
    python3 scripts/make_fault_stations.py tpv205         # one problem
    python3 scripts/make_fault_stations.py --check        # files up to date?

Scope: Chunhui Zhao, 2026-10-16 (TPV205/TPV104 plan rev 1.3, decision Q4);
TPV26 2026-10-20 (rheology plan rev 1, phase R0).
"""
import argparse
import os
import sys

# (strike km, depth km) in the order the files list them
_TPV102_POINTS = [(0, 3), (0, 7.5), (0, 12), (9, 7.5), (12, 3), (12, 12),
                  (-9, 7.5), (-12, 3), (-12, 12)]
_TPV205_POINTS = [(-12, 0), (-12, 7.5), (-7.5, 0), (-7.5, 7.5), (-4.5, 0), (-4.5, 7.5),
                  (0, 0), (0, 3), (0, 7.5), (0, 12), (4.5, 0), (4.5, 7.5),
                  (7.5, 0), (7.5, 7.5), (12, 0), (12, 7.5)]
_TPV26_POINTS = [(-15, 10), (-5, 0), (-5, 5), (-5, 10), (-5, 15), (0, 10), (5, 10),
                 (10, 10), (15, 0), (15, 5), (15, 10), (15, 15)]
# TPV26/27 Part 9: (fault-normal km, strike km) on the surface; + = far side = +Y
_TPV26_OFF_POINTS = [(3, -5), (-3, -5), (3, 5), (-3, 5), (3, 15), (-3, 15)]


def body_name(n, s):
    """SCEC off-fault names: body<normal x10, 3 digits>st<strike x10, 3 digits>dp000."""
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


PROBLEMS = {
    "tpv102": ("SCEC TPV102 on-fault stations", [(flt_name(s, d), s, d) for s, d in _TPV102_POINTS]),
    "tpv104": ("SCEC TPV104 on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV102_POINTS]),
    "tpv205": ("SCEC TPV205 (TPV5) on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV205_POINTS]),
    "tpv26": ("SCEC TPV26/TPV27 on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV26_POINTS]),
}


def _num(v):
    v = float(v) + 0.0          # no "-0"
    return f"{v:.17g}"


def render(problem):
    title, rows = PROBLEMS[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand.",
           "# Global frame of the mesh: X along strike, Y fault normal, Z up (m); fault on Y = 0.",
           "name,x,y,z"]
    for name, s, d in rows:
        out.append(f"{name},{_num(1000.0 * s)},{_num(0.0)},{_num(0.0 - 1000.0 * d)}")
    return "\n".join(out) + "\n"


# off-fault (surface) stations: problem -> (title, [(name, strike km, normal km)])
OFF_FAULT = {
    "tpv26": ("SCEC TPV26/TPV27 off-fault surface stations",
              [(body_name(n, s), s, n) for n, s in _TPV26_OFF_POINTS]),
}


def render_off(problem):
    title, rows = OFF_FAULT[problem]
    out = [f"# {title}, written by scripts/make_fault_stations.py; do not edit by hand.",
           "# Global frame of the mesh: X along strike, Y fault normal (+ = far side), Z up (m).",
           "name,x,y,z"]
    for name, s, n in rows:
        out.append(f"{name},{_num(1000.0 * s)},{_num(1000.0 * n)},{_num(0.0)}")
    return "\n".join(out) + "\n"


def target(root, problem):
    return os.path.join(root, problem, "stations", f"{problem}_fault_stations.csv")


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
        jobs = [(target(a.out_root, p), render(p), len(PROBLEMS[p][1]))]
        if p in OFF_FAULT:
            jobs.append((target_off(a.out_root, p), render_off(p), len(OFF_FAULT[p][1])))
        for path, text, n in jobs:
            path = os.path.normpath(path)
            if a.check:
                same = os.path.exists(path) and open(path).read() == text
                print(f"{'ok     ' if same else 'DIFFERS'} {path}")
                bad += 0 if same else 1
                continue
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                f.write(text)
            print(f"wrote {path} ({n} stations)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
