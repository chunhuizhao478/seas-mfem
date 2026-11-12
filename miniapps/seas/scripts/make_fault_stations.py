#!/usr/bin/env python3
"""Write the on-fault station files read by the Nitsche dynamic-rupture driver.

One CSV per problem, <problem>/stations/<problem>_fault_stations.csv:

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
  tpv6    2007RuthRalphletter2.pdf, on-fault stations (five points; the
          same positions as the DRDG3D files in tpv6/benchmark_data/)
  tpv7    the same five points as TPV6
  tpv31   TPV31_32_Description_v03.pdf, on-fault stations (thirty points:
          strike 0, 6, 12 km at ten depths; the same positions as the
          SeisSol and EQdyna files in tpv31/benchmark_data/)
TPV104, TPV205, TPV6, TPV7 and TPV31 stations are named like the SCEC
reference files (x2_<strike km>_x3_<depth km>).

Usage (from miniapps/seas):
    python3 scripts/make_fault_stations.py                # all problems
    python3 scripts/make_fault_stations.py tpv205         # one problem
    python3 scripts/make_fault_stations.py --check        # files up to date?

Scope: Chunhui Zhao, 2026-10-16 (TPV205/TPV104 plan rev 1.3, decision Q4); TPV6, TPV7
and TPV31 added 2026-10-19 (TPV6/TPV7/TPV31 plan rev 1).
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
_TPV6_POINTS = [(0, 0), (-12, 0), (12, 0), (-12, 7.5), (12, 7.5)]
_TPV31_POINTS = [(s, d) for s in (0, 6, 12) for d in (0, 0.2, 0.5, 1, 2.4, 3, 5, 7.5, 10, 12)]


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
    "tpv6": ("SCEC TPV6 on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV6_POINTS]),
    "tpv7": ("SCEC TPV7 on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV6_POINTS]),
    "tpv31": ("SCEC TPV31 on-fault stations", [(scec_name(s, d), s, d) for s, d in _TPV31_POINTS]),
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
