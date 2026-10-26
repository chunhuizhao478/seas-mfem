#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# scec_bp5_export.py — fully dynamic BP5 cycle: the SCEC upload files of the
# BP5 description (section 4) from a seas_fd_cycle_driver run with
# [cycle.benchmark] (cycle/bench_output.hpp).  Standard library only.
#
#   python3 cycle/scec_bp5_export.py RUN_DIR PREFIX OUT_DIR \
#       --modeler "Chunhui Zhao" --element-size "1000 m" --order 2 \
#       [--problem "SEAS Benchmark BP5-FD"] [--version COMMIT] [--event 1] \
#       [--region-x2 -34000 34000] [--region-x3 0 20000]
#
# Writes into OUT_DIR (made if missing; existing files of these names are
# replaced):
#   fltst_strk*dp*.dat   section 4.1, from RUN_DIR/PREFIX_fltst_*.txt
#   blkst_strk*.dat      section 4.2, from RUN_DIR/PREFIX_blkst_*.txt
#   global.dat           section 4.3, from PREFIX_source.txt
#   catalog.dat          section 4.4, from PREFIX_events.txt
#   slip_2_depth.dat ... section 4.5, from PREFIX_profile_<name>_<q>.txt
#                        (profile names "depth" and "strike")
#   rupture.dat          section 4.6, from PREFIX_rupture_event<event>.txt
# and prints one line per file.  Numbers are written as the description
# recommends for C users: 21.13E for time (and the off-fault velocities),
# 14.6E for the rest.
#
# Rows: the time series keep only rows whose time is larger than the
# previous kept row's (the description wants increasing time; a regime
# switch writes the switch time twice, the dynamic end state and the QD
# start).  The station, off-fault and source files of a run have the same
# rows, so they keep the same times.  rupture.dat keeps the fault DOFs in
# the region (default: the domain A of section 4.6) and one row per
# position: the fault DOFs are per face, so a node shared by faces appears
# once per face; the earliest time is kept.  The converter changes no
# value otherwise: the run writes the components, signs and units of the
# description (slip_2, slip_3 = kernel strike, dip, as the station files).

import argparse
import datetime
import glob
import math
import os
import re
import sys

T_FMT = "%21.13E"
V_FMT = "%14.6E"


def read_table(path):
    """Comment lines and the rows of a RowFile (cycle/fd_cycle_output.hpp)."""
    header, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            s = line.split()
            if s:
                rows.append([float(x) for x in s])
    return header, rows


def commit_of(header):
    for h in header:
        m = re.search(r"commit (\S+)", h)
        if m:
            return m.group(1)
    return "unknown"


def increasing(rows):
    out, last = [], -math.inf
    for r in rows:
        if r[0] > last:
            out.append(r)
            last = r[0]
    return out, len(rows) - len(out)


def step_stats(rows):
    dts = [b[0] - a[0] for a, b in zip(rows, rows[1:])]
    return (min(dts), max(dts)) if dts else (0.0, 0.0)


class Exporter:
    def __init__(self, a):
        self.a = a
        self.date = datetime.datetime.now().strftime("%Y/%m/%d")
        self.version = a.version

    def head(self, f, extra):
        f.write("# This is the file header:\n")
        f.write("# problem=%s\n" % self.a.problem)
        f.write("# code=%s\n" % self.a.code)
        f.write("# version=%s\n" % self.version)
        f.write("# modeler=%s\n" % self.a.modeler)
        f.write("# date=%s\n" % self.date)
        f.write("# element_size=%s\n" % self.a.element_size)
        f.write("# polynomial_degree=%d\n" % self.a.order)
        for e in extra:
            f.write("# %s\n" % e)

    def src(self, name):
        return os.path.join(self.a.run_dir, self.a.prefix + name)

    def dst(self, name):
        return os.path.join(self.a.out_dir, name)

    def series(self, src, dst, location, columns, fields, fmts):
        header, rows = read_table(src)
        if not self.version:
            self.version = commit_of(header)
        rows, dropped = increasing(rows)
        if rows and len(rows[0]) != len(fields):
            sys.exit("%s: %d columns, expected %d" % (src, len(rows[0]), len(fields)))
        lo, hi = step_stats(rows)
        with open(dst, "w") as f:
            self.head(f, ["location=%s" % location,
                          "minimum_time_step=%.6E" % lo,
                          "maximum_time_step=%.6E" % hi,
                          "num_time_steps=%d" % len(rows)] +
                      ["Column #%d = %s" % (i + 1, c) for i, c in enumerate(columns)] +
                      ["The line below lists the names of the data fields"])
            f.write(" ".join(fields) + "\n")
            f.write("# Here is the time-series data.\n")
            for r in rows:
                f.write(" ".join(fm % v for fm, v in zip(fmts, r)) + "\n")
        print("%s: %d rows (%d repeated times dropped) from %s" %
              (os.path.basename(dst), len(rows), dropped, os.path.basename(src)))

    def on_fault(self):
        files = sorted(glob.glob(self.src("_fltst_*.txt")))
        for src in files:
            name = os.path.basename(src)[len(self.a.prefix) + 1:-4]
            m = re.match(r"fltst_strk([+-]\d+)dp([+-]?\d+)$", name)
            if not m:
                sys.exit("%s: not a BP5 on-fault station name" % name)
            loc = "on fault, %dkm along strike, %dkm depth" % (int(m.group(1)), int(m.group(2)))
            self.series(src, self.dst(name + ".dat"), loc,
                        ["Time (s)", "Slip_2 (m)", "Slip_3 (m)", "Slip_rate_2 (log10 m/s)",
                         "Slip_rate_3 (log10 m/s)", "Shear_stress_2 (MPa)",
                         "Shear_stress_3 (MPa)", "State (log10 s)"],
                        ["t", "slip_2", "slip_3", "slip_rate_2", "slip_rate_3",
                         "shear_stress_2", "shear_stress_3", "state"],
                        [T_FMT] + [V_FMT] * 7)
        return len(files)

    def off_fault(self):
        files = sorted(glob.glob(self.src("_blkst_*.txt")))
        for src in files:
            name = os.path.basename(src)[len(self.a.prefix) + 1:-4]
            m = re.match(r"blkst_strk([+-]\d+)fn([+-]\d+)dp([+-]?\d+)$", name)
            if not m:
                sys.exit("%s: not a BP5 off-fault station name" % name)
            loc = ("off fault, %dkm along strike, %dkm away from the fault, %dkm depth" %
                   (int(m.group(1)), int(m.group(2)), int(m.group(3))))
            self.series(src, self.dst(name + ".dat"), loc,
                        ["Time (s)", "Displacement_1 (m)", "Displacement_2 (m)",
                         "Displacement_3 (m)", "Velocity_1 (m/s)", "Velocity_2 (m/s)",
                         "Velocity_3 (m/s)"],
                        ["t", "disp_1", "disp_2", "disp_3", "vel_1", "vel_2", "vel_3"],
                        [T_FMT] + [V_FMT] * 3 + [T_FMT] * 3)
        return len(files)

    def source(self):
        src = self.src("_source.txt")
        self.series(src, self.dst("global.dat"),
                    "VW patch + transition zone (|x2| <= %g m, %g <= x3 <= %g m)" %
                    (max(abs(self.a.region_x2[0]), abs(self.a.region_x2[1])),
                     self.a.region_x3[0], self.a.region_x3[1]),
                    ["Time (s)", "Max slip rate (log10 m/s)", "Moment rate (N-m/s)"],
                    ["t", "max_slip_rate", "moment_rate"], [T_FMT, V_FMT, V_FMT])

    def catalog(self):
        src = self.src("_events.txt")
        header, rows = read_table(src)
        with open(self.dst("catalog.dat"), "w") as f:
            self.head(f, ["Column #1 = Event number",
                          "Column #2 = Event start time (s)",
                          "Column #3 = Event end time (s)",
                          "Column #4 = Rupture area (m^2)",
                          "Column #5 = Stress avg at start time (MPa)",
                          "Column #6 = Stress avg at end time (MPa)",
                          "Column #7 = Avg slip in rupture area (m)",
                          "The line below lists the names of the data fields"])
            f.write("event_no t_start t_end rup_area avg_stress_start avg_stress_end avg_slip\n")
            f.write("# Here are the data\n")
            for r in rows:
                f.write("%d %s %s %s\n" % (int(round(r[0])), T_FMT % r[1], T_FMT % r[2],
                                           " ".join(V_FMT % v for v in r[3:])))
        print("catalog.dat: %d events from %s" % (len(rows), os.path.basename(src)))

    def profiles(self):
        n = 0
        for prof, coord, label in (("depth", "x3", "Depth (m)"), ("strike", "x2", "Strike (m)")):
            for q, what in (("slip_2", "Horizontal slip (Slip 2) (m)"),
                            ("slip_3", "Vertical slip (Slip 3) (m)"),
                            ("stress_2", "Horizontal shear stress (Stress 2) (MPa)"),
                            ("stress_3", "Vertical shear stress (Stress 3) (MPa)")):
                src = self.src("_profile_%s_%s.txt" % (prof, q))
                if not os.path.exists(src):
                    continue
                header, rows = read_table(src)
                first, rest = rows[0], rows[1:]
                if first[0] != 0.0 or first[1] != 0.0:
                    sys.exit("%s: the first row is not 0 0 and the positions" % src)
                rest, dropped = increasing(rest)
                npt = len(first) - 2
                dst = self.dst("%s_%s.dat" % (q, prof))
                with open(dst, "w") as f:
                    self.head(f, ["Row #1 = %s with two zeros first" % label,
                                  "Column #1 = Time (s)",
                                  "Column #2 = Max slip rate (log10 m/s)",
                                  "Columns #3-%d = %s" % (npt + 2, what),
                                  "The line below lists the names of the data fields"])
                    f.write("%s\nt\nmax_slip_rate\n%s\n" % (coord, q))
                    f.write("# Here are the data\n")
                    f.write(" ".join([V_FMT % 0.0, V_FMT % 0.0] +
                                     [V_FMT % v for v in first[2:]]) + "\n")
                    for r in rest:
                        f.write(" ".join([T_FMT % r[0], V_FMT % r[1]] +
                                         [V_FMT % v for v in r[2:]]) + "\n")
                print("%s: %d rows x %d points (%d repeated times dropped)" %
                      (os.path.basename(dst), len(rest), npt, dropped))
                n += 1
        return n

    def rupture(self):
        src = self.src("_rupture_event%d.txt" % self.a.event)
        header, rows = read_table(src)
        tol = 1e-6 * (abs(self.a.region_x2[0]) + abs(self.a.region_x2[1]) +
                      abs(self.a.region_x3[0]) + abs(self.a.region_x3[1]) + 1.0)
        best = {}
        for r in rows:
            x2, x3, t = r[3], r[4], r[5]
            if not (self.a.region_x2[0] - tol <= x2 <= self.a.region_x2[1] + tol and
                    self.a.region_x3[0] - tol <= x3 <= self.a.region_x3[1] + tol):
                continue
            key = (round(x2, 3), round(x3, 3))
            if key not in best or t < best[key][2]:
                best[key] = (x2, x3, t)
        pts = sorted(best.values(), key=lambda p: (p[1], p[0]))
        with open(self.dst("rupture.dat"), "w") as f:
            self.head(f, ["Event %d; first time slip rate >= the run's rupture_V "
                          "(1 mm/s); 1.0E+09 = never" % self.a.event,
                          "Column #1 = x2 (m)", "Column #2 = x3 (m)", "Column #3 = time (s)",
                          "The line below lists the names of the data fields"])
            f.write("x2 x3 t\n")
            f.write("# Here are the data\n")
            for x2, x3, t in pts:
                f.write("%s %s %s\n" % (T_FMT % x2, T_FMT % x3, T_FMT % t))
        print("rupture.dat: %d positions (of %d fault DOF rows) from %s" %
              (len(pts), len(rows), os.path.basename(src)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_dir")
    p.add_argument("prefix")
    p.add_argument("out_dir")
    p.add_argument("--modeler", required=True)
    p.add_argument("--element-size", required=True, help='e.g. "1000 m"')
    p.add_argument("--order", type=int, required=True, help="polynomial degree")
    p.add_argument("--problem", default="SEAS Benchmark BP5-FD")
    p.add_argument("--code", default="seas-mfem (seas_fd_cycle_driver)")
    p.add_argument("--version", default="", help="default: the commit in the run's headers")
    p.add_argument("--event", type=int, default=1, help="event of rupture.dat")
    p.add_argument("--region-x2", type=float, nargs=2, default=[-34000.0, 34000.0])
    p.add_argument("--region-x3", type=float, nargs=2, default=[0.0, 20000.0])
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    e = Exporter(a)
    if e.on_fault() == 0:
        sys.exit("no %s_fltst_*.txt in %s" % (a.prefix, a.run_dir))
    e.off_fault()
    e.source()
    e.catalog()
    e.profiles()
    if os.path.exists(e.src("_rupture_event%d.txt" % a.event)):
        e.rupture()
    else:
        print("rupture.dat: not written (no %s_rupture_event%d.txt yet)" % (a.prefix, a.event))


if __name__ == "__main__":
    main()
