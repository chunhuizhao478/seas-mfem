#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# scec_bp3_export.py — 2D in-plane design rev 3, phase E: the SCEC upload
# files of the BP3-QD description (section 4) from a seas_fd_cycle_driver run
# of a config/fd_cycle/bp3_qd_*.toml deck (station_sigma_n and
# [cycle.benchmark] with the "dip" profile and profile_sigma_n).  Standard
# library only.
#
#   python3 cycle/scec_bp3_export.py RUN_DIR PREFIX OUT_DIR \
#       --modeler "Chunhui Zhao" --element-size "25 m" --order 1 \
#       --scenario "dip 60, thrust" [--problem "SEAS Benchmark BP3-QD"] \
#       [--version COMMIT] [--domain "..."]
#
# Writes into OUT_DIR (made if missing; existing files of these names are
# replaced):
#   fltst_dp*.dat        section 4.1, from RUN_DIR/PREFIX_fltst_dp*.txt
#   srfst_fn*.dat        section 4.2, from RUN_DIR/PREFIX_srfst_fn*.txt
#   slip.dat, shear_stress.dat, normal_stress.dat
#                        section 4.3, from PREFIX_profile_dip_{slip_3,
#                        stress_3,sigma_n}.txt
# and prints one line per file.  Numbers are written as the description
# recommends for C users: 21.13E for time, 14.6E for the rest.
#
# Columns: on the BP3 decks (up = [0, 1], ref_normal = [1, 0]) the kernel's
# dip slot is BP3's delta (eq. 6) and its tau, so the converter takes
#   slip = slip_dip, slip_rate = log10 |V_dip|, shear_stress = tau_dip,
#   normal_stress = sigma_n (compression positive), state = log10 theta
# and refuses a file whose strike slot is not zero (a deck with another
# frame).  The surface stations are written as the run gives them: the
# decks' volume_axes e1 = (1, 0), e2 = (0, -1) make disp_2 and vel_2
# positive downward, as section 4.2 asks.  Rows keep only times larger than
# the previous kept row's (the description wants increasing time).

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
        f.write("# element size=%s\n" % self.a.element_size)
        f.write("# polynomial degree=%d\n" % self.a.order)
        f.write("# scenario=%s\n" % self.a.scenario)
        for e in extra:
            f.write("# %s\n" % e)

    def src(self, name):
        return os.path.join(self.a.run_dir, self.a.prefix + name)

    def dst(self, name):
        return os.path.join(self.a.out_dir, name)

    def write_series(self, src, dst, location, columns, fields, rows, dropped):
        lo, hi = step_stats(rows)
        with open(dst, "w") as f:
            self.head(f, ["location=%s" % location,
                          "minimum time step=%.6E" % lo,
                          "maximum time step=%.6E" % hi,
                          "num time steps=%d" % len(rows)] +
                      ["Column #%d = %s" % (i + 1, c) for i, c in enumerate(columns)] +
                      ["The line below lists the names of the data fields"])
            f.write(" ".join(fields) + "\n")
            f.write("# Here is the time-series data.\n")
            for r in rows:
                f.write(" ".join([T_FMT % r[0]] + [V_FMT % v for v in r[1:]]) + "\n")
        print("%s: %d rows (%d repeated times dropped) from %s" %
              (os.path.basename(dst), len(rows), dropped, os.path.basename(src)))

    def read_rows(self, src, ncol):
        header, rows = read_table(src)
        if not self.version:
            self.version = commit_of(header)
        rows, dropped = increasing(rows)
        if rows and len(rows[0]) != ncol:
            sys.exit("%s: %d columns, expected %d" % (src, len(rows[0]), ncol))
        return rows, dropped

    def on_fault(self):
        files = sorted(glob.glob(self.src("_fltst_dp*.txt")))
        for src in files:
            name = os.path.basename(src)[len(self.a.prefix) + 1:-4]
            m = re.match(r"fltst_dp(\d{3})$", name)
            if not m:
                sys.exit("%s: not a BP3 on-fault station name" % name)
            # t, slip_s, slip_d, log10 V_s, log10 V_d, tau_s, tau_d, log10 theta, sigma_n
            rows, dropped = self.read_rows(src, 9)
            for r in rows:
                if r[1] != 0.0 or r[5] != 0.0:
                    sys.exit("%s: strike slip %g or strike stress %g at t = %g: not a BP3 "
                             "deck's frame (the slip is in the dip slot)" %
                             (src, r[1], r[5], r[0]))
            out = [[r[0], r[2], r[4], r[6], r[8], r[7]] for r in rows]
            loc = "on fault, %gkm down-dip distance" % (int(m.group(1)) / 10.0)
            self.write_series(src, self.dst(name + ".dat"), loc,
                              ["Time (s)", "Slip (m)", "Slip rate (log10 m/s)",
                               "Shear stress (MPa)", "Normal stress (MPa)",
                               "State (log10 s)"],
                              ["t", "slip", "slip_rate", "shear_stress", "normal_stress",
                               "state"], out, dropped)
        return len(files)

    def surface(self):
        files = sorted(glob.glob(self.src("_srfst_fn*.txt")))
        for src in files:
            name = os.path.basename(src)[len(self.a.prefix) + 1:-4]
            m = re.match(r"srfst_fn([+-])(\d{2})$", name)
            if not m:
                sys.exit("%s: not a BP3 surface station name" % name)
            rows, dropped = self.read_rows(src, 5)
            d = int(m.group(2))
            loc = ("on surface, %s%dkm distance off-fault" % (m.group(1), d) if d else
                   "on surface, 0%s km (%s side of the fault trace)" %
                   (m.group(1), "+x" if m.group(1) == "+" else "-x"))
            self.write_series(src, self.dst(name + ".dat"), loc,
                              ["Time (s)", "Displacement 1 (m)", "Displacement 2 (m)",
                               "Velocity 1 (m/s)", "Velocity 2 (m/s)"],
                              ["t", "disp_1", "disp_2", "vel_1", "vel_2"], rows, dropped)
        return len(files)

    def profiles(self):
        n = 0
        zero = self.src("_profile_dip_slip_2.txt")
        if os.path.exists(zero):
            header, rows = read_table(zero)
            for r in rows[1:]:
                if any(v != 0.0 for v in r[2:]):
                    sys.exit("%s: strike slip %g at t = %g: not a BP3 deck's frame" %
                             (zero, max(abs(v) for v in r[2:]), r[0]))
        for q, out, what in (("slip_3", "slip", "Slip (m)"),
                             ("stress_3", "shear_stress", "Shear stress (MPa)"),
                             ("sigma_n", "normal_stress", "Normal stress (MPa)")):
            src = self.src("_profile_dip_%s.txt" % q)
            if not os.path.exists(src):
                sys.exit("%s: missing (the BP3 decks write it)" % src)
            header, rows = read_table(src)
            if not self.version:
                self.version = commit_of(header)
            first, rest = rows[0], rows[1:]
            if first[0] != 0.0 or first[1] != 0.0:
                sys.exit("%s: the first row is not 0 0 and the positions" % src)
            rest, dropped = increasing(rest)
            npt = len(first) - 2
            dst = self.dst(out + ".dat")
            with open(dst, "w") as f:
                f.write("# This is the file header:\n")
                f.write("# problem=%s\n" % self.a.problem)
                f.write("# modeler=%s\n" % self.a.modeler)
                f.write("# date=%s\n" % self.date)
                f.write("# code=%s\n" % self.a.code)
                f.write("# code version=%s\n" % self.version)
                f.write("# element size=%s\n" % self.a.element_size)
                f.write("# polynomial degree=%d\n" % self.a.order)
                f.write("# scenario=%s\n" % self.a.scenario)
                f.write("# Row #1 = Distance down dip (m) with two zeros first\n")
                f.write("# Column #1 = Time (s)\n")
                f.write("# Column #2 = Max slip rate (log10 m/s)\n")
                f.write("# Columns #3-%d = %s\n" % (npt + 2, what))
                f.write("# Computational domain size: %s\n" % self.a.domain)
                f.write("# The line below lists the names of the data fields\n")
                f.write("xd\nt max_slip_rate %s\n" % out)
                f.write("# Here are the data\n")
                f.write(" ".join([V_FMT % 0.0, V_FMT % 0.0] +
                                 [V_FMT % v for v in first[2:]]) + "\n")
                for r in rest:
                    f.write(" ".join([T_FMT % r[0], V_FMT % r[1]] +
                                     [V_FMT % v for v in r[2:]]) + "\n")
            print("%s: %d rows x %d points (%d repeated times dropped) from %s" %
                  (os.path.basename(dst), len(rest), npt, dropped, os.path.basename(src)))
            n += 1
        return n


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_dir")
    p.add_argument("prefix")
    p.add_argument("out_dir")
    p.add_argument("--modeler", required=True)
    p.add_argument("--element-size", required=True, help='e.g. "25 m"')
    p.add_argument("--order", type=int, required=True, help="polynomial degree")
    p.add_argument("--scenario", required=True, help='e.g. "dip 60, thrust"')
    p.add_argument("--problem", default="SEAS Benchmark BP3-QD")
    p.add_argument("--code", default="seas-mfem (seas_fd_cycle_driver)")
    p.add_argument("--version", default="", help="default: the commit in the run's headers")
    p.add_argument("--domain", default="depth 400 km, surface half-width max(400, 400 cot dip) "
                   "km, sides parallel to the fault (bp3/mesh/bp3.geo)")
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    e = Exporter(a)
    if e.on_fault() == 0:
        sys.exit("no %s_fltst_dp*.txt in %s" % (a.prefix, a.run_dir))
    e.surface()
    e.profiles()


if __name__ == "__main__":
    main()
