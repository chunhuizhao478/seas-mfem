#!/usr/bin/env python3
"""Summarize the QD speed timing runs of one or more job folders.

    python3 jobs/vista_qd_speed/qd_speed_summary.py DIR [DIR ...]

DIR is a folder written by jobs/vista_qd_speed/qd_speed_timing.sbatch
(bp5/results_qd_speed/timing_nN_rR_cC_pP_jobJ).  For every runK_MODE
subfolder the end-of-run [STAGE-TIMERS] table (stage_timers_final.txt, from
seas_bp5_full --stage-timers) gives the loop wall per accepted step and each
part's time; gemv_bw.txt gives the layout's bandwidth.  One line per run:

  layout, mode, exit code, s/step, stages/step, then per part the slowest
  rank's ms per call (solve, slip_rhs + dirichlet_rhs, traction, friction)
  and the share of the loop wall, matrix_setup (assembly + factorization)
  and the GEMV bandwidth per node.

Times are the slowest rank's (the table's "max s" column), as the table
prints them.  A run without a table (stopped early) is listed with its exit
code only.
"""
import os
import re
import sys

HEAD = re.compile(r"\[STAGE-TIMERS\] end of run: (\d+) accepted steps, (\d+) "
                  r"stages, loop wall ([0-9.]+) s \(([0-9.eE+-]+) s/step\), "
                  r"(\d+) ranks")
SETUP = re.compile(r"\[STAGE-TIMERS\] matrix_setup .*: max ([0-9.]+) s")
ROW = re.compile(r"\[STAGE-TIMERS\] (\S+)\s+(\d+)\s+([0-9.eE+-]+)\s+"
                 r"([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+"
                 r"([0-9.]+)$")
GEMV = re.compile(r"GEMV fp64 .* ([0-9.]+) GB/s per node")


def parse_table(path):
    t = {"parts": {}}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            m = HEAD.search(line)
            if m:
                t["steps"] = int(m.group(1))
                t["stages"] = int(m.group(2))
                t["wall"] = float(m.group(3))
                t["s_per_step"] = float(m.group(4))
                t["ranks"] = int(m.group(5))
                continue
            m = SETUP.search(line)
            if m:
                t["setup"] = float(m.group(1))
                continue
            m = ROW.search(line)
            if m:
                t["parts"][m.group(1)] = {
                    "calls": int(m.group(2)), "max_s": float(m.group(3)),
                    "ms_call": float(m.group(5)), "share": float(m.group(7))}
    return t if "steps" in t else None


def run_dirs(root):
    out = []
    for name in os.listdir(root):
        m = re.fullmatch(r"run(\d+)_(\w+)", name)
        if m and os.path.isdir(os.path.join(root, name)):
            out.append((int(m.group(1)), m.group(2), os.path.join(root, name)))
    return sorted(out)


def main(argv):
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    cols = ("layout", "mode", "exit", "s/step", "st/step", "solve ms",
            "rhs ms", "trac ms", "fric ms", "solve %", "face %", "setup s",
            "GB/s/node")
    print("  ".join(f"{c:>10}" for c in cols))
    for root in argv:
        layout = os.path.basename(os.path.normpath(root))
        m = re.match(r"timing_(n\d+_r\d+_c\d+)", layout)
        layout = m.group(1) if m else layout
        bw = ""
        gpath = os.path.join(root, "gemv_bw.txt")
        if os.path.exists(gpath):
            m = GEMV.search(open(gpath).read())
            bw = m.group(1) if m else ""
        for _, mode, d in run_dirs(root):
            ec = "?"
            if os.path.exists(os.path.join(d, "exit_code")):
                ec = open(os.path.join(d, "exit_code")).read().strip()
            tpath = os.path.join(d, "stage_timers_final.txt")
            t = parse_table(tpath) if os.path.exists(tpath) else None
            if t is None:
                print("  ".join(f"{v:>10}" for v in (layout, mode, ec)))
                continue
            p = t["parts"]

            def ms(name):
                return p[name]["ms_call"] if name in p else 0.0

            def share(*names):
                return sum(p[n]["share"] for n in names if n in p)
            # The matrix path times B*slip and D*g inside slip_rhs and
            # dirichlet_rhs; the loops time them there too.
            rhs = ms("slip_rhs") + ms("dirichlet_rhs")
            vals = (layout, mode, ec, f"{t['s_per_step']:.4f}",
                    f"{t['stages'] / max(t['steps'], 1):.2f}",
                    f"{ms('solve'):.2f}", f"{rhs:.2f}", f"{ms('traction'):.2f}",
                    f"{ms('friction'):.2f}", f"{share('solve'):.1f}",
                    f"{share('slip_expand', 'slip_rhs', 'dirichlet_rhs', 'traction'):.1f}",
                    f"{t.get('setup', 0.0):.1f}", bw)
            print("  ".join(f"{v:>10}" for v in vals))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
