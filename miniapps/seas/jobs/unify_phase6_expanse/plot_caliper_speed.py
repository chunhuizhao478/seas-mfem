#!/usr/bin/env python3
"""Speed-vs-polynomial-order graph from the Phase-C Caliper region reports.

Usage:  python3 plot_caliper_speed.py [report_dir]

Scans `report_dir` (default: this script's directory) for the Phase-C
outputs named `{tpv205|tpv102}_p{1|2|3}_cali_<jobid>.cali-region-report.txt`
(the runtime-report channel written by the
`*_npsweep_cali_expanse.sbatch` jobs), extracts per-region average and
maximum time per rank, and writes:

  caliper_speed_vs_order.csv   (problem, P, region, avg_s, max_s)
  caliper_speed_vs_order.png   one panel per problem: avg time/rank vs P
                               for the headline regions (log-y), with the
                               Max/Avg imbalance shown as error caps.

runtime-report table format (calc.inclusive,region.count,profile.mpi):
indented region-path column followed by numeric columns; the first three
numeric columns are Min/Max/Avg time per rank.  Parsing is tolerant: rows
that do not match are skipped.

Regions plotted (when present):
  seas::WaveOperator::AdvanceADER            whole ADER step
  seas::WaveOperator::ComputeADERFaceFluxRHS interior face flux (corrector)
  seas::WaveOperator::ComputeADERSharedFaceFluxRHS  seam face flux
  seas::WaveOperator::ComputeVolumeRHS       volume derivative
  seas::spatial_dyn::friction_substep        friction iterator (when marked)

Standard library + numpy + matplotlib.
"""

import csv
import glob
import os
import re
import sys

import numpy as np

HEADLINE = [
    "AdvanceADER",
    "ComputeADERFaceFluxRHS",
    "ComputeADERSharedFaceFluxRHS",
    "ComputeVolumeRHS",
    "friction_substep",
]

FNAME_RE = re.compile(
    r"(?P<prob>tpv205|tpv102)_p(?P<p>[123])_cali_\d+\.cali-region-report\.txt$")


def parse_report(path):
    """Return {region_leaf: (avg_s, max_s)} from a runtime-report table."""
    out = {}
    for line in open(path, errors="replace"):
        # region rows: possibly-indented path then >=3 float columns
        m = re.match(r"^(\s*)([A-Za-z_][\w:<>~ .\-]*?)\s+((?:[0-9.eE+\-]+\s+)+)",
                     line.rstrip("\n"))
        if not m:
            continue
        nums = m.group(3).split()
        try:
            vals = [float(x) for x in nums]
        except ValueError:
            continue
        if len(vals) < 3:
            continue
        leaf = m.group(2).strip().split("::")[-1].split("/")[-1].strip()
        mn, mx, avg = vals[0], vals[1], vals[2]
        # keep the LARGEST entry per leaf (a leaf can appear at several
        # tree depths; the headline number is the dominant one)
        if leaf not in out or avg > out[leaf][0]:
            out[leaf] = (avg, mx)
    return out


def main(argv):
    d = argv[1] if len(argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    data = {}   # (prob, P) -> {leaf: (avg, max)}
    for f in sorted(glob.glob(os.path.join(d, "*.cali-region-report.txt"))):
        m = FNAME_RE.search(os.path.basename(f))
        if not m:
            continue
        key = (m.group("prob"), int(m.group("p")))
        data[key] = parse_report(f)
        print(f"parsed {os.path.basename(f)}: {len(data[key])} regions")
    if not data:
        print(f"no {{prob}}_p{{P}}_cali_*.cali-region-report.txt files under {d}")
        return 1

    csv_path = os.path.join(d, "caliper_speed_vs_order.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["problem", "P", "region", "avg_s_per_rank",
                    "max_s_per_rank"])
        for (prob, p) in sorted(data):
            for leaf, (avg, mx) in sorted(data[(prob, p)].items()):
                w.writerow([prob, p, leaf, f"{avg:.6g}", f"{mx:.6g}"])
    print("wrote", csv_path)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable — CSV written, skipping the PNG.")
        return 0

    probs = sorted({k[0] for k in data})
    fig, axes = plt.subplots(1, len(probs), figsize=(6.5 * len(probs), 5),
                             squeeze=False)
    for ax, prob in zip(axes[0], probs):
        ps = sorted(p for (pr, p) in data if pr == prob)
        for leaf in HEADLINE:
            avg = [data[(prob, p)].get(leaf, (np.nan, np.nan))[0] for p in ps]
            mx = [data[(prob, p)].get(leaf, (np.nan, np.nan))[1] for p in ps]
            if all(np.isnan(avg)):
                continue
            avg = np.array(avg); mx = np.array(mx)
            ax.errorbar(ps, avg, yerr=[np.zeros_like(avg), mx - avg],
                        marker="o", capsize=4, label=leaf)
        ax.set_xlabel("polynomial order P (ADER O = P+1)")
        ax.set_ylabel("time per rank [s]  (cap = max rank)")
        ax.set_title(f"{prob} unified substep — Caliper speed vs order")
        ax.set_xticks(ps)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(d, "caliper_speed_vs_order.png")
    fig.savefig(png, dpi=150)
    print("wrote", png)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
