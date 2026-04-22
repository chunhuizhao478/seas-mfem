#!/usr/bin/env python3
"""
TPV102 v9.2.0 plan §4.5 (RANK-1 diagnostic) — sigma_n rank consistency.

Drives `seas_tpv102_driver` at 1-rank and 4-rank decompositions with
IDENTICAL inputs (mesh, order, tfinal), then compares sigma_n(t) at
every fault station at a fixed sampling time.  Classifies H-V92-M
(MPI partition-boundary bit-disagreement on shared fault faces) per
the decision matrix in plan §4.5.

Usage
-----
    python3 test_sigma_n_rank_consistency.py \
        --driver /path/to/seas_tpv102_driver \
        --mesh   tpv102/mesh/tpv102_1000m.msh \
        --tfinal 4.0 \
        --order  1 \
        --cfl    0.5 \
        --output-root tpv102/r_v92_r4_rank_consistency

Guarded defaults
----------------
Default mesh is the 1000 m local mesh, NOT the 200 m production mesh.
Feedback memory "TPV102 no local reproducer runs on production mesh"
forbids running the 200 m driver locally.  The rank-consistency test
can still be informative at 1000 m (partition cuts still intersect
the fault, just at fewer places); if you intend to use the 200 m mesh
you must invoke this script from a Frontera sbatch per user directive
(and override `--mesh`).

Output
------
Writes `sigma_n_rank_consistency_report.json` to the --output-root.
Exit code 0 on PASS (ranks agree within --tol MPa), 1 on FAIL.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults — match plan §4.5
# ---------------------------------------------------------------------------
DEFAULT_MESH    = "tpv102/mesh/tpv102_1000m.msh"
DEFAULT_TFINAL  = 4.0
DEFAULT_ORDER   = 1
DEFAULT_CFL     = 0.5
DEFAULT_TCHECK  = [2.0, 3.5]        # match plan §4.5 sampling points
DEFAULT_TOL_MPA = 0.1                # 100 kPa tolerance per plan §4.5

STATIONS = [
    "flt_0_3", "flt_0_7.5", "flt_0_12",
    "flt_9_7.5", "flt_n9_7.5",
    "flt_12_3", "flt_12_12",
    "flt_n12_3", "flt_n12_12",
]

# Column indices in `tpv_station_*.dat` (0-based, from driver output header):
#   0 time, 1 slip1, 2 slip2, 3 V1, 4 V2, 5 tau1, 6 tau2, 7 sigma_n, 8 log10_theta
COL_TIME    = 0
COL_SIGMA_N = 7


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--driver",      default="./seas_tpv102_driver",
                   help="Path to compiled seas_tpv102_driver binary.")
    p.add_argument("--mesh",        default=DEFAULT_MESH)
    p.add_argument("--tfinal",      type=float, default=DEFAULT_TFINAL)
    p.add_argument("--order",       type=int,   default=DEFAULT_ORDER)
    p.add_argument("--cfl",         type=float, default=DEFAULT_CFL)
    p.add_argument("--tol-mpa",     type=float, default=DEFAULT_TOL_MPA,
                   help="Tolerance on |sigma_n(1-rank) - sigma_n(4-rank)| in MPa.")
    p.add_argument("--t-check",     type=float, nargs="+", default=DEFAULT_TCHECK,
                   help="Sampling times (s) at which to compare sigma_n.")
    p.add_argument("--mpiexec",     default="mpirun",
                   help="MPI launcher command (mpirun / srun).")
    p.add_argument("--output-root", default="tpv102/r_v92_r4_rank_consistency")
    p.add_argument("--ranks-b",     type=int, default=4,
                   help="Second rank count (default 4 — do NOT exceed 14 "
                        "per local-oversubscription memo).")
    p.add_argument("--skip-run",    action="store_true",
                   help="Skip the driver invocations (useful for unit-testing "
                        "the post-processing on existing run directories).")
    p.add_argument("--run-a",       default=None,
                   help="Override directory for run A (low-rank / baseline). "
                        "When set, implies --skip-run for run A.")
    p.add_argument("--run-b",       default=None,
                   help="Override directory for run B (high-rank / production). "
                        "When set, implies --skip-run for run B.")
    return p.parse_args()


# ---------------------------------------------------------------------------
def run_driver(args, ranks, run_dir):
    """Run seas_tpv102_driver with `ranks` MPI ranks."""
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = []
    if ranks > 1:
        cmd.extend([args.mpiexec, "-np", str(ranks)])
    cmd.extend([
        args.driver,
        "--mesh",   args.mesh,
        "--order",  str(args.order),
        "--tfinal", str(args.tfinal),
        "--cfl",    str(args.cfl),
        "--output-dir",    str(run_dir),
        "--output-prefix", "tpv",
    ])
    print(f"[{ranks} rank(s)] {' '.join(cmd)}")
    log_path = run_dir / "run.log"
    with open(log_path, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              check=False)
    if proc.returncode != 0:
        print(f"  driver failed with rc={proc.returncode}")
        print(f"  tail of {log_path}:")
        try:
            with open(log_path) as fh:
                print("\n".join(fh.readlines()[-40:]))
        except OSError:
            pass
        raise SystemExit(1)


def load_station(run_dir, station):
    """Read (time, sigma_n) from tpv_station_<station>.dat."""
    fpath = run_dir / f"tpv_station_{station}.dat"
    times, sig = [], []
    if not fpath.is_file():
        return times, sig
    with open(fpath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) <= COL_SIGMA_N:
                continue
            times.append(float(parts[COL_TIME]))
            sig.append(float(parts[COL_SIGMA_N]))
    return times, sig


def interpolate(times, vals, t):
    """Linear interpolation of `vals` at `t` (times must be sorted)."""
    if not times:
        return None
    if t <= times[0]:  return vals[0]
    if t >= times[-1]: return vals[-1]
    # Bisect
    lo, hi = 0, len(times) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if times[mid] <= t: lo = mid
        else:               hi = mid
    frac = (t - times[lo]) / (times[hi] - times[lo])
    return vals[lo] + frac * (vals[hi] - vals[lo])


# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    root_dir = Path(args.output_root)
    run_a    = Path(args.run_a) if args.run_a else root_dir / "rank1"
    run_b    = Path(args.run_b) if args.run_b else root_dir / f"rank{args.ranks_b}"

    if not args.skip_run and args.run_a is None:
        run_driver(args, 1, run_a)
    if not args.skip_run and args.run_b is None:
        run_driver(args, args.ranks_b, run_b)

    # Compare sigma_n at each station at each t-check value.
    report = {
        "mesh":          args.mesh,
        "tfinal":        args.tfinal,
        "order":         args.order,
        "cfl":           args.cfl,
        "tol_mpa":       args.tol_mpa,
        "t_check":       args.t_check,
        "ranks_a":       1,
        "ranks_b":       args.ranks_b,
        "stations":      {},
        "worst_diff_mpa": 0.0,
        "classification": None,
    }

    worst = 0.0
    worst_at = {}
    for st in STATIONS:
        times_a, sig_a = load_station(run_a, st)
        times_b, sig_b = load_station(run_b, st)
        st_entry = {"t_check": {}}
        for t in args.t_check:
            sa = interpolate(times_a, sig_a, t)
            sb = interpolate(times_b, sig_b, t)
            if sa is None or sb is None:
                st_entry["t_check"][str(t)] = {
                    "sigma_n_rank1":  sa,
                    "sigma_n_rank":   sb,
                    "diff_mpa":       None,
                    "note":           "missing station data",
                }
                continue
            diff_mpa = (sa - sb) * 1e-6
            st_entry["t_check"][str(t)] = {
                "sigma_n_rank1":  sa,
                "sigma_n_rank":   sb,
                "diff_mpa":       diff_mpa,
            }
            if abs(diff_mpa) > worst:
                worst = abs(diff_mpa)
                worst_at = {"station": st, "t": t}
        report["stations"][st] = st_entry

    report["worst_diff_mpa"] = worst
    report["worst_location"] = worst_at

    # Plan §4.5 decision matrix ------------------------------------------
    #   |Δσ_n| <= tol_mpa at all stations → ranks agree (H-V92-M eliminated
    #     at this mesh / tfinal; may need longer tfinal or production
    #     mesh to trigger partition-dependent bug).
    #   |Δσ_n| >  tol_mpa at any station → H-V92-M CONFIRMED: shared-fault
    #     ctor bit disagreement on `elem1_on_plus` or canonical normal.
    if worst <= args.tol_mpa:
        report["classification"] = "RANKS_AGREE"
        verdict = "PASS"
    else:
        report["classification"] = "RANKS_DIFFER_H_V92_M_CONFIRMED"
        verdict = "FAIL"

    report_path = root_dir / "sigma_n_rank_consistency_report.json"
    root_dir.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report: {report_path}")
    print(f"Worst diff: {worst:.4f} MPa at {worst_at}")
    print(f"Verdict  : {verdict}")

    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
