#!/usr/bin/env python3
"""Turn G0's two floor measurements at one order and rank count into the
gate tolerances of the GPU plan.

floor_N is the larger of bp5_gate_check.py's max error for the repeat run
and for the perturbed run (each against the golden).  For every phase,
TOL = max(tol_phase, 10 floor_N); a TOL above the 1e-6 ceiling is flagged,
since the plan investigates such a spread before using a gate at N ranks.
Runs on Horizon add G1's Vista-to-Horizon difference to floor_N later.
"""
import argparse
import json
import math
import sys

TOL_PHASE = {
    "G1_host_syncs_cpu": 1e-12,   # G1 host syncs on Vista: no arithmetic change
    "G2_operators_cpu": 1e-10,    # G2: operators instead of face loops
    "G1_first_horizon_run": 1e-8,  # plus G3-G5 (before the Horizon floor is added)
    "G3_to_G5_gpu": 1e-8,
}
CEILING = 1e-6

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--order", type=int, required=True)
ap.add_argument("--ranks", type=int, required=True)
ap.add_argument("--repeat", required=True)
ap.add_argument("--perturbed", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

rep = json.load(open(a.repeat))
per = json.load(open(a.perturbed))
problems = [f"repeat: {p}" for p in rep["failures"]] + \
           [f"perturbed: {p}" for p in per["failures"]]
floor = max(rep["max_error"], per["max_error"])
tol = {k: max(v, 10.0 * floor) for k, v in TOL_PHASE.items()}
out = {
    "order": a.order, "ranks": a.ranks,
    "repeat_max_error": rep["max_error"], "repeat_worst": rep["worst"],
    "repeat_identical_files": rep["errors"]["identical_files"],
    "perturbed_max_error": per["max_error"], "perturbed_worst": per["worst"],
    "floor": floor, "tol": tol, "ceiling": CEILING,
    "above_ceiling": 10.0 * floor > CEILING,
    "problems": problems,
}
json.dump(out, open(a.out, "w"), indent=1)
print(f"G0 floor p={a.order} ranks={a.ranks}: repeat {rep['max_error']:.3e} "
      f"({rep['errors']['identical_files']}/{rep['errors']['files']} files "
      f"byte-identical), perturbed {per['max_error']:.3e}, floor {floor:.3e}")
for k, v in tol.items():
    print(f"  TOL {k}: {v:.3e}")
if out["above_ceiling"]:
    print(f"  WARNING: 10 x floor = {10 * floor:.3e} exceeds the {CEILING:.0e} "
          "ceiling; investigate before using a gate at this rank count")
for p in problems:
    print(f"  PROBLEM: {p}")
sys.exit(1 if problems or not math.isfinite(floor) else 0)
