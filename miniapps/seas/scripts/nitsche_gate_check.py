#!/usr/bin/env python3
"""Validate a Nitsche dynamic-rupture run's station files against a golden run.

GPU dynamic-rupture plan D0 gate (jobs/horizon/check_dyn_gate.sh), the
Nitsche counterpart of bp5/bp5_gate_check.py, whose reader and relative L2
error it imports.  seas_nitsche_dyn_driver writes one file per station,
<prefix>_station_<name>.dat (dynamic/fault_station_file.hpp), with the
columns

    t slip1 slip2 V1 V2 tau1 tau2 sigma_n state

The files are checked first:

  * the run and the golden hold the same, non-empty set of stations (the
    _copies.dat files are not part of the gate), and --expect-stations N,
    when given, is the number of stations;
  * every non-comment line parses into 9 finite numbers;
  * times increase strictly;
  * the row count matches the golden's, and every row's time matches the
    golden's within 1e-12 of the golden's time span (the step size is fixed,
    so the times agree to the printed digits; the largest difference is
    reported).

It then computes each station field's relative L2 error over the rows.  The
two components of slip, of slip rate and of shear traction are measured
against the record of the whole vector at that station, not of the component
alone:

    error(tau1) = ||run.tau1 - golden.tau1||_2 / ||(golden.tau1, golden.tau2)||_2

and likewise for slip1, slip2, V1, V2 and tau2, while sigma_n and state are
measured against their own record, ||run - golden||_2 / ||golden||_2 (the
formula of bp5_gate_check.relative_l2).  A component that is round-off on its
own (the dip shear traction tau1 of TPV205 is 1e-4 to 1e4 Pa against a shear
record of 1e9 Pa) would otherwise fail any gate on a change at round-off.
One more guard: when a station's golden record of the vector (or field) is
below 1e-6 of the largest record over all stations (a station the rupture has
not reached, whose slip and slip rate are zero), the error is taken relative
to 1e-6 of that largest record instead.  Such entries are counted in the
report.

With --tol, every error must be finite and <= TOL.  Without --tol the script
only reports (used to measure floors).  The report also says how many files
are byte-identical, for bitwise checks.  The files hold 11 significant digits,
so differences below about 1e-11 relative to a value do not show.

Exit code: 0 pass, 1 fail, 2 usage error.  Standard library only.
"""

import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bp5"))
import bp5_gate_check as gc  # noqa: E402  (read_table, relative_l2, FileProblem)

FIELDS = ["slip1", "slip2", "V1", "V2", "tau1", "tau2", "sigma_n", "state"]
NCOLS = 1 + len(FIELDS)
# column(s) whose joint record is the denominator of each field's error
GROUPS = {"slip1": (1, 2), "slip2": (1, 2), "V1": (3, 4), "V2": (3, 4),
          "tau1": (5, 6), "tau2": (5, 6), "sigma_n": (7,), "state": (8,)}
NEAR_ZERO = 1e-6     # station record / largest station record of the group
TIME_RTOL = 1e-12    # time difference / golden time span


def station_names(d, prefix):
    """Station names from <prefix>_station_<name>.dat, _copies.dat left out."""
    head = os.path.join(d, f"{prefix}_station_")
    names = []
    for path in glob.glob(glob.escape(head) + "*.dat"):
        name = path[len(head):-len(".dat")]
        if not name.endswith("_copies"):
            names.append(name)
    return sorted(names)


def check_increasing(path, rows):
    for i in range(1, len(rows)):
        if not rows[i][0] > rows[i - 1][0]:
            raise gc.FileProblem(f"{path}: time does not increase at data row "
                                 f"{i + 1} ({rows[i - 1][0]!r} -> "
                                 f"{rows[i][0]!r})")


def l2(x):
    return math.sqrt(sum(v * v for v in x))


def compare(run_dir, run_prefix, golden_dir, golden_prefix,
            expect_stations=None):
    """Validate all files, then compute every error.

    Returns (problems, result): problems lists invalid files and mismatched
    layouts (empty when all is valid), result holds the errors.
    """
    problems = []
    result = {"stations": {}, "identical_files": 0, "files": 0,
              "near_zero_entries": 0, "max_time_diff": 0.0}
    run_names = station_names(run_dir, run_prefix)
    gold_names = station_names(golden_dir, golden_prefix)
    if not gold_names:
        problems.append(f"no station files {golden_dir}/{golden_prefix}"
                        "_station_*.dat in the golden")
    if run_names != gold_names:
        only_run = sorted(set(run_names) - set(gold_names))
        only_gold = sorted(set(gold_names) - set(run_names))
        problems.append(f"station sets differ: only in run {only_run}, "
                        f"only in golden {only_gold}")
    if expect_stations is not None and len(gold_names) != expect_stations:
        problems.append(f"golden has {len(gold_names)} stations, expected "
                        f"{expect_stations}")

    tables = {}
    for name in sorted(set(run_names) & set(gold_names)):
        loaded = {}
        for side, d, p in (("run", run_dir, run_prefix),
                           ("golden", golden_dir, golden_prefix)):
            path = os.path.join(d, f"{p}_station_{name}.dat")
            try:
                rows = gc.read_table(path, NCOLS)
                check_increasing(path, rows)
                loaded[side] = rows
            except gc.FileProblem as e:
                problems.append(str(e))
        if len(loaded) != 2:
            continue
        result["files"] += 1
        with open(os.path.join(run_dir, f"{run_prefix}_station_{name}.dat"),
                  "rb") as fa, \
             open(os.path.join(golden_dir,
                               f"{golden_prefix}_station_{name}.dat"),
                  "rb") as fb:
            if fa.read() == fb.read():
                result["identical_files"] += 1
        run, gold = loaded["run"], loaded["golden"]
        if len(run) != len(gold):
            problems.append(f"station {name}: row count {len(run)} != golden "
                            f"{len(gold)}")
            continue
        span = gold[-1][0] - gold[0][0]
        dt_max = max(abs(a[0] - b[0]) for a, b in zip(run, gold))
        result["max_time_diff"] = max(result["max_time_diff"], dt_max)
        if dt_max > TIME_RTOL * max(span, abs(gold[-1][0])):
            problems.append(f"station {name}: times differ from the golden by "
                            f"up to {dt_max!r}")
            continue
        tables[name] = (run, gold)

    def record(gold, cols):
        return math.sqrt(sum(r[c] * r[c] for r in gold for c in cols))

    scale = {cols: max((record(gold, cols) for _, gold in tables.values()),
                       default=0.0) for cols in set(GROUPS.values())}
    for name, (run, gold) in tables.items():
        errs = {}
        for k, field in enumerate(FIELDS, start=1):
            sim = [r[k] for r in run]
            ref = [r[k] for r in gold]
            cols = GROUPS[field]
            floor = NEAR_ZERO * scale[cols]
            den = record(gold, cols)
            if den >= floor and cols == (k,):
                errs[field] = gc.relative_l2(sim, ref)
                continue
            if den < floor:
                result["near_zero_entries"] += 1
                den = floor
            num = l2([a - b for a, b in zip(sim, ref)])
            if den < 1e-30:   # the group is zero at every station
                errs[field] = 0.0 if num < 1e-30 else 1e30
            else:
                errs[field] = num / den
        result["stations"][name] = errs
    return problems, result


def all_errors(result):
    for st, errs in result["stations"].items():
        for field, e in errs.items():
            yield st, field, e


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--run-prefix", required=True)
    ap.add_argument("--golden-dir", required=True)
    ap.add_argument("--golden-prefix", required=True)
    ap.add_argument("--tol", type=float, default=None,
                    help="gate tolerance; omit to report errors only")
    ap.add_argument("--expect-stations", type=int, default=None,
                    help="number of stations the golden must hold")
    ap.add_argument("--json", default=None, help="write the errors here")
    args = ap.parse_args(argv)

    problems, result = compare(args.run_dir, args.run_prefix,
                               args.golden_dir, args.golden_prefix,
                               args.expect_stations)
    errors = list(all_errors(result))
    finite = [e for _, _, e in errors if math.isfinite(e)]
    max_err = max(finite) if finite else float("nan")
    worst = max(errors, key=lambda x: x[2] if math.isfinite(x[2])
                else float("inf")) if errors else None

    failures = list(problems)
    for st, field, e in errors:
        if not math.isfinite(e):
            failures.append(f"{st} {field}: non-finite error {e!r}")
        elif args.tol is not None and e > args.tol:
            failures.append(f"{st} {field}: error {e:.6e} > tol "
                            f"{args.tol:.6e}")
    if not errors and not failures:
        failures.append("no station was compared")

    print(f"nitsche_gate_check: run {args.run_dir}/{args.run_prefix} vs "
          f"golden {args.golden_dir}/{args.golden_prefix}")
    print(f"  stations compared: {len(result['stations'])}, fields: "
          f"{len(FIELDS)}, near-zero records: {result['near_zero_entries']}")
    print(f"  byte-identical files: {result['identical_files']} / "
          f"{result['files']}")
    print(f"  largest time difference: {result['max_time_diff']:.3e} s")
    if worst is not None:
        print(f"  max relative L2 error: {max_err:.6e} ({worst[0]} / "
              f"{worst[1]})")
    if args.tol is not None:
        print(f"  tolerance: {args.tol:.6e}")
    for f in failures:
        print(f"  FAIL: {f}")
    if args.tol is None:
        verdict = "INVALID FILES" if failures else "REPORT ONLY"
    else:
        verdict = "FAIL" if failures else "PASS"
    print(f"  GATE: {verdict}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump({"run_dir": args.run_dir, "run_prefix": args.run_prefix,
                       "golden_dir": args.golden_dir,
                       "golden_prefix": args.golden_prefix,
                       "tol": args.tol, "max_error": max_err,
                       "worst": list(worst[:2]) if worst else None,
                       "failures": failures, "verdict": verdict,
                       "errors": result}, f, indent=1)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
