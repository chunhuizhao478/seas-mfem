#!/usr/bin/env python3
"""Validate a BP5 run's station and global files against a golden run.

GPU plan G0 gate (jobs/horizon/check_gate.sh).  seas_bp5_full's built-in
comparator (--regression-tolerance) is not enough on its own: it skips a
missing reference, drops data lines it cannot parse, holds the reference
constant outside its time range, and passes a NaN error.  This script checks
the files themselves first:

  * all 10 BP5 stations (io/bp5_benchmark_output.hpp:821-830) and the global
    file exist in both the run and the golden;
  * every non-comment line parses into the expected number of columns
    (8 for a station, 2 for the global file) and every value is finite;
  * times increase strictly, with two exceptions that the driver's output
    produces: the second row may repeat the first row's time (with
    --petsc-ts the files get the initial write and the step-0 monitor
    write), and the last row may repeat the time of the row before (the
    end-of-run write; its values can differ from the step's row at
    round-off when the operator was evaluated again in between, as
    --check-accepted-state does).  Those rows are compared like any other;
  * start time and row count match the golden's exactly.

End times are compared with a tolerance, not exactly: the step size is
adaptive, so any run that is not bitwise identical (two runs of one binary
on several ranks already differ at round-off) ends at a slightly different
time.  The end time's relative difference is reported with the field errors
(as global / end_time) and is held to the same tolerance.

It then recomputes each field's relative L2 error with the comparator's
formula (tests/verification/bp5_verification_full.cpp, InterpolateOnto and
RelativeL2Error: the golden interpolated linearly onto the run's times, held
constant outside its range, ||run - golden||_2 / ||golden||_2), and the event
count (upward crossings of log10(V_max) = -3 in the global file).

With --tol, every error must be finite and <= TOL and the event counts must
match.  Without --tol the script only reports (used to measure floors).
Identical files give errors of order 1e-17, not 0, from the interpolation's
round-off (the comparator's too); the report also says how many files are
byte-identical, for bitwise checks.

Exit code: 0 pass, 1 fail, 2 usage error.  Standard library only.
"""

import argparse
import json
import math
import os
import sys

# io/bp5_benchmark_output.hpp:821-830
STATIONS = [
    "fltst_strk-36dp+00", "fltst_strk-16dp+00", "fltst_strk+00dp+00",
    "fltst_strk+16dp+00", "fltst_strk+36dp+00", "fltst_strk-24dp+10",
    "fltst_strk-16dp+10", "fltst_strk+00dp+10", "fltst_strk+16dp+10",
    "fltst_strk+00dp+22",
]
STATION_FIELDS = ["slip_strike", "slip_dip", "log10_V_strike", "log10_V_dip",
                  "tau_strike", "tau_dip", "log10_state"]
EVENT_THRESHOLD_LOG10 = -3.0  # CountEvents() default in the driver


class FileProblem(Exception):
    pass


def read_table(path, ncols):
    """Return the rows of a SCEC text file as lists of floats.

    Every line that is not blank and does not start with '#' must hold
    exactly ncols finite numbers.
    """
    if not os.path.isfile(path):
        raise FileProblem(f"missing: {path}")
    rows = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) != ncols:
                raise FileProblem(f"{path}:{lineno}: {len(parts)} columns, "
                                  f"expected {ncols}")
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                raise FileProblem(f"{path}:{lineno}: unparseable line")
            for v in vals:
                if not math.isfinite(v):
                    raise FileProblem(f"{path}:{lineno}: non-finite value")
            rows.append(vals)
    if not rows:
        raise FileProblem(f"{path}: no data rows")
    return rows


def check_times(path, rows):
    """Times increase strictly, except that the second and the last row may
    repeat the time of the row before (see the module docstring)."""
    for i in range(1, len(rows)):
        if rows[i][0] > rows[i - 1][0]:
            continue
        edge_repeat = (i in (1, len(rows) - 1)
                       and rows[i][0] == rows[i - 1][0])
        if not edge_repeat:
            raise FileProblem(f"{path}: time does not increase at data row "
                              f"{i + 1} ({rows[i - 1][0]!r} -> {rows[i][0]!r})")


def interpolate_onto(x_ref, y_ref, x_target):
    """InterpolateOnto() from bp5_verification_full.cpp, line for line."""
    out = [0.0] * len(x_target)
    n_ref = len(x_ref)
    j = 0
    for i, xt in enumerate(x_target):
        if xt <= x_ref[0]:
            out[i] = y_ref[0]
            continue
        if xt >= x_ref[-1]:
            out[i] = y_ref[-1]
            continue
        while j < n_ref - 2 and x_ref[j + 1] < xt:
            j += 1
        frac = (xt - x_ref[j]) / (x_ref[j + 1] - x_ref[j])
        out[i] = y_ref[j] + frac * (y_ref[j + 1] - y_ref[j])
    return out


def relative_l2(sim, ref):
    """RelativeL2Error() from bp5_verification_full.cpp."""
    if len(sim) != len(ref) or not sim:
        return 1e30
    num = math.sqrt(sum((a - b) ** 2 for a, b in zip(sim, ref)))
    den = math.sqrt(sum(b * b for b in ref))
    if den < 1e-30:
        return 0.0 if num < 1e-30 else 1e30
    return num / den


def count_events(log10_vmax, threshold=EVENT_THRESHOLD_LOG10):
    count, above = 0, False
    for v in log10_vmax:
        if v > threshold and not above:
            count += 1
            above = True
        elif v <= threshold:
            above = False
    return count


def compare(run_dir, run_prefix, golden_dir, golden_prefix):
    """Validate all files, then compute every error.

    Returns (problems, result) where problems is a list of strings (empty
    when every file is valid) and result holds the errors.
    """
    problems = []
    result = {"stations": {}, "global": None, "events": None,
              "identical_files": 0, "files": 0}

    def load(kind, name, ncols):
        tables = {}
        for side, d, p in (("run", run_dir, run_prefix),
                           ("golden", golden_dir, golden_prefix)):
            path = os.path.join(d, f"{p}_{name}.txt")
            try:
                rows = read_table(path, ncols)
                check_times(path, rows)
                tables[side] = rows
            except FileProblem as e:
                problems.append(str(e))
        if len(tables) != 2:
            return None
        result["files"] += 1
        with open(os.path.join(run_dir, f"{run_prefix}_{name}.txt"),
                  "rb") as fa, \
             open(os.path.join(golden_dir, f"{golden_prefix}_{name}.txt"),
                  "rb") as fb:
            if fa.read() == fb.read():
                result["identical_files"] += 1
        run, gold = tables["run"], tables["golden"]
        for what, a, b in (("start time", run[0][0], gold[0][0]),
                           ("row count", len(run), len(gold))):
            if a != b:
                problems.append(f"{kind} {name}: {what} {a!r} != golden "
                                f"{b!r}")
        return run, gold

    for st in STATIONS:
        loaded = load("station", st, 8)
        if loaded is None:
            continue
        run, gold = loaded
        t_run = [r[0] for r in run]
        t_gold = [r[0] for r in gold]
        errs = {}
        for k, field in enumerate(STATION_FIELDS, start=1):
            ref = interpolate_onto(t_gold, [r[k] for r in gold], t_run)
            errs[field] = relative_l2([r[k] for r in run], ref)
        result["stations"][st] = errs

    loaded = load("global", "global", 2)
    if loaded is not None:
        run, gold = loaded
        t_run = [r[0] for r in run]
        ref = interpolate_onto([r[0] for r in gold], [r[1] for r in gold],
                               t_run)
        t_end_run, t_end_gold = run[-1][0], gold[-1][0]
        result["global"] = {
            "log10_Vmax": relative_l2([r[1] for r in run], ref),
            "end_time": (abs(t_end_run - t_end_gold) / abs(t_end_gold)
                         if t_end_gold != 0.0 else abs(t_end_run)),
        }
        result["events"] = {"run": count_events([r[1] for r in run]),
                            "golden": count_events([r[1] for r in gold])}
    return problems, result


def all_errors(result):
    for st, errs in result["stations"].items():
        for field, e in errs.items():
            yield st, field, e
    if result["global"] is not None:
        for field, e in result["global"].items():
            yield "global", field, e


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--run-prefix", required=True)
    ap.add_argument("--golden-dir", required=True)
    ap.add_argument("--golden-prefix", required=True)
    ap.add_argument("--tol", type=float, default=None,
                    help="gate tolerance; omit to report errors only")
    ap.add_argument("--json", default=None, help="write the errors here")
    args = ap.parse_args(argv)

    problems, result = compare(args.run_dir, args.run_prefix,
                               args.golden_dir, args.golden_prefix)
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
    ev = result["events"]
    if ev is not None and ev["run"] != ev["golden"]:
        failures.append(f"event count {ev['run']} != golden {ev['golden']}")

    print(f"bp5_gate_check: run {args.run_dir}/{args.run_prefix} vs golden "
          f"{args.golden_dir}/{args.golden_prefix}")
    print(f"  stations compared: {len(result['stations'])} / {len(STATIONS)}"
          f", global: {'yes' if result['global'] is not None else 'no'}")
    if ev is not None:
        print(f"  events: run {ev['run']}, golden {ev['golden']}")
    print(f"  byte-identical files: {result['identical_files']} / "
          f"{result['files']}")
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
