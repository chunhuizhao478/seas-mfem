#!/usr/bin/env python3
"""Unit tests for bp5/bp5_gate_check.py (GPU plan G0 gate).

Run from miniapps/seas:  python3 tests/unit/test_bp5_gate_check.py
Builds a synthetic golden (10 stations + global file, the layout
seas_bp5_full writes with --petsc-ts --write-every-step) in a temporary
directory and checks that the gate passes, fails or flags each case.
"""
import contextlib
import io
import json
import math
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bp5"))
import bp5_gate_check as gc  # noqa: E402

TIMES = [0.0, 0.0, 0.01, 0.03, 0.06, 0.1, 0.15, 0.21, 0.21]


def station_rows(times, k):
    rows = []
    for t in times:
        rows.append([t, 1e-9 * t * (k + 1), 1e-12 * t, -9.0 + 0.1 * t,
                     -20.0 + t, 13.0 + 0.01 * t + k, 1e-5 * (1 + t), 8.1 + t])
    return rows


def global_rows(times):
    # log10 V_max: starts above the -3 event threshold (BP5 nucleation),
    # so the event count is 1.
    return [[t, -2.0 - t] for t in times[1:]]


def write(path, rows):
    with open(path, "w") as f:
        f.write("# header line\n# t slip ...\n")
        for r in rows:
            f.write(" ".join(f"{v:.12E}" for v in r) + "\n")


def make_run(d, prefix, times=TIMES, edit=None):
    os.makedirs(d, exist_ok=True)
    for k, st in enumerate(gc.STATIONS):
        rows = station_rows(times, k)
        if edit:
            rows = edit(st, rows)
        if rows is not None:
            write(os.path.join(d, f"{prefix}_{st}.txt"), rows)
    rows = global_rows(times)
    if edit:
        rows = edit("global", rows)
    if rows is not None:
        write(os.path.join(d, f"{prefix}_global.txt"), rows)


class GateCheckTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gate_check_test_")
        self.gold = os.path.join(self.tmp, "gold")
        make_run(self.gold, "g")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def gate(self, edit=None, tol=None, times=TIMES):
        run = os.path.join(self.tmp, "run")
        shutil.rmtree(run, ignore_errors=True)
        make_run(run, "r", times=times, edit=edit)
        js = os.path.join(self.tmp, "out.json")
        args = ["--run-dir", run, "--run-prefix", "r", "--golden-dir",
                self.gold, "--golden-prefix", "g", "--json", js]
        if tol is not None:
            args += ["--tol", str(tol)]
        with contextlib.redirect_stdout(io.StringIO()):
            rc = gc.main(args)
        with open(js) as f:
            return rc, json.load(f)

    def test_identical_passes(self):
        rc, d = self.gate(tol=1e-12)
        self.assertEqual(rc, 0)
        self.assertEqual(d["verdict"], "PASS")
        self.assertEqual(d["errors"]["identical_files"], 11)
        self.assertEqual(len(d["errors"]["stations"]), 10)
        self.assertLess(d["max_error"], 1e-15)

    def test_report_only(self):
        rc, d = self.gate()
        self.assertEqual(rc, 0)
        self.assertEqual(d["verdict"], "REPORT ONLY")

    def test_small_change_measured(self):
        def edit(name, rows):
            if name == "fltst_strk+00dp+10":
                for r in rows:
                    r[5] *= 1 + 1e-9
            return rows
        rc, d = self.gate(edit=edit, tol=1e-12)
        self.assertEqual(rc, 1)
        self.assertEqual(d["worst"], ["fltst_strk+00dp+10", "tau_strike"])
        self.assertAlmostEqual(d["max_error"] / 1e-9, 1.0, places=2)
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 0)

    def test_missing_station(self):
        def edit(name, rows):
            return None if name == "fltst_strk-24dp+10" else rows
        rc, d = self.gate(edit=edit)
        self.assertEqual(rc, 1)
        self.assertEqual(d["verdict"], "INVALID FILES")
        self.assertEqual(len(d["errors"]["stations"]), 9)

    def test_missing_global(self):
        rc, d = self.gate(edit=lambda n, r: None if n == "global" else r,
                          tol=1e-6)
        self.assertEqual(rc, 1)

    def test_nan_value(self):
        def edit(name, rows):
            if name == "fltst_strk+16dp+00":
                rows[4][3] = float("nan")
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("non-finite" in f for f in d["failures"]))

    def test_short_line(self):
        def edit(name, rows):
            if name == "fltst_strk+16dp+00":
                rows[3] = rows[3][:7]
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("columns" in f for f in d["failures"]))

    def test_time_goes_back(self):
        def edit(name, rows):
            if name == "global":
                rows[3][0] = rows[2][0]
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("does not increase" in f for f in d["failures"]))

    def test_end_rows_repeat_time(self):
        # The end-of-run row repeats the last step's time; its values may
        # differ at round-off.
        def edit(name, rows):
            if name != "global":
                rows[-1][6] *= 1 + 1e-13
            return rows
        rc, d = self.gate(edit=edit, tol=1e-10)
        self.assertEqual(rc, 0)

    def test_row_count(self):
        times = TIMES[:-2] + [0.18] + TIMES[-2:]
        rc, d = self.gate(times=times, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("row count" in f for f in d["failures"]))

    def test_end_time_tolerance(self):
        times = TIMES[:-2] + [0.21 * (1 + 1e-10)] * 2
        rc, d = self.gate(times=times, tol=1e-12)
        self.assertEqual(rc, 1)
        e = d["errors"]["global"]["end_time"]
        self.assertAlmostEqual(e / 1e-10, 1.0, places=3)
        rc, d = self.gate(times=times, tol=1e-8)
        self.assertEqual(rc, 0)

    def test_event_count(self):
        def edit(name, rows):
            if name == "global":
                rows[3][1] = -4.0   # drops below -3 and rises again
            return rows
        rc, d = self.gate(edit=edit, tol=1.0)
        self.assertEqual(rc, 1)
        self.assertTrue(any("event count" in f for f in d["failures"]))

    def test_interpolation_matches_driver(self):
        # InterpolateOnto: held constant outside the reference range.
        out = gc.interpolate_onto([0.0, 1.0, 2.0], [0.0, 10.0, 0.0],
                                  [-1.0, 0.5, 1.5, 3.0])
        self.assertEqual(out, [0.0, 5.0, 5.0, 0.0])
        self.assertEqual(gc.relative_l2([1.0], [0.0]), 1e30)
        self.assertTrue(math.isclose(gc.relative_l2([3.0, 4.0], [3.0, 3.0]),
                                     1.0 / math.sqrt(18.0)))


if __name__ == "__main__":
    unittest.main(verbosity=1)
