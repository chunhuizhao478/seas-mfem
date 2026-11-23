#!/usr/bin/env python3
"""Unit tests for scripts/nitsche_gate_check.py (GPU dynamic-rupture plan D0).

Run from miniapps/seas:  python3 tests/unit/test_nitsche_gate_check.py
Builds a synthetic golden (station files in the layout seas_nitsche_dyn_driver
writes, dynamic/fault_station_file.hpp) in a temporary directory and checks
that the gate passes, fails or flags each case.
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import nitsche_gate_check as nc  # noqa: E402

STATIONS = ["st000dp075", "st045dp075", "st120dp030", "far"]
DT = 0.005
NROWS = 41


def station_rows(name, k, nrows=NROWS):
    rows = []
    for i in range(nrows):
        t = i * DT
        if name == "far":
            # not reached: zero slip and slip rate, round-off V
            slip, v = 0.0, 1e-22 * (i % 3)
        else:
            slip, v = 0.1 * t * (k + 1), 0.5 + t
        rows.append([t, slip, 0.01 * slip, v, 0.01 * v, 7.0e7 - 1e5 * t,
                     1e5 * k, -1.2e8, slip])
    return rows


def write(path, rows):
    with open(path, "w") as f:
        f.write("# fault station ...\n# Columns: time slip1 slip2 V1 V2 "
                "tau1 tau2 sigma_n slip_path\n# 1 = dip, 2 = strike\n")
        for r in rows:
            f.write(" ".join(f"{v:.10e}" for v in r) + "\n")


def make_run(d, prefix, edit=None, nrows=NROWS):
    os.makedirs(d, exist_ok=True)
    for k, st in enumerate(STATIONS):
        rows = station_rows(st, k, nrows)
        if edit:
            rows = edit(st, rows)
        if rows is not None:
            write(os.path.join(d, f"{prefix}_station_{st}.dat"), rows)
    # a copies file, which the gate leaves out
    write(os.path.join(d, f"{prefix}_station_{STATIONS[0]}_copies.dat"),
          [[0.0, 1.0, 2.0]])


class NitscheGateCheckTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="nitsche_gate_check_test_")
        self.gold = os.path.join(self.tmp, "gold")
        make_run(self.gold, "g")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def gate(self, edit=None, tol=None, nrows=NROWS, expect=None):
        run = os.path.join(self.tmp, "run")
        shutil.rmtree(run, ignore_errors=True)
        make_run(run, "r", edit=edit, nrows=nrows)
        js = os.path.join(self.tmp, "out.json")
        args = ["--run-dir", run, "--run-prefix", "r", "--golden-dir",
                self.gold, "--golden-prefix", "g", "--json", js]
        if tol is not None:
            args += ["--tol", str(tol)]
        if expect is not None:
            args += ["--expect-stations", str(expect)]
        with contextlib.redirect_stdout(io.StringIO()):
            rc = nc.main(args)
        with open(js) as f:
            return rc, json.load(f)

    def test_identical_passes(self):
        rc, d = self.gate(tol=1e-12, expect=4)
        self.assertEqual(rc, 0)
        self.assertEqual(d["verdict"], "PASS")
        self.assertEqual(d["errors"]["identical_files"], 4)
        self.assertEqual(sorted(d["errors"]["stations"]), sorted(STATIONS))
        self.assertEqual(d["max_error"], 0.0)

    def test_report_only(self):
        rc, d = self.gate()
        self.assertEqual(rc, 0)
        self.assertEqual(d["verdict"], "REPORT ONLY")

    def test_expect_stations(self):
        rc, d = self.gate(tol=1e-6, expect=5)
        self.assertEqual(rc, 1)
        self.assertTrue(any("expected 5" in f for f in d["failures"]))

    def test_small_change_measured(self):
        def edit(name, rows):
            if name == "st045dp075":
                for r in rows:
                    r[5] *= 1 + 1e-8
            return rows
        rc, d = self.gate(edit=edit, tol=1e-10)
        self.assertEqual(rc, 1)
        self.assertEqual(d["worst"], ["st045dp075", "tau1"])
        self.assertAlmostEqual(d["max_error"] / 1e-8, 1.0, places=2)
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 0)

    def test_small_component_against_vector(self):
        # tau2 is 1e5 Pa against tau1 = 7e7 Pa: a relative 1e-8 change of
        # tau2 alone is 1.4e-11 of the shear record
        def edit(name, rows):
            if name == "st045dp075":
                for r in rows:
                    r[6] *= 1 + 1e-8
            return rows
        rc, d = self.gate(edit=edit, tol=1e-10)
        self.assertEqual(rc, 0)
        e = d["errors"]["stations"]["st045dp075"]["tau2"]
        self.assertAlmostEqual(e / (1e-8 * 1e5 / 7e7), 1.0, places=2)

    def test_near_zero_record(self):
        # round-off on a station the rupture has not reached: relative to
        # its own record the error is O(1); relative to 1e-6 of the field's
        # largest record it is tiny
        def edit(name, rows):
            if name == "far":
                for r in rows:
                    r[3] = 3e-22
            return rows
        rc, d = self.gate(edit=edit, tol=1e-10)
        self.assertEqual(rc, 0)
        self.assertGreater(d["errors"]["near_zero_entries"], 0)
        self.assertLess(d["errors"]["stations"]["far"]["V1"], 1e-10)

    def test_zero_field_everywhere(self):
        # sigma_n is zero at every station in the golden; any change fails
        def edit(name, rows):
            for r in rows:
                r[7] = 0.0 if name != "far" else 1e-3
            return rows
        os.makedirs(self.gold, exist_ok=True)
        make_run(self.gold, "g", edit=lambda n, r: [x[:7] + [0.0] + x[8:]
                                                     for x in r])
        rc, d = self.gate(edit=edit, tol=1.0)
        self.assertEqual(rc, 1)
        self.assertEqual(d["errors"]["stations"]["far"]["sigma_n"], 1e30)

    def test_missing_station(self):
        rc, d = self.gate(edit=lambda n, r: None if n == "st120dp030" else r)
        self.assertEqual(rc, 1)
        self.assertEqual(d["verdict"], "INVALID FILES")
        self.assertTrue(any("station sets differ" in f for f in d["failures"]))
        self.assertEqual(len(d["errors"]["stations"]), 3)

    def test_no_golden(self):
        shutil.rmtree(self.gold)
        os.makedirs(self.gold)
        rc, d = self.gate(tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("no station files" in f for f in d["failures"]))

    def test_nan_value(self):
        def edit(name, rows):
            if name == "st000dp075":
                rows[4][3] = float("nan")
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("non-finite" in f for f in d["failures"]))

    def test_short_line(self):
        def edit(name, rows):
            if name == "st000dp075":
                rows[3] = rows[3][:8]
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("columns" in f for f in d["failures"]))

    def test_time_repeats(self):
        def edit(name, rows):
            if name == "st000dp075":
                rows[-1][0] = rows[-2][0]
            return rows
        rc, d = self.gate(edit=edit, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("does not increase" in f for f in d["failures"]))

    def test_row_count(self):
        rc, d = self.gate(nrows=NROWS - 1, tol=1e-6)
        self.assertEqual(rc, 1)
        self.assertTrue(any("row count" in f for f in d["failures"]))

    def test_time_shift(self):
        def edit(name, rows):
            for r in rows:
                r[0] *= 1 + 1e-9
            return rows
        rc, d = self.gate(edit=edit, tol=1.0)
        self.assertEqual(rc, 1)
        self.assertTrue(any("times differ" in f for f in d["failures"]))


if __name__ == "__main__":
    unittest.main(verbosity=1)
