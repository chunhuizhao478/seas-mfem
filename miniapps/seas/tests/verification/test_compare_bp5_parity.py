#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# test_compare_bp5_parity.py — unit tests for the Phase-7 parity harness
# (tests/verification/compare_bp5_parity.py).  Locally runnable (no mesh/run):
# validates the comparison + parse logic, incl. the R-701 NaN-passes bug.
#
#   pytest tests/verification/test_compare_bp5_parity.py
#   (or: python3 tests/verification/test_compare_bp5_parity.py)

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_bp5_parity import compare_pair, parse_station, NCOL  # noqa: E402

FIN = [0.0, 0.0, 0.0, -9.0, -30.0, 30.0, 0.0, 1.5]   # one finite row


def test_identical_passes():
    a = [FIN[:], [1e6] + FIN[1:]]
    b = [r[:] for r in a]
    ok, _, _, na, nb = compare_pair(a, b, rtol=1e-6, atol=1e-9)
    assert ok is True and na == nb


def test_R701_nan_fails():
    # b has a NaN in log10V_strike (col 3); pre-fix this PASSED (NaN > atol is False).
    a = [FIN[:]]
    bad = FIN[:]
    bad[3] = float("nan")
    ok, _, _, _, _ = compare_pair(a, [bad], rtol=1e-6, atol=1e-9)
    assert ok is False


def test_R701_inf_fails():
    a = [FIN[:]]
    bad = FIN[:]
    bad[4] = float("-inf")
    ok, _, _, _, _ = compare_pair(a, [bad], rtol=1e-6, atol=1e-9)
    assert ok is False


def test_R704_log_column_abs_tol():
    # A small log10V diff (0.5) must fail under log_atol=1e-6 ...
    a = [FIN[:]]
    b = FIN[:]
    b[3] = FIN[3] + 0.5
    ok, _, _, _, _ = compare_pair(a, [b], rtol=1e-6, atol=1e-9, log_atol=1e-6)
    assert ok is False
    # ... but pass under a loose log_atol.
    ok2, _, _, _, _ = compare_pair(a, [b], rtol=1e-6, atol=1e-9, log_atol=1.0)
    assert ok2 is True


def test_R705_time_misalignment_fails():
    a = [[0.0] + FIN[1:]]
    b = [[5.0] + FIN[1:]]   # same data, different time -> not comparable
    ok, _, _, _, _ = compare_pair(a, b, rtol=1e-6, atol=1e-9)
    assert ok is False


def test_R706_corrupt_row_raises(tmp_path):
    p = tmp_path / "spatial_seas_fltst.txt"
    p.write_text(
        "time slip_strike slip_dip log10V_strike log10V_dip "
        "tau_strike tau_dip log10state\n"
        "0 0 0 -9 -30 30 0 1.5\n"
        "1e6 GARBAGE 0 -9 -30 30 0 1.5\n")
    try:
        parse_station(str(p))
    except ValueError:
        return
    assert False, "expected ValueError on a corrupt data row"


def test_row_count_mismatch_fails():
    a = [FIN[:], [1e6] + FIN[1:]]
    b = [FIN[:]]
    ok, _, _, na, nb = compare_pair(a, b, rtol=1e-6, atol=1e-9)
    assert ok is False and na != nb


if __name__ == "__main__":
    import tempfile
    n = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            if "tmp_path" in fn.__code__.co_varnames:
                import pathlib
                with tempfile.TemporaryDirectory() as d:
                    fn(pathlib.Path(d))
            else:
                fn()
            n += 1
            print("PASS", name)
    print("\n{} tests passed".format(n))
