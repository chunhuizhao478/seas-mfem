#!/usr/bin/env python3
"""Unit tests for project_friction_to_fault.py (pure functions).

    cd project_7.0_alternative/friction/code && pytest -q test_project_friction_to_fault.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from project_friction_to_fault import (
    read_profile_csv,
    eval_profile,
    transition_depth_km,
    DEFAULT_PARAM_A_CSV,
    DEFAULT_PARAM_AMB_CSV,
)


# ----------------------------------------------------------------------
# read_profile_csv
# ----------------------------------------------------------------------
def test_read_profile_csv_roundtrip(tmp_path: Path):
    p = tmp_path / "prof.csv"
    p.write_text("-0.5, 0\n0.0, 11\n0.5, 60\n")
    depth, val = read_profile_csv(p)
    assert depth.tolist() == [0.0, 11.0, 60.0]
    assert val.tolist() == [-0.5, 0.0, 0.5]


def test_read_profile_csv_skips_blank_and_comment(tmp_path: Path):
    p = tmp_path / "prof.csv"
    p.write_text("# header\n\n0.01, 0\n0.02, 10\n")
    depth, val = read_profile_csv(p)
    assert depth.tolist() == [0.0, 10.0]
    assert val.tolist() == [0.01, 0.02]


def test_read_profile_csv_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        read_profile_csv(tmp_path / "nope.csv")


def test_read_profile_csv_too_few_rows(tmp_path: Path):
    p = tmp_path / "prof.csv"
    p.write_text("0.01, 0\n")
    with pytest.raises(ValueError, match="need >= 2 data rows"):
        read_profile_csv(p)


def test_read_profile_csv_non_increasing_depth(tmp_path: Path):
    p = tmp_path / "prof.csv"
    p.write_text("0.01, 0\n0.02, 10\n0.03, 10\n")  # duplicate depth
    with pytest.raises(ValueError, match="strictly increasing"):
        read_profile_csv(p)


def test_read_profile_csv_non_numeric(tmp_path: Path):
    p = tmp_path / "prof.csv"
    p.write_text("0.01, 0\nfoo, 10\n")
    with pytest.raises(ValueError, match="non-numeric"):
        read_profile_csv(p)


# ----------------------------------------------------------------------
# eval_profile (flat-clamped piecewise-linear)
# ----------------------------------------------------------------------
def test_eval_profile_linear_interior():
    d = np.array([0.0, 10.0])
    v = np.array([0.0, 1.0])
    got = eval_profile(d, v, np.array([0.0, 5.0, 10.0]))
    assert np.allclose(got, [0.0, 0.5, 1.0])


def test_eval_profile_flat_clamp_outside():
    d = np.array([5.0, 10.0])
    v = np.array([-0.5, 0.5])
    # below x[0] clamps to v[0]; above x[-1] clamps to v[-1]
    got = eval_profile(d, v, np.array([0.0, 100.0]))
    assert got[0] == pytest.approx(-0.5)
    assert got[1] == pytest.approx(0.5)


# ----------------------------------------------------------------------
# transition_depth_km
# ----------------------------------------------------------------------
def test_transition_depth_simple_root():
    depth = np.array([0.0, 11.0, 60.0])
    amb = np.array([-0.5, 0.0, 0.5])
    assert transition_depth_km(depth, amb) == pytest.approx(11.0)


def test_transition_depth_interior_root():
    depth = np.array([0.0, 8.0, 60.0])
    amb = np.array([-0.2, -0.1, 0.4])  # crosses 0 between 8 and 60
    root = transition_depth_km(depth, amb)
    # linear root: 8 + 0.1/(0.5) * 52 = 8 + 10.4 = 18.4
    assert root == pytest.approx(18.4)


def test_transition_depth_no_crossing_returns_none():
    depth = np.array([0.0, 10.0, 20.0])
    amb = np.array([-0.3, -0.2, -0.1])  # all VW, never crosses
    assert transition_depth_km(depth, amb) is None


# ----------------------------------------------------------------------
# Acceptance: the committed 11 km CSVs really put the transition at 11 km
# ----------------------------------------------------------------------
def test_committed_csvs_transition_at_11km():
    amb_depth, amb_val = read_profile_csv(DEFAULT_PARAM_AMB_CSV)
    t = transition_depth_km(amb_depth, amb_val)
    assert t == pytest.approx(11.0, abs=1e-6)


def test_committed_a_profile_all_positive():
    _d, a_val = read_profile_csv(DEFAULT_PARAM_A_CSV)
    assert np.all(a_val > 0.0)


def test_b_is_positive_over_fault_depth_range():
    # b(depth) = a(depth) - (a-b)(depth) must stay > 0 across the meshed
    # fault depths (0..~16.5 km) for a physical rate-and-state law.
    a_d, a_v = read_profile_csv(DEFAULT_PARAM_A_CSV)
    amb_d, amb_v = read_profile_csv(DEFAULT_PARAM_AMB_CSV)
    depths = np.linspace(0.0, 16.5, 200)
    a = eval_profile(a_d, a_v, depths)
    amb = eval_profile(amb_d, amb_v, depths)
    b = a - amb
    assert np.all(b > 0.0)
