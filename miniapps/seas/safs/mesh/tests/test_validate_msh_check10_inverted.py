"""Unit tests for the inverted-tet HARD-fail in
validate_msh.check_10_tet_quality.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_validate_msh_check10_inverted.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import validate_msh  # noqa: E402


def _well_shaped_tet_points():
    """A regular-ish tet with positive signed volume."""
    return np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=np.float64)


def test_check10_passes_on_well_shaped_tet():
    points = _well_shaped_tet_points()
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    r = validate_msh.check_10_tet_quality(
        points, tets, tris=None, ttags=None)
    assert r.passed, f"expected PASS; got {r}"
    assert r.metrics["n_inverted_tets"] == 0
    assert r.metrics["worst_inverted_tet_idx"] == -1


def test_check10_fails_on_single_inverted_tet():
    """Swap two vertices to flip orientation → signed_vol < 0."""
    points = _well_shaped_tet_points()
    # (0,2,1,3) is the same regular tet with vertices 1<->2 swapped,
    # which inverts the signed volume.
    tets = np.array([[0, 2, 1, 3]], dtype=np.int64)
    r = validate_msh.check_10_tet_quality(
        points, tets, tris=None, ttags=None)
    assert not r.passed
    assert r.metrics["n_inverted_tets"] == 1
    assert r.metrics["worst_inverted_tet_idx"] == 0
    assert r.metrics["worst_inverted_tet_signed_vol_m3"] < 0.0
    assert "inverted tet" in r.message
    assert "DG solver" in r.message  # informative error


def test_check10_fails_on_zero_volume_degenerate_tet():
    """All four vertices collinear → signed_vol == 0 ⇒ FAIL."""
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    r = validate_msh.check_10_tet_quality(
        points, tets, tris=None, ttags=None)
    assert not r.passed
    assert r.metrics["n_inverted_tets"] == 1


def test_check10_inverted_tet_count_aggregates():
    """Mix of N positive + M inverted: count must be M and the worst
    must be the most-negative signed_vol."""
    pts = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 2.0],   # gives a much-more-negative signed_vol
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],   # positive
        [0, 2, 1, 3],   # negative (small)
        [0, 2, 1, 4],   # negative (larger magnitude due to z=-2)
    ], dtype=np.int64)
    r = validate_msh.check_10_tet_quality(
        pts, tets, tris=None, ttags=None)
    assert not r.passed
    assert r.metrics["n_inverted_tets"] == 2
    # The worst (most-negative) is index 2 (signed_vol ~ -1/3 vs -1/6).
    assert r.metrics["worst_inverted_tet_idx"] == 2
    assert (r.metrics["worst_inverted_tet_signed_vol_m3"]
            < -1.0 / 6.0 - 1e-9)


def test_check10_inverted_short_circuits_before_sliver_check():
    """Even if global γ stats look fine on |vol|, an inverted tet must
    HARD-fail BEFORE the sliver-fraction gate runs (so the message
    points at the inversion, not at the sliver count)."""
    pts = _well_shaped_tet_points()
    # Single inverted tet; both metrics aggregate cleanly otherwise.
    tets = np.array([[0, 2, 1, 3]], dtype=np.int64)
    r = validate_msh.check_10_tet_quality(
        pts, tets, tris=None, ttags=None)
    assert not r.passed
    # Message must reference the inversion, NOT the sliver gate.
    assert "inverted" in r.message.lower()
    assert "sliver" not in r.message.lower()
