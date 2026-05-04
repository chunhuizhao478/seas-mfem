"""Unit tests for validate_msh.py check_12_surface_closure and
check_13_fault_orientation (REVIEW.md R-001 / R-002).

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_validate_msh_check12_check13.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import validate_msh  # noqa: E402


# ---------------------------------------------------------------------------
# R-001: check_12_surface_closure
# ---------------------------------------------------------------------------
def test_R001_check12_detects_surface_hole():
    """A tet whose every face has no tagged tri must FAIL check_12."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0.5, 0.5, -1.0]], dtype=np.float64)
    tris = np.array([], dtype=np.int64).reshape(0, 3)
    ttags = np.array([], dtype=np.int32)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    r = validate_msh.check_12_surface_closure(points, tris, ttags, tets)
    assert not r.passed
    assert r.metrics["n_unlabeled_holes"] == 4


def test_R001_check12_passes_on_closed_mesh():
    """All 4 faces tagged → check_12 PASS, 0 unlabeled holes."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0.5, 0.5, -1.0]], dtype=np.float64)
    tris = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
                    dtype=np.int64)
    ttags = np.array([5, 4, 1, 2], dtype=np.int32)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    r = validate_msh.check_12_surface_closure(points, tris, ttags, tets)
    assert r.passed
    assert r.metrics["n_unlabeled_holes"] == 0
    assert r.metrics["n_bdry_faces"] == 4


def test_R001_check12_partial_holes():
    """3 faces tagged + 1 missing → 1 hole reported, FAIL."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0.5, 0.5, -1.0]], dtype=np.float64)
    tris = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3]], dtype=np.int64)
    ttags = np.array([5, 4, 1], dtype=np.int32)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    r = validate_msh.check_12_surface_closure(points, tris, ttags, tets)
    assert not r.passed
    assert r.metrics["n_unlabeled_holes"] == 1


def test_R001_check12_no_tets_fails():
    """Empty mesh fails check_12."""
    points = np.zeros((0, 3), dtype=np.float64)
    tris = np.zeros((0, 3), dtype=np.int64)
    ttags = np.zeros(0, dtype=np.int32)
    tets = np.zeros((0, 4), dtype=np.int64)
    r = validate_msh.check_12_surface_closure(points, tris, ttags, tets)
    assert not r.passed


# ---------------------------------------------------------------------------
# R-002: check_13_fault_orientation
# ---------------------------------------------------------------------------
def test_R002_check13_detects_flipped_winding():
    """Two adjacent fault tris that traverse a shared edge in the SAME
    direction must FAIL check_13 with n_winding_flips_needed == 1.

    Tri-A=(0,1,2) traverses edge (1,2) in direction 1→2.
    Tri-B=(1,2,3) ALSO traverses edge (1,2) in direction 1→2.
    Same direction across the shared edge ⇒ inconsistent winding
    (their outward normals lie on opposite sides of the surface).
    """
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
    ttags = np.array([100, 100], dtype=np.int32)
    r = validate_msh.check_13_fault_orientation(points, tris, ttags)
    assert not r.passed
    assert r.metrics["n_winding_flips_needed"] == 1


def test_R002_check13_passes_on_consistent_winding():
    """Two adjacent tris that traverse a shared edge in OPPOSITE
    directions PASS check_13.

    Tri-A=(0,1,2) traverses edge (2,0) in direction 2→0.
    Tri-B=(0,2,3) traverses edge (0,2) in direction 0→2.
    Opposite directions across the shared edge ⇒ consistent winding.
    """
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    ttags = np.array([100, 100], dtype=np.int32)
    r = validate_msh.check_13_fault_orientation(points, tris, ttags)
    assert r.passed
    assert r.metrics["n_winding_flips_needed"] == 0
    assert r.metrics["n_components"] == 1


def test_R002_check13_reports_nonmanifold_edges_as_metric():
    """A fault edge shared by 3 fault tris is non-manifold (>2).

    Per REVIEW.md R-002: non-manifold fault edges are unavoidable for
    branching multi-fault SAFS surfaces (e.g., MJVS + SAF along a shared
    polyline).  check_13 must REPORT the count as a metric but PASS as
    long as winding flips == 0 — non-manifold is "soft warn", not fail.
    """
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0, 0, 1], [0, -1, 0]], dtype=np.float64)
    # All three tris share edge (0,1).  Tri 0 traverses 0->1; tris 1
    # and 2 traverse 1->0 (opposite to tri 0).  Per the BFS algorithm
    # in check_13, both neighbors are reached from tri 0 with
    # OPPOSITE direction on the shared edge, so n_flips_needed == 0;
    # the (1,2) pair is not visited because BFS marks them visited
    # before they can revisit each other.  Non-manifold-edge count is
    # still 1 (edge (0,1) shared by 3 tris).
    tris = np.array([[0, 1, 2], [1, 0, 3], [1, 0, 4]], dtype=np.int64)
    ttags = np.array([100, 100, 100], dtype=np.int32)
    r = validate_msh.check_13_fault_orientation(points, tris, ttags)
    # Non-manifold edge is reported but does not fail the check.
    assert r.metrics["n_nonmanifold_edges"] >= 1
    assert r.passed, (
        "non-manifold edges with 0 winding flips must be soft-warn, "
        "not hard fail (REVIEW.md R-002 spec)")


def test_R002_check13_no_fault_tris_fails():
    """No tag-100 triangles → FAIL."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    tris = np.array([[0, 1, 2]], dtype=np.int64)
    ttags = np.array([5], dtype=np.int32)
    r = validate_msh.check_13_fault_orientation(points, tris, ttags)
    assert not r.passed


def test_R002_check13_passes_on_box_only_mesh():
    """No fault tris and no other tris also fails (no fault to check)."""
    # A consistently-wound 2-tri patch (NOT tagged 100) is irrelevant.
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    ttags = np.array([5, 5], dtype=np.int32)
    r = validate_msh.check_13_fault_orientation(points, tris, ttags)
    assert not r.passed  # no fault tris → fail
