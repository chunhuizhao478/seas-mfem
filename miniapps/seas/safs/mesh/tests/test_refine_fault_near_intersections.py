"""Unit tests for refine_fault_near_intersections.py.

Covers REVIEW_intersection_refinement_investigation.md:
  R-001: polyline-vertex triangles MUST be markable for refinement.
  R-004: a parametric helper that asserts γ_min between two .msh
         files is not bit-identical (used by integration runs to
         catch "refinement did nothing" outcomes).

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_refine_fault_near_intersections.py -v
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the parent (mesh/) importable.
_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import refine_fault_near_intersections as rfn  # noqa: E402


# ---------------------------------------------------------------------------
# R-001: polyline-vertex triangles get marked
# ---------------------------------------------------------------------------
def _two_triangles_sharing_edge() -> tuple[
        np.ndarray, np.ndarray, np.ndarray, set[tuple[int, int, int]]]:
    """Two triangles sharing edge (a, b).
    Vertices a and b are flagged as polyline vertices (shared with a
    second fault not modelled here).  Both triangles have a, b as
    vertices, so their centroids include the shared edge midpoint.
    """
    snap_m = 0.001
    # a, b on shared edge; c, d are interior apexes.
    V = np.asarray([
        [0.0, 0.0, 0.0],   # a (polyline)
        [1.0, 0.0, 0.0],   # b (polyline)
        [0.5, 1.0, 0.0],   # c (interior)
        [0.5, -1.0, 0.0],  # d (interior)
    ], dtype=np.float64)
    T = np.asarray([
        [0, 1, 2],   # tri 0: (a, b, c)
        [1, 0, 3],   # tri 1: (b, a, d)
    ], dtype=np.int64)

    polyline_pts_xyz = V[:2].copy()  # a and b
    polyline_keys = {
        rfn._snap_key(V[0], snap_m),
        rfn._snap_key(V[1], snap_m),
    }
    return V, T, polyline_pts_xyz, polyline_keys


def test_R001_polyline_adjacent_triangles_are_marked_by_default():
    """The R-001 fix: polyline-bordering triangles MUST be marked when
    `include_polyline_triangles=True` (the default after the fix).
    Both triangles share polyline vertices a and b; both should mark.
    """
    V, T, pts, keys = _two_triangles_sharing_edge()
    snap_m = 0.001
    # Use a band radius large enough to cover both triangle centroids.
    band = 5.0

    marked = rfn.mark_triangles_in_band(
        V, T, pts, keys, snap_m, band,
        include_polyline_triangles=True)
    assert marked.tolist() == [True, True], (
        f"R-001 fix: both polyline-adjacent triangles must be marked; "
        f"got {marked.tolist()}")


def test_R001_legacy_behaviour_excludes_polyline_triangles():
    """When the user opts into the pre-R-001 behaviour, polyline-
    vertex triangles are excluded — preserved for regression
    reproducibility.
    """
    V, T, pts, keys = _two_triangles_sharing_edge()
    snap_m = 0.001
    band = 5.0

    marked = rfn.mark_triangles_in_band(
        V, T, pts, keys, snap_m, band,
        include_polyline_triangles=False)
    # Both triangles have polyline vertices → both excluded.
    assert marked.tolist() == [False, False], (
        f"legacy mode: triangles with any polyline vertex must be "
        f"excluded; got {marked.tolist()}")


def test_R001_default_mark_function_arg_is_inclusive():
    """The function default must be the FIXED behaviour (include
    polyline-vertex triangles) so callers that don't update get
    the bug fix automatically.
    """
    V, T, pts, keys = _two_triangles_sharing_edge()
    marked_default = rfn.mark_triangles_in_band(V, T, pts, keys, 0.001, 5.0)
    marked_explicit = rfn.mark_triangles_in_band(
        V, T, pts, keys, 0.001, 5.0,
        include_polyline_triangles=True)
    np.testing.assert_array_equal(marked_default, marked_explicit)


# ---------------------------------------------------------------------------
# R-006: legacy alias exists
# ---------------------------------------------------------------------------
def test_R006_legacy_function_name_alias():
    """`collect_polyline_endpoints` must remain importable as an
    alias of `collect_polyline_vertices` for backward compatibility.
    """
    assert hasattr(rfn, "collect_polyline_endpoints")
    assert hasattr(rfn, "collect_polyline_vertices")
    assert rfn.collect_polyline_endpoints is rfn.collect_polyline_vertices


# ---------------------------------------------------------------------------
# R-004: bit-identity guard helper
# ---------------------------------------------------------------------------
def assert_gamma_min_changed(refined_value: float, baseline_value: float,
                              *, require_improvement: bool = True) -> None:
    """Assert that refinement actually moved the worst-tet quality.

    A bit-identical γ_min means HXT produced the same near-degenerate
    tet on the same unchanged geometry — strong evidence that the
    refinement did not reach the polyline-adjacent vertex set.

    Use this in integration tests / CI runs that compare a refined
    mesh against an unrefined baseline.

    Args:
      refined_value:  γ_min from the refined run.
      baseline_value: γ_min from the unrefined baseline run.
      require_improvement: if True (default), also assert the refined
        γ_min strictly improved.  Set False to allow regression while
        still failing on bit-identity.

    Raises:
      AssertionError: if the two values are bit-identical.
    """
    if refined_value == baseline_value:
        raise AssertionError(
            f"refinement did not change worst-tet quality: both runs "
            f"report γ_min = {refined_value!r} bit-identically.  HXT "
            f"produced the same degenerate tet on unchanged input "
            f"geometry near the polyline.  The refinement reached a "
            f"region that does not drive the worst sliver.  See "
            f"REVIEW_intersection_refinement_investigation.md R-004.")
    if require_improvement and refined_value <= baseline_value:
        raise AssertionError(
            f"refinement made γ_min worse (or no improvement): "
            f"refined={refined_value!r}, baseline={baseline_value!r}.  "
            f"Either the refinement is mis-targeted or it introduced "
            f"new degenerate tets elsewhere.")


def _read_gamma_min(report_path: Path) -> float:
    """Parse `gamma_min:` from a validate_msh.py validation_report.txt."""
    txt = report_path.read_text()
    for line in txt.splitlines():
        s = line.strip()
        if s.startswith("gamma_min:"):
            return float(s.split(":", 1)[1].strip())
    raise ValueError(f"no `gamma_min:` line found in {report_path}")


def test_R004_helper_flags_bit_identical_gamma_min():
    """The helper used by integration runs must reject bit-identical
    γ_min — the empirical signature of "refinement did nothing"
    documented in REVIEW R-004.
    """
    bit_identical = 1.8015246456285393e-10
    with pytest.raises(AssertionError, match="did not change"):
        assert_gamma_min_changed(bit_identical, bit_identical)


def test_R004_helper_accepts_an_actual_improvement():
    """An improvement (refined > baseline) is the success signal."""
    assert_gamma_min_changed(refined_value=0.05, baseline_value=1.8e-10)


def test_R004_helper_flags_regression_by_default():
    """Refinement that strictly worsens γ_min should fail — that's
    a different bug class than "no change" but equally bad.
    """
    with pytest.raises(AssertionError, match="worse"):
        assert_gamma_min_changed(refined_value=1e-12,
                                  baseline_value=1.8e-10)


def test_R004_helper_can_skip_improvement_check():
    """Some debugging workflows want to detect bit-identity but
    tolerate regression — `require_improvement=False`.
    """
    # Regression (refined < baseline) — should pass with the flag.
    assert_gamma_min_changed(refined_value=1e-12,
                              baseline_value=1.8e-10,
                              require_improvement=False)
    # But bit-identity still fails.
    with pytest.raises(AssertionError, match="did not change"):
        assert_gamma_min_changed(refined_value=1.8e-10,
                                  baseline_value=1.8e-10,
                                  require_improvement=False)


def test_R004_read_gamma_min_from_validation_report(tmp_path):
    """Parser must extract `gamma_min:` from the validator's text
    output format.
    """
    rpt = tmp_path / "validation_report.txt"
    rpt.write_text("\n".join([
        "validate_msh.py report — somefile.msh",
        "================================================",
        "[FAIL] 10_tet_quality: 78 sliver tet(s) ...",
        "        gamma_min: 1.8015246456285393e-10",
        "        gamma_mean: 0.81",
    ]))
    assert _read_gamma_min(rpt) == 1.8015246456285393e-10


def test_R004_existing_5b_vs_5d_runs_were_bit_identical():
    """Regression sanity: confirm the historical empirical observation
    behind R-004.  The 5b (no surface refinement) and 5d (with
    refinement applied) runs reported the IDENTICAL γ_min.  This test
    only runs if both validation reports exist; otherwise it skips.
    """
    repo_root = _HERE.parent.parent.parent.parent.parent  # safs/mesh/tests/ → repo
    candidates_baseline = [
        repo_root / "miniapps/seas/safs/mesh/output/"
                    "newset_5b_garnetfirst_2000/output/validation_report.txt",
    ]
    candidates_refined = [
        repo_root / "miniapps/seas/safs/mesh/output/"
                    "newset_5d_surface_refined/output/validation_report.txt",
    ]
    base = next((p for p in candidates_baseline if p.exists()), None)
    refn = next((p for p in candidates_refined if p.exists()), None)
    if base is None or refn is None:
        pytest.skip(f"missing validator outputs (base={base}, refn={refn})")
    g_base = _read_gamma_min(base)
    g_refn = _read_gamma_min(refn)
    # The R-004 finding: these were observed to be bit-identical.
    # The helper should flag that as "refinement did nothing".
    if g_base == g_refn:
        with pytest.raises(AssertionError, match="did not change"):
            assert_gamma_min_changed(g_refn, g_base)
    else:
        # If a future re-run produced different values, R-004 is no
        # longer reproducing — flag for human review.
        pytest.skip(
            f"5b vs 5d γ_min differ (base={g_base!r}, refn={g_refn!r}); "
            f"R-004 historical scenario no longer reproduces, but the "
            f"bit-identity guard helper itself is exercised by the "
            f"`flags_bit_identical_gamma_min` test above.")
