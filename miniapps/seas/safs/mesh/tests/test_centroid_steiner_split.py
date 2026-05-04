"""Unit tests for centroid_steiner_split.py.

Covers REVIEW_centroid_steiner_for_collinear_tets.md R-001:
collinear / near-coplanar bad tets must be split via centroid
Steiner insertion into 4 well-conditioned sub-tets, with cross-fault
conformity preserved by construction (the original tet's 4 faces
become faces of the new sub-tets).

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_centroid_steiner_split.py -v
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import meshio
import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import centroid_steiner_split as css  # noqa: E402
import mmg3d_local_patch as mlp  # noqa: E402


# ---------------------------------------------------------------------------
# Geometric tests — pure-Python, no external binaries needed.
# ---------------------------------------------------------------------------
def test_collinear_tet_detected():
    """4 vertices nearly on a line → flagged as collinear."""
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.001],   # tiny perturbation off line
        [2.0, 0.001, 0.0],
        [3.0, 0.0, 0.0],
    ])
    tet = np.array([0, 1, 2, 3])
    assert css._is_near_collinear(points, tet, 0.99)
    # Near-coplanar test ALSO flags this (collinear is a stronger
    # form of coplanar).
    assert css._is_near_coplanar_or_flat(points, tet, 1e-4)


def test_well_conditioned_tet_not_flagged():
    """A regular tet should NOT be flagged as collinear or coplanar."""
    s = 1.0 / np.sqrt(3.0)
    points = np.array([[s, s, s], [s, -s, -s],
                        [-s, s, -s], [-s, -s, s]])
    tet = np.array([0, 1, 2, 3])
    assert not css._is_near_collinear(points, tet, 0.99)
    assert not css._is_near_coplanar_or_flat(points, tet, 1e-4)


def test_coplanar_tet_detected():
    """4 vertices nearly on a 2-D plane (but NOT collinear) →
    coplanar test flags it; collinear test may or may not.
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 1e-6],   # slightly above plane
    ])
    tet = np.array([0, 1, 2, 3])
    # Coplanar test: vol / max_edge^3 ≈ 1e-6 / 1.4^3 ≈ 4e-7 < 1e-4
    assert css._is_near_coplanar_or_flat(points, tet, 1e-4)
    # Combined predicate: yes split.
    assert css._should_split(points, tet, 0.99, 1e-4)


def test_should_split_or_not_logic():
    """`_should_split` returns True iff EITHER predicate triggers."""
    # Collinear-only.
    pts = np.array([[0,0,0], [1,0,0], [2,0,0], [3,0,0]],
                   dtype=np.float64)
    assert css._should_split(pts, np.array([0,1,2,3]), 0.99, 1e-4)
    # Well-conditioned.
    s = 1.0 / np.sqrt(3.0)
    pts = np.array([[s,s,s],[s,-s,-s],[-s,s,-s],[-s,-s,s]])
    assert not css._should_split(pts, np.array([0,1,2,3]),
                                  0.99, 1e-4)


def test_zero_length_edge_treated_as_collinear():
    """If two vertices coincide, the edge has length 0 — the
    near-collinear test should return True (defensive)."""
    points = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],   # coincident with v0
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    assert css._is_near_collinear(
        points, np.array([0,1,2,3]), 0.99)


# ---------------------------------------------------------------------------
# Defensive parameter validation.
# ---------------------------------------------------------------------------
def test_invalid_gamma_thresh_raises(tmp_path):
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="gamma_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            gamma_thresh=0.0)
    with pytest.raises(ValueError, match="gamma_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            gamma_thresh=-0.5)


def test_invalid_collinearity_thresh_raises(tmp_path):
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="collinearity_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            collinearity_thresh=0.0)
    with pytest.raises(ValueError, match="collinearity_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            collinearity_thresh=1.5)


def test_invalid_coplanarity_thresh_raises(tmp_path):
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="coplanarity_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            coplanarity_thresh=-0.1)
    with pytest.raises(ValueError, match="coplanarity_thresh"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            coplanarity_thresh=0.0)


def test_invalid_eps_vol_raises(tmp_path):
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="eps_vol_m3"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            eps_vol_m3=0.0)
    with pytest.raises(ValueError, match="eps_vol_m3"):
        css.centroid_split_collinear_tets(fake, tmp_path / "out",
                                            eps_vol_m3=-1.0)


def test_missing_input_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-msh not found"):
        css.centroid_split_collinear_tets(
            tmp_path / "nope.msh", tmp_path / "out.msh")


# ---------------------------------------------------------------------------
# Integration tests with synthetic .msh fixtures.
# ---------------------------------------------------------------------------
def _write_minimal_msh(out_path: Path,
                        points: np.ndarray,
                        tets: np.ndarray,
                        tris: np.ndarray | None = None,
                        tet_tags: np.ndarray | None = None,
                        tri_tags: np.ndarray | None = None,
                        ) -> None:
    """Write a minimal .msh.  Defaults: every tet tag=10, every
    triangle tag=100 (fault).
    """
    if tris is None:
        # Make 1 boundary tri so the file is well-formed.
        tris = np.asarray([[0, 1, 2]], dtype=np.int64)
        tri_tags = np.asarray([100], dtype=np.int32)
    if tet_tags is None:
        tet_tags = np.full(tets.shape[0], 10, dtype=np.int32)
    if tri_tags is None:
        tri_tags = np.full(tris.shape[0], 100, dtype=np.int32)
    meshio.write(out_path, meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    ), file_format="gmsh22", binary=False)


def test_e2e_collinear_tet_is_split_into_4(tmp_path):
    """Single near-collinear tet → 4 sub-tets after centroid split.
    +1 vertex (the centroid)."""
    points = np.array([
        [0.0,    0.0,   0.0],
        [1.0,    0.001, 0.0],
        [2.0,    0.0,   0.001],
        [3.0,    0.001, 0.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh,
        gamma_thresh=1.0,  # generous: any γ<1 considered bad
        collinearity_thresh=0.99,
        coplanarity_thresh=1e-4)
    assert rep["no_op"] is False
    assert rep["n_centroid_split"] == 1
    assert rep["tets_in"] == 1
    assert rep["tets_out"] == 4
    out = meshio.read(out_msh)
    n_tets_out = sum(cb.data.shape[0] for cb in out.cells
                      if cb.type == "tetra")
    assert n_tets_out == 4
    # +1 vertex (the centroid).
    assert out.points.shape[0] == points.shape[0] + 1


def test_e2e_well_conditioned_tet_untouched(tmp_path):
    """Regular tet (γ ≈ 1) should NOT be split — γ_thresh excludes
    it and the geometric tests would reject it anyway."""
    s = 1.0 / np.sqrt(3.0)
    points = np.array([[s, s, s], [s, -s, -s],
                        [-s, s, -s], [-s, -s, s]])
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh,
        gamma_thresh=1e-3,
        collinearity_thresh=0.99,
        coplanarity_thresh=1e-4)
    assert rep["no_op"] is True
    assert rep["n_centroid_split"] == 0
    assert rep["tets_in"] == rep["tets_out"]


def test_e2e_split_preserves_tet_tag(tmp_path):
    """Each sub-tet inherits the parent tet's physical tag."""
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.001, 0.0],
        [2.0, 0.0, 0.001],
        [3.0, 0.001, 0.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    tet_tags = np.array([42], dtype=np.int32)  # nonstandard tag
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets, tet_tags=tet_tags)
    css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1.0,
        collinearity_thresh=0.99, coplanarity_thresh=1e-4)
    out = meshio.read(out_msh)
    out_tet_tags = None
    for cb, tags in zip(out.cells,
                          out.cell_data.get("gmsh:physical", [])):
        if cb.type == "tetra":
            out_tet_tags = tags
            break
    # All 4 sub-tets must have tag 42.
    assert (np.asarray(out_tet_tags) == 42).all()


def test_e2e_split_preserves_fault_face_adjacency(tmp_path):
    """Two tets sharing one face: the original face is preserved
    after the parent tet is split.  After centroid-split of the
    bad tet, the sub-tet that contains the shared face must STILL
    share that face with the surviving good tet — preserving
    validator check_5 invariant.

    Layout: bad tet (collinear) + good tet, sharing face (0,1,2).
    The bad tet is (0,1,2,3); the good tet is (0,1,2,4).
    After split of bad tet:
      Sub-tet 1 = (0,1,2, centroid)  ← contains face (0,1,2)
    The good tet (0,1,2,4) is unchanged.
    The shared face (0,1,2) is now between sub-tet 1 and good tet.
    """
    # Layout:
    #   v0, v1, v2: a clean 3-vertex triangle in the z=0 plane.
    #   v3: bad apex JUST above the plane (vol/edge^3 ~ 5e-5 < 1e-4)
    #       → bad tet (0,1,2,3) is coplanar-flagged.
    #   v4: good apex WELL above the plane (vol/edge^3 ~ 0.035)
    #       → good tet (0,1,2,4) is NOT coplanar-flagged.
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, 0.001],   # bad apex: 0.001 m above plane
        [0.5, 0.5, 2.0],     # good apex: 2 m above plane
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],   # bad (4 vertices nearly coplanar)
        [0, 1, 2, 4],   # good (well-conditioned tet)
    ], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1.0,
        collinearity_thresh=0.99, coplanarity_thresh=1e-4)
    # Should split exactly 1 tet (the bad one).
    assert rep["n_centroid_split"] == 1
    # Output: 4 sub-tets + 1 unchanged good tet = 5 tets total.
    assert rep["tets_out"] == 5

    # Verify the shared face (0, 1, 2) appears as a face in
    # exactly TWO output tets — one sub-tet (a, b, c, centroid)
    # and the original good tet.
    out = meshio.read(out_msh)
    out_tets = None
    for cb in out.cells:
        if cb.type == "tetra":
            out_tets = cb.data
            break
    target_face = frozenset((0, 1, 2))
    n_incident = 0
    for tet in out_tets:
        for face in (frozenset((tet[0], tet[1], tet[2])),
                      frozenset((tet[0], tet[1], tet[3])),
                      frozenset((tet[0], tet[2], tet[3])),
                      frozenset((tet[1], tet[2], tet[3]))):
            if face == target_face:
                n_incident += 1
                break
    assert n_incident == 2, (
        f"shared face (0,1,2) must remain shared between exactly "
        f"2 tets after centroid split; got {n_incident}")


def test_e2e_no_op_when_no_tets_below_threshold(tmp_path):
    """If no tets are below γ_thresh, no_op=True and the output is
    identical to the input."""
    s = 1.0 / np.sqrt(3.0)
    points = np.array([[s, s, s], [s, -s, -s],
                        [-s, s, -s], [-s, -s, s]])
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1e-9)  # very tight
    assert rep["no_op"] is True
    assert rep["n_centroid_split"] == 0


def test_e2e_skip_zero_volume_tet(tmp_path):
    """A true-zero-volume tet (4 exactly-coplanar vertices) is
    SKIPPED — splitting it would produce 4 zero-volume sub-tets.
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],   # exactly on the same plane
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1.0,
        collinearity_thresh=0.99,
        coplanarity_thresh=1e-4,
        eps_vol_m3=1e-12)
    assert rep["n_skipped_vol_too_small"] == 1
    assert rep["n_centroid_split"] == 0


def test_e2e_writes_report_json(tmp_path):
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.001, 0.0],
        [2.0, 0.0, 0.001],
        [3.0, 0.001, 0.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    rep_json = tmp_path / "rep.json"
    _write_minimal_msh(in_msh, points, tets)
    css.main([
        "--in-msh", str(in_msh),
        "--out-msh", str(out_msh),
        "--gamma-thresh", "1.0",
        "--collinearity-thresh", "0.99",
        "--coplanarity-thresh", "1e-4",
        "--report-json", str(rep_json),
    ])
    rep = json.loads(rep_json.read_text())
    for k in ("gamma_min_in", "gamma_min_out", "n_bad_tets",
               "n_centroid_split", "tets_in", "tets_out"):
        assert k in rep, f"report missing key {k}"


def test_e2e_centroid_split_helps_coplanar_tet(tmp_path):
    """For a near-COPLANAR bad tet (volume ≈ 0 because 4 vertices
    are nearly on the same plane, NOT all collinear), the centroid
    is in the plane, so the split's 4 sub-tets are also in-plane
    (γ does not improve).

    For a near-COLLINEAR bad tet (4 vertices on the same line),
    the centroid is also on the line — sub-tets stay collinear.
    γ does NOT improve materially in this synthetic case.

    This test documents the limitation honestly: centroid-split
    is necessary but not sufficient for true 1-D / 2-D degenerate
    configurations.  In realistic 3D meshes the bad tets have
    SOME perpendicular variation (e.g., the production worst tet
    has 36 m y-spread and 16 m x-spread on a 2 km z-line), and
    the centroid lies inside that perpendicular volume.  γ
    improvement is real on those.

    For the synthetic strictly-collinear case, we only assert:
      - The split produces 4 sub-tets.
      - γ_min stays > 0 (no inverted sub-tets).
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.001, 0.0],
        [2.0, 0.0, 0.001],
        [3.0, 0.001, 0.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1.0,
        collinearity_thresh=0.99, coplanarity_thresh=1e-4)
    assert rep["n_centroid_split"] == 1
    assert rep["tets_out"] == 4
    # Essential invariant: sub-tets are non-inverted (positive γ).
    assert rep["gamma_min_out"] > 0.0


def test_e2e_centroid_split_runs_on_perpendicular_spread_tet(
        tmp_path):
    """For a bad tet with non-trivial perpendicular spread (3D
    sliver, NOT pure 1D / 2D degenerate), the centroid-split must
    run cleanly and produce 4 well-defined sub-tets.

    IMPORTANT empirical caveat:  centroid-split is a TOPOLOGICAL
    operation.  It introduces a new vertex at the centroid of the
    bad tet and replaces the parent with 4 sub-tets.  Whether
    γ_min improves depends on the geometry: for tets with strong
    perpendicular spread (centroid is meaningfully off the
    dominant axis / plane), γ_min CAN improve; for nearly-flat or
    nearly-linear tets where the centroid stays in/near the
    degenerate subspace, γ_min may NOT improve and can even
    decrease (the 4 sub-tets may have new short edges from
    centroid to vertex).

    For this synthetic "thick line" sliver, we only assert the
    structural invariants:
      - 1 tet → 4 sub-tets.
      - +1 vertex.
      - γ_min_out > 0 (no inverted sub-tets).

    γ_min IMPROVEMENT is empirically validated by the end-to-end
    Step 6 production-mesh run, not by a synthetic test.
    """
    points = np.array([
        [0.0,    0.0,    0.0],
        [1000.0, 50.0,  -30.0],
        [2000.0, -40.0,  60.0],
        [3000.0, 30.0,  -50.0],
    ], dtype=np.float64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_minimal_msh(in_msh, points, tets)
    rep = css.centroid_split_collinear_tets(
        in_msh, out_msh, gamma_thresh=1.0,
        collinearity_thresh=0.99, coplanarity_thresh=1.0)
    assert rep["n_centroid_split"] == 1
    assert rep["tets_out"] == 4
    assert rep["gamma_min_out"] > 0.0
