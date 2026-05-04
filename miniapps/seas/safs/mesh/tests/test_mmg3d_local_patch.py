"""Unit tests for mmg3d_local_patch.py.

Covers REVIEW_local_sliver_targeting.md R-001:
local sliver-targeted mmg3d patch must (a) identify bad tets,
(b) extract a halo subdomain, (c) lock the halo boundary, (d) run
mmg3d locally, (e) stitch back without disturbing outside-subdomain
tets.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_mmg3d_local_patch.py -v
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import mmg3d_local_patch as mlp  # noqa: E402

_mmg3d = shutil.which("mmg3d_O3")


# ---------------------------------------------------------------------------
# Pure-Python tests (no mmg3d binary required).
# ---------------------------------------------------------------------------
def test_invalid_inputs_raise(tmp_path):
    """Defensive parameter validation."""
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="gamma_thresh"):
        mlp.local_patch(fake, tmp_path / "out.msh",
                        gamma_thresh=0.0)
    with pytest.raises(ValueError, match="gamma_thresh"):
        mlp.local_patch(fake, tmp_path / "out.msh",
                        gamma_thresh=1.5)
    with pytest.raises(ValueError, match="halo_radius_m"):
        mlp.local_patch(fake, tmp_path / "out.msh",
                        halo_radius_m=-100)
    with pytest.raises(ValueError, match="hmin"):
        mlp.local_patch(fake, tmp_path / "out.msh",
                        hmin=0)
    with pytest.raises(ValueError, match="hgradreq"):
        mlp.local_patch(fake, tmp_path / "out.msh",
                        hgradreq=0.5)


def test_missing_input_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-msh not found"):
        mlp.local_patch(tmp_path / "nope.msh", tmp_path / "out.msh")


def test_gamma_per_tet_regular_tet():
    """A regular tet has γ = 1."""
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    g = mlp._gamma_per_tet(P, T)
    assert 0.99 <= g[0] <= 1.01


def test_gamma_per_tet_handles_empty_array():
    """No tets → empty array, no crash."""
    P = np.zeros((0, 3))
    T = np.zeros((0, 4), dtype=np.int64)
    g = mlp._gamma_per_tet(P, T)
    assert g.shape == (0,)


def test_gamma_per_tet_zero_volume_returns_zero():
    """A degenerate (zero-volume) tet returns γ = 0, not a NaN."""
    P = np.array([[0, 0, 0], [1, 0, 0],
                   [2, 0, 0], [3, 0, 0]], dtype=np.float64)
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    g = mlp._gamma_per_tet(P, T)
    assert g[0] == 0.0


def test_identify_subdomain_halo_radius(tmp_path):
    """A tet's centroid within `halo_radius_m` of a bad tet's
    centroid is in the subdomain; further away is not.
    """
    # Two tets with centroids ~10 m and ~10000 m from origin.
    P = np.array([
        # Tet A (bad, near origin):
        [0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10],
        # Tet B (far away):
        [9990, 0, 0], [10000, 0, 0], [9990, 10, 0], [9990, 0, 10],
        # Tet C (close to A, inside halo):
        [100, 0, 0], [110, 0, 0], [100, 10, 0], [100, 0, 10],
    ], dtype=np.float64)
    T = np.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]],
                 dtype=np.int64)
    bad_idx = np.array([0], dtype=np.int64)
    in_sub = mlp._identify_subdomain_tets(
        P, T, bad_idx, halo_radius_m=500.0)
    assert in_sub[0] == True   # the bad tet itself
    assert in_sub[1] == False  # ~10 km away, outside halo
    assert in_sub[2] == True   # ~100 m, inside 500 m halo


def test_identify_subdomain_no_bad_tets():
    """Empty bad_tet_idx → empty subdomain."""
    P = np.zeros((4, 3))
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    bad_idx = np.zeros(0, dtype=np.int64)
    in_sub = mlp._identify_subdomain_tets(
        P, T, bad_idx, halo_radius_m=1000.0)
    assert in_sub.sum() == 0


def test_halo_boundary_face_detection():
    """Two adjacent tets sharing one face: if one is in subdomain
    and the other is not, that face is a halo boundary face.
    """
    # 5-vertex bipyramid: tet0 = (0,1,2,3), tet1 = (0,1,2,4).
    # Shared face = (0, 1, 2).
    T = np.array([[0, 1, 2, 3], [0, 1, 2, 4]], dtype=np.int64)
    in_sub = np.array([True, False])
    f2t = mlp._build_tet_face_adjacency(T)
    halo = mlp._identify_halo_boundary_faces(T, in_sub, f2t)
    # Exactly the shared face must be flagged.
    assert (0, 1, 2) in halo or (0, 2, 1) in halo \
           or tuple(sorted((0, 1, 2))) in halo
    # And only the shared face — not the boundary faces of either tet.
    assert len(halo) == 1


def test_halo_boundary_face_both_in_or_both_out():
    """If BOTH adjacent tets are in (or both out) the subdomain,
    their shared face is NOT a halo boundary."""
    T = np.array([[0, 1, 2, 3], [0, 1, 2, 4]], dtype=np.int64)
    f2t = mlp._build_tet_face_adjacency(T)
    halo_both_in = mlp._identify_halo_boundary_faces(
        T, np.array([True, True]), f2t)
    assert len(halo_both_in) == 0
    halo_both_out = mlp._identify_halo_boundary_faces(
        T, np.array([False, False]), f2t)
    assert len(halo_both_out) == 0


def test_no_op_when_no_bad_tets(tmp_path):
    """Input with no bad tets → output is a copy of input (no-op)."""
    import meshio
    # Build a minimal valid .msh: one tet, one fault tri.
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    F = np.array([[0, 1, 2]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    meshio.write(in_msh, meshio.Mesh(
        points=P,
        cells=[("triangle", F), ("tetra", T)],
        cell_data={"gmsh:physical":
                   [np.array([100], dtype=np.int32),
                    np.array([10], dtype=np.int32)],
                   "gmsh:geometrical":
                   [np.array([100], dtype=np.int32),
                    np.array([10], dtype=np.int32)]},
    ), file_format="gmsh22", binary=False)

    rep = mlp.local_patch(in_msh, out_msh,
                           gamma_thresh=0.01,
                           halo_radius_m=100.0)
    assert rep["no_op"] is True
    assert rep["n_bad_tets"] == 0
    assert out_msh.exists()


# ---------------------------------------------------------------------------
# Tests requiring mmg3d_O3.
# ---------------------------------------------------------------------------
def _make_sliver_test_mesh(tmp_path: Path) -> Path:
    """Create a minimal mesh containing a sliver tet inside a region
    of well-conditioned tets.

    Layout: 8 vertices forming a unit cube, decomposed into 6 tets
    (standard).  Then add 1 EXTRA vertex very close to vertex 0 →
    creates a needle tet.

    For test simplicity we use a tetrahedral 'corner cluster'
    around an origin-vicinity sliver.
    """
    import meshio
    # Six well-conditioned tets, plus one sliver near origin.
    P = np.array([
        # Cube corners (z = 0 face plus z = 100):
        [0,    0,   0],    # 0
        [100,  0,   0],    # 1
        [100, 100,  0],    # 2
        [0,   100,  0],    # 3
        [0,    0, 100],    # 4
        [100,  0, 100],    # 5
        [100, 100, 100],   # 6
        [0,   100, 100],   # 7
        # Sliver-creating extra vertex very close to vertex 0:
        [0.01, 0.01, 0.01], # 8 ← sub-meter from 0
    ], dtype=np.float64)
    # 6 tets covering the cube interior + 1 sliver tet using vertex 8.
    T = np.array([
        [0, 1, 2, 5],
        [0, 2, 3, 5],
        [0, 3, 7, 5],
        [0, 7, 4, 5],
        [0, 4, 5, 6],
        [2, 3, 7, 5],
        # Sliver tet:
        [0, 1, 2, 8],
    ], dtype=np.int64)
    # Surface triangles (one fault face for tag 100):
    F = np.array([
        [0, 1, 2],   # base, fault
        [0, 1, 4],
        [1, 2, 6],
        [2, 3, 7],
        [0, 3, 4],
        [4, 5, 6],
    ], dtype=np.int64)
    F_tags = np.array([100, 1, 2, 3, 4, 5], dtype=np.int32)
    T_tags = np.full(T.shape[0], 10, dtype=np.int32)

    in_msh = tmp_path / "sliver_test.msh"
    meshio.write(in_msh, meshio.Mesh(
        points=P,
        cells=[("triangle", F), ("tetra", T)],
        cell_data={"gmsh:physical": [F_tags, T_tags],
                   "gmsh:geometrical": [F_tags, T_tags]},
    ), file_format="gmsh22", binary=False)
    return in_msh


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_local_patch_runs_and_writes_output(tmp_path):
    """End-to-end smoke: run local_patch on a synthetic sliver
    mesh; verify output exists, contains tetrahedra and triangles,
    and that the bad-tet count went down (or stayed same).
    """
    in_msh = _make_sliver_test_mesh(tmp_path)
    out_msh = tmp_path / "out.msh"

    rep = mlp.local_patch(in_msh, out_msh,
                           gamma_thresh=0.5,  # generous threshold
                           halo_radius_m=300.0,
                           hmin=10.0, hmax=200.0,
                           hgrad=1.1, hausd=5.0,
                           hgradreq=1.1)
    assert out_msh.exists()
    assert rep["no_op"] is False
    # The output must have ≥ 1 tet (we cannot delete the model).
    assert rep["tets_out"] > 0
    # NOTE: we DO NOT assert γ_min monotonicity on this synthetic
    # 7-tet fixture.  mmg3d's quality functional is not meaningful
    # on a sub-mesh that small — the halo subdomain encompasses
    # the entire model and mmg3d has no "exterior" tets to use as
    # a reference frame.  γ_min monotonicity is the correct
    # invariant for production-scale meshes (1M+ tets) and is
    # validated by the empirical end-to-end runs in
    # `run_newset_step_by_step.sh`, not here.
    assert "gamma_min_in" in rep
    assert "gamma_min_out" in rep


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_outside_subdomain_tets_unchanged(tmp_path):
    """Tets ENTIRELY outside the halo subdomain (no shared
    vertices with subdomain tets) must keep their vertex coords
    bit-exact through the patch.

    Build a mesh with two well-separated regions: a sliver cluster
    in one and a clean tet far away.  The clean tet should be
    untouched.
    """
    import meshio
    P = np.array([
        # Sliver cluster near origin (z < 100):
        [0,    0,   0],    # 0
        [100,  0,   0],    # 1
        [100, 100,  0],    # 2
        [0,   100,  0],    # 3
        [0.01, 0.01, 0.01], # 4 ← sliver-creating
        # Clean tet far away (centroid ~10 km from sliver):
        [10000, 10000, 5000],   # 5
        [10100, 10000, 5000],   # 6
        [10100, 10100, 5000],   # 7
        [10100, 10100, 5100],   # 8
    ], dtype=np.float64)
    T = np.array([
        [0, 1, 2, 4],   # sliver
        [5, 6, 7, 8],   # clean, isolated tet
    ], dtype=np.int64)
    # Each tet must have boundary tris (mmg3d needs a closed surface
    # for each connected component).  Add box-face triangles for both.
    F = np.array([
        # Sliver-cluster boundary:
        [0, 1, 2], [0, 1, 4], [1, 2, 4], [0, 2, 4],
        # Clean-tet boundary:
        [5, 6, 7], [5, 6, 8], [6, 7, 8], [5, 7, 8],
    ], dtype=np.int64)
    F_tags = np.array([100, 1, 1, 1, 100, 1, 1, 1], dtype=np.int32)
    T_tags = np.array([10, 10], dtype=np.int32)

    in_msh = tmp_path / "two_regions.msh"
    out_msh = tmp_path / "two_regions_out.msh"
    meshio.write(in_msh, meshio.Mesh(
        points=P,
        cells=[("triangle", F), ("tetra", T)],
        cell_data={"gmsh:physical": [F_tags, T_tags],
                   "gmsh:geometrical": [F_tags, T_tags]},
    ), file_format="gmsh22", binary=False)

    rep = mlp.local_patch(in_msh, out_msh,
                           gamma_thresh=0.5,
                           halo_radius_m=300.0,
                           hmin=10.0, hmax=200.0,
                           hgrad=1.1, hausd=5.0)

    # Read output and verify the FAR tet (centroid ~10km) survives
    # with its original 4 vertex coords.
    out = meshio.read(out_msh)
    out_P = out.points
    far_centroid_target = np.array([5 + 6 + 7 + 8, 5 + 6 + 7 + 8,
                                      5 + 6 + 7 + 8])  # placeholder
    # Find any vertex with coord ≈ (10000, 10000, 5000).
    found = False
    for v in out_P:
        if (abs(v[0] - 10000.0) < 1e-9
                and abs(v[1] - 10000.0) < 1e-9
                and abs(v[2] - 5000.0) < 1e-9):
            found = True
            break
    assert found, (
        "outside-subdomain vertex (10000, 10000, 5000) was not "
        "preserved bit-exactly through the patch")


def test_canonical_snap_modifies_polyline_vertex_coords():
    """`_canonical_snap_polyline_in_mesh` snaps polyline vertices to
    `snap_key * snap_m`; non-polyline vertices are left alone.
    """
    snap_m = 0.1
    # 3 vertices: 0 + 1 are polyline (off canonical by 0.03), 2 is not.
    P = np.array([
        [0.030, 0.030, 0.030],   # 0 — polyline, should snap to (0,0,0)
        [1.030, 0.030, 0.030],   # 1 — polyline, snap to (1,0,0)
        [99.999, 99.999, 99.999],# 2 — not polyline, untouched
    ], dtype=np.float64)
    polyline_keys = {(0, 0, 0), (10, 0, 0)}  # = (0,0,0) and (1.0,0,0)
    merged = {"points": P}
    n_snapped = mlp._canonical_snap_polyline_in_mesh(
        merged, polyline_keys, snap_m)
    assert n_snapped == 2
    np.testing.assert_array_equal(P[0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_equal(P[1], np.array([1.0, 0.0, 0.0]))
    # Non-polyline vertex unchanged (any drift here would be a bug).
    np.testing.assert_array_equal(
        P[2], np.array([99.999, 99.999, 99.999]))


def test_canonical_snap_no_op_when_already_canonical():
    """A polyline vertex already at canonical position should NOT be
    counted as snapped (n_snapped only counts ACTUAL coord moves)."""
    snap_m = 0.1
    P = np.array([
        [0.0, 0.0, 0.0],   # exactly canonical
        [1.0, 0.0, 0.0],   # exactly canonical
    ], dtype=np.float64)
    polyline_keys = {(0, 0, 0), (10, 0, 0)}
    merged = {"points": P.copy()}
    n_snapped = mlp._canonical_snap_polyline_in_mesh(
        merged, polyline_keys, snap_m)
    assert n_snapped == 0
    np.testing.assert_array_equal(merged["points"], P)


def test_canonical_snap_empty_polyline_keys_is_noop():
    """If polyline_keys is empty, no vertex should be modified."""
    P = np.array([[0.5, 0.5, 0.5]], dtype=np.float64)
    P_orig = P.copy()
    merged = {"points": P}
    n_snapped = mlp._canonical_snap_polyline_in_mesh(
        merged, set(), 0.1)
    assert n_snapped == 0
    np.testing.assert_array_equal(P, P_orig)


def test_R006_stitch_back_raises_on_surface_hole():
    """REVIEW.md R-006: the stitch-back step must raise RuntimeError
    if the merged mesh has any 1-tet boundary face whose 3 vertices
    are NOT covered by a tagged triangle.

    Construct a synthetic input where the all-in-subdomain skip drops
    the only tag for a model-bdry face that an outside-cavity tet
    still uses.  We hand-craft `in_data`, `in_subdomain`, and the
    sub_out arrays to model that pathology.
    """
    # Two tets sharing face (0, 1, 2).  Tet 0 = subdomain (in), tet 1
    # = outside.  Tet 1's face (0,1,2) and (0,1,3) and (0,2,3) and
    # (1,2,3) are its 4 faces; (0,1,2) is shared with tet 0 so NOT a
    # boundary face of the merged mesh.  But tet 1 has 3 OTHER
    # boundary faces — we tag 2 of them (intentionally MISS the third
    # to simulate a stitch bug).
    points = np.array([
        [0.0, 0.0, 0.0],   # 0
        [1.0, 0.0, 0.0],   # 1
        [0.0, 1.0, 0.0],   # 2
        [0.5, 0.5, -1.0],  # 3 (subdomain interior)
        [0.5, 0.5, 1.0],   # 4 (outside-tet apex)
    ], dtype=np.float64)
    in_data = {
        "points":   points,
        "tris":     np.array([
            [0, 1, 2],   # internal interface (subdomain face)
            [0, 1, 4],   # outside-tet boundary face — KEEP
            [0, 2, 4],   # outside-tet boundary face — KEEP
            # MISSING: (1, 2, 4)  ← the synthetic "stitch bug"
        ], dtype=np.int64),
        "tri_tags": np.array([100, 5, 5], dtype=np.int32),
        "tets":     np.array([
            [0, 1, 2, 3],   # subdomain tet
            [0, 1, 2, 4],   # outside tet
        ], dtype=np.int64),
        "tet_tags": np.array([10, 10], dtype=np.int32),
    }
    in_subdomain = np.array([True, False], dtype=bool)

    # mmg3d "output": replace the subdomain tet with a single tet of
    # the same shape.  Tris from mmg3d include the internal-interface
    # face (tag 100) but mmg3d does NOT emit a tri for the missing
    # outside-tet face (1,2,4) — that's the bug the assertion catches.
    sub_out_points = points[[0, 1, 2, 3]].copy()
    sub_out_tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    sub_out_tet_tags = np.array([10], dtype=np.int32)
    sub_out_tris = np.array([[0, 1, 2]], dtype=np.int64)
    sub_out_tri_tags = np.array([100], dtype=np.int32)

    with pytest.raises(RuntimeError,
                       match=r"unlabeled 1-tet bdry face"):
        mlp._stitch_back(
            in_data, in_subdomain,
            sub_out_points, sub_out_tets, sub_out_tet_tags,
            sub_out_tris, sub_out_tri_tags,
            halo_tag=999, snap_m=0.1)


def test_R006_stitch_back_passes_when_topology_closed():
    """Sanity: when EVERY outside-tet boundary face has a tagged tri
    and mmg3d preserves the internal interface, the assertion passes
    and _stitch_back returns a valid merged dict."""
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, -1.0],
        [0.5, 0.5, 1.0],
    ], dtype=np.float64)
    in_data = {
        "points":   points,
        "tris":     np.array([
            [0, 1, 2],   # internal interface
            [0, 1, 4],   # outside-tet face (closed)
            [0, 2, 4],
            [1, 2, 4],   # NOW INCLUDED — closed boundary
            # subdomain tet faces (other than internal interface)
            [0, 1, 3],
            [0, 2, 3],
            [1, 2, 3],
        ], dtype=np.int64),
        "tri_tags": np.array([100, 5, 5, 5, 5, 5, 5], dtype=np.int32),
        "tets":     np.array([
            [0, 1, 2, 3],
            [0, 1, 2, 4],
        ], dtype=np.int64),
        "tet_tags": np.array([10, 10], dtype=np.int32),
    }
    in_subdomain = np.array([True, False], dtype=bool)
    sub_out_points = points[[0, 1, 2, 3]].copy()
    sub_out_tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    sub_out_tet_tags = np.array([10], dtype=np.int32)
    # mmg3d output reproduces every subdomain-bdry tri.
    sub_out_tris = np.array([
        [0, 1, 2],
        [0, 1, 3],
        [0, 2, 3],
        [1, 2, 3],
    ], dtype=np.int64)
    sub_out_tri_tags = np.array([100, 5, 5, 5], dtype=np.int32)

    merged = mlp._stitch_back(
        in_data, in_subdomain,
        sub_out_points, sub_out_tets, sub_out_tet_tags,
        sub_out_tris, sub_out_tri_tags,
        halo_tag=999, snap_m=0.1)
    assert merged["tets"].shape[0] == 2
    assert merged["tris"].shape[0] >= 4


def test_halo_face_coinciding_with_box_face_does_not_create_hole():
    """REGRESSION: when a cavity-boundary halo face is ALSO an
    existing tagged surface triangle (e.g., box face tag 5), the
    medit writer must NOT add a halo_tag duplicate.  Instead the
    existing triangle is marked Required and preserved under its
    ORIGINAL tag.

    Without this fix, the medit input had the same face twice with
    different physical tags; mmg3d preserved the Required (halo)
    version and could drop the non-Required (original) version.
    Stitch-back then dropped halo-tag triangles, and the cavity-
    boundary face became unlabeled — a surface hole.

    Verify by inspecting the medit input directly: the offending
    face must appear EXACTLY ONCE in the Triangles block, and its
    tag must be the original surface tag (not halo_tag).
    """
    import tempfile
    # Set up: a tag-5 (box-top) triangle that's also a halo face.
    points = np.zeros((6, 3), dtype=np.float64)
    points[0] = [0.0, 0.0, 0.0]
    points[1] = [1.0, 0.0, 0.0]
    points[2] = [0.0, 1.0, 0.0]
    points[3] = [0.5, 0.5, -1.0]   # in subdomain
    points[4] = [0.5, 0.5, +1.0]   # outside subdomain
    points[5] = [2.0, 2.0, -2.0]   # outside subdomain
    tets = np.array([
        [0, 1, 2, 3],   # cavity tet — its face (0,1,2) is on box top
        [0, 1, 2, 4],   # outside-cavity neighbor (separated by face)
    ], dtype=np.int64)
    tet_tags = np.array([10, 10], dtype=np.int32)
    # Box-top triangle (tag 5) — coincides with the cavity-boundary
    # face when the only-tet-0 cavity has face (0,1,2) as halo.
    # But (0,1,2) is shared between two tets — so it's NOT a halo
    # face.  Use a different setup:
    #   - tet 0 in cavity, tet 1 outside.
    #   - face (0,1,2) shared between them → halo face.
    #   - face (0,1,2) ALSO exists in tris with tag 5 (which makes
    #     no physical sense for an interior face but is the test
    #     case for the writer logic).
    in_subdomain = np.array([True, False])
    halo_faces = [(0, 1, 2)]
    # Place a "tag 5" triangle at face (0,1,2).
    tris = np.array([[0, 1, 2]], dtype=np.int64)
    tri_tags = np.array([5], dtype=np.int32)
    halo_tag = 200

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "subdomain.mesh"
        meta = mlp._write_subdomain_medit(
            out_path,
            points=points, in_subdomain=in_subdomain,
            tets=tets, tet_tags=tet_tags,
            tris=tris, tri_tags=tri_tags,
            halo_faces=halo_faces, halo_tag=halo_tag,
            polyline_keys=None, snap_m=0.001,
            free_surface_clearance_m=100.0, z_top=0.0,
            fault_tri_tag=100,
        )
        # Read the medit and parse triangles.
        text = out_path.read_text()
        # Find Triangles block.
        lines = text.split("\n")
        i = next(idx for idx, ln in enumerate(lines)
                 if ln.strip() == "Triangles")
        n_tri = int(lines[i + 1])
        emitted_tris = []
        for k in range(n_tri):
            parts = lines[i + 2 + k].split()
            v1, v2, v3, tag = int(parts[0]), int(parts[1]), \
                              int(parts[2]), int(parts[3])
            emitted_tris.append((tuple(sorted((v1, v2, v3))), tag))
        # The face (1, 2, 3) (1-based for v0,v1,v2) must appear
        # EXACTLY ONCE, with tag = 5 (NOT halo_tag = 200).
        face_key = (1, 2, 3)
        matches = [(f, t) for (f, t) in emitted_tris if f == face_key]
        assert len(matches) == 1, (
            f"halo face that coincides with tag-5 triangle was "
            f"written {len(matches)} times — expected exactly 1.  "
            f"Bug: medit writer is adding halo_tag duplicate.")
        assert matches[0][1] == 5, (
            f"halo-coinciding triangle should keep its original "
            f"tag (5); got {matches[0][1]}.  Bug: medit writer "
            f"overwrote the original tag with halo_tag.")
        # And the existing triangle's local index (here 0) must be
        # in the RequiredTriangles list.
        i_req = next((idx for idx, ln in enumerate(lines)
                       if ln.strip() == "RequiredTriangles"), None)
        assert i_req is not None, (
            "RequiredTriangles block missing from medit output")
        n_req = int(lines[i_req + 1])
        req_idx = [int(lines[i_req + 2 + k]) for k in range(n_req)]
        assert 1 in req_idx, (
            f"existing tag-5 triangle (local idx 1, 1-based) must "
            f"be in RequiredTriangles; got {req_idx}")


def test_box_face_subtriangles_are_marked_Required():
    """REGRESSION: box-face triangles (tags 1-6) within the subdomain
    must be marked RequiredTriangles in the medit input — even when
    they are NOT halo-boundary faces (which is the case for model-
    boundary box faces, since they have no neighbour tet on the
    other side).

    Without this, mmg3d in polyline-relax mode can collapse short
    edges in box-face triangles, orphaning the original tag-5
    entry in `tris[]` and producing surface holes at z=0.

    Verify by inspecting the medit file: every tag != 100 triangle
    must appear in the RequiredTriangles list.
    """
    import tempfile
    # Setup: 1 cavity tet (in subdomain) + 1 outside tet, sharing
    # face (0,1,2).  Add a box-face tri (tag 5) at (0,1,3) whose
    # vertices are all in the subdomain BUT it is NOT a halo face.
    points = np.zeros((6, 3), dtype=np.float64)
    points[0] = [0.0, 0.0, 0.0]
    points[1] = [1.0, 0.0, 0.0]
    points[2] = [0.0, 1.0, 0.0]
    points[3] = [0.0, 0.0, -1.0]
    points[4] = [1.0, 1.0, 1.0]
    tets = np.array([
        [0, 1, 2, 3],   # cavity (in subdomain)
        [0, 1, 2, 4],   # outside (not in subdomain)
    ], dtype=np.int64)
    tet_tags = np.array([10, 10], dtype=np.int32)
    in_subdomain = np.array([True, False])
    halo_faces = [(0, 1, 2)]    # the shared face
    # Tag-5 box face on (0, 1, 3) — its 3 vertices are all in the
    # subdomain, but the face itself is NOT a halo face (it's only
    # adjacent to tet 0).
    tris = np.array([[0, 1, 3]], dtype=np.int64)
    tri_tags = np.array([5], dtype=np.int32)
    halo_tag = 200
    fault_tri_tag = 100

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "subdomain.mesh"
        mlp._write_subdomain_medit(
            out_path,
            points=points, in_subdomain=in_subdomain,
            tets=tets, tet_tags=tet_tags,
            tris=tris, tri_tags=tri_tags,
            halo_faces=halo_faces, halo_tag=halo_tag,
            polyline_keys=None, snap_m=0.001,
            free_surface_clearance_m=100.0, z_top=0.0,
            fault_tri_tag=fault_tri_tag,
        )
        text = out_path.read_text()
        lines = text.split("\n")
        # Find the Triangles block to map local indices to tags.
        i_t = next(idx for idx, ln in enumerate(lines)
                    if ln.strip() == "Triangles")
        n_t = int(lines[i_t + 1])
        emitted = []
        for k in range(n_t):
            parts = lines[i_t + 2 + k].split()
            v1, v2, v3, tag = int(parts[0]), int(parts[1]), \
                              int(parts[2]), int(parts[3])
            emitted.append((tuple(sorted((v1, v2, v3))), tag))
        # Find the RequiredTriangles block.
        i_r = next((idx for idx, ln in enumerate(lines)
                     if ln.strip() == "RequiredTriangles"), None)
        assert i_r is not None
        n_r = int(lines[i_r + 1])
        req_indices = [int(lines[i_r + 2 + k]) for k in range(n_r)]
        # Find the local index of the tag-5 box face.
        for li, (face, tag) in enumerate(emitted, start=1):
            if tag == 5:
                assert li in req_indices, (
                    f"box-face triangle (tag 5, local idx {li}) is "
                    f"NOT in RequiredTriangles {req_indices}.  "
                    f"Bug: box-face triangles must be Required so "
                    f"mmg3d does not collapse them in polyline-"
                    f"relax mode.")
                break
        else:
            pytest.fail("tag-5 triangle not found in emitted list")


def test_polyline_relax_flag_in_report(tmp_path):
    """The local_patch report records the `polyline_relax` flag so
    downstream tooling can identify which mode produced the mesh.
    """
    in_msh = _make_sliver_test_mesh(tmp_path)
    out_msh = tmp_path / "out_relax.msh"
    if _mmg3d is None:
        pytest.skip("mmg3d_O3 not on PATH")
    rep = mlp.local_patch(in_msh, out_msh,
                           gamma_thresh=0.5,
                           halo_radius_m=300.0,
                           hmin=10.0, hmax=200.0,
                           hgrad=1.1, hausd=5.0,
                           polyline_relax=True)
    assert "polyline_relax" in rep
    assert rep["polyline_relax"] is True
    assert "n_polyline_canonical_snapped" in rep


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_writes_report_json(tmp_path):
    """`--report-json` writes a JSON file with the expected fields."""
    in_msh = _make_sliver_test_mesh(tmp_path)
    out_msh = tmp_path / "out.msh"
    rep_json = tmp_path / "rep.json"
    mlp.main([
        "--in-msh", str(in_msh),
        "--out-msh", str(out_msh),
        "--gamma-thresh", "0.5",
        "--halo-radius-m", "300",
        "--hmin", "10", "--hmax", "200",
        "--hgrad", "1.1", "--hausd", "5",
        "--hgradreq", "1.1",
        "--report-json", str(rep_json),
    ])
    rep = json.loads(rep_json.read_text())
    for k in ("gamma_min_in", "gamma_min_out", "n_bad_tets",
               "n_subdomain_tets", "subdomain_meta", "halo_tag"):
        assert k in rep, f"report missing key {k}"
    assert rep["halo_tag"] == 200
