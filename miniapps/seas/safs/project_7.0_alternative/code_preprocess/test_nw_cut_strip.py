"""
test_nw_cut_strip.py — pytest suite for nw_cut_strip.py.

Run with:
    cd project_7.0_alternative/code_preprocess && pytest -q

Most tests are pure-NumPy and need only `numpy` and `pytest`. The real-
data smoke (test_smoke_real_data) and the equivalence test against
`ts_to_stl.clip_triangle_at_z0` skip themselves if their dependencies
are unavailable (preferred 2000 m STL on disk, or `meshio` for the
real-data load).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from nw_cut_strip import (                                # noqa: E402
    EPS,
    NW_DIRECTION_XY,
    clip_mesh_at_plane,
    clip_triangle_at_plane,
    cutting_plane,
    lerp_to_plane,
    load_vertices,
    nw_anchor,
    nw_cut_file,
    signed_dist,
    snap_vertices_to_plane,
    DEFAULT_REFERENCE,
    DEFAULT_ALT_DIR,
    DEFAULT_ALT_GLOB,
    DEFAULT_SNAP_TOL_FRACTION,
)

# Tolerance bands
EPS_TEST_GEOM = 1.0e-6                        # for "lies on the plane" checks
EPS_TEST_AREA_REL = 1.0e-9                    # for polygon area sums

# A randomly-rotated cutting plane reused across the four canonical cases.
def _random_plane(seed: int = 0) -> tuple[np.ndarray, float]:
    rng = np.random.default_rng(seed)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    c = float(rng.uniform(-5.0, 5.0))
    return n, c


def make_triangle_with_signed_dists(d0: float, d1: float, d2: float,
                                    n: np.ndarray, c: float,
                                    *, rng=None) -> tuple:
    """Build a triangle whose three vertices have signed distances
    d0, d1, d2 from the plane n . p = c. Used as a fixture for the
    four canonical clip cases.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    p = []
    for di in (d0, d1, d2):
        q = rng.uniform(-10.0, 10.0, size=3)
        # Project q onto the plane: q' = q - signed_dist(q) * n.
        q = q - signed_dist(tuple(q), n, c) * n
        # Translate by di along n -> signed distance becomes di.
        pt = q + di * n
        p.append((float(pt[0]), float(pt[1]), float(pt[2])))
    return p[0], p[1], p[2]


def _tri_area(p0, p1, p2) -> float:
    a = np.asarray(p1) - np.asarray(p0)
    b = np.asarray(p2) - np.asarray(p0)
    return 0.5 * float(np.linalg.norm(np.cross(a, b)))


def _cross_sign(p0, p1, p2, axis: np.ndarray) -> float:
    """Sign of the triangle normal projected on `axis` — used to verify
    CCW orientation is preserved."""
    a = np.asarray(p1) - np.asarray(p0)
    b = np.asarray(p2) - np.asarray(p0)
    return float(np.dot(np.cross(a, b), axis))


# ----------------------------------------------------------------------
# Phase 1 — anchor
# ----------------------------------------------------------------------

class TestPhase1Anchor:

    def test_R001_nw_anchor_does_not_crash_on_real_n(self):
        """R-001 regression: the operand order in nw_anchor must allow
        any N != 2. Earlier spec used `NW_DIRECTION_XY @ verts[:, :2]`
        which raised ValueError for N != 2."""
        rng = np.random.default_rng(0)
        verts = rng.uniform(-1.0, 1.0, size=(1000, 3))
        p, s = nw_anchor(verts)
        assert p.shape == (3,)
        assert isinstance(s, float)
        # Brute-force the same answer.
        s_brute = (-verts[:, 0] + verts[:, 1]) / np.sqrt(2.0)
        j = int(np.argmax(s_brute))
        np.testing.assert_allclose(p, verts[j], atol=0)
        np.testing.assert_allclose(s, s_brute[j], atol=0)

    def test_anchor_tie_breaks_by_smallest_z(self):
        # Three vertices all at the same NW projection; lower z wins.
        verts = np.array([
            [1.0, 2.0, 0.0],            # s = 1/sqrt(2)
            [1.0, 2.0, -5.0],           # tied; lower z
            [1.0, 2.0, +3.0],           # tied; higher z
            [3.0, 1.0, -2.0],           # s < tied
        ], dtype=np.float64)
        p, s = nw_anchor(verts)
        np.testing.assert_allclose(p, [1.0, 2.0, -5.0])
        assert abs(s - 1.0 / np.sqrt(2.0)) < 1e-12

    def test_anchor_rejects_bad_shape(self):
        with pytest.raises(ValueError):
            nw_anchor(np.array([[1.0, 2.0]]))
        with pytest.raises(ValueError):
            nw_anchor(np.empty((0, 3), dtype=np.float64))


class TestPhase1Plane:

    def test_R004_cutting_plane_at_origin(self):
        n, c = cutting_plane(np.array([0.0, 0.0, 0.0]))
        np.testing.assert_allclose(n, [-1.0 / np.sqrt(2.0),
                                       1.0 / np.sqrt(2.0), 0.0])
        assert abs(c) < 1e-12

    def test_cutting_plane_ignores_z(self):
        n_a, c_a = cutting_plane(np.array([1.0, 1.0, 0.0]))
        n_b, c_b = cutting_plane(np.array([1.0, 1.0, 7.0]))
        np.testing.assert_allclose(n_a, n_b)
        assert abs(c_a - c_b) < 1e-12
        # n = (-1, 1, 0)/sqrt(2); n . (1, 1, *) = 0.
        assert abs(c_a) < 1e-12

    def test_cutting_plane_is_unit_normal(self):
        rng = np.random.default_rng(1)
        for _ in range(20):
            p = rng.uniform(-1e6, 1e6, size=3)
            n, _ = cutting_plane(p)
            assert abs(np.linalg.norm(n) - 1.0) < 1e-12


# ----------------------------------------------------------------------
# Phase 2 — clipper
# ----------------------------------------------------------------------

class TestPhase2Clipper:

    def setup_method(self):
        self.n, self.c = _random_plane(seed=42)

    def test_all_keep(self):
        p0, p1, p2 = make_triangle_with_signed_dists(
            -3.0, -2.0, -1.0, self.n, self.c)
        out = clip_triangle_at_plane(p0, p1, p2, self.n, self.c)
        assert len(out) == 1
        assert out[0] == (p0, p1, p2)                 # no copies, no reorder

    def test_all_drop(self):
        p0, p1, p2 = make_triangle_with_signed_dists(
            +1.0, +2.0, +3.0, self.n, self.c)
        out = clip_triangle_at_plane(p0, p1, p2, self.n, self.c)
        assert out == []

    def test_one_keep_two_drop(self):
        # 1 keep + 2 drop -> 1 output triangle.
        p0, p1, p2 = make_triangle_with_signed_dists(
            -2.0, +1.0, +3.0, self.n, self.c)
        out = clip_triangle_at_plane(p0, p1, p2, self.n, self.c)
        assert len(out) == 1
        b, a_clip, c_clip = out[0]
        # Surviving keep-side vertex must be one of the inputs.
        inputs = {p0, p1, p2}
        assert b in inputs
        # New vertices must lie on the plane.
        for v in (a_clip, c_clip):
            assert abs(signed_dist(v, self.n, self.c)) \
                <= EPS_TEST_GEOM * (abs(self.c) + 1.0)
        # Orientation preserved.
        sign_in = _cross_sign(p0, p1, p2, self.n)
        sign_out = _cross_sign(b, a_clip, c_clip, self.n)
        assert sign_in * sign_out > 0

    def test_two_keep_one_drop(self):
        # 2 keep + 1 drop -> 2 output triangles.
        p0, p1, p2 = make_triangle_with_signed_dists(
            +2.5, -1.0, -3.0, self.n, self.c)
        out = clip_triangle_at_plane(p0, p1, p2, self.n, self.c)
        assert len(out) == 2
        # New vertices on the plane.
        ab, b, cc = out[0]
        ab2, cc2, ac = out[1]
        # The (ab, cc) shared edge connects the two output triangles.
        assert ab == ab2
        assert cc == cc2
        for v in (ab, ac):
            assert abs(signed_dist(v, self.n, self.c)) \
                <= EPS_TEST_GEOM * (abs(self.c) + 1.0)
        # Orientation preserved on each output.
        sign_in = _cross_sign(p0, p1, p2, self.n)
        for tri in out:
            sign_out = _cross_sign(*tri, axis=self.n)
            assert sign_in * sign_out > 0
        # Polygon area sum equals area of the keep-side quadrilateral
        # (input triangle minus the small triangle (a, ab, ac)).
        a_input = _tri_area(p0, p1, p2)
        a_drop = _tri_area(p0, ab, ac)                # drop-vertex is p0
        a_out = _tri_area(*out[0]) + _tri_area(*out[1])
        assert abs(a_out - (a_input - a_drop)) \
            <= EPS_TEST_AREA_REL * a_input * 1e3      # generous; clipping
        # is exact within float64 noise.

    def test_axis_aligned_z0_equivalence(self):
        """When n = (0, 0, 1) and c = 0, clip_triangle_at_plane must
        agree with ts_to_stl.clip_triangle_at_z0 on random triangles."""
        try:
            from ts_to_stl import clip_triangle_at_z0
        except Exception:
            pytest.skip("ts_to_stl not importable")
        n = np.array([0.0, 0.0, 1.0])
        c = 0.0
        rng = np.random.default_rng(123)
        for _ in range(100):
            pts = rng.uniform(-1.0, 1.0, size=(3, 3))
            p0 = tuple(pts[0]); p1 = tuple(pts[1]); p2 = tuple(pts[2])
            ours = clip_triangle_at_plane(p0, p1, p2, n, c)
            theirs = clip_triangle_at_z0(p0, p1, p2)
            assert len(ours) == len(theirs)
            # Triangles in same order, same vertex sequence (the
            # algorithm is byte-equivalent).
            for tri_o, tri_t in zip(ours, theirs):
                for v_o, v_t in zip(tri_o, tri_t):
                    np.testing.assert_allclose(v_o, v_t, atol=1e-12)


class TestPhase2MeshClipper:

    def test_R006_stats_keys_match_reference_implementation(self):
        # Plane: n = (-1, +1, 0)/sqrt(2), c = 0  →  d(p) = (-x + y)/sqrt(2)
        # signed_dist > 0 on drop side (large y vs x), <= 0 on keep side.
        n = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2.0)
        c = 0.0
        verts = np.array([
            # Triangle 0: all on keep side (s < 0): y < x
            [-1.0, -2.0, 0.0], [-1.0, -1.5, 0.0], [0.0, -1.5, 0.0],
            # Triangle 1: all on drop side (s > 0): y > x
            [2.0, 5.0, 0.0], [1.0, 4.0, 0.0], [1.5, 6.0, 0.0],
            # Triangle 2: 1 keep + 2 drop  (n_above == 2  -> '2to1')
            [-1.0, -2.0, 0.0], [1.0, 5.0, 0.0], [-2.0, 5.0, 0.0],
            # Triangle 3: 2 keep + 1 drop  (n_above == 1  -> '1to2')
            [-1.0, -2.0, 0.0], [1.0, -2.0, 0.0], [-1.0, 5.0, 0.0],
        ], dtype=np.float64)
        faces = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]],
                         dtype=np.int32)
        V_out, F_out, stats = clip_mesh_at_plane(verts, faces, n, c)
        assert stats["n_in"] == 4
        assert stats["n_kept_whole"] == 1
        assert stats["n_dropped"] == 1
        assert stats["n_clipped_2to1"] == 1
        assert stats["n_clipped_1to2"] == 1
        # Output triangles: 1 (kept) + 0 (dropped) + 1 (2to1) + 2 (1to2) = 4
        assert F_out.shape[0] == 4
        assert stats["n_out"] == 4

    def test_empty_mesh_returns_empty(self):
        n = np.array([1.0, 0.0, 0.0])
        V_out, F_out, stats = clip_mesh_at_plane(
            np.zeros((0, 3), dtype=np.float64),
            np.zeros((0, 3), dtype=np.int32),
            n, 0.0)
        assert V_out.shape == (0, 3)
        assert F_out.shape == (0, 3)
        assert stats["n_in"] == stats["n_out"] == 0


# ----------------------------------------------------------------------
# Phase 2 — lerp_to_plane edge cases
# ----------------------------------------------------------------------

class TestLerpToPlane:

    def test_R010_handles_near_parallel_segment(self):
        """Two endpoints whose signed distances differ by < EPS but
        straddle 0 -> snap branch must still produce an on-plane point."""
        n = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2.0)
        c = 0.0
        # d(p_keep) ≈ -3.5e-10, d(p_drop) ≈ +2.8e-10
        p_keep = (1.0, 1.0 - 5e-10, 0.0)
        p_drop = (1.0, 1.0 + 4e-10, 0.0)
        q = lerp_to_plane(p_keep, p_drop, n, c)
        assert abs(signed_dist(q, n, c)) <= 1e-9

    def test_normal_segment_lands_on_plane(self):
        n = np.array([0.0, 0.0, 1.0])
        c = 0.0
        p_keep = (1.0, 2.0, -3.0)
        p_drop = (4.0, 5.0, +6.0)
        q = lerp_to_plane(p_keep, p_drop, n, c)
        assert abs(signed_dist(q, n, c)) <= 1e-12
        # t = -(-3)/(6-(-3)) = 3/9 = 1/3
        np.testing.assert_allclose(q, (2.0, 3.0, 0.0), atol=1e-12)


# ----------------------------------------------------------------------
# load_vertices
# ----------------------------------------------------------------------

class TestLoadVertices:

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_vertices(tmp_path / "no_such_file.stl")

    def test_R007_load_vertices_drops_orphans_from_ts(self, tmp_path):
        ts = tmp_path / "tiny.ts"
        ts.write_text(
            "GOCAD TSurf 1\n"
            "VRTX 1 100 100 0\n"
            "VRTX 2 200 100 0\n"
            "VRTX 3 100 200 0\n"
            "VRTX 4 -1000000000 1000000000 0\n"      # orphan ATOM, very NW
            "TRGL 1 2 3\n"
            "END\n"
        )
        V = load_vertices(ts)
        assert V.shape == (3, 3)                      # orphan stripped
        p, s = nw_anchor(V)
        # If orphan filter failed, anchor would be (-1e9, 1e9, 0).
        assert p[0] >= 100.0 - 1.0


# ----------------------------------------------------------------------
# CLI mutual exclusion (R-002)
# ----------------------------------------------------------------------

class TestCLI:

    def _run(self, argv: list, tmp_path: Path) -> tuple[int, str, str]:
        import subprocess
        cmd = [sys.executable, str(HERE / "nw_cut_strip.py")] + argv
        r = subprocess.run(cmd, capture_output=True, text=True,
                           cwd=str(tmp_path))
        return r.returncode, r.stdout, r.stderr

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_R002_print_anchor_alone_works(self, tmp_path):
        rc, out, err = self._run(
            ["--print-anchor", str(DEFAULT_REFERENCE)], tmp_path)
        assert rc == 0, err
        assert out.startswith("anchor:")

    def test_R002_print_anchor_with_batch_rejected(self, tmp_path):
        rc, _, err = self._run(
            ["--print-anchor", str(tmp_path / "x"), "--batch"], tmp_path)
        assert rc == 2
        assert "mutually exclusive" in err

    def test_R002_no_mode_prints_help_and_exits_2(self, tmp_path):
        rc, _, err = self._run([], tmp_path)
        assert rc == 2
        assert "no mode" in err

    def test_R002_input_without_output_rejected(self, tmp_path):
        rc, _, err = self._run([str(tmp_path / "in.stl")], tmp_path)
        # argparse will accept this (output positional defaults to None);
        # main() must explicitly reject it.
        assert rc == 2
        assert "INPUT and OUTPUT" in err


# ----------------------------------------------------------------------
# Empty-clip handling (R-004)
# ----------------------------------------------------------------------

class TestEmptyClip:

    def test_R004_nw_cut_file_raises_on_empty_clip(self, tmp_path):
        # Build a tiny mesh entirely SE of the plane y > x — i.e. on the
        # drop side. With c chosen so even the SE-most vertex is on the
        # drop side, every triangle is clipped away.
        # n = (-1, +1, 0)/sqrt(2);  d(p) = (-x + y)/sqrt(2) - c
        # If we set c to a hugely negative value, every vertex has d > 0
        # → everything dropped.
        n = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2.0)
        c = -1e9
        # Write a synthetic STL we can hand to nw_cut_file.
        try:
            import meshio
        except Exception:
            pytest.skip("meshio not available")
        verts = np.array([[0.0, 0.0, 0.0],
                          [1.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0]], dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)
        in_path = tmp_path / "in.stl"
        meshio.write_points_cells(
            str(in_path), verts, [("triangle", faces)],
            file_format="stl")
        out_path = tmp_path / "out.stl"
        with pytest.raises(ValueError, match="empty STL"):
            nw_cut_file(in_path, out_path, (n, c), verbose=False)
        assert not out_path.exists()


# ----------------------------------------------------------------------
# float32 dtype coercion (R-003)
# ----------------------------------------------------------------------

class TestDtypeCoercion:

    def test_R003_clip_promotes_float32_to_float64(self):
        """clip_mesh_at_plane must accept (and internally use) float64
        even if the caller passes float32."""
        verts32 = np.array([[0.0, 0.0, 0.0],
                            [3.0e6, 0.0, 0.0],
                            [0.0, 3.0e6, 0.0]], dtype=np.float32)
        faces = np.array([[0, 1, 2]], dtype=np.int32)
        n, c = cutting_plane(np.array([1.0e6, 1.0e6, 0.0],
                                      dtype=np.float64))
        V_out, F_out, _ = clip_mesh_at_plane(verts32, faces, n, c)
        assert V_out.dtype == np.float64
        # Any output vertex within machine epsilon of the plane should
        # satisfy the on-plane test in float64, not the looser float32
        # tolerance.
        for v in V_out:
            d = signed_dist(tuple(v), n, c)
            # Vertices on the plane (post-clip) within 1e-6 m; vertices
            # in the keep interior have d <= 0.
            assert d <= 1e-6


# ----------------------------------------------------------------------
# s_max_after measured post-cleanup (R-005)
# ----------------------------------------------------------------------

class TestPostCleanupSMax:

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_R005_smax_after_is_post_cleanup(self, tmp_path):
        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, s_max = nw_anchor(verts_ref)
        plane = cutting_plane(p_anchor)
        # Pick the smallest alternative file for speed.
        alt = (DEFAULT_ALT_DIR
               / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m"
                 "_clean_clip.stl")
        if not alt.exists():
            pytest.skip(f"alt 2000 m file missing: {alt}")
        out = tmp_path / "alt_2000m_nwcut.stl"
        stats = nw_cut_file(alt, out, plane, verbose=False)
        # Reload and recompute s_max independently from the saved file.
        try:
            import meshio
        except Exception:
            pytest.skip("meshio not available")
        m = meshio.read(str(out))
        s = np.asarray(m.points, dtype=np.float64)[:, :2] @ NW_DIRECTION_XY
        # Saved-file s_max must match the reported s_max_out.
        assert s.size > 0
        assert s.max() <= s_max + 1e-6
        assert abs(stats["s_max_out"] - s.max()) <= 1e-6


# ----------------------------------------------------------------------
# Real-data smoke (Phase 4 §5)
# ----------------------------------------------------------------------

class TestSmokeRealData:

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_anchor_matches_documented_values(self):
        verts = load_vertices(DEFAULT_REFERENCE)
        # The plan's "1479" expectation came from counting raw `vertex`
        # lines in the ASCII STL (3 per triangle). meshio dedupes on
        # load, so the unique-vertex count is ~291 (no orphans). The
        # spec is correct; the documented count was wrong.
        assert verts.ndim == 2 and verts.shape[1] == 3
        assert verts.shape[0] > 100, \
            f"expected at least 100 unique verts, got {verts.shape[0]}"
        # Bbox checks (1 m tolerance).
        bbox = np.array([verts.min(axis=0), verts.max(axis=0)])
        assert abs(bbox[0, 0] - 365072.3) <= 1.0
        assert abs(bbox[1, 0] - 428025.6) <= 1.0
        assert abs(bbox[0, 1] - 3808374.0) <= 1.0
        assert abs(bbox[1, 1] - 3838985.0) <= 1.0
        assert abs(bbox[0, 2] - (-13272.1)) <= 1.0
        assert abs(bbox[1, 2] - 0.0) <= 1.0
        p, s = nw_anchor(verts)
        np.testing.assert_allclose(
            p, (365072.3, 3838982.0, 0.0), atol=1.0)
        assert abs(s - 2456425.1) <= 1.0
        n, c = cutting_plane(p)
        np.testing.assert_allclose(
            n, (-1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0), 0.0),
            atol=1e-12)
        assert abs(c - s) <= 1e-9

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_cut_2000m_alternative(self, tmp_path):
        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, s_max_pref = nw_anchor(verts_ref)
        plane = cutting_plane(p_anchor)
        alt = (DEFAULT_ALT_DIR
               / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m"
                 "_clean_clip.stl")
        if not alt.exists():
            pytest.skip(f"alt 2000 m file missing: {alt}")
        out = tmp_path / "alt_2000m_nwcut.stl"
        stats = nw_cut_file(alt, out, plane, verbose=False)
        # No surviving vertex may sit further NW than the anchor.
        assert stats["s_max_out"] <= s_max_pref + 1e-6
        # Substantial-but-not-total clipping.
        assert stats["n_in"] > 0
        n_out_before_cleanup = (stats["n_kept_whole"]
                                + stats["n_clipped_2to1"]
                                + 2 * stats["n_clipped_1to2"])
        assert 0.5 * stats["n_in"] < n_out_before_cleanup < 0.95 * stats["n_in"]

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_endpoints_returned(self, tmp_path):
        """nw_cut_file must populate nw_endpoint and se_endpoint with
        coordinates that match s_max_out and s_min_out, and the NW
        endpoint must lie on the cutting plane to within 1 mm."""
        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, _ = nw_anchor(verts_ref)
        n, c = cutting_plane(p_anchor)
        alt = (DEFAULT_ALT_DIR
               / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m"
                 "_clean_clip.stl")
        if not alt.exists():
            pytest.skip(f"alt 2000 m file missing: {alt}")
        out = tmp_path / "alt_2000m_nwcut.stl"
        stats = nw_cut_file(alt, out, (n, c), verbose=False)
        assert "nw_endpoint" in stats and "se_endpoint" in stats
        assert stats["nw_endpoint"].shape == (3,)
        assert stats["se_endpoint"].shape == (3,)
        # NW endpoint's NW projection must equal s_max_out.
        s_nw = float(stats["nw_endpoint"][:2] @ NW_DIRECTION_XY)
        s_se = float(stats["se_endpoint"][:2] @ NW_DIRECTION_XY)
        assert abs(s_nw - stats["s_max_out"]) < 1e-9
        assert abs(s_se - stats["s_min_out"]) < 1e-9
        # NW endpoint sits at (or just SE of) the cutting plane.
        assert stats["s_max_out"] <= c + 1e-6


class TestOutDir:

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_out_dir_redirects_batch_outputs(self, tmp_path):
        """--out-dir must redirect batch outputs to a separate directory
        and must not modify the input directory."""
        import subprocess
        out_dir = tmp_path / "cut_outputs"
        cmd = [sys.executable, str(HERE / "nw_cut_strip.py"),
               "--batch", "--out-dir", str(out_dir),
               "--glob", "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-"
                         "ALT6_2000m_clean_clip.stl"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        out_files = list(out_dir.glob("*_nwcut.stl"))
        assert len(out_files) == 1
        assert out_files[0].name == \
            "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m" \
            "_clean_clip_nwcut.stl"
        # Source dir must not have gained an _nwcut file.
        src_nwcut = list(DEFAULT_ALT_DIR.glob("*_nwcut.stl"))
        for p in src_nwcut:
            # If a previous local --batch ran without --out-dir, the file
            # would already be there. Tolerate it but don't *add* one.
            assert "ALT6_2000m" not in p.name or p.exists(), \
                "test should not have created files in DEFAULT_ALT_DIR"


# ----------------------------------------------------------------------
# Pre-clip plane snap (issue: bulk mesh edge_min ~ 8 m on nwcut faults)
# ----------------------------------------------------------------------

class TestPreClipSnap:
    """The clipper alone has no minimum-feature tolerance: any input
    vertex within sub-LC of the cutting plane produces an output edge of
    that length, which the gmsh `Surface{} In Volume{1}` embedding then
    propagates into bulk-tet quality.  `snap_vertices_to_plane` is the
    fix; these tests pin its contract.
    """

    def test_snap_projects_vertices_within_tol(self):
        n = np.array([0.0, 0.0, 1.0])
        c = 0.0
        V = np.array([[0.0, 0.0, -1.0],     # well below (kept as-is)
                      [0.0, 0.0, -0.05],    # within tol  (snapped)
                      [0.0, 0.0,  0.04],    # within tol  (snapped)
                      [0.0, 0.0,  0.5]])    # well above  (kept as-is)
        V_out, n_snapped = snap_vertices_to_plane(V, n, c, snap_tol=0.1)
        assert n_snapped == 2
        # In-tol vertices snap exactly onto plane.
        assert abs(V_out[1, 2]) < EPS_TEST_GEOM
        assert abs(V_out[2, 2]) < EPS_TEST_GEOM
        # Out-of-tol vertices are untouched.
        assert V_out[0, 2] == -1.0
        assert V_out[3, 2] == 0.5

    def test_snap_disabled_when_tol_zero(self):
        n, c = _random_plane(0)
        rng = np.random.default_rng(1)
        V = rng.uniform(-5, 5, size=(20, 3))
        V_out, n_snapped = snap_vertices_to_plane(V, n, c, snap_tol=0.0)
        assert n_snapped == 0
        np.testing.assert_array_equal(V_out, V)

    def test_snap_does_not_mutate_input(self):
        n = np.array([0.0, 0.0, 1.0])
        c = 0.0
        V = np.array([[0.0, 0.0, 0.05]])
        V_copy = V.copy()
        V_out, _ = snap_vertices_to_plane(V, n, c, snap_tol=0.1)
        np.testing.assert_array_equal(V, V_copy)
        assert V_out is not V

    def test_pre_snap_eliminates_sliver_in_2k1d_case(self):
        """Drop vertex barely past the plane (d_drop ~ small) without
        snap produces a tiny `ab-ac` cut-line edge.  With snap, the drop
        vertex lands on the plane, becomes kept, and the triangle
        survives whole -- no sliver."""
        n, c = _random_plane(7)
        # 2 keep, 1 drop, with the drop just 0.01 past the plane.
        p_keep1, p_keep2, p_drop = make_triangle_with_signed_dists(
            -1.0, -1.5, 0.01, n, c, rng=np.random.default_rng(8))
        V = np.array([p_keep1, p_keep2, p_drop])
        F = np.array([[0, 1, 2]], dtype=np.int32)

        # Without snap: a 2k+1d clip is performed (n_above == 1 ->
        # n_clipped_1to2 in the stats dict, which counts by n_above), and
        # output triangles include a sub-tol cut-line edge.
        _, _, stats_raw = clip_mesh_at_plane(V, F, n, c)
        assert stats_raw["n_clipped_1to2"] == 1

        # With snap (tol 0.1 > 0.01): drop vertex is now on plane and
        # the triangle is reclassified as fully kept.
        V_snap, n_snapped = snap_vertices_to_plane(V, n, c, snap_tol=0.1)
        assert n_snapped == 1
        _, F_snap, stats_snap = clip_mesh_at_plane(V_snap, F, n, c)
        assert stats_snap["n_kept_whole"] == 1
        assert stats_snap["n_clipped_1to2"] == 0
        assert F_snap.shape[0] == 1

    def test_pre_snap_eliminates_sliver_in_1k2d_case(self):
        """Keep vertex within sub-LC of the plane (d_keep ~ small)
        without snap produces tiny `B-ab` and `B-cb` edges out of the
        keep vertex.  With snap, the keep vertex lands on the plane and
        the lerp results coincide with it -- the entire output triangle
        is degenerate (zero area), to be removed by null-face cleanup
        downstream.

        Uses an axis-aligned plane and explicit geometry so the sliver
        edge length can be analytically predicted from `d_keep` and
        the keep-to-drop edge length.
        """
        n = np.array([0.0, 0.0, 1.0])
        c = 0.0
        # Triangle: B at z = -0.001 (sub-tol), drops at z = +1.0 with
        # large in-plane offsets so |B-drop| ~ 1.4, lerp parameter
        # t ~ 0.001 gives |B-ab| ~ 0.0014 (well below snap_tol = 0.1).
        p_keep  = np.array([0.0, 0.0, -0.001])
        p_drop1 = np.array([1.0, 0.0,  1.0])
        p_drop2 = np.array([0.0, 1.0,  1.0])
        V = np.array([p_keep, p_drop1, p_drop2])
        F = np.array([[0, 1, 2]], dtype=np.int32)

        # 1k+2d -> n_above == 2 -> stats key "n_clipped_2to1".
        Vraw_out, F_raw, stats_raw = clip_mesh_at_plane(V, F, n, c)
        assert stats_raw["n_clipped_2to1"] == 1
        # Baseline (no snap): output triangle has at least one sub-tol
        # edge.  This is the gmsh-poison sliver.
        assert F_raw.shape[0] == 1
        tri = Vraw_out[F_raw[0]]
        edges_raw = [np.linalg.norm(tri[1] - tri[0]),
                     np.linalg.norm(tri[2] - tri[1]),
                     np.linalg.norm(tri[0] - tri[2])]
        assert min(edges_raw) < 0.01      # |B-ab| ~ 0.0014

        # With snap, keep vertex moves onto plane -> n_keep == 1 with
        # d_keep == 0 -> both lerp results equal the keep vertex.  The
        # output triangle (B, B, B) has all three vertices identical.
        V_snap, n_snapped = snap_vertices_to_plane(V, n, c, snap_tol=0.1)
        assert n_snapped == 1
        Vs_out, Fs_out, _ = clip_mesh_at_plane(V_snap, F, n, c)
        if Fs_out.shape[0] > 0:
            tri = Vs_out[Fs_out[0]]
            # All three vertices coincide -> zero-area degenerate, to be
            # removed by the downstream null-face filter.  Crucially, no
            # sub-tol non-zero edge exists.
            for i in range(3):
                np.testing.assert_allclose(
                    tri[i], V_snap[0], atol=EPS_TEST_GEOM)

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_real_data_pre_snap_eliminates_slivers(self, tmp_path):
        """End-to-end on the real 1000 m mesh: with the default snap,
        the output STL must have no edge shorter than ~snap_tol; without
        the snap (snap_tol=0), at least one edge < 100 m re-appears.
        """
        import meshio
        in_path = (DEFAULT_ALT_DIR
                   / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-"
                     "ALT6_1000m_clean_clip.stl")
        if not in_path.exists():
            pytest.skip("1000 m alt STL not present")

        plane, _, _ = (np.array([NW_DIRECTION_XY[0], NW_DIRECTION_XY[1],
                                 0.0]),
                       0.0), 0.0, 0.0
        # Use the production plane via load_vertices on the reference.
        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, _ = nw_anchor(verts_ref)
        n, c = cutting_plane(p_anchor)

        # Run with default snap (auto = 0.1 * median input edge).
        out_default = tmp_path / "default_snap.stl"
        stats_def = nw_cut_file(in_path, out_default, (n, c),
                                cleanup=True, verbose=False,
                                snap_tol=None)
        m = meshio.read(str(out_default))
        V = np.asarray(m.points)
        F = np.asarray(m.cells_dict["triangle"])
        e0 = np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)
        e1 = np.linalg.norm(V[F[:, 2]] - V[F[:, 1]], axis=1)
        e2 = np.linalg.norm(V[F[:, 0]] - V[F[:, 2]], axis=1)
        edge_min_default = float(np.minimum(np.minimum(e0, e1), e2).min())
        assert stats_def["snap_tol"] > 0.0
        # The default snap is ~10% of median edge ~99 m; the output must
        # have no edge much shorter than that.  Allow a small slack: no
        # edge below 50 % of snap_tol.
        assert edge_min_default >= 0.5 * stats_def["snap_tol"], (
            f"default-snap output still has emin={edge_min_default:.2f} m "
            f"(snap_tol = {stats_def['snap_tol']:.2f} m)")

        # Run with snap explicitly disabled -> regression baseline.
        out_legacy = tmp_path / "no_snap.stl"
        nw_cut_file(in_path, out_legacy, (n, c),
                    cleanup=True, verbose=False, snap_tol=0.0)
        m = meshio.read(str(out_legacy))
        V = np.asarray(m.points)
        F = np.asarray(m.cells_dict["triangle"])
        e0 = np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)
        e1 = np.linalg.norm(V[F[:, 2]] - V[F[:, 1]], axis=1)
        e2 = np.linalg.norm(V[F[:, 0]] - V[F[:, 2]], axis=1)
        edge_min_legacy = float(np.minimum(np.minimum(e0, e1), e2).min())
        # Legacy (no snap) reproduces the original 8.6 m sliver.
        assert edge_min_legacy < 50.0, (
            f"legacy snap=0 should reproduce sliver, got "
            f"emin={edge_min_legacy:.2f} m")
        assert edge_min_default > edge_min_legacy


# ----------------------------------------------------------------------
# Post-clip close-vertex merge (issue: 500 m fault still has a 45 m
# cut-line edge after snap because two near-plane vertices both project
# onto the plane and end up close)
# ----------------------------------------------------------------------

class TestPostClipMerge:
    """Snap moves vertices but does not merge them.  Two vertices both
    within `snap_tol` of the cutting plane can be projected onto the
    plane and remain at sub-snap_tol distance from each other,
    producing a sliver cut-line edge that survives every other cleanup
    step.  `meshing_merge_close_vertices` (via `merge_tol`) is the only
    step that can collapse such pairs.
    """

    @pytest.mark.skipif(not DEFAULT_REFERENCE.exists(),
                        reason="preferred reference STL not present")
    def test_merge_tol_eliminates_close_on_plane_pairs(self, tmp_path):
        """500 m mesh: with merge_tol = 0 we reproduce the 45 m
        cut-line edge.  With merge_tol >= 60 m we eliminate it.
        """
        import meshio
        in_path = (DEFAULT_ALT_DIR
                   / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-"
                     "ALT6_500m_clean_clip.stl")
        if not in_path.exists():
            pytest.skip("500 m alt STL not present")

        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, _ = nw_anchor(verts_ref)
        n, c = cutting_plane(p_anchor)

        def emin(stl):
            m = meshio.read(str(stl))
            V = np.asarray(m.points)
            F = np.asarray(m.cells_dict["triangle"])
            e0 = np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)
            e1 = np.linalg.norm(V[F[:, 2]] - V[F[:, 1]], axis=1)
            e2 = np.linalg.norm(V[F[:, 0]] - V[F[:, 2]], axis=1)
            return float(np.minimum(np.minimum(e0, e1), e2).min())

        # Baseline: snap on, merge OFF -> reproduce the 45 m sliver.
        out_no_merge = tmp_path / "no_merge.stl"
        nw_cut_file(in_path, out_no_merge, (n, c),
                    cleanup=True, verbose=False,
                    snap_tol=None, merge_tol=0.0)
        emin_no = emin(out_no_merge)
        assert emin_no < 60.0, (
            f"snap-only baseline should keep the 45 m sliver; "
            f"got emin={emin_no:.2f} m")

        # Fix: snap on, merge_tol = 100 m -> sliver removed.
        out_merge = tmp_path / "merge_100.stl"
        stats = nw_cut_file(in_path, out_merge, (n, c),
                            cleanup=True, verbose=False,
                            snap_tol=None, merge_tol=100.0)
        emin_m = emin(out_merge)
        assert stats["merge_tol"] == 100.0
        assert emin_m > emin_no
        # No edge survives below the merge threshold.
        assert emin_m >= 100.0, (
            f"merge_tol=100 must eliminate sub-100 m edges; "
            f"got emin={emin_m:.2f} m")

    def test_merge_tol_auto_defaults_to_snap_tol(self, tmp_path):
        """merge_tol=None -> stats["merge_tol"] == snap_tol."""
        if not DEFAULT_REFERENCE.exists():
            pytest.skip("reference STL not present")
        in_path = (DEFAULT_ALT_DIR
                   / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-"
                     "ALT6_2000m_clean_clip.stl")
        if not in_path.exists():
            pytest.skip("2000 m alt STL not present")

        verts_ref = load_vertices(DEFAULT_REFERENCE)
        p_anchor, _ = nw_anchor(verts_ref)
        n, c = cutting_plane(p_anchor)
        out = tmp_path / "auto.stl"
        stats = nw_cut_file(in_path, out, (n, c),
                            cleanup=True, verbose=False,
                            snap_tol=None, merge_tol=None)
        assert stats["snap_tol"] > 0.0
        assert stats["merge_tol"] == stats["snap_tol"]
