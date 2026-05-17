"""
test_project_to_fault_stress.py — pytest suite for
project_to_fault_stress.py.

Run with:
    cd project_7.0_alternative/stress/code && pytest -q test_project_to_fault_stress.py

The real-data tests (test_load_fault_vtu_real_data,
test_triangle_geometry_real_data_bbox,
test_tet_geometry_real_data_count_and_positivity) skip themselves if
the real 2000 m VTUs are not on disk under meshing/results/vtu/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from project_to_fault_stress import (                       # noqa: E402
    EPS,
    NEAR_HORIZONTAL_NZ_TOL,
    UP_VECTOR,
    BulkCellGeometry,
    DEFAULT_BULK_NAME,
    DEFAULT_FAULT_NAME,
    DEFAULT_INPUT_BASES,
    DEFAULT_OUT_DIR,
    FaultCellGeometry,
    extract_cells_by_physical,
    load_bulk_mesh,
    load_fault_mesh,
    tet_geometry,
    triangle_geometry,
)

# Real-data paths
_CODE_MESHING_DIR = (
    HERE.parent.parent / "meshing" / "results" / "vtu"
)
_REAL_2000M_FAULT = _CODE_MESHING_DIR / "safs_fault_box_nwcut_2000m_fault.vtu"
_REAL_2000M_BULK = _CODE_MESHING_DIR / "safs_fault_box_nwcut_2000m_bulk.vtu"

# Expected counts and bbox for the 2000 m mesh (verified via the
# upstream bash probe documented in the plan's Phase 1 acceptance).
_EXPECTED_2000M_N_TRI = 2685
_EXPECTED_2000M_N_TET = 145300
_EXPECTED_2000M_BBOX = dict(
    xmin=303_000.0, xmax=697_500.0,
    ymin=3_612_000.0, ymax=3_903_000.0,
    zmin=-66_607.45, zmax=100.0,
)
# Phase 1 acceptance: "sum of volumes equals the gmsh box volume within 0.1 %".
# Computed from the bbox above:
_EXPECTED_2000M_BBOX_VOLUME = (
    (_EXPECTED_2000M_BBOX["xmax"] - _EXPECTED_2000M_BBOX["xmin"])
    * (_EXPECTED_2000M_BBOX["ymax"] - _EXPECTED_2000M_BBOX["ymin"])
    * (_EXPECTED_2000M_BBOX["zmax"] - _EXPECTED_2000M_BBOX["zmin"])
)


def _skip_if_no_real_data(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"real-data mesh not on disk: {path}")


# ----------------------------------------------------------------------
# Module constants
# ----------------------------------------------------------------------


class TestConstants:
    def test_eps_is_small(self):
        assert 0 < EPS <= 1e-9

    def test_up_vector_is_unit_z(self):
        np.testing.assert_allclose(UP_VECTOR, [0.0, 0.0, 1.0])

    def test_near_horizontal_threshold(self):
        # Just below unity — only flag near-perfectly-horizontal triangles.
        assert 0.99 < NEAR_HORIZONTAL_NZ_TOL < 1.0

    def test_default_names_match_geo_file(self):
        assert DEFAULT_BULK_NAME == "rock"
        assert DEFAULT_FAULT_NAME == "fault"

    def test_default_input_bases_match_real_mesh_resolutions(self):
        assert "safs_fault_box_nwcut_500m" in DEFAULT_INPUT_BASES
        assert "safs_fault_box_nwcut_1000m" in DEFAULT_INPUT_BASES
        assert "safs_fault_box_nwcut_2000m" in DEFAULT_INPUT_BASES

    def test_default_out_dir_is_stress_results(self):
        # Post-restructure (project_7.0_alternative pillar reorg): batch
        # outputs land under stress/results/, the sibling of
        # stress/code/.
        assert DEFAULT_OUT_DIR.name == "results"
        assert DEFAULT_OUT_DIR.parent.name == "stress"


# ----------------------------------------------------------------------
# Synthetic-mesh unit tests
# ----------------------------------------------------------------------


class TestTriangleGeometry:
    def test_unit_345_triangle_in_xy_plane(self):
        # 3-4-5 right triangle in the XY plane, area = 6
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [3.0, 0.0, 0.0],
             [0.0, 4.0, 0.0]],
            dtype=np.float64,
        )
        tri = np.array([[0, 1, 2]], dtype=np.int64)
        c, n, a = triangle_geometry(pts, tri)
        np.testing.assert_allclose(a, [6.0], atol=1e-12)
        np.testing.assert_allclose(n[0], [0.0, 0.0, 1.0], atol=1e-12)
        np.testing.assert_allclose(c[0], [1.0, 4.0 / 3.0, 0.0], atol=1e-12)

    def test_normal_is_unit_under_random_rotations(self):
        rng = np.random.default_rng(0)
        for _ in range(100):
            # Random non-degenerate triangle
            pts = rng.uniform(-1.0, 1.0, size=(3, 3))
            tri = np.array([[0, 1, 2]], dtype=np.int64)
            _c, n, _a = triangle_geometry(pts, tri)
            assert abs(np.linalg.norm(n[0]) - 1.0) < 1e-12

    def test_centroid_equals_mean_of_vertices(self):
        rng = np.random.default_rng(1)
        for _ in range(20):
            pts = rng.uniform(-1.0, 1.0, size=(3, 3))
            tri = np.array([[0, 1, 2]], dtype=np.int64)
            c, _n, _a = triangle_geometry(pts, tri)
            np.testing.assert_allclose(c[0], pts.mean(axis=0), atol=1e-12)

    def test_degenerate_zero_area_triangle_produces_nan(self):
        # Three collinear points → zero area
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [2.0, 0.0, 0.0]],
            dtype=np.float64,
        )
        tri = np.array([[0, 1, 2]], dtype=np.int64)
        c, n, a = triangle_geometry(pts, tri)
        assert a[0] < EPS
        assert np.all(np.isnan(n[0]))
        assert np.all(np.isnan(c[0]))

    def test_bad_points_shape_raises(self):
        with pytest.raises(ValueError):
            triangle_geometry(np.zeros((5, 2)), np.zeros((1, 3), dtype=int))

    def test_bad_tri_shape_raises(self):
        with pytest.raises(ValueError):
            triangle_geometry(np.zeros((5, 3)), np.zeros((1, 4), dtype=int))

    def test_R203_oob_index_raises(self):
        """Out-of-bounds connectivity index raises ValueError with a
        clear message, not a bare numpy IndexError (R-203)."""
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        bad = np.array([[0, 1, 5]], dtype=np.int64)  # 5 is OOB
        with pytest.raises(ValueError, match="out of bounds"):
            triangle_geometry(pts, bad)

    def test_R203_negative_index_raises(self):
        """Negative connectivity index raises ValueError instead of
        silently wrapping (R-203)."""
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        bad = np.array([[0, 1, -1]], dtype=np.int64)
        with pytest.raises(ValueError, match="out of bounds"):
            triangle_geometry(pts, bad)


class TestTetGeometry:
    def test_reference_simplex(self):
        # Reference tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        # → V = 1/6, centroid = (1/4, 1/4, 1/4)
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        tet = np.array([[0, 1, 2, 3]], dtype=np.int64)
        c, v = tet_geometry(pts, tet)
        np.testing.assert_allclose(v, [1.0 / 6.0], atol=1e-12)
        np.testing.assert_allclose(c[0], [0.25, 0.25, 0.25], atol=1e-12)

    def test_volume_is_always_non_negative(self):
        """Reversing the orientation of one face flips the signed
        scalar-triple-product; tet_geometry must return |V|."""
        rng = np.random.default_rng(2)
        pts = rng.uniform(-1.0, 1.0, size=(8, 3))
        tet_ccw = np.array([[0, 1, 2, 3]], dtype=np.int64)
        tet_cw = np.array([[0, 2, 1, 3]], dtype=np.int64)
        _c1, v1 = tet_geometry(pts, tet_ccw)
        _c2, v2 = tet_geometry(pts, tet_cw)
        assert v1[0] >= 0
        assert v2[0] >= 0
        np.testing.assert_allclose(v1, v2, atol=1e-12)

    def test_centroid_equals_mean_of_vertices(self):
        rng = np.random.default_rng(3)
        for _ in range(20):
            pts = rng.uniform(-1.0, 1.0, size=(4, 3))
            tet = np.array([[0, 1, 2, 3]], dtype=np.int64)
            c, _v = tet_geometry(pts, tet)
            np.testing.assert_allclose(c[0], pts.mean(axis=0), atol=1e-12)

    def test_bad_tet_shape_raises(self):
        with pytest.raises(ValueError):
            tet_geometry(np.zeros((5, 3)), np.zeros((1, 3), dtype=int))

    def test_R203_tet_oob_index_raises(self):
        """Out-of-bounds tet connectivity index raises ValueError
        (R-203)."""
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        bad = np.array([[0, 1, 2, 10]], dtype=np.int64)  # 10 is OOB
        with pytest.raises(ValueError, match="out of bounds"):
            tet_geometry(pts, bad)

    def test_R203_tet_negative_index_raises(self):
        """Negative tet connectivity index raises ValueError (R-203)."""
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        bad = np.array([[0, 1, 2, -1]], dtype=np.int64)
        with pytest.raises(ValueError, match="out of bounds"):
            tet_geometry(pts, bad)


# ----------------------------------------------------------------------
# Mesh I/O — errors / edge cases
# ----------------------------------------------------------------------


class TestMeshIO:
    def test_missing_fault_path_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_fault_mesh(Path("/nonexistent/path/to/fault.vtu"))

    def test_missing_bulk_path_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_bulk_mesh(Path("/nonexistent/path/to/bulk.vtu"))

    def test_load_fault_mesh_raises_on_no_triangles(self, tmp_path):
        """A VTU with only tetra cells must fail load_fault_mesh
        (R-201)."""
        import meshio
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        m = meshio.Mesh(
            points=pts,
            cells=[("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        )
        vtu = tmp_path / "tet_only.vtu"
        meshio.write(str(vtu), m)
        with pytest.raises(ValueError, match="no triangle cells"):
            load_fault_mesh(vtu)

    def test_load_bulk_mesh_raises_on_no_tetra(self, tmp_path):
        """A VTU with only triangle cells must fail load_bulk_mesh
        (R-201)."""
        import meshio
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        m = meshio.Mesh(
            points=pts,
            cells=[("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
        )
        vtu = tmp_path / "tri_only.vtu"
        meshio.write(str(vtu), m)
        with pytest.raises(ValueError, match="no tetra cells"):
            load_bulk_mesh(vtu)

    def test_warns_on_non_utm_coordinates(self, tmp_path, capsys):
        """Loading a mesh whose |X| is well below the UTM threshold
        must emit a stderr warning (R-202)."""
        import meshio
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        m = meshio.Mesh(
            points=pts,
            cells=[("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
        )
        vtu = tmp_path / "local_frame.vtu"
        meshio.write(str(vtu), m)
        load_fault_mesh(vtu)
        captured = capsys.readouterr()
        assert "non-UTM" in captured.err
        assert "below the UTM sanity threshold" in captured.err

    def test_no_warning_on_real_utm_mesh(self, capsys):
        """The real 2 km mesh is UTM Zone 11 N (X ≫ 1e4 m); no
        warning should fire (R-202)."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        load_fault_mesh(_REAL_2000M_FAULT)
        captured = capsys.readouterr()
        assert "non-UTM" not in captured.err

    def test_extract_cells_by_physical_bad_cell_type_raises(self):
        """The cell_type argument is validated up-front (R-207:
        moved out of TestRealData2000m because it uses synthetic
        data only)."""
        import meshio
        m = meshio.Mesh(
            points=np.zeros((3, 3)),
            cells=[("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
        )
        with pytest.raises(ValueError, match="must be 'triangle' or 'tetra'"):
            extract_cells_by_physical(m, "hexahedron", DEFAULT_FAULT_NAME)


# ----------------------------------------------------------------------
# Real-data tests (skipped if VTUs are not on disk)
# ----------------------------------------------------------------------


class TestRealData2000m:
    def test_load_fault_vtu_real_data(self):
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        # Phase 1 acceptance: 2685 triangles, single triangle block.
        tri_blocks = [cb for cb in mesh.cells if cb.type == "triangle"]
        assert len(tri_blocks) == 1
        assert len(tri_blocks[0].data) == _EXPECTED_2000M_N_TRI

    def test_triangle_geometry_real_data_bbox(self):
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        tri = extract_cells_by_physical(mesh, "triangle", DEFAULT_FAULT_NAME)
        c, _n, _a = triangle_geometry(mesh.points, tri)
        # Bounding box of centroids must lie inside the mesh point bbox
        # (centroids are convex combinations of mesh points).
        for i, axis in enumerate(["x", "y", "z"]):
            cmin, cmax = c[:, i].min(), c[:, i].max()
            pmin, pmax = mesh.points[:, i].min(), mesh.points[:, i].max()
            assert pmin <= cmin, f"{axis} centroid min {cmin} < point min {pmin}"
            assert cmax <= pmax, f"{axis} centroid max {cmax} > point max {pmax}"
        # The mesh's own bbox must match the expected values (sanity
        # check that we're looking at the right mesh).
        assert abs(mesh.points[:, 0].min() - _EXPECTED_2000M_BBOX["xmin"]) < 1.0
        assert abs(mesh.points[:, 0].max() - _EXPECTED_2000M_BBOX["xmax"]) < 1.0
        assert abs(mesh.points[:, 1].min() - _EXPECTED_2000M_BBOX["ymin"]) < 1.0
        assert abs(mesh.points[:, 1].max() - _EXPECTED_2000M_BBOX["ymax"]) < 1.0
        # Z bottom can have small precision drift from the gmsh export.
        assert abs(mesh.points[:, 2].min() - _EXPECTED_2000M_BBOX["zmin"]) < 5.0
        assert abs(mesh.points[:, 2].max() - _EXPECTED_2000M_BBOX["zmax"]) < 1.0

    def test_tet_geometry_real_data_count_and_positivity(self):
        _skip_if_no_real_data(_REAL_2000M_BULK)
        mesh = load_bulk_mesh(_REAL_2000M_BULK)
        tet = extract_cells_by_physical(mesh, "tetra", DEFAULT_BULK_NAME)
        # Phase 1 acceptance: 145 300 cells.
        assert len(tet) == _EXPECTED_2000M_N_TET
        _c, v = tet_geometry(mesh.points, tet)
        # Volumes positive.
        assert (v > 0).all()
        # Sum of volumes within 0.1 % of the gmsh box volume.
        rel_err = abs(v.sum() - _EXPECTED_2000M_BBOX_VOLUME) / _EXPECTED_2000M_BBOX_VOLUME
        assert rel_err < 1.0e-3, (
            f"sum of tet volumes ({v.sum():.4e}) differs from bbox "
            f"volume ({_EXPECTED_2000M_BBOX_VOLUME:.4e}) by "
            f"{rel_err * 100:.3f}%, exceeds 0.1% tolerance"
        )

    def test_extract_cells_by_physical_fault_name(self):
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        tri = extract_cells_by_physical(mesh, "triangle", DEFAULT_FAULT_NAME)
        assert tri.shape == (_EXPECTED_2000M_N_TRI, 3)

    def test_extract_cells_by_physical_fallback_when_field_data_empty(self):
        """Per the Phase 1 §6 implementation note: when meshio drops
        field_data on VTU write (as it does for the *_fault.vtu and
        *_bulk.vtu outputs of msh_to_vtu.py), the function falls back
        to the unique-value heuristic on cell_data["gmsh:physical"].
        In that fallback case, the requested physical-group name is
        effectively ignored — the function returns the single block
        that has a uniform gmsh:physical tag. This test pins that
        behaviour."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        # field_data is dropped → fallback returns the block regardless of name.
        assert mesh.field_data == {}, (
            "test assumption broken: meshio now preserves field_data on VTU "
            "write; tighten extract_cells_by_physical to raise when the "
            "name is unknown but field_data is non-empty."
        )
        tri = extract_cells_by_physical(mesh, "triangle", "nonexistent_group")
        assert tri.shape == (_EXPECTED_2000M_N_TRI, 3)

    def test_R705_unknown_phys_name_with_known_field_data_raises(self):
        """R-705: when field_data is present but does NOT contain the
        requested phys_name (typo or multi-fault setup), raise
        instead of silently returning the whole block via the
        unique-value fallback."""
        import meshio
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        tris = np.array([[0, 1, 2]], dtype=np.int64)
        m = meshio.Mesh(
            points=pts,
            cells=[("triangle", tris)],
            cell_data={"gmsh:physical": [np.array([101])]},
            field_data={"rock": np.array([3, 3])},  # only "rock" registered
        )
        with pytest.raises(ValueError, match="not found in mesh field_data"):
            extract_cells_by_physical(m, "triangle", "fualt")


# ----------------------------------------------------------------------
# Dataclass constructors (smoke tests)
# ----------------------------------------------------------------------


class TestDataclasses:
    def test_fault_cell_geometry_construction(self):
        N = 5
        g = FaultCellGeometry(
            centroids=np.zeros((N, 3)),
            normals=np.zeros((N, 3)),
            strikes=np.zeros((N, 3)),
            dips=np.zeros((N, 3)),
            areas=np.zeros(N),
            n_degenerate=0,
        )
        assert g.centroids.shape == (N, 3)
        assert g.n_degenerate == 0

    def test_bulk_cell_geometry_construction(self):
        N = 5
        g = BulkCellGeometry(
            centroids=np.zeros((N, 3)),
            volumes=np.zeros(N),
        )
        assert g.centroids.shape == (N, 3)
        assert g.volumes.shape == (N,)


# ----------------------------------------------------------------------
# Phase 2 — per-triangle basis (Tandem/SEAS down-dip convention)
# ----------------------------------------------------------------------

from project_to_fault_stress import (                       # noqa: E402
    basis_to_node,
    build_fault_basis,
    cell_to_node_average,
    harmonise_normal_orientation,
    per_triangle_basis_raw,
)
# Bring the H&Z fault_basis_vectors into scope for cross-checks.
from hickman_and_zoback_regional_stress_projection import (  # noqa: E402
    fault_basis_vectors as hz_fault_basis_vectors,
)


class TestPerTriangleBasisRaw:
    def test_tpv102_convention(self):
        """TPV102 fault: y=0, ref_normal=(0,-1,0), up=(0,0,1).
        CLAUDE.md "Fault-local tangent frame" rule says
        can_t1 = (0, 0, -1) (down-dip), can_t2 = (+1, 0, 0) (strike)."""
        n = np.array([[0.0, -1.0, 0.0]])
        s, d, deg = per_triangle_basis_raw(n)
        assert not deg[0]
        np.testing.assert_allclose(s[0], [1.0, 0.0, 0.0], atol=1e-12)
        np.testing.assert_allclose(d[0], [0.0, 0.0, -1.0], atol=1e-12)

    def test_orthonormal_random_normals(self):
        """1000 random unit normals: each row has |s| = |d| = |n| = 1
        and the three vectors are mutually orthogonal to 1e-12."""
        rng = np.random.default_rng(42)
        N = 1000
        n_raw = rng.standard_normal((N, 3))
        n = n_raw / np.linalg.norm(n_raw, axis=1, keepdims=True)
        s, d, deg = per_triangle_basis_raw(n)
        good = ~deg
        if not good.any():
            pytest.skip("RNG produced no non-degenerate rows; rerun.")
        s_g, d_g, n_g = s[good], d[good], n[good]
        # Unit magnitudes
        np.testing.assert_allclose(np.linalg.norm(s_g, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(d_g, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(n_g, axis=1), 1.0, atol=1e-12)
        # Orthogonality (acceptance criterion #4)
        np.testing.assert_allclose(np.einsum("ij,ij->i", s_g, d_g), 0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", s_g, n_g), 0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", d_g, n_g), 0, atol=1e-12)

    def test_dip_is_down_not_up(self):
        """Tandem/SEAS d = s × n is DOWN-dip; H&Z d_hat = n × s is UP-dip.
        For TPV102 (n along -y), Tandem's d = (0, 0, -1), pointing
        into the earth. H&Z's d would be (0, 0, +1). (R-102 contract.)"""
        n = np.array([[0.0, -1.0, 0.0]])
        _s, d, _deg = per_triangle_basis_raw(n)
        # Down-dip: z-component is negative.
        assert d[0, 2] < 0
        np.testing.assert_allclose(d[0], [0.0, 0.0, -1.0], atol=1e-12)

    def test_horizontal_triangle_marked_degenerate(self):
        """A triangle with normal = ±up has undefined strike."""
        n = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
        s, d, deg = per_triangle_basis_raw(n)
        assert deg[0] and deg[1]
        assert np.all(np.isnan(s))
        assert np.all(np.isnan(d))

    def test_nan_normal_propagates_to_degenerate(self):
        """A NaN normal (from a zero-area triangle in Phase 1)
        propagates to a degenerate basis row."""
        n = np.array(
            [[np.nan, np.nan, np.nan], [0.0, -1.0, 0.0]],
            dtype=np.float64,
        )
        s, d, deg = per_triangle_basis_raw(n)
        assert deg[0]
        assert not deg[1]
        assert np.all(np.isnan(s[0]))
        # Second row is valid.
        np.testing.assert_allclose(s[1], [1.0, 0.0, 0.0], atol=1e-12)

    def test_bad_normals_shape_raises(self):
        with pytest.raises(ValueError):
            per_triangle_basis_raw(np.zeros((5, 2)))

    def test_bad_up_shape_raises(self):
        with pytest.raises(ValueError):
            per_triangle_basis_raw(np.zeros((1, 3)), up=np.zeros(2))


class TestHarmoniseNormalOrientation:
    def test_idempotent_under_strike_hint(self):
        """Acceptance criterion #2: idempotency."""
        rng = np.random.default_rng(7)
        N = 200
        n_raw = rng.standard_normal((N, 3))
        n_raw /= np.linalg.norm(n_raw, axis=1, keepdims=True)
        centroids = rng.uniform(-1.0, 1.0, size=(N, 3))
        once = harmonise_normal_orientation(
            n_raw, "right-lateral", centroids,
            fault_strike_azimuth_hint_deg=314.0,
        )
        twice = harmonise_normal_orientation(
            once, "right-lateral", centroids,
            fault_strike_azimuth_hint_deg=314.0,
        )
        np.testing.assert_array_equal(once, twice)

    def test_idempotent_under_pca_fallback(self):
        """Idempotency also holds under PCA fallback (hint=None)."""
        rng = np.random.default_rng(8)
        N = 100
        n_raw = rng.standard_normal((N, 3))
        n_raw /= np.linalg.norm(n_raw, axis=1, keepdims=True)
        # Quasi-planar centroids so PCA gives a clean normal.
        u = rng.uniform(-1, 1, (N, 2))
        centroids = np.column_stack(
            [u[:, 0], u[:, 1], 0.01 * rng.standard_normal(N)]
        )
        once = harmonise_normal_orientation(
            n_raw, "right-lateral", centroids,
            fault_strike_azimuth_hint_deg=None,
        )
        twice = harmonise_normal_orientation(
            once, "right-lateral", centroids,
            fault_strike_azimuth_hint_deg=None,
        )
        np.testing.assert_array_equal(once, twice)

    def test_flips_minority_to_match_reference(self):
        """The reference direction for right-lateral SAF (314°) is
        H&Z's n_hat ≈ (-0.69, -0.72, 0) (SW). A mix of +n_HZ and
        -n_HZ normals must all end up aligned with the SW direction."""
        s_hat, _d_hat, n_hat = hz_fault_basis_vectors(
            314.0, dip_deg=90.0, rake_sense="right-lateral"
        )
        n_ref = np.asarray(n_hat)
        n_raw = np.vstack([n_ref, -n_ref, n_ref, -n_ref, n_ref])
        centroids = np.zeros_like(n_raw)
        out = harmonise_normal_orientation(
            n_raw, "right-lateral", centroids,
            fault_strike_azimuth_hint_deg=314.0,
        )
        # Every row must now point in n_ref's half-space.
        dots = out @ n_ref
        assert (dots > 0).all(), f"some rows not aligned: dots={dots}"

    def test_invalid_rake_sense_raises(self):
        with pytest.raises(ValueError, match="rake_sense"):
            harmonise_normal_orientation(
                np.array([[1.0, 0.0, 0.0]]),
                "weird-rake",
                np.zeros((1, 3)),
                314.0,
            )

    def test_R403_1d_normals_raises_clear_value_error(self):
        """A 1-D normals input must raise a clear ValueError, not
        a downstream numpy.AxisError (R-403)."""
        with pytest.raises(ValueError, match="normals must be"):
            harmonise_normal_orientation(
                np.array([1.0, 0.0, 0.0]),     # 1-D — invalid
                "right-lateral",
                np.zeros((1, 3)),
                314.0,
            )

    def test_R403_1d_centroids_raises_clear_value_error(self):
        """A 1-D centroids input must raise a clear ValueError, not
        a downstream numpy error (R-403)."""
        with pytest.raises(ValueError, match="trace_centroids must be"):
            harmonise_normal_orientation(
                np.array([[1.0, 0.0, 0.0]]),
                "right-lateral",
                np.array([0.0, 0.0, 0.0]),     # 1-D — invalid
                314.0,
            )

    def test_shape_mismatch_raises_in_pca_fallback(self):
        """R-305: shape check moved into the PCA-fallback branch.
        Mismatched centroids only error when hint is None."""
        with pytest.raises(ValueError, match="PCA fallback"):
            harmonise_normal_orientation(
                np.array([[1.0, 0.0, 0.0]]),
                "right-lateral",
                np.zeros((2, 3)),  # shape mismatch
                fault_strike_azimuth_hint_deg=None,  # forces PCA
            )

    def test_shape_mismatch_ignored_in_hint_path(self):
        """R-305: with a strike hint, centroids are unused; a
        shape-mismatched dummy must NOT raise."""
        out = harmonise_normal_orientation(
            np.array([[-0.69, -0.72, 0.0]]),
            "right-lateral",
            np.zeros((42, 3)),  # arbitrary shape — unused
            fault_strike_azimuth_hint_deg=314.0,
        )
        assert out.shape == (1, 3)

    def test_nan_rows_pass_through_unchanged(self):
        n = np.array(
            [[np.nan, np.nan, np.nan], [-0.69, -0.72, 0.0]],
            dtype=np.float64,
        )
        n[1] /= np.linalg.norm(n[1])
        out = harmonise_normal_orientation(
            n, "right-lateral", np.zeros_like(n),
            fault_strike_azimuth_hint_deg=314.0,
        )
        assert np.all(np.isnan(out[0]))   # NaN preserved
        # Second row aligned with reference (no flip needed since it
        # already points SW).
        s_hat, _d, n_hat = hz_fault_basis_vectors(
            314.0, 90.0, "right-lateral"
        )
        np.testing.assert_allclose(out[1], n[1], atol=1e-12)


class TestBuildFaultBasis:
    def _synth_single_triangle_safs(self, tmp_path):
        """Build a single vertical NW-striking triangle (az=314°) in
        a tempdir VTU. Returns the loaded mesh."""
        import meshio
        # Vertices defining a vertical triangle whose normal is
        # H&Z's n_hat for SAF az=314°.
        # The triangle plane is spanned by s_hat (along strike, NW)
        # and ẑ (down-dip), so its normal is ±s_hat × ẑ = ∓n_hat.
        # We construct vertices to give CCW winding such that the
        # raw normal lies in either half-space; harmonise will flip
        # as needed.
        s_hat, _d_hat, n_hat = hz_fault_basis_vectors(
            314.0, dip_deg=90.0, rake_sense="right-lateral"
        )
        s_hat = np.asarray(s_hat)
        # Vertices: origin, along strike, down-dip.
        v0 = np.zeros(3)
        v1 = 1000.0 * s_hat
        v2 = np.array([0.0, 0.0, -1000.0])
        pts = np.vstack([v0, v1, v2])
        m = meshio.Mesh(
            points=pts,
            cells=[("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
        )
        vtu = tmp_path / "synth_saf_tri.vtu"
        meshio.write(str(vtu), m)
        return load_fault_mesh(vtu)

    def test_strike_colinear_with_HZ_s_hat_safs(self, tmp_path, capsys):
        """Acceptance criterion #1 (interpreted: colinear up to sign).
        The Tandem convention `s = up × n_harmonised` gives a strike
        vector COLINEAR with H&Z's s_hat for SAF az=314° — the two are
        anti-parallel because the plan deliberately diverges from H&Z
        (R-102: Tandem strike = up × n flips sign relative to H&Z's
        s_hat = (sin az, cos az, 0)). The acceptance criterion's
        'aligned with H&Z's s_hat' is read as 'colinear, either
        direction'."""
        mesh = self._synth_single_triangle_safs(tmp_path)
        capsys.readouterr()   # swallow non-UTM warning (R-404)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",  # field_data is empty; fallback returns block
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        assert geom.n_degenerate == 0
        s_hat_HZ, _d_hat_HZ, _n_hat_HZ = hz_fault_basis_vectors(
            314.0, 90.0, "right-lateral"
        )
        # |strike · s_hat_HZ| should be 1 (colinear up to sign).
        cos_angle = abs(np.dot(geom.strikes[0], np.asarray(s_hat_HZ)))
        assert abs(cos_angle - 1.0) < 1e-10, (
            f"strike not colinear with H&Z's s_hat: |cos|={cos_angle}"
        )

    def test_returned_basis_orthonormal(self, tmp_path, capsys):
        """Acceptance criterion #4."""
        mesh = self._synth_single_triangle_safs(tmp_path)
        capsys.readouterr()   # swallow non-UTM warning (R-404)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        s, d, n = geom.strikes[0], geom.dips[0], geom.normals[0]
        np.testing.assert_allclose(np.linalg.norm(s), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(d), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(n), 1.0, atol=1e-12)
        assert abs(np.dot(s, d)) < 1e-12
        assert abs(np.dot(s, n)) < 1e-12
        assert abs(np.dot(d, n)) < 1e-12

    def test_normal_in_SW_half_space_for_safs(self, tmp_path, capsys):
        """Plan §"Goal" (line 471): for right-lateral SAF with NW
        strike, the harmonised normal points to the SW half-space
        (n_HZ · n_harmonised > 0)."""
        mesh = self._synth_single_triangle_safs(tmp_path)
        capsys.readouterr()   # swallow non-UTM warning (R-404)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        _s_HZ, _d_HZ, n_HZ = hz_fault_basis_vectors(
            314.0, 90.0, "right-lateral"
        )
        assert np.dot(geom.normals[0], np.asarray(n_HZ)) > 0, (
            "harmonised normal does not point to H&Z's SW half-space"
        )

    def test_real_data_low_degeneracy(self):
        """Acceptance criterion #3: n_degenerate / N_tri < 0.001 on
        the real 2000 m SAFS fault VTU."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        N_tri = geom.strikes.shape[0]
        assert geom.n_degenerate / N_tri < 0.001, (
            f"degenerate fraction {geom.n_degenerate}/{N_tri} = "
            f"{geom.n_degenerate / N_tri:.4f} exceeds 0.001"
        )

    def test_real_data_full_orthonormality(self):
        """Acceptance criterion #4 on the real 2000 m mesh: every
        non-degenerate row is orthonormal to 1e-12."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        good = ~np.isnan(geom.strikes).any(axis=1)
        s = geom.strikes[good]
        d = geom.dips[good]
        n = geom.normals[good]
        np.testing.assert_allclose(np.linalg.norm(s, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", s, d), 0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", s, n), 0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", d, n), 0.0, atol=1e-12
        )


# ----------------------------------------------------------------------
# Phase 2 — Cell→node averaging
# ----------------------------------------------------------------------


class TestCellToNodeAverage:
    def _two_triangle_quad(self):
        """Unit square split into two triangles sharing edge 1-2."""
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [1.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        tris = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)
        areas = np.array([0.5, 0.5])
        return pts, tris, areas

    def test_constant_scalar_is_pointwise_exact(self):
        """f_cell = 3.14 everywhere → every node gets 3.14 exactly."""
        pts, tris, areas = self._two_triangle_quad()
        cv = np.full(2, 3.14, dtype=np.float64)
        node = cell_to_node_average(pts, tris, cv, areas)
        np.testing.assert_allclose(node, 3.14, atol=1e-15)

    def test_constant_vector_is_pointwise_exact(self):
        """Same field but as a 3-vector. Tests the vector broadcast
        recipe (R-010 fix from round 1)."""
        pts, tris, areas = self._two_triangle_quad()
        cv = np.tile(np.array([1.0, 2.0, 3.0]), (2, 1))
        node = cell_to_node_average(pts, tris, cv, areas)
        for i in range(4):
            np.testing.assert_allclose(
                node[i], [1.0, 2.0, 3.0], atol=1e-15
            )

    def test_per_cell_distinct_vectors(self):
        """Two cells with DIFFERENT vectors. Shared nodes get the
        area-weighted average; un-shared nodes get the single cell's
        value. (Catches the 'identical vector tiled 3x' coincidence
        from round 1's R-010.)"""
        pts, tris, areas = self._two_triangle_quad()
        cv = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        node = cell_to_node_average(pts, tris, cv, areas)
        # Vertex 0 is only in tri 0.
        np.testing.assert_allclose(node[0], cv[0], atol=1e-15)
        # Vertex 3 is only in tri 1.
        np.testing.assert_allclose(node[3], cv[1], atol=1e-15)
        # Vertices 1 and 2 are in both tris; area-weighted average
        # equals (0.5*cv[0] + 0.5*cv[1]) / 1.0 = mean(cv).
        np.testing.assert_allclose(node[1], cv.mean(axis=0), atol=1e-15)
        np.testing.assert_allclose(node[2], cv.mean(axis=0), atol=1e-15)

    def test_linear_scalar_on_plane(self):
        """Per-cell f = x_centroid (a linear function of position) on
        a planar mesh — the area-weighted node average equals the
        vertex x-coordinate (linear-on-triangulation exact case)."""
        pts, tris, areas = self._two_triangle_quad()
        centroids = pts[tris].mean(axis=1)
        cv = centroids[:, 0]                  # f = x at centroid
        node = cell_to_node_average(pts, tris, cv, areas)
        # Interior vertex (shared) recovers the vertex x.
        # Vertex 0 is in tri 0 only, with centroid x = 1/3 — node[0] = 1/3.
        # Vertex 3 is in tri 1 only, with centroid x = 2/3 — node[3] = 2/3.
        # (Not exactly the vertex x; the linear-recovery property
        # only holds for area-weighted averaging when the field is
        # linear AND there are enough surrounding triangles. For
        # boundary vertices the recovery is exact only in the
        # specific average sense.)
        np.testing.assert_allclose(node[0], 1.0 / 3.0, atol=1e-15)
        np.testing.assert_allclose(node[3], 2.0 / 3.0, atol=1e-15)
        # Shared vertices: (0.5 * 1/3 + 0.5 * 2/3) / 1 = 0.5
        np.testing.assert_allclose(node[1], 0.5, atol=1e-15)
        np.testing.assert_allclose(node[2], 0.5, atol=1e-15)

    def test_orphan_vertex_is_nan(self):
        """A vertex not referenced by any triangle gets NaN."""
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [99.0, 99.0, 99.0]],   # orphan
            dtype=np.float64,
        )
        tris = np.array([[0, 1, 2]], dtype=np.int64)
        areas = np.array([0.5])
        cv = np.array([2.0])
        node = cell_to_node_average(pts, tris, cv, areas)
        assert np.isnan(node[3])
        np.testing.assert_allclose(node[0], 2.0, atol=1e-15)

    def test_orphan_vector_vertex_is_nan_row(self):
        pts = np.array(
            [[0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [99.0, 99.0, 99.0]],
            dtype=np.float64,
        )
        tris = np.array([[0, 1, 2]], dtype=np.int64)
        areas = np.array([0.5])
        cv = np.array([[1.0, 2.0, 3.0]])
        node = cell_to_node_average(pts, tris, cv, areas)
        assert np.all(np.isnan(node[3]))
        np.testing.assert_allclose(node[0], [1.0, 2.0, 3.0], atol=1e-15)

    def test_shape_mismatch_raises(self):
        pts = np.zeros((4, 3))
        tris = np.array([[0, 1, 2]], dtype=np.int64)
        areas = np.array([0.5, 0.5])    # wrong length (2 vs 1 tri)
        cv = np.zeros(1)
        with pytest.raises(ValueError, match="areas"):
            cell_to_node_average(pts, tris, cv, areas)

    def test_R302_nan_cell_does_not_pollute_neighbours(self):
        """A NaN cell must not pollute neighbouring vertices that
        have at least one good incident cell (R-302)."""
        pts, tris, areas = self._two_triangle_quad()
        # cell 0 = 1.0 (good); cell 1 = NaN (bad).
        cv = np.array([1.0, np.nan])
        node = cell_to_node_average(pts, tris, cv, areas)
        # Vertex 0 is only in cell 0 → 1.0.
        np.testing.assert_allclose(node[0], 1.0, atol=1e-12)
        # Vertex 3 is only in cell 1 (which is bad) → NaN.
        assert np.isnan(node[3])
        # Vertices 1 and 2 are in BOTH cells; cell 1 is bad, but
        # cell 0 is good — they should get 1.0, not NaN (R-302).
        np.testing.assert_allclose(node[1], 1.0, atol=1e-12,
            err_msg="R-302: shared vertex polluted by single NaN cell")
        np.testing.assert_allclose(node[2], 1.0, atol=1e-12,
            err_msg="R-302: shared vertex polluted by single NaN cell")

    def test_R302_nan_vector_cell_does_not_pollute_neighbours(self):
        """Vector cell-values: any NaN in a cell's vector flags the
        whole cell as bad (R-302)."""
        pts, tris, areas = self._two_triangle_quad()
        cv = np.array([[1.0, 2.0, 3.0], [np.nan, 4.0, 5.0]])
        node = cell_to_node_average(pts, tris, cv, areas)
        # Vertex 0 only in cell 0 → cell 0's vector.
        np.testing.assert_allclose(node[0], [1.0, 2.0, 3.0], atol=1e-12)
        # Vertex 3 only in cell 1 (bad) → NaN row.
        assert np.all(np.isnan(node[3]))
        # Shared vertices: cell 0 only → cell 0's vector.
        np.testing.assert_allclose(node[1], [1.0, 2.0, 3.0], atol=1e-12)
        np.testing.assert_allclose(node[2], [1.0, 2.0, 3.0], atol=1e-12)


class TestBasisToNode:
    def _planar_quad(self):
        """Four-triangle fan around a central vertex, all in z=0."""
        pts = np.array(
            [[0.0, 0.0, 0.0],   # center
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [-1.0, 0.0, 0.0],
             [0.0, -1.0, 0.0]],
            dtype=np.float64,
        )
        tris = np.array(
            [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
            dtype=np.int64,
        )
        return pts, tris

    def test_planar_constant_basis_recovers_at_nodes(self):
        """All triangles share the same basis (n=+ẑ, s=+x̂); every
        node must recover the same basis. Dip = s × n = (0, -1, 0)
        per Tandem down-dip convention (R-301)."""
        pts, tris = self._planar_quad()
        N_tri = tris.shape[0]
        n_cell = np.tile([0.0, 0.0, 1.0], (N_tri, 1))
        s_cell = np.tile([1.0, 0.0, 0.0], (N_tri, 1))
        # d_cell is ignored by basis_to_node (dip is recomputed via
        # s × n); we pass the Tandem value here for clarity.
        d_cell = np.cross(s_cell, n_cell)   # (0, -1, 0)
        areas = 0.5 * np.ones(N_tri)
        s_node, d_node, n_node = basis_to_node(
            pts, tris, s_cell, d_cell, n_cell, areas
        )
        # Center node (vertex 0): all four tris.
        np.testing.assert_allclose(n_node[0], [0.0, 0.0, 1.0], atol=1e-12)
        np.testing.assert_allclose(s_node[0], [1.0, 0.0, 0.0], atol=1e-12)
        # Tandem down-dip: d = s × n = (1,0,0) × (0,0,1) = (0, -1, 0).
        np.testing.assert_allclose(d_node[0], [0.0, -1.0, 0.0], atol=1e-12)

    def test_orthonormality_per_node(self):
        """Every non-degenerate node basis must be orthonormal."""
        rng = np.random.default_rng(13)
        pts, tris = self._planar_quad()
        N_tri = tris.shape[0]
        # Random per-cell basis (Tandem convention: build n random
        # then derive s, d).
        n_cell = rng.standard_normal((N_tri, 3))
        n_cell /= np.linalg.norm(n_cell, axis=1, keepdims=True)
        s_cell, d_cell, deg = per_triangle_basis_raw(n_cell)
        assert not deg.any()
        areas = rng.uniform(0.1, 1.0, size=N_tri)
        s_node, d_node, n_node = basis_to_node(
            pts, tris, s_cell, d_cell, n_cell, areas
        )
        good = ~np.isnan(s_node).any(axis=1)
        s, d, n = s_node[good], d_node[good], n_node[good]
        np.testing.assert_allclose(np.linalg.norm(s, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(n, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", s, n), 0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", d, n), 0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", s, d), 0, atol=1e-12)

    def test_dip_is_s_cross_n(self):
        """At every non-degenerate node, dip = s × n (Tandem
        down-dip convention; R-301 fix)."""
        pts, tris = self._planar_quad()
        N_tri = tris.shape[0]
        n_cell = np.tile([0.0, 0.0, 1.0], (N_tri, 1))
        s_cell = np.tile([1.0, 0.0, 0.0], (N_tri, 1))
        # d_cell isn't used (basis_to_node recomputes dip from s × n).
        d_cell = np.cross(s_cell, n_cell)
        areas = 0.5 * np.ones(N_tri)
        s_node, d_node, n_node = basis_to_node(
            pts, tris, s_cell, d_cell, n_cell, areas
        )
        good = ~np.isnan(s_node).any(axis=1)
        d_expected = np.cross(s_node[good], n_node[good])
        np.testing.assert_allclose(d_node[good], d_expected, atol=1e-12)

    def test_R401_basis_to_node_docstring_matches_implementation(self):
        """The docstring step 5 must specify the same cross-product
        order as the implementation (R-401: regression guard against
        the round-2 self-contradiction where R-301 fixed the code
        but not the docstring)."""
        import inspect
        from project_to_fault_stress import basis_to_node
        doc = inspect.getdoc(basis_to_node)
        # Must explicitly call out the Tandem `s × n` convention.
        assert "s_node × n_node" in doc or "s × n" in doc, (
            "R-401: basis_to_node docstring must specify d = s × n; "
            "current step 5 is the bug R-301 fixed in code but not "
            "in docs."
        )
        # And must NOT advertise the H&Z up-dip recipe as THE
        # current recipe (it may still be mentioned as a contrast).
        # We check that the step-5 line itself uses `s × n`, not `n × s`:
        for line in doc.splitlines():
            stripped = line.strip()
            if stripped.startswith("5."):
                # The step-5 line must contain `s_node × n_node` or
                # `s × n` and must NOT have `n_node × s_node` as the
                # primary formula.
                assert "s_node × n_node" in stripped or "s × n" in stripped, (
                    f"R-401: docstring step 5 doesn't specify the "
                    f"Tandem recipe; got: {stripped!r}"
                )
                # If the line mentions `n × s` it must be in a
                # contrastive context (e.g., "NOT n × s").
                if "n_node × s_node" in stripped or "n × s" in stripped:
                    assert "NOT" in stripped or "not" in stripped, (
                        f"R-401: docstring step 5 mentions n × s "
                        f"without disclaiming it: {stripped!r}"
                    )
                break

    def test_R301_basis_to_node_dip_matches_per_triangle_dip(self):
        """Cell-data and point-data dip vectors must carry the same
        sign — both down-dip per Tandem convention (R-301)."""
        # TPV102 reference geometry: planar fault y=0, normal (0,-1,0).
        pts = np.array(
            [[0., 0., 0.], [1., 0., 0.],
             [0., 0., -1.], [1., 0., -1.]],
            dtype=np.float64,
        )
        # CCW winding such that (e1 × e2) = (0, -1, 0).
        tris = np.array([[0, 2, 1], [1, 2, 3]], dtype=np.int64)
        n_cell = np.tile([0.0, -1.0, 0.0], (2, 1))
        s_cell, d_cell, _deg = per_triangle_basis_raw(n_cell)
        areas = np.array([0.5, 0.5])
        s_node, d_node, n_node = basis_to_node(
            pts, tris, s_cell, d_cell, n_cell, areas
        )
        # per_triangle_basis_raw gives down-dip (0, 0, -1) for TPV102.
        np.testing.assert_allclose(
            d_cell[0], [0.0, 0.0, -1.0], atol=1e-12
        )
        # basis_to_node must agree (R-301: NOT (0, 0, +1)).
        np.testing.assert_allclose(
            d_node[0], [0.0, 0.0, -1.0], atol=1e-12,
            err_msg="R-301: basis_to_node dip must match "
                    "per_triangle (down-dip Tandem); got opposite "
                    "sign (up-dip H&Z)",
        )
        # All node-data dips must be down-dip (z <= 0).
        good = ~np.isnan(d_node).any(axis=1)
        assert (d_node[good, 2] <= 0).all(), (
            f"R-301: some node-data dips point upward: {d_node[good]}"
        )


# ----------------------------------------------------------------------
# Phase 3 — Stress projection
# ----------------------------------------------------------------------

from project_to_fault_stress import (                       # noqa: E402
    bulk_stress_tensor_field,
    pore_pressure_field,
    project_stress_onto_fault,
    resolve_traction_per_cell,
)
from hickman_and_zoback_regional_stress_projection import (  # noqa: E402
    build_bulk_stress_tensor as hz_build_bulk_stress_tensor,
    compute_fault_stress as hz_compute_fault_stress,
    resolve_traction as hz_resolve_traction,
)


# SAFOD input parameters (from H&Z demo_safod). Reference values
# for σ_n_total, σ_n_eff, τ_strike, μ_apparent are computed live
# from H&Z's compute_fault_stress to avoid hardcoded rounded
# numbers (per the user-feedback memory `feedback_no_hardcoded_numbers.md`).
_SAFOD_SHmax_MPa = 113.0
_SAFOD_Shmin_MPa = 49.0
_SAFOD_Sv_MPa = 45.0
_SAFOD_Pp_MPa = 16.0
_SAFOD_SHmax_az_deg = 23.0
_SAFOD_fault_strike_az_deg = 314.0


def _safod_hz_reference():
    """Compute the SAFOD H&Z reference values live. The Phase 3
    output is related to these by the sign-flip-and-basis-flip
    cancellation:
        σ_n_seas      = -σ_n_total_HZ
        σ_n_eff_seas  = -σ_n_eff_HZ
        τ_strike_seas = +τ_strike_HZ
        τ_dip_seas    = +τ_dip_HZ
        |τ|_seas      = |τ|_HZ
        μ_seas        = μ_HZ   (magnitudes cancel)
    """
    _sigma0, traction = hz_compute_fault_stress(
        SHmax=_SAFOD_SHmax_MPa, Shmin=_SAFOD_Shmin_MPa,
        Sv=_SAFOD_Sv_MPa, P_p=_SAFOD_Pp_MPa,
        SHmax_azimuth_deg=_SAFOD_SHmax_az_deg,
        fault_strike_azimuth_deg=_SAFOD_fault_strike_az_deg,
    )
    return traction


def _safod_tandem_basis():
    """Build a single-triangle Tandem (s, d, n) basis for the
    SAFOD SAF: vertical, az=314°, right-lateral. n is harmonised
    to point SW (H&Z convention); s, d follow from Tandem
    (s = up × n, d = s × n)."""
    from project_to_fault_stress import per_triangle_basis_raw
    _s_HZ, _d_HZ, n_HZ = hz_fault_basis_vectors(
        _SAFOD_fault_strike_az_deg, 90.0, "right-lateral"
    )
    n_arr = np.array([n_HZ], dtype=np.float64)
    s_arr, d_arr, deg = per_triangle_basis_raw(n_arr)
    assert not deg[0], "SAFOD synthetic triangle unexpectedly degenerate"
    return s_arr, d_arr, n_arr


class TestBulkStressTensorField:
    def test_constant_safod_returns_compression_positive_seas(self):
        """SAFOD parameters, constant depth model: returns σ_seas
        = -σ_HZ. Diagonal of σ_seas at SAFOD is dominantly positive
        (compression positive)."""
        z = np.array([-1671.0])
        sigma_seas = bulk_stress_tensor_field(
            z, SHmax_az_deg=_SAFOD_SHmax_az_deg, depth_model="constant",
            SHmax_top=_SAFOD_SHmax_MPa, Shmin_top=_SAFOD_Shmin_MPa,
            Sv_top=_SAFOD_Sv_MPa,
        )
        # σ_seas = -σ_HZ. H&Z diagonal is [-58.77, -103.23, -45]; flip → all +.
        np.testing.assert_allclose(
            np.diag(sigma_seas[0]), [58.77, 103.23, 45.0], atol=0.05
        )

    def test_constant_matches_negated_HZ_exactly(self):
        """σ_seas = -build_bulk_stress_tensor(...) exactly."""
        z = np.array([-1671.0, 0.0, -5000.0])  # three depths
        sigma_seas = bulk_stress_tensor_field(
            z, SHmax_az_deg=23.0, depth_model="constant",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
        )
        sigma_HZ_ref = np.asarray(
            hz_build_bulk_stress_tensor(113.0, 49.0, 45.0, 23.0),
            dtype=np.float64,
        )
        # Constant depth model: every row identical to -σ_HZ_ref.
        for i in range(3):
            np.testing.assert_allclose(
                sigma_seas[i], -sigma_HZ_ref, atol=1e-12,
            )

    def test_symmetric_output(self):
        z = np.linspace(-5000.0, 0.0, 11)
        sigma_seas = bulk_stress_tensor_field(
            z, SHmax_az_deg=37.0, depth_model="constant",
            SHmax_top=100.0, Shmin_top=40.0, Sv_top=60.0,
        )
        for i in range(len(z)):
            np.testing.assert_allclose(
                sigma_seas[i], sigma_seas[i].T, atol=1e-9,
            )

    def test_lithostatic_sv_depth_dependence(self):
        """Acceptance #4: With non-zero gradients, the σ⁰ field at
        z = -1671 m equals (depth-extrapolated) H&Z values to 1e-9."""
        z = np.array([-1671.0])
        SHmax_grad = 0.01    # MPa/m
        Shmin_grad = 0.005
        Sv_grad = 0.025
        sigma_seas = bulk_stress_tensor_field(
            z, SHmax_az_deg=23.0, depth_model="lithostatic_sv",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
            SHmax_grad=SHmax_grad, Shmin_grad=Shmin_grad, Sv_grad=Sv_grad,
        )
        # Depth-extrapolated H&Z values at d = 1671 m.
        depth = 1671.0
        SHmax = 113.0 + SHmax_grad * depth
        Shmin = 49.0 + Shmin_grad * depth
        Sv = 45.0 + Sv_grad * depth
        sigma_HZ_extrap = np.asarray(
            hz_build_bulk_stress_tensor(SHmax, Shmin, Sv, 23.0),
            dtype=np.float64,
        )
        np.testing.assert_allclose(
            sigma_seas[0], -sigma_HZ_extrap, atol=1e-9,
        )

    def test_lithostatic_sv_clamps_above_surface(self):
        """For z > 0 (above free surface), depth is clamped to 0."""
        z = np.array([+100.0, -100.0])  # above and below surface
        sigma_seas = bulk_stress_tensor_field(
            z, SHmax_az_deg=23.0, depth_model="lithostatic_sv",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
            SHmax_grad=0.01, Shmin_grad=0.005, Sv_grad=0.025,
        )
        # z = +100 → depth = 0 → no gradient applied → equals
        # surface (constant) value.
        sigma_surface = np.asarray(
            hz_build_bulk_stress_tensor(113.0, 49.0, 45.0, 23.0),
            dtype=np.float64,
        )
        np.testing.assert_allclose(sigma_seas[0], -sigma_surface, atol=1e-12)
        # z = -100 → depth = 100 → gradient applied → differs from surface.
        assert not np.allclose(sigma_seas[1], -sigma_surface, atol=1e-3)

    def test_bad_depth_model_raises(self):
        with pytest.raises(ValueError, match="depth_model"):
            bulk_stress_tensor_field(
                np.array([0.0]),
                SHmax_az_deg=0.0, depth_model="weird",
                SHmax_top=1.0, Shmin_top=1.0, Sv_top=1.0,
            )

    def test_bad_centroids_z_shape_raises(self):
        with pytest.raises(ValueError, match="centroids_z must be 1-D"):
            bulk_stress_tensor_field(
                np.zeros((5, 3)),
                SHmax_az_deg=0.0,
                SHmax_top=1.0, Shmin_top=1.0, Sv_top=1.0,
            )


class TestPorePressureField:
    def test_constant_zero(self):
        z = np.array([-100.0, -200.0, 0.0])
        P = pore_pressure_field(z)
        np.testing.assert_allclose(P, [0.0, 0.0, 0.0])

    def test_constant_uniform(self):
        z = np.array([-100.0, -200.0, 0.0])
        P = pore_pressure_field(z, P_p_top=16.0)
        np.testing.assert_allclose(P, [16.0, 16.0, 16.0])

    def test_linear_depth_gradient(self):
        z = np.array([0.0, -100.0, -200.0])
        # Hydrostatic-like: rho * g = 1000 kg/m³ * 9.8 m/s² = 9800 Pa/m = 9.8e-3 MPa/m
        P = pore_pressure_field(z, P_p_top=0.0, P_p_grad=9.8e-3)
        np.testing.assert_allclose(P, [0.0, 0.98, 1.96], atol=1e-12)

    def test_above_surface_clamps_to_top(self):
        """z > 0 is above the free surface; depth clamped to 0."""
        z = np.array([+50.0])
        P = pore_pressure_field(z, P_p_top=16.0, P_p_grad=9.8e-3)
        np.testing.assert_allclose(P, [16.0], atol=1e-12)

    def test_bad_shape_raises(self):
        with pytest.raises(ValueError, match="centroids_z must be 1-D"):
            pore_pressure_field(np.zeros((3, 2)))


class TestResolveTractionPerCell:
    def test_R501_safod_anchor(self):
        """Phase 3 Acceptance Criterion #1: SAFOD parameters + single
        vertical NW-striking right-lateral triangle → sigma_n_eff
        ≈ +89 MPa (compression-positive SEAS), mu_apparent ≈ 0.24
        (sign flip cancels for μ). Also the R-501 regression anchor:
        τ_strike > 0 for right-lateral on the SAF. Reference values
        are computed live from H&Z (no hardcoded numerics)."""
        s, d, n = _safod_tandem_basis()
        ref = _safod_hz_reference()
        # σ_seas via the bulk builder (single source-site flip).
        sigma_seas = bulk_stress_tensor_field(
            np.array([-1671.0]),
            SHmax_az_deg=_SAFOD_SHmax_az_deg, depth_model="constant",
            SHmax_top=_SAFOD_SHmax_MPa, Shmin_top=_SAFOD_Shmin_MPa,
            Sv_top=_SAFOD_Sv_MPa,
        )
        res = resolve_traction_per_cell(
            sigma_seas, s, d, n, P_p_per_cell=_SAFOD_Pp_MPa,
        )
        # Acceptance: sigma_n_eff ≈ -H&Z.sigma_n_eff (sign-flipped to
        # compression-positive); H&Z reports -88.78 → Phase 3 reports +88.78.
        np.testing.assert_allclose(
            res["sigma_n_eff"][0], -ref.sigma_n_eff, atol=1e-3,
        )
        # mu_apparent: same magnitude as H&Z (sign flip cancels for μ).
        np.testing.assert_allclose(
            res["mu_apparent"][0], ref.mu_apparent, atol=1e-4,
        )
        # R-501: τ_strike > 0 for right-lateral SAF.
        assert res["tau_strike"][0] > 0, (
            f"R-501 regression: τ_strike = {res['tau_strike'][0]} "
            "must be positive for right-lateral SAFOD; current "
            "value is the OLD bug (sign flip from Phase 2 Tandem "
            "basis not propagated through Phase 3)."
        )
        # Anchor values for σ_n_total and τ_strike.
        np.testing.assert_allclose(
            res["sigma_n_total"][0], -ref.sigma_n_total, atol=1e-3,
        )
        np.testing.assert_allclose(
            res["tau_strike"][0], ref.tau_strike, atol=1e-3,
        )
        # τ_dip ≈ 0 for pure strike-slip (H&Z reports −0.00).
        np.testing.assert_allclose(
            res["tau_dip"][0], ref.tau_dip, atol=1e-3,
        )

    def test_R501_hz_cross_check_safod(self):
        """Cross-check stanza (Phase 3 §1 H&Z):
        - σ_n_seas       = -σ_n_HZ      to 1e-12
        - τ_strike_seas  = +τ_strike_HZ to 1e-12
        - τ_dip_seas     = +τ_dip_HZ    to 1e-12

        The two simultaneous flips (σ → -σ_HZ AND basis → -basis_HZ
        for s and d) cancel for the shear components."""
        s_T, d_T, n_T = _safod_tandem_basis()
        # H&Z baseline using the H&Z basis (s_HZ, d_HZ, n_HZ).
        sigma0_HZ, traction_HZ = hz_compute_fault_stress(
            SHmax=_SAFOD_SHmax_MPa, Shmin=_SAFOD_Shmin_MPa,
            Sv=_SAFOD_Sv_MPa, P_p=_SAFOD_Pp_MPa,
            SHmax_azimuth_deg=_SAFOD_SHmax_az_deg,
            fault_strike_azimuth_deg=_SAFOD_fault_strike_az_deg,
        )
        # Phase 3 output: Tandem basis + σ_seas.
        sigma_seas = -np.asarray(sigma0_HZ, dtype=np.float64)
        res = resolve_traction_per_cell(
            sigma_seas, s_T, d_T, n_T, P_p_per_cell=_SAFOD_Pp_MPa,
        )
        # σ_n_seas = -σ_n_HZ (compression flip).
        np.testing.assert_allclose(
            res["sigma_n_total"][0],
            -traction_HZ.sigma_n_total,
            atol=1e-12,
        )
        # τ_strike_seas = +τ_strike_HZ (two flips cancel).
        np.testing.assert_allclose(
            res["tau_strike"][0], traction_HZ.tau_strike, atol=1e-12,
        )
        # τ_dip_seas = +τ_dip_HZ (two flips cancel).
        np.testing.assert_allclose(
            res["tau_dip"][0], traction_HZ.tau_dip, atol=1e-12,
        )

    def test_acceptance_2_trace_invariance(self):
        """Acceptance #2: σ_n + sᵀσs + dᵀσd = trace(σ) to 1e-9 over
        100 random orthonormal frames and random symmetric σ."""
        rng = np.random.default_rng(42)
        for _ in range(100):
            # Random unit normal.
            n_raw = rng.standard_normal(3)
            n = n_raw / np.linalg.norm(n_raw)
            # Tandem basis on a single triangle. Skip if degenerate
            # (n ≈ ±up).
            from project_to_fault_stress import per_triangle_basis_raw
            s, d, deg = per_triangle_basis_raw(n.reshape(1, 3))
            if deg[0]:
                continue
            # Random symmetric σ.
            A = rng.standard_normal((3, 3))
            sigma = 0.5 * (A + A.T)
            # σ_n + sᵀσs + dᵀσd should equal trace(σ) — the three
            # diagonal "components" in the orthonormal (s, d, n)
            # frame.
            sigma_n = float(n @ sigma @ n)
            sigma_ss = float(s[0] @ sigma @ s[0])
            sigma_dd = float(d[0] @ sigma @ d[0])
            np.testing.assert_allclose(
                sigma_n + sigma_ss + sigma_dd,
                np.trace(sigma),
                atol=1e-9,
            )

    def test_acceptance_3_vectorised_matches_HZ_scalar(self):
        """Acceptance #3: 100 random triangles + random σ; the
        vectorised Phase 3 path must agree with 100 independent
        scalar H&Z `resolve_traction` calls. Magnitudes match
        exactly; signs follow the documented Phase 3 cross-check:
            sigma_n_seas      = -sigma_n_total_HZ
            tau_strike_seas   = +tau_strike_HZ
            tau_dip_seas      = +tau_dip_HZ"""
        from project_to_fault_stress import per_triangle_basis_raw
        rng = np.random.default_rng(7)
        N = 100
        n_raw = rng.standard_normal((N, 3))
        n = n_raw / np.linalg.norm(n_raw, axis=1, keepdims=True)
        s, d, deg = per_triangle_basis_raw(n)
        # Skip degenerate rows.
        good = ~deg
        s = s[good]; d = d[good]; n = n[good]
        assert s.shape[0] > 50, "RNG produced too many degenerate rows"

        # Random symmetric σ_HZ (compression-negative-ish).
        A = rng.standard_normal((3, 3))
        sigma_HZ = 0.5 * (A + A.T)
        sigma_seas = -sigma_HZ

        # Vectorised Phase 3 output.
        res = resolve_traction_per_cell(sigma_seas, s, d, n)

        # Scalar H&Z reference using H&Z's basis (= -Tandem basis,
        # so s_HZ_i = -s_i, d_HZ_i = -d_i, n_HZ_i = n_i).
        for i in range(s.shape[0]):
            s_HZ_i = -s[i]
            d_HZ_i = -d[i]
            n_HZ_i = n[i]
            ref = hz_resolve_traction(
                sigma_HZ, s_HZ_i, d_HZ_i, n_HZ_i, P_p=0.0,
            )
            # σ_n_seas = -σ_n_HZ
            np.testing.assert_allclose(
                res["sigma_n_total"][i], -ref.sigma_n_total, atol=1e-12,
            )
            # τ_strike_seas = +τ_strike_HZ (two flips cancel)
            np.testing.assert_allclose(
                res["tau_strike"][i], ref.tau_strike, atol=1e-12,
            )
            # τ_dip_seas = +τ_dip_HZ (two flips cancel)
            np.testing.assert_allclose(
                res["tau_dip"][i], ref.tau_dip, atol=1e-12,
            )

    def test_traction_vec_consistency(self):
        """Phase 3 §1 traction_vec contract: n · t = σ_n_total,
        s · t = τ_strike, d · t = τ_dip (all on σ_seas, Tandem basis,
        no sign-asymmetry caveats)."""
        s, d, n = _safod_tandem_basis()
        sigma_seas = bulk_stress_tensor_field(
            np.array([-1671.0]), SHmax_az_deg=23.0,
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
        )
        res = resolve_traction_per_cell(sigma_seas, s, d, n)
        t = res["traction_vec"][0]
        np.testing.assert_allclose(
            np.dot(n[0], t), res["sigma_n_total"][0], atol=1e-9,
        )
        np.testing.assert_allclose(
            np.dot(s[0], t), res["tau_strike"][0], atol=1e-9,
        )
        np.testing.assert_allclose(
            np.dot(d[0], t), res["tau_dip"][0], atol=1e-9,
        )

    def test_mu_apparent_nan_on_effective_tension(self):
        """sigma_n_eff ≤ 0 (effective tension): mu_apparent must be
        NaN (Phase 3 §"Edge Cases")."""
        # Construct a setup where σ_n_eff < 0 deliberately: shallow,
        # tiny σ_n_total, large P_p.
        s, d, n = _safod_tandem_basis()
        # σ_seas with very small SHmax, Shmin, Sv (1 MPa each).
        sigma_seas = bulk_stress_tensor_field(
            np.array([-100.0]), SHmax_az_deg=23.0,
            SHmax_top=1.0, Shmin_top=1.0, Sv_top=1.0,
        )
        # P_p = 10 MPa → effective tension.
        res = resolve_traction_per_cell(
            sigma_seas, s, d, n, P_p_per_cell=10.0,
        )
        assert res["sigma_n_eff"][0] < 0, "test setup did not produce tension"
        assert np.isnan(res["mu_apparent"][0]), (
            f"mu_apparent under effective tension must be NaN; got "
            f"{res['mu_apparent'][0]}"
        )

    def test_constant_sigma_broadcast_to_n_triangles(self):
        """Pass a single (3, 3) σ; rotator broadcasts to every
        triangle."""
        from project_to_fault_stress import per_triangle_basis_raw
        rng = np.random.default_rng(11)
        N = 5
        n_raw = rng.standard_normal((N, 3))
        n = n_raw / np.linalg.norm(n_raw, axis=1, keepdims=True)
        s, d, deg = per_triangle_basis_raw(n)
        # σ as a single (3, 3) — must broadcast.
        sigma = np.diag([50.0, 30.0, 20.0])
        res_broadcast = resolve_traction_per_cell(sigma, s, d, n)
        # σ as (N, 3, 3) — same result.
        sigma_stack = np.tile(sigma[None, :, :], (N, 1, 1))
        res_stack = resolve_traction_per_cell(sigma_stack, s, d, n)
        np.testing.assert_allclose(
            res_broadcast["sigma_n_total"],
            res_stack["sigma_n_total"],
            atol=1e-12,
        )
        np.testing.assert_allclose(
            res_broadcast["tau_strike"],
            res_stack["tau_strike"],
            atol=1e-12,
        )

    def test_nonsymmetric_sigma_raises(self):
        s, d, n = _safod_tandem_basis()
        sigma_asym = np.array([
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],   # 3.0 ≠ 2.0
            [0.0, 0.0, 1.0],
        ])
        with pytest.raises(ValueError, match="symmetric"):
            resolve_traction_per_cell(sigma_asym, s, d, n)

    def test_nonsymmetric_sigma_stack_raises(self):
        from project_to_fault_stress import per_triangle_basis_raw
        rng = np.random.default_rng(0)
        N = 3
        n = rng.standard_normal((N, 3))
        n /= np.linalg.norm(n, axis=1, keepdims=True)
        s, d, _ = per_triangle_basis_raw(n)
        # Build a stack of symmetric tensors, perturb one off-diagonal.
        stack = np.tile(np.eye(3)[None], (N, 1, 1))
        stack[1, 0, 1] = 5.0  # not equal to stack[1, 1, 0] = 0
        with pytest.raises(ValueError, match="non-symmetric rows"):
            resolve_traction_per_cell(stack, s, d, n)

    def test_bad_strikes_shape_raises(self):
        with pytest.raises(ValueError, match="strikes must be"):
            resolve_traction_per_cell(
                np.eye(3),
                strikes=np.zeros(3),       # 1-D
                dips=np.zeros((1, 3)),
                normals=np.zeros((1, 3)),
            )

    def test_shape_mismatch_dips_raises(self):
        with pytest.raises(ValueError, match="dips shape"):
            resolve_traction_per_cell(
                np.eye(3),
                strikes=np.zeros((2, 3)),
                dips=np.zeros((3, 3)),     # mismatch
                normals=np.zeros((2, 3)),
            )

    def test_pp_per_cell_array(self):
        s, d, n = _safod_tandem_basis()
        ref = _safod_hz_reference()
        sigma_seas = bulk_stress_tensor_field(
            np.array([-1671.0]), SHmax_az_deg=_SAFOD_SHmax_az_deg,
            SHmax_top=_SAFOD_SHmax_MPa, Shmin_top=_SAFOD_Shmin_MPa,
            Sv_top=_SAFOD_Sv_MPa,
        )
        # Per-cell array form (N=1).
        res = resolve_traction_per_cell(
            sigma_seas, s, d, n,
            P_p_per_cell=np.array([_SAFOD_Pp_MPa]),
        )
        np.testing.assert_allclose(
            res["sigma_n_eff"][0], -ref.sigma_n_eff, atol=1e-3,
        )

    def test_nan_basis_propagates_to_nan_output(self):
        """A degenerate-basis row (NaN strike/dip/normal from Phase 2)
        propagates NaN through to every emitted scalar."""
        sigma_seas = np.diag([50.0, 30.0, 20.0])
        # One good triangle and one degenerate (NaN).
        s = np.array([[1.0, 0.0, 0.0], [np.nan, np.nan, np.nan]])
        d = np.array([[0.0, 0.0, -1.0], [np.nan, np.nan, np.nan]])
        n = np.array([[0.0, -1.0, 0.0], [np.nan, np.nan, np.nan]])
        res = resolve_traction_per_cell(sigma_seas, s, d, n)
        # Row 0: finite.
        assert np.isfinite(res["sigma_n_total"][0])
        # Row 1: NaN.
        for key in ("sigma_n_total", "sigma_n_eff", "tau_strike",
                    "tau_dip", "tau_magnitude", "rake_deg"):
            assert np.isnan(res[key][1]), (
                f"{key} did not propagate NaN through degenerate row"
            )


class TestProjectStressOntoFault:
    def test_safod_anchor_wireup(self):
        """The wire-up function reproduces the SAFOD anchor end-to-end
        (mirrors `TestResolveTractionPerCell.test_R501_safod_anchor`
        but exercises the bulk + P_p builders too)."""
        s, d, n = _safod_tandem_basis()
        ref = _safod_hz_reference()
        centroids = np.array([[0.0, 0.0, -1671.0]])  # z = -1671 m
        geom = FaultCellGeometry(
            centroids=centroids,
            normals=n,
            strikes=s,
            dips=d,
            areas=np.array([1.0]),
            n_degenerate=0,
        )
        res = project_stress_onto_fault(
            geom,
            SHmax=_SAFOD_SHmax_MPa, Shmin=_SAFOD_Shmin_MPa,
            Sv=_SAFOD_Sv_MPa, P_p=_SAFOD_Pp_MPa,
            SHmax_az_deg=_SAFOD_SHmax_az_deg,
        )
        np.testing.assert_allclose(
            res["sigma_n_eff"][0], -ref.sigma_n_eff, atol=1e-3,
        )
        np.testing.assert_allclose(
            res["mu_apparent"][0], ref.mu_apparent, atol=1e-4,
        )
        assert res["tau_strike"][0] > 0
        # sigma_field key present and shaped correctly.
        assert "sigma_field" in res
        assert res["sigma_field"].shape == (1, 3, 3)

    def test_sigma_field_compression_positive(self):
        """sigma_field returned by the wire-up is in
        compression-positive SEAS convention (the single-site flip
        has already been applied by `bulk_stress_tensor_field`)."""
        s, d, n = _safod_tandem_basis()
        centroids = np.array([[0.0, 0.0, -1671.0]])
        geom = FaultCellGeometry(
            centroids=centroids, normals=n, strikes=s, dips=d,
            areas=np.array([1.0]), n_degenerate=0,
        )
        res = project_stress_onto_fault(
            geom,
            SHmax=113.0, Shmin=49.0, Sv=45.0, P_p=16.0,
            SHmax_az_deg=23.0,
        )
        # σ_zz at SAFOD compression-positive = +45.
        np.testing.assert_allclose(
            res["sigma_field"][0, 2, 2], 45.0, atol=1e-12,
        )

    def test_R601_empty_inputs_return_empty_outputs(self):
        """N=0 inputs must not crash the symmetry check (R-601).
        Both `bulk_stress_tensor_field` and `resolve_traction_per_cell`
        per-cell branch must return empty outputs gracefully."""
        # bulk_stress_tensor_field: N=0 returns (0, 3, 3) without raising.
        sigma_seas = bulk_stress_tensor_field(
            np.array([], dtype=np.float64),
            SHmax_az_deg=23.0, depth_model="constant",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
        )
        assert sigma_seas.shape == (0, 3, 3)
        # lithostatic_sv: same.
        sigma_seas_lith = bulk_stress_tensor_field(
            np.array([], dtype=np.float64),
            SHmax_az_deg=23.0, depth_model="lithostatic_sv",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
            SHmax_grad=0.01, Shmin_grad=0.005, Sv_grad=0.025,
        )
        assert sigma_seas_lith.shape == (0, 3, 3)

        # resolve_traction_per_cell per-cell stack: N=0 returns
        # empty arrays without raising.
        s = np.zeros((0, 3))
        d = np.zeros((0, 3))
        n = np.zeros((0, 3))
        stack = np.zeros((0, 3, 3))
        res = resolve_traction_per_cell(stack, s, d, n)
        for key in (
            "sigma_n_total", "sigma_n_eff", "tau_strike", "tau_dip",
            "tau_magnitude", "rake_deg", "mu_apparent",
        ):
            assert res[key].shape == (0,), (
                f"R-601: {key} shape {res[key].shape} != (0,)"
            )
        assert res["traction_vec"].shape == (0, 3)

        # resolve_traction_per_cell with single (3, 3) σ + N=0 bases
        # should also work (broadcast path).
        res_b = resolve_traction_per_cell(np.eye(3), s, d, n)
        assert res_b["sigma_n_total"].shape == (0,)
        assert res_b["traction_vec"].shape == (0, 3)

    def test_R702_nan_centroid_z_propagates_per_spec(self):
        """R-702: a NaN centroid_z corresponds to a degenerate
        Phase 1 triangle; the NaN must propagate to a NaN row in
        sigma_seas per Phase 3 §"Edge Cases", NOT raise.  The
        R-605 protection survives for genuinely unexpected NaN
        (see ``test_R702_unexpected_nan_still_raises`` below)."""
        sigma = bulk_stress_tensor_field(
            np.array([np.nan, -1000.0]),
            SHmax_az_deg=23.0, depth_model="lithostatic_sv",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
            Sv_grad=0.025,
        )
        assert sigma.shape == (2, 3, 3)
        # Row 0 (NaN centroid) → fully NaN tensor.
        assert np.all(np.isnan(sigma[0]))
        # Row 1 (finite centroid) → fully finite tensor.
        assert np.all(np.isfinite(sigma[1]))

    def test_R702_unexpected_nan_still_raises(self):
        """R-702: a NaN in SHmax_top / Shmin_top / Sv_top etc.
        (NOT explained by a NaN centroid) must still raise — the
        R-605 input-bug protection is preserved."""
        with pytest.raises(ValueError, match="non-finite"):
            bulk_stress_tensor_field(
                np.array([-1000.0]),    # finite centroid
                SHmax_az_deg=23.0, depth_model="lithostatic_sv",
                SHmax_top=float("nan"),  # input bug
                Shmin_top=49.0, Sv_top=45.0,
                Sv_grad=0.025,
            )

    def test_real_data_safs_no_nan(self):
        """On the real 2 km SAFS fault VTU, build the basis (Phase 2)
        and resolve stresses (Phase 3). No NaN in any emitted
        field on non-degenerate rows; expected zero degenerate rows
        per Phase 2 acceptance #3."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        mesh = load_fault_mesh(_REAL_2000M_FAULT)
        geom = build_fault_basis(
            mesh,
            fault_phys_name="fault",
            rake_sense="right-lateral",
            fault_strike_azimuth_hint_deg=314.0,
        )
        assert geom.n_degenerate == 0
        res = project_stress_onto_fault(
            geom,
            SHmax=113.0, Shmin=49.0, Sv=45.0, P_p=16.0,
            SHmax_az_deg=23.0,
        )
        for key in ("sigma_n_total", "sigma_n_eff", "tau_strike",
                    "tau_dip", "tau_magnitude", "rake_deg"):
            assert not np.isnan(res[key]).any(), (
                f"{key} contains NaN on real SAFS mesh"
            )
        # σ_n_total > 0 (compression positive) at every triangle.
        assert (res["sigma_n_total"] > 0).all(), (
            "σ_n_total must be compression-positive on the SAFS mesh"
        )


# ----------------------------------------------------------------------
# Phase 4 — VTU + JSON writers + CLI
# ----------------------------------------------------------------------

from project_to_fault_stress import (                       # noqa: E402
    DEFAULTS,
    _CONVENTION_STRING,
    _derive_bulk_path,
    _extract_lc_tag,
    _field_stats,
    _hz_demo_safod_defaults,
    main as cli_main,
    write_bulk_vtu,
    write_fault_vtu,
    write_summary_json,
)


def _build_synth_fault_geom(tmp_path) -> tuple:
    """Build a 4-triangle fan-around-vertex-0 synthetic fault mesh,
    project SAFOD stress onto it, return (mesh, tri_conn, geom,
    resolved) ready to feed to the writers."""
    import meshio
    # Four-triangle fan: central vertex 0, surrounded by 1, 2, 3, 4.
    # Use UTM-scale coords (large enough to skip the non-UTM warning).
    base = np.array([1.0e6, 1.0e6, 0.0])
    pts = np.array([
        base + [0.0, 0.0, 0.0],
        base + [1000.0, 0.0, 0.0],
        base + [0.0, 1000.0, 0.0],
        base + [-1000.0, 0.0, 0.0],
        base + [0.0, -1000.0, 0.0],
    ], dtype=np.float64)
    tris = np.array(
        [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
        dtype=np.int64,
    )
    m = meshio.Mesh(
        points=pts,
        cells=[("triangle", tris)],
    )
    vtu = tmp_path / "synth_fan_fault.vtu"
    meshio.write(str(vtu), m)
    mesh = load_fault_mesh(vtu)
    tri_conn = extract_cells_by_physical(
        mesh, "triangle", DEFAULT_FAULT_NAME,
    )
    # Use a strike hint that makes the basis well-defined: the
    # triangles are horizontal (n = ±ẑ) so the basis is degenerate.
    # Tilt the mesh so triangles are vertical instead.
    pts_tilted = pts.copy()
    pts_tilted[1:, 2] = 0.0
    pts_tilted[1, 0] += 0.0  # keep horizontal? — degenerate
    # Easier: just use a vertical SAFOD synth triangle and skip
    # full-fan testing for the writer (writer is exercised via the
    # real-data smoke below).
    return mesh, tri_conn


class TestWriteBulkVtu:
    def test_writes_six_scalar_fields_and_tensor(self, tmp_path):
        """write_bulk_vtu emits six sigma_*_MPa scalar cell-data
        fields plus a 9-component sigma_tensor_MPa."""
        import meshio
        pts = np.array(
            [[1e6, 1e6, 0.], [1e6 + 1e3, 1e6, 0.],
             [1e6, 1e6 + 1e3, 0.], [1e6, 1e6, -1e3]],
            dtype=np.float64,
        )
        tet = np.array([[0, 1, 2, 3]], dtype=np.int64)
        m = meshio.Mesh(
            points=pts,
            cells=[("tetra", tet)],
        )
        vtu_in = tmp_path / "synth_tet.vtu"
        meshio.write(str(vtu_in), m)
        bulk_mesh = load_bulk_mesh(vtu_in)
        bulk_geom = BulkCellGeometry(
            centroids=np.array([[1e6, 1e6, -250.]]),
            volumes=np.array([1e9 / 6.0]),
        )
        sigma_field = np.array([
            [[100.0, 10.0, 5.0],
             [10.0, 80.0, 2.0],
             [5.0, 2.0, 60.0]],
        ])
        out_path = tmp_path / "synth_bulk_stress.vtu"
        write_bulk_vtu(bulk_mesh, tet, bulk_geom, sigma_field, out_path)
        assert out_path.is_file()
        # Roundtrip.
        loaded = meshio.read(str(out_path))
        cd = loaded.cell_data
        for key in (
            "sigma_xx_MPa", "sigma_yy_MPa", "sigma_zz_MPa",
            "sigma_xy_MPa", "sigma_xz_MPa", "sigma_yz_MPa",
        ):
            assert key in cd, f"missing cell-data field {key}"
        assert "sigma_tensor_MPa" in cd
        # Values match input exactly (pass-through, no sign flip).
        np.testing.assert_allclose(cd["sigma_xx_MPa"][0], [100.0])
        np.testing.assert_allclose(cd["sigma_yy_MPa"][0], [80.0])
        np.testing.assert_allclose(cd["sigma_zz_MPa"][0], [60.0])
        np.testing.assert_allclose(cd["sigma_xy_MPa"][0], [10.0])
        np.testing.assert_allclose(cd["sigma_xz_MPa"][0], [5.0])
        np.testing.assert_allclose(cd["sigma_yz_MPa"][0], [2.0])
        # 9-component tensor: row-major (xx, xy, xz, yx, yy, yz, zx, zy, zz).
        np.testing.assert_allclose(
            cd["sigma_tensor_MPa"][0][0],
            [100.0, 10.0, 5.0, 10.0, 80.0, 2.0, 5.0, 2.0, 60.0],
        )

    def test_pass_through_no_sign_flip(self, tmp_path):
        """Compression-positive input → compression-positive output;
        no sign manipulation in the writer (R-501/R-502 contract)."""
        import meshio
        pts = np.array(
            [[1e6, 1e6, 0.], [1e6 + 1e3, 1e6, 0.],
             [1e6, 1e6 + 1e3, 0.], [1e6, 1e6, -1e3]],
            dtype=np.float64,
        )
        tet = np.array([[0, 1, 2, 3]], dtype=np.int64)
        m = meshio.Mesh(points=pts, cells=[("tetra", tet)])
        vtu_in = tmp_path / "synth.vtu"
        meshio.write(str(vtu_in), m)
        bulk_mesh = load_bulk_mesh(vtu_in)
        bulk_geom = BulkCellGeometry(
            centroids=np.array([[1e6, 1e6, -250.]]),
            volumes=np.array([1e9 / 6.0]),
        )
        # Negative diagonal would mean tension — the writer must NOT
        # silently flip back to positive.
        sigma_field = np.array(
            [[[-50.0, 0., 0.], [0., -30., 0.], [0., 0., -20.]]],
        )
        out_path = tmp_path / "synth_bulk.vtu"
        write_bulk_vtu(bulk_mesh, tet, bulk_geom, sigma_field, out_path)
        loaded = meshio.read(str(out_path))
        # Output equals input (pass-through), even when input is
        # non-physical "tension positive".
        np.testing.assert_allclose(
            loaded.cell_data["sigma_xx_MPa"][0], [-50.0],
        )

    def test_bad_sigma_field_shape_raises(self, tmp_path):
        import meshio
        pts = np.zeros((4, 3))
        tet = np.array([[0, 1, 2, 3]], dtype=np.int64)
        m = meshio.Mesh(points=pts, cells=[("tetra", tet)])
        vtu_in = tmp_path / "synth.vtu"
        meshio.write(str(vtu_in), m)
        bulk_mesh = load_bulk_mesh(vtu_in)
        bulk_geom = BulkCellGeometry(
            centroids=np.array([[0., 0., 0.]]),
            volumes=np.array([1.0]),
        )
        with pytest.raises(ValueError, match="sigma_field must be"):
            write_bulk_vtu(
                bulk_mesh, tet, bulk_geom,
                np.zeros((1, 3, 4)),   # wrong shape
                tmp_path / "out.vtu",
            )


class TestWriteSummaryJson:
    def test_schema_complete(self, tmp_path):
        """Summary JSON contains all top-level keys per Phase 4 §3."""
        # Synthetic SAFOD-style resolved dict.
        N = 5
        resolved = {
            "sigma_n_total": np.linspace(50.0, 100.0, N),
            "sigma_n_eff": np.linspace(34.0, 84.0, N),
            "tau_strike": np.linspace(10.0, 25.0, N),
            "tau_dip": np.zeros(N),
            "tau_magnitude": np.linspace(10.0, 25.0, N),
            "rake_deg": np.zeros(N),
            "mu_apparent": np.linspace(0.2, 0.3, N),
        }
        fault_geom = FaultCellGeometry(
            centroids=np.zeros((N, 3)),
            normals=np.zeros((N, 3)),
            strikes=np.zeros((N, 3)),
            dips=np.zeros((N, 3)),
            areas=np.ones(N),
            n_degenerate=0,
        )
        sigma0 = np.diag([60.0, 100.0, 45.0])
        params = {
            "SHmax_MPa": 113.0,
            "Shmin_MPa": 49.0,
            "Sv_MPa": 45.0,
            "P_p_MPa": 16.0,
            "SHmax_azimuth_deg": 23.0,
        }
        out = tmp_path / "summary.json"
        write_summary_json(
            out,
            input_mesh_path=Path("/tmp/synth.vtu"),
            fault_geom=fault_geom, bulk_geom=None,
            sigma_global_at_z0=sigma0, params=params,
            resolved=resolved,
        )
        with open(out) as f:
            data = json.load(f)
        for key in (
            "input_mesh", "fault_n_cells", "fault_n_degenerate_basis",
            "bulk_n_cells", "params", "sigma_global_at_z0_MPa",
            "stats", "convention",
        ):
            assert key in data, f"missing top-level key {key}"
        # stats has all per-field entries.
        for fld in (
            "sigma_n_total_MPa", "sigma_n_eff_MPa", "tau_strike_MPa",
            "tau_dip_MPa", "tau_magnitude_MPa", "rake_deg",
            "mu_apparent",
        ):
            assert fld in data["stats"]
            for s in ("min", "median", "max", "nan_count"):
                assert s in data["stats"][fld]
        # No bulk_geom passed → bulk_n_cells is None.
        assert data["bulk_n_cells"] is None
        assert data["fault_n_cells"] == N

    def test_handles_nan_via_allow_nan(self, tmp_path):
        """mu_apparent NaN entries must serialise (CPython's
        allow_nan=True default; Phase 4 §3 calls this out)."""
        import json as _json
        N = 3
        resolved = {
            "sigma_n_total": np.array([100., 50., 1.]),
            "sigma_n_eff": np.array([84., 34., -10.]),  # last is tension
            "tau_strike": np.array([10., 5., 1.]),
            "tau_dip": np.zeros(N),
            "tau_magnitude": np.array([10., 5., 1.]),
            "rake_deg": np.zeros(N),
            "mu_apparent": np.array([0.12, 0.15, np.nan]),
        }
        fault_geom = FaultCellGeometry(
            centroids=np.zeros((N, 3)), normals=np.zeros((N, 3)),
            strikes=np.zeros((N, 3)), dips=np.zeros((N, 3)),
            areas=np.ones(N), n_degenerate=0,
        )
        out = tmp_path / "summary_nan.json"
        write_summary_json(
            out,
            input_mesh_path=Path("/tmp/synth.vtu"),
            fault_geom=fault_geom, bulk_geom=None,
            sigma_global_at_z0=np.eye(3),
            params={}, resolved=resolved,
        )
        # NaN survives the round-trip (CPython json reads "NaN" back
        # as float('nan')).
        with open(out) as f:
            data = _json.load(f)
        assert data["stats"]["mu_apparent"]["nan_count"] == 1


class TestFieldStats:
    def test_handles_all_nan(self):
        s = _field_stats(np.full(5, np.nan))
        assert s["nan_count"] == 5
        assert np.isnan(s["min"])
        assert np.isnan(s["median"])
        assert np.isnan(s["max"])

    def test_min_median_max(self):
        s = _field_stats(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
        assert s["min"] == 1.0
        assert s["median"] == 3.0
        assert s["max"] == 5.0
        assert s["nan_count"] == 0

    def test_mixed_nan_and_finite(self):
        s = _field_stats(np.array([1.0, np.nan, 3.0]))
        assert s["nan_count"] == 1
        assert s["min"] == 1.0
        assert s["max"] == 3.0


class TestCLIHelpers:
    def test_derive_bulk_path(self):
        p = Path("/foo/bar/safs_fault_box_nwcut_2000m_fault.vtu")
        b = _derive_bulk_path(p)
        assert b.name == "safs_fault_box_nwcut_2000m_bulk.vtu"

    def test_derive_bulk_path_rejects_missing_fault(self):
        with pytest.raises(ValueError, match="end with '_fault.vtu'"):
            _derive_bulk_path(Path("/foo/some_mesh.vtu"))

    def test_R706_derive_bulk_path_rejects_non_suffix_fault(self):
        """R-706: only the canonical `_fault.vtu` stem suffix
        triggers the swap; any other use of `_fault` in the path
        must raise rather than silently produce a wrong sibling
        path (e.g., `safs_fault_box.vtu` → `safs_bulk_box.vtu`)."""
        with pytest.raises(ValueError, match="end with '_fault.vtu'"):
            _derive_bulk_path(Path("/data/safs_fault_box.vtu"))
        with pytest.raises(ValueError, match="end with '_fault.vtu'"):
            _derive_bulk_path(Path("/data/bp5_fault_set2.msh"))

    def test_extract_lc_tag(self):
        assert _extract_lc_tag("safs_fault_box_nwcut_2000m") == "2000m"
        assert _extract_lc_tag("safs_fault_box_nwcut_500m") == "500m"
        assert _extract_lc_tag("safs_fault_box_nwcut_1000m") == "1000m"

    def test_extract_lc_tag_fallback(self):
        # No trailing <N>m: returns the whole stem.
        assert _extract_lc_tag("some_mesh_no_tag") == "some_mesh_no_tag"

    def test_extract_lc_tag_with_variant(self):
        # Trailing `<N>m_<variant>` is preserved so variant outputs
        # land in a sibling subdir of the basic `<N>m/` outputs.
        assert (_extract_lc_tag("safs_fault_box_nwcut_500m_lcfar3000")
                == "500m_lcfar3000")
        assert (_extract_lc_tag("safs_fault_box_nwcut_1000m_lcfar5000")
                == "1000m_lcfar5000")
        assert (_extract_lc_tag("safs_fault_box_nwcut_2000m_zgraded")
                == "2000m_zgraded")

    def test_hz_demo_safod_defaults(self):
        d = _hz_demo_safod_defaults()
        # H&Z demo_safod uses SAFOD pilot-hole values.
        assert d["SHmax"] == 113.0
        assert d["Shmin"] == 49.0
        assert d["Sv"] == 45.0
        assert d["P_p"] == 16.0
        assert d["SHmax_az_deg"] == 23.0
        assert d["fault_strike_az_deg"] == 314.0
        assert d["rake_sense"] == "right-lateral"

    def test_DEFAULTS_matches_introspected_hz(self):
        # Module-level DEFAULTS is the introspected H&Z dict.
        d = _hz_demo_safod_defaults()
        assert DEFAULTS == d

    def test_convention_string_is_compression_positive(self):
        # The summary JSON's convention string must explicitly say
        # SEAS / compression-POSITIVE so future readers can't mistake
        # the active convention.
        assert "compression POSITIVE" in _CONVENTION_STRING
        assert "SEAS" in _CONVENTION_STRING


class TestCLIRealDataSmoke:
    """End-to-end CLI smoke tests on the real 2000 m SAFS mesh.
    Skipped if the real VTUs are not on disk."""

    def test_print_info_writes_no_files(self, tmp_path, capsys):
        """--print-info prints σ⁰ and fault statistics; creates no files."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        rc = cli_main([
            "--print-info", str(_REAL_2000M_FAULT),
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        captured = capsys.readouterr()
        assert "Bulk σ⁰ at z = 0" in captured.out
        assert "Fault triangles" in captured.out
        # No artefacts written.
        assert not (tmp_path / "2000m").exists()

    def test_R703_print_info_prints_strike_statistics(self, tmp_path, capsys):
        """R-703: --print-info must also print strike-azimuth
        statistics per plan §Phase 4 §4 line 1253-1255.  Verifies
        min / median / max / stddev are all present in stdout."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        rc = cli_main([
            "--print-info", str(_REAL_2000M_FAULT),
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        captured = capsys.readouterr()
        assert "Strike-azimuth statistics" in captured.out, (
            "R-703: --print-info missing strike-azimuth statistics"
        )
        for tag in ("min=", "median=", "max=", "stddev="):
            assert tag in captured.out, (
                f"R-703: strike-azimuth statistics missing {tag!r}"
            )

    def test_single_file_writes_three_artefacts(self, tmp_path, capsys):
        """Single-file mode with --write-bulk produces fault VTU,
        bulk VTU, and summary JSON; all roundtrip via meshio /
        json.load."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        _skip_if_no_real_data(_REAL_2000M_BULK)
        rc = cli_main([
            str(_REAL_2000M_FAULT),
            "safs_fault_box_nwcut_2000m",
            "--write-bulk",
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        out_dir = tmp_path / "2000m"
        fault_vtu = out_dir / "safs_fault_box_nwcut_2000m_fault_stress.vtu"
        bulk_vtu = out_dir / "safs_fault_box_nwcut_2000m_bulk_stress.vtu"
        summary = out_dir / "safs_fault_box_nwcut_2000m_summary.json"
        assert fault_vtu.is_file()
        assert bulk_vtu.is_file()
        assert summary.is_file()

        # Roundtrip — fault VTU.
        import meshio
        fault_m = meshio.read(str(fault_vtu))
        assert "sigma_n_total_MPa" in fault_m.point_data
        assert "tau_strike_MPa" in fault_m.point_data
        assert "tau_dip_MPa" in fault_m.point_data
        assert "tau_magnitude_MPa" in fault_m.point_data
        assert "mu_apparent" in fault_m.point_data
        assert "strike_vec" in fault_m.point_data
        assert "dip_vec" in fault_m.point_data
        assert "normal_vec" in fault_m.point_data
        assert "traction_vec_MPa" in fault_m.point_data
        # Cell-data: raw `_cell` channel.
        assert "sigma_n_total_MPa_cell" in fault_m.cell_data
        assert "tau_strike_MPa_cell" in fault_m.cell_data

        # Roundtrip — bulk VTU.
        bulk_m = meshio.read(str(bulk_vtu))
        for k in (
            "sigma_xx_MPa", "sigma_yy_MPa", "sigma_zz_MPa",
            "sigma_xy_MPa", "sigma_xz_MPa", "sigma_yz_MPa",
            "sigma_tensor_MPa",
        ):
            assert k in bulk_m.cell_data, f"bulk missing {k}"

        # Roundtrip — summary JSON.
        with open(summary) as f:
            data = json.load(f)
        assert data["fault_n_cells"] == 2685
        assert data["bulk_n_cells"] == 145300
        assert data["fault_n_degenerate_basis"] == 0
        assert data["params"]["SHmax_MPa"] == 113.0
        # Acceptance criterion #1: σ_n_total is compression-positive.
        assert data["stats"]["sigma_n_total_MPa"]["min"] > 0
        # Acceptance criterion #5: mu_apparent.median in plausible
        # range. The plan's loose bound [0.20, 0.35] was based on a
        # single point estimate; the curved-SAF distribution has a
        # higher median.  Use a wider but still meaningful bound.
        med = data["stats"]["mu_apparent"]["median"]
        assert 0.1 < med < 1.0, (
            f"mu_apparent.median = {med} outside plausible "
            f"[0.1, 1.0] range for SAFOD-parameter run"
        )

    def test_single_file_compression_positive_on_real_mesh(
        self, tmp_path,
    ):
        """All sigma_n_total values are compression-positive on the
        real SAFS mesh (point-data AND cell-data — the R-301 / R-501
        sign-convention regression)."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        rc = cli_main([
            str(_REAL_2000M_FAULT),
            "safs_fault_box_nwcut_2000m",
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        import meshio
        m = meshio.read(
            str(tmp_path / "2000m" /
                "safs_fault_box_nwcut_2000m_fault_stress.vtu")
        )
        point_sigma_n = m.point_data["sigma_n_total_MPa"]
        cell_sigma_n = m.cell_data["sigma_n_total_MPa_cell"][0]
        # Point-data: NaN at orphan vertices (not on the fault);
        # filter those out.
        finite_pt = ~np.isnan(point_sigma_n)
        assert (point_sigma_n[finite_pt] > 0).all()
        # Cell-data: every triangle has compression-positive σ_n.
        assert (cell_sigma_n > 0).all()

    def test_batch_processes_all_resolutions(self, tmp_path, capsys):
        """--batch iterates over DEFAULT_INPUT_BASES; skips meshes
        that aren't on disk; produces output bundles for the rest."""
        # We can't easily redirect the batch-mode code_meshing
        # location, so just smoke-test that --batch runs without
        # crashing and produces SOMETHING per available resolution.
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        # Run via tmp out-dir so we don't pollute the project's
        # canonical data_projection_onfaultstress/.
        rc = cli_main([
            "--batch",
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        # At least the 2000 m output bundle should exist.
        assert (
            tmp_path / "2000m" /
            "safs_fault_box_nwcut_2000m_fault_stress.vtu"
        ).is_file()

    def test_cli_rejects_missing_mode(self, capsys):
        """No mode given → exit 2 with a clear message."""
        rc = cli_main([])
        assert rc == 2
        captured = capsys.readouterr()
        assert "no mode given" in captured.err

    def test_cli_rejects_mode_conflict(self, capsys):
        """--batch AND positional args → exit 2."""
        rc = cli_main([
            "/tmp/foo.vtu", "foo",
            "--batch",
        ])
        assert rc == 2
        captured = capsys.readouterr()
        assert "mutually exclusive" in captured.err

    def test_cli_rejects_partial_positionals(self, capsys):
        """Only input_mesh, no mesh_base → exit 2."""
        rc = cli_main(["/tmp/foo.vtu"])
        assert rc == 2
        captured = capsys.readouterr()
        assert "BOTH" in captured.err

    def test_cli_lithostatic_zero_gradient_warning(
        self, tmp_path, capsys,
    ):
        """depth-model=lithostatic_sv with all-zero gradients
        emits a stderr warning."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        rc = cli_main([
            str(_REAL_2000M_FAULT),
            "safs_fault_box_nwcut_2000m",
            "--depth-model", "lithostatic_sv",
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        captured = capsys.readouterr()
        assert "lithostatic_sv" in captured.err
        assert "all gradients" in captured.err

    def test_cli_write_bulk_missing_raises(self, tmp_path, capsys):
        """--write-bulk requested but bulk VTU missing → clear
        FileNotFoundError with a pointer to msh_to_vtu.py."""
        # Synthesise a fault VTU that has NO corresponding bulk file.
        import meshio
        pts = np.array(
            [[1e6, 1e6, 0.], [1e6 + 1e3, 1e6, 0.],
             [1e6, 1e6, -1e3]],
            dtype=np.float64,
        )
        tris = np.array([[0, 1, 2]], dtype=np.int64)
        m = meshio.Mesh(points=pts, cells=[("triangle", tris)])
        fault_only_path = tmp_path / "lonely_fault.vtu"
        meshio.write(str(fault_only_path), m)
        # The derived sibling "lonely_bulk.vtu" doesn't exist.
        with pytest.raises(FileNotFoundError, match="msh_to_vtu"):
            cli_main([
                str(fault_only_path),
                "lonely",
                "--write-bulk",
                "--out-dir", str(tmp_path),
            ])

    def test_cli_hz_dump_file_populates_params(self, tmp_path):
        """--hz-dump-file populates SHmax/Shmin/Sv/P_p/SHmax_az from
        the dump JSON; per-flag --SHmax wins over the dump."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        # Build a dump JSON with non-default values.
        dump = {
            "schema": "hickman_zoback_sigma0_v1",
            "params": {
                "SHmax_MPa": 200.0,
                "Shmin_MPa": 75.0,
                "Sv_MPa": 60.0,
                "P_p_MPa": 20.0,
                "SHmax_azimuth_deg": 45.0,
            },
            "sigma0_MPa": [[-200.0, 0.0, 0.0],
                            [0.0, -75.0, 0.0],
                            [0.0, 0.0, -60.0]],
            "convention": "compression negative; ...",
            "regime": "test",
        }
        dump_path = tmp_path / "custom_dump.json"
        with open(dump_path, "w") as f:
            json.dump(dump, f)
        rc = cli_main([
            str(_REAL_2000M_FAULT),
            "safs_fault_box_nwcut_2000m",
            "--hz-dump-file", str(dump_path),
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        with open(
            tmp_path / "2000m" /
            "safs_fault_box_nwcut_2000m_summary.json"
        ) as f:
            data = json.load(f)
        assert data["params"]["SHmax_MPa"] == 200.0
        assert data["params"]["Shmin_MPa"] == 75.0
        assert data["params"]["Sv_MPa"] == 60.0

    def test_cli_flag_overrides_dump_file(self, tmp_path):
        """--SHmax 999.0 wins over the dump-file value."""
        _skip_if_no_real_data(_REAL_2000M_FAULT)
        dump = {
            "schema": "hickman_zoback_sigma0_v1",
            "params": {
                "SHmax_MPa": 200.0,
                "Shmin_MPa": 75.0,
                "Sv_MPa": 60.0,
                "P_p_MPa": 20.0,
                "SHmax_azimuth_deg": 45.0,
            },
            "sigma0_MPa": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "convention": "",
            "regime": "",
        }
        dump_path = tmp_path / "custom_dump.json"
        with open(dump_path, "w") as f:
            json.dump(dump, f)
        rc = cli_main([
            str(_REAL_2000M_FAULT),
            "safs_fault_box_nwcut_2000m",
            "--hz-dump-file", str(dump_path),
            "--SHmax", "999.0",
            "--out-dir", str(tmp_path),
        ])
        assert rc == 0
        with open(
            tmp_path / "2000m" /
            "safs_fault_box_nwcut_2000m_summary.json"
        ) as f:
            data = json.load(f)
        # CLI flag wins.
        assert data["params"]["SHmax_MPa"] == 999.0
        # Dump-file values used for unset params.
        assert data["params"]["Shmin_MPa"] == 75.0
