"""test_build_stress_safs.py — Phase 5 pytest companion.

Covers:
- ``build_uniform_utm_grid`` (Phase 5 §3): monotone-increasing axes,
  enclosure guarantee, malformed-bbox raises.
- ``evaluate_stress_field_on_grid`` (Phase 5 §4): six-field output,
  shape, MPa→Pa unit conversion (no sign flip), x/y constancy,
  symmetry of the (Nz, 3, 3) source stack.
- ``derive_field_bounds`` (Phase 5 §5): rule semantics, every cell in
  ``[v_min, v_max]``, ``v_min < v_max`` (writer-side guard).
- ``build_stress_safs`` (Phase 5 §6): end-to-end round-trip — write the
  sidecar, reload via h5py, assert schema-v1 attrs, exact six-field
  set, and the grid-strictly-encloses-mesh-bbox acceptance.
- CLI (Phase 5 §7): synth-mesh smoke and real-data smoke (skipped
  when the 500/1000/2000 m bulk VTUs are not on disk).

Run via::

    cd project_7.0_alternative/stress/code && pytest -q test_build_stress_safs.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import meshio
import numpy as np
import pytest


_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from build_stress_safs import (                              # noqa: E402
    DEFAULT_BOUNDS_FACTOR,
    DEFAULT_GRID_DX_M,
    DEFAULT_GRID_PAD_M,
    _FIELD_NAMES,
    _mesh_bbox,
    _union_bbox,
    build_stress_safs,
    build_uniform_utm_grid,
    derive_field_bounds,
    evaluate_stress_field_on_grid,
    main as cli_main,
)
from project_to_fault_stress import bulk_stress_tensor_field  # noqa: E402


# Real-data fixtures — the SAFS bulk VTUs that ship with the repo.
_CODE_MESHING_DIR = (
    _THIS_DIR.parent.parent / "meshing" / "results" / "vtu"
)
_REAL_500M_BULK = (
    _CODE_MESHING_DIR / "safs_fault_box_nwcut_500m_bulk.vtu"
)
_REAL_1000M_BULK = (
    _CODE_MESHING_DIR / "safs_fault_box_nwcut_1000m_bulk.vtu"
)
_REAL_2000M_BULK = (
    _CODE_MESHING_DIR / "safs_fault_box_nwcut_2000m_bulk.vtu"
)


def _skip_if_no_real_data(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"real-data fixture not on disk: {path}")


# ---------------------------------------------------------------------
# Synthetic-mesh fixtures
# ---------------------------------------------------------------------


def _write_synth_bulk_vtu(
    tmp_path: Path,
    name: str = "synth_bulk.vtu",
    bbox: tuple[float, float, float, float, float, float] = (
        1.0e6, 1.0e6 + 1.0e4,        # x in [1.0e6, 1.01e6]
        2.0e6, 2.0e6 + 1.0e4,        # y in [2.0e6, 2.01e6]
        -3.0e3, 1.0e2,               # z in [-3000, 100]
    ),
) -> Path:
    """Write a single-tet bulk VTU whose vertex cloud spans ``bbox``
    on all axes (the four corners of the tet are chosen so that
    the bbox is reached on every axis)."""
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    pts = np.array(
        [
            [xmin, ymin, zmin],   # 0: low corner
            [xmax, ymin, zmin],   # 1: extends x
            [xmin, ymax, zmin],   # 2: extends y
            [xmin, ymin, zmax],   # 3: extends z
            [xmax, ymax, zmax],   # 4: high corner
        ],
        dtype=np.float64,
    )
    # Two non-degenerate tets covering all five points.
    tets = np.array(
        [[0, 1, 2, 3], [1, 2, 3, 4]],
        dtype=np.int64,
    )
    m = meshio.Mesh(
        points=pts,
        cells=[("tetra", tets)],
    )
    p = tmp_path / name
    meshio.write(str(p), m)
    return p


# ---------------------------------------------------------------------
# Module constants — regression guards
# ---------------------------------------------------------------------


class TestModuleConstants:
    def test_R801_symmetry_tol_matches_phase3(self):
        """R-801: Phase 5's _SIGMA_SYMMETRY_TOL_PA must equal Phase 3's
        _SIGMA_SYMMETRY_TOL_MPA after the MPa → Pa unit conversion."""
        import build_stress_safs as m5
        from project_to_fault_stress import _SIGMA_SYMMETRY_TOL_MPA
        np.testing.assert_allclose(
            m5._SIGMA_SYMMETRY_TOL_PA,
            _SIGMA_SYMMETRY_TOL_MPA * 1.0e6,
            rtol=0.0, atol=0.0,
        )

    def test_R802_default_out_under_stress_results(self):
        """R-802 (post-restructure): Phase 5's DEFAULT_OUT must resolve
        under stress/results/, the canonical sibling of stress/code/
        after the project_7.0_alternative pillar reorganisation."""
        from build_stress_safs import DEFAULT_OUT
        assert DEFAULT_OUT.name == "stress_safs.h5"
        assert DEFAULT_OUT.parent.name == "results"
        assert DEFAULT_OUT.parent.parent.name == "stress"


# ---------------------------------------------------------------------
# build_uniform_utm_grid
# ---------------------------------------------------------------------


class TestBuildUniformUtmGrid:
    def test_axes_are_strictly_monotone_increasing(self):
        bbox = (1.0e6, 1.05e6, 2.0e6, 2.05e6, -5.0e4, 1.0e3)
        x, y, z = build_uniform_utm_grid(bbox, dx_m=1000.0, pad_m=2000.0)
        for axis in (x, y, z):
            assert axis.size >= 2
            assert np.all(np.diff(axis) > 0)
            assert axis.dtype == np.float64

    def test_strictly_encloses_bbox_by_pad(self):
        bbox = (1.0e6, 1.05e6, 2.0e6, 2.05e6, -5.0e4, 1.0e3)
        dx, pad = 1000.0, 2000.0
        x, y, z = build_uniform_utm_grid(bbox, dx_m=dx, pad_m=pad)
        # Grid extends below bbox.min - pad and above bbox.max + pad.
        assert x[0] <= bbox[0] - pad
        assert x[-1] >= bbox[1] + pad
        assert y[0] <= bbox[2] - pad
        assert y[-1] >= bbox[3] + pad
        assert z[0] <= bbox[4] - pad
        assert z[-1] >= bbox[5] + pad

    def test_uniform_spacing(self):
        bbox = (0.0, 1.0e3, 0.0, 1.0e3, 0.0, 1.0e3)
        # bbox doesn't satisfy "subsurface" convention but the helper
        # is generic — only checks min < max per axis.
        bbox = (1.0e6, 1.0e6 + 1.0e3, 2.0e6, 2.0e6 + 1.0e3, -1.0e3, 1.0e2)
        dx = 100.0
        x, y, z = build_uniform_utm_grid(bbox, dx_m=dx, pad_m=0.0)
        np.testing.assert_allclose(np.diff(x), dx, atol=1e-9)
        np.testing.assert_allclose(np.diff(y), dx, atol=1e-9)
        np.testing.assert_allclose(np.diff(z), dx, atol=1e-9)

    def test_bad_dx_raises(self):
        bbox = (1.0e6, 1.05e6, 2.0e6, 2.05e6, -5.0e4, 1.0e3)
        with pytest.raises(ValueError, match="dx_m must be > 0"):
            build_uniform_utm_grid(bbox, dx_m=0.0, pad_m=1.0)
        with pytest.raises(ValueError, match="dx_m must be > 0"):
            build_uniform_utm_grid(bbox, dx_m=-1.0, pad_m=1.0)

    def test_negative_pad_raises(self):
        bbox = (1.0e6, 1.05e6, 2.0e6, 2.05e6, -5.0e4, 1.0e3)
        with pytest.raises(ValueError, match="pad_m must be >= 0"):
            build_uniform_utm_grid(bbox, dx_m=1000.0, pad_m=-1.0)

    def test_malformed_bbox_raises(self):
        # x: min >= max
        with pytest.raises(ValueError, match="min < max on every axis"):
            build_uniform_utm_grid(
                (1.0e6, 1.0e6, 2.0e6, 2.1e6, -1.0e3, 1.0),
                dx_m=1000.0, pad_m=0.0,
            )
        # z: min >= max
        with pytest.raises(ValueError, match="min < max on every axis"):
            build_uniform_utm_grid(
                (1.0e6, 1.1e6, 2.0e6, 2.1e6, 1.0, -1.0),
                dx_m=1000.0, pad_m=0.0,
            )


# ---------------------------------------------------------------------
# evaluate_stress_field_on_grid
# ---------------------------------------------------------------------


class TestEvaluateStressFieldOnGrid:
    def _grid(self):
        x = np.array([1.0e6, 1.0e6 + 1e3, 1.0e6 + 2e3], dtype=np.float64)
        y = np.array([2.0e6, 2.0e6 + 1e3], dtype=np.float64)
        z = np.array([-3.0e3, -1.5e3, 0.0], dtype=np.float64)
        return x, y, z

    def test_returns_exactly_six_fields(self):
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        assert set(f.keys()) == set(_FIELD_NAMES)
        # Insertion order must match the schema-v1 canonical order.
        assert tuple(f.keys()) == _FIELD_NAMES

    def test_field_shapes(self):
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        expected = (x.size, y.size, z.size)
        for k, arr in f.items():
            assert arr.shape == expected, f"{k} shape {arr.shape}"
            assert arr.dtype == np.float64

    def test_constant_in_x_and_y(self):
        """Homogeneous σ⁰ → field varies only with z; identical across
        every (i, j) at fixed k."""
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        for name, arr in f.items():
            # arr[i, j, k] must equal arr[0, 0, k] for all i, j.
            ref = arr[0, 0, :]
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    np.testing.assert_allclose(
                        arr[i, j, :], ref, atol=0.0, rtol=0.0,
                        err_msg=f"{name}[{i}, {j}, :] differs from [0,0,:]",
                    )

    def test_mpa_to_pa_unit_conversion_no_sign_flip(self):
        """The field at a single z slice must equal
        ``bulk_stress_tensor_field(z) * 1e6`` (Pa) with NO sign flip
        (R-501 / R-502)."""
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        sigma_MPa = bulk_stress_tensor_field(
            z,
            SHmax_az_deg=23.0, depth_model="constant",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
        )
        sigma_Pa = sigma_MPa * 1.0e6
        # sigma_xx[i, j, k] == sigma_Pa[k, 0, 0].
        for i, j in ((0, 0), (1, 0), (0, 1)):
            np.testing.assert_allclose(
                f["sigma_xx"][i, j, :], sigma_Pa[:, 0, 0],
            )
            np.testing.assert_allclose(
                f["sigma_yy"][i, j, :], sigma_Pa[:, 1, 1],
            )
            np.testing.assert_allclose(
                f["sigma_zz"][i, j, :], sigma_Pa[:, 2, 2],
            )
            np.testing.assert_allclose(
                f["sigma_xy"][i, j, :], sigma_Pa[:, 0, 1],
            )
            np.testing.assert_allclose(
                f["sigma_yz"][i, j, :], sigma_Pa[:, 1, 2],
            )
            np.testing.assert_allclose(
                f["sigma_xz"][i, j, :], sigma_Pa[:, 0, 2],
            )

    def test_diagonal_is_compression_positive(self):
        """SAFOD parameters: every sigma_xx / yy / zz cell positive
        (compression POSITIVE SEAS)."""
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        for k in ("sigma_xx", "sigma_yy", "sigma_zz"):
            assert (f[k] > 0).all(), (
                f"{k} contains non-positive cells; compression POSITIVE "
                f"SEAS convention violated"
            )

    def test_R803_symmetry_on_function_output(self):
        """R-803: the six per-component (Nx, Ny, Nz) arrays returned
        by ``evaluate_stress_field_on_grid`` must satisfy
        ``sigma_xy == sigma_yx``, ``sigma_xz == sigma_zx``,
        ``sigma_yz == sigma_zy`` at every (i, j, k) — verified
        against a fresh ``bulk_stress_tensor_field`` rotation as
        ground truth."""
        x, y, z = self._grid()
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        sigma_Pa = bulk_stress_tensor_field(
            z,
            SHmax_az_deg=23.0, depth_model="constant",
            SHmax_top=113.0, Shmin_top=49.0, Sv_top=45.0,
        ) * 1.0e6
        # Every (i, j) slice of the output equals the upper-triangle
        # of σ_seas at that k. Asserting on the function's output
        # (not on the source) catches a regression in this module's
        # broadcast / dict-construction logic.
        for i in range(f["sigma_xy"].shape[0]):
            for j in range(f["sigma_xy"].shape[1]):
                np.testing.assert_allclose(
                    f["sigma_xy"][i, j, :], sigma_Pa[:, 1, 0],
                    atol=1e-3,
                )
                np.testing.assert_allclose(
                    f["sigma_xz"][i, j, :], sigma_Pa[:, 2, 0],
                    atol=1e-3,
                )
                np.testing.assert_allclose(
                    f["sigma_yz"][i, j, :], sigma_Pa[:, 2, 1],
                    atol=1e-3,
                )

    def test_lithostatic_sv_depth_dependence(self):
        """With non-zero gradients, sigma_zz at z = -1000 m exceeds
        sigma_zz at z = 0."""
        x = np.array([1.0e6], dtype=np.float64)
        y = np.array([2.0e6], dtype=np.float64)
        z = np.array([-1.0e3, 0.0], dtype=np.float64)
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            depth_model="lithostatic_sv", Sv_grad=0.025,
        )
        # sigma_zz at -1000 m = (45 + 0.025*1000) MPa = 70 MPa = 7.0e7 Pa
        np.testing.assert_allclose(
            f["sigma_zz"][0, 0, 0], 7.0e7, atol=1.0,
        )
        # sigma_zz at 0 m = 45 MPa = 4.5e7 Pa
        np.testing.assert_allclose(
            f["sigma_zz"][0, 0, 1], 4.5e7, atol=1.0,
        )

    def test_R805_lithostatic_SHmax_grad_propagates(self):
        """R-805: non-zero ``SHmax_grad`` must propagate to
        depth-varying ``sigma_xx``.  With SHmax aligned to +x
        (geological azimuth 90° cw from north = east), the rotated
        tensor's xx-entry equals ``SHmax(d)`` directly — no
        off-axis contribution from Shmin/Sv."""
        x = np.array([1.0e6], dtype=np.float64)
        y = np.array([2.0e6], dtype=np.float64)
        z = np.array([-2.0e3, 0.0], dtype=np.float64)
        f = evaluate_stress_field_on_grid(
            x, y, z,
            SHmax=100.0, Shmin=50.0, Sv=40.0,
            SHmax_az_deg=90.0,          # SHmax along +x (east)
            depth_model="lithostatic_sv",
            SHmax_grad=0.030,
        )
        # SHmax at z = -2000: 100 + 0.030 * 2000 = 160 MPa = 1.6e8 Pa
        np.testing.assert_allclose(
            f["sigma_xx"][0, 0, 0], 1.6e8, atol=10.0,
        )
        # SHmax at z = 0: 100 MPa = 1.0e8 Pa
        np.testing.assert_allclose(
            f["sigma_xx"][0, 0, 1], 1.0e8, atol=10.0,
        )

    def test_R804_asymmetric_source_raises(self, monkeypatch):
        """R-804: when ``bulk_stress_tensor_field`` returns
        asymmetric tensors (a regression in Phase 3),
        ``evaluate_stress_field_on_grid`` must raise rather than
        silently emit a non-symmetric σ field to the sidecar."""
        import build_stress_safs as m

        def _asym_field(centroids_z, **kwargs):
            Nz = centroids_z.size
            out = np.zeros((Nz, 3, 3), dtype=np.float64)
            # Symmetric diagonal but DELIBERATELY asymmetric
            # off-diagonal large enough to clear the 1e-3 Pa
            # tolerance after the MPa → Pa conversion.
            out[:, 0, 0] = 1.0
            out[:, 1, 1] = 1.0
            out[:, 2, 2] = 1.0
            out[:, 0, 1] = 0.5
            out[:, 1, 0] = 0.0     # asymmetry = 0.5 MPa = 5e5 Pa
            return out

        monkeypatch.setattr(m, "bulk_stress_tensor_field", _asym_field)
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 1.0])
        z = np.array([-1.0e3, 0.0])
        with pytest.raises(ValueError, match="not symmetric"):
            m.evaluate_stress_field_on_grid(
                x, y, z,
                SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            )

    def test_bad_axis_shapes_raise(self):
        with pytest.raises(ValueError, match="1-D"):
            evaluate_stress_field_on_grid(
                np.zeros((2, 2)), np.zeros(2), np.zeros(2),
                SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            )

    def test_empty_axis_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            evaluate_stress_field_on_grid(
                np.zeros(0), np.zeros(2), np.zeros(2),
                SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            )


# ---------------------------------------------------------------------
# derive_field_bounds
# ---------------------------------------------------------------------


class TestDeriveFieldBounds:
    def test_positive_only_field(self):
        f = {"a": np.array([10.0, 20.0, 30.0])}
        b = derive_field_bounds(f, safety_factor=1.2)
        # min = 10 / 1.2 = 8.333..., max = 30 * 1.2 = 36
        v_min, v_max, units = b["a"]
        assert v_min < 10.0 < 30.0 < v_max
        assert units == "Pa"

    def test_negative_only_field(self):
        f = {"a": np.array([-30.0, -20.0, -10.0])}
        b = derive_field_bounds(f, safety_factor=1.2)
        v_min, v_max, units = b["a"]
        # min = -30 * 1.2 = -36, max = -10 / 1.2 = -8.333...
        assert v_min < -30.0 < -10.0 < v_max
        assert units == "Pa"

    def test_mixed_sign_field(self):
        f = {"a": np.array([-10.0, 0.0, 20.0])}
        b = derive_field_bounds(f, safety_factor=1.2)
        v_min, v_max, _ = b["a"]
        # min = -10 * 1.2 = -12, max = 20 * 1.2 = 24
        assert v_min < -10.0
        assert v_max > 20.0

    def test_bounds_strictly_enclose_data(self):
        rng = np.random.default_rng(0)
        f = {"a": rng.standard_normal(100) * 1.0e7}
        b = derive_field_bounds(f, safety_factor=1.5)
        v_min, v_max, _ = b["a"]
        assert v_min < float(f["a"].min())
        assert v_max > float(f["a"].max())
        assert v_min < v_max

    def test_constant_zero_field_widened(self):
        """Identically-zero field would collapse the bounds to [0, 0];
        the helper widens to ±1 Pa so the writer's strict-inequality
        guard passes."""
        f = {"a": np.zeros(5)}
        b = derive_field_bounds(f, safety_factor=1.2)
        v_min, v_max, _ = b["a"]
        assert v_min < 0.0 < v_max

    def test_units_third_element(self):
        f = {"a": np.array([1.0, 2.0])}
        b = derive_field_bounds(f)
        assert len(b["a"]) == 3
        assert b["a"][2] == "Pa"

    def test_bad_safety_factor_raises(self):
        with pytest.raises(ValueError, match="safety_factor must be > 1.0"):
            derive_field_bounds({"a": np.array([1.0])}, safety_factor=1.0)


# ---------------------------------------------------------------------
# _mesh_bbox / _union_bbox
# ---------------------------------------------------------------------


class TestMeshBboxHelpers:
    def test_mesh_bbox_synth(self, tmp_path):
        bbox = (1.0e6, 1.0e6 + 1e4, 2.0e6, 2.0e6 + 1e4, -3.0e3, 1.0e2)
        p = _write_synth_bulk_vtu(tmp_path, bbox=bbox)
        got = _mesh_bbox(p)
        assert got == pytest.approx(bbox, abs=1e-6)

    def test_mesh_bbox_missing_raises(self):
        with pytest.raises(FileNotFoundError, match="mesh not found"):
            _mesh_bbox(Path("/nonexistent/path/mesh.vtu"))

    def test_union_bbox_pads_outward(self):
        bboxes = [
            (0., 1., 0., 1., 0., 1.),
            (0.5, 2., -1., 0.5, -1., 0.5),
        ]
        u = _union_bbox(bboxes)
        assert u == (0., 2., -1., 1., -1., 1.)

    def test_union_bbox_empty_raises(self):
        with pytest.raises(ValueError, match="empty bbox list"):
            _union_bbox([])


# ---------------------------------------------------------------------
# build_stress_safs — end-to-end
# ---------------------------------------------------------------------


class TestBuildStressSafs:
    def test_synth_roundtrip(self, tmp_path):
        """Full pipeline on a synthetic mesh: write sidecar, reload,
        confirm schema-v1 attrs and six-field set."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh_path],
            out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0,
            SHmax_az_deg=23.0,
            dx_m=1000.0, pad_m=2000.0,
            mesh_tag="synth",
        )
        assert out_path.is_file()
        with h5py.File(out_path, "r") as h5:
            attrs = {
                k: (v.decode() if isinstance(v, bytes) else v)
                for k, v in h5.attrs.items()
            }
            # Required schema-v1 attrs.
            assert attrs["schema_version"] == "data_projection_v1"
            assert attrs["crs"] == "EPSG:32611"
            assert attrs["units"] == "m"
            assert attrs["z_positive"] == "elevation"
            # Optional attrs we set.
            assert attrs["mesh_tag"] == "synth"
            assert "source" in attrs
            assert "source_crs" in attrs
            # Six exact fields, no more, no less.
            assert set(h5["fields"].keys()) == set(_FIELD_NAMES)
            assert len(h5["fields"].keys()) == 6
            # Each field has min_value / max_value / units attrs.
            for name in _FIELD_NAMES:
                ds = h5["fields"][name]
                assert "min_value" in ds.attrs
                assert "max_value" in ds.attrs
                assert "units" in ds.attrs
                assert (
                    ds.attrs["units"].decode()
                    if isinstance(ds.attrs["units"], bytes)
                    else ds.attrs["units"]
                ) == "Pa"

    def test_grid_strictly_encloses_mesh_bbox(self, tmp_path):
        """Acceptance criterion #4: grid bbox encloses each mesh bbox
        by at least pad_m on every face."""
        bbox = (1.0e6, 1.0e6 + 1e4, 2.0e6, 2.0e6 + 1e4, -3.0e3, 1.0e2)
        mesh_path = _write_synth_bulk_vtu(tmp_path, bbox=bbox)
        out_path = tmp_path / "stress_safs.h5"
        pad = 2500.0
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            dx_m=500.0, pad_m=pad, mesh_tag="synth",
        )
        with h5py.File(out_path, "r") as h5:
            x = h5["grid"]["x"][:]
            y = h5["grid"]["y"][:]
            z = h5["grid"]["z"][:]
        assert x[0] <= bbox[0] - pad
        assert x[-1] >= bbox[1] + pad
        assert y[0] <= bbox[2] - pad
        assert y[-1] >= bbox[3] + pad
        assert z[0] <= bbox[4] - pad
        assert z[-1] >= bbox[5] + pad

    def test_empty_mesh_paths_raises(self, tmp_path):
        with pytest.raises(ValueError, match="mesh_paths is empty"):
            build_stress_safs(
                mesh_paths=[], out_path=tmp_path / "x.h5",
                SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            )

    def test_lithostatic_zero_gradient_warning(self, tmp_path, capsys):
        """depth_model=lithostatic_sv with all-zero gradients prints a
        stderr warning (mirrors Phase 3 §Edge Cases)."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            depth_model="lithostatic_sv",
            dx_m=1000.0, pad_m=2000.0, mesh_tag="synth",
        )
        captured = capsys.readouterr()
        assert "lithostatic_sv" in captured.err
        assert "all gradients" in captured.err

    def test_union_bbox_when_multiple_meshes(self, tmp_path):
        """Two meshes with disjoint x-extents: the grid must contain
        the UNION."""
        m1 = _write_synth_bulk_vtu(
            tmp_path, name="m1.vtu",
            bbox=(0.0, 1.0e3, 0.0, 1.0e3, -1.0e3, 0.0),
        )
        m2 = _write_synth_bulk_vtu(
            tmp_path, name="m2.vtu",
            bbox=(5.0e3, 6.0e3, 0.0, 1.0e3, -1.0e3, 0.0),
        )
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[m1, m2], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            dx_m=500.0, pad_m=1000.0, mesh_tag="synth",
        )
        with h5py.File(out_path, "r") as h5:
            x = h5["grid"]["x"][:]
        assert x[0] <= 0.0 - 1000.0
        assert x[-1] >= 6.0e3 + 1000.0

    def test_creates_missing_parent_dir(self, tmp_path):
        """Writer creates the destination directory if missing
        (delegated to sidecar.write_sidecar)."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "nested" / "deep" / "stress_safs.h5"
        assert not out_path.parent.exists()
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
        )
        assert out_path.is_file()


# ---------------------------------------------------------------------
# Cross-phase contract (Phase 5 -> Phase 6 C++ StressField3D)
# ---------------------------------------------------------------------


class TestPhase5Phase6Contract:
    """R-902: pin the cross-phase field-name + canonical-order
    contract.  Phase 6's C++ `StressField3D` reader hard-codes the
    six field names (sigma_xx, sigma_yy, sigma_zz, sigma_xy,
    sigma_yz, sigma_xz) in schema-v1 canonical order; this test
    verifies the Phase 5 writer emits exactly those names so a
    field-name drift between phases is caught at Python-test time.
    """

    # Canonical order — must match the StressField3D's component
    # ordering at io/stress_field_3d.cpp:component_names_() and
    # build_stress_safs._FIELD_NAMES.
    _PHASE6_EXPECTED_FIELDS = (
        "sigma_xx", "sigma_yy", "sigma_zz",
        "sigma_xy", "sigma_yz", "sigma_xz",
    )

    def test_R902_phase5_field_names_match_phase6_reader(self, tmp_path):
        """Phase 5 sidecar must contain exactly the six field names
        expected by Phase 6's StressField3D reader."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            dx_m=1000.0, pad_m=2000.0, mesh_tag="synth",
        )
        with h5py.File(out_path, "r") as h5:
            on_disk = set(h5["fields"].keys())
        expected = set(self._PHASE6_EXPECTED_FIELDS)
        assert on_disk == expected, (
            f"R-902: Phase 5 sidecar field names {sorted(on_disk)} "
            f"do not match Phase 6 StressField3D's expected set "
            f"{sorted(expected)}.  A field-name drift between phases "
            f"would make StressField3D's ctor abort at runtime."
        )

    def test_R902_phase5_field_units_are_Pa(self, tmp_path):
        """Phase 6's StressField3D contract requires Pa-valued
        components (R-501/R-502).  Verify every sigma_* dataset
        carries the 'Pa' units attribute Phase 5 sets."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            dx_m=1000.0, pad_m=2000.0, mesh_tag="synth",
        )
        with h5py.File(out_path, "r") as h5:
            for name in self._PHASE6_EXPECTED_FIELDS:
                ds = h5["fields"][name]
                units = ds.attrs["units"]
                if isinstance(units, bytes):
                    units = units.decode()
                assert units == "Pa", (
                    f"R-902: Phase 5 emitted {name!r} with units "
                    f"{units!r}, expected 'Pa'."
                )

    def test_R902_schema_v1_attrs_present(self, tmp_path):
        """Phase 6's StressField3D ctor delegates to DataField3D for
        the schema-v1 attribute checks (schema_version, crs,
        z_positive).  Verify Phase 5 emits those attrs."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[mesh_path], out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            dx_m=1000.0, pad_m=2000.0, mesh_tag="synth",
        )
        with h5py.File(out_path, "r") as h5:
            attrs = {
                k: (v.decode() if isinstance(v, bytes) else v)
                for k, v in h5.attrs.items()
            }
        assert attrs["schema_version"] == "data_projection_v1"
        assert attrs["crs"] == "EPSG:32611"
        assert attrs["z_positive"] == "elevation"


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


class TestCLI:
    def test_cli_synth_smoke(self, tmp_path):
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        out_path = tmp_path / "stress_safs.h5"
        rc = cli_main([
            "--meshes", str(mesh_path),
            "--out", str(out_path),
            "--SHmax", "113.0",
            "--Shmin", "49.0",
            "--Sv", "45.0",
            "--SHmax-az", "23.0",
        ])
        assert rc == 0
        assert out_path.is_file()

    def test_cli_missing_required_flag_raises_systemexit(self, tmp_path):
        """argparse raises SystemExit when a required flag is missing."""
        mesh_path = _write_synth_bulk_vtu(tmp_path)
        with pytest.raises(SystemExit):
            cli_main([
                "--meshes", str(mesh_path),
                # Missing --SHmax, --Shmin, --Sv, --SHmax-az.
            ])

    def test_cli_error_path_returns_nonzero(self, tmp_path, capsys):
        """A FileNotFoundError on a missing mesh maps to rc=1 with a
        stderr message (matches build_velocity_cvmh.py)."""
        rc = cli_main([
            "--meshes", "/nonexistent/mesh.vtu",
            "--out", str(tmp_path / "x.h5"),
            "--SHmax", "113.0", "--Shmin", "49.0", "--Sv", "45.0",
            "--SHmax-az", "23.0",
        ])
        assert rc == 1
        captured = capsys.readouterr()
        assert "ERROR" in captured.err


# ---------------------------------------------------------------------
# Real-data smoke
# ---------------------------------------------------------------------


class TestRealDataSmoke:
    def test_real_2000m_roundtrip(self, tmp_path):
        """End-to-end on the real 2 km bulk VTU.  Verifies acceptance
        criteria #1 (schema-v1 attrs) and #2 (exactly six sigma_*
        fields)."""
        _skip_if_no_real_data(_REAL_2000M_BULK)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[_REAL_2000M_BULK],
            out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            mesh_tag="safs_fault_box_nwcut",
        )
        assert out_path.is_file()
        with h5py.File(out_path, "r") as h5:
            assert set(h5["fields"].keys()) == set(_FIELD_NAMES)
            x = h5["grid"]["x"][:]
            y = h5["grid"]["y"][:]
            z = h5["grid"]["z"][:]
            attrs = {
                k: (v.decode() if isinstance(v, bytes) else v)
                for k, v in h5.attrs.items()
            }
        assert attrs["schema_version"] == "data_projection_v1"
        assert attrs["crs"] == "EPSG:32611"
        # The grid encloses the real mesh bbox by at least pad_m.
        mesh_bbox = _mesh_bbox(_REAL_2000M_BULK)
        assert x[0] <= mesh_bbox[0] - DEFAULT_GRID_PAD_M
        assert x[-1] >= mesh_bbox[1] + DEFAULT_GRID_PAD_M
        assert y[0] <= mesh_bbox[2] - DEFAULT_GRID_PAD_M
        assert y[-1] >= mesh_bbox[3] + DEFAULT_GRID_PAD_M
        assert z[0] <= mesh_bbox[4] - DEFAULT_GRID_PAD_M
        assert z[-1] >= mesh_bbox[5] + DEFAULT_GRID_PAD_M

    def test_real_500_1000_2000_union(self, tmp_path):
        """Acceptance criterion #4: grid encloses the SMALLEST (500 m)
        mesh by at least pad_m when all three resolutions are unioned.
        Skipped if any of the three real meshes is absent."""
        for p in (_REAL_500M_BULK, _REAL_1000M_BULK, _REAL_2000M_BULK):
            _skip_if_no_real_data(p)
        out_path = tmp_path / "stress_safs.h5"
        build_stress_safs(
            mesh_paths=[_REAL_500M_BULK, _REAL_1000M_BULK, _REAL_2000M_BULK],
            out_path=out_path,
            SHmax=113.0, Shmin=49.0, Sv=45.0, SHmax_az_deg=23.0,
            mesh_tag="safs_fault_box_nwcut",
        )
        assert out_path.is_file()
        with h5py.File(out_path, "r") as h5:
            x = h5["grid"]["x"][:]
            y = h5["grid"]["y"][:]
            z = h5["grid"]["z"][:]
        bbox_500 = _mesh_bbox(_REAL_500M_BULK)
        assert x[0] <= bbox_500[0] - DEFAULT_GRID_PAD_M
        assert x[-1] >= bbox_500[1] + DEFAULT_GRID_PAD_M
        assert y[0] <= bbox_500[2] - DEFAULT_GRID_PAD_M
        assert y[-1] >= bbox_500[3] + DEFAULT_GRID_PAD_M
        assert z[0] <= bbox_500[4] - DEFAULT_GRID_PAD_M
        assert z[-1] >= bbox_500[5] + DEFAULT_GRID_PAD_M
