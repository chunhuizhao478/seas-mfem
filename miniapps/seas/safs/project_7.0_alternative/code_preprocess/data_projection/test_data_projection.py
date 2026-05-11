"""Test suite for data_projection (Phase 2). Test catalog T-2-1..T-2-12.

Run with::

    cd safs/project_7.0_alternative/code_preprocess && \\
        pytest -q data_projection/test_data_projection.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
if str(PARENT) not in sys.path:
    sys.path.insert(0, str(PARENT))

from data_projection import (                          # noqa: E402
    BBoxContainmentError,
    SCHEMA_VERSION,
    assert_grid_contains_mesh,
    geographic_to_utm11n,
    grid_bbox,
    grid_contains_bbox,
    mesh_msh_bbox,
    read_cvmh_ascii,
    read_sidecar_attrs,
    utm11n_to_geographic,
    write_sidecar,
)


# ---------------------------------------------------------------------------
# Synthetic CVM-H slice fixture
# ---------------------------------------------------------------------------

def _write_synthetic_cvmh_slice(path: Path, *,
                                lon_axis: np.ndarray,
                                lat_axis: np.ndarray,
                                depth_m: float,
                                fields: dict,
                                spacing_deg: float = None) -> None:
    """Write a synthetic CVM-H ASCII slice file conforming to the CVM-H
    `# Lon,Lat,Vp(m/s),Vs(m/s),Density(kg/m^3)` format. The first
    sample row is the SW corner (lon[0], lat[0]); the last is the NE
    corner (lon[-1], lat[-1]). Values inside ``fields`` must be 2-D
    arrays of shape (Nlon, Nlat).

    ``spacing_deg`` defaults to the actual lon-step of ``lon_axis``;
    callers can override only if they want to test a header/spacing
    mismatch.
    """
    n_lon = lon_axis.size
    n_lat = lat_axis.size
    n_total = n_lon * n_lat
    if spacing_deg is None:
        spacing_deg = float(lon_axis[1] - lon_axis[0]) if n_lon > 1 else 0.0
    with open(path, "w") as f:
        f.write(f"# Title: synthetic slice at {depth_m} m\n")
        f.write("# CVM(abbr): synthetic\n")
        f.write(f"# Data_type: {','.join(fields)}\n")
        f.write(f"# Depth(m): {depth_m}\n")
        f.write(f"# Spacing(degree): {spacing_deg}\n")
        f.write(f"# Lon_pts: {n_lon}\n")
        f.write(f"# Lat_pts: {n_lat}\n")
        f.write(f"# Total_pts: {n_total}\n")
        f.write(f"# Lat1: {lat_axis[0]}\n")
        f.write(f"# Lon1: {lon_axis[0]}\n")
        f.write(f"# Lat2: {lat_axis[-1]}\n")
        f.write(f"# Lon2: {lon_axis[-1]}\n")
        cols = ["Lon", "Lat"]
        for n in fields:
            cols.append(f"{n}(synthetic)")
        f.write(f"# {','.join(cols)}\n")
        for i in range(n_lon):
            for j in range(n_lat):
                row = [f"{lon_axis[i]}", f"{lat_axis[j]}"]
                for arr in fields.values():
                    row.append(f"{float(arr[i, j])}")
                f.write(",".join(row) + "\n")


@pytest.fixture
def synthetic_cvmh_dir(tmp_path):
    """Two-slice synthetic CVM-H archive at depths 0, 1000 m."""
    lon_axis = np.array([-118.0, -117.5, -117.0], dtype=np.float64)
    lat_axis = np.array([33.5, 34.0], dtype=np.float64)
    # Vp grows with depth, Vs is constant, density grows with lon.
    fields_d0 = {
        "Vp":      4000.0 + 0.0 * np.outer(np.arange(3), np.ones(2)),
        "Vs":      2500.0 + 0.0 * np.outer(np.arange(3), np.ones(2)),
        "Density": 2400.0 + 50.0 * np.outer(np.arange(3), np.ones(2)),
    }
    fields_d1 = {
        "Vp":      4500.0 + 0.0 * np.outer(np.arange(3), np.ones(2)),
        "Vs":      2700.0 + 0.0 * np.outer(np.arange(3), np.ones(2)),
        "Density": 2450.0 + 50.0 * np.outer(np.arange(3), np.ones(2)),
    }
    _write_synthetic_cvmh_slice(
        tmp_path / "velocity_raw_0.0km.bp",
        lon_axis=lon_axis, lat_axis=lat_axis,
        depth_m=0.0, fields=fields_d0)
    _write_synthetic_cvmh_slice(
        tmp_path / "velocity_raw_0.1km.bp",
        lon_axis=lon_axis, lat_axis=lat_axis,
        depth_m=1000.0, fields=fields_d1)
    return tmp_path


# ---------------------------------------------------------------------------
# T-2-1 — orientation
# ---------------------------------------------------------------------------

def test_T_2_1_read_cvmh_orientation(synthetic_cvmh_dir):
    paths = sorted(synthetic_cvmh_dir.glob("velocity_raw_*.bp"))
    stack = read_cvmh_ascii(paths)
    assert stack.lon_grid.shape == (3,)
    assert stack.lat_grid.shape == (2,)
    assert stack.depths_m.shape == (2,)
    assert pytest.approx(stack.lon_grid[0]) == -118.0
    assert pytest.approx(stack.lon_grid[-1]) == -117.0
    assert pytest.approx(stack.lat_grid[0]) == 33.5
    assert pytest.approx(stack.lat_grid[-1]) == 34.0
    # Density grows with lon: 2400 at lon=-118, 2450 at lon=-117.5,
    # 2500 at lon=-117. After read_cvmh_ascii, density is sorted lon
    # ascending so density[0,:,k] = 2400 (lon=-118), density[-1,:,k]
    # = 2500 (lon=-117).
    assert np.allclose(stack.fields["density"][0, :, 0], 2400.0)
    assert np.allclose(stack.fields["density"][-1, :, 0], 2500.0)


# ---------------------------------------------------------------------------
# T-2-2 — filename trap (depth from header, not filename)
# ---------------------------------------------------------------------------

def test_T_2_2_filename_depth_label_is_ignored(tmp_path):
    """File named `velocity_raw_5km.bp` but with `# Depth(m): 50000`
    should yield ``depths_m[0] == 50000.0``. (The CVM-H archive in
    production uses misleading ``X.Ykm`` labels equal to depth_km/10.)"""
    lon_axis = np.array([-118.0, -117.5], dtype=np.float64)
    lat_axis = np.array([33.5, 34.0], dtype=np.float64)
    fields = {
        "Vp": 5000.0 * np.ones((2, 2)),
        "Vs": 3000.0 * np.ones((2, 2)),
        "Density": 2700.0 * np.ones((2, 2)),
    }
    p = tmp_path / "velocity_raw_5km.bp"  # filename label = "5km"
    _write_synthetic_cvmh_slice(
        p, lon_axis=lon_axis, lat_axis=lat_axis,
        depth_m=50000.0, fields=fields)
    stack = read_cvmh_ascii([p])
    assert stack.depths_m.size == 1
    assert pytest.approx(stack.depths_m[0]) == 50000.0


# ---------------------------------------------------------------------------
# T-2-3 — mismatched lon grids across slices is an error
# ---------------------------------------------------------------------------

def test_T_2_3_mismatched_grids(tmp_path):
    fields = {
        "Vp": 5000.0 * np.ones((2, 2)),
        "Vs": 3000.0 * np.ones((2, 2)),
        "Density": 2700.0 * np.ones((2, 2)),
    }
    _write_synthetic_cvmh_slice(
        tmp_path / "velocity_raw_a.bp",
        lon_axis=np.array([-118.0, -117.5]),
        lat_axis=np.array([33.5, 34.0]),
        depth_m=0.0, fields=fields)
    _write_synthetic_cvmh_slice(
        tmp_path / "velocity_raw_b.bp",
        lon_axis=np.array([-118.5, -118.0]),    # different
        lat_axis=np.array([33.5, 34.0]),
        depth_m=1000.0, fields=fields)
    paths = sorted(tmp_path.glob("velocity_raw_*.bp"))
    with pytest.raises(ValueError, match="lon_grid different"):
        read_cvmh_ascii(paths)


# ---------------------------------------------------------------------------
# T-2-4 — CRS round-trip
# ---------------------------------------------------------------------------

def test_T_2_4_crs_round_trip():
    rng = np.random.default_rng(42)
    lon = rng.uniform(-119.0, -116.0, size=1024)
    lat = rng.uniform(32.5, 35.0, size=1024)
    x, y = geographic_to_utm11n(lon, lat)
    lon2, lat2 = utm11n_to_geographic(x, y)
    np.testing.assert_allclose(lon2, lon, atol=1e-7)
    np.testing.assert_allclose(lat2, lat, atol=1e-7)


# ---------------------------------------------------------------------------
# T-2-5 — CRS calibration (anchor lon=-118.1677, lat=33.3428)
# ---------------------------------------------------------------------------

def test_T_2_5_crs_anchor_calibration():
    """Verifies a published-corner CVM-H anchor maps to the expected UTM
    location. Cross-checked against pyproj's published EPSG transform.
    """
    lon, lat = np.array([-118.1677]), np.array([33.3428])
    x, y = geographic_to_utm11n(lon, lat)
    # Reference values produced via pyproj 3.7.2 + EPSG database; tolerance
    # 0.1 m to detect any silent transformer-config breakage.
    expected_x = 391339.10822640534
    expected_y = 3689899.487805032
    assert abs(float(x[0]) - expected_x) < 0.1
    assert abs(float(y[0]) - expected_y) < 0.1


# ---------------------------------------------------------------------------
# T-2-6 / T-2-7 — bbox containment guards
# ---------------------------------------------------------------------------

def test_T_2_6_bbox_contains_happy():
    x = np.array([0.0, 5.0, 10.0])
    y = np.array([0.0, 10.0])
    z = np.array([-100.0, 0.0])
    mesh = {"xmin": 1.0, "xmax": 9.0, "ymin": 0.0, "ymax": 10.0,
            "zmin": -50.0, "zmax": 0.0}
    assert grid_contains_bbox(x, y, z, mesh)
    assert_grid_contains_mesh(x, y, z, mesh)   # must not raise


def test_T_2_7_bbox_contains_sad():
    x = np.array([0.0, 5.0, 10.0])
    y = np.array([0.0, 10.0])
    z = np.array([-100.0, 0.0])
    mesh = {"xmin": -1.0, "xmax": 9.0, "ymin": 0.0, "ymax": 10.0,
            "zmin": -50.0, "zmax": 0.0}   # x violates
    assert not grid_contains_bbox(x, y, z, mesh)
    with pytest.raises(BBoxContainmentError) as excinfo:
        assert_grid_contains_mesh(x, y, z, mesh)
    msg = str(excinfo.value)
    # Both bbox tables must appear in the message.
    assert "data grid bbox" in msg
    assert "mesh bbox" in msg


# ---------------------------------------------------------------------------
# T-2-8 — mesh_msh_bbox round-trip via meshio
# ---------------------------------------------------------------------------

def test_T_2_8_mesh_msh_bbox(tmp_path):
    """Build a small synthetic gmsh-format mesh and confirm
    mesh_msh_bbox reads its bbox correctly."""
    import meshio
    points = np.array([
        [0.0, 0.0, 0.0],
        [1000.0, 0.0, 0.0],
        [1000.0, 500.0, 0.0],
        [0.0, 500.0, 0.0],
        [500.0, 250.0, -100.0],
    ], dtype=np.float64)
    # One tetrahedron; topology is irrelevant to the bbox test.
    cells = [("tetra", np.array([[0, 1, 2, 4]], dtype=np.int64))]
    p = tmp_path / "tiny.msh"
    meshio.write_points_cells(str(p), points, cells, file_format="gmsh")
    bbox = mesh_msh_bbox(p)
    assert pytest.approx(bbox["xmin"]) == 0.0
    assert pytest.approx(bbox["xmax"]) == 1000.0
    assert pytest.approx(bbox["ymin"]) == 0.0
    assert pytest.approx(bbox["ymax"]) == 500.0
    assert pytest.approx(bbox["zmin"]) == -100.0
    assert pytest.approx(bbox["zmax"]) == 0.0


# ---------------------------------------------------------------------------
# T-2-9 — sidecar I/O round-trip
# ---------------------------------------------------------------------------

def test_T_2_9_sidecar_round_trip(tmp_path):
    import h5py
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 2.0, 3.0])
    z = np.array([-3.0, -2.0, -1.0, 0.0])
    Vp = 4000.0 * np.ones((4, 4, 4))
    fields = {"Vp": Vp}
    field_bounds = {"Vp": (2000.0, 9000.0, "m/s")}
    out = tmp_path / "round_trip.h5"
    write_sidecar(out, x, y, z, fields, attrs={}, field_bounds=field_bounds)
    attrs = read_sidecar_attrs(out)
    assert attrs["schema_version"] == SCHEMA_VERSION
    assert attrs["crs"] == "EPSG:32611"
    assert attrs["units"] == "m"
    assert attrs["z_positive"] == "elevation"
    with h5py.File(out, "r") as h5:
        np.testing.assert_array_equal(h5["grid/x"][...], x)
        np.testing.assert_array_equal(h5["grid/y"][...], y)
        np.testing.assert_array_equal(h5["grid/z"][...], z)
        np.testing.assert_array_equal(h5["fields/Vp"][...], Vp)
        ds_attrs = dict(h5["fields/Vp"].attrs)
        # h5py stores strings as bytes by default; coerce to str.
        units = ds_attrs["units"]
        units = units.decode() if isinstance(units, bytes) else units
        assert units == "m/s"
        assert pytest.approx(float(ds_attrs["min_value"])) == 2000.0
        assert pytest.approx(float(ds_attrs["max_value"])) == 9000.0


# ---------------------------------------------------------------------------
# T-2-10 — NaN guard
# ---------------------------------------------------------------------------

def test_T_2_10_writer_rejects_nan(tmp_path):
    x = np.array([0.0, 1.0])
    y = np.array([0.0, 1.0])
    z = np.array([-1.0, 0.0])
    Vp = 4000.0 * np.ones((2, 2, 2))
    Vp[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="contains 1 NaN cell"):
        write_sidecar(
            tmp_path / "with_nan.h5", x, y, z,
            fields={"Vp": Vp},
            attrs={},
            field_bounds={"Vp": (2000.0, 9000.0, "m/s")})


# ---------------------------------------------------------------------------
# T-2-11 — per-field min guard (Vs floor)
# ---------------------------------------------------------------------------

def test_T_2_11_writer_rejects_subfloor_vs(tmp_path):
    x = np.array([0.0, 1.0])
    y = np.array([0.0, 1.0])
    z = np.array([-1.0, 0.0])
    Vs = 2500.0 * np.ones((2, 2, 2))
    Vs[1, 1, 1] = 100.0   # sub-seafloor / water-saturated
    with pytest.raises(ValueError) as excinfo:
        write_sidecar(
            tmp_path / "with_seafloor.h5", x, y, z,
            fields={"Vs": Vs},
            attrs={},
            field_bounds={"Vs": (1500.0, 5000.0, "m/s")})
    assert "1 cell(s) below min=1500" in str(excinfo.value)


# ---------------------------------------------------------------------------
# T-2-12 — full pipeline smoke (synthetic, deterministic)
# ---------------------------------------------------------------------------

def test_T_2_12_build_pipeline_synthetic(tmp_path, synthetic_cvmh_dir):
    """End-to-end: read synthetic CVM-H slices, build sidecar with the
    real driver, verify schema attrs and per-field min/max ranges.
    """
    import h5py
    from data_projection.build_velocity_cvmh import build_velocity_sidecar

    out = tmp_path / "synthetic.h5"
    # The synthetic fixture uses Vs = 2500..2700, well within [1500, 5000].
    build_velocity_sidecar(
        raw_dir=synthetic_cvmh_dir,
        out_path=out,
        grid_dx=10000.0,    # ~0.1 deg in metres
        vp_min_mps=2000, vp_max_mps=9000,
        vs_min_mps=1500, vs_max_mps=5000,
        rho_min_kgm3=2000, rho_max_kgm3=3500,
        mesh_msh=None,
        verbose=False)
    with h5py.File(out, "r") as h5:
        assert h5.attrs["schema_version"] == SCHEMA_VERSION
        for name in ("Vp", "Vs", "density"):
            assert name in h5["fields"], f"field '{name}' missing"
            arr = h5[f"fields/{name}"][...]
            v_min = float(h5[f"fields/{name}"].attrs["min_value"])
            v_max = float(h5[f"fields/{name}"].attrs["max_value"])
            assert (arr >= v_min).all()
            assert (arr <= v_max).all()
            assert not np.isnan(arr).any()


# ---------------------------------------------------------------------------
# Bonus: build_velocity_sidecar G-3 (mesh containment) refuses too-small
# ---------------------------------------------------------------------------

def test_build_pipeline_refuses_when_mesh_too_big(tmp_path,
                                                   synthetic_cvmh_dir):
    """G-3: writer-side mesh-bbox containment guard fires when the mesh
    extends beyond the data grid."""
    import meshio
    from data_projection.build_velocity_cvmh import build_velocity_sidecar

    # Build a mesh that spans far beyond the synthetic data hull
    # (which is roughly UTM x ~ [40 km, 100 km], y ~ [3.7 M, 3.8 M]).
    points = np.array([
        [-1_000_000.0, -1_000_000.0,  -100.0],
        [ 2_000_000.0, -1_000_000.0,  -100.0],
        [ 2_000_000.0,  5_000_000.0,  -100.0],
        [-1_000_000.0,  5_000_000.0,    50.0],
    ], dtype=np.float64)
    cells = [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))]
    msh = tmp_path / "huge.msh"
    meshio.write_points_cells(str(msh), points, cells, file_format="gmsh")

    out = tmp_path / "should_not_be_written.h5"
    with pytest.raises(BBoxContainmentError):
        build_velocity_sidecar(
            raw_dir=synthetic_cvmh_dir, out_path=out, grid_dx=10000.0,
            vp_min_mps=2000, vp_max_mps=9000,
            vs_min_mps=1500, vs_max_mps=5000,
            rho_min_kgm3=2000, rho_max_kgm3=3500,
            mesh_msh=msh, verbose=False)
    assert not out.exists()


# ---------------------------------------------------------------------------
# R-006 — extend_z_top_m clones the surface slice
# ---------------------------------------------------------------------------

def test_R006_extend_z_top_clones_surface_slice(tmp_path,
                                                 synthetic_cvmh_dir):
    """With extend_z_top_m=500, the sidecar's z axis grows by exactly
    one entry, the entry equals z[-2] + 500, and every field's last
    z slab equals its second-to-last (clone of surface)."""
    import h5py
    from data_projection.build_velocity_cvmh import build_velocity_sidecar
    out = tmp_path / "ext.h5"
    build_velocity_sidecar(
        raw_dir=synthetic_cvmh_dir, out_path=out, grid_dx=10000.0,
        vp_min_mps=1000, vp_max_mps=9000,
        vs_min_mps=100,  vs_max_mps=5000,
        rho_min_kgm3=1000, rho_max_kgm3=3500,
        extend_z_top_m=500.0, mesh_msh=None, verbose=False)
    with h5py.File(out, "r") as h5:
        z = h5["grid/z"][...]
        assert z.size >= 2
        assert pytest.approx(z[-1] - z[-2]) == 500.0
        for name in ("Vp", "Vs", "density"):
            arr = h5[f"fields/{name}"][...]
            assert arr.shape[2] == z.size
            np.testing.assert_array_equal(
                arr[:, :, -1], arr[:, :, -2],
                err_msg=f"top slab of '{name}' must clone surface")


# ---------------------------------------------------------------------------
# R-007 — extra source field raises
# ---------------------------------------------------------------------------

def test_R007_extra_source_field_raises(tmp_path):
    """If the CVM-H source has a fourth column not in the canonical
    map (e.g. Qp), the build must fail loudly, not silently drop it."""
    from data_projection.build_velocity_cvmh import build_velocity_sidecar
    lon_axis = np.array([-118.0, -117.5, -117.0], dtype=np.float64)
    lat_axis = np.array([33.5, 34.0], dtype=np.float64)
    fields = {
        "Vp":      4000.0 * np.ones((3, 2)),
        "Vs":      2500.0 * np.ones((3, 2)),
        "Density": 2400.0 * np.ones((3, 2)),
        "Qp":      150.0  * np.ones((3, 2)),
    }
    p = tmp_path / "velocity_raw_0.0km.bp"
    _write_synthetic_cvmh_slice(
        p, lon_axis=lon_axis, lat_axis=lat_axis,
        depth_m=0.0, fields=fields)
    with pytest.raises(ValueError) as exc:
        build_velocity_sidecar(
            raw_dir=tmp_path, out_path=tmp_path / "out.h5",
            grid_dx=10000.0,
            vp_min_mps=1000, vp_max_mps=9000,
            vs_min_mps=100,  vs_max_mps=5000,
            rho_min_kgm3=1000, rho_max_kgm3=3500,
            mesh_msh=None, verbose=False)
    assert "qp" in str(exc.value).lower()


# ---------------------------------------------------------------------------
# R-002 — paraview companion file written next to the .h5 sidecar
# ---------------------------------------------------------------------------

def test_R002_paraview_companion(tmp_path, synthetic_cvmh_dir):
    """With paraview_export=True, a .vti or .vtr companion is written
    next to the .h5 and contains the three named field arrays."""
    from data_projection.build_velocity_cvmh import build_velocity_sidecar
    out_h5 = tmp_path / "with_pv.h5"
    build_velocity_sidecar(
        raw_dir=synthetic_cvmh_dir, out_path=out_h5,
        grid_dx=10000.0,
        vp_min_mps=1000, vp_max_mps=9000,
        vs_min_mps=100,  vs_max_mps=5000,
        rho_min_kgm3=1000, rho_max_kgm3=3500,
        mesh_msh=None, paraview_export=True, verbose=False)
    companion = out_h5.with_suffix(".vti")
    if not companion.exists():
        companion = out_h5.with_suffix(".vtr")
    assert companion.exists(), "expected companion .vti or .vtr next to .h5"
    text = companion.read_text()
    for token in ("Vp", "Vs", "density",
                  "AppendedData", "VTKFile"):
        assert token in text, f"companion file lacks '{token}'"
