"""pytest suite for convert_cvm_to_asagi.py (pure numpy + h5py + netCDF4).

Run:  python3 -m pytest -q test_convert_cvm_to_asagi.py
"""

import h5py
import numpy as np
import pytest

from convert_cvm_to_asagi import (
    main,
    read_asagi_netcdf,
    read_sidecar,
    resample_z,
    trilinear_sample,
    validate_axes,
    velocities_to_moduli,
    write_asagi_netcdf,
)

FLOAT32_EPS = np.finfo(np.float32).eps


def write_synthetic_sidecar(path, x, y, z, vp, vs, rho):
    """Minimal schema-v1 sidecar with canonical root attrs."""
    with h5py.File(path, "w") as h5:
        h5.attrs["schema_version"] = "data_projection_v1"
        h5.attrs["crs"] = "EPSG:32611"
        h5.attrs["units"] = "m"
        h5.attrs["z_positive"] = "elevation"
        h5.attrs["created_at"] = "2026-06-09T00:00:00Z"
        h5.attrs["mesh_tag"] = "synthetic"
        grid = h5.create_group("grid")
        grid.create_dataset("x", data=np.asarray(x, dtype=np.float64))
        grid.create_dataset("y", data=np.asarray(y, dtype=np.float64))
        grid.create_dataset("z", data=np.asarray(z, dtype=np.float64))
        flds = h5.create_group("fields")
        flds.create_dataset("Vp", data=vp)
        flds.create_dataset("Vs", data=vs)
        flds.create_dataset("density", data=rho)


def test_validate_axes_rejects_nonuniform_x():
    x = np.array([0.0, 1.0, 2.5, 3.0])  # non-uniform
    y = np.array([0.0, 1.0, 2.0])
    z = np.array([0.0, 1.0, 5.0])  # non-uniform z is allowed
    with pytest.raises(ValueError, match="axis 'x' is not uniform"):
        validate_axes(x, y, z)
    # uniform x/y with non-uniform z passes
    validate_axes(np.array([0.0, 1.0, 2.0]), y, z)


def test_validate_axes_rejects_decreasing():
    good = np.array([0.0, 1.0, 2.0])
    bad = np.array([0.0, 2.0, 1.0])
    with pytest.raises(ValueError, match="not strictly increasing"):
        validate_axes(good, good, bad)


def test_moduli_constant_model_parity():
    # The verified constant medium: mu = lambda = 32 GPa, rho = 2670.
    # Exact velocities: Vs = sqrt(mu/rho), Vp = sqrt((lambda + 2 mu)/rho).
    mu_ref, lam_ref, rho_ref = 3.2e10, 3.2e10, 2670.0
    vs = np.sqrt(mu_ref / rho_ref)            # ~3461.94 m/s
    vp = np.sqrt((lam_ref + 2 * mu_ref) / rho_ref)  # ~5996.25 m/s
    rho, mu, lam = velocities_to_moduli(
        np.full((2, 2, 2), vp), np.full((2, 2, 2), vs),
        np.full((2, 2, 2), rho_ref))
    assert np.allclose(mu, mu_ref, rtol=1e-12)
    assert np.allclose(lam, lam_ref, rtol=1e-12)
    assert np.allclose(rho, rho_ref, rtol=1e-12)


def test_moduli_rejects_lambda_nonpositive():
    # Vp < sqrt(2) * Vs  =>  lambda <= 0  =>  hard error.
    vp = np.full((2, 2, 2), 1000.0)
    vs = np.full((2, 2, 2), 800.0)
    rho = np.full((2, 2, 2), 2000.0)
    with pytest.raises(ValueError, match="lambda <= 0"):
        velocities_to_moduli(vp, vs, rho)


def test_resample_z_exact_at_source_levels():
    nx, ny = 3, 2
    z_src = np.array([-100.0, -50.0, -20.0, 0.0])
    rng = np.random.default_rng(7)
    field = rng.uniform(1.0, 2.0, (nx, ny, z_src.size))
    z_out = np.array([-100.0, -50.0, -20.0, 0.0])  # all source levels
    out = resample_z(field, z_src, z_out)
    assert np.allclose(out, field, rtol=0, atol=1e-15)
    # subset of source levels, plus a midpoint
    out2 = resample_z(field, z_src, np.array([-50.0, -35.0, 0.0]))
    assert np.allclose(out2[:, :, 0], field[:, :, 1], atol=1e-15)
    assert np.allclose(out2[:, :, 2], field[:, :, 3], atol=1e-15)
    mid = 0.5 * (field[:, :, 1] + field[:, :, 2])
    assert np.allclose(out2[:, :, 1], mid, atol=1e-15)


def test_resample_z_rejects_extrapolation():
    field = np.ones((2, 2, 3))
    z_src = np.array([-10.0, -5.0, 0.0])
    with pytest.raises(ValueError, match="extrapolation"):
        resample_z(field, z_src, np.array([-11.0, 0.0]))


def test_netcdf_roundtrip_affine_field(tmp_path):
    # Affine fields are reproduced exactly by trilinear interpolation, so
    # the only error is float32 storage rounding (<= ~10 eps relative).
    x = np.linspace(0.0, 4000.0, 5)
    y = np.linspace(0.0, 3000.0, 4)
    z = np.linspace(-6000.0, 0.0, 7)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    rho = 2000.0 + 0.01 * xx + 0.02 * yy - 0.005 * zz
    mu = 3.0e10 + 1e5 * xx - 2e5 * yy + 3e5 * zz
    lam = 2.5e10 - 1e5 * xx + 1e5 * yy - 2e5 * zz
    nc = tmp_path / "affine.nc"
    write_asagi_netcdf(nc, x, y, z, rho, mu, lam, "float32",
                       {"converter": "test"})

    xg, yg, zg, fields = read_asagi_netcdf(nc)
    assert np.array_equal(xg, x) and np.array_equal(yg, y)
    assert np.array_equal(zg, z)

    rng = np.random.default_rng(11)
    pts = np.column_stack([rng.uniform(x[0], x[-1], 200),
                           rng.uniform(y[0], y[-1], 200),
                           rng.uniform(z[0], z[-1], 200)])
    analytic = {
        "rho": 2000.0 + 0.01 * pts[:, 0] + 0.02 * pts[:, 1]
               - 0.005 * pts[:, 2],
        "mu": 3.0e10 + 1e5 * pts[:, 0] - 2e5 * pts[:, 1]
              + 3e5 * pts[:, 2],
        "lambda": 2.5e10 - 1e5 * pts[:, 0] + 1e5 * pts[:, 1]
                  - 2e5 * pts[:, 2],
    }
    for name in ("rho", "mu", "lambda"):
        got = trilinear_sample(xg, yg, zg, fields[name], pts)
        rel = np.abs(got - analytic[name]) / np.abs(analytic[name])
        assert rel.max() <= 10 * FLOAT32_EPS, (
            f"{name}: max rel err {rel.max():.3e}")


def test_netcdf_layout_compound_zyx(tmp_path):
    # The on-disk layout must be data(z, y, x) with compound members
    # rho/mu/lambda — verified against raw netCDF4 introspection, the
    # same layout as the SeisSol Training reference files.
    import netCDF4

    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 2.0])
    z = np.array([-1.0, 0.0])
    shape = (x.size, y.size, z.size)
    # Encode the index in the value: val = ix*100 + iy*10 + iz.
    ix, iy, iz = np.meshgrid(np.arange(4), np.arange(3), np.arange(2),
                             indexing="ij")
    marker = (ix * 100 + iy * 10 + iz).astype(np.float64)
    ones = np.ones(shape)
    nc = tmp_path / "layout.nc"
    write_asagi_netcdf(nc, x, y, z, marker, ones, ones, "float64", {})
    with netCDF4.Dataset(nc) as ds:
        var = ds.variables["data"]
        assert var.dimensions == ("z", "y", "x")
        assert set(var.dtype.names) == {"rho", "mu", "lambda"}
        raw = var[:]
        # raw is (nz, ny, nx): entry [iz, iy, ix] must be the marker.
        assert raw["rho"][1, 2, 3] == pytest.approx(3 * 100 + 2 * 10 + 1)
        assert raw["rho"][0, 0, 0] == pytest.approx(0.0)


def test_read_sidecar_rejects_missing_field(tmp_path):
    path = tmp_path / "broken.h5"
    x = np.array([0.0, 1.0])
    write_synthetic_sidecar(path, x, x, x,
                            np.full((2, 2, 2), 6000.0),
                            np.full((2, 2, 2), 3500.0),
                            np.full((2, 2, 2), 2700.0))
    with h5py.File(path, "a") as h5:
        del h5["fields/Vs"]
    with pytest.raises(ValueError, match="fields/Vs"):
        read_sidecar(path)


def test_read_sidecar_rejects_wrong_crs(tmp_path):
    path = tmp_path / "crs.h5"
    x = np.array([0.0, 1.0])
    write_synthetic_sidecar(path, x, x, x,
                            np.full((2, 2, 2), 6000.0),
                            np.full((2, 2, 2), 3500.0),
                            np.full((2, 2, 2), 2700.0))
    with h5py.File(path, "a") as h5:
        h5.attrs["crs"] = "EPSG:4326"
    with pytest.raises(ValueError, match="crs"):
        read_sidecar(path)


def test_main_end_to_end_constant_medium(tmp_path, monkeypatch, capsys):
    # Full CLI run on a synthetic constant-medium sidecar covering the
    # mesh bbox: the self-check must PASS (exit code 0) and the file must
    # reproduce mu = lambda = 32 GPa.
    import convert_cvm_to_asagi as conv

    mu_ref, lam_ref, rho_ref = 3.2e10, 3.2e10, 2670.0
    vs = np.sqrt(mu_ref / rho_ref)
    vp = np.sqrt((lam_ref + 2 * mu_ref) / rho_ref)
    x = np.linspace(0.0, 10000.0, 6)
    y = np.linspace(0.0, 8000.0, 5)
    z = np.array([-5000.0, -3000.0, -2000.0, -500.0, -250.0, 0.0])
    shape = (x.size, y.size, z.size)
    sidecar = tmp_path / "const.h5"
    write_synthetic_sidecar(sidecar, x, y, z,
                            np.full(shape, vp), np.full(shape, vs),
                            np.full(shape, rho_ref))
    # Point the self-check bbox inside the synthetic grid.
    monkeypatch.setattr(conv, "MESH_BBOX",
                        ((500.0, 9500.0), (500.0, 7500.0),
                         (-4500.0, -100.0)))
    out = tmp_path / "const.nc"
    rc = main(["--sidecar", str(sidecar), "--out", str(out),
               "--z-min", "-5000", "--z-max", "0", "--dz", "250"])
    assert rc == 0
    assert "PASS" in capsys.readouterr().out
    _, _, zg, fields = read_asagi_netcdf(out)
    assert zg.size == 21  # -5000 .. 0 @ 250
    assert np.allclose(fields["mu"], mu_ref, rtol=1e-6)
    assert np.allclose(fields["lambda"], lam_ref, rtol=1e-6)


def test_main_warns_on_missed_source_level(tmp_path, monkeypatch, capsys):
    # Source has a level at -750 that a dz=500 grid from -2000 cannot
    # represent: the converter must warn (model smoothed there).
    import convert_cvm_to_asagi as conv

    x = np.linspace(0.0, 3000.0, 4)
    z = np.array([-2000.0, -1000.0, -750.0, -500.0, 0.0])
    shape = (4, 4, 5)
    sidecar = tmp_path / "kink.h5"
    write_synthetic_sidecar(sidecar, x, x, z,
                            np.full(shape, 6000.0),
                            np.full(shape, 3500.0),
                            np.full(shape, 2700.0))
    monkeypatch.setattr(conv, "MESH_BBOX",
                        ((100.0, 2900.0), (100.0, 2900.0),
                         (-1900.0, -100.0)))
    rc = main(["--sidecar", str(sidecar), "--out", str(tmp_path / "k.nc"),
               "--z-min", "-2000", "--z-max", "0", "--dz", "500"])
    out = capsys.readouterr().out
    assert "WARNING" in out and "-750" in out
    assert rc == 0  # constant medium: smoothing is exact anyway


def test_main_rejects_z_range_outside_source(tmp_path):
    x = np.array([0.0, 1000.0, 2000.0])
    z = np.array([-1000.0, -500.0, 0.0])
    shape = (3, 3, 3)
    sidecar = tmp_path / "small.h5"
    write_synthetic_sidecar(sidecar, x, x, z,
                            np.full(shape, 6000.0),
                            np.full(shape, 3500.0),
                            np.full(shape, 2700.0))
    with pytest.raises(SystemExit):
        main(["--sidecar", str(sidecar), "--out", str(tmp_path / "o.nc"),
              "--z-min", "-2000", "--z-max", "0"])


def test_vs_floor_lifts_vp(tmp_path, monkeypatch, capsys):
    import convert_cvm_to_asagi as conv

    # Medium with a slow layer: Vs = 200 m/s. With --vs-floor 500 the
    # written mu must reflect Vs >= 500 and lambda must stay positive.
    x = np.linspace(0.0, 4000.0, 4)
    z = np.array([-3000.0, -1000.0, 0.0])
    shape = (4, 4, 3)
    vs = np.full(shape, 200.0)
    vp = np.full(shape, 1000.0)
    rho = np.full(shape, 2000.0)
    sidecar = tmp_path / "slow.h5"
    write_synthetic_sidecar(sidecar, x, x, z, vp, vs, rho)
    monkeypatch.setattr(conv, "MESH_BBOX",
                        ((100.0, 3900.0), (100.0, 3900.0),
                         (-2900.0, -100.0)))
    out = tmp_path / "slow.nc"
    rc = main(["--sidecar", str(sidecar), "--out", str(out),
               "--z-min", "-3000", "--z-max", "0", "--dz", "500",
               "--vs-floor", "500"])
    assert rc == 0
    _, _, _, fields = read_asagi_netcdf(out)
    mu_expected = 2000.0 * 500.0**2
    assert np.allclose(fields["mu"], mu_expected, rtol=1e-6)
    assert (fields["lambda"] > 0).all()
