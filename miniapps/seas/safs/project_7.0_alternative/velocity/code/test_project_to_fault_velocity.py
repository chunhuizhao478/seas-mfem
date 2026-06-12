#!/usr/bin/env python3
"""Unit tests for project_to_fault_velocity.py.

Run:
    cd project_7.0_alternative/velocity/code && \
        pytest -q test_project_to_fault_velocity.py
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import meshio
import numpy as np
import pytest

from project_to_fault_velocity import (
    DERIVED_FIELDS,
    SIDECAR_FIELDS,
    VelocityGrid,
    count_out_of_bounds,
    load_velocity_sidecar,
    main,
    sample_fields_at_points,
    trilinear_sample,
    write_fault_velocity_vtu,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _linear_field(gx, gy, gz, a, b, c, d):
    """F[i,j,k] = a + b·x + c·y + d·z on the tensor grid."""
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    return a + b * X + c * Y + d * Z


def _make_sidecar(path: Path, gx, gy, gz, fields: dict) -> Path:
    with h5py.File(str(path), "w") as h5:
        g = h5.create_group("grid")
        g.create_dataset("x", data=np.asarray(gx, dtype=np.float64))
        g.create_dataset("y", data=np.asarray(gy, dtype=np.float64))
        g.create_dataset("z", data=np.asarray(gz, dtype=np.float64))
        f = h5.create_group("fields")
        for name, arr in fields.items():
            f.create_dataset(name, data=np.asarray(arr, dtype=np.float64))
        h5.attrs["schema_version"] = "data_projection_v1"
        h5.attrs["crs"] = "EPSG:32611"
    return path


# Non-uniform axes (z deliberately non-uniform like the real sidecar).
GX = np.array([0.0, 1000.0, 2000.0, 3500.0])
GY = np.array([0.0, 800.0, 1600.0])
GZ = np.array([-5000.0, -1000.0, -500.0, 0.0, 100.0])


@pytest.fixture
def sidecar(tmp_path):
    """Sidecar with linear Vp/Vs/density fields (trilinear-exact)."""
    fields = {
        "Vp": _linear_field(GX, GY, GZ, 5000.0, 0.1, 0.05, -0.2),
        "Vs": _linear_field(GX, GY, GZ, 3000.0, 0.05, 0.02, -0.1),
        "density": _linear_field(GX, GY, GZ, 2500.0, 0.01, 0.0, -0.01),
    }
    return _make_sidecar(tmp_path / "vel.h5", GX, GY, GZ, fields)


def _make_msh(path: Path) -> Path:
    """Tiny Gmsh v2.2 mesh: one vertical fault triangle pair
    (physical surface 'fault', tag 101) + one tet ('rock', tag 1),
    inside the GX×GY×GZ box."""
    points = np.array(
        [
            # fault quad (vertical plane y = 400)
            [500.0, 400.0, -2000.0],
            [1500.0, 400.0, -2000.0],
            [1500.0, 400.0, -200.0],
            [500.0, 400.0, -200.0],
            # extra vertex for the tet
            [1000.0, 1200.0, -1000.0],
        ]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("tetra", np.array([[0, 1, 2, 4]])),
    ]
    cell_data = {
        "gmsh:physical": [
            np.array([101, 101]),
            np.array([1]),
        ],
        "gmsh:geometrical": [
            np.array([1, 1]),
            np.array([1]),
        ],
    }
    mesh = meshio.Mesh(
        points=points, cells=cells, cell_data=cell_data,
        field_data={"fault": np.array([101, 2]),
                    "rock": np.array([1, 3])},
    )
    meshio.write(str(path), mesh, file_format="gmsh22", binary=False)
    return path


# ----------------------------------------------------------------------
# trilinear_sample
# ----------------------------------------------------------------------


class TestTrilinearSample:
    def test_exact_on_linear_field(self):
        F = _linear_field(GX, GY, GZ, 10.0, 0.3, -0.7, 0.11)
        rng = np.random.default_rng(42)
        P = np.column_stack(
            [
                rng.uniform(GX[0], GX[-1], 200),
                rng.uniform(GY[0], GY[-1], 200),
                rng.uniform(GZ[0], GZ[-1], 200),
            ]
        )
        got = trilinear_sample(GX, GY, GZ, F, P)
        want = 10.0 + 0.3 * P[:, 0] - 0.7 * P[:, 1] + 0.11 * P[:, 2]
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-9)

    def test_grid_nodes_reproduced(self):
        F = _linear_field(GX, GY, GZ, 1.0, 2.0, 3.0, 4.0)
        P = np.array([[GX[i], GY[j], GZ[k]]
                      for i in range(len(GX))
                      for j in range(len(GY))
                      for k in range(len(GZ))])
        got = trilinear_sample(GX, GY, GZ, F, P)
        want = F.ravel()
        np.testing.assert_allclose(got, want, rtol=1e-12)

    def test_out_of_bounds_clamped_to_boundary(self):
        F = _linear_field(GX, GY, GZ, 0.0, 1.0, 0.0, 0.0)
        P = np.array(
            [
                [GX[0] - 5000.0, GY[1], GZ[1]],   # below x range
                [GX[-1] + 5000.0, GY[1], GZ[1]],  # above x range
            ]
        )
        got = trilinear_sample(GX, GY, GZ, F, P)
        np.testing.assert_allclose(got, [GX[0], GX[-1]], rtol=1e-12)

    def test_nan_point_yields_nan(self):
        F = _linear_field(GX, GY, GZ, 1.0, 0.0, 0.0, 0.0)
        P = np.array([[np.nan, GY[1], GZ[1]], [GX[1], GY[1], GZ[1]]])
        got = trilinear_sample(GX, GY, GZ, F, P)
        assert np.isnan(got[0])
        assert np.isfinite(got[1])

    def test_shape_validation(self):
        F = _linear_field(GX, GY, GZ, 1.0, 0.0, 0.0, 0.0)
        with pytest.raises(ValueError, match=r"P must be \(N, 3\)"):
            trilinear_sample(GX, GY, GZ, F, np.zeros((4, 2)))
        with pytest.raises(ValueError, match="F shape"):
            trilinear_sample(GX, GY, GZ, F[:-1], np.zeros((1, 3)))


class TestCountOutOfBounds:
    def test_counts_only_finite_outside(self):
        P = np.array(
            [
                [GX[1], GY[1], GZ[1]],            # inside
                [GX[0] - 1.0, GY[1], GZ[1]],      # outside x
                [GX[1], GY[-1] + 1.0, GZ[1]],     # outside y
                [np.nan, GY[1], GZ[1]],           # NaN — not counted
            ]
        )
        assert count_out_of_bounds(GX, GY, GZ, P) == 2


# ----------------------------------------------------------------------
# load_velocity_sidecar
# ----------------------------------------------------------------------


class TestLoadVelocitySidecar:
    def test_roundtrip(self, sidecar):
        grid = load_velocity_sidecar(sidecar)
        np.testing.assert_array_equal(grid.x, GX)
        np.testing.assert_array_equal(grid.y, GY)
        np.testing.assert_array_equal(grid.z, GZ)
        assert set(grid.fields) == set(SIDECAR_FIELDS)
        assert grid.fields["Vp"].shape == (len(GX), len(GY), len(GZ))
        assert grid.attrs["schema_version"] == "data_projection_v1"

    def test_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            load_velocity_sidecar(tmp_path / "nope.h5")

    def test_missing_field(self, tmp_path):
        p = _make_sidecar(
            tmp_path / "bad.h5", GX, GY, GZ,
            {"Vp": np.zeros((len(GX), len(GY), len(GZ)))},
        )
        with pytest.raises(KeyError, match="fields/Vs"):
            load_velocity_sidecar(p)

    def test_non_monotone_axis(self, tmp_path):
        gx_bad = GX.copy()
        gx_bad[1], gx_bad[2] = gx_bad[2], gx_bad[1]
        p = _make_sidecar(
            tmp_path / "bad.h5", gx_bad, GY, GZ,
            {n: np.zeros((len(GX), len(GY), len(GZ)))
             for n in SIDECAR_FIELDS},
        )
        with pytest.raises(ValueError, match="strictly increasing"):
            load_velocity_sidecar(p)

    def test_shape_mismatch(self, tmp_path):
        p = _make_sidecar(
            tmp_path / "bad.h5", GX, GY, GZ,
            {
                "Vp": np.zeros((len(GX), len(GY), len(GZ) - 1)),
                "Vs": np.zeros((len(GX), len(GY), len(GZ))),
                "density": np.zeros((len(GX), len(GY), len(GZ))),
            },
        )
        with pytest.raises(ValueError, match="grid shape"):
            load_velocity_sidecar(p)


# ----------------------------------------------------------------------
# sample_fields_at_points (derived fields)
# ----------------------------------------------------------------------


class TestSampleFields:
    def test_derived_mu_and_ratio(self, sidecar):
        grid = load_velocity_sidecar(sidecar)
        P = np.array([[1200.0, 700.0, -900.0]])
        out = sample_fields_at_points(grid, P)
        vp = 5000.0 + 0.1 * 1200.0 + 0.05 * 700.0 - 0.2 * (-900.0)
        vs = 3000.0 + 0.05 * 1200.0 + 0.02 * 700.0 - 0.1 * (-900.0)
        rho = 2500.0 + 0.01 * 1200.0 - 0.01 * (-900.0)
        np.testing.assert_allclose(out["Vp_m_per_s"], [vp], rtol=1e-12)
        np.testing.assert_allclose(out["Vs_m_per_s"], [vs], rtol=1e-12)
        np.testing.assert_allclose(
            out["density_kg_per_m3"], [rho], rtol=1e-12
        )
        np.testing.assert_allclose(
            out["mu_GPa"], [rho * vs ** 2 * 1e-9], rtol=1e-12
        )
        np.testing.assert_allclose(
            out["vp_vs_ratio"], [vp / vs], rtol=1e-12
        )

    def test_ratio_nan_on_nonpositive_vs(self):
        nx, ny, nz = len(GX), len(GY), len(GZ)
        grid = VelocityGrid(
            x=GX, y=GY, z=GZ,
            fields={
                "Vp": np.full((nx, ny, nz), 4000.0),
                "Vs": np.zeros((nx, ny, nz)),
                "density": np.full((nx, ny, nz), 2700.0),
            },
        )
        out = sample_fields_at_points(
            grid, np.array([[GX[1], GY[1], GZ[1]]])
        )
        assert np.isnan(out["vp_vs_ratio"][0])
        assert out["mu_GPa"][0] == 0.0


# ----------------------------------------------------------------------
# VTU writer + end-to-end CLI
# ----------------------------------------------------------------------


class TestWriteFaultVtu:
    def test_shape_validation(self, tmp_path):
        mesh = meshio.Mesh(
            points=np.zeros((3, 3)),
            cells=[("triangle", np.array([[0, 1, 2]]))],
        )
        with pytest.raises(ValueError, match="point field"):
            write_fault_velocity_vtu(
                mesh, np.array([[0, 1, 2]]),
                {"Vp_m_per_s": np.zeros(2)},   # wrong: N_pt = 3
                {}, tmp_path / "o.vtu",
            )
        with pytest.raises(ValueError, match="cell field"):
            write_fault_velocity_vtu(
                mesh, np.array([[0, 1, 2]]),
                {},
                {"Vp_m_per_s": np.zeros(2)},   # wrong: N_tri = 1
                tmp_path / "o.vtu",
            )


class TestEndToEnd:
    def test_msh_to_fault_and_bulk_vtu(self, sidecar, tmp_path):
        msh = _make_msh(tmp_path / "tiny.msh")
        out_dir = tmp_path / "out"
        rc = main(
            [
                str(msh),
                "--sidecar", str(sidecar),
                "--write-bulk",
                "--out-dir", str(out_dir),
            ]
        )
        assert rc == 0

        fault_vtu = out_dir / "tiny_fault_velocity.vtu"
        bulk_vtu = out_dir / "tiny_bulk_velocity.vtu"
        summary = out_dir / "tiny_summary.json"
        assert fault_vtu.is_file()
        assert bulk_vtu.is_file()
        assert summary.is_file()

        # Fault VTU: point data exact against the analytic linear
        # field at the 4 fault vertices; orphan vertex (the tet apex,
        # index 4) is NaN.
        m = meshio.read(str(fault_vtu))
        assert m.cells[0].type == "triangle"
        assert m.cells[0].data.shape == (2, 3)
        pts = m.points
        for out_name in list(SIDECAR_FIELDS.values()) + list(
            DERIVED_FIELDS
        ):
            assert out_name in m.point_data, out_name
            assert f"{out_name}_cell" in m.cell_data, out_name
        vp = m.point_data["Vp_m_per_s"]
        want = (5000.0 + 0.1 * pts[:, 0] + 0.05 * pts[:, 1]
                - 0.2 * pts[:, 2])
        np.testing.assert_allclose(vp[:4], want[:4], rtol=1e-9)
        assert np.isnan(vp[4])
        # gmsh:physical carried as the fault tag on every cell.
        np.testing.assert_array_equal(
            m.cell_data["gmsh:physical"][0], [101, 101]
        )

        # Cell data exact at the centroids.
        cents = pts[m.cells[0].data].mean(axis=1)
        want_c = (5000.0 + 0.1 * cents[:, 0] + 0.05 * cents[:, 1]
                  - 0.2 * cents[:, 2])
        np.testing.assert_allclose(
            m.cell_data["Vp_m_per_s_cell"][0], want_c, rtol=1e-9
        )

        # Bulk VTU: one tet, centroid-sampled cell data only.
        mb = meshio.read(str(bulk_vtu))
        assert mb.cells[0].type == "tetra"
        tc = mb.points[mb.cells[0].data].mean(axis=1)
        want_b = (5000.0 + 0.1 * tc[:, 0] + 0.05 * tc[:, 1]
                  - 0.2 * tc[:, 2])
        np.testing.assert_allclose(
            mb.cell_data["Vp_m_per_s_cell"][0], want_b, rtol=1e-9
        )
        assert not mb.point_data  # bulk writer is cell-data only

        # Summary JSON schema.
        with open(summary) as f:
            payload = json.load(f)
        assert payload["fault_n_cells"] == 2
        assert payload["bulk_n_cells"] == 1
        assert payload["fault_n_points"] == 4
        assert payload["n_clamped_to_grid_bbox"]["fault_points"] == 0
        stats = payload["stats"]
        assert "fault_cell" in stats and "fault_point" in stats
        assert "bulk_cell" in stats
        assert stats["fault_point"]["Vp_m_per_s"]["nan_count"] == 0

    def test_missing_fault_group_errors(self, sidecar, tmp_path):
        msh = _make_msh(tmp_path / "tiny.msh")
        with pytest.raises(ValueError, match="physical group"):
            main(
                [
                    str(msh),
                    "--sidecar", str(sidecar),
                    "--fault-name", "no_such_group",
                    "--out-dir", str(tmp_path / "out"),
                ]
            )
