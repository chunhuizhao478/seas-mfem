"""
Unit tests for `verify_onfault_stress.py` (Phase 8 of
PLAN_onfaultstress.md).

The verifier reads VTU artefacts and the analytic prediction, computes
L∞ / L2 errors, and exits non-zero on tolerance violation.  These
tests build synthetic VTUs in-memory via meshio and exercise the four
public verifier functions plus the CLI.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pytest

_THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS))
import project_to_fault_stress as ptfs
import verify_onfault_stress as vof


# ----------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------
@pytest.fixture
def safod_params() -> dict:
    """Demo SAFOD analytic parameters (from H&Z section)."""
    return {
        "SHmax": 113.0,
        "Shmin": 49.0,
        "Sv": 45.0,
        "SHmax_az_deg": 23.0,
        "depth_model": "constant",
        "SHmax_grad": 0.0,
        "Shmin_grad": 0.0,
        "Sv_grad": 0.0,
        "fault_strike_azimuth_hint": 314.0,  # SAF NW strike
        "rake_sense": "right-lateral",
    }


def _write_bulk_vtu_with_constant_sigma(
    path: Path,
    params: dict,
    n_cells: int = 8,
    perturb: dict | None = None,
):
    """Build a tiny tetrahedral mesh and write it as a VTU with
    the expected six sigma_*_Pa cell-data fields.

    The recorded values come from the analytic predictor evaluated at
    each cell centroid (so a clean read should pass with tol_pa = 1).
    `perturb` is a {field_name: delta_pa} map that adds a delta to the
    field's first cell, simulating a corrupted output.
    """
    import meshio

    rng = np.random.default_rng(42)
    pts = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [1.0, 1.0, -1.0],
    ])
    # Connect 5 tetrahedra inside the unit cube.
    tets = np.array([
        [0, 1, 2, 3],
        [1, 4, 2, 7],
        [1, 5, 3, 7],
        [2, 6, 3, 7],
        [1, 2, 3, 7],
    ])[:n_cells // 1 if n_cells <= 5 else 5]
    # Compute cell centroids.
    centroids = np.mean(pts[tets], axis=1)
    # Analytic σ_seas at each centroid.
    sigma_pa = vof._bulk_sigma_seas_pa(centroids[:, 2], params)
    cell_data = {
        "sigma_xx_Pa": [sigma_pa[:, 0, 0].copy()],
        "sigma_yy_Pa": [sigma_pa[:, 1, 1].copy()],
        "sigma_zz_Pa": [sigma_pa[:, 2, 2].copy()],
        "sigma_xy_Pa": [sigma_pa[:, 0, 1].copy()],
        "sigma_yz_Pa": [sigma_pa[:, 1, 2].copy()],
        "sigma_xz_Pa": [sigma_pa[:, 0, 2].copy()],
    }
    if perturb:
        for name, delta in perturb.items():
            cell_data[name][0][0] += float(delta)

    meshio.write_points_cells(
        str(path), pts, [("tetra", tets)],
        cell_data=cell_data,
    )
    return centroids


def _write_fault_vtu_with_predicted_components(
    path: Path,
    params: dict,
    n_tri: int = 4,
):
    """Build a tiny triangular fault VTU near z=0 and write the
    expected cell-data + point-data fields populated from the
    analytic prediction.
    """
    import meshio

    # Vertical fault triangles (strike along x, dip down -z). These
    # are NOT horizontal so per_triangle_basis_raw produces a proper
    # (s, d, n) basis.
    pts = np.array([
        [0.0, 0.0,  0.0],     # vertex 0
        [1.0, 0.0,  0.0],     # vertex 1
        [0.0, 0.0, -1.0],     # vertex 2
        [1.0, 0.0, -1.0],     # vertex 3
        [0.5, 0.0, -2.0],     # vertex 4
    ])
    tris = np.array([
        [0, 1, 3],
        [0, 3, 2],
        [2, 3, 4],
        [3, 1, 4],
    ])[:n_tri]
    centroids = np.mean(pts[tris], axis=1)
    # Reproduce the verifier's fault geometry to compute the expected
    # cell-data values.
    e1 = pts[tris[:, 1]] - pts[tris[:, 0]]
    e2 = pts[tris[:, 2]] - pts[tris[:, 0]]
    cross = np.cross(e1, e2)
    cross_norm = np.linalg.norm(cross, axis=1, keepdims=True)
    areas = 0.5 * cross_norm[:, 0]
    n_raw = cross / np.where(cross_norm > 0.0, cross_norm,
                             np.ones_like(cross_norm))
    try:
        n_oriented = ptfs.harmonise_normal_orientation(
            n_raw, params["rake_sense"], centroids,
            fault_strike_azimuth_hint_deg=params["fault_strike_azimuth_hint"],
        )
    except Exception:
        n_oriented = n_raw
    s_cell, d_cell, _ = ptfs.per_triangle_basis_raw(n_oriented)
    sigma_pa = vof._bulk_sigma_seas_pa(centroids[:, 2], params)
    Sn = np.einsum("nij,nj->ni", sigma_pa, n_oriented)

    sigma_n_cell = np.einsum("ni,ni->n", n_oriented, Sn) / 1.0e6   # MPa
    tau_s_cell   = np.einsum("ni,ni->n", s_cell, Sn) / 1.0e6
    tau_d_cell   = np.einsum("ni,ni->n", d_cell, Sn) / 1.0e6

    # Point-data: re-derive vertex basis via basis_to_node.
    s_node, d_node, n_node = ptfs.basis_to_node(
        pts, tris, s_cell, d_cell, n_oriented, areas,
    )
    sigma_pa_pts = vof._bulk_sigma_seas_pa(pts[:, 2], params)
    Sn_pts = np.einsum("nij,nj->ni", sigma_pa_pts, n_node)
    sigma_n_pt = np.einsum("ni,ni->n", n_node, Sn_pts) / 1.0e6
    tau_s_pt   = np.einsum("ni,ni->n", s_node, Sn_pts) / 1.0e6
    tau_d_pt   = np.einsum("ni,ni->n", d_node, Sn_pts) / 1.0e6

    cell_data = {
        "sigma_n_total_MPa": [sigma_n_cell],
        "tau_strike_MPa":    [tau_s_cell],
        "tau_dip_MPa":       [tau_d_cell],
    }
    point_data = {
        "sigma_n_total_MPa": sigma_n_pt,
        "tau_strike_MPa":    tau_s_pt,
        "tau_dip_MPa":       tau_d_pt,
    }
    meshio.write_points_cells(
        str(path), pts, [("triangle", tris)],
        point_data=point_data, cell_data=cell_data,
    )
    return centroids


# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------
class TestPhase8VerifierAPI:
    """Exercise the four public verifier functions."""

    def test_bulk_cell_data_clean(self, tmp_path, safod_params):
        vtu = tmp_path / "bulk.vtu"
        _write_bulk_vtu_with_constant_sigma(vtu, safod_params)
        results = vof.verify_bulk_cell_data(vtu, safod_params)
        assert len(results) == 6, "six sigma_* components verified"
        for r in results:
            assert r.passed, f"{r.field_name} L_inf={r.l_inf_err}"
            assert r.l_inf_err < vof.DEFAULT_TOL_PA_BULK

    def test_bulk_cell_data_detects_corruption(self, tmp_path, safod_params):
        vtu = tmp_path / "bulk_bad.vtu"
        # Add a 1 GPa offset to the first cell of sigma_xx_Pa.
        _write_bulk_vtu_with_constant_sigma(
            vtu, safod_params,
            perturb={"sigma_xx_Pa": 1.0e9},
        )
        results = vof.verify_bulk_cell_data(vtu, safod_params)
        bad = [r for r in results if not r.passed]
        assert len(bad) == 1, "exactly one field is flagged"
        assert bad[0].field_name == "sigma_xx_Pa"
        assert bad[0].l_inf_err >= 1.0e9 - 1.0

    def test_fault_cell_data_clean(self, tmp_path, safod_params):
        vtu = tmp_path / "fault.vtu"
        _write_fault_vtu_with_predicted_components(vtu, safod_params)
        results = vof.verify_fault_cell_data(vtu, safod_params)
        assert len(results) == 3
        for r in results:
            assert r.passed, f"{r.field_name} L_inf={r.l_inf_err}"

    def test_fault_point_data_clean(self, tmp_path, safod_params):
        vtu = tmp_path / "fault_pt.vtu"
        _write_fault_vtu_with_predicted_components(vtu, safod_params)
        results = vof.verify_fault_point_data(vtu, safod_params)
        assert len(results) == 3
        for r in results:
            assert r.passed, f"{r.field_name} L_inf={r.l_inf_err}"

    def test_verification_result_to_dict(self):
        r = vof.VerificationResult(
            field_name="x",
            l_inf_err=1.0, l_inf_loc=np.array([1.0, 2.0, 3.0]),
            l2_err=0.5,
            observed_min=0.0, observed_max=1.0,
            passed=True, tol_abs=1.0, tol_rel=1.0e-6,
        )
        d = r.to_dict()
        assert d["field_name"] == "x"
        assert d["l_inf_loc"] == [1.0, 2.0, 3.0]


class TestPhase8CLI:
    """Exercise the verifier CLI."""

    def test_cli_clean_bulk_exits_zero(self, tmp_path, safod_params):
        vtu = tmp_path / "bulk.vtu"
        _write_bulk_vtu_with_constant_sigma(vtu, safod_params)
        params_json = tmp_path / "params.json"
        with open(params_json, "w") as fh:
            json.dump(safod_params, fh)
        rc = vof.main([
            "--bulk-vtu", str(vtu),
            "--params", str(params_json),
        ])
        assert rc == 0

    def test_cli_corrupt_bulk_exits_nonzero(self, tmp_path, safod_params):
        vtu = tmp_path / "bulk_bad.vtu"
        _write_bulk_vtu_with_constant_sigma(
            vtu, safod_params,
            perturb={"sigma_xx_Pa": 1.0e9},
        )
        params_json = tmp_path / "params.json"
        with open(params_json, "w") as fh:
            json.dump(safod_params, fh)
        rc = vof.main([
            "--bulk-vtu", str(vtu),
            "--params", str(params_json),
        ])
        assert rc != 0

    def test_cli_writes_report_json(self, tmp_path, safod_params):
        vtu = tmp_path / "bulk.vtu"
        _write_bulk_vtu_with_constant_sigma(vtu, safod_params)
        params_json = tmp_path / "params.json"
        with open(params_json, "w") as fh:
            json.dump(safod_params, fh)
        report = tmp_path / "report.json"
        rc = vof.main([
            "--bulk-vtu", str(vtu),
            "--params", str(params_json),
            "--report-json", str(report),
        ])
        assert rc == 0
        assert report.exists()
        data = json.loads(report.read_text())
        assert isinstance(data, list)
        assert all("field_name" in r for r in data)


class TestPhase8Tolerances:
    """Tolerance budget defaults (per plan §2086-2090)."""

    def test_default_tolerances_exist(self):
        assert vof.DEFAULT_TOL_PA_CELL == 1.0
        assert vof.DEFAULT_TOL_PA_NODE == 1.0e3
        assert vof.DEFAULT_TOL_PA_BULK == 1.0
        assert vof.DEFAULT_TOL_PA_H1_PROJ == 1.0e3
        assert vof.DEFAULT_REL_TOL == 1.0e-6


class TestPhase8MixedCellBlocks:
    """R-102 regression: mixed cell-block bulk VTUs must not crash.

    meshio enforces cell_data alignment with cell blocks (each entry
    must match its block size), so the realistic mixed-block scenario
    is: tetra + triangle blocks where sigma_*_Pa is defined on both.
    The verifier must concatenate centroids and cell_data across both
    blocks without shape mismatch (R-102: _cell_centroids_with_field
    aligns with _cell_data_array's concatenation).
    """

    def test_R102_mixed_cell_blocks_concatenated(self, tmp_path,
                                                  safod_params):
        import meshio

        # Tetra block (4 verts) + triangle block (3 verts).
        pts = np.array([
            [0.0, 0.0,  0.0],
            [1.0, 0.0,  0.0],
            [0.0, 1.0,  0.0],
            [0.0, 0.0, -1.0],
            [2.0, 0.0,  0.0],
            [3.0, 0.0,  0.0],
            [2.0, 0.0, -1.0],
        ])
        cells = [
            ("tetra",    np.array([[0, 1, 2, 3]])),
            ("triangle", np.array([[4, 5, 6]])),
        ]
        # Centroids per block.
        c_tet = np.mean(pts[cells[0][1][0]], axis=0)
        c_tri = np.mean(pts[cells[1][1][0]], axis=0)
        # Predicted σ_xx at each centroid.
        sigma_tet = vof._bulk_sigma_seas_pa(np.array([c_tet[2]]),
                                             safod_params)[0, 0, 0]
        sigma_tri = vof._bulk_sigma_seas_pa(np.array([c_tri[2]]),
                                             safod_params)[0, 0, 0]
        cell_data = {
            "sigma_xx_Pa": [
                np.array([sigma_tet]),
                np.array([sigma_tri]),
            ],
        }
        out = tmp_path / "mixed.vtu"
        meshio.write_points_cells(
            str(out), pts, cells, cell_data=cell_data,
        )
        # Must not raise a shape mismatch; both cells must verify.
        results = vof.verify_bulk_cell_data(out, safod_params)
        assert len(results) == 1
        assert results[0].field_name == "sigma_xx_Pa"
        assert results[0].passed
