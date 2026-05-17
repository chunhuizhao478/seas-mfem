"""Tests for build_size_field.py — Phase 2 of fault_zone_projection_plan
(v3) PLAN.md.

R-2-G-T7: build_size_field unit (PLAN.md §Phase 2 → Acceptance Criteria).

Run with::

    cd safs/project_7.0_alternative/velocity/code && \\
        pytest -q data_projection/test_build_size_field.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from build_size_field import (                                      # noqa: E402
    build_pos,
    compute_size_field,
    write_pos_scalar_point,
)
from sidecar import write_sidecar                                   # noqa: E402


# ---------------------------------------------------------------------
# Pure-numerics tests
# ---------------------------------------------------------------------

def test_constant_field_yields_lc_far():
    """Constant Vs ≡ 2500 → ‖∇Vs‖ ≡ 0 → LC ≡ lc_far everywhere."""
    gx = np.linspace(0, 10, 5)
    gy = np.linspace(0, 10, 5)
    gz = np.linspace(0, 10, 5)
    F = np.full((5, 5, 5), 2500.0)
    lc = compute_size_field(F, gx, gy, gz,
                            lc_near=1500.0, lc_far=10000.0,
                            lc_min=500.0, alpha=8.0,
                            smooth_sigma=0.0)
    assert lc.shape == F.shape
    assert np.all(lc == 10000.0), (
        f"constant field must yield uniform lc_far; got "
        f"min={lc.min()}, max={lc.max()}")


def test_linear_field_yields_uniform_lc():
    """Linear Vs(x,y,z) = 100x → ‖∇Vs‖ uniform → after norm=1
    everywhere → LC uniform = lc_far / (1 + alpha)."""
    gx = np.linspace(0, 10, 6)
    gy = np.linspace(0, 10, 6)
    gz = np.linspace(0, 10, 6)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    F = 100.0 * X
    alpha = 8.0
    lc = compute_size_field(F, gx, gy, gz,
                            lc_near=1500.0, lc_far=10000.0,
                            lc_min=500.0, alpha=alpha,
                            smooth_sigma=0.0)
    expected = 10000.0 / (1.0 + alpha)
    # np.gradient gives slightly less-clean numerics at the edges
    # (edge_order=2); the bulk should still match.
    bulk = lc[1:-1, 1:-1, 1:-1]
    np.testing.assert_allclose(bulk, expected, rtol=1.0e-6)


def test_lc_min_floor_clamps():
    """Sharp delta on a single voxel → ∇F spike → after normalisation
    LC at that voxel saturates at lc_min."""
    gx = np.linspace(0, 10, 7)
    gy = np.linspace(0, 10, 7)
    gz = np.linspace(0, 10, 7)
    F = np.zeros((7, 7, 7))
    F[3, 3, 3] = 1000.0     # delta
    lc = compute_size_field(F, gx, gy, gz,
                            lc_near=1500.0, lc_far=10000.0,
                            lc_min=500.0, alpha=8.0,
                            smooth_sigma=0.0)
    # The delta's gradient is the global max; after norm it sets
    # neighbour voxels to lc_far / 9 ≈ 1111 m, still above lc_min.
    # Bumping alpha up enough makes lc_min bind: alpha >= 19.
    lc_strict = compute_size_field(F, gx, gy, gz,
                                   lc_near=1500.0, lc_far=10000.0,
                                   lc_min=500.0, alpha=19.0,
                                   smooth_sigma=0.0)
    assert lc_strict.min() == pytest.approx(500.0, abs=1.0)
    # The high-alpha case must NOT push values above lc_far.
    assert lc_strict.max() <= 10000.0 + 1.0e-9


def test_zero_alpha_disables_refinement():
    """alpha = 0 → LC ≡ lc_far regardless of gradient."""
    gx = np.linspace(0, 10, 6)
    gy = np.linspace(0, 10, 6)
    gz = np.linspace(0, 10, 6)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    F = 100.0 * X + 200.0 * Y - 300.0 * Z
    lc = compute_size_field(F, gx, gy, gz,
                            lc_near=1500.0, lc_far=10000.0,
                            lc_min=500.0, alpha=0.0,
                            smooth_sigma=0.0)
    np.testing.assert_allclose(lc, 10000.0, rtol=0, atol=0)


# ---------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------

def test_negative_alpha_raises():
    F = np.zeros((4, 4, 4))
    g = np.linspace(0, 1, 4)
    with pytest.raises(ValueError, match="alpha"):
        compute_size_field(F, g, g, g,
                           lc_near=1500.0, lc_far=10000.0,
                           lc_min=500.0, alpha=-1.0)


def test_lc_min_above_lc_far_raises():
    F = np.zeros((4, 4, 4))
    g = np.linspace(0, 1, 4)
    with pytest.raises(ValueError, match="lc_min .* > lc_far"):
        compute_size_field(F, g, g, g,
                           lc_near=1500.0, lc_far=10000.0,
                           lc_min=20000.0, alpha=8.0)


def test_non_monotone_axis_raises():
    F = np.zeros((4, 4, 4))
    g = np.array([0.0, 1.0, 1.0, 2.0])    # not strictly increasing
    with pytest.raises(ValueError, match="strictly monotone"):
        compute_size_field(F, g, g, g,
                           lc_near=1500.0, lc_far=10000.0,
                           lc_min=500.0, alpha=8.0)


def test_F_shape_mismatch_raises():
    F = np.zeros((4, 4, 5))
    g = np.linspace(0, 1, 4)
    with pytest.raises(ValueError, match="shape"):
        compute_size_field(F, g, g, g,
                           lc_near=1500.0, lc_far=10000.0,
                           lc_min=500.0, alpha=8.0)


# ---------------------------------------------------------------------
# .pos round-trip via meshio
# ---------------------------------------------------------------------

def test_pos_roundtrip_via_text(tmp_path: Path):
    """write_pos_scalar_point produces a parseable gmsh PostView; the
    SP records carry the (x, y, z, value) we asked for."""
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([0.0, 1.0])
    z = np.array([-1.0, 0.0])
    vals = np.arange(3 * 2 * 2, dtype=np.float64).reshape(3, 2, 2)
    out = tmp_path / "view.pos"
    write_pos_scalar_point(out, x, y, z, vals, view_name="test_view")

    text = out.read_text()
    assert text.startswith('View "test_view" {')
    assert text.rstrip().endswith("};")
    # 3*2*2 = 12 SP records.
    assert text.count("SP(") == 12
    # First record: x=0, y=0, z=-1, value=vals[0,0,0]=0.
    assert "SP(0.000000,0.000000,-1.000000){0.000000};" in text
    # Last record: x=2, y=1, z=0, value=vals[2,1,1]=11.
    assert "SP(2.000000,1.000000,0.000000){11.000000};" in text


def test_voxel_stride_drops_to_minimum(tmp_path: Path):
    """A stride that drops the axis below 4 voxels must error out,
    per build_size_field docstring."""
    # Build a tiny sidecar.
    x = np.linspace(0.0, 1.0, 4)        # only 4 entries
    y = np.linspace(0.0, 1.0, 4)
    z = np.linspace(-1.0, 0.0, 4)
    Vs = np.full((4, 4, 4), 2500.0)
    p = tmp_path / "tiny.h5"
    write_sidecar(p, x, y, z, fields={"Vs": Vs}, attrs={},
                  field_bounds={"Vs": (0.0, 5000.0, "m/s")})
    out_pos = tmp_path / "tiny.pos"
    with pytest.raises(ValueError, match="requires >= 4"):
        build_pos(sidecar=p, out_pos=out_pos, voxel_stride=3,
                  verbose=False)


def test_build_pos_full_loop(tmp_path: Path):
    """End-to-end: write a synthetic sidecar, run build_pos, confirm
    the .pos has the expected number of SP records and records-per-axis
    matches the downsampled grid."""
    x = np.linspace(0.0, 100.0, 12)
    y = np.linspace(0.0, 100.0, 12)
    z = np.linspace(-50.0, 50.0, 12)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    Vs = 2500.0 + 100.0 * X / 100.0
    p = tmp_path / "sc.h5"
    write_sidecar(p, x, y, z, fields={"Vs": Vs}, attrs={},
                  field_bounds={"Vs": (0.0, 5000.0, "m/s")})
    out_pos = tmp_path / "sc.pos"
    meta = build_pos(sidecar=p, out_pos=out_pos, voxel_stride=2,
                     alpha=8.0, smooth_sigma=0.0, verbose=False)
    # Stride 2 over 12 axis entries → indices [0,2,4,6,8,10] (6),
    # plus index 11 appended because the last original entry must be
    # preserved → 7 voxels per axis; 7^3 = 343 SP records.
    assert meta["axes_n"] == (7, 7, 7)
    text = out_pos.read_text()
    assert text.count("SP(") == 343
