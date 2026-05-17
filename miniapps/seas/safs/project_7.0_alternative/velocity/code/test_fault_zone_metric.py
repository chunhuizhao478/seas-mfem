"""Tests for fault_zone_metric.py — Phase 1 of fault_zone_projection_plan
(v3) PLAN.md.

Test catalog (T-1-1 .. T-1-4) per PLAN.md §Phase 1 → Acceptance Criteria.

Run with::

    cd safs/project_7.0_alternative/velocity/code && \\
        pytest -q data_projection/test_fault_zone_metric.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fault_zone_metric import (                                     # noqa: E402
    BandRow,
    DEFAULT_BANDS,
    _p1_proj_at_cubature,
    _p2_basis,
    _p2_proj_at_cubature,
    _split_conn_by_offset,
    per_tet_l2_by_band,
)
from sidecar import write_sidecar                                   # noqa: E402

from plot_comparison_with_dg0 import TET_CUB_BC, TET_CUB_W          # noqa: E402


# ---------------------------------------------------------------------
# Synthetic helpers
# ---------------------------------------------------------------------

def _write_unit_cube_stl(path: Path) -> None:
    """One degenerate STL triangle far below the cube so distance-to-
    fault is always large; bands collapse to the [3000, inf) bucket
    when we don't care about distance stratification."""
    # ASCII STL: a single tiny triangle at y=1e8, far enough that all
    # cubature points fall in the last band.
    path.write_text(
        "solid synthetic\n"
        "facet normal 0 1 0\n"
        "  outer loop\n"
        "    vertex 0 100000000 0\n"
        "    vertex 1 100000000 0\n"
        "    vertex 0 100000000 1\n"
        "  endloop\n"
        "endfacet\n"
        "endsolid synthetic\n"
    )


def _write_unit_cube_p1_vtu(path: Path, vert_vals: np.ndarray,
                            field_name: str = "Vs") -> np.ndarray:
    """Write an MFEM-style P1 vtu of a unit cube split into 5 tets.
    Returns the (8, 3) vertex coordinates."""
    pts = np.array([
        [0.0, 0.0, 0.0],   # 0
        [1.0, 0.0, 0.0],   # 1
        [1.0, 1.0, 0.0],   # 2
        [0.0, 1.0, 0.0],   # 3
        [0.0, 0.0, 1.0],   # 4
        [1.0, 0.0, 1.0],   # 5
        [1.0, 1.0, 1.0],   # 6
        [0.0, 1.0, 1.0],   # 7
    ], dtype=np.float64)
    # 5-tet decomposition of a unit cube.
    tets = np.array([
        [0, 1, 2, 5],
        [0, 2, 3, 7],
        [0, 4, 5, 7],
        [2, 5, 6, 7],
        [0, 2, 5, 7],
    ], dtype=np.int64)
    n_tets = tets.shape[0]
    conn_flat = tets.ravel()
    offsets = np.cumsum([4] * n_tets)
    types = np.full(n_tets, 10, dtype=np.int64)        # VTK_TETRA = 10

    def _arr(name, vals, dtype="Float64", n_components=1):
        flat = " ".join(f"{float(v):.17g}" for v in np.asarray(vals).ravel())
        nc = f' NumberOfComponents="{n_components}"' if n_components > 1 else ""
        return (f'<DataArray type="{dtype}" Name="{name}"'
                f' format="ascii"{nc}>\n{flat}\n</DataArray>')

    pts_arr = _arr("Points", pts.ravel(), n_components=3)
    conn_arr = _arr("connectivity", conn_flat, dtype="Int64")
    offs_arr = _arr("offsets", offsets, dtype="Int64")
    types_arr = _arr("types", types, dtype="UInt8")
    field_arr = _arr(field_name, vert_vals)

    xml = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian">\n'
        '  <UnstructuredGrid>\n'
        f'    <Piece NumberOfPoints="{pts.shape[0]}" '
        f'NumberOfCells="{n_tets}">\n'
        f'      <Points>{pts_arr}</Points>\n'
        f'      <Cells>{conn_arr}{offs_arr}{types_arr}</Cells>\n'
        f'      <PointData Scalars="{field_name}">{field_arr}</PointData>\n'
        '    </Piece>\n'
        '  </UnstructuredGrid>\n'
        '</VTKFile>\n'
    )
    path.write_text(xml)
    return pts


def _write_unit_cube_p2_vtu(path: Path, node_vals: np.ndarray,
                            field_name: str = "Vs"):
    """Write an MFEM P2 vtu of a unit cube split into 5 tets, with a
    10-DOF connectivity per tet.  ``node_vals`` is indexed by the
    P2 DOF id used in the connectivity array.

    Returns ``(pts, p2_nodes, conn_p2)`` so tests can sanity-check the
    DOF ↔ coordinate mapping.
    """
    corners = np.array([
        [0.0, 0.0, 0.0],   # 0
        [1.0, 0.0, 0.0],   # 1
        [1.0, 1.0, 0.0],   # 2
        [0.0, 1.0, 0.0],   # 3
        [0.0, 0.0, 1.0],   # 4
        [1.0, 0.0, 1.0],   # 5
        [1.0, 1.0, 1.0],   # 6
        [0.0, 1.0, 1.0],   # 7
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 5],
        [0, 2, 3, 7],
        [0, 4, 5, 7],
        [2, 5, 6, 7],
        [0, 2, 5, 7],
    ], dtype=np.int64)
    n_tets = tets.shape[0]

    # Build a unique edge-mid catalogue per (corner, corner) pair across
    # all tets (corner pair as a sorted tuple, value = (Nv-1) + edge_id).
    edge_index: dict[tuple[int, int], int] = {}
    extra_pts: list[np.ndarray] = []
    conn_p2 = np.zeros((n_tets, 10), dtype=np.int64)
    edge_def = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    for ti, tet in enumerate(tets):
        # Corners first.
        conn_p2[ti, :4] = tet
        for k, (i, j) in enumerate(edge_def):
            ci, cj = int(tet[i]), int(tet[j])
            key = (min(ci, cj), max(ci, cj))
            if key not in edge_index:
                mid = 0.5 * (corners[ci] + corners[cj])
                edge_index[key] = corners.shape[0] + len(extra_pts)
                extra_pts.append(mid)
            conn_p2[ti, 4 + k] = edge_index[key]

    pts = np.vstack([corners] + ([np.stack(extra_pts)] if extra_pts else []))
    if node_vals.shape[0] != pts.shape[0]:
        raise ValueError(
            f"node_vals must have one entry per P2 DOF; got "
            f"{node_vals.shape[0]} for pts={pts.shape[0]}")

    conn_flat = conn_p2.ravel()
    offsets = np.cumsum([10] * n_tets)
    # VTK_QUADRATIC_TETRA = 24 (used by MFEM's high-order vtu output).
    types = np.full(n_tets, 24, dtype=np.int64)

    def _arr(name, vals, dtype="Float64", n_components=1):
        flat = " ".join(f"{float(v):.17g}" for v in np.asarray(vals).ravel())
        nc = f' NumberOfComponents="{n_components}"' if n_components > 1 else ""
        return (f'<DataArray type="{dtype}" Name="{name}"'
                f' format="ascii"{nc}>\n{flat}\n</DataArray>')

    pts_arr = _arr("Points", pts.ravel(), n_components=3)
    conn_arr = _arr("connectivity", conn_flat, dtype="Int64")
    offs_arr = _arr("offsets", offsets, dtype="Int64")
    types_arr = _arr("types", types, dtype="UInt8")
    field_arr = _arr(field_name, node_vals)

    xml = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian">\n'
        '  <UnstructuredGrid>\n'
        f'    <Piece NumberOfPoints="{pts.shape[0]}" '
        f'NumberOfCells="{n_tets}">\n'
        f'      <Points>{pts_arr}</Points>\n'
        f'      <Cells>{conn_arr}{offs_arr}{types_arr}</Cells>\n'
        f'      <PointData Scalars="{field_name}">{field_arr}</PointData>\n'
        '    </Piece>\n'
        '  </UnstructuredGrid>\n'
        '</VTKFile>\n'
    )
    path.write_text(xml)
    return pts, conn_p2


def _write_box_sidecar(path: Path, *,
                       Vs_fn,
                       x: np.ndarray = None,
                       y: np.ndarray = None,
                       z: np.ndarray = None) -> None:
    """Write a sidecar over a coarse box that strictly contains the
    unit-cube fixture and evaluates the supplied callable on the grid.
    """
    if x is None: x = np.linspace(-2.0, 3.0, 6)
    if y is None: y = np.linspace(-2.0, 3.0, 6)
    if z is None: z = np.linspace(-2.0, 3.0, 6)
    Vs = np.empty((x.size, y.size, z.size), dtype=np.float64)
    for ix, xv in enumerate(x):
        for iy, yv in enumerate(y):
            for iz, zv in enumerate(z):
                Vs[ix, iy, iz] = float(Vs_fn(xv, yv, zv))
    write_sidecar(
        path, x, y, z,
        fields={"Vs": Vs},
        attrs={},
        field_bounds={"Vs": (Vs.min() - 1.0, Vs.max() + 1.0, "m/s")})


# ---------------------------------------------------------------------
# T-1-1 — P1 linear field projects exactly
# ---------------------------------------------------------------------

def test_T_1_1_p1_linear_exact(tmp_path: Path):
    """A linear sidecar Vs(x,y,z) = 100x + 200y + 300z + 1000 captured
    at vertices is reproduced exactly by H1-P1 in-tet barycentric
    interpolation; per-tet L² ≤ 1e-6."""
    def Vs_fn(x, y, z): return 100.0 * x + 200.0 * y + 300.0 * z + 1000.0

    sidecar = tmp_path / "linear.h5"
    _write_box_sidecar(sidecar, Vs_fn=Vs_fn)

    vtu = tmp_path / "cube_p1.vtu"
    pts = _write_unit_cube_p1_vtu(
        vtu, vert_vals=np.array([Vs_fn(*p) for p in [
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1],
        ]]))

    stl = tmp_path / "fake.stl"
    _write_unit_cube_stl(stl)

    rows = per_tet_l2_by_band(
        vtu, sidecar, stl,
        z_target=0.5, half_thick=0.6,
        bands=((0.0, 1.0e9),))    # one band, all tets in it
    assert len(rows) == 1
    assert rows[0].n_tets == 5
    assert rows[0].l2_proj < 1.0e-6, (
        f"P1 must be exact on a linear field; got L²={rows[0].l2_proj}")


# ---------------------------------------------------------------------
# T-1-2 — constant field => L² = 0
# ---------------------------------------------------------------------

def test_T_1_2_constant_zero(tmp_path: Path):
    """Constant Vs ≡ 2500 → reported L² is identically zero."""
    def Vs_fn(x, y, z): return 2500.0

    sidecar = tmp_path / "constant.h5"
    _write_box_sidecar(sidecar, Vs_fn=Vs_fn)
    vtu = tmp_path / "cube_p1.vtu"
    _write_unit_cube_p1_vtu(vtu, vert_vals=np.full(8, 2500.0))
    stl = tmp_path / "fake.stl"
    _write_unit_cube_stl(stl)

    rows = per_tet_l2_by_band(
        vtu, sidecar, stl,
        z_target=0.5, half_thick=0.6,
        bands=((0.0, 1.0e9),))
    assert rows[0].n_tets == 5
    assert rows[0].l2_proj == 0.0


# ---------------------------------------------------------------------
# T-1-3 — P2 captures linear field exactly
# ---------------------------------------------------------------------

def test_T_1_3_p2_linear_exact(tmp_path: Path):
    """H1-P2 on a P1 linear sidecar — zero per-tet L² (P2 captures
    P1 exactly).  This is the test that pins the MFEM P2 edge ordering
    in fault_zone_metric._P2_EDGES.
    """
    def Vs_fn(x, y, z): return 100.0 * x + 200.0 * y + 300.0 * z + 1000.0

    sidecar = tmp_path / "linear.h5"
    _write_box_sidecar(sidecar, Vs_fn=Vs_fn)

    vtu = tmp_path / "cube_p2.vtu"
    # Build the P2 catalogue once with placeholder node values to learn
    # the unique-node count + coordinate map; then re-write with the
    # actual Vs values at each node coordinate.
    n_nodes = _count_unique_p2_nodes()
    pts, conn_p2 = _write_unit_cube_p2_vtu(
        vtu, node_vals=np.zeros(n_nodes))
    node_vals = np.array([Vs_fn(*p) for p in pts])
    _write_unit_cube_p2_vtu(vtu, node_vals=node_vals)

    stl = tmp_path / "fake.stl"
    _write_unit_cube_stl(stl)

    rows = per_tet_l2_by_band(
        vtu, sidecar, stl,
        z_target=0.5, half_thick=0.6,
        bands=((0.0, 1.0e9),))
    assert rows[0].n_tets == 5
    assert rows[0].l2_proj < 1.0e-6, (
        f"P2 must be exact on a P1 field; got L²={rows[0].l2_proj}")


def _count_unique_p2_nodes() -> int:
    """Build the same P2 catalogue as ``_write_unit_cube_p2_vtu`` and
    count the resulting unique nodes (corners + unique edge mids)."""
    tets = np.array([
        [0, 1, 2, 5],
        [0, 2, 3, 7],
        [0, 4, 5, 7],
        [2, 5, 6, 7],
        [0, 2, 5, 7],
    ], dtype=np.int64)
    edge_def = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    edge_index: dict[tuple[int, int], int] = {}
    n_extra = 0
    for tet in tets:
        for (i, j) in edge_def:
            ci, cj = int(tet[i]), int(tet[j])
            key = (min(ci, cj), max(ci, cj))
            if key not in edge_index:
                edge_index[key] = 8 + n_extra
                n_extra += 1
    return 8 + n_extra


# ---------------------------------------------------------------------
# T-1-4 — baseline reproduction on real SAFS data (gated)
# ---------------------------------------------------------------------

ROOT = Path("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/"
            "safs/project_7.0_alternative")


@pytest.mark.parametrize("tag,expected_global", [
    ("500m",  148.0),
    ("1000m", 207.0),
    ("2000m", 238.0),
])
def test_T_1_4_baseline_global_matches_v3(tag, expected_global):
    """The CLI's GLOBAL slab L² reproduces the v3 §0 table within ±1
    m/s on each baseline mesh.  Skipped when the SAFS sidecar / mesh
    fixtures are absent (CI may not have them)."""
    vtu = (ROOT / "velocity" / "results" / "preview"
           / f"projected_velocity_{tag}" / "Cycle000000"
           / "proc000000.vtu")
    sidecar = ROOT / "velocity" / "results" / "velocity_safs.h5"
    stl = (ROOT / "meshing" / "results" / "stl_nwcut"
           / f"SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_{tag}"
             "_clean_clip_nwcut.stl")
    if not (vtu.is_file() and sidecar.is_file() and stl.is_file()):
        pytest.skip(f"SAFS fixtures for {tag} not present")
    rows = per_tet_l2_by_band(vtu, sidecar, stl,
                              z_target=-1000.0, half_thick=200.0)
    n_all = sum(r.n_tets for r in rows
                if r.n_tets > 0 and not math.isnan(r.l2_proj))
    sumsq = sum(r.n_tets * r.l2_proj ** 2 for r in rows
                if r.n_tets > 0 and not math.isnan(r.l2_proj))
    global_l2 = math.sqrt(sumsq / n_all)
    assert abs(global_l2 - expected_global) < 1.0, (
        f"{tag}: global L²={global_l2:.2f} m/s; expected "
        f"{expected_global} ± 1 (v3 §0 reproduction)")


# ---------------------------------------------------------------------
# Unit-level kernels
# ---------------------------------------------------------------------

def test_p1_proj_at_cubature_matches_einsum():
    """_p1_proj_at_cubature is the same einsum that
    plot_comparison_with_dg0.py uses."""
    rng = np.random.default_rng(0)
    vert_vals = rng.uniform(-1, 1, size=20)
    conn = rng.integers(0, 20, size=(7, 4)).astype(np.int64)
    out = _p1_proj_at_cubature(vert_vals, conn)
    expected = np.einsum("qi,ti->tq", TET_CUB_BC, vert_vals[conn])
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_p2_basis_partition_of_unity():
    """The 10 P2 basis functions sum to 1 at every interior point and
    match the corner-only sum at vertices (just one nonzero phi)."""
    bc_centroid = np.array([[0.25, 0.25, 0.25, 0.25]])
    phi = _p2_basis(bc_centroid)
    # At the centroid each corner phi = 0.25 * (0.5 - 1) = -0.125; each
    # edge phi = 4 * 0.25 * 0.25 = 0.25. Sum = 4 * (-0.125) + 6 * 0.25
    # = -0.5 + 1.5 = 1.0.
    assert abs(phi.sum() - 1.0) < 1.0e-12

    # At each corner exactly one corner-phi = 1, all others 0.
    for i in range(4):
        bc = np.zeros((1, 4)); bc[0, i] = 1.0
        phi = _p2_basis(bc)
        assert abs(phi[0, i] - 1.0) < 1.0e-12
        # All other 9 must be 0
        zero_mask = np.ones(10, dtype=bool); zero_mask[i] = False
        assert np.allclose(phi[0, zero_mask], 0.0, atol=1.0e-12)


def test_per_voxel_midpoint_residuals_linear_exact(tmp_path: Path):
    """Metric A: sample at sidecar voxel midpoints. For a linear field
    captured exactly by H1-P1 vertex projection, residuals at midpoints
    must also be ~0 (linear is captured at every off-grid point too)."""
    from fault_zone_metric import (
        per_voxel_midpoint_residuals)
    from sidecar import write_sidecar

    def Vs_fn(x, y, z): return 100.0 * x + 200.0 * y + 300.0 * z + 1000.0
    x = np.linspace(-2.0, 3.0, 6)
    y = np.linspace(-2.0, 3.0, 6)
    z = np.linspace(-2.0, 3.0, 6)
    Vs = np.empty((6, 6, 6))
    for ix, xv in enumerate(x):
        for iy, yv in enumerate(y):
            for iz, zv in enumerate(z):
                Vs[ix, iy, iz] = Vs_fn(xv, yv, zv)
    sc = tmp_path / "linear.h5"
    write_sidecar(sc, x, y, z,
                  fields={"Vs": Vs}, attrs={},
                  field_bounds={"Vs": (Vs.min()-1, Vs.max()+1, "m/s")})

    vtu = tmp_path / "cube.vtu"
    pts = _write_unit_cube_p1_vtu(
        vtu, vert_vals=np.array([
            Vs_fn(*p) for p in [
                [0,0,0], [1,0,0], [1,1,0], [0,1,0],
                [0,0,1], [1,0,1], [1,1,1], [0,1,1],
            ]]))

    stl = tmp_path / "fake.stl"
    _write_unit_cube_stl(stl)

    samples = per_voxel_midpoint_residuals(
        vtu, sc, stl, z_target=0.5, half_thick=0.6,
        bands=((0.0, 1.0e9),))
    if samples.res.size == 0:
        pytest.skip("no voxel midpoints in slab")
    # P1 captures linear exactly even at off-grid midpoints — and
    # trilinear truth of a linear field is also exact.  Both agree.
    rms = float(np.sqrt(np.mean(samples.res ** 2)))
    assert rms < 1.0e-9, (
        f"metric-A linear field must be ~exact; got RMS={rms}")


def test_grid_point_truth_is_bit_identical(tmp_path: Path):
    """per_grid_point_residuals must use the raw sidecar voxel values
    F[i,j,k] as truth, not a trilinear interpolation. Verified by
    constructing a sidecar with a non-smooth field that disagrees
    with itself between voxels (a single-voxel spike), then sampling
    the projection at the spike's voxel coordinate.

    The 'truth' returned must equal F[i,j,k] to bit (==), and the
    'proj' must be the linear-in-tet interpolation of vertex values
    (NOT the spike-aware truth)."""
    from fault_zone_metric import per_grid_point_residuals

    # Sidecar: zero everywhere except one voxel.
    x = np.linspace(-2.0, 3.0, 6)
    y = np.linspace(-2.0, 3.0, 6)
    z = np.linspace(-2.0, 3.0, 6)
    Vs = np.zeros((6, 6, 6), dtype=np.float64)
    # Put the spike at a voxel that lies inside the unit-cube fixture
    # at (0.5, 0.5, 0.5).  That maps to (i, j, k) = (3, 3, 3) since
    # x[3] = 1.0 (closest > 0.5)? Let's check: x = [-2, -1, 0, 1, 2, 3]
    # So x[2] = 0, x[3] = 1, gz[2] = 0, gz[3] = 1.  No exact 0.5
    # match, but the voxel at (3, 3, 3) = (1, 1, 1) lies on a corner
    # of the unit cube so it's still in/on the mesh.  Make spike at
    # (i, j, k) = (3, 3, 3).
    Vs[3, 3, 3] = 7777.7
    from sidecar import write_sidecar
    sc = tmp_path / "spike.h5"
    write_sidecar(sc, x, y, z,
                  fields={"Vs": Vs}, attrs={},
                  field_bounds={"Vs": (-1.0, 9999.0, "m/s")})

    # Mesh: P1 unit cube; vertex Vs values from trilinear-eval-at-vertex
    # of the spike sidecar.  All vertices land at (0/1, 0/1, 0/1)
    # which are sidecar grid corners; Vs at vertex (1,1,1) = spike.
    vtu = tmp_path / "cube.vtu"
    pts = _write_unit_cube_p1_vtu(
        vtu, vert_vals=np.array([
            Vs[2, 2, 2], Vs[3, 2, 2], Vs[3, 3, 2], Vs[2, 3, 2],
            Vs[2, 2, 3], Vs[3, 2, 3], Vs[3, 3, 3], Vs[2, 3, 3],
        ]))

    stl = tmp_path / "fake.stl"
    _write_unit_cube_stl(stl)

    # Sample at every sidecar voxel that lies on the slab z=0.5.
    # gz = [-2,-1,0,1,2,3]; pick z_target=0.5, half_thick=0.6 -> only
    # gz[3] = 1 hits.  But gz[3]=1 is exactly at a vertex so the locator
    # may snap either side.  Use z_target=0.999, half_thick=0.05 to
    # uniquely select gz[3].
    samples = per_grid_point_residuals(
        vtu, sc, stl,
        z_target=0.999, half_thick=0.05,
        bands=((0.0, 1.0e9),))

    # The spike grid point is at (gx[3], gy[3], gz[3]) = (1, 1, 1) which
    # is corner v6 of the cube.  truth at that point must equal Vs[3,3,3]
    # = 7777.7 *bit-exactly* (since we read it straight from F).
    spike_idx = np.flatnonzero(
        np.isclose(samples.xyz[:, 0], 1.0) &
        np.isclose(samples.xyz[:, 1], 1.0) &
        np.isclose(samples.xyz[:, 2], 1.0))
    assert spike_idx.size == 1, (
        f"expected one sample at (1,1,1); got {spike_idx.size}")
    assert samples.truth[spike_idx[0]] == 7777.7, (
        "truth at the spike voxel must be the raw sidecar value, "
        "not a trilinearly-smoothed average")


def test_split_conn_by_offset():
    """P1 (4 DOFs/elem) and P2 (10 DOFs/elem) tet meshes are detected
    correctly; mixed-cell-type meshes raise."""
    # P1: 3 tets
    conn = np.arange(12, dtype=np.int64)
    offs = np.array([4, 8, 12], dtype=np.int64)
    dofs, conn_2d = _split_conn_by_offset(conn, offs)
    assert dofs == 4
    assert conn_2d.shape == (3, 4)

    # P2: 2 tets
    conn = np.arange(20, dtype=np.int64)
    offs = np.array([10, 20], dtype=np.int64)
    dofs, conn_2d = _split_conn_by_offset(conn, offs)
    assert dofs == 10
    assert conn_2d.shape == (2, 10)

    # Mixed: error
    conn = np.arange(14, dtype=np.int64)
    offs = np.array([4, 14], dtype=np.int64)
    with pytest.raises(ValueError, match="non-homogeneous"):
        _split_conn_by_offset(conn, offs)
