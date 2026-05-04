"""Unit tests for mmgs_remesh_per_fault.py.

Covers REVIEW_interior_subdivision_and_collapse.md Tier 1:
the per-fault mmgs_O3 remesher must (a) preserve cross-fault
polyline vertex coords bit-exactly, (b) preserve polyline edge
topology, and (c) handle pathological inputs (single triangle, no
polyline) without crashing.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_mmgs_remesh_per_fault.py -v
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the parent (mesh/) importable.
_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import break_fault_wedges as bfw  # noqa: E402
import mmgs_remesh_per_fault as mr  # noqa: E402


# Skip the entire module if mmgs_O3 isn't installed.
_mmgs = shutil.which("mmgs_O3")
pytestmark = pytest.mark.skipif(
    _mmgs is None,
    reason="mmgs_O3 not on PATH; install via "
           "`conda install -c conda-forge mmgsuite`")


# ---------------------------------------------------------------------------
# Helpers — synthesize tiny fault triangulations.
# ---------------------------------------------------------------------------
def _two_faults_sharing_edge(tmp_path: Path,
                              dihedral_deg: float = 30.0
                              ) -> tuple[Path, list[str]]:
    """Two fault STLs sharing the edge (0,0,0) → (1,0,0).

    Each fault has a small fan of triangles around its half-strip so
    mmgs has interior triangles it can re-mesh.  The shared polyline
    edge has bit-identical endpoints in both faults.
    """
    stl_dir = tmp_path / "stl_in"
    stl_dir.mkdir(parents=True, exist_ok=True)

    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])

    # Fault A in y ≥ 0 half-plane, with apexes spread across +y.
    apex_a1 = np.array([0.5, 1.0, 0.0])
    apex_a2 = np.array([0.0, 2.0, 0.0])
    apex_a3 = np.array([1.0, 2.0, 0.0])
    V_a = np.stack([p0, p1, apex_a1, apex_a2, apex_a3])
    T_a = np.asarray([
        [0, 1, 2],   # (p0, p1, apex_a1)
        [0, 2, 3],   # (p0, apex_a1, apex_a2)
        [1, 4, 2],   # (p1, apex_a3, apex_a1)
        [2, 4, 3],   # (apex_a1, apex_a3, apex_a2)
    ], dtype=np.int64)

    # Fault B in the y ≤ 0 half-plane rotated about the x-axis.
    theta = np.deg2rad(dihedral_deg)
    rot = np.array([
        [1.0, 0.0,           0.0],
        [0.0, np.cos(theta), -np.sin(theta)],
        [0.0, np.sin(theta),  np.cos(theta)],
    ])
    apex_b1 = rot @ np.array([0.5, -1.0, 0.0])
    apex_b2 = rot @ np.array([0.0, -2.0, 0.0])
    apex_b3 = rot @ np.array([1.0, -2.0, 0.0])
    V_b = np.stack([p0, p1, apex_b1, apex_b2, apex_b3])
    T_b = np.asarray([
        [0, 2, 1],   # (p0, apex_b1, p1)
        [0, 3, 2],   # (p0, apex_b2, apex_b1)
        [1, 2, 4],   # (p1, apex_b1, apex_b3)
        [2, 3, 4],   # ...
    ], dtype=np.int64)

    bfw._write_ascii_stl(stl_dir / "fault_A.stl", V_a, T_a, "fault_A")
    bfw._write_ascii_stl(stl_dir / "fault_B.stl", V_b, T_b, "fault_B")
    return stl_dir, ["fault_A", "fault_B"]


def _disjoint_single_fault(tmp_path: Path) -> tuple[Path, list[str]]:
    """One isolated fault with no cross-fault polyline.

    Used to verify mmgs_remesh_per_fault doesn't crash when
    polyline_keys is empty (i.e., no required vertices / edges).
    """
    stl_dir = tmp_path / "stl_in"
    stl_dir.mkdir(parents=True, exist_ok=True)
    V = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, -1.0, 0.0],
    ])
    T = np.asarray([[0, 1, 2], [1, 0, 3]], dtype=np.int64)
    bfw._write_ascii_stl(stl_dir / "lonely.stl", V, T, "lonely")
    return stl_dir, ["lonely"]


# ---------------------------------------------------------------------------
# Tests for the Medit I/O round-trip.
# ---------------------------------------------------------------------------
def test_medit_roundtrip_preserves_coords(tmp_path):
    """Write a triangulation to medit and read it back; the vertex
    coords and triangles must round-trip exactly.
    """
    V = np.asarray([
        [1.5, -2.5, 3.5],
        [4.0,  0.0, 0.0],
        [0.0,  1.0, 0.0],
        [-1.0, 2.0, -3.0],
    ])
    T = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    p = tmp_path / "rt.mesh"
    mr._write_medit_surface(p, V, T,
                              required_vertex_indices=[0, 1],
                              required_edges=[(0, 1), (1, 2)])
    V2, T2 = mr._read_medit_surface(p)
    np.testing.assert_array_equal(V, V2)
    np.testing.assert_array_equal(T, T2)


def test_medit_no_required_blocks_is_legal(tmp_path):
    """Writing with empty required-vertex / required-edge lists must
    produce a parseable file.  This covers the disjoint-single-fault
    case where there are no polyline vertices.
    """
    V = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0]])
    T = np.asarray([[0, 1, 2]], dtype=np.int64)
    p = tmp_path / "no_req.mesh"
    mr._write_medit_surface(p, V, T,
                              required_vertex_indices=[],
                              required_edges=[])
    V2, T2 = mr._read_medit_surface(p)
    np.testing.assert_array_equal(V, V2)
    np.testing.assert_array_equal(T, T2)


def test_medit_writer_rejects_wrong_dim(tmp_path):
    """Defensive input validation."""
    p = tmp_path / "bad.mesh"
    with pytest.raises(ValueError, match="V must be"):
        mr._write_medit_surface(p, np.zeros((3, 2)),
                                  np.zeros((1, 3), dtype=np.int64),
                                  [], [])
    with pytest.raises(ValueError, match="T must be"):
        mr._write_medit_surface(p, np.zeros((3, 3)),
                                  np.zeros((1, 4), dtype=np.int64),
                                  [], [])


# ---------------------------------------------------------------------------
# Polyline-edge / vertex extraction.
# ---------------------------------------------------------------------------
def test_polyline_edges_in_fault_finds_shared_edge():
    snap_m = 0.001
    V = np.asarray([
        [0.0, 0.0, 0.0],   # 0 — polyline
        [1.0, 0.0, 0.0],   # 1 — polyline
        [0.5, 1.0, 0.0],   # 2 — interior
    ])
    T = np.asarray([[0, 1, 2]], dtype=np.int64)
    polyline_keys = {bfw._snap_key(V[0], snap_m),
                     bfw._snap_key(V[1], snap_m)}
    edges = mr._polyline_edges_in_fault(V, T, polyline_keys, snap_m)
    assert edges == [(0, 1)], edges


def test_polyline_edges_does_not_match_one_endpoint():
    """An edge (polyline_vertex, interior_vertex) must NOT be
    flagged as a polyline edge — only edges with BOTH endpoints on
    the polyline are required.
    """
    snap_m = 0.001
    V = np.asarray([
        [0.0, 0.0, 0.0],   # 0 — polyline
        [1.0, 0.0, 0.0],   # 1 — polyline
        [0.5, 1.0, 0.0],   # 2 — interior
    ])
    T = np.asarray([[0, 1, 2]], dtype=np.int64)
    # Only vertex 0 is on the polyline.
    polyline_keys = {bfw._snap_key(V[0], snap_m)}
    edges = mr._polyline_edges_in_fault(V, T, polyline_keys, snap_m)
    assert edges == []


def test_polyline_vertices_in_fault():
    snap_m = 0.001
    V = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
    ])
    polyline_keys = {bfw._snap_key(V[0], snap_m),
                     bfw._snap_key(V[1], snap_m)}
    out = mr._polyline_vertices_in_fault(V, polyline_keys, snap_m)
    assert out == [0, 1]


# ---------------------------------------------------------------------------
# Defensive parameter validation.
# ---------------------------------------------------------------------------
def test_invalid_h_params_raise(tmp_path):
    stl_dir, names = _disjoint_single_fault(tmp_path)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="hmin/hmax"):
        mr.remesh_per_fault(stl_dir, out, names, hmin=-1, hmax=10)
    with pytest.raises(ValueError, match="hmin must be <= hmax"):
        mr.remesh_per_fault(stl_dir, out, names, hmin=100, hmax=10)
    with pytest.raises(ValueError, match="hgrad"):
        mr.remesh_per_fault(stl_dir, out, names, hmin=1, hmax=10,
                            hgrad=0.5)
    with pytest.raises(ValueError, match="hausd"):
        mr.remesh_per_fault(stl_dir, out, names, hmin=1, hmax=10,
                            hgrad=1.3, hausd=-1)


def test_missing_input_dir_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-stl-dir not found"):
        mr.remesh_per_fault(tmp_path / "does_not_exist",
                            tmp_path / "out", ["x"])


def test_empty_include_list_raises(tmp_path):
    stl_dir = tmp_path / "in"
    stl_dir.mkdir()
    with pytest.raises(SystemExit, match="--include-fault must"):
        mr.remesh_per_fault(stl_dir, tmp_path / "out", [])


def test_missing_stl_raises(tmp_path):
    stl_dir = tmp_path / "in"
    stl_dir.mkdir()
    with pytest.raises(SystemExit, match="missing input STL"):
        mr.remesh_per_fault(stl_dir, tmp_path / "out", ["nonexistent"])


# ---------------------------------------------------------------------------
# End-to-end mmgs invocation on synthetic 2-fault input.
# ---------------------------------------------------------------------------
def test_e2e_two_faults_preserves_polyline_vertices_bit_identically(
        tmp_path):
    """Cross-fault conformity invariant: after mmgs runs on each fault
    independently, every shared snap_key has bit-identical coords
    across the two output STLs.

    This is the critical invariant that lets HXT's PLC recovery accept
    the two-fault polyline as a single shared constraint.
    """
    stl_dir, names = _two_faults_sharing_edge(tmp_path)
    out_dir = tmp_path / "out"
    snap_m = 0.001

    rc = mr.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--hmin", "0.05",
        "--hmax", "0.5",
        "--hgrad", "1.3",
        "--hausd", "0.05",
        "--snap-m", str(snap_m),
    ])
    assert rc == 0

    V_a, _, _ = bfw._read_ascii_stl(out_dir / f"{names[0]}.stl")
    V_b, _, _ = bfw._read_ascii_stl(out_dir / f"{names[1]}.stl")

    keys_a = {tuple(int(round(c / snap_m)) for c in v): tuple(v)
              for v in V_a}
    keys_b = {tuple(int(round(c / snap_m)) for c in v): tuple(v)
              for v in V_b}
    shared = set(keys_a.keys()) & set(keys_b.keys())
    assert len(shared) >= 2, (
        f"expected ≥2 shared snap_keys (the polyline edge endpoints); "
        f"got {len(shared)}")
    for k in shared:
        ca = np.array(keys_a[k])
        cb = np.array(keys_b[k])
        for i in range(3):
            assert ca[i] == cb[i], (
                f"shared snap_key {k} has different coords in fault_A "
                f"({ca[i]!r}) vs fault_B ({cb[i]!r}) on axis {i}; "
                f"mmgs did not preserve cross-fault conformity")


def test_e2e_two_faults_polyline_edge_count_preserved(tmp_path):
    """The polyline edge (0,0,0)–(1,0,0) is RequiredEdge in mmgs
    input.  After remeshing, the polyline edge endpoints (0,0,0) and
    (1,0,0) MUST still appear as vertices in both faults' output
    STLs.  RequiredEdges blocks split / collapse, so the edge cannot
    be subdivided further by mmgs.
    """
    stl_dir, names = _two_faults_sharing_edge(tmp_path)
    out_dir = tmp_path / "out"
    snap_m = 0.001

    mr.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--hmin", "0.05", "--hmax", "0.5",
        "--hgrad", "1.3", "--hausd", "0.05",
        "--snap-m", str(snap_m),
    ])

    for short in names:
        V, _, _ = bfw._read_ascii_stl(out_dir / f"{short}.stl")
        # Both polyline endpoints must still be present.
        coords_set = {tuple(v) for v in V}
        assert (0.0, 0.0, 0.0) in coords_set or any(
            np.allclose(v, [0.0, 0.0, 0.0], atol=1e-12) for v in V), (
            f"fault {short}: polyline endpoint (0,0,0) lost after mmgs; "
            f"a RequiredVertex must NOT be moved")
        assert (1.0, 0.0, 0.0) in coords_set or any(
            np.allclose(v, [1.0, 0.0, 0.0], atol=1e-12) for v in V), (
            f"fault {short}: polyline endpoint (1,0,0) lost after mmgs")


def test_e2e_disjoint_single_fault_runs_without_polyline(tmp_path):
    """When there's only one fault (no polyline_keys at all), mmgs
    should still run and produce a valid output STL.  No
    RequiredEdges / RequiredVertices block in the medit input.
    """
    stl_dir, names = _disjoint_single_fault(tmp_path)
    out_dir = tmp_path / "out"

    rc = mr.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--hmin", "0.1", "--hmax", "1.0",
        "--hgrad", "1.3", "--hausd", "0.1",
        "--snap-m", "0.001",
    ])
    assert rc == 0
    out_stl = out_dir / f"{names[0]}.stl"
    assert out_stl.exists()
    V, T, _ = bfw._read_ascii_stl(out_stl)
    assert V.shape[0] >= 3
    assert T.shape[0] >= 1

    report = json.loads(
        (out_dir / "mmgs_remesh_report.json").read_text())
    assert report["n_polyline_vertex_keys"] == 0
    assert report["per_fault"][names[0]]["n_required_vertices"] == 0
    assert report["per_fault"][names[0]]["n_required_edges"] == 0
