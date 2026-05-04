"""Unit tests for dedup_coplanar_facets.py.

Covers the Step-6-unblocking fix for autorefine R-501: when two
faults are physically adjacent, CGAL 6.1 autorefine emits
overlapping triangulations whose pairs of triangles share an edge
AND lie on the same plane.  The dedup pass must remove one of each
such pair before HXT sees the input.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_dedup_coplanar_facets.py -v
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import break_fault_wedges as bfw  # noqa: E402
import dedup_coplanar_facets as dcf  # noqa: E402


# ---------------------------------------------------------------------------
# Defensive parameter validation.
# ---------------------------------------------------------------------------
def test_missing_input_dir_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-stl-dir not found"):
        dcf.dedup_coplanar(tmp_path / "missing", tmp_path / "out", ["x"])


def test_empty_include_list_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    with pytest.raises(SystemExit, match="--include-fault must"):
        dcf.dedup_coplanar(in_dir, tmp_path / "out", [])


def test_invalid_snap_m_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    with pytest.raises(ValueError, match="snap_m"):
        dcf.dedup_coplanar(in_dir, tmp_path / "out", ["x"], snap_m=-1)


def test_invalid_coplanar_tol_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    with pytest.raises(ValueError, match="coplanar_tol"):
        dcf.dedup_coplanar(in_dir, tmp_path / "out", ["x"],
                            coplanar_tol=-0.1)
    with pytest.raises(ValueError, match="coplanar_tol"):
        dcf.dedup_coplanar(in_dir, tmp_path / "out", ["x"],
                            coplanar_tol=2.0)


def test_missing_stl_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    with pytest.raises(SystemExit, match="missing input STL"):
        dcf.dedup_coplanar(in_dir, tmp_path / "out", ["nope"])


# ---------------------------------------------------------------------------
# Core algorithmic tests.
# ---------------------------------------------------------------------------
def _write_two_fault_overlap(tmp_path: Path
                              ) -> tuple[Path, list[str]]:
    """Build two faults A and B that BOTH triangulate the SAME planar
    quad (0,0,0)-(1,0,0)-(1,1,0)-(0,1,0) into two triangles, but with
    DIFFERENT diagonal cuts:

      Fault A: (0,0,0)-(1,0,0)-(1,1,0) and (0,0,0)-(1,1,0)-(0,1,0)
               — diagonal from (0,0) to (1,1)
      Fault B: (0,0,0)-(1,0,0)-(0,1,0) and (1,0,0)-(1,1,0)-(0,1,0)
               — diagonal from (1,0) to (0,1)

    The two triangulations are coplanar.  Each triangle in fault A has
    a coplanar partner in fault B sharing 2 vertices.  The dedup pass
    should drop one of each pair, leaving each fault with 1 triangle
    (or some equivalent reduction).
    """
    stl_dir = tmp_path / "stl_in"
    stl_dir.mkdir(parents=True, exist_ok=True)

    # Use shared coords so snap-key matching works.
    V_a = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    T_a = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)

    V_b = V_a.copy()
    T_b = np.asarray([[0, 1, 3], [1, 2, 3]], dtype=np.int64)

    bfw._write_ascii_stl(stl_dir / "fault_A.stl", V_a, T_a, "fault_A")
    bfw._write_ascii_stl(stl_dir / "fault_B.stl", V_b, T_b, "fault_B")
    return stl_dir, ["fault_A", "fault_B"]


def test_coplanar_overlap_pair_is_detected(tmp_path):
    """The two-fault-overlap fixture must produce ≥ 2 coplanar
    overlap pairs (each fault A triangle pairs with a fault B
    triangle along a shared edge)."""
    stl_dir, names = _write_two_fault_overlap(tmp_path)
    out_dir = tmp_path / "out"
    rc = dcf.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    assert rc == 0
    rep = json.loads(
        (out_dir / "dedup_coplanar_report.json").read_text())
    # Two faults each with 2 triangles, all 4 coplanar; non-manifold
    # edges of the unit square at the boundary edges.  Expected
    # detected pairs: at LEAST 1 (the shared diagonal edge).
    assert rep["n_coplanar_pairs_detected"] >= 1, (
        f"expected ≥1 coplanar pair; got "
        f"{rep['n_coplanar_pairs_detected']}")


def test_loser_is_from_fault_listed_later(tmp_path):
    """Selection rule: the triangle dropped must come from the fault
    that is LATER in the include-fault list.  Verify by listing
    fault_B first and seeing that fault_A's triangles get dropped.
    """
    stl_dir, names = _write_two_fault_overlap(tmp_path)
    out_dir = tmp_path / "out"
    # Reverse order: B first, A second.
    dcf.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[1],  # B first
        "--include-fault", names[0],  # A second
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    rep = json.loads(
        (out_dir / "dedup_coplanar_report.json").read_text())
    # Fault A (listed second) must have at least 1 triangle dropped;
    # fault B (listed first) must have 0 dropped.
    assert rep["per_fault"][names[0]]["n_dropped_coplanar"] >= 1, (
        f"fault_A (listed later) should have ≥ 1 coplanar drop; "
        f"got {rep['per_fault'][names[0]]['n_dropped_coplanar']}")
    assert rep["per_fault"][names[1]]["n_dropped_coplanar"] == 0, (
        f"fault_B (listed first) should be untouched; got "
        f"{rep['per_fault'][names[1]]['n_dropped_coplanar']} drops")


def test_non_coplanar_triangles_are_kept(tmp_path):
    """Two faults sharing an edge at a NON-zero dihedral must have
    NO triangles dropped.  This is a regression guard against the
    dedup pass dropping legitimate cross-fault polyline triangles.
    """
    stl_dir = tmp_path / "stl_in"
    stl_dir.mkdir(parents=True, exist_ok=True)
    # Fault A in y=0 plane.  Fault B rotated 30° about the x-axis.
    V_a = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
    ])
    T_a = np.asarray([[0, 1, 2]], dtype=np.int64)
    theta = np.deg2rad(30)
    rot = np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta),  np.cos(theta)],
    ])
    apex_b = rot @ np.array([0.5, 1.0, 0.0])
    V_b = np.stack([V_a[0], V_a[1], apex_b])
    T_b = np.asarray([[0, 1, 2]], dtype=np.int64)
    bfw._write_ascii_stl(stl_dir / "A.stl", V_a, T_a, "A")
    bfw._write_ascii_stl(stl_dir / "B.stl", V_b, T_b, "B")

    out_dir = tmp_path / "out"
    dcf.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", "A",
        "--include-fault", "B",
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    rep = json.loads(
        (out_dir / "dedup_coplanar_report.json").read_text())
    assert rep["n_coplanar_pairs_detected"] == 0, (
        f"non-coplanar pair (30° dihedral) must NOT be detected; "
        f"got {rep['n_coplanar_pairs_detected']} pairs")
    for name in ("A", "B"):
        assert rep["per_fault"][name]["n_dropped_coplanar"] == 0, (
            f"non-coplanar pair must not produce drops in {name}")


def test_disjoint_faults_have_no_pairs(tmp_path):
    """Two faults that share NO vertices must produce zero pairs."""
    stl_dir = tmp_path / "stl_in"
    stl_dir.mkdir(parents=True, exist_ok=True)
    V_a = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0],
    ])
    V_b = np.asarray([
        [10.0, 10.0, 10.0], [11.0, 10.0, 10.0], [10.5, 11.0, 10.0],
    ])
    T = np.asarray([[0, 1, 2]], dtype=np.int64)
    bfw._write_ascii_stl(stl_dir / "A.stl", V_a, T, "A")
    bfw._write_ascii_stl(stl_dir / "B.stl", V_b, T, "B")

    out_dir = tmp_path / "out"
    dcf.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", "A",
        "--include-fault", "B",
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    rep = json.loads(
        (out_dir / "dedup_coplanar_report.json").read_text())
    assert rep["n_coplanar_pairs_detected"] == 0
    for name in ("A", "B"):
        assert rep["per_fault"][name]["T_out"] == 1, (
            f"disjoint fault {name} must keep its 1 triangle")


def test_dedup_idempotent(tmp_path):
    """Running dedup twice must give the same output as running it
    once (no remaining coplanar overlaps after the first pass).
    """
    stl_dir, names = _write_two_fault_overlap(tmp_path)
    out_dir1 = tmp_path / "out1"
    out_dir2 = tmp_path / "out2"
    dcf.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir1),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    dcf.main([
        "--in-stl-dir", str(out_dir1),
        "--out-stl-dir", str(out_dir2),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--snap-m", "0.001",
        "--coplanar-tol", "1e-3",
    ])
    rep2 = json.loads(
        (out_dir2 / "dedup_coplanar_report.json").read_text())
    # Second pass should find 0 new pairs (idempotent).
    assert rep2["n_coplanar_pairs_detected"] == 0, (
        f"second pass must find 0 pairs (idempotent); got "
        f"{rep2['n_coplanar_pairs_detected']}")
    # Triangle counts must match between the two passes.
    for name in names:
        rep1 = json.loads(
            (out_dir1 / "dedup_coplanar_report.json").read_text())
        assert (rep1["per_fault"][name]["T_out"]
                == rep2["per_fault"][name]["T_out"]), (
            f"idempotence violation for {name}")
