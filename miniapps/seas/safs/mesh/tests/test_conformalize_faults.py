"""Phase 2 acceptance tests for conformalize_faults.py.

Phase 4 of PLAN_cgal_corefine.md: this file exercises the LEGACY
Python conformalizer.  The production path is now
`tools/corefine_faults` (CGAL).  This test module is annotated with
`pytest.mark.legacy` (see `mesh/tests/pytest.ini`) and runs only
when invoked explicitly via `pytest -m legacy`.  See
`mesh/tests/README.md` for the 90-day retention clock policy
(P-013).

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_conformalize_faults.py -v -m legacy
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.legacy

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import fault_intersect as fi                        # noqa: E402
import conformalize_faults as cf                    # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _square_in_y0(half: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    V = np.array([
        (-half, 0.0, -half),
        (+half, 0.0, -half),
        (+half, 0.0, +half),
        (-half, 0.0, +half),
    ], dtype=np.float64)
    T = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    return V, T


def _square_in_x0(half: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    V = np.array([
        (0.0, -half, -half),
        (0.0, +half, -half),
        (0.0, +half, +half),
        (0.0, -half, +half),
    ], dtype=np.float64)
    T = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    return V, T


def _write_simple_stl(path: Path, name: str,
                       V: np.ndarray, T: np.ndarray) -> None:
    with open(path, "w") as fh:
        fh.write(f"solid {name}\n")
        for tri in T:
            a, b, c = V[tri[0]], V[tri[1]], V[tri[2]]
            n = np.cross(b - a, c - a)
            nrm = np.linalg.norm(n)
            n = n / nrm if nrm > 0 else np.array([0.0, 0.0, 1.0])
            fh.write(f"  facet normal {n[0]:.17e} {n[1]:.17e} {n[2]:.17e}\n")
            fh.write("    outer loop\n")
            for v in (a, b, c):
                fh.write(
                    f"      vertex {v[0]:.17e} {v[1]:.17e} {v[2]:.17e}\n"
                )
            fh.write("    endloop\n  endfacet\n")
        fh.write(f"endsolid {name}\n")


# ---------------------------------------------------------------------------
# Synthetic 90°-cross — primary acceptance test
# ---------------------------------------------------------------------------
def test_two_perpendicular_squares_conformalize(tmp_path):
    """Synthetic acceptance: two 90°-crossing 1m squares. After
    conformalize, polyline edge appears as exactly one new edge in
    each fault; the interior-crossing gate reports 0; manifold gate
    PASS; polyline edge-coincidence PASS."""
    V_A, T_A = _square_in_y0(1.0)
    V_B, T_B = _square_in_x0(1.0)
    in_dir = tmp_path / "in"; in_dir.mkdir()
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write_simple_stl(in_dir / "sqA.stl", "sqA", V_A, T_A)
    _write_simple_stl(in_dir / "sqB.stl", "sqB", V_B, T_B)
    rep = cf.conformalize(in_dir, out_dir, ["sqA", "sqB"])

    # Output STLs exist.
    assert (out_dir / "sqA.stl").exists()
    assert (out_dir / "sqB.stl").exists()
    assert (out_dir / "triangle_to_fault.json").exists()
    assert (out_dir / "intersection_report.json").exists()

    # All 4 gates passed (the conformalize call would have raised otherwise).
    pair_key = "sqA__x__sqB"
    gates = rep["pairs"][pair_key]["gates"]
    assert gates["manifold_A"] == "PASS" or gates["manifold_A"].startswith("PASS_WITH")
    assert gates["manifold_B"] == "PASS" or gates["manifold_B"].startswith("PASS_WITH")
    assert gates["polyline_edge_coincidence"] == "PASS"
    assert gates["interior_crossing_only"] == "PASS"

    # Per-fault counts: each square split from 2 → 6 triangles
    # (canonical 3 children per parent, 2 parents per square).
    assert rep["pairs"][pair_key]["post_split_n_tri_A"] == 6
    assert rep["pairs"][pair_key]["post_split_n_tri_B"] == 6
    # Vertices: 4 originals + 3 polyline pierces = 7 per square.
    assert rep["per_fault"]["sqA"]["n_vertices"] == 7
    assert rep["per_fault"]["sqB"]["n_vertices"] == 7

    # triangle_to_fault.json schema-valid.
    t2f = json.loads((out_dir / "triangle_to_fault.json").read_text())
    assert t2f["schema_version"] == 1
    assert t2f["faults"]["sqA"]["n_triangles"] == 6
    assert t2f["faults"]["sqB"]["n_triangles"] == 6
    assert t2f["faults"]["sqA"]["range"] == [0, 6]
    assert t2f["faults"]["sqB"]["range"] == [6, 12]
    assert t2f["n_total_triangles"] == 12


def test_refine_uniform_midpoint_works():
    """Knob C: midpoint refinement reduces max edge below target.
    10m square split into 2 tris with 28.3m diagonal → with target 5m
    needs 3 rounds (28.3 → 14.1 → 7.07 → 3.54 ≤ 5).  4^3 × 2 = 128 tris."""
    V, T = _square_in_y0(10.0)
    target = 5.0
    V_r, T_r = fi.refine_uniform_midpoint(V, T, target)
    a, b, c = V_r[T_r[:, 0]], V_r[T_r[:, 1]], V_r[T_r[:, 2]]
    max_edge = float(np.maximum.reduce([
        np.linalg.norm(b - a, axis=1),
        np.linalg.norm(c - b, axis=1),
        np.linalg.norm(a - c, axis=1),
    ]).max())
    assert max_edge <= target, f"max edge {max_edge:.2f} > target {target}"
    assert T_r.shape[0] == 128, T_r.shape[0]


def test_refine_below_target_is_no_op():
    """Knob C: input already at target → unchanged."""
    V, T = _square_in_y0(0.5)         # 0.5m square, max edge ≈ 0.71m
    V_r, T_r = fi.refine_uniform_midpoint(V, T, target_edge_length_m=10.0)
    assert V_r.shape[0] == V.shape[0]
    assert T_r.shape[0] == T.shape[0]


def test_conformalize_with_knob_c_synthetic(tmp_path):
    """End-to-end: Knob C on a tilted-cross test.  Squares tilted
    slightly off-perpendicular so refined-mesh edge midpoints don't
    fall exactly on the other plane (which would trigger the
    ambiguity skip without gmpy2).  Refinement to target=2.5m,
    then conformalize, then verify all 4 gates pass."""
    # Square A in the y=ε plane (slightly off y=0).
    V_A = np.array([
        [-10, 0.5, -10], [10, 0.5, -10], [10, 0.5, 10], [-10, 0.5, 10],
    ], dtype=np.float64)
    T_A = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    # Square B in the x=ε plane.
    V_B = np.array([
        [0.3, -10, -10], [0.3, 10, -10], [0.3, 10, 10], [0.3, -10, 10],
    ], dtype=np.float64)
    T_B = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    in_dir = tmp_path / "in"; in_dir.mkdir()
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write_simple_stl(in_dir / "sqA.stl", "sqA", V_A, T_A)
    _write_simple_stl(in_dir / "sqB.stl", "sqB", V_B, T_B)
    rep = cf.conformalize(
        in_dir, out_dir, ["sqA", "sqB"],
        target_edge_length_m=2.5,
    )
    pair_key = "sqA__x__sqB"
    assert pair_key in rep["pairs"], (
        f"expected pair {pair_key} in report; got {list(rep['pairs'])}"
    )
    gates = rep["pairs"][pair_key]["gates"]
    assert gates["manifold_A"] == "PASS" or gates["manifold_A"].startswith("PASS_WITH")
    assert gates["manifold_B"] == "PASS" or gates["manifold_B"].startswith("PASS_WITH")
    assert gates["polyline_edge_coincidence"] == "PASS"
    assert gates["interior_crossing_only"] == "PASS"
    # Each square refined from 2 → 128 tris (3 rounds of midpoint),
    # then more children added by polyline insertion.
    assert rep["pairs"][pair_key]["post_split_n_tri_A"] >= 128
    assert rep["pairs"][pair_key]["post_split_n_tri_B"] >= 128


def test_disjoint_pair_no_op(tmp_path):
    """Two faults far apart → no crossings → no per-pair report entry,
    but per-fault output exists and the input mesh is preserved
    (triangle count unchanged)."""
    V_A, T_A = _square_in_y0(1.0)
    V_B = V_A.copy() + np.array([100.0, 100.0, 100.0])
    T_B = T_A.copy()
    in_dir = tmp_path / "in"; in_dir.mkdir()
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write_simple_stl(in_dir / "fA.stl", "fA", V_A, T_A)
    _write_simple_stl(in_dir / "fB.stl", "fB", V_B, T_B)
    rep = cf.conformalize(in_dir, out_dir, ["fA", "fB"])
    assert rep["pairs"] == {}
    # Per-fault triangle count unchanged (2 each).
    assert rep["per_fault"]["fA"]["n_triangles"] == 2
    assert rep["per_fault"]["fB"]["n_triangles"] == 2


def test_single_fault_passthrough(tmp_path):
    """One fault, no pair to process → conformalize is a no-op
    pass-through.  Verifies the smoke-test invariant that the existing
    single-fault pipeline (Mill Creek smoke) is not regressed by the
    new step."""
    V_A, T_A = _square_in_y0(1.0)
    in_dir = tmp_path / "in"; in_dir.mkdir()
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write_simple_stl(in_dir / "fA.stl", "fA", V_A, T_A)
    rep = cf.conformalize(in_dir, out_dir, ["fA"])
    assert rep["pairs"] == {}
    assert rep["per_fault"]["fA"]["n_triangles"] == 2


# ---------------------------------------------------------------------------
# Validation gate primitives — direct unit tests
# ---------------------------------------------------------------------------
def test_is_manifold_no_t_junctions_clean_mesh():
    """A clean 2-triangle mesh sharing one edge has no T-junctions."""
    V = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]], dtype=np.float64)
    T = np.array([[0,1,2],[0,2,3]], dtype=np.int64)
    ok, bad = cf.is_manifold_no_t_junctions(V, T)
    assert ok
    assert bad == []


def test_is_manifold_no_t_junctions_detects_t_junction():
    """A T-junction: triangle T_a is split with a midpoint vertex on
    edge (u, v); neighbour T_b has the same edge undivided.  The
    midpoint vertex sits on T_b's edge — that's the T-junction."""
    V = np.array([
        [0, 0, 0],     # 0  u
        [2, 0, 0],     # 1  v
        [1, 1, 0],     # 2  apex of T_a
        [1, -1, 0],    # 3  apex of T_b
        [1, 0, 0],     # 4  midpoint on edge (0, 1)
    ], dtype=np.float64)
    T = np.array([
        [0, 4, 2],     # T_a left half
        [4, 1, 2],     # T_a right half
        [0, 1, 3],     # T_b — undivided, contains vertex 4 on its edge
    ], dtype=np.int64)
    ok, bad = cf.is_manifold_no_t_junctions(V, T)
    assert not ok, "expected manifold gate to detect the T-junction"
    assert any({0, 1} == set(e) for e in bad), bad


def test_is_manifold_detects_three_fold_edge():
    """Three triangles sharing one edge — non-manifold."""
    V = np.array([
        [0,0,0],[1,0,0],[0,1,0],[0,-1,0],[0,0,1],
    ], dtype=np.float64)
    T = np.array([
        [0,1,2],
        [0,1,3],
        [0,1,4],
    ], dtype=np.int64)
    ok, bad = cf.is_manifold_no_t_junctions(V, T)
    assert not ok
    assert (0, 1) in bad


def test_verify_polyline_in_mesh_positive():
    """A polyline whose vertices are exactly mesh vertices and whose
    consecutive pairs are mesh edges → gate PASS."""
    V = np.array([[0,0,0],[1,0,0],[2,0,0]], dtype=np.float64)
    T = np.array([[0,1,2]], dtype=np.int64)        # one degenerate-ish tri
    pl = fi.Polyline(
        points=[(0,0,0),(1,0,0),(2,0,0)],
        closed=False,
        fault_a="A", fault_b="B",
    )
    assert cf.verify_polyline_in_mesh(V, T, pl)


def test_verify_polyline_in_mesh_missing_edge():
    """Polyline vertex pair is NOT an edge of the mesh → gate FAIL."""
    V = np.array([[0,0,0],[1,0,0],[5,0,0]], dtype=np.float64)
    T = np.array([[0,1,2]], dtype=np.int64)
    pl = fi.Polyline(
        points=[(0,0,0),(5,0,0)],     # vertex pair (0, 5) is in T,
                                        # but the path between them
                                        # via mesh edges goes via vertex 1
        closed=False,
        fault_a="A", fault_b="B",
    )
    # Actually edge (0, 5) IS in T (triangle is [0, 1, 2] = indices to V's
    # rows; edges are (0,1), (1,2), (2,0) — i.e., V[0]-V[2] is an edge
    # in mesh terms despite being collinear with V[1].  So this test
    # should PASS for these data; flip to a definitive miss:
    V2 = np.array([[0,0,0],[1,0,0],[2,0,0],[3,0,0]], dtype=np.float64)
    T2 = np.array([[0,1,2]], dtype=np.int64)
    pl2 = fi.Polyline(
        points=[(0,0,0),(3,0,0)],     # vertex 3 is not in T2
        closed=False,
        fault_a="A", fault_b="B",
    )
    # vertex 3 in V2 is at distance 0 from polyline (3,0,0) point — match
    # within tol — so we'll find the polyline vertex.  But edge (0, 3)
    # is NOT in the triangle, so the gate should report False.
    assert not cf.verify_polyline_in_mesh(V2, T2, pl2)


def test_verify_polyline_in_mesh_tolerance():
    """Mesh-vertex round-off below match_tol_m is OK; large offset is not."""
    V_clean = np.array([[0,0,0],[1,0,0]], dtype=np.float64)
    T = np.array([[0,1,0]], dtype=np.int64)
    pl = fi.Polyline(
        points=[(1e-5, 0, 0), (1.0 - 1e-5, 0, 0)],   # ~10 µm off
        closed=False, fault_a="A", fault_b="B",
    )
    assert cf.verify_polyline_in_mesh(V_clean, T, pl, match_tol_m=2e-2)
    pl_far = fi.Polyline(
        points=[(0.5, 0, 0), (1.5, 0, 0)],         # 0.5m off → far above tol
        closed=False, fault_a="A", fault_b="B",
    )
    assert not cf.verify_polyline_in_mesh(V_clean, T, pl_far,
                                            match_tol_m=2e-2)


# ---------------------------------------------------------------------------
# fault_intersect.dedup_mesh
# ---------------------------------------------------------------------------
def test_dedup_mesh_collapses_coincident_vertices():
    """STL-style 3-slot-per-tri input → deduplicated mesh with shared
    vertex indices for adjacent triangles."""
    V = np.array([
        [0,0,0],[1,0,0],[0,1,0],     # tri 0
        [0,0,0],[0,1,0],[1,1,0],     # tri 1, shares (0,0,0) and (0,1,0)
    ], dtype=np.float64)
    T = np.array([[0,1,2],[3,4,5]], dtype=np.int64)
    V_out, T_out = fi.dedup_mesh(V, T, snap_m=1e-6)
    assert V_out.shape == (4, 3)
    # The shared edge appears in both triangles.
    edges_set: set[tuple[int, int]] = set()
    for tri in T_out:
        a, b, c = (int(x) for x in tri)
        for u, v in ((a, b), (b, c), (c, a)):
            edges_set.add((min(u, v), max(u, v)))
    # In the deduped 4-vertex mesh, vertex (0,0,0) is index ?, (0,1,0)
    # is some index — there should be a shared edge between them.
    # Easier check: edge count = 5 (3+3 per tri minus 1 shared).
    assert len(edges_set) == 5


def test_dedup_mesh_preserves_unique_vertices():
    """Already-unique input → vertex count unchanged (3 in, 3 out).
    np.unique may reorder rows, so we don't compare T_out to T directly;
    instead we verify the geometry round-trips: V_out[T_out] equals
    V[T] (modulo row order)."""
    V = np.array([[0,0,0],[1,0,0],[0,1,0]], dtype=np.float64)
    T = np.array([[0,1,2]], dtype=np.int64)
    V_out, T_out = fi.dedup_mesh(V, T, snap_m=1e-6)
    assert V_out.shape == (3, 3)
    # Verify each original triangle's geometry is preserved.
    for ti in range(T.shape[0]):
        orig_coords = sorted(tuple(V[i].tolist()) for i in T[ti])
        out_coords = sorted(tuple(V_out[i].tolist()) for i in T_out[ti])
        assert orig_coords == out_coords, (orig_coords, out_coords)


# ---------------------------------------------------------------------------
# tri_tri_intersect_3d_interior_only
# ---------------------------------------------------------------------------
def test_interior_only_returns_empty_after_conformalize(tmp_path):
    """Post-conformalize integration test for tri_tri_intersect_3d_interior_only:
    after the 90°-cross is conformalized, the post-split STLs share
    polyline-edge endpoints (vertex matches across faults).  The raw
    `tri_tri_intersect_3d` would report those shared edges as full-edge
    crossings; `tri_tri_intersect_3d_interior_only` must filter them
    out and report 0.

    A purely-synthetic standalone test (two coplanar-edge triangles)
    triggers the gmpy2-absent ambiguity skip, which makes both raw and
    interior return empty — that case can't distinguish a working
    filter from a broken one.  The integration test does."""
    V_A, T_A = _square_in_y0(1.0)
    V_B, T_B = _square_in_x0(1.0)
    in_dir = tmp_path / "in"; in_dir.mkdir()
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write_simple_stl(in_dir / "sqA.stl", "sqA", V_A, T_A)
    _write_simple_stl(in_dir / "sqB.stl", "sqB", V_B, T_B)
    cf.conformalize(in_dir, out_dir, ["sqA", "sqB"])

    V_A2_raw, T_A2_raw = fi._read_stl_arrays(out_dir / "sqA.stl")
    V_B2_raw, T_B2_raw = fi._read_stl_arrays(out_dir / "sqB.stl")
    V_A2, T_A2 = fi.dedup_mesh(V_A2_raw, T_A2_raw, snap_m=1e-3)
    V_B2, T_B2 = fi.dedup_mesh(V_B2_raw, T_B2_raw, snap_m=1e-3)

    interior = fi.tri_tri_intersect_3d_interior_only(
        V_A2, T_A2, V_B2, T_B2, "sqA", "sqB"
    )
    assert interior == [], (
        f"interior-crossing gate found {len(interior)} residual "
        f"interior crossings after conformalization"
    )


# ---------------------------------------------------------------------------
# Real-fixture acceptance — Mill × SBMT-SAF
# ---------------------------------------------------------------------------
_ALL8_STL_DIR = (
    Path(__file__).resolve().parent.parent
    / "output" / "all8_2000m" / "stl"
)


@pytest.mark.skipif(
    not (_ALL8_STL_DIR / "safs_sbmt_saf.stl").exists(),
    reason="all-8 fixture not present; run run_all8_2000m.sh first",
)
def test_real_mill_x_sbmt_saf_conformalize(tmp_path):
    """PLAN Phase 2 real-fixture acceptance.

    Mill Creek strand × SBMT San Andreas at 2000m: pre-split
    794 + 1346 = 2140 triangles; post-split count grows by roughly
    2 × N_pierce_total (each polyline vertex inserted as a Steiner
    point on each side adds ≈ 2 children).

    All four gates must PASS:
      - manifold_A, manifold_B
      - polyline_edge_coincidence
      - interior_crossing_only

    Runtime ≤ 60 s on a laptop.
    """
    out_dir = tmp_path
    t0 = time.perf_counter()
    rep = cf.conformalize(
        _ALL8_STL_DIR, out_dir,
        ["safs_sbmt_millcreek", "safs_sbmt_saf"],
        clearance_m=100.0,
    )
    elapsed = time.perf_counter() - t0
    assert elapsed <= 60.0, f"conformalize took {elapsed:.2f} s (budget 60 s)"

    pair_key = "safs_sbmt_millcreek__x__safs_sbmt_saf"
    assert pair_key in rep["pairs"]
    pair = rep["pairs"][pair_key]
    # With knob C enabled (target_edge_length_m=1000), real CFM data
    # may produce a small number of T-junctions from cross-parent
    # projection round-off on non-coplanar adjacent triangles.  These
    # are tolerated as warnings (HXT meshes through them); we accept
    # both PASS and PASS_WITH_<N>_T_JUNCTIONS as success here.
    for k in ("manifold_A", "manifold_B"):
        assert (pair["gates"][k] == "PASS"
                or pair["gates"][k].startswith("PASS_WITH")), (
            f"gate {k}: {pair['gates'][k]}"
        )
    assert pair["gates"]["polyline_edge_coincidence"] == "PASS"
    assert pair["gates"]["interior_crossing_only"] == "PASS"

    # Triangle count grew on both sides.
    assert pair["post_split_n_tri_A"] > pair["pre_split_n_tri_A"]
    assert pair["post_split_n_tri_B"] > pair["pre_split_n_tri_B"]

    # Output files written.
    assert (out_dir / "safs_sbmt_millcreek.stl").exists()
    assert (out_dir / "safs_sbmt_saf.stl").exists()
    assert (out_dir / "triangle_to_fault.json").exists()
    assert (out_dir / "intersection_report.json").exists()
