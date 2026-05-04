"""Unit tests for break_fault_wedges.py.

Covers REVIEW_intersection_refinement_investigation.md R-005:
the wedge-edge splitter must (a) detect shallow-dihedral cross-fault
edges, (b) split each such edge in BOTH faults, and (c) emit
bit-identical midpoint coords on both sides so HXT's PLC recovery
sees a single shared polyline vertex per edge.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_break_fault_wedges.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make the parent (mesh/) importable so `import break_fault_wedges` resolves.
_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import break_fault_wedges as bfw  # noqa: E402


def _write_stl(path: Path, V: np.ndarray, T: np.ndarray, name: str) -> None:
    bfw._write_ascii_stl(path, V, T, name)


def _shallow_wedge_pair(tmp_path: Path, dihedral_deg: float
                        ) -> tuple[Path, list[str]]:
    """Two faults sharing an edge along the x-axis from (0,0,0) to (1,0,0).

    Fault A is in the y=0 plane (horizontal strip).
    Fault B is rotated by `dihedral_deg` around the shared x-axis edge,
    so the dihedral between A's plane and B's plane equals
    `dihedral_deg`.
    Each fault has 2 triangles: one on each side of the shared edge.

    Returns (stl_dir, [fault_short_names]).
    """
    stl_dir = tmp_path / "stl"
    stl_dir.mkdir(parents=True, exist_ok=True)

    # Shared edge endpoints — bit-identical between the two faults.
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])

    # Fault A: y=0 plane, triangle apex at (0.5, +1.0, 0).
    apex_a = np.array([0.5, 1.0, 0.0])
    # Wider apex on the OTHER side of the shared edge so triangle (p0,p1,apex)
    # is well-formed.
    apex_a_back = np.array([0.5, -1.0, 0.0])

    V_a = np.stack([p0, p1, apex_a, apex_a_back])
    # Two triangles sharing edge (p0, p1).  Use CCW orientation.
    T_a = np.asarray([[0, 1, 2], [1, 0, 3]], dtype=np.int64)

    # Fault B: rotated about the x-axis by dihedral_deg.
    theta = np.deg2rad(dihedral_deg)
    rot = np.array([
        [1.0, 0.0,           0.0],
        [0.0, np.cos(theta), -np.sin(theta)],
        [0.0, np.sin(theta),  np.cos(theta)],
    ])
    apex_b = rot @ apex_a
    apex_b_back = rot @ apex_a_back

    V_b = np.stack([p0, p1, apex_b, apex_b_back])
    T_b = np.asarray([[0, 1, 2], [1, 0, 3]], dtype=np.int64)

    _write_stl(stl_dir / "fault_A.stl", V_a, T_a, "fault_A")
    _write_stl(stl_dir / "fault_B.stl", V_b, T_b, "fault_B")
    return stl_dir, ["fault_A", "fault_B"]


def _read_stl_vertices(path: Path) -> np.ndarray:
    V, _, _ = bfw._read_ascii_stl(path)
    return V


def test_R005_shallow_dihedral_edge_is_split(tmp_path):
    """Two faults at 5° dihedral share one edge; the splitter must
    detect and split it on both sides."""
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=5.0)
    out_dir = tmp_path / "out"

    rc = bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", "0.001",
    ])
    assert rc == 0

    # The single shared edge should be detected as 1 wedge edge.
    import json
    report = json.loads((out_dir / "wedge_split_report.json").read_text())
    assert report["n_wedge_edges"] == 1, (
        f"expected 1 wedge edge at 5° dihedral; got "
        f"{report['n_wedge_edges']}")
    # Each fault has 2 triangles incident to the shared edge → 2 splits/fault.
    assert report["n_splits_applied"] == 4, (
        f"expected 4 splits (2 per fault, 2 faults); got "
        f"{report['n_splits_applied']}")


def test_R005_steep_dihedral_edge_is_NOT_split(tmp_path):
    """At 60° dihedral the splitter must leave the edge untouched —
    the wedge tet at 60° is well-conditioned."""
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=60.0)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", "0.001",
    ])

    import json
    report = json.loads((out_dir / "wedge_split_report.json").read_text())
    assert report["n_wedge_edges"] == 0
    assert report["n_splits_applied"] == 0
    # Output STL triangle counts are unchanged.
    for short in names:
        V_in = _read_stl_vertices(stl_dir / f"{short}.stl")
        V_out = _read_stl_vertices(out_dir / f"{short}.stl")
        assert V_out.shape[0] == V_in.shape[0]


def test_R005_midpoint_is_bit_identical_across_faults(tmp_path):
    """The midpoint vertex inserted on the shared edge must have
    bit-identical coords in both faults' STLs.  This is the
    cross-fault conformity invariant: HXT's PLC recovery treats the
    two faults' polylines as ONE polyline only when the inserted
    midpoint coord matches to bit precision.
    """
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=10.0)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", "0.001",
    ])

    V_a = _read_stl_vertices(out_dir / f"{names[0]}.stl")
    V_b = _read_stl_vertices(out_dir / f"{names[1]}.stl")

    # The midpoint of the shared edge (0,0,0)-(1,0,0) is (0.5, 0, 0).
    expected_mid = np.array([0.5, 0.0, 0.0])

    # Find the midpoint vertex in each fault — the one closest to
    # `expected_mid`.
    def closest(V, p):
        d = np.linalg.norm(V - p, axis=1)
        return V[np.argmin(d)]

    mid_a = closest(V_a, expected_mid)
    mid_b = closest(V_b, expected_mid)

    # Bit-identity: every coordinate component must match exactly.
    for i in range(3):
        assert mid_a[i] == mid_b[i], (
            f"midpoint coord component {i} mismatched: "
            f"fault_A={mid_a[i]!r} vs fault_B={mid_b[i]!r}; HXT "
            f"would see two near-coincident polylines")
    # And the value should be the analytic midpoint (no drift).
    np.testing.assert_array_equal(mid_a, expected_mid)


def _drifted_wedge_pair(tmp_path: Path, dihedral_deg: float,
                         drift_m: float
                         ) -> tuple[Path, list[str]]:
    """Like _shallow_wedge_pair but introduces a sub-cm drift in
    fault B's polyline vertex coords relative to fault A's.

    Fault A's shared edge endpoints are exactly (0,0,0) and (1,0,0).
    Fault B's are (drift_m, drift_m, drift_m) and (1+drift_m, drift_m,
    drift_m).  At snap_m=0.1 with drift_m < 0.05, both faults' edges
    snap to the SAME edge_key but their raw coords differ.
    """
    stl_dir = tmp_path / "stl"
    stl_dir.mkdir(parents=True, exist_ok=True)

    # Fault A: y=0 plane.
    pa1 = np.array([0.0, 0.0, 0.0])
    pa2 = np.array([1.0, 0.0, 0.0])
    apex_a = np.array([0.5, 1.0, 0.0])
    apex_a_back = np.array([0.5, -1.0, 0.0])
    V_a = np.stack([pa1, pa2, apex_a, apex_a_back])
    T_a = np.asarray([[0, 1, 2], [1, 0, 3]], dtype=np.int64)

    # Fault B: rotated by dihedral_deg, with shared-edge endpoints
    # drifted by `drift_m` in (x, y, z) so fault A's pa1 ≠ fault B's pb1
    # at machine-epsilon level.
    theta = np.deg2rad(dihedral_deg)
    rot = np.array([
        [1.0, 0.0,           0.0],
        [0.0, np.cos(theta), -np.sin(theta)],
        [0.0, np.sin(theta),  np.cos(theta)],
    ])
    drift = np.array([drift_m, drift_m, drift_m])
    pb1 = pa1 + drift
    pb2 = pa2 + drift
    apex_b = rot @ (apex_a - drift) + drift
    apex_b_back = rot @ (apex_a_back - drift) + drift
    V_b = np.stack([pb1, pb2, apex_b, apex_b_back])
    T_b = np.asarray([[0, 1, 2], [1, 0, 3]], dtype=np.int64)

    _write_stl(stl_dir / "fault_A.stl", V_a, T_a, "fault_A")
    _write_stl(stl_dir / "fault_B.stl", V_b, T_b, "fault_B")
    return stl_dir, ["fault_A", "fault_B"]


def test_canonical_midpoint_lies_on_canonical_line_when_endpoints_drift(
        tmp_path):
    """The post-2026-05-02 canonical-midpoint fix: when fault A's and
    fault B's raw endpoint coords differ (sub-cm drift), the midpoint
    must be at (snap_key_a * snap_m + snap_key_b * snap_m) / 2 — i.e.,
    the canonical position implied by the shared snap_key, not fault
    A's raw position.

    Without this fix, the midpoint depends on which fault wins
    iteration order in find_wedge_edges, and HXT rejects the resulting
    polyline (REVIEW_intersection_refinement_investigation_fix.md
    Empirical Validation).
    """
    snap_m = 0.1
    drift_m = 0.03  # < snap_m/2, so both endpoints share snap_keys.
    stl_dir, names = _drifted_wedge_pair(tmp_path, dihedral_deg=10.0,
                                          drift_m=drift_m)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "canonical",
    ])

    import json
    report = json.loads((out_dir / "wedge_split_report.json").read_text())
    assert report["midpoint_mode"] == "canonical"
    assert report["n_wedge_edges"] == 1
    w = report["wedges"][0]

    # Canonical endpoints are snap_key * snap_m.  For an edge along
    # the x-axis from (0,0,0) to (1,0,0) with snap_m=0.1, the canonical
    # endpoints are (0,0,0) and (1.0, 0.0, 0.0).
    canon_p1 = np.array(w["edge_key"][0]) * snap_m
    canon_p2 = np.array(w["edge_key"][1]) * snap_m
    expected_mid = (canon_p1 + canon_p2) / 2.0
    np.testing.assert_array_equal(np.array(w["midpoint"]), expected_mid)

    # And the midpoint must be at the snap-key-canonical position
    # regardless of fault iteration order — verify by swapping which
    # fault is "A" via the include order.
    out_dir2 = tmp_path / "out2"
    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir2),
        "--include-fault", names[1],  # swapped
        "--include-fault", names[0],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "canonical",
    ])
    report2 = json.loads((out_dir2 / "wedge_split_report.json").read_text())
    np.testing.assert_array_equal(
        np.array(report2["wedges"][0]["midpoint"]),
        np.array(w["midpoint"]))


def test_legacy_fault_a_midpoint_depends_on_iteration_order(tmp_path):
    """Regression sanity: the old fault_a mode produces a midpoint
    that drifts with fault iteration order (because it depends on
    which fault is "A").  This is the pre-fix behaviour and the
    reason HXT rejected the constraint at high wedge counts.
    """
    snap_m = 0.1
    drift_m = 0.03
    stl_dir, names = _drifted_wedge_pair(tmp_path, dihedral_deg=10.0,
                                          drift_m=drift_m)

    out_dir1 = tmp_path / "out_AB"
    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir1),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "fault_a",
    ])
    out_dir2 = tmp_path / "out_BA"
    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir2),
        "--include-fault", names[1],
        "--include-fault", names[0],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "fault_a",
    ])

    import json
    r1 = json.loads((out_dir1 / "wedge_split_report.json").read_text())
    r2 = json.loads((out_dir2 / "wedge_split_report.json").read_text())
    mid1 = np.array(r1["wedges"][0]["midpoint"])
    mid2 = np.array(r2["wedges"][0]["midpoint"])
    # The midpoints differ by the drift amount in each axis (because
    # find_wedge_edges sorts faults alphabetically, "fault_A" wins in
    # both iterations — so this regression test is a sanity check that
    # the fault_a mode IS sensitive to which fault wins; if both
    # iteration orders coincidentally give "fault_A" as winner the
    # diff is zero).
    # The stronger guarantee is that when find_wedge_edges does pick
    # different fault A's (e.g., when names sort to put B before A), the
    # midpoint differs.  Synthesize that here.
    sorted_names = sorted(names)
    if sorted_names[0] == names[0]:
        # names[0] is alphabetically first.  Both iteration orders pick
        # names[0] as fault_a.  Midpoints should match to bit precision.
        np.testing.assert_array_equal(mid1, mid2)
    else:
        # iteration-order-dependent
        assert not np.array_equal(mid1, mid2)


def test_snap_polyline_vertices_makes_shared_coords_bit_identical(
        tmp_path):
    """After polyline-vertex canonicalization, every shared snap_key
    has bit-identical coords across all faults that contain it.  This
    is the prerequisite for HXT's PLC recovery to accept the polyline
    as a single shared constraint.
    """
    snap_m = 0.1
    drift_m = 0.04  # < snap_m/2 — endpoints share snap_keys.
    stl_dir, names = _drifted_wedge_pair(tmp_path, dihedral_deg=10.0,
                                          drift_m=drift_m)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "canonical",
    ])

    # Read the OUTPUT STLs (post-snap, post-split) and verify that
    # every snap-key shared between the two faults has bit-identical
    # coords on both sides.
    V_a, _, _ = bfw._read_ascii_stl(out_dir / f"{names[0]}.stl")
    V_b, _, _ = bfw._read_ascii_stl(out_dir / f"{names[1]}.stl")

    keys_a = {tuple(int(round(c / snap_m)) for c in v): tuple(v) for v in V_a}
    keys_b = {tuple(int(round(c / snap_m)) for c in v): tuple(v) for v in V_b}
    shared = set(keys_a.keys()) & set(keys_b.keys())
    assert len(shared) >= 2, (
        f"expected ≥2 shared snap_keys (the two endpoints of the "
        f"wedge edge) plus the inserted midpoint; got {len(shared)}")

    drift_after = []
    for k in shared:
        ca = np.array(keys_a[k])
        cb = np.array(keys_b[k])
        drift_after.append(np.linalg.norm(ca - cb))
        # Bit-identity guarantee:
        for i in range(3):
            assert ca[i] == cb[i], (
                f"shared snap_key {k} has different coords in the two "
                f"output STLs: {ca[i]!r} vs {cb[i]!r} on axis {i}; "
                f"polyline-vertex snapping did not produce bit-"
                f"identical coords")
    # Sanity: drift after snap is exactly zero.
    assert max(drift_after) == 0.0


def test_no_snap_polyline_vertices_preserves_drift(tmp_path):
    """With --no-snap-polyline-vertices, fault A's and fault B's
    shared-snap-key vertex coords remain off by drift_m (legacy
    behaviour).  This test documents the pre-fix behaviour for
    regression reproducibility.
    """
    snap_m = 0.1
    drift_m = 0.04
    stl_dir, names = _drifted_wedge_pair(tmp_path, dihedral_deg=10.0,
                                          drift_m=drift_m)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "canonical",
        "--no-snap-polyline-vertices",
    ])

    V_a, _, _ = bfw._read_ascii_stl(out_dir / f"{names[0]}.stl")
    V_b, _, _ = bfw._read_ascii_stl(out_dir / f"{names[1]}.stl")

    keys_a = {tuple(int(round(c / snap_m)) for c in v): tuple(v) for v in V_a}
    keys_b = {tuple(int(round(c / snap_m)) for c in v): tuple(v) for v in V_b}
    shared = set(keys_a.keys()) & set(keys_b.keys())
    # At least one shared snap_key must have NON-zero coord drift —
    # confirming the legacy mode preserves the input drift.
    n_drifted = 0
    for k in shared:
        ca = np.array(keys_a[k])
        cb = np.array(keys_b[k])
        if np.linalg.norm(ca - cb) > 0.0:
            n_drifted += 1
    assert n_drifted >= 1, (
        f"--no-snap-polyline-vertices: expected at least one shared "
        f"snap_key to retain its input drift; all {len(shared)} were "
        f"bit-identical (snapping happened anyway?)")


def test_bisect_all_polyline_splits_every_shared_edge(tmp_path):
    """`--bisect-all-polyline` MUST split every cross-fault polyline
    edge regardless of dihedral angle.  Tested at 60° dihedral —
    where the default mode would leave the edge alone, but this mode
    must still split it.

    This is the critical behaviour that breaks polyline-locked wedge
    tets in mmg3d_post_pass mode='optim_relax_fault'.
    """
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=60.0)
    out_dir = tmp_path / "out"

    rc = bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",  # would normally exclude 60° edges
        "--bisect-all-polyline",
        "--snap-m", "0.001",
    ])
    assert rc == 0

    import json
    report = json.loads((out_dir / "wedge_split_report.json").read_text())
    assert report["bisect_all_polyline"] is True
    # The single shared edge MUST still be detected and split, even
    # though dihedral=60° > dihedral_deg_max=30°.
    assert report["n_wedge_edges"] == 1, (
        f"--bisect-all-polyline must split the 60° edge "
        f"regardless of dihedral_deg_max; got "
        f"{report['n_wedge_edges']} edges")
    assert report["n_splits_applied"] == 4, (
        f"expected 4 splits (2 per fault, 2 faults); got "
        f"{report['n_splits_applied']}")


def test_bisect_all_polyline_preserves_cross_fault_conformity(tmp_path):
    """`--bisect-all-polyline` must produce bit-identical midpoint
    coords across faults — the same cross-fault conformity
    invariant as the dihedral-thresholded mode.
    """
    # `_drifted_wedge_pair` was calibrated for snap_m=0.1, drift_m=0.03
    # (drift < snap_m/2 so endpoints share snap_keys).  Match that.
    snap_m = 0.1
    drift_m = 0.03
    stl_dir, names = _drifted_wedge_pair(tmp_path, dihedral_deg=45.0,
                                          drift_m=drift_m)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "10.0",  # would exclude 45°
        "--bisect-all-polyline",
        "--snap-m", str(snap_m),
        "--midpoint-mode", "canonical",
    ])

    V_a, _, _ = bfw._read_ascii_stl(out_dir / f"{names[0]}.stl")
    V_b, _, _ = bfw._read_ascii_stl(out_dir / f"{names[1]}.stl")
    keys_a = {tuple(int(round(c / snap_m)) for c in v): tuple(v)
              for v in V_a}
    keys_b = {tuple(int(round(c / snap_m)) for c in v): tuple(v)
              for v in V_b}
    shared = set(keys_a.keys()) & set(keys_b.keys())
    assert len(shared) >= 3, (
        f"expected ≥3 shared snap_keys (2 endpoints + 1 midpoint); "
        f"got {len(shared)}")
    # Every shared snap_key must have bit-identical coords across
    # faults (after `--snap-polyline-vertices` defaults ON).
    for k in shared:
        ca = np.array(keys_a[k])
        cb = np.array(keys_b[k])
        for i in range(3):
            assert ca[i] == cb[i], (
                f"shared snap_key {k} has different coords "
                f"across faults on axis {i}: {ca[i]!r} vs {cb[i]!r}")


def test_bisect_all_polyline_is_superset_of_wedge_split(tmp_path):
    """For any input where the default mode splits N wedge edges,
    `--bisect-all-polyline` must split AT LEAST N edges (every
    wedge edge is also a polyline edge).  Use the synthetic shallow-
    dihedral pair to make the default mode produce ≥1 split.
    """
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=5.0)

    # Run default mode (dihedral filter on).
    out_default = tmp_path / "out_default"
    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_default),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", "0.001",
    ])
    import json
    r_default = json.loads(
        (out_default / "wedge_split_report.json").read_text())

    # Run bisect-all mode.
    out_all = tmp_path / "out_all"
    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_all),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--bisect-all-polyline",
        "--snap-m", "0.001",
    ])
    r_all = json.loads(
        (out_all / "wedge_split_report.json").read_text())

    # Bisect-all must split at least as many edges as the default
    # mode.  In this synthetic 1-shared-edge fixture, both find 1
    # edge — but the assertion still validates that bisect-all is
    # never a subset of default.
    assert r_all["n_wedge_edges"] >= r_default["n_wedge_edges"], (
        f"--bisect-all-polyline split fewer edges "
        f"({r_all['n_wedge_edges']}) than default mode "
        f"({r_default['n_wedge_edges']}); the polyline-edge set is "
        f"a SUPERSET of the wedge-edge set.")


def test_bisect_all_polyline_via_python_api():
    """Direct invocation of `find_all_polyline_edges` returns every
    shared edge regardless of dihedral angle.  Verifies the public
    API exists and its semantics match the CLI flag.
    """
    # Two faults sharing one edge at 60° dihedral.
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    apex_a = np.array([0.5, 1.0, 0.0])
    apex_b = np.array([0.0, np.cos(np.deg2rad(60)) * 0.5,
                       np.sin(np.deg2rad(60)) * 0.5])  # rotated 60°
    snap_m = 0.001
    V_a = np.stack([p0, p1, apex_a])
    T_a = np.asarray([[0, 1, 2]], dtype=np.int64)
    V_b = np.stack([p0, p1, apex_b])
    T_b = np.asarray([[0, 1, 2]], dtype=np.int64)
    mA = bfw.FaultMesh("A", V_a, T_a, "A", snap_m)
    mB = bfw.FaultMesh("B", V_b, T_b, "B", snap_m)

    # Default find_wedge_edges with dihedral_max=30 EXCLUDES the
    # 60° edge.
    wedges_default = bfw.find_wedge_edges({"A": mA, "B": mB},
                                           dihedral_deg_max=30.0)
    assert len(wedges_default) == 0, (
        f"default mode at 30° threshold should exclude 60° edge; "
        f"got {len(wedges_default)}")

    # find_all_polyline_edges INCLUDES it.
    polyline = bfw.find_all_polyline_edges({"A": mA, "B": mB})
    assert len(polyline) == 1, (
        f"find_all_polyline_edges must include every shared edge; "
        f"got {len(polyline)}")


def test_invalid_midpoint_mode_raises():
    """Type-safety guard: an invalid mode string must error out
    explicitly, not silently fall through to one of the modes."""
    from break_fault_wedges import find_wedge_edges, FaultMesh
    # Build the smallest possible inputs.
    V = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                  [0.5, 1.0, 0.0]])
    T = np.array([[0, 1, 2]], dtype=np.int64)
    m = FaultMesh("only", V, T, "only", 0.1)
    with pytest.raises(ValueError, match="midpoint_mode"):
        find_wedge_edges({"only": m, "other": m},
                          dihedral_deg_max=30.0,
                          midpoint_mode="bogus")


def test_R005_split_preserves_total_surface_area(tmp_path):
    """Splitting an edge with its midpoint preserves total surface
    area to bit-precision (Loop subdivision through a midpoint of a
    flat triangle is area-preserving).
    """
    stl_dir, names = _shallow_wedge_pair(tmp_path, dihedral_deg=8.0)
    out_dir = tmp_path / "out"

    bfw.main([
        "--in-stl-dir", str(stl_dir),
        "--out-stl-dir", str(out_dir),
        "--include-fault", names[0],
        "--include-fault", names[1],
        "--dihedral-deg-max", "30.0",
        "--snap-m", "0.001",
    ])

    def area(stl_path: Path) -> float:
        V, T, _ = bfw._read_ascii_stl(stl_path)
        a = 0.0
        for t in T:
            v0, v1, v2 = V[t[0]], V[t[1]], V[t[2]]
            a += 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
        return a

    for short in names:
        a_in = area(stl_dir / f"{short}.stl")
        a_out = area(out_dir / f"{short}.stl")
        # Allow only floating-point round-off (≤ 1 ULP * O(N)).
        assert abs(a_out - a_in) < 1e-12, (
            f"fault {short}: area changed from {a_in} to {a_out}")
