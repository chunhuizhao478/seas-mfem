"""Unit tests for orient_fault_surface.py (REVIEW.md R-003).

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_orient_fault_surface.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import orient_fault_surface  # noqa: E402


def test_R003_orient_fixes_flipped_winding():
    """A pair of fault tris with same-direction shared edge: the
    second tri's winding must be flipped, recorded in stats."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)  # 1 flip
    ttags = np.array([100, 100], dtype=np.int32)
    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)
    assert stats["n_flipped"] == 1
    # Vertex set unchanged but ordering flipped.
    assert frozenset(new_tris[1]) == frozenset((1, 2, 3))
    assert tuple(new_tris[1]) != tuple(tris[1])


def test_R003_orient_idempotent():
    """If the input already has consistent winding, no flips occur."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    ttags = np.array([100, 100], dtype=np.int32)
    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)
    assert stats["n_flipped"] == 0
    np.testing.assert_array_equal(new_tris, tris)


def test_R003_orient_running_twice_is_idempotent():
    """Running the orient pass on its own output produces no new flips."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
    ttags = np.array([100, 100], dtype=np.int32)
    new_tris1, stats1 = orient_fault_surface._orient_fault(
        tris, ttags, points)
    new_tris2, stats2 = orient_fault_surface._orient_fault(
        new_tris1, ttags, points)
    assert stats2["n_flipped"] == 0
    np.testing.assert_array_equal(new_tris1, new_tris2)


def test_R003_orient_passes_through_box_tris_unchanged():
    """Non-fault triangles are passed through without modification."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                      dtype=np.float64)
    tris = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
    ttags = np.array([5, 6], dtype=np.int32)  # Box tags, not fault.
    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)
    assert stats["n_flipped"] == 0
    assert stats["n_fault_tris"] == 0
    np.testing.assert_array_equal(new_tris, tris)


def test_R003_orient_handles_multiple_components():
    """Two disconnected fault patches → 2 components in stats."""
    points = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0],
        [10, 0, 0], [11, 0, 0], [11, 1, 0],
    ], dtype=np.float64)
    tris = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int64)
    ttags = np.array([100, 100], dtype=np.int32)
    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)
    assert stats["n_components"] == 2
    assert stats["n_flipped"] == 0


def test_R003_orient_roundtrips_via_msh(tmp_path):
    """End-to-end .msh read → orient → write → read produces a mesh
    where check_13 passes."""
    import meshio
    import validate_msh

    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0.5, 0.5, -1.0], [0.5, 0.5, -2.0]],
                      dtype=np.float64)
    # Two fault tris with mismatched winding + a tet to satisfy
    # the writer's cell layout.
    fault_tris = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.int64)
    fault_tags = np.array([100, 100], dtype=np.int32)
    tets = np.array([[0, 1, 4, 5]], dtype=np.int64)
    tet_tags = np.array([10], dtype=np.int32)

    in_path = tmp_path / "in.msh"
    out_path = tmp_path / "out.msh"
    mesh = meshio.Mesh(
        points=points,
        cells=[("triangle", fault_tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [fault_tags, tet_tags],
                   "gmsh:geometrical": [fault_tags, tet_tags]},
    )
    meshio.write(in_path, mesh, file_format="gmsh22", binary=False)

    rc = orient_fault_surface.main(
        ["--in-msh", str(in_path), "--out-msh", str(out_path)])
    assert rc == 0

    # Re-read and check_13 should pass.
    pts2, tris2, ttags2, tets2, _ = orient_fault_surface._read_mesh(
        out_path)
    r = validate_msh.check_13_fault_orientation(pts2, tris2, ttags2)
    assert r.passed


def test_R003_orient_does_not_propagate_through_x_junction():
    """REGRESSION: at a non-manifold (count != 2) edge — typical at a
    fault-fault X-junction in branching SAFS geometry — orient must
    NOT propagate orientation across.  Otherwise two crossing fault
    sheets get forced into the same winding convention and one of
    them ends up back-faced (visible as triangular gaps in ParaView
    even though every topology check passes).

    Setup: two perpendicular fault sheets sharing edge (0, 1).
      Sheet A: tris (0, 1, 2) and (0, 1, 3) on one side  (z=0 plane)
      Sheet B: tris (0, 1, 4) and (0, 1, 5) on another side (y=0 plane)
    The shared edge (0, 1) has count=4 — an X-junction.

    Sheet A's tris have OPPOSITE windings on (0, 1) — they are the two
    halves of one manifold sheet — so internally consistent.
    Sheet B's tris also have opposite windings — internally consistent.
    But sheet A's (0, 1) direction differs from sheet B's: forcing
    propagation would flip half of sheet B (or A).
    """
    points = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],            # 0, 1 (shared)
        [0.5, 1.0, 0.0], [0.5, -1.0, 0.0],           # 2, 3 in z=0 plane
        [0.5, 0.0, 1.0], [0.5, 0.0, -1.0],           # 4, 5 in y=0 plane
    ], dtype=np.float64)

    # Sheet A: (0,1,2) traverses (0,1); (0,1,3) reversed traverses (1,0)
    # Sheet B: (0,1,4) traverses (0,1); (0,1,5) reversed traverses (1,0)
    # Both sheets are internally consistent (0 flips needed within
    # each sheet via manifold-only BFS).
    tris = np.array([
        [0, 1, 2],   # sheet A, (0,1)
        [1, 0, 3],   # sheet A, (1,0)  - opposite, manifold-consistent
        [0, 1, 4],   # sheet B, (0,1)
        [1, 0, 5],   # sheet B, (1,0)  - opposite, manifold-consistent
    ], dtype=np.int64)
    ttags = np.array([100, 100, 100, 100], dtype=np.int32)

    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)

    # Manifold-only BFS sees count==4 at edge {0,1} as a barrier;
    # each sheet's two tris are connected only through their
    # respective {0,2}/{1,2}/{0,3}/{1,3} manifold edges within the
    # shared-edge subgraph — but those edges have count=1 here,
    # so each tri is its own component.
    # Either way, NO flips should be needed and components >= 2.
    assert stats["n_flipped"] == 0, (
        f"manifold-only BFS must not flip across non-manifold "
        f"edges; got {stats['n_flipped']} flips")
    assert stats["n_nonmanifold_edges"] >= 1, (
        "the X-junction edge {0,1} must be reported as non-manifold")
    # All 4 tris remain in their original winding.
    assert (new_tris == tris).all(), (
        "no tri should have been flipped at the X-junction")


def test_R003_orient_repairs_pair_inconsistency_at_x_junction():
    """REGRESSION: at a count==4 X-junction, two crossing fault sheets
    each have their own pair of tris.  Each pair MUST be propagated
    INTERNALLY (one flip if they disagree on the edge direction), but
    the two pairs must NOT cross-contaminate each other.  Earlier
    versions of this code either over-propagated (12k false flips at
    SAFS X-junctions) or under-propagated (legitimate within-pair
    flips skipped).
    """
    points = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],            # 0, 1 (shared)
        [0.5, 1.0, 0.0], [0.5, -1.0, 0.0],           # 2, 3  (sheet A)
        [0.5, 0.0, 1.0], [0.5, 0.0, -1.0],           # 4, 5  (sheet B)
    ], dtype=np.float64)

    # Sheet A: tris (0,1,2) and (0,1,3) — BOTH traverse (0,1).
    # That is an internal manifold-pair WINDING ERROR within sheet A:
    # one of them should be flipped to (1,0,3).
    # Sheet B: tris (0,1,4) and (1,0,5) — already correct (opposite).
    tris = np.array([
        [0, 1, 2],   # sheet A, (0,1)  - same dir as sheet A's other tri
        [0, 1, 3],   # sheet A, (0,1)  - WRONG; should be (1,0,3)
        [0, 1, 4],   # sheet B, (0,1)
        [1, 0, 5],   # sheet B, (1,0) - already opposite
    ], dtype=np.int64)
    ttags = np.array([100, 100, 100, 100], dtype=np.int32)

    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)

    # Sheet A's tri 1 must be flipped (winding (0,1,3) → (0,3,1))
    # Sheet B must remain untouched.
    assert stats["n_flipped"] == 1, (
        f"expected exactly 1 flip (sheet A's bad tri); got "
        f"{stats['n_flipped']}")
    # New row 1 should now traverse (1,0) on edge {0,1} — i.e.
    # original (0,1,3) became (0,3,1).
    assert tuple(new_tris[1]) == (0, 3, 1), (
        f"sheet A's tri 1 should have winding (0,3,1); got "
        f"{tuple(new_tris[1])}")
    # Sheet B's tris untouched.
    assert tuple(new_tris[2]) == (0, 1, 4)
    assert tuple(new_tris[3]) == (1, 0, 5)


def test_R003_orient_handles_t_junction_count3():
    """At a count==3 T-junction, the 3 tris should partition into
    1 dihedral pair (the through-going fault) + 1 terminator
    (the fault that ends here).  Orientation propagates within the
    pair; the terminator gets no propagation across this edge.
    """
    points = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],            # 0, 1 (shared)
        [0.5, 1.0, 0.0], [0.5, -1.0, 0.0],           # 2, 3  through-going
        [0.5, 0.0, 1.0],                              # 4    terminator
    ], dtype=np.float64)

    # Through-going pair: (0,1,2) and (0,1,3).  Same direction on
    # edge {0,1} — needs flip.
    # Terminator: (0,1,4) — orthogonal sheet ending here, no
    # antiparallel partner among the remaining tris.
    tris = np.array([
        [0, 1, 2],
        [0, 1, 3],   # SAME direction — must be flipped
        [0, 1, 4],   # terminator (orthogonal sheet)
    ], dtype=np.int64)
    ttags = np.array([100, 100, 100], dtype=np.int32)

    new_tris, stats = orient_fault_surface._orient_fault(
        tris, ttags, points)

    assert stats["n_flipped"] == 1, (
        f"expected 1 flip on through-going pair; got "
        f"{stats['n_flipped']}")
    # Pair member: tri 1 must now traverse (1,0)
    assert tuple(new_tris[1]) == (0, 3, 1), (
        f"through-going pair's tri 1 should be (0,3,1); got "
        f"{tuple(new_tris[1])}")
    # Terminator untouched
    assert tuple(new_tris[2]) == (0, 1, 4)
    assert stats["n_nonmanifold_edges"] >= 1
