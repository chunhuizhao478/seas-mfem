"""test_subdivide_fault_stl.py — invariants of the 1->4 midpoint subdivision.

The whole reason to subdivide (rather than CGAL-remesh) to reach ~215 m is that
a midpoint split is *exactly* shape-preserving and geometry-preserving, so the
`_triq` fault's 26.59 deg worst min-angle survives unchanged while the edge
length halves.  These tests pin that contract on a tiny deterministic mesh
(fast, no file dependency) and, when present, on the real produced STL.

Run:
    conda activate pythonenv
    pytest -q test_subdivide_fault_stl.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import meshio
import remesh_fault_stl as R
from subdivide_fault_stl import subdivide_1to4, _max_edge_incidence


# A two-triangle patch with: one edge (v0-v1) on z = 0 (a "trace" edge), one
# shared interior edge (v0-v2, 2-incident), and non-degenerate shape so angle
# checks are meaningful.  z <= 0 everywhere (the mesher's free-surface side).
PTS = np.array([
    [0.0, 0.0, 0.0],     # v0  on z=0
    [1.0, 0.0, 0.0],     # v1  on z=0
    [1.0, 1.0, -1.0],    # v2
    [0.0, 1.0, -1.0],    # v3
])
TRIS = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)


def _z0_trace_count(pts, tris, tol=R.EPS_Z):
    return R.z0_trace_segment_count(np.asarray(pts, float),
                                    np.asarray(tris, np.int64), tol)


def test_counts_tris_x4_verts_plus_unique_edges():
    """tris -> 4x; new verts = old verts + #unique_edges (midpoints shared)."""
    p2, t2 = subdivide_1to4(PTS, TRIS)
    assert t2.shape[0] == 4 * TRIS.shape[0]
    # input unique edges: 01,12,02,23,03 -> 5 (the shared 0-2 counted once)
    assert p2.shape[0] == PTS.shape[0] + 5


def test_midpoint_shared_on_interior_edge_is_manifold():
    """The shared interior edge's midpoint is welded (one vertex, not two), so
    the surface stays 2-manifold (max edge incidence == 2)."""
    p2, t2 = subdivide_1to4(PTS, TRIS)
    assert _max_edge_incidence(t2) == 2
    # exactly one new vertex sits at the 0-2 edge midpoint (0.5,0.5,-0.5)
    mid02 = np.array([0.5, 0.5, -0.5])
    n_at_mid = int((np.linalg.norm(p2 - mid02, axis=1) < 1e-12).sum())
    assert n_at_mid == 1, f"interior-edge midpoint not welded (found {n_at_mid})"


def test_angles_preserved_exactly():
    """1->4 midpoint split yields 4 triangles each similar to the parent, so the
    angle distribution is unchanged (this is the property CGAL-to-250 lacks)."""
    a_in, q_in = R.tri_metrics(PTS, TRIS)
    p2, t2 = subdivide_1to4(PTS, TRIS)
    a_out, q_out = R.tri_metrics(p2, t2)
    assert abs(float(a_out.min()) - float(a_in.min())) < 1e-9
    assert abs(float(q_out.min()) - float(q_in.min())) < 1e-9
    # every child angle must already exist in the parent's angle set
    parent_angles = np.sort(np.unique(np.round(a_in, 6)))
    child_angles = np.sort(np.unique(np.round(a_out, 6)))
    assert np.allclose(parent_angles, child_angles, atol=1e-6)


def test_z0_trace_doubles_and_stays_on_plane():
    """z=0 trace edges double (1 -> 2) and every trace midpoint is exactly on
    z = 0 (so run_z0cut_meshing.py's surface-rupture invariant is preserved)."""
    assert _z0_trace_count(PTS, TRIS) == 1
    p2, t2 = subdivide_1to4(PTS, TRIS)
    assert _z0_trace_count(p2, t2) == 2
    assert float(p2[:, 2].max()) <= R.EPS_Z          # nothing lifted above 0
    # the new midpoint of the z=0 edge (0-1) is bit-exactly on z=0
    mid01 = np.array([0.5, 0.0, 0.0])
    assert (np.linalg.norm(p2 - mid01, axis=1) < 1e-12).any()


def test_geometry_preserved_planar_patch():
    """All new vertices are midpoints of existing edges, hence lie on the input
    facets: a planar input stays planar (zero geometric drift)."""
    # planar triangle in the plane x + y + z = 0
    pts = np.array([[0.0, 0.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 0.0, 1.0]])
    tris = np.array([[0, 1, 2]], dtype=np.int64)
    p2, _ = subdivide_1to4(pts, tris)
    resid = p2[:, 0] + p2[:, 1] + p2[:, 2]           # plane eq residual
    assert np.allclose(resid, 0.0, atol=1e-12)


def test_levels_two_is_x16():
    """--levels 2 (applied as two passes) -> tris x16, edges halved twice."""
    p1, t1 = subdivide_1to4(PTS, TRIS)
    p2, t2 = subdivide_1to4(p1, t1)
    assert t2.shape[0] == 16 * TRIS.shape[0]
    assert _max_edge_incidence(t2) == 2


def test_rejects_bad_shape():
    """Non-(M,3) connectivity is a hard error, not a silent reshape."""
    with pytest.raises(ValueError):
        subdivide_1to4(PTS, np.array([[0, 1, 2, 3]], dtype=np.int64))


# --- tie-in to the real produced artifacts (skip if not yet built) -----------
STL_TRIQ = HERE / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut_triq.stl"
STL_SUBDIV = HERE / "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_250m_clean_clip_nwcut_triqsubdiv.stl"


def test_real_subdiv_stl_halves_edge_preserves_shape():
    """The produced ~215 m STL must be 4x the _triq tri count, ~half the median
    edge, the SAME worst min-angle, and reach z = 0 as a single component."""
    if not (STL_TRIQ.is_file() and STL_SUBDIV.is_file()):
        pytest.skip("triq / triqsubdiv STL not present (run subdivide first)")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pin, tin = R._read_triangles(meshio.read(str(STL_TRIQ)))
        pout, tout = R._read_triangles(meshio.read(str(STL_SUBDIV)))
    a_in, _ = R.tri_metrics(pin, tin)
    a_out, _ = R.tri_metrics(pout, tout)
    assert tout.shape[0] == 4 * tin.shape[0]
    r = R.median_edge_length(pout, tout) / R.median_edge_length(pin, tin)
    assert 0.45 < r < 0.55, f"median edge ratio {r:.3f} != ~0.5"
    assert abs(float(a_out.min()) - float(a_in.min())) < 1e-3
    assert float(a_out.min()) >= 25.0
    assert float(pout[:, 2].max()) <= R.DEFAULT_Z_SNAP_TOL
    assert R.n_connected_components(pout, tout) == 1
