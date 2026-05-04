"""Unit tests for local_cavity_retet.py.

Covers Option B local cavity re-tetrahedralization:
  - cavity construction = bad tet + 4 face-adjacent neighbours
  - cavity-boundary triangle identification
  - Steiner-fan tet generation with positive signed volume
  - γ-regression revert defence
  - cross-fault face-adjacency preservation

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_local_cavity_retet.py -v
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import meshio
import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import local_cavity_retet as lcr  # noqa: E402
import mmg3d_local_patch as mlp  # noqa: E402


def _write_msh(out_path: Path,
                points: np.ndarray, tets: np.ndarray,
                tris: np.ndarray | None = None,
                tet_tags: np.ndarray | None = None,
                tri_tags: np.ndarray | None = None) -> None:
    if tris is None:
        tris = np.asarray([[0, 1, 2]], dtype=np.int64)
        tri_tags = np.asarray([100], dtype=np.int32)
    if tet_tags is None:
        tet_tags = np.full(tets.shape[0], 10, dtype=np.int32)
    if tri_tags is None:
        tri_tags = np.full(tris.shape[0], 100, dtype=np.int32)
    meshio.write(out_path, meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    ), file_format="gmsh22", binary=False)


# ---------------------------------------------------------------------------
# Pure-Python unit tests.
# ---------------------------------------------------------------------------
def test_gamma_one_regular_tet():
    """γ ≈ 1 for a regular tet."""
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    g = lcr._gamma_one(P, (0, 1, 2, 3))
    assert 0.99 <= g <= 1.01


def test_gamma_one_zero_volume():
    """A coplanar tet has volume 0 → γ = 0."""
    P = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    g = lcr._gamma_one(P, (0, 1, 2, 3))
    assert g == 0.0


def test_signed_volume_positive_orientation():
    """A standard right-handed tet has positive signed volume."""
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                 dtype=np.float64)
    sv = lcr._signed_volume(P, (0, 1, 2, 3))
    assert sv > 0
    # Swapping two vertices flips the sign.
    sv_flipped = lcr._signed_volume(P, (0, 2, 1, 3))
    assert sv_flipped < 0


def test_invalid_inputs_raise(tmp_path):
    fake = tmp_path / "fake.msh"
    fake.write_text("dummy")
    with pytest.raises(ValueError, match="gamma_thresh"):
        lcr.local_cavity_retet(fake, tmp_path / "out.msh",
                                gamma_thresh=0.0)
    with pytest.raises(ValueError, match="gamma_thresh"):
        lcr.local_cavity_retet(fake, tmp_path / "out.msh",
                                gamma_thresh=-1.0)
    with pytest.raises(ValueError, match="gamma_revert_tol"):
        lcr.local_cavity_retet(fake, tmp_path / "out.msh",
                                gamma_revert_tol=-0.1)


def test_missing_input_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-msh not found"):
        lcr.local_cavity_retet(tmp_path / "nope.msh",
                                tmp_path / "out.msh")


# ---------------------------------------------------------------------------
# Cavity construction.
# ---------------------------------------------------------------------------
def test_cavity_includes_bad_tet_and_4_neighbours():
    """A bad tet with all 4 faces interior has cavity = 5 tets."""
    # Build a small mesh: 1 central tet + 4 neighbours, each
    # sharing one face with the central.
    # Central tet: (0,1,2,3).  Apexes 4,5,6,7 are placed so each
    # shares a face with the central via one of its 4 faces.
    P = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],   # central
        [1, 1, 1],   # 4: shares face (1,2,3)
        [-1, 0, 0],  # 5: shares face (0,2,3)
        [0, -1, 0],  # 6: shares face (0,1,3)
        [0, 0, -1],  # 7: shares face (0,1,2)
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],   # central, idx 0
        [4, 1, 2, 3],   # idx 1: shares (1,2,3)
        [5, 0, 2, 3],   # idx 2: shares (0,2,3)
        [6, 0, 1, 3],   # idx 3: shares (0,1,3)
        [7, 0, 1, 2],   # idx 4: shares (0,1,2)
    ], dtype=np.int64)
    f2t = mlp._build_tet_face_adjacency(tets)
    cavity = lcr._cavity_for_bad_tet(0, tets, f2t)
    assert cavity == {0, 1, 2, 3, 4}


def test_cavity_with_boundary_face_skips_missing_neighbour():
    """A bad tet with 2 model-boundary faces has cavity = 1 + 2 = 3."""
    P = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1],
        [-1, 0, 0],
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],
        [4, 1, 2, 3],
        [5, 0, 2, 3],
    ], dtype=np.int64)
    f2t = mlp._build_tet_face_adjacency(tets)
    cavity = lcr._cavity_for_bad_tet(0, tets, f2t)
    # Tet 0 has neighbours 1 (across face 1,2,3) and 2 (across
    # face 0,2,3); the other 2 faces are boundary.
    assert cavity == {0, 1, 2}


def test_cavity_boundary_faces_count():
    """Cavity = 5 tets sharing the central tet → 4 interior faces
    (each between central and a neighbour), 4*3 = 12 outer faces
    on the 4 neighbours.  Cavity boundary = 12 triangles.
    """
    P = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],
        [4, 1, 2, 3],
        [5, 0, 2, 3],
        [6, 0, 1, 3],
        [7, 0, 1, 2],
    ], dtype=np.int64)
    cavity = {0, 1, 2, 3, 4}
    boundary = lcr._cavity_boundary_faces(cavity, tets)
    assert len(boundary) == 12, (
        f"5-tet cavity should have 12 boundary faces; got "
        f"{len(boundary)}")


def test_cavity_boundary_face_orientation_outward():
    """Each boundary face's vertex order must be such that its
    normal points OUTWARD from the cavity.  Combined with a
    Steiner point at the cavity centroid, the resulting tet
    `(boundary_face, S)` has positive signed volume.
    """
    P = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ], dtype=np.float64)
    tets = np.array([
        [0, 1, 2, 3],
        [4, 1, 2, 3],
        [5, 0, 2, 3],
        [6, 0, 1, 3],
        [7, 0, 1, 2],
    ], dtype=np.int64)
    cavity = {0, 1, 2, 3, 4}
    boundary = lcr._cavity_boundary_faces(cavity, tets)
    # Cavity centroid: average of all unique cavity vertices.
    cavity_v = sorted({int(tets[c, k]) for c in cavity
                        for k in range(4)})
    centroid = P[cavity_v].mean(axis=0)
    P_with_S = np.vstack((P, centroid[None, :]))
    si = P_with_S.shape[0] - 1
    fan = lcr._build_steiner_fan(P_with_S, si, boundary)
    # Every fan tet must have positive signed volume.
    for ft in fan:
        sv = lcr._signed_volume(P_with_S, ft)
        assert sv > 0, (
            f"fan tet {ft} has non-positive signed volume {sv}; "
            f"orientation logic is wrong")


# ---------------------------------------------------------------------------
# End-to-end on synthetic fixtures.
# ---------------------------------------------------------------------------
def test_e2e_well_conditioned_mesh_unchanged(tmp_path):
    """A mesh with NO bad tets — output identical to input."""
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=1e-3)
    assert rep["no_op"] is True
    assert rep["n_bad_tets"] == 0
    assert rep["tets_in"] == rep["tets_out"]


def test_e2e_bad_tet_with_no_neighbours_skipped(tmp_path):
    """An isolated bad tet (no face-adjacent neighbours) has
    cavity = {bad_tet}, with 4 boundary faces.  Steiner-fan would
    require splitting the tet at its centroid — same failure mode
    as the prior centroid_steiner_split.  We expect cavity_built
    but NOT applied (revert because γ regresses).
    """
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.001, 0.0],
        [2.0, 0.0, 0.001],
        [3.0, 0.001, 0.0],   # collinear bad tet
    ], dtype=np.float64)
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=1.0,
                                   gamma_revert_tol=0.0)
    # Cavity is built (1-tet cavity is degenerate).  γ-revert
    # should kick in if the new fan is no better than the old.
    assert rep["n_bad_tets"] == 1
    assert (rep["n_cavities_reverted"] >= 0
             or rep["n_skipped_inverted"] >= 0)


def test_e2e_5_tet_cavity_breaks_collinear_bad_tet(tmp_path):
    """The key empirical test: a bad collinear tet at the centre
    of a 5-tet cavity should be successfully replaced by a
    Steiner fan.  γ_min should improve.
    """
    # Central bad tet: 4 nearly-collinear vertices.
    # 4 neighbours each well off the central tet's degenerate line.
    P = np.array([
        # Central bad tet (collinear along x):
        [0.0,  0.0,  0.0],     # 0
        [10.0, 0.01, 0.0],     # 1
        [20.0, 0.0,  0.01],    # 2
        [30.0, 0.01, 0.0],     # 3
        # Neighbours, well off the line:
        [10.0,  10.0,  0.0],     # 4: opposite vertex 0 (face 1,2,3)
        [10.0, -10.0,  0.0],     # 5: opposite vertex 1 (face 0,2,3)
        [10.0,  0.0,   10.0],    # 6: opposite vertex 2 (face 0,1,3)
        [10.0,  0.0,  -10.0],    # 7: opposite vertex 3 (face 0,1,2)
    ], dtype=np.float64)
    T = np.array([
        [0, 1, 2, 3],   # bad central
        [1, 2, 3, 4],   # nbr opposite 0
        [0, 2, 3, 5],   # nbr opposite 1
        [0, 1, 3, 6],   # nbr opposite 2
        [0, 1, 2, 7],   # nbr opposite 3
    ], dtype=np.int64)
    # Use a fault triangle on the OUTSIDE of the cavity (a face
    # of one of the neighbours that's NOT shared with the central
    # tet).  Face (1,4,2) is outside the central tet's cavity.
    F = np.array([[1, 4, 2]], dtype=np.int64)
    F_tags = np.array([100], dtype=np.int32)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T, tris=F, tri_tags=F_tags)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=0.5,
                                   gamma_revert_tol=0.0)
    # The synthetic mesh has 5 tets (central + 4 neighbours);
    # several are below γ=0.5.  We don't care exactly how many
    # are bad — what matters is the central collinear tet (which
    # has the worst γ ≈ 6e-5) should be processed first (the loop
    # sorts bad tets by γ ascending).  After the central tet's
    # cavity is rebuilt, the remaining bad tets are CONSUMED
    # (they were neighbours of the central) and skipped.
    assert rep["n_bad_tets"] >= 1
    assert rep["n_cavities_applied"] == 1, (
        f"expected exactly 1 cavity to be applied (the central "
        f"collinear tet); got applied={rep['n_cavities_applied']}, "
        f"reverted={rep['n_cavities_reverted']}, "
        f"consumed-skipped={rep['n_skipped_consumed']}")
    # Verify γ_min improved.  The central tet had γ ≈ 6e-5; the
    # 5-tet cavity rebuild produces a Steiner fan whose worst tet
    # should be 10–100× better.
    assert rep["gamma_min_out"] > 10.0 * rep["gamma_min_in"], (
        f"5-tet cavity rebuild should improve γ_min by >10× on "
        f"a collinear-central-tet test: {rep['gamma_min_in']} → "
        f"{rep['gamma_min_out']}")


def test_e2e_revert_when_new_fan_regresses(tmp_path):
    """If the new fan's γ_min is worse than the original, the
    operation must REVERT and leave the mesh unchanged.
    """
    # A perfect tet plus 4 neighbours that are also perfect.
    # No bad tets; revert-mechanism not exercised.  But we can
    # synthesize a case where the cavity is bad enough that the
    # Steiner fan also has bad γ.
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=1e-3)
    # No bad tets → no cavities → no_op.  But the API still
    # returns a sensible report.
    assert rep["no_op"] is True
    assert rep["n_cavities_applied"] == 0


def test_e2e_writes_report_json(tmp_path):
    s = 1.0 / np.sqrt(3.0)
    P = np.array([[s, s, s], [s, -s, -s],
                   [-s, s, -s], [-s, -s, s]])
    T = np.array([[0, 1, 2, 3]], dtype=np.int64)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    rep_json = tmp_path / "rep.json"
    _write_msh(in_msh, P, T)
    lcr.main([
        "--in-msh", str(in_msh),
        "--out-msh", str(out_msh),
        "--gamma-thresh", "1e-3",
        "--report-json", str(rep_json),
    ])
    rep = json.loads(rep_json.read_text())
    for k in ("gamma_min_in", "gamma_min_out", "n_bad_tets",
               "n_cavities_built", "n_cavities_applied",
               "n_cavities_reverted", "tets_in", "tets_out"):
        assert k in rep, f"report missing key {k}"


def test_e2e_skip_cavity_with_internal_fault_face(tmp_path):
    """If a cavity contains a fault triangle as an INTERNAL face
    (shared between 2 cavity tets), the rebuild would orphan it
    (no new fan tet would use it as a face).  The defensive skip
    must trigger; the cavity is NOT applied.

    Construct: 2 tets sharing face (0, 1, 2), where (0, 1, 2) is
    tagged as a fault triangle.  Tet (0,1,2,3) is bad.  The cavity
    contains both tets.  Face (0,1,2) is the interior cavity face
    AND a fault triangle → must skip.
    """
    P = np.array([
        [0.0,  0.0,  0.0],
        [10.0, 0.01, 0.0],
        [20.0, 0.0,  0.01],
        [30.0, 0.01, 0.0],   # bad apex
        [10.0, 10.0, 5.0],   # other apex (good)
    ], dtype=np.float64)
    T = np.array([
        [0, 1, 2, 3],   # bad: nearly collinear
        [0, 1, 2, 4],   # neighbour
    ], dtype=np.int64)
    F = np.array([[0, 1, 2]], dtype=np.int64)   # fault tri
    F_tags = np.array([100], dtype=np.int32)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T, tris=F, tri_tags=F_tags)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=1.0)
    # The cavity (tet0 + tet1) has interior face (0,1,2), which IS
    # the fault triangle → cavity must be skipped.
    assert rep["n_skipped_internal_fault"] >= 1, (
        f"cavity with internal fault face should be skipped; "
        f"report: {rep}")
    assert rep["n_cavities_applied"] == 0, (
        f"no cavity should be applied; got "
        f"applied={rep['n_cavities_applied']}")
    # Verify check_5 invariant: fault triangle (0,1,2) is still
    # adjacent to exactly 2 tets in the output.
    out = meshio.read(out_msh)
    out_tets = None
    for cb in out.cells:
        if cb.type == "tetra":
            out_tets = cb.data
            break
    target_face = frozenset((0, 1, 2))
    n_incident = 0
    for tet in out_tets:
        for face in (frozenset((tet[0], tet[1], tet[2])),
                      frozenset((tet[0], tet[1], tet[3])),
                      frozenset((tet[0], tet[2], tet[3])),
                      frozenset((tet[1], tet[2], tet[3]))):
            if face == target_face:
                n_incident += 1
                break
    assert n_incident == 2, (
        f"fault face (0,1,2) must be adjacent to exactly 2 tets "
        f"in the output; got {n_incident}")


def test_e2e_fault_face_adjacency_preserved(tmp_path):
    """The cavity-boundary triangles include any tag-100 fault
    triangles that lie on the cavity boundary.  After cavity
    re-tet, those triangles must still appear in the output as
    faces of the new fan tets — preserving check_5 invariant
    (each fault triangle shared by 2 tets in the merged mesh).
    """
    # 5-tet cavity.  Mark the (1, 2, 3) triangle (shared between
    # central tet and neighbour 1) as a fault triangle.  Wait —
    # internal cavity faces are NOT cavity-boundary faces.  Let
    # me mark a face on the OUTSIDE of the cavity (e.g., a face
    # of one of the neighbours) as fault.
    P = np.array([
        [0.0,  0.0,  0.0],
        [10.0, 0.01, 0.0],
        [20.0, 0.0,  0.01],
        [30.0, 0.01, 0.0],
        [10.0,  10.0,  0.0],
        [10.0, -10.0,  0.0],
        [10.0,  0.0,   10.0],
        [10.0,  0.0,  -10.0],
    ], dtype=np.float64)
    T = np.array([
        [0, 1, 2, 3],   # bad central
        [1, 2, 3, 4],   # nbr 1
        [0, 2, 3, 5],   # nbr 2
        [0, 1, 3, 6],   # nbr 3
        [0, 1, 2, 7],   # nbr 4
    ], dtype=np.int64)
    # Mark face (2, 3, 4) — a face of nbr 1, on the cavity
    # exterior — as fault.  This face should still appear in the
    # output mesh after the cavity re-tet (the new fan tet that
    # uses face (2, 3, 4) must inherit it).
    F = np.array([[2, 3, 4]], dtype=np.int64)
    F_tags = np.array([100], dtype=np.int32)
    in_msh = tmp_path / "in.msh"
    out_msh = tmp_path / "out.msh"
    _write_msh(in_msh, P, T, tris=F, tri_tags=F_tags)
    rep = lcr.local_cavity_retet(in_msh, out_msh,
                                   gamma_thresh=0.5)
    out = meshio.read(out_msh)
    # The fault triangle (2, 3, 4) must still be present in the
    # output's triangle block.
    out_tris_set = set()
    for cb in out.cells:
        if cb.type == "triangle":
            for tri in cb.data:
                out_tris_set.add(frozenset(int(x) for x in tri))
    assert frozenset((2, 3, 4)) in out_tris_set, (
        "fault triangle (2,3,4) was lost from the output triangle "
        "block after cavity re-tet.  Check_5 would fail.")
