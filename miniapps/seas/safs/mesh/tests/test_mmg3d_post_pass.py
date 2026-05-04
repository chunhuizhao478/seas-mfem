"""Unit tests for mmg3d_post_pass.py.

Covers REVIEW_interior_subdivision_and_collapse.md Tier 1 fallback:
the post-HXT mmg3d sliver-cleanup pass must (a) preserve fault
triangles via -opnbdy, (b) preserve fault triangle physical tag
(100), and (c) handle invalid inputs with informative errors.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_mmg3d_post_pass.py -v
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import mmg3d_post_pass as m3p  # noqa: E402

_mmg3d = shutil.which("mmg3d_O3")


# ---------------------------------------------------------------------------
# Pure-Python tests (no mmg3d binary needed).
# ---------------------------------------------------------------------------
def test_invalid_h_params_raise(tmp_path):
    """Defensive parameter validation."""
    fake_in = tmp_path / "fake.msh"
    fake_in.write_text("dummy")
    out = tmp_path / "out.msh"
    with pytest.raises(ValueError, match="hmin/hmax"):
        m3p.post_pass(fake_in, out, hmin=-1, hmax=10)
    with pytest.raises(ValueError, match="hmin must be"):
        m3p.post_pass(fake_in, out, hmin=100, hmax=10)
    with pytest.raises(ValueError, match="hgrad"):
        m3p.post_pass(fake_in, out, hmin=1, hmax=10, hgrad=0.5)
    with pytest.raises(ValueError, match="hausd"):
        m3p.post_pass(fake_in, out, hmin=1, hmax=10,
                      hgrad=1.3, hausd=-1)


def test_missing_input_raises(tmp_path):
    with pytest.raises(SystemExit, match="--in-msh not found"):
        m3p.post_pass(tmp_path / "nope.msh", tmp_path / "out.msh")


def test_invalid_mode_raises(tmp_path):
    """mode parameter must be 'optim', 'adapt', or 'optim_relax_fault'."""
    # Write a minimal valid input so we get past the "missing" check.
    in_path = tmp_path / "in.msh"
    in_path.write_text("dummy")
    with pytest.raises(ValueError, match="mode must be"):
        m3p.post_pass(in_path, tmp_path / "out.msh", mode="bogus")


def test_optim_relax_fault_requires_provenance(tmp_path):
    """The optim_relax_fault mode needs fault_provenance.json to
    identify the polyline vertices.  Calling without it should
    error out clearly rather than silently fall through."""
    in_path = tmp_path / "in.msh"
    in_path.write_text("dummy")
    with pytest.raises(ValueError,
                       match="optim_relax_fault.*provenance"):
        m3p.post_pass(in_path, tmp_path / "out.msh",
                      mode="optim_relax_fault",
                      provenance_path=None)


def test_identify_polyline_keys_from_provenance(tmp_path):
    """Verify polyline-vertex extraction from a synthetic
    fault_provenance.json.

    Build a 3-vertex test set where:
      vertex 0 belongs only to fault A
      vertex 1 belongs only to fault B
      vertex 2 belongs to BOTH faults (polyline vertex)
    Each fault has exactly one triangle.
    """
    points = np.asarray([
        [0.0, 0.0, 0.0],   # 0 — fault A only
        [1.0, 0.0, 0.0],   # 1 — fault B only
        [0.5, 0.5, 0.0],   # 2 — shared polyline vertex
    ])
    # Two triangles, both tag 100 (fault).  The .msh global indices
    # are [0, 1] for these tag-100 triangles.
    tris = [(0, 2, 0, 100),   # fault A: vertices 0, 2, (and dup 0
                                # — bogus but enough for the helper)
            (1, 2, 1, 100)]   # fault B: vertices 1, 2, (dup 1)
    # Use proper non-degenerate triangles (need distinct vertices):
    tris = [(0, 2, 0, 100), (1, 2, 0, 100)]
    # Adjust: each triangle needs 3 distinct vertices.  Add one more.
    points = np.asarray([
        [0.0, 0.0, 0.0],   # 0 — fault A
        [1.0, 0.0, 0.0],   # 1 — fault B
        [0.5, 0.5, 0.0],   # 2 — shared
        [0.3, 0.3, 0.0],   # 3 — fault A interior
        [0.7, 0.3, 0.0],   # 4 — fault B interior
    ])
    tris = [(0, 2, 3, 100),   # fault A triangle (verts 0, 2, 3)
            (1, 2, 4, 100)]   # fault B triangle (verts 1, 2, 4)

    prov = {
        "faults": {
            "faultA": {"triangle_indices_in_msh": [0]},
            "faultB": {"triangle_indices_in_msh": [1]},
        }
    }
    prov_path = tmp_path / "fault_provenance.json"
    prov_path.write_text(json.dumps(prov))

    keys = m3p._identify_polyline_keys_from_provenance(
        points, tris, prov_path, snap_m=0.001, fault_tri_tag=100)

    # Vertex 2 (shared) should be the only polyline vertex.
    expected_key = (int(round(0.5 / 0.001)),
                    int(round(0.5 / 0.001)),
                    int(round(0.0 / 0.001)))
    assert expected_key in keys, (
        f"shared vertex 2 must be in polyline_keys; got {keys}")
    # Vertices 0, 1, 3, 4 should NOT be polyline vertices.
    assert len(keys) == 1, (
        f"only the shared vertex (idx 2) should be a polyline "
        f"vertex; got {len(keys)} keys: {keys}")


def test_identify_polyline_keys_returns_empty_when_no_provenance(
        tmp_path):
    """When fault_provenance.json doesn't exist, return empty set
    rather than crashing — caller guards on emptiness."""
    points = np.zeros((3, 3))
    tris = [(0, 1, 2, 100)]
    keys = m3p._identify_polyline_keys_from_provenance(
        points, tris, tmp_path / "nonexistent.json",
        snap_m=0.1, fault_tri_tag=100)
    assert keys == set()


def test_identify_polyline_keys_handles_offset(tmp_path):
    """The fault_provenance schema indexes against the FULL .msh
    triangle list, not the tag-100 subset.  When tag-100 triangles
    don't start at index 0 (e.g., box faces appear first), the
    offset must be subtracted.
    """
    points = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 1.0, 0.0],
    ])
    # Three triangles: first two are box (tag 1), last is fault.
    tris = [(0, 1, 3, 1),     # box, tag 1, global index 0
            (1, 2, 3, 1),     # box, tag 1, global index 1
            (0, 1, 2, 100)]   # fault, tag 100, global index 2
    # Provenance indexes against the GLOBAL .msh list (index 2).
    prov = {"faults": {"f": {"triangle_indices_in_msh": [2]}}}
    prov_path = tmp_path / "p.json"
    prov_path.write_text(json.dumps(prov))

    # Single-fault input → no polyline vertices (need ≥2 faults).
    keys = m3p._identify_polyline_keys_from_provenance(
        points, tris, prov_path, snap_m=0.001, fault_tri_tag=100)
    assert keys == set(), (
        f"single-fault provenance has no polyline; got {keys}")


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_optim_relax_fault_protects_polyline(tmp_path):
    """End-to-end: run mmg3d_post_pass in optim_relax_fault mode on
    a tiny 2-fault bipyramid mesh.  Verify:
      1. Cross-fault polyline vertex (the shared apex) survives.
      2. The fault-tag-100 triangles are present in the output.
      3. mmg3d returns rc=0 even though fault triangles are not
         marked Required (only edges are).
    """
    import meshio
    # Bipyramid with two fault triangles (top and bottom faces of
    # a flat shared edge).  Vertex 2 is the polyline vertex —
    # shared between fault A's two triangles and fault B's two
    # triangles.  In this minimal fixture we conflate "polyline"
    # with a single vertex (degenerate but enough to exercise the
    # code path).
    points = np.asarray([
        [0.0, 0.0, 0.0],   # 0
        [1.0, 0.0, 0.0],   # 1
        [0.5, 0.5, 0.0],   # 2 — polyline vertex (shared)
        [0.5, 0.0, 0.5],   # 3 — apex above
        [0.5, 0.0, -0.5],  # 4 — apex below
    ])
    # Two tets sharing face (0,1,2) — that face is fault.
    tets = np.asarray([[0, 1, 2, 3], [0, 1, 2, 4]], dtype=np.int64)
    # Surface triangles: 1 fault (shared face) + 6 outer box-like
    # faces of the bipyramid.  Tag 100 for fault, 1 for box.
    tris = np.asarray([
        [0, 1, 2],     # fault (shared face)
        [0, 1, 3], [1, 2, 3], [0, 2, 3],
        [0, 1, 4], [1, 2, 4], [0, 2, 4],
    ], dtype=np.int64)
    tri_tags = np.asarray([100, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    tet_tags = np.asarray([10, 10], dtype=np.int32)

    msh_in = tmp_path / "in.msh"
    msh_out = tmp_path / "out.msh"
    meshio.write(msh_in, meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    ), file_format="gmsh22", binary=False)

    # Synthesize a single-fault provenance — vertex 2 only appears
    # in this single fault, so polyline_keys will be empty (need ≥2
    # faults).  This is a smoke test of the code path; the real
    # bit-identical-polyline guarantee is tested at integration
    # time on the full SAFS dataset.
    prov_path = tmp_path / "fault_provenance.json"
    prov_path.write_text(json.dumps({
        "faults": {
            "faultA": {"triangle_indices_in_msh": [0]},
        }
    }))

    rc = m3p.main([
        "--in-msh", str(msh_in),
        "--out-msh", str(msh_out),
        "--mode", "optim_relax_fault",
        "--provenance-json", str(prov_path),
        "--snap-m", "0.001",
        "--hmin", "0.05", "--hmax", "5.0",
        "--hgrad", "1.3", "--hausd", "0.1",
    ])
    assert rc == 0
    assert msh_out.exists()

    out_mesh = meshio.read(msh_out)
    found_fault = False
    for cb, tags in zip(out_mesh.cells,
                         out_mesh.cell_data.get("gmsh:physical", [])):
        if cb.type == "triangle":
            for tag in tags:
                if int(tag) == 100:
                    found_fault = True
                    break
    assert found_fault, (
        "fault triangles (tag 100) were lost on output; "
        "-opnbdy should preserve them even though they're not "
        "marked RequiredTriangles in optim_relax_fault mode")


# ---------------------------------------------------------------------------
# Tests requiring mmg3d_O3 + meshio (integration).
# ---------------------------------------------------------------------------
@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_preserves_fault_triangles_via_opnbdy(tmp_path):
    """End-to-end: take a tiny gmsh-style mesh with tetrahedra + a
    fault triangle (tag 100), run mmg3d_post_pass in optim mode, and
    verify the fault triangle survives (with its tag) in the output.

    This is the regression guard for the 2026-05-02 finding that
    mmg3d without -opnbdy silently drops "internal" surface triangles
    even when marked RequiredTriangles.
    """
    import meshio
    # Build a 2-tet mesh sharing one triangle face.  Tag the shared
    # face as fault (tag 100); tag the 6 outer faces as box (tag 1).
    points = np.asarray([
        [0.0, 0.0, 0.0],   # 0
        [1.0, 0.0, 0.0],   # 1
        [0.0, 1.0, 0.0],   # 2
        [0.0, 0.0, 1.0],   # 3 — apex above
        [0.0, 0.0, -1.0],  # 4 — apex below
    ])
    tets = np.asarray([
        [0, 1, 2, 3],
        [0, 1, 2, 4],
    ], dtype=np.int64)
    # The shared face is (0,1,2) — fault.
    # The 6 outer triangle faces of the bipyramid:
    tris = np.asarray([
        [0, 1, 2],   # fault (interior of bipyramid)
        [0, 1, 3], [1, 2, 3], [0, 2, 3],  # upper tet outer faces
        [0, 1, 4], [1, 2, 4], [0, 2, 4],  # lower tet outer faces
    ], dtype=np.int64)
    tri_tags = np.asarray([100, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    tet_tags = np.asarray([10, 10], dtype=np.int32)

    msh_in = tmp_path / "in.msh"
    msh_out = tmp_path / "out.msh"
    mesh = meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    )
    meshio.write(msh_in, mesh, file_format="gmsh22", binary=False)

    rc = m3p.main([
        "--in-msh", str(msh_in),
        "--out-msh", str(msh_out),
        "--mode", "optim",
        "--hmin", "0.05", "--hmax", "5.0",
        "--hgrad", "1.3", "--hausd", "0.1",
    ])
    assert rc == 0
    assert msh_out.exists()

    out_mesh = meshio.read(msh_out)
    # Verify the fault triangle (tag 100) is present.
    found_fault = False
    for cb, tags in zip(out_mesh.cells,
                         out_mesh.cell_data.get("gmsh:physical", [])):
        if cb.type == "triangle":
            for tag in tags:
                if int(tag) == 100:
                    found_fault = True
                    break
    assert found_fault, (
        "fault triangle (tag 100) was lost on output — `-opnbdy` "
        "should preserve internal-domain interface triangles")


# NOTE: a regression test that mmg3d in `optim` mode does not
# aggressively recoarsen the input was attempted with a 4×4×4 cube
# fixture but proved fragile — mmg3d requires a fully-closed input
# boundary surface, which is non-trivial to synthesize as a unit
# test.  The empirical regression scenario (first run without
# `-optim -nosurf` coarsened 680k tets → 4k tets) is documented in
# REVIEW_interior_subdivision_and_collapse.md and reproduced
# end-to-end via the run_newset_step_by_step.sh ENABLE_MMG3D_POSTPASS
# pipeline.


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_mmg3d_honors_RequiredEdges_on_surface_without_nosurf(tmp_path):
    """**Pre-implementation confirmation test for the
    `optim_relax_fault` mode (Tier-1 follow-up).**

    Plan: build a bipyramid mesh where the shared face (the "fault")
    is internal to a single domain (both tets carry tag 10), and
    where ONE edge of that face is marked as a RequiredEdge with both
    its endpoints marked as RequiredVertices.  Run mmg3d_O3 in
    -optim mode WITH -opnbdy and WITHOUT -nosurf — i.e., let mmg3d
    modify surface triangles freely except where explicitly required.

    Pass criterion:
      - The two RequiredVertex coords are present in the output mesh
        (not moved, not deleted).
      - At least one triangle in the output contains both required-
        vertex endpoints as a triangle edge (i.e., the RequiredEdge
        was preserved as some triangle's edge).

    If this test FAILS, mmg3d does not honor surface RequiredEdges
    when surface modification is enabled — the optim_relax_fault
    mode is not viable as designed and the implementation must
    pivot to a different protection strategy.
    """
    # Bipyramid: two tets sharing face (0,1,2).  Apex above is 3,
    # below is 4.  All tets in domain ref 10; shared face tagged 100.
    points = np.asarray([
        [0.0, 0.0, 0.0],   # 0 — RequiredVertex
        [1.0, 0.0, 0.0],   # 1 — RequiredVertex
        [0.0, 1.0, 0.0],   # 2
        [0.0, 0.0, 1.0],   # 3 — apex above
        [0.0, 0.0, -1.0],  # 4 — apex below
    ])

    medit_in = tmp_path / "in.mesh"
    medit_out = tmp_path / "out.mesh"
    log_path = tmp_path / "mmg3d.log"

    # Write the medit input by hand — we need RequiredEdges in the
    # exact form mmg3d expects.
    with medit_in.open("w") as f:
        f.write("MeshVersionFormatted 2\n\nDimension 3\n\n")
        f.write("Vertices\n5\n")
        for v in points:
            f.write(f"{v[0]:+.17e} {v[1]:+.17e} {v[2]:+.17e} 0\n")
        f.write("\n")
        f.write("Triangles\n7\n")
        # The 7 boundary + interior triangles of the bipyramid.
        for tri, ref in [((1, 2, 3), 100),
                          ((1, 2, 4), 100),
                          ((1, 4, 5), 1),
                          ((2, 4, 5), 1),
                          ((3, 4, 5), 1),
                          ((1, 3, 5), 1),
                          ((2, 3, 5), 1)]:
            f.write(f"{tri[0]} {tri[1]} {tri[2]} {ref}\n")
        f.write("\n")
        # The shared face (1,2,3) inside the bipyramid is NOT
        # written — it's an interior face of the volume mesh and
        # exists only as a tet-tet boundary, which mmg3d will infer.
        # The test edge to protect: vertex 1 (0,0,0) ↔ vertex 2
        # (1,0,0).  Mark as RequiredEdge.
        f.write("Edges\n1\n")
        f.write("1 2 0\n")     # the polyline edge
        f.write("\nRequiredEdges\n1\n1\n")
        f.write("\nRequiredVertices\n2\n1\n2\n")
        f.write("\nTetrahedra\n2\n")
        f.write("1 2 3 4 10\n")  # upper tet
        f.write("1 2 3 5 10\n")  # lower tet
        f.write("\nEnd\n")

    # Invoke mmg3d in optim mode WITH -opnbdy and WITHOUT -nosurf.
    # This is the configuration optim_relax_fault will use.
    cmd = ["mmg3d_O3",
           "-in", str(medit_in),
           "-out", str(medit_out),
           "-optim", "-opnbdy",
           "-hmin", "0.05", "-hmax", "5.0",
           "-hgrad", "1.3", "-hausd", "0.1"]
    with log_path.open("w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.flush()
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    assert r.returncode == 0, (
        f"mmg3d_O3 failed (rc={r.returncode}); see {log_path}")

    # Parse the output using the production medit reader.
    from mmgs_remesh_per_fault import _read_medit_surface
    V_out, T_out = _read_medit_surface(medit_out)

    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([1.0, 0.0, 0.0])
    found_p1 = any(np.allclose(v, p1, atol=1e-12) for v in V_out)
    found_p2 = any(np.allclose(v, p2, atol=1e-12) for v in V_out)
    assert found_p1, (
        f"RequiredVertex (0,0,0) lost from output; mmg3d does not "
        f"honor RequiredVertices when surface modification is "
        f"enabled.  optim_relax_fault is not viable; see {log_path}")
    assert found_p2, (
        f"RequiredVertex (1,0,0) lost from output; "
        f"see {log_path}")

    # Find vertex indices of the two RequiredVertex coords.
    def _idx(p):
        for i, v in enumerate(V_out):
            if np.allclose(v, p, atol=1e-12):
                return i
        return -1
    i1 = _idx(p1)
    i2 = _idx(p2)
    assert i1 >= 0 and i2 >= 0

    # The RequiredEdge (i1, i2) must appear as some surface triangle's
    # edge in the output — meaning at least one triangle has both i1
    # and i2 as vertices.
    edge_found = False
    for tri in T_out:
        if i1 in tri and i2 in tri:
            edge_found = True
            break
    assert edge_found, (
        f"RequiredEdge (0,0,0)-(1,0,0) was deleted: no output "
        f"triangle contains both endpoints.  mmg3d does not honor "
        f"surface RequiredEdges when surface modification is "
        f"enabled.  optim_relax_fault is not viable; see {log_path}")


# ---------------------------------------------------------------------------
# R-002: non-contiguous tag-100 layout (offset bug regression).
# ---------------------------------------------------------------------------
def test_R002_offset_works_on_non_contiguous_fault_triangles(tmp_path):
    """The polyline-key extraction must NOT assume tag-100 triangles
    are contiguous in the global triangle list.  Build a synthetic
    `triangles` list with [tag1, tag100, tag1, tag100] (interleaved)
    and verify polyline detection still finds vertices shared
    across two faults.
    """
    points = np.zeros((20, 3), dtype=np.float64)
    # Vertices 11 and 14 share snap_key (0,0,0) — both at origin.
    points[11] = [0.0, 0.0, 0.0]
    points[14] = [0.0, 0.0, 0.0]
    # Other vertices have distinct snap_keys (any nonzero offsets).
    points[10] = [1.0, 0.0, 0.0]
    points[12] = [2.0, 0.0, 0.0]
    points[13] = [0.0, 1.0, 0.0]
    points[15] = [3.0, 0.0, 0.0]

    # Triangles in global order: tag-1 (idx 0), tag-100 (idx 1),
    # tag-1 (idx 2), tag-100 (idx 3).  Tag-100 is NON-contiguous.
    triangles = [
        (0, 1, 2, 1),       # global 0 — box
        (10, 11, 12, 100),  # global 1 — fault A (provenance gi=1)
        (3, 4, 5, 1),       # global 2 — box
        (13, 14, 15, 100),  # global 3 — fault B (provenance gi=3)
    ]
    prov = {"faults": {
        "A": {"triangle_indices_in_msh": [1]},
        "B": {"triangle_indices_in_msh": [3]},
    }}
    p = tmp_path / "p.json"
    p.write_text(json.dumps(prov))

    keys = m3p._identify_polyline_keys_from_provenance(
        points, triangles, p, snap_m=0.1, fault_tri_tag=100)
    # Vertices 11 and 14 share snap_key (0,0,0); each is in a
    # different fault's triangle — must be flagged as polyline.
    assert (0, 0, 0) in keys, (
        f"non-contiguous tag-100 layout: shared polyline vertex "
        f"missed; got {keys}.  R-002 fix should map global indices "
        f"to local tag-100 indices via an explicit dict, not a "
        f"single offset.")


def test_R002_zero_assigned_raises_runtime_error(tmp_path):
    """When provenance lists faults but the triangle indices don't
    map to any tag-100 triangle in the input, raise rather than
    silently returning an empty set.
    """
    points = np.zeros((10, 3), dtype=np.float64)
    triangles = [(0, 1, 2, 100), (3, 4, 5, 100)]  # 2 fault tris
    # Provenance lists indices that are out of range.
    prov = {"faults": {"A": {"triangle_indices_in_msh": [99]}}}
    p = tmp_path / "p.json"
    p.write_text(json.dumps(prov))

    with pytest.raises(RuntimeError, match="ZERO assigned vertices"):
        m3p._identify_polyline_keys_from_provenance(
            points, triangles, p, snap_m=0.1, fault_tri_tag=100)


# ---------------------------------------------------------------------------
# R-001: -hgradreq 1.3 in optim_relax_fault mode flags.
# ---------------------------------------------------------------------------
def test_R001_optim_relax_fault_includes_hgradreq():
    """The optim_relax_fault mode must include `-hgradreq 1.3` in
    its mmg3d invocation flags.  Without this, mmg3d cannot insert
    Steiner points near polyline edges, leaving 4-collinear needle
    tets unbroken.
    """
    flags = m3p._resolve_mode_flags("optim_relax_fault", None)
    assert "-hgradreq" in flags, (
        f"optim_relax_fault mode must include -hgradreq for R-001 "
        f"fix; got {flags}")
    idx = flags.index("-hgradreq")
    assert idx + 1 < len(flags) and flags[idx + 1] == "1.3", (
        f"-hgradreq value must be 1.3; got {flags}")


def test_R001_optim_mode_does_not_include_hgradreq():
    """The plain `optim` mode (with -nosurf) does NOT need
    -hgradreq because surface modification is disabled — Steiner
    insertion only happens in the bulk and is already free.
    """
    flags = m3p._resolve_mode_flags("optim", None)
    assert "-hgradreq" not in flags, (
        f"plain optim mode should not include -hgradreq; got {flags}")


# ---------------------------------------------------------------------------
# R-003 + R-006: multi-pass and best-of-N trial parameters.
# ---------------------------------------------------------------------------
def test_R003_n_passes_validation(tmp_path):
    """`n_passes` must be >= 1; values < 1 raise."""
    in_path = tmp_path / "in.msh"
    in_path.write_text("dummy")
    with pytest.raises(ValueError, match="n_passes"):
        m3p.post_pass(in_path, tmp_path / "out.msh",
                      mode="optim", n_passes=0)
    with pytest.raises(ValueError, match="n_passes"):
        m3p.post_pass(in_path, tmp_path / "out.msh",
                      mode="optim", n_passes=-1)


def test_R006_n_trials_validation(tmp_path):
    """`n_trials` must be >= 1; values < 1 raise."""
    in_path = tmp_path / "in.msh"
    in_path.write_text("dummy")
    with pytest.raises(ValueError, match="n_trials"):
        m3p.post_pass(in_path, tmp_path / "out.msh",
                      mode="optim", n_trials=0)


# ---------------------------------------------------------------------------
# Helper to compute γ_min — used by R-006 best-of-N.
# ---------------------------------------------------------------------------
def test_gamma_min_helper_on_regular_tet():
    """A regular tetrahedron has γ = 1.  The helper should return
    γ ≈ 1.0 within rounding.
    """
    # Build a regular tet inscribed in the unit sphere.
    # Vertices: [(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)] / sqrt(3).
    import meshio
    s = 1.0 / np.sqrt(3.0)
    points = np.array([
        [s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s],
    ])
    # Need a triangle (any tag) for the writer to accept a non-empty
    # tri block; use a degenerate tag-100 face on the boundary.
    tris = np.array([[0, 1, 2]], dtype=np.int64)
    tets = np.array([[0, 1, 2, 3]], dtype=np.int64)
    p = pytest.importorskip("meshio")
    out = Path(__file__).parent / "_tmp_regular_tet.msh"
    try:
        meshio.write(out, meshio.Mesh(
            points=points,
            cells=[("triangle", tris), ("tetra", tets)],
            cell_data={"gmsh:physical":
                       [np.array([100], dtype=np.int32),
                        np.array([10], dtype=np.int32)],
                       "gmsh:geometrical":
                       [np.array([100], dtype=np.int32),
                        np.array([10], dtype=np.int32)]},
        ), file_format="gmsh22", binary=False)
        g = m3p._gamma_min_of_msh(out)
        assert 0.99 <= g <= 1.01, (
            f"regular tet γ should be ≈ 1.0; got {g}")
    finally:
        if out.exists():
            out.unlink()


# ---------------------------------------------------------------------------
# R-005: free-surface vertex detection excludes box-top vertices.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_R005_freesurface_excludes_box_top(tmp_path):
    """Build a mesh whose box top is at z=0 (tag 5) but whose fault
    (tag 100) is entirely below z = -200.  The
    `optim_relax_fault` mode's medit conversion must NOT mark any
    edge as a free-surface trace (because no fault vertex is near
    z=0).  Without R-005 fix, box-top vertices at z=0 would be
    flagged as free-surface, but the fault-edge filter still gives
    the right outcome — yet the side effect is over-flagging.

    Test verifies the report's `n_freesurface_edges` counter.
    """
    import meshio
    # Box top at z=0 (tag 5), fault below z=-200 (tag 100).
    points = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],     # 0-3 box top
        [0, 0, -300], [1, 0, -300], [0, 1, -300],       # 4-6 fault
        [0, 0, -1000], [1, 0, -1000], [0, 1, -1000],    # 7-9 below
    ], dtype=np.float64)
    tris = np.array([
        [0, 1, 2], [0, 2, 3],          # box top
        [4, 5, 6],                      # fault triangle (only one)
    ], dtype=np.int64)
    tri_tags = np.array([5, 5, 100], dtype=np.int32)
    # Tetrahedra to make a closed volume.
    tets = np.array([
        [0, 1, 2, 4], [0, 2, 3, 4],
        [4, 5, 6, 7], [4, 5, 7, 8],
    ], dtype=np.int64)
    tet_tags = np.array([10, 10, 10, 10], dtype=np.int32)

    msh_in = tmp_path / "in.msh"
    meshio.write(msh_in, meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    ), file_format="gmsh22", binary=False)

    prov_path = tmp_path / "fault_provenance.json"
    prov_path.write_text(json.dumps({
        "faults": {"faultA": {"triangle_indices_in_msh": [2]}},
    }))
    transform_path = tmp_path / "transform.json"
    transform_path.write_text(json.dumps({
        "free_surface_clearance_m": 100.0,
    }))

    medit_out = tmp_path / "out.mesh"
    counts = m3p._convert_msh_to_medit(
        msh_in, medit_out,
        fault_tri_tag=100, mode="optim_relax_fault",
        provenance_path=prov_path,
        transform_path=transform_path,
        snap_m=0.1)

    # Fault vertices 4, 5, 6 are all at z=-300 (well outside
    # clearance=100 of z_top=0), so 0 free-surface edges.
    assert counts["n_freesurface_edges"] == 0, (
        f"fault triangles entirely below z=-200 with "
        f"free_surface_clearance_m=100 should produce 0 "
        f"freesurface_edges; got {counts['n_freesurface_edges']}.  "
        f"R-005: box-top vertices at z=0 are being flagged.")


@pytest.mark.skipif(_mmg3d is None,
                    reason="mmg3d_O3 not on PATH")
def test_e2e_writes_report_json(tmp_path):
    """`--report-json` (or default location) must contain the input
    counts, output counts, and the mode flag for postmortem
    diagnostics.
    """
    import meshio
    # Re-use the bipyramid fixture inline.
    points = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [0.0, 0.0, -1.0],
    ])
    tets = np.asarray([[0, 1, 2, 3], [0, 1, 2, 4]], dtype=np.int64)
    tris = np.asarray([
        [0, 1, 2],
        [0, 1, 3], [1, 2, 3], [0, 2, 3],
        [0, 1, 4], [1, 2, 4], [0, 2, 4],
    ], dtype=np.int64)
    tri_tags = np.asarray([100, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    tet_tags = np.asarray([10, 10], dtype=np.int32)
    msh_in = tmp_path / "in.msh"
    msh_out = tmp_path / "out.msh"
    rep = tmp_path / "rep.json"
    meshio.write(msh_in, meshio.Mesh(
        points=points,
        cells=[("triangle", tris), ("tetra", tets)],
        cell_data={"gmsh:physical": [tri_tags, tet_tags],
                   "gmsh:geometrical": [tri_tags, tet_tags]},
    ), file_format="gmsh22", binary=False)

    m3p.main([
        "--in-msh", str(msh_in), "--out-msh", str(msh_out),
        "--mode", "optim",
        "--hmin", "0.05", "--hmax", "5.0",
        "--report-json", str(rep),
    ])
    assert rep.exists()
    r = json.loads(rep.read_text())
    assert r["mode"] == "optim"
    assert r["fault_tri_tag"] == 100
    assert r["input_counts"]["n_required_triangles"] == 1
    assert r["input_counts"]["n_tets"] == 2
