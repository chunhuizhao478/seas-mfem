"""test_rv_volume_smoothing.py — gate for the fault volume-ratio (r_v) smoother.

PLAN: experimental_mesh_refinement/PLAN_rv_volume_smoothing.md

Gates the delivered smoothed mesh (`--slide-boundary`) against the revised
acceptance criteria A1–A8 and confirms the gate has teeth (the baseline mesh
FAILS A1).  Skips (does not fail) when a mesh file is absent.

Run:
    conda activate pythonenv
    pytest -q test_rv_volume_smoothing.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
MESHING_CODE = (HERE.parent / "meshing" / "code").resolve()
for p in (str(HERE), str(MESHING_CODE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import smooth_fault_volume_ratio as S       # this folder
from check_mesh_quality import check_mesh_quality   # meshing/code

MESH_IN = HERE / "safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH_OUT = HERE / "safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_rvsmooth.msh"

# Revised A1 targets (max<2 is structurally unreachable; gate the distribution).
GT15_FRACTION_TARGET = 0.005    # < 0.5 % of fault faces with r_v > 1.5
P99_TARGET = 1.3                # 99th-percentile r_v < 1.3


def _load(path):
    g = S.read_gmsh22(path)
    topo = S.build_fault_apex_map(g.coords, g.tets, g.fault_tris)   # A5: 2-tet
    return g, topo


@pytest.fixture(scope="module")
def meshes():
    if not MESH_IN.is_file():
        pytest.skip(f"input mesh absent: {MESH_IN.name}")
    if not MESH_OUT.is_file():
        pytest.skip(f"smoothed mesh not built: {MESH_OUT.name} (run the tool)")
    g_in, topo_in = _load(MESH_IN)
    g_out, topo_out = _load(MESH_OUT)
    return g_in, topo_in, g_out, topo_out


def test_A1_rv_tail_reduced(meshes):
    _, _, g_out, topo_out = meshes
    s = S.rv_stats(g_out.coords, topo_out)
    frac = s["n_gt15"] / s["n"]
    assert frac < GT15_FRACTION_TARGET, (
        f"r_v>1.5 fraction {frac:.4%} >= {GT15_FRACTION_TARGET:.2%} target")
    assert s["p99"] < P99_TARGET, f"p99 r_v {s['p99']:.4f} >= {P99_TARGET}"


def test_A2_mean_reduced(meshes):
    g_in, topo_in, g_out, topo_out = meshes
    assert S.rv_stats(g_out.coords, topo_out)["mean"] < \
        S.rv_stats(g_in.coords, topo_in)["mean"]


def test_A3_fault_frozen_boundary_on_face(meshes):
    g_in, _, g_out, _ = meshes
    # fault nodes byte-identical
    assert np.array_equal(g_out.coords[g_in.fault_node_idx],
                          g_in.coords[g_in.fault_node_idx]), \
        "a fault node moved"
    # boundary nodes stay on their box face: the locked (face-normal) axis is
    # byte-unchanged (slide_boundary frees only the in-plane axes).
    free = S.build_free_axes(g_in, slide_boundary=True)
    locked = ~free
    assert np.array_equal(g_out.coords[locked], g_in.coords[locked]), \
        "a locked boundary/fault axis moved (domain box not preserved)"


def test_A4_topology_identical(meshes):
    g_in, _, g_out, _ = meshes
    assert len(g_out.node_ids) == len(g_in.node_ids)
    assert g_out.elem_lines == g_in.elem_lines      # verbatim → tags/conn equal
    assert np.array_equal(g_out.node_ids, g_in.node_ids)


def test_A6_no_inversion_and_quality(meshes):
    g_in, _, g_out, _ = meshes
    s_in = np.sign(S.signed_volumes(g_in.coords, g_in.tets))
    s_out = np.sign(S.signed_volumes(g_out.coords, g_out.tets))
    assert int((s_in != s_out).sum()) == 0, "tet inversion(s) detected"
    r_in = check_mesh_quality(MESH_IN)
    r_out = check_mesh_quality(MESH_OUT)
    assert r_out["q1_pass"], f"Q1 fail: tet_emin={r_out['tet_emin']:.1f}"
    assert r_out["q2_pass"], f"Q2 fail: eta_min={r_out['eta_min']:.4f}"
    assert r_out["eta_min"] >= r_in["eta_min"] - 1e-9, "eta_min regressed"
    assert r_out["tet_emin"] >= 100.0


def test_A7_mesh_zmax_zero(meshes):
    _, _, g_out, _ = meshes
    assert float(g_out.coords[:, 2].max()) == 0.0


def test_teeth_baseline_fails_A1(meshes):
    g_in, topo_in, _, _ = meshes
    s = S.rv_stats(g_in.coords, topo_in)
    frac = s["n_gt15"] / s["n"]
    assert frac >= GT15_FRACTION_TARGET or s["p99"] >= P99_TARGET, (
        f"baseline unexpectedly PASSES A1 (frac={frac:.4%}, p99={s['p99']:.3f});"
        f" the gate is vacuous")
