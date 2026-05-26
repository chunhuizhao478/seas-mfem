"""test_fault_triangle_quality.py — Phase 2 fault-triangle quality gate.

PLAN: meshing/docs/PLAN_fault_triangle_quality_2026-05-26.md (Phase 2).

Gates the *embedded* fault triangulation (Physical Surface 101) of the
remeshed mesh against the plan targets:

    worst min-angle >= ANGLE_TARGET  (default 25 deg)
    tri_qmin        >= TRIQ_TARGET   (default 0.5)

and confirms the gate has teeth by asserting the *baseline* mesh FAILS it.

Skips (does not fail) when a mesh file is absent, so the suite passes on a
fresh checkout / before Phase 2 has been run.

Run:
    conda activate pythonenv
    pytest -q test_fault_triangle_quality.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# --- locate the experimental tool (this dir) + the existing checker -----------
HERE = Path(__file__).resolve().parent
MESHING_CODE = (HERE.parent / "meshing" / "code").resolve()
for p in (str(HERE), str(MESHING_CODE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from check_mesh_quality import _extract_blocks  # existing (meshing/code)
from remesh_fault_stl import tri_metrics        # this folder (min_angle + tri_q)

import meshio

ANGLE_TARGET = 25.0
TRIQ_TARGET = 0.5

MESH_NEW = HERE / "safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh"
MESH_OLD = (HERE.parent / "meshing" / "results" / "msh"
            / "safs_fault_box_nwcut_500m_lcfar3000_z0embed.msh")


def _fault_quality(msh_path: Path) -> tuple[float, float, int]:
    """Return (worst_min_angle_deg, tri_qmin, n_below_target) for Physical
    Surface 101 of a mesh."""
    m = meshio.read(str(msh_path))
    _, fault_tri, pts = _extract_blocks(m)
    if fault_tri is None or fault_tri.shape[0] == 0:
        raise AssertionError(f"{msh_path.name}: no Physical Surface 101 tris")
    min_angle, tri_q = tri_metrics(np.asarray(pts, float),
                                   np.asarray(fault_tri, np.int64))
    n_below = int(((min_angle < ANGLE_TARGET) | (tri_q < TRIQ_TARGET)).sum())
    return float(min_angle.min()), float(tri_q.min()), n_below


def test_new_mesh_meets_fault_quality_target():
    """The remeshed mesh's fault triangulation clears the angle/q targets."""
    if not MESH_NEW.is_file():
        pytest.skip(f"mesh not built: {MESH_NEW.name} (run Phase 2 first)")
    worst_angle, qmin, n_below = _fault_quality(MESH_NEW)
    assert worst_angle >= ANGLE_TARGET, (
        f"worst fault min-angle {worst_angle:.2f} deg < {ANGLE_TARGET} deg "
        f"({n_below} tris below target)")
    assert qmin >= TRIQ_TARGET, (
        f"fault tri_qmin {qmin:.4f} < {TRIQ_TARGET} "
        f"({n_below} tris below target)")


def test_baseline_mesh_fails_gate_has_teeth():
    """Sanity: the un-remeshed baseline mesh FAILS the same gate (so a pass on
    the new mesh is meaningful, not vacuous)."""
    if not MESH_OLD.is_file():
        pytest.skip(f"baseline mesh absent: {MESH_OLD.name}")
    worst_angle, qmin, _ = _fault_quality(MESH_OLD)
    assert worst_angle < ANGLE_TARGET or qmin < TRIQ_TARGET, (
        f"baseline unexpectedly PASSES the gate "
        f"(worst_angle={worst_angle:.2f}, qmin={qmin:.4f}); the gate is "
        f"vacuous and needs tightening")
