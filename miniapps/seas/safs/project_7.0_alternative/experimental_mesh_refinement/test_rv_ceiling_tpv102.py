"""test_rv_ceiling_tpv102.py — hard r_v <= 1.5 ceiling gate for the TPV102
SCEC benchmark mesh.

Zhang et al. (2023) report spatial-spike oscillations (SSOs) when the Eq.(18)
fault volume ratio r_v > 1.5.  For the TPV102 benchmark we enforce a HARD guard:
the smoothed mesh must have NO fault face with r_v > 1.5, while still passing the
project bulk gates Q1 (tet_emin >= 100 m) and Q2 (eta_min > 0.1) and keeping the
fault geometry and domain box exactly.

The smoothed mesh is produced by:
    smooth_fault_volume_ratio.py --in tpv102_200m.msh --out tpv102_200m_rvsmooth.msh \
        --fault-phys 3 --rock-phys 1 --boundary-phys 1,5 \
        --slide-boundary --rv-ceiling 1.5 --max-sweeps 200

Skips (does not fail) when the meshes are absent.

Run:  conda activate pythonenv && pytest -q test_rv_ceiling_tpv102.py
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

import smooth_fault_volume_ratio as S
from check_mesh_quality import check_mesh_quality

RV_CEILING = 1.5
TPV102_TAGS = dict(fault_phys=3, rock_phys=1, boundary_phys=(1, 5))


def _find(name: str) -> Path | None:
    """Locate a TPV102 mesh in the worktree or a plain checkout."""
    repo = HERE.parents[4]
    for pat in (f".claude/worktrees/*/miniapps/seas/tpv102/mesh/{name}",
                f"miniapps/seas/tpv102/mesh/{name}"):
        hits = sorted(repo.glob(pat))
        if hits:
            return hits[0]
    return None


@pytest.fixture(scope="module")
def meshes():
    mi = _find("tpv102_200m.msh")
    mo = _find("tpv102_200m_rvsmooth.msh")
    if mi is None or mo is None:
        pytest.skip("TPV102 baseline / rvsmooth mesh not found")
    gi = S.read_gmsh22(mi, **TPV102_TAGS)
    go = S.read_gmsh22(mo, **TPV102_TAGS)
    return mi, mo, gi, go


def test_hard_ceiling_no_face_above_1p5(meshes):
    _, _, _, go = meshes
    topo = S.build_fault_apex_map(go.coords, go.tets, go.fault_tris)
    s = S.rv_stats(go.coords, topo)
    assert s["n_gt15"] == 0, (
        f"{s['n_gt15']} fault face(s) have r_v > {RV_CEILING} "
        f"(max={s['max']:.4f}); hard ceiling violated")
    assert s["max"] <= RV_CEILING + 1e-9, f"max r_v {s['max']:.4f} > {RV_CEILING}"


def test_bulk_quality_gates(meshes):
    mi, mo, _, _ = meshes
    ri, ro = check_mesh_quality(mi), check_mesh_quality(mo)
    assert ro["q1_pass"], f"Q1 fail: tet_emin={ro['tet_emin']:.1f}"
    assert ro["q2_pass"], f"Q2 fail: eta_min={ro['eta_min']:.4f}"
    assert ro["mesh_zmax"] == 0.0


def test_fault_and_box_preserved(meshes):
    _, _, gi, go = meshes
    assert np.array_equal(go.coords[gi.fault_node_idx],
                          gi.coords[gi.fault_node_idx]), "a fault node moved"
    free = S.build_free_axes(gi, slide_boundary=True)
    locked = ~free
    assert np.array_equal(go.coords[locked], gi.coords[locked]), \
        "a locked boundary/fault axis moved (domain box not preserved)"
    assert go.elem_lines == gi.elem_lines, "topology changed"


def test_no_inverted_tets(meshes):
    _, _, gi, go = meshes
    s_in = np.sign(S.signed_volumes(gi.coords, gi.tets))
    s_out = np.sign(S.signed_volumes(go.coords, go.tets))
    assert int((s_in != s_out).sum()) == 0, "tet inversion(s) detected"


def test_teeth_baseline_violates_ceiling(meshes):
    _, _, gi, _ = meshes
    topo = S.build_fault_apex_map(gi.coords, gi.tets, gi.fault_tris)
    s = S.rv_stats(gi.coords, topo)
    assert s["max"] > RV_CEILING, (
        f"baseline max r_v {s['max']:.4f} <= {RV_CEILING}; gate is vacuous")
