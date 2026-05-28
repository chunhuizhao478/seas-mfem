#!/usr/bin/env python3
"""verify_250m_mesh.py — one-pass quality verification of the 250 m triqsubdiv
mesh.  Reads the (large) mesh ONCE and reports every gate that matters for the
dynamic-rupture driver:

  * Joe-Liu eta gates: Q1 (min tet edge; sub-km floor 50 m) and Q2 (eta>0.1).
  * Fault-triangle shape gate: worst min-angle >= 25 deg, tri_qmin >= 0.5
    (the property the old 250 m subdiv mesh lacked: 9.1 deg trace needles).
  * eta<0.1 sliver tets and HOW MANY SIT ON THE FAULT (>=3 fault-surface
    vertices) — these are the dynamic-rupture blow-up seed
    (DEBUG_100m_mesh_fault_slivers_2026-05-27.md).  This is the decisive check.
  * mesh_zmax (surface-rupture invariant) and fault embedding sanity.

Usage (conda activate pythonenv):
    python verify_250m_mesh.py [MESH.msh]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MESHING = (HERE.parent / "meshing" / "code").resolve()
for p in (str(HERE), str(MESHING)):
    if p not in sys.path:
        sys.path.insert(0, p)

import warnings
warnings.simplefilter("ignore")
import meshio
from check_mesh_quality import _extract_blocks, _tet_metrics, QBINS
from remesh_fault_stl import tri_metrics

DEFAULT_MESH = (HERE
                / "safs_fault_box_nwcut_250m_lcfar3000_z0embed_triqsubdiv.msh")
ANGLE_TARGET = 25.0
TRIQ_TARGET = 0.5
Q1_FLOOR = 50.0      # sub-km floor (DEBUG_z0cut_fine_resolution.md)
Q2_FLOOR = 0.1


def main(argv: list[str] | None = None) -> int:
    msh = Path(argv[0]) if argv else DEFAULT_MESH
    if not msh.is_file():
        print(f"ERROR: mesh not found: {msh}", file=sys.stderr)
        return 1
    print(f"reading {msh.name} ({msh.stat().st_size/1e6:.0f} MB) ...",
          flush=True)
    m = meshio.read(str(msh))
    tetra, fault_tri, pts = _extract_blocks(m)
    if fault_tri is None or fault_tri.shape[0] == 0:
        print("ERROR: no Physical Surface 101 (fault) in mesh", file=sys.stderr)
        return 1
    mesh_zmax = float(pts[:, 2].max())
    print(f"n_tet={tetra.shape[0]:,}  n_fault_tri={fault_tri.shape[0]:,}  "
          f"mesh_zmax={mesh_zmax:.6g}", flush=True)

    # ---- bulk Joe-Liu eta gates ----
    tet_e, eta = _tet_metrics(tetra, pts)
    tet_emin = float(tet_e.min())
    eta_min = float(eta.min())
    q1 = tet_emin >= Q1_FLOOR
    q2 = eta_min > Q2_FLOOR
    hist, _ = np.histogram(eta, bins=QBINS)
    print("\n--- bulk (Joe-Liu eta) ---")
    print(f"  tet_emin = {tet_emin:.2f} m   "
          f"Q1@{Q1_FLOOR:.0f}: {'PASS' if q1 else 'FAIL'}   "
          f"(Q1@100: {'PASS' if tet_emin >= 100 else 'FAIL — expected, ~56 m fault edges'})")
    print(f"  eta_min  = {eta_min:.4f}   eta_med = {np.median(eta):.4f}   "
          f"Q2(eta>0.1): {'PASS' if q2 else 'FAIL'}")
    print(f"  eta hist {[round(b,1) for b in QBINS]}: {hist.tolist()}")

    # ---- fault-triangle shape gate (verbatim-embedded from the subdiv STL) ----
    ang, q = tri_metrics(np.asarray(pts, float),
                         np.asarray(fault_tri, np.int64))
    fault_worst = float(ang.min())
    tri_qmin = float(q.min())
    gate = (fault_worst >= ANGLE_TARGET) and (tri_qmin >= TRIQ_TARGET)
    print("\n--- fault triangles (Physical Surface 101) ---")
    print(f"  worst min-angle = {fault_worst:.2f} deg   "
          f"tris<25 = {int((ang < ANGLE_TARGET).sum())}")
    print(f"  tri_qmin = {tri_qmin:.4f}   tri_q<0.5 = {int((q < TRIQ_TARGET).sum())}")
    print(f"  fault gate (>=25 deg & q>=0.5): {'PASS' if gate else 'FAIL'}")

    # ---- DECISIVE: eta<0.1 slivers, and how many are ON the fault ----
    sl = np.where(eta < Q2_FLOOR)[0]
    print("\n--- eta<0.1 sliver location (blow-up seed check) ---")
    print(f"  eta<0.1 sliver tets: {sl.size}")
    on_fault = 0
    if sl.size:
        fn = np.zeros(pts.shape[0], dtype=bool)
        fn[np.unique(fault_tri)] = True
        nfault = fn[tetra[sl]].sum(axis=1)        # fault verts per sliver tet
        on_fault = int((nfault >= 3).sum())
        z = pts[tetra[sl]].mean(axis=1)[:, 2]
        print(f"    ON fault (>=3 fault verts, a fault FACE): {on_fault}")
        print(f"    touching fault (>=1 fault vert)         : {int((nfault>=1).sum())}")
        print(f"    sliver centroid z: min={z.min():.0f} "
              f"med={np.median(z):.0f} max={z.max():.0f} m")

    # ---- verdict ----
    ok = (q1 and q2 and gate and mesh_zmax <= 1e-6 and on_fault == 0)
    print("\n=== VERDICT ===")
    print(f"  Q1@50={q1}  Q2={q2}  fault_gate={gate}  "
          f"mesh_zmax=0:{mesh_zmax <= 1e-6}  on_fault_slivers={on_fault}")
    print(f"  {'ALL GATES PASS' if ok else 'GATE FAILURE — inspect above'}")
    return 0 if ok else 4


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
