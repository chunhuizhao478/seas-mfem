"""check_mesh_quality.py — bulk + fault mesh-quality verifier.

PLAN_mesh_quality.md Phase 4 — verification harness.

Computes the same Joe-Liu eta and triangle-q metrics used in the project's
quality charts and enforces two hard gates:

  Q1: min bulk tet edge >= min_edge_floor_m   (default 100 m)
  Q2: min Joe-Liu eta   >  min_eta_floor      (default 0.1)

Metric definitions (see PLAN_mesh_quality.md):
  Joe-Liu eta   = 12 * (3V)^(2/3) / sum(edge_length^2)   in [0, 1]
  Triangle q    = 4*sqrt(3)*A / sum(edge_length^2)       in [0, 1]
  Histogram bins: [0.0, 0.1, 0.3, 0.5, 0.7, 1.0+1e-9]

Usage:
    # programmatic
    from check_mesh_quality import check_mesh_quality
    report = check_mesh_quality(Path('mesh.msh'))

    # CLI: exits 0 iff every mesh passes both Q1 and Q2.
    python check_mesh_quality.py PATH [PATH ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import meshio
import numpy as np


QBINS = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0 + 1.0e-9]


def _extract_blocks(m: meshio.Mesh):
    """Return (tetra_array, fault_tri_array, points) for a mesh.

    Raises ValueError if the mesh has no tetra block.
    fault_tri_array may be None if Physical Surface 101 is absent.
    """
    pts = m.points
    if "gmsh:physical" not in m.cell_data:
        raise ValueError(
            "mesh has no 'gmsh:physical' cell data; cannot identify "
            "the fault Physical Surface 101.  Was it produced by "
            "run_nwcut_meshing.py?")
    phys = m.cell_data["gmsh:physical"]
    tetra = None
    fault_tri = None
    fault_blocks = []
    for cb, tags in zip(m.cells, phys):
        if cb.type == "tetra":
            tetra = cb.data
        elif cb.type == "triangle":
            tag_arr = np.asarray(tags)
            mask = tag_arr == 101
            if mask.any():
                fault_blocks.append(cb.data[mask])
    if tetra is None or tetra.shape[0] == 0:
        raise ValueError("mesh has no tetra block (or it is empty)")
    if fault_blocks:
        # If multiple blocks carry tag 101, take their union (matches
        # the convention used by project_to_fault_stress.py's reader).
        fault_tri = (np.concatenate(fault_blocks, axis=0)
                     if len(fault_blocks) > 1 else fault_blocks[0])
    return tetra, fault_tri, pts


def _tet_metrics(tetra: np.ndarray, pts: np.ndarray):
    """Return (edge_lengths (N,6), joe_liu_eta (N,)) for the tetra block."""
    epairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    A = pts[tetra[:, [p[0] for p in epairs]]]
    B = pts[tetra[:, [p[1] for p in epairs]]]
    tet_e = np.linalg.norm(A - B, axis=2)
    p0 = pts[tetra[:, 0]]
    p1 = pts[tetra[:, 1]]
    p2 = pts[tetra[:, 2]]
    p3 = pts[tetra[:, 3]]
    V = np.abs(np.einsum('ij,ij->i', p1 - p0,
                         np.cross(p2 - p0, p3 - p0))) / 6.0
    # Joe-Liu eta = 12 * (3V)^(2/3) / sum(L^2); regular tet -> 1.
    sumL2 = (tet_e ** 2).sum(axis=1)
    # Guard against zero-volume / zero-edge degeneracies so we report
    # eta = 0 rather than crashing (sliver detector wants to see them).
    with np.errstate(divide='ignore', invalid='ignore'):
        eta = np.where(sumL2 > 0.0,
                       12.0 * np.cbrt((3.0 * V) ** 2) / sumL2,
                       0.0)
    return tet_e, eta


def _tri_metrics(fault_tri: np.ndarray, pts: np.ndarray):
    """Return (edge_lengths (N,3), areas (N,), q (N,))."""
    tpairs = [(0, 1), (1, 2), (2, 0)]
    Ta = pts[fault_tri[:, [p[0] for p in tpairs]]]
    Tb = pts[fault_tri[:, [p[1] for p in tpairs]]]
    tri_e = np.linalg.norm(Ta - Tb, axis=2)
    v0 = pts[fault_tri[:, 0]]
    v1 = pts[fault_tri[:, 1]]
    v2 = pts[fault_tri[:, 2]]
    tri_a = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    sumL2 = (tri_e ** 2).sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        tri_q = np.where(sumL2 > 0.0,
                         (4.0 * np.sqrt(3.0) * tri_a) / sumL2,
                         0.0)
    return tri_e, tri_a, tri_q


def check_mesh_quality(
    msh_path: Path,
    min_edge_floor_m: float = 100.0,
    min_eta_floor: float = 0.1,
) -> dict:
    """Return a structured quality report for a gmsh mesh.

    Parameters
    ----------
    msh_path : Path
        Path to a .msh file produced by run_nwcut_meshing.py.
    min_edge_floor_m : float
        Q1 floor: minimum bulk tet edge length is required to be >= this.
    min_eta_floor : float
        Q2 floor: minimum Joe-Liu tet quality is required to be strictly
        greater than this.

    Returns
    -------
    dict with keys:
        n_tet, n_tri_fault, mesh_zmax,
        tet_emin, tet_emed, tet_emax,
        eta_min, eta_med, eta_mean, eta_hist (5-bucket counts),
        tri_emin, tri_emed, tri_emax,
        tri_qmin, tri_qmed, tri_qmean, tri_q_hist,
        q1_pass: bool   # tet_emin >= min_edge_floor_m
        q2_pass: bool   # eta_min  >  min_eta_floor
        violations: list[str]

    Triangle fields fall back to None if Physical Surface 101 is absent.

    Raises
    ------
    FileNotFoundError if msh_path does not exist.
    ValueError if the mesh has no tetra block.
    """
    p = Path(msh_path)
    if not p.is_file():
        raise FileNotFoundError(f"mesh not found: {p}")
    m = meshio.read(str(p))
    tetra, fault_tri, pts = _extract_blocks(m)

    tet_e, eta = _tet_metrics(tetra, pts)
    eta_hist, _ = np.histogram(eta, bins=QBINS)

    tet_emin = float(tet_e.min())
    tet_emax = float(tet_e.max())
    tet_emed = float(np.median(tet_e))
    eta_min = float(eta.min())
    eta_med = float(np.median(eta))
    eta_mean = float(eta.mean())

    report = {
        "n_tet": int(tetra.shape[0]),
        "n_tri_fault": int(fault_tri.shape[0]) if fault_tri is not None
                        else None,
        "mesh_zmax": float(pts[:, 2].max()),
        "tet_emin": tet_emin,
        "tet_emed": tet_emed,
        "tet_emax": tet_emax,
        "eta_min": eta_min,
        "eta_med": eta_med,
        "eta_mean": eta_mean,
        "eta_hist": eta_hist.tolist(),
        "tri_emin": None, "tri_emed": None, "tri_emax": None,
        "tri_qmin": None, "tri_qmed": None, "tri_qmean": None,
        "tri_q_hist": None,
    }
    if fault_tri is not None:
        tri_e, tri_a, tri_q = _tri_metrics(fault_tri, pts)
        tri_q_hist, _ = np.histogram(tri_q, bins=QBINS)
        report["tri_emin"] = float(tri_e.min())
        report["tri_emed"] = float(np.median(tri_e))
        report["tri_emax"] = float(tri_e.max())
        report["tri_qmin"] = float(tri_q.min())
        report["tri_qmed"] = float(np.median(tri_q))
        report["tri_qmean"] = float(tri_q.mean())
        report["tri_q_hist"] = tri_q_hist.tolist()

    # Q1 / Q2 gates per plan: strict semantics.
    #   Q1: tet_emin >= min_edge_floor_m   (allow equality)
    #   Q2: eta_min  >  min_eta_floor      (strict, no equality)
    q1_pass = bool(tet_emin >= float(min_edge_floor_m))
    q2_pass = bool(eta_min > float(min_eta_floor))
    violations = []
    if not q1_pass:
        violations.append(
            f"Q1 FAIL: tet_emin={tet_emin:.2f} m < "
            f"{min_edge_floor_m:.0f} m floor")
    if not q2_pass:
        n_below = int((eta <= float(min_eta_floor)).sum())
        violations.append(
            f"Q2 FAIL: eta_min={eta_min:.4f} <= {min_eta_floor:.2f} "
            f"floor ({n_below} sliver tet(s))")
    report["q1_pass"] = q1_pass
    report["q2_pass"] = q2_pass
    report["violations"] = violations
    return report


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def _format_table(rows: list[dict]) -> str:
    """Markdown table of per-mesh Q1/Q2 results."""
    head = (
        "| mesh | n_tet | tet_emin [m] | eta_min | Q1 | Q2 |\n"
        "|---|---:|---:|---:|:-:|:-:|"
    )
    lines = [head]
    for r in rows:
        lines.append(
            f"| {r['name']} | {r['report']['n_tet']:,} | "
            f"{r['report']['tet_emin']:.2f} | "
            f"{r['report']['eta_min']:.4f} | "
            f"{'✓' if r['report']['q1_pass'] else '✗'} | "
            f"{'✓' if r['report']['q2_pass'] else '✗'} |"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("meshes", nargs="+", type=Path,
                    help="one or more gmsh .msh files to verify")
    ap.add_argument("--min-edge-floor-m", type=float, default=100.0,
                    help="Q1 floor on bulk tet min edge (default: %(default)s)")
    ap.add_argument("--min-eta-floor", type=float, default=0.1,
                    help="Q2 floor on Joe-Liu eta (default: %(default)s)")
    args = ap.parse_args(argv)

    rows = []
    rc = 0
    for p in args.meshes:
        try:
            r = check_mesh_quality(p, args.min_edge_floor_m,
                                   args.min_eta_floor)
        except (FileNotFoundError, ValueError) as exc:
            print(f"ERROR: {p}: {exc}", file=sys.stderr)
            rc = 1
            continue
        rows.append({"name": p.name, "report": r})
        if r["violations"]:
            for v in r["violations"]:
                print(f"  {p.name}: {v}", file=sys.stderr)
            rc = 1

    if rows:
        print(_format_table(rows))
    return rc


if __name__ == "__main__":
    sys.exit(main())
