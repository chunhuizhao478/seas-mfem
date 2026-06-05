"""check_mesh_quality.py — bulk + per-surface mesh-quality verifier for
GOCAD-exported Abaqus .inp meshes (project_7.0_preferred).

Adapted from project_7.0_alternative/meshing/code/check_mesh_quality.py.
The metric definitions (Joe-Liu eta, triangle q, histogram bins, Q1/Q2
gates) are IDENTICAL to that gmsh-based version; only the input layer is
changed so it consumes the GOCAD Abaqus .inp (C3D4 tets + *SURFACE
element-side sets) that this directory's gocad_inp_to_vtu.py also reads.

Two differences from the alternative (gmsh) version, both forced by the
input rather than choice:

  1. Input format.  A gmsh .msh tags the fault with Physical Surface 101
     in `gmsh:physical` cell data.  An Abaqus .inp has no such tag; faults
     and boundaries are `*SURFACE, TYPE=ELEMENT` blocks pointing at
     `*ELSET ..._S{1..4}` element-side groups.  We expand those into
     triangles with the same Abaqus C3D4 face convention used by
     gocad_inp_to_vtu.py and report q PER named surface (every *SURFACE),
     plus an aggregate over all surface triangles.

  2. Scale.  These meshes can have tens of millions of tets, so the bulk
     tet metric is computed in element-axis CHUNKS to bound peak memory
     (the alternative version's single fully-vectorized call allocates
     (N,6,3) intermediates and would exhaust RAM at ~47M tets).  The
     per-chunk call uses the unchanged _tet_metrics().

Metric definitions (unchanged from the alternative version):
  Joe-Liu eta   = 12 * (3V)^(2/3) / sum(edge_length^2)   in [0, 1]
  Triangle q    = 4*sqrt(3)*A / sum(edge_length^2)        in [0, 1]
  Histogram bins: [0.0, 0.1, 0.3, 0.5, 0.7, 1.0+1e-9]

Hard gates (bulk tets), strict semantics as in the original:
  Q1: min bulk tet edge >= min_edge_floor_m   (default 100 m)
  Q2: min Joe-Liu eta   >  min_eta_floor       (default 0.1)

Usage:
    conda activate pythonenv

    # programmatic
    from check_mesh_quality import check_inp_mesh_quality
    report = check_inp_mesh_quality(Path('SAFtopo_test1km.inp'))

    # CLI: exits 0 iff the bulk passes both Q1 and Q2.
    python check_mesh_quality.py PATH.inp [PATH.inp ...]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import meshio
import numpy as np


QBINS = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0 + 1.0e-9]
QBIN_LABELS = ["<0.1", "0.1-0.3", "0.3-0.5", "0.5-0.7", ">0.7"]

# Abaqus C3D4 face -> (0-based) node indices.  Identical to the table in
# gocad_inp_to_vtu.py (the standard Abaqus C3D4 face->node convention).
FACE_NODES = {
    1: (0, 1, 2),   # S1 opposite node 4
    2: (0, 3, 1),   # S2 opposite node 3
    3: (1, 3, 2),   # S3 opposite node 2
    4: (2, 3, 0),   # S4 opposite node 1
}


# ---------------------------------------------------------------------
# Metric kernels — kept verbatim from the alternative version so the
# quality numbers are directly comparable across the two projects.
# ---------------------------------------------------------------------

def _tet_metrics(tetra: np.ndarray, pts: np.ndarray):
    """Return (edge_lengths (N,6), joe_liu_eta (N,)) for a tetra block."""
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


def _tri_metrics(tri: np.ndarray, pts: np.ndarray):
    """Return (edge_lengths (N,3), areas (N,), q (N,)) for a triangle block."""
    tpairs = [(0, 1), (1, 2), (2, 0)]
    Ta = pts[tri[:, [p[0] for p in tpairs]]]
    Tb = pts[tri[:, [p[1] for p in tpairs]]]
    tri_e = np.linalg.norm(Ta - Tb, axis=2)
    v0 = pts[tri[:, 0]]
    v1 = pts[tri[:, 1]]
    v2 = pts[tri[:, 2]]
    tri_a = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    sumL2 = (tri_e ** 2).sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        tri_q = np.where(sumL2 > 0.0,
                         (4.0 * np.sqrt(3.0) * tri_a) / sumL2,
                         0.0)
    return tri_e, tri_a, tri_q


# ---------------------------------------------------------------------
# Abaqus .inp surface extraction — same approach as gocad_inp_to_vtu.py.
# ---------------------------------------------------------------------

def _build_block_offsets(m: meshio.Mesh) -> list[int]:
    offsets = [0]
    for cb in m.cells:
        offsets.append(offsets[-1] + len(cb.data))
    return offsets


def _cell_set_global_tet_ids(m: meshio.Mesh, offsets: list[int],
                             name: str) -> np.ndarray:
    """Global tet indices (into the concatenated tet array) for a meshio
    cell set.  Non-tetra blocks in the set are dropped (GOCAD .inp is
    tet-only)."""
    if name not in m.cell_sets:
        return np.array([], dtype=np.int64)
    pieces = []
    for blk_i, ids in enumerate(m.cell_sets[name]):
        if ids is None or len(ids) == 0:
            continue
        if m.cells[blk_i].type != "tetra":
            continue
        pieces.append(np.asarray(ids, dtype=np.int64) + offsets[blk_i])
    return np.concatenate(pieces) if pieces else np.array([], dtype=np.int64)


def _extract_surface_triangles(m: meshio.Mesh, offsets: list[int],
                               all_tetra: np.ndarray,
                               surface_name: str) -> np.ndarray | None:
    """Concat a *SURFACE's S1..S4 element-side ELSETs into one (N,3)
    triangle array.  Returns None if the surface has no triangles."""
    tris = []
    for sn, nodes in FACE_NODES.items():
        gids = _cell_set_global_tet_ids(m, offsets, f"{surface_name}_S{sn}")
        if len(gids) == 0:
            continue
        tris.append(all_tetra[gids][:, list(nodes)])
    return np.concatenate(tris, axis=0) if tris else None


def _surface_names(m: meshio.Mesh) -> list[str]:
    """Surface base names, derived from the `..._S{1..4}` element-side
    ELSET keys meshio loaded into cell_sets (strip the trailing _S<n>).

    Using the cell_set keys avoids re-reading the multi-GB .inp text just
    to grep `*SURFACE` NAME= lines; every named surface owns at least one
    `_S<n>` ELSET, so the base-name set is identical.
    """
    bases = []
    seen = set()
    side = re.compile(r"^(.*)_S[1-4]$")
    for key in m.cell_sets:
        mobj = side.match(key)
        if not mobj:
            continue
        base = mobj.group(1)
        if base not in seen:
            seen.add(base)
            bases.append(base)
    return bases


# ---------------------------------------------------------------------
# Quality report
# ---------------------------------------------------------------------

def _stats(arr: np.ndarray) -> tuple[float, float, float]:
    """(min, median, max) as floats."""
    return float(arr.min()), float(np.median(arr)), float(arr.max())


def check_inp_mesh_quality(
    inp_path: Path,
    min_edge_floor_m: float = 100.0,
    min_eta_floor: float = 0.1,
    chunk: int = 2_000_000,
) -> dict:
    """Structured bulk + per-surface quality report for a GOCAD .inp mesh.

    Parameters
    ----------
    inp_path : Path
        Path to the Abaqus .inp (C3D4 tets, *SURFACE element-side sets).
    min_edge_floor_m : float
        Q1 floor: minimum bulk tet edge length must be >= this (default 100 m).
    min_eta_floor : float
        Q2 floor: minimum Joe-Liu eta must be strictly > this (default 0.1).
    chunk : int
        Number of tetra processed per metric chunk (memory bound).

    Returns
    -------
    dict with keys:
        n_node, n_tet, mesh_zmin, mesh_zmax,
        tet_emin, tet_emed, tet_emax,
        eta_min, eta_med, eta_mean, eta_hist (5-bucket counts),
        n_below_eta_floor,
        q1_pass, q2_pass, violations,
        surfaces: list of per-*SURFACE dicts
            {name, n_tri, tri_emin, tri_emed, tri_emax,
             q_min, q_med, q_mean, q_hist},
        surf_all: aggregate-over-all-surface-triangles dict (or None).

    Raises
    ------
    FileNotFoundError if inp_path does not exist.
    ValueError       if the mesh has no tetra block.
    """
    p = Path(inp_path)
    if not p.is_file():
        raise FileNotFoundError(f"mesh not found: {p}")
    if chunk <= 0:
        raise ValueError(f"chunk must be positive, got {chunk}")

    m = meshio.read(str(p))
    pts = np.asarray(m.points, dtype=np.float64)
    if pts.size == 0:
        raise ValueError(f"{p}: mesh has no points")

    tetra_blocks = [cb.data for cb in m.cells if cb.type == "tetra"]
    if not tetra_blocks or sum(b.shape[0] for b in tetra_blocks) == 0:
        raise ValueError(f"{p}: mesh has no tetra block (or it is empty)")
    all_tetra = (np.concatenate(tetra_blocks, axis=0)
                 if len(tetra_blocks) > 1 else tetra_blocks[0])
    n_tet = int(all_tetra.shape[0])

    # --- bulk tet metrics, chunked to bound peak memory ---------------
    eta_all = np.empty(n_tet, dtype=np.float64)
    edges_all = np.empty((n_tet, 6), dtype=np.float32)
    for s in range(0, n_tet, chunk):
        e = min(s + chunk, n_tet)
        tet_e, eta = _tet_metrics(all_tetra[s:e], pts)
        eta_all[s:e] = eta
        edges_all[s:e] = tet_e.astype(np.float32)

    tet_emin, tet_emed, tet_emax = _stats(edges_all)
    eta_min = float(eta_all.min())
    eta_med = float(np.median(eta_all))
    eta_mean = float(eta_all.mean())
    eta_hist, _ = np.histogram(eta_all, bins=QBINS)
    n_below = int((eta_all <= float(min_eta_floor)).sum())

    # --- per-*SURFACE triangle metrics --------------------------------
    offsets = _build_block_offsets(m)
    surfaces = []
    all_surf_tris = []
    for name in _surface_names(m):
        tris = _extract_surface_triangles(m, offsets, all_tetra, name)
        if tris is None or len(tris) == 0:
            continue
        tri_e, _tri_a, tri_q = _tri_metrics(tris, pts)
        tri_emin, tri_emed, tri_emax = _stats(tri_e)
        q_min, q_med, q_max = _stats(tri_q)
        q_hist, _ = np.histogram(tri_q, bins=QBINS)
        # Drop trailing _C3D4 for display, matching gocad_inp_to_vtu.py.
        disp = name[:-5] if name.endswith("_C3D4") else name
        surfaces.append({
            "name": disp,
            "n_tri": int(tris.shape[0]),
            "tri_emin": tri_emin, "tri_emed": tri_emed, "tri_emax": tri_emax,
            "q_min": q_min, "q_med": q_med, "q_mean": float(tri_q.mean()),
            "q_hist": q_hist.tolist(),
        })
        all_surf_tris.append(tris)

    surf_all = None
    if all_surf_tris:
        agg = np.concatenate(all_surf_tris, axis=0)
        tri_e, _tri_a, tri_q = _tri_metrics(agg, pts)
        tri_emin, tri_emed, tri_emax = _stats(tri_e)
        q_min, q_med, q_max = _stats(tri_q)
        q_hist, _ = np.histogram(tri_q, bins=QBINS)
        surf_all = {
            "name": "ALL surfaces",
            "n_tri": int(agg.shape[0]),
            "tri_emin": tri_emin, "tri_emed": tri_emed, "tri_emax": tri_emax,
            "q_min": q_min, "q_med": q_med, "q_mean": float(tri_q.mean()),
            "q_hist": q_hist.tolist(),
        }

    # --- Q1 / Q2 gates (strict semantics, as the original) ------------
    q1_pass = bool(tet_emin >= float(min_edge_floor_m))
    q2_pass = bool(eta_min > float(min_eta_floor))
    violations = []
    if not q1_pass:
        violations.append(
            f"Q1 FAIL: tet_emin={tet_emin:.2f} m < "
            f"{min_edge_floor_m:.0f} m floor")
    if not q2_pass:
        violations.append(
            f"Q2 FAIL: eta_min={eta_min:.4f} <= {min_eta_floor:.2f} "
            f"floor ({n_below} sliver tet(s))")

    return {
        "n_node": int(pts.shape[0]),
        "n_tet": n_tet,
        "mesh_zmin": float(pts[:, 2].min()),
        "mesh_zmax": float(pts[:, 2].max()),
        "tet_emin": tet_emin, "tet_emed": tet_emed, "tet_emax": tet_emax,
        "eta_min": eta_min, "eta_med": eta_med, "eta_mean": eta_mean,
        "eta_hist": eta_hist.tolist(),
        "n_below_eta_floor": n_below,
        "q1_pass": q1_pass, "q2_pass": q2_pass, "violations": violations,
        "surfaces": surfaces, "surf_all": surf_all,
    }


# ---------------------------------------------------------------------
# CLI / reporting
# ---------------------------------------------------------------------

def _fmt_hist(hist: list[int]) -> str:
    return "  ".join(f"{lab}:{n:,}" for lab, n in zip(QBIN_LABELS, hist))


def _print_report(name: str, r: dict, min_edge_floor_m: float,
                  min_eta_floor: float) -> None:
    print(f"\n========== {name} ==========")
    print(f"nodes      : {r['n_node']:,}")
    print(f"bulk tets  : {r['n_tet']:,}")
    print(f"z range    : [{r['mesh_zmin']:.1f}, {r['mesh_zmax']:.1f}] m")
    print("\n-- bulk tetrahedra --")
    print(f"edge [m]   : min {r['tet_emin']:.2f}   median {r['tet_emed']:.2f}"
          f"   max {r['tet_emax']:.2f}")
    print(f"Joe-Liu eta: min {r['eta_min']:.4f}   median {r['eta_med']:.4f}"
          f"   mean {r['eta_mean']:.4f}")
    print(f"eta hist   : {_fmt_hist(r['eta_hist'])}")
    q1 = "PASS" if r["q1_pass"] else "FAIL"
    q2 = "PASS" if r["q2_pass"] else "FAIL"
    print(f"Q1 (emin >= {min_edge_floor_m:.0f} m) : {q1}")
    print(f"Q2 (eta_min > {min_eta_floor:.2f})   : {q2}"
          f"   ({r['n_below_eta_floor']:,} tet(s) <= floor)")
    for v in r["violations"]:
        print(f"  !! {v}", file=sys.stderr)

    if r["surfaces"]:
        print("\n-- per *SURFACE (triangles) --")
        print("| surface | n_tri | edge min [m] | edge med [m] | "
              "q_min | q_med | q_mean |")
        print("|---|---:|---:|---:|---:|---:|---:|")
        rows = list(r["surfaces"])
        if r["surf_all"] is not None:
            rows.append(r["surf_all"])
        for s in rows:
            print(f"| {s['name']} | {s['n_tri']:,} | "
                  f"{s['tri_emin']:.2f} | {s['tri_emed']:.2f} | "
                  f"{s['q_min']:.4f} | {s['q_med']:.4f} | {s['q_mean']:.4f} |")
    else:
        print("\n-- per *SURFACE (triangles) -- : none found")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("meshes", nargs="+", type=Path,
                    help="one or more GOCAD Abaqus .inp files to verify")
    ap.add_argument("--min-edge-floor-m", type=float, default=100.0,
                    help="Q1 floor on bulk tet min edge (default: %(default)s)")
    ap.add_argument("--min-eta-floor", type=float, default=0.1,
                    help="Q2 floor on Joe-Liu eta (default: %(default)s)")
    ap.add_argument("--chunk", type=int, default=2_000_000,
                    help="tetra processed per metric chunk (default: %(default)s)")
    args = ap.parse_args(argv)

    rc = 0
    for p in args.meshes:
        try:
            r = check_inp_mesh_quality(p, args.min_edge_floor_m,
                                       args.min_eta_floor, args.chunk)
        except (FileNotFoundError, ValueError) as exc:
            print(f"ERROR: {p}: {exc}", file=sys.stderr)
            rc = 1
            continue
        _print_report(p.name, r, args.min_edge_floor_m, args.min_eta_floor)
        if r["violations"]:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
