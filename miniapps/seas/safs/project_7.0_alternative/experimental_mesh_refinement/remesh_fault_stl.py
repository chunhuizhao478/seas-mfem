#!/usr/bin/env python3
"""remesh_fault_stl.py — surface-remesh a cut SAFS fault STL to lift fault
triangle quality (min-angle / tri_q), then snap the top boundary to z = 0.

PLAN: meshing/docs/PLAN_fault_triangle_quality_2026-05-26.md (Phase 1).

Why this exists: gmsh embeds the cut fault STL *verbatim* as the Physical
Surface 101 triangulation (STL tri count == fault tri count), so fault
triangle quality == input STL quality.  The cut STL carries ~185 near-zero
clip needles along the z = 0 trace (ts_to_stl.py per-triangle clip).  This
tool isotropically remeshes the STL surface (MMG mmgs_O3), holding the open
free boundary (z = 0 trace + NW-cut + deep perimeter) and bounding the
geometric drift via mmgs -hausd, then re-snaps the top-boundary nodes to
exactly z = 0 so run_z0cut_meshing.py's free-surface invariant + trace embed
still hold.  Re-mesh the output with the EXISTING run_z0cut_meshing.py to a
new .msh (Phase 2 of the plan).

Engine: MMG mmgs_O3 5.8.x (primary).  --engine gmsh / pymeshlab are Phase-3
fallbacks (stubs here — implement only if mmgs misses the target).

Usage (conda activate pythonenv):
    python remesh_fault_stl.py INPUT.stl OUTPUT.stl
    python remesh_fault_stl.py INPUT.stl OUTPUT.stl --target-edge 500 --hausd 25
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import meshio
import numpy as np

DEFAULT_HAUSD = 25.0            # [m] max deviation from the input surface (mmgs -hausd)
DEFAULT_HGRAD = 1.3             # size gradation (mmgs -hgrad)
DEFAULT_RIDGE_ANGLE = 45.0     # [deg] ridge-detection angle (mmgs -ar)
DEFAULT_Z_SNAP_TOL = 1.0e-3    # [m] |z| < tol on the top boundary -> snap to exact 0
EPS_Z = 1.0e-6                 # [m] "on the z = 0 plane" (matches run_z0cut_meshing.py)
CORNER_ANGLE_TARGET = 25.0     # [deg] min-angle target for the surface gate report
DEFAULT_CGAL_ITERS = 8         # isotropic_remeshing iterations (flip+relax regularize)
# CGAL isotropic_remeshing only rebuilds edges outside [4/5, 4/3]*target; a
# target AT the current median is a near-no-op.  0.8*median lands the current
# edge distribution in the active band so the needles get regularized.  This
# is a tuning *fraction* (derived from the mesh's own median), not a hardcoded
# absolute size.  Empirically 0.8 clears worst-angle >= 25 / tri_q >= 0.6 on
# the 500 m SAFS fault; smaller fractions over-refine and re-seed boundary
# slivers (see the parameter sweep in IMPLEMENTATION_REPORT.md).
DEFAULT_CGAL_TARGET_FRACTION = 0.8


# ----------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------

def _read_triangles(mesh: meshio.Mesh) -> tuple[np.ndarray, np.ndarray]:
    """Return (points (N,3) float, tris (M,3) int) for the first triangle block."""
    tri = None
    for cb in mesh.cells:
        if cb.type == "triangle":
            tri = cb.data
            break
    if tri is None or np.asarray(tri).shape[0] == 0:
        raise ValueError("mesh has no (non-empty) triangle block")
    return np.asarray(mesh.points, dtype=float), np.asarray(tri, dtype=np.int64)


def write_medit(mesh_path: Path, pts: np.ndarray, tris: np.ndarray,
                required_edges: np.ndarray | None = None,
                required_verts: np.ndarray | None = None) -> None:
    """Write a Medit (.mesh) file, optionally with RequiredEdges /
    RequiredVertices sections so mmgs pins those entities (does not move,
    split, collapse, or swap them).

    NOTE (deviation from plan step 3): the plan wrote the .mesh with meshio,
    but meshio cannot emit Medit RequiredEdges/RequiredVertices, which are
    essential to keep the z = 0 trace exactly in-plane through mmgs (a free
    mmgs run lifts the trace up to ~45 m off z = 0).  STL is still *read*
    with meshio; only the medit *write* is explicit here.
    Triangle/edge indices are written 1-based per the Medit convention.
    """
    out: list[str] = ["MeshVersionFormatted 2", "Dimension 3", "", "Vertices",
                      str(pts.shape[0])]
    for x, y, z in pts:
        out.append(f"{x:.10g} {y:.10g} {z:.10g} 0")
    if required_edges is not None and required_edges.shape[0] > 0:
        out += ["", "Edges", str(required_edges.shape[0])]
        for a, b in required_edges:
            out.append(f"{int(a) + 1} {int(b) + 1} 1")
        out += ["", "RequiredEdges", str(required_edges.shape[0])]
        out += [str(k + 1) for k in range(required_edges.shape[0])]
    out += ["", "Triangles", str(tris.shape[0])]
    for a, b, c in tris:
        out.append(f"{int(a) + 1} {int(b) + 1} {int(c) + 1} 0")
    if required_verts is not None and required_verts.shape[0] > 0:
        out += ["", "RequiredVertices", str(required_verts.shape[0])]
        out += [str(int(v) + 1) for v in required_verts]
    out += ["", "End", ""]
    Path(mesh_path).write_text("\n".join(out))


def _z0_boundary(points: np.ndarray, tris: np.ndarray,
                 tol: float = EPS_Z) -> tuple[np.ndarray, np.ndarray]:
    """Return (z0_edges (K,2), z0_nodes (P,)) — the free-boundary edges whose
    both endpoints lie on the z = 0 plane, and the unique nodes they touch."""
    be = _boundary_edges(tris)
    if be.size == 0:
        return np.empty((0, 2), int), np.empty((0,), int)
    on = (np.abs(points[be[:, 0], 2]) < tol) & (np.abs(points[be[:, 1], 2]) < tol)
    z0e = be[on]
    return z0e, np.unique(z0e.reshape(-1)) if z0e.size else np.empty((0,), int)


def stl_to_medit(stl: Path, mesh_path: Path, z_tol: float = EPS_Z) -> dict:
    """Read an STL with meshio (vertices merged on read), find the z = 0 trace
    boundary, and write a Medit .mesh that pins that trace (RequiredEdges +
    RequiredVertices) so mmgs keeps it exactly in-plane.  Returns a stats dict."""
    m = meshio.read(str(stl))
    pts, tris = _read_triangles(m)
    z0e, z0n = _z0_boundary(pts, tris, z_tol)
    write_medit(mesh_path, pts, tris, required_edges=z0e, required_verts=z0n)
    return {
        "n_verts": int(pts.shape[0]),
        "n_tris": int(tris.shape[0]),
        "n_z0_trace_edges": int(z0e.shape[0]),
        "n_z0_trace_nodes": int(z0n.shape[0]),
        "zmin": float(pts[:, 2].min()),
        "zmax": float(pts[:, 2].max()),
        "bbox": (float(pts[:, 0].min()), float(pts[:, 0].max()),
                 float(pts[:, 1].min()), float(pts[:, 1].max()),
                 float(pts[:, 2].min()), float(pts[:, 2].max())),
    }


def _read_medit_tris(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Minimal Medit (.mesh) reader for the Vertices + Triangles blocks only.

    mmgs output carries extra sections (Edges, Ridges, RequiredEdges,
    Corners, Normals); this parser ignores them so we never depend on
    meshio's medit reader handling every mmgs section.  Triangle indices are
    converted from Medit's 1-based to 0-based.
    """
    lines = Path(path).read_text().split("\n")
    pts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    i = 0
    n = len(lines)
    while i < n:
        tok = lines[i].strip()
        if tok == "Vertices":
            cnt = int(lines[i + 1].split()[0])
            i += 2
            for k in range(cnt):
                p = lines[i + k].split()
                pts.append((float(p[0]), float(p[1]), float(p[2])))
            i += cnt
            continue
        if tok == "Triangles":
            cnt = int(lines[i + 1].split()[0])
            i += 2
            for k in range(cnt):
                p = lines[i + k].split()
                tris.append((int(p[0]) - 1, int(p[1]) - 1, int(p[2]) - 1))
            i += cnt
            continue
        i += 1
    if not pts or not tris:
        raise ValueError(f"medit file {path} has no Vertices/Triangles blocks")
    return np.asarray(pts, dtype=float), np.asarray(tris, dtype=np.int64)


def _read_remeshed(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read mmgs Medit output (Vertices + Triangles).

    NOTE (deviation from plan step 5): mmgs output always carries Medit
    sections meshio's reader rejects (RequiredEdges, Ridges, Normals,
    Tangents), so we use the explicit `_read_medit_tris` parser rather than
    meshio for the mmgs output.  meshio is still used to *read the input STL*.
    """
    return _read_medit_tris(path)


def _snap_and_write_stl(pts: np.ndarray, tris: np.ndarray, stl_out: Path,
                        z_snap_tol: float) -> dict:
    """Snap near-z=0 nodes to exact z = 0, assert the free-surface invariant,
    and write a binary STL.  Both engines pin the trace to z = 0 (mmgs via
    RequiredEdges, CGAL via constrained-border relax), so this snap only
    cleans float roundoff; an above-band node is a hard error."""
    on_top = (pts[:, 2] > -z_snap_tol) & (pts[:, 2] <= z_snap_tol)
    n_snapped = int(on_top.sum())
    pts[on_top, 2] = 0.0
    zmax = float(pts[:, 2].max())
    if zmax > z_snap_tol:
        raise ValueError(
            f"remesh lifted the surface above z = 0 (max z = {zmax:.6g} m > "
            f"z_snap_tol = {z_snap_tol:g} m); this would break the "
            f"free-surface invariant.  Reduce the target / check the input STL.")
    pts[:, 2] = np.where(np.abs(pts[:, 2]) < z_snap_tol, 0.0, pts[:, 2])
    meshio.write(str(stl_out),
                 meshio.Mesh(pts, [("triangle", tris)]),
                 file_format="stl")
    return {
        "n_verts": int(pts.shape[0]),
        "n_tris": int(tris.shape[0]),
        "n_snapped_to_z0": n_snapped,
        "zmin": float(pts[:, 2].min()),
        "zmax": float(pts[:, 2].max()),
        "points": pts,
        "tris": tris,
    }


def medit_to_stl_snapped(mesh_out: Path, stl_out: Path,
                         z_snap_tol: float) -> dict:
    """Read the remeshed Medit mesh (mmgs) and finalize to a snapped STL."""
    pts, tris = _read_remeshed(mesh_out)
    return _snap_and_write_stl(pts, tris, stl_out, z_snap_tol)


def write_off(off_path: Path, pts: np.ndarray, tris: np.ndarray) -> None:
    """Write a double-precision OFF (CGAL input) from merged arrays."""
    meshio.write(str(off_path),
                 meshio.Mesh(np.asarray(pts, float), [("triangle", tris)]),
                 file_format="off")


# ----------------------------------------------------------------------------
# Geometry / quality metrics (vectorized; formulas match
# locate_fault_slivers.py and check_mesh_quality.py)
# ----------------------------------------------------------------------------

def median_edge_length(points: np.ndarray, tris: np.ndarray) -> float:
    """Median 3D edge length over the unique (deduplicated) edge set.

    Matches PLAN_retriangulate.reference_edge_length semantics: each undirected
    edge counted once.
    """
    e = np.stack([tris[:, [0, 1, 2]], tris[:, [1, 2, 0]]], axis=2)  # (M,3,2)
    e = np.sort(e, axis=2)                                          # (min,max) per edge
    uniq = np.unique(e.reshape(-1, 2), axis=0)
    L = np.linalg.norm(points[uniq[:, 0]] - points[uniq[:, 1]], axis=1)
    return float(np.median(L))


def tri_metrics(points: np.ndarray, tris: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (min_angle_deg (M,), tri_q (M,)).

    tri_q = 4*sqrt(3)*A / sum(edge^2)  in [0,1]  (1 = equilateral).
    """
    p = points[tris]  # (M,3,3)

    def _ang(u, v):
        nu = np.linalg.norm(u, axis=1)
        nv = np.linalg.norm(v, axis=1)
        c = np.einsum("ij,ij->i", u, v) / (nu * nv + 1e-30)
        return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))

    a0 = _ang(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    a1 = _ang(p[:, 2] - p[:, 1], p[:, 0] - p[:, 1])
    a2 = _ang(p[:, 0] - p[:, 2], p[:, 1] - p[:, 2])
    min_angle = np.minimum(np.minimum(a0, a1), a2)

    e0 = np.linalg.norm(p[:, 1] - p[:, 0], axis=1)
    e1 = np.linalg.norm(p[:, 2] - p[:, 1], axis=1)
    e2 = np.linalg.norm(p[:, 0] - p[:, 2], axis=1)
    sumL2 = e0 ** 2 + e1 ** 2 + e2 ** 2
    area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tri_q = np.where(sumL2 > 0.0, 4.0 * np.sqrt(3.0) * area / sumL2, 0.0)
    return min_angle, tri_q


def _boundary_edges(tris: np.ndarray) -> np.ndarray:
    """Return (K,2) array of undirected boundary edges (used by exactly one tri)."""
    e = np.stack([tris[:, [0, 1, 2]], tris[:, [1, 2, 0]]], axis=2)  # (M,3,2)
    e = np.sort(e, axis=2)                                          # (min,max) per edge
    uniq, counts = np.unique(e.reshape(-1, 2), axis=0, return_counts=True)
    return uniq[counts == 1]


def n_boundary_edges_per_tri(tris: np.ndarray) -> np.ndarray:
    """Return (M,) count of boundary (free) edges on each triangle."""
    be = _boundary_edges(tris)
    bset = set(map(tuple, be.tolist()))
    out = np.zeros(tris.shape[0], dtype=int)
    for m, t in enumerate(tris):
        c = 0
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            if (min(int(x), int(y)), max(int(x), int(y))) in bset:
                c += 1
        out[m] = c
    return out


def z0_trace_segment_count(points: np.ndarray, tris: np.ndarray,
                           tol: float = EPS_Z) -> int:
    """Number of boundary edges with both endpoints on the z = 0 plane."""
    be = _boundary_edges(tris)
    if be.size == 0:
        return 0
    za = np.abs(points[be[:, 0], 2]) < tol
    zb = np.abs(points[be[:, 1], 2]) < tol
    return int((za & zb).sum())


def n_connected_components(points: np.ndarray, tris: np.ndarray) -> int:
    """Number of connected components of the surface (vertex graph from tri edges)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = points.shape[0]
    rows = np.concatenate([tris[:, 0], tris[:, 1], tris[:, 2]])
    cols = np.concatenate([tris[:, 1], tris[:, 2], tris[:, 0]])
    g = coo_matrix((np.ones(rows.shape[0]), (rows, cols)), shape=(n, n))
    _, labels = connected_components(g, directed=False)
    used = np.unique(tris)
    return int(np.unique(labels[used]).size)


# ----------------------------------------------------------------------------
# Hausdorff spot-check (sampled point-to-triangle distance)
# ----------------------------------------------------------------------------

def _closest_dist_point_tris(p: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Distance from point p (3,) to each triangle in T (k,3,3).

    Vectorized closest-point-on-triangle (Ericson, Real-Time Collision
    Detection).  The seven Voronoi regions partition the plane, so applying
    the region masks over the interior default (in any order) is exact.
    """
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    ab, ac = b - a, c - a
    ap = p[None, :] - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p[None, :] - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p[None, :] - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = va + vb + vc
    safe = np.where(np.abs(denom) > 0.0, denom, 1.0)
    v = vb / safe
    w = vc / safe
    res = a + v[:, None] * ab + w[:, None] * ac      # interior default
    mA = (d1 <= 0) & (d2 <= 0)
    res[mA] = a[mA]
    mB = (d3 >= 0) & (d4 <= d3)
    res[mB] = b[mB]
    mC = (d6 >= 0) & (d5 <= d6)
    res[mC] = c[mC]
    with np.errstate(divide="ignore", invalid="ignore"):
        vab = np.where((d1 - d3) != 0, d1 / (d1 - d3), 0.0)
    mAB = (va <= 0) & (d1 >= 0) & (d3 <= 0)
    res[mAB] = (a + vab[:, None] * ab)[mAB]
    with np.errstate(divide="ignore", invalid="ignore"):
        wac = np.where((d2 - d6) != 0, d2 / (d2 - d6), 0.0)
    mAC = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    res[mAC] = (a + wac[:, None] * ac)[mAC]
    den_bc = (d4 - d3) + (d5 - d6)
    with np.errstate(divide="ignore", invalid="ignore"):
        wbc = np.where(den_bc != 0, (d4 - d3) / den_bc, 0.0)
    mBC = (vc <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    res[mBC] = (b + wbc[:, None] * (c - b))[mBC]
    return np.linalg.norm(p[None, :] - res, axis=1)


def sampled_hausdorff(out_pts: np.ndarray, in_pts: np.ndarray,
                      in_tris: np.ndarray, n_sample: int = 2000,
                      k: int = 24, seed: int = 0) -> float:
    """Sampled one-sided Hausdorff: max over a random sample of output vertices
    of the distance to the nearest input triangle (candidates pruned by a
    cKDTree on input-triangle centroids)."""
    from scipy.spatial import cKDTree
    cen = in_pts[in_tris].mean(axis=1)
    tree = cKDTree(cen)
    rng = np.random.default_rng(seed)
    m = out_pts.shape[0]
    idx = rng.choice(m, size=min(n_sample, m), replace=False)
    P = out_pts[idx]
    kk = min(k, cen.shape[0])
    _, nbr = tree.query(P, k=kk)
    nbr = np.atleast_2d(nbr.reshape(P.shape[0], kk))
    worst = 0.0
    for i in range(P.shape[0]):
        d = _closest_dist_point_tris(P[i], in_pts[in_tris[nbr[i]]])
        worst = max(worst, float(d.min()))
    return worst


# ----------------------------------------------------------------------------
# Engine: mmgs
# ----------------------------------------------------------------------------

def run_mmgs(mesh_in: Path, mesh_out: Path, args) -> None:
    """Invoke mmgs_O3 for an isotropic surface remesh.  Raises SystemExit on a
    missing binary or a non-zero mmgs exit (mmgs stderr/stdout surfaced)."""
    mmgs = args.mmgs_bin or shutil.which("mmgs_O3")
    if not mmgs:
        # Fall back to a binary sitting next to the active interpreter
        # (so `pythonenv/bin/python remesh_fault_stl.py` works without PATH).
        cand = Path(sys.executable).parent / "mmgs_O3"
        if cand.exists():
            mmgs = str(cand)
    if not mmgs or not Path(mmgs).exists():
        print("ERROR: mmgs_O3 not found.  Expected on PATH (conda activate "
              "pythonenv) or pass --mmgs-bin "
              "/Users/chunhuizhao/miniforge/envs/pythonenv/bin/mmgs_O3",
              file=sys.stderr)
        raise SystemExit(2)
    cmd = [
        str(mmgs),
        "-in", str(mesh_in),
        "-out", str(mesh_out),
        "-hsiz", repr(float(args.target_edge)),
        "-hausd", repr(float(args.hausd)),
        "-hgrad", repr(float(args.hgrad)),
        "-ar", repr(float(args.ridge_angle)),
        "-v", "5" if args.verbose else "1",
    ]
    print("  mmgs:", " ".join(cmd))
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        print("ERROR: mmgs_O3 failed (exit %d).  This usually means a "
              "non-manifold or self-intersecting input surface." % exc.returncode,
              file=sys.stderr)
        print("----- mmgs stdout -----\n" + (exc.stdout or ""), file=sys.stderr)
        print("----- mmgs stderr -----\n" + (exc.stderr or ""), file=sys.stderr)
        raise SystemExit(exc.returncode)
    if args.verbose:
        print(res.stdout)
    if not mesh_out.exists():
        # mmgs sometimes appends ".o.mesh"; check and adopt it.
        alt = mesh_out.with_suffix(".o.mesh")
        if alt.exists():
            alt.replace(mesh_out)
        else:
            print(f"ERROR: mmgs reported success but no output at {mesh_out}",
                  file=sys.stderr)
            raise SystemExit(3)


def run_cgal(off_in: Path, off_out: Path, target_edge: float, iters: int,
             cgal_bin: Path | None) -> None:
    """Invoke the compiled CGAL isotropic-remesh tool (constrained borders,
    collapse + relax) for OFF in/out.  Raises SystemExit on a missing binary
    or a non-zero exit (stdout/stderr surfaced)."""
    binp = cgal_bin or (Path(__file__).resolve().parent / "remesh_cgal"
                        / "build" / "remesh_fault_cgal")
    if not Path(binp).exists():
        print(f"ERROR: CGAL binary not built: {binp}\n  Build it in the "
              f"cgal-61 env (see remesh_cgal/CMakeLists.txt header), or pass "
              f"--cgal-bin.", file=sys.stderr)
        raise SystemExit(2)
    cmd = [str(binp), str(off_in), str(off_out),
           repr(float(target_edge)), str(int(iters))]
    print("  cgal:", " ".join(cmd))
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: remesh_fault_cgal failed (exit {exc.returncode}).",
              file=sys.stderr)
        print("----- cgal stdout -----\n" + (exc.stdout or ""), file=sys.stderr)
        print("----- cgal stderr -----\n" + (exc.stderr or ""), file=sys.stderr)
        raise SystemExit(exc.returncode)
    print("    " + (res.stdout or "").strip().replace("\n", "\n    "))
    if not off_out.exists():
        print(f"ERROR: CGAL reported success but no output at {off_out}",
              file=sys.stderr)
        raise SystemExit(3)


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------

def report_surface_quality(points: np.ndarray, tris: np.ndarray,
                           angle_target: float = CORNER_ANGLE_TARGET) -> dict:
    """Compute + print the surface quality gate; classify sub-target tris into
    geometry-forced corners (>=2 boundary edges) vs clip needles/interior."""
    min_angle, tri_q = tri_metrics(points, tris)
    worst = float(min_angle.min())
    qmin = float(tri_q.min())
    bad = min_angle < angle_target
    n_bad = int(bad.sum())
    nbe = n_boundary_edges_per_tri(tris) if n_bad else np.zeros(tris.shape[0], int)
    corner = int((bad & (nbe >= 2)).sum())
    needle = n_bad - corner

    print(f"  fault/surface tris       : {tris.shape[0]}")
    print(f"  worst min-angle [deg]    : {worst:.2f}  (target >= {angle_target:.0f})")
    print(f"  tri_qmin                 : {qmin:.4f}  (target >= 0.5)")
    print(f"  tri_q median             : {float(np.median(tri_q)):.4f}")
    print(f"  tris min-angle < {angle_target:.0f} deg : {n_bad}")
    if n_bad:
        print(f"    geometry-forced corners (>=2 bnd edges): {corner}")
        print(f"    clip needles / interior (<=1 bnd edge) : {needle}")
        order = np.argsort(min_angle)
        cen = points[tris].mean(axis=1)
        print("    worst (min_angle, #bnd, centroid x,y,z):")
        for m in order[:min(8, n_bad)]:
            cx, cy, cz = cen[m]
            print(f"      {min_angle[m]:5.2f}  {nbe[m]}  "
                  f"({cx:.0f},{cy:.0f},{cz:.0f})")
    return {
        "worst_min_angle": worst, "tri_qmin": qmin,
        "tri_qmed": float(np.median(tri_q)),
        "n_bad": n_bad, "n_corner": corner, "n_needle": needle,
        "pass_angle": worst >= angle_target, "pass_q": qmin >= 0.5,
    }


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="input cut fault STL (reaches z=0)")
    ap.add_argument("output", type=Path, help="output remeshed STL (new file)")
    ap.add_argument("--engine", choices=("cgal", "mmgs", "gmsh", "pymeshlab"),
                    default="cgal",
                    help="remesh engine (default: cgal — the only one that "
                         "clears worst-angle>=25 / tri_q>=0.5; mmgs gets close "
                         "but cannot collapse the pinned z=0 trace segments)")
    ap.add_argument("--target-edge", type=float, default=None,
                    help="target edge [m]; default = input median 3D edge "
                         "length (mmgs -hsiz) or 0.8*median (cgal)")
    ap.add_argument("--cgal-iters", type=int, default=DEFAULT_CGAL_ITERS,
                    help="cgal isotropic_remeshing iterations (default %(default)s)")
    ap.add_argument("--cgal-bin", type=Path, default=None,
                    help="path to the compiled remesh_fault_cgal "
                         "(default: remesh_cgal/build/remesh_fault_cgal)")
    ap.add_argument("--hausd", type=float, default=DEFAULT_HAUSD,
                    help="mmgs -hausd: max deviation from the input surface [m] "
                         "(default %(default)s)")
    ap.add_argument("--hgrad", type=float, default=DEFAULT_HGRAD,
                    help="mmgs -hgrad size gradation (default %(default)s)")
    ap.add_argument("--ridge-angle", type=float, default=DEFAULT_RIDGE_ANGLE,
                    help="mmgs -ar ridge-detection angle [deg] (default %(default)s)")
    ap.add_argument("--z-snap-tol", type=float, default=DEFAULT_Z_SNAP_TOL,
                    help="snap top-boundary nodes with |z|<tol to exact 0 [m] "
                         "(default %(default)s)")
    ap.add_argument("--mmgs-bin", type=Path, default=None,
                    help="path to mmgs_O3 (default: shutil.which)")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="dir for intermediate .mesh files (default: tempdir)")
    ap.add_argument("--keep-temp", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    if not args.input.is_file():
        print(f"ERROR: input STL not found: {args.input}", file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.engine in ("gmsh", "pymeshlab"):
        raise NotImplementedError(
            f"--engine {args.engine} is an unused Phase-3 fallback (A/B); the "
            f"cgal engine meets the target.  See "
            f"PLAN_fault_triangle_quality_2026-05-26.md Phase 3.")

    workdir = Path(tempfile.mkdtemp(prefix="remesh_fault_")) \
        if args.workdir is None else args.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    # Read the input arrays once; reuse for target-edge, z=0 in-count, Hausdorff.
    ip, it = _read_triangles(meshio.read(str(args.input)))
    in_zmax = float(ip[:, 2].max())
    print(f"[1/4] read STL     : {args.input.name}")
    print(f"      in: {ip.shape[0]} verts, {it.shape[0]} tris, "
          f"z=[{float(ip[:,2].min()):.1f}, {in_zmax:.1f}]")
    if abs(in_zmax) > EPS_Z:
        print(f"WARNING: input STL zmax={in_zmax:.6g} != 0; run_z0cut_meshing.py "
              f"requires fault_zmax == 0.", file=sys.stderr)
    if args.target_edge is None:
        med = median_edge_length(ip, it)
        args.target_edge = (med if args.engine == "mmgs"
                            else round(DEFAULT_CGAL_TARGET_FRACTION * med, 3))
        print(f"      median edge = {med:.1f} m -> target_edge = "
              f"{args.target_edge:.1f} m ({args.engine})")

    if args.engine == "mmgs":
        mesh_in = workdir / "fault_in.mesh"
        mesh_out = workdir / "fault_out.mesh"
        in_stats = stl_to_medit(args.input, mesh_in, z_tol=EPS_Z)
        print(f"      pinned z=0 trace: {in_stats['n_z0_trace_edges']} edges, "
              f"{in_stats['n_z0_trace_nodes']} nodes (RequiredEdges/Vertices)")
        if in_stats["n_z0_trace_edges"] == 0:
            print("WARNING: no z=0 trace edges found to pin; mmgs may move the "
                  "top boundary off-plane.", file=sys.stderr)
        print(f"[2/4] mmgs remesh  : -hsiz {args.target_edge:.1f} "
              f"-hausd {args.hausd:.1f} -hgrad {args.hgrad} -ar {args.ridge_angle}")
        run_mmgs(mesh_in, mesh_out, args)
        print(f"[3/4] Medit -> STL : snap |z|<{args.z_snap_tol:g} to 0")
        out_stats = medit_to_stl_snapped(mesh_out, args.output, args.z_snap_tol)
    else:  # cgal
        off_in = workdir / "fault_in.off"
        off_out = workdir / "fault_out.off"
        write_off(off_in, ip, it)
        print(f"[2/4] CGAL remesh  : isotropic target={args.target_edge:.1f} "
              f"iters={args.cgal_iters} (constrained borders, collapse+relax)")
        run_cgal(off_in, off_out, args.target_edge, args.cgal_iters, args.cgal_bin)
        cp, ct = _read_triangles(meshio.read(str(off_out)))
        print(f"[3/4] OFF -> STL   : snap |z|<{args.z_snap_tol:g} to 0")
        out_stats = _snap_and_write_stl(cp, ct, args.output, args.z_snap_tol)
    pts, tris = out_stats["points"], out_stats["tris"]
    print(f"      out: {out_stats['n_verts']} verts, {out_stats['n_tris']} tris, "
          f"z=[{out_stats['zmin']:.1f}, {out_stats['zmax']:.1f}], "
          f"snapped {out_stats['n_snapped_to_z0']} nodes to z=0")

    # ---- structural checks (must hold for run_z0cut_meshing.py) ----
    ncomp = n_connected_components(pts, tris)
    z0_out = z0_trace_segment_count(pts, tris)
    z0_in = z0_trace_segment_count(ip, it)
    print(f"[4/4] checks")
    print(f"      connected components     : {ncomp}  (want 1)")
    print(f"      z=0 trace segments out/in: {z0_out} / {z0_in}")
    if ncomp != 1:
        print(f"WARNING: output STL has {ncomp} connected components; "
              f"run_z0cut_meshing.py needs exactly one fault surface.",
              file=sys.stderr)
    if z0_out < 1:
        print("WARNING: output STL has NO z=0 free-boundary edges; the trace "
              "embed will fail (empty volume).", file=sys.stderr)
    n_tris_in = int(it.shape[0])
    ratio = out_stats["n_tris"] / max(1, n_tris_in)
    if ratio < 0.5 or ratio > 2.0:
        print(f"WARNING: tri count changed by {ratio:.2f}x "
              f"({n_tris_in} -> {out_stats['n_tris']}); check "
              f"--target-edge / iters.", file=sys.stderr)

    # ---- geometry fidelity ----
    haus = sampled_hausdorff(pts, ip, it)
    print(f"      sampled Hausdorff out->in: {haus:.2f} m  (budget {args.hausd:.1f} m)")
    if haus > args.hausd * 1.5:
        print(f"WARNING: sampled Hausdorff {haus:.2f} m exceeds 1.5x the "
              f"--hausd budget; the surface may have drifted.", file=sys.stderr)

    # ---- surface quality gate (report; the 2 NW-tip corners may not clear) ----
    print("  --- surface quality gate ---")
    q = report_surface_quality(pts, tris)

    print(f"\nwrote {args.output}")
    print("NEXT (Phase 2): re-mesh with the EXISTING mesher, e.g.")
    print("  python <meshing>/code/run_z0cut_meshing.py "
          f"--stl {args.output} --out <experimental>/<name>_triq.msh")
    if not args.keep_temp and args.workdir is None:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print(f"(intermediate Medit files kept in {workdir})")

    # exit 0 if the remesh produced a structurally valid STL; the angle target
    # is reported (the 2 geometry-forced NW-tip corners are expected to remain
    # < target — see the plan).  Hard-fail only on broken structure.
    if ncomp != 1 or z0_out < 1 or out_stats["zmax"] > args.z_snap_tol:
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
