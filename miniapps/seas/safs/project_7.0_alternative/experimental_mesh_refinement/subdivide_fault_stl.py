#!/usr/bin/env python3
"""subdivide_fault_stl.py — geometry-preserving 1->4 midpoint subdivision of a
cut fault STL.  Halves the edge length WITHOUT changing surface shape, the z = 0
trace, or the triangle min-angle distribution.

Why this exists: refining the near-fault resolution below the `_triq` median
(430 m) cannot be done by re-running the CGAL remesh.  A direct CGAL target
<= 250 m re-seeds sub-25 deg trace needles (the IMPLEMENTATION_REPORT.md CGAL
sweep: target 250 m -> worst-angle 23.3 deg, 4 bad tris; 200 m -> 21.5 deg).
Those needles are exactly the shallow-trace blow-up seed the `_triq` pipeline
removed, so CGAL cannot be pushed finer without giving them back.

A 1->4 midpoint split sidesteps that entirely: each parent triangle becomes 4
sub-triangles, every one *similar* to the parent (the central medial triangle
is the parent scaled by 1/2), so **every angle is preserved exactly** — the
worst min-angle stays 26.59 deg — while each edge halves (430 m -> ~215 m).
Midpoints lie on the parent facets, so the surface and the z = 0 trace are
geometrically exact (sampled/exact Hausdorff = 0).  The result is consumed by
run_z0cut_meshing.py exactly like any other cut fault STL (it reaches z = 0 as
an open free boundary, one connected component, watertight interior).

Midpoints are shared between the two triangles incident on each interior edge
(keyed by the sorted vertex-index pair), so the output stays watertight:
interior edges stay 2-incident, boundary/trace edges stay 1-incident and simply
double in count.

Usage (conda activate pythonenv):
    python subdivide_fault_stl.py INPUT_triq.stl OUTPUT_triqsubdiv.stl
    python subdivide_fault_stl.py INPUT.stl OUTPUT.stl --levels 1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import meshio
import numpy as np

# Reuse the project's canonical metric + STL-write helpers so the subdivided
# surface is measured with the EXACT formulas used by the quality gate
# (remesh_fault_stl.py lives in this same directory).
import remesh_fault_stl as R


def subdivide_1to4(points: np.ndarray,
                   tris: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One 1->4 midpoint subdivision of a triangle mesh.

    Returns (new_points, new_tris).  Midpoints are deduplicated per undirected
    edge (sorted vertex-index pair) so the surface stays watertight: every
    interior edge keeps exactly two incident triangles and no new non-manifold
    edge is created.  Sub-triangle winding matches the parent (outward normal
    preserved): for parent (a, b, c) with midpoints m_ab, m_bc, m_ca the four
    children are (a, m_ab, m_ca), (b, m_bc, m_ab), (c, m_ca, m_bc),
    (m_ab, m_bc, m_ca).
    """
    points = np.asarray(points, dtype=float)
    tris = np.asarray(tris, dtype=np.int64)
    if tris.ndim != 2 or tris.shape[1] != 3:
        raise ValueError(f"tris must be (M,3); got {tris.shape}")
    n_v = points.shape[0]
    n_t = tris.shape[0]

    # Unique undirected edges -> one new midpoint vertex each.  e is laid out as
    # [ (0,1) edges of all tris ; (1,2) edges ; (2,0) edges ], so the inverse
    # map reshapes to (3, M) rows = [m_ab, m_bc, m_ca].
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]],
                       axis=0)
    e = np.sort(e, axis=1)
    uniq, inv = np.unique(e, axis=0, return_inverse=True)
    inv = np.asarray(inv).reshape(-1)
    mid = 0.5 * (points[uniq[:, 0]] + points[uniq[:, 1]])
    new_points = np.vstack([points, mid])

    m = inv.reshape(3, n_t).T + n_v          # (M,3): cols = m_ab, m_bc, m_ca
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    m_ab, m_bc, m_ca = m[:, 0], m[:, 1], m[:, 2]
    new_tris = np.empty((4 * n_t, 3), dtype=np.int64)
    new_tris[0::4] = np.stack([a, m_ab, m_ca], axis=1)
    new_tris[1::4] = np.stack([b, m_bc, m_ab], axis=1)
    new_tris[2::4] = np.stack([c, m_ca, m_bc], axis=1)
    new_tris[3::4] = np.stack([m_ab, m_bc, m_ca], axis=1)
    return new_points, new_tris


def _max_edge_incidence(tris: np.ndarray) -> int:
    """Maximum number of triangles sharing any single undirected edge.

    2 for a clean manifold interior; 1 only edges are boundary; >2 means a
    non-manifold edge (would break the volume mesh)."""
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]],
                       axis=0)
    e = np.sort(e, axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return int(counts.max())


def report(points: np.ndarray, tris: np.ndarray, label: str) -> dict:
    """Print + return the surface metrics that must hold for the mesher."""
    med = R.median_edge_length(points, tris)
    ang, q = R.tri_metrics(points, tris)
    z0e, z0n = R._z0_boundary(points, tris)
    be = R._boundary_edges(tris)
    ncomp = R.n_connected_components(points, tris)
    n_bad = int((ang < R.CORNER_ANGLE_TARGET).sum())
    stats = {
        "n_verts": int(points.shape[0]),
        "n_tris": int(tris.shape[0]),
        "median_edge": med,
        "min_angle": float(ang.min()),
        "p10_angle": float(np.percentile(ang, 10)),
        "med_angle": float(np.median(ang)),
        "tri_qmin": float(q.min()),
        "tri_qmed": float(np.median(q)),
        "n_below_25": n_bad,
        "n_z0_edges": int(z0e.shape[0]),
        "n_boundary_edges": int(be.shape[0]),
        "n_components": ncomp,
        "zmin": float(points[:, 2].min()),
        "zmax": float(points[:, 2].max()),
    }
    print(f"  [{label}] verts={stats['n_verts']} tris={stats['n_tris']} "
          f"components={ncomp}")
    print(f"  [{label}] median_edge={med:.1f} m  "
          f"min_angle={stats['min_angle']:.2f} (p10={stats['p10_angle']:.2f}, "
          f"med={stats['med_angle']:.2f})  tri_qmin={stats['tri_qmin']:.4f}")
    print(f"  [{label}] tris<{R.CORNER_ANGLE_TARGET:.0f}deg={n_bad}  "
          f"z0_trace_edges={stats['n_z0_edges']}  "
          f"boundary_edges={stats['n_boundary_edges']}  "
          f"z=[{stats['zmin']:.1f},{stats['zmax']:.1f}]")
    return stats


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="input cut fault STL (reaches z=0)")
    ap.add_argument("output", type=Path, help="output subdivided STL (new file)")
    ap.add_argument("--levels", type=int, default=1,
                    help="number of 1->4 subdivisions (default 1: edge x0.5, "
                         "tris x4; 2: edge x0.25, tris x16)")
    ap.add_argument("--z-snap-tol", type=float, default=R.DEFAULT_Z_SNAP_TOL,
                    help="snap top-boundary nodes with |z|<tol to exact 0 [m] "
                         "(default %(default)s)")
    args = ap.parse_args(argv)

    if not args.input.is_file():
        print(f"ERROR: input STL not found: {args.input}", file=sys.stderr)
        return 1
    if args.levels < 1:
        print(f"ERROR: --levels must be >= 1 (got {args.levels})",
              file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] read STL     : {args.input.name}")
    ip, it = R._read_triangles(meshio.read(str(args.input)))
    in_zmax = float(ip[:, 2].max())
    if abs(in_zmax) > R.EPS_Z:
        print(f"WARNING: input STL zmax={in_zmax:.6g} != 0; "
              f"run_z0cut_meshing.py requires fault_zmax == 0.", file=sys.stderr)
    in_stats = report(ip, it, "in")

    print(f"[2/4] subdivide    : {args.levels} level(s) 1->4")
    pts, tris = ip, it
    for lv in range(args.levels):
        pts, tris = subdivide_1to4(pts, tris)
        print(f"      level {lv + 1}: {tris.shape[0]} tris, "
              f"{pts.shape[0]} verts")

    # All input z <= 0, so midpoints are z <= 0; the snap only cleans float
    # roundoff on the exact-0 trace nodes (their midpoints are bit-exact 0).
    print(f"[3/4] snap + write : |z|<{args.z_snap_tol:g} -> 0, binary STL")
    out_stats_io = R._snap_and_write_stl(pts, tris, args.output, args.z_snap_tol)
    pts, tris = out_stats_io["points"], out_stats_io["tris"]
    out_stats = report(pts, tris, "out")

    # ---- structural / invariant checks (must hold for run_z0cut_meshing.py) ----
    print(f"[4/4] checks")
    ok = True
    max_inc = _max_edge_incidence(tris)
    if max_inc > 2:
        print(f"ERROR: output has a non-manifold edge (max edge incidence "
              f"{max_inc} > 2); subdivision corrupted the surface.",
              file=sys.stderr)
        ok = False
    if out_stats["n_components"] != in_stats["n_components"]:
        print(f"ERROR: component count changed "
              f"{in_stats['n_components']} -> {out_stats['n_components']}; "
              f"subdivision must not split/merge the surface.", file=sys.stderr)
        ok = False
    if out_stats["n_components"] != 1:
        print(f"WARNING: output has {out_stats['n_components']} connected "
              f"components; run_z0cut_meshing.py needs exactly one fault "
              f"surface.", file=sys.stderr)
    # 1->4 doubles every edge: tris x4^L, z=0 trace edges x2^L, all boundary
    # edges x2^L.  Verify the exact ratios (a strong corruption tripwire).
    f = 4 ** args.levels
    s = 2 ** args.levels
    if out_stats["n_tris"] != f * in_stats["n_tris"]:
        print(f"ERROR: tri count {out_stats['n_tris']} != {f}x input "
              f"{in_stats['n_tris']}.", file=sys.stderr)
        ok = False
    if out_stats["n_z0_edges"] != s * in_stats["n_z0_edges"]:
        print(f"ERROR: z=0 trace edges {out_stats['n_z0_edges']} != {s}x input "
              f"{in_stats['n_z0_edges']}.", file=sys.stderr)
        ok = False
    if out_stats["n_boundary_edges"] != s * in_stats["n_boundary_edges"]:
        print(f"ERROR: boundary edges {out_stats['n_boundary_edges']} != {s}x "
              f"input {in_stats['n_boundary_edges']}.", file=sys.stderr)
        ok = False
    # Midpoint subdivision is exactly shape-preserving: the worst min-angle must
    # not move (beyond float noise).  This is the whole reason to subdivide
    # instead of CGAL-remeshing to 250 m.
    d_ang = abs(out_stats["min_angle"] - in_stats["min_angle"])
    if d_ang > 1.0e-3:
        print(f"ERROR: worst min-angle moved by {d_ang:.4g} deg "
              f"({in_stats['min_angle']:.4f} -> {out_stats['min_angle']:.4f}); "
              f"subdivision should preserve angles exactly.", file=sys.stderr)
        ok = False
    if out_stats["zmax"] > args.z_snap_tol:
        print(f"ERROR: output zmax={out_stats['zmax']:.6g} > z_snap_tol; "
              f"would break the free-surface invariant.", file=sys.stderr)
        ok = False

    print(f"      max edge incidence       : {max_inc}  (want 2)")
    print(f"      min-angle preserved       : "
          f"{in_stats['min_angle']:.2f} -> {out_stats['min_angle']:.2f} deg "
          f"(delta {d_ang:.2e})")
    print(f"      median edge               : "
          f"{in_stats['median_edge']:.1f} -> {out_stats['median_edge']:.1f} m")

    print(f"\nwrote {args.output}")
    print("NEXT: re-mesh with the EXISTING mesher (proven sliver-free recipe):")
    print(f"  python ../meshing/code/run_z0cut_meshing.py \\")
    print(f"      --stl {args.output.name} \\")
    print(f"      --out <name>_triqsubdiv.msh \\")
    print(f"      --lc-near 250 --surface-buffer-size 1500 "
          f"--surface-buffer-depth 1000")

    return 0 if ok else 4


if __name__ == "__main__":
    sys.exit(main())
