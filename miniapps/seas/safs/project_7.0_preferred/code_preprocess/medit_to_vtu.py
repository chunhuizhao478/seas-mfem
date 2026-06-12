#!/usr/bin/env python3
"""medit_to_vtu.py — Split a CGAL Mesh_3 MEDIT (.mesh) output produced by
mesh_volume.cpp (with output_to_medit(rebind=true, show_patches=true,
all_cells=true)) into:

    <out_base>_bulk.vtu  — rock tetrahedra (single 'rock' attribute=1):
                              subdomain 1 (rock-side-A) ∪ subdomain 2 (rock-
                              side-B) with edge_min ≥ MIN_EDGE, plus subdomain
                              0 "cavity filler" cells that pass the same
                              edge_min ≥ MIN_EDGE AND q_iso ≥ Q_MIN gate.
                           Excluded: any tet with edge_min < MIN_EDGE; subdomain-
                           0 Delaunay convex-hull buffer cells whose quality
                           does not meet Q_MIN.
    <out_base>_fault.vtu — all surface triangles whose patch ref is NOT the box
                           ref (identified automatically as the patch with the
                           most triangles).

Cell data attached:
    bulk : 'attribute' = 1 for every kept tet (single MFEM attribute);
           'subdomain' = original MEDIT ref (0/1/2) for debug;
           'q_iso'     = normalized isoperimetric tet quality;
           'edge_min'  = shortest tet edge [m].
    fault: 'patch'     = MEDIT ref (per-fault patch id).

Filter rationale:
    - Polyhedral_complex_mesh_domain_3 with N faults all using subdomain pair
      (1,2) leaves a fraction of cells unlabeled (subdomain 0) due to
      ambiguous adjacency around the polyline intersections.  Most of those
      "subdomain 0" cells are GOOD interior tets that just happened to be
      near an ambiguous patch.  A smaller fraction are slivers / Delaunay
      buffer cells outside the box envelope.
    - The two-gate filter (q_iso ≥ Q_MIN AND edge_min ≥ MIN_EDGE) cleanly
      separates the "good interior cells" from the "bad buffer cells", and
      empirically gives 100 % point-in-tet coverage of the box interior on
      the SAFS 2000 m fixture.

Usage:
    python medit_to_vtu.py <input.mesh> <out_base>
        [--manifest manifest.json]
        [--min-edge 100.0]                 # hard floor on edge_min
        [--region0-q-min 0.3]              # subdomain 0 quality gate
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import meshio


def tet_iso_q(P):
    """Normalized isoperimetric tet quality, vectorized over (N,4,3) tet array.
    Matches check_msh_quality.py:tet_iso_q.  q ≈ 1 for regular tet, 0 for
    fully degenerate."""
    a, b, c, d = P[:,0], P[:,1], P[:,2], P[:,3]
    V = np.abs(np.einsum('ij,ij->i', b - a, np.cross(c - a, d - a))) / 6.0
    A = (np.linalg.norm(np.cross(b - a, c - a), axis=1) +
         np.linalg.norm(np.cross(b - a, d - a), axis=1) +
         np.linalg.norm(np.cross(c - a, d - a), axis=1) +
         np.linalg.norm(np.cross(c - b, d - b), axis=1)) / 2.0
    q = np.zeros(len(P))
    good = (V > 1e-30) & (A > 0)
    q[good] = (12.0 * (3.0 * V[good]) ** (2.0/3.0)) / A[good] / (
              3.0 ** (1.0/3.0) * 6.0 ** (2.0/3.0)
              ) * 6.0 ** (2.0/3.0) * 3.0 ** (1.0/3.0) / 6.0 ** (2.0/3.0)
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_mesh", type=Path)
    ap.add_argument("out_base",   type=Path)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--min-edge", type=float, default=100.0,
                    help="hard floor on tet edge_min [m] (default 100; tets "
                         "with any edge below this are dropped)")
    ap.add_argument("--region0-q-min", type=float, default=0.3,
                    help="quality threshold for keeping subdomain-0 cells "
                         "(default 0.3 — matches the project bar; 0.05 keeps "
                         "marginal cells; 0.0 keeps everything except sub-edge-"
                         "min ones)")
    args = ap.parse_args()

    if not args.input_mesh.is_file():
        print(f"error: input not found: {args.input_mesh}", file=sys.stderr); return 1

    if args.manifest is None:
        project = args.input_mesh.parent.parent
        args.manifest = project / "data_corefined" / "manifest.json"
    if not args.manifest.is_file():
        print(f"error: manifest not found: {args.manifest}", file=sys.stderr); return 1
    m = json.loads(args.manifest.read_text())
    print(f"input        : {args.input_mesh}")
    print(f"manifest     : {args.manifest}")
    print(f"min-edge     : {args.min_edge} m")
    print(f"region0-q-min: {args.region0_q_min}")

    print(f"\nloading {args.input_mesh} ...")
    msh = meshio.read(str(args.input_mesh))
    pts = msh.points
    print(f"  {len(pts)} vertices")

    tetra_blocks = []; tet_refs = []
    tri_blocks   = []; tri_refs = []
    for cb, refs in zip(msh.cells, msh.cell_data['medit:ref']):
        if cb.type == "tetra":
            tetra_blocks.append(cb.data); tet_refs.append(refs)
        elif cb.type == "triangle":
            tri_blocks.append(cb.data); tri_refs.append(refs)
    if not tetra_blocks:
        print("error: no tetrahedra in input", file=sys.stderr); return 1
    tets = np.concatenate(tetra_blocks); refs_t = np.concatenate(tet_refs)
    tris = np.concatenate(tri_blocks) if tri_blocks else np.zeros((0,3), int)
    refs_s = np.concatenate(tri_refs) if tri_refs else np.zeros((0,), int)
    print(f"  {len(tets)} tets, {len(tris)} triangles")
    print(f"  tet ref histogram: {dict(sorted(Counter(refs_t.tolist()).items()))}")
    print(f"  tri ref histogram: {dict(sorted(Counter(refs_s.tolist()).items()))}")

    # ----- per-tet edge_min and isoperimetric quality -----
    P = pts[tets]
    e = np.stack([
        np.linalg.norm(P[:,0]-P[:,1], axis=1),
        np.linalg.norm(P[:,0]-P[:,2], axis=1),
        np.linalg.norm(P[:,0]-P[:,3], axis=1),
        np.linalg.norm(P[:,1]-P[:,2], axis=1),
        np.linalg.norm(P[:,1]-P[:,3], axis=1),
        np.linalg.norm(P[:,2]-P[:,3], axis=1),
    ], axis=1)
    emin = e.min(axis=1)
    qiso = tet_iso_q(P)

    # ----- F1 filter (rock + region 0 with q≥region0_q_min, all edge≥min-edge) -----
    is_rock     = (refs_t == 1) | (refs_t == 2)
    is_buffer_q = (refs_t == 0) & (qiso >= args.region0_q_min)
    keep_mask   = (is_rock | is_buffer_q) & (emin >= args.min_edge)

    n_rock_keep   = (is_rock & (emin >= args.min_edge)).sum()
    n_rock_drop_e = (is_rock & (emin <  args.min_edge)).sum()
    n_buffer_keep = (is_buffer_q & (emin >= args.min_edge)).sum()
    n_buffer_drop_q = ((refs_t == 0) & (qiso < args.region0_q_min)).sum()
    n_buffer_drop_e = ((refs_t == 0) & (qiso >= args.region0_q_min) & (emin < args.min_edge)).sum()

    print(f"\nbulk filter F1 (edge_min≥{args.min_edge}, region-0 q≥{args.region0_q_min}):")
    print(f"  KEPT: {keep_mask.sum()} tets")
    print(f"    rock (sub 1+2)              : {n_rock_keep}")
    print(f"    buffer (sub 0, q≥threshold) : {n_buffer_keep}")
    print(f"  DROPPED:")
    print(f"    rock with edge_min<{args.min_edge:.0f}: {n_rock_drop_e}")
    print(f"    buffer with q<{args.region0_q_min}     : {n_buffer_drop_q}")
    print(f"    buffer with edge_min<{args.min_edge:.0f}: {n_buffer_drop_e}")

    bulk_tets = tets[keep_mask]
    bulk_emin = emin[keep_mask]
    bulk_qiso = qiso[keep_mask]
    bulk_subd = refs_t[keep_mask].astype(np.int32)
    # Single MFEM attribute for the whole rock region.
    bulk_attr = np.ones(len(bulk_tets), dtype=np.int32)

    print()
    print(f"  bulk q_iso   : min={bulk_qiso.min():.4f}  med={np.median(bulk_qiso):.4f}")
    print(f"  bulk edge    : min={bulk_emin.min():.2f} m  med={np.median(bulk_emin):.0f} m")
    print(f"  q≥0.05 frac  : {(bulk_qiso >= 0.05).mean():.4f}")
    print(f"  q≥0.30 frac  : {(bulk_qiso >= 0.30).mean():.4f}")

    bulk_path = args.out_base.with_name(args.out_base.name + "_bulk.vtu")
    bulk_path.parent.mkdir(parents=True, exist_ok=True)
    bulk = meshio.Mesh(points=pts, cells=[("tetra", bulk_tets)],
                       cell_data={
                           "attribute": [bulk_attr],   # MFEM-friendly single tag
                           "subdomain": [bulk_subd],   # debug: 0/1/2 origin
                           "q_iso":     [bulk_qiso],
                           "edge_min":  [bulk_emin],
                       })
    bulk.write(str(bulk_path), binary=True)
    print(f"\nwrote {bulk_path}")

    # ----- fault: triangles with patch ref != box ref -----
    if len(tris) == 0:
        print("warning: no surface triangles"); return 0

    # Identify the BOX patch ref by largest triangle count
    # (the box has 6 faces each ~LC_NEAR-triangulated; faults are smaller).
    counts = Counter(refs_s.tolist())
    box_ref = max(counts, key=counts.get)
    print(f"\nidentified box surface ref: {box_ref} ({counts[box_ref]} triangles)")
    fault_mask = refs_s != box_ref
    fault_tris = tris[fault_mask]
    fault_refs = refs_s[fault_mask]
    fault_counts = Counter(fault_refs.tolist())
    print(f"fault refs (per-fault triangle counts): {dict(sorted(fault_counts.items()))}")
    print(f"total fault triangles: {len(fault_tris)}")

    if len(fault_tris) == 0:
        print("warning: no fault triangles found"); return 0

    fault_path = args.out_base.with_name(args.out_base.name + "_fault.vtu")
    fault = meshio.Mesh(points=pts, cells=[("triangle", fault_tris)],
                        cell_data={"patch": [fault_refs]})
    fault.write(str(fault_path), binary=True)
    print(f"wrote {fault_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
