#!/usr/bin/env python3
"""locate_slivers.py — Identify sliver tets in the bulk mesh and write
a visualization VTU + summary table grouped by nearest fault.

Outputs:
    code_meshing/safs_multifault_box_<R>m_slivers.vtu
        — only the sliver tets (q_iso < threshold) with cell data
          'q_iso', 'edge_min', 'depth_z', 'nearest_fault_patch'
    stdout: per-fault sliver count + xy/z distribution

This file is meant for the user to inspect in ParaView and decide which
fault traces (if any) to collapse in the input STL.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import meshio
from scipy.spatial import cKDTree


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bulk",  type=Path,
                    default=project / "code_meshing" / "safs_multifault_box_2000m_bulk.vtu")
    ap.add_argument("--fault", type=Path,
                    default=project / "code_meshing" / "safs_multifault_box_2000m_fault.vtu")
    ap.add_argument("--out",   type=Path, default=None,
                    help="output sliver VTU (default = <bulk>_slivers.vtu)")
    ap.add_argument("--q-threshold",   type=float, default=0.05,
                    help="q_iso below which a tet is a sliver (default 0.05)")
    ap.add_argument("--edge-threshold", type=float, default=100.0,
                    help="edge_min below which a tet is a sliver (default 100 m)")
    args = ap.parse_args()

    if not args.bulk.is_file():
        print(f"error: bulk not found: {args.bulk}", file=sys.stderr); return 1
    if not args.fault.is_file():
        print(f"error: fault not found: {args.fault}", file=sys.stderr); return 1
    if args.out is None:
        args.out = args.bulk.with_name(args.bulk.stem.replace("_bulk", "") + "_slivers.vtu")

    bulk  = meshio.read(str(args.bulk))
    fault = meshio.read(str(args.fault))
    pts   = bulk.points
    tets  = bulk.cells_dict['tetra']
    qiso     = bulk.cell_data.get('q_iso',    [None])[0]
    edge_min = bulk.cell_data.get('edge_min', [None])[0]
    if qiso is None or edge_min is None:
        # Recompute if missing.
        P = pts[tets]
        e = np.stack([
            np.linalg.norm(P[:,0]-P[:,1], axis=1),
            np.linalg.norm(P[:,0]-P[:,2], axis=1),
            np.linalg.norm(P[:,0]-P[:,3], axis=1),
            np.linalg.norm(P[:,1]-P[:,2], axis=1),
            np.linalg.norm(P[:,1]-P[:,3], axis=1),
            np.linalg.norm(P[:,2]-P[:,3], axis=1),
        ], axis=1)
        edge_min = e.min(axis=1)
        a, b, c, d = P[:,0], P[:,1], P[:,2], P[:,3]
        V = np.abs(np.einsum('ij,ij->i', b-a, np.cross(c-a, d-a))) / 6.0
        A = (np.linalg.norm(np.cross(b-a, c-a), axis=1)
             + np.linalg.norm(np.cross(b-a, d-a), axis=1)
             + np.linalg.norm(np.cross(c-a, d-a), axis=1)
             + np.linalg.norm(np.cross(c-b, d-b), axis=1)) / 2.0
        qiso = np.zeros(len(P))
        good = (V > 1e-30) & (A > 0)
        qiso[good] = (12.0 * (3.0 * V[good]) ** (2.0/3.0)) / A[good] / (
                     3.0 ** (1.0/3.0) * 6.0 ** (2.0/3.0)
                     ) * 6.0 ** (2.0/3.0) * 3.0 ** (1.0/3.0) / 6.0 ** (2.0/3.0)

    # Sliver mask.
    sliver_mask = (qiso < args.q_threshold) | (edge_min < args.edge_threshold)
    n_sliv = int(sliver_mask.sum())
    print(f"\n=== Sliver detection ===")
    print(f"  bulk tets       : {len(tets)}")
    print(f"  q_iso threshold : {args.q_threshold}")
    print(f"  edge threshold  : {args.edge_threshold} m")
    print(f"  sliver tets     : {n_sliv}  ({n_sliv/len(tets)*100:.3f} %)")
    if n_sliv == 0:
        print("  (no slivers — nothing to do)")
        return 0

    sliver_tets = tets[sliver_mask]
    sliver_q    = qiso[sliver_mask]
    sliver_em   = edge_min[sliver_mask]
    centroids   = pts[sliver_tets].mean(axis=1)

    # ----- nearest fault patch for each sliver -----
    fault_tris    = fault.cells_dict['triangle']
    fault_patches = fault.cell_data['patch'][0]
    # Per-fault centroids (we'll find nearest patch per sliver via centroid kd-tree).
    fault_centroids = fault.points[fault_tris].mean(axis=1)
    tree_f = cKDTree(fault_centroids)
    nearest_dist, nearest_tri = tree_f.query(centroids, k=1)
    nearest_patch = fault_patches[nearest_tri]

    # ----- summary -----
    print(f"\n=== Sliver location ===")
    print(f"  centroid z range: [{centroids[:,2].min():.0f}, {centroids[:,2].max():.0f}] m")
    print(f"  distance to nearest fault triangle:")
    print(f"    min={nearest_dist.min():.2f}m  median={np.median(nearest_dist):.0f}m  "
          f"max={nearest_dist.max():.0f}m")
    print(f"  slivers within 100 m of fault: {(nearest_dist < 100).sum()}")
    print(f"  slivers within 500 m of fault: {(nearest_dist < 500).sum()}")
    print(f"  slivers within 5 km of fault : {(nearest_dist < 5000).sum()}")

    print(f"\n=== Slivers grouped by nearest fault patch ===")
    # Decode patch id to fault basename if available.
    print(f"  (patch IDs match fault.vtu 'patch' cell data)")
    cnt = Counter(nearest_patch.tolist())
    for pid, n in sorted(cnt.items()):
        msk = nearest_patch == pid
        z_msk = centroids[msk][:, 2]
        d_msk = nearest_dist[msk]
        q_msk = sliver_q[msk]
        em_msk = sliver_em[msk]
        print(f"  patch {pid}: {n:4d} slivers  "
              f"z=[{z_msk.min():7.0f}, {z_msk.max():7.0f}] m  "
              f"d_to_fault: med={np.median(d_msk):4.0f}m  "
              f"q_iso: min={q_msk.min():.2e} med={np.median(q_msk):.2e}  "
              f"edge_min: min={em_msk.min():.2e}m")

    # ----- depth bands -----
    print(f"\n=== Slivers by depth band ===")
    bands = [(0, "z >= 0"),
             (-100, "0 > z >= -100"),
             (-1000, "-100 > z >= -1km"),
             (-5000, "-1km > z >= -5km"),
             (-15000, "-5km > z >= -15km"),
             (-100000, "z < -15km"),]
    z_arr = centroids[:, 2]
    for i in range(len(bands)-1):
        lo, hi = bands[i+1][0], bands[i][0]
        n = ((z_arr >= lo) & (z_arr < hi)).sum() if i > 0 else (z_arr >= hi).sum()
        print(f"  {bands[i][1]:>20}: {n:4d}")
    n_below = (z_arr < bands[-1][0]).sum()  # below -15km
    print(f"  {bands[-1][1]:>20}: {n_below:4d}")

    # ----- print worst N slivers with coords -----
    worst_idx = np.argsort(sliver_q)[:20]
    print(f"\n=== Worst 20 slivers (sorted by q_iso) ===")
    print(f"  {'rank':>4}  {'q_iso':>11}  {'edge_min':>9}  "
          f"{'patch':>5}  {'d_fault':>7}  {'centroid (UTM)':>32}")
    for r, i in enumerate(worst_idx):
        c = centroids[i]
        print(f"  {r+1:>4}  {sliver_q[i]:>11.3e}  {sliver_em[i]:>9.2f}  "
              f"{int(nearest_patch[i]):>5d}  {nearest_dist[i]:>7.0f}  "
              f"({c[0]:>10.0f}, {c[1]:>10.0f}, {c[2]:>8.0f})")

    # ----- write per-sliver CSV for easy inspection -----
    csv_path = args.out.with_suffix('.csv')
    with open(csv_path, 'w') as f:
        f.write("idx,q_iso,edge_min,patch,d_fault,cx,cy,cz\n")
        for i in range(n_sliv):
            c = centroids[i]
            f.write(f"{i},{sliver_q[i]:.6e},{sliver_em[i]:.3f},"
                    f"{int(nearest_patch[i])},{nearest_dist[i]:.2f},"
                    f"{c[0]:.3f},{c[1]:.3f},{c[2]:.3f}\n")
    print(f"\nwrote {csv_path}")

    # ----- write sliver-only VTU for ParaView -----
    sliv_mesh = meshio.Mesh(
        points=pts,
        cells=[("tetra", sliver_tets)],
        cell_data={
            "q_iso":               [sliver_q],
            "edge_min":            [sliver_em],
            "depth_z":             [centroids[:, 2]],
            "dist_to_fault":       [nearest_dist],
            "nearest_fault_patch": [nearest_patch.astype(np.int32)],
        },
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sliv_mesh.write(str(args.out), binary=True)
    print(f"\nwrote {args.out}")
    print(f"  → load in ParaView alongside _bulk.vtu and _fault.vtu;")
    print(f"    color slivers by 'nearest_fault_patch' to see which fault each sits next to,")
    print(f"    or by 'depth_z' / 'q_iso' for spatial / quality distribution.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
