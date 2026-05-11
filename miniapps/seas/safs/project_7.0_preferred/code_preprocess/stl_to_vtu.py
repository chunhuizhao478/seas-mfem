#!/usr/bin/env python3
"""stl_to_vtu.py — Convert ASCII STL fault meshes to VTU for visualization.

Reads a directory of `*.stl` files, writes a sibling `*.vtu` per file
with per-triangle quality (q_tri = 4*sqrt(3)*A / sum_e^2) and edge-length
arrays attached.

Usage:
    conda activate pythonenv
    python stl_to_vtu.py [--in-dir DIR] [--out-dir DIR]
"""

import argparse
import sys
from pathlib import Path

import meshio
import numpy as np


def tri_quality(verts: np.ndarray) -> np.ndarray:
    """Per-triangle quality q = 4*sqrt(3)*A / (e_a^2 + e_b^2 + e_c^2)."""
    a = verts[:, 0]
    b = verts[:, 1]
    c = verts[:, 2]
    ab2 = ((b - a) ** 2).sum(axis=1)
    bc2 = ((c - b) ** 2).sum(axis=1)
    ca2 = ((a - c) ** 2).sum(axis=1)
    sum_sq = ab2 + bc2 + ca2
    cross = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cross, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = 4.0 * np.sqrt(3.0) * area / np.where(sum_sq > 0, sum_sq, 1.0)
    return np.where(sum_sq > 0, q, 0.0)


def tri_edges(verts: np.ndarray) -> np.ndarray:
    """Min, median, max edge length per triangle."""
    a, b, c = verts[:, 0], verts[:, 1], verts[:, 2]
    e1 = np.linalg.norm(b - a, axis=1)
    e2 = np.linalg.norm(c - b, axis=1)
    e3 = np.linalg.norm(a - c, axis=1)
    edges = np.stack([e1, e2, e3], axis=1)
    return edges  # shape (N, 3)


def convert(stl_path: Path, vtu_path: Path) -> dict:
    m = meshio.read(str(stl_path))
    tris = m.cells_dict.get("triangle")
    if tris is None or len(tris) == 0:
        raise RuntimeError(f"no triangles in {stl_path}")
    pts = m.points
    tri_verts = pts[tris]                    # (N, 3, 3)
    q = tri_quality(tri_verts)
    edges = tri_edges(tri_verts)
    edge_min = edges.min(axis=1)
    edge_max = edges.max(axis=1)

    cell_data = {
        "q_tri":     [q],
        "edge_min":  [edge_min],
        "edge_max":  [edge_max],
    }
    out_mesh = meshio.Mesh(
        points=pts,
        cells=[("triangle", tris)],
        cell_data=cell_data,
    )
    meshio.write(str(vtu_path), out_mesh, file_format="vtu")

    return {
        "n_verts": int(len(pts)),
        "n_tris":  int(len(tris)),
        "q_min":   float(q.min()),
        "q_median":float(np.median(q)),
        "q_max":   float(q.max()),
        "edge_min":float(edge_min.min()),
        "edge_median": float(np.median(np.concatenate([edges[:,0], edges[:,1], edges[:,2]]))),
        "edge_max":float(edge_max.max()),
    }


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", type=Path,
                    default=project / "data_corefined",
                    help="dir of *.stl to convert")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output dir for *.vtu (default: same as in-dir)")
    ap.add_argument("--pattern", default="*.stl",
                    help="glob pattern (default '*.stl')")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.out_dir is None:
        args.out_dir = args.in_dir
    args.out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(args.in_dir.glob(args.pattern))
    if not files:
        print(f"error: no files match {args.pattern} in {args.in_dir}",
              file=sys.stderr)
        return 1

    print(f"{'file':<55} {'V':>6} {'F':>6} {'edge_min':>9} {'edge_med':>9} {'edge_max':>9} {'q_min':>7} {'q_med':>7}")
    for stl in files:
        vtu = args.out_dir / (stl.stem + ".vtu")
        s = convert(stl, vtu)
        print(f"{stl.name[:55]:<55} {s['n_verts']:>6} {s['n_tris']:>6}"
              f" {s['edge_min']:>9.1f} {s['edge_median']:>9.1f} {s['edge_max']:>9.1f}"
              f" {s['q_min']:>7.4f} {s['q_median']:>7.4f}")

    print(f"\nwrote {len(files)} VTU file(s) to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
