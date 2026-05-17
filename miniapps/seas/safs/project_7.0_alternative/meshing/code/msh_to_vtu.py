#!/usr/bin/env python3
"""
msh_to_vtu.py — Split a gmsh .msh into two .vtu files:
    <out_base>_bulk.vtu   — all tetrahedra in the named bulk physical group
    <out_base>_fault.vtu  — all triangles in the named fault physical group

Cell-data attached to the output:
    `gmsh:physical`     — physical-group tag of the source cell
    `gmsh:geometrical`  — geometric entity tag (when present)
    `quality`           — minimum scaled Jacobian (only for the bulk VTU; a
                          rough quality proxy ~ ratio of inradius to
                          circumradius scaled to [0, 1])

Both VTUs share the full original point cloud so ParaView's "Append
Datasets" can recombine them losslessly if needed.

Usage:
    python msh_to_vtu.py INPUT.msh OUTPUT_BASE
        [--bulk-name rock] [--fault-name fault]
        [--no-quality]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import meshio


def tet_quality(points, conn):
    """Compute a per-tetrahedron quality score in [0, 1].

    Uses (3 * 6^(2/3) * V^(2/3)) / sum(face_area) — equivalent to the
    "isoperimetric" tet quality (1 = regular tet, 0 = degenerate).  This
    is fast, robust, and what ParaView reports as "Aspect Frobenius
    inverse" for tets.
    """
    p = points[conn]                   # (N_tet, 4, 3)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    vol = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0

    def tri_area(u, v, w):
        return 0.5 * np.linalg.norm(np.cross(v - u, w - u), axis=1)

    a012 = tri_area(p[:, 0], p[:, 1], p[:, 2])
    a013 = tri_area(p[:, 0], p[:, 1], p[:, 3])
    a023 = tri_area(p[:, 0], p[:, 2], p[:, 3])
    a123 = tri_area(p[:, 1], p[:, 2], p[:, 3])
    s = a012 + a013 + a023 + a123

    # Regular tet of edge L: V = L^3/(6*sqrt(2)), face area = sqrt(3)/4 * L^2,
    # sum = sqrt(3) L^2.  ratio  V^(2/3) / sum = 1 / (6^(2/3) * sqrt(3))
    # Multiply by that constant -> regular tet returns 1.
    norm = 6.0 ** (2.0 / 3.0) * np.sqrt(3.0)
    q = np.where(s > 0, norm * vol ** (2.0 / 3.0) / s, 0.0)
    return q


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input_msh", type=Path, help="input gmsh .msh file")
    ap.add_argument("output_base", type=Path,
                    help="output path prefix; produces "
                         "<base>_bulk.vtu and <base>_fault.vtu")
    ap.add_argument("--bulk-name", default="rock",
                    help="name of the 3-D physical group to extract as bulk")
    ap.add_argument("--fault-name", default="fault",
                    help="name of the 2-D physical group to extract as fault")
    ap.add_argument("--no-quality", action="store_true",
                    help="do not compute the per-tet quality field")
    args = ap.parse_args()

    if not args.input_msh.is_file():
        print(f"error: input not found: {args.input_msh}", file=sys.stderr)
        return 1

    print(f"Reading {args.input_msh} ...")
    mesh = meshio.read(args.input_msh)
    print(f"  points: {len(mesh.points)}")
    print(f"  cell blocks:")
    for cb in mesh.cells:
        print(f"    {cb.type}: {len(cb.data)}")

    field_data = mesh.field_data or {}
    name_to_tag = {name: tag for name, (tag, _dim) in field_data.items()}
    bulk_tag = name_to_tag.get(args.bulk_name)
    fault_tag = name_to_tag.get(args.fault_name)
    if bulk_tag is None:
        print(f"error: physical group {args.bulk_name!r} not found in mesh "
              f"(have: {list(name_to_tag)})", file=sys.stderr)
        return 1
    if fault_tag is None:
        print(f"error: physical group {args.fault_name!r} not found in mesh "
              f"(have: {list(name_to_tag)})", file=sys.stderr)
        return 1
    print(f"  physical: {args.bulk_name!r} -> tag {bulk_tag}, "
          f"{args.fault_name!r} -> tag {fault_tag}")

    phys = mesh.cell_data.get("gmsh:physical")
    geom = mesh.cell_data.get("gmsh:geometrical")
    if phys is None:
        print("error: mesh has no 'gmsh:physical' cell data; was the .msh "
              "written with -format msh4 / has Physical groups defined?",
              file=sys.stderr)
        return 1

    bulk_blocks = []
    bulk_phys = []
    bulk_geom = []
    fault_blocks = []
    fault_phys = []
    fault_geom = []

    for i, cb in enumerate(mesh.cells):
        cell_phys = phys[i]
        cell_geom = geom[i] if geom is not None else None
        if cb.type == "tetra":
            mask = cell_phys == bulk_tag
            if mask.any():
                bulk_blocks.append(meshio.CellBlock("tetra", cb.data[mask]))
                bulk_phys.append(cell_phys[mask])
                if cell_geom is not None:
                    bulk_geom.append(cell_geom[mask])
        elif cb.type == "triangle":
            mask = cell_phys == fault_tag
            if mask.any():
                fault_blocks.append(meshio.CellBlock("triangle", cb.data[mask]))
                fault_phys.append(cell_phys[mask])
                if cell_geom is not None:
                    fault_geom.append(cell_geom[mask])

    if not bulk_blocks:
        print(f"error: no tetrahedra in physical group {args.bulk_name!r}",
              file=sys.stderr)
        return 1
    if not fault_blocks:
        print(f"error: no triangles in physical group {args.fault_name!r}",
              file=sys.stderr)
        return 1

    n_tets = sum(len(b.data) for b in bulk_blocks)
    n_fault = sum(len(b.data) for b in fault_blocks)
    print(f"  bulk  cells: {n_tets} tets")
    print(f"  fault cells: {n_fault} triangles")

    def edge_length_stats(label, blocks, kind):
        if kind == "tetra":
            edges = np.array([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
        else:
            edges = np.array([(0, 1), (1, 2), (2, 0)])
        all_lens = []
        for cb in blocks:
            p = mesh.points[cb.data]
            for ei, ej in edges:
                all_lens.append(np.linalg.norm(p[:, ei] - p[:, ej], axis=1))
        L = np.concatenate(all_lens)
        print(f"  {label} edge length: "
              f"min={L.min():.2f}  median={np.median(L):.2f}  "
              f"mean={L.mean():.2f}  max={L.max():.2f}")
        return L

    L_tet = edge_length_stats("bulk", bulk_blocks, "tetra")
    L_fault = edge_length_stats("fault", fault_blocks, "triangle")

    bulk_cell_data = {"gmsh:physical": bulk_phys}
    if bulk_geom:
        bulk_cell_data["gmsh:geometrical"] = bulk_geom
    if not args.no_quality:
        q_blocks = []
        for cb in bulk_blocks:
            q_blocks.append(tet_quality(mesh.points, cb.data))
        q_all = np.concatenate(q_blocks)
        bulk_cell_data["quality"] = q_blocks
        print(f"  tet quality: min={q_all.min():.4f}  "
              f"mean={q_all.mean():.4f}  max={q_all.max():.4f}")
        n_bad = int((q_all < 0.1).sum())
        n_med = int(((q_all >= 0.1) & (q_all < 0.3)).sum())
        n_ok = int((q_all >= 0.3).sum())
        print(f"             q<0.1: {n_bad} ({100*n_bad/len(q_all):.2f}%)  "
              f"0.1<=q<0.3: {n_med} ({100*n_med/len(q_all):.2f}%)  "
              f"q>=0.3: {n_ok}")

    fault_cell_data = {"gmsh:physical": fault_phys}
    if fault_geom:
        fault_cell_data["gmsh:geometrical"] = fault_geom

    bulk_mesh = meshio.Mesh(
        points=mesh.points, cells=bulk_blocks, cell_data=bulk_cell_data
    )
    fault_mesh = meshio.Mesh(
        points=mesh.points, cells=fault_blocks, cell_data=fault_cell_data
    )

    bulk_path = args.output_base.with_name(args.output_base.name + "_bulk.vtu")
    fault_path = args.output_base.with_name(args.output_base.name + "_fault.vtu")
    bulk_path.parent.mkdir(parents=True, exist_ok=True)
    bulk_mesh.write(bulk_path)
    fault_mesh.write(fault_path)
    print(f"Wrote {bulk_path} ({bulk_path.stat().st_size / 1024 / 1024:.2f} MB)")
    print(f"Wrote {fault_path} ({fault_path.stat().st_size / 1024 / 1024:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
