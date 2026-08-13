#!/usr/bin/env python3
"""puml_to_msh.py -- deployed SeisSol PUML/HDF5 -> Gmsh v2.2 ASCII with SAFS tags.

BC codes in the PUML byte-packed `boundary` -> SAFS physical tags:
    3 (dynamic rupture) -> 101 (fault)     [deduped: one triangle per face;
                                            interior fault faces appear on BOTH tets]
    1 (free surface)    -> 102 (top)
    5 (absorbing)       -> 104 (sides+bottom; msh_to_puml maps 104 -> 5)
Volume tets -> physical 1.  Point order is preserved 1:1 (metric .sol relies on it).
"""
import argparse
import sys
from pathlib import Path

import h5py
import meshio
import numpy as np

# SeisSol PUML face k -> local tet vertices
FACE_VERTS = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
BC_TO_TAG = {3: 101, 1: 102, 5: 104}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("puml", type=Path)
    ap.add_argument("out_msh", type=Path)
    args = ap.parse_args(argv)

    with h5py.File(args.puml, "r") as f:
        pts = f["geometry"][:].astype(np.float64)
        tets = f["connect"][:].astype(np.int64)
        boundary = f["boundary"][:].astype(np.int64)
    print(f"puml: {len(pts):,} nodes, {len(tets):,} tets")

    tris_by_tag = {t: [] for t in BC_TO_TAG.values()}
    seen_fault = set()
    for k in range(4):
        codes = (boundary >> (8 * k)) & 0xFF
        for bc, tag in BC_TO_TAG.items():
            e = np.flatnonzero(codes == bc)
            if not len(e):
                continue
            tri = tets[e][:, list(FACE_VERTS[k])]
            if bc == 3:                        # interior faces: dedupe both-side tags
                keep = []
                for row in tri:
                    key = tuple(sorted(row.tolist()))
                    if key not in seen_fault:
                        seen_fault.add(key)
                        keep.append(row)
                if keep:
                    tris_by_tag[tag].append(np.asarray(keep, dtype=np.int64))
            else:
                tris_by_tag[tag].append(tri)

    cells = [("tetra", tets)]
    physical = [np.ones(len(tets), dtype=int)]
    for tag in sorted(tris_by_tag):
        blocks = tris_by_tag[tag]
        if not blocks:
            continue
        tri = np.concatenate(blocks)
        cells.append(("triangle", tri))
        physical.append(np.full(len(tri), tag, dtype=int))
        print(f"  tag {tag}: {len(tri):,} triangles")
    if not seen_fault:
        print("ERROR: no fault (BC=3) faces found", file=sys.stderr)
        return 1

    mesh = meshio.Mesh(points=pts, cells=cells,
                       cell_data={"gmsh:physical": physical,
                                  "gmsh:geometrical": [p.copy() for p in physical]})
    args.out_msh.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(str(args.out_msh), mesh, file_format="gmsh22", binary=False)
    print(f"wrote {args.out_msh}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
