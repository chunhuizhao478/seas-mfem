#!/usr/bin/env python
"""reorient_negative_tets.py — flip tetrahedra with negative signed volume to
positive orientation (swap local vertices 2 and 3), in a Gmsh v2.2 .msh.

WHY: `remove_sliver_tets.py`'s edge-removal flips can leave a handful of tets
listed in NEGATIVE-volume orientation (valid geometry, wrong winding).  Such a
.msh reads as "inverted" and, fed to SeisSol directly, an inverted-tet mesh
blows the bulk energy up (Mw 88, 1e271 -> NaN at t~1s).  `msh_to_puml.py`
already re-orients when writing the PUML, so the HDF5 is safe regardless; this
tool makes the intermediate .msh / _safstags.msh consistent (0 inverted) too.

Swapping vertices 2<->3 negates the signed volume and leaves the tet's four
FACES (as node-sets) unchanged, so fault/boundary face tags and embedding are
preserved exactly.  Verifies 0 negative-volume tets on output.

Usage:
    conda activate pythonenv
    python reorient_negative_tets.py IN.msh OUT.msh
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import fastmsh as meshio   # memory-lean drop-in: meshio.read OOMs at 48.9M tets


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", type=Path)
    ap.add_argument("outfile", type=Path)
    ap.add_argument("--binary", action="store_true",
                    help="gmsh22 BINARY intermediate")
    args = ap.parse_args(argv)

    m = meshio.read(str(args.infile))
    pts = np.asarray(m.points, dtype=np.float64)
    n_flipped = 0
    for ci, cb in enumerate(m.cells):
        if cb.type != "tetra":
            continue
        T = np.asarray(cb.data).copy()
        a, b, c, d = pts[T[:, 0]], pts[T[:, 1]], pts[T[:, 2]], pts[T[:, 3]]
        V = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
        inv = V < 0
        T[inv, 2], T[inv, 3] = T[inv, 3], T[inv, 2].copy()
        m.cells[ci] = meshio.CellBlock("tetra", T)
        n_flipped += int(inv.sum())
        # verify
        a, b, c, d = pts[T[:, 0]], pts[T[:, 1]], pts[T[:, 2]], pts[T[:, 3]]
        V2 = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
        n_neg = int((V2 <= 0).sum())
        if n_neg:
            print(f"ERROR: {n_neg} tets still <= 0 volume after re-orient "
                  "(degenerate?)", file=sys.stderr)
            return 2

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(str(args.outfile), m, file_format="gmsh22", binary=args.binary)
    print(f"re-oriented {n_flipped} negative-volume tets -> {args.outfile} "
          f"(0 inverted)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
