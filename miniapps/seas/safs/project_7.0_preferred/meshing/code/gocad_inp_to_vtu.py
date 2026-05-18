"""gocad_inp_to_vtu.py — split a SKUA-GOCAD-exported Abaqus .inp mesh
into one bulk VTU + one VTU per named *SURFACE so each fault / boundary
can be opened and toggled independently in ParaView.

Background
----------
GOCAD's `*SURFACE, TYPE=ELEMENT, NAME=<surface>` blocks point at one or
more `*ELSET, ELSET=<surface>_S{n}` element-side groups, where Sn is the
1-based Abaqus C3D4 face id (S1..S4).  meshio reads these into
``Mesh.cell_sets`` but does NOT expand them into triangle blocks, so
ParaView's native Abaqus reader sees the bulk tets but not the named
fault surfaces.  This script does the expansion using the standard
Abaqus C3D4 face -> node convention:

    S1: opposite node 4 -> (1, 2, 3)   -> 0-based (0, 1, 2)
    S2: opposite node 3 -> (1, 4, 2)   -> 0-based (0, 3, 1)
    S3: opposite node 2 -> (2, 4, 3)   -> 0-based (1, 3, 2)
    S4: opposite node 1 -> (3, 4, 1)   -> 0-based (2, 3, 0)

Output layout (next to the input ``.inp`` unless --out-dir is given):

    <stem>_bulk.vtu              all tets
    <stem>_<surface>.vtu         one per *SURFACE (trailing _C3D4 dropped
                                  from filename)

Usage
-----
::

    conda activate pythonenv
    python gocad_inp_to_vtu.py PATH_TO.inp [--out-dir DIR]

Recipe to inspect in ParaView
-----------------------------
1. File -> Open -> select all <stem>_*.vtu at once -> Apply.
2. Bulk: Representation = Wireframe (or Surface + Opacity ~ 0.15).
3. Each fault VTU: Representation = Surface, Coloring = Solid Color
   (pick a distinct color per fault).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import meshio
import numpy as np


# Abaqus C3D4 face -> (0-based) node indices.
# Reference: Abaqus Analysis User's Guide section on element libraries
# (the same convention every modern FE preprocessor exports).
FACE_NODES = {
    1: (0, 1, 2),   # S1 opposite node 4
    2: (0, 3, 1),   # S2 opposite node 3
    3: (1, 3, 2),   # S3 opposite node 2
    4: (2, 3, 0),   # S4 opposite node 1
}


def _build_block_offsets(m: meshio.Mesh) -> list[int]:
    offsets = [0]
    for cb in m.cells:
        offsets.append(offsets[-1] + len(cb.data))
    return offsets


def _cell_set_global_tet_ids(m: meshio.Mesh, offsets: list[int],
                              name: str) -> np.ndarray:
    """Return global tet indices (into the concatenated tet array) for
    a meshio cell set.  Only ``tetra`` entries are kept; non-tetra
    blocks in the cell set are silently dropped (GOCAD .inp tets only).
    """
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
    """For a *SURFACE NAME=<surface_name>, concat its S1..S4 ELSETs into
    one (N, 3) triangle array.  Returns None if no triangles were found
    (e.g. the surface name doesn't match any ELSET)."""
    tris = []
    for sn, nodes in FACE_NODES.items():
        es = f"{surface_name}_S{sn}"
        gids = _cell_set_global_tet_ids(m, offsets, es)
        if len(gids) == 0:
            continue
        tris.append(all_tetra[gids][:, list(nodes)])
    return np.concatenate(tris, axis=0) if tris else None


def _surface_names_from_inp_text(text: str) -> list[str]:
    """Return every NAME= from `*SURFACE, TYPE=ELEMENT, NAME=<name>`
    lines, in file order.  meshio does not expose *SURFACE definitions,
    so we re-parse the raw text."""
    return re.findall(r"\*SURFACE,\s*TYPE=ELEMENT,\s*NAME=(\S+)", text)


def split_inp(inp_path: Path, out_dir: Path | None = None,
              verbose: bool = True) -> dict:
    """Split an Abaqus .inp into bulk + per-surface VTUs.

    Parameters
    ----------
    inp_path : Path
        Path to the .inp file (GOCAD-style: C3D4 tets, *SURFACE blocks).
    out_dir : Path, optional
        Output directory (default: same dir as the input file).
    verbose : bool
        If True, print a per-surface summary while writing.

    Returns
    -------
    dict with keys ``bulk`` (Path) and ``surfaces`` (dict[name -> Path]).

    Raises
    ------
    FileNotFoundError if the input doesn't exist.
    ValueError       if the input has no tetra blocks.
    """
    inp_path = Path(inp_path)
    if not inp_path.is_file():
        raise FileNotFoundError(f"input .inp not found: {inp_path}")
    out_dir = Path(out_dir) if out_dir is not None else inp_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = inp_path.stem

    m = meshio.read(str(inp_path))
    pts = m.points
    tetra_blocks = [(i, cb.data) for i, cb in enumerate(m.cells)
                    if cb.type == "tetra"]
    if not tetra_blocks:
        raise ValueError(f"{inp_path}: no tetra blocks found")
    all_tetra = np.concatenate([d for _, d in tetra_blocks], axis=0)
    offsets = _build_block_offsets(m)

    # 1) Bulk: all tets, no cell_sets (keeps file small and ParaView fast).
    bulk_path = out_dir / f"{stem}_bulk.vtu"
    meshio.write(bulk_path,
                 meshio.Mesh(points=pts, cells=[("tetra", all_tetra)]),
                 file_format="vtu")
    if verbose:
        print(f"bulk -> {bulk_path.name}  "
              f"({bulk_path.stat().st_size/1e6:.1f} MB, "
              f"{len(all_tetra):,} tets)")

    # 2) Per *SURFACE: extract triangles, write to a separate VTU using
    #    only the points that surface actually references.
    surface_names = _surface_names_from_inp_text(inp_path.read_text())
    surface_paths: dict[str, Path] = {}
    if verbose:
        print(f"\nextracting {len(surface_names)} surface(s):")
    for n in surface_names:
        tris = _extract_surface_triangles(m, offsets, all_tetra, n)
        if tris is None or len(tris) == 0:
            if verbose:
                print(f"  {n:60s}  0 tris (skipped)")
            continue
        used = np.unique(tris.ravel())
        remap = -np.ones(pts.shape[0], dtype=np.int64)
        remap[used] = np.arange(len(used))
        tris_local = remap[tris]
        sub_stem = n[:-5] if n.endswith("_C3D4") else n  # drop trailing _C3D4
        out = out_dir / f"{stem}_{sub_stem}.vtu"
        meshio.write(out, meshio.Mesh(points=pts[used],
                                       cells=[("triangle", tris_local)]),
                     file_format="vtu")
        surface_paths[sub_stem] = out
        if verbose:
            print(f"  {sub_stem:60s}  {len(tris):>7,} tris  "
                  f"({out.stat().st_size/1e6:5.2f} MB)")

    return {"bulk": bulk_path, "surfaces": surface_paths}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inp", type=Path,
                    help="Path to the GOCAD-exported Abaqus .inp file")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output directory (default: same as .inp)")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress per-surface summary")
    args = ap.parse_args(argv)
    try:
        split_inp(args.inp, out_dir=args.out_dir, verbose=not args.quiet)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
