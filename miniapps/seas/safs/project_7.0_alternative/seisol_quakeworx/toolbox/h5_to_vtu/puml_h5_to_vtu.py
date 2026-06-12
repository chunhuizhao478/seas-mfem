#!/usr/bin/env python3
"""Convert a SeisSol PUML mesh (.puml.h5) to ParaView-readable VTU files.

Outputs (into --out-dir):
    <stem>_bulk.vtu   all tetrahedra; cell data: group (material region),
                      n_fault_faces (number of BC-3 faces of the tet, handy
                      for thresholding near-fault elements)
    <stem>_fault.vtu  deduplicated dynamic-rupture (BC 3) triangles
    with --all-bcs additionally one surface VTU per other boundary code
    present (e.g. <stem>_bc1_free_surface.vtu, <stem>_bc5_absorbing.vtu).

PUML format (pumgen output):
    geometry  (nVert, 3) float64
    connect   (nElem, 4) vertex ids
    group     (nElem,)   material region id
    boundary  (nElem,)   per-face boundary codes packed into one integer:
                         8 bits/face for boundary-format 'i32',
                         16 bits/face for 'i64'; or already-unpacked (nElem,4).

Face k of a tet uses local vertices (PUML Numbering<TETRAHEDRON>):
    face 0: (1,0,2)   face 1: (0,1,3)   face 2: (1,2,3)   face 3: (2,0,3)
Verified against safs_mesh.puml.h5: with this ordering all 34,096 BC-1
(free-surface) faces lie exactly on the z=0 plane; the lexicographic
ordering misplaces 1,856 of them. The script re-runs a relaxed version of
this sanity check on every conversion and aborts on mismatch.

Boundary codes (SeisSol convention): 0 interior, 1 free surface,
3 dynamic rupture (fault), 5 absorbing.

Usage:
    python3 puml_h5_to_vtu.py <mesh.puml.h5> [--out-dir DIR] [--all-bcs]

Requires only numpy + h5py (no vtk/meshio).
"""
import argparse
import os
import struct
import sys

import h5py
import numpy as np

# PUML Numbering<TETRAHEDRON>::facevertices() — see module docstring.
FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))

BC_NAMES = {1: "free_surface", 3: "fault", 5: "absorbing"}

VTK_TETRA = 10
VTK_TRIANGLE = 5


# ---------------------------------------------------------------- VTU writer
def write_vtu(path, points, cells, cell_type, cell_data=None):
    """Write one XML VTU (appended raw binary, little endian, UInt64 headers).

    points: (nPts,3) float64; cells: (nCells,k) int; cell_data: {name: (nCells,)}.
    """
    points = np.ascontiguousarray(points, dtype="<f8")
    cells = np.ascontiguousarray(cells, dtype="<i8")
    n_pts, n_cells = len(points), len(cells)
    k = cells.shape[1]
    offsets = np.arange(1, n_cells + 1, dtype="<i8") * k
    types = np.full(n_cells, cell_type, dtype="u1")
    cell_data = cell_data or {}

    blocks = []   # (xml name, vtk type string, ncomp, raw bytes)
    blocks.append(("Points:Points", "Float64", 3, points.tobytes()))
    blocks.append(("Cells:connectivity", "Int64", 1, cells.tobytes()))
    blocks.append(("Cells:offsets", "Int64", 1, offsets.tobytes()))
    blocks.append(("Cells:types", "UInt8", 1, types.tobytes()))
    cd_blocks = []
    for name, arr in cell_data.items():
        arr = np.ascontiguousarray(arr, dtype="<i4")
        cd_blocks.append((name, "Int32", 1, arr.tobytes()))

    # appended-data offsets: each block is <UInt64 nbytes><raw bytes>
    raw_all = blocks + cd_blocks
    offs, pos = [], 0
    for _, _, _, raw in raw_all:
        offs.append(pos)
        pos += 8 + len(raw)

    def da(tag_name, vtype, ncomp, off):
        comp = f' NumberOfComponents="{ncomp}"' if ncomp > 1 else ""
        return (f'        <DataArray type="{vtype}" Name="{tag_name}"{comp}'
                f' format="appended" offset="{off}"/>\n')

    xml = []
    xml.append('<?xml version="1.0"?>\n')
    xml.append('<VTKFile type="UnstructuredGrid" version="1.0" '
               'byte_order="LittleEndian" header_type="UInt64">\n')
    xml.append('  <UnstructuredGrid>\n')
    xml.append(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n')
    xml.append('      <Points>\n')
    xml.append(da("Points", "Float64", 3, offs[0]))
    xml.append('      </Points>\n')
    xml.append('      <Cells>\n')
    xml.append(da("connectivity", "Int64", 1, offs[1]))
    xml.append(da("offsets", "Int64", 1, offs[2]))
    xml.append(da("types", "UInt8", 1, offs[3]))
    xml.append('      </Cells>\n')
    if cd_blocks:
        xml.append('      <CellData>\n')
        for (name, vtype, ncomp, _), off in zip(cd_blocks, offs[4:]):
            xml.append(da(name, vtype, ncomp, off))
        xml.append('      </CellData>\n')
    xml.append('    </Piece>\n')
    xml.append('  </UnstructuredGrid>\n')
    xml.append('  <AppendedData encoding="raw">\n_')

    with open(path, "wb") as fh:
        fh.write("".join(xml).encode())
        for _, _, _, raw in raw_all:
            fh.write(struct.pack("<Q", len(raw)))
            fh.write(raw)
        fh.write(b"\n  </AppendedData>\n</VTKFile>\n")


# ------------------------------------------------------------- PUML reading
def unpack_boundary(bnd, fmt):
    """Return (nElem, 4) int array of per-face BC codes."""
    bnd = np.asarray(bnd)
    if bnd.ndim == 2:            # already unpacked
        return bnd.astype(np.int64)
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if bits is None:
        sys.exit(f"ERROR: unsupported boundary-format {fmt!r} (expected i32/i64)")
    b = bnd.astype(np.int64)
    mask = (1 << bits) - 1
    return np.stack([(b >> (bits * k)) & mask for k in range(4)], axis=1)


def compact_surface(geom, tris):
    """Renumber triangle vertices to a compact point array."""
    used, inv = np.unique(tris.ravel(), return_inverse=True)
    return geom[used], inv.reshape(tris.shape)


def extract_bc_triangles(geom, conn, bc, code):
    """All faces with the given BC code, deduplicated (fault faces are
    flagged on both adjacent tets), as compacted (points, tris)."""
    tri_list = [conn[bc[:, k] == code][:, FACE_VERTS[k]] for k in range(4)]
    tris = np.concatenate(tri_list, axis=0)
    if len(tris) == 0:
        return None, None, 0
    _, keep = np.unique(np.sort(tris, axis=1), axis=0, return_index=True)
    n_raw = len(tris)
    tris = tris[np.sort(keep)]   # keep first occurrence, original orientation
    pts, tris_c = compact_surface(geom, tris)
    return pts, tris_c, n_raw


def sanity_check_face_ordering(geom, conn, bc):
    """Free-surface (BC 1) faces must be flat at z = z_max; a wrong face
    ordering picks wrong triangles of the boundary tets and breaks this."""
    sel = bc == 1
    if not sel.any():
        return  # mesh without a free surface: nothing to check
    ztop = geom[:, 2].max()
    zspan = geom[:, 2].max() - geom[:, 2].min()
    tol = 1e-6 * max(zspan, 1.0)
    bad = 0
    for k in range(4):
        tri = conn[sel[:, k]][:, FACE_VERTS[k]]
        bad += int(np.sum(~np.all(np.abs(geom[tri][:, :, 2] - ztop) < tol, axis=1)))
    if bad:
        sys.exit(f"ERROR: face-ordering sanity check failed — {bad} free-surface "
                 f"faces not on the z={ztop} plane. The PUML face numbering of "
                 "this file does not match FACE_VERTS; do not trust the output.")


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mesh", help="input .puml.h5 (pumgen output)")
    ap.add_argument("--out-dir", default=".", help="output directory (default: cwd)")
    ap.add_argument("--all-bcs", action="store_true",
                    help="also write one surface VTU per non-fault BC code")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    stem = os.path.basename(args.mesh)
    for suf in (".h5", ".puml"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]

    with h5py.File(args.mesh, "r") as f:
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
        group = f["group"][:].astype(np.int32)
        fmt = f.attrs.get("boundary-format", "i32")
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        bc = unpack_boundary(f["boundary"][:], fmt)

    print(f"mesh: {len(conn)} tets, {len(geom)} vertices, "
          f"boundary-format {fmt}, groups {np.unique(group).tolist()}")
    codes, counts = np.unique(bc, return_counts=True)
    print("face BC codes:", {int(c): int(n) for c, n in zip(codes, counts) if c})

    sanity_check_face_ordering(geom, conn, bc)

    # bulk: all tets
    n_fault_faces = np.sum(bc == 3, axis=1).astype(np.int32)
    bulk_path = os.path.join(args.out_dir, f"{stem}_bulk.vtu")
    write_vtu(bulk_path, geom, conn, VTK_TETRA,
              {"group": group, "n_fault_faces": n_fault_faces})
    print(f"wrote {bulk_path}  ({len(conn)} tets)")

    # fault: BC-3 triangles, deduplicated
    pts, tris, n_raw = extract_bc_triangles(geom, conn, bc, 3)
    if tris is None:
        print("no BC-3 (fault) faces in this mesh — fault VTU skipped")
    else:
        fault_path = os.path.join(args.out_dir, f"{stem}_fault.vtu")
        write_vtu(fault_path, pts, tris, VTK_TRIANGLE,
                  {"bc": np.full(len(tris), 3, np.int32)})
        print(f"wrote {fault_path}  ({len(tris)} triangles, "
              f"deduplicated from {n_raw} flagged faces)")

    if args.all_bcs:
        for code in codes:
            code = int(code)
            if code in (0, 3):
                continue
            pts, tris, n_raw = extract_bc_triangles(geom, conn, bc, code)
            name = BC_NAMES.get(code, f"code{code}")
            path = os.path.join(args.out_dir, f"{stem}_bc{code}_{name}.vtu")
            write_vtu(path, pts, tris, VTK_TRIANGLE,
                      {"bc": np.full(len(tris), code, np.int32)})
            print(f"wrote {path}  ({len(tris)} triangles)")


if __name__ == "__main__":
    main()
