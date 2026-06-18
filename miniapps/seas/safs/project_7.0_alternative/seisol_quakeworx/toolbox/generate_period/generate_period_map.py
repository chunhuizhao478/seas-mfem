#!/usr/bin/env python3
"""Per-cell resolved-period map (grid size / local Vs) on the free surface.

Reads the SeisSol free-surface output mesh (safs-surface.xdmf -> raw
connect/geometry .bin), matches every surface triangle to its owning
boundary tetrahedron in the PUML mesh (.puml.h5) to get the true element
size, samples the local shear-wave speed Vs = sqrt(mu/rho) from the ASAGI
material grid (safs_material_cvm.nc, trilinear — the same `linear`
interpolation the easi !ASAGI block requests), and writes one triangle VTU
with cell data:

    h_tri_max_edge   [m]  longest edge of the surface triangle
    h_tet_max_edge   [m]  longest edge of the attached tetrahedron
                          (the mesh size that controls resolution; equals
                          h_tri_max_edge when --no-puml)
    h_ip             [m]  h_tet_max_edge / (degree+1): nominal spacing of
                          the degree+1 interpolation/integration points
                          along the longest edge (--degree, default 3 =
                          the CONVERGENCE_ORDER-4 QuakeWorx build).
                          NOTE: SeisSol's modal DG has no literal nodal
                          set; this is the standard SEM-style effective-
                          resolution normalisation (uniform average).
    vs               [m/s] Vs at the attached tet's centroid (surface-
                          triangle centroid when --no-puml)
    period           [s]  h_tet_max_edge / vs  — element-crossing time
    period_ip        [s]  h_ip / vs            — sub-point-crossing time

The two period fields are the SAME criterion under different
normalisations — do not double-count.  Minimum resolved period:

    T_min = epw     * period      with epw     ~ 2 elements/wavelength (O4)
          = ppw_pts * period_ip   with ppw_pts ~ 2*(degree+1) = 8 sub-
                                       points/wavelength (identical number)

This script deliberately leaves the epw/ppw factor out of both fields.

The material .nc is located via the `file:` entry of the easi material
yaml (resolved next to the yaml), or directly via --material-nc.

Defaults point at the v2_0_0_RSSRW case; everything is overridable:

    conda run -n pythonenv python3 generate_period_map.py
    ... --surface-xdmf X --puml Y --material-yaml Z --out OUT.vtu

Requires only numpy + h5py (no vtk/meshio/netCDF4) — the .nc is read as
plain HDF5.  PUML face numbering follows ../h5_to_vtu/puml_h5_to_vtu.py
(verified there against this very mesh: all 34,096 BC-1 faces flat at
z=0) and the same free-surface flatness sanity check is re-run here.
"""
import argparse
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# PUML Numbering<TETRAHEDRON>::facevertices() — see puml_h5_to_vtu.py.
FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))
TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
VTK_TRIANGLE = 5

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
QW_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # seisol_quakeworx/
DEF_XDMF = os.path.join(QW_DIR, "safs_seisol_v2_0_0_RSSRW_output", "safs-surface.xdmf")
DEF_PUML = os.path.join(QW_DIR, "safs_seisol_v2_0_0_RSSRW", "safs_mesh.puml.h5")
DEF_YAML = os.path.join(QW_DIR, "safs_seisol_v2_0_0_RSSRW", "safs_material_cvm.yaml")
DEF_OUT = os.path.join(SCRIPT_DIR, "safs_surface_period.vtu")


# ---------------------------------------------------------------- VTU writer
def write_vtu(path, points, cells, cell_type, cell_data=None):
    """One XML VTU (appended raw binary, little endian, UInt64 headers).

    Same writer as puml_h5_to_vtu.py, generalised so cell_data values may
    be float32/float64/int32 (dtype decides the DataArray type).
    """
    points = np.ascontiguousarray(points, dtype="<f8")
    cells = np.ascontiguousarray(cells, dtype="<i8")
    n_pts, n_cells = len(points), len(cells)
    k = cells.shape[1]
    offsets = np.arange(1, n_cells + 1, dtype="<i8") * k
    types = np.full(n_cells, cell_type, dtype="u1")
    cell_data = cell_data or {}

    vtk_type = {"f4": "Float32", "f8": "Float64", "i4": "Int32", "i8": "Int64"}

    blocks = [
        ("Points", "Float64", 3, points.tobytes()),
        ("connectivity", "Int64", 1, cells.tobytes()),
        ("offsets", "Int64", 1, offsets.tobytes()),
        ("types", "UInt8", 1, types.tobytes()),
    ]
    cd_blocks = []
    for name, arr in cell_data.items():
        arr = np.ascontiguousarray(arr)
        code = arr.dtype.str.lstrip("<>|=")
        if code not in vtk_type:
            arr = arr.astype("<f8")
            code = "f8"
        arr = arr.astype("<" + code)
        cd_blocks.append((name, vtk_type[code], 1, arr.tobytes()))

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


# --------------------------------------------- SeisSol surface-output reader
def read_surface_xdmf(xdmf_path):
    """Return (points (nV,3) f8, tris (nT,3) i8) from the first Uniform grid.

    SeisSol writes Format="Binary" raw little-endian files referenced
    relative to the .xdmf; topology/geometry are identical across steps.
    """
    base = os.path.dirname(os.path.abspath(xdmf_path))
    root = ET.parse(xdmf_path).getroot()
    grid = root.find(".//Grid[@GridType='Uniform']")
    if grid is None:
        sys.exit(f"ERROR: no Uniform grid found in {xdmf_path}")

    def read_item(item, np_kind):
        dims = tuple(int(d) for d in item.get("Dimensions").split())
        prec = int(item.get("Precision", "8"))
        if item.get("Format") != "Binary":
            sys.exit(f"ERROR: expected Format=Binary, got {item.get('Format')!r}")
        dtype = f"<{np_kind}{prec}"
        path = os.path.join(base, item.text.strip())
        arr = np.fromfile(path, dtype=dtype)
        if arr.size != np.prod(dims):
            sys.exit(f"ERROR: {path}: {arr.size} values, expected {np.prod(dims)}")
        return arr.reshape(dims)

    topo = grid.find("Topology")
    if topo.get("TopologyType") != "Triangle":
        sys.exit(f"ERROR: expected Triangle topology, got {topo.get('TopologyType')!r}")
    tris = read_item(topo.find("DataItem"), "i").astype(np.int64)
    pts = read_item(grid.find("Geometry").find("DataItem"), "f").astype(np.float64)
    return pts, tris


# ------------------------------------------------------------- PUML reading
def unpack_boundary(bnd, fmt):
    bnd = np.asarray(bnd)
    if bnd.ndim == 2:
        return bnd.astype(np.int64)
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if bits is None:
        sys.exit(f"ERROR: unsupported boundary-format {fmt!r} (expected i32/i64)")
    b = bnd.astype(np.int64)
    mask = (1 << bits) - 1
    return np.stack([(b >> (bits * k)) & mask for k in range(4)], axis=1)


def puml_free_surface_tets(puml_path):
    """For each BC-1 (free-surface) face: face centroid, owning-tet max
    edge length, and owning-tet centroid.  Returns three aligned arrays."""
    with h5py.File(puml_path, "r") as f:
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
        fmt = f.attrs.get("boundary-format", "i32")
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        bc = unpack_boundary(f["boundary"][:], fmt)

    # Same flatness sanity check as puml_h5_to_vtu.py — guards the
    # FACE_VERTS numbering this matching depends on.
    ztop = geom[:, 2].max()
    tol = 1e-6 * max(geom[:, 2].max() - geom[:, 2].min(), 1.0)

    face_cents, tet_h, tet_cents = [], [], []
    for k in range(4):
        sel = bc[:, k] == 1
        if not sel.any():
            continue
        tets = conn[sel]                                  # (n,4) vertex ids
        tri = tets[:, FACE_VERTS[k]]                      # (n,3)
        tri_xyz = geom[tri]                               # (n,3,3)
        if not np.all(np.abs(tri_xyz[:, :, 2] - ztop) < tol):
            sys.exit("ERROR: PUML face-ordering sanity check failed (BC-1 "
                     "faces not flat at z=ztop) — do not trust the matching.")
        tet_xyz = geom[tets]                              # (n,4,3)
        edges = np.stack([np.linalg.norm(tet_xyz[:, a] - tet_xyz[:, b], axis=1)
                          for a, b in TET_EDGES], axis=1)
        face_cents.append(tri_xyz.mean(axis=1))
        tet_h.append(edges.max(axis=1))
        tet_cents.append(tet_xyz.mean(axis=1))
    if not face_cents:
        sys.exit(f"ERROR: no BC-1 (free-surface) faces in {puml_path}")
    return (np.concatenate(face_cents), np.concatenate(tet_h),
            np.concatenate(tet_cents))


def match_faces(query_cents, face_cents, round_dec=3):
    """Index of each query centroid in face_cents (coordinate match,
    rounded to `round_dec` decimals = mm).  -1 where unmatched."""
    keys = {c.tobytes(): i
            for i, c in enumerate(np.round(face_cents, round_dec))}
    out = np.full(len(query_cents), -1, dtype=np.int64)
    for i, c in enumerate(np.round(query_cents, round_dec)):
        out[i] = keys.get(c.tobytes(), -1)
    return out


# --------------------------------------------------------- material sampling
def resolve_material_nc(yaml_path):
    """Pull the !ASAGI `file:` entry out of the easi yaml (custom tags make
    yaml.safe_load unusable) and resolve it next to the yaml."""
    with open(yaml_path) as fh:
        for line in fh:
            line = line.split("#", 1)[0]
            m = re.search(r"\bfile:\s*(\S+)", line)
            if m:
                return os.path.join(os.path.dirname(os.path.abspath(yaml_path)),
                                    m.group(1))
    sys.exit(f"ERROR: no 'file:' entry found in {yaml_path}")


def sample_vs(nc_path, pts):
    """Trilinear Vs = sqrt(mu/rho) at pts (n,3), clamped to the grid
    (ASAGI clamps to the boundary the same way).  Axes x,y,z are 1-D;
    `data` is a (z,y,x) compound with fields rho/mu/lambda."""
    with h5py.File(nc_path, "r") as f:
        ax = {k: f[k][:].astype(np.float64) for k in ("x", "y", "z")}
        rho = f["data"].fields("rho")[...].astype(np.float64)
        mu = f["data"].fields("mu")[...].astype(np.float64)

    idx, frac = {}, {}
    for k, q in (("x", pts[:, 0]), ("y", pts[:, 1]), ("z", pts[:, 2])):
        a = ax[k]
        qc = np.clip(q, a[0], a[-1])
        i = np.clip(np.searchsorted(a, qc, side="right") - 1, 0, len(a) - 2)
        idx[k] = i
        frac[k] = (qc - a[i]) / (a[i + 1] - a[i])

    def trilin(vol):
        ix, iy, iz = idx["x"], idx["y"], idx["z"]
        tx, ty, tz = frac["x"], frac["y"], frac["z"]
        v = np.zeros(len(pts))
        for dz in (0, 1):
            for dy in (0, 1):
                for dx in (0, 1):
                    w = ((tz if dz else 1 - tz) * (ty if dy else 1 - ty)
                         * (tx if dx else 1 - tx))
                    v += w * vol[iz + dz, iy + dy, ix + dx]
        return v

    return np.sqrt(trilin(mu) / trilin(rho))


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--surface-xdmf", default=DEF_XDMF,
                    help="SeisSol free-surface output .xdmf")
    ap.add_argument("--puml", default=DEF_PUML,
                    help=".puml.h5 volume mesh (owning-tet element size)")
    ap.add_argument("--no-puml", action="store_true",
                    help="skip the PUML; use surface-triangle size + centroid")
    ap.add_argument("--material-yaml", default=DEF_YAML,
                    help="easi material yaml whose !ASAGI file: names the .nc")
    ap.add_argument("--material-nc", default=None,
                    help="ASAGI .nc directly (overrides --material-yaml)")
    ap.add_argument("--degree", type=int, default=3,
                    help="polynomial degree p of the SeisSol build "
                         "(CONVERGENCE_ORDER - 1; default 3 for the O4 "
                         "QuakeWorx app) — sets h_ip = h/(p+1)")
    ap.add_argument("--out", default=DEF_OUT, help="output .vtu path")
    args = ap.parse_args()

    pts, tris = read_surface_xdmf(args.surface_xdmf)
    tri_xyz = pts[tris]                                   # (nT,3,3)
    tri_cents = tri_xyz.mean(axis=1)
    e = [np.linalg.norm(tri_xyz[:, a] - tri_xyz[:, b], axis=1)
         for a, b in ((0, 1), (1, 2), (2, 0))]
    h_tri = np.max(e, axis=0)
    print(f"surface: {len(tris)} triangles, {len(pts)} vertices, "
          f"z in [{pts[:,2].min():.1f}, {pts[:,2].max():.1f}] m")
    print(f"h_tri_max_edge [m]: min {h_tri.min():.1f}  median "
          f"{np.median(h_tri):.1f}  max {h_tri.max():.1f}")

    h_tet = h_tri.copy()
    sample_pts = tri_cents.copy()
    if not args.no_puml:
        face_cents, tet_h_all, tet_cents_all = puml_free_surface_tets(args.puml)
        m = match_faces(tri_cents, face_cents)
        n_miss = int(np.sum(m < 0))
        print(f"PUML matching: {len(m) - n_miss}/{len(m)} surface triangles "
              f"matched to BC-1 faces ({len(face_cents)} in mesh)")
        if n_miss > 0.01 * len(m):
            sys.exit(f"ERROR: {n_miss} unmatched triangles (>1%) — surface "
                     "output and PUML mesh do not correspond; re-run with "
                     "--no-puml to use surface-triangle sizes instead.")
        ok = m >= 0
        h_tet[ok] = tet_h_all[m[ok]]
        sample_pts[ok] = tet_cents_all[m[ok]]
        print(f"h_tet_max_edge [m]: min {h_tet.min():.1f}  median "
              f"{np.median(h_tet):.1f}  max {h_tet.max():.1f}")

    nc_path = args.material_nc or resolve_material_nc(args.material_yaml)
    print(f"material: {nc_path}")
    vs = sample_vs(nc_path, sample_pts)
    period = h_tet / vs
    n_pts_edge = args.degree + 1
    h_ip = h_tet / n_pts_edge
    period_ip = period / n_pts_edge
    print(f"vs [m/s]: min {vs.min():.1f}  median {np.median(vs):.1f}  "
          f"max {vs.max():.1f}")
    print(f"period h/vs [s]: min {period.min():.4f}  median "
          f"{np.median(period):.4f}  max {period.max():.4f}")
    print(f"period_ip (h/{n_pts_edge})/vs [s] (p={args.degree}): "
          f"min {period_ip.min():.4f}  median {np.median(period_ip):.4f}  "
          f"max {period_ip.max():.4f}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_vtu(args.out, pts, tris, VTK_TRIANGLE, {
        "h_tri_max_edge": h_tri.astype(np.float32),
        "h_tet_max_edge": h_tet.astype(np.float32),
        "h_ip": h_ip.astype(np.float32),
        "vs": vs.astype(np.float32),
        "period": period.astype(np.float32),
        "period_ip": period_ip.astype(np.float32),
    })
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
