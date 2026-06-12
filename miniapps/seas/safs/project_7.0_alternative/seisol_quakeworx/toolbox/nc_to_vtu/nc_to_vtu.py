#!/usr/bin/env python3
"""Convert an ASAGI material netCDF (rectilinear grid) to a ParaView VTU.

Input layout (as written by convert_cvm_to_asagi.py, netCDF4 = HDF5,
read here with h5py — no netCDF4/vtk/meshio dependency):
    x (nx,), y (ny,), z (nz,)   axis coordinates, metres
    data (nz, ny, nx)           compound (rho, mu, lambda), SI units
(a file with separate plain datasets rho/mu/lambda is also accepted)

Output: one VTU of hexahedral cells on the rectilinear grid with POINT data
    rho_kg_m3, mu_GPa, lambda_GPa            (raw fields)
    Vs_m_s, Vp_m_s, poisson_ratio            (derived; --no-derived to skip)
Values live on grid nodes, so point data + voxel cells lets ParaView
interpolate exactly like easi/ASAGI `interpolation: linear`.

The full SAFS grid (264 x 195 x 181) gives a ~660 MB VTU; use
--stride N to subsample every Nth node per axis (e.g. --stride 2
≈ 1/8 the size) and/or --no-derived for a lighter file.

Usage:
    python3 nc_to_vtu.py <material.nc> [--out-dir DIR] [--stride N]
        [--no-derived]
"""
import argparse
import os
import struct
import sys

import h5py
import numpy as np

VTK_HEXAHEDRON = 12  # universally supported (meshio chokes on VOXEL=11)

_VTK_TYPE = {np.dtype("<f8"): "Float64", np.dtype("<f4"): "Float32",
             np.dtype("<i8"): "Int64", np.dtype("<i4"): "Int32",
             np.dtype("u1"): "UInt8"}


def write_vtu(path, points, cells, cell_type, point_data):
    """XML VTU, appended raw binary, little endian, UInt64 headers."""
    points = np.ascontiguousarray(points)
    cells = np.ascontiguousarray(cells)
    n_pts, n_cells = len(points), len(cells)
    offsets = (np.arange(1, n_cells + 1, dtype="<i8") * cells.shape[1])
    if cells.dtype == np.dtype("<i4"):
        offsets = offsets.astype("<i4")
    types = np.full(n_cells, cell_type, dtype="u1")

    # connectivity MUST be a flat single-component array — declaring it
    # with NumberOfComponents=k makes vtkXMLUnstructuredGridReader fail
    # ("could not be created with one component"); conn_flat is used for
    # BOTH the payload and the header declaration below
    conn_flat = cells.reshape(-1)
    blocks = [("Points", points), ("connectivity", conn_flat),
              ("offsets", offsets), ("types", types)]
    pd = list(point_data.items())
    raw_all = blocks + pd
    offs, pos = [], 0
    for _, a in raw_all:
        offs.append(pos)
        pos += 8 + a.nbytes

    def da(name, arr, off):
        vtype = _VTK_TYPE[arr.dtype]
        ncomp = 1 if arr.ndim == 1 else arr.shape[1]
        comp = f' NumberOfComponents="{ncomp}"' if ncomp > 1 else ""
        return (f'        <DataArray type="{vtype}" Name="{name}"{comp}'
                f' format="appended" offset="{off}"/>\n')

    xml = ['<?xml version="1.0"?>\n',
           '<VTKFile type="UnstructuredGrid" version="1.0" '
           'byte_order="LittleEndian" header_type="UInt64">\n',
           '  <UnstructuredGrid>\n',
           f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n',
           '      <Points>\n', da("Points", points, offs[0]),
           '      </Points>\n', '      <Cells>\n',
           da("connectivity", conn_flat, offs[1]),
           da("offsets", offsets, offs[2]),
           da("types", types, offs[3]), '      </Cells>\n']
    if pd:
        xml.append('      <PointData>\n')
        for (name, arr), off in zip(pd, offs[4:]):
            xml.append(da(name, arr, off))
        xml.append('      </PointData>\n')
    xml += ['    </Piece>\n', '  </UnstructuredGrid>\n',
            '  <AppendedData encoding="raw">\n_']
    with open(path, "wb") as fh:
        fh.write("".join(xml).encode())
        for _, a in raw_all:
            fh.write(struct.pack("<Q", a.nbytes))
            fh.write(a.tobytes())
        fh.write(b"\n  </AppendedData>\n</VTKFile>\n")


def load_nc(path, stride):
    """Return xs, ys, zs and {rho, mu, lambda} (nz, ny, nx) float64."""
    with h5py.File(path, "r") as f:
        xs = f["x"][::stride]
        ys = f["y"][::stride]
        zs = f["z"][::stride]
        fields = {}
        if "data" in f and f["data"].dtype.names:  # compound ASAGI var
            data = f["data"][::stride, ::stride, ::stride]
            for name in ("rho", "mu", "lambda"):
                if name not in data.dtype.names:
                    sys.exit(f"ERROR: compound var lacks field {name!r} "
                             f"(has {data.dtype.names})")
                fields[name] = data[name].astype(np.float64)
        else:                                      # separate plain datasets
            for name in ("rho", "mu", "lambda"):
                if name not in f:
                    sys.exit(f"ERROR: dataset {name!r} not found in {path}")
                fields[name] = f[name][::stride, ::stride, ::stride] \
                    .astype(np.float64)
        attrs = {k: (v.decode() if isinstance(v, bytes) else v)
                 for k, v in f.attrs.items()}
    nz, ny, nx = fields["rho"].shape
    if (len(zs), len(ys), len(xs)) != (nz, ny, nx):
        sys.exit(f"ERROR: axis lengths ({len(zs)},{len(ys)},{len(xs)}) do "
                 f"not match data shape ({nz},{ny},{nx}); expected "
                 "data[z, y, x] layout")
    return xs, ys, zs, fields, attrs


def hex_grid(xs, ys, zs):
    """Points (n,3) f32 in x-fastest order + hexahedron connectivity i32."""
    nx, ny, nz = len(xs), len(ys), len(zs)
    pts = np.empty((nz, ny, nx, 3), dtype="<f4")
    pts[..., 0] = xs[None, None, :]
    pts[..., 1] = ys[None, :, None]
    pts[..., 2] = zs[:, None, None]
    pts = pts.reshape(-1, 3)

    iz, iy, ix = np.meshgrid(np.arange(nz - 1), np.arange(ny - 1),
                             np.arange(nx - 1), indexing="ij")
    base = ((iz * ny + iy) * nx + ix).ravel()
    # VTK_HEXAHEDRON vertex order: bottom face counter-clockwise
    # (-x-y, +x-y, +x+y, -x+y), then the same on the top face
    corner = np.array([0, 1, nx + 1, nx,
                       nx * ny, nx * ny + 1, nx * ny + nx + 1, nx * ny + nx])
    conn = (base[:, None] + corner[None, :]).astype("<i4")
    return pts, conn


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("nc", help="ASAGI material netCDF (netCDF4/HDF5)")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--stride", type=int, default=1,
                    help="subsample every Nth grid node per axis (default 1)")
    ap.add_argument("--no-derived", action="store_true",
                    help="skip Vs/Vp/Poisson-ratio derived fields")
    args = ap.parse_args()
    if args.stride < 1:
        sys.exit("ERROR: --stride must be >= 1")

    os.makedirs(args.out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.nc))[0]
    if args.stride > 1:
        stem += f"_s{args.stride}"

    xs, ys, zs, fld, attrs = load_nc(args.nc, args.stride)
    nx, ny, nz = len(xs), len(ys), len(zs)
    print(f"grid: {nx} x {ny} x {nz} = {nx*ny*nz} nodes "
          f"(stride {args.stride}); x [{xs[0]:.0f}, {xs[-1]:.0f}] "
          f"y [{ys[0]:.0f}, {ys[-1]:.0f}] z [{zs[0]:.0f}, {zs[-1]:.0f}] m")
    if "units" in attrs:
        print("units:", attrs["units"])

    rho, mu, lam = fld["rho"], fld["mu"], fld["lambda"]
    point_data = {
        "rho_kg_m3": rho.reshape(-1).astype("<f4"),
        "mu_GPa": (mu / 1e9).reshape(-1).astype("<f4"),
        "lambda_GPa": (lam / 1e9).reshape(-1).astype("<f4"),
    }
    if not args.no_derived:
        with np.errstate(invalid="ignore", divide="ignore"):
            point_data["Vs_m_s"] = np.sqrt(mu / rho) \
                .reshape(-1).astype("<f4")
            point_data["Vp_m_s"] = np.sqrt((lam + 2.0 * mu) / rho) \
                .reshape(-1).astype("<f4")
            point_data["poisson_ratio"] = (lam / (2.0 * (lam + mu))) \
                .reshape(-1).astype("<f4")
    for name, arr in point_data.items():
        bad = int((~np.isfinite(arr)).sum())
        tag = f"  [{bad} non-finite]" if bad else ""
        print(f"  {name:14s} min {np.nanmin(arr):12.4g}  "
              f"max {np.nanmax(arr):12.4g}{tag}")

    pts, conn = hex_grid(xs, ys, zs)
    out = os.path.join(args.out_dir, f"{stem}.vtu")
    write_vtu(out, pts, conn, VTK_HEXAHEDRON, point_data)
    print(f"wrote {out}  ({len(pts)} points, {len(conn)} hexahedra, "
          f"{os.path.getsize(out)/1e6:.0f} MB)")


if __name__ == "__main__":
    main()
