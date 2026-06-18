#!/usr/bin/env python3
"""generate_velocity_nc_from_raw.py — raw CVM ASCII slices -> ASAGI NetCDF,
in one file.

Single-file combination of the two-stage pipeline that produced
safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc:

  Stage 1 (was velocity/code/build_velocity_cvmh.py + raw_readers.py + crs.py)
    - read every velocity_raw_*.bp ASCII slice (one horizontal slice per
      depth; lon/lat regular grid; columns lon,lat,vp,vs,density)
    - reproject samples EPSG:4326 -> UTM 11N (EPSG:32611, pyproj)
    - build the inscribed axis-aligned UTM grid at --grid-dx spacing
      (rounded inward so every node is inside the source hull)
    - per slice, per field: scipy LinearNDInterpolator (linear,
      Delaunay-based) onto the UTM grid; NaN cells fail loud
    - flip depth -> elevation (z negative down) and clone the surface
      slice +--extend-z-top metres up (covers mesh PAD_TOP)

  Stage 2 (was safs_seisol_v2_0_0_RSSRW/convert_cvm_to_asagi.py)
    - moduli AT THE SOURCE NODES:  mu = rho Vs^2,
      lambda = rho (Vp^2 - 2 Vs^2)   (fails if lambda <= 0 anywhere)
    - linear z-resample onto the uniform output axis (ASAGI needs
      equidistant axes; conversion-before-resampling keeps the stored
      grid the exact piecewise-linear interpolant of node moduli)
    - write COARDS NetCDF4: x/y/z float64 + compound data(z,y,x)
      {rho; mu; lambda}, the layout SeisSol's easi !ASAGI reader expects
    - fixed-seed round-trip self-check (trilinear NetCDF samples vs an
      independent evaluation of the source stack); non-zero exit on FAIL

Production parameters of the shipped safs_material_cvm.nc are the
defaults below: --grid-dx 1500, --extend-z-top 100, --z-min -45000,
--z-max 100 (axis tops out at 0), --dz 250, --dtype float32, no
vs-floor.

Usage (from the folder containing toolbox/):
    python3 toolbox/generate_velocity_nc_from_raw/generate_velocity_nc_from_raw.py
        [--raw-dir raw_data/multiscale_statewise_cvm]
        [--out toolbox/generate_velocity_nc_from_raw/safs_material_cvm.nc]
        [--write-sidecar PATH.h5] [--mesh-puml PATH.puml.h5]
        [--compare-to existing.nc]

Requires: numpy, scipy, pyproj, netCDF4 (and h5py only for
--write-sidecar / --mesh-puml).  All defaults are resolved relative to
this script, so the folder tree stays relocatable.
"""
from __future__ import annotations

import argparse
import datetime
import os
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent          # the quakeworx case-folder root
DEFAULT_RAW_DIR = _ROOT / "raw_data" / "multiscale_statewise_cvm"
DEFAULT_OUT = _HERE / "safs_material_cvm.nc"

# ---- self-check constants (identical to convert_cvm_to_asagi.py) ----
# SAFS SeisSol mesh bounding box (safs_mesh.puml.h5 geometry, UTM 11N, m)
MESH_BBOX = ((314840.1, 672329.5), (3642278.0, 3888985.0), (-41607.61, 0.0))
SELF_CHECK_SEED = 20260609
SELF_CHECK_NPTS = 1000
SELF_CHECK_MEDIAN_TOL = 1.0e-6
SELF_CHECK_MAX_TOL = 5.0e-3

_REQUIRED_HEADER_KEYS = ("Depth(m)", "Lon_pts", "Lat_pts", "Total_pts",
                         "Lat1", "Lon1", "Spacing(degree)")


# ======================================================================
# Stage 1a — raw .bp slice reader  (from velocity/code/raw_readers.py)
# ======================================================================
def _parse_header(path):
    header, column_line = {}, None
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if not line.startswith("#"):
                break
            stripped = line[1:].strip()
            if stripped.lower().startswith("lon") and "," in stripped:
                column_line = stripped
                continue
            if ":" in stripped:
                key, _, value = stripped.partition(":")
                header[key.strip()] = value.strip()
    if column_line is None:
        raise ValueError(f"{path}: missing column-spec header line")
    header["_column_line"] = column_line
    missing = [k for k in _REQUIRED_HEADER_KEYS if k not in header]
    if missing:
        raise ValueError(f"{path}: missing header keys {missing}")
    return header


def read_raw_slices(raw_dir, verbose=True):
    """All velocity_raw_*.bp in raw_dir -> (lon, lat, depths_m, fields)
    with fields {vp, vs, density} -> (n_lon, n_lat, n_z) float64,
    depth ascending."""
    paths = sorted(Path(raw_dir).glob("velocity_raw_*.bp"))
    if not paths:
        raise FileNotFoundError(f"no velocity_raw_*.bp in {raw_dir}")
    if verbose:
        print(f"stage 1: {len(paths)} raw slices in {raw_dir}")

    slices, field_names = [], None
    for p in paths:
        h = _parse_header(p)
        depth_m = float(h["Depth(m)"])
        n_lon, n_lat = int(h["Lon_pts"]), int(h["Lat_pts"])
        n_total = int(h["Total_pts"])
        if n_lon * n_lat != n_total:
            raise ValueError(f"{p}: Lon_pts*Lat_pts != Total_pts")
        cols = []
        for tok in h["_column_line"].split(","):
            tok = tok.strip()
            paren = tok.find("(")
            cols.append((tok[:paren] if paren >= 0 else tok).strip().lower())
        fld_cols = [c for c in cols if c not in ("lon", "lat")]
        if field_names is None:
            field_names = fld_cols
        elif fld_cols != field_names:
            raise ValueError(f"{p}: field columns {fld_cols} != {field_names}")

        body = np.loadtxt(str(p), delimiter=",", comments="#")
        if body.shape != (n_total, len(cols)):
            raise ValueError(f"{p}: body shape {body.shape} != "
                             f"({n_total}, {len(cols)})")
        # corner pinning: first row = SW (Lon1, Lat1); last row = NE
        # derived from spacing (the archive's Lon2/Lat2 header values
        # are buggy by up to ~3 sample widths -> ignored)
        lon1, lat1 = float(h["Lon1"]), float(h["Lat1"])
        sp = float(h["Spacing(degree)"])
        if (abs(body[0, 0] - lon1) > 1e-3 or abs(body[0, 1] - lat1) > 1e-3):
            raise ValueError(f"{p}: first row is not the (Lon1, Lat1) corner")
        if (abs(body[-1, 0] - (lon1 + (n_lon - 1) * sp)) > 1e-3 or
                abs(body[-1, 1] - (lat1 + (n_lat - 1) * sp)) > 1e-3):
            raise ValueError(f"{p}: last row is not the derived NE corner")

        lon_u = np.unique(np.round(body[:, 0], 6))
        lat_u = np.unique(np.round(body[:, 1], 6))
        if lon_u.size != n_lon or lat_u.size != n_lat:
            raise ValueError(f"{p}: duplicate/noisy lon-lat values")
        li = np.searchsorted(lon_u, np.round(body[:, 0], 6))
        ti = np.searchsorted(lat_u, np.round(body[:, 1], 6))
        flds = {}
        for ci, name in enumerate(cols):
            if name in ("lon", "lat"):
                continue
            arr = np.full((n_lon, n_lat), np.nan)
            arr[li, ti] = body[:, ci]
            if np.isnan(arr).any():
                raise ValueError(f"{p}: NaN holes after reshape ('{name}')")
            flds[name] = arr
        slices.append((depth_m, lon_u, lat_u, flds))

    ref_lon, ref_lat = slices[0][1], slices[0][2]
    for d, lo, la, _ in slices[1:]:
        if not (np.array_equal(lo, ref_lon) and np.array_equal(la, ref_lat)):
            raise ValueError(f"slice at depth {d} m on a different lon/lat "
                             "grid")
    slices.sort(key=lambda s: s[0])
    depths = np.array([s[0] for s in slices])
    if np.any(np.diff(depths) <= 0):
        raise ValueError("duplicate depths in slice set")
    fields = {n: np.empty((ref_lon.size, ref_lat.size, depths.size))
              for n in field_names}
    for k, (_, _, _, fl) in enumerate(slices):
        for n in field_names:
            fields[n][:, :, k] = fl[n]
    if sorted(field_names) != ["density", "vp", "vs"]:
        raise ValueError(f"expected fields vp/vs/density; got {field_names}")
    return ref_lon, ref_lat, depths, fields


# ======================================================================
# Stage 1b — UTM grid + resample  (from velocity/code/build_velocity_cvmh.py)
# ======================================================================
def geographic_to_utm11n(lon, lat):
    try:
        from pyproj import Transformer
    except ImportError:
        sys.exit("ERROR: pyproj is required (pip install pyproj)")
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
    x, y = tr.transform(np.asarray(lon).ravel(), np.asarray(lat).ravel())
    return (np.asarray(x).reshape(np.shape(lon)),
            np.asarray(y).reshape(np.shape(lat)))


def build_utm_grid(lon, lat, grid_dx):
    """Inscribed axis-aligned UTM bbox of the projected source quad,
    rounded INWARD to grid_dx multiples (every node inside the hull)."""
    LON, LAT = np.meshgrid(lon, lat, indexing="ij")
    X, Y = geographic_to_utm11n(LON, LAT)
    xmin = np.ceil(float(X[0, :].max()) / grid_dx) * grid_dx
    xmax = np.floor(float(X[-1, :].min()) / grid_dx) * grid_dx
    ymin = np.ceil(float(Y[:, 0].max()) / grid_dx) * grid_dx
    ymax = np.floor(float(Y[:, -1].min()) / grid_dx) * grid_dx
    if xmin > xmax or ymin > ymax:
        raise ValueError(f"grid_dx={grid_dx} too large for source extent")
    x_axis = np.arange(xmin, xmax + 0.5 * grid_dx, grid_dx)
    y_axis = np.arange(ymin, ymax + 0.5 * grid_dx, grid_dx)
    return x_axis, y_axis, X, Y


def resample_to_utm(lon, lat, depths, fields, x_axis, y_axis, X_pts, Y_pts,
                    extend_z_top_m, verbose=True):
    """LinearNDInterpolator per slice/field onto the UTM grid; flip depth
    -> elevation; clone the surface slice extend_z_top_m up.
    Returns (z_axis ascending elevation, {Vp,Vs,density}->(nx,ny,nz))."""
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay

    pts_xy = np.stack([X_pts.ravel(), Y_pts.ravel()], axis=1)
    tri = Delaunay(pts_xy)        # one triangulation, reused per slice
    XX, YY = np.meshgrid(x_axis, y_axis, indexing="ij")
    target = np.stack([XX.ravel(), YY.ravel()], axis=1)

    nx, ny, nz = x_axis.size, y_axis.size, depths.size
    out = {}
    for name in fields:
        arr = np.empty((nx, ny, nz))
        for k in range(nz):
            interp = LinearNDInterpolator(tri, fields[name][:, :, k].ravel(),
                                          fill_value=np.nan)
            arr[:, :, k] = interp(target).reshape(nx, ny)
        n_nan = int(np.isnan(arr).sum())
        if n_nan:
            raise ValueError(
                f"{n_nan} NaN cells in '{name}' after resample — the UTM "
                "grid extends beyond the source hull")
        if verbose:
            print(f"  resampled {name}: ({nx}, {ny}, {nz}), "
                  f"range [{arr.min():.6g}, {arr.max():.6g}]")
        out[name] = arr

    z_axis = np.sort(-depths)             # elevation ascending
    flipped = {n: a[:, :, ::-1].copy() for n, a in out.items()}
    if extend_z_top_m > 0.0:
        z_axis = np.concatenate([z_axis, [z_axis[-1] + extend_z_top_m]])
        for n in flipped:
            flipped[n] = np.concatenate(
                [flipped[n], flipped[n][:, :, -1:]], axis=2)
    canonical = {"vp": "Vp", "vs": "Vs", "density": "density"}
    return z_axis, {canonical[n]: a for n, a in flipped.items()}


# ======================================================================
# Stage 2 — moduli + z-resample + NetCDF  (from convert_cvm_to_asagi.py)
# ======================================================================
def velocities_to_moduli(vp, vs, rho):
    mu = rho * vs**2
    lam = rho * (vp**2 - 2.0 * vs**2)
    if not np.all(mu > 0):
        raise ValueError(f"mu must be positive; min {mu.min():.6g}")
    if not np.all(lam > 0):
        raise ValueError(
            f"lambda <= 0 at {(lam <= 0).sum()} node(s); requires "
            f"Vp > sqrt(2) Vs (min Vp/Vs = {(vp / vs).min():.4f})")
    return rho, mu, lam


def resample_z(field, z_src, z_out):
    if z_out[0] < z_src[0] or z_out[-1] > z_src[-1]:
        raise ValueError("z_out outside z_src (no extrapolation)")
    idx = np.clip(np.searchsorted(z_src, z_out, side="right") - 1,
                  0, z_src.size - 2)
    w = np.clip((z_out - z_src[idx]) / (z_src[idx + 1] - z_src[idx]), 0., 1.)
    return field[:, :, idx] * (1.0 - w) + field[:, :, idx + 1] * w


def write_asagi_netcdf(path, x, y, z, rho, mu, lam, dtype, attrs):
    try:
        import netCDF4
    except ImportError:
        sys.exit("ERROR: netCDF4 is required (pip install netCDF4)")
    np_dtype = np.dtype(dtype)
    material_t = np.dtype([("rho", np_dtype), ("mu", np_dtype),
                           ("lambda", np_dtype)])
    with netCDF4.Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.createDimension("x", x.size)
        ds.createDimension("y", y.size)
        ds.createDimension("z", z.size)
        for name, ax in (("x", x), ("y", y), ("z", z)):
            v = ds.createVariable(name, "f8", (name,))
            v[:] = ax
        mtype = ds.createCompoundType(material_t, "material")
        data = ds.createVariable("data", mtype, ("z", "y", "x"))
        buf = np.empty((z.size, y.size, x.size), dtype=material_t)
        buf["rho"] = np.transpose(rho, (2, 1, 0)).astype(np_dtype)
        buf["mu"] = np.transpose(mu, (2, 1, 0)).astype(np_dtype)
        buf["lambda"] = np.transpose(lam, (2, 1, 0)).astype(np_dtype)
        data[:] = buf
        for key, val in attrs.items():
            ds.setncattr(key, val)


def read_asagi_netcdf(path):
    import netCDF4
    with netCDF4.Dataset(str(path), "r") as ds:
        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        y = np.asarray(ds.variables["y"][:], dtype=np.float64)
        z = np.asarray(ds.variables["z"][:], dtype=np.float64)
        raw = ds.variables["data"][:]
        fields = {n: np.transpose(np.asarray(raw[n], dtype=np.float64),
                                  (2, 1, 0))
                  for n in ("rho", "mu", "lambda")}
    return x, y, z, fields


def trilinear_sample(xg, yg, zg, field, pts):
    def locate(g, c):
        i = np.clip(np.searchsorted(g, c, side="right") - 1, 0, g.size - 2)
        return i, np.clip((c - g[i]) / (g[i + 1] - g[i]), 0.0, 1.0)
    ix, wx = locate(xg, pts[:, 0])
    iy, wy = locate(yg, pts[:, 1])
    iz, wz = locate(zg, pts[:, 2])
    out = np.zeros(len(pts))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                wgt = ((wx if dx else 1 - wx) * (wy if dy else 1 - wy)
                       * (wz if dz else 1 - wz))
                out += wgt * field[ix + dx, iy + dy, iz + dz]
    return out


def self_check(nc_path, x, y, z, fields, smoothed_band_zmin=None):
    """Trilinear NetCDF samples vs independent source-stack evaluation at
    SELF_CHECK_NPTS fixed-seed points in the SAFS mesh bbox.

    A coarse --dz that misses shallow source z-levels deliberately
    smooths the model there (documented production trade-off at
    dz = 250: the 4 levels at -100..-400 m, worst ~28% on mu in the
    weathered zone — see safs_material_cvm.yaml).  When
    ``smoothed_band_zmin`` is given, max-tolerance excursions whose
    sample points ALL lie above that elevation are reported but do not
    fail the check; the median tolerance (which catches transposition /
    ordering / resampling bugs at every depth) always applies.
    """
    rng = np.random.default_rng(SELF_CHECK_SEED)
    pts = np.column_stack([rng.uniform(lo, hi, SELF_CHECK_NPTS)
                           for (lo, hi) in MESH_BBOX])
    xg, yg, zg, nc_fields = read_asagi_netcdf(nc_path)
    rho_n, mu_n, lam_n = velocities_to_moduli(
        fields["Vp"], fields["Vs"], fields["density"])
    ref = {"rho": rho_n, "mu": mu_n, "lambda": lam_n}
    rel = np.concatenate([
        np.abs(trilinear_sample(xg, yg, zg, nc_fields[n], pts)
               - trilinear_sample(x, y, z, ref[n], pts))
        / np.abs(trilinear_sample(x, y, z, ref[n], pts))
        for n in ("rho", "mu", "lambda")])
    med, mx = float(np.median(rel)), float(rel.max())
    note = ""
    max_ok = mx <= SELF_CHECK_MAX_TOL
    if not max_ok and smoothed_band_zmin is not None:
        bad_z = np.tile(pts[:, 2], 3)[rel > SELF_CHECK_MAX_TOL]
        if bad_z.size and bad_z.min() >= smoothed_band_zmin:
            max_ok = True
            note = (f" [max-err excursions confined to the dz-smoothed "
                    f"shallow band z >= {smoothed_band_zmin:.0f} m — "
                    f"expected at this --dz]")
    return med, mx, (med <= SELF_CHECK_MEDIAN_TOL and max_ok), note


# ======================================================================
# Optional extras
# ======================================================================
def write_sidecar_h5(path, x, y, z, fields, source_names):
    """Intermediate schema-v1 sidecar (same layout the two-stage pipeline
    stored under velocity/results/<version>/velocity_safs.h5)."""
    import h5py
    units = {"Vp": "m/s", "Vs": "m/s", "density": "kg/m^3"}
    with h5py.File(path, "w") as h5:
        h5.attrs.update({
            "schema_version": "data_projection_v1",
            "crs": "EPSG:32611", "units": "m", "z_positive": "elevation",
            "created_at": datetime.datetime.now(datetime.timezone.utc)
                          .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "source": ",".join(source_names), "source_crs": "EPSG:4326",
            "mesh_tag": "",
        })
        g = h5.create_group("grid")
        for n, ax in (("x", x), ("y", y), ("z", z)):
            g.create_dataset(n, data=ax)
        f = h5.create_group("fields")
        for n, arr in fields.items():
            d = f.create_dataset(n, data=arr)
            d.attrs["units"] = units[n]
            d.attrs["min_value"] = float(arr.min()) - 1.0
            d.attrs["max_value"] = float(arr.max()) + 1.0


def check_mesh_containment(puml_path, x, y, z):
    """G-3 guard: the data grid must contain the PUML mesh bbox."""
    import h5py
    with h5py.File(puml_path, "r") as f:
        geom = f["geometry"]
        lo = geom[:].min(axis=0)
        hi = geom[:].max(axis=0)
    for k, (name, ax) in enumerate((("x", x), ("y", y), ("z", z))):
        if lo[k] < ax[0] or hi[k] > ax[-1]:
            raise ValueError(
                f"mesh bbox {name} [{lo[k]:.1f}, {hi[k]:.1f}] not contained "
                f"in data grid [{ax[0]:.1f}, {ax[-1]:.1f}]")
    print(f"mesh guard  : PUML bbox contained in data grid "
          f"({os.path.basename(str(puml_path))})")


def compare_to_reference(out_path, ref_path):
    """Element-wise comparison of two ASAGI NetCDFs (axes + fields)."""
    x1, y1, z1, f1 = read_asagi_netcdf(out_path)
    x2, y2, z2, f2 = read_asagi_netcdf(ref_path)
    ok = True
    for n, (a, b) in (("x", (x1, x2)), ("y", (y1, y2)), ("z", (z1, z2))):
        same = a.shape == b.shape and np.array_equal(a, b)
        ok &= same
        print(f"compare     : axis {n}: "
              f"{'identical' if same else 'DIFFERENT'} ({a.size} vs {b.size})")
    if ok:
        for n in ("rho", "mu", "lambda"):
            d = np.abs(f1[n] - f2[n])
            scale = np.abs(f2[n]).max()
            print(f"compare     : {n:6s} max abs diff {d.max():.6g} "
                  f"(rel {d.max() / scale:.3e})"
                  f"{'  [bit-identical]' if d.max() == 0 else ''}")
            ok &= bool(d.max() <= 1e-6 * scale)
    print(f"compare     : {'MATCH' if ok else 'MISMATCH'} vs {ref_path}")
    return ok


# ======================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR,
                    help="directory of velocity_raw_*.bp slices "
                    "(default: raw_data/multiscale_statewise_cvm)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="output ASAGI NetCDF (default: next to this script)")
    ap.add_argument("--grid-dx", type=float, default=1500.0,
                    help="UTM grid spacing, m (default %(default)s)")
    ap.add_argument("--extend-z-top", type=float, default=100.0,
                    help="clone the surface slice this many m up "
                    "(default %(default)s = mesh PAD_TOP)")
    ap.add_argument("--z-min", type=float, default=-45000.0)
    ap.add_argument("--z-max", type=float, default=100.0)
    ap.add_argument("--dz", type=float, default=250.0,
                    help="uniform output z spacing, m (default %(default)s "
                    "= production value; 50 captures every source level)")
    ap.add_argument("--dtype", choices=("float32", "float64"),
                    default="float32")
    ap.add_argument("--write-sidecar", type=Path, default=None,
                    help="also write the intermediate schema-v1 sidecar h5")
    ap.add_argument("--mesh-puml", type=Path, default=None,
                    help="optional PUML mesh for the bbox containment guard")
    ap.add_argument("--compare-to", type=Path, default=None,
                    help="compare the produced NetCDF against an existing "
                    "reference .nc (e.g. the shipped safs_material_cvm.nc)")
    args = ap.parse_args(argv)

    # ---- stage 1 ----
    lon, lat, depths, raw_fields = read_raw_slices(args.raw_dir)
    src_names = sorted(p.name for p in Path(args.raw_dir)
                       .glob("velocity_raw_*.bp"))
    x_axis, y_axis, X_pts, Y_pts = build_utm_grid(lon, lat, args.grid_dx)
    print(f"utm grid    : x [{x_axis[0]:.0f}, {x_axis[-1]:.0f}] "
          f"({x_axis.size}) | y [{y_axis[0]:.0f}, {y_axis[-1]:.0f}] "
          f"({y_axis.size}) @ dx {args.grid_dx:g} m")
    z_axis, fields = resample_to_utm(lon, lat, depths, raw_fields,
                                     x_axis, y_axis, X_pts, Y_pts,
                                     args.extend_z_top)
    print(f"sidecar grid: {x_axis.size} x {y_axis.size} x {z_axis.size} "
          f"(z [{z_axis[0]:.0f}, {z_axis[-1]:.0f}], non-uniform)")
    if args.mesh_puml is not None:
        check_mesh_containment(args.mesh_puml, x_axis, y_axis, z_axis)
    if args.write_sidecar is not None:
        write_sidecar_h5(args.write_sidecar, x_axis, y_axis, z_axis,
                         fields, src_names)
        print(f"sidecar     : wrote {args.write_sidecar}")

    # ---- stage 2 ----
    if args.z_min < z_axis[0] or args.z_max > z_axis[-1]:
        sys.exit(f"ERROR: --z-min/--z-max outside source z range "
                 f"[{z_axis[0]}, {z_axis[-1]}]")
    z_out = np.arange(args.z_min, args.z_max + 0.5 * args.dz, args.dz)
    if z_out[-1] > args.z_max:
        z_out = z_out[:-1]
    in_range = z_axis[(z_axis >= args.z_min) & (z_axis <= z_out[-1])]
    k = np.rint((in_range - args.z_min) / args.dz)
    missed = in_range[~np.isclose(args.z_min + k * args.dz, in_range,
                                  atol=1e-6 * args.dz)]
    if missed.size:
        print(f"WARNING     : {missed.size} source z-level(s) off the "
              f"output grid (smoothed): "
              f"{', '.join(f'{v:.0f}' for v in missed)}")

    rho_n, mu_n, lam_n = velocities_to_moduli(
        fields["Vp"], fields["Vs"], fields["density"])
    rho = resample_z(rho_n, z_axis, z_out)
    mu = resample_z(mu_n, z_axis, z_out)
    lam = resample_z(lam_n, z_axis, z_out)

    # record the raw dir relative to the case-folder root when inside it
    # (keeps uploaded artifacts free of machine-specific absolute paths)
    raw_rel = os.path.relpath(os.path.abspath(args.raw_dir), _ROOT)
    attrs = {
        "source_raw_dir": (str(args.raw_dir) if raw_rel.startswith("..")
                           else raw_rel),
        "source_slices": ",".join(src_names),
        "source_crs": "EPSG:4326",
        "crs": "EPSG:32611",
        "z_positive": "elevation",
        "units": "SI: rho kg/m^3, mu/lambda Pa, axes m",
        "grid_dx": float(args.grid_dx),
        "extend_z_top": float(args.extend_z_top),
        "dz": float(args.dz),
        "vs_floor": "none",
        "converter": "generate_velocity_nc_from_raw.py (single-file "
                     "combination of build_velocity_cvmh.py + "
                     "convert_cvm_to_asagi.py)",
        "converted_at": datetime.datetime.now(datetime.timezone.utc)
                        .strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    write_asagi_netcdf(args.out, x_axis, y_axis, z_out, rho, mu, lam,
                       args.dtype, attrs)
    print(f"written     : {args.out} "
          f"({os.path.getsize(args.out) / 1e6:.1f} MB, {args.dtype}); "
          f"output z [{z_out[0]:.0f}, {z_out[-1]:.0f}] @ dz {args.dz:g} "
          f"({z_out.size} levels)")
    vs_eff, vp_eff = np.sqrt(mu / rho), np.sqrt((lam + 2 * mu) / rho)
    print(f"  rho [{rho.min():.6g}, {rho.max():.6g}] kg/m^3 | "
          f"mu [{mu.min():.6g}, {mu.max():.6g}] Pa")
    print(f"  Vs  [{vs_eff.min():.6g}, {vs_eff.max():.6g}] m/s | "
          f"Vp [{vp_eff.min():.6g}, {vp_eff.max():.6g}] m/s")

    band_zmin = float(missed.min() - args.dz) if missed.size else None
    med, mx, passed, note = self_check(args.out, x_axis, y_axis, z_axis,
                                       fields, smoothed_band_zmin=band_zmin)
    print(f"self-check  : {SELF_CHECK_NPTS} pts (seed {SELF_CHECK_SEED}) "
          f"median rel err {med:.3e} (tol {SELF_CHECK_MEDIAN_TOL:g}), "
          f"max rel err {mx:.3e} (tol {SELF_CHECK_MAX_TOL:g}) -> "
          f"{'PASS' if passed else 'FAIL'}{note}")

    ok = passed
    if args.compare_to is not None:
        ok &= compare_to_reference(args.out, args.compare_to)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
