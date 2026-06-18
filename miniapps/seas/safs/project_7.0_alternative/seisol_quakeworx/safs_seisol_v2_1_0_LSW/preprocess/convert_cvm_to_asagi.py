#!/usr/bin/env python3
"""convert_cvm_to_asagi.py — CVM sidecar HDF5 -> SeisSol/easi ASAGI NetCDF.

Reads a schema-v1 velocity sidecar (rectilinear UTM 11N grid of Vp/Vs/density,
non-uniform z), converts velocities to elastic moduli AT THE SOURCE GRID NODES

    mu     = rho * Vs^2
    lambda = rho * (Vp^2 - 2 * Vs^2)

then resamples the moduli along z onto a uniform axis (ASAGI requires
equidistant spacing per dimension).  Conversion happens before resampling so
the stored grid is exactly the piecewise-linear interpolant of the node
moduli — the same function SeisSol's linear ASAGI lookup evaluates.  With
the default dz = 50 m every source z-level in range is an exact output node
(all SAFS sidecar levels are multiples of 50 m), so no source structure is
smoothed away; a coarser --dz that misses source levels triggers a warning.

Output is a COARDS NetCDF with compound variable material{rho;mu;lambda}
laid out data(z, y, x), matching the layout SeisSol's ASAGI reader expects
(reference: SeisSol/Training sulawesi/3dvel_Sulawesi.nc and
SeisSol/preprocessing/science/generating_ASAGI_file.py).

A round-trip self-check always runs after writing: trilinear samples of the
written NetCDF are compared against an independent evaluation of the source
sidecar at 1000 fixed-seed points inside the SAFS mesh bounding box.  Exit
status is non-zero if the self-check fails.
"""

from __future__ import annotations

import argparse
import datetime
import os
import sys
from pathlib import Path

import h5py
import numpy as np

SCHEMA_VERSION = "data_projection_v1"
CANONICAL_CRS = "EPSG:32611"
CANONICAL_UNITS = "m"
CANONICAL_Z_POSITIVE = "elevation"

# SAFS SeisSol mesh bounding box (safs_mesh.puml.h5 geometry, UTM 11N, m).
# Used only by the round-trip self-check to draw sample points.
MESH_BBOX = (
    (314840.1, 672329.5),
    (3642278.0, 3888985.0),
    (-41607.61, 0.0),
)
SELF_CHECK_SEED = 20260609
SELF_CHECK_NPTS = 1000
SELF_CHECK_MEDIAN_TOL = 1.0e-6
SELF_CHECK_MAX_TOL = 5.0e-3


def read_sidecar(path):
    """Read a schema-v1 sidecar.

    Returns (x, y, z, fields, attrs) where fields maps
    {'Vp','Vs','density'} -> (nx,ny,nz) float64 arrays.
    Raises ValueError naming the missing/violating item on any schema error.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        attrs = {k: (v.decode() if isinstance(v, bytes) else v)
                 for k, v in h5.attrs.items()}
        required_attrs = {
            "schema_version": SCHEMA_VERSION,
            "crs": CANONICAL_CRS,
            "units": CANONICAL_UNITS,
            "z_positive": CANONICAL_Z_POSITIVE,
        }
        for key, expected in required_attrs.items():
            if key not in attrs:
                raise ValueError(
                    f"read_sidecar: '{path}' is missing required root "
                    f"attribute '{key}' (expected '{expected}')")
            if str(attrs[key]) != expected:
                raise ValueError(
                    f"read_sidecar: '{path}' attr '{key}' is "
                    f"'{attrs[key]}', expected '{expected}'")
        axes = {}
        for ax in ("x", "y", "z"):
            dset = f"grid/{ax}"
            if dset not in h5:
                raise ValueError(
                    f"read_sidecar: '{path}' is missing dataset '{dset}'")
            axes[ax] = np.asarray(h5[dset][:], dtype=np.float64)
        fields = {}
        for name in ("Vp", "Vs", "density"):
            dset = f"fields/{name}"
            if dset not in h5:
                raise ValueError(
                    f"read_sidecar: '{path}' is missing dataset '{dset}'")
            fields[name] = np.asarray(h5[dset][:], dtype=np.float64)

    expected_shape = (axes["x"].size, axes["y"].size, axes["z"].size)
    for name, arr in fields.items():
        if arr.shape != expected_shape:
            raise ValueError(
                f"read_sidecar: field '{name}' shape {arr.shape} != "
                f"grid shape {expected_shape}")
        if np.isnan(arr).any():
            raise ValueError(
                f"read_sidecar: field '{name}' contains NaN "
                f"(schema v1 forbids NaN)")
    return axes["x"], axes["y"], axes["z"], fields, attrs


def validate_axes(x, y, z):
    """x and y must be uniform (equidistant) and ascending; z strictly
    ascending (may be non-uniform — it gets resampled)."""
    for name, ax in (("x", x), ("y", y), ("z", z)):
        if ax.ndim != 1 or ax.size < 2:
            raise ValueError(
                f"validate_axes: axis '{name}' must be 1-D with >= 2 "
                f"points; got shape {ax.shape}")
        d = np.diff(ax)
        if not np.all(d > 0):
            raise ValueError(
                f"validate_axes: axis '{name}' is not strictly increasing")
    for name, ax in (("x", x), ("y", y)):
        d = np.diff(ax)
        if not np.allclose(d, d[0], rtol=1e-9, atol=1e-6):
            raise ValueError(
                f"validate_axes: axis '{name}' is not uniform "
                f"(spacing range [{d.min():.6g}, {d.max():.6g}]); ASAGI "
                f"requires equidistant spacing and this converter only "
                f"resamples z")


def resample_z(field, z_src, z_out):
    """Linear resampling of a (nx,ny,nz_src) field along z onto z_out.

    z_out must lie within [z_src[0], z_src[-1]] (interpolation only).
    Returns (nx, ny, len(z_out)) float64.
    """
    field = np.asarray(field, dtype=np.float64)
    z_src = np.asarray(z_src, dtype=np.float64)
    z_out = np.asarray(z_out, dtype=np.float64)
    if field.ndim != 3 or field.shape[2] != z_src.size:
        raise ValueError(
            f"resample_z: field shape {field.shape} incompatible with "
            f"z_src size {z_src.size}")
    if z_out[0] < z_src[0] or z_out[-1] > z_src[-1]:
        raise ValueError(
            f"resample_z: z_out [{z_out[0]}, {z_out[-1]}] extends outside "
            f"z_src [{z_src[0]}, {z_src[-1]}]; extrapolation is not allowed")
    idx = np.searchsorted(z_src, z_out, side="right") - 1
    idx = np.clip(idx, 0, z_src.size - 2)
    z0 = z_src[idx]
    z1 = z_src[idx + 1]
    w = np.clip((z_out - z0) / (z1 - z0), 0.0, 1.0)
    return field[:, :, idx] * (1.0 - w) + field[:, :, idx + 1] * w


def velocities_to_moduli(vp, vs, rho):
    """(Vp, Vs, density) -> (rho, mu, lambda), SI units.

    Raises ValueError if any resulting mu or lambda is non-positive.
    """
    vp = np.asarray(vp, dtype=np.float64)
    vs = np.asarray(vs, dtype=np.float64)
    rho = np.asarray(rho, dtype=np.float64)
    mu = rho * vs**2
    lam = rho * (vp**2 - 2.0 * vs**2)
    if not np.all(mu > 0):
        raise ValueError(
            f"velocities_to_moduli: mu must be positive everywhere; "
            f"min(mu) = {mu.min():.6g} (check Vs > 0, rho > 0)")
    if not np.all(lam > 0):
        n_bad = int((lam <= 0).sum())
        raise ValueError(
            f"velocities_to_moduli: lambda <= 0 at {n_bad} point(s) "
            f"(min {lam.min():.6g}); requires Vp > sqrt(2)*Vs everywhere "
            f"(min Vp/Vs = {(vp / vs).min():.4f})")
    return rho, mu, lam


def write_asagi_netcdf(path, x, y, z, rho, mu, lam, dtype, attrs):
    """Write the ASAGI COARDS NetCDF.

    Coordinate variables x/y/z are always float64; the compound
    material{rho;mu;lambda} members use `dtype` ('float32' or 'float64').
    Field arrays come in as (nx,ny,nz) and are stored as data(z,y,x).
    """
    import netCDF4

    np_dtype = np.dtype(dtype)
    if np_dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
        raise ValueError(
            f"write_asagi_netcdf: dtype must be float32 or float64, "
            f"got '{dtype}'")
    nx, ny, nz = x.size, y.size, z.size
    for name, arr in (("rho", rho), ("mu", mu), ("lambda", lam)):
        if arr.shape != (nx, ny, nz):
            raise ValueError(
                f"write_asagi_netcdf: field '{name}' shape {arr.shape} "
                f"!= ({nx}, {ny}, {nz})")

    material_t = np.dtype([("rho", np_dtype), ("mu", np_dtype),
                           ("lambda", np_dtype)])
    with netCDF4.Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.createDimension("x", nx)
        ds.createDimension("y", ny)
        ds.createDimension("z", nz)
        vx = ds.createVariable("x", "f8", ("x",))
        vy = ds.createVariable("y", "f8", ("y",))
        vz = ds.createVariable("z", "f8", ("z",))
        vx[:] = x
        vy[:] = y
        vz[:] = z
        mtype = ds.createCompoundType(material_t, "material")
        data = ds.createVariable("data", mtype, ("z", "y", "x"))
        buf = np.empty((nz, ny, nx), dtype=material_t)
        buf["rho"] = np.transpose(rho, (2, 1, 0)).astype(np_dtype)
        buf["mu"] = np.transpose(mu, (2, 1, 0)).astype(np_dtype)
        buf["lambda"] = np.transpose(lam, (2, 1, 0)).astype(np_dtype)
        data[:] = buf
        for key, val in attrs.items():
            ds.setncattr(key, val)


def read_asagi_netcdf(path):
    """Read back a file written by write_asagi_netcdf.

    Returns (x, y, z, {'rho','mu','lambda'} -> (nx,ny,nz) float64).
    """
    import netCDF4

    with netCDF4.Dataset(str(path), "r") as ds:
        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        y = np.asarray(ds.variables["y"][:], dtype=np.float64)
        z = np.asarray(ds.variables["z"][:], dtype=np.float64)
        raw = ds.variables["data"][:]
        fields = {}
        for name in ("rho", "mu", "lambda"):
            arr = np.asarray(raw[name], dtype=np.float64)  # (nz, ny, nx)
            fields[name] = np.transpose(arr, (2, 1, 0))    # -> (nx, ny, nz)
    return x, y, z, fields


def trilinear_sample(xg, yg, zg, field, pts):
    """Trilinear interpolation of field (nx,ny,nz) on a rectilinear grid
    (axes may be non-uniform) at pts (N,3).  Points must lie inside the
    grid bounding box."""
    pts = np.asarray(pts, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"trilinear_sample: pts must be (N,3); "
                         f"got {pts.shape}")
    for k, g in enumerate((xg, yg, zg)):
        lo, hi = pts[:, k].min(), pts[:, k].max()
        if lo < g[0] or hi > g[-1]:
            raise ValueError(
                f"trilinear_sample: points along axis {k} span "
                f"[{lo:.6g}, {hi:.6g}], outside grid [{g[0]:.6g}, "
                f"{g[-1]:.6g}]")

    def locate(g, c):
        i = np.searchsorted(g, c, side="right") - 1
        i = np.clip(i, 0, g.size - 2)
        w = np.clip((c - g[i]) / (g[i + 1] - g[i]), 0.0, 1.0)
        return i, w

    ix, wx = locate(np.asarray(xg, dtype=np.float64), pts[:, 0])
    iy, wy = locate(np.asarray(yg, dtype=np.float64), pts[:, 1])
    iz, wz = locate(np.asarray(zg, dtype=np.float64), pts[:, 2])
    out = np.zeros(pts.shape[0], dtype=np.float64)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                wgt = ((wx if dx else 1.0 - wx)
                       * (wy if dy else 1.0 - wy)
                       * (wz if dz else 1.0 - wz))
                out += wgt * field[ix + dx, iy + dy, iz + dz]
    return out


def _apply_vs_floor(vp, vs, vs_floor):
    """Clamp Vs up to vs_floor, then lift Vp to keep lambda > 0."""
    vs = np.maximum(vs, vs_floor)
    vp = np.maximum(vp, np.sqrt(2.0) * 1.001 * vs)
    return vp, vs


def self_check(nc_path, x_src, y_src, z_src, fields_src, vs_floor):
    """Round-trip check: trilinear samples of the written NetCDF vs an
    independent evaluation of the source sidecar, at SELF_CHECK_NPTS
    fixed-seed points in MESH_BBOX.
    Returns (median_rel_err, max_rel_err, passed).

    The reference converts velocities to moduli at the SOURCE grid nodes
    (the vs-floor, if any, applied at the nodes first — same as the
    converter) and interpolates trilinearly using the source's
    NON-UNIFORM z axis.  Because the converter resamples node moduli onto
    a uniform z grid that (at the default dz) contains every source
    level, both paths evaluate the same piecewise-linear function: any
    transposition, member-order, dimension-order, or resampling bug shows
    up as a large error, while a correct file agrees to float32 rounding.
    """
    rng = np.random.default_rng(SELF_CHECK_SEED)
    pts = np.column_stack([
        rng.uniform(lo, hi, SELF_CHECK_NPTS) for (lo, hi) in MESH_BBOX
    ])

    xg, yg, zg, nc_fields = read_asagi_netcdf(nc_path)

    vp_n, vs_n = fields_src["Vp"], fields_src["Vs"]
    if vs_floor is not None:
        vp_n, vs_n = _apply_vs_floor(vp_n, vs_n, vs_floor)
    rho_n, mu_n, lam_n = velocities_to_moduli(
        vp_n, vs_n, fields_src["density"])
    ref_axes = (x_src, y_src, z_src)

    ref = {"rho": rho_n, "mu": mu_n, "lambda": lam_n}
    rel_errs = []
    for name in ("rho", "mu", "lambda"):
        got = trilinear_sample(xg, yg, zg, nc_fields[name], pts)
        want = trilinear_sample(*ref_axes, ref[name], pts)
        rel_errs.append(np.abs(got - want) / np.abs(want))
    rel = np.concatenate(rel_errs)
    med = float(np.median(rel))
    mx = float(rel.max())
    passed = (med <= SELF_CHECK_MEDIAN_TOL) and (mx <= SELF_CHECK_MAX_TOL)
    return med, mx, passed


def _plot_check(nc_path, x_src, y_src, z_src, fields_src, out_png):
    """Optional side-by-side mu slice (z = -1000 m): NetCDF vs source."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"plot-check skipped: matplotlib unavailable ({exc})")
        return
    z_slice = -1000.0
    xg, yg, zg, nc_fields = read_asagi_netcdf(nc_path)
    xx, yy = np.meshgrid(xg, yg, indexing="ij")
    pts = np.column_stack([xx.ravel(), yy.ravel(),
                           np.full(xx.size, z_slice)])
    mu_nc = trilinear_sample(xg, yg, zg, nc_fields["mu"], pts)
    _, mu_n, _ = velocities_to_moduli(
        fields_src["Vp"], fields_src["Vs"], fields_src["density"])
    mu_src = trilinear_sample(x_src, y_src, z_src, mu_n, pts)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, dat, title in ((axes[0], mu_nc, "mu from NetCDF"),
                           (axes[1], mu_src, "mu from sidecar")):
        im = ax.pcolormesh(xg, yg, dat.reshape(xg.size, yg.size).T,
                           shading="auto")
        ax.set_title(f"{title} (z = {z_slice:.0f} m)")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"plot-check written: {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sidecar", required=True,
                    help="path to schema-v1 velocity sidecar HDF5")
    ap.add_argument("--out", default="safs_material_cvm.nc",
                    help="output ASAGI NetCDF path (default: %(default)s)")
    ap.add_argument("--z-min", type=float, default=-45000.0,
                    help="bottom of the uniform output z axis, m "
                         "(default: %(default)s)")
    ap.add_argument("--z-max", type=float, default=100.0,
                    help="top of the uniform output z axis, m "
                         "(default: %(default)s)")
    ap.add_argument("--dz", type=float, default=50.0,
                    help="uniform output z spacing, m (default: "
                         "%(default)s — captures every SAFS sidecar "
                         "z-level exactly). If dz does not divide "
                         "z_max - z_min, the last level is the largest "
                         "z_min + k*dz <= z_max. A dz that misses "
                         "source levels smooths the model there "
                         "(warning issued).")
    ap.add_argument("--vs-floor", type=float, default=None,
                    help="optional Vs clamp (m/s) applied at the source "
                         "grid nodes before moduli conversion; Vp is "
                         "lifted to sqrt(2)*1.001*Vs where needed "
                         "(default: off)")
    ap.add_argument("--dtype", choices=("float32", "float64"),
                    default="float32",
                    help="storage dtype of the material members "
                         "(default: %(default)s)")
    ap.add_argument("--plot-check", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="write a PNG comparing a mu slice from the "
                         "NetCDF vs the sidecar")
    args = ap.parse_args(argv)

    x, y, z, fields, src_attrs = read_sidecar(args.sidecar)
    validate_axes(x, y, z)

    if args.z_min < z[0] or args.z_max > z[-1]:
        ap.error(
            f"--z-min/--z-max [{args.z_min}, {args.z_max}] must lie within "
            f"the source z range [{z[0]}, {z[-1]}] (interpolation only, "
            f"no extrapolation)")
    if args.dz <= 0:
        ap.error(f"--dz must be positive; got {args.dz}")
    if args.z_min >= args.z_max:
        ap.error(f"--z-min ({args.z_min}) must be < --z-max ({args.z_max})")
    if args.vs_floor is not None and args.vs_floor <= 0:
        ap.error(f"--vs-floor must be positive; got {args.vs_floor}")

    z_out = np.arange(args.z_min, args.z_max + 0.5 * args.dz, args.dz)
    if z_out[-1] > args.z_max:  # guard float drift past z_max
        z_out = z_out[:-1]

    print(f"source grid : {x.size} x {y.size} x {z.size} "
          f"(z in [{z[0]:.0f}, {z[-1]:.0f}], non-uniform)")
    print(f"output grid : {x.size} x {y.size} x {z_out.size} "
          f"(z in [{z_out[0]:.0f}, {z_out[-1]:.0f}] @ dz = {args.dz:g})")

    # Source z-levels in range that the uniform grid cannot represent are
    # smoothed over by the resampling — warn so this is a deliberate choice.
    in_range = z[(z >= args.z_min) & (z <= z_out[-1])]
    k = np.rint((in_range - args.z_min) / args.dz)
    missed = in_range[~np.isclose(args.z_min + k * args.dz, in_range,
                                  atol=1e-6 * args.dz)]
    if missed.size:
        print(f"WARNING     : {missed.size} source z-level(s) not on the "
              f"output grid (model smoothed there): "
              f"{', '.join(f'{v:.0f}' for v in missed)}")

    vp_n, vs_n = fields["Vp"], fields["Vs"]
    if args.vs_floor is not None:
        n_clamped = int((vs_n < args.vs_floor).sum())
        vp_n, vs_n = _apply_vs_floor(vp_n, vs_n, args.vs_floor)
        print(f"vs-floor    : {args.vs_floor:g} m/s "
              f"({n_clamped} of {vs_n.size} source nodes clamped)")
    rho_n, mu_n, lam_n = velocities_to_moduli(vp_n, vs_n,
                                              fields["density"])

    rho = resample_z(rho_n, z, z_out)
    mu = resample_z(mu_n, z, z_out)
    lam = resample_z(lam_n, z, z_out)
    for name, arr in (("rho", rho), ("mu", mu), ("lambda", lam)):
        if np.isnan(arr).any():
            raise ValueError(
                f"resampled field '{name}' contains NaN — source data or "
                f"resampling is corrupt")

    out_attrs = {
        "source_sidecar": os.path.abspath(args.sidecar),
        "source_created_at": str(src_attrs.get("created_at", "unknown")),
        "source_mesh_tag": str(src_attrs.get("mesh_tag", "unknown")),
        "crs": CANONICAL_CRS,
        "z_positive": CANONICAL_Z_POSITIVE,
        "units": "SI: rho kg/m^3, mu/lambda Pa, axes m",
        "dz": float(args.dz),
        "vs_floor": (float(args.vs_floor) if args.vs_floor is not None
                     else "none"),
        "converter": "convert_cvm_to_asagi.py",
        "converted_at": datetime.datetime.now(datetime.timezone.utc)
                        .strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    write_asagi_netcdf(args.out, x, y, z_out, rho, mu, lam,
                       args.dtype, out_attrs)
    size_mb = os.path.getsize(args.out) / 1e6

    vp_eff = np.sqrt((lam + 2.0 * mu) / rho)
    vs_eff = np.sqrt(mu / rho)
    print(f"written     : {args.out} ({size_mb:.1f} MB, dtype {args.dtype})")
    print( "field       :        min            max")
    print(f"  rho       : {rho.min():14.6g} {rho.max():14.6g}  kg/m^3")
    print(f"  mu        : {mu.min():14.6g} {mu.max():14.6g}  Pa")
    print(f"  lambda    : {lam.min():14.6g} {lam.max():14.6g}  Pa")
    print(f"  Vs (impl.): {vs_eff.min():14.6g} {vs_eff.max():14.6g}  m/s")
    print(f"  Vp (impl.): {vp_eff.min():14.6g} {vp_eff.max():14.6g}  m/s")

    med, mx, passed = self_check(args.out, x, y, z, fields,
                                 args.vs_floor)
    print(f"self-check  : {SELF_CHECK_NPTS} pts (seed {SELF_CHECK_SEED}) "
          f"median rel err {med:.3e} (tol {SELF_CHECK_MEDIAN_TOL:g}), "
          f"max rel err {mx:.3e} (tol {SELF_CHECK_MAX_TOL:g})")
    print(f"self-check  : {'PASS' if passed else 'FAIL'}")

    if args.plot_check:
        _plot_check(args.out, x, y, z, fields,
                    Path(args.out).with_suffix("").name + "_check.png")

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
