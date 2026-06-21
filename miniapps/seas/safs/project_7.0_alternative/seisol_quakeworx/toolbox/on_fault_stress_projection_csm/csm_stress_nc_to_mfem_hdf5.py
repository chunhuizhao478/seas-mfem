#!/usr/bin/env python3
"""Convert a SeisSol-convention CSM initial-stress netCDF into an MFEM-SEAS
`StressField3D` HDF5 sidecar (schema `data_projection_v1`).

Why this exists
---------------
The SAFS v3-ALT MFEM port (`PLAN_mixed_flux_seissol_port_2026-06-19.md`, Phase 1)
runs the SAME spatially-varying CSM stress field that the SeisSol QuakeWorx case
reads through ASAGI.  SeisSol reads `safs_stress_csm_d12.nc` (compression-NEGATIVE,
effective — hydrostatic P_p and the gate-removal overpressure already baked in).
MFEM's `seas_spatial_dyn_driver` reads a different on-disk format (the HDF5
`data_projection_v1` schema consumed by `io/stress_field_3d.cpp` /
`io/data_field_3d.cpp`) and uses the OPPOSITE sign convention.

Sign convention (THE load-bearing detail)
------------------------------------------
  * SeisSol / the nc : compression-NEGATIVE (geophysics), effective Pa.
  * MFEM             : compression-POSITIVE (`sigma_n > 0 == compression`,
                       miniapps/seas/CLAUDE.md), and `StressField3D::Evaluate`
                       is a STRICT PASS-THROUGH (no sign flip in C++ — the flip
                       MUST happen here, in the writer).
  => emit  sigma_MFEM[ij] = -s_seisol[ij]  for ALL SIX components.

There is NO special-casing of `sigma_xy`: the global `-1` flip is applied to
every component, and MFEM's per-DOF fault basis (`field_coefficient.cpp:472`,
`tau2 = t2 . (S.n)`) handles the on-fault rotation / rake.  The strike rake is
invariant to the normal orientation because `strike = up x n` is bilinear in n
(verified in REVIEW_mixed_flux_port_plan_2026-06-19.md).

Pore pressure
-------------
The nc field is EFFECTIVE (overpressure baked in spatially), so the MFEM config
MUST set `[pore_pressure] P_p_pa = 0.0` — `ProjectFaultPreStress`
(`field_coefficient.cpp:489-503`) subtracts P_p UNCONDITIONALLY; a non-zero P_p
here would double-subtract.  This script does not touch P_p; it only flips sign.

Usage
-----
  python3 csm_stress_nc_to_mfem_hdf5.py \
      --in  ../../safs_seisol_v3_0_0_LSW_ALT/safs_stress_csm_d12.nc \
      --out <repo>/.../stress/results/.../safs_stress_csm_d12_mfem.h5 \
      [--pad-m 0.0] \
      [--validate csm_yhsm2013_stress_on_safs_mesh_B_fault_stress.vtu]

The emitted `.h5` is large and gitignored — document its provenance and scp it to
Frontera alongside `velocity_safs.h5`.
"""

import argparse
import sys

import numpy as np

# Source (nc) compound-field names -> MFEM (h5) canonical dataset names.
# Order is the StressField3D canonical order (stress_field_3d.cpp:47-52).
NAME_MAP = [
    ("s_xx", "sigma_xx"),
    ("s_yy", "sigma_yy"),
    ("s_zz", "sigma_zz"),
    ("s_xy", "sigma_xy"),
    ("s_yz", "sigma_yz"),
    ("s_xz", "sigma_xz"),
]

HYPO_DEFAULT = (606971.0, 3707270.0, -4965.62)  # SAFS hypocenter, UTM 11N (m)


# ---------------------------------------------------------------------------
# netCDF read
# ---------------------------------------------------------------------------
def load_csm_nc(path):
    """Read the SeisSol-convention CSM stress nc.

    Returns (x, y, z, comps) where x/y/z are 1-D float64 axes (strictly
    increasing) and comps maps the SIX MFEM dataset names to (Nx, Ny, Nz)
    float64 arrays already SIGN-FLIPPED to compression-positive effective Pa.
    """
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit(
            "csm_stress_nc_to_mfem_hdf5: needs the `netCDF4` Python package to "
            "read the compound nc (conda activate pythonenv). Import error: %s" % exc
        )

    ds = Dataset(str(path), "r")
    try:
        ds.set_always_mask(False)
    except Exception:
        pass

    for axis in ("x", "y", "z"):
        if axis not in ds.variables:
            raise SystemExit(
                "csm_stress_nc_to_mfem_hdf5: nc '%s' is missing coordinate "
                "variable '%s' (expected x, y, z)." % (path, axis)
            )
    if "data" not in ds.variables:
        raise SystemExit(
            "csm_stress_nc_to_mfem_hdf5: nc '%s' is missing the compound 'data' "
            "variable (expected fields s_xx..s_xz)." % path
        )

    x = np.asarray(ds.variables["x"][:], dtype=np.float64)
    y = np.asarray(ds.variables["y"][:], dtype=np.float64)
    z = np.asarray(ds.variables["z"][:], dtype=np.float64)

    # data is a compound array shaped (nz, ny, nx) with float32 fields.
    raw = ds.variables["data"][:]
    src_names = set(raw.dtype.names or ())
    missing = [s for s, _ in NAME_MAP if s not in src_names]
    if missing:
        raise SystemExit(
            "csm_stress_nc_to_mfem_hdf5: nc 'data' compound is missing fields %s "
            "(has %s)." % (missing, sorted(src_names))
        )

    nx, ny, nz = len(x), len(y), len(z)
    comps = {}
    for src, dst in NAME_MAP:
        a = np.asarray(raw[src], dtype=np.float64)          # (nz, ny, nx)
        if a.shape != (nz, ny, nx):
            raise SystemExit(
                "csm_stress_nc_to_mfem_hdf5: field '%s' has shape %s, expected "
                "(nz,ny,nx)=%s." % (src, a.shape, (nz, ny, nx))
            )
        a = np.transpose(a, (2, 1, 0))                      # -> (nx, ny, nz)
        # SIGN FLIP: compression-negative effective -> compression-positive.
        comps[dst] = -a
    ds.close()
    return x, y, z, comps


# ---------------------------------------------------------------------------
# optional edge padding (nearest-edge extrapolation, == ASAGI clamp behaviour)
# ---------------------------------------------------------------------------
def pad_grid(x, y, z, comps, pad_m):
    """Extend every axis outward by `pad_m` metres, replicating the edge slice.

    StressField3D HARD-ABORTS on any out-of-bbox fault DOF (data_field_3d.cpp:616).
    The nc is built to cover the fault, but MFEM's strict OOB check is less
    forgiving than ASAGI's clamp.  This reproduces ASAGI's nearest-edge clamp by
    one (or more) replicated cells per face so a fault DOF exactly on / slightly
    past the nc boundary does not abort.  Default pad_m=0 keeps the nc extent.
    """
    if pad_m <= 0.0:
        return x, y, z, comps

    def extend_axis(ax):
        d0 = ax[1] - ax[0]
        d1 = ax[-1] - ax[-2]
        n_lo = int(np.ceil(pad_m / d0))
        n_hi = int(np.ceil(pad_m / d1))
        lo = ax[0] - d0 * np.arange(n_lo, 0, -1)
        hi = ax[-1] + d1 * np.arange(1, n_hi + 1)
        return np.concatenate([lo, ax, hi]), n_lo, n_hi

    xx, nxl, nxh = extend_axis(x)
    yy, nyl, nyh = extend_axis(y)
    zz, nzl, nzh = extend_axis(z)
    out = {}
    for k, a in comps.items():
        a = np.pad(a, ((nxl, nxh), (nyl, nyh), (nzl, nzh)), mode="edge")
        out[k] = a
    return xx, yy, zz, out


# ---------------------------------------------------------------------------
# HDF5 write (data_projection_v1)
# ---------------------------------------------------------------------------
def write_mfem_hdf5(path, x, y, z, comps):
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "csm_stress_nc_to_mfem_hdf5: needs the `h5py` Python package to write "
            "the sidecar. Import error: %s" % exc
        )

    for ax, lab in ((x, "x"), (y, "y"), (z, "z")):
        if not np.all(np.diff(ax) > 0.0):
            raise SystemExit(
                "csm_stress_nc_to_mfem_hdf5: axis '%s' is not strictly increasing "
                "(StressField3D requires monotone-increasing axes)." % lab
            )
    nx, ny, nz = len(x), len(y), len(z)

    with h5py.File(str(path), "w") as f:
        # Root attributes required by DataField3D (data_field_3d.cpp:213-240).
        f.attrs["schema_version"] = "data_projection_v1"
        f.attrs["crs"] = "EPSG:32611"
        f.attrs["units"] = "m"
        f.attrs["z_positive"] = "elevation"
        f.attrs["source"] = "csm_stress_nc_to_mfem_hdf5.py (sign-flipped from CSM nc)"
        f.attrs["sign_convention"] = "compression-POSITIVE effective Pa (MFEM)"

        g = f.create_group("grid")
        g.create_dataset("x", data=x.astype(np.float64))
        g.create_dataset("y", data=y.astype(np.float64))
        g.create_dataset("z", data=z.astype(np.float64))

        fields = f.create_group("fields")
        for _, dst in NAME_MAP:
            a = np.ascontiguousarray(comps[dst].reshape((nx, ny, nz)), dtype=np.float64)
            if np.isnan(a).any():
                raise SystemExit(
                    "csm_stress_nc_to_mfem_hdf5: field '%s' contains NaN; the v1 "
                    "schema forbids NaN." % dst
                )
            vmin = float(a.min())
            vmax = float(a.max())
            if not (vmin < vmax):
                # A genuinely constant component (e.g. all-zero s_yz/s_xz) would
                # trip the strict reader check; nudge max so min<max holds.
                vmax = vmin + 1.0
            d = fields.create_dataset(dst, data=a)
            d.attrs["units"] = "Pa"
            d.attrs["min_value"] = float(a.min())
            d.attrs["max_value"] = vmax


# ---------------------------------------------------------------------------
# trilinear sampling (mirrors DataField3D's trilinear path for self-checks)
# ---------------------------------------------------------------------------
def trilinear(x, y, z, field, px, py, pz):
    def locate(ax, p):
        if p < ax[0] or p > ax[-1]:
            raise ValueError(
                "sample point %g outside axis bbox [%g, %g]" % (p, ax[0], ax[-1])
            )
        i = int(np.searchsorted(ax, p) - 1)
        i = min(max(i, 0), len(ax) - 2)
        t = (p - ax[i]) / (ax[i + 1] - ax[i])
        return i, t

    ix, tx = locate(x, px)
    iy, ty = locate(y, py)
    iz, tz = locate(z, pz)
    c = 0.0
    for dx in (0, 1):
        wx = tx if dx else 1.0 - tx
        for dy in (0, 1):
            wy = ty if dy else 1.0 - ty
            for dz in (0, 1):
                wz = tz if dz else 1.0 - tz
                c += wx * wy * wz * field[ix + dx, iy + dy, iz + dz]
    return c


def tensor_at(x, y, z, comps, p):
    s = {dst: trilinear(x, y, z, comps[dst], p[0], p[1], p[2]) for _, dst in NAME_MAP}
    return np.array(
        [
            [s["sigma_xx"], s["sigma_xy"], s["sigma_xz"]],
            [s["sigma_xy"], s["sigma_yy"], s["sigma_yz"]],
            [s["sigma_xz"], s["sigma_yz"], s["sigma_zz"]],
        ]
    )


# ---------------------------------------------------------------------------
# --validate against the projected on-fault VTU
# ---------------------------------------------------------------------------
def validate_against_vtu(x, y, z, comps, vtu_path):
    """Sample the sign-flipped field at each fault-facet centroid, project onto
    the facet normal, and compare sigma_n_eff / |tau| to the VTU's per-cell
    arrays.  Reports median / p95 relative error.  The driver's --print-derived
    + cycle-0 fault dump remain the authoritative acceptance gate."""
    try:
        import meshio

        m = meshio.read(vtu_path)
        pts = np.asarray(m.points, dtype=np.float64)
        tris = None
        for cb in m.cells:
            if cb.type in ("triangle", "triangle3"):
                tris = np.asarray(cb.data, dtype=np.int64)
                break
        if tris is None:
            raise SystemExit("validate: VTU has no triangle cells.")
        cd = m.cell_data
        def cell_arr(name):
            for k, v in cd.items():
                if k == name:
                    return np.concatenate([np.asarray(a) for a in v])
            # some writers store as point_data on cells; fall back
            if name in m.point_data:
                return np.asarray(m.point_data[name])
            return None
        ref_sn = cell_arr("sigma_n_eff_MPa_cell")
        ref_tau = cell_arr("tau_magnitude_MPa_cell")
    except ImportError:
        print("  [validate] skipped: `meshio` not importable (conda activate pythonenv).")
        return

    if ref_sn is None or ref_tau is None:
        print("  [validate] VTU lacks sigma_n_eff_MPa_cell / tau_magnitude_MPa_cell; "
              "skipping projection comparison.")
        return

    cen = pts[tris].mean(axis=1)
    v1 = pts[tris[:, 1]] - pts[tris[:, 0]]
    v2 = pts[tris[:, 2]] - pts[tris[:, 0]]
    nrm = np.cross(v1, v2)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)

    sn = np.full(len(cen), np.nan)
    tau = np.full(len(cen), np.nan)
    n_oob = 0
    for i in range(len(cen)):
        try:
            S = tensor_at(x, y, z, comps, cen[i])
        except ValueError:
            n_oob += 1
            continue
        n = nrm[i]
        t = S @ n
        snv = float(n @ t)                       # compression-positive normal stress
        tauv = float(np.linalg.norm(t - snv * n))
        sn[i] = abs(snv) / 1.0e6                  # -> MPa, magnitude
        tau[i] = tauv / 1.0e6

    ok = ~np.isnan(sn)
    def rel(a, b):
        m = ok & (np.abs(b) > 1.0)
        r = np.abs(a[m] - b[m]) / np.abs(b[m])
        return np.median(r), np.percentile(r, 95)
    sn_med, sn_p95 = rel(sn, np.asarray(ref_sn, dtype=np.float64))
    tau_med, tau_p95 = rel(tau, np.asarray(ref_tau, dtype=np.float64))
    print("  [validate] facets=%d  OOB=%d" % (len(cen), n_oob))
    print("  [validate] sigma_n_eff rel-err  median=%.2e  p95=%.2e  (target <=5e-3)"
          % (sn_med, sn_p95))
    print("  [validate] |tau|       rel-err  median=%.2e  p95=%.2e  (target <=1e-2)"
          % (tau_med, tau_p95))
    if n_oob:
        print("  [validate] WARNING: %d centroids OOB — increase --pad-m or the nc "
              "extent, else the driver will abort on the full mesh." % n_oob)


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_nc", required=True,
                    help="source SeisSol-convention CSM stress netCDF")
    ap.add_argument("--out", dest="out_h5", required=True,
                    help="output MFEM StressField3D HDF5 sidecar")
    ap.add_argument("--pad-m", type=float, default=0.0,
                    help="extend every grid axis outward by this many metres via "
                         "edge-slice replication (default 0 = keep nc extent)")
    ap.add_argument("--validate", dest="validate_vtu", default=None,
                    help="projected on-fault VTU to compare sigma_n_eff/|tau| against")
    ap.add_argument("--hypo", type=float, nargs=3, default=list(HYPO_DEFAULT),
                    metavar=("X", "Y", "Z"),
                    help="hypocenter sample point for the sanity print")
    args = ap.parse_args(argv)

    print("csm_stress_nc_to_mfem_hdf5: reading %s" % args.in_nc)
    x, y, z, comps = load_csm_nc(args.in_nc)
    print("  grid: nx=%d ny=%d nz=%d  x[%g,%g] y[%g,%g] z[%g,%g]"
          % (len(x), len(y), len(z), x[0], x[-1], y[0], y[-1], z[0], z[-1]))

    if args.pad_m > 0.0:
        x, y, z, comps = pad_grid(x, y, z, comps, args.pad_m)
        print("  padded to nx=%d ny=%d nz=%d (pad_m=%g)"
              % (len(x), len(y), len(z), args.pad_m))

    # Hypocenter sanity print (compression-positive effective tensor, MPa).
    hp = tuple(args.hypo)
    try:
        S = tensor_at(x, y, z, comps, hp) / 1.0e6
        print("  hypocenter (%g,%g,%g) effective tensor [MPa, comp+]:" % hp)
        print("    sxx=%.2f syy=%.2f szz=%.2f sxy=%.2f syz=%.2f sxz=%.2f"
              % (S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[1, 2], S[0, 2]))
    except ValueError as exc:
        print("  hypocenter sample skipped: %s" % exc)

    print("csm_stress_nc_to_mfem_hdf5: writing %s" % args.out_h5)
    write_mfem_hdf5(args.out_h5, x, y, z, comps)
    print("  wrote 6 fields (sigma_xx..sigma_xz), compression-POSITIVE effective Pa.")
    print("  REMINDER: the MFEM config must set [pore_pressure] P_p_pa = 0.0 "
          "(field is already effective).")

    if args.validate_vtu:
        print("csm_stress_nc_to_mfem_hdf5: validating against %s" % args.validate_vtu)
        validate_against_vtu(x, y, z, comps, args.validate_vtu)
    return 0


if __name__ == "__main__":
    sys.exit(main())
