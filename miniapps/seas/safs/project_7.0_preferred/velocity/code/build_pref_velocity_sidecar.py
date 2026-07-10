#!/usr/bin/env python3
"""Build the PREFERRED MFEM velocity sidecar (velocity_safs.h5, data_projection_v1)
from the SAME material nc SeisSol consumes (safs_material_cvm.nc, compound
{rho, mu, lambda}).  MFEM and SeisSol therefore sample identical material.

  Vp = sqrt((lambda + 2 mu) / rho),  Vs = sqrt(mu / rho),  density = rho
Far-field mesh extent beyond this grid is handled by the gated OOBPolicy::Clamp
(ASAGI-style), so the grid only needs to cover the fault/data hull.

Phase 4 of
  miniapps/seas/safs/project_7.0_preferred/document/
    PLAN_thermal_case2_mixedflux_port_2026-07-08.md

The source nc used to be HARD-CODED to
    {ALT}/seisol_quakeworx/safs_seisol_v3_0_0_LSW_PREFERRED/safs_material_cvm.nc
which no longer exists (`seisol_quakeworx` moved to `safs/seisol_quakeworx/` on
2026-06-26) and which, in any case, was the v3.0.0 CVM cube (297 x 196 x 194) --
NOT the statewide cube (452 x 361 x 194) the v3_4_1 decks ship.  Both `--in` and
`--out` are now REQUIRED so a stale default can never silently produce a sidecar
that disagrees with the deck SeisSol actually ran.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]          # .../miniapps/seas
ALT = REPO / "safs" / "project_7.0_alternative"
_SIDECAR_DIR = ALT / "velocity" / "code"
if not (_SIDECAR_DIR / "sidecar.py").is_file():
    # R-004: `parents[4]` is a positional assumption about this file's depth.
    # If the script is moved or copied (which is how these builders get used),
    # sys.path.insert would silently succeed and `import sidecar` would either
    # raise a bare ModuleNotFoundError or, worse, bind to SOME OTHER sidecar.py
    # on sys.path and emit a wrong-schema file.  Fail here, naming the path.
    raise SystemExit(
        "%s: cannot find the shared schema-v1 writer at\n"
        "    %s\n"
        "This script locates it relative to its own path "
        "(parents[4] must be miniapps/seas).\n"
        "Run it from its home in the repo, or fix REPO above."
        % (__file__.rsplit("/", 1)[-1], _SIDECAR_DIR / "sidecar.py")
    )
sys.path.insert(0, str(_SIDECAR_DIR))
from sidecar import write_sidecar  # noqa: E402


def field_zyx_to_xyz(data, name):
    a = np.asarray(data[name], dtype=np.float64)   # (z, y, x)
    if np.ma.isMaskedArray(a):
        a = a.filled(np.nan)
    return np.ascontiguousarray(np.transpose(a, (2, 1, 0)))   # -> (x, y, z)


def bnd(a, units):
    """Pad bounds slightly so edge cells satisfy strict min<=cell<=max."""
    lo, hi = float(a.min()), float(a.max())
    pad = max(1.0, 1e-6 * (hi - lo))
    return (lo - pad, hi + pad, units)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Convert the SeisSol ASAGI CVM material nc into an MFEM "
                    "data_projection_v1 velocity sidecar (Vp/Vs/density).")
    ap.add_argument("--in", dest="src", required=True, type=Path,
                    help="safs_material_cvm.nc (ASAGI, compound "
                         "data{rho, mu, lambda})")
    ap.add_argument("--out", dest="dst", required=True, type=Path,
                    help="destination velocity_safs.h5")
    ap.add_argument("--mesh-tag", default="safv4_deep19km_sub500m_flattop_"
                                          "nwtrim_fixed_opt",
                    help="informational `mesh_tag` root attribute")
    args = ap.parse_args(argv)

    if not args.src.is_file():
        raise SystemExit("build_pref_velocity_sidecar: --in '%s' not found"
                         % args.src)

    try:
        import netCDF4
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "build_pref_velocity_sidecar: needs the `netCDF4` Python package. "
            "Import error: %s" % exc)

    ds = netCDF4.Dataset(str(args.src))
    try:
        for axis in ("x", "y", "z"):
            if axis not in ds.variables:
                raise SystemExit(
                    "build_pref_velocity_sidecar: '%s' has no '%s' coordinate "
                    "variable; is it an ASAGI grid?" % (args.src, axis))
        if "data" not in ds.variables:
            raise SystemExit(
                "build_pref_velocity_sidecar: '%s' has no 'data' variable "
                "(expected an ASAGI compound {rho, mu, lambda})." % args.src)

        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        y = np.asarray(ds.variables["y"][:], dtype=np.float64)
        z = np.asarray(ds.variables["z"][:], dtype=np.float64)
        data = ds.variables["data"][:]              # compound, shape (z, y, x)

        members = getattr(data.dtype, "names", None) or ()
        for want in ("rho", "mu", "lambda"):
            if want not in members:
                raise SystemExit(
                    "build_pref_velocity_sidecar: compound 'data' is missing "
                    "member '%s'; has %s." % (want, list(members)))

        rho = field_zyx_to_xyz(data, "rho")
        mu = field_zyx_to_xyz(data, "mu")
        lam = field_zyx_to_xyz(data, "lambda")
    finally:
        ds.close()

    print(f"grid: nx={x.size} ny={y.size} nz={z.size}")
    print(f"  x[{x.min():.1f},{x.max():.1f}] y[{y.min():.1f},{y.max():.1f}] "
          f"z[{z.min():.1f},{z.max():.1f}]")

    # Physicality guard BEFORE deriving velocities (report, don't silently fill).
    for nm, arr in (("rho", rho), ("mu", mu), ("lambda", lam)):
        nnan = int(np.isnan(arr).sum())
        print(f"  {nm}: min={np.nanmin(arr):.4g} max={np.nanmax(arr):.4g} "
              f"NaN={nnan}")
        if nnan:
            idx = np.argwhere(np.isnan(arr))[:5]
            raise SystemExit(
                f"FATAL: {nnan} NaN cells in '{nm}' (e.g. {idx.tolist()}); "
                "ASAGI input should be gap-free — investigate before writing "
                "sidecar.")
    if np.nanmin(rho) <= 0 or np.nanmin(mu) < 0:
        raise SystemExit("FATAL: non-physical rho<=0 or mu<0 in material nc.")
    # Vp needs lambda + 2 mu > 0 (P-wave modulus); a zero/negative modulus would
    # yield NaN velocities that write_sidecar would then reject far downstream.
    pmod = lam + 2.0 * mu
    if np.nanmin(pmod) <= 0:
        raise SystemExit(
            "FATAL: non-physical P-wave modulus lambda + 2*mu <= 0 "
            f"(min = {np.nanmin(pmod):.4g}) in material nc.")

    Vp = np.sqrt(pmod / rho)
    Vs = np.sqrt(mu / rho)
    density = rho
    for nm, arr in (("Vp", Vp), ("Vs", Vs), ("density", density)):
        print(f"  {nm}: min={arr.min():.4g} max={arr.max():.4g} "
              f"NaN={int(np.isnan(arr).sum())}")

    fields = {"Vp": Vp, "Vs": Vs, "density": density}
    bounds = {"Vp": bnd(Vp, "m/s"), "Vs": bnd(Vs, "m/s"),
              "density": bnd(density, "kg/m^3")}
    attrs = {
        "schema_version": "data_projection_v1",
        "crs": "EPSG:32611",
        "units": "m",
        "z_positive": "elevation",
        "source": (f"{args.src.name} (muscal CVM, the SeisSol ASAGI material) "
                   "-> Vp/Vs/density via build_pref_velocity_sidecar.py"),
        "mesh_tag": args.mesh_tag,
    }
    args.dst.parent.mkdir(parents=True, exist_ok=True)
    write_sidecar(args.dst, x, y, z, fields, attrs, bounds)
    print(f"WROTE {args.dst}")
    print(f"  size = {os.path.getsize(args.dst)/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
