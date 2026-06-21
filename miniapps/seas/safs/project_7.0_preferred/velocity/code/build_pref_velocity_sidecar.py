#!/usr/bin/env python3
"""Build the PREFERRED MFEM velocity sidecar (velocity_safs.h5, data_projection_v1)
from the SAME material nc SeisSol consumes (safs_material_cvm.nc, compound
{rho, mu, lambda}).  MFEM and SeisSol therefore sample identical material.

  Vp = sqrt((lambda + 2 mu) / rho),  Vs = sqrt(mu / rho),  density = rho
Far-field mesh extent beyond this grid is handled by the gated OOBPolicy::Clamp
(ASAGI-style), so the grid only needs to cover the fault/data hull.
"""
import sys
import numpy as np
import netCDF4

REPO = "/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver"
SEAS = f"{REPO}/miniapps/seas"
ALT = f"{SEAS}/safs/project_7.0_alternative"
PREF = f"{SEAS}/safs/project_7.0_preferred"
sys.path.insert(0, f"{ALT}/velocity/code")
from sidecar import write_sidecar  # noqa: E402

NC = f"{ALT}/seisol_quakeworx/safs_seisol_v3_0_0_LSW_PREFERRED/safs_material_cvm.nc"
OUT = f"{PREF}/velocity/results/multiscale_statewise_cvm/velocity_safs.h5"

ds = netCDF4.Dataset(NC)
x = np.asarray(ds.variables["x"][:], dtype=np.float64)
y = np.asarray(ds.variables["y"][:], dtype=np.float64)
z = np.asarray(ds.variables["z"][:], dtype=np.float64)
data = ds.variables["data"][:]              # compound, shape (z, y, x)

def field_zyx_to_xyz(name):
    a = np.asarray(data[name], dtype=np.float64)   # (z, y, x)
    if np.ma.isMaskedArray(a):
        a = a.filled(np.nan)
    return np.transpose(a, (2, 1, 0))              # -> (x, y, z)

rho = field_zyx_to_xyz("rho")
mu = field_zyx_to_xyz("mu")
lam = field_zyx_to_xyz("lambda")
print(f"grid: nx={x.size} ny={y.size} nz={z.size}")
print(f"  x[{x.min():.1f},{x.max():.1f}] y[{y.min():.1f},{y.max():.1f}] z[{z.min():.1f},{z.max():.1f}]")

# Physicality guard BEFORE deriving velocities (report, don't silently fill).
for nm, arr in (("rho", rho), ("mu", mu), ("lambda", lam)):
    nnan = int(np.isnan(arr).sum())
    print(f"  {nm}: min={np.nanmin(arr):.4g} max={np.nanmax(arr):.4g} NaN={nnan}")
    if nnan:
        idx = np.argwhere(np.isnan(arr))[:5]
        sys.exit(f"FATAL: {nnan} NaN cells in '{nm}' (e.g. {idx.tolist()}); "
                 "ASAGI input should be gap-free — investigate before writing sidecar.")
if np.nanmin(rho) <= 0 or np.nanmin(mu) < 0:
    sys.exit("FATAL: non-physical rho<=0 or mu<0 in material nc.")

Vp = np.sqrt((lam + 2.0 * mu) / rho)
Vs = np.sqrt(mu / rho)
density = rho
for nm, arr in (("Vp", Vp), ("Vs", Vs), ("density", density)):
    print(f"  {nm}: min={arr.min():.4g} max={arr.max():.4g} NaN={int(np.isnan(arr).sum())}")

# Pad bounds slightly so edge cells satisfy strict min<=cell<=max in the reader.
def bnd(a, units):
    lo, hi = float(a.min()), float(a.max())
    pad = max(1.0, 1e-6 * (hi - lo))
    return (lo - pad, hi + pad, units)

fields = {"Vp": Vp, "Vs": Vs, "density": density}
bounds = {"Vp": bnd(Vp, "m/s"), "Vs": bnd(Vs, "m/s"),
          "density": bnd(density, "kg/m^3")}
attrs = {
    "schema_version": "data_projection_v1",
    "crs": "EPSG:32611",
    "units": "m",
    "z_positive": "elevation",
    "source": ("safs_material_cvm.nc (PREFERRED muscal CVM, the SeisSol ASAGI "
               "material) -> Vp/Vs/density via build_pref_velocity_sidecar.py"),
    "mesh_tag": "safv4_deep_500m_opt",
}
import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
write_sidecar(OUT, x, y, z, fields, attrs, bounds)
print(f"WROTE {OUT}")
print(f"  size = {os.path.getsize(OUT)/1e6:.1f} MB")
