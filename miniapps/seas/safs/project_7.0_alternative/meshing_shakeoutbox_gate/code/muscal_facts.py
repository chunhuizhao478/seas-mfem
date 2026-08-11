#!/usr/bin/env python3
"""muscal_facts.py -- MUSCAL's own grid and Vs, and how it differs from the deck nc.

The deck's `safs_material_cvm.nc` is a RESAMPLED product (1500 m lateral,
250 m uniform vertical).  MUSCAL is the source model, on 0.01 deg with a
depth axis that is fine near the surface.  Whether the frequency gate should be
judged against the source or against the resampled file is not a detail: the
deck nc's z = 0 bin bottoms out at Vs 155 m/s and is 3.78x slower than the bin
250 m below it, and that single step is what sets the entire collar's price.

Prints both grids side by side and, on the shared footprint, the near-surface
Vs each one reports -- including MUSCAL averaged over the deck's 250 m bin,
which is the like-for-like comparison.
"""
import argparse

import h5py
import numpy as np

MUSCAL = "/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc"
CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--muscal", default=MUSCAL)
    ap.add_argument("--cvm", default=CVM)
    a = ap.parse_args()

    with h5py.File(a.muscal, "r") as f:
        lon, lat, dep = f["longitude"][:], f["latitude"][:], f["depth"][:]
        print("MUSCAL")
        print(f"  longitude {len(lon):>5}  {lon.min():.4f} .. {lon.max():.4f}   "
              f"d {np.diff(lon).min():.4f} .. {np.diff(lon).max():.4f} deg")
        print(f"  latitude  {len(lat):>5}  {lat.min():.4f} .. {lat.max():.4f}   "
              f"d {np.diff(lat).min():.4f} .. {np.diff(lat).max():.4f} deg")
        print(f"  depth     {len(dep):>5}  {dep.min():,.1f} .. {dep.max():,.1f}")
        d = np.sort(dep)
        print(f"  depth axis head: {np.array2string(d[:16], precision=1, max_line_width=110)}")
        print(f"  d(depth) min {np.diff(d).min():,.1f}  max {np.diff(d).max():,.1f}")
        # near-surface levels only -- reading the whole cube is 1.4 GB
        order = np.argsort(dep)
        print("\n  Vs by MUSCAL level (shallowest 12):")
        print(f"    {'depth':>9}{'finite':>12}{'min':>8}{'p1':>8}{'p10':>8}{'p50':>8}{'p90':>8}")
        for k in order[:12]:
            v = f["vs"][k]
            m = np.isfinite(v) & (v > 0)
            if not m.any():
                print(f"    {dep[k]:>9,.1f}   all NaN")
                continue
            vv = v[m]
            print(f"    {dep[k]:>9,.1f}{m.sum():>12,}{vv.min():>8,.0f}"
                  f"{np.percentile(vv,1):>8,.0f}{np.percentile(vv,10):>8,.0f}"
                  f"{np.percentile(vv,50):>8,.0f}{np.percentile(vv,90):>8,.0f}")

    with h5py.File(a.cvm, "r") as f:
        cz = f["z"][:]
        dd = f["data"]
        mu = np.asarray(dd["mu"]).astype(np.float64)
        rho = np.asarray(dd["rho"]).astype(np.float64)
    vs = np.sqrt(np.maximum(mu, 0.0) / rho)
    print("\nDECK CVM (what the gate currently reads, and what ASAGI serves)")
    print(f"  z {len(cz)}  {cz.min():,.0f} .. {cz.max():,.0f}   uniform dz "
          f"{np.diff(np.sort(cz)).min():,.0f}")
    print(f"    {'z':>9}{'min':>8}{'p1':>8}{'p10':>8}{'p50':>8}{'p90':>8}")
    for zq in (0, -250, -500, -1000):
        k = int(np.argmin(np.abs(cz - zq)))
        v = vs[k][np.isfinite(vs[k]) & (vs[k] > 0)]
        print(f"    {cz[k]:>9,.0f}{v.min():>8,.0f}{np.percentile(v,1):>8,.0f}"
              f"{np.percentile(v,10):>8,.0f}{np.percentile(v,50):>8,.0f}"
              f"{np.percentile(v,90):>8,.0f}")


if __name__ == "__main__":
    main()
