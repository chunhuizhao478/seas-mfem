#!/usr/bin/env python3
"""Facts about the CVM grid that the frequency gate is judged against.

Everything downstream (the size field, the LEB target, the affordability
argument) depends on the CVM's actual lattice and on how fast Vs recovers with
depth, so measure it once and quote it rather than re-deriving it per script.
"""
import argparse

import h5py
import numpy as np

CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cvm", default=CVM)
    a = ap.parse_args()

    with h5py.File(a.cvm, "r") as f:
        x, y, z = f["x"][:], f["y"][:], f["z"][:]
        d = f["data"]
        mu = np.asarray(d["mu"]).astype(np.float64)
        rho = np.asarray(d["rho"]).astype(np.float64)
    vs = np.sqrt(np.maximum(mu, 0.0) / rho)

    print(f"file {a.cvm}")
    print(f"  x {len(x):>5}  {x.min():>12,.0f} .. {x.max():>12,.0f}   "
          f"dx {np.diff(x).min():,.1f} .. {np.diff(x).max():,.1f}")
    print(f"  y {len(y):>5}  {y.min():>12,.0f} .. {y.max():>12,.0f}   "
          f"dy {np.diff(y).min():,.1f} .. {np.diff(y).max():,.1f}")
    print(f"  z {len(z):>5}  {z.min():>12,.0f} .. {z.max():>12,.0f}")
    print(f"  data shape {vs.shape}  (z, y, x)")

    print("\n  z axis (top 20, descending from the top of the model):")
    zs = z[np.argsort(-z)][:20]
    print("   ", np.array2string(zs, precision=0, max_line_width=100))
    dz = np.diff(np.sort(z))
    print(f"  dz min {dz.min():,.0f}  max {dz.max():,.0f}")

    print("\n  Vs by level -- this is what makes the gate self-referential:")
    print(f"    {'z [m]':>9}{'min':>8}{'p1':>8}{'p10':>8}{'p50':>8}{'p90':>8}")
    for zq in (3250, 0, -125, -250, -500, -1000, -2000, -4000, -8000, -20000):
        k = int(np.argmin(np.abs(z - zq)))
        v = vs[k]
        v = v[np.isfinite(v) & (v > 0)]
        if v.size == 0:
            print(f"    {z[k]:>9,.0f}   (no finite data)")
            continue
        print(f"    {z[k]:>9,.0f}{v.min():>8,.0f}{np.percentile(v,1):>8,.0f}"
              f"{np.percentile(v,10):>8,.0f}{np.percentile(v,50):>8,.0f}"
              f"{np.percentile(v,90):>8,.0f}")

    # The ratio between the surface bin and the one below is the size of the
    # trap: refining a surface cell drags its barycentre from the second bin
    # into the first, and the required dx drops by this factor.
    k0 = int(np.argmin(np.abs(z - 0.0)))
    order = np.argsort(-z)
    below = order[list(order).index(k0) + 1] if list(order).index(k0) + 1 < len(order) else k0
    v0, v1 = vs[k0], vs[below]
    m = np.isfinite(v0) & np.isfinite(v1) & (v0 > 0) & (v1 > 0)
    r = v1[m] / v0[m]
    print(f"\n  step from z={z[k0]:,.0f} to z={z[below]:,.0f}: "
          f"Vs ratio median {np.median(r):.2f}, p90 {np.percentile(r,90):.2f}, max {r.max():.2f}")


if __name__ == "__main__":
    main()
