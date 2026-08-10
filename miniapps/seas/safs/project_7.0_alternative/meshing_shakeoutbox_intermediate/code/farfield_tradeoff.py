#!/usr/bin/env python3
"""Where does the collar spend its cells, and what would a UNIFORM far field cost?

Two questions, both measured rather than argued:

  1. how the current collar's 22.27M cells are distributed with depth -- i.e.
     how much of it is the thin slow surface layer;
  2. for a UNIFORM far-field size h, the cell count and the resulting
     resolved-frequency compliance over the collar footprint.

The tension to expose: the 0.5 Hz p3 gate is f = 0.75*Vs/dx with Vs at the
element barycentre, and the CVM's near-surface Vs in the collar runs 186-560
m/s.  Holding the gate AT the free surface therefore forces dx <= 280-840 m
over 238,386 km2, which is exactly what makes the collar expensive.  A uniform
far-field size cannot do both; this quantifies the choice.
"""
import sys
import numpy as np
import h5py

CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")
GATE = 0.6667
TARGET_E = (72000.0, 786000.0)
TARGET_N = (3524000.0, 4007000.0)
CORNERS = np.array([[362120., 3996867.], [235495., 3777546.],
                    [599606., 3567327.], [726231., 3786648.]])
Z_BOT = -40000.0
CH = 2_000_000


def order_rect(C):
    c = C.mean(0)
    return C[np.argsort(np.arctan2(C[:, 1] - c[1], C[:, 0] - c[0]))]


def inside_rect(X, Y, C):
    C = order_rect(C)
    ins = np.ones(X.shape, bool)
    for i in range(4):
        a, b = C[i], C[(i + 1) % 4]
        e = b - a
        ins &= (e[0] * (Y - a[1]) - e[1] * (X - a[0])) >= 0.0
    return ins


with h5py.File(CVM, "r") as f:
    gx, gy, gz = f["x"][:], f["y"][:], f["z"][:]
    d = f["data"]
    VS = np.sqrt(np.asarray(d["mu"]).astype(np.float64)
                 / np.asarray(d["rho"]).astype(np.float64))

X, Y = np.meshgrid(gx, gy, indexing="xy")
in_box = (X >= TARGET_E[0]) & (X <= TARGET_E[1]) & (Y >= TARGET_N[0]) & (Y <= TARGET_N[1])
collar2d = in_box & ~inside_rect(X, Y, CORNERS)
ncol = int(collar2d.sum())
area = ncol * 2.25
print(f"collar footprint: {ncol:,} CVM columns = {area:,.0f} km2\n")

# ---- 1. where the current collar's cells actually are ----------------------
if len(sys.argv) > 1:
    with h5py.File(sys.argv[1], "r") as f:
        P = f["geometry"][:]
        conn = f["connect"]
        nt = conn.shape[0]
        nt0 = int(sys.argv[2])
        edges = np.array([0, -100, -250, -500, -1000, -2000, -4000, -10000, -40001])
        cnt = np.zeros(len(edges) - 1, np.int64)
        for s0 in range(nt0, nt, CH):
            T = conn[s0:min(s0 + CH, nt)].astype(np.int64)
            zb = P[T].mean(1)[:, 2]
            for k in range(len(edges) - 1):
                cnt[k] += int(((zb <= edges[k]) & (zb > edges[k + 1])).sum())
            del T, zb
    tot = cnt.sum()
    print(f"CURRENT collar: {tot:,} cells by barycentre depth")
    cum = 0
    for k in range(len(edges) - 1):
        cum += cnt[k]
        print(f"   {edges[k]:>7,} .. {edges[k+1]:>7,} m : {cnt[k]:>12,}  "
              f"({100*cnt[k]/tot:5.1f} %)   cumulative {100*cum/tot:5.1f} %")
    print()

# ---- 2. uniform far-field size: cost vs compliance -------------------------
# cells ~ 6 * volume / h^3 is the standard tet-per-cube count; the collar's
# own measured density is used as a sanity anchor where available.
print("UNIFORM far-field size h -> cost and 0.5 Hz p3 compliance")
print(f"{'h [m]':>7}{'cells':>14}{'vs 22.27M':>11}   "
      f"{'volume frac with Vs/h >= 0.6667':>34}")
vol_km3 = area * (0.0 - Z_BOT) / 1e3
kz = (gz >= Z_BOT) & (gz <= 0.0)
zsel = gz[kz]
V3 = VS[kz][:, collar2d]                       # (nz, ncol)
w = np.ones(len(zsel))                          # equal 250 m slabs
for h in (500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000):
    cells = 6.0 * vol_km3 * 1e9 / h ** 3
    ok = (V3 / h >= GATE)
    frac = float((ok * w[:, None]).sum() / (w.sum() * V3.shape[1]))
    # and how deep the failure reaches
    bad_z = zsel[(~ok).any(1)]
    print(f"{h:>7,}{cells:>14,.0f}{cells/22.27e6:>10.2f}x   {100*frac:>32.2f} %"
          f"   fails above z = {bad_z.max() if len(bad_z) else 0:,.0f} m")

print("\nVs over the collar footprint, by depth slab:")
for zt in (0, -125, -250, -500, -1000, -2000, -5000):
    k = int(np.argmin(np.abs(gz - zt)))
    s = VS[k][collar2d]
    print(f"   z {gz[k]:>7,.0f} m : Vs min {s.min():>7,.0f}  p10 {np.percentile(s,10):>7,.0f}"
          f"  med {np.median(s):>7,.0f}   -> dx for 0.5 Hz: "
          f"min {s.min()/GATE:>7,.0f}  med {np.median(s)/GATE:>7,.0f} m")
