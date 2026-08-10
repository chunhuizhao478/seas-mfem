#!/usr/bin/env python3
"""Minimum-gain far-field size profile h(z) for a resolved-frequency gate.

The gate is f = Vs/dx with Vs NEAREST-GRID AT THE ELEMENT BARYCENTRE.  That
makes the sizing self-referential near the surface: a cell of vertical extent h
sitting on the free surface has its barycentre at roughly z = -h/3, so the
constraint is

    h <= Vs(-h/3) / gate                                    (*)

Chasing Vs(0) instead -- as a pointwise field would -- drives h down, which
drags the barycentre INTO the slowest CVM bin and tightens the constraint
further.  That feedback is what made the gate-driven collar cost 22.27M cells
and STILL miss the gate on 5.5 % of them, while a uniform 2,500 m collar cost
13.4M and missed 0.315 %.

This solves (*) by fixed-point iteration on a lateral PERCENTILE of Vs, so the
profile is set by the bulk of the far field rather than by its slowest column
(which is a thin skin and cannot be resolved at any affordable cost).  Deeper
cells are then allowed to grow to whatever Vs permits, which is where the
saving is: uniform 2,500 m wastes cells below ~2 km, where Vs 2,800-3,300 m/s
permits 3.5-5 km.
"""
import argparse
import numpy as np
import h5py

CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")
TARGET_E = (72000.0, 786000.0)
TARGET_N = (3524000.0, 4007000.0)
CORNERS = np.array([[362120., 3996867.], [235495., 3777546.],
                    [599606., 3567327.], [726231., 3786648.]])
Z_BOT = -40000.0


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--pct", type=float, default=10.0,
                    help="lateral Vs percentile the profile is sized to")
    ap.add_argument("--h-min", type=float, default=400.0)
    ap.add_argument("--h-max", type=float, default=6000.0)
    a = ap.parse_args()

    with h5py.File(CVM, "r") as f:
        gx, gy, gz = f["x"][:], f["y"][:], f["z"][:]
        d = f["data"]
        VS = np.sqrt(np.asarray(d["mu"]).astype(np.float64)
                     / np.asarray(d["rho"]).astype(np.float64))
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    box = ((X >= TARGET_E[0]) & (X <= TARGET_E[1])
           & (Y >= TARGET_N[0]) & (Y <= TARGET_N[1]))
    collar = box & ~inside_rect(X, Y, CORNERS)
    kz = (gz >= Z_BOT) & (gz <= 0.0)
    z = gz[kz]
    prof = np.percentile(VS[kz][:, collar], a.pct, axis=1)     # Vs_pct(z)
    area_km2 = int(collar.sum()) * 2.25
    print(f"collar {area_km2:,.0f} km2   gate {a.gate}   Vs percentile p{a.pct:g}")

    def vs_at(zq):
        return np.interp(np.clip(zq, z.min(), z.max()), z, prof)

    # fixed point on h <= Vs(-h/3)/gate, marching down from the surface
    print(f"\n{'z [m]':>9}{'Vs_p':>9}{'h allowed':>11}{'cum cells':>13}")
    zc, cells, rows = 0.0, 0.0, []
    while zc > Z_BOT:
        h = a.h_max
        for _ in range(60):                       # fixed point
            h_new = np.clip(vs_at(zc - h / 3.0) / a.gate, a.h_min, a.h_max)
            if abs(h_new - h) < 1.0:
                h = h_new
                break
            h = 0.5 * (h + h_new)
        h = float(min(h, zc - Z_BOT if zc - Z_BOT > 0 else h))
        slab = min(h, zc - Z_BOT)
        cells += 6.0 * area_km2 * 1e6 * slab / h ** 3
        rows.append((zc, float(vs_at(zc - h / 3.0)), h, cells))
        zc -= slab
    for zc, v, h, c in rows[:14]:
        print(f"{zc:>9,.0f}{v:>9,.0f}{h:>11,.0f}{c:>13,.0f}")
    if len(rows) > 14:
        print(f"      ... {len(rows)-14} more layers")
    print(f"\nlayers {len(rows)}   TOTAL far-field cells ~ {cells:,.0f}")
    print(f"  vs uniform 2,500 m : {6.0*area_km2*1e6*40000/2500**3:,.0f}")
    print(f"  vs the gate-driven collar shipped earlier : 22,266,601")


if __name__ == "__main__":
    main()
