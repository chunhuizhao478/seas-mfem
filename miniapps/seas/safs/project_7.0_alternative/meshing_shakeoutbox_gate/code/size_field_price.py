#!/usr/bin/env python3
"""size_field_price.py -- what does actually holding the gate cost, per column?

The gate reads Vs NEAREST-GRID AT THE BARYCENTRE, and the CVM's z = 0 bin is
3.78x slower (median) than the bin 250 m below it.  So for a cell of vertical
extent h resting on the flat free surface, whose barycentre sits at
z = -beta*h, the constraint is

    h <= Vs(-beta*h) / gate                                        (*)

which is NOT monotone in h.  Making a cell smaller drags its barycentre UP into
the slow surface bin and tightens its own constraint.  (*) therefore admits two
disjoint families of solutions in a slow column:

  * a FINE branch,  h <= Vs(0)/gate          (barycentre stays in the z=0 bin);
  * a COARSE branch, h large enough that -beta*h reaches competent rock.

Between them is a FORBIDDEN ZONE where no cell size passes.  That is the whole
story of this mesh: the shipped uniform 2,500 m collar sits on the coarse
branch, which is why it beat a much finer gate-driven collar 29x.  The residual
failures are exactly the columns whose coarse branch does not exist.

This script scans h per column and prices every option, so the build target is
chosen from measured numbers instead of from a percentile guess.
"""
import argparse

import h5py
import numpy as np

CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")
TARGET_E = (72000.0, 786000.0)
TARGET_N = (3524000.0, 4007000.0)
# the ALT parent footprint: a rectangle rotated ~30 deg to the SAF strike
CORNERS = np.array([[362120., 3996867.], [235495., 3777546.],
                    [599606., 3567327.], [726231., 3786648.]])
Z_BOT = -40000.0
CELL_KM2 = 2.25            # one CVM column footprint, 1500 x 1500 m


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
    ap.add_argument("--cvm", default=CVM)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--beta", type=float, default=0.25,
                    help="barycentre depth / max edge for the SHALLOWEST cell of a "
                         "layer; 0.25 is the conservative (worst-cell) value")
    ap.add_argument("--h-min", type=float, default=150.0)
    ap.add_argument("--h-max", type=float, default=3000.0)
    ap.add_argument("--region", choices=("collar", "parent", "box"), default="collar")
    a = ap.parse_args()

    with h5py.File(a.cvm, "r") as f:
        gx, gy, gz = f["x"][:], f["y"][:], f["z"][:]
        d = f["data"]
        VS = np.sqrt(np.maximum(np.asarray(d["mu"]).astype(np.float64), 0.0)
                     / np.asarray(d["rho"]).astype(np.float64))
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    box = ((X >= TARGET_E[0]) & (X <= TARGET_E[1])
           & (Y >= TARGET_N[0]) & (Y <= TARGET_N[1]))
    inpar = inside_rect(X, Y, CORNERS)
    mask = {"collar": box & ~inpar, "parent": box & inpar, "box": box}[a.region]
    ncol = int(mask.sum())
    area_km2 = ncol * CELL_KM2

    # Vs profile per selected column, ordered top-down over [Z_BOT, 0]
    kz = (gz >= Z_BOT) & (gz <= 0.0)
    z = gz[kz]
    o = np.argsort(-z)
    z = z[o]                                    # descending: 0, -250, ...
    V = VS[kz][o][:, mask]                      # (nz, ncol)
    nz = len(z)
    print(f"region {a.region}: {ncol:,} CVM columns = {area_km2:,.0f} km2 "
          f"({100*mask.mean():.1f} % of the grid)")
    print(f"gate {a.gate}  beta {a.beta}  h in [{a.h_min:,.0f}, {a.h_max:,.0f}] m\n")

    def vs_at(zq, col):
        """nearest-grid Vs at depth zq for every column (vectorised over columns)."""
        k = np.abs(z[:, None] - zq[None, :]).argmin(0) if np.ndim(zq) else \
            int(np.abs(z - zq).argmin())
        return V[k, col] if np.ndim(zq) else V[k, col]

    # ---- 1. surface layer: scan h and find the feasible branches -----------
    hs = np.arange(a.h_min, a.h_max + 1.0, 25.0)
    zk = np.abs(z[:, None] - (-a.beta * hs)[None, :]).argmin(0)     # (nh,) bin index
    feas = V[zk, :] >= a.gate * hs[:, None]                          # (nh, ncol)

    any_feas = feas.any(0)
    hmax_feas = np.where(any_feas, hs[np.where(feas, hs[:, None], -1).argmax(0)], np.nan)
    # fine branch = feasible sizes contiguous from h_min upward
    fine_ok = feas[0]
    fine_top = np.full(ncol, np.nan)
    run = np.ones(ncol, bool)
    for i in range(len(hs)):
        run &= feas[i]
        fine_top[run] = hs[i]
    # coarse branch = the largest feasible h, when it is NOT part of the fine run
    coarse = np.where(any_feas & ~np.isclose(np.nan_to_num(fine_top, nan=-1), hmax_feas),
                      hmax_feas, np.nan)

    print("  SURFACE LAYER -- feasible cell sizes per column")
    print(f"    columns with NO feasible h at all      : {int((~any_feas).sum()):>8,} "
          f"({100*(~any_feas).mean():5.2f} %)  -> {int((~any_feas).sum())*CELL_KM2:,.0f} km2")
    print(f"    columns with a fine branch (h>={a.h_min:,.0f})   : {int(fine_ok.sum()):>8,} "
          f"({100*fine_ok.mean():5.2f} %)")
    print(f"    columns with a coarse branch           : {int(np.isfinite(coarse).sum()):>8,} "
          f"({100*np.isfinite(coarse).mean():5.2f} %)")
    for q in (1, 10, 50, 90):
        print(f"      largest feasible surface h  p{q:<2d} : "
              f"{np.nanpercentile(hmax_feas, q):>8,.0f} m")
    ft = fine_top[np.isfinite(fine_top)]
    if ft.size:
        print(f"      fine-branch ceiling Vs(0)/gate  p1 {np.percentile(ft,1):,.0f}  "
              f"p10 {np.percentile(ft,10):,.0f}  med {np.median(ft):,.0f} m")

    # how does the shipped uniform 2,500 m collar fare here?
    for h0 in (1500.0, 2000.0, 2500.0, 3000.0):
        i = int(np.abs(hs - h0).argmin())
        print(f"    uniform h = {h0:,.0f} m  -> surface layer passes in "
              f"{100*feas[i].mean():5.2f} % of columns "
              f"({int((~feas[i]).sum())*CELL_KM2:,.0f} km2 failing)")

    # ---- 2. full column march, cheapest feasible mesh ----------------------
    # Greedy top-down: at each z_top take the LARGEST feasible h (cheapest cell).
    print("\n  FULL COLUMN -- greedy cheapest feasible size field")
    cells = np.zeros(ncol)
    ztop = np.zeros(ncol)
    stuck = np.zeros(ncol, bool)
    nlayer = np.zeros(ncol, int)
    hsurf = np.full(ncol, np.nan)
    for it in range(400):
        live = (ztop > Z_BOT) & ~stuck
        if not live.any():
            break
        cols = np.nonzero(live)[0]
        zt = ztop[cols]
        # candidate barycentre depth for each (h, column)
        best = np.full(len(cols), np.nan)
        for i in range(len(hs) - 1, -1, -1):
            h = hs[i]
            need = np.isnan(best)
            if not need.any():
                break
            zb = zt[need] - a.beta * h
            k = np.abs(z[:, None] - zb[None, :]).argmin(0)
            ok = V[k, cols[need]] >= a.gate * h
            idx = np.nonzero(need)[0][ok]
            best[idx] = h
        nf = np.isnan(best)
        stuck[cols[nf]] = True
        good = ~nf
        c, h = cols[good], best[good]
        slab = np.minimum(h, ztop[c] - Z_BOT)
        cells[c] += 6.0 * CELL_KM2 * 1e6 * slab / h ** 3
        if it == 0:
            hsurf[c] = h
        ztop[c] -= slab
        nlayer[c] += 1
    print(f"    columns that could not be filled : {int(stuck.sum()):,} "
          f"({100*stuck.mean():.2f} %)")
    tot = np.nansum(cells[~stuck])
    print(f"    layers per column  med {int(np.median(nlayer)):,}  max {int(nlayer.max()):,}")
    print(f"    surface h chosen   p1 {np.nanpercentile(hsurf,1):,.0f}  "
          f"p10 {np.nanpercentile(hsurf,10):,.0f}  med {np.nanmedian(hsurf):,.0f}  "
          f"p90 {np.nanpercentile(hsurf,90):,.0f} m")
    print(f"    ANALYTIC far-field cells ~ {tot:,.0f}")
    print(f"      (the shipped uniform 2,500 m collar analytic estimate is "
          f"{6.0*area_km2*1e6*abs(Z_BOT)/2500**3:,.0f} for a MEASURED 13,405,498 -- "
          f"the build runs ~4.6x the estimate, so read these as ~4.6x too)")
    print(f"    => expect roughly {4.6*tot:,.0f} real cells")

    # ---- 3. price of ONLY fixing the columns the coarse branch misses ------
    print("\n  TARGETED FIX -- keep 2,500 m where it works, fine-branch elsewhere")
    i25 = int(np.abs(hs - 2500.0).argmin())
    bad = ~feas[i25]
    print(f"    columns needing the fine branch : {int(bad.sum()):,} "
          f"({100*bad.mean():.2f} %) = {int(bad.sum())*CELL_KM2:,.0f} km2")
    hb = fine_top[bad]
    hb = np.where(np.isfinite(hb), hb, a.h_min)
    # a fine surface skin of thickness = 1 cell, then grade back to 2,500 m
    skin = 6.0 * CELL_KM2 * 1e6 / hb ** 2
    print(f"      their fine ceiling  p1 {np.percentile(hb,1):,.0f}  "
          f"p10 {np.percentile(hb,10):,.0f}  med {np.median(hb):,.0f} m")
    print(f"      one-cell fine skin costs ~ {skin.sum():,.0f} analytic "
          f"(~{4.6*skin.sum():,.0f} real) cells")


if __name__ == "__main__":
    main()
