#!/usr/bin/env python3
"""beta_stats.py -- the barycentre-depth / max-edge ratio, per block and per size.

beta decides whether a free-surface cell's barycentre lands in the CVM's slow
z = 0 bin or in competent rock 250 m down, and the collar's price swings by
almost 2x on it (28.6M cells at beta 0.294 vs 46.5M at 0.124).  The first
measurement mixed the parent's 390 m surface cells with the collar's 2,265 m
ones, which is not the statistic a collar size field needs.

Also validates the beta-free alternative: pooling Vs over the cell's own
VERTICAL extent is a lower bound on Vs at ANY barycentre inside the cell, so
`dx <= Vs_pooled/gate` is a sufficient condition that no beta assumption can
break.  Report how conservative that is against the gate actually measured.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))
from collar_lib import VsGrid, tet_edge_lengths  # noqa: E402

CH = 1_000_000
CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--cvm", default=CVM)
    ap.add_argument("--first-tet", type=int, default=0)
    ap.add_argument("--gate", type=float, default=0.8)
    a = ap.parse_args()

    vs = VsGrid(a.cvm)
    with h5py.File(a.mesh, "r") as f:
        G = f["geometry"][:]
        nt = f["connect"].shape[0]
    print(f"{Path(a.mesh).name}: tets [{a.first_tet:,}, {nt:,})", flush=True)

    beta, bdx, bpass = [], [], []
    npool_fail = ngate_fail = nsurf = 0
    with h5py.File(a.mesh, "r") as f:
        conn = f["connect"]
        for s0 in range(a.first_tet, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            P = G[T]
            zmax = P[:, :, 2].max(1)
            m = zmax > -1e-6                      # rests on the flat free surface
            if not m.any():
                del T, P
                continue
            Pm = P[m]
            dmax = tet_edge_lengths(G, T[m]).max(1)
            bary = Pm.mean(1)
            b = (-bary[:, 2]) / dmax
            # measured gate
            r = vs.at(bary) / dmax
            # beta-free sufficient condition: min Vs over the cell's z extent
            zlo = Pm[:, :, 2].min(1)
            vpool = np.array([vs.column_min(bary[i:i + 1, :2], zlo[i], 0.0)[0]
                              for i in range(0, 0)]) if False else None
            nsurf += int(m.sum())
            ngate_fail += int((r < a.gate).sum())
            beta.append(b.astype(np.float32))
            bdx.append(dmax.astype(np.float32))
            bpass.append((r >= a.gate))
            del T, P, Pm

    B = np.concatenate(beta)
    D = np.concatenate(bdx)
    OK = np.concatenate(bpass)
    print(f"\n  free-surface-resting cells in this block: {len(B):,}   "
          f"gate {a.gate} failures {int((~OK).sum()):,}")
    print("\n  beta = |barycentre z| / max edge")
    print(f"    ALL          p1 {np.percentile(B,1):.3f}  p10 {np.percentile(B,10):.3f}  "
          f"med {np.median(B):.3f}  p90 {np.percentile(B,90):.3f}")
    print("\n  by cell size:")
    print(f"    {'dx band':<16}{'n':>12}{'beta p1':>10}{'p10':>8}{'med':>8}{'p90':>8}"
          f"{'fail':>10}")
    bands = [(0, 250), (250, 500), (500, 1000), (1000, 1500), (1500, 2000),
             (2000, 2500), (2500, 3000), (3000, 1e9)]
    for lo, hi in bands:
        m = (D >= lo) & (D < hi)
        if not m.any():
            continue
        print(f"    {lo:>6,}-{hi if hi<1e9 else 0:<9,.0f}{int(m.sum()):>12,}"
              f"{np.percentile(B[m],1):>10.3f}{np.percentile(B[m],10):>8.3f}"
              f"{np.median(B[m]):>8.3f}{np.percentile(B[m],90):>8.3f}"
              f"{int((~OK[m]).sum()):>10,}")

    # the barycentre depth that matters: does it clear the z=0 bin (|z|<125)?
    zb = B * D
    print(f"\n  barycentre depth of surface cells: "
          f"in the CVM z=0 bin (|z| < 125 m): {int((zb < 125).sum()):,} "
          f"({100*(zb<125).mean():.1f} %)")
    for lo, hi in ((0, 125), (125, 375), (375, 875), (875, 1e9)):
        m = (zb >= lo) & (zb < hi)
        if m.any():
            print(f"    bary {lo:>5,} .. {hi if hi<1e9 else 0:<7,.0f} m : "
                  f"{int(m.sum()):>10,}  failures {int((~OK[m]).sum()):>9,} "
                  f"({100*(~OK[m]).mean():5.2f} %)")


if __name__ == "__main__":
    main()
