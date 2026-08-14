#!/usr/bin/env python3
"""lts_cost.py -- clustered rate-2 LTS cost of a PUML mesh.

The project rule (learned on the s=30 cut, 2026-07-27): NEVER judge a mesh
change by dt_min or by tet count alone.  A handful of tiny cells can drop
dt_min by 33x while costing only +2.7% under SeisSol's rate-2 LTS, because
those cells sit alone in the fastest cluster.  The decision number is

    cost = sum_i 1 / dt_cluster(i)          [element-updates per simulated second]

dt_i = CFL * 2 * r_insphere,i / ((2p+1) * vp,i)   (SeisSol ADER-DG)
cluster k = floor(log2(dt_i/dt_min)),  dt_cluster = dt_min * 2^k

vp is read from the deck's CVM NetCDF by nearest grid node at the element
barycenter -- the same nearest-grid convention the resolution gate uses.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import h5py

FACE_VERTS = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("puml", type=Path)
    ap.add_argument("--cvm", type=Path, required=True)
    ap.add_argument("--order", type=int, default=4, help="SeisSol ORDER (p+1)")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--max-clusters", type=int, default=0,
                    help="0 = unlimited (SeisSol default)")
    ap.add_argument("--label", default=None)
    a = ap.parse_args(argv)

    with h5py.File(a.puml, "r") as f:
        V = f["geometry"][:].astype(np.float64)
        T = f["connect"][:].astype(np.int64)
    P = V[T]
    d = P[:, 1:] - P[:, 0:1]
    vol = np.abs(np.einsum('ij,ij->i', np.cross(d[:, 0], d[:, 1]), d[:, 2])) / 6.0
    A = np.zeros(len(T))
    for f_ in FACE_VERTS:
        Q = P[:, list(f_)]
        A += 0.5 * np.linalg.norm(np.cross(Q[:, 1] - Q[:, 0], Q[:, 2] - Q[:, 0]), axis=1)
    r_in = 3.0 * vol / A

    with h5py.File(a.cvm, "r") as f:
        gx, gy, gz = f["x"][:], f["y"][:], f["z"][:]
        D = f["data"][:]
    bc = P.mean(axis=1)
    ix = np.clip(np.searchsorted(gx, bc[:, 0]), 0, len(gx) - 1)
    iy = np.clip(np.searchsorted(gy, bc[:, 1]), 0, len(gy) - 1)
    iz = np.clip(np.searchsorted(gz, bc[:, 2]), 0, len(gz) - 1)
    rho = D["rho"][iz, iy, ix].astype(np.float64)
    mu = D["mu"][iz, iy, ix].astype(np.float64)
    lam = D["lambda"][iz, iy, ix].astype(np.float64)
    vp = np.sqrt((lam + 2.0 * mu) / rho)

    p = a.order - 1
    dt = a.cfl * 2.0 * r_in / ((2 * p + 1) * vp)
    dtmin = float(dt.min())
    k = np.floor(np.log2(dt / dtmin)).astype(np.int64)
    if a.max_clusters:
        k = np.minimum(k, a.max_clusters - 1)
    dtc = dtmin * (2.0 ** k)
    cost = float((1.0 / dtc).sum())
    gts = len(T) / dtmin

    print(f"=== LTS cost: {a.label or a.puml.name} ===")
    print(f"  tets {len(T):,}   ORDER {a.order} (p={p})   CFL {a.cfl}")
    print(f"  vp  min {vp.min():,.0f}  med {np.median(vp):,.0f}  max {vp.max():,.0f} m/s")
    print(f"  r_insphere min {r_in.min():.3f}  med {np.median(r_in):,.1f} m")
    print(f"  dt_min {dtmin*1e6:,.3f} us   dt_med {np.median(dt)*1e6:,.1f} us")
    print(f"  clusters {int(k.max())+1}")
    nb = np.bincount(k)
    for i, n in enumerate(nb):
        if n:
            print(f"    cluster {i:2d}: {n:10,d} tets   dt {dtmin*2**i*1e6:9,.2f} us"
                  f"   cost {n/(dtmin*2**i):.4e}")
    print(f"  GTS cost {gts:.4e}   LTS cost {cost:.4e}   (LTS speedup {gts/cost:.2f}x)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
