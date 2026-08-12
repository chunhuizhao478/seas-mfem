#!/usr/bin/env python3
"""price_build.py -- what a from-scratch ShakeOut-box heavy mesh costs.

Design field, on a lattice over the whole domain box:

    h(x,y,z) = min( fault_band_envelope(d) , gate_ceiling(x,y,z) )

then Lipschitz-limited (|grad h| <= hgrad-1) so every size change is gradual by
construction, then integrated as K * sum(6 dV / h^3) with K measured for this
pipeline (3.66, calibrate_density.py).

The fault-band envelope is MEASURED off the deployed mesh, not assumed
(measure_envelope.py): ~110 m out to 1 km at every depth, then 244 / 357 / 518 /
761 m at 1.5 / 2 / 3 / 5 km.

The gate ceiling is the same self-consistent MUSCAL construction used to close
the gates -- the largest h whose own barycentre still measures Vs/dx >= gate.

Prints the full-spec price and a COARSE-BAND price, because the full spec is a
single-shot fill of ~10^8 cells that no 36 GB machine can run: the staged route
is to fill at a coarse band and recover the band by red/LEB refinement, which is
how this lineage was actually built (fb200 -> refine2 = fault edge/4).
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_gate" / "code"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))

TARGET_E = (72000.0, 786000.0)
TARGET_N = (3524000.0, 4007000.0)
Z_BOT = -40000.0
LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
BC_DR = 3
CH = 1_000_000
KFAC = 3.66

# measured envelope: (distance from fault [m], median max edge [m])
ENV_D = np.array([0.0, 400.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 8000.0])
ENV_H = np.array([110.0, 120.0, 220.0, 250.0, 360.0, 520.0, 760.0, 1200.0])


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def fault_centroids(mesh):
    out = []
    with h5py.File(mesh, "r") as f:
        G = f["geometry"][:]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                hit |= face_code(Bc, s) == BC_DR
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            for s in range(4):
                m = face_code(Bs, s) == BC_DR
                if m.any():
                    out.append(G[Tc[m][:, list(LOCAL_FACES[s])]].mean(1))
    return np.vstack(out)


def lipschitz(H, dx, dy, dzk, slope, iters=14):
    H = H.copy()
    dz = np.asarray(dzk, float).reshape(-1, 1, 1)
    for _ in range(iters):
        b = H.sum()
        H[1:] = np.minimum(H[1:], H[:-1] + slope * dz)
        H[:-1] = np.minimum(H[:-1], H[1:] + slope * dz)
        H[:, 1:] = np.minimum(H[:, 1:], H[:, :-1] + slope * dy)
        H[:, :-1] = np.minimum(H[:, :-1], H[:, 1:] + slope * dy)
        H[:, :, 1:] = np.minimum(H[:, :, 1:], H[:, :, :-1] + slope * dx)
        H[:, :, :-1] = np.minimum(H[:, :, :-1], H[:, :, 1:] + slope * dx)
        if abs(H.sum() - b) < 1e-7 * abs(b):
            break
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fault-mesh", required=True)
    ap.add_argument("--gate", type=float, default=0.8)
    ap.add_argument("--hgrad", type=float, default=1.3)
    ap.add_argument("--h-max", type=float, default=3000.0)
    ap.add_argument("--band-scale", type=float, action="append",
                    help="multiply the measured envelope by this (repeatable). "
                         "1 = the deployed design; 4 = a coarse base fill that a "
                         "red refinement then brings back to spec")
    ap.add_argument("--lat", type=float, default=1000.0, help="lattice spacing")
    a = ap.parse_args()
    scales = a.band_scale or [1.0, 2.0, 4.0]

    from muscal_vs import MuscalVs
    from collar_lib import VsGrid
    CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_"
           "THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_"
           "deep40km/safs_material_cvm.nc")
    mus = MuscalVs(backup=VsGrid(CVM))

    gx = np.arange(TARGET_E[0], TARGET_E[1] + 1, a.lat)
    gy = np.arange(TARGET_N[0], TARGET_N[1] + 1, a.lat)
    z = -mus.dep[(mus.dep >= 0) & (mus.dep <= -Z_BOT)]
    nz, ny, nx = len(z), len(gy), len(gx)
    print(f"lattice {nz} x {ny} x {nx} = {nz*ny*nx:,} nodes "
          f"({a.lat:,.0f} m lateral)", flush=True)

    print("[fault] building distance field…", flush=True)
    C = fault_centroids(a.fault_mesh)
    print(f"[fault] {len(C):,} facet centroids", flush=True)
    tree = cKDTree(C)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    D = np.empty((nz, ny, nx), np.float32)
    Q = np.empty((ny * nx, 3))
    Q[:, 0] = X.ravel()
    Q[:, 1] = Y.ravel()
    for k, zz in enumerate(z):
        Q[:, 2] = zz
        d, _ = tree.query(Q, distance_upper_bound=ENV_D[-1] * 3, workers=-1)
        D[k] = np.where(np.isfinite(d), d, ENV_D[-1] * 3).reshape(ny, nx)
    print("[gate] building the self-consistent ceiling…", flush=True)

    # gate ceiling: largest h in the contiguous feasible run (see field3d.py)
    V = np.empty((nz, ny, nx), np.float32)
    P = np.empty((ny * nx, 3))
    P[:, 0] = X.ravel()
    P[:, 1] = Y.ravel()
    for k, zz in enumerate(z):
        P[:, 2] = zz
        V[k] = mus.at(P).reshape(ny, nx)
    hs = np.arange(100.0, a.h_max + 1, 50.0)
    HG = np.full((nz, ny, nx), np.nan, np.float32)
    run = np.ones((nz, ny, nx), bool)
    for h in hs:
        zb = z[:, None] - 0.294 * h
        kk = np.abs(z[None, :] - zb).argmin(1)
        run &= V[kk] >= a.gate * h
        HG[run] = h
    HG[~np.isfinite(HG)] = 100.0

    edges = np.concatenate(([z[0]], 0.5 * (z[1:] + z[:-1]), [z[-1]]))
    thk = np.abs(np.diff(edges))
    dvk = (a.lat * a.lat * thk).reshape(-1, 1, 1)
    vol = float((dvk * np.ones((1, ny, nx))).sum())
    print(f"[vol] {vol/1e9:,.0f} km3\n")

    print(f"{'band scale':>12}{'band h @0':>11}{'raw cells':>16}{'graded cells':>16}")
    for sc in scales:
        HB = np.interp(D, ENV_D, ENV_H * sc).astype(np.float32)
        H = np.minimum(HB, HG)
        n_raw = float((6.0 * dvk / H.astype(np.float64) ** 3).sum()) * KFAC
        HL = lipschitz(H, a.lat, a.lat, thk, a.hgrad - 1.0)
        n_lim = float((6.0 * dvk / HL.astype(np.float64) ** 3).sum()) * KFAC
        print(f"{sc:>12.1f}{ENV_H[0]*sc:>11,.0f}{n_raw:>16,.0f}{n_lim:>16,.0f}")

    print(f"\n  deployed heavy, for reference: 160,825,939 tets")
    print(f"  (hgrad {a.hgrad}, K={KFAC} measured, envelope from measure_envelope.py)")


if __name__ == "__main__":
    main()
