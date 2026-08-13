#!/usr/bin/env python3
"""surface_field.py -- the 2-D size field for the PLC's TOP surface.

h_surf(x,y) = min( h_fault(d_trace) , Vs_surface / gate , hmax ), graded to |grad h| <= g.

The gate term is what makes the free surface expensive and it is not optional: at
1 Hz for p5 the requirement is Vs/dx >= 0.8, and MUSCAL's surface Vs falls to
~200 m/s in the basins, so h must drop to ~250 m over large areas regardless of
distance from the fault.  Vs is pooled over the top 150 m rather than sampled at
depth 0, because a surface triangle belongs to cells that extend downward.
"""
import numpy as np, netCDF4 as ncdf
from pyproj import Transformer
from scipy.spatial import cKDTree

GATE, HMAX, G = 0.8, 5000.0, 0.15
PD = np.array([0, 125, 250, 500, 1000, 2000, 4000, 7500, 15000], float)
PH = np.array([115, 115, 141, 215, 262, 529, 700, 1000, 1500], float)
DXY = 750.0

def build(corners, chains_P, steep_P=None, muscal="/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc"):
    W, N, E, S = corners
    u = (S - W) / np.linalg.norm(S - W); v = (N - W) / np.linalg.norm(N - W)
    L1, L2 = np.linalg.norm(S - W), np.linalg.norm(N - W)
    ns, nt = int(L1 / DXY) + 1, int(L2 / DXY) + 1
    ss = np.linspace(0, L1, ns); tt = np.linspace(0, L2, nt)
    SS, TT = np.meshgrid(ss, tt, indexing="ij")
    XY = W[None, None, :] + SS[..., None] * u + TT[..., None] * v
    xy = XY.reshape(-1, 2)
    d, _ = cKDTree(chains_P[:, :2]).query(xy, k=1, distance_upper_bound=25000.0, workers=-1)
    d = np.where(np.isfinite(d), d, 25000.0)
    hf = np.interp(d, PD, PH, right=PH[-1])
    hf = np.where(d > 15000, PH[-1] + (d - 15000) * 0.20, hf)
    hf = np.where(d >= 24999.0, np.inf, hf)
    # The trace is at z = 0.00 at the 1st, 50th AND 99th percentile, but 37 of 6,575
    # nodes (0.56 %, in 4 clusters) excurse to -49.9..+25.8 m.  A 115 m triangle cannot
    # follow a 50 m notch, and the fault triangles just beneath then poke through the
    # free surface -- that is the entire self-intersection (one 0.5 km patch, measured).
    # Refine locally there; the cost is a few thousand triangles.
    if steep_P is not None and len(steep_P):
        ds, _ = cKDTree(steep_P[:, :2]).query(xy, k=1, distance_upper_bound=600.0)
        hf = np.minimum(hf, np.where(np.isfinite(ds), 40.0 + 0.15 * np.nan_to_num(ds, posinf=0.0), np.inf))
        print(f"[surf field] steep-trace refinement at {len(steep_P)} nodes "
              f"({int(np.isfinite(ds).sum()):,} grid points affected)")
    m = ncdf.Dataset(muscal)
    mlon = m["longitude"][:].data.astype(np.float64); mlat = m["latitude"][:].data.astype(np.float64)
    mdep = m["depth"][:].data.astype(np.float64)
    fwd = Transformer.from_crs("EPSG:32611", "EPSG:4326", always_xy=True)
    lo, la = fwd.transform(xy[:, 0], xy[:, 1])
    def near(ax, val):
        i = np.clip(np.searchsorted(ax, val), 0, len(ax) - 1); im = np.maximum(i - 1, 0)
        return np.where(np.abs(ax[im] - val) < np.abs(ax[i] - val), im, i).astype(np.int32)
    ii, jj = near(mlon, lo), near(mlat, la)
    kz = [int(near(mdep, np.array([float(z)]))[0]) for z in (0, 50, 100, 150)]
    vs = np.full(len(xy), np.inf, np.float32)
    for k in sorted(set(kz)):
        layer = np.asarray(m["vs"][k], np.float32)[jj, ii]
        vs = np.fmin(vs, np.where(np.isfinite(layer), layer, np.inf))
    hg = np.where(np.isfinite(vs), vs / GATE, HMAX)
    H = np.minimum(np.minimum(hf, hg), HMAX).reshape(ns, nt)
    for _ in range(400):
        b = H.copy()
        np.minimum(H[1:], H[:-1] + G * DXY, out=H[1:])
        np.minimum(H[:-1], H[1:] + G * DXY, out=H[:-1])
        np.minimum(H[:, 1:], H[:, :-1] + G * DXY, out=H[:, 1:])
        np.minimum(H[:, :-1], H[:, 1:] + G * DXY, out=H[:, :-1])
        if np.abs(H - b).max() < 1.0: break
    tri = float((L1 * L2) / np.mean(H ** 2) / (np.sqrt(3) / 4))
    print(f"[surf field] grid {ns}x{nt} @ {DXY:.0f} m   h {H.min():.1f}..{H.max():.1f} m "
          f"(med {np.median(H):.0f})")
    print(f"[surf field] gate-limited fraction: {100*np.mean(hg.reshape(ns,nt) <= np.minimum(hf.reshape(ns,nt),HMAX)):.1f} %"
          f"   ~{tri/1e6:.2f} M top triangles expected")
    return dict(H=H, ss=ss, tt=tt, u=u, v=v, W=W, L1=L1, L2=L2, ns=ns, nt=nt)

if __name__ == "__main__":
    Q = np.load("build_tmp/domain_corners_utm.npy")
    fs = np.load("build_tmp/fault_surface.npz")
    ch = np.load("build_tmp/trace_chains.npz")
    idx = np.unique(np.concatenate([ch[k] for k in ch.files]))
    tp = fs["P"][idx]
    steep = tp[np.abs(tp[:, 2]) > 5.0]
    print(f"[surf field] {len(steep)} steep trace nodes (|z|>5 m) of {len(tp):,}")
    out = build(Q, tp, steep_P=steep)
    np.savez_compressed("build_tmp/surface_field.npz", **{k: v for k, v in out.items()})
    print("[out] build_tmp/surface_field.npz")
