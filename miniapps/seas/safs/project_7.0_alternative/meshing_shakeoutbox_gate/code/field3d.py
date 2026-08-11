#!/usr/bin/env python3
"""field3d.py -- the gate-feasible 3-D size field, and what it costs.

THE DESIGN RULE, and why it is safe

The gate is not monotone in h (see size_field_price.py): shrinking a surface
cell drags its barycentre into the CVM's slow z = 0 bin and tightens its own
constraint, so the feasible set of a slow column can be two disjoint intervals
with a FORBIDDEN ZONE between them.  A mesher's gradient limiter only ever
REDUCES sizes, so a field that parks slow columns on the upper (coarse) branch
gets smeared straight through the forbidden zone at the branch boundary.

So the field is defined as the CONTIGUOUS feasible ceiling from h_min upward:

    H(x,y,z) = sup { h : every h' in [h_min, h] satisfies h' <= Vs(z - beta h')/gate }

which is downward-closed by construction.  Any h <= H is feasible, so gradient
limiting -- and any undershoot by the mesher -- keeps the cell feasible.  That
is the property the previous gate-driven collar lacked.

beta is the barycentre-depth / max-edge ratio, MEASURED on the shipped collar
by census_box.py (median 0.294, p10 0.124, p1 0.083), not assumed.

Outputs the field as npz on the CVM lattice, plus an integrated cell count.
"""
import argparse
from pathlib import Path

import h5py
import numpy as np

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


def feasible_ceiling(V, z, gate, beta, hs):
    """H(z,y,x): contiguous feasible ceiling from hs[0] upward.

    V is (nz, ny, nx) on the DESCENDING z axis `z`.  For a cell whose top sits
    at level k, size h, the barycentre is z[k] - beta*h and the gate wants
    Vs(that depth) >= gate*h.
    """
    nz, ny, nx = V.shape
    H = np.full((nz, ny, nx), np.nan, np.float32)
    run = np.ones((nz, ny, nx), bool)
    for h in hs:
        zb = z[:, None] - beta * h                       # (nz,1) barycentre depth
        k = np.abs(z[None, :] - zb).argmin(1)            # nearest level for each top level
        ok = V[k] >= gate * h                            # (nz,ny,nx)
        run &= ok
        H[run] = h
    return H


def lipschitz_limit(H, dx, dy, dz, slope, iters=12):
    """H <- min over neighbours of (H_nb + slope*distance); 6-neighbour sweeps.

    This is the discrete form of |grad H| <= slope, i.e. a gmsh/mmg `hgrad` of
    (1 + slope).  Sweeping alternately in +/- along each axis converges in a few
    passes on a field this smooth.

    `dz` may be a scalar or a per-gap array (MUSCAL's depth axis is 50 m near
    the surface and 1,000 m deep, so a single spacing would be wrong by 20x).
    """
    H = H.copy()
    dzc = np.asarray(dz, float)
    if dzc.ndim:
        dzc = dzc.reshape(-1, 1, 1)
    for _ in range(iters):
        before = H.sum()
        H[1:, :, :] = np.minimum(H[1:, :, :], H[:-1, :, :] + slope * dzc)
        H[:-1, :, :] = np.minimum(H[:-1, :, :], H[1:, :, :] + slope * dzc)
        H[:, 1:, :] = np.minimum(H[:, 1:, :], H[:, :-1, :] + slope * dy)
        H[:, :-1, :] = np.minimum(H[:, :-1, :], H[:, 1:, :] + slope * dy)
        H[:, :, 1:] = np.minimum(H[:, :, 1:], H[:, :, :-1] + slope * dx)
        H[:, :, :-1] = np.minimum(H[:, :, :-1], H[:, :, 1:] + slope * dx)
        if abs(H.sum() - before) < 1e-6 * abs(before):
            break
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cvm", default=CVM)
    ap.add_argument("--gate", type=float, default=0.8)
    ap.add_argument("--beta", type=float, default=0.294)
    ap.add_argument("--h-min", type=float, default=150.0)
    ap.add_argument("--h-max", type=float, default=4000.0)
    ap.add_argument("--hgrad", type=float, default=1.3)
    ap.add_argument("--step", type=float, default=25.0)
    ap.add_argument("--muscal", action="store_true",
                    help="build the field from MUSCAL (source model, 50 m depth "
                         "step near the surface) instead of the deck's resampled "
                         "250 m nc")
    ap.add_argument("--kfac", type=float, default=3.66,
                    help="measured cells-per-requested-size factor for this pipeline")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    with h5py.File(a.cvm, "r") as f:
        cx, cy, cz = f["x"][:], f["y"][:], f["z"][:]
        d = f["data"]
        VS = np.sqrt(np.maximum(np.asarray(d["mu"]).astype(np.float64), 0.0)
                     / np.asarray(d["rho"]).astype(np.float64)).astype(np.float32)

    # The domain box reaches E 786 km and S 3,524 km, but the CVM grid stops at
    # E 705 km / N 3,543 km -- roughly a fifth of the collar lies OUTSIDE the
    # CVM hull.  Build the field on a lattice covering the whole box and sample
    # Vs with nearest-node clamping, which is exactly what VsGrid (and ASAGI at
    # runtime) does there.  Integrating over the CVM's own footprint instead
    # silently drops 21 % of the collar volume.
    gx = np.arange(TARGET_E[0], TARGET_E[1] + 1.0, cx[1] - cx[0])
    gy = np.arange(TARGET_N[0], TARGET_N[1] + 1.0, cy[1] - cy[0])
    kz = (cz >= Z_BOT) & (cz <= 0.0)
    z = cz[kz]
    o = np.argsort(-z)
    z = z[o]

    def near(axis, q):
        i = np.clip(np.searchsorted(axis, q), 0, len(axis) - 1)
        j = np.clip(i - 1, 0, len(axis) - 1)
        return np.where(np.abs(axis[j] - q) <= np.abs(axis[i] - q), j, i)

    X, Y = np.meshgrid(gx, gy, indexing="xy")
    collar = ~inside_rect(X, Y, CORNERS)            # the lattice IS the box

    if a.muscal:
        # MUSCAL is the SOURCE model: 0.01 deg laterally and a 50 m depth step
        # through the top 500 m.  Use its own vertical axis -- resampling it
        # onto the deck file's 250 m lattice would reintroduce exactly the
        # binning artefact this is meant to remove.
        import sys as _s
        _s.path.insert(0, str(Path(__file__).resolve().parent))
        _s.path.insert(0, str(Path(__file__).resolve().parents[2]
                              / "meshing_shakeoutbox_intermediate" / "code"))
        from muscal_vs import MuscalVs
        from collar_lib import VsGrid as _VG
        mus = MuscalVs(backup=_VG(a.cvm))
        z = -mus.dep[(mus.dep >= 0.0) & (mus.dep <= -Z_BOT)]     # descending, <= 0
        V = np.empty((len(z), len(gy), len(gx)), np.float32)
        PT = np.empty((len(gy) * len(gx), 3))
        PT[:, 0] = X.ravel()
        PT[:, 1] = Y.ravel()
        for k, zz in enumerate(z):
            PT[:, 2] = zz
            V[k] = mus.at(PT).reshape(len(gy), len(gx))
        mus.report()
    else:
        ix, iy = near(cx, gx), near(cy, gy)
        V = VS[kz][o][:, iy[:, None], ix[None, :]]  # (nz, ny, nx), z descending

    hs = np.arange(a.h_min, a.h_max + 1.0, a.step)
    print(f"lattice {V.shape} (z,y,x)   z {z[0]:,.0f} .. {z[-1]:,.0f}")
    print(f"gate {a.gate}  beta {a.beta}  h in [{a.h_min:,.0f},{a.h_max:,.0f}]  "
          f"hgrad {a.hgrad}", flush=True)

    H = feasible_ceiling(V, z, a.gate, a.beta, hs)
    nan = ~np.isfinite(H)
    if nan.any():
        print(f"  WARNING: {int(nan[:, collar].sum()):,} collar lattice nodes have NO "
              f"feasible h >= {a.h_min:,.0f} m; clamped to h_min")
        H[nan] = a.h_min

    Hc = H[:, collar]
    print("\n  raw feasible ceiling over the collar:")
    probes = [("z=0", 0.0), ("z=-100", -100.0), ("z=-250", -250.0),
              ("z=-500", -500.0), ("z=-1000", -1000.0), ("z=-2000", -2000.0),
              ("z=-5000", -5000.0)]
    for lab, zq in probes:
        k = int(np.abs(z - zq).argmin())
        if True:
            v = Hc[k]
            print(f"    {lab:<9} p1 {np.percentile(v,1):>7,.0f}  p10 {np.percentile(v,10):>7,.0f}"
                  f"  med {np.median(v):>7,.0f}  p90 {np.percentile(v,90):>7,.0f}")

    gaps = np.abs(np.diff(z))                    # gap k is between level k and k+1
    HL = lipschitz_limit(H, gx[1] - gx[0], gy[1] - gy[0], gaps, a.hgrad - 1.0)
    HLc = HL[:, collar]
    print(f"\n  after Lipschitz limiting (slope {a.hgrad-1:.2f}):")
    probes = [("z=0", 0.0), ("z=-100", -100.0), ("z=-250", -250.0),
              ("z=-500", -500.0), ("z=-1000", -1000.0), ("z=-2000", -2000.0),
              ("z=-5000", -5000.0)]
    for lab, zq in probes:
        k = int(np.abs(z - zq).argmin())
        if True:
            v = HLc[k]
            print(f"    {lab:<9} p1 {np.percentile(v,1):>7,.0f}  p10 {np.percentile(v,10):>7,.0f}"
                  f"  med {np.median(v):>7,.0f}  p90 {np.percentile(v,90):>7,.0f}")

    # cells = K * integral of 6/h^3 dV over the collar, on the lattice.
    # K is MEASURED on the shipped collar by calibrate_density.py: a uniform
    # 2,500 m request produced 13,405,498 cells in 9,535,430 km3, i.e. 3.66x
    # the ideal-lattice count for the REQUESTED size.
    # per-level slab thickness (MUSCAL's axis is not uniform)
    edges = np.concatenate(([z[0]], 0.5 * (z[1:] + z[:-1]), [z[-1]]))
    thk = np.abs(np.diff(edges))
    dvk = ((gx[1] - gx[0]) * (gy[1] - gy[0]) * thk).reshape(-1, 1)
    dv = float(np.mean(dvk))
    vol = float(dvk.sum()) * int(collar.sum())
    n_raw = float((6.0 * dvk / Hc.astype(np.float64) ** 3).sum())
    n_lim = float((6.0 * dvk / HLc.astype(np.float64) ** 3).sum())
    print(f"\n  collar volume on this lattice {vol/1e9:,.0f} km3 "
          f"(shipped collar measured 9,535,430 km3)")
    print(f"  ANALYTIC collar cells   raw field {n_raw:,.0f}   limited {n_lim:,.0f}")
    print(f"  x{a.kfac} measured build factor -> raw {a.kfac*n_raw:,.0f}   "
          f"limited {a.kfac*n_lim:,.0f}")
    print(f"  shipped uniform 2,500 m collar, measured: 13,405,498")

    # where the cost sits
    print(f"\n  cost by depth band (limited field, x{a.kfac}):")
    cuts = [0.0, -125.0, -500.0, -1000.0, -2000.0, -5000.0, Z_BOT - 1.0]
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        k = np.nonzero((z <= lo) & (z > hi))[0]
        if not len(k):
            continue
        band = (6.0 * dvk[k] / HLc[k].astype(np.float64) ** 3).sum() * a.kfac
        print(f"    z {lo:>8,.0f} .. {hi:>9,.0f} m ({len(k):>3} levels): "
              f"{band:>14,.0f} cells")

    if a.out:
        np.savez_compressed(a.out, H=HL.astype(np.float32), x=gx, y=gy, z=z,
                            collar=collar, gate=a.gate, beta=a.beta, hgrad=a.hgrad)
        print(f"\n  field -> {a.out}")


if __name__ == "__main__":
    main()
