"""plot_old_vs_new_baseline.py — focused 3-panel verification.

Compare, at z = -1000 m:
  Panel 1 — sidecar raw F[i, j, k=27] (truth, on the 1500 m grid).
  Panel 2 — OLD H1-P1 baseline mesh slice
            (lc_far=10000, dist_outer=40000).
  Panel 3 — NEW H1-P1 baseline mesh slice
            (lc_far=5000, dist_outer=80000) — refined ramp.

The "slice" is built by, for every tet that straddles z=-1000 m,
sampling the FE-linear value at the tet's xy-centroid lifted to
z=-1000 m.  Linear barycentric on the 4 corner DOFs.

The stdout summary reports per-distance-band residuals (|GF slice -
sidecar truth|) so the OLD vs NEW improvement is quantitative.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)
sys.path.insert(0, HERE)
sys.path.insert(0, "/tmp")
from parse_mfem_vtu import parse_vtu                                # noqa: E402
import pyproj                                                       # noqa: E402

from plot_comparison_with_dg0 import (                              # noqa: E402
    load_stl_vertices,
    stl_surface_trace_xy,
    trilinear,
)


ROOT = ("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/"
        "project_7.0_alternative")
SIDECAR    = f"{ROOT}/data_projected/velocity_safs.h5"
PREVIEW_OLD = (f"{ROOT}/data_projected/preview/"
               "projected_velocity_500m/Cycle000000/proc000000.vtu")
PREVIEW_NEW = (f"{ROOT}/data_projected/preview/"
               "projected_velocity_500m_d80km/Cycle000000/proc000000.vtu")
PREVIEW_NEW2 = (f"{ROOT}/data_projected/preview/"
                "projected_velocity_500m_lcfar3000/Cycle000000/proc000000.vtu")
PREVIEW_NEW3 = (f"{ROOT}/data_projected/preview/"
                "projected_velocity_500m_zgraded/Cycle000000/proc000000.vtu")
STL = (f"{ROOT}/data_cutnwfault/"
       "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m"
       "_clean_clip_nwcut.stl")
OUT = f"{ROOT}/data_projected/comparison_old_vs_new_baseline_z1km.png"


TRANS = pyproj.Transformer.from_crs("EPSG:32611", "EPSG:4326",
                                     always_xy=True)


def utm_to_lonlat(x_m, y_m):
    return TRANS.transform(x_m, y_m)


def tet_slice(vtu_path, gx, gy, gz, Vs, z_target=-1000.0):
    """Read .vtu, return (xc, yc, gf_at_slice, truth_at_slice, n_cells)
    for every tet that straddles z=z_target.  gf_at_slice is the FE-
    linear barycentric value at (xc, yc, z_target).  truth_at_slice
    is sidecar trilinear at the same point."""
    m = parse_vtu(vtu_path)
    pts = np.ascontiguousarray(m["points"], dtype=np.float64)
    conn = m["conn"].astype(np.int64)
    offs = m["offsets"].astype(np.int64)
    Vs_dof = m["point_data"]["Vs"]
    n_per = int(offs[0])
    n_cells = offs.size
    corners_idx = conn.reshape(n_cells, n_per)[:, :4]
    V = pts[corners_idx]
    F = Vs_dof[corners_idx]
    z = V[:, :, 2]
    sel = (z.min(axis=1) <= z_target) & (z.max(axis=1) >= z_target)
    Vt = V[sel]; Ft = F[sel]
    xc = Vt[:, :, 0].mean(axis=1)
    yc = Vt[:, :, 1].mean(axis=1)
    qpts = np.stack([xc, yc, np.full(xc.size, z_target)], axis=1)
    Amat = np.stack([Vt[:, 1] - Vt[:, 0],
                     Vt[:, 2] - Vt[:, 0],
                     Vt[:, 3] - Vt[:, 0]], axis=2)
    rhs = (qpts - Vt[:, 0])[:, :, None]
    bc = np.linalg.solve(Amat, rhs).squeeze(-1)
    b0 = 1.0 - bc.sum(axis=1)
    gf = (b0 * Ft[:, 0] + bc[:, 0] * Ft[:, 1]
          + bc[:, 1] * Ft[:, 2] + bc[:, 2] * Ft[:, 3])
    truth = trilinear(gx, gy, gz, Vs, qpts)
    z_height = z[sel].max(axis=1) - z[sel].min(axis=1)
    return xc, yc, gf, truth, n_cells, z_height


def print_residual_table(name, xc, yc, gf, truth, z_height, trace):
    from scipy.spatial import cKDTree
    tree = cKDTree(trace[:, :2])
    d, _ = tree.query(np.column_stack([xc, yc]), k=1)
    diff = gf - truth
    print(f"\n=== {name}  ({gf.size} tets straddle z=-1000 m) ===")
    print(f"  {'dist band (m)':<15} {'n_tets':>7}  "
          f"{'mean|d|':>9} {'med|d|':>9} {'p95|d|':>9} "
          f"{'max|d|':>9}  {'tet ht p95':>10}")
    for lo, hi in [(0, 500), (500, 1500), (1500, 3000),
                   (3000, 10000), (10000, 30000),
                   (30000, 80000), (80000, 1e18)]:
        mask = (d >= lo) & (d < hi)
        if mask.sum() == 0: continue
        ad = np.abs(diff[mask]); zh = z_height[mask]
        hi_s = "inf" if hi > 1e9 else f"{int(hi)}"
        print(f"  [{int(lo):>5}, {hi_s:>5}):  {mask.sum():>7}  "
              f"{ad.mean():>9.1f} {np.median(ad):>9.1f} "
              f"{np.percentile(ad,95):>9.1f} {ad.max():>9.1f}  "
              f"{np.percentile(zh,95):>10.0f}")
    print(f"  {'GLOBAL':<15} {gf.size:>7}  "
          f"{np.mean(np.abs(diff)):>9.1f} {np.median(np.abs(diff)):>9.1f} "
          f"{np.percentile(np.abs(diff),95):>9.1f} "
          f"{np.max(np.abs(diff)):>9.1f}")


def main():
    z_target = -1000.0
    print(f"reading sidecar {SIDECAR}")
    with h5py.File(SIDECAR, "r") as f:
        gx_s = f["grid/x"][...]; gy_s = f["grid/y"][...]
        gz_s = f["grid/z"][...]; Vs = f["fields/Vs"][...]
    k_sc = int(np.argmin(np.abs(gz_s - z_target)))
    print(f"  sidecar slice k={k_sc}, gz[k]={gz_s[k_sc]:.1f} m")
    II, JJ = np.meshgrid(np.arange(gx_s.size), np.arange(gy_s.size),
                          indexing="ij")
    II = II.ravel(); JJ = JJ.ravel()
    sx_m = gx_s[II]; sy_m = gy_s[JJ]; sv = Vs[II, JJ, k_sc]

    print(f"reading OLD vtu: {PREVIEW_OLD}")
    oxc, oyc, ogf, otrue, on_c, oz_h = tet_slice(
        PREVIEW_OLD, gx_s, gy_s, gz_s, Vs)
    print(f"  OLD: {on_c} cells, {ogf.size} straddle z=-1000 m")

    print(f"reading NEW vtu: {PREVIEW_NEW}")
    if not os.path.isfile(PREVIEW_NEW):
        print(f"  !!! {PREVIEW_NEW} not found — projection not yet run")
        return
    nxc, nyc, ngf, ntrue, nn_c, nz_h = tet_slice(
        PREVIEW_NEW, gx_s, gy_s, gz_s, Vs)
    print(f"  NEW: {nn_c} cells, {ngf.size} straddle z=-1000 m")

    has_new2 = os.path.isfile(PREVIEW_NEW2)
    if has_new2:
        print(f"reading NEW2 vtu: {PREVIEW_NEW2}")
        n2xc, n2yc, n2gf, n2true, n2_c, n2z_h = tet_slice(
            PREVIEW_NEW2, gx_s, gy_s, gz_s, Vs)
        print(f"  NEW2: {n2_c} cells, {n2gf.size} straddle z=-1000 m")
    else:
        print(f"  (NEW2 {PREVIEW_NEW2} not yet projected — skipping panel)")

    has_new3 = os.path.isfile(PREVIEW_NEW3)
    if has_new3:
        print(f"reading NEW3 (z-graded) vtu: {PREVIEW_NEW3}")
        n3xc, n3yc, n3gf, n3true, n3_c, n3z_h = tet_slice(
            PREVIEW_NEW3, gx_s, gy_s, gz_s, Vs, z_target=z_target)
        print(f"  NEW3: {n3_c} cells, {n3gf.size} straddle z={z_target:.0f} m")
    else:
        print(f"  (NEW3 {PREVIEW_NEW3} not yet projected — skipping panel)")

    # Fault trace overlay
    stl_v = load_stl_vertices(STL)
    trace = stl_surface_trace_xy(stl_v, z_target)
    flon, flat = utm_to_lonlat(trace[:, 0], trace[:, 1])
    if len(flon) > 5000:
        s = len(flon) // 5000 + 1
        flon = flon[::s]; flat = flat[::s]

    # Residual stats
    print_residual_table("OLD (lc_far=10000, dist_outer=40000)",
                          oxc, oyc, ogf, otrue, oz_h, trace)
    print_residual_table("NEW (lc_far=5000, dist_outer=80000)",
                          nxc, nyc, ngf, ntrue, nz_h, trace)
    if has_new2:
        print_residual_table("NEW2 (lc_far=3000, dist_outer=80000)",
                              n2xc, n2yc, n2gf, n2true, n2z_h, trace)
    if has_new3:
        print_residual_table("NEW3 (z-graded + dense sidecar)",
                              n3xc, n3yc, n3gf, n3true, n3z_h, trace)

    # ---- Plot ---------------------------------------------------------
    s_lon, s_lat = utm_to_lonlat(sx_m, sy_m)
    o_lon, o_lat = utm_to_lonlat(oxc, oyc)
    n_lon, n_lat = utm_to_lonlat(nxc, nyc)
    if has_new2:
        n2_lon, n2_lat = utm_to_lonlat(n2xc, n2yc)
    if has_new3:
        n3_lon, n3_lat = utm_to_lonlat(n3xc, n3yc)
    extent = [s_lon.min(), s_lon.max(), s_lat.min(), s_lat.max()]

    cmap_v = plt.get_cmap("jet")
    vmin, vmax = 0.0, 3.5

    ncol = 3 + int(has_new2) + int(has_new3)
    fig = plt.figure(figsize=(6 * ncol, 6.5))
    gs = fig.add_gridspec(1, ncol, wspace=0.25,
                          left=0.04, right=0.99,
                          top=0.83, bottom=0.13)

    def setup(ax, title):
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("auto")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        ax.scatter(flon, flat, s=0.06, c="k", alpha=0.5)

    ax0 = fig.add_subplot(gs[0])
    s0 = ax0.scatter(s_lon, s_lat, s=2, c=sv / 1000.0,
                     cmap=cmap_v, vmin=vmin, vmax=vmax,
                     edgecolors="none")
    setup(ax0, f"Sidecar raw F[i,j,k={k_sc}] @ z={gz_s[k_sc]:.0f} m\n"
               f"{sv.size} pts on 1500 m grid, NO interpolation")
    fig.colorbar(s0, ax=ax0, fraction=0.046, label="Vs (km/s)")

    ax1 = fig.add_subplot(gs[1])
    s1 = ax1.scatter(o_lon, o_lat, s=2, c=ogf / 1000.0,
                     cmap=cmap_v, vmin=vmin, vmax=vmax,
                     edgecolors="none")
    setup(ax1, f"OLD baseline H1-P1 slice @ z=-1000 m\n"
               f"lc_far=10000, dist_outer=40000  "
               f"({on_c} cells, {ogf.size} straddle)")
    fig.colorbar(s1, ax=ax1, fraction=0.046, label="Vs (km/s)")

    ax2 = fig.add_subplot(gs[2])
    s2 = ax2.scatter(n_lon, n_lat, s=2, c=ngf / 1000.0,
                     cmap=cmap_v, vmin=vmin, vmax=vmax,
                     edgecolors="none")
    setup(ax2, f"NEW baseline H1-P1 slice @ z=-1000 m\n"
               f"lc_far=5000, dist_outer=80000  "
               f"({nn_c} cells, {ngf.size} straddle)")
    fig.colorbar(s2, ax=ax2, fraction=0.046, label="Vs (km/s)")

    if has_new2:
        ax3 = fig.add_subplot(gs[3])
        s3 = ax3.scatter(n2_lon, n2_lat, s=2, c=n2gf / 1000.0,
                         cmap=cmap_v, vmin=vmin, vmax=vmax,
                         edgecolors="none")
        setup(ax3, f"NEW2 baseline H1-P1 slice @ z=-1000 m\n"
                   f"lc_far=3000, dist_outer=80000  "
                   f"({n2_c} cells, {n2gf.size} straddle)")
        fig.colorbar(s3, ax=ax3, fraction=0.046, label="Vs (km/s)")

    if has_new3:
        col = 3 + int(has_new2)
        ax4 = fig.add_subplot(gs[col])
        s4 = ax4.scatter(n3_lon, n3_lat, s=2, c=n3gf / 1000.0,
                         cmap=cmap_v, vmin=vmin, vmax=vmax,
                         edgecolors="none")
        setup(ax4,
              f"NEW3 z-graded + dense sidecar @ z={z_target:.0f} m\n"
              f"top500m=500, [-3km,-500m]=1000, deeper=2-3km  "
              f"({n3_c} cells, {n3gf.size} straddle)")
        fig.colorbar(s4, ax=ax4, fraction=0.046, label="Vs (km/s)")

    fig.suptitle(
        f"OLD vs NEW baseline H1-P1: does refined ramp + lower lc_far "
        f"cut the 1-2.7 km/s far-field interpolation errors?",
        fontsize=11)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
