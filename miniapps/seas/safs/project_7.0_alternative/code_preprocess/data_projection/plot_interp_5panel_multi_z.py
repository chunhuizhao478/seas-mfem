"""plot_interp_5panel_multi_z.py — 5-panel continuous-field comparison
at multiple z slices.

Renders sidecar + 4 mesh GF slices as continuous interpolated fields
(NOT scatter dots) on a fine 2D xy grid:

  Panel 1  Sidecar:    trilinear of velocity_safs.h5 at (X, Y, z_target)
  Panel 2  OLD:        FE P1 GF, baseline mesh, evaluated at z=z_target
  Panel 3  NEW:        FE P1 GF, lc_far=5000  mesh
  Panel 4  NEW2:       FE P1 GF, lc_far=3000  mesh
  Panel 5  NEW3:       FE P1 GF, z-graded     mesh + dense sidecar

Each mesh panel is rendered by:
  (i) For every tet that straddles z=z_target, compute the FE-linear
      barycentric value at the tet's xy-centroid lifted to z=z_target.
  (ii) Feed that scattered cloud into LinearNDInterpolator over a 2D
       Delaunay and resample on a regular xy grid -> imshow.

Generates one PNG per requested depth in
  data_projected/comparison_5panel_interp_z<DEPTH>km.png
"""
from __future__ import annotations
import os, sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT); sys.path.insert(0, HERE); sys.path.insert(0, "/tmp")
from parse_mfem_vtu import parse_vtu                           # noqa: E402
import pyproj                                                  # noqa: E402
from plot_comparison_with_dg0 import (                         # noqa: E402
    load_stl_vertices, stl_surface_trace_xy, trilinear,
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

DEPTHS_KM   = [0, 1, 2, 3, 4, 5, 6]   # km below surface (positive numbers)
DISPLAY_DX_M = 1000.0                  # 2D display grid spacing (m)

TRANS = pyproj.Transformer.from_crs("EPSG:32611", "EPSG:4326",
                                     always_xy=True)


def utm_to_lonlat(x_m, y_m):
    return TRANS.transform(x_m, y_m)


def tet_slice_scatter(vtu_path, z_target):
    """Return (xc, yc, gf_at_slice) for every tet that straddles z=z_target.
    gf is FE-linear barycentric interp of corner DOFs at (xc, yc, z_target)."""
    m = parse_vtu(vtu_path)
    pts  = np.ascontiguousarray(m["points"], dtype=np.float64)
    conn = m["conn"].astype(np.int64)
    offs = m["offsets"].astype(np.int64)
    Vs_dof = m["point_data"]["Vs"]
    n_per = int(offs[0]); n_cells = offs.size
    corners_idx = conn.reshape(n_cells, n_per)[:, :4]
    V = pts[corners_idx]; F = Vs_dof[corners_idx]
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
    return xc, yc, gf


def scatter_to_grid(xc, yc, vals, X, Y):
    """LinearNDInterpolator on (xc, yc, vals) evaluated at (X, Y)."""
    P = np.column_stack([xc, yc])
    _, uidx = np.unique(P, axis=0, return_index=True)
    interp = LinearNDInterpolator(P[uidx], vals[uidx], fill_value=np.nan)
    return interp(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)


def build_grid(gx, gy):
    """Fine xy grid spanning the sidecar bbox at DISPLAY_DX_M spacing."""
    xs = np.arange(gx[0], gx[-1] + DISPLAY_DX_M, DISPLAY_DX_M)
    ys = np.arange(gy[0], gy[-1] + DISPLAY_DX_M, DISPLAY_DX_M)
    X, Y = np.meshgrid(xs, ys, indexing="xy")  # imshow-friendly
    return xs, ys, X, Y


def render_one_depth(z_target_m, gx, gy, gz, Vs_grid, X, Y,
                     extent_lon, extent_lat, flon, flat,
                     out_path):
    z_km_pos = -z_target_m / 1000.0
    print(f"\n=== rendering z = {z_target_m:.0f} m ===")

    # Panel 1: sidecar trilinear at (X, Y, z_target)
    print("  sidecar trilinear ...")
    qxyz = np.column_stack([X.ravel(), Y.ravel(),
                             np.full(X.size, z_target_m)])
    sidecar_field = trilinear(gx, gy, gz, Vs_grid, qxyz).reshape(X.shape)

    # Panels 2-5: mesh slices via scatter -> grid
    panels = []
    for label, vtu in [
        ("OLD baseline (lc_far=10000)", PREVIEW_OLD),
        ("NEW (lc_far=5000, dist_outer=80000)", PREVIEW_NEW),
        ("NEW2 (lc_far=3000, dist_outer=80000)", PREVIEW_NEW2),
        ("NEW3 z-graded + dense sidecar", PREVIEW_NEW3),
    ]:
        if not os.path.isfile(vtu):
            print(f"  SKIP missing: {vtu}")
            panels.append((label, None))
            continue
        print(f"  reading {label} ...")
        xc, yc, gf = tet_slice_scatter(vtu, z_target_m)
        print(f"    {gf.size} tets straddle z={z_target_m:.0f} m; gridding ...")
        if gf.size < 4:
            panels.append((label, np.full(X.shape, np.nan)))
            continue
        field = scatter_to_grid(xc, yc, gf, X, Y)
        panels.append((label, field))

    # Convert grid X, Y from UTM -> lon/lat for display
    LON, LAT = utm_to_lonlat(X, Y)

    # ---- Plot ----
    cmap_v = plt.get_cmap("jet")
    vmin, vmax = 0.0, 3.5
    fig = plt.figure(figsize=(30, 6.0))
    gs = fig.add_gridspec(1, 5, wspace=0.25,
                          left=0.04, right=0.99, top=0.84, bottom=0.13)

    def draw(ax, field, title):
        # imshow on lon/lat extent — note imshow assumes uniform grid
        im = ax.imshow(field / 1000.0, origin="lower",
                       extent=[extent_lon[0], extent_lon[1],
                                extent_lat[0], extent_lat[1]],
                       aspect="auto",
                       cmap=cmap_v, vmin=vmin, vmax=vmax,
                       interpolation="bilinear")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        if flon is not None:
            ax.scatter(flon, flat, s=0.06, c="k", alpha=0.6)
        fig.colorbar(im, ax=ax, fraction=0.046, label="Vs (km/s)")

    draw(fig.add_subplot(gs[0]), sidecar_field,
         f"Sidecar trilinear @ z={z_target_m:.0f} m\n(continuous truth)")
    for i, (label, field) in enumerate(panels, start=1):
        ax = fig.add_subplot(gs[i])
        if field is None:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes); continue
        draw(ax, field, f"{label}\n@ z={z_target_m:.0f} m (FE-interp)")

    fig.suptitle(
        f"Continuous-field comparison @ z = {z_target_m:.0f} m: sidecar "
        f"(trilinear truth) vs 4 mesh GFs (FE-linear interpolated)",
        fontsize=12)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    print(f"reading sidecar {SIDECAR}")
    with h5py.File(SIDECAR, "r") as f:
        gx = f["grid/x"][...].astype(np.float64)
        gy = f["grid/y"][...].astype(np.float64)
        gz = f["grid/z"][...].astype(np.float64)
        Vs = f["fields/Vs"][...].astype(np.float64)

    xs, ys, X, Y = build_grid(gx, gy)
    LON, LAT = utm_to_lonlat(X, Y)
    extent_lon = (LON.min(), LON.max())
    extent_lat = (LAT.min(), LAT.max())
    print(f"  display grid: {X.shape} ({xs.size} x {ys.size}) at {DISPLAY_DX_M} m")

    # Fault trace overlay
    stl_v = load_stl_vertices(STL)
    flon = flat = None
    for z_show in [-1000.0]:
        trace = stl_surface_trace_xy(stl_v, z_show)
        if trace is not None and len(trace):
            flon_t, flat_t = utm_to_lonlat(trace[:, 0], trace[:, 1])
            if len(flon_t) > 5000:
                s = len(flon_t) // 5000 + 1
                flon_t = flon_t[::s]; flat_t = flat_t[::s]
            flon, flat = flon_t, flat_t
            break

    for d_km_pos in DEPTHS_KM:
        z_target = -float(d_km_pos) * 1000.0
        out = (f"{ROOT}/data_projected/"
               f"comparison_5panel_interp_z{int(d_km_pos)}km.png")
        render_one_depth(z_target, gx, gy, gz, Vs, X, Y,
                          extent_lon, extent_lat, flon, flat, out)


if __name__ == "__main__":
    main()
