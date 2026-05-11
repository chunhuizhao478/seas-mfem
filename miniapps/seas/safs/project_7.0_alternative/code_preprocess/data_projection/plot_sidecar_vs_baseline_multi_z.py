"""plot_sidecar_vs_baseline_multi_z.py — 3×3 multi-z multi-field
comparison.  Generates one PNG per (field, z) pair:

For each field in {Vs, Vp, density}, for each depth z, the figure has
9 panels in a 3×3 grid:

    row 0:  sidecar truth | baseline value     | lc_far=5000 value
    row 1:  lc_far=3000 value | z-graded value | err baseline
    row 2:  err lc_far=5000   | err lc_far=3000 | err z-graded

All panels rendered as continuous interpolated fields (LinearNDInterp
on tet-slice scatter for mesh GFs; trilinear on the sidecar).  White =
outside the data hull / outside the sidecar bbox.

Error panel = (mesh GF) - (sidecar trilinear at the same xy, z).
"""
from __future__ import annotations
import os, sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator

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
SIDECAR = f"{ROOT}/data_projected/velocity_safs.h5"
PREVIEW_BASELINE = (f"{ROOT}/data_projected/preview/"
                     "projected_velocity_500m/Cycle000000/proc000000.vtu")
PREVIEW_LCFAR5K = (f"{ROOT}/data_projected/preview/"
                    "projected_velocity_500m_d80km/Cycle000000/proc000000.vtu")
PREVIEW_LCFAR3K = (f"{ROOT}/data_projected/preview/"
                    "projected_velocity_500m_lcfar3000/Cycle000000/proc000000.vtu")
PREVIEW_ZGRADED = (f"{ROOT}/data_projected/preview/"
                    "projected_velocity_500m_zgraded/Cycle000000/proc000000.vtu")
STL = (f"{ROOT}/data_cutnwfault/"
       "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m"
       "_clean_clip_nwcut.stl")

DEPTHS_KM = [0, 1, 2, 3, 4, 5, 7.5, 10, 15, 20]

# Per-field config: (sidecar dataset key, .vtu PointData key, vmax for
# value panel, abs-max for error panel, units string).
FIELD_CONFIG = {
    "Vs":      ("fields/Vs",      "Vs",      3500.0, 1000.0, "m/s"),
    "Vp":      ("fields/Vp",      "Vp",      8000.0, 2000.0, "m/s"),
    "density": ("fields/density", "density", 3300.0,  500.0, "kg/m^3"),
}

DISPLAY_DX_M = 1000.0     # imshow grid spacing

TRANS = pyproj.Transformer.from_crs("EPSG:32611", "EPSG:4326",
                                     always_xy=True)


def utm_to_lonlat(x_m, y_m):
    return TRANS.transform(x_m, y_m)


def tet_slice_scatter(vtu_path, z_target, field_key):
    """For every tet straddling z=z_target, return (xc, yc, gf) where
    gf is FE-linear barycentric interp of the requested PointData
    field at (xc, yc, z_target)."""
    m = parse_vtu(vtu_path)
    pts  = np.ascontiguousarray(m["points"], dtype=np.float64)
    conn = m["conn"].astype(np.int64)
    offs = m["offsets"].astype(np.int64)
    if field_key not in m["point_data"]:
        raise KeyError(
            f"{vtu_path} has no PointData '{field_key}'; available: "
            f"{list(m['point_data'].keys())}")
    fld = m["point_data"][field_key]
    n_per = int(offs[0]); n_cells = offs.size
    corners_idx = conn.reshape(n_cells, n_per)[:, :4]
    V = pts[corners_idx]; F = fld[corners_idx]
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
    P = np.column_stack([xc, yc])
    _, uidx = np.unique(P, axis=0, return_index=True)
    interp = LinearNDInterpolator(P[uidx], vals[uidx], fill_value=np.nan)
    return interp(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)


def render_one_field_z(field_name, h5_key, vtu_key,
                        vmax_val, vmax_err, units,
                        z_target_m,
                        gx, gy, gz, F_raw,
                        X, Y, extent_lon, extent_lat,
                        flon, flat,
                        panel_specs, out_path):
    print(f"\n--- field={field_name}  z={z_target_m:.0f} m ---")

    # Panel 0: sidecar trilinear at (X, Y, z_target)
    print("  sidecar trilinear ...")
    qxyz = np.column_stack([X.ravel(), Y.ravel(),
                             np.full(X.size, z_target_m)])
    sidecar_field = trilinear(gx, gy, gz, F_raw, qxyz).reshape(X.shape)

    # Panels 1..N: each mesh's tet-slice for the requested field
    mesh_fields = []
    for name, vtu in panel_specs:
        print(f"  {name} tet-slice + grid ...")
        try:
            xc, yc, gf = tet_slice_scatter(vtu, z_target_m, vtu_key)
        except KeyError as e:
            print(f"    SKIP: {e}")
            mesh_fields.append((name, None, 0))
            continue
        print(f"    {gf.size} tets straddle; gridding ...")
        if gf.size < 4:
            mesh_fields.append((name, np.full(X.shape, np.nan), 0))
            continue
        field = scatter_to_grid(xc, yc, gf, X, Y)
        mesh_fields.append((name, field, gf.size))

    # ---- Plot ----
    cmap_v = plt.get_cmap("jet")
    cmap_v.set_bad(color="white")
    cmap_e = plt.get_cmap("RdBu_r")
    cmap_e.set_bad(color="white")
    vmin = 0.0

    fig = plt.figure(figsize=(18, 16))
    gs = fig.add_gridspec(3, 3, wspace=0.20, hspace=0.30,
                          left=0.04, right=0.99,
                          top=0.92, bottom=0.05)

    def draw_value(ax, field, title):
        im = ax.imshow(field, origin="lower",
                       extent=[extent_lon[0], extent_lon[1],
                                extent_lat[0], extent_lat[1]],
                       aspect="auto", cmap=cmap_v,
                       vmin=vmin, vmax=vmax_val,
                       interpolation="bilinear")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Longitude (deg)"); ax.set_ylabel("Latitude (deg)")
        if flon is not None:
            ax.scatter(flon, flat, s=0.06, c="k", alpha=0.6)
        fig.colorbar(im, ax=ax, fraction=0.046,
                     label=f"{field_name} ({units})")

    def draw_error(ax, mesh_field, title):
        if mesh_field is None:
            ax.text(0.5, 0.5, "(missing)", ha="center", va="center",
                    transform=ax.transAxes); return
        diff = mesh_field - sidecar_field
        im = ax.imshow(diff, origin="lower",
                       extent=[extent_lon[0], extent_lon[1],
                                extent_lat[0], extent_lat[1]],
                       aspect="auto", cmap=cmap_e,
                       vmin=-vmax_err, vmax=vmax_err,
                       interpolation="bilinear")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Longitude (deg)"); ax.set_ylabel("Latitude (deg)")
        if flon is not None:
            ax.scatter(flon, flat, s=0.06, c="k", alpha=0.6)
        fig.colorbar(im, ax=ax, fraction=0.046,
                     label=f"Δ{field_name} ({units})")

    # Row 0: sidecar | baseline | lc_far=5000
    draw_value(fig.add_subplot(gs[0, 0]), sidecar_field,
               f"Sidecar trilinear (truth) {field_name} @ z={z_target_m:.0f} m")
    if len(mesh_fields) >= 1:
        n, f, k = mesh_fields[0]
        draw_value(fig.add_subplot(gs[0, 1]), f,
                   f"{n}\n{field_name} @ z={z_target_m:.0f} m  ({k} tets)")
    if len(mesh_fields) >= 2:
        n, f, k = mesh_fields[1]
        draw_value(fig.add_subplot(gs[0, 2]), f,
                   f"{n}\n{field_name} @ z={z_target_m:.0f} m  ({k} tets)")

    # Row 1: lc_far=3000 | z-graded | err baseline
    if len(mesh_fields) >= 3:
        n, f, k = mesh_fields[2]
        draw_value(fig.add_subplot(gs[1, 0]), f,
                   f"{n}\n{field_name} @ z={z_target_m:.0f} m  ({k} tets)")
    if len(mesh_fields) >= 4:
        n, f, k = mesh_fields[3]
        draw_value(fig.add_subplot(gs[1, 1]), f,
                   f"{n}\n{field_name} @ z={z_target_m:.0f} m  ({k} tets)")
    if len(mesh_fields) >= 1:
        n, f, _ = mesh_fields[0]
        draw_error(fig.add_subplot(gs[1, 2]), f,
                   f"err: {n}\n(GF - sidecar)")

    # Row 2: err lc=5000 | err lc=3000 | err z-graded
    if len(mesh_fields) >= 2:
        n, f, _ = mesh_fields[1]
        draw_error(fig.add_subplot(gs[2, 0]), f, f"err: {n}\n(GF - sidecar)")
    if len(mesh_fields) >= 3:
        n, f, _ = mesh_fields[2]
        draw_error(fig.add_subplot(gs[2, 1]), f, f"err: {n}\n(GF - sidecar)")
    if len(mesh_fields) >= 4:
        n, f, _ = mesh_fields[3]
        draw_error(fig.add_subplot(gs[2, 2]), f, f"err: {n}\n(GF - sidecar)")

    fig.suptitle(
        f"{field_name} @ z = {z_target_m:.0f} m  —  "
        f"5 value panels + 4 error panels (mesh - sidecar) "
        f"on a continuous 1 km display grid",
        fontsize=13)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    print(f"reading sidecar {SIDECAR}")
    with h5py.File(SIDECAR, "r") as f:
        gx = f["grid/x"][...].astype(np.float64)
        gy = f["grid/y"][...].astype(np.float64)
        gz = f["grid/z"][...].astype(np.float64)
        sidecar_fields = {fname: f[h5key][...].astype(np.float64)
                          for fname, (h5key, *_) in FIELD_CONFIG.items()}

    xs = np.arange(gx[0], gx[-1] + DISPLAY_DX_M, DISPLAY_DX_M)
    ys = np.arange(gy[0], gy[-1] + DISPLAY_DX_M, DISPLAY_DX_M)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    LON, LAT = utm_to_lonlat(X, Y)
    extent_lon = (LON.min(), LON.max())
    extent_lat = (LAT.min(), LAT.max())
    print(f"  display grid: {X.shape}")

    panel_specs = [
        ("Baseline H1-P1 (lc_far=10000)", PREVIEW_BASELINE),
        ("lc_far=5000",                   PREVIEW_LCFAR5K),
        ("lc_far=3000",                   PREVIEW_LCFAR3K),
        ("z-graded + dense sidecar",      PREVIEW_ZGRADED),
    ]
    panel_specs = [(name, p) for name, p in panel_specs if os.path.isfile(p)]
    print(f"  mesh panels: {[n for n, _ in panel_specs]}")

    stl_v = load_stl_vertices(STL)

    for d_km_pos in DEPTHS_KM:
        z_target = -float(d_km_pos) * 1000.0
        # fault trace at z_target (or fallback)
        trace = stl_surface_trace_xy(stl_v, z_target)
        if trace is None or len(trace) == 0:
            trace = stl_surface_trace_xy(stl_v, -1000.0)
        flon = flat = None
        if trace is not None and len(trace):
            flon, flat = utm_to_lonlat(trace[:, 0], trace[:, 1])
            if len(flon) > 5000:
                s = len(flon) // 5000 + 1
                flon = flon[::s]; flat = flat[::s]

        depth_str = (f"{int(d_km_pos)}"
                     if d_km_pos == int(d_km_pos)
                     else f"{d_km_pos}".replace(".", "p"))

        for fname, (h5key, vtu_key, vmax_val, vmax_err, units) in \
                FIELD_CONFIG.items():
            out = (f"{ROOT}/data_projected/"
                   f"compare_3x3_{fname}_z{depth_str}km.png")
            render_one_field_z(fname, h5key, vtu_key,
                                vmax_val, vmax_err, units,
                                z_target, gx, gy, gz,
                                sidecar_fields[fname],
                                X, Y, extent_lon, extent_lat,
                                flon, flat,
                                panel_specs, out)


if __name__ == "__main__":
    main()
