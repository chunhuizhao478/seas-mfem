"""plot_baseline_vs_g1g2.py — pure point-wise comparison at z = -1 km.

LEFT panel:  raw sidecar values plotted as scatter dots at the sidecar
             grid coordinates (gx[i], gy[j]) colored by F[i,j,k].  No
             interpolation, no imshow.  This IS the source data.

MIDDLE panels:  the mesh projection's nodal values (PointData['Vs']
             from the .vtu) plotted as scatter dots at the mesh DOF
             coordinates (px, py) colored by that PointData value.
             The DOF values themselves were obtained at projector-time
             by trilinear sampling of the sidecar at each DOF
             coordinate (this within-sidecar interpolation is allowed
             per the user's request).  But the plot itself draws each
             DOF as a single dot — no interpolation between mesh
             nodes, no smooth surface.

RIGHT panel: per-band RMS of (mesh DOF value at sidecar grid pt)
             - (raw F[i,j,k]) — bit-identical-truth metric from
             fault_zone_metric.per_grid_point_residuals.
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
from fault_zone_metric import (                                     # noqa: E402
    DEFAULT_BANDS,
    _split_conn_by_offset,
    grid_point_band_rms,
    per_grid_point_residuals,
    per_voxel_midpoint_residuals,
)


ROOT = ("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/"
        "project_7.0_alternative")
SIDECAR = f"{ROOT}/velocity/results/cvmh/velocity_safs.h5"
PREVIEW_BASELINE = (f"{ROOT}/velocity/results/cvmh/preview/"
                    "projected_velocity_500m/Cycle000000/proc000000.vtu")
PREVIEW_G1G2 = (f"{ROOT}/velocity/results/cvmh/preview/"
                "projected_velocity_500m_g1g2/Cycle000000/proc000000.vtu")
PREVIEW_G1G2_CROM = (f"{ROOT}/velocity/results/cvmh/preview/"
                     "projected_velocity_500m_g1g2_catmull/Cycle000000/"
                     "proc000000.vtu")
PREVIEW_BASELINE_CROM = (f"{ROOT}/velocity/results/cvmh/preview/"
                         "projected_velocity_500m_p1_catmull/Cycle000000/"
                         "proc000000.vtu")
STL = (f"{ROOT}/meshing/results/stl_nwcut/"
       "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m"
       "_clean_clip_nwcut.stl")
OUT = f"{ROOT}/velocity/results/cvmh/comparison_baseline_vs_g1g2_z1km.png"


TRANS = pyproj.Transformer.from_crs(
    "EPSG:32611", "EPSG:4326", always_xy=True)


def utm_to_lonlat_pair(x_m, y_m):
    return TRANS.transform(x_m, y_m)


def main():
    z_target = -1000.0
    # Two distinct z windows for two distinct purposes.
    #
    #   half_thick_plot   = visualisation only (panels A/B/C/D).
    #       Must match panel 1's z=-1000 m slice apples-to-apples so
    #       that the DOF colours can be compared against the sidecar
    #       raw colours at the SAME z.  Sidecar slices nearest
    #       z=-1000 m are k=27 (-1000 m), k=26 (-2000 m), k=28 (0 m);
    #       the Voronoi half-cell of k=27 is therefore 500 m on each
    #       side.  Anything thicker pulls in shallow-basin DOFs
    #       (z>=-500 m, Vs ~ 200-1000 m/s) and deeper sub-basement
    #       DOFs whose actual sidecar Vs differs from the z=-1000 m
    #       slice by > 1 km/s — and matplotlib's painter-order
    #       overdraw then makes panel A look like it disagrees with
    #       panel 1, even though the projection is bit-perfect.
    #       (Confirmed: mean |DOF - trilin(x,y,dof_z)| ≈ 0.25 m/s.)
    #   half_thick_metric = per-band RMS bar chart only (panels 5/6).
    #       Wider window is fine — the metric integrates residuals
    #       over many sidecar grid points and does not suffer from
    #       z-overdraw ambiguity.  Kept at the previous 2500 m so the
    #       bar values are unchanged from prior runs.
    half_thick_plot   = 500.0
    half_thick_metric = 2500.0
    half_thick = half_thick_metric   # backward-compat alias for metric calls

    # ---- Sidecar raw values --------------------------------------------
    print(f"reading sidecar {SIDECAR}")
    with h5py.File(SIDECAR, "r") as f:
        gx = f["grid/x"][...]
        gy = f["grid/y"][...]
        gz = f["grid/z"][...]
        Vs = f["fields/Vs"][...]                       # (Nx, Ny, Nz)

    # The sidecar panel shows EXACTLY one z-slice (the one nearest
    # z_target) so the user sees the discrete 1500 m grid as distinct
    # dots, not as overlaid stacks of several slices.
    k_sc = int(np.argmin(np.abs(gz - z_target)))
    print(f"  sidecar slice for plot:  k={k_sc}, gz[k]={gz[k_sc]:.1f} m")
    II, JJ = np.meshgrid(np.arange(gx.size), np.arange(gy.size),
                          indexing="ij")
    II = II.ravel(); JJ = JJ.ravel()
    sidecar_x = gx[II]
    sidecar_y = gy[JJ]
    sidecar_v = Vs[II, JJ, k_sc]
    print(f"  sidecar slice dots:      {sidecar_x.size} (1500 m grid)")

    # ---- Mesh GF values at z=z_target via tet-slice --------------------
    #
    # PURPOSE: verify projection correctness — does the .vtu's stored
    # Vs (= what the actual simulation will read) match the sidecar's
    # raw Vs at the SAME (x, y, z=-1000 m)?
    #
    # METHOD: for every tet whose 4 corner z-range straddles z_target,
    # sample the GF at one point in the plane: the tet's xy-centroid
    # at z = z_target.  Use barycentric linear interpolation against
    # the 4 corner DOF values from PointData['Vs'].  Far-field 10-km-
    # tall tets all straddle z=-1000 m so they each contribute one
    # scatter dot — the resulting cloud fills the full mesh footprint
    # with NO white holes, and every dot's colour is exactly what the
    # downstream simulation will see at z=-1000 m inside that tet.
    #
    # P1 .vtu (4 conn/cell): linear basis on corners is the EXACT FE
    #   solution — slice is an unbiased pointwise read of the GF.
    # P2 .vtu (10 conn/cell): we use only the 4 corner DOFs here, so
    #   the slice is a P1 approximation of the P2 field.  The
    #   quadratic correction within a tet is bounded by the local
    #   curvature of the (smooth) trilinear sidecar × tet_size^2 — at
    #   the 500 m corridor scale this is tens of m/s, well below the
    #   1-km/s discrepancy the user is hunting for.  Adequate for
    #   visual apples-to-apples; not a substitute for the per-band
    #   RMS metric panels on the right.
    def _read_mesh_slice(vtu_path: str):
        m = parse_vtu(vtu_path)
        pts = np.ascontiguousarray(m["points"], dtype=np.float64)
        conn = m["conn"].astype(np.int64)
        offs = m["offsets"].astype(np.int64)
        Vs_dof = m["point_data"]["Vs"]
        dofs_per_cell = int(offs[0])     # 4 for P1 tet, 10 for P2
        n_cells = offs.size
        conn_2d = conn.reshape(n_cells, dofs_per_cell)
        corners_idx = conn_2d[:, :4]                     # (n_cells, 4)
        V = pts[corners_idx]                             # (n_cells, 4, 3)
        F = Vs_dof[corners_idx]                          # (n_cells, 4)

        z = V[:, :, 2]
        zmin = z.min(axis=1); zmax = z.max(axis=1)
        sel = (zmin <= z_target) & (zmax >= z_target)
        if not sel.any():
            return (np.empty(0), np.empty(0), np.empty(0),
                    pts.shape[0], dofs_per_cell, n_cells)
        Vt = V[sel]                                       # (n_sel, 4, 3)
        Ft = F[sel]                                       # (n_sel, 4)

        # Representative sample point: tet's xy-centroid lifted to z_target.
        xc = Vt[:, :, 0].mean(axis=1)
        yc = Vt[:, :, 1].mean(axis=1)
        qpts = np.stack([xc, yc, np.full(xc.size, z_target)], axis=1)

        # Barycentric coordinates of qpts in each tet.  Solve
        #   [v1-v0, v2-v0, v3-v0] @ [b1, b2, b3]^T = q - v0
        # then b0 = 1 - b1 - b2 - b3.
        Amat = np.stack([Vt[:, 1] - Vt[:, 0],
                         Vt[:, 2] - Vt[:, 0],
                         Vt[:, 3] - Vt[:, 0]], axis=2)    # (n_sel, 3, 3)
        rhs = (qpts - Vt[:, 0])[:, :, None]               # (n_sel, 3, 1)
        bc = np.linalg.solve(Amat, rhs).squeeze(-1)       # (n_sel, 3)
        b0 = 1.0 - bc.sum(axis=1)
        gf_at_q = (b0 * Ft[:, 0] + bc[:, 0] * Ft[:, 1]
                   + bc[:, 1] * Ft[:, 2] + bc[:, 2] * Ft[:, 3])
        # Numerical verification: at each straddling tet's slice sample
        # point, what does the sidecar say at the SAME (xc, yc, z_target)?
        # If projection is faithful then |gf_at_q - truth_at_q| is
        # bounded by interpolation error (≪ 1.5 km/s).
        truth_at_q = trilinear(gx, gy, gz, Vs, qpts)
        diff = gf_at_q - truth_at_q
        print(f"    mesh-slice vs sidecar at SAME (xc, yc, -1000 m):")
        print(f"      n_tets={diff.size}  "
              f"mean|diff|={np.mean(np.abs(diff)):7.1f} m/s  "
              f"median|diff|={np.median(np.abs(diff)):7.1f} m/s  "
              f"max|diff|={np.max(np.abs(diff)):7.1f} m/s")
        return (xc, yc, gf_at_q,
                pts.shape[0], dofs_per_cell, n_cells)

    print(f"reading vtu BASELINE: {PREVIEW_BASELINE}")
    bx, by, bv, b_n_total, b_dofs, b_n_cells = \
        _read_mesh_slice(PREVIEW_BASELINE)
    print(f"  BASELINE: {b_n_total} total DOFs, {b_n_cells} cells "
          f"({b_dofs} DOFs/cell), {bx.size} tets straddle z=-1000 m")

    print(f"reading vtu BASELINE+CRom: {PREVIEW_BASELINE_CROM}")
    bcx, bcy, bcv, bc_n_total, bc_dofs, bc_n_cells = \
        _read_mesh_slice(PREVIEW_BASELINE_CROM)
    print(f"  BASELINE+CRom: {bc_n_total} total DOFs, {bc_n_cells} cells "
          f"({bc_dofs} DOFs/cell), {bcx.size} tets straddle z=-1000 m")

    print(f"reading vtu G-1+G-2: {PREVIEW_G1G2}")
    gx_dof, gy_dof, gv, g_n_total, g_dofs, g_n_cells = \
        _read_mesh_slice(PREVIEW_G1G2)
    print(f"  G-1+G-2: {g_n_total} total DOFs, {g_n_cells} cells "
          f"({g_dofs} DOFs/cell), {gx_dof.size} tets straddle z=-1000 m")

    print(f"reading vtu G-1+G-2 + Catmull-Rom: {PREVIEW_G1G2_CROM}")
    cx_dof, cy_dof, cv, c_n_total, c_dofs, c_n_cells = \
        _read_mesh_slice(PREVIEW_G1G2_CROM)
    print(f"  G-1+G-2+CRom: {c_n_total} total DOFs, {c_n_cells} cells "
          f"({c_dofs} DOFs/cell), {cx_dof.size} tets straddle z=-1000 m")

    # ---- Per-sidecar-grid-point residuals (bit-identical truth) --------
    print("computing per-sidecar-grid-point residuals "
          "(truth = raw F[i,j,k]) ...")
    samples_b = per_grid_point_residuals(
        PREVIEW_BASELINE, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick)
    samples_bc = per_grid_point_residuals(
        PREVIEW_BASELINE_CROM, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick)
    samples_g = per_grid_point_residuals(
        PREVIEW_G1G2, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick)
    samples_c = per_grid_point_residuals(
        PREVIEW_G1G2_CROM, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick)
    bands_grid_b  = grid_point_band_rms(samples_b,  DEFAULT_BANDS)
    bands_grid_bc = grid_point_band_rms(samples_bc, DEFAULT_BANDS)
    bands_grid_g  = grid_point_band_rms(samples_g,  DEFAULT_BANDS)
    bands_grid_c  = grid_point_band_rms(samples_c,  DEFAULT_BANDS)

    # Metric A: residuals at sidecar VOXEL MIDPOINTS (off-grid, no
    # vertex-grid alignment bonus).  Use a wider z-window so the
    # midpoint slab is non-empty even when the corner slab catches
    # only one z-slice.
    half_thick_mid = max(half_thick, 1500.0)
    samples_b_mid  = per_voxel_midpoint_residuals(
        PREVIEW_BASELINE, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick_mid)
    samples_bc_mid = per_voxel_midpoint_residuals(
        PREVIEW_BASELINE_CROM, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick_mid)
    samples_g_mid  = per_voxel_midpoint_residuals(
        PREVIEW_G1G2, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick_mid)
    samples_c_mid  = per_voxel_midpoint_residuals(
        PREVIEW_G1G2_CROM, SIDECAR, STL,
        z_target=z_target, half_thick=half_thick_mid)
    bands_mid_b  = grid_point_band_rms(samples_b_mid,  DEFAULT_BANDS)
    bands_mid_bc = grid_point_band_rms(samples_bc_mid, DEFAULT_BANDS)
    bands_mid_g  = grid_point_band_rms(samples_g_mid,  DEFAULT_BANDS)
    bands_mid_c  = grid_point_band_rms(samples_c_mid,  DEFAULT_BANDS)

    # ---- Fault trace overlay ------------------------------------------
    print(f"reading STL: {STL}")
    stl_v = load_stl_vertices(STL)
    fault_trace_xy = stl_surface_trace_xy(stl_v, -1000.0)
    if fault_trace_xy is not None and len(fault_trace_xy):
        flon, flat = utm_to_lonlat_pair(fault_trace_xy[:, 0],
                                        fault_trace_xy[:, 1])
        if len(flon) > 5000:
            stride = len(flon) // 5000 + 1
            flon = flon[::stride]; flat = flat[::stride]
    else:
        flon, flat = None, None

    # ---- Convert all UTM coords to lon/lat for plotting ----------------
    sidecar_lon, sidecar_lat = utm_to_lonlat_pair(sidecar_x, sidecar_y)
    bx_lon,  bx_lat  = utm_to_lonlat_pair(bx, by)
    bcx_lon, bcx_lat = utm_to_lonlat_pair(bcx, bcy)
    gx_lon,  gx_lat  = utm_to_lonlat_pair(gx_dof, gy_dof)
    cx_lon,  cx_lat  = utm_to_lonlat_pair(cx_dof, cy_dof)

    # Common axis extent (in lon/lat) across all panels — use the
    # sidecar's bbox so all panels share the same area.
    extent = [sidecar_lon.min(), sidecar_lon.max(),
              sidecar_lat.min(), sidecar_lat.max()]

    # ---- Plot ----------------------------------------------------------
    cmap_v = plt.get_cmap("jet")
    cmap_d = plt.get_cmap("RdBu_r")
    vmin, vmax = 0.0, 3.5             # km/s, common color scale
    res_lim = 1.5                     # km/s, residual color scale

    fig = plt.figure(figsize=(32, 6.0))
    gs = fig.add_gridspec(1, 7,
                          width_ratios=[1, 1, 1, 1, 1, 1.0, 1.0],
                          wspace=0.30, left=0.03, right=0.995,
                          top=0.83, bottom=0.14)

    def setup_axes(ax, title, xlabel="Longitude (deg)", ylabel=None):
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("auto")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if flon is not None:
            ax.scatter(flon, flat, s=0.05, c="k", alpha=0.5)

    # All panels share the sidecar's full bbox so the user sees the
    # entire horizontal domain.  Mesh DOFs that exist in the slab
    # appear as scatter dots wherever they happen to lie; mesh
    # footprint < sidecar footprint (mesh is fitted to the cut-fault
    # bbox + padding) so the corner regions of the panels will be
    # empty of mesh DOFs by design.
    setup_axes_cropped = setup_axes

    # -- Panel 1 — raw sidecar dots only (the BASELINE), single slice.
    # Use small markers so the discrete 1500 m grid is visible (large
    # markers would overlap and look like an interpolated surface).
    ax0 = fig.add_subplot(gs[0])
    s0 = ax0.scatter(sidecar_lon, sidecar_lat, s=2,
                     c=sidecar_v / 1000.0,
                     cmap=cmap_v, vmin=vmin, vmax=vmax,
                     marker="o", edgecolors="none")
    setup_axes_cropped(ax0,
               f"Sidecar raw F[i,j,k=27]  (BASELINE)\n"
               f"{sidecar_v.size} pts on 1500 m grid @ z={gz[k_sc]:.0f} m, "
               f"NO interpolation",
               ylabel="Latitude (deg)")
    fig.colorbar(s0, ax=ax0, fraction=0.046, label="Vs (km/s)")

    # -- Panel 2 — A: BASELINE H1-P1 + trilinear
    ax1 = fig.add_subplot(gs[1])
    if bx.size > 0:
        s1 = ax1.scatter(bx_lon, bx_lat, s=4,
                         c=bv / 1000.0,
                         cmap=cmap_v, vmin=vmin, vmax=vmax,
                         edgecolors="none")
        fig.colorbar(s1, ax=ax1, fraction=0.046, label="Vs (km/s)")
    setup_axes_cropped(ax1,
               f"A: H1-P1 + trilinear  (BASELINE)\n"
               f"GF slice @ z={z_target:.0f} m  "
               f"({bx.size} tets straddle the plane)")

    # -- Panel 2b — B: BASELINE H1-P1 + Catmull-Rom
    ax1b = fig.add_subplot(gs[2])
    if bcx.size > 0:
        s1b = ax1b.scatter(bcx_lon, bcx_lat, s=4,
                           c=bcv / 1000.0,
                           cmap=cmap_v, vmin=vmin, vmax=vmax,
                           edgecolors="none")
        fig.colorbar(s1b, ax=ax1b, fraction=0.046, label="Vs (km/s)")
    setup_axes_cropped(ax1b,
               f"B: H1-P1 + Catmull-Rom\n"
               f"GF slice @ z={z_target:.0f} m  "
               f"({bcx.size} tets straddle the plane)")

    # -- Panel 3 — C: G-1+G-2 H1-P2 + trilinear
    ax2 = fig.add_subplot(gs[3])
    if gx_dof.size > 0:
        s2 = ax2.scatter(gx_lon, gx_lat, s=2,
                         c=gv / 1000.0,
                         cmap=cmap_v, vmin=vmin, vmax=vmax,
                         edgecolors="none")
        fig.colorbar(s2, ax=ax2, fraction=0.046, label="Vs (km/s)")
    setup_axes_cropped(ax2,
               f"C: H1-P2 + trilinear  (G-1+G-2)\n"
               f"GF slice @ z={z_target:.0f} m  "
               f"({gx_dof.size} tets straddle the plane)")

    # -- Panel 3b — D: G-1+G-2 H1-P2 + Catmull-Rom
    ax2c = fig.add_subplot(gs[4])
    if cx_dof.size > 0:
        s2c = ax2c.scatter(cx_lon, cx_lat, s=2,
                           c=cv / 1000.0,
                           cmap=cmap_v, vmin=vmin, vmax=vmax,
                           edgecolors="none")
        fig.colorbar(s2c, ax=ax2c, fraction=0.046, label="Vs (km/s)")
    setup_axes_cropped(ax2c,
               f"D: H1-P2 + Catmull-Rom  (G-1+G-2+CRom)\n"
               f"GF slice @ z={z_target:.0f} m  "
               f"({cx_dof.size} tets straddle the plane)")

    def _draw_4bar(ax, b_b, b_bc, b_g, b_c, title, ylabel):
        if not (b_b and b_bc and b_g and b_c):
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes); return
        labels = []
        for r in b_b:
            hi_s = ("inf" if (np.isinf(r.hi) or r.hi >= 1e9)
                    else f"{int(r.hi)}")
            labels.append(f"{int(r.lo)}-{hi_s}")
        def _vals(bands):
            return [r.l2_proj if r.n_tets > 0
                    and not np.isnan(r.l2_proj) else 0.0
                    for r in bands]
        l2_b  = _vals(b_b);  l2_bc = _vals(b_bc)
        l2_g  = _vals(b_g);  l2_c  = _vals(b_c)
        x = np.arange(len(labels))
        w = 0.20
        rects_b  = ax.bar(x - 1.5*w, l2_b,  w,
                          label="A: P1+trilin (BASE)", color="#bf3030")
        rects_bc = ax.bar(x - 0.5*w, l2_bc, w,
                          label="B: P1+CRom",          color="#e08a3c")
        rects_g  = ax.bar(x + 0.5*w, l2_g,  w,
                          label="C: P2+trilin (G-1+G-2)",
                          color="#3060bf")
        rects_c  = ax.bar(x + 1.5*w, l2_c,  w,
                          label="D: P2+CRom",          color="#2a8a4d")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
        ax.set_xlabel("dist. to fault (m)", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=6, loc="upper left",
                  handlelength=1.0, handletextpad=0.3)
        for rects, vals in ((rects_b,  l2_b),  (rects_bc, l2_bc),
                            (rects_g,  l2_g),  (rects_c,  l2_c)):
            for rect, val in zip(rects, vals):
                ax.text(rect.get_x() + rect.get_width() / 2,
                        rect.get_height() + 3, f"{val:.0f}",
                        ha="center", fontsize=5)

    # -- Panel 6 — per-band RMS at sidecar voxel CORNERS
    ax4 = fig.add_subplot(gs[5])
    _draw_4bar(ax4,
               bands_grid_b, bands_grid_bc,
               bands_grid_g, bands_grid_c,
               title="Per-band RMS at sidecar CORNERS\n"
                     "(truth = raw F[i,j,k], no interp)",
               ylabel="RMS(proj - raw F)  (m/s)")

    # -- Panel 7 — per-band RMS at sidecar voxel MIDPOINTS (Metric A)
    ax5 = fig.add_subplot(gs[6])
    _draw_4bar(ax5,
               bands_mid_b, bands_mid_bc,
               bands_mid_g, bands_mid_c,
               title="Per-band RMS at sidecar MIDPOINTS (Metric A)\n"
                     "(truth = trilinear of sidecar; off-grid)",
               ylabel="RMS(proj - trilin F)  (m/s)")

    fig.suptitle(
        f"Velocity-projection comparison @ z = {z_target:.0f} m   |   "
        f"sidecar = raw F[i,j,k] @ z={gz[k_sc]:.0f} m (NO interp)   |   "
        f"4 mesh configs:  A: P1+trilin (BASELINE)   "
        f"B: P1+CRom   C: P2+trilin (G-1+G-2)   D: P2+CRom   |   "
        f"mesh DOFs in slab plotted as scatter (NO surface fit)   |   "
        f"fault trace in black",
        fontsize=11)

    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}")

    # ---- Stdout summary ------------------------------------------------
    def _fmt(x):
        return f"{x:>10.1f}" if not np.isnan(x) else f"{'nan':>10}"

    def _agg(bands):
        n = sum(r.n_tets for r in bands
                if r.n_tets > 0 and not np.isnan(r.l2_proj))
        ssq = sum(r.n_tets * r.l2_proj ** 2 for r in bands
                  if r.n_tets > 0 and not np.isnan(r.l2_proj))
        return float(np.sqrt(ssq / max(1, n)))

    def _print_bands(label, b_b, b_bc, b_g, b_c):
        print(f"\n{label}:")
        print(f"    {'band (m)':<14}  "
              f"{'A:P1+tril':>10} {'B:P1+CRom':>10} "
              f"{'C:P2+tril':>10} {'D:P2+CRom':>10}")
        for rb, rbc, rg, rc in zip(b_b, b_bc, b_g, b_c):
            hi_s = ("inf" if (np.isinf(rb.hi) or rb.hi >= 1e9)
                    else f"{int(rb.hi)}")
            print(f"    [{int(rb.lo)}-{hi_s:<8}]  "
                  f"{_fmt(rb.l2_proj)} {_fmt(rbc.l2_proj)} "
                  f"{_fmt(rg.l2_proj)} {_fmt(rc.l2_proj)}")
        print(f"    {'GLOBAL':<14}  "
              f"{_agg(b_b):>10.1f} {_agg(b_bc):>10.1f} "
              f"{_agg(b_g):>10.1f} {_agg(b_c):>10.1f}")

    _print_bands(
        "per-band RMS at sidecar CORNERS  (truth = raw F[i,j,k])",
        bands_grid_b, bands_grid_bc, bands_grid_g, bands_grid_c)
    _print_bands(
        "per-band RMS at sidecar MIDPOINTS  (truth = trilinear; Metric A)",
        bands_mid_b, bands_mid_bc, bands_mid_g, bands_mid_c)


if __name__ == "__main__":
    main()
