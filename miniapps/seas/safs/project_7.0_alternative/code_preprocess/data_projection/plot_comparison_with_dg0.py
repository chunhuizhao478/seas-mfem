"""
Generate an updated comparison PNG for the velocity-projection feature.

Layout (per mesh resolution), 5 panels stacked horizontally:

  1. Sidecar (CVM-H raw) at z = -1 km, smooth trilinear.
  2. Current H1-P1 mesh, rendered with TRUE in-tet linear interpolation
     between mesh vertex values on the z = -1 km slab (this faithfully
     reproduces what ParaView shows and exposes the "streak" artefact
     near sharp gradients in the source).
  3. Proposed F-2 (DG-L²(0)) candidate: per-tet centroid sample of the
     trilinear sidecar.  Each tet renders as one constant Vs colour.
  4. Diff (H1 - sidecar) with the cut-fault STL trace overlaid in black,
     and global RMS in the title.
  5. Cell-mean residual L²(Vs) stratified by tet-centroid distance to
     the fault, as a bar chart (H1-P1 today vs F-2 candidate).

The DG-L²(0) panel is synthesised by sampling the trilinear sidecar at
each tet centroid -- exactly what `FieldProjector::ProjectDG0` will
write on real builds (see fault_zone_projection_plan_v2.md §2).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import LinearNDInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/tmp")
from parse_mfem_vtu import parse_vtu  # noqa: E402

import pyproj

ROOT = ("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/"
        "project_7.0_alternative")
SIDECAR = f"{ROOT}/data_projected/velocity_safs.h5"
PREVIEW = f"{ROOT}/data_projected/preview"
STL_DIR = f"{ROOT}/data_cutnwfault"
OUT_DIR = f"{ROOT}/data_projected"

STL = {
    "500m":  f"{STL_DIR}/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut.stl",
    "1000m": f"{STL_DIR}/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_1000m_clean_clip_nwcut.stl",
    "2000m": f"{STL_DIR}/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_clean_clip_nwcut.stl",
}

TRANS = pyproj.Transformer.from_crs(
    "EPSG:26911", "EPSG:4326", always_xy=True)


def utm_to_lonlat_pair(x_m, y_m):
    return TRANS.transform(x_m, y_m)


def load_stl_vertices(path):
    """Return per-triangle vertex array (Ntri, 3 vertices, 3 coords)."""
    with open(path, "rb") as f:
        f.read(80)
        try:
            n = int(np.frombuffer(f.read(4), dtype=np.uint32)[0])
            payload = f.read(50 * n)
            if len(payload) == 50 * n and n > 0:
                buf = np.frombuffer(payload, dtype=np.uint8).reshape(n, 50)
                return np.frombuffer(
                    np.ascontiguousarray(buf[:, 12:48]),
                    dtype=np.float32).reshape(n, 3, 3).astype(np.float64)
        except Exception:
            pass
    triangles = []
    cur = []
    with open(path, "r", errors="ignore") as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith("vertex"):
                cur.append([float(t) for t in ln.split()[1:4]])
                if len(cur) == 3:
                    triangles.append(cur)
                    cur = []
    return np.array(triangles)


def stl_surface_trace_xy(stl_v, z_target=-1000.0, half_thick=600.0):
    """Return the x-y trace of the cut fault surface near z=z_target by
    keeping triangles whose z-range includes the slab.  Used to plot a
    line on the diff panel."""
    zmin = stl_v[:, :, 2].min(axis=1)
    zmax = stl_v[:, :, 2].max(axis=1)
    sel = (zmax >= z_target - half_thick) & (zmin <= z_target + half_thick)
    centroids = stl_v[sel].mean(axis=1)
    return centroids[:, :2]  # x, y in UTM


def trilinear(gx, gy, gz, F, P):
    nx, ny, nz = len(gx), len(gy), len(gz)
    px = np.clip(np.searchsorted(gx, P[:, 0]) - 1, 0, nx - 2)
    py = np.clip(np.searchsorted(gy, P[:, 1]) - 1, 0, ny - 2)
    pz = np.clip(np.searchsorted(gz, P[:, 2]) - 1, 0, nz - 2)
    fx = np.clip((P[:, 0] - gx[px]) / (gx[px + 1] - gx[px]), 0, 1)
    fy = np.clip((P[:, 1] - gy[py]) / (gy[py + 1] - gy[py]), 0, 1)
    fz = np.clip((P[:, 2] - gz[pz]) / (gz[pz + 1] - gz[pz]), 0, 1)
    c00 = F[px, py, pz] * (1 - fx) + F[px + 1, py, pz] * fx
    c10 = F[px, py + 1, pz] * (1 - fx) + F[px + 1, py + 1, pz] * fx
    c01 = F[px, py, pz + 1] * (1 - fx) + F[px + 1, py, pz + 1] * fx
    c11 = F[px, py + 1, pz + 1] * (1 - fx) + F[px + 1, py + 1, pz + 1] * fx
    return ((c00 * (1 - fy) + c10 * fy) * (1 - fz)
            + (c01 * (1 - fy) + c11 * fy) * fz)


def render_h1_at_z_via_tet_slice(pts, conn, vert_vals,
                                 xq_u, yq_u, z_target):
    """Slice the H1-P1 GF at z = z_target by, for each tet that
    straddles the plane, computing the value at (x_c, y_c, z_target)
    via barycentric linear interpolation against the four vertex
    Vs values.  The resulting (x_c, y_c, Vs) point cloud is then fed
    to a 2-D LinearNDInterpolator -- this is a faithful proxy for
    ParaView's tet-locator-based slicing.

    Tets whose z extent does not straddle z_target are skipped.
    """
    p = pts
    z = p[:, 2]
    zmin = z[conn].min(axis=1)
    zmax = z[conn].max(axis=1)
    sel = (zmin <= z_target + 1.0) & (zmax >= z_target - 1.0)
    if sel.sum() == 0:
        return np.full((len(yq_u), len(xq_u)), np.nan)

    tet = conn[sel]
    n_t = tet.shape[0]
    V = p[tet]                # (n_t, 4, 3)
    F = vert_vals[tet]        # (n_t, 4)

    # x,y of the tet's xy-centroid; z fixed to z_target
    xc = V[:, :, 0].mean(axis=1)
    yc = V[:, :, 1].mean(axis=1)
    qpts = np.stack([xc, yc, np.full(n_t, z_target)], axis=1)

    # Barycentric coordinates of qpts in each tet:
    #   solve [v1-v0, v2-v0, v3-v0] @ [b1,b2,b3]^T = q - v0
    A = np.stack([V[:, 1] - V[:, 0],
                  V[:, 2] - V[:, 0],
                  V[:, 3] - V[:, 0]], axis=2)  # (n_t, 3, 3)
    rhs = (qpts - V[:, 0])[:, :, None]                          # (n_t, 3, 1)
    bc = np.linalg.solve(A, rhs).squeeze(-1)                     # (n_t, 3)
    b0 = 1.0 - bc.sum(axis=1)
    Vs_at_q = (b0 * F[:, 0] + bc[:, 0] * F[:, 1]
               + bc[:, 1] * F[:, 2] + bc[:, 2] * F[:, 3])

    # Reduce to a 2D scattered cloud and interpolate on the (xq, yq) grid.
    P = np.column_stack([xc, yc])
    _, uidx = np.unique(P, axis=0, return_index=True)
    interp = LinearNDInterpolator(P[uidx], Vs_at_q[uidx],
                                  fill_value=np.nan)
    XX, YY = np.meshgrid(xq_u, yq_u, indexing="xy")
    return interp(np.column_stack([XX.ravel(),
                                   YY.ravel()])).reshape(XX.shape)


# Backwards-compat alias used below
render_h1_in_tet_linear = render_h1_at_z_via_tet_slice


# Tet 4-point Gauss cubature (degree-of-precision 2) on the reference
# tet vol = 1/6.  Reference simplex vertices: (0,0,0), (1,0,0), (0,1,0),
# (0,0,1).  Barycentric coords of the 4 cubature points (Stroud T3:5-1):
_GA = (5.0 - np.sqrt(5.0)) / 20.0    # ~0.1382
_GB = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0  # ~0.5854
TET_CUB_BC = np.array([
    [_GB, _GA, _GA, _GA],
    [_GA, _GB, _GA, _GA],
    [_GA, _GA, _GB, _GA],
    [_GA, _GA, _GA, _GB],
])
TET_CUB_W = np.full(4, 0.25)  # equal weights, exact for any quadratic
                              # polynomial on a tet


def tet_cell_mean_via_cubature(V_xyz, gx, gy, gz, F):
    """Approximate the cell mean of the trilinear sidecar over each tet
    using 4-point Gauss cubature.  V_xyz: (n_t, 4, 3) tet vertices."""
    n_t = V_xyz.shape[0]
    # 4 cubature points per tet in physical coords:
    qpts = np.einsum("qi,tij->tqj", TET_CUB_BC, V_xyz)  # (n_t, 4, 3)
    qflat = qpts.reshape(-1, 3)
    vals = trilinear(gx, gy, gz, F, qflat).reshape(n_t, 4)
    return (vals * TET_CUB_W).sum(axis=1)


def render_dg0_at_z(pts, conn, gx, gy, gz, F,
                    xq_u, yq_u, z_target, mode="centroid"):
    """Per-tet projected value, rendered via 2-D nearest-tet-centroid
    Voronoi.

    mode:
       'centroid' — single trilinear-sidecar evaluation at the tet
                    centroid.  Cheapest; F-2 v2 §2.1 baseline.
       'volavg'   — 4-point Gauss tet cubature (exact for trilinear).
                    What an L2-projection of the sidecar onto L²(0)
                    would write.  Recommended default once we move
                    past plot demos.
    """
    ec = pts[conn].mean(axis=1)
    mask = (ec[:, 2] >= z_target - 800) & (ec[:, 2] <= z_target + 800)
    if mask.sum() == 0:
        mask = (ec[:, 2] >= z_target - 1500) & (ec[:, 2] <= z_target + 1500)
    ec_xy = ec[mask, :2]
    if mode == "centroid":
        cell_vals = trilinear(gx, gy, gz, F, ec[mask])
    elif mode == "volavg":
        V_xyz = pts[conn[mask]]                                 # (n, 4, 3)
        cell_vals = tet_cell_mean_via_cubature(V_xyz, gx, gy, gz, F)
    else:
        raise ValueError(f"unknown DG-0 render mode: {mode}")
    tree = cKDTree(ec_xy)
    XX, YY = np.meshgrid(xq_u, yq_u, indexing="xy")
    Q = np.column_stack([XX.ravel(), YY.ravel()])
    _, idx = tree.query(Q, k=1)
    return cell_vals[idx].reshape(YY.shape)


def stratified_metrics_h1_vs_dg0(pts, conn, vert_vals,
                                 gx, gy, gz, F, stl_xy,
                                 z_target, half_thick=200.0):
    """For each tet on the slab, compute the **per-tet integrated L²**
    of the residual `proj(x) - f_source(x)` for each of the three
    projections {H1-P1, DG-0 centroid, DG-0 vol-avg}, evaluated at
    the same 4 Gauss tet-cubature points so the comparison is fair.

    Returns (lo, hi, n_tets, l2_h1, l2_dg0_cent, l2_dg0_vol).
    """
    ec = pts[conn].mean(axis=1)
    mask = (ec[:, 2] >= z_target - half_thick) \
         & (ec[:, 2] <= z_target + half_thick)
    if mask.sum() == 0:
        return None
    tet = conn[mask]
    V_xyz = pts[tet]                          # (n, 4, 3)
    Vverts = vert_vals[tet]                   # (n, 4)
    n_t = V_xyz.shape[0]

    # 4 cubature points in physical coords + barycentric-weighted
    # H1-P1 evaluation at those points.
    qpts = np.einsum("qi,tij->tqj", TET_CUB_BC, V_xyz)    # (n, 4, 3)
    h1_at_q = np.einsum("qi,ti->tq", TET_CUB_BC, Vverts)  # (n, 4)
    f_truth = trilinear(gx, gy, gz, F,
                        qpts.reshape(-1, 3)).reshape(n_t, 4)

    cent  = trilinear(gx, gy, gz, F, V_xyz.mean(axis=1))   # (n,)
    cmean = (f_truth * TET_CUB_W).sum(axis=1)              # (n,)

    res_h1 = h1_at_q - f_truth
    res_c  = cent[:, None]  - f_truth
    res_v  = cmean[:, None] - f_truth

    l2_h1 = np.sqrt((res_h1**2 * TET_CUB_W).sum(axis=1))
    l2_c  = np.sqrt((res_c **2 * TET_CUB_W).sum(axis=1))
    l2_v  = np.sqrt((res_v **2 * TET_CUB_W).sum(axis=1))

    if stl_xy is None:
        return None
    tree = cKDTree(stl_xy)
    d, _ = tree.query(V_xyz.mean(axis=1)[:, :2], k=1)
    bands = [(0, 500), (500, 1500), (1500, 3000), (3000, np.inf)]
    out = []
    for lo, hi in bands:
        sel = (d >= lo) & (d < hi)
        if sel.sum() == 0:
            out.append((lo, hi, 0, 0.0, 0.0, 0.0))
            continue
        out.append((lo, hi, int(sel.sum()),
                    float(np.sqrt(np.mean(l2_h1[sel]**2))),
                    float(np.sqrt(np.mean(l2_c [sel]**2))),
                    float(np.sqrt(np.mean(l2_v [sel]**2)))))
    return out


def make_figure_for(tag: str):
    pvtu = (f"{PREVIEW}/projected_velocity_{tag}/"
            f"Cycle000000/proc000000.vtu")
    if not os.path.exists(pvtu):
        print(f"[skip] {tag}: vtu not found")
        return
    print(f"\n=== {tag} ===")

    with h5py.File(SIDECAR, "r") as f:
        gx = f["grid/x"][...]
        gy = f["grid/y"][...]
        gz = f["grid/z"][...]
        Vs = f["fields/Vs"][...]

    m = parse_vtu(pvtu)
    pts = m["points"]
    conn = m["conn"].reshape(-1, 4)
    Vs_h1 = m["point_data"]["Vs"]

    stl_v = (load_stl_vertices(STL[tag])
             if STL.get(tag) and os.path.exists(STL[tag]) else None)
    stl_xy = (stl_v.mean(axis=1)[:, :2] if stl_v is not None else None)
    fault_trace_xy = (stl_surface_trace_xy(stl_v, -1000.0)
                      if stl_v is not None else None)

    z_target = -1000.0
    nx, ny = len(gx), len(gy)
    xq_u = np.linspace(gx.min(), gx.max(), 360)
    yq_u = np.linspace(gy.min(), gy.max(), 280)
    XX, YY = np.meshgrid(xq_u, yq_u, indexing="xy")

    # Sidecar at z_target (sample on the same query grid for fair diff)
    P3 = np.column_stack([XX.ravel(), YY.ravel(),
                          np.full(XX.size, z_target)])
    Vs_sc = trilinear(gx, gy, gz, Vs, P3).reshape(XX.shape)

    Vs_h1_grid = render_h1_at_z_via_tet_slice(pts, conn, Vs_h1,
                                              xq_u, yq_u, z_target)
    Vs_dg0c_grid = render_dg0_at_z(pts, conn, gx, gy, gz, Vs,
                                   xq_u, yq_u, z_target,
                                   mode="centroid")
    Vs_dg0v_grid = render_dg0_at_z(pts, conn, gx, gy, gz, Vs,
                                   xq_u, yq_u, z_target,
                                   mode="volavg")

    # Mask cells outside the mesh footprint (NaN from LinearND)
    mesh_mask = ~np.isnan(Vs_h1_grid)

    diff_h1   = np.where(mesh_mask, Vs_h1_grid   - Vs_sc, np.nan)
    diff_dg0c = np.where(mesh_mask, Vs_dg0c_grid - Vs_sc, np.nan)
    diff_dg0v = np.where(mesh_mask, Vs_dg0v_grid - Vs_sc, np.nan)

    bands = stratified_metrics_h1_vs_dg0(pts, conn, Vs_h1,
                                         gx, gy, gz, Vs, stl_xy,
                                         z_target)

    # UTM -> lon/lat for the imshow extent
    lon_lo, lat_lo = utm_to_lonlat_pair(xq_u.min(), yq_u.min())
    lon_hi, lat_hi = utm_to_lonlat_pair(xq_u.max(), yq_u.max())
    extent = [lon_lo, lon_hi, lat_lo, lat_hi]
    if fault_trace_xy is not None and len(fault_trace_xy):
        flon, flat = utm_to_lonlat_pair(fault_trace_xy[:, 0],
                                        fault_trace_xy[:, 1])
        # subsample for plotting
        if len(flon) > 5000:
            stride = len(flon) // 5000 + 1
            flon = flon[::stride]
            flat = flat[::stride]
    else:
        flon, flat = None, None

    cmap_v = plt.get_cmap("jet")
    cmap_d = plt.get_cmap("RdBu_r")
    vmin, vmax = 0.0, 3.5         # km/s
    diff_lim = 1.5                # km/s

    fig = plt.figure(figsize=(26, 5.6))
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.75],
                          wspace=0.30, left=0.03, right=0.995,
                          top=0.83, bottom=0.14)

    def imshow_map(ax, Z, **kw):
        return ax.imshow(Z, origin="lower", extent=extent,
                         aspect="auto", **kw)

    ax0 = fig.add_subplot(gs[0])
    h0 = imshow_map(ax0, Vs_sc / 1000.0,
                    cmap=cmap_v, vmin=vmin, vmax=vmax)
    ax0.set_title(f"Sidecar (CVM-H raw) at z = -1000 m\n"
                  f"({Vs.shape[0]} x {Vs.shape[1]} UTM grid)",
                  fontsize=10)
    ax0.set_ylabel("Latitude (deg)")
    ax0.set_xlabel("Longitude (deg)")
    if flon is not None:
        ax0.scatter(flon, flat, s=0.05, c="k", alpha=0.5)
    fig.colorbar(h0, ax=ax0, fraction=0.046, label="Vs (km/s)")

    ax1 = fig.add_subplot(gs[1])
    h1 = imshow_map(ax1, np.where(mesh_mask, Vs_h1_grid, np.nan) / 1000.0,
                    cmap=cmap_v, vmin=vmin, vmax=vmax)
    ax1.set_title(f"Current H1-P1 ({tag} mesh)\n"
                  f"in-tet linear interp", fontsize=10)
    ax1.set_xlabel("Longitude (deg)")
    if flon is not None:
        ax1.scatter(flon, flat, s=0.05, c="k", alpha=0.5)
    fig.colorbar(h1, ax=ax1, fraction=0.046, label="Vs (km/s)")

    ax2 = fig.add_subplot(gs[2])
    h2 = imshow_map(ax2, np.where(mesh_mask, Vs_dg0c_grid, np.nan) / 1000.0,
                    cmap=cmap_v, vmin=vmin, vmax=vmax)
    ax2.set_title(f"F-2 DG-L2(0) centroid ({tag})\n"
                  f"single-point sample / cell", fontsize=10)
    ax2.set_xlabel("Longitude (deg)")
    if flon is not None:
        ax2.scatter(flon, flat, s=0.05, c="k", alpha=0.5)
    fig.colorbar(h2, ax=ax2, fraction=0.046, label="Vs (km/s)")

    ax2b = fig.add_subplot(gs[3])
    h2b = imshow_map(ax2b, np.where(mesh_mask, Vs_dg0v_grid, np.nan) / 1000.0,
                     cmap=cmap_v, vmin=vmin, vmax=vmax)
    ax2b.set_title(f"F-2 DG-L2(0) vol-avg ({tag})\n"
                   f"4-pt Gauss tet cubature", fontsize=10)
    ax2b.set_xlabel("Longitude (deg)")
    if flon is not None:
        ax2b.scatter(flon, flat, s=0.05, c="k", alpha=0.5)
    fig.colorbar(h2b, ax=ax2b, fraction=0.046, label="Vs (km/s)")

    ax3 = fig.add_subplot(gs[4])
    h3 = imshow_map(ax3, diff_h1 / 1000.0,
                    cmap=cmap_d, vmin=-diff_lim, vmax=diff_lim)
    rms_h1   = float(np.sqrt(np.nanmean((diff_h1   / 1000.0) ** 2)))
    rms_dg0c = float(np.sqrt(np.nanmean((diff_dg0c / 1000.0) ** 2)))
    rms_dg0v = float(np.sqrt(np.nanmean((diff_dg0v / 1000.0) ** 2)))
    ax3.set_title(
        f"Diff (H1 - sidecar)\n"
        f"H1 RMS={rms_h1:.3f}  cent={rms_dg0c:.3f}  "
        f"volavg={rms_dg0v:.3f}  (km/s)",
        fontsize=9)
    ax3.set_xlabel("Longitude (deg)")
    if flon is not None:
        ax3.scatter(flon, flat, s=0.05, c="k", alpha=0.5)
    fig.colorbar(h3, ax=ax3, fraction=0.046,
                 label="Vs - sidecar (km/s)")

    ax4 = fig.add_subplot(gs[5])
    if bands is not None:
        labels = [f"{lo}-{hi if hi < 1e9 else 'inf'}"
                  for (lo, hi, _, _, _, _) in bands]
        l2_h1   = [v[3] for v in bands]
        l2_dg0c = [v[4] for v in bands]
        l2_dg0v = [v[5] for v in bands]
        x = np.arange(len(labels))
        w = 0.27
        b1 = ax4.bar(x - w, l2_h1,   w, label="H1-P1",
                     color="#bf3030")
        b2 = ax4.bar(x,     l2_dg0c, w, label="DG-0 cent",
                     color="#3060bf")
        b3 = ax4.bar(x + w, l2_dg0v, w, label="DG-0 vol-avg",
                     color="#2a8a4d")
        ax4.set_xticks(x)
        ax4.set_xticklabels(labels, rotation=20, ha="right",
                            fontsize=8)
        ax4.set_xlabel("dist. to fault (m)", fontsize=8)
        ax4.set_ylabel("per-tet L2 |proj - source|  (m/s)",
                       fontsize=9)
        ax4.set_title("Per-tet integrated L2\n(4-pt cubature)",
                      fontsize=9)
        ax4.legend(fontsize=7, loc="upper right")
        for rects, vals in ((b1, l2_h1), (b2, l2_dg0c), (b3, l2_dg0v)):
            for rect, val in zip(rects, vals):
                ax4.text(rect.get_x() + rect.get_width() / 2,
                         rect.get_height() + 3, f"{val:.0f}",
                         ha="center", fontsize=6)

    fig.suptitle(
        f"Velocity-projection comparison @ z = -1 km, mesh = {tag}\n"
        f"plan: fault_zone_projection_plan_v2.md   |   fault trace overlaid in black",
        fontsize=12)

    out_path = (f"{OUT_DIR}/comparison_sidecar_vs_mesh_"
                f"{tag}_z1km_v2.png")
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")

    if bands is not None:
        print("  per-tet integrated L2 |proj - source|  (m/s):")
        print(f"    {'band (m)':<14} {'n_tets':>6} {'H1':>8}"
              f" {'DG0-cent':>10} {'DG0-vol':>10}")
        for lo, hi, n, l2_h, l2_c, l2_v in bands:
            band_str = (f"[{lo}-"
                        f"{('inf' if hi >= 1e9 else int(hi))}]")
            print(f"    {band_str:<14} {n:>6} {l2_h:>8.1f}"
                  f" {l2_c:>10.1f} {l2_v:>10.1f}")
        # Aggregate global per-tet L2:
        n_all = sum(b[2] for b in bands)
        if n_all > 0:
            wsum = lambda i: sum(b[2] * b[i]**2 for b in bands)
            print(f"    {'GLOBAL':<14} {n_all:>6} "
                  f"{np.sqrt(wsum(3)/n_all):>8.1f} "
                  f"{np.sqrt(wsum(4)/n_all):>10.1f} "
                  f"{np.sqrt(wsum(5)/n_all):>10.1f}")
    print(f"  IMG-RMS (km/s, query-grid render):  "
          f"H1={rms_h1:.4f}  DG0-cent={rms_dg0c:.4f}  "
          f"DG0-volavg={rms_dg0v:.4f}")


def main():
    for tag in ("500m", "1000m", "2000m"):
        try:
            make_figure_for(tag)
        except Exception as e:
            print(f"[error] {tag}: {e}")
            raise


if __name__ == "__main__":
    main()
