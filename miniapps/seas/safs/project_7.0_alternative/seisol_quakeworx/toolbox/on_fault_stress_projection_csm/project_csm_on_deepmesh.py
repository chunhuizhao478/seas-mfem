#!/usr/bin/env python3
"""project_csm_on_deepmesh.py — apply the SAME CSM C1 k=2.39 on-fault stress
projection used for safs_seisol_v3 to a DIFFERENT fault mesh (a Gmsh v2.2 .msh,
physical-surface fault tag 103), and regenerate the strike-distance vs depth
mu_app / S-ratio figure (the csm_S_ratio_const_mu_strike_depth.png layout) so we
can check whether the SAME gates appear on the deeper mesh.

Reuses project_csm_stress_to_vtu.py verbatim for the physics (read_csm_csv,
geographic_to_utm11n, linear CSM tensor interpolation, magnitudes_C1,
build_tensor_from_axes, triangle_geometry, harmonise_normals, tandem_basis,
resolve_tractions) so the projection is byte-for-byte the production config
(axes=csm, sampling=linear, closure=ratio k=2.39, hydrostatic Pp, MUSCAL Sv).
The ONLY new code is the Gmsh-.msh fault reader (the production tool reads PUML).

Outputs (into this toolbox folder):
  csm_C1k239_deepmesh_strike_depth.png   the 2-panel figure
  csm_C1k239_deepmesh_fault_stress.vtu   per-facet mu_app/S-ratio for ParaView
  csm_C1k239_deepmesh_gates.json         detected gate bands + comparison

Run with the env that has scipy+pyproj+matplotlib+netCDF4 (pythonenv):
  /Users/chunhuizhao/miniforge/envs/pythonenv/bin/python project_csm_on_deepmesh.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_csm_stress_to_vtu import (  # noqa: E402
    read_csm_csv, csm_tensors_tension, interpolate_csm_field,
    geographic_to_utm11n, muscal_sv_total_profile, sv_total_at, pore_pressure,
    magnitudes_C1, build_tensor_from_axes, triangle_geometry,
    harmonise_normals, tandem_basis, resolve_tractions, write_vtu, EPS)

# ---- config: identical physics to safs_seisol_v3_0_0_LSW -------------------
ALT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))  # project_7.0_alternative
CSM_CSV = os.path.join(ALT, "seisol_quakeworx/raw_data/"
                       "yang_and_hauksson_orientation/CSM_data_1781290479811.csv")
MAT_NC = os.path.join(ALT, "seisol_quakeworx/safs_seisol_v2_1_0_RSSRW/"
                      "safs_material_cvm.nc")
MESH = os.path.join(ALT, "seisol_quakeworx/safs_seisol_v2_1_0_RSSRW/"
                    "preprocess/safs_mesh_deep.msh")
FAULT_TAG = 103
K_RATIO = 2.39
MU_S, MU_D = 0.47, 0.10
STRIKE_AZ = 314.0
S_ARREST = 1.77                                # S below this propagates
HYPO = np.array([606971.0, 3707270.0, -4965.62])
# v3 reference gates (PUML mesh) for the side-by-side comparison
V3_GATES = [(-20, -10, "SE-tip", 94), (45, 65, "San Gorgonio", 95),
            (180, 190, "", 81), (220, 245, "", 94)]


def read_gmsh_fault(path, fault_tag):
    """Parse a Gmsh v2.2 ASCII mesh; return (pts, tris) for the triangles whose
    physical tag == fault_tag, deduplicated and with a compacted point array
    (same (pts, tris) contract as project_csm_stress_to_vtu.extract_fault)."""
    with open(path) as fh:
        line = fh.readline()
        while line and not line.startswith("$Nodes"):
            line = fh.readline()
        nn = int(fh.readline())
        coords = np.empty((nn, 3))
        nid = np.empty(nn, dtype=np.int64)
        for i in range(nn):
            p = fh.readline().split()
            nid[i] = int(p[0])
            coords[i] = (float(p[1]), float(p[2]), float(p[3]))
        # node-id -> row (ids may be 1-based / non-contiguous)
        id2row = {int(k): i for i, k in enumerate(nid)}
        while line and not line.startswith("$Elements"):
            line = fh.readline()
        ne = int(fh.readline())
        tri = []
        for _ in range(ne):
            p = fh.readline().split()
            if int(p[1]) != 2:                 # 2 = 3-node triangle
                continue
            ntag = int(p[2])
            if int(p[3]) != fault_tag:          # physical tag is the first tag
                continue
            n = p[3 + ntag:3 + ntag + 3]
            tri.append((id2row[int(n[0])], id2row[int(n[1])], id2row[int(n[2])]))
    tris = np.asarray(tri, dtype=np.int64)
    if len(tris) == 0:
        sys.exit(f"ERROR: no triangles with physical tag {fault_tag} in {path}")
    n_raw = len(tris)
    _, keep = np.unique(np.sort(tris, axis=1), axis=0, return_index=True)
    tris = tris[np.sort(keep)]
    used, inv = np.unique(tris.ravel(), return_inverse=True)
    print(f"fault tag {fault_tag}: {n_raw} triangles -> {len(tris)} unique "
          f"(dedup removed {n_raw - len(tris)})")
    return coords[used], inv.reshape(tris.shape)


def detect_gates(s, zkm, S_ratio, barrier, bin_km=5.0,
                 seis=(-12.0, -3.0), frac_thr=0.70):
    """Auto-detect arrest gates: bin along strike, compute the seismogenic-band
    fraction of facets that ARREST (pure barrier OR S>=S_ARREST), and return
    contiguous runs of bins above frac_thr as (s_lo, s_hi, peak_pct)."""
    m = (zkm > seis[0]) & (zkm < seis[1])
    arrest = barrier | (S_ratio >= S_ARREST)
    edges = np.arange(np.floor(s.min()), np.ceil(s.max()) + bin_km, bin_km)
    frac = np.full(len(edges) - 1, np.nan)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        sel = m & (s >= a) & (s < b)
        if sel.sum() >= 10:
            frac[i] = arrest[sel].mean()
    gates, i = [], 0
    hot = frac >= frac_thr
    while i < len(hot):
        if hot[i]:
            j = i
            while j + 1 < len(hot) and hot[j + 1]:
                j += 1
            gates.append((float(edges[i]), float(edges[j + 1]),
                          float(np.nanmax(frac[i:j + 1]) * 100)))
            i = j + 1
        else:
            i += 1
    return gates, edges, frac


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # ---- fault geometry ----
    pts, tris = read_gmsh_fault(MESH, FAULT_TAG)
    cent, nraw, areas = triangle_geometry(pts, tris)
    normals, _ = harmonise_normals(nraw, STRIKE_AZ)
    strikes, dips, degen = tandem_basis(normals)
    print(f"facets {len(tris)}; degenerate basis {int(degen.sum())}; "
          f"centroid z [{cent[:,2].min():.0f}, {cent[:,2].max():.0f}] m")

    # ---- depth-dependent magnitudes (MUSCAL Sv - hydrostatic Pp) ----
    depth = np.maximum(-cent[:, 2], 0.0)
    dg, sg, _ = muscal_sv_total_profile(MAT_NC)
    Sv_total = sv_total_at(depth, dg, sg)
    Pp = pore_pressure(depth, "hydrostatic")
    Sv_eff = Sv_total - Pp
    nbad = int(np.sum(Sv_eff <= 0))
    if nbad:
        print(f"warning: {nbad} facets with Sv_eff<=0 (shallow); clamped to EPS")
        Sv_eff = np.maximum(Sv_eff, EPS)

    # ---- CSM orientation + shape ratio: LINEAR tensor interpolation (production) ----
    lon, lat, _dep, S, *_ = read_csm_csv(CSM_CSV)
    cx, cy = geographic_to_utm11n(lon, lat)
    S_facet, n_fb = interpolate_csm_field(cx, cy, S, cent[:, 0], cent[:, 1])
    if n_fb:
        print(f"warning: {n_fb}/{len(tris)} facets outside CSM hull -> NN fallback")
    T = csm_tensors_tension(S_facet)
    w, U = np.linalg.eigh(T)
    denom = w[:, 2] - w[:, 0]
    R = np.clip((w[:, 1] - w[:, 0]) / np.where(denom > EPS, denom, np.nan),
                0.0, 1.0)

    # ---- C1 k=2.39 closure -> tensor -> on-fault tractions ----
    sig1, sig2, sig3 = magnitudes_C1(Sv_eff, R, K_RATIO)
    sigma = build_tensor_from_axes(U, sig1, sig2, sig3)
    r = resolve_tractions(sigma, strikes, dips, normals, Pp)
    mu = r["mu_apparent"]
    sn = r["sigma_n_eff"]

    # ---- const-mu S-ratio (mu_s=0.47, mu_d=0.10) ----
    barrier = mu < MU_D                              # pure barrier: never weakens
    with np.errstate(divide="ignore", invalid="ignore"):
        S_ratio = np.where(mu > MU_D, (MU_S - mu) / (mu - MU_D), np.nan)

    # ---- strike distance s (hypocenter origin, NW positive) ----
    su = np.array([np.sin(np.radians(STRIKE_AZ)), np.cos(np.radians(STRIKE_AZ))])
    s = ((cent[:, :2] - HYPO[:2]) @ su) / 1000.0
    zkm = cent[:, 2] / 1000.0

    fin = np.isfinite(mu)
    print(f"\nmu_app   : min {np.nanmin(mu):.3f} med {np.nanmedian(mu):.3f} "
          f"max {np.nanmax(mu):.3f}")
    print(f"strike s : SE end {s.min():.0f} km, hypocenter 0, NW end {s.max():.0f} km")
    seis = (zkm > -12) & (zkm < -3) & fin
    print(f"seismogenic mu_app median {np.nanmedian(mu[seis]):.3f}; "
          f"pure-barrier facets {int(barrier.sum())}")

    gates, edges, frac = detect_gates(s, zkm, S_ratio, barrier)
    print(f"\n=== DETECTED GATES (deep mesh; seismogenic arrest fraction >=70%) ===")
    for a, b, pct in gates:
        print(f"  s = {a:6.0f} .. {b:6.0f} km   peak arrest {pct:.0f}%")
    print("=== v3 reference gates (PUML mesh) ===")
    for a, b, nm, pct in V3_GATES:
        print(f"  s = {a:6.0f} .. {b:6.0f} km   {nm:13s} {pct}%")

    # ---------------------------------------------------------------- figure
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    order = np.argsort(zkm)                         # deep first so shallow on top
    sc0 = ax0.scatter(s[order], zkm[order], c=mu[order], s=3,
                      cmap="viridis", vmin=0.0, vmax=0.45, linewidths=0)
    ax0.set_title("mu_app (CSM C1 k=2.39 stress) on the DEEP SAFS fault "
                  "(safs_mesh_deep, fault tag 103)")
    ax0.set_ylabel("depth (km)")
    fig.colorbar(sc0, ax=ax0, label="mu_app", pad=0.01)

    Scl = np.clip(S_ratio, 0, 3.5)
    sc1 = ax1.scatter(s[order], zkm[order], c=Scl[order], s=3,
                      cmap="RdYlGn_r", vmin=0.0, vmax=3.5, linewidths=0)
    ax1.scatter(s[barrier], zkm[barrier], c="k", s=4, label="pure barrier",
                linewidths=0)
    ax1.set_title("S-ratio = (mu_s-mu_app)/(mu_app-mu_d), mu_s=0.47 mu_d=0.10 : "
                  "GREEN S<1.77 propagates, RED arrests (GATE)")
    ax1.set_ylabel("depth (km)")
    ax1.set_xlabel(f"strike distance s (km)   [SE end ~ {s.min():.0f}, "
                   f"hypocenter 0, NW end ~ {s.max():.0f}]")
    ax1.legend(loc="lower right", markerscale=2, framealpha=0.9)
    fig.colorbar(sc1, ax=ax1, label="S-ratio (clip 3.5)", pad=0.01)

    zlo, zhi = zkm.min() - 0.5, zkm.max() + 0.5
    for ax in (ax0, ax1):
        ax.axvline(0.0, color="blue", lw=1)
        ax.set_ylim(zlo, zhi)
        for a, b, pct in gates:
            for axx in (ax0, ax1):
                axx.add_patch(Rectangle((a, zlo), b - a, zhi - zlo, fill=False,
                                        edgecolor="red", lw=1.8, ls="--"))
        ax.margins(x=0.01)
    ax1.text(0.5, zlo + 0.6, "hypocenter", color="blue", fontsize=8)
    for a, b, pct in gates:
        ax0.text((a + b) / 2, zhi - 0.3, f"gate\ns={a:.0f}..{b:.0f} km ({pct:.0f}%)",
                 color="darkred", ha="center", va="top", fontsize=8)

    fig.tight_layout()
    png = os.path.join(HERE, "csm_C1k239_deepmesh_strike_depth.png")
    fig.savefig(png, dpi=130)
    print(f"\nwrote {png}")

    # ---- VTU for ParaView ----
    vtu = os.path.join(HERE, "csm_C1k239_deepmesh_fault_stress.vtu")
    write_vtu(vtu, pts, tris, 5, cell_data={
        "mu_apparent_cell": mu, "sigma_n_eff_MPa_cell": sn,
        "tau_magnitude_MPa_cell": r["tau_magnitude"],
        "S_ratio_cell": np.nan_to_num(S_ratio, nan=-1.0),
        "pure_barrier_cell": barrier.astype(np.int32),
        "strike_distance_km_cell": s, "csm_R_cell": R,
        "Sv_eff_MPa_cell": Sv_eff})
    print(f"wrote {vtu}")

    out = {"mesh": MESH, "fault_tag": FAULT_TAG, "n_facets": int(len(tris)),
           "k_ratio": K_RATIO, "mu_s": MU_S, "mu_d": MU_D,
           "strike_s_range_km": [float(s.min()), float(s.max())],
           "mu_app_seismogenic_median": float(np.nanmedian(mu[seis])),
           "pure_barrier_facets": int(barrier.sum()),
           "detected_gates_s_km": [[a, b, p] for a, b, p in gates],
           "v3_reference_gates_s_km": [[a, b, nm, p] for a, b, nm, p in V3_GATES]}
    jpath = os.path.join(HERE, "csm_C1k239_deepmesh_gates.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {jpath}")


if __name__ == "__main__":
    main()
