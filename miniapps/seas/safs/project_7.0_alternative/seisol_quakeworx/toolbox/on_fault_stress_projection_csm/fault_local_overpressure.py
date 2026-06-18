#!/usr/bin/env python3
"""fault_local_overpressure.py — fault-local graded pore-pressure overpressure
for the SAFS CSM gate-removal (PLAN_gate_removal_csm_2026-06-17.md, Phase 2).

Approach (2) in the gate-removal task: at a restraining-bend gate the fault is
over-clamped, so the apparent friction mu_app = tau/sigma_n_eff collapses toward
mu_d and the rupture arrests.  Raising the LOCAL pore pressure reduces the
effective normal stress (Fulton & Saffer 2009 overpressured low-permeability
fault core), which raises mu_app there WITHOUT changing tau (pore pressure is
isotropic): subtract an isotropic DeltaPp*I from the EFFECTIVE stress tensor.

  mu_app(DeltaPp) = tau / (sigma_n_eff - DeltaPp)        (tau unchanged)

WHY A 2D (s,z) FIELD, NOT A 1D DEPTH PROFILE.  The feasibility sweep
(gate_removal_sweep.py) proved a per-facet graded overpressure removes all gates,
but a single DeltaPp per depth applied across the whole strike band OVER-RAISES
the well-oriented facets (they go supercritical, mu_app > mu_s).  So the carrier
is a GRADED field on a 2D (strike-distance s, depth) grid fine enough that each
cell is orientation-homogeneous.  Per cell the overpressure is the smaller of
(a) what brings the worst facet to the propagation floor and (b) what keeps the
best facet just below mu_s and above the tension / lithostatic caps — so it is
NEVER supercritical by construction (it under-raises a heterogeneous cell rather
than over-raise it).

The field is computed ONCE on the fault (where facet normals exist), saved as an
.npz (s_centers km, depth_centers m, DeltaPp[ns,nz] MPa + metadata), and then
sampled at arbitrary points (fault facet centroids for the VTU/verification, or
volume-grid nodes for the ASAGI nc) by `delta_pp_mpa` — the SINGLE shared
carrier, so the nc and the projection agree under the self-check.

Plain-text math only.  Reuses no project physics (pure geometry + the caps).
"""
import os

import numpy as np

EPS = 1.0e-9


def floor_mu(mu_s, mu_d, s_target):
    """mu_app needed for S = (mu_s-mu_app)/(mu_app-mu_d) < s_target."""
    return (mu_s + s_target * mu_d) / (1.0 + s_target)


def strike_distance_km(x, y, strike_az_deg, hypo_xy):
    """Along-strike distance (km) from the hypocenter, NW positive (matches the
    projection tools: su = (sin az, cos az))."""
    su = np.array([np.sin(np.radians(strike_az_deg)),
                   np.cos(np.radians(strike_az_deg))])
    q = np.stack([np.asarray(x).ravel() - hypo_xy[0],
                  np.asarray(y).ravel() - hypo_xy[1]], axis=1)
    return (q @ su) / 1000.0


def _cosine_taper(d, width):
    """0 at d<=0, 1 at d>=width, cosine ramp between (d = distance INTO the band
    from its edge, width = taper width, same units)."""
    if width <= 0:
        return np.where(d > 0, 1.0, 0.0)
    t = np.clip(d / width, 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(np.pi * t))


def build_overpressure_field(s_km, depth_m, sigma_n_eff, tau, Sv_total, Pp_h,
                             seis_mask, gate_bands_km, mu_s, strike_az_deg, hypo_xy,
                             mu_d=0.10, s_target=1.7, lambda_max=0.9,
                             sn_floor_mpa=0.5, mu_s_margin=0.005, s_bin_km=1.0,
                             z_bin_m=250.0, s_taper_km=1.0, z_taper_m=500.0,
                             seis_zkm=(-12.0, -3.0), core_target_mu=None,
                             meta=None):
    """Build the 2D (s, depth) graded overpressure grid from a PROJECTED fault.

    All per-facet arrays are length-N, compression-positive MPa (sigma_n_eff,
    tau, Sv_total, Pp_h).  gate_bands_km = list of (s_lo, s_hi) along-strike
    bands to overpressure.  Returns a dict:
        s_centers  (ns,)  strike-distance bin centers (km)
        z_centers  (nz,)  DEPTH bin centers (m, >=0)
        DPp        (ns,nz) overpressure (MPa, >=0), 0 outside the bands
        meta       dict (recipe params, for provenance + sampling)
    Cells are filled only from seismogenic facets inside a gate band; a
    cosine taper rolls DPp to 0 over s_taper_km / z_taper_m at the edges.
    """
    floor = floor_mu(mu_s, mu_d, s_target)
    # cores are aimed at `aim` (>= floor) so volume-grid blur + the conservative
    # neighborhood cap still leave them above the floor after baking; the
    # subcritical cap below independently bounds every facet < mu_s, so aiming
    # higher never risks pre-slip.  Default aim = floor.
    aim = float(core_target_mu) if core_target_mu else floor
    aim = min(max(aim, floor), mu_s - mu_s_margin)
    s_km = np.asarray(s_km, float)
    depth_m = np.asarray(depth_m, float)
    # ---- grid extent: union of gate bands (+ taper), seismogenic depth (+ taper)
    s_lo = min(b[0] for b in gate_bands_km) - s_taper_km
    s_hi = max(b[1] for b in gate_bands_km) + s_taper_km
    zd_lo = -seis_zkm[1] * 1000.0 - z_taper_m       # shallow edge (m depth)
    zd_hi = -seis_zkm[0] * 1000.0 + z_taper_m       # deep edge (m depth)
    zd_lo = max(zd_lo, 0.0)
    s_edges = np.arange(s_lo, s_hi + 0.5 * s_bin_km, s_bin_km)
    z_edges = np.arange(zd_lo, zd_hi + 0.5 * z_bin_m, z_bin_m)
    s_centers = 0.5 * (s_edges[:-1] + s_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    ns, nz = len(s_centers), len(z_centers)
    DPp = np.zeros((ns, nz))

    # facets eligible to DRIVE the field: seismogenic + inside a band
    in_band = np.zeros(len(s_km), bool)
    for a, b in gate_bands_km:
        in_band |= (s_km >= a) & (s_km < b)
    drive = seis_mask & in_band & np.isfinite(sigma_n_eff) & np.isfinite(tau)
    if not drive.any():
        raise ValueError("build_overpressure_field: no seismogenic facets in "
                         f"the gate bands {gate_bands_km}")

    si = np.clip(np.searchsorted(s_edges, s_km) - 1, 0, ns - 1)
    zi = np.clip(np.searchsorted(z_edges, depth_m) - 1, 0, nz - 1)
    cap_tension = sigma_n_eff - sn_floor_mpa
    cap_litho = lambda_max * Sv_total - Pp_h
    need_f = np.clip(sigma_n_eff - tau / aim, 0.0, None)               # to aim
    subcap_f = np.clip(sigma_n_eff - tau / (mu_s - mu_s_margin), 0.0, None)

    # per-cell: need (worst facet -> floor) and allow = the largest DeltaPp that
    # keeps the cell's best facet below mu_s and above the tension/lithostatic
    # caps.  allow = +inf for empty cells (no constraint).
    need_grid = np.zeros((ns, nz))
    allow_grid = np.full((ns, nz), np.inf)
    for i in range(ns):
        for j in range(nz):
            sel = drive & (si == i) & (zi == j)
            if not sel.any():
                continue
            need_grid[i, j] = float(np.max(need_f[sel]))
            subcap = float(np.min(subcap_f[sel]))
            cap = float(min(np.min(cap_tension[sel]), np.min(cap_litho[sel])))
            allow_grid[i, j] = max(min(subcap, cap), 0.0)
    # Neighborhood (3x3) MIN of `allow`: a bilinear sample at any facet is a
    # convex combination of the 4 surrounding cell-center values, all within the
    # 3x3 of the facet's own cell, so capping each cell's DeltaPp at the 3x3-min
    # allow guarantees the interpolated value never exceeds that facet's own
    # subcritical/tension cap -> NEVER supercritical, NEVER tension (proof in
    # the test_build_never_supercritical_and_roundtrips unit test).
    allow_nb = allow_grid.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            shifted = np.full((ns, nz), np.inf)
            si0, si1 = max(di, 0), ns + min(di, 0)
            di0, di1 = max(-di, 0), ns + min(-di, 0)
            sj0, sj1 = max(dj, 0), nz + min(dj, 0)
            dj0, dj1 = max(-dj, 0), nz + min(-dj, 0)
            shifted[si0:si1, sj0:sj1] = allow_grid[di0:di1, dj0:dj1]
            allow_nb = np.minimum(allow_nb, shifted)
    filled = np.isfinite(allow_grid)
    n_hetero = int(np.sum((need_grid > allow_nb + 1e-9) & filled))
    DPp = np.clip(np.minimum(need_grid, np.where(np.isfinite(allow_nb),
                                                 allow_nb, need_grid)), 0.0, None)

    # cosine tapers -> DPp smoothly to 0 at the band/depth edges (no stress jump)
    # along-strike: distance into the nearest band; depth: into the seismogenic
    sw = np.zeros(ns)
    for a, b in gate_bands_km:
        din = np.minimum(s_centers - a, b - s_centers)        # >0 inside band
        sw = np.maximum(sw, _cosine_taper(din + s_taper_km, s_taper_km)
                        * (din > -s_taper_km))
    zw = np.minimum(
        _cosine_taper(z_centers - zd_lo, z_taper_m),
        _cosine_taper(zd_hi - z_centers, z_taper_m))
    DPp *= sw[:, None] * zw[None, :]

    md = {"recipe_mu_s": float(mu_s), "mu_d": float(mu_d),
          "s_target": float(s_target), "lambda_max": float(lambda_max),
          "floor": float(floor), "core_target_mu_aim": float(aim),
          "s_taper_km": float(s_taper_km), "z_taper_m": float(z_taper_m),
          "sn_floor_mpa": float(sn_floor_mpa),
          "gate_bands_km": [list(map(float, b)) for b in gate_bands_km],
          "n_hetero_cells_under_raised": int(n_hetero),
          "s_bin_km": float(s_bin_km), "z_bin_m": float(z_bin_m),
          "DPp_max_MPa": float(DPp.max()), "DPp_mean_nonzero_MPa":
          float(DPp[DPp > 0].mean()) if (DPp > 0).any() else 0.0,
          "strike_az_deg": float(strike_az_deg),
          "hypo_xy": [float(hypo_xy[0]), float(hypo_xy[1])]}
    if meta:
        md.update(meta)
    return {"s_centers": s_centers, "z_centers": z_centers, "DPp": DPp,
            "meta": md}


def save_overpressure_field(path, field):
    """Save the (s,z) field as .npz (s_centers, z_centers, DPp, meta-as-json)."""
    import json
    np.savez(path, s_centers=field["s_centers"], z_centers=field["z_centers"],
             DPp=field["DPp"], meta_json=json.dumps(field["meta"]))


def load_overpressure_field(path):
    """Load an .npz overpressure field written by save_overpressure_field."""
    import json
    if not os.path.exists(path):
        raise FileNotFoundError(f"overpressure field not found: {path}")
    z = np.load(path, allow_pickle=False)
    return {"s_centers": z["s_centers"], "z_centers": z["z_centers"],
            "DPp": z["DPp"], "meta": json.loads(str(z["meta_json"]))}


def delta_pp_mpa(field, x, y, z, strike_az_deg=None, hypo_xy=None):
    """Bilinearly sample DeltaPp (MPa, >=0) at points (x,y,z) [m].  Returns 0
    outside the (s,depth) grid.  Used for BOTH fault facet centroids and
    volume-grid nodes so the projection and the nc carry the identical field.
    strike_az / hypo default to the values stored in the field metadata at
    build time (so all consumers use the SAME strike frame)."""
    if strike_az_deg is None:
        strike_az_deg = field["meta"]["strike_az_deg"]
    if hypo_xy is None:
        hypo_xy = field["meta"]["hypo_xy"]
    sc = field["s_centers"]
    zc = field["z_centers"]
    DPp = field["DPp"]
    s_km = strike_distance_km(x, y, strike_az_deg, hypo_xy)
    depth = np.maximum(-np.asarray(z, float).ravel(), 0.0)

    def locate(g, q):
        if len(g) < 2:
            return np.zeros(len(q), int), np.zeros(len(q)), np.ones(len(q), bool)
        i = np.clip(np.searchsorted(g, q) - 1, 0, len(g) - 2)
        w = (q - g[i]) / (g[i + 1] - g[i])
        inside = (q >= g[0]) & (q <= g[-1])
        return i, np.clip(w, 0.0, 1.0), inside

    si, ws, in_s = locate(sc, s_km)
    zj, wz, in_z = locate(zc, depth)
    out = ((1 - ws) * (1 - wz) * DPp[si, zj]
           + ws * (1 - wz) * DPp[np.minimum(si + 1, len(sc) - 1), zj]
           + (1 - ws) * wz * DPp[si, np.minimum(zj + 1, len(zc) - 1)]
           + ws * wz * DPp[np.minimum(si + 1, len(sc) - 1),
                           np.minimum(zj + 1, len(zc) - 1)])
    out = np.where(in_s & in_z, out, 0.0)
    return np.maximum(out, 0.0).reshape(np.shape(z))


def apply_overpressure_to_tensor(sigma_comp_pos, dpp_mpa):
    """Subtract the isotropic overpressure from a compression-POSITIVE tensor
    stack (N,3,3) MPa, in place-safe: sigma_eff_new = sigma_eff - DeltaPp*I.
    Off-diagonals (shear) are untouched so tau is exact."""
    out = sigma_comp_pos.copy()
    for i in range(3):
        out[:, i, i] -= dpp_mpa
    return out
