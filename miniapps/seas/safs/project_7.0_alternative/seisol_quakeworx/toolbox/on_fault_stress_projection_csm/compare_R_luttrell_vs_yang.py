#!/usr/bin/env python3
"""compare_R_luttrell_vs_yang.py

Compare the CSM stress-shape ratio
    R = (S2 - S3) / (S1 - S3)        (Gephart & Forsyth 1984, tension-positive)
between two Community Stress Models, sampled on the SAFS dynamic-rupture fault:

  * Luttrell-2017   (raw_data/luttrell_differential_stress)  -- meaningful
                     orientation AND magnitude
  * YHSM-2013       (raw_data/yang_and_hauksson_orientation)  -- orientation
                     only (the model V3 actually uses for the C1 k=2.39 stress)

R is column 13 of both CSV exports.  Both models are depth-invariant, so R on
the fault is a function of map position only; we sample R(lon,lat) onto each
fault facet centroid (nearest-neighbour in UTM 11N, the projection pipeline's
convention).

The "large S region" is read straight from the same const-mu S-ratio field the
user's strike-depth figure was built from
(safs_seisol_v3_0_0_LSW/csm_S_ratio_const_mu_fault.vtu : S_ratio_cell /
mu_app_cell).  We highlight where the fault arrests (S >= S_arrest) and compare
R there vs. the propagating remainder.

Outputs (into --outdir):
  compare_R_luttrell_vs_yang_strike_depth.png   2D scatter, both models + dR
  compare_R_luttrell_vs_yang_along_strike.png   1D R(s) curves + dR + S overlay
  compare_R_luttrell_vs_yang_summary.json       statistics, incl. large-S region
"""
import argparse
import json
import os
import re
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_csm_stress_to_vtu import (  # noqa: E402
    COL_LON, COL_LAT, COL_R, load_puml, sanity_check_face_ordering,
    extract_fault, triangle_geometry, geographic_to_utm11n)

_VTK_RD = {"Float64": "<f8", "Float32": "<f4", "Int64": "<i8",
           "Int32": "<i4", "UInt8": "u1"}


def read_R_lonlat(path):
    """Parse (lon, lat, R) from a CSM csv; drop rows with non-finite R."""
    lon, lat, R = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\n").split(",")
            try:
                lo = float(row[COL_LON].strip().strip('"'))
                la = float(row[COL_LAT].strip().strip('"'))
                rs = row[COL_R].strip()
                rv = float(rs) if rs not in ("", "NaN", "nan") else np.nan
            except (ValueError, IndexError):
                continue
            lon.append(lo); lat.append(la); R.append(rv)
    lon, lat, R = map(np.asarray, (lon, lat, R))
    good = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(R)
    if not good.any():
        sys.exit(f"ERROR: no finite R rows in {path}")
    return lon[good], lat[good], R[good]


def read_vtu_cell_arrays(path, names):
    with open(path, "rb") as fh:
        raw = fh.read()
    marker = b'<AppendedData encoding="raw">\n_'
    mk = raw.find(marker)
    if mk < 0:
        sys.exit(f"ERROR: {path} is not the expected appended-raw VTU")
    start = mk + len(marker)
    hdr = raw[:mk].decode("latin1")
    arrs = re.findall(
        r'<DataArray type="([^"]+)" Name="([^"]+)"'
        r'(?: NumberOfComponents="(\d+)")? format="appended" offset="(\d+)"/>',
        hdr)
    out = {}
    for typ, nm, nc, off in arrs:
        if nm not in names:
            continue
        off = int(off); nc = int(nc) if nc else 1
        nbytes = struct.unpack_from("<Q", raw, start + off)[0]
        buf = raw[start + off + 8:start + off + 8 + nbytes]
        a = np.frombuffer(buf, dtype=_VTK_RD[typ])
        out[nm] = a.reshape(-1, nc) if nc > 1 else a
    missing = [n for n in names if n not in out]
    if missing:
        sys.exit(f"ERROR: {path} missing cell arrays {missing}")
    return out


def nearest_sample(cx, cy, vals, qx, qy):
    """Nearest-neighbour lookup of `vals` defined at (cx,cy) onto (qx,qy)."""
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        sys.exit("ERROR: scipy is required (conda env 'pythonenv')")
    tree = cKDTree(np.stack([cx, cy], axis=1))
    d, idx = tree.query(np.stack([qx, qy], axis=1), k=1)
    return vals[idx], d


def strike_coordinate(cx, cy, hypo_xy, strike_az_deg=314.0, n_bins=160):
    """Arc-length strike coordinate s (km) for each facet centroid.

    Build a smooth fault-trace polyline by binning centroids along the PCA
    primary axis, then assign each facet the cumulative arc length of its
    nearest polyline vertex.  Origin at the hypocenter, NW positive."""
    P = np.stack([cx, cy], axis=1)
    mean = P.mean(axis=0)
    Pc = P - mean
    _, _, Vt = np.linalg.svd(Pc, full_matrices=False)
    u = Vt[0]                                   # primary (strike) axis
    nw = np.array([np.sin(np.radians(strike_az_deg)),
                   np.cos(np.radians(strike_az_deg))])
    if u @ nw < 0:
        u = -u
    t = Pc @ u
    # ordered trace polyline: bin by t-quantile, mean position per bin
    edges = np.quantile(t, np.linspace(0, 1, n_bins + 1))
    edges[-1] += 1e-6
    bidx = np.clip(np.digitize(t, edges) - 1, 0, n_bins - 1)
    trace = []
    for b in range(n_bins):
        sel = bidx == b
        if sel.any():
            trace.append(P[sel].mean(axis=0))
    trace = np.asarray(trace)
    seg = np.linalg.norm(np.diff(trace, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])         # m
    # nearest trace vertex for each facet and for the hypocenter
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        sys.exit("ERROR: scipy is required (conda env 'pythonenv')")
    tt = cKDTree(trace)
    _, fi = tt.query(P, k=1)
    s = arc[fi]
    _, hi = tt.query(np.asarray(hypo_xy)[None, :], k=1)
    s0 = arc[hi[0]]
    return (s - s0) / 1000.0                                # km, origin=hypo


def field_stats(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return dict(n=0)
    return dict(n=int(a.size), min=float(a.min()), p05=float(np.percentile(a, 5)),
                median=float(np.median(a)), mean=float(a.mean()),
                p95=float(np.percentile(a, 95)), max=float(a.max()))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    raw = os.path.normpath(os.path.join(here, "..", "..", "raw_data"))
    v3 = os.path.normpath(os.path.join(here, "..", "..",
                                       "safs_seisol_v3_0_0_LSW"))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--luttrell-csv",
                    default=os.path.join(raw, "luttrell_differential_stress",
                                         "CSM_data_1781537253009.csv"))
    ap.add_argument("--yang-csv",
                    default=os.path.join(raw, "yang_and_hauksson_orientation",
                                         "CSM_data_1781290479811.csv"))
    ap.add_argument("--mesh",
                    default=os.path.join(here, "..", "..",
                                         "safs_seisol_v2_0_0_RSSRW",
                                         "safs_mesh.puml.h5"))
    ap.add_argument("--sratio-vtu",
                    default=os.path.join(v3, "csm_S_ratio_const_mu_fault.vtu"))
    ap.add_argument("--hypocenter", type=float, nargs=3,
                    default=[606971.0, 3707270.0, -4965.62])
    ap.add_argument("--s-arrest", type=float, default=1.77,
                    help="S-ratio arrest threshold (green<S propagates)")
    ap.add_argument("--seis-ztop", type=float, default=-3000.0)
    ap.add_argument("--seis-zbot", type=float, default=-15000.0)
    ap.add_argument("--outdir", default=v3)
    args = ap.parse_args()

    # --- fault geometry (extract_fault order == the S-ratio VTU order) ---
    geom, conn, bc = load_puml(args.mesh)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    centroids, _n, _a = triangle_geometry(pts, tris)
    n_facets = len(tris)
    cxg, cyg, czg = centroids[:, 0], centroids[:, 1], centroids[:, 2]
    z = czg

    # --- mu_app / S_ratio (the field the user's figure was built from) ---
    cd = read_vtu_cell_arrays(args.sratio_vtu,
                              ["S_ratio_cell", "mu_app_cell", "barrier_cell"])
    if len(cd["mu_app_cell"]) != n_facets:
        sys.exit(f"ERROR: S-ratio VTU has {len(cd['mu_app_cell'])} facets, "
                 f"mesh fault has {n_facets} -- mismatch")
    S = cd["S_ratio_cell"]
    mu_app = cd["mu_app_cell"]
    barrier = cd["barrier_cell"].astype(bool)

    # --- facet centroids -> lon/lat?  No: CSM is in lon/lat -> UTM, sample
    #     nearest in UTM (model depth-invariant, 2D). ---
    loL, laL, RL = read_R_lonlat(args.luttrell_csv)
    loY, laY, RY = read_R_lonlat(args.yang_csv)
    cxL, cyL = geographic_to_utm11n(loL, laL)
    cxY, cyY = geographic_to_utm11n(loY, laY)
    R_lut, dL = nearest_sample(cxL, cyL, RL, cxg, cyg)
    R_yng, dY = nearest_sample(cxY, cyY, RY, cxg, cyg)
    dR = R_lut - R_yng

    # --- strike coordinate (arc length, origin hypocenter, NW positive) ---
    s = strike_coordinate(cxg, cyg, np.asarray(args.hypocenter[:2]))

    # --- region masks ---
    seis = (z < args.seis_ztop) & (z > args.seis_zbot) & (~barrier)
    largeS = seis & np.isfinite(S) & (S >= args.s_arrest)     # arrest band
    propag = seis & np.isfinite(S) & (S < args.s_arrest)      # propagating
    pure_barrier = seis & (mu_app <= 0.12)                    # mu_app<=mu_d

    # --- console summary ---
    print(f"facets: {n_facets};  nearest-sample max dist: "
          f"Luttrell {dL.max():.0f} m, Yang {dY.max():.0f} m")
    print(f"strike span: SE {s.min():+.1f} km  hypo 0  NW {s.max():+.1f} km")
    print(f"\nR over ALL fault facets:")
    print(f"  Luttrell : {field_stats(R_lut)}")
    print(f"  Yang     : {field_stats(R_yng)}")
    print(f"  dR=L-Y   : {field_stats(dR)}")
    for name, m in (("seismogenic band", seis),
                    (f"LARGE-S arrest band (S>={args.s_arrest})", largeS),
                    (f"propagating band (S<{args.s_arrest})", propag),
                    ("pure-barrier (mu_app<=0.12)", pure_barrier)):
        if m.sum() == 0:
            print(f"\n[{name}] : 0 facets"); continue
        print(f"\n[{name}]  ({int(m.sum())} facets, "
              f"{100*m.sum()/n_facets:.1f}%)")
        print(f"  R_Luttrell : {field_stats(R_lut[m])}")
        print(f"  R_Yang     : {field_stats(R_yng[m])}")
        print(f"  dR = L-Y   : {field_stats(dR[m])}")

    # --- agreement metrics + targeted arrest windows ---
    def agree(mask):
        d = dR[mask]; d = d[np.isfinite(d)]
        a = R_lut[mask]; b = R_yng[mask]
        ok = np.isfinite(a) & np.isfinite(b)
        cc = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else None
        return dict(n=int(mask.sum()), rms_dR=float(np.sqrt(np.mean(d**2))),
                    mean_abs_dR=float(np.mean(np.abs(d))),
                    p95_abs_dR=float(np.percentile(np.abs(d), 95)),
                    max_abs_dR=float(np.max(np.abs(d))),
                    frac_absdR_gt_0p05=float(np.mean(np.abs(d) > 0.05)),
                    frac_absdR_gt_0p10=float(np.mean(np.abs(d) > 0.10)),
                    corr_R=cc)

    print("\n=== AGREEMENT (R_Luttrell vs R_Yang) ===")
    for nm, m in (("all", np.ones(n_facets, bool)), ("seismogenic", seis),
                  ("large_S_arrest", largeS)):
        ag = agree(m)
        print(f"  [{nm}] corr={ag['corr_R']:.4f}  RMS(dR)={ag['rms_dR']:.4f}  "
              f"mean|dR|={ag['mean_abs_dR']:.4f}  p95|dR|={ag['p95_abs_dR']:.3f}"
              f"  |dR|>0.05: {100*ag['frac_absdR_gt_0p05']:.1f}%  "
              f"|dR|>0.10: {100*ag['frac_absdR_gt_0p10']:.1f}%")

    windows = {"SE_arrest_gate_s35_65": (35.0, 65.0),
               "central_hi_R_s80_115": (80.0, 115.0),
               "NW_arrest_s185_260": (185.0, 260.0)}
    print("\n=== TARGETED STRIKE WINDOWS (seismogenic band) ===")
    win_summary = {}
    for nm, (lo, hi) in windows.items():
        m = seis & (s >= lo) & (s < hi)
        if m.sum() == 0:
            continue
        rl, ry = field_stats(R_lut[m]), field_stats(R_yng[m])
        ag = agree(m)
        print(f"  [{nm}]  ({int(m.sum())} facets)  "
              f"R_Lut med {rl['median']:.3f}  R_Yang med {ry['median']:.3f}  "
              f"dR med {field_stats(dR[m])['median']:+.3f}  "
              f"RMS(dR) {ag['rms_dR']:.3f}  max|dR| {ag['max_abs_dR']:.3f}")
        win_summary[nm] = {"s_window_km": [lo, hi], "R_luttrell": rl,
                           "R_yang": ry, "dR_L_minus_Y": field_stats(dR[m]),
                           "agreement": ag}

    summary = {
        "n_facets": n_facets,
        "nearest_max_dist_m": {"luttrell": float(dL.max()),
                               "yang": float(dY.max())},
        "strike_span_km": {"SE": float(s.min()), "NW": float(s.max())},
        "s_arrest": args.s_arrest,
        "agreement": {nm: agree(m) for nm, m in
                      (("all", np.ones(n_facets, bool)),
                       ("seismogenic", seis), ("large_S_arrest", largeS))},
        "targeted_windows": win_summary,
        "regions": {}}
    for name, m in (("all", np.ones(n_facets, bool)),
                    ("seismogenic", seis), ("large_S_arrest", largeS),
                    ("propagating", propag), ("pure_barrier", pure_barrier)):
        summary["regions"][name] = {
            "n": int(m.sum()),
            "R_luttrell": field_stats(R_lut[m]),
            "R_yang": field_stats(R_yng[m]),
            "dR_L_minus_Y": field_stats(dR[m])}
    os.makedirs(args.outdir, exist_ok=True)
    jpath = os.path.join(args.outdir,
                         "compare_R_luttrell_vs_yang_summary.json")
    with open(jpath, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {jpath}")

    # --- plots ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zk = z / 1000.0
    # large-S strike extent for shading
    def shade_bands(ax, mask, color, label):
        if mask.sum() == 0:
            return
        ss = np.sort(s[mask])
        # crude contiguous-run shading
        gaps = np.where(np.diff(ss) > 3.0)[0]
        starts = np.concatenate([[0], gaps + 1])
        ends = np.concatenate([gaps, [len(ss) - 1]])
        first = True
        for a, b in zip(starts, ends):
            ax.axvspan(ss[a], ss[b], color=color, alpha=0.12,
                       label=label if first else None)
            first = False

    # Figure 1: strike-depth scatter, both models + difference
    fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True, sharey=True)
    sc0 = axes[0].scatter(s, zk, c=R_yng, s=3, cmap="viridis",
                          vmin=0.0, vmax=1.0)
    axes[0].set_title("R = (S2-S3)/(S1-S3)  --  YHSM-2013 (Yang & Hauksson) "
                      "[model V3 uses]")
    fig.colorbar(sc0, ax=axes[0], label="R")
    sc1 = axes[1].scatter(s, zk, c=R_lut, s=3, cmap="viridis",
                          vmin=0.0, vmax=1.0)
    axes[1].set_title("R = (S2-S3)/(S1-S3)  --  Luttrell-2017 "
                      "(differential stress)")
    fig.colorbar(sc1, ax=axes[1], label="R")
    vmax = np.nanpercentile(np.abs(dR), 98)
    sc2 = axes[2].scatter(s, zk, c=dR, s=3, cmap="RdBu_r",
                          vmin=-vmax, vmax=vmax)
    axes[2].set_title(f"dR = R_Luttrell - R_Yang   "
                      f"(red: Luttrell more strike-slip / higher R)")
    fig.colorbar(sc2, ax=axes[2], label="dR")
    for ax in axes:
        shade_bands(ax, largeS, "red", f"S>={args.s_arrest} (arrest)")
        ax.axhline(args.seis_ztop / 1000.0, ls=":", c="grey", lw=0.8)
        ax.axhline(args.seis_zbot / 1000.0, ls=":", c="grey", lw=0.8)
        ax.set_ylabel("depth (km)")
        ax.legend(loc="lower right", fontsize=8)
    axes[-1].set_xlabel("strike distance s (km)  [SE - , hypocenter 0, NW +]")
    fig.tight_layout()
    p1 = os.path.join(args.outdir,
                      "compare_R_luttrell_vs_yang_strike_depth.png")
    fig.savefig(p1, dpi=110)
    print(f"wrote {p1}")

    # Figure 2: along-strike 1D, seismogenic-band binned medians
    fig2, ax2 = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    sbin = np.arange(np.floor(s.min()), np.ceil(s.max()) + 2, 2.0)
    mid = 0.5 * (sbin[:-1] + sbin[1:])

    def binned(vals, mask):
        out = np.full(len(mid), np.nan)
        for i in range(len(mid)):
            sel = mask & (s >= sbin[i]) & (s < sbin[i + 1]) & np.isfinite(vals)
            if sel.any():
                out[i] = np.median(vals[sel])
        return out

    ax2[0].plot(mid, binned(R_yng, seis), "-", color="C0",
                label="R Yang (YHSM-2013)")
    ax2[0].plot(mid, binned(R_lut, seis), "-", color="C3",
                label="R Luttrell-2017")
    ax2[0].set_ylabel("R (seismogenic-band median)")
    ax2[0].set_ylim(0, 1)
    ax2[0].legend(loc="upper right")
    ax2[0].set_title("Stress-shape ratio R along strike (median over the "
                     "seismogenic band)")
    ax2[1].plot(mid, binned(dR, seis), "-", color="purple",
                label="dR = Luttrell - Yang")
    ax2[1].axhline(0, color="k", lw=0.6)
    ax2[1].set_ylabel("dR (median)")
    ax2[1].set_xlabel("strike distance s (km)  [SE - , hypocenter 0, NW +]")
    ax2[1].legend(loc="upper right")
    for ax in ax2:
        shade_bands(ax, largeS, "red", f"S>={args.s_arrest} (arrest)")
        ax.legend(loc="upper right", fontsize=8)
    fig2.tight_layout()
    p2 = os.path.join(args.outdir,
                      "compare_R_luttrell_vs_yang_along_strike.png")
    fig2.savefig(p2, dpi=110)
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
