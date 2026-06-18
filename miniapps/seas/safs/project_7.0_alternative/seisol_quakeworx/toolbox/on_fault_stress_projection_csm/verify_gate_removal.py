#!/usr/bin/env python3
"""verify_gate_removal.py — end-to-end gate-removal check (PLAN Phase 4).

Samples the BAKED nc at the V3 PUML fault facet centroids EXACTLY as SeisSol
will (trilinear, BaseDRInitializer rotateStressToFaultCS is just n.sigma.n),
computes mu_app and the seismic ratio S per facet, and compares BEFORE
(production nc, mu_s_before) vs AFTER (baked recipe nc, mu_s_after) per gate:
  - %facets with S < S_TARGET in the seismogenic band,
  - pure-barrier count (mu_app < mu_d),
  - subcriticality (max mu_app over active set < mu_s_after),
  - tension (min sigma_n_eff over all facets > 0).
Writes a markdown result + a before/after strike-depth S figure.  This is the
real test of what SeisSol receives; S<S_TARGET is the analytic propagation
proxy (a SeisSol dynamic-rupture run on Frontera is the ultimate confirmation).
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_csm_stress_to_vtu import (  # noqa: E402
    load_puml, sanity_check_face_ordering, extract_fault, triangle_geometry,
    harmonise_normals, tandem_basis, resolve_tractions)
from csm_stress_to_asagi import trilinear_sample  # noqa: E402
from gate_removal_sweep import (  # noqa: E402
    GATES_KM, SEIS_ZKM, ACTIVE_ZKM, MU_D, S_TARGET, STRIKE_AZ, HYPO,
    DEFAULT_MESH)

CIDX = {"s_xx": (0, 0), "s_yy": (1, 1), "s_zz": (2, 2),
        "s_xy": (0, 1), "s_yz": (1, 2), "s_xz": (0, 2)}


def sample_nc_on_fault(ncpath, cent, strikes, dips, normals):
    """Sample the comp-negative-Pa nc at facet centroids, back-convert to comp+
    MPa, resolve tractions (P_p baked in -> P_p=0)."""
    smp = trilinear_sample(ncpath, cent[:, 0], cent[:, 1], cent[:, 2])
    sig = np.zeros((len(cent), 3, 3))
    for f, (a, b) in CIDX.items():
        sig[:, a, b] = -smp[f] / 1.0e6
        sig[:, b, a] = sig[:, a, b]
    return resolve_tractions(sig, strikes, dips, normals, 0.0)


def s_ratio(mu, mu_s, mu_d):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(mu > mu_d, (mu_s - mu) / (mu - mu_d), np.inf)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mesh", default=DEFAULT_MESH)
    ap.add_argument("--nc-before",
                    default=os.path.join(HERE, "..", "..",
                    "safs_seisol_v3_0_0_LSW", "safs_stress_csm.nc"))
    ap.add_argument("--nc-after",
                    default=os.path.join(HERE, "..", "..",
                    "safs_seisol_v3_0_0_LSW", "safs_stress_csm_B.nc"))
    ap.add_argument("--mu-s-before", type=float, default=0.47)
    ap.add_argument("--mu-s-after", type=float, default=0.318)
    ap.add_argument("--recipe", default="B")
    ap.add_argument("--out", default=os.path.join(
        HERE, "GATE_REMOVAL_RESULT_B_2026-06-17.md"))
    ap.add_argument("--fig", default=os.path.join(
        HERE, "gate_removal_B_strike_depth.png"))
    args = ap.parse_args()

    geom, conn, bc = load_puml(args.mesh)
    sanity_check_face_ordering(geom, conn, bc)
    pts, tris = extract_fault(geom, conn, bc)
    cent, nraw, areas = triangle_geometry(pts, tris)
    normals, _ = harmonise_normals(nraw, STRIKE_AZ)
    strikes, dips, degen = tandem_basis(normals)
    su = np.array([np.sin(np.radians(STRIKE_AZ)), np.cos(np.radians(STRIKE_AZ))])
    s_km = ((cent[:, :2] - HYPO[:2]) @ su) / 1000.0
    zkm = cent[:, 2] / 1000.0
    seis = (zkm > SEIS_ZKM[0]) & (zkm < SEIS_ZKM[1])
    active = (zkm < ACTIVE_ZKM[1]) & (zkm > ACTIVE_ZKM[0])

    rb = sample_nc_on_fault(args.nc_before, cent, strikes, dips, normals)
    ra = sample_nc_on_fault(args.nc_after, cent, strikes, dips, normals)
    mub, mua = rb["mu_apparent"], ra["mu_apparent"]
    finb = np.isfinite(mub) & (rb["sigma_n_eff"] > 0)
    fina = np.isfinite(mua) & (ra["sigma_n_eff"] > 0)
    Sb = s_ratio(mub, args.mu_s_before, MU_D)
    Sa = s_ratio(mua, args.mu_s_after, MU_D)

    # ---- global checks on the AFTER field ----
    act_a = active & fina
    max_mu_after = float(np.nanmax(mua[act_a]))
    subcrit = max_mu_after < args.mu_s_after
    min_sn_after = float(np.nanmin(ra["sigma_n_eff"][fina]))
    min_sn_seis_after = float(np.nanmin(ra["sigma_n_eff"][fina & seis]))

    L = []
    w = L.append
    w(f"# Gate-removal result — recipe {args.recipe} (end-to-end, baked nc)\n")
    w(f"Sampled the baked nc at {len(tris)} PUML fault facet centroids (as "
      f"SeisSol does) and computed S = (mu_s - mu_app)/(mu_app - mu_d).\n")
    w(f"- BEFORE: `{os.path.basename(args.nc_before)}` with mu_s={args.mu_s_before}"
      f" (production V3).")
    w(f"- AFTER : `{os.path.basename(args.nc_after)}` with mu_s={args.mu_s_after}"
      f" (recipe {args.recipe}: lower k + fault-local overpressure).")
    w(f"- mu_d={MU_D}, S_TARGET={S_TARGET}, seismogenic z in {SEIS_ZKM} km.\n")
    w("## Per-gate seismogenic S < 1.7 (before -> after)")
    w("| gate | n_seis | mu_app med B->A | pure-barriers B->A | %S<1.7 B->A | "
      "residual S>=1.7 facets (after) |")
    w("|---|---|---|---|---|---|")
    all_removed = True
    rows = []
    for gname, (a, b) in GATES_KM.items():
        g = (s_km >= a) & (s_km < b) & seis
        nb = int((g & finb).sum())
        if nb == 0:
            w(f"| {gname} | 0 | - | - | - | (not on PUML mesh) |")
            continue
        ga, gb = g & fina, g & finb
        barr_b = int((mub[gb] < MU_D).sum())
        barr_a = int((mua[ga] < MU_D).sum())
        pctb = 100.0 * np.mean(Sb[gb] < S_TARGET)
        pcta = 100.0 * np.mean(Sa[ga] < S_TARGET)
        resid = int((Sa[ga] >= S_TARGET).sum())
        if pcta < 95.0:
            all_removed = False
        w(f"| {gname} | {int(ga.sum())} | {np.nanmedian(mub[gb]):.3f}->"
          f"{np.nanmedian(mua[ga]):.3f} | {barr_b}->{barr_a} | "
          f"{pctb:.0f}%->{pcta:.0f}% | {resid} |")
        rows.append((gname, a, b, pcta, resid))
    w("")
    w("## Global checks on the AFTER field")
    w(f"- Subcriticality: max mu_app (active z in {ACTIVE_ZKM} km) = "
      f"{max_mu_after:.4f}  {'<' if subcrit else '>='}  mu_s={args.mu_s_after}"
      f"  -> {'PASS (no t=0 pre-slip)' if subcrit else 'FAIL'}.")
    w(f"- Tension: min sigma_n_eff (all facets) = {min_sn_after:.2f} MPa; "
      f"min in seismogenic band = {min_sn_seis_after:.2f} MPa  -> "
      f"{'PASS (>0)' if min_sn_after > 0 else 'FAIL (<=0, SeisSol clamp)'}.")
    w(f"- Overall: {'ALL TARGETED GATES REMOVED (>=95% S<1.7)' if all_removed else 'SOME GATES NOT FULLY REMOVED — see residual column'}.")
    w("")
    w("NOTE: S<1.7 is the analytic propagation proxy (Das & Aki). A SeisSol "
      "dynamic-rupture run on Frontera is the ultimate confirmation that the "
      "rupture crosses; this static check is necessary, not sufficient.\n")
    with open(args.out, "w") as fh:
        fh.write("\n".join(L))
    print("\n".join(L))
    print(f"\nwrote {args.out}")

    # ---- before/after figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        fig, axes = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
        for ax, (Sv, mus, tag, fin) in zip(axes, [
                (Sb, args.mu_s_before, "BEFORE (production, mu_s=%.2f)" %
                 args.mu_s_before, finb),
                (Sa, args.mu_s_after, "AFTER recipe %s (mu_s=%.3f + local "
                 "overpressure)" % (args.recipe, args.mu_s_after), fina)]):
            Scl = np.clip(np.where(np.isfinite(Sv), Sv, 99), 0, 3.5)
            order = np.argsort(zkm)
            sc = ax.scatter(s_km[order], zkm[order], c=Scl[order], s=3,
                            cmap="RdYlGn_r", vmin=0, vmax=3.5, linewidths=0)
            ax.set_title(f"S-ratio  {tag}  (GREEN S<1.7 propagates, RED arrests)")
            ax.set_ylabel("depth (km)")
            for gname, (a, b) in GATES_KM.items():
                if b > s_km.max():
                    continue
                ax.add_patch(Rectangle((a, zkm.min()), b - a,
                             zkm.max() - zkm.min(), fill=False,
                             edgecolor="blue", lw=1.5, ls="--"))
            fig.colorbar(sc, ax=ax, label="S-ratio (clip 3.5)", pad=0.01)
        axes[1].set_xlabel("strike distance s (km)  [SE ~ -18, hypocenter 0, "
                           "NW ~ 266]")
        fig.tight_layout()
        fig.savefig(args.fig, dpi=130)
        print(f"wrote {args.fig}")
    except Exception as e:
        print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main()
