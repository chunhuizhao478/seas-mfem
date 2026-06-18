#!/usr/bin/env python3
"""gate_removal_sweep.py — feasibility sweep for removing the SAFS CSM arrest
gates (S-ratio S >= S_TARGET) by COMBINING two levers:

  Lever 1 (closure k): lower the regional differential-stress ratio k so the
    GLOBAL apparent-friction ceiling drops, which is the only thing that lets
    the constant static friction mu_s drop without t=0 pre-slip.
  Lever 2 (fault-local overpressure): subtract an isotropic, fault-local excess
    pore pressure DeltaPp from the EFFECTIVE stress tensor at the gate, raising
    mu_app = tau/sigma_n_eff there (Fulton & Saffer 2009 overpressured
    restraining-bend core), capped sub-lithostatic (lambda = Pp/Sv <= LAMBDA_MAX).

This tool is PURE ANALYSIS: it re-projects the production CSM C1 field (closure
ratio, hydrostatic Pp, MUSCAL Sv, axes=csm orientation, linear sampling) once
per closure ratio k, applies the graded fault-local overpressure model, and
writes a markdown report comparing recipes B (k=1.8) and C (k=1.5) across all
four gates.  It modifies NO production artifact (no .nc, no .yaml, no VTU).

Physics is reused VERBATIM from project_csm_stress_to_vtu.py /
project_csm_on_deepmesh.py; only the closure k and the (post-closure) isotropic
overpressure are varied here.  See PLAN_gate_removal_csm_2026-06-17.md.

Plain-text math only:
  floor(mu_s)  = (mu_s + S_TARGET*mu_d) / (1 + S_TARGET)   -> mu_app >= floor <=> S < S_TARGET
  DeltaPp_req  = sigma_n_eff - tau/mu_target               (>=0 when mu_app < mu_target)
  cap_tension  = sigma_n_eff - SN_FLOOR                    (SeisSol anti-tension clamp)
  cap_litho    = LAMBDA_MAX*Sv_total - Pp_hydro            (sub-lithostatic Pp cap)
  mu_app(DeltaPp) = tau / (sigma_n_eff - DeltaPp)          (tau unchanged: Pp isotropic)

Run under the env with scipy+pyproj+netCDF4 (pythonenv):
  /Users/chunhuizhao/miniforge/envs/pythonenv/bin/python gate_removal_sweep.py
"""
import argparse
import csv as csv_mod
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_csm_stress_to_vtu import (  # noqa: E402
    read_csm_csv, csm_tensors_tension, interpolate_csm_field,
    geographic_to_utm11n, muscal_sv_total_profile, sv_total_at, pore_pressure,
    magnitudes_C1, build_tensor_from_axes, triangle_geometry,
    harmonise_normals, tandem_basis, resolve_tractions, load_puml,
    sanity_check_face_ordering, extract_fault, EPS)
from project_csm_on_deepmesh import read_gmsh_fault  # noqa: E402

# ----------------------------------------------------------- frozen constants
ALT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))  # project_7.0_alternative
# V3 LSW runs on the PUML mesh (parameters.par: MeshFile='safs_mesh.puml.h5');
# the deep .msh was only used to CHARACTERISE the gates.  Phase-0 decision:
# measure S on the mesh SeisSol actually reads -> PUML by default.
DEFAULT_MESH = os.path.join(
    ALT, "seisol_quakeworx/safs_seisol_v3_0_0_LSW/safs_mesh.puml.h5")
CSM_CSV = os.path.join(ALT, "seisol_quakeworx/raw_data/"
                       "yang_and_hauksson_orientation/CSM_data_1781290479811.csv")
MAT_NC = os.path.join(
    ALT, "seisol_quakeworx/safs_seisol_v3_0_0_LSW/safs_material_cvm.nc")
GMSH_FAULT_TAG = 103                         # only used for a .msh mesh

GATES_KM = {"SanGorgonio": (40.0, 71.0), "g181": (181.0, 186.0),
            "g221": (221.0, 246.0), "g266": (266.0, 271.0)}
SEIS_ZKM = (-12.0, -3.0)                     # seismogenic propagation band
ACTIVE_ZKM = (-15.0, -0.3)                   # subcriticality band (no surface
                                             # slivers, no deep locked band)
MU_D = 0.10
S_TARGET = 1.7
SN_FLOOR_MPA = 0.5                           # min sigma_n_eff after overpressure
STRIKE_AZ = 314.0
HYPO = np.array([606971.0, 3707270.0, -4965.62])
DEFAULT_K_LIST = [2.39, 2.2, 2.0, 1.8, 1.6, 1.5]
RECIPES = {"B": 1.8, "C": 1.5}              # name -> closure k
DEPTH_BIN_KM = 1.0
MIN_BIN_FACETS = 10


def floor_mu(mu_s):
    """mu_app needed for S < S_TARGET (S = (mu_s-mu_app)/(mu_app-mu_d))."""
    return (mu_s + S_TARGET * MU_D) / (1.0 + S_TARGET)


def load_fault(mesh_path):
    """(pts, tris) for the fault, supporting PUML .h5 (BC-3) or Gmsh .msh
    (physical tag 103) — matches whatever mesh SeisSol reads."""
    if mesh_path.endswith(".msh"):
        return read_gmsh_fault(mesh_path, GMSH_FAULT_TAG)
    geom, conn, bc = load_puml(mesh_path)
    sanity_check_face_ordering(geom, conn, bc)
    return extract_fault(geom, conn, bc)


def pct(a, ps):
    if a.size == 0:
        return [float("nan")] * len(ps)
    return [float(np.nanpercentile(a, p)) for p in ps]


class Field:
    """k-independent fault geometry + CSM orientation/shape, projected once."""

    def __init__(self, mesh_path):
        pts, tris = load_fault(mesh_path)
        cent, nraw, areas = triangle_geometry(pts, tris)
        normals, _ = harmonise_normals(nraw, STRIKE_AZ)
        strikes, dips, degen = tandem_basis(normals)
        self.cent, self.normals, self.strikes, self.dips = \
            cent, normals, strikes, dips
        self.n_facets = len(tris)
        self.n_degen = int(degen.sum())
        depth = np.maximum(-cent[:, 2], 0.0)
        dg, sg, _ = muscal_sv_total_profile(MAT_NC)
        self.Sv_total = sv_total_at(depth, dg, sg)
        self.Pp_h = pore_pressure(depth, "hydrostatic")
        self.Sv_eff = np.maximum(self.Sv_total - self.Pp_h, EPS)
        lon, lat, _d, S, *_ = read_csm_csv(CSM_CSV)
        cx, cy = geographic_to_utm11n(lon, lat)
        S_facet, self.n_fb = interpolate_csm_field(
            cx, cy, S, cent[:, 0], cent[:, 1])
        T = csm_tensors_tension(S_facet)
        w, self.U = np.linalg.eigh(T)
        den = w[:, 2] - w[:, 0]
        self.R = np.clip((w[:, 1] - w[:, 0]) /
                         np.where(den > EPS, den, np.nan), 0.0, 1.0)
        su = np.array([np.sin(np.radians(STRIKE_AZ)),
                       np.cos(np.radians(STRIKE_AZ))])
        self.s_km = ((cent[:, :2] - HYPO[:2]) @ su) / 1000.0
        self.zkm = cent[:, 2] / 1000.0
        self.seis = (self.zkm > SEIS_ZKM[0]) & (self.zkm < SEIS_ZKM[1])

    def project(self, k):
        """Per-k mu_app, tau, sigma_n_eff (comp+ MPa) + the active mask."""
        sig1, sig2, sig3 = magnitudes_C1(self.Sv_eff, self.R, k)
        sigma = build_tensor_from_axes(self.U, sig1, sig2, sig3)
        r = resolve_tractions(sigma, self.strikes, self.dips,
                              self.normals, self.Pp_h)
        mu, tau, sn = r["mu_apparent"], r["tau_magnitude"], r["sigma_n_eff"]
        active = ((self.zkm < ACTIVE_ZKM[1]) & (self.zkm > ACTIVE_ZKM[0])
                  & np.isfinite(mu) & (sn > 1.0))
        return mu, tau, sn, active


def caps(field, sn, lambda_max):
    """Per-facet DeltaPp caps: tension (sigma_n_eff floor) and lithostatic."""
    cap_tension = sn - SN_FLOOR_MPA
    cap_litho = lambda_max * field.Sv_total - field.Pp_h
    return np.minimum(cap_tension, cap_litho)


def graded_gate(field, mu, tau, sn, active, gate_mask, mu_s, lambda_max):
    """Graded per-facet overpressure that lifts each gate seismogenic facet to
    floor(mu_s).  Returns a dict of metrics + the per-facet DeltaPp/mu_new."""
    floor = floor_mu(mu_s)
    band = gate_mask & field.seis & active
    need = band & (mu < floor)
    cap = caps(field, sn, lambda_max)
    DPp_req = np.where(need, sn - tau / floor, 0.0)
    unfix = need & (DPp_req > cap)
    DPp = np.clip(np.minimum(DPp_req, cap), 0.0, None)
    mu_new = np.where(DPp > 0.0, tau / np.maximum(sn - DPp, EPS), mu)
    cross = band & np.isfinite(mu_new) & (mu_new >= floor - 1e-9)
    lam_after = np.where(band, (field.Pp_h + DPp) / np.maximum(field.Sv_total,
                         EPS), np.nan)
    lifted = need & ~unfix
    return {
        "floor": floor, "n_band": int(band.sum()), "n_need": int(need.sum()),
        "n_unfix": int(unfix.sum()), "unfix_mask": unfix,
        "cross_frac": float(cross[band].mean()) if band.any() else float("nan"),
        "DPp_lifted": DPp[lifted], "lam_band": lam_after[band],
        "mu_new": mu_new, "band": band, "DPp": DPp,
    }


def depth_envelope(field, mu, tau, sn, active, gate_mask, mu_s, lambda_max):
    """A production-realizable lambda_op(z) depth profile for one gate band:
    per 1 km bin, the just-enough DeltaPp_need(z), and whether it stays under
    the subcritical cap (no band facet pushed to mu_s) and the tension/litho
    caps.  Returns feasibility + the resulting %S<1.7 and max mu_app."""
    floor = floor_mu(mu_s)
    band = gate_mask & field.seis & active
    cap = caps(field, sn, lambda_max)
    edges = np.arange(SEIS_ZKM[0], SEIS_ZKM[1] + DEPTH_BIN_KM, DEPTH_BIN_KM)
    DPp_profile = np.zeros(field.n_facets)     # per-facet DeltaPp from the envelope
    feasible = True
    bins_used = 0
    rows = []
    for zlo, zhi in zip(edges[:-1], edges[1:]):
        sel = band & (field.zkm >= zlo) & (field.zkm < zhi)
        nsel = int(sel.sum())
        if nsel < MIN_BIN_FACETS:
            continue
        bins_used += 1
        dpp_need = float(np.max(np.clip(sn[sel] - tau[sel] / floor, 0.0, None)))
        dpp_sub = float(np.min(np.clip(sn[sel] - tau[sel] / mu_s, 0.0, None)))
        dpp_cap = float(np.min(cap[sel]))
        ok = dpp_need <= min(dpp_sub, dpp_cap)
        feasible = feasible and ok
        DPp_profile[sel] = min(dpp_need, dpp_sub, dpp_cap) if ok else dpp_need
        rows.append((0.5 * (zlo + zhi), nsel, dpp_need, dpp_sub, dpp_cap, ok))
    mu_env = np.where((DPp_profile > 0) & band,
                      tau / np.maximum(sn - DPp_profile, EPS), mu)
    cross = band & np.isfinite(mu_env) & (mu_env >= floor - 1e-9)
    # subcriticality of the WHOLE active set after the envelope bake
    mu_full = np.where((DPp_profile > 0) & active,
                       tau / np.maximum(sn - DPp_profile, EPS), mu)
    max_mu_after = float(np.nanmax(mu_full[active])) if active.any() else float("nan")
    return {
        "feasible": bool(feasible and bins_used > 0),
        "bins_used": bins_used, "rows": rows,
        "cross_frac": float(cross[band].mean()) if band.any() else float("nan"),
        "max_mu_after": max_mu_after, "subcritical": max_mu_after < mu_s,
    }


def mu_s_for_k(field, k):
    """mu_s = max(mu_app over active) + 0.02, rounded UP to 3 decimals."""
    mu, _, _, active = field.project(k)
    mx = float(np.nanmax(mu[active])) if active.any() else float("nan")
    return float(np.ceil((mx + 0.02) * 1000.0) / 1000.0), mx


def build_report(field, k_list, recipes, lambda_max, mesh_path):
    L = []
    w = L.append
    w("# Gate-removal feasibility sweep (CSM C1) — recipes B vs C\n")
    w(f"_Generated by gate_removal_sweep.py on the V3 production mesh._\n")
    w("")
    w("## Setup")
    w(f"- Mesh: `{os.path.relpath(mesh_path, ALT)}` "
      f"({field.n_facets} fault facets, {field.n_degen} degenerate-basis, "
      f"{field.n_fb} CSM-hull NN-fallback).")
    w(f"- Physics: production CSM C1 (closure ratio, hydrostatic Pp, MUSCAL "
      f"Sv, axes=csm, linear sampling); only k and a post-closure isotropic "
      f"fault-local overpressure are varied.")
    w(f"- mu_d = {MU_D} (fixed); S_TARGET = {S_TARGET}; "
      f"LAMBDA_MAX = {lambda_max} (Pp/Sv cap); SN_FLOOR = {SN_FLOOR_MPA} MPa.")
    w(f"- Seismogenic band z in {SEIS_ZKM} km; active (subcriticality) band "
      f"z in {ACTIVE_ZKM} km; strike origin = hypocenter, NW positive.")
    w(f"- Gate bands (km): " +
      ", ".join(f"{n} {v}" for n, v in GATES_KM.items()) + ".")
    w("")
    w("S < 1.7 per facet <=> mu_app >= floor(mu_s) = (mu_s + 1.7*mu_d)/2.7.")
    w("Overpressure raises mu_app via DeltaPp = sigma_n_eff - tau/mu_target "
      "(tau unchanged; Pp isotropic), capped by tension and lambda<=" +
      f"{lambda_max}.")
    w("")

    # ---- k-ceiling table -------------------------------------------------
    w("## Lever k: the global mu_app ceiling (pins mu_s)")
    w("`mu_s` must exceed the fault-wide max mu_app or well-oriented facets "
      "pre-slip at t=0. Lowering k lowers that ceiling.\n")
    w("| k | median mu_app | max mu_app (active) | min feasible mu_s (max+0.02) |")
    w("|---|---|---|---|")
    krows = []
    for k in k_list:
        mu, _, _, active = field.project(k)
        md = float(np.nanmedian(mu[active]))
        mx = float(np.nanmax(mu[active]))
        ms = float(np.ceil((mx + 0.02) * 1000.0) / 1000.0)
        krows.append((k, md, mx, ms))
        w(f"| {k:.2f} | {md:.3f} | {mx:.4f} | {ms:.3f} |")
    w("")

    # ---- per-recipe per-gate analysis -----------------------------------
    csv_rows = []
    verdicts = {}
    for rname, k in recipes.items():
        mu, tau, sn, active = field.project(k)
        mu_s, mx = mu_s_for_k(field, k)
        w(f"## Recipe {rname}: k = {k:.2f}, mu_s = {mu_s:.3f} "
          f"(max mu_app {mx:.4f}), floor = {floor_mu(mu_s):.3f}")
        w("| gate | n_seis | need lift | UNFIXABLE | %S<1.7 graded | "
          "%S<1.7 envelope | DPp med/p90/max MPa | lambda med/max | "
          "env-feasible | max mu_app after | subcritical |")
        w("|---|---|---|---|---|---|---|---|---|---|---|")
        removed_graded, removed_env = [], []
        for gname, (a, b) in GATES_KM.items():
            gmask = (field.s_km >= a) & (field.s_km < b)
            g = graded_gate(field, mu, tau, sn, active, gmask, mu_s, lambda_max)
            e = depth_envelope(field, mu, tau, sn, active, gmask, mu_s,
                               lambda_max)
            d50, d90, d100 = pct(g["DPp_lifted"], (50, 90, 100))
            l50, l100 = pct(g["lam_band"], (50, 100))
            w(f"| {gname} | {g['n_band']} | {g['n_need']} | {g['n_unfix']} | "
              f"{100*g['cross_frac']:.1f}% | {100*e['cross_frac']:.1f}% | "
              f"{d50:.1f}/{d90:.1f}/{d100:.1f} | "
              f"{l50:.2f}/{l100:.2f} | {'Y' if e['feasible'] else 'N'} | "
              f"{e['max_mu_after']:.3f} | {'Y' if e['subcritical'] else 'N'} |")
            if g["cross_frac"] >= 0.95 and g["n_unfix"] == 0:
                removed_graded.append(gname)
            if e["feasible"] and e["cross_frac"] >= 0.95 and e["subcritical"]:
                removed_env.append(gname)
            # barrier cores list
            if g["n_unfix"]:
                idx = np.where(g["unfix_mask"])[0]
                pts_ = [(round(float(field.s_km[i]), 1),
                         round(float(field.zkm[i]), 1)) for i in idx[:8]]
                w(f"  - {gname}: {g['n_unfix']} UNFIXABLE barrier cores "
                  f"(s,z km, first 8): {pts_}")
            csv_rows.append([rname, k, mu_s, gname, g["n_band"], g["n_need"],
                             g["n_unfix"], round(100*g["cross_frac"], 1),
                             round(100*e["cross_frac"], 1), e["feasible"],
                             round(e["max_mu_after"], 4), e["subcritical"]])
        w("")
        verdicts[rname] = (k, mu_s, removed_graded, removed_env)
        w(f"**Recipe {rname} verdict:** graded removes "
          f"{removed_graded if removed_graded else 'NONE'}; depth-envelope "
          f"removes {removed_env if removed_env else 'NONE'} "
          f"(>=95% facets S<1.7, subcritical, lambda<={lambda_max}).")
        w("")

    # ---- recommendation --------------------------------------------------
    w("## Recommendation (B vs C)")
    kb, msb, gb, eb = verdicts["B"]
    kc, msc, gc, ec = verdicts["C"]
    w(f"- Recipe B (k={kb:.2f}, mu_s={msb:.3f}): graded removes {gb}; "
      f"envelope removes {eb}.")
    w(f"- Recipe C (k={kc:.2f}, mu_s={msc:.3f}): graded removes {gc}; "
      f"envelope removes {ec}.")
    w("- Trade-off: B keeps k nearer the Cajon Pass strong-crust value "
      "(Zoback & Healy 1992, k~3 for mu~0.6) and needs MORE overpressure; "
      "C weakens the regional field below that constraint but needs the LEAST "
      "overpressure. Both keep mu_s constant (no facet-varying friction). "
      "Pick by how much regional-stress weakening vs how much fault-zone "
      "overpressure is acceptable; the tables above quantify both.")
    w("")
    return "\n".join(L), krows, csv_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mesh", default=DEFAULT_MESH,
                    help="fault mesh (PUML .h5 [V3 default] or Gmsh .msh)")
    ap.add_argument("--out", default=os.path.join(
        HERE, "GATE_REMOVAL_SWEEP_2026-06-17.md"))
    ap.add_argument("--csv-out", default=os.path.join(
        HERE, "GATE_REMOVAL_SWEEP_2026-06-17.csv"))
    ap.add_argument("--k-list", default=",".join(str(k) for k in DEFAULT_K_LIST))
    ap.add_argument("--recipes", default="B,C")
    ap.add_argument("--lambda-max", type=float, default=0.9)
    args = ap.parse_args()

    k_list = [float(x) for x in args.k_list.split(",")]
    recipes = {r: RECIPES[r] for r in args.recipes.split(",")}
    print(f"projecting on {args.mesh}")
    field = Field(args.mesh)
    print(f"facets {field.n_facets}  degenerate {field.n_degen}  "
          f"NN-fallback {field.n_fb}")

    report, krows, csv_rows = build_report(
        field, k_list, recipes, args.lambda_max, args.mesh)
    with open(args.out, "w") as fh:
        fh.write(report)
    with open(args.csv_out, "w", newline="") as fh:
        wcsv = csv_mod.writer(fh)
        wcsv.writerow(["recipe", "k", "mu_s", "gate", "n_seis", "n_need",
                       "n_unfix", "pct_S<1.7_graded", "pct_S<1.7_env",
                       "env_feasible", "max_mu_after", "subcritical"])
        wcsv.writerows(csv_rows)
    print(f"wrote {args.out}")
    print(f"wrote {args.csv_out}")

    # ---- acceptance echo -------------------------------------------------
    mx239 = next((r[2] for r in krows if abs(r[0] - 2.39) < 1e-9), None)
    print(f"\nACCEPTANCE CHECK:")
    if mx239 is not None:
        ok239 = abs(mx239 - 0.449) <= 0.001
        print(f"  max mu_app(k=2.39) = {mx239:.4f}  "
              f"(expect 0.449 +/- 0.001: {'PASS' if ok239 else 'CHECK'})")
    print(f"  (B/C San Gorgonio unfixable counts are in the report tables)")


if __name__ == "__main__":
    main()
