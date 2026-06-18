# Implementation Plan: Remove SAFS CSM arrest gates (S < 1.7) by combining a closure-k stress reduction with fault-local graded overpressure

> NOTE ON LOCATION: the code-plan skill nominally writes `PLAN.md` to the repo
> root.  This repo's convention is dated `PLAN_*.md` documents colocated with
> the toolbox they govern (e.g. `PLAN_csm_stress_magnitudes_2026-06-12.md`), so
> this plan lives in the toolbox folder instead.  `/code-implement` should be
> pointed at THIS file.

## Overview
A SE-nucleated linear-slip-weakening rupture on the SAFS CSM C1 stress field
arrests at the San Gorgonio gate (strike s = 40-71 km) and at three more
seismogenic gates (s ~ 181-186, 221-246, 266-271 km) because the apparent
friction `mu_app = tau/sigma_n_eff` collapses toward `mu_d` there, driving the
seismic ratio `S = (mu_s - mu_app)/(mu_app - mu_d)` far above the propagation
limit (~1.7).  This plan removes those gates (target `S < 1.7` per facet in the
seismogenic band) by COMBINING two levers, then bakes the result into the
SeisSol-facing artifacts:

- Lever 1 (closure k): lower the regional differential-stress ratio `k` so the
  GLOBAL `mu_app` ceiling drops, which is the only thing that lets the constant
  `mu_s` drop without t=0 pre-slip.
- Lever 2 (fault-local overpressure): subtract an isotropic, fault-local excess
  pore pressure `DeltaPp(s,z)` from the EFFECTIVE stress tensor at the gate,
  raising `mu_app` there (Fulton & Saffer 2009 overpressured restraining-bend
  core), capped sub-lithostatic (Hubbert-Rubey `lambda = Pp/Sv <= 0.9`).

Two recipes are carried in parallel and compared on real numbers before any
production artifact is regenerated:
- Recipe B (moderate): k ~ 1.8, mu_s ~ 0.30.
- Recipe C (aggressive): k ~ 1.5, mu_s ~ 0.22 ("global mu_app ~ 0.2").

## Background facts already established (do NOT re-derive)
All verified on the DEEP mesh (132,333 fault-tag-103 facets) with the
production C1 projection (closure=ratio, hydrostatic Pp, MUSCAL Sv, axes=csm,
sampling=linear); see `/tmp/gate_diag.py`, `/tmp/gate_diag2.py`.

1. `mu_app` is SCALE-INVARIANT in `Sv_eff`: in C1, sig2 = Sv_eff,
   sig3 = Sv_eff/((1-R)k+R), sig1 = k*sig3, so the whole EFFECTIVE tensor is
   proportional to Sv_eff and `mu_app = tau/sigma_n_eff` does not depend on the
   Sv_eff magnitude.  => Reducing Sv_eff by adding pore pressure is a NO-OP on
   mu_app.  Lever 2 must instead subtract DeltaPp*I from the tensor AFTER the
   closure, which is isotropic and therefore changes sigma_n_eff but not tau.

2. Subcriticality is HARD: `mu_s` must exceed the fault-wide max `mu_app`
   (active set z in (-15 km, -0.3 km), sigma_n_eff > 1 MPa) or those
   well-oriented facets slip at t=0.  At k=2.39 the active max is 0.4493, which
   is exactly why production pins mu_s = 0.47.  Lever-k controls this ceiling:
   ```
   k     median mu_app   max mu_app(active)   min feasible mu_s (= max + 0.02)
   2.39    0.296           0.4493              0.469
   2.20    0.274           0.4044              0.424
   2.00    0.246           0.3536              0.374
   1.80    0.214           0.2981              0.318
   1.60    0.176           0.2371              0.257
   1.50    0.154           0.2041              0.224
   ```

3. GLOBAL overpressure is counterproductive: the per-MPa fractional gain
   d(ln mu_app)/dPp = 1/sigma_n_eff is LARGEST at the low-sigma_n_eff
   well-oriented facets, so a whole-fault overpressure raises the ceiling
   faster than the clamped gate and breaks subcriticality first.  Overpressure
   MUST be fault-local (concentrated where mu_app < floor).

4. S < 1.7 requires, with mu_d = 0.10 fixed:
   ```
   floor(mu_s) = (mu_s + 1.7*mu_d) / (1 + 1.7) = (mu_s + 0.17) / 2.7
   mu_s=0.47 -> floor 0.237 ;  mu_s=0.318 -> floor 0.181 ;  mu_s=0.224 -> floor 0.146
   ```
   A facet crosses iff `mu_app >= floor(mu_s)`.

5. To raise a facet from `mu_app = tau/sigma_n_eff` to a target `mu_t`, the
   required isotropic overpressure is
   ```
   DeltaPp = sigma_n_eff - tau/mu_t       (>= 0 when mu_app < mu_t)
   ```
   Caps (all must hold):
   - tension cap   : DeltaPp <= sigma_n_eff - SN_FLOOR   (SN_FLOOR = 0.5 MPa;
     SeisSol hardcodes the anti-tension clamp min(0, sigma_n) -> NaN past it).
   - lithostatic cap: DeltaPp <= (LAMBDA_MAX - lambda_hydro(z)) * Sv_total(z),
     LAMBDA_MAX = 0.9, lambda_hydro = Pp_hydro/Sv_total ~ 0.37-0.40.
   A facet is "unfixable" if `sigma_n_eff - tau/floor` exceeds the binding cap.

6. Diagnostic feasibility result (graded, per-facet, San Gorgonio seismogenic
   band, mu_d=0.10):
   - Recipe A (mu_s=0.47, k=2.39): 0 unfixable but lambda up to 0.90 (median
     0.63) — exceeds Fulton & Saffer support; REJECTED in favor of B/C.
   - Recipe B (k=1.8, mu_s=0.30 region): max mu_app ~ 0.298 subcritical;
     floor 0.174; 0 unfixable; less overpressure.
   - Recipe C (k=1.5, mu_s=0.22 region): max mu_app ~ 0.204; floor ~0.146;
     0 unfixable; least overpressure.

## Sign and convention reference (load before touching code)
- SEAS/projection internal: compression POSITIVE, MPa; `sigma_n_eff = n.sigma.n`
  from the EFFECTIVE tensor; `tau = |shear of sigma.n|`; `mu_app = tau/sigma_n_eff`.
- nc / SeisSol: compression NEGATIVE, Pa; stored diag `s_xx = -sigma_xx_comp+ * 1e6`.
- Overpressure baking (isotropic, raises mu_app):
  - projection (comp+ MPa): `sigma_facet[:,i,i] -= DeltaPp` for i in {0,1,2}.
  - nc (comp- Pa): `comps[s_xx|s_yy|s_zz] += DeltaPp_MPa * 1e6`  (less negative
    = less compressive); off-diagonals (s_xy,s_yz,s_xz) UNCHANGED so tau is exact.
- Strike coordinate (must match the projection): with STRIKE_AZ = 314 deg,
  `su = (sin(az), cos(az))`, `s_km = ((x,y) - HYPO[:2]) . su / 1000`,
  `HYPO = (606971, 3707270, -4965.62)`, NW positive.

## Constraints
- Interface: reuse the existing projection functions in
  `project_csm_stress_to_vtu.py` verbatim (read_csm_csv, csm_tensors_tension,
  interpolate_csm_field, csm_axes_and_shape, magnitudes_C1,
  build_tensor_from_axes, tandem_basis, harmonise_normals, resolve_tractions,
  muscal_sv_total_profile, sv_total_at, pore_pressure, write_vtu).  Do NOT fork
  the physics.
- Byte-identical default: every modification to a SHARED producer
  (`project_csm_stress_to_vtu.py`, `project_csm_on_deepmesh.py`,
  `csm_stress_to_asagi.py`) MUST default the overpressure OFF and reproduce the
  current artifacts bit-for-bit when off (regression contract; existing
  `*_summary.json`, `safs_stress_csm.nc` unchanged).
- mu_d = 0.10, d_c = 2.5 fixed.  Deep locked band z in (-20000,-15000) m
  (mu_s = 1e6) unchanged.
- SeisSol anti-tension clamp: sigma_n_eff must stay > 0.5 MPa on every facet
  after baking.  Memory: `seissol-normal-stress-clamp-and-rssrw-deep-nan`.
- Local-only verification: run the Python toolbox in
  `/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python`.  NEVER run a full
  SAFS SeisSol mesh locally (Frontera only) — validate by re-projection +
  nc self-check.  Memory: `feedback-no-local-mesh-runs`.
- SAF-specific references only (Fulton & Saffer 2009; Lachenbruch & Sass;
  Hardebeck & Hauksson; Zoback & Healy / Cajon Pass for k).  No imported
  parameters from other earthquakes.  Memory: `feedback-saf-references-only`.
- Plain-text math only in all docs/outputs (no LaTeX).  Memory:
  `feedback-no-latex-plain-text-math`.

## Phase 0: Confirm the production mesh and freeze the gate set

### Goal
The exact fault mesh the V3 LSW production run consumes is known, and the gate
bands + seismogenic window are frozen as named constants, so Phases 1-4 all
measure the same facets.

### Files to Create
- none (investigation only; record findings inline in Phase 1's tool header).

### Detailed Requirements
1. Determine which mesh the V3 run uses.  The nc self-check comment in
   `safs_initial_stress.yaml` cites "60,394 fault facet centroids" =
   `safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5` (PUML), while the gate analysis
   used `safs_seisol_v2_1_0_RSSRW/preprocess/safs_mesh_deep.msh` (132,333
   facets, Gmsh v2.2).  Resolve by inspecting the V3 case folder's
   `parameters.par` / mesh reference.  Record the answer.
   - DECISION RULE: the FEASIBILITY tool (Phase 1) and the final S-gate
     VERIFICATION (Phase 4) MUST run on whichever mesh SeisSol actually reads
     for V3.  If the deep mesh is the run mesh, use it everywhere; if the PUML
     mesh is, use it (it is shallower, ~16.5 km, still spans the gates).  The
     nc is mesh-independent; only the where-we-measure-S changes.
2. Freeze constants (used by all phases):
   ```
   GATES_KM = {"SanGorgonio": (40.0, 71.0), "g181": (181.0, 186.0),
               "g221": (221.0, 246.0), "g266": (266.0, 271.0)}
   SEIS_ZKM = (-12.0, -3.0)          # seismogenic propagation band
   ACTIVE_ZKM = (-15.0, -0.3)        # subcriticality test band (excludes
                                     # surface slivers + the deep locked band)
   MU_D = 0.10 ;  S_TARGET = 1.7 ;  LAMBDA_MAX = 0.9 ;  SN_FLOOR_MPA = 0.5
   STRIKE_AZ = 314.0 ; HYPO = (606971.0, 3707270.0, -4965.62)
   ```

### Acceptance Criteria
- [ ] The V3 production fault mesh is identified and written into the Phase-1
      tool header; the same mesh path is used for Phase 1 and Phase 4.
- [ ] Gate bands and depth windows above are encoded as module constants.

### Dependencies
- Depends on: nothing.  Required by: Phases 1-4.

## Phase 1: Feasibility sweep tool + B-vs-C report (DECISION GATE)

### Goal
A single reusable tool projects the field once per k, applies the graded
fault-local overpressure model, and emits a markdown report quantifying, for
recipes B and C across ALL four gates, whether every gate can reach S < 1.7
while the fault stays subcritical and sub-lithostatic — so the user can pick B
or C on real numbers.

### Files to Create
- `gate_removal_sweep.py` (toolbox) — the sweep/feasibility tool.
- `GATE_REMOVAL_SWEEP_2026-06-17.md` (toolbox) — generated report (the tool
  writes it; checked in as the decision artifact).

### Detailed Requirements
1. Reuse the deep-mesh projection path from `project_csm_on_deepmesh.py`
   (`read_gmsh_fault`) and the physics from `project_csm_stress_to_vtu.py`.
   Factor the per-k projection into:
   ```
   def project_field(k):
       # returns dict with mu_app, tau, sigma_n_eff, Sv_total, Sv_eff, Pp_hydro,
       # s_km, zkm, masks  -- all per facet, comp+ MPa.
   ```
   The orientation interpolation (`interpolate_csm_field`) and geometry are
   k-independent; compute them ONCE and only re-run `magnitudes_C1` +
   `build_tensor_from_axes` + `resolve_tractions` per k (cheap).
2. Graded per-facet feasibility, given (k, mu_s, gate_mask):
   ```
   floor = (mu_s + S_TARGET*MU_D) / (1 + S_TARGET)
   need  = gate_mask & SEIS & ACTIVE & (mu_app < floor)
   DPp_req      = sigma_n_eff - tau/floor                       # to reach floor
   cap_tension  = sigma_n_eff - SN_FLOOR_MPA
   cap_litho    = (LAMBDA_MAX - Pp_hydro/Sv_total) * Sv_total   # = LAMBDA_MAX*Sv_total - Pp_hydro
   cap          = minimum(cap_tension, cap_litho)
   unfixable    = need & (DPp_req > cap)
   DPp_applied  = clip(DPp_req, 0, cap)                         # per facet (graded)
   mu_app_new   = tau / maximum(sigma_n_eff - DPp_applied, EPS)
   ```
3. Depth-profile (production-realizable) feasibility, given (k, mu_s, gate band):
   bin the seismogenic gate-band facets by depth (e.g. 1 km bins).  Per bin z:
   ```
   DPp_need(z)        = max over band facets at z of clip(sigma_n_eff - tau/floor, 0, inf)
   DPp_cap_subcrit(z) = min over band facets at z of clip(sigma_n_eff - tau/mu_s, 0, inf)
   DPp_cap(z)         = min over band facets at z of minimum(cap_tension, cap_litho)
   feasible(z)        = DPp_need(z) <= min(DPp_cap_subcrit(z), DPp_cap(z))
   DPp_profile(z)     = DPp_need(z)  (the just-enough envelope)  if feasible(z)
   ```
   A gate is "envelope-feasible" iff feasible(z) holds for every populated bin.
   Report the per-bin DPp_profile and the resulting lambda_op(z) =
   (Pp_hydro(z) + DPp_profile(z)) / Sv_total(z).
4. Subcriticality after baking: apply DPp_applied (graded) and, separately,
   DPp_profile (envelope) to the band facets, recompute mu_app over the FULL
   active set, and assert `max(mu_app) < mu_s`.  (Off-band facets unchanged;
   band facets must not exceed mu_s.)
5. Recipes to evaluate (with the auto mu_s rule `mu_s = max_mu_app(k) + 0.02`,
   rounded up to 3 decimals):
   - Recipe B: k = 1.8.
   - Recipe C: k = 1.5.
   - Plus a frontier sweep k in {2.39, 2.2, 2.0, 1.8, 1.6, 1.5} for context.
6. Report contents (`GATE_REMOVAL_SWEEP_2026-06-17.md`), plain text:
   - The k-ceiling table (fact 2 above), regenerated from this run.
   - For B and C, a per-gate table: n_seis facets, floor, %facets reaching
     S<1.7 (graded and envelope), #unfixable, DPp median/p90/max, lambda
     median/max, envelope-feasible (Y/N), max mu_app(active), subcritical (Y/N).
   - A one-line verdict per recipe: "removes {which gates}; needs lambda up to X".
   - A recommendation (B vs C) with the trade-off stated (B keeps k closer to
     the Cajon Pass strong-crust value; C is weaker but needs least Pp).

### Interfaces
- CLI: `gate_removal_sweep.py [--mesh PATH] [--out REPORT.md]
  [--k-list 2.39,2.2,2.0,1.8,1.6,1.5] [--recipes B,C] [--lambda-max 0.9]`.
- Pure-analysis: writes only the report + optional per-recipe CSV; does NOT
  modify any production artifact.

### Edge Cases to Handle
- Pure-barrier cores with tau so small that `sigma_n_eff - tau/floor > cap`:
  count as unfixable, list their (s,z), and state in the report that they
  remain barriers (do NOT silently drop them from the %).
- Depth bins with < 10 facets: skip with a logged note (avoid noisy envelopes).
- A gate where DPp_need(z) > DPp_cap_subcrit(z) (band too heterogeneous for a
  single envelope): mark envelope-infeasible and recommend either narrowing the
  band or the per-facet graded bake (Phase 2 supports both).

### Acceptance Criteria
- [ ] `gate_removal_sweep.py` runs under pythonenv on the production mesh and
      writes `GATE_REMOVAL_SWEEP_2026-06-17.md`.
- [ ] The report shows, for B and C, per-gate %S<1.7, #unfixable, lambda<=0.9
      envelope feasibility, and subcriticality.
- [ ] The k-ceiling table reproduces fact 2 (max mu_app at k=2.39 = 0.449 +/-
      0.001) — sanity that the projection matches the established baseline.
- [ ] Reproduces the established San Gorgonio graded result (0 unfixable for
      both B and C) — regression against the diagnostic.

### DECISION GATE
After Phase 1, STOP and report B-vs-C to the user; proceed to Phase 2/3 only
with the user's chosen recipe (k, mu_s) and the per-gate lambda_op envelopes.

### Dependencies
- Depends on: Phase 0.  Required by: Phases 2-4 (supplies the chosen recipe).

## Phase 2: Fault-local overpressure baking in the projection (default OFF)

### Goal
The projection physics can subtract a fault-local graded overpressure
`DeltaPp(s,z)` from the effective tensor, controlled by a config; with the
config absent it is a no-op and all existing artifacts remain byte-identical.

### Files to Create
- `fault_local_overpressure.py` (toolbox) — the overpressure field model.

### Files to Modify
- `project_csm_stress_to_vtu.py` — add an optional overpressure step after the
  effective tensor is built, before `resolve_tractions`.
- `project_csm_on_deepmesh.py` — pass the overpressure config through so gate
  re-detection (Phase 4) sees the baked field.
- `csm_stress_to_asagi.py` — apply the same overpressure to the volume grid
  before `write_stress_asagi` (Phase 3 uses this).

### Detailed Requirements
1. `fault_local_overpressure.py` defines:
   ```
   def load_overpressure_config(path) -> dict | None
       # JSON: {"gates": [{"name","s_lo_km","s_hi_km","taper_km",
       #                    "z_profile": [[depth_m, DPp_MPa], ...]}],
       #        "lambda_max": 0.9}
       # depth profile is piecewise-linear in DPp(MPa) vs depth_m, 0 outside range.

   def delta_pp_mpa(cfg, x, y, z, strike_az, hypo, sv_total_fn) -> np.ndarray
       # For points (x,y,z) [m] returns DeltaPp [MPa] >= 0:
       #   s_km = ((x,y)-hypo[:2]) . su / 1000
       #   for each gate: w = cosine_taper(s_km, s_lo, s_hi, taper_km) in [0,1]
       #                  dpp = interp(depth=-z, z_profile)   (0 outside profile)
       #                  acc = max(acc, w*dpp)   (gates don't overlap; max is safe)
       #   clip so (Pp_hydro(z)+DeltaPp) <= lambda_max * sv_total(z)   (final guard)
       #   returns acc
   ```
   `cosine_taper` = 1 inside [s_lo,s_hi], 0.5*(1-cos) ramp over taper_km at each
   edge, 0 beyond.  No facet-normal needed (works for volume nodes AND facet
   centroids), which is why the carrier is (s,z) not per-facet.
2. In `project_csm_stress_to_vtu.py`:
   - add `--overpressure-config PATH` (default None).
   - after `sigma_facet = build_tensor_from_axes(...)` (and the optional taper),
     if cfg is not None:
     ```
     dpp = delta_pp_mpa(cfg, centroids[:,0], centroids[:,1], centroids[:,2],
                        args.strike_hint_az, HYPO_used, sv_total_fn)
     for i in (0,1,2):
         sigma_facet[:, i, i] -= dpp          # comp+ MPa, isotropic
     ```
   - record in the summary `params`: `overpressure_config` path + a stats block
     (DeltaPp field_stats, lambda_after field_stats) when active; ABSENT when
     inactive (so the constant/None summary stays byte-identical).
   - HYPO for the projection main: the projection currently has no hypocenter
     origin baked in; thread `HYPO` and `STRIKE_AZ` from the same constants used
     by the deep-mesh tool so s_km is identical across tools.
3. In `csm_stress_to_asagi.py`:
   - add `--overpressure-config PATH` (default None).
   - inside the grid loop, after `sigma = build_tensor_from_axes(...)` per
     z-slice, if cfg: compute `dpp = delta_pp_mpa(cfg, GX.ravel(), GY.ravel(),
     full(zz), ...)` (per (x,y) at this z) and subtract from sigma diag in
     comp+ MPa BEFORE the comp-negative store; equivalently add dpp*1e6 to the
     comp-negative `comps[s_xx|s_yy|s_zz]`.  Off-diagonals unchanged.
   - extend the nc `attrs` with `overpressure: <config summary>` when active.
4. In `project_csm_on_deepmesh.py`: add `--overpressure-config PATH`, thread it
   into the projection block and into the gate detection so the re-detected
   gates.json reflects the baked field.

### Interfaces
- One config schema shared by all three producers (`load_overpressure_config`).
- `delta_pp_mpa(...)` is the single source of the DeltaPp field; the projection
  and the nc builder MUST call it identically so the self-check (which
  re-projects the nc) closes.

### Edge Cases to Handle
- cfg None -> all producers byte-identical to current outputs (assert in tests).
- Surface z-slice where sv_eff <= 0 (nc builder already zeroes it): DeltaPp must
  also be 0 there (no negative effective stress).
- Final lambda guard: even if a z_profile entry over-specifies DPp, clamp so
  Pp <= lambda_max*Sv_total (never drive sigma_n_eff below SN_FLOOR).

### Acceptance Criteria
- [ ] With no `--overpressure-config`, `project_csm_stress_to_vtu.py` and
      `csm_stress_to_asagi.py` reproduce the current VTU/nc byte-for-byte
      (md5 compare against a pre-change copy).
- [ ] Unit check: with a 1-gate config, on the band facets sigma_n_eff drops by
      exactly the sampled DeltaPp and tau is unchanged to 1e-9 MPa (isotropy).
- [ ] Unit check: `delta_pp_mpa` evaluated at facet centroids vs at the nearest
      volume-grid nodes agree to within the grid-resolution tolerance (so the
      nc self-check will close).

### Dependencies
- Depends on: Phase 1 (chosen recipe + envelopes -> the config file).
- Required by: Phase 3.

## Phase 3: Regenerate production artifacts for the chosen recipe

### Goal
The V3 LSW case folder carries a gate-removed stress nc and a matching friction
yaml for the chosen recipe, with provenance updated and the self-check passing.

### Files to Create
- `overpressure_<recipe>_2026-06-17.json` (toolbox) — the chosen config from
  Phase 1 (gate bands + lambda_op depth profiles + lambda_max).

### Files to Modify
- `../../safs_seisol_v3_0_0_LSW/safs_stress_csm.nc` — regenerated via
  `csm_stress_to_asagi.py --k-ratio <chosen> --overpressure-config <json>`.
- `../../safs_seisol_v3_0_0_LSW/safs_fault.yaml` — set `mu_s` (LuaMap) to the
  chosen value; keep the deep locked band (z in (-20000,-15000) -> 1e6); update
  the DESIGN NOTE to describe k, mu_s, and the fault-local overpressure (and that
  the gate is now EXPECTED TO PROPAGATE).
- `../../safs_seisol_v3_0_0_LSW/safs_initial_stress.yaml` — update the header
  provenance (new k, "+ fault-local overpressure at the gate(s), lambda<=0.9");
  update the !ConstantMap fallback ONLY if k changed the regional tensor (note:
  fallback is unreachable in practice; leave numbers but note they are the
  hydrostatic-k=2.39 tensor unless the user wants them updated).

### Detailed Requirements
1. Generate the nc on the SAME grid as production (dx=1000, dz=250, box
   xmin..zmax defaults) with `--k-ratio <chosen>` and `--overpressure-config`.
2. The nc self-check must PASS against a FRESH overpressure-baked projection VTU
   (regenerate `csm_yhsm2013_stress_on_safs_mesh_fault_stress.vtu` with the same
   k + config first, so ground truth matches).  Do NOT check a baked nc against
   the old hydrostatic VTU.
3. `safs_fault.yaml` mu_s: use the recipe's mu_s (= max mu_app(k) + 0.02 margin,
   from the Phase-1 report).  Verify the comment's "> fault-wide max mu_app"
   line is updated with the new number.
4. Keep mu_d=0.10, d_c=2.5, cohesion=0, nucleation block UNCHANGED unless the
   Phase-1 report shows the hypocenter facet's mu_app moved (k change lowers it;
   re-state tau_s, SE, dtau_nuc checks in the comment with the new numbers — the
   nucleation patch must still nucleate: r_os > Lnuc/2).

### Edge Cases to Handle
- If k changed, the hypocenter mu_app drops (e.g. 0.296 -> lower); confirm the
  20 MPa nucleation over-stress still satisfies r_os = R*sqrt(ln(dtau/SE)) >
  Lnuc/2 with the new sigma_n / SE; if not, FLAG (do not silently retune
  nucleation — that is a separate user decision).
- nc size/grid unchanged so downstream paths and ASAGI config still resolve.

### Acceptance Criteria
- [ ] `csm_stress_to_asagi.py ... --k-ratio <chosen> --overpressure-config <json>`
      prints `SELF-CHECK: PASS` against the freshly baked VTU.
- [ ] Max nc-interpolated mu_app over the active set < chosen mu_s (printed).
- [ ] `safs_fault.yaml` mu_s matches the recipe; deep barrier + mu_d + nucleation
      intact; DESIGN NOTE updated to "expected to propagate past San Gorgonio".

### Dependencies
- Depends on: Phase 1 (recipe), Phase 2 (baking).  Required by: Phase 4.

## Phase 4: Validation + harden the self-check

### Goal
The gate removal is demonstrated end-to-end on the production mesh (S<1.7 across
the targeted gates, subcritical, sigma_n_eff>0), the figure/JSON show the gates
gone, and the nc self-check no longer silently passes on a corrupt nc.

### Files to Modify
- `project_csm_on_deepmesh.py` — run with `--overpressure-config` to regenerate
  `csm_C1k239_deepmesh_gates.json` + `..._strike_depth.png` (renamed/suffixed
  for the recipe, e.g. `csm_gateremoved_<recipe>_deepmesh_*`); keep the original
  baseline artifacts.
- `csm_stress_to_asagi.py` — apply REVIEW_nc fixes P-001 (out-of-grid coverage
  abort), P-002 (finiteness abort before write), P-005 (gate on p95 not median).

### Detailed Requirements
1. Re-detect gates on the baked field and assert: for each TARGETED gate, the
   seismogenic arrest fraction (S>=1.7 OR pure-barrier) drops below the
   detection threshold (gate no longer flagged), and per-facet `S < 1.7` holds
   for >= 95% of seismogenic band facets (report the residual % and list any
   unfixable barrier cores explicitly).
2. Subcriticality: assert `max(mu_app over active set) < mu_s` on the baked field.
3. Tension: assert `min(sigma_n_eff over all facets) > SN_FLOOR_MPA` on the baked
   field (or, for the always-tiny surface slivers, > 0 — match the existing
   `Sv_eff<=0 clamp` policy and document it).
4. Self-check hardening (REVIEW_nc):
   - P-002: before `write_stress_asagi`, abort if any STRESS_FIELDS sample is
     non-finite (count + offending grid-node count).
   - P-001: after sampling, abort if any fault facet centroid is outside the nc
     grid box (coverage), with the deepest-facet vs z-floor message.
   - P-005: gate PASS on `percentile(rel,95) < tol` (retune `--selfcheck-tol`).
5. Produce a short `GATE_REMOVAL_RESULT_<recipe>_2026-06-17.md`: before/after S
   stats per gate, the lambda_op profiles used, subcriticality margin, and the
   figure path.

### Edge Cases to Handle
- A targeted gate with residual unfixable barrier cores (deep, highly clamped):
  report it; the rupture may still cross if the cores are isolated (< process
  zone Lb wide) — note this is a SeisSol-run question, out of local scope.
- "All seismogenic gates" stretch goal: if g266 cannot be made envelope-feasible
  under lambda<=0.9, report which gates ARE removed and stop (do not raise
  lambda past 0.9 without asking).

### Acceptance Criteria
- [ ] Baked deep-mesh re-detection shows San Gorgonio removed (>=95% seismogenic
      facets S<1.7), and reports status of the other three gates.
- [ ] Subcriticality and tension asserts pass on the baked field.
- [ ] Hardened `csm_stress_to_asagi.py` aborts (non-zero) on a deliberately
      under-covered grid (`--zmin -12000`) and on an injected NaN (P-001/P-002
      tests).
- [ ] `GATE_REMOVAL_RESULT_<recipe>_2026-06-17.md` written with before/after.

### Dependencies
- Depends on: Phase 3.  Required by: nothing.

## Testing Strategy
- Phase 1: regression vs the established diagnostic (max mu_app(k) table;
  San Gorgonio 0-unfixable for B and C).  No new physics — same projection.
- Phase 2: (a) byte-identical OFF test (md5 of VTU + nc vs pre-change);
  (b) isotropy test (tau unchanged, sigma_n_eff drops by DeltaPp);
  (c) field-consistency test (delta_pp_mpa at facet vs grid nodes).
- Phase 3: nc self-check PASS against the freshly baked VTU; nc-interpolated max
  mu_app < mu_s.
- Phase 4: gate-removal assertions on the baked deep-mesh re-detection;
  P-001/P-002/P-005 self-check unit tests (abort on under-coverage / NaN /
  localized error).
- All Python runs use `/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python`.

## Risk Assessment
- Volume-grid vs facet mismatch: the overpressure carrier is a smooth (s,z)
  field, but S is measured per facet.  Risk: the coarse nc grid (dx=1000) blurs
  the envelope so a facet ends up below floor or above mu_s after sampling.
  Detection: Phase 4 re-projects the BAKED nc (not the analytic field) and
  re-checks S and subcriticality per facet.  Mitigation: tighten the depth
  profile or, if a band is too heterogeneous (DPp_need > DPp_cap_subcrit),
  fall back to the per-facet graded bake on a refined near-fault grid.
- Over-raising within a gate band: a uniform-along-strike envelope can push the
  already-high-mu_app band facets toward mu_s.  Detection: the Phase-1
  DPp_cap_subcrit(z) check and Phase-4 subcriticality assert.  Mitigation:
  narrow the band taper or lower lambda_op(z) to the subcritical cap.
- k change ripples into nucleation: lowering k lowers the hypocenter mu_app and
  sigma_n, changing the nucleation energy budget.  Detection: Phase 3 re-states
  the r_os vs Lnuc check.  Mitigation: FLAG for user decision; do not retune
  nucleation silently.
- Fallback !ConstantMap tensor in safs_initial_stress.yaml is the old k=2.39
  hydrostatic tensor; if k changes it is inconsistent with the ASAGI field.
  Practically unreachable (mesh containment), but note it; update only on
  request.
- Two meshes in play (deep 132k vs PUML 60k).  Phase 0 resolves which is the
  run mesh; using the wrong one would validate S on facets SeisSol never sees.
