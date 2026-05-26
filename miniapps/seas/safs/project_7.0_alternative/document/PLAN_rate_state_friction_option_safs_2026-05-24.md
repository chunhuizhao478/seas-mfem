# Implementation Plan: selectable LSW / rate-and-state friction in the SAFS spatial driver — 2026-05-24

## Overview

Make **linear slip-weakening (LSW, existing)** and **rate-and-state (RS, aging-law,
"like TPV102")** two user-selectable friction options in `seas_spatial_dyn_driver`,
switched by the existing TOML field `[meta].law = "slip_weakening" | "rate_state"`. The
LSW path must stay **byte-exact** when `law = slip_weakening` (the RS branch is inert in
that case). Reference machinery is this repo's own dynamic RS implementation (TPV102 aging
law); **Tandem is explicitly out of scope** here — the SAFS spatial path is fully-dynamic
ADER-DG modelled on the TPV* family, not the BP5 quasi-dynamic path.

This is a **wiring** task: the RS solver, state evolution, ADER flux, dispatch, config
schema, and per-DOF resolver already exist. Two pieces are genuinely new (an equilibrium-ψ
seed and a nucleation-callback overload on the RS iterator). **No implementation until this
plan is approved.**

## What already exists (reused verbatim — DO NOT reimplement)

| Piece | Location | Status |
|---|---|---|
| Config enum `FrictionLawKind{SlipWeakening,RateState}` | `spatial/code/spatial_friction.hpp:70` | done |
| `RateStateBlock` (f_0,V_0,eta_auto,a,b,Dc,V_init,sigma_n + spatial rules) | `spatial_friction.hpp:241-253` | done |
| `cfg.rate_state` optional + `cfg.law` | `spatial_friction.hpp:258,269` | done |
| `RateStatePerDOFParams{a,b,Dc,V_init,f_0,V_0,eta,sigma_n_eff}` | `spatial_friction.hpp:293-296` | done |
| `ResolveRateState(...)` per-DOF resolver (+validation, ParMesh & Mesh overloads) | `spatial_friction.cpp:957-1115` | **fully implemented** |
| `InitializeFaultDOFs_Spatial_RS(...)` per-DOF DOFData seeding | `dynamic/spatial_setup.hpp:302-346` | done (but seeds `psi=0` — see Gap 4) |
| DOFData RS fields `a`, `Dc`, `psi` | `dynamic/fault_face_flux.hpp:47-50` | done |
| RS solver `SolveSlipRatePsi` (Brent), `FrictionCoefficientPsi`, `InitialStatePsi` | `friction/dieterich_ruina.hpp` | done |
| Aging-law state evolution `AgingLawPsi(b,V0,f0)` | `friction/state_evolution.hpp:177` | done |
| RS dynamic substep iterator `Tpv102SubStepIterator(FaultFaceFlux&, const AgingLawPsi&)` | `dynamic/tpv102_substep_iterator.hpp:68-137` | done |
| RS ADER fault flux `EvaluateADER(...)` | `dynamic/fault_face_flux.cpp:633` | done |
| Dispatch `FaultFrictionLaw::RateAndState → EvaluateADER` (interior + shared) | `dynamic/wave_operator.inl:3916, 4887` | done |

## Constraints

- **LSW byte-exactness**: when `law = slip_weakening`, every change below is a no-op
  branch; the existing LSW code path is untouched.
- **No new friction physics**: reuse `AgingLawPsi` + `dieterich_ruina` + `EvaluateADER`
  as-is. The aging law is SCEC Eq. (2): `dθ/dt = 1 − Vθ/Dc`, ψ-space
  `dψ/dt = (b·V0/Dc)·[exp((f0−ψ)/b) − V/V0]` (`AgingLawPsi::Rate`,
  `state_evolution.hpp:190`).
- **Global b, V0, f0** (the `AgingLawPsi` ctor is scalar). Per-DOF `a`, `Dc`,
  `sigma_n_eff`, `tau_pre`, `psi` vary via DOFData. SAFS `[friction.rate_state]` per-DOF
  spatial rules for `b/V_0/f_0` are NOT honoured in this first cut (asserted away — see
  Gap 2). This matches TPV102, which varies only `a` per-DOF.
- **Nucleation reuse**: the existing `gradual_overstress` shear-stress accumulator
  (writes `tau{1,2}_nuc`) is reused for RS; the RS trial traction consumes
  `tau{1,2}_0 + tau{1,2}_nuc + trial` identically to LSW. No new nucleation kind.
- **No hardcoded constants** (project rule): all RS params come from the TOML
  `[friction.rate_state]` block / per-DOF resolver.

## Design decisions (request explicit confirmation)

- **D1 — state law = aging (`AgingLawPsi`, TPV102).** Per the user's "like tpv102".
  Slip-law / SRW (TPV104) is deferred. *Confirm.*
- **D2 — global b/V0/f0; per-DOF a/Dc.** First cut restricts `[friction.rate_state]`
  spatial rules to `a`, `Dc` (+ `sigma_n` via stress projection). If the TOML carries
  `b`/`V_0`/`f_0` spatial rules, the driver **aborts** with a clear message (rather than
  silently ignoring them). *Confirm, or request per-DOF b/V0/f0 (larger change — would
  need `AgingLawPsi` to take per-DOF vectors).*
- **D3 — equilibrium initialisation (single-phase).** Seed per-DOF
  `ψ_i = InitialStatePsi(|τ_pre_i|, V_init_i, σ_n_eff_i, η_s_i, a_i)` so the fault starts
  at steady state for `V_init` (no iterative 4-phase). *Confirm.*
- **D4 — nucleation = `gradual_overstress` (shear).** Reused for RS via a new callback
  overload on `Tpv102SubStepIterator`. *Confirm (vs. an RS-native V/state perturbation,
  deferred).*

---

## Phase 1 — config acceptance + dispatch scaffolding (driver only)

### Goal
The driver accepts `law = rate_state` and threads a single `is_lsw` boolean through the
existing setup; no behaviour change for LSW.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`

### Change 1.1 — remove the `rate_state` rejection (lines 714-719)
```diff
   const bool is_lsw =
      (cfg.law == spatial::FrictionLawKind::SlipWeakening);
-  MFEM_VERIFY(is_lsw,
-              "spatial_dyn_driver: only [meta].law = \"slip_weakening\" "
-              "is supported in this commit; rate_state path is a "
-              "deferred follow-up (plan §Phase 4 Edge Cases).");
+  // D2 first-cut guard: rate_state honours per-DOF a/Dc + global b/V0/f0 only.
+  if (!is_lsw)
+  {
+     MFEM_VERIFY(cfg.rate_state.has_value(),
+                 "spatial_dyn_driver: [meta].law=rate_state but the "
+                 "[friction.rate_state] block is absent in TOML.");
+     // Reject per-DOF b/V_0/f_0 spatial rules (not yet supported; AgingLawPsi
+     // is global in b/V0/f0).  a/Dc/sigma_n spatial rules ARE allowed.
+     for (const auto &r : cfg.rate_state->spatial)
+     {
+        MFEM_VERIFY(r.field != "b" && r.field != "V_0" && r.field != "f_0",
+                    "spatial_dyn_driver: per-DOF spatial rule on '" << r.field
+                    << "' is not supported (b/V_0/f_0 are global in this cut); "
+                    "use the [friction.rate_state] default instead.");
+     }
+  }
```
> NOTE: the exact field name (`r.field`) must be confirmed against `SpatialRule`
> (`spatial_friction.hpp:189-201`) during implementation — if the struct names the target
> differently (e.g. `r.param`), use that. This is the one identifier to verify before
> coding.

### Acceptance Criteria
- [ ] `law = slip_weakening` run is byte-identical to current `main` (the new branch is
      not entered).
- [ ] `law = rate_state` with only `a`/`Dc` spatial rules passes the guard; with a
      `b`/`V_0`/`f_0` spatial rule it aborts with the message above.

---

## Phase 2 — per-DOF parameter resolve + DOFData init + equilibrium ψ (driver)

### Goal
Under `law = rate_state`, resolve per-DOF RS params, seed DOFData via the RS initializer,
and overwrite `psi` with the steady-state equilibrium value.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` (§11 resolve, §14 init)

### Change 2.1 — branch the resolver (current lines 1146-1152)
```diff
-  MFEM_VERIFY(cfg.slip_weakening.has_value(),
-              "spatial_dyn_driver: [meta].law=slip_weakening but the "
-              "[friction.slip_weakening] block is absent in TOML.");
   spatial::SpatialFrictionResolver resolver;
-  const spatial::SlipWeakeningPerDOFParams lsw =
-     resolver.ResolveSlipWeakening(*cfg.slip_weakening,
-                                   dof_coords_3d, dof_to_attr);
+  spatial::SlipWeakeningPerDOFParams lsw;     // populated iff is_lsw
+  spatial::RateStatePerDOFParams    rs;       // populated iff !is_lsw
+  if (is_lsw)
+  {
+     MFEM_VERIFY(cfg.slip_weakening.has_value(),
+                 "spatial_dyn_driver: [meta].law=slip_weakening but the "
+                 "[friction.slip_weakening] block is absent in TOML.");
+     lsw = resolver.ResolveSlipWeakening(*cfg.slip_weakening,
+                                         dof_coords_3d, dof_to_attr);
+  }
+  else
+  {
+     // ResolveRateState signature (spatial_friction.cpp:1090) — confirm the
+     // exact arg list during implementation; it needs dof_coords_3d,
+     // dof_to_elem, dof_to_attr, and the projected sigma_n for eta=auto.
+     rs = resolver.ResolveRateState(*cfg.rate_state, dof_coords_3d,
+                                    dof_to_elem, dof_to_attr,
+                                    /*sigma_n_total_per_dof=*/geom.sigma_n_per_dof(),
+                                    material);
+  }
```
> The precise `ResolveRateState` parameter list must be read from
> `spatial_friction.cpp:1090-1100` at implementation time (the MFEM_VERIFYs at 976-982
> reference `dof_coords_3d`, `dof_to_elem`, `dof_to_attr`, `sigma_n_total_per_dof`).

### Change 2.2 — branch DOFData init (current lines 1217-1226)
```diff
   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
-     spatial::InitializeFaultDOFs_Spatial<ParMesh>(
-        dof_data, num_fault_total, dof_to_elem, material, pmesh,
-        lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
-        dummy_T_forced, dummy_t0_decay,
-        dof_ips);
+     if (is_lsw)
+     {
+        spatial::InitializeFaultDOFs_Spatial<ParMesh>(
+           dof_data, num_fault_total, dof_to_elem, material, pmesh,
+           lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
+           dummy_T_forced, dummy_t0_decay,
+           dof_ips);
+     }
+     else
+     {
+        spatial::InitializeFaultDOFs_Spatial_RS<ParMesh>(
+           dof_data, num_fault_total, dof_to_elem, material, pmesh,
+           rs, geom.GetTauPre(), geom.sigma_n_per_dof());
+        // Gap 4 — equilibrium psi: InitializeFaultDOFs_Spatial_RS seeds psi=0
+        // (spatial_setup.hpp:341).  Overwrite with the steady-state value so
+        // the fault starts at V_init in equilibrium.
+        SeedEquilibriumPsi_RS(dof_data, rs);   // new helper, Change 2.3
+     }
   }
```

### Change 2.3 — NEW equilibrium-ψ helper (file-local, in the anonymous namespace)
```cpp
// Seed per-DOF state variable psi from steady state at V_init so the RS fault
// starts in equilibrium: tau_pre = sigma_n_eff * f(V_init, psi) + eta*V_init.
// Uses the existing InitialStatePsi (friction/dieterich_ruina.hpp).  No
// iterative 4-phase init (TPV102-style single equilibrium seed, D3).
static void SeedEquilibriumPsi_RS(std::vector<DOFData> &dof_data,
                                  const spatial::RateStatePerDOFParams &rs)
{
   for (size_t i = 0; i < dof_data.size(); ++i)
   {
      DOFData &d = dof_data[i];
      const real_t tau0 = std::sqrt(d.tau1_0 * d.tau1_0 + d.tau2_0 * d.tau2_0);
      d.psi = mfem::seas::DieterichRuina::InitialStatePsi(
                 tau0, rs.V_init(i), d.sigma_n0, d.eta_s, d.a);
   }
}
```
> Confirm the exact `InitialStatePsi` signature/namespace (`dieterich_ruina.hpp`) and the
> DOFData impedance field name (`eta_s`) at implementation time. `d.sigma_n0` is the
> effective normal stress already seeded by `InitializeFaultDOFs_Spatial_RS` from
> `sigma_n_eff`.

### Edge Cases
- `num_fault_total == 0` (rank with no fault DOFs): both branches skip — unchanged.
- `V_init`, `sigma_n_eff`, `a` already validated `> 0` by `ResolveRateState`
  (`spatial_friction.cpp:1068-1076`), so `InitialStatePsi` inputs are well-posed.

### Acceptance Criteria
- [ ] RS run: every `dof_data[i].a/.Dc/.psi` is finite and `psi` ≠ 0.
- [ ] A unit test (Phase 5) asserts `tau_pre ≈ sigma_n_eff·f(V_init,psi)+eta·V_init` per DOF
      after seeding (equilibrium residual < 1e-6 relative).

---

## Phase 3 — friction-law dispatch + iterator selection + time-loop helper

### Goal
Route the RS run through `FaultFrictionLaw::RateAndState`, the `Tpv102SubStepIterator`, and
an RS time-loop helper.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`
- `dynamic/tpv102_substep_iterator.hpp` / `.cpp` (new callback overload — Phase 4)

### Change 3.1 — branch the friction law (current line 971)
```diff
-  wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
+  wave.SetFaultFrictionLaw(is_lsw ? FaultFrictionLaw::LSW
+                                  : FaultFrictionLaw::RateAndState);
```

### Change 3.2 — RS iterator + state evolution (parallel to current lines 1736-1742)
```cpp
// Built only when !is_lsw.  AgingLawPsi is GLOBAL in b/V0/f0 (D2): take them
// from the [friction.rate_state] defaults (resolver guarantees the spatial
// rules don't touch them).
std::unique_ptr<AgingLawPsi> rs_state_evo;
std::unique_ptr<Tpv102SubStepIterator> rs_iterator;
if (!is_lsw)
{
   rs_state_evo = std::make_unique<AgingLawPsi>(
      cfg.rate_state->b_default,
      cfg.rate_state->V_0_default,
      cfg.rate_state->f_0_default);
   rs_iterator = std::make_unique<Tpv102SubStepIterator>(fault_flux, *rs_state_evo);
   const int O = std::max(1, cfg.numerics.ader_order);
   std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
   std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
   rs_iterator->SetSubSteps(deltaT, weights);
   rs_iterator->SetDiagNumLocalFaultQPs(wave.GetNumLocalFaultQPs());
}
```
> The existing `Tpv205SubStepIterator substep_iterator(...)` block (1736-1746) stays, built
> only when `is_lsw` (wrap it in `if (is_lsw) { … }`).

### Change 3.3 — RS time-loop helper (NEW, parallel to `AdvanceADERWithSubStep_Spatial`)
`AdvanceADERWithSubStep_Spatial` (line 365) takes `Tpv205SubStepIterator&` concretely and
calls the **callback overload** of `AdvanceWithSubStepStates` (line 439). The
`Tpv102SubStepIterator` has a different signature (no callback; ends in
`FrictionSolver::Method`). Add a sibling helper:

```cpp
// RS twin of AdvanceADERWithSubStep_Spatial: identical predictor/EvaluateBulk/
// imposed-state plumbing (lines 384-460), but drives Tpv102SubStepIterator and
// the gradual_overstress callback added in Phase 4.
void AdvanceADERWithSubStep_Spatial_RS(
   WaveOperator<ParMesh> &wave,
   Tpv102SubStepIterator &iterator,
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const Vector &Q, real_t dt_step, int ader_order, real_t t_step_start,
   Vector &Q_new,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   /* … lines 377-435 verbatim (SetSubSteps scale, tau_nodes,
      ComputeADERSubStepStates, EvaluateBulkAtFaultQPsCanonical, I_imp alloc) … */
   if (n_total_fault_qps > 0)
   {
      iterator.AdvanceWithSubStepStates(dof_data, fault_coords,
                                        Q_pointwise_plus, Q_pointwise_minus,
                                        dt_step, t_step_start,
                                        I_imp_plus_flat.data(),
                                        I_imp_minus_flat.data(),
                                        nuc_callback);   // Phase-4 overload
   }
   /* … lines 448-460 verbatim (ImposedGuard, SetSubStepFaultImposedStates,
      AdvanceADER) … */
}
```
> Rationale for a parallel helper rather than templating: the two iterators' public
> `AdvanceWithSubStepStates` signatures differ (callback vs `FrictionSolver::Method`), so a
> template would not unify cleanly. ~40 lines, mostly identical; the shared body (384-435,
> 448-460) is a candidate for a small `static` helper to avoid duplication, but
> straight duplication is acceptable for review clarity. *Decide during implementation.*

### Change 3.4 — branch the time-loop call (current line 1893)
```diff
-     AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
-                                    fault_coords, Q, dt_step, ader_order,
-                                    t_step_start, Q_new, nuc_cb);
+     if (is_lsw)
+     {
+        AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
+                                       fault_coords, Q, dt_step, ader_order,
+                                       t_step_start, Q_new, nuc_cb);
+     }
+     else
+     {
+        AdvanceADERWithSubStep_Spatial_RS(wave, *rs_iterator, dof_data,
+                                          fault_coords, Q, dt_step, ader_order,
+                                          t_step_start, Q_new, nuc_cb);
+     }
```
(`nuc_cb`, lines 1753-1761, is unchanged — it calls `ApplyGradualOverstressIncrement`,
which writes `tau{1,2}_nuc` and is friction-law-agnostic.)

### Acceptance Criteria
- [ ] LSW run unchanged (the `is_lsw` branches select the existing code).
- [ ] RS run dispatches `EvaluateADER` (verify via the `FaultFrictionLaw::RateAndState`
      arm at `wave_operator.inl:3916`).

---

## Phase 4 — gradual_overstress nucleation callback on the RS iterator (NEW iterator code)

### Goal
`Tpv102SubStepIterator::AdvanceWithSubStepStates` must accept a per-sub-step
`nuc_callback` (as `Tpv205SubStepIterator` already does) so the SAFS `gradual_overstress`
shear increment fires before each sub-step's friction+ψ solve.

### Files to Modify
- `dynamic/tpv102_substep_iterator.hpp` — add the overload declaration
- `dynamic/tpv102_substep_iterator.cpp` — implement it

### Change 4.1 — new overload declaration (`tpv102_substep_iterator.hpp`, after line 127)
```cpp
/// Callback overload (SAFS spatial driver): `nuc_callback(t_sub_end, dt_sub)`
/// fires ONCE per ADER sub-step BEFORE the per-QP friction+psi pipeline, so a
/// gradual_overstress increment into DOFData::tau{1,2}_nuc is visible to the
/// trial traction.  Mirrors Tpv205SubStepIterator's callback overload.
void AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro, real_t t_macro_start,
   real_t *I_imp_plus_flat, real_t *I_imp_minus_flat,
   const std::function<void(real_t, real_t)> &nuc_callback,
   FrictionSolver::Method method = FrictionSolver::Method::NewtonRaphsonStable);
```

### Change 4.2 — implementation (`tpv102_substep_iterator.cpp`)
The existing (non-callback) `AdvanceWithSubStepStates` already runs the per-sub-step loop
(friction solve → analytic ψ update via `UpdateStateAnalytic` aging law → imposed-state
accumulate). The new overload is that loop with **one added line per sub-step**: call
`nuc_callback(t_sub_end, dt_sub)` at the sub-step endpoint, immediately before the per-QP
friction pipeline — exactly where `Tpv205SubStepIterator` fires its callback. Preferred
implementation: factor the existing loop body into a private method and have both overloads
call it (the callback overload passing a non-empty callback; the original passing
`[](real_t,real_t){}`), to avoid divergence.

> The exact insertion point + the per-sub-step `t_sub_end`/`dt_sub` derivation must mirror
> `tpv205_substep_iterator.cpp`'s callback overload (read it at implementation time so the
> nucleation-endpoint convention matches bit-for-bit).

### Edge Cases
- `nuc_callback` empty (`[](real_t,real_t){}`) ⇒ identical to the non-callback overload
  (must be byte-exact — guard with a test).
- Nucleation disabled (`cfg.nucleation.enabled == false`): the driver's `nuc_cb`
  early-returns; the increment is a no-op.

### Acceptance Criteria
- [ ] Calling the new overload with an empty callback is byte-identical to the existing
      `AdvanceWithSubStepStates` (no nucleation regression for native TPV102).
- [ ] With the SAFS `gradual_overstress` callback, `tau{1,2}_nuc` accumulates the smoothStep
      increment over `[0,T_nuc]` and the RS rupture nucleates.

---

## Phase 5 — tests

### Files to Create
- `tests/unit/test_spatial_rs_equilibrium_seed.cpp` — Phase-2 ψ seed: build a synthetic
  `RateStatePerDOFParams` + DOFData, run `SeedEquilibriumPsi_RS`, assert the equilibrium
  residual `|τ_pre − (σ_n·f(V_init,ψ)+η·V_init)| / σ_n < 1e-8` per DOF.
- `tests/unit/test_tpv102_iterator_nuc_callback_identity.cpp` — Phase-4: the new callback
  overload with an empty callback reproduces the existing overload bit-for-bit on a 2-QP
  fixture (guards the refactor).

### Files to Modify
- `tests/unit/test_spatial_friction_config.cpp` — add a `law = "rate_state"` parse case +
  the per-DOF `b/V_0/f_0` spatial-rule rejection (Change 1.1 guard). (Verify the existing
  file already covers RS parse; if so, extend only the rejection case.)
- `Makefile` — wire the two new unit tests + add to `test` aggregate and
  `test-v92-regression-gates` (per the R-002 lesson: a guard that never runs is not a
  guard).

### Integration (Frontera, not local)
- Small SAFS RS smoke: `law = rate_state`, short `tfinal`, on the Dc2 mesh. Expectations:
  pre-nucleation the fault sits at `V ≈ V_init` (locked steady state), nucleation triggers
  an RS rupture, `V_substep_max == V_max` (R-004 is law-agnostic — `EvaluateBulkAtFaultQPsCanonical`
  is shared), `max_slip` physical. A dedicated `…_rate_state_…toml` config + sbatch
  (clone the Dc2 LSW config, swap `[meta].law` + add `[friction.rate_state]`).

### Regression
- `make test-v92-regression-gates` green (LSW byte-exactness + the new RS guards).
- A `law = slip_weakening` run diffs bit-identical against pre-change (the `is_lsw`
  branches are inert).

---

## Risk Assessment

- **R-A (shared-fault RS + R-001 dual-solve).** On shared faces the macro routine
  *discards* the iterator's `I_imp` and **re-solves** `EvaluateADER` (RS) on `I/dt`
  (`wave_operator.inl:4887`, the same R-001 structure we hit for LSW). For RS this re-solve
  reads `data.psi` but the **ψ time-integration** lives in the iterator's per-sub-step
  `UpdateStateAnalytic`. The shared re-solve must consume a consistent ψ, or shared-face RS
  QPs will weaken on a ψ history the applied flux never evolved — the RS analog of R-001.
  **Mitigation:** add a shared-fault RS check to the smoke run (`R-101`-style cross-rank ψ
  consistency); treat R-001 as still-open for RS even though it's correctness-moot for LSW
  post-R-004. *This is the top risk and may surface a follow-up.*
- **R-B (R-004).** Already fixed and law-agnostic (`EvaluateBulkAtFaultQPsCanonical`), so RS
  inherits the correct neighbour read — no action.
- **R-C (equilibrium seed).** A wrong `InitialStatePsi` makes the fault creep or lock at
  t=0. Guarded by the Phase-5 equilibrium-residual test.
- **R-D (global b/V0/f0).** If a SAFS RS config needs spatial `b/V0/f0`, the Phase-1 guard
  aborts loudly (no silent wrong physics); lifting it is a scoped follow-up
  (per-DOF `AgingLawPsi`).
- **R-E (nucleation amplitude for RS).** `gradual_overstress` was tuned as a shear kick for
  LSW; for RS the same Δτ may over/under-drive nucleation. The smoke run calibrates; not a
  code risk.

## Summary of every code change (for approval)

| # | File | Change | New/Edit |
|---|---|---|---|
| 1.1 | `spatial_dyn_driver.cpp:714-719` | replace rate_state rejection with the D2 per-DOF-rule guard | edit |
| 2.1 | `spatial_dyn_driver.cpp:1146-1152` | branch resolve: `ResolveSlipWeakening` vs `ResolveRateState` | edit |
| 2.2 | `spatial_dyn_driver.cpp:1217-1226` | branch init: `InitializeFaultDOFs_Spatial` vs `_RS` + seed call | edit |
| 2.3 | `spatial_dyn_driver.cpp` (anon ns) | `SeedEquilibriumPsi_RS` helper | **new** |
| 3.1 | `spatial_dyn_driver.cpp:971` | `SetFaultFrictionLaw(is_lsw ? LSW : RateAndState)` | edit |
| 3.2 | `spatial_dyn_driver.cpp:1736-1746` | build `AgingLawPsi` + `Tpv102SubStepIterator` when `!is_lsw`; wrap LSW iterator in `if(is_lsw)` | edit |
| 3.3 | `spatial_dyn_driver.cpp` (anon ns) | `AdvanceADERWithSubStep_Spatial_RS` helper | **new** |
| 3.4 | `spatial_dyn_driver.cpp:1893` | branch the time-loop call on `is_lsw` | edit |
| 4.1 | `tpv102_substep_iterator.hpp:~127` | `nuc_callback` overload declaration | **new** |
| 4.2 | `tpv102_substep_iterator.cpp` | implement the callback overload (factor loop body) | **new** |
| 5.* | 2 new unit tests + `test_spatial_friction_config.cpp` + `Makefile` | tests + wiring | new/edit |
| 5.cfg | new `…_rate_state_…toml` + sbatch | RS smoke config (clone Dc2 LSW, swap law) | **new** |

**No code will be written until this plan is approved.** Open confirmations: D1 (aging
law), D2 (global b/V0/f0 + abort on spatial rules), D3 (single-phase equilibrium seed), D4
(reuse gradual_overstress) — and acknowledgement of risk R-A (shared-fault RS / R-001).
