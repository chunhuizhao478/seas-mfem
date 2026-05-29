# Code Review: Phase 14 — Explicit Runge–Kutta (RK4/RK45) time integrator (2026-05-29)

> Supersedes the prior [DIAG-SIGN] speckle-instrumentation review (recoverable via git history).
> Scope is the Phase 14 RK time-stepper implementation against the plan.

## Review Scope
- Plan: `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md` (§ Phase 14, lines 3031–3372)
- Files reviewed:
  - `miniapps/seas/dynamic/rk_time_stepper.hpp` (tableau + `AdvanceRKCoupled_Spatial`)
  - `miniapps/seas/dynamic/rk_time_stepper.cpp` (tableau factories, `ValidateTableau`, `PsiRate`)
  - `miniapps/seas/tests/unit/test_rk_time_stepper.cpp`
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (CLI/guards/banner ~530–650, ~755–830; dt ~1758–1810; tableau select + time-loop branch ~2470–2660)
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (`TimeIntegratorKind`, `RkCflFactor`, `CflSafetyFactor`, parser ~1093–1101)
  - `miniapps/seas/dynamic/wave_operator.{hpp,inl}` (`ComputeMaxDt` RK-aware branch ~5530–5599, `SetCflRkAware`, Mult fault paths ~2607–2620 / ~3245–3262)
  - `miniapps/seas/dynamic/nucleation_method.hpp`, `dynamic/spatial_nucleation.{hpp,cpp}` (`ApplyGradualOverstress{Absolute,Increment}` + compact-circular variants)
  - `miniapps/seas/dynamic/fault_face_flux.{hpp,cpp}` (`Evaluate`/`WriteBackState`, ψ-stateless tripwire)
  - `miniapps/seas/Makefile` (RK object + test target)
- Domain context: `CLAUDE.md`, project memory (Phase-14 RK plan, mixed-flux/RK note, CFL-calibration note), and the cited reference template `git:8461c67 drivers/tpv102_driver.cpp` (the proven coupled-RK4 loop).
- Build/run done in this review: `seas_test_rk_time_stepper` builds and **20/20 tests pass**; `drivers/spatial_dyn_driver.o` compiles cleanly (against the main-repo MFEM at the same commit `30bf3e3`).

## Verified correct (not findings — so the fix agent does not re-litigate)
- RK4 and DP45 Butcher coefficients (validator + tests pass); DP45 FSAL (`a[6]==b`, `b[6]=0`); `ValidateTableau` identities (Σb=1, Σbhat=1, c=Σa, strict-lower-triangular).
- `PsiRate` aging/SRW arms equal the underlying kernels bit-for-bit (T4). The aging arm correctly uses the **ψ-space** rate `AgingLawPsi::Rate` (not the θ-space `1−Vθ/Dc` the plan text loosely wrote), consistent with `DOFData::psi` being the ψ-space state the rest of the RS infra integrates and with the ADER `UpdateStateAnalytic`. `SlipLawSRWPsi` ctor arg order matches (`muW ← f_w_default`); per-QP `V_w(m)/a(m)` used; production-mode set.
- The (Q, ψ, slip) **stage coupling order is byte-faithful to the 8461c67 reference** for RK4: stage-ψ inputs (`ψ_n+½dt·psi_k0`, `ψ_n+½dt·psi_k1`, `ψ_n+dt·psi_k2`) and the final ψ/slip combine (`ψ_n+dt/6·(psi_k0+2psi_k1+2psi_k2+psi_k3)`, `slip += dt·Σ b_i V_k`) match term-for-term.
- Absolute-nucleation design correctly avoids the double-apply an increment form would incur on repeated-abscissa stages (RK4 stages 2&3 at t+dt/2; DP45 stages 6&7 at t+dt) — `ApplyAbsolute` SETs, is idempotent (T5).
- Coupling is wired: driver binds the same `dof_data` the Mult fault solve writes (`SetFaultDOFData(&dof_data)` at `spatial_dyn_driver.cpp:1595`); the Mult fault solve reads `data.psi` + `tau*_nuc` (via `ComputeStageState`) and writes `slip_rate/V1/V2` (`WriteBackState`) while remaining ψ-stateless (tripwire `fault_face_flux.cpp:432`); the Mult fault solve uses the default Brent method (matches production / CLAUDE.md).
- Guards present and correct: rk ⇒ `!is_lsw`, rk ⇒ `interior_flux==Scalar`, rk ⇒ `rate_state.has_value()`; `--time-integrator` parse (CLI + TOML) with abort on unknown; Makefile object + link + test target complete.
- ADER byte-exactness preserved: `cfl_rk_aware_` defaults false ⇒ `ComputeMaxDt` mixed-flux factors stay 1.0/0.9/0.4; ADER `dt` uses unchanged `CflSafetyFactor`; the RK path is additive/CLI-gated.

## Findings

### [R-001] MODERATE [rk_time_stepper.hpp:204-230 / AdvanceRKCoupled_Spatial] — RK4 reports fault observables at the wrong (non-endpoint) state; header comment provably false

**Category:** BUG (also a DEVIATION from the cited reference)

**Description:**
After the macro-step, `AdvanceRKCoupled_Spatial` leaves `dof_data[m].{slip_rate, V1, V2, tau1_corr, tau2_corr, sigma_n_corr}` at the values written by the **last stage's** `wave.Mult` (the combine loop overwrites only `psi`, `slip1`, `slip2`). The header comment (lines 227–230) justifies this as:

> `// dof_data[m].{slip_rate,V1,V2,tau*_corr,sigma_n_corr} retain the final`
> `// stage's values (c_{s-1}=1 for both RK4 and FSAL DP45 ⇒ endpoint) ...`

That claim is **false for RK4**. `c_{s-1}=1` is the stage *time*, not the stage *state*. For RK4 the last stage's input is `Q^(3) = Q + dt·k2`, which is **not** the combined endpoint `Q_new = Q + dt/6·(k0+2k1+2k2+k3)`, and the last stage's ψ is `ψ_n + dt·psi_k2`, not `psi_new`. Verified numerically on `y'=−0.7y`, dt=0.1: RK4 last-stage state `0.932364250` vs endpoint `0.932393834` (diff `2.96e-5`); DP45 last-stage `==` endpoint to `0.0e0` (FSAL). So on the `--time-integrator rk4` path the reported instantaneous fault observables correspond to an intermediate predictor state, not `Q(t+dt)`.

This is exactly the half-step phase-lag the project's own reference fixed: `git:8461c67 drivers/tpv102_driver.cpp` (round-9 `R-V92-K01`, comment: *"the values correspond to t_n + dt/2 — a half-step phase lag. Semantically wrong for instantaneous observables regardless of magnitude"*) adds an explicit **endpoint re-evaluation** — one extra `wave.Mult(Q, k_endpoint)` at `Q(t+dt)` with `psi_new` and nucleation at `t+dt` — so `dof_data` observables are self-consistent with the output time. The new stepper omits that re-eval, reintroducing the bug for RK4. **DP45 (the production integrator) is correct here for free** via FSAL. The state evolution `(Q, ψ, slip)` is correct for both schemes; only the *reported* RK4 instantaneous observables are wrong.

**Trigger:**
`--time-integrator rk4` (a user-selectable run option, not only the bulk validation harness) on any rate-and-state fault run.

**Actual behavior:**
`dof_data[m].slip_rate` (→ driver `V_max` at `spatial_dyn_driver.cpp:2655`, which gates ParaView write cadence and the blow-up monitor) and `V1/V2/tau*_corr/sigma_n_corr` (→ ParaView fault fields `:2526-2553`) are RK4 stage-4 predictor values, phase-lagged from `Q(t+dt)` by O(dt) (worse near a fast rupture front).

**Expected behavior:**
Reported fault observables should be the endpoint values at `(Q_new, psi_new, nuc(t+dt))`, matching the reference's `R-V92-K01` endpoint re-eval and matching what DP45 already produces.

**Suggested fix:** after the final combine, re-evaluate the fault observables at the endpoint when the last stage is not already the endpoint (non-FSAL). `psi` is already `psi_new` from the combine, so no extra ψ write is needed:
```diff
       dof_data[m].psi    = psi_new;
       dof_data[m].slip1 += slip1_inc;
       dof_data[m].slip2 += slip2_inc;
-      // dof_data[m].{slip_rate,V1,V2,tau*_corr,sigma_n_corr} retain the final
-      // stage's values (c_{s-1}=1 for both RK4 and FSAL DP45 ⇒ endpoint), so
-      // the driver's V_max diagnostic + ParaView read end-of-step observables.
    }
+
+   // R-001: only an FSAL tableau (a[s-1]==b) evaluates its last stage AT the
+   // endpoint; for non-FSAL (RK4) the last stage is the predictor Q+dt·k_{s-2},
+   // so dof_data observables would be phase-lagged from Q(t+dt) (the 8461c67
+   // R-V92-K01 bug). Re-evaluate at (Q_new, psi_new, nuc(t+dt)) so the driver's
+   // V_max + ParaView read self-consistent end-of-step values. (DP45 is FSAL ⇒
+   // its last stage already ran at the endpoint, so this is skipped.)
+   bool last_stage_is_endpoint = (s >= 1);
+   for (int j = 0; j < s; ++j)
+   { if (tab.a[s-1][j] != tab.b[j]) { last_stage_is_endpoint = false; break; } }
+   if (n > 0 && nuc == nuc /*always*/ && !last_stage_is_endpoint)
+   {
+      // psi already == psi_new from the combine above.
+      if (nuc) { nuc->ApplyAbsolute(dof_data, t_step_start + dt_step); }
+      Vector k_scratch(height);
+      wave.Mult(Q_new, k_scratch);   // overwrites slip_rate/V1/V2/tau*_corr/sigma_n_corr
+      for (int m = 0; m < n; ++m)
+      {
+         if (dof_data[m].slip_rate > dof_data[m].slip_rate_substep_max)
+         { dof_data[m].slip_rate_substep_max = dof_data[m].slip_rate; }
+      }
+   }
```
(Alternative, if RK4 is to stay strictly a bulk-validation stepping-stone: reword the comment to say the rk4-path fault observables are non-endpoint and document it; but the driver currently exposes rk4 as a full run option, so the re-eval is the safer fix.)

**Test case:** (runnable standalone, reuses the test file's tableau structures; pins the root cause without a fault fixture)
```cpp
// Add to test_rk_time_stepper.cpp.  The last RK stage's INPUT state equals the
// combined endpoint ONLY for FSAL (DP45), not for RK4 — so the header claim
// "c_{s-1}=1 ⇒ endpoint" is false for RK4. After the fix, the stepper must do an
// endpoint re-eval for any tableau where this returns false.
static bool LastStageIsEndpoint(const RKTableau& tab)
{
   const real_t lambda = -0.7, dt = 0.1, y0 = 1.0;
   const int s = tab.stages;
   std::vector<real_t> k(s, 0.0);
   real_t last_stage = y0;
   for (int i = 0; i < s; ++i)
   {
      real_t ys = y0;
      for (int j = 0; j < i; ++j) { ys += dt * tab.a[i][j] * k[j]; }
      k[i] = lambda * ys;                 // f(y_stage)
      if (i == s - 1) { last_stage = ys; }
   }
   real_t y_end = y0;
   for (int i = 0; i < s; ++i) { y_end += dt * tab.b[i] * k[i]; }
   return std::abs(last_stage - y_end) <= 1e-12 * std::abs(y_end);
}
void test_R001_rk4_last_stage_not_endpoint()
{
   TEST_TRUE(!LastStageIsEndpoint(MakeRK4Tableau()),
             "RK4 last stage is NOT the endpoint (driver needs endpoint re-eval)");
   TEST_TRUE(LastStageIsEndpoint(MakeDormandPrinceRK45Tableau()),
             "DP45 last stage IS the endpoint (FSAL)");
}
```

---

### [R-002] LOW [rk_time_stepper.cpp:145-184 / PsiRate] — friction-law objects reconstructed per fault DOF per RK stage

**Category:** QUALITY (performance; deviation from the reference idiom)

**Description:**
`PsiRate` constructs a fresh `AgingLawPsi` (aging arm) or a fresh `SlipLawSRWPsi` + `SetProductionMode()` (SRW arm) on **every call** — once per fault DOF per RK stage (`s × n` constructions per macro-step). The cited reference (`git:8461c67 tpv102_driver.cpp:825`) builds the law **once** outside the time loop. The law objects are stateless value-holders over global scalars, so reconstructing them is pure overhead on the hottest loop (production SAFS fault: O(10⁴–10⁵) DOFs × thousands of steps × 7 stages). Correctness is unaffected (T4 pins it); flagged only because it diverges from the proven idiom and lands on the per-stage hot path the plan's own risk table already calls out.

**Trigger:** any rk4/rk45 fault run.

**Actual behavior:** law object allocated inside the innermost `PsiRate` call.

**Expected behavior:** build the law once and reuse across DOFs/stages/steps.

**Suggested fix:** build the law objects once in `AdvanceRKCoupled_Spatial` (dispatched on `rs_cfg.state_evolution`, before the stage loop) and thread them into `PsiRate`, leaving `PsiRate` as the per-DOF `Rate`/`Rate_SRW` arithmetic only. Non-blocking.

**Test case:** N/A (performance; T4 already locks the numerics).

---

### [R-003] LOW [tests/unit/test_rk_time_stepper.cpp] — the coupled (Q, ψ, slip) fault path and the np>1 shared-fault path are never exercised

**Category:** QUALITY (test-coverage gap vs the plan's acceptance criteria)

**Description:**
Every test exercises either the **bulk-only** path (T3, no fault DOFs) or helpers in isolation (T4 `PsiRate`, T5 `ApplyAbsolute`). Nothing exercises `AdvanceRKCoupled_Spatial` with an **actual fault** — the heart of Phase 14.2/14.3/14.5 (stage-ψ-write → `Mult`/`Evaluate` → `PsiRate` → b-weighted ψ/slip combine; per-stage `ApplyAbsolute` reaching the fault solve). So the most error-prone part of the feature is unverified locally. Specifically:
- Phase 14.2 AC ("rk* + mixed=none reproduces the ADER-O2 upwind baseline") has no serial smoke proxy.
- Phase 14.5 AC ("np>1: a fault crossing a rank seam is conservation-consistent on the RK path, verified explicitly, **not inherited from ADER**") has **no test at all**, and the plan warns the ADER R-1600/R-1601 shared-fault workaround "does NOT transfer." (Code reading is reassuring — the Mult shared-fault path feeds canonical-frame Q "bit-identical across ranks", `wave_operator.inl:3245-3262`, so equal V ⇒ equal `PsiRate` ⇒ rank-consistent ψ — but this is asserted, not tested.)

**Trigger:** N/A (absence of coverage).

**Actual behavior:** coupled-fault and parallel-fault RK behavior untested.

**Expected behavior:** a serial single-fault one-step smoke + an np2 shared-fault ψ-consistency test.

**Suggested fix (new tests):**
```cpp
// (serial) one coupled RK4 step on a small cube + planar fault fixture:
//   - seed equilibrium psi; run AdvanceRKCoupled_Spatial(RK4) for one dt;
//   - assert dof_data[m].psi moved by ≈ dt·PsiRate(stage avg) (sign + O(dt));
//   - assert slip{1,2} == dt·Σ b_i V{1,2}_k (RK-consistent integral);
//   - after the R-001 fix: assert dof_data[m].slip_rate == slip_rate from a
//     fresh wave.Mult(Q_new) (endpoint self-consistency).
// (parallel np=2) fault crossing the rank seam:
//   - one RK45 step per rank; Allgather shared-DOF psi/slip_rate;
//     assert bit-identical across the seam (Phase 14.5 AC).
```

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — RK4 non-endpoint fault observables + false comment)
- Low issues: 2 (R-002 per-call law construction; R-003 missing coupled/parallel-fault tests)
- Plan compliance: **FULL** for Phases 14.0–14.5 as implemented (RK4 + DP45 tableaux, selector/CLI/parser/guards, `PsiRate`, absolute nucleation, RK-aware CFL, mixed-flux wiring, Makefile + test). Phase 14.6 (adaptive) correctly deferred (`bhat` stored, unused). The one substantive deviation is R-001: the reference's endpoint re-eval was not carried over for the non-FSAL RK4 path.
- Verdict: **PASS WITH FIXES** — the integrator's *state* evolution is correct and the production DP45 path is endpoint-consistent. Fix R-001 so the rk4 path reports correct fault observables; address R-002/R-003 as quality follow-ups.

## Unreviewed Areas
- **Runtime stability of the RK CFL calibration** (`RkCflFactor=3/(2N+1)`; `cfl_mixed_flux_factor` RK arm 0.6/0.7) — in-source comment and project memory both state these are a *starting calibration validated on Frontera*, not by a local unit test. Effective CFLs at N=1/cfl=0.5 are None=0.5, Adjacent=0.30, AllContinuous=0.35 (math checks vs DRDG3D anchors); the "runs to t_final without runaway" claim (Phase 14.5 AC) cannot be verified locally (no full-mesh runs per project policy).
- **ADER byte-exactness** end-to-end (`make test` green + scalar dry-run bit-identical) was reasoned-through (additive, CLI-gated, default `cfl_rk_aware_=false`) but not re-run here.
- **`ApplyGradualOverstressCompactCircularAbsolute` / `InstantaneousOverstressCircular::ApplyAbsolute`** read and look correct (strike-only set; idempotent) but are not hit by T5 (Gaussian only). The instantaneous kind is LSW-only (TPV31) ⇒ unreachable on the RS-only RK path ⇒ benign.
- Single-precision build (`MFEM_USE_SINGLE`) not exercised.
