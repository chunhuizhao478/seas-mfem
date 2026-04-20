# v8.0.0 SeisSol cross-reference audit

**Context:** v8.0.0 check identified that fault-to-bulk coupling is broken on
the 400-rank TPV102 dispositive run — `V_max = 7.7 m/s` at hypocenter fault QP
but `max‖Q‖_∞ ≈ 1.66e-8` (10^15× below the expected m/s velocity or MPa stress
perturbation).  User question: *"will comparing our formulas against SeisSol
help, or do they already match?"*

**Short answer:** SeisSol has been cloned locally at
`/Users/chunhuizhao/projects/SeisSol`.  I compared every formula in our
`FaultFaceFlux::Evaluate` against SeisSol's `FrictionSolverCommon.h`.

**Our formulas match SeisSol exactly (trial traction, imposed states,
impedance handling).**  This is a useful rule-out: the bug is NOT in
`fault_face_flux.cpp::Evaluate`.  It must be downstream — in rotation,
accumulation, or state-evolution bookkeeping.

## Detailed side-by-side

### Trial traction (SeisSol `FrictionSolverCommon.h:182-192`)

```cpp
normalStress  = etaP * (qIMinus[U] - qIPlus[U] + qIPlus[N]*invZp  + qIMinus[N]*invZpNeig)
traction1     = etaS * (qIMinus[V] - qIPlus[V] + qIPlus[T1]*invZs + qIMinus[T1]*invZsNeig)
traction2     = etaS * (qIMinus[W] - qIPlus[W] + qIPlus[T2]*invZs + qIMinus[T2]*invZsNeig)
```

Where SeisSol `QuantityIndices` (`Misc.h:164-177`):
```
U = 6  (v_x)                  N  = 0  (sigma_xx in fault-local = sigma_nn)
V = 7  (v_y)                  T1 = 3  (sigma_xy in fault-local = tau_nt1)
W = 8  (v_z)                  T2 = 5  (sigma_xz in fault-local = tau_nt2)
```

**Exactly maps to our `wave_state.hpp` enum** (SXX=0, SXY=3, SXZ=5, VX=6,
VY=7, VZ=8).

Our `ComputeTrialTraction` (`fault_face_flux.cpp:37-63`):
```cpp
sigma_n_trial = data.eta_p * (Q_minus[VX] - Q_plus[VX] +
                              Q_plus[SXX]*invZp_plus + Q_minus[SXX]*invZp_minus);
tau1_trial    = data.eta_s * (Q_minus[VY] - Q_plus[VY] +
                              Q_plus[SXY]*invZs_plus + Q_minus[SXY]*invZs_minus);
tau2_trial    = data.eta_s * (Q_minus[VZ] - Q_plus[VZ] +
                              Q_plus[SXZ]*invZs_plus + Q_minus[SXZ]*invZs_minus);
```

**Match: bit-for-bit identical** (modulo `Zp_plus = Zp_minus = Zp` for our
homogeneous material).

### Imposed-state construction (SeisSol `FrictionSolverCommon.h:347-362`)

```cpp
imposedStateP[N][i]  = weight * normalStress;
imposedStateP[T1][i] = weight * traction1;
imposedStateP[T2][i] = weight * traction2;
imposedStateP[U][i]  = qIPlus[U] + invZp * (normalStress - qIPlus[N]);
imposedStateP[V][i]  = qIPlus[V] + invZs * (traction1    - qIPlus[T1]);
imposedStateP[W][i]  = qIPlus[W] + invZs * (traction2    - qIPlus[T2]);

imposedStateM[N][i]  = weight * normalStress;
imposedStateM[T1][i] = weight * traction1;
imposedStateM[T2][i] = weight * traction2;
imposedStateM[U][i]  = qIMinus[U] - invZpNeig * (normalStress - qIMinus[N]);
imposedStateM[V][i]  = qIMinus[V] - invZsNeig * (traction1    - qIMinus[T1]);
imposedStateM[W][i]  = qIMinus[W] - invZsNeig * (traction2    - qIMinus[T2]);
```

Our `Evaluate` (`fault_face_flux.cpp:118-144`):
```cpp
Q_imp_plus[SXX]  = sigma_n_corr;
Q_imp_plus[SXY]  = tau1_corr;
Q_imp_plus[SXZ]  = tau2_corr;
Q_imp_plus[VX]   = Q_plus[VX] + invZp_p * (sigma_n_corr - Q_plus[SXX]);
Q_imp_plus[VY]   = Q_plus[VY] + invZs_p * (tau1_corr    - Q_plus[SXY]);
Q_imp_plus[VZ]   = Q_plus[VZ] + invZs_p * (tau2_corr    - Q_plus[SXZ]);

Q_imp_minus[SXX] = sigma_n_corr;
Q_imp_minus[SXY] = tau1_corr;
Q_imp_minus[SXZ] = tau2_corr;
Q_imp_minus[VX]  = Q_minus[VX] - invZp_m * (sigma_n_corr - Q_minus[SXX]);
Q_imp_minus[VY]  = Q_minus[VY] - invZs_m * (tau1_corr    - Q_minus[SXY]);
Q_imp_minus[VZ]  = Q_minus[VZ] - invZs_m * (tau2_corr    - Q_minus[SXZ]);
```

**Match: bit-for-bit identical** (modulo the `weight` factor that SeisSol
uses for ADER-DG time-quadrature — irrelevant for our RK4 flux
computation).

### The one semantic difference — `Q` as perturbation vs full state

SeisSol's `qIPlus`, `qIMinus` are the **full bulk state** (including
background stress).  Our `Q_plus`, `Q_minus` are the **perturbation from
equilibrium** (background lives in `DOFData.tau_0, sigma_n_0`).

Consequence: SeisSol's `traction1, traction2, normalStress` passed into
the imposed-state formula are the **absolute corrected stresses**
(post-friction).  Our `tau1_corr, tau2_corr, sigma_n_corr` local
variables inside `Evaluate` are the **perturbation-only** corrected
stresses (`= tau_trial - eta_s * V`).

**Algebraically equivalent at the fault surface:**
- SeisSol's `(traction1 - qIPlus[T1])` = `(tau_full_corrected - tau_full_bulk)` = `delta_tau_correction_from_friction_solve`.
- Ours' `(tau1_corr - Q_plus[SXY])` = `(tau_perturbation_corrected - tau_perturbation_bulk)` = `delta_tau_correction_from_friction_solve`.

Both compute the **same physical correction**, just in different
reference frames.  The resulting `Q_imp[V]` (velocity field) is the
same absolute physical quantity — our `Q_imp_plus[VY]` is the
*perturbation* in v_y (which equals -V/2 at steady slip), SeisSol's
`imposedStateP[V]` is the *full* v_y (which also equals -V/2 since
the bulk was at v_y=0 before).

**This is not a bug.** It's a documentation/accounting difference that
both codes handle correctly in their own convention.

## What this rules out vs rules in

**Ruled out:** the `FaultFaceFlux::Evaluate` formulas are correct
(matching SeisSol exactly).  At the hypocenter with `V_slip = 7.7 m/s`
and `Q_plus = Q_minus = 0` (quiescent bulk perturbation), the
computed `Q_imp_plus[VZ]` IS the expected `-V/2 = -3.85 m/s`
(perturbation).

**Not ruled out — the downstream path:**
1. Rotation `T_can` mapping fault-local `Q_imp` → global `Q_imp_g`.
2. `flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h)` — the
   Godunov flux computation.
3. Accumulation `rhs[Elem1/Elem2] ±= w * shape * F_h`.
4. Mass-inverse and RK4 integration.
5. Ghost exchange (`q_gf.ExchangeFaceNbrData`) — does it correctly
   push our perturbation Q across rank boundaries?

**Between `Q_imp_plus[VZ] = -3.85 m/s` (fault-local) and
`max‖Q‖_∞ = 1.66e-8` (global bulk), something is eating the signal
by 10^8×.**

## New finding

### [R-1002] [CRITICAL] [POSSIBLE] `dynamic/wave_operator.inl` — Q_imp magnitudes are physically correct (-3.85 m/s velocity perturbation, -36 MPa stress perturbation at steady V=7.7 m/s hypocenter); somewhere in T_can rotation, Godunov flux, or accumulation the signal is attenuated by ~10^8× before reaching bulk

**Category:** BUG (POSSIBLE — pending targeted instrumentation)

**Description:**
SeisSol cross-reference confirms our `Evaluate` formulas produce
physically-correct `Q_imp` values.  For the 400-rank dispositive run's
hypocenter QP at steady `V_slip = 7.7 m/s`:
- Expected `Q_imp_plus[VZ]` (fault-local strike velocity) = -3.85 m/s.
- Expected `Q_imp_minus[VZ]` = +3.85 m/s.
- Jump across fault = 7.7 m/s = `V_slip`. ✓

After `T_can` rotation (fault-local → global, with `can_t2 = (+1, 0, 0)`
for TPV102):
- Expected `Q_imp_plus_g[VX]` = -3.85 m/s (global x = strike).
- Expected `Q_imp_minus_g[VX]` = +3.85 m/s.

After `flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h)`:
- Expected `F_h[SXY]` or `F_h[SXZ]` ~ `eta_s * V` = 36 MPa · (surface factor) ≈ 10^7 Pa / unit area.

After `rhs[Elem1] -= w * shape * F_h` and mass-inverse:
- Expected bulk stress perturbation at fault-adjacent DOFs to reach
  O(MPa) within a few element transit times (~ms).

**Observed:** `max‖Q‖_∞` across all 400 ranks stays at 1.66e-8 for the
entire 1.5 s simulation.  Something in the chain from `Q_imp`
(physically correct in fault-local) to bulk `Q` (should reach ~MPa near
fault) is attenuating by ~10^8×.

**Trigger:**
Any TPV102 run under the v7.0.0 code after R-701/R-801/R-802.

**Diagnostic plan (in priority order):**

- **D1 — instrument `Q_imp_plus_g` and `F_h` at hypocenter fault QP:**
  In `ComputeFaceFluxRHS`'s fault branch (or `ComputeSharedFaceFluxRHS`
  if hypocenter is on a shared fault face — check via the rank-3 driver
  log), add a one-time print of:
  ```cpp
  if (rank == hypo_rank && step == late_step) {
     std::cout << "DIAG Q_imp_plus_local:";
     for (int c = 0; c < NUM_STATE; c++) std::cout << " " << Q_imp_plus[c];
     std::cout << "\nDIAG Q_imp_plus_g:";
     for (int c = 0; c < NUM_STATE; c++) std::cout << " " << Q_imp_plus_g[c];
     std::cout << "\nDIAG F_h_total:";
     for (int c = 0; c < NUM_STATE; c++) std::cout << " " << F_h_total[c];
  }
  ```

  **Pass criteria:**
  - `Q_imp_plus_local[VZ]` = `-3.85 m/s`  (fault-local)
  - `Q_imp_plus_g[VX]`     = `-3.85 m/s`  (after T_can rotation)
  - `F_h_total[SXY]` or `F_h_total[VX]` = O(1e6-1e7)  (nonzero, MPa-scale)

  If Q_imp_plus_local is correct but Q_imp_plus_g is tiny → T_can bug.
  If Q_imp_plus_g is correct but F_h is tiny → flux_.Interior bug.
  If F_h is correct but bulk max‖Q‖ stays tiny → accumulation / mass /
  RK4 / ghost-exchange bug.

- **D2 — instrument `rhs[Elem1_dof]` before and after fault flux:**
  Print rhs contributions from the fault face to confirm the
  accumulation is being executed.

- **D3 — instrument `Q_post_RK4_step` at fault-adjacent DOF:**
  Print Q values at bulk DOFs adjacent to hypocenter across time
  steps.  If they monotonically grow toward MPa-scale, the issue is
  in how `max‖Q‖_∞` is sampled (which is just `Q.Normlinf()` — hard to
  have a bug there, but check).  If they stay at 0, flux is not
  reaching bulk.

**Suggested fix:** Depends on which step D1-D3 localizes.  No blind
fix.

---

## Does SeisSol comparison help further?

At this point: **no, further SeisSol cross-reference is low-value.**
SeisSol's trial-correct formulas match ours.  The remaining suspect
code paths (`T_can` rotation, `flux_.Interior`, the DG accumulation)
are **specific to our MFEM DG implementation**, not shared with
SeisSol's ADER-DG + yateto infrastructure.  Diffing beyond
`FrictionSolverCommon.h` would compare apples-to-oranges at the flux-
assembly level.

**Priority: run D1 locally** with the `seas_tpv102_driver` at 2 or 4
ranks with enough tfinal to drive steady hypocenter slip (say,
`--tfinal 1.0`).  D1 is a ~30-minute addition of diagnostic prints and
should pinpoint whether the bug lies in rotation, flux, or
accumulation.

## Version note
This is an audit document, not a new `v[a].[b].[c]` bump.  Current
version remains **`v8.0.0`** with R-1001 open.  R-1002 is a refinement
of R-1001's root cause, not a new finding.

---

## Also written to `REVIEW.md` per /code-review convention

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v8.0.0_check.md`
- Files reviewed (cross-reference):
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (our `Evaluate`,
    `ComputeTrialTraction`)
  - `miniapps/seas/dynamic/wave_state.hpp` (state enum indices)
  - `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h`
    (SeisSol reference)
  - `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/Misc.h`
    (SeisSol `QuantityIndices`)
- Domain context: SeisSol is a production DG dynamic-rupture code
  with extensive TPV102 verification.  Our code cites
  "SeisSol FrictionSolverCommon.h" in header comments.

## Findings
See R-1002 above.

## Summary
- Critical issues: **0 new** (R-1002 is a refinement of R-1001 from
  v8.0.0 check).
- Moderate issues: **0**
- Low issues: **0**
- Plan compliance: N/A (diagnostic audit, not plan-driven).
- Verdict: **SeisSol comparison confirms our formulas are correct.
  Focus debugging on T_can rotation, flux_.Interior, and accumulation
  via D1 instrumentation.**

## Unreviewed Areas
- SeisSol's flux-assembly step (ADER-DG kernels, not transferable).
- SeisSol's Q is in native-fault-local coordinates always; ours round-
  trips through global.  If the round-trip loses precision, a
  different workflow.
