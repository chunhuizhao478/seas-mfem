# Implementation Plan: TPV205 LSW-Native Fields & TPV205-Aware ADER Fault Dispatch

**Date:** 2026-04-27
**Owner:** SEAS-MFEM TPV205 dynamic-rupture follow-up
**Author:** code-plan agent (under user authority — see REVIEW.md R-016)
**Consumed by:** /code-implement, /code-fix, /code-review

## Overview

Replace the TPV205 `DOFData` field repurposing (`data.a ← μ_s`,
`data.psi ← μ_d`, `data.Dc ← d_c`) with **LSW-native fields** on `DOFData`,
add a **TPV205-aware ADER dispatch** in the wave operator, and remove the
`np > 1` abort proposed as a stop-gap by REVIEW round-3 (R-016). The redesign
restores correct rupture physics across MPI rank boundaries while preserving
TPV102/TPV104 byte-for-byte. Every phase is additive: at the end of each
phase the project still builds and the existing test suite passes.

## Constraints

### Interface constraints — what cannot change
- `FaultFaceFlux::Evaluate`, `EvaluateTotal`, `EvaluateADER`, `EvaluateADERTotal`,
  `ComputeStageState`, `BuildImposedState`, `WriteBackState` signatures and
  semantics. TPV102/TPV104 depend on these byte-for-byte.
- `DOFData` field layout and offsets for all PRE-EXISTING fields. New
  fields **append** to the struct; pre-existing fields stay where they
  are. Any reordering would break the `R-V92-H07` psi-invariance assertion
  in `Evaluate` (which captures `data.psi` by value at entry), the
  TPV104 substep iterator's reads of `data.a / data.psi / data.Dc`, and
  the `EvaluateTotal` zero-pre-stress contract assertion.
- `WaveOperator<MeshType>` template public API — only **add** the
  `SetFaultFrictionLaw` setter and `FaultFrictionLaw` enum. Do not
  rename or remove members.
- The R-1600/R-1601 frame-mismatch fallback for **rate-and-state** shared
  faces stays exactly as it is. Only the LSW branch is new.

### Dependency constraints — what must be used
- `LSWFrictionCoefficient_TPV205` (in `dynamic/tpv205_friction.hpp`) —
  the canonical μ(δ) helper. The R-002 strength-barrier short-circuit
  lives there; do not duplicate.
- `SolveLSW_TPV205` — the canonical closed-form V solve. The R-003
  barrier short-circuit lives there; do not duplicate.
- `FaultFaceFlux::ComputeTrialTraction` (static, pure) — must be the
  trial-traction source for both the iterator and the new
  `EvaluateADER_LSW`.
- `FaultFaceFlux::BuildImposedState` (member, pure on `EvalStageState`) —
  must be the imposed-state builder for both the iterator and
  `EvaluateADER_LSW`.
- `FaultFaceFlux::WriteBackState` (member) — must be the post-call
  data-writeback path for `EvaluateADER_LSW` so `data.tau*_corr /
  sigma_n_corr / V*` carry TOTAL physical traction in the same convention
  as `EvaluateADER`.

### Convention constraints — patterns to follow
- BP5/Tandem canonical fault-local frame: `tangent1 = dip,
  tangent2 = strike` (per `miniapps/seas/CLAUDE.md` "Sign Conventions").
  The new fields hold material parameters only — frame-independent.
- Sign convention: `sigma_n > 0` = compression. Already encoded in
  `SolveLSW_TPV205`; preserved by the new dispatch.
- `real_t` is `double` in production builds and may be `float` under
  `-DMFEM_USE_SINGLE`. All literals comparing against the new fields
  use `static_cast<real_t>(...)` or `0.0` (which is fine — implicit
  promotion is OK for zero defaults).
- Defensive design: every change is **additive** — default behaviour
  (`FaultFrictionLaw::RateAndState`) is byte-identical to pre-change
  builds. The new LSW dispatch is opt-in via the new setter.

### Numerical / performance constraints
- ADER `EvaluateADER_LSW` must be O(dt²)-consistent with the iterator's
  per-sub-step LSW solve at O = 1 (the trivial coincidence case where
  `tau_node = dt/2`). Proof outline (mirrors the pre-existing
  `EvaluateADER` derivation): for a smooth Q(τ) on [0, dt],
  `Q̄ = ∫Q dτ / dt = Q(dt/2) + O(dt²)` (Simpson on smooth Q); LSW closed
  form is Lipschitz in (τ_total, σ_n) so the imposed Q is O(dt²)-consistent
  with feeding Q at the midpoint. The wrap pattern `Q̄ = I/dt; Q_imp =
  LSW_solve(Q̄); I_imp = Q_imp · dt` is identical to the rate-and-state
  variant.
- The new dispatch branch in `wave_operator.inl` must not regress the
  per-call cost on the rate-and-state path: a single `if (law_ ==
  FaultFrictionLaw::LSW)` branch is below the noise floor of one
  Brent iteration.

---

## Phase 1: Add fields + dispatch enum + EvaluateADER_LSW (no callers yet)

### Goal
After this phase the code compiles, all existing tests pass byte-identically,
and three new entry points exist but are **never called** by any existing
code path: (a) `DOFData::lsw_mu_s / lsw_mu_d / lsw_d_c` fields with zero
defaults, (b) `FaultFrictionLaw` enum + `WaveOperator::SetFaultFrictionLaw`
setter (with default `RateAndState`), (c)
`FaultFaceFlux::EvaluateADER_LSW(...)`.

### Files to Modify

- `dynamic/fault_face_flux.hpp` — append three LSW-native fields to `DOFData`;
  declare `EvaluateADER_LSW`.
- `dynamic/fault_face_flux.cpp` — implement `EvaluateADER_LSW` mirroring
  `EvaluateADER`'s I/dt wrap, but using the LSW closed-form on the
  time-averaged Q.
- `dynamic/wave_operator.hpp` — declare `enum class FaultFrictionLaw {
  RateAndState, LSW };` and add a setter + private member
  `fault_friction_law_` with default `RateAndState`.

### Files to Create
None.

### Detailed Requirements

**1.1 — DOFData field append (`dynamic/fault_face_flux.hpp`).**

Insert the following block in `DOFData`, **after** the existing
`tau1_corr / tau2_corr / sigma_n_corr` block (around the current
end of the struct, just before any `#ifdef SEAS_DIAG_FAULT_FLUX
bool diag_print` field). Keep these fields *appended* so no existing
struct offset moves.

```cpp
   // -----------------------------------------------------------------------
   // Linear slip-weakening (LSW) parameters — separate slot from the
   // rate-and-state `a / psi / Dc` so that:
   //   1. A code path that reads `data.a / data.psi / data.Dc` for
   //      rate-and-state friction (Evaluate, EvaluateTotal, the TPV104
   //      iterator, every Brent / Newton call) cannot accidentally
   //      consume LSW values and produce off-by-orders-of-magnitude
   //      strengths (REVIEW R-016).
   //   2. A code path that reads these LSW fields cannot accidentally
   //      consume rate-and-state values: zero defaults make
   //      `LSWFrictionCoefficient_TPV205(δ, 0, 0, 0)` deterministic.
   //
   // Populated only by InitializeFaultDOFs_TPV205; pre-stress / nucleation /
   // impedance fields above are shared with rate-and-state callers.
   real_t lsw_mu_s = 0.0;   ///< LSW static friction μ_s (≥ mu_s_barrier ⇒ barrier QP)
   real_t lsw_mu_d = 0.0;   ///< LSW dynamic friction μ_d
   real_t lsw_d_c  = 0.0;   ///< LSW slip-weakening critical distance d_c [m]
```

The `SEAS_DIAG_FAULT_FLUX bool diag_print` field, if present, must remain
the LAST field so `sizeof(DOFData)` only changes by `3 * sizeof(real_t)`
in builds without the diag flag, and by the same in builds with it.

**1.2 — Declare `EvaluateADER_LSW` in `dynamic/fault_face_flux.hpp`.**

Insert immediately after the existing `EvaluateADERTotal` declaration:

```cpp
   /// @brief LSW counterpart to `EvaluateADER` — time-integrated Riemann
   /// solve for the linear slip-weakening law (SCEC TPV5 / TPV205).
   ///
   /// I-form wrap pattern identical to `EvaluateADER`:
   ///   Q̄± = I±/dt, run the LSW closed form on Q̄±, scale I_imp = dt·Q_imp.
   /// The friction physics is the SCEC TPV5 §7-11 closed form:
   ///   μ(δ) via `LSWFrictionCoefficient_TPV205`,
   ///   V_abs / V1 / V2 / τ*_corr via `SolveLSW_TPV205`,
   /// reading the LSW-native fields `data.lsw_mu_s / lsw_mu_d / lsw_d_c`.
   /// `data.a / data.psi / data.Dc` are NOT read.
   ///
   /// On return `data.{slip_rate, V1, V2, tau1_corr, tau2_corr,
   /// sigma_n_corr}` carry the time-averaged values over [t_n, t_n+dt],
   /// matching `EvaluateADER`'s station-output convention; `data.tau*_corr`
   /// and `data.sigma_n_corr` are TOTAL (= pre-stress + nuc + trial-scale
   /// corrected) per `WriteBackState`. `data.psi` is NOT touched.
   ///
   /// `data.slip{1,2}` accumulation is the SAME single-step forward-Euler
   /// step the iterator applies — `data.slip{1,2} += V{1,2} * dt` —
   /// so a one-shot driver path stays slip-consistent with the substep
   /// iterator path at O = 1.  (Multi-substep callers should use the
   /// iterator instead; `EvaluateADER_LSW` is a one-shot per-call helper.)
   ///
   /// @param[in,out] data   Per-DOF state.
   /// @param[in]  I_plus    Time-integrated + side state (NUM_STATE).
   /// @param[in]  I_minus   Time-integrated − side state (NUM_STATE).
   /// @param[in]  dt        Time step (> 0).
   /// @param[out] I_imp_plus   Time-integrated imposed + state (NUM_STATE).
   /// @param[out] I_imp_minus  Time-integrated imposed − state (NUM_STATE).
   void EvaluateADER_LSW(DOFData &data,
                         const real_t *I_plus, const real_t *I_minus,
                         real_t dt,
                         real_t *I_imp_plus, real_t *I_imp_minus) const;
```

No `FrictionSolver::Method` parameter — LSW has no root finder.

**1.3 — Implement `EvaluateADER_LSW` in `dynamic/fault_face_flux.cpp`.**

Append at the end of the file, before the closing namespace braces. Do
NOT touch any other function body in this file.

```cpp
// ---------------------------------------------------------------------------
// REVIEW R-016: LSW counterpart to EvaluateADER.  Wraps the I/dt → I·dt
// pattern around the closed-form LSW solve so the wave operator's
// shared-fault dispatch path (wave_operator.inl R-1600 fallback) can
// route TPV205 through correct physics instead of Brent on the
// rate-and-state law.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADER_LSW(DOFData &data,
                                     const real_t *I_plus,
                                     const real_t *I_minus,
                                     real_t dt,
                                     real_t *I_imp_plus,
                                     real_t *I_imp_minus) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADER_LSW: dt must be > 0, got " << dt);

   // Homogeneous-material check (mirrors Evaluate / EvaluateTotal).  The
   // v9.0.0 Pelties-9 per-side flux assumes A_plus == A_minus.
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             static_cast<real_t>(1e-12) * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "FaultFaceFlux::EvaluateADER_LSW: bimaterial fault face "
               "(Zp_plus=" << data.Zp_plus << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus << " Zs_minus=" << data.Zs_minus
               << ").  Extend per-side handling before running this configuration.");

   // Step 0: Q̄± = I±/dt.
   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = static_cast<real_t>(1.0) / dt;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      Q_avg_plus[c]  = I_plus[c]  * inv_dt;
      Q_avg_minus[c] = I_minus[c] * inv_dt;
   }

   // Step 1: trial traction (pure helper).
   EvalStageState s;
   ComputeTrialTraction(data, Q_avg_plus, Q_avg_minus,
                        s.sigma_n_trial, s.tau1_trial, s.tau2_trial);

   // Step 2: total traction (TPV205 has zero nucleation channels by
   // contract; pre-stress lives in DOFData).
   s.sigma_n_total = data.sigma_n0 + data.sigma_n_nuc + s.sigma_n_trial;
   s.tau1_total    = data.tau1_0   + data.tau1_nuc    + s.tau1_trial;
   s.tau2_total    = data.tau2_0   + data.tau2_nuc    + s.tau2_trial;
   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
                              + s.tau2_total * s.tau2_total);

   // Step 3: μ_eff(δ) at the slip magnitude carried in DOFData (the
   // iterator updates δ each substep; this one-shot call sees the
   // current cumulative slip).  Reads ONLY the LSW-native fields.
   const real_t delta = std::sqrt(data.slip1 * data.slip1
                                  + data.slip2 * data.slip2);
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       data.lsw_mu_s,
                                                       data.lsw_mu_d,
                                                       data.lsw_d_c);

   // Step 4: closed-form LSW solve (R-002/R-003 barrier short-circuits live
   // inside SolveLSW_TPV205).  Sets s.V_abs, s.V1, s.V2, s.tau{1,2}_corr.
   SolveLSW_TPV205(s.tau1_trial, s.tau2_trial,
                   s.tau1_total, s.tau2_total,
                   s.sigma_n_total, data.eta_s,
                   mu_eff,
                   s.V_abs, s.V1, s.V2,
                   s.tau1_corr, s.tau2_corr);

   // σ_n is unaffected by friction — TRIAL-scale value matches the
   // rate-and-state path's CompleteFromVabs convention.
   s.sigma_n_corr = s.sigma_n_trial;

   // Step 5: forward-Euler slip accumulation, identical to the iterator's
   // per-substep update at O = 1.  Multi-substep callers route through the
   // iterator instead and do NOT call EvaluateADER_LSW.
   data.slip1 += s.V1 * dt;
   data.slip2 += s.V2 * dt;

   // Step 6: imposed Q-state (Eq. 11-12 of the FaultFaceFlux pipeline).
   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   BuildImposedState(data, s, Q_avg_plus, Q_avg_minus,
                     Q_imp_plus, Q_imp_minus);

   // Step 7: write back V/slip_rate/τ*_corr/σ_n_corr to DOFData
   // (TOTAL physical traction — pre + nuc + trial-scale corrected).
   WriteBackState(data, s);

   // Step 8: rescale Q_imp back to time-integrated form.
   for (int c = 0; c < NUM_STATE; ++c)
   {
      I_imp_plus[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus[c] = Q_imp_minus[c] * dt;
   }
}
```

Note: this function intentionally DOES include slip accumulation
(matching iterator's `StepOneQP_`). The wave-operator dispatch (Phase 2)
must call this **at most once per macro-step per QP** — if a future
substep loop runs `EvaluateADER_LSW` per substep without disabling
slip accumulation, the slip would over-integrate. The Phase 2
dispatch always runs it once per macro-step in the shared-fault
branch, so this is correct as-is.

**1.4 — `FaultFrictionLaw` enum + setter in `dynamic/wave_operator.hpp`.**

Add immediately above the existing `WaveOperator` class definition
(near line 70, alongside the other `enum class` declarations like
`FreeSurfaceBCMode` and `MixedFluxMode`):

```cpp
/// REVIEW R-016: friction-law tag the driver sets at init so the wave
/// operator's fault dispatch (interior + shared) can route to the
/// correct ADER closure.  Default `RateAndState` keeps TPV102 / TPV104
/// byte-identical; `LSW` routes through `FaultFaceFlux::EvaluateADER_LSW`
/// for TPV205.
enum class FaultFrictionLaw { RateAndState, LSW };
```

Inside `class WaveOperator`, immediately after `SetFaultFlux` /
`GetFaultFlux` (around current line 112), add:

```cpp
   /// REVIEW R-016: select the ADER fault dispatch.  Default
   /// `FaultFrictionLaw::RateAndState` runs `fault_flux_->EvaluateADER` /
   /// `EvaluateADERTotal`; `FaultFrictionLaw::LSW` runs
   /// `fault_flux_->EvaluateADER_LSW` and is required by TPV205.
   void SetFaultFrictionLaw(FaultFrictionLaw law)
   { fault_friction_law_ = law; }
   FaultFrictionLaw GetFaultFrictionLaw() const
   { return fault_friction_law_; }
```

Add to the private members (around current line 602, alongside
`fault_flux_`):

```cpp
   FaultFrictionLaw fault_friction_law_ = FaultFrictionLaw::RateAndState;
```

### Interfaces

- `EvaluateADER_LSW(DOFData &, const real_t*, const real_t*, real_t,
  real_t*, real_t*) const` — new; signature matches `EvaluateADER`
  minus the `FrictionSolver::Method` parameter (no root finder).
- `WaveOperator::SetFaultFrictionLaw(FaultFrictionLaw)` — new; setter
  called once at init by drivers.
- `enum class FaultFrictionLaw` — new; two values, `RateAndState`
  (default) and `LSW`.

### Edge Cases to Handle

- `dt ≤ 0` in `EvaluateADER_LSW`: `MFEM_VERIFY` aborts loudly (matches
  `EvaluateADER`).
- Bimaterial fault face (`Zp_plus != Zp_minus` or `Zs_plus != Zs_minus`):
  `MFEM_VERIFY` aborts loudly (matches `EvaluateADER`).
- All-zero LSW fields (legacy `DOFData` constructed by a non-TPV205
  caller):
  - `LSWFrictionCoefficient_TPV205(δ, 0, 0, 0)` enters the
    `mu_s >= 0.5 * mu_s_barrier` test, which is false (0 < 5000),
    so falls through. `δ <= 0` returns `mu_s = 0`. `δ > 0 && δ >= d_c = 0`
    returns `mu_d = 0`. So `mu_eff = 0` ⇒ `tau_strength = 0` ⇒ `V_abs =
    |τ_total|/η_s` (unconstrained sliding). This is the wrong
    physics for any non-LSW caller, but the dispatch in Phase 2
    only routes here when `fault_friction_law_ == LSW` — which
    requires the driver to have explicitly set the law. Still:
    add an `MFEM_VERIFY(data.lsw_mu_s > 0 || data.lsw_mu_d > 0 ||
    data.lsw_d_c > 0, ...)` at the top of `EvaluateADER_LSW` to
    catch the misuse early.

### Acceptance Criteria
- [ ] `make seas_test_fault_face_flux -j8` passes byte-identical to
  pre-change.
- [ ] `make seas_test_tpv104_substep_iterator -j8` passes byte-identical.
- [ ] `make seas_tpv102_driver -j8` and `make seas_tpv104_driver -j8`
  link without warnings (excluding the harmless `-rpath` ld warning).
- [ ] `nm seas_tpv102_driver | grep EvaluateADER_LSW` shows the new
  symbol present (it links into the binary even if not called).
- [ ] No existing test changes its output by even one bit.

### Dependencies
- Depends on: nothing (additive only).
- Required by: Phase 2.

---

## Phase 2: Wire dispatch enum into wave_operator.inl interior + shared-fault paths

### Goal
After this phase, the wave operator's fault branch routes to
`fault_flux_->EvaluateADER_LSW(...)` when `fault_friction_law_ ==
FaultFrictionLaw::LSW`. Default (`RateAndState`) is byte-identical to
pre-change. **No driver wires the LSW law yet** — TPV102 / TPV104 must
still produce identical output.

### Files to Modify

- `dynamic/wave_operator.inl` — add two additive `if (fault_friction_law_
  == FaultFrictionLaw::LSW)` branches:
  1. Interior fault path, around line 3614 (`fault_flux_->EvaluateADER`
     call inside the `else` branch of the substep-side-channel
     dispatch).
  2. Shared-fault R-1600 fallback, around line 4518
     (`fault_flux_->EvaluateADER` call after the SHARED FALLBACK
     comment block).

### Files to Create
None.

### Detailed Requirements

**2.1 — Interior fault path branch (around line 3614).**

Replace the existing block:

```cpp
                  else
                  {
                     // v9.4.0 Commit 3: fluctuation-Q ADER dispatch;
                     // has_bulk_bg_ already asserted at the top of
                     // ComputeADERFaceFluxRHS (Q_bg = 0 is valid).
                     fault_flux_->EvaluateADER(fdata,
                                               I_plus_local, I_minus_local,
                                               dt,
                                               I_imp_plus, I_imp_minus);
                  }
```

with:

```cpp
                  else
                  {
                     // v9.4.0 Commit 3: fluctuation-Q ADER dispatch;
                     // has_bulk_bg_ already asserted at the top of
                     // ComputeADERFaceFluxRHS (Q_bg = 0 is valid).
                     // REVIEW R-016: route LSW callers (TPV205) through
                     // EvaluateADER_LSW; the default RateAndState path
                     // is byte-identical to pre-change.
                     if (fault_friction_law_ == FaultFrictionLaw::LSW)
                     {
                        fault_flux_->EvaluateADER_LSW(
                           fdata,
                           I_plus_local, I_minus_local,
                           dt,
                           I_imp_plus, I_imp_minus);
                     }
                     else
                     {
                        fault_flux_->EvaluateADER(fdata,
                                                  I_plus_local, I_minus_local,
                                                  dt,
                                                  I_imp_plus, I_imp_minus);
                     }
                  }
```

**2.2 — Shared-fault R-1600 fallback branch (around line 4518).**

Replace the existing block:

```cpp
                  // SHARED FALLBACK: this branch always runs the inline
                  // EvaluateADER regardless of substep_I_imp_*_flat_.
                  fault_flux_->EvaluateADER(fdata,
                                            I_plus_local, I_minus_local,
                                            dt,
                                            I_imp_plus, I_imp_minus);
                  (void)substep_I_imp_plus_flat_;
                  (void)substep_I_imp_minus_flat_;
                  (void)substep_n_total_fault_qps_;
```

with:

```cpp
                  // SHARED FALLBACK: this branch always runs the inline
                  // ADER closure regardless of substep_I_imp_*_flat_
                  // (R-1600/R-1601 frame-mismatch on shared QPs).
                  // REVIEW R-016: dispatch on the friction-law tag —
                  // RateAndState keeps the original Brent path
                  // byte-identical (TPV102/TPV104); LSW (TPV205) runs
                  // the closed-form solver.  Without this branch, TPV205
                  // shared-fault QPs at np > 1 silently consume LSW
                  // values via Brent and stall the rupture front.
                  if (fault_friction_law_ == FaultFrictionLaw::LSW)
                  {
                     fault_flux_->EvaluateADER_LSW(
                        fdata,
                        I_plus_local, I_minus_local,
                        dt,
                        I_imp_plus, I_imp_minus);
                  }
                  else
                  {
                     fault_flux_->EvaluateADER(fdata,
                                               I_plus_local, I_minus_local,
                                               dt,
                                               I_imp_plus, I_imp_minus);
                  }
                  (void)substep_I_imp_plus_flat_;
                  (void)substep_I_imp_minus_flat_;
                  (void)substep_n_total_fault_qps_;
```

**2.3 — Audit no other `fault_flux_->Evaluate*(` call sites (≤ ~10).**

Run `grep -n "fault_flux_->Evaluate" miniapps/seas/dynamic/wave_operator.inl`
and confirm all sites that participate in the ADER fault dispatch are
covered. Sites at lines ~2419, ~3058 are the **non-ADER**
`fault_flux_->Evaluate(...)` (single-stage RK4 dispatch) — those are
NOT touched in this phase. They are not on TPV205's hot path
(TPV205 is ADER-only per `tpv205_driver.cpp` line 1539-1544).
**Document this in the diff as a comment** so future TPV205 RK4
wiring knows to extend dispatch.

### Interfaces
None changed — only internal dispatch.

### Edge Cases to Handle

- `fault_friction_law_ == FaultFrictionLaw::RateAndState` (default):
  branch falls through to existing `EvaluateADER` — byte-identical.
- `fault_flux_ == nullptr`: pre-existing path already segfaults; no new
  guard added (keeping behaviour symmetric with existing code).

### Acceptance Criteria
- [ ] `make all -j8` succeeds without new warnings.
- [ ] **TPV102 binary byte-identical** (compare `md5sum
  seas_tpv102_driver` pre/post; or, if compiler-flag rebuild changes
  hash, run `seas_test_ader_tpv102_smoke` and confirm output
  identical to a pre-change golden).
- [ ] **TPV104 binary byte-identical** (run
  `seas_test_tpv104_substep_iterator_parity`,
  `seas_test_tpv104_substep_one_shot_parity`,
  `seas_test_tpv104_substep_dispatch_parity`;
  all must pass with `Δ = 0` to floating-point bit precision).
- [ ] BP5 build still links (`make seas_bp5_full -j8`) and
  `make test-bp5-smoke` passes.
- [ ] `grep -n "fault_friction_law_" miniapps/seas/dynamic/wave_operator.inl`
  returns exactly two hits (the two new `if` branches).

### Dependencies
- Depends on: Phase 1 (uses `EvaluateADER_LSW` and `FaultFrictionLaw`).
- Required by: Phase 3.

---

## Phase 3: TPV205 driver / setup / iterator / station / ParaView migrate to LSW-native fields

### Goal
After this phase, every TPV205 code path reads / writes the new LSW-native
fields (`d.lsw_mu_s / lsw_mu_d / lsw_d_c`) instead of the rate-and-state
slots (`d.a / d.psi / d.Dc`). The driver registers `FaultFrictionLaw::LSW`
with the wave operator and removes the round-3 R-016 stop-gap (if it was
applied; if not, this phase establishes the same outcome by other means).
TPV205 4-rank smoke (4 ADER-O2 steps) reproduces V_max ≈ 0.0779 m/s, and
TPV102 / TPV104 are still byte-identical (they don't touch the new fields).

### Files to Modify

- `dynamic/tpv205_setup.hpp`
  - `InitializeFaultDOFs_TPV205` — populate `d.lsw_mu_s / lsw_mu_d /
    lsw_d_c`; explicitly zero `d.a / d.psi / d.Dc`.
  - `TPV205StationWriter::WriteStep` — read LSW-native fields; replace
    inline μ_eff with a single `LSWFrictionCoefficient_TPV205` call
    (closes R-017).
- `dynamic/tpv205_substep_iterator.cpp`
  - `Tpv205SubStepIterator::StepOneQP_` — read `d.lsw_mu_s / lsw_mu_d /
    lsw_d_c` instead of `d.a / d.psi / d.Dc`.
- `drivers/tpv205_driver.cpp`
  - Wire `wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW)` after
    `SetFaultFlux`.
  - ParaView setup loop (`pv_local_a / pv_local_Dc` populating, around
    line 1731) — read `d.lsw_mu_s / d.lsw_d_c`.
  - Per-step ParaView state channel (`pv_local_state(i) =
    LSWFrictionCoefficient_TPV205(delta, d.a, d.psi, d.Dc)`, around
    line 1931) — read `d.lsw_mu_s / d.lsw_mu_d / d.lsw_d_c`.
  - Remove the R-016 round-3 abort guard (if present) and replace
    the `nprocs > 1` banner with the corrected message.

### Files to Create
None.

### Detailed Requirements

**3.1 — `InitializeFaultDOFs_TPV205` (`dynamic/tpv205_setup.hpp`).**

Locate the existing block:

```cpp
      // LSW friction parameters (DOFData field repurposing).
      d.a   = ComputeMuS_TPV205(along_strike, down_dip);   // μ_s
      d.psi = ComputeMuD_TPV205(along_strike, down_dip);   // μ_d
      d.Dc  = ComputeDc_TPV205(along_strike, down_dip);    // d_c
```

Replace with:

```cpp
      // R-016: LSW-native fields.  Populates the new dedicated slots so
      // every LSW reader (iterator, station writer, ParaView, the new
      // EvaluateADER_LSW) consumes them directly — no field repurposing,
      // no rate-and-state code path can accidentally read these as a /
      // psi / Dc.
      d.lsw_mu_s = ComputeMuS_TPV205(along_strike, down_dip);  // μ_s
      d.lsw_mu_d = ComputeMuD_TPV205(along_strike, down_dip);  // μ_d
      d.lsw_d_c  = ComputeDc_TPV205(along_strike, down_dip);   // d_c

      // Defensive: zero the rate-and-state slots.  TPV205 must NEVER
      // run rate-and-state.  If a future regression sends TPV205 data
      // into Evaluate / EvaluateTotal, the FrictionSolver guard sees
      // a == 0 and aborts loudly rather than computing nonsense
      // (REVIEW R-016 was caused by Brent silently running on
      // repurposed values).
      d.a   = 0.0;
      d.psi = 0.0;
      d.Dc  = 0.0;
```

**3.2 — `TPV205StationWriter::WriteStep` (`dynamic/tpv205_setup.hpp`).**

Replace the inline μ_eff block at lines ~303-314 with:

```cpp
         const DOFData &d = dof_data[idx];

         // R-017 / R-016: route through the canonical helper so the
         // strength-barrier short-circuit applies uniformly across
         // station traces, ParaView, and the friction solve.  Reads
         // the LSW-native fields directly.
         const real_t delta = std::sqrt(d.slip1 * d.slip1
                                        + d.slip2 * d.slip2);
         const real_t mu_eff =
            LSWFrictionCoefficient_TPV205(delta,
                                          d.lsw_mu_s, d.lsw_mu_d,
                                          d.lsw_d_c);
```

Remove the now-dead `mu_s_qp / mu_d_qp / d_c_qp` locals. Add the include
`#include "tpv205_friction.hpp"` at the top of `tpv205_setup.hpp` if
not already present.

**3.3 — `Tpv205SubStepIterator::StepOneQP_` (`dynamic/tpv205_substep_iterator.cpp`).**

Replace lines ~99-105:

```cpp
   const real_t delta = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
   const real_t mu_s_qp = d.a;     // μ_s stored in `a`
   const real_t mu_d_qp = d.psi;   // μ_d stored in `psi`
   const real_t d_c_qp  = d.Dc;
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       mu_s_qp,
                                                       mu_d_qp,
                                                       d_c_qp);
```

with:

```cpp
   // R-016: read LSW-native fields directly — no field repurposing.
   const real_t delta = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       d.lsw_mu_s,
                                                       d.lsw_mu_d,
                                                       d.lsw_d_c);
```

**3.4 — `drivers/tpv205_driver.cpp`: wire `SetFaultFrictionLaw`.**

After the existing line:

```cpp
   wave.SetFaultFlux(&fault_flux);
```

(around line 1428 — locate via grep), add:

```cpp
   // R-016: tell the wave operator to dispatch the LSW closed-form on
   // the fault branch (interior + shared-fault fallback).  Without this,
   // shared-fault QPs at np > 1 fall through to fault_flux_->EvaluateADER
   // (Brent on rate-and-state), which produces wrong physics for LSW
   // and stalls the rupture front at MPI rank boundaries.
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
```

**3.5 — `drivers/tpv205_driver.cpp`: ParaView setup loop (around line 1731).**

Replace:

```cpp
      for (int i = 0; i < num_fault_total; i++)
      {
         pv_local_a(i)  = dof_data[i].a;
         pv_local_Dc(i) = dof_data[i].Dc;
         pv_local_x2(i) = fault_coords[i](0);
         pv_local_x3(i) = std::abs(fault_coords[i](2));
      }
```

with:

```cpp
      for (int i = 0; i < num_fault_total; i++)
      {
         // R-016: ParaView "a" / "Dc" channels carry LSW μ_s / d_c
         // (mirrors station-writer convention).  TPV104's fault output
         // uses d.a = direct-effect parameter; TPV205 reuses the SAME
         // channel slot for the LSW μ_s.
         pv_local_a(i)  = dof_data[i].lsw_mu_s;
         pv_local_Dc(i) = dof_data[i].lsw_d_c;
         pv_local_x2(i) = fault_coords[i](0);
         pv_local_x3(i) = std::abs(fault_coords[i](2));
      }
```

**3.6 — `drivers/tpv205_driver.cpp`: per-step ParaView state channel (around line 1928-1933).**

Replace:

```cpp
         {
            const real_t delta = std::sqrt(d.slip1 * d.slip1
                                           + d.slip2 * d.slip2);
            pv_local_state(i) = mfem::seas::LSWFrictionCoefficient_TPV205(
                                   delta, d.a, d.psi, d.Dc);
         }
```

with:

```cpp
         {
            // R-016: read LSW-native fields directly.
            const real_t delta = std::sqrt(d.slip1 * d.slip1
                                           + d.slip2 * d.slip2);
            pv_local_state(i) = mfem::seas::LSWFrictionCoefficient_TPV205(
                                   delta,
                                   d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
         }
```

**3.7 — `drivers/tpv205_driver.cpp`: remove the R-016 abort guard.**

If round-3 R-016's stop-gap was applied (an `MPI_Abort` or
`MFEM_ABORT` inside `if (nprocs > 1)`), remove it.

Replace any remaining `nprocs > 1` banner block with:

```cpp
   if (nprocs > 1 && rank == 0)
   {
      std::cout << "[tpv205_driver] MPI dispatch (nprocs=" << nprocs
                << "): shared-fault QPs use FaultFrictionLaw::LSW "
                << "dispatch in wave_operator (REVIEW R-016 fix; "
                << "shared-fault closed-form LSW solve via "
                << "FaultFaceFlux::EvaluateADER_LSW).\n";
   }
```

Confirm via `grep -n "MPI_Abort\|MFEM_ABORT" miniapps/seas/drivers/tpv205_driver.cpp`
that no R-016-tagged abort remains.

### Interfaces
- `InitializeFaultDOFs_TPV205` signature unchanged.
- `TPV205StationWriter::WriteStep` signature unchanged.
- `Tpv205SubStepIterator::StepOneQP_` signature unchanged.

### Edge Cases to Handle

- `disable_nucleation == true`: the post-init override
  `for (auto &d : dof_data) { d.tau2_0 = TPV205Params::tau_back; }`
  does NOT touch the LSW fields, so μ_s / μ_d / d_c stay correct
  (including the strength-barrier mask). No change needed.
- `num_fault_total == 0`: `dof_data` is empty; all loops are no-ops.
  No path reads `d.lsw_*`. ✓
- A barrier QP with `lsw_mu_s = TPV205Params::mu_s_barrier`:
  `LSWFrictionCoefficient_TPV205` and `SolveLSW_TPV205` short-circuit
  to V = 0. ParaView "state" channel renders 10000 (sentinel; users
  should clamp the colour ramp).

### Acceptance Criteria
- [ ] `make seas_tpv205_driver -j8` succeeds.
- [ ] `make seas_tpv102_driver -j8` and `make seas_tpv104_driver -j8`
  still succeed (no DOFData-offset breakage).
- [ ] `seas_test_tpv104_substep_iterator_parity`,
  `seas_test_tpv104_substep_one_shot_parity` pass byte-identical.
- [ ] `seas_test_tpv102_local`, `seas_test_tpv104_smoke` pass.
- [ ] TPV205 dry-run banner (`./seas_tpv205_driver --dry-run`) shows
  "Friction solver: LSW closed-form ..." (R-004 fix preserved).
- [ ] TPV205 single-rank 4-step smoke reproduces V_max ≈ 0.0779 m/s
  at step 0 (within 1% — the analytic LSW initial slip rate at the
  nucleation patch).
- [ ] TPV205 4-rank 4-step smoke reproduces V_max ≈ 0.0779 m/s on the
  rank that owns the nucleation patch. (Same as round-2; the iterator
  still owns interior fault QPs; the wire-up adds dispatch-routing
  for SHARED QPs but the smoke test doesn't yet stress them.)
- [ ] `grep -n "data\.a\|data\.psi\|data\.Dc" miniapps/seas/dynamic/tpv205_*.{hpp,cpp} miniapps/seas/drivers/tpv205_driver.cpp miniapps/seas/config/tpv205_params.hpp`
  returns no hits in non-comment lines (sanity grep).
- [ ] No `MPI_Abort` or `MFEM_ABORT` referencing R-016 remains in
  `drivers/tpv205_driver.cpp`.

### Dependencies
- Depends on: Phase 2.
- Required by: Phase 4.

---

## Phase 4: Tests — parity at O = 1 + np > 1 rupture-front MPI smoke

### Goal
After this phase the project has two new tests pinning the redesign:
(a) `seas_test_tpv205_evaluate_ader_lsw_parity` — proves
`FaultFaceFlux::EvaluateADER_LSW` agrees with
`Tpv205SubStepIterator::StepOneQP_` at O = 1 (pointwise Q == time-average Q),
and (b) `seas_test_tpv205_mpi_rupture_crossing` — np = 4 smoke that
forces the rupture front to cross at least one rank boundary and
asserts the rank-boundary V matches the analytic LSW value.

### Files to Create
- `tests/unit/test_tpv205_evaluate_ader_lsw_parity.cpp` — new.
- `tests/verification/test_tpv205_mpi_rupture_crossing.cpp` — new.

### Files to Modify
- `Makefile` — add the two new test targets (mirror the existing
  `seas_test_tpv104_substep_iterator_parity` and
  `seas_test_tpv104_smoke` Makefile stanzas).
- `tpv205/jobs/` (NEW directory or `jobs/tpv205/`) — sbatch script
  for Frontera np = 4 smoke (per project convention; not required
  for local CI).

### Detailed Requirements

**4.1 — `test_tpv205_evaluate_ader_lsw_parity.cpp`.**

The test constructs a single fault QP with TPV205 LSW parameters and
asserts that `FaultFaceFlux::EvaluateADER_LSW(data, I±, dt, ...)`
agrees with `Tpv205SubStepIterator::AdvanceWithSubStepStates(...)` at
O = 1 (single sub-step covering the whole macro-step, Q_pointwise[0] =
I/dt, weights = {1.0}).

Parameters:
- `data.lsw_mu_s = 0.677, data.lsw_mu_d = 0.525, data.lsw_d_c = 0.4`
  (interior rupture-area QP; not barrier).
- `data.sigma_n0 = 120e6, data.tau1_0 = 0, data.tau2_0 = 81.6e6`
  (nucleation patch).
- `data.tau1_nuc = 0, data.tau2_nuc = 0, data.sigma_n_nuc = 0`.
- `data.a = data.psi = data.Dc = 0` (rate-and-state slots zeroed).
- Impedances from `TPV205Params`.
- `dt = 1e-3` (1 ms — matches typical TPV205 macro-step).
- `Q_plus = Q_minus = 0` (initial fluctuation-Q rest state); thus
  `I_plus = I_minus = 0`.

Expected:
- After `EvaluateADER_LSW`: `data.slip_rate ≈ 0.0779 m/s` (analytic
  LSW initial slip rate at the nucleation patch).
- After running the iterator with the same `data` (re-initialised
  to the same state) and `dt`: `data.slip_rate ≈ 0.0779 m/s` and
  `Q_imp_plus / Q_imp_minus` match the `EvaluateADER_LSW` outputs to
  within `1e-12 * max(|component|, 1)`.

```cpp
// Sketch of the test body — full implementation per project's existing
// test patterns (see tests/unit/test_tpv104_substep_iterator_parity.cpp).
TEST(TPV205, R016_EvaluateADER_LSW_AgreesWithIteratorAtO1) {
   FaultFaceFlux flux(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   DOFData data_a, data_b;
   InitOneNucleationPatchQP(data_a);   // populates lsw_mu_s/mu_d/d_c, sigma_n0, tau2_0
   data_b = data_a;

   const real_t dt = 1e-3;
   real_t I_plus[NUM_STATE] = {0}, I_minus[NUM_STATE] = {0};
   real_t I_imp_a_p[NUM_STATE], I_imp_a_m[NUM_STATE];
   real_t I_imp_b_p[NUM_STATE], I_imp_b_m[NUM_STATE];

   // Path A: EvaluateADER_LSW (one-shot)
   flux.EvaluateADER_LSW(data_a, I_plus, I_minus, dt, I_imp_a_p, I_imp_a_m);

   // Path B: iterator at O=1
   Tpv205SubStepIterator iter(flux);
   iter.SetSubSteps({dt}, {1.0});
   std::vector<DOFData> dof_data{data_b};
   std::vector<Vector> coords(1, Vector(3));  // dummy coords
   coords[0] = 0.0;
   std::vector<std::vector<real_t>> Qp(1, std::vector<real_t>(NUM_STATE, 0.0));
   std::vector<std::vector<real_t>> Qm(1, std::vector<real_t>(NUM_STATE, 0.0));
   iter.AdvanceWithSubStepStates(dof_data, coords, Qp, Qm, dt, 0.0,
                                 I_imp_b_p, I_imp_b_m);

   // Assert agreement
   const real_t tol = 1e-12;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      EXPECT_NEAR(I_imp_a_p[c], I_imp_b_p[c],
                  tol * std::max(std::abs(I_imp_a_p[c]), real_t(1)));
      EXPECT_NEAR(I_imp_a_m[c], I_imp_b_m[c],
                  tol * std::max(std::abs(I_imp_a_m[c]), real_t(1)));
   }
   EXPECT_NEAR(data_a.slip_rate, data_b.slip_rate, 1e-12);
   // Analytic check: V = (81.6 - 0.677*120) MPa / eta_s = 0.36e6 / 4.624e6
   const real_t analytic = (81.6e6 - 0.677*120e6) / TPV205Params::eta_s;
   EXPECT_NEAR(data_a.slip_rate, analytic, 1e-3 * analytic);  // 0.1 % tol
}
```

**4.2 — `test_tpv205_mpi_rupture_crossing.cpp`.**

A `MFEM_USE_MPI`-gated test. Builds a small TPV205-style mesh (the
existing `tpv2053d_200m.msh` is too large for a unit test; either
generate a 2 km × 2 km × 2 km cube with a fault at y=0 in
`SetUp()` via `seas_generate_mesh`, or check in a tiny mesh under
`tpv205/mesh/tpv205_unittest_2km.msh`). Partitions across np = 2 ranks
so the fault is split at x = 0 (passing through the nucleation patch).
Runs 4 ADER-O2 steps (dt ≈ 1 ms). Expected V_max on each rank that
owns part of the nucleation patch matches the analytic 0.0779 m/s
to within 1 %.

Failure mode this test catches: if the wave-operator dispatch
regresses to the rate-and-state Brent on shared QPs, V_max on the
rank-boundary side drops to ≈ 0 and the rupture front stalls. The
test would catch that as `V_max < 0.5 * analytic`.

```cpp
TEST(TPV205, R016_RuptureCrossesMPIRankBoundary) {
   int rank, nprocs;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   ASSERT_EQ(nprocs, 2) << "Run with mpirun -np 2.";

   // Build / load a 2 km cube mesh with fault at y=0, partitioned at x=0.
   //   ... (use seas_generate_mesh helper or a checked-in tiny .msh)
   ParMesh pmesh = LoadTinyTPV205Mesh(MPI_COMM_WORLD);

   // Force partition split at x = 0: the nucleation patch (centred at
   // x = 0, depth = 7.5 km — but for the tiny test, scale to fit) sits
   // exactly on the rank boundary.

   WaveOperator<ParMesh> wave(pmesh, /*order*/1,
                              TPV205Params::lambda, TPV205Params::mu,
                              TPV205Params::rho, MakeBC());
   FaultFaceFlux fault_flux(TPV205Params::rho, TPV205Params::cp, TPV205Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);   // <-- the test of interest

   // Init dof_data, iterator, Q = 0, run 4 steps.
   //   ... (mirror tpv205_driver.cpp's init block)

   // Step 0 V_max check.
   const real_t analytic = (81.6e6 - 0.677*120e6) / TPV205Params::eta_s;
   real_t V_max_local = 0.0;
   for (auto &d : dof_data) {
      // Only check QPs inside the nucleation patch.
      if (InNucleationPatch(...)) {
         V_max_local = std::max(V_max_local, d.slip_rate);
      }
   }
   real_t V_max_global = 0.0;
   MPI_Allreduce(&V_max_local, &V_max_global, 1, MPI_DOUBLE, MPI_MAX,
                 MPI_COMM_WORLD);

   EXPECT_NEAR(V_max_global, analytic, 0.05 * analytic)  // 5% tol
       << "Rupture stalled at MPI rank boundary — likely R-016 regression.";

   // Per-rank V_max: every rank that owns part of the nucleation patch
   // must see V_max within 5% of analytic.  This catches the case where
   // ONE rank gets 0.0779 m/s and the other gets 0 because shared QPs
   // ran the wrong dispatch.
   if (RankOwnsPatchEdge(rank)) {
      EXPECT_GT(V_max_local, 0.5 * analytic)
          << "Rank " << rank << " saw V_max=" << V_max_local
          << " (analytic=" << analytic << ") — rupture stalled at boundary.";
   }
}
```

**4.3 — Makefile additions.**

Mirror the existing TPV104 substep-iterator-parity Makefile stanza.
Names:

```makefile
TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC = tests/unit/test_tpv205_evaluate_ader_lsw_parity.cpp
TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ = $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC:.cpp=.o)
TEST_TPV205_MPI_RUPTURE_CROSSING_SRC     = tests/verification/test_tpv205_mpi_rupture_crossing.cpp
TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ     = $(TEST_TPV205_MPI_RUPTURE_CROSSING_SRC:.cpp=.o)

seas_test_tpv205_evaluate_ader_lsw_parity: $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ) \
    $(TPV205_SHARED_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ) \
	    $(TPV205_SHARED_OBJS) $(MFEM_LIBS)

seas_test_tpv205_mpi_rupture_crossing: $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) \
    $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) \
    $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) \
	    $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) \
	    $(PML_LAYER_OBJ) $(MFEM_LIBS)

# Compile rules
$(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ): %.o: $(SRC)%.cpp \
    $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
$(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ): %.o: $(SRC)%.cpp \
    $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
```

Add both to the `test:` aggregate target if one exists for TPV205,
or add a `test-tpv205:` target mirroring `test-tpv104`.

### Interfaces
None changed — only test-target additions.

### Edge Cases to Handle

- Test 4.1 (`AdvanceWithSubStepStates`'s contract requires
  `fault_coords.size() == dof_data.size()`): supply a 1-element
  coord vector with `Vector(3)` zeros. Verified to pass the iterator's
  size check.
- Test 4.2 partition: ensure the partition explicitly puts the nucleation
  patch on the rank boundary; relying on ParMETIS may not always
  bisect at x = 0. Use the `--partition-file` mechanism with a
  hand-built sidecar partition: ranks 0/1 split at x = 0.
- Test 4.2 mesh: if the tiny test mesh isn't checked in, build it
  on the fly with `seas_generate_mesh` in the test's `SetUp()`.

### Acceptance Criteria
- [ ] `make seas_test_tpv205_evaluate_ader_lsw_parity -j8` builds.
- [ ] `./seas_test_tpv205_evaluate_ader_lsw_parity` passes (parity within 1e-12).
- [ ] `make seas_test_tpv205_mpi_rupture_crossing -j8` builds.
- [ ] `mpirun -np 2 ./seas_test_tpv205_mpi_rupture_crossing` passes
  on Frontera dev queue (or the Frontera sbatch helper if local
  MPI cannot oversubscribe per `feedback_no_local_oversubscribe`).
- [ ] After **deliberately reverting Phase 2's branch** (i.e.,
  forcing the wave op back to the pre-change `EvaluateADER` for shared
  QPs), `seas_test_tpv205_mpi_rupture_crossing` FAILS — confirming
  the test catches the R-016 regression.

### Dependencies
- Depends on: Phase 3 (driver wires `SetFaultFrictionLaw`).
- Required by: nothing.

---

## Testing Strategy

| Phase | Tests run                                                                                                | Expected outcome                                  |
|-------|----------------------------------------------------------------------------------------------------------|---------------------------------------------------|
| 1     | Existing TPV102 / TPV104 unit + driver smoke tests                                                       | Byte-identical pass                               |
| 2     | TPV102 / TPV104 substep parity tests; BP5 smoke                                                          | Byte-identical pass                               |
| 3     | TPV205 dry-run banner; TPV205 4-rank 4-step smoke                                                        | V_max ≈ 0.0779 m/s; banner says "lsw-closed-form" |
| 4     | New TPV205 parity test; new TPV205 MPI rupture-crossing test                                             | Both pass; revert-Phase-2 dry-run fails MPI test  |

Reference solutions:
- Analytic LSW initial slip rate at the SCEC TPV5 nucleation patch:
  `V = (τ_nuc − μ_s · σ_n) / η_s = (81.6 − 81.24) MPa / 4.624 MPa·s/m =
  0.0778 m/s` (matches the implementer's smoke result).
- DRDG3D reference traces under `tpv205/benchmark_data/DRDG3D_*` are
  the long-form benchmark; they're consumed by separate convergence
  studies, not this plan's tests.

## Risk Assessment

| Risk | How to detect | Mitigation                                                                                                              |
|------|---------------|-------------------------------------------------------------------------------------------------------------------------|
| **DOFData layout breakage**: appending fields shifts an offset some other code depends on. | TPV102/TPV104 parity tests fail with non-zero `Δ`. | Append at the end (after `tau1_corr / tau2_corr / sigma_n_corr` and before `#ifdef SEAS_DIAG_FAULT_FLUX bool diag_print`); never reorder existing fields.  Run `seas_test_fault_face_flux` immediately after Phase 1.                                  |
| **`EvaluateADER_LSW` slip double-counting** if the wave operator dispatch ends up calling it more than once per macro-step per QP (e.g., interior path branch + a future per-substep loop both fire). | Iterator parity test 4.1 fails with `slip > expected`. | Phase 2 dispatches per macro-step only.  Add an `MFEM_ASSERT` later if a second dispatch site is added.                  |
| **Frontera resource gate**: project rule `feedback_frontera_approval` requires explicit user approval before any sbatch.  | The plan mentions Frontera in test 4.2.            | Phase 4 only WRITES the sbatch script; running it is gated on user approval per the existing convention.                |
| **Bit-rot of `Tpv205SubStepIterator::Advance`** (REVIEW R-011) leaves a stale code path that reads `d.lsw_*` only after Phase 3.  | Compile error if `Advance` body still touches `d.a / d.psi / d.Dc`. | Phase 3 explicitly searches `dynamic/tpv205_substep_iterator.cpp` for the `Advance` body and migrates its reads — same translation as `StepOneQP_`.  If `Advance` is deleted instead (option in R-011 fix), this risk evaporates. |
| **Existing `_check.md` / `_fix.md` documents drifting from this plan.** | Round-4 review confusion. | Phase 3 updates `REVIEW.md` to mark R-016 / R-017 as FIXED with cross-reference to this plan's filename.                |
| **`SEAS_DIAG_FAULT_FLUX` build alters `DOFData` size**: the diag flag adds a `bool diag_print` field.  Appending the LSW fields changes total size by 24 bytes; with the diag flag, total grows by the same 24 bytes plus the existing `bool` (no realignment issue on x86_64).  | Compilation passes; `sizeof(DOFData)` is consistent across translation units (single compile flag enforced by the build system). | No action; just be aware that the size check must be done with the same `SEAS_DIAG_FAULT_FLUX` setting.                |

## Roll-back triggers

| Phase | Trigger                                                                                                    | Roll-back action                                                                              |
|-------|------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| 1     | Any TPV102 / TPV104 unit test fails with non-zero Δ vs pre-change baseline.                                | `git revert` Phase-1 commit; investigate which existing reader iterates DOFData by offset.    |
| 2     | TPV104 parity test fails or BP5 smoke fails.                                                               | `git revert` Phase-2 commit; the dispatch enum default `RateAndState` should make this impossible. |
| 3     | TPV205 4-rank 4-step smoke V_max no longer matches 0.0779 m/s within 1 %.                                  | `git revert` Phase-3 commit; check that `SetFaultFrictionLaw(LSW)` is wired BEFORE the time loop. |
| 4     | New tests fail to build due to test-fixture mesh problems.                                                 | Skip MPI rupture-crossing test as `DISABLED_` until the tiny mesh is checked in / scripted.   |

## Post-merge follow-ups (out of scope for this plan)

- TPV205 RK4 dispatch path (the non-ADER `fault_flux_->Evaluate(...)` sites
  at wave_operator.inl ~L2419 and ~L3058) is currently unwired for LSW.
  TPV205 driver hard-locks to ADER, so this is dormant. If a future
  TPV205 RK4 driver lands, extend the dispatch the same way.
- Delete `Tpv205SubStepIterator::Advance` once Phase 4 confirms no
  caller touches it (REVIEW R-011 close-out).
- Consider promoting the LSW-native fields to a separate `LSWParams`
  POD and pointing `DOFData` at it via a discriminated union (would
  bound DOFData growth at one pointer, but breaks the no-touch-extreme
  list further; defer.)
