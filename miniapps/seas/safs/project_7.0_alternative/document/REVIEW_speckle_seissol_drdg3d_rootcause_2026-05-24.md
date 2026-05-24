# Code Review: SAFS speckle runaway — DEBUG doc audit + SeisSol/drdg3d comparison — 2026-05-24

## Review Scope
- **Reviewed document:** `safs/project_7.0_alternative/document/DEBUG_speckle_normal_velocity_jump_2026-05-23.md`
- **SAFS source reviewed:**
  - `dynamic/tpv205_substep_iterator.cpp` (`StepOneQP_`, `AdvanceWithSubStepStates`)
  - `dynamic/fault_face_flux.cpp` (`ComputeTrialTraction`, `EvaluateADER_LSW`, `BuildImposedState`)
  - `dynamic/wave_operator.inl` (`ComputeADERSubStepStates`, `EvaluateBulkAtFaultQPsCanonical`, `ComputeADERFaceFluxRHS`, `ComputeADERSharedFaceFluxRHS`, ctor ghost-GF construction)
  - `dynamic/wave_operator.hpp` (ghost-GF members)
  - `drivers/spatial_dyn_driver.cpp` (`AdvanceADERWithSubStep_Spatial`)
- **Reference implementations compared (read by sub-agents, file:line verified):**
  - SeisSol `/Users/chunhuizhao/projects/SeisSol` — `FrictionLaws/`, `Kernels/DynamicRupture.cpp`, `Solver/TimeStepping/TimeCluster.cpp`
  - drdg3d `/Users/chunhuizhao/projects/drdg3d` — `mod_wave.F90`, `mod_eqns.F90`, `seis3d.F90`, `mod_exchange.F90`
- **Domain context:** project `CLAUDE.md`, prior `REVIEW_DEBUG_speckle_normal_velocity_jump`, auto-memory `project_safs_vn_leak_not_frame.md`.

> **Note on file placement.** I did **not** overwrite the repo-root `REVIEW.md` — it currently holds the in-flight diagnostics-code review (`R-001…R-005`) paired with `REVIEW_fix.md`. This investigation is a distinct root-cause study and is filed alongside the existing speckle `REVIEW_*.md` chain in `…/document/`. Move it into `REVIEW.md` if you want the `/code-fix` agent to consume it directly.

---

## The reference invariant (what both SeisSol and drdg3d guarantee)

Both references enforce **one slip rate, one flux**:

> The slip rate computed by the friction solve is used simultaneously for (a) the radiation-damping traction that is fluxed back to the bulk, **and** (b) the cumulative slip that weakens μ. They are the *same number* at the *same point and time*. There is no second, independently-evaluated solve.

- **SeisSol:** one `slipRateMagnitude = (absoluteTraction − strength)·invEtaS` (`FrictionLaws/CpuImpl/LinearSlipWeakening.h:77-86`) drives the corrected flux traction (`:89-96`), the directional slip (`:99-102`), **and** the slip-weakening state `accumulatedSlipMagnitude` (`:185-186`). The corrected per-time-node tractions are time-integrated with `timeWeights` into a single `imposedState` (`FrictionLaws/FrictionSolverCommon.h:329-363`) applied once by `nodalFlux` (`Kernels/LinearCK/Neighbor.cpp:93-98`). MPI-boundary ("Copy") fault faces use the **identical** friction/flux code as interior ("Interior") faces — only the source of the neighbor's predictor DOFs (ghost storage) differs (`TimeCluster.cpp:243-244`, `CellLocalMatrices.cpp:529-536`). **No partition-dependent flux path.**
- **drdg3d:** one `vv1/vv2` (`mod_wave.F90:963-964`) both modifies the tangential flux `fstar(2,3,8,9) -= eta*vv` (`:921-924`) and accumulates `slip1/slip2` in the RK loop (`seis3d.F90:275-352`). MPI fault faces use the generic ghost path with no special-casing (`mod_wave.F90:445-467`, `mod_exchange.F90:47-115`). Additionally drdg3d **welds the normal direction**: only tangential components of the flux are touched; `fstar(4)` (normal velocity) and `fstar(1)` (normal traction) keep their continuous Godunov values, so a normal-velocity jump (opening) **cannot arise by construction**.

The SAFS dynamic Godunov trial traction itself is correct and matches both references term-for-term (`fault_face_flux.cpp:62-64` ≡ SeisSol `FrictionSolverCommon.h:182-192` ≡ drdg3d `mod_wave.F90:675-701`). The deviation is **not** in the flux formula; it is in **how many times, and on which state, the solve is run.**

---

## Findings

### [R-001] CRITICAL [wave_operator.inl:4856-4895 / tpv205_substep_iterator.cpp:125-126] — Shared faces run TWO decoupled solves; slip-weakening and applied flux use different slip rates (single-slip-rate invariant violated)

**Category:** DEVIATION / BUG

**Description:**
On a **shared (MPI-partition-boundary) fault face**, SAFS runs the friction solve **twice on two different bulk states**, and the slip that weakens μ comes from the solve whose flux is **thrown away**:

1. **Iterator (predictor) solve** — `AdvanceADERWithSubStep_Spatial` (`spatial_dyn_driver.cpp:439`) calls `iterator.AdvanceWithSubStepStates(dof_data, …)` over **all** fault QPs incl. shared. For each sub-step it solves LSW on the **pointwise CK predictor** `Q_pointwise_±[o]` and accumulates `d.slip1 += s.V1·dt_sub` (`tpv205_substep_iterator.cpp:125-126`) — the **sole** slip-accumulation site in the code (confirmed `fault_face_flux.cpp:815-824`: "slip evolution is owned exclusively by the iterator"). It also builds an imposed state into `I_imp_*`.
2. **Macro (corrector) solve** — `AdvanceADER` → `ComputeADERSharedFaceFluxRHS` **discards** the iterator's `I_imp_*` for shared QPs (`wave_operator.inl:4856-4858, 4893-4895`) and re-solves via `EvaluateADER_LSW` on the **time-integral** `I_{plus,minus}_local/dt` (`:4878-4884`). This solve **reads** `data.slip1/slip2` for μ (`fault_face_flux.cpp:794`) but does **not** write it; its corrected traction is the flux actually applied to the bulk.

So on shared faces: `μ ← δ(predictor V)` while `flux ← V(time-integral)`. The two slip rates are unrelated. SeisSol/drdg3d make this impossible — the rate that weakens μ *is* the rate that is fluxed (and damped). This is the structural enabler of the runaway: the iterator's predictor V can diverge arbitrarily from the macro V, accumulate unbounded `max_slip` (the symptom in §4 of the DEBUG doc), and drive μ→μ_d off a slip history that the applied flux never reflects. Interior faces do **not** have this defect — `ComputeADERFaceFluxRHS:3871-3884` applies the iterator's own `I_imp_*` as the flux, so there the invariant holds (one solve, like SeisSol).

This is the same R-1601 fallback the DEBUG doc references, but the doc treats it as a *stability workaround*; the reference comparison shows it is the **mechanism**, not a side issue.

**Trigger:** Any shared fault face under dynamic slip at `np > 1`. (At `np = 1` `fault_shared_faces_` is empty and the defect cannot fire — consistent with the runaway being observed only at `np > 1`.)

**Actual behavior:** μ on shared QPs is weakened by the predictor iterator's slip; the applied bulk flux is a separate macro-dt solve. The predictor's `max_slip`/V runs away (→4.6e6 m, →1.1e7 m/s) invisibly to the bounded macro output (V_max≈5–7 m/s, σ_n≈47 MPa).

**Expected behavior (reference):** one slip rate per shared QP drives both the slip-weakening state and the flux applied to both bulk sides, exactly as on interior faces and exactly as SeisSol's Copy-layer faces.

**Suggested fix:** eliminate the dual solve on shared faces. Preferred (matches SeisSol interior parity and the SAFS interior path): make `ComputeADERSharedFaceFluxRHS` consume the iterator's per-sub-step `substep_I_imp_*` for shared QPs, the same way `ComputeADERFaceFluxRHS:3871-3884` does for interior QPs, once the frame-convention reconciliation that R-1601 deferred is done (the predictor and macro **already use the identical `sign_flipped` frame and `elem1_on_plus` routing** — verified below in R-003 — so the original R-1601 "frame mismatch" rationale at `:4838-4854` no longer holds and should be re-tested):
```diff
- // SHARED FALLBACK: always run inline ADER closure, discard substep_I_imp_*
- fault_flux_->EvaluateADER_LSW(fdata, I_plus_local, I_minus_local, dt, I_imp_plus, I_imp_minus);
- (void)substep_I_imp_plus_flat_; (void)substep_I_imp_minus_flat_;
+ // Consume the iterator's per-sub-step imposed state on shared QPs too
+ // (single-slip-rate invariant; mirrors ComputeADERFaceFluxRHS:3871-3884).
+ if (substep_I_imp_plus_flat_ && substep_I_imp_minus_flat_) {
+    const real_t *sp = substep_I_imp_plus_flat_  + dof_idx*NUM_STATE;
+    const real_t *sm = substep_I_imp_minus_flat_ + dof_idx*NUM_STATE;
+    for (int c=0;c<NUM_STATE;c++){ I_imp_plus[c]=sp[c]; I_imp_minus[c]=sm[c]; }
+ } else {
+    fault_flux_->EvaluateADER_LSW(fdata, I_plus_local, I_minus_local, dt, I_imp_plus, I_imp_minus);
+ }
```
If R-1601 must stay for now, the minimum correctness patch is to **stop weakening μ off the discarded solve**: on shared QPs, accumulate `d.slip` from the **macro** V (the one that is actually fluxed), not the predictor V — i.e. move the shared-QP slip accumulation into `EvaluateADER_LSW` and skip it in the iterator for `i >= diag_num_local_fault_qps_`. That restores the invariant `μ ← δ(applied V)` even with the dual dispatch.

**Test case (C++/MFEM unit test, the stack here — adapt `tests/unit/`):**
```cpp
// test_shared_slip_rate_invariant.cpp
// On a 2-rank fault with one shared face, after one macro step the slip
// increment recorded for a shared QP must equal V_applied * dt, where
// V_applied is the slip rate of the flux deposited into the bulk RHS —
// NOT the predictor iterator's V. Today they differ on shared faces.
TEST(SharedFault, SlipWeakeningRateEqualsFluxedRate) {
   auto [dof, V_applied, dt] = RunOneMacroStepSharedQP(/*np=*/2);
   double delta_slip = dof.slip1 - dof.slip1_before;          // accumulated
   ASSERT_NEAR(delta_slip, V_applied * dt, 1e-12 * std::abs(V_applied*dt));
}
```

---

### [R-002] HIGH [DEBUG doc §"Root-cause investigation" H-A/H-B] — The doc's own `[MACRO]` data already EXCLUDES H-A (bulk pumping); H-A/H-B are not "indiscriminable"

**Category:** ASSUMPTION (analysis error in the reviewed document)

**Description:**
The doc states H-A (bulk pumping) and H-B (predictor/ghost divergence) are "**not yet discriminable from this log**." They are discriminable, and the log already **rules out H-A**. The doc reports two facts that are jointly inconsistent with bulk pumping:

- **Fact 1:** within one macro step the predictor `sn_vjump` is **flat** (o=0→o=1: −4.170e7→−4.219e7, +1%). ⇒ the predictor jump is dominated by `D(0)=Q` (the base state), not by a large-τ CK term.
- **Fact 2:** the macro time-integral sides are **both bounded all run**: `vn_plus∈[−0.016,0.495]`, `vn_minus∈[−0.123,0.231]` m/s.

`I/dt` is, by construction, the macro-step **time-average** of the same predictor field, and both `ComputeADERSubStepStates` (pointwise) and `ComputeADERTimeIntegrated` (integral) are the **identical element-local CK recursion** of the **same** input `Q` (verified: `wave_operator.inl:1355-1404` vs the integrated variant; the factorial phases differ but are algebraically equal, as the `:1381-1394` comment proves). Therefore:

> If the predictor is flat at value `P` across the step (Fact 1), then `I/dt ≈ P`. H-A claims the *shared bulk field* `Q` grows secularly; but both paths read that same `Q`, so if `Q` grew, **Fact 2 (`I/dt`) would grow too**. It does not. ⇒ the shared bulk field is **not** pumping; the divergence lives in the **predictor/ghost evaluation path itself** (H-B).

Moreover, H-A is mechanically impossible on shared faces for an independent reason: the iterator's imposed flux on shared faces is **discarded** (R-001), so the iterator cannot pump the shared bulk; only the (bounded) macro flux is applied there.

**Consequence for the test plan:** DEBUG doc tests #5 (tighten dt) and the §7 primary fix ("make the frame geometrically correct") are **dead ends** — both target mechanisms already excluded. The remaining work is to localize *within H-B* (self vs nbr), which the doc's own "decisive next measurement" (print `Q_self[VX]` / `Q_nbr[VX]` separately) does — see R-004 for the predicted result.

**Trigger:** N/A (analytical).

**Suggested fix:** in the DEBUG doc, replace "two surviving sub-hypotheses, not yet discriminable" with "H-A excluded by Fact 1 + Fact 2; locus is the predictor/ghost evaluation (H-B); remaining question is self-side vs neighbor-side." Drop test #5 and the §7 frame-fix as the primary direction.

**Test case:** the discriminator is already shipped — the `[FRAME]` block at `wave_operator.inl:2258-2278` has `Q_self` and `Q_nbr` in scope; add two fields `self_vx=%+.6e nbr_vx=%+.6e` (`Q_self[VX]`, `Q_nbr[VX]`) and rerun. Assertion of the expected result is R-004.

---

### [R-003] MODERATE [DEBUG doc §TL;DR + §7] — Document is internally inconsistent: TL;DR/§7 still prescribe a frame fix that the doc itself refuted; the frame is provably not the cause

**Category:** DEVIATION (doc inconsistency that misdirects the fix)

**Description:**
The doc's TL;DR (lines 16-20) and §7 "Fix direction" (lines 347-360) still name the **fault frame** ("make the canonical fault frame on the curvilinear shared face geometrically correct") as the *primary* fix, while the doc's own 2026-05-23 UPDATE (lines 108-121) and the prior `REVIEW_DEBUG_speckle…` already refuted the frame mechanism. The reference comparison closes this definitively:

- The predictor path (`EvaluateBulkAtFaultQPsCanonical`, shared branch) and the macro path (`ComputeADERSharedFaceFluxRHS`) build the canonical frame **identically**: both negate `{normal,tangent1,tangent2}` on `qpd.sign_flipped` (`wave_operator.inl:2199-2204` vs `:4805-4810`) and route ± identically on `elem1_on_plus` (`:2302-2305` vs `:4831-4834`). **A shared frame cannot produce a predictor-vs-macro divergence** — both paths use the very same `can_n/can_t1/can_t2` and the same `Tinv_can`.
- The frame is orthonormal by construction (`fault_basis.hpp` cross products; `test_fault_basis_qp_orthonormality` passes), so `V·can_n ≡ 0` and friction slip cannot inject a normal jump — already established in auto-memory `project_safs_vn_leak_not_frame.md`.

Leaving the frame fix in §7 will send the `/code-fix` agent to re-derive `ComputeOrientedFrame` — wasted effort on a refuted cause.

**Suggested fix:** rewrite §7 to point at R-001 (eliminate the shared-face dual solve / restore the single-slip-rate invariant) as the primary fix, and R-004 (the predictor ghost-exchange divergence) as the trigger to localize. Remove "frame accuracy" from the fix list.

**Test case:** `test_frame_diag_projection_identity` (already in `tests/unit/`) demonstrates the frame identity holds; no new test — this is a doc correction.

---

### [R-004] MODERATE [POSSIBLE] [wave_operator.inl:2025-2059] — Prime concrete suspect: the R-1601 `vdim=NUM_STATE` byNODES batched ghost exchange is the ONLY material code difference between the bounded macro path and the runaway predictor path

**Category:** BUG (POSSIBLE — locus high-confidence, specific defect not yet proven)

**Description:**
After R-001/R-002/R-003, the runaway is localized to the **neighbor (ghost) side of the per-sub-step predictor evaluation**. Comparing the two shared-face ghost exchanges line-by-line, they are identical *except* in the exchange mechanism:

| | Macro `ComputeADERSharedFaceFluxRHS` | Predictor `EvaluateBulkAtFaultQPsCanonical` |
|---|---|---|
| Ghost GF | scalar `ghost_gf_` (vdim=1) | `ghost_gf_full_state_` (vdim=NUM_STATE, **byNODES**) |
| Exchange | `NUM_STATE` per-component `ExchangeFaceNbrData` (`:4639-4650`) | **one** batched `ExchangeFaceNbrData` (`:2037`) |
| Unpack | `nbr_data[c] = src` per component (well-tested) | assumes `FaceNbrData()` layout `src[c*n_face_nbr_dofs + i]` (`:2050-2058`) |
| Read | `nbr_data[c][nbr_idx*ndof_per_el_+i]` (`:4772`) | `nbr_data[c][nbr_idx*ndof_per_el_+i]` (`:2170`) |

The macro per-component scalar exchange is the **pre-R-1601** mechanism and is battle-tested (it is what the predictor path *also* used before R-1601 — see the `wave_operator.hpp:997-1004` comment). R-1601 replaced **only the predictor** exchange with the batched byNODES path to save `8·O` collectives; the macro kept the scalar path "because of the deep-copy correctness guard" (`:1007-1009`). The runaway is on the predictor path; the bounded path kept the old exchange.

Why this is the suspect, quantitatively: at onset (Fact 1, flat-in-τ), both `I/dt` (macro) and `Q_pointwise` (predictor) reduce to `≈ D(0)_nbr` = the neighbor's base `Q`. The doc measures them differing by ~10× at t=0.477 (predictor `sn_vjump≈−4.17e7` Pa ⇒ `dv_n≈−5.2` m/s; macro `vn` bounded ≤~0.5 m/s) and by 2.6e8 at end. Two CK expansions of the **same** base `Q` cannot differ by 10× — so the two paths are **reading different neighbor DOF values**, i.e. an exchange/layout inconsistency on the byNODES path. The byNODES `FaceNbrData` layout assumption (`src[c*n_face_nbr_dofs+i]` with the face-neighbor element ordering matching the scalar space's `nbr_idx*ndof_per_el_`) is plausible by inspection but is a separate, less-exercised MFEM code path; a mismatch here produces exactly the observed signature (**self side agrees, neighbor side diverges, present from onset, grows as the decoupled μ-feedback saturates**).

I cannot prove the byNODES layout is wrong from static inspection alone (it *looks* correct), so this is flagged POSSIBLE — but it is the single highest-probability concrete defect and is cheap to test decisively.

**Trigger:** `np > 1` shared fault face, every sub-step.

**Suggested fix (decisive test first, then fix):**
1. **A/B the exchange (fastest):** temporarily replace the predictor's batched exchange with the macro's per-component scalar exchange (loop `c`, pack `ghost_gf_`, `ExchangeFaceNbrData`, copy) — i.e. revert to pre-R-1601 mechanism — and rerun the Dc2 case. If the runaway vanishes, R-1601's byNODES batching is the bug; then fix the layout/ordering rather than reverting (to keep the `8·O` collective saving).
2. If confirmed, the fix is to derive the neighbor index through the **vector** space's face-neighbor vdof map (`pfes_full_state_->GetFaceNbrElementVDofs` / `DofToVDof` for byNODES) instead of assuming `c*n_face_nbr_dofs + nbr_idx*ndof_per_el_ + i`, OR keep byNODES only if a unit test (below) proves bit-identity with the scalar exchange.

**Test case (deterministic, CI-able — no Frontera run needed):**
```cpp
// test_ghost_exchange_bynodes_vs_scalar.cpp  (2 MPI ranks, 1 shared face)
// Fill Q[c*ndof_total + dof] = encode(rank, c, dof). Run BOTH exchanges on the
// SAME Q and assert the neighbor values read at every (c, nbr_idx, i) match.
for (int c=0;c<NUM_STATE;c++)
  for (i over nbr element dofs)
    ASSERT_EQ(nbr_data_bynodes[c][nbr_idx*ndof_per_el_+i],
              nbr_data_scalar [c][nbr_idx*ndof_per_el_+i]);  // must be bit-identical
```

---

### [R-005] LOW [DEBUG doc §5 R-005d / §3] — Free-slide-under-tension matches both references (doc correct), but the reason the references don't blow up is the structural coupling SAFS breaks on shared faces — sharpen the fix direction accordingly

**Category:** QUALITY (correct conclusion, incomplete reasoning that affects the fix)

**Description:**
The doc's refutation of "change the friction law" (§5 R-005d) is **correct**: SeisSol (`LinearSlipWeakening.h:151-160`, `strength=−cohesion−μ·min(σ_n,0)`) and drdg3d (`mod_wave.F90:837`, `Tau_str=μ·max(0,−Tau_n)+C0`) both free-slide under tension, so changing SAFS's LSW would deviate *from* the references. But the doc does not state *why* the references stay bounded under the same free-slide law, which matters for the fix:
- In both references the **same** slip rate that frees under tension is immediately returned as the radiation-damping flux `−η·V` to the bulk (SeisSol `:89-96`; drdg3d `:921-924`), which **damps** the velocity jump that caused the tension. SAFS **discards** that damping flux on shared faces (R-001) — the freeing solve has no feedback to the bulk it destabilizes.
- drdg3d additionally **welds the normal** (only tangential flux modified), so the opening that triggers free-slide cannot occur at all.

**Suggested fix:** none to code here; this confirms the R-001 fix direction (restore the slip-rate↔flux coupling on shared faces) and offers a defense-in-depth option borrowed from drdg3d: weld the normal-velocity jump in the shared imposed state (`BuildImposedState` already welds it for the iterator path — `fault_face_flux.cpp:278,286` — so a normal-welded shared corrector is consistent with the rest of the code). Note this is *defense in depth*, not a substitute for R-001: the DEBUG doc's cap run 7747835 already showed a no-opening cap alone is necessary-but-not-sufficient.

**Test case:** covered by R-001's invariant test.

---

## The bug, identified

**Root cause (definitive):** SAFS violates the SeisSol/drdg3d **single-slip-rate invariant on shared (MPI-boundary) fault faces** (R-001). The slip that weakens μ is accumulated by the per-sub-step **predictor** iterator, while the flux applied to the bulk is a **separate** macro-dt solve on the time-integral, and the iterator's own (normal-welded, radiation-damping) imposed state is **discarded** there. Interior faces are correct; only shared faces are decoupled (the R-1601 fallback).

**Trigger (high-confidence locus, specific defect flagged POSSIBLE):** the predictor's neighbor-side normal-velocity-jump evaluation diverges from the macro's from the first sub-step (R-002 proves the divergence is in the predictor/ghost evaluation, not the bulk field or the frame). The single material code difference on that path is the **R-1601 `vdim=NUM_STATE` byNODES batched ghost exchange** (R-004), which replaced the battle-tested per-component scalar exchange that the bounded macro path still uses.

**Why the symptom looks the way it does:** the predictor's runaway V integrates into `max_slip` (→4.6e6 m) and pushes μ→μ_d; once μ saturates at μ_d (δ>d_c), the macro flux stays bounded (hence bounded V_max≈5 m/s and compressive σ_n in the output), while the predictor — decoupled and undamped on shared faces — keeps integrating its corrupted neighbor jump. That is exactly the DEBUG doc §4 "iterator-vs-output decoupling," now explained mechanistically rather than just observed.

**Recommended fix order:**
1. **R-004 A/B test** (revert predictor exchange to per-component scalar for one run). Cheap, decisive on the trigger.
2. **R-001** (eliminate the shared-face dual solve; consume the iterator's `I_imp_*` on shared QPs, or at minimum accumulate shared-QP slip from the *applied* macro V). This restores the reference invariant and is the durable fix.
3. **R-003** doc correction (drop the frame fix from §7).
4. **R-005** optional defense-in-depth (weld normal jump in the shared corrector).

---

## Summary
- Critical issues: 1 (R-001)
- High/Moderate issues: 3 (R-002 high, R-003 moderate, R-004 moderate/possible)
- Low issues: 1 (R-005)
- Document-hypothesis verdict: **PARTIALLY CORRECT** — the doc correctly localized to "per-sub-step predictor/ghost path on shared faces (R-1303/R-1601)" and correctly refuted the frame/tension/friction-law fixes, but (a) still prescribes the refuted frame fix in §7 (R-003), (b) wrongly calls H-A/H-B indiscriminable when its own data excludes H-A (R-002), and (c) does not identify the single-slip-rate-invariant violation (R-001) that is the actual root design flaw the reference comparison exposes.
- Verdict: **FAIL — must fix R-001 before the shared-fault dynamic path is physical; run the R-004 A/B test to confirm the trigger.**

## Appendix A — R-001 in equations + roadmap

### A.0 The per-QP friction primitive `F`

Every solve (iterator or macro) calls the **same** primitive. Given the two bulk traces `Q⁺,Q⁻` at a fault QP (in the canonical frame: normal velocity `v_n`, tangential velocity `v_t`, normal stress `σ_n`, shear `τ`) and the current accumulated slip `δ`, it computes (`fault_face_flux.cpp:ComputeTrialTraction` + `SolveLSW_TPV205`):

```
(1) Godunov trial traction (radiation form):
      σ_n,trial = η_p (v_n⁻ − v_n⁺)  +  η_p (σ_n⁺/Z_p⁺ + σ_n⁻/Z_p⁻)     [the "opening" channel]
      τ_trial   = η_s (v_t⁻ − v_t⁺)  +  η_s (τ⁺/Z_s⁺  + τ⁻/Z_s⁻)
(2) σ_n,tot = σ_n0 + σ_n,trial ;   τ_tot = τ_0 + τ_trial
(3) μ(δ)    = μ_s − (μ_s − μ_d)·min(δ/d_c, 1)                            [LSW weakening]
(4) τ_str   = μ(δ)·max(σ_n,tot, 0)                                       [free-slide when σ_n,tot ≤ 0]
(5) V       = max(0, (|τ_tot| − τ_str)/η_s)                              [slip rate]
(6) τ_corr  = τ_tot − η_s·V                                             [APPLIED shear traction]
(7) Q_imp±  : shear ← τ_corr, normal welded so [[v_n]]_imp = 0           [BuildImposedState]
```

Write this compactly as **`(V, Q_imp) = F(Q⁺, Q⁻; δ)`**.

**The brake lives in (6).** The traction actually pushed onto the bulk is `τ_corr = τ_tot − η_s·V`. The `−η_s·V` term is radiation damping: it is *anti-parallel* to the slip, so it removes exactly the momentum that the velocity jump in (1) injected. Step (7) does the same for the normal channel by welding `[[v_n]]=0`. **A solve is only dissipative if its `Q_imp` is actually applied to the bulk.**

### A.1 Reference (SeisSol/drdg3d) and SAFS **interior** — one solve, stable

Loop over time-quadrature nodes `o` (SeisSol) / RK stages (drdg3d) / sub-steps (SAFS interior):

```
for o:  (V_o, Q_imp,o) = F( Q̃⁺(τ_o), Q̃⁻(τ_o); δ )      // Q̃ = pointwise predictor at τ_o
        δ      ← δ + V_o·Δt_o                            // (a) weakens μ      ── uses V_o
        Î_imp  ← Î_imp + w_o·Q_imp,o                      // (b) flux           ── uses SAME V_o
apply Î_imp to BOTH bulk sides.
```

`V_o` appears in **both** (a) and (b). So the slip that weakens μ is exactly the slip whose damping flux `−η_s·V_o` is fed back to the bulk. The map "jump → V → damping flux → smaller jump" is a **negative feedback** (gain < 1 by construction). SAFS interior is literally this: `tpv205_substep_iterator.cpp:125-126` does (a), and `ComputeADERFaceFluxRHS:3871-3884` applies the iterator's own `Î_imp` as the flux. ✔ stable.

### A.2 SAFS **shared** — two solves, the invariant breaks

On a shared face the code runs `F` **twice on different states**, and crosses the wires:

```
SOLVE A  (iterator, on the POINTWISE predictor)   tpv205_substep_iterator.cpp / StepOneQP_
   for o: (V_oᴬ, Q_imp,oᴬ) = F( Q̃⁺(τ_o), Q̃⁻(τ_o); δ )
          δ ← δ + V_oᴬ·Δt_o          ←── weakens μ           (slip uses Vᴬ)
          Î_impᴬ accumulated ……………… DISCARDED   wave_operator.inl:4856-4895

SOLVE B  (macro, on the TIME-INTEGRAL)            ComputeADERSharedFaceFluxRHS
   (Vᴮ, Q_impᴮ) = F( Ī⁺/dt, Ī⁻/dt; δ )            reads δ for μ, does NOT write δ  (fault_face_flux.cpp:815-824)
   apply Q_impᴮ to bulk           ←── flux uses Vᴮ
```

Collapsing the two:

```
        μ  ←  δ  ←  ∫ Vᴬ dt          (PREDICTOR slip rate)
      flux ←  Q_impᴮ  ←  Vᴮ          (TIME-INTEGRAL slip rate)
```

`Vᴬ ≠ Vᴮ` in general, and **the damping `−η_s·Vᴬ` from Solve A is thrown away** (its `Î_impᴬ` is discarded). So the predictor's opening jump `J ≡ ṽ_n⁻ − ṽ_n⁺` has **no brake**: nothing feeds `−η_s·Vᴬ` back to oppose it. Meanwhile Solve B's own brake `−η_s·Vᴮ` keeps the *time-integral* jump `Ī/dt` bounded — which is why `[MACRO] vn`, `V_max`, and output `σ_n` all look healthy while `max_slip = ∫Vᴬdt` explodes.

### A.3 Why it actually runs away — the conjunction with R-004

Crucial subtlety: **the decoupling alone is not sufficient.** If the predictor and the time-integral agreed (`Q̃ ≈ Ī/dt`, which holds to `O(dt²)` for the *same* base state — both are the identical element-local CK expansion, `wave_operator.inl:1355-1404`), then `F` is called on nearly equal arguments and `Vᴬ ≈ Vᴮ`; δ would track the realized macro slip and μ would stay consistent — **no runaway**. The runaway needs a **seed**: `J^pred ≠ J^macro`. The doc measured them differing ~10× already at onset (t=0.477) and `2.6e8` at the end — a discrepancy that two CK expansions of the *same* `Q` cannot produce, hence an **exchange/evaluation divergence on the neighbor side** (R-004, prime suspect = the R-1601 byNODES batched exchange).

So the runaway = **spark × no-brake**:

```
   R-004 (spark)         :  J^pred diverges from J^macro on the neighbor side
   R-001 (missing brake) :  the predictor solve's damping −η_s·Vᴬ is discarded; δ integrates Vᴬ undamped
   ─────────────────────────────────────────────────────────────────────────────
   ⇒ closed loop:   J^pred ↑  →  σ_n,trial = η_p·J^pred collapses (eq.1)
                              →  τ_str → 0 (free slide, eq.4-5)  →  Vᴬ ↑
                              →  δ ↑  →  μ → μ_d (eq.3)
                              →  macro releases full stress (eq.6 with small μ)
                              →  (no damping on the predictor channel) J^pred ↑↑   [doc gain ≈1.25–1.32 / macro step]
```

The loop is **sign-indefinite** (doc §3, cap run 7747835): for an *opening* QP (`J<0`) σ_n goes tensile and frees the slide; for a *closing* QP (`J>0`) σ_n→+TPa but `τ_tot` grows even faster, so `V=(|τ_tot|−τ_str)/η_s` still blows up. Both blow up precisely because the channel that drives δ has lost its `−η_s·V` brake.

**Worked numbers at the seed** (doc §2.2, qp 480, `η_p = ½ρc_p = 8.0e6 Pa·s/m`, `σ_n0 = 49 MPa`):
- `J^pred = −5.2 m/s` ⇒ `σ_n,trial = η_p·J^pred ≈ 8.0e6·(−5.2) = −4.16e7 Pa = −41.6 MPa` (matches `sn_vjump = −4.17e7`); bulk term `sn_sterm ≈ 0`.
- `σ_n,tot = 49 − 41.7 = +6.0 MPa` — *just* above collapse. Tensile onset threshold: `J_crit = −σ_n0/η_p = −49e6/8.0e6 ≈ −6.1 m/s`. A small further increase in `|J^pred|` drives `σ_n,tot < 0` ⇒ `τ_str = 0` ⇒ free slide ⇒ `Vᴬ = |τ_tot|/η_s` huge ⇒ δ runaway. On interior faces the same `−η_s·Vᴬ` would have pushed `J` back below threshold; on shared faces it is discarded.

### A.4 Roadmap (fix order)

| Step | Action | Restores | Durable? |
|---|---|---|---|
| **0** | **R-004 A/B test** — swap the predictor's R-1601 byNODES exchange (`wave_operator.inl:2025-2059`) for the macro's per-component scalar exchange (`:4639-4650`) for one Dc2 run. Runaway vanishes ⇒ exchange is the spark. | confirms trigger | — (diagnostic) |
| **1** | **R-001 primary** — on shared QPs, apply the iterator's `Î_imp` as the flux, like interior (`:3871-3884`); delete the `EvaluateADER_LSW`-on-`I/dt` shared branch. Re-test the R-1601 frame rationale (now moot: predictor & macro already share `sign_flipped`+`elem1_on_plus`). | single-solve invariant `μ ← δ(applied V)`; the `−η_s·V` brake **and** the `[[v_n]]=0` weld on shared faces | ✔ durable, matches references |
| **1′** | **R-001 minimal** (if R-1601 must stay) — accumulate shared-QP slip from the **macro** `Vᴮ` (move the `δ +=` into `EvaluateADER_LSW`, skip it in the iterator for `i ≥ diag_num_local_fault_qps_`). | `μ ← δ(applied V)` only | partial (two solves remain) |
| **2** | **R-004 fix** — once step 0 confirms, fix the byNODES face-neighbor vdof mapping (or keep the scalar exchange). Guard with the 2-rank bit-identity unit test (R-004). | correct neighbor predictor | ✔ |
| **3** | **R-003** — rewrite DEBUG doc §7: drop the frame fix; point at R-001. | doc consistency | ✔ |
| **4** | **R-005** (optional) — weld `[[v_n]]=0` in the shared corrector too (drdg3d-style), as defense-in-depth. | normal channel can't open | ✔ belt-and-suspenders |

**Why either of step 1 or step 2 stops the runaway** (and why step 1 is the durable one): step 2 removes the *spark* (`Vᴬ≈Vᴮ`, δ consistent); step 1 restores the *brake* (the predictor's `−η_s·Vᴬ` is fed back, so the loop is dissipative even if the exchange is still imperfect — you get a bounded answer). Do both: step 1 for physical correctness of the coupling, step 2 for a correct neighbor state.

## Appendix B — Round-2 re-audit: updated DEBUG doc + `test_ghost_exchange_bynodes_vs_scalar` (2026-05-24)

Re-review of `DEBUG_speckle_normal_velocity_jump_2026-05-23.md` after it adopted the R-001/R-004 framing and added the R-004 confirmation test. **I rebuilt nothing and ran the committed-but-untracked test directly.**

### Verified correct (reproduced, not taken on faith)
- **The R-004 test is sound and genuinely RED.** `mpirun -np 2 ./seas_test_ghost_exchange_bynodes_vs_scalar` → `FAILED … 1836 mismatched values` (918/rank). The test builds `pfes` (vdim=1) and `pfes_full` (vdim=NUM_STATE, `Ordering::byNODES`) **exactly as the ctor does** (`wave_operator.inl:149,160-164`), fills a `1e6·rank+1e3·c+i` fingerprint, runs both exchanges on the same data, and compares the production unpack `src_full[c*n_fn+j]` against the battle-tested per-component scalar exchange. The comparison faithfully reproduces the production read (`wave_operator.inl:2050-2058, 2170`). A RED result means the byNODES unpack reads the wrong neighbor value — not a test artifact.
- **The doc's specific example is exact.** Doc: "`(c=0,j=27)` returns `encode(c=1,i=108)` instead of `encode(c=0,i=135)`." Observed (rank 1): `byNODES=1108 scalar=135` ⇒ `1108 = encode(rank0,c=1,i=108)`, `135 = encode(rank0,c=0,i=135)`. The byNODES layout returns a **wrong component _and_ wrong dof** — `src[c*n_fn+j]` is not the FaceNbrData layout for a byNODES vector space. R-004 (POSSIBLE in the Round-1 review) is now **CONFIRMED**; the doc reports it accurately.
- **Fix order R-004→R-001 is correct and improves on the Round-1 roadmap.** The doc's rationale (consuming `I_imp` on shared faces while the predictor ghost is still scrambled would feed the corruption straight into the bulk → the `tau=4e28→SIGABRT` R-1601 was added to prevent) is right. Round-1 §A.4 listed R-001 primary before R-004 fix; **adopt the doc's order** (R-004 first). Either fix alone bounds the *runaway*, but R-004 must precede the "consume `I_imp`" variant of R-001.
- **The §7/§8 "SUPERSEDED" redirects away from the frame are correct** (frame built identically in predictor `:2199` and macro `:4805`).

### New findings on the updated document

#### [D-001] MODERATE — Blast radius understated: R-004 is not SAFS-only
The scrambled byNODES exchange is in `EvaluateBulkAtFaultQPsCanonical`'s **shared branch** (`wave_operator.inl:2025-2059`), which is called by **all four** substep drivers — `grep` confirms `tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`, `spatial_dyn_driver.cpp` all route through `ComputeADERSubStepStates`/`AdvanceWithSubStepStates`. So **every** substep driver running at `np>1` with shared *fault* faces reads scrambled neighbor predictor state on those faces. It has not blown up TPV runs because (a) the **applied flux** on shared faces uses the correct scalar exchange (`ghost_gf_`, `:4639-4650` — the test confirms that path is right), so only the iterator's shared-QP *slip accumulation* consumes the scramble (R-001 decoupling), and (b) TPV102/104 are rate-and-state (no free-slide-under-tension), so a corrupted predictor yields a bounded-but-wrong V rather than a runaway. **Implication the doc omits:** any `np>1` substep result on a mesh with shared fault faces has *suspect shared-fault slip / state-weakening* since R-1601 landed — SAFS is the catastrophic case, not the only affected one. The doc should state this scope so the fix and re-validation cover TPV102/104/205, not just SAFS.

#### [D-002] MODERATE — Unreconciled tension: the confirmed all-component scramble vs the retained "velocity-channel-only" claim
The test shows the scramble corrupts **stress** components too (`c=0`/`SXX` reads `c=1`'s value). But §2.2/§2.4/§5 still assert the contamination is **velocity-only** (`sn_sterm ≈ 0`, "−3.28 Pa", "only the velocity channel is contaminated, not the bulk stress") and use it to *refute* "corrupted ghost bulk stress." These cannot both be taken at face value: `sn_sterm = η_p(σ_n⁺/Z_p⁺ + σ_n⁻/Z_p⁻)` reads the **neighbor** normal stress `σ_n⁺`; its near-perfect cancellation to −3.28 Pa requires `σ_n⁺ ≈ −σ_n⁻`, which a scrambled `σ_n⁺` would destroy (note `σ_n/Z_p ~ 3 m/s`, so `η_p·σ_n/Z_p ~ 24 MPa` per side — *not* negligible unless the two sides genuinely cancel). Either the seed-time `sn_sterm≈0` is an early-time/coincidental-cancellation artifact (stresses small before the rupture builds), or the scramble's effect on the canonical normal stress is being mis-attributed. **Action:** the `[SLIP]` trace already prints `sn_sterm` — check whether it grows once bulk stresses build; reconcile §2.2/§2.4/§5 with R-004 (the "stress channel is clean" argument in §5 is now suspect). Not fatal to the R-001/R-004 conclusion, but a live loose end in a doc that is now the canonical record.

#### [D-003] LOW — Refuted narrative left inline (readability/safety for the fix agent)
§2.3, §2.4, §3, §6, and §7 still present the frame-leak / tension story in full, current-tense, with only bolted-on "SUPERSEDED" banners. A `/code-fix` agent (or a new reader) skimming §2.4 ("~24° effective `can_n` error") or §7 ("make the canonical fault frame geometrically correct") can act on the refuted direction. Move the refuted sections into a clearly-fenced "Superseded history" appendix, or prefix each with a one-line `⛔ REFUTED — see §8` so they cannot be mistaken for current guidance.

#### [D-004] LOW — Citation drift on the slip accumulator
The sole slip-accumulation site is `tpv205_substep_iterator.cpp:125-126` (`d.slip1 += s.V1*dt_sub`). The TL;DR R-008 block (line 63) and §4 (line 343) cite `:122`, which is actually `s.sigma_n_corr = s.sigma_n_trial;`. §8 correctly cites `:125`. Normalize all references to `:125-126` — the doc is now the root-cause record and feeds a fix agent.

#### [D-005] LOW — R-004 fix is under-specified; prefer the proven revert, and confirm the guard runs in the aggregate
§8.5(1) offers "read through the vector space's face-neighbour vdof map" *or* "revert to per-component scalar." The scramble is **not** a simple byVDIM transpose (the test shows component *and* dof mixing), so a hand-rolled byNODES stride is high-risk. Recommend: **make the proven per-component scalar exchange the fix** (correctness over the `8·O`-collective saving); if perf later demands batching, derive indices from MFEM's `pfes_full->GetFaceNbrElementVDofs(...)`, never a literal `c*n_fn+j`. The test target exists (`Makefile:3518`, `test-ghost-exchange-bynodes-vs-scalar`) but is an untracked file and I did not confirm it is in the parallel regression aggregate — wire it in so it actually gates (it must flip RED→GREEN with the fix).

### Round-2 verdict
The document's **core conclusion is now correct and empirically confirmed** (R-001 dual-solve + R-004 byNODES scramble, fix order R-004→R-001). Remaining issues are scope (D-001), one analytical loose end (D-002), and documentation hygiene (D-003/D-004/D-005) — none overturn the diagnosis. **PASS WITH FIXES.**

## Appendix C — Review of the applied R-004 fix (2026-05-24)

Adversarial review of the committed R-004 fix (`git diff miniapps/seas/dynamic/wave_operator.inl`), the updated guard test, and the Makefile wiring. **I built and ran everything below**, I did not take the diff on faith.

### What the fix does
Instead of the recommended revert-to-scalar, the implementer kept the **batched** R-1601 exchange (one collective) and fixed only the **unpack**: the pre-fix slab read `nbr_data[c][nbr_idx*ndof_per_el_+i]` (which assumed a component-major `src[c*n_fn+j]` layout) is replaced by a layout-agnostic read through MFEM's own vdof map:
```cpp
Array<int> nbr_vdofs;
pfes_full_state_->GetFaceNbrElementVDofs(nbr_idx, nbr_vdofs);   // byNODES ⇒ nbr_vdofs[c*ndof+i]
...
s_nbr += shape2(i) * nbr_all[nbr_vdofs[c * ndof_per_el_ + i]];
```

### Verified correct (with evidence I generated)
- **Guard test GREEN.** `mpirun -np 2 ./seas_test_ghost_exchange_bynodes_vs_scalar` → `fix_mismatch=0` at every (element, component, dof); the test now compares the *production* read pattern `src_full[fv[c*ndof_e+i]]` against the scalar-exchange ground truth `scalar_nbr[c][sv[i]]` (not the abandoned slab formula), and informationally confirms the old formula still gets `918` reads/rank wrong. The test was correctly updated to follow the fix, and is wired into the `test:` aggregate (Makefile:940) **and** `test-v92-regression-gates`.
- **Production recompiles clean.** Forced rebuild of `dynamic/wave_operator.o` (touch + make) — no errors; no dangling `nbr_data` reference left behind.
- **No regression.** `test-interior-vs-shared-branch-live` (np=2) → `worst relative |k_serial − k_parallel| = 0.0` over all 9 components, PASSED. The shared-fault flux still matches the interior/serial branch to machine zero.
- **Architecturally sound choice.** Keeping the batched exchange + fixing the unpack preserves the R-1600 single-collective contract and avoids the `8·O`-collective regression a scalar revert would cost. The exchange (line 2037) is untouched; only the read changed.
- **byte-exact at np=1** (the shared branch is gated on `n_shared_total>0`, so it never executes serially).

### Findings

#### [R4F-001] MODERATE [wave_operator.inl:EvaluateBulkAtFaultQPsCanonical] — release-critical indexing invariants are `MFEM_ASSERT` (debug-only) and the stride uses the global `ndof_per_el_` instead of the element-local `ndof2`
**Category:** ASSUMPTION / ROBUSTNESS

**Description:** The new read `nbr_all[nbr_vdofs[c * ndof_per_el_ + i]]` is correct **only if** (a) `nbr_vdofs.Size() == NUM_STATE*ndof2` and (b) `ndof_per_el_ == ndof2` (so the `c * ndof_per_el_` stride matches the vdof array's actual `[c*ndof2 + i]` byNODES layout). Both invariants are guarded by **`MFEM_ASSERT`** (`ndof2==ndof_per_el_` at the pre-existing line ~2103; the new `nbr_vdofs.Size()` check just below), which is **compiled out in a release (`NDEBUG`) MFEM build**. If either ever fails (a heterogeneous-order ghost), the index silently lands on the **wrong slot or out of bounds — re-introducing exactly the kind of scramble R-004 just fixed, with no diagnostic.** This is the *same* release-safety gap the codebase already decided to close: the macro path's analogous ndof2 check was deliberately promoted ASSERT→VERIFY under **R-1508** (`ComputeADERSharedFaceFluxRHS:4761`, comment: "a Release build still fails loud … rather than silently producing wrong shape evaluations"). The fix should follow that precedent. Separately, the stride should use the element-local `ndof2` (already computed at ~:2102) — which is what the *guard test itself uses* (`fv[c*ndof_e+i]`, `ndof_e = sv.Size()`) — rather than the global member `ndof_per_el_`, removing the hidden dependence on invariant (b) entirely.

**Trigger:** Release build (`MFEM_DEBUG=NO`) with any future heterogeneous-order fault-adjacent ghost element; or a refactor that changes `ndof_per_el_`'s meaning.

**Actual behavior:** silent wrong-slot read (no abort) in release if the homogeneity invariant breaks.

**Expected behavior:** fail loud in release; index self-consistently with the vdof array's own per-element size.

**Suggested fix:**
```diff
-            Array<int> nbr_vdofs;
-            pfes_full_state_->GetFaceNbrElementVDofs(nbr_idx, nbr_vdofs);
-            MFEM_ASSERT(nbr_vdofs.Size() == NUM_STATE * ndof2,
-                        "EvaluateBulkAtFaultQPsCanonical: neighbour vdof count "
-                        "!= NUM_STATE * ndof2 (R-004 byNODES map).");
+            Array<int> nbr_vdofs;
+            pfes_full_state_->GetFaceNbrElementVDofs(nbr_idx, nbr_vdofs);
+            MFEM_VERIFY(nbr_vdofs.Size() == NUM_STATE * ndof2,   // release-checked, R-1508
+                        "EvaluateBulkAtFaultQPsCanonical: neighbour vdof count "
+                        << nbr_vdofs.Size() << " != NUM_STATE * ndof2 = "
+                        << NUM_STATE * ndof2 << " (R-004 byNODES map).");
```
```diff
-                     s_nbr += shape2(i)
-                              * nbr_all[nbr_vdofs[c * ndof_per_el_ + i]];
+                     s_nbr += shape2(i)
+                              * nbr_all[nbr_vdofs[c * ndof2 + i]];   // stride = element-local ndof2
```

**Test case:** the existing `test_ghost_exchange_bynodes_vs_scalar` stays GREEN after the change (its stride is already `ndof_e == ndof2`), proving the stride swap is a zero-risk equivalence on the homogeneous space; the `MFEM_VERIFY` then fires at runtime under a heterogeneous ghost instead of scrambling silently. (A dedicated heterogeneous-order test is out of scope — the codebase universally assumes uniform L2 order — so the guard is the VERIFY + the unchanged GREEN test.)

#### [R4F-002] LOW [POSSIBLE] [wave_operator.inl:EvaluateBulkAtFaultQPsCanonical] — signed/encoded vdofs not decoded; correctness tied to "FEC is always L2"
**Category:** ASSUMPTION

**Description:** `GetFaceNbrElementVDofs` can return **negative** indices encoding orientation sign-flips for oriented FE spaces; `nbr_all[nbr_vdofs[...]]` indexes directly with no `FiniteElementSpace::DecodeDof`/sign handling. For the **L2 (DG)** collection the wave operator uses, dofs are never sign-encoded, so this is **safe today** (and consistent with the test, which also doesn't decode). But it silently couples the fix's correctness to "the space is L2"; a future vector/oriented space would make `nbr_all[negative]` an OOB read.

**Trigger:** reuse of this path with a sign-encoding FE space (not currently the case).

**Suggested fix (defensive, no-op for L2) or document the assumption:**
```diff
-                     s_nbr += shape2(i)
-                              * nbr_all[nbr_vdofs[c * ndof2 + i]];
+                     {
+                        const int vd = nbr_vdofs[c * ndof2 + i];  // L2: always >= 0
+                        s_nbr += shape2(i) * nbr_all[vd >= 0 ? vd : -1 - vd];
+                     }
```
Or add a one-line comment asserting the L2 (no-sign) assumption next to the read.

#### [R4F-003] LOW [wave_operator.inl:EvaluateBulkAtFaultQPsCanonical] — unnecessary full-vector deep copy per sub-step; the stated rationale no longer holds
**Category:** QUALITY / PERF

**Description:** `Vector nbr_all(src.Size()); nbr_all = src;` deep-copies the entire FaceNbrData (`NUM_STATE * num_face_nbr_dofs`) on every call (O per macro step). The comment justifies it as "pre-R-1601 alias-safety (subsequent `q_gf_full` reads could invalidate the storage)" — but the new code never re-reads or re-exchanges `q_gf_full` inside the loop, so `src` (the `FaceNbrData()` reference) stays valid for the whole loop. The copy is dead defensive cost.

**Suggested fix:**
```diff
-         Vector nbr_all(src.Size());
-         nbr_all = src;
+         // src (FaceNbrData) stays valid for the whole loop below — q_gf_full is
+         // not re-exchanged or modified — so alias it directly (no per-call copy).
+         const Vector &nbr_all = src;
```
(If a defensive copy is still desired, keep it but correct the comment — the "subsequent reads invalidate" rationale is stale.)

### Round-3 verdict
The R-004 fix is **correct and verified** (guard test GREEN, recompiles, branch-equivalence rel 0, batched collective preserved). The vdof-map approach is a legitimate alternative to the scalar revert and is arguably better (no collective regression). Findings are **hardening only** — no correctness defect found in the fix. **PASS WITH (minor) FIXES**: apply R4F-001 (release-safety + self-consistent stride, with the R-1508 precedent) before relying on this in production builds; R4F-002/R4F-003 are optional.

> Per-finding note for the fix agent: these are appended here (the speckle review chain) rather than the repo-root `REVIEW.md`, which holds the in-flight diagnostics handoff. Move R4F-001 into `REVIEW.md` if you want it consumed by `/code-fix` directly.

## Unreviewed Areas
- The exact MFEM internal `FaceNbrData()` layout for a `vdim=NUM_STATE` byNODES `ParGridFunction` was reasoned about but not traced into MFEM source — this is why R-004 is POSSIBLE, not CRITICAL. The R-004 unit test resolves it deterministically.
- `ComputeADERSubStepStates` / `ComputeADERTimeIntegrated` factorial recursions were verified algebraically equivalent (`:1381-1394`) but not numerically diffed at the seed QP — out of scope for static review.
- The cross-rank reconcile (`ExchangeAndPairSharedFaultQPs`, Pass 2 at `:4670-4676`) was not audited; it forces the *output* fields equal post-solve and does not change the slip-weakening decoupling, so it is downstream of R-001.
