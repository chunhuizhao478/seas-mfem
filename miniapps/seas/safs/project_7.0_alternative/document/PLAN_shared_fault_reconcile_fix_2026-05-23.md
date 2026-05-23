# Implementation Plan: method-invariant shared-fault cross-rank reconcile

## Overview
On a shared (cross-rank) fault face the per-QP friction state (`DOFData`) is
computed **redundantly** on the two ranks that own the two sides. Their inputs
(the bulk Q at the QP) are NOT bit-identical — the same physical QP is
interpolated through two different element parametrisations (`shape1∘Loc1` vs
`shape2∘Loc2`), so they differ at ~1e-14. R-701 deleted the v5 R-501
owner-broadcast on the (false) premise that the canonical frame makes the inputs
bit-identical. The Frontera XRANK trace (job 7747036) proved this is the SAFS
blow-up: at the slip-onset kink `V_abs=max(0,(|τ|−τ_str)/η_s)` the 1e-14 tips one
rank to slip and the other to lock → full desync → blow-up; the R-101 guard
aborts first (t≈0.455 s).

This plan delivers, in one document, the three pieces requested:
1. **Reproducible local test** — a deterministic np=2 unit test that injects a
   1-ULP cross-rank input difference at the slip-onset threshold and reproduces
   the desync (RED before the fix), for BOTH friction laws.
2. **General (method-invariant) fix** — make the shared-fault `DOFData` and
   imposed state **single-valued** across the two ranks (owner computes,
   broadcasts; non-owner adopts), for ALL friction laws, before flux assembly.
3. **Correctness guard** — the cross-rank consistency verify, kept as a runtime
   guard (method-invariant), expected to report `max_rel_diff == 0` after the
   fix, so any future regression aborts loudly.

Source of truth for the diagnosis: `spatial_dynamic_rupture_speckle_blowup_2026-05-22.md`.
Design review (findings R-001..R-006): `spatial_dynamic_rupture_reconcile_review_2026-05-23.md`.
This document SUPERSEDES the fix sections of
`PLAN_dip_channel_v1_cross_rank_stability_2026-05-22.md` (covariant-decomposition
direction — dropped as a symptom-treatment).

## Constraints
- **Method-invariant (review R-001):** the reconcile applies to LSW,
  LSW_ForcedRupture, AND RateState in the shared-fault ADER path. No
  per-law / per-rake / kink-proximity gating (that re-introduces
  problem-specificity).
- **Regression contract (REVISED, R-001) — needs the byte-exact relaxation that
  the method-invariant directive implies:** any cross-rank reconcile makes the
  non-owner adopt the owner's value, changing TPV102/104 by ~1e-14, so
  bit-exactness vs the pre-fix binary is impossible. New contract:
  TPV102/104/205 + BP5 **physically-exact** (`worst_rel ≤ 1e-13` vs pre-fix) AND
  now **cross-rank bit-identical** (R-101 `max_rel_diff == 0`). The latter is a
  strict correctness improvement over the current cross-rank-inconsistent-at-1e-14
  baseline. (Overrides the CLAUDE.md "TPV/BP5 byte-exact" rule per the directive;
  recorded here.)
- **MPI-collective-safe (review R-005, R-1600 deadlock class):** every rank with
  shared fault faces participates consistently; ranks without must not deadlock.
- **Files Requiring Extreme Care** (`CLAUDE.md`): `wave_operator.inl/.hpp`,
  `fault_face_flux.cpp`. Full verification suite required.
- **Owner rule (review R-003):** owner = the **lower MPI rank** of the two
  sharing the face (globally unique, deterministic, geometry-independent —
  immune to the θ≈90° degenerate band that `elem1_on_plus` is fragile in).
  `elem1_on_plus` stays ONLY for the +/- flux role.
- **No hardcoded constants** (derive tolerances from η_s/σ_n/μ); **Brent** for
  rate-state friction (unchanged); the reconcile does not touch the solvers.

## Background math (what each rank computes, and the seed)
Per shared fault QP, in the bit-identical canonical frame `(can_n, can_t1,
can_t2)` (frame is static + bit-identical — confirmed; see diagnosis):
```
I_self_can = Tinv_can · I_self     I_nbr_can = Tinv_can · I_nbr
(Q_plus, Q_minus) = elem1_on_plus ? (I_self_can, I_nbr_can) : (I_nbr_can, I_self_can)
EvaluateADER_LSW / EvaluateADER : (Q_plus,Q_minus) → DOFData{V1,V2,τ*_corr,σn_corr,slip*}, I_imp_plus/minus(canonical)
```
Cross-rank, `I_self_can(A) == I_nbr_can(B)` only to ~11 digits (interpolation FP,
~1e-14). The friction map is Lipschitz away from the kink (rate-state: bounded
1e-14 divergence — invisible but real) and NON-smooth at `|τ|=τ_str` (LSW: 1e-14
→ O(1) branch split). The fix removes the divergence at the source: the
shared-QP `DOFData` + `I_imp` are computed once (owner) and copied — so both
ranks hold identical values regardless of the 1e-14 input gap.

---

## Phase 1: Reproducible local test (injection-based, np=2, method-invariant)

### Goal
A standalone np=2 unit test that DETERMINISTICALLY reproduces the cross-rank
desync at the slip-onset threshold by injecting a 1-ULP input difference, and
asserts the two ranks' `DOFData` are bit-identical — RED before Phase 2, GREEN
after — for BOTH RateState and LSW (proving method-invariance, review R-001).

### Files to Create
- `tests/unit/test_shared_fault_reconcile_cross_rank.cpp` — the reproducer.

### Files to Modify
- `Makefile` — 5 entries mirroring the tilted-test wiring (SRC at ~323, OBJ at
  ~623, link target `seas_test_shared_fault_reconcile_cross_rank` with the same
  object set as `seas_test_rupture_tilted_fault_serial_vs_parallel` (Makefile
  ~1671), compile rule (~2719), run target `test-shared-fault-reconcile-cross-rank`
  invoking `$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2`). NOT in `make test`.
- `dynamic/fault_face_flux.hpp` (TEST-ONLY hook, behind `#ifdef SEAS_TEST_INTERNAL`
  which already exists): add an optional static debug perturbation
  `static real_t s_seas_test_qplus_xz_perturb_ulp = 0.0;` consumed in
  `ComputeTrialTraction` to add `s*ULP*Q_plus[SXZ]` to `Q_plus[SXZ]` ONLY when
  set — so the test can inject a controlled, rank-specific 1-ULP difference
  without faking the whole bulk. If `SEAS_TEST_INTERNAL` is not defined the hook
  compiles out (production byte-exact). (Alternative if a fault-flux hook is
  undesirable: inject by perturbing one rank's bulk Q DOF directly in the test
  before the step — preferred if it reproduces; decide at implementation.)

### Detailed Requirements
1. **Fixture.** Reuse `BuildTwoTetFaultMesh()` (planar y=0, 2 tets) and the
   `cy<0→rank0 / cy≥0→rank1` partition from
   `test_rupture_multistep_serial_vs_parallel.cpp`. One shared fault face, 3 QPs.
2. **Drive to the slip-onset threshold.** Set per-DOF LSW params and a uniform
   shear load such that `|τ_total|` crosses `μ_s·σ_n` during the run: derive the
   load from `TPV205Params` so the QP is AT the kink at some step (assert the
   run reaches a step where `0 < V_abs < V_small` on at least one rank, i.e. the
   onset). For the rate-state leg use `TPV102Params` total-stress setup.
3. **Inject a 1-ULP cross-rank seed.** On rank 0 only, set the test perturbation
   so its `Q_plus[SXZ]` differs from rank 1's matching input by exactly 1 ULP
   (`std::nextafter`). This emulates the interpolation seed deterministically.
4. **Run** `kNSteps` ADER-2 steps through the threshold (≥ a few steps past
   onset).
5. **Assert (the captured bug).** After each step, gather the shared-QP `DOFData`
   from both ranks (reuse `VerifySharedFaultDOFDataConsistency(tol, &wr, &wf,
   /*abort=*/false)` from the diagnostic commit) and record `worst_rel`. The
   test asserts `worst_rel == 0` (bit-identical) at every step.
   - **Before Phase 2:** FAILS — LSW leg desyncs to `worst_rel≈1.0`; rate-state
     leg shows `worst_rel≈1e-14` (bounded but nonzero — proves the defect is
     method-invariant, just sub-threshold for rate-state).
   - **After Phase 2:** PASSES — `worst_rel == 0` for both laws.
6. **Run BOTH laws** in the same executable (loop over {RateState, LSW}); the
   headline assertion is per-law. This is what proves method-invariance.

### Interfaces
- Reuse `VerifySharedFaultDOFDataConsistency(real_t, double*, int*, bool)`
  (added in commit a9bd4d2).
- `static real_t FaultFaceFlux::s_seas_test_qplus_xz_perturb_ulp` (TEST-ONLY).

### Edge Cases to Handle
- np≠2 → SKIPPED, return 77.
- A run where the QP never reaches onset → assert nucleation/onset occurred,
  else FAIL "load did not cross the slip-onset threshold".

### Acceptance Criteria
- [ ] `make test-shared-fault-reconcile-cross-rank` builds + runs at np=2.
- [ ] Before Phase 2: RED (LSW `worst_rel→O(1)`; rate-state `worst_rel~1e-14`>0).
- [ ] After Phase 2: GREEN — `worst_rel == 0` for BOTH laws, every step.
- [ ] Production build (`SEAS_TEST_INTERNAL` undefined): the perturbation hook
      compiles out; TPV/BP5 unaffected.

### Dependencies
- Depends on: nothing. Required by: Phase 2 (oracle), Phase 4.

---

## Phase 2: General method-invariant reconcile (the fix)

### Goal
After the friction solve on a shared fault face, the two ranks hold
**bit-identical** `DOFData` and assemble their side's flux from the **owner's**
imposed state — for ALL friction laws — so the slip/lock decision is made once
and copied. The desync is impossible regardless of the 1e-14 input gap.

### Files to Modify
- `dynamic/wave_operator.inl` `ComputeADERSharedFaceFluxRHS` (LSW branch
  ~4671-4797 AND the rate-state `else` at ~4779): restructure the shared-fault QP
  handling into **two passes** with a reconcile between them (review R-002).
- `dynamic/wave_operator.hpp`: declare the reconcile helper + a per-QP assembly
  buffer struct.
- (Reuse) the R-101 record gather/pair logic from
  `VerifySharedFaultDOFDataConsistency` (`:5470+`): factor it into a shared
  helper so the reconcile and the verify use ONE proven matcher (DRY).

### Detailed Requirements
1. **Factor the cross-rank pairing** out of `VerifySharedFaultDOFDataConsistency`
   into:
   ```cpp
   // Gathers one record per local shared-fault QP (face-vertex-key + qp_idx +
   // rank + payload[NPAY]) via MPI_Allgatherv, pairs each QP across the two
   // ranks by (face-key, qp_idx), and invokes cb(local_idx, peer_payload,
   // peer_rank) for each local QP that has a peer.  Returns #unpaired.
   int WaveOperator<MeshType>::ExchangeAndPairSharedFaultQPs(
        int npay, const std::vector<double>& local_payload,
        const std::function<void(int local_qp, const double* peer_payload,
                                 int peer_rank)>& cb) const;
   ```
   The verify calls it with `npay=8` (the 8 fields) and a compare callback; the
   reconcile calls it with `npay = 8 + 2*NUM_STATE` (DOFData + I_imp_plus_can +
   I_imp_minus_can) and an overwrite callback. Identical face-key construction
   (sorted global vertex IDs) as today.
2. **Per-QP assembly buffer.** Define
   ```cpp
   struct SharedFaultQPAssembly {
      int    dof_offset1, ndof;            // local Elem1 assembly target
      bool   elem1_on_plus;
      real_t w;                            // ip.weight * |J_F|
      real_t can_n[3], can_t1[3], can_t2[3];
      real_t shape1[MAX_NDOF];             // fe1 shape at the QP
      real_t I_imp_plus_can[NUM_STATE];    // owner-reconcilable
      real_t I_imp_minus_can[NUM_STATE];
      int    dof_idx;                      // into fault_dof_data_
   };
   ```
3. **Pass 1 (compute, no assembly).** Loop shared fault QPs exactly as today
   through the friction dispatch (LSW / LSW_ForcedRupture / rate-state) to fill
   `fault_dof_data_[dof_idx]` and `I_imp_plus/minus` (canonical frame). Store the
   `SharedFaultQPAssembly` for each QP. Do NOT assemble `rhs` yet.
4. **Reconcile.** Build the payload per QP = the 8 DOFData fields +
   `I_imp_plus_can` + `I_imp_minus_can`. Call `ExchangeAndPairSharedFaultQPs`.
   In the callback, if `peer_rank < my_rank_` (peer is owner), OVERWRITE this
   rank's `fault_dof_data_[dof_idx]` (8 fields) and the buffer's
   `I_imp_plus_can / I_imp_minus_can` with the peer's payload. (Owner = lower
   rank, R-003.) `MFEM_VERIFY` every local shared-fault QP was paired (else the
   classification is asymmetric — fail loud).
5. **Pass 2 (assemble from reconciled state).** For each buffered QP, rotate the
   (possibly-overwritten) `I_imp_*_can` to global via `T_can` (rebuilt from the
   buffered `can_*`; bit-identical on both ranks) and assemble
   `rhs[c·ndof_total_ + dof_offset1 + i] += assemble_sign · w · shape1[i] · F_h_side[c]`
   exactly as the current inline code, but consuming the **reconciled** imposed
   state (review R-002, R-004).
6. **Method-invariant:** steps 3-5 run for every friction law; no gate. The
   rate-state `else` branch is included.
7. **Collective safety (R-005):** `ExchangeAndPairSharedFaultQPs` uses the same
   global `MPI_Allreduce(any_shared)` short-circuit + `MPI_Allgatherv` as the
   verify (every rank participates; ranks with no shared fault QPs contribute a
   zero-length buffer). No participation depends on local fault presence.

### Interfaces
- `int ExchangeAndPairSharedFaultQPs(int npay, const std::vector<double>&,
   const std::function<void(int,const double*,int)>&) const;` (new, in
   `wave_operator.hpp`, parallel-only).
- `struct SharedFaultQPAssembly` (new, file-scope in the .inl).

### Edge Cases to Handle
- Rank with shared fault faces but a QP whose peer never appears (mesh
  classification asymmetry) → `MFEM_VERIFY` fail (same policy as the verify's
  `n_unpaired>0`).
- A shared face where `elem1_on_plus` is the same on both ranks (the θ≈90°
  side-fix edge) → owner-by-rank still picks exactly one owner, so the reconcile
  is well-defined even if the +/- role were momentarily wrong (defense in depth).
- Serial build / np=1 → no shared faces → the two-pass degenerates to the
  current single-rank path (no exchange); byte-exact.

### Acceptance Criteria
- [ ] Phase-1 test GREEN for both laws (`worst_rel == 0` every step).
- [ ] TPV102/104/205 + BP5: physically-exact `worst_rel ≤ 1e-13` vs pre-fix AND
      R-101 `max_rel_diff == 0` (cross-rank bit-identical — was ~1e-14).
- [ ] The local fault/TPV suite (fault-basis trio, shared-fault role/dof-data
      consistency, interior-flux-path ×3, godunov identity, tpv102
      locked/absorbing/pepper/ader-smoke, multistep & tilted serial-vs-parallel,
      sign-flipped truth table) all green.
- [ ] np=2,4,8 smoke (incl. a partition where some ranks have NO shared fault
      faces) — no deadlock (R-005).

### Dependencies
- Depends on: Phase 1. Required by: Phase 4.

---

## Phase 3: Correctness guard (runtime, method-invariant)

### Goal
A permanent, method-invariant runtime guard that the shared-fault `DOFData` is
cross-rank-consistent, so any future regression (a missed field in the payload,
an exchange bug, a new law that bypasses the reconcile) aborts loudly instead of
silently desyncing.

### Files to Modify
- `dynamic/wave_operator.inl` `VerifySharedFaultDOFDataConsistency`: (a) refactor
  it onto the shared `ExchangeAndPairSharedFaultQPs` helper from Phase 2 (one
  matcher); (b) update the stale R-501 wording in the abort message to point at
  the Phase-2 reconcile; (c) keep the default abort tol at `1e-10` but note that
  POST-fix the expected `max_rel_diff` is `0` (the guard now also catches any
  nonzero divergence as a real regression, not roundoff).
- `drivers/spatial_dyn_driver.cpp` (and tpv104/tpv205 drivers if they call it):
  keep the existing gated cadence (`step==0 || (nucleation.enabled && t≤T_nuc_s
  && step%100==0)`); add a CLI/env to raise the cadence for debugging
  (`SEAS_VERIFY_XRANK_EVERY=N`), default unchanged.

### Detailed Requirements
1. The guard MUST remain method-invariant (it already checks all shared QPs
   regardless of law — preserve).
2. Post-fix the guard is the regression sentinel: a single `max_rel_diff > 0`
   (beyond a derived ULP floor) is a real bug, not noise. Document this in the
   verify header and the abort message.
3. The guard's per-call cost (Allgatherv) is unchanged; its cadence is unchanged.
   (The per-substep reconcile in Phase 2 is the always-on consistency mechanism;
   this guard is the periodic audit.)

### Acceptance Criteria
- [ ] Post-fix Dc2 Frontera run: the guard reports `max_rel_diff == 0` through
      and past the old t≈0.455 s abort point; no abort.
- [ ] Inject a deliberate regression (skip one field in the Phase-2 payload) →
      the guard aborts on that field. (A negative test in the Phase-1 unit test:
      with the reconcile disabled via a test switch, the guard trips.)

### Dependencies
- Depends on: Phase 2. Required by: Phase 4.

---

## Phase 4: Full regression + Frontera re-validation

### Goal
End-to-end: the Dc2 mesh advances past t=1.0 s (no R-101 abort, no blow-up);
the fix is confirmed config-agnostic.

### Detailed Requirements
1. Local regression battery (Phase-2 AC list) green; document unchanged
   pre-existing failures (adjacent-triangle pepper-bug, r101 missing-precondition,
   macOS MUMPS Bus error).
2. Re-run the `SEAS_DIAG_XRANK` Dc2 sbatch: at the onset QP both ranks now show
   identical V1/V2/τ*_corr through t=0.478 s; R-101 passes; run reaches tfinal
   (or well past t=1.0 s).
3. Re-run (or have the user run) the physical `D_c=1.0` baseline config to
   confirm the fix is not tuned to Dc2.

### Acceptance Criteria
- [ ] Local regression green (modulo documented pre-existing).
- [ ] Frontera Dc2 AND D_c=1.0: no R-101 abort, no blow-up; `V_max` peaks then
      decreases (rupture, not runaway).

### Dependencies
- Depends on: Phase 2, Phase 3.

## Testing Strategy
- **Primary local oracle:** Phase-1 injection test (np=2, both laws) — RED→GREEN,
  proves method-invariance and that the reconcile absorbs the seed.
- **Guard:** Phase-3 R-101 verify — post-fix `max_rel_diff==0`; negative test
  (disable reconcile) trips it.
- **Regression (physically-exact + cross-rank-identical):** TPV102/104/205 + BP5
  via the listed targets; capture pre-fix outputs, assert `≤1e-13` AND R-101
  `max_rel_diff==0`. Use the stash-rebuild-baseline method to quantify the
  ~1e-14 TPV change and confirm it is the reconcile (non-owner→owner), not a
  physics change.
- **MPI robustness:** np=2/4/8, incl. a no-shared-fault-on-some-ranks partition.

## Risk Assessment
- **Deadlock (R-005, R-1600 class):** the new Allgatherv must be reached by every
  rank unconditionally (mirror the verify's `any_shared` short-circuit). Detect:
  np=4/8 smoke hangs.
- **Per-substep Allgatherv cost (R-006):** correctness-first; if hot on 400r,
  migrate `ExchangeAndPairSharedFaultQPs` to point-to-point over the face-nbr
  topology (same payload). Do NOT add a kink-proximity gate. Detect: step-rate
  vs the pre-fix run.
- **Two-pass refactor regressions:** buffering the assembly context risks an
  indexing slip (dof_offset1/shape1/w). Mitigate: the interior-fault path is
  untouched; the tilted + multistep + adjacent-triangle tests cover shared-face
  assembly; assert Pass-2 reproduces the pre-refactor RHS when the reconcile is a
  no-op (ranks already agree — force via a serial-equivalent np=1).
- **Byte-exact relaxation:** the ~1e-14 TPV change overrides a CLAUDE.md
  non-negotiable; recorded (review R-001) per the method-invariant directive.
- **Tricky existing code:** `ComputeADERSharedFaceFluxRHS` (4671-4797, the
  R-1601 inline fallback + the assembly), `VerifySharedFaultDOFDataConsistency`
  (5364+, the matcher being factored), the canonical-frame reconstruction
  (4695-4709, must stay bit-identical across ranks for `T_can` rotation in
  Pass 2).
