# Code Review: Phase 4 + Phase 5 of the interior/shared fault unification (2026-07-10)

## Review Scope

- Plan: `miniapps/seas/document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md`
  (§Phase 4 "Forced rupture uses one clock", §Phase 5 "Collapse the duplicated
  interior/shared fault code into one routine", and the "PHASE 4 + PHASE 5
  IMPLEMENTED" banner at the top of that file).
- Files reviewed:
  - `dynamic/wave_operator.hpp` (helper declaration :599, FRAME BIT CONTRACT doc,
    R-1003/R-1601 History note)
  - `dynamic/wave_operator.inl` (helper definition :3940; interior call site :4361;
    shared call site :5293; RK pair; diag remaps; stage-avg `_qq` block :4604)
  - `dynamic/friction_substep_iterator.hpp` (`WaveOpLaw()` doc)
  - `drivers/spatial_dyn_driver.cpp` (R-002 Option-A guard removal :860)
  - `tests/unit/test_forced_rupture_iterator_parity.cpp` (`Test_Seam_Interior_Mu_Identical`)
  - `tpv26/REVIEW_tpv26_phase2_3_2026-07-09.md` (R-002 Fix Status row)
  - Collaterals greped for staleness: `dynamic/fault_face_flux.{hpp,cpp}`,
    `spatial/code/spatial_friction.cpp`, `dynamic/spatial_nucleation.{hpp,cpp}`.
- Domain context: `miniapps/seas/CLAUDE.md`, root `CLAUDE.md` (never revert prior fixes
  without citing debug docs),
  `debug_document/tpv104_debug_document/R1601_root_cause_2026-07-10.md`,
  `document/REVIEW_phase2_3_unify_2026-07-10.md` (prior round),
  `tests/parallel/test_fault_frame_bit_census.cpp` (Phase-0 census).
- Baseline for diffs: all session changes are uncommitted, so `git show HEAD`
  (`4e98bd0`) is the pre-Phase-0 state; Phase 2/3 deltas were subtracted using the
  Phase-2/3 review doc.

## Refactor-fidelity verification performed by this review (all clean)

These are the review's own measurements, recorded so the fix round does not redo them:

1. **RK pair byte-identity (hunt item 2).** `sed`-extracted lines 2676–3939 of the
   working `wave_operator.inl` (both RK functions `ComputeFaceFluxRHS` +
   `ComputeSharedFaceFluxRHS`) against lines 2676–3912 of `git show HEAD`: the ONLY
   diff is the Phase-5 helper's leading comment block appended after the shared-RK
   function. The hand-restored RK fault block is **byte-identical** to its original.
2. **Helper vs original interior block (hunt item 1).** Full-function diff of
   `ComputeADERFaceFluxRHS` (HEAD 3913–4996 vs working 4049–5024): frame-build
   ternaries, `Tinv_can` rotation loop, `elem1_on_plus` routing, substep-buffer gate
   condition, LSW v_imp recovery, and the `LSW_ForcedRupture → LSW → RateAndState`
   dispatch (incl. `VerifyForcedRuptureTimeReady` + `GetTime()`) are op-for-op
   identical inside the helper. The only structural changes: (a) routing is now an
   array **copy** into `I_plus_can/I_minus_can` instead of pointer aliases —
   value-identical, no FP change; (b) `T_can` is built at the call site AFTER the
   helper instead of alongside `Tinv_can` before it — `GodunovFlux::BuildRotation`
   and `BuildRotationInverse` are static, stateless pure functions of
   `(can_n, can_t1, can_t2)` (`godunov_flux.cpp:221/263`), so build order cannot
   change any matrix entry; (c) the `SEAS_DIAG_TPV104_FAULT_BASIS` print now runs
   after the helper instead of between frame build and rotation — stderr-only, no
   state read that changed meaning.
3. **Shared branch.** Full-function diff of `ComputeADERSharedFaceFluxRHS` shows
   exactly: the R-005 `MFEM_VERIFY(dt > 0)` at entry, the R-1601-fallback block
   replaced by the helper call (Phase 2+5), the documented [MACRO] comment rewrite,
   the diag remaps, and the Phase-3 ordering-contract + reconcile assertion. Nothing
   else. The qa deferred-buffering and Pass-2 scatter are untouched.
4. **XRANK/C-1s diag remap correctness (hunt item 3).** Helper routes
   `I_plus_can = I_self_can` iff `elem1_on_plus` (`wave_operator.inl:3989`), so the
   remap `self = (elem1_on_plus ? I_plus_can : I_minus_can)` /
   `nbr = (elem1_on_plus ? I_minus_can : I_plus_can)` reproduces the original
   `I_self_can`/`I_nbr_can` values exactly in both branch polarities. Correct.
5. **Frame-bit contract.** Interior passes `!elem1_on_plus`, shared passes
   `qpd.sign_flipped`; the census-mandated non-unification is honored and documented
   in the header. The `should_negate_frame_qq` build at :4604 belongs to the
   `FaultEvalStageAvgMode::Trial/Theta` two-pass averaging algorithm — structurally
   different (ComputeStageState/CompleteFromTrial, no I_imp), unchanged from HEAD,
   legitimately outside Phase-5 scope (see R-105 for the doc consequence).
6. **Phase 4(a) GetTime reachability.** Working tree has exactly ONE physics
   `GetTime()` on the fault path — the helper's inline one-shot
   `LSW_ForcedRupture` arm (:4033); the other four uses (:5365, :5418, :5433,
   :5493) are env-gated diag prints. With the buffer installed the one-shot arm is
   unreachable. Verified.
7. **Phase 4(c) guard removal.** The Option-A `t0_s > 0 && nprocs > 1` MFEM_VERIFY
   is gone from the driver (replaced by an explanatory comment at :860-864); the
   `t0_s >= 0` validation survives TWICE (`spatial/code/spatial_friction.cpp:1494`
   parser + `dynamic/spatial_nucleation.cpp:464` resolver). No test sets `t0_s = 0`
   with `nprocs > 1` expecting an abort (hunt item 7 clean; the only test use is
   `test_forced_rupture_resolver.cpp:76` with `t0_s = 0.5`, serial).
8. **Release-build asserts (hunt item 5).** The helper's `MFEM_ASSERT`s (dt>0,
   dof_idx range) compile out in release, but BOTH ADER entries carry
   `MFEM_VERIFY(dt > 0)` (interior at function top, shared added by R-005) and both
   call sites bound-check `dof_idx` before calling. Belt-and-suspenders — acceptable;
   no finding.
9. **Test reruns on the reviewed tree:** `test_forced_rupture_iterator_parity`
   68/68; `test_friction_substep_iterator_parity` 39/39. Pre-refactor reference
   trajectories still on disk at
   `/Users/chunhuizhao/.claude/jobs/e2a416ea/tmp/refs_pre_phase5/` (4 files).

## Review-generated evidence for the Phase-4 acceptance criterion

The plan's Phase-4 acceptance asks for an np=1 vs np=2 forced-rupture comparison "on
a fault forced across a seam". The implementation gate did not include a multi-rank
**forced-rupture** run (its np=2 run was TPV104 rate-and-state), so this review ran
one. Configs derived from `tpv26/configs/tpv26_spatial_forcedrupture_smoke.toml`
(t0_s = 0.05, and a t0_s = 0 step-ramp variant — the configuration the removed guard
used to FORBID at np > 1):

| run | V_max_global | notes |
|---|---|---|
| np=1, t0_s=0.05 | 0.199649 | matches the implementation gate value |
| np=2, t0_s=0.05 | 0.199649 | all 16 station files **byte-identical** to np=1 — but `shared = 0`: METIS put no seam on the fault, so np=2 does not exercise the seam |
| np=4, t0_s=0.05 | 0.199649 | **6 of 261 fault DOFs on shared faces** (the seam case). 14/16 station files byte-identical to np=1; 2 files differ in ONE line, ONE value, LAST printed digit (max_rel 2.3e-11) — partition FP noise, far below any rupture-time-relevant scale. R-101 init check: 3 shared pairs matched, max_rel_diff = 0; zero reconcile aborts. |
| np=1 vs np=2, t0_s=0 | 0.252328 both | formerly-guarded config runs clean and np-consistently; station files byte-identical |

Verdict on the criterion: **substantively demonstrated** (np≤2 exact, np=4 seam case
equal to the last printed digit). The np≥4 last-digit noise is consistent with the
known, pre-existing np≥4 attractor scope and does not indicate a seam-clock lag.

## Findings

### [R-101] MODERATE [tests/unit/test_forced_rupture_iterator_parity.cpp:150 Test_Seam_Interior_Mu_Identical] — the replacement test is vacuous: it compares one pure function against itself

**Category:** BUG (test) / DEVIATION

**Description:**
`Test_Seam_Interior_Mu_Identical` computes `mu_interior` and `mu_seam` by calling
`spatial::LSWFrictionCoefficient_ForcedRupture` **twice with bit-identical
arguments** (lines 163–168) and asserts their equality at tol 0.0. Two calls to the
same pure function with the same arguments are equal by the C++ abstract machine; the
assertion cannot fail under ANY change to the friction code, the wave operator, the
iterator, or the dispatch. Plan Phase 4 Detailed Requirement 2 asked for a test that
*guards* the one-clock contract ("Replace the bounded-lag test with ... seam and
interior QPs at equal (δ, T_forced, t0) produce bit-identical μ"). As written the
test provides zero regression value while its name and comment claim coverage of the
"unified substep dispatch" — the dangerous part: a future regression that
reintroduces a seam clock (e.g. someone re-adding an inline shared arm reading
`GetTime()`) keeps this test green.

**Trigger:**
Any regression whatsoever; the test passes unconditionally.

**Actual behavior:**
`f(x) == f(x)` asserted 40 times.

**Expected behavior:**
A test that fails if the μ contract the unified dispatch depends on changes. The
honest unit-level guard is a pinned closed-form reference: an independent
reimplementation of the spec composition (`barrier short-circuit; f1 = clamp(δ/d_c);
f2 = forced ramp; μ = μ_s + (μ_d − μ_s)·max(f1,f2)`) compared at tol 0.0. That pins
the algebra the ONE shared call site (post-Phase-5 helper → iterator `StepOneQP_`)
consumes, and any conscious change to the composition must update the pin. (The
seam==interior property itself is structural post-Phase-5 — one call site — and is
enforced at integration level by `test_shared_fault_substep_parity_np2` +
the tet2x2 np=2==np=1 harness; the test comment should say so.)

**Suggested fix:**

```diff
-            const real_t mu_interior =
-               spatial::LSWFrictionCoefficient_ForcedRupture(
-                  delta, kMuS, kMuD, kDc, t, T, t0);
-            const real_t mu_seam =
-               spatial::LSWFrictionCoefficient_ForcedRupture(
-                  delta, kMuS, kMuD, kDc, t, T, t0);
-            TEST_NEAR(mu_seam, mu_interior, 0.0,
-                      "seam == interior mu, bit-exact (same helper, same "
-                      "t_sub_end; unified substep dispatch)");
+            // Independent closed-form reference (mirrors the spec composition
+            // of spatial_friction.hpp:841 op-for-op).  Seam == interior is
+            // structural post-unify-Phase-5 (ONE call site inside
+            // ComputeFaultQPImposedStatesCanonical_ / StepOneQP_); what CAN
+            // regress is the composition itself — pin it bit-exactly.
+            real_t f1;
+            if (delta <= 0.0)      { f1 = 0.0; }
+            else if (delta >= kDc) { f1 = 1.0; }
+            else                   { f1 = delta / kDc; }
+            real_t f2;
+            if (t < T)             { f2 = 0.0; }
+            else if (t < T + t0)   { f2 = (t - T) / t0; }
+            else                   { f2 = 1.0; }
+            const real_t factor = (f1 > f2) ? f1 : f2;
+            const real_t mu_ref = kMuS + (kMuD - kMuS) * factor;
+            const real_t mu_impl =
+               spatial::LSWFrictionCoefficient_ForcedRupture(
+                  delta, kMuS, kMuD, kDc, t, T, t0);
+            TEST_NEAR(mu_impl, mu_ref, 0.0,
+                      "one-clock mu contract: impl == closed-form spec "
+                      "reference, bit-exact (unified substep dispatch pins "
+                      "the SINGLE composition both seam and interior consume)");
```

Also rename the function to `Test_OneClock_Mu_Composition_Pinned` (or keep the name
but fix the comment block at :142-149 to state the structural argument instead of
implying two paths are compared).

**Test case:** the rewritten test IS the test. It fails if anyone edits the
`max(f1,f2)` composition, the barrier short-circuit, or the ramp branches; it passes
bit-exactly today because the reference mirrors `spatial_friction.hpp:841-892`
op-for-op (verified against the source in this review).

---

### [R-102] MODERATE [drivers/spatial_dyn_driver.cpp:1983-2004] — the R-008 shared-face nucleation warning describes the RETIRED R-1601 behavior and now misinforms every multi-rank run

**Category:** DEVIATION (stale-reference sweep miss — hunt item 6)

**Description:**
The warning block reads: "The wave operator's shared-face EvaluateADER_LSW call
reads DOFData::tau{1,2}_nuc AFTER the Phase N per-substep iterator has accumulated
the full smoothStep increment for the macrostep, so those DOFs see the perturbation
as an end-of-macrostep step rather than a smooth ramp (1st-order time-accuracy
degradation)". That was TRUE under the R-1601 fallback (shared QPs = inline one-shot
over macro dt, after the iterator). Unify-plan Phase 2 retired exactly that: shared
QPs now consume the iterator's per-substep buffer, so they see the SAME per-substep
ramp as interior QPs. The warning now (a) asserts a numerical-accuracy defect that
no longer exists on the production substep path, (b) tells users to tighten dt for
no reason, and (c) hardcodes "the gradual_overstress perturbation" even when
`[nucleation] kind = "forced_rupture"` (observed: this review's np=4 forced-rupture
run printed the gradual_overstress text). This is user-facing misinformation about
the exact defect the plan just dissolved.

**Trigger:**
Any np > 1 run with `nucleation.enabled` and ≥1 fault DOF on a shared face — e.g.
`mpirun -np 4 ./seas_spatial_dyn_driver --config tpv26/configs/tpv26_spatial_forcedrupture_smoke.toml`
(reproduced in this review; meanwhile the np=4 station traces match np=1 to the last
printed digit, directly refuting the warned degradation).

**Actual behavior:**
```
[spatial_dyn] WARNING: 6 of 261 fault DOFs (2.29885%) live on shared faces and will
see the gradual_overstress perturbation as an end-of-macrostep step rather than a
smooth ramp.  Tighten dt (smaller macrostep) to reduce the 1st-order time-accuracy
error on shared faces.
```

**Expected behavior:**
The step-degradation claim applies only to the ADER ONE-SHOT path (no substep buffer
installed). On the substep path (production for LSW/RS/forced-rupture) shared and
interior QPs are per-substep identical (unify plan Phase 2; proven by the parity
oracle and the byte-exactness gate). The message should either be dropped on the
substep path or reworded to an informational shared-DOF count.

**Suggested fix:**

```diff
-   // R-008: warn when a non-trivial fraction of fault DOFs live on
-   // shared faces.  The wave operator's shared-face EvaluateADER_LSW
-   // call reads DOFData::tau{1,2}_nuc AFTER the Phase N per-substep
-   // iterator has accumulated the full smoothStep increment for the
-   // macrostep, so those DOFs see the perturbation as an end-of-
-   // macrostep step rather than a smooth ramp (1st-order time-
-   // accuracy degradation).  Quantify and warn so the user can
-   // tighten dt or accept the trade-off.
+   // R-008 (HISTORY): pre-unify (R-1601 fallback era) shared fault QPs ran an
+   // inline one-shot over the macro dt AFTER the per-substep iterator, so they
+   // saw nucleation perturbations as an end-of-macrostep step.  Unify-plan
+   // Phase 2 (PLAN_unify_interior_shared_fault_substep_2026-07-09.md) retired
+   // that: shared QPs consume the iterator's per-substep buffer and see the
+   // identical per-substep ramp as interior QPs.  Keep an INFORMATIONAL
+   // shared-DOF count (useful for partition diagnostics); the accuracy
+   // warning no longer applies on the substep path.
    if (cfg.nucleation.enabled && num_shared_global > 0 && rank == 0)
    {
       const real_t shared_frac = (num_fault_global > 0)
          ? (static_cast<real_t>(num_shared_global)
             / static_cast<real_t>(num_fault_global))
          : 0.0;
-      std::cout << "[spatial_dyn] WARNING: " << num_shared_global
-                << " of " << num_fault_global << " fault DOFs ("
-                << (100.0 * shared_frac) << "%) live on shared faces "
-                << "and will see the gradual_overstress perturbation as "
-                << "an end-of-macrostep step rather than a smooth ramp.  "
-                << "Tighten dt (smaller macrostep) to reduce the "
-                << "1st-order time-accuracy error on shared faces.\n";
+      std::cout << "[spatial_dyn] INFO: " << num_shared_global
+                << " of " << num_fault_global << " fault DOFs ("
+                << (100.0 * shared_frac) << "%) live on shared (rank-"
+                << "boundary) faces.  Since unify-plan Phase 2 these consume "
+                << "the same per-substep friction buffer as interior DOFs "
+                << "(identical nucleation ramp; no end-of-macrostep step).\n";
    }
```

**Test case:**
```bash
# After the fix: the np=4 forced-rupture smoke must not print "WARNING ... smooth ramp"
mpirun -np 4 ./seas_spatial_dyn_driver --config tpv26/configs/tpv26_spatial_forcedrupture_smoke.toml \
  | grep -c "end-of-macrostep step"   # expect 0
# and the byte-level np-consistency evidence stands:
# stations np=4 vs np=1 differ at most in the last printed digit (verified 2026-07-10).
```

---

### [R-103] MODERATE [dynamic/wave_operator.inl:3956 + plan banner] — Phase-5 acceptance criterion "(INV) asserted inside the single routine" silently dropped

**Category:** DEVIATION

**Description:**
Plan Phase 5 acceptance includes: "The (INV) invariant is asserted inside the single
routine, so it covers both face classes", where (INV) is the geometric orientation
invariant `dot(can_n, centroid_plus − centroid_minus) < 0` (plan :429). The helper
asserts only `dt > 0` and the `dof_idx` range; no (INV) check exists anywhere in the
working tree. The banner documents the frame-bit-parameter deviation but does NOT
mention dropping (INV). There are real mitigating facts — (INV) was defined inside
the WITHDRAWN Phase 1 (its bug hypothesis was falsified by the census), and the
helper has no access to element centroids (asserting there requires new plumbing) —
but an acceptance criterion silently dropped is exactly the pattern this review
exists to catch. It matters doubly here because (INV) is the natural diagnostic for
the OPEN np≥4 attractor whose prime suspect is an interior-vs-shared frame-ORIENTATION
flip relative to the DOFData baked-basis seeding: an (INV)-style orientation check at
both call sites is precisely how that investigation will start.

**Trigger:**
Read the helper (:3956-3963) and grep: `grep -n "centroid" dynamic/wave_operator.inl`
→ no (INV) assertion. Read the banner (plan :62-85) → no mention of the dropped
criterion.

**Actual behavior:**
Criterion unmet, deviation undocumented.

**Expected behavior:**
Either the assertion exists (debug/env-gated, at the call sites where `e1/e2`
respectively `ftr` + face-neighbor data are in scope), or the banner explicitly
records the deferral and re-targets (INV) at the attractor investigation.

**Suggested fix (recommended: option B, documentation — option A is new plumbing
beyond a mechanical fix and touches the hot loop):**

```diff
 # ✅ PHASE 4 + PHASE 5 IMPLEMENTED (2026-07-10, same session)
 ...
   The env-gated `[XRANK]`/`[C-1s]`/`[MACRO]` diagnostics were remapped to the helper's routed
   outputs.  An initial blanket text-replace accidentally hit the RK path
   (`ComputeFaceFluxRHS`) — caught by the compiler, reverted verbatim; RK remains untouched.
+- **Dropped criterion (documented deviation):** the Phase-5 acceptance line "the (INV)
+  invariant is asserted inside the single routine" is NOT implemented.  (INV)
+  (`dot(can_n, centroid_plus − centroid_minus) < 0`) was specified inside the WITHDRAWN
+  Phase 1; the helper has no element-centroid access, so asserting it there needs new
+  plumbing.  (INV) is deliberately DEFERRED to the np≥4 attractor investigation, where a
+  call-site orientation census (which HAS e1/e2 / face-neighbor context) is the first
+  planned probe — the attractor's prime suspect is exactly an interior-vs-shared frame-
+  orientation flip that an (INV)-style check would expose.
```

**Test case:**
```bash
# Documentation fix: the banner must acknowledge the dropped criterion.
grep -c "INV" miniapps/seas/document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md
# The Phase 4+5 banner section must contain an (INV) deferral entry (currently 0 mentions
# in the banner; the string appears only in the withdrawn-Phase-1/Phase-5 sections below).
```
(If option A is chosen instead: on the 2-tet fixture, tampering the frame bit via
`TamperSharedFaultElem1OnPlus`-style test hooks must abort the orientation check —
mirror of the Phase-1 acceptance sketch. That is a new parallel test, not a
mechanical fix.)

---

### [R-104] LOW [dynamic/fault_face_flux.hpp:452, dynamic/fault_face_flux.cpp:853] — present-tense references to the retired R-1601 fallback as the rationale for the slip-ownership rule

**Category:** QUALITY (stale reference — hunt item 6)

**Description:**
Both comments justify the (still correct and load-bearing) rule "EvaluateADER_LSW
must not accumulate slip" with: "The wave operator's R-1601 shared-fault fallback
re-invokes EvaluateADER_LSW AFTER the iterator on shared QPs; if this function also
accumulated slip, shared-fault QPs at np > 1 would double-count". The fallback no
longer exists (Phase 2); the inline dispatch now runs only when the substep buffer is
NOT installed, i.e. when the iterator did not run — the described double-count
scenario is unreachable. A maintainer who greps for the fallback will find nothing
and may conclude the rule is dead and delete it. The RULE must stay (one-shot-path
correctness + division of ownership); only the RATIONALE text is stale.

**Suggested fix** (same rewording at both sites):

```diff
-   /// R-001 (final review): `data.slip{1,2}` is NOT touched by this
-   /// call — slip evolution is owned exclusively by
-   /// `Tpv205SubStepIterator::StepOneQP_`, which integrates `slip{1,2}
-   /// += V{1,2} * dt_sub` once per sub-step over every fault QP.  The
-   /// wave operator's R-1601 shared-fault fallback re-invokes
-   /// `EvaluateADER_LSW` AFTER the iterator on shared QPs; if this
-   /// function also accumulated slip, shared-fault QPs at np > 1 would
-   /// double-count slip and the rupture front would accelerate
-   /// artificially across MPI rank boundaries.  See REVIEW R-001.
+   /// R-001 (final review): `data.slip{1,2}` is NOT touched by this
+   /// call — slip evolution is owned exclusively by
+   /// `Tpv205SubStepIterator::StepOneQP_`, which integrates `slip{1,2}
+   /// += V{1,2} * dt_sub` once per sub-step over every fault QP.
+   /// HISTORY: under the retired R-1601 shared-fault fallback (superseded
+   /// by unify-plan Phase 2, see R1601_root_cause_2026-07-10.md) this
+   /// function ran AFTER the iterator on shared QPs, so accumulating slip
+   /// here double-counted at np > 1.  The rule REMAINS load-bearing on the
+   /// one-shot (no-buffer) path and preserves the iterator's exclusive
+   /// ownership of slip.  See REVIEW R-001.
```

(No test case — comment-only change; LOW.)

---

### [R-105] LOW [plan banner :64-75 + Phase-5 acceptance] — banner enumeration of remaining frame-build sites is wrong; the literal grep acceptance metric is unmet/meaningless

**Category:** QUALITY (doc accuracy)

**Description:**
Two inaccuracies that will mislead the (open) np≥4 attractor investigation, which
must reason precisely about "which code builds canonical frames":
1. Banner: "Verified no GetTime() ... the two physics uses live only in the inline
   one-shot arms" — post-Phase-5 there is exactly ONE physics use (the helper's
   one-shot arm, `wave_operator.inl:4033`).
2. Banner: "`should_negate_frame` expressions in the corrector: 6 → 2 (the untouched
   RK pair)" and plan acceptance "grep -c 'should_negate_frame' drops from 6 to ≤ 2".
   Actual remaining frame-build decision sites: producer
   `EvaluateBulkAtFaultQPsCanonical` interior loop (:2221, `!elem1_on_plus`), RK
   interior `ComputeFaceFluxRHS` (:2972, `!elem1_on_plus`), RK shared
   `ComputeSharedFaceFluxRHS` (:3646, `qpd.sign_flipped` inline — never spelled
   `should_negate_frame`), stage-avg `_qq` block (:4604, `!elem1_on_plus`), and the
   helper itself (:3967). A literal `grep -c` counts 19 lines (HEAD: 16) — the
   metric as written INCREASED and never measured what was meant. What Phase 5
   actually delivered (and what should be stated): the ADER corrector's two NODAL
   fault branches now share ONE frame-build/dispatch site (the helper).

**Suggested fix:** correct the banner text to the enumeration above and restate the
acceptance line as "the ADER corrector nodal fault branches contain zero frame-build
sites outside `ComputeFaultQPImposedStatesCanonical_` (remaining sites: producer
:2221, RK pair :2972/:3646, stage-avg :4604 — all outside Phase-5 scope)". (No test
case — doc-only; LOW.)

---

### [R-106] LOW [dynamic/wave_operator.hpp:599] — private-convention helper (`..._` suffix) declared in the public section

**Category:** QUALITY

**Description:**
`ComputeFaultQPImposedStatesCanonical_` carries the trailing-underscore private
naming convention and mutates `fdata` / dispatches physics, but is declared in the
`public:` section (last access specifier before :599 is `public:` at :126; the
class's `private:` starts at :1287). Nothing outside the class calls it (both call
sites are member functions); leaving it public invites external callers to bypass
the fault-branch guards (dof-offset lookup, R-101 flag check) that the call sites
perform.

**Suggested fix:** move the declaration (with its doc block) below the `private:`
specifier at :1287 — zero behavioral change, both call sites are members. Verify
with a clean rebuild.

(No test case — access-control move; LOW. If any test calls it directly under
`SEAS_TEST_INTERNAL`, keep it public and instead add a doc line "private by
convention; public only for test access" — grep found no such caller.)

---

## Summary

- Critical issues: **0**
- Moderate issues: **3** (R-101 vacuous test, R-102 stale user-facing warning
  contradicting Phase 2, R-103 silently dropped (INV) acceptance criterion)
- Low issues: **3** (R-104 stale R-1601 rationale ×2 sites, R-105 banner/metric
  accuracy, R-106 access specifier)
- Plan compliance: **PARTIAL** — Phase 4 items (a)-(e) all delivered and now
  evidence-backed (this review supplied the missing np-multi forced-rupture seam
  run); Phase 5's pure-refactor contract is PROVEN (RK byte-identity, op-exact
  extraction, full bit-exactness gate + review reruns), with two documented
  deviations (frame-bit parameter — justified by the census) and one UNdocumented
  deviation ((INV) dropped, R-103).
- Verdict: **PASS WITH FIXES** — no correctness defect found in the shipped
  numerics; all three MODERATEs are guard/communication defects (a test that cannot
  fail, a warning that asserts a dissolved defect, an unrecorded dropped criterion).
  Fix before closing the goal.

## Fix Status (applied 2026-07-10, same session)

| ID | Severity | Status | What was done |
|---|---|---|---|
| R-101 | MODERATE | **FIXED** | Test rewritten as `Test_OneClock_Mu_Composition_Pinned`: closed-form spec reference (f1/f2/max composition mirrored op-for-op) compared at tol 0.0 across the same 40-point grid; comment now states the structural seam==interior argument. 68/68. |
| R-102 | MODERATE | **FIXED** | Driver message WARNING→INFO per the suggested diff; stale end-of-macrostep/tighten-dt claim removed, HISTORY comment cites unify plan Phase 2 + this review. Verified: np=4 smoke prints the new INFO line and zero occurrences of the old text; V_max unchanged. |
| R-103 | MODERATE | **FIXED (option B)** | (INV) dropped-criterion deferral entry added to the plan's Phase 4+5 banner, re-targeting (INV) at the np≥4 attractor investigation (call-site orientation census). |
| R-104 | LOW | **FIXED** | Both `fault_face_flux.{hpp,cpp}` comments reworded to HISTORY (retired R-1601 fallback, root-cause doc cited); the slip-ownership rule and its one-shot-path justification retained. |
| R-105 | LOW | **FIXED** | Banner corrected: ONE physics `GetTime()` (helper one-shot arm :4033); frame-build site enumeration replaces the ill-posed "6 → 2" grep metric. |
| R-106 | LOW | **FIXED** | Helper declaration + doc moved below `private:` (next to `ComputeADERFaceFluxRHS`), with a "callers must perform the fault-branch guards" note. Clean rebuild. |

**Post-fix verification (all green):** driver + both parity tests rebuilt
(`SEAS_FS_LIB=""`, worktree MFEM overrides); forced-rupture parity 68/68; substep
parity 39/39; forced-rupture smoke np=1 V_max_global = 0.199649 and np=4 = 0.199649
(6 shared fault DOFs, zero reconcile aborts, stale warning text absent, new INFO line
present); TPV26 smoke np=1 V_max_global = 0.185543; all 32 np=1/np=4 station files
**byte-identical** to the review-run references (the fixes are output-neutral except
the reworded rank-0 INFO line).

## Unreviewed Areas

- The RK pair's internals beyond byte-identity to HEAD (unchanged code, out of
  Phase-5 scope by plan decision).
- `ComputeADERFaceFluxRHS_CachedInterior_` (:804) and the `--face-cache` path —
  untouched by this change set (verified no diff vs HEAD in that region).
- The stage-avg `_qq` branch's numerics (:4550-4720) — unchanged from HEAD;
  only its existence/scope was assessed (R-105).
- Phase 6 (gold regeneration) and the np≥4 attractor — explicitly out of scope,
  tracked in the plan's Remaining section.
- Makefile deltas — belong to Phases 0-3 (already reviewed in
  `REVIEW_phase2_3_unify_2026-07-10.md`); Phases 4+5 required no Makefile change
  (confirmed: the FR-parity and substep-parity targets rebuilt and ran green).
