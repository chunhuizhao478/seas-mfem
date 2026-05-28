# Code Review: Phase 5 — unified method-oriented sub-step iterator — 2026-05-27

> NOTE: the pre-existing `REVIEW.md` (TPV205, 2026-04-27) and
> `REVIEW_phase4_drename_2026-05-27.md` are unrelated and left intact. Point
> `/code-fix` at THIS file for the Phase-5 findings.

## Review Scope
- **Plan:** `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`, §Phase 5 (incl. the 2026-05-27 ⚠️ deviation block).
- **Worktree:** `/Users/chunhuizhao/projects/seas-mfem-safs` (branch `safs`). The `/code-fix` agent must operate here.
- **Files reviewed:** `friction/state_policies.hpp`, `dynamic/friction_substep_iterator.{hpp,cpp}`, `dynamic/friction_iterator_factory.cpp`, `tests/unit/test_friction_substep_iterator_parity.cpp`, `tests/unit/test_friction_iterator_factory.cpp`, `Makefile` (Phase-5 additions).
- **Domain context:** root + miniapps/seas `CLAUDE.md`; the three oracle iterators `dynamic/tpv{205,102,104}_substep_iterator.*`; the implementer completion report.

## Independent verification performed
- `RunSubSteps_` confirmed to index `Qp_per_substep[o]`/`Qm_per_substep[o]`, `time_weights_[o]`, `deltaT_[o]`, and to set `t_sub_cursor = t_macro_start` then `t_sub_end = t_sub_cursor + dt_sub` — **correct** (no indexing bug).
- LSW `StepOneQP_` is a faithful line-by-line lift of `tpv205_substep_iterator.cpp:75-142`; the RS `step_fn` matches the tpv102/tpv104 per-QP body (ComputeStageState → slip → ψ(policy) → BuildImposedState → WriteBack(last)).
- `SEAS_HEADERS → DYNAMIC_HEADERS` includes `fault_face_flux.hpp`, so the new objects rebuild on a `DOFData` change (no stale-ABI risk).
- Test outcomes (this session): `seas_test_friction_substep_iterator_parity` **27/27**; `seas_test_friction_iterator_factory` **10/10**; `seas_spatial_dyn_driver` builds EXIT=0.

**No CRITICAL product bug found.** The unified iterators reproduce the standalone oracles bit-for-bit on the tested fixtures. The findings below are test-robustness gaps that weaken the "bit-for-bit" acceptance guarantee, plus two low-severity items.

## Findings

### [R-001] MODERATE [tests/unit/test_friction_substep_iterator_parity.cpp:MakeQField/MakeQuadrature/main] — parity fixture is degenerate (constant Q per sub-step, uniform weights, t0=0) → cannot catch sub-step-indexing or t_macro_start regressions

**Category:** EDGE_CASE (test coverage)

**Description:**
The parity fixture feeds the **same** `Q` to every sub-step (`Qp.assign(O, qp)`), **uniform** weights (`weights.assign(O, 1.0/O)`), and `t_macro_start = 0.0`. With these degenerate inputs the unified and standalone iterators produce identical output **even if** a unified-iterator regression read `Qp[0]` instead of `Qp[o]`, used `time_weights_[0]` instead of `[o]`, or dropped the `t_macro_start` offset — because all sub-steps see identical `Q`/weights and the nucleation time base is 0. The headline acceptance criterion is "reproduces … **bit-for-bit**"; the current test verifies that only on inputs that can't distinguish the per-sub-step paths.

**Trigger:** A future edit to `RunSubSteps_` that mis-indexes the sub-step (`[0]` vs `[o]`) or omits `t_macro_start`.

**Actual behavior:** Test passes regardless (degenerate inputs mask the difference).

**Expected behavior:** Test fails when the unified iterator diverges from the standalone on per-sub-step-varying inputs.

**Suggested fix:** vary `Q` per sub-step, use non-uniform weights, and a non-zero `t0`.
```diff
-void MakeQField(int n, int O, std::vector<std::vector<real_t>> &Qp,
-                std::vector<std::vector<real_t>> &Qm)
-{
-   std::vector<real_t> qp(static_cast<size_t>(NUM_STATE) * n, 0.0);
-   std::vector<real_t> qm(static_cast<size_t>(NUM_STATE) * n, 0.0);
-   for (int i = 0; i < n; ++i)
-   {
-      qp[static_cast<size_t>(i) * NUM_STATE + VZ] = +1.0e-3 * (1.0 + 0.1 * i);
-      qm[static_cast<size_t>(i) * NUM_STATE + VZ] = -1.0e-3 * (1.0 + 0.1 * i);
-   }
-   Qp.assign(O, qp);
-   Qm.assign(O, qm);
-}
+void MakeQField(int n, int O, std::vector<std::vector<real_t>> &Qp,
+                std::vector<std::vector<real_t>> &Qm)
+{
+   Qp.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
+   Qm.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
+   for (int o = 0; o < O; ++o)
+      for (int i = 0; i < n; ++i)
+      {
+         // Per-sub-step AND per-QP variation so [o] vs [0] indexing matters.
+         const real_t s = 1.0e-3 * (1.0 + 0.1 * i) * (1.0 + 0.3 * o);
+         Qp[o][static_cast<size_t>(i) * NUM_STATE + VZ] = +s;
+         Qm[o][static_cast<size_t>(i) * NUM_STATE + VZ] = -s;
+      }
+}
```
```diff
-void MakeQuadrature(int O, real_t dt_macro,
-                    std::vector<real_t> &deltaT, std::vector<real_t> &weights)
-{
-   deltaT.assign(O, dt_macro / static_cast<real_t>(O));
-   weights.assign(O, 1.0 / static_cast<real_t>(O));
-}
+void MakeQuadrature(int O, real_t dt_macro,
+                    std::vector<real_t> &deltaT, std::vector<real_t> &weights)
+{
+   // Non-uniform sub-steps + weights (still Σδt = dt_macro, Σw = 1) so a
+   // [0]-vs-[o] weight/deltaT mis-index is observable.
+   deltaT.assign(O, 0.0); weights.assign(O, 0.0);
+   real_t dsum = 0.0, wsum = 0.0;
+   for (int o = 0; o < O; ++o)
+   { deltaT[o] = (1.0 + 0.5 * o); weights[o] = (1.0 + 0.25 * o);
+     dsum += deltaT[o]; wsum += weights[o]; }
+   for (int o = 0; o < O; ++o)
+   { deltaT[o] = deltaT[o] / dsum * dt_macro; weights[o] /= wsum; }
+}
```
```diff
-   const real_t t0       = 0.0;
+   const real_t t0       = 3.7;   // non-zero so t_macro_start propagation is tested
```

**Test case:**
```python
def test_R001_per_substep_varying_Q_still_bit_parity():
    # With Q varying per sub-step, non-uniform weights, t0 != 0:
    # the unified iterator must STILL match the standalone bit-for-bit
    # (max|diff| == 0) for LSW/aging/SRW at O in {1,2,3}.
    # Pre-fix: the degenerate fixture would pass even with a [o]->[0] bug.
    assert parity_max_diff(O=2, vary_substep=True, nonuniform_w=True, t0=3.7) == 0.0
```

---

### [R-002] MODERATE [tests/unit/test_friction_substep_iterator_parity.cpp:main] — parity is not guarded against vacuous (degenerate) output

**Category:** EDGE_CASE (test coverage)

**Description:**
Every assertion is `std == uni` (`MaxDofVecDiff`/`MaxVecDiff == 0`). Nothing asserts the friction solve produced **non-trivial** output. If a fixture change drove `V_abs ≈ 0` everywhere (locked fault, no slip), the comparison would pass on all-zero / unchanged DOFData without exercising `ComputeStageState`/`SolveLSW_TPV205`/the ψ-update at all — a green test that proves nothing.

**Trigger:** A fixture/material change that yields zero slip (e.g., trial traction below frictional strength on every QP).

**Actual behavior:** Bit-parity passes vacuously; the friction kernels are never meaningfully exercised.

**Expected behavior:** The test asserts the standalone produced a meaningful state change (non-zero slip/velocity, ψ moved) before trusting the bit-parity.

**Suggested fix:** add a non-triviality gate per iterator on the standalone result.
```diff
         TEST_EQ0(MaxDofVecDiff(dof_std, dof_uni), "LSW DOFData bit-parity");
         TEST_EQ0(MaxVecDiff(Ip_std, Ip_uni),      "LSW I_imp_plus bit-parity");
         TEST_EQ0(MaxVecDiff(Im_std, Im_uni),      "LSW I_imp_minus bit-parity");
+        // Guard against vacuous parity: the solve must have moved the state.
+        TEST_ASSERT_NONTRIVIAL(std::abs(dof_std[0].V2) > 0.0 ||
+                               std::abs(dof_std[0].slip2) > 0.0,
+                               "LSW fixture exercises a non-trivial friction solve");
```
(add an analogous `TEST_ASSERT_NONTRIVIAL` for aging and SRW, and a small
`#define TEST_ASSERT_NONTRIVIAL(cond,msg)` that increments counters like `TEST_EQ0`).

**Test case:**
```python
def test_R002_fixture_is_nontrivial():
    # The standalone LSW/aging/SRW runs must change V2 or slip2 from the
    # fixture's initial 0; otherwise the bit-parity assertions are vacuous.
    assert standalone_lsw_run().V2 != 0.0 or standalone_lsw_run().slip2 != 0.0
```

---

### [R-003] LOW [dynamic/friction_iterator.hpp] — Phase-2 adapters are now dead/divergent code (factory no longer returns them)

**Category:** QUALITY

**Description:**
After re-pointing `MakeFrictionIterator` to `LinearSlipWeakeningIterator` / `RateStateAgingIterator`, the Phase-2 adapters `LswFrictionIterator` and `RateStateAgingFrictionIterator` (in `friction_iterator.hpp`) are no longer constructed by the factory or referenced by the factory test. They still exist and still wrap the **standalone** iterators — so the codebase now has two LSW paths (`LswFrictionIterator → Tpv205`, `LinearSlipWeakeningIterator → unified`) and two aging paths. A future caller that picks `LswFrictionIterator` would silently get the standalone path, not the production one.

**Trigger:** A new caller (or test) constructing `LswFrictionIterator`/`RateStateAgingFrictionIterator` expecting the production path.

**Actual behavior:** Two parallel implementations of the same law coexist with no compile-time signal that one is retired.

**Expected behavior:** Either remove the now-unused adapters, or mark them `[[deprecated]]` with a comment pointing to the unified classes.

**Suggested fix (minimal — deprecate rather than delete, to avoid breaking `test_advance_interface_compiles` if it references them):**
```diff
-class LswFrictionIterator : public IFrictionIterator
+/// DEPRECATED (Phase 5): the factory now returns LinearSlipWeakeningIterator
+/// (unified). Retained only for the Phase-2 compile test; do not use in new code.
+class LswFrictionIterator : public IFrictionIterator
```
(and the same on `RateStateAgingFrictionIterator`). Confirm no remaining
non-test references with `grep -rn 'LswFrictionIterator\|RateStateAgingFrictionIterator' --include='*.cpp'`.

---

### [R-004] LOW [dynamic/friction_substep_iterator.cpp:LinearSlipWeakeningIterator::Advance] — env-gated [SLIP] trace is untested and carries an unguarded division

**Category:** QUALITY

**Description:**
The `[SLIP]` trace lifted from tpv205 computes `sn_sterm = d.eta_p * (Qp_i[SXX]/d.Zp_plus + Qm_i[SXX]/d.Zp_minus)`. It only runs under `SEAS_DIAG_SLIP` + `V_abs > thr`, so the parity test (env unset) never exercises it, and a `d.Zp_plus == 0` fixture would divide by zero. This is identical to the oracle (not a new defect) and diagnostic-only, but it is now duplicated in a second location and untested here.

**Trigger:** Running the unified LSW path with `SEAS_DIAG_SLIP=1` and a spiking QP where `d.Zp_plus == 0`.

**Actual behavior:** Untested; would divide by zero on a degenerate impedance (same as the oracle).

**Expected behavior:** Either covered by an env-on smoke, or annotated that `d.Zp_plus/Zp_minus > 0` is a precondition (the impedances are always positive in practice).

**Suggested fix:** add a brief comment documenting the `Zp_plus/Zp_minus > 0` precondition at the trace site (no behavior change), or a follow-up env-on diagnostic test. No production-path change required.

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 degenerate parity inputs; R-002 vacuous-parity guard)
- Low issues: 2 (R-003 dead Phase-2 adapters; R-004 untested/unguarded LSW trace)
- **Plan compliance:** FULL for the product code (unified `SubStepIteratorBase`/`RateStateSubStepIterator<Policy>`/`LinearSlipWeakeningIterator`, both policies, factory re-point, all documented deviations sound). PARTIAL on the *acceptance evidence*: the parity test passes but under-verifies the "bit-for-bit" claim (R-001/R-002).
- **Verdict:** PASS WITH FIXES — no source-code bug blocks proceeding; strengthen the parity test (R-001, R-002) so the bit-for-bit guarantee is genuinely enforced, then optionally R-003/R-004 cleanup.

## Unreviewed Areas
- **SAFS+RS / LSW full MPI smoke** (Phase-5 acceptance criterion #2: "bit-identical to Phase-3 results") — not run (heavyweight; local MPI bus-errors). The 27/27 iterator-level bit-parity + factory test + driver build are strong proxies, but the end-to-end smoke is unverified.
- **Factory SRW dispatch** — intentionally gated until Phase 6 (`state_evolution` selector + `rs->V_w`); the `RateStateSlipLawSrwIterator` is exercised only by the parity test, not by any production config yet.
- **Non-`safs` worktree** — the rename/Phase-5 work lives on branch `safs`; `system/spatial_dyn_driver` is unaffected.
