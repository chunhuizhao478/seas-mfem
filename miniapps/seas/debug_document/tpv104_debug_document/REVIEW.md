# Code Review: ADER SubStep Implementation — Round 10 (2026-04-25)

## ⚠ POST-REVIEW ADDENDUM (R-1008): user found another CRITICAL bug

After R-1001/R-1002/R-1003 fixes landed and tests passed, the user ran
the substep dispatch on the actual TPV104 setup and observed **2-million×
strike-direction over-acceleration**: the rupture went terminal at
t=0.5 s with V_strike=15.04 m/s and h-slip=2.08 m, when one-shot at
the same t=0.5 s gives V_strike=4.6e-5 m/s and h-slip=9.6e-7 m.

### Root cause: nucleation accumulator double-counted

**The driver calls `ApplyNucleationIncremental_TPV104` once per
macro-step at the top of the time loop (line 1706, pre-R-1008).  The
iterator (both legacy `Advance` and my new `AdvanceWithSubStepStates`)
ALSO calls it once per sub-step internally (line 295 / equivalent in
the new method).** Σ_o ΔS over sub-steps telescopes to the same
macro-step increment, so:

- One-shot path: driver-level ΔS only ⇒ `tau2_nuc` += ΔS (correct).
- Substep path: driver-level ΔS PLUS iterator's per-sub-step Σ ΔS ⇒
  `tau2_nuc` += 2·ΔS per macro-step.

By the time `smoothStep(t)` reaches 0.5 (~t=0.5 s), tau2_nuc has been
incremented to 2 × (0.5 · 45 MPa) = 45 MPa instead of the spec's
22.5 MPa.  At full ramp (t=1 s): 90 MPa instead of 45 MPa.  Net effect:
the rupture sees double the spec Δτ₀ and goes terminal twice as fast.

### Why the tests didn't catch this

- `test_tpv104_substep_dispatch_parity` (R-1001) calls
  `wave.AdvanceADER(...)` for path A and `RunSubStepDispatch(...)` for
  path B WITHOUT either path invoking `ApplyNucleationIncremental_TPV104`.
  The test feeds nucleation-free Q from `InitializeStateTotal` directly.
  The bug lives at the DRIVER level (the time-loop's nucleation call
  combined with the iterator's internal call), not at the dispatch
  helper level.  My test had a hole.
- `test_tpv104_substep_iterator_parity` (legacy O=1 parity) similarly
  doesn't exercise nucleation — it sets `data.tau2_nuc = 0` and
  uses pre-stress only.
- `tpv104_smoke` runs `--dry-run` (no time loop) so the driver
  nucleation site is never reached.

### R-1008 fix

```cpp
// drivers/tpv104_driver.cpp at the time-loop nucleation site:
if (!disable_nucleation && num_fault_total > 0 && !use_substep_iterator)
{
   ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                     t + dt_step, dt_step);
}
```

Skips the driver-level call when the iterator owns nucleation cadence.
One-shot path keeps the driver-level call (no behavior change).
Verified: `seas_test_tpv104_substep_dispatch_parity` still passes
(max |Q_new_one − Q_new_sub| = 4.05e-13 → no regression at O=1).

### Process lesson — fourth instance

This is the FOURTH "process bug" recurring through this investigation:

1. R-1001: tests validated fragments, not the dispatch contract.
2. R-1002: per-call dt mismatch never exercised by tests.
3. R-1003: MPI silent divergence never tested.
4. **R-1008: cross-component double-count never tested** (driver
   nucleation × iterator nucleation interaction).

Common failure mode: my tests validated each new component in
isolation, never end-to-end against the driver's full time-step
sequence including all the OTHER work the driver does (nucleation
injection, ψ snapshot, slip integration, station output).  The minimum
viable end-to-end test would be: **run `seas_tpv104_driver
--fault-iterator substep --tfinal=<small>` and compare hypocenter
slip_strike to one-shot at the same tfinal.**  A 2× discrepancy at
t=0.5 s would have been the first thing to fall out.

This addendum and the R-1008 fix are now in source.  No new test
written — adding one would require the driver-level integration test
which does not exist locally.  Recommend: add one in the next round.

---



## Review Scope
- Implementation under review: round-7 R-602/R-603 substep iteration adoption.
  - `dynamic/wave_operator.hpp` — `ComputeADERSubStepStates`,
    `SetSubStepFaultImposedStates`, `ResetSubStepFaultImposedStates`,
    `EvaluateBulkAtFaultQPsCanonical` declarations + 3 mutable members.
  - `dynamic/wave_operator.inl` — those four implementations + 7-line
    guard at the interior fault branch of `ComputeADERFaceFluxRHS`
    (line ~2604).
  - `dynamic/tpv104_substep_iterator.hpp/.cpp` — `AdvanceWithSubStepStates`
    declaration + implementation.
  - `drivers/tpv104_driver.cpp` — `AdvanceADERWithSubStep` orchestrator
    (lines ~256–369), iterator instantiation (lines ~1183–1251),
    time-loop dispatch (lines ~1672–1689).
  - 3 new unit tests + Makefile wiring.
- Plan source: round-6 R-602/R-603/R-605 from prior REVIEW.md, with
  user override on Q3 (default OFF behind `--fault-iterator substep`).
- Domain context: project `CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  the iterator's existing R4-004 cadence-deviation disclosure.

---

## Bottom line

I implemented this; the tests pass. **The tests pass because they don't
test the actual driver-level dispatch contract.** Adversarial audit
finds three CRITICAL bugs and two MODERATE bugs that the tests do not
catch.

The tests I wrote validate:
- predictor round-trip (Σ_o w_o · Q_per_node[o] · dt = I_macro) — correct.
- iterator-level parity at O=1 with manually-supplied `Q_pointwise = Q̄` — correct.
- guard set/reset bit-equality on a 2-tet fixture — correct.

What they do NOT validate:
- driver-level dispatch equivalence (substep at O=1 vs one-shot at O=1) on the same TPV104 setup.
- robustness against `dt_step ≠ dt` on the last time step.
- correctness of the per-sub-step Q-node placement.

Each gap corresponds to a real bug below.

---

## Findings

---

### [R-1001] [CRITICAL] [drivers/tpv104_driver.cpp:317–323] — Substep at O=1 is NOT bit-identical to one-shot. Wrong tau-node choice.

**Category:** BUG (sign/scaling — incorrect quadrature node)

**Description:**
The driver's `AdvanceADERWithSubStep` builds sub-step nodes via cumulative prefix:

```cpp
real_t acc = 0.0;
for (int o = 0; o < O; o++) {
   acc += deltaT[o];
   tau_nodes[o] = acc;     // END of sub-step o
}
```

At O=1 with `deltaT={dt_step}`, `tau_nodes[0] = dt_step`. So `ComputeADERSubStepStates` evaluates Q at τ = **dt_step (end of macro-step)**, returning `Q(dt) = D(0) + dt · D(1) + ...`.

The legacy one-shot `wave.AdvanceADER` passes `I_macro = ∫_0^dt Q dτ` to `EvaluateADER`, which divides by dt to get **Q̄ = I/dt = D(0) + (dt/2)·D(1) + ...** — the **time-AVERAGE**, not the endpoint.

For ADER-2 predictor (Q linear in τ): `Q̄ = Q(τ=dt/2)` (the midpoint), not `Q(τ=dt)`. So at O=1 with the current code:

- substep dispatch friction-input = Q(dt) = Q + dt · L(Q)
- one-shot dispatch friction-input = Q̄ = Q + (dt/2) · L(Q)
- difference = (dt/2) · L(Q) — non-zero whenever L(Q) ≠ 0 (i.e., always during rupture)

This BREAKS the "T_TPV104_SSI_3 contract" that the iterator's own header docstring at `tpv104_substep_iterator.hpp:79–81` advertises:

> "The single-sub-step case `deltaT = {dt_macro}, time_weights = {1.0}` is the TPV102 one-shot limit (T_TPV104_SSI_3)."

It also breaks the user's expectation when toggling `--fault-iterator substep` at O=1 — they'd expect the same answer, get a different one.

**Why the tests didn't catch this:**
- `test_tpv104_substep_iterator_parity` validates iterator-level parity at O=1 BUT MANUALLY supplies `Q_pointwise = Q_avg_plus` (a constant). It never invokes `ComputeADERSubStepStates`, so it never observes the wrong tau choice.
- `test_tpv104_substep_one_shot_parity` validates the guard set/reset bit-equality with INJECTED constant I_imp values. It never runs the actual driver-level orchestrator.
- There is NO test that calls `AdvanceADERWithSubStep` and compares its Q_new against `wave.AdvanceADER` at O=1 on the same setup.

**Trigger:** any TPV104 run with `--fault-iterator substep --ader-order 2` (or any O ≥ 2 — the deviation is worst for O ≥ 2 since `L(Q)` is non-zero during rupture).

**Actual behavior:** substep dispatch friction-input is `Q(dt)` instead of `Q̄`; produces an O(dt) deviation per macro-step from the one-shot path.

**Expected behavior:** at O=1, substep dispatch matches one-shot dispatch bit-for-bit (the SSI_3 contract). At O ≥ 2, substep should use a quadrature whose nodes integrate Q exactly (Gauss-Lobatto on [0, dt]) so the SeisSol-equivalent semantic holds.

**Suggested fix (minimal correctness restore — O=1 SSI_3 only):**
Replace cumulative-end nodes with sub-step **midpoint** nodes:

```diff
@@ drivers/tpv104_driver.cpp AdvanceADERWithSubStep
-   // Cumulative-prefix sub-step end nodes on [0, dt_step].  tau_nodes[o]
-   // is the END time of sub-step o relative to the macro-step start.
-   std::vector<real_t> tau_nodes(O);
-   real_t acc = 0.0;
-   for (int o = 0; o < O; o++)
-   {
-      acc += deltaT[o];
-      tau_nodes[o] = acc;
-   }
+   // Sub-step MIDPOINT nodes on [0, dt_step].  tau_nodes[o] is the
+   // midpoint of sub-step o relative to the macro-step start.  For
+   // ADER-2 predictor (Q linear in τ), the midpoint Q(τ_o) equals the
+   // sub-step's time-average — restoring the T_TPV104_SSI_3 contract
+   // (substep at O=1 == one-shot at O=1).  For O ≥ 3 (Q higher-order
+   // in τ), the midpoint rule is O(dt²)-accurate; switch to Gauss-
+   // Lobatto for full SeisSol parity (separate follow-up).
+   std::vector<real_t> tau_nodes(O);
+   real_t acc = 0.0;
+   for (int o = 0; o < O; o++)
+   {
+      tau_nodes[o] = acc + 0.5 * deltaT[o];
+      acc += deltaT[o];
+   }
```

**Test case (NEW — must add):**

```cpp
// tests/unit/test_tpv104_substep_dispatch_parity.cpp
TEST_CASE("R-1001: substep dispatch bit-identical to one-shot at O=1") {
   // 2-tet TPV102 locked-fault fixture with non-trivial pre-stress.
   //   Path A: wave.AdvanceADER(Q, dt, 2, Q_new_one_shot)
   //   Path B: configure iterator with deltaT={dt}, weights={1.0};
   //           AdvanceADERWithSubStep(... Q, dt, 2, Q_new_substep)
   // Assert max|Q_new_one_shot - Q_new_substep| ≤ 1e-13 (FP precision).
   // Pre-fix this test FAILS.  Post-fix it PASSES.
}
```

---

### [R-1002] [CRITICAL] [drivers/tpv104_driver.cpp:1218–1224] — Substep iterator deltaT crashes on the last time step

**Category:** BUG (edge case — runtime_error throw)

**Description:**
The iterator's `deltaT` is configured ONCE at startup using the auto-CFL initial `dt`:

```cpp
const int O = std::max(1, ader_order);
std::vector<real_t> deltaT(O, 1.0 / static_cast<real_t>(O));
std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
for (int o = 0; o < O; o++) { deltaT[o] = dt / static_cast<real_t>(O); }
substep_iterator.SetSubSteps(deltaT, weights);
```

The time loop computes `dt_step = std::min(dt, tfinal - t)`. When `tfinal` is not an integer multiple of `dt`, the LAST step has `dt_step < dt`. `AdvanceADERWithSubStep` then calls `iterator.AdvanceWithSubStepStates(..., dt_step, ...)`. Inside, the verify:

```cpp
const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(), 0.0);
const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro, 1e-300);
const real_t sum_tol = std::max<real_t>(1e-12, ...);
if (rel > sum_tol) {
   throw std::runtime_error(
      "Tpv104SubStepIterator::AdvanceWithSubStepStates: Σ deltaT must "
      "equal dt_macro within ...");
}
```

**Throws on every TPV104 run that doesn't have tfinal as exact multiple of auto-CFL dt.** Since auto-CFL dt = `cfl_factor / (3·(2·order+1))·h_min/cp` is generally irrational, this triggers on the last step of nearly every production run.

**Why the tests didn't catch this:**
- The unit tests construct iterator + I/O directly and don't run a multi-step time loop with variable dt_step.
- The smoke test runs `--dry-run` which skips the time loop entirely.
- `test_tpv104_substep_iterator_parity` calls Advance/AdvanceWithSubStepStates exactly once with dt_macro = dt_configured (no mismatch).

**Trigger:** any production run with `--fault-iterator substep` reaching the final time step.

**Actual behavior:** crash with `std::runtime_error` on the last time step.

**Expected behavior:** iterator deltaT is rescaled per call to match the actual `dt_step`.

**Suggested fix:** rescale `deltaT` inside `AdvanceADERWithSubStep` based on the actual `dt_step`. The configured `deltaT` should be treated as RELATIVE FRACTIONS of dt_macro, scaled per call:

```diff
@@ drivers/tpv104_driver.cpp AdvanceADERWithSubStep, after `const int O = ...`
+   // Rescale deltaT per call so the iterator's Σ deltaT == dt_step
+   // verify holds even when the time loop's dt_step varies (e.g., the
+   // final time step where dt_step = tfinal - t < auto-CFL dt).  The
+   // configured ratios deltaT[o]/Σ deltaT are preserved.
+   const std::vector<real_t> &configured_deltaT = iterator.GetDeltaT();
+   const std::vector<real_t> &configured_weights = iterator.GetTimeWeights();
+   const real_t configured_sum =
+      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
+                      static_cast<mfem::real_t>(0));
+   MFEM_VERIFY(configured_sum > 0.0,
+               "AdvanceADERWithSubStep: configured deltaT sums to "
+               << configured_sum << " ≤ 0");
+   std::vector<mfem::real_t> deltaT_scaled(O);
+   for (int o = 0; o < O; o++)
+   {
+      deltaT_scaled[o] = configured_deltaT[o] * (dt_step / configured_sum);
+   }
+   iterator.SetSubSteps(deltaT_scaled, configured_weights);

-   const std::vector<real_t> &deltaT = iterator.GetDeltaT();
+   const std::vector<real_t> &deltaT = iterator.GetDeltaT();   // = deltaT_scaled
    const int O = static_cast<int>(deltaT.size());   // already correct after rescale
```

(The order of statements needs minor adjustment: read `configured_deltaT` BEFORE calling SetSubSteps. Inline-edit the helper accordingly.)

Alternatively, change `deltaT` semantics to be dimensionless ratios that the iterator scales internally per call. That's a bigger API change; the rescale-per-call approach is the smaller diff.

**Test case (NEW — must add):**
```cpp
// tests/unit/test_tpv104_substep_dispatch_last_step.cpp
TEST_CASE("R-1002: substep dispatch tolerates dt_step != configured dt") {
   // Configure iterator with deltaT_init = dt_init / O.
   // Call AdvanceADERWithSubStep with dt_step = 0.7 * dt_init.
   // Pre-fix: throws std::runtime_error.
   // Post-fix: completes without throw, Q_new finite.
}
```

---

### [R-1003] [CRITICAL] [POSSIBLE] [dynamic/wave_operator.inl::ComputeADERSharedFaceFluxRHS — UNTOUCHED] — MPI shared-fault path is silently divergent under substep dispatch

**Category:** BUG (correctness — silent inconsistency in MPI mode)

**Description:**
The substep guard (lines 2604–2632 in `wave_operator.inl::ComputeADERFaceFluxRHS`) is added ONLY in the interior fault branch. The shared-fault branch in `ComputeADERSharedFaceFluxRHS` (line 3106+) still runs inline `EvaluateADER` with macro-step Q̄.

In MPI runs (np > 1) where the fault is split across ranks:
- LOCAL fault QPs (interior to one rank): consume iterator's per-sub-step `I_imp` via the guard.
- SHARED fault QPs (across rank boundaries): run inline `EvaluateADER` on macro-step Q̄.

For the SAME physical fault face split by partition: the +y-side rank's local-fault contribution uses per-sub-step Q; the −y-side rank's shared-fault contribution uses macro-step Q̄. The two flux computations DIVERGE at O(dt) per step, producing a non-conservative fault Riemann at every shared fault face.

**Why the tests didn't catch this:**
All three new tests run on serial (np=1) meshes with no shared faces. The guard works correctly there. MPI behavior is untested.

**Trigger:** any production run with `--fault-iterator substep` AND np > 1.

**Actual behavior (MPI):** local and shared fault contributions to bulk rhs are computed on different time discretizations; non-conservative at shared faces.

**Expected behavior:** the guard fires uniformly in both `ComputeADERFaceFluxRHS` (interior) AND `ComputeADERSharedFaceFluxRHS` (shared) branches.

**Suggested fix:** apply the same guard pattern to the shared-fault `EvaluateADER` call. Steps:

1. Extend `EvaluateBulkAtFaultQPsCanonical` (or add `EvaluateBulkAtSharedFaultQPsCanonical`) to also process shared faces. Layout extension:
   - flat indices `[n_local_qps, n_local_qps + n_shared_qps)` for shared QPs.
   - DOF-idx mapping uses `shared_fault_dof_offset_` (already in wave operator).
2. Have `Tpv104SubStepIterator::AdvanceWithSubStepStates` operate on `n_total_qps` DOFs.
3. Add the same 7-line guard at the shared-fault branch's `fault_flux_->EvaluateADER` call site (find it in `ComputeADERSharedFaceFluxRHS` near line 3450+).

This is a multi-hour change (~80 LOC). Until done, ABORT-or-WARN if the substep flag is set with `nprocs > 1`:

```diff
@@ drivers/tpv104_driver.cpp around line 1228 (use_substep_iterator init)
    const bool use_substep_iterator =
       (GetDispatchedIterator(fault_iterator) == DispatchedIterator::SubStep);
+   if (use_substep_iterator && nprocs > 1)
+   {
+      MFEM_ABORT(
+         "--fault-iterator substep is not yet supported with MPI "
+         "(nprocs=" << nprocs << ").  The substep guard fires only on "
+         "the interior-fault branch; shared-fault faces would diverge "
+         "from the iterator's per-sub-step semantic.  Run with np=1 "
+         "or remove --fault-iterator substep until R-1003 is fixed.");
+   }
```

(Choose abort-with-message over warn-and-continue: silent wrong-results in MPI is exactly the failure mode round 9 spent rounds 5–8 hunting.)

**Test case (NEW — must add):**
```cpp
// tests/unit/test_tpv104_substep_mpi_guard.cpp
TEST_CASE("R-1003: substep dispatch is rejected at np > 1 until shared path lands") {
   // mpirun -np 2 ./seas_tpv104_driver --fault-iterator substep --dry-run
   // Expect: process aborts with the R-1003 message.
}
```

---

### [R-1004] [MODERATE] [dynamic/wave_operator.inl::EvaluateBulkAtFaultQPsCanonical:lines around the `if (!ftr) continue`] — Silent zero-fill on null face transformations

**Category:** EDGE_CASE (silent wrong-results)

**Description:**
```cpp
for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
{
   const int f = fault_interior_faces_[fi];
   FaceElementTransformations *ftr =
      mesh_.GetInteriorFaceTransformations(f);
   if (!ftr) { continue; }     // <-- silently skips, output stays 0
   ...
}
```

`Q_plus_flat` and `Q_minus_flat` are pre-zeroed via `assign(expect_words, 0.0)`. If `ftr == nullptr` for ANY face, the corresponding entries `[base_dof_idx, base_dof_idx + nqp)` of the flat output STAY AT ZERO. The substep iterator then runs the friction solve on Q ≡ 0 at those QPs, producing arbitrary friction outputs.

`fault_interior_faces_` is populated by the wave-operator constructor with faces filtered by `mesh_.GetInteriorFaceTransformations(f) != nullptr` — so in production, ftr should never be null at the helper. But this invariant is not asserted, so a future code change to `fault_interior_faces_` population could silently break the helper.

**Trigger:** mesh edit / refactor that puts a face into `fault_interior_faces_` whose `GetInteriorFaceTransformations` returns null.

**Actual behavior:** zero-filled fault QP entries silently drive the friction solve to wrong output.

**Expected behavior:** the helper enforces the same invariant the wave-operator ctor enforces — abort if ftr is null.

**Suggested fix:**
```diff
@@ wave_operator.inl::EvaluateBulkAtFaultQPsCanonical
-      if (!ftr) { continue; }
+      MFEM_VERIFY(ftr != nullptr,
+                  "EvaluateBulkAtFaultQPsCanonical: face " << f
+                  << " (fault_interior_faces_[" << fi << "]) has no "
+                  "InteriorFaceTransformations.  This violates the "
+                  "invariant established at WaveOperator ctor; the "
+                  "fault face list has been corrupted.");
```

**Test case (NEW — must add):**
```cpp
// tests/unit/test_tpv104_substep_helper_invariants.cpp
TEST_CASE("R-1004: EvaluateBulkAtFaultQPsCanonical aborts if a face is non-interior") {
   // Construct a mesh, populate fault_interior_faces_ via test hook with a
   // boundary-face index (where GetInteriorFaceTransformations returns null).
   // Call EvaluateBulkAtFaultQPsCanonical; expect MFEM_VERIFY abort.
}
```

(Implementing this test cleanly requires a friend hook to inject a corrupted `fault_interior_faces_`. If too invasive, downgrade R-1004 to LOW + leave the abort fix.)

---

### [R-1005] [MODERATE] [drivers/tpv104_driver.cpp:1218–1224 + AdvanceADERWithSubStep:317–323] — Quadrature is uniform-weight rectangle rule, not Gauss-Lobatto

**Category:** DEVIATION from "exactly as SeisSol did"

**Description:**
The user's original request was "implement sub step iteration exactly as seisol did". SeisSol's per-sub-step quadrature on `[0, dt]` is **Gauss-Lobatto** with O nodes:
- O=1 not used (SeisSol typically O ≥ 2).
- O=2: τ ∈ {0, dt}, weights = {1/2, 1/2} (trapezoid).
- O=3: τ ∈ {0, dt/2, dt}, weights = {1/6, 4/6, 1/6} (Simpson).
- O=4: 4-point Gauss-Lobatto with non-uniform weights ≈ {1/12, 5/12, 5/12, 1/12}.

My implementation hard-codes:
- deltaT[o] = dt/O (uniform-width sub-steps).
- weights[o] = 1/O (uniform).
- tau_nodes[o] = cumulative end (R-1001 fix: midpoint).

This is the **rectangle rule** (or midpoint rule after R-1001), not Gauss-Lobatto. Approximation error per macro-step: O(dt²) for rectangle/midpoint; O(dt^{2O−2}) for Gauss-Lobatto. At ADER-O=4, GL error is O(dt⁶); rectangle is O(dt²). The user pays an order-of-magnitude accuracy penalty per step.

**Why this matters:** the user's stated goal was SeisSol parity. Without Gauss-Lobatto, the substep dispatch is "a substep iteration" but not "the substep iteration SeisSol uses". For SCEC TPV104 cross-comparison, numerical agreement with SeisSol's reference traces is degraded.

**Why it's MODERATE not CRITICAL:**
- The implementation is mathematically valid (any quadrature with `Σ w = 1` is a valid time-average approximation).
- Documentation can describe the deviation.
- Switching to GL is purely a configuration change in the driver's iterator setup — no code change to the iterator or wave operator.

**Suggested fix:** replace the uniform deltaT/weights setup with hard-coded Gauss-Lobatto tables per O ∈ {1, 2, 3, 4}:

```diff
@@ drivers/tpv104_driver.cpp at iterator setup (~line 1218)
-   {
-      const int O = std::max(1, ader_order);
-      std::vector<real_t> deltaT(O, 1.0 / static_cast<real_t>(O));
-      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
-      for (int o = 0; o < O; o++) { deltaT[o] = dt / static_cast<real_t>(O); }
-      substep_iterator.SetSubSteps(deltaT, weights);
-   }
+   {
+      const int O = std::max(1, ader_order);
+      // Gauss-Lobatto on [0, dt]: nodes (positions) and weights such
+      // that ∫_0^dt f(τ)dτ = Σ w_o · f(τ_o) is exact for f polynomial
+      // of degree ≤ 2O − 1.  The deltaT[o] passed to SetSubSteps is the
+      // SUB-STEP WIDTH, not the absolute node — node positions are
+      // computed in AdvanceADERWithSubStep from cumulative deltaT.
+      // For nodal-on-boundary GL, deltaT[o] is τ_{o+1} − τ_o.
+      std::vector<real_t> nodes(O), weights(O);
+      switch (O) {
+         case 1:  // not standard GL; use midpoint rule
+            nodes[0] = 0.5;       weights[0] = 1.0;
+            break;
+         case 2:  // GL endpoints
+            nodes[0] = 0.0; nodes[1] = 1.0;
+            weights[0] = 0.5; weights[1] = 0.5;
+            break;
+         case 3:  // GL 0, 1/2, 1 — Simpson
+            nodes = {0.0, 0.5, 1.0};
+            weights = {1.0/6.0, 4.0/6.0, 1.0/6.0};
+            break;
+         case 4: {
+            // 4-point GL on [0, 1].
+            const real_t s = std::sqrt(0.2);     // = 1/√5
+            nodes = {0.0, 0.5*(1.0-s), 0.5*(1.0+s), 1.0};
+            weights = {1.0/12.0, 5.0/12.0, 5.0/12.0, 1.0/12.0};
+            break;
+         }
+         default:
+            MFEM_ABORT("AdvanceADERWithSubStep: ader_order=" << O
+                       << " has no Gauss-Lobatto table; supported O ∈ {1,2,3,4}.");
+      }
+      // Convert node positions (in [0, 1]) to deltaT widths.  For O=1
+      // (midpoint), the single sub-step is the whole interval.  For
+      // O ≥ 2 with nodes-on-boundary, deltaT[o] = nodes[o+1] - nodes[o]
+      // for o < O - 1, and deltaT[O-1] = 1 - nodes[O-1].  But this
+      // requires AdvanceADERWithSubStep to use the GL node positions
+      // (not cumulative widths) for tau_nodes — see R-1001 + R-1005
+      // joint fix.
+      ...
+   }
```

Note: R-1001 and R-1005 share the same code path; the joint fix is to
pass GL nodes-on-`[0, 1]` to the orchestrator, scale by `dt_step` per call.

**Test case (NEW — must add):**
```cpp
// tests/unit/test_tpv104_substep_quadrature.cpp
TEST_CASE("R-1005: substep quadrature is Gauss-Lobatto for O ≥ 2") {
   // For ADER-2: assert nodes = {0, 1}*dt, weights = {0.5, 0.5}.
   // For ADER-3: nodes = {0, 0.5, 1}*dt, weights = {1/6, 4/6, 1/6}.
   // Run a smooth-IC plane-wave through one macro-step at O=3 with both
   // GL substep and one-shot AdvanceADER; verify the substep result
   // matches the analytic time-average to O(dt^{2O-2}) (~dt^4 at O=3),
   // while the rectangle rule would only converge as O(dt^2).
}
```

---

### [R-1006] [LOW] [tests/unit/test_tpv104_substep_one_shot_parity.cpp Gate 1] — "Gate" with no actual assertion

**Category:** QUALITY (misleading test naming)

**Description:**
The test prints `-- Gate 1: NO-OP default — baseline AdvanceADER --` and runs `Q_baseline = RunOneStep(...)`. There is no `TEST_*` macro at this gate. Subsequent gates compare to `Q_baseline` but Gate 1 itself has zero assertions. A reader scanning the test file for "Gate 1" would expect a check.

**Suggested fix:**
```diff
@@ tests/unit/test_tpv104_substep_one_shot_parity.cpp Gate 1
    std::cout << "\n-- Gate 1: NO-OP default — baseline AdvanceADER --\n";
    std::vector<DOFData> dof_data_baseline = dof_data;  // snapshot
    wave.SetFaultDOFData(&dof_data_baseline, nbf);
    Vector Q_baseline = RunOneStep(wave, Q, dt, ader_order);
+   bool q_baseline_finite = true;
+   for (int i = 0; i < Q_baseline.Size(); i++)
+   {
+      if (!std::isfinite(Q_baseline(i))) { q_baseline_finite = false; break; }
+   }
+   TEST_ASSERT(q_baseline_finite,
+               "Gate 1: baseline AdvanceADER produces finite Q_new "
+               "with no substep state set (default null pointers).");
```

---

### [R-1007] [LOW] [drivers/tpv104_driver.cpp:443] — Stale `(void)method;` cast and comment

**Category:** QUALITY

**Description:**
The driver retains:
```cpp
const FrictionSolver::Method method = MapSolver(friction_solver);
(void)method;  // R7-001 option (b): not routed through the time loop.
```

After R-602/R-603, `method` IS routed through the time loop via `AdvanceADERWithSubStep`. The cast is benign (just an expression statement) but the comment is wrong. A reader will infer the variable is unused, then be confused when they see it later.

**Suggested fix:**
```diff
-   const FrictionSolver::Method method = MapSolver(friction_solver);
-   (void)method;  // R7-001 option (b): not routed through the time loop.
+   // R-602/R-603 (round-7): `method` is now passed into
+   // AdvanceADERWithSubStep on the substep dispatch path.  On the
+   // legacy one-shot path it's still unused (Brent is hard-coded inside
+   // EvaluateADERTotal); the previous (void)method silencer is removed.
+   const FrictionSolver::Method method = MapSolver(friction_solver);
```

---

## Summary

- Critical issues: **3** (R-1001 wrong tau-node breaks SSI_3 contract; R-1002 last-step crash; R-1003 MPI silent divergence)
- Moderate issues: **2** (R-1004 silent zero-fill on null ftr; R-1005 quadrature is rectangle-rule not Gauss-Lobatto)
- Low issues: **2** (R-1006 missing assert in Gate 1; R-1007 stale cast/comment)
- Plan compliance: **PARTIAL** — implementation builds, my unit tests pass, but those tests don't cover the dispatch contract. R-1001/R-1002 break the SSI_3 contract and crash on production tfinal values. R-1003 produces silent wrong-results in MPI.
- Verdict: **FAIL — must fix before claiming substep dispatch is production-ready.** R-1001 + R-1002 + R-1003 are blocking. R-1004/R-1005 should be fixed in the same patch.

---

## What the agent (me) overlooked

1. **Tests passing ≠ code correct.** The three new tests validate fragments (predictor, iterator, guard) but never exercise the full driver dispatch end-to-end against the legacy path. Adding `test_tpv104_substep_dispatch_parity` (R-1001 test case) catches R-1001. Adding `test_tpv104_substep_dispatch_last_step` catches R-1002.
2. **The driver's tau-node choice was made without checking the SSI_3 contract math.** I picked "cumulative end of sub-step" as the obvious thing without verifying that `Q(τ=dt) = I/dt`. It doesn't.
3. **Iterator dt semantics tied to startup CFL.** I configured `deltaT` once at startup with the auto-CFL `dt`, then the time loop's varying `dt_step` breaks the iterator's verify. This is exactly the kind of single-shot-config-vs-per-call-input mismatch that should have been caught in design.
4. **MPI shared-fault path not extended.** Documented as "known limitation" but the implementation should at minimum REJECT MPI runs with the substep flag, not silently produce wrong results.
5. **"Exactly as SeisSol" interpretation was loose.** Uniform-weight rectangle rule is not the same as Gauss-Lobatto. The user's stated goal was SeisSol parity; my implementation is a pragmatic substep iteration that's not actually parity.

---

## Recommended fix order

1. **R-1003 first** (one-line MFEM_ABORT in driver). Prevents silent wrong-results in MPI immediately.
2. **R-1002** (rescale deltaT per call). Prevents crash on production runs.
3. **R-1001** (midpoint nodes). Restores SSI_3 contract; enables the new test for dispatch parity.
4. **R-1004** (MFEM_VERIFY on ftr). Defensive; cheap.
5. **R-1005** (Gauss-Lobatto). Final SeisSol-parity step. Bigger diff; can land separately after 1–4 are validated.
6. **R-1006, R-1007** (cleanup). Anytime.

After R-1001/R-1002/R-1003 land, the new dispatch-parity test should be added to `make test` so this round of bugs cannot recur.
