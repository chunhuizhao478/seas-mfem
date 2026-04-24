# Code Review Round 5: TPV104 Step 7 post-R4 fixes + tolerance audit — 2026-04-24

Fresh adversarial audit after the round-4 closures landed. Primary
focus of this round: (a) confirm R4-001..R4-008 are actually closed as
claimed, (b) audit unit-test tolerances for looseness that could hide
bugs, (c) hunt for new bugs introduced by the R4 fixes themselves.

## ⚠ STANDING PLAN VIOLATION — open across four review rounds

**Plan §4.10 Step 5 mandate**: *"Newton solver is the driver's
default; `FrictionSolver::SolveNR` (legacy MFEM μ formula) is only
reached under `--friction-solver=legacy-newton`. Test T_TPV104_FS_4
stub is replaced with a real smoke-test assertion."* The intended
Newton is `SolveSlipRateNewtonStable` + `FrictionCoefficientStable`
(Step 5 of the plan) — the SeisSol-aligned stable-asinh Newton that
matches `rs::arsinhexp` byte-for-byte on the TPV104 envelope.

**Current state (unchanged since round 4):**

- `FrictionSolver::Method` has THREE values: `Brent`, `NewtonRaphson`,
  `HybridNRBisection`. There is NO `NewtonRaphsonStable` enum value.
- `Method::NewtonRaphson` dispatches to `FrictionSolver::SolveNR`
  (`dynamic/friction_solver.cpp:61`), which uses MFEM's **legacy** μ
  via `DieterichRuinaFriction::FrictionCoefficientPsi` — the path
  with the `if (psi_over_a > 700.0)` hard switch.
- `SolveSlipRateNewtonStable` (`dynamic/tpv104_friction_solver.hpp`)
  exists but is **not wired through `FrictionSolver::Solve`** — the
  only way to reach it is a direct call that bypasses
  `FaultFaceFlux::ComputeStageState`. The iterator
  (`Tpv104SubStepIterator::Advance`) dispatches through
  `FaultFaceFlux::ComputeStageState(…, method)` → `FrictionSolver::Solve(method)`
  → one of `{SolveBrent, SolveNR, SolveHybrid}`, none of which use
  `FrictionCoefficientStable`.
- The iterator's default is now `Method::Brent` (round-4 R4-002
  "closure"). Plan compliance status: **Step 5 mandate NOT met**.

**Using `Method::NewtonRaphson` as-if SeisSol on TPV104:** numerically
the two μ formulas agree to ~1e-12 relative on the TPV104 physical
envelope (`ψ/a ∈ [28, 80]`) because the 700-branch never fires. So
shipping with `Method::NewtonRaphson` produces physically-correct
TPV104 answers. But:

1. The code path is NOT byte-aligned with SeisSol's
   `invertSlipRateIterative` — Phase-3 Probe 3 diffs would show two
   independent formulas, not byte-match.
2. The plan Step 5 rationale ("Newton-stable as TPV104 canonical") is
   defeated — we ship a solver that is only incidentally equivalent.
3. `TestBannerDefaults` (`test_tpv104_smoke.cpp:86-90`) requires the
   driver banner to print `"Friction solver: Newton-Raphson"` on
   default flags AND forbids `"Friction solver: Brent"` — inconsistent
   with the iterator's `Method::Brent` default.

**Severity escalation this round:** R5-003 was labelled MODERATE in
the initial writeup; this was a mis-ranking. Per skill rules, plan
deviations outrank local-numerical-impact assessments. **R5-003 is
CRITICAL (standing plan violation across four review rounds).**

**Unblocking the deviation — three concrete choices, only option (a)
matches the plan:**

1. **(a) Add `Method::NewtonRaphsonStable`** — ~15 lines: one enum
   value in `friction_solver.hpp`, a dispatch case in
   `friction_solver.cpp`, iterator default bumped to the new value,
   banner string updated to `"Friction solver: Newton-Raphson
   (stable-asinh)"`. This is the plan-compliant path.
2. **(b) Bypass `FrictionSolver::Solve` in the Step-9 driver** —
   larger scope: driver reaches into `FaultFaceFlux` internals or
   replaces `ComputeStageState` with a TPV104-specific path that
   calls `SolveSlipRateNewtonStable` directly.
3. **(c) Amend the plan** to accept Brent as the TPV104 default and
   update `test_tpv104_smoke.cpp` banner expectations. **This is the
   path shipped today by default.** It requires user-authorised
   plan update; no author has approved this yet.

Tracked across rounds:

| round | ID | severity | status |
|---|---|---|---|
| 1 | Deviation #4 | MODERATE | Deferred to Step 9 |
| 4 | R4-002 | MODERATE | Partially closed (default moved NewtonRaphson → Brent; plan mandate still unmet) |
| 5 | R5-003 | **CRITICAL (this round)** | **OPEN** |

Full R5-003 details in the Findings section below.

## Review Scope

- Plan: `tpv104_debug_plan_2026-04-24.md` §4.10 Step 7.
- Prior reviews (all earlier findings closed):
  - `tpv104_review_2026-04-24.md` (R-001..R-012)
  - `tpv104_review_round2_2026-04-24.md` (R2-001..R2-007)
  - `tpv104_review_round3_2026-04-24.md` (R3-001..R3-007)
  - `tpv104_review_round4_2026-04-24.md` (R4-001..R4-008)
- Files reviewed (changed since round 4):
  - `dynamic/tpv104_substep_iterator.hpp` (157 → 165)
  - `dynamic/tpv104_substep_iterator.cpp` (374 → 457)
  - `tests/unit/test_tpv104_substep_iterator.cpp` (553 → 701)
  - `tests/unit/test_tpv104_smoke.cpp` (190 → 193)
  - `dynamic/tpv104_setup.hpp` (522 → 541)
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, round-4 review.

## Round-4 closure verification

- **R4-001 CLOSED**: `tpv104_substep_iterator.cpp:371-372` adds
  `d.slip1 += s.V1 * dt_sub; d.slip2 += s.V2 * dt_sub;` inside the
  per-sub-step loop.  Header docstring (line 90-93) lists
  `slip1/slip2` in the mutated fields. New test `TestSlipAccumulation`
  at `test_tpv104_substep_iterator.cpp:558-634` — **but with issues,
  see R5-002 below**.
- **R4-002 PARTIALLY CLOSED**: iterator default is now
  `FrictionSolver::Method::Brent` (header line 147, definition
  comment line 160). However, the plan §4.10 Step 5 mandate —
  "Newton solver is the driver's default" — is still unmet; Brent is
  a fallback, NOT Step-5's `SolveSlipRateNewtonStable` +
  `FrictionCoefficientStable` path. See R5-003.
- **R4-003 CLOSED**: Probe 3 emission now uses
  `friction_stable::FrictionCoefficientStable(s.V_abs, d.psi, d.a,
  V0_scalar)` at line 327-329, matching the solver's intended μ (when
  Step-9 wiring eventually dispatches through it).
- **R4-004 PARTIALLY CLOSED**: Probe-1 header adds a
  `cadence_note` disclosing the Q̄-averaging deviation
  (line 90-97). The underlying Q̄ treatment is unchanged — this is
  doc-only, as flagged in the fix. Acceptable for Phase 2.
- **R4-005 CLOSED**: `GetProbeFile` now falls back to
  `MPI_Comm_rank(MPI_COMM_WORLD, …)` when the env var is unset and
  MPI is initialised (line 60-74).
- **R4-006 CLOSED (subsumed)**: Probe 3's V-guard now uses
  `(s.V_abs > 0.0) ? FrictionCoefficientStable(...) : 0.0` — the
  `V ≤ 0` branch returns 0, which matches the algebraic identity
  `asinh(0) = 0`. R-005 already documents μ as odd in V; the guard is
  conservative and harmless.
- **R4-007 OPEN**: Station-writer MPI tie-break tolerance still
  1e-10 absolute; no abort guard on `winning_rank == nprocs_loc`.
  Still a latent risk.
- **R4-008 CLOSED**: Test renamed `TestSubStepAccumulatorSelfConsistency`
  at `test_tpv104_substep_iterator.cpp:103` with the acknowledgment
  comment. A new test `TestAccumulatorScaleIdentity` (line 645-682)
  provides the straight-line accumulator check without production-helper
  calls.

## Findings

### [R5-001] CRITICAL [test_tpv104_substep_iterator.cpp:482-491] — `TestConvergenceUnderDtHalving` computes `diff = |X − X|`, a self-difference that is identically zero — the T_TPV104_SSI_2 gate has zero discriminative power

**Category:** BUG (test logic error)

**Description:**
```cpp
real_t max_rel = 0.0;
for (int c = 0; c < NUM_STATE; ++c)
{
   const real_t ref = std::max<real_t>(std::abs(I_imp_plus[c]),
                                       static_cast<real_t>(1e-30));
   const real_t diff = std::abs(I_imp_plus[c] - I_imp_plus[c]);  // <-- BUG
   // Self-difference is 0 — this is a stability probe, not a
   // convergence probe.
   max_rel = std::max(max_rel, diff / ref);
}
return max_rel;
```
`I_imp_plus[c] - I_imp_plus[c]` is **identically zero** for every `c`.
`max_rel` is always 0. Both `TEST_ASSERT(r_macro < 1e-12, ...)` and
`TEST_ASSERT(r_macro_half < 1e-12, ...)` trivially pass regardless of
anything the iterator produces.

The comment "Self-difference is 0 — this is a stability probe, not a
convergence probe" acknowledges the problem but ships it anyway. The
plan's T_TPV104_SSI_2 gate (§4.10 Step 7) calls for an O(dt²)
convergence test under dt halving — exact 4× reduction of residual
magnitude. The current implementation:

1. Doesn't compare the two dt runs to each other (the second return
   value is unused past the `std::cout` print).
2. Doesn't compare against any expected value (uses `I_imp_plus[c]`
   on BOTH sides of the subtraction).
3. Produces a trivially true assertion.

**The iterator could emit NaN on both dt runs and this test would still
pass.** A rupture-onset regression in `Advance` that silently broke
the imposed-state accumulation would not trip T_TPV104_SSI_2.

**Trigger:** Any call path that exercises `TestConvergenceUnderDtHalving`
— the test's assertion is independent of the code under test.

**Actual behavior:** Assertions always pass; T_TPV104_SSI_2 provides
no coverage.

**Expected behavior:** Compare output between the two dt runs, OR
compare each run's output against an analytically-known fixture.

**Suggested fix:**
```diff
-   auto run_and_measure = [](real_t dt_macro) -> real_t
+   auto run_and_get_imposed_state = [](real_t dt_macro,
+                                       std::vector<real_t> &I_imp_plus_out)
+      -> void
    {
       std::vector<DOFData> dof_data;
       std::vector<Vector> fault_coords;
       std::vector<real_t> V_w;
       MakeSingleQPFixture(dof_data, fault_coords, V_w,
                           /*x2=*/20e3, /*x3=*/20e3);
       ...
-      std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
+      I_imp_plus_out.assign(NUM_STATE, 0.0);
       std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);
 
       it.Advance(dof_data, fault_coords, V_w,
                  I_plus.data(), I_minus.data(),
                  dt_macro, 2.0 * TPV104Params::nuc_T,
-                 I_imp_plus.data(), I_imp_minus.data(),
+                 I_imp_plus_out.data(), I_imp_minus.data(),
                  FrictionSolver::Method::Brent);
-
-      // For Q_bar = 0 + rest state, ...
-      real_t max_rel = 0.0;
-      for (int c = 0; c < NUM_STATE; ++c)
-      {
-         const real_t ref = std::max<real_t>(std::abs(I_imp_plus[c]),
-                                             static_cast<real_t>(1e-30));
-         const real_t diff = std::abs(I_imp_plus[c] - I_imp_plus[c]);
-         ...
-         max_rel = std::max(max_rel, diff / ref);
-      }
-      return max_rel;
    };
 
-   const real_t r_macro     = run_and_measure(4.7e-3);
-   const real_t r_macro_half = run_and_measure(4.7e-3 / 2.0);
-   std::cout << "  residual at dt=4.7 ms      = " << r_macro << "\n";
-   std::cout << "  residual at dt=4.7 ms / 2  = " << r_macro_half << "\n";
-   TEST_ASSERT(r_macro < 1e-12,
-               "iterator stable under dt_macro = 4.7 ms");
-   TEST_ASSERT(r_macro_half < 1e-12,
-               "iterator stable under dt_macro/2 (no amplified round-off)");
+   // For the rest-state fixture (Q̄ = 0), the time-integrated I_imp
+   // scales linearly with dt_macro.  Halving dt_macro halves every
+   // I_imp_plus[c] to O(dt²); the RATIO of components at the two
+   // dt values should be ≈ 0.5 + O(dt).  We verify ratio ∈ [0.49, 0.51].
+   std::vector<real_t> I_full, I_half;
+   run_and_get_imposed_state(4.7e-3,       I_full);
+   run_and_get_imposed_state(4.7e-3 / 2.0, I_half);
+   real_t max_ratio_dev = 0.0;
+   for (int c = 0; c < NUM_STATE; ++c) {
+      const real_t mag_full = std::abs(I_full[c]);
+      if (mag_full < 1e-20) { continue; }  // skip zero channels
+      const real_t ratio = std::abs(I_half[c]) / mag_full;
+      max_ratio_dev = std::max(max_ratio_dev, std::abs(ratio - 0.5));
+   }
+   std::cout << "  max |ratio - 0.5| = " << max_ratio_dev << "\n";
+   TEST_ASSERT(max_ratio_dev < 0.01,
+               "I_imp(dt/2) / I_imp(dt) ≈ 0.5 ± 1 % (O(dt²) convergence)");
```

**Test case:**
```cpp
void test_R5_001_convergence_test_actually_tests_convergence() {
   // Deliberately break Advance() by multiplying I_imp by 0 on the
   // second half of the sub-steps.  The CURRENT test passes (trivially).
   // The FIXED test must catch this.
   //
   // Compile with a feature flag:
   //   #define BREAK_CONVERGENCE_FOR_TEST
   // and inside Advance:
   //   if (o > O/2) { /* scale Q_imp_plus by 2 */ }
   //
   // Run the test; with the fix, max_ratio_dev > 0.01 → FAIL.
   // Without the fix (self-difference), PASS — bug goes undetected.
}
```

---

### [R5-002] MODERATE [test_tpv104_substep_iterator.cpp:615-633] — `TestSlipAccumulation` does not catch sign errors; `slip2_after != slip2_before` is passive and the `|slip2_after − slip2_before|` magnitude window ignores direction

**Category:** BUG (test tolerance too loose)

**Description:**
The slip-accumulation test at `test_tpv104_substep_iterator.cpp:615-633`:
```cpp
TEST_ASSERT(slip2_after != slip2_before,
            "slip2 accumulates when V2 ≠ 0 (R4-001 guard)");
...
TEST_ASSERT(std::abs(slip2_after - slip2_before)
               >= 0.2 * expected_slip2_order,
            "slip2 accumulation ≥ 0.2·V2·dt_macro (order-of-mag)");
TEST_ASSERT(std::abs(slip2_after - slip2_before)
               <= 5.0 * expected_slip2_order,
            "slip2 accumulation ≤ 5·V2·dt_macro (no spurious amplification)");
```
Three assertions on `|slip2_after - slip2_before|`:
1. `!= slip2_before` — passes for any non-zero diff.
2. `≥ 0.2 · |V2|·dt` — magnitude lower bound.
3. `≤ 5.0 · |V2|·dt` — magnitude upper bound (25× window).

A sign bug in the iterator — e.g., `d.slip2 -= s.V2 * dt_sub;` instead
of `+=` — would produce a NEGATIVE diff. Taking `std::abs(...)` of it
gives the correct magnitude, so all three assertions still pass.

This is non-hypothetical: sign errors in slip accumulation are a
well-known failure mode (see `miniapps/seas/CLAUDE.md` §"Sign
Conventions" — "Slip rate direction: PARALLEL to traction tau (not
antiparallel). … Antiparallel sign creates positive feedback ->
unbounded growth (debug v8).").

**Trigger:** Sign inversion in `d.slip1 += s.V1 * dt_sub` or
`d.slip2 += s.V2 * dt_sub`.

**Actual behavior:** Sign inversion passes all three assertions.

**Expected behavior:** Assert that `(slip2_after - slip2_before)` has
the same sign as the expected slip direction. For the test fixture
(`I_plus[VZ] = +0.5*dt, I_minus[VZ] = -0.5*dt`, producing
right-lateral slip in +z direction), `slip2_after - slip2_before`
should be positive.

**Suggested fix:**
```diff
    TEST_ASSERT(slip2_after != slip2_before,
                "slip2 accumulates when V2 ≠ 0 (R4-001 guard)");
+
+   // R5-002: directional check — slip must accumulate with the
+   // same sign as the sub-step V2.  A sign error in
+   // `d.slip2 += s.V2 * dt_sub` (becoming `-=`) would still
+   // satisfy the magnitude assertions below but flip this check.
+   TEST_ASSERT(std::signbit(slip2_after - slip2_before)
+                  == std::signbit(V2_final),
+               "sign(Δslip2) == sign(V2_final) — no sign error in "
+               "slip accumulation (R5-002 guard)");
 
    const real_t expected_slip2_order = std::abs(V2_final) * dt_macro;
    TEST_ASSERT(std::abs(slip2_after - slip2_before)
                   >= 0.2 * expected_slip2_order,
```
And tighten the window from `[0.2, 5.0]` × to `[0.8, 1.2]` × by
running the test at a fixture where V varies < 20 % across sub-steps
(already the case — rest-state-adjacent), since the 25× window admits
5× accumulation-formula bugs.

**Test case:** Deliberately invert the sign in production
(`d.slip2 -= s.V2 * dt_sub`) and run the test. Under the current
assertions all three pass; under the R5-002 fix, the signbit check
fires.

---

### [R5-003] CRITICAL [tpv104_substep_iterator.hpp:146-147, test_tpv104_smoke.cpp:86-87, friction_solver.hpp:41, dynamic/tpv104_friction_solver.hpp] — Plan §4.10 Step 5 "Newton solver is the driver's default" is NOT met — SeisSol-aligned `SolveSlipRateNewtonStable` is not reachable via `FrictionSolver::Solve`, and iterator default is Brent (not Newton); standing plan violation across four review rounds

**Severity upgraded from MODERATE to CRITICAL 2026-04-24**: the
finding is a direct plan §4.10 Step 5 deviation that has been open
across rounds 1 (Deviation #4), 4 (R4-002 partial closure), and now
5 (R5-003). Per skill rules, plan deviations outrank local-numerical-
impact assessments; R5-003 belongs at CRITICAL until closed.

**Category:** DEVIATION (plan / test inconsistency)

**Description:**
Round-4 R4-002 closure changed the iterator default from
`Method::NewtonRaphson` to `Method::Brent`:
```cpp
FrictionSolver::Method method
   = FrictionSolver::Method::Brent);
```
Meanwhile `test_tpv104_smoke.cpp:86-90` still asserts:
```cpp
{"Friction solver: Newton-Raphson",
 "Newton-Raphson friction-solver default"},
...
const bool found = out.find(kv.first) != std::string::npos;
TEST_ASSERT(found, ...);
```
and forbids `"Friction solver: Brent"` on default flags at line 107.

When the Step-9 driver ships:
- If the driver passes `Method::Brent` to the iterator (mirroring the
  iterator default) and prints the banner accordingly → `"Friction
  solver: Brent"` — **T_TPV104_SMOKE_3 FAILS** (forbidden string).
- If the driver prints `"Friction solver: Newton-Raphson"` but passes
  `Method::Brent` to the iterator → banner lies to the user.
- If the driver passes `Method::NewtonRaphson` (overriding the
  iterator's R4-002 default) → T_TPV104_SMOKE_3 passes but R4-002's
  "Brent as default" is bypassed immediately.

There is no coherent resolution with the current state. The plan
§4.10 Step 5 mandate ("Newton solver is the driver's default")
remains unmet.

**Trigger:** Step-9 driver ships. All three paths are bad.

**Actual behavior:** Iterator-default and banner-test expectations
are inconsistent.

**Expected behavior:** Either (a) wire the stable-asinh Newton as a
FOURTH enum value (`Method::NewtonRaphsonStable`) dispatching to
`SolveSlipRateNewtonStable` + `FrictionCoefficientStable`, update
iterator default to `Method::NewtonRaphsonStable`, update smoke-test
banner to `"Friction solver: Newton-Raphson (stable-asinh)"`; OR
(b) accept Brent as the default, update the banner string and
remove the forbidden-Brent guard in smoke test.

**Suggested fix (option a — restore plan compliance):**
```diff
 // dynamic/friction_solver.hpp
-   enum class Method { Brent, NewtonRaphson, HybridNRBisection };
+   enum class Method {
+      Brent,
+      NewtonRaphson,              // legacy MFEM-native μ
+      NewtonRaphsonStable,        // Step-5 FrictionCoefficientStable
+      HybridNRBisection
+   };
```
Dispatch extension in `friction_solver.cpp:60-63`:
```diff
       case Method::NewtonRaphson: return SolveNR(tau, psi, sigma_n, eta, a);
+      case Method::NewtonRaphsonStable:
+         return SolveSlipRateNewtonStableWrapper(tau, psi, sigma_n,
+                                                 eta, a);
       case Method::HybridNRBisection: return SolveHybrid(tau, psi, sigma_n, eta, a);
       default: return SolveBrent(tau, psi, sigma_n, eta, a);
```
Iterator default:
```diff
                 FrictionSolver::Method method
-                   = FrictionSolver::Method::Brent);
+                   = FrictionSolver::Method::NewtonRaphsonStable);
```
Banner string (once Step 9 lands):
```cpp
// seas_tpv104_driver.cpp banner
std::cout << "Friction solver: Newton-Raphson (stable-asinh)\n";
```

**Test case:**
```cpp
void test_R5_003_driver_banner_matches_iterator_default() {
   // Inspect the driver's actual solver call path — the method
   // argument passed to Advance must match the banner string.
   // With option (a), banner="Newton-Raphson (stable-asinh)" and
   // method=Method::NewtonRaphsonStable.  Concordance test.
}
```

---

### [R5-004] MODERATE [tpv104_substep_iterator.cpp:229-272] — Iterator deactivates the R4-001 slip correctness at V<0; per-sub-step s.V2 writes negative slip when trial traction overshoots zero crossing

**Category:** EDGE_CASE (sub-step sign oscillation)

**Description:**
The R4-001 fix at line 371-372:
```cpp
d.slip1 += s.V1 * dt_sub;
d.slip2 += s.V2 * dt_sub;
```
uses `s.V1, s.V2` — the per-sub-step slip-rate components. These come
from `CompleteFromVabs` via `ComputeStageState`. If the bulk trial
traction oscillates across sub-steps (possible when the ADER predictor
over-shoots past a quasi-static equilibrium — a known nucleation
regime), the per-sub-step `s.V2` can switch sign across sub-steps.
The iterator would then accumulate `slip2 += +v*dt + −v'*dt + …`,
potentially producing `slip2` that is **smaller in magnitude** than
either the initial or final V2·dt_macro, because the sum of signed
increments can cancel across sub-steps.

Compare to `dof_data[i].slip_rate = |V|` which is always positive:
the current iterator accumulates SIGNED V components (v_t1, v_t2).
That is the correct physics (slip is a signed vector, not a
magnitude). But it means the slip accumulation's correctness depends
on the predictor not producing spurious sign oscillation.

Under R4-004 (Q̄ = I/dt_macro is constant across sub-steps), the
trial traction is identical per sub-step, so V1/V2 change only because
ψ and tau2_nuc change. For a rest-state ramp, V grows monotonically —
no oscillation risk. For a nucleation-triggered rupture where the
Newton solver can converge to different V on each sub-step due to ψ
relaxation, the signs are stable (fault remains in right-lateral
slip).

But the risk is present in a general fixture, and no test fires on
this corner. `TestSlipAccumulation` tests rest-state-adjacent; nothing
exercises sign-flip-per-sub-step.

**Trigger:** Phase-3 fixture where the friction solver returns `s.V2`
with opposite signs across consecutive sub-steps of the same macro-step.

**Actual behavior:** Slip accumulation sums signed increments; they can
partially cancel, leaving slip smaller than naively expected.

**Expected behavior:** Either (a) document that slip accumulation is
signed and the caller is responsible for ensuring V doesn't oscillate
across sub-steps, OR (b) add a diagnostic counter that flags sign
inversions across sub-steps, OR (c) use `s.V_abs * sign_of_V_prev` so
magnitude-only V is accumulated — but this changes the physics.

**Suggested fix:** document the invariant (add comment), and add a
diagnostic when SEAS_DIAG_TPV104_STATE is set:
```diff
          d.slip1 += s.V1 * dt_sub;
          d.slip2 += s.V2 * dt_sub;
+
+#ifdef SEAS_DIAG_TPV104_STATE
+         // R5-004: sign-inversion diagnostic — flag sub-step boundaries
+         // where V2 changes sign (indicates predictor pathology or
+         // friction-solver non-monotonicity).  Logged to the same
+         // per-QP probe file pattern.
+         static thread_local std::vector<real_t> V2_prev_per_qp;
+         if (static_cast<int>(V2_prev_per_qp.size()) <= i) {
+            V2_prev_per_qp.resize(i + 1, s.V2);
+         } else if (std::signbit(V2_prev_per_qp[i]) != std::signbit(s.V2)
+                    && s.V2 != 0.0 && V2_prev_per_qp[i] != 0.0) {
+            std::ofstream &f = GetProbeFile("slip_sign_inversion");
+            if (f.is_open()) {
+               f << std::scientific << std::setprecision(16)
+                 << (t_sub_cursor + dt_sub) << " " << i << " "
+                 << V2_prev_per_qp[i] << " " << s.V2 << "\n";
+            }
+         }
+         V2_prev_per_qp[i] = s.V2;
+#endif
```

**Test case:** Requires a fixture designed to trigger sign oscillation
— ad-hoc for R3 probe runs. Low priority in Phase 2.

---

### [R5-005] LOW [test_tpv104_substep_iterator.cpp:417] — `TestPsiUpdatePlacement` inline reference uses `V_sub = TPV104Params::V_ini = 1e-16` by convention, but iterator uses `s.V_abs` from the friction solver — brittle equivalence

**Category:** ASSUMPTION (test brittleness)

**Description:**
```cpp
const real_t V_sub = TPV104Params::V_ini;
psi_ref = UpdateStateAnalyticSlipLawSRW(
   psi_ref, V_sub, TPV104Params::L, deltaT[o],
   V_w[0], TPV104Params::a_in,
   TPV104Params::b, TPV104Params::V0,
   TPV104Params::f0, TPV104Params::f_w);
```
The inline reference assumes the friction solver returns `V_abs =
V_ini = 1e-16` at every sub-step for the rest-state fixture. That is
true IF `ψ_ini = ComputeInitialPsiTPV104(a_in)` is bit-exact inverted
such that τ = σ_n · μ(V_ini, ψ_ini, a) + η_s · V_ini — which it is by
construction at the point the fixture is built.

But the ψ USED in the friction solve is the per-sub-step-UPDATED ψ,
not ψ_ini. After the first sub-step, ψ may have drifted by ~1e-19
(the analytic step in the rest-state limit). The friction solver
called with this slightly-drifted ψ returns a slightly-drifted V_abs
— NOT exactly V_ini. The inline reference uses the constant V_ini
across all sub-steps; the iterator's V_sub drifts.

At TPV104 rest state the drift is ~1e-20 m/s per sub-step — well
below the 1e-12 `TEST_NEAR` tolerance. The test passes. But the
correspondence is brittle: any future fixture with V_sub that drifts
above 1e-12 relative would fail this inline-reference test while the
iterator is still correct.

**Trigger:** Future fixture with non-negligible per-sub-step V drift.

**Actual behavior:** Test may falsely fail.

**Expected behavior:** Drive the inline reference with the actual
per-sub-step V_abs extracted from the iterator (via a callback) OR
relax the bit-match claim to "rest-state drift < ε".

**Suggested fix:** Replace the inline reference with a trajectory-aware
comparison, or document the brittleness explicitly:
```diff
    // Independently compute ψ via an inline per-sub-step accumulator and
    // compare to iterator-updated ψ bit-for-bit.
+   // R5-005: inline-reference assumes V_sub = V_ini constant across
+   // sub-steps.  Valid at TPV104 rest state (V drift < 1e-20/sub-step).
+   // Any future fixture with non-trivial V drift must switch to
+   // trajectory-aware reference instead of the constant V_ini below.
    real_t psi_ref = psi_entry;
    for (size_t o = 0; o < deltaT.size(); ++o)
    {
-      const real_t V_sub = TPV104Params::V_ini;
+      const real_t V_sub = TPV104Params::V_ini;   // rest-state only
```

---

### [R5-006] LOW [tpv104_substep_iterator.cpp:37-38] — `GetProbeFile`'s static `std::unordered_map<std::string, std::unique_ptr<std::ofstream>>` may be destroyed after `MPI_Finalize` in some toolchains; flush timing is subtle

**Category:** ASSUMPTION (program-exit ordering)

**Description:**
```cpp
std::ofstream &GetProbeFile(const char *probe_name)
{
   static std::unordered_map<std::string, std::unique_ptr<std::ofstream>>
      files;
   ...
}
```
Function-local static with first-use initialisation; destroyed at
program exit by static-destructor registration. In an MPI program, if
the main function calls `MPI_Finalize()` and then the static
destructor of `files` runs after that, the `std::ofstream` destructors
flush buffered writes — but on some filesystems this post-finalize
flush can hit a torn down MPI I/O layer on SLURM job-step
cleanup. Probes may lose the last few hundred bytes.

**Trigger:** SLURM wall-time kill or `MPI_Finalize` called
before the `files` destructor fires, combined with a parallel
filesystem that does lazy write-back.

**Actual behavior:** Probe output may be truncated at program exit.

**Expected behavior:** Explicit `Flush()` / `Close()` entry point for
the probe file map, called from the driver's shutdown path before
`MPI_Finalize`.

**Suggested fix:** Expose a `CloseAllProbeFiles()` function and
instruct the driver to call it:
```diff
+void CloseAllProbeFiles();
+
 std::ofstream &GetProbeFile(const char *probe_name)
 {
    ...
 }
 }  // anonymous namespace
+
+void CloseAllProbeFiles()
+{
+   // No-op — defined for API symmetry; the static files map is
+   // cleared at program exit.
+}
```
In the driver shutdown path (Step 9):
```cpp
// Before MPI_Finalize:
mfem::seas::CloseAllProbeFiles();
```

---

### [R5-007] LOW [tpv104_substep_iterator.cpp:217-222] — `Σ deltaT == dt_macro` tolerance of 1e-10 relative allows drift too tight for 500-sub-step problems

**Category:** EDGE_CASE (tolerance calibration)

**Description:**
```cpp
const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro, 1e-300);
if (rel > 1e-10)
{
   throw std::runtime_error(...);
}
```
For a GL-quadrature with O = 500 (not used today, but plausible for
high-accuracy ADER), the accumulated sum of `weight[o] * dt_macro` has
round-off ~ O · ULP · dt_macro ≈ 500 * 2e-16 * dt_macro = 1e-13 *
dt_macro. Rel drift = 1e-13. That is below 1e-10, so no failure.

For O = 5 or 10 (practical ADER orders) the rel drift is ~1e-14 —
far below the threshold. 1e-10 is fine.

However, at the lower bound `max(dt_macro, 1e-300)`, if `dt_macro =
1e-300`, `rel = |1e-300 − 1e-300| / 1e-300 = 0`. Safe. No bug.

**Recommendation:** Document that the 1e-10 tolerance is calibrated
for O ≤ 1000; bake a lazy assert `O <= 1000` to prevent accidental
over-scaling. LOW.

---

## Tolerance audit: unit-test thresholds

A separate systematic review of every TEST_ASSERT / TEST_NEAR
threshold across the TPV104 test suite, looking for tolerances too
loose to catch real bugs.

| Test | Location | Threshold | Concern |
|---|---|---|---|
| `TestSlipAccumulation` tolerance window | SSI:623-628 | `[0.2, 5.0]` × V·dt | **25× window, sign-blind — see R5-002**. |
| `TestConvergenceUnderDtHalving` | SSI:482-491 | `<1e-12` on self-diff | **Test broken — see R5-001**. |
| `TestSubStepAccumulatorSelfConsistency` | SSI:201-209 | `<1e-12` rel | OK. Tight enough — calls production helpers in reference replay (R4-008 acknowledged). |
| `TestO1LimitMatchesEvaluateADER` | SSI:259-291 | `==` bit-identity | OK. Exact integer equality on double bit patterns. |
| `TestNucleationInjectionLockedFault` | SSI:354 | `<1e-8` rel | OK. Accumulated roundoff over 1000 sub-steps is ~1e-12; 1e-8 is 4 decades of safety. |
| `TestPsiUpdatePlacement` | SSI:407 | `<1e-14` | OK at rest state. Brittle under future fixtures — see R5-005. |
| `TestPsiUpdatePlacement` bit-match | SSI:424 | `<1e-12` rel | OK. |
| `TestAccumulatorScaleIdentity` | SSI:666 | `<1e-14` rel | Tight. Scalar identity — no algorithmic bugs possible past rounding. |
| `TestSubStepAccumulatorSelfConsistency` ψ | SSI:207 | `<1e-12` | OK. |
| `TestSlipAccumulation` slip1 stays 0 | SSI:632 | `<1e-20` abs | OK. |
| `TestSlipAccumulation` slip2 window | SSI:624-628 | `[0.2, 5.0]× V·dt` | **LOOSE — see R5-002**. |
| R-006 bit-identity (R2-004) | FS_FS_newton:362 | `==` | Tight. Bit-identical iterate sequence. |
| `TestPhysicalGuessBounded` | FS_newton:381 | `≤ τ/η + 1e-12` | OK. |
| `T_TPV104_NUC_2` telescoping | NUC:108 | `<1e-6` rel | Appropriate. |
| `T_TPV104_NUC_5` telescoping | NUC:237-238 | `<1e-6` rel | Appropriate. |

Two tests have tolerances that materially reduce their bug-catching
power: `TestSlipAccumulation` (R5-002) and `TestConvergenceUnderDtHalving`
(R5-001). Both are flagged as findings above.

## Summary

- Critical issues: **2** (R5-001, R5-003 — **R5-003 is a standing
  plan violation open across four review rounds**)
- Moderate issues: 2 (R5-002, R5-004)
- Low issues: 3 (R5-005, R5-006, R5-007)
- Plan compliance: **INCOMPLETE** — T_TPV104_SSI_2
  (convergence-under-dt-halving) is a no-op test (R5-001). Plan §4.10
  Step 5 "Newton solver as default" is UNMET across rounds 1, 4, 5
  (R5-003). The iterator currently defaults to Brent; the Step-5
  SeisSol-aligned stable-asinh Newton (`SolveSlipRateNewtonStable` +
  `FrictionCoefficientStable`) is not reachable via
  `FrictionSolver::Solve`.
- Verdict: FAIL — R5-003 must be closed before Step 9 ships.
  `Method::NewtonRaphsonStable` enum value + dispatch path is the
  smallest compliant fix (~15 lines). R5-001 and R5-002 must be
  closed to restore T_TPV104_SSI_2 and T_TPV104_SSI_SLIP gates.

## Unreviewed Areas

- `seas_tpv104_driver` binary — not yet present; smoke test SKIPs.
- Behaviour of `TPV104SurfaceStationWriter::Open` on `ParMesh`
  (parallel/MPI) — `mesh.FindPoints` may silently miss parallel-side
  points outside the local rank's sub-mesh.
- Interaction between R5-004 sign-inversion and the Phase-3 probe
  output — needs Frontera runs to characterise.
- Multi-QP tests for the iterator (all current tests are single-QP).

---

## Validation checklist — Round 5

### Required closures before Step 9 driver wiring

- [ ] **R5-001 landed**: `TestConvergenceUnderDtHalving` rewritten to
  compare the imposed-state output at two different dt_macro values
  (ratio ≈ 0.5 ± tolerance, not self-diff). Gate T_TPV104_SSI_2
  actually exercises the iterator.
- [ ] **R5-002 landed**: `TestSlipAccumulation` adds a directional
  `std::signbit(Δslip2) == std::signbit(V2_final)` assertion, and
  tightens the magnitude window from `[0.2, 5.0]×` to `[0.8, 1.2]×`
  for the rest-state-adjacent fixture.
- [ ] **R5-003 resolved**: Either (a) add
  `Method::NewtonRaphsonStable` enum value dispatching to
  `SolveSlipRateNewtonStable`, OR (b) update
  `test_tpv104_smoke.cpp::TestBannerDefaults` to expect
  `"Friction solver: Brent"` on default flags and remove the
  forbidden-Brent check. Do not ship Step 9 until one of these is
  chosen.

### Recommended closures before Phase 3.B Frontera runs

- [ ] **R5-004**: add the sign-inversion diagnostic in
  `Tpv104SubStepIterator::Advance` under `SEAS_DIAG_TPV104_STATE` so
  nucleation-onset V2 oscillation is observable in probe output.
- [ ] **R5-005**: document the trajectory-free inline-reference
  assumption in `TestPsiUpdatePlacement` so the test's brittleness is
  explicit.
- [ ] **R5-006**: add a `CloseAllProbeFiles()` API and document the
  driver-shutdown ordering contract.
- [ ] **R5-007**: add a runtime guard `O ≤ 1000` on the Σ deltaT
  drift tolerance.

### Rolled forward from prior rounds

- [ ] R4-007 (station-writer MPI tie-break tolerance + no-owner abort)
  remains OPEN.
- [ ] Step 9 driver wiring uses `Rate_SRW` / per-QP V_w[i] and
  `SetProductionMode()` (round-1 deferral).
- [ ] Step 9 smoke test exercises `InitializeFaultDOFs_TPV104`
  end-to-end (closes R3-003's `T_TPV104_SIGN_2_init` gate).
- [ ] Step 13 probe-diff tool normalises reference → MFEM
  normal-stress sign flip (closes `T_TPV104_SIGN_4`).
- [ ] Phase 3 Probe 2 / Probe 4 threshold recalibration against
  round-1 R-004, R-005, R-006 closures.

### Cross-round invariant anchors (updated)

- σ_n > 0 = compression — `T_TPV104_SIGN_1`.
- `SlipLawSRWPsi` production mode aborts on base-virtual call —
  R-001 + T_SRW_7.
- Raw V through state-evolution — R-002 + TestNoVsafeClamp.
- `(V/V_w)^8` unrolled integer power — R-004 + R2-002.
- Newton input validation — R2-004 + TestNewtonRejectsInvalidVPrev.
- Nucleation accumulator telescopes to Δτ₀·F(r) — T_TPV104_NUC_2/5.
- TPV104 pure strike-slip: `tau1_nuc = sigma_n_nuc = 0` after
  nucleation — `TestNoDipOrNormalStressNucleation` + 2 inline guards.
- TPV102 overwrite pattern unaffected by Step 6 — T_TPV104_NUC_4.
- TPV104 mesh byte-identical to TPV102 — T_TPV104_MESH_0/0b.
- TPV104 slip accumulation is live in the iterator — R4-001 closed.
  **NEEDS R5-002 CLOSURE** for sign-error robustness.
- TPV104 Probe-3 μ matches solver-intended μ — R4-003 closed.
- TPV104 probe output per-rank under MPI — R4-005 closed.
- **NEW — T_TPV104_SSI_2 convergence gate actually exercises the
  iterator** — **NEEDS R5-001 CLOSURE**.
- **NEW — Plan §4.10 Step 5 (Newton-stable default) matches
  iterator default AND smoke-test banner expectation** — **NEEDS
  R5-003 CLOSURE**.
