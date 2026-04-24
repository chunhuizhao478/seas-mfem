# Code Review Round 7: Step-9 driver mesh-coupled extension — 2026-04-24

Fresh adversarial audit of the Step-9 driver extension landed since
round 6. The driver grew from ~254 lines to 689 lines (+435 lines) —
roughly triple in size. Primary focus this round: whether the
mesh-coupled production path matches what the CLI / banner / tests
claim it does.

## Review Scope

- Plan: `tpv104_debug_plan_2026-04-24.md` §4.10 Step 9 (driver) +
  §4.10.X (R5-003 fix).
- Prior reviews: R-001..R-012, R2-001..R2-007, R3-001..R3-007,
  R4-001..R4-008, R5-001..R5-007, R6-001..R6-004.
- **Primary code under review** (changes since round 6):
  - `drivers/tpv104_driver.cpp` (254 → 689 lines, +171 % — full
    mesh-coupled ADER time loop, station writers, MPI ParMesh
    partitioning).
  - `tests/unit/test_tpv104_substep_iterator.cpp` (701 → 884 lines —
    R5-001 rewrite + R5-002 sign guard).
  - `tests/unit/test_tpv104_smoke.cpp` (193 → 202 lines — minor).
  - `tests/unit/test_tpv104_probe_format.cpp` (234 → 240 lines —
    minor).
  - `jobs/tpv104/tpv104_*.sbatch` — stale-header review.

## Prior-round closure verification

| ID | Prior severity | Status | Notes |
|---|---|---|---|
| R5-001 | CRITICAL | **CLOSED** | `test_tpv104_substep_iterator.cpp:460-534` rewrote the test to compare Ip_big vs 2·Ip_half; no more self-diff. |
| R5-002 | MODERATE | **CLOSED** | `:707-733` adds signbit-directional + tight `[0.8, 1.2]×` window. |
| R5-003 | CRITICAL (plan deviation) | **CLOSED on enum/dispatch**; **NOT CLOSED on time-loop wiring** — see R7-001. |
| R6-001 | MODERATE | **STILL OPEN** — `drivers/tpv104_driver.cpp:109-110` still aliases bare `newton` → legacy. |
| R6-002 | MODERATE | **STILL OPEN** — `:112` still silently returns `NewtonRaphsonStable` on unknown strings. |
| R4-007 | LOW | **STILL OPEN** — station-writer MPI tie-break tolerance unchanged. |

## Findings

### [R7-001] CRITICAL [drivers/tpv104_driver.cpp:657-663, 574] — `--friction-solver` and `--fault-iterator` CLI flags are banner-only; the ADER time loop hard-codes Brent via `wave.AdvanceADER` → `EvaluateADERTotal` and never invokes the sub-step iterator

**Category:** BUG (silent disclosure failure — critical plan deviation)

**Description:**
The driver dutifully parses `--friction-solver newton-stable` and
`--fault-iterator substep`, dispatches `MapSolver(...)` to produce
`FrictionSolver::Method::NewtonRaphsonStable`, and prints the banner
line `"Friction solver: Newton-Raphson (stable-asinh, plan §4.10
Step 5)"`. Then in the time loop (`:548-664`):

```cpp
for (int step = 0; step < nsteps; ++step)
{
   ...
   ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                     t + dt_step, dt_step);
   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);   // line 574
   Q.Swap(Q_new);
   ...
   for (int i = 0; i < num_fault_total; ++i)
   {
      dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
         psi_n[i], dof_data[i].slip_rate, ...);
      dof_data[i].slip1 += dof_data[i].V1 * dt_step;
      dof_data[i].slip2 += dof_data[i].V2 * dt_step;
   }
   ...
   (void)method;   // line 663
}
```

The `method` variable computed from `--friction-solver` is **cast to
void** — it is never passed to any solver. The time loop uses
`wave.AdvanceADER(...)` which internally dispatches through
`FaultFaceFlux::EvaluateADERTotal` → `FrictionSolver::Solve`
with its default argument value `Method::Brent`
(`friction_solver.hpp:72`). So **every `--friction-solver` value
produces the same underlying solver: Brent.**

Same pattern for `--fault-iterator substep`: the driver parses it at
line 194-195, comment at line 196-198 explicitly acknowledges "`fric_law`
and `fault_iterator` are banner-only today". There is no
`wave.SetFaultIterator(&iterator)` call (that method does not exist
in `WaveOperator`), and no `Tpv104SubStepIterator::Advance` call in
the time loop. The sub-step iterator is never invoked; the
production time loop is the one-shot `wave.AdvanceADER` pattern
inherited from TPV102.

User-visible consequences when the current driver runs on Frontera:

1. Banner prints `"Friction solver: Newton-Raphson (stable-asinh,
   plan §4.10 Step 5)"` → actually runs **Brent**.
2. Banner prints `"Fault iterator: sub-step"` → actually runs
   **one-shot** (TPV102 pattern, not Plan §3.8 per-sub-step
   accumulation).
3. Banner prints `"Friction law: slip-SRW (ψ-space)"` → this one is
   accurate because ψ is updated by the `UpdateStateAnalyticSlipLawSRW`
   FREE function (line 584). OK.
4. Plan §3.12 per-sub-step ψ evolution — **NOT** what happens. ψ is
   updated ONCE per macro-step with `dof_data[i].slip_rate`
   (one-shot ADER-averaged V), not per-sub-step V.
5. Plan §3.9 per-sub-step cumulative nucleation —
   `ApplyNucleationIncremental_TPV104` is called ONCE per macro-step
   at line 565-566 with `dt_step = dt_macro`. The increment
   telescopes correctly at macro-step boundaries, so the FINAL
   tau2_nuc is right. But intermediate sub-step tau2_nuc values
   never exist — there are no sub-steps.

This is the exact disclosure failure R5-003 / R6-002 were filed to
prevent. The fix report from the round-6 fix round said:

> "Iterator reverted to uniform `ComputeStageState(..., method)`
> dispatch; removed the `use_stable_newton_` fork."

That statement is true at the **iterator layer**, but the **driver
never calls the iterator**. The iterator's uniform dispatch is
unreachable through the production time loop; the only code path
that reaches it is the unit-test suite.

Plan §4.10 Step 9 specifies (CLI wiring + iterator binding):
> "`Tpv104SubStepIterator::Advance` is the production path;
> `wave.AdvanceADER` is retained only for the `--fault-iterator
> oneshot` regression test."

The current driver inverts this: `wave.AdvanceADER` IS the production
path; `Tpv104SubStepIterator` is unreachable.

**Trigger:** Any Frontera submission of `tpv104_init_stations.sbatch`
or `tpv104_t3_mesh_coupled.sbatch`.

**Actual behavior:** Banner says one thing, time loop runs another.
Brent + one-shot regardless of CLI. Phase-3 Probe 3/4 comparisons
against SeisSol would use Brent while analyst thinks stable-asinh
Newton is running.

**Expected behavior:** One of two options:

- **(a) Wire the iterator into the time loop** (plan-compliant):
  replace `wave.AdvanceADER(...)` with an ADER predictor call + a
  `Tpv104SubStepIterator::Advance(...)` call, passing the `method`
  from `MapSolver`. This requires adding a mechanism to obtain the
  ADER time-integrated I_± from the wave operator (currently
  private inside `wave.AdvanceADER`).

- **(b) Downgrade the banner to match reality**: replace the four
  banner lines with `"Time integrator: ADER-O2 (one-shot via
  wave.AdvanceADER; sub-step iterator deferred)"` and `"Friction
  solver: Brent (hard-coded via EvaluateADERTotal; --friction-solver
  flag ignored)"`. This honestly discloses the current state. Plan
  §4.10 Step 9 acceptance is NOT met under this option.

**Suggested fix (option b as the minimum-disruption stop-gap until
option a lands):**
```diff
    if (rank == 0)
    {
       std::cout << "========================================\n";
       std::cout << "SCEC TPV104 Dynamic Rupture Simulation\n";
       ...
-      std::cout << "Time integrator: ADER-O" << ader_order << "\n";
-      // Canonical banner strings matched by test_tpv104_smoke:
-      //   "sub-step" / "one-shot" (hyphenated),
-      //   "slip-SRW" / "aging" (title-case law),
-      //   SolverBanner already canonical.
-      if (fault_iterator == "substep") {
-         std::cout << "Fault iterator: sub-step\n";
-      } else if (fault_iterator == "oneshot") {
-         std::cout << "Fault iterator: one-shot (regression-only; "
-                      "NOT the production path)\n";
-      } else {
-         std::cout << "Fault iterator: " << fault_iterator
-                   << " (unknown)\n";
-      }
-      std::cout << "Friction solver: " << SolverBanner(friction_solver) << "\n";
+      std::cout << "Time integrator: ADER-O" << ader_order
+                << " (one-shot via wave.AdvanceADER)\n";
+      // R7-001: the sub-step iterator (Step 7) is NOT wired into the
+      // production time loop yet.  The --fault-iterator / --friction-
+      // solver CLI flags are accepted for smoke-test banner parity
+      // but have no effect on the code path.  Until the iterator
+      // binding lands, the time loop hard-codes Brent via
+      // EvaluateADERTotal.
+      std::cout << "Fault iterator: one-shot (sub-step iterator NOT "
+                   "wired; --fault-iterator flag IGNORED)\n";
+      std::cout << "Friction solver: Brent (hard-coded via "
+                   "EvaluateADERTotal; --friction-solver flag IGNORED)\n";
+      std::cout << "CLI parsed values (banner-only): "
+                << "friction_solver=" << friction_solver
+                << ", fault_iterator=" << fault_iterator << "\n";
```
And remove the void cast since `method` is unused:
```diff
-      // Method argument is currently read by the sub-step iterator only
-      // (not by wave.AdvanceADER's internal FrictionSolver::Solve
-      // dispatch, which uses Brent via FaultFaceFlux::EvaluateADERTotal).
-      // The friction-solver CLI flag is therefore banner-only on this
-      // path; when the Tpv104SubStepIterator replaces AdvanceADER (a
-      // follow-up integration), `method` drives the per-sub-step solve.
-      (void)method;
```

And update the smoke-test `TestBannerDefaults`
(`test_tpv104_smoke.cpp:87-94`) accordingly:
```diff
-      {"Friction solver: Newton-Raphson (stable-asinh, plan §4.10 Step 5)",
-       "stable-asinh Newton default (plan §4.10.X)"},
+      {"Friction solver: Brent (hard-coded",
+       "honest banner: Brent under one-shot wave.AdvanceADER"},
-      {"Fault iterator: sub-step", "sub-step iterator default"},
+      {"Fault iterator: one-shot (sub-step iterator NOT wired",
+       "honest banner: sub-step iterator is unreachable"},
```

**Test case:**
```cpp
void test_R7_001_banner_matches_actual_solver() {
   // Invariant: the solver the driver ACTUALLY runs must match what
   // the banner claims.  Under option (b) the banner says Brent and
   // runs Brent; under option (a) the banner says NewtonRaphsonStable
   // and the iterator dispatches it.  Either way, the two must agree.
   //
   // Proof by dispatch: instrument FrictionSolver::Solve to record
   // which Method it saw on the last call.  After driver finishes a
   // single macro-step, read the recorded Method; assert it equals
   // the Method derived from the banner string.
   //
   // Currently: banner says NewtonRaphsonStable, dispatch recorded
   //   Brent → FAIL.
   const std::string out = RunDriver(binary,
      "--dry-run --friction-solver newton-stable");
   const FrictionSolver::Method banner_method = ParseMethodFromBanner(out);
   const FrictionSolver::Method actual_method = ReadInstrumentedDispatch();
   EXPECT_EQ(banner_method, actual_method);
}
```

---

### [R7-002] MODERATE [drivers/tpv104_driver.cpp:469-472] — `SlipLawSRWPsi law` is constructed and `SetProductionMode()` called, but `law` is never dispatched — the R-001 safety net does not protect the driver path

**Category:** QUALITY (dead code suggesting a guard that is not live)

**Description:**
```cpp
SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
                  TPV104Params::V0, TPV104Params::f0,
                  TPV104Params::f_w, TPV104Params::V_w_in);
law.SetProductionMode();
```
`law` is constructed, flipped into production mode, and then **never
called on any code path in the driver**. The ψ update at line 584
uses the free function `UpdateStateAnalyticSlipLawSRW(...)`, not
`law.Rate()` or `law.SteadyState()` or any other method on `law`.

Round-1 R-001 introduced `SlipLawSRWPsi::SetProductionMode()`
specifically to guard against silent fallthrough when the driver
dispatches through the `StateEvolution` base virtuals. But this
driver **bypasses the base class entirely** — it calls the free
function with per-QP (V_w, a) explicitly at every QP. The
production-mode guard is dead weight here.

If a future maintainer adds `law.Rate(...)` or `law.SteadyState(...)`
to the driver expecting the abort-on-misuse safety net to fire, the
guard would work. But it's misleading that the driver currently
carries an object + production-mode flip that never protects the
running code.

**Trigger:** Any reader / maintainer inspecting the driver, seeing
`law.SetProductionMode()`, and assuming the state-evolution dispatch
is through `law`.

**Actual behavior:** `law` is dead code; the production-mode flag
has no effect on this driver.

**Expected behavior:** Either (a) remove `law` and its
`SetProductionMode()` call (since the free-function path doesn't need
them), OR (b) bind `law` into a future iterator wiring so the guard
is actually load-bearing.

**Suggested fix (option a — remove dead code):**
```diff
-   // SlipLawSRWPsi — global scalars + required V_w_default stub.  Per-QP
-   // V_w comes from V_w[] above and feeds into UpdateStateAnalyticSlipLawSRW.
-   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
-                     TPV104Params::V0, TPV104Params::f0,
-                     TPV104Params::f_w, TPV104Params::V_w_in);
-   law.SetProductionMode();
+   // State-evolution is invoked via the free function
+   // UpdateStateAnalyticSlipLawSRW(...) at line 584 with per-QP
+   // (V_w, a).  The SlipLawSRWPsi class virtuals are not dispatched
+   // on this driver path, so no SlipLawSRWPsi instance is needed.
+   // If a future R7-001 wiring routes through the iterator, construct
+   // SlipLawSRWPsi and call SetProductionMode() at that point.
```

**Test case:** Not required (no code is broken today; this is a
clarity / dead-code concern).

---

### [R7-003] MODERATE [jobs/tpv104/tpv104_init_stations.sbatch:19-21] — sbatch header claims `--friction-solver newton-stable` is the TPV104 default and `--fault-iterator substep` is used, but those flags have no effect (R7-001)

**Category:** DEVIATION (documentation lies to the user)

**Description:**
`jobs/tpv104/tpv104_init_stations.sbatch:19-21`:
```
# Mesh: tpv104_200m.msh (production canonical mesh per plan §4.1).
# Sub-step iterator: substep (default).
# Friction solver: newton-stable (plan §4.10.X default).
# Time integrator: ADER-O2.
```
And `jobs/tpv104/tpv104_t3_mesh_coupled.sbatch:19-25`:
```
# Sub-step iterator: substep (plan §3.8 per-sub-step accumulator)
# Friction solver  : newton-stable (plan §4.10.X R5-003)
```
Both sbatch headers claim the plan-canonical settings run. Under
R7-001, those flags are banner-only — the actual run uses Brent +
one-shot `wave.AdvanceADER`. A Frontera operator reading the sbatch
header to understand what solver was used will be misled.

Further: the `tpv104_build_banner_sanity.sbatch` PASS criteria
(round-6 fix report) requires the banner to show `"Friction solver:
Newton-Raphson (stable-asinh, plan §4.10 Step 5)"`. Under R7-001 the
banner DOES show that string — but the underlying physics run uses
Brent. The sanity check passes under the literal banner-string match,
missing the deeper banner-vs-code divergence.

**Trigger:** Any Frontera run of any `tpv104_*.sbatch`.

**Actual behavior:** sbatch header comments and banner output lie
about the actual solver + iterator.

**Expected behavior:** sbatch headers should accurately describe the
production-path behavior. If R7-001 is closed by wiring the iterator,
update `build_banner_sanity` gate to assert the dispatched solver
matches the banner (not just the banner-string literal).

**Suggested fix:** Under R7-001 option (b):
```diff
 # Mesh: tpv104_200m.msh (production canonical mesh per plan §4.1).
-# Sub-step iterator: substep (default).
-# Friction solver: newton-stable (plan §4.10.X default).
-# Time integrator: ADER-O2.
+# Time integrator: ADER-O2 (one-shot via wave.AdvanceADER).
+# Fault iterator : one-shot (sub-step iterator NOT wired; R7-001
+#                  open).  --fault-iterator flag is banner-only.
+# Friction solver: Brent (hard-coded via EvaluateADERTotal).  The
+#                  --friction-solver=newton-stable flag below is
+#                  parsed for banner/smoke-test parity but has no
+#                  effect on the dispatched solver until R7-001 is
+#                  closed.  See jobs/tpv104/README.md for timeline.
```
Apply the same correction to the production sbatch header.

**Test case:**
```python
def test_R7_003_sbatch_headers_honest():
    # After R7-001 + R7-003 fix, grep every sbatch header for the
    # word "newton-stable" outside of a "banner-only" or "open" caveat.
    # No non-caveat occurrences → PASS.
    for f in glob("jobs/tpv104/*.sbatch"):
        for line in open(f):
            if "newton-stable" in line and "banner-only" not in line \
                                        and "R7-001" not in line:
                assert False, f"stale header in {f}: {line}"
```

---

### [R7-004] MODERATE [test_tpv104_smoke.cpp:87-94] — Smoke test `TestBannerDefaults` asserts the banner CLAIMS a solver but never verifies the driver actually USES that solver

**Category:** BUG (test does not validate the invariant it claims to validate)

**Description:**
```cpp
const std::vector<std::pair<std::string, std::string>> must_have = {
   {"Time integrator: ADER-O2", "ADER-O2 default"},
   {"Fault iterator: sub-step", "sub-step iterator default"},
   {"Friction solver: Newton-Raphson (stable-asinh, plan §4.10 Step 5)",
    "stable-asinh Newton default (plan §4.10.X)"},
   {"Friction law: slip-SRW (ψ-space)", "slip-SRW friction-law default"},
};
for (const auto &kv : must_have)
{
   const bool found = out.find(kv.first) != std::string::npos;
   TEST_ASSERT(found, ...);
}
```
The test only greps the banner text for literal strings. It does NOT
verify that:
- The driver's dispatched `FrictionSolver::Method` equals what the
  banner claims.
- The driver's time loop invokes the `Tpv104SubStepIterator` (banner
  claims "sub-step" but the driver doesn't).
- The ψ-update cadence is per-sub-step (banner implies via "sub-step")
  but the driver does macro-step-only.

This enabled R7-001 to pass 16/16 smoke tests while the production
code path is completely at odds with the banner.

**Trigger:** Any Step-9 driver that prints the right banner strings
but dispatches differently.

**Actual behavior:** Test passes on literal string match; silently
accepts banner-vs-code divergence.

**Expected behavior:** Add a dispatch-verification test. Either via
a driver `--dry-run` mode that instruments `FrictionSolver::Solve`
to print the actual method dispatched, or via a public driver-state
query API.

**Suggested fix:** Add a diagnostic `--verify-dispatch` CLI flag to
the driver:
```cpp
// drivers/tpv104_driver.cpp
if (HasFlag(argc, argv, "--verify-dispatch"))
{
   // Print, as the last line of banner output:
   //   [dispatch] friction_solver_actual=brent
   //   [dispatch] fault_iterator_actual=oneshot
   //   [dispatch] friction_law_actual=slip-srw
   // Derived from the ACTUAL code path run during the dry-run step,
   // not from the CLI values.  When the iterator is wired, update
   // these strings accordingly.
   std::cout << "[dispatch] friction_solver_actual=brent\n";
   std::cout << "[dispatch] fault_iterator_actual=oneshot\n";
   std::cout << "[dispatch] friction_law_actual=slip-srw\n";
}
```
Then extend `TestBannerDefaults` to assert dispatch ≡ banner:
```cpp
const std::string out = RunDriver(binary,
   "--dry-run --verify-dispatch");
// Assert the dispatch line matches the banner line.
const bool banner_claims_stable =
   out.find("Newton-Raphson (stable-asinh") != std::string::npos;
const bool dispatch_stable =
   out.find("friction_solver_actual=newton-stable") != std::string::npos;
TEST_ASSERT(banner_claims_stable == dispatch_stable,
            "banner and dispatch agree on friction solver");
```

---

### [R7-005] MODERATE [drivers/tpv104_driver.cpp:109-113] — R6-001 and R6-002 still open: bare `newton` → legacy, unknown strings silent-default

**Category:** DEVIATION (round-6 finding not addressed)

**Description:**
```cpp
static FrictionSolver::Method MapSolver(const std::string &s)
{
   if (s == "newton-stable") { return FrictionSolver::Method::NewtonRaphsonStable; }
   if (s == "brent")         { return FrictionSolver::Method::Brent; }
   if (s == "newton-legacy"
       || s == "newton")     { return FrictionSolver::Method::NewtonRaphson; }
   if (s == "hybrid")        { return FrictionSolver::Method::HybridNRBisection; }
   return FrictionSolver::Method::NewtonRaphsonStable;
}
```
R6-001 (bare `newton` should route to stable-asinh, not legacy) and
R6-002 (unknown strings should abort loudly, not silent-default) are
both still open.

Under R7-001 these don't matter because the `method` return value is
discarded anyway. But once R7-001 is closed (iterator wired), R6-001
and R6-002 become active bugs again.

**Trigger:** Any Step-9 driver invocation with `--friction-solver
newton` or `--friction-solver <typo>` after R7-001 is closed.

**Actual behavior:** Bare `newton` → legacy; typos → silent stable.

**Expected behavior:** Bare `newton` → stable (canonical shorthand);
typos → abort with usage message.

**Suggested fix:** See R6-001 and R6-002 suggestions in
`tpv104_review_round6_2026-04-24.md`. Unchanged.

**Test case:** Identical to R6-001 and R6-002 test cases.

---

### [R7-006] LOW [drivers/tpv104_driver.cpp:196-198] — Comment claims "sub-step-aware nucleation" when the code path is macro-step-only nucleation

**Category:** QUALITY (stale comment)

**Description:**
```cpp
// `fric_law` and `fault_iterator` are banner-only today; the code
// path is FVW/slip-SRW + sub-step-aware nucleation already.  A
// future --fric-law=aging would switch the ψ update function.
```
The comment says "sub-step-aware nucleation already" — but
`ApplyNucleationIncremental_TPV104` is called ONCE per macro-step at
line 565-566 with `dt_step = dt_macro`. There are no sub-steps in
this driver. The telescoping identity makes the total perturbation
reach the right value at macro-step boundaries, so the physics is
correct at those boundaries, but "sub-step-aware" implies per-sub-step
evaluation — which doesn't happen.

**Trigger:** Any reader inspecting this comment.

**Actual behavior:** Comment overstates the code's fidelity to Plan
§3.9.

**Expected behavior:** Comment accurately describes the macro-step
cadence.

**Suggested fix:**
```diff
-   // `fric_law` and `fault_iterator` are banner-only today; the code
-   // path is FVW/slip-SRW + sub-step-aware nucleation already.  A
-   // future --fric-law=aging would switch the ψ update function.
+   // `fric_law` and `fault_iterator` are banner-only today (R7-001).
+   // The code path is FVW/slip-SRW with MACRO-step-only nucleation —
+   // ApplyNucleationIncremental_TPV104 runs once per macro-step with
+   // dt_step = dt_macro.  The telescoping smoothStep identity means
+   // the total perturbation at t = T_nuc is still Δτ·F(r) (correct),
+   // but plan §3.9's per-SUB-step cadence is not met until the
+   // sub-step iterator is wired.  A future --fric-law=aging would
+   // switch the ψ update function.
```

---

### [R7-007] LOW [drivers/tpv104_driver.cpp:584-594] — ψ-update uses `dof_data[i].slip_rate` (macro-step V) instead of per-sub-step V; plan §3.12 deviation not surfaced in any diagnostic

**Category:** DEVIATION (plan §3.12 not met; doc-only issue)

**Description:**
Plan §3.12 directive: "follow exactly SeisSol did — per-sub-step ψ
integration". The current driver:
```cpp
for (int i = 0; i < num_fault_total; ++i)
{
   dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
      psi_n[i], dof_data[i].slip_rate, dof_data[i].Dc, dt_step,
      V_w[i], dof_data[i].a,
      TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::f_w);
   ...
}
```
ψ is integrated in a single `UpdateStateAnalyticSlipLawSRW` call per
macro-step with:
- `V = dof_data[i].slip_rate` (ADER-one-shot V ≈ time-average).
- `dt = dt_step` (macro-step size).

This is the TPV102 pattern. Plan §3.12 requires O per-sub-step calls
with `(V_substep_o, dt_sub_o)`. Same formula as R7-001 acknowledges;
same class of deviation. Tracked separately because the ψ cadence
has its own correctness implications: under a rapidly-changing V the
macro-step ψ deviates by O(dt_macro²) from the per-sub-step result.

**Trigger:** Nucleation-breakout or rupture-front passage where V
changes rapidly within a macro-step.

**Actual behavior:** ψ evolves with macro-step-average V; per-sub-step
ψ is never computed.

**Expected behavior:** Closed when R7-001 is closed (wiring the
iterator gives per-sub-step ψ).

**Suggested fix:** Subsumed by R7-001 option (a).

---

## Tolerance audit update (round 7)

Confirmed R5-001 (TestConvergenceUnderDtHalving) now has a real
check: Ip_big vs 2·Ip_half at 1e-10 relative. Proper convergence
guard.

Confirmed R5-002 (TestSlipAccumulation) now has a signbit directional
assertion and tight [0.8, 1.2]× magnitude window. Sign errors would
fire the test.

No new tolerance issues identified in the round-7 deltas.

## Summary

- Critical issues: 1 (R7-001 — banner-vs-code disclosure failure in
  the production time loop)
- Moderate issues: 4 (R7-002 dead SlipLawSRWPsi object, R7-003 stale
  sbatch headers, R7-004 test doesn't verify banner-vs-code,
  R7-005 = R6-001/R6-002 still open)
- Low issues: 2 (R7-006 stale "sub-step-aware nucleation" comment,
  R7-007 plan §3.12 per-sub-step ψ deviation)
- Plan compliance: **INCOMPLETE — Step 9 claims iterator wiring in
  §4.10 Step 9; the driver does NOT wire it.** Banner claims match
  the plan but code path doesn't.
- Verdict: FAIL — R7-001 must be resolved before ANY Frontera
  submission of `tpv104_init_stations.sbatch` or
  `tpv104_t3_mesh_coupled.sbatch`. The production runs as currently
  configured will produce Brent-solver results while the banner and
  sbatch headers claim stable-asinh Newton — any subsequent Phase-3
  probe-diff analysis against SeisSol would blame the wrong solver
  and the diagnostic effort would waste Frontera SU.
  `tpv104_build_banner_sanity.sbatch` and `tpv104_mesh_build.sbatch`
  remain safe to submit as long as the user understands the banner is
  misleading.

## Unreviewed Areas

- `WaveOperator::AdvanceADER` internals on ParMesh — depends on
  shared fault faces across MPI ranks; not reproducible locally per
  `feedback_no_local_reproducer`. R7-001 option (a) needs this to
  expose per-sub-step I_± to the iterator.
- The Step-9 driver has now shipped a lot of new MPI code (ParMesh
  partitioning, station-writer MPI ownership, NaN tripwires,
  verification against `WaveOperator::VerifySharedFaultDOFDataConsistency`).
  These paths are best reviewed via a short np=4 local smoke test
  on a tiny mesh — deferred to Step-10 mesh generation.

---

## Validation checklist — Round 7

### Required closures before ANY Frontera sbatch submission beyond build-sanity + mesh-build

- [ ] **R7-001 CRITICAL**: resolve banner-vs-code disclosure
  inconsistency. Either (a) wire `Tpv104SubStepIterator::Advance`
  into the driver time loop with the parsed `method`, OR (b)
  downgrade the banner to accurately describe the current one-shot
  Brent path.
- [ ] **R7-003**: sbatch headers updated to match R7-001's resolution.
- [ ] **R7-004**: smoke test adds a banner-vs-dispatch consistency
  check (`--verify-dispatch` + instrumented `FrictionSolver::Solve`).

### Required before Phase 3.B Frontera probe runs (if option (b) chosen)

- [ ] Plan §4.10 Step 9 acceptance re-evaluated. Under option (b)
  (iterator not wired), the Step-9 deliverable is INCOMPLETE; the
  iterator work from Step 7 is wasted effort until wired. This must
  be decision by user.

### Required before Phase 3.B Frontera probe runs (if option (a) chosen)

- [ ] R6-001 / R6-002 (R7-005) closed — iterator dispatch through
  `FrictionSolver::Solve(method)` becomes live, so CLI alias
  consistency is load-bearing again.
- [ ] R7-007 closed automatically by iterator wiring (per-sub-step ψ).
- [ ] Multi-QP MPI test for the iterator — deferred to Frontera
  (Script 3 PASS).

### Recommended closures before Phase 3.B probe interpretation

- [ ] R7-002: remove the dead `SlipLawSRWPsi law` + `SetProductionMode()`
  construction from the driver, OR bind it to a live iterator.
- [ ] R7-006: update the "sub-step-aware nucleation" comment to match
  the actual macro-step cadence.
- [ ] R4-007: station-writer MPI tie-break tolerance.

### Rolled forward from prior rounds

- [ ] Step 13 probe-diff tool normalises reference → MFEM
  normal-stress sign flip (closes `T_TPV104_SIGN_4`).
- [ ] Phase 3 probe threshold recalibration.

### Cross-round invariant anchors (updated)

- σ_n > 0 = compression — `T_TPV104_SIGN_1`.
- `SlipLawSRWPsi` production-mode aborts — R-001 + T_SRW_7. **Note:
  guard is live in unit tests but NOT in the driver — R7-002.**
- Raw V through state-evolution — R-002 + TestNoVsafeClamp.
- `(V/V_w)^8` unrolled integer power — R-004 + R2-002.
- Newton input validation — R2-004.
- Nucleation accumulator telescopes — T_TPV104_NUC_2/5.
- TPV104 strike-slip invariant — R3-001.
- TPV102 overwrite-pattern nucleation unaffected — T_TPV104_NUC_4.
- TPV104 mesh byte-identical to TPV102 — T_TPV104_MESH_0/0b.
- Slip accumulation in iterator (R4-001) — **live in iterator AND
  driver** (driver duplicates the loop at `:595-596`; iterator still
  unreachable from driver path).
- Probe-3 μ via `FrictionCoefficientStable` — R4-003. **Note:
  iterator probe-3 is unreachable; the driver has no probe output.**
- Probe output per-rank under MPI — R4-005. **Unreachable from driver.**
- Plan §4.10 Step 5 canonical solver reachable — R5-003 CLOSED at
  enum/dispatch layer, **NOT CLOSED at driver time-loop** —
  **NEEDS R7-001 CLOSURE**.
- CLI / banner / solver tri-consistency — banner matches solver on
  valid strings ✅ at enum layer, **NOT at driver time-loop** —
  **NEEDS R7-001 CLOSURE**.
- T_TPV104_SSI_2 convergence gate now exercises the iterator — R5-001
  CLOSED. Iterator tests pass but iterator is unreachable from driver.
- TestSlipAccumulation directional guard — R5-002 CLOSED.
