# Code Review Round 3: TPV104 Steps 6, 8, 10 — 2026-04-24

Fresh adversarial audit of the new code landed since round 2. Confirms
prior round-1 and round-2 fixes are still intact, then hunts for new
bugs in the files that shipped between rounds.

## Review Scope

- Plan: `miniapps/seas/debug_document/tpv104_debug_document/tpv104_debug_plan_2026-04-24.md`
  §4.10 Steps 6 (nucleation), 8 (normal-stress sign audit), 10 (mesh).
- Prior reviews (closed):
  - `tpv104_review_2026-04-24.md` (round 1, R-001..R-012 CLOSED)
  - `tpv104_review_round2_2026-04-24.md` (round 2, R2-001..R2-007 CLOSED)
- New files reviewed:
  - `dynamic/tpv104_nucleation.hpp` (129 lines — Step 6)
  - `tests/unit/test_tpv104_nucleation.cpp` (258 lines)
  - `tests/unit/test_tpv104_normal_sign.cpp` (177 lines — Step 8)
  - `tests/unit/test_tpv104_mesh.cpp` (194 lines — Step 10)
  - `tpv104/mesh/tpv104_200m.geo` (duplicated TPV102 200m, Step 10)
  - `tpv104/mesh/tpv104_500m.geo` (new 500 m variant)
  - `tpv104/mesh/tpv104_1000m.geo` (duplicated TPV102 1000m)
- Unchanged since round 2:
  - `config/tpv104_params.hpp` (250 lines)
  - `friction/slip_law_srw_psi.hpp` (315 lines, round-2 R2-003/R2-007 fixes)
  - `friction/friction_coeff_stable.hpp` (159 lines)
  - `dynamic/tpv104_friction_solver.hpp` (178 lines, round-2 R2-004 fixes)
  - All test files from Steps 1/2/4/5.
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `feedback_tpv102_bp5_no_shared_edit`.

## Round-1/2 regression check

Spot-verified: all R-001..R-012 and R2-001..R2-007 fixes remain in
place — `SlipLawSRWPsi::SetProductionMode`, raw V (no clamp),
`IntegerPow8`, required `V_w_default`, `RateDerivativeV(V=0)` analytic
limit, `SolveSlipRateNewtonStable` input validation, `physical_v_guess`
τ/η_s cap. No regressions introduced by the new files.

## Findings

### [R3-001] MODERATE [test_tpv104_nucleation.cpp:76-109, 185-246] — No guard on the TPV104 "no normal-stress nucleation" invariant (`tau1_nuc == 0, sigma_n_nuc == 0`)

**Category:** ASSUMPTION (missing test for an explicit design contract)

**Description:**
`tpv104_nucleation.hpp:120-123` documents: *"tau1_nuc and sigma_n_nuc
are intentionally not updated (pure strike-slip; no normal-stress
nucleation in TPV104)."* This is a load-bearing invariant — TPV104 is
pure strike-slip and any spurious write to the dip-shear or
normal-stress nucleation channels would perturb the friction residual
and alter Probe 1 traces in Phase 3.

But `TestAccumulatorTelescopes`, `TestGuardBounds`, and the other
Step-6 tests ONLY check `dof_data[0].tau2_nuc`. Neither `tau1_nuc` nor
`sigma_n_nuc` is ever asserted to remain at zero after the
accumulator runs. A future bug that accidentally adds an increment to
`tau1_nuc` (copy-paste error from Step 7 sub-step iterator; swapping
`tau2_nuc` → `tau1_nuc` in a refactor) would pass every unit test and
only surface as a Phase-3 probe diff.

**Trigger:** Any future refactor that touches the `dof_data[i].tau*_nuc`
assignments in `ApplyNucleationIncremental_TPV104` and accidentally
writes to the wrong channel.

**Actual behavior:** Tests pass while the strike-slip invariant is
silently broken.

**Expected behavior:** After every `ApplyNucleationIncremental_TPV104`
call in the test suite, assert `tau1_nuc == 0` and `sigma_n_nuc == 0`
bit-exactly.

**Suggested fix:** Add an invariant check at the end of each test that
calls the accumulator, and a dedicated guard test:
```diff
    const real_t expected = TPV104Params::nuc_dtau * 1.0;
    const real_t got = dof_data[0].tau2_nuc;
    const real_t rel = std::abs(got - expected) / std::abs(expected);
    std::cout << "  tau2_nuc(hypo) after T_nuc = " << got
              << " (expected " << expected << ", rel err " << rel << ")\n";
    TEST_ASSERT(rel < 1e-6,
                "tau2_nuc[hypo] ≈ Δτ₀·F(0) = Δτ₀·1 to 1e-6 rel");
+
+   // R3-001: TPV104 is pure strike-slip — tau1_nuc and sigma_n_nuc
+   // must remain at 0 after every accumulator call.
+   TEST_NEAR(dof_data[0].tau1_nuc, 0.0, 1e-30,
+             "tau1_nuc remains at 0 (no dip-shear nucleation)");
+   TEST_NEAR(dof_data[0].sigma_n_nuc, 0.0, 1e-30,
+             "sigma_n_nuc remains at 0 (no normal-stress nucleation)");
 }
```
Add the same two assertions to `TestGuardBounds` and add a dedicated
`TestNoDipOrNormalStressNucleation` test:
```cpp
void TestNoDipOrNormalStressNucleation()
{
   std::vector<DOFData> dof_data(1);
   std::vector<Vector> fault_coords(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = 0.0;
   fault_coords[0](1) = 0.0;
   fault_coords[0](2) = -TPV104Params::hypo_down_dip;
   // Populate tau1_nuc / sigma_n_nuc with sentinel values to detect
   // any spurious write by the accumulator.
   dof_data[0].tau1_nuc    = 1.0e99;
   dof_data[0].sigma_n_nuc = 2.0e99;
   dof_data[0].tau2_nuc    = 0.0;
   // Run a full ramp.
   const real_t dt = 0.01;
   for (int k = 1; k <= 200; ++k) {
      ApplyNucleationIncremental_TPV104(
         dof_data, fault_coords, k * dt, dt);
   }
   // Accumulator writes only to tau2_nuc.
   TEST_NEAR(dof_data[0].tau1_nuc,    1.0e99, 0.0,
             "tau1_nuc untouched (sentinel preserved)");
   TEST_NEAR(dof_data[0].sigma_n_nuc, 2.0e99, 0.0,
             "sigma_n_nuc untouched (sentinel preserved)");
}
```

**Test case:** The new test above demonstrates the hole. Under the
current production code (which correctly skips the two channels), the
sentinels stay pristine. Under a buggy refactor that writes to
tau1_nuc, the sentinel is overwritten and the assertion fires.

---

### [R3-002] [MODERATE] [POSSIBLE] [test_tpv104_nucleation.cpp:93] — `int(nuc_T / dt)` may truncate to 999 instead of 1000, placing `TestAccumulatorTelescopes` at the 1e-6 tolerance boundary

**Category:** EDGE_CASE (floating-point rounding + hardcoded tolerance)

**Description:**
```cpp
const real_t dt = 1e-3;
const int n_steps = static_cast<int>(TPV104Params::nuc_T / dt);
```
`TPV104Params::nuc_T = 1.0`. The expected n_steps is 1000. But
`1e-3` in IEEE double is not exactly representable — its true value
is 0.001000000000000000020816681711721685... When you compute
`1.0 / 0.001000000000000000020816...`, the result is
999.99999999999997..., which `static_cast<int>` truncates to **999**.

With `n_steps = 999`, the loop runs `k = 1..999`, final `t_end =
0.999`. The telescoping sum yields `smoothStep(0.999, 1.0) −
smoothStep(0, 1.0) = smoothStep(0.999, 1.0)`.

Evaluating `smoothStep(0.999, 1.0)`:
- `tau = -0.001`
- `tau² / (0.999 · (0.999 − 2)) = 1e-6 / -1.000999 ≈ -9.99e-7`
- `exp(-9.99e-7) ≈ 0.9999990`

Final `tau2_nuc = 45e6 · 0.9999990 = 44999955`, expected 45e6. Relative
error = 45 / 45e6 = **1e-6** — exactly at the test's tolerance. Rounding
one ULP either way can flip PASS↔FAIL.

On some toolchains `double(1.0) / double(1e-3)` rounds up to exactly
1000.0 (depending on the x87 vs SSE rounding mode); on others it
rounds down. The test is seed-stable but toolchain-fragile.

**Trigger:** Any compiler / libm combination that rounds `1.0/1e-3`
below 1000.0.

**Actual behavior:** Potentially flaky test at the 1e-6 boundary.

**Expected behavior:** Either (a) use an integer n_steps literal and
derive dt from it (`const int n = 1000; const real_t dt = nuc_T / n;`),
which guarantees k*dt == nuc_T bit-exactly when n divides evenly, OR
(b) extend the final iterate to `k = n_steps + 1` so the `t_end ≥
nuc_T` endpoint is always reached.

**Suggested fix (option a — explicit integer step count):**
```diff
-   const real_t dt = 1e-3;
-   const int n_steps = static_cast<int>(TPV104Params::nuc_T / dt);
+   // R3-002: derive dt from an integer step count to avoid
+   // truncation via static_cast<int>(1.0 / 1e-3) → 999 under some
+   // IEEE rounding paths.
+   constexpr int n_steps = 1000;
+   const real_t dt = TPV104Params::nuc_T / n_steps;

    for (int k = 1; k <= n_steps; ++k)
    {
       const real_t t_end = k * dt;
       ApplyNucleationIncremental_TPV104(dof_data, fault_coords, t_end, dt);
    }
```

**Test case:**
```cpp
void test_R3_002_integer_step_count_hits_T_nuc() {
   // With 1000 steps of nuc_T/1000, final t_end = 1000 * (nuc_T/1000)
   // MAY NOT equal nuc_T bit-exactly (accumulated rounding) but will
   // land within 1 ULP.  n_steps truncation cannot drop a step.
   constexpr int n_steps = 1000;
   const real_t dt = 1.0 / n_steps;
   const real_t t_final = n_steps * dt;
   // Both 1.0 and t_final must satisfy smoothStep(t_final, 1.0) ≥ 1 - 1e-15.
   const real_t ss = SmoothStep_TPV104(t_final, 1.0);
   ASSERT_GT(ss, 1.0 - 1e-15);
}
```

---

### [R3-003] MODERATE [test_tpv104_normal_sign.cpp:71-93] — `TestSigmaN0Propagation` tests field storage, not production init (`InitializeFaultDOFs_TPV104` not yet wired)

**Category:** DEVIATION (test title implies production coverage the test does not provide)

**Description:**
The plan (§4.10 Step 8) specifies T_TPV104_SIGN_2 as: *"TPV104 sign
propagation: at `t = 0`, `dof_data[i].sigma_n0 ==
TPV104Params::sigma_n > 0` at every QP."* The implemented test
manually populates DOFData:
```cpp
DOFData d;
d.sigma_n0 = TPV104Params::sigma_n;
d.tau1_0   = 0.0;
d.tau2_0   = TPV104Params::tau_ini;
```
It does NOT call any `InitializeFaultDOFs_TPV104` function (which
doesn't exist yet — Step 3 deferred to Step 9 per round-1 checklist).
So the test verifies only that `DOFData` can STORE a positive value,
not that any production init function writes a positive value. A
Step-9 bug that, say, computes `d.sigma_n0 = -TPV104Params::sigma_n`
would not be caught by this test.

The test's docstring acknowledges the deferral (*"Since Step 3 … has
not shipped yet, we test the structural property directly"*), but the
gate name `T_TPV104_SIGN_2` in the plan implies end-to-end coverage
that is not yet delivered.

**Trigger:** Step 9's `InitializeFaultDOFs_TPV104` introduces a sign
bug; existing tests miss it.

**Actual behavior:** Test passes with zero coverage of the production
init sign propagation.

**Expected behavior:** Either (a) rename/mark the current test
"T_TPV104_SIGN_2_storage" and defer the full sign-propagation check
to Step 9 (as `test_tpv104_smoke.cpp` or in `test_tpv104_setup.cpp`),
OR (b) ship a placeholder `InitializeFaultDOFs_TPV104` as a thin
Step-9 prerequisite and wire it into the test now.

**Suggested fix (option a):** adjust the acceptance matrix so plan
§4.10 Step 8 T_TPV104_SIGN_2 is explicitly split into a storage-gate
(closed here) and a production-init-gate (deferred to Step 9).
Update the test docstring:
```diff
 // ----------------------------------------------------------------------------
-// T_TPV104_SIGN_2 — TPV104 DOFData sigma_n0 propagation.
-// Since Step 3 (InitializeFaultDOFs_TPV104) has not shipped yet, we
-// test the structural property directly: a manually-populated DOFData
-// with `sigma_n0 = TPV104Params::sigma_n` reports the positive value
-// and the fault-flux computation treats it as compression.
+// T_TPV104_SIGN_2_storage — DOFData storage invariant only.
+//
+// The FULL T_TPV104_SIGN_2 sign-propagation gate (initialise DOFData
+// via production code and verify σ_n0 > 0 at every QP) requires
+// Step 3's InitializeFaultDOFs_TPV104, which has not shipped.
+// Move the end-to-end gate to test_tpv104_setup.cpp when Step 3
+// lands and mark T_TPV104_SIGN_2_init as the Step-9 gate.
 // ----------------------------------------------------------------------------
```

**Test case:**
```cpp
void test_R3_003_production_init_sets_positive_sigma_n0() {
   // Expected at Step 3 landing:
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   PopulateFaultQPCoords(fault_coords, /* test mesh */);
   InitializeFaultDOFs_TPV104(dof_data, fault_coords.size(), fault_coords);
   for (const auto &d : dof_data) {
      EXPECT_GT(d.sigma_n0, 0.0);
      EXPECT_NEAR(d.sigma_n0, TPV104Params::sigma_n, 1e-3);
   }
}
```
(Currently cannot compile — `InitializeFaultDOFs_TPV104` undefined.)

---

### [R3-004] LOW [tpv104_nucleation.hpp:108-110] — `dS == 0.0` exact-equality fast-path skips slightly-negative roundoff drift

**Category:** EDGE_CASE (potential tiny-drift accumulation)

**Description:**
```cpp
const real_t dS = SmoothStepIncrement_TPV104(t_substep_end, dt_substep,
                                             TPV104Params::nuc_T);
// Fast-path: before nucleation begins or after it completes, the
// increment is exactly 0 and we can skip the per-QP loop.
if (dS == 0.0) { return; }
```
`SmoothStep_TPV104` is monotonically non-decreasing by construction,
so mathematically `dS ≥ 0`. But floating-point roundoff can produce
a **slightly negative** dS at the boundaries (e.g., both `smoothStep`
calls return the same very-close-to-1 value via different code paths;
their subtraction yields, say, -2e-17). The exact-equality guard
misses this, and the per-QP loop runs with `dS < 0`, subtracting a
tiny negative from `tau2_nuc`. Over many calls, this could accumulate
a small reverse drift.

**Trigger:** Pathological `(t, dt, t0)` where smoothStep computed via
the ramp formula yields a value bit-different from 1 by one ULP at
both endpoints.

**Actual behavior:** Slightly negative dS survives the equality check
and feeds into the accumulator.

**Expected behavior:** Skip when dS is non-positive, or when
|dS| < ε_abs where ε_abs reflects roundoff.

**Suggested fix:**
```diff
-   // Fast-path: before nucleation begins or after it completes, the
-   // increment is exactly 0 and we can skip the per-QP loop.
-   if (dS == 0.0) { return; }
+   // Fast-path: clamp tiny non-positive increments to 0 (guards
+   // against sub-ULP negative drift from smoothStep roundoff at the
+   // ramp boundaries).  R3-004 review fix.
+   if (dS <= 0.0) { return; }
```

**Test case:** (LOW severity — drift magnitude is sub-ULP, no
repro in the current test envelope).

---

### [R3-005] LOW [test_tpv104_mesh.cpp:57-68, 112-123] — Mesh byte-identical tests hardcode header-line cutoffs (30 for 200 m, 36 for 1000 m)

**Category:** QUALITY (test fragility)

**Description:**
```cpp
auto tail_from_30 = [](const std::string &path) -> std::string { ... };
auto tail_from_36 = [](const std::string &path) -> std::string { ... };
```
The two line-number constants are hand-picked to match the CURRENT
header structure (200 m: 28-line comment block + 1 blank; 1000 m:
34-line comment block + 1 blank). If a maintainer adds a comment line
to the TPV104 200 m header but not the TPV102 200 m header (or vice
versa), the tail extraction silently slides past or stops short of
the code section, and the comparison would either miss a real
divergence or falsely flag one.

**Trigger:** Editing the header of either `.geo` file without
updating the test's hardcoded `lineno >= N` constants.

**Actual behavior:** Test passes or fails based on comment-block
length, not geometry content.

**Expected behavior:** Auto-locate the first non-comment, non-blank
line in each file and align from there, OR use a sentinel string like
`"lc = "` as the alignment anchor.

**Suggested fix:**
```diff
-   auto tail_from_30 = [](const std::string &path) -> std::string
+   // R3-005 review fix: locate the first code line ("lc = ...")
+   // and extract from there.  Immune to header-length changes.
+   auto tail_from_code = [](const std::string &path) -> std::string
    {
       std::ifstream in(path);
       std::string line;
       std::ostringstream oss;
-      int lineno = 0;
+      bool in_code = false;
       while (std::getline(in, line))
       {
-         ++lineno;
-         if (lineno >= 30) { oss << line << "\n"; }
+         // Strip trailing \r (handle CRLF line endings).
+         if (!line.empty() && line.back() == '\r') { line.pop_back(); }
+         if (!in_code && line.rfind("lc = ", 0) == 0) { in_code = true; }
+         if (in_code) { oss << line << "\n"; }
       }
       return oss.str();
    };
```
Apply the same pattern to the 1000 m test and remove `tail_from_36`.

**Test case:** manually add a blank comment line to the top of
`tpv104_200m.geo` only, and verify the test still reports byte-identity
against TPV102 (it should, because both tails now align on `lc = `).

---

### [R3-006] LOW [test_tpv104_mesh.cpp:144] — `"lc_fault = 500;"` literal string match fails if whitespace drifts

**Category:** QUALITY (brittle test)

**Description:**
```cpp
TEST_ASSERT(body.find("lc_fault = 500;") != std::string::npos,
            "500m variant sets lc_fault = 500");
```
If the `.geo` file is ever reformatted (e.g., `lc_fault=500;` without
spaces, or `lc_fault  = 500;` with tab-aligned spacing), the string
search returns `npos` and the test fails. Gmsh parses all these
whitespace variants identically.

**Trigger:** Any whitespace-preserving reformat of the `.geo` file.

**Actual behavior:** False negative.

**Expected behavior:** Use a regex or tokenised match.

**Suggested fix (minimal):** match the variable name and value
separately:
```diff
-   TEST_ASSERT(body.find("lc_fault = 500;") != std::string::npos,
-               "500m variant sets lc_fault = 500");
+   // R3-006: tolerate whitespace variation around `=` and `;`.
+   const size_t name = body.find("lc_fault");
+   TEST_ASSERT(name != std::string::npos,
+               "500m variant declares lc_fault");
+   const size_t eq = body.find('=', name);
+   const size_t semi = body.find(';', eq);
+   TEST_ASSERT(eq != std::string::npos && semi != std::string::npos,
+               "lc_fault = <value>; found");
+   std::string val = body.substr(eq + 1, semi - eq - 1);
+   // strip whitespace
+   val.erase(std::remove_if(val.begin(), val.end(), ::isspace),
+             val.end());
+   TEST_ASSERT(val == "500",
+               "500m variant sets lc_fault to 500");
```

---

### [R3-007] LOW [tpv104_nucleation.hpp:55-60] — `SmoothStepIncrement_TPV104` does not validate inputs; NaN / inf dt silently produces a nonsense increment

**Category:** ASSUMPTION (silent NaN propagation)

**Description:**
```cpp
inline real_t SmoothStepIncrement_TPV104(real_t current_time, real_t dt,
                                         real_t t0)
{
   return SmoothStep_TPV104(current_time, t0)
          - SmoothStep_TPV104(current_time - dt, t0);
}
```
If `dt` is NaN: `current_time - NaN = NaN`. `SmoothStep(NaN, t0)`:
`NaN <= 0` → false, `NaN < t0` → false, so it enters the
return-1.0 branch. Increment = `SmoothStep(current_time, t0) − 1.0`,
which for `current_time < t0` is a finite negative number — silently
fed into the accumulator as a spurious decrement.

If `dt` is +inf: `current_time - inf = -inf`. `SmoothStep(-inf, t0)`:
`-inf <= 0` is true → returns 0. Increment = `SmoothStep(current_time,
t0) − 0 = full ramp value`, inflating a single sub-step to the full
perturbation.

Production usage never passes NaN / inf (solver provides a finite
substep dt), but the function has no documented precondition.

**Trigger:** Any caller passing non-finite dt (a buggy ADER
predictor, a misaligned time-index calculation).

**Actual behavior:** Silent numerical garbage.

**Expected behavior:** Document the finite-dt precondition, OR assert
it via `MFEM_ASSERT(std::isfinite(dt) && dt > 0.0, ...)`.

**Suggested fix:**
```diff
 inline real_t SmoothStepIncrement_TPV104(real_t current_time, real_t dt,
                                          real_t t0)
 {
+   // R3-007: finite-dt precondition.  ApplyNucleationIncremental_TPV104
+   // assumes a strictly positive sub-step size; NaN or non-positive dt
+   // would silently corrupt the accumulator.
+   MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
+               "SmoothStepIncrement_TPV104: dt must be finite and "
+               "positive; got dt = " << dt);
    return SmoothStep_TPV104(current_time, t0)
           - SmoothStep_TPV104(current_time - dt, t0);
 }
```

---

## Summary

- Critical issues: 0
- Moderate issues: 3 (R3-001, R3-002, R3-003)
- Low issues: 4 (R3-004, R3-005, R3-006, R3-007)
- Plan compliance: PARTIAL — Steps 6, 8, 10 all shipped. Step 8's
  T_TPV104_SIGN_2 is weaker than the plan specifies (R3-003); Step 10's
  T_TPV104_MESH_1/2/3 are explicitly deferred to Step 9 (expected).
- Verdict: PASS WITH FIXES — R3-001 should be closed before Phase 3
  probe runs so the strike-slip invariant is actively guarded. R3-002
  should be closed now because it is one-line and the test is flaky at
  the tolerance boundary. R3-003 should be closed by a nomenclature
  update so the acceptance matrix distinguishes the storage gate from
  the production-init gate. R3-004..R3-007 are polish items that do
  not block Step 9 landing.

## Unreviewed Areas

- The physical `.msh` files generated from the `.geo` templates are
  not shipped and their count invariants (T_TPV104_MESH_1) remain
  DEFERRED to the Step 10 Phase-2 acceptance matrix. Gate only
  activates under the `pythonenv` conda environment with a running
  gmsh.
- `tpv104_nucleation.hpp`'s call site inside a future Step-7 sub-step
  iterator: we cannot verify correct integration cadence without the
  iterator code. The current tests exercise the accumulator in
  isolation.
- `BP5Params` default member values (T_TPV104_SIGN_1 relies on
  `BP5Params bp5;` default-constructing `sigma_n = 25e6`). This was not
  independently verified in this round; if a maintainer changes the
  BP5 config defaults, the test will catch the drift at runtime but
  not at compile time.

## Validation checklist — Round 3

Carry forward into Step 9 / Phase 3 planning.

### Required closures before Phase 3 probe runs

- [x] **R3-001 landed**: `tau1_nuc` and `sigma_n_nuc` explicitly
  asserted at 0 after every `ApplyNucleationIncremental_TPV104` call in
  the test suite. New dedicated `TestNoDipOrNormalStressNucleation`
  test with sentinel values that survive the full ramp.
  **CLOSED 2026-04-24** —
  `test_tpv104_nucleation.cpp:271-321` adds
  `TestNoDipOrNormalStressNucleation` with `tau1_nuc = 1e99`,
  `sigma_n_nuc = 2e99` sentinels that must be BIT-EXACTLY preserved
  after 200 sub-steps. The dip/normal invariants are also asserted
  inline after the accumulator calls in `TestAccumulatorTelescopes`
  (`:121-123`) and `TestGuardBounds` (`:264-266`).
- [x] **R3-002 landed**: `TestAccumulatorTelescopes` uses
  `constexpr int n_steps = 1000; const real_t dt = nuc_T / n_steps;`
  to avoid the `int(1.0/1e-3) = 999` truncation path. Re-run under
  clang, gcc, and icc to confirm no ULP-boundary flakiness.
  **CLOSED 2026-04-24** —
  `test_tpv104_nucleation.cpp:92-101` now derives `dt` from an
  integer step count: `constexpr int n_steps = 1000; const real_t
  dt = TPV104Params::nuc_T / n_steps;`. Multi-toolchain re-run
  deferred to CI but no longer at the 1e-6 boundary.

### Recommended closures before Step 9 landing

- [x] **R3-003 nomenclature update**: split T_TPV104_SIGN_2 into
  "_storage" (closed now) and "_init" (deferred to Step 9). Update the
  Phase-2 acceptance matrix accordingly.
  **CLOSED 2026-04-24** —
  `test_tpv104_normal_sign.cpp:65-87` renames the test to
  `T_TPV104_SIGN_2_storage`, docstring explicitly splits the gate
  into `_storage` (closed) and `_init` (deferred to Step 9).
  `:186-187` adds a deferred placeholder `T_TPV104_SIGN_2_init`
  stub that emits no TEST_ASSERT (matches R-012 convention).
- [x] **R3-004**: clamp `dS <= 0.0` in `ApplyNucleationIncremental_TPV104`
  (guards sub-ULP negative drift at ramp boundaries).
  **CLOSED 2026-04-24** —
  `tpv104_nucleation.hpp:117-124` replaces `if (dS == 0.0) { return; }`
  with `if (dS <= 0.0) { return; }`, fast-pathing both zero and
  sub-ULP negative increments.
- [x] **R3-005**: replace `tail_from_30` / `tail_from_36` line-number
  cutoffs in `test_tpv104_mesh.cpp` with a sentinel-string alignment
  anchor (`"lc = "`). Also strip trailing `\r` on read to handle CRLF.
  **CLOSED 2026-04-24** —
  `test_tpv104_mesh.cpp:42-60` introduces sentinel-anchor tail
  extraction: the reader locates the first line beginning `"lc = "`
  and returns everything from that line onward. Trailing `\r` is
  stripped via `if (!line.empty() && line.back() == '\r') {
  line.pop_back(); }`. Both 200m and 1000m byte-identical tests
  (`:80, :124`) now use the sentinel form; the hardcoded 30 / 36
  literals are removed.
- [x] **R3-006**: replace literal `"lc_fault = 500;"` match with a
  whitespace-tolerant name+value parse.
  **CLOSED 2026-04-24** —
  `test_tpv104_mesh.cpp:146-161` parses `lc_fault`'s value by
  locating `"lc_fault"`, then `=`, then `;`, stripping whitespace
  from the captured value via
  `std::remove_if(val.begin(), val.end(),
                  [](unsigned char c){ return std::isspace(c); })`.
- [x] **R3-007**: add `MFEM_ASSERT(std::isfinite(dt) && dt > 0.0, ...)`
  at the top of `SmoothStepIncrement_TPV104`.
  **CLOSED 2026-04-24** —
  `tpv104_nucleation.hpp:56-65` adds:
  ```cpp
  /// R3-007 (review round 3): `dt` must be finite and strictly positive.
  MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
              "SmoothStepIncrement_TPV104: dt must be finite and "
              "positive; got dt = " << dt);
  ```

## Round-3 closure verdict

All seven R3 findings are closed in code AND guarded by dedicated
tests or documentation anchors. The "TPV104 pure strike-slip, no
normal-stress nucleation" invariant — flagged in the round-3 summary
as the single most-pressing Phase-3 gate — now has three
independent guards: the sentinel test (`TestNoDipOrNormalStressNucleation`),
the inline `TEST_NEAR(tau1_nuc, 0.0, 1e-30)` after `TestAccumulatorTelescopes`,
and the same inline check after `TestGuardBounds`.

### Rolled forward from prior rounds

- [ ] Step 9 driver wiring uses `Rate_SRW` / per-QP V_w[i] and
  `SetProductionMode()` (round-1 deferral).
- [ ] Step 9 smoke test exercises `InitializeFaultDOFs_TPV104`
  end-to-end (closes R3-003's production-init gate).
- [ ] Step 10 Phase-2 acceptance matrix tracks T_TPV104_MESH_1/2/3
  as running under the `pythonenv` conda env; add a Makefile target
  `make tpv104-mesh` for .msh production.
- [ ] Step 13 probe-diff tool normalises the SeisSol → MFEM
  normal-stress sign flip (closes T_TPV104_SIGN_4).
- [ ] Phase 3 Probe 2 / Probe 4 threshold recalibration against the
  round-1 R-004, R-005, R-006 closures.

### Cross-round invariant anchors

The following invariants are now guarded at multiple levels and should
survive any future refactor:

- σ_n > 0 = compression (BP5, TPV102, TPV104 all agree) —
  `T_TPV104_SIGN_1`.
- `SlipLawSRWPsi` production mode aborts on base-virtual call —
  R-001 + T_SRW_7.
- Raw V passes through state-evolution formulas (no V_safe clamp) —
  R-002 + TestNoVsafeClamp.
- `(V/V_w)^8` via unrolled integer power, not `std::pow` —
  R-004 + R2-002 + TestIntegerPowerUnrolled.
- Newton input validation throws on invalid inputs —
  R2-004 + TestNewtonRejectsInvalidVPrev.
- Nucleation accumulator telescopes to Δτ₀·F(r) over [0, T_nuc] —
  T_TPV104_NUC_2 + T_TPV104_NUC_5.
- TPV104 pure strike-slip: tau1_nuc = sigma_n_nuc = 0 after
  nucleation calls — **NEEDS R3-001 CLOSURE**.
- TPV102 overwrite-pattern nucleation remains bit-unaffected by
  Step 6 —  T_TPV104_NUC_4.
- TPV104 mesh geometry byte-identical to TPV102 equivalents (modulo
  header) — T_TPV104_MESH_0/0b.
