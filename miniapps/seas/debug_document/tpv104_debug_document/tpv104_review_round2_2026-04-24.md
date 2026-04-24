# Code Review Round 2: TPV104 Steps 1, 2, 4, 5 — 2026-04-24

Fresh adversarial audit (all three passes re-executed from scratch). Confirms
the prior R-001..R-012 fixes from `tpv104_review_2026-04-24.md` landed, then
hunts for **new** bugs introduced by those fixes or missed in the first pass.

## Review Scope

- Plan: `miniapps/seas/debug_document/tpv104_debug_document/tpv104_debug_plan_2026-04-24.md`
  §4.10 Steps 1, 2, 4, 5.
- Prior review (closed): `tpv104_review_2026-04-24.md`.
- Files reviewed (current state):
  - `miniapps/seas/config/tpv104_params.hpp` (250 lines)
  - `miniapps/seas/friction/slip_law_srw_psi.hpp` (290 lines)
  - `miniapps/seas/friction/friction_coeff_stable.hpp` (159 lines, renamed
    from `friction_coeff_seissol.hpp` — public API renamed to
    `friction_stable::` + `FrictionCoefficientStable*`)
  - `miniapps/seas/dynamic/tpv104_friction_solver.hpp` (128 lines, public
    entry renamed to `SolveSlipRateNewtonStable`)
  - `miniapps/seas/tests/unit/test_tpv104_params.cpp` (211 lines)
  - `miniapps/seas/tests/unit/test_slip_law_srw_psi.cpp` (490 lines)
  - `miniapps/seas/tests/unit/test_friction_coeff_stable.cpp` (253 lines,
    renamed + namespace updated)
  - `miniapps/seas/tests/unit/test_tpv104_friction_solver_newton.cpp`
    (373 lines)
- SeisSol cross-reference:
  - `src/DynamicRupture/FrictionLaws/RateAndStateCommon.h`
  - `src/DynamicRupture/FrictionLaws/CpuImpl/FastVelocityWeakeningLaw.h`
  - `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h`

## Round-1 findings status

- R-001 (production-mode guard on base virtuals) — CLOSED. `production_mode_`
  flag + `throw std::runtime_error` in each base virtual. T_SRW_7 verifies
  the abort and R-011 guard.
- R-002 (drop `V_safe` clamp) — CLOSED. Raw `V` used everywhere. Dedicated
  `TestNoVsafeClamp` at `V = 1e-48` passes bit-identically with the inline
  reference.
- R-003 (standalone reference lambdas) — CLOSED. Both `test_slip_law_srw_psi`
  and `test_friction_coeff_stable` use fully inlined `reference_*` helpers
  with no calls into production.
- R-004 (unrolled integer 8th-power) — CLOSED. `IntegerPow8(x) = (x²)²·(x²)²`.
- R-005 (drop `V ≤ 0` early-return in μ) — CLOSED. `TestNoVNonPositiveGuard`
  verifies μ is odd in V.
- R-006 (unclamped Newton first-guess) — CLOSED. `TestFirstGuessUnclamped`
  runs `V_prev = 1e-50` and asserts bit-identity with the reference.
- R-007 (physically-informed Newton initial guess) — **NEW BUG INTRODUCED**.
  See R2-001 below.
- R-008 (`W` docstring says "full VW depth") — CLOSED.
- R-009 (`ComputeInitialPsiTPV104` uses logsinh) — CLOSED. Small-`a` stability
  test added.
- R-010 (narrower T_SRW_1 envelope) — CLOSED.
- R-011 (T_SRW_7 checks R-001 trap) — CLOSED.
- R-012 (T_TPV104_FS_4 stub removal) — CLOSED. No unconditional
  `TEST_ASSERT(true)`.

## Findings

### [R2-001] CRITICAL [test_tpv104_friction_solver_newton.cpp:138-147] — `physical_v_guess` massively overestimates V at the high-τ end of the test envelope, poisoning Newton convergence

**Category:** BUG (the R-007 fix replaces a mediocre guess with a worse one)

**Description:**
The R-007 "fix" replaced the old initial-guess formula
`V_guess = tau / (σ_n·f_0 + η_s)` with a derivation that assumes the large
`asinh`-argument regime:
```cpp
const real_t c = tau / (sigma_n * a);
const real_t arg = c - psi / a;
return V0 * std::exp(arg);
```
This derivation IGNORES the radiation-damping term `η·V` in the balance
τ = σ_n·μ(V) + η·V. For mid-range τ, `η·V` dominates and `V ≈ τ/η`, so the
asymptotic `V₀·exp((τ − σ_n·ψ)/(σ_n·a))` overshoots by many orders of
magnitude.

Concrete numbers from the test envelope (`tau ∈ [1e4, 1e8]`, `psi ∈ [0.3,
1.0]`, `a ∈ [0.008, 0.02]`, σ_n = 1.2e8, η_s = 4.625e6):

| sample | τ    | ψ   | a     | `physical_v_guess`  | true rest-state V | overshoot factor |
|---|---|---|---|---|---|---|
| worst | 1e8 | 0.3 | 0.008 | `1e-6·exp(66.7)` ≈ **9e22** m/s | ≈ 21.6 m/s (η-dominated) | **10²¹** |
| typical hot | 5e7 | 0.5 | 0.01 | `1e-6·exp(−8.3)` ≈ 2.5e-10 m/s | ≈ 6 m/s | — |
| typical cold | 1e5 | 0.7 | 0.015 | `1e-6·exp(−46.2)` ≈ 1e-26 m/s | ≈ 1e-26 m/s | OK |

Starting Newton at V = 9e22 m/s when the true root is 21.6 m/s: the first
Newton step computes g = -η·V at V = 9e22 ≈ -4e29; dg ≈ -1 (in the large-V
regime dμ/dV → 0). step = g/dg ≈ 9e22. V_new = max(almostZero, 0) = 1e-45,
a 67-decade collapse. From V = 1e-45 the Newton step-size is ~g/dg = 21/large
≈ tiny (due to huge dμ/dV at sub-normal V), so Newton stalls and cannot
converge in 60 iterations.

**Empirical consequence:** the test's `non_conv < N / 100 = 50` tolerance is
fragile — a single RNG seed that draws 51+ samples with `τ > 5e7` will fail.
Current seed `20260424u` happens to produce < 50 non-convergences, so the
test is green, but the guard is seed-dependent.

**Trigger:** Any future change to the RNG seed, N (5000 → 10000), or the
τ range, OR a compiler that rounds `exp(66)` slightly differently.

**Actual behavior:** Non-convergence rate is a function of the RNG seed,
not of the solver correctness.

**Expected behavior:** Initial guess that is bounded above by the
physically-plausible steady-state V, regardless of τ.

**Suggested fix:**
```diff
 static real_t physical_v_guess(real_t tau, real_t psi, real_t sigma_n,
-                               real_t a, real_t V0)
+                               real_t a, real_t V0, real_t eta_s)
 {
-   const real_t c = tau / (sigma_n * a);
-   // 2·V0·exp(-ψ/a)·sinh(c).  For large c, sinh(c) ≈ exp(c)/2;
-   // use the stable form exp(c - ψ/a) to avoid overflow when ψ/a and
-   // c are both large.
-   const real_t arg = c - psi / a;
-   return V0 * std::exp(arg);
+   // Bound the asymptotic asinh-inversion guess by the radiation-damping
+   // limit V ≤ τ/η_s, which is the overstressed-fault saturation.  For
+   // low τ where the asinh asymptote is valid, the exp term is small
+   // and unchanged; for high τ where η·V dominates, we cap at τ/η.
+   const real_t c   = tau / (sigma_n * a);
+   const real_t arg = c - psi / a;
+   const real_t asymptotic = V0 * std::exp(arg);
+   const real_t damping_cap = tau / eta_s;
+   return std::min(asymptotic, damping_cap);
 }
```
And update the call site:
```diff
-      const real_t V_guess = physical_v_guess(tau, psi, kSigmaN_TPV, a,
-                                              kV0_TPV);
+      const real_t V_guess = physical_v_guess(tau, psi, kSigmaN_TPV, a,
+                                              kV0_TPV, kEta_TPV);
```

**Test case:**
```cpp
void test_R2_001_physical_guess_bounded_at_high_tau() {
   // At τ=1e8, ψ=0.3, a=0.008: asymptotic guess overshoots by 20+ decades.
   const real_t g = physical_v_guess(1e8, 0.3, 1.2e8, 0.008, 1e-6, 4.625e6);
   // After the fix, the guess is capped by τ/η = 21.6 m/s.
   EXPECT_LE(g, 1e8 / 4.625e6 + 1e-12);
   EXPECT_GT(g, 0.0);

   // At low τ, the asymptotic form is untouched.
   const real_t g_low = physical_v_guess(1e4, 0.7, 1.2e8, 0.015, 1e-6, 4.625e6);
   const real_t expected_low = 1e-6 * std::exp(1e4/(1.2e8*0.015) - 0.7/0.015);
   EXPECT_NEAR(g_low, expected_low, 1e-30);
}
```

---

### [R2-002] MODERATE [test_slip_law_srw_psi.cpp:394-434] — `TestIntegerPowerUnrolled` does not actually assert that production uses the unrolled form instead of `std::pow`

**Category:** BUG (test is toothless against the regression it claims to guard)

**Description:**
The test prints a diff between `r8_unrolled = r4·r4` and `r8_pow = std::pow(r, 8.0)`
but never **asserts** that they differ (or that production matches the
unrolled form *specifically*). The only active `TEST_ASSERT` is:
```cpp
TEST_ASSERT(BitIdentical(got, expected),
            "production PsiSS_SRW uses unrolled integer power (R-004)");
```
where `expected` is computed inline using the same unrolled form that
production uses. So if a future refactor swaps production from
`IntegerPow8(x)` to `std::pow(x, 8.0)` AND the compiler produces
bit-identical output for `std::pow(r, 8.0)` and `r4·r4` (which it does on
some stdlibs via a specialisation), the test would pass with both forms —
the R-004 guard silently evaporates.

Even worse: on toolchains where `std::pow(r, 8.0) ≠ r4·r4` (most glibc
builds), the test would STILL pass because the `expected` side uses the
unrolled form and production matches it — the diff is relative to the
unrolled reference, not to `std::pow`.

**Trigger:** Any future refactor that accidentally reintroduces
`std::pow(x, 8.0)` in `PsiSS_SRW` or `UpdateStateAnalyticSlipLawSRW`.

**Actual behavior:** The test name claims an R-004 guard but the assertion
does not verify the production-vs-std::pow distinction.

**Expected behavior:** Explicitly assert `r8_unrolled != r8_pow` (skip the
test if they happen to match on this toolchain), AND assert production
matches the unrolled form (already done), AND assert production does NOT
bit-match an inline `std::pow`-based reference.

**Suggested fix:**
```diff
 void TestIntegerPowerUnrolled()
 {
    std::cout << "\n[R-004] (V/V_w)^8 matches unrolled integer power\n";
 
    const real_t V = 0.11, V_w = 0.1;
    const real_t r = V / V_w;
    const real_t r2 = r * r;
    const real_t r4 = r2 * r2;
    const real_t r8_unrolled = r4 * r4;
    const real_t r8_pow      = std::pow(r, 8.0);
 
-   std::cout << "  r^8 unrolled = " << r8_unrolled << "\n"
-             << "  r^8 std::pow = " << r8_pow << "\n"
-             << "  diff ULP      = "
-             << std::abs(r8_unrolled - r8_pow) / r8_unrolled << "\n";
+   // If the toolchain specialises std::pow(x, 8.0) to the unrolled form,
+   // the R-004 distinction is moot on THIS toolchain — but we still must
+   // confirm production is correct on toolchains where they DO differ.
+   const bool toolchain_distinguishes = (r8_unrolled != r8_pow);
+   std::cout << "  toolchain distinguishes std::pow from unrolled: "
+             << toolchain_distinguishes << "\n";
 
    // ... existing expected recomputation ...
 
    const real_t got = SlipLawSRWPsi::PsiSS_SRW(V_sample, V_w, a, b,
                                                V0, f0, muW);
    TEST_ASSERT(BitIdentical(got, expected),
                "production PsiSS_SRW uses unrolled integer power (R-004)");
+
+   // Compute an "anti-reference" that uses std::pow and assert that
+   // production does NOT match it on toolchains where the two forms
+   // differ.  This catches the regression where production is swapped
+   // back to std::pow(x, 8.0).
+   if (toolchain_distinguishes) {
+      const real_t denom_pow = std::pow(1.0 + r8_pow, 1.0 / 8.0);
+      const real_t f_ss_pow  = muW + (f_LV - muW) / denom_pow;
+      const real_t x_ref     = 2.0 * V0 / V_sample;
+      const real_t c_ref     = f_ss_pow / a;
+      const real_t signCp    = (c_ref >= 0) ? 1.0 : -1.0;
+      const real_t absCp     = std::abs(c_ref);
+      const real_t anti_expected =
+         a * (absCp + std::log(x_ref / 2.0 * -signCp * std::expm1(-2.0 * absCp)));
+      TEST_ASSERT(!BitIdentical(got, anti_expected),
+                  "production PsiSS_SRW does NOT bit-match std::pow-based "
+                  "reference (regression guard)");
+   }
 }
```

**Test case:** the above test body, manually exercised by replacing
`IntegerPow8` with `std::pow(x, 8.0)` in production — the new
`!BitIdentical(got, anti_expected)` assertion must then FAIL, proving the
guard is live.

---

### [R2-003] MODERATE [slip_law_srw_psi.hpp:208-214] — `RateDerivativeV` produces `NaN` at `V = 0`

**Category:** BUG (latent NaN in a test-fixture-mode API)

**Description:**
```cpp
real_t RateDerivativeV(real_t V, real_t psi, real_t L) const override
{
   ...
   const real_t eps = std::max(static_cast<real_t>(1e-8) * std::abs(V),
                               static_cast<real_t>(1e-18));
   const real_t ps_plus  = SteadyState_SRW(V + eps, V_w_default_, a_);
   const real_t ps_mid   = SteadyState_SRW(V,       V_w_default_, a_);
   const real_t dps_dV   = (ps_plus - ps_mid) / eps;
   return -(1.0 / L) * (psi - ps_mid) + (V / L) * dps_dV;
}
```
With the R-002 fix, `SteadyState_SRW(V = 0, ...)` no longer clamps V. So:
- `log(V/V0) = log(0) = −∞`
- `f_LV = max(0, f0 − (b − a)·(−∞)) = max(0, +∞) = +∞`
- `f_ss = muW + (+∞ − muW) / (1 + 0^8)^(1/8) = +∞`
- `psi_ss = a · logsinh(inf, inf/a)`. `logsinh(x, c)` with `c = inf`:
  `expm1(−2·inf) = −1`, `log(x/2 · 1) = log(x·0.5·1) = log(x/2)` where
  `x = 2·V0/V = 2·V0/0 = +∞`. So `log(+∞) = +∞`. `psi_ss = +∞`.

Then `ps_mid = +∞` and `ps_plus = finite` (since V+eps = 1e-18 > 0). So
`dps_dV = (finite − ∞) / eps = −∞`. Finally:
`-(V/L)·dps_dV = -(0/L)·(−∞) = 0·(−∞) = NaN`.

T_SRW_7 happens to call `RateDerivativeV(V = 0.2, ...)` where V > 0, so
the NaN case is never exercised.

**Trigger:** Any caller passing `V = 0` (or sub-normal V) into the base
virtual. The current test uses V=0.2 and never hits this.

**Actual behavior:** Returns NaN silently, propagates into any caller
using the virtual for implicit-integrator diagnostics.

**Expected behavior:** Either return a sensible value (the analytic
closed-form of `−(1/L)·(ψ − ψ_ss(0))` and `0·dψ_ss/dV = 0` for V = 0),
OR guard with a loud abort.

**Suggested fix:**
```diff
 real_t RateDerivativeV(real_t V, real_t psi, real_t L) const override
 {
    if (production_mode_)
    {
       throw std::runtime_error(
          "SlipLawSRWPsi::RateDerivativeV base-virtual called in "
          "production mode.");
    }
+   // At V = 0 the (V/L)·dψ_ss/dV term vanishes algebraically even
+   // though dψ_ss/dV → -∞; the first term -(1/L)·(ψ - ψ_ss(0)) is
+   // also singular because ψ_ss(0) = +∞.  The analytic limit as
+   // V → 0⁺ is -(V/L)·(dψ_ss/dV) → 0 (locked fault contributes no
+   // Rate derivative).  Return that limit directly to avoid NaN.
+   if (V == 0.0)
+   {
+      return -(1.0 / L) * (psi - SteadyState_SRW(
+         std::numeric_limits<real_t>::min(), V_w_default_, a_));
+   }
    const real_t eps = std::max(static_cast<real_t>(1e-8) * std::abs(V),
                                static_cast<real_t>(1e-18));
    ...
 }
```
Alternatively, document the `V ≠ 0` precondition loudly and add a
runtime assertion.

**Test case:**
```cpp
void test_R2_003_rate_derivative_at_zero_V() {
   SlipLawSRWPsi law(0.01, 0.014, 1e-6, 0.6, 0.1, /*V_w_default=*/0.1);
   const real_t d = law.RateDerivativeV(0.0, /*psi=*/0.5, /*L=*/0.4);
   EXPECT_TRUE(std::isfinite(d));
}
```

---

### [R2-004] MODERATE [tpv104_friction_solver.hpp:69-123] — `SolveSlipRateNewtonStable` does not validate `V_prev`; negative V_prev corrupts the first iterate

**Category:** EDGE_CASE (silent propagation of invalid input)

**Description:**
```cpp
// R-006: first guess unclamped — match reference.
real_t V = V_prev;
...
for (it = 0; it < max_iter; ++it)
{
   const real_t mu = friction_stable::FrictionCoefficientStable(
                        V, psi, a, V0);
   g = -inv_eta_s * (sigma_n_abs * mu - tau_abs) - V;
   ...
}
```
With the R-005 fix, `FrictionCoefficientStable` is an odd function of V:
`μ(−V) = −μ(V)`. If `V_prev < 0` (accidentally, e.g. from a buggy caller
passing un-initialised memory):
- `μ(V_prev) < 0`, so `σ_n · μ < 0`.
- `g = -(σ_n·μ − τ)/η − V_prev = (|σ_n·μ| + τ)/η − V_prev > 0` (since
  V_prev < 0, subtracting makes it larger).
- Newton step `V − g/dg` with `dg < 0`: step is negative, so V_new =
  max(almostZero, V_old - negative) = max(almostZero, V_old + positive).
  Might land at almostZero or oscillate.

Either way, the physically meaningless V_prev < 0 propagates at least one
iteration with negative μ — polluting the friction law and potentially
the convergence trajectory.

**Trigger:** Caller bug that passes an uninitialised or intermediate-step
negative V into the solver. Unlikely in the current Phase-2 tests
(T_TPV104_FS_3 draws V_prev from `std::pow(10.0, logV(rng))` with `logV ∈
[-18, -6]`, so always positive). Real risk is in Step 9 driver wiring
where the per-QP `slip_rate` field is read without validation.

**Actual behavior:** Negative V_prev silently propagates for ≥ 1
iteration.

**Expected behavior:** Either (a) assert V_prev ≥ 0 at entry (loud
failure), OR (b) clamp V_prev = max(V_prev, almostZero) with a comment
documenting the R-006 tension.

**Suggested fix (option a — loud failure):**
```diff
 inline real_t SolveSlipRateNewtonStable(
    real_t tau_abs, real_t psi, real_t sigma_n, real_t eta_s,
    real_t a, real_t V0,
    real_t V_prev,
    int max_iter = 60,
    real_t tol   = 1e-8,
    int *iterations = nullptr,
    bool *has_converged = nullptr)
 {
+   MFEM_ASSERT(std::isfinite(V_prev) && V_prev >= 0.0,
+               "SolveSlipRateNewtonStable: V_prev must be finite and "
+               "non-negative; got V_prev = " << V_prev);
+   MFEM_ASSERT(std::isfinite(tau_abs) && tau_abs >= 0.0,
+               "SolveSlipRateNewtonStable: tau_abs must be finite and "
+               "non-negative; got " << tau_abs);
+   MFEM_ASSERT(std::isfinite(sigma_n),
+               "SolveSlipRateNewtonStable: sigma_n must be finite; got "
+               << sigma_n);
+   MFEM_ASSERT(std::isfinite(eta_s) && eta_s > 0.0,
+               "SolveSlipRateNewtonStable: eta_s must be finite and "
+               "positive; got " << eta_s);
    if (tau_abs <= 0.0)
    ...
```

**Test case:**
```cpp
void test_R2_004_newton_rejects_negative_V_prev() {
   int it = 0; bool conv = false;
   bool aborted = false;
   try {
      SolveSlipRateNewtonStable(40e6, 0.5, 120e6, 4.625e6, 0.01, 1e-6,
                                /*V_prev=*/-1e-12, 60, 1e-8, &it, &conv);
   } catch (...) { aborted = true; }
   EXPECT_TRUE(aborted);  // MFEM_ASSERT fires
}
```

---

### [R2-005] LOW [test_tpv104_friction_solver_newton.cpp:233-283] — `TestStandaloneByteMatch` threshold relaxed from bit-identical to 1e-13 rel without test-level regression coverage

**Category:** QUALITY (unjustified weakening of the R-003 guard)

**Description:**
T_TPV104_FS_3 was originally supposed to verify bit-identity between the
production Newton loop and a standalone inline reference (R-003 in the
prior review). The current version accepts 1e-13 relative + matching
iter/conv:
```cpp
// The production and the standalone reference use algebraically
// identical arithmetic but the compiler emits them as distinct
// instruction sequences (inlined helpers vs straight-line).  Over
// the Newton iterate sequence, ULP-level round-off differences
// accumulate.  Accept ≤ 1e-13 relative and matching iter/conv.
// The point of R-003 is to catch FORMULA drift, not bit-identity.
```
This justification is reasonable (ULP drift from instruction-sequence
reordering is real), but the 1e-13 bound is pulled out of a hat —
there is no test that **measures** the actual drift and fails loudly if
it grows. A future refactor that introduces a subtle formula drift of
~1e-14 would be accepted, eroding the guard over time.

**Trigger:** Slow formula drift via rounding ordering changes across
refactors.

**Actual behavior:** The 1e-13 threshold is fixed; no trend tracking.

**Expected behavior:** Either (a) print `max_rel` across all 1000
samples and require it stays below a visible ratchet (e.g., the current
measured max), or (b) add a separate bit-identity test on a narrower
sub-envelope where instruction-ordering is not an issue (e.g., single
Newton iteration on a canonical rest state).

**Suggested fix:**
```diff
    std::cout << "  max rel err = " << max_rel << " over " << N << " samples\n";
    TEST_ASSERT(ok == 1000,
                "All 1000 Newton solves match standalone reference to 1e-13 "
                "rel with matching iter/conv");
+
+   // R2-005 ratchet: record the observed max-rel-err.  If it creeps
+   // above 1e-14 on this toolchain the test will still pass at 1e-13
+   // but fail at the ratchet — catching slow formula drift.
+   TEST_ASSERT(max_rel < 1e-14,
+               "max rel err stays below the 1e-14 ratchet (catches slow "
+               "formula drift across refactors)");
 }
```
Note: the ratchet value (1e-14) should match the actually-observed
max_rel on the reference toolchain; run the test once and bake in
0.5× that number so future refactors must stay at least as good.

**Test case:** the test itself serves as the test; no additional
fixture needed.

---

### [R2-006] LOW [tpv104_friction_solver.hpp:69-90] — τ or σ_n of `NaN` / `+inf` silently passes through the early-return checks

**Category:** EDGE_CASE (silent NaN propagation)

**Description:**
```cpp
if (tau_abs <= 0.0) { ... return 0.0; }
if (sigma_n <= 0.0) { ... return (eta_s > 0.0) ? (tau_abs / eta_s) : 0.0; }
```
`NaN <= 0.0` evaluates to `false`. So NaN inputs bypass the guard and
feed into the Newton loop, producing NaN output. The caller sees
`has_converged = true` (eventually — `|NaN| < tol` is false, so we
exit with `has_converged = |NaN| < tol = false`).

Actually — in the failure path, `has_converged = std::abs(g) < tol`
where `g = NaN`, so `abs(NaN) = NaN`, and `NaN < tol = false`. So
`has_converged = false` as expected, and the returned V is NaN. That
is the "correct" propagation, so not wrong per se — but a finite-input
caller can get a silent NaN back if any of `(tau, psi, sigma_n, eta_s,
a, V0, V_prev)` is NaN.

**Trigger:** Any NaN in the caller's friction inputs. Production path
should never have NaN, but a buggy ADER predictor could produce one.

**Actual behavior:** Silent NaN propagation + `has_converged = false`.

**Expected behavior:** Either (a) assert `std::isfinite` on all inputs
(matching R2-004's recommendation), OR (b) document the NaN passthrough
contract.

**Suggested fix:** See R2-004 — the `MFEM_ASSERT(std::isfinite(...))`
guards would cover this case too. If kept separate:
```diff
+   if (!std::isfinite(tau_abs) || !std::isfinite(sigma_n) ||
+       !std::isfinite(eta_s)   || !std::isfinite(a)       ||
+       !std::isfinite(V0)      || !std::isfinite(V_prev)  ||
+       !std::isfinite(psi))
+   {
+      if (iterations)    { *iterations = 0; }
+      if (has_converged) { *has_converged = false; }
+      return std::numeric_limits<real_t>::quiet_NaN();
+   }
    if (tau_abs <= 0.0) { ... return 0.0; }
```
**Do not fix** (LOW) if R2-004's assertion lands — it subsumes this.

---

### [R2-007] LOW [slip_law_srw_psi.hpp:156, 223] — `SlipLawSRWPsi` default `V_w_default = 1.0` is the strengthening value, not the VW-core value

**Category:** QUALITY (subtle default that biases test fixtures toward the wrong regime)

**Description:**
```cpp
SlipLawSRWPsi(real_t a_scalar, real_t b_scalar, real_t V0_scalar,
              real_t f0_scalar, real_t muW_scalar,
              real_t V_w_default = 1.0)
```
The default `V_w_default = 1.0` corresponds to `V_w_out` (outside the VW
core), i.e. the strengthening regime. A test fixture that constructs
`SlipLawSRWPsi(a, b, V0, f0, muW)` without specifying V_w_default
silently uses the strengthening friction. If the test author expects
"default = TPV104 core physics", they will silently get wrong behavior
(rupture cannot nucleate in strengthening).

None of the current tests are tripped by this (T_SRW_1 random V_w, T_SRW_3
passes V_w explicitly, T_SRW_7 uses `kVw_TPV = 0.1`), but it is a
landmine for Step 9 onward when the driver constructs a single shared
`SlipLawSRWPsi` and per-QP `V_w[i]` overrides.

**Trigger:** Any future test or driver construction
`SlipLawSRWPsi(scalars...)` without an explicit V_w_default.

**Actual behavior:** Silent fallthrough to V_w = 1.0.

**Expected behavior:** Either (a) remove the default argument so every
construction MUST pass V_w_default explicitly, OR (b) pick a default
that matches TPV104's VW core (0.1), documented loudly.

**Suggested fix (option a):**
```diff
    SlipLawSRWPsi(real_t a_scalar, real_t b_scalar, real_t V0_scalar,
                  real_t f0_scalar, real_t muW_scalar,
-                 real_t V_w_default = 1.0)
+                 real_t V_w_default)
       : a_(a_scalar), b_(b_scalar), V0_(V0_scalar),
```
Existing call sites will break at compile time, flagging which tests
need an explicit V_w_default. Update each.

---

## Summary

- Critical issues: 1 (R2-001)
- Moderate issues: 3 (R2-002, R2-003, R2-004)
- Low issues: 3 (R2-005, R2-006, R2-007)
- Plan compliance: PARTIAL — all prior Steps 1/2/4/5 fixes landed cleanly
  except the R-007 "physically-informed V_guess" which introduces a new
  large-τ overshoot bug. The other new findings are latent edge cases
  (NaN/negative inputs, test toothlessness) that do not block Step 9 but
  should be cleared before Phase 3 production runs.
- Verdict: PASS WITH FIXES — R2-001 must be fixed before T_TPV104_FS_1 is
  accepted as a stable Phase-2 gate (currently seed-dependent). R2-002
  must be fixed to make the R-004 guard actually guard. R2-003 and R2-004
  should be fixed before Step 9 wires `SlipLawSRWPsi` into the driver;
  otherwise a NaN or negative V from the ADER predictor can silently
  corrupt friction state.

## Validation checklist going into Step 9

- [x] R2-001: `physical_v_guess` capped by `τ/η_s`; T_TPV104_FS_1 passes
  with three different RNG seeds (`20260424`, `1`, `42`).
  **CLOSED 2026-04-24** —
  `test_tpv104_friction_solver_newton.cpp:144-152` caps with
  `std::min(asymptotic, tau / eta_s)`; dedicated
  `TestPhysicalGuessBounded` at `:370-395` asserts high-τ ≤ τ/η_s and
  low-τ = uncapped asymptotic.
- [x] R2-002: `TestIntegerPowerUnrolled` includes the `!BitIdentical(got,
  anti_expected)` regression guard.
  **CLOSED 2026-04-24** —
  `test_slip_law_srw_psi.cpp:460-484` builds an anti-reference with
  `std::pow(r, 8.0)` on toolchains that distinguish it from the
  unrolled form, and asserts production does **not** byte-match it.
- [x] R2-003: `RateDerivativeV(V=0, ...)` returns finite; new unit test.
  **CLOSED 2026-04-24** —
  `slip_law_srw_psi.hpp:227-232` adds an analytic-limit branch using
  `V_proxy = std::numeric_limits<real_t>::min()` for ψ_ss evaluation;
  `test_slip_law_srw_psi.cpp:357-379` (`TestRateDerivativeVAtZero`)
  asserts finite output at `V = 0` and `V = 1e-30`.
- [x] R2-004: `SolveSlipRateNewtonStable` aborts on `V_prev < 0` OR on
  any non-finite input; new unit test covers the assertion.
  **CLOSED 2026-04-24** —
  `tpv104_friction_solver.hpp:85-126` throws `std::runtime_error` on
  each of `tau_abs < 0 | !finite`, `!finite(sigma_n)`,
  `eta_s ≤ 0 | !finite`, `a ≤ 0 | !finite`, `V0 ≤ 0 | !finite`,
  `!finite(psi)`, and `V_prev < 0 | !finite`.
  `TestNewtonRejectsInvalidVPrev` at `:400-456` covers negative V_prev,
  NaN V_prev, NaN τ, and the V_prev = 0 accept case.
- [x] R2-005: ratchet test tracks max-rel-err across refactors.
  **CLOSED 2026-04-24** —
  `test_tpv104_friction_solver_newton.cpp:292-299` adds
  `TEST_ASSERT(max_rel < 3e-14, …)` as a drift monitor beside the
  1e-13 ok-threshold. Ratchet value calibrated from the observed
  ~1.16e-14 on clang-19 + O2.
- [x] R2-006: NaN propagation policy documented (subsumed by R2-004 if
  `std::isfinite` guards land).
  **CLOSED 2026-04-24** — subsumed by R2-004: every non-finite input
  now throws loudly at solver entry; NaN cannot reach the Newton loop.
- [x] R2-007: `V_w_default` default value removed from the class
  constructor OR changed to 0.1 (VW-core value) with a loud comment.
  **CLOSED 2026-04-24** — `slip_law_srw_psi.hpp:163-169` makes
  `V_w_default` a **required** constructor argument (default = 1.0
  removed). All test fixtures updated to pass `kVw_TPV = 0.1`
  explicitly.

## Round-2 closure verdict

All seven R2 findings are closed in code AND guarded by dedicated unit
tests. Remaining deferred items below are Step-9 or Phase-3 work and do
not affect Steps 1/2/4/5 acceptance.

## Deferred (out of Round-2 scope)

The following round-1 checklist items remain open pending downstream
work; not closeable at the current phase:

- Step 9 driver wiring (uses `Rate_SRW` with per-QP V_w[i] and a[i],
  calls `SetProductionMode()`, verifies V_w[i] side-channel lifetime,
  CLI default = Newton). Tracked by Step 9 smoke tests
  (`test_tpv104_smoke.cpp`, not yet shipped).
- Phase 3.B threshold recalibrations on Probes 2 and 4 — Phase 3 work.
- Phase 3 probe-output interpretation rules — Phase 3 work.

## Residual low-severity observations (new, 2026-04-24)

Noted during checklist closure; none block Step 9 landing. Filed for
future-attention during Phase 3 / Step 9.

1. **Reference-vs-production behavior at τ < 0.** Production
   `SolveSlipRateNewtonStable` throws on `tau_abs < 0` (R2-004);
   test-local `reference_newton` (`:95-128`) silently returns 0.0 via
   the `tau_abs <= 0.0` branch. `TestStandaloneByteMatch` never draws
   negative τ (range `[1e5, 1e8]`), so the divergence is invisible.
   Harmless at present — document if the reference ever serves a wider
   diagnostic role.
2. **`TestPhysicalGuessBounded` does not cover the mid-τ crossover**
   where `asymptotic == damping_cap`. The `std::min` branch is trivial
   but a one-line crossover check would harden the coverage.
3. **Error messages use `const char* + std::to_string(...)`** — depends
   on the `operator+(const char*, const std::string&)` overload from
   `<string>`. Standard C++; compiles on all supported toolchains.
   Not a bug, just noted.
4. **`TestPhysicalGuessBounded` and `TestRateDerivativeVAtZero` do not
   verify interaction with production mode.** `SetProductionMode()`
   correctly guards these entry points (R-001 check is ordered before
   R2-003's V=0 branch in `RateDerivativeV`); `T_SRW_7` already proves
   the production-mode throw for base virtuals. No additional test
   needed.

## Unreviewed areas

- `friction_stable::ArsinhExp` stable-branch math identity: verified on
  paper (`asinh(x·eᶜ) = sign(x)·(c + log(|x| + √(x² + e⁻²ᶜ)))`), but no
  unit test covers the case `c < 0` with large |x| that routes to the
  stable branch. Current `TestReferenceByteMatch` envelope (`V ∈ [1e-18,
  10]`, `ψ ∈ [0.1, 2.0]`, `a ∈ [0.005, 0.04]`) has `c = ψ/a ∈ [2.5, 400]`
  — all positive. A sample with negative c would exercise
  `exp_neg = 1 / cExp` on the stable branch.
- `ComputeA_TPV104` / `ComputeVw_TPV104` with inputs outside the test's
  `|x| ≤ 20 km, z ≤ 20 km` range: the `Boxcar_TPV104` function asserts
  `ax >= W + w → return 0.0` — no issue up to `double` overflow. Not a
  concern.
- Build/Makefile integration (Step 12) has not shipped; tests' ability to
  locate `friction/friction_coeff_stable.hpp` at multiple working
  directories (`in.open("friction/...")` → `../friction/...` →
  `../../friction/...`) is opportunistic and may fail under Step 12's
  final build-tree layout. Flag for Step 12 landing.
- `slip_law_srw_psi.hpp`'s `SlipLawSRWPsi::InProductionMode()` accessor
  is defined but never tested; no risk but noted for completeness.
