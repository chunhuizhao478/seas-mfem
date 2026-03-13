# Phase 4 Quick Fixes

**Version 3.0** - Updated 2025-02-06 after all v2.0 documentation fixes applied

This document provides actionable fixes for discrepancies identified in the Phase 4 debug report.

---

## Status Summary

### Completed Fixes (v1.0)

| Fix | Priority | Status |
|-----|----------|--------|
| Fix 1: Remove double time increment in driver pseudocode | **HIGH** | COMPLETED |
| Fix 2: psi to theta notation in Phase 4 docs | MEDIUM | COMPLETED |
| Fix 3: ComputeRHS signature in docs | MEDIUM | COMPLETED |
| Fix 4: SetInitialCondition flow in docs | MEDIUM | COMPLETED |
| Fix 5: File organization section | LOW | COMPLETED |
| Fix 6: DomainOpType alias | LOW | COMPLETED |
| Fix 7: Traction formula comment | LOW | COMPLETED |
| Fix 8: Brent to Newton in coupling diagram | LOW | COMPLETED |

### Completed Fixes (v2.0)

| Fix | Priority | Status |
|-----|----------|--------|
| Fix 9: Main driver constructor signatures | **MEDIUM** | COMPLETED |
| Fix 10: Remove GetTime() and UpdateDomainAndTraction() from doc | LOW | COMPLETED |
| Fix 11: Add GetDomain() and GetFault() to doc | LOW | COMPLETED |
| Fix 12: Init() return value in doc | LOW | COMPLETED |
| Fix 13: Add state size check to doc | LOW | COMPLETED |
| Fix 14: params_ access pattern in doc | LOW | COMPLETED |
| Fix 15: AdaptiveTimeStepper missing API in doc | LOW | COMPLETED |
| Fix 16: Test file paths (integration to unit) | LOW | COMPLETED |
| Fix 17: Test framework and tolerances in doc | LOW | COMPLETED |

### Remaining Fixes (v3.0)

| Fix | Priority | Status |
|-----|----------|--------|
| Fix 18: Constructor body in doc overview section | **MEDIUM** | Pending |
| Fix 19: BP2 short test macro and tolerance | **MEDIUM** | Pending |
| Fix 20: Document all 9 test functions or add note | **MEDIUM** | Pending |
| Fix 21: GetParams() ordering in SetInitialCondition | LOW | Pending |
| Fix 22: MFEM_VERIFY error message detail | LOW | Pending |
| Fix 23: ComputeRHS return type note | LOW | Pending |
| Fix 24: Domain constructor arg count comment | LOW | Pending |
| Fix 25: VerifyStressEquilibrium tol parameter | LOW | Pending |
| Fix 26: Add destructor to doc class | LOW | Pending |
| Fix 27: t_final literal to named constant | LOW | Pending |
| Fix 28: BP2 test missing upper-bound slip check | LOW | Pending |

---

## Fix 18: Constructor Body in Doc Overview Section (MEDIUM)

**File:** `document/phase4_seas_operator.md`

**Current (lines 131-134):**
```cpp
SEASQuasiDynamicOperator(DomainOpType *domain,
                          RateStateFaultOperator<MeshType> *fault)
    : TimeDependentOperator(fault->StateSize()),
      domain_(domain), fault_(fault) {}
```

**Fixed (Option A -- expand body):**
```cpp
SEASQuasiDynamicOperator(DomainOpType *domain,
                          RateStateFaultOperator<MeshType> *fault)
    : TimeDependentOperator(fault->StateSize()),
      domain_(domain), fault_(fault)
{
    MFEM_VERIFY(domain_ != nullptr, "Domain operator must not be null");
    MFEM_VERIFY(fault_ != nullptr, "Fault operator must not be null");
    u_gf_ = std::make_unique<GridFuncType>(&domain_->GetFESpace());
    *u_gf_ = 0.0;
    slip_.SetSize(fault_->NumNodes());
    traction_.SetSize(fault_->NumNodes());
}
```

**Fixed (Option B -- add note):**
```cpp
SEASQuasiDynamicOperator(DomainOpType *domain,
                          RateStateFaultOperator<MeshType> *fault)
    : TimeDependentOperator(fault->StateSize()),
      domain_(domain), fault_(fault)
{
    // Null-pointer checks, u_gf_ allocation, work vector sizing
    // See Implementation Details section for full constructor body
}
```

**Why:** The implementation constructor allocates `u_gf_`, sizes work vectors (`slip_`, `traction_`), and adds null-pointer guards. The empty `{}` body in the overview would result in a non-functional operator if used as-is.

---

## Fix 19: BP2 Short Test Macro and Tolerance (MEDIUM)

**File:** `document/phase4_seas_operator.md`

**Current (line 402):**
```cpp
TEST_REL_NEAR(V_initial, params.V_init, 0.05, "Initial V ~ V_init");
```

**Fixed:**
```cpp
TEST_ASSERT(std::abs(V_initial - params.V_init) / params.V_init < 0.1,
            "Initial V_max ~ V_init");
```

**Why:** `test_bp2_short.cpp` does not define `TEST_REL_NEAR` (only `TEST_ASSERT` and `TEST_NEAR`). The tolerance is 10% (0.1), not 5% (0.05).

---

## Fix 20: Document All 9 Test Functions or Add Note (MEDIUM)

**File:** `document/phase4_seas_operator.md`

### Option A: Add a test summary note

**Add after line 385 (before Short BP2 Simulation section):**
```
// Full test suite includes 9 functions with 29 assertions total:
// TestOperatorConstruction, TestInitialCondition, TestSteadyStateSlip,
// TestStressBalanceMaintained, TestSlipConservation, TestThetaPositive,
// TestAdaptiveTimeStepper, TestRK4WithAdaptiveDt, TestMultDeterministic
//
// Only key physics tests are shown below. See source for full suite.
```

### Also fix the function name (line 375):

**Current:**
```cpp
void TestStressBalanceDuringStepping() {
```

**Fixed:**
```cpp
void TestStressBalanceMaintained() {
```

**Why:** The doc only shows 3 of 9 test functions and one has the wrong name. This understates the actual test coverage.

---

## Fix 21: GetParams() Ordering in SetInitialCondition (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (lines 264-274):**
```cpp
    // Verify stress equilibrium
    real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
    MFEM_VERIFY(eq_error < 1e-6,
                "Initial stress equilibrium error too large: " << eq_error);

    // Verify initial slip rate is close to V_init (10% tolerance for coarse meshes)
    const BP2Params &params = fault_->GetParams();
    real_t V_rel_err = ...
```

**Fixed:**
```cpp
    const BP2Params &params = fault_->GetParams();

    // Verify stress equilibrium
    real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
    MFEM_VERIFY(eq_error < 1e-6,
                "Initial stress equilibrium error too large: " << eq_error);

    // Verify initial slip rate is close to V_init (10% tolerance for coarse meshes)
    real_t V_rel_err = ...
```

**Why:** The implementation retrieves `GetParams()` before the equilibrium check (seas_operator.hpp line 162). Functionally equivalent but matches the implementation ordering.

---

## Fix 22: MFEM_VERIFY Error Message Detail (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (lines 273-274):**
```cpp
    MFEM_VERIFY(V_rel_err < 0.1,
                "Initial V_max differs from V_init by " << V_rel_err * 100 << "%");
```

**Fixed:**
```cpp
    MFEM_VERIFY(V_rel_err < 0.1,
                "Initial V_max = " << V_max
                << " differs from V_init = " << params.V_init
                << " by " << V_rel_err * 100 << "%");
```

**Why:** The implementation includes actual values for easier debugging.

---

## Fix 23: ComputeRHS Return Type Note (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (line 230):**
```cpp
    fault_->ComputeRHS(traction_, state, rate);
```

**Fixed:**
```cpp
    // ComputeRHS returns V_max (discarded here; use GetMaxSlipRate() instead)
    fault_->ComputeRHS(traction_, state, rate);
```

**Why:** `ComputeRHS` returns `real_t` (maximum slip rate). The return value is intentionally discarded in `Mult()` since `GetMaxSlipRate()` provides the same information, but this is not obvious without the comment.

---

## Fix 24: Domain Constructor Arg Count Comment (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (line 445):**
```cpp
    // Create domain operator (5 args: mesh, order, mu, Vp, Wf)
```

**Fixed:**
```cpp
    // Create domain operator (mesh, order, mu, Vp, Wf; 6th arg DGMethod defaults to IP)
```

**Why:** `AntiplaneDomainOperator` has 6 parameters: 4 required (`mesh`, `order`, `mu`, `Vp`) and 2 optional (`Wf=40.0e3`, `method=DGMethod::IP`).

---

## Fix 25: VerifyStressEquilibrium tol Parameter (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (line 265):**
```cpp
    real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
```

**Fixed:**
```cpp
    // VerifyStressEquilibrium has optional 3rd arg: tol (default 1e-10)
    real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
```

**Why:** The method signature is `VerifyStressEquilibrium(traction, state, tol=1e-10)`. The 2-arg call works using the default, but documenting the optional parameter helps users who may want a different tolerance.

---

## Fix 26: Add Destructor to Doc Class (LOW)

**File:** `document/phase4_seas_operator.md`

**Add after line 175 (after constructor declaration):**
```cpp
    /// Destructor
    ~SEASQuasiDynamicOperator() override = default;
```

**Why:** The implementation has an explicit defaulted destructor (seas_operator.hpp line 61). Important for classes with `unique_ptr` members.

---

## Fix 27: t_final Literal to Named Constant (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (line 396):**
```cpp
    real_t t_final = 1.0 * 3.15576e7;  // 1 year in seconds
```

**Fixed:**
```cpp
    real_t t_final = BP2Params::seconds_per_year;  // 1 year in seconds
```

**Why:** The implementation uses the named constant `BP2Params::seconds_per_year` defined in `config/bp2_params.hpp`. Using the constant avoids magic numbers and ensures consistency if the value changes.

---

## Fix 28: BP2 Test Missing Upper-Bound Slip Check (LOW)

**File:** `document/phase4_seas_operator.md`

**Current (line 415):**
```cpp
    // 5. Mean slip is positive
    TEST_ASSERT(mean_slip > 0.0, "Mean slip positive after 1 year");
```

**Fixed:**
```cpp
    // 5. Mean slip is positive and bounded
    TEST_ASSERT(mean_slip > 0.0, "Mean slip is positive after 1 year");
    TEST_ASSERT(mean_slip < 10.0 * expected_slip,
                "Mean slip within 10x of V_init*t");
```

**Why:** The implementation has both a lower-bound and upper-bound check on mean slip. The doc only shows the lower bound.

---

## Verification After Fixes

After applying documentation fixes, rebuild and run tests to confirm no regression:

```bash
conda activate mfem-dev
cd miniapps/seas
make seas_test_quasi_dynamic seas_test_bp2_short
./seas_test_quasi_dynamic    # Expected: 29/29 pass
./seas_test_bp2_short        # Expected: 8/8 pass
```

No code changes are required -- all remaining fixes (18-28) are documentation-only.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-02-06 | Initial: 8 actionable fixes (1 HIGH, 3 MEDIUM, 4 LOW) |
| 2.0 | 2025-02-06 | v1.0 fixes completed; 9 new remaining fixes (1 MEDIUM, 8 LOW) |
| 3.0 | 2025-02-06 | v2.0 fixes completed; 11 new remaining fixes (3 MEDIUM, 8 LOW) |
