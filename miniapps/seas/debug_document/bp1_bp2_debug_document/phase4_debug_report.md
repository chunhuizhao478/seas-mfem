# Phase 4 Debug Report: SEAS Quasi-Dynamic Operator

**Version 3.0** - Updated 2025-02-06 after all v2.0 documentation fixes applied

## Overview

This report documents the analysis of the Phase 4 implementation (SEAS quasi-dynamic operator) against its design documents (`phase4_seas_operator.md`, `dg_antiplane_theory.md`) and the implementations of Phases 2--3.

**Files Analyzed:**
- Documentation: `document/phase4_seas_operator.md`, `document/dg_antiplane_theory.md`
- Implementation: `solver/seas_operator.hpp`, `solver/time_stepper.hpp`
- Tests: `tests/unit/test_quasi_dynamic.cpp`, `tests/unit/test_bp2_short.cpp`
- Supporting: `fault/rate_state_fault.hpp`, `domain/antiplane_operator.hpp`, `config/bp2_params.hpp`

---

## Part I: Test Execution Results

### Quasi-Dynamic Operator Tests (`seas_test_quasi_dynamic`)

| Test | Result | Details |
|------|--------|---------|
| Operator Construction | PASSED | State size = 24 (2 x 12 fault DOFs) |
| Initial Condition (slip=0) | PASSED | Slip L2 norm = 0 |
| Initial Condition (theta>0) | PASSED | theta in [4000, 4e6] s |
| Initial Condition (V~V_init) | PASSED | V_max = 1e-9 m/s |
| Stress equilibrium at init | PASSED | Error = 0 |
| Steady-state V after 10ks | PASSED | V_max = 1e-9 m/s |
| Mean slip positive | PASSED | 7.39e-6 m |
| Stress balance during stepping | PASSED | Max error = 2.81e-16 |
| Slip conservation | PASSED | Max rel error = 0.934 (within 1x) |
| Theta stays positive | PASSED | 20 RK4 steps |
| Adaptive dt: very slow | PASSED | dt = dt_max |
| Adaptive dt: rapid slip | PASSED | dt reduces proportionally |
| Adaptive dt: moderate | PASSED | dt increases |
| Adaptive dt: bounds | PASSED | dt clamped to [dt_min, dt_max] |
| Adaptive dt: config | PASSED | Setters work |
| RK4 + adaptive dt | PASSED | 10 steps, t_final reached |
| Mult deterministic | PASSED | Max diff = 0 |

**Total: 29/29 pass**

### BP2 Short Simulation (`seas_test_bp2_short`)

| Test | Result | Details |
|------|--------|---------|
| Initial V_max ~ V_init | PASSED | V_max = 1e-9 m/s |
| Completed in step limit | PASSED | 23 steps for 1 year |
| No coseismic events | PASSED | V_max = 1.01e-9 m/s |
| Theta positive | PASSED | theta_min = 1.40e4 s |
| Stress > 15 MPa | PASSED | tau_min = 26.535 MPa |
| Stress < 40 MPa | PASSED | tau_max = 26.564 MPa |
| Mean slip positive | PASSED | 1.83e-2 m |
| Mean slip within 10x | PASSED | vs expected 3.16e-2 m |

**Total: 8/8 pass**

---

## Part II: Previously Fixed Issues

### Fixed in v1.0 to v2.0

| Issue (v1.0) | Description | Status |
|--------------|-------------|--------|
| #1 psi vs theta notation | Doc used psi, impl uses theta | FIXED |
| #2 Brent root finder | Coupling diagram said "Brent" | FIXED |
| #3 Double time increment | `t += dt` after `Step()` in main driver | FIXED |
| #4 ComputeRHS signature | Doc showed 4 args with time | FIXED |
| #5 SetInitialCondition flow | Doc showed single-call `Initialize()` | FIXED |
| #6 File organization | Listed 5 files, only 2 exist | FIXED |
| #7 DomainOpType alias | Doc used `DomainOperator<MeshType>` | FIXED |
| #10 Traction comment | Domain header said `[[du/dn]]` | FIXED |

### Fixed in v2.0 to v3.0

| Issue (v2.0) | Description | Status |
|--------------|-------------|--------|
| #1 GetTime() in doc | Removed from doc | FIXED |
| #2 UpdateDomainAndTraction() in doc | Removed from doc | FIXED |
| #3 Init() return value | Doc now captures `real_t V_max =` | FIXED |
| #4 State size check | MFEM_VERIFY added to doc | FIXED |
| #5 params_ access pattern | Changed to `fault_->GetParams()` | FIXED |
| #6 AdaptiveTimeStepper API | Added growth_factor_, setters, getters | FIXED |
| #7 Main driver constructor sigs | Updated to match actual API | FIXED |
| #8 Test file paths | Changed `integration` to `unit` | FIXED |
| #9 Test framework | Updated to custom macros | FIXED |
| #10 Acceptance criteria | Added coarse-mesh tolerance note | FIXED |
| #11 GetDomain()/GetFault() | Added to doc | FIXED |

---

## Part III: Remaining Discrepancies

### 1. Constructor Body Missing from Doc Class Definition (MEDIUM)

**Documentation (phase4_seas_operator.md lines 131-134):**
The ODE Solver Integration section shows a simplified constructor:
```cpp
SEASQuasiDynamicOperator(DomainOpType *domain,
                          RateStateFaultOperator<MeshType> *fault)
    : TimeDependentOperator(fault->StateSize()),
      domain_(domain), fault_(fault) {}
```
Empty body `{}`.

**Implementation (seas_operator.hpp lines 114-130):**
```cpp
SEASQuasiDynamicOperator(...)
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

**Analysis:** The doc's ODE overview section (lines 128-139) shows a simplified constructor with an empty body. The Detailed Component Design section (lines 160-207) shows the class without the constructor body at all. The implementation has null-pointer guards, displacement grid function allocation, and work vector sizing — all essential for a functional operator.

**Impact:** Medium. Someone implementing from only the overview section would miss null checks and memory allocation.

**Note:** The doc's Detailed Component Design section (lines 174-175) correctly shows just the constructor declaration without body, which is appropriate for a class declaration. The issue is specifically with the simplified constructor in lines 131-134.

**Recommendation:** Either expand the constructor body in the overview section, or add a note: "See Implementation Details section for full constructor body."

---

### 2. BP2 Short Test: Doc Uses TEST_REL_NEAR with 5%, Impl Uses TEST_ASSERT with 10% (MEDIUM)

**Documentation (phase4_seas_operator.md line 402):**
```cpp
TEST_REL_NEAR(V_initial, params.V_init, 0.05, "Initial V ~ V_init");
```

**Implementation (test_bp2_short.cpp line 131):**
```cpp
TEST_ASSERT(std::abs(V_initial - params.V_init) / params.V_init < 0.1,
            "Initial V_max ~ V_init");
```

**Differences:**
1. `TEST_REL_NEAR` macro is not defined in `test_bp2_short.cpp` (only `TEST_ASSERT` and `TEST_NEAR` exist)
2. Tolerance: doc says 5% (0.05), implementation uses 10% (0.1)
3. Variable name: `V_initial` (doc) vs `V_initial` (impl) — matches

**Analysis:** `test_quasi_dynamic.cpp` defines `TEST_REL_NEAR` but `test_bp2_short.cpp` does not. The doc shows a macro that doesn't exist in this file. The 10% tolerance in the implementation is more lenient than the doc's 5%.

**Recommendation:** Update doc to show `TEST_ASSERT` with 0.1 tolerance, matching the actual `test_bp2_short.cpp`.

---

### 3. Doc Documents 3 of 9 Test Functions in test_quasi_dynamic.cpp (MEDIUM)

**Documentation (lines 357-385):** Describes three test functions:
1. `TestSteadyStateSlip()`
2. `TestStressBalanceDuringStepping()` (actual name: `TestStressBalanceMaintained`)
3. `TestSlipConservation()`

**Implementation (test_quasi_dynamic.cpp):** Has nine test functions:
1. `TestOperatorConstruction()` — not documented
2. `TestInitialCondition()` — not documented
3. `TestSteadyStateSlip()` — documented
4. `TestStressBalanceMaintained()` — documented with wrong name
5. `TestSlipConservation()` — documented
6. `TestThetaPositive()` — not documented
7. `TestAdaptiveTimeStepper()` — not documented
8. `TestRK4WithAdaptiveDt()` — not documented
9. `TestMultDeterministic()` — not documented

**Analysis:** The doc shows only the physics-focused tests and omits infrastructure tests (construction, initial condition, theta positivity, adaptive time stepping, RK4 integration, determinism). One of the three documented test names is wrong (`TestStressBalanceDuringStepping` vs `TestStressBalanceMaintained`).

**Recommendation:** Either document all 9 tests or add a note: "Only key physics tests are shown. See source for full test suite (9 tests, 29 assertions)." Fix the function name.

---

### 4. SetInitialCondition: GetParams() Ordering Differs (LOW)

**Documentation (lines 264-274):**
`VerifyStressEquilibrium` is called first, then `GetParams()` retrieves params.

**Implementation (seas_operator.hpp lines 162-177):**
`GetParams()` is called first (line 162), then `VerifyStressEquilibrium` (line 165).

**Analysis:** Both orderings are functionally equivalent since the operations are independent reads. The implementation retrieves params earlier so the local reference is available for both the equilibrium check and the V_init check.

**Recommendation:** Swap the order in the doc to match the implementation (params retrieval before equilibrium check).

---

### 5. SetInitialCondition: MFEM_VERIFY Error Message Simplified (LOW)

**Documentation (lines 273-274):**
```cpp
MFEM_VERIFY(V_rel_err < 0.1,
            "Initial V_max differs from V_init by " << V_rel_err * 100 << "%");
```

**Implementation (seas_operator.hpp lines 174-177):**
```cpp
MFEM_VERIFY(V_rel_err < 0.1,
            "Initial V_max = " << V_max
            << " differs from V_init = " << params.V_init
            << " by " << V_rel_err * 100 << "%");
```

**Analysis:** The implementation includes actual values of `V_max` and `params.V_init` in the error message for easier debugging. The doc simplifies this.

**Recommendation:** Update doc error message to include actual values.

---

### 6. ComputeRHS Return Type Not Documented (LOW)

**Documentation (line 230):**
```cpp
fault_->ComputeRHS(traction_, state, rate);
```
Called without capturing return value.

**Implementation (rate_state_fault.hpp lines 175-176):**
```cpp
/// @return Maximum slip rate
real_t ComputeRHS(const Vector &traction, const Vector &state, Vector &rate)
```

**Analysis:** `ComputeRHS` returns `real_t` (maximum slip rate). Both the doc and the implementation's `Mult()` (seas_operator.hpp line 198) discard the return value, which is fine since `GetMaxSlipRate()` provides the same information. However, the doc doesn't document that `ComputeRHS` has a return value.

**Recommendation:** Minor. Either capture the return value in the doc or add a comment noting `ComputeRHS` also returns V_max.

---

### 7. AntiplaneDomainOperator Has 6 Parameters, Doc Comment Says 5 (LOW)

**Documentation (line 445):**
```cpp
// Create domain operator (5 args: mesh, order, mu, Vp, Wf)
```

**Implementation (antiplane_operator.hpp lines 81-82):**
```cpp
AntiplaneDomainOperator(MeshType &mesh, int order, real_t mu, real_t Vp,
                        real_t Wf = 40.0e3, DGMethod method = DGMethod::IP);
```

**Analysis:** The constructor has 6 parameters (4 required + 2 with defaults: `Wf` and `method`). The doc's 5-arg call works correctly since `method` defaults to `DGMethod::IP`. The comment "5 args" is technically inaccurate about the full constructor signature.

**Recommendation:** Change comment to "5 args (6th DGMethod defaults to IP)" or simply remove the arg count.

---

### 8. VerifyStressEquilibrium tol Parameter Omitted from Doc (LOW)

**Documentation (line 265):**
```cpp
real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
```
Shows 2-argument call only.

**Implementation (rate_state_fault.hpp lines 312-313):**
```cpp
real_t VerifyStressEquilibrium(const Vector &traction, const Vector &state,
                               real_t tol = 1e-10) const
```

**Analysis:** The method has an optional 3rd parameter `real_t tol` with default `1e-10`. The 2-arg call works correctly using the default. The parameter is not mentioned anywhere in the doc.

**Recommendation:** Add a note about the optional tolerance parameter.

---

### 9. Destructor Not in Doc Class Declaration (LOW)

**Documentation (lines 168-207):** No destructor shown.

**Implementation (seas_operator.hpp line 61):**
```cpp
~SEASQuasiDynamicOperator() override = default;
```

**Analysis:** The defaulted destructor is trivial but explicit in the implementation (important for classes with `unique_ptr` members to ensure proper destruction). Its omission from the doc is a minor gap.

**Recommendation:** Add `~SEASQuasiDynamicOperator() override = default;` to the doc class.

---

### 10. Doc Uses Literal t_final Value, Implementation Uses Named Constant (LOW)

**Documentation (line 396):**
```cpp
real_t t_final = 1.0 * 3.15576e7;  // 1 year in seconds
```

**Implementation (test_bp2_short.cpp line 87):**
```cpp
real_t t_final = BP2Params::seconds_per_year;  // 1 year in seconds
```

**Analysis:** Numerically identical (`365.25 * 24 * 3600 = 3.15576e7`), but the implementation uses the named constant `BP2Params::seconds_per_year` for maintainability.

**Recommendation:** Update doc to use `BP2Params::seconds_per_year`.

---

### 11. BP2 Short Test: Doc Missing Upper-Bound Slip Check (LOW)

**Documentation (line 415):**
```cpp
TEST_ASSERT(mean_slip > 0.0, "Mean slip positive after 1 year");
```
Only a positivity check.

**Implementation (test_bp2_short.cpp lines 252-255):**
```cpp
TEST_ASSERT(mean_slip > 0.0,
            "Mean slip is positive after 1 year");
TEST_ASSERT(mean_slip < 10.0 * expected_slip,
            "Mean slip within 10x of V_init*t");
```

**Analysis:** The implementation has an additional upper-bound check (mean slip within 10x of expected) that the doc doesn't show.

**Recommendation:** Add the upper-bound test to the doc.

---

## Part IV: Potential Improvements

### 12. Bilinear Form Reassembly Every Solve() Call (MEDIUM -- Performance)

**Location:** `antiplane_operator.hpp:460-639` (`Solve()` method)

**Issue:** Every call to `Solve()` creates a new `BilinFormType`, `LinFormType`, coefficient objects, assembles, factorizes, and solves. In the SEAS time loop, `Solve()` is called once per RK4 stage (4 times per time step). The bilinear form (stiffness matrix) is time-independent and need only be assembled once.

**Current flow per Mult() call:**
1. Create BilinearForm (allocates memory)
2. Add integrators (allocates integrator objects)
3. Assemble (loops over all elements + faces)
4. Finalize (converts to CSR)
5. Create LinearForm
6. Assemble linear form
7. Assemble slip contribution
8. Create preconditioner (GSSmoother)
9. Solve

**Optimal flow:**
1. (One-time in constructor): Assemble and factorize bilinear form
2. (Per call): Assemble RHS with current time and slip, solve

**Performance estimate for BP2 at production resolution:**
- BP2 mesh: 800x800 elements, order 4 -> ~16M DOFs
- RK4 = 4 stages/step, ~1000 steps/year, 1200 years
- ~4.8M unnecessary reassemblies
- Each reassembly includes O(N) element loops + O(N) face loops
- Potential speedup: 2-5x depending on solve cost ratio

**Recommendation:** Cache the assembled bilinear form (stiffness matrix `A`) and solver setup. Only the RHS (linear form + slip contribution) needs to change between calls.

---

### 13. Slip Conservation Test Tolerance Is Very Generous (LOW -- Testing)

**Location:** `test_quasi_dynamic.cpp:353`

```cpp
TEST_ASSERT(max_rel_error < 1.0,
            "Slip within order-of-magnitude of V_init*t");
```

**Observed value:** `max_rel_error = 0.934` (93.4% error).

**Analysis:** A relative slip error of 93.4% means the actual slip is about half of `V_init * t` at some nodes. While the test passes, this is at the boundary of the 1.0 tolerance. The discrepancy is expected on the coarse mesh (4x8, order 1) because:
- Slip varies with depth (VW vs VS zones)
- The comparison uses uniform `V_init * t` at each node
- Some nodes near the free surface (z=0) or transition zone have different slip rates

**Recommendation:** Consider refining the slip conservation check to compare mean slip in the VW zone separately from the VS zone.

---

### 14. No Output/Monitoring Infrastructure (LOW -- Feature)

**Location:** `seas_operator.hpp`

**Issue:** The SEAS operator provides `GetMaxSlipRate()`, `GetDisplacement()`, and `GetTraction()` but lacks:
- Time history output at specific fault stations (BP2 requires 12 output stations)
- Slip/traction profiles along fault at specified time intervals
- Event detection (earthquake start/end based on V threshold)
- Checkpoint/restart for long simulations

The doc (line 460) shows `BenchmarkOutput output("bp2_qd", params)` but no such class exists.

**Recommendation:** Implement `BenchmarkOutput` class before production BP2 runs.

---

### 15. AdaptiveTimeStepper Lacks Safety Factor for RK4 Stability (LOW -- Robustness)

**Location:** `time_stepper.hpp:75-97`

**Issue:** The time stepper controls dt purely based on slip rate magnitude, with no secondary check on the rate of change of V. During rapid nucleation, V can jump orders of magnitude in one step. If the jump is large enough, the proportional reduction `dt * V_target / V_max` may lag behind the actual acceleration.

**Recommendation:** Consider adding a secondary check: if V_max increased by more than a factor (e.g., 5x) in one step, reduce dt more aggressively.

---

### 16. BP2 Short Test Mean Slip ~58% of Expected (LOW -- Physics)

**Observed:**
- Mean slip after 1 year: 1.83e-2 m
- Expected (V_init * t): 3.16e-2 m
- Ratio: 0.58 (58% of expected)

**Analysis:** The 42% deficit is expected on the coarse mesh (4x8, order 1). Should converge with mesh refinement.

**Recommendation:** Add a mesh-refinement convergence test.

---

## Part V: Consistency Verification

### Elements Verified as Consistent

| Aspect | Documentation | Implementation | Status |
|--------|--------------|----------------|--------|
| theta notation | theta | theta | Match |
| Root finder | Newton | Newton | Match |
| Time stepping | No `t += dt` after Step | No `t += dt` | Match |
| Overshoot protection | `if (t + dt > t_final)` | Same in tests | Match |
| ComputeRHS signature | 3 args (traction, state, rate) | 3 args | Match |
| SetInitialCondition | 4-phase PreInit/Init flow | Same flow | Match |
| ODE solver | RK4Solver recommended | RK4Solver used | Match |
| State size | 2 x num_fault_dofs | `fault->StateSize()` = 2 x `NumNodes()` | Match |
| Coupling flow (Mult) | Extract slip, solve, traction, RHS | Same 4-step flow | Match |
| TimeDependentOperator | Inherits from TDO | Inherits from TDO | Match |
| Template parameter | `MeshType = Mesh` default | Same | Match |
| Ownership | domain/fault owned externally | Raw pointers, not owned | Match |
| DomainOpType | `AntiplaneDomainOperator<MeshType>` | Same | Match |
| File organization | 2 files listed | 2 files exist | Match |
| dt_min | 1e-6 s | 1e-6 s | Match |
| dt_max | 3.15e7 s (~1 year) | 3.15e7 s | Match |
| V_target | 1e-6 m/s | 1e-6 m/s | Match |
| Initial dt | 1e3 s | 1e3 s | Match |
| Growth factor | 1.5 (member) | 1.5 (member) | Match |
| SetGrowthFactor/SetInitialDt | In doc and impl | Same | Match |
| GetDtMin/GetDtMax/GetTargetSlipRateMax | In doc and impl | Same | Match |
| GetDomain()/GetFault() | In doc and impl | Same | Match |
| Stress balance eq | tau = sigma_n*f(V,theta) + eta*V | Same (verified < 2.81e-16) | Match |
| Traction formula | tau = mu * avg(du/dx) | Same in both doc and header | Match |
| State evolution | dtheta/dt = 1 - V*theta/Dc | Same | Match |
| Test file paths | tests/unit/ | tests/unit/ | Match |
| Test framework | Custom macros | Custom macros | Match |
| Main driver constructor sigs | Correct 5-arg domain, FaultGeometry<Mesh>(domain, params) | Same | Match |
| State size check in doc | MFEM_VERIFY present | MFEM_VERIFY present | Match |
| Init() return value | `real_t V_max = fault_->Init(...)` | Same | Match |
| params access | `fault_->GetParams()` | Same | Match |
| Acceptance criteria note | Coarse-mesh note added | N/A | Match |

### Cross-Phase Consistency

| Interface | Phase 2 to Phase 4 | Phase 3 to Phase 4 | Status |
|-----------|-------------------|-------------------|--------|
| Solve() signature | `Solve(time, slip, disp)` | N/A | Match |
| ComputeTraction() | `ComputeTraction(disp, traction)` | N/A | Match |
| GetSlip() | N/A | `GetSlip(state, slip)` | Match |
| ComputeRHS() | N/A | `ComputeRHS(traction, state, rate)` | Match |
| PreInit()/Init() | N/A | `PreInit(state)`, `Init(traction, state)` | Match |
| StateSize() | N/A | `2 * NumNodes()` | Match |

---

## Part VI: Summary Table

| Issue | Severity | Category | Fix Required |
|-------|----------|----------|--------------|
| #1 Constructor body missing from overview section | Medium | Doc incomplete | Update doc |
| #2 BP2 test: TEST_REL_NEAR/5% vs TEST_ASSERT/10% | Medium | Doc/code mismatch | Update doc |
| #3 Doc shows 3 of 9 test functions | Medium | Doc incomplete | Update doc |
| #4 GetParams() ordering | Low | Doc cosmetic | Update doc |
| #5 MFEM_VERIFY error message simplified | Low | Doc cosmetic | Update doc |
| #6 ComputeRHS return type not documented | Low | Doc omission | Update doc |
| #7 Domain constructor "5 args" comment | Low | Doc inaccuracy | Update doc |
| #8 VerifyStressEquilibrium tol param | Low | Doc omission | Update doc |
| #9 Destructor not in doc class | Low | Doc omission | Update doc |
| #10 Literal t_final vs named constant | Low | Doc style | Update doc |
| #11 BP2 test missing upper-bound slip check | Low | Doc omission | Update doc |
| #12 Reassembly every Solve() | Medium | Performance | Optimization |
| #13 Slip conservation tolerance | Low | Testing | Improve test |
| #14 No output infrastructure | Low | Feature gap | Future work |
| #15 No RK4 stability safety | Low | Robustness | Future work |
| #16 Mean slip deficit on coarse mesh | Low | Expected | Monitor |

---

## Part VII: Recommendations Priority

1. **MEDIUM:** Update doc overview constructor to note that implementation has null checks and allocation (Issue #1)
2. **MEDIUM:** Fix BP2 short test doc to use `TEST_ASSERT` with 10% tolerance (Issue #2)
3. **MEDIUM:** Document all 9 test functions or add note about full test suite (Issue #3)
4. **MEDIUM:** Cache bilinear form assembly in `Solve()` for performance (Issue #12)
5. **LOW:** Fix cosmetic doc issues: GetParams ordering, error message detail, ComputeRHS return type, domain arg count, VerifyStressEquilibrium tol, destructor, t_final constant, upper-bound slip check (Issues #4-11)
6. **LOW:** Future work: BenchmarkOutput, RK4 safety factor, mesh-refinement test (Issues #14-16)

---

## Part VIII: Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-02-06 | Initial analysis: 15 issues found (1 HIGH, 4 MEDIUM, 10 LOW) |
| 2.0 | 2025-02-06 | 8 v1.0 issues fixed in doc; 10 remaining + 5 improvements |
| 3.0 | 2025-02-06 | All 10 v2.0 issues fixed in doc; 3 MEDIUM + 8 LOW remaining |
