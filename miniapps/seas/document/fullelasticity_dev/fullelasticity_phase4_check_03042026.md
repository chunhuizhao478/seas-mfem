# Phase 4 Implementation Review: SEASQuasiDynamicOperator Wiring and Integration

**Date**: 2026-03-04
**Last updated**: 2026-03-04 (revision 3 — all bugs and code quality issues resolved)
**Reviewing**: `fullelasticity_phase4_plan_03022026.md` vs. current implementation
**Files checked**: `solver/seas_operator.hpp`, `solver/time_stepper.hpp`,
`tests/unit/test_bp5_integration.cpp`, `tests/unit/test_quasi_dynamic.cpp`,
`fault/rate_state_fault.hpp`, `domain/elasticity_operator.hpp`,
`fault/fault_geometry.hpp`, `Makefile`
**References**: Tandem `app/form/SeasQDOperator.cpp`, Tandem `app/form/SeasQDOperator.h`,
Phase 4 plan, SCEC BP5-QD spec

---

## 1. Overall Assessment

Phase 4 implementation is complete and correct. **All bugs and code quality
issues have been resolved** across revisions 1–3:

- **§2.1** DormandPrinceRK45 diagnostic now uses configurable `state_per_node_`
  (default 2 for backward compatibility, settable via `SetStatePerNode(3)` for BP5).
- **§2.2** `PBP5SEASOp` alias now correctly guarded by `#ifdef MFEM_USE_MPI`.
- **§3.1** All type aliases from the plan are now defined: `BP2SEASOp`,
  `BP5DomainOp`, `BP5SEASOp`, and `PBP5SEASOp`.
- **§4.2** `TestBP5ShortRK45Run` now calls `rk45.SetStatePerNode(3)`.
- **§4.3** Clarifying comment added for V_max initialization flow.

**0 open bugs, 0 open code quality issues.** Only deferred items remain
(conditional solve optimization, parallel diagnostic enhancement — both INFO).

**Cumulative changes across all revisions**:
- `time_stepper.hpp`: Added `SetStatePerNode()` method, `state_per_node_` member
  (default 2), diagnostic now uses `state_per_node_` for DOF/variable identification
- `seas_operator.hpp`: Added `BP2SEASOp`, `PBP5SEASOp` type aliases;
  `PBP5SEASOp` wrapped in `#ifdef MFEM_USE_MPI`; clarifying comment on V_max
- `test_bp5_integration.cpp`: Added `rk45.SetStatePerNode(3)` call

---

## 2. Bugs / Potential Bugs

### 2.1 ~~[LOW] DormandPrinceRK45 diagnostic hardcodes 2-component state layout~~

**Status**: COMPLETED

**File**: `solver/time_stepper.hpp:177, 421-424, 504`

The diagnostic now uses a configurable `state_per_node_` member:

```cpp
void SetStatePerNode(int spn) { state_per_node_ = spn; }  // line 177
```

```cpp
int dof = worst_idx / state_per_node_;                     // line 421
int comp = worst_idx % state_per_node_;                    // line 422
const char *comp_name = (comp == state_per_node_ - 1)      // line 423
                           ? "(theta/psi)" : "(slip)";     // line 424
```

```cpp
int state_per_node_ = 2;  // line 504 — default for BP2 backward compat
```

Default value of 2 maintains backward compatibility with existing BP2 code.
For BP5, callers set `rk45.SetStatePerNode(3)` before time stepping.

### 2.2 ~~[LOW] PBP5SEASOp alias lacks `#ifdef MFEM_USE_MPI` guard~~

**Status**: COMPLETED

**File**: `solver/seas_operator.hpp:228-233`

The `PBP5SEASOp` alias is now correctly wrapped in an MPI guard:

```cpp
#ifdef MFEM_USE_MPI
// Parallel BP5 type alias
using PBP5SEASOp = SEASQuasiDynamicOperator<ParMesh,
                      ElasticityDomainOperator<ParMesh>,
                      RateStateFaultOperator<ParMesh, 2>>;
#endif
```

This is consistent with `seas_types.hpp` which also guards `ParMesh`
references with `#ifdef MFEM_USE_MPI`.

---

## 3. Deviations from Plan

### 3.1 ~~Missing BP2 and parallel type aliases~~

**Status**: COMPLETED

`seas_operator.hpp:220-230` now defines all aliases from the plan:

```cpp
// BP2 type alias (uses default template arguments)
using BP2SEASOp = SEASQuasiDynamicOperator<Mesh>;            // line 221

// BP5 type aliases
using BP5DomainOp = ElasticityDomainOperator<Mesh>;          // line 224
using BP5SEASOp   = SEASQuasiDynamicOperator<Mesh, BP5DomainOp, BP5FaultOp>;  // line 225

// Parallel BP5 type alias
using PBP5SEASOp = SEASQuasiDynamicOperator<ParMesh,         // line 228-230
                      ElasticityDomainOperator<ParMesh>,
                      RateStateFaultOperator<ParMesh, 2>>;
```

`BP2FaultOp` and `BP5FaultOp` remain in `rate_state_fault.hpp:798-799`
(appropriate placement alongside the class they alias).

### 3.2 More integration tests than planned (positive deviation)

**Plan**: 2 tests (TestBP5FullStackInit, TestBP5ShortRun) plus optional
scalar-vector consistency.

**Code**: 7 tests in `test_bp5_integration.cpp`:
1. TestBP5FullStackConstruction
2. TestBP5SetInitialCondition
3. TestBP5MultConsistency
4. TestBP5StateRoundTrip
5. TestBP5ShortRK4Run
6. TestBP5ShortRK45Run (not in plan)
7. TestBP5StressBalanceDuringTimeStep (not in plan)

This is a positive deviation. The additional tests for RK45 integration and
stress balance during time stepping provide stronger verification.

---

## 4. Code Quality Issues

### 4.1 [INFO] No conditional solve optimization

**File**: `solver/seas_operator.hpp` (Mult, lines 195-218)

Tandem's `SeasQDOperator` has an `update_internal_state()` method that can
skip the domain solve when state hasn't changed since the last RHS evaluation.
This is important for output/checkpoint operations that need displacement or
traction without recomputing.

The MFEM implementation always performs the full solve in `Mult()`. This is
correct and simpler. The optimization is not needed until long-running
simulations with frequent output become a bottleneck.

**Recommendation**: No action needed now. Consider adding a caching mechanism
(dirty flag on state) if profiling shows the domain solve is the bottleneck
during output steps.

### 4.2 ~~[INFO] TestBP5ShortRK45Run does not call SetStatePerNode(3)~~

**Status**: COMPLETED

**File**: `tests/unit/test_bp5_integration.cpp:412`

The test now correctly calls `rk45.SetStatePerNode(3)`:

```cpp
DormandPrinceRK45 rk45;
rk45.SetAbsTol(1e-7);
rk45.SetDt(1e3);
rk45.SetDtMax(1e6);
rk45.SetStatePerNode(3);  // BP5: 3 components per node
rk45.Init(*fix.seas_op);
```

### 4.3 ~~[INFO] SetInitialCondition: V_max from Init() immediately overwritten~~

**Status**: COMPLETED

**File**: `solver/seas_operator.hpp:166`

Clarifying comment added:

```cpp
// V_max is recomputed below after the verification re-solve
real_t V_max = fault_->Init(traction_, state);
```

### 4.4 [INFO] DormandPrinceRK45 diagnostic worst-DOF inaccuracy in parallel

**File**: `time_stepper.hpp:369-384, 419-434`

In parallel L-infinity mode, `worst_idx` is tracked locally per rank
(line 376), but `err_norm` is globally reduced via `GlobalMax` (line 384).
Only rank 0 prints the diagnostic (line 419), so it reports its own local
worst DOF alongside the global `err_norm`. If the globally worst DOF is on
a non-root rank, the printed DOF index and per-DOF values are misleading.

**Impact**: None. This only affects verbose diagnostic output (`SetVerbose(true)`)
in parallel runs. The error norm, accept/reject decisions, and dt adjustments
are all correct. The same pattern exists in the 2-norm path.

**Recommendation**: No action needed for Phase 4. If parallel debugging
becomes important, consider using `MPI_MAXLOC` to identify the global
worst rank and DOF together.

---

## 5. Correctness Verification Against Plan and Tandem

### 5.1 Template generalization: CORRECT

`seas_operator.hpp:48-50`:
```cpp
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneDomainOperator<MeshType>,
          typename FaultOpType = RateStateFaultOperator<MeshType>>
class SEASQuasiDynamicOperator : public TimeDependentOperator
```

Matches plan exactly. Default arguments maintain backward compatibility
with existing BP2 code that uses `SEASQuasiDynamicOperator<Mesh>`.

### 5.2 Work vector sizing: CORRECT

`seas_operator.hpp:144-146`:
```cpp
slip_.SetSize(fault_->SlipSize());         // N for BP2, 2N for BP5
traction_.SetSize(fault_->TractionSize()); // N for BP2, 2N for BP5
```

`TimeDependentOperator` base size set from `fault->StateSize()` at line 134.
Matches plan specification.

### 5.3 Mult() coupling flow: CORRECT

`seas_operator.hpp:195-218`:
```
1. fault_->GetSlip(state, slip_)          ✓ Extract slip from state
2. domain_->Solve(t, slip_, *u_gf_)       ✓ Solve domain with slip BC
3. domain_->ComputeTraction(...)           ✓ Get fault traction
4. fault_->ComputeRHS(traction_, state, rate)  ✓ Compute rates
```

Matches both plan and Tandem's `rhs()` pattern. The branching between scalar
and vector happens inside the domain and fault operators, not here.

### 5.4 SetInitialCondition(): CORRECT

`seas_operator.hpp:149-192` follows Tandem's two-phase init:

```
Phase 1: fault_->PreInit(state)                       ✓ slip=0, psi=placeholder
Phase 2: GetSlip → Solve → ComputeTraction            ✓ domain solve with zero slip
Phase 3: fault_->Init(traction_, state)                ✓ compute psi from equilibrium
Phase 4: Re-solve → ComputeRHS → VerifyStressEquilibrium  ✓ verify
```

The MFEM implementation adds an explicit verification phase (re-solve and
check stress equilibrium) that Tandem does not have. This is a good addition
for debugging.

### 5.5 Backward compatibility: CORRECT

Existing tests use `SEASQuasiDynamicOperator<Mesh>` which defaults to
`<Mesh, AntiplaneDomainOperator<Mesh>, RateStateFaultOperator<Mesh>>` — the
original BP2 configuration. No existing code needs modification.

Verified by plan's "Regression" section: all existing test suites pass with
the template changes.

### 5.6 TimeDependentOperator base: CORRECT

`seas_operator.hpp:134`: `TimeDependentOperator(fault->StateSize())`

The operator height equals the fault state size (2N for BP2, 3N for BP5),
which is the correct ODE system dimension.

### 5.7 DormandPrinceRK45 numerical correctness: CORRECT

The DOPRI5(4) integrator at `time_stepper.hpp:130-577` is numerically correct:
- Butcher tableau coefficients match standard DOPRI5(4) values ✓
- FSAL property correctly implemented (k_[6] → k_[0]) ✓
- Error estimation uses weighted norm (L-inf or RMS) ✓
- NaN detection at each stage with aggressive rejection ✓
- Parallel error norm reduction via MPIContext ✓
- PI controller with PETSc-style safety factor ✓
- Configurable `state_per_node_` for diagnostic output ✓

### 5.8 AdaptiveTimeStepper: CORRECT

`time_stepper.hpp:39-128` implements slip-rate-based adaptive time stepping
with the standard V_target strategy. Not used in BP5 integration tests but
available for BP5 driver programs.

### 5.9 Tandem coupling pattern: CORRECT

Compared against Tandem's `SeasQDOperator.cpp`:
- MFEM's `Mult()` ↔ Tandem's `rhs()`: same flow (skip ghost exchange in serial) ✓
- MFEM's `SetInitialCondition()` ↔ Tandem's `initial_condition()`: same two-phase pattern ✓
- MFEM stores displacement as GridFunction ↔ Tandem uses linear solver's `x()` ✓
- MFEM uses `GetMaxSlipRate()` ↔ Tandem tracks VMax in friction operator ✓

### 5.10 Makefile: CORRECT

- `TEST_BP5_INTEGRATION_SRC` at line 56 ✓
- `TEST_BP5_INTEGRATION_OBJ` at line 77 ✓
- `seas_test_bp5_integration` build target at lines 309-310 ✓
- `test-bp5-integration` phony target at lines 460-461 ✓
- `.PHONY` includes `test-bp5-integration` at line 226 ✓
- Added to `SEQ_MINIAPPS` at line 135 ✓
- Added to `test:` target at line 404 ✓
- Object file rule at lines 389-391 depends on `$(SEAS_HEADERS)` ✓

### 5.11 Type aliases: CORRECT

All aliases from the plan are defined:
- `BP2FaultOp`, `BP5FaultOp` in `rate_state_fault.hpp:798-799` ✓
- `BP2SEASOp` in `seas_operator.hpp:221` ✓
- `BP5DomainOp`, `BP5SEASOp` in `seas_operator.hpp:224-225` ✓
- `PBP5SEASOp` in `seas_operator.hpp:230-232` (guarded by `#ifdef MFEM_USE_MPI`) ✓

---

## 6. Test Coverage Assessment

### 6.1 Integration Tests (test_bp5_integration.cpp)

| Test | What It Verifies | Status |
|------|-----------------|--------|
| TestBP5FullStackConstruction | Template instantiation, StateSize=3N, SlipSize=2N, pointer consistency | ✓ |
| TestBP5SetInitialCondition | 4-phase init completes, V_max>0, slip=0, traction finite | ✓ |
| TestBP5MultConsistency | Rate vector finite, \|V\|>0 per node, deterministic (two calls identical) | ✓ |
| TestBP5StateRoundTrip | GetSlip=2N, GetTheta=N positive, traction=2N, displacement finite | ✓ |
| TestBP5ShortRK4Run | 10 RK4 steps, state finite, V_max in (0,1), slip accumulated, theta>0 | ✓ |
| TestBP5ShortRK45Run | 5 accepted DOPRI5(4) steps within 50 attempts, state finite, dt>0 | ✓ |
| TestBP5StressBalanceDuringTimeStep | Stress equilibrium < 1e-4 after 3 RK4 steps | ✓ |

### 6.2 Test Fixture Quality

The `BP5IntegrationFixture` struct correctly constructs the full BP5 stack:
1. 3D hex mesh (2x2x2 elements, 50km x 60km x 40km) with boundary attributes ✓
2. `ElasticityDomainOperator<Mesh>` with IP method ✓
3. `FaultGeometry<Mesh>` constructed with BP5Params ✓
4. `DieterichRuinaFriction` with BP5 constants ✓
5. `AgingLawPsi` with BP5 parameters ✓
6. `RateStateFaultOperator<Mesh, 2>` (BP5 vector fault) ✓
7. `BP5SEASOp` (full SEAS operator using type alias) ✓

The graceful skip on `nf == 0` (no fault faces found) prevents test crashes
on mesh configurations that don't produce fault-crossing faces.

### 6.3 Remaining Test Gaps

| Gap | Priority | Notes |
|-----|----------|-------|
| Scalar-vector consistency test | LOW | Plan mentions as optional. Would require comparing antiplane 2D result against 3D elasticity with purely antiplane loading. Non-trivial to set up correctly. |
| AdaptiveTimeStepper with BP5 | LOW | Existing BP2 tests cover AdaptiveTimeStepper. BP5 uses DormandPrinceRK45 directly. |
| SetDisplacement() round-trip | INFO | Checkpoint feature, likely tested in checkpoint test suite. |
| Parallel BP5 integration (PBP5SEASOp) | INFO | Planned for parallel phase, not Phase 4 scope. |

---

## 7. Summary Table

| Item | Status | Severity | Action |
|------|--------|----------|--------|
| DormandPrinceRK45 diagnostic format | **COMPLETED** | ~~LOW~~ | SetStatePerNode() added, default 2 for backward compat |
| PBP5SEASOp MPI guard | **COMPLETED** | ~~LOW~~ | Wrapped in `#ifdef MFEM_USE_MPI` |
| Missing BP2SEASOp, PBP5SEASOp aliases | **COMPLETED** | ~~INFO~~ | All aliases defined in seas_operator.hpp |
| Test SetStatePerNode(3) | **COMPLETED** | ~~INFO~~ | Added to TestBP5ShortRK45Run |
| V_max clarifying comment | **COMPLETED** | ~~INFO~~ | Comment added at seas_operator.hpp:166 |
| No conditional solve optimization | Deferred | INFO | No action needed now |
| Parallel diagnostic worst-DOF | Deferred | INFO | Only affects verbose output in MPI runs |
| Template generalization | Correct | -- | -- |
| Work vector sizing | Correct | -- | -- |
| Mult() coupling flow | Correct | -- | -- |
| SetInitialCondition() flow | Correct | -- | -- |
| Backward compatibility | Correct | -- | -- |
| TimeDependentOperator base | Correct | -- | -- |
| DormandPrinceRK45 numerics | Correct | -- | -- |
| AdaptiveTimeStepper | Correct | -- | -- |
| Tandem coupling pattern match | Correct | -- | -- |
| Type aliases | Correct | -- | -- |
| Makefile | Correct | -- | -- |
| Integration tests (7 tests) | Correct | -- | More than planned (positive) |
| Test fixture (BP5IntegrationFixture) | Correct | -- | -- |
| Scalar-vector consistency | Not tested | LOW | Optional per plan |
