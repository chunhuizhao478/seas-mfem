# Phase 4: SEASQuasiDynamicOperator Wiring and Integration

**Date**: 2026-03-02
**Status**: Template changes COMPLETE (compiles with all existing tests); Integration testing REMAINING
**Prerequisite**: Phase 3 complete (RateStateFaultOperator generalized)

---

## Overview

Phase 4 wires the full BP5 stack together through `SEASQuasiDynamicOperator`. The operator was already templated on `DomainOpType` from Phase 1. Phase 4 adds `FaultOpType` as a third template parameter so that the SEAS operator can work with both scalar (BP2) and vector (BP5) fault operators.

---

## Implementation

**File**: `solver/seas_operator.hpp` (EXTENDED)

### Template Generalization

```cpp
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneDomainOperator<MeshType>,
          typename FaultOpType = RateStateFaultOperator<MeshType>>
class SEASQuasiDynamicOperator : public TimeDependentOperator
{
public:
   SEASQuasiDynamicOperator(DomainOpType *domain,
                             FaultOpType *fault,
                             MPIContext *mpi_ctx = nullptr);

   void SetInitialCondition(Vector &state);
   void Mult(const Vector &state, Vector &rate) const override;
   // ...

private:
   DomainOpType *domain_;
   FaultOpType *fault_;
   mutable Vector slip_;
   mutable Vector traction_;
};
```

### Work Vector Sizing

Work vectors are sized using the fault operator's query methods:

```cpp
SEASQuasiDynamicOperator(...) : TimeDependentOperator(fault->StateSize()), ...
{
   slip_.SetSize(fault_->SlipSize());         // N for BP2, 2N for BP5
   traction_.SetSize(fault_->TractionSize()); // N for BP2, 2N for BP5
}
```

### Mult() Flow (Unchanged Logic)

The coupling flow in each `Mult()` call is unchanged:

```
1. fault_->GetSlip(state, slip_)          Extract slip from state
2. domain_->Solve(t, slip_, *u_gf_)       Solve domain with slip BC
3. domain_->ComputeTraction(*u_gf_, slip_, traction_)  Get fault traction
4. fault_->ComputeRHS(traction_, state, rate)          Compute rates
```

This works for both BP2 and BP5 because:
- `GetSlip` extracts `SlipComponents` per node (1 or 2)
- `Solve` takes a slip vector of the correct size (N or 2N)
- `ComputeTraction` returns traction of the same size
- `ComputeRHS` handles scalar/vector branching internally via `if constexpr`

### SetInitialCondition() Flow (Unchanged Logic)

Two-phase initialization following Tandem:

```
1. fault_->PreInit(state)                     Set slip=0, psi=placeholder
2. fault_->GetSlip(state, slip_)              Get zero slip
3. domain_->Solve(0.0, slip_, *u_gf_)         Solve with zero slip
4. domain_->ComputeTraction(*u_gf_, slip_, traction_)  Initial traction
5. fault_->Init(traction_, state)              Compute psi from equilibrium
6. Verify: re-solve, ComputeRHS, check V_max ≈ V_init
```

---

## BP5 Type Aliases

```cpp
// Full BP5 stack
using BP5DomainOp = ElasticityDomainOperator<Mesh>;
using BP5FaultOp  = RateStateFaultOperator<Mesh, 2>;
using BP5SEASOp   = SEASQuasiDynamicOperator<Mesh, BP5DomainOp, BP5FaultOp>;

// Existing BP2 stack (unchanged, uses defaults)
using BP2DomainOp = AntiplaneDomainOperator<Mesh>;
using BP2FaultOp  = RateStateFaultOperator<Mesh, 1>;
using BP2SEASOp   = SEASQuasiDynamicOperator<Mesh, BP2DomainOp, BP2FaultOp>;
// or equivalently: SEASQuasiDynamicOperator<Mesh>

// Parallel variants
using PBP5SEASOp = SEASQuasiDynamicOperator<ParMesh,
                      ElasticityDomainOperator<ParMesh>,
                      RateStateFaultOperator<ParMesh, 2>>;
```

---

## Integration Testing (REMAINING)

### Test: Full BP5 Initialization

**File**: `tests/unit/test_bp5_integration.cpp` (TO CREATE)

```cpp
void TestBP5FullStackInit() {
   // 1. Create 3D mesh with fault at x1=0
   Mesh mesh = CreateTestMesh3D(2, 2, 2, 50e3, 60e3, 40e3);

   // 2. Construct domain operator
   BP5Params params;
   ElasticityDomainOperator<Mesh> domain(mesh, 1,
      params.lambda(), params.mu(), params.Vp, params.Wf, params.lf);

   // 3. Construct fault geometry
   FaultGeometry<Mesh> geom(domain, params);

   // 4. Construct fault operator
   DieterichRuinaFriction friction;
   AgingLawPsi evolution(params.f0, params.b, params.V0);
   RateStateFaultOperator<Mesh, 2> fault(&geom, &friction, &evolution, params);

   // 5. Construct SEAS operator
   SEASQuasiDynamicOperator<Mesh,
      ElasticityDomainOperator<Mesh>,
      RateStateFaultOperator<Mesh, 2>> seas(&domain, &fault);

   // 6. Initialize
   Vector state(fault.StateSize());
   seas.SetInitialCondition(state);

   // 7. Verify V_max ≈ max(V_init components)
   TEST_ASSERT(seas.GetMaxSlipRate() > 0, "V_max > 0 after init");
}
```

### Test: Short Time Integration

```cpp
void TestBP5ShortRun() {
   // Construct full stack (as above)
   // Run 10 RK4 steps
   RK4Solver ode;
   ode.Init(seas);
   real_t dt = 1e4;  // 10,000 seconds
   for (int i = 0; i < 10; i++) {
      ode.Step(state, t, dt);
   }
   // Verify state is finite and V_max is reasonable
}
```

### Test: Scalar-Vector Consistency

For pure antiplane loading on a 3D mesh (slip only in x2 direction, no x3 slip), the 3D elasticity result should match the 2D antiplane result to within discretization error.

---

## Regression

All existing test suites pass with the template changes:
- `seas_test_quasi_dynamic` — uses default template `SEASQuasiDynamicOperator<Mesh>` → backward compatible
- `seas_test_fault_operator` — uses `RateStateFaultOperator<Mesh>` (SlipComponents=1 default) → unchanged
- `seas_test_antiplane` — full BP2 pipeline → unchanged

---

## Files Modified

| File | Change |
|------|--------|
| `solver/seas_operator.hpp` | Added `FaultOpType` template parameter, sized work vectors via `fault_->SlipSize()` and `fault_->TractionSize()` |

## Files To Create (Integration Tests)

| File | Purpose |
|------|---------|
| `tests/unit/test_bp5_integration.cpp` | Full BP5 stack initialization + short time integration |

---

## Key Design Decisions

1. **Three template parameters**: `MeshType`, `DomainOpType`, `FaultOpType`. Default values maintain full backward compatibility with existing code.

2. **No runtime dispatch**: Unlike the original plan's suggestion of querying `NumSlipComponents()` at runtime, the template approach gives compile-time resolution with zero overhead.

3. **Work vectors sized from fault operator**: `slip_.SetSize(fault_->SlipSize())` instead of hardcoded sizes — works for both scalar (N) and vector (2N) cases.

4. **Unchanged coupling flow**: `Mult()` and `SetInitialCondition()` have identical logic for BP2 and BP5. The branching happens inside the domain and fault operators, not in the SEAS operator.

---

## Remaining Work Summary

| Item | Status | Notes |
|------|--------|-------|
| Template generalization | COMPLETE | `FaultOpType` template parameter added |
| Work vector sizing | COMPLETE | Uses `SlipSize()`, `TractionSize()` |
| Backward compatibility | COMPLETE | All 13 test suites pass |
| Integration test (full stack init) | TODO | `test_bp5_integration.cpp` |
| Integration test (short run) | TODO | RK4 time stepping with BP5 |
| Parallel integration test | TODO | `ParMesh` + MPI |
| BP5 driver program | TODO | `bp5/bp5_verification.cpp` (Phase 5) |
| BP5 output / I/O | TODO | On-fault/off-fault time series (Phase 5) |
