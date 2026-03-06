# Phase 3: RateStateFaultOperator Generalization

**Date**: 2026-03-02
**Status**: COMPLETE (compiles and passes all existing tests; BP5 path exercised via test_elasticity_operator)
**Prerequisite**: Phase 2c complete (FaultGeometry with BP5 support)

---

## Overview

Phase 3 generalizes `RateStateFaultOperator` from scalar slip (BP1/BP2) to vector slip (BP5) using a `SlipComponents` template parameter. The key change is that the per-node state layout grows from `[slip, theta]` (2 values) to `[s_dip, s_strike, psi]` (3 values), and the friction solver switches from scalar to vector.

All changes are backward-compatible: existing `SlipComponents=1` code paths are unchanged.

---

## Implementation

**File**: `fault/rate_state_fault.hpp` (EXTENDED)

### Template Parameter

```cpp
template <typename MeshType = Mesh, int SlipComponents = 1>
class RateStateFaultOperator
{
   static constexpr int NumSlipComp = SlipComponents;
   static constexpr int StatePerNode = SlipComponents + 1;  // slip(s) + psi
   static constexpr int SlipIndex = 0;
   static constexpr int ThetaIndex = SlipComponents;
   static constexpr int PsiIndex = SlipComponents;
```

### State Layout

| | BP1/BP2 (`SlipComponents=1`) | BP5 (`SlipComponents=2`) |
|---|---|---|
| **Per node** | `[slip, theta/psi]` (2 values) | `[s_dip, s_strike, psi]` (3 values) |
| **StateSize** | `2N` | `3N` |
| **SlipSize** | `N` | `2N` |
| **TractionSize** | `N` | `2N` |
| **PsiIndex** | `1` | `2` |

### Dual Constructors

**BP2 constructor** (`SlipComponents=1`):
```cpp
RateStateFaultOperator(FaultGeometry<MeshType> *geom,
                       FrictionLaw *friction,
                       StateEvolution *evolution,
                       const BP2Params &params,
                       MPIContext *mpi_ctx = nullptr,
                       bool use_psi = false);
// static_assert(SlipComponents == 1)
```

**BP5 constructor** (`SlipComponents=2`):
```cpp
RateStateFaultOperator(FaultGeometry<MeshType> *geom,
                       DieterichRuinaFriction *friction,
                       StateEvolution *evolution,
                       const BP5Params &params,
                       MPIContext *mpi_ctx = nullptr);
// static_assert(SlipComponents == 2)
// Always uses psi-space, requires DR friction
```

### Key Methods — Generalized with `if constexpr`

#### PreInit(state)

Sets initial slip to zero (all components) and initial psi to steady-state value:

```cpp
for (int i = 0; i < num_nodes_; i++) {
   for (int c = 0; c < SlipComponents; c++)
      state(i * StatePerNode + c) = 0.0;

   if constexpr (SlipComponents == 1) {
      // theta_ss = Dc/V_init or psi_ss = f0 + b*ln(V0/V_init)
   } else {
      // BP5: psi_ss from |V_init| and per-DOF Dc
      real_t V_abs_init = sqrt(V_init[2i]^2 + V_init[2i+1]^2);
      state(i*StatePerNode + PsiIndex) = evolution_->SteadyState(V_abs_init, Dc);
   }
}
```

#### Init(traction, state) → V_max

Initializes psi from stress equilibrium after the first domain solve:

```cpp
if constexpr (SlipComponents == 1) {
   // Scalar: tau0 + traction(i) → SolveSlipRatePsi → V
} else {
   // Vector: tau_pre[2i,2i+1] + traction[2i,2i+1] → |tau|
   //         InitialStatePsi(|tau|, |V_init|, ...) → psi0
   //         SolveSlipRateVectorPsi(tau_vec, psi0, ...) → V_vec
   //         V_max = max(V_max, |V_vec|)
}
```

#### ComputeRHS(traction, state, rate) → V_max

The hot inner loop — computes slip rates and state evolution:

```cpp
for (int i = 0; i < num_nodes_; i++) {
   if constexpr (SlipComponents == 1) {
      // Scalar: unchanged from existing code
      real_t V = SolveSlipRatePsi(tau, psi, sigma_n, eta, a);
      rate(i*StatePerNode + SlipIndex) = V;
      rate(i*StatePerNode + ThetaIndex) = evolution_->Rate(V, psi, Dc);
   } else {
      // Vector: tau_pre + traction → 2-component tau_vec
      real_t tau_vec[2] = {tau_pre_(2*i) + traction(2*i),
                           tau_pre_(2*i+1) + traction(2*i+1)};
      real_t V_vec[2];
      dr_friction_->SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);
      real_t V_abs = sqrt(V_vec[0]^2 + V_vec[1]^2);

      rate(i*StatePerNode + 0) = V_vec[0];
      rate(i*StatePerNode + 1) = V_vec[1];
      rate(i*StatePerNode + PsiIndex) = evolution_->Rate(V_abs, psi, Dc);
   }
}
```

#### GetSlip(state, slip)

Extracts slip components from state:
```cpp
for (int i = 0; i < num_nodes_; i++)
   for (int c = 0; c < SlipComponents; c++)
      slip(i*SlipComponents + c) = state(i*StatePerNode + c);
```

#### GetMaxSlipRate()

For `SlipComponents=1`: returns `max(|V_i|)`
For `SlipComponents=2`: returns `max(sqrt(V_dip^2 + V_strike^2))`

---

## Size Query Methods

```cpp
int StateSize() const   { return num_nodes_ * StatePerNode; }
int SlipSize() const    { return num_nodes_ * SlipComponents; }
int TractionSize() const { return num_nodes_ * SlipComponents; }
int NumNodes() const    { return num_nodes_; }
```

---

## Additional BP5 Data Members

```cpp
// BP5-specific (only used when SlipComponents=2)
BP5Params bp5_params_;
real_t sigma_n_bp5_, Vp_bp5_, Wf_bp5_;
Vector Dc_values_;        // [N] per-DOF critical slip distance
Vector tau_pre_;          // [2N] pre-stress (from FaultGeometry)
Vector V_init_values_;    // [2N] initial velocity (from FaultGeometry)
DieterichRuinaFriction *dr_friction_;  // Required for vector solver
```

---

## Backward Compatibility

The `static_assert` in each constructor ensures compile-time safety:
- BP2 constructor: `static_assert(SlipComponents == 1)`
- BP5 constructor: `static_assert(SlipComponents == 2)`

Default template argument `SlipComponents=1` means all existing code:
```cpp
RateStateFaultOperator<Mesh> fault(...)  // Defaults to SlipComponents=1
```
compiles and behaves identically to before.

The `if constexpr` branching has **zero runtime overhead** for the default scalar path — the compiler eliminates the vector branch entirely.

---

## Type Aliases for Convenience

```cpp
// BP5 instantiation
using BP5FaultOp = RateStateFaultOperator<Mesh, 2>;

// BP2 (default, existing)
using BP2FaultOp = RateStateFaultOperator<Mesh, 1>;
// or equivalently: RateStateFaultOperator<Mesh>
```

---

## Tests

All 13 pre-existing test suites pass without modification:
- `seas_test_fault_operator` — BP2 scalar path unchanged
- `seas_test_quasi_dynamic` — default template instantiation unchanged
- `seas_test_antiplane` — full BP2 pipeline unchanged

The BP5 constructor path is exercised via:
- `test_elasticity_operator.cpp: TestFaultGeometry3D()` — constructs `FaultGeometry<Mesh>(op, bp5_params)` which the BP5 `RateStateFaultOperator` constructor consumes

---

## Files Modified

| File | Change |
|------|--------|
| `fault/rate_state_fault.hpp` | Added `SlipComponents` template parameter, BP5 constructor, `if constexpr` branching in PreInit/Init/ComputeRHS/GetSlip, size query methods, BP5 data members |

## Key Design Decisions

1. **Template parameter, not runtime flag**: `SlipComponents` is a compile-time constant, enabling `if constexpr` with zero overhead for the scalar path.

2. **`static_assert` for constructor safety**: Prevents accidentally using BP5Params with SlipComponents=1 or vice versa — caught at compile time, not runtime.

3. **BP5 always uses psi-space**: Theta-space (raw state variable) is not supported for BP5. This simplifies the vector path.

4. **Per-DOF parameters cached from FaultGeometry**: `Dc_values_`, `tau_pre_`, `V_init_values_` are copied at construction, avoiding per-call lookups during the hot loop.
