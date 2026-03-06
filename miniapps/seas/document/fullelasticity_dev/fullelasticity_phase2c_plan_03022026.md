# Phase 2c: FaultGeometry 3D Extension

**Date**: 2026-03-02
**Status**: COMPLETE (tested via test_elasticity_operator.cpp, TestFaultGeometry3D)
**Prerequisite**: Phase 2b complete (ElasticityDomainOperator providing 2D fault coordinates)

---

## Overview

Phase 2c extends `FaultGeometry` to support 3D (BP5) fault surfaces with spatially varying parameters. The existing BP2 constructor handles 1D depth-dependent parameters `a(z)`, `eta(z)`. The new BP5 constructor handles 2D spatially varying parameters `a(x2, x3)`, `L(x2, x3)`, `tau_pre(x2, x3)`, `V_init(x2, x3)`.

---

## Implementation

**File**: `fault/fault_geometry.hpp` (EXTENDED)

### New Constructor for BP5

```cpp
FaultGeometry(DomainOperator<MeshType> &domain_op, const BP5Params &params,
              MPIContext *mpi_ctx = nullptr)
   : bp5_params_(params), mpi_ctx_(mpi_ctx), is_bp5_(true)
{
   num_fault_dofs_ = domain_op.GetNumFaultDOFs();
   // ... parallel setup ...

   // Get 2D fault coordinates from domain operator
   domain_op.GetFaultCoords2D(coords_x2_, coords_x3_);

   // Store depths for compatibility
   depths_.SetSize(num_fault_dofs_);
   for (int i = 0; i < num_fault_dofs_; i++)
      depths_(i) = coords_x3_(i);

   // Precompute per-DOF parameters using BP5 2D functions
   ComputeBP5Params();
}
```

### Per-DOF Parameter Computation

`ComputeBP5Params()` precomputes the following arrays at construction time:

```cpp
void ComputeBP5Params() {
   a_values_.SetSize(num_fault_dofs_);
   eta_values_.SetSize(num_fault_dofs_);
   dc_values_.SetSize(num_fault_dofs_);
   tau_pre_.SetSize(2 * num_fault_dofs_);
   V_init_vec_.SetSize(2 * num_fault_dofs_);

   for (int i = 0; i < num_fault_dofs_; i++) {
      real_t x2 = coords_x2_(i), x3 = coords_x3_(i);

      a_values_(i)  = bp5_params_.a_of_x2_x3(x2, x3);
      eta_values_(i) = bp5_params_.eta();
      dc_values_(i)  = bp5_params_.Dc_of_x2_x3(x2, x3);

      real_t tau[2];
      bp5_params_.tau0_vec(x2, x3, tau);
      tau_pre_(2*i)     = tau[0];
      tau_pre_(2*i + 1) = tau[1];

      real_t Vi[2];
      bp5_params_.V_init_vec(x2, x3, Vi);
      V_init_vec_(2*i)     = Vi[0];
      V_init_vec_(2*i + 1) = Vi[1];
   }
}
```

### New Accessors

| Method | Returns | Layout |
|--------|---------|--------|
| `GetDcValues()` | Critical slip distance per DOF | `[N]` scalar |
| `GetTauPre()` | Pre-stress vector per DOF | `[2N]`: `[τ_dip_0, τ_strike_0, ...]` |
| `GetVInit()` | Initial velocity per DOF | `[2N]`: `[V_dip_0, V_strike_0, ...]` |
| `GetCoordsX2()` | Along-strike coordinate | `[N]` scalar |
| `GetCoordsX3()` | Depth coordinate | `[N]` scalar |
| `IsBP5()` | Whether BP5 constructor was used | `bool` |
| `GetBP5Params()` | BP5 parameter reference | `const BP5Params&` |

### New Data Members

```cpp
// BP5-specific
BP5Params bp5_params_;
bool is_bp5_ = false;
Vector coords_x2_, coords_x3_;
Vector dc_values_;
Vector tau_pre_;        // [2*N]: (tau_dip, tau_strike) per DOF
Vector V_init_vec_;     // [2*N]: (V_dip, V_strike) per DOF
```

### Backward Compatibility

The existing BP2 constructor and all BP2-specific methods (`GetParams()`, `GetVWDepth()`, `GetVSDOFs()`, etc.) are unchanged. The BP5 path is selected at construction time via the `is_bp5_` flag. All pre-existing tests pass without modification.

---

## Data Flow

```
ElasticityDomainOperator
  ├── GetNumFaultDOFs()       → num_fault_dofs_
  ├── GetFaultCoords2D(x2,x3) → coords_x2_, coords_x3_
  └── GetFaultDepths(depths)   → depths_

FaultGeometry (BP5 constructor)
  ├── Calls BP5Params::a_of_x2_x3(x2, x3)  → a_values_
  ├── Calls BP5Params::eta()                 → eta_values_ (constant)
  ├── Calls BP5Params::Dc_of_x2_x3(x2, x3)  → dc_values_
  ├── Calls BP5Params::tau0_vec(x2, x3)      → tau_pre_ [2N]
  └── Calls BP5Params::V_init_vec(x2, x3)   → V_init_vec_ [2N]
```

---

## Spatially Varying Parameters (BP5 vs BP2)

| Parameter | BP2 (1D) | BP5 (2D) |
|-----------|----------|----------|
| a | `a(z)` — depth only | `a(x2, x3)` — VW core, transition, VS zones |
| b | Constant (0.019) | Constant (0.03) |
| Dc (L) | Constant (0.008 m) | `L(x2, x3)` — 0.14 m default, 0.13 m nucleation |
| eta | Constant `μ/(2cs)` | Constant `μ/(2cs)` |
| sigma_n | Constant (50 MPa) | Constant (25 MPa) |
| tau_pre | Scalar `tau0()` | Vector `tau0_vec(x2, x3)` [2-component] |
| V_init | Constant `V_init` | Vector `V_init_vec(x2, x3)` [2-component] |

---

## Parallel Considerations

For parallel execution with `ParMesh`:
- Each rank computes parameters for its local fault DOFs only
- `GatherToRoot` and `GatherToRootDedup` generalize from 1D depth to the existing depth-based approach (using `coords_x3_` as depth)
- Future: Full 2D `(x2, x3)` deduplication may be needed for non-Cartesian partitioning

---

## Tests

Tested via `tests/unit/test_elasticity_operator.cpp`, `TestFaultGeometry3D()`:

1. DOF count matches domain operator
2. `IsBP5()` returns true
3. `GetAValues()` has correct size `[N]`
4. `GetEtaValues()` has correct size `[N]`, values match `params.eta()`
5. `GetDcValues()` has correct size `[N]`
6. `GetTauPre()` has correct size `[2N]`

---

## Files Modified

| File | Change |
|------|--------|
| `fault/fault_geometry.hpp` | Added BP5 constructor, ComputeBP5Params(), new accessors, BP5 data members |

## Key Design Notes

- **No new files**: Only extends existing `fault_geometry.hpp`
- **Constructor overloading**: BP2 vs BP5 distinguished by parameter type (`BP2Params` vs `BP5Params`)
- **Single-pass parameter computation**: All per-DOF parameters computed once at construction, cached for time stepping
- **2-component vectors**: `tau_pre_` and `V_init_vec_` use interleaved layout `[dip_0, strike_0, dip_1, strike_1, ...]` consistent with `RateStateFaultOperator` state layout
