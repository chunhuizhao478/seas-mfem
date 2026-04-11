# BP5 Frontera Fix: zeroIn Bracket Failure at Scale

## Root Cause Analysis

The Frontera run (100 ranks, 62992-element Gmsh mesh) aborts at the **first RK45 time step** (dt=2.17s) with:
```
zeroIn: F(a) and F(b) must have different signs.
a=0 F(a)=5.15e+10 b=1.11e+04 F(b)=7.63e-06
```

### Chain of Events

1. **Initial dt too large**: dt_init = 0.5 * L_nuc / V_nuc = 2.17s
2. **RK45 stages amplify**: DormandPrince has 6 stages; stages 1-3 produce increasing V
3. **Nucleation instability**: V_nuc=0.03 m/s in VW zone → stress drop from slip → |tau| increases → V increases (earthquake nucleation feedback)
4. **State variable collapse**: High V causes rapid psi decrease (dpsi/dt ∝ -V/V0, V/V0 = 1e7 at V=10 m/s)
5. **Friction vanishes**: When psi << 0, exp(psi/a) → 0, so f(V,psi) → 0
6. **Bracket failure**: V_hi = tau/eta. With f≈0: F(V_hi) = -sigma_n*f ≈ -0 ≈ +7.6e-6 (floating point error makes it slightly positive instead of negative)
7. **MFEM_VERIFY aborts**: zeroIn fails, MPI_Abort kills all ranks before RK45 can reject the step

### Why It Works Locally but Fails at Scale

This is NOT a parallel-specific bug. The same dynamics would occur on any mesh with a nucleation zone at V_nuc=0.03 m/s when dt is too large. It passes locally because:
- The tiny 4×2×1 test mesh has only 2-4 fault DOFs
- The test only does a single Mult() + one Euler step (dt=0.1s), not full RK45

### Tandem Comparison

Tandem has the same bracket [0, tau/eta] and also throws on bracket failure. However, Tandem uses PETSc's time stepper which can catch exceptions and retry with smaller dt. Our MFEM_VERIFY calls MPI_Abort with no recovery.

## Fix Plan

### Fix 1: Make Friction Solver Robust (CRITICAL)

**File**: `friction/dieterich_ruina.hpp`

**Change in `SolveSlipRatePsi()`** (lines 312-336):
- Before calling `zeroIn`, check if F(V_hi) ≥ 0 (bracket invalid)
- If invalid: return V = tau/eta (the physically correct limit when friction → 0)
- This allows the RK45 error estimator to compute a large error and reject the step

```cpp
real_t SolveSlipRatePsi(real_t tau, real_t psi, real_t sigma_n,
                        real_t eta, real_t a,
                        int *iterations = nullptr) const
{
   if (sigma_n <= 0.0) { /* existing code */ }
   if (tau <= 0.0) { return 0.0; }

   real_t V_lo = 0.0;
   real_t V_hi = tau / eta;

   auto residual = [&](real_t V) -> real_t {
      return tau - sigma_n * FrictionCoefficientPsi(V, psi, a) - eta * V;
   };

   real_t Fa = residual(V_lo);  // = tau > 0
   real_t Fb = residual(V_hi);  // should be -sigma_n*f < 0

   // Handle degenerate bracket: when psi is very negative, f → 0,
   // so F(V_hi) ≈ 0 and may be slightly positive due to floating point.
   // In this regime, friction is negligible and V ≈ tau/eta.
   if (Fb >= 0.0)
   {
      if (iterations) { *iterations = 0; }
      return V_hi;  // V = tau/eta (no friction resistance)
   }

   real_t V = zeroIn(V_lo, V_hi, residual);
   if (iterations) { *iterations = 0; }
   return V;
}
```

### Fix 2: Reduce Initial Time Step (IMPORTANT)

**File**: `tests/verification/bp5_verification_full.cpp`

**Change**: Reduce the initial dt multiplier from 0.5 to 0.01:
```cpp
// Old: dt_init = min(1000, 0.5 * L_nuc / max(V_init, 1e-20))  → 2.17s
// New: dt_init = min(1000, 0.01 * L_nuc / max(V_init, 1e-20)) → 0.043s
real_t dt_init = std::min(1e3, 0.01 * params.L_nuc /
                          std::max(V_init, 1e-20));
```

This gives dt_init ≈ 0.043s, which is 50× smaller. The adaptive RK45 controller will quickly ramp up during the interseismic period. This matches the typical practice in SEAS codes of starting conservatively.

### Fix 3: Add Diagnostic Warning (NICE TO HAVE)

When the bracket is degenerate, log a warning with the relevant state:
```cpp
if (Fb >= 0.0)
{
   // Log warning (once per step, not per DOF)
   return V_hi;
}
```

## Files Changed

1. `friction/dieterich_ruina.hpp` — Robust bracket handling in `SolveSlipRatePsi()`
2. `tests/verification/bp5_verification_full.cpp` — Smaller initial dt

## Testing Plan

1. Rebuild locally: `seas_test_elasticity_operator`, `seas_test_bp5_parallel_smoke`, `seas_bp5_full`
2. Run all existing tests (should still pass — the fix only affects degenerate cases)
3. Submit to Frontera: the simulation should now survive the initial transient and begin normal time stepping
