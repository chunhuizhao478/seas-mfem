# Fix BP5 Newton Solver Failure in DieterichRuinaFriction

## Context

The BP5 simulation aborts at step ~82 (t≈20s) with Newton solver failure in `SolveSlipRatePsi()`. The failing parameters are:
- `tau=1.42e7, psi=0.789, sigma_n=2.5e7, eta=4.62e6, a=0.004`
- `psi/a = 197.25` → `exp(psi/a) ≈ 7e85`

**Root cause**: With large `psi/a`, even at `V_min = 1e-30`, the friction coefficient is so large that `sigma_n * f(V_min) > tau`. The true solution requires `V ≈ 1e-86`, but the Newton solver's lower bound `V_lo = 1e-30` prevents convergence. The solver gets stuck at V_lo indefinitely.

### Detailed Numerical Analysis

At the failing node (rank 25), the residual function `F(V) = sigma_n * f(V, psi) + eta * V - tau` evaluates as:

| V value | sinh_arg = (V/2V0)*exp(psi/a) | f = a*asinh(sinh_arg) | sigma_n*f | F(V) |
|---------|-------------------------------|----------------------|-----------|------|
| V0 = 1e-6 | 0.5 * 7e85 ≈ 3.5e85 | ≈ 0.79 | 19.8e6 | +5.6e6 (too high) |
| V_min = 1e-30 | 5e-25 * 7e85 ≈ 3.5e61 | ≈ 0.57 | 14.25e6 | +0.06e6 (still positive!) |
| V ≈ 1e-86 | ≈ 1.0 | ≈ 0.004*asinh(1) ≈ 0.0035 | ≈ 0.088e6 | ≈ -14.1e6 (sign change) |
| V = 0 | 0 | 0 | 0 | -14.2e6 (negative) |

The true solution lies between 1e-86 and 1e-30, far below the current V_min floor.

## Current Implementation (Broken)

**File**: `miniapps/seas/friction/dieterich_ruina.hpp`, lines 337-383

The current `SolveSlipRatePsi()` uses Newton-Raphson with:
- `V_lo = V_min_ = 1e-30` (too high — prevents reaching the true solution)
- `V_hi = min(tau/eta, 100.0)`
- Initial guess: `V = V0 = 1e-6`
- When Newton step would go below V_lo, it bisects: `V_new = 0.5*(V + V_lo)`
- This traps the solver at V_lo where F > 0 for all iterations → failure

## Tandem's Approach (Reference)

**File**: `/Users/chunhuizhao/projects/tandem/app/localoperator/DieterichRuinaAgeing.h`, lines 81-120

Tandem uses **Brent's method** (not Newton):
- Bracket: `[0, tau/eta]`
- Residual: `R(V) = tauAbs - sigma_n * f(V, psi) - eta * V`
- At V=0: R(0) = tauAbs > 0 (always positive)
- At V=tau/eta: R = -sigma_n * f(...) < 0 (always negative)
- Guaranteed sign change → guaranteed convergence

**Brent's method implementation**: `/Users/chunhuizhao/projects/tandem/src/util/Zero.cpp`
- Hybrid: inverse quadratic interpolation + bisection fallback
- NaN/Inf detection at each step
- Convergence: `|xm| <= 2*eps*|b| + tol/2`

Key design difference: Tandem **never clamps V to a minimum** — it lets V reach exactly 0.0, which is valid because `0.0 * exp(psi/a) = 0.0` in IEEE 754 (as long as exp(psi/a) is finite, i.e., psi/a < 710).

## Proposed Fix

### Replace Newton-Raphson with Brent's method in `SolveSlipRatePsi()`

**File to modify**: `miniapps/seas/friction/dieterich_ruina.hpp`

### Changes:

1. **Add a private `zeroIn()` method** — port Tandem's Brent solver from `Zero.cpp`

2. **Rewrite `SolveSlipRatePsi()`**:
   ```cpp
   real_t SolveSlipRatePsi(real_t tau, real_t psi, real_t sigma_n,
                           real_t eta, real_t a, int *iterations = nullptr) const
   {
      if (sigma_n <= 0.0) { /* existing tension handling */ }

      // Bracket: [0, tau/eta]
      real_t V_lo = 0.0;
      real_t V_hi = tau / eta;  // When friction=0, all stress goes to damping

      // Residual: R(V) = tau - sigma_n * f(V, psi) - eta * V
      // R(0) = tau > 0, R(tau/eta) = -sigma_n * f < 0  → valid bracket
      auto residual = [&](real_t V) -> real_t {
         real_t f = FrictionCoefficientPsi(V, psi, a);
         return tau - sigma_n * f - eta * V;
      };

      real_t V = zeroIn(V_lo, V_hi, residual, 1e-12);
      return V;
   }
   ```

3. **Fix `FrictionCoefficientPsi()` for V=0**: Remove `V = std::max(V, V_min_)` clamping, or add special handling:
   ```cpp
   real_t FrictionCoefficientPsi(real_t V, real_t psi, real_t a) const
   {
      if (V == 0.0) return 0.0;  // f(0, psi) = a * asinh(0) = 0
      real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(psi / a);
      return a * std::asinh(sinh_arg);
   }
   ```

4. **Also update `SolveSlipRate()` (theta-space)** with the same Brent approach for consistency.

### Edge Case: exp(psi/a) Overflow (psi/a > 710)

If psi/a exceeds ~710, `exp(psi/a)` overflows to inf, and `0.0 * inf = NaN`. Add protection:
```cpp
if (V == 0.0 || (V < 1e-300 && psi/a > 700.0)) return 0.0;
```

For non-zero V with large psi/a, use log-space evaluation:
```
f = a * asinh(exp(log(V/(2V0)) + psi/a))
```
For large argument x: `asinh(x) ≈ log(2x) = log(2) + log(x)`, so:
```
f ≈ a * (log(2) + log(V/(2V0)) + psi/a) = a*log(2) + a*log(V/(2V0)) + psi
```

## Verification

1. Build: `conda activate mfem-dev && make` in the seas build directory
2. Run BP5 — should pass step 82 without Newton failure
3. Check that slip rate V at the previously-failing node is physically reasonable (very small, ~1e-86, consistent with a locked fault segment at depth with high normal stress)
4. Verify overall simulation behavior matches expected BP5 benchmark results
