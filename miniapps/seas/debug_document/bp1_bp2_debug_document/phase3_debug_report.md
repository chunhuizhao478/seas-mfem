# Phase 3 Debug Report: Fault Operator

**Version 4.0** - Updated 2025-02-06 after thorough re-examination with Tandem/SCEC verification

## Overview

This report documents discrepancies between the Phase 3 documentation (`dg_antiplane_theory.md`, `phase3_fault_operator.md`) and the current implementation in the fault, friction, and config directories.

**Files Analyzed:**
- Documentation: `dg_antiplane_theory.md`, `phase3_fault_operator.md`
- Implementation: `fault/fault_geometry.hpp`, `fault/rate_state_fault.hpp`, `friction/friction_law.hpp`, `friction/dieterich_ruina.hpp`, `friction/state_evolution.hpp`, `config/bp2_params.hpp`

---

## Part IX: Unit Test Verification (NEW)

### Test Execution Results

**Friction Law Tests (`seas_test_friction_law`):**
- Total tests: 32
- Passed: 32
- Failed: 0

**State Evolution Tests (`seas_test_state_evolution`):**
- Total tests: 23
- Passed: 23
- Failed: 0

**TOTAL: 55/55 tests pass**

### Physical Reasonableness Verification

All test results are physically reasonable and match expected values:

| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| τ₀ (pre-stress) | ~26.546 MPa | 26.5461 MPa | ✅ CORRECT |
| θ₀ at z=0 | ~4000 s (log₁₀ ≈ 3.6) | 4000 s (log₁₀ = 3.60206) | ✅ CORRECT |
| Characteristic timescale | ~46 days | 46.3 days | ✅ CORRECT |
| Newton iterations (typical) | <50 | 0-24 | ✅ EFFICIENT |
| Velocity-weakening (a<b) | f_ss decreases with V | ✅ Verified | ✅ CORRECT |
| Velocity-strengthening (a>b) | f_ss increases with V | ✅ Verified | ✅ CORRECT |
| Aging law locked | dθ/dt ≈ 1 | 1.0 ± 1e-10 | ✅ CORRECT |
| Slip law locked | dθ/dt ≈ 0 | <1e-8 | ✅ CORRECT |
| State convergence | Approaches θ_ss | Within 10% | ✅ CORRECT |
| Analytical derivatives | Match numerical | Within 1e-5 | ✅ CORRECT |

---

## Critical Discrepancies

### 1. Root Finding Algorithm (MAJOR)

| Aspect | Documentation | Implementation | Impact |
|--------|--------------|----------------|--------|
| Algorithm | **Brent's method** (bisection + interpolation hybrid) | **Newton-Raphson** with bounds checking | Medium - both converge, Newton is faster |
| Reference | Tandem `src/util/Zero.cpp` lines 13-100 | `dieterich_ruina.hpp` lines 151-230 | Newton converges in 0-24 iterations |
| Bracketing | Uses interval [0, τ/η] | Uses interval [V_min, min(τ/η, 100)] | **MATCHES** Tandem's approach |

**Tandem Implementation (DieterichRuinaAgeing.h lines 97-98):**
```cpp
double a = 0.0;
double b = tauAbs / eta;  // Upper bound = τ/η
```

**Our Implementation (dieterich_ruina.hpp lines 173-179):**
```cpp
// Bounds for bracketing (slip rate should be positive)
real_t V_lo = V_min_;
// Physical upper bound: when friction = 0, tau = eta * V
// V_max = tau / eta (with fallback if eta is very small)
// Following Tandem's approach for tighter bracketing
real_t V_hi = (eta > 1e-6) ? (tau / eta) : 100.0;
V_hi = std::min(V_hi, 100.0);  // Cap at 100 m/s for safety
```

**Analysis:** The implementation now correctly uses τ/η as the upper bound, matching Tandem's approach. The key difference is the root-finding algorithm:
- **Tandem**: Uses Brent's method (guaranteed convergence, ~30-50 iterations)
- **MFEM SEAS**: Uses Newton-Raphson (faster, 0-24 iterations, requires derivative)

Both approaches are valid. Newton-Raphson is faster because we already compute df/dV for implicit time stepping. The bounds checking ensures robustness.

**Recommendation:** Document this deliberate design choice in phase3_fault_operator.md

---

### 2. State Variable Representation

| Aspect | Documentation | Implementation | Mathematically Equivalent? |
|--------|--------------|----------------|---------------------------|
| Variable | ψ (psi) natural log form | θ (theta) physical age | **Yes** |
| Friction formula | f = a · asinh[(V/2V₀) · exp(ψ/a)] | f = a · asinh[(V/2V₀) · exp((f₀ + b·ln(V₀θ/Dc))/a)] | Yes |
| State evolution | dψ/dt = (b·V₀/L) · [exp((f₀-ψ)/b) - V/V₀] | dθ/dt = 1 - V·θ/Dc | Yes |

**Relationship:** ψ = f₀ + b·ln(V₀·θ/Dc)

**Documentation (dg_antiplane_theory.md Section 6, phase3_fault_operator.md lines 69-107):**
Uses ψ (psi) form throughout, with the ageing law:
```
dψ/dt = (b·V₀/L) · [exp((f₀ - ψ)/b) - V/V₀]
```

**Implementation (state_evolution.hpp line 87-90):**
Uses θ (theta) form directly:
```cpp
real_t Rate(real_t V, real_t theta, real_t Dc) const override
{
   return 1.0 - V * theta / Dc;
}
```

**Analysis:** Both formulations are mathematically equivalent. The implementation uses the more intuitive "physical age" interpretation (θ has units of seconds), which is actually what Tandem uses internally despite the documentation using ψ notation. The implementation is **correct** but the documentation should be clearer about this equivalence.

**Recommendation:** Update documentation to:
1. Clarify the relationship between ψ and θ
2. Explicitly note that implementation uses θ form
3. Show both forms are equivalent

---

### 3. Missing Root Finder File

| Aspect | Documentation | Implementation | Status |
|--------|--------------|----------------|--------|
| File | `fault/root_finder.hpp` | Does not exist | **MISSING** |
| Location | Listed in phase3_fault_operator.md file organization | Not created | File structure inconsistent |

**Documentation (phase3_fault_operator.md lines 751-760):**
```
fault/
├── fault_geometry.hpp       # Fault surface geometry
├── friction_law.hpp         # Dieterich-Ruina Ageing law
├── root_finder.hpp          # Brent's method root finder  <-- MISSING
├── rate_state_fault.hpp     # Rate-state fault operator
...
```

**Current Implementation Structure:**
```
fault/
├── fault_geometry.hpp       # EXISTS
├── rate_state_fault.hpp     # EXISTS
friction/
├── friction_law.hpp         # EXISTS
├── dieterich_ruina.hpp      # EXISTS (contains Newton solver)
├── state_evolution.hpp      # EXISTS
```

**Analysis:** The documentation specifies a separate `root_finder.hpp` file with Brent's method implementation, but this file was never created. The slip rate solving is instead embedded in `dieterich_ruina.hpp` using Newton-Raphson.

**Recommendation:** Either:
1. Create `root_finder.hpp` with Brent's method as documented, OR
2. Update documentation to reflect current file organization

---

### 4. Critical Slip Distance Value (Dc vs L)

| Aspect | phase3_fault_operator.md | bp2_params.hpp | SCEC BP2 Spec |
|--------|-------------------------|----------------|---------------|
| Parameter name | L (Dc) = 0.008 m | Dc = 0.004 m | Dc = 4 mm = 0.004 m |

**Documentation (phase3_fault_operator.md line 45):**
```
| Characteristic length | L (Dc) | 0.008 | m |
```

**Implementation (bp2_params.hpp line 62):**
```cpp
/// Critical slip distance [m] (key BP2 parameter, differs from BP1)
real_t Dc = 0.004;
```

**SCEC BP2 Specification:**
```
Dc = 4 mm = 0.004 m
```

**Analysis:** The documentation lists L = 0.008 m, but the implementation correctly uses Dc = 0.004 m per the SCEC BP2 specification. The documentation value is **INCORRECT**.

**Verified from Tandem bp1.lua (line 38-40):**
```lua
function BP1:L(x, y)
    return 0.008  -- BP1 uses L = 0.008 m
end
```

**Note:** The 0.008 m value is from **BP1** specification. The documentation incorrectly used the BP1 value instead of the BP2 value (0.004 m). This is likely a copy-paste error from BP1 documentation.

**Recommendation:** Fix documentation to use Dc = 0.004 m (or 4 mm) to match both the implementation and SCEC BP2 specification.

---

### 5. Slip Rate Solver Upper Bound - **[FIXED]**

| Aspect | Documentation | Implementation | Status |
|--------|--------------|----------------|--------|
| Upper bound | V_max = τ/η (physically motivated) | V_hi = (eta > 1e-6) ? (tau/eta) : 100.0 | **FIXED** |
| Source | Tandem Zero.cpp | dieterich_ruina.hpp lines 175-179 | Matches Tandem |

**Current Implementation (dieterich_ruina.hpp lines 175-179):**
```cpp
// Physical upper bound: when friction = 0, tau = eta * V
// V_max = tau / eta (with fallback if eta is very small)
// Following Tandem's approach for tighter bracketing
real_t V_hi = (eta > 1e-6) ? (tau / eta) : 100.0;
V_hi = std::min(V_hi, 100.0);  // Cap at 100 m/s for safety
```

**Analysis:** The implementation now uses the physically-motivated upper bound τ/η as documented, with:
- Fallback to 100 m/s if η is very small (avoiding division by near-zero)
- Safety cap at 100 m/s (seismic upper limit)

This matches Tandem's approach and provides tighter bracketing for faster Newton convergence.

**Status:** ✅ FIXED - Implementation now matches documentation.

---

## Minor Discrepancies

### 6. Friction Law File Organization

| Documentation | Implementation |
|--------------|----------------|
| `fault/friction_law.hpp` | `friction/friction_law.hpp` |
| Single friction directory | Separate `friction/` directory |

The implementation uses a separate `friction/` directory, which is actually better organization than the documented `fault/` subdirectory.

**Recommendation:** Update documentation to reflect the `friction/` directory structure.

---

### 7. State Layout Comments

**Documentation (phase3_fault_operator.md line 499):**
```
/// State layout: [slip_0, psi_0, slip_1, psi_1, ...]
```

**Implementation (rate_state_fault.hpp lines 38-42):**
```cpp
/// State layout: [slip_0, theta_0, slip_1, theta_1, ...]
/// - Each fault node has 2 state variables: slip and theta (state variable)
/// - slip: accumulated fault slip [m]
/// - theta: state variable from rate-and-state friction law [s]
```

**Analysis:** Documentation uses `psi` but implementation uses `theta`. Both refer to the state variable, but the naming is inconsistent. Implementation is clearer with units.

**Recommendation:** Harmonize naming in documentation to use `theta` with clear units.

---

### 8. DieterichRuinaFriction Class Name

| Documentation | Implementation |
|--------------|----------------|
| `DieterichRuinaAgeing` | `DieterichRuinaFriction` |

**Documentation (phase3_fault_operator.md line 406):**
```cpp
class DieterichRuinaAgeing {
```

**Implementation (dieterich_ruina.hpp line 38):**
```cpp
class DieterichRuinaFriction : public FrictionLaw
```

**Analysis:** The implementation name is more generic and better since the class handles the friction law, not the evolution law. The evolution law is correctly separated in `state_evolution.hpp`.

**Recommendation:** Update documentation to use `DieterichRuinaFriction`.

---

## Consistent Elements (No Discrepancies)

The following elements are correctly implemented per documentation:

### Pre-stress τ₀ Computation
- Documentation: ~26.546 MPa
- Implementation: `bp2_params.hpp::tau0()` correctly computes this value

### Radiation Damping
- Documentation: η = √(μ·ρ) / 2
- Implementation: η = μ / (2·cs) in `bp2_params.hpp::eta()`
- These are equivalent: μ = ρ·cs² → μ/(2cs) = ρ·cs/2 = √(ρμ)/2 ✓

### Depth-Dependent a(z) Profile
- Correctly implemented in `bp2_params.hpp::a_of_z()`
- Matches Tandem bp1.lua specification

### BP2 Parameters
All material properties match SCEC specification:
- ρ = 2670 kg/m³ ✓
- cs = 3464 m/s ✓
- σ_n = 50 MPa ✓
- a₀ = 0.010 ✓
- amax = 0.025 ✓
- b = 0.015 ✓
- V₀ = 10⁻⁶ m/s ✓
- f₀ = 0.6 ✓

### State Evolution Laws
- AgingLaw: dθ/dt = 1 - V·θ/Dc correctly implemented
- SlipLaw: dθ/dt = -V·θ/Dc · ln(V·θ/Dc) correctly implemented

---

## Summary Table

| Issue | Severity | Status | Fix Required |
|-------|----------|--------|--------------|
| Root finding algorithm (Brent vs Newton) | Medium | Discrepancy | Optional - document choice |
| State variable (ψ vs θ notation) | Low | Equivalent | Clarify in docs |
| Missing root_finder.hpp | Low | Missing file | Update docs or create file |
| Dc value (0.008 vs 0.004) | **High** | Doc error | Fix documentation |
| Slip rate upper bound | Low | **FIXED** | ✅ No action needed |
| File organization | Low | Different | Update docs |
| State layout naming | Low | Inconsistent | Update docs |
| Class name | Low | Different | Update docs |

---

## Recommendations Priority

1. **HIGH:** Fix Dc value in documentation from 0.008 to 0.004 m
2. **MEDIUM:** Decide on root finding algorithm and update docs or code
3. **LOW:** Update documentation for consistent naming (θ vs ψ, class names, file organization)

---

## Part X: Implementation Verification Summary

### What's Working Correctly

1. **Friction Coefficient Computation**
   - Dieterich-Ruina regularized form correctly implemented
   - Direct effect (∂f/∂V > 0) verified
   - Evolution effect (∂f/∂θ > 0) verified
   - Velocity-weakening/strengthening behavior correct

2. **Slip Rate Solver (Newton-Raphson)**
   - Converges efficiently (0-24 iterations)
   - Handles edge cases: zero eta, negative sigma_n
   - Uses physically-motivated bounds [V_min, min(τ/η, 100)]

3. **Initial State Computation**
   - Correctly inverts friction law
   - Satisfies stress balance: τ₀ = σₙ·f(V_init, θ₀) + η·V_init
   - Handles small 'a' values without overflow

4. **State Evolution Laws**
   - Aging law: dθ/dt = 1 - V·θ/Dc ✓
   - Slip law: dθ/dt = -(V·θ/Dc)·ln(V·θ/Dc) ✓
   - Both converge to steady state θ_ss = Dc/V
   - Key difference preserved: Aging heals when locked, Slip does not

5. **BP2 Parameters**
   - τ₀ ≈ 26.546 MPa (matches SCEC spec)
   - η = μ/(2·cs) correctly computed
   - a(z) depth profile correct

### Items Still Requiring Documentation Updates

| Item | Action Required |
|------|-----------------|
| Dc value | Change 0.008 → 0.004 in phase3_fault_operator.md |
| ψ vs θ notation | Add clarification section to dg_antiplane_theory.md |
| Newton vs Brent | Document Newton choice in phase3_fault_operator.md |
| File organization | Update docs to reflect friction/ directory |

---

## Part XI: Tandem Reference Implementation Comparison (NEW)

### Verified Against Tandem Source Files

| File | Tandem | MFEM SEAS | Match |
|------|--------|-----------|-------|
| Friction law | `DieterichRuinaAgeing.h` | `dieterich_ruina.hpp` | ✅ Yes |
| Root finder | `util/Zero.cpp` | (inline in dieterich_ruina.hpp) | ⚠️ Different algorithm |
| State evolution | `DieterichRuinaAgeing.h:122-125` | `state_evolution.hpp` | ✅ Equivalent |

### Key Tandem Implementation Details Verified

**1. Friction Function F(σ_n, V, ψ)** (DieterichRuinaAgeing.h:157-161):
```cpp
double F(...) const {
    auto a = p_[index].get<A>();
    double e = exp(psi / a);
    double f = a * asinh((V / (2.0 * cp_.V0)) * e);
    return snAbs * f;
}
```
✅ Our implementation uses equivalent θ form

**2. Tension Handling** (DieterichRuinaAgeing.h:92-95):
```cpp
if (snAbs <= 0.0) {
    snAbs = 0.0;
    V = tauAbs / eta;  // Viscous sliding
}
```
✅ Our implementation matches (dieterich_ruina.hpp:157-167)

**3. Upper Bound for Slip Rate** (DieterichRuinaAgeing.h:97-98):
```cpp
double b = tauAbs / eta;  // V_max = τ/η
```
✅ Our implementation matches (dieterich_ruina.hpp:178)

**4. State Evolution (ψ form)** (DieterichRuinaAgeing.h:122-125):
```cpp
double state_rhs(...) const {
    double myL = p_[index].get<L>();
    return cp_.b * cp_.V0 / myL * (exp((cp_.f0 - psi) / cp_.b) - V / cp_.V0);
}
```
✅ Equivalent to our θ form: dθ/dt = 1 - V·θ/Dc

### Summary of Tandem Comparison

| Aspect | Tandem | MFEM SEAS | Status |
|--------|--------|-----------|--------|
| Friction formula | f = a·asinh[(V/2V₀)·exp(ψ/a)] | f = a·asinh[(V/2V₀)·exp((f₀+b·ln(V₀θ/Dc))/a)] | ✅ Equivalent |
| State variable | ψ (natural log form) | θ (physical age) | ✅ Equivalent |
| Root finder | Brent's method | Newton-Raphson | ⚠️ Different (faster) |
| Upper bound | τ/η | min(τ/η, 100) | ✅ Matches |
| Tension case | V = τ/η | V = τ/η | ✅ Matches |
| State evolution | dψ/dt = (bV₀/L)[exp((f₀-ψ)/b) - V/V₀] | dθ/dt = 1 - Vθ/Dc | ✅ Equivalent |

---

## Part XII: Latest Verification Summary (2025-02-06)

### Test Results Confirmation

All 55 Phase 3 unit tests pass:

```
=== Friction Law Tests ===
Total tests: 32 | Passed: 32 | Failed: 0

=== State Evolution Tests ===
Total tests: 23 | Passed: 23 | Failed: 0
```

### Physical Values Verified Against SCEC BP2 Specification

| Parameter | SCEC BP2 Spec | Implementation | Status |
|-----------|---------------|----------------|--------|
| ρ (density) | 2670 kg/m³ | 2670 kg/m³ | ✅ Match |
| cs (shear wave speed) | 3464 m/s | 3464 m/s | ✅ Match |
| σ_n (normal stress) | 50 MPa | 50 MPa | ✅ Match |
| Dc (critical slip distance) | 4 mm = 0.004 m | 0.004 m | ✅ Match |
| V₀ (reference velocity) | 10⁻⁶ m/s | 10⁻⁶ m/s | ✅ Match |
| f₀ (reference friction) | 0.6 | 0.6 | ✅ Match |
| a₀ (VW zone) | 0.010 | 0.010 | ✅ Match |
| amax (VS zone) | 0.025 | 0.025 | ✅ Match |
| b (state evolution) | 0.015 | 0.015 | ✅ Match |
| H (VW depth) | 15 km | 15000 m | ✅ Match |
| h (transition width) | 3 km | 3000 m | ✅ Match |
| V_p (plate rate) | 10⁻⁹ m/s | 10⁻⁹ m/s | ✅ Match |

### Computed Values Verified

| Computed Value | Expected | Actual | Status |
|----------------|----------|--------|--------|
| μ (shear modulus) | ~32 GPa | 32.04 GPa | ✅ Correct |
| η (radiation damping) | ~4.63 MPa·s/m | 4.625 MPa·s/m | ✅ Correct |
| τ₀ (pre-stress) | ~26.5 MPa | 26.5461 MPa | ✅ Correct |
| θ₀ at z=0 (initial state) | ~4000 s | 4000 s | ✅ Correct |
| θ_ss (steady state) | Dc/V = 4×10⁶ s | Converges correctly | ✅ Correct |

### Tandem Compatibility Verified

| Feature | Tandem | MFEM SEAS | Compatibility |
|---------|--------|-----------|---------------|
| Friction coefficient | f = a·asinh(arg) | Same formula | ✅ Identical |
| Tension handling | V = τ/η | V = τ/η | ✅ Identical |
| Upper bound | τ/η | min(τ/η, 100) | ✅ Compatible |
| State evolution (θ form equivalent) | dψ/dt formula | dθ/dt = 1 - Vθ/Dc | ✅ Equivalent |
| Initial state (psi_init equivalent) | ψ = a·ln(...) | θ = (Dc/V₀)·exp(...) | ✅ Equivalent |

### No New Issues Found

This re-examination confirms:
1. All unit tests pass with physically reasonable results
2. Implementation matches SCEC BP2 specification
3. Implementation is functionally equivalent to Tandem
4. All previously identified issues remain accurately documented
5. No regression from previous verification
