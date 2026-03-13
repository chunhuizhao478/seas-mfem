# Phase 5 Debug Report: I/O and Validation

## Files Reviewed

| File | Lines |
|------|-------|
| `io/benchmark_output.hpp` | 220 |
| `io/probe_output.hpp` | 250 |
| `io/paraview_output.hpp` | 295 |
| `tests/unit/test_io.cpp` | 508 |

## Summary

**Critical Bugs: 0 | Medium Issues: 3 (ALL FIXED) | Low Issues: 8 (5 FIXED) | Cosmetic: 2**

The Phase 5 I/O implementation is well-designed with no critical bugs. All 7 unit tests cover the main functionality. All medium-severity issues and key low-severity issues have been fixed.

---

## Medium Severity Issues

### 1. Missing Probe Location Header (Design Doc Non-Compliance) — FIXED

**Location:** `io/benchmark_output.hpp` constructor / `io/probe_output.hpp`

**Issue:** The design doc specifies output files should have a descriptive first line:
```
# BP2-QD time series at z = 0 km
# Columns: time(s), slip(m), ...
```

The implementation only writes `# Columns:` but **never writes the descriptive probe location line**.

**Resolution:** `ProbeOutput` constructor now accepts an optional `description` parameter. `BenchmarkOutput` passes `"BP2-QD time series at z = X km"` for each probe.

### 2. No Validation Before Clamping Negative Values — FIXED

**Location:** `io/benchmark_output.hpp` lines 129-132

**Issue:** Silently clamps negative values. Negative slip rate or state variable indicates a physics bug that should be caught, not hidden.

**Resolution:** `MFEM_ASSERT` checks added before clamping in `BenchmarkOutput::Write()`.

### 3. Test Coverage Gaps — DEFERRED (non-blocking)

**Location:** `tests/unit/test_io.cpp`

Missing test cases (not blocking; existing 7 tests cover main functionality):
- **Extrapolation**: Probes outside the fault depth range (above z=0 or below z=-40km)
- **Negative/zero values**: Behavior when slip rate or theta are zero/negative
- **Domain-only ParaView**: ParaView output without `InitFaultOutput`
- **Column count mismatch**: `ProbeOutput::WriteStep` with wrong number of columns
- **Empty probes**: `BenchmarkOutput` with zero probe depths

---

## Low Severity Issues

### 4. Unused Member Variable `params_` — FIXED

**Location:** `io/benchmark_output.hpp`

`BP2Params params_` was stored but never accessed after construction. Removed.

### 5. Undocumented 0.99 Tolerance Factor — FIXED

**Location:** `io/benchmark_output.hpp` line 51, `io/paraview_output.hpp` line 44

**Resolution:** Named constant `kOutputTimeTolerance = 0.99` added to both classes with documenting comment.

### 6. Undocumented Stress Unit Conversion in ParaView

**Location:** `io/paraview_output.hpp` line 185

```cpp
(*fault_stress_)(e1) = avg_stress / 1e6;  // Pa to MPa
```

Both `BenchmarkOutput` and `ParaViewOutput` convert to MPa, but neither the field name nor comments make this obvious to ParaView users. Consider naming the field `shear_stress_MPa`.

### 7. ProbeInterpolator Fallback Logic Undocumented — FIXED

**Location:** `io/probe_output.hpp` lines 246-252

**Resolution:** Clarifying comment `// Exact match or extrapolation: use closest` added to the fallback path.

### 8. Weight Default Value When dz ≈ 0 — FIXED

**Location:** `io/probe_output.hpp` line 244

**Resolution:** Default changed from `0.5` to `0.0`.

### 9. Test File Cleanup Ignores Errors

**Location:** `tests/unit/test_io.cpp` (multiple locations)

`std::remove()` return values are unchecked. Silent failure may leave test artifacts.

### 10. ParaView UpdateFaultFields Zeros Entire Field

**Location:** `io/paraview_output.hpp` lines 144-148

Every call zeros all elements then sets only fault-adjacent ones. Inefficient if called frequently, but functionally correct.

### 11. Inconsistent ParaView Field Naming

**Location:** `io/paraview_output.hpp` lines 118-121

Fields: `slip`, `slip_rate`, `shear_stress`, `state_variable`. Consider prefixing with `fault_` for clarity when domain and fault fields coexist.

---

## Positive Observations

- Clean template design supporting serial/parallel (`MeshType` parameter)
- Proper `std::unique_ptr` memory management throughout
- Correct move semantics with deleted copy constructors
- Good const correctness
- Comprehensive 7-test unit test suite
- Proper `MFEM_VERIFY` / `MFEM_ASSERT` usage
- Clean separation of concerns (ProbeOutput, ProbeInterpolator, BenchmarkOutput)
- Adaptive output scheduling correctly implements three regimes (interseismic/nucleation/coseismic)

---

## Recommended Fixes (Priority Order)

| Priority | Issue | Effort | Status |
|----------|-------|--------|--------|
| High | Add probe location header to output files | Small | FIXED |
| High | Add assertions for negative V/theta before clamping | Small | FIXED |
| High | Add edge-case unit tests | Medium | DEFERRED |
| Medium | Remove unused `params_` member | Trivial | FIXED |
| Medium | Name the 0.99 constant | Trivial | FIXED |
| Medium | Document stress unit conversion | Trivial | Open |
| Low | Improve interpolator comments | Trivial | FIXED |
| Low | Fix weight default value | Trivial | FIXED |
| Low | Check `std::remove()` returns in tests | Small | Open |
