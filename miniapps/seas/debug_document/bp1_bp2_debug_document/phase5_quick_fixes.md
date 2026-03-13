# Phase 5 Quick Fixes

## Fix 1: Add Assertions for Negative Physical Values — COMPLETE

**File:** `io/benchmark_output.hpp`, lines 129-132

**Applied:** `MFEM_ASSERT` guards added before clamping in `Write()` method.

---

## Fix 2: Remove Unused `params_` Member — COMPLETE

**File:** `io/benchmark_output.hpp`

**Applied:** `params_` member and `params_(params)` initializer removed. `BP2Params` kept in constructor signature for API compatibility.

---

## Fix 3: Name the Output Tolerance Constant — COMPLETE

**File:** `io/benchmark_output.hpp` line 51, `io/paraview_output.hpp` line 44

**Applied:** `static constexpr real_t kOutputTimeTolerance = 0.99` added to both classes with documenting comment.

---

## Fix 4: Add Probe Location Header — COMPLETE

**File:** `io/benchmark_output.hpp` lines 80-84, `io/probe_output.hpp` lines 43-55

**Applied:** `ProbeOutput` constructor accepts optional `description` parameter. `BenchmarkOutput` passes `"BP2-QD time series at z = X km"` for each probe.

---

## Fix 5: Clarify Interpolator Weight Default — COMPLETE

**File:** `io/probe_output.hpp` line 244

**Applied:** Default changed from `0.5` to `0.0`.
