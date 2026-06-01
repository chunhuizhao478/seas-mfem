# Fix Report: Lever-3 review R-001..R-005

## Summary
- Findings addressed: **5 of 5**.
- Files modified: `dynamic/wave_operator.inl`, `tests/unit/test_wave_operator_spatial_derivative.cpp`,
  `tests/unit/bench_ader_hotpath.cpp`, `Makefile`.
- Files added: `tests/parallel/test_wave_operator_cached_parallel.cpp`.
- Tests added: 6 serial assertions (heterogeneous ×2, Mult ×2, R-003 ×2) + 1
  parallel test (np=2/np=4).
- Test suite: **PASS** — serial `seas_test_wave_operator_spatial_derivative`
  **54/54**; parallel cached-equivalence np=2 **1.39e-15**, np=4 **1.40e-15**
  (both exit 0); regression green incl. `bimaterial_wave_operator_parity` 9/9.

## Changes Made
1. **R-001 (heterogeneous material)** → `test_…spatial_derivative.cpp`:
   `TestCachedEquivalence` now also builds a `BimaterialWaveOperator` from
   **position-varying `FunctionCoefficient`s** (per-element-varying `A_d^e`, as
   TPV31's `depth_profile_1d`) and asserts cached `ComputeVolumeRHS` == OnTheFly.
   Result: **3.1e-16 / 4.4e-16** (orders 2/3). This exercises the per-element
   flux application that the constant-material case could not.
2. **R-002 (Mult + ParMesh coverage)** →
   - **Mult (RK path):** added a serial `Mult` cached-vs-OnTheFly check
     (`5.5e-16 / 1.3e-15`), so the cached `ComputeVolumeRHS` reached via `Mult`
     is now covered (the ADER tests do not exercise `Mult`).
   - **ParMesh:** new standalone `tests/parallel/test_wave_operator_cached_parallel.cpp`
     (+ `seas_test_wave_operator_cached_parallel` target) — partitioned ParMesh,
     per-rank `ApplySpatialDerivative` + `ComputeVolumeRHS` cached==OnTheFly,
     `MPI_Allreduce(MAX)`. np=2 **1.39e-15**, np=4 **1.40e-15**.
3. **R-003 (combined 2× budget)** → `test_…spatial_derivative.cpp`: assert
   `SetDerivMode(Cached)` SUCCEEDS at a budget == `2×ElementDerivativeCacheBytes`
   (a too-low guard would `MFEM_ABORT` before the assertion). The abort path
   itself remains inspection-only (no death-test facility).
4. **R-004 (homogeneous-order guard)** → `wave_operator.inl`: added
   `MFEM_VERIFY(ndof == ndof_per_el_, …)` to the cached `ComputeVolumeRHS`
   branch AND the cached `ApplySpatialDerivative` branch (the deferred Lever-1
   R-004), for symmetry with the OnTheFly kernels. The homogeneous regression
   set still passes (guard never trips on a homogeneous mesh).
5. **R-005 (bench `rhs_vol`)** → `bench_ader_hotpath.cpp`: zero `rhs_vol` at the
   start of each `measure()` (outside the timed block) and add it to the final
   checksum.

## Deviation (documented, per the fix process)
- **R-002 ParMesh was put in a NEW file, not the existing
  `test_parallel_wave_operator.cpp`.** While adding the parallel check there, I
  found that file's P4/P5 are **pre-existing-broken at the worktree base
  (e0a99f0)**: they construct `WaveOperator<ParMesh>` with `bc.absorbing_attrs={5}`
  but never set `bc.fault_attr=0` (BoundaryConfig defaults to 3) NOR call
  `SetAbsorbingBackground`, so `Mult` aborts on `fault_attr=3` then on
  `has_bulk_bg_` — at both np=1 and np=2. These are **unrelated to Lever 3**
  (OnTheFly `Mult`, untouched by this work). Per the fix process I did NOT fix
  the pre-existing P4/P5 (noted here instead); I reverted my edits to that file
  and delivered R-002b as a standalone test that does not depend on the broken
  `Mult` path.

## Verification
- [x] R-001: heterogeneous (varying `A_d^e`) volume RHS cached==OnTheFly ≤4.4e-16.
- [x] R-002: Mult (serial) cached==OnTheFly ≤1.3e-15; ParMesh np=2/np=4 ≤1.4e-15.
- [x] R-003: `SetDerivMode(Cached)` succeeds at the combined 2× budget.
- [x] R-004: homogeneous guard in both cached branches; regression unaffected.
- [x] R-005: `rhs_vol` zeroed per measure + in checksum.
- [x] Gates intact: Lever-2 bit-exact parity still 0.000e+00; all ≤1e-12
      equivalence checks green; default `OnTheFly` byte-exact (parity 9/9).

## Unresolved Findings
- None of R-001..R-005. (Out of scope, noted not fixed: the pre-existing P4/P5
  setup bugs in `test_parallel_wave_operator.cpp`.)

## New Tests
- heterogeneous volume RHS, Mult-path, R-003 budget (serial) —
  `test_wave_operator_spatial_derivative.cpp` (now 54/54).
- `test_wave_operator_cached_parallel.cpp` — ParMesh cached==OnTheFly (R-002b).

## Ready for Re-Review: YES
