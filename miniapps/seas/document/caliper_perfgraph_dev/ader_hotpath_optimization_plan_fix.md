# Fix Report: ADER hot-path optimization — Step 0 (fold R-001…R-004) + Lever 1

## Summary
- Findings addressed: R-002, R-004 implemented in code; R-001, R-003 folded into
  the plan (their code lands with Lever 2, not in this Lever-1 step).
- Files modified: `dynamic/wave_operator.hpp`, `dynamic/wave_operator.inl`,
  `tests/unit/test_wave_operator_spatial_derivative.cpp`, `Makefile`,
  `document/caliper_perfgraph_dev/ader_hotpath_optimization_plan_2026-05-31.md`
  (+ regenerated `.pdf`), `tests/unit/bench_ader_hotpath.cpp`.
- Files added: `dynamic/elem_derivative_cache.hpp`.
- Tests added: 6 assertions (scalar + bimaterial equivalence × orders {2,3},
  + R-004 byte-count × {2,3}) in the spatial-derivative test.
- Test suite: **PASS.** Equivalence to ≤ **2.5e-15** relative (budget 1e-12);
  regression set all green incl. bimaterial parity 9/9 (default path byte-exact).
- Strategy: flagged in-place + isolated helper, `DerivMode` default `OnTheFly`
  (per `[[ader-opt-landing-strategy]]`).

## Step 0 — R-001…R-004 folded into the plan
Added a "Review guardrails" block at the top of §4 + targeted edits: Lever 1
gained the `DerivMode` landing note and the **R-004** memory-budget guard; Lever 2
gained the **R-001** (unconditional / no fault-less skip → deadlock) and **R-003**
(preserve both factorial-weight phases) warnings; §5 split Lever 1 acceptance into
unit-equivalence (≤1e-12) + end-to-end **physical** metrics, and gated Lever 2
parity against a **Lever-1-ON** baseline (**R-002**). PDF regenerated (7 pp.).

## Changes Made (Lever 1)
1. **New `dynamic/elem_derivative_cache.hpp`** — `BuildElementDerivativeOperators`
   builds per-element `D_d^e = M_e⁻¹ K_d^e` using the SAME `IntRules.Get(geom,
   2*order_)` + `CalcShape`/`CalcPhysDShape` as the on-the-fly kernel, then
   `Mult(Minv, Kd, D)`. Plus `ElementDerivativeCacheBytes(ne, ndof)` for the
   R-004 guard. Geometry-only ⇒ one cache serves scalar + bimaterial.
2. **`wave_operator.hpp`** — `enum class DerivMode { OnTheFly, Cached }`;
   `SetDerivMode` / `GetDerivMode` / `SetDerivCacheBudgetBytes`; private members
   `deriv_mode_` (default `OnTheFly`), `elem_deriv_op_`, `deriv_cache_budget_bytes_`
   (1 GiB default); file-scope include of the helper (the `.inl` is included
   inside `namespace mfem::seas`, so it cannot `#include`).
3. **`wave_operator.inl`** — `ApplySpatialDerivative` gains a leading
   `if (deriv_mode_ == Cached)` branch (per-element per-component `dQ_c = D·Q_c`
   mat-vec) that early-returns; the original quadrature kernel is the untouched
   `else`. `SetDerivMode` definition builds/frees the cache and enforces the
   **R-004** budget (`MFEM_VERIFY` with the MB figure) before allocating.
4. **`Makefile`** — `seas_test_wave_operator_spatial_derivative` link line gains
   the bimaterial deps (`GODUNOV_FLUX_POOL`/`HETEROGENEOUS_MATERIAL`/
   `SPATIAL_FRICTION`) for the BimaterialWaveOperator pass.
5. **`bench_ader_hotpath.cpp`** — upgraded to time BOTH modes and print the
   OnTheFly/Cached speedup (this is the Lever-1 measurement).

**Default-build invariance:** with the default `OnTheFly` the new branch is never
taken and `elem_deriv_op_` stays empty — the emitted hot path is the original
code. Confirmed by the byte-exact regression set below.

## Verification
- [x] **R-002 equivalence:** `Cached` vs `OnTheFly` ApplySpatialDerivative — scalar
      1.1e-15 / 2.5e-15 (order 2/3), bimaterial 1.5e-15 / 2.1e-15. ≤ 1e-12. ✓
- [x] **R-004 guard:** `ElementDerivativeCacheBytes` == 3·ne·ndof²·8 exactly
      (test), and `SetDerivMode(Cached)` aborts fail-loud above the budget (code;
      death-path not unit-tested — harness has no death tests). ✓
- [x] **Default byte-exact:** `seas_test_wave_operator` 23/23,
      `seas_test_bimaterial_wave_operator_parity` 9/9,
      `seas_test_ader_ck_predictor` 11/11,
      `seas_test_ader_linear_wave_equivalence` 4/4,
      `seas_test_wave_operator_spatial_derivative` 34/34. ✓
- [x] **R-001 / R-003:** folded into the plan; these gate the Lever-2 code (the
      merged recursion), not Lever 1 — Lever 1 does not touch the CK recursion.
- [x] **Measured speedup (Lever 1 alone, `make bench`):**

  | region | p2/O3 | p3/O3 |
  |---|---|---|
  | ApplySpatialDerivative (×3 dir) | **13.8×** | **18.3×** |
  | ComputeADERSubStepStates | 12.5× | 17.4× |
  | AdvanceADER (rec#2 + corrector) | 2.78× | 3.32× |
  | **macro-step aggregate** | **4.07×** | **5.10×** |

  Exceeds the plan's 3–8×-on-the-kernel estimate. (Build the worktree per
  `[[worktree-build-mfem-dir-override]]`; `make clean` per
  `[[makefile-no-header-deps-stale-o]]`.)

## Unresolved Findings
- **R-001, R-003** — not code-fixed here by design: they constrain the **Lever 2**
  merge (the redundant-recursion removal), which was explicitly the *next* step.
  Folded into the plan so the Lever-2 implementation honors them.
- **R-005** (p2/p3 doc clarity) — the benchmark now reports `ndof/elem`; confirm
  the profiled run's `[mesh].order` when quoting one memory number.
- **R-006** (hoist the per-step Allreduce) — Lever 4.1, not in Lever-1 scope.

## New Tests
- `TestCachedEquivalence` (in `test_wave_operator_spatial_derivative.cpp`) —
  covers R-002 (scalar + bimaterial, orders 2/3) and R-004 (byte count).

## Ready for Re-Review: YES (Lever 1)
Default builds are byte-identical; the `Cached` path is validated to round-off
and measured at 4–5× macro-step. Next: Lever 2 (merged CK recursion) under
R-001/R-003.
