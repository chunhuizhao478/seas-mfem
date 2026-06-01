# Implementation Complete: Lever 3 (cached volume-RHS operator)

## What was implemented
- **`BuildElementVolumeOperators`** (`dynamic/elem_derivative_cache.hpp`) — builds
  per-element geometry-only operators `S_d^e[i,m] = Σ_q w_q ∂_dφ_i(x_q) φ_m(x_q)`
  (the weak-derivative stiffness; `K_d^e` transposed), using the SAME
  `IntRules.Get(geom, 2*order_)` + `CalcShape`/`CalcPhysDShape` as the on-the-fly
  `ComputeVolumeRHS`.
- **`elem_volume_op_`** cache member (`wave_operator.hpp`), built/freed alongside
  `elem_deriv_op_` by `SetDerivMode`; empty under the default `OnTheFly`.
- **`Cached` branch in `ComputeVolumeRHS`** (`wave_operator.inl`): replaces the
  per-QP quadrature with `rhs_c += Σ_d A_d^e · (S_d^e Q)` — the geometry `S_d^e`
  is cached, the per-element flux matrices `A_d^e` (`FluxForElem_(e)`) are applied
  at run time, so one cache serves the scalar and bimaterial operators.
  Accumulates into `rhs` (`+=`), same contract as the OnTheFly path.
- **`SetDerivMode`** now builds both caches and the **R-004 budget** accounts for
  the combined `2 × 3·ne·ndof²·8 B` (D_d^e + S_d^e).
- **Equivalence test** + **bench** timing.

## Decisions made
- **Gated under the SAME `Cached` flag** (`--deriv-cache`): the flag now enables
  Lever 1 (`ApplySpatialDerivative`) AND Lever 3 (`ComputeVolumeRHS`) — both are
  "precomputed geometry operators". No new flag. Lever 2 (`--shared-ck-recursion`)
  is independent.
- **Separate build function** (not reusing Lever 1's `K_d^e` via transpose):
  Lever 1 stores `D_d^e = M⁻¹K_d^e`, not `K_d^e`, so there is nothing to transpose;
  a dedicated one-time builder is lower-risk and keeps Lever 1 untouched. The
  extra setup pass is negligible (built once, not per step).
- **Machine-eps, not bit-exact** (REVIEW R-002): the cached `Σ_q` folds into
  `S_d^e`, re-associating round-off vs the OnTheFly per-QP accumulation.
- **PML loop in `AdvanceADER` NOT cached** — the plan defers it (inert for TPV31,
  no PML); left as the OnTheFly quadrature.

## Known limitations
- **Smaller speedup than Lever 1** (3.7× vs 13–18×): `ComputeVolumeRHS`'s cached
  form still applies the per-element 9×9 flux matrices `A_d^e` (`Σ_k A(c,k)…`),
  which is irreducible real work; only the geometry/quadrature is cached.
- The `Mult` (RK) path also calls `ComputeVolumeRHS`, so it gets the cached path
  too when `--deriv-cache` is set — covered by the equivalence test, opt-in only.

## Files changed
- [modified] `dynamic/elem_derivative_cache.hpp` — `BuildElementVolumeOperators`.
- [modified] `dynamic/wave_operator.hpp` — `elem_volume_op_` member.
- [modified] `dynamic/wave_operator.inl` — `Cached` branch in `ComputeVolumeRHS`;
  `SetDerivMode` builds/frees the volume cache + combined budget.
- [modified] `tests/unit/test_wave_operator_spatial_derivative.cpp` —
  `VolumeRHSCachedVsOnTheFlyMaxRelErr` + 4 assertions.
- [modified] `tests/unit/bench_ader_hotpath.cpp` — `ComputeVolumeRHS` timing row.

## Testing
- **Lever 3 equivalence:** cached vs OnTheFly `ComputeVolumeRHS` — scalar
  2.4e-16/5.5e-16, bimaterial 4.0e-16/5.6e-16 (orders 2/3), ≤ 1e-12. The
  recursion-level test (`AdvanceADER`, which calls `ComputeVolumeRHS`) still
  passes at 1.4e-15, and Lever 2 bit-exact parity stays **0.000e+00** at Cached.
- `seas_test_wave_operator_spatial_derivative` **48/48** (was 44; +4).
- Regression all green: `wave_operator` 23/23,
  `bimaterial_wave_operator_parity` **9/9** (default byte-exact),
  `ader_ck_predictor` 11/11, `ader_linear_wave_equivalence` 4/4,
  `ader_tpv102_smoke` 4/4.
- **Speedup (`make bench`, p2/O3):** `ComputeVolumeRHS` **3.70×** (53.9→14.6 ms);
  Cached `AdvanceADER` 173→132 ms; **combined macro-step 4.63× → 6.00×**.
  Frontera projection: volume RHS ~80 s → ~22 s; re-profile to confirm.
