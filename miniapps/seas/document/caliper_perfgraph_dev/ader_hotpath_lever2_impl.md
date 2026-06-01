# Implementation Complete: Lever 2 (shared CK recursion)

## What was implemented
- **Merged recursion** `WaveOperator::ComputeADERSubStepStatesAndIntegral`
  (`wave_operator.{hpp,inl}`) — runs the Cauchy-Kovalevskaya recursion ONCE and
  accumulates BOTH the substep nodal Taylor states (weights τ^k/k!) and the time
  integral I (weights dt^{k+1}/(k+1)!), each with its OWN factorial phase
  (R-003: substep denom=k+1, integral denom=k+2 — NOT unified).
- **`AdvanceADER` overload** — added `const Vector *I_precomputed = nullptr`;
  when non-null it reuses the supplied I instead of recomputing
  `ComputeADERTimeIntegrated`, skipping the redundant second recursion. Default
  null ⇒ existing 4-arg call sites unchanged.
- **Driver wiring** (`spatial_dyn_driver.cpp`) — `--shared-ck-recursion` opt-in
  (default OFF). `AdvanceADERWithSubStep_Spatial` gains a `use_shared_ck` param;
  when set it calls the merged routine (producing `Q_per_node` + `shared_I`) then
  `AdvanceADER(…, &shared_I)`. The branch is taken identically on **all ranks**
  (same CLI value), so `Q_per_node` still feeds the per-substep
  `ExchangeFaceNbrData` collective — **R-001** (no rank-conditional skip).
- **Parity test** (`test_wave_operator_spatial_derivative.cpp`) +
  **bench** measurement (`bench_ader_hotpath.cpp`).

## Decisions made
- **Option (a) from the plan** (extend the substep routine to also emit I + an
  AdvanceADER overload), not (b) (a private `ComputeCKStates`). Reason: (a) maps
  1:1 onto the existing driver call sites and keeps both originals intact for the
  parity gate; minimal surface.
- **`--shared-ck-recursion` default OFF**, matching the staged opt-in pattern of
  Lever 1's `--deriv-cache`. Lever 2 is byte-identical so default-on would also be
  safe, but opt-in keeps blast radius minimal until Frontera-validated. The two
  flags are independent and compose (`--deriv-cache --shared-ck-recursion`).
- **Reused `ck_substep_*_buf_`** for the merged routine's D(k) scratch (it
  subsumes the substep predictor); bit-identical regardless of buffer choice.

## Known limitations
- **Lever 2's standalone benefit shrinks once Lever 1 is on.** Measured macro-step
  (no-fault micro-bench): merge alone is **1.70×/1.73×** in OnTheFly mode but only
  **1.15×/1.14×** in Cached mode (p2/p3) — because with the cache the recursion is
  already cheap and the corrector (`ComputeVolumeRHS` + `ComputeADERFaceFluxRHS`)
  dominates. The corrector is Lever 3/4 territory, untouched here.
- On Frontera the real corrector (with the fault face flux, ~23% and 4×
  imbalanced) is a larger share than in the no-fault bench, so the end-to-end
  speedup will be below the bench's combined figure — confirm by re-profiling.
- The driver end-to-end (`--shared-ck-recursion`) path is not run locally
  (no-local-full-mesh); verified by compile + the bit-exact parity unit test.

## Files changed
- [modified] `dynamic/wave_operator.hpp` — `ComputeADERSubStepStatesAndIntegral`
  decl; `AdvanceADER` `I_precomputed` param.
- [modified] `dynamic/wave_operator.inl` — merged routine; `AdvanceADER` reuse-I.
- [modified] `drivers/spatial_dyn_driver.cpp` — `--shared-ck-recursion`;
  `AdvanceADERWithSubStep_Spatial(use_shared_ck)`.
- [modified] `tests/unit/test_wave_operator_spatial_derivative.cpp` —
  `TestSharedCKParity` (bit-exact, both DerivModes, orders {2,3,4}).
- [modified] `tests/unit/bench_ader_hotpath.cpp` — merged macro-step timing +
  Lever-1/Lever-2/combined speedup report.

## Testing
- **Lever 2 parity (the gate):** merged == separate **bit-for-bit (0.000e+00)**
  for `Q_per_node` AND `Q_new`, DerivMode ∈ {OnTheFly, Cached}, orders {2,3,4}.
- `seas_test_wave_operator_spatial_derivative` **44/44** (was 38; +6 parity).
- Plan §5 regression set all green: `wave_operator` 23/23,
  `bimaterial_wave_operator_parity` **9/9** (default byte-exact),
  `ader_ck_predictor` 11/11, `ader_linear_wave_equivalence` 4/4,
  `ader_tpv102_smoke` 4/4.
- Driver builds + links; `--shared-ck-recursion` present in the binary.
- **Speedup (`make bench`):**

  | | p2/O3 | p3/O3 |
  |---|---|---|
  | Lever 1 (cache) on 2-recursion macro-step | 4.02× | 5.22× |
  | Lever 2 (merge) — OnTheFly / Cached | 1.70× / 1.15× | 1.73× / 1.14× |
  | **Combined (OnTheFly 2-rec → Cached merged)** | **4.63×** | **5.93×** |
