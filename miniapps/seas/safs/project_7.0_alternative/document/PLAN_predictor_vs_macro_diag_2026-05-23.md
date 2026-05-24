# Implementation Plan: `[MACRO]` predictor-vs-macro diagnostic (test (B)) — 2026-05-23

## Overview

Add a single runtime-env-gated `[MACRO]` stderr trace inside the shared-fault **macro
solve** (`WaveOperator::ComputeADERSharedFaceFluxRHS`, Pass-1 QP loop,
`wave_operator.inl:4073-4185`) that prints, at the seed shared QP, the macro solve's
normal-traction decomposition on the **time-integrated state** `Q_avg = I/dt`. Compared
against the iterator's per-sub-step `[SLIP] sn_vjump` / `[FRAME] dv_n` at the same QP, it
settles **R-008** from `REVIEW_DEBUG_speckle_normal_velocity_jump_2026-05-23.md`: is the
runaway-driving normal-velocity jump born in the **per-sub-step predictor / ghost path**
(iterator collapses, macro bounded) or is it present in the integrated state too?

This is a **diagnostic only**; byte-exact when the env switch is unset.

> **AS-BUILT CORRECTION (post-review R-001, REVIEW.md 2026-05-23).** The first
> implementation mistakenly inserted `[MACRO]` into `ComputeADERFaceFluxRHS` (the **interior**
> fault-flux routine), which never processes the shared seed QP. It has been relocated to
> **`ComputeADERSharedFaceFluxRHS`**, immediately after the `EvaluateADER_LSW` dispatch (~the
> `[XRANK]` seed-locator). The shared path has **no** `states[qq]`/`Q_avg_*_qq` array — the
> decomposition is computed from the canonical time-integrated `I_{plus,minus}_local · (1/dt)`
> (exactly what `EvaluateADER_LSW` consumes), and the printed totals are the **written**
> `fdata.sigma_n_corr`/`fdata.slip_rate` (strictly better than the interior draft's pre-face-
> average `states[qq]`). Sections below that reference `states[qq]`/`Q_avg_*_qq`/`:4184`
> describe the superseded interior draft; the shared-path variables are the as-built.

## Background (the contrast it measures)

Two code paths solve friction on the same shared fault QP from **different inputs**:

| Path | Input state | Produces | Trace |
|---|---|---|---|
| **Iterator** (`tpv205_substep_iterator`, via `EvaluateBulkAtFaultQPsCanonical`) | per-sub-step **predictor** `Q̃±` (CK extrapolation + R-1303/R-1601 ghost coverage) | the **collapsing** σ_n (slip accumulator) | `[SLIP]` `sn_vjump`, `[FRAME]` `dv_n` |
| **Macro solve** (`ComputeADERSharedFaceFluxRHS`) | time-integrated `Q_avg = I/dt` (`:4175-4176`) → `ComputeStageState` (`:4181`) | the **bounded** σ_n **written to output** (`WriteBackState :4292`) | **(this plan) `[MACRO]`** |

§4's observation (`fault.vtkhdf` σ_n bounded ~+49 MPa while the iterator collapses) is only
consistent if the divergence is **not** QP-static (same frame, same QP, same
`FaultFaceFlux`): the difference must be the **input state** (per-sub-step predictor vs
integrated). `[MACRO]` measures the macro side so the two can be compared directly.

The macro normal-traction decomposition mirrors `ComputeTrialTraction` Eq. 7a
(`fault_face_flux.cpp:62-64`):
```
sn_vjump_macro = eta_p·(Q_avg_minus[VX] − Q_avg_plus[VX])
sn_sterm_macro = eta_p·(Q_avg_plus[SXX]/Zp⁺ + Q_avg_minus[SXX]/Zp⁻)
sn_vjump_macro + sn_sterm_macro == states[qq].sigma_n_trial   (self-check, printed)
sigma_n_total_macro = sigma_n0 + sigma_n_nuc + sigma_n_trial  (= states[qq].sigma_n_total)
```

**Read:** `[MACRO] sigma_n_total` stays ≈ +49 MPa (bounded, matching the output min 30.4 MPa)
while `[SLIP] sigma_n_tot` at the same QP collapses ⇒ the opening is a **per-sub-step
predictor / ghost** phenomenon (R-1303/R-1601), confirming R-008 and re-pointing the fix
there. If `[MACRO]` *also* collapses ⇒ the opening is in the integrated state too (not
predictor-specific) — a different locus.

## Constraints

- **Byte-exact when `SEAS_DIAG_MACRO` is unset.** Only a gated `fprintf` is added in Pass 1;
  no change to `Q_avg`, `states`, `BuildImposedState`, `WriteBackState`, or `rhs`. Preserves
  the TPV205/BP5 regression (this routine is shared with TPV* via the template).
- **Runtime env-gated** (mirror `[XRANK]`/`[FRAME]`/`[SLIP]` — `std::getenv`, static parse),
  **not** compile-time `SEAS_DIAG_FAULT_FLUX`.
- **No struct changes.**
- **Single insertion point**, shared-fault macro-solve Pass 1 only.
- **Reuse the seed locator** `SEAS_DIAG_FRAME_XYZ`/`SEAS_DIAG_FRAME_R` (already in the Dc2
  sbatch for `[FRAME]`) so the seed coordinate is set once; new enable `SEAS_DIAG_MACRO`.

## Phase 1: the `[MACRO]` trace

### Files to Modify
- `dynamic/wave_operator.inl` — gated `[MACRO]` trace in `ComputeADERSharedFaceFluxRHS`
  Pass-1 loop, **immediately after** `ComputeStageState` (after `:4184`, inside the
  `for (qq…)` body, before `}` at `:4185`). Headers already present (`<cstdio>/<cstdlib>/
  <cmath>` at `:10-14`; `getenv` already used).
- `jobs/safs/spatial_dyn_resolution_Dc2_8N_400r_dev_2hr_safs.sbatch` — `export
  SEAS_DIAG_MACRO="${SEAS_DIAG_MACRO:-1}"` + echo line (reuses `SEAS_DIAG_FRAME_XYZ/_R`).

### Files to Create
- `tests/unit/test_macro_diag_decomposition_identity.cpp` — asserts
  `sn_vjump_macro + sn_sterm_macro == FaultFaceFlux::ComputeTrialTraction(...).sigma_n_trial`
  (the trace's decomposition faithfully reproduces the production normal channel).

### Detailed Requirements

1. **Env gates** (function-local `static const` lambdas, parsed once):
   - `macro_diag` ← `getenv("SEAS_DIAG_MACRO")` truthy.
   - target `(x,y,z)` ← `getenv("SEAS_DIAG_FRAME_XYZ")` `"%lf,%lf,%lf"` (default
     `607518,3706359,-4543`); radius ← `getenv("SEAS_DIAG_FRAME_R")` (default `300`),
     squared. (Shared seed locator with `[FRAME]`.)

2. **Per-QP gate**: only when `macro_diag`, compute `ftr->Face->Transform(ip_qq, xq)`; emit
   iff `‖xq − target‖² ≤ r²`.

3. **Quantities** (all in-scope at `:4184`: `Q_avg_plus_qq`, `Q_avg_minus_qq`,
   `states[qq]`, `fdata_const`, `dof_idx_per_qp[qq]`, `ip_qq`, `ftr`, `my_rank_`,
   `GetTime()`):
   - `sn_vjump_macro = fdata_const.eta_p·(Q_avg_minus_qq[VX] − Q_avg_plus_qq[VX])`
   - `sn_sterm_macro = fdata_const.eta_p·(Q_avg_plus_qq[SXX]/fdata_const.Zp_plus +
     Q_avg_minus_qq[SXX]/fdata_const.Zp_minus)`
   - print also `Q_avg_plus_qq[VX]`, `Q_avg_minus_qq[VX]`, `states[qq].sigma_n_trial`,
     `states[qq].sigma_n_total`, `states[qq].V_abs`.

4. **Output** — one `fprintf(stderr,…)+fflush`:
   ```
   [MACRO] qp=%d c=(%.1f,%.1f,%.1f) rank=%d t=%.6e vn_plus=%+.6e vn_minus=%+.6e
           sn_vjump=%+.6e sn_sterm=%+.6e sigma_n_trial=%+.6e sigma_n_tot=%+.6e V_abs=%+.6e
   ```
   `sn_vjump+sn_sterm` must equal `sigma_n_trial` on every line (in-situ self-check).

### Edge Cases
- Switch unset ⇒ no coords, no print, byte-exact.
- Serial build ⇒ `my_rank_` member is 0; compiles.
- Malformed `SEAS_DIAG_FRAME_XYZ`/`_R` ⇒ keep defaults.
- Target off this rank's partition ⇒ no output here; the owning rank prints.
- **Sign caveat (document in the trace comment):** the macro `±` routing uses the R-101
  geometric flag (`!elem1_on_plus_per_qp`, `:4119-4135`) while the iterator's predictor path
  uses `qpd.sign_flipped`; the two may differ by a global sign. Compare **magnitudes /
  boundedness**, not the sign of `sn_vjump`.

### Acceptance Criteria
- [ ] `SEAS_DIAG_MACRO` unset ⇒ byte-exact: `make test-ader-interior-vs-shared-branch-live`
      (routes through `ComputeADERSharedFaceFluxRHS`) still rel=0 across orders {2,3,4}.
- [ ] `make seas_spatial_dyn_driver` compiles.
- [ ] `test_macro_diag_decomposition_identity` passes: `sn_vjump+sn_sterm == sigma_n_trial`
      from the real `ComputeTrialTraction` (≤ 16 ULP), incl. `Zp⁺≠Zp⁻`.
- [ ] On Frontera with `SEAS_DIAG_MACRO=1`: `[MACRO]` lines at the seed coords; compare
      `[MACRO] sigma_n_tot` vs `[SLIP] sigma_n_tot` at the same QP/nearest time.

### Dependencies
- Depends on: nothing (self-contained). Complements the shipped `[FRAME]`/`[SLIP]`.
- Required by: settling R-008 (predictor vs integrated locus).

## Testing Strategy
1. Local compile (`make seas_spatial_dyn_driver`).
2. `make test-macro-diag-decomposition-identity` (new) — decomposition faithfulness.
3. `make test-ader-interior-vs-shared-branch-live` — byte-exact (inert when unset).
4. Frontera: `SEAS_DIAG_MACRO=1` short run; grep `[MACRO]` + `[SLIP]`; compare σ_n_tot.

## Risk Assessment
- **`ftr->Face->Transform` cost** — gated behind `macro_diag`; zero in production.
- **Pass-1 vs post-average value** — `states[qq].sigma_n_trial` is the **per-QP** value here
  (Pass 2 face-averaging at `:4190-4259` runs later); per-QP is the right apples-to-apples
  vs the iterator's per-QP `sn_vjump`. (The written output is the face-average; the per-QP
  trial is what we want for the locus question.)
- **`Q_avg = I/dt` semantics** — it is the ADER time-integrated state, not raw base Q; the
  comparison is "integrated-state jump bounded?" vs "per-sub-step-predictor jump collapsing?"
  — exactly the R-008 discriminator. Documented in the trace comment.
