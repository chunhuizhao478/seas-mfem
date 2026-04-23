# TPV102 Nucleation Code Review

Date: 2026-04-22

Context:
- User-observed Frontera run shows `V_max` staying near `1e-7 m/s` through `t ~ 1.7 s`
- Run banner says `State representation: total-Q (pre-stress baked into Q)` and `Time integrator: ADER-O(2)`
- Goal: compare current MFEM nucleation handling against SeisSol and identify the likely bug suppressing nucleation

## Findings

### 1. High: total-Q nucleation is injected into the evolving bulk state `Q`, but the state tracker assumes that nucleation remains persistently embedded in `Q` across time steps

MFEM’s current total-Q driver does this:

- zeroes fault pre-stress in `DOFData`, so nucleation can no longer act through `tau2_0`
  - `drivers/tpv102_driver.cpp:498-503`
- declares that total-Q nucleation will instead be written into bulk `Q[SXY]`
  - `drivers/tpv102_driver.cpp:521-537`
- for ADER, applies nucleation once per step with
  - `ApplyNucleationTotal(Q, ..., t + dt/2)`
  - `drivers/tpv102_driver.cpp:1086-1100`
- then advances the PDE and swaps in the new bulk state
  - `wave.AdvanceADER(Q, dt_step, ader_order, Q_new);`
  - `Q.Swap(Q_new);`
  - `drivers/tpv102_driver.cpp:1111-1112`

The core problem is in the contract of `ApplyNucleationTotal()`:

- it stores a running `dtau_applied_plus/minus`
  - `dynamic/tpv102_setup_total.hpp:545-562`
- each call only adds the *increment*
  - `delta = dtau_new - dtau_applied`
  - `dynamic/tpv102_setup_total.hpp:568-576`
- then updates the tracker so the next call assumes the old nucleation is still physically present in `Q`
  - `dynamic/tpv102_setup_total.hpp:680-705`

That assumption is not valid once `Q` has been advanced by the wave equation and fault coupling and then replaced by `Q_new`. The PDE step does not preserve a clean decomposition

`Q = background prestress + persistent nucleation prestress + dynamic perturbation`

at those nodal DOFs. After `AdvanceADER(...)`, the previously injected nucleation contribution may have been partially radiated, mixed, or altered by the bulk/fault update, but the tracker still treats it as if it were intact and therefore only adds a tiny increment on the next step.

This is exactly where MFEM diverges from SeisSol’s logic.

In SeisSol, nucleation is applied to a separate persistent fault prestress field:

- `BaseFrictionLaw::evaluate()` calls `adjustInitialStress(...)` before each friction update
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:101-116`
- `adjustInitialStress(...)` increments `initialStressInFaultCS`
  - `SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:422-447`

Crucially, `initialStressInFaultCS` is not the evolving bulk wave state. It is dedicated fault storage. So SeisSol preserves the full nucleation history exactly where the friction solve uses it. MFEM’s current total-Q implementation does not.

Why this can explain the Frontera symptom:

- if the full nucleation amplitude is not being maintained step to step,
- the fault only sees a weak residual perturbation,
- and `V_max` will creep upward very slowly instead of triggering a dynamic nucleation patch.

This is the strongest code-level explanation I see for the reported `V_max` log.

### 2. Medium: the total-Q migration removed the old `DOFData` nucleation path, but the replacement is not semantically equivalent to SeisSol

In the old fluctuation path, nucleation modified `DOFData.tau2_0` directly:

- `dynamic/tpv102_setup.hpp:135-147`

That old path is intentionally disabled for total-Q:

- driver zeroes `sigma_n0/tau1_0/tau2_0`
  - `drivers/tpv102_driver.cpp:498-503`
- `EvaluateTotal()` asserts those fields must remain zero
  - `dynamic/fault_face_flux.cpp:258-274`
- `EvaluateTotal()` computes friction entirely from `Q_plus/Q_minus`
  - `dynamic/fault_face_flux.cpp:294-337`

So in total-Q mode, the only way nucleation can work is if the bulk-state injection is correct.

Comparing to SeisSol:

- SeisSol’s friction law still forms total tractions from
  - `initialStressInFaultCS + faultStresses`
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:135-138`
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:222-225`

That means SeisSol preserves a clean split:

- dynamic wave contribution from the predictor
- persistent initial/nucleation stress on the fault

MFEM’s current total-Q path collapses both into a single evolving bulk `Q`, then tries to maintain the nucleation component with a delta tracker. That is a much weaker construction and, based on the current behavior, likely the wrong one.

### 3. Medium: ADER midpoint freezing of nucleation is a secondary discrepancy from SeisSol, but it is probably not the main reason nucleation fails so badly

MFEM ADER path:

- applies nucleation once at `t + dt/2`
  - `drivers/tpv102_driver.cpp:1086-1100`

SeisSol:

- applies nucleation inside the per-time-index friction loop
  - `SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:98-121`

So SeisSol updates nucleation at each sub-time-step, while MFEM freezes it at a midpoint value over the whole ADER step.

This is a real accuracy difference, but by itself it should cause an `O(dt^2)` timing/amplitude discrepancy, not a nearly complete failure to nucleate. I would treat this as a secondary issue, not the primary root cause.

## Comparison Summary

SeisSol nucleation semantics:

1. Keep fault initial stress in dedicated storage (`initialStressInFaultCS`)
2. Increment that storage at each friction sub-step with `adjustInitialStress(...)`
3. Form total tractions as `initialStressInFaultCS + dynamic faultStresses`
4. Never rely on the evolving bulk predictor state to remember the nucleation history

Current MFEM total-Q semantics:

1. Bake background prestress into bulk `Q`
2. Zero `DOFData` prestress
3. Inject nucleation deltas into bulk `Q[SXY]`
4. Advance the PDE and overwrite `Q` with `Q_new`
5. Assume the previous nucleation amplitude is still present in `Q`, so only inject the delta next step

That step 5 assumption is the weak link.

## Bottom Line

My main conclusion is:

- the current total-Q nucleation implementation is not preserving nucleation in a persistent fault-side state the way SeisSol does
- instead, it writes nucleation into the evolving bulk state and then only applies future increments
- that is very likely why the Frontera run shows `V_max` barely moving instead of properly nucleating

## Recommended Fix Direction

Best fix direction:

1. Stop treating the evolving bulk `Q` as the authoritative store for cumulative nucleation history.
2. Restore a persistent fault-side nucleation/prestress quantity, analogous to SeisSol’s `initialStressInFaultCS`, and combine it with the dynamic fault traction inside the friction solve.
3. If total-Q must be kept, then recompute or re-impose the full current nucleation field each step from a separate persistent source, rather than using `delta = dtau_new - dtau_applied` against the already-evolved `Q`.

In other words: the SeisSol comparison suggests the right abstraction is a persistent fault prestress source, not a delta patch written into the evolving bulk solution vector.
