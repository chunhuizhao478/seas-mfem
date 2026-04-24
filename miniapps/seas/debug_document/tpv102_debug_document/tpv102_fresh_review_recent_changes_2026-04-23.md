# TPV102 Fresh Review of Recent Changes and Next-Step Instrumentation

Date: 2026-04-23

## Findings

1. The recent split-prestress refactor is directionally closer to SeisSol, but it does **not** fix the live `test_adjacent_triangle_fault_uniformity` failure.
2. The current failure sequence does **not** support “fault prestress stored in total-Q” as the sole root cause anymore.
3. The strongest live clue is this: on the 8-triangle fixture, the fault QPs are still uniform at step 0, but the bulk DOFs adjacent to the fault are already patterned after the first update. The fault spread appears only later.
4. That shifts the highest-priority suspicion from “wrong nucleation storage” or “wrong friction input contract” to the **fault-to-bulk update path**, especially the ADER fault integration shortcut and the face-flux deposition/lift path.
5. SeisSol differs from MFEM not only in prestress storage, but also in **how imposed fault states are time-integrated and injected back into the bulk**. That difference is now a higher-value target than the earlier face-averaging idea.

## What I Reviewed

Recent MFEM changes reviewed:

- [fault_face_flux.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:111)
- [fault_face_flux.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:441)
- [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1119)
- [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1957)
- [tpv102_driver.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/drivers/tpv102_driver.cpp)
- [tpv102_setup.hpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/tpv102_setup.hpp)
- [tpv102_setup_total.hpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/tpv102_setup_total.hpp)
- [test_adjacent_triangle_fault_uniformity.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_adjacent_triangle_fault_uniformity.cpp)
- [test_persistent_nuc_prestress_channel.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_persistent_nuc_prestress_channel.cpp)
- [test_tpv102_pepper_reproducer.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_tpv102_pepper_reproducer.cpp)

Relevant SeisSol paths reviewed:

- [BaseFrictionLaw.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:72)
- [BaseFrictionLaw.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:105)
- [BaseFrictionLaw.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:158)
- [FrictionSolverCommon.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:144)
- [FrictionSolverCommon.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:295)
- [FrictionSolverCommon.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:422)
- [RateAndState.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:135)
- [RateAndState.h](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:222)
- [Neighbor.cpp](/Users/chunhuizhao/projects/SeisSol/src/Kernels/LinearCK/Neighbor.cpp:93)
- [CellLocalMatrices.cpp](/Users/chunhuizhao/projects/SeisSol/src/Initializer/CellLocalMatrices.cpp:804)

## Tests Run

- `make -C miniapps/seas test-adjacent-triangle-fault-uniformity-serial`
  Result: **failed**
- `make -C miniapps/seas test-persistent-nuc-prestress-channel`
  Result: **passed 20/20**
- `make -C miniapps/seas test-tpv102-pepper-reproducer`
  Result: **passed 6/6**

## Fresh Read of the Recent Changes

### 1. The split-prestress fault contract itself looks internally correct

The recent changes in [FaultFaceFlux::Evaluate()](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:116) are reasonable and are the right direction for SeisSol-style fault input semantics:

- total traction for friction is now
  `static prestress + persistent nucleation + dynamic trial traction`
- imposed states still use the **trial-scale corrected traction**, which is consistent with SeisSol’s split between `faultStresses` / `tractionResults` and `initialStressInFaultCS`
- the isolated regression [test_persistent_nuc_prestress_channel.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_persistent_nuc_prestress_channel.cpp) passes, including the important “no tau2_nuc leak into the Riemann jump” check

I do **not** see the current failure as “persistent nucleation channel is still wrong”.

### 2. The production path now matches SeisSol better on representation, but not on full operator flow

MFEM now dispatches the fault through fluctuation-Q:

- RK4 path calls [fault_flux_->Evaluate(...)](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1129)
- ADER path calls [fault_flux_->EvaluateADER(...)](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1960)

This is closer to SeisSol than the previous total-Q path because SeisSol also uses:

- dynamic traces from `qInterpolatedPlus/Minus`
- plus persistent prestress/nucleation in `initialStressInFaultCS`

But SeisSol does **more** than that:

- it computes `faultStresses` from per-substep `qInterpolated±` arrays
- it updates persistent prestress via `adjustInitialStress(...)`
- it computes `tractionResults`
- then it builds **time-integrated imposed states** in [postcomputeImposedStateFromNewStress(...)](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:295)
- then those imposed states are injected into the bulk through precomputed `fluxSolver` / `nodalFlux` kernels in [Neighbor.cpp](/Users/chunhuizhao/projects/SeisSol/src/Kernels/LinearCK/Neighbor.cpp:93) and [CellLocalMatrices.cpp](/Users/chunhuizhao/projects/SeisSol/src/Initializer/CellLocalMatrices.cpp:804)

MFEM still differs in two important ways:

- ADER uses the shortcut in [EvaluateADER()](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:441): divide `I` by `dt`, call `Evaluate()` once on the average state, multiply back by `dt`
- fault-to-bulk deposition is done on the fly by per-QP runtime loops in [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1155) and [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1976), not by SeisSol-style precomputed face kernels

## Most Important Live Observation from the Failing Adjacent-Triangle Test

The current split-prestress branch still fails, but **how** it fails matters:

### Step 0

At step 0, the test prints:

- all 24 fault QPs have identical `V2`, `tau2_corr`, `sigma_n_corr`
- `tau1_corr` is exactly zero everywhere

So the fault solve is initially symmetric under the current split-prestress path.

### After the first update

At step 1:

- the fault outputs are still uniform
- but the bulk DOFs in the fault-adjacent tetrahedra are already patterned
- `SXY` shows alternating per-tet values
- `SXZ` is no longer zero and differs strongly between symmetry-equivalent adjacent tets

That is the most useful current fact.

It means the first observable asymmetry is appearing in the **bulk update next to the fault**, not in the initial per-QP friction state.

### Later steps

By step 2:

- `slip_rate` spread becomes nonzero
- `tau2_corr` spread becomes nonzero but remains tiny
- `tau1_corr` becomes nonzero first, initially at machine-noise scale

By step 19:

- `tau1_corr` spread grows to about `2.236`
- `slip_rate` spread grows to about `4.56e-5`
- `tau2_corr` spread remains only `1.74e-6`
- `sigma_n_corr` spread remains only `2.05e-6`

This pattern is not “everything is equally contaminated”. The dominant visible growth is in the **shear-1 / cross-channel response**.

## Why This Weakens the Old “Total-Q ULP Seed Is the Root Cause” Story

The older ULP-seed story is still plausible as one numerical ingredient, but it is no longer sufficient as the root-cause explanation for the current branch:

- [test_tpv102_pepper_reproducer.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_tpv102_pepper_reproducer.cpp) still shows that a uniform large scalar can reconstruct with a 1-ULP per-QP difference
- but the production path no longer reconstructs 120 MPa background prestress from bulk `Q` at the fault
- the adjacent-triangle test still fails anyway

So:

- the total-Q reconstruction issue may have been a real amplifier on the old branch
- but the **current** failing mechanism survives the split-prestress change
- therefore it cannot be the whole bug

This is the main reason I would not treat the v9.4.0 `tpv102_pepper_reproducer` gate as sufficient evidence for the architectural direction by itself.

## Fresh Comparison with SeisSol

### What MFEM now matches

MFEM now matches SeisSol on the following high-level contract:

- friction takes dynamic trial traction from the bulk traces
- background prestress and nucleation live in persistent fault storage
- total traction for the friction law is `persistent + dynamic`

That part is good.

### What MFEM still does differently in ways that matter for this failure

#### A. ADER fault integration is not SeisSol-style

SeisSol:

- keeps `qInterpolatedPlus/Minus` for multiple time substeps
- updates nucleation inside the time-substep loop
- computes `faultStresses` and `tractionResults` per substep
- builds imposed states by time-weighted accumulation in [postcomputeImposedStateFromNewStress(...)](/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:295)

MFEM ADER:

- receives time-integrated `I_plus/I_minus`
- divides by `dt`
- calls `Evaluate()` once on the average state in [EvaluateADER()](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:441)
- multiplies the imposed state back by `dt`

That is a materially different operator, not just a storage difference.

#### B. Fault-to-bulk injection is runtime and per-QP in MFEM, precomputed in SeisSol

SeisSol:

- uses precomputed `fluxSolver` / `nodalFlux` kernels for dynamic rupture
- the face-to-volume injection path is matrix-based and fixed

MFEM:

- rotates at runtime
- calls `flux_.Interior(...)` at runtime
- deposits with explicit loops
  `rhs +=/-= w * shape(i) * F_h[c]`
  in [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1155) and [wave_operator.inl](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl:1976)

That runtime deposition path is exactly where the adjacent-triangle test first shows nonuniformity in the bulk.

#### C. SeisSol keeps the fault path in point arrays; MFEM reconstructs and re-lifts face data every call

This is not about face averaging. SeisSol is still pointwise on the face. The difference is that SeisSol’s pointwise path is built around precomputed `QInterpolated`, `imposedState`, and `fluxSolver` arrays, whereas MFEM rebuilds face-local traces, rotations, and lifts on each face visit.

## Independent Assessment of the Current Test Suite

### Useful tests

- `test_persistent_nuc_prestress_channel`
  Value: proves the split-prestress fault contract itself is sensible.

- `test_adjacent_triangle_fault_uniformity`
  Value: currently the best live reproducer because it exposes the first nonuniform bulk step on a multi-triangle planar fault.

### Tests that are useful but currently over-interpreted

- `test_tpv102_pepper_reproducer`
  Value: proves one ULP-scale interpolation effect exists on this platform.
  Limitation: does **not** prove that effect is the active cause of the current split-prestress failure.

- the various two-tet pepper diagnostics
  Value: useful for local operator probing.
  Limitation: many of them are tuned to one specific hypothesis and miss the multi-triangle symmetry breaking that the adjacent-triangle test exposes.

## What I Think the Next Instrumentation Should Do

The next tests should focus on the **first nonuniform step**, not on long-time pepper summaries.

### Priority 1: First-step bulk symmetry audit on the 8-triangle fixture

Extend [test_adjacent_triangle_fault_uniformity.cpp](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_adjacent_triangle_fault_uniformity.cpp) with a dedicated “step-0 → step-1 symmetry dump” mode that records, for selected symmetry-equivalent triangles:

- `Q_plus_local`, `Q_minus_local`
- `Q_imp_plus`, `Q_imp_minus`
- `F_h_plus`, `F_h_minus`
- the lifted per-element contribution `w * shape * F_h`
- `rhs` before inverse mass
- `rhs` after inverse mass
- resulting `Q_new` on the fault-adjacent tets

Goal:

- determine the first stage where two symmetry-equivalent triangles diverge
- distinguish whether the divergence starts in:
  - canonical trace reconstruction
  - fault flux solve
  - flux-to-rhs lift
  - inverse mass / full update

### Priority 2: One-step ADER vs one-step RK4 uniformity comparison on the same 8-triangle mesh

Add a companion test:

- same 8-triangle fixture
- same uniform persistent nucleation
- one single step only
- compare:
  - ADER-2
  - RK4 with the same split-prestress fault path

If the first symmetry break is much stronger in ADER than RK4, the focus should move to [EvaluateADER()](/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/fault_face_flux.cpp:441) and the ADER face-flux path first.

This is a higher-value discriminator than more generic perturbation tests.

### Priority 3: SeisSol-style imposed-state reference test for one face

Build a test-only reference operator for one face that follows the SeisSol pattern more closely:

- per-substep `qInterpolated±`
- per-substep `faultStresses`
- per-substep `tractionResults`
- time-weighted imposed-state accumulation

Then compare that reference imposed state against MFEM’s `EvaluateADER(I/dt)` shortcut on the same symmetric inputs.

Goal:

- determine whether the first asymmetry is already present in MFEM’s ADER fault-imposed state, before any face-lift to bulk DOFs

### Priority 4: Fault-to-bulk lift reference test

Build a test that bypasses friction entirely:

- prescribe identical imposed states on symmetry-equivalent fault triangles
- run only the MFEM fault-to-bulk deposition path
- compare the resulting adjacent-tet rhs / `Q_new`

If this alone breaks symmetry, the remaining bug is in the lift/deposition path, not in friction.

### Priority 5: Demote the current ULP reproducer from gate to diagnostic

Keep `test_tpv102_pepper_reproducer`, but do not use it as the primary gate for the fix direction.

It is a platform-dependent numerical clue, not a sufficient explanation of the current live failure.

## Proposed Next Step Order

1. Keep the split-prestress refactor in place for now.
   It is not the final fix, but it is still closer to SeisSol and the isolated fault tests pass.

2. Do **not** jump to another structural fix yet.
   The current evidence is not specific enough to justify a new architectural change.

3. Add the first-step instrumentation to `test_adjacent_triangle_fault_uniformity`.
   This is the most direct route to the actual bug location.

4. Add the one-step ADER-vs-RK4 comparison on the same fixture.
   This is the fastest way to test whether the remaining bug is primarily ADER-specific.

5. If ADER is implicated, compare MFEM `EvaluateADER` against a SeisSol-style time-integrated imposed-state reference next.

## Bottom Line

The recent code changes improved the fault representation and likely fixed a real total-Q/nucleation issue, but they did **not** remove the adjacent-triangle symmetry break. The fresh live evidence says the remaining problem appears first in the bulk update next to the fault, after an initially uniform fault solve. That makes the highest-value next instrumentation:

- first-step bulk symmetry tracing on the 8-triangle fixture
- ADER-vs-RK4 one-step comparison
- a SeisSol-style imposed-state reference comparison

Those three tests should identify whether the remaining bug lives in:

- MFEM’s ADER fault integration shortcut
- the fault-to-bulk face lift / deposition path
- or a smaller residual mismatch in the per-face operator construction

