# TPV102 MFEM vs SeisSol Code-Flow Comparison

Date: 2026-04-22

Scope:
- Compare SeisSol (`/Users/chunhuizhao/projects/SeisSol`) and current MFEM TPV102 path (`/Users/chunhuizhao/projects/seas-mfem/miniapps/seas`)
- Focus on overall bulk/fault code flow
- Highlight differences that could plausibly create the observed fault-surface "pepper" pattern

## Executive Summary

My current read is:

1. The strongest architectural difference is not the constitutive law itself. Both codes compute trial fault tractions from plus/minus fault-side states and then apply friction.
2. The highest-risk MFEM area for a pepper pattern is the fault-bulk coupling and face bookkeeping around the fault, especially the repeated per-QP reconstruction of canonical frames, side ordering, and separate interior/shared fault branches in `wave_operator.inl`.
3. The MFEM bulk solver could still be involved, but the fault coupling path is the first place I would look for a cell-scale, triangle-local artifact.
4. The difference between SeisSol ADER and MFEM RK4 is important for parity, but by itself it is not a good explanation for a spatial pepper pattern. It is much better at explaining temporal lag/bias than cellwise speckle.
5. The fact that SeisSol also keeps initial stress in dedicated fault storage means "MFEM uses separate pre-stress bookkeeping" is not, by itself, a sufficient root-cause explanation.

## Compared Code Paths

MFEM:
- `dynamic/wave_operator.inl`
- `dynamic/fault_face_flux.cpp`
- `dynamic/tpv102_setup.hpp`
- `drivers/tpv102_driver.cpp`
- `dynamic/godunov_flux.hpp`

SeisSol:
- `src/Proxy/KernelHost.cpp`
- `src/Kernels/DynamicRupture.cpp`
- `src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h`
- `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h`
- `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h`
- `src/DynamicRupture/Initializer/RateAndStateInitializer.cpp`

## Overall Flow: MFEM

Top-level bulk step:
- `WaveOperator::Mult()` does `ComputeVolumeRHS()`, then `ComputeFaceFluxRHS()`, then `ComputeSharedFaceFluxRHS()` in parallel runs. See `dynamic/wave_operator.inl:483-499`.

Fault coupling in MFEM:
- Fault faces are split into two separate collections: interior fault faces and shared fault faces. See `dynamic/wave_operator.inl:258-320`.
- On an interior fault face, MFEM:
  - interpolates `Q_self` and `Q_nbr` at the face quadrature point,
  - reconstructs a canonical fault basis,
  - rotates both states into canonical coordinates,
  - decides which side is `Q_plus` and `Q_minus`,
  - calls `FaultFaceFlux::Evaluate()` or `EvaluateTotal()`,
  - rotates imposed states back to global coordinates,
  - applies per-side Godunov fluxes back to the bulk RHS.
  See `dynamic/wave_operator.inl:1055-1188`.
- The shared-face path repeats the same basic logic in a second branch, with separate shared-face metadata and rank-side handling. See `dynamic/wave_operator.inl:1608-1665`.

MFEM fault friction solve:
- `FaultFaceFlux::ComputeTrialTraction()` computes `sigma_n_trial`, `tau1_trial`, `tau2_trial` directly from `Q_plus/Q_minus` and per-side impedances. See `dynamic/fault_face_flux.cpp:39-64`.
- `FaultFaceFlux::Evaluate()` then:
  - adds pre-stress from `DOFData`,
  - solves the friction equation,
  - computes corrected tractions,
  - constructs imposed states,
  - writes `slip_rate`, `V1`, `V2`, `tau*_corr`, `sigma_n_corr` back into `DOFData`.
  See `dynamic/fault_face_flux.cpp:70-223`.

MFEM TPV102 state setup:
- TPV102 initializes the bulk `Q` to zero and stores initial fault pre-stress in `DOFData` (`sigma_n0`, `tau1_0`, `tau2_0`). See `dynamic/tpv102_setup.hpp:76-147`.
- Nucleation is applied by mutating `DOFData.tau2_0` before each RK4 stage. See `dynamic/tpv102_setup.hpp:123-147` and `drivers/tpv102_driver.cpp:1175-1236`.

MFEM time stepping:
- TPV102 advances the bulk state with explicit RK4, sampling fault outputs at each stage and then doing an endpoint `wave.Mult(Q, ...)` refresh so instantaneous fault observables match `Q(t+dt)`. See `drivers/tpv102_driver.cpp:1175-1340`.

## Overall Flow: SeisSol

Top-level bulk step:
- SeisSol uses a predictor/integral split:
  - ADER predictor on each cell: `ProxyKernelHostAder::run()` in `src/Proxy/KernelHost.cpp:38-60`
  - local volume integral: `ProxyKernelHostLocal::run()` in `src/Proxy/KernelHost.cpp:140-149`
  - neighbor coupling: `ProxyKernelHostNeighbor::run()` in `src/Proxy/KernelHost.cpp:151-209`

Dynamic rupture entry:
- Dynamic rupture uses a dedicated Godunov/DR path on fault faces, not an inlined branch inside the generic face loop. See `ProxyKernelHostGodunovDR::run()` in `src/Proxy/KernelHost.cpp:249-275`.

SeisSol fault-side state construction:
- `DynamicRupture::spaceTimeInterpolation()` evaluates predicted bulk states at time quadrature points and rotates them to the fault-local frame using precomputed face data (`TinvT`, `plusSide`, `minusSide`, `faceRelation`). See `src/Kernels/DynamicRupture.cpp:51-100`.

SeisSol friction pipeline:
- `BaseFrictionLaw::evaluate()` is face-local and time-quadrature-local:
  - precompute trial stresses from `qInterpolatedPlus/Minus`,
  - run friction updates,
  - postcompute imposed states from the new tractions.
  See `src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:49-168`.
- Trial fault stresses are precomputed by `precomputeStressFromQInterpolated()`. See `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:143-222`.
- Imposed states are constructed by `postcomputeImposedStateFromNewStress()`. See `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:287-365`.
- Rate-and-state uses stored `initialStressInFaultCS` plus the dynamic trial tractions, then updates slip, slip rate, and corrected tractions in one face-local loop. See `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:117-250`.

SeisSol initialization:
- Initial fault stress and RS state live in dedicated per-face/per-point storage such as `InitialStressInFaultCS`, `StateVariable`, `SlipRate1`, `SlipRate2`. See `src/DynamicRupture/Initializer/RateAndStateInitializer.cpp:38-90`.

## Key Architectural Differences

### 1. MFEM embeds fault physics inside the generic face loop; SeisSol has a dedicated DR pipeline

MFEM:
- The generic face loop decides at runtime whether a given face QP is a regular face, boundary face, interior fault face, or shared fault face.
- Fault logic is mixed into the same `ComputeFaceFluxRHS()` / `ComputeSharedFaceFluxRHS()` traversal.

SeisSol:
- Dynamic rupture has its own dedicated face-local data structures and kernels.
- The DR pipeline is naturally separated into:
  - interpolate/rotate fault-side states,
  - compute trial stresses,
  - solve friction,
  - reconstruct imposed states.

Why it matters for pepper:
- MFEM has more moving parts per face QP and more opportunities for a one-face bookkeeping inconsistency.
- SeisSol’s flow is more standardized and less dependent on ad hoc runtime reconstruction at each face visit.

### 2. MFEM has two separate fault branches: interior and shared

MFEM:
- Interior fault branch: `dynamic/wave_operator.inl:1055-1188`
- Shared fault branch: `dynamic/wave_operator.inl:1608-1665`

SeisSol:
- Uses dedicated DR face data and side metadata once in the DR layer, rather than duplicated interior/shared logic in the main solver flow.

Why it matters for pepper:
- A triangle-scale pepper pattern is exactly the kind of artifact that can be produced if interior faces and shared faces are almost, but not exactly, consistent in frame construction, side ordering, or field packing.
- This remains my highest-value MFEM code difference relative to SeisSol.

### 3. MFEM repeatedly reconstructs the canonical frame and `Q_plus/Q_minus` mapping in the face loop

MFEM:
- Reconstructs `can_n`, `can_t1`, `can_t2`, builds `T_can/Tinv_can`, rotates states, and chooses plus/minus at every fault QP visit. See `dynamic/wave_operator.inl:1092-1160` and `1608-1652`.

SeisSol:
- Uses precomputed DR face information (`DRFaceInformation`, `DRGodunovData`, `TinvT`) and a dedicated interpolation kernel. See `src/Kernels/DynamicRupture.cpp:51-100`.

Why it matters for pepper:
- If the pepper is due to a sign flip, a tangent swap, or a one-branch-only mismatch, MFEM’s runtime reconstruction path is the prime suspect.
- SeisSol’s precomputed face metadata reduces that risk substantially.

### 4. MFEM uses mutable `DOFData` as a side-channel beside bulk `Q`

MFEM:
- `DOFData` holds pre-stress, state variable, slip, slip rate, corrected tractions, and output fields.
- `FaultFaceFlux::Evaluate()` mutates `DOFData` every call. See `dynamic/fault_face_flux.cpp:205-211`.
- TPV102 nucleation also mutates `DOFData.tau2_0` stage by stage. See `dynamic/tpv102_setup.hpp:135-147`.

SeisSol:
- Uses dedicated DR layer arrays for `InitialStressInFaultCS`, `StateVariable`, `SlipRate1/2`, `Traction1/2`, `ImposedStatePlus/Minus`. See `src/DynamicRupture/FrictionLaws/FrictionSolver.cpp:23-41` and `src/DynamicRupture/Initializer/RateAndStateInitializer.cpp:38-90`.

Why it matters for pepper:
- A separate mutable side-channel is not wrong, but it creates more opportunities for stale or partially-overwritten per-face state.
- This was already real for the old RK4 output bug. Even after that fix, `DOFData` remains a larger bookkeeping surface than SeisSol’s dedicated DR storage.

### 5. Time integration differs, but this is probably not the main pepper mechanism

MFEM:
- Explicit RK4 on bulk `Q` and fault `psi`, with stage sampling and endpoint refresh. See `drivers/tpv102_driver.cpp:1175-1340`.

SeisSol:
- Time quadrature inside the ADER/DR path:
  - `spaceTimeInterpolation()` gives per-time-point fault states,
  - friction is solved over those time points,
  - imposed states are integrated with `timeWeights`.
  See `src/DynamicRupture/FrictionLaws/CpuImpl/BaseFrictionLaw.h:98-166` and `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:287-365`.

Why it matters:
- This is a strong parity difference.
- But it is better at explaining phase/bias differences than a checkerboard or pepper artifact tied to specific triangles.

## What Is Probably Not the Root Cause

### 1. "MFEM uses separate initial stress storage, SeisSol does not"

This is not accurate enough to be a root-cause claim.

SeisSol also stores initial stress separately in `InitialStressInFaultCS` and adds it inside the friction law:
- `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:135-138`
- `src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:222-225`
- `src/DynamicRupture/Initializer/RateAndStateInitializer.cpp:53-87`

So "separate prestress bookkeeping" alone is not enough to explain the pepper.

### 2. "The constitutive RS law itself differs so much that it explains pepper"

The friction-law structure is actually very similar:
- trial stress from plus/minus states
- add initial stress
- solve for slip rate
- compute corrected tractions
- reconstruct imposed states

That makes a pure constitutive-law bug less likely than a data-flow or mapping bug.

### 3. "RK4 vs ADER alone explains pepper"

Not a strong match for the observed symptom.

Temporal integration differences usually create:
- lag
- bias
- amplitude drift
- smoother field distortion

They are a weaker explanation for cell-scale spatial speckle on the fault.

## Most Plausible Pepper Mechanisms in Current MFEM

### H1. Interior/shared fault branch mismatch

Why plausible:
- MFEM has duplicated fault coupling logic in `ComputeFaceFluxRHS()` and `ComputeSharedFaceFluxRHS()`.
- A subtle mismatch in one branch can produce triangle-local artifacts without completely destroying the run.
- This is the cleanest code-structure difference versus SeisSol.

Code hotspots:
- `dynamic/wave_operator.inl:1055-1188`
- `dynamic/wave_operator.inl:1608-1665`

### H2. Canonical-basis / side-ordering inconsistency at fault QPs

Why plausible:
- MFEM reconstructs `can_n`, `can_t1`, `can_t2`, and `elem1_on_plus` dynamically.
- A sign error or tangent swap would directly perturb `Q_plus/Q_minus`, trial tractions, and imposed states.
- That kind of error can look like a triangle-wise pepper field on fault outputs.

Code hotspots:
- `dynamic/wave_operator.inl:1092-1134`
- `dynamic/wave_operator.inl:1608-1631`

### H3. Fault-side `Q` trace reconstruction mismatch rather than friction-law failure

Why plausible:
- `FaultFaceFlux::Evaluate()` itself is algebraically simple once `Q_plus/Q_minus` are correct.
- If `Q_plus/Q_minus` are slightly wrong per face, all fault observables inherit that error together: `tau`, `sigma_n`, `V`, and downstream slip/state evolution.
- This matches the multi-field pepper symptom better than a bug in only `sigma_n`.

Code hotspots:
- interpolation and rotation in `dynamic/wave_operator.inl:887-913`, `1119-1128`, `1608-1617`
- trial traction in `dynamic/fault_face_flux.cpp:111-170`

### H4. `DOFData` synchronization / ownership / shared-face pairing issue

Why plausible:
- MFEM still has explicit face-to-DOF tracking, shared-face offsets, and fault-surface output synchronization logic.
- A pairing or ownership mismatch can create triangle-local output corruption even if the bulk solve is mostly correct.

Code hotspots:
- fault-face list and basis setup: `dynamic/wave_operator.inl:258-320`
- shared-face output packing: `dynamic/wave_operator.inl:2680-2775`

### H5. Bulk solver issue away from the fault

Why plausible:
- If future bulk outputs show the same pattern in `Q[SYY]` / `Q[SXY]` adjacent to the fault, then the pepper is upstream of the friction solve.

Why not first:
- Right now the most distinctive MFEM-vs-SeisSol complexity is in the fault face bookkeeping, not in the pure volume derivative path.
- `WaveOperator::Mult()` volume part is structurally much simpler than the fault-face branches. See `dynamic/wave_operator.inl:483-499` and `831-925`.

## Ranking

Most likely:
1. MFEM fault-face bookkeeping mismatch between interior/shared paths
2. MFEM canonical frame / plus-minus assignment inconsistency
3. MFEM `DOFData` face ownership or pairing issue

Possible but secondary:
4. Upstream bulk trace interpolation near the fault
5. Nucleation timing/bookkeeping interaction with fault storage

Less likely as primary pepper cause:
6. Rate-and-state constitutive algebra itself
7. RK4 vs ADER time integration alone
8. Separate prestress storage alone

## Practical Conclusions

### Conclusion 1

The SeisSol comparison points much more strongly to MFEM fault-bulk coupling complexity than to the RS constitutive law.

### Conclusion 2

If the pepper survives the corrected RK4 endpoint-refresh run, I would investigate the fault coupling path before blaming the whole bulk solver.

### Conclusion 3

Implementing I-06 may still help robustness and parity, but based on this code comparison it does not look like the single most targeted fix for a pepper artifact.

## Recommended Next Checks

1. Add coarse bulk output near the fault for `Q[SYY]`, `Q[SXY]`, `Q[SXZ]`, and velocity, then compare the bulk field against the fault pepper at the same time.
2. Instrument MFEM to dump `Q_plus_local`, `Q_minus_local`, `can_n`, `can_t1`, `can_t2`, and `elem1_on_plus` for a few noisy triangles on both interior and shared paths.
3. Compare the same physical fault triangle on both ranks across a partition seam and verify bitwise or near-bitwise agreement in the canonical inputs to `FaultFaceFlux::Evaluate()`.
4. Run a locked-fault or no-slip fixture through the same mesh and see whether the pepper remains in trial tractions. If yes, the bug is upstream of the nonlinear friction solve.
5. If bulk fields are smooth but fault outputs are peppered, focus entirely on `wave_operator.inl` fault-path rotation, side ordering, and DOFData mapping.

## Bottom Line

From the SeisSol comparison, the best current explanation for the pepper pattern is:

- not "MFEM bulk elastodynamics are obviously wrong everywhere"
- not "the RS friction formula is fundamentally different"
- but "MFEM has a much more fragile fault-face data path than SeisSol, and that path is exactly where a triangle-scale pepper artifact would most naturally be generated"
