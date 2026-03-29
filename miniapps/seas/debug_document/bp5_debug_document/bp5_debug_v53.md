# BP5 Debug v53: p=1 Fault Discretization Mismatch with Tandem

**Date**: 2026-03-28
**Status**: ROOT CAUSE IDENTIFIED, NO CODE CHANGES YET
**Previous**: v52 (traction coherence, BLR residual amplification, dip contamination)
**Branch**: `feature/elasticity`

---

## 1. Motivation

After the recent reverts, the previously suspected `dt_init` / v53 startup-path issue is no
longer the leading explanation for the BP5 `p=1`, near-fault 1000 m mismatch.

Two earlier hypotheses were also checked and are no longer the primary suspects:

1. **The IP `3x` penalty factor is wrong at `p=1`**  
   This is **false**. A term-by-term code comparison shows the `3x` factor is the geometric
   normalization required for MFEM's tetrahedral `CalcOrtho/Weight()` convention to match
   Tandem's physical `A/V` penalty coefficient.

2. **MFEM boundary loading misses Tandem's center-line `Vp*t` branch**  
   This is also **false**. MFEM has separate non-fault `Y=0` interior/shared Dirichlet logic
   that reproduces the Tandem center-line jump construction.

The remaining question was: if penalty and boundary loading are not the source, what still
makes MFEM `p=1` disagree with Tandem `p=1`?

The answer is now clear:

> **Current MFEM `p=1` does not solve the same discrete fault problem as Tandem `p=1`.**

MFEM collapses each triangular fault face to **one constant fault DOF** at `p=1`, while
Tandem uses **three nodal fault DOFs per triangular face** at `p=1`.

That changes the sampled BP5 problem data (`a`, `Dc/L`, `V_init`, `tau_pre`) near the sharp
nucleation-zone boundary and makes the comparison non-like-for-like.

---

## 2. Main Finding

### 2.1 Tandem `p=1` Uses Three Fault Nodes per Triangle

Tandem's rate-and-state fault space is always built from
`NodalRefElement<DomainDimension - 1>(PolynomialDegree)`.

Source:

```cpp
auto RateAndStateBase::Space() -> NodalRefElement<DomainDimension - 1u> {
    return NodalRefElement<DomainDimension - 1u>(
        PolynomialDegree, WarpAndBlendFactory<DomainDimension - 1u>(), ALIGNMENT);
}
```

For 3D elasticity, the fault is 2D, so at `PolynomialDegree = 1` a triangular fault face has
**three nodal basis functions / three nodal samples**.

Tandem then assigns friction/problem data at **every fault node**:

```cpp
for (std::size_t index = 0; index < num_nodes; ++index) {
    auto params = pfun(fault_.storage()[index].template get<Coords>());
    law_.set_params(index, params);
}
```

So Tandem `p=1` is already a **multi-DOF fault discretization**.

### 2.2 MFEM `p=1` Uses One Face-Averaged Fault DOF

MFEM explicitly hard-codes:

- `p=1` -> `face_order = 0`
- `face_order = 0` -> `nbf = 1`
- `nbf = 1` means one constant fault DOF per face

Source:

```cpp
// p=1: face_order=0, nbf=1 -> backward compatible (face average)
int face_fe_order = (method_ == DGMethod::IP && order_ >= 2) ? order_ : 0;
face_quad_ = std::make_unique<FaceQuadrature>(face_fe_order, ...);
nbf_per_face_ = face_quad_->NumBasisFunctions();
```

and

```cpp
///   - p=1 (vol_order=1): face_order=0, nbf=1 -> proven stable, backward compatible
///   - p≥2 (vol_order≥2): face_order=vol_order, nbf=(p+1)(p+2)/2 -> matches Tandem
```

So current MFEM only matches Tandem's fault discretization at `p>=2`, not at `p=1`.

### 2.3 Consequence for BP5 Initialization

MFEM evaluates BP5 parameter fields at its fault DOFs. At `p=1`, that means one sample per
face, effectively at the face centroid:

```cpp
// At nbf=1 (p=1): single centroid point -> same as old code
const IntegrationRule &nir = face_quad_->GetNodalRule();
...
fault_x2_(i * nbf + kk) = coords(0);
fault_x3_(i * nbf + kk) = -coords(2);
```

Then:

```cpp
for (int i = 0; i < num_fault_dofs_; i++)
{
   real_t x2 = coords_x2_(i);
   real_t x3 = coords_x3_(i);

   a_values_(i) = bp5_params_.a_of_x2_x3(x2, x3);
   dc_values_(i) = bp5_params_.Dc_of_x2_x3(x2, x3);
   bp5_params_.tau0_vec(x2, x3, tau);
   bp5_params_.V_init_vec(x2, x3, Vi);
}
```

Therefore, on each `p=1` fault triangle:

- Tandem samples BP5 data at **3 face nodes**
- MFEM samples BP5 data at **1 face centroid**

This is not a small implementation detail. It changes the discrete nucleation patch itself.

---

## 3. Why This Matters for the 1000 m BP5 `p=1` Case

BP5 has a **sharp rectangular nucleation zone**:

- depth: `hs + ht <= x3 <= hs + ht + H`
- strike: `-l/2 <= x2 <= -l/2 + w`

In Tandem, those tests are applied to each fault node. In current MFEM `p=1`, they are
applied only at the face centroid.

This means a face that straddles the nucleation boundary can be treated as:

- **partially inside** in Tandem (`3` nodal samples, potentially mixed values), but
- **fully inside or fully outside** in MFEM (`1` centroid sample only)

That directly changes:

1. `Dc = L_nuc` vs `L0`
2. `V_init = V_nuc` vs `Vp`
3. `tau_pre`
4. local friction state through `psi` equilibrium initialization

The result is a different discrete seed geometry and different initial traction/slip-rate
balance near the nucleation front.

For the 1000 m mesh, this is a strong candidate for why Tandem `p=1` nucleates while current
MFEM `p=1` does not.

---

## 4. Secondary Configuration Difference: `bp5_outside`

There is also a smaller, but real, BP5 configuration mismatch:

- Tandem's stock BP5 case uses `scenario = "bp5_outside"`
- `bp5_outside = BP5.new({eps=1e-3})`
- MFEM currently uses exact inclusion logic with no `eps` offset

This affects points lying exactly on the nucleation-zone boundary.

However, this is **not** the primary issue found here. The larger inconsistency is the fault
discretization mismatch itself:

- Tandem `p=1`: 3 nodal samples / face
- MFEM `p=1`: 1 centroid sample / face

---

## 5. Conclusion

The current BP5 `p=1` mismatch is best explained by a **discrete model mismatch**, not by a
penalty-coefficient mismatch.

### 5.1 What Is Now Disproved

1. **`3x` penalty is wrong for `p=1`**  
   Disproved. The `3x` factor is required to match Tandem's IP penalty coefficient.

2. **MFEM is missing Tandem's `Y=0` center-line loading branch**  
   Disproved. That logic exists in the interior/shared Dirichlet assembly.

3. **The reverted `bp5_v53_p1_exact_7618935.out` startup path is still the active blocker**  
   No longer supported after the reverts.

### 5.2 What Is Now the Leading Root Cause

> **MFEM `p=1` fault discretization is not Tandem `p=1` fault discretization.**

Current MFEM uses one face-averaged fault DOF per triangle at `p=1`.  
Tandem uses three nodal fault DOFs per triangle at `p=1`.

As a result, MFEM and Tandem are not sampling the same BP5 fields on the same discrete fault
space, so the observed nucleation mismatch is unsurprising.

---

## 6. Proposed Fix

### 6.1 Goal

Achieve **exact Tandem-style fault sampling at `p=1`**:

- triangular fault face
- polynomial degree `p=1`
- **3 nodal fault DOFs per face**
- BP5 parameters evaluated at those 3 face nodes
- slip interpolated from those 3 DOFs to quadrature points
- traction projected back onto those 3 DOFs

In short:

> **MFEM `p=1` should use the same nodal fault space as Tandem `p=1`, not a face average.**

### 6.2 Concrete Implementation Target

For the IP elasticity path:

1. Remove the special-case reduction
   - current: `p=1 -> face_order=0`
   - target: `p=1 -> face_order=1`

2. Construct `FaceQuadrature` with `face_order = order_` for `p=1` as well
   - for triangles, this gives `nbf = (1+1)(1+2)/2 = 3`

3. Keep per-DOF BP5 parameter evaluation exactly as already written
   - once `nbf_per_face_ = 3`, the existing `GetFaultCoords2D()` and
     `ComputeBP5Params()` machinery will naturally sample all three nodes

4. Keep the fault traction/slip machinery in nodal form
   - interpolate nodal slip to quad points
   - compute traction at quad points
   - Galerkin project traction back to nodal DOFs

This is already the design of the multi-DOF fault path. The required fix is mainly to make
that path active for `p=1`, not only for `p>=2`.

### 6.3 Exact Matching Target

To match Tandem as closely as possible, the intended final state is:

- `p=1` fault faces use **3 nodal samples per face**
- node locations follow the same reference-face nodal set as Tandem's
  `NodalRefElement<2>(1)`
- BP5 parameter evaluation is nodal, not centroidal
- `bp5_outside`-style inclusion tolerance should also be considered after the nodal fix

The nodal fault-space change is the first-priority correction.

---

## 7. Recommended Validation Plan After the Fix

After implementing the `p=1 -> 3 DOF / face` change:

1. Print `nbf_per_face_` at startup for `p=1`
   - expected: `3`

2. Dump a few fault-face coordinates at the nucleation boundary
   - confirm three distinct nodal points per face, not one centroid

3. Compare nucleation-zone classification counts
   - number of DOFs with `Dc = 0.13`
   - number of DOFs with `V_init = 0.01`
   - compare with Tandem

4. Re-run BP5 near-fault 1000 m, `p=1`
   - compare station nucleation timing against Tandem

5. Only after this, revisit secondary differences
   - `bp5_outside` vs exact inclusion
   - any remaining strike/dip traction differences

---

## 8. Bottom Line

The current BP5 `p=1` comparison is not apples-to-apples.

The critical mismatch is:

- **Tandem `p=1`: 3 fault nodes per triangular face**
- **MFEM `p=1`: 1 face-average DOF per triangular face**

The proposed fix is to make MFEM `p=1` use the same nodal fault discretization as Tandem,
with **three node samples per element face exactly**.

---

## 9. Implemented in Code

### 9.1 Fix 1: p=1 IP Fault Space Now Uses 3 Nodal DOFs per Triangle

This fix has now been implemented.

Current MFEM IP behavior:

- `face_fe_order = order_` for all IP orders, including `p=1`
- triangular fault faces at `p=1` therefore use `nbf = 3`
- BP5 parameters are evaluated at those 3 face nodes

This removes the previous MFEM-only special case:

- old: `p=1 -> face_order=0 -> nbf=1`
- new: `p=1 -> face_order=1 -> nbf=3`

Unit-test status after the change:

- `make test-elasticity-operator`
- result: **300 passed, 0 failed**

### 9.2 Fix 2: BP5 Nucleation Inclusion Now Matches Tandem `bp5_outside`

This fix has also now been implemented.

MFEM previously used exact rectangular inclusion:

```cpp
x3 >= hs+ht && x3 <= hs+ht+H &&
x2 >= -l/2 && x2 <= -l/2+w
```

Tandem stock BP5 uses:

- `scenario = "bp5_outside"`
- `eps = 1e-3`

so boundary points are classified with a small outward tolerance.

MFEM now exposes the same behavior through:

- `BP5Params::nucleation_eps = 1e-3` by default
- CLI override: `--nucleation-eps`
- `--nucleation-eps 0.0` reproduces the old exact behavior

Unit-test status after the change:

- `make test-bp5-params`
- result: **123 passed, 0 failed**

Verification build status:

- `make seas_bp5_full`
- result: build succeeds

---

## 10. Current Outcome After Both Fixes

### 10.1 What Improved

The `p=1` nodal-fault fix clearly improved the initial nucleation behavior.

At the nucleation-front station `strk-24dp+10`, MFEM and Tandem are now close in the
early phase:

- at `t ~ 1 s`, both are near `log10(V_strike) ~ -1.91`
- at `t ~ 5 s`, both are near `log10(V_strike) ~ -1.86`

So the local seed behavior is no longer the main failure mode.

### 10.2 What Still Fails

The run still decays relative to Tandem one station ahead of the front.

At `strk-16dp+10`:

- MFEM near `t ~ 1 s`: `log10(V_strike) = -9.056`
- Tandem near `t ~ 1 s`: `log10(V_strike) = -8.975`

At `t ~ 5 s`:

- MFEM: `log10(V_strike) = -9.139`
- Tandem: `log10(V_strike) = -8.567`

So the remaining mismatch appears in **stress transfer / front propagation**, not in the
local nucleation seed itself.

### 10.3 `bp5_outside` Fix Verdict

The `bp5_outside` inclusion fix was worth implementing for consistency, but:

> **It is not sufficient to fix the remaining mismatch.**

The user has already reported that this change is **not helping** the overall `p=1`
Tandem comparison. The current evidence agrees: the main remaining error is no longer
boundary-node classification of the seed region.

---

## 11. Updated Diagnosis

After the two implemented fixes:

1. `3x` penalty mismatch  
   still **disproved**

2. `p=1` face-average vs 3-node fault discretization mismatch  
   **fixed**

3. `bp5_outside` vs exact nucleation inclusion mismatch  
   **fixed**, but **not enough**

The leading remaining problem is now:

> **MFEM nucleates locally, but the early propagation/loading transfer toward adjacent
> stations is weaker than Tandem.**

The most defensible remaining target is now:

- traction recovery / stress transfer beyond the seed

The previously considered alternatives are weaker:

- `elastic_sigma_n` was already addressed and rejected in earlier BP5 debug notes
- BLR vs exact MUMPS does not change the observed decay in the user's current runs
- `p=1` IP serial-vs-parallel shared-face inconsistency is not supported by the new test below

---

## 12. New Check: BP5 IP Serial-vs-Parallel Shared-Face Diagnostic

To test whether the production mismatch comes from the parallel shared-face IP path,
a new diagnostic was added to:

- `miniapps/seas/tests/parallel/test_bp5_parallel_smoke.cpp`

The new check:

- uses `bp5/mesh/reference/bp5_tandem_coarse.msh`
- runs BP5 with `p=1`, `DGMethod::IP`
- compares serial vs parallel on the same reference mesh
- checks both:
  - shared-face duplicate agreement before deduplication
  - serial vs parallel fault fields after deduplication
- evaluates two states:
  - post-initialization equilibrium state
  - a deterministic nonuniform slip perturbation state

Verification status:

- `make test-bp5-smoke`
- result: **passes on both `np=2` and `np=4`**

Observed diagnostic values:

- `np=4`, initialization:
  - duplicate shared-face keys: `3`
  - max duplicate relative spread in strike traction / strike slip rate / strike RHS:
    `0 / 0 / 0`
- `np=4`, perturbed slip state:
  - duplicate shared-face keys: `3`
  - max duplicate relative spread:
    about `1.7e-14 / 0 / 0`
  - serial-vs-parallel relative difference in strike traction / strike slip rate /
    strike RHS:
    about `9.3e-8 / 6.7e-11 / 6.7e-11`

Verdict:

> **The BP5 `p=1` IP shared-face parallel path is not showing a meaningful
> serial-vs-parallel inconsistency on the reference coarse mesh.**

So the current Tandem mismatch is **unlikely** to be primarily an MPI shared-face sign
or assembly bug.

---

## 13. New Check: Dip Contamination Still Persists in `v53`

The new `p=1` fault-space fix improved local nucleation at `strk-24dp+10`, but the
current `v53` station comparison still shows excessive dip contamination at the next
station ahead of the front, `strk-16dp+10`.

Comparing MFEM `v53` against Tandem:

At approximately `0.1 s`:

- MFEM:
  - `tau_dip = -4.35e-4 MPa`
- Tandem:
  - `tau_dip = -4.77e-5 MPa`
- ratio:
  - MFEM dip traction magnitude is about **9.1x** Tandem

At approximately `1.0 s`:

- MFEM:
  - `tau_dip = -4.33e-3 MPa`
- Tandem:
  - `tau_dip = -2.52e-4 MPa`
- ratio:
  - MFEM dip traction magnitude is about **17.2x** Tandem

At approximately `5.0 s`:

- MFEM:
  - `tau_dip = -2.86e-2 MPa`
- Tandem:
  - `tau_dip = +3.33e-3 MPa`
- ratio:
  - MFEM dip traction magnitude is still about **8.6x** Tandem

The strike traction stays much closer over the same interval, while dip traction
diverges strongly.

The dip slip-rate ratio shows the same pattern:

- at `~1.0 s`, `|V_dip / V_strike|`:
  - MFEM: about `2.2e-4`
  - Tandem: about `1.6e-5`
  - MFEM is about **14x** larger

Verdict:

> **The remaining `p=1` mismatch is still most consistent with traction
> contamination, especially excess dip traction ahead of the nucleation front,
> not with the already-fixed seed discretization issue.**

---

## 14. Current Best Next Step

The previous proposal to use `--traction-stress-only` as a likely fix was too strong.
That run has now been checked directly and it blows up very quickly.

Observed behavior from the user's full run with `--traction-stress-only`:

- earthquake begins effectively at `t=0`
- `V_max` grows rapidly instead of decaying
- the run becomes strongly unstable

Correct interpretation:

> **The penalty correction cannot simply be removed.**
>
> It is providing essential stabilization, even though the full traction still
> appears too contaminated to match Tandem.

So the right next step is not "stress-only as the solution." It is:

- **measure the stress part and correction part separately at the BP5 stations**

### 14.1 New Diagnostic Implemented

A new station-level traction decomposition diagnostic has been implemented:

- `ElasticityDomainOperator::ComputeTractionComponents(...)`
  now returns:
  - total traction
  - stress contribution
  - correction contribution
- `BP5BenchmarkOutput` / `ParallelBP5BenchmarkOutput`
  can now write separate station files for this decomposition
- BP5 driver flag:
  - `--diag-station-traction-decomp`

This writes files of the form:

- `<prefix>_tracdec_<station>.txt`

with columns:

- `time(s)`
- `tau_stress_strike(MPa)`
- `tau_stress_dip(MPa)`
- `tau_corr_strike(MPa)`
- `tau_corr_dip(MPa)`
- `tau_total_strike(MPa)`
- `tau_total_dip(MPa)`

The decomposition is evaluated only on actual BP5 output writes, so it does not
change the solve path.

### 14.2 Updated Best Isolation Run

- rerun the same `p=1`, near-fault 1000 m case with:
  - `--diag-station-traction-decomp`

Reason:

- it preserves the actual stable run
- it exposes whether the remaining mismatch at `strk-16dp+10` comes mainly from:
  - the stress part,
  - the correction part,
  - or their balance

---

## 15. Result From Station-Level Traction Decomposition

The new decomposition files have now been checked directly from the user's run:

- `bp5_v53a_ip_p1_tracdec_tracdec_fltst_strk-24dp+10.txt`
- `bp5_v53a_ip_p1_tracdec_tracdec_fltst_strk-16dp+10.txt`
- `bp5_v53a_ip_p1_tracdec_tracdec_fltst_strk+00dp+10.txt`

### 15.1 Main Finding

At the critical near-front station `strk-16dp+10`, the correction term is not a
small adjustment.

It dominates both strike and dip traction from the beginning of the run.

Examples from `strk-16dp+10`:

At about `0.106 s`:

- strike:
  - `tau_stress = 9.87e-04 MPa`
  - `tau_corr   = 4.06e-03 MPa`
  - `tau_total  = 5.05e-03 MPa`
  - correction / stress magnitude = about `4.1`

- dip:
  - `tau_stress = -2.86e-04 MPa`
  - `tau_corr   = -1.03e-03 MPa`
  - `tau_total  = -1.31e-03 MPa`
  - correction / stress magnitude = about `3.6`

At about `0.984 s`:

- strike:
  - `tau_stress = 9.29e-03 MPa`
  - `tau_corr   = 2.85e-02 MPa`
  - `tau_total  = 3.78e-02 MPa`
  - correction / stress magnitude = about `3.1`

- dip:
  - `tau_stress = -1.99e-03 MPa`
  - `tau_corr   = -9.32e-03 MPa`
  - `tau_total  = -1.13e-02 MPa`
  - correction / stress magnitude = about `4.7`

At about `4.99 s`:

- strike:
  - `tau_stress = 4.20e-02 MPa`
  - `tau_corr   = 1.15e-01 MPa`
  - `tau_total  = 1.57e-01 MPa`
  - correction / stress magnitude = about `2.75`

- dip:
  - `tau_stress = -7.21e-03 MPa`
  - `tau_corr   = -4.23e-02 MPa`
  - `tau_total  = -4.95e-02 MPa`
  - correction / stress magnitude = about `5.87`

So at `strk-16dp+10`, about `75%` to `85%` of the total traction comes from the
correction term rather than from the physical stress term.

### 15.2 Contrast With Other Stations

At `strk-24dp+10`:

- early strike traction is still mostly stress-driven
- for example, at about `0.106 s`, strike correction / stress is only about `0.06`
- but dip is already strongly correction-dominated

At `strk+00dp+10`:

- both components are much smaller
- strike correction / stress is only about `0.19` to `0.26`
- dip correction / stress is about `0.62` to `0.65`
- so the correction is present, but not dominating the way it does at `-16 km`

### 15.3 Interpretation

This sharpens the diagnosis:

> **The remaining `p=1` mismatch is primarily a front-local traction-balance
> problem.**

More specifically:

- `strk-24dp+10` nucleation is reasonably stress-driven in strike
- but one station ahead, `strk-16dp+10`, the correction term overwhelms the
  physical stress term in both strike and dip
- that is consistent with the observed mismatch against Tandem at the front

This also explains why `--traction-stress-only` is not a fix:

- without the correction term, the run blows up
- so the correction is providing essential stabilization
- but in the current formulation it is also too large near the front, so the
  stable total traction is not Tandem-like

### 15.4 Updated Conclusion

The most defensible current statement is:

> **The remaining bug is not primarily in MPI/shared-face handling, not in BLR,
> and not in the already-fixed `p=1` fault discretization.**
>
> **It is in the composition of the total traction passed to friction near the
> rupture front, where the correction term dominates too strongly.**

---

## 16. Next Diagnostic: Jump Residual At Stations

The next missing quantity is the actual fault jump residual that feeds the
penalty correction:

- `R = [[u]] - delta`

The traction decomposition alone shows that the correction dominates near
`strk-16dp+10`, but it does not yet tell us whether that happens because:

- the residual `R` itself is large there, or
- `R` is modest and the penalty scaling is amplifying it too strongly.

### 16.1 Code Change

A new station-level diagnostic has now been added:

- driver flag:
  - `--diag-station-jump-residual`

- new output files:
  - `<prefix>_jumpres_<station>.txt`

- columns:
  - `time(s)`
  - `jump_res_strike(m)`
  - `jump_res_dip(m)`
  - `jump_res_mag(m)`

Implementation notes:

- `ElasticityDomainOperator` now exposes
  - `ComputeTractionDiagnostics(...)`
  which returns:
  - total traction
  - stress part
  - correction part
  - jump residual in local `[dip, strike]` basis

- the BP5 output path and parallel wrapper now support writing the new station
  files on the same schedule as the traction-decomposition files

- the initial half-finished residual implementation had scope bugs in the IP
  shared-face/interior code; those have been repaired and the executable now
  builds again

### 16.2 Verification

Build verification completed:

- `PATH=/Users/chunhuizhao/miniforge/envs/mfem-dev/bin:$PATH make -C miniapps/seas seas_bp5_full`

### 16.3 Purpose

This diagnostic is intended to answer the next concrete question:

> **At `strk-16dp+10`, is the oversized correction caused mainly by a large
> displacement-jump residual, or by the penalty multiplying a residual that is
> already small enough that Tandem would not be correction-dominated there?**

That is the next decision point for narrowing the remaining mismatch.
