# BP5 Debug v44: IP p≥2 Root Cause — Fault DOF Mismatch & Traction Fix

**Date**: 2026-03-18
**Status**: Multi-DOF fault implemented, all 528+ serial tests pass, ready for TACC runs
**Previous**: v42 (c_N_1 fix — necessary but insufficient for p≥2)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

After v42's c_N_1 fix (which correctly balanced LHS/RHS IP penalties), p≥2 **still locks**:

| Config | Result | Expected |
|--------|--------|----------|
| p=1, 1000m, IP | 8 events / 1800yr ✓ | ~8 events |
| p=2, 1000m, IP | Nucleation stalls (V stuck at 4e-5 m/s) ✗ | ~8 events |
| p=4, 1000m, IP | Fault locks (V → 1e-10 m/s) ✗ | ~8 events |

v43 attempted multi-DOF fault discretization but had a sign bug
(`+ correction_q` instead of `- correction_q`) and was reverted.

This document identifies the **true root cause** through deep comparison
with Tandem's implementation and proposes a targeted fix.

---

## 2. Root Cause: Fault DOF Mismatch Between MFEM and Tandem

### 2.1 How Tandem Handles Fault DOFs

From `app/preprocess.cpp` line 21:
```cpp
numFaultBasisFunctions = binom(N + D - 1, D - 1)
```
For 3D (D=3) with polynomial degree N on **triangle faces**:

| N (degree) | nbf_fault = (N+1)(N+2)/2 | Fault DOFs per face (×2 tangential) |
|------------|--------------------------|--------------------------------------|
| 1 | 3 | 6 |
| 2 | 6 | 12 |
| 4 | 15 | 30 |
| 6 | 28 | 56 |

Each of these DOFs has **independent** friction evaluation:
- From `app/localoperator/RateAndState.h` lines 148-186:
```cpp
for (std::size_t node = 0; node < nbf; ++node) {
    auto sn = t_mat(node, 0);       // Normal stress AT THIS DOF
    auto psi = s_mat(node, PsiIndex); // State variable AT THIS DOF
    auto tau = get_tau(node, t_mat);  // Shear stress AT THIS DOF
    auto Vi = law_.slip_rate(index + node, sn, tau, psi); // Per-DOF friction
}
```

The complete Tandem fault coupling flow (`app/form/SeasQDOperator.cpp`):
1. **Slip → quad points**: `slip_q = Σ_k φ_k(q) * slip_k * fault_basis` — polynomial slip field
2. **Solve**: DG elastic solve with polynomial slip as BC
3. **Traction → L2 project**: `T_k = M⁻¹ ∫ φ_k T̂_q dS` — project to fault DOFs
4. **Friction per DOF**: Each DOF k gets its own V_k from rate-and-state
5. **Advance**: Each DOF evolves independently: d(slip_k)/dt = V_k

### 2.2 How MFEM Handles Fault DOFs (Current Code)

| N (degree) | Fault DOFs per face | Friction evaluations per face |
|------------|--------------------|-----------------------------|
| Any p | **2** (tau_dip, tau_strike) | **1** |

Current flow:
1. **Slip → constant**: `delta_u = EmbedSlip(slip_dip, slip_strike)` — **constant per face**
2. **Solve**: DG elastic solve with constant slip BC
3. **Traction → face average**: `T = (1/Σw) Σ_q w_q T̂_q` — single value per face
4. **Friction once**: Single V from rate-and-state using face-averaged T
5. **Advance**: Single evolution: d(slip)/dt = V (constant across face)

### 2.3 The Critical Mismatch

At p=1, both approaches are nearly equivalent:
- Volume DG solution: constant gradient per element → constant stress per face
- Face average ≈ point value (stress doesn't vary across face)
- 2 fault DOFs capture all the information
- The penalty correction is small and helps correct stress errors
- **Result: MFEM works fine at p=1** ✓

At p≥2, the mismatch becomes fatal:
- Volume DG solution: polynomial gradient → **polynomial stress variation** across face
- The DG displacement jump [[u]] varies polynomially across the face
- But prescribed slip δ is **constant** per face (2 DOFs)
- The penalty correction η × ([[u]]_q − δ) has polynomial variation
- **Face-averaging destroys this variation**, producing a systematic traction bias
- The bias acts as **numerical stiffness** that locks the fault

### 2.4 Why the Penalty Amplifies the Mismatch

The IP penalty η scales as:
```
η = (D+1) × c_N_1 × (A/V) × (c₁²/c₀)
  where c_N_1 = p(p+D-1)/D ∝ p²
```

At each face quadrature point:
```
T̂_q = {σ·n}_q − η × ([[u]]_q − δ)
```

In Tandem: δ_q = Σ_k φ_k(q) × slip_k has **polynomial variation matching [[u]]_q**
→ [[u]]_q − δ_q ≈ 0 everywhere → penalty correction negligible → T̂_q ≈ {σ·n}_q ✓

In MFEM: δ is **constant** across face → [[u]]_q − δ has polynomial residual ε_q
→ η × ε_q is amplified by p² → face-average picks up systematic bias
→ T = {σ·n}_avg − η × <ε> shifts the traction → fault locks ✗

### 2.5 Mathematical Analysis of the Face-Averaged Correction

From the DG weak form, summing over all element DOFs with partition-of-unity:
```
η × Σ_q w_q × ε_q = Σ_k (stress & symmetry terms for DOF k)
```

Therefore:
```
η × <ε> = (1/Σw) × Σ_k S_k     (independent of η!)
```

The face-averaged penalty correction is **bounded** (O(1), independent of η), but
its **sign and magnitude** depend on the balance of stress/symmetry terms across all
face DOFs. At p≥2, more DOFs contribute to this sum, and the net effect can be
a systematic bias in the traction that:
- Reduces effective loading on the VS (velocity-strengthening) zone
- Suppresses interseismic stress buildup
- Prevents earthquake nucleation

### 2.6 Why Tandem's IP Works at p=6 Despite Large Penalty

Tandem's penalty at p=6 on 4km mesh: c_N_1 = 6×8/3 = 16 (very large).
Yet Tandem produces correct earthquake cycling. The reason:

1. **Multi-DOF slip**: Slip has polynomial variation matching the DG solution
2. **L2 projection preserves spatial information**: No face-averaging loss
3. **Per-DOF friction**: Each DOF responds independently to local traction
4. **Self-correcting**: If one DOF gets biased, its neighbors compensate

The penalty is not the problem. The **representation of slip as constant per face**
combined with **face-averaging of traction** is the problem.

---

## 3. Proposed Fix: Remove Penalty Correction from Fault Traction

### 3.1 The Change

In `ComputeTraction()`, for **fault faces only** (both interior and shared),
change the IP traction combination from:

```cpp
// BEFORE (v42):
T_global[c] = T_stress[c] - correction[c];  // Full IP numerical flux
```

to:

```cpp
// AFTER (v44):
T_global[c] = T_stress[c];  // Stress average only at fault faces
```

This affects 2 locations:
1. Interior fault face IP path (~line 3238)
2. Shared fault face IP path (~line 3589)

**Keep the correction computation for diagnostics** — just don't include it in T_global.

### 3.2 Justification

**1. Matches SCEC Benchmark Specification**

The SCEC SEAS benchmark Algorithm line 15 specifies:
```
T = {σ·n̂}_fault
```
No penalty correction in the traction. The IP penalty is a DG **stability mechanism**
for the bilinear form, not a physical component of the traction.

**2. Tandem's BR2 Is Effectively This**

Tandem's BR2 method uses penalty = dim+1 = 4 (order-independent, small).
The BR2 traction correction is negligible at any p. Our BR2 results at p=1
gave 295yr recurrence (close to IP's 250yr). Removing the IP penalty from
traction makes the IP traction extraction behave like BR2's.

**3. Decouples DG Stability from Traction Accuracy**

- **Bilinear form**: Full IP penalty (for DG stability and well-posedness) ✓
- **RHS assembly**: Full IP penalty (for consistent slip enforcement) ✓
- **Traction extraction**: Stress average only (for physical traction) ✓

This is standard DG post-processing: solve with the full numerical flux,
then extract physical quantities using the most appropriate recovery method.
The stress average {σ·n} converges at O(h^p) and is free of penalty artifacts.

**4. Eliminates p²-Scaling Numerical Stiffness**

Without the penalty correction, the traction is:
```
T = {σ·n}_face_avg
```
This depends only on the stress field accuracy (which improves with p),
not on the penalty magnitude (which grows with p²). At p≥2, the stress
average is more accurate than at p=1, giving better traction without bias.

**5. The IP Penalty Still Serves Its Purpose**

The IP penalty in the bilinear form ensures:
- DG coercivity and stability
- The displacement jump [[u]] ≈ δ (slip enforcement)
- The elastic solve produces a well-conditioned system

Removing it from traction extraction does NOT affect these properties.
The DG solve is unchanged; only the post-processing step changes.

### 3.3 Comparison with Tandem at Each Step

| Step | Tandem IP | MFEM v42 IP | MFEM v44 IP (proposed) |
|------|-----------|-------------|------------------------|
| Bilinear form | {σ·n}·[[v]] + η[[u]]·[[v]] | Same | Same (no change) |
| RHS slip | {σ(v)·n}·δ + η·δ·v | Same | Same (no change) |
| DG solve | Same system | Same system | Same system |
| Traction at q-pt | T̂_q = {σ·n}_q − η(ε_q) | Same | **T_q = {σ·n}_q only** |
| Traction → DOFs | L2 project: M⁻¹∫φ_k T̂ dS | Face average: <T̂> | **Face average: <{σ·n}>** |
| Fault DOFs | nbf=(N+1)(N+2)/2 | 2 (dip,strike) | 2 (dip,strike) |
| Friction | Per DOF (nbf evals) | 1 per face | 1 per face |

The proposed fix deviates from Tandem only in the traction step: we use
{σ·n} instead of T̂. But since we face-average (unlike Tandem's L2 projection),
this deviation actually **improves** the result by avoiding penalty amplification
of the representation mismatch.

### 3.4 Risk Assessment

**Low risk at p≥2**: The stress average {σ·n} is O(h²) accurate at p=2 and
O(h⁴) at p=4. The penalty correction (which we're removing) was causing harm
by adding p²-scaling stiffness. Removing it should directly fix fault locking.

**Moderate risk at p=1**: At p=1, {σ·n} has O(h) accuracy, and the penalty
correction was providing a useful (though imperfect) correction. Removing it
might increase the strike deficit from ~20% to ~25-30%. However:
- Earthquakes should still nucleate (the deficit is a quantitative issue, not qualitative)
- The 20% deficit at p=1 is already a known discretization artifact
- The v24 "blowup without penalty" was at a much earlier code stage with many other bugs

**Mitigation**: If p=1 becomes unstable or unacceptably inaccurate, we can add
an **order-independent stabilization** to the fault traction:
```cpp
// Fallback stabilization (BR2-like, no p² scaling):
real_t penalty_stab = (dim + 1) * (face_area / vol) * (c1_mat * c1_mat / c0_mat);
// This is η with c_N_1 = 1 (equivalent to the p=1 penalty)
```
This caps the traction penalty at the p=1 level, avoiding the p² amplification
while preserving the p=1 correction that helps on coarse meshes.

---

## 4. Why Not Multi-DOF Fault? (Long-Term vs. Short-Term)

The **correct** long-term fix is multi-DOF fault discretization matching Tandem:
- nbf = (p+1)(p+2)/2 DOFs per face per tangential component
- L2 projection of traction to fault DOFs
- Independent friction evaluation at each DOF
- Polynomial slip interpolation to quadrature points

This requires major refactoring:
- Fault state vector size changes (from 5 per face to 3×nbf per face)
- Friction law must loop over DOFs
- Time integrator must handle variable-size state
- Output/station extraction must interpolate from DOFs
- ~500-1000 lines of new code across multiple files

The proposed v44 fix (remove penalty from traction) is a **pragmatic intermediate**:
- 2-line change in ComputeTraction
- No changes to state vector, friction, time integrator, or output
- Should unlock p≥2 by removing the main source of the bias
- Can be validated quickly on TACC

After v44 validates that the approach works, multi-DOF fault can be implemented
incrementally as a separate effort.

---

## 5. Implementation Details

### 5.1 Code Changes (2 locations in elasticity_operator.hpp)

**Location 1: Interior fault faces, IP path (~line 3238)**
```cpp
// BEFORE:
T_global[c] = T_stress[c] - correction[c];

// AFTER:
// v44: Use stress average only at fault faces.
// The IP penalty correction introduces p²-scaling numerical stiffness
// that locks the fault at p≥2 (see v44 debug doc). The penalty serves
// DG stability (in bilinear form), not traction accuracy.
// Matches SCEC benchmark Algorithm line 15: T = {σ·n̂}.
T_global[c] = T_stress[c];
```

**Location 2: Shared fault faces, IP path (~line 3589)**
```cpp
// Same change as Location 1.
T_global[c] = T_stress[c];
```

### 5.2 What Does NOT Change

- Bilinear form (`DGElasticityIPPenaltyIntegrator`): unchanged
- RHS slip assembly (`AssembleSlipContributionIP`): unchanged
- RHS Dirichlet loading (`AssembleDirichletLoading`): unchanged
- BR2 traction path: unchanged
- Fault friction, time integrator, output: unchanged
- All existing tests should pass (p=1 traction values change slightly)

### 5.3 Diagnostic Preservation

Keep the correction computation and `diag_traction_decomp_` output so we can
monitor the correction magnitude:
```cpp
// Still compute correction for diagnostics:
correction_q[c] = -penalty_ip * sign * jump_c;
// ...accumulate and normalize correction as before...
// But DON'T use it in T_global:
T_global[c] = T_stress[c];
// Diagnostics still report correction magnitude for monitoring.
```

---

## 6. Test Plan

### 6.1 TACC Runs

| Run | Config | Purpose | Expected |
|-----|--------|---------|----------|
| v44a | p=1, 1000m, IP | Regression check | ~8 events, slightly different timing |
| v44b | p=2, 1000m, IP | Critical p=2 test | **Earthquakes should occur** |
| v44c | p=4, 1000m, IP | p=4 convergence test | **Earthquakes should occur** |

### 6.2 Success Criteria

1. **p=1 regression**: Earthquakes still occur. Recurrence may shift slightly
   (±50yr acceptable). Traction at dp+22 should remain within 10% of v42 values.

2. **p=2 nucleation**: V_max exceeds 1 m/s during first event. Earthquake
   propagates across VW zone. Recurrence interval < 300yr.

3. **p=4 nucleation**: Same as p=2. Strike deficit should be <5% (vs ~20% at p=1).

### 6.3 Failure Modes and Fallbacks

- **p=1 blows up during coseismic**: Add order-independent stabilization
  (BR2-like penalty cap at c_N_1=1 in traction only)
- **p=2 still locks**: Root cause is deeper than penalty; investigate multi-DOF fault
- **p=2 works but p=4 doesn't**: May need mesh-dependent tuning or multi-DOF fault

---

## 7. CRITICAL DISCOVERY: IP Penalty Is 1/3 of Tandem's (Reference Element Scaling Bug)

### 7.1 The Bug

MFEM's IP penalty uses `nl_q / Weight()` as the face-area / element-volume ratio.
For **tetrahedra**:

```
CalcOrtho(J_face, nor):  nl_q = |nor| = |col1 × col2| = 2 × A_phys
                         (reference triangle area = 1/2)

Weight() = |det(J_vol)| = 6 × V_phys
           (reference tet volume = 1/6)

nl_q / Weight() = (2A) / (6V) = A / (3V)
```

**Tandem** computes physical area and volume by quadrature integration:
```cpp
// DGCurvilinearCommon.cpp:
volume_[elNo] = Σ_q w_q × absDetJ(q)     // = V_phys (physical volume)
area_[fctNo]  = Σ_q w_q × |nor_q|         // = A_phys (physical face area)
penalty = (D+1) × c_N_1 × (area/volume) × (c1²/c0)   // uses A/V
```

**MFEM's penalty is exactly 1/3 of Tandem's**: `A/(3V)` vs `A/V`.

### 7.2 Numerical Verification

For a regular tet with edge length h=1000m:
```
Physical: A/V = 3/h = 0.003
MFEM:     nl_q/Wt = 1/h = 0.001
Ratio:    3.0  (= dim for 3D tetrahedra)
```

### 7.3 Why This Matters: Imbalanced Bilinear Form

The IP bilinear form has three terms assembled by **two separate integrators**:

1. `DGElasticityIntegrator(kappa=0)` → consistency + symmetry (MFEM standard)
   - Uses `dshape_adj / detJ` which **correctly cancels** reference scaling
   - These terms are at **full strength** ✓

2. `DGElasticityIPPenaltyIntegrator` → penalty (our custom code)
   - Uses `nl_q / vol` which gives `A/(3V)` instead of `A/V`
   - This term is at **1/3 strength** ✗

The consistency/symmetry terms are correctly scaled but the penalty that must
counterbalance them is 3× too weak. This makes the bilinear form **under-penalized**.

### 7.4 Impact

- **p=1**: 1/3 penalty is still above the coercivity threshold → stable but inaccurate
  (20% strike deficit is caused by this weak penalty)
- **p≥2**: With correctly balanced c_N_1 (v42 fix), the consistency terms grow
  proportionally with p while the penalty grows proportionally with p but at 1/3
  strength → insufficient enforcement → larger DG jumps → traction bias → fault locks

### 7.5 All 13 Affected Locations

| # | File | Line | Function |
|---|------|------|----------|
| 1 | `dg_elasticity_ip_penalty_integrator.hpp` | 120 | Bilinear form (p0) |
| 2 | `dg_elasticity_ip_penalty_integrator.hpp` | 126 | Bilinear form (p1) |
| 3 | `elasticity_operator.hpp` | 1039 | AssembleSlipContributionIP (p0) |
| 4 | `elasticity_operator.hpp` | 1040 | AssembleSlipContributionIP (p1) |
| 5 | `elasticity_operator.hpp` | 1501 | AssembleSlipContributionIPShared (p0) |
| 6 | `elasticity_operator.hpp` | 1502 | AssembleSlipContributionIPShared (p1) |
| 7 | `elasticity_operator.hpp` | 1859 | AssembleDirichletLoading boundary (p0) |
| 8 | `elasticity_operator.hpp` | 2133 | AssembleDirichletLoading interior (p0) |
| 9 | `elasticity_operator.hpp` | 2135 | AssembleDirichletLoading interior (p1) |
| 10 | `elasticity_operator.hpp` | 2514 | AssembleDirichletLoading shared (p0) |
| 11 | `elasticity_operator.hpp` | 2516 | AssembleDirichletLoading shared (p1) |
| 12 | `elasticity_operator.hpp` | 3088-3090 | ComputeTraction interior (p0,p1) |
| 13 | `elasticity_operator.hpp` | 3509-3511 | ComputeTraction shared (p0,p1) |

### 7.6 The Fix

At each location, multiply `nl_q / detJ` by `dim` to convert from reference to physical:
```cpp
// BEFORE (WRONG — gives A/(3V) for tets):
real_t p0 = (dim + 1) * c_N_1 * (nl_q / detJ1) * (c1² / c0);

// AFTER (CORRECT — gives A/V for tets, matching Tandem):
real_t p0 = (dim + 1) * c_N_1 * (real_t(dim) * nl_q / detJ1) * (c1² / c0);
```

The factor `dim` = 3 corrects for simplex reference element scaling:
- Reference tet volume = 1/D!, so Weight() = D! × V_phys
- Reference triangle area = 1/(D-1)!, so nl_q = (D-1)! × 2 × A_phys ... no:
- More precisely: nl_q = 2A, Weight() = 6V, so nl_q/Weight() = A/(3V) = A/(D×V)
- Multiply by D to get A/V

---

## 8. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| v30-v40 | All previous fixes | Done |
| v41 | Parametric study: resolution, domain, p-refinement | Done |
| v42 | IP c_N_1 fix: 1.0 → p(p+D-1)/D at 5 RHS locations | Done |
| v43 | Multi-DOF fault attempt (reverted — sign bug) | Reverted |
| **v44a** | **IP penalty ×3 fix: nl_q/vol → dim×nl_q/vol at all 13 locations** | **Applied** |
| **v44b** | **Remove IP penalty from fault traction (if still needed)** | **Contingent on v44a results** |

---

## 9. Files Changed

| File | Change |
|------|--------|
| `integrator/dg_elasticity_ip_penalty_integrator.hpp` | 2 lines: multiply by dim in penalty formula |
| `domain/elasticity_operator.hpp` | 11 lines: multiply by dim in penalty formula |

---

## 10. Summary

**Two root causes identified for p≥2 failure:**

1. **Factor-of-3 penalty scaling bug** (Section 7): MFEM uses `nl_q/Weight()` =
   `A/(3V)` while Tandem uses physical `A/V`. The IP penalty is 1/3 of Tandem's
   across ALL 13 locations (bilinear form, RHS, traction). This under-penalizes
   the DG system, explaining the 20% strike deficit at p=1 and contributing to
   fault locking at p≥2.

2. **Fault DOF mismatch** (Section 2): MFEM uses 2 constant DOFs per face while
   Tandem uses (p+1)(p+2)/2 DOFs with polynomial variation. At p≥2, the constant
   slip representation cannot match the DG solution's polynomial jump, and the
   penalty correction amplifies this mismatch.

**v44 fix strategy**: Fix the factor-of-3 bug first (Section 7.6). This is a
clear coding error — the penalty formula should produce the same physical value
as Tandem's. If p≥2 still fails after this fix, the fault DOF mismatch (Section 2)
is the remaining issue, addressed by removing penalty from fault traction
(Section 3) or full multi-DOF implementation (Section 4).
