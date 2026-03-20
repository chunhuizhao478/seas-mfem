# BP5 Debug v47: IP Penalty ×3 Correction — Reference Element Scaling Fix

**Date**: 2026-03-20
**Status**: Analysis complete, implementing fix
**Previous**: v46 (Tandem initialization defaults + multi-DOF slip indexing fix)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

v46 results with Tandem initialization defaults show:
- **p=1**: Closeup match to benchmark is reasonable, but later events diverge
  (timing drift, slip offset, stress divergence by event 3-4)
- **p=2**: Closeup match is good initially, but diverges MORE than p=1 at later
  times (larger timing drift, bigger dip-slip deviation, starting earlier)

Both orders show systematic drift from Tandem reference. p=2 should be MORE
accurate than p=1 (higher spatial resolution), but shows WORSE long-term behavior.
This points to a systematic error that scales with polynomial order.

---

## 2. Root Cause: IP Penalty at 1/3 of Tandem's Value

### 2.1 The Mathematical Proof

The IP penalty formula in both codes is:
```
p(K) = (D+1) × c_N_1 × (A_face / V_elem) × (c₁² / c₀)
```

The difference is in how `A_face / V_elem` is computed:

**Tandem** (Elasticity.cpp:278-290, DGCurvilinearCommon.cpp:55-94):
```
area_[fctNo] = Σ_q w_q × |normal_q|     (integrated physical face area)
volume_[elNo] = Σ_q w_q × |det(J_q)|    (integrated physical element volume)
A/V = area_[fctNo] / volume_[elNo]       → physical ratio
```

**MFEM** (elasticity_operator.hpp:1077, ip_penalty_integrator.hpp:120):
```
nl_q = CalcOrtho(Jacobian()).Norml2()    (face Jacobian determinant at quad point)
vol  = Trans.Elem1->Weight()             (element Jacobian determinant)
A/V  = nl_q / vol                        → reference-scaled ratio
```

For linear (flat-faced) tetrahedra with constant Jacobians:
```
nl_q = |J_face| = 2 × A_phys / (ref triangle area) = 2 × A_phys / (1/2) = ...
```

More precisely:
- Physical face area: `A_phys = nl_q × (reference face area) = nl_q × 1/(D-1)!`
  → For D=3: `A_phys = nl_q × 1/2`
- Physical element volume: `V_phys = Weight() × (reference tet volume) = Weight() × 1/D!`
  → For D=3: `V_phys = Weight() × 1/6`

Therefore:
```
A_phys / V_phys = (nl_q × 1/2) / (Weight() × 1/6)
                = (nl_q / Weight()) × (1/2) / (1/6)
                = (nl_q / Weight()) × 3
                = (nl_q / Weight()) × D
```

**MFEM's `nl_q / Weight()` = (1/D) × (A_phys / V_phys) = (1/3) × Tandem's value.**

The correction factor is exactly D = dim = 3 for 3D.

### 2.2 Why This Factor Arises

The reference simplex measures are:
```
ref_triangle_area = 1/(D-1)! = 1/2!  = 1/2    (for D=3 faces)
ref_tet_volume    = 1/D!     = 1/3!  = 1/6    (for D=3 elements)
```

MFEM's `CalcOrtho` returns the face Jacobian determinant (maps from ref → physical),
and `Weight()` returns the element Jacobian determinant. The ratio of these determinants
maps ref_face/ref_vol, NOT physical area/volume:
```
nl_q / Weight() = (A_phys / ref_face) / (V_phys / ref_vol)
                = (A_phys / V_phys) × (ref_vol / ref_face)
                = (A_phys / V_phys) × (1/6) / (1/2)
                = (A_phys / V_phys) / 3
```

General D: `nl_q / Weight() = (A_phys / V_phys) × [1/D! / 1/(D-1)!] = (A_phys / V_phys) / D`

### 2.3 Why It Affects p=2 More Than p=1

Both p=1 and p=2 have the same 1/3 penalty deficit. But the impact differs:

1. **Coercivity threshold scales as p²**: `c_N_1 = p(p+D-1)/D` gives:
   - p=1: c_N_1 = 1.0
   - p=2: c_N_1 = 8/3 ≈ 2.67

   The consistency terms grow proportionally with c_N_1 while the penalty is
   c_N_1/3. At p=2, the ratio of consistency-to-penalty is further from
   the designed balance point.

2. **Multi-DOF amplification**: At p=2, each face has 6 DOFs with polynomial
   variation. The weak penalty allows larger DG jumps, which create polynomial
   stress variations across the face. These interact with the multi-DOF
   friction evaluation, causing systematic traction bias.

3. **Error accumulation**: Both orders accumulate timing errors over earthquake
   cycles. p=2 starts with slightly less spatial error (higher resolution) but
   the weaker-than-designed penalty causes larger penalty correction errors,
   which compound faster over time.

---

## 3. Why the v45 Revert Was Wrong

### 3.1 The v45 "Balance" Argument — Refuted

v45 (bp5_debug_v45.md, Section 3) argued:

> "Applying ×3 to only the penalty shifts the consistency-penalty balance,
> producing a different (not better) DG solution."

This argument has two flaws:

**Flaw 1: Consistency terms ARE already at physical strength.**
v45's own table (Section 3.2) shows:
| Term | v42 (no ×3) | v44 (with ×3) | Tandem |
|------|-------------|---------------|--------|
| Consistency | 1.0× | 1.0× | 1.0× |
| Symmetry | 1.0× | 1.0× | 1.0× |
| Penalty | 0.33× | 1.0× | 1.0× |

The v44 column (with ×3) matches Tandem EXACTLY. The v42 column has an
IMBALANCED system (consistency at 1.0× but penalty at 0.33×). Claiming
that v44 creates an "imbalance" is backwards — v42 is the imbalanced one.

**Flaw 2: "Both are valid SIPG" misses the point.**
Both 1/3 penalty and full penalty are technically stable SIPG if above
the coercivity threshold. But we want to MATCH Tandem's solution, not
just any valid SIPG solution. Tandem uses full penalty, so we must too.

### 3.2 The v44e "Results Got Worse" — Misleading Context

v45 Section 6 states that v44e (×3 fix on clean v42 baseline) showed
"p=1 earthquake timing shifted significantly" and concluded ×3 was harmful.

This conclusion was flawed because:

1. **v44e had the slip indexing bug** (only discovered in v46). While this
   bug doesn't affect p=1 (nbf=1, indices coincide), it means the entire
   v44e codebase was untested at p≥2, making it impossible to evaluate
   the ×3 fix where it matters most.

2. **Comparison was against Tandem p=4 reference**, not Tandem p=1. A p=1
   MFEM solution should NOT match a p=4 reference — the discretization
   error is ~20%. The "better match" of v42 was **compensating error**:
   the 1/3 penalty error happened to offset the spatial discretization
   error, giving a fortuitously good match at p=1.

3. **"Shifted" ≠ "worse"**. The v44e results were different from v42, as
   expected (different penalty = different solution). The question is which
   converges faster as p increases — and that was never tested because
   the slip indexing bug made p≥2 crash.

### 3.3 Why v47 Can Properly Evaluate the ×3 Fix

v47 has all prerequisites that v44 lacked:
1. ✅ Slip indexing bug fixed (v46)
2. ✅ Multi-DOF fault working at p=2 (v46)
3. ✅ Tandem initialization defaults (v46 Phase 2)
4. ✅ Clean penalty fix without other confounding changes

The test: p=1 and p=2 with ×3 penalty should BOTH move closer to Tandem.
If p=2 matches better than p=1, that confirms correct convergence behavior.

---

## 4. Detailed Code Comparison: Tandem vs MFEM

### 4.1 Tandem Penalty (Physical A/V)

**File**: `tandem/app/localoperator/Elasticity.cpp:278-290`
```cpp
void Elasticity::prepare_penalty(std::size_t fctNo, FacetInfo const& info,
                                 LinearAllocator<double>&) {
    auto const p = [&](int side) {
        const auto [c0, c1] = stiffness_tensor_bounds(info.up[side]);
        constexpr double c_N_1 = InverseInequality<Dim>::trace_constant(PolynomialDegree - 1);
        return (Dim + 1) * c_N_1 * (area_[fctNo] / volume_[info.up[side]]) * (c1 * c1 / c0);
    };
    if (info.up[0] != info.up[1]) {
        penalty_[fctNo] = (p(0) + p(1)) / 4.0;
    } else {
        penalty_[fctNo] = p(0);
    }
}
```

Where:
- `area_[fctNo]` = integrated face area (DGCurvilinearCommon.cpp:90-94)
- `volume_[elNo]` = integrated element volume (DGCurvilinearCommon.cpp:55-59)
- These are PHYSICAL values (not reference-scaled)

**Assembly** (Elasticity.cpp:422-443, kernels/elasticity.py:141-145):
```python
# penalty term in assembleSurface:
c20 * w_q * E_q[x] * L_q[y]
# where L_q (IP lift) = E_q * delta * nl_q  (line 122-123)
# → penalty * Σ_q w_q * nl_q * φ_i * φ_j = penalty * ∫_face φ_i φ_j ds
```

### 4.2 MFEM Penalty (Reference-Scaled A/(3V))

**File**: `seas/integrator/dg_elasticity_ip_penalty_integrator.hpp:116-132`
```cpp
real_t vol1 = Trans.Elem1->Weight();   // = |det(J_elem)| = D! × V_phys = 6V
real_t p0 = (dim_ + 1) * c_N_1 * (nl_q / vol1) * (c1 * c1 / c0);
//                                      ^^^^^^^^
//                                      nl_q / Weight() = A/(3V), NOT A/V
```

**After fix** (v47):
```cpp
real_t vol1 = Trans.Elem1->Weight();
// Physical A/V ratio: correct for ref element measure ratio (1/(D-1)!) / (1/D!) = D
real_t p0 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / vol1) * (c1 * c1 / c0);
//                                  ^^^^^^^^^^^^
//                                  multiply by dim = 3 to get physical A/V
```

### 4.3 All 15 Locations Requiring the Fix

In `integrator/dg_elasticity_ip_penalty_integrator.hpp`:
| Line | Context |
|------|---------|
| 120 | Bilinear form: p0 (elem1) |
| 126 | Bilinear form: p1 (elem2) |

In `domain/elasticity_operator.hpp`:
| Line | Context |
|------|---------|
| 1077 | AssembleSlipContributionIP: p0 (fault interior faces) |
| 1078 | AssembleSlipContributionIP: p1 (fault interior faces) |
| 1523 | AssembleFarFieldDirichletRHS_IP: p0 (skeleton faces) |
| 1524 | AssembleFarFieldDirichletRHS_IP: p1 (skeleton faces) |
| 1878 | AssembleBoundaryDirichletRHS_IP: p0 (boundary faces) |
| 2152 | AssembleSkeletonFarFieldRHS_IP: p0 (skeleton Dirichlet) |
| 2154 | AssembleSkeletonFarFieldRHS_IP: p1 (skeleton Dirichlet) |
| 2533 | AssembleSkeletonFarFieldRHS_IP_v2: p0 |
| 2535 | AssembleSkeletonFarFieldRHS_IP_v2: p1 |
| 3096 | ComputeTraction (interior): p0 |
| 3098 | ComputeTraction (interior): p1 |
| 3544 | ComputeTraction (shared): p0 |
| 3546 | ComputeTraction (shared): p1 |

### 4.4 Comparison: v44 ×3 Fix vs v47 ×3 Fix

| Aspect | v44 | v47 |
|--------|-----|-----|
| Slip indexing | Bug present (fi instead of fi*nbf+kk) | Fixed (v46) |
| Initialization | SCEC defaults (delta_tau=138kPa, V_nuc=0.03) | Tandem defaults (no delta_tau, V_nuc=0.01) |
| Multi-DOF | Partial, untested at p=2 (MUMPS issues) | Working, tested at p=1 and p=2 |
| Penalty fix locations | 13 locations | 15 locations (includes both ComputeTraction variants) |
| Other changes | Also removed penalty from traction | Pure penalty scaling only |

---

## 5. Tandem Inverse Inequality Reference

**File**: `tandem/src/form/InverseInequality.h:27-28`
```cpp
constexpr static double trace_constant(unsigned N) {
    return (N + 1) * (N + D) / static_cast<double>(D);
}
```

Called with `PolynomialDegree - 1`, so for polynomial degree p:
```
c_N_1 = trace_constant(p - 1) = p × (p + D - 1) / D
```

For D=3:
| p | c_N_1 |
|---|-------|
| 1 | 1 × 3 / 3 = 1.0 |
| 2 | 2 × 4 / 3 = 2.667 |
| 4 | 4 × 6 / 3 = 8.0 |
| 6 | 6 × 8 / 3 = 16.0 |

Our code matches exactly (elasticity_operator.hpp:1076):
```cpp
real_t c_N_1 = order_ * (order_ + dim - 1.0) / dim;
```

---

## 6. Fix Implementation

### 6.1 The Fix

At all 15 locations, multiply `nl_q / detJ` (or `face_area / vol`) by `dim`:

```cpp
// BEFORE (1/3 of physical):
real_t p0 = (dim + 1) * c_N_1 * (nl_q / detJ1) * (c1_mat * c1_mat / c0_mat);

// AFTER (correct physical A/V):
real_t p0 = (dim + 1) * c_N_1 * (real_t(dim) * nl_q / detJ1) * (c1_mat * c1_mat / c0_mat);
```

### 6.2 Why This Is Self-Consistent

The penalty must be consistent between:
1. **Bilinear form** (stiffness matrix K): determines DG solution
2. **RHS assembly** (slip BC, Dirichlet BC): provides load vector
3. **Traction recovery** (ComputeTraction): post-processes solution

All 15 locations use the same formula, so applying ×dim everywhere maintains
the LHS-RHS-traction consistency that v42 established.

The consistency/symmetry terms (MFEM's `DGElasticityIntegrator`) use
`adj(J)/det(J) × nor` which correctly evaluates to the physical integral
(the det(J) factors cancel). These terms need NO modification.

### 6.3 Expected Impact

With correct penalty:
1. **DG jumps decrease**: stronger enforcement → solution jumps [[u]] are smaller
2. **Traction accuracy improves**: penalty correction is smaller relative to stress
3. **p-convergence improves**: penalty now scales correctly with p, matching Tandem's
   designed coercivity margin
4. **p=1 results will shift**: the compensating error is removed, exposing true
   discretization error at p=1. This is EXPECTED and CORRECT behavior.
5. **p=2 should improve significantly**: the multi-DOF fix (v46) + correct penalty
   (v47) should produce the designed convergence behavior

---

## 7. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v30-v42 | All previous fixes through c_N_1 | Done |
| v44 | Penalty ×3 (first attempt) — reverted in v45 due to flawed analysis | Reverted |
| v45 | Analysis claiming ×3 is wrong + slip indexing bug discovery | Done |
| v46 | Multi-DOF slip indexing fix + Tandem initialization defaults | Done |
| **v47** | **Penalty ×3 re-applied with corrected analysis and clean codebase** | **In Progress** |
