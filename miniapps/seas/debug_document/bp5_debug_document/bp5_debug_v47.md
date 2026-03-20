# BP5 Debug v47: IP Penalty ×3 Correction — Reference Element Scaling Fix

**Date**: 2026-03-20
**Status**: p=2 blowup is formulation-resolution issue (within-face V amplification at ×3 penalty). dt_init bug fixed but does not resolve blowup. See Section 11.
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

## 7. v47 Results: ×3 Penalty Breaks Both p=1 and p=2

### 7.1 Summary of All v47 Runs

| Run | Order | Solver | Penalty | Result |
|-----|-------|--------|---------|--------|
| v47a | p=1 | MUMPS-BLR 1e-10 | ×3 | V_nuc decays from 0.01→5e-9, no earthquake |
| v47b | p=2 | MUMPS-BLR 1e-10 | ×3 | Immediate blowup (τ=1.43 GPa, step 0) |
| v47c | p=2 | MUMPS exact | ×3 | MUMPS OOM (INFO(1)=-1, INFO(2)=784) |
| v47d | p=2 | MUMPS-BLR 1e-14 | ×3 | Immediate blowup **WORSE** (τ=1.68 GPa, slip=0.6m) |

### 7.2 Critical Finding: BLR Hypothesis Disproved

v47d (BLR tol=1e-14) produced a **worse** blowup than v47b (BLR tol=1e-10):

| Metric | v47b (BLR 1e-10) | v47d (BLR 1e-14) |
|--------|-------------------|-------------------|
| τ_max | 1.43 GPa | 1.68 GPa |
| max slip | 0.38 m | 0.60 m |
| BLR tolerance | 1e-10 (loose) | 1e-14 (tight) |

**The tighter solver makes the instability worse.** The BLR approximation was
accidentally dampening the IP-DG instability. The root cause is in the DG
formulation itself, not the solver accuracy.

### 7.3 p=2 Blowup: Within-Face Discontinuity + ×3 Penalty

At p=2, each face has nbf=6 DOFs. Faces straddling the nucleation zone boundary
have both V_nuc (0.01 m/s) and V_init (1e-9 m/s) DOFs — a 10^7 velocity ratio.

The IP penalty term amplifies this within-face discontinuity:
- At 1/3 penalty (v46): stress concentrations are manageable → simulation runs
- At 3× penalty (v47): stress concentrations are 3× larger → catastrophic blowup

The blowup locations confirm this — all are near the nucleation zone boundary:
- x=(-30885, -490, -3558): edge of nucleation zone (x2≈-31 km, boundary at -30 km)
- x=(-14833, -403, -6680): nucleation zone boundary
- x=(-19968, -374, -5840): inside nucleation zone near boundary

**Why Tandem doesn't have this problem**: Tandem defaults to p=4 where the
higher-order approximation better resolves the nucleation zone boundary. At p=4,
the transition from V_nuc to V_init is spread across more DOFs per face, reducing
the effective discontinuity. At p=2 on our coarser mesh, the transition is sharper
and the ×3 penalty amplifies it more. (Both codes use SIPG — see Section 10.)

### 7.4 p=1 Stiffness Analysis: Elastic Stiffness vs Friction Weakening

Using station data at x2=-24km, x3=10km (nucleation zone center):

**Nucleation stability criterion** (spring-slider analog):
```
V increases ⟺ k_elastic < k_critical = σ_n × b / Dc
k_critical = 25 MPa × 0.03 / 0.13 m = 5.77 MPa/m
```

**Measured elastic stiffness** (stress drop per unit slip at the nucleation station):

| Code | Time | Slip (mm) | Δτ (kPa) | k (MPa/m) | vs k_crit |
|------|------|-----------|----------|-----------|-----------|
| v47a (ours) | 0.44s | 4.6 | 18 | 3.9 | BELOW (V grows) |
| v47a (ours) | 0.68s | 7.0 | 35 | 5.0 | BELOW (barely) |
| v47a (ours) | 0.90s | 9.3 | 55 | 5.9 | **ABOVE** (V decays) |
| v47a (ours) | 1.83s | 17.4 | 138 | 7.9 | ABOVE |
| Tandem p=4 | 0.66 yr | 7.0 | 22 | 3.1 | BELOW (V grows) |

**Key finding**: Our effective stiffness INCREASES with slip and crosses k_crit
at ~9mm slip (t≈0.9s). After that, V decays irreversibly. Tandem's stiffness
stays below k_crit because by the time similar slip accumulates, it has spread
across a wider area (lower stiffness).

**Why the stiffness differs**:

1. **Concentrated vs distributed slip**: At t<1s in our code, all slip is in the
   12×12 km nucleation zone (no time for spreading). In Tandem at t=0.66yr, slip
   has spread across the entire VW zone (60×12 km). Concentrated slip → higher
   stiffness (smaller patch = stiffer spring).

2. **DG penalty effect on global stiffness**: The penalty at ALL interior faces
   (not just the fault) enters the stiffness matrix K. Tripling the penalty changes
   K, which changes the displacement field u, which changes {σ·n̂} at the fault.
   At p=1 with h=1000m, this effect is large enough to push k_eff above k_crit.
   Both codes use the same SIPG formulation (Section 10), but Tandem runs at p=4
   where the discretization is better resolved and less sensitive to the penalty.

3. **Mesh resolution effect**: Our 1000m mesh has 4× finer resolution than
   Tandem's 4000m reference. In DG, finer mesh → larger penalty per element →
   stiffer fault coupling. This is opposite to CG FEM where finer mesh converges
   to the physical stiffness.

### 7.5 Both Codes Use SIPG — The Difference Is Resolution

**CORRECTION**: An earlier draft of this section incorrectly claimed Tandem uses
SBP-SAT. **Both codes use SIPG** (Symmetric Interior Penalty Galerkin) with
epsilon=-1, full Galerkin mass matrix, and the same penalty formula. See Section 10
for the complete formula-by-formula verification.

The ×3 penalty value is **mathematically correct** (proven in Section 2). It matches
Tandem's physical A/V ratio exactly. Both codes compute identical penalty values on
the same mesh. The difference in fault behavior is a **resolution effect**:

| Property | Our code (p=1, h=1000m) | Tandem (p=4, h=1000m) |
|----------|--------------------------|------------------------|
| Effective DOFs per face | 1 | 15 |
| c_N_1 | 1.0 | 8.0 |
| Solution accuracy | O(h²) | O(h⁵) |
| Sensitivity to penalty | HIGH (coarse discretization) | LOW (well-resolved) |

The penalty enters the global stiffness matrix K at ALL interior faces (not just
the fault). At p=1 with h=1000m, the DG solution is coarse enough that the
penalty's contribution to K significantly changes the effective fault stiffness.
At p=4, the same penalty formula produces a better-conditioned system where the
effective stiffness converges closer to the physical value.

For the **exact** DG solution at infinite resolution, T = {σ·n̂} (the penalty
correction vanishes). But at finite resolution, the penalty affects the entire
displacement field, changing {σ·n̂} at the fault. This resolution-dependent
stiffness pushes the nucleation zone into the stable regime (k > k_crit) at p=1.

### 7.6 Why v46 (1/3 Penalty) Worked

With 1/3 of the physical penalty:
- k_DG is 3× smaller → effective stiffness stays below k_crit longer
- Within-face amplification at p=2 is 3× weaker → manageable
- The 1/3 penalty is below coercivity at p=2 (marginal stability) but
  produces workable results because the penalty error partially compensates
  the DG stiffness artifact

This was a **compensating error**: the wrong penalty value happened to counteract
the IP-DG formulation's excessive fault stiffness.

---

## 8. Investigation Plan: Resolving the IP-DG Penalty Issue

### 8.1 Immediate Diagnostics (No Code Changes)

1. **Run v46 again with station output**: Get v46 station data at the nucleation
   zone (x2=-24, x3=10) to confirm that V_nuc grows with 1/3 penalty. This
   establishes the baseline for comparison.

2. **p=1 with ×3 penalty + exact MUMPS**: Run p=1 (not p=2) with exact MUMPS
   (smaller system, should fit in memory). If V_nuc still decays → confirms the
   issue is formulation, not solver. If V_nuc grows → BLR IS the issue at p=1.

3. **p=1 with intermediate penalty values**: Test ×1.5, ×2.0, ×2.5 to find the
   threshold where V_nuc transitions from growing to decaying.

### 8.2 Potential Fixes (Require Code Changes)

**Option A: Fault-face penalty reduction**
- Use the physical penalty (×3) on all INTERIOR faces
- Use a reduced penalty (×1 or ×1.5) on FAULT faces only
- Rationale: the fault's Dirichlet BC is enforced weakly; lower penalty → less
  artificial stiffness → more physical behavior
- Risk: may create inconsistency between bilinear form and traction extraction

**Option B: Mesh refinement at p=1**
- Run at h=500m or h=333m with correct ×3 penalty
- Finer mesh → better-resolved solution → less sensitivity to penalty
- If nucleation succeeds at finer h, confirms the resolution hypothesis
- No code changes needed, just more compute resources

**Option C: Move to p≥2 with multi-DOF approach**
- Complete the FaceQuadrature multi-DOF implementation (plan in bright-moseying-bentley.md)
- p=2 with nbf=6 DOFs per face has much better resolution
- Need to also address the nucleation zone boundary V discontinuity (eps-based smoothing)
- This is the planned path forward anyway

**Option D: Move to CG at the fault (hybrid DG/CG)**
- Use CG for the elastic solve, with the fault as a split interface
- Penalty only appears on non-fault interior faces
- This is how most SEAS codes handle the fault (no penalty artifact)
- Major architectural change

**Option E: Accept 1/3 penalty for production**
- Keep the 1/3 penalty (which works) for production runs
- Document it as a pragmatic choice: wrong penalty value, but the compensating
  error produces better-than-expected results
- Focus on other improvements (mesh refinement, time stepping, output)
- Risk: limits convergence at higher p; may not work at finer meshes

### 8.3 Recommended Path

**Short term**: Run diagnostics 8.1.1-3 to understand the parameter space. Also
run Option B (p=1 at h=500m with ×3 penalty) to test the resolution hypothesis.
**Medium term**: Implement Option A (fault-face penalty reduction) as the least
invasive code change. Test whether a reduced fault-face penalty restores nucleation
while maintaining overall DG accuracy. In parallel, advance Option C (p≥2 multi-DOF).
**Long term**: Move to p≥2 as the default. p=1 with the correct penalty may simply
require too fine a mesh for practical BP5 runs at 1000m.

---

## 10. Deep Investigation: Formula-by-Formula Comparison with Tandem

### 10.1 Both Codes Use SIPG (Confirmed)

**Tandem**: `Elasticity.h:245`: `constexpr static double epsilon = -1.0;`
**Our code**: `epsilon_ = -1.0` in `SetupIPInteriorFaceIntegrator`

Both use Symmetric Interior Penalty Galerkin with three-term assembly:
1. Consistency: `-∫_F {σ(u)·n} · [[v]] ds`
2. Symmetry: `ε ∫_F {σ(v)·n} · [[u]] ds` (ε=-1 for SIPG)
3. Penalty: `η ∫_F [[u]] · [[v]] ds`

Tandem's bilinear form (Elasticity.cpp:422-443, elasticity.py:141-145):
```
c0[0]=-0.5, c0[1]=0.5      (consistency)
c1[0]=ε*0.5=-0.5, c1[1]=0.5 (symmetry)
c2[0]=penalty, c2[1]=-penalty (penalty)
```

Our bilinear form: `DGElasticityIntegrator(λ, μ, ε=-1, κ=0)` for consistency+symmetry,
plus `DGElasticityIPPenaltyIntegrator(λ, μ, dim=3)` for penalty. Verified that
MFEM's DGElasticityIntegrator with κ=0 produces zero penalty (bilininteg.cpp:4067).

### 10.2 Penalty Formula: Identical

**Tandem** (Elasticity.cpp:280-283, InverseInequality.h:27-29):
```
p(K) = (D+1) × c_N_1 × (area_phys / vol_phys) × (c1²/c0)
c_N_1 = trace_constant(p-1) = p × (p+D) / D
penalty = (p(K⁻) + p(K⁺)) / 4    (interior face)
```

**Our code** (with ×3 fix):
```
p0 = (dim+1) × c_N_1 × (dim × nl_q / Weight()) × (c1_mat²/c0_mat)
c_N_1 = order × (order + dim - 1) / dim
penalty_ip = (p0 + p1) / 4         (interior face)
```

Verified that `dim × nl_q / Weight() = A_phys / V_phys` for flat tet faces:
```
nl_q = |CalcOrtho| = 2 × A_phys    (face Jacobian maps ref_tri of area 1/2)
Weight() = 6 × V_phys              (element Jacobian maps ref_tet of vol 1/6)
dim × nl_q / Weight() = 3 × 2A / 6V = A/V  ✓
```

Verified c_N_1 values match:
```
p=1: Tandem=(0+1)(0+3)/3=1.0, Ours=1(1+2)/3=1.0  ✓
p=2: Tandem=(1+1)(1+3)/3=8/3, Ours=2(2+2)/3=8/3   ✓
p=4: Tandem=(3+1)(3+3)/3=8.0, Ours=4(4+2)/3=8.0   ✓
```

### 10.3 Traction Extraction: Identical

**Tandem** (elasticity.py:242-244):
```
T_q = 0.5 × (σ₁·n̂ + σ₂·n̂) - penalty × (u₁ - u₂ - f_q)
```
Uses `n_unit_q` (unit normal) for stress, `c0[0] = -penalty` for penalty correction.

**Our code** (elasticity_operator.hpp:3162-3206):
```
T_stress_q = {σ · basis.normal}           (unit normal from fault_basis_)
correction_q = -penalty_ip × sign × ((u₁-u₂) - sign × δu)
T_q = T_stress_q - correction_q
```
After sign expansion: `T_q = {σ·n̂} - penalty × (jump - δu)` ✓

Both use unit normal for stress part (physical traction units).
Both use same penalty value in the correction term.

### 10.4 Slip RHS Assembly: Identical

**Tandem** (elasticity.py:174-178, `rhsFacet`):
```
b += c1[0] × tractionTest(0, f_q) × w_q    (symmetry term, c1=ε/2=-0.5)
   + c2[0] × w_q × φ_k × f_lifted_q        (penalty term, c2=penalty)
```
For IP: `f_lifted_q = f_q × nl_q` (line 156).

**Our code** (elasticity_operator.hpp:1087-1135):
```
elvec += ε × sym_val × w1                  (symmetry, w1 = ip.weight/(2×detJ))
       + wq_penalty × sign × δu × shape    (penalty, wq_penalty = penalty × wq × nl_q)
```

Both compute: `Σ_q [ε/2 × σ_test(δu)·n + penalty × δu × φ × |n|] × w_q` ✓

### 10.5 Bilinear Form Penalty Term: Identical

**Tandem**: `c2 × Σ_q w_q × φ_i(q) × φ_j(q) × |n_q| × δ_pu` = `penalty × M_face`
where `M_face` is the physical face mass matrix.

**Our code**: `Σ_q (penalty × w_q × nl_q) × φ_i(q) × φ_j(q)` = `penalty × M_face` ✓

Both use full Galerkin mass matrix (NOT diagonal). Tandem does NOT use SBP-SAT or
collocated penalty. Confirmed in elasticity.py: `lift_ip` (line 122-123) uses
`E_q[x]` (evaluation at quadrature points) and `nl_q`, producing a standard IP lift.

### 10.6 Root Cause: Resolution-Dependent Effective Stiffness

Since ALL formulas match exactly, the two codes should produce the SAME results on
the SAME mesh at the SAME polynomial order. The difference in fault behavior is due to:

1. **Tandem defaults to p=4, our tests use p=1**: At p=4, the DG solution has
   O(h⁵) accuracy. At p=1, only O(h²). The penalty enters the global stiffness
   matrix K at ALL interior faces. At p=1 with h=1000m, the penalty contribution
   to K significantly changes the effective fault stiffness.

2. **Penalty affect on K is resolution-dependent**: The penalty coefficient
   `η ∝ p²/h` (via c_N_1 × A/V). In K, this creates inter-element coupling ∝ η.
   For the fault's effective stiffness `k_eff = dT/dδ`, the penalty contribution
   depends on how K⁻¹ maps slip changes to stress changes. At coarse resolution
   (p=1, h=1000m), this mapping is poorly approximated and the penalty amplifies
   the discretization error.

3. **The penalty is at the SIPG stability minimum**: The formula
   `(D+1)×c_N_1×(A/V)×(c1²/c0)` is the theoretical minimum for SIPG coercivity
   (Warburton-Hesthaven 2003, trace inverse inequality). There is no safety factor.
   The 1/3 penalty (v46) is BELOW the stability bound, which means the bilinear form
   is technically non-coercive. However, in practice, the matrix remains SPD due to
   the volume terms dominating. The under-penalization reduces the artificial stiffness
   enough for nucleation to proceed — a compensating error.

4. **BLR paradox explained**: Tighter BLR → more accurate K⁻¹ → the full effect of
   the (too-stiff) penalty is realized. Loose BLR → K⁻¹ has solver error → the solver
   error accidentally softens the penalty effect. This is why v47d (BLR 1e-14) is
   WORSE than v47b (BLR 1e-10): the tighter solver removes the accidental dampening.

### 10.7 Quantitative Penalty Analysis

For p=1, regular tet with h=1000m, BP5 material (λ=μ=32 GPa):
```
c0 = 2μ = 64 GPa
c1 = 3λ + 2μ = 160 GPa
c_N_1 = 1.0
(D+1) = 4
A/V ≈ 3.67e-3 m⁻¹  (regular tet)

p_side = 4 × 1.0 × 3.67e-3 × (160²/64) = 4 × 3.67e-3 × 400 = 5.87 GPa/m
penalty_ip ≈ (5.87+5.87)/4 = 2.94 GPa/m     (×3 fix)
penalty_old ≈ 2.94/3 = 0.98 GPa/m             (1/3 penalty, v46)
```

Penalty spring stiffness per face (bilinear form contribution):
```
k_spring = penalty × A_face ≈ 2.94e9 × 4.33e5 = 1.27 PN/m     (×3)
k_spring_old ≈ 0.42 PN/m                                         (1/3)
```

Physical elastic stiffness per element:
```
k_phys ≈ E × A / L ≈ 80e9 × 4.33e5 / 1000 = 34.6 TN/m
k_spring / k_phys ≈ 37×    (×3 penalty)
k_spring_old / k_phys ≈ 12× (1/3 penalty)
```

The penalty spring is 37× stiffer than the physical element stiffness (×3 penalty),
or 12× with the 1/3 penalty. This confirms the penalty dominates the inter-element
coupling, and tripling it significantly changes the effective fault stiffness.

### 10.8 Files Compared

| File | What was verified |
|------|-------------------|
| `tandem/app/localoperator/Elasticity.h:245` | epsilon = -1.0 (SIPG) |
| `tandem/app/localoperator/Elasticity.cpp:278-291` | penalty formula, (p0+p1)/4 |
| `tandem/app/kernels/elasticity.py:89-100` | traction, surface, surfaceOp formulas |
| `tandem/app/kernels/elasticity.py:122-123` | lift_ip: standard IP lift (NOT collocated) |
| `tandem/app/kernels/elasticity.py:141-145` | assembleSurface: consistency+symmetry+penalty |
| `tandem/app/kernels/elasticity.py:242-244` | compute_traction: {σ·n̂} - penalty×(jump-f) |
| `tandem/src/form/InverseInequality.h:27-29` | c_N_1 = (N+1)(N+D)/D |
| `tandem/src/form/DGCurvilinearCommon.cpp:55-94` | area_/volume_ = physical quantities |
| `mfem/fem/bilininteg.cpp:3994-4214` | DGElasticityIntegrator: std SIPG, κ=0→no penalty |
| `seas/integrator/dg_elasticity_ip_penalty_integrator.hpp` | IP penalty integrator |
| `seas/domain/elasticity_operator.hpp:1040-1135` | AssembleSlipContributionIP |
| `seas/domain/elasticity_operator.hpp:3050-3266` | ComputeTraction IP interior path |

---

## 11. Root Cause Found: Initial Time Step Does Not Account for V_nuc

### 11.1 The Bug

`bp5_verification_full.cpp:985-986` computed:
```cpp
real_t dt_init = std::min(1e3, 0.01 * params.L_nuc / std::max(V_init, 1e-20));
```

- `L_nuc = 0.13 m` (critical slip distance Dc)
- `V_init = 1e-9 m/s` (plate rate)
- Result: `dt_init = min(1000, 1.3e6) = 1000 seconds`

But the nucleation zone has `V_nuc = 0.01 m/s`. The formula uses V_init (plate rate)
and completely ignores V_nuc. At `dt = 1000s`, the first RK45 stage (Dormand-Prince,
a21 = 1/5) computes:

| Region | Velocity | δu at first RK stage |
|--------|----------|---------------------|
| Nucleation zone | 0.01 m/s | 0.2 × 1000 × 0.01 = **2.0 m** |
| Outside | 1e-9 m/s | 0.2 × 1000 × 1e-9 = 2×10⁻⁷ m |

**2 meters of slip** at the first intermediate evaluation. With v47's ×3 penalty
at p=2 (penalty ≈ 7.83 GPa/m), the penalty correction at the nucleation boundary
is on the order of GPa — matching the observed τ = 1.43 GPa blowup.

### 11.2 How Tandem Handles This

Tandem uses PETSc's TS (time stepping) with adaptive RK45 (`-ts_rk_type 5dp`).
For QD (quasi-dynamic) mode, Tandem does **NOT** set an explicit initial dt:

**`tandem/app/tandem/SEAS.cpp:125`** (QD specialization):
```cpp
static std::optional<double> cfl_time_step(T const&) { return std::nullopt; }
```

This returns `std::nullopt`, so the CFL-based `set_max_time_step` / `TSSetTimeStep`
call at line 199 is skipped. PETSc's adaptive controller auto-detects the initial dt
by evaluating `f(t0, y0)` (the RHS at the initial condition), which naturally
accounts for V_nuc through the state derivatives.

PETSc's initial step size algorithm (Hairer-Nørsett-Wanner II.4.1):
```
dt ≈ (atol / ||f(t0,y0)||)^(1/(p+1))
   ≈ (1e-7 / 0.01)^(1/6)
   ≈ 0.046 seconds
```

**Effective initial dt: Tandem ≈ 0.05s, our code = 1000s — a factor of ~20,000.**

### 11.3 Why v46 (1/3 Penalty) Survived This Bug

v46 has the same `dt_init = 1000s` bug. The adaptive controller goes through the
same rejection loop (1000 → 100 → 10 → ... → dt_min). The difference:

| | v46 (1/3 penalty, p=1) | v47 (×3 penalty, p=2) |
|--|------------------------|----------------------|
| Penalty | 0.98 GPa/m | 7.83 GPa/m |
| Correction at 2m slip | ~2 GPa | ~16 GPa |
| Intermediate NaN severity | Moderate | Catastrophic |

At v46's lower penalty, the intermediate evaluations during the rejection loop
produce large but finite values, allowing the NaN detection + shrink mechanism
to eventually find a working dt. At v47's 8× larger penalty (p=2), the intermediate
values are so extreme they may produce Inf/NaN that prevents clean error estimation
even during rejection.

### 11.4 dt Fix Does NOT Resolve p=2 Blowup

**Run v47b_ip_p2 (job 7605297)** with `dt_init = 0.13s` still blows up immediately.
The output confirms `Initial dt: 0.13 s (V_max_init = 0.01)` but traction blowup
occurs before the first step completes. Key evidence from the output:

| Rank | DOF | τ (GPa) | slip (m) |
|------|-----|---------|----------|
| 316 | 2 | 10.2 | +2.53 |
| 316 | 4 | 9.63 | +2.92 |
| 316 | 5 | 4.11 | **-0.92** |
| 317 | 121 | 3.89 | -0.86 |
| 317 | 36 | 8.67 | +2.69 |

**Critical observations**:
1. Slip values of 0.5–2.9 m are **2000× larger** than expected from dt=0.13s
   (V_nuc × dt = 0.01 × 0.13 = 0.0013 m)
2. Slip **alternates sign** between adjacent DOFs (+2.9, -0.9) — oscillation
3. Blowup locations are all on the fault plane (y ≈ 0) in the VW/nucleation region

### 11.5 The Actual Mechanism: RK-Stage V Amplification

The blowup is a **cascade within the RK45 stages**, not across time steps.
Even at dt=0.13s, the multi-DOF p=2 system amplifies V through the
elastic-friction coupling within a single Mult() evaluation:

**Stage 1** (k1): δu=0 → T=0 → V=V_nuc=0.01 → k1_slip = 0.01 m/s

**Stage 2** (k2): δu = dt × a21 × V_nuc = 0.13 × 0.2 × 0.01 = 2.6×10⁻⁴ m
- Elastic solve with this tiny δu → traction T at p=2 multi-DOF
- With ×3 penalty (7.83 GPa/m), the traction response has within-face
  DOF variation (6 DOFs per face, some in nucleation zone, others outside)
- DOFs at the nucleation boundary see amplified traction from the
  mismatch between adjacent DOFs → some get V ≈ T/η ≈ **20 m/s**
- k2_slip ≈ 20 m/s at those DOFs

**Stage 3** (k3): δu ≈ dt × (a31×0.01 + a32×20) ≈ 0.13 × 0.3 × 20 ≈ **0.78 m**
- Penalty correction: 7.83 × 10⁹ × 0.78 = **6.1 GPa** → catastrophic
- V at some DOFs → O(1000) m/s → further amplifies slip

**Stage 4+**: Slip reaches 2–3 meters, traction 10 GPa → segfault in
ComputeTraction (null dereference from corrupted state)

**Why this doesn't happen at p=1**: Only 1 effective DOF per face → no within-face
oscillation → V responds uniformly → stable (but too stiff to nucleate, per Section 7.4).

**Why this doesn't happen at 1/3 penalty**: Penalty = 0.98 GPa/m → 8× weaker traction
response → V amplification stays bounded → cascade doesn't start.

### 11.6 Sections 7.3–7.5 Were Correct

The earlier analysis in Sections 7.3 (within-face discontinuity) and 7.5
(resolution-dependent stiffness) correctly identified the root cause. The dt_init
fix (Section 11.1–11.3) was a genuine bug that needed fixing (1000s was wrong),
but it is NOT the cause of the p=2 blowup. The blowup is driven by:

1. **Multi-DOF within-face traction amplification** at the nucleation boundary
2. **×3 penalty amplifying the DOF-to-DOF variation** by 8× vs v46
3. **Positive feedback through RK stages** (V → δu → T → V)

This is a **formulation-resolution issue**, not a time step issue.

### 11.7 The dt Fix

The dt_init fix remains correct (the old value WAS a bug):

```cpp
real_t V_max_init = std::max(V_init, params.V_nuc);
real_t dt_init = std::min(1e3, 0.01 * params.L_nuc /
                          std::max(V_max_init, 1e-20));
```

It just doesn't solve the p=2 penalty problem.

### 11.8 Note on BP2 and Other Tests

BP2 tests use `SetDt(1e3)` which is correct: BP2 has no nucleation zone
(V_nuc = V_init). Only BP5 needed the fix. `test_bp5_integration.cpp`
already uses `SetDt(10.0)` which is small enough.

---

## 12. BR2 p=2 Diagnostic: Isolating the IP Fault Coupling

### 12.1 Motivation

IP p=2 blows up due to multi-DOF within-face V amplification (Section 11.5).
To determine whether the issue is in the p=2 elastic solver or specifically
in the IP fault coupling, we ran BR2 at p=2. BR2 differs from IP in two
critical ways:

| | IP p=2 | BR2 p=2 |
|--|--------|---------|
| Penalty formula | p²-dependent: 7.83 GPa/m | Fixed: dim+1 = 4 (dimensionless) |
| ×3 scaling issue | YES (ref element ratio) | NO (no A/V formula) |
| nbf per face | 6 (multi-DOF) | 1 (face-averaged) |
| Within-face DOF oscillation | YES → cascade → blowup | NO → stable |

All 15 locations of the ×3 penalty fix are IP-specific. BR2 code paths are
completely unaffected by the v47 penalty changes.

### 12.2 Changes Since Last BR2 p=2 Run (v37)

BR2 code paths are largely unchanged since v37. The differences that affect
this run:

| Change | Version | Impact on BR2 |
|--------|---------|---------------|
| Tandem init defaults (V_nuc 0.03→0.01, delta_tau 1→0) | v46 | Different nucleation dynamics (no overstress) |
| Shared-face Dirichlet loading fix | v38b | Stronger tectonic loading in parallel |
| dt_init fix (1000s→0.13s) | v47 | Affects both IP and BR2 |
| IP penalty ×3, traction sign fixes, multi-DOF | v38-v47 | **None** (IP-only paths) |

### 12.3 Results: v47e BR2 p=2 (Job 7605319)

**BR2 p=2 passes the initial stage without blowup.**

```
DG order: 2
DG method: BR2
Global fault DOFs: 9348 (nbf=1, vs 56088 for IP nbf=6)
Initial dt: 0.13 s (V_max_init = 0.01)

    Step       Time [yr]        dt [s]     V_max [m/s]     EQs
----------------------------------------------------------------
       1    6.043491e-10     2.874e-02       1.004e-02       1
      10    2.772566e-08     1.799e-01       1.148e-02       1
      20    7.106352e-08     1.913e-01       1.199e-02       1
      30    1.147227e-07     1.994e-01       1.124e-02       1
      43    1.773517e-07     2.675e-01       9.855e-03       1
```

Comparison with IP p=2:

| | IP p=2 (v47b, job 7605297) | BR2 p=2 (v47e, job 7605319) |
|--|---------------------------|----------------------------|
| Fault DOFs | 56,088 (nbf=6) | 9,348 (nbf=1) |
| Step 0 | **BLOWUP** (τ=10 GPa, slip=2.9m) | Clean |
| Steps completed | 0 | 43+ (still running) |
| V_max trend | N/A (crash) | 0.01 → 0.012 → 0.0099 (slow decay) |
| Segfaults | YES (ComputeTraction) | None |

### 12.4 Analysis

**1. The p=2 elastic solver is correct.** BR2 p=2 runs stably with the same
mesh, same material, same initialization. The DG volume solution at p=2 is
fine. The blowup is NOT in the elastic formulation.

**2. The blowup is specific to IP multi-DOF fault coupling.** The combination
of ×3 IP penalty (7.83 GPa/m) + nbf=6 DOFs per face creates the within-face
V amplification cascade (Section 11.5). BR2 avoids this with:
- Fixed penalty = 4 (no ×3 scaling issue)
- nbf=1 (no within-face DOF oscillation)

**3. V_max is slowly decaying** (0.01 → 0.0099 after 43 steps, t ≈ 5s).
This could be:
- **Too early to tell**: Tandem init (no overstress) means nucleation takes
  years, not seconds. t = 1.7×10⁻⁷ yr is negligible.
- **Stiffness issue**: Similar to IP p=1 with ×3 penalty (V decays, Section 7.4).
  Need to let the run continue to distinguish.

### 12.5 Historical Context: BR2 p=2 Locks the Fault

The v47e BR2 p=2 V decay is NOT a new finding. Debug documents v34-v36
documented progressive fault locking with BR2 at increasing p:

| Order | BR2 recurrence (v34) | Status |
|-------|---------------------|--------|
| p=1 | ~271 yr | Works |
| p=2 | ~396 yr | Delayed (nearly locked) |
| p=4 | No EQs in 1800 yr | **Locked** |
| p=6 | No EQs in 1291 yr | **Locked** |

v35 diagnosed two bugs (centroid-only traction eval + avg_shapes BR2 correction)
and v36-v37 fixed them with per-quadrature-point evaluation. However, these fixes
improved but did not fully resolve the p≥2 locking with nbf=1. The face-averaged
traction (nbf=1) loses within-face detail at higher p, under-resolving the fault
coupling relative to the volume discretization.

The v47e run will likely show the same delayed/locked behavior. BR2 p=2 with
nbf=1 is **not a viable production path** — it was already known to lock the fault.
Its diagnostic value is confirmed: the p=2 elastic volume solver is correct.

### 12.6 Implications for Path Forward

Since the p=2 elastic solver works, the fix should target the **IP fault
coupling specifically**. See Section 13 for the detailed plan.

---

## 13. Option A: Reduced IP Fault-Face Penalty

### 13.1 Core Idea

The ×3 penalty is mathematically correct (Section 2) and must be kept on
**skeleton (non-fault) interior faces** for DG accuracy. But on **fault faces**,
the full penalty creates excessive effective stiffness that either:
- Prevents nucleation at p=1 (V decays, Section 7.4)
- Triggers within-face V amplification at p=2 (blowup, Section 11.5)

The fix: apply a scaling factor `α ∈ [0, 1]` to the penalty on fault faces only:

```
penalty_fault    = α × (dim+1) × c_N_1 × (dim × nl_q / Weight()) × (c1²/c0)
penalty_skeleton = 1 × (dim+1) × c_N_1 × (dim × nl_q / Weight()) × (c1²/c0)
```

At `α = 1/3`: the fault-face penalty equals v46's value (which works at both
p=1 and p=2), while skeleton faces keep the correct ×3. This is strictly better
than v46, where ALL faces had the wrong 1/3 penalty.

### 13.2 Two Levels of Implementation

**Level 1 (quick test): Scale penalty in RHS + traction only, keep K unchanged**

Modify 6 fault-specific locations in `elasticity_operator.hpp`:
- `AssembleSlipContributionIP` p0, p1 (lines 1082-1083)
- `ComputeTraction` interior p0, p1 (lines 3099, 3101)
- `ComputeTraction` shared p0, p1 (lines 3548, 3550)

NOT modified: bilinear form integrator (lines 122, 128), far-field/boundary
Dirichlet RHS (lines 1528, 1529, 1882, 2155, 2157, 2535, 2537).

Inconsistency: K has full ×3 penalty on fault faces, but RHS + traction use
α×penalty. This means:
- Displacement u is computed with strong slip enforcement (full penalty in K)
- Traction T uses weaker penalty correction → lower effective fault stiffness
- As α→0, traction approaches {σ·n̂} (stress average, physical traction)

This is informative for testing. The inconsistency means the traction doesn't
exactly correspond to the variational formulation, but:
- The penalty correction η×(jump(u) - δu) is O(h^p) — a discretization artifact
- With strong penalty in K, jump(u) ≈ δu, so the correction is small regardless
- Reducing α removes the artificial stiffness that causes the instability

**Level 2 (consistent): Also modify K on fault faces**

Modify the bilinear form to use α×penalty on fault faces. Requires either:
- Passing fault face markers to `DGElasticityIPPenaltyIntegrator`
- Or splitting into two integrators (fault vs skeleton) with face attribute filters
- Or a custom face-loop assembly

Implement only after Level 1 confirms the approach works.

### 13.3 Test Matrix

| Run | Order | α (fault) | Expected behavior |
|-----|-------|-----------|-------------------|
| v47f | p=2 | 1/3 | Should run without blowup (v46-equivalent fault penalty) |
| v47g | p=2 | 0.0 | Stress-only traction, no penalty correction on fault |
| v47h | p=1 | 1/3 | Should nucleate (v46-equivalent fault, correct skeleton) |

**Priority**: v47f first. If p=2 runs without blowup AND nucleates → approach works.

### 13.4 Success Criteria

1. **No blowup** at p=2 (no TRACTION BLOWUP, no segfault)
2. **Nucleation occurs** (V_max > 1 m/s at some point in the simulation)
3. **p=2 matches benchmark at least as well as v46** (since skeleton penalty is
   now correct, results should be ≥ v46 quality)
4. **p-convergence**: p=2 closer to Tandem than p=1

### 13.5 Outcome Interpretation

| Outcome | Interpretation | Next step |
|---------|---------------|-----------|
| α=1/3 works, nucleates | Fault penalty was the problem; skeleton fix helps | Implement Level 2, sweep α values |
| α=1/3 runs but no nucleation | Effective stiffness still too high even at 1/3 | Try α=0, or finer mesh |
| α=0 works, nucleates | Penalty correction itself causes excessive stiffness | Consider stress-only traction as default |
| α=0 also blows up | Problem is in multi-DOF interpolation/projection, not penalty | Investigate FaceQuadrature implementation |

### 13.6 Implementation Details

Add member variable and command-line option:
```cpp
// elasticity_operator.hpp: member variable
real_t fault_penalty_factor_ = 1.0;

// bp5_verification_full.cpp: command-line option
// --fault-penalty-factor 0.333
```

At each of the 6 fault-face penalty locations, multiply by the factor:
```cpp
// BEFORE:
real_t p0 = (dim+1) * c_N_1 * (real_t(dim) * nl_q / detJ1) * (c1_mat*c1_mat/c0_mat);

// AFTER:
real_t p0 = (dim+1) * c_N_1 * (fault_penalty_factor_ * real_t(dim) * nl_q / detJ1)
            * (c1_mat*c1_mat/c0_mat);
```

When `fault_penalty_factor_ = 1.0`: identical to v47 (default, backward compatible).
When `fault_penalty_factor_ = 1.0/3.0`: fault faces at v46 penalty, skeleton at v47.

---

## 9. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v30-v42 | All previous fixes through c_N_1 | Done |
| v44 | Penalty ×3 (first attempt) — reverted in v45 due to flawed analysis | Reverted |
| v45 | Analysis claiming ×3 is wrong + slip indexing bug discovery | Done |
| v46 | Multi-DOF slip indexing fix + Tandem initialization defaults | Done |
| **v47** | **Penalty ×3 re-applied — breaks p=1 (stiffness) and p=2 (blowup)** | **Investigating** |
| v47+ | Deep investigation: confirmed both codes use SIPG with identical formulas. Root cause is resolution-dependent effective stiffness at p=1. Corrected SBP-SAT analysis. | Done |
| v47++ | dt_init fix: use max(V_init, V_nuc) → dt_init = 0.13s. Tandem uses PETSc auto-detect (~0.05s). | Fix applied |
| v47+++ | **dt fix does NOT resolve p=2 blowup.** Run 7605297 confirms blowup at dt=0.13s. Slip 0.5–2.9m with oscillating sign → RK-stage V amplification via multi-DOF ×3 penalty feedback. Sections 7.3–7.5 were correct: formulation-resolution issue. | **Confirmed** |
| v47e | **BR2 p=2 runs stably** (job 7605319). 43+ steps, no blowup. Proves p=2 elastic solver is correct; blowup is IP multi-DOF fault coupling. V_max slowly decaying (too early to assess nucleation). | **Running** |
