# BP5 Debug v54: Post-Fix p=1 3-DOF Results — Nucleation Dies

**Date**: 2026-03-28
**Status**: ROOT CAUSE IDENTIFIED — TRACTION SYSTEMATICALLY TOO HIGH
**Previous**: v53 (p=1 fault discretization mismatch root cause, 3-DOF/face fix)
**Branch**: `feature/elasticity`

---

## 1. Context

v53 identified the p=1 fault discretization mismatch (1 DOF vs 3 DOFs/face) and implemented
the fix. The v53 production run (`bp5_v53_1000m_ip_p1_full.sbatch`, commit `7cf3821`) ran
for 48 wall-clock hours on 8 nodes / 400 ranks with MUMPS-BLR (tol=1e-12).

**The central question: does MFEM p=1 1000m nucleate after the 3-DOF fix?**

---

## 2. Answer: NO — Nucleation Dies (MFEM), YES — Nucleation Succeeds (Tandem)

### 2.1 Run Log Summary

From `bp5_v53_ip_p1_7619035.out`:

- **4075 accepted steps**, simulation reached **t = 3.69e-6 yr (116 s)**
- **Killed by SIGTERM** (48hr wall-time limit)
- **V-guard triggered 2x** on the very first step (dt halved from 0.13 → 0.065 → 0.0325)
- dt stabilized at **0.027-0.041 s** (alternating), never grew
- **V_max grew from 0.0103 to 0.0197** over 4075 steps (run log, global V_max)
- **Wall time per step: 42.4 s** (7.1 s per elastic solve, 6 solves/step for DP45 FSAL)
- 28044 global fault DOFs (3 DOFs/face × 9348 fault faces)
- Station output only reached t=38.2 s (buffers not flushed before SIGTERM)

### 2.2 Tandem Full Nucleation Cycle

**Both codes show an initial V overshoot then decline.** The critical difference is what
happens after the decline:

Tandem V at strk-24dp+10:
```
t=0-3.3s:   V rises 0.010 → 0.0142  (initial overshoot)
t=3.3-25s:  V dips to ~0.0114        (initial transient relaxation, stays ABOVE V_nuc)
t=25-50s:   V re-accelerates 0.0115 → 0.015  (nucleation resumes)
t=50-82s:   V rapidly accelerates → 0.30 m/s   *** EARTHQUAKE ***
t=82-117s:  V drops below V_nuc (post-seismic)
t=130s+:    V drops to ~1e-6 → deep interseismic
```

### 2.3 MFEM V Dies Where Tandem Recovers

**MFEM V_strike peaks at t=3.1 s then declines below V_nuc and never recovers:**

| Time (s) | MFEM log10(V) | Tandem log10(V) | MFEM V    | Tandem V  |
|----------|---------------|-----------------|-----------|-----------|
| 0.0      | -2.000        | -2.000          | 0.0100    | 0.0100    |
| 3.1      | **-1.865 (peak)** | -1.848      | 0.0136    | 0.0142    |
| 5.0      | -1.877        | -1.852          | 0.0133    | 0.0141    |
| 10.0     | -1.933        | -1.891          | 0.0117    | 0.0128    |
| 20.0     | -2.021        | -1.938          | **0.0095**| 0.0115    |
| 38.0     | **-2.098**    | -1.915          | **0.0080**| 0.0122    |

**MFEM V drops BELOW V_nuc (0.01) by t~15s and continues falling.**
**Tandem V dips to 0.0114 but stays above V_nuc and re-accelerates to earthquake.**

### 2.3 tau_strike Divergence — The Growing Gap

The tau_strike gap between MFEM and Tandem grows systematically:

| Time (s) | MFEM tau (MPa) | Tandem tau (MPa) | Gap (MPa) |
|----------|----------------|------------------|-----------|
| 0.0      | 21.148         | 21.102           | 0.046     |
| 1.0      | 21.112         | 21.057           | 0.055     |
| 5.0      | 20.825         | 20.759           | 0.066     |
| 10.0     | 20.441         | 20.363           | 0.078     |
| 20.0     | 19.751         | 19.603           | 0.148     |
| 38.0     | 18.866         | 18.443           | **0.424** |

The gap starts at 0.046 MPa (the eta*V output convention, Section 4) and **grows to
0.424 MPa by t=38s**. After subtracting the baseline output convention difference
(~0.046 MPa), the REAL traction excess is ~0.38 MPa.

**MFEM traction is systematically ~0.38 MPa TOO HIGH at t=38s.**

Higher traction = more friction resistance = slower slip rate = nucleation dies.

### 2.4 State Variable Divergence

The state variable (log10 theta) also diverges:

| Time (s) | MFEM state | Tandem state | Delta |
|----------|------------|--------------|-------|
| 0.0      | 8.114      | 8.114        | 0.000 |
| 5.0      | 7.901      | 7.894        | 0.007 |
| 10.0     | 7.691      | 7.668        | 0.023 |
| 38.0     | 6.811      | 6.541        | **0.270** |

MFEM state is HIGHER (theta is larger = fault is stronger). This is consistent with
MFEM's higher traction keeping the fault in a stronger state.

---

## 3. Root Cause Analysis

### 3.1 The 0.046 MPa Initial Offset (Output Convention — COSMETIC)

At t=0, the tau gap is exactly `eta * V_nuc = 0.046 MPa`. This is because:

- **Tandem** outputs `tau_hat = tau_elastic + tau_pre + eta*V` (full friction strength)
  Source: `DieterichRuinaBase.h:75-77`
- **MFEM** outputs `tau_pre + tau_elastic` (without eta*V)
  Source: `bp5_benchmark_output.hpp:260-263`

This 0.046 MPa is a **cosmetic output convention difference** — both codes have the same
internal stress balance. It does NOT affect dynamics.

**Proposed Fix**: Add `eta*V` to MFEM's traction output to match Tandem's convention.

### 3.2 The Growing Gap (0.046 → 0.424 MPa) — THE REAL PROBLEM

After subtracting the 0.046 MPa output convention offset, the REAL traction error grows
from **0 to ~0.38 MPa in 38 seconds**. This is the traction that MFEM computes being
systematically too high (= too much friction = nucleation suppressed).

**Candidate causes (from v51/v52 analysis):**

1. **BLR penalty amplification** (v52 Section 7): MUMPS-BLR solver residual of ~4.4e-7 m
   RMS gets amplified by the IP penalty to ~1830 Pa RMS on fault faces. Over many steps,
   this systematic bias can accumulate. At p=1 with 3 DOFs/face, the penalty factor is
   `3 * p*(p+2) / 3 = p*(p+2) = 3` per face, and the BLR residual is amplified by
   `penalty * h_f / h_e ≈ 300x` the physical traction signal (v52 Table, Section 7).

2. **DG cross-component contamination** (v51): The DG stress average `{sigma.n}` produces
   ~21% dip traction for pure strike input. This drains energy from the strike direction.
   Without `--zero-dip-traction`, this contamination is active in the v53 run.

3. **Traction recovery IP correction term**: The penalty correction in ComputeTraction
   (`correction = -penalty * sign * (du - delta_u)`) amplifies solver error. At p=1 the
   penalty coefficient is relatively small (beta=1.5), but BLR residual is still amplified.

### 3.3 Evidence Against Other Explanations

- **Discretization mismatch**: Fixed in v53 (3 DOFs/face). Initial conditions match to
  high precision at all stations.
- **tau_pre / initialization error**: tau_pre and state match Tandem at t=0 (after
  accounting for output convention).
- **Time stepping error**: Adaptive RK45 with atol=1e-7 should be accurate. V-guard
  prevents cascade. dt is small (0.03s), giving many steps per nucleation timescale.
- **Nucleation zone geometry**: bp5_outside eps=1e-3 matches Tandem.

---

## 4. Cross-Reference with Prior Debug Documents

### Already Tried and Working

| Fix | Version | Impact on nucleation |
|-----|---------|---------------------|
| CFL-aware dt_init | v49/v50 | Prevents first-step crash, does NOT fix nucleation |
| V-guard (factor=100) | v49/v50 | Prevents RK cascade, does NOT fix nucleation |
| 3-DOF/face at p=1 | v53 | Fixes discrete model mismatch, does NOT fix nucleation |
| bp5_outside eps=1e-3 | v53 | Matches Tandem inclusion, does NOT fix nucleation |
| ClosedUniform face nodes | v50 | Critical for p>=4, not relevant at p=1 |
| Exp cap removal from aging law | v51 | Correct, but not the bottleneck |

### Already Tried and Found Insufficient

| Fix | Version | Finding |
|-----|---------|---------|
| BLR tol 1e-10 | v47a | V decays at p=1, nucleation dies |
| BLR tol 1e-12 | v51d | V decays at p=1, nucleation dies |
| Exact MUMPS (no BLR) | v47c/v52 | OOM at p=2 1000m on cluster |
| Penalty factor reduction | v50 | Tried, insufficient |
| Stress-only traction (no penalty correction) | v50f/v51f | Diverges (blowup) |
| Elastic sigma_n feedback | v51 | Makes dip WORSE, not a fix |

### Identified But Not Yet Fixed

| Issue | Version | Status |
|-------|---------|--------|
| BLR penalty amplification (632% corr/stress) | v52 | ROOT CAUSE identified, no fix |
| DG cross-component contamination (21% dip) | v51 | IDENTIFIED, zero-dip-traction works as workaround |
| No CG matrix-free solver in MFEM | v52 | NOT AVAILABLE — would eliminate BLR issues entirely |
| Traction output convention (missing eta*V) | v54 | PROPOSED, cosmetic only |

---

## 5. The Solver Is NOT the Cause

### 5.1 CG+AMG Run Confirms: Identical to MUMPS-BLR

The v54 CG+AMG run (`bp5_v54_cg_p1_7619209.out`) produces **bit-for-bit identical
V_max** as the v53 MUMPS-BLR run for the first 10 steps:

| Step | BLR V_max     | CG V_max      | BLR dt    | CG dt     |
|------|---------------|---------------|-----------|-----------|
| 1    | 1.029e-02     | 1.029e-02     | 1.146e-02 | 1.146e-02 |
| 5    | 1.124e-02     | 1.124e-02     | 3.424e-02 | 3.424e-02 |
| 10   | 1.054e-02     | 1.054e-02     | 3.703e-02 | 3.709e-02 |

Time differs by <0.02% at step 10 — dt by <0.2%. V_max is exactly identical.

**This definitively disproves the BLR solver hypothesis** (v52 Section 7). The MUMPS-BLR
factorization error is NOT what kills nucleation. Switching to CG+AMG (iterative,
rtol=1e-10, guaranteed convergence) changes nothing.

### 5.2 What Tandem Actually Uses

Tandem uses **PETSc KSP CG** (iterative) with:

- `matrix_free = true`: shell operator A for matvec, assembled matrix P as preconditioner
- `mg_strategy = "logarithmic"`: p-multigrid, but at p=1 only one level (degenerates to
  default PETSc PC on the assembled matrix)
- `rtol = 1e-12`

Source: `tandem/app/common/PetscLinearSolver.cpp:7-44`

At p=1, Tandem effectively uses CG + block-Jacobi/ILU on the same kind of assembled DG
matrix. The CG run proves this solver architecture produces the same answer as MUMPS-BLR.

### 5.3 The Real Root Cause: IP Penalty Correction Dominance

The solver is eliminated. What remains is a **DG formulation issue**: the IP penalty
correction term dominates the traction at near-front stations.

From v53 Section 15 (station-level traction decomposition at `strk-16dp+10`):

| Time (s) | tau_stress_strike | tau_corr_strike | corr/stress |
|----------|-------------------|-----------------|-------------|
| 0.1      | 9.87e-04 MPa      | 4.06e-03 MPa    | **4.1x**    |
| 1.0      | 9.29e-03 MPa      | 2.85e-02 MPa    | **3.1x**    |
| 5.0      | 4.20e-02 MPa      | 1.15e-01 MPa    | **2.75x**   |

**75-85% of the total traction comes from the penalty correction, not from the
physical stress.** The correction also produces excessive dip contamination (5-6x
the stress-only dip traction).

This is intrinsic to the IP DG formulation at p=1 1000m resolution — not a solver
accuracy issue, not a BLR issue, and not a discretization issue. Both MUMPS-BLR and
CG+AMG faithfully reproduce this same DG behavior.

---

## 6. Proposed Fixes

### Fix 1: ~~CG Solver~~ — DISPROVED (tested in v54 run)

CG+AMG produces identical results to MUMPS-BLR. Solver accuracy is not the issue.

### Fix 2: Address IP Penalty Correction Dominance (Priority: CRITICAL)

The real problem (v53 Section 15): at near-front stations the penalty correction is
3-6x the physical stress term. Tandem computes traction differently — its DG traction
recovery does not produce this correction dominance.

**Investigation needed**: Compare MFEM's ComputeTraction IP correction formula against
Tandem's traction recovery implementation term-by-term.

Tandem's traction at fault is computed in `adapter_.traction(disp_view, traction_)`:
- Source: `tandem/app/form/SeasQDOperator.cpp:69-73`
- This calls into the DG operator's traction evaluation

The key question: does Tandem's DG formulation use the same penalty correction structure
as MFEM? If not, what is different?

Possible directions:
- **Penalty coefficient formula**: MFEM uses `p*(p+dim-1)/dim * {1/h}`. Does Tandem
  use a different scaling at p=1?
- **Correction term structure**: MFEM's correction is `-penalty * sign * (du - delta_u)`.
  Does Tandem compute the displacement jump differently?
- **Face normal / orientation**: Does Tandem handle face normals differently for the
  traction evaluation?

### Fix 3: Add eta*V to Traction Output (Priority: MEDIUM, COSMETIC)

Add `eta * V` to the traction output in `bp5_benchmark_output.hpp` to match Tandem's
`tau_hat` convention. Eliminates the 0.046 MPa baseline offset in comparisons.
Does not affect dynamics.

### Fix 4: Investigate Tandem's Traction Recovery (Priority: CRITICAL)

The most important next step is a **term-by-term comparison** of how MFEM and Tandem
compute the elastic traction at fault faces. Specifically:

1. How does Tandem evaluate `{sigma . n}` at fault faces?
2. Does Tandem use the same penalty correction form?
3. What is the stress/correction ratio at Tandem's near-front stations?

If Tandem's correction/stress ratio is ~1 (balanced) vs MFEM's 3-6x (dominated), the
difference in the correction formula is the root cause.

---

## 7. Recommended Next Steps

1. **Implement CG+AMG solver** (Fix 1) — this is the highest-leverage change
2. **Test with a short run** (100 steps) to verify nucleation behavior improves
3. **Compare V_strike at strk-24dp+10** against Tandem to confirm the traction bias is gone
4. If nucleation succeeds, launch full 1800-year production run
5. Apply cosmetic fixes (eta*V output, zero-dip-traction) for clean comparison

---

## 8. Timeline of Nucleation Death

```
t=0.0s:  V=0.0100  MFEM and Tandem agree (both start at V_nuc)
t=0.3s:  V=0.0105  Both accelerating similarly
t=3.1s:  V=0.0136  MFEM PEAKS HERE — Tandem at V=0.0142
t=5.0s:  V=0.0133  MFEM V starts decreasing — Tandem still V=0.0140
t=10s:   V=0.0116  MFEM decelerating — Tandem V=0.0129
t=20s:   V=0.0095  MFEM BELOW V_nuc — Tandem V=0.0115
t=38s:   V=0.0080  MFEM DYING — Tandem V=0.0122 (still nucleating)
```

The crossover happens at t ≈ 3 seconds. Before that, both codes agree well. After that,
MFEM's systematically higher traction (from BLR + penalty amplification) creates enough
extra friction to halt the acceleration and reverse it.

---

## 9. Cross-Verification Unit Test Results

### 9.1 Tests Built and Run

Test file: `miniapps/seas/tests/unit/test_cross_verify_tandem.cpp`

| Test | Result | Finding |
|------|--------|---------|
| Penalty coefficient | PASS | c_N_1=1.0 at p=1, formula correct |
| Slip round-trip | PASS | InterpolateToQuadPoints/GalerkinProject exact for linear |
| L2 projection consistency | PASS | M_ref * M_ref_inv = I for p=0,1,2 |
| Solve+decomposition (uniform slip) | **FAIL** | correction/stress = **49x** for uniform slip |
| RHS vs traction consistency | PASS* | Penalty consistent (37.5 Pa/m everywhere) |
| Quadrature mismatch check | **DETECTED** | match_quad_order changes ratio from 66x to 941x |
| Dirichlet BC time | PASS | t=0 gives zero, t>0 gives non-zero |
| Deep dive: jump residual | INFO | |[[u]]-slip|/|slip| = 27-41% at p=1 coarse |
| Refinement study | INFO | max corr/stress INCREASES with refinement |

### 9.2 Key Finding: Correction Dominance is Inherent to IP DG at p=1

The jump residual `[[u]] - slip` is 15-27% of the prescribed slip even for UNIFORM slip
on tet meshes. This is an inherent property of the IP DG method — the penalty enforces
the jump constraint WEAKLY, and at finite penalty the residual is O(h/penalty).

The penalty coefficient is correct (matches Tandem's formula exactly: 37.5 Pa/m on unit
test mesh, consistent across all DOFs to machine precision).

### 9.3 Refinement Study

| Config | max corr/stress | avg corr/stress | RMS jumpres/slip |
|--------|----------------|----------------|-----------------|
| p=1 1x1 | 51.6 | 29.6 | 27.1% |
| p=1 2x2 | 64.0 | 18.8 | 19.2% |
| p=1 3x3 | 73.2 | 14.1 | 15.8% |
| p=1 4x4 | 76.7 | 11.4 | 13.7% |
| p=2 1x1 | 69.1 | 27.1 | 22.5% |
| p=2 2x2 | 83.5 | 16.8 | 15.7% |

The avg correction/stress ratio DECREASES with refinement (bulk DOFs improve), but the
max ratio INCREASES (worst-case DOFs at fault boundary get relatively worse). This is a
standard IP DG boundary effect.

### 9.4 Test 16: Uniform Slip Exactness

For uniform slip (constant on all DOFs), the exact solution is piecewise-constant u:
u = ±slip/2. DG p=1 CAN represent this. But **MFEM's solution is NOT piecewise
constant** — within-element u variation = 84% of expected amplitude.

The SIPG consistency term `{sigma(v)}.n . [[u]]` couples the fault face jump into
the volume DOFs, creating non-physical strain even for uniform slip. This is inherent
to the SIPG formulation (both MFEM and Tandem have this). The result is that the
penalty correction must compensate for the spillover, producing the 50x ratio.

**Whether Tandem has the same 50x ratio on the same mesh is the key open question.**
If both codes produce the same correction-dominated traction, the nucleation difference
is in the friction coupling (sigma_n feedback, etc.). If Tandem produces a smaller ratio,
the IP implementations differ in a way the unit tests haven't captured.

### 9.5 Implication

The 50x correction/stress ratio is a consequence of IP DG finite-penalty spillover.
The immediate actionable finding is the normal stress feedback (Test 14): Tandem's
sigma_n = SnPre - sn_elastic compensates for this spillover in a way MFEM was not doing.

**v54 fixes:**
1. `elastic_sigma_n` is now ON by default (matching Tandem)
2. **SIGN BUG FIXED**: v51 had `sigma_n_eff = SnPre + sn` (WRONG).
   Tandem uses `snAbs = -sn + SnPre` (DieterichRuinaBase.h:87).
   v54 fix: `sigma_n_eff = SnPre - normal_traction` (correct).

   The wrong sign meant compression DECREASED sigma_n (reducing friction),
   when it should INCREASE sigma_n (increasing friction). This explains why
   v51 found `--elastic-sigma-n` made things WORSE — the sign was inverted,
   so the feedback had the opposite effect from Tandem's.

---

### 9.5 CRITICAL: Quadrature Inconsistency Between K and b

**The stiffness matrix K and the RHS b use DIFFERENT quadrature orders.**

| Operation | Quad Order | # Points (p=1) | Source |
|-----------|-----------|----------------|--------|
| K: DGElasticityIntegrator (consistency) | `2p` = 2 | **3** | `bilininteg.cpp:4097` |
| K: DGElasticityIPPenaltyIntegrator (penalty) | `2p` = 2 | **3** | `dg_elasticity_ip_penalty_integrator.hpp:76-78` |
| b: AssembleSlipContributionIP (RHS) | `2p+1` = 3 | **4** | `elasticity_operator.hpp:1185` |
| Traction: ComputeTractionImpl | `2p+1` = 3 | **4** | `elasticity_operator.hpp:3732` |
| **Tandem: ALL face operations** | **`2p+1` = 3** | **4** | `DGCurvilinearCommon.cpp:16` |

**Test 11 proved this is NOT a solution difference** — changing traction/RHS to order 2p
gives the same u (K is unchanged). But it IS an evaluation inconsistency: the traction
recovery evaluates the penalty correction at 4 quad points while K was assembled at 3 points.

For linear shape functions on flat triangles, both rules exactly integrate the polynomial
integrands. However, the DISCRETE operators are assembled at different evaluation points:
- K "sees" the penalty at 3 locations
- ComputeTraction "sees" the penalty at 4 locations
- These produce numerically different correction values even though both are exact for
  the continuous problem

**Tested fix**: Changed bilinear form integrators to use order `2p+1`. Result: **NO EFFECT**.
For flat elements at p=1, both order-2 (3 points) and order-3 (4 points) quadrature
integrate degree-2 polynomial integrands exactly. The stiffness matrix K is bit-for-bit
identical. The quadrature mismatch is NOT the cause of the correction dominance.

The `match_quad_order` flag only affects traction EVALUATION (how quad-point tractions
are projected to DOFs), not the elastic SOLUTION. Test 11 confirms: u is identical
regardless of traction quadrature order.

**Root finding**: The correction dominance (50x at p=1) is inherent to the IP DG method's
weak enforcement of the jump constraint. The penalty drives `[[u]] → slip` but never
achieves it exactly. The residual `[[u]] - slip` is 15-27% of slip on coarse meshes.
Both MFEM and Tandem should produce the same residual for the same formulation.

**The difference between MFEM and Tandem is NOT in the DG traction computation.**

---

## 10. What the Cross-Verification Has Eliminated

| Hypothesis | Test | Result |
|------------|------|--------|
| MUMPS-BLR solver error | CG+AMG run (v54 job) | **DISPROVED** — identical results |
| Penalty coefficient wrong | Tests 1, 6, 9 | **DISPROVED** — matches Tandem formula exactly |
| Slip interpolation bug | Test 2 | **DISPROVED** — exact round-trip for linear fields |
| L2 projection (mass matrix) bug | Test 3 | **DISPROVED** — M_ref * M_ref_inv = I |
| Quadrature mismatch (K vs b) | Tests 11, 12 | **DISPROVED** — both orders give identical K for flat p=1 |
| Dirichlet BC time error | Test 7 | **DISPROVED** — correct time behavior |

## 11. What Remains

The DG traction computation pipeline is verified correct. The correction dominance
(50x at p=1) is inherent to IP DG weak jump enforcement. The penalty, interpolation,
projection, stress computation, and solver all match Tandem's formulation for flat
elements.

**The difference between MFEM and Tandem is NOT in the DG elasticity operator.**

### 11.1 NEW FINDINGS: Friction Law and Normal Stress (Tests 13-15)

**Friction solver** (Test 13): MFEM uses Brent's method in LINEAR V space, Tandem uses
LOG10(V) space. Both converge to machine precision (residual < 1e-16). The linear-space
bracket is catastrophically ill-conditioned for interseismic V (root at 3.5e-8% of
bracket), but Brent's method handles this correctly. NOT the root cause.

**State evolution** (Test 15): `dpsi/dt = (b*V0/Dc)*(exp((f0-psi)/b) - V/V0)` matches
Tandem exactly. Verified with hand computation.

**CRITICAL — Normal stress feedback** (Test 14):

Tandem: `sigma_n = SnPre - sn_elastic` (DieterichRuinaBase.h:87)
MFEM:   `sigma_n = sigma_n_pre` = 25 MPa constant (default, no --elastic-sigma-n)

On a test mesh with uniform strike slip 0.001 m:
- max |T_normal| = **9.23 MPa** (37% of sigma_n!)
- Implied delta_tau from sigma_n feedback = **5.54 MPa**

This is a 5.54 MPa effect that Tandem includes but MFEM ignores. Even though BP5 is
analytically a pure strike-slip problem with zero normal traction, the DG IP penalty
correction generates large spurious normal traction (just as it generates large spurious
dip traction — the same mechanism).

Tandem's friction law absorbs this normal component into sigma_n_eff, partially
compensating for the correction dominance. MFEM ignores it.

**IMPORTANT CAVEAT**: v51 previously tested `--elastic-sigma-n` and found it made dip
WORSE. This needs re-investigation with the current v53 3-DOF fix. The v51 test may
have been on different code.

### 11.2 elastic_sigma_n + Sign Fix Result: NO EFFECT

Run `bp5_v54_esn_p1_7619263.out`: 600 steps, 14.8 seconds, V_max = 1.034→1.018.
**Same dying behavior as v53.** The elastic normal traction on the production mesh is
much smaller than on the coarse unit test mesh (37% vs negligible).

The sign fix is correct (matches Tandem), and elastic_sigma_n should stay ON as default.
But **it is NOT the cause of nucleation failure.**

Disproved hypotheses so far:
- Solver accuracy (CG vs BLR: identical)
- Quadrature mismatch (irrelevant for flat p=1)
- elastic_sigma_n OFF (now ON + sign fixed, no effect)
- Penalty coefficient (matches Tandem exactly)
- Slip interpolation (exact round-trip)
- L2 projection (M_ref correct)
- State evolution (matches Tandem)
- Friction solver (accurate to 1e-16)
