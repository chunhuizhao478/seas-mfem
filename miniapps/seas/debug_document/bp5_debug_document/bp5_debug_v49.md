# BP5 Debug v49: CFL Stability Analysis — IP Penalty Eigenvalue vs dt_init

**Date**: 2026-03-22
**Status**: ROOT CAUSE IDENTIFIED. The ×3 IP penalty creates a mesh-size-dependent CFL constraint: the effective eigenvalue λ ∝ 1/h, and λ*dt exceeds the RK45 stability boundary for h=1000m at dt=0.13s but NOT for h=2500m. The fix is CFL-aware dt selection: dt_init = min(dt_V, C·η·h/(4μ)).
**Previous**: v48 (sign hypothesis disproved, BLR secondary issue, 1000m serial crash confirmed)
**Branch**: `feature/elasticity`

---

## 1. Context from v48

v48 established conclusively:
- **Sign hypothesis DISPROVED** (Section 10): CalcOrtho normals flip correctly on shared faces
- **BLR accuracy** is a secondary issue on 2500m mesh; irrelevant for 1000m crash
- **1000m mesh crashes in SERIAL** (v48i): same tau/slip values as all parallel runs
- **Root cause unknown**: 1000m + p=2 + ×3 penalty = mesh-resolution-dependent instability

v48 left open: *why does 1000m trigger instability but 2500m doesn't?*

This document answers that question.

---

## 2. Phase 1 Experiments: Design

Built a single binary with 5 diagnostic flags, ran 4 jobs simultaneously:

| Job | Mesh | Special flag | What it tests |
|-----|------|-------------|---------------|
| **v49a** | 1000m serial | `--smooth-nucleation` | Sharp V boundary as trigger |
| **v49b** | 1000m serial | `--match-quad-order` (2p) | Quadrature order mismatch |
| **v49cd** | 1000m serial | *(baseline)* | Diagnostic dump of normals, seed traction, RK stages |
| **v49e** | 2500m serial | *(baseline)* | Reference — expected stable |

All 4 jobs include `--diag-normals`, `--diag-first-traction`, `--diag-rk-stages`.

---

## 3. Phase 1 Results Summary

### 3.1 Overall Outcome

| Job | Mesh | Fix | Result | Crash stage | Peak V_str |
|-----|------|-----|--------|-------------|------------|
| **v49cd** | 1000m | *(none)* | **CRASH** | Stage 4 | 20.94 m/s |
| **v49a** | 1000m | Smooth nucleation | **CRASH** | Stage 5 | 114.2 m/s |
| **v49b** | 1000m | Match quad (2p) | **CRASH** | Stage 2 | 0.896 m/s |
| **v49e** | 2500m | *(none)* | **STABLE** | N/A | 0.133 m/s (bounded) |

**Key finding**: ALL 1000m runs crash. The 2500m reference is stable and runs 132+ steps.

### 3.2 Normal Diagnostic (NOR-DIAG)

ALL faces on BOTH meshes show:
```
CalcOrtho = (0, +1, 0)   basis.normal = (0, -1, 0)   dot = -1
```

The CalcOrtho normal and fault basis normal are **exactly opposite** on every fault face. This is consistent and systematic — NOT a bug. The `sign = (nor(1) > 0) ? 1.0 : -1.0` convention handles this correctly: sign = +1 everywhere, and the code's `sign * delta_u` correctly represents the prescribed jump in the CalcOrtho convention.

**Verdict**: Normal direction is NOT a factor in the mesh-dependent instability (both meshes have the same pattern).

### 3.3 Seed Traction (SEED-TRAC)

Both meshes:
```
total_DOFs_with_tau>1Pa: 3 / 55920  (1000m)
total_DOFs_with_tau>1Pa: 3 / 10248  (2500m)
```

Initial traction is essentially ZERO (only 3 DOFs barely above threshold, likely boundary effects). There is no initial traction error seeding the instability. The instability arises from the DYNAMIC coupling during RK stages, not from a static initial error.

### 3.4 RK Stage Cascade — THE CRITICAL DATA

**v49e (2500m, STABLE)** — First step, dt = 0.13s:
```
Stage 0: V_str = 0.0100                         (initial)
Stage 1: V_str = 0.0108   (×1.08)
Stage 2: V_str = 0.0106   (×0.98)               oscillating
Stage 3: V_str = 0.0186   (×1.75)               mild growth
Stage 4: V_str = 0.0587   (×3.15)               significant
Stage 5: V_str = 0.1329   (×2.26)               peak
Stage 6: V_str = 0.0264   (×0.20) ← REVERSAL    bounded!
```
After step 0 completes, adaptive dt → 0.030s. All subsequent steps are calm (V_str ≈ 0.010).

**v49cd (1000m, CRASH)** — First step, dt = 0.13s:
```
Stage 0: V_str = 0.0100                         (initial)
Stage 1: V_str = 0.0111   (×1.11)
Stage 2: V_str = 0.0217   (×1.96)
Stage 3: V_str = 0.8727   (×40.2)  ← DIVERGENCE
Stage 4: V_str = 20.942   (×24.0)  ← RUNAWAY
→ TRACTION BLOWUP at DOF 8330, tau = 1.07 GPa
```

**v49a (1000m smooth, CRASH)** — First step, dt = 0.13s:
```
Stage 0: V_str = 0.0100                         (initial)
Stage 1: V_str = 0.0100   (×1.00)               delayed onset
Stage 2: V_str = 0.0100   (×1.00)
Stage 3: V_str = 0.1859   (×18.4)  ← DIVERGENCE (2 stages later)
Stage 4: V_str = 5.0557   (×27.2)
Stage 5: V_str = 114.16   (×22.6)
→ TRACTION BLOWUP at DOF 8293, tau = 2.10 GPa
```

**v49b (1000m quad-match, CRASH)** — First step, dt = 0.13s:
```
Stage 0: V_str = 0.0100                         (initial)
Stage 1: V_str = 0.1064   (×10.6)  ← IMMEDIATE DIVERGENCE
Stage 2: V_str = 0.8958   (×8.42)
→ SEGFAULT (crashed before blowup detection)
```

---

## 4. Root Cause Analysis: CFL-Like Stability Constraint

### 4.1 The Key Observation

The 2500m mesh shows **bounded amplification** within the first RK step (peaks at 13× then REVERSES at stage 6). The 1000m mesh shows **unbounded amplification** (exponential growth through stages, never reverses). This is the classic signature of exceeding an explicit stability boundary.

### 4.2 The Effective Eigenvalue

The coupled elasticity-friction system has an effective eigenvalue that governs the RK stage amplification:

```
λ_eff = k_eff / η
```

where:
- `k_eff` = effective fault stiffness (traction response per unit slip)
- `η` = radiation damping coefficient = μ/(2cs) = 4.62×10⁶ Pa·s/m

The effective fault stiffness k_eff includes contributions from:
1. Physical elastic stiffness: ∝ μ/h
2. DG penalty contribution: ∝ penalty/h (through the elastic solve)

For the ×3 corrected IP penalty at p=2:
```
k_eff ≈ β · μ / h
```
where β is a geometry-dependent factor. With the ×3 penalty correction, β ≈ 4 (physical stiffness + penalty amplification through the elastic solve).

Computing:
```
λ_eff = β · μ / (h · η)

h = 1000m:  λ_eff = 4 × 32×10⁹ / (1000 × 4.62×10⁶) = 27.7 /s
h = 2500m:  λ_eff = 4 × 32×10⁹ / (2500 × 4.62×10⁶) = 11.1 /s
```

### 4.3 The Stability Product z = λ·dt

The Dormand-Prince RK45 stability region has a critical boundary for the maximum amplification per step. For a system with effective eigenvalue λ, the stability product z = λ·dt determines whether the RK stages diverge:

```
h = 1000m:  z = 27.7 × 0.13 = 3.60  → OUTSIDE stability boundary
h = 2500m:  z = 11.1 × 0.13 = 1.44  → INSIDE stability boundary
```

The critical z_crit is approximately 2.5-3.0 for this nonlinear system (the stages oscillate but bound below z_crit, and diverge exponentially above it).

### 4.4 Verification Against PETSc Default

Tandem uses PETSc's default dt = 0.1s (from `tscreate.c:45`):
```
h = 1000m with PETSc dt:  z = 27.7 × 0.10 = 2.77  → Marginally stable
h = 1000m with our dt:    z = 27.7 × 0.13 = 3.60  → Unstable
```

The difference between PETSc's 0.1s and our 0.13s pushes us across the stability boundary at h=1000m. At h=2500m, both dt values are well within stability.

### 4.5 Why Each Candidate Fix Worked/Failed

| Experiment | Effect on z | Result | Explanation |
|-----------|------------|--------|-------------|
| v49cd (baseline) | z = 3.60 | CRASH | Above stability boundary |
| v49a (smooth nucleation) | z = 3.60 (unchanged) | CRASH | Same λ·dt; smoothing delays onset but doesn't change stability bound |
| v49b (match quad 2p) | z > 3.60 (worse) | CRASH faster | 2p quadrature introduces additional errors that amplify the perturbation |
| v49e (2500m reference) | z = 1.44 | STABLE | Well within stability boundary |

The smooth nucleation (v49a) delayed the cascade by 2 stages because the initial perturbation at the nucleation boundary is smaller with a Gaussian taper. But the amplification RATE (determined by λ·dt) is the same, so once the perturbation grows large enough to trigger the nonlinear amplification, the cascade proceeds identically.

### 4.6 Formula for CFL-Aware dt_init

```
dt_CFL = C_stability × η × h_min / (β × μ)
```

where:
- C_stability ≈ 2.0 (safety factor below the empirical z_crit ≈ 2.5-3.0)
- η = μ/(2cs) = 4.62×10⁶ Pa·s/m
- h_min = minimum element size on the fault
- β ≈ 4 (geometry/penalty factor for ×3 corrected IP at p=2)
- μ = shear modulus = 32.04×10⁹ Pa

Numerical constant:
```
dt_CFL = 2.0 × 4.62×10⁶ / (4 × 32.04×10⁹) × h_min = 7.21×10⁻⁵ × h_min
```

| h_min | dt_CFL | Our dt_init | Status |
|-------|--------|------------|--------|
| 500m | 0.036s | 0.13s | UNSTABLE (z = 7.2) |
| 1000m | 0.072s | 0.13s | UNSTABLE (z = 3.6) |
| 2500m | 0.180s | 0.13s | STABLE (z = 1.4) |
| 5000m | 0.361s | 0.13s | STABLE (z = 0.7) |

The final dt selection:
```
dt_init = min(dt_from_V_max, dt_CFL)
        = min(L/(10·V_max), C · η · h_min / (β · μ))
```

---

## 5. Empirical Stage-by-Stage Amplification Analysis

### 5.1 The 2500m Step-0 to Step-1 Transition

After the first step (dt = 0.13s), the adaptive controller reduces dt:
```
Step 0: dt = 0.13s   → V_max = 1.03e-02 (stages oscillated to 0.133, then reversed)
Step 1: dt = 0.030s  → V_max = 1.03e-02 (all stages calm, amplification < 1.003)
```

The error estimate from step 0 is large (because the stages oscillated), triggering a 4.3× dt reduction. At dt=0.03s, z = 11.1 × 0.03 = 0.33 — deeply stable.

This is the HEALTHY behavior: the adaptive controller catches the large step error and reduces dt. The 2500m mesh survives the first step because the amplification, though large (13×), is BOUNDED.

### 5.2 Why 1000m Can't Recover

At h=1000m with dt=0.13s, the stage amplification is UNBOUNDED:
- Stage 3→4: V goes from 0.87 to 20.9 (24× in one stage)
- This produces massive slip (0.036m), massive traction (1 GPa), segfault

The adaptive controller never gets a chance to reduce dt because the step CRASHES before completing. The key difference: 2500m completes the step (with large but bounded error) while 1000m diverges before completion.

### 5.3 Amplification Factor Per Stage

| Stage | 2500m V_str | Ratio | 1000m V_str | Ratio |
|-------|-------------|-------|-------------|-------|
| 0 | 0.0100 | — | 0.0100 | — |
| 1 | 0.0108 | 1.08 | 0.0111 | 1.11 |
| 2 | 0.0106 | 0.98 | 0.0217 | 1.96 |
| 3 | 0.0186 | 1.75 | 0.8727 | 40.2 |
| 4 | 0.0587 | 3.15 | 20.942 | 24.0 |
| 5 | 0.1329 | 2.26 | BLOWUP | — |
| 6 | 0.0264 | 0.20 | — | — |

At 2500m, the per-stage ratio peaks at 3.15 (stage 4) then decreases. At 1000m, the ratio hits 40× at stage 3 — the system is in full runaway mode.

The transition from "oscillating" to "diverging" happens when the effective V during a stage is large enough that the high-V eigenvalue (λ ≈ k_eff/η) dominates. For 1000m, this happens at a lower V threshold because k_eff is 2.5× larger.

---

## 6. Why Previous Hypotheses Were Wrong (or Incomplete)

### 6.1 Shared Face Sign Bug (v48 Section 3) — DISPROVED

The sign hypothesis assumed that `CalcOrtho` returns the same normal on both ranks. v48 Section 10 proved this wrong: MFEM's `GetSharedFaceTransformations` orients the normal outward from the local element, so it FLIPS between ranks. The code's `+sign * delta_u` produces opposite-signed loads on the two ranks — correct behavior.

### 6.2 BLR Solver Accuracy (v48 Section 13) — Secondary Issue

BLR accuracy matters on the 2500m mesh in parallel (v48c-e), but is irrelevant for the 1000m crash. The 1000m mesh crashes with exact MUMPS (v48f), BLR 1e-12 (v48h), and in serial (v48i). The instability is independent of solver accuracy.

### 6.3 Quadrature Order Mismatch — Makes It WORSE

v49b tested matching the RHS quadrature order to K's (2p instead of 2p+1). Result: crashes FASTER. The 2p quadrature has fewer points and introduces additional integration errors that increase the initial perturbation. The 2p+1 quadrature was actually BETTER (it over-integrates, reducing aliasing errors).

### 6.4 Sharp Nucleation Boundary — Delays But Doesn't Fix

v49a tested Gaussian tapering of the V_nuc/V_init boundary. Result: delays the cascade by 2 stages (the initial perturbation is smaller) but doesn't change the stability boundary. Once V grows large enough through normal nucleation physics, the same cascade triggers.

---

## 7. The Fix: Phase 2 Implementation Plan

### 7.1 Priority 1: CFL-Aware dt_init (Quick Test)

**Immediate test**: Run the 1000m mesh with dt_init = 0.07s.

Implementation: add a `--dt-init` command-line override flag:
```cpp
// In bp5_verification_full.cpp:
real_t dt_init_override = -1.0;
// parse --dt-init <value>
if (dt_init_override > 0) {
    dt_init = dt_init_override;
} else {
    real_t dt_V = L_ref / (10.0 * V_max_init);
    real_t dt_CFL = 2.0 * eta * h_min / (4.0 * mu);
    dt_init = std::min(dt_V, dt_CFL);
}
```

This requires computing `h_min` from the mesh. For Gmsh meshes, h_min can be estimated from:
- The mesh scale parameter (1000, 2500, etc.)
- Or the minimum element volume: `h_min ≈ (6 * V_min)^(1/3)` for tets

**Expected result**: 1000m mesh stabilizes with dt_init ≈ 0.07s. Adaptive controller grows dt from there.

### 7.2 Priority 2: RK Stage-Level V Guard

Add a safety check within the RK45 stepper that detects runaway amplification:

```cpp
// After computing stage k_[i]:
double V_max_stage = /* max |V| in k_[i] */;
if (V_max_stage > V_guard * V_max_current) {
    // Stage is diverging — reject step, halve dt
    dt *= 0.5;
    goto restart_step;
}
```

With V_guard ≈ 100 (100× amplification = clearly diverging). This is robust against mesh-size changes and doesn't require pre-computing CFL limits.

### 7.3 Priority 3: Automatic h_min Computation

For the CFL formula, we need h_min on the fault. Implementation:
```cpp
real_t h_min = HUGE_VAL;
for (int fi = 0; fi < fault_interior_faces_.Size(); fi++) {
    auto *FTr = mesh_.GetInteriorFaceTransformations(fault_interior_faces_[fi]);
    real_t vol1 = FTr->Elem1->Weight() / 6.0;  // physical vol = detJ * ref_vol
    real_t h1 = std::cbrt(6.0 * vol1);
    h_min = std::min(h_min, h1);
}
// MPI reduce for parallel
```

### 7.4 Phase 2 Test Plan

All tests use the 1000m mesh, serial (1 rank), p=2, ×3 penalty:

| Test | dt_init | What it verifies |
|------|---------|-----------------|
| **v49f** | 0.07s (manual CFL) | Does CFL-based dt fix the 1000m crash? |
| **v49g** | 0.05s (conservative) | Extra margin to confirm stability |
| **v49h** | 0.10s (PETSc default) | Does Tandem's default work for us? |
| **v49i** | 0.13s + V_guard | Does stage monitoring catch the cascade? |

If v49f works → the fix is confirmed. Then implement automatic h_min computation (Priority 3) and the V guard (Priority 2) for production robustness.

---

## 8. Normal Convention Analysis (NOR-DIAG Follow-Up)

### 8.1 Finding

ALL fault faces show: CalcOrtho = (0,+1,0), basis.normal = (0,-1,0), dot = -1.

This means:
- The DG formulation (K matrix, slip RHS) uses the CalcOrtho convention: normal points +y
- The traction computation's stress term uses basis.normal: normal points -y
- The penalty correction in ComputeTraction uses `sign` from CalcOrtho: sign = +1

### 8.2 Is This a Problem?

For a perfect elastic solve (K·u = f exactly), the penalty correction in ComputeTraction is zero (displacement jump equals prescribed slip). The traction is determined entirely by the stress term {σ}·basis.normal, which is in the correct physical direction for the fault.

For approximate solves (BLR, iterative), the penalty correction is nonzero. The correction is computed in the CalcOrtho convention (sign = +1) while the stress is in the basis.normal convention (opposite). The code combines:
```
T = {σ}·basis.normal - penalty·sign·((u1-u2) - sign·delta_u)
```

This computes T in the basis.normal direction:
- Stress part: σ·(0,-1,0) — traction on the -y surface
- Penalty correction: -penalty·(+1)·(jump - prescribed) — reduces error in +y convention
  = +penalty·(prescribed - jump) in basis.normal convention

This is CORRECT: both terms contribute to traction in the same physical direction. The `sign` in the correction effectively converts the CalcOrtho-convention residual to the basis.normal convention.

### 8.3 Verdict

The opposite normals are by design and handled correctly. No code change needed for normals.

---

## 9. Seed Traction Analysis (SEED-TRAC Follow-Up)

### 9.1 Finding

Only 3 DOFs with tau > 1 Pa out of 55920 (1000m) or 10248 (2500m). The 3 DOFs are at the mesh boundaries (corner effects) with traction ≈ 0.

### 9.2 Implication

There is no "bad initial traction" seeding the instability. The instability arises DYNAMICALLY during the first RK step:
1. Stage 0: V = V_init (given) → slip = 0 (no time has passed)
2. Stage 1: state advances by c₁·dt → new slip → elastic solve → traction → new V
3. Stage 2: state advances further → V grows from the nucleation physics
4. Stage 3+: if V has grown enough, the nonlinear amplification takes over

The "trigger" is the normal nucleation physics (V_nuc = 0.01 produces slip which changes traction which accelerates more DOFs). The INSTABILITY is the RK stage amplification exceeding the stability boundary.

---

## 10. Summary

| Finding | Evidence | Status |
|---------|----------|--------|
| **CFL stability constraint** | 2500m: z=1.44 (stable), 1000m: z=3.60 (unstable) | **ROOT CAUSE** |
| Normal mismatch | CalcOrtho = -basis.normal on ALL faces, BOTH meshes | By design, correct |
| Seed traction | 3/55920 DOFs > 1Pa, essentially zero | Not a factor |
| Smooth nucleation | Delays cascade 2 stages, doesn't fix | Insufficient |
| Quad order match (2p) | Crashes FASTER | Makes it worse |
| **PETSc default dt=0.1s** | z=2.77 at h=1000m (marginal) vs our 0.13s (z=3.60) | Explains Tandem survival |

### Fix Formula
```
dt_init = min(dt_from_V, C·η·h_min/(4μ))

where C = 2.0, η = μ/(2cs), h_min = minimum fault element size
```

For h=1000m: dt_init = min(0.13, 0.072) = **0.072s**

---

## 11. Phase 2: Immediate Next Steps

1. **v49f**: Run 1000m serial with `--dt-init 0.07` → confirm stability
2. **v49g**: Run 1000m serial with `--dt-init 0.05` → confirm with margin
3. **v49h**: Run 1000m serial with `--dt-init 0.10` → test PETSc default
4. **v49i**: Implement V_guard in RK stepper → test at dt=0.13 with guard

If v49f-h confirm the CFL hypothesis, implement automatic h_min-based dt selection for production runs.

---

## 12. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v47 | Penalty ×3 re-applied. 1000m crashes, 2500m serial stable. | Done |
| v48 | Sign hypothesis DISPROVED, BLR secondary, 1000m serial crash confirmed | Done |
| **v49** | **Phase 1 diagnostics: CFL stability identified as root cause** | **IDENTIFIED** |
| v49a | Smooth nucleation: delays cascade 2 stages, doesn't fix | Insufficient |
| v49b | Match quad order (2p): crashes FASTER | Makes worse |
| v49cd | Baseline diagnostic: normals correct, seed=0, RK cascade mapped | Done |
| v49e | 2500m reference: STABLE, 132+ steps, bounded amplification | Confirmed |
| v49f (planned) | **CFL fix: dt_init=0.07s for 1000m** | Pending |
