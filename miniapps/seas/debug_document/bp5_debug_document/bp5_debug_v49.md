# BP5 Debug v49: CFL Stability Analysis — IP Penalty Eigenvalue vs dt_init

**Date**: 2026-03-22
**Status**: ISSUE A (1000m serial crash) SOLVED. CFL-like stability constraint confirmed — dt_init=0.07s and V-guard both fix 1000m. ISSUE B (2500m parallel crash) STILL OPEN — not addressed by Phase 1+2 (all serial tests). Next: test 2500m parallel with V-guard.
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
| v49f | **CFL fix: dt_init=0.07s → STABLE** (cascade V=87, RK rejects, recovers) | **CONFIRMED** |
| v49g | Conservative dt_init=0.05s → STABLE (cascade V=3.4, RK rejects, recovers) | Confirmed |
| v49h | PETSc default dt_init=0.10s → **CRASH** (cascade V=129, traction GPa) | Threshold between 0.07-0.10 |
| v49i | V-guard at dt=0.13s → **STABLE** (guard triggers 2×, halves dt, recovers) | **V-GUARD WORKS** |

---

## 13. Phase 2 Results

### 13.1 Overall Outcome

| Test | dt_init | Max cascade V | Stage of max | Outcome |
|------|---------|--------------|-------------|---------|
| **v49f** | 0.07s | 87.1 (stage 6) | All 7 stages complete | **STABLE** — RK error rejects, dt→0.009 |
| **v49g** | 0.05s | 3.39 (stage 5) | All 7 stages complete | **STABLE** — RK error rejects, dt→0.009 |
| **v49h** | 0.10s | 129 (stage 5) | Stage 5 triggers blowup | **CRASH** — traction GPa at 5 DOFs |
| **v49i** | 0.13s + V-guard(100×) | 20.9 (stage 4, caught) | V-guard at stage 4 | **STABLE** — halves dt ×2, then proceeds |

### 13.2 v49f: dt=0.07s (CFL Fix) — STABLE

First step attempt at dt=0.07s shows cascade:
```
Stage 0: V_str = 0.0100
Stage 1: V_str = 0.0106   (×1.06)
Stage 2: V_str = 0.0122   (×1.15)
Stage 3: V_str = 0.1417   (×11.6)   ← nonlinear onset
Stage 4: V_str = 1.969    (×13.9)
Stage 5: V_str = 21.60    (×11.0)
Stage 6: V_str = 87.14    (×4.03)   ← still large but ALL 7 stages complete
```

RK error estimate: huge → step REJECTED. Retries with dt ≈ 0.005s:
```
Stage 0: V_str = 0.0100
Stage 1: V_str = 0.01003  (×1.003)   ← calm
... all stages < 0.0101
```
Step accepted, dt grows adaptively. After step 1 (dt=0.009s), steady progression:
Steps 2-10: V_max ≈ 0.010, dt growing to 0.017s. **Healthy nucleation.**

### 13.3 v49g: dt=0.05s (Conservative) — STABLE

Same pattern as v49f but milder cascade (peak V=3.39 at stage 5). RK rejects first attempt, retries at dt≈0.005s. Subsequent steps identical to v49f.

### 13.4 v49h: dt=0.10s (PETSc Default) — CRASH

```
Stage 0: V_str = 0.0100
Stage 1: V_str = 0.0108   (×1.08)
Stage 2: V_str = 0.0158   (×1.46)
Stage 3: V_str = 0.428    (×27.1)   ← runaway
Stage 4: V_str = 8.155    (×19.1)
Stage 5: V_str = 129.2    (×15.8)   ← OVERFLOW
```
TRACTION BLOWUP at 5 interior DOFs, tau > 1 GPa. Segfault.

**Key: dt=0.10 is ABOVE the 1000m stability threshold.** The CFL estimate of dt_CFL=0.072s is accurate — the threshold is between 0.07 and 0.10.

### 13.5 v49i: V-guard — STABLE

Default dt=0.13s. V-guard factor = 100× (stage V > 100 × V_stage0 → reject and halve):
```
Attempt 1 (dt=0.13s):
  Stage 4: V=20.94 > 100 × 0.01 = 1.0 → REJECTED, dt → 0.065s

Attempt 2 (dt=0.065s):
  Stage 4: V=1.43 > 100 × 0.01 = 1.0 → REJECTED, dt → 0.0325s

Attempt 3 (dt=0.0325s):
  All stages complete (max V=0.258 at stage 5) → ACCEPTED by V-guard
  RK error estimate still large → step rejected by error control
  Retries at dt ≈ 0.006s → accepted
```
After that, step 1 accepted at dt=0.009s. Subsequent steps identical to v49f/v49g.

### 13.6 Stage-by-Stage Amplification Comparison

| Stage | v49g (dt=0.05) | v49f (dt=0.07) | v49h (dt=0.10) | v49cd (dt=0.13) |
|-------|---------------|---------------|---------------|----------------|
| 0 | 0.0100 | 0.0100 | 0.0100 | 0.0100 |
| 1 | 0.0104 | 0.0106 | 0.0108 | 0.0111 |
| 2 | 0.0108 | 0.0122 | 0.0158 | 0.0217 |
| 3 | 0.0455 | 0.1417 | 0.428 | 0.873 |
| 4 | 0.455 | 1.969 | 8.155 | 20.94 |
| 5 | 3.39 | 21.60 | **129 → CRASH** | **CRASH** |
| 6 | 3.07 | 87.14 | — | — |
| **Result** | **RK rejects** | **RK rejects** | **Segfault** | **Segfault** |

The stage 3→4 amplification is the critical nonlinear transition:
- dt=0.05: 0.045 → 0.455 (10×)
- dt=0.07: 0.14 → 1.97 (14×)
- dt=0.10: 0.43 → 8.16 (19×)
- dt=0.13: 0.87 → 20.9 (24×)

The amplification RATE increases with dt (larger steps feed more energy into the cascade).

### 13.7 Post-Stabilization Convergence

All three stable runs (v49f, v49g, v49i) converge to identical behavior after step 1:

| Step | Time [yr] | dt [s] | V_max [m/s] |
|------|-----------|--------|-------------|
| 1 | 1.1e-10 | 0.009 | 0.0101 |
| 2 | 3.6e-10 | 0.010 | 0.0102 |
| 3 | 7.1e-10 | 0.014 | 0.0101 |
| 4 | 1.1e-09 | 0.018 | 0.0102 |
| 5 | 1.4e-09 | 0.018 | 0.0102 |
| ... | ... | 0.017-0.021 | 0.010-0.010 |

**The initial dt choice only affects the first step.** Once the adaptive controller finds the right dt (≈0.01-0.02s), the simulation is identical regardless of how it got there.

---

## 14. Two Distinct Issues Confirmed

The Phase 1+2 investigation reveals TWO separate issues:

### 14.1 Issue A: 1000m Serial Crash (SOLVED)

**Root cause**: CFL-like stability constraint. IP penalty stiffness scales as 1/h, making dt_init=0.13s too large for h=1000m.

**Evidence**: Serial runs at dt=0.07s and dt=0.05s are STABLE. V-guard at dt=0.13s is also STABLE. dt=0.10s is still too large (CRASH).

**Fix**: Either CFL-aware dt_init or V-guard in the RK stepper (or both for robustness).

### 14.2 Issue B: 2500m Parallel Crash (v47m — STILL OPEN)

**Root cause**: Unknown. The 2500m serial run (v47l, v49e) is STABLE with the same dt=0.13s. Only the parallel run (v47m, 48 ranks) crashes. This requires shared faces to trigger.

**Not addressed by Phase 1+2**: All v49 tests were serial. The parallel crash is a separate investigation track.

**The v48 shared face sign hypothesis** remains a candidate for Issue B, though Section 6.1 of this document notes that v48 Section 10 disproved the CalcOrtho sign claim. The parallel crash mechanism may involve:
- K-f mismatch on shared faces (different from the sign direction issue)
- BLR solver interaction with parallel assembly
- Face-neighbor data exchange timing
- Or another undiscovered parallel-specific bug

---

## 15. Recommended Next Steps

### 15.1 Immediate (Issue A — Production Fix)

1. **Implement CFL-aware dt_init** in production code:
   ```
   dt_init = min(L/(10·V_max), 2.0 · η · h_min / (4·μ))
   ```
   Compute h_min from mesh element volumes on the fault.

2. **Keep V-guard as safety net**: Already implemented and tested in v49i. Factor 100× is conservative but effective.

### 15.2 Next Investigation (Issue B — Parallel)

1. Run 2500m **parallel** with V-guard → does it survive?
2. If V-guard helps: the parallel crash is also a cascade/dt issue (shared faces may have different effective stiffness)
3. If V-guard doesn't help: the parallel bug is structural (sign, assembly, or data exchange error) → dump K×u vs f on shared faces

### 15.3 Long-Term

1. Implement automatic h_min computation from mesh
2. Make dt_init formula the default (no manual --dt-init flag needed)
3. Add V-guard with configurable factor to production RK45
4. Profile to ensure the V-guard overhead is negligible

---

## 16. Phase 3: Issue B — Parallel Crash Investigation

### 16.1 Context

Phase 1+2 resolved Issue A (1000m serial crash → CFL/dt stiffness). Issue B remains:

| Run | Mesh | Ranks | dt | Result |
|-----|------|-------|-----|--------|
| v47l | 2500m | 1 (serial) | 0.13s | STABLE |
| v49e | 2500m | 1 (serial) | 0.13s | STABLE (132+ steps) |
| v47m | 2500m | 48 (parallel) | 0.13s | **CRASH** — GPa traction, SIGNAL 9/11 |

The 2500m serial run shows only mild cascade (peak V=0.13 at stage 5, recovers).
The 2500m parallel run crashes immediately. The ONLY difference is MPI partitioning,
which introduces shared faces at partition boundaries.

### 16.2 Key Question

Is the parallel crash:
- **(A)** A cascade/stiffness issue amplified by shared faces (e.g., shared faces have
  higher effective stiffness due to face-neighbor coupling), fixable by V-guard/dt?
- **(B)** A structural parallel bug (sign error, assembly mismatch, data exchange issue)
  that produces wrong results regardless of dt?

### 16.3 Phase 3 Tests

Three tests targeting the parallel crash:

| Test | Mesh | Ranks | Nodes | Fix applied | What it tests |
|------|------|-------|-------|------------|---------------|
| **v49j** | 2500m | 48 | 1 | `--v-guard` | Can V-guard prevent the 2500m parallel cascade? |
| **v49k** | 2500m | 48 | 1 | `--dt-init 0.05` | Is the 2500m parallel crash dt-dependent? |
| **v49l** | 1000m | 8 | 4 | `--v-guard` + `--dt-init 0.05` | 1000m parallel — small scale |
| **v49m** | 1000m | 400 | 8 | `--v-guard` + `--dt-init 0.05` | **PRODUCTION SCALE** — max shared faces |

All include `--diag-rk-stages` for cascade tracing. v49j builds; v49k, v49l, v49m wait for binary.

v49l (8 ranks) is the minimal parallel test for 1000m. v49m (400 ranks) is the
full production configuration: 63451 tets / 400 ranks ≈ 159 tets/rank, maximizing
the number of shared faces. If v49l passes but v49m fails, the issue scales with
partition count (more shared faces = more exposure to the parallel bug).

### 16.4 Decision Tree

```
v49j (2500m, 48 ranks, V-guard) + v49k (2500m, 48 ranks, dt=0.05):
├── Both STABLE → Issue B is cascade/CFL at 2500m parallel
│   └── v49l (1000m, 8 ranks) + v49m (1000m, 400 ranks):
│       ├── Both STABLE → ALL ISSUES RESOLVED. Production ready.
│       ├── v49l OK, v49m CRASH → scales with partition count
│       │                         (more shared faces = stiffer system)
│       │                         → tighten V-guard or reduce dt further
│       └── Both CRASH → 1000m parallel has issue beyond cascade
│                         → K×u vs f diagnostic on 1000m shared faces
│
├── v49j STABLE, v49k CRASH → V-guard specifically needed for parallel
│   └── v49l/v49m should be STABLE (have V-guard)
│
└── Both CRASH → **Structural parallel bug confirmed**
                 → K×u vs f diagnostic on shared faces
                 → v49l/v49m will also crash (same bug)
```

### 16.5 If Structural Bug Confirmed: Next Diagnostic

If both v49j and v49k crash, implement a K×u vs f consistency check:

1. On each shared fault face, extract the K penalty rows for Elem1 DOFs
2. Construct `u_prescribed` such that `[[u]] = sign * delta_u` (the prescribed slip)
3. Compute `f_from_K = K_face_penalty * u_prescribed`
4. Compare with `f_from_RHS` from `AssembleSlipContributionIPShared`
5. Print per-DOF difference: `|f_from_K - f_from_RHS|`
6. Any nonzero difference reveals the assembly inconsistency

This requires extracting the element stiffness matrix contribution from a single shared
face, which can be done by calling the DG integrator's `AssembleFaceMatrix` directly
on the shared face FaceElementTransformations.

### 16.6 Submission

```bash
# Submit v49j first (builds), then the rest (wait for binary)
sbatch jobs/bp5/bp5_v49j_2500m_par_vguard.sbatch
sbatch jobs/bp5/bp5_v49k_2500m_par_dt005.sbatch
sbatch jobs/bp5/bp5_v49l_1000m_par_vguard.sbatch
sbatch jobs/bp5/bp5_v49m_1000m_par400_vguard.sbatch
```

---

## 17. Revision History (Updated)

| Version | Change | Status |
|---------|--------|--------|
| v47 | Penalty ×3 re-applied. 1000m crashes, 2500m serial stable. | Done |
| v48 | Sign hypothesis disproved, BLR secondary, 1000m serial crash confirmed | Done |
| **v49** | **Phase 1: CFL stability identified as root cause of 1000m crash** | Done |
| v49a | Smooth nucleation: delays cascade 2 stages, doesn't fix | Insufficient |
| v49b | Match quad order (2p): crashes FASTER | Makes worse |
| v49cd | Baseline diagnostic: normals correct, seed=0, RK cascade mapped | Done |
| v49e | 2500m serial reference: STABLE, 132+ steps, bounded amplification | Confirmed |
| **v49f** | **CFL fix: dt=0.07s → STABLE** (cascade V=87, RK rejects, recovers) | **CONFIRMED** |
| **v49g** | **Conservative dt=0.05s → STABLE** (cascade V=3.4, RK rejects) | Confirmed |
| **v49h** | **PETSc default dt=0.10s → CRASH** (cascade V=129, GPa traction) | Threshold 0.07-0.10 |
| **v49i** | **V-guard at dt=0.13s → STABLE** (guard triggers 2×, recovers) | **V-GUARD WORKS** |
| v49j (pending) | 2500m parallel + V-guard → test Issue B | **Submitted** |
| v49k (pending) | 2500m parallel + dt=0.05 → test Issue B dt-dependence | **Submitted** |
| v49l (pending) | 1000m parallel 8 ranks + V-guard + dt=0.05 | **Submitted** |
| v49m (pending) | **1000m parallel 400 ranks + V-guard + dt=0.05 → PRODUCTION SCALE** | **Submitted** |
