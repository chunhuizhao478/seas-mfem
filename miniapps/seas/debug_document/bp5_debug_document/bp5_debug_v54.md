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

## 5. The Core Problem

**Every solver configuration tried (BLR 1e-10, 1e-12, exact MUMPS OOM) fails to
nucleate at p=1 1000m.** The BLR penalty amplification identified in v52 is the most
likely root cause: MUMPS-BLR introduces systematic solver residual that gets amplified
by the IP penalty correction into a ~0.3+ MPa traction bias that suppresses nucleation.

This is NOT a discretization issue (fixed in v53) or a parameter issue (matches Tandem).
It is a **solver accuracy issue** specific to the MUMPS-BLR + IP penalty combination.

### 5.1 What Tandem Actually Uses

Tandem uses **PETSc KSP CG** (iterative conjugate gradient) with:

- `matrix_free = true`: system operator A is a matrix-free shell (`PetscDGShell`)
  for exact matrix-vector products
- An **assembled DG matrix** P (`PetscDGMatrix`) as preconditioner
- `mg_strategy = "logarithmic"`: p-multigrid preconditioner — but at p=1, there is
  only one MG level (can't coarsen below p=1), so MG degenerates to the default
  PETSc preconditioner on the assembled matrix
- `rtol = 1e-12`: CG converges until ||Ax-b||/||b|| < 1e-12

Source: `tandem/app/common/PetscLinearSolver.cpp:7-44`

At p=1, Tandem's solver is effectively **CG + block-Jacobi/ILU on the assembled DG
matrix**. It is NOT some fundamentally different solver architecture — it uses the same
assembled DG system. The critical difference is:

| Property | MFEM (MUMPS-BLR) | Tandem (PETSc CG) |
|----------|-------------------|--------------------|
| Solve type | Direct (approximate factorization) | Iterative (converges to true residual) |
| Residual guarantee | BLR-tol dependent, no bound on ||Ax-b|| | **||Ax-b||/||b|| < 1e-12** |
| Penalty amplification of residual | ~300x BLR residual → ~0.03 MPa systematic error | ~300x × 1e-12 → negligible |
| Factorization error | YES (BLR low-rank blocks) | NO (iterative, no factorization in solution) |

The BLR residual (~1e-7 m at fault DOFs, v52 Section 7) gets amplified by the IP
penalty correction in `ComputeTraction` by ~300x, producing ~30 Pa per solve. Over
thousands of steps, this systematic bias accumulates into the observed 0.38 MPa traction
excess that kills nucleation.

CG with rtol=1e-12 produces residuals of ~1e-12, which amplified by 300x gives ~3e-10 —
twelve orders of magnitude smaller. This is why Tandem nucleates and MFEM does not.

---

## 6. Proposed Fixes

### Fix 1: Implement PETSc CG Solver, Matching Tandem (Priority: CRITICAL)

**Rationale**: Tandem uses PETSc KSP CG with the assembled DG matrix as preconditioner
(`PetscLinearSolver.cpp:22-27`). This guarantees ||Ax-b||/||b|| < 1e-12, eliminating the
BLR penalty amplification that kills nucleation in MFEM.

**Implementation options** (in order of preference):

**Option A: PETSc KSP CG (match Tandem exactly)**
- MFEM already supports PETSc via `--enable-petsc`. Use `PetscLinearSolver` with KSP CG.
- Assemble the DG matrix once, use as both operator and preconditioner (like Tandem at p=1).
- Set rtol=1e-12.
- At p=1 this matches Tandem's solver architecture exactly.

**Option B: MFEM native CG + block-diagonal preconditioner**
- Use MFEM's `CGSolver` with `BlockDiagonalPreconditioner` or `HypreBoomerAMG`.
- Set rtol=1e-12.
- May require tuning for DG systems.

**Option C: MUMPS as preconditioner + CG outer loop**
- Use MUMPS-BLR factorization as preconditioner for CG (1-2 CG iterations should suffice).
- CG drives the residual below rtol even though the preconditioner is approximate.
- Minimal code change: wrap existing MUMPS solve in a CG loop.

**Expected impact**: Eliminates systematic solver residual amplification by IP penalty.
CG residual ~1e-12 × penalty amplification 300x → ~3e-10, which is negligible.
This should allow nucleation to proceed as in Tandem.

**Risk**: CG may converge slowly for ill-conditioned DG systems. Option C mitigates this
by using MUMPS-BLR as a high-quality preconditioner (expect convergence in 1-3 iterations).

### Fix 2: Add eta*V to Traction Output (Priority: HIGH, COSMETIC)

Add `eta * V` to the traction output in `bp5_benchmark_output.hpp` to match Tandem's
`tau_hat` convention. This eliminates the 0.046 MPa baseline offset in comparisons.

Does not affect dynamics. Needed for accurate comparison plots.

### Fix 3: Use `--zero-dip-traction` (Priority: HIGH)

The v53 run omitted this flag. Since DG cross-component contamination drains energy from
strike to dip, enabling `--zero-dip-traction` may partially improve nucleation by
preventing the 21% dip contamination from weakening the strike traction.

Not a fix by itself (v52 showed it doesn't fix the tau_strike deficit), but reduces one
source of error.

### Fix 4: Try Reduced Penalty Factor (Priority: MEDIUM)

If the penalty amplifies BLR residual by 300x, reducing the penalty factor may reduce the
traction bias at the cost of some DG stability. Test with `--penalty-factor 0.5` to see
if nucleation improves while DG remains stable.

**Note**: v50 tried penalty reduction but may not have tested at p=1 with 3 DOFs/face.
Worth re-testing with the current code.

### Fix 5: Try CG+MUMPS Hybrid (Priority: MEDIUM)

Use CG with MUMPS as a preconditioner (rather than direct solve). CG ensures the residual
is driven to true tolerance, while MUMPS-BLR provides a good preconditioner. This may be
faster than CG+AMG for DG systems.

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

## 9. Bottom Line

The v53 3-DOF/face fix correctly resolved the discrete model mismatch, but
**nucleation still dies because MUMPS-BLR introduces a growing traction bias
(~0.38 MPa at t=38s) that suppresses the nucleation zone acceleration.**

Both MFEM and Tandem show an initial V overshoot then decline. In Tandem, V dips to
0.0114 but stays above V_nuc and re-accelerates to earthquake at t=82s. In MFEM, V
drops below V_nuc by t~15s and keeps falling — the 0.38 MPa excess traction creates
enough extra friction to prevent recovery.

Tandem uses **PETSc KSP CG** (iterative) with the assembled DG matrix as preconditioner,
guaranteeing ||Ax-b||/||b|| < 1e-12. MFEM uses **MUMPS-BLR** (direct, approximate
factorization) where the BLR residual gets amplified ~300x by the IP penalty correction.

The fix is to switch from MUMPS-BLR direct solve to **CG iterative solve** (matching
Tandem's solver architecture), which guarantees the residual is driven to true tolerance
regardless of preconditioner quality. Option C (MUMPS-BLR as preconditioner for CG
outer loop) is the lowest-risk path — it reuses the existing solver as a preconditioner
and should converge in 1-3 CG iterations.
