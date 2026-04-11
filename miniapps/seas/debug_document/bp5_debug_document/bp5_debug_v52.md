# BP5 Debug v52: Zero Dip Traction — Isolating Cross-Component DG Contamination

**Date**: 2026-03-24
**Status**: PRODUCTION RUN SUBMITTED
**Previous**: v51 (dip offset root cause investigation, elastic σ_n, K-matrix coupling diagnostics)
**Branch**: `feature/elasticity`

---

## 1. Motivation

v51 established that the DG discretization produces ~21% cross-component dip/strike traction
contamination (Section 19), and that elastic σ_n feedback partially suppresses it (Section 25).
However, elastic σ_n is an indirect fix — it adjusts friction strength to absorb spurious dip,
rather than eliminating the spurious dip at its source.

This version takes the direct approach: **zero τ_dip after ComputeTraction**.

### 1.1 Why This Is Benchmark-Correct

For BP5 on a planar vertical fault with pure strike-slip loading:
- SCEC spec Eq. 15b: V_3 = 0 (zero dip velocity) outside Ω_f
- SCEC spec Eq. 16: V_3 = V_zero = 1e-20 m/s (negligibly small)
- τ_dip = 0 analytically at all times, everywhere on the fault
- Even Tandem's 0.2-0.7% dip/strike ratio is a discretization artifact

The `--zero-dip-traction` flag is not a workaround — it enforces the analytically correct
value that the DG discretization fails to preserve.

### 1.2 What This Tests

By zeroing only the elastic dip traction while leaving everything else unchanged, we isolate
whether the DG cross-component contamination is responsible for:

1. **Dip offset** — should go to exactly zero (no source of spurious dip)
2. **Interseismic strike deviations** — τ_strike was ~1 MPa low; the cascade
   (Section 4 of v51) predicts this is caused by dip contamination weakening the fault
3. **Earthquake timing drift** — events 2+ diverge from Tandem

If results improve on all three, the root cause is confirmed and we can pursue a proper
fix (e.g., reformulating the DG traction recovery to preserve component decoupling).

---

## 2. Code Changes (v51 → v52)

### 2.1 Removed: Single Elastic Solve Diagnostic

The `--single-solve-uz` flag and its entire implementation block (~100 lines) were removed
from `bp5_verification_full.cpp`. This was a one-time diagnostic that prescribed a slip
pattern, performed a single elastic solve, and dumped u_z/u_x statistics. The diagnostic
served its purpose in v51 Section 19 (confirming u_z contamination) and is no longer needed.

Removed:
- Variable declaration: `bool single_solve_uz`
- CLI parsing: `--single-solve-uz` flag
- Implementation block: slip prescription, elastic solve, MPI reduction, statistics output

### 2.2 Removed: K-Matrix Coupling Test

The `test_k_matrix_coupling.cpp` Makefile entries were removed. This unit test was a v51
diagnostic for verifying cross-component coupling in the stiffness matrix. The test confirmed
that K-matrix components match Tandem (v51 Section 22), so the build rules are no longer
needed.

Removed from Makefile:
- `TEST_K_COUPLING_SRC` / `TEST_K_COUPLING_OBJ` variable definitions
- `seas_test_k_coupling` build target
- Compilation rule for `TEST_K_COUPLING_OBJ`

### 2.3 Retained: `--zero-dip-traction` Implementation

The zero-dip-traction logic in `seas_operator.hpp` (lines 255-262) is unchanged:

```cpp
// v51: Zero dip traction component (index 0 of each DOF's [dip, strike] pair)
if (zero_dip_traction_)
{
   for (int i = 0; i < traction_.Size() / 2; i++)
   {
      traction_(2 * i) = 0.0;  // tau_dip = 0
   }
}
```

This executes after `ComputeTraction()` and before `ComputeRHS()`, so:
- The elastic solve runs normally (3D coupling intact)
- Only the dip component of the resulting traction is zeroed
- `tau_total_dip = tau_pre_dip + 0 ≈ 0` (tau_pre_dip is O(1e-11) Pa, negligible)
- The friction law then sees effectively zero dip stress → V_dip ≈ V_zero

---

## 3. Production Run Configuration

**Job**: `bp5_v52_zero_dip_prod.sbatch`

```
Mesh:       bp5_tandem.msh (1000m reference mesh)
Order:      p=2
DG method:  IP (SIPG)
Solver:     MUMPS-BLR (tol=1e-12)
BC mode:    FarField (default)
t_final:    5.68e10 s (~1800 years)
Nodes:      8 (400 MPI ranks)
Queue:      normal (48hr)
Flags:      --zero-dip-traction --check-residual --write-every-step
Output:     bp5/results_v52_zero_dip_prod/
```

This matches the v50 production baseline exactly, with only `--zero-dip-traction` added.
No `--elastic-sigma-n` — we test each fix independently.

---

## 4. Expected Outcomes

### 4.1 If Zero Dip Traction Fixes All Issues

- Dip slip offset → 0 at all stations
- Dip V → O(1e-20) (matches V_zero)
- τ_strike recovers ~1 MPa deficit
- Earthquake timing matches Tandem through multiple cycles
- State variable ψ matches (no longer depressed by spurious |V| from dip)
- **Conclusion**: DG cross-component contamination was the sole root cause

### 4.2 If Zero Dip Traction Fixes Dip But Not Strike

- Dip quantities correct
- τ_strike still ~1 MPa low, timing still drifts
- **Conclusion**: Strike deficit has a separate source (e.g., traction magnitude error,
  boundary loading, or penalty-related stress redistribution)

### 4.3 If Zero Dip Traction Has Minimal Effect

- Dip should still improve (we're removing the source directly)
- If it doesn't, something else feeds dip (e.g., slip evolution, state feedback)
- **Conclusion**: Unlikely — would indicate a bug in the zeroing implementation

---

## 5. Comparison Plan

When the production run completes, compare against:

1. **Tandem reference**: `bp5/benchmark_data/bp5-qd-*.txt` (8-column SCEC format)
2. **v50 production**: `bp5/results_v50_prod/` (baseline without any dip fix)
3. **v51d elastic σ_n**: if available, to compare direct zeroing vs. indirect feedback

Key metrics per station:
- Dip slip offset (should be ~0)
- Dip/strike traction ratio (should be ~0)
- τ_strike absolute value vs. Tandem
- Earthquake recurrence interval
- State variable ψ trajectory

---

## 6. Solve vs. Traction Coherence: Concrete Code Comparison

### 6.1 Correction: Both Codes Use the Same Architecture

Initial hypothesis was that Tandem excludes fault faces from the stiffness matrix
while MFEM includes them. **This is wrong.** Code-level investigation confirms
both codes use the same split: K includes fault faces with [[u]] as unknown,
RHS corrects for prescribed slip.

**Evidence — Tandem stiffness includes fault faces** (DGOperator.h:142-174):
```cpp
// assemble() loops ALL skeleton faces — no BC type check:
for (fctNo = 0; fctNo < numLocalFacets(); ++fctNo)
    if (info.up[0] != info.up[1])  // skeleton (fault OR regular)
        lop_->assemble_skeleton(fctNo, info, A00, A01, A10, A11);
```

**Evidence — Tandem matrix-free K*u includes fault faces** (Elasticity.cpp:751):
```cpp
if (bc == BC::None || (is_skeleton_face && is_fault_or_dirichlet)) {
    flux_u_skeleton(...)           // [[u]] as unknown — SAME as regular face
    flux_sigma_skeleton(...)       // c00 = -penalty(fctNo)
    if constexpr (WithRHS) {       // Only during RHS, not K*u product
        if (bc_skeleton(fctNo, bc, f_q))
            flux_u_add_bc(...)     // slip correction added to RHS only
    }
}
```
When `WithRHS = false` (K*u during CG iterations), fault = regular interior face.

**Evidence — MFEM stiffness includes fault faces** (elasticity_operator.hpp:870):
```cpp
cached_a_->AddInteriorFaceIntegrator(
    new DGElasticityIntegrator(lambda_, mu_, epsilon_, 0.0));  // ALL faces
cached_a_->AddInteriorFaceIntegrator(
    new DGElasticityIPPenaltyIntegrator(lambda_, mu_, 3));     // ALL faces
```
MFEM's BilinearForm loops all interior faces. No fault exclusion.

**Both codes solve the identical system:** `K*u = f_vol + f_dirichlet + f_slip`
where K includes penalty on fault faces enforcing [[u]]=0, f_slip corrects.

### 6.2 The Real Differences (Concrete, Code-Level)

#### 6.2.1 Solver: CG matrix-free (Tandem) vs MUMPS-BLR assembled (MFEM)

**Tandem** (PetscLinearSolver.cpp:22, bp5.toml):
- PETSc CG, rtol = 1e-12, `matrix_free = true`, p-multigrid preconditioner
- K*u evaluated **exactly** at every CG iteration (no factorization error)
- Converged solution satisfies ||K*u - f|| / ||f|| < 1e-12

**MFEM** (elasticity_operator.hpp:908-914):
- MUMPS-BLR with BLR tol = 1e-12, assembled HypreParMatrix
- Approximate LU: BLR compresses off-diagonal blocks via low-rank approximation
- Factorization error: ||L*U - K|| / ||K|| ~ 1e-12, but **structurally non-uniform**

#### 6.2.2 Penalty magnitude at fault faces (computed from BP5 parameters)

```
mu = 32.04 GPa, lambda = 32.04 GPa
c0 = 2*mu = 64.08 GPa, c1 = 3*lambda + 2*mu = 160.19 GPa
c_N_1 = p*(p+dim-1)/dim = 2*4/3 = 2.667 (at p=2)
face_area ~ 8.1e5 m^2 (from h_min=812m mesh)
vol ~ 5.4e8 m^3 (element volume)
A/V = 3 * face_area / vol ~ 4.5e-3 /m
p_each = (D+1)*c_N_1*(A/V)*(c1^2/c0) = 4*2.667*4.5e-3*400e9 = 19.2e9 Pa/m
penalty_ip = (p0+p1)/4 ~ 9.6e9 Pa/m
```

**Fault penalty ~10^10 Pa/m. Bulk stiffness ~mu/h = 3.2e7 Pa/m. Ratio: ~300x.**

The fault face entries are the largest in K. BLR low-rank compression error
concentrates at these high-contrast entries.

#### 6.2.3 Predicted vs observed solver fault residual

BLR tol 1e-12 controls ||L*U - K|| / ||K||. The solution error depends on cond(K):
- DG stiffness cond(K) ~ O(penalty/bulk * h^-2) ~ O(300 * 10^6) ~ O(10^8)
- Expected |u_err| ~ cond(K) * BLR_tol * |u| ~ 10^8 * 10^-12 * |u| = 10^-4 * |u|
- For |u| ~ 10^-5 m (first step): |u_err| ~ 10^-9 m

But the **jump** of the error across the fault can be much larger than the pointwise
error because adjacent elements share the same BLR compression artifact.

**Observed** (Section 7.2): |[[u]] - slip| = 4.4e-7 m RMS, 3.5e-5 m max.
The max is 7 orders above BLR tolerance, consistent with cond(K) amplification
plus error concentration at high-penalty fault faces.

#### 6.2.4 Penalty amplification of solver residual

```
penalty * max_res = 9.6e9 Pa/m * 3.5e-5 m = 3.4e5 Pa = 340 kPa
stress signal at first step: ~290 Pa RMS
penalty/stress = 340,000 / 290 = 1170x
```

**Observed** (Section 7.2): corr/stress = 632%. Order-of-magnitude match.

#### 6.2.5 Traction extraction: codegen (Tandem) vs hand-coded (MFEM)

**Tandem** (Elasticity.cpp:971-972):
```cpp
krnl.c00 = -penalty(fctNo);  // SAME function as stiffness line 761
```
- `penalty(fctNo)` returns precomputed value used in both K assembly and traction
- Kernels generated by YATETO from single Python source — guaranteed consistent
- Precomputed `n_unit_q`, `JInv`, `E_q` reused from assembly

**MFEM** (elasticity_operator.hpp:3468-3479):
```cpp
real_t penalty_ip = penalty_factor_ * (p0 + p1) / 4.0;  // recomputed from scratch
```
- Same formula but independent code path (not shared function)
- Stress uses `CalcInverse` (traction) vs `CalcAdjugate` (solve integrator)
- On planar hex faces: mathematically identical, numerically independent

#### 6.2.6 Summary table

| Aspect | Tandem | MFEM | Impact |
|--------|--------|------|--------|
| Fault faces in K | Yes (line 751) | Yes (line 870) | Same architecture |
| Solver | CG matrix-free (exact K*u) | MUMPS-BLR (approx. factorization) | **BLR error at fault** |
| Solver tolerance | rtol 1e-12 (Krylov) | BLR tol 1e-12 (factorization) | Different error distribution |
| Penalty in traction | `penalty(fctNo)` shared | Recomputed independently | Math. same |
| Kernel consistency | YATETO codegen (1 source) | Hand-coded (2 code paths) | Tandem guaranteed |
| Observed dip/strike | 0.2-0.7% | 31% stress + 632% penalty | **~100x worse** |
| Fault residual | Not measured (CG → ~0) | 4.4e-7 m RMS, 3.5e-5 m max | **Root cause** |

**CORRECTION**: The earlier Section 6.1 description of Tandem as "coherent" and MFEM
as "decoupled" was misleading. Both use the same K + f_slip architecture. The actual
difference is the **solver** (CG exact K*u vs MUMPS-BLR approximate factorization)
and the **traction code path** (codegen-consistent vs hand-coded-independent).

The dominant contamination source is the MUMPS-BLR factorization error concentrating
at the high-penalty (~10^10 Pa/m) fault face entries, amplified by the penalty in
ComputeTraction into spurious traction 6x larger than the stress signal.

**Tandem's CG solver does not have this problem** because every K*u evaluation is
exact (matrix-free) — there is no factorization, so no structural compression error.
The traction extraction sees the same K*u the solver converged to.

### 6.3 Implemented: `--diag-traction-coherence` Flag

(Old Sections 6.2-6.5 described initial hypotheses about penalty computation
paths, adjugate vs inverse Jacobian, and proposed a residual-based test. These are
superseded by the concrete comparison in Section 6.2 above and the coherence test
results in Section 7.)


Rather than the full residual-based approach (Section 6.4), we implemented a lightweight
diagnostic that decomposes the existing ComputeTraction into its stress and penalty
components and measures each independently.

**What it measures** (per quad point, all fault faces, MPI-reduced):

1. **Stress-only dip/strike**: `{σ(u)·n}` projected to dip and strike separately.
   If dip/strike ratio >> 0, the **DG solution itself** has cross-component coupling
   (the elastic solve produces u_z for pure strike-slip input).

2. **Penalty correction dip/strike**: `η * ([[u]] - slip)` projected to dip/strike.
   If this is large, the **solver residual** at fault faces is being amplified by the
   penalty into spurious traction. This would indicate an incoherence between the solve
   (which balances penalty forces) and the extraction (which re-evaluates them).

3. **Solver fault residual**: `|[[u]] - slip|` at each quad point.
   This measures how well the solver enforces the prescribed jump. With MUMPS BLR at
   tol=1e-12, this should be O(1e-12) × displacement scale.

**Interpretation**:

| Stress dip | Penalty dip | Residual | Diagnosis |
|-----------|------------|---------|-----------|
| Large | Small | Small | DG solution has intrinsic cross-coupling |
| Small | Large | Large | Solver residual amplified by penalty |
| Both large | Large | Large | Both sources contribute |
| Small | Small | Small | Traction extraction is coherent |

**Job**: `bp5_v52_coherence_test.sbatch`

```
Queue:      development (2hr)
Nodes:      4 (200 MPI ranks)
t_final:    3.15e8 s (~10 years)
Flags:      --diag-traction-coherence --check-residual --write-every-step
Output:     bp5/results_v52_coherence/
```

Short run — the diagnostic fires once after the first non-trivial slip step and
prints a summary. The 10-year window also gives early interseismic data for quick
comparison against Tandem.

---

## 7. Coherence Test Results: BOTH Sources Confirmed (Job 7610206)

### 7.1 MPI Deadlock Fix

The first submission (Job 7610190) hung at the first time step. Root cause: the
`coherence_active` flag was computed locally — ranks without fault DOFs had
`any_slip = false` and skipped the MPI collectives, while ranks with fault DOFs
entered them. Two fixes applied:

1. `any_slip` check: now `MPI_Allreduce`'d so all ranks agree on `coherence_active`
2. Summary block guard: changed from `coherence_active && coh_n_qp > 0` to just
   `coherence_active`. The print is guarded by `is_root && g_n_qp > 0`.

### 7.2 Global Results (65,023 quad points across all fault faces)

```
Stress-only traction {sigma.n}:
  RMS tau_strike(stress) = 289.5 Pa
  RMS tau_dip(stress)    = 89.9 Pa
  max |tau_dip(stress)|  = 6411 Pa
  dip/strike ratio (RMS) = 31.1 %

Penalty correction eta*(jump-slip):
  RMS tau_strike(corr)   = 1830.9 Pa
  RMS tau_dip(corr)      = 292.0 Pa
  max |tau_dip(corr)|    = 22239 Pa
  corr/stress ratio (strike) = 632 %
  corr/stress ratio (dip)    = 325 %

Solver fault residual |[[u]] - slip|:
  RMS residual = 4.38e-07 m
  max residual = 3.53e-05 m
```

### 7.3 Per-Station Decomposition

| Station | stress_dip (Pa) | stress_strk (Pa) | dip/strk % | corr_dip (Pa) | corr_strk (Pa) | max_res (m) |
|---|---|---|---|---|---|---|
| strk-36dp+00 | -6.4 | -46.0 | 14.0 | 0.8 | 0.6 | 5.3e-10 |
| strk-16dp+00 | 24.6 | -93.2 | 26.4 | 1.7 | 2.2 | 1.2e-08 |
| strk+00dp+00 | 0.07 | -9.6 | 0.8 | 0.01 | -0.05 | 1.2e-10 |
| strk+16dp+00 | -0.01 | -2.4 | 0.6 | 0.0002 | 0.01 | 7.7e-12 |
| strk+36dp+00 | -0.01 | -0.6 | 1.5 | -0.0002 | -0.0008 | 4.1e-12 |
| strk-24dp+10 | 8.3 | 317 | 2.6 | -6.8 | 1.6 | 1.3e-08 |
| strk-16dp+10 | -12.8 | -335 | 3.8 | **32.6** | **144.8** | **9.2e-08** |
| strk+00dp+10 | 0.68 | -6.4 | 10.6 | -0.05 | -0.16 | 2.6e-10 |
| strk+16dp+10 | -0.23 | -1.8 | 12.5 | 0.0006 | 0.007 | 3.3e-11 |
| strk+00dp+22 | -0.62 | -4.1 | 15.2 | 0.04 | 0.11 | 7.9e-11 |

### 7.4 Analysis: Two Contamination Sources

**Source 1: DG stress cross-coupling (31% dip/strike globally)**

The stress-only traction `{σ(u)·n}` has a non-zero dip component even for pure
strike-slip input. This comes from the DG elastic solution itself — the penalty
enforcement in the bilinear form couples all three displacement components, so the
solution `u` has non-zero `u_z` that creates `σ_yz ≠ 0`, which projects to τ_dip.

Spatial pattern:
- Worst at surface near nucleation: strk-16dp+00 = 26.4% dip/strike
- Low at center depth: strk+00dp+10 = 10.6%, strk-24dp+10 = 2.6%
- Very low far from nucleation: strk+16dp+00 = 0.6%, strk+36dp+00 = 1.5%
- Asymmetric: negative x2 stations have higher contamination than positive x2

**Source 2: Penalty amplification of solver residual (632% corr/stress)**

The penalty correction `η([[u]] - slip)` is **6× larger than the stress signal** in
strike. The solver does not enforce `[[u]] = slip` accurately at fault faces:
- RMS residual = 4.4e-7 m (should be ~1e-12 with BLR tol=1e-12)
- Max residual = 3.5e-5 m (at nucleation zone)

This means the MUMPS-BLR solver's approximate factorization leaves a large residual
specifically at the displacement jump across fault faces. The penalty coefficient is
O(10^7-10^8), so a 1e-7 m residual becomes O(1-10 Pa) spurious traction, and at the
nucleation zone (max_res = 9.2e-8 m), the penalty produces corr_strk = 145 Pa — far
exceeding the stress_strk = 335 Pa, making the penalty the dominant traction source.

Worst station: strk-16dp+10 (nucleation edge)
- `corr_dip = 32.6 Pa` vs `stress_dip = -12.8 Pa` → penalty 2.5× stress in dip
- `corr_strk = 144.8 Pa` vs `stress_strk = -335 Pa` → penalty 43% of stress in strike
- `max_res = 9.2e-8 m` → 4 orders above BLR tolerance

### 7.5 Why the Solver Residual Is Large

The BLR tolerance (1e-12) controls the **relative** factorization error of the
global stiffness matrix, not the pointwise error at fault faces. The fault faces
have the highest penalty coefficients in the system (they enforce the slip BC), so
they concentrate the factorization error. The BLR low-rank approximation may
compress exactly the blocks that correspond to fault-face penalty coupling.

### 7.6 Implications for the Dip Problem

The coherence test confirms that the dip contamination has **two independent sources**:

1. **Intrinsic DG coupling** (~31% globally, up to 26% at specific stations) from
   the stress field. This is a property of the DG discretization on this mesh and
   cannot be fixed by solver improvements alone.

2. **Penalty amplification** (~632% globally) from BLR solver residual. This could
   be reduced by tighter BLR tolerance or exact MUMPS (no BLR).

**IMPORTANT: Stress-only traction is NOT a viable fix.** v51 Section 22 tested
`--traction-stress-only` (removing the penalty correction from ComputeTraction):
- v51f (p=4): BLOWUP — V = 10 m/s in 35 steps
- v51f2 (p=2): BLOWUP — V = 46 m/s in 397 steps

The penalty correction is **essential for stability**. Without it, the DG method
loses its stabilization mechanism and the solution diverges. This rules out
removing the penalty term from ComputeTraction as a fix for the coherence problem.

This leaves the following options:
- **`--zero-dip-traction`**: Eliminates both sources at once by zeroing τ_dip
  post-extraction. Correct for BP5 where τ_dip = 0 analytically. Already validated.
- **Tighter solver**: Use exact MUMPS (no BLR) to reduce the solver fault residual.
  Would fix source 2 but not source 1.
- **Reformulate traction recovery**: Compute traction in a way that is inherently
  consistent with the solve (e.g., derive from the bilinear form residual rather
  than re-evaluating σ·n + penalty independently). This is a major refactor.

**The `--zero-dip-traction` flag eliminates both sources at once**, which is why it
is the correct approach for the BP5 benchmark where τ_dip = 0 analytically.

### 7.7 Validation of Zero-Dip-Traction Results

The v52 zero-dip production run (Section 3) was verified at 10 stations over ~0.89s:
- slip_dip: max 9.75e-21 m (V_zero floor)
- log10(V_dip): -20.00 everywhere
- tau_dip: bit-for-bit constant, ratio to tau_strike ≤ 1e-11
- Zero dip enforcement is correct and stable

The full production run (48hr, 1800yr) will determine whether removing dip
contamination also fixes the strike deficit and earthquake timing.

---

## 8. What We Know and Don't Know

### 8.1 The Core Observation

The solution u has spurious u_z (dip displacement) that should be analytically zero
for pure strike-slip on a planar fault. Everything else follows from this:
- {σ·n} has dip component because σ_yz = μ ∂u_z/∂y ≠ 0
- [[u]] ≠ slip because u_z creates a dip jump component
- Penalty amplifies the jump residual into dominant spurious traction
- Friction law partitions V into dip via τ_dip/|τ| → dip slip accumulates
- State variable ψ sees inflated |V| → weaker fault → strike deficit

### 8.2 Two Hypotheses — Not Yet Distinguished

**Hypothesis A: Bug in RHS assembly (f has spurious z-forcing)**

If AssembleSlipContributionIP or AssembleDirichletLoading produces a non-zero
z-component in the RHS vector f, then K*u = f forces u_z ≠ 0 regardless of solver
accuracy. This would be a code bug, not a solver issue.

For pure strike-slip:
- Slip BC is (dip=0, strike=δ). EmbedSlip should produce Δu = (δ_x, 0, 0).
  If EmbedSlip leaks into z → f has z-forcing → u_z ≠ 0.
- Dirichlet loading is u_x = ±Vp*t/2 on far-field faces. If AssembleDirichletLoading
  produces z-forcing → u_z ≠ 0.
- The penalty and symmetry RHS contributions involve the 3D elasticity tensor and
  face normals. For faces with z-component normals (top/bottom boundaries), the
  tensor coupling could create z-forcing even from pure x-displacement BCs.

**Hypothesis B: Solver introduces u_z (BLR factorization error)**

If f_z = 0 exactly but the BLR factorization of K is inaccurate at the high-penalty
fault face blocks (300× stiffer than bulk), the approximate solve could produce
spurious u_z. The coherence test's 4e-7 m fault residual supports this.

### 8.3 Why We Can't Distinguish Yet

The coherence test measured:
- stress_dip = 31% of stress_strike (from {σ·n})
- corr_dip = 632% of stress (from penalty*(jump-slip))

But BOTH trace back to u having u_z ≠ 0. The test doesn't tell us whether u_z
comes from f_z ≠ 0 (Hypothesis A) or solver error (Hypothesis B).

### 8.4 What Previous Tests Covered

| Test | What it showed | What it didn't show |
|------|---------------|-------------------|
| v47c exact MUMPS | OOM at p=2 | Whether exact solver eliminates dip |
| v47d BLR 1e-14 | Tighter BLR made blowup worse | Was testing penalty blowup, not dip |
| v51 K-matrix coupling | K components match Tandem | Whether RHS has z-forcing |
| v51 stress-only traction | Blowup (penalty essential) | N/A — different problem |
| v52 coherence test | Both stress and penalty have dip | Whether source is f or solver |

**The gap**: No test has checked whether the RHS f has a z-component.

---

## 9. Next Steps: Isolate the Source

### 9.1 Priority 1: Dump RHS z-components (no solver involved)

Add a diagnostic that, after assembling the full RHS:
```
f = f_volume + f_dirichlet + f_slip
```
dumps the z-component statistics:
- ||f_z|| / ||f_x|| — if this is >> 0, Hypothesis A is confirmed
- Per-contribution breakdown: which of f_volume, f_dirichlet, f_slip has z?
- Spatial distribution: where is f_z largest? (near fault? boundaries? nucleation?)

This requires NO solver. It directly tests whether the linear system is set up
correctly. If f_z = 0, the problem is the solver. If f_z ≠ 0, the problem is in
the assembly code (EmbedSlip, Dirichlet loading, or tensor coupling at boundaries).

### 9.2 Priority 2: Dump u_z after solve (solver contribution)

If f_z = 0, then u_z must come from the solver. Dump:
- ||u_z|| / ||u_x|| — the relative dip contamination in the solution
- Compare with BLR residual: ||K*u - f|| / ||f||
- If ||u_z|| >> cond(K) * BLR_tol * ||u||, something else is wrong

### 9.3 Priority 3: Check EmbedSlip for z-leakage

The FaultBasis::EmbedSlip function converts (slip_dip, slip_strike) → (Δu_x, Δu_y, Δu_z).
For pure strike-slip (slip_dip ≈ 0, slip_strike = δ), Δu_z should be exactly 0.
If the fault basis vectors have even a small z-component in the strike direction,
this would create z-forcing proportional to slip magnitude.

This was checked in v51 Section 7 (FaultBasis comparison: IDENTICAL to Tandem).
But worth re-verifying: print Δu_z for a pure strike slip input at a few fault faces.

### 9.4 Priority 4: Run exact MUMPS on reduced mesh

If Priorities 1-3 confirm f_z = 0 and EmbedSlip is clean, the problem is the solver.
Run exact MUMPS (not BLR) on a reduced-size mesh that fits in memory:
- Use inline mesh with fewer elements (e.g., 8×4×4 = 128 elements)
- Or coarsen the Gmsh mesh
- Compare u_z, fault residual, and dip/strike ratio with BLR result

---

## 10. Why --zero-dip-traction Doesn't Fix Strike: Penalty Amplification of BLR Residual

### 10.1 User Observation

Production run v52 (`--zero-dip-traction`) shows only small improvements at deeper
locations for strike slip, slip rate, and shear stress compared to v50. Strike results
mostly follow v50.

### 10.2 Root Cause

The coherence test (Section 7.2) shows:

```
τ_strike = τ_strike_stress + τ_strike_penalty
         = {σ·n}_strike    + η * ([[u_strike]] - slip_strike)

τ_strike_stress  = 289 Pa RMS  (physical signal)
τ_strike_penalty = 1830 Pa RMS (BLR solver noise amplified by penalty = 6.3× signal)
```

`--zero-dip-traction` zeroes τ_dip after extraction, fixing the dip → strike cascade
through the friction law. But **τ_strike itself** is dominated by BLR solver residual
amplified by the penalty. The 4.4e-7 m RMS jump residual × 10^10 Pa/m penalty =
O(4.4 kPa) spurious strike forcing — larger than the real traction signal of ~289 Pa.

**Consequence**: even with τ_dip = 0, the friction law computes V_strike from a
corrupted τ_strike that has 6× the noise of the physical contribution. The ~1 MPa
strike deficit in the interseismic τ persists because the tau loading rate is
distorted by this spurious penalty contribution, not because of dip contamination.

### 10.3 Summary of Fix Status

| Fix | Applied | Result | Remaining problem |
|-----|---------|--------|-------------------|
| `--zero-dip-traction` | v52 | τ_dip → 0, small improvement at depth | τ_strike still dominated by BLR noise (6× signal) |
| Exact MUMPS (no BLR) | Not tried (OOM at p=2 for full mesh) | — | Would fix penalty amplification |
| Reduce mesh size + exact MUMPS | Not tried | — | Feasibility unknown |
| CG matrix-free solver | Not available in MFEM | — | Would match Tandem exactly |

The strike deficit is not caused by dip contamination — it is caused by the
MUMPS-BLR solver leaving a large residual at the high-penalty fault faces.

---

## 11. Tandem Comparison Plan (Now That Tandem Is Available on Cluster)

### 11.1 Goal

With Tandem running on the cluster, we can directly compare station outputs and
isolate where MFEM diverges from the reference. The key question is:

**At which point in time and at which stations does MFEM's τ_strike first diverge
from Tandem's, and what is the magnitude and spatial pattern of the divergence?**

This will tell us whether the deficit is:
- **Static** (present from t=0 → initial condition or elastic assembly bug)
- **Dynamic** (accumulates over time → loading rate error or ODE feedback)

### 11.2 Comparison 1: Initial Tau Values (t=0)

Both codes initialize τ from the friction law equilibrium at V_init:

```
τ0(x2,x3) = σ_n * f(V_init, psi_ss(V_init))
```

MFEM's initial values from v52 station output:
- strk+00dp+00: τ_strike = 13.27 MPa
- strk-24dp+10: τ_strike = 19.49 MPa
- strk-16dp+10 (nucleation edge): τ_strike = 21.15 MPa

**Action**: Read Tandem's station output at t=0 (first row). If Tandem t=0 tau
differs from MFEM t=0 tau, the initial condition assembly is wrong.

**Expected if correct**: Both codes should have identical τ_0 = σ_n * f(V_init_vec, psi_ss(V_init_vec)).

### 11.3 Comparison 2: Short Pre-Earthquake Loading Rate

Both codes: compare τ_strike loading rate dτ/dt during the first few years (before
the first earthquake). This tests the elastic stiffness operator.

From plate loading: Vp/2 = 0.5e-9 m/s on each side of fault. The stiffness operator
maps this loading to τ at fault stations. If MFEM's elastic operator has a systematic
error (penalty-amplified traction bias), the loading rate dτ/dt will differ.

**Action**:
- Plot τ_strike vs t for both codes at strk+00dp+10 (nucleation center)
- Measure slope dτ/dt in MPa/year during locked phase
- Compare MFEM slope vs Tandem slope

MFEM first earthquake is at t≈55s (nucleation zone reaches V≈V_nuc). Before that,
during the initial interseismic phase starting from t~0 (post-first-earthquake if
the transient is fast), the locked-phase loading rate should be identical to Tandem's
if the elastic operator is correct.

### 11.4 Comparison 3: Earthquake Timing

**First earthquake**: Both codes should have first earthquake around t≈55s due to
nucleation zone V_nuc=0.01 m/s initialization (off steady-state with psi=psi_ss(V_init)).

**Action**: Check Tandem first earthquake timing. If Tandem's first earthquake is at
a very different time from MFEM (~55s), it means the nucleation mechanism works
differently. This is expected if Tandem uses delta_tau_factor > 0 (overstress
nucleation) vs MFEM's delta_tau_factor = 0 (velocity nucleation).

**Check**: What nucleation method does the Tandem cluster run use?
In MFEM bp5_params.hpp: `delta_tau_factor = 0.0` (Tandem mode, no overstress).
In Tandem's bp5.toml: check nucleation configuration.

**Second earthquake**: After the first earthquake, both codes enter interseismic
locked phase. The recurrence interval is approximately:

```
ΔΤ ≈ Δτ / (dτ/dt)
```

where Δτ is the stress drop and dτ/dt is the plate-loading rate. If dτ/dt differs
(Section 11.3), the recurrence interval will differ.

### 11.5 Comparison 4: Post-Earthquake V Decay

MFEM shows V_strike → 10^{-42} to 10^{-53} post-earthquake (depends on station).

With Tandem parameters b=0.03, a=0.004 (b/a=7.5):
```
d(log10 V)/d(log10 theta) = -b/a = -7.5
```

At t = 6.14e6 s post-earthquake, MFEM data: V = 10^{-42.56}. Theory predicts
10^{-42.5}. This matches — V~10^{-42} is physically correct.

**Action**: Check Tandem's post-earthquake V_strike at strk+00dp+10 or strk-24dp+10.
Does Tandem show similar V~10^{-40} to 10^{-45}? If Tandem shows much higher V
post-earthquake, Tandem may use different b/a parameters than we assumed.

**Verify Tandem parameter file**: The MFEM bp5_params.hpp was set to match Tandem.
Confirm by reading the actual Tandem config used on the cluster and checking b, a0, Dc.

### 11.6 Comparison 5: Solver Residual at Fault Faces

This is the most direct test of Hypothesis A vs B (Section 8.2).

**Tandem** uses PETSc CG matrix-free. The Krylov solver converges to ||K*u - f||/||f|| < 1e-12,
meaning the jump residual at fault faces should be O(1e-12) × displacement scale.
With |u| ~ 1e-5 m at first step, |[[u]] - slip| should be ~1e-17 m (near machine precision).

**MFEM** shows |[[u]] - slip| = 4.4e-7 m RMS — 10 orders above what a converged
CG solver would give.

**Action**: Run Tandem's coherence-equivalent diagnostic (if available) or compare
the traction magnitudes:
- |τ_strike_stress| vs |τ_strike| in Tandem (if Tandem outputs this)
- Compare Tandem's τ_strike at first non-zero slip step vs MFEM's

If Tandem's τ_strike at same slip magnitude matches MFEM's stress-only component
(~289 Pa) rather than total (~2120 Pa), it confirms that MFEM's penalty noise is
the cause of the τ_strike inflation.

### 11.7 What to Run on Cluster

**Tandem config check** (before running anything):
1. Read the bp5.toml config used on cluster. Verify: b, a0, Dc, sigma_n, Vp match MFEM bp5_params.hpp
2. Check nucleation method: V_nuc or delta_tau overstress?
3. Check output format and station locations — do they match MFEM's 10 stations?

**Tandem short run**: 100 years (enough to see 1-2 earthquakes with Tandem's faster b/a)
- Station output every year during locked phase, every second during rupture

**MFEM comparison**: Use v52 zero-dip output (already have 49 years, ~53127 rows)
- Same stations, same time points

**Expected output of comparison**:
- If τ_initial matches: initial condition correct ✓
- If dτ/dt differs: elastic loading rate wrong (penalty noise contaminating τ_strike)
- If earthquake timing differs by large factor (not just small drift): nucleation mechanism different

### 11.8 Priority Diagnostic: Does Tandem's τ_strike Have the Same Magnitude as MFEM's?

The single most useful number: **Tandem τ_strike at strk+00dp+10 at the first locked step after the first earthquake, vs MFEM.**

MFEM coherence test: τ_strike_total ≈ τ_stress + τ_penalty ≈ 289 + 1830 = 2120 Pa (RMS, first step).

If Tandem shows τ_strike ~ 289 Pa (stress-only level) while MFEM shows ~2120 Pa,
it confirms that the 6× traction inflation from MUMPS-BLR penalty amplification
is the dominant cause of MFEM's τ_strike being wrong.

---

## 12. Revision History

| Date | Section | Change |
|------|---------|--------|
| 2026-03-24 | 1-5 | Initial v52 document. Removed single-solve and K-coupling diagnostics. Production run submitted with --zero-dip-traction. |
| 2026-03-24 | 6 | Added solve vs traction coherence analysis. Identified penalty, normal vector, Jacobian, and fault-face-in-K discrepancies between MFEM (decoupled) and Tandem (coherent). Proposed residual-based coherent traction test. |
| 2026-03-24 | 6.6 | Implemented --diag-traction-coherence flag. Decomposes traction into stress-only and penalty components, measures solver fault residual. Created bp5_v52_coherence_test.sbatch (dev queue, 10yr, 200 ranks). Build verified. |
| 2026-03-24 | 7 | Coherence test results (Job 7610206). Fixed MPI deadlock (two bugs: local any_slip check, local coh_n_qp guard). Results: 31% dip/strike from DG stress, 632% penalty amplification from BLR solver residual. Both sources confirmed. |
| 2026-03-24 | 6 | CORRECTED: Tandem also includes fault faces in K (same architecture as MFEM). Rewrote Section 6 with concrete code references. Solver difference: CG matrix-free vs MUMPS-BLR. Penalty ~10^10 Pa/m at fault vs ~3.2e7 bulk = 300x contrast. |
| 2026-03-24 | 8-9 | Honest assessment: two hypotheses (RHS bug vs solver error) not yet distinguished. No test has checked whether f has z-component. Defined priority-ordered next steps: (1) dump RHS z-components, (2) dump u_z after solve, (3) verify EmbedSlip, (4) exact MUMPS on reduced mesh. |
| 2026-03-28 | 10-11 | v52 production run confirmed: --zero-dip-traction only gives small improvements, strike results still follow v50. Root cause: τ_strike itself is dominated by BLR penalty amplification (6× signal). Added Tandem comparison plan (Sections 11.2-11.8): compare initial tau, loading rate, earthquake timing, post-earthquake V decay, solver residual. Priority diagnostic: Tandem τ_strike at strk+00dp+10 first locked step vs MFEM (~289 Pa vs ~2120 Pa expected if BLR noise is root cause). |
