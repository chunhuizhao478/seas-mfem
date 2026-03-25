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

## 6. Solve vs. Traction Coherence: MFEM Decoupled, Tandem Coherent

### 6.1 The Coherence Problem

Tandem and MFEM handle elastic solve and traction extraction differently. This is a
potential source of the dip contamination and strike traction error.

**Tandem: Coherent pipeline.** The elastic solve and traction extraction use the
**exact same DG flux formula**:

```
Tandem solve:     surfaceOp kernel (elasticity.py:105-110)
  traction_q = 0.5 * (σ(u⁺)·n + σ(u⁻)·n)     [consistency]
  + penalty * (jump - slip_bc)                   [stabilization]

Tandem traction:  compute_traction kernel (elasticity.py:242-248)
  traction_q = 0.5 * (σ(u⁺)·n + σ(u⁻)·n)     [SAME consistency]
  + c0 * (jump - slip_bc)                        [SAME penalty]
```

Both kernels call the identical `traction(x, normal)` function (lines 89-91), which
computes the Cauchy stress σ·n = λ(∇·u)n + μ(∇u + ∇uᵀ)·n. The penalty coefficients
are the same. The displacement basis evaluations (`E_q`) are the same. The slip BC
(`f_q`) is the same. **The traction Tandem extracts is exactly the flux the solver used.**

The code flow (SeasQDOperator.rhs(), lines 34-40):
```cpp
solve(time, state);              // K*u = f  (surfaceOp kernel)
update_traction(state);          // T = flux(u) at fault  (compute_traction kernel)
friction_->rhs(time, traction_, state, result);
```

The `update_traction` call (SeasQDOperator.cpp:69-74) goes through
`AdapterOperator::traction()` which calls `traction_skeleton()` /
`traction_boundary()` on each fault face, using the solved displacement. These methods
(Elasticity.cpp:948-1017) invoke the `compute_traction` kernel — the same formula as
the bilinear form's face integral.

**MFEM: Decoupled pipeline.** The elastic solve and traction extraction use
**different code paths with independently computed operators**:

```
MFEM solve:
  Stiffness K:   MFEM's DGElasticityIntegrator (kappa=0) [consistency + symmetry]
                 + DGElasticityIPPenaltyIntegrator        [penalty]
                 → assembled into HypreParMatrix
  RHS f:         AssembleSlipContributionIP()              [hand-coded slip RHS]
                 + AssembleDirichletLoading()               [hand-coded BC RHS]
  Solve:         K * u = f  via MUMPS/CG

MFEM traction:  ComputeTraction() — entirely separate hand-coded computation
  T_stress_q:    ∇u evaluated manually at quad points, averaged, multiplied by C:n
  correction_q:  penalty * (u_jump - slip), penalty recomputed from scratch
  T = T_stress - correction, projected onto fault basis
```

### 6.2 Specific Discrepancies Identified

#### 6.2.1 Penalty computation path

**Solve** (AssembleSlipContributionIP, line 1200-1208):
- `nor` = `CalcOrtho(FTr->Jacobian())` at each quad point (varies with q)
- `nl_q = nor.Norml2()` — the face Jacobian determinant
- penalty uses `nl_q` at each quad point
- Applied weight: `wq_penalty = penalty_ip * ip.weight * nl_q`

**ComputeTraction** (line 3416-3427):
- `nor` = `CalcOrtho(FTr->Jacobian())` at **face centroid only** (line 3392)
- `face_area = nor.Norml2()` — face Jacobian at centroid
- penalty uses `face_area` (constant, not per-quad-point)
- Applied per quad point but with centroid-derived penalty

For non-planar faces or non-uniform Jacobians, these differ. On the hex meshes used
for BP5, faces are planar, so this should be exact. But it's still a latent bug for
general meshes.

#### 6.2.2 Stress computation: adjugate vs inverse Jacobian

**Solve** (line 1185-1196):
- Uses `CalcAdjugate` → `dshape_adj = dshape_ref * adjJ`
- Weight: `w = ip.weight / (2.0 * detJ)` which divides out the adjugate's detJ factor
- Net: `dshape_phys = dshape_ref * J⁻¹` (mathematically equivalent)

**ComputeTraction** (line 3396-3506):
- Uses `CalcInverse` → `dshape_phys = dshape_ref * J⁻¹` directly
- No detJ factor in weights (just `ip.weight`, accumulated into `sum_wq`)
- Projects via `GalerkinProject` which divides by face mass matrix

These are mathematically equivalent for well-conditioned elements. But the adjugate
path avoids division by detJ, while the inverse path requires it. If detJ is small
(thin elements), the inverse path may have different rounding characteristics.

#### 6.2.3 Normal vector: centroid vs per-quad-point

**Solve**: `nor` is recomputed at every quad point via `CalcOrtho(FTr->Jacobian())`
at the current integration point. For the stress-related terms (symmetry/consistency),
the normal varies with the quad point.

**ComputeTraction**: `nor` is computed once at the face centroid (line 3392) and
`basis.normal` (from FaultBasis, constant per face) is used for stress dot product
(line 3534): `T_stress_q[ci] += stress_ij * basis.normal[cj]`.

For planar faces, `CalcOrtho` returns the same direction at every quad point (magnitude
may differ due to Jacobian parameterization, but direction is constant). `basis.normal`
is a unit vector. So the traction direction is consistent. But the **magnitude scaling**
differs:

- Solve uses `nor` (not unit — includes face area factor)
- ComputeTraction uses `basis.normal` (unit vector) for the stress·n dot product

This is compensated by the different quadrature weight handling, but it means the two
code paths must be carefully balanced. Any imbalance creates a traction that doesn't
match what the solver "expects."

#### 6.2.4 Fault face exclusion from stiffness matrix

**Critical**: The MFEM stiffness matrix `K` is assembled by MFEM's
`DGElasticityIntegrator` which loops over **all interior faces**, including fault faces.
This means `K` treats fault faces as regular DG interfaces with consistency + symmetry +
penalty terms on the displacement jump.

The slip contribution `AssembleSlipContributionIP` then adds the **RHS correction** for
the prescribed jump `[[u]] = slip` on fault faces.

In the solve, the stiffness matrix enforces `K*u = f` where `K` includes the penalty
for `[[u]] = 0` on fault faces (as if no slip), and `f` includes the correction for the
prescribed slip. The solution `u` satisfies:

```
K_regular * u + penalty * (u_jump) = f_vol + penalty * slip + symmetry * slip
```

After solving, `u_jump ≈ slip` (enforced by the penalty). **ComputeTraction** then
re-evaluates the flux formula to extract the traction.

In Tandem, the stiffness matrix is **matrix-free** — the same kernel (`surfaceOp`)
handles both the matrix-vector product and the slip BC simultaneously. There's no
separation between "matrix treats fault as no-slip" and "RHS corrects for slip."
The `set_slip()` call modifies what the kernel sees as the prescribed jump, so the
matrix-vector product itself incorporates the slip BC coherently.

### 6.3 Why This Matters for the Dip Error

The decoupled approach creates a subtle issue: **the traction ComputeTraction extracts
may not exactly equal the traction the solver implicitly balanced against.**

If the solve balances forces using one penalty and normal computation, but
ComputeTraction evaluates with slightly different quadrature/normals, the extracted
traction has an error proportional to the difference. For pure strike-slip:

- The strike traction error is O(ε) × τ_strike — a small relative error
- The dip traction should be exactly zero, but the error is O(ε) × τ_strike
- This creates a **spurious dip traction** of order ε × 20 MPa ≈ O(0.01-1 MPa)

This is consistent with the 0.1-1.2 MPa dip offsets observed in v51 Section 3.

### 6.4 Proposed Test: Coherent Traction Extraction

To test whether the solve/traction decoupling is the root cause, we can implement a
**coherent traction test** that extracts traction using the same operator used in the
solve, rather than the separate ComputeTraction code.

#### Test Design: Residual-Based Traction

Instead of re-evaluating σ·n and penalty separately in ComputeTraction, compute
traction from the **residual of the bilinear form** restricted to fault faces:

```
T_fault = (K * u - f_non_fault) restricted to fault DOFs
```

This is exactly what the solver balanced against. If `K * u = f_total`, then:
```
f_total = f_volume + f_dirichlet + f_slip_fault
```
and:
```
K * u = f_non_fault + f_slip_fault
```
The traction on fault faces is embedded in the K*u product on fault-adjacent DOFs.

Concretely, for each fault face DOF `i`:
```
traction_i = (K * u)(i) - f_non_fault(i)
```

This is automatically coherent with the solve because it uses the **same K matrix**.

#### Implementation Outline

```cpp
void ComputeTractionCoherent(const GridFuncType &displacement,
                              const Vector &slip_bc, real_t time,
                              Vector &traction)
{
   // 1. Compute K * u (matrix-vector product)
   Vector Ku(displacement.Size());
   cached_Ah_.Mult(displacement, Ku);

   // 2. Compute f_non_fault = f_volume + f_dirichlet (no slip contribution)
   LinFormType b_no_slip(fes_.get());
   b_no_slip.Assemble();
   AssembleDirichletLoading(b_no_slip, time);
   // Do NOT add AssembleSlipContributionIP

   // 3. Residual on fault DOFs = K*u - f_non_fault
   //    This residual encodes the traction the solver balanced
   Vector residual(displacement.Size());
   residual = Ku;
   residual -= b_no_slip;

   // 4. Extract fault-face residual and project to (dip, strike) via FaultBasis
   //    The residual at fault DOFs is the integrated traction × test function
   //    Need to "un-integrate" via inverse face mass matrix (L2 projection)
   for each fault face fi:
      extract residual at face DOFs
      apply inverse face mass matrix
      project to local (dip, strike) frame
}
```

#### What This Tests

If the coherent traction extraction produces smaller dip contamination than the current
ComputeTraction, it confirms the decoupled approach is a significant error source.

The test should be run on the same mesh/parameters as the production run (1000m, p=2, IP)
and compare:
1. Current `ComputeTraction` dip/strike ratio per DOF
2. Coherent residual-based dip/strike ratio per DOF

Expected: coherent extraction has O(solver_tol) dip contamination, while current has
O(21%) from the re-evaluation inconsistency.

### 6.5 Alternative: Match ComputeTraction to Solve Exactly

Instead of the residual approach, we could fix ComputeTraction to use exactly the same
computation as the solve:

1. Use `CalcOrtho` per quad point (not centroid-only) for penalty
2. Use adjugate Jacobian (not inverse) to match rounding
3. Use the same quadrature rule as the bilinear form integrator
4. Include the symmetry term in the traction (currently ComputeTraction may omit it)

This is harder to verify than the residual approach but avoids the need for inverse mass
matrices.

### 6.6 Implemented: `--diag-traction-coherence` Flag

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

## 7. Revision History

| Date | Section | Change |
|------|---------|--------|
| 2026-03-24 | 1-5 | Initial v52 document. Removed single-solve and K-coupling diagnostics. Production run submitted with --zero-dip-traction. |
| 2026-03-24 | 6 | Added solve vs traction coherence analysis. Identified penalty, normal vector, Jacobian, and fault-face-in-K discrepancies between MFEM (decoupled) and Tandem (coherent). Proposed residual-based coherent traction test. |
| 2026-03-24 | 6.6 | Implemented --diag-traction-coherence flag. Decomposes traction into stress-only and penalty components, measures solver fault residual. Created bp5_v52_coherence_test.sbatch (dev queue, 10yr, 200 ranks). Build verified. |
