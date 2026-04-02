# BP5 Debug v57: MPI Core-Count Dependence of Solution

**Date:** 2026-04-01
**Issue:** Simulation results depend on number of MPI ranks (100 vs 400)
**Severity:** Systematic divergence growing over time (~1% at step 1, ~7% by step 410)

## 1. Problem Statement

Two runs of the same BP5 v57 configuration (1000m, IP, p=1, PETSc TS RK45, MUMPS direct solver) produce different results depending on MPI rank count:

| Job | Ranks | Output |
|-----|-------|--------|
| 7627262 | 100 | `bp5_v57_pts_p1_mumps` |
| 7627480 | 400 | `bp5_v57_pts_p1_mumps_nowrite` |

Both use the same mesh (`bp5_tandem_exact.msh`), same parameters, same solver, and same PETSc TS configuration (`petsc_ts_rk45_tandem.cfg`).

### Observed Divergence

| Time (yr) | V_max (400 ranks) | V_max (100 ranks) | Difference |
|-----------|--------------------|--------------------|------------|
| 3.17e-10  | 1.097e-02          | 1.108e-02          | ~1.0%      |
| 7.3e-08   | 1.374e-02          | 1.399e-02          | ~1.8%      |
| 3.02e-07  | 1.106e-02          | 1.188e-02          | **~7.4%**  |

The divergence grows systematically over time.

## 2. Proven Facts (Frontera Diagnostic Data)

Targeted diagnostics in the first non-trivial RK45 stage:

| Quantity | 100 ranks | 400 ranks | Interpretation |
|----------|-----------|-----------|----------------|
| `||slip||_inf` | 2.00e-05 | 2.00e-05 | identical input |
| `sum|K_ij|` | 8.119568279813302e+21 | 8.119568279813346e+21 | K entry values match (5e-15 rel) |
| `K_nnz` | 43,434,416 | 43,443,554 | **NNZ differs by 9,138** |
| `||b||_2` | 2.324035e+11 | 2.323174e+11 | **RHS differs 0.037%** |
| `||u||_inf` | 6.845e-05 | 6.410e-05 | **displacement differs 6.4%** |
| P = identity | yes (both) | yes (both) | DG prolongation correct |
| global rows/cols | 754,488 | 754,488 | same global DOF count |
| interior faces (sum) | 114,808 | 106,618 | fewer with more ranks |
| shared faces (sum) | 19,626 | 36,006 | more with more ranks |
| interior + shared/2 | **124,621** | **124,621** | face accounting correct |
| shared fault faces | 4 | 55 | very few shared fault faces |
| boundary elements | 12,954 | 12,954 | same |

**What is established:**
- The slip (fault state input) is identical.
- The stiffness matrix entry VALUES are the same to machine precision.
- The stiffness matrix NNZ differs by 9,138 (structural zeros, not value differences).
- The RHS vector differs by 0.037%.
- The displacement differs by 6.4% (amplified from RHS).
- Face counting is correct: no faces lost or double-counted.
- The prolongation P is identity (correct for DG).

**What is NOT yet established:**
- The exact mechanism producing the 0.037% RHS difference.
- Whether the NNZ difference contributes to the displacement error or is benign.
- Whether the reparameterization hypothesis is correct (see Section 4).

## 3. Discarded Hypotheses

### 3.1 Shared Fault DOF Duplication (Discarded)

**Hypothesis:** Both ranks independently evolve fault state (psi, slip, V) for
shared fault faces, causing divergent evolution.

**Fix attempted:** Owned-fault layout with ghost synchronization (commits f94d000, 978155d).

**Result:** No effect. Only 4 shared fault face instances exist across 100 ranks.
`global_owned_dofs = global_local_dofs = 27936` — the state vector is identical.

### 3.2 DOF-Within-Face Permutation (Discarded)

**Hypothesis:** DOF ordering within shared faces differs between ranks, causing
communicated values to land at wrong physical locations.

**Fix attempted:** Canonical DOF permutation sorted by global vertex ID (commit 978155d).

**Result:** No effect. The permutation changes the owned↔local mapping, but with
essentially no shared fault faces, no communication occurs.

### 3.3 KeepNbrBlock and skip_zeros (Discarded)

**Fix attempted:** `KeepNbrBlock(true)` (commit 6cf06f8) and `Assemble(skip_zeros=0)` (commit cad9a31).

**Result:** No effect on HypreParMatrix NNZ or values.

## 4. Working Hypothesis: Shared-Face Parameterization Inconsistency

### 4.1 Hypothesis

MFEM's DG assembly splits face processing between `BilinearForm::Assemble()`
(interior faces) and `ParBilinearForm::AssembleSharedFaces()` (shared faces).
For the same physical face processed as interior on one partition vs shared on
another, the `FaceElementTransformations` may differ because:

1. The face vertex ordering in `faces[FaceNo]->GetVertices()` is element-dependent
2. `CalcOrtho(Trans.Jacobian(), nor)` gives different normals (`eltrans.cpp`)
3. `Trans.GetElement1IntPoint()` / `GetElement2IntPoint()` map face reference
   points to different physical locations via Loc1/Loc2 (`eltrans.cpp:606`)

The integrator (`dg_elasticity_ip_combined_integrator.hpp:45`) uses both
`CalcOrtho` and `GetElement1/2IntPoint`, so any inconsistency propagates
into the face matrix and RHS.

For polynomial integrands on planar faces with exact quadrature, the integral
should be parameterization-independent (the DG formulation is orientation-invariant).
But the slip RHS interpolates fault DOF values that can have sub-element
discontinuities (nucleation zone `eps=0.001m << h=1000m`), which may not be
exactly representable and could sample differently under different parameterizations.

### 4.2 Supporting Evidence

- The `sum|K_ij|` matches to 5e-15, consistent with the DG bilinear form being
  orientation-invariant for polynomial integrands.
- The NNZ difference (9,138) suggests some structural entries are zero under one
  parameterization but O(ε) under another.
- The 0.037% RHS difference is order-of-magnitude consistent with a small number
  of shared fault faces (~4 out of ~9,312) having different slip interpolation.
- The displacement error (6.4%) is consistent with κ ≈ 160 amplification.

### 4.3 What Remains Unproven

The hypothesis has NOT been verified by the decisive micro-test:

> Assemble the same triangular face as serial interior and as parallel shared.
> Compare the local face matrix and slip RHS after allowing only the expected
> block swap / side permutation. Show the mismatch, then apply a canonical
> Face+Loc1+Loc2 reparameterization and show the mismatch disappears.

The small hex/tet test mesh does not reproduce the bug (no shared fault faces
with METIS partitioning on a 48-element mesh). The Frontera data supports the
hypothesis but does not isolate the mechanism.

### 4.4 Proposed Fix (Pending Micro-Test Verification)

Following Tandem's sorted-simplex approach, canonicalize the FULL face reference
mapping for shared faces before any face integrator call:

1. **Face transformation** — permute PointMat columns to canonical vertex order
2. **Loc1** — compose reference-face permutation into face→elem1 mapping
3. **Loc2** — compose reference-face permutation into face→elem2 mapping

All three must be consistent. A PointMat-only patch would leave `CalcOrtho` and
`GetElement1/2IntPoint` on different parameterizations.

### 4.5 Key MFEM Code Paths

- Face transformation: `mesh.cpp:540-561` — PointMat from face vertices
- Face-to-element mapping: `eltrans.cpp:606` — eip1/eip2 through Loc1/Loc2
- Interior face assembly: `bilinearform.cpp:672-696`
- Shared face assembly: `pbilinearform.cpp:229-272`
- Slip RHS: `elasticity_operator.hpp:AssembleSlipContributionIPShared`

## 5. Next Step: Micro-Test

Build a test that:
1. Constructs a mesh with a known triangular fault face
2. Processes that face as serial interior → records face matrix F_serial and slip RHS b_serial
3. Processes the same face as parallel shared (2 ranks) → records face matrix F_shared and slip RHS b_shared
4. Compares F_serial vs F_shared (after expected block swap) and b_serial vs b_shared
5. If they differ: applies Face+Loc1+Loc2 canonicalization and shows the mismatch disappears

This proves or disproves the hypothesis before implementing the full wrapper.

## 6. Verification Jobs

Three SBATCH scripts for Frontera validation (post-fix):

| Script | Ranks | Nodes | Purpose |
|--------|-------|-------|---------|
| `bp5_v57_mpi_test_100ranks.sbatch` | 100 | 2 | Baseline |
| `bp5_v57_mpi_test_200ranks.sbatch` | 200 | 4 | Intermediate |
| `bp5_v57_mpi_test_400ranks.sbatch` | 400 | 8 | Original scale |

## 7. References

- Tandem sorted simplex: `/Users/chunhuizhao/projects/tandem/src/mesh/Simplex.h:15-37`
- Tandem facetParam: `/Users/chunhuizhao/projects/tandem/src/geometry/Curvilinear.cpp:333-346`
- MFEM face transformation: `/Users/chunhuizhao/projects/mfem/mesh/mesh.cpp:540-561`
- MFEM face-to-element mapping: `/Users/chunhuizhao/projects/mfem/fem/eltrans.cpp:606`
- MFEM shared face assembly: `/Users/chunhuizhao/projects/mfem/fem/pbilinearform.cpp:229-272`
- Original runs: job 7627262 (100 ranks), job 7627480 (400 ranks)
- Diagnostic runs: jobs 7627560, 7627563, 7627565, 7627580, 7627582
