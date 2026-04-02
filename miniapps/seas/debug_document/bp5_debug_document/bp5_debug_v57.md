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

The divergence grows systematically over time, ruling out MUMPS floating-point sensitivity (which would give ~1e-12 level differences).

## 2. Root Cause Analysis

### 2.1 Tandem Reference: Canonical Fault DOF Ordering

Tandem (`src/form/BoundaryMap.cpp:43-51`) sorts fault faces by a **global face ID** (`facets.l2cg(fctNo)`), ensuring deterministic, partition-independent DOF ordering:

```cpp
std::sort(theFctNos.begin(), theFctNos.end(),
    [&rank, &facets](auto const& a, auto const& b) {
        return std::make_tuple(!(rank <= a.first), a.first, facets.l2cg(a.second)) <
               std::make_tuple(!(rank <= b.first), b.first, facets.l2cg(b.second));
    });
```

Additionally, Tandem uses `ScatterPlan` to synchronize shared fault data between ranks.

### 2.2 MFEM-SEAS: No Ownership or Synchronization

In MFEM-SEAS (`elasticity_operator.hpp:580-605`), fault faces are appended in local mesh iteration order with no global canonical ordering:

```cpp
for (int f = 0; f < num_faces; f++) {
    if (IsFaultFace3D(FTr))
        fault_interior_faces_.Append(f);   // :588
}
for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++) {
    if (IsFaultFace3DShared(sf))
        fault_shared_faces_.Append(sf);    // :605
}
```

**Critical consequence:** For shared fault faces, both ranks independently maintain and evolve their own copy of the fault DOFs (psi, slip, V). There is no mechanism to keep these copies synchronized.

### 2.3 Why This Causes Divergence

1. Different MPI partitioning → different faces become "interior" vs "shared"
2. For shared fault faces, both ranks run the rate-state friction law independently on their own copies
3. The face `FaceElementTransformations` differ between ranks (different elem1), causing the nodal rule to map DOF indices to different physical coordinates
4. Different physical coordinates → different friction parameters (a, V_init, tau_pre) at each DOF index
5. Different parameters → different state evolution → different slip values contributed to the elasticity RHS from each side of the shared face
6. The inconsistency feeds back through the elasticity solve and traction recovery, amplifying over time

## 3. Implemented Fix (v57 Owned-Fault Layout)

### 3.1 Architecture: Owned vs Ghost Fault DOFs

The fix introduces an **owned-fault / ghost-fault** pattern, analogous to Tandem's `ScatterPlan`:

- **Owned fault DOFs:** Interior faces are always owned. For shared faces, the rank with the lowest rank number among sharing ranks is the owner.
- **Ghost fault DOFs:** Non-owner ranks store ghost copies that are populated from the owner before each domain solve.
- **State vector:** Contains only owned DOFs (no duplication across ranks).

### 3.2 Key Code Changes

**`elasticity_operator.hpp` — `BuildOwnedFaultLayout()`:**
- Uses `mesh_.GetGlobalFaceIndices()` to determine global face IDs
- `MPI_Allgather` of shared face GIDs to determine ownership (lowest rank wins)
- Builds `SharedFaultCommBlock` structures for MPI send/recv

**`elasticity_operator.hpp` — `RestrictToOwnedFault()` / `ExpandOwnedToLocalFault()`:**
- `RestrictToOwnedFault`: Extracts owned DOF subset from full local view
- `ExpandOwnedToLocalFault`: Scatters owned DOFs to full local view + MPI communication of owned values to ghost ranks via `MPI_Isend`/`MPI_Irecv`

**`fault_geometry.hpp`:**
- Uses `GetNumOwnedFaultDOFs()` for all parameter/state array sizing
- Restricts coordinates and depths to owned view before computing friction parameters

**`seas_operator.hpp` — `Mult()` and `SetInitialCondition()`:**
- Before domain solve: `ExpandOwnedToLocalFault(slip → local_slip)` (ghost sync)
- After traction recovery: `RestrictToOwnedFault(local_traction → traction)` (owned only)
- Rate-state operates on owned DOFs exclusively

**`test_parallel_elasticity.cpp`:**
- New `test_owned_fault_layout` verifying global owned count equals serial unique count and ghost values match after sync

### 3.3 Data Flow

```
State vector (owned DOFs only)
    │
    ├─ GetSlip() → owned slip
    │       │
    │       ▼
    │   ExpandOwnedToLocalFault() ← MPI ghost sync
    │       │
    │       ▼
    │   local_slip (full local view)
    │       │
    │       ├─ Solve(local_slip) → displacement
    │       │
    │       ├─ ComputeTraction(displacement, local_slip) → local_traction
    │       │
    │       ▼
    │   RestrictToOwnedFault() → owned traction
    │       │
    │       ▼
    └─ ComputeRHS(traction, state) → rate (owned DOFs)
```

## 4. DOF Permutation Fix (v57b)

### 4.1 Problem: First fix attempt had no effect

The initial owned-fault fix (Section 3) produced byte-for-byte identical results to the pre-fix runs because the communicated DOF values landed at wrong physical locations on the receiving rank. MFEM's shared face vertex ordering is element-dependent (`faces[FaceNo]->GetVertices()` differs between ranks), so DOF 0 on the owner maps to a different physical vertex than DOF 0 on the ghost.

### 4.2 Solution: Canonical DOF permutation (Tandem sorted-simplex)

Following Tandem's `Simplex.h` convention, we sort face vertices by ascending global vertex ID to produce a canonical ordering. A per-face permutation maps between MFEM's local DOF order and canonical order.

**Implementation in `BuildOwnedFaultLayout()`:**
1. Get global vertex IDs via `mesh_.GetGlobalVertexIndices()`
2. For each fault face, get face vertices via `mesh_.GetFaceVertices()`
3. Sort (global_vertex_id, mfem_local_index) pairs by global ID
4. Store `canonical_to_local_perm_[face][k]` = MFEM local DOF for canonical DOF k
5. For p=1: DOF k = vertex k, so vertex permutation = DOF permutation

**Applied at three points:**
- `owned_fault_dof_to_local_dof_` includes the permutation (used by `RestrictToOwnedFault` and owned-face expansion)
- `ExpandOwnedToLocalFault` ghost recv applies receiver's `canonical_to_local_perm_`
- `ExpandOwnedToLocalFault` ghost send packs data in canonical order (from owned state which is canonical)

### 4.3 Result: No effect

Both the owned-fault fix and the canonical DOF permutation produced byte-for-byte
identical results to the pre-fix runs. The root cause is NOT shared fault faces
(only 4 instances across 100 ranks — negligible).

## 5. True Root Cause: Shared-Face Parameterization in DG Assembly

### 5.1 Diagnostic Data (Frontera 100 vs 400 ranks)

| Quantity | 100 ranks | 400 ranks | Interpretation |
|----------|-----------|-----------|----------------|
| `||slip||_inf` | 2.00e-05 | 2.00e-05 | identical input |
| `sum|K_ij|` | 8.1196e+21 | 8.1196e+21 | K values same (5e-15 rel) |
| `K_nnz` | 43,434,416 | 43,443,554 | structural zeros differ by 9,138 |
| `||b||_2` | 2.3240e+11 | 2.3232e+11 | **RHS differs 0.037%** |
| `||u||_inf` | 6.85e-05 | 6.41e-05 | **displacement differs 6.4%** |

### 5.2 Mechanism

MFEM's DG assembly splits face processing:
- `BilinearForm::Assemble()` → interior faces (both elements local)
- `ParBilinearForm::AssembleSharedFaces()` → shared faces (elem1 local, elem2 remote)

For the SAME physical face, the face transformation differs because:
1. The face vertex ordering in `faces[FaceNo]->GetVertices()` is element-dependent
2. `CalcOrtho(Trans.Jacobian(), nor)` gives different normals
3. `Trans.GetElement1IntPoint()` / `GetElement2IntPoint()` map face reference
   points to different physical locations via Loc1/Loc2

For polynomial integrands on planar faces with exact quadrature, the integral
is parameterization-independent. But:
- The **slip RHS** interpolates fault DOF values that can have sub-element
  discontinuities (nucleation zone boundary with `nucleation_eps=0.001m << h=1000m`).
  Different quadrature point locations sample this discontinuity differently.
- The **K matrix** has 9,138 structural zeros that differ (entries exactly zero for
  one parameterization but O(ε) for another), explaining NNZ mismatch but not values.

The 0.037% RHS error is consistent with ~4 shared fault faces out of 9,312 total
(0.04%). Amplified by condition number κ ≈ 160 → 6.4% displacement error.

### 5.3 Key MFEM Code Paths

- Face transformation: `mesh.cpp:540-561` — PointMat columns from `faces[FaceNo]->GetVertices()`
- Face-to-element mapping: `eltrans.cpp:606` — eip1/eip2 through Loc1/Loc2
- Interior face assembly: `bilinearform.cpp:672-696` — GetInteriorFaceTransformations
- Shared face assembly: `pbilinearform.cpp:229-272` — GetSharedFaceTransformations
- Slip RHS: `elasticity_operator.hpp:AssembleSlipContributionIPShared`

### 5.4 Attempted Fixes (No Effect)

| Fix | Commit | Effect |
|-----|--------|--------|
| Owned-fault DOF layout | f94d000 | None — only 4 shared fault faces |
| Canonical DOF permutation | 978155d | None — permutation doesn't change face parameterization |
| `KeepNbrBlock(true)` | 6cf06f8 | None — doesn't change HypreParMatrix after RAP |
| `Assemble(skip_zeros=0)` | cad9a31 | None — no exact zeros in face matrix |

### 5.5 Required Fix: Full Shared-Face Reparameterization

Following Tandem's sorted-simplex approach, canonicalize the FULL face reference
mapping for shared faces before any face integrator call:

1. **Face transformation** — permute PointMat columns to canonical vertex order
2. **Loc1** — compose reference-face permutation into face→elem1 mapping
3. **Loc2** — compose reference-face permutation into face→elem2 mapping

This must be a general wrapper applied in `ParBilinearForm::AssembleSharedFaces()`
and in all custom shared-face code (slip RHS, traction recovery). NOT a PointMat-only
patch — Loc1/Loc2 control the element integration point mapping and must be consistent.

## 6. Verification Jobs

Three SBATCH scripts created for Frontera validation:

| Script | Ranks | Nodes | Purpose |
|--------|-------|-------|---------|
| `bp5_v57_mpi_test_100ranks.sbatch` | 100 | 2 | Baseline |
| `bp5_v57_mpi_test_200ranks.sbatch` | 200 | 4 | Intermediate |
| `bp5_v57_mpi_test_400ranks.sbatch` | 400 | 8 | Original scale |

All use identical parameters: 1000m mesh, IP, p=1, PETSc TS RK45, MUMPS exact, `--tandem-time-stepping`, `--checkpoint-interval 0`. Output directories indicate the rank count for easy comparison.

## 6. References

- Tandem `BoundaryMap.cpp` (fault face canonical sorting): `/Users/chunhuizhao/projects/tandem/src/form/BoundaryMap.cpp:43-51`
- Tandem `ScatterPlan` (ghost synchronization): `/Users/chunhuizhao/projects/tandem/src/parallel/ScatterPlan.h`
- MFEM shared face orientation: `/Users/chunhuizhao/projects/mfem/mesh/pmesh.cpp:2322-2547`
- Original runs: job 7627262 (100 ranks), job 7627480 (400 ranks)
