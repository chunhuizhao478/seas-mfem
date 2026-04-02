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

## 4. Known Limitation: Missing DOF Permutation

### 4.1 MFEM Shared Face Orientation

MFEM's `GetSharedFaceTransformations` returns face transformations where the face vertex ordering depends on the local element (elem1). The same reference integration point can map to **different physical coordinates** on the two ranks sharing a face.

Evidence from MFEM source (`mesh/pmesh.cpp`):
- `shared_trias[sf].v` stores canonical (global-vertex-ID-sorted) ordering
- `faces[FaceNo]->GetVertices()` uses local element-dependent ordering
- Orientation info is tracked in `Elem1Inf`/`Elem2Inf` via `GetTriOrientation()`

### 4.2 Impact on the Fix

The `ExpandOwnedToLocalFault` communication copies DOF values directly by index (`kk=0,1,2`) without applying orientation-based permutation. If the face vertex ordering differs between owner and ghost, DOF values are placed at wrong physical locations on the ghost rank.

**Expected impact:** For flat faces with parameters that are nearly constant within a face, the error is negligible. For faces crossing parameter transition zones (VW/VS boundary, nucleation zone edge), the error could be significant but smaller than the original duplicate-evolution bug.

### 4.3 Validation Plan

Run the same BP5 v57 configuration at three rank counts (100, 200, 400) with the fix applied:

- If V_max agrees to ~1e-10 → permutation issue is negligible
- If ~0.1-1% difference persists → need to add `GetTriOrientation`-based DOF permutation

## 5. Verification Jobs

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
