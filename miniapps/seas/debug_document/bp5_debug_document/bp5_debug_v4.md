# BP5 Debug v4: Missing Shared Face Support in ElasticityDomainOperator

## Problem

After adding the DG penalty correction to `ComputeTraction` (bp5_debug_v3), BP5 still blows up on Frontera (56 MPI ranks):
```
zeroIn: F(a) and F(b) must have different signs. a=0 F(a)=2.8e10 b=6053 F(b)=3.8e-6
```
F(a)=2.8e10 means traction ~ 28 GPa -- wildly incorrect. Happens on ranks 26, 29.

## Root Cause

The elasticity operator completely ignores shared faces in parallel. When the mesh
is partitioned across MPI ranks, fault faces on partition boundaries become "shared
faces" but the elasticity operator only tracks `fault_interior_faces_` (faces local
to a single rank). Any fault face on an MPI boundary:

1. **Not detected** -- `SetupFaultInfo()` only calls `GetInteriorFaceTransformations()`
2. **Not included in slip RHS** -- `AssembleSlipContributionIP/BR2()` only loops over `fault_interior_faces_`
3. **Not included in traction** -- `ComputeTraction()` only loops over `fault_interior_faces_`
4. **No FaultBasis** -- `FaultBasis::Compute()` uses `GetInteriorFaceTransformations()`
5. **No depth/coords** -- `GetFaultDepths()`, `GetFaultCoords2D()` miss shared faces

Result: fault DOFs on shared faces get zero traction and zero slip, causing the
friction solver to produce wildly incorrect slip rates and blowup.

The antiplane operator (`antiplane_operator.hpp`) already handles shared faces
correctly -- it has `fault_shared_faces_`, `IsFaultFaceShared()`,
`AssembleSlipContributionIPShared()`, `AssembleSlipContributionBR2Shared()`, and
shared-face loops in `ComputeTraction()`.

## Fix

### Files Modified

1. **`domain/elasticity_operator.hpp`** -- Main changes:
   - Added `fault_shared_faces_` member
   - Added `IsFaultFace3DShared()` method
   - Modified `SetupFaultInfo()` to detect shared faces
   - Added `AssembleSlipContributionIPShared()` for IP method
   - Added `AssembleSlipContributionBR2Shared()` for BR2 method
   - Updated `Solve()` to call shared-face slip assembly
   - Extended `ComputeTraction()` for shared faces (with DG penalty correction)
   - Updated `GetFaultDepths()` and `GetFaultCoords2D()` for shared faces
   - Added `GetFaultSharedFaces()` accessor

2. **`fault/fault_basis.hpp`** -- Added `AppendSharedFaces()` method:
   - Uses `GetSharedFaceTransformations()` instead of `GetInteriorFaceTransformations()`
   - Appends to existing `basis_` array computed by `Compute()`

3. **`tests/parallel/test_parallel_elasticity.cpp`** -- New parallel test file with 6 tests:
   - `test_shared_fault_detection`: Verifies shared fault faces detected
   - `test_serial_parallel_fault_dof_count`: Global count = interior + shared/2 == serial
   - `test_parallel_zero_slip_traction`: Zero slip -> zero traction
   - `test_serial_parallel_traction_consistency`: Max traction matches serial
   - `test_parallel_traction_bounded`: No NaN/Inf/blowup
   - `test_parallel_solve_with_slip`: Non-zero displacement

4. **`Makefile`** -- Added build target for parallel elasticity test

### Key Implementation Details

- Shared faces: only Elem1 is local; Elem2 is a face-neighbor accessed via
  `pfes->GetFaceNbrFE(FTr->Elem2No - mesh_.GetNE())`
- DOF values for face-neighbor elements obtained from `par_u.FaceNbrData()`
  after `par_u.ExchangeFaceNbrData()`
- Mass matrix inverses for face-neighbor elements already precomputed in
  `PrecomputeMassInverse()` (indexed by `ne + i`)
- Slip assembly for shared faces: only Elem1 contributions added to local RHS
  (the other rank handles its Elem1 independently)
- BR2 lifting: two-sided (both Elem1 and Elem2 contribute to the lifting)
  but only Elem1 RHS assembled locally

## Verification

- All 77 serial tests pass
- All 6 new parallel tests pass with 2 and 4 ranks
- Serial-parallel traction consistency: rel_err < 1e-11
