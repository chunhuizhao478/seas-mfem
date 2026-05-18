# Fix Report — REVIEW.md round 2 (R-101 through R-105)

## Summary
- Findings addressed: **5 of 5** — all LOW-severity round-2 findings fixed
- Files modified: **4**
  - `miniapps/seas/fault/fault_geometry.hpp` (R-101 docstring + R-103 fallback t2)
  - `miniapps/seas/domain/elasticity_operator_traction.inl` (R-104 dedup)
  - `miniapps/seas/drivers/project_stress_to_mesh.cpp` (R-105 pre-flight)
  - `miniapps/seas/safs/.../verify_onfault_stress.py` (R-102 mixed-block)
- Tests added: **1** new (`test_R102_mixed_cell_blocks_concatenated`)
- Test suite: **PASS — 116 / 116 C++ + 183 / 183 Python**

## Changes Made

### R-101 [LOW] — Updated stale docstring in `ComputePerDOFCoordsAndBasis_`

**File modified:** `miniapps/seas/fault/fault_geometry.hpp`

The class-level docstring previously described `t2_i ← n_i × t1_i`
(the pre-R-001 buggy convention). Updated to reflect the post-fix
Gram-Schmidt-of-input-t2 algorithm, with an inline note explaining
WHY the cross-product form is unsafe (loses FaultBasis sign-flip) and
referencing the R-001 contract.

### R-102 [LOW] — Mixed cell-block resilience in Phase 8 verifier

**File modified:** `verify_onfault_stress.py`, `test_verify_onfault_stress.py`

Added `_cell_centroids_with_field(m, name)` that aligns centroids with
the cell blocks where the named field is actually defined. Updated
`_cell_data_array` to skip `None` / empty entries. Both
`verify_bulk_cell_data` and `verify_fault_cell_data` now use the new
helper inside the per-field loop, so mixed-cell-block meshes
(tetra + triangle in a bulk VTU) verify correctly instead of crashing
on shape mismatch.

`verify_fault_cell_data` additionally falls back to a clear stderr
diagnostic when the fault VTU's per-field cell count does not match
the triangle-only fault geometry (the per-cell basis is undefined on
non-triangle blocks, so verification is skipped with a warning).

Added regression test `test_R102_mixed_cell_blocks_concatenated`:
mixed tetra+triangle bulk VTU with `sigma_xx_Pa` populated on both
blocks must verify successfully against the analytic prediction.

### R-103 [LOW] — Re-derive t2 from cross product in the t1-degeneracy fallback

**File modified:** `miniapps/seas/fault/fault_geometry.hpp`

In the degenerate-t1 fallback (when the projected input t1 collapses
below 1e-12), the code now sets `used_t1_fallback = true` and re-derives
`t2 = n × t1_fallback` instead of Gram-Schmidt-projecting the (also
likely degenerate) input t2. The FaultBasis sign convention is
undefined for this DOF anyway, so sign loss is acceptable and the
cross product cannot collapse to zero by construction (n and
t1_fallback are unit and mutually orthogonal).

This makes the fallback path actually usable in the degenerate input
case it was meant to handle — previously the subsequent
`t2_len > 1e-12` MFEM_VERIFY would have aborted on the degenerate
input. The fallback remains unreachable in practice (FaultBasis aborts
earlier on the same degeneracy), so no test regression is required.

### R-104 [LOW] — Extracted per-face helpers from `GetFaultDOFCoords3D` / `GetFaultDOFBasis`

**File modified:** `miniapps/seas/domain/elasticity_operator_traction.inl`

Extracted two anonymous-namespace inline helpers (`WriteFaceDOFCoords3D_`,
`WriteFaceDOFBasis_`) that contain the per-face DOF iteration body
previously duplicated between the interior-face and shared-face loops
of each function. Each loop body now calls one helper. The interior /
shared branch retains its own face accessor lookup (`GetInteriorFaceTransformations`
vs `GetSharedFaceTransformations`) but the per-DOF math lives in
exactly one place per accessor.

Phase 6.A tests still pass — the refactor preserves the exact
arithmetic and storage pattern (verified via T_6A_3..T_6A_6 which
check orthonormality and FaultBasis-frame equivalence element-wise).

### R-105 [LOW, POSSIBLE] — Pre-flight bbox check before filesystem mutation

**File modified:** `miniapps/seas/drivers/project_stress_to_mesh.cpp`

Added a pre-flight block (after FE space construction, before any
ParaViewDataCollection setup) that:
1. Computes the parallel mesh bbox via per-rank min/max + MPI_Allreduce.
2. Opens the sidecar via `StressField3D probe(sidecar_path)`.
3. Calls `probe.ContainsBBox(...)`; on failure prints a clear error
   listing both bboxes and exits 1 cleanly via MPI_Finalize().

The downstream `ProjectStress` call still performs the per-component
containment check (defence in depth). Phase 7 round-trip test still
passes.

## Verification

- [x] R-101: docstring matches implementation; future maintainer
  reading the doc cannot be misled back to the buggy convention.
- [x] R-102: `_cell_centroids_with_field` aligns sizes; new test
  `test_R102_mixed_cell_blocks_concatenated` exercises a mixed
  tetra+triangle bulk VTU and passes.
- [x] R-103: fallback now sets `t2 = n × t1_fallback` and the
  subsequent Gram-Schmidt projection on the orthonormal frame is a
  no-op modulo FP.
- [x] R-104: helper extracted; T_6A_3..T_6A_6 still pass to 1e-12.
- [x] R-105: pre-flight bbox check added; Phase 7 round-trip test
  still passes; containment-failure path now aborts before
  filesystem mutation.

## Test Suite Results

C++ (linked against /Users/chunhuizhao/projects/seas-mfem/libmfem.a +
openblas; run with DYLD_LIBRARY_PATH):

| Suite | Tests | Passed |
|-------|-------|--------|
| Phase 6 Tranche 1 (`seas_test_stress_field_3d`) | 56 | 56 |
| Phase 6.A (`seas_test_fault_dof_basis`) | 18 | 18 |
| Phase 6 §4 (`seas_test_project_fault_prestress`) | 12 | 12 |
| Phase 6 §5 (`seas_test_compute_safs_params`) | 13 | 13 |
| Phase 6 §6 (`seas_test_safs_mode_wiring`) | 5 | 5 |
| Phase 6 §7 (`seas_test_safs_stress_config`) | 9 | 9 |
| Phase 7 (`seas_test_project_stress_to_mesh`) | 3 | 3 |
| **Total C++** | **116** | **116** |

Python (`pytest` in `pythonenv`):

| Suite | Tests | Passed |
|-------|-------|--------|
| `test_project_to_fault_stress.py` (Phase 1-4) | 129 | 129 |
| `test_build_stress_safs.py` (Phase 5) | 44 | 44 |
| `test_verify_onfault_stress.py` (Phase 8, +1 R-102 test) | 10 | 10 |
| **Total Python** | **183** | **183** |

## New Tests
- `test_R102_mixed_cell_blocks_concatenated` (Python) — exercises the
  R-102 fix; mixed tetra+triangle bulk VTU with `sigma_xx_Pa` on both
  blocks must verify without shape mismatch.

## Ready for Re-Review: YES

All round-2 LOW findings are addressed; total test count is 299 with
zero failures. The codebase passes the full suite of plan-defined
Phase 6 (Tranches 1, 2, 3), Phase 7, and Phase 8 unit tests.
