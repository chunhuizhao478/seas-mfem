# Known Disabled Tests

Tests that are present in the source tree but short-circuit at runtime
with a SKIP message.  Tracked here so coverage gaps are searchable.

## test_io.cpp::TestParaViewCombinedOutput

**Disabled:** 2026-05-17 (REVIEW.md R-006)

**Why:** The Phase 5 ParaViewOutput fault API (`InitFaultOutput`,
`UpdateFaultFields`) was refactored in the Phase 6 split-bulk-solutions
work into `InitFaultOutputBP5(fault_interior_faces, fault_shared_faces,
nbf_per_face)` and per-component `UpdateFaultFieldsBP5(slip, slip_rate,
traction, state, normal_stress)` using BP5-style 2-component (dip +
strike) layouts.  The test was written for antiplane (1-component slip
on a `BP2MeshGenerator` mesh), which does not map naturally onto the
BP5 fault layout.

**Coverage left intact by related tests:**
- `test_fault_surface_vtkhdf.cpp` — single-file VTKHDF fault output round trip.
- `test_fault_surface_vtkhdf_mpi.cpp` — np=4 gather-merge correctness.
- `test_fault_surface_vtu_binary.cpp` — gather-to-rank-0 binary VTU.
- `test_paraview_schedule_cap.cpp` — adaptive snapshot cap behavior.
- `test_io.cpp::TestParaViewOutput` — single ParaView PVD/VTKHDF file is created (under HDF5 build).

**Coverage gap that has no current replacement:**
- "Single combined `<prefix>.pvd` (or `<prefix>/volume.vtkhdf`) per
  driver — domain and fault fields land in the same collection" was
  the only test that asserted both end up in one collection.  All
  current tests verify domain-only OR fault-only output, not the
  combined case.

**Replacement plan (when prioritized):**
- Rewrite using a BP5 gmsh-generated 3D mesh (e.g.
  `bp5/mesh/bp5_2000m.msh`) and `ElasticityDomainOperator<Mesh>` so
  the InitFaultOutputBP5 API can be exercised with realistic
  inputs.  This adds gmsh-generated mesh data as a test prerequisite,
  which is why it has not been done as part of the merge unblock.

## (No other disabled tests at this time.)
