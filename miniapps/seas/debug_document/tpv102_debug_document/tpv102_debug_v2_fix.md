# Fix Report: tpv102_debug_v2

Follow-up to `tpv102_debug_v2_check.md`.  All nine code-level findings
addressed; four dev-queue sbatch scripts and one MPI regression test
added per the check-doc recommendations.

## Summary

- Findings addressed: 9 / 9 (R-101 through R-109)
- Files modified:
  - `dynamic/wave_operator.hpp`
  - `dynamic/wave_operator.inl`
  - `drivers/tpv102_driver.cpp`
  - `io/paraview_output.hpp`
  - `Makefile`
- Files added:
  - `tests/parallel/test_r101_shared_fault.cpp` (new MPI regression test)
  - `jobs/tpv102/tpv102_1000m_p1_1.5s_4rank_dev.sbatch`
  - `jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`
  - `jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`
  - `jobs/tpv102/tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch`
- Tests added: 1 new test source file (two test cases).
- Test suite: PASS
  - `seas_test_fault_face_flux` — 19 / 19
  - `seas_test_tpv102_setup`    — 24 / 24
  - `seas_test_wave_operator`   — 17 / 17
  - `seas_test_parallel_wave_operator` (4 ranks) — 5 / 5
  - `seas_test_tpv102_local`    — 16 / 16
  - `seas_test_r101_shared_fault` (2 ranks) — SKIPPED (see below)

## Changes Made

### R-101 [CRITICAL] Verify shared-fault DOFData consistency at runtime
- `dynamic/wave_operator.hpp` — new method
  `VerifySharedFaultDOFDataConsistency(real_t tol=1e-10)` (public const).
- `dynamic/wave_operator.inl` — implementation: packs every shared-fault
  QP's `(cx, cy, cz, tau1_corr, V1, psi)` into a buffer, `MPI_Allgatherv`s
  across ranks, sorts entries by centroid, and asserts bit-equality on
  paired entries.  On mismatch, calls `MFEM_ABORT` with the offending
  centroid + field + delta.  All ranks participate in the collective
  even when locally empty (fixed a latent deadlock in the first draft).
- `drivers/tpv102_driver.cpp::main` — invoke the check once after the
  first RK4 step so any invariant failure is loud and fatal instead of
  silent drift.

### R-102 [MODERATE] Hard-abort on unresolved peer_rank
- `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` — replaced
  the `(peer_rank < 0) || my_rank < peer_rank` owner check with
  `MFEM_VERIFY(shared_face_peer_[sf].resolved, ...)` +
  `MFEM_VERIFY(peer_rank != my_rank_, ...)`.  A ctor resolution failure
  now aborts with a descriptive message instead of silently returning
  to the pre-R-001 buggy path.

### R-103 [MODERATE] Direct stdlib includes
- `drivers/tpv102_driver.cpp` — added direct `#include <limits>`,
  `#include <vector>`, `#include <cstddef>` (for `offsetof`).

### R-104 [MODERATE] Single-shot schedule gate
- `io/paraview_output.hpp` — new `CommitSchedule(cycle, time, V_max)`
  that advances `last_write_time_` without re-evaluating the gate.
- `drivers/tpv102_driver.cpp::paraview_write` — replaced the
  `PeekShouldWrite`-then-`Save`/`ShouldWrite` double-evaluation with
  `PeekShouldWrite` → (pack + write via ForceSave or CommitSchedule) →
  `WriteFaultSurfaceVTU`.  The gate is now evaluated exactly once per
  call and the mutation (advance last_write_time_) is always explicit.

### R-105 [MODERATE] Named MinDist struct + static_asserts for MPI_DOUBLE_INT
- `drivers/tpv102_driver.cpp::main` — replaced the anonymous
  `struct { double d; int r; }` with a named `MinDist` struct that
  carries `static_assert(offsetof(...) == 0)` and
  `static_assert(offsetof(r) == sizeof(double))`, catching any future
  layout deviation at compile time.  Also eliminates the silent
  `real_t → double` widening in the prior initialiser.

### R-106 [LOW] Non-collapsing qnorm watch list
- `drivers/tpv102_driver.cpp::main` — replaced the one-sided
  `{hypo_rank, hypo+1, hypo+4, nprocs-1}` with bilateral neighbours
  plus rank 0 and `nprocs-1`.  For `hypo_rank` at either end of the
  range the watch now still shows 3–4 ranks including the far side,
  so the diagnostic stays useful.  Verified at runtime: the 4-rank
  coarse run shows `r0=1.8e-27 r1=2.0e-10 r2=6.9e-24 r3=1.2e-12`
  at t=1.0 s (hypo_rank=1).

### R-107 [LOW] Separate `resolved` from `peer_rank`
- `dynamic/wave_operator.hpp` — replaced `std::vector<int>
  shared_face_peer_rank_` with `std::vector<SharedFacePeer>
  shared_face_peer_` where `SharedFacePeer { bool resolved; int
  peer_rank; }`.  R-102's `MFEM_VERIFY` checks `resolved` explicitly,
  distinguishing "ctor lookup failed" from "peer_rank == 0".

### R-108 [LOW] const correctness
- Addressed by R-104's `PeekShouldWrite` (const) + `CommitSchedule`
  (non-const) split.  Read and write are now lexically distinct.

### R-109 [LOW] Cache my_rank_ in the operator
- `dynamic/wave_operator.hpp` — new `int my_rank_ = 0` member.
- `dynamic/wave_operator.inl` (ctor) — populated once from
  `pmesh.GetMyRank()`.
- `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` — use
  `my_rank_` instead of re-querying per Mult call.

## New Tests

`tests/parallel/test_r101_shared_fault.cpp` — 2-rank MPI test that
loads `tpv102/mesh/tpv102_1000m.msh`, runs one RK4 stage and then
ten RK4 steps, and calls `VerifySharedFaultDOFDataConsistency` after
each.  Pass criterion: the runtime check does not abort.  The test
prints a per-rank fault-layout summary so SKIP conditions are
explicit.

Makefile targets:
- `make seas_test_r101_shared_fault`
- `make test-r101-shared-fault` (runs with `-n 2`)

## Unresolved Findings / Discussion

### R-101 unit-test caveat — MFEM mesh-partitioning observation

During test development I discovered that MFEM does NOT classify
TPV102-style fault faces (generated by Gmsh's `BooleanFragments`,
which duplicates vertices on either side of the fault surface) as
*shared* faces after `ParMesh` partitioning.  Each rank's portion of
the fault appears in its `GetNSharedFaces()` count as zero and in
its `GetNBE()` list as local interior-BE faces.

Concretely, at 2, 4, and 8 MPI ranks on `tpv102_1000m.msh`:
```
rank 0: local=2307 shared=0
rank 1: local=2631 shared=0
(global total: 4938 local fault QPs, 0 shared fault QPs)
```

Implication for R-001:
- The `ComputeSharedFaceFluxRHS` shared-fault code path is NOT
  reached for these meshes — the fault is handled entirely by
  `ComputeFaceFluxRHS`'s interior-fault branch.
- R-001's (+,-) swap is defensive code for a code path that never
  fires on the current meshes, and thus cannot be the cause of the
  v1 "rupture-stops-at-partition-seam" symptom.
- The v1 symptom was therefore (most likely) fixed by **R-005**
  (deep-copy of `nbr_data`) in the ordinary cross-rank bulk-wave
  path, not by R-001.
- The R-101 runtime diagnostic is still valuable as a long-range
  guard: if a future mesh format or MFEM version classifies fault
  faces as shared, the diagnostic will fire on every real run.

The R-101 unit test consequently SKIPs (reporting "no shared fault
faces on this partitioning") on all tested mesh + rank combinations.
Fully exercising the shared-fault path would require a mesh with a
fault surface that is NOT duplicated at the vertex level — not
currently available in the repo.

### Graded-mesh antisymmetry test (v2 rec #3)
Not added as a separate test: the existing `test_parallel_wave_operator`
covers generic cross-rank propagation, and the 4-rank TPV102 coarse run
(reproduced in this session, see v1 fix report) exercises a graded-mesh
partition in practice.  Adding a dedicated synthetic graded-mesh test
is a nice-to-have but duplicates coverage.

### Visual ParaView smoke (v2 rec #4)
Shipped as `jobs/tpv102/tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch`.
The sbatch comments describe the exact visual-inspection procedure.

## Verification

- [x] R-101: runtime diagnostic implemented; driver calls it after step 0;
      all ranks participate in the collective.  Dedicated 2-rank test
      written but SKIPs with current meshes (see caveat above).
- [x] R-102: MFEM_VERIFY abort on unresolved peer; second VERIFY on
      `peer_rank != my_rank_`.
- [x] R-103: `<limits>`, `<vector>`, `<cstddef>` included directly.
- [x] R-104: single `PeekShouldWrite` gate; `ForceSave` / `CommitSchedule`
      perform the advance; `WriteFaultSurfaceVTU` only runs when the gate
      passed.
- [x] R-105: named `MinDist` with two `static_assert`s.
- [x] R-106: bilateral watch list; r0 / r(nprocs-1) always present.
- [x] R-107: `SharedFacePeer { bool resolved; int peer_rank; }` struct.
- [x] R-108: const-correct split via R-104's new helper.
- [x] R-109: `my_rank_` member populated once in ctor.

## Frontera Dev-Queue Workflow (v2 rec #1–#3 + visual smoke)

Run in this order; do NOT submit the 400-rank job until both cheaper
jobs pass cleanly.

1. `sbatch jobs/tpv102/tpv102_1000m_p1_1.5s_4rank_dev.sbatch` (~15 min)
2. `sbatch jobs/tpv102/tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch`
   (~15 min, optional visual smoke — useful if R-101 diagnostic SKIPs
   and you want a sanity check via ParaView)
3. `sbatch jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch` (~2 hr)
4. `sbatch jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch` (~1.5 hr
   within the 2 hr dev-queue window) — this is the dispositive test
   reproducing the v1 symptom configuration.

Each job writes `--debug-qnorm` output; look for:
- `[qnorm:watch]` showing ||Q||_∞ ramping up on the hypocenter rank
  first, then its neighbours within ~100 ms of nucleation.
- `[R-101 check]` line: either "consistency OK: N pairs matched"
  or the diagnostic is dormant (shared-fault branch not reached).
  An `MFEM_ABORT("R-101 shared-fault DOFData consistency FAILED…")`
  would be a hard-fail.

## Ready for Re-Review: YES
