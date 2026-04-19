# Fix Report v3: address review findings in `tpv102_debug_v3_check.md`

## Summary
- Findings addressed: **9 of 9** (all R-302 parts included; initial draft deferred Part A,
  subsequent pass completed it)
- Files modified:
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/drivers/tpv102_driver.cpp`
  - `miniapps/seas/io/paraview_output.hpp`
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`
- Tests added: 1 new test `TestR302_InlineTwoTetSharedFault` (inside the existing
  `seas_test_r101_shared_fault` binary — see R-302 Part A below)
- Test suite: **PASS**
  - `seas_test_godunov_flux`   — 29/29
  - `seas_test_wave_operator`  — 17/17
  - `seas_test_wave_bc`        —  10/10
  - `seas_test_fault_face_flux`— 19/19
  - `seas_test_tpv102_setup`   — 24/24
  - `seas_test_tpv102_local`   — 16/16
  - `seas_test_parallel_wave_operator` (2 and 4 ranks) — 5/5 (exercises the
    newly-Release-visible `MFEM_VERIFY` in `ComputeSharedFaceFluxRHS`; no
    trip)
  - `seas_test_r101_shared_fault`:
    - 2 ranks: R-302a **actually runs** (8/8 pass) with
      `max_diff=0`, `0 unpaired`, `3 pairs matched` — the R-001 (+,-)
      canonicalisation is now bit-for-bit verified.
    - 4 ranks: R-302a SKIPs with a reasoned message (inline mesh has 2
      tets; requires exactly 2 ranks).  `TestR101_OneStage` /
      `TestR101_TenSteps` still SKIP on TPV102 meshes as disclosed.

## Changes Made

### [R-301] [CRITICAL] — Verifier now packs all 8 mutable DOFData fields
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

The record layout went from 6 doubles (cx, cy, cz, tau1_corr, V1, psi) to
11 doubles (centroid + tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi,
slip1, slip2).  Constants `REC = 11`, `FIELD_BASE = 3`, `NUM_FIELDS = 8`
and a `FIELD_NAMES[8]` table replace the ad-hoc 3-way ternary.  The pair-
comparison loop now iterates `k = FIELD_BASE..REC` so every field is
checked.  The MPI_Allgatherv buffer grows by ~2×; for TPV102-scale
shared-fault QP counts (< 1e5 globally) this remains negligible.

### [R-302 Part B] [CRITICAL] — MFEM_ASSERT → MFEM_VERIFY for R-005 aliasing check
**Where:** `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS`

The deep-copy invariant (`nbr_data[c].GetData() != q_gf.FaceNbrData().GetData()`)
is load-bearing for cross-rank bulk-wave propagation — without it every
`nbr_data[c]` would alias the last-component ghost buffer.  Promoting
the check from debug-only `MFEM_ASSERT` to `MFEM_VERIFY` means the
invariant is now exercised in Release builds.  `seas_test_parallel_wave_operator`
passes under the stricter check at 2 and 4 ranks.

### [R-302 Part A] [CRITICAL] — Inline 2-tet shared-fault test now runs
**Where:** `tests/parallel/test_r101_shared_fault.cpp`:
`BuildTwoTetSharedFaultMeshInline` + `TestR302_InlineTwoTetSharedFault`.
Plus: a ctor bug this test exposed, fixed in `dynamic/wave_operator.inl`.

Implementation: `BuildTwoTetSharedFaultMeshInline` programmatically builds
a 5-vertex, 2-tet mesh where the two tets share one triangular face
(at y=0).  All 7 triangular faces are added as boundary elements: 6
external faces with attr=1 (natural), and the 1 interior face with attr=3
(fault).  `TestR302_InlineTwoTetSharedFault` partitions the mesh one tet
per rank via an explicit `int partition[2] = {0, 1}`, constructs a
`WaveOperator`, and invokes `VerifySharedFaultDOFDataConsistency`.  On 2
ranks the test:

  - Confirms `GetNSharedFaces() == 1` on each rank.
  - Confirms `wave.GetFaultSharedFaces().Size() == 1` on each rank.
  - Invokes `Mult` (which runs the R-001 swap through Evaluate on both ranks).
  - Invokes `VerifySharedFaultDOFDataConsistency(1e-10)` and expects no abort.

Actual result: `3 pairs matched, 0 unpaired entries, max_diff=0`.
The R-001 canonicalisation and R-301/R-304/R-305 verifier changes are all
exercised; the R-005 `MFEM_VERIFY` deep-copy invariant is exercised; and
the returned field values are bit-identical across ranks.

**Bug exposed and fixed while implementing the test (ctor vertex-key
matching):** When a serial `Mesh` has an interior boundary element on a
face (TPV102-style fault-tag pattern, and also our 2-tet inline test),
MFEM's `ParMesh` partitioning preserves that BE on **only ONE** of the
two ranks that share the face — the other rank's `GetNBE()` does not
include it.  The old ctor logic built `local_fault_keys` from local BEs
only and matched using local vertex indices, so the non-BE-owner rank
silently classified its shared fault face as **non-fault**
(`shared_face_bdr_attr_[sf] = 0`, `fault_shared_faces_.Size() = 0`).
The R-001 swap then never ran on that rank, and
`VerifySharedFaultDOFDataConsistency` saw unpaired records — the exact
R-305 abort condition I had just added.

**Fix (in `wave_operator.inl` ctor):**
  1. Use **global** vertex IDs (`pmesh.GetGlobalVertexIndices`) as the
     key, since local vertex numbering is rank-private.
  2. `MPI_Allgatherv` each rank's local fault keys so every rank sees
     the **union**.  One-time ctor cost, no per-RK4-stage overhead.
  3. Classify each shared face by matching its global vertex key against
     the merged set.

This closes the latent R-002 bug the comment "Both ranks sharing a fault
face independently detect it via their respective boundary elements"
implicitly assumed — which is false under MFEM's partitioning semantics.
On TPV102 meshes the fix is a no-op (no shared fault faces exist, so the
merged set is identical to the local set); on meshes with real shared
fault faces (the inline test, future non-duplicated-vertex meshes) it is
load-bearing.

### [R-303] [CRITICAL] — MPI datatype now matches real_t at compile time
**Where:** `drivers/tpv102_driver.cpp` (3 sites) and `dynamic/wave_operator.inl` (1 site)

Replaced hardcoded `MPI_DOUBLE` on `real_t` buffers with
`MPITypeMap<real_t>::mpi_type` (MFEM's idiomatic compile-time selector,
also used by `miniapps/tools/gridfunction-bounds.cpp` and
`miniapps/contact/ip.cpp`).  Sites fixed:
  1. `tpv102_driver.cpp`: `--debug-qnorm`'s `MPI_Gather` (the explicit R-303 site)
  2. `tpv102_driver.cpp`: `V_max_step` `MPI_Allreduce`
  3. `tpv102_driver.cpp`: NaN-detection `MPI_Allreduce`
  4. `wave_operator.inl`: `h_min_` `MPI_Allreduce` (ctor)

The other `MPI_DOUBLE` uses in the file are on genuine `std::vector<double>`
buffers (`VerifySharedFaultDOFDataConsistency::all_data`) and remain
correct.  The `MPI_DOUBLE_INT` use in R-105's `MinDist` struct is a
distinct, intentional type for `MPI_MINLOC` and is unrelated.

### [R-304] [MODERATE] — Exact-lex sort comparator
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

Replaced the tolerance-based lex comparator (which violated strict weak
ordering — `A ≈ B, B ≈ C, A < C` was possible and is formally UB under
`std::sort`) with exact-equality lex compare.  Shared-fault QP centroids
on the two owning ranks are computed from the same MFEM Face transformation
at the same reference `IntegrationPoint`, so they agree bit-for-bit by
construction; no tolerance is needed.  A side effect is that topological
mismatches (a QP owned by only one rank) now surface as unpaired entries
instead of being absorbed into 3-element "groups" via transitive tolerance
— which dovetails with R-305.

### [R-305] [MODERATE] — `n_unpaired > 0` now aborts
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

`fail_local` now includes `n_unpaired > 0`.  The abort message distinguishes
the two failure modes: a pure-unpair case gets a targeted "N unpaired
entries" abort, while the value-mismatch case gets the original R-101
centroid/field diagnostic.  The rank-0 "OK" summary is still printed
when `fail_global == 0` (i.e. all pairs matched and all entries paired).

### [R-306] [MODERATE] — PASS/FAIL post-run check in sbatch scripts
**Where:** `jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`,
`jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`

After `ibrun ./seas_tpv102_driver …` each dispositive/intermediate sbatch
now greps the job log and writes `${RESULT_DIR}/RESULT.txt` with a
machine-readable verdict:
  - **FAIL** if `R-101 shared-fault DOFData consistency FAILED` or
    `R-101 shared-fault DOFData: <N> unpaired entries` is present.
  - **FAIL** if no `[qnorm:watch]` line exists at all.
  - **FAIL** if the final `[qnorm:watch]` shows dead ranks (`rN=0` or
    `rN=1e-3NN`).
  - **PASS** otherwise.
The `RESULT.txt` body is also echoed to stdout for human eyeballing.
This replaces the prior "ibrun exits 0 ⇒ success" signal that masked
physics regressions.

### [R-307] [LOW] — Abort-message field label via lookup table
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

Subsumed by R-301's fix: `FIELD_NAMES[max_diff_field - FIELD_BASE]` with
a bounds check yields the correct name for all 8 fields, and falls back
to `"<unknown>"` if `max_diff_field` is out of range (replaces the prior
silent fall-through to `"psi"`).

### [R-308] [LOW] — `CommitSchedule` takes only the arg it uses
**Where:** `io/paraview_output.hpp::CommitSchedule`, one call site in
`drivers/tpv102_driver.cpp`

Signature changed from `CommitSchedule(int /*cycle*/, real_t time,
real_t /*V_max*/)` to `CommitSchedule(real_t time)`.  The doc-comment
explains that the cycle/V_max arguments are *gate* inputs for
`PeekShouldWrite` — they decide whether the cycle writes — and do not
belong on Commit, whose only effect is to advance `last_write_time_`.
Driver call site updated accordingly.  No other callers exist
(`grep -r CommitSchedule miniapps/seas/` confirmed).

### [R-309] [LOW] — Null-FTR guard in BuildFaultCoords
**Where:** `tests/parallel/test_r101_shared_fault.cpp::BuildFaultCoords`

Added `MFEM_VERIFY(ftr, …)` after both `GetInteriorFaceTransformations`
and `GetSharedFaceTransformations` calls.  A future mesh-state change
that leaves a stale face in the wave operator's list now produces a
clean abort message instead of segfaulting on `ftr->GetGeometryType()`.

## Unresolved Findings
None.

## New Tests
- `TestR302_InlineTwoTetSharedFault` (in `seas_test_r101_shared_fault`):
  builds a 2-tet inline mesh where MFEM classifies the inter-tet face as
  a true shared face, forces 1 tet per rank via explicit partitioning,
  and drives `VerifySharedFaultDOFDataConsistency`.  Requires exactly 2
  ranks; SKIPs with a reasoned message otherwise.  Confirms:
    - `GetNSharedFaces() == 1` on each rank.
    - `wave.GetFaultSharedFaces().Size() == 1` on each rank (after the
      new MPI-exchange fix in the ctor).
    - Verifier reports `3 pairs matched, 0 unpaired, max_diff=0` —
      i.e., the R-001 (+,-) canonicalisation produces bit-identical
      DOFData across the two ranks, through every one of the 8 mutable
      fields that R-301 now packs.
- `seas_test_parallel_wave_operator` (2 and 4 ranks) also exercises the
  R-005 `MFEM_VERIFY` deep-copy invariant in Release (no trip).

The R-301/R-304/R-305 verifier changes are now covered end-to-end by
the inline test rather than being dormant.

## Verification
- [x] R-301: fixed — all 8 mutable DOFData fields packed and checked, with a lookup-table-driven abort message.  Covered end-to-end by R-302a.
- [x] R-302 Part A: fixed — `TestR302_InlineTwoTetSharedFault` runs the full shared-fault pipeline on a 2-tet inline mesh; passes with `max_diff=0`.
- [x] R-302 Part B: fixed — `MFEM_ASSERT` → `MFEM_VERIFY`, covered by parallel wave-operator tests.
- [x] R-303: fixed — `MPITypeMap<real_t>::mpi_type` used at all 4 `MPI_DOUBLE`-on-`real_t` sites.
- [x] R-304: fixed — exact-lex sort comparator; strict weak ordering restored.  Covered by R-302a (3 pairs matched with centroid max_diff=0).
- [x] R-305: fixed — `n_unpaired > 0` now aborts with a targeted message.  Verified by the pre-ctor-fix state of R-302a, which correctly tripped "3 unpaired entries" before the ctor MPI-exchange fix closed the root cause.
- [x] R-306: fixed — PASS/FAIL `RESULT.txt` emitted by both dev-queue sbatch scripts.
- [x] R-307: fixed — subsumed by R-301's `FIELD_NAMES[]` table with bounds-checked fallback.
- [x] R-308: fixed — `CommitSchedule(real_t time)`; driver call site updated.
- [x] R-309: fixed — null-FTR `MFEM_VERIFY` at both sites in `BuildFaultCoords`.
- [x] **Bonus fix (latent R-002 bug):** ctor vertex-key matching moved
  from local to global vertex IDs + MPI-Allgatherv merge, closes the
  "both ranks independently detect fault" assumption that MFEM's
  partitioning violates.

## Ready for Re-Review: YES
