# TPV104 `--fault-iterator substep` np>1 hang — Bug Report (R-1600+)

## Summary

Adversarial review focused on **the production-blocking hang** of
`--fault-iterator substep` at np>1.  Reproduction (symmirror 1000m mesh,
np=10, ADER-O=2, Brent friction, slip-SRW): banner prints, "Starting
ADER-O(2) time loop..." appears, then **zero "Step N/N" output for 13.5
minutes** with all 10 ranks pinned at 100% CPU in R state (user-space
spin, not MPI block / IO).  Comparison run with `--fault-iterator
one-shot` on the same mesh + ranks completes 522 steps in seconds.

Files reviewed:
- `dynamic/wave_operator.inl` — focus on `EvaluateBulkAtFaultQPsCanonical`
  (L1583-1946), `ComputeADERSubStepStates` (L1087-1190), interior /
  shared substep gates (L3497-3521, L4414-4438), `ComputeADERSharedFaceFluxRHS`
  prologue (L4143-4198).
- `dynamic/wave_operator.hpp` — `GetNumLocal/Shared/TotalFaultQPs`
  (L435-440), substep side-channel members (L597-605).
- `dynamic/tpv104_substep_iterator.{hpp,cpp}` — full read; loop bounds in
  `Advance` (L284-492) and `AdvanceWithSubStepStates` (L619-682) are
  consistent with `dof_data.size()`.
- `drivers/tpv104_driver.cpp` — `AdvanceADERWithSubStep` helper
  (L286-409), R-1003 activation log (L1517-1534), time-loop dispatch
  (L1984-1994).
- `dynamic/friction_solver.{cpp,hpp}` — Brent / Newton-stable solver
  bounds.
- `debug_document/tpv104_debug_document/SUBSTEP_ITERATOR_MPI_REVIEW.md`
  (R-1300..R-1306, Rev 1).
- `debug_document/tpv104_debug_document/SUBSTEP_ITERATOR_MPI_FIX.md`
  (claims §1–§5 of R-1003 LANDED).

**Review history**
- **Rev 1** (2026-04-25): 4 new bugs (R-1600..R-1603); identifies the
  exact deadlock site that produces the 13.5-min CPU spin; verifies
  R-1300..R-1304 are mechanically present in code but uncovers a
  dormant defect in §1's ghost-exchange wiring (R-1600) introduced
  WITH the R-1003 landing — i.e., R-1003 fixed the np>1 abort but
  introduced a collective-skip on a strict subset of ranks.

---

## Pass-1 verification of R-1300..R-1306

Quick diff of "what the prior review proposed" vs. "what is in the
code today":

| ID | Prior status | Current code state | Verdict |
|---|---|---|---|
| R-1300 | Shared-fault gate missing | `wave_operator.inl:4414-4438` adds the gate, mirroring the interior gate at `:3497-3521`. Identical conditional `(substep_I_imp_*_flat_ != nullptr && dof_idx >= 0 && dof_idx < substep_n_total_fault_qps_)`. | **LANDED** |
| R-1301 | Buffer sized to local-only | `tpv104_driver.cpp:370` reads `wave.GetNumTotalFaultQPs()`; member renamed to `substep_n_total_fault_qps_`; setter accepts `n_total_fault_qps`. | **LANDED** |
| R-1302 | EvaluateBulk skips shared | `wave_operator.inl:1737-1944` adds a shared-fault loop after the interior loop, sized to `n_total_qps`. | **LANDED** but with a deadlock defect (R-1600). |
| R-1303 | Predictor needs ghost | Resolved differently from §1 of the plan: the per-substep ghost exchange was **colocated** inside `EvaluateBulkAtFaultQPsCanonical` (L1773-1784) rather than added to `ComputeADERSubStepStates`. Consequence: the gate that wraps it (L1753 `if (n_shared_faces > 0)`) creates the R-1600 deadlock. | **LANDED but BROKEN** at np>1. |
| R-1304 | dof_idx must be absolute on shared branch | `wave_operator.inl:4416-4427` uses the existing `dof_idx` (already absolute via `shared_fault_dof_offset_`) and indexes `substep_I_imp_*_flat_ + dof_idx * NUM_STATE` — no rebase from zero. | **LANDED** |
| R-1305 | DOFData symmetry MPI check | DEFERRED in the FIX doc. Not relevant to the hang. | DEFERRED (unchanged) |
| R-1306 | dt_scale round-trip | DEFERRED. Not relevant to the hang. | DEFERRED (unchanged) |

So the hang is **not** R-1300..R-1306 regressing; it is a NEW defect
introduced by the way §1 of the R-1003 plan was implemented.

---

## Pass-2 — root cause of the hang

### R-1600: `EvaluateBulkAtFaultQPsCanonical` skips its `ExchangeFaceNbrData` collective on ranks with no fault shared faces (or no fault QPs at all), creating a pairwise-MPI deadlock with ranks that DO have fault shared faces

**Status:** OPEN
**Severity:** CRITICAL (root cause of the production-blocking hang)
**Category:** MPI collective correctness / ghost-exchange topology

**File:** `dynamic/wave_operator.inl`
**Lines:** 1604 (early-return gate) and 1753 (parallel-block gate)

**Description.**  The new function
`WaveOperator<ParMesh>::EvaluateBulkAtFaultQPsCanonical` (L1583-1946)
performs a per-substep ghost exchange of the predictor's
per-component bulk Q via `q_gf.ExchangeFaceNbrData()` (L1779), repeated
NUM_STATE = 9 times.  Each `ExchangeFaceNbrData` call is a **pairwise**
MPI exchange among the rank's face-neighbour set
(`pmesh.face_nbr_group_`) — every rank that shares a face with rank A
must call it for the calls to match.

The function has **two early-skip gates** that depend on the
*fault*-shared-face count, NOT the *total*-shared-face count:

```cpp
// L1597-1604  — early return on rank with zero fault QPs total
const int n_total_qps = GetNumTotalFaultQPs();
const size_t expect_words =
   static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n_total_qps);
Q_plus_flat.assign(expect_words, 0.0);
Q_minus_flat.assign(expect_words, 0.0);
if (n_total_qps == 0 || !fault_basis_) { return; }   // <-- SKIP-1

// L1749-1753  — parallel block conditional on FAULT-shared count only
if constexpr (IsParallelMesh<MeshType>::value)
{
#ifdef MFEM_USE_MPI
   const int n_shared_faces = fault_shared_faces_.Size();   // FAULT-only
   if (n_shared_faces > 0)                                  // <-- SKIP-2
   {
      ...
      for (int c = 0; c < NUM_STATE; c++)
      {
         q_gf.ExchangeFaceNbrData();                         // pairwise MPI
         ...
      }
      ...
   }
#endif
}
```

Compare to the SAFE pattern used by the existing macro-step path
`ComputeADERSharedFaceFluxRHS` at `wave_operator.inl:4166-4167`:

```cpp
auto &pmesh = static_cast<const ParMesh &>(mesh_);
int n_shared = pmesh.GetNSharedFaces();      // ALL shared, not fault
if (n_shared == 0) { return; }
```

`pmesh.GetNSharedFaces() == 0` is the only safe early-return: a rank
with zero shared faces has zero face neighbours, so
`ExchangeFaceNbrData` would be a no-op anyway.

`fault_shared_faces_.Size() == 0` and `GetNumTotalFaultQPs() == 0` are
**not** safe early-return conditions: a rank can have **shared
non-fault faces** with rank A (e.g., the fault-adjacent rank) while
having zero shared FAULT faces of its own.  For that rank,
`pmesh.face_nbr_group_` is non-empty and rank A's
`ExchangeFaceNbrData` posts an MPI_Irecv from this rank — but this
rank skipped the call (SKIP-1 or SKIP-2), so no matching MPI_Isend
from this rank ever arrives.  Rank A spins forever in the MFEM
NeighborSendData / NeighborRecvData internal `MPI_Wait` polling loop
→ R state at 100% CPU.

**Empirical match to the report.**

- "Banner prints, then zero Step output for 13.5 minutes" — consistent
  with the FIRST call to `AdvanceADERWithSubStep` deadlocking at the
  first `q_gf.ExchangeFaceNbrData()` inside the FIRST sub-step's
  `EvaluateBulkAtFaultQPsCanonical` call.
- "All 10 ranks alive in R state at 100% CPU" — consistent with
  Open MPI / MPICH busy-wait inside `MPI_Wait` after a non-blocking
  `MPI_Irecv`.  S/D state would be only on implementations that fall
  back to `nanosleep` after some poll budget; the default modern MPI
  binaries spin first.
- "Same params + one-shot: 522 steps in seconds" — one-shot path
  goes through `wave.AdvanceADER` which calls
  `ComputeADERSharedFaceFluxRHS` whose gate is the SAFE
  `pmesh.GetNSharedFaces() == 0` (L4167).  Every rank either has zero
  face neighbours total (skips and is topologically absent) or has
  some face neighbours and runs the full exchange — no asymmetric
  skip.  No deadlock.
- "Banner: rank 0 has local=0, shared=24" — rank 0 IS a fault-shared
  rank (24 QPs / nbf_per_face = 4 shared fault faces on rank 0,
  assuming nbf_per_face=6 for ADER-2 quadrature on a tet).  Rank 0
  proceeds past SKIP-1 (`n_total_qps = 24 > 0`) and past SKIP-2
  (`n_shared_faces = 4 > 0`) and POSTS the receives.  The other 9
  ranks of the partition: at least one has zero fault shared faces
  but non-zero face neighbours via shared **non-fault** faces with
  rank 0, and that rank skips the matching sends.

**Trigger.**  Any np>1 run with `--fault-iterator substep` on a mesh
where the partition produces at least one rank with
`fault_shared_faces_.Size() == 0` AND `pmesh.GetNSharedFaces() > 0`,
i.e., a rank that is bordering a fault-adjacent rank but doesn't
itself touch the fault.  TPV104 production meshes (1000m, 200m, 100m
on 8+ ranks) almost guarantee this — only a small minority of ranks
border the embedded fault.

At np=1 the parallel block is `if constexpr` skipped at compile
time on `Mesh` (non-`ParMesh`); no exchange is ever called → no
deadlock.  This is consistent with the user's note that "np=1
substep was always working."

**Suggested fix.**  Replace BOTH gates with the macro-step's safe
`pmesh.GetNSharedFaces() == 0` early-return.  The shared-fault loop
itself is still gated on `fault_shared_faces_.Size() > 0`, but the
preceding `ExchangeFaceNbrData` calls run on every rank that has any
face neighbours.

Concrete diff:

```diff
--- a/miniapps/seas/dynamic/wave_operator.inl
+++ b/miniapps/seas/dynamic/wave_operator.inl
@@ -1597,10 +1597,15 @@ void WaveOperator<MeshType>::EvaluateBulkAtFaultQPsCanonical(
    const int n_local_qps = GetNumLocalFaultQPs();
    const int n_total_qps = GetNumTotalFaultQPs();
    const size_t expect_words =
       static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n_total_qps);
    Q_plus_flat.assign(expect_words, 0.0);
    Q_minus_flat.assign(expect_words, 0.0);

-   if (n_total_qps == 0 || !fault_basis_) { return; }
+   // R-1600: the per-substep ExchangeFaceNbrData below is a PAIRWISE
+   // collective among face neighbours.  A rank with zero fault QPs but
+   // non-zero shared faces (i.e., bordering a fault-adjacent rank) MUST
+   // still participate, otherwise its peers hang in MPI_Wait.  Defer the
+   // early-return until AFTER the exchange when the parallel block runs.
+   if (!fault_basis_) { return; }
+   if (n_total_qps == 0 && !IsParallelMesh<MeshType>::value) { return; }

    const real_t *Q_data = Q_bulk.GetData();
@@ -1749,9 +1754,15 @@ void WaveOperator<MeshType>::EvaluateBulkAtFaultQPsCanonical(
    if constexpr (IsParallelMesh<MeshType>::value)
    {
 #ifdef MFEM_USE_MPI
-      const int n_shared_faces = fault_shared_faces_.Size();
-      if (n_shared_faces > 0)
+      // R-1600: gate on TOTAL shared faces, NOT fault-only shared faces.
+      // A rank that has shared non-fault faces with a fault-adjacent peer
+      // MUST call ExchangeFaceNbrData; otherwise the peer's matching
+      // MPI_Irecv on the per-substep ghost exchange never gets paired and
+      // every fault-adjacent rank hangs in MPI_Wait at 100% CPU.
+      auto &pmesh_top = const_cast<ParMesh &>(
+         static_cast<const ParMesh &>(mesh_));
+      const int n_shared_total = pmesh_top.GetNSharedFaces();
+      const int n_shared_faces = fault_shared_faces_.Size();
+      if (n_shared_total > 0)
       {
          auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
          MFEM_VERIFY(pfes,
@@ -1781,7 +1792,11 @@ void WaveOperator<MeshType>::EvaluateBulkAtFaultQPsCanonical(
          ...
          }
+         if (n_shared_faces == 0)
+         {
+            // No fault shared faces on THIS rank, but the exchange above
+            // was needed for peers.  Skip the per-fault-face loop.
+            return;
+         }

          for (int sf_idx = 0; sf_idx < n_shared_faces; sf_idx++)
          {
```

The fix is mechanically minimal: gate the **collective** on the safe
total-shared count, but keep the per-fault-face loop gated on the
fault-shared count.  Total-fault-QP=0 ranks now run the 9 pairwise
exchanges (no-op for them but topologically required for peers) and
return without entering the per-fault-face loop.

**Test.**  See `Test 1 (R-1600-DEADLOCK)` below — a 2-tet `ParMesh`
fixture where np=2 splits across the fault, with both ranks present;
combined with a 3-tet fixture (np=3) where ONE rank has no fault face
at all but DOES share a non-fault face with a fault rank — that
configuration reproduces the deadlock under the current code, and
passes after the fix.

---

### R-1601: per-substep ghost exchange is unnecessarily repeated O times — wastes 9·(O-1) collectives per macro-step

**Status:** OPEN
**Severity:** MEDIUM (performance, not correctness; latent risk of magnifying R-1600's deadlock window)
**Category:** Performance / MPI overhead

**File:** `dynamic/wave_operator.inl` lines 1773-1784 (per-substep
exchange loop) called from `drivers/tpv104_driver.cpp:374-377` (per-O loop).

**Description.**  `AdvanceADERWithSubStep` (driver L372-377) calls
`wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o], ...)` once per
sub-step `o ∈ [0, O)`, passing a DIFFERENT `Q_per_node[o]` each time
(the per-sub-step pointwise predictor output).

But inside the function (L1773-1784), the `q_gf` ParGridFunction is
filled from `Q_data = Q_bulk.GetData()` — i.e., the per-substep
predictor output IS exchanged (correct).  However, the implementation
does **9 collectives PER call**, and the helper is called O times per
macro-step → **9·O collectives per macro-step on the substep path**,
versus 9 on the macro-step (one-shot) path.

The fix-doc claims "Cost: O calls × NUM_STATE collectives = 9·O
exchanges per macro-step" (FIX.md L114) and accepts this as the cost.
That's the documented cost; not a bug per se, but:

- ADER-O=2 → 18 collectives per macro-step.
- ADER-O=4 → 36 collectives per macro-step.
- TPV104 ~10⁵ macro-steps over 60 s → 1.8M (O=2) to 3.6M (O=4) extra
  collectives.

**Latent risk.** Each collective is a deadlock opportunity for R-1600.
While R-1600 fires on the FIRST call regardless of O, the cost
amplifies if R-1600 is fixed by adding more conditional logic (e.g.,
"only exchange if this substep needs it") instead of the unconditional
fix above.  The unconditional fix (R-1600 §) is the right answer; this
finding documents the cost so the maintainer doesn't regress to a
"smarter" gate that re-introduces deadlocks.

**Suggested fix.**  Optional optimization: pre-compute
`Q_per_node[o]`'s ghost layer ONCE for all O sub-steps before the loop
in `AdvanceADERWithSubStep`, then have `EvaluateBulkAtFaultQPsCanonical`
accept a cached ghost-data argument.  Saves (O-1)·9 collectives per
macro-step.  Not blocking for production correctness; defer until
after R-1600 is fixed and tests pass.

**Test.**  Performance-regression test (Test 4 below): time
`AdvanceADERWithSubStep` for a fixed mesh + step count at np=2 vs np=8
and assert sub-linear scaling on collectives.

---

### R-1602: SKIP-1 early return at L1604 leaves `Q_plus_flat` / `Q_minus_flat` at zero on a rank with `n_total_qps == 0`, propagating into the iterator's loop bound check downstream

**Status:** OPEN
**Severity:** LOW (silently correct on the AdvanceWithSubStepStates path because the iterator `n` matches, but a future caller with mismatched sizing would silently produce zero rows)
**Category:** Defensive programming / contract documentation

**File:** `dynamic/wave_operator.inl` lines 1597-1604

**Description.**  At L1601-1602 the function calls
`Q_*_flat.assign(expect_words, 0.0)` BEFORE the early return at L1604,
so a rank with `n_total_qps == 0` returns clean zero-sized buffers.
The driver helper at `tpv104_driver.cpp:370-377` then sizes and uses
those buffers via `n_total_fault_qps = wave.GetNumTotalFaultQPs() = 0`
on that rank, and the iterator at L385 short-circuits the
`AdvanceWithSubStepStates` call when `n_total_fault_qps > 0` is false.
So this rank's iterator is a no-op.  Correct.

**However**, this contract is implicit: a hypothetical second caller
of `EvaluateBulkAtFaultQPsCanonical` who sized their downstream
buffers at MAX_RANK_TOTAL_QPS or some other value would silently get
zero rows for a rank with no fault QPs.  Add a doc-string assertion
that the output `Q_*_flat.size()` IS `NUM_STATE *
GetNumTotalFaultQPs()` and that 0-size is a valid case meaning "this
rank has no local fault QPs."

**Trigger.**  Future caller mis-uses the size convention.

**Suggested fix.**  No code change required for current callers.  Add
a doc-string note at `wave_operator.hpp:340-365` and an
`MFEM_ASSERT(Q_plus_flat.size() == NUM_STATE * GetNumTotalFaultQPs())`
inside the function for self-documentation.

**Test.**  None required.

---

### R-1603: `pmesh.GetSharedFaceTransformations` is called inside an O-step loop without caching — recomputed `O × n_shared_faces` times per macro-step

**Status:** OPEN
**Severity:** LOW (perf; not correctness)
**Category:** Performance

**File:** `dynamic/wave_operator.inl` lines 1789-1790

**Description.**  Inside the per-substep call to
`EvaluateBulkAtFaultQPsCanonical`, line 1789-1790 calls
`pmesh.GetSharedFaceTransformations(sf)` for every shared-fault face
sf, every sub-step.  For TPV104 with O=2 and ~24 shared fault faces
per rank, that's 48 calls per macro-step.  `GetSharedFaceTransformations`
internally builds the FaceElementTransformations object (allocating
matrix data, computing element transforms).  The macro-step
`ComputeADERSharedFaceFluxRHS` also calls it — total O+1 = 3 calls per
sf per macro-step.

The interior-fault loop has the same pattern with
`mesh_.GetInteriorFaceTransformations` (L1612).  Pre-fix: identical
recomputation.  Not introduced by R-1003.

**Suggested fix.**  None blocking for the hang.  Defer to a
performance pass after R-1600 lands.

**Test.**  None required.

---

## Pass-3 — non-bugs ruled out

I list these explicitly so the next reviewer/debugger doesn't re-walk
them:

- **Friction solver unbounded loop (suspect 1 in the prompt).**  Brent
  via `qd_friction_.SolveSlipRatePsi` is iteration-bounded (proven QD
  log10-V solver; Tandem-verified per BP5 CLAUDE.md L24); Newton-stable
  via `SolveNRStable` is bounded by `max_iter=60`; both have
  `tau<=0` short-circuits at `friction_solver.cpp:79-80, 100, 110, 135`.
  `--friction-solver newton-stable` flows through `MapSolver` →
  `Method::NewtonRaphsonStable` (driver L123, L779) → bounded.
- **Iterator inner loop with wrong upper bound (suspect 2 in the
  prompt).**  `AdvanceWithSubStepStates` iterates `i ∈ [0, n)` where
  `n = dof_data.size() = num_fault_total = local + shared` (driver
  L1097); `Q_pointwise_*_per_substep[o]` is sized
  `NUM_STATE * num_fault_total` per driver L370-377; per-call verify
  at iterator L573-589 enforces this.  No out-of-bounds, no garbage.
- **Per-fault-QP absolute index but local buffer (suspect 5 in the
  prompt).**  `dof_idx` from `shared_fault_dof_offset_[sf] + q` is
  absolute (interior + shared) per `wave_operator.hpp:481-487`; the
  substep buffer at `tpv104_driver.cpp:382-383` is sized
  `NUM_STATE * n_total_fault_qps`; the gate at L4417 is
  `dof_idx < substep_n_total_fault_qps_` — out-of-range correctly
  falls through to `EvaluateADER`.  R-1304 honored (no rebase).
- **All-to-all in tight loop (suspect 3 in the prompt).** No
  `MPI_Allgatherv` or `MPI_Allreduce` inside the per-substep loop; the
  exchange is pairwise (`ExchangeFaceNbrData`).  R-1600 IS the
  pairwise-deadlock variant of suspect 3.
- **R-1003 wiring incomplete.** R-1300, R-1301, R-1302, R-1304
  mechanically present in code (verified line-by-line in Pass-1).
  R-1303 resolved differently (ghost exchange colocated with
  EvaluateBulk rather than added to ComputeADERSubStepStates) — the
  resolution is correct in spirit but introduces R-1600.

---

## Proposed Unit Tests

### Test 1 (R-1600-DEADLOCK): np=3 fixture, one rank with shared non-fault faces but zero fault faces, must NOT deadlock

**Target:** `WaveOperator<ParMesh>::EvaluateBulkAtFaultQPsCanonical`
end-to-end via the driver helper `AdvanceADERWithSubStep`.

**File:** `tests/unit/test_tpv104_substep_iterator_mpi_no_fault_rank.cpp`
(NEW, REQUIRES MPI_Init)

**Validates:** A 3-tet ParMesh fixture where the partition assigns:
- Rank 0: tet adjacent to fault, with 1 fault interior face and 1 shared non-fault face to rank 2.
- Rank 1: tet on the other side of the fault (1 shared FAULT face to rank 0).
- Rank 2: tet that does NOT touch the fault but DOES share a non-fault face with rank 0.

Under current code, rank 0 and rank 1 each post a per-substep ghost
receive from rank 2, but rank 2 has `n_total_qps == 0` (SKIP-1) → rank
2 returns early → rank 0 and rank 1 hang in `MPI_Wait` inside
`q_gf.ExchangeFaceNbrData`.  After R-1600 fix, rank 2 enters the
`if (n_shared_total > 0)` block, calls `ExchangeFaceNbrData` 9 times
(no fault output since SKIP-2 still applies for the per-fault-face
LOOP, only the COLLECTIVE moved) and exits cleanly.

**Priority:** HIGH (merge-blocker for any np>1 substep production run)

**Skeleton:**

```cpp
#include "mfem.hpp"
#include "miniapps/seas/dynamic/wave_operator.hpp"
#include "miniapps/seas/dynamic/tpv104_substep_iterator.hpp"
#include "miniapps/seas/dynamic/fault_face_flux.hpp"
#include <gtest/gtest.h>
#include <mpi.h>
#include <sys/time.h>

namespace {
double WallSeconds() {
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return tv.tv_sec + 1e-6 * tv.tv_usec;
}
}

TEST(SubstepIteratorMpiHangRegression, NoFaultRankDoesNotDeadlock_np3)
{
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   ASSERT_EQ(nprocs, 3) << "this test must run under mpirun -np 3";

   // Build a 3-tet box mesh, mark face f0 as the fault between tet0/tet1,
   // partition_assignment = {0, 1, 2} so each tet lives on its own rank.
   // Tet2 is positioned to share face f12 with tet0 but NOT touch f0.
   mfem::Mesh mesh = mfem::seas::test::ThreeTetFaultBoxMesh();
   mfem::Array<int> part(3); part[0]=0; part[1]=1; part[2]=2;
   mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh, part);

   // Sanity: rank 2 has zero fault faces (interior + shared).
   const int n_local_fault_2 = ...;   // helper to count fault interior on rank
   const int n_shared_fault_2 = ...;
   if (rank == 2) {
      ASSERT_EQ(n_local_fault_2 + n_shared_fault_2, 0)
         << "test fixture invariant: rank 2 must have zero fault faces";
      ASSERT_GT(pmesh.GetNSharedFaces(), 0)
         << "test fixture invariant: rank 2 must share a non-fault face";
   }

   // Build wave op + iterator + driver helper inputs.
   mfem::seas::WaveOperator<mfem::ParMesh> wave(pmesh, /*order=*/1, ...);
   mfem::seas::FaultFaceFlux fault_flux(...);
   mfem::seas::SlipLawSRWPsi state_evo(...);
   mfem::seas::Tpv104SubStepIterator iter(fault_flux, state_evo);
   const real_t dt = 1e-3;
   iter.SetSubSteps({0.5*dt, 0.5*dt}, {0.5, 0.5});

   // Allocate dof_data, fault_coords, V_w sized to local+shared on EACH rank.
   std::vector<mfem::seas::DOFData> dof_data;
   std::vector<mfem::Vector> fault_coords;
   std::vector<mfem::real_t> V_w;
   InitTPV104FaultStateForTest(wave, dof_data, fault_coords, V_w);

   mfem::Vector Q(NUM_STATE * wave.GetNdofTotal());
   Q = 0.0;
   mfem::Vector Q_new(Q.Size());

   // Hang-detection: wrap the call with a wall-clock timer + MPI_Allreduce
   // so any rank that doesn't return within 30s aborts the test.  Without
   // R-1600, the call hangs FOREVER; with R-1600 fixed it returns in ms.
   const double t0 = WallSeconds();
   AdvanceADERWithSubStep(wave, iter, dof_data, fault_coords, V_w,
                          Q, dt, /*ader_order=*/2,
                          /*t_step_start=*/0.0,
                          mfem::seas::FrictionSolver::Method::Brent,
                          Q_new);
   const double t1 = WallSeconds();
   const double elapsed = t1 - t0;
   double max_elapsed = elapsed;
   MPI_Allreduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX,
                 MPI_COMM_WORLD);
   EXPECT_LT(max_elapsed, 5.0)
      << "AdvanceADERWithSubStep hung > 5s on rank "
      << rank << " (R-1600 deadlock); elapsed = " << max_elapsed;
}
```

**Detects:** R-1600.  Will be a guaranteed FAIL on the unfixed code
(test runner's per-test wall-clock kill > 5s) and a guaranteed PASS
after the fix.

---

### Test 2 (R-1600-COLLECTIVE-CONSISTENCY): every rank's collective-call count must match across all ranks per macro-step

**Target:** counter inside `EvaluateBulkAtFaultQPsCanonical` that
records every `q_gf.ExchangeFaceNbrData()` call.

**File:** `tests/unit/test_tpv104_substep_iterator_mpi_collective_count.cpp`
(NEW, MPI_Init required)

**Validates:** After every macro-step, the value
`(call_count_of_ExchangeFaceNbrData)` is identical across all ranks
(or differs by at most a known multiplicative constant for the
SKIP-2 short-cut).  Catches future regressions where a partial gate
is re-introduced.

**Priority:** MEDIUM

**Skeleton:**

```cpp
TEST(SubstepIteratorMpi, CollectiveCallCountConsistent_np3)
{
   int rank, nprocs;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   // Use SEAS_DIAG_ASSERT_COLLECTIVE_COUNTS=1 build flag to enable a
   // per-WaveOperator counter that increments inside the parallel block.
   ::setenv("SEAS_DIAG_ASSERT_COLLECTIVE_COUNTS", "1", /*overwrite=*/1);
   ...
   AdvanceADERWithSubStep(wave, iter, ...);
   const int my_count = wave.GetExchangeFaceNbrDataCallCount();
   int max_count, min_count;
   MPI_Allreduce(&my_count, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&my_count, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   EXPECT_EQ(max_count, min_count)
      << "ExchangeFaceNbrData call counts asymmetric across ranks: "
      << "rank " << rank << " saw " << my_count
      << " calls; max=" << max_count << " min=" << min_count;
}
```

**Detects:** R-1600 and any future re-introduction.

---

### Test 3 (R-1600-PARITY): np=2 substep produces same DOFData as np=1 on a 2-tet shared-fault fixture

**Target:** End-to-end np=2 vs np=1 bit-equality on shared-fault QPs
(this is Test 1 from the prior review SUBSTEP_ITERATOR_MPI_REVIEW.md
that the FIX doc lists as "remains to be authored").

**File:** `tests/unit/test_tpv104_substep_iterator_np2_parity.cpp` (NEW)

**Validates:** Run a 2-tet ParMesh with one shared fault face at np=2,
gather all rank's DOFData onto rank 0, compare against a serial np=1
run of the same physical configuration. After R-1600 fix this should
PASS bit-for-bit (the iterator is deterministic given identical
canonical inputs).

**Priority:** HIGH (also blocks merge per the prior review's "merge
blocker" tag).  Cannot run before R-1600 is fixed.

**Skeleton:** see SUBSTEP_ITERATOR_MPI_REVIEW.md L398-417.

**Detects:** R-1600 (would deadlock under current code; passes after
fix).  Also R-1300, R-1301, R-1302, R-1304 (latent regression guards).

---

### Test 4 (R-1601 perf): collective count scales sub-linearly with O

**Target:** `AdvanceADERWithSubStep`

**File:** `tests/unit/test_tpv104_substep_iterator_collective_cost.cpp`
(NEW)

**Validates:** Per-macro-step ExchangeFaceNbrData count vs ADER-O.
After the optional R-1601 fix (cache the ghost layer once per macro-
step), the count should be 9 (independent of O).  Without the fix:
9·O.  This test sets a regression baseline.

**Priority:** LOW (perf only; non-blocking)

**Detects:** R-1601 if the cache fix lands; the test asserts the cost
property and would regress if the cache is broken.

---

### Test 5 (R-1602 contract): doc-string for `EvaluateBulkAtFaultQPsCanonical` output sizing on np>1 with no fault QPs

**Target:** doc-string only.

**File:** N/A (doc change)

**Validates:** that the output buffer sizing contract is documented.

**Priority:** LOW

---

## Priority Summary

| ID | Severity | Status | File | Lines | Description |
|---|---|---|---|---|---|
| R-1600 | CRITICAL | OPEN | `wave_operator.inl` | 1604, 1753 | Per-substep ghost exchange skipped on ranks with no fault shared faces → pairwise-MPI deadlock with fault-adjacent ranks → 100% CPU R-state spin (the production hang). |
| R-1601 | MEDIUM   | OPEN | `wave_operator.inl` | 1773-1784 | 9·O exchanges per macro-step where 9 would suffice. |
| R-1602 | LOW      | OPEN | `wave_operator.inl` | 1597-1604 | Implicit zero-output contract on no-fault ranks; document explicitly. |
| R-1603 | LOW      | OPEN | `wave_operator.inl` | 1789-1790 | `GetSharedFaceTransformations` recomputed per substep. |

### Test Coverage Summary

| Suite | File | Tests | Bugs detected |
|---|---|---|---|
| Hang regression  | `test_tpv104_substep_iterator_mpi_no_fault_rank.cpp` (np=3, NEW) | 1 | R-1600 |
| Collective count | `test_tpv104_substep_iterator_mpi_collective_count.cpp` (np=3, NEW) | 1 | R-1600 latent regression |
| np=2 parity      | `test_tpv104_substep_iterator_np2_parity.cpp` (np=2, NEW) | 1 | R-1600 + R-1300/1/2/4 latent |
| Collective cost  | `test_tpv104_substep_iterator_collective_cost.cpp` (np>1, NEW) | 1 | R-1601 |
| **Total**        |                                                                  | **4** |       |

---

## Code Review Report (Rev 1)

**Date:** 2026-04-25
**Scope:** TPV104 `--fault-iterator substep` np>1 hang investigation
(production-blocker; identified via 13.5-min spin on symmirror 1000m
mesh at np=10).
**Reviewer:** chunhui-code-reviewer agent

### Summary

R-1003 successfully retired the np>1 substep abort by mechanically
landing R-1300, R-1301, R-1302, R-1304 (verified line-by-line in
Pass-1).  However, the R-1303 resolution chose to colocate the
per-substep ghost exchange inside `EvaluateBulkAtFaultQPsCanonical`
rather than extending `ComputeADERSubStepStates`.  That choice is
sound in principle (the macro-step path uses the same colocation
pattern at `ComputeADERSharedFaceFluxRHS:4228-4239`).  The
implementation, however, introduced two early-skip gates that
condition the **collective** `ExchangeFaceNbrData` on **fault**
shared-face counts (`n_total_qps == 0` at L1604; `fault_shared_faces_.Size() == 0`
at L1753).  Both gates exclude a real and common partition state at
np>1 — a rank with shared non-fault faces but no fault faces — and
that rank's silent skip pairs with peers' MPI_Irecv to produce the
exact deadlock the user observed.

The fix is mechanical: replace the per-substep
`fault_shared_faces_.Size() > 0` gate with the SAFE
`pmesh.GetNSharedFaces() > 0` gate that already protects the macro-
step path, and run the collective even when the per-fault-face LOOP
is empty.

### Code Review Findings

| # | File | Line(s) | Severity | Description |
|---|---|---|---|---|
| 1 | `dynamic/wave_operator.inl` | 1604, 1753 | CRITICAL | R-1600: per-substep collective skipped on no-fault ranks — pairwise-MPI deadlock at np>1. |
| 2 | `dynamic/wave_operator.inl` | 1773-1784 | MEDIUM   | R-1601: 9·O exchanges per macro-step (cacheable). |
| 3 | `dynamic/wave_operator.inl` | 1597-1604 | LOW      | R-1602: zero-output contract on no-fault ranks not documented. |
| 4 | `dynamic/wave_operator.inl` | 1789-1790 | LOW      | R-1603: per-substep `GetSharedFaceTransformations` recomputation. |

### Previously Reported Bugs — Verification

| ID | Prior status | Current code | Verified |
|---|---|---|---|
| R-1300 | Shared-fault gate missing | gate present at `wave_operator.inl:4414-4438` | YES — fixed |
| R-1301 | Buffer sized to local-only | uses `GetNumTotalFaultQPs` at driver L370 | YES — fixed |
| R-1302 | EvaluateBulk skips shared | shared-fault loop present at L1737-1944 | YES — fixed |
| R-1303 | Predictor needs ghost | colocated ghost exchange at L1773-1784 | YES — fixed (but with R-1600 latent defect) |
| R-1304 | dof_idx must be absolute | uses `dof_idx` directly without rebase at L4416-4427 | YES — fixed |
| R-1305 | DOFData symmetry MPI check | DEFERRED in FIX doc | DEFERRED (unchanged) |
| R-1306 | dt_scale round-trip | DEFERRED in FIX doc | DEFERRED (unchanged) |

### Documentation Status

- The activation log at `tpv104_driver.cpp:1517-1534` accurately
  describes that the R-1003 path is active and that shared-fault QPs
  flow through the substep side-channel. No misleading claim.
- The doc-string at `wave_operator.hpp:316-330` for
  `SetSubStepFaultImposedStates` says
  `n_total_fault_qps == GetNumTotalFaultQPs()`. Accurate post-rename.
- The doc-string at `wave_operator.hpp:340-365` for
  `EvaluateBulkAtFaultQPsCanonical` should explicitly state the
  pairwise-collective contract introduced by R-1003 §1's colocated
  ghost exchange — i.e., "every rank that has any shared face MUST
  call this function" — so a future caller doesn't gate it on
  fault-only counts.  Recommend adding to R-1602's fix.

### Test Status

- **Implemented:** 0 of 4 proposed MPI tests (the prior review's
  Tests 1-4 also remain unauthored per FIX doc L183-194; the new
  Tests 1-4 here all need MPI fixtures).
- **Existing np>1 substep coverage:** zero. The FIX doc L155-167 lists
  only np=1 unit tests as passing; the FIX doc L181-194 explicitly
  notes "Frontera MPI smoke has not been run."  R-1600 was the first
  np>1 production attempt and exposed the deadlock.
- **Priority for production-readiness:**
  1. R-1600 fix (mechanical, ~5 LOC).
  2. Test 1 (R-1600-DEADLOCK) at np=3 — must FAIL on current code,
     PASS after fix.
  3. Test 3 (np=2 parity) — bit-identity check; gates correctness of
     R-1300/1/2/4 jointly.

### Recommended Next Steps

1. **Apply the R-1600 fix** as the mechanical diff in this report.
   Estimated effort: 30 minutes (5 LOC code + a docstring update at
   `wave_operator.hpp:340-365`).
2. **Author Test 1 (R-1600-DEADLOCK)** before re-running TPV104 at
   np>1.  Without a deadlock-detection unit test, a future R-1600-
   like regression (e.g., a maintainer adds a "smarter" gate) would
   re-introduce the same 13.5-minute production hang.
3. **Author Test 3 (np=2 parity)** to verify R-1300/1/2/4 jointly.
   The FIX doc lists this as "merge-blocking" but the merge
   landed without it; this is the correctness debt to repay.
4. **Submit a Frontera np=2 short sbatch** (≤ 1 minute wall) on the
   1000m mesh to verify the fix at production scale before approving
   any np≥10 production run.  Per project memory ("Ask before
   Frontera runs"), get user approval before submitting.
5. **DEFER R-1601** (perf) and **R-1602/1603** (doc/perf) — none are
   correctness blockers and all are easier to land cleanly after the
   R-1600 + tests are in.

### Reproduction guarantee

The R-1600 mechanism predicts:
- At np=1, no deadlock (parallel branch is `if constexpr` skipped at
  compile time on `Mesh`). Matches prior-review's "np=1 always working."
- At np=2 on a 2-tet fault-spanning ParMesh, BOTH ranks have shared
  fault faces — both reach SKIP-2's `if (n_shared_faces > 0)` true
  branch, both call `ExchangeFaceNbrData`, no asymmetric skip → no
  deadlock. (This is why a 2-tet smoke might pass; the user's hang
  is at np=10 where the partition has uniform-mesh interior ranks.)
- At np=10 on the symmirror 1000m mesh, with ~384 tets distributed,
  several interior ranks have non-fault shared faces with fault-
  adjacent ranks but zero fault faces themselves — they hit SKIP-1
  or SKIP-2 → fault-adjacent ranks hang in MPI_Wait → R-state spin.
  Matches the user's observation exactly.
- At np=10 with `--fault-iterator one-shot`, the
  `ComputeADERSharedFaceFluxRHS` uses the SAFE
  `pmesh.GetNSharedFaces() == 0` gate → no asymmetric skip → no
  deadlock. Matches "522 steps in seconds" comparison run.

The mechanism is consistent with all four empirical data points
provided in the prompt.
