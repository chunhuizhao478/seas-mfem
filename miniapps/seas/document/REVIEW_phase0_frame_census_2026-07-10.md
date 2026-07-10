# Code Review: Phase 0 — fault frame-bit census — 2026-07-10

## Review Scope
- **Plan:** `miniapps/seas/document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md`, Phase 0
- **Files reviewed:**
  - `tests/parallel/test_fault_frame_bit_census.cpp` (new)
  - `Makefile` (new vars, obj rule, link rule, `test-fault-frame-bit-census`)
- **Domain context:** `CLAUDE.md` (sign conventions, antiparallel-traction blow-up),
  `dynamic/wave_operator.inl` (FaultBasis construction, R-1601 fallback),
  `fault/fault_basis.hpp`, `debug_document/tpv104_debug_document/SUBSTEP_ITERATOR_MPI_FIX.md`.
- **Context:** the census RAN and **falsified** the plan's root cause. Phase 1 is withdrawn.
  Defect (1) — shared fault QPs get a one-shot macro-`dt` solve while interior QPs get O sub-steps —
  **survives** and now has no oracle. That materially raises the severity of R-001.

---

## Findings

### [R-001] MODERATE [Phase 0 deliverable (b)] — The RED parity test was never written; the surviving defect has no oracle

**Category:** DEVIATION

**Description:**
Phase 0 required **two** deliverables. Only (a) was implemented.
- (a) the census diagnostic — DONE (`test_fault_frame_bit_census.cpp`).
- (b) `tests/parallel/test_shared_fault_substep_parity_np2.cpp` — **MISSING.**

Phase 0's acceptance criteria list "The np=1 vs np=2 comparison FAILS today (this is the point)".
No such test exists. This was tolerable while the frame-bit hypothesis was live, because Phase 1 was
going to be the fix. Now that the hypothesis is dead, **(b) is the only artifact that pins the
defect which actually survives**, and it is the oracle any future fix must be measured against.

**Trigger:** `make test-parallel` — nothing exercises the interior-vs-shared flux divergence.

**Actual behavior:** The repository asserts nothing about defect (1).

**Expected behavior:** A test that demonstrates the sub-step buffer is honored on interior fault
faces and *ignored* on shared fault faces.

**Suggested fix:** A full np=1-vs-np=2 physics comparison needs a macro-step harness and
cross-partition QP matching. A **strictly cheaper and more decisive** test exercises the exact
mechanism: install a *sentinel* sub-step buffer and observe whether the corrector consumes it.

```cpp
// tests/parallel/test_shared_fault_substep_parity_np2.cpp  (new)
//   Install a sentinel I_imp buffer, run AdvanceADER, compare Q_new against a
//   run with NO buffer installed.
//     np=1 (fault face INTERIOR): Q_new MUST differ  -> buffer consumed.
//     np=2 (fault face SHARED)  : Q_new is IDENTICAL -> buffer DISCARDED (R-1601).
// The second assertion is the defect, characterised.  A third, env-gated
// assertion (SEAS_EXPECT_SHARED_FAULT_PARITY=1) demands the buffer BE consumed on
// shared faces; it fails today and is the RED oracle for any future fix.
```

Reuse `BuildTwoTetSharedFault()` / `PartitionByYSign()` from the census (extract them to a shared
header — see R-006). Required public API: `WaveOperator::SetFaultDOFData`, `SetFaultFlux`,
`SetFaultFrictionLaw(FaultFrictionLaw::LSW)`, `SetSubStepFaultImposedStates`,
`ResetSubStepFaultImposedStates`, `AdvanceADER`. Seed `DOFData` with non-zero
`lsw_mu_s/lsw_mu_d/lsw_d_c/sigma_n0/tau2_0` and impedances, or `EvaluateADER_LSW`'s own
`MFEM_VERIFY` ("all LSW-native fields are zero") will abort.

**Test case:**
```cpp
// Sentinel-buffer discrimination.  `flat` = NUM_STATE * GetNumTotalFaultQPs().
std::vector<real_t> sentinel(flat, 1.0e-3);
Vector Q_no_buf(Q.Size()), Q_with_buf(Q.Size());
wave.ResetSubStepFaultImposedStates();
wave.AdvanceADER(Q, dt, O, Q_no_buf);
wave.SetSubStepFaultImposedStates(sentinel.data(), sentinel.data(), n_total_fault_qps);
wave.AdvanceADER(Q, dt, O, Q_with_buf);
const real_t d = MaxAbsDiff(Q_no_buf, Q_with_buf);
if (nprocs == 1) { ASSERT(d > 0.0,  "interior fault QPs CONSUME the substep buffer"); }
else             { ASSERT(d == 0.0, "shared fault QPs DISCARD it (R-1601 fallback) -- DEFECT"); }
// env-gated RED oracle:
if (getenv("SEAS_EXPECT_SHARED_FAULT_PARITY") && nprocs == 2)
{ ASSERT(d > 0.0, "RED: shared fault QPs SHOULD consume the substep buffer"); }
```

---

### [R-002] MODERATE [test_fault_frame_bit_census.cpp:RunCensus] — Short-circuit `&&` leaves the second canonical normal unnormalized, and C2 runs anyway on garbage

**Category:** BUG

**Description:**
```cpp
c.have_shared = Normalize3(c.can_n_from_sign_flipped)
                && Normalize3(c.can_n_from_elem1_on_plus);
```
`&&` short-circuits. If the first vector is degenerate, `Normalize3` returns `false` and the second
is **never normalized**. Worse, `main()` computes the C2 dot products **unconditionally** — they are
guarded by `if (g_nprocs == 2)`, not by `have_shared`. A failed `CENSUS_ASSERT(c.have_shared, ...)`
does not stop execution, so the two `dot` assertions then compare zero/unnormalized vectors and emit
misleading failures instead of a clean "degenerate face" diagnosis.

**Trigger:** any shared fault face whose per-QP canonical normals sum to (near) zero — e.g. a
sign-bimodal face, which is exactly the pathology this census exists to detect. The census would
then report a *wrong* answer on the one input it most needs to get right.

**Actual behavior:** `can_n_from_elem1_on_plus` stays unnormalized; `dot_e1` is meaningless; the
C2 assertions fail for the wrong reason.

**Expected behavior:** Normalize both unconditionally; skip C2 (with an explicit diagnostic) when
either normal is degenerate.

**Suggested fix:**
```diff
-      if (si == 0 && !bd.qp_data.empty())
-      {
-         for (int d = 0; d < 3; ++d)
-         {
-            c.can_n_from_sign_flipped[d]  = acc_sf[d];
-            c.can_n_from_elem1_on_plus[d] = acc_e1[d];
-         }
-         c.have_shared = Normalize3(c.can_n_from_sign_flipped)
-                         && Normalize3(c.can_n_from_elem1_on_plus);
-      }
+      if (si == 0 && !bd.qp_data.empty())
+      {
+         for (int d = 0; d < 3; ++d)
+         {
+            c.can_n_from_sign_flipped[d]  = acc_sf[d];
+            c.can_n_from_elem1_on_plus[d] = acc_e1[d];
+         }
+         // Evaluate BOTH (no short-circuit): a degenerate face must not leave
+         // one vector silently unnormalized.
+         const bool ok_sf = Normalize3(c.can_n_from_sign_flipped);
+         const bool ok_e1 = Normalize3(c.can_n_from_elem1_on_plus);
+         c.have_shared = ok_sf && ok_e1;
+      }
```
and in `main()`:
```diff
-      CENSUS_ASSERT(c.have_shared, "np=2: shared-face canonical normals built");
+      CENSUS_ASSERT(c.have_shared, "np=2: shared-face canonical normals built");
+      if (!c.have_shared)
+      {
+         if (g_rank == 0)
+         { std::cerr << "  [C2] SKIPPED: degenerate face-averaged normal.\n"; }
+      }
+      else
+      {
+         ... existing Allgather + dot assertions ...
+      }
```

**Test case:**
```cpp
// R-002: a degenerate accumulator must not silently leave the 2nd vector unnormalized.
real_t a[3] = {0,0,0}, b[3] = {0,2,0};
const bool ok_a = Normalize3(a);
const bool ok_b = Normalize3(b);       // must still run
ASSERT(!ok_a && ok_b, "both Normalize3 calls evaluated");
ASSERT(std::abs(b[1] - 1.0) < 1e-15, "second vector WAS normalized");
```

---

### [R-003] MODERATE [test_fault_frame_bit_census.cpp:RunCensus] — Shared basis-index formula is correct but never exercised by the fixture

**Category:** ASSUMPTION

**Description:**
```cpp
const int basis_idx = int_faces.Size() + si;   // shared face si
```
This is **correct** — `wave_operator.inl:358` calls `fault_basis_->Compute(mesh, fault_interior_faces_, ...)`
and then appends shared faces (`AppendSharedFaces`), so basis indices run interior-first. But the
fixture has `interior=1, shared=0` at np=1 and `interior=0, shared=1` at np=2. **The offset term
`int_faces.Size()` is always zero when it is used.** A wrong formula (e.g. `si` alone, or
`shr_faces.Size() + si`) would pass this test identically. The census's core measurement rests on an
index convention it does not verify.

**Trigger:** any mesh where a rank owns both interior and shared fault faces — i.e. every real mesh.

**Actual behavior:** untested; a future change to FaultBasis ordering silently corrupts the census.

**Expected behavior:** assert the convention, and exercise a mixed fixture.

**Suggested fix (1) — cheap assertion, do this at minimum:**
```diff
    const FaultBasis *fb = op.GetFaultBasis();
    MFEM_VERIFY(fb != nullptr, "census: WaveOperator has no FaultBasis");
+   // FaultBasis is built interior-first (wave_operator.inl:358 Compute(fault_interior_faces_))
+   // then AppendSharedFaces.  The shared basis index below depends on that layout.
+   MFEM_VERIFY(fb->NumFaces() == int_faces.Size() + shr_faces.Size(),
+               "census: FaultBasis face count (" << fb->NumFaces()
+               << ") != interior (" << int_faces.Size() << ") + shared ("
+               << shr_faces.Size() << ") -- basis-index convention changed");
```
(hoist `int_faces`/`shr_faces` above this check).

**Suggested fix (2) — exercise the mixed case:** add a second fixture using
`BuildFourHexTwoFault()` (two y=0 fault faces) with a partition that leaves **one fault face
interior and one shared** on rank 0, so `int_faces.Size() > 0` *and* `shr_faces.Size() > 0`
simultaneously. Assert `c.n_interior_faces > 0 && c.n_shared_faces > 0` on that rank.

**Test case:**
```cpp
// R-003: on a rank owning BOTH face classes, the shared basis index must be offset.
ASSERT(c.n_interior_faces > 0 && c.n_shared_faces > 0, "mixed fixture: both classes present");
ASSERT(fb->NumFaces() == c.n_interior_faces + c.n_shared_faces, "interior-first basis layout");
// and C1 on the interior subset must still be 0, while the shared subset disagrees on one rank.
```

---

### [R-004] MODERATE [Makefile] — The STOP-gate test is in no aggregate, so it never runs

**Category:** DEVIATION

**Description:**
`test-fault-frame-bit-census` exists but appears in **no** aggregate target. `Makefile:5007` defines
`test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel
test-parallel-fault test-shared-fault-reconcile-cross-rank test-bimaterial-deriv-cache-parity` — the
census is absent, and it is not in `SEQ_MINIAPPS` either. Plan Phase 0 required registration.

This matters more than ordinary test hygiene: the census's C2 assertions are the **regression guard**
that stops a future contributor from "unifying" the frame bit and reintroducing the `4e28` overflow.
A guard that never executes is not a guard.

**Trigger:** `make test-parallel`, or any CI invocation.

**Suggested fix:**
```diff
-test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel test-parallel-fault test-shared-fault-reconcile-cross-rank test-bimaterial-deriv-cache-parity
+test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel \
+               test-parallel-fault test-shared-fault-reconcile-cross-rank \
+               test-bimaterial-deriv-cache-parity test-fault-frame-bit-census
```
Note `seas_test_bimaterial_seam_fault_np2` is likewise absent from `SEQ_MINIAPPS`, so parallel tests
are conventionally **not** listed there. Registering in `test-parallel` is the correct and sufficient
action; do **not** add it to `SEQ_MINIAPPS`.

**Test case:** `make -n test-parallel | grep -c seas_test_fault_frame_bit_census` returns ≥ 1.

---

### [R-005] LOW [test_fault_frame_bit_census.cpp:main] — Pass/test counters are per-rank but the failure count is globally reduced

**Category:** BUG

**Description:** The summary prints rank 0's local `g_num_passed` / `g_num_tests` alongside a
globally reduced `failed_global`. A failure occurring only on rank 1 yields the self-contradictory
line `Phase 0 census: 5 / 5 passed, 1 failed`. The exit code is correct; the report is not.

**Trigger:** any rank-asymmetric failure (e.g. the C1 disagreement lives on exactly one rank, so
this is a realistic path).

**Suggested fix:**
```diff
    int failed_global = 0;
    MPI_Allreduce(&g_num_failed, &failed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
+   int passed_global = 0, tests_global = 0;
+   MPI_Allreduce(&g_num_passed, &passed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
+   MPI_Allreduce(&g_num_tests,  &tests_global,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (g_rank == 0)
    {
-      std::cout << "  Phase 0 census: " << g_num_passed << " / " << g_num_tests
+      std::cout << "  Phase 0 census: " << passed_global << " / " << tests_global
                 << " passed, " << failed_global << " failed\n";
```

---

### [R-006] LOW [test_fault_frame_bit_census.cpp] — Fixture duplicated verbatim; unused include

**Category:** QUALITY

**Description:** `BuildTwoTetSharedFault()` and `PartitionByYSign()` are copy-pasted from
`tests/parallel/test_bimaterial_seam_fault_np2.cpp`. R-001's new test will need them a third time.
`#include "../../dynamic/wave_state.hpp"` is unused.

**Suggested fix:** extract to `tests/parallel/shared_fault_fixtures.hpp` (header-only) and include it
from all three tests. Drop the unused include. Do this as part of R-001 so the duplication does not
triple.

---

### [R-007] LOW [test_fault_frame_bit_census.cpp:main] — C2 silently skipped at np > 2, yielding a vacuous pass

**Category:** EDGE_CASE

**Description:** `if (g_nprocs == 2)` guards C2. Running `mpirun -np 4 ./seas_test_fault_frame_bit_census`
exits 0 having asserted almost nothing. The `test-fault-frame-bit-census` target only invokes np=1
and np=2, but a developer probing at np=4 gets false reassurance.

**Suggested fix:**
```diff
+   MFEM_VERIFY(g_nprocs == 1 || g_nprocs == 2,
+               "census fixture is defined for np=1 (interior fault face) or "
+               "np=2 (shared fault face); got np=" << g_nprocs);
```

---

## Summary
- Critical issues: **0**
- Moderate issues: **4** (R-001, R-002, R-003, R-004)
- Low issues: **3** (R-005, R-006, R-007)
- **Plan compliance: PARTIAL** — Phase 0 (a) done and green; (b) **missing**; (c) partially done
  (built + run target, but registered in no aggregate); (d) not applicable to (a), pending on (b).
- **Verdict: PASS WITH FIXES.** The census is sound enough that its headline result stands: it
  measured `dot = +1` for `sign_flipped` and `dot = -1` for `!elem1_on_plus`, on a face where
  `have_shared` was true and `int_faces.Size() == 0` (so R-002 and R-003 could not have corrupted
  it). The falsification of the plan's root cause is therefore **safe to rely on**. R-001 must be
  closed before any further phase, because defect (1) now has no oracle.

## Unreviewed Areas
- The census was run only on the 2-tet fixture. `tpv104_symmirror_1000m.msh` and
  `tpv2053d_200m.msh` are gitignored and **absent locally**, so the plan's intended census meshes
  were never exercised. The y-mirror/FP-bimodality question is therefore untested; only the
  rank-locality question was answered.
- `wave_operator.inl`'s RHS scatter of `I_imp_plus/I_imp_minus` onto Elem1/Elem2 was not audited.
  It is the remaining candidate location for the true cause of the R-1601 overflow.
- No production code was changed in Phase 0 (`git status` clean under `dynamic/ drivers/ spatial/
  fault/`), so no regression surface was introduced.
