# Code Review: 2026-07-17 — MPI-3 shared-memory sidecar windows (Phases 0–2)

> Reviews the implemented code against
> `document/code_optimization_dev/PLAN_sidecar_mpi_shared_memory_2026-07-17.md`.
> The previous review (RK45+LSW, 2026-05-29) is archived at
> `REVIEW_rk45_lsw_impl_2026-05-29.md`.

## Review Scope
- Plan: `document/code_optimization_dev/PLAN_sidecar_mpi_shared_memory_2026-07-17.md`
- Files reviewed: `io/data_field_3d.{hpp,cpp}`, `io/stress_field_3d.{hpp,cpp}`,
  `spatial/code/spatial_velocity.{hpp,cpp}`, `spatial/code/spatial_stress.{hpp,cpp}`,
  `drivers/spatial_dyn_driver.cpp` (flag/guard/3 call sites),
  `tests/unit/test_data_field_3d_shared_mem.cpp`, `Makefile` (test wiring)
- Domain context: `miniapps/seas/CLAUDE.md` (no-refactor/no-TODO rules, MPI patterns),
  `general/error.cpp:182` (MFEM_ABORT → global MPI_Abort, verified),
  driver MPIContext comment (`MPI_Comm_free`-after-finalize abort hazard)
- Verification performed during review: re-read all final hunks; recompiled
  `io/data_field_3d.o` with warnings visible (none); ran
  `make test-data-field-3d-shared-mem` (np=2: 12/12, np=4: 24/24) and the five
  regression targets (all pass); **ran a purpose-built probe** that destroys a
  shared-mode `DataField3D` AFTER `MPI_Finalize` (the driver's real lifetime) —
  exit 0, no crash.

## Findings

### [R-001] [MODERATE] [io/data_field_3d.cpp:~DataField3D] — Whole-run windows are never freed before MPI_Finalize (plan §7 wording unachievable; guard is correct but untested in-tree)

**Category:** DEVIATION / ASSUMPTION

**Description:**
Plan §7's mitigation says "ensure sidecar readers are destroyed before
`MPI_Finalize`". That is NOT achievable for the two long-lived owners:
`vel_bundle` and `rs_sidecar` (+ the resolver's captured `shared_ptr` copies)
are main-scope locals, and the driver calls `MPI_Finalize()`
(`drivers/spatial_dyn_driver.cpp:4006` and the dry-run exit) before main's
locals unwind. They also cannot be `reset()` earlier: the R-008 lifetime
contract requires `vel_bundle` to outlive `wave_ptr`, which itself outlives
finalize. The implementation therefore skips `MPI_Win_unlock_all`/`MPI_Win_free`
when `MPI_Finalized()` is true. Per the MPI standard, finalizing with live
windows is formally erroneous; OpenMPI tolerates it (session-dir cleanup), and
the review probe confirms no crash — but the in-tree test suite never exercises
this exact path (the unit test's readers are all scoped and freed pre-finalize).

**Trigger:**
Any production run with `--sidecar-shared-mem`: material (3 windows) and
friction (2+) windows are alive at `MPI_Finalize`; their dtors run after.

**Actual behavior:**
Dtor detects `MPI_Finalized` and skips the collective free; OS reclaims the
shm segment at process exit. Verified by probe (exit 0). The transient stress
windows (6) are freed properly inside `apply_csm_impl`.

**Expected behavior:**
Same runtime behavior, but (a) the guard path must be covered by a permanent
regression test, and (b) plan §7 must state the real mitigation (guarded skip)
instead of the false "destroyed before finalize" claim.

**Suggested fix:**
1. Add the probe as a permanent single-rank test file
   `tests/unit/test_data_field_3d_shared_mem_finalize.cpp` (fixture writer +
   `new DataField3D(path, "F", Abort, MPI_COMM_SELF)` → `MPI_Finalize()` →
   `delete` → print OK), with target:
```diff
 test-data-field-3d-shared-mem: seas_test_data_field_3d_shared_mem
 	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_data_field_3d_shared_mem
 	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_data_field_3d_shared_mem
+	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_data_field_3d_shared_mem_finalize
```
2. In plan §7, replace the row "Window lifetime vs MPI_Finalize | … ensure
   sidecar readers are destroyed before MPI_Finalize …" with the actual
   mechanism: "long-lived readers outlive finalize by design; dtor skips the
   collective free under `MPI_Finalized()` (probe-tested); transient stress
   windows free normally."

**Test case:**
```cpp
// tests/unit/test_data_field_3d_shared_mem_finalize.cpp (np=1)
int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);
   write_min_fixture(path);                       // 2x2x2 grid, field "F"=2.0
   auto* f = new mfem::seas::DataField3D(path, "F",
                mfem::seas::OOBPolicy::Abort, MPI_COMM_SELF);
   if (f->Evaluate(0.5,0.5,0.5) != 2.0) { return 1; }
   MPI_Finalize();
   delete f;                                      // must not crash (guard path)
   std::puts("dtor-after-finalize OK");
   return 0;
}
```
(Ran during review as an ad-hoc probe: prints `eval=2`, `dtor-after-finalize
OK`, exit 0.)

---

### [R-002] [MODERATE] [tests/unit/test_data_field_3d_shared_mem.cpp] — No coverage of the shared-mode FAILURE path (loading-rank validation abort must kill all node ranks, not hang them)

**Category:** EDGE_CASE (test gap)

**Description:**
The load-bearing safety claim of the design — "a NaN/out-of-range/H5Dread
failure on node-rank-0 MPI_Aborts the whole job while peers wait at the publish
barrier" — is argued from `general/error.cpp:182` but never tested. The
per-rank reader has fork()-based abort tests (`test_data_field_3d.cpp`); the
shared-mode reader has none. fork() inside an MPI process is not reliable, so
the same idiom cannot be reused directly.

**Trigger:**
A sidecar with a NaN or out-of-range cell (or truncated dataset) loaded with
`--sidecar-shared-mem` on np≥2. Expected: whole job aborts promptly. A
regression that broke the abort-before-barrier ordering (e.g. moving
validation after the publish barrier, or an MFEM build where MFEM_ABORT does
not reach MPI_Abort) would instead hang every non-root rank at
`MPI_Barrier(node_comm_)` — a wall-clock-eating deadlock on a cluster.

**Actual behavior:**
Untested; correctness rests on code reading only.

**Expected behavior:**
An automated np=2 test proves: bad fixture → nonzero exit for the whole
`mpirun`, within a timeout (no hang).

**Suggested fix:**
Add a hidden self-test mode to the existing test binary and drive it from the
Makefile with a bounded, negated invocation:
```diff
 # in tests/unit/test_data_field_3d_shared_mem.cpp main(), before fixtures:
+   if (argc > 1 && std::string(argv[1]) == "--nan-fixture-abort-child")
+   {
+      // rank 0 writes a fixture whose "F" contains a NaN cell; ALL ranks
+      // then construct shared-mode readers.  EXPECTED: global abort.
+      std::string p = make_tmp_path("nanchild");
+      if (g_world_rank == 0) { write_nan_fixture(p); }
+      p = bcast_path(p, MPI_COMM_WORLD);
+      MPI_Barrier(MPI_COMM_WORLD);
+      DataField3D bad(p, "F", OOBPolicy::Abort, node_comm);   // aborts here
+      return 0;   // NOT reached
+   }
```
```diff
 test-data-field-3d-shared-mem: seas_test_data_field_3d_shared_mem
 	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_data_field_3d_shared_mem
 	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_data_field_3d_shared_mem
+	@echo "--- expecting ABORT (NaN fixture, shared mode) ---"
+	! timeout 60 $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 \
+	    ./seas_test_data_field_3d_shared_mem --nan-fixture-abort-child \
+	    > /dev/null 2>&1
```
(`! timeout 60 …` asserts nonzero exit AND bounds a hang at 60 s. On macOS
without `timeout`, use `gtimeout` from coreutils or guard with
`command -v timeout`; fall back to the un-timed `!` form.)

**Test case:** the Makefile stanza above IS the test (asserts abort, bounds
hang).

---

### [R-003] [LOW] [io/data_field_3d.cpp:load_] — Missing defensive checks: loading rank's `base` null, and `disp_unit` from shared_query unchecked

**Category:** ASSUMPTION

**Description:**
After `MPI_Win_allocate_shared` returns `MPI_SUCCESS` on the loading rank,
`base` is assumed non-null (standard-guaranteed for size>0, but a null deref
here would be a silent SEGV inside the converted writes). On peers, the
returned `qdisp` (disp_unit) is ignored; a mismatch would indicate a
heterogeneous build.

**Trigger:**
Buggy/exotic MPI implementation; mixed-ABI node. Not reachable on OpenMPI 4.x.

**Actual behavior:** no check.

**Expected behavior:** abort with a precise message.

**Suggested fix:**
```diff
       if (err != MPI_SUCCESS)
       {
          ...
       }
+      if (loads && base == nullptr)
+      {
+         H5Dclose(did);
+         H5Fclose(file);
+         MFEM_ABORT("DataField3D: MPI_Win_allocate_shared returned a null "
+                    "base pointer for field '" << field_name << "'");
+      }
       if (!loads)
       {
          MPI_Aint qsize = 0;
          int qdisp = 0;
          err = MPI_Win_shared_query(shared_win_, 0, &qsize, &qdisp, &base);
          ...
+         if (qdisp != static_cast<int>(sizeof(real_t)))
+         {
+            H5Dclose(did);
+            H5Fclose(file);
+            MFEM_ABORT("DataField3D: shared window disp_unit " << qdisp
+                       << " != sizeof(real_t) " << sizeof(real_t)
+                       << " for field '" << field_name << "'");
+         }
```

**Test case:** not demonstrable on a conforming MPI (defensive only) — hence LOW.

---

### [R-004] [LOW] [tests/unit/test_data_field_3d_shared_mem.cpp:27-42] — `std::array` used without `#include <array>`

**Category:** QUALITY

**Description:**
`sample_points()` returns `std::vector<std::array<real_t, 3>>` but the file
never includes `<array>`; it compiles via transitive inclusion from
`data_field_3d.hpp` (which includes `<array>`), so a future header cleanup
there would break this test.

**Trigger:** header hygiene change in an included header.

**Actual behavior:** compiles by luck of transitivity.

**Expected behavior:** self-sufficient includes.

**Suggested fix:**
```diff
 #include <cmath>
+#include <array>
 #include <cstdio>
```

**Test case:** n/a (compile-time; LOW).

---

## Fix round 1 (2026-07-17, /code-fix) — ALL FINDINGS RESOLVED
- **R-001 FIXED**: permanent np=1 regression
  `tests/unit/test_data_field_3d_shared_mem_finalize.cpp` added + wired into
  `test-data-field-3d-shared-mem`; plan §7 row reworded to the real mechanism
  (Finalized-guard skip; R-008 ordering forbids early reset).  Run: prints
  `dtor-after-finalize OK`, exit 0.
- **R-002 FIXED**: `--nan-fixture-abort-child` mode added to the shared-mem
  test (clean-exit-0 on a BROKEN guard so the negated Makefile check trips);
  Makefile stanza asserts nonzero mpirun exit AND distinguishes a hang
  (timeout/gtimeout exit 124 → FAIL) from the expected abort.  Run:
  `OK: aborted as expected (mpirun exit 1)`.
- **R-003 FIXED**: `loads && base == nullptr` abort after
  `MPI_Win_allocate_shared`; `qdisp != sizeof(real_t)` abort after
  `MPI_Win_shared_query` (both with field-named messages).
- **R-004 FIXED**: `#include <array>` added.
- Post-fix verification: shared-mem target green (np=2 12/12, np=4 24/24,
  finalize OK, NaN-abort OK); full classic suite unchanged (data-projection
  ALL PASSED, prestress 12/12, friction-sidecar 40/40, velocity-bundle 6/6,
  safs-params 25/25); driver relinks.

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 deviation/coverage, R-002 failure-path coverage)
- Low issues: 2 (R-003 defensive checks, R-004 include hygiene)
- Plan compliance: FULL for the implemented scope (Phases 0–2 + §6 gates;
  Phase 3 explicitly deferred to Expanse per the plan's own gating — plan
  status line updated accordingly). One §7 mitigation was reworded in code
  (Finalized-guard instead of impossible destroy-before-finalize) — flagged as
  R-001 rather than silently accepted.
- Verdict: **PASS WITH FIXES** — apply R-001/R-002 (coverage + doc wording)
  and the two LOW patches; no algorithmic changes required.

## Positive verification performed (not findings)
- Flag OFF ⇒ classic path: all five pre-existing sidecar test targets pass
  unchanged (data-projection roll-up, project-fault-prestress 12/12,
  spatial-friction-sidecar 40/40, spatial-velocity-bundle 6/6,
  compute-safs-params 25/25).
- Shared ⇒ per-rank bit-identity + cross-rank payload identity: np=2 12/12,
  np=4 24/24 (the corner sweep samples EVERY stored cell, so the §6(b)
  "checksum" requirement is fully covered, not sampled).
- Collective-ordering audit: all three loaders are config-driven and execute
  in identical order on every rank (velocity → stress → friction); node-comm
  collectives complete before the next world collective; no interleaving
  deadlock found.
- Copy/move safety: `DataField3D` copy deleted; no in-repo code copies or
  moves it (StressField3D members constructed in place; bundles hold
  `unique_ptr`); every consumer (driver, 2 standalone projector tools, all
  test binaries) recompiles and links.
- dtor-after-finalize probe: exit 0 (see R-001).

## Unreviewed Areas
- On-cluster behavior (Lustre HDF5 metadata storm at 128 ranks/node, OpenMPI
  4.1.x on Expanse, real ~253 MB windows) — Phase 3 scope, needs the cluster.
- Serial (non-MPI) MFEM build of the io/ files: code paths are `#ifdef`-clean
  by inspection, but no serial build exists in this environment to compile.
