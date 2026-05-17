# Code Review: PETSc TS restart plan (2026-05-16)

**Date:** 2026-05-16

**Plan reviewed:** `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`

**Scope:** Adversarial audit of the plan document itself (NOT yet-implemented code). Goal is to find specification bugs that would cause the next agent to implement broken code if they followed the plan literally. Every finding cross-references the actual codebase at HEAD (commit `407f456`) to verify the plan's claims about line numbers, ordering, and existing-code interactions.

## Review Scope
- Plan: `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md` (566 lines)
- Cross-referenced source files:
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp` (lines 1091-1099 early-abort to delete, 1572-1591 R-003 padded-shrink, 2222-2358 PETSc init + V1 restart, 2400-2426 post-Run, 2585-2638 checkpoint write sites).
  - `miniapps/seas/io/checkpoint.hpp` (V1 format, lines 56-130 Write, 151-235 Read).
  - `miniapps/seas/io/paraview_output.hpp` (lines 366-368 `GetTotalSnapshotsWritten`, 454 `last_write_time_(-1e30)`, 1501-1516 ShouldWrite cadence state, 1747-1752 `CommitOnceAtCycle`).
  - `linalg/petsc.cpp` (PetscODESolver::Run line 4357, Step line 4319).
- Domain context: `miniapps/seas/CLAUDE.md`, existing R-001..R-006 / R-104..R-106 reviews, the round-3 FIX_REPORT.

## Findings

### [R-301] [CRITICAL] [plan §"Files to Modify" + §"Format version detection"] — Plan contradicts itself about the checkpoint header tag

**Category:** DEVIATION (internal inconsistency)

**Description:**
The plan gives **two mutually incompatible specifications** for how the file-format version is encoded:

1. Line 69 ("Files to Modify → checkpoint.hpp") says:
   > "bump the file-format tag from `SEAS_CHECKPOINT_V1` to `SEAS_CHECKPOINT_V2`. The reader recognises both tags: V1 is read as today; V2 has an additional trailing block."

2. Line 154 ("Format version detection") says:
   > "ReadCheckpoint does not need to change its tag check — the V1 path STILL reads `SEAS_CHECKPOINT_V1` as today. The header tag itself remains `SEAS_CHECKPOINT_V1` for backwards compatibility."

These are mutually exclusive choices. The implementer would have to guess which one to apply. Each leads to different code:

- **Option A (bump header):** `WriteCheckpoint` writes `SEAS_CHECKPOINT_V2` at the top; `ReadCheckpoint` accepts EITHER `V1` or `V2`. V1 files are still readable; V2 files use the new tag. Forward-incompatible (old reader rejects V2 file).
- **Option B (keep header, append trailing block):** `WriteCheckpoint` writes `SEAS_CHECKPOINT_V1` unchanged; the new `WritePetscTSCheckpoint` appends a `PETSC_TS_V2` trailing block in append mode. Old reader stops at the V1 block and ignores the trailing extension. Fully backwards compatible.

The plan's actual code listings at lines 116-125 (writer body) and 150 (reader hint) clearly use **Option B** (no header tag change; new tag is the trailing `PETSC_TS_V2`). So line 69 is a documentation error — it should NOT say "bump the file-format tag from V1 to V2".

**Trigger:**
Implementer reads "Files to Modify" first (line 69), modifies the V1 header tag, then later reaches the "Format version detection" section (line 154) and is confused. May ship Option A (broken backwards compat) instead of the intended Option B.

**Actual behavior (if Option A is mistakenly chosen):**
- Old `seas_bp5_full` binaries can no longer read checkpoints written by the new binary (because they see `V2` and bail out at the existing `read_tag("SEAS_CHECKPOINT_V1")` MFEM_VERIFY).
- Operators with mixed-version environments (e.g., re-running an older job on Frontera while the build cluster has a newer binary) hit cryptic checkpoint-parse errors.

**Expected behavior:**
Option B is consistent with the plan's stated "Existing V1 checkpoint files (written by current production runs) must still be loadable by the new reader (forward-compatible)" constraint at line 29 AND with the writer body at line 119 (`out << "PETSC_TS_V2\n";` is OUTSIDE the V1 header section). Lock in Option B everywhere.

**Suggested fix:**

```diff
 ### Files to Modify

-- `miniapps/seas/io/checkpoint.hpp` — bump the file-format tag from `SEAS_CHECKPOINT_V1` to `SEAS_CHECKPOINT_V2`. The reader recognises both tags: V1 is read as today (no TS-internal fields); V2 has an additional trailing block read by `ReadPetscTSCheckpoint`. The writer always writes V2.
+- `miniapps/seas/io/checkpoint.hpp` — **NO CHANGES** to the file-format tag or to the existing `WriteCheckpoint`/`ReadCheckpoint` signatures.  The V2 extension is appended AFTER the V1 block by the new `WritePetscTSCheckpoint` function (defined in `petsc_ts_checkpoint.hpp`); the new `ReadPetscTSCheckpoint` reads the trailing `PETSC_TS_V2` block.  Old V1-only readers stop at the V1 block and silently ignore the V2 trailing block — full backwards compatibility.
```

**Test case:**
```cpp
// In test_bp5_petsc_ts_restart, add as sub-test 1.5:
void test_R301_v1_header_unchanged()
{
   // Write a V2 checkpoint (V1 + V2 blocks).  Read the file header
   // with the existing read_tag — must succeed at "SEAS_CHECKPOINT_V1".
   WriteCheckpoint(prefix, ...);
   WritePetscTSCheckpoint(prefix, ...);
   std::ifstream in(CheckpointFilename(prefix, 0));
   std::string tag;
   in >> tag;
   TEST_ASSERT(tag == "SEAS_CHECKPOINT_V1",
               "R-301: header tag must NOT bump to V2; backwards "
               "compatibility requires V1 header + appended V2 block");
}
```

---

### [R-302] [CRITICAL] [plan §"Files to Modify" item 2 + §"Detailed Requirements" item 4] — Phase 1 restart block placed at WRONG line in the driver; cross-check `std::abs(t - ts_t) < ...` compares uninitialised `t` against `ts_t`

**Category:** BUG (ordering)

**Description:**
The plan instructs the implementer to insert the V2 restart block at line 73:

> "Add in the PetscTS init block (after `petsc_ode->Init(...)` around line 2269): if `!restart_prefix.empty()`, call `ReadPetscTSCheckpoint(...)`."

And the code listing at line 169-209 has this cross-check at line 185:

```cpp
MFEM_VERIFY(std::abs(t - ts_t) < 1e-12 * std::abs(t),
            "Checkpoint inconsistency: V1 time=" << t
            << " differs from V2 time=" << ts_t);
```

This compares `t` (the V1 time) to `ts_t` (the V2 time). But `t` is only set by the existing `ReadCheckpoint` call which is at **line 2352** in the driver (not line 2269). At line ~2269 (where the plan wants the V2 block), `t` is still the default `t = 0.0` from line 2328.

Verified by grep:
```
$ grep -n "ReadCheckpoint" miniapps/seas/tests/verification/bp5_verification_full.cpp
2352:      bool ok = ReadCheckpoint(restart_prefix, t, restart_dt, ...
```

So if the implementer follows the plan literally, the V2 block runs BEFORE V1 has loaded `t`. The cross-check `std::abs(0 - ts_t) < 1e-12 * std::abs(0) = 0` always fails (since ts_t > 0 for any non-zero checkpoint time and 1e-12 * 0 = 0). MFEM_VERIFY aborts on every restart attempt with a misleading "Checkpoint inconsistency: V1 time=0 differs from V2 time=..." error.

**Trigger:**
Any `seas_bp5_full --restart PREFIX --petsc-ts` invocation after the plan is implemented as literally specified.

**Actual behavior:**
Driver aborts with "Checkpoint inconsistency: V1 time=0 differs from V2 time=..." on every restart.

**Expected behavior:**
The V2 block must run AFTER the existing V1 restart block at lines 2346-2358, so `t`, `current_dt`, `state`, etc. are populated before the cross-check.

**Suggested fix:**
Replace the "after `petsc_ode->Init(...)` around line 2269" placement with "after the existing V1 restart block at line ~2358, immediately before the time-stepping loop at line ~2390":

```diff
 ### Files to Modify

 - `miniapps/seas/tests/verification/bp5_verification_full.cpp`:
   - **Delete** the `--restart` + `--petsc-ts` early-abort at lines 1091-1099.
-  - **Add** in the PetscTS init block (after `petsc_ode->Init(...)` around line 2269): if `!restart_prefix.empty()`, call `ReadPetscTSCheckpoint(restart_prefix, t_from_ckpt, dt_next_from_ckpt, step_from_ckpt, rejections_from_ckpt, snapshots_so_far_from_ckpt, &mpi)`. Then `TSSetTime(ts, t_from_ckpt)`, `TSSetTimeStep(ts, dt_next_from_ckpt)`, `TSSetStepNumber(ts, step_from_ckpt)`. Set the monitor context's `num_seismic_events` / `in_seismic_event` from the V1 fields (already done by current restart path at lines 2346-2358).
+  - **Add** the V2 restart block IMMEDIATELY AFTER the existing V1 restart block (which currently lives at lines 2346-2358 and is conditioned on `!restart_prefix.empty()`).  The V2 block runs ONLY when `use_petsc_ts && !restart_prefix.empty()` and calls `ReadPetscTSCheckpoint`, then `TSSetTime(ts, ts_t)`, `TSSetTimeStep(ts, ts_dt_next)`, `TSSetStepNumber(ts, ts_step)`.  **Ordering is load-bearing**: by this point `t`, `current_dt`, `state`, `num_seismic_events`, `in_seismic_event` have ALL been populated by the V1 `ReadCheckpoint` at line 2352; the cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` is valid.  Do NOT place the V2 block inside the PetscTS init at line 2269 — `t` is still 0 there and the cross-check would always fail.
```

And update the code listing at lines 169-209 to reflect this — the comment "AFTER `petsc_ode->Init(...)` but BEFORE `petsc_ode->Run(...)`" should change to "AFTER the existing V1 restart block at lines 2346-2358".

**Test case:**
```cpp
// In test_bp5_petsc_ts_restart, sub-test 3 — the cross-check must pass:
void test_R302_v1_v2_cross_check_passes()
{
   // Run scenario A, checkpoint, restart, run scenario B.
   // The MFEM_VERIFY in the V2 restart block must NOT fire.
   // If R-302 isn't fixed, this test aborts with "Checkpoint
   // inconsistency: V1 time=0 differs from V2 time=...".
}
```

---

### [R-303] [MODERATE] [plan §"Detailed Requirements" item 4, line 199] — `step_rejections` assignment in the V2 restart block is immediately overwritten by existing post-Run code; previous-run rejections are silently lost

**Category:** BUG (ordering / dead store)

**Description:**
The plan's V2 restart block at line 199 assigns:
```cpp
step_rejections = ts_rejections;       // for the summary line
```

But the existing post-Run code at `bp5_verification_full.cpp:2414-2422` overwrites `step_rejections`:
```cpp
{
   petsc::TS ts = *petsc_ode;
   PetscInt ts_steps = 0;
   TSGetStepNumber(ts, &ts_steps);
   step = static_cast<int>(ts_steps);
   PetscInt rejects = 0;
   TSGetStepRejections(ts, &rejects);
   step_rejections = static_cast<int>(rejects);  // <-- OVERWRITES line-199's value
}
```

PETSc's `TSGetStepRejections` returns the rejection count tracked by the TS instance. `TSSetStepNumber` (which the V2 block calls per plan) does NOT reset the rejection counter — verified via PETSc 3.15 source `src/ts/interface/ts.c:TSSetStepNumber` only sets `ts->steps`. So `ts->reject` is whatever it was at PetscODESolver construction (0). After Run, `TSGetStepRejections` returns the rejections incurred during THIS Run only — NOT the cumulative count across the V1+V2 windows.

Net result: the final summary line at `bp5_verification_full.cpp:2668` (`step_rejections`) reports ONLY the new run's rejections, losing the saved `ts_rejections` from the checkpoint. The plan's "for the summary line" comment is misleading because the assignment is a dead store.

**Trigger:**
Any restart run. The summary line at end of run reports a too-low rejection count.

**Actual behavior:**
Summary reports rejections from the post-restart run only; pre-checkpoint rejections are lost.

**Expected behavior:**
Summary reports the cumulative rejection count across the whole simulation, OR explicitly labels the count as "post-restart only" with the previous count documented separately.

**Suggested fix:**
Save the pre-restart count in a separate variable that the summary code accumulates:

```diff
       step_rejections = ts_rejections;       // for the summary line
+      // R-303: this is the PRE-restart rejection count; the post-Run
+      // code at line 2414-2422 will overwrite `step_rejections` with
+      // the post-restart count.  Save the pre-restart count in a
+      // separate accumulator so the final summary reports the sum.
+      int restart_rejections_carryover = ts_rejections;
```

And in the post-Run code (lines 2414-2422):

```diff
       PetscInt rejects = 0;
       TSGetStepRejections(ts, &rejects);
-      step_rejections = static_cast<int>(rejects);
+      // R-303: TSGetStepRejections returns THIS-Run's rejections only;
+      // accumulate with the pre-restart count for the summary line.
+      step_rejections = restart_rejections_carryover
+                      + static_cast<int>(rejects);
```

The plan should add this requirement and update the listing at line 169-209.

**Test case:**
```cpp
void test_R303_rejections_accumulate_across_restart()
{
   // Run scenario A for 200 steps; record final step_rejections (call it R_A).
   // Run scenario B (100 steps + checkpoint + restart + 100 steps); record
   // final step_rejections (call it R_B).  Assert R_A == R_B.
   // Pre-fix: R_B equals only the post-restart rejection count, less than R_A.
   // Post-fix: R_B equals the cumulative count, matching R_A.
}
```

---

### [R-304] [MODERATE] [plan §"Numerical/behavioural constraints" + §"Risk Assessment"] — V2 checkpoint schema is missing `last_write_time_`, `last_v_max_`, and `current_regime_`; first ParaView ShouldWrite after restart will mis-fire

**Category:** EDGE_CASE (state not carried across restart)

**Description:**
The plan saves `total_snapshots_written_` (renamed `paraview_snapshots`) in the V2 block via `SetTotalSnapshotsWritten` (plan line 247-258). This handles the cap budget correctly. But the `ParaViewOutput::AdaptiveSchedule` state machine has THREE other per-instance fields that drive write decisions and are NOT in the V2 schema:

| Field                  | Default value     | What it controls                                                              |
|------------------------|-------------------|-------------------------------------------------------------------------------|
| `last_write_time_`     | `-1e30` (line 454)| `ShouldWrite` compares `time - last_write_time_ >= dt_out * tol`              |
| `last_v_max_`          | `0.0`             | Backs the single-arg `CommitSchedule(time)` shim                              |
| `current_regime_`      | `0`               | State-machine hysteresis (`NextRegime(V, prev_regime)`)                       |

After the V2 restart loads `total_snapshots_written_` but leaves `last_write_time_ = -1e30`, the first `pv_out->ShouldWrite(...)` call at the post-restart time `T` computes `T - (-1e30) = huge >> dt_out * 0.99` and WRITES UNCONDITIONALLY. This produces an extra snapshot that the fresh-run trajectory would not produce. With `max_total_snapshots = 5000` and 100s of restarts in a long campaign, the user accumulates an unbounded number of spurious snapshots.

Similarly, `current_regime_ = 0` on restart means hysteresis state is lost — if the pre-checkpoint regime was `2` (coseismic) with `hysteresis_factor > 1`, the V_max value that should keep the system in regime 2 (because it's above `V_co_exit = V_coseismic / hysteresis_factor`) instead drops to regime 0 (because we started fresh at regime 0). The first paraview_write after restart fires with the wrong regime classification.

The plan's risk table at line 522 mentions `last_committed_cycle_` (which actually doesn't need restoration — the dedup is cycle-relative, not absolute) but completely misses the three fields above that DO need restoration.

**Trigger:**
Any restart of a run that has emitted at least one paraview snapshot before checkpointing.

**Actual behavior:**
- First ShouldWrite after restart always returns true, producing an extra snapshot at `t = saved_t`.
- Hysteresis state lost; regime classification may transition incorrectly on the first post-restart step.

**Expected behavior:**
The V2 schema includes `last_write_time_`, `last_v_max_`, and `current_regime_`. `SetTotalSnapshotsWritten` is supplemented by a `RestoreScheduleState(last_write_time, last_v_max, current_regime)` setter that restores all three.

**Suggested fix:**
Extend the V2 schema and `WritePetscTSCheckpoint` / `ReadPetscTSCheckpoint` signatures:

```diff
 ///   PETSC_TS_V2
 ///   petsc_ts_time     <real_t>
 ///   petsc_ts_dt_next  <real_t>
 ///   petsc_ts_step     <int>
 ///   petsc_ts_rejections <int>
 ///   paraview_snapshots <int>
+///   paraview_last_write_time   <real_t>
+///   paraview_last_v_max        <real_t>
+///   paraview_current_regime    <int>
```

And in `ParaViewOutput`, alongside `SetTotalSnapshotsWritten`, add:

```cpp
/// @brief R-304: restore the adaptive-schedule state from a checkpoint.
///
/// Pairs with `SetTotalSnapshotsWritten`.  Use ONLY at restart time,
/// before any Save/ShouldWrite call.  Without this restoration the
/// first ShouldWrite after restart fires unconditionally (because
/// `last_write_time_` defaults to -1e30) and the regime state machine
/// resets to interseismic.
void RestoreScheduleState(real_t last_write_time, real_t last_v_max,
                          int current_regime)
{
   last_write_time_ = last_write_time;
   last_v_max_      = last_v_max;
   current_regime_  = current_regime;
}
```

And on the write side (plan line 222-238), add accessors:
```cpp
real_t last_write_time = pv_out ? pv_out->GetLastWriteTime() : -1e30;
real_t last_v_max      = pv_out ? pv_out->GetLastVMax()      : 0.0;
int    current_regime  = pv_out ? pv_out->GetCurrentRegime() : 0;
```

(The getters need to be added to `ParaViewOutput` as well — currently `last_write_time_` etc. are private.)

**Test case:**
```cpp
void test_R304_first_post_restart_should_write_does_not_fire()
{
   // Run A: ShouldWrite at t=0 (fires), t=0.5*dt_inter (does not fire),
   //        t=1.0*dt_inter (fires).  Checkpoint at t=1.0*dt_inter
   //        immediately after the fire.
   // Run B: restart from the checkpoint at t=1.0*dt_inter.
   //        Call ShouldWrite at t=1.001*dt_inter (immediately after
   //        the saved time, well within tol).
   // Pre-fix: ShouldWrite returns true (extra snapshot!).
   // Post-fix: ShouldWrite returns false (last_write_time_ properly
   //           restored, so time - last_write_time_ ≈ 0 < tol).
   TEST_ASSERT(!pv.ShouldWrite(step_after_restart, 1.001 * dt_inter, V),
               "R-304: first ShouldWrite after restart must NOT fire "
               "if t is within tol of the saved last_write_time");
}
```

---

### [R-305] [MODERATE] [plan §"Testing Strategy" sub-test 3] — Procedure says "Run PetscTS for K=200 steps" but `PetscODESolver::Run` doesn't take a step count; PETSc adaptive RK45 makes "200 steps" non-deterministic

**Category:** ASSUMPTION (test design is under-specified)

**Description:**
The plan's most important test (Sub-test 3, lines 410-435) is the fresh-vs-restart trajectory agreement test. The procedure says:

> "2. Run scenario A (fresh): create state, SetInitialCondition, run PetscTS for K=200 steps. Save final state_A, final t_A, final dt_A."
> "3. Run scenario B (restart-mid-way): a. ... b. Run PetscTS for K/2=100 steps."

But the `PetscODESolver::Run(Vector &x, real_t &t, real_t &dt, real_t t_final)` API (linalg/petsc.cpp:4357) takes a `t_final`, not a step count. Adaptive RK45 with `atol = 1e-7` produces non-deterministic step counts because the controller's dt depends on local error magnitudes; two seemingly-identical runs can take a slightly different number of steps to reach the same `t_final`. So "run for 200 steps" is ambiguous — it could mean:

a. Set `t_final = t_init + 200 * dt_init` and call `Run` once. Adaptive controller may take more or fewer than 200 steps.
b. Call `PetscODESolver::Step(state, t, dt)` (linalg/petsc.cpp:4319) 200 times in a loop.
c. Set `MaxIter = 200` on the TS and call Run with `t_final = LARGE`. Run exits when 200 steps complete.

The implementer would have to pick one without guidance. Option (a) breaks the "fresh vs restart same step count" assumption. Option (b) requires loop bookkeeping. Option (c) requires `petsc_ode->SetMaxIter(200)` and may interact with the existing `SetMaxIter(max_steps)` at `bp5_verification_full.cpp:2268`.

The "Save final t_A" + "Save final t_B" comparison at step 4 (`|t_A - t_B| < 1e-12`) assumes both runs end at the same `t`. Under option (a), they DO (both stop at `t_final`). Under option (b)/(c), they MAY differ because the adaptive controller picks slightly different dt sequences in the two runs even before the restart.

**Trigger:**
Implementer reads the test description, picks one of the three interpretations, ships the test. If they pick (b) or (c) under default PETSc behaviour, the `t_A == t_B` assertion may pass OR fail depending on whether the trajectories happen to align.

**Actual behavior:**
Test is ambiguous; under (b)/(c) the equality `t_A == t_B` may not hold even when the trajectory IS correct.

**Expected behavior:**
The plan specifies exactly which API to use, and the assertions match that API's guarantees.

**Suggested fix:**
Use approach (c) with `TSSetMaxSteps` to fix the step count exactly, and call `Run` with `t_final = HUGE_VAL`. Then `TSGetStepNumber` after Run reports exactly N steps (assuming no NaN). Then `t_A == t_B` holds because both runs took exactly N steps from the same initial condition along the same trajectory.

Actually the most robust approach is **option (a) but state in terms of t_final, not step count**:

```diff
 **Sub-test 3: Fresh-vs-restart trajectory agreement (the core correctness test).**

 This is the most important test. Procedure (in pseudocode):

 ```
-1. Construct a minimal BP5-shaped problem (inline mesh, P=1 IP, ~1000
-   fault DOFs, atol=1e-7, rtol=1e-50).  Construct seas_op, fault_op.
-2. Run scenario A (fresh): create state, SetInitialCondition, run
-   PetscTS for K=200 steps.  Save final state_A, final t_A, final dt_A.
-3. Run scenario B (restart-mid-way):
-     a. Create state, SetInitialCondition (identical to A).
-     b. Run PetscTS for K/2=100 steps.
-     c. Call WriteCheckpoint + WritePetscTSCheckpoint to /tmp/test_restart.
+1. Construct a minimal BP5-shaped problem (inline mesh, P=1 IP, ~1000
+   fault DOFs, atol=1e-7, rtol=1e-50).  Construct seas_op, fault_op.
+   Pick T_full and T_mid such that T_mid < T_full.  Suggested:
+   T_full = 100 * dt_init, T_mid = 50 * dt_init.
+2. Run scenario A (fresh): create state, SetInitialCondition, call
+   petsc_ode->Run(state, t, dt, T_full).  Save final state_A, t_A=T_full,
+   final dt_A.  PetscTS uses TSSolve with TS_EXACTFINALTIME_MATCHSTEP
+   (set by MFEM at petsc.cpp:4365), so t_A == T_full exactly.
+3. Run scenario B (restart-mid-way):
+     a. Create state, SetInitialCondition (identical to A).
+     b. Call petsc_ode->Run(state, t, dt, T_mid).
+     c. Call WriteCheckpoint + WritePetscTSCheckpoint.
      d. Destroy everything; re-construct fresh.
      e. ReadCheckpoint + ReadPetscTSCheckpoint, restore t, dt, state, TS.
-     g. Run PetscTS for the remaining K/2=100 steps.
+     g. Call petsc_ode->Run(state, t, dt, T_full).
      h. Save final state_B, final t_B, final dt_B.
 4. Assert:
-   - |t_A - t_B| < 1e-12  (time should match exactly — both runs
-     terminate at the same K-th step's accepted t).
+   - t_A == t_B == T_full bit-exactly (TS_EXACTFINALTIME_MATCHSTEP).
    - max_i |state_A(i) - state_B(i)| <= atol + rtol * max_i |state_A(i)|.
 ```
```

**Test case:**
The corrected test description above IS the test specification.

---

### [R-306] [LOW] [plan §"Files to Modify" item 4 vs. §"Detailed Requirements" item 6] — `SetTotalSnapshotsWritten` body is specified twice with different implementations

**Category:** QUALITY (specification ambiguity)

**Description:**
Two different bodies for `SetTotalSnapshotsWritten`:

- Line 78 (summary in "Files to Modify"):
  ```cpp
  void SetTotalSnapshotsWritten(int n) { total_snapshots_written_ = n; }
  ```
  No input validation.

- Line 254-257 (detailed requirement):
  ```cpp
  void SetTotalSnapshotsWritten(int n)
  {
     total_snapshots_written_ = std::max(n, 0);
  }
  ```
  Clamps negative values to 0.

The detailed version is clearly the intended implementation (the doc-comment at line 248-253 mentions the clamp), but a quick-reading implementer might cut-and-paste the summary version and skip the clamp.

**Trigger:**
Restart from a corrupted checkpoint that contains a negative `paraview_snapshots` value (file truncation, partial write, manual edit). The non-clamping version silently sets `total_snapshots_written_ = -5` (say); subsequent `CommitOnceAtCycle` increments it to -4, -3, ..., 0, 1, ..., effectively giving the user `max_total_snapshots + |negative|` extra writes before the cap engages.

**Actual behavior:**
With the non-clamping version, a corrupted checkpoint silently extends the cap budget.

**Expected behavior:**
Both summary and detailed sections show the clamping version.

**Suggested fix:**
```diff
 - `miniapps/seas/io/paraview_output.hpp`:
-  - **Add** a public setter `void SetTotalSnapshotsWritten(int n) { total_snapshots_written_ = n; }`. Doc-comment: "Used by `--restart` to restore the snapshot counter; should only be called once, immediately after construction, before any `Save`/`ShouldWrite` call. Negative values are clamped to 0."
+  - **Add** a public setter `SetTotalSnapshotsWritten(int n)` whose body clamps negative inputs: `total_snapshots_written_ = std::max(n, 0);`.  Full body and doc-comment are in §"Detailed Requirements" item 6 below.  Negative values are clamped because a corrupted V2 block (e.g., partial write) could otherwise quietly extend the cap budget by abs(n).
```

**Test case:**
```cpp
void test_R306_setter_clamps_negative()
{
   ParaViewOutput<Mesh> pv("/tmp/test_R306", smesh, 1);
   pv.SetTotalSnapshotsWritten(-5);
   TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 0,
               "R-306: SetTotalSnapshotsWritten must clamp negative "
               "inputs to 0 to defend against corrupted checkpoints");
}
```

---

### [R-307] [LOW] [plan §"Phase 3" line 367] — `TSRKSetStageVectors` is referenced but does not exist in upstream PETSc; only the proposed `TSRKGetStageVectors` is mentioned

**Category:** ASSUMPTION (Phase-3 dependency unspecified)

**Description:**
Phase 3 line 359 proposes adding `TSRKGetStageVectors` to upstream PETSc:

> "**NEW upstream PETSc PR** required: add `PetscErrorCode TSRKGetStageVectors(TS ts, PetscInt *nstages, Vec **Y)` to `src/ts/impls/explicit/rk/rk.c`. This is the only public API gap."

But line 367 then says:

> "On the READ side: load the values, allocate a PETSc Vec of the right size, and set `Y[6]` (or equivalent) on the active TS via `TSRKSetStageVectors`."

`TSRKSetStageVectors` is a DIFFERENT function (Set, not Get) that ALSO doesn't exist in upstream PETSc and is NOT in the proposed PR list. The plan needs both Get (to extract on checkpoint) AND Set (to restore on restart). Currently the plan only flags Get as a PR.

**Trigger:**
Implementer of Phase 3 reads line 367, tries to call `TSRKSetStageVectors`, gets a linker error. Goes back to the plan and finds no spec for Set. Has to invent the API.

**Actual behavior (if implementer just adds the Set function silently):**
Likely correct but un-reviewed PETSc patch shipped alongside the Get patch.

**Expected behavior:**
Phase 3 explicitly lists BOTH the Get and Set additions as upstream PETSc PRs (or both as MFEM-side private-API workarounds).

**Suggested fix:**
```diff
 ### Files to Modify

 - `miniapps/seas/io/petsc_ts_checkpoint.hpp` — extend V2 block (or bump to V3) with `petsc_ts_fsal_stage_size N` + N values for the DP5 FSAL stage k_7 = k_1 vector.

-- **NEW upstream PETSc PR** required: add `PetscErrorCode TSRKGetStageVectors(TS ts, PetscInt *nstages, Vec **Y)` to `src/ts/impls/explicit/rk/rk.c`. This is the only public API gap.
+- **NEW upstream PETSc PR** required: add BOTH of these to `src/ts/impls/explicit/rk/rk.c`:
+  * `PetscErrorCode TSRKGetStageVectors(TS ts, PetscInt *nstages, Vec **Y)` — for the checkpoint WRITE side.
+  * `PetscErrorCode TSRKSetStageVectors(TS ts, PetscInt nstages, Vec *Y)` — for the checkpoint READ side.  Restores the FSAL stage so the first step after restart uses the inherited K[6] instead of re-evaluating it.
+  Both functions are new public API.  Without `TSRKSetStageVectors`, the Phase 3 read side cannot complete and Phase 3 collapses to (essentially) Phase 1 behaviour.
```

No test case needed — this is a Phase-3 spec issue and Phase 3 is blocked-on-upstream regardless.

---

## Summary
- Critical issues: 2 (R-301, R-302)
- Moderate issues: 3 (R-303, R-304, R-305)
- Low issues: 2 (R-306, R-307)
- Plan compliance: **PARTIAL** — the plan is well-structured and covers the right concepts (V1/V2 layout, R-003 interaction, cap-counter preservation), but contains two CRITICAL self-contradictions that would cause the next agent to ship broken code if followed literally, plus three MODERATE specification gaps that would yield subtly-wrong behaviour or under-specified tests.
- Verdict: **FAIL — must fix R-301, R-302, R-303, R-304, R-305 before the next agent implements.** R-301 (header tag contradiction) and R-302 (V2 block placement) are the highest-priority because they would cause the implementation to fail at restart time and at compile/runtime respectively. R-303 (rejection accumulator) and R-304 (missing schedule state) cause silently-wrong behaviour in production. R-305 (under-specified test) reduces the value of the regression test. R-306 (clamping ambiguity) and R-307 (Phase 3 Set API) can be fixed when the relevant section is touched.

## Unreviewed Areas
- The plan's Phase 2 claim that TSADAPTBASIC has no multi-step history was NOT verified against PETSc 3.15 source. The plan acknowledges this gap in the Risk Assessment table (line 517: "Read PETSc 3.15 source ... to verify"). If the claim is wrong, Phase 2 has hidden requirements.
- The Risk Assessment item at line 519 (`TSSetStepNumber` not resetting `_TSMonitorSet` log) was not verified against PETSc source. If the claim is wrong, monitor `step` arguments after restart would not be correctly contiguous with pre-checkpoint `step`s, potentially confusing the I/O cadence.
- The plan's interaction between V2 restart and the existing checkpoint-restart code at `bp5_verification_full.cpp:2346-2358` (V1 path) was reviewed only for ordering (R-302); the plan does NOT specify what happens if `restart_prefix` is set but `use_petsc_ts` is false on a V2-only file. The current code would happily read the V1 fields and ignore the V2 trailer — which IS the right behaviour, but the plan should make this explicit (probably an "Edge Cases" bullet alongside the existing "V1 file fed to a V2-aware run" bullet at line 275).
