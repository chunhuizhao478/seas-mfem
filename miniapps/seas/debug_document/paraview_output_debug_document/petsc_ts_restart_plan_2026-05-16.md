# Implementation Plan — PETSc TS state serialization for `--restart` support in `seas_bp5_full`

**Date:** 2026-05-16
**Branch:** `feature/paraview-compaction` (or a new branch off it)
**Target file (per FIX_REPORT note "Option C"):** make `--restart` + `--petsc-ts` work instead of aborting.

## Overview

Today the BP5 driver rejects the combination `--restart PREFIX` + `--petsc-ts` with a hard error (`bp5_verification_full.cpp:1091-1099`) because the existing per-rank checkpoint only saves MFEM-side state (slip+psi vector, displacement, traction, slip_rate, FSAL k0). PETSc TS internal state — TS time/dt/step number, TSAdapt history, and TSRK stage memory — is NOT in the checkpoint, so a naive restart would silently lose adaptive-dt history and (for DP5) the FSAL stage.

This plan delivers a `--restart` + `--petsc-ts` combination that works correctly at **RK45 tolerance** (i.e., the restarted trajectory matches a fresh-run trajectory to within `atol + rtol * |y|`), while explicitly NOT promising bit-exact reproducibility. Bit-exact restart requires PETSc-internal access to TSAdapt and TSRK that PETSc does not expose in standard releases; the gap is documented and tracked as future work.

The plan is broken into three phases:

| Phase | Scope                                          | Deliverable                                    | MFEM/PETSc patches needed?     |
|-------|------------------------------------------------|------------------------------------------------|--------------------------------|
| 1     | Tolerance-correct restart via `TSSetTime` / `TSSetTimeStep` / pre-PlaceMemory state | Drop the early-abort; checkpoint format bumped to V2; tests verify trajectory agreement within atol | None — uses only public PETSc TS API |
| 2     | TSAdapt state extraction so the controller's "next dt" matches across restart | One extra `petsc_ts_dt_next` field in the checkpoint; restart calls `TSSetTimeStep(dt_next)` instead of guessing | None — uses public `TSGetTimeStep` (already saved); requires verification that the "next" dt is the one used |
| 3     | Bit-exact restart (FSAL stage preserved across the seam) | Saves the DP5 FSAL stage vector; restart skips the first-stage re-evaluation | Yes — MFEM patch to expose TSRK stage memory, OR upstream PETSc PR to add `TSRKGetStageVectors` |

Phase 1 is the production deliverable. Phases 2 and 3 are documented as follow-ups with clear interfaces so a future maintainer can pick them up without re-discovering the analysis.

## Constraints

**Interface constraints (must not change):**

- `WriteCheckpoint` and `ReadCheckpoint` signatures in `miniapps/seas/io/checkpoint.hpp` are stable; existing non-PetscTS callers must keep working. New TS-related fields go in a separate function (`WritePetscTSCheckpoint` / `ReadPetscTSCheckpoint`) called only when `use_petsc_ts == true`.
- The per-rank checkpoint filename pattern `{prefix}_checkpoint_r{rank}.txt` stays unchanged. The new TS-side fields are appended INSIDE that file (extending V1 → V2), not in a separate file, so that a single `--restart PREFIX` argument loads everything.
- Existing V1 checkpoint files (written by current production runs) must still be loadable by the new reader (forward-compatible). The new reader detects V1 vs V2 by the header tag.
- The driver-side `--restart` and `--checkpoint-interval` CLI flags keep their current names and semantics.
- The R-003 padded-shrink workaround at `bp5_verification_full.cpp:1572-1591` must continue to function — the restart path must not regress the zero-fault-DOF-rank invariant.

**Dependency constraints:**

- Phase 1 uses ONLY publicly-documented PETSc TS API: `TSSetTime`, `TSSetTimeStep`, `TSSetStepNumber`, `TSGetTime`, `TSGetTimeStep`, `TSGetStepNumber`, `TSGetStepRejections`. All are stable across PETSc 3.15+ (Frontera ships PETSc 3.15).
- `MFEM::PetscODESolver` already exposes `operator petsc::TS()` for raw TS access (`linalg/petsc.hpp:979`), used existing-style at `bp5_verification_full.cpp:2270-2291` and `2414-2422`. No new MFEM wrapper required for Phase 1.
- Plan-text checkpoint format (one value per line, 17-digit scientific) is preserved. Binary checkpoint format is out of scope.

**Convention constraints:**

- Naming follows existing seas pattern: snake_case functions, `seas::` namespace, `mfem::seas::` in headers.
- Floating-point precision: `std::setprecision(17)` + `std::scientific` (matches existing).
- MFEM_VERIFY for unrecoverable errors; `return false` for missing-file kinds of failures (matches `ReadCheckpoint`).
- Tests follow the `seas_test_*` pattern (e.g., `seas_test_bp5_petsc_ts_restart`), wired into the Makefile next to `seas_test_bp5_petsc_ts_zero_fault_rank` (the existing R-003 pattern).

**Numerical / behavioural constraints:**

- Phase 1: restarted trajectory must match a fresh-run trajectory at the same final time to within `atol + rtol * |y|` (BP5 atol = 1e-7, rtol = 1e-50 — so effectively `1e-7`). This is the standard "correct-within-tolerance" restart contract for adaptive ODE solvers.
- The cap-aware ParaView schedule (`SnapshotCapAwareInterval`) must continue to function across restart — the snapshot counter (`total_snapshots_written_`) must be saved/restored, otherwise restart resets the cap budget and a long-running job would over-write past K = 5000.
- The earthquake-detection state machine (`num_seismic_events`, `in_seismic_event`) is already saved/restored by V1; behaviour preserved.
- The R-003 padded-shrink at `bp5_verification_full.cpp:1572-1591` must run BEFORE the restart's `ReadCheckpoint` call (i.e., construct `state` with the padded-shrink trick, then `ReadCheckpoint` `SetSize(n_from_file)` — which is a shrink, capacity preserved). See R-105 doc-comment.

---

## Phase 1: Tolerance-correct PETSc TS restart

### Goal

After this phase, the combination `--restart PREFIX --petsc-ts` is accepted by `seas_bp5_full` and produces a trajectory whose long-term behaviour matches a fresh run to within RK45 atol+rtol tolerance. The early-abort at `bp5_verification_full.cpp:1091-1099` is gone.

### Files to Create

- `miniapps/seas/io/petsc_ts_checkpoint.hpp` — new header with `WritePetscTSCheckpoint(prefix, t, dt_next, step, rejections, snapshots_so_far, mpi)` and `ReadPetscTSCheckpoint(prefix, t, dt_next, step, rejections, snapshots_so_far, mpi)`. Format: appended to the existing per-rank file as a `PETSC_TS_V2` trailing block. Both functions are no-ops on non-MFEM_USE_PETSC builds (compile-time `#ifdef`).

- `miniapps/seas/tests/unit/test_bp5_petsc_ts_restart.cpp` — new MPI integration test. Drives a synthetic short BP5 run, checkpoints mid-run, restarts, and asserts trajectory agreement within atol+rtol. See Testing Strategy.

### Files to Modify

- `miniapps/seas/io/checkpoint.hpp` — **NO CHANGES** to the file-format tag or to the existing `WriteCheckpoint` / `ReadCheckpoint` signatures.  The V2 extension is appended AFTER the V1 block by the new `WritePetscTSCheckpoint` function (defined in `petsc_ts_checkpoint.hpp`); the new `ReadPetscTSCheckpoint` reads the trailing `PETSC_TS_V2` block.  Old V1-only readers stop at the V1 block and silently ignore the V2 trailing block — full backwards compatibility.  See §"Format version detection" below for the canonical statement.

- `miniapps/seas/tests/verification/bp5_verification_full.cpp`:
  - **Delete** the `--restart` + `--petsc-ts` early-abort at lines 1091-1099.
  - **Add** the V2 restart block IMMEDIATELY AFTER the existing V1 restart block at lines 2346-2358.  **Ordering is load-bearing**: the V1 `ReadCheckpoint` call at line 2352 populates `t`, `current_dt`, `state`, `num_seismic_events`, `in_seismic_event` from the V1 fields; the V2 block then reads `ts_t`, `ts_dt_next`, `ts_step`, `ts_rejections`, `pv_snapshots` from the trailing `PETSC_TS_V2` block and calls `TSSetTime(ts, ts_t)`, `TSSetTimeStep(ts, ts_dt_next)`, `TSSetStepNumber(ts, ts_step)`.  The cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` is only meaningful BECAUSE `t` has been populated by the preceding V1 read.  **Do NOT place the V2 block inside the PetscTS init at line ~2269** — `t` is still 0 there, the cross-check would always fail, and every restart attempt would abort with a misleading "V1 time=0 differs from V2 time=..." error.
  - **Modify** the checkpoint-writing call sites at lines 624-628 (in the monitor) and 2591 (mid-run) and 2635 (final) to ALSO call `WritePetscTSCheckpoint` immediately after `WriteCheckpoint`. The fields come from: `TSGetTime`, `TSGetTimeStep`, `TSGetStepNumber`, `TSGetStepRejections` on the `ts` handle; `snapshots_so_far = pv_out ? pv_out->GetTotalSnapshotsWritten() : 0`; per R-304, also `pv_out->GetLastWriteTime()`, `pv_out->GetLastVMax()`, `pv_out->GetCurrentRegime()` (new accessors — see item 6 below).
  - **Add** in the `paraview_write` lambda (or in the PetscTS init block): after `pv_out` is constructed, if `!restart_prefix.empty() && use_petsc_ts`, restore the snapshot counter via a new `pv_out->SetTotalSnapshotsWritten(snapshots_so_far_from_ckpt)` setter (added in next bullet).

- `miniapps/seas/io/paraview_output.hpp`:
  - **Add** a public setter `SetTotalSnapshotsWritten(int n)` whose body clamps negative inputs: `total_snapshots_written_ = std::max(n, 0);`.  Doc-comment: "Used by `--restart` to restore the snapshot counter; should only be called once, immediately after construction, before any `Save`/`ShouldWrite` call.  Negative values are clamped to 0 to defend against corrupted/truncated V2 checkpoints."  Full body in §"Detailed Requirements" item 6 below.
  - **Add** a public setter `RestoreScheduleState(real_t last_write_time, real_t last_v_max, int current_regime)` that restores the three private fields from a V2 checkpoint.  Required by R-304 to prevent the first post-restart ShouldWrite from firing unconditionally and to preserve regime hysteresis across the seam.  Full body in §"Detailed Requirements" item 6 below.
  - **Add** three public read-only accessors: `GetLastWriteTime()`, `GetLastVMax()`, `GetCurrentRegime()`.  Used by the WRITE side to snapshot the schedule state at checkpoint time.  Full bodies in §"Detailed Requirements" item 6 below.

- `miniapps/seas/Makefile` — add build rules for `seas_test_bp5_petsc_ts_restart`, following the pattern of `seas_test_bp5_petsc_ts_zero_fault_rank` at lines 232, 446, 1070-1072, 1660-1662. Add the new MPI-only target `test-bp5-petsc-ts-restart: seas_test_bp5_petsc_ts_restart` invoking `$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_bp5_petsc_ts_restart`. Do NOT add to the `test:` umbrella (mirrors the R-104 decision — MPI required).

### Detailed Requirements

**1. `WritePetscTSCheckpoint` signature and body.**

```cpp
namespace mfem { namespace seas {

#ifdef MFEM_USE_PETSC
/// @brief Append PETSc TS internal state to an already-written V1 checkpoint.
///
/// Must be called AFTER WriteCheckpoint on the same prefix.  Opens the
/// file in append mode and writes the PETSC_TS_V2 trailing block:
///   PETSC_TS_V2
///   petsc_ts_time              <real_t>
///   petsc_ts_dt_next           <real_t>
///   petsc_ts_step              <int>
///   petsc_ts_rejections        <int>
///   paraview_snapshots         <int>
///   paraview_last_write_time   <real_t>   (R-304)
///   paraview_last_v_max        <real_t>   (R-304)
///   paraview_current_regime    <int>      (R-304)
///
/// `petsc_ts_dt_next` is the dt PETSc will use for the NEXT step, as
/// returned by TSGetTimeStep AFTER the just-completed step.  This is
/// the value we need to feed back via TSSetTimeStep on restart.
///
/// `paraview_last_write_time`, `paraview_last_v_max`,
/// `paraview_current_regime` capture the ParaViewOutput adaptive
/// schedule state at checkpoint time so the restart's first
/// ShouldWrite call does not fire unconditionally (the default
/// `last_write_time_ = -1e30` would otherwise make any `time - last`
/// look huge and trigger an extra snapshot at t=saved_t).  See R-304.
inline void WritePetscTSCheckpoint(const std::string &prefix,
                                   real_t t, real_t dt_next,
                                   int step, int rejections,
                                   int paraview_snapshots,
                                   real_t paraview_last_write_time,
                                   real_t paraview_last_v_max,
                                   int paraview_current_regime,
                                   const MPIContext *mpi);
#endif

}}
```

Body: appends to `CheckpointFilename(prefix, rank)` (the same file as `WriteCheckpoint`).

```cpp
std::ofstream out(filename, std::ios::app);
out << std::setprecision(17) << std::scientific;
out << "PETSC_TS_V2\n";
out << "petsc_ts_time " << t << "\n";
out << "petsc_ts_dt_next " << dt_next << "\n";
out << "petsc_ts_step " << step << "\n";
out << "petsc_ts_rejections " << rejections << "\n";
out << "paraview_snapshots " << paraview_snapshots << "\n";
out << "paraview_last_write_time " << paraview_last_write_time << "\n";
out << "paraview_last_v_max " << paraview_last_v_max << "\n";
out << "paraview_current_regime " << paraview_current_regime << "\n";
```

**2. `ReadPetscTSCheckpoint` signature and body.**

```cpp
#ifdef MFEM_USE_PETSC
/// @brief Read the PETSC_TS_V2 trailing block from a checkpoint file.
///
/// Must be called AFTER ReadCheckpoint (which reads the V1 prefix).
/// Returns true if the trailing PETSC_TS_V2 block was found and parsed;
/// returns false (with sensible defaults written into the out-params)
/// if the file is V1-only (no trailing block).  Aborts on a malformed
/// V2 block.
///
/// Caller MUST treat a `false` return as a hard error when
/// `--petsc-ts` is active, because the V1-only fall-back ("restart
/// the MFEM-side state but let PETSc TS start fresh from dt_init")
/// produces a misleading trajectory.  The driver enforces this with
/// MFEM_VERIFY (see §"Driver-side restart-with-PetscTS flow" below).
///
/// @return true if PETSc TS state was loaded; false if file is V1
/// (no PetscTS fields).
inline bool ReadPetscTSCheckpoint(const std::string &prefix,
                                  real_t &t, real_t &dt_next,
                                  int &step, int &rejections,
                                  int &paraview_snapshots,
                                  real_t &paraview_last_write_time,
                                  real_t &paraview_last_v_max,
                                  int &paraview_current_regime,
                                  const MPIContext *mpi);
#endif
```

Body: reopens the file, skips past the V1 block (uses a temporary inner reader for the V1 fields, discarding values — OR, simpler, reads the whole file into a string and seeks for the `PETSC_TS_V2` tag). Recommended: just open the file, read tags one at a time until `PETSC_TS_V2` is found OR EOF is hit; on EOF return false; on found, parse the 5 fields.

**3. Format version detection.**

`ReadCheckpoint` does not need to change its tag check — the V1 path STILL reads `SEAS_CHECKPOINT_V1` as today (the V2 extension is appended after V1, not a replacement of the header). The header tag itself remains `SEAS_CHECKPOINT_V1` for backwards compatibility. The TS-side fields are gated by their own `PETSC_TS_V2` tag inside the same file.

This means a V2-aware writer produces a file that an old V1-only reader can still parse (it just stops at the FSAL fields and leaves the TS extension unread). And a new V2-aware reader handling an old V1 file finds no `PETSC_TS_V2` tag and returns false from `ReadPetscTSCheckpoint`, allowing the caller to fall back to "restart but PETSc TS starts fresh" (which is the pre-fix behaviour — abort, OR optionally proceed with dt = dt_init).

**4. Driver-side restart-with-PetscTS flow.**

Replace the early-abort block with:

```cpp
// bp5_verification_full.cpp ~ line 1091
// (DELETE the early-abort block.)
```

Add a new restart block IMMEDIATELY AFTER the V1 restart block (`if (!restart_prefix.empty()) { ... bool ok = ReadCheckpoint(...) ... }` at lines 2346-2358).  Both blocks are gated on `!restart_prefix.empty()`; the V2 block has the additional `&& use_petsc_ts` gate.  Placing the V2 block AFTER V1 is required (R-302): the cross-check at line 11 below depends on `t` being populated by the V1 `ReadCheckpoint`.

```cpp
#ifdef MFEM_USE_PETSC
if (use_petsc_ts && !restart_prefix.empty())
{
   real_t ts_t = 0.0, ts_dt_next = 0.0;
   int ts_step = 0, ts_rejections = 0, pv_snapshots = 0;
   real_t ts_last_write_time = -1e30;        // R-304
   real_t ts_last_v_max       = 0.0;         // R-304
   int    ts_current_regime   = 0;           // R-304
   const bool have_ts_state = ReadPetscTSCheckpoint(
      restart_prefix, ts_t, ts_dt_next, ts_step, ts_rejections,
      pv_snapshots, ts_last_write_time, ts_last_v_max,
      ts_current_regime, &mpi);
   MFEM_VERIFY(have_ts_state,
               "--restart with --petsc-ts requires a V2 checkpoint "
               "(written by a build >= 2026-05-XX).  V1 checkpoints "
               "do not contain PETSc TS state; cannot continue.");

   // V1 ReadCheckpoint already set `t`, `current_dt`, and `state`
   // from the V1 fields.  Cross-check vs. V2 fields:
   MFEM_VERIFY(std::abs(t - ts_t) < 1e-12 * std::abs(t),
               "Checkpoint inconsistency: V1 time=" << t
               << " differs from V2 time=" << ts_t);

   petsc::TS ts = *petsc_ode;
   PetscErrorCode ierr;
   ierr = TSSetTime(ts, ts_t);            PCHKERRQ(ts, ierr);
   ierr = TSSetTimeStep(ts, ts_dt_next);  PCHKERRQ(ts, ierr);
   ierr = TSSetStepNumber(ts, ts_step);   PCHKERRQ(ts, ierr);

   if (pv_out)
   {
      pv_out->SetTotalSnapshotsWritten(pv_snapshots);
      // R-304: also restore the schedule state (last_write_time,
      // last_v_max, current_regime) so the first ShouldWrite after
      // restart does not fire unconditionally and the regime
      // state-machine carries hysteresis across the seam.
      pv_out->RestoreScheduleState(ts_last_write_time,
                                   ts_last_v_max,
                                   ts_current_regime);
   }
   // R-303: PRE-restart rejection count.  The post-Run code at
   // line 2414-2422 will OVERWRITE `step_rejections` with the
   // post-restart count (TSGetStepRejections returns this-Run's
   // count only, not the cumulative).  Save the pre-restart count
   // here in a separate accumulator so the final summary line
   // reports the sum across the restart seam.
   restart_rejections_carryover = ts_rejections;

   if (mpi.IsRoot())
   {
      std::cout << "PETSc TS restart: t=" << ts_t << " s, "
                << "dt_next=" << ts_dt_next << " s, step=" << ts_step
                << ", rejections=" << ts_rejections << ", "
                << "paraview_snapshots=" << pv_snapshots << "\n";
   }
}
#endif
```

Add a local declaration `int restart_rejections_carryover = 0;` next to the existing `int step_rejections = 0;` at `bp5_verification_full.cpp:2223`.

Modify the post-Run rejection-accumulation code at `bp5_verification_full.cpp:2414-2422`:

```diff
       PetscInt rejects = 0;
       TSGetStepRejections(ts, &rejects);
-      step_rejections = static_cast<int>(rejects);
+      // R-303: TSGetStepRejections returns THIS-Run's rejections
+      // only (it isn't reset by TSSetStepNumber but it ISN'T
+      // pre-populated from the checkpoint either — see PETSc
+      // 3.15 `src/ts/interface/ts.c`).  Accumulate with the
+      // pre-restart count so the summary reports the sum.
+      step_rejections = restart_rejections_carryover
+                      + static_cast<int>(rejects);
```

**5. Checkpoint-writing call sites.**

Three sites currently call `WriteCheckpoint`:

- Monitor (line 619-628): every `checkpoint_interval` steps
- Mid-run (line 2588-2594): the non-PetscTS path's per-step checkpointing
- Final (line 2632-2638): once at end-of-run

The monitor and final sites are reached on the PetscTS path. After each `WriteCheckpoint(...)` call (in the monitor and final blocks), add an immediately-following:

```cpp
#ifdef MFEM_USE_PETSC
if (use_petsc_ts && petsc_ode)
{
   petsc::TS ts = *petsc_ode;
   PetscReal ts_dt_next_q;
   PetscInt ts_step_q, ts_rejections_q;
   TSGetTimeStep(ts, &ts_dt_next_q);
   TSGetStepNumber(ts, &ts_step_q);
   TSGetStepRejections(ts, &ts_rejections_q);
   const int    pv_snap       = pv_out ? pv_out->GetTotalSnapshotsWritten() : 0;
   const real_t pv_last_write = pv_out ? pv_out->GetLastWriteTime()        : -1e30;
   const real_t pv_last_vmax  = pv_out ? pv_out->GetLastVMax()             :  0.0;
   const int    pv_regime     = pv_out ? pv_out->GetCurrentRegime()        :  0;
   WritePetscTSCheckpoint(full_prefix, time, ts_dt_next_q,
                          static_cast<int>(ts_step_q),
                          static_cast<int>(ts_rejections_q),
                          pv_snap, pv_last_write, pv_last_vmax,
                          pv_regime, mon->mpi);
}
#endif
```

The mid-run site (non-PetscTS path) does NOT need the TS extension.

**6. `SetTotalSnapshotsWritten` + `RestoreScheduleState` setters and three new getters.**

Add to `paraview_output.hpp` near `GetTotalSnapshotsWritten` (existing read-only getter at line 368):

```cpp
/// @brief Restore the snapshot counter from a checkpoint.
///
/// Use ONLY at restart time, before any Save/ShouldWrite call.  The
/// cap-aware schedule needs the counter to survive restart, otherwise
/// the cap budget resets and a restarted long run will over-write past
/// `max_total_snapshots`.  Clamps negative values to 0 (defends
/// against a corrupted/truncated V2 block).
void SetTotalSnapshotsWritten(int n)
{
   total_snapshots_written_ = std::max(n, 0);
}

/// @brief R-304: Restore the adaptive-schedule state from a checkpoint.
///
/// Pairs with `SetTotalSnapshotsWritten`.  Use ONLY at restart time,
/// before any Save/ShouldWrite call.  Without this restoration the
/// first ShouldWrite after restart fires unconditionally (because
/// `last_write_time_` defaults to -1e30 and `time - (-1e30)` always
/// exceeds dt_out * tol) and the regime state machine resets to
/// interseismic regardless of where the pre-checkpoint trajectory was.
void RestoreScheduleState(real_t last_write_time, real_t last_v_max,
                          int current_regime)
{
   last_write_time_ = last_write_time;
   last_v_max_      = last_v_max;
   current_regime_  = current_regime;
}

/// @brief R-304: read-only accessors for the schedule state.
/// Used by `WritePetscTSCheckpoint` to snapshot the state at
/// checkpoint time.  These were previously private; exposing them
/// is necessary for the restart machinery.
real_t GetLastWriteTime() const { return last_write_time_; }
real_t GetLastVMax()      const { return last_v_max_; }
int    GetCurrentRegime() const { return current_regime_; }
```

### Interfaces

**New exposed functions:**

- `mfem::seas::WritePetscTSCheckpoint(prefix, t, dt_next, step, rejections, pv_snapshots, pv_last_write_time, pv_last_v_max, pv_current_regime, mpi)` — header-only, inline. (R-304: includes the three pv schedule-state fields.)
- `mfem::seas::ReadPetscTSCheckpoint(prefix, t, dt_next, step, rejections, pv_snapshots, pv_last_write_time, pv_last_v_max, pv_current_regime, mpi)` — header-only, inline, returns `bool`. (R-304: includes the three pv schedule-state out-params.)
- `mfem::seas::ParaViewOutput::SetTotalSnapshotsWritten(int)` — public method, clamps negative input; no behaviour change for callers that don't restart.
- `mfem::seas::ParaViewOutput::RestoreScheduleState(real_t, real_t, int)` — public method, restores the schedule's last-write-time, last-v-max, and regime; no behaviour change for callers that don't restart. (R-304)
- `mfem::seas::ParaViewOutput::GetLastWriteTime()`, `GetLastVMax()`, `GetCurrentRegime()` — three new public const accessors, used by `WritePetscTSCheckpoint` to snapshot the schedule state. (R-304)

**Connection to existing code:**

- `WriteCheckpoint` continues to be called as today; `WritePetscTSCheckpoint` follows immediately on the PetscTS path. The two write to the same file.
- `ReadCheckpoint` continues to be called as today; `ReadPetscTSCheckpoint` is called AFTER it on the PetscTS path, returning `false` for legacy V1 files.

### Edge Cases to Handle

- **V1 file fed to a V2-aware run.** `ReadPetscTSCheckpoint` returns `false`. The V2 driver path treats this as a hard error with `MFEM_VERIFY` and a clear message ("V1 checkpoint does not contain PETSc TS state — re-write with a build >= 2026-05-XX or use --no-petsc-ts to restart on the MFEM time-stepper"). Decision: do NOT silently fall back to "restart but TS starts fresh" — the user explicitly asked for `--petsc-ts` and a silent dt-reset would produce a misleading trajectory.

- **Rank-count mismatch.** Already caught by `ReadCheckpoint` (file_num_ranks vs current size) at `checkpoint.hpp:196-201`. No new check needed for the V2 trailing block (V2 is per-rank like V1).

- **Missing PetscTS checkpoint file but V1 present and `--petsc-ts` set.** `ReadCheckpoint` returns `true` (V1 fields), `ReadPetscTSCheckpoint` returns `false`. Driver aborts (same as the V1-only case above).

- **Zero-fault-DOF rank restart.** `ReadCheckpoint` calls `state.SetSize(0)`. The R-003 padded-shrink at line 1572-1591 must run BEFORE `ReadCheckpoint`, so capacity is already >= 1; the SetSize(0) shrink preserves the allocation. Verified by an extension of `seas_test_bp5_petsc_ts_zero_fault_rank` (Phase 1 test 5 below).

- **`dt_next == 0`.** Can happen if TSGetTimeStep is called BEFORE the first step or after `TSSetConvergedReason(TS_DIVERGED_*)`. Detect on read: if `dt_next <= 0`, fall back to `dt_init`. Log a one-line warning on rank 0.

- **`paraview_snapshots > max_total_snapshots`.** Means the previous run already exhausted the cap. Restart should preserve this state — `SetTotalSnapshotsWritten(n)` accepts any non-negative `n`. The cap-exhausted branch will fire on the first `ShouldWrite` call after restart (per `RecomputeIntervalForCap`'s `remaining_budget <= 0` path) which is the correct behaviour.

- **State Vector size mismatch between V1 and V2 sections.** Can't happen — they refer to the same vector. But if somehow the V2 `petsc_ts_time` differs from V1 `time` (file corruption, partial-write), the cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` catches it.

### Acceptance Criteria

- [ ] `seas_bp5_full --restart PREFIX --petsc-ts` no longer aborts with "ERROR: --restart is not supported with --petsc-ts".
- [ ] A V1 checkpoint (written by an old build) read by a V2-aware driver with `--petsc-ts` aborts with a clear message naming the missing V2 fields.
- [ ] A V1 checkpoint read by a V2-aware driver WITHOUT `--petsc-ts` (i.e., MFEM time stepper path) still works exactly as today (V1 backwards compatibility).
- [ ] A V2 checkpoint contains both the V1 block and the `PETSC_TS_V2` trailing block, in that order, in the same file `{prefix}_checkpoint_r{rank}.txt`.
- [ ] Test `seas_test_bp5_petsc_ts_restart` passes: a short synthetic run that checkpoints mid-way and restarts produces a final state matching a fresh-run final state to within atol+rtol*|y|.
- [ ] Test `seas_test_bp5_petsc_ts_zero_fault_rank` continues to pass (R-003 invariant preserved through restart).
- [ ] Existing schedule-cap (23/23) and kinematics (31/31) tests continue to pass (no regression on non-restart paths).
- [ ] `make seas_bp5_full seas_tpv102_driver seas_tpv104_driver seas_tpv205_driver` builds clean (no warnings).
- [ ] The bp5_phase6_paraview_zfp_normal_48hr.sbatch can be re-submitted (after a one-time edit to enable restart on second submission via `--restart`) and resume cleanly.

### Dependencies

- Depends on: nothing (uses only existing public PETSc TS API and existing MFEM PetscODESolver).
- Required by: Phase 2.

---

## Phase 2: TSAdapt-aware dt_next preservation

### Goal

After this phase, the dt PETSc uses for the FIRST step after restart is bit-exact equal to the dt the fresh run would have used at the same `t`, eliminating the small (<1%) trajectory divergence Phase 1 leaves due to the controller "warming up" from a single just-set dt instead of an inherited error history.

### Files to Modify

- `miniapps/seas/io/petsc_ts_checkpoint.hpp` — extend the V2 block with an optional `petsc_ts_adapt_state` field (currently the I-controller has no history, so this is a no-op record reserved for future PI/PID upgrades). Bump tag to `PETSC_TS_V3` only when an actual extension lands; for Phase 2 the V2 tag is reused since the existing `petsc_ts_dt_next` already captures the I-controller's full state.

- `miniapps/seas/tests/unit/test_bp5_petsc_ts_restart.cpp` — tighten the tolerance assertion (Phase 1 was atol + rtol; Phase 2 expects the first-step dt to match bit-exactly).

### Detailed Requirements

For the `TSADAPTBASIC` adapter that the bp5_phase6 sbatch uses (`tests/verification/petsc_ts_rk45_tandem.cfg`: `-ts_type rk -ts_rk_type 5dp`), the controller is an I-only controller with no internal history beyond the dt it just produced. So Phase 1's `petsc_ts_dt_next` field is already sufficient — Phase 2 is a NO-OP IMPLEMENTATION for the BP5 production configuration.

The only Phase 2 work is the TEST tightening: prove that the first-step dt is bit-exact.

```cpp
// In seas_test_bp5_petsc_ts_restart, add:
TEST_ASSERT(dt_first_step_restart == dt_first_step_fresh,
            "Phase 2: first-step dt after restart must be bit-exact "
            "equal to the fresh-run dt at the same simulated time "
            "(TSADAPTBASIC is an I-controller; the saved dt_next IS "
            "the full controller state)");
```

If a future PR changes the adapter to TSADAPTGLEE or TSADAPTHISTORY (both of which have multi-step history), Phase 2 must be expanded to save that history. The plan flags this in the test's comment.

### Acceptance Criteria

- [ ] The first-step dt after restart equals the fresh-run dt at the same `t`, bit-exactly (under TSADAPTBASIC).
- [ ] A test that switches to a hypothetical multi-history adapter would FAIL the bit-exact assertion (manual procedure documented).

### Dependencies

- Depends on: Phase 1.
- Required by: Phase 3 (optional).

---

## Phase 3: FSAL stage preservation (bit-exact restart)

### Goal

After this phase, the restarted trajectory is BIT-EXACTLY identical to a fresh-run trajectory for the same wall-clock window, including the first step after restart (which currently re-evaluates the first DP5 stage instead of inheriting it via FSAL).

### Files to Modify

- `miniapps/seas/io/petsc_ts_checkpoint.hpp` — extend V2 block (or bump to V3) with `petsc_ts_fsal_stage_size N` + N values for the DP5 FSAL stage k_7 = k_1 vector.

- **NEW upstream PETSc PR** required: add BOTH of these to `src/ts/impls/explicit/rk/rk.c`:
  * `PetscErrorCode TSRKGetStageVectors(TS ts, PetscInt *nstages, Vec **Y)` — for the checkpoint WRITE side.  Returns the internal RK stage memory.
  * `PetscErrorCode TSRKSetStageVectors(TS ts, PetscInt nstages, Vec *Y)` — for the checkpoint READ side.  Restores the FSAL stage (Y[6] in DP5) so the first step after restart uses the inherited K[6] instead of re-evaluating the RHS at the seam.
  Both are new public API.  Without `TSRKSetStageVectors`, the Phase 3 read side cannot complete and Phase 3 collapses back to (essentially) Phase 1 behaviour for the first step after restart.

- **Alternative MFEM-only path:** add BOTH `PetscODESolver::GetRKStageVectors(int &n, Vec **Y)` AND `PetscODESolver::SetRKStageVectors(int n, Vec *Y)` to `linalg/petsc.{hpp,cpp}` that reach into the TSRK private headers (`#include "petsc/private/tsimpl.h"`) to extract and restore the stage memory.  This avoids both upstream PRs but creates a private-API dependency that breaks across PETSc versions.

### Detailed Requirements

1. On the WRITE side: call `TSRKGetStageVectors` (the new upstream-PR API, OR the MFEM private-API equivalent — see "Files to Modify" above) on the `ts` handle, extract the FSAL stage vector (for DP5, this is `Y[6]`), and append its values to the V2/V3 checkpoint as `petsc_ts_fsal_stage` + N + values.

2. On the READ side: load the values, allocate a PETSc Vec of the right size, and set `Y[6]` (or equivalent) on the active TS via the matching `TSRKSetStageVectors` (the new upstream-PR API for the READ side, OR the MFEM private-API equivalent).

3. The interaction with non-DP5 RK types (rk2, rk3, rk4, ...) needs an early-out: only DP5 uses FSAL. For other types, skip the FSAL save/restore. Detect via `TSRKGetType(ts, &rk_type)` and compare against `TSRK5DP`.

### Acceptance Criteria

- [ ] The restarted trajectory's first step state-vector update is bit-exactly equal to the fresh-run state-vector update at the same `t`.
- [ ] The integrated total snapshot count over a 1000-step restart-from-step-500 matches a 1000-step fresh run to the byte.
- [ ] Non-DP5 TS types still work (no FSAL save/restore attempted).

### Dependencies

- Depends on: Phase 2.
- Required by: nothing.
- **Blocked on:** PETSc upstream PR acceptance OR the MFEM private-header workaround being deemed acceptable.

---

## Testing Strategy

### Phase 1 — `seas_test_bp5_petsc_ts_restart` (NEW, MPI-only)

Mirrors the `seas_test_bp5_petsc_ts_zero_fault_rank` skeleton (its own `main`, no `MPI_Init` in the schedule-cap test). Requires `MFEM_USE_PETSC`. Run via `mpirun -np 2 ./seas_test_bp5_petsc_ts_restart`. Skips with INFO message on serial / single-rank.

**Sub-test 1: V1 backwards compatibility (no PetscTS).**

```cpp
// Write a V1 checkpoint with the existing WriteCheckpoint (TS field absent).
// Read it back with the new V2-aware ReadCheckpoint.  Verify all V1
// fields round-trip exactly.  ReadPetscTSCheckpoint returns false.
```

**Sub-test 2: V2 round-trip.**

```cpp
// Write V1 + V2 (WriteCheckpoint + WritePetscTSCheckpoint).
// Read both back.  Verify all 12 fields (V1: t, dt, step, num_eq,
// in_eq, state, displacement, traction, slip_rate, fsal_init, k0;
// V2: t, dt_next, step, rejections, paraview_snapshots) round-trip
// exactly (bit-exact for the integer/string fields; 17-digit
// equality for the floats — which IS bit-exact in this format).
```

**Sub-test 3: Fresh-vs-restart trajectory agreement (the core correctness test).**

This is the most important test. Express the run windows as ABSOLUTE end-times (`T_mid`, `T_full`) and use the existing `PetscODESolver::Run(state, t, dt, t_final)` API.  This is the only formulation that gives a deterministic stopping point:
- `Run` calls `TSSolve` with `TS_EXACTFINALTIME_MATCHSTEP` set by MFEM (linalg/petsc.cpp:4365), so the trajectory ALWAYS lands exactly at `t = t_final` (bit-exact under FP).
- "Number of accepted steps" is NOT a deterministic primitive under adaptive RK45 (controller dt is dependent on local error magnitudes), so phrasing the test in terms of step counts would make `t_A == t_B` non-deterministic.

Procedure (in pseudocode):

```
1. Construct a minimal BP5-shaped problem (inline mesh, P=1 IP, ~1000
   fault DOFs, atol=1e-7, rtol=1e-50).  Construct seas_op, fault_op.
   Pick T_full and T_mid such that 0 < T_mid < T_full.  Suggested:
   dt_init = 0.01 s, T_mid = 1.0 s, T_full = 2.0 s (small enough to
   run in <1 s wall-clock per scenario; long enough for the adaptive
   controller to take many real steps).
2. Run scenario A (fresh): create state, SetInitialCondition (deterministic
   initial state), call petsc_ode->Run(state, t, dt, T_full).  Save final
   state_A, final t_A.  By TS_EXACTFINALTIME_MATCHSTEP, t_A == T_full
   bit-exactly.
3. Run scenario B (restart-mid-way):
     a. Create state, SetInitialCondition (identical to A).
     b. Call petsc_ode->Run(state, t, dt, T_mid).  After this, t == T_mid.
     c. Call WriteCheckpoint + WritePetscTSCheckpoint to /tmp/test_restart.
     d. Destroy the PetscODESolver instance, the seas_op, fault_op,
        and state Vector (full teardown — verifies the checkpoint
        round-trip is robust to memory recycling).
     e. Re-construct everything fresh from the same problem-setup code.
     f. Apply the R-003 padded-shrink to state, then call
        ReadCheckpoint + ReadPetscTSCheckpoint.  Apply TSSetTime,
        TSSetTimeStep, TSSetStepNumber from the V2 fields, and
        SetTotalSnapshotsWritten + RestoreScheduleState on pv_out
        (per R-304).
     g. Call petsc_ode->Run(state, t, dt, T_full).
     h. Save final state_B, final t_B.
4. Assert:
   - t_A == T_full and t_B == T_full bit-exactly (both runs honour
     TS_EXACTFINALTIME_MATCHSTEP).  This is the deterministic
     stopping primitive — NOT "200 steps".
   - max_i |state_A(i) - state_B(i)| <= atol + rtol * max_i |state_A(i)|.
     With atol=1e-7, rtol=1e-50, the bound is ~1e-7.  This is the
     standard "correct-within-tolerance" restart contract for
     adaptive ODE solvers — Phase 1 is NOT bit-exact (see Phase 3
     for the bit-exact target).
```

**Sub-test 4: Snapshot counter survives restart.**

```cpp
// Set up a small ParaViewOutput with max_total_snapshots = 10.  Run a
// short scenario A that produces 7 snapshots.  Checkpoint.  Restart.
// Run more snapshots; assert the total never exceeds 10 (the cap
// budget was preserved).
```

**Sub-test 5: Zero-fault-DOF rank restart.**

```cpp
// Use a 2-rank partitioning where rank 1 has 0 fault DOFs.  Restart
// from a V2 checkpoint on this partitioning.  Assert that rank 1's
// state.GetMemory().Empty() is FALSE after ReadCheckpoint (because
// the padded-shrink ran first), and that the PetscTS Run completes
// without the R-003 abort.
```

**Sub-test 6: Rank-count mismatch.**

```cpp
// Write a 2-rank checkpoint.  Try to read it on a 1-rank run.
// Assert MFEM_VERIFY fires with "num_ranks mismatch" message.
```

**Sub-test 7: V1 file rejected when --petsc-ts.**

```cpp
// Write a V1-only checkpoint (no PetscTS extension).  Try to start
// the V2 driver with --restart + --petsc-ts.  Assert MFEM_VERIFY
// fires with the documented "V1 checkpoint does not contain PETSc
// TS state" message.
```

### Phase 2 — Tighten Sub-test 3 to bit-exact first-step dt

```cpp
// In sub-test 3 above, after step (g) runs ONE step, capture the dt
// it just used.  In a separate run, recompute scenario A but stop
// after step K/2 + 1 and capture the dt that step used.  Assert
// dt_step_K/2_plus_1_restart == dt_step_K/2_plus_1_fresh (bit-exact).
```

### Phase 3 — Tighten Sub-test 3 to bit-exact state

```cpp
// In sub-test 3, replace
//   max_i |state_A(i) - state_B(i)| <= atol + rtol * max_i |state_A(i)|
// with
//   for (i = 0..state_A.Size()-1)
//      assert state_A(i) == state_B(i)  // bit-exact
```

### Build / test wiring

```
$ cd miniapps/seas
$ make seas_test_bp5_petsc_ts_restart        # Phase 1+
$ mpirun -np 2 ./seas_test_bp5_petsc_ts_restart
  Sub-test 1 (V1 backwards compat) ... PASS
  Sub-test 2 (V2 round-trip) ........ PASS
  Sub-test 3 (fresh-vs-restart) ..... PASS (within atol)
  Sub-test 4 (snapshot counter) ..... PASS
  Sub-test 5 (zero-fault-rank) ...... PASS
  Sub-test 6 (rank mismatch) ........ PASS
  Sub-test 7 (V1 + --petsc-ts) ...... PASS
  === Summary: 7 / 7 passed; 0 failed ===
```

The new test is MPI-only. Per R-104 / FIX_REPORT decision, it stays OUT of `make test` umbrella (mpirun not portable). It DOES get its own target `test-bp5-petsc-ts-restart` and a one-line README note in `regime_budget_2026-05-16.md` that operators should run it as a pre-submit step.

---

## Risk Assessment

### What could go wrong

| Risk                                                              | Detection                                          | Mitigation                                                              |
|-------------------------------------------------------------------|----------------------------------------------------|-------------------------------------------------------------------------|
| TSADAPTBASIC has hidden multi-step history not captured by `TSGetTimeStep` | Sub-test 3 fails (trajectory divergence > atol)    | Read PETSc 3.15 source `src/ts/adapt/impls/basic/adaptbasic.c` to verify; if history exists, add to V2 schema |
| FSAL stage NOT inherited produces > 1% trajectory drift over many steps | Long-window variant of Sub-test 3 (10k steps not 200) | Phase 3 fix; until then document as expected behaviour                  |
| `TSSetStepNumber` doesn't reset PETSc's internal `_TSMonitorSet` log; output cycle counts drift | Manual check: monitor callback's `step` arg vs. the saved counter | Test that the monitor fires with `step == saved_step + 1` after restart |
| Append-mode writing of V2 block races with WriteCheckpoint's close | Always serial (one rank, one file at a time); no race | None needed                                                             |
| ParaView snapshot files overwritten on restart (cycle numbers reset) | Sub-test would catch — check `<output_dir>/fault.vtkhdf` doesn't lose timesteps | Append-mode for VTKHDF or rename strategy — out of scope for Phase 1; document as known limitation |
| Restart breaks the cap-aware `last_committed_cycle_` dedup        | Sub-test 4 catches double-bump                     | Save `last_committed_cycle_` in V2; restore in `SetTotalSnapshotsWritten` |

### Known tricky areas

1. **`bp5_verification_full.cpp:1572-1591` — R-003 padded-shrink.** Must execute BEFORE `ReadCheckpoint`. The current restart block at lines 2346-2358 already calls `ReadCheckpoint(prefix, ..., state, ...)` where `state` is the padded-shrink'd Vector from line 1572. Verify the new V2-aware restart path preserves this ordering.

2. **Monitor callback's `petsc_mon_ctx.snapshots_so_far` is not currently a field.** Adding `SetTotalSnapshotsWritten` to `ParaViewOutput` is enough; the monitor doesn't need to know about the counter restoration.

3. **`TSGetTimeStep` returns the NEXT-step dt, not the JUST-COMPLETED dt** (per the comment at `petsc.cpp:4345`). This is the dt we want to save for restart — it's the controller's prediction for the next step. Correct.

4. **`TSSetTime` on a TS that was already initialised — does it work?** Yes, PETSc allows it. From PETSc docs: "TSSetTime: Sets the time. Logically Collective on TS. Must be called before TSSolve, OR can be called inside a monitor to manually advance time. Usually called only at problem setup." Calling it AFTER `TSSetFromOptions` but BEFORE `TSSolve` (which is what `petsc_ode->Run()` does internally) is safe.

5. **`TSSetSolution` vs. PlaceMemory.** MFEM's `PetscODESolver::Run` calls `X->PlaceMemory(x.GetMemory(), true)` then `TSSolve(ts, X->x)`. The Vec `X->x` aliases the user's state Vector. After restart, we have a fresh `state` Vector (with the loaded values from `ReadCheckpoint`). The next `petsc_ode->Run(state, ...)` will PlaceMemory afresh, picking up the loaded values. So we do NOT need to call `TSSetSolution` explicitly — the PlaceMemory in Run does it.

6. **Phase 3 FSAL — needs a PETSc PR.** Document this clearly. The Phase 3 work is BLOCKED on either (a) upstream PETSc accepting `TSRKGetStageVectors` or (b) the MFEM-side workaround being accepted.

### Out of scope

- Restart across DIFFERENT MPI rank counts (re-partitioning the mesh + redistributing the state). The existing `ReadCheckpoint` rank-count check forbids this; the V2 extension does not relax it.
- Restart across DIFFERENT mesh files. The existing checkpoint format does not record the mesh; rest is the operator's responsibility.
- Binary (non-text) checkpoint format. The plain-text format is preserved for debuggability; binary would be a separate plan.
- Snapshot-file (`.vtkhdf`) append-mode for restart. The VTKHDF writer currently overwrites; on restart the fault.vtkhdf for the second window would clobber the first. A workaround for now: rename the previous `fault.vtkhdf` to `fault_part1.vtkhdf` manually before submitting the restart, OR set `--output-dir` to a NEW directory on restart. Both are documented in the sbatch comment.

---

## How to start implementing (for the next agent)

1. Read this plan top to bottom.
2. Confirm with `git log --oneline -5` that you're on `feature/paraview-compaction` at or after commit `407f456` ("BP5 Phase 6 ParaView: fix ResetMemory crash + cap rework + per-regime cadence").
3. Implement Phase 1 in a SINGLE atomic commit; do not bleed it into Phase 2 or 3.
4. Run `make seas_test_bp5_petsc_ts_restart` + `mpirun -np 2 ./seas_test_bp5_petsc_ts_restart` and confirm 7/7 PASS before marking Phase 1 done.
5. Run `make test-paraview-schedule-cap test-bp5-petsc-ts-zero-fault-rank test-kinematics-field-set` to confirm no regression on the previously-fixed bugs.
6. Submit a new test sbatch for Frontera (a derivative of `bp5_phase6_paraview_zfp_normal_48hr.sbatch` that runs for ~2 h with `--checkpoint-interval 1000`) to gather a V2 checkpoint. Then submit a follow-up `--restart` job and verify the trajectory continues cleanly.
7. Defer Phases 2 and 3; they are documented for a later pass. The R-203 cleanup (round-3 review) is unrelated to this plan but should be in place; verify `RecomputeIntervalForCap` does NOT take a `regime` parameter (post-R-203 cleanup state).

---

## See also

- `regime_budget_2026-05-16.md` (same directory) — empirical write-count budget for BP5 production runs.
- `FIX_REPORT.md` (project root) — round-3 fixes that landed on this branch.
- `REVIEW.md` (project root) — round-3 review findings already addressed.
- `miniapps/seas/io/checkpoint.hpp` — existing V1 checkpoint format.
- `miniapps/seas/tests/verification/bp5_verification_full.cpp` lines 1091-1099 (the abort to delete), 2346-2358 (V1 restart path that the V2 extension hooks into).
- `linalg/petsc.cpp` lines 4357-4394 (`PetscODESolver::Run` — the PlaceMemory/ResetMemory round-trip that V2 plugs into).
