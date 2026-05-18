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
  - **Add** the V2 restart block IMMEDIATELY AFTER the existing V1 restart block at lines 2346-2358.  **Ordering is load-bearing**: the V1 `ReadCheckpoint` call at line 2352 populates `t`, `current_dt`, `state`, `num_seismic_events`, `in_seismic_event` from the V1 fields; the V2 block then reads the 10 PETSC_TS_V2 fields, cross-checks `t == ts_t`, calls `TSSetStepNumber(ts, ts_step)`, and writes `t = ts_t; current_dt = ts_dt_next;` so that PetscODESolver::Run picks up the V2 values when it calls TSSetTime/TSSetTimeStep on entry (R-002: the explicit TSSetTime / TSSetTimeStep calls have been REMOVED because `Run` overwrites them).  The cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` is only meaningful BECAUSE `t` has been populated by the preceding V1 read.  **Do NOT place the V2 block inside the PetscTS init at line ~2269** — `t` is still 0 there, the cross-check would always fail, and every restart attempt would abort with a misleading "V1 time=0 differs from V2 time=..." error.
  - **Modify** the PetscTS-path checkpoint-writing call sites — the monitor (lines 619-628) and the final checkpoint in main (lines 2631-2641) — to ALSO call `WritePetscTSCheckpoint` immediately after `WriteCheckpoint`.  The mid-run site at line 2591 is in the non-PetscTS `else` branch and does NOT get the V2 extension (R-009).  The fields come from: `TSGetTime`, `TSGetTimeStep`, `TSGetStepNumber`, `TSGetStepRejections` on the `ts` handle (monitor uses its `ts` parameter directly; main uses `*petsc_ode`); `snapshots_so_far = pv_out ? pv_out->GetTotalSnapshotsWritten() : 0`; per R-304, also `pv_out->GetLastWriteTime()`, `pv_out->GetLastVMax()`, `pv_out->GetCurrentRegime()`; per R-004, `pv_out->GetLastCommittedCycle()`; per R-006, `pv_out->GetLastVolumeWriteTime()` (new accessors — see item 6 below).

  - **Extend** `BP5MonitorCtx` (defined at `bp5_verification_full.cpp:483-513`) so the monitor callback can reach the data the V2 checkpoint needs.  The monitor sees ONLY its `void *ctx` (→ `BP5MonitorCtx*`) and the `TS ts` parameter — `pv_out`, `use_petsc_ts`, `petsc_ode`, and `restart_rejections_carryover` are all local to `main` and OUT OF SCOPE in the monitor (R-001).  Add two members at the end of the struct:
    ```cpp
    // V2 PETSc-TS restart support (R-001 / R-005).
    seas::ParaViewOutput<ParMesh> *pv_out = nullptr;     // may be nullptr
    int restart_rejections_carryover = 0;                 // R-303 / R-005
    ```
    And initialise both in main, alongside the existing `petsc_mon_ctx.* = ...` block at lines 2296-2313:
    ```cpp
    petsc_mon_ctx.pv_out = pv_out.get();                   // R-001
    petsc_mon_ctx.restart_rejections_carryover = restart_rejections_carryover;  // R-005
    ```
    The monitor-side V2 snippet (see §5 below) then dereferences `mon->pv_out` and `mon->restart_rejections_carryover` instead of the out-of-scope main locals, and the TS handle comes from the callback's `ts` parameter directly (no `*petsc_ode` indirection).
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
///   petsc_ts_time                  <real_t>
///   petsc_ts_dt_next               <real_t>
///   petsc_ts_step                  <int>
///   petsc_ts_rejections            <int>  (cumulative across restarts; R-005)
///   paraview_snapshots             <int>
///   paraview_last_write_time       <real_t>   (R-304)
///   paraview_last_v_max            <real_t>   (R-304)
///   paraview_current_regime        <int>      (R-304)
///   paraview_last_committed_cycle  <int>      (R-004)
///   paraview_last_volume_write_time <real_t>  (R-006)
///
/// `petsc_ts_dt_next` is the dt PETSc will use for the NEXT step, as
/// returned by TSGetTimeStep AFTER the just-completed step.  This is
/// the value we need to feed back via TSSetTimeStep on restart.
///
/// `petsc_ts_rejections` is the CUMULATIVE rejection count across the
/// entire restart chain — the WRITE-side caller is responsible for
/// adding `restart_rejections_carryover` + `TSGetStepRejections(ts)`
/// before passing the sum here.  Saving the per-run value alone would
/// lose all prior restarts' rejections on a 2+-link chain (R-005).
///
/// `paraview_last_write_time`, `paraview_last_v_max`,
/// `paraview_current_regime` capture the ParaViewOutput adaptive
/// schedule state at checkpoint time so the restart's first
/// ShouldWrite call does not fire unconditionally (the default
/// `last_write_time_ = -1e30` would otherwise make any `time - last`
/// look huge and trigger an extra snapshot at t=saved_t).  See R-304.
///
/// `paraview_last_committed_cycle` (R-004) is the dedup key used by
/// ParaViewOutput::CommitSchedule's same-step guard.  Without it,
/// `last_committed_cycle_` stays at INT_MIN on restart and the first
/// post-restart CommitSchedule cannot dedup against the (-INT_MIN)
/// default; the snapshot counter then over-bumps by 1 per restart event.
///
/// `paraview_last_volume_write_time` (R-006) preserves the independent
/// volume-PV cadence across restart.  Without it the first
/// ForceSaveImpl after restart fires the volume save unconditionally
/// because `last_volume_write_time_` defaults to -1e30 and
/// `time - (-1e30) ~ +inf` always exceeds `volume_pv_dt_ * tol`.
inline void WritePetscTSCheckpoint(const std::string &prefix,
                                   real_t t, real_t dt_next,
                                   int step, int rejections,
                                   int paraview_snapshots,
                                   real_t paraview_last_write_time,
                                   real_t paraview_last_v_max,
                                   int paraview_current_regime,
                                   int paraview_last_committed_cycle,
                                   real_t paraview_last_volume_write_time,
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
out << "paraview_last_committed_cycle "
    << paraview_last_committed_cycle << "\n";        // R-004
out << "paraview_last_volume_write_time "
    << paraview_last_volume_write_time << "\n";      // R-006

// R-012: mirror the WriteCheckpoint pattern — explicit close, then
// barrier — so the "after this returns, all ranks have committed"
// contract holds for the V2 trailing block too.  Without it, rank N
// may still be appending when rank 0 advances to the next
// checkpoint iteration and TRUNC-opens the file.
out.close();
if (mpi) { mpi->Barrier(); }
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
                                  int &paraview_last_committed_cycle,
                                  real_t &paraview_last_volume_write_time,
                                  const MPIContext *mpi);
#endif
```

Body: reopens the file, scans past the V1 block by reading whitespace-separated tokens until the `PETSC_TS_V2` tag is found or EOF is hit (numeric tokens in the V1 body cannot match the alpha-numeric tag, so the scan is robust).  On EOF, return `false` (the file is V1-only).  On a found tag, parse the 10 fields IN THE ORDER WRITTEN by §1 above — `petsc_ts_time`, `petsc_ts_dt_next`, `petsc_ts_step`, `petsc_ts_rejections`, `paraview_snapshots`, `paraview_last_write_time`, `paraview_last_v_max`, `paraview_current_regime`, `paraview_last_committed_cycle`, `paraview_last_volume_write_time` (R-008 + R-004 + R-006).  Each field uses the same `read_tag(expected_name); in >> value;` pattern as `ReadCheckpoint`.  Abort via `MFEM_VERIFY` on any tag mismatch.

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
   int    ts_last_committed_cycle =
             std::numeric_limits<int>::min();   // R-004
   real_t ts_last_volume_write_time = -1e30; // R-006
   const bool have_ts_state = ReadPetscTSCheckpoint(
      restart_prefix, ts_t, ts_dt_next, ts_step, ts_rejections,
      pv_snapshots, ts_last_write_time, ts_last_v_max,
      ts_current_regime, ts_last_committed_cycle,
      ts_last_volume_write_time, &mpi);
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
   // R-002: PetscODESolver::Run() unconditionally calls
   // TSSetTime(ts, t) and TSSetTimeStep(ts, dt) on entry
   // (linalg/petsc.cpp:4362-4363).  An explicit TSSetTime /
   // TSSetTimeStep call HERE would be silently overwritten on
   // the next call to `petsc_ode->Run(state, t, current_dt,
   // t_final)`.  Instead, update the C++ `t` and `current_dt`
   // variables that Run() reads on entry — those are the
   // load-bearing ones.  Only TSSetStepNumber survives Run()
   // (Run never resets the step counter), so it stays.
   ierr = TSSetStepNumber(ts, static_cast<PetscInt>(ts_step));  // R-010
   PCHKERRQ(ts, ierr);

   // V1 ReadCheckpoint already set `t = ts_t_v1`; the cross-check
   // above guarantees ts_t == t, so the reassignment is a no-op
   // today.  But `current_dt` was set to V1's restart_dt, which
   // may diverge from ts_dt_next in any future change that adds a
   // CFL clamp or dt_init override.  Make V2 authoritative.
   t          = ts_t;
   current_dt = ts_dt_next;

   // R-003: TSGetTimeStep can return 0 if the checkpoint was
   // written before the first accepted step or just after a
   // TSSetConvergedReason(TS_DIVERGED_*).  Fall back to dt_init
   // and log on rank 0 so PETSc has a non-zero starting dt.
   if (current_dt <= 0.0)
   {
      if (mpi.IsRoot())
      {
         std::cout << "PETSc TS restart: V2 ts_dt_next was "
                   << current_dt << " <= 0; falling back to "
                   << "dt_init = " << dt_init << " s\n";
      }
      current_dt = dt_init;
   }

   if (pv_out)
   {
      pv_out->SetTotalSnapshotsWritten(pv_snapshots);
      // R-304: also restore the schedule state (last_write_time,
      // last_v_max, current_regime) so the first ShouldWrite after
      // restart does not fire unconditionally and the regime
      // state-machine carries hysteresis across the seam.  R-006:
      // RestoreScheduleState now also takes last_volume_write_time
      // so the volume-PV cadence survives restart (without this,
      // the first ForceSaveImpl after restart fires the volume save
      // unconditionally because last_volume_write_time_ defaults to
      // -1e30 and time - (-1e30) ~ +inf).
      pv_out->RestoreScheduleState(ts_last_write_time,
                                   ts_last_v_max,
                                   ts_current_regime,
                                   ts_last_volume_write_time);
      // R-004: restore the dedup cycle key so CommitSchedule's
      // same-step guard works correctly on the first post-restart
      // commit (otherwise total_snapshots_written_ over-bumps by 1
      // per restart event).
      pv_out->SetLastCommittedCycle(ts_last_committed_cycle);
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

- Monitor (line 619-628): every `checkpoint_interval` steps  → V2 extension required
- Mid-run (line 2588-2594): the non-PetscTS path's per-step checkpointing  → NO V2 extension (R-009)
- Final (line 2632-2638): once at end-of-run  → V2 extension required on the PetscTS path

The monitor and final sites need DIFFERENT snippets because they have DIFFERENT identifiers in scope:

- the monitor sees only its `TS ts` parameter and `mon = static_cast<BP5MonitorCtx*>(ctx)`, so it MUST go through `mon->pv_out` / `mon->mpi` / `mon->restart_rejections_carryover` (R-001) and use the parameter `ts` directly (no `petsc_ode` indirection);
- the final-checkpoint site is in `main`, where `petsc_ode`, `use_petsc_ts`, `pv_out`, and `restart_rejections_carryover` are all in scope as locals.

**Monitor site** — replace the existing `WriteCheckpoint(...)` call at lines 619-628 with the existing call PLUS the following V2 trailing block immediately after:

```cpp
#ifdef MFEM_USE_PETSC
{
   // Note: the monitor is installed ONLY when --petsc-ts is active,
   // so the use_petsc_ts gate is structurally true and elided.
   // `ts` is the callback's TS parameter — no *petsc_ode indirection
   // (petsc_ode is local to main and OUT OF SCOPE here; R-001).
   PetscReal ts_dt_next_q;
   PetscInt  ts_step_q, ts_rejections_q;
   TSGetTimeStep(ts, &ts_dt_next_q);
   TSGetStepNumber(ts, &ts_step_q);
   TSGetStepRejections(ts, &ts_rejections_q);
   const int    pv_snap          = mon->pv_out ? mon->pv_out->GetTotalSnapshotsWritten() : 0;
   const real_t pv_last_write    = mon->pv_out ? mon->pv_out->GetLastWriteTime()         : -1e30;
   const real_t pv_last_vmax     = mon->pv_out ? mon->pv_out->GetLastVMax()              :  0.0;
   const int    pv_regime        = mon->pv_out ? mon->pv_out->GetCurrentRegime()         :  0;
   const int    pv_last_commit   = mon->pv_out ? mon->pv_out->GetLastCommittedCycle()    : std::numeric_limits<int>::min();   // R-004
   const real_t pv_last_vol_time = mon->pv_out ? mon->pv_out->GetLastVolumeWriteTime()   : -1e30;                              // R-006
   // R-005: save the CUMULATIVE rejection count, not this-run's
   // count, so chained restarts do not lose prior runs' rejections.
   const int    cum_rejects      = mon->restart_rejections_carryover
                                   + static_cast<int>(ts_rejections_q);
   WritePetscTSCheckpoint(mon->full_prefix, time, ts_dt_next_q,
                          static_cast<int>(ts_step_q),
                          cum_rejects,
                          pv_snap, pv_last_write, pv_last_vmax,
                          pv_regime, pv_last_commit, pv_last_vol_time,
                          mon->mpi);
}
#endif
```

**Final site** — replace the existing `WriteCheckpoint(...)` call at lines 2632-2641 with the existing call PLUS the following V2 trailing block immediately after:

```cpp
#ifdef MFEM_USE_PETSC
if (use_petsc_ts && petsc_ode)
{
   petsc::TS ts = *petsc_ode;
   PetscReal ts_dt_next_q;
   PetscInt  ts_step_q, ts_rejections_q;
   TSGetTimeStep(ts, &ts_dt_next_q);
   TSGetStepNumber(ts, &ts_step_q);
   TSGetStepRejections(ts, &ts_rejections_q);
   const int    pv_snap          = pv_out ? pv_out->GetTotalSnapshotsWritten() : 0;
   const real_t pv_last_write    = pv_out ? pv_out->GetLastWriteTime()         : -1e30;
   const real_t pv_last_vmax     = pv_out ? pv_out->GetLastVMax()              :  0.0;
   const int    pv_regime        = pv_out ? pv_out->GetCurrentRegime()         :  0;
   const int    pv_last_commit   = pv_out ? pv_out->GetLastCommittedCycle()    : std::numeric_limits<int>::min();   // R-004
   const real_t pv_last_vol_time = pv_out ? pv_out->GetLastVolumeWriteTime()   : -1e30;                              // R-006
   // R-005: save the CUMULATIVE rejection count, not this-run's
   // count, so chained restarts do not lose prior runs' rejections.
   const int    cum_rejects      = restart_rejections_carryover
                                   + static_cast<int>(ts_rejections_q);
   WritePetscTSCheckpoint(full_prefix, t, ts_dt_next_q,
                          static_cast<int>(ts_step_q),
                          cum_rejects,
                          pv_snap, pv_last_write, pv_last_vmax,
                          pv_regime, pv_last_commit, pv_last_vol_time,
                          &mpi);
}
#endif
```

The mid-run site (non-PetscTS path) does NOT need the TS extension.

**6. `SetTotalSnapshotsWritten` + `RestoreScheduleState` + `SetLastCommittedCycle` setters and five new getters.**

Add to `paraview_output.hpp` near `GetTotalSnapshotsWritten` (existing read-only getter at line 418):

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

/// @brief R-304 + R-006: Restore the adaptive-schedule state from a
/// checkpoint.
///
/// Pairs with `SetTotalSnapshotsWritten`.  Use ONLY at restart time,
/// before any Save/ShouldWrite call.  Without this restoration the
/// first ShouldWrite after restart fires unconditionally (because
/// `last_write_time_` defaults to -1e30 and `time - (-1e30)` always
/// exceeds dt_out * tol) and the regime state machine resets to
/// interseismic regardless of where the pre-checkpoint trajectory was.
///
/// R-006: `last_volume_write_time` preserves the independent volume-PV
/// cadence across restart — without it the first ForceSaveImpl after
/// restart fires the volume save unconditionally for the same -1e30
/// reason as above.
void RestoreScheduleState(real_t last_write_time, real_t last_v_max,
                          int current_regime,
                          real_t last_volume_write_time)
{
   last_write_time_        = last_write_time;
   last_v_max_             = last_v_max;
   // R-007: clamp to the documented [0,2] range to defend against a
   // corrupted or forward-incompatible V2 block.  NextRegime's switch
   // has a default-0 branch that would self-recover on the next call,
   // but the assignment itself should not silently install an invalid
   // sentinel.
   current_regime_         = std::clamp(current_regime, 0, 2);
   last_volume_write_time_ = last_volume_write_time;          // R-006
}

/// @brief R-004: restore the dedup-key cycle counter used by
/// CommitSchedule's same-step guard.
///
/// Use ONLY at restart time, alongside SetTotalSnapshotsWritten and
/// RestoreScheduleState.  Without this restoration the first
/// post-restart CommitSchedule cannot dedup against the (-INT_MIN)
/// default sentinel and `total_snapshots_written_` over-bumps by 1
/// per restart event (this leaks budget for long-running cap-aware
/// jobs that restart many times).
void SetLastCommittedCycle(int cycle) { last_committed_cycle_ = cycle; }

/// @brief R-304 + R-004 + R-006: read-only accessors for the schedule
/// state.  Used by `WritePetscTSCheckpoint` to snapshot the state at
/// checkpoint time.  These fields were previously private; exposing
/// them is necessary for the restart machinery.
real_t GetLastWriteTime()       const { return last_write_time_; }
real_t GetLastVMax()            const { return last_v_max_; }
int    GetCurrentRegime()       const { return current_regime_; }
int    GetLastCommittedCycle()  const { return last_committed_cycle_; }      // R-004
real_t GetLastVolumeWriteTime() const { return last_volume_write_time_; }    // R-006
```

### Interfaces

**New exposed functions:**

- `mfem::seas::WritePetscTSCheckpoint(prefix, t, dt_next, step, rejections, pv_snapshots, pv_last_write_time, pv_last_v_max, pv_current_regime, pv_last_committed_cycle, pv_last_volume_write_time, mpi)` — header-only, inline.  10 fields (5 PetscTS + 5 paraview-schedule).  `rejections` is the CUMULATIVE count across the restart chain (R-005).
- `mfem::seas::ReadPetscTSCheckpoint(prefix, t, dt_next, step, rejections, pv_snapshots, pv_last_write_time, pv_last_v_max, pv_current_regime, pv_last_committed_cycle, pv_last_volume_write_time, mpi)` — header-only, inline, returns `bool`.  Same 10 fields as out-params.  R-004 + R-006: extended schema vs the original R-304 draft.
- `mfem::seas::ParaViewOutput::SetTotalSnapshotsWritten(int)` — public method, clamps negative input; no behaviour change for callers that don't restart.
- `mfem::seas::ParaViewOutput::RestoreScheduleState(real_t, real_t, int, real_t)` — public method, restores the schedule's last-write-time, last-v-max, regime, and last-volume-write-time; no behaviour change for callers that don't restart.  R-007: clamps `regime` to [0,2].  R-006: 4th param added.
- `mfem::seas::ParaViewOutput::SetLastCommittedCycle(int)` — public method, restores the dedup cycle key (R-004).  No behaviour change for callers that don't restart.
- `mfem::seas::ParaViewOutput::GetLastWriteTime()`, `GetLastVMax()`, `GetCurrentRegime()`, `GetLastCommittedCycle()`, `GetLastVolumeWriteTime()` — five public const accessors, used by `WritePetscTSCheckpoint` to snapshot the schedule state.  R-004 / R-006 added the last two.

**Connection to existing code:**

- `WriteCheckpoint` continues to be called as today; `WritePetscTSCheckpoint` follows immediately on the PetscTS path. The two write to the same file.
- `ReadCheckpoint` continues to be called as today; `ReadPetscTSCheckpoint` is called AFTER it on the PetscTS path, returning `false` for legacy V1 files.

### Edge Cases to Handle

- **V1 file fed to a V2-aware run.** `ReadPetscTSCheckpoint` returns `false`. The V2 driver path treats this as a hard error with `MFEM_VERIFY` and a clear message ("V1 checkpoint does not contain PETSc TS state — re-write with a build >= 2026-05-XX or use --no-petsc-ts to restart on the MFEM time-stepper"). Decision: do NOT silently fall back to "restart but TS starts fresh" — the user explicitly asked for `--petsc-ts` and a silent dt-reset would produce a misleading trajectory.

- **Rank-count mismatch.** Already caught by `ReadCheckpoint` (file_num_ranks vs current size) at `checkpoint.hpp:196-201`. No new check needed for the V2 trailing block (V2 is per-rank like V1).

- **Missing PetscTS checkpoint file but V1 present and `--petsc-ts` set.** `ReadCheckpoint` returns `true` (V1 fields), `ReadPetscTSCheckpoint` returns `false`. Driver aborts (same as the V1-only case above).

- **Zero-fault-DOF rank restart.** `ReadCheckpoint` calls `state.SetSize(0)`. The R-003 padded-shrink at line 1572-1591 must run BEFORE `ReadCheckpoint`, so capacity is already >= 1; the SetSize(0) shrink preserves the allocation. Verified by an extension of `seas_test_bp5_petsc_ts_zero_fault_rank` (Phase 1 test 5 below).

- **`dt_next == 0`.** Can happen if TSGetTimeStep is called BEFORE the first step or after `TSSetConvergedReason(TS_DIVERGED_*)`. Detect in the V2 driver block (§4): if `current_dt <= 0` after `current_dt = ts_dt_next`, fall back to `dt_init` and log a one-line warning on rank 0.  R-003 implementation lives in §4; do NOT push this guard inside `ReadPetscTSCheckpoint` (the reader is value-blind and should faithfully return the stored value).

- **`paraview_snapshots > max_total_snapshots`.** Means the previous run already exhausted the cap. Restart should preserve this state — `SetTotalSnapshotsWritten(n)` accepts any non-negative `n`. The cap-exhausted branch will fire on the first `ShouldWrite` call after restart (per `RecomputeIntervalForCap`'s `remaining_budget <= 0` path) which is the correct behaviour.

- **State Vector size mismatch between V1 and V2 sections.** Can't happen — they refer to the same vector. But if somehow the V2 `petsc_ts_time` differs from V1 `time` (file corruption, partial-write), the cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` catches it.

### Acceptance Criteria

- [ ] `seas_bp5_full --restart PREFIX --petsc-ts` no longer aborts with "ERROR: --restart is not supported with --petsc-ts".
- [ ] A V1 checkpoint (written by an old build) read by a V2-aware driver with `--petsc-ts` aborts with a clear message naming the missing V2 fields.
- [ ] A V1 checkpoint read by a V2-aware driver WITHOUT `--petsc-ts` (i.e., MFEM time stepper path) still works exactly as today (V1 backwards compatibility).
- [ ] A V2 checkpoint contains both the V1 block and the `PETSC_TS_V2` trailing block, in that order, in the same file `{prefix}_checkpoint_r{rank}.txt`.  The V2 trailing block has exactly 10 fields in the order documented in §1 (R-008): `petsc_ts_time`, `petsc_ts_dt_next`, `petsc_ts_step`, `petsc_ts_rejections` (cumulative; R-005), `paraview_snapshots`, `paraview_last_write_time`, `paraview_last_v_max`, `paraview_current_regime`, `paraview_last_committed_cycle` (R-004), `paraview_last_volume_write_time` (R-006).
- [ ] Test `seas_test_bp5_petsc_ts_restart` passes: a short synthetic run that checkpoints mid-way and restarts produces a final state matching a fresh-run final state to within atol+rtol*|y|.
- [ ] Test `seas_test_bp5_petsc_ts_zero_fault_rank` continues to pass (R-003-padded-shrink invariant preserved through restart — NOTE: this is the OLDER R-003 from the round-1 review, not the R-003 dt_next-guard in this plan).
- [ ] R-002 regression: corrupting V1.dt while leaving V2.ts_dt_next correct, the restart still produces the correct trajectory (proves V2 actually drives the post-restart dt rather than being silently overwritten by `PetscODESolver::Run`).
- [ ] R-003 regression: a V2 checkpoint with `ts_dt_next = 0` triggers the `dt_init` fallback path and prints the documented rank-0 warning.
- [ ] R-004 regression: `total_snapshots_written_` after restart matches the fresh-run count (no over-bump by 1 per restart).
- [ ] R-005 regression: in a 3-link restart chain (A → B → C), `step_rejections` reported at the end of C is the SUM of A+B+C rejections.
- [ ] R-006 regression: with `--volume-pv-dt 0.1`, the number of volume snapshots in window `[T_mid, T_full]` matches between a fresh run and a restart-from-T_mid run (no extra snapshot at the seam).
- [ ] R-007 regression: a V2 checkpoint with `paraview_current_regime = 7` is loaded successfully and the in-memory `current_regime_` is clamped to 0.
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
     f. Apply the R-003-padded-shrink (round-1; the older R-003) to
        `state`, then call ReadCheckpoint + ReadPetscTSCheckpoint.
        Apply `TSSetStepNumber(ts, ts_step)` (the only one Run does
        not overwrite, per R-002), and set `t = ts_t; current_dt =
        ts_dt_next; if (current_dt <= 0) current_dt = dt_init;`
        (per R-002 + R-003), then `pv_out->SetTotalSnapshotsWritten`,
        `pv_out->RestoreScheduleState(..., ts_last_volume_write_time)`
        (R-006 fourth arg), and `pv_out->SetLastCommittedCycle(...)`
        (R-004).
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

R-011: dropped from this test binary.  Inside an `mpirun -np 2` launch, switching from 2 ranks to 1 rank requires spawning a `system("mpirun -np 1 ./helper")` subprocess or leaking into MPIContext internals — both fragile.  The V1 `ReadCheckpoint` already catches rank-count mismatch at `checkpoint.hpp:196-201` via `MFEM_VERIFY(file_num_ranks == size, ...)`, and V2 inherits this because the V2 trailing block is per-rank like V1 (mismatch is detected before the V2 reader is even invoked).  The existing `seas_test_bp5_petsc_ts_zero_fault_rank` already exercises the rank-count guard implicitly; no new dedicated test for V2 is added.

**Sub-test 7: V1 file rejected when --petsc-ts.**

```cpp
// Write a V1-only checkpoint (no PetscTS extension).  Try to start
// the V2 driver with --restart + --petsc-ts.  Assert MFEM_VERIFY
// fires with the documented "V1 checkpoint does not contain PETSc
// TS state" message.
```

**Sub-test 8: `dt_next == 0` fallback (R-003).**

```cpp
// Write a V2 checkpoint with ts_dt_next = 0.0 (synthesise by hand,
// or fish out a real one from a forced TS_DIVERGED checkpoint).
// Restart with that file; capture rank-0 stdout.  Assert:
//   (a) the run completes without aborting;
//   (b) rank-0 log contains "V2 ts_dt_next was 0 <= 0; falling back
//       to dt_init";
//   (c) the first post-restart dt equals dt_init.
```

**Sub-test 9: Multi-restart rejection accumulation (R-005).**

```cpp
// Three-link chain: run A from t=0 to T_A, checkpoint, restart to
// T_B, checkpoint, restart to T_C.  Read PETSc rejection counts
// from each run's summary line.  Assert:
//   (a) n_AB >= n_A and n_ABC >= n_AB (monotonicity);
//   (b) the V2 rejections field in prefix_B equals n_AB (cumulative
//       up to that checkpoint), NOT just n_B alone.
```

**Sub-test 10: Volume-PV cadence survives restart (R-006).**

```cpp
// Set --volume-pv-dt 0.1, --paraview-volume-vtu (or hdf5).  Run
// scenario A from t=0 to 2.0 s, checkpoint at 1.0 s.  Restart from
// 1.0 s, run to 2.0 s.  Assert:
//   snap_count_in_window(prefix_A, 1.0, 2.0)
// == snap_count_in_window(prefix_B, 1.0, 2.0)
// (i.e. no extra snapshot at the post-restart seam from
//  last_volume_write_time_ defaulting to -1e30).
```

**Sub-test 11: Regime clamp on V2 corruption (R-007).**

```cpp
// Hand-write a V2 trailing block with paraview_current_regime = 7.
// Restart; assert ReadPetscTSCheckpoint returns true (the value is
// faithfully read) AND pv_out->GetCurrentRegime() == 0 (the setter
// clamped it).  Sanity-check that the next ShouldWrite call doesn't
// take the invalid-regime path.
```

**Sub-test 12: First post-restart commit doesn't over-bump snapshot counter (R-004).**

```cpp
// Run scenario A producing N writes, checkpoint with
// total_snapshots_written_ = N and last_committed_cycle_ = C.
// Restart, then call CommitSchedule at exactly the saved
// last_write_time_ (within tol).  Assert:
//   pv_out->GetTotalSnapshotsWritten() == N  (NOT N+1, i.e. the
//   dedup correctly recognised "same step as last commit" thanks
//   to SetLastCommittedCycle(C) being called in the V2 restore).
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
  Sub-test 1  (V1 backwards compat) ............ PASS
  Sub-test 2  (V2 round-trip) .................. PASS
  Sub-test 3  (fresh-vs-restart) ............... PASS (within atol)
  Sub-test 4  (snapshot counter) ............... PASS
  Sub-test 5  (zero-fault-rank) ................ PASS
  Sub-test 6  (rank mismatch — see R-011 note) . SKIP (covered by V1 path)
  Sub-test 7  (V1 + --petsc-ts) ................ PASS
  Sub-test 8  (dt_next=0 fallback / R-003) ..... PASS
  Sub-test 9  (multi-restart rejects / R-005) .. PASS
  Sub-test 10 (volume-PV cadence / R-006) ...... PASS
  Sub-test 11 (regime clamp / R-007) ........... PASS
  Sub-test 12 (commit-cycle dedup / R-004) ..... PASS
  === Summary: 11 / 11 passed; 0 failed; 1 skipped (sub-test 6) ===
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
| Append-mode writing of V2 block races with the NEXT WriteCheckpoint's TRUNC-open on a different rank | Each rank owns its file; intra-rank order is sequential | R-012: WritePetscTSCheckpoint ends with explicit `close()` + `Barrier()`, mirroring WriteCheckpoint, so the "all ranks committed before next checkpoint" invariant holds for the V2 block too |
| ParaView snapshot files overwritten on restart (cycle numbers reset) | Sub-test would catch — check `<output_dir>/fault.vtkhdf` doesn't lose timesteps | Append-mode for VTKHDF or rename strategy — out of scope for Phase 1; document as known limitation |
| Restart breaks the cap-aware `last_committed_cycle_` dedup        | Sub-test 12 catches the over-bump-by-1             | R-004: save `paraview_last_committed_cycle` as the 9th V2 field; new `pv_out->SetLastCommittedCycle(...)` setter restores it at V2 read time |

### Known tricky areas

1. **`bp5_verification_full.cpp:1572-1591` — R-003 padded-shrink.** Must execute BEFORE `ReadCheckpoint`. The current restart block at lines 2346-2358 already calls `ReadCheckpoint(prefix, ..., state, ...)` where `state` is the padded-shrink'd Vector from line 1572. Verify the new V2-aware restart path preserves this ordering.

2. **Monitor callback's `petsc_mon_ctx.snapshots_so_far` is not currently a field.** Adding `SetTotalSnapshotsWritten` to `ParaViewOutput` is enough; the monitor doesn't need to know about the counter restoration.

3. **`TSGetTimeStep` returns the NEXT-step dt, not the JUST-COMPLETED dt** (per the comment at `petsc.cpp:4345`). This is the dt we want to save for restart — it's the controller's prediction for the next step. Correct.

4. **`TSSetTime` / `TSSetTimeStep` on a TS that was already initialised — does it work?  R-002 caveat.**  Yes, PETSc allows it, BUT `MFEM::PetscODESolver::Run()` unconditionally calls `TSSetTime(ts, t)` and `TSSetTimeStep(ts, dt)` on entry (`linalg/petsc.cpp:4362-4363`), using the C++ `t` / `dt` args you pass in.  So any explicit `TSSetTime(ts, ts_t)` or `TSSetTimeStep(ts, ts_dt_next)` you call BEFORE `Run` is silently overwritten the moment `Run` starts.  The load-bearing thing to update at restart time is the C++ `t` and `current_dt` variables themselves (NOT just the TS internal state).  `TSSetStepNumber` is the only of the three that `Run` does NOT overwrite, so it stays in the V2 restart block.  See R-002 in REVIEW.md round 4 for the full diagnosis.

5. **`TSSetSolution` vs. PlaceMemory.** MFEM's `PetscODESolver::Run` calls `X->PlaceMemory(x.GetMemory(), true)` then `TSSolve(ts, X->x)`. The Vec `X->x` aliases the user's state Vector. After restart, we have a fresh `state` Vector (with the loaded values from `ReadCheckpoint`). The next `petsc_ode->Run(state, ...)` will PlaceMemory afresh, picking up the loaded values. So we do NOT need to call `TSSetSolution` explicitly — the PlaceMemory in Run does it.

6. **Phase 3 FSAL — needs a PETSc PR.** Document this clearly. The Phase 3 work is BLOCKED on either (a) upstream PETSc accepting `TSRKGetStageVectors` or (b) the MFEM-side workaround being accepted.

### Out of scope

- Restart across DIFFERENT MPI rank counts (re-partitioning the mesh + redistributing the state). The existing `ReadCheckpoint` rank-count check forbids this; the V2 extension does not relax it.
- Restart across DIFFERENT mesh files. The existing checkpoint format does not record the mesh; rest is the operator's responsibility.
- Binary (non-text) checkpoint format. The plain-text format is preserved for debuggability; binary would be a separate plan.
- Snapshot-file (`.vtkhdf`) append-mode for restart. The VTKHDF writer currently overwrites; on restart the fault.vtkhdf for the second window would clobber the first. A workaround for now: rename the previous `fault.vtkhdf` to `fault_part1.vtkhdf` manually before submitting the restart, OR set `--output-dir` to a NEW directory on restart. Both are documented in the sbatch comment.

- **TPV104 restart support — V1 ONLY (Phase 4 delivered 2026-05-17; V3 pending).**  REVOKED partial deferral: the user instructed "you must proceed with tpv104 as well" and a TPV104 V1 checkpoint port was delivered:
  - `io/tpv104_checkpoint.hpp` — `TPV104_CHECKPOINT_V1` format (Q + dof_data dynamic fields), MPIContext* + raw-MPI overloads.
  - `drivers/tpv104_driver.cpp` — `--restart` + `--checkpoint-interval` CLI flags, safety-check block, restart-load + write call sites.
  - `tests/unit/test_tpv104_checkpoint.cpp` — 6 implemented sub-tests covering V1 round-trip (#1), wrong-format/dof-size/Q-size guards (#4/#5/#7 — SKIP, MFEM_VERIFY non-catchable in-process; opt-in subprocess runtime check at #8 via `SEAS_TEST_RUNTIME_QSIZE_CHECK=1`), driver+header source-grep covering R-001..R-007 + R-101/R-102/R-104 (#6), and cross-overload byte-identity round-trip (#9, gates R-007 dedup).
  - `jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch` — dev queue, 8N×400r, two-phase verification.

  **What V1 carries:** time `t`, dt, step, bulk wave field `Q` (NUM_STATE × ndof_total), 9 dynamic fault-DOF fields (psi, slip_rate, V1, V2, slip1, slip2, tau1_nuc, tau2_nuc, sigma_n_nuc).
  **What V1 does NOT carry (deferred to V3):**
  - The SECONDARY `pv_bulk_out` collection's schedule state (regime, last_write_time).  V2 today carries the PRIMARY `pv_out` collection only.  TPV104 with `--paraview-bulk-dt > 0` AND regime-adaptive cadence on the bulk collection will silently restart that collection from default regime — the driver emits a `WARNING: --restart resumes the PRIMARY ParaView collection's schedule state ... V3 will carry both collections.` at restart-load time.
  - DOFData static fields (impedances, a, Dc, prestress, LSW params) — re-initialised by `InitializeFaultDOFs_TPV104`; do not need to round-trip.
  - DOFData corrected-traction fields (tau1_corr, tau2_corr, sigma_n_corr) — recomputed on first post-restart step.

  **V3 (future work, NOT in this plan):** extend the V2 trailing block to carry **per-collection** ParaView schedule state.  The V2 layout serializes one collection's 10 fields; V3 needs `n_collections` × {regime, last_write_time, last_v_max, last_committed_cycle, last_volume_write_time, snapshots, ...}.  When V3 lands, both the TPV104 driver warning and this Out-of-scope entry must be retracted.

---

## How to start implementing (for the next agent)

1. Read this plan top to bottom.
2. Confirm with `git log --oneline -5` that you're on `feature/paraview-compaction` at or after commit `407f456` ("BP5 Phase 6 ParaView: fix ResetMemory crash + cap rework + per-regime cadence").
3. Implement Phase 1 in a SINGLE atomic commit; do not bleed it into Phase 2 or 3.
4. Run `make seas_test_bp5_petsc_ts_restart` + `mpirun -np 2 ./seas_test_bp5_petsc_ts_restart` and confirm 11/11 PASS + 1 SKIP (sub-test 6, per R-011) before marking Phase 1 done.
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
