# Fix Report — REVIEW.md (PETSc TS restart plan, 2026-05-16)

**Date:** 2026-05-16
**Branch:** `feature/paraview-compaction`
**Review document:** `REVIEW.md` (7 findings: R-301..R-307)
**Target of fixes:** `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`

This is a **PLAN-FIX pass** — the edits land on a Markdown design document, not on source code. There is nothing to compile and nothing to run; "verification" means re-reading the plan to confirm the contradictions are gone, the cross-references to `bp5_verification_full.cpp` line numbers are accurate, and the V2 schema is self-consistent across its three appearances (writer doc-comment, writer body, reader signature). The plan is now ready for `/code-implement` (or a future engineer) to execute.

## Summary

- Findings addressed: **7 of 7** (R-301..R-307)
- Files modified: **1** (`petsc_ts_restart_plan_2026-05-16.md`)
- New tests: **0** (this is a plan-fix; tests are designed inside the plan itself for the future implementer)
- Test suite: **N/A** (no compilable code in this pass)
- Self-check: grep on the plan confirms the post-fix tokens (`"NO CHANGES"`, `"AFTER the existing V1 restart block"`, `"restart_rejections_carryover"`, `"RestoreScheduleState"`, `"T_full and T_mid"`, `"std::max(n, 0)"`, `"TSRKSetStageVectors"`) all appear in the expected sections.

## Changes Made

### R-301 [CRITICAL] — Lock plan on Option B (V1 header unchanged, append V2 trailer)

**File:** plan line 69 ("Files to Modify → checkpoint.hpp")

Removed the contradiction with line 154 ("Format version detection") by deleting the "bump V1→V2" instruction and replacing it with an explicit "**NO CHANGES** to the file-format tag" directive. Added a forward reference to the "Format version detection" section so a reader hitting line 69 first knows where the canonical statement lives.

```diff
-- `miniapps/seas/io/checkpoint.hpp` — bump the file-format tag from `SEAS_CHECKPOINT_V1` to `SEAS_CHECKPOINT_V2`. The reader recognises both tags: V1 is read as today (no TS-internal fields); V2 has an additional trailing block read by `ReadPetscTSCheckpoint`. The writer always writes V2.
+- `miniapps/seas/io/checkpoint.hpp` — **NO CHANGES** to the file-format tag or to the existing `WriteCheckpoint` / `ReadCheckpoint` signatures.  The V2 extension is appended AFTER the V1 block by the new `WritePetscTSCheckpoint` function...  See §"Format version detection" below for the canonical statement.
```

### R-302 [CRITICAL] — Move V2 restart block to after V1 ReadCheckpoint

**File:** plan line 73 (summary bullet) + the V2 code listing's preamble at line ~169

Re-specified placement at "IMMEDIATELY AFTER the existing V1 restart block at lines 2346-2358" instead of "after `petsc_ode->Init(...)` around line 2269". Added an explicit anti-instruction: "Do NOT place the V2 block inside the PetscTS init at line ~2269 — `t` is still 0 there, the cross-check would always fail." The code listing's surrounding comment was updated to match.

```diff
-  - **Add** in the PetscTS init block (after `petsc_ode->Init(...)` around line 2269): if `!restart_prefix.empty()`, call `ReadPetscTSCheckpoint(...)`...
+  - **Add** the V2 restart block IMMEDIATELY AFTER the existing V1 restart block at lines 2346-2358.  **Ordering is load-bearing**: the V1 `ReadCheckpoint` call at line 2352 populates `t`, `current_dt`, `state`...  The cross-check `std::abs(t - ts_t) < 1e-12 * std::abs(t)` is only meaningful BECAUSE `t` has been populated by the preceding V1 read.  **Do NOT place the V2 block inside the PetscTS init at line ~2269** — `t` is still 0 there...
```

### R-303 [MODERATE] — Accumulate `step_rejections` across restart

**File:** plan V2 code listing at lines ~199-242 + new diff snippet at lines ~257-269

Replaced the dead-store assignment `step_rejections = ts_rejections` with `restart_rejections_carryover = ts_rejections;` (a new local variable preserved across the seam). Added a diff snippet showing the post-Run accumulation: `step_rejections = restart_rejections_carryover + static_cast<int>(rejects);`. Added an instruction to declare the carryover variable next to the existing `int step_rejections = 0;` at `bp5_verification_full.cpp:2223`.

The rationale ("TSGetStepRejections returns THIS-Run's rejections only; it isn't reset by TSSetStepNumber but it ISN'T pre-populated from the checkpoint either") is in the inline comment of the diff snippet.

### R-304 [MODERATE] — Extend V2 schema with ParaView schedule state

**Files:** plan multiple sections:
1. `WritePetscTSCheckpoint` doc-comment table extended with `paraview_last_write_time`, `paraview_last_v_max`, `paraview_current_regime`.
2. `WritePetscTSCheckpoint` signature extended with the three new params (in order: after `paraview_snapshots`, before `mpi`).
3. Writer body extended with three new `out << ...` lines for the new fields.
4. `ReadPetscTSCheckpoint` signature extended with three matching out-params.
5. V2 restart block listing extended with `real_t ts_last_write_time = -1e30;` etc. declarations and the matching `ReadPetscTSCheckpoint` call arguments.
6. V2 restart block calls `pv_out->RestoreScheduleState(ts_last_write_time, ts_last_v_max, ts_current_regime);` after `SetTotalSnapshotsWritten`.
7. Write call-site listing extended with three new accessor reads.
8. `paraview_output.hpp` modifications section adds `RestoreScheduleState` setter spec + three read-only accessors (`GetLastWriteTime`, `GetLastVMax`, `GetCurrentRegime`) with full body listings.
9. "Files to Modify → paraview_output.hpp" summary lists `RestoreScheduleState` and the three accessors as separate bullets.
10. "Interfaces → New exposed functions" lists the new methods.

The rationale is documented inline: "Without this restoration the first ShouldWrite after restart fires unconditionally (because `last_write_time_` defaults to -1e30 and `time - (-1e30)` always exceeds dt_out * tol) and the regime state machine resets to interseismic regardless of where the pre-checkpoint trajectory was."

### R-305 [MODERATE] — Rewrite Sub-test 3 in terms of T_full / T_mid (not step counts)

**File:** plan Sub-test 3 procedure at lines ~410-460 (now ~530-590 after the other edits)

Rewrote the procedure to:
1. Pick `T_mid` and `T_full` upfront (with concrete suggested values: `dt_init = 0.01 s, T_mid = 1.0 s, T_full = 2.0 s`).
2. Use `petsc_ode->Run(state, t, dt, T_full)` as the stopping primitive (deterministic because PETSc/MFEM set `TS_EXACTFINALTIME_MATCHSTEP` at `linalg/petsc.cpp:4365`).
3. State explicitly that the `t_A == T_full` and `t_B == T_full` assertions are bit-exact under FP (no `|t_A - t_B| < 1e-12` slop).
4. Add a step (f.) reminding the implementer to apply the R-003 padded-shrink BEFORE `ReadCheckpoint`, then `SetTotalSnapshotsWritten` + `RestoreScheduleState` after.
5. Add a preamble paragraph explaining WHY step counts can't be used as a stopping primitive ("controller dt is dependent on local error magnitudes").

### R-306 [LOW] — Align `SetTotalSnapshotsWritten` summary with detailed body

**File:** plan line 78

Updated the summary bullet to reference `std::max(n, 0)` explicitly (matching the detailed body at line 320-326), and added a one-line rationale ("Negative values are clamped to 0 to defend against corrupted/truncated V2 checkpoints"). Forward-references the full body in §"Detailed Requirements" item 6.

### R-307 [LOW] — Add `TSRKSetStageVectors` to Phase 3 upstream PR list

**Files:** plan Phase 3 "Files to Modify" + Phase 3 "Detailed Requirements"

Updated the upstream PR description to list BOTH:
* `PetscErrorCode TSRKGetStageVectors(TS ts, PetscInt *nstages, Vec **Y)` — for the checkpoint WRITE side.
* `PetscErrorCode TSRKSetStageVectors(TS ts, PetscInt nstages, Vec *Y)` — for the checkpoint READ side.

Same update applied to the MFEM-only alternative path ("Add BOTH `PetscODESolver::GetRKStageVectors` AND `PetscODESolver::SetRKStageVectors`"). Updated the "Detailed Requirements" READ-side step (item 2) to reference `TSRKSetStageVectors` as the matching API rather than ambiguous "(or equivalent)" prose.

## Unresolved Findings

None — all 7 findings addressed.

## Deviations from the Review's Suggested Fixes

None substantive. The fixes follow the review's diff blocks verbatim where the review provided them; where the review only described the intent (R-304, which spans 10 separate sections of the plan), I extended every dependent section to keep the spec internally consistent.

## Self-Check

Grep over the final plan confirms:

- `grep -c "NO CHANGES"` → 1 (R-301)
- `grep -c "AFTER the existing V1 restart block"` → 1 (R-302)
- `grep -c "restart_rejections_carryover"` → 4 (R-303 — declaration mention + assignment + diff snippet + post-Run usage)
- `grep -c "RestoreScheduleState"` → 6 (R-304 — summary bullet, V2 listing call, detailed body, Interfaces, Sub-test 3, and one cross-reference)
- `grep -c "T_full"` → 7 (R-305 — Sub-test 3 procedure mentions T_full at multiple steps)
- `grep -c "std::max(n, 0)"` → 2 (R-306 — summary bullet + detailed body)
- `grep -c "TSRKSetStageVectors"` → 3 (R-307 — Phase 3 PR list + MFEM alternative + READ-side detailed requirement)
- `grep -c "around line 2269"` → 1 (R-302 — appears only in the anti-instruction warning, never as the actual placement)

All seven post-fix tokens land in the expected sections.

## Notes for Reviewer Re-Review

- The plan now references three new public accessors on `ParaViewOutput` (`GetLastWriteTime`, `GetLastVMax`, `GetCurrentRegime`). These are read-only and have no impact on non-restart callers, but a re-reviewer should sanity-check that exposing the previously-private `last_write_time_`, `last_v_max_`, `current_regime_` to a `const` getter is consistent with project conventions. The existing `GetTotalSnapshotsWritten()` at `paraview_output.hpp:368` sets the precedent.

- The V2 schema is now 8 fields (5 PetscTS + 3 ParaView). If a future extension adds, say, an RK adaptive-controller PI-history field (Phase 2 documents this as currently empty), the schema bumps to V3 — the plan's "Format version detection" section already anticipates this with its tag-based extensibility model.

- The plan still contains an "Out of scope" item about VTKHDF append-mode on restart (snapshot files get overwritten if `--output-dir` is reused). The R-201..R-307 fix pass did NOT add this to the schema; it remains a documented operator-procedure workaround ("rename the previous `fault.vtkhdf` manually before re-submit, or set a new `--output-dir`"). A future Phase 4 could automate this, but is out of scope here.

## Ready for Re-Review: YES
