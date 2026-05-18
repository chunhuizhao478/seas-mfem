# Paraview-Compaction Merge — Round-2 Review Audit (2026-05-17)

## Context

This document tracks the **second round** of adversarial code review on the
`feature/safs-quasi-dynamic` branch after the `feature/paraview-compaction` merge
(`0ce73a4`) and the follow-up fix sweep (`1a805ef seas: address REVIEW.md R-001..R-009
findings`).

- **Branch:** `feature/safs-quasi-dynamic`
- **Merge commit:** `0ce73a4 Merge feature/paraview-compaction into feature/safs-quasi-dynamic`
- **Round-1 review:** previous `REVIEW.md` (R-001..R-009 — all addressed by `1a805ef`)
- **Round-2 review (this round):** `REVIEW.md` (R-201..R-207) — written by the
  /code-review pass invoked 2026-05-17.

The round-1 fixes are kept; this round audits the post-fix code as a *fresh* artifact
rather than verifying each round-1 finding.

## Merge Verification

| Check | Result | Evidence |
|-------|--------|----------|
| Conflict markers in working tree | None | `git status --porcelain` clean; recent commits show no `===` / `<<<` resolutions left over. |
| MFEM core changes applied | Applied | `mesh/vtkhdf.{cpp,hpp}` and `fem/datacollection.{cpp,hpp}` carry the HDFCompression selector + ZFP filter dispatch (Phase 2d.2 of the paraview-compaction plan). |
| MFEM rebuild required by merge | Done | macOS rebuild on `MFEM_USE_HDF5=YES + MFEM_USE_H5Z_ZFP=YES` succeeded; first cold pass surfaced 5 pre-existing latent issues (commit `6e45aaf seas: unblock full test sweep`). |
| Round-1 R-001..R-009 addressed | Yes | `1a805ef` commit message + diff cross-check; R-001 split-counter + SAFS guard, R-003 base default sizing, R-004 inline dedup removed, R-006 SKIP→`KNOWN_DISABLED_TESTS.md`, R-007 invariant-token banner check all applied. |
| Paraview tests present after merge | Yes | `seas_test_fault_surface_vtkhdf{,_mpi,_zfp}`, `seas_test_paraview_schedule_cap`, `seas_test_paraview_hdf_smoke`, `seas_test_paraview_rank0_warning_gate`, `seas_test_volume_hdf_compression`, `seas_test_vtkhdf_zfp` all built. |
| SAFS smoke tests present | Yes | `seas_test_safs_mode_wiring`, `seas_test_compute_safs_params`, `seas_test_safs_stress_config`, `seas_test_project_fault_prestress` all built. |
| Driver paraview wiring | TPV only | `tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp` wire `seas::ParaViewOutput`; `seas_driver.cpp` (BP5) does not (per CLAUDE.md §"Phase 4 deferred deviation (R-313)"). |

**Conclusion:** Merge is mechanically correct.  No conflict residue, no missing files,
no broken includes.  Round-1 fixes are present and behaviourally compatible.

## Round-2 Findings — Summary

| ID | Severity | Site | Category | Status |
|----|----------|------|----------|--------|
| R-201 | MODERATE | `drivers/tpv1{02,04,205}_driver.cpp` paraview_write else-branch | BUG (latent) | open — see /code-fix |
| R-202 | MODERATE | `trace/face_trace_logger.hpp:409-416` | BUG (latent) | open — see /code-fix |
| R-203 | MODERATE | `fault/rate_state_fault.hpp:214-223` | ASSUMPTION | open — see /code-fix |
| R-204 | LOW | `io/paraview_output.hpp:2253-2268` | QUALITY (performance) | open — partial /code-fix |
| R-205 | LOW | `fault/fault_geometry.hpp:840-845` | ASSUMPTION (maintainability) | open — see /code-fix |
| R-206 | LOW | `fault/rate_state_fault.hpp:156 / fault_geometry_safs.inl:58` | ASSUMPTION (contract) | open — see /code-fix |
| R-207 | LOW | `io/fault_vtu_binary.hpp:431-441` | ASSUMPTION (impl-defined hash) | open — see /code-fix |

Full description and suggested diffs in `/Users/chunhuizhao/projects/seas-mfem-safs/REVIEW.md`.

## SAFS-Specific Risk

The SAFS-mode wiring (Phase 6 §6 of `PLAN_onfaultstress.md`) is the primary deliverable
of `feature/safs-quasi-dynamic`.  The round-2 audit identified the following SAFS
hazards introduced or unaddressed by the merge:

1. **R-203 — t1-fallback sign-flip.**  The round-1 R-001 guard inside `SetSAFSMode`
   only rejects zero-normal fallbacks.  t1-degenerate fallbacks produce a basis whose
   sign convention is undefined per Tandem-FaultBasis semantics; a sidecar projection
   at such a DOF can end up *antiparallel* to the elastic traction, creating the
   positive-feedback regime CLAUDE.md flags (debug v8) as causing unbounded growth.
   The runtime warning at `fault_geometry.hpp:1000-1007` mentions only the zero-normal
   silent-zero hazard, not the t1 sign-flip risk.
2. **R-206 — `tau_pre_` caching order dependency.**  The `RateStateFaultOperator` BP5
   ctor copies `geom_->tau_pre_` into a member at construction time.
   `ComputeSAFSParams` overwrites the geom's `tau_pre_` with the sidecar projection.
   The canonical order (`ctor → ComputeSAFSParams → SetSAFSMode`) is correct; any
   reversal silently caches sidecar values into the BP5-mode read site.

Neither hazard breaks existing passing tests (`T_66_1..T_66_5`); both are latent
contracts that would surface on the first SAFS run that hits a real-mesh t1-fallback
DOF or a test that reorders construction.

## Verification Strategy for the /code-fix Pass

The `/code-fix` agent should consume `REVIEW.md` and apply the diffs in finding ID
order.  After each fix:

- **R-201:** Re-build and run `seas_test_paraview_schedule_cap`,
  `seas_test_paraview_hdf_smoke`, `seas_test_paraview_rank0_warning_gate`.  None
  currently exercise hysteresis>1, so they will pass — the fix is validated by
  inspection + the new `test_R201_fault_only_commit_uses_v_max.cpp` test.
- **R-202:** Add `test_R202_antiplane_tracer_refuses.cpp` (no existing test covers this
  composition).  Re-build all `seas_test_*` targets — should be a no-op for
  current code (no Antiplane+tracer composition is built).
- **R-203:** Re-run `seas_test_safs_mode_wiring`; the existing `T_66_5` positive
  control must keep passing on the healthy BP5 fixture (NumZeroNormalFallbacks ==
  NumT1Fallbacks == 0).  Add `T_66_6` negative-control (deferred — requires synthetic
  fault geometry helper not present in current test).
- **R-204:** No functional test; visual verification of `<output_dir>/FaultSurface/`
  layout after a small TPV104 dry-run with `--paraview-fault-vtu` and several writes.
  Atomicity tested by interrupting an actual run (manual, deferred).
- **R-205:** No directly observable test (the early-return path isn't exercised twice
  in any current code path).  Annotated comment suffices.
- **R-206:** Add `test_R206_construct_after_compute_safs_aborts.cpp` if test harness
  permits — currently SAFS tests in `test_compute_safs_params.cpp` use canonical order
  so no regression.
- **R-207:** Re-run `seas_test_fault_surface_vtkhdf_mpi` (cross-rank gather with field
  filter) on N=2 ranks; the hash check must still pass.  The FNV-1a replacement is
  deterministic, so the test outcome is unchanged.

## Risk Acceptance for Production SAFS Sweep

The round-2 audit does **not** find a critical bug that blocks SAFS development.  The
two SAFS-relevant findings (R-203, R-206) are LATENT — they will only fire on a future
run that either (a) builds a fault mesh whose `GetFaultDOFBasis` produces t1-degenerate
DOFs, or (b) reverses the construct-vs-ComputeSAFSParams order.

The current SAFS test fixture (`bp5_1000m.msh` per `T_66_5`) is reported to have 0
zero-normal and 0 t1-degenerate fallbacks, so the existing positive control passes.
Production SAFS meshes (the curvilinear CFM-derived meshes in
`miniapps/seas/safs/CFM_data_step/`) **have not been audited** against the t1-fallback
threshold in this round; the operator's `NumT1Fallbacks()` accessor should be logged
on the first ApplySAFSMode call to confirm 0 on production meshes BEFORE relying on
SAFS-mode results.

## Files Created / Modified by This Round

| File | Purpose |
|------|---------|
| `REVIEW.md` (root) | Round-2 findings (R-201..R-207) — full diffs and test cases |
| `miniapps/seas/debug_document/paraview_output_debug_document/round2_review_audit_2026-05-17.md` | This document — round-2 audit history |

No source files are modified by the review pass itself; the /code-fix invocation will
apply the diffs proposed in `REVIEW.md`.
