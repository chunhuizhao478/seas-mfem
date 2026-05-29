# Code Review: Phase 10 — TPV31 (round 2, post-fix fresh adversarial pass) — 2026-05-29

## Review Scope
- Plan: `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md` §Phase 10
- This is a **fresh round-2 review** after the /code-fix pass that addressed round-1 R-001…R-006.
  All three passes were re-run from scratch on the changed files (not just a fix-checklist):
  `drivers/spatial_dyn_driver.cpp`, `dynamic/nucleation_factory.{hpp,cpp}`,
  `dynamic/spatial_nucleation.{hpp,cpp}`, `dynamic/tpv205_substep_iterator.cpp`,
  `spatial/code/spatial_friction.cpp`, `tpv31/configs/tpv31.toml`, and the two TPV31 unit tests.
- Verification performed this round (`conda activate mfem-dev`, worktree MFEM_DIR override):
  - `seas_spatial_dyn_driver` **builds clean** with R-001/R-002/R-003.
  - `seas_test_tpv31_canonical_rotation` **52/52** (incl. new T-8 sign regression).
  - `seas_test_tpv31_nucleation_and_cohesion` **23/23** (incl. new I-6 µ-scaling, C-1z axis-z, C-1b region).
  - `seas_test_nucleation_factory` **22/22** (signature change is backward-compatible).
  - `seas_test_tpv31_station_writer` **36/36**, `seas_test_heterogeneous_material` **all pass**.

## Round-1 findings — verification
- **[R-001] FIXED & VERIFIED.** `spatial_dyn_driver.cpp:1424-1437` now passes `-dp.sigma_xy_per_mu`
  (with a comment tying it to the `constant_tensor` D3.1 negation at `:1345`).  New test T-8
  (`test_tpv31_canonical_rotation.cpp`) projects the source through the no-flip rule and asserts
  `tau_strike = +30 MPa` (right-lateral), and asserts the **un-negated** value yields a negative
  (left-lateral) traction — so a regression that removes the negation fails loudly. ✅
- **[R-002] FIXED & VERIFIED.** `InstantaneousOverstressCircularSpec.mu_ref_pa` added + parsed
  (`spatial_friction.cpp:1322`), resolver scales by `mu(point)/mu_ref` when a non-empty
  `mu_at_xyz` + `mu_ref_pa>0` are supplied (`spatial_nucleation.cpp`), driver builds the µ
  callback and threads it (`spatial_dyn_driver.cpp:1505-1523`), `tpv31.toml` sets
  `mu_ref_pa = 32.03812032e9`.  New test I-6 verifies scale 1 / scale 2 / disabled-when-empty /
  disabled-when-mu_ref=0. The default (0) keeps every non-TPV31 path byte-identical. ✅
- **[R-003] FIXED & VERIFIED (build).** `TPV31StationWriter` is now constructed + `Open`ed
  (MPI overload) before the time loop, `WriteStep` called every step, `Close` after — gated on
  `cfg.problem.tag == "tpv31"` (`spatial_dyn_driver.cpp`).  Driver compiles and links.  Produces
  the `_station_*.dat` files `tpv31/visualize_results.py` consumes.  (End-to-end station output
  not run locally — needs the 50 m mesh, which is not built; per CLAUDE.md "no local full-mesh
  runs".) ✅ build-verified; ⏳ run-verification deferred to the cluster smoke (Phase 10 AC #2).
- **[R-004] FIXED & VERIFIED.** T-8 (projection sign), I-6 (µ-scaling), C-1z (production axis-z
  depth=max(0,-z) + above-surface cap clamp), C-1b (z-bounds region restriction, distinguishable
  from the floor clamp via a distinct default).  All pass. ✅
- **[R-005] FIXED.** `tpv205_substep_iterator.cpp` now documents that the standalone oracle is
  cohesion-unaware by design (C0≡0); comment-only, no behavior change. ✅
- **[R-006] FIXED.** `test_tpv31_nucleation_and_cohesion.cpp` C-1 comment corrected (the y=5000
  DOF is zeroed by the floor clamp, NOT by a y-bound the Depth-kind rule does not apply); the real
  region-restriction case is now C-1b. ✅

## New Findings (introduced or surfaced by the fixes)

### [R-007] [LOW] [spatial_dyn_driver.cpp — TPV31StationWriter::Open on restart] — restart truncates the SCEC station traces

**Category:** EDGE_CASE

**Description:** `TPV31StationWriter::Open` opens each `.dat` with `std::ofstream::open(fname)`
(truncate). The driver calls `Open` unconditionally for `tpv31` runs, including on restart
(`restart_prefix` non-empty); only the *initial* `WriteStep(t_initial,…)` is gated off on restart.
So a restarted TPV31 run **erases the pre-restart station history** and resumes the trace from
`t_restart + dt` (the `t_restart` instant itself is also skipped). ParaView/HDF5 fault output is
unaffected; only the per-station `.dat` traces lose their pre-restart rows.

**Trigger:** `seas_spatial_dyn_driver tpv31.toml --restart <cp>` after a checkpoint.

**Actual:** station `.dat` files truncated and restarted mid-run; `visualize_results.py` then plots
only the post-restart segment.

**Expected:** continuation appends to the existing traces (or the writer is told it is a restart).

**Suggested fix:** lowest-risk option — only open the writer on a fresh run, and on restart either
(a) open in append mode, or (b) document that TPV31 validation runs must complete without restart.
Minimal code change (keeps the hrs-ref-identical writer header untouched):
```cpp
   // Open (and seed t_initial) only on a fresh run; restart would truncate the
   // pre-restart traces.  TODO(append-on-restart): add an append-mode Open if a
   // checkpointed TPV31 validation run is needed.
   if (cfg.problem.tag == "tpv31" && restart_prefix.empty())
   {
      tpv31_stations = std::make_unique<TPV31StationWriter>();
      ... Open ...; tpv31_stations->WriteStep(cfg.time.t_initial, dof_data);
   }
```
This is safe because TPV31 SCEC validation is a single 15 s run (checkpoints are a crash-safety
net, not a workflow step). NOT fixed in this pass to avoid an API change / behavior choice without
the user's call; flagged for decision. Severity LOW (no wrong physics; affects only restart trace
continuity).

**Test case:** integration-level (run a few steps, checkpoint, restart, assert the `.dat` contains
rows with `t < t_restart`) — deferred with the cluster smoke; not unit-testable cheaply.

---

### [R-008] [LOW] [POSSIBLE] [spatial_dyn_driver.cpp:1505-1521 — nuc µ callback for non-depth-profile heterogeneous material] — instantaneous µ-scaling silently disabled for a sidecar material

**Category:** ASSUMPTION

**Description:** The driver builds `nuc_mu_at_xyz` only for `Mode::Constant` (returns `mu_const`) or
a `depth_profile_wrapper` (Coefficient). For a **sidecar (CVM) heterogeneous** material the
callback stays empty, so an instantaneous-circular nucleation with `mu_ref_pa>0` would silently
fall back to uniform amplitude (no µ-scaling). This is consistent with the depth-proportional
*stress* path, which instead **aborts** for sidecar (`spatial_dyn_driver.cpp:1418`). TPV31 uses
`depth_profile_1d`, so this never triggers for the only config that sets `mu_ref_pa>0`; hence
POSSIBLE/LOW. If a future sidecar config combined instantaneous nucleation with µ-scaling it would
get an unscaled patch with no warning.

**Suggested fix (defer):** when `cfg.nucleation.kind == InstantaneousOverstressCircular` AND
`instantaneous_circular.mu_ref_pa > 0` AND no coordinate-mu lookup is available, abort with a clear
message (mirror the depth-proportional-stress abort), rather than silently dropping the scaling.
Not fixed now (no such config exists; would add an abort path with no test coverage).

---

## Byte-exactness / regression safety (verified this round)
- **Non-TPV31 configs unchanged.** R-001 touches only `stress.kind="depth_proportional"` (TPV31
  only). R-002's callback is forwarded only to the instantaneous resolver and only scales when
  `mu_ref_pa>0` (default 0); Gaussian/compact-circular/static ignore it. R-003 is gated on the
  TPV31 problem tag. So TPV205/102/104 and SAFS paths are byte-identical (confirmed:
  `nucleation_factory` 22/22 unchanged, the 3-arg `MakeNucleation` call sites still compile/pass).
- **No dangling reference** in the µ callback: `nuc_mu_at_xyz` captures `eval` by reference but is
  consumed synchronously inside `MakeNucleation` (the resolver evaluates it per-DOF and returns
  amplitudes); `depth_profile_wrapper` outlives the call.

## Summary
- Critical issues: 0
- Moderate issues: 0
- Low issues: 2 (R-007 restart truncation; R-008 sidecar µ-scaling edge — both deferred-by-design)
- Plan compliance: **FULL** for the reviewable+buildable scope — material (req 1), depth-proportional
  stress sign (req 2/3), cohesion (req 3), µ-scaled nucleation (req 2), and station output wiring
  (req 5) are all in place and unit-verified.
- Verdict: **PASS.** Round-1 R-001 (CRITICAL) is fixed and locked by a regression test; R-002/R-003
  (MODERATE) are implemented and (for R-002) unit-verified. The two new LOW findings are edge cases
  that do not affect a standard single-run TPV31 validation and are flagged for a product decision.

## Unreviewed Areas (unchanged from round 1)
- End-to-end TPV31 run (needs the 50 m mesh; not built locally — only `tpv31_50m.geo` is present —
  and "no local full-mesh runs" per CLAUDE.md). Phase 10 AC #1 (`--verify-dispatch`) and AC #2
  (8-rank smoke writing station output that matches the SCEC overlays) are cluster-side.
- Phase 9/13 bimaterial wave operator, PML, RK stepper — out of Phase-10 scope.
