# Code Review: Phase 5 (QD time loop) + Phase 6 (checkpoint/restart) — 2026-06-04

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` §Phase 5, §Phase 6.
- Files reviewed:
  - Phase 5: `drivers/spatial_seas_driver.cpp` (friction → fault-op → §B `ComputeParams` dispatch → SAFS → RK45 loop), `spatial/code/spatial_friction.{hpp,cpp}` (bp5_analytic flag), `tests/unit/test_spatial_seas_qd_coupling_smoke.cpp`, `config/safs_qd/bp5_phase5_50step.toml`, `jobs/safs_qd/spatial_seas_bp5_50step.sbatch`.
  - Phase 6: `io/seas_qd_checkpoint.hpp`, `tests/unit/test_seas_qd_checkpoint.cpp`, `drivers/spatial_seas_driver.cpp` (restart/checkpoint wiring + restart-safety).
- Domain context: CLAUDE.md (4-phase init, RK45 SetDt/Step, sign rules, extreme-care files), the plan (R-001/002/003/005/009), project memory ([[project_qd_safs_coupling_needs_clean_fault_basis]], [[no-local-reproducer]]), the implementer completion reports.
- Method: three adversarial passes; verified the RK45 `Step` dt-out vs `GetDt()` semantics (`time_stepper.hpp:230` + accept path) and the `geom` accessor set directly; build + the two new unit tests run.

## Findings

### [R-601] MODERATE [BUG] drivers/spatial_seas_driver.cpp:(RK45 loop) — mid-loop checkpoint persists the wrong dt (breaks restart reproduction)

**Category:** BUG

**Description:**
`DormandPrinceRK45::Step` sets its OUT `dt = dt_` at entry (the trial = the step actually taken) and, on accept, updates the internal `dt_` to the NEXT trial (`time_stepper.hpp:230` + the accept/grow path). So after a step, `dt_taken` is the step just taken while `ode.GetDt()` is the next trial. For a restart to reproduce the un-checkpointed trajectory, the RK45 must be seeded with the NEXT trial dt (what the un-restarted run uses for its next step). The **mid-loop** checkpoint writes `dt_taken`:
```
if (ckpt_every > 0 && step % ckpt_every == 0) { write_checkpoint(t, dt_taken, step); }
```
while the **final** checkpoint correctly writes `ode.GetDt()`. On restart, `ode.SetDt(ck.dt)` is then seeded with `dt_taken` (the already-taken step) instead of the next trial, so the first post-restart step uses a different dt than the un-restarted run's corresponding step → the trajectories diverge at the RK45-tolerance level (~1e-7), NOT to round-off. This breaks the Phase-6 acceptance ("restart reproduces the trajectory to round-off for ≥10 steps").

**Trigger:** restart from a mid-loop checkpoint (any `--checkpoint-every N` where the run is checkpointed before the final step).

**Actual behavior:** restart seeds dt = the accepted step of the checkpointed step → diverges from the un-restarted run.

**Expected behavior:** restart seeds dt = the next trial (`ode.GetDt()` after the checkpointed step) → bit-reproduces.

**Suggested fix:** use `ode.GetDt()` in the mid-loop checkpoint too (match the final checkpoint):
```diff
-         if (ckpt_every > 0 && step % ckpt_every == 0)
-         {
-            write_checkpoint(t, dt_taken, step);
-         }
+         if (ckpt_every > 0 && step % ckpt_every == 0)
+         {
+            write_checkpoint(t, ode.GetDt(), step);
+         }
```

**Test case:**
```
// test_R601_restart_reproduces_dt (FRONTERA, BP5 mesh — needs SAFS coupling):
//   run A: --checkpoint-every 5 --max-steps 20  → checkpoint at step 5.
//   run B: --restart <A's prefix> --max-steps 20 (15 more steps).
//   assert run B's state at each step == run A's state to round-off (1e-14),
//   NOT just to RK45 tolerance.  Pre-fix: diverges at ~1e-7 after the restart.
```

---

### [R-602] MODERATE [POSSIBLE] drivers/spatial_seas_driver.cpp:(restart branch) — persisted static params are never validated on restart (config drift undetected)

**Category:** ASSUMPTION / EDGE_CASE

**Description:**
The checkpoint persists per-DOF static params (a, Dc, eta, V_init, σ_n, τ_pre), but the restart path reads them into `ck` and uses ONLY `ck.state`/`ck.t`/`ck.step`/`ck.dt`. The driver re-derives geom (re-runs `SetRateStatePerDOF` + `ComputeParams`) before the restart branch, so the persisted params are dead on read AND a config change between the original run and the restart (e.g., a different `a_default`) is silently applied: the restart continues the OLD evolved state under the NEW friction params → a physically wrong continuation with no warning. The plan persists these params precisely so a restart can validate against config drift.

**Trigger:** restart with a `[friction.rate_state]` or `[stress]` value changed from the original run.

**Actual behavior:** new params + old state, silently.

**Expected behavior:** abort if the re-derived geom params disagree with the checkpoint's.

**Suggested fix:** after the restart restore (geom is already re-derived), compare and abort on mismatch:
```cpp
auto check_param = [&](const char *nm, const Vector &re, const Vector &ck_v){
   MFEM_VERIFY(re.Size() == ck_v.Size(), "restart: " << nm << " size drift");
   real_t e = 0.0;
   for (int i = 0; i < re.Size(); i++) e = std::max(e, std::abs(re(i)-ck_v(i)));
   MFEM_VERIFY(e <= 1e-12*(1.0+std::abs(re.Normlinf())),
               "restart: " << nm << " differs from checkpoint (config drift) max="<<e);
};
check_param("a", geom.GetAValues(), ck.a);
check_param("Dc", geom.GetDcValues(), ck.Dc);
check_param("eta", geom.GetEtaValues(), ck.eta);
check_param("V_init", geom.GetVInit(), ck.V_init);
check_param("sigma_n", geom.sigma_n_per_dof(), ck.sigma_n);
check_param("tau_pre", geom.GetTauPre(), ck.tau_pre);
```

**Test case:**
```
// test_R602_restart_rejects_config_drift (can be a unit test on the header +
//   a stub geom): write a checkpoint with a=0.010; "re-derive" a=0.012;
//   the validation must abort.  (Local: exercise check_param directly.)
```

---

### [R-603] MODERATE [POSSIBLE] drivers/spatial_seas_driver.cpp:(§B dispatch) — no σ_n consistency guard (R-003)

**Category:** ASSUMPTION

**Description:**
The plan (R-003) requires the stress-source σ_n to equal the TOTAL σ_n the resolver used in Phase 3. The §B comment documents this but nothing enforces it. For `FaultLocalPrestress`, `[stress].sigma_n_pa` and `[friction.rate_state].sigma_n_default` are independent config values; a mismatch is only caught later as a `SetInitialCondition` equilibrium failure (`> 1e-6`), with a confusing message that points at equilibrium rather than the σ_n mismatch.

**Trigger:** a config where `[stress].sigma_n_pa != [friction.rate_state].sigma_n_default` (FaultLocalPrestress path).

**Suggested fix:** add an explicit guard in the `FaultLocalPrestress` branch:
```cpp
MFEM_VERIFY(std::abs(cfg.stress.sigma_n_pa - cfg.rate_state->sigma_n_default)
            <= 1e-9 * (std::abs(cfg.stress.sigma_n_pa) + 1.0),
            "spatial_seas: [stress].sigma_n_pa (" << cfg.stress.sigma_n_pa
            << ") must equal [friction.rate_state].sigma_n_default ("
            << cfg.rate_state->sigma_n_default << ") for a consistent "
            "equilibrium (plan R-003).");
```
(Verify the exact field path for the rate-state σ_n default in `RateStateConfig`.)

**Test case:**
```
// test_R603_sigma_n_mismatch_rejected: a config with sigma_n_pa=50e6 and
//   sigma_n_default=40e6 must abort at the dispatch with the σ_n message,
//   not proceed to an equilibrium failure.  (Needs the driver / a config stub.)
```

---

### [R-604] LOW drivers/spatial_seas_driver.cpp:(bp5_analytic branch) — fault basis col 0 used for all DOFs without a planarity assert

**Category:** ASSUMPTION

**Description:**
The bp5_analytic branch reads `n/t1/t2` from column 0 of `geom.fault_dof_basis()` and passes that single basis to `Bp5AnalyticStressSource`, which applies it uniformly to all DOFs. This is correct ONLY for a planar fault (all columns equal up to a global sign — verified for BP5 in `test_spatial_seas_bp5_analytic`). On a non-planar fault it would silently use the wrong basis for every non-col-0 DOF. There is no assert that the basis columns are uniform.

**Suggested fix:** assert col-uniformity (up to a global sign) before using col 0, or document that bp5_analytic requires a planar fault:
```cpp
// (after extracting col 0) sanity-check a few columns share the normal direction:
for (int j = 1; j < std::min(geom.NumFaultDOFs(), 8); ++j) {
   const real_t dot = Bsis(0,0)*Bsis(0,j) + Bsis(1,0)*Bsis(1,j) + Bsis(2,0)*Bsis(2,j);
   MFEM_VERIFY(std::abs(std::abs(dot) - 1.0) < 1e-9,
               "spatial_seas: [stress].bp5_analytic requires a planar fault "
               "(column "<<j<<" normal differs from column 0).");
}
```

---

### [R-605] LOW drivers/spatial_seas_driver.cpp:write_checkpoint — aborts if output_dir does not exist

**Category:** EDGE_CASE

**Description:**
`ckpt_prefix = cfg.output.output_dir + "/spatial_seas"`; if `output_dir` does not exist, the first `WriteSpatialSeasCheckpoint` hits `MFEM_VERIFY(out.good())` and aborts. The driver does not create `output_dir`.

**Suggested fix:** create the directory once on rank 0 before the loop (the dynamic driver does this), e.g. `std::filesystem::create_directories(cfg.output.output_dir)` guarded by `rank == 0` + a barrier, or document that `--output-dir` must pre-exist.

---

### [R-606] LOW spatial/code/spatial_friction.cpp + driver — `[stress].bp5_analytic=true` silently overrides `[stress].kind`

**Category:** QUALITY

**Description:**
When `bp5_analytic=true`, the dispatch's first `if` wins and `[stress].kind` is ignored. A config setting both `kind="depth_proportional"` and `bp5_analytic=true` silently runs bp5_analytic. Harmless but surprising.

**Suggested fix:** warn once on rank 0 if `bp5_analytic` is set together with a non-default `kind`, or document the precedence in the StressSpec comment (the comment already says "reaches … WITHOUT extending kind" — add "and takes precedence over [stress].kind").

---

## Summary
- Critical issues: **0**.
- Moderate issues: **3** — R-601 (mid-loop checkpoint dt breaks restart reproduction — the strongest find), R-602 (restart doesn't validate persisted params → config drift), R-603 (no σ_n consistency guard).
- Low issues: **3** — R-604 (bp5_analytic planarity), R-605 (output_dir not created), R-606 (bp5_analytic precedence).
- Plan compliance: **PARTIAL** — Phase 5 wiring + Phase 6 checkpoint are implemented and (locally) verified; the SAFS init/stepping + restart-reproduces runs are Frontera (documented). Phase 6 stations/probe + ParaView (fault VTU) are NOT yet wired (implementer-documented remaining work; Phase 7 needs station output).
- Verdict: **PASS WITH FIXES** — R-601 must be fixed before the Frontera restart-reproduction test can pass; R-602/R-603 are cheap robustness guards; R-604/605/606 are LOW. No blocker for the local build/tests (both new unit tests pass), but R-601 blocks the Phase-6 checkpoint acceptance.

## Unreviewed Areas
- Phase 6 stations/probe + ParaView (fault VTU): not yet implemented (out of scope this round).
- Behavior on the BP5 mesh (SAFS coupling init/stepping, restart reproduction): Frontera — not locally runnable ([[project_qd_safs_coupling_needs_clean_fault_basis]]).
- The Frontera config/sbatch (`bp5_phase5_50step.toml` / `spatial_seas_bp5_50step.sbatch`): schema-matched to working configs but not executed (production mesh).
