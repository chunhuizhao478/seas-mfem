# Code Review: [DIAG-SIGN] speckle instrumentation (LSW triq-mesh test) — 2026-05-26

> Supersedes the prior σ_n-strength-floor review (recoverable via git history).
> Scope is the speckle diagnostic added in commit 5909cfe.

## Review Scope
- Plan: no formal PLAN.md; spec is the conversational debug task + commit
  `5909cfe` ("instrument the diagnostics ... capture the bug location") and
  `safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_sliver_blowup_2026-05-25.md`.
- Files reviewed (commit 5909cfe):
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` ([DIAG-SIGN] block + reset)
  - `miniapps/seas/dynamic/fault_face_flux.hpp` (DOFData.sigma_n_substep_min)
  - `miniapps/seas/dynamic/tpv205_substep_iterator.cpp` (LSW sub-step update)
  - `miniapps/seas/dynamic/wave_operator.inl` (shared-fault sub-step update)
  - `miniapps/seas/jobs/safs/spatial_dyn_slipweakening_nocap_triq_8N_400r_dev_2hr_safs.sbatch`
- Domain context: `miniapps/seas/CLAUDE.md` (σ_n>0 compression; no hardcoded
  constants), memory `project_safs_vn_leak_not_frame.md` (speckle locus =
  per-sub-step predictor/ghost `[[v_n]]` leak, R-1303/R-1601).

## Verified correct (not findings)
- Scale consistency: `data.sigma_n_corr` (`fault_face_flux.cpp:322`) and
  `s.sigma_n_total` (`:91`/`:182`) are both `σ_n0 + σ_n_nuc + trial` → the
  end-of-step vs sub-step comparison is like-for-like.
- `sigma_n_total` IS populated on both the LSW iterator and the shared-fault
  RS-style path (`fault_face_flux.cpp:182`), so shared QPs do not feed a
  spurious `0` into the `min` (first hypothesis tested — refuted).
- Sign: tracking `min` σ_n (most negative = most tensile) and `< floor` is the
  correct engagement test; nocap (`floor<0`) correctly reports DISABLED and
  skips `n_below_floor` while still counting `n_tensile`.

## Findings

### [R-001] [MODERATE] spatial_dyn_driver.cpp:[DIAG-SIGN] — only ONE speckle spot localized per step

**Category:** BUG (capture-completeness)

**Description:**
The screenshot shows 3–4 simultaneous speckle clusters, but `[DIAG-SIGN-DOF]`
reports only the single GLOBAL most-tensile DOF (`MPI_MINLOC` → one rank, one
`argmin_ss_local`). Other concurrent spots are invisible (counted in
`n_tensile`, never localized). The stated goal is to "capture the bug
location"; with multiple spots this captures at most one location per step.

**Trigger:** ≥2 fault DOFs tensile in the same macro-step (the observed case).

**Actual behavior:** one `[DIAG-SIGN-DOF]` line (global worst) per diag step.

**Expected behavior:** localize all (or top-K) tensile DOFs so every speckle
cluster's coordinates are recorded.

**Suggested fix:** every rank dumps its local tensile DOFs (capped), instead of
only the global MINLOC winner.
```diff
-         if (rank == ss_out.r && argmin_ss_local >= 0)
-         {
-            const DOFData &d = dof_data[argmin_ss_local];
-            ... single global-worst dump ...
-         }
+         // Dump EVERY local DOF whose sub-step σ_n is tensile/below floor,
+         // capped, so all concurrent speckle spots are localized.
+         {
+            const real_t sign_thr = (sn_floor >= 0.0) ? sn_floor : 0.0;
+            int dumped = 0;
+            for (int i = 0; i < num_fault_total && dumped < 32; ++i)
+            {
+               const DOFData &d = dof_data[i];
+               if (d.sigma_n_substep_min < 1.0e299 &&
+                   d.sigma_n_substep_min < sign_thr)
+               {
+                  std::cout << "[DIAG-SIGN-DOF] rank " << rank << " dof " << i
+                            << " xyz=(" << dof_coords_3d(3*i) << ","
+                            << dof_coords_3d(3*i+1) << ","
+                            << dof_coords_3d(3*i+2) << ")"
+                            << " sigma_n_substep_min=" << d.sigma_n_substep_min
+                            << " sigma_n_corr(end)=" << d.sigma_n_corr
+                            << " V_substep_max=" << d.slip_rate_substep_max
+                            << "\n";
+                  ++dumped;
+               }
+            }
+         }
```

**Test case:**
```python
def test_R001_all_tensile_spots_localized():
    # 2 engineered tensile DOFs (or parse a known multi-spot log).
    lines = [l for l in run_log
             if l.startswith("[DIAG-SIGN-DOF]") and at_step(l, "t=3.")]
    assert len({xyz_of(l) for l in lines}) >= n_known_tensile_spots
```

---

### [R-002] [MODERATE] sbatch + [DIAG-SIGN] — captures WHERE σ_n is tensile but not WHY (opening vs bulk); `SEAS_DIAG_SLIP=0`

**Category:** DEVIATION / capture-completeness

**Description:**
The documented mechanism is the per-sub-step `[[v_n]]` OPENING leak; the
existing slip trace decomposes `sigma_n_trial = sn_vjump + sn_sterm` (opening
velocity jump vs bulk normal stress, `tpv205_substep_iterator.cpp:~427`).
`[DIAG-SIGN]` tracks only the aggregate `sigma_n_total` min — it localizes the
tensile DOF but cannot say whether the cause is `sn_vjump` (predictor/ghost
opening) or `sn_sterm` (bulk). The new nocap sbatch leaves that trace OFF
(`SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"`), so the run will not capture the
root-cause split at the located DOF.

**Trigger:** running the nocap sbatch as committed.

**Actual behavior:** localizes tensile DOFs; `sn_vjump`/`sn_sterm` never print.

**Expected behavior:** also emit the opening-vs-bulk decomposition at/near the
located speckle DOFs.

**Suggested fix:** enable the slip diagnostic with a threshold below the ~13 m/s
front peak so it fires on the speckle, not just the blow-up:
```diff
-export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-0}"
-export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-10.0}"
+export SEAS_DIAG_SLIP="${SEAS_DIAG_SLIP:-1}"
+export SEAS_DIAG_SLIP_VTHR="${SEAS_DIAG_SLIP_VTHR:-15.0}"
```
(Alternative: fold `sn_vjump`/`sn_sterm` into `[DIAG-SIGN-DOF]` so the split
prints at the located DOF without the full per-QP slip trace.)

**Test case:**
```python
def test_R002_causal_decomposition_emitted():
    assert any("sn_vjump" in l and "sn_sterm" in l for l in run_log)
```

---

### [R-003] [MODERATE] [POSSIBLE] spatial_dyn_driver.cpp:2136 — sub-step min reset every macro-step but read only at step%100 ‖ V>10

**Category:** BUG (capture-completeness)

**Description:**
`sigma_n_substep_min` is reset to `1e300` at the TOP of **every** macro-step
(:2136), but `[DIAG-SIGN]` reads it only when
`step % 100 == 0 || V_max_step > 10.0`. For ~99/100 macro-steps the sub-step
minimum is computed then overwritten unread. A tensile sub-step transient while
global `V_max < 10` (the speckle's FIRST onset) on a non-%100 step is invisible
— exactly the "did a sub-step go tensile before the output recovered?" regime.
POSSIBLE because once the front reaches V≈13 (~t=3 s) every step is sampled, so
the *sustained* window is covered; only the earliest pre-front onset is at risk.

**Trigger:** tensile sub-step at a step with `step%100!=0 && V_max<=10`.

**Actual behavior:** that step's min is reset away before any diag reads it.

**Expected behavior:** the tensile minimum should persist across the diag
interval (or be scanned every step).

**Suggested fix:** keep a driver-scope running interval-min, reset only after a
print (not every macro-step). Update it each macro-step from
`dof_data[i].sigma_n_substep_min`; report+reset in the `[DIAG-SIGN]` print.
```diff
   for (int i = 0; i < num_fault_total; ++i)
   {
      dof_data[i].slip_rate_substep_max = 0.0;
      dof_data[i].sigma_n_substep_min = 1.0e300;
   }
+  // (declared once, outside the time loop)
+  //   static real_t diag_iv_sn_min = 1.0e300; static int diag_iv_argmin = -1;
+  // updated every macro-step from dof_data[i].sigma_n_substep_min, reported
+  // and reset inside the rank-0 [DIAG-SIGN] print.
```
(If the maintainer accepts the V>10 gate covers the window of interest,
downgrade to LOW and document the pre-onset gap.)

**Test case:**
```python
def test_R003_transient_on_unsampled_step_captured():
    # Inject tensile sub-step at step 1050 (not %100, V<10), recover by 1051.
    assert diag_sign_at(1100).interval_sn_min <= injected_tensile_value
```

---

### [R-004] [LOW] multiple files — magic sentinel `1.0e300`/`1.0e299` instead of `std::numeric_limits`

**Category:** QUALITY / ASSUMPTION

**Description:**
`1.0e300` (DOFData default in `fault_face_flux.hpp`, driver reset `:2136`, loop
inits) and the guard `< 1.0e299` are hardcoded magic numbers — CLAUDE.md forbids
them; codebase idiom is `std::numeric_limits<real_t>::max()` (e.g. ComputeMaxDt).
Also fragile under `MFEM_USE_SINGLE` (`float` max ≈ 3.4e38): `1.0e300f` → `+inf`,
which works only by inf semantics.

**Suggested fix:**
```diff
-   real_t sigma_n_substep_min = 1.0e300;
+   real_t sigma_n_substep_min = std::numeric_limits<real_t>::max();
```
```diff
-         dof_data[i].sigma_n_substep_min = 1.0e300;  // [DIAG-SIGN] tensile tracker
+         dof_data[i].sigma_n_substep_min =
+            std::numeric_limits<real_t>::max();  // [DIAG-SIGN] tensile tracker
```
```diff
-            if (snss < 1.0e299)   // a sub-step value was recorded this macro step
+            if (snss < std::numeric_limits<real_t>::max())  // value recorded
```
(plus the matching driver `1.0e300` inits; ensure `<limits>` is included in
`fault_face_flux.hpp`).

**Test case:** N/A (covered by a single-precision build pass).

---

### [R-005] [LOW] [POSSIBLE] wave_operator.inl:4328 — shared-QP `sigma_n_total` is pre-average

**Category:** EDGE_CASE

**Description:**
On the shared-fault path `ComputeStageState` sets `sigma_n_total` from the
rank-LOCAL trial (`:182`), then the code averages `sigma_n_trial`
(`states[qq].sigma_n_trial = sn_avg`) WITHOUT recomputing `sigma_n_total`. So
`:4328` min-tracks the pre-average `sigma_n_total`, while end-of-step
`data.sigma_n_corr` for shared QPs uses the averaged path — slightly different
bases. Shared QPs are ~0.089 % of DOFs (162/182055), so negligible aggregate,
but a shared-QP speckle is reported with a marginally off σ_n.

**Suggested fix:**
```diff
+                  const real_t sn_tot_avg = fdata_qq.sigma_n0
+                                          + fdata_qq.sigma_n_nuc
+                                          + states[qq].sigma_n_trial;
                   fdata_qq.sigma_n_substep_min =
-                     std::min(fdata_qq.sigma_n_substep_min,
-                              states[qq].sigma_n_total);
+                     std::min(fdata_qq.sigma_n_substep_min, sn_tot_avg);
```
(Confirm `states[qq].sigma_n_trial` is the averaged value at this line.)

**Test case:**
```python
def test_R005_shared_qp_consistency():
    # On a >=2-rank run, quiescent shared-fault DOF:
    assert abs(shared_dof.sigma_n_substep_min - shared_dof.sigma_n_corr) < tol
```

---

### [R-006] [LOW] spatial_dyn_driver.cpp:[DIAG-SIGN] — every-step printing after V>10 bloats the log

**Category:** QUALITY

**Description:**
Once `V_max_step>10` (~t=3 s) the diag block fires every step through the event
(~30k steps), adding ≥2 `[DIAG-SIGN]*` lines/step (more with R-001's fix) atop
`[DIAG]`/`[DIAG-ONSET]`. The 7751832 log was already 8.9 MB with slip off; this
can multiply it and bury the signal.

**Suggested fix:** print the aggregate line on a coarser cadence or only when
tensile:
```diff
-         if (rank == 0)
+         if (rank == 0 && (step % 100 == 0 || n_ss_tensile_g > 0))
            { std::cout << "[DIAG-SIGN] step " << step << ... }
```

**Test case:** N/A (log-volume/quality).

---

## Summary
- Critical issues: 0
- Moderate issues: 3 (R-001 multi-spot; R-002 no opening/bulk decomposition +
  slip diag off; R-003 pre-onset sampling gap)
- Low issues: 3 (R-004 magic sentinel; R-005 shared-QP scale; R-006 log volume)
- Plan compliance: PARTIAL — the instrument correctly/consistently tracks
  per-sub-step tensile σ_n and floor engagement (scale + population checks
  pass) and WILL localize the worst tensile DOF once V>10. But "capture the bug
  location" is only partly met: one spot/step (R-001), no opening-vs-bulk cause
  (R-002), possible pre-front onset miss (R-003).
- Verdict: PASS WITH FIXES — usable to confirm tensile σ_n and locate the worst
  spot; apply R-001/R-002 (and decide R-003) before relying on it to fully
  localize/diagnose the speckle.

## Unreviewed Areas
- Live-run numerical behavior (Frontera unreachable; not executed). R-001/R-003
  capture claims should be re-checked against the first nocap+diag run.
- `tpv102`/`tpv104` byte-exact regression not re-run; the new field is
  transient/non-serialized and only min-updated, so byte-exactness is expected
  but unverified — confirm with `make test`.
