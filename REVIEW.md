# Code Review: Adaptive ParaView Output — Implementation Diff (2026-04-20)

## TL;DR

Implementation of `PLAN.md` (Adaptive ParaView Output for BP5 production). Fresh adversarial audit of the modified files. **All eight prior-round findings (R-P01..R-P08 in the plan review) are addressed in the code.** The implementer even improved on the plan's R-P07 fix by introducing a two-arg `CommitSchedule(time, V_max)` overload (plan specified only a one-arg shim update) — the two-arg form atomically advances `last_write_time_`, `current_regime_`, and `last_v_max_`, and the BP5 fault-only driver path uses it correctly.

One **MODERATE** finding: the TPV102 driver's `pv_no_domain` fault-only path still calls the single-arg `CommitSchedule(time)`. It is *safe today* (TPV102 always configures `output_every_n_steps` or `fixed_dt`, so the adaptive branch is never hit and `current_regime_` is irrelevant for scheduling), but the 1-arg shim's docstring explicitly warns that `PeekShouldWrite`-gated callers should use the 2-arg overload. The new docstring identifies TPV102 by name as a safe 1-arg caller — that claim is load-bearing on the "always step-based" property and would silently break if TPV102 ever enabled adaptive hysteresis.

Two **LOW** findings: a boundary-stickiness divergence at exact-threshold V values when hyst=1 (purely theoretical in floating point), and a missed opportunity to warn the user when hyst/v-threshold flags are silently ignored by step-based mode.

**Verdict: PASS WITH FIXES — R-I01 should be applied for forward-safety; R-I02 and R-I03 are LOW. The production Frontera submission is unblocked by this review.**

## Review Scope

- Plan: `/Users/chunhuizhao/projects/seas-mfem/PLAN.md` (475 lines)
- Files reviewed (working-tree diffs against HEAD):
  - `miniapps/seas/io/paraview_output.hpp` (+213 / -43)
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp` (+126 / -9)
  - `miniapps/seas/tests/unit/test_io.cpp` (+187 / -1)
  - `miniapps/seas/jobs/bp5/bp5_system_update_full_production.sbatch` (+18 / -1)
- Cross-referenced callers: `miniapps/seas/drivers/tpv102_driver.cpp` (unchanged), `miniapps/seas/tests/verification/bp2_serial_smoke.cpp` (unchanged, pre-existing API call pattern)
- Domain context: `CLAUDE.md` (repo root), `miniapps/seas/CLAUDE.md`, `miniapps/seas/ARCHITECTURE.md` (by reference).
- Memory constraints honored: no BP5 source rename, no local 1000m reproducer, Frontera approval still gated.

## Prior-Round Findings — Verification

| ID | Severity | Status | Evidence |
|----|----------|--------|----------|
| R-P01 | CRITICAL | **FIXED** | `bp5_verification_full.cpp:1654-1711` fault-only path populates all 5 `pv_local_*` vectors via 5 `ExpandOwnedToLocalFault` calls (plus `pv_local_normal_stress = 0.0` fallback), before `WriteFaultSurfaceVTU`. |
| R-P02 | MODERATE | **FIXED** | `bp5_verification_full.cpp:1565-1575` override block is placed inside the `if (use_paraview)` guard (:1510), between construction and the IC-write call at :1769. |
| R-P03 | MODERATE | **FIXED** | `paraview_output.hpp:134-151` case 2 drop semantics are documented; `test_io.cpp:541-573` `TestParaViewCoseismicDropSkipsNucleation` verifies the intentional skip of the [V_nu_exit, V_nu_enter] window from regime 2. |
| R-P04 | LOW | **FIXED** | `paraview_output.hpp:155-170` `Validate()` is a one-shot method; `bp5_verification_full.cpp:1573` calls it once after CLI overrides; `test_io.cpp:654-667` exercises it. |
| R-P05 | LOW | **FIXED** | `test_io.cpp:468` now uses `BP5Params::seconds_per_year` (matching the production constant in `paraview_output.hpp`). |
| R-P06 | LOW | **NOOP** | Flag-parsing style kept as plain `if`-chain per existing convention; not worth changing. No diff regression. |
| R-P07 | CRITICAL | **FIXED+IMPROVED** | `paraview_output.hpp:815-825` adds two-arg `CommitSchedule(time, V_max)` that advances regime; `bp5_verification_full.cpp:1698` uses the two-arg form; `test_io.cpp:594-634` `TestParaViewFaultOnlyAdvancesRegime` regression-tests the scenario. The fix goes beyond the plan's "CommitRegime" suggestion by atomizing the state update. |
| R-P08 | LOW | **FIXED** | `paraview_output.hpp:127,128,131` entry comparisons use strict `>` (not `>=`); `test_io.cpp:523-538` `TestParaViewLegacyBoundaryPreserved` locks in the V=1e-3 / V=1e-6 legacy semantics. |

## New Findings

---

### [R-I01] [MODERATE] [`drivers/tpv102_driver.cpp:728`] — TPV102 `pv_no_domain` path uses the deprecated 1-arg `CommitSchedule(time)` while the shim docstring warns PeekShouldWrite-gated callers to use the 2-arg form

**Category:** DEVIATION (API CONSISTENCY + LATENT TRAP)

**Description:**

The two-arg `CommitSchedule(time, V_max)` was introduced in this change to atomically advance `last_write_time_`, `current_regime_`, and `last_v_max_` after a `PeekShouldWrite`-gated fault-only write. The 1-arg shim was preserved at `paraview_output.hpp:839` with an explicit warning in its docstring (lines 830-838):

> Back-compat single-arg shim preserved for existing call sites (**TPV102 driver** and legacy BP5 paths). Delegates to the two-arg form using the last V_max recorded by Save/ShouldWrite; this is correct for callers that either (a) use the default hysteresis_factor = 1.0 (regime is stateless in V), or (b) call Save/ShouldWrite just before this. Callers that want hysteresis to advance through a PeekShouldWrite-gated path MUST call the two-arg overload instead.

TPV102 is identified by name as a safe 1-arg caller. The claim depends on two load-bearing properties:

1. **Property A**: TPV102 always sets `output_every_n_steps > 0` OR `fixed_dt > 0`, so the adaptive branch in `PeekShouldWrite` never runs. Verified: `drivers/tpv102_driver.cpp:634-644` unconditionally sets one of these two — the `else { pv_out->output_every_n_steps = output_interval; }` fallback at `:644` is NON-optional.
2. **Property B**: TPV102 never configures `hysteresis_factor > 1`. Verified: TPV102 has no CLI plumbing for the adaptive schedule.

Both properties hold today, so TPV102 is safe. BUT:

- A future TPV102 diagnostic that enables `sched.hysteresis_factor = 10` at runtime (e.g., to probe rupture ring-down) would silently defeat hysteresis — `last_v_max_` stays at 0.0 because `PeekShouldWrite` is const and never updates it, and 1-arg `CommitSchedule(time)` delegates to `CommitSchedule(time, last_v_max_=0.0)` → `NextRegime(0.0, 0) = 0` → `current_regime_` pinned at 0 forever. The bug would appear out of nowhere: no compile error, no abort, no warning.
- The inconsistency also weakens the new API's self-documenting property: "all PeekShouldWrite-gated writes should use 2-arg CommitSchedule" is the design invariant, and leaving TPV102 on the 1-arg path is exactly the kind of exception that accumulates into confusing legacy-vs-modern code.

**Trigger:** Not today. Future regression when anyone adds adaptive scheduling to TPV102.

**Actual behavior:** Safe by accident (property A above).

**Expected behavior:** Mirror the BP5 fault-only pattern. The two-arg `CommitSchedule(time, V_max)` is already in the header; switching the TPV102 call site is a one-line change that eliminates the latent trap.

**Suggested fix:**

```diff
 drivers/tpv102_driver.cpp:724-729
       if (pv_no_domain)
       {
          // Fault-surface PVD only; advance the schedule ourselves.
-         pv_out->CommitSchedule(time);
+         // Use the 2-arg overload so hysteresis (if ever enabled) advances
+         // correctly in this PeekShouldWrite-gated path.  Matches the BP5
+         // fault-only pattern at bp5_verification_full.cpp:1698.
+         pv_out->CommitSchedule(time, V_max);
       }
```

While at it, also update the docstring on the 1-arg shim to remove TPV102 from the "safe caller" list (since TPV102 will no longer use it):

```diff
 io/paraview_output.hpp:830-838
    /// Back-compat single-arg shim preserved for existing call sites
-   /// (TPV102 driver and legacy BP5 paths).  Delegates to the two-arg
+   /// (legacy BP2 smoke + bp2_serial_smoke test path).  Delegates to the
+   /// two-arg form using the last V_max recorded by Save/ShouldWrite;
    /// ...
```

**Test case:**

```cpp
// tests/unit/test_io.cpp — add regression test
void TestParaViewTPV102StyleFaultOnlyAdvancesRegime()
{
   std::cout << "\n=== Test: TPV102-Style Fault-Only Advances Regime ===\n";
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 1, Element::HEXAHEDRON,
                                     2.0, 2.0, 1.0);
   ParaViewOutput<Mesh> pv("test_pv_tpv102_style", mesh, 1);
   auto &s = pv.GetSchedule();
   s.v_coseismic       = 1e-3;
   s.v_nucleation      = 1e-6;
   s.hysteresis_factor = 10.0;        // NOT default — simulate hypothetical TPV102 adaptive use
   s.dt_coseismic      = 0.5;
   s.dt_nucleation     = 10.0;
   s.dt_interseismic   = 1.0 * BP5Params::seconds_per_year;
   s.Validate();

   // Drive the TPV102 fault-only pattern: PeekShouldWrite then CommitSchedule.
   // With the current 1-arg call this fails (last_v_max_ = 0 → regime stuck at 0).
   // With the 2-arg call (R-I01 fix) this passes.
   TEST_ASSERT(pv.PeekShouldWrite(0, 0.0, 2e-3), "step 0: IC write fires");
   pv.CommitSchedule(0.0, 2e-3);    // R-I01 fix: use 2-arg

   TEST_ASSERT(pv.PeekShouldWrite(1, 0.5, 2e-3), "step 1: dt_co elapsed");
   pv.CommitSchedule(0.5, 2e-3);

   TEST_ASSERT(pv.PeekShouldWrite(2, 1.0, 2e-4),
               "step 2: hysteresis keeps dt_coseismic even as V dropped to 2e-4");
}
```

This test, run against the current TPV102 driver's 1-arg path (substitute `pv.CommitSchedule(time)` — no 2nd arg), would fail. Run against the R-I01 fix, it passes.

---

### [R-I02] [LOW] [`paraview_output.hpp:NextRegime case 2`] — Boundary stickiness at exact-V-threshold values with hysteresis_factor=1 diverges from legacy `OutputInterval` semantics for stateful Save() flows

**Category:** ASSUMPTION (DOCUMENTATION CLAIM OVERSTATED)

**Description:**

The header comment at `paraview_output.hpp:109-111` claims:

> Exit comparisons are strict <, so an exact-threshold V stays in its current regime.

This implementation choice, combined with strict `>` for entry (good, per R-P08), creates a *stickiness* at V exactly equal to `v_coseismic` (default 1e-3) from prev=2 that the legacy `OutputInterval` never had:

- Legacy `OutputInterval(V=1e-3)`: `V_max > 1e-3` is false → returns nucleation (`1.0s`).
- New `NextRegime(V=1e-3, prev=2)`: case 2, `V < V_co_exit=1e-3` is false (strict <) → stays at 2 (coseismic, dt=0.01s).
- New `NextRegime(V=1e-3, prev=0 or 1)`: returns 1 (nucleation) via strict `>` check — matches legacy.

So the static helper (`OutputInterval`, which internally starts from `prev=0`) matches legacy, but the stateful `Save()` flow with `prev=2` DIVERGES at V=1e-3 exactly. The `TestParaViewLegacyBoundaryPreserved` test only covers the static-helper path, not the stateful path.

Plan constraint line 33 claims byte-compatibility for BP5/BP2 smoke tests. In floating-point arithmetic, V_max never hits exactly 1e-3, so the divergence is purely theoretical — tests pass and production runs are unaffected.

**Trigger:** A test / production run where V_max arrives at exactly `v_coseismic` from a prior-step regime 2. Requires exact floating-point equality, which does not occur in practice.

**Actual behavior:** New stateful path gets one extra write at dt_coseismic cadence at the boundary.

**Expected behavior:** Either (a) strictly match legacy by using `<=` on the exit check when prev=2 (symmetric with legacy's `V_max > 1e-3` — which excludes V=1e-3 from coseismic), or (b) keep the current behavior and clarify the docstring that byte-compat applies to the static helper, not the stateful Save flow.

**Suggested fix (option b, minimal):** Update the header comment at `paraview_output.hpp:107-111` to be precise:

```diff
       /// Entry comparisons are STRICT > (not >=) to preserve byte
       /// compatibility with the legacy OutputInterval which used
-      /// V_max > 1e-3 / V_max > 1e-6.  Exit comparisons are strict <,
-      /// so an exact-threshold V stays in its current regime.
+      /// V_max > 1e-3 / V_max > 1e-6.  Exit comparisons are strict <,
+      /// so an exact-threshold V stays in its current regime — this
+      /// introduces stickiness in the STATEFUL Save() path at V exactly
+      /// equal to a threshold from prev=2 that the (stateless) legacy
+      /// OutputInterval did not have.  Byte-compat is preserved for the
+      /// STATIC helper OutputInterval(V) (which always starts prev=0),
+      /// not for the stateful Save() flow.  In practice V never hits
+      /// exactly 1e-3 in floating-point, so the divergence is theoretical.
```

**Test case:** N/A — divergence only manifests at bit-exact threshold values which do not arise in any real run.

---

### [R-I03] [LOW] [`bp5_verification_full.cpp:1576-1608`] — Log line silently omits the active hysteresis / V-threshold settings when step-based or fixed-dt mode short-circuits the adaptive schedule

**Category:** QUALITY (OBSERVABILITY)

**Description:**

If the user passes `--paraview-every 100 --paraview-hyst 10 --paraview-v-co 1e-3`, the adaptive schedule is silently bypassed (step-based takes precedence per `Save()`'s `output_every_n_steps` short-circuit at `paraview_output.hpp:742-750`), but the CLI override block at `bp5_verification_full.cpp:1565-1574` still applies the hyst / V overrides. The log line at `:1580-1607` prints:

```
  ParaView output: ON (every 100 steps)
```

— WITHOUT echoing the hyst=10 or v_co=1e-3 values. The user's mental model ("I enabled hysteresis") does not match the runtime behavior ("hysteresis is a no-op, every 100 steps writes unconditionally"). No warning, no diagnostic.

**Trigger:** User combines step-based or fixed-dt flags with adaptive-schedule overrides.

**Actual behavior:** Overrides silently applied to `sched.*` fields but never used.

**Expected behavior:** Either (a) warn when the mixed combination is detected, or (b) unconditionally echo the adaptive-schedule values in the log so the user sees what's ignored.

**Suggested fix:**

```diff
 bp5_verification_full.cpp:1580-1607
       if (mpi.IsRoot())
       {
+         const auto &s = pv_out->GetSchedule();
          if (paraview_step_interval > 0)
          {
             std::cout << "  ParaView output: ON (every "
                       << paraview_step_interval << " steps)\n";
+            if (pv_hyst > 0 || pv_v_co > 0 || pv_v_nu > 0 ||
+                pv_dt_co > 0 || pv_dt_nu > 0 || pv_dt_inter > 0)
+            {
+               std::cout << "    WARN: adaptive flags (hyst/v/dt) are "
+                         "IGNORED in step-based mode.\n";
+            }
          }
          else if (paraview_dt > 0.0)
          {
             std::cout << "  ParaView output: ON (every "
                       << paraview_dt / BP5Params::seconds_per_year
                       << " yr)\n";
+            if (pv_hyst > 0 || pv_v_co > 0 || pv_v_nu > 0 ||
+                pv_dt_co > 0 || pv_dt_nu > 0 || pv_dt_inter > 0)
+            {
+               std::cout << "    WARN: adaptive flags (hyst/v/dt) are "
+                         "IGNORED in fixed-dt mode.\n";
+            }
          }
          else
          {
-            const auto &s = pv_out->GetSchedule();
             std::cout << "  ParaView output: ON (adaptive schedule"
```

**Test case:** Manual — run `seas_bp5_full --paraview-every 100 --paraview-hyst 10` and grep the stdout for "WARN". Not worth a unit test.

---

## Summary

- **Critical issues: 0** (all prior CRITICAL findings R-P01, R-P07 addressed with evidence).
- **Moderate issues: 1** (R-I01 — TPV102 1-arg/2-arg inconsistency, safe today but latent).
- **Low issues: 2** (R-I02 boundary stickiness doc-claim; R-I03 silent adaptive-flag override in step/dt modes).
- **Plan compliance: FULL** — all 6 phases implemented; all 8 prior-round findings (R-P01..R-P08) addressed; one positive deviation (R-P07 fix via atomic 2-arg `CommitSchedule(time, V_max)` is stronger than the plan's `CommitRegime` suggestion).
- **Memory-constraint compliance: FULL** — no BP5 source rename (additive changes only to shared `paraview_output.hpp`), no local 1000m runs in tests, Frontera submission still requires approval, `--paraview-fault-only` disk-footprint rationale documented inline in the production sbatch.
- **Verdict: PASS WITH FIXES — R-I01 should be applied for forward-safety (symmetry + remove latent trap); R-I02 and R-I03 are LOW and can defer.** The production Frontera submission is unblocked by this review.

## Unreviewed Areas

- **`WriteFaultSurfaceVTU` with empty local fault set** (ranks with zero fault faces in 400-rank partitioning). Existing code at `paraview_output.hpp:410-620` is presumed to handle via the `has_fault_output_` guard, not re-audited here.
- **`ForceSave(cycle, time)` in adaptive mode**: does not update `current_regime_` or `last_v_max_` (`paraview_output.hpp:912-917` only advances `last_write_time_`). If someone uses PeekShouldWrite + ForceSave in adaptive mode, hysteresis state is stale. Not triggered by this PR; noted for future.
- **`bp2_serial_smoke.cpp:416-417`** uses `UpdateFaultFields(...)` (no BP5 suffix) and `Save(step, t)` (2-arg), neither of which matches the current header API. This is pre-existing breakage (not introduced by this PR, not modified). Out of scope for this review; worth a separate cleanup ticket.
- **Checkpoint-restart interaction with `current_regime_`**: on restart, regime resets to 0 (constructor default) — the post-restart first few writes could miscategorize if the restart happens mid-rupture. Same-severity note as in the prior plan review; not in this PR's scope.
- **`TestParaViewOutputInterval` + `TestParaViewLegacyBoundaryPreserved` combined coverage**: tests the static helper only. The stateful Save flow's regime-stickiness at exact V thresholds (R-I02) is not directly covered. Theoretical risk only.
