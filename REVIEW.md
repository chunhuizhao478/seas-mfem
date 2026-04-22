# Code Review: Claim verification — "fault observables written as RK4 stage averages"

**Date:** 2026-04-22 (round-9; targeted claim check)
**Reviewer:** code-review agent
**User ask:** *"Can you confirm this is still true? — High: the driver is writing fault observables as RK4 stage averages, not as values re-evaluated from the final state Q(t+dt). In `.../tpv102_driver.cpp:944`, `tau1_corr`, `tau2_corr`, and `sigma_n_corr` are overwritten with `(k1 + 2k2 + 2k3 + k4)/6`. Then `.../tpv102_setup.hpp:289` writes those fields directly to station output, and `.../tpv102_driver.cpp:1051` does that immediately after the RK4 step. That is a semantic bug for instantaneous outputs: stage-averaged traction/stress is not the same thing as traction/stress at Q^{n+1}. This is the strongest code-level explanation I see for the inflated sigma_n signal."*

## Review Scope

- `miniapps/seas/drivers/tpv102_driver.cpp` — RK4 loop line 944 (for-loop header), the averaging assignments at lines 971-973, and the station write at line 1051-1053.
- `miniapps/seas/dynamic/tpv102_setup.hpp` — `TPV102StationWriter::WriteStep` at line 289.
- Prior review rounds — round-3 R-V92-E02 (H-V92-K hypothesis), round-4 F01/F02/F03 (magnitude dimensional analysis), round-7 R-V92-I01 (_k4 discriminator broken).

## Claim verification

### Part 1 — "tpv102_driver.cpp:944 overwrites `tau*_corr` and `sigma_n_corr` with (k1 + 2k2 + 2k3 + k4)/6"

**VERIFIED WITH ONE CAVEAT.** Line 944 is the FOR-LOOP HEADER:
```cpp
for (int i = 0; i < num_fault_total; i++)
```

The actual assignments are at lines 971-973 (inside that loop):
```cpp
dof_data[i].tau1_corr    = (t1c_k1[i] + 2*t1c_k2[i] + 2*t1c_k3[i] + t1c_k4[i]) / 6.0;
dof_data[i].tau2_corr    = (t2c_k1[i] + 2*t2c_k2[i] + 2*t2c_k3[i] + t2c_k4[i]) / 6.0;
dof_data[i].sigma_n_corr = (snc_k1[i] + 2*snc_k2[i] + 2*snc_k3[i] + snc_k4[i]) / 6.0;
```
Also overwrites `V1`, `V2`, `slip_rate`, plus `slip1 += V1·dt`, `slip2 += V2·dt`.  The `t1c_k_i`, `t2c_k_i`, `snc_k_i` are stage-i SNAPSHOTS of `dof_data[i].tau1_corr/tau2_corr/sigma_n_corr` captured immediately after the corresponding `wave.Mult` call at each stage (lines 825-826, 850-851, 866-867, 898-899 in the current source).  So yes — the final `tau*_corr` and `sigma_n_corr` values written to `dof_data` at step end are the Butcher-weighted Simpson-average of stage-wise post-friction outputs.

Inline comment at line 970 confirms the intent:
```cpp
// R-001 fix: RK4-weighted corrected tractions for consistent station output
```

### Part 2 — "tpv102_setup.hpp:289 writes those fields directly to station output"

**VERIFIED.** `WriteStep` at line 289 reads from `dof_data[idx]` (the station's owner DOF) and writes to the `.dat` file at lines 305-314:
```cpp
files_[s] << std::scientific << std::setprecision(10)
          << t << " "
          << d.slip1 << " " << d.slip2 << " "
          << d.V1 << " "    << d.V2    << " "
          << d.tau1_corr << " " << d.tau2_corr << " "
          << d.sigma_n_corr << " " << log10_theta << "\n";
```
`d = dof_data[idx]` is assigned at line 298.

So `tau1_corr`, `tau2_corr`, `sigma_n_corr` in the `.dat` files ARE the step-end Simpson-averaged values — NOT re-evaluated at `Q(t+dt)`.

### Part 3 — "tpv102_driver.cpp:1051 writes immediately after the RK4 step"

**VERIFIED.** Line 1051-1053:
```cpp
if (step % output_interval == 0 || step == nsteps - 1)
{
   station_writer.WriteStep(t, dof_data);
   surface_writer.WriteStep(t, Q);
   ...
}
```
Timeline within one loop iteration:
1. Stages 1-4 Mult calls + stage-i captures into `t1c_k_i`, `t2c_k_i`, `snc_k_i` buffers.
2. Q update with `(k1+2k2+2k3+k4)/6` (line 906-909).
3. `t += dt_step` (line 911).
4. _k4 stage-4 snapshot to VTU buffers (line 919-929; current work, see round-7 R-V92-I01 for the overwrite bug).
5. RK4 averaging of `dof_data[i].{psi, V1, V2, tau1_corr, tau2_corr, sigma_n_corr, slip1, slip2}` (lines 944-974, including the claim's lines 971-973).
6. V_max tracking, diagnostic blocks.
7. Output write at line 1051 via `station_writer.WriteStep(t, dof_data)`.

So the station write uses `dof_data` whose tau/σ_n fields were just Simpson-averaged.  The `t` label written to the `.dat` file is `t_{n+1}` (post-increment at step 3), but the values are evaluated at `t_n + dt/2` (Simpson mean of stage-i samples at stage times `t_n`, `t_n+dt/2`, `t_n+dt/2`, `t_n+dt`).

## Claim status: confirmed as a real bug, severity downgraded

### The factual observation is CORRECT

All three line references are accurate (modulo line 944 being the for-loop header, not the assignment itself).  The station `.dat` output at time `t` is the Simpson-mean of the stage-wise post-friction values over `[t_n, t_n+dt]`, i.e. the midpoint-like approximation of `tau/σ_n` at `t_n + dt/2`.  This IS a semantic mismatch between the time label and the value.

### The severity label "HIGH — strongest explanation for the inflated σ_n" is OVERSTATED

Quantitative argument (mirrors round-4 REVIEW R-V92-F03 dimensional analysis):

- The Simpson-mean differs from the endpoint value by approximately:
  `|tau(t_n+dt/2) − tau(t_n+dt)| ≈ (dt/2) · |∂tau/∂t|`
- At peak rupture, `|∂σ_n/∂t|` ~ O(10 GPa/s) locally near the rupture tip.  At `dt ~ 10⁻⁴ s`, the per-step phase-lag offset is `5·10⁻⁵ · 10¹⁰ = 5·10⁵ Pa = 0.5 MPa`.
- **This is a per-step SNAPSHOT offset, not a cumulative drift.**  Each step's output is an independent sample with an independent midpoint bias; the biases do NOT add.
- Observed σ_n deviation at `flt_0_7.5` is `~98 MPa average (late)`, `218 MPa peak`.
- A ~0.5 MPa phase-lag cannot explain a 98 MPa sustained drift — off by 200×.

### Additional factor the original claim missed: the averaged values do NOT feed back into the simulation

The claim implicitly worries that stage-averaged `tau1_corr` etc. contaminates the next step.  Inspection of the code shows they DO NOT:

- At the start of the next RK4 step, stage-1 `wave.Mult(Q, k1)` internally calls `fault_flux_->Evaluate(fdata, Q_plus, Q_minus, ...)`.  `Evaluate` computes `sigma_n_trial`, `tau1_trial`, `tau2_trial` from the CURRENT Q's velocity and stress jumps (Pelties eq. 7), then runs the friction Brent solve from the current `fdata.psi` and `fdata.sigma_n0/tau1_0/tau2_0`.  It OVERWRITES `fdata.sigma_n_corr`, `tau1_corr`, `tau2_corr`, `V1`, `V2` with the stage-1 friction-corrected values.
- The previous step's Simpson-averaged values are IMMEDIATELY replaced before anything reads them.
- Only `fdata.psi` (integrated per F01+F02 via coupled RK4 on psi — round-4 fix) and `fdata.slip1/slip2` (accumulated correctly via `slip += V_avg·dt` — Simpson-exact) survive step-to-step.

**Therefore:**
- The averaging is a **pure output transformation** — simulation dynamics are unaffected.
- The phase-lag is a cosmetic / visualization bug, not a driver of the pathology.
- H-V92-K (round-3 hypothesis, round-4 F03 bounded at ~1% of observed pathology) is the correct classification; this claim's "HIGH / strongest explanation" severity is not supported by dimensional analysis.

## Findings

### [R-V92-K01] [MODERATE] [tpv102_driver.cpp:971-973 + tpv102_setup.hpp:305-314] — Station output shows Simpson-mean tau/σ_n at time `t_n + dt/2` but labels the sample time as `t = t_{n+1}`; half-step phase lag in observables

**Category:** BUG (semantic — time label ≠ value time)

**Description:**
Lines 971-973 overwrite `dof_data[i].tau1_corr/tau2_corr/sigma_n_corr` with the RK4-Butcher-weighted Simpson mean of stage-i post-friction snapshots.  The station writer at `tpv102_setup.hpp:289-314` writes these values alongside the current time `t` (post-increment).  The Simpson mean approximates the time-averaged value over `[t_n, t_{n+1}]`, which equals the midpoint value `tau(t_n + dt/2)` to O(dt⁴) for smooth tau(τ).  But the output file labels it as `t = t_{n+1}`.  Half-step phase lag.

**Trigger:**
Any station .dat output during a rupture where `|∂tau/∂t|` or `|∂σ_n/∂t|` are non-trivial.

**Magnitude (dimensional bound):**
Per-sample bias ≈ `dt/2 · |∂tau/∂t|`.  For `dt = 10⁻⁴` s and peak rupture-tip `|∂σ_n/∂t| ~ 10 GPa/s`: bias ≈ 0.5 MPa per sample.  NOT cumulative.

**Impact on observed 218 MPa σ_n peak:** bias is 200× smaller than the observed deviation.  This bug alone cannot explain the primary pathology.

**Suggested fix:**
Two options, pick one:

**(A) Label the time correctly.**  If keeping the Simpson-mean semantics (useful for time-averaged validation), write `t - dt/2` in the `.dat` file:

```diff
--- a/miniapps/seas/dynamic/tpv102_setup.hpp
+++ b/miniapps/seas/dynamic/tpv102_setup.hpp
@@ -305,1 +305,5 @@
-         files_[s] << std::scientific << std::setprecision(10)
+         // R-V92-K01 (round-9 REVIEW): tau*_corr and sigma_n_corr
+         // stored in dof_data are the RK4-Butcher-weighted Simpson
+         // mean over [t-dt, t], which approximates the midpoint value
+         // tau(t - dt/2).  Write the midpoint time to match.
+         files_[s] << std::scientific << std::setprecision(10)
-                   << t << " "
+                   << (t - 0.5 * last_dt_) << " "    // Simpson midpoint
                    << d.slip1 << " "
```

(requires plumbing `last_dt_` through the writer — minor API change.)

**(B) Re-evaluate tau/σ_n at Q(t+dt) after the RK4 step.**  Add a final post-averaging `fault_flux_->Evaluate` call on the updated Q (not on stage-i buffers) to produce endpoint tractions:

```diff
--- a/miniapps/seas/drivers/tpv102_driver.cpp
+++ b/miniapps/seas/drivers/tpv102_driver.cpp
@@ -974,0 +974,20 @@
+#ifdef SEAS_OUTPUT_ENDPOINT_TRACTIONS
+      // R-V92-K01 (round-9 REVIEW): replace the Simpson-mean tau*_corr
+      // and sigma_n_corr with values re-evaluated at Q(t+dt).  This
+      // removes the half-step phase lag between the time label and the
+      // output value.  Gated on a compile flag so the default behaviour
+      // (Simpson mean) is unchanged pending broader review.
+      {
+         Vector k_unused(Q.Size());
+         wave.Mult(Q, k_unused);                  // triggers Evaluate at Q(t+dt)
+         // dof_data[i].{tau1_corr, tau2_corr, sigma_n_corr, V1, V2,
+         //              slip_rate} are now endpoint values.
+         // Slip accumulators have already been updated with the
+         // Simpson-mean V_avg on line 968-969 (that IS the correct
+         // RK4 integral of V over the step — do NOT re-accumulate).
+      }
+#endif
```

(Adds one extra `wave.Mult` per output-interval step — modest cost.)

**Recommended:** Option (A) — because endpoint re-evaluation via option (B) has its own issues (Mult has side effects beyond the friction state, and the extra friction Brent call adds iteration count).  Option (A) is a single-line fix that makes the claim false; the station output then correctly labels `t = t_n + dt/2`, aligning value and label.

**Test case:**
```python
def test_R_V92_K01_station_time_label_matches_value_time():
    # Synthetic TPV102 run where tau evolves linearly: tau(τ) = a + b·τ.
    # After one RK4 step over [t_n, t_n+dt], dof_data.tau*_corr is set to
    # the Simpson mean tau_avg = a + b·(t_n + dt/2).
    #
    # Before fix: station .dat row reads "t_{n+1}  tau_avg = a + b·(t_n+dt/2)"
    #             — the tuple (t_label, tau_value) misrepresents tau(t_label)
    #             = a + b·t_{n+1} by 0.5·b·dt.
    #
    # After fix (A): station .dat row reads "t_n+dt/2  tau_avg" — self-
    #                consistent.
    run_one_step_linear_tau_fixture(a=100e6, b=1e9, dt=1e-4)
    for row in load_station_dat("flt_0_7.5.dat"):
        expected_tau_from_label = a + b * row.t
        # Within one step: should match to O(dt^4) given linear input.
        assert abs(row.tau1 - expected_tau_from_label) < 1e3, \
            f"phase lag: tau({row.t}) = {row.tau1}, expected {expected_tau_from_label}"
```

---

### [R-V92-K02] [LOW] [tpv102_driver.cpp:970 comment] — In-source comment "RK4-weighted corrected tractions for consistent station output" understates the semantic shift

**Category:** QUALITY (readability / misleading comment)

**Description:**
Line 970:
```cpp
// R-001 fix: RK4-weighted corrected tractions for consistent station output
```

"Consistent station output" suggests "the output is in sync with the simulation state".  In fact the output is in sync with the MIDPOINT of the step, while the time label is in sync with the ENDPOINT.  A clearer comment would surface the tradeoff:

```diff
-         // R-001 fix: RK4-weighted corrected tractions for consistent station output
+         // R-001 fix: RK4-weighted corrected tractions for consistent station output.
+         // NOTE (R-V92-K01): the Butcher (1,2,2,1)/6 weights produce the
+         // SIMPSON TIME-AVERAGE of tau*_corr over [t_n, t_{n+1}], i.e.
+         // tau*_corr(t_n + dt/2) to O(dt^4), NOT the endpoint value
+         // tau*_corr(t_{n+1}).  The half-step phase lag between the
+         // output value and the output time label is an open R-V92-K01
+         // bug — see REVIEW round-9 for the fix options.
         dof_data[i].tau1_corr = (t1c_k1[i] + 2*t1c_k2[i] + 2*t1c_k3[i] + t1c_k4[i]) / 6.0;
```

---

## Verdict on the user's claim

- **Factual part ("lines 944, 289, 1051"):** **CONFIRMED** with the minor caveat that 944 is the for-loop header and the actual assignments are 971-973.
- **Severity label "HIGH":** **DOWNGRADED to MODERATE.**  The bug is real (half-step phase lag between output value and output time label) but cannot account for the 218 MPa σ_n peak:
  - Per-step bias bound: ~0.5 MPa (at dt=10⁻⁴ s, |∂σ_n/∂t|~10 GPa/s).  200× smaller than the observed signal.
  - Non-cumulative: each sample has an independent midpoint bias.
  - No simulation feedback: `Evaluate` at stage 1 of the NEXT step immediately overwrites the averaged values, so dynamics are unaffected.
- **"Strongest code-level explanation for the inflated σ_n signal":** **NO.**  This is H-V92-K, which round-4 R-V92-F03 dimensional analysis bounded at ~1% of the observed pathology.  The primary driver remains in the flux-path (H-V92-G / H-V92-U) — see round-5 R-V92-G02 next-step tree.

## Recommended next step

1. Apply R-V92-K01 fix (option A — single-line change in `WriteStep`).  Closes the claim as a legitimate bug while correctly scoping its magnitude.
2. Continue with the round-7 R-V92-I01 `_k4` VTU fix (still blocking the H-V92-K empirical check).
3. Continue with the round-5 R-V92-G01 dt-halving experiment, which is the decisive discriminator for time-integration vs flux-path primary.

The user's claim identifies a bug that should be fixed.  It should NOT be used as a reason to skip R-V92-I01 or the dt-halving test; those address the primary pathology, which this bug does not explain.

## Summary

- **Critical:** 0
- **Moderate:** 1 (R-V92-K01)
- **Low:** 1 (R-V92-K02)
- **Plan compliance:** N/A (claim verification, not plan audit).
- **Verdict:** **PASS WITH FIXES.**  The claim identifies a real but narrower bug than characterized; severity label "HIGH / primary cause" is not supported by dimensional analysis.  Fix R-V92-K01, continue with the primary-pathology diagnostics.

## Unreviewed areas

- **Alternative fix path:** moving the `fault_flux_->Evaluate` call that produces `tau*_corr` out of `wave.Mult` and into a dedicated driver-level "endpoint friction" call after the RK4 Q update.  Would eliminate R-V92-K01 at the source but is a larger refactor.  Not audited.
- **Whether the `surface_writer` at line 1054 has the same issue** — `surface_writer.WriteStep(t, Q)` takes Q (the endpoint state), not `dof_data`, so it is at t_{n+1} correctly.  No bug there.
- **Checkpoint/restart consistency** — if the simulation checkpoints `dof_data` at `t` after the averaging, restart reads back the Simpson-mean values and uses them to initialize the next step.  Round-1 of the next segment's stage-1 Mult overwrites them immediately, so no data corruption.  Not a fresh bug.
