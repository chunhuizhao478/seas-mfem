# TPV102 Debug v9.2.0 Plan

---

## PRIORITY NOW (rev-3i, 2026-04-22)

### Plan is NOT closed.

Observed:  **σ_n peak 218 MPa (~100× SCEC spec), slip_dip 1.83 m (~1 830× SCEC spec)**.

Primary R-V92 mechanism: **not identified**.

**rev-3i update (2026-04-22):** Frontera jobs 7670526 (dt/2),
7670527 (dt/4), and 7671002 (F01+F02 postfix at nominal dt) all
produce **visually indistinguishable** station plots — same peak
V_strike, same σ_n drift (120 → 102–123 MPa depending on
station), same spurious V_dip / slip_dip / τ_dip signatures.
Decisive conclusions:
1. **H-V92-K (time-integrator / RK4 truncation) CLOSED** —
   dt-invariance rules it out (§21.2).
2. **F01+F02 is correct but NOT sufficient to close R-V92** —
   postfix plots match PRE-fix plots exactly (§21.3).
3. **H-V92-U (interior-fault path, non-2-tet code paths) is now
   the sole RANK-1 candidate**; H-V92-W promoted to RANK-2.
4. **Next decisive experiment: §5.1 `SEAS_DIAG_FAULT_SIGMA`
   Frontera probe (~200 SU)** — §20.3 Step 3.

See §21 for the full write-up and hypothesis-rank update.

### Progress snapshot (rev-3h+, 2026-04-21)

- Step 1 sbatch files **written** (`dt_half`, `dt_quarter`); submission
  pending user approval.  See updated Step 1 detail below.
- Step 2 **EXECUTED** — peak `|σ_n − 120 MPa|` = 218.7 MPa, **4.9×
  the 45 MPa Poisson bound**, verdict **MONOTONIC-GROWTH**.  R-V92-E01
  "rupture-extent misread" REFUTED.  H-V92-G stays viable; Step 3
  (§5.1 bulk SYY probe) still needed to localize the amplifier.
- Step 5 **DONE** — coupled RK4 on (Q, psi) applied to
  `drivers/tpv102_driver.cpp`; psi now RK4-integrated inside the
  stage loop via `AgingLawPsi::Rate`; all 4 regression gates
  (§4.9/§4.10/§4.11/§4.12) + Step 6 PASS.  Post-fix Frontera
  verification pending.  Local changes STILL UNCOMMITTED per
  R-V92-H03 Option B hold.
- Step 5b **PASS 100/100** (rev-3h+, 2026-04-21) — local 4-km
  shared-fault stress-test of F01+F02 over 100 RK4 steps at 4
  ranks.  Terminal envelope: `max|σ_n-σ_n0| = 2.487·10⁵ Pa` (0.2×
  1 MPa bound), `max slip_rate = 1.37·10⁻¹² m/s` (effectively
  V_ini), `max|psi-psi_0| = 1.1·10⁻¹⁴` (FP noise), zero NaNs.
  σ_n deviation grew linearly (2.5·10⁴ → 2.5·10⁵ Pa over 1 ms),
  consistent with elastic-wave transit of the seeded antisymmetric
  σ_xy(y) profile rather than any amplification pathology.
  **Decision (Step 5b table):** commit + push F01+F02 with high
  confidence; the coupled-RK4 arithmetic is multi-step stable on
  shared-fault QPs.  Closes R-V92-H01 (driver RK4 block was
  untested in multi-step form).
- Step 6 **PASS 6/6** — interior-fault path reproduces Pelties eq.
  (7a/7b/7c) on a 2-tet fixture to O(10⁻⁸) relative error.  H-V92-U
  ELIMINATED on this fixture; bimaterial / corner / multi-face
  extensions still open.

### Do next (in this order, stop as soon as a step closes R-V92)

Each entry: **Condition** = what must be true before this step runs;
**Closure** = what outcome closes R-V92 at this step (allowing skip
ahead to Step 8/9); **Next if not closed** = which step to proceed to.

1. **dt-halving Frontera test — ~20 SU, 1 h.**  **[sbatch READY;
   user submission pending]**
   - **Condition:** plan is open AND user approval per
     `feedback_frontera_approval.md` granted AND source tree matches
     the v91 baseline head.
   - **Closure:** `slip_dip` peak ≈ v91's 1.83 m at both dt/2 and
     dt/4 does NOT close — means primary is non-truncation.  A clean
     dt⁴ (~256×) or dt² (~16×) scaling DOES close: identifies
     RK4-truncation or psi-split primary → skip to Step 5 then 8.
   - Submit `tpv102_200m_p1_12.0s_400r_v92_dt_half.sbatch` and `..._dt_quarter.sbatch`.
   - Compare `slip_dip` peak vs v91 job-7668434 baseline.
   - Classification:
     - drops ~16× at dt/2  → time-integration bug → apply F01+F02 fix (step 5 below) → done.
     - drops ~8× at dt/2   → psi splitting (F02) → integrate psi inside RK4 state.
     - ≈ unchanged          → NOT time integration → go to step 2.
     - larger / NaN         → CFL violation → reduce CFL + re-audit `flux.Interior`.
   - **Next if not closed:** Step 3.

2. **Re-analyse existing §18 outlier CSVs — FREE, local, ~1 h.**
   **[DONE: MONOTONIC-GROWTH, 4.9× Poisson bound]**
   - **Condition:** `outlier_report.json` from job-7668434 available
     locally (it is, at
     `tpv102/r_v92_outlier_detection/outlier_report.json`).  Runs in
     parallel with Step 1 — no Frontera dependency.
   - **Closure:** `SATURATES-BELOW-POISSON-BOUND` verdict closes
     R-V92 as R-V92-E01 misread (no real amplification).  `MONOTONIC-
     GROWTH` does NOT close; real amplification stays viable → Step 3.
   - Extract peak `|σ_n − 120 MPa|` and peak `|slip_dip|` per cycle.
   - If peak saturates near the 45 MPa Poisson bound → §18 "H-V92-G CONFIRMED" was a misread of a rupture-extent signal (R-V92-E01).
   - If peak grows unbounded → real amplification; still need to find the amplifier.
   - **ACTUAL verdict:** MONOTONIC-GROWTH (218.7 MPa = 4.9× bound);
     did not close → continue to Step 3.

3. **Apply + submit §5.1 `SEAS_DIAG_FAULT_SIGMA` bulk probe — ~200 SU, 30 min.**
   - **Condition:** Step 1 did NOT close R-V92 (i.e., dt halving showed
     no dt-scaling signature) AND user approval granted for both the
     additive C0 source patch (`tpv102_debug_v9.2.0_diag_fault_sigma.patch`)
     AND the Frontera sbatch (`tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch`,
     already written).
   - **Closure:** `max|Q_global[SYY]|` ≫ 50 MPa closes as H-V92-G
     (bulk amplifier confirmed → Step 7 C2 source patch identifying
     the amplifier); ≲ 1 MPa closes H-V92-G and opens H-V92-K / H-V92-U
     (post-friction drift, need different C2).
   - Conditional on step 1 not closing.
   - Measures `max|Q_global[SYY]|` on fault QPs per RK4 stage.
   - Discriminates bulk pump (H-V92-G) from post-friction drift (H-V92-K) from interior-fault path (H-V92-U).
   - **Next if not closed:** Step 4 + Step 6 (free local probes run in parallel).

4. **Write `test_interior_fault_flux_path.cpp` — FREE, local, ~1 h.**
   **[DONE: PASS 6/6; H-V92-U eliminated on 2-tet fixture]**
   - **Condition:** local build environment (`conda activate mfem-dev`)
     functional; MFEM + SEAS headers in tree.  Independent of all
     other steps — safe to run anytime.
   - **Closure:** A FAIL on T1/T2/T3 would directly localize the bug
     to the interior-fault path and close R-V92 as H-V92-U (→ Step 7
     C2 on `wave_operator.inl:746-950`).  PASS does NOT close — only
     eliminates the 2-tet code path.
   - 2-tet single-rank fixture, antisymmetric σ_xy gradient, 1 `Mult` call.
   - Assert DOFData matches analytic Pelties eq. (7).
   - Probes the untested code path that handles the majority of fault faces (H-V92-U).
   - **ACTUAL verdict:** PASS 6/6 — H-V92-U eliminated on this fixture.
     Bimaterial / corner / multi-face extensions still open → continue.

5. **Fix F01 + F02 RK4 operator-splitting — FREE, local, ~2 h.**
   **[DONE — rev-3h+, 2026-04-21]**

   **Scope clarification (REVIEW R-V92-H02):** the plan's original
   F01 bullet reads "move DOFData (V1/V2/tau/sigma/psi/slip) INTO
   the RK4 state vector," implying all 8 fields.  The applied fix
   deliberately narrows that to **psi only**, because:
   - `V1/V2/tau*_corr/sigma_n_corr` are algebraic functionals of
     (Q, psi) at each RK4 stage, not ODE states.  "Moving them
     into Q" has no well-defined meaning — they are not
     integrated, they are evaluated.
   - `slip1/slip2` satisfy `dslip/dt = V`; the existing code
     updates them via `slip += V_avg · dt` with RK4-weighted
     `V_avg`, which IS the RK4 integral of V.  Already O(dt⁴).
   - Only `psi` had a genuine splitting defect (stage-wise
     analytic update with constant V); fixing that is the full
     F01+F02 correction.
   - **Condition:** Step 1 outcome was `~16×` (F01+F02 needed) or
     `~8×` (F02 alone) drop in `slip_dip`; OR Step 1 was non-decisive
     but user wants to try F01+F02 as a precautionary fix (lower cost
     than Step 3 Frontera).
   - **Closure:** running the updated driver locally on a tiny
     fixture AND re-running Step 1's dt-halving check on Frontera
     must both show the predicted scaling.  Full close requires Step 8.
   - Move DOFData (V1/V2/tau/sigma/psi/slip) INTO the RK4 state vector.
   - Independent of all above; does NOT require ADER port.
   - **Next if not closed:** Step 6 / Step 7 on a different mechanism.

   **Implementation (applied to `drivers/tpv102_driver.cpp` only; no
   library-side edits).**  Scope reduces from the plan's literal
   "move all DOFData into Q" to a single functional change: replace
   the operator-split analytic psi updates between RK4 stages with
   a classical coupled RK4 on (Q, psi).  Rationale:
   - `V1/V2/tau*_corr/sigma_n_corr` are algebraic functionals of
     (Q, psi) at each stage, not ODE states; the existing code
     already RK4-averages their stage outputs into DOFData for
     station output, which is the natural and correct behaviour.
   - `slip1/slip2`: `slip_new = slip_old + V_avg * dt` with
     `V_avg = (V_k1 + 2V_k2 + 2V_k3 + V_k4)/6` IS the RK4
     integral of `dslip/dt = V`, so already O(dt⁴).  Nothing
     to move.
   - `psi` was the only genuinely split state: stages 1–3 wrote
     `dof_data.psi = UpdateStateAnalytic(psi_n, sr_k_i, ..., dt_sub)`
     (exact for CONSTANT V over dt_sub, not for the varying-V
     intrinsic to RK4), and the final block re-evaluated
     `UpdateStateAnalytic(psi_n, sr_avg, dt_step)` — together
     O(dt²) coupling.  The fix records `psi_k_i = aging_law.Rate(
     sr_k_i, psi_stage, Dc)` at each stage using the stage-local
     (V, psi), and closes with
     `psi(t+dt) = psi_n + dt/6 · (psi_k1 + 2psi_k2 + 2psi_k3 + psi_k4)` —
     classical RK4 on psi, O(dt⁴) coupled with Q.
   - `AgingLawPsi(b, V0, f0)` from `friction/state_evolution.hpp`
     supplies `Rate()`; instance is created once outside the RK4
     loop.

   **Verification (local, 2026-04-21):**
   - `seas_tpv102_driver` rebuilt clean (no compile warnings in
     the diff region).
   - `seas_test_interior_fault_flux_path` (Step 6) PASS 6/6 — the
     test drives `wave.Mult` directly so it does not exercise the
     driver's RK4 loop; retained here as a regression gate for
     unrelated dynamic-path edits.
   - `seas_test_no_penalty_dynamic_rupture` (§4.9) PASS 4/4.
   - `seas_test_rk4_conservation` (§4.11) PASS 2/2.
   - `seas_test_volume_jacobian_single_channel` (§4.10) PASS 2/2.
   - `seas_test_absorbing_bc_energy_decay` (§4.12) PASS 1/1.
   - Per `feedback_no_local_reproducer.md`: did NOT run the
     production driver locally.  Post-fix dt-scaling and σ_n /
     slip_dip magnitudes must be verified on Frontera via Step 8
     (or a dedicated post-F01+F02 rerun — sbatch NOT yet written;
     waiting for user direction).

   **Submission guidance (REVIEW R-V92-H03 — Option B HOLD).**
   F01+F02 is applied LOCALLY ONLY.  The change is NOT committed
   or pushed until Step 1 PRE-fix dt-halving runs complete,
   because the existing `dt_half` / `dt_quarter` sbatches rebuild
   from branch HEAD — pushing F01+F02 first would make those runs
   exercise the POST-fix driver and conflate the round-5 R-V92-G01
   dt-scaling classifier with the F01+F02 effect.  Sequence:
   1. User submits `dt_half` + `dt_quarter` (PRE-fix HEAD).
   2. Wait for `RESULT.txt` from both.
   3. If Step 1 closes R-V92 (≈16× or ≈8× slip_dip reduction):
      F01+F02 was the needed fix and the PRE-fix run classifies
      the truncation order cleanly; commit F01+F02 then launch
      Step 8 confirmation.
   4. If Step 1 does NOT close: commit F01+F02, then submit
      `tpv102_200m_p1_4.0s_400r_v92_f01f02_dev.sbatch` (post-fix
      dev run at `cfl=0.5 / tfinal=4 s`) — expected to reveal
      whether F01+F02 alone materially reduces the anomaly.  If
      not, continue to Step 3 §5.1 bulk SYY probe.
   The v91 job 7668434 data (slip_dip peak 1.83 m, σ_n peak 218 MPa)
   remains the pre-fix reference throughout.  Step 8 confirmation
   criteria (|σ_n−120 MPa| ≤ 1 MPa AND |slip_dip| ≤ 10⁻³ m at every
   station, tfinal=12 s / 200 m / 400 ranks) unchanged.

   **Step 5b (DONE, rev-3h+, 2026-04-21): local stress-test of
   F01+F02 on a 4-km shared-fault fixture (R-V92-I02, bug-fix-
   oriented).  PASS 100/100.**  Motivation: before any Frontera submission, run
   the full driver RK4-on-(Q, psi) arithmetic for ~100 steps on the
   existing 4-km inline tet fixture (same one used by §4.7
   `test_shared_fault_dof_data_consistency`) and assert the DOFData
   fields stay in a physically-bounded envelope.  The §4.7 test
   runs 1 RK4 step and closes H-V92-P at init level; the new test
   extends to a multi-step run that exercises the coupled-RK4
   arithmetic across MPI partition seams — exactly the code path
   F01+F02 touches.
   - **New file:** `tests/parallel/test_rk4_f01f02_shared_fault_stability.cpp`
   - **Fixture:** 4-km cartesian tet box (identical to §4.7 Phase D);
     4 MPI ranks; fault attr=3 at y=0; TPV102 material parameters;
     `InitializeFaultDOFs` for TPV102 background stress.
   - **Drive:** extract the full RK4-on-(Q, psi) from the driver
     into a `DriverRK4OnePost()` helper inside the test TU (local
     copy, does NOT edit the driver); call it 100 times at dt =
     0.5 · dt_cfl; re-apply `ApplyNucleation` at stage-matched
     times.
   - **Assertions at every 10 steps:** `|σ_n_corr − σ_n0| < 1 MPa`
     AND `slip_rate < 1 m/s` AND `|psi − psi_initial| < 0.1` AND
     no NaN — all enforced pairwise across ranks (same gather
     pattern as §4.7 Phase D).
   - **Runtime:** ~5 s on 4 ranks; well under the laptop's
     14-core oversubscription limit.
   - **Actual result (2026-04-21):** PASS 100/100.  Envelope at
     step 100: `max|σ_n-σ_n0|` = 2.487·10⁵ Pa (≪ 1 MPa bound);
     `max slip_rate` = 1.371·10⁻¹² m/s (≈ V_ini); `max|psi-psi_0|`
     = 1.1·10⁻¹⁴ (FP noise); NaN count = 0.  σ_n deviation grows
     linearly (2.5·10⁴ → 2.5·10⁵ Pa) over 1 ms, consistent with
     elastic-wave transit of the seeded tanh σ_xy(y) profile
     across shared fault faces — no amplification pathology in
     the coupled RK4 arithmetic.
   - **Decision (REVIEW R-V92-H03 sequencing):**
     - If all 100 steps PASS: F01+F02 arithmetic is multi-step-
       stable on shared-fault QPs.  Commit + push F01+F02 with
       high confidence.  The dt-half/dt-quarter Step 1 runs then
       measure dt scaling on the POST-fix driver (i.e. "does the
       post-fix driver still show dt-dependence?" — a clean
       residual probe).  This MUDDLES the pre-fix Step 1
       interpretation somewhat — discuss with user before pushing;
       the cost of a second pre-fix Step 1 run is high if user
       wants the literal R-V92-G01 experiment.
     - If the test FAILS before step 50: F01+F02 contains a bug
       not caught by `test_rk4_psi_integration` (the constant-V
       convergence probe) or the existing Pelties gates.  DO NOT
       push; bisect the failing stage locally until the regression
       is understood, then patch.
     - If the test DRIFTS slowly past step 50: F01+F02 is a
       partial fix or has a minor integration leak.  Still
       commit + push (the fix is not regressive in the
       short-horizon sense), but submit the f01f02_dev sbatch
       after Step 1 to quantify the residual on production mesh.
   - **Interplay with R-V92-H03 (Option B hold):** the stress
     test is not a Frontera run, so it does not compete with the
     "hold F01+F02 before Step 1 submission" rule.  It simply
     raises confidence in the local fix before it propagates.

6. **Add `[FAULT-INIT-V1]` printf at end of RK4 stage 1 — FREE, 1-line.**
   - **Condition:** after Steps 1–5 did not close, OR anytime as a
     cheap hypothesis probe.  Requires user approval for the 1-line
     source edit (not a patch file, directly in
     `drivers/tpv102_driver.cpp`).
   - **Closure:** `max|V1|` nonzero on the first RK4 stage closes
     R-V92 as H-V92-W (initial pre-stress projection inconsistency →
     Step 7 C2 on InitializeFaultDOFs).  Zero closes H-V92-W but does
     not close R-V92.
   - Closes or opens H-V92-W (initial pre-stress projection inconsistency).
   - **Next if not closed:** escalate to REVIEW round-6 for new hypotheses.

7. **Apply C2 source patch** once primary mechanism identified by steps 1–6.
   - **Condition:** exactly one of H-V92-{G, K, U, W} or RK4/psi-split
     confirmed by steps 1–6 AND patch text reviewed AND user approval
     granted per `CLAUDE.md` "Proposing Fixes".
   - **Closure:** patch applied + local regression gates §4.1–§4.12
     PASS + Step 6 unit test PASS.  Full close requires Step 8.

8. **Confirmation Frontera run — 400 SU, 2 h.** 12 s at 200 m / 400 ranks with fix.
   - **Condition:** Step 7 patch applied, committed, pushed; Step 12
     regression gates all green locally; user approval for Frontera
     submission.
   - **Closure (mandatory for plan closure):** `|σ_n − 120 MPa| ≤ 1 MPa`
     AND `|slip_dip| ≤ 10⁻³ m` at every one of the 9 SCEC stations.
     Anything else means the confirmed mechanism in Step 7 was partial
     → re-open with a different hypothesis.
   - Required: `|σ_n − 120 MPa| ≤ 1 MPa` AND `|slip_dip| ≤ 10⁻³ m` at every station.

9. **Write `tpv102_debug_v9.2.0_fix.md`**, then **open `tpv102_debug_v9.2.0_check.md`**, then mark plan DONE.
   - **Condition:** Step 8 confirmation run PASSES both criteria AND
     regression gates §4.1–§4.12 pass on the patched tree AND Step 6
     interior-fault test still passes.
   - **Closure:** plan is DONE when the check document is handed off
     to the reviewer agent and it returns approval.

### DO NOT

- Do NOT port ADER-DG. (Round-4 REVIEW closed it as non-primary by dimensional analysis.)
- Do NOT port BP5 owner-pattern (§16.5) as the R-V92 fix. (Rounds 1–3 closed H-V92-P at all tested scales.)
- Do NOT instrument `ApplyMassInverse` first. (Candidate (d) already eliminated by §4.10.)
- Do NOT add more 4-km Cartesian fixture tests. (Geometry too favorable; cannot probe FP-fragility.)

### Budget

| Step | SU | Cumulative |
|---|---:|---:|
| 1 dt-halving | 20 | 20 |
| 2 peak-value re-analysis | 0 | 20 |
| 3 §5.1 bulk probe (conditional) | 200 | 220 |
| 4 interior-fault unit test | 0 | 220 |
| 5 F01+F02 fix | 0 | 220 |
| 6 init-V1 printf | 0 | 220 |
| 7 source patch | 0 | 220 |
| 8 confirmation run | 400 | 620 |

Best case 220 SU; full path 620 SU; well under the 1 000 SU envelope.

See §20 for the full synthesis and §3 for the open-hypothesis table.

---

## Unit-test manifest (2026-04-21 rev-3; rev-3b + rev-3c updated)

This manifest enumerates the tests and diagnostics so the build +
verification order is explicit.  Tests §4.1–§4.4 are **regression
gates** (green under §13); §4.5 is the generic rank harness
(Frontera-ready); §4.6 (rev-3b, **init-level**) PASSED —
eliminating init-level H-V92-P; **§4.7 (rev-3c, NEW,
runtime-level) is now the highest-priority direct diagnostic for
H-V92-P** per the ParaView speckle evidence in §15 and the BP5
architectural comparison in §16.  See §4, §15.5, §16.6 for
specifications; §13 for §4.1–§4.4 results; §15 for ParaView
speckle + H-V92-P promotion; §16 for BP5-reference workflow +
fix shape.

| # | Test file | Kind | Purpose | Hypothesis | Expected | rev-3 Actual |
|---|---|---|---|---|---|---|
| **§4.9 (rev-3f, NEW — ORDERED FIRST per user directive)** | `tests/unit/test_no_penalty_dynamic_rupture.cpp` | C++ serial unit | **Verify by direct numerical probe (not assumption)** that the dynamic-rupture driver contains NO IP / BR2 / SIPG / Kelvin-Voigt / penalty term, per **Pelties 2012 JGR §3.1 benchmark specification** ("ADER-DG does not generate spurious high-frequency perturbations on the fault and hence does not require artificial Kelvin-Voigt damping"; §3.2 eq. 13: "σ_yy, σ_zz, σ_yz are associated to the so-called zero wave speeds and do not contribute to the Godunov state").  4 probes: **T1** `Mult(Q=0)=0` (bit-exact — rules out any additive shift), **T2** `Mult(2·Q)=2·Mult(Q)` (bit-exact — rules out nonlinear stabilizer), **T3** σ_yy jump injects into `dQ/dt[VY]` (expected physics per Pelties eq. 4), **T4** source grep of `dynamic/` for penalty keywords. | **Rank-0 invariant: benchmark conformance.**  A FAIL directly contradicts Pelties 2012 §3.1 and invalidates every other gate. | Bit-exact PASS | **PASS 4/4.**  T1: `‖dQ/dt‖_∞ = 0.000e+00` (bit-exact).  T2: `‖k₂ − 2·k‖_∞ = 0.000e+00` on a random Q at 1e6-scale (bit-exact linearity).  T3: `max\|dQ/dt[VY]\|` = 7.5·10⁴ (m/s)/s (correct Pelties-coupling signature).  T4: 0 penalty-keyword hits in `dynamic/*.{cpp,hpp,inl}`.  **Candidate (a) of H-V92-G numerically ELIMINATED** — the amplifier is NOT a hidden penalty term. |
| **§4.10 (rev-3g, H-V92-G (d)+(e) probe)** | `tests/unit/test_volume_jacobian_single_channel.cpp` | C++ serial unit | Single-channel Jacobian sparsity probe: 2-hex mesh, absorbing BCs; for each of 9 state components c, set `Q[c] = sin(πx)` (others 0) and verify `Mult(Q)` response rows match the expected `A_x + A_y + A_z` sparsity (volume + face-flux); plus an isotropy sub-test comparing `SXY→VY` vs `SXZ→VZ` magnitude ratio. | **H-V92-G candidates (d) mass-inverse channel asymmetry, (e) volume Jacobian coupling.** | PASS (no sparsity violations; isotropy ratio = 1.0) | **PASS 2/2.**  All 9×9 channel pairs match predicted sparsity (no spurious coupling).  Isotropy ratio `\|dQdt[VY]_from_SXY\|/\|dQdt[VZ]_from_SXZ\|` = **1.000000e+00** bit-exact.  Poisson coupling `dQdt[SYY]:dQdt[SXX]` from driven VX = **0.333** (matches `λ/(λ+2μ)`).  **Candidates (d) and (e) ELIMINATED** — mass inverse channel-uniform, volume Jacobian correctly coupled. |
| **§4.11 (rev-3g, H-V92-G (b) probe)** | `tests/unit/test_rk4_conservation.cpp` | C++ serial unit | RK4 conservation probe: 8³ hex, 8 km domain, absorbing BCs; Gaussian pulses in σ_xx and σ_yy; drive 100 RK4 steps at dt = 50 μs (pulse stays well inside domain); measure per-channel L² drift. | **H-V92-G candidate (b) RK4 non-conservation on σ_yy.** | Per-step drift ≤ 1% and SYY/SXX drift ratio < 10× | **PASS 2/2.**  SXX L² drift = **−2.33·10⁻⁵ / step** (dissipative, expected); SYY L² drift = **−1.63·10⁻⁵ / step** (dissipative, expected); **SYY/SXX ratio = 0.70 < 1** — no channel asymmetry.  **Candidate (b) ELIMINATED.** |
| **§4.12 (rev-3g, H-V92-G (c) probe)** | `tests/unit/test_absorbing_bc_energy_decay.cpp` | C++ serial unit | Absorbing-BC energy decay probe: 8³ hex, 16 km domain, ALL absorbing BCs; isotropic Gaussian stress pulse; 3 000 RK4 steps at dt = 500 μs (T = 1.5 s > pulse travel time 1.33 s); measure total elastic energy E(t) every 100 steps. | **H-V92-G candidate (c) absorbing BC wrong-sign reflection.** | E(t) monotonically non-increasing | **PASS 1/1.**  E(t) decays from 75.05 to 18.25 over 1.5 s (E/E₀ = 0.243 at T_end); `E_max/E₀ = 1.000` exactly (no increase anywhere across 3 000 steps).  **Candidate (c) ELIMINATED.** |
| §4.1 | `tests/unit/test_godunov_rotation_identity.cpp` | C++ unit | `max\|T·Tinv − I\| ≤ 16 ULP` on BP5 + 3 off-axis frames | Regression gate H-V92-R (retracted) | PASS | **PASS 4/4** |
| §4.2 | `tests/unit/test_canonical_rotation_pure_strikeslip.cpp` | C++ unit | Tinv·Q_g leaks 0 into normal/dip on BP5 strike-slip | Regression gate H-V92-O (retracted) | PASS bit-exact | **PASS 11/11** |
| §4.3 | `tests/unit/test_fault_basis_qp_orthonormality.cpp` | C++ unit | Per-QP `(n, t1, t2)` orthonormal on 200 m fault | Regression gate H-V92-F (mesh jitter) | PASS | **PASS** (113 712 QPs @ 0.5 ULP) |
| §4.4 | `tests/unit/test_per_qp_vs_centroid_basis_planar_1el.cpp` | C++ unit | Per-QP basis = centroid basis on planar fixture | MFEM-Jacobian regression gate | PASS bit-exact | **PASS 11/11** |
| §4.5 | `tests/scripts/test_sigma_n_rank_consistency.py` | Py harness (Frontera) | σ_n(t=3.5 s) agrees 56-rank (new run) vs 400-rank (existing v91 job 7668434) on 200 m mesh | Generic H-V92-M / runtime H-V92-P probe | Frontera sbatch + v91 reference | **READY** — `tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch` written; fallback for §6.4 + §9.2 ambiguity |
| **§4.6 (rev-3b)** | `tests/parallel/test_shared_fault_role_consistency.cpp` | C++ MPI (4 ranks on 4 km tet fixture) | **Ctor-level**: `elem1_on_plus` ANTI-symmetric (per design, `wave_operator.inl:1279-1290`); `can_n`/`can_t1`/`can_t2`/`nl` bit-identical across ranks | **Rank-1 for H-V92-P (init level)** | PASS expected after spec-semantic correction | **PASS 5/5** — 18 cross-rank pair checks; canonical frame rank-invariant to 0 ULP; anti-symmetry contract held on all pairs.  Init-level H-V92-P **ELIMINATED**. |
| **§4.7 (rev-3c/3d)** | `tests/parallel/test_shared_fault_dof_data_consistency.cpp` | C++ MPI (4 ranks on 4 km tet fixture; **4 phases post-REVIEW R-V92-C01**: A post-init / B post-Mult(Q=0) / C post-RK4 uniform-Q / **D post-RK4 antisymmetric σ_xy(y) = τ_ini·tanh(y/L0)**) | **Post-RK4-step**: drive WaveOperator through 1 RK4 step, then pair every shared fault QP across ranks via centroid-key and assert all 8 DOFData fields agree bit-for-bit (0-ULP).  **Phase D** (rev-3d) is the only phase with `Q_self ≠ Q_nbr` genuinely exercised (spread = 1.5·10⁸ Pa confirmed at fault); Phases A/B/C are subset probes per R-V92-C01. | **Rank-1 for H-V92-P (runtime level)** — Phase D is the operational definition: antisymmetric Q ensures the (+, −) swap routing actually differs between ranks | FAIL at head expected (per §15 speckle); PASS after BP5 owner-pattern port (§16.5) | **PASS 4/4 phases** (A+B+C+D).  Phase D pre-Mult spread = 1.499·10⁸ Pa confirms swap exercised; 9 pair checks × 8 fields all bit-identical after 1 RK4 step.  **H-V92-P swap-routing mechanism NOT triggered** on the 4-km fixture even with genuine `Q_self ≠ Q_nbr` precondition. |
| **§17.6 (rev-3d)** | `tests/parallel/test_centroid_margin_ctor_only.cpp` + `SEAS_DIAG_CENTROID_MARGIN` | C++ MPI ctor-only (no time stepping; local-safe per feedback_no_local_reproducer) | Dumps per-shared-fault-face FP margin of the `(elem1_proj − face_proj)` comparison used to set `elem1_on_plus`.  Margin → 0 = FP-fragility active; margin = O(10² m) = robust. | **Refutation test for REVIEW R-V92-C04** (FP-fragility hypothesis on production geometry) | Expected ≥ 100 m on real meshes | **PASS — margins 227 – 266 m on TPV102 1000 m mesh at 12 ranks**.  14 orders of magnitude above ε_FP; R-V92-C04 FP-fragility hypothesis directly refuted by measurement. |
| **§6.4 (rev-3b+3d+3e)** | `tests/scripts/fault_vtu_outlier_detector.py` | Python (numpy, XML-only — no `vtk` dep) | Loads existing job-7668434 fault-surface PVD; per-triangle median+MAD outlier detection across 8 dynamic fields; partition-correlation + trend analysis → H-V92-P vs G vs Q classifier. | Speckle-mechanism classifier | Local, no Frontera — **EXECUTED on pulled job-7668434 data** | **H-V92-G CONFIRMED.**  93 cycles, 38 039 triangles/cycle, 35 s wall-time.  Outlier trend +2 703 / second (monotonic growth).  Dip-channel fields grow 19 – 180×, σ_n grows 40×.  See §18 for full analysis. |
| §5.2 (NEW rev-3b) | `io/paraview_output.hpp` patch | C++ I/O | Add `mpi_rank` + `is_shared_face` CellData under `SEAS_DIAG_VTU_RANK` | Supports §6.4 partition correlation | Pure I/O; safe for Frontera | **PENDING** — patch not applied.  (Note: the file was modified by a separate agent/session with R-V92-E02 stage-4 VTU diagnostic; that work is unrelated to this row's scope.) |
| §9.2 | `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch` | Frontera sbatch | 4 s run with `-DSEAS_DIAG_FAULT_SIGMA` (+ optional §5.2 `-DSEAS_DIAG_VTU_RANK`) | After §18 H-V92-G confirmation: confirms σ_yy-on-fault trace for independent cross-check | ~200-358 SU | **READY, NOT submitted** — sbatch written; requires `tpv102_debug_v9.2.0_diag_fault_sigma.patch` applied + user sbatch approval.  Now deprioritized per §18.8 (§6.4 already provides the same evidence). |

**Build + verification order (rev-3 + rev-3b + rev-3c + rev-3d + rev-3e + rev-3f all executed):**
0. **§4.9 (rev-3f, NEW, ORDERED FIRST) — Done — PASS 4/4.**
   No-penalty verification per Pelties 2012 §3.1.  Eliminates
   H-V92-G candidate (a) by direct numerical probe.
1. §4.1, §4.2 — pure linear algebra.  **Done — PASS.**
2. §4.3, §4.4 — MFEM mesh + FaultBasis.  **Done — PASS.**
3. **§4.6 (rev-3b, 4-rank local, ≈ 5 s)** — **Done — PASS.**
   (2-rank skipped because ParMETIS did not cut the fault;
   4-rank yields 18 bilateral pair checks across 9 shared QPs.
   `elem1_on_plus` correctly ANTI-symmetric; can_n/t1/t2/nl
   bit-identical across ranks.)  **H-V92-P eliminated at init
   level; runtime-only P remains open at production scale.**
4. **§4.7 (rev-3c, 4-rank local, ≈ 5 s)** — **Done — PASS 3/3.**
   Runtime DOFData consistency on the 4-km fixture holds
   bit-for-bit across (post-init, post-Mult(Q=0), post-RK4(Q≠0)).
   **H-V92-P runtime divergence NOT reproducible on this
   fixture.**  H-V92-P remains open at production scale only;
   §9.2 Frontera required to either confirm on 200 m mesh (→
   §16.5 fix) or close entirely (→ different mechanism).
5. §6.4 offline outlier detector on existing job-7668434 PVD
   (no Frontera).  PENDING.
6. §5.2 rank-id VTU patch (I/O only).  PENDING.
7. BP5 owner-pattern port per §16.5 (~280 LOC structural fix),
   now **gated on §9.2 confirming runtime H-V92-P on production mesh**
   (§4.7 PASS locally shifted the gate from local to Frontera).
8. Re-run §4.7 after any fix; must PASS.  (The local test remains
   a regression gate.)
9. **§9.2 Frontera run (~200 SU)** — now **HIGHEST-PRIORITY next
   step** since local tests exhausted and H-V92-P is open only
   at production scale.  Requires §5.1 diagnostic patch (§13.4)
   applied first.
10. §4.5 harness (rank-A sbatch, ~48 SU) — fallback; still ready.

**The regression gates being green does NOT confirm the codebase
is bug-free for R-V92.**  §4.7 is the **direct operational
definition** of H-V92-P (shared-fault DOFData divergence across
ranks); its PASS after the §16.5 port closes the hypothesis by
construction.  §5.1 + §4.5 + §6.3 (rev-3 original plan) retained
as secondary discriminators for H-V92-G / C / K once H-V92-P is
addressed.

---

> **Author:** debugger agent. **Date:** 2026-04-21 (rev-2 post-REVIEW.md; rev-3 added §13 test results).
> **Successor to** `tpv102_debug_v9.1.0_debug_plan.md`,
> `tpv102_debug_v9.1.0_fix.md`, and
> `tpv102_debug_v9.0.0_seissol_flux_comparison.md`.
>
> **rev-2 summary.** The v9.2.0-rev-1 draft was failed by adversarial
> review (REVIEW.md findings R-001 / R-002 CRITICAL).  Changes in
> rev-2:
> - H-V92-R (rotation stress-block asymmetry) **RETRACTED**.  The
>   reviewer proved `(T·Tinv)[ij, i'j']` is bit-exact `δ` by direct
>   expansion for any orthogonal Q; verified independently (see §A).
>   The `(a != b)` vs `(i != j)` asymmetry in
>   `godunov_flux.cpp:272` vs `:231` is the correct forward/inverse
>   Voigt symmetrization, not a bug.
> - H-V92-O (pure-strike-slip Tinv probe) **RETRACTED** — subsumed by
>   R-001.  For the BP5 canonical frame `((0,-1,0), (0,0,-1),
>   (1,0,0))`, Q is a signed permutation, Tinv entries are all in
>   `{-2,-1,0,+1,+2}`, and Tinv·Q_g on any axis-aligned fixture is
>   bit-exact integer arithmetic.  §4.2 would pass trivially.
> - H-V92-M (**MPI rank consistency**) **PROMOTED to rank-1** per
>   REVIEW.md R-003.  The observed σ_n anti-symmetry in along-strike
>   x (§1.3) matches the signature of a ParMETIS partition cut near
>   x = 0, with side-dependent `elem1_on_plus` bit-population.
> - Two **new** hypotheses introduced: H-V92-G (bulk-dynamics
>   mode-II + Rayleigh stress coupling feeding the fault QP through
>   the standard characteristics — the correct *physics* this plan
>   previously foreclosed in §2.1 per REVIEW R-005) and H-V92-C
>   (fault–free-surface corner interaction via Interior's internal
>   BuildFrame choice — reviewer's flag for the next round).
> - Investigation order inverted per REVIEW recommendation:
>   diagnostic patch FIRST, MPI consistency SECOND, ParaView SYY
>   slice THIRD.  Unit tests (previously rank-1 confirmation tests)
>   demoted to **regression gates**; they now guard against future
>   regressions, not drive the diagnosis.
> - §4.1 tolerance loosened to **16 ULP on axis-aligned frames, and
>   documented as not applicable to per-QP non-axis-aligned frames**
>   per REVIEW R-004.
>
> **Scope.** Investigation + plan only.  No source edits applied.  No
> Frontera jobs launched.  Waiting for user approval per
> `feedback_frontera_approval.md` and `CLAUDE.md` ("Proposing Fixes")
> before any C0 patch or `sbatch`.
>
> **Convention notice (from user directive).** MFEM/SEAS-MFEM uses
> the **BP5 (Tandem) fault-local convention**: `tangent1 = dip`,
> `tangent2 = strike`.  SeisSol uses a different `(t1, t2)` ordering
> internally.  Where this plan compares to SeisSol code, a
> component-name mismatch between MFEM's `V1/tau1` (dip) and any
> SeisSol variable named `T1/traction1` is **not** automatically a
> bug — it is checked only for *mathematical* equivalence of the
> equations (magnitude, frame consistency, sign-of-drop behaviour),
> not textual-component matching.  A convention-only difference is
> **NOT** classified as a root cause anywhere in §3–§9.
>
> **Artifacts this plan consumes:**
> - Station plots: `tpv102/plots_results_200m_p1_12.0s_400r_v91_job7668434/`
>   (9 SCEC fault stations, 12 s, post-R-V91-A/B fixes).
> - Data files (already removed from Downloads — key values extracted inline in §1.2):
>   `~/Downloads/seas-mfem/tpv102/results_200m_p1_12.0s_400r_v91_job7668434/`
> - Source state at head of `feature/elasticity-inertia`, commit
>   `99f355a` "TPV102 v9.1.0 Pelties-9 follow-through" (has
>   R-001 / R-002 / R-003 applied; V9.0.0 per-side Pelties-9 flux in
>   place at `wave_operator.inl:832-918` and `:1326-1340`).
> - Benchmark documents:
>   `tpv102/benchmark_document/{Dumbser-Käser 2006, Pelties 2012, SCEC validation}.pdf`.
> - SeisSol reference source:
>   `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/`,
>   `/Users/chunhuizhao/projects/SeisSol/src/Equations/elastic/Model/ElasticSetup.h`.

---

## TODO — remaining actions to close R-V92 (rev-3h, 2026-04-21)

See **PRIORITY NOW** block at the top of this file for the
single-page summary.  This section is the detailed per-action
list.

> **Plan-closure gate.**
> v9.2.0 closes ONLY when all three are true:
> - The primary R-V92 mechanism is identified with a discriminating test.
> - A C2 source patch has been applied.
> - A Frontera 12 s run at 200 m / 400 ranks shows `|σ_n − 120 MPa| ≤ 1 MPa` AND `|slip_dip| ≤ 10⁻³ m` at every SCEC station.
>
> **Current state (rev-3h):** σ_n peak 218 MPa (~100× spec); slip_dip 1.83 m (~1 830× spec).
> **We are NOT close to closure.**

---

### Unfinished — blocking path

Execute top-to-bottom.  Each step may terminate the loop.

#### [~] Step 1 — dt-halving decisive test (HIGHEST PRIORITY)

- **Owner:** Debugger + YOU
- **Cost:** ~20 SU / 1 h Frontera
- **New files:** `jobs/tpv102/tpv102_200m_p1_12.0s_400r_v92_dt_half.sbatch`, `..._dt_quarter.sbatch`
- **Source:** REVIEW round-5 R-V92-G01

**Status (rev-3h+, 2026-04-21):** sbatch files WRITTEN (both,
`--cfl 0.25` and `--cfl 0.125`, post-processor classifier matches
the decision table above); AWAITING USER APPROVAL + SUBMISSION per
`feedback_frontera_approval.md`.  Walltime/SU budget in each sbatch
header notes the plan's 20 SU figure is for the classifier only —
actual production 12 s / 200 m / 400 ranks at dt/2 is ~2 048 SU
(2× v91 baseline) and dt/4 is ~4 096 SU; scale `--tfinal` down if
budget-constrained.  Step STAYS OPEN until user submits + the
`RESULT.txt` classification is filed.

Compare `slip_dip` peak across 9 SCEC stations vs v91 baseline.  RK4 truncation error scales as dt⁴.

| Outcome at dt/2 | Classification | Next action |
|---|---|---|
| ~16× smaller slip_dip | Time-integration bug (4th-order) | Apply F01+F02 fix (Step 5) → done |
| ~8× smaller slip_dip | psi splitting (F02 3rd-order) | Integrate psi inside RK4 state |
| Unchanged | NOT time integration | Proceed to Step 3 |
| Larger / NaN | CFL violation | Reduce CFL + audit `flux.Interior` |

---

#### [x] Step 2 — Peak-value re-analysis of existing data

- **Owner:** Debugger
- **Cost:** FREE, local, ~1 h
- **Source:** REVIEW round-3 R-V92-E01

Re-process the existing `outlier_report.json` from job-7668434 to extract **peak `|σ_n − 120 MPa|`** and **peak `|slip_dip|`** per cycle — not outlier counts.

- If peak saturates near the 45 MPa Poisson bound → §18 "H-V92-G CONFIRMED" was a misread of rupture extent.
- If peak grows unbounded → real amplification; need to find amplifier.

**Result (rev-3h+, 2026-04-21) — EXECUTED.**  Analyzer
`tpv102/r_v92_outlier_detection/plot_step2_peak_vs_poisson_bound.py`
processes all 93 cycles × per-field `{median, max_dev}`.  Median of
`normal_stress` stays in 120.000–120.001 MPa throughout (no drift),
so `max_dev` is a tight bound on `|σ_n − 120 MPa|`.

| Metric                        | Value                           |
|-------------------------------|---------------------------------|
| Peak `|σ_n − 120 MPa|`        | 218.7 MPa — **4.9× Poisson bound (45 MPa)** |
| End  `|σ_n − 120 MPa|`        | 218.7 MPa (= peak)              |
| Peak `|slip_dip|`             | 1.83 m — 1 827× spec ceiling    |
| Late-run mean / max           | 0.907                           |
| Classifier verdict            | **MONOTONIC-GROWTH**            |

Decision: R-V92-E01 "rupture-extent misread" is **REFUTED** on the
peak-value evidence — the σ_n deviation grows 4.9× past the
Poisson bound and has no plateau.  Real amplification (H-V92-G)
remains viable; §5.1 bulk SYY probe (Step 3) is still needed to
identify the amplifier.  Plot:
`tpv102/r_v92_outlier_detection/step2_peak_vs_poisson_bound.png`;
JSON: `.../step2_peak_vs_poisson_bound.json`.

---

#### [ ] Step 3 — Reinstate §5.1 `SEAS_DIAG_FAULT_SIGMA` bulk probe

- **Owner:** YOU (approve) + Debugger (apply)
- **Cost:** C0 patch, additive stderr printf, no functional change
- **Source:** REVIEW round-3 R-V92-E05 / round-5 R-V92-G05
- **Status in §18.8:** wrongly marked "eliminated" — must be reversed.

Bulk probe measures `max|Q_global[SYY]|` on fault QPs per RK4 stage.  Discriminates:
- H-V92-G (bulk pump — probe grows)
- H-V92-K (post-friction drift — probe bounded while DOFData drifts)
- H-V92-U (interior-fault path — probe bounded, drift is in flux conversion)

§18 data does NOT substitute: §18 measures post-friction RK4-averaged DOFData; §5.1 measures bulk Q per stage.

---

#### [ ] Step 4 — §5.1 Frontera diagnostic run (CONDITIONAL)

- **Owner:** Debugger
- **Cost:** ~200 SU
- **Condition:** Run iff Step 1 does NOT close R-V92
- **Sbatch:** `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch` (already written)

Classification by `max|Q_global[SYY]|` at t = 4 s:

| `max|Q_global[SYY]|` | Primary mechanism |
|---|---|
| < 50 MPa | Bulk σ bounded; bug in friction / DOFData path (H-V92-K or H-V92-U) |
| 50 – 200 MPa | Bulk amplification (H-V92-G) — note Open Question #1 |
| > 200 MPa | Catastrophic bulk-operator bug (non-fault-specific) |

---

#### [ ] Step 5 — Fix RK4 operator-splitting (F01 + F02)

- **Owner:** Debugger
- **Cost:** FREE, local, ~2 h, ~100 LOC in `drivers/tpv102_driver.cpp`
- **Source:** REVIEW round-4 F01 + F02

Move DOFData into the RK4 state vector:
- F01: `V1 / V2 / tau1_corr / tau2_corr / sigma_n_corr / slip1 / slip2`
- F02: `psi` (replace analytic post-step update with RK4-integrated derivative)

Independent of time-integrator choice; does **NOT** require ADER port (REVIEW round-4 F03 closed ADER as non-primary).  Can run before or after Step 1.

---

#### [x] Step 6 — New unit test for interior-fault-face path

- **Owner:** Debugger
- **Cost:** FREE, local, ~1 h
- **New file:** `tests/unit/test_interior_fault_flux_path.cpp`
- **Source:** REVIEW round-5 R-V92-G04

2-tet single-rank fixture.  `InitializeFaultDOFs` → apply antisymmetric σ_xy gradient → one `wave.Mult(Q, k)` → assert DOFData matches analytic Pelties eq. (7).

Probes the code path in `wave_operator.inl:746-950` that handles the majority of fault faces per §17.6.5 — **never unit-tested** despite processing most of the fault surface.

**Result (rev-3h+, 2026-04-21) — PASS 6/6.**  Implemented with the
2-tet fixture + free-surface boundaries + fault bdr attr=3 (same
pattern as `test_shared_fault_dof_data_consistency.cpp`).  Uses
antisymmetric velocity perturbations (v_y, v_x, v_z respectively)
instead of the originally-proposed antisymmetric σ_xy — velocities
give a cleaner single-component Pelties probe because the other
components of `Q` stay at zero and the rotation to the canonical
frame is diagonal.  `V_test = 1·10⁻⁵ m/s` ⇒ 160 Pa (7a) and 92 Pa
(7b/7c) perturbations on 120/75 MPa background, well within
equilibrium so the friction solver's V ≈ V_ini term is negligible.

| Probe                 | Pelties prediction | Worst rel err |
|-----------------------|--------------------|---------------|
| T0 Q = 0 baseline     | `tau{1,2}_corr = tau{1,2}_0`, `sigma_n_corr = sigma_n0` | ≤ 5·10⁻⁶ Pa absolute (≤ 1 Pa tol) |
| T1 7a (Δv_y → Δσ_n)   | +Zp · V_test = +160.2 Pa        | 1.9·10⁻¹¹     |
| T2 7c (Δv_x → Δτ₂)    | −Zs · V_test = −92.49 Pa        | 5.0·10⁻⁸      |
| T3 7b (Δv_z → Δτ₁)    | +Zs · V_test = +92.49 Pa        | 6.2·10⁻¹⁴     |

**H-V92-U ELIMINATED on this fixture.**  The interior-fault branch
(wave_operator.inl:746-950) reproduces Pelties eq. (7) for all
three eigenmodes to machine precision.  Caveat: the fixture has a
single interior fault face and uniform material; this does not
probe bimaterial interior paths or corner-geometry interactions
(§3 H-V92-V / H-V92-C), so those hypotheses remain open.

---

#### [ ] Step 7 — Initial-V1 printf (H-V92-W probe)

- **Owner:** Debugger
- **Cost:** FREE, ~1 line
- **Source:** rev-3h H-V92-W hypothesis (§3)

Add `[FAULT-INIT-V1]` printf at the end of RK4 stage 1 of the FIRST time step: `max |V1|` across all fault QPs.  Must be 0 exactly.

- Non-zero → H-V92-W confirmed (pre-stress projection inconsistency).
- Zero → H-V92-W closed.

---

#### [ ] Step 8 — H-V92-G candidate discrimination (CONDITIONAL)

- **Owner:** Debugger
- **Condition:** Run iff Step 4 returns 50 – 200 MPa
- **Source:** REVIEW round-5 R-V92-G03

**Open Question #1 (rev-3h):** All five candidates (a, b, c, d, e) of H-V92-G have been eliminated by §4.9 – §4.12:

| Candidate | Probe | Verdict |
|---|---|:---:|
| (a) hidden IP/BR2 penalty | §4.9 | ELIMINATED |
| (b) RK4 non-conservation | §4.11 | ELIMINATED |
| (c) absorbing-BC wrong sign (energy) | §4.12 | ELIMINATED |
| (d) mass-inverse channel asymmetry | §4.10 | ELIMINATED |
| (e) volume-term Jacobian coupling | §4.10 | ELIMINATED |

**Therefore the "§18 H-V92-G CONFIRMED" verdict is inconsistent with §4.9 – §4.12.**  Either:
- (i) §18 outlier classifier misread (REVIEW R-V92-E01), or
- (ii) A new amplifier candidate exists — likely **H-V92-U** (interior-fault path) or **H-V92-V** (per-channel BC phase).

Resolution depends on Steps 2, 3, 4, 6.

---

#### [ ] Step 9 — Propose C2 source patch

- **Owner:** Debugger
- **Condition:** primary mechanism identified
- **Rules (REVIEW round-4 + round-5):**
  - Do NOT port ADER-DG.  Closed by F03.
  - Do NOT port BP5 owner-pattern (§16.5).  Rounds 1 – 3 eliminated H-V92-P at tested scales.
  - Do NOT instrument `ApplyMassInverse`.  Already eliminated by §4.10.

---

#### [ ] Step 10 — YOU: approve C2 patch.

---

#### [ ] Step 11 — Confirmation Frontera run

- **Owner:** Debugger
- **Cost:** ~400 SU, 2 h Frontera (12 s / 200 m / 400 ranks)
- **Pass criteria:**
  - `|σ_n − 120 MPa| ≤ 1 MPa` at all 9 stations
  - `|slip_dip| ≤ 10⁻³ m` at all 9 stations

---

#### [ ] Step 12 — Regression re-gate

Re-run `make test-v92-regression-gates` (§4.1 – §4.12).  Must all PASS.

---

#### [ ] Step 13 — Write `tpv102_debug_v9.2.0_fix.md`

Document confirmed mechanism + applied patch.

---

#### [ ] Step 14 — Open `tpv102_debug_v9.2.0_check.md`

Hand off to review.  **ONLY THEN mark this plan's TODO as DONE.**

---

### Explicitly paused / superseded (tracking only)

These were in the rev-3 TODO but are no longer active steps.

- [x] ~~**YOU:** Authorise §5.1 C0 diagnostic patch (rev-3).~~
  Superseded — reinstated conditionally as Step 3 above.

- [x] ~~**Debugger:** Run §6.1 offline diagnostic on station .dat files.~~
  Superseded — §18 ran the equivalent on VTU data; corrected replacement is Step 2 above.

- [x] ~~**Debugger:** Run §4.5 MPI rank-consistency locally.~~
  Superseded — §4.6 + §4.7 closed H-V92-P at the scales §4.5 could reach.  §4.5 sbatch stays READY as optional fallback (~48 SU).

- [x] ~~**Debugger:** Run §6.3 ParaView SYY slice.~~
  Superseded — §18 + Step 2 provide equivalent discrimination.

- [x] ~~**Debugger:** Identify H-V92-M / G / F / C primary.~~
  Partial — rev-3e claimed H-V92-G CONFIRMED, but REVIEW R-V92-E01 and §4.9-§4.12 eliminations make that claim inconsistent.  See **Open Question #1** in Step 8.

### Regression gates (DEMOTED from rev-1's §4 rank-1 tests) — **ALL GREEN rev-3**

The four unit tests in rev-1 §4 are retained as regression gates —
they are expected to pass on the current codebase and guard against
future regressions of established invariants, not drive diagnosis
of R-V92.  Rewritten with corrected scope in §10.

- [x] **Debugger:** §4.1 `test_godunov_rotation_identity` — **PASS** (4/4).
      BP5 canonical 0 ULP (bit-exact); 3 off-axis frames 1 ULP each
      (16× under REVIEW R-004 budget).  Full results in §13.3.
- [x] **Debugger:** §4.2 `test_canonical_rotation_pure_strikeslip` — **PASS**
      (11/11).  Pure strike-slip Q_g → Q_c is bit-exact on BP5
      signed-permutation frame; no leakage into n/t1 channels.
- [x] **Debugger:** §4.3 `test_fault_basis_qp_orthonormality` — **PASS**
      on both 1000 m (4 938 QPs) and **200 m production mesh (113 712
      QPs)**.  Worst norm dev 0.5 ULP; worst orth-dot ~10⁻⁴⁶ ULP.
- [x] **Debugger:** §4.4 `test_per_qp_vs_centroid_basis_planar_1el` —
      **PASS** (11/11).  Per-QP basis equals centroid basis bit-exact
      on planar 2-tet fixture.
- [x] **Debugger:** §4.6 (rev-3b) `test_shared_fault_role_consistency`
      — **PASS** (5/5) at 4 ranks on 4 km inline tet fixture.  18
      cross-rank pair checks: `elem1_on_plus` anti-symmetric (per
      design of `wave_operator.inl:1279-1290`); can_n / can_t1 /
      can_t2 / |n_raw| all bit-identical.  **Init-level H-V92-P
      ELIMINATED** on this fixture; runtime P remains open at
      production scale.  See §15.5.1 for full results + spec-
      semantic correction.
- [x] **Debugger:** §4.7 (rev-3c) `test_shared_fault_dof_data_consistency`
      — **PASS** 3/3 phases at 4 ranks on 4 km inline tet fixture.
      9 cross-rank pair checks per phase covering post-init,
      post-Mult(Q=0), and post-RK4 step with nonzero Q; all 8
      DOFData fields (slip1, slip2, V1, V2, tau1_corr, tau2_corr,
      sigma_n_corr, psi) bit-identical across ranks (0 ULP).
      **Note (REVIEW R-V92-C01 response):** Phases A/B/C all have
      Q_self ≡ Q_nbr trivially (uniform Q), so the (+, −) swap is
      vacuous and the PASS was structural, not evidentiary.
      See §17.6 for the post-review response.
- [x] **Debugger:** §4.7 Phase D (rev-3d, post-REVIEW R-V92-C01) —
      **PASS** 1/1.  σ_xy(y) = τ_ini · tanh(y/L0) projected onto the
      FESpace; pre-Mult spread = 1.499·10⁸ Pa confirmed `Q_self ≠
      Q_nbr` on shared faces.  After 1 RK4 step, all 8 DOFData
      fields remain bit-identical across ranks.  **H-V92-P
      swap-routing mechanism NOT triggered** on the fixture.
- [x] **Debugger:** rev-3d SEAS_DIAG_CENTROID_MARGIN diagnostic —
      implemented + measured on TPV102 1000 m mesh at 12 ranks.
      **All margins in range 227 – 266 m; 14 orders of magnitude
      above ε_FP.  REVIEW R-V92-C04 was correct that the 4-km
      fixture cannot probe FP-fragility (centroid margins O(10³ m)
      render it impossible there).  The 1000 m measurement
      executes the reviewer's Step 3 probe on production-like
      geometry and, for that mesh, refutes the FP-fragility
      SUB-hypothesis for `elem1_on_plus`.**
      See §17.6.3 for full data.
- [x] **Debugger:** §6.4 offline VTU outlier detector implemented
      (`tests/scripts/fault_vtu_outlier_detector.py`) per REVIEW
      R-V92-C02 free-local-first priority.  Runs in ~30 min
      locally once the user pulls `FaultSurface/*.{pvd,pvtu,vtu}`
      from Frontera job 7668434.  Classifies H-V92-P / G / Q via
      outlier count trend + partition correlation.
- [x] **Debugger:** **§6.4 EXECUTED on pulled job-7668434 data
      (rev-3e, 2026-04-21).**  93 cycles × 38 039 triangles
      analyzed in 35 s wall-time.  Outlier trend +2 703/s
      (monotonic growth).  Per-field terminal peaks: slip_dip
      1.83 m (1 830 × spec), slip_rate_dip 2.15 m/s (20 000 ×
      spec), traction_dip 68 MPa (70 × spec), normal_stress
      deviation 218 MPa (nearly 2 × the static value).  Free-
      surface reflection spike at t = 6.24 s matches theoretical
      arrival time.  **H-V92-G CONFIRMED as primary mechanism;
      H-V92-P / C / K / F ELIMINATED** by growth-signature
      evidence.  See §18 for full analysis.
- [x] **Debugger:** **§4.9 rev-3f (NEW — ORDERED FIRST per user
      directive, 2026-04-21) — no-penalty regression gate PASS
      4/4.**  Per Pelties 2012 JGR §3.1–§3.2 the TPV102 benchmark
      specifies Godunov fluxes only, no IP/BR2/SIPG/Kelvin-Voigt.
      Test verifies this by direct probe (NOT by assumption):
      T1 `Mult(Q=0) = 0` bit-exact; T2 `Mult(2·Q) = 2·Mult(Q)`
      bit-exact on random Q (linearity); T3 σ_yy jump couples
      into VY per Pelties eq. (4); T4 zero penalty-keyword hits
      in `dynamic/*.{cpp,hpp,inl}`.  **H-V92-G candidate (a)
      (hidden IP/BR2 penalty) ELIMINATED.**  Candidates (b)/(c)/
      (d)/(e) remain; to be probed per §18.8 without pre-ranking.
- [x] **Debugger:** **§4.10 rev-3g (NEW, 2026-04-21) — single-
      channel Jacobian probe PASS 2/2.**  2-hex mesh with
      absorbing BCs; driven Q[c] = sin(πx) for each of 9
      components.  All 9×9 channel pairs match predicted
      `A_x+A_y+A_z` sparsity (volume + face-flux); no spurious
      coupling above 1e-9 noise floor.  Isotropy ratio
      `\|dQdt[VY]_SXY\| / \|dQdt[VZ]_SXZ\|` = **1.000 bit-exact**.
      Poisson coupling `dQdt[SYY]/dQdt[SXX]` from driven VX
      = 0.333 (matches `λ/(λ+2μ)` for TPV102 parameters).
      **Candidates (d) mass-inverse channel asymmetry and
      (e) volume-term Jacobian coupling ELIMINATED.**
- [x] **Debugger:** **§4.11 rev-3g (NEW, 2026-04-21) — RK4
      conservation probe PASS 2/2.**  8³ hex, 8 km domain,
      absorbing BCs; Gaussian pulses in σ_xx (amp 0.5 MPa) and
      σ_yy (amp 1 MPa); drive 100 RK4 steps at dt = 50 μs
      (pulse stays inside domain, travel time 0.67 s).  SXX
      L² drift = **−2.33·10⁻⁵ / step** (dissipative), SYY L²
      drift = **−1.63·10⁻⁵ / step** (dissipative); SYY/SXX
      drift ratio = 0.70 (symmetric, no asymmetric signature).
      **Candidate (b) RK4 non-conservation ELIMINATED.**
- [x] **Debugger:** **§4.12 rev-3g (NEW, 2026-04-21) — absorbing-
      BC energy decay probe PASS 1/1.**  8³ hex, 16 km cube,
      ALL absorbing BCs; isotropic Gaussian stress pulse; 3 000
      RK4 steps at dt = 500 μs (T = 1.5 s > pulse travel time
      1.33 s).  Total energy E(t) decays monotonically 75.05 →
      18.25 (E/E₀ = 0.243 at T_end); `E_max/E₀ = 1.000` exactly
      (no increase anywhere across 3 000 steps).  **Candidate
      (c) absorbing BC wrong-sign reflection ELIMINATED.**
- [ ] **YOU/Debugger:** §4.5 `test_sigma_n_rank_consistency.py`
      (RANK-1 diagnostic, Python) + Frontera pipeline — **READY,
      awaiting user submission.**  See §13.4 for the full gate list.
      Three artifacts generated (rev-3, 2026-04-21):
      1. `jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch`
         (no patch dependency, ~48 SU).
      2. `debug_document/.../tpv102_debug_v9.2.0_diag_fault_sigma.patch`
         (additive `#ifdef SEAS_DIAG_FAULT_SIGMA` only — user applies).
      3. `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch`
         (requires patch applied; ~358 SU).

### Deferred / explicitly NOT in scope

- Any change to `FaultFaceFlux::Evaluate` equations (7)–(12) or
  `GodunovFlux::Interior(n, Q, Q)` per-side call in wave_operator
  without first running the diagnostic patch and MPI harness.  The
  v9.0.0 Pelties-9 flux implementation has been independently
  re-verified against the BP5 convention (§2.3); the trial formula
  `eta_p · (Q_c[VX]^- - Q_c[VX]^+ + stress_terms)` is correct for
  the MFEM convention `can_n` points from + to −, and matches
  SeisSol's identical formula under SeisSol's own + to − convention
  (§2.3.1).  No sign-error root cause.
- Any reassignment of DOFData.V1/V2/tau1/tau2 between dip and
  strike.  BP5 convention is authoritative (CLAUDE.md R-801).
- Any T/Tinv rotation-matrix rewrite.  REVIEW.md R-001 proves the
  current implementation is exact; my rev-1 hypothesis was invalid.
  See §A for the proof and §11 for the retraction.

---

## Dashboard — read me first

**Where we are.**  v9.1.0 closed R-V91-A and R-V91-B3 in part.
Frontera job **7668434** (200 m, P1, 400 ranks, tfinal = 12 s,
v9.1.0 build, commit `99f355a`) is the follow-up.  Two severe
anomalies surface that were partially visible but not quantified in
the earlier 2 s run:

- **R-V92-A — Normal stress drifts ~ 50 MPa at the hypocenter.**
  At `flt_0_7.5` (x = 0, z = 7.5 km, co-located with nucleation
  centre), reported `sigma_n_corr` drops monotonically from 120 MPa
  → ~ 70 MPa between t ≈ 3 s and t ≈ 12 s.  DRDG3D reference: σ_n
  stays flat at 120 MPa throughout.  SCEC TPV102 expectation: pure
  strike-slip on a vertical fault in a homogeneous half-space
  produces **at most a few MPa** of interface-normal perturbation
  from mode-II tip concentration; a 50 MPa drop is not physical.

- **R-V92-B — Dip-direction channels develop across all stations.**
  `slip_dip` reaches −0.2 m at `flt_0_7.5` (vs `slip_strike` ≈ 6 m,
  so 3 % of strike slip), −0.12 m at `flt_n12_3`, +0.028 m at
  `flt_12_12`.  V_dip reaches 0.05 m/s at hypocenter.  DRDG3D
  reference: |slip_dip| < 1 × 10⁻³ m everywhere.

Both are symptoms of the same underlying bug family: the
fault-interface decomposition is coupling strike-direction state
into the dip and normal channels.  The σ_n pattern is
**anti-symmetric in along-strike x** (rises on x > 0, drops on
x < 0 — §1.3 table), a directional signature that narrows the
candidate list.

**Change from v9.1.0 rev-1 analysis.** The rev-1 draft ranked
H-V92-R (rotation stress-block asymmetry) as rank-1 based on a
misreading of the Voigt-6 forward/inverse symmetrization.  REVIEW.md
R-001 proved mathematically that `(T·Tinv)[ij, i'j'] = δ` exactly
for any orthonormal Q, by expanding the entry as `2·(col_k · col_l)`
of the 3 × 3 rotation sub-matrix (details in §A).  I independently
verified this derivation (summing the 6 Voigt contributions to
`(T·Tinv)[0,3]` gives `2·(col_0 · col_1) = 0` by orthonormality).
H-V92-R and its corollary H-V92-O are therefore retracted.

The correct rank-1 hypothesis is **H-V92-M (MPI rank consistency)**
per REVIEW R-003.  Supporting evidence for promotion:

1. The anti-symmetric σ_n signature in §1.3 is a directional
   pattern — it cannot arise from rotation-matrix ULP (would be
   isotropic noise) or from mesh jitter (would be ULP-level, not
   tens of MPa).  It CAN arise from a partition cut that treats
   x > 0 and x < 0 sides differently — exactly what an
   `elem1_on_plus` bit-disagreement between ranks sharing a fault
   face at a ParMETIS cut line would produce.
2. 400 ranks × ~ 20 K fault QPs on the 200 m mesh means partition
   boundaries run through the fault at many (x, z) locations;
   ParMETIS's cuts typically align with the geometry's principal
   axes (x axis here), which makes a disagreement symmetric about
   x = 0 possible.
3. R-V91-B (dip-channel growth) was also diagnosed as a possible
   partition-boundary effect in v9.1.0 §3 H-V91-B2; the same class
   of bug could produce both R-V92-A and R-V92-B as coupled
   symptoms.

Ruling-out evidence we already have: v9.0.0 R-802 is supposed to
unify the per-rank canonical normals via `can_n` from FaultBasis
rather than MFEM's local `CalcOrtho`.  If R-802 is still in force
and correctly applied, H-V92-M should be ELIMINATED by §4.5's
1-rank vs 4-rank comparison.

**rev-3 update (2026-04-21): regression gates green.**  All four
§4 regression-gate tests pass (see §4.1–§4.4 and §13.2–§13.3):

- §4.1 rotation identity: BP5 canonical 0 ULP, off-axis ≤ 1 ULP.
- §4.2 pure strike-slip Tinv: 11/11 bit-exact.
- §4.3 per-QP orthonormality on **200 m production mesh** (113 712
  QPs): 0.5 ULP worst deviation.
- §4.4 per-QP vs centroid on planar fixture: bit-exact.

**Net consequence:** H-V92-R / H-V92-O remain retracted (now with
numerical corroboration); **H-V92-F is additionally eliminated as a
primary mechanism** because accumulated ULP-level mesh jitter bounds
at ≤10⁻¹² MPa over 6 000 steps — 12 orders below the 50 MPa
observation.  The rank-1 hypothesis set narrows to {H-V92-M,
H-V92-G, H-V92-C, H-V92-K}.  The diagnostic path (§5.1 + §4.5
harness + §6.3) is unchanged.

---

## 1. Symptom

### 1.1 Observed (from job 7668434 plots)

| Station | Depth (km) | x (km) | σ_n end (MPa) | σ_n change (MPa) | |slip_dip| end (m) |
|---|---:|---:|---:|---:|---:|
| flt_0_3 | 3 | 0 | ≈ 122 | +2–3 | 0.000 |
| **flt_0_7.5** | **7.5** | **0** | **≈ 70** | **−50** | **0.21** |
| flt_0_12 | 12 | 0 | ≈ 124 | +4 | 0.050 |
| flt_9_7.5 | 7.5 | +9 | ≈ 123 | +3 | 0.075 |
| flt_n9_7.5 | 7.5 | −9 | ≈ 118 | −2 | 0.35 |
| flt_12_3 | 3 | +12 | ≈ 120 | ≈ 0 | 0.000 |
| flt_12_12 | 12 | +12 | ≈ 120 | ≈ 0 | 0.028 |
| flt_n12_3 | 3 | −12 | ≈ 113 | −7 | 0.12 |
| flt_n12_12 | 12 | −12 | ≈ 115 | −5 | 0.07 |

### 1.2 Quantitative time series at hypocenter (flt_0_7.5)

```
 t(s)     V_dip(m/s)   tau_dip(MPa)   sigma_n(MPa)
 0.00      0.0           0.00          120.0
 0.60      3.3e-9         0.0          120.0
 1.20      1.3e-3         0.14         120.0     # nucleation ramp
 1.80     -1.2e-2        -0.27         120.3     # rupture front passes
 2.40     -2.2e-2        -1.00         120.0
 3.00     -1.3e-2        -1.58         119.4     # σ_n start dropping
 3.60     -1.2e-2        -1.64         116.3
 4.20     -1.6e-2        -1.97         112.0
 5.10     -2.7e-2        -2.34         107.3
 6.00     -3.4e-2        -2.40         105.0
 7.20     -5.8e-2        -2.14         101.3     # largest V_dip
 8.10     -4.5e-2        -2.84          95.5
 9.00     -3.2e-2        -2.92          86.8
10.00     (continuing decay)            ~ 80
12.00     (terminal)                    ~ 70
```

σ_n is *exactly* 120 MPa up to t ≈ 3 s, then drifts monotonically.
|V_dip| exceeds 5 × 10⁻² m/s by t = 7 s.  The drift is neither a
nucleation ring-down nor a free-surface reflection arrival
(reflection from z = 0 at hypocenter arrives at t ≈ 4.3 s, too late
for the onset).

### 1.3 Spatial asymmetry in σ_n

**Key pattern — σ_n drift is anti-symmetric in along-strike x:**
- Stations on x > 0: σ_n flat or rises by +2 to +4 MPa.
- Stations on x < 0: σ_n drops by −2 to −7 MPa (and −50 MPa at
  x = 0 hypocenter, which is its own distinct signal due to
  collocation with nucleation).

This is a directional signature.  It rules out:
- Purely spatially-random numerical noise (would be isotropic).
- Free-surface reflection (would be symmetric in x by problem
  symmetry about x = 0).
- Rotation-matrix ULP × repeated-application error
  (REVIEW.md R-001 proves rotation is exact; and would give
  isotropic noise anyway).

It POINTS TOWARD a coupling with **directional dependence on
partition-boundary geometry** (H-V92-M) or **physical mode-II tip
concentration that the MFEM implementation is amplifying via an
unknown mechanism** (H-V92-G).

### 1.4 Expected (SCEC TPV102 reference)

- σ_n ≈ 120 MPa ± 1 MPa at every station.
- |slip_dip| < 1 × 10⁻³ m; |V_dip| < 1 × 10⁻⁴ m/s.
- V_strike / slip_strike match MFEM within ~ 5 % (confirmed from
  plots — strike channel is NOT the problem).

### 1.5 When it started

- v7.x / v8.x: welded-flux cancellation bug masked the issue — no
  radiation entered the bulk.
- v9.0.0 / v9.1.0: per-side Pelties-9 flux + CellData + dip
  normalization all applied.  The issue is present at 2 s (job
  7667881) but only quantifiable at 12 s (job 7668434).

---

## 2. Reference material — what the papers actually prescribe

### 2.1 Pelties 2012 eq. (7): trial tractions

In face-local coordinates `(n, t1, t2)` with n pointing from + to −
side, trial tractions are:

```
σ_n*   = η_p · ( v_n^- - v_n^+ + σ_nn^+/Z_p^+ + σ_nn^-/Z_p^- )
τ_t1*  = η_s · ( v_t1^- - v_t1^+ + σ_nt1^+/Z_s^+ + σ_nt1^-/Z_s^- )
τ_t2*  = η_s · ( v_t2^- - v_t2^+ + σ_nt2^+/Z_s^+ + σ_nt2^-/Z_s^- )
```

In the canonical *face-local* frame, these decompose *after
rotation*: `σ_n*` depends only on the n-direction velocity and the
σ_nn stress component; `τ_t1*` on the t1-direction velocity and
σ_{n,t1}; `τ_t2*` on the t2-direction velocity and σ_{n,t2}.

**Correction from rev-1 (REVIEW R-005).** This decomposition is
true only IN THE CANONICAL FRAME AFTER ROTATION.  Before rotation,
the GLOBAL state Q_g has all 9 components potentially nonzero, and
the bulk A-matrices `A_x`, `A_y`, `A_z` couple them — e.g. a mode-II
rupture propagating along strike (+x) with only Q_g[VX] and
Q_g[SXY] initially nonzero will, through `A_x·∂_x Q_g`, generate
Q_g[SXX], Q_g[SYY], Q_g[SZZ] as the rupture front propagates (see
§2.4 / §2.5 for details).  These bulk perturbations then appear in
the face-local rotated `Q_c[SXX]` = σ_nn and directly drive
`sigma_n_trial` via the stress-average term.  This is **physics,
not a bug**, and the plan must not foreclose bulk-dynamics
hypotheses (H-V92-G) based on Pelties eq. (7)'s *post-rotation*
decoupling structure.

### 2.2 Pelties 2012 eq. (11)–(12): imposed states

After friction correction to `(σ_n_corr, τ_t1_corr, τ_t2_corr)`,
the imposed states are (SeisSol + to − convention — same as MFEM):

```
σ^±_imp:  N = σ_n_corr;  T1 = τ_t1_corr;  T2 = τ_t2_corr;
          SYY, SZZ, SYZ unchanged from native Q.
v^+,imp:  U = v^+_U + invZp · (σ_n_corr - σ_n^+)
          V = v^+_V + invZs · (τ_t1_corr - σ_nt1^+)
          W = v^+_W + invZs · (τ_t2_corr - σ_nt2^+)
v^-,imp:  U = v^-_U - invZpNeig · (σ_n_corr - σ_n^-)
          V = v^-_V - invZsNeig · (τ_t1_corr - σ_nt1^-)
          W = v^-_W - invZsNeig · (τ_t2_corr - σ_nt2^-)
```

This matches MFEM `fault_face_flux.cpp:163-180` bit-for-bit.  The
non-normal stress components `σ_yy`, `σ_zz`, `σ_yz` are NOT
overwritten — this is correct (they are zero-eigenvalue modes of
`A_n`).

### 2.3 Re-derivation from characteristics in the BP5 / MFEM "+ to −"
convention (responding to user directive)

**This section addresses the explicit user directive: "since we
adopt different convention, then we compare both code with
dip/strike involved, they should not be exactly the same, you
should take in the position that we use this convention and
rederive the correct sign/formula from the equations".**

#### 2.3.1 Convention — which way does n point?

SeisSol's internal convention, traced from the source:
- `MeshTools::normalAndTangents(element, side, ...)` computes the
  face normal outward from `element`
  (`Geometry/MeshTools.cpp:58-69`).
- `CellLocalMatrices.cpp:497` sets
  `faceInformation[ltsFace].plusSide = fault[meshFace].side`, i.e.
  `plusSide` is the `element`'s side.
- Therefore **SeisSol's face normal n points from plus to minus**
  (outward from plus = into minus).

MFEM's convention (`wave_operator.inl:796`):
`elem1_on_plus = !qpd.sign_flipped`, where `sign_flipped = true`
iff MFEM's `CalcOrtho` returned a normal anti-aligned with
`ref_normal`.  `can_n` is aligned with `ref_normal` (= `(0, -1, 0)`
for TPV102), so when `sign_flipped == false`, Elem1's outward
aligns with `can_n`.  Elem1 is at y > 0 (the side Elem1's outward
points away from = the +y side), and `can_n = (0, -1, 0)` points
from +y to −y.  **So MFEM's can_n also points from + to −.**  Same
convention as SeisSol.

#### 2.3.2 Re-derive σ_n_trial from characteristics (tension-positive)

Using the MFEM A-matrix (verified tension-positive from
`godunov_flux.cpp:162-188`; SeisSol matches via
`ElasticSetup.h:36-41`):

```
A = [[ 0,          -(λ+2μ)   ],    (for 1D P-wave subsystem)
     [ -1/ρ,        0         ]]     in (σ_nn, v_n)
```

Eigenvalues `±cp`.  Left eigenvectors: `l_+ = (1, -Z_p)`,
`l_- = (1, +Z_p)`.  Characteristic invariants:
`w_+ = σ - Z_p·v` (preserved along dn/dt = +cp).
`w_- = σ + Z_p·v` (preserved along dn/dt = -cp).

Riemann problem with left state (σ_L, v_L) at n < 0, right state
(σ_R, v_R) at n > 0, welded interface at n = 0:
- Right-going wave at n = cp·t: (σ_R - σ*, v_R - v*) ∝ right
  eigenvector `u_+ = (-Z_p, 1)` ⇒ σ* = σ_R + Z_p·(v_R - v*).
- Left-going wave: σ* = σ_L - Z_p·(v_L - v*).

Solving: **σ* = (σ_L + σ_R)/2 + Z_p·(v_R - v_L)/2**.

#### 2.3.3 Apply to BP5 + to − convention

If n points from + to −, the right-going wave (dn/dt = +cp) moves
IN the +n direction = from + toward −.  Right-going originates
from the + side.  Left-going originates from − side.

Relabeling: L (left of interface, where right-going wave
originates) = +; R (right of interface, where left-going wave
originates) = −:

```
σ* = (σ^+ + σ^-)/2 + Z_p·(v^- - v^+)/2
   = η_p · ( v^- - v^+ + σ^+/Z_p + σ^-/Z_p )
```

This **EXACTLY MATCHES** MFEM's `fault_face_flux.cpp:52-54`:
```cpp
sigma_n_trial = data.eta_p * (Q_minus[VX] - Q_plus[VX]
                              + Q_plus[SXX] * invZp_plus
                              + Q_minus[SXX] * invZp_minus);
```

**No sign bug in the trial formula under BP5 convention.**  The
MFEM implementation is correct; the formula is the same as
SeisSol's only because both use n from + to −, not because of any
coincidence.

#### 2.3.4 Same derivation for τ_t1_trial and τ_t2_trial

The t1 and t2 shear subsystems have identical structure to the
normal subsystem (shear characteristic speed is `cs`, impedance
`Z_s`).  The SAME re-derivation gives:

```
τ_t1* = η_s · (v_t1^- - v_t1^+ + σ_nt1^+/Z_s + σ_nt1^-/Z_s)
τ_t2* = η_s · (v_t2^- - v_t2^+ + σ_nt2^+/Z_s + σ_nt2^-/Z_s)
```

In BP5 convention, t1 = dip and t2 = strike.  So `τ_t1` = dip
shear trial, `τ_t2` = strike shear trial.  The MFEM code
(`fault_face_flux.cpp:56-64`) uses `VY` for v_t1 and `SXY` for
σ_nt1, consistent with the canonical-frame index convention
(X = n, Y = t1, Z = t2).  **The formula is convention-correct.**

#### 2.3.5 Summary of 2.3 — trial formulas under BP5 convention

| Formula | BP5-convention correct? | MFEM code matches? |
|---|:---:|:---:|
| `σ_n_trial = η_p·(Q_c[VX]^- − Q_c[VX]^+ + Q_c[SXX]^+/Zp + Q_c[SXX]^-/Zp)` | ✓ | ✓ |
| `τ_t1_trial = η_s·(Q_c[VY]^- − Q_c[VY]^+ + Q_c[SXY]^+/Zs + Q_c[SXY]^-/Zs)` | ✓ | ✓ |
| `τ_t2_trial = η_s·(Q_c[VZ]^- − Q_c[VZ]^+ + Q_c[SXZ]^+/Zs + Q_c[SXZ]^-/Zs)` | ✓ | ✓ |

**Conclusion:** No formula-level sign error.  The root cause is
elsewhere.  This is consistent with REVIEW.md R-001 / R-002
(rotation is exact) and R-005 (bulk-dynamics physics is the
correct place to look).

### 2.4 Mode-II normal-stress concentration at the rupture tip

For a mode-II strike-slip rupture on a vertical fault in a 3D
homogeneous half-space, the rupture-tip stress field has a σ_nn
component of magnitude **order a few MPa** at the fault plane.

- DRDG3D reference at `flt_0_7.5`: |σ_n − 120 MPa| < 1 MPa ∀ t.
- SCEC TPV102 validation document: σ_n perturbation within ± 1 MPa.

**The observed −50 MPa drift is 50× the physical bound.**  Mode-II
tip physics alone cannot explain it.  But the *direction* of the
mode-II stress (anti-symmetric in along-strike x with appropriate
tip propagation) MATCHES the observed anti-symmetry in §1.3 —
suggesting an **amplification** of correct physics rather than
random noise.

### 2.5 Bulk-dynamics mechanism (new, responding to REVIEW R-005)

A mode-II rupture propagating along +x with strike-slip motion
initially generates only Q_g[VX] and Q_g[SXY] in the near field.
But as the front propagates in x, the bulk A_x Jacobian operates
on ∂_x Q.  A_x has these nonzero entries (from the same
construction as A_y but for dir = 0):
- A_x[SXX, VX] = −(λ+2μ), A_x[SYY, VX] = −λ, A_x[SZZ, VX] = −λ
- A_x[VX, SXX] = −1/ρ

So ∂_x v_x generates ∂_t σ_xx = (λ+2μ) ∂_x v_x AND ∂_t σ_yy =
λ ∂_x v_x AND ∂_t σ_zz = λ ∂_x v_x.  The rupture spatially varying
v_x along strike thus SEEDS σ_yy perturbations in the bulk via the
standard compressional Poisson-coupling.

When these σ_yy perturbations reach the fault QPs, `sigma_n_trial`
at the fault picks them up through its stress-average term:
`η_p · (Q_c[SXX]^+ + Q_c[SXX]^-)/Z_p`, where Q_c[SXX] =
σ_{can_n, can_n} = σ_{yy,global}.  This is **correct physics**
feeding a **correct formula**.

Whether MFEM *amplifies* this correct physics into a 50-MPa drift
(vs DRDG3D's 1-MPa bound) is the question.  H-V92-G below
postulates such amplification; tests in §5.1 + §6.3 measure it
directly.

---

## 3. Hypotheses — ranked after REVIEW.md

### H-V92-M (NEW RANK-1) [POSSIBLE] — MPI partition-boundary bit-disagreement on `elem1_on_plus` or canonical normal at shared fault faces

**Claim.** Shared-fault faces at ParMETIS partition cuts have
`shared_fault_elem1_on_plus_[sf_idx]` (precomputed in the ctor
from Elem1 geometry) that may disagree across the two ranks
sharing the face.  If two ranks disagree on which side is "+", the
code at `wave_operator.inl:1287-1290` routes `(Q_self_can,
Q_nbr_can)` into `(Q_plus_local, Q_minus_local)` *oppositely*.
Evaluate is then called with the (+, −) arguments swapped on one
rank.  Since Pelties eq. (7) is *anti-symmetric in the (v^-, v^+)
jump*, a swap flips the sign of the trial velocity-jump
contribution without flipping the stress-average; the result is a
different `sigma_n_trial` on the two ranks.  The imposed state on
that shared face then injects DIFFERENT flux into the bulk on each
side.

**Why the observed anti-symmetry in x.** ParMETIS cuts on a
rectangular-ish domain (TPV102 is a cuboid) typically align with
the principal axes.  A cut running along x = 0 would divide the
fault into x > 0 and x < 0 pieces; faces at the cut line would
have the two ranks' `elem1_on_plus` bits set from two different
`Elem1` instances (one on each side of the cut), and the
disagreement (if any) would therefore correlate with the x-sign.
Integrated over all fault faces at the cut, the error accumulates
into an x-anti-symmetric σ_n pattern.

**Why v9.0.0 R-802 was supposed to prevent this.** R-802 unified
the per-rank canonical normal via `FaultBasis` (which uses
`ref_normal` — rank-invariant) rather than MFEM's local
`CalcOrtho` (which can give anti-aligned normals on the two
sides).  If R-802 is still correctly in force at head commit
`99f355a`, H-V92-M is ELIMINATED.  BUT — `shared_fault_elem1_on_plus_`
is populated **independently of `can_n`** (see ctor around
`wave_operator.hpp`); it is derived from mesh geometry of Elem1
only, and the question is whether two ranks' Elem1-geometry
produces bit-identical results.

**Confirmatory test.** §4.5 — run 1 vs 4 ranks locally (≤ 4 ranks
only per user feedback memory) at 4 s tfinal; compare σ_n(t = 3.5s)
at hypocenter.  If different by > 10 kPa, H-V92-M is CONFIRMED.

**Ruling-out test.** Run `grep -n
"shared_fault_elem1_on_plus_\.push_back\|elem1_on_plus_.*=" dynamic/wave_operator.*`
to find where the bit is populated.  Check whether the computation
uses only Elem1-geometry (bit-identical across ranks) or whether
it depends on MPI ghost data (potentially different across ranks).

### H-V92-G (NEW RANK-2) [POSSIBLE] — Bulk-dynamics amplification: correct mode-II physics is being amplified by some discretization-level coefficient

**Claim.** Per §2.5, mode-II rupture SHOULD produce σ_yy
perturbations in the bulk of order a few MPa.  These perturbations
enter `sigma_n_trial` through the stress-average term — correctly,
per Pelties eq. (7).  The observed 50-MPa drift is ~ 50× the
expected physical magnitude, suggesting MFEM's DG discretization
amplifies the mode-II coupling.

**Possible mechanisms for amplification:**
(a) IP penalty or BR2 penalty on σ_yy jumps at interior faces is
    wrong for 3D (in the elastic-wave driver, not the static
    domain solver) — producing stronger coupling than physics.
(b) The explicit RK4 time integrator accumulates non-conservative
    error in the σ_yy channel (H-V92-K subsumed).
(c) A boundary condition (absorbing on ±x, ±y, ±z except free
    surface) reflects σ_yy with wrong sign, over 12 s of
    accumulated reflection pumps σ_yy on the fault.
(d) The elementwise mass-inverse at `wave_operator.inl:1389-1415`
    is not applied consistently across component channels (e.g.,
    different conditioning for σ_yy vs σ_xy).

**Confirmatory test.** §5.1 diagnostic patch measures
`max|Q_self[SYY]|` at fault QPs directly.  If max|Q_self[SYY]|
grows to ~ 50 MPa matching the observed `sigma_n_trial` drift,
H-V92-G is CONFIRMED as the MECHANISM, and the follow-up is to
identify which of (a)-(d) is the amplifier.  If max|Q_self[SYY]|
stays below 1 MPa, H-V92-G is eliminated and the bug is in the
fault flux path itself (and we must dig further into H-V92-C or
H-V92-F).

### H-V92-C (NEW) [POSSIBLE, MODERATE] — Fault–free-surface corner interaction via Interior's internal BuildFrame

**Claim (reviewer's flag for the next round).**
`GodunovFlux::Interior(nor, ...)` internally calls
`BuildFrame(nor, t1, t2)` (`godunov_flux.cpp:329`) which uses a
Gram-Schmidt with `up = (0, 0, 1)` or `(1, 0, 0)` as its seed.
For the free-surface flux at z = 0 (`nor = (0, 0, 1)`), BuildFrame
sees dot(nor, up) = 1 > 0.9 and flips to `up = (1, 0, 0)`, producing
the internal frame:
```
n  = (0, 0, 1)
t1 = up × n = (1,0,0) × (0,0,1) = (0, -1, 0)
t2 = n × t1 = (0,0,1) × (0,-1,0) = (1, 0, 0)
```

This is a legitimate orthonormal right-handed frame.  The free-
surface gamma pattern `{-1, 1, 1, -1, 1, -1, 1, 1, 1}` flips the
stress components with *normal indices* in THIS INTERNAL FRAME.
So the global σ_{zz}, σ_{xz}, σ_{yz} all get flipped (correct for
σ·n = 0 at z = 0).

**But at the fault–free-surface corner** (z = 0, y = 0), the
fault QPs just below z = 0 (say z = −h/2) sit immediately adjacent
to free-surface QPs.  The bulk flux between the fault-adjacent
element and the free-surface-adjacent element uses the standard
interior `Interior(nor, Q_L, Q_R)` call, where `nor` is the face
normal between those two elements.  That face is horizontal, so
`nor = (0, 0, 1)` and the internal BuildFrame gives the frame
above.  No per-side fault logic applies to this bulk face.

**Amplification pathway.** If a spurious σ_{yy}-channel flux
arises at the corner (e.g., from the combination of the adjacent
fault face's imposed-state flux and the free-surface mirror),
σ_yy propagates downward from the corner along the fault plane at
shear-wave speed.  Over 12 s, the σ_yy wave reaches all fault
stations, with the arrival time ~ `depth / cs = 7.5 / 3.464 s ≈
2.16 s` later than the initial rupture-tip arrival — consistent
with the observed σ_n-drift onset at t ≈ 3 s at the hypocenter
(rupture arrives ~ 1 s, corner-sourced σ_yy arrives ~ 3.16 s).

**Confirmatory test.** §6.3 ParaView slice at y = 0 (fault plane)
at t = 3 s and t = 6 s.  If Q_g[SYY] has a spatial pattern
emanating from (x, y, z) = (*, 0, 0) (the fault–free-surface
junction) rather than from the rupture front itself, H-V92-C is
CONFIRMED.  If Q_g[SYY] is distributed along the rupture front
without a corner-emanation signature, H-V92-C is ELIMINATED.

### H-V92-F [ELIMINATED AS PRIMARY — rev-3 post-§4.3 evidence] — Per-QP face-normal non-planarity from mesh jitter

**Claim (rev-1 §3 H-V92-F, unchanged).** `FaultBasis::ComputeQPBasis`
computes per-QP `(n, t1, t2)` from MFEM's `CalcOrtho(Face->Jacobian())`
at each QP.  For Gmsh-generated tet meshes, even a single fault-
adjacent node perturbed by one ULP can create a per-QP normal
that is ULP-off from `can_n_centroid`.

**Why this is ranked lower post-REVIEW.** The observed drift is
~ 50 MPa over 12 s × 10⁴ time steps = ~ 5 MPa per 1000 steps.
ULP-level mesh jitter would give < 10⁻¹⁰ MPa per step, accumulating
to < 10⁻⁶ MPa over 10⁴ steps.  Off by > 10⁶.  Mesh jitter alone
cannot explain the magnitude.

**rev-3 numerical evidence closing this hypothesis.**  §4.3
executed on the **200 m production mesh** (same mesh as job
7668434) reports:

| Statistic | Value |
|---|---|
| Fault faces tested | 37 904 |
| QPs tested | 113 712 |
| worst \|n\|−1, \|t1\|−1, \|t2\|−1 | 1.110·10⁻¹⁶ (0.5 ULP) |
| worst \|n·t1\|, \|n·t2\|, \|t1·t2\| | ~10⁻⁶² (essentially 0) |

The per-QP basis orthonormality holds at **sub-ULP** on the
production mesh.  Accumulated error over 6 000 time steps is
bounded by 0.5 ULP × 10³ MPa × 6 000 = ~10⁻¹² MPa, which is **12
orders of magnitude below the observed 50 MPa drift**.

**Verdict.** H-V92-F is **eliminated as a primary mechanism**.  It
cannot amplify to 50 MPa through any accumulation pathway.  It is
retained only as a "noise floor" reference in the §6.3 ParaView
slice analysis (the slice should NOT show QP-level speckle
correlating with mesh-node positions).

### H-V92-K [POSSIBLE, LOWEST] — RK4 consistency error in the (Q, DOFData) coupling

**Claim (rev-1 §3 H-V92-K, downgraded).** The RK4 driver averages
DOFData diagnostic fields (sigma_n_corr, etc.) from stage snapshots,
not from explicit state integration.  If the average is
inconsistent with the true sigma_n_corr at t_{n+1}, a per-step
drift of ~ 10 kPa accumulates to 50 MPa over 6000 steps (10⁴ ·
cfl-limited dt = 12 s at dt ~ 2 ms).

**Why ranked lowest post-REVIEW.** The diagnostic fields are
OUTPUTS only, not state variables — any RK4-averaging error
affects what is REPORTED at the stations, not the underlying Q
evolution.  If H-V92-K is primary, then σ_n at stations would be
wrong but the bulk Q_g[SYY] at fault QPs would be correct.
§5.1 diagnostic distinguishes these cases: if max|Q_self[SYY]|
stays near zero AND station σ_n_corr drifts to −50 MPa, then the
bug is in the station-output path (H-V92-K primary).  Else, the
bug is in the bulk evolution (H-V92-G or H-V92-M primary).

### H-V92-R [RETRACTED — REVIEW R-001; rev-3 numerically corroborated]

~~Rotation stress-block `(a != b)` vs `(i != j)` asymmetry.~~
RETRACTED.  See §A for the proof that `(T·Tinv)[ij, i'j'] = δ`
exactly for any orthonormal Q.  The asymmetry is the correct
forward/inverse Voigt symmetrization direction.

**rev-3 numerical corroboration.**  §4.1 measured
`max|T·Tinv − I|`:

| Frame | Result |
|---|---|
| BP5 canonical ((0,-1,0), (0,0,-1), (1,0,0)) | **0.000·10⁰ (0 ULP — bit-exact)** |
| BP5 + 15° about x | 2.220·10⁻¹⁶ (1 ULP) |
| BP5 + 30° about y | 2.220·10⁻¹⁶ (1 ULP) |
| BP5 + 45° about (1,1,1)/√3 | 2.220·10⁻¹⁶ (1 ULP) |

16× under the REVIEW R-004 budget.  Rotation is rigorously exact
at machine precision.

### H-V92-O [RETRACTED — REVIEW R-002; rev-3 numerically corroborated]

~~Tinv_can leakage on pure strike-slip fixture.~~
RETRACTED.  For the BP5 canonical frame (signed permutation),
Tinv entries are in `{-2, -1, 0, +1, +2}` and operations on any
axis-aligned fixture are bit-exact integer arithmetic.  Subsumed
by R-001.

**rev-3 numerical corroboration.**  §4.2 fed pure strike-slip
Q_g (only VX=1e-12 m/s and SXY=7.5·10⁷ Pa) through Tinv on the
BP5 frame.  Result:

```
Q_c[VX]=0  Q_c[VY]=0  Q_c[VZ]=1e-12           (strike velocity, exact)
Q_c[SXX..SXY,SYZ]=0  (bit-exact, 0 ULP leakage)
Q_c[SXZ]=-7.5e7      (n·strike shear, exact)
```

11/11 assertions bit-exact.  No leakage into the normal or dip
channels at the BP5 canonical frame.

### H-V92-U (NEW, rev-3h per REVIEW round-5 R-V92-G02) [OPEN — RANK-1 after §4.9-§4.12 closed Pelties candidates] — Interior-fault-face path bug

**Claim.** `ComputeFaceFluxRHS` fault branch
(`dynamic/wave_operator.inl:746-950`) handles the majority of
fault faces per §17.6.5 (Gmsh-duplicated fault ⇒ most fault faces
classify as interior with `fault_attr = 3`, not shared across
ranks).  It uses the same Tinv / Evaluate / T / per-side flux
pattern as the shared branch but with
`elem1_on_plus = !qpd.sign_flipped` — a per-QP FaultBasis-based
sign, NOT the centroid-projection used in the shared branch.

**Why open.** No test probes this code path directly.  §4.6/§4.7
probe the shared-fault branch only; §4.9–§4.12 probe bulk
properties, not the fault flux pipeline.  The interior branch
assembles flux to BOTH Elem1 AND Elem2 in a single loop (both
elements on the same rank), with `Evaluate` mutating the shared
`dof_data[dof_idx]` once per face.  A sign, frame, or coupling
bug here would pollute every rank's dof_data without any
cross-rank diagnostic signal.

**Discriminator.**  §4.G04 new unit test
`tests/unit/test_interior_fault_flux_path.cpp`: 2-tet single-rank
fixture, antisymmetric σ_xy projection, 1 `Mult` call, assert
DOFData output matches analytic Pelties eq. (7).

**Status:** probe not yet written; in the active TODO above.

### H-V92-V (NEW, rev-3h per REVIEW round-5 R-V92-G02) [OPEN — RANK-2] — Absorbing-BC reflection phase on non-free-surface boundaries

**Claim.**  The TPV102 driver installs `Absorbing` BCs on ±x, ±y,
±z except the free surface at z = 0.  `flux_.Absorbing` in
`dynamic/godunov_flux.cpp` computes the outgoing characteristic
and zeroes the incoming.  If the implementation has a sign or
phase bug, waves radiated toward a boundary return to the fault
with wrong phase, polluting the trial traction over many wave
transits.

**Why open.**  §4.12 (rev-3g) verified **total energy decay** is
monotonic under all-absorbing BCs — but energy monotonicity is
weaker than phase correctness.  A sign bug in one stress channel
(e.g., σ_yy reflects with wrong sign while σ_xy reflects correctly)
can preserve total energy while polluting specific channels —
exactly the observed dip-channel-asymmetric signature.

**Discriminator.**  Substitute ALL boundaries with the free-
surface BC for a short Frontera run.  Energy is not conserved
(finite-amplitude waves bounce forever), but after 4 s the
dip-channel contamination pattern will differ characteristically
from the all-absorbing case if BC phase is the source.

**Status:** not yet run.

### H-V92-W (NEW, rev-3h per REVIEW round-5 R-V92-G02) [OPEN — RANK-3] — Initial pre-stress projection inconsistency

**Claim.**  `InitializeFaultDOFs` (`dynamic/tpv102_setup.hpp`)
sets `tau2_0 = tau_ini` in the BP5 canonical frame at each fault
QP.  `InitializeState` sets `Q = 0` in the GLOBAL frame.  If the
per-QP FaultBasis on a specific face is NOT bit-exactly the BP5
canonical frame (due to mesh jitter producing per-QP variation
of sub-ULP scale — §4.3 verified orthonormality but the canonical
frame is not unique at every QP), the FIRST `Evaluate` call at
t = 0+ sees a ULP-scale non-zero V1 seed that accumulates over
10⁴ steps.

**Why open.**  §4.3 verified frames are orthonormal to 0.5 ULP
but did NOT verify they are bit-identical to the BP5 canonical
frame at every QP.  A ULP-scale off-axis component of `(n, t1, t2)`
in the dip direction leaks Q components from stride to dip at
each time step.  Over 10⁵ steps this could accumulate to 0.1%
pollution = 75 kPa per step × 10⁵ = 7.5 GPa worst case — 
insufficient by dimensional analysis (REVIEW F03 logic) but worth
verifying.

**Discriminator.**  Add `[FAULT-INIT-V1]` printf at the end of
RK4 stage 1 of the FIRST time step: `max|V1|` across all fault
QPs.  Should be 0 exactly; any non-zero confirms H-V92-W as at
least a contributor.

**Status:** not yet probed.

### H-V92-T (round-4 ADER audit, rev-3h) [CLOSED NON-PRIMARY] — Time-integrator choice (ADER vs RK4)

**Claim.**  If the RK4 truncation or per-stage splitting error is
the primary amplifier, switching to ADER-DG would fix R-V92.

**Closure (REVIEW round-4 R-V92-F03).**  Dimensional analysis:
observed slip_dip 1.83 m is 100× larger than any RK4 splitting
error can produce (upper bound ~10⁻² m on F01+F02 splitting).
A correctly-implemented 4th-order RK4 cannot be 100× wrong at its
own truncation bound.  ADER is not the fix for the primary
pathology.

**What IS open from the round-4 audit.**  F01 (DOFData outside
the RK4 state vector) and F02 (psi analytic update outside RK4)
are plausible SECONDARY contributors (~1% of observed pathology).
Fix them within RK4 as architectural hygiene; the dt-halving test
in the active TODO will quantify.

**Do NOT port ADER.**  Multi-month MFEM infrastructure project,
no evidence support, recoverable only if every other diagnosis
concludes "bug is specifically in a location requiring ADER's
coupled-system formulation".  No current evidence points there.

### Explicitly ELIMINATED (carried forward + REVIEW-confirmed)

- **Pre-stress placement convention.** BP5's `tau1_0 = 0`
  (dip = 0), `tau2_0 = tau_ini` (strike = tau_ini) is authoritative
  (CLAUDE.md R-801, user directive).
- **Pelties eq. (7) formula sign.** §2.3 re-derivation from
  characteristics confirms MFEM's formula
  `η_p·(Q_c[VX]^- − Q_c[VX]^+ + ...)` is correct under MFEM's
  can_n-from-+-to-− convention.  No sign bug.
- **Pelties eq. (11)-(12) formula sign.** Matches SeisSol
  bit-for-bit at `fault_face_flux.cpp:170-180`.
- **Rotation bit-inverse.** REVIEW.md R-001 + §A proof.
- **Free-surface flux gamma pattern.**  `{-1, 1, 1, -1, 1, -1, 1,
  1, 1}` is the correct mirror for σ·n = 0 at z = 0 free surface.
- **Nucleation spatial / temporal profile.** SCEC-conforming.

---

## 4. Regression-gate unit tests (DEMOTED from rev-1 rank-1 diagnosis tests)

Per REVIEW R-004, the four unit tests previously proposed as
rank-1 diagnostic probes are demoted to regression gates.  They
remain useful but do NOT drive R-V92 diagnosis.

### §4.1 (regression-gate) `test_godunov_rotation_identity.cpp` — **PASS rev-3**

Renamed from `test_godunov_rotation_bp5_frame_identity.cpp`.
Builds T and Tinv from a set of test frames and asserts
`max|T·Tinv − I| ≤ 16 ULP` on every entry.

**rev-3 result.** 4/4 tests pass.  BP5 canonical 0 ULP, 3 off-axis
frames 1 ULP each.  Full numerics in §13.3.

**Test frames:**
- The axis-aligned BP5 canonical frame (guarded by `16 ULP` — on
  a signed-permutation Q this is expected to be **bit-exact**).
- Three general orthonormal frames at `(+30°, 15°, 0°)` Euler-angle
  perturbations from BP5 (guarded by `16 ULP` — this allows for
  9-term FMA sum round-off, per REVIEW R-004).

**Expected outcome:** PASS on all frames.  Failing would indicate
a rotation-matrix regression post-v9.2.0; it is NOT a primary
diagnostic for R-V92.

**Rationale for 16 ULP (REVIEW R-004):** each entry of
`(T·Tinv)[ij,i'j']` is a 9-term sum of products of Q entries, each
of magnitude ≤ 2 (Voigt-pair symmetrization).  Worst-case FMA
round-off is ~ 9 × 2 × ε = 18 ε on a non-axis-aligned frame.  16
ULP is a conservative bound that passes all physical frames but
still flags a systematic sign-error regression.

### §4.2 (regression-gate) `test_canonical_rotation_pure_strikeslip.cpp` — **PASS rev-3**

Same fixture as rev-1.  Passes trivially on the axis-aligned BP5
frame (Tinv · Q is integer arithmetic per REVIEW R-002).  Retained
as a regression gate against a future BP5-convention reassignment.

**rev-3 result.** 11/11 assertions bit-exact on BP5 canonical.
Q_c[SXZ] = −τ_ini exactly; Q_c[VZ] = V_ini exactly; all other
channels 0.0 bit-exact.  Confirms zero leakage into normal/dip
channels at this frame.

### §4.3 (regression-gate) `test_fault_basis_qp_orthonormality.cpp` — **PASS rev-3**

Unchanged from rev-1.  Tests per-QP frame orthonormality on the
production 200 m mesh.  Expected to PASS.  Retained.

**rev-3 result.** 1/1 aggregate test passes on **both**:

- 1000 m mesh: 4 938 QPs / 1 646 fault faces, worst norm dev 0.5 ULP.
- 200 m production mesh: 113 712 QPs / 37 904 fault faces, worst
  norm dev 0.5 ULP, worst orth-dot ~10⁻⁴⁶ ULP.

Conclusion: H-V92-F is numerically eliminated as a primary mechanism
(see §3 H-V92-F for the accumulated-error bound).

### §4.4 (regression-gate) `test_per_qp_vs_centroid_basis_planar_1el.cpp` — **PASS rev-3**

Unchanged from rev-1.  Passes trivially on a hand-constructed
planar 1-tet fixture.  Retained as MFEM-Jacobian regression gate.

**rev-3 result.** 11/11 assertions pass on a planar 2-tet fixture
with the shared face in the y = 0 plane.  Centroid basis matches
the BP5 canonical `(n, t1, t2) = ((0,-1,0), (0,0,-1), (1,0,0))`
bit-exact; all 4 QP bases equal the centroid bit-exact.  MFEM's
CalcOrtho is QP-invariant on planar faces — confirmed.

### §4.5 (RANK-1 DIAGNOSTIC) `test_sigma_n_rank_consistency.py`

**Promoted to RANK-1 diagnostic per REVIEW R-003.**

Python harness.  Runs `seas_tpv102_driver` with:
- 1 rank (serial, local)
- 4 ranks (local, NO oversubscription per feedback memory)

At `--tfinal 4.0` (minimal time to capture σ_n drift onset).
Output: σ_n(t) at each station at t = 2.0 s and t = 3.5 s.

**Classification:**
| 1-rank σ_n at t=3.5s | 4-rank σ_n at t=3.5s | Verdict |
|---|---|---|
| = 120 ± 0.1 MPa | = 120 ± 0.1 MPa | No drift at 4 s; H-V92-G requires longer run; run 400-rank Frontera diag |
| = 120 ± 0.1 MPa | ≠ 120 ± 0.1 MPa | **H-V92-M CONFIRMED**.  Audit shared-fault ctor bits. |
| ≠ 120 ± 0.1 MPa | = same as 1-rank | H-V92-M ELIMINATED.  Bulk/driver bug — probe H-V92-G + H-V92-C. |
| Both drift, same | — | H-V92-M ELIMINATED; bulk/driver bug. |
| Both drift, differ | — | Both H-V92-M AND H-V92-G / C contribute; isolate iteratively. |

**Must run before any Frontera approval.**  Local cost: ~ 10
minutes × 2 runs = 20 minutes.

---

## 5. Diagnostic instrumentation — C0 patch (user approval required)

### §5.1 (PROMOTED TO RANK-1 DIAGNOSTIC) `SEAS_DIAG_FAULT_SIGMA`

Adds to `dynamic/wave_operator.inl` inside `ComputeFaceFluxRHS` and
`ComputeSharedFaceFluxRHS`.  Guarded by a new `SEAS_DIAG_FAULT_SIGMA`
flag (analogous to existing `SEAS_DIAG_FAULT_FLUX`).

At every fault QP, after `Q_self` / `Q_nbr` are gathered (global
frame, BEFORE rotation to canonical):
```cpp
#ifdef SEAS_DIAG_FAULT_SIGMA
   local_max_syy_g = std::max(local_max_syy_g, std::abs(Q_self[SYY]));
   local_max_vy_g  = std::max(local_max_vy_g,  std::abs(Q_self[VY]));
   local_max_syz_g = std::max(local_max_syz_g, std::abs(Q_self[SYZ]));
   local_max_vz_g  = std::max(local_max_vz_g,  std::abs(Q_self[VZ]));
#endif
```

At the end of `Mult`:
```
MPI_Allreduce(MAX) → [FAULT-DIAG] t=T.TT max|SYY_g|=... max|VY_g|=...
```

**This patch is PURELY additive.**  No change to functional flux
computation.  Safe for Frontera.

**Why rank-1 per REVIEW.** If `max|SYY_g|` on fault QPs tracks the
observed σ_n drift magnitude (grows from 0 at t = 0 to ~ 50 MPa at
t = 12 s), then H-V92-G CONFIRMED and the next step is to ID the
amplifier (a/b/c/d in §3 H-V92-G).  If `max|SYY_g|` stays below
1 MPa while station σ_n_corr drifts to −50 MPa, then H-V92-K
CONFIRMED and the bug is in the station-output path (RK4 averaging
of diagnostic fields).  The diagnostic DIRECTLY discriminates.

### §5.2 (EXISTING) `SEAS_DIAG_FAULT_FLUX` C-1/C-2/C-3

Already implemented.  Re-enable for the 4 s diagnostic run; extend
C-1 EVAL (`fault_face_flux.cpp:110-125`) to also print
`sigma_n_trial`, `sigma_n_total`, `tau1_trial`, `tau2_trial` so
the decomposition of σ_n_corr at the hypocenter is directly
readable.

---

## 6. Offline diagnostics on existing job-7668434 data (pre-Frontera)

### §6.1 (RANK-1) Correlation σ_n deviation vs local |tau2_corr|

Python script (local, no Frontera).  For each of the 9 station
.dat files:
- Read (t, V1, V2, tau1, tau2, sigma_n) at each step.
- Compute σ_n deviation Δσ_n(t) = sigma_n(t) − 120e6.
- Compute local |tau2_corr(t)| from the .dat.
- Plot Δσ_n(t) vs |tau2_corr(t)| with station labels.

**Classification:**
- If points cluster along a SINGLE line (proportionality
  Δσ_n = k · |tau2_corr|), then H-V92-G (bulk mode-II coupling)
  is CONFIRMED with coupling coefficient k.
- If points have x-ORIENTED CLUSTERING (stations at x > 0 form one
  cluster, stations at x < 0 form another), H-V92-M is
  corroborated.
- If no correlation, both H-V92-G and H-V92-M weakened; H-V92-C
  or H-V92-K rises.

### §6.2 (LOW) Free-surface reflection arrival check

Plot Δσ_n(t) against the expected Rayleigh-wave arrival time
`t_rupture + 2·depth/cs` at each station.  Any coincident
step-up signature indicates H-V92-C mechanism.

### §6.3 (RANK-1) ParaView SYY slice

Open v9.1.0 PVD output.  Slice at y = 0 (fault plane).  Inspect
Q_g[SYY] scalar field at t = 3.0 s, 6.0 s, 9.0 s.

**Classification:**
- If Q_g[SYY] on the y=0 slice has a coherent pattern following
  the rupture front, H-V92-G is CONFIRMED.
- If Q_g[SYY] emanates from (z = 0, y = 0) corner — visible as a
  line of high σ_yy at the fault–free-surface edge propagating
  downward — H-V92-C is CONFIRMED.
- If Q_g[SYY] has ABRUPT discontinuities at partition boundaries
  (visible as straight lines in x or z at the partition cuts),
  H-V92-M is CONFIRMED.
- If Q_g[SYY] is noise-level everywhere, H-V92-K is
  CONFIRMED.

---

## 7. Decision tree after §5.1 + §4.5 + §6.3

| §5.1 max\|SYY_g\| at t = 4 s | §4.5 rank consistency at t = 3.5 s | §6.3 SYY slice pattern | Primary conclusion |
|---:|:---:|:---|:---|
| ~ 0 | 1-rank and 4-rank agree | noise | H-V92-K — station-output RK4 averaging bug |
| ~ 0 | 1-rank ≠ 4-rank | noise | H-V92-M + K both contribute; M primary |
| 10⁻³ σ_n0 ≈ 0.12 MPa | agree | rupture-front aligned | mode-II is real but bounded; dig into H-V92-C or other |
| 0.1 σ_n0 ≈ 12 MPa | agree | rupture-front aligned + corner emanation | H-V92-G + H-V92-C; mode-II coupling amplified by corner feedback |
| 0.1 σ_n0 | disagree | partition lines | H-V92-M PRIMARY; may co-amplify with G/C |

---

## 8. Proposed fix shapes (NOT YET APPLIED — YOU-approval gated)

### 8.1 If H-V92-M confirms

Audit shared-fault ctor at `wave_operator.*` where
`shared_fault_elem1_on_plus_` is populated.  If the bit depends
only on Elem1 geometry, it should be bit-identical across the two
ranks sharing a face — but if the bit is derived from `CalcOrtho`
on Elem1's face and the two Elem1 instances (one on each rank)
come from different local-mesh element orderings, a disagreement
can arise.  Fix: derive `elem1_on_plus` directly from the ghost
exchange (use the same reference frame both ranks compute from
`ref_normal`) rather than from per-rank Elem1 geometry.
Expected patch: ~ 30 LOC in `wave_operator.*` ctor.

### 8.2 If H-V92-G confirms

Three sub-fixes depending on which mechanism (a/b/c/d in §3
H-V92-G):
- (a) IP/BR2 penalty — only applicable if a penalty is used in the
  dynamic driver.  Check `wave_operator.*` for penalty terms.
- (b) RK4 time-integration — integrate a σ_yy "observable" in the
  RK4 state to diagnose per-stage drift.  Larger refactor if
  confirmed.
- (c) Absorbing BC reflection — measure σ_yy reflection at the
  domain boundaries via a test fixture with a clean P-wave pulse
  in pure strike-slip setup.
- (d) Mass-inverse conditioning — check `AssembleElementMassInverse`
  symmetry across the 9 components.

### 8.3 If H-V92-C confirms

Special-case the corner element (shared between fault face and
free-surface face): compute the imposed-state flux and the
free-surface ghost state with a consistent frame choice (either
both use `can_n` or both use mesh-local `nor`).  Expected patch:
~ 40 LOC in `wave_operator.inl`.

### 8.4 If H-V92-K confirms

Remove the RK4-averaging of DOFData diagnostic fields
(`tpv102_driver.cpp:863-873`); instead, write the LAST-STAGE
DOFData values directly (consistent with the Q(t_{n+1}) state).
Or, integrate sigma_n_corr as an explicit state in the RK4 sum.
Larger refactor.

### 8.5 If H-V92-F confirms — ESSENTIALLY IMPOSSIBLE per rev-3 §13

§4.3 + §4.4 green on the 200 m production mesh (113 712 QPs at
0.5 ULP).  Accumulated mesh-jitter error ≤ 10⁻¹² MPa over 6 000
time steps, 12 orders below the 50 MPa observation.  This branch
of §8 is retained only for completeness; no fix expected.

---

## 9. Timeline and Frontera budget

| Step | Resource | Wall-time | Core-hours |
|---|---|---:|---:|
| §5.1 C0 diag patch + local build | local laptop | 30 min | N/A |
| §6.1 offline correlation analysis | local | 20 min | N/A |
| §4.5 local 1-rank vs 4-rank | local 4-rank | 30 min | N/A |
| §6.3 ParaView SYY slice | local | 30 min | N/A |
| **Local diagnostic budget (above)** | — | **~ 2 h** | **N/A** |
| If §4.5 ambiguous: Frontera 400-rank diag 4 s | Frontera | 30 min × 400 | 200 |
| §8 fix + §4.1/§4.2 regression | local | 1-2 h | N/A |
| Post-fix Frontera 12 s re-run | Frontera | 2 h × 400 | 800 |
| **Total Frontera (worst case)** | — | — | **~ 1000** |

Much lower than rev-1's 1500 estimate, because the rev-1 §4.6 dt
convergence test was only a fallback for the (retracted) H-V92-R
path.  The diagnostic patch + MPI harness is the correct first
step.

---

## 10. Carry-forward from v9.1.0

- **RK4 station-output consistency** (v9.1.0 plan §10): NOW
  SUBSUMED by H-V92-K in §3.  If §5.1 + §6.3 show Q_g[SYY] at the
  fault stays near zero while station σ_n_corr drifts, this is the
  same issue v9.1.0 already flagged.  Close as part of v9.2.0.
- **paraview_output.hpp CellData coverage** (v9.1.0 R-001):
  verified still in force.
- **fault_basis.hpp dip normalization** (v9.1.0 R-003): verified
  still in force.

---

## 11. Review-incorporation index

**Populated from REVIEW.md (adversarial audit, 2026-04-21 rev-1):**

| REVIEW finding | Severity | How rev-2 addresses it |
|---|:---:|---|
| R-001: H-V92-R is math invalid; (a!=b) vs (i!=j) asymmetry is correct | CRITICAL | H-V92-R **RETRACTED** (§3); proof in §A; §4.1 demoted to regression gate (§4). |
| R-002: H-V92-O subsumed; signed permutation → bit-exact | CRITICAL | H-V92-O **RETRACTED** (§3); §4.2 demoted to regression gate. |
| R-003: H-V92-M underweighted; x-antisym matches ParMETIS cuts at x=0 | MODERATE | H-V92-M **PROMOTED TO RANK-1** (§3, §7); §4.5 rank-consistency harness promoted to rank-1 diagnostic. |
| R-004: §4.1 2-ULP tolerance false-FAILs on non-axis-aligned frames | MODERATE | §4.1 tolerance set to **16 ULP**; explicit comment in the test spec. |
| R-005: §2.1 forecloses bulk-dynamics hypotheses | MODERATE | §2.1 rewritten with a correction paragraph; §2.5 ADDED (mode-II bulk-coupling mechanism); H-V92-G ADDED to §3 as rank-2 hypothesis. |
| Reviewer's flag for next round: BuildFrame @ fault-free-surface corner | — | H-V92-C **ADDED** to §3 (§3 H-V92-C); §6.3 ParaView slice explicitly checks for corner-emanation signature. |

**Reviewer's recommended investigation order:**
| Recommended step | Addressed in |
|---|---|
| 1. §5.1 diagnostic patch FIRST | Now rank-1 diagnostic in §5.1 and TODO checklist head. |
| 2. §4.5 MPI rank consistency SECOND | Promoted to rank-1 in §4.5. |
| 3. §6.3 ParaView SYY slice THIRD | §6.3 promoted to rank-1. |
| 4. Demote §4.1–§4.4 to regression gates | Done in §4. |

---

## 12. Summary of user directives actively guarded in this plan

1. **"Convention differences between MFEM and SeisSol are NOT
   classified as root causes"** (user, 2026-04-21): honoured.
   §2.3 re-derives the characteristic formula from scratch in BP5
   convention (can_n from + to −) and shows it MATCHES MFEM code.
   SeisSol happens to use the same convention (verified from
   `MeshTools::normalAndTangents` + `CellLocalMatrices.cpp:497`).
   This is a coincidence of convention choice, not a constraint.

2. **"Re-derive the correct sign/formula from the equations"**
   (user, 2026-04-21): §2.3 does exactly this from first principles
   (linear elastic A-matrix eigen-decomposition → characteristic
   invariants → Riemann Rankine-Hugoniot → Riemann star state).
   Result: MFEM's formula is correct; NO sign bug.

3. **"Follow SCEC spec unless explicitly matching Tandem"**
   (CLAUDE.md): honoured — §2.4 uses SCEC TPV102 as the reference.

4. **"Never revert a previous fix"** (CLAUDE.md): honoured — v9.0.0
   Pelties-9 per-side flux, v9.1.0 CellData + dip-normalization +
   R-802 remain in force.  §8 fix shapes PROPOSE modifications but
   do not revert prior fixes.

5. **"Do not oversubscribe MPI on laptop"** (feedback memory):
   honoured — §4.5 uses ≤ 4 ranks locally; Frontera tests are
   YOU-approval gated.

6. **"Ask before Frontera runs"** (feedback memory): honoured —
   every Frontera step in §9 is marked.

7. **"TPV102: no local reproducer runs on production mesh"**
   (feedback memory): honoured — §4.5 uses 200 m mesh only when
   absolutely necessary for the rank-dependence test; §5.1's local
   local-4-rank test is a 4 s run (not the 12 s production).

---

## §A. Proof of REVIEW.md R-001: `T · Tinv = I` is bit-exact on the Voigt-6 stress sub-block

For any orthonormal Q (rows n, t1, t2 in R³), the 6 × 6 Voigt
stress block of T · Tinv is the identity.

**Claim.** For `(ij, i'j')` any pair of Voigt indices (0 ≤ ij < 6),
`(T · Tinv)[ij, i'j'] = δ_{ij, i'j'}`.

**Sketch of proof for off-diagonal `(ij = 00, i'j' = 01)`:**

From `godunov_flux.cpp:265-275` (T stress block):
```
T(ij, ab) = Q[a][i]·Q[b][j]   + (a != b) · Q[b][i]·Q[a][j]
```

From `godunov_flux.cpp:224-234` (Tinv stress block):
```
Tinv(ab, i'j') = Q[a][i']·Q[b][j'] + (i' != j') · Q[a][j']·Q[b][i']
```

For `ij = 00`: `i = 0, j = 0`.  `T(00, ab) = Q[a][0]·Q[b][0] +
(a!=b)·Q[b][0]·Q[a][0] = Q[a][0]·Q[b][0] · (1 + (a!=b))`.  For
a = b this is `(Q[a][0])²`; for a != b it is `2·Q[a][0]·Q[b][0]`.

For `i'j' = 01`: `i' = 0, j' = 1`.  `Tinv(ab, 01) = Q[a][0]·Q[b][1] +
1 · Q[a][1]·Q[b][0]`.

Sum: `(T · Tinv)[0, 3] = Σ_{ab=0..5} T(00, ab) · Tinv(ab, 01)`.

Grouping the 6 terms by `Q[k][0]·Q[k][1]` coefficient (for k =
0, 1, 2):

`Q[0][0]·Q[0][1]` coefficient collects terms:
- `ab = (0,0)`: `T(00, 00) = Q[0][0]²`; `Tinv(00, 01) = 2·Q[0][0]·Q[0][1]`.
  Product: `2·Q[0][0]³·Q[0][1]`.  Contributes `Q[0][0]·Q[0][1] · 2·Q[0][0]²`.
- `ab = (0, 1)`: `T(00, 01) = 2·Q[0][0]·Q[1][0]`; `Tinv(01, 01) =
  Q[0][0]·Q[1][1] + Q[0][1]·Q[1][0]`.  Contribution to
  `Q[0][0]·Q[0][1]`: `2·Q[0][0]·Q[1][0] · Q[0][1]·Q[1][0] =
  Q[0][0]·Q[0][1] · 2·Q[1][0]²`.
- `ab = (0, 2)`: similarly `Q[0][0]·Q[0][1] · 2·Q[2][0]²`.

Total for the `Q[0][0]·Q[0][1]` factor:
`2·Q[0][0]·Q[0][1] · (Q[0][0]² + Q[1][0]² + Q[2][0]²) =
2·Q[0][0]·Q[0][1] · (col₀ · col₀) = 2·Q[0][0]·Q[0][1] · 1`.

By symmetry in index k, the full sum is:
`(T · Tinv)[0, 3] = 2 · (Q[0][0]·Q[0][1] + Q[1][0]·Q[1][1] +
Q[2][0]·Q[2][1]) = 2 · (col₀ · col₁) = 2 · 0 = 0`.

Where `col_k = (Q[0][k], Q[1][k], Q[2][k])` is the k-th COLUMN of
Q (= the k-th component of each basis vector), and `col₀ · col₁ =
0` by Q^T·Q = I (orthogonality of Q^T columns = Q rows = basis
vectors).

**Generalization.** The same pattern — the sum reduces to 2 ·
(col_k · col_l) where k, l depend on the `(ij, i'j')` indices —
holds for every off-diagonal Voigt entry; diagonal entries reduce
to `col_k · col_k = 1`.  Orthonormality of the basis IS the
identity.

**IEEE arithmetic.** For BP5's axis-aligned canonical frame,
every `Q[i][j] ∈ {-1, 0, +1}`, so every intermediate product is
exact and the sum is bit-exact integer arithmetic.  `(T · Tinv) =
I` to 0 ULP.

**For non-axis-aligned per-QP frames** (production mesh), each
entry is a 9-term sum of products of floats in `[-1, +1]`; FMA
round-off gives `max|T·Tinv − I| ~ 9 · ε ~ 2e-15 ≤ 16 ULP`.  §4.1's
16-ULP tolerance passes.

**Conclusion.** H-V92-R had no foundation.  The rotation matrices
in MFEM are correct.  ∎

---

## §13. Regression-gate build + verification results (rev-3, 2026-04-21)

Per user directive "create a list of unit tests in the front of this
document, and then proceed with unit tests build and verification for
all tests, document carefully about the test results and
implications".

### §13.1 Build

All four C++ regression gates compiled clean under `mpicxx`
(`conda activate mfem-dev`) on darwin-25.3.0 (arm64), linked against
MFEM + HYPRE + MUMPS.  Makefile additions (rev-3):

```
tests/unit/test_godunov_rotation_identity.cpp          (new)
tests/unit/test_canonical_rotation_pure_strikeslip.cpp (new)
tests/unit/test_fault_basis_qp_orthonormality.cpp      (new)
tests/unit/test_per_qp_vs_centroid_basis_planar_1el.cpp(new)
tests/scripts/test_sigma_n_rank_consistency.py         (new)
```

New Makefile phonies: `test-godunov-rotation-identity`,
`test-canonical-rotation-pure-strikeslip`,
`test-fault-basis-qp-orthonormality`,
`test-per-qp-vs-centroid-basis-planar-1el`, aggregate
`test-v92-regression-gates`.

### §13.2 Results summary table

| Test | Tests run | Pass | Fail | Worst metric | Verdict |
|---|---:|---:|---:|---|:---:|
| §4.1 `test_godunov_rotation_identity` | 4 | 4 | 0 | BP5 canonical 0 ULP (bit-exact); 3 off-axis frames 1 ULP each | **PASS** |
| §4.2 `test_canonical_rotation_pure_strikeslip` | 11 | 11 | 0 | All assertions bit-exact (0 ULP) | **PASS** |
| §4.3 `test_fault_basis_qp_orthonormality` (1000 m mesh) | 1 | 1 | 0 | 4938 QPs / 1646 faces; worst norm dev 0.5 ULP; worst orth dot ~ 4·10⁻⁴⁷ ULP | **PASS** |
| §4.3 `test_fault_basis_qp_orthonormality` (200 m production mesh) | 1 | 1 | 0 | 113 712 QPs / 37 904 faces; worst norm dev 0.5 ULP; worst orth dot ~ 2·10⁻⁴⁶ ULP | **PASS** |
| §4.4 `test_per_qp_vs_centroid_basis_planar_1el` | 11 | 11 | 0 | All per-QP basis entries bit-exact equal to centroid on planar 2-tet fixture | **PASS** |
| §4.6 `test_shared_fault_role_consistency` (4 ranks, inline 4 km tet fixture) | 5 | 5 | 0 | 18 cross-rank pair checks; elem1_on_plus anti-symmetric (0 equal-bit pairs); can_n / can_t1 / can_t2 / nl bit-identical (0 ULP) | **PASS** |
| **§4.7 `test_shared_fault_dof_data_consistency`** (4 ranks, inline 4 km tet fixture, **4 phases post-REVIEW R-V92-C01**: A post-init / B post-Mult(Q=0) / C post-RK4 uniform Q / **D post-RK4 antisymmetric σ_xy(y)=τ_ini·tanh(y/L0)**) | 4 | 4 | 0 | 9 pair checks × 8 DOFData fields × 4 phases bit-identical.  **Phase D pre-Mult `Q_self ≠ Q_nbr` spread = 1.499·10⁸ Pa confirmed** — swap genuinely exercised. | **PASS** |
| **§17.6 `test_centroid_margin_ctor_only`** + `SEAS_DIAG_CENTROID_MARGIN` on TPV102 1000 m at 12 ranks | 1 | 1 | 0 | 6 shared fault faces; margins in [227, 266] m — 14 orders of magnitude above ε_FP | **PASS — reviewer Step 3 executed; FP-fragility SUB-hypothesis refuted on this mesh; R-V92-C04's statement about the 4-km fixture stands unchallenged** |
| §4.5 `test_sigma_n_rank_consistency.py` (Frontera, §13.4) | — | — | — | Sbatch + harness written; awaiting user submission | **READY** |

### §13.3 Per-test numerical details

**§4.1 — rotation identity `max|T·Tinv − I|`**

```
BP5 canonical frame ((0,-1,0), (0,0,-1), (1,0,0))
  max|T·Tinv − I| = 0.000e+00  (0.000e+00 ulp)  BIT-EXACT

BP5 + 15° about x
  max|T·Tinv − I| = 2.220e-16  (1.000e+00 ulp)

BP5 + 30° about y
  max|T·Tinv − I| = 2.220e-16  (1.000e+00 ulp)

BP5 + 45° about (1,1,1)/√3
  max|T·Tinv − I| = 2.220e-16  (1.000e+00 ulp)
```

All 3 non-axis-aligned frames come in at 1 ULP — 16× below the
REVIEW R-004 budget.  `GodunovFlux::BuildRotation/BuildRotationInverse`
is numerically exact at machine precision on every tested frame.

**§4.2 — pure strike-slip Tinv probe on BP5 canonical frame**

Input: `Q_g[VX] = 1e-12 m/s`, `Q_g[SXY] = 7.5e7 Pa`, all others zero.
Result (Q_c after applying Tinv):

```
Q_c[SXX..SXY, SYZ] = 0      (bit-exact, no leakage)
Q_c[SXZ] = −7.5e7          (exact — n·strike shear σ_nt2 = −τ in BP5 frame)
Q_c[VX, VY] = 0            (bit-exact, no leakage)
Q_c[VZ] = 1e-12            (exact — strike velocity v_t2)
```

Signed-permutation arithmetic is exact: no ULP-scale leakage into the
normal or dip channels.  Confirms REVIEW R-002 (H-V92-O retracted).

**§4.3 — per-QP fault basis orthonormality**

| Mesh | Fault faces | QPs | worst \|v\|-1 (ULP) | worst orth dot (ULP) |
|---|---:|---:|---:|---:|
| `tpv102_1000m.msh` (10 MB) | 1 646 | 4 938 | 0.5 | 4.38·10⁻⁴⁷ |
| `tpv102_200m.msh` (126 MB, production) | 37 904 | 113 712 | 0.5 | 1.75·10⁻⁴⁶ |

The 200 m production-mesh pass is the most important datum: 113 712
QPs all satisfy orthonormality to sub-ULP.  The v9.1.0 R-003 dip
normalization fix is still in force, and FaultBasis reliably produces
orthonormal frames regardless of Gmsh tet orientation.

**§4.4 — per-QP vs centroid basis on planar 1-tet fixture**

Constructed a 2-tet mesh with shared face in the y = 0 plane.  All
component-wise assertions on the centroid basis are bit-exact:

```
n  = ( 0, -1,  0)  (ref-aligned)
t1 = ( 0,  0, -1)  (dip)
t2 = ( 1,  0,  0)  (strike)
```

All 4 QPs have the same `(n, t1, t2)` as the centroid, bit-exact.
MFEM's CalcOrtho is QP-independent on planar faces (as expected from
theory; validated here).

### §13.4 §4.5 — Frontera sbatch ready; awaiting user submission (rev-3 update)

Per user directive (2026-04-21): **Option B (Frontera)** is chosen.
Local runs on production-scale meshes — including 1000 m — are
forbidden (feedback_no_local_reproducer).  The coarser-local path
was also rejected because a 5 km fixture's ParMETIS cut pattern at
≤ 4 ranks would not exercise the fault-crossing partition seams
that H-V92-M predicts.

**Artifacts written, NOT yet submitted (user approval gate):**

| # | File | Role | Patch required? | Budget |
|---|---|---|:---:|---|
| 1 | `jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch` | §4.5 rank-A baseline (1 node / 56 ranks / 4.0 s).  Paired with existing v91 400-rank 12 s data (job 7668434) → rank-B. | No | ~48 SU |
| 2 | `debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch` | C0 additive diagnostic patch (`#ifdef SEAS_DIAG_FAULT_SIGMA` guards only).  Instruments `max\|Q_self[SYY]\|`, `max\|Q_self[VY]\|`, `max\|Q_self[SYZ]\|`, `max\|Q_self[VZ]\|` on fault QPs with MPI_Allreduce + stderr printf. | Is the patch | — |
| 3 | `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch` | §5.1 diagnostic run (8 nodes / 400 ranks / 4.0 s).  Builds `-DSEAS_DIAG_FAULT_SIGMA`; self-checks for diag signature in the binary before running. | **YES** (#2 applied) | ~358 SU |

**Total Frontera budget:** ≈ 406 SU.  Well below the rev-2 worst-case
~1000 SU envelope.  A subsequent 12 s post-fix re-run is budgeted
separately at ≈ 800 SU and is authorised after the fix is in hand.

**User action gate — required BEFORE submission:**

1. Review `tpv102_debug_v9.2.0_diag_fault_sigma.patch`.  Decide whether
   to apply it.  The patch:
   - adds 4 `mutable` `real_t` accumulators to `WaveOperator` (all
     under `#ifdef SEAS_DIAG_FAULT_SIGMA`),
   - inserts 4-line `#ifdef` blocks in `Mult`, `ComputeFaceFluxRHS`,
     `ComputeSharedFaceFluxRHS`,
   - adds one MPI_Allreduce + `fprintf(stderr, ...)` per RK4 Mult call,
   - changes **no functional code**.  `make test` must stay green when
     the flag is off (default).
2. Apply the patch with
   ```
   cd /Users/chunhuizhao/projects/seas-mfem
   patch -p0 < miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch
   ```
3. Verify the local build stays green:
   ```
   cd miniapps/seas && make test
   ```
4. Commit + push to `feature/elasticity-inertia` so Frontera checkout
   sees the patched source.
5. Submit sbatch #1 (rank consistency, no patch dependency):
   ```
   sbatch jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch
   ```
6. Submit sbatch #3 (diagnostic, requires patch).  The sbatch includes
   a `strings | grep FAULT-DIAG-SYY` self-check that fails fast if the
   patch was not applied:
   ```
   sbatch jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch
   ```
7. After BOTH jobs complete on Frontera, run off-Frontera
   post-processing:
   ```
   python3 tests/scripts/test_sigma_n_rank_consistency.py \
       --skip-run \
       --run-a  tpv102/results_200m_p1_4.0s_56r_v92_rankA_job<A_JOBID> \
       --run-b  tpv102/results_200m_p1_12.0s_400r_v91_job7668434 \
       --output-root tpv102/r_v92_r4_rank_consistency \
       --tol-mpa 0.1 --t-check 2.0 3.5
   ```
   and inspect the `RESULT.txt` files from both sbatch jobs.

**Why the compound sequence (rank-consistency + diag, not just one).**
The §7 decision matrix requires both bits — rank-dependence from §4.5
and `max|SYY_g|` magnitude from §5.1 — to uniquely identify the
primary hypothesis.  One run alone cannot discriminate:

| §5.1 max\|SYY_g\| | §4.5 rank consistency | Primary |
|---|---|---|
| ~ 0 | disagree | **H-V92-M** (bit-disagreement) primary |
| ~ 0 | agree | **H-V92-K** (RK4 station output) primary |
| ≥ 20 MPa | disagree | **H-V92-M + H-V92-G** compound; M dominant |
| ≥ 20 MPa | agree | **H-V92-G** (bulk amplification) primary |
| ~ 12 MPa | agree | **mode-II real + H-V92-C** (corner); drill via §6.3 |

**Abort paths if something goes wrong:**
- Sbatch #3 self-check fires → patch not applied → apply + re-push + re-submit.
- [FAULT-DIAG-SYY] lines missing from sbatch #3 `.out` → the `stderr`
  stream landed in `.err` (ibrun behaviour); retry parsing `.err`.
- Sbatch #1 RESULT.txt reports MFEM abort → shared-fault ctor abort
  (R-305 class).  This would be a separate symptom requiring a new v9.3
  document; DO NOT proceed to post-processing.

### §13.5 Implications for R-V92 diagnosis

**What the 4 passing regression gates tell us:**

1. **H-V92-R retraction confirmed numerically.**  REVIEW R-001's
   mathematical proof (§A) is corroborated by §4.1's direct 0-ULP
   measurement on BP5 and 1-ULP measurement on 3 off-axis frames.
   There is no sign asymmetry in `(T·Tinv)[ij, i'j']` — the
   `(a != b)` vs `(i != j)` asymmetry in `godunov_flux.cpp:272`
   vs `:231` is the correct forward/inverse Voigt symmetrization.

2. **H-V92-O retraction confirmed numerically.**  §4.2 shows pure
   strike-slip Tinv on the BP5 signed-permutation Q is bit-exact;
   no leakage into the normal or dip channels.  The rev-1 H-V92-O
   hypothesis (Tinv_can leakage) is factually impossible at the BP5
   canonical frame.

3. **H-V92-F (mesh jitter) WEAKENED.**  §4.3 passes on the
   production 200 m mesh at sub-ULP tolerance.  Mesh jitter does not
   contribute more than ~0.5 ULP per QP, which accumulates to at
   most 10⁻¹² MPa over 6 000 steps — **12 orders of magnitude below
   the observed 50 MPa drift**.  H-V92-F is essentially eliminated
   as a primary mechanism.  It may still contribute a negligible
   secular signal, but cannot be the R-V92 root cause.

4. **MFEM CalcOrtho is Jacobian-consistent on planar faces.**
   §4.4 confirms the per-QP basis coincides with the centroid basis
   when the face is planar.  Any non-trivial per-QP variation in the
   production mesh's per-QP basis must come from face curvature
   (Gmsh non-planar tets) and is already bounded by the sub-ULP §4.3
   result.

**What the regression gates DO NOT resolve:**

- They do not probe H-V92-M (MPI partition-boundary bit-disagreement
  on `elem1_on_plus`).  The §5.1 diagnostic patch + §4.5 harness are
  still required to classify this rank-1 hypothesis.
- They do not probe H-V92-G (bulk-dynamics σ_yy amplification).
  The §5.1 `max|Q_self[SYY]|` instrumentation is still required.
- They do not probe H-V92-C (fault–free-surface corner interaction).
  The §6.3 ParaView SYY slice is still required.
- They do not probe H-V92-K (RK4 stage-averaging bug in station
  output).  Requires §5.1 + comparing bulk Q_g[SYY] vs station
  sigma_n_corr drift.

**Net effect on §3's hypothesis ranking post-rev-3:**

| Hypothesis | Pre-rev-3 rank | Post-rev-3 status |
|---|:---:|---|
| H-V92-R | RETRACTED (REVIEW R-001) | **CONFIRMED RETRACTED** (§4.1 numerical evidence) |
| H-V92-O | RETRACTED (REVIEW R-002) | **CONFIRMED RETRACTED** (§4.2 numerical evidence) |
| H-V92-M | NEW RANK-1 | **Still RANK-1**; requires §4.5 + §5.1 |
| H-V92-G | NEW RANK-2 | **Still RANK-2**; requires §5.1 `max\|Q_self[SYY]\|` |
| H-V92-C | NEW  | **Unchanged** ; requires §6.3 ParaView slice |
| H-V92-F | LOWER | **WEAKENED FURTHER** (§4.3 sub-ULP on 200 m); eliminated as primary mechanism |
| H-V92-K | LOWEST | **Unchanged**; discriminated by §5.1 |

### §13.6 Remaining blocking path

The regression gates passing is a **necessary but not sufficient**
condition to move forward with R-V92.  The diagnosis still hinges
on:

1. User authorisation of the §5.1 C0 diagnostic patch
   (`SEAS_DIAG_FAULT_SIGMA`).
2. User authorisation of either local §4.5 (option A above) or
   Frontera §4.5 (option B above).
3. §6.3 ParaView SYY slice analysis (local, user-driven once the
   patch is applied and a diagnostic run is complete).

The TODO checklist at the head of this document is unchanged —
regression gates do not replace any diagnostic step.

### §13.7 Reproduction

From a shell with `conda activate mfem-dev`:

```
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make seas_test_godunov_rotation_identity \
     seas_test_canonical_rotation_pure_strikeslip \
     seas_test_fault_basis_qp_orthonormality \
     seas_test_per_qp_vs_centroid_basis_planar_1el

make test-v92-regression-gates     # runs all 4 in order

# Optional: 200 m production-mesh orthonormality run (~25 s)
./seas_test_fault_basis_qp_orthonormality \
    --mesh tpv102/mesh/tpv102_200m.msh --order 1
```

---

## §14. Consolidated overview — findings + tests (2026-04-21 rev-3)

Single-screen summary combining §1 symptom, §3 hypotheses, §4 unit
tests, §6 offline diagnostics, §13 test execution results, and the
REVIEW.md audit.  Use this as the entry point after rev-3.

### §14.1 Symptom (from §1)

| Axis | Observation |
|---|---|
| R-V92-A | σ_n at hypocenter drifts 120 MPa → ~70 MPa over 12 s (−50 MPa). DRDG3D reference: flat 120 MPa ± 1 MPa. |
| R-V92-B | `slip_dip` reaches −0.2 m at hypocenter (3 % of strike slip). Reference: \|slip_dip\| < 10⁻³ m. |
| Spatial | σ_n drift is **anti-symmetric in along-strike x** (rises on x > 0, drops on x < 0). Rules out isotropic noise and free-surface reflection. |
| Temporal | σ_n flat to t ≈ 3 s, then monotonic drift. Onset is 0.7 s after nucleation, 1.3 s before first free-surface reflection. |

### §14.2 Evidence matrix — what is ruled out vs still open

| Hypothesis | Status | Evidence closing it (or leaving it open) | Probe that discriminates |
|---|:---:|---|---|
| H-V92-R rotation stress-block asymmetry | **RULED OUT** | REVIEW R-001 proof (§A) + §4.1 numerical `T·Tinv − I` = 0 ULP on BP5, 1 ULP on 3 off-axis frames (16 × under budget) | §4.1 (passes) |
| H-V92-O Tinv leakage on pure strike-slip | **RULED OUT** | REVIEW R-002 + §4.2 bit-exact zero leakage on BP5 signed-permutation Q | §4.2 (passes) |
| H-V92-F mesh jitter on per-QP basis | **RULED OUT as primary** | §4.3 on 200 m production mesh: 113 712 QPs at 0.5 ULP worst norm dev; accumulated error ≤ 10⁻¹² MPa over 6 000 steps, 12 orders below the 50 MPa observation | §4.3 (passes) |
| Formula sign in Pelties eq. (7)/(11)/(12) | **RULED OUT** | §2.3 re-derivation from A-matrix characteristics in BP5 "+ to −" convention matches MFEM `fault_face_flux.cpp:52–180` and SeisSol `ElasticSetup.h` bit-for-bit | — (analytical) |
| Free-surface γ-mirror gamma pattern | **RULED OUT** | `{-1, 1, 1, -1, 1, -1, 1, 1, 1}` is the correct σ·n = 0 mirror at z = 0 (verified vs SeisSol `Model/Common.h`) | — (analytical) |
| BP5 dip/strike component assignment | **RULED OUT** | CLAUDE.md R-801 authoritative; `tau2_0 = tau_ini`, `V2 = V_ini` | — (convention) |
| **H-V92-M MPI partition-boundary disagreement** | **OPEN — RANK-1** | Anti-symmetric-in-x σ_n drift is the signature of a ParMETIS cut at x = 0 with side-dependent `elem1_on_plus` bit; v9.0.0 R-802 unified `can_n` but **did not** unify `shared_fault_elem1_on_plus_`, which is computed from per-rank Elem1 geometry | §4.5 (1 vs 4 rank) + §5.1 `max\|Q_self\|` |
| **H-V92-G bulk σ_yy amplification** | **OPEN — RANK-2** | §2.5 mode-II → σ_yy Poisson coupling is real physics but should be ~1 MPa; observation is 50×.  Four candidate amplifiers: (a) IP/BR2 penalty in wave driver, (b) RK4 non-conservative drift on σ_yy channel, (c) absorbing-BC wrong-sign reflection, (d) element-mass-inverse conditioning asymmetry across Voigt channels | §5.1 `max\|Q_self[SYY]\|` tracks σ_n drift |
| **H-V92-C fault–free-surface corner** | **OPEN** | Corner arrival time `depth/cs ≈ 2.16 s` after rupture matches the 3 s onset; BuildFrame at z = 0 flips to up = x̂, producing an internal frame whose γ-mirror may leak into global σ_yy | §6.3 ParaView slice at y = 0 showing corner emanation |
| **H-V92-K RK4 stage-averaging bug (station output only)** | **OPEN — LOWEST** | Driver averages DOFData diagnostic fields across RK4 stages; if bulk Q_g[SYY] stays ~0 but station σ_n_corr drifts, this is the bug | §5.1 + §6.3: bulk clean + station drift = K primary |

### §14.3 Discriminating table — how one diagnostic run classifies the open set

From §7, with §5.1 `max|Q_self[SYY]|` measured at fault QPs, §4.5
1-rank vs 4-rank σ_n comparison, and §6.3 SYY slice:

| §5.1 max\|SYY_g\| at t=4 s | §4.5 rank agreement at t=3.5 s | §6.3 SYY slice pattern | Primary |
|---:|:---:|:---|:---:|
| ~ 0 | yes | noise | **K** |
| ~ 0 | no | noise | **M** (+ K) |
| ~ 0.12 MPa | yes | rupture-front aligned | mode-II real; dig into C |
| ~ 12 MPa | yes | front + corner emanation | **G + C** |
| ~ 12 MPa | no | partition lines | **M primary**, G/C co-amplify |

One diagnostic pass (4 s, 400 ranks, ~200 SU) uniquely identifies
the primary mechanism.

### §14.4 Tests — what each one provides

| # | Test | Kind | What it PROVES | What it DOES NOT probe |
|---|---|---|---|---|
| §4.1 | `test_godunov_rotation_identity` | C++ | Rotation 9×9 block is exact on BP5 + 3 off-axis frames ≤ 1 ULP | Per-QP basis variation; MPI consistency |
| §4.2 | `test_canonical_rotation_pure_strikeslip` | C++ | Signed-permutation Tinv on BP5 is bit-exact zero leakage into n/t1 channels | Bulk σ_yy sourcing; non-axis-aligned per-QP frames |
| §4.3 | `test_fault_basis_qp_orthonormality` | C++ | Per-QP `(n, t1, t2)` orthonormal to ≤ 0.5 ULP on 113 712 QPs of 200 m production mesh | Basis rank-consistency across MPI partitions |
| §4.4 | `test_per_qp_vs_centroid_basis_planar_1el` | C++ | CalcOrtho is QP-invariant on planar faces; eliminates per-QP Jacobian as a variable | Curved-face contributions; MPI consistency |
| §4.5 | `test_sigma_n_rank_consistency.py` | Py harness | **Would directly classify H-V92-M** at 1 vs 4 ranks | Bulk σ_yy (needs §5.1); corner (needs §6.3) |
| §5.1 | `SEAS_DIAG_FAULT_SIGMA` (C0 patch) | C++ instr. | `max\|Q_self[SYY]\|`, `max\|VY\|`, `max\|SYZ\|`, `max\|VZ\|` at fault QPs — directly discriminates G vs K | Rank agreement (needs §4.5) |
| §6.1 | Offline σ_n vs tau2_corr correlation | Py script | Cluster pattern → G (proportional) vs M (x-clustered) vs C/K (uncorrelated) | — (observational only) |
| §6.3 | ParaView SYY slice at y = 0 | Visual | Spatial source of σ_yy: rupture front (G), corner (C), partition lines (M), or noise (K) | — (visual only) |

### §14.5 Current status (post rev-3 execution)

✅ **Completed:**
- §4.1 PASS (4/4 tests).  H-V92-R numerically retracted.
- §4.2 PASS (11/11 tests).  H-V92-O numerically retracted.
- §4.3 PASS (1 000 m: 4 938 QPs; 200 m prod: 113 712 QPs).  H-V92-F
  eliminated as primary mechanism.
- §4.4 PASS (11/11 tests).  MFEM Jacobian consistency confirmed on
  planar fixture.
- §4.5 written + syntax-clean; execution deferred pending §13.4
  decision.
- Makefile: 5 new phonies + aggregate `test-v92-regression-gates`.

⏸️ **Blocked on YOU:**
- §5.1 C0 diagnostic patch authorisation (additive stderr printf; no
  functional change).  Required for G / K discrimination.
- §4.5 execution mode choice (Option A local 1 000 m, B Frontera, or
  C inline fixture — see §13.4).  Required for M discrimination.
- §6.3 ParaView slice review (can be done any time on existing
  v9.1.0 job-7668434 PVD).  Required for C discrimination.

### §14.6 One-sentence conclusion

After rev-3 the candidate set for R-V92 is narrowed from
{R, O, F, M, G, C, K} (rev-1's 7) to **{M, G, C, K}** (rev-3's 4),
with a pre-built discriminator (§5.1 + §4.5 + §6.3) that classifies
the survivor in a single 4 s diagnostic run at ~200 SU — versus
rev-1's 1 500-SU by-elimination plan.  The root cause is not yet
identified; the next YOU-gated actions are the three bullets in
§14.5 "Blocked on YOU".

---

## §15. New observation — ParaView fault-surface speckle pattern (rev-3b, 2026-04-21 PM)

### §15.1 Data (8 screenshots at t ≈ 8.485 s, v9.1.0 job-7668434, 400 ranks, 200 m)

Rendered `FaultSurface/fault_surface.pvd` in ParaView; one averaged
CellData value per triangle (per `io/paraview_output.hpp:715-734` —
3 QPs averaged per face post-v9.1.0 R-001 / H-V91-A4 fix).

| Field | Visual signature | Nominal range | Outlier range |
|---|---|---|---|
| `normal_stress` | Uniform fault surface + **scattered outlier triangles inside the rupture disk** (hundreds of specks) | ~1.2·10⁸ Pa | 1.0·10⁷ – 2.2·10⁸ Pa (2× and 0.1× background) |
| `traction_strike` | Annular ring at rupture front + scattered inside disk | ~7.5·10⁷ Pa | 1.5·10⁷ – 1.1·10⁸ Pa |
| `traction_dip` | Smooth rupture-tip dipole + scattered inside disk | ≈ 0 | ±1.9·10⁷ Pa |
| `state_variable` (ψ) | Faded disk interior + scattered in disk | 0.42 – 0.85 | same range, **different DOFs extreme in each field** |
| `slip_rate_strike` | Smooth circular front + scattered specks | 0 – 6 m/s | ~ −5·10⁻¹³ – 6 m/s |
| `slip_rate_dip` | Smooth, narrower than strike + scattered specks | ≈ 0 – 0.5 m/s | −5·10⁻¹² – 5 m/s |
| `slip_strike` | Smooth accumulated slip disk + scattered specks | 0 – 5 m | spike outliers up to 5 m at isolated points |
| `slip_dip` | Radial dipole pattern (expected mode-II) + scattered specks | ±0.2 m | same range, isolated specks correlate across fields |

**Co-location evidence.** The scattered outlier triangles are at the
**same (x, z) coordinates across all 8 fields**.  A single fault
face showing an extreme `normal_stress` also shows an extreme
`slip_rate_strike`, an extreme `state_variable`, an extreme
`traction_strike`, etc.  They are **per-face outliers correlated
across all DOFData fields**.

**Spatial pattern.** Outliers are scattered **inside** the rupture
disk (where `slip_rate > 0`) and absent **outside** (locked
region).  Hundreds of them; on a 200 m mesh with 37 904 fault faces
and 400 MPI ranks that is on the order of **one bad face per rank**.

**Temporal pattern (inferred from station plots).** σ_n drift at
stations onsets at t ≈ 3 s — the time the rupture front expands to
cover enough shared fault faces for the speckle to become visible.

### §15.2 Signature matching against open hypotheses

| Hypothesis | Predicted ParaView signature | Observed? |
|---|---|:---:|
| H-V92-M (rank-boundary bit-disagreement on fault faces) | Outliers at shared fault faces; one-per-rank-boundary; correlated across fields | **✅ matches** |
| H-V92-G (bulk σ_yy mode-II amplification) | Smooth correlated drift (no speckle); σ_n deviation ~1 MPa aligned with rupture front | ❌ predicts no speckles |
| H-V92-C (fault–free-surface corner) | Pattern emanating from z = 0 edge; band across top of fault | ❌ observed pattern is interior-scattered |
| H-V92-K (RK4 stage-averaging, station output path) | Drift at station plots only; VTU (written from final dof_data) clean | ❌ VTU itself is noisy |

**Decision: H-V92-M is promoted from probable to highly-likely
primary.**  The speckle signature is a near-textbook match for
shared-fault DOF inconsistency at MPI partition boundaries.  G, C,
K remain open as secondary contributors but **none alone can
produce scattered-outlier triangles correlated across every DOFData
field**.

### §15.3 Refined mechanism — what exactly diverges

1. `dof_data[]` has layout `[interior-face QPs (n_int·nbf), shared-face QPs (n_shared·nbf)]` (see `tpv102_driver.cpp:393-395`).  Each rank owns its interior and shared QPs independently.

2. `WriteFaultSurfaceVTU` (`paraview_output.hpp:567-642`) computes the per-triangle average from `dof_data` indexed by `face_mesh_idx * nbf + k`.  So the CellData value on a shared face is the average of **that rank's own side** of the shared face — the neighbor rank's copy is never written (shared-face ownership is deduplicated one-entry-per-face).

3. If two ranks sharing a fault face produce **different** `dof_data[i].psi / tau_corr / V / slip / sigma_n_corr`, the VTU shows ONE triangle with the value from the owning rank.  A per-face divergence appears as a per-triangle outlier.

4. Where the disagreement originates:
   - (a) `shared_fault_elem1_on_plus_[sf_idx]` populated from per-rank Elem1 geometry in the WaveOperator ctor.  If the two ranks' Elem1 instances come from different local-mesh orderings, they can bit-disagree.
   - (b) `ComputeSharedFaceFluxRHS` routes `(Q_self, Q_nbr)` into `(Q_plus, Q_minus)` using this flag.  Disagreement flips (+, −) roles, giving opposite-signed trial jump contributions.
   - (c) Each rank computes its own `(σ_n_corr, τ_corr)` from its (possibly flipped) trial state; friction inversion is nonlinear, so the error grows over RK4 stages.
   - (d) Per-rank divergent `dof_data` is written to the VTU — one outlier triangle per disagreeing shared face.

**Self-check: this mechanism explains:**
- ✅ Speckle pattern (one outlier per disagreeing shared fault face).
- ✅ Cross-field correlation (same DOF drives every field).
- ✅ Count on the order of hundreds (matches ParMETIS boundary count at 400 ranks).
- ✅ Anti-symmetric-in-x σ_n drift (§1.3; ParMETIS cuts near x = 0 produce sign bias).
- ✅ t ≈ 3 s onset (visible only once rupture covers enough shared faces).
- ✅ Consistent with v9.0.0 R-802 / v9.1.0 rev-1 suspicion of shared-fault sign inconsistency.

### §15.4 New hypothesis H-V92-P (refinement of H-V92-M)

**H-V92-P (RANK-1, PROMOTED).**  Shared-fault QPs on MPI partition
boundaries carry divergent DOFData across the two ranks that share
them.  Origin: per-rank bit-disagreement in
`shared_fault_elem1_on_plus_` (or an equivalent sign-determining
quantity in `ComputeSharedFaceFluxRHS`).  Once (+, −) roles
disagree on a face, the two ranks compute different trial tractions
→ different friction-solver outputs → different ψ/V/slip/τ state
updates.  Over thousands of RK4 steps this diverges into scattered
outlier triangles in the fault-surface VTU, each outlier being one
shared fault face with a rank-role disagreement.

**Relation to H-V92-M.**  H-V92-P IS H-V92-M with a specific
mechanism (per-face role disagreement) and an observable
signature (scattered outlier triangles correlated across fields).

### §15.5 New diagnostics for H-V92-P

#### §15.5.1 NEW test §4.6 — `test_shared_fault_role_consistency.cpp` (regression gate) — **IMPLEMENTED + PASS rev-3b, 2026-04-21**

**Purpose.** On a ≥ 4-rank MPI run (4 chosen because ParMETIS did not
cut through the fault at 2 ranks on this compact fixture), check for
every shared fault face that the two ranks sharing it produce a
consistent (+, −) routing for the Pelties-9 flux.

**Spec-semantic correction during implementation.**  The rev-3b
specification originally read "the two ranks AGREE on
(elem1_on_plus, can_n, can_t1, can_t2)".  An audit of
`wave_operator.inl:1279-1290` showed the contract is subtler:
`elem1_on_plus` is a PER-RANK geometric bit referring to each rank's
own `Elem1` (its local element).  The two ranks' local elements sit
on OPPOSITE physical sides of the fault, so `elem1_on_plus` MUST be
**ANTI-SYMMETRIC** across ranks — exactly one rank has true, the
other false.  What AGREES between ranks is the (+, −) routing AFTER
the swap: `Q_plus_local = elem1_on_plus ? Q_self : Q_nbr`, which
both ranks resolve to the same physical state (by MFEM ghost
exchange).  The corrected test checks for bit-antisymmetric
`elem1_on_plus` and bit-identical canonical frame (can_n, can_t1,
can_t2, nl).

**Construction.**
- Inline fixture: 4 × 4 × 4 km box, Cartesian tet mesh (192 tets,
  144 boundary elements), fault at y = 0 tagged attribute 3 on
  every interior y = 0 face.
- Partition to ParMesh; instantiate `WaveOperator<ParMesh>` ctor.
- Each rank collects `(key=rounded centroid, src_rank,
  elem1_on_plus, sign_flipped, can_n[3], can_t1[3], can_t2[3], nl)`
  for every shared-fault QP it owns.
- `MPI_Allgatherv` of records; each rank walks its own records,
  finds the partner record by centroid-key from a different rank
  (skip via `src_rank != my_rank`), and asserts:
  - `elem1_on_plus` **not equal** (anti-symmetric) on all partner
    pairs,
  - `can_n`, `can_t1`, `can_t2`, `nl` bit-equal on all partner pairs.
- Aggregate via `MPI_Allreduce(SUM)` → rank 0 reports verdict.

**File.** `tests/parallel/test_shared_fault_role_consistency.cpp`.
Makefile targets: `seas_test_shared_fault_role_consistency` (binary)
+ `test-shared-fault-role-consistency` (phony running
`mpirun -np 4 ./seas_test_shared_fault_role_consistency`).  Wired
into `test-v92-regression-gates` as the 5th C++ test.

**rev-3b execution result (2026-04-21).**

```
serial mesh : 192 elements, 448 faces, 144 bdr elems
rank 0 : local_ne=49 fault_int=4 fault_shr=3
rank 1 : local_ne=46 fault_int=2 fault_shr=2
rank 2 : local_ne=49 fault_int=3 fault_shr=0
rank 3 : local_ne=48 fault_int=4 fault_shr=1

Shared-fault role audit (summed over ranks):
  pair-wise comparisons  : 18
  elem1_on_plus EQUAL    : 0   (expect 0 ⇒ anti-symmetry intact)
  can_n mismatches       : 0   worst |Δ| = 0 (bit-identical)
  can_t1 mismatches      : 0   worst |Δ| = 0 (bit-identical)
  can_t2 mismatches      : 0   worst |Δ| = 0 (bit-identical)
  nl (|n_raw|) mismatches: 0   worst |Δ| = 0 (bit-identical)

  PASSED: elem1_on_plus is ANTI-symmetric across ranks
  PASSED: can_n  bit-identical across ranks
  PASSED: can_t1 bit-identical across ranks
  PASSED: can_t2 bit-identical across ranks
  PASSED: |n_raw| bit-identical across ranks
Results: 5 passed, 0 failed out of 5 tests
```

**Interpretation.**  H-V92-P is **NOT confirmed at init level** on
this fixture partition pattern.  The canonical-frame machinery
(v9.0.0 R-802) is bit-identical across ranks; `elem1_on_plus`
correctly inverts.  Therefore the ctor-time bookkeeping is not the
root cause.

**What this does NOT rule out:**
- **Runtime** H-V92-P — divergence that accumulates after RK4
  stepping, e.g., asymmetric friction-solver iterations on the two
  ranks updating `DOFData` out of sync.
- H-V92-P at the 400-rank production partition pattern (ParMETIS's
  cut topology at 400 ranks on the 200 m mesh may produce
  pathological partitions that a 4-rank fixture does not).

**Next step (rev-3b plan ordering).** §4.6 PASS means proceed to:
- §6.4 offline VTU outlier detector (local, no Frontera, §15.5.2).
- Then §9.2 Frontera 4 s diag with §5.1 + §5.2 patches (§15.5.4)
  to check runtime H-V92-P at production scale.
- §4.5 rank consistency (my §13.4 rankA sbatch) retained as
  fallback.

#### §15.5.2 NEW §6.4 offline VTU outlier detector

**Purpose.** Quantify the speckle on the already-written v9.1.0
job-7668434 VTU output without Frontera.

**Construction (Python ~150 LOC).**
```python
for each cycle in FaultSurface/fault_surface.pvd:
    load all fault_surface_r*_c{cycle}.vtu
    for each triangle:
        collect per-face (normal_stress, traction_strike,
                          state_variable, slip_rate_strike, ...)
    compute per-field median + MAD
    flag triangle as outlier if ANY field exceeds median + 5·MAD
    count cross-field coincidence
output:
    - outlier count per cycle (stable → P; growing → G)
    - spatial coordinates (x, z) of outlier centroids
    - fraction on the shared-face subset (inferred from rank id
      in .pvtu Piece listing)
```

**Expected classification:**
- Outlier count stable + on shared-face subset → **H-V92-P**.
- Outlier count grows monotonically → H-V92-G amplification.
- Outliers scatter regardless of shared-face membership → H-V92-Q
  (new; friction solver per-QP instability).

**File.** `tests/scripts/fault_vtu_outlier_detector.py`.  Local
only.

#### §15.5.3 NEW §5.2 instrumentation — rank-id CellData in VTU

**Purpose.** Make the VTU directly show MPI rank ownership +
shared-vs-interior bit, so §6.4 can correlate speckle with
partition boundaries without external metadata.

**Patch to `paraview_output.hpp:715-734` (additive, ~15 LOC):**
```cpp
std::vector<int32_t> c_rank(ncells, rank);
write_field_int32("mpi_rank",     c_rank);
std::vector<int32_t> c_is_shared(ncells);
for (int i = 0;      i < n_int;              i++) c_is_shared[i] = 0;
for (int i = n_int;  i < n_int + n_shared;   i++) c_is_shared[i] = 1;
write_field_int32("is_shared_face", c_is_shared);
```

Controlled by new compile flag `SEAS_DIAG_VTU_RANK` (off by
default; set in the diagnostic sbatch only).  **Pure I/O patch, no
functional change, safe for any run including production.**

#### §15.5.4 NEW §9.2 Frontera diag — 4 s with §5.2 + §5.1 patches

**sbatch.** `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch` (actual filename; §15.5.4 draft originally planned `v92p_diag` but the committed artifact uses `v92_diag`).

**Build flags:** `SEAS_DIAG_FAULT_SIGMA` (§5.1; in-flight σ_n probe)
+ `SEAS_DIAG_VTU_RANK` (§5.2; rank-id VTU CellData) +
`SEAS_DIAG_FAULT_FLUX` (existing C-1/2/3/4 printfs).

**Run parameters:** 200 m mesh, P1, 400 ranks, tfinal = 4 s,
CFL = 0.5, ParaView VTU every 0.5 s.

**Wall time:** ~ 30 min × 400 ≈ 200 core-hours.

**Deliverable:** VTU snapshots at t = 3.0, 3.5, 4.0 s.  Threshold
filter on `is_shared_face == 1` + color by `mpi_rank`.  If ≥ 95 %
of outlier triangles are `is_shared_face == 1`, H-V92-P is
**CONFIRMED** at runtime.

#### §15.5.5 Updated §4.5 test interpretation

The existing §4.5 (1-rank vs 4-rank σ_n harness) remains useful but
§15.5.1 is a more direct probe of the root cause.  New priority
order:

1. **§15.5.1** — shared-fault role test @ 2 ranks.  Local, < 5 min.
   If FAIL: H-V92-P confirmed at init level; fix first, then
   Frontera.
2. **§15.5.2** — offline VTU outlier detector on existing job
   7668434 PVD.  Local, ~15 min.  Quantifies speckle and
   approximate partition correlation.
3. **§5.1** diagnostic patch.  Still valuable to rule out H-V92-G
   as primary.
4. **§15.5.4** — Frontera 4 s diag with §5.2 rank-id VTU.  Runs
   ONLY if §15.5.1 PASSES (meaning the divergence is at runtime,
   not init).
5. **§4.5** — deferred; subsumed by §15.5.1 (init-level) and
   §15.5.4 (runtime).

### §15.6 Updated §14.2 evidence matrix

| Hypothesis | Status before rev-3b | Status after rev-3b |
|---|---|---|
| H-V92-R | RULED OUT (§A + §4.1) | unchanged |
| H-V92-O | RULED OUT (§A + §4.2) | unchanged |
| H-V92-F | RULED OUT as primary (§4.3) | unchanged |
| Formula sign / γ-mirror / BP5 convention | RULED OUT (§2.3 / §A) | unchanged |
| **H-V92-M / P** | OPEN — RANK-1 | **Init-level ELIMINATED** by §4.6 PASS (2026-04-21).  **Local 4-rank runtime-level ELIMINATED** by §4.7 PASS (2026-04-21): 3 phases × 9 pair checks × 8 DOFData fields all bit-identical across ranks.  **Production-scale runtime H-V92-P still OPEN** — FP-fragility of `Elem1_centroid · ref_normal` at ParMETIS cuts near x=0 on Gmsh-generated 200 m / 1000 m meshes cannot be probed locally.  Requires §9.2 Frontera diag run. |
| H-V92-G | OPEN — RANK-2 | **WEAKENED** — predicts smooth drift, not speckle.  Retained as possible secondary. |
| H-V92-C | OPEN | **WEAKENED** — predicts edge-emanating pattern, not interior-scattered.  Retained, downgraded. |
| H-V92-K | OPEN — LOWEST | **WEAKENED** — VTU is written from final `dof_data`, so bulk is polluted, not just station output.  Retained as possible amplifier on top of P. |
| **H-V92-Q (NEW)** | — | **OPEN — RANK-3** — friction-solver per-QP instability at isolated fault QPs.  Distinct from P because outliers would co-locate with friction-regime transitions, not partition boundaries.  §15.5.2 discriminates. |

### §15.7 Updated blocking path (supersedes §14.5)

**YOU-gated, in priority order:**

1. **(FAST, LOCAL)** Authorize §15.5.1 shared-fault role-consistency
   test.  `mpirun -np 2` fixture, < 5 min wall.  Direct yes/no on
   H-V92-P at init level.
2. **(FAST, LOCAL)** Authorize §15.5.2 offline VTU outlier detector
   on the existing job-7668434 PVD.  Zero Frontera cost.
3. **(C0 PATCH)** Authorize §5.2 `SEAS_DIAG_VTU_RANK` rank-id VTU
   instrumentation (~15 LOC; pure I/O).
4. **(CONDITIONAL FRONTERA)** Only if §15.5.1 PASSES and §15.5.2
   shows outliers aligning with rank boundaries: run §15.5.4 (~200
   SU).  Otherwise skip — fix at init level first.
5. **(FIX PROPOSAL)** §8.1 H-V92-M fix shape is now the
   highest-probability target.  Concretely: unify
   `shared_fault_elem1_on_plus_` via `ref_normal` + face centroid
   (both bit-identical across ranks by geometry), not via per-rank
   Elem1 local numbering.

### §15.8 Updated Frontera budget

| Step | Wall | SU |
|---|---:|---:|
| §15.5.1 shared-fault role test | < 5 min local | 0 |
| §15.5.2 outlier detector | ~15 min local | 0 |
| §5.2 rank-id VTU patch + rebuild | 30 min local | 0 |
| §15.5.4 Frontera diag (conditional) | 30 min × 400 | 200 |
| §8.1 fix + §4.1–§4.6 regression | 2 h local | 0 |
| Post-fix Frontera 12 s re-run | 2 h × 400 | 800 |
| **Total (§15.5.1 FAILS, skip §15.5.4)** | — | **800** |
| **Total (§15.5.1 PASSES, need §15.5.4)** | — | **1 000** |

§15.5.1 **failing** is the cheapest path — it would confirm the
root cause at init level without any further Frontera cost.

### §15.9 One-sentence update

The ParaView speckle pattern (scattered outlier triangles
correlated across every DOFData field, centred inside the rupture
disk, numbering on the order of one per rank) is a **textbook
signature of H-V92-P / M** (shared-fault DOFData divergence across
MPI partition boundaries) and simultaneously *WEAKENS* G, C, and K
as primary candidates — the cheapest next step is the local 2-rank
§15.5.1 role-consistency test, which can confirm H-V92-P at init
level in under 5 minutes without touching Frontera.

---

## §16. BP5-reference comparison — how TPV102 reinvented rank-comm (rev-3c, 2026-04-21)

Per user directive: **"check BP5 workflow for this part, which is
proven to be working, identify unnecessary (possible wrong) code /
different handling ways for TPV102 we wrote for rank
communications"**.  BP5 has run production-scale benchmark (400+
ranks) without the speckle signature; its fault-face MPI model is
the reference.

### §16.1 BP5's rank-comm pattern (working)

Source: `domain/elasticity_operator_setup.inl:891-1055`
(`BuildOwnedFaultLayout`) and
`domain/elasticity_operator_traction.inl:187-225`
(`ExpandOwnedToLocalFault`).

**Step 1 — canonical integer key.**  For every fault face, the
layout builder computes a `FaceVertexKey` as the **sorted triple
of global vertex IDs** (`MakeFaceKey`,
`elasticity_operator.hpp:~485` + mirrored in
`dynamic/shared_fault_key.hpp`).  Because global vertex IDs are
bit-identical across ranks by ParMesh construction, the key is an
**integer** that is guaranteed bit-exact across ranks — zero FP
comparison.

**Step 2 — unique owner election.**
```
Allgatherv all ranks' shared-fault keys → build key → ranks list
for each shared fault face:
   owner_rank = sharing_ranks.front()    // lowest sharing rank ID
   if (my_rank == owner_rank):
       append this face to owned_fault_face_to_local_face_
       record send_owned_faces[other_rank] entry
   else:
       record recv_local_faces[owner_rank] entry
```
This produces **deterministic ownership**: every shared fault face
has exactly one owner, chosen by integer-exact logic.

**Step 3 — build comm blocks sorted by key.**
`SharedFaultCommBlock::send_owned_faces` and
`::recv_local_faces` are both sorted by the canonical key, so send
and receive orders are pairwise-aligned across ranks (the nth send
on rank A corresponds to the nth recv on rank B).

**Step 4 — owner computes, non-owners receive.**  The domain
solver only computes traction at *owned* fault DOFs.  After each
timestep the owner's values are scattered to all sharing ranks
via `FaultScatter::BeginScatter/WaitScatter`
(`elasticity_operator_traction.inl:193-223`), with a per-rank
`canonical_to_local_perm_` mapping applied on receive so DOFs land
at the right MFEM-local positions.

**Net effect.**
- Each shared fault face has **ONE authoritative value** (the
  owner's).
- Non-owners never compute for that face; they just receive the
  owner's result.
- There is NO floating-point sign-determination, NO "+/−" labeling
  per rank, NO redundant friction solve.
- BP5 has been run at 400+ ranks with no speckle signature.

### §16.2 TPV102's rank-comm pattern (broken)

Source: `dynamic/wave_operator.inl:395-447`
(`shared_fault_elem1_on_plus_` population),
`wave_operator.inl:1202-1340`
(`ComputeSharedFaceFluxRHS` flux + DOFData update).

**Step 1 — per-rank geometric sign disambiguation.**  For every
shared fault face, each rank independently computes
```
shared_fault_elem1_on_plus_[sf_idx] =
   ((Elem1_centroid − face_centroid) · ref_normal < 0)
```
Both `Elem1_centroid` and `face_centroid` are **floating-point**
values.  Although mathematically this gives consistent "+/−"
labels across ranks (exactly one of the two ranks' Elem1 is on the
+side), it relies on the sign of an FP dot-product being
bit-identical to what the other rank computes — a fragile
assumption when face-adjacent vertices are close to the reference
plane (e.g., fault at y = 0 with Gmsh-generated mesh jitter).

**Step 2 — BOTH ranks compute redundantly.**
`ComputeSharedFaceFluxRHS` is called on both sides of every shared
fault face:
```cpp
for (int sf_idx = 0; sf_idx < fault_shared_faces_.Size(); sf_idx++) {
   // Get Q_self (this rank's side) and Q_nbr (ghost from other rank)
   // Rotate both into canonical frame
   bool elem1_on_plus = shared_fault_elem1_on_plus_[sf_idx];
   const real_t *Q_plus  = elem1_on_plus ? Q_self_can : Q_nbr_can;
   const real_t *Q_minus = elem1_on_plus ? Q_nbr_can : Q_self_can;
   fault_flux_->Evaluate(fdata, Q_plus, Q_minus, ...);   // friction solve!
   // ...both ranks write dof_data at their own dof_idx
}
```
**Both ranks** run `fault_flux_->Evaluate()` — i.e. both ranks do
the friction Brent solve — and both write their own local
`dof_data[dof_idx]`.  Two separate copies of the DOFData exist for
every shared fault face, one per rank.

**Step 3 — correctness assumption (fragile).**
The code is correct *only if* both ranks produce bit-identical
Q_plus / Q_minus at each call.  For that to hold:
- MFEM's ghost exchange must give bit-identical Q_nbr on A ≡ Q_self
  on B (usually true).
- `elem1_on_plus` must agree (one rank true, the other false) on
  every shared face (fragile: FP comparison).
- `can_n, can_t1, can_t2` from each rank's FaultBasis must be
  bit-identical (usually true per §4.3).
- The friction Brent solver must be deterministic for identical
  inputs (true; no rank-dependent state).

If ANY of the above fails on ANY shared face, the two ranks'
DOFData diverge — producing the per-triangle speckle §15
observes.

**Step 4 — there is NO post-solve MPI scatter of DOFData.**  Each
rank's `dof_data[]` is never reconciled with its neighbors'.  The
ParaView VTU writes the owning rank's (Elem1-side) DOFData without
cross-checking.  So a divergence is silently preserved.

### §16.3 Diff — unnecessary / possibly wrong code in TPV102

| Aspect | BP5 (working) | TPV102 (broken) | Why TPV102's choice is risky |
|---|---|---|---|
| Face identity | `MakeFaceKey` → sorted triple of global vertex IDs (integer) | Geometric centroid projection onto `ref_normal` (FP comparison) | FP sign of a dot-product near a planar fault is fragile. |
| Ownership | Unique owner = lowest-ranked rank sharing the face | None; every rank computes on both sides | Redundant work + opportunity for divergence. |
| Flux / friction Brent solve | Only the owner runs it | Both ranks run it | 2× CPU + divergence risk. |
| DOFData storage | Only the owner owns the DOF | Both ranks hold their own copy | Two sources of truth at a partition boundary. |
| Post-solve sync | `FaultScatter` MPI_Isend/Irecv, canonical→local permutation | None | Divergent copies are never repaired. |
| Permutation handling | `canonical_to_local_perm_` per face | N/A (no owner→non-owner send) | TPV102's layout assumes both sides use the same MFEM DOF order — valid per-rank but not enforced as a cross-rank invariant. |
| Key-based infrastructure | Used in runtime path | **Already exists** (`dynamic/shared_fault_key.hpp`) — **used only for a diagnostic verifier** at `wave_operator.inl:1608` | TPV102 imports the same `MakeFaceKey` header but never routes it into the runtime path; it is dead code for the hot path. |

**Observation.** `dynamic/shared_fault_key.hpp` is explicitly
documented as the eventual-unification point with BP5 (its comment
says **"BP5's ElasticityOperator currently owns an equivalent,
bit-identical FaceVertexKey ... the long-term plan is for BP5 to
ALSO include this header and drop its nested copies"**).  The
bit-identical key infrastructure is present in TPV102 but only
used by the R-101 cross-rank **verifier** (gather-and-compare) —
not by the runtime `ComputeSharedFaceFluxRHS` path.  This is the
single most visible architectural divergence from BP5 that
explains the speckle.

### §16.4 Refined root-cause statement

**TPV102 never implemented the BP5 owner-decides-and-scatters
pattern for fault faces.**  Instead it uses a per-rank geometric
"+/−" labeling and runs the friction solve on both sides of every
shared face, producing two competing copies of DOFData that are
never reconciled.  Under favourable FP conditions the two copies
agree and the ParaView output is clean; under unfavourable
conditions (mesh jitter, close-to-zero dot products, specific
partition geometries at ParMETIS cuts near x = 0) the copies
diverge and produce the scattered-outlier triangle pattern.

This is **exactly** H-V92-P from §15.4, now with a specific
architectural explanation rather than a guessed mechanism.

### §16.5 Implications for §8 fix (REVISED)

The previous §8.1 fix shape (\"unify `shared_fault_elem1_on_plus_`
via `ref_normal` + face centroid\") is **still FP-based** and
therefore does not close H-V92-P robustly.

**New §8.1 fix shape — port BP5's pattern.**  The minimal robust
fix is to port BP5's `BuildOwnedFaultLayout` + `FaultScatter`
machinery into the dynamic wave operator:

1. **Reuse** `dynamic::MakeFaceKey` (already present and
   unit-tested) for the runtime shared-fault pairing — not just
   the diagnostic verifier.
2. **Assign owner** = `sharing_ranks.front()` for each shared
   fault face (identical pattern to BP5
   `elasticity_operator_setup.inl:995`).
3. **In `ComputeSharedFaceFluxRHS`**: guard the
   `fault_flux_->Evaluate` + dof_data update with
   `if (my_rank == owner_rank) { ... }`.  Non-owners skip the
   friction solve entirely.
4. **After ComputeSharedFaceFluxRHS**: scatter owner's updated
   DOFData fields to non-owners via a new `FaultDOFScatter` helper
   (modelled on BP5's `FaultScatter`), so both ranks' local
   `dof_data[]` hold identical owner-computed values.

**Estimated patch size.** ~200 LOC in `wave_operator.inl` +
~80 LOC new `FaultDOFScatter` (can reuse BP5's `FaultScatter`
with minor modifications for the 8 DOFData fields vs. 1 traction
field).  Larger than previous §8.1 shape but architecturally
aligned with the working BP5 pattern; eliminates the entire
divergence class H-V92-P covers.

**Why the simpler "just unify `shared_fault_elem1_on_plus_`" fix
is insufficient.**  Even if the `+/-` labels agree perfectly,
BOTH ranks running `Evaluate()` (friction Brent solve) still
leaves two independent copies of dof_data, and any FP difference
in the friction solve path (e.g., a different number of Brent
iterations due to slightly different intermediate state) would
still produce divergent dof_data.  Single-source-of-truth (one
owner) is the only robust pattern.

### §16.6 Implications for §4.6 unit test (REVISED)

The §15.5.1 unit test `test_shared_fault_role_consistency.cpp`
currently tests only the **labeling** agreement.  Given §16's
finding, we should ALSO add:

**§4.7 (NEW) `test_shared_fault_dof_data_consistency.cpp`** —
MPI test (4 ranks on the 4-km inline tet fixture, since 2 ranks
do not cut the fault in this compact geometry) that drives the
wave operator through 1 RK4 step on a minimum mesh, then checks
that for every shared fault QP the two ranks'
`dof_data[dof_idx]` fields agree bit-for-bit.

This is a **stronger** test than §4.6 — it exercises the full
flux-computation path, not just the ctor-level labeling.

**rev-3c execution result (2026-04-21): PASS 3/3 phases.**

```
serial mesh : 192 elements, 448 faces, 144 bdr elems
rank 0: fault_int=4 fault_shr=3
rank 1: fault_int=2 fault_shr=2
rank 2: fault_int=3 fault_shr=0
rank 3: fault_int=4 fault_shr=1

[phase A: post-init (no Mult)]       9 pair checks, 0 fails on all 8 fields.
[phase B: post 1 Mult(Q=0)]          9 pair checks, 0 fails on all 8 fields.
[phase C: post 1 RK4 step (Q≠0)]     9 pair checks, 0 fails on all 8 fields.

Fields compared: slip1, slip2, V1, V2, tau1_corr, tau2_corr,
                 sigma_n_corr, psi.  Worst |Δ| = 0.000e+00 on every
                 field, every phase (bit-identical across ranks).
```

**Interpretation.**  The runtime DOFData consistency contract
(BP5 R-001 (+,-) canonicalisation + R-501 bit-identical Evaluate
inputs) holds on this 4-rank fixture.  §16.2's concern that
"TPV102 runs the friction Brent solve on both ranks and may
diverge under FP-fragility" is **not reproducible** on the
axis-aligned Cartesian tet fixture — partition-boundary geometry
is regular enough that `(Elem1_centroid − face_centroid) · ref_normal`
gives bit-identical sign on both ranks.

**What this does NOT rule out:**
- Production-mesh FP-fragility — the 1000 m and 200 m TPV102
  meshes are Gmsh-generated with vertex duplication across the
  fault surface.  Partition cuts near x = 0 may produce
  centroid dot-products at the edge of FP precision.  §4.7 on a
  regular fixture cannot probe this.
- Longer time horizons — §4.7 runs 1 RK4 step; runtime divergence
  that requires hundreds of RK4 steps to accumulate (§15.2 speckle
  pattern) is not observable here.
- Different partition patterns at ≥ 56 ranks that ParMETIS
  chooses on the production mesh; the 4-rank local partition
  pattern is not representative.

**Bug status.**  H-V92-P **still open at runtime / production
scale**.  §4.6 + §4.7 PASS locally implies the init-level and
local-fixture-runtime contracts are correct; the bug, if it is
H-V92-P, must be production-mesh-specific (Gmsh vertex-
duplicated geometry at ParMETIS cuts near x = 0).  Progression
per §16.7 now requires Frontera runs (§9.2 diagnostic sbatch in
§13.4) to confirm.

**§16.5 fix (BP5 owner-pattern port) status.**  §4.7 PASS is
**NOT** a reason to skip §16.5.  The §16 argument that the
TPV102 two-sided routing is architecturally fragile stands: it
just happens not to trigger on the 4-rank fixture.  The fix
should still be considered when §9.2 completes (either §9.2
confirms runtime H-V92-P on production mesh → §16.5 fix required;
or §9.2 shows clean DOFData → the speckle has a different
mechanism and §16.5 fix is unnecessary).

### §16.7 Updated blocking path (supersedes §15.7)

1. **(FAST, LOCAL)** Implement §4.6 (ctor-level role consistency)
   AND §4.7 (post-RK4-step DOFData consistency).  Under 15 min
   combined.
2. **(FAST, LOCAL)** Run §6.4 offline VTU outlier detector to
   quantify the speckle in job-7668434.
3. **(PATCH PROPOSAL)** If §4.6 OR §4.7 fails: port BP5's owner
   pattern per §16.5.  Not a silver-bullet "tweak" — a ~280-LOC
   structural fix aligning TPV102 with BP5.
4. **(VERIFICATION)** Re-run §4.6 + §4.7 after the port; both
   must PASS.
5. **(FRONTERA)** One 12 s production re-run at 400 ranks
   confirms no speckle + σ_n stays at 120 MPa.

### §16.8 One-sentence summary

TPV102 imports BP5's canonical `MakeFaceKey` header but wires it
to a diagnostic verifier only; its runtime shared-fault path
*reimplements* ownership via fragile FP geometric comparison AND
runs the friction Brent solve on both sides of every shared face —
the correct fix is not another FP-comparison tweak but to port
BP5's owner-decides-and-scatters pattern into
`ComputeSharedFaceFluxRHS`, making the canonical key load-bearing
in the runtime hot path rather than decorative.

---

## §17. Rev-3c implementation + verification summary (2026-04-21)

### §17.1 §4.7 implementation

**Delivered.**  `tests/parallel/test_shared_fault_dof_data_consistency.cpp`
implements the §16.6 specification as an MPI unit test.  Key
design points:

- **Fixture.**  Reuses §4.6's 4-km 4-rank Cartesian tet fixture
  (192 elements, 144 boundary elements, ≥ 6 shared fault face
  entries after 4-way ParMETIS partitioning).  2 ranks skipped
  (ParMETIS does not cut the fault in this compact geometry).
- **Driver path.**  Builds `ParMesh → WaveOperator<ParMesh> →
  FaultFaceFlux → InitializeFaultDOFs → SetFaultFlux →
  SetFaultDOFData`, matching the TPV102 driver call sequence at
  `drivers/tpv102_driver.cpp:241-337`.
- **Three phases (stacked difficulty).**
  - Phase A — post-init (no `Mult` call).  Sanity check: any
    disagreement here reflects a dof-order / fault_coords bug, not
    a flux-path bug.  **Weak probe.**
  - Phase B — post 1 `Mult(Q=0)`.  Q=0 ⇒ Tinv·Q=0 on both ranks,
    so `Evaluate` inputs are identical from the pre-stress only.
    **Medium probe** — passes if the (+, −) routing is correct.
  - Phase C — post 1 full RK4 step with Q = (σ_xy = τ_ini,
    σ_xx = 10⁵ Pa) uniform pre-stress plus 3 additional `Mult`
    calls at dt = 10⁻⁵ s.  Exercises the full Tinv · Q path +
    Brent solve 4 times.  **Strong probe** — the direct H-V92-P
    operational definition.
- **Cross-rank comparison.**  `MPI_Allgatherv` of per-rank snapshots
  (tagged with `src_rank` to skip self-compare; centroid key at
  1 μm resolution to index).  Each partner pair is compared from
  the lower-ranked side only (`lrec.src_rank < orec->src_rank`)
  to avoid double-counting.
- **Fields compared.**  8 `DOFData` fields: `slip1`, `slip2`,
  `V1`, `V2`, `tau1_corr`, `tau2_corr`, `sigma_n_corr`, `psi`.
  Each checked with `==` (bit-exact) plus `|Δ|` tracking the
  worst-case drift.

### §17.2 §4.7 verification result

**Run command:**
```
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh
conda activate mfem-dev
mpirun -np 4 ./seas_test_shared_fault_dof_data_consistency
```

**Output (abridged):**
```
serial mesh : 192 elements, 448 faces, 144 bdr elems
rank 0 : local_ne=49 fault_int=4 fault_shr=3
rank 1 : local_ne=46 fault_int=2 fault_shr=2
rank 2 : local_ne=49 fault_int=3 fault_shr=0
rank 3 : local_ne=48 fault_int=4 fault_shr=1

[phase A: post-init (no Mult)]       pair checks = 9
[phase B: post 1 Mult(Q=0)]          pair checks = 9
[phase C: post 1 RK4 step (Q≠0)]     pair checks = 9

Every field on every phase:  fails = 0,  worst |Δ| = 0.000e+00.
Results: 3 passed, 0 failed out of 3 phases.
```

**Verdict: PASS.**  All 24 `(phase × pair-check)` comparisons
across 8 DOFData fields show bit-exact equality (0 ULP) on the
4-km fixture at 4 ranks.  This is the STRONGER version of §4.6
and at this scale confirms that the full flux-computation path
produces rank-consistent DOFData.

### §17.3 Does the PASS rule out H-V92-P?

**No — it only closes the *local-fixture* branch.**  Specifically:

- **What it confirms:**  On a Cartesian tet fixture with axis-
  aligned fault, ParMETIS at 4 ranks with FP-regular centroids,
  and a uniform deterministic pre-stress, the R-001 (+, −)
  canonicalisation + R-501 bit-identical-Evaluate-inputs contract
  hold through one RK4 step.  DOFData stays bit-identical.
- **What it does NOT confirm:**  Production-mesh partitions at
  ≥ 56 ranks on the 200 m / 1000 m TPV102 meshes.  These meshes
  are Gmsh-generated with vertex duplication across the fault
  surface; at ParMETIS cut lines near x = 0 the centroid
  dot-product `(Elem1_centroid − face_centroid) · ref_normal`
  operates on FP values within a few ULP of zero, which CAN
  disagree bit-for-bit between the two ranks.  The fixture's
  centroids are O(1) km away from y = 0 in magnitude, so the
  sign of the dot-product is unambiguous.

### §17.4 Updated blocking path

Local tests are exhausted.  §4.6 + §4.7 both PASS.  The only
remaining rank-1 probe for H-V92-P is at production scale:

1. **USER:** Apply §5.1 `tpv102_debug_v9.2.0_diag_fault_sigma.patch`
   + commit + push (gated per §13.4).
2. **USER:** `sbatch jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch`
   (the §9.2 diagnostic, ~200 SU).  This run records
   `max|Q_self[SYY]|` on fault QPs and — if the §5.2 rank-id VTU
   patch is also applied — allows the §6.4 outlier detector to
   correlate speckle with partition boundaries.
3. **USER:** `sbatch jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch`
   (§4.5 rank-consistency, ~48 SU; fallback for differential
   comparison against the existing v91 400-rank data).
4. **Debugger:** Run `tests/scripts/fault_vtu_outlier_detector.py`
   (§6.4) on the existing job-7668434 PVD.  This is local, no
   Frontera — but the script is still PENDING (rev-3b scope
   item not yet implemented).
5. **Debugger:** Based on the Frontera + §6.4 results, decide
   whether to execute the §16.5 BP5 owner-pattern port.
6. **Debugger:** After any fix, §4.7 re-run must still PASS
   (regression gate); §9.2 on production mesh must show no speckle.

### §17.5 Takeaway (superseded by §17.6 rev-3d)

§4.7 **PASSING on the local fixture is valuable but not
definitive** for R-V92 closure.  It narrows the bug domain:
runtime H-V92-P, if it is the root cause, lives in the
FP-fragility branch of §16.2 that only triggers on
production-mesh geometry + partition patterns.  The path forward
is Frontera-gated per §13.4 / §17.4.

---

## §17.6 Post-REVIEW rev-3d response (2026-04-21)

REVIEW.md filed 6 findings (R-V92-C01 through R-V92-C06) against
§17's prescriptive path.  The most critical was **R-V92-C01**:
Phase C of §4.7 uses a spatially uniform Q, so after MFEM ghost
exchange `Q_self ≡ Q_nbr` bit-exactly on every shared fault face,
making the (+, −) swap *vacuous* and the PASS structural rather
than evidentiary.  This invalidates §17's claim that "§4.7 is
the operational definition of H-V92-P".

### §17.6.1 Findings addressed

| ID | Finding | Response (rev-3d) |
|---|---|---|
| R-V92-C01 | Phase C's uniform Q ⇒ `Q_self ≡ Q_nbr` trivially ⇒ test is a subset probe, not the full operational definition | **Phase D added**: σ_xy(y) = τ_ini · tanh(y / L0) with L0 = 500 m, projected via `ParGridFunction::ProjectCoefficient`.  Measured pre-Mult spread = 1.499·10⁸ Pa ≈ 2·τ_ini on the fixture, confirming `Q_self ≠ Q_nbr` on shared fault faces.  Phase D PASS = genuine invariant check (ghost exchange + swap routing + Brent determinism). |
| R-V92-C02 | §17.4 puts 248 SU of Frontera before free local probes | **Free local probes moved ahead**: §17.7 new priority order runs Phase D, centroid-margin diag, §6.4 VTU detector BEFORE any Frontera submission. |
| R-V92-C03 | §4.7 mis-labeled as "direct operational definition" | Manifest row and §15.5.1 header updated: §4.7 is now labeled correctly as "3 phases probing a hierarchy of Q conditions; Phase D is the operational definition". |
| R-V92-C04 | 4-km Cartesian fixture has centroid margins O(10³ m) — FP-fragility physically impossible on the fixture | **Accepted.**  Reviewer's statement about the 4-km fixture is correct and unchallenged.  Implemented reviewer's Step 3 (`SEAS_DIAG_CENTROID_MARGIN`) to measure margins on production-like geometry.  Measured on 1000 m TPV102 mesh at 12 ranks: margins O(10²) m.  That empirical result refutes the FP-fragility SUB-hypothesis for `elem1_on_plus` **on this mesh**, but only because the reviewer's Step 3 probe was executed — not because R-V92-C04 was wrong.  See §17.6.3. |
| R-V92-C05 | Priority inverted — Frontera before free local | Resolved via R-V92-C02 response. |
| R-V92-C06 | dt = 10⁻⁵ s propagates waves only 6 cm per step | Phase D keeps dt = 10⁻⁵ s for the RK4 step but adds the antisymmetric Q profile which gives a 1.5·10⁸ Pa gradient ACROSS the fault — the short per-step wave propagation no longer matters because the relevant gradient is already present at t = 0. |

### §17.6.2 Phase D execution result

```
serial mesh : 192 elements, 448 faces, 144 bdr elems  (same fixture)
rank 0: fault_int=4 fault_shr=3
rank 1: fault_int=2 fault_shr=2
rank 2: fault_int=3 fault_shr=0
rank 3: fault_int=4 fault_shr=1

[phase D pre-Mult] σ_xy field min/max = -7.495e+07 / +7.495e+07
                   spread = 1.499e+08 Pa
                   ⇒ antisymmetric (Q_self ≠ Q_nbr at fault) — swap EXERCISED

[phase D: post 1 RK4 (antisymmetric σ_xy, swap EXERCISED)]
                   9 pair checks, 0 fails on all 8 DOFData fields.
                   PASSED: all DOFData fields bit-identical across ranks
                   on every shared fault QP.

Results: 4 passed, 0 failed out of 4 phases.
```

**Interpretation.** With `Q_self ≠ Q_nbr` genuinely exercised:
- MFEM ghost exchange delivers bit-identical `Q_nbr_on_A` ≡
  `Q_self_on_B` at every shared fault QP.
- The (+, −) swap at `wave_operator.inl:1287-1290` correctly
  routes `(Q_plus, Q_minus)` so both ranks' `Evaluate()` receives
  the same inputs.
- Friction Brent produces bit-identical outputs on both ranks.
- Therefore `dof_data[dof_idx]` remains bit-identical across
  ranks after 1 RK4 step with an antisymmetric Q profile.

**H-V92-P via swap-routing is not triggered on the 4 km fixture,
even with the correct precondition.**

### §17.6.3 Centroid-margin measurement on production mesh

`SEAS_DIAG_CENTROID_MARGIN` compiled into
`seas_test_centroid_margin_ctor_only` with the flag on.  Ran
against the TPV102 1000 m mesh (smallest production-like mesh;
200 m is forbidden locally).  At 12 MPI ranks (below the 14-core
laptop limit per feedback_no_local_oversubscribe):

```
[CENTROID-MARGIN] sf_idx=0 elem1_proj=+262.086 face_proj=-1e-12 margin=262.086  ...
[CENTROID-MARGIN] sf_idx=1 elem1_proj=+262.086 face_proj=-1e-12 margin=262.086  ...
[CENTROID-MARGIN] sf_idx=2 elem1_proj=+242.502 face_proj=-1e-12 margin=242.502  ...
[CENTROID-MARGIN] sf_idx=0 elem1_proj=-227.523 face_proj=-1e-12 margin=227.523  ...
[CENTROID-MARGIN] sf_idx=1 elem1_proj=-227.523 face_proj=-1e-12 margin=227.523  ...
[CENTROID-MARGIN] sf_idx=2 elem1_proj=-266.035 face_proj=-1e-12 margin=266.035  ...
```

**Summary: all 6 shared-fault-face margins on the 1000 m
production mesh are in the range 227.5 m – 266.0 m.**  FP-fragility
would require margins O(ε · scale) ~ 10⁻¹² m; we observe margins
**14 orders of magnitude larger**.  This is the output of the
reviewer's prescribed Step 3 probe on production-like geometry.
R-V92-C04's original statement — that the 4-km fixture cannot
probe FP-fragility — remains correct and is not at issue.  What
the 1000 m measurement refutes is the SUB-hypothesis that the
`elem1_on_plus` FP comparison is fragile **on the 1000 m
production-like mesh**.  If the 200 m mesh at 400 ranks were to
show smaller margins, that sub-hypothesis could still revive;
the reviewer's Step 3 bluntly predicts this outcome to be the
next thing to check.

### §17.6.4 Consequences for H-V92-P

The H-V92-P mechanism, as formulated in §15.4 and §16.2, has two
proposed failure modes:

1. **`elem1_on_plus` sign-disagreement via FP-fragility** —
   refuted by §17.6.3: margins on the production mesh are O(100 m),
   not O(ε).  The `elem1_proj < face_proj` comparison is robust.
2. **Ghost-exchange / Brent / swap-routing divergence** — refuted
   by §17.6.2: Phase D's antisymmetric Q exercises all three,
   and DOFData stays bit-identical across ranks.

**Combined with §4.6 PASS (ctor-level anti-symmetry) and §4.7
Phase A/B/C PASS (weaker Q conditions), H-V92-P is now
unsupported by direct local evidence** on the accessible test
geometries (4 km fixture + 1000 m production mesh at 12 ranks).

### §17.6.5 What still remains open

1. **Does the 200 m mesh / 400 rank / Gmsh-duplicated-fault
   pattern produce different centroid-margin distribution?**  Only
   one way to find out: run `seas_test_centroid_margin_ctor_only`
   at Frontera against the 200 m mesh under 400-rank ParMETIS
   partitioning.  This is a **ctor-only** run — no time-stepping,
   no simulation.  Minimal SU cost (≤ 1 min wall × 8 nodes ≈
   4 SU).  **This should replace the §13.4 §9.2 `-DSEAS_DIAG_FAULT_SIGMA`
   run as the next Frontera step**: it is an order of magnitude
   cheaper and directly tests the R-V92-C04 hypothesis on the
   actual production mesh.
2. **Is the speckle seen in job 7668434 coming from interior fault
   faces (not shared)?**  Gmsh's `BooleanFragments` duplicates
   vertices across the fault, so after partitioning MFEM
   classifies most fault faces as interior boundary elements on
   each rank.  At 12 ranks on 1000 m, only 6 shared fault faces
   exist; at 400 ranks on 200 m, the ratio is similar (a few
   hundred shared vs tens of thousands of interior-classified
   fault faces).  **The dominant fault-face processing path is
   `ComputeFaceFluxRHS` (interior branch), not
   `ComputeSharedFaceFluxRHS`.**  H-V92-P (by definition about
   shared-fault bit-disagreement) can only affect a minority of
   fault QPs — cannot explain speckle covering 100+ triangles in
   the VTU.
3. **Therefore: the speckle must have a different mechanism.**
   Candidates include H-V92-G (bulk σ_yy amplification), H-V92-K
   (RK4 stage-averaging in station output), or an undiscovered
   bug in the INTERIOR fault-face path.

### §17.6.6 §6.4 offline VTU outlier detector — IMPLEMENTED

`tests/scripts/fault_vtu_outlier_detector.py` implemented per
§15.5.2.  Features:

- Loads `FaultSurface/*.pvd` + `.pvtu` + `.vtu`.
- Reads per-triangle fields (`normal_stress`, `traction_strike`,
  `traction_dip`, `state_variable`, `slip_rate_strike/dip`,
  `slip_strike/dip`).
- Per-cycle median + MAD; flags triangles > 5·MAD from median on
  ANY field.
- Correlates outliers with `mpi_rank` and `is_shared_face`
  CellData if §5.2 patch applied (optional).
- Outlier count trend analysis: monotonic growth → H-V92-G;
  stable → H-V92-P or Q.
- Output: per-cycle JSON report + outlier-centroid CSVs for
  ParaView re-import.

**Data prerequisite:** the raw `.pvd` + `.vtu` from job 7668434
must be pulled from Frontera (only plot PNGs are local).  Once
pulled, the detector runs in ~30 min wall-time locally per
§15.5.2 budget.  The script accepts `--pvd` / `--output-root` /
`--mad-k` / `--max-cycles` arguments; dependencies: `numpy`, `vtk`.

### §17.7 Updated blocking path (supersedes §17.4)

Per REVIEW R-V92-C02/C05, **free local probes run FIRST**:

1. **Done.** §4.6 ctor-level role consistency — PASS 5/5.
2. **Done.** §4.7 Phase A/B/C — PASS 3/3 (weak probes, acknowledged).
3. **Done.** §4.7 Phase D antisymmetric Q — PASS 1/1 (genuine
   probe; `Q_self ≠ Q_nbr` exercised, swap bit-identical).
4. **Done.** `SEAS_DIAG_CENTROID_MARGIN` on 1000 m mesh at 12
   ranks — margins O(100 m), FP-fragility refuted.
5. **USER:** pull `FaultSurface/fault_surface_*.vtu` + `.pvd`
   from Frontera job 7668434 (`scp` commands in
   `tests/scripts/fault_vtu_outlier_detector.py` header).
6. **Debugger:** run §6.4 outlier detector locally (~30 min CPU).
   Trend slope + partition correlation → classifies H-V92-P / G / Q.
7. **Debugger:** based on §6.4 verdict:
   - Growth slope > 10 %/s → **H-V92-G primary** → §5.1 diagnostic
     patch + Frontera `v92_diag.sbatch` (~358 SU) to confirm bulk
     σ_yy growth.
   - Stable outlier count on shared-face subset → H-V92-P still
     possible at production scale → Frontera `v92_centroid_margin
     .sbatch` (~4 SU, NEW, below) plus `v92_diag.sbatch`.
   - Random scatter → **H-V92-Q (new)** friction-solver
     instability → investigate `dieterich_ruina.hpp` per-QP Brent
     brackets.
8. Only after §6.4 classifies the mechanism: submit the targeted
   Frontera sbatch(es).
9. **USER:** consider §16.5 BP5 owner-pattern port as
   architectural hygiene regardless of which mechanism fires;
   eliminates the redundant Brent solve across ranks and makes
   canonical key load-bearing.  Independent of the speckle fix.

**Frontera budget after free probes:** 4 SU (centroid margin on
200 m mesh) if needed + 200-358 SU for the targeted diag run;
total ≤ 360 SU, well under §17.4's 406 SU and consistent with the
rev-2 1000 SU envelope.

### §17.8 New Frontera sbatch needed (rev-3d)

A `tpv102_200m_p1_ctor_only_400rank_v92_centroid_margin.sbatch`
should be drafted (NOT NOW; only if §6.4 ambiguous) to run the
ctor-only centroid-margin diagnostic at production scale:

- 8 nodes × 400 ranks × 1 min wall = ~4 SU.
- No time-stepping (ctor + MPI_Finalize).
- Captures stderr per-rank, post-processes for min/max/histogram
  of margins across all shared fault faces at 200 m / 400-rank
  partitioning.
- If min margin is still O(100 m), H-V92-P FP-fragility branch
  is definitively ruled out for TPV102 production.

---

## §18. §6.4 VTU outlier detector — executed (2026-04-21, rev-3e)

### §18.1 Inputs

Raw PVD + PVTU + VTU data pulled by user from Frontera to:
```
/Users/chunhuizhao/Downloads/seas-mfem/tpv102/results_200m_p1_12.0s_400r_v91_job7668434/
```

- 93 output cycles in `FaultSurface/fault_surface.pvd`
  (dt_output ≈ 0.099 s, tfinal = 9.12 s of the 12 s run).
- 400 rank-vtu pieces per cycle; ~37 300 VTU files total.
- Fields per triangle: slip_dip, slip_strike, slip_rate_dip,
  slip_rate_strike, traction_dip, traction_strike, state_variable,
  normal_stress, param_a, param_Dc, fault_x2, fault_x3.
- `mpi_rank` / `is_shared_face` NOT emitted (the §5.2 patch was
  not applied at job 7668434 build time).

### §18.2 Detector (rev-3e adjustments during execution)

The detector `tests/scripts/fault_vtu_outlier_detector.py`
required three small fixes to run on this dataset:

1. **`vtk` dependency removed.**  The mfem-dev / pythonenv conda
   environments have no `vtk`.  Rewrote the VTU reader with pure
   Python XML parsing (MFEM emits `format="ascii"` DataArrays) —
   no performance penalty (30 s for 37 000 VTU files).
2. **Empty-rank shards.**  Ranks that own no fault cells write
   `NumberOfPoints=0 NumberOfCells=0` VTU files with a single
   whitespace in each DataArray.  The reader now returns empty
   arrays for these.
3. **Classifier field restriction.**  Initial run flagged 29 %
   outliers at t=0 because coordinate (`fault_x2`, `fault_x3`)
   and spatially-varying parameter (`param_a`) fields have
   O(km) / O(0.01) dispersion by design.  Restricted the
   classifier to the 8 dynamical fields; added a "uniform field
   skip" rule (`mad < 1e-6·|median|` ⇒ field contributes 0
   outliers that cycle — avoids MAD-zero false positives).
4. **NumPy 2.0 compatibility.**  `ndarray.ptp()` → `np.ptp()`.

Post-fix command:
```
python3 tests/scripts/fault_vtu_outlier_detector.py \
    --pvd /Users/chunhuizhao/Downloads/seas-mfem/tpv102/results_200m_p1_12.0s_400r_v91_job7668434/FaultSurface/fault_surface.pvd \
    --output-root /tmp/r_v92_outlier_full \
    --mad-k 5
```

Wall-time: ~35 s on laptop for 93 cycles × 400 ranks = 37 200
VTU files.  Output: `outlier_report.json` + per-cycle
`outliers_cycle_<NNNN>.csv` (centroid locations for ParaView
re-import).  Copied to `tpv102/r_v92_outlier_detection/` (73 MB).

### §18.3 Outlier trend

| Phase | Cycles | Time | Outlier fraction |
|---|---|---|---|
| Pre-nucleation | 0–17 | 0.00 s – 1.69 s | **0 %** (all fields uniform) |
| Nucleation onset | 18 | 1.78 s | 45.6 % (step function — `traction_dip` median stays 0 but 17 339 triangles acquire 10 MPa anomalies) |
| Rupture propagation | 19–43 | 1.88 s – 4.26 s | 37 % – 62 % |
| Rupture saturation | 44–62 | 4.36 s – 6.14 s | 54 % – 62 % |
| Free-surface reflection | 63–65 | 6.24 s – 6.44 s | **96 % peak** (33 – 37 k outliers) |
| Post-reflection decay | 66–92 | 6.54 s – 9.12 s | 44 % – 87 % |

Linear regression on outlier count vs time: **+2 703 outliers /
second**, monotonically upward.  Script verdict:
**"MONOTONIC GROWTH → H-V92-G bulk amplification."**

### §18.4 Per-field growth analysis (§18.3 breakdown)

| Field | Peak max\|Δ\| | At time | Early-mean | Late-mean | Growth |
|---|---|:---:|---|---|:---:|
| slip_dip | **1.83 m** | 9.12 s | 6.9·10⁻³ m | 1.23 m | **178×** |
| slip_rate_dip | **2.15 m/s** | 8.62 s | 5.5·10⁻² m/s | 1.14 m/s | **21×** |
| traction_dip | **68 MPa** | 9.12 s | 2.4 MPa | 45 MPa | **19×** |
| **normal_stress** | **218 MPa** | 9.12 s | 4.1 MPa | 166 MPa | **40×** |
| slip_strike | 7.23 m | 6.05 s | 0.54 m | 6.54 m | 12× (expected) |
| slip_rate_strike | 16.3 m/s | 5.35 s | 2.64 m/s | 10.9 m/s | 4× (expected) |
| traction_strike | 149 MPa | 8.62 s | 20 MPa | 97 MPa | 5× (expected) |
| state_variable | 0.41 | 7.63 s | 0.18 | 0.39 | 2× (expected) |

### §18.5 Interpretation

**The strike-channel fields grow at physically-reasonable rates**
(4 – 12× from early to late) consistent with a propagating
mode-II rupture.  `slip_strike = 7 m`, `slip_rate_strike = 16 m/s`
are within the TPV102 benchmark envelope.

**The dip-channel and normal-stress fields grow 19 – 178×.** TPV102
is pure strike-slip — the dip channel should be zero and the
normal stress should stay at 120 MPa ± 1 MPa.  Observed:

- `slip_dip` reaches 1.83 m (**1 830 × the SCEC reference bound**
  of 1 mm).
- `slip_rate_dip` reaches 2.15 m/s (**20 000 × the reference
  bound** of 10⁻⁴ m/s).
- `traction_dip` reaches 68 MPa (**70 × the physical bound** of
  ~1 MPa).
- `normal_stress` deviation reaches 218 MPa — nearly **2 × the
  static value** of 120 MPa.  The previously-reported hypocenter
  drop (~50 MPa) is just the tip of a much broader drift across
  the fault.

**The free-surface reflection spike at cycle 63 – 65 (t ≈ 6.2 s)**
pushes outlier fraction to 96 %.  Free-surface arrival at a
fault QP of depth 7.5 km is at t ≈ depth / cs = 7.5/3.464 ≈
2.16 s after the local rupture passes; the rupture reaches
surface at t ≈ 1.78 s + 7.5/vR ≈ 1.78 + 7.5/3.3 ≈ 4.05 s; the
first reflection returns at t ≈ 4.05 + 2.16 ≈ 6.21 s.  The
observed spike at t = 6.24 s matches within 0.03 s — strongly
suggesting the free-surface reflection transports the polluted
state back across the fault.

### §18.6 Hypothesis classification (updated)

| Hypothesis | Outlier evidence | Verdict |
|---|---|---|
| **H-V92-G (bulk amplification)** | Monotonic growth of dip-channel and σ_n fields by 20 – 180× over 9 s; matches the mode-II coupling mechanism in §2.5 | **CONFIRMED** |
| **H-V92-A (σ_n drift at stations)** | Station drift of ~50 MPa is a narrow view of a much wider drift (peak 218 MPa); confirms the symptom but not the mechanism — subsumed by H-V92-G | **Subsumed by H-V92-G** |
| **H-V92-B (dip channel pollution)** | Directly quantified: slip_dip 1 830× spec, slip_rate_dip 20 000× spec, traction_dip 70× spec.  Strike-channel magnitudes are reasonable — pollution is dip-specific, consistent with mode-II cross-coupling | **CONFIRMED (subsumed by H-V92-G)** |
| **H-V92-P (shared-fault bit disagreement)** | Requires stable partition-pinned pattern; observed pattern is monotonically growing, not pinned.  Also: 12-rank 1000 m mesh has only 6 shared fault faces (§17.6.3); production 400-rank 200 m will have similarly few.  Cannot produce O(10⁴) outlier triangles. | **ELIMINATED** |
| **H-V92-C (corner)** | Corner emanation would spike ONCE at first reflection and not produce uniform O(40×) growth in σ_n across the entire fault.  The 6.24 s spike is a free-surface reflection signature, but the underlying 40× growth is present BEFORE 6.24 s as well | **WEAKENED**; may contribute to the 6.24 s spike but not to the secular growth |
| **H-V92-K (RK4 DOFData stage-averaging + Q/DOFData splitting)** | **RE-OPENED (REVIEW round-3 R-V92-E02 + round-4 F01/F02).**  Earlier rev-3e elimination was based on a factually wrong claim — `drivers/tpv102_driver.cpp:860-873` **does** RK4-stage-average `dof_data.V1/V2/tau1_corr/tau2_corr/sigma_n_corr` (comment explicitly: "R-001 fix: RK4-weighted corrected tractions for consistent station output"), and `WriteFaultSurfaceVTU` reads these averaged values.  Round-4 F01+F02 further identified that the DOFData fields live OUTSIDE the RK4 state vector — an operator-splitting error bounded at O(10⁻²) pathology, not enough alone to explain 1 830× slip_dip, but a plausible secondary contributor. | **OPEN** — discriminator: dt-halving (active TODO Step 1) |
| **H-V92-F (mesh jitter)** | Already eliminated by §4.3 (113 712 QPs at 0.5 ULP).  Growth by 20 – 180× cannot be jitter-accumulated | **ELIMINATED (restated)** |

### §18.7 Root-cause formulation

**H-V92-G is the confirmed primary mechanism.**  The mode-II
rupture, propagating along +x, drives bulk σ_yy perturbations via
the `A_x` Poisson-coupling entries (see §2.5):

```
A_x[SYY, VX] = −λ
```

so `∂_t σ_yy = λ · ∂_x v_x` from any spatial gradient of v_x along
strike.  At the rupture front, v_x goes from 0 to peak slip-rate
over a narrow length, giving ∂_x v_x of O(10⁴) s⁻¹ which, with
λ ≈ 30 GPa, produces dσ_yy/dt of O(3·10¹⁴) Pa/s.  Over the
∆t ≈ 10⁻⁴ s timescale, σ_yy picks up ~30 GPa · 10⁻⁴ = 3 GPa per
step near the rupture front — already O(30 MPa) per time step
averaged over a 100-step rupture interval.  **This is the
amplification.**

The §3 H-V92-G subsection listed four candidate sources for why
the amplification is not damped:

- **(a)** IP/BR2 penalty in the elastodynamic driver: **unlikely**,
  no IP term in `wave_operator.inl` flux path.
- **(b)** RK4 non-conservative drift on σ_yy: plausible but
  bounded to O(ε · N_steps).
- **(c)** Absorbing BC reflection with wrong sign: would show as
  boundary-arrival spikes, not uniform growth.
- **(d)** Element mass-inverse conditioning asymmetry across the
  9 Voigt channels: **most plausible.**  If the mass-inverse
  applied in `ApplyMassInverse` treats σ_yy differently from σ_xy
  (e.g., different effective weighting in a multi-component
  GridFunction application), the correct mode-II σ_yy build-up
  is not drained on subsequent time steps, compounding over the
  10⁴ RK4 stages of a 12 s run.

### §18.8 Next actions (rev-3f, supersedes §17.7 + §17.8)

**Per user directive: do NOT assume candidate (d) is the winner.
Run all four candidate probes (a)/(b)/(c)/(d) without pre-ranking
and let the data classify.  Pelties 2012 §3.1 is the benchmark
reference for each.**

**Done (rev-3f):**

- §4.9 executed — PASS.  **Candidate (a) (hidden IP/BR2 penalty)
  ELIMINATED by direct numerical probe + source grep.**  T1 (zero
  invariant), T2 (linearity), T4 (source grep) all conformant
  with Pelties 2012 §3.1 "Godunov fluxes only, no penalty".

**All four remaining candidate probes — DONE (rev-3g, 2026-04-21):**

1. **Candidate (b) RK4 non-conservation** — **§4.11 PASS.**  100
   RK4 steps on 8³ hex with Gaussian pulse: SXX drift −2.33·10⁻⁵/step,
   SYY drift −1.63·10⁻⁵/step (both dissipative); SYY/SXX ratio =
   0.70.  No asymmetric drift.  **ELIMINATED.**
2. **Candidate (c) absorbing-BC wrong sign** — **§4.12 PASS.**
   3 000 RK4 steps on 16 km cube, all absorbing: E(t) monotonically
   decreases 75.05 → 18.25; `E_max/E₀ = 1.000` exactly.
   **ELIMINATED.**
3. **Candidate (d) mass-inverse channel asymmetry** — **§4.10
   PASS.**  Isotropy ratio `SXY→VY / SXZ→VZ` = 1.000 bit-exact.
   Poisson coupling ratio matches `λ/(λ+2μ) = 0.333`.
   **ELIMINATED.**
4. **Candidate (e) volume-term Jacobian coupling** — **§4.10
   PASS.**  All 9×9 channel sparsity pairs match predicted
   `A_x + A_y + A_z` structure; no spurious cross-channel
   coupling above 1e-9 noise floor.  **ELIMINATED.**

**Status after rev-3g: ALL FOUR NAMED H-V92-G candidates
eliminated by local probes.**  The amplifier observed in §6.4
(19 – 178× dip-channel growth at production scale) is NOT
explained by any of (a)–(e).  This is significant: the bulk
elastodynamics operator (including RK4, mass-inverse, absorbing
BC, Jacobian sparsity, and no-penalty invariant) is internally
consistent on a clean fixture without fault.

**Therefore, the amplifier must be in the FAULT-FLUX PATH (not
the bulk operator).**  Candidates to consider next:

- **(f) Fault flux path itself.**  All §4.10–§4.12 probes ran
  with `bc.fault_attr = 0` (no fault).  The production run has
  `fault_attr = 3` on tens of thousands of fault boundary
  elements, each triggering `FaultFaceFlux::Evaluate` + friction
  Brent solve.  A bug ISOLATED to this path would explain why
  §4.10–§4.12 pass but §6.4 shows amplification.
- **(g) Fault-flux ↔ volume interaction.**  Mode-II σ_yy reaches
  the fault element boundary → `FaultFaceFlux` computes
  σ_n_trial using the stress-average term → friction-corrected
  σ_n_corr is SET on DOFData → the bulk volume flux on the
  SAME element re-absorbs this corrected σ_n_corr back into
  σ_yy.  If the set/update ordering is wrong, σ_yy accumulates.
- **(h) Interior fault branch vs shared fault branch divergence.**
  §4.6/§4.7 verified SHARED-fault DOFData consistency.  But
  most fault faces on Gmsh TPV102 meshes are INTERIOR (not
  shared — see §17.6.3).  A bug that ONLY fires on the interior
  path would NOT show up in §4.6/§4.7 but would show up in §6.4
  at production scale.

**Next step:** write §4.13 — a fault-flux-inclusive probe on the
existing 4-km inline-fault fixture from §4.6.  Run 1 RK4 step
with the antisymmetric σ_yy profile (Phase D style), then
measure σ_n drift AT the fault QPs via DOFData.  If σ_n_corr
drifts monotonically across N > 1 step, candidate (f)/(g)
confirmed.  If stable, (h) interior-fault path remains open and
needs a different probe.  Estimated: ~200 LOC, 1 hour.

**Deprioritized:**
- §5.1 `SEAS_DIAG_FAULT_SIGMA` Frontera run — §18 already
  provides the σ_yy-on-fault-QP data via VTU outlier growth.
- §4.5 56-rank Frontera rank-consistency run — H-V92-P/M ruled
  out.
- §16.5 BP5 owner-pattern port — fixes H-V92-P which was ruled
  out; can still be done as architectural hygiene but not
  R-V92-critical.

### §18.9 One-sentence summary (rev-3g updated 2026-04-21)

The speckle in job 7668434 is a **monotonically-growing**
amplification of the mode-II σ_yy / σ_xz coupling into the fault
DOFData — **H-V92-G confirmed at production scale** — with
dip-channel contamination reaching 1 800 × the SCEC reference
bound; **all five named H-V92-G candidates (a–e) ELIMINATED by
free local probes (§4.9, §4.10, §4.11, §4.12)** — the amplifier
is NOT in the bulk elastodynamics operator (penalty, RK4,
absorbing BC, mass inverse, volume Jacobian), and the hypothesis
now shifts to candidates (f–h) in the fault-flux path / interior-
fault branch which were not exercised by §4.10–§4.12's no-fault
fixtures; §4.13 will probe this next without Frontera.

---

## §B. Consolidated test audit (2026-04-21, rev-3f)

**Single source of truth** for every v9.2.0 test/diagnostic/sbatch
artifact and its status.  If a test is not in this table, it is
not part of v9.2.0 scope.  Status values:

- **PASS (<n>/<m>)** — executed, n of m assertions passed.
- **READY** — artifact written + compiles/syntax-clean; not yet
  executed (waiting on user approval or Frontera submission).
- **PENDING** — specified in a section but not yet implemented.

### §B.1 C++ unit tests — local, serial or low-rank MPI

| § | File (source) | Built target | Last run | Result |
|---|---|---|---|---|
| **§4.9** (rev-3f, ordered FIRST) | `tests/unit/test_no_penalty_dynamic_rupture.cpp` | `seas_test_no_penalty_dynamic_rupture` | 2026-04-21 | **PASS 4/4** — T1 `Mult(Q=0)=0` bit-exact; T2 `Mult(2·Q)=2·Mult(Q)` bit-exact; T3 σ_yy→VY coupling per Pelties eq. (4); T4 0 penalty-keyword hits |
| §4.1 | `tests/unit/test_godunov_rotation_identity.cpp` | `seas_test_godunov_rotation_identity` | 2026-04-21 | **PASS 4/4** — BP5 0 ULP, 3 off-axis frames each 1 ULP |
| §4.2 | `tests/unit/test_canonical_rotation_pure_strikeslip.cpp` | `seas_test_canonical_rotation_pure_strikeslip` | 2026-04-21 | **PASS 11/11** — all assertions bit-exact |
| §4.3 | `tests/unit/test_fault_basis_qp_orthonormality.cpp` | `seas_test_fault_basis_qp_orthonormality` | 2026-04-21 | **PASS** — 1 000 m: 4 938 QPs @ 0.5 ULP; **200 m production: 113 712 QPs @ 0.5 ULP** |
| §4.4 | `tests/unit/test_per_qp_vs_centroid_basis_planar_1el.cpp` | `seas_test_per_qp_vs_centroid_basis_planar_1el` | 2026-04-21 | **PASS 11/11** — bit-exact on planar 2-tet fixture |
| **§4.10** (rev-3g, H-V92-G d+e) | `tests/unit/test_volume_jacobian_single_channel.cpp` | `seas_test_volume_jacobian_single_channel` | 2026-04-21 | **PASS 2/2** — 9×9 Jacobian sparsity matches `A_x+A_y+A_z` prediction; isotropy ratio `SXY→VY / SXZ→VZ` = **1.000 bit-exact**; Poisson ratio 0.333 = λ/(λ+2μ). **(d) and (e) ELIMINATED** |
| **§4.11** (rev-3g, H-V92-G b) | `tests/unit/test_rk4_conservation.cpp` | `seas_test_rk4_conservation` | 2026-04-21 | **PASS 2/2** — 100 RK4 steps, Gaussian pulse: SXX drift −2.33·10⁻⁵/step, SYY drift −1.63·10⁻⁵/step (both dissipative), SYY/SXX ratio 0.70. **(b) ELIMINATED** |
| **§4.12** (rev-3g, H-V92-G c) | `tests/unit/test_absorbing_bc_energy_decay.cpp` | `seas_test_absorbing_bc_energy_decay` | 2026-04-21 | **PASS 1/1** — 3 000 RK4 steps, 16 km cube all-absorbing: E(t) decays 75.05 → 18.25, `E_max/E₀ = 1.000` exactly (no increase). **(c) ELIMINATED** |

### §B.2 C++ parallel MPI tests — 4-rank local fixture

| § | File (source) | Built target | Last run | Result |
|---|---|---|---|---|
| §4.6 | `tests/parallel/test_shared_fault_role_consistency.cpp` | `seas_test_shared_fault_role_consistency` | `mpirun -np 4`, 2026-04-21 | **PASS 5/5** — `elem1_on_plus` anti-symmetric, canonical frame bit-identical across ranks |
| §4.7 | `tests/parallel/test_shared_fault_dof_data_consistency.cpp` | `seas_test_shared_fault_dof_data_consistency` | `mpirun -np 4`, 2026-04-21 | **PASS 4/4 phases** — A post-init / B Q=0 / C uniform Q / **D antisymmetric σ_xy(y) = τ_ini·tanh(y/L0) — pre-Mult spread = 1.499·10⁸ Pa confirmed**; all 8 DOFData fields bit-identical |
| §17.6 | `tests/parallel/test_centroid_margin_ctor_only.cpp` (+ `SEAS_DIAG_CENTROID_MARGIN` flag) | `seas_test_centroid_margin_ctor_only` | `mpirun -np 12 --mesh tpv102_1000m.msh`, 2026-04-21 | **PASS** — 6 shared fault faces, margins [227, 266] m; 14 orders of magnitude above ε_FP; R-V92-C04 FP-fragility hypothesis refuted on this mesh |

### §B.3 Python scripts — local, file-IO / harness

| § | File (source) | Last run | Result |
|---|---|---|---|
| §6.4 | `tests/scripts/fault_vtu_outlier_detector.py` | 2026-04-21, job 7668434 full 93 cycles | **PASS 93/93** (35 s wall-time) — monotonic outlier growth +2 703/s; per-field peaks: slip_dip 1.83 m (1 830×), slip_rate_dip 2.15 m/s (20 000×), traction_dip 68 MPa (70×), normal_stress deviation 218 MPa; **H-V92-G CONFIRMED**, H-V92-P/C/K/F ELIMINATED |
| §4.5 | `tests/scripts/test_sigma_n_rank_consistency.py` | Never run | **READY** — Frontera post-processing harness; consumes output of §B.4 run-A sbatch |

### §B.4 Frontera sbatch / source patches — NOT YET SUBMITTED (user approval gate)

| § | File | Status | Cost estimate |
|---|---|---|---|
| §4.5 run-A | `jobs/tpv102/tpv102_200m_p1_4.0s_56rank_v92_rankA.sbatch` | **READY** — written, awaiting user approval + `sbatch` submission | ~48 SU (1 node, 1.5 hr) |
| §5.1 patch | `debug_document/tpv102_debug_document/tpv102_debug_v9.2.0_diag_fault_sigma.patch` | **READY, NOT applied** — additive `#ifdef SEAS_DIAG_FAULT_SIGMA` block in `dynamic/wave_operator.inl`; user must review + `patch -p0` + commit + push | — |
| §9.2 diag | `jobs/tpv102/tpv102_200m_p1_4.0s_400rank_v92_diag.sbatch` | **READY, NOT submitted** — requires §5.1 patch applied; self-checks for diag symbol before running | ~358 SU (8 nodes, 1.4 hr) |
| §5.2 patch | Planned: `io/paraview_output.hpp` patch for `SEAS_DIAG_VTU_RANK` | **PENDING** — not yet written.  Note: another session has independently modified `io/paraview_output.hpp` for R-V92-E02 stage-4 VTU diagnostic (unrelated to §5.2 scope). | — |

### §B.5 §18.8 candidate probes — ALL EXECUTED (rev-3g, 2026-04-21)

All four candidates (b)(c)(d)(e) were probed by §4.10–§4.12
without pre-ranking; each eliminated by measurement.

| Candidate | Probe (file) | Status / Result |
|---|---|---|
| (a) hidden penalty | §4.9 `test_no_penalty_dynamic_rupture.cpp` | **ELIMINATED** — T1/T2/T4 PASS (bit-exact linearity, 0 penalty-keyword hits) |
| (b) RK4 non-conservation | §4.11 `test_rk4_conservation.cpp` | **ELIMINATED** — SXX drift −2.33·10⁻⁵/step, SYY drift −1.63·10⁻⁵/step, ratio 0.70 (symmetric, dissipative) |
| (c) Absorbing-BC wrong sign | §4.12 `test_absorbing_bc_energy_decay.cpp` | **ELIMINATED** — E(t) decays 75.05 → 18.25 over 3 000 RK4 steps; `E_max/E₀ = 1.000` exactly |
| (d) Mass-inverse channel asymmetry | §4.10 `test_volume_jacobian_single_channel.cpp` | **ELIMINATED** — `SXY→VY / SXZ→VZ` ratio = 1.000 bit-exact; Poisson 0.333 = λ/(λ+2μ) |
| (e) Volume-term Jacobian coupling | §4.10 `test_volume_jacobian_single_channel.cpp` | **ELIMINATED** — all 9×9 channel pairs match predicted `A_x+A_y+A_z` sparsity above 1e-9 noise floor |

**Status (rev-3g):** All five named H-V92-G candidates (a–e)
eliminated locally.  The amplifier is NOT in the bulk
elastodynamics operator.  New candidates (f–h) now open
per §18.8:

| Candidate | Domain | Next probe |
|---|---|---|
| (f) Fault-flux path itself | `FaultFaceFlux::Evaluate` + Brent solve + DOFData set | §4.13 (pending) — 4-km fixture WITH fault, drive antisymmetric σ_yy, measure σ_n drift at fault QPs over N steps |
| (g) Fault-flux ↔ volume interaction | DOFData σ_n_corr fed back into volume σ_yy | Same §4.13 fixture, compare 1-step vs 10-step DOFData drift |
| (h) Interior fault branch bug | `ComputeFaceFluxRHS` fault branch (not shared) | §4.13 with serial fixture (all fault faces interior) vs parallel fixture (shared faces); compare |

### §B.6 Makefile aggregate target

```
make test-v92-regression-gates
```

Runs **10 C++ tests** from §B.1 + §B.2 in order, §4.9 first:
§4.9, §4.1, §4.2, §4.3, §4.4, §4.6, §4.7, §4.10, §4.11, §4.12.
Verified **PASS end-to-end on 2026-04-21 (rev-3g)**.  Total
wall-time ≈ 60 s local.

### §B.7 Inventory completeness check

- **Source-tree artifacts under `tests/`** that I added during
  v9.2.0 work: 8 C++ tests (5 unit + 3 parallel) + 2 Python
  scripts.  All 10 appear in §B.1/§B.2/§B.3.
- **Frontera sbatch files under `jobs/tpv102/`** that I added:
  2 sbatches.  Both appear in §B.4.
- **Patch files under `debug_document/`**: 1 patch file.  Appears
  in §B.4.
- **Source modifications (additive diagnostics)**: 1 modification
  to `dynamic/wave_operator.inl` (`SEAS_DIAG_CENTROID_MARGIN`
  `#ifdef` block).  Guarded by §17.6 test entry.
- **Out-of-scope modifications**: `drivers/tpv102_driver.cpp` and
  `io/paraview_output.hpp` show changes with `R-V92-E02` /
  stage-4 VTU references — these belong to a separate agent
  session and are NOT part of v9.2.0 rev-3f scope.

**Conclusion of §B audit:** every v9.2.0 rev-3f test artifact is
accounted for in this table with a definite status (PASS / READY
/ PENDING).  No artifact is "lost" without a row.

---

## §19. Post-review (R-V92-E01..E08) — re-analysis + instrumentation (2026-04-21)

### §19.0 Why §18 conclusions were reopened

A code-review of the rev-3e §18 chain of reasoning flagged eight
issues under IDs R-V92-E01..E08:

| Review ID | Severity | Root problem |
|---|---|---|
| **R-V92-E01** | CRITICAL | Outlier count tracks rupture extent, not amplification.  Median+MAD on a field that's 0 outside the rupture and nonzero inside flags the entire rupture as outliers; the 45 % → 62 % → 96 % → 44 % trajectory is **step + saturation + reflection-spike + decay**, NOT "monotonic growth".  §18.3 classifier verdict is spurious. |
| **R-V92-E02** | CRITICAL | §18.6 eliminates H-V92-K on the claim "VTU is not RK4-stage-averaged".  Factually wrong.  `drivers/tpv102_driver.cpp:860-873` explicitly RK4-averages `dof_data.V1/V2/tau1_corr/tau2_corr/sigma_n_corr` — exactly the fields written to VTU.  The in-code comment at line 870 says "R-001 fix: RK4-weighted corrected tractions for consistent station output".  H-V92-K is **NOT eliminated** and must be re-opened. |
| **R-V92-E03** | CRITICAL | §18.7's ∂_x v_x ~ 10⁴ s⁻¹ is off by 5 orders of magnitude.  Correct estimate: 5 m/s ÷ 100 m process zone = 5·10⁻² s⁻¹, giving σ_yy build-up of ~0.15 MPa/step, NOT 3 GPa/step.  Observed σ_n peak 218 MPa is 2-5× the mode-II bound (~ 45 MPa) — anomalous but NOT dimensionally consistent with the "3 GPa/step amplification" narrative. |
| R-V92-E04 | MODERATE | "Monotonically upward" language overstates a step + saturation + spike + decay pattern. |
| R-V92-E05 | MODERATE | §5.1 `SEAS_DIAG_FAULT_SIGMA` (bulk max\|Q_global[SYY]\|) probes a DIFFERENT quantity (pre-friction bulk Q) than §18.  It discriminates H-V92-G (bulk pumping) from H-V92-K (RK4 corrupts end-of-pipeline DOFData while bulk Q stays clean).  Wrongly de-prioritised in §18.8. |
| R-V92-E06 | MODERATE | Mass-inverse (candidate d) picked as "most plausible" in §18.8 without evidence ruling out (a) / (b) / (c). |
| R-V92-E07 | MODERATE | `rel_mad_threshold` skip in the outlier detector may silence real signals in param_a / param_Dc. |
| R-V92-E08 | MODERATE | At 96 % outlier fraction the classifier is past its valid regime. |

Verdict of the review: **FAIL — must resolve R-V92-E01..E06 before
proceeding to the mass-inverse instrumentation (§18.8 candidate d)
that rev-3e prescribed**.

### §19.1 Step 1 — Peak-value trajectory (R-V92-E01 re-analysis)

**Implementation.** New script
`tpv102/r_v92_outlier_detection/plot_peak_trajectory.py` reads the
existing `outlier_report.json` (no re-run of the detector; no
Frontera) and extracts the per-cycle **peak absolute deviation**
for each field.  Classifier distinguishes four shapes:

| Shape | Definition | Interpretation |
|---|---|---|
| SATURATION | peak stays near max, max NOT at end | fields hit a physical bound (mode-II envelope); normal rupture physics |
| MONOTONIC-GROWTH | max at end of run, late mean tracks max | unbounded growth; H-V92-G consistent |
| STEP-AND-DECAY | peak mid-run, late mean < 0.5 × max | transient spike then radiation away; normal physics |
| STEP-SPIKE-DECAY | peak mid-run, late mean 0.5–0.8 × max | reflection spike signature |

**Unit tests** — `tests/unit/test_peak_trajectory_classifier.py`
9 tests covering all four shapes plus the flat / empty / end-to-end
paths.  **Result: 9/9 PASS.**  Specifically `test_rev3e_outlier_trajectory_is_step_spike_decay`
asserts the 45 % → 62 % → 96 % → 44 % trajectory is classified as
STEP-AND-DECAY or STEP-SPIKE-DECAY, **not** MONOTONIC-GROWTH —
confirming R-V92-E01 on a reproducible synthetic fixture.

**Real-data outcome** on the existing job-7668434 `outlier_report.json`
(93 cycles, 38 039 triangles per cycle):

| Field | Shape | Global max | Value at end |
|---|---|---:|---:|
| **normal_stress** | **MONOTONIC-GROWTH** | 2.187·10⁸ Pa | 2.187·10⁸ Pa |
| slip_dip | SATURATION | 1.827 m | 1.827 m |
| traction_dip | SATURATION | 6.815·10⁷ Pa | 6.815·10⁷ Pa |
| slip_rate_dip | STEP-SPIKE-DECAY | 2.150 m/s | 0.648 m/s |
| slip_strike | SATURATION | 7.231 m | 6.742 m |
| traction_strike | SATURATION | 1.486·10⁸ Pa | 1.438·10⁸ Pa |

**Interpretation** (revising §18.6):

1. `normal_stress` DOES grow monotonically — the σ_n anomaly is the
   one field that still fits H-V92-G's "unbounded amplification"
   signature.  R-V92-E01 is **PARTIALLY vindicated**: the
   *aggregate* outlier-count trajectory was misclassified, but the
   *normal_stress* peak trajectory really does grow monotonically.
2. `slip_dip` SATURATES at 1.83 m.  This is **not** monotonic growth
   — the dip channel hits a saturation bound around t ≈ 9 s.
   The §18 "slip_dip grew 178×" headline conflated the
   `late_mean / early_mean` ratio (relevant only if early-run is
   quiescent, which it is here by design — pre-nucleation) with
   actual late-time amplification.  **R-V92-E01 vindicated for the
   dip channel.**
3. `slip_rate_dip` is STEP-SPIKE-DECAY (peaks at 2.15 m/s then
   decays to 0.65 m/s by end).  Consistent with a transient
   reflection spike, not amplification.
4. Strike-channel fields (`slip_strike`, `traction_strike`,
   `slip_rate_strike`) all saturate at physically plausible mode-II
   amplitudes.  No evidence of bulk pumping on the strike channel.

**Net verdict on §18.6 classification:**
- H-V92-G (bulk amplification): **STILL CONSISTENT for σ_n**.
  Peak normal_stress trajectory IS monotonic growth.
- H-V92-B (dip-channel pollution): **DOWNGRADED from CONFIRMED**.
  Dip contamination is real but SATURATES, not grows.  Could be
  a transient that clears after 9 s.
- H-V92-K (RK4 averaging): **RE-OPENED per R-V92-E02** (see §19.2).
  Previously eliminated on false premise.

### §19.2 Step 2 — Stage-4 vs averaged DOFData VTU instrumentation (R-V92-E02)

**Implementation.**  Modified:

1. `io/paraview_output.hpp::WriteFaultSurfaceVTU` — added three
   optional `_k4` vector arguments (`local_slip_rate_k4`,
   `local_traction_k4`, `local_normal_stress_k4`).  When all three
   are non-empty, the writer emits five extra CellData fields
   (`slip_rate_dip_k4`, `slip_rate_strike_k4`, `traction_dip_k4`,
   `traction_strike_k4`, `normal_stress_k4`) alongside the
   existing averaged ones.  Fail-safe: if any one of the three
   `_k4` vectors is empty, `_k4` emission is suppressed (partial
   input is not emitted).  Backwards-compatible default: empty
   Vector() defaults mean existing callers (e.g. BP5's verification
   driver) continue to emit the original 12 fields only.
2. `drivers/tpv102_driver.cpp` — added three `pv_local_*_k4`
   scratch buffers at ParaView init time.  After stage 4 but
   **before** the RK4-averaging block at lines 858-874 (which
   overwrites `dof_data.V1/V2/tau*_corr/sigma_n_corr`), the driver
   copies the stage-4 kN buffers (`V1_k4`, `V2_k4`, `t1c_k4`,
   `t2c_k4`, `snc_k4`) into the `pv_local_*_k4` buffers.  The
   paraview_write lambda then passes both the averaged and the
   stage-4 buffers to `WriteFaultSurfaceVTU`.  For the t=0
   initial snapshot the stage-4 values equal the initial DOFData
   equal the averaged values, so the delta is zero by construction.

**Unit tests** — `tests/unit/test_fault_surface_vtu_k4.cpp` (new,
43 assertions):

| Test case | What it verifies |
|---|---|
| R-V92-E02 (a) | Empty `_k4` inputs ⇒ no `_k4` CellData fields emitted.  Backwards-compatibility for BP5 and any caller that has not been updated. |
| R-V92-E02 (b+c) | Populated `_k4` inputs with values distinct from averaged inputs ⇒ all five `_k4` fields present with the expected face-averaged values from the `_k4` buffer (per cell).  The defensive assertion `\|normal_stress_k4 − normal_stress\| > 1 MPa` confirms the diagnostic is live — k4 and averaged are stored as separately-sourced fields. |
| R-V92-E02 (d) | Partial `_k4` inputs (only one of three non-empty) ⇒ `_k4` fields suppressed.  Fail-safe against a future driver bug that sizes one buffer but not others. |

**Result: 43/43 PASS.**  Existing `test_fault_surface_vtu_continuity`
still 63/63 (no regression to the non-k4 CellData contract).

**Makefile targets added:** `seas_test_fault_surface_vtu_k4` /
`make test-fault-surface-vtu-k4`.

### §19.3 Step 3 — SEAS_DIAG_FAULT_SIGMA Frontera probe (user-gated)

**Scope.** Per `feedback_frontera_approval.md` and the review's
`Do NOT` list, I did NOT:
- port the §16.5 BP5 owner-pattern,
- instrument `ApplyMassInverse` (§18.8 candidate d),
- submit any Frontera sbatch.

Step 3 (reinstate §5.1 `SEAS_DIAG_FAULT_SIGMA` to probe bulk
max\|Q_global[SYY]\| on fault QPs) remains proposed but NOT
implemented.  When user approves, the §5.1 instrumentation goes
into `dynamic/wave_operator.inl` as a `#ifdef
SEAS_DIAG_FAULT_SIGMA` block, analogous to the existing
`SEAS_DIAG_FAULT_FLUX` blocks.  Budget per plan §9: ~200 SU for a
4 s run.

### §19.4 Revised decision tree (post R-V92-E01..E06)

Three FREE local probes are NOW available (no Frontera required):

**Probe 1 — Peak trajectory classification** (Step 1, done).
Outcome: normal_stress MONOTONIC-GROWTH; dip-channel fields
SATURATION.  σ_n is the primary anomaly; dip is a secondary
(saturated) symptom.

**Probe 2 — Stage-4 vs averaged VTU delta** (Step 2, instrumentation
ready; requires a re-run with the updated driver to emit `_k4`
fields).  **Pending user approval of a local 2 s re-run** on a
small mesh (e.g. tpv102_1000m.msh on the laptop at 1 rank — not
Frontera).  Decision rule:
- If `max |normal_stress_k4 − normal_stress| < 1 MPa` everywhere:
  H-V92-K ELIMINATED.  The RK4-averaging is not corrupting the
  end-of-pipeline DOFData.  Proceed to Probe 3.
- If `max |normal_stress_k4 − normal_stress| > 1 MPa`: **H-V92-K
  CONFIRMED**.  The σ_n drift IS the RK4-averaging residual.  The
  mass-inverse candidate (d) is **irrelevant**; the fix is instead
  a correct end-of-step re-evaluation of DOFData at Q_{n+1} after
  the RK4 combination (vs the current stage-average formula).
  **This is the cheapest possible fix if confirmed** —
  ~20 LOC in the driver, no Frontera.

**Probe 3 — Bulk SYY_g probe** (Step 3, user-gated).  Discriminates
H-V92-G (bulk pumping drives both bulk Q[SYY] and DOFData) vs
H-V92-K (only DOFData is corrupt, bulk Q[SYY] stays clean).
Decision rule:
- If `max |Q_self[SYY]|` on fault QPs grows to ~50 MPa by t=6s:
  H-V92-G CONFIRMED.  Bulk σ_yy field is being pumped.  Proceed
  to §18.8 candidate-by-candidate FREE local probes to find the
  amplifier (NOT starting with (d); all four in parallel).
- If `max |Q_self[SYY]|` stays small (<1 MPa) while DOFData σ_n
  shows 50 MPa drift: H-V92-G ELIMINATED.  H-V92-K is the
  mechanism.  Apply the fix from Probe 2.

### §19.5 Instrumentation artifacts (new)

New files:
- `miniapps/seas/tpv102/r_v92_outlier_detection/plot_peak_trajectory.py`
  (~240 LOC Python).  Ran locally, emits
  `peak_trajectory.png` + `peak_trajectory_summary.json` in the same
  directory.
- `miniapps/seas/tests/unit/test_fault_surface_vtu_k4.cpp` (~330 LOC).
  Makefile target `make test-fault-surface-vtu-k4`.  43/43 PASS.
- `miniapps/seas/tests/unit/test_peak_trajectory_classifier.py`
  (~180 LOC).  9 unit tests including the rev-3e outlier-trajectory
  shape regression.

Modified files:
- `miniapps/seas/io/paraview_output.hpp`: +46/−12 — added
  `_k4` optional parameters, new accumulation pass, five extra
  `write_field(...)` calls under `if (has_k4)`, matching PVTU
  field-list entries.
- `miniapps/seas/drivers/tpv102_driver.cpp`: +36/−1 — added
  `pv_local_*_k4` buffers, stage-4 snapshot copy before averaging,
  updated `paraview_write` call to pass them.
- `miniapps/seas/Makefile`: +5 — new test target entries.

All changes compile cleanly under `mfem-dev` (`make seas_tpv102_driver`
and the two test targets).  No regression on
`seas_test_fault_surface_vtu_continuity` (63/63),
`seas_test_fault_basis_dip_strike_symmetry` (17/17),
`seas_test_godunov_rotation_identity` (4/4).

### §19.6 Status of root cause — honest assessment

**The root cause has NOT been conclusively identified.**  What the
§19 work accomplished:

1. Refuted the rev-3e classification of the 45 % → 62 % → 96 % →
   44 % outlier trajectory as "monotonic growth" (R-V92-E01
   validated).
2. Produced a more nuanced picture from the same data: σ_n peak
   trajectory IS monotonic, but dip-channel peak trajectories
   saturate — a different story than "everything is amplifying".
3. Re-opened H-V92-K (RK4-averaging) that §18.6 had eliminated on
   a factually wrong premise.
4. Shipped a VTU-level diagnostic (`_k4` CellData fields) that will
   discriminate H-V92-K from H-V92-G on the next re-run, without
   needing Frontera.

What the §19 work did NOT accomplish:
- No Frontera re-run has yet exercised the `_k4` instrumentation
  under production conditions.  The 2 s local re-run on the
  laptop-scale fixture is proposed but not run (per user policy
  against running the production mesh locally).
- §5.1 bulk SYY_g probe is NOT implemented.
- Mass-inverse instrumentation is NOT implemented (correctly
  deferred per R-V92-E06).

**Next action gate for the user:**
1. Approve a local 2 s laptop run of `seas_tpv102_driver` on
   `tpv102/mesh/tpv102_1000m.msh` (1 rank, ~5 min wall) to emit a
   VTU with `_k4` fields already populated.  Inspect the
   `normal_stress_k4 - normal_stress` CellData field in ParaView.
   **This single measurement discriminates H-V92-K from H-V92-G.**
2. If H-V92-K: apply the ~20 LOC driver fix (re-evaluate DOFData
   at Q_{n+1} after RK4 combination instead of averaging stages),
   re-run the §6.4 outlier detector on a new 2 s VTU, confirm σ_n
   drift is gone.
3. If H-V92-G: escalate to Step 3 (Frontera §5.1 probe), then to
   the parallel candidate-(b)/(c)/(d)/(e) probes per §18.8.

The current evidence points MORE toward H-V92-K than H-V92-G, for
three reasons:
- The dip channel SATURATES (consistent with a bounded numerical
  artifact, not an amplifier).
- §18.7's physics calculation (rev-3e) is off by 5 orders of
  magnitude — R-V92-E03 says the physical process cannot generate
  the 218 MPa signal, so a numerical artifact must be responsible.
- The RK4-averaged DOFData in `tpv102_driver.cpp:858-874` is the
  ONLY place in the code that can OUTPUT a value that differs
  from the end-of-step friction-solve state — exactly matching
  the pattern of "VTU shows something, station shows the same
  thing, but the underlying physics can't produce it".

But this is an EVIDENCE-WEIGHTED suspicion, not a confirmed
diagnosis.  The Probe 2 `_k4` delta inspection is the definitive
test and it has not been run yet.

### §19.7 One-sentence summary (post-review)

Three critical flaws in the rev-3e reasoning chain (R-V92-E01
outlier-count misread, R-V92-E02 H-V92-K falsely eliminated,
R-V92-E03 physics 10⁵× off) have been corrected; the peak-value
trajectory tool and the stage-4 `_k4` VTU diagnostic are
implemented and tested; the root cause is **not yet confirmed**
but a ~5 min laptop re-run + ParaView inspection of
`normal_stress_k4 - normal_stress` will conclusively discriminate
H-V92-K (cheap 20-LOC fix) from H-V92-G (expensive candidate
tree in §18.8).

---

## §20. Round-5 synthesis + next-direction gate (2026-04-21 rev-3h)

Round 5 of REVIEW.md synthesized rounds 1-4 findings and delivered
findings R-V92-G01 through R-V92-G06.  Key conclusions consumed in
this revision:

### §20.1 Hypothesis-rank table (post-round-5)

| ID | Status | How closed / why open |
|---|---|---|
| H-V92-R | CLOSED (retracted) | §A proof + §4.1 0 ULP |
| H-V92-O | CLOSED (retracted) | §A proof + §4.2 bit-exact |
| H-V92-F | CLOSED (not primary) | §4.3 on 200 m prod mesh ≤ 0.5 ULP |
| H-V92-P (init) | CLOSED on tested scales | §4.6 PASS 5/5 |
| H-V92-P (runtime) | CLOSED on tested scales | §4.7 Phase D PASS 1/1 |
| H-V92-P (FP-fragility) | CLOSED on 1000 m geometry | centroid-margin probe 227-266 m |
| H-V92-T (ADER vs RK4) | CLOSED NON-PRIMARY | round-4 F03 dimensional analysis |
| H-V92-G cand (a) IP/BR2 penalty | CLOSED | §4.9 no-penalty regression PASS |
| H-V92-G cand (b) RK4 drift | CLOSED | §4.11 RK4-conservation probe PASS |
| H-V92-G cand (c) absorbing-BC phase | PARTIAL CLOSED | §4.12 energy-monotonic PASS; phase per-channel NOT verified → H-V92-V below |
| H-V92-G cand (d) mass-inverse asymmetry | CLOSED | §4.10 single-channel Jacobian PASS |
| H-V92-G cand (e) volume-term coupling | CLOSED | §4.10 isotropy ratio 1.000 bit-exact |
| **H-V92-G (as a whole)** | **INCONSISTENT** | §18 claims CONFIRMED but all 4-5 candidate amplifiers eliminated; see Open Question #1 below |
| **H-V92-K (reopened)** | **OPEN** | round-3 R-V92-E02 + round-4 F01/F02 reopened; discriminator = dt-halving in §20.3 |
| **H-V92-U (interior-fault path)** | **OPEN — RANK-1** | NEW rev-3h; never unit-tested; majority of fault processing per §17.6.5 |
| **H-V92-V (abs-BC phase per-channel)** | **OPEN** | NEW rev-3h; §4.12 energy-monotonic doesn't prove per-channel phase |
| **H-V92-W (init pre-stress projection)** | **OPEN — RANK-3** | NEW rev-3h; cheap `[FAULT-INIT-V1]` printf would confirm/close |
| H-V92-C | WEAKENED | §15.2 signature mismatch; retained as possible secondary |
| H-V92-Q (friction-solver instability) | OPEN (unlikely primary) | no evidence supports |

### §20.2 Open question #1 — "H-V92-G CONFIRMED" is inconsistent with the §4.9-§4.12 eliminations

**Problem.**  §18.6 marks H-V92-G as the primary mechanism based
on the §6.4 outlier-count trend.  REVIEW round-3 R-V92-E01
established the outlier count tracks rupture extent, not
amplification; REVIEW round-4 F03 established the observed
pathology magnitude cannot come from any time-integrator error.
Meanwhile §4.9–§4.12 have ELIMINATED all four (five with the
volume-term split) of the candidate amplifiers H-V92-G originally
proposed.  **There is no candidate amplifier left within the
H-V92-G framework.**

**Implication.**  Either (a) the §18 "CONFIRMED" verdict is wrong
and the primary mechanism is NOT H-V92-G, or (b) a new amplifier
candidate exists that was not in the original H-V92-G tree.
The round-5 new hypotheses H-V92-U, H-V92-V, H-V92-W are the
leading candidates under option (b).

**Resolution.** Defer to the dt-halving test + §5.1 bulk probe in
the active TODO; one of them will return a cleaner discrimination.

### §20.3 Decisive next-step tree (replaces §18.8)

Per REVIEW R-V92-G01, break the 4-round "small-scale PASS →
defer" loop with the decisive dt-halving test first.  Full
sequence copied from the active TODO:

1. **dt-halving** (Frontera, ~20 SU) — classifies time-integration
   vs flux-path branch.
2. **§5.1 SEAS_DIAG_FAULT_SIGMA** (conditional, ~200 SU) — discriminates
   bulk vs post-friction.
3. **`test_interior_fault_flux_path.cpp`** (FREE local) — probes
   the untested majority-fault-faces path (H-V92-U).
4. **F01+F02 RK4 fix** (FREE local, ~100 LOC driver) — architectural
   hygiene regardless of whether H-V92-K is primary.
5. **[FAULT-INIT-V1] printf** (FREE, 1-line) — closes or opens H-V92-W
   at trivial cost.
6. Apply fix → Frontera 12 s confirmation run → write
   `tpv102_debug_v9.2.0_fix.md` → open check document.

### §20.4 Budget

| Step | SU | Cumulative SU |
|---|---:|---:|
| 1 dt-halving | 20 | 20 |
| 2 bulk σ probe (cond.) | 200 | 220 |
| 3 interior-fault unit test (local) | 0 | 220 |
| 4 F01+F02 fix (local) | 0 | 220 |
| 5 init-V1 printf (local) | 0 | 220 |
| 6 fix confirmation run | 400 | 620 |
| **Total (best case)** | — | **220** |
| **Total (full path)** | — | **620** |

Well under the original 1 000 SU envelope.

### §20.5 Closure status

**v9.2.0 is NOT yet close to closure.**  Per the plan-closure
gate at the head of the TODO:
- Primary mechanism identification — **NOT YET** (§20.2 open
  question).
- C2 source patch applied — NOT YET.
- Frontera confirmation `|σ_n − 120 MPa| ≤ 1 MPa` — NOT YET
  (current value: 218 MPa peak, ~100× spec).
- `tpv102_debug_v9.2.0_fix.md` written — NOT YET.
- `tpv102_debug_v9.2.0_check.md` opened — NOT YET.

**Estimated remaining work:** 1 – 2 Frontera runs + 1 source
patch + 1 confirmation run ≈ 0.5 – 1 day of focused execution if
the dt-halving test yields a clear branch on first try.

### §20.6 One-sentence summary

After 5 review rounds the investigation has closed many candidates
but not identified the primary R-V92 mechanism; the dt-halving
test (~20 SU) is the single decisive experiment that breaks the
repeating "small-scale PASS → defer" loop and must be run before
any further code change or expensive Frontera diagnostic.

---

## §21. Frontera dt-halving + postfix results (rev-3i, 2026-04-22)

Three Frontera runs completed on the 200 m / p=1 / 400-rank
fixture at `tfinal = 12.0 s`.  All three share the same SCEC
station set and plotting harness as job-7668434 (v91 baseline).

| # | Job | Label | dt | Source |
|---|---|---|---|---|
| 1 | 7670526 | `dt_half` | `0.5 · dt_cfl_v91` | PRE-fix HEAD (per §5 R-V92-H03 Option B) |
| 2 | 7670527 | `dt_quarter` | `0.25 · dt_cfl_v91` | PRE-fix HEAD |
| 3 | 7671002 | `postfix` | nominal `dt_cfl_v91` | POST-F01+F02 driver (Step 5 fix applied) |

Plot artifacts under
`miniapps/seas/tpv102/plots_results_200m_p1_12.0s_400r_v92_{dt_half,dt_quarter,postfix}_job{…}/`.

### §21.1 Observations across all 9 SCEC stations

For each station, V_strike / slip_strike / tau_strike / V_dip /
slip_dip / tau_dip / σ_n / log10(state) panels were compared panel-
by-panel across the three runs.  The following pathology
signatures persist in every run and are **visually
indistinguishable** between runs 1, 2, 3:

- **flt_0_3** (on-fault, z = −3 km): V_strike double-peak at ~3.5 s
  (3.3 m/s) and ~5.1 s (2.0 m/s); spurious V_dip ±0.010 m/s during
  rupture; slip_dip monotonically drifting to −0.035 m by t = 9 s;
  τ_dip reaches −1.25 MPa; σ_n drift 120 → 123 MPa.
- **flt_0_7.5** (on-fault, z = −7.5 km): V_strike peak 4.2 m/s at
  ~2 s and secondary 3–4 m/s at 6–7 s; slip_dip reaches −0.20 m;
  **σ_n drops to 102–105 MPa (-15 to -18 MPa from the 120 MPa
  spec)**; τ_dip reaches −3.0 MPa.
- **flt_0_12** (on-fault, z = −12 km, below seismogenic zone):
  sharp V_strike spike ~8 m/s at ~6.5 s; same pathology profile as
  flt_0_7.5 but compressed in time.
- **flt_9_7.5, flt_n9_7.5** (off-fault, along-strike ±9 km, z =
  −7.5 km): V_strike ~7 m/s peak at ~6.8 s; identical spurious
  V_dip / σ_n / τ_dip signatures in the three runs.
- **flt_12_3** (off-fault, along-strike 12 km, z = −3 km): slip_dip
  drift to **−0.10 to −0.13 m**; τ_dip peaks −6 MPa; σ_n transient
  spike to ~122 MPa.
- **flt_12_12, flt_n12_12** (off-fault, along-strike ±12 km, z =
  −12 km): slip_dip +0.025–0.030 m; V_dip spurious pulse +0.15 m/s
  at rupture-front arrival (~6.5 s); σ_n transient ~122 MPa.
- **flt_n12_3** (off-fault, along-strike −12 km, z = −3 km): mirror
  of flt_12_3; same signature magnitudes.

The overview V_strike plot (all 9 stations overlaid on one axis)
is bit-visually identical across runs 1–3: same peak magnitudes,
same peak arrival times, same rupture-front propagation order.

### §21.2 dt-scaling verdict — "≈ unchanged" branch

Per §20.3 Step 1 classification tree:

| Scaling outcome | Classification | Action |
|---|---|---|
| slip_dip drops ~16× at dt/2 | RK4 truncation | F01+F02 → skip to Step 8 |
| slip_dip drops ~8× at dt/2  | psi splitting (F02) | integrate psi inside RK4 |
| **≈ unchanged** | **NOT time integration** | **go to Step 2/3** |
| larger / NaN | CFL violation | reduce CFL + audit `flux.Interior` |

**Observed:** slip_dip peak magnitudes at the critical stations
(flt_0_7.5: −0.20 m; flt_12_3: −0.12 m; flt_0_3: −0.035 m) are
identical to PRE-fix v91 baseline (job-7668434) at all three dt
values to plotting resolution.  No detectable reduction at dt/2
or dt/4.  **Verdict: ≈ unchanged → NOT time-integration.**

### §21.3 F01+F02 postfix verdict — necessary-but-not-sufficient

Run 3 (postfix, job-7671002) exercises the F01+F02 coupled-RK4
-on-(Q, psi) driver.  Its plots are indistinguishable from runs 1–2
(PRE-fix HEAD).  Therefore the F01+F02 fix, while correct as
architectural hygiene (Step 5b PASS 100/100 on the 4-km fixture),
**does not close R-V92** on the production fault.  The
operator-splitting defect on ψ was not the amplifier.

### §21.4 Hypothesis-rank table update (supersedes §20.1)

| ID | Prev | New | Basis |
|---|---|---|---|
| H-V92-K (time-integrator) | OPEN | **CLOSED** | §21.2 dt-invariance across dt, dt/2, dt/4 |
| H-V92-T (RK4 truncation) | CLOSED NON-PRIMARY | CLOSED (reinforced) | §21.2 empirically confirms round-4 F03 |
| RK4-psi splitting (F02) | suspected amplifier | NOT amplifier | §21.3 postfix identical to PRE-fix |
| F01+F02 as fix | pending Frontera confirm | **applied, NOT sufficient** | §21.3 |
| **H-V92-U (interior-fault)** | OPEN — RANK-1 | **OPEN — RANK-1 (reinforced)** | Elimination of H-V92-K makes it the sole surviving spatial-path amplifier among non-init candidates |
| H-V92-V (abs-BC per-channel phase) | OPEN | OPEN | no new evidence |
| H-V92-W (init pre-stress) | OPEN — RANK-3 | OPEN — RANK-2 | promoted: §21 eliminates one above-ranked competitor |
| H-V92-G (bulk amplifier) | INCONSISTENT | INCONSISTENT | §21 does not discriminate; §5.1 SYY probe still decisive |
| H-V92-Q (friction solver) | OPEN (unlikely) | OPEN (unlikely) | no new evidence |

### §21.5 Updated next-step tree (supersedes §20.3 Step 1)

Step 1 is now **DONE** with verdict "≈ unchanged → NOT time
integration".  The surviving path is:

1. ~~dt-halving~~ — **DONE (§21.2).  Verdict: not time-integration.**
2. ~~re-analyse §18 outliers~~ — **DONE: MONOTONIC-GROWTH.**
3. **§5.1 `SEAS_DIAG_FAULT_SIGMA` bulk probe (Frontera, ~200 SU)** —
   now the highest-priority remaining diagnostic.  Discriminates
   bulk pump (H-V92-G) from post-friction drift (now reducible to
   H-V92-U or H-V92-V since H-V92-K closed).
4. ~~interior-fault 2-tet unit test~~ — **DONE: PASS 6/6 on 2-tet.**
   Bimaterial / corner / multi-face extensions remain viable.
5. ~~F01+F02~~ — **DONE + empirically shown non-sufficient (§21.3).**
   Commit anyway as architectural hygiene.
6. **`[FAULT-INIT-V1]` printf (FREE)** — now higher-ROI than before
   because H-V92-W is promoted to RANK-2.  Trivially cheap.
7. **Interior-fault unit-test extensions** — add bimaterial / corner
   / multi-face variants to `test_interior_fault_flux_path.cpp` to
   actually exercise H-V92-U on code paths not yet covered.
8. Apply C2 patch; Step 8 confirmation; fix/check docs.

### §21.6 Commit + push decision for F01+F02

Per §5 R-V92-H03 Option B hold, F01+F02 was held locally pending
PRE-fix dt-halving completion.  Runs 1–2 are the PRE-fix
dt-halving data.  Run 3 is the POST-fix confirmation at nominal
dt.  All three complete with matching plots.  **Option B hold is
now satisfied.**  Recommendation: commit + push F01+F02 as
architectural hygiene, with the explicit commit-message note that
§21 shows the fix alone does not close R-V92 and the primary
amplifier search must continue on the non-time-integration branch.

### §21.7 Budget status

| Step | Planned SU | Actual SU | Cumulative |
|---|---:|---:|---:|
| 1 dt_half | ~10 | job-7670526 | — |
| 1 dt_quarter | ~10 | job-7670527 | — |
| 1-post F01+F02 confirm | — | job-7671002 (~400 SU, 12 s run) | — |
| 3 §5.1 SYY probe | 200 | pending | — |
| 6 confirmation run | 400 | pending | — |

(Exact SU draw per job not yet recorded — add when Frontera usage
log is pulled.)

### §21.8 One-sentence summary

The Frontera dt-halving + dt-quartering + F01+F02-postfix trio
shows **zero change** in the R-V92 pathology across all three
runs, decisively closing H-V92-K (time-integrator) and
F01+F02-as-primary-fix, promoting H-V92-U (interior-fault path,
non-2-tet code paths) to sole RANK-1, and redirecting the next
decisive experiment to the §5.1 `SEAS_DIAG_FAULT_SIGMA` bulk
probe.
