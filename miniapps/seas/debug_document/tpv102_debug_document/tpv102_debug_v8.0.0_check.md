# Code Review v8.0.0: rupture pinned at hypocenter — fault-to-bulk coupling broken after R-801 convention flip

**Versioning note:** Previous `v7.0.0`.  This round found **one CRITICAL bug**
(rupture fails to propagate away from the hypocenter; fault slip is not
radiating bulk waves at the correct magnitude).  Bumps `[a]`.
**New version: `v8.0.0`**.

## TL;DR — v7.1.0 was premature; do not apply the RESULT.txt regex fix

The Frontera 400-rank 200m dispositive run (job 7665297) **physically
failed** in a way the v7.0.0 test suite was blind to:

- Hypocenter fault QP reaches `V_max = 7.7 m/s` (nucleation → breakaway
  works **locally**).
- Every one of the nine SCEC fault stations **off** the hypocenter
  remains at `V2 = V_ini = 1e-12 m/s` and `tau2 = tau_ini = 75 MPa`
  (pristine initial state) at `tfinal = 1.5 s`.
- Maximum `‖Q‖_∞` across all 88.7M bulk DOFs = `1.66e-8`.  For a
  physical 7.7 m/s slipping fault, bulk velocity near the fault
  should be O(m/s) and stress perturbations O(10 MPa = 1e7 Pa).
  **~10^15 orders of magnitude too small.**

The qnorm watch going "green" (r4, r0 climbing by orders of magnitude)
is not a sign of healthy propagation — it's low-amplitude numerical
leakage from a broken coupling.  The R-101 verifier, R-302a, R-501a,
T-R801, and T-R802 tests all pass because they check *cross-rank
DOFData consistency* and *zero-sum flux*, not physical coupling
magnitude.  The fault-to-bulk coupling is broken and none of the
tests can detect it.

**The v7.1.0 RESULT.txt regex fix would mask this by turning the FAIL
into a PASS.  Do NOT apply v7.1.0 R-901.**

## Review Scope
- Log: `tpv102_200m_400r_v2_7665297.out` (full run, 5268 steps).
- Station files under
  `/Users/chunhuizhao/Downloads/seas-mfem/tpv102/results_200m_p1_1.5s_400r_v2_job7665297/`.
- Code path that may contain the root cause:
  - `miniapps/seas/dynamic/fault_face_flux.cpp::Evaluate` — trial-correct
    Riemann solver.
  - `miniapps/seas/dynamic/wave_operator.inl` — interior-fault branch
    (lines 705-844) and shared-fault branch (lines 964-1213), both use
    the canonical-frame reconstruction from `FaultBasisQPData::sign_flipped`.
  - `miniapps/seas/dynamic/godunov_flux.cpp::Interior` — Godunov
    upwind; internal frame construction via `BuildFrame`.
  - `miniapps/seas/dynamic/tpv102_setup.hpp::InitializeFaultDOFs` +
    `ApplyNucleation` — now write `tau2_0 = tau_ini, V2 = V_ini`
    (BP5 convention, component 2 = strike).
- v7.0.0 audit findings that are now superseded.

## Empirical evidence

### 1. Fault station data — only the hypocenter ruptures

At `tfinal = 1.5 s`, column V2 (strike slip rate, BP5 convention):

| Station | Position (strike, down-dip) km | Distance from hypo | V2 (m/s) |
|---|---|---|---|
| **flt_0_7.5** | **(0, 7.5)** | **0 km (HYPOCENTER)** | **7.699** |
| flt_0_3 | (0, 3) | 4.5 km up-dip | 1.0e-12 (V_ini) |
| flt_0_12 | (0, 12) | 4.5 km down-dip | 1.0e-12 |
| flt_9_7.5 | (9, 7.5) | 9 km along-strike | 1.0e-12 |
| flt_n9_7.5 | (-9, 7.5) | 9 km along-strike | 1.0e-12 |
| flt_12_3, flt_12_12, flt_n12_3, flt_n12_12 | ±12km, ±3/12 km | 13 km | 1.0e-12 |

**Every station except the hypocenter is at numerically-exact initial
V_ini.**  Time-series shows the off-hypocenter stations' `V2`
*decrease* slightly from `V_ini = 1.00e-12` to `9.9999999868e-13`
over 1.5 s — consistent with plate-rate state variable relaxation on
a locked fault QP that receives no dynamic waves.  **No rupture front
arrives anywhere off the hypocenter.**

### 2. Hypocenter V_max growth + saturation

```
t=0.3s  V_max=2.1e-08  (slow slip, nucleation loading)
t=0.7s  V_max=1.7e-02  (approaching V_nuc)
t=0.9s  V_max=2.1e-01  (breakaway at hypo)
t=1.1s  V_max=2.4      (acceleration)
t=1.2s  V_max=7.70     (peak)
t=1.3s  V_max=7.69871  (saturated, bit-identical)
t=1.4s  V_max=7.69871
t=1.5s  V_max=7.69871
```

Bit-identical V_max at 7.69871 over 1000+ time steps.  Originally
I interpreted this as "steady dynamic propagation"; in light of the
station data, the correct interpretation is **the hypocenter QP
converged to its friction-law steady-state slip rate** for the
imposed `tau_ini + peak_dtau` + equilibrium psi, and the rupture
does not propagate away to any other QP.

### 3. Bulk qnorm vs fault slip — the smoking gun

| Quantity | Observed | Expected for V_max=7.7 m/s |
|---|---|---|
| Max \|Q\|_∞ (any component, any DOF) | `1.66e-08` | `O(1)` m/s (velocity) or `O(1e7)` Pa (stress) |
| Gap | | **~10^15 orders of magnitude** |

Physics of dynamic rupture: a fault slipping at `V = 7.7 m/s` with
shear modulus `mu ≈ 32 GPa` and shear wave speed `c_s ≈ 3464 m/s`
radiates body waves with characteristic stress amplitude
`sigma ≈ rho * c_s * V/2 ≈ 2670 * 3464 * 3.85 ≈ 3.6e7 Pa = 36 MPa`
and particle velocity `V/2 ≈ 3.85 m/s`.  These are the bulk state
perturbations **adjacent to the fault**.  The observed bulk maximum
is `1.66e-8`, about **10^15× below** the physically-expected value.

This is not a resolution artifact (`Λ_dyn/h ≈ 0.8` at 200m is
borderline, not obliterating).  It is evidence that **the fault
flux is not injecting bulk wave energy at physically-correct
magnitude**.

### 4. Why all v7.0.0 tests pass anyway

- R-101 verifier: compares `DOFData` across ranks for shared QPs.
  DOFData is identical on both ranks (405 pairs at 3 ULP).  The
  verifier says nothing about whether the DOFData values themselves
  are the physically-correct result of the fault flux.  **Blind to
  this bug.**
- T-R302a: `Q=0` initial, one Mult.  Tests pairing only.  **Blind.**
- T-R501a: nonzero Q, 4 Mults, pair-consistency check.  Q is a
  prescribed sinusoid, not an evolved rupture.  **Blind.**
- T-R801 (Parts A, B): asserts `|V2| > |V1|` after one Mult on
  `Q=0`.  With `tau2_0 = tau_ini` (strike) the V2 component IS the
  larger one even if propagation is broken.  **Blind.**
- T-R802: `Σ_ranks k_shared_face ≈ 0` at `Q=0`.  Tests conservation
  of flux to machine zero when there's essentially nothing to flux.
  **Blind.**

The test suite verifies *cross-rank consistency* and *zero-sum flux*
but **not** *physical magnitude of fault-to-bulk coupling*.  This is
a genuine gap in the v7.0.0 test strategy.

## Findings

### [R-1001] [CRITICAL] Fault-to-bulk coupling fails to radiate bulk waves at physical magnitude; rupture pinned at hypocenter after R-801 convention flip

**Category:** BUG (physics-level)

**Description:**
After the R-801 Option A refactor (v7.0.0) unified the whole TPV102
pipeline on BP5's canonical fault-local frame (`tangent1 = dip,
tangent2 = strike`), the hypocenter's nucleation zone ruptures
correctly *locally* (V_max reaches physical steady state), but no
physically-meaningful bulk wave energy reaches any off-hypocenter
fault QP.  Bulk `Q` globally has `max‖Q‖_∞ = 1.66e-08` — 10^15× below
what a physically-correct fault → bulk coupling would produce for a
V=7.7 m/s slipping fault.

Since the trial traction at neighboring fault QPs is computed from
their local `Q_self, Q_nbr` (and these are essentially zero), they
never exceed yield and never rupture.  The result is a pinned-at-
hypocenter slow-slip nucleus that reaches friction-law equilibrium
and stops — not a propagating dynamic rupture.

**Trigger:**
Every TPV102 run under the v7.0.0 code.  The bug may predate v7.0.0
but is newly *unmasked* at Frontera scale; earlier runs (4-rank
1000m, 8-rank 1000m) were at Q-amplitudes too small to
discriminate, and shorter physical time windows (tfinal = 0.1 – 0.5 s
on the smoke tests) never reached full breakaway anyway.  The v7.0.0
audit's 4-rank smoke run at `--tfinal 0.1` completed without abort
precisely because the bug doesn't show up until after breakaway at
~t=1s.

**Actual behavior:**
See empirical section above.

**Expected behavior:**
At `tfinal = 1.5 s` on a 200m 400-rank TPV102 configuration:
- Fault stations at ±4.5 km (flt_0_3, flt_0_12) should show arrival
  of the rupture front by ~t = distance/c_s ≈ 1.3 s, with V2
  climbing into O(1 m/s).
- Bulk `max‖Q‖_∞` should be O(MPa) for stress perturbations near the
  rupture front.
- R-402 regex's "r399 dead" check should NOT be the failure mode —
  the correct failure is *stations never slipped*.

**Hypotheses for root cause (in order of likelihood):**

1. **H1:** R-801 Option A moved strike from `tau1_0` to `tau2_0` in
   `InitializeFaultDOFs` and `ApplyNucleation`, but the rotation
   logic in `Evaluate` → fault-local frame → back to global may still
   be routing the strike response through the wrong tangent axis.
   Steady-state hypocenter slip works (Brent solver is driven by
   `Theta = sqrt(tau1_total² + tau2_total²)` which is axis-swap-
   invariant), but the bulk-radiated component (`Q_imp[VZ]` vs
   `Q_imp[VY]`) might be tiny or cancellation-prone under the new
   convention.  **Probable.**

2. **H2:** The R-802 `accum_sign = elem1_on_plus ? +1 : -1` for the
   shared-fault flux is correct for conservation, but in combination
   with how the interior-fault branch handles the non-canonical
   `nor` accumulation there could be a sign-cancellation.  A fault
   face is EITHER interior-local (both elements on one rank) OR
   shared; the pairs are mutually exclusive.  But if both branches
   write into `rhs` for different faces with flipped conventions,
   bulk perturbation at a given DOF sums to near-zero.  **Possible.**

3. **H3:** Pre-R-801 the interior-fault branch used
   `GodunovFlux::BuildFrame` (t1=strike for TPV102) which was
   COINCIDENTALLY correct because the driver was ALSO writing strike
   into tau1_0.  R-801 unified everything on BP5 convention, but the
   **Q_imp reconstruction formula** (`Q_imp_plus[VY] = Q_plus[VY] +
   invZs * (tau1_corr - Q_plus[SXY])`) was not re-audited.  If that
   formula has a subtle sign or index error that was self-consistent
   under the OLD convention but not the new one, strike-component
   bulk velocities come out near zero.  **Possible.**

4. **H4:** Not a code bug — 200m mesh really is too coarse to
   propagate TPV102 rupture when `Λ_dyn/h ≈ 0.8`.  But: the v1
   debug documents and prior Tandem benchmarks show propagation
   IS achievable at 200m.  And the 10^15× magnitude gap is beyond
   what under-resolution could explain — under-resolution produces
   jagged or noise-limited rupture, not total absence of bulk-wave
   radiation.  **Unlikely but cannot be ruled out without comparing
   to a Tandem 200m run.**

**Required diagnostic before choosing a fix:**
Before attempting a code fix, do ONE of the following to isolate the
root cause:

- **D1:** Write a unit test that drives `FaultFaceFlux::Evaluate` on
  a SINGLE fault QP with hand-picked `DOFData` representing steady
  TPV102 hypocenter state (`tau2_0 = tau_ini + peak_dtau`, other
  fields at equilibrium).  Check that `Q_imp_plus[VZ]` (strike
  component of imposed bulk velocity, fault-local) is the expected
  `~ tau2_corr / Zs ≈ 6-7 m/s`.  If `Q_imp_plus[VZ]` is tiny,
  `Evaluate` itself is broken under BP5 convention.  If it's
  physical, the bug is downstream (in the rotation or accumulation).

- **D2:** Instrument `ComputeFaceFluxRHS`'s interior-fault branch to
  dump `Q_imp_plus_g, Q_imp_minus_g, F_h_total` for the hypocenter
  fault QP at one Mult call.  Compare the post-rotation global-frame
  `Q_imp_plus_g[VX]` (= strike component of imposed bulk velocity in
  global) against the expected ~m/s magnitude.  If tiny in global
  frame, the rotation is wrong.

- **D3:** Revert JUST `tpv102_setup.hpp` to pre-R-801 (strike in
  `tau1_0`/V1) and re-run.  If rupture then propagates: the bug is
  in the BP5 convention unification (H1/H3).  If it still doesn't:
  the bug is in R-701/R-802 (H2).

**Suggested fix:** None yet — pending root-cause isolation via D1/D2/D3.

**Test case (to be added regardless of fix):**
```cpp
// tests/parallel/test_r101_shared_fault.cpp — new CRITICAL test:
// On the inline 2-tet mesh or a 4-tet extension, drive the driver's
// RK4 loop for enough steps to reach breakaway on the single fault
// QP at the "hypocenter" (= the one with peak delta_tau).  Then
// query a fault QP 2-4 elements AWAY along strike and assert its V2
// has grown by at least 3 orders of magnitude from V_ini by the end
// of the simulation.
//
// This is the "propagation" test the v7.0.0 suite was missing.
//
// Alternative: assert max‖Q‖_∞ > 1e-3 Pa anywhere in the bulk after
// breakaway.  Physically, a 0.1 m/s slipping fault radiates MPa-
// scale stress perturbations into the adjacent bulk; any v-scale or
// MPa-scale amplitude is qualitatively sufficient to confirm
// coupling.
TEST_R1001_RUPTURE_PROPAGATES_AWAY_FROM_HYPO:
   // ... build 4-tet inline mesh with two fault QPs "along strike"
   // ... initialize one as hypocenter (tau2_0 = tau_ini + peak_dtau)
   // ... run ~1000 steps of the driver's RK4 loop with nucleation
   // Assertion: the NON-hypocenter fault QP's V2 > 1e-9 m/s by end.
```

This test would have caught the v7.0.0 regression locally.  Its
absence is the test-strategy gap.

---

## Summary
- Critical issues: **1** (R-1001 — fault-to-bulk coupling silently
  broken; rupture pinned at hypocenter, does not propagate; root
  cause between R-701/R-801/R-802).
- Moderate issues: **0**.
- Low issues: **0**.
- Plan compliance: **the v7.0.0 fix chain passed its own test
  suite**, but those tests are insufficient to detect the observed
  failure mode.  Test-strategy gap, not a plan deviation.
- **Verdict: FAIL — do NOT proceed with v7.1.0 or any further
  Frontera runs until R-1001 root cause is isolated and fixed.**

## What the v7.1.0 plan should now become

- **Reject v7.1.0 R-901 RESULT.txt regex fix.**  The FAIL verdict
  is correct in outcome even if its stated reason (r399 dead)
  misleads.  Fixing the regex would silence a correct alarm.
- **Preserve the existing RESULT.txt regex** as-is until R-1001 is
  fixed; a false-positive FAIL is preferable to a false-positive
  PASS.

## Version bump
- Pre: `v7.0.0`.
- Found: 1 CRITICAL, 0 MODERATE, 0 LOW.  Bumps `[a]`.
- **Post: `v8.0.0`**.  Doc: `tpv102_debug_v8.0.0_check.md`.

## Recommended next steps

1. **Apply diagnostic D3 first** (revert `tpv102_setup.hpp`
   component assignment to pre-R-801 values) — lowest-effort
   probe.  30 minutes locally.
2. Based on D3 result, apply D1 or D2 to pinpoint whether the
   failure is in `Evaluate`, rotation, or accumulation.
3. Add the R-1001 propagation unit test to
   `test_r101_shared_fault.cpp`, failing under current v7.0.0 code,
   to have a local fail-case before attempting a fix.
4. Fix R-1001 with the minimum-change approach informed by D1/D2.
5. Re-run the full local test suite + a local 4-rank 1000m smoke +
   re-submit the 400-rank dispositive.  Success criteria: at least
   3 of 9 off-hypocenter fault stations show `V2 > 1e-3 m/s` by
   `tfinal = 1.5 s`, and `max‖Q‖_∞ > 1e3 Pa` anywhere.

## Postmortem on v7.0.0 audit

The v7.0.0 audit I wrote said "no blockers remaining, code is
Frontera-ready."  That claim was based on:
- All tests green locally.
- A 4-rank 1000m `tfinal=0.1` smoke run completing without abort.

Both signals were technically correct but insufficient:
- The tests don't probe propagation — they probe cross-rank
  consistency.
- A `tfinal=0.1` run doesn't reach breakaway, so any post-breakaway
  bug is invisible.

**Lesson:** for any dynamic-rupture refactor, before declaring
"Frontera-ready" require a **local end-to-end test that runs past
breakaway and asserts off-hypocenter slip**.  Even at 1000m
resolution, 4 ranks × 1.5 s would produce enough signal to detect
whether bulk waves carry stress to neighboring QPs.  The v7.0.0
audit did not include such a test; adding R-1001 closes this gap
retroactively.

## Unreviewed Areas
- Whether the same bug exists for BP5 (quasi-dynamic, elliptic solve).
  R-801 Option A unified BP5 and TPV102 on the same convention; if
  BP5's tests continue to pass (461/461 locally), BP5 is either
  self-consistent under the convention or its tests don't probe the
  failure mode either.  Worth a separate audit.
- Whether the PRE-R-801 code (v6 or earlier) would have produced a
  propagating rupture on the 400-rank dispositive.  We know from
  the Frontera job 7664983 that v7.0.0-pre-R-801 aborted at the
  R-101 verifier; we never saw that code run to completion.  So we
  don't have a "known-good" reference.  Tandem on the same mesh
  would be the cleanest external reference.
