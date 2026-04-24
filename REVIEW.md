# Code Review: 2026-04-23 Round 10 — §K Round-9 Frozen-Friction + dt-Refinement Verdict

## Review Scope
- Plan: `phase1_arm1_findings_2026-04-23.md` §K (round-9 deliverable, three-option verdict).
- New artifacts:
  - `phase2e_r9_freezeC_kuhn.txt`, `phase2e_r9_freezeC_d4.txt` (raw logs).
- Domain context: REVIEW round 9 R-001..R-006, prior rounds 1-8,
  CLAUDE.md (Tandem reference), MEMORY rules.

## Bottom-line up front

**§K is the most epistemically honest deliverable in 9 rounds.**  The
implementer correctly admits no single-line localization has been
found, identifies the seed-hunting pattern as exhausted, and presents
three options for the user to decide.

**But the §K options are NOT equivalent in information gain.**  Two
critical empirical anomalies in §K's data + one option-mislabeling
mean the user, if presented the three options as currently framed,
would likely pick a path that doesn't disambiguate the underlying
question.

The recommendation is straightforward: **Option X (Tandem benchmark)
first**, defer Y/Z.  Option X is the only one that can disprove the
"diffuse-noise floor" hypothesis — and disproving it is the only way
to redirect the investigation toward an actual fix.  Y can wait;
Z (accept the floor + adjust gate) is Direction D in disguise, which
the user already rejected.

---

## Findings (Round 10)

### [R-001] [CRITICAL] [phase2e_r9_freezeC_*.txt / dt-exponent -0.17 anomaly] — Pepper WORSENS with smaller dt; this fits neither "sub-resolved physics" nor "per-step injection" cleanly; needs explanation before being used as evidence

**Category:** ASSUMPTION (data interpretation may be wrong)

**Description:**
§K reports dt-refinement result:
> dt exponent = -0.17 on BOTH fixtures → NOT sub-resolved physics.

Pepper exponent -0.17 means `tau1_corr ∝ dt^(-0.17)`.  As dt → 0
(more refined), pepper INCREASES weakly.  The implementer correctly
notes this rules out sub-resolved physics (which would give positive
exponent — pepper decreases as dt → 0, since the under-resolved
instability frequency scales with dt).

But the negative exponent is unusual and the implementer doesn't
explain WHY pepper grows with refinement.  Two distinct mechanisms
predict negative exponents, with very different interpretations:

(M1) **Per-step ULP injection accumulating linearly in N=t/dt.**  If
the source injects ULP-magnitude noise per step regardless of dt,
total noise after fixed simulation time `t` scales as N·ULP = t/dt·ULP
→ pepper ∝ dt^(-1).  Pure per-step injection gives exponent -1.

(M2) **Numerical drift in the friction-solver tolerance.**  Brent's
method tolerance is typically a fixed value (~1e-8); per-iteration
work is the same regardless of dt.  Smaller dt = more friction calls
(N times), so cumulative drift scales as N · drift_per_call.  Same
exponent -1.

(M3) **Mixed mechanism.**  Per-step injection (M1, exponent -1)
combined with amplitude saturation at large N (exponent +1) gives a
net intermediate exponent.  -0.17 is consistent with such a mix where
the saturation barely beats the injection.

(M4) **Dimensional artifact.**  The dt-refinement test runs to fixed
SIMULATION TIME or fixed STEP COUNT?  The two give different scaling.
At fixed time t = 20·dt_baseline, halving dt doubles N but keeps t
constant.  Pepper ∝ dt^(-0.17) at fixed t.  At fixed N=20 steps
(varying simulation time t = 20·dt), pepper would have a very
different scaling.  The implementer doesn't specify which.

The -0.17 exponent is INFORMATIVE but incompletely interpreted:
- It does rule out "sub-resolved physics" (would need positive
  exponent).
- It does NOT confirm "code bug outside friction and outside ψ"
  cleanly.  M1, M2, and M3 all predict the observed sign without
  being a code bug per se.
- M2 in particular suggests the friction solver IS contributing
  per-call drift, contradicting §K's matrix-branch verdict.

**Trigger:** The user picks Option Y (authorize FREEZE-A/B) on the
basis of "amplifier ruled out, must be friction mechanism (FREEZE-A/B
domain)"; FREEZE-A/B tests then return mixed verdicts because the
underlying mechanism is per-step injection scaling with N (M1/M2),
not a "single-line friction bug".

**Actual behavior:** §K reports `-0.17 → NOT sub-resolved physics`
and stops there.  The negative-but-not-near-(-1) exponent is not
explained.

**Expected behavior:** Disambiguate M1 vs M2 vs M3 vs M4 before
relying on the dt-scaling result as Phase 2C-fault evidence.  Two
specific tests:

1. **Verify dt-scaling test design.**  Document whether dt-refinement
   runs to fixed t or fixed N.  If fixed-N, repeat with fixed-t (or
   vice versa) to identify which mechanism is sensitive to which
   parameter.
2. **Test friction-call count vs pepper.**  Run the pepper guard at
   N steps but with the friction solver called once per step (current)
   vs once every k steps (skip friction on intermediate steps).  If
   pepper scales with friction-call count, M2 (friction-solver per-call
   drift) is confirmed.  If it doesn't, M1 (per-step injection from
   another source) stands.

**Suggested fix (documentation only):**
```diff
@@ phase1_arm1_findings_2026-04-23.md §K dt-refinement section
- dt exponent = -0.17 on BOTH fixtures → NOT sub-resolved physics.
+ dt exponent = -0.17 on BOTH fixtures.
+
+ Sign analysis:
+ - Positive exponent (pepper drops as dt → 0): would indicate
+   sub-resolved physics.  RULED OUT.
+ - Exponent near -1: would indicate per-step injection (per-call
+   ULP noise accumulating linearly in N).
+ - Observed exponent -0.17: between zero and -1.  Compatible with
+   either (a) per-step injection partially saturated at large N,
+   (b) mixed injection + amplification mechanism, or (c) test-design
+   artifact (fixed-time vs fixed-step scaling).
+
+ Open: which scaling protocol was used (fixed t = 20·dt_baseline,
+ varying N; vs fixed N=20, varying t)?  Re-run with the alternate
+ protocol to disambiguate.
+
+ Open: friction-call-count test (run pepper guard with friction
+ called every k steps) — distinguishes per-call drift (friction
+ contributes ULP each invocation) from per-step injection (some
+ other source contributes per step regardless of friction calls).
```

**Test case:**
```cpp
// tests/unit/test_R001_round10_dt_friction_call_count.cpp
TEST(R001Round10, PepperScalesWithFrictionCallCountNotJustDt) {
   for (int friction_skip : {1, 2, 4}) {  // call friction every k steps
      auto pepper = RunPepperGuardWithFrictionSkip("kuhn",
                                                     friction_skip);
      std::cout << "friction skip=" << friction_skip
                << ": tau1_corr=" << pepper.tau1_corr << "\n";
   }
   // If pepper drops by ~k× when friction is called 1/k as often,
   // friction is contributing per-call drift (M2).  If invariant,
   // M1 stands.
}
```

---

### [R-002] [CRITICAL] [phase2e_r9_freezeC_*.txt / FREEZE-C ratio = 1.0000 exactly] — Both fixtures report ratio = exactly 1.0000; either ψ has no effect (genuine but surprising) OR the FREEZE-C implementation is a no-op (test bug)

**Category:** ASSUMPTION (verification needed)

**Description:**
§K reports:
> FREEZE-C ratio = 1.0000 on BOTH fixtures → ψ rate-state feedback
> NOT amplifier.

Ratio = exactly 1.0000 (4 decimal places, no perturbation) on BOTH
fixtures is suspicious.  If FREEZE-C made ANY change to the simulation
state, the result should differ by SOMETHING — even at machine
precision.  The exact-1.0000 result has three explanations:

(A) **Genuine: ψ literally doesn't change between steps in this 20-step
    test, so freezing it is a no-op.**  Possible if the rupture drive
    is short enough that ψ updates by < 1 ULP per step.  But: under
    `tau2_nuc = nuc_dtau` rupture drive over 20 steps with dt=5e-5,
    `psi` should evolve measurably (rate-state evolution at high V is
    fast).  Verify by printing ψ before and after the 20-step run.

(B) **Test bug: the FREEZE-C implementation is a no-op.**  If the
    implementation didn't actually freeze ψ (e.g., wrong field hooked,
    or freeze was applied to a copy that gets reset), the test ran
    standard physics and reported the standard result vs itself →
    ratio = 1.0000.  The "verdict: ψ NOT amplifier" would then be
    based on a non-test.

(C) **Pepper is fully determined by quantities ψ doesn't affect.**
    Genuine result: ψ enters the friction solve as
    `theta = exp(psi/a)/(2V0)`, but if Theta dominates (high V,
    rupture phase), the V_abs solution is essentially independent of
    `theta` for this regime.  Possible at high V where the
    `eta_s·V` term dominates `strength·f(V)`.

The implementer's verification report doesn't show the ψ values to
prove (A) or rule out (B).  An exact 1.0000 result calls for direct
verification.

**Trigger:** The user picks Option Y or Z based on the "ψ NOT
amplifier" verdict; the verdict is actually a no-op test; round 10
or 11 discovers ψ IS the amplifier and the diagnostic was broken.

**Actual behavior:** Ratio = 1.0000 reported as evidence; no
verification of whether ψ values actually differed.

**Expected behavior:** Print ψ at step 0 and step 20 in BOTH baseline
AND FREEZE-C runs.  If baseline ψ at step 20 differs from step 0 (ψ
DID evolve), and FREEZE-C ψ at step 20 EQUALS step 0 (freeze took
effect), the result is genuine.  If baseline ψ doesn't evolve, the
test is too short.  If FREEZE-C ψ ALSO evolved, the freeze didn't
take.

**Suggested fix:**
Add to the test:
```diff
@@ tests/unit/test_r9_frozen_friction_C.cpp main()
+   // R-002 round 10 verification: confirm FREEZE-C actually froze ψ.
+   // If baseline ψ_step20 == ψ_step0 (no evolution), test is too short.
+   // If FREEZE-C ψ_step20 != ψ_step0 (freeze did NOT take), the test
+   // is broken and the "ψ not amplifier" verdict is unsupported.
+   std::cout << "  baseline ψ at step 0:  " << dof_data_baseline[0].psi << "\n";
+   std::cout << "  baseline ψ at step 20: " << dof_data_baseline[0].psi_after_20 << "\n";
+   std::cout << "  FREEZE-C ψ at step 0:  " << dof_data_freezeC[0].psi << "\n";
+   std::cout << "  FREEZE-C ψ at step 20: " << dof_data_freezeC[0].psi_after_20 << "\n";
+   ASSERT_NE(dof_data_baseline[0].psi_after_20, dof_data_baseline[0].psi)
+      << "Baseline ψ did not evolve — test is too short to exercise ψ feedback.";
+   ASSERT_EQ(dof_data_freezeC[0].psi_after_20, dof_data_freezeC[0].psi)
+      << "FREEZE-C ψ did evolve — freeze did NOT take effect.";
```

**Test case:** The verification print itself.

---

### [R-003] [CRITICAL] [phase1_arm1_findings_2026-04-23.md §K Option Y framing] — Option Y is described as "authorize FREEZE-A/B flux-layer modifications", but FREEZE-A/B are TEST-LEVEL interventions; they do not require flux-layer freeze unblock

**Category:** BUG (option mislabeled, may falsely block diagnostic)

**Description:**
§K presents Option Y:
> Option Y (freeze unblock required): authorize FREEZE-A/B flux-layer
> modifications.

Per round-9 R-003's specification, FREEZE-A and FREEZE-B are
implemented as test-level wrappers around the friction-solve call
site:
- FREEZE-A: skip the `fault_flux_->Evaluate(...)` call in
  `tpv102_driver.cpp`'s RK4 loop; substitute `tau*_corr = tau*_trial,
  V = 0`.
- FREEZE-B: cache the step-1 `dof_data` after one Evaluate call;
  reuse cached values at subsequent steps without calling Evaluate.

Neither modifies `FaultFaceFlux::Evaluate`, `wave_operator.inl`'s
fault dispatch, `PrecomputedFaceFluxes`, or any production flux-layer
code.  Both modify ONLY:
- `tpv102_driver.cpp` (test-driver scope; not production), AND/OR
- A new test binary that calls the production code differently.

The flux-layer freeze (per §C R-006) covers production code changes
to `wave_operator.inl`, `PrecomputedFaceFluxes`, `FaultFaceFlux`, etc.
It does NOT cover test-level wrappers around production functions.

By mislabeling Option Y as requiring freeze unblock, §K artificially
blocks a freeze-allowed diagnostic.  The user might decline Option Y
on the basis "I don't want to authorize flux-layer changes" when in
fact Option Y is identical in scope to the round-9 FREEZE-C work that
was already done.

**Trigger:** User reads "freeze unblock required" and declines Option
Y; instead picks Option X or Z; the next-most-informative diagnostic
is skipped on a process technicality.

**Actual behavior:** Option Y mislabeled.  Round-9 already implemented
FREEZE-C without unblock authorization; FREEZE-A/B are no different in
scope.

**Expected behavior:** Re-classify Option Y as freeze-allowed, same
as Options X and Z.

**Suggested fix:**
```diff
@@ phase1_arm1_findings_2026-04-23.md §K options
- Option Y (freeze unblock required): authorize FREEZE-A/B flux-layer
- modifications.
+ Option Y (freeze-allowed): run FREEZE-A and FREEZE-B variants.
+ Both are test-level wrappers around the friction-solve call site
+ (analogous in scope to round-9's FREEZE-C); neither modifies
+ production flux-layer code.  The flux-layer freeze (per §C R-006)
+ applies to production code changes; test-level wrappers around
+ production-code call sites are diagnostic and freeze-allowed.
+
+ User authorization is only required if FREEZE-A/B verdicts INDICATE
+ a fix to FaultFaceFlux source code.  Until then, the diagnostic
+ runs are freeze-allowed.
```

**Test case:** N/A (process clarification).

---

### [R-004] [MODERATE] [phase1_arm1_findings_2026-04-23.md §K Option X / Y / Z framing] — Three options presented as user choice; they're not equivalent in information gain; Option X (Tandem) should be the strong recommendation

**Category:** DEVIATION (decision framing)

**Description:**
§K presents three options as a user-decision menu.  The framing
treats them as roughly equivalent paths.  They are not:

| Option | Informs | Decisive? | Info gain |
|---|---|---|---|
| X (Tandem benchmark) | Whether ANY MFEM-specific bug exists at this resolution | YES — Tandem clean → bug exists; Tandem dirty → diffuse-noise floor confirmed | HIGH |
| Y (FREEZE-A/B) | Whether friction MECHANISM contributes (vs ψ which is ruled out) | Partially — verdicts could be ambiguous (per round-9 R-001 mixed-regime concern) | MODERATE |
| Z (accept floor + adjust gate) | Nothing new | NO — commits to "accept" without further investigation | ZERO |

Option X has the unique property of being able to **disprove** the
"diffuse-noise floor / no localizable bug" hypothesis.  If Tandem at
the same 2×2×2 fixture with same ψ-update parameters and same dt
shows tau1_corr near ULP, then:
- The diffuse-noise hypothesis is FALSE.
- A localizable bug DOES exist in MFEM SEAS specifically.
- The next investigation step is to identify where MFEM and Tandem
  diverge.

If Tandem shows ~Pa-scale pepper at the same fixture, then:
- The diffuse-noise hypothesis is supported.
- The pepper is a discretization-floor artifact common to both
  codes.
- Direction D applies (with cross-code validation).

Option X is the **only** option that can falsify §K's diffuse-noise
hypothesis.  Option Y at best refines the candidate space within
the assumption that a bug exists.  Option Z assumes no bug exists.

§K should strongly recommend X.  Y can run in parallel if compute
budget allows.  Z is Direction D — already user-rejected.

**Trigger:** User picks Option Y or Z without running X first; round
10 spends its budget refining within an unverified hypothesis.

**Suggested fix:**
```diff
@@ phase1_arm1_findings_2026-04-23.md §K options
- Verdict: §K ends with THREE OPTIONS requiring user decision:
- - Option X (freeze-allowed, recommended): benchmark against Tandem at same fixture resolution.
- - Option Y (freeze unblock required): authorize FREEZE-A/B flux-layer modifications.
- - Option Z (freeze-allowed, pragmatic): accept the 2–3 Pa pepper floor at 2×2×2 fixture and adjust unit-test gate tolerance.
+ Verdict: §K's diffuse-noise hypothesis is consistent with current
+ data but is NOT proven.  The three options have very different
+ information-gain profiles:
+
+ STRONGLY RECOMMENDED: Option X (Tandem benchmark, freeze-allowed,
+ ~1 day).  Tandem at the same 2×2×2 / dt / ψ params is the only
+ test that can disprove the diffuse-noise hypothesis.
+ - If Tandem clean (tau1_corr ≤ ULP) → bug exists in MFEM SEAS;
+   round 10 targets the MFEM-vs-Tandem code-path diff.
+ - If Tandem dirty (tau1_corr ~ Pa) → diffuse-noise hypothesis
+   confirmed by independent reference implementation.  Consider
+   Direction D (which user previously rejected, but with new
+   cross-code evidence).
+
+ DEFERRED: Option Y (FREEZE-A/B variants).  Useful if Option X
+ confirms an MFEM-specific bug; redundant if Option X confirms
+ diffuse-noise.  Run after X.
+
+ NOT RECOMMENDED: Option Z (accept floor + adjust gate).  This is
+ Direction D, which the user explicitly rejected.  Re-presenting it
+ as "pragmatic" without new evidence does not change the user's
+ prior decision.
```

**Test case:** N/A (recommendation framing).

---

### [R-005] [MODERATE] [phase1_arm1_findings_2026-04-23.md §K diffuse-noise hypothesis] — "No uniquely localizable bug exists" is unfalsifiable as currently stated; Option X is the only proposed test that can falsify it

**Category:** ASSUMPTION (epistemic framing)

**Description:**
§K's honest assessment:
> the pepper signature is consistent with DIFFUSE per-step numerical
> noise (NOT a single-line bug).  Rounds 5–9 have each named a
> target rescinded by the next round — pattern suggests no uniquely
> localizable bug exists.

This is reasonable inductive reasoning from 5 rounds of failed
localization.  But "no uniquely localizable bug exists" is not
falsifiable in the strict Popperian sense — absence of evidence is
not evidence of absence.  The five rounds show "we couldn't find one
with these methods", not "one doesn't exist".

Possible alternative hypotheses §K doesn't address:
- (H1) The bug IS uniquely localizable but at a layer Arm 1 / Phase
  2D static + amplifier diagnostics don't probe (e.g., MFEM internal
  element-local mass-matrix lumping, FE-collection orientation
  tables, MPI ghost-buffer state).
- (H2) The bug is in a code path that's BIT-EXACT under the test
  fixtures (which is why the diagnostics never see it) but
  ULP-noisy under production-scale conditions.  Not measurable on
  the 2×2×2 fixture but real at TPV102 scale.
- (H3) The bug is a missing term (not an incorrect term) — e.g., a
  consistent-mass-matrix correction that should be applied but isn't.
  Diagnostics that compare existing terms can't see a missing term.

Option X (Tandem benchmark) is the only proposed test that can
disprove the diffuse-noise hypothesis without needing to enumerate
H1/H2/H3.  Tandem at the same fixture with bit-exact-equivalent
inputs would either show pepper (confirms diffuse-noise) or not
(falsifies it).

**Suggested fix:**
```diff
@@ phase1_arm1_findings_2026-04-23.md §K honest assessment
- Honest assessment: the pepper signature is consistent with DIFFUSE
- per-step numerical noise (NOT a single-line bug).  Rounds 5–9 have
- each named a target rescinded by the next round — pattern suggests
- no uniquely localizable bug exists.
+ Honest assessment: the pepper signature is consistent with DIFFUSE
+ per-step numerical noise.  This hypothesis is supported by 5 rounds
+ of failed localization (rounds 5-9 each named a target rescinded
+ by the next round) but is NOT proven.
+
+ Alternative hypotheses not yet ruled out:
+ - The bug is at a code layer Arm 1 / Phase 2D doesn't probe
+   (mass-matrix lumping, FE-collection orientation tables, MPI
+   ghost state).
+ - The bug is bit-exact-clean on the 2×2×2 fixture but ULP-noisy
+   at production scale.
+ - The bug is a MISSING term (consistent-mass correction, etc.) —
+   diagnostics that compare existing terms cannot see this.
+
+ Option X (Tandem benchmark) is the unique test that can falsify
+ the diffuse-noise hypothesis without enumerating these
+ alternatives.  Recommend X first.
```

---

### [R-006] [LOW] [phase1_arm1_findings_2026-04-23.md §K user-decision framing] — After 5 rounds of failed seed-hunting, asking the user to choose without a clear recommendation pushes ownership of an investigative decision to a non-investigator

**Category:** QUALITY (process)

**Description:**
§K's closing:
> Please indicate which option to pursue for round 10: X (Tandem
> benchmark), Y (authorize FREEZE-A/B), or Z (accept floor + adjust
> gate).

After 5 rounds of investigation, the implementer has the deepest
context for which option will yield the most information.  Asking
the user to choose without a recommendation:
- Defers an investigative decision to a non-investigator.
- Spreads decision overhead across the user's time budget.
- Invites the user to pick based on time/cost considerations rather
  than information-gain.

A better framing: implementer recommends Option X with reasoning;
user can override if Tandem isn't installed/accessible OR if the
user wants to commit to Y/Z directly.

Per CLAUDE.md, Tandem is a primary reference for SEAS implementation
(`Always refer to Tandem code at /Users/chunhuizhao/projects/tandem`).
The Tandem benchmark is institutionally established as the validation
path; recommending Option X aligns with project guidance.

**Suggested fix:**
```diff
@@ phase1_arm1_findings_2026-04-23.md §K closing
- Please indicate which option to pursue for round 10: X (Tandem
- benchmark), Y (authorize FREEZE-A/B), or Z (accept floor + adjust
- gate).
+ Recommended: Option X (Tandem benchmark).  Tandem is the project's
+ established SEAS reference (per CLAUDE.md), is locally available
+ at /Users/chunhuizhao/projects/tandem, and Option X is the unique
+ test that can falsify the diffuse-noise hypothesis.  ~1 day cost.
+
+ Override paths (if user prefers not Option X):
+ - Y (FREEZE-A/B): run as supplemental diagnostic; redundant if X
+   resolves but useful as parallel evidence.
+ - Z (accept floor): user already rejected Direction D; not
+   recommended unless cross-code validation (Option X) confirms
+   floor-is-floor.
+
+ Awaiting user confirmation to proceed with Option X.
```

---

## Summary
- Critical issues: **3** (R-001 dt exponent unexplained; R-002
  FREEZE-C ratio = exactly 1.0000 needs verification; R-003 Option
  Y mislabeled as flux-layer-modification).
- Moderate issues: **2** (R-004 three options not equivalent in
  info gain; R-005 diffuse-noise hypothesis unfalsifiable without X).
- Low issues: **1** (R-006 user-decision framing).
- Plan compliance: **PARTIAL** — round-9 R-001..R-006 verifications
  closed; new findings are about how §K interprets and presents
  the round-9 verdict.
- Verdict: **PASS WITH FIXES — §K should be amended per R-002 (verify
  FREEZE-C took effect) and R-003 (re-classify Option Y as freeze-
  allowed) before the user is asked to choose.  Strong recommendation
  per R-004: Option X first.**

## Unreviewed Areas
- The `test_r9_frozen_friction_C.cpp` source itself — only the verdict
  was inspected; the FREEZE-C implementation (where ψ is reset) was
  not directly read.  R-002 verification depends on this.
- The dt-refinement test design (fixed-time vs fixed-step protocol)
  — needed to apply R-001's M1/M2/M3/M4 disambiguation.
- Tandem's TPV102 setup — would benefit from a quick read of
  `/Users/chunhuizhao/projects/tandem` to confirm the 2×2×2 fixture
  is reproducible there before Option X is selected.

---

## Suggested Next Step

The investigation has reached a real decision point.  The implementer
correctly admits no single-line localization has been found and
presents three options.  But the options are mis-equivalent and one
is mis-labeled — the user, presented as-is, would likely make a
suboptimal choice.

### Round 10 Step 0 (~1 hour, freeze-allowed) — fix §K presentation

Apply REVIEW round 10 R-002 (verify FREEZE-C took effect) and R-003
(re-classify Option Y as freeze-allowed) and R-004 (recommend X
strongly).  These are documentation-only fixes; no test re-run
needed.

### Round 10 Step 1 (~1 day, freeze-allowed) — execute Option X

**Tandem benchmark at the same 2×2×2 fixture:**

1. Set up Tandem with the same TPV102 parameters (per CLAUDE.md
   guidance).
2. Configure Tandem to use the SAME mesh as the test fixture
   (`Mesh::MakeCartesian3D(2, 2, 2, TETRAHEDRON, 1000.0, 1000.0,
   1000.0)`).  If Tandem doesn't accept arbitrary meshes, use the
   closest equivalent (typically Tandem uses simplicial meshes from
   Gmsh).
3. Run Tandem for 20 steps with `dt = 5e-5`, uniform `tau2_nuc =
   nuc_dtau` rupture drive.
4. Measure tau1_corr spread across fault QPs at step 19.
5. Compare to MFEM SEAS's 2.236 Pa.

**Verdict matrix:**

| Tandem tau1_corr | MFEM tau1_corr | Conclusion |
|---|---|---|
| ≤ 1e-6 (clean) | 2.24 (dirty) | Bug exists in MFEM specifically; round-11 targets MFEM-vs-Tandem code-path diff |
| ~1 Pa (similar) | 2.24 | Discretization floor common to both implementations; supports diffuse-noise hypothesis with cross-code evidence; Direction D becomes defensible |
| ~10× MFEM (much worse) | 2.24 | Tandem has a different bug; not informative for MFEM |
| Cannot run on equivalent fixture | — | Option X infeasible; fall back to Option Y |

### Round 10 Step 2 (depends on Step 1 verdict)

**If Tandem clean (MFEM-specific bug):**
- Round 11 inspects the code-path diff.  Tandem uses (per CLAUDE.md
  references) the same Pelties Eq. 7 trial-traction formula and same
  Brent friction solve.  Diff is in the DG operator (face flux,
  basis), the time integrator, OR auxiliary infrastructure.
- Use git-diff-style enumeration of the differences.
- Pick the highest-leverage diff for round 11 instrumentation.

**If Tandem dirty (cross-code floor):**
- Direction D becomes defensible with cross-code evidence.  User can
  re-evaluate the prior rejection.
- Alternative: discretization refinement.  Run pepper guard at
  4×4×4 fixture (8× refinement); check if pepper drops below the
  v9.4.0 §11 1e-10 gate.  If yes, the production mesh (which is
  much finer than 2×2×2) may be already below floor.

**If Tandem infeasible:**
- Fall back to Option Y (FREEZE-A/B), but only after R-002 verifies
  the round-9 FREEZE-C wasn't a no-op.

### What NOT to do

- Do not present §K's three options to the user without applying
  R-002, R-003, R-004 fixes first.  The user, presented as-is, may
  decline Option Y on the false-positive "freeze unblock required"
  label, or pick Option Z without realizing it's Direction D.
- Do not commit to Direction D without cross-code evidence.  The
  user's prior rejection of Direction D was made without Tandem
  data; with cross-code confirmation the rejection might be revised,
  but it should be the user's revised decision, not an implementer
  re-presentation of a previously-rejected option.
- Do not run round 11 on a layer-level localization.  After 5 rounds
  of failed seed-naming, the next round must produce either a
  cross-code verdict (X), a clean FREEZE-A/B disambiguation (Y), or
  an explicit Direction-D decision (Z).  No more "we localized to a
  layer".
