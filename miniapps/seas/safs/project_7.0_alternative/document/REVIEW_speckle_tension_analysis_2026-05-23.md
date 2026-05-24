# Review of the tensile-σ_n speckle analysis — 2026-05-23

## Scope
Adversarial audit of the agent analysis that concludes the tensile `σ_n` runaway is "the
iterator's dynamic Godunov normal traction from a **corrupted shared-QP predictor Q**" and proposes
to next read `EvaluateBulkAtFaultQPsCanonical` for a stale/mis-rotated ghost `Q⁻`.

**Code claims verified (the analysis is accurate here):**
- `spatial_setup.hpp:86` `d.sigma_n0 = sigma_n_eff`; `:88` `sigma_n_nuc = 0`. ✓
- `s.sigma_n_total = sigma_n0 + sigma_n_nuc + s.sigma_n_trial` (`StepOneQP_:88`,
  `EvaluateADER_LSW:785`). ✓
- `σ_n_trial = η_p·(Q⁻[VX] − Q⁺[VX] + Q⁺[SXX]/Z_p⁺ + Q⁻[SXX]/Z_p⁻)` (`fault_face_flux.cpp:62-64`). ✓
- friction never changes `σ_n` (`σ_n_corr = σ_n_trial`, `:813` / `tpv205_friction.hpp` comment). ✓
- η_p ≈ 0.5·ρ·c_p ≈ 8×10⁶, η_p/Z_p = 0.5. ✓

**What I agree with:** the *mechanism* — the tension lives entirely in `σ_n_trial` (the bulk
Godunov normal traction), not in friction/`σ_n0`/`σ_n_nuc`; and it is a **positive feedback**
(V↑ → radiated normal traction more tensile → `τ_str→0` → V↑). Both are correct and important.

**What I do NOT agree with:** the *diagnosis* (corrupted ghost `Q`) and the *next step* (hunt the
ghost path). The conclusion overreaches the evidence and, more importantly, **excludes the most
likely cause** — the physical, uncapped dip→normal feedback (plan B.2 fact 3 / R-005) — by framing
the question as a binary "ghost-Q vs frame-error." Findings below.

---

## Findings

### [R-001] CRITICAL [analysis diagnosis] — "corrupted ghost Q" contradicts the Phase-0 finding (ghost DOFs bit-exact) and the physical-feedback evidence; the most likely cause (uncapped normal-channel feedback) is excluded by a false binary

**Category:** ASSUMPTION (wrong root-cause attribution)

**Description:**
The analysis frames the cause as ghost-`Q` corruption "rather than a frame error" — a binary that
omits the third, best-supported option: **the −43 MPa seed `σ_n_trial` is a genuine dynamic
normal-stress change (dip slip on the 32° dipping fault radiating into the normal channel), present
on both ranks, amplified to −9 TPa by the uncapped positive feedback** the analysis itself
identified. Three pieces of evidence point here and away from ghost corruption:
1. **Phase 0 already measured the ghost path bit-exact** (job 7747304: "raw +side DOFs are
   bit-identical across the two ranks at every sub-step → rules out ghost-exchange + bulk drift";
   the only cross-rank gap was ~1e-14 *interpolation*). A −43 MPa ghost error is ~15 orders above
   that gap. The analysis's hypothesis contradicts this prior finding and does not reconcile with it.
2. **A 1e-14 interpolation gap cannot produce −43 MPa at the seed.** So the −43 MPa is not a
   cross-rank artifact; both ranks compute ≈−43 MPa (agreeing to 1e-14). That is the signature of a
   *physical* dynamic traction, not a corrupted input.
3. **Plan B.2 predicted exactly this:** on a dipping fault, dip slip couples into the
   **unconstrained** normal channel (no no-opening cap; `σ_n_corr = σ_n_trial`, friction caps only
   shear), "this is why the original blow-up reached `σ_n = −124 GPa` tension." A self-sustaining
   tensile runaway with no cap is the B.2 mechanism, single-rank, identical on both ranks.

**Trigger:** any strongly-slipping QP on the dipping fault inside the rupture area; the radiated
normal traction goes tensile and, uncapped, feeds back.

**Actual behavior (of the analysis):** concludes "corrupted shared-QP predictor Q" and proposes to
audit the ghost path next.

**Expected behavior:** treat the uncapped normal-channel feedback (R-005 / B.2 fact 3) as the
leading hypothesis; the ghost path is a lower-prior branch already weakened by Phase 0.

**Suggested fix (decisive experiment — cheaper and more conclusive than the ghost read):** test the
amplifier directly. Add the no-opening response (zero tangential V when `σ_n_total ≤ 0`) behind an
env gate and rerun:
```cpp
// tpv205_friction.hpp SolveLSW_TPV205, after sigma_n_pos:
if (sigma_n_total <= 0.0) {           // env-gated test of the no-opening hypothesis
   V_abs = 0.0; V1 = 0.0; V2 = 0.0;
   tau1_corr = tau1_trial; tau2_corr = tau2_trial; return;
}
```
If the runaway stops, the cause is the uncapped feedback (regardless of the seed's origin) and the
ghost hunt is moot. If it persists, *then* the seed origin matters and the ghost read is warranted.
(Check the SCEC TPV205 spec / SeisSol convention before keeping this — it changes physics; the
branch is inert for `σ_n>0` so the TPV205 byte-exact regression is preserved.)

**Test case:**
```
# Frontera (or a reduced local mini-SAFS): Dc2 run with SEAS_TENSION_NOOPEN=1.
# Assert: max_slip stays O(10 m), sigma_n_total never < 0 at the worst DOF.
# Compare against the baseline (runaway) to confirm the cap is what bounds it.
```

---

### [R-002] CRITICAL [analysis premise] — `worst_rel = 0` does NOT establish a "single-rank cause"; it is forced by the reconcile, so it cannot distinguish a physical (symmetric) cause from a propagated ghost (asymmetric) one

**Category:** ASSUMPTION (the symmetry inference is invalid)

**Description:**
The chain "`worst_rel = 0` ⇒ not cross-rank ⇒ single-rank cause" (carried from the prior framing
into this analysis) is unsound. The reconcile **forces** the 8 output fields (incl. `sigma_n_corr`,
`slip1/2`, `V1/2`) bit-identical across ranks regardless of whether the *inputs* agree. So
`worst_rel = 0` is guaranteed post-reconcile and is consistent with BOTH (a) a genuinely symmetric
physical cause AND (b) a cross-rank ghost corruption whose boss value is then copied to the peer.
It therefore cannot be used as evidence for either. The `[SLIP]` trace prints only `my_rank`'s
iterator value (−9 TPa), so it likewise does not establish whether the two ranks' iterator
`σ_n_trial` agree.

**Trigger:** inferring causation from the post-reconcile consistency metric.

**Actual behavior:** the analysis treats "both ranks hold the same unphysical state" as proof of a
single-rank cause; it is a reconcile artifact.

**Expected behavior:** to test cross-rank symmetry of the *input*, compare BOTH ranks' iterator
`σ_n_trial` (pre-reconcile) at the same shared QP.

**Suggested fix (instrument both ranks pre-reconcile):**
```cpp
// In the [SLIP] trace, the value is already per-rank (prints my_rank_). Run the
// Dc2 sbatch's existing per-rank split (grep 'rank=N') on the [SLIP] lines for the
// SAME (x,y,z) QP and diff the two ranks' sigma_n_tot/sigma_n_trial at equal t,o:
//   identical (to ~1e-12)  => symmetric/physical (R-001 feedback) — ghost is NOT the cause;
//   differ by O(MPa)        => asymmetric — THEN the ghost path (R-001 branch) is warranted.
```

**Test case:**
```
# Re-split the [SLIP] trace by rank (the sbatch already splits XRANK by rank);
# assert the two owning ranks' sigma_n_trial at the shared QP agree to ~1e-12
# (=> physical) or diverge (=> ghost). This is the test the analysis skipped.
```

---

### [R-003] MODERATE [analysis: the −9 TPa argument] — conflates the runaway RESULT with the seed CAUSE

**Category:** BUG (cause/effect inversion)

**Description:**
The analysis uses the t=1.0 value (`σ_n_trial ≈ −9 TPa`) to argue the iterator's predictor `Q` is
"corrupted." But it earlier (correctly) states this is a positive feedback — so −9 TPa is the
*consequence* of the runaway (slip already radiating garbage into `Q`), not evidence about the
*input* that started it. −9 TPa implies a v-jump of ~1.1×10⁶ m/s (≈1000·c_p) — i.e. the state is
already blown up; reasoning about its provenance is uninformative. Only the **seed** (−43 MPa at
t=0.477) bears on causation, and even there the DOF is already slipping at 12.6 m/s (well into the
event, not the true onset).

**Suggested fix:** restrict all causal claims to the earliest pre-runaway samples; explicitly label
the −9 TPa as downstream of the feedback, not as input evidence.

**Test case:** trace `σ_n_trial` from the *first* sub-step it goes tensile at the worst DOF; assert
the magnitude argument (R-004) is evaluated there, not at t=1.0.

---

### [R-004] MODERATE [analysis: the "43% leak too large for a frame error" magnitude argument] — rests on the v-jump interpretation the analysis admits it cannot confirm

**Category:** ASSUMPTION

**Description:**
The "−5.4 m/s normal jump = 43% of the 12.6 m/s tangential ⇒ too large for a ~25° frame error ⇒
ghost corruption" argument assumes the entire −43 MPa is the **v-jump** term. The analysis itself
notes the log "only prints the summed `σ_n_total`" and cannot split the v-jump term from the
`Q[SXX]` stress term. If the −43 MPa is dominated by the **stress** term (`Q⁺[SXX]/Z_p⁺ +
Q⁻[SXX]/Z_p⁻`), the "43% leak" framing is irrelevant: a −43 MPa dynamic normal-stress perturbation
at a strongly-rupturing point on a dipping fault is physically ordinary (tens of MPa dynamic Δσ_n
is expected from dip-normal coupling), not evidence of corruption. So the headline inference is
contingent on an unmeasured split and should not be stated as a conclusion.

**Suggested fix (the decomposition the analysis correctly wants — do this, but read it neutrally):**
```cpp
// In ComputeTrialTraction or the trace, print the two contributions separately:
const real_t sn_vjump  = data.eta_p * (Q_minus[VX] - Q_plus[VX]);
const real_t sn_stress  = data.eta_p * (Q_plus[SXX]*invZp_plus + Q_minus[SXX]*invZp_minus);
// trace: sigma_n_trial = sn_vjump + sn_stress;  print both.
// stress-dominated => physical dynamic Δσ_n (R-001), NOT a v-jump "leak".
```

**Test case:** assert the decomposed trace shows whether `sn_vjump` or `sn_stress` dominates at the
seed; the analysis's conclusion only holds if `sn_vjump` dominates AND its sign/magnitude is
inconsistent with the local frame tilt.

---

### [R-005] MODERATE [POSSIBLE] [analysis: implicit "the runaway DOF is a shared QP"] — never verified; the entire ghost hypothesis requires it

**Category:** ASSUMPTION

**Description:**
The "iterator vs shared-path" decoupling and the ghost-`Q` hypothesis only exist if the runaway DOF
is a **shared** fault QP (has a ghost neighbor and an inline `ComputeADERSharedFaceFluxRHS` solve).
For an **interior** fault QP there is no ghost and no shared-path solve — the iterator is the only
solve, and "the shared-path stays bounded" has no meaning. The decisive partition overlay
(`vtkProcessId` on the fault) was never run (it was flagged as the open test two rounds ago). If the
worst DOF is interior, the ghost-`Q` diagnosis is categorically inapplicable.

**Suggested fix:** color the fault by `vtkProcessId` / partition in ParaView and check whether the
runaway/speckle DOFs sit on rank seams (shared) or in rank interiors (interior). One render.

**Test case:**
```
# ParaView: load fault.vtkhdf, color by vtkProcessId; overlay the max-slip DOFs.
# Assert seam-coincident (=> shared, ghost hypothesis admissible) vs interior
# (=> ghost hypothesis void; cause is single-rank physics/predictor).
```

---

### [R-006] MODERATE [analysis: "iterator corrupted, shared-path clean"] — the comparison is confounded three ways

**Category:** ASSUMPTION

**Description:**
The inference "iterator `σ_n_trial → −9 TPa` but output stays ~49 MPa ⇒ the iterator's predictor Q
is corrupted, the shared-path's isn't" is confounded:
1. **Different inputs.** The iterator uses the ADER **sub-step predictor** states
   (`Q_pointwise_*` at midpoint nodes); the shared-path inline solve uses the **macro
   time-integrated** `I`. They are expected to differ — and the predictor can overshoot at an
   under-resolved front (R-006/R-007 of the debug review) far more than the macro-integrated state.
   That alone explains a wild iterator value with a bounded macro value — no ghost corruption needed.
2. **The "~49 MPa output" is the RECONCILED value.** `sigma_n_corr` is one of the 8 reconciled
   fields; post-reconcile both ranks show the boss's value, which **masks** whatever the non-boss
   shared-path solve produced. So "the shared-path isn't corrupted" is not established — it is
   hidden by the reconcile.
3. **Same-DOF not confirmed.** Whether the −9 TPa trace and the ~49 MPa "output" are read at the
   *same* DOF is unstated; if not, the decoupling compares two different points.

**Suggested fix:** compare, at one identified shared QP and one sub-step, the iterator's
`σ_n_trial` against the shared-path inline solve's `σ_n_trial` **pre-reconcile** (not the reconciled
output), and separately against the macro-integrated `I` decomposition — to separate "predictor
overshoot" from "ghost corruption."

**Test case:**
```
# One shared QP, one macro step: log (a) iterator sigma_n_trial (substep predictor),
# (b) inline-solve sigma_n_trial PRE-reconcile, (c) the reconciled output.
# If (a) wild but (b) bounded with the SAME ghost Q => predictor overshoot, not ghost.
# If (a)≈(b) both wild => not a substep-vs-macro artifact.
```

---

## Summary
- Critical issues (in the analysis): 2 — R-001 (ghost-Q diagnosis excludes the better-supported
  uncapped-feedback cause and contradicts Phase-0 bit-exact-ghost), R-002 (`worst_rel=0` does not
  imply single-rank; it's a reconcile artifact).
- Moderate issues: 4 — R-003 (cause/effect inversion of the −9 TPa), R-004 (magnitude argument
  rests on an unconfirmed v-jump split), R-005 (unverified shared-vs-interior), R-006 (confounded
  iterator-vs-shared comparison).
- **Verdict: PARTIAL agreement.** The mechanism (tension in `σ_n_trial`; positive feedback) is
  correct and well-verified. The **diagnosis (corrupted ghost Q) is not supported** and is likely
  the wrong target: the evidence (physical −43 MPa seed consistent across ranks per Phase 0; an
  uncapped normal channel; a dipping fault) points to the **B.2 dip→normal uncapped feedback**
  (R-005 of the debug review), a single-rank physics bug — not a ghost-exchange defect.
- **Recommended next step (replaces "go read the ghost path"):** run the two cheap, decisive tests
  first — (1) the no-opening-cap experiment (R-001 fix: does capping `V` under tension stop the
  runaway?) and (2) the `σ_n_trial` v-jump/stress decomposition + the two-rank pre-reconcile diff
  (R-002/R-004/R-006) + the `vtkProcessId` overlay (R-005). These either confirm the uncapped-
  feedback cause (most likely) or, only if the runaway is genuinely cross-rank/asymmetric, justify
  the ghost-path read the analysis proposed. Reading `EvaluateBulkAtFaultQPsCanonical` first risks
  a long hunt for a defect Phase 0 already found absent.

## Unreviewed Areas
- The actual `[SLIP]` log values (`+6 MPa`, `−9.4 TPa`, `12.6 m/s`) are taken from the analysis's
  quotation; not independently re-read from the run logs.
- `EvaluateBulkAtFaultQPsCanonical` ghost/self path was NOT read here (the analysis's proposed
  target) — deliberately, because R-001/R-002 argue it is the wrong first read; revisit only if the
  two-rank pre-reconcile diff (R-002) shows genuine asymmetry.
- Whether `Q[SXX]` in the canonical frame is tension- or compression-positive at the bulk→fault
  interface was not re-derived; the friction path treats `σ_n>0` as compression (`max(σ_n,0)`),
  which is internally consistent, but the bulk-stress sign convention feeding `σ_n_trial` should be
  confirmed if the decomposition (R-004) shows the stress term dominating.
