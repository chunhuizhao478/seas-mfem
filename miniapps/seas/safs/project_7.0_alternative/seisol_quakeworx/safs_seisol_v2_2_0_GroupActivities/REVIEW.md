# Code Review: GROUP_ACTIVITY_PLAN.md (rate-and-state SAFS activity) — 2026-06-16

## Review Scope
- Plan: `safs_seisol_v2_2_0_GroupActivities/GROUP_ACTIVITY_PLAN.md`
- Files reviewed: the SeisSol decks/yaml embedded in the plan (`safs_fault_rs.yaml`,
  `parameters_aging.par`/`parameters_slip.par`), the size budget, and the ParaView recipes.
- Domain context consulted: SeisSol source `v1.3.1-1760-g49bdd63e4`
  (`DRParameters.{h,cpp}`, `Factory.cpp`, `RateAndStateInitializer.cpp`, `BaseDRInitializer.cpp`,
  `FaultWriter.cpp`, `EnergyOutput.cpp`, `ParameterReader.h`); existing repo decks
  (`tpv104/`, `safs_seisol_v2_1_0_RSSRW/`); project `CLAUDE.md`.

Note: the "code" here is a plan whose embedded SeisSol decks are intended to be copied verbatim
into runnable files, so the review targets deck correctness and the load-bearing numeric claims.
Test cases are given as shell/verification checks (the artifact is not Python).

## Findings

### [R-001] [MODERATE] parameters_aging.par/parameters_slip.par (plan line 207) — `GPwise = 1` is not a parameter in the cited SeisSol version

**Category:** DEVIATION / BUG

**Description:**
The decks set `GPwise = 1` in `&DynamicRupture` and the plan claims it "resolve[s] a(z)/b(z)
Gauss-point-wise". `GPwise` does not exist anywhere in the cited reference source
(`grep -rin gpwise src/` returns nothing). In this version, fault parameters are ALWAYS evaluated
at fault Gauss points via `FaultGPGenerator` (`BaseDRInitializer.cpp:176`), so the flag is both
invalid and unnecessary. Unknown keys are handled by `warnUnknown` (`ParameterReader.h:141`) — a
warning, not an error — so it is non-fatal, but the plan presents an obsolete key as functional.
The repo's `tpv104/parameters.par` also carries `GPwise`, which is consistent with it being a
legacy key from an older SeisSol.

**Trigger:** Any run using these decks against SeisSol v1.3.1-1760 (the version the plan says it
was "verified against").

**Actual behavior:** SeisSol warns "unknown parameter GPwise" (or silently ignores it); a(z)/b(z)
are evaluated GP-wise regardless. The plan asserts the flag is doing necessary work — false.

**Expected behavior:** No obsolete key; the plan should state GP-wise evaluation is the default.
If the QuakeWorx SeisSol is actually OLDER and reads `GPwise`, that must be stated explicitly
(the plan's precision-unknown caveat already admits the target version is uncertain — this same
uncertainty applies to `GPwise`).

**Suggested fix:**
```diff
- GPwise = 1                                       ! resolve a(z)/b(z) Gauss-point-wise
+ ! (GP-wise fault-parameter evaluation is the DEFAULT in SeisSol >= v1.3; no flag needed.
+ !  If the target QuakeWorx SeisSol is older and requires it, add: GPwise = 1)
```
Also delete the `GPwise = 1` mention from the Phase-2 deck comment block and the inline
"resolve a(z)/b(z) Gauss-point-wise" rationale.

**Test case:**
```bash
# Demonstrates the key is not recognized by the cited source:
test ! "$(grep -rin 'gpwise' /Users/chunhuizhao/projects/SeisSol/src/)" \
  && echo "PASS: GPwise absent from source -> obsolete/no-op"
# And confirm GP evaluation is the default path:
grep -q "FaultGPGenerator" /Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/Initializer/BaseDRInitializer.cpp \
  && echo "PASS: fault params evaluated GP-wise by default"
```

---

### [R-002] [MODERATE] [POSSIBLE] safs_fault_rs.yaml (plan lines 151-152) + line 19 — default nucleation is inherited from the SRW case and is likely too weak to PROPAGATE plain RS, contradicting "ready-to-run"

**Category:** ASSUMPTION / DEVIATION

**Description:**
The nucleation (`Tnuc_s = 24 MPa`, `R = 6000 m`) and `RS_iniSlipRate1 = 1d-12` are copied verbatim
from `safs_seisol_v2_1_0_RSSRW` (FL=103). Those values were tuned so that a STRONG-velocity-
weakening rupture (tiny critical length a_c ~ 65 m) runs away. With strong weakening REMOVED
(FL=3/4), the governing scale is the rate-and-state nucleation length
`L_nuc = mu*Dc/((b-a)*sn)` ~= 9 km (CVM) / ~12.5 km (constant), while the forced overstress
footprint is only `r_os = R*sqrt(ln(24/22.6))` ~= 1.5 km << L_nuc. The plan itself documents this
as the top risk, yet line 19 still advertises the decks as "ready-to-run" and Phase 2 presents the
24 MPa / 6000 m values as the defaults. The deliverable as written may produce only a confined
nucleation patch — which silently defeats task [2]'s rupture-speed map and task [3]'s PGV-vs-Vr
correlation (both need a propagating rupture with spatial extent), even though M0/P0 remain
trivially "verifiable" on a near-point source.

**Trigger:** Running `parameters_aging.par`/`parameters_slip.par` as-is, especially on the
constant material (largest L_nuc).

**Actual behavior:** Energy CSV `seismic_moment` likely plateaus right after `t_0` (rupture does
not escape the nucleation patch); Vr/PGV maps are near-degenerate.

**Expected behavior:** A propagating event so all three tasks have signal; or an explicit
"smoke-test-then-tune; NOT ready as-is" framing rather than "ready-to-run".

**Suggested fix:** Make the safer nucleation the default AND downgrade the readiness claim:
```diff
-      local R  = 6000.0
-      return { Tnuc_s = 24.0e6 * math.exp(-r2 / (R*R)) }
+      local R  = 8000.0          -- enlarged from 6000 (SRW value) so the forced region
+                                 -- approaches L_nuc for plain RS (FL=3/4); re-confirm by smoke test
+      return { Tnuc_s = 30.0e6 * math.exp(-r2 / (R*R)) }
```
```diff
- This is a post-processing/pedagogy deliverable plus ready-to-run SeisSol input decks.
+ This is a post-processing/pedagogy deliverable plus SeisSol input decks that MUST pass the
+ Phase-6 smoke test (rupture nucleates AND propagates) before distribution — plain RS without
+ strong weakening is prone to confined, non-propagating events; the decks ship with enlarged
+ nucleation and may still need the Risk-Assessment tuning.
```

**Test case:**
```bash
# After a constant-material run, the rupture must propagate, not stall at the patch:
python3 - <<'PY'
import csv
m=[(float(r['time']),float(r['measurement'])) for r in csv.DictReader(open('output/safs-energy.csv'))
   if r['variable']=='seismic_moment']
m0_at_t0   = next(v for t,v in m if t>=1.0)     # end of nucleation ramp t_0=1 s
m0_final   = m[-1][1]
assert m0_final > 5*m0_at_t0, f"rupture did not propagate past nucleation: {m0_at_t0:.2e} -> {m0_final:.2e}"
print("PASS: moment grows well beyond the nucleation patch")
PY
```

---

### [R-003] [MODERATE] WORKSHEET STF recipe (plan line 329) + Phase 4 — `x 3.2e10` STF/M0 is constant-material-only and is silently wrong when the same recipe is reused on CVM

**Category:** BUG

**Description:**
Phase 3 task 3 builds the moment-rate as `Integrate(ASl)(t) x 3.2e10 = M0(t)`. The constant
`3.2e10` is the homogeneous rigidity and is correct ONLY for the constant material. Phase 4 tells
participants to "recompute ... exactly as in Phase 3" on the CVM, and flags the M0 caveat for the
single final M0 (task 2) — but does NOT carry that caveat to the STF (task 3). On CVM,
`M0(t) = integral(mu(x)*slip(t) dA) != 3.2e10 * integral(slip dA)`, so a participant following the
Phase-3 STF recipe on CVM produces a wrong moment-rate curve and will (incorrectly) conclude
SeisSol disagrees.

**Trigger:** A group applies the Phase-3 STF steps to their CVM run.

**Actual behavior:** ParaView STF on CVM is off by the spatially-varying mu weighting; mismatch vs
SeisSol `seismic_moment(t)` is misread as an error.

**Expected behavior:** STF on CVM either uses SeisSol's `seismic_moment(t)` directly (no ParaView
mu needed) or the advanced resampled-mu integral; the `x 3.2e10` shortcut is labeled
constant-material-only.

**Suggested fix:**
```diff
- 3. **Moment-rate / STF** — `Integrate Variables` -> `Plot Data Over Time` of `ASl` x 3.2e10 = M0(t);
-    finite-difference for the rate. Truth: CSV `seismic_moment(t)` differenced at 0.25 s. Expect
-    matching shape; ParaView coarser (2 s).
+ 3. **Moment-rate / STF (CONSTANT material only for the `x 3.2e10` form)** — `Integrate Variables`
+    -> `Plot Data Over Time` of `ASl` x 3.2e10 = M0(t); finite-difference for the rate. Truth: CSV
+    `seismic_moment(t)` differenced at 0.25 s. Expect matching shape; ParaView coarser (2 s).
+    ON CVM: the `x 3.2e10` factor is INVALID (mu varies); use SeisSol `seismic_moment(t)` directly,
+    or the Phase-4 resampled-mu integral. Reusing `x 3.2e10` on CVM gives a wrong STF.
```
Add the same one-line warning to Phase 4 requirement 3.

**Test case:**
```bash
# On a CVM run, the constant-mu STF endpoint must NOT match SeisSol's moment (proves the caveat):
# integral(ASl)*3.2e10  vs  CSV seismic_moment(final)  -> expected to DIFFER by the mu-weighting.
echo "Manual: confirm ParaView Integrate(ASl)_final * 3.2e10 != energy-CSV seismic_moment_final on CVM"
```

---

### [R-004] [LOW] Constraints (plan lines 63-64) — "per-job < 1 GB" assumes the QuakeWorx limit is per-download, not per-user-total

**Category:** ASSUMPTION

**Description:**
The budget proves each of the 4 runs is ~0.47 GB (double). The plan states the 1 GB limit is
per-job. If QuakeWorx instead enforces a per-user or per-project storage quota, a group running
both materials (2 runs ~= 0.94 GB) or the whole class (4 runs ~= 1.9 GB) could exceed it. The
original user constraint ("does it exceed 1GB") was stated for a single download; per-job is the
reasonable reading but is unverified.

**Trigger:** Accumulated storage across a group's/class's runs on a per-user-quota system.

**Suggested fix:**
```diff
- laptop). Each of the 4 runs is its own job and must individually clear 1 GB.
+ laptop). Each of the 4 runs is its own job and individually clears 1 GB (~0.47 GB double).
+ VERIFY whether the 1 GB limit is per-job (assumed here) or a per-user/per-project quota; if the
+ latter, delete each run's download before fetching the next, or coarsen intervals (see levers).
```

**Test case:** n/a (policy verification, not code).

---

### [R-005] [LOW] parameters_aging.par `&Output` (plan line 249) — dead `refinement = 1` under `Format = 10` is confusable with the fault-refinement budget warning

**Category:** QUALITY

**Description:**
`&Output` carries `refinement = 1` (volume-output refinement), inert because `Format = 10` disables
volume output. Harmless, but the plan elsewhere warns that "fault `refinement = 1` ... blows past
1 GB" (lines 293, 451) — a reader can conflate the inert volume `refinement = 1` with the fault
`&Elementwise refinement` (correctly 0) and either panic or "fix" the wrong knob.

**Trigger:** Reader/maintainer editing the deck after seeing the budget warning.

**Suggested fix:**
```diff
- TimeInterval = 5.0
- refinement = 1
+ TimeInterval = 5.0                              ! moot (Format=10)
+ refinement = 1                                  ! VOLUME refinement; moot under Format=10. NOT the
+                                                 ! fault &Elementwise refinement (that is 0, below).
```

**Test case:** n/a (clarity).

---

### [R-006] [LOW] Mw constant consistency (plan line 327) — keep `-6.07` everywhere; an earlier value (`-9.05`/`-6.03`) was used in discussion

**Category:** QUALITY

**Description:**
The plan's `Mw = (2/3)*log10(M0) - 6.07` matches SeisSol (`EnergyOutput.cpp:654,675`) and is
correct. Flagging only to prevent the downstream WORKSHEET/ANSWER_KEY (not yet written) from
reintroducing the alternative `Mw = (2/3)(log10 M0 - 9.05)` form (= -6.03) that appeared in the
design discussion; the two differ by ~0.04 magnitude units.

**Trigger:** Authoring WORKSHEET.md / ANSWER_KEY.md from memory.

**Suggested fix:** No change to the plan. When writing the worksheet/key, use `-6.07` to match
SeisSol's own terminal/CSV-derived Mw.

**Test case:**
```bash
grep -q "2.0 / 3.0 \* std::log10" /Users/chunhuizhao/projects/SeisSol/src/ResultWriter/EnergyOutput.cpp \
 && grep -q "6.07" /Users/chunhuizhao/projects/SeisSol/src/ResultWriter/EnergyOutput.cpp \
 && echo "PASS: SeisSol uses (2/3)log10(M0) - 6.07"
```

---

## Summary
- Critical issues: 0
- Moderate issues: 3  (R-001 GPwise obsolete; R-002 nucleation too weak for plain RS vs "ready-to-run"; R-003 CVM STF mu-weighting)
- Low issues: 3  (R-004 per-job quota assumption; R-005 dead volume refinement; R-006 Mw constant consistency)
- Plan compliance: PARTIAL — decks are internally consistent and the SeisSol formulas/sizes are
  verified correct, but R-001 (obsolete key) and R-002 (untuned nucleation under a "ready-to-run"
  label) are deviations from a truly turnkey deliverable.
- Verdict: PASS WITH FIXES — apply R-001/R-002/R-003 before distributing; R-002's smoke test is the
  gating check.

## Verified-correct (audited, no issue found)
- FL mapping: aging=3, slip=4, both `<NoTP>`, shared initializer (`Factory.cpp:64,93`). One yaml,
  two FL lines — correct.
- RS required keys all present and case-insensitively matched; `rs_muw`/`rs_srW` correctly OMITTED
  for FL=3/4; spatial `rs_b` honored (`RateAndStateInitializer.cpp:138`); fallback `RS_b=0.0168`
  equals the hypocenter `a+0.004` — consistent.
- a-b sign: shallow `amb=-0.004` => `b=a+0.004` => VW (a-b<0) — correct.
- Constant material: `mu=lambda=3.2e10` => Poisson 0.25, Vs=3462 m/s, Vp=5996 m/s; harmonic-mean
  mu = 3.2e10 => `M0 = 3.2e10*P0` exact — correct (plan's Vp supersedes the archive yaml comment's
  erroneous 5.39 km/s).
- Size budget: fault geometry written ONCE (`FaultWriter.cpp:71-88`, `addSyncBuffer` at init),
  6 fault components (SRs,SRd,Vr,ASl,PSR,RT), 6 surface fields (v1..u3); ~0.47 GB/run double
  arithmetic checks out; surface forced to legacy `surfacevtkorder=-1` and fault `refinement=0`.
- ParaView P0: `Integrate Variables` on CELL data (fault output is nCells-sized => cell data) gives
  sum(value*area) == SeisSol potency formula — correct; units m^3 (P0), N*m (M0).
- Nucleation mechanism `Tnuc_s` is friction-law-independent (`BaseDRInitializer`) — valid for FL=3/4.

## Unreviewed Areas
- WORKSHEET.md and ANSWER_KEY.md do not exist yet (plan Phases 3-6 deliverables) — not reviewable;
  R-003 and R-006 are pre-emptive guards for when they are written.
- The actual deck FILES are not yet materialized (only embedded in the plan); this review audits the
  embedded content. Re-review after the files are created.
- QuakeWorx SeisSol build version/precision is unknown; R-001 (GPwise) and the single/double budget
  split both hinge on it and could not be verified against the deployed app.
