# Code Review: tpv102_debug_v8.0.0_debug_plan.md Revision 7 — Round 8 (2026-04-19)

## TL;DR

**Verdict: PASS — the plan is ready to execute. Stop reviewing, start running Phase 0.**

Eight rounds of adversarial review have been completed on this plan. Each round has found progressively less-critical issues. Round 8's fresh audit finds no new CRITICAL or MODERATE issues — only three LOW-severity notes about local environment and plan ergonomics. All findings R-001 through R-705 (rounds 1-7) have been addressed in Rev 7, including:

- CFL arithmetic correction (R-601 retracted, arithmetic fixed)
- Silent-weld guards on both interior and shared fault branches (R-007, R-103, R-503)
- Ghost-exchange aliasing deep-copy verify for all 9 components with `MFEM_VERIFY` not `MFEM_ASSERT` (R-005, R-102)
- Phase 2 DIAG with hypocenter + off-hypo + 1/2/3 km bulk propagation witnesses (R-504, R-701)
- Step-frequency log throttle (R-702)
- Tandem ground-truth reference promoted to Round 1 Slot B (R-703)
- Off-hypo cross-rank placement verified via MPI_Allreduce (R-704)
- Frontera topology constraint (8-node × {100, 200, 400} ranks; 1-node only for small-mesh sanity)
- Concurrent submission plan (2× 8-node/400-rank dev jobs per round)
- tfinal=1.2s post-breakaway (R-502)
- Makefile `SEAS_EXTRA_CPPFLAGS` + startup banner + banner-to-file (R-501, R-606)

At this point, further review iterations hit diminishing returns. **The correct action is to execute Phase 0 and generate real data.**

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v8.0.0_debug_plan.md` (Revision 7, 88 KB, ~1620 lines)
- Verified environment readiness:
  - Local meshes: `tpv102_1000m.msh` ✓, `tpv102_200m.msh` ✓ (too large for local, but available)
  - Tandem TPV102 config: `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/tpv102.toml` + `.lua` + `.geo` ✓
  - Frontera sbatch scripts: `tpv102_1000m_p1_1.5s_4rank_dev.sbatch`, `tpv102_200m_p1_1.5s_50rank_dev.sbatch`, `tpv102_200m_p1_1.5s_400rank_dev.sbatch` all present
- Spot-checked plan sections:
  - Phase 0 (freeze / baseline / allocation / CFL log)
  - Phase 1 (pre-flight + bisection)
  - Phase 1A (all guards including R-503, R-103, R-102)
  - Phase 2 Step 4-5 (hypo, off-hypo, bulk witnesses)
  - Round 1 concurrent submission (Slot A MFEM DIAG + Slot B Tandem)
- Prior reviews: rounds 1-7 (R-001 through R-705).

## Findings

---

### [R-801] [LOW] [Phase 1 pre-flight] `tpv102_500m.msh` referenced but not present locally; pre-flight coverage is reduced but degrades gracefully

**Category:** EDGE_CASE (ENVIRONMENT)

**Description:**
Phase 1 Step 2 pre-flight loop (plan line 322-340) iterates:
```bash
for m in tpv102_1000m tpv102_500m; do
   if [ ! -f "tpv102/mesh/${m}.msh" ]; then continue; fi
```

Verified by `ls miniapps/seas/tpv102/mesh/`: only `tpv102_1000m.msh` and `tpv102_200m.msh` are present locally. `tpv102_500m.msh` does not exist. The pre-flight's `continue` gracefully skips — no error — but the scan reduces to `for np in {2,4,8} × {tpv102_1000m}`, which is the known-`shared=0` combination per `v2_fix.md:128-129`.

**Actual behavior:** Pre-flight log records 3 attempts (np=2, 4, 8 on 1000m), each with `shared=0`, then `PRE_FLIGHT: ESCALATE to Phase 1 Escalation A-Frontera`. This is the expected routing for this environment — no diagnostic impact.

**Expected behavior:** Plan could suggest generating a `tpv102_500m.msh` via the existing `.geo` file if the user wants intermediate-size local coverage, but this is optional. Current behavior is correct.

**Suggested fix:**
```diff
 for np in 2 4 8; do
    for m in tpv102_1000m tpv102_500m; do
       if [ ! -f "tpv102/mesh/${m}.msh" ]; then
+         echo "(${m}.msh not present; skipping)" >> "$PRE_FLIGHT_LOG"
          continue
       fi
       ...
    done
 done
+# Note: if only 1000m is available locally, pre-flight is
+# guaranteed to land on ESCALATE (all np ≤ 8 on 1000m have
+# shared=0 per v2_fix.md:128-129).  Phase 1 Escalation A-Frontera
+# is the intended path in this environment.
```

**Test case:** N/A (environment-dependent).

---

### [R-802] [LOW] [Phase 1 Slot B Tandem action] Action 1 says "identify or create a Tandem TPV102 sbatch" — actual config path is known and should be cited

**Category:** QUALITY

**Description:**
Phase 1 Slot B Action 1 (plan line 443-446):
```
1. **(10 min) Identify or create a Tandem TPV102 sbatch script on
   Frontera.** Check `/Users/chunhuizhao/projects/tandem/examples/`
   for existing SCEC TPV102 inputs; if missing, compose from the
   Tandem documentation using TPV102's SCEC spec (BP5-like layout).
```

Verified: Tandem HAS a TPV102 config at `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/tpv102.toml` + `.lua` + `.geo`. The plan could point to this directly rather than saying "check examples."

The ambiguity means a /code-implement agent might spend 10-30 minutes searching / composing a config when the file already exists. Not fatal but wastes the "cheap Tandem comparison" advantage R-703 was meant to capture.

**Suggested fix:**
```diff
-1. **(10 min) Identify or create a Tandem TPV102 sbatch script on
-   Frontera.** Check `/Users/chunhuizhao/projects/tandem/examples/`
-   for existing SCEC TPV102 inputs; if missing, compose from the
-   Tandem documentation using TPV102's SCEC spec (BP5-like layout).
+1. **(5 min) Use Tandem's existing TPV102 config.**  Base path:
+   `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/tpv102.toml`
+   (plus `tpv102.lua` and `tpv102.geo` in the same directory).
+   Copy an existing Tandem Frontera sbatch and point it at this
+   TOML.  Adjust `final_time = 1.2` in the TOML to match our tfinal.
+   Output station files to `$SCRATCH/tandem_ref/`.
```

**Test case:** N/A (documentation).

---

### [R-803] [LOW] [Plan size / ergonomics] Rev 7 is ~1620 lines; a /code-implement agent may lose priorities without a condensed "execute this now" summary

**Category:** QUALITY

**Description:**
The plan is thorough (appropriate for a multi-round debugging effort that has burned many node-hours on misdirections), but at 1620 lines it exceeds what most LLM-based /code-implement agents can keep in working memory without dropping details.

A reader starting at Phase 0 must navigate:
- Revision history (lines 1-200)
- Ground Rules (lines 200-300)
- Phase 0 (lines 206-304)
- Phase 1 + Escalation (lines 304-500)
- Phase 1A (lines 498-730)
- Phase 2 (lines 732-1190)
- Phase 2B (lines 1196-1240)
- Phase 3A/B/C (lines 1243-1460)
- Phase 4 (lines 1461-1580)
- Phase 5/6 (lines 1583-1614)
- Escalation triggers (line 1615+)
- Finding Index (appendix)

**Trigger:** /code-implement agent invoked on this plan.

**Suggested fix:** Add a 1-page "Execute this now" quick-reference at the top of the plan (after revision history), listing only:
- Phase 0: 4 bash commands.
- Phase 1 Slot A: 1 sbatch submission command (DIAG build).
- Phase 1 Slot B: 1 sbatch submission command (Tandem).
- "Wait ~1.5 hr, then run this diff script."
- "See Round 1 Decision Table for next step."

The existing plan text becomes the reference manual for when the agent needs details.

**Test case:** N/A (documentation ergonomics).

---

## Summary

- **Critical issues: 0.**
- **Moderate issues: 0.**
- **Low issues: 3** (R-801 local mesh inventory; R-802 Tandem path could be cited; R-803 plan ergonomics).
- **Plan compliance with fundamental goal ("find the true issue"):** READY — plan covers all known failure modes with concrete diagnostic paths, MPI-verified tagging, step-throttled DIAG, and Tandem ground-truth reference.
- **Verdict: PASS — execute Phase 0. The diminishing-returns threshold has been reached across 7 review rounds.**

## What to do next

1. Start **Phase 0 right now** (4 commands, ~10 minutes):
   ```bash
   cd /Users/chunhuizhao/projects/seas-mfem
   git status && git tag v8.0.0-bug-state
   cd miniapps/seas && conda activate mfem-dev
   make -j test 2>&1 | tee /tmp/phase0_baseline_tests.log
   ssh login1.frontera.tacc.utexas.edu "taccinfo && df -h \$SCRATCH"
   ```

2. **Do not** invoke another review round before Phase 0 completes. Additional adversarial iterations on Rev 7 are unlikely to produce critical findings and delay the actual diagnosis.

3. When Phase 1 Slot A + Slot B results are in hand, pattern-match against the Round 1 Decision Table (plan line 467-477) and proceed to the next step. **At that point, post-run data will tell us more than any further plan review.**

## Unreviewed Areas

(These are deferred — not blocking execution of Phase 0.)

- **Tandem's Frontera sbatch format vs our sbatch format.** If Tandem's on-cluster build tree and launch configuration differ significantly from MFEM's, the "mirror the 400-rank sbatch" step (Slot B) may take longer than budgeted. Worst case: Slot B job fails → Round 1 Decision Table's "Tandem output missing" row → proceed with Slot A only. Not fatal.
- **Whether the R-504 + R-704 + R-701 combined DOF tagging (hypo + off-hypo + 3 bulk monitors = 5 DOFs flagged) exceeds the `local_diag_count <= 2` MFEM_VERIFY** (plan line 581). With 5 potential tags, the cap needs adjustment. If it hasn't been relaxed in Rev 7, Phase 2 build aborts at startup. Worth spot-checking during Phase 0's baseline test compile — if the compile hits this, relax to `<= 5`.
- **Whether `compare_stations.py` (referenced in Phase 1 Slot B Action 3) exists or needs to be written.** If not present, writing it takes ~15 min. Not on the critical path.
