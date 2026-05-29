# Fix Report: Caliper perfgraph PLAN — round 1

Source review: `caliper_perfgraph_check.md` (R-001 … R-005).
Artifact fixed: `caliper_perfgraph_plan.md` (the review found defects *in the
plan*; this round edits the plan, not C++).

## Summary
- Findings addressed: **5 of 5**
- Files modified: `document/caliper_perfgraph_dev/caliper_perfgraph_plan.md`
- Tests added: n/a in the unit-test sense (the artifact is a markdown plan); the
  review's two verification *checks* (region-presence greps for R-001/R-002) are
  now embedded in the plan's recipe + Phase-2 acceptance criteria.
- Test suite: n/a (no compilable artifact this round). Consistency verified by
  grep (all 5 fix markers present; 0 stale phrasings; code fences balanced).

## Changes Made
1. **R-001 (MODERATE) — empty-perfgraph build gap.** Phase 3 step 3 now reads
   `make clean && make seas_spatial_dyn_driver` with an inline note that
   skipping the clean is the #1 cause of an empty perfgraph (stale no-op
   objects). Applied the reviewer's diff verbatim + the explanatory sentence.
2. **R-002 (MODERATE) — friction solve not annotated.** Added a first-class
   region for the friction sub-step solve in three places: the Phase 2 Goal
   (now `… → {friction_substep, WaveOperator regions}`), a new **item 3** under
   Phase 2 "Files to Modify" (bracket `iterator.Advance(...)` at
   `spatial_dyn_driver.cpp:~449` with `MFEM_PERF_BEGIN/END
   ("seas::spatial_dyn::friction_substep")`), and the Phase 2 acceptance
   criteria. I went slightly beyond the reviewer's single suggestion by
   documenting BOTH the call-site BEGIN/END form (primary) and the
   exception-safe per-iterator `MFEM_PERF_SCOPE` alternative, with an explicit
   "pick ONE, not both (double-counting)" instruction — the reviewer noted the
   alternative but did not warn against applying both.
3. **R-003 (LOW) — masked exit code.** The verification loop now sets `rc=0`
   before the loop, `rc=1` on failure, and ends with
   `[ "$rc" -eq 0 ] || { echo "VERIFICATION FAILED"; exit 1; }`. Applied the
   reviewer's diff verbatim.
4. **R-004 (LOW) — Phase-2 coverage understated.** Added a mandatory "Sole
   behavioral gate" bullet to the Phase 2 acceptance criteria stating there is
   NO unit test over `spatial_dyn_driver.cpp` and that the line-by-line
   `git diff` inspection is the real guard. Applied the reviewer's intent.
5. **R-005 (LOW) — merge target unspecified.** Added a "Merge target
   (REQUIRED)" bullet to Phase 3 "Files to Modify": the worktree branch must
   merge into `system/spatial_dyn_driver` (where the `TPV31_PERFGRAPH` hook +
   p1/p2 jobs live), else the annotations and their runner are separated.

### Cascading consistency fixes (not separate findings, required by R-002)
Adding the friction `MFEM_PERF_BEGIN/END` pair made two earlier "only
`MFEM_PERF_SCOPE`" statements inaccurate; updated both so the reviewer's
"reject any non-macro hunk" rule still matches what the plan actually asks for:
- Constraints bullet "Zero numerical/perf impact …" now lists `MFEM_PERF_SCOPE`
  *and* the one `MFEM_PERF_BEGIN/END` pair as the permitted edits.
- Risk Assessment "Accidental logic edit …" mitigation updated likewise.
(Phase 1 acceptance criteria intentionally left as "MFEM_PERF_SCOPE lines only"
— Phase 1 touches `wave_operator.inl`/`bimaterial_wave_operator.inl`, which are
SCOPE-only; the BEGIN/END pair is a Phase 2 driver edit.)

## Unresolved Findings
- None. All 5 findings resolved in this pass.

## Deviations from suggested fixes
- R-002: enhanced (not contradicted) the reviewer's fix — added the
  "pick ONE approach" guard and the dual call-site/per-iterator options. The
  reviewer's primary call-site BEGIN/END suggestion is the recommended path.
- No reviewer suggestion was found to be wrong; all were applied as-is or
  enhanced for plan consistency.

## New Tests
- None addable at the plan stage. The R-001 and R-002 verification *checks*
  (`grep 'seas::WaveOperator::Mult'` and `grep 'friction_substep'` against a
  `runtime-report` capture) from the review now live in the plan as the
  Caliper-build acceptance criteria — they become executable once the
  annotations + a Caliper build exist (Frontera step).

## Ready for Re-Review: YES
The plan now self-consistently specifies: per-step + macro-step + friction +
wave-operator regions; a clean Caliper build recipe; a fail-loud local
unit-test gate with the explicit Phase-2 diff guard; and a pinned merge target.
Next step is `/code-implement` against the corrected
`caliper_perfgraph_plan.md`.
