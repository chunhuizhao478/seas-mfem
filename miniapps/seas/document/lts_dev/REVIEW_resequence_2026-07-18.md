# Plan-logic review: LTS phase resequencing (rev 5) — 2026-07-18

> Written to a dedicated file (not the root `REVIEW.md`, which the plan references
> for the P-001…P-024 findings) so that record is preserved.

## Scope
- Plan: `document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md` (rev 5 resequencing).
- Question: is the re-ordered phase sequence **logically sound** — no forward
  dependency, no acceptance criterion that can't hold in its new phase, no
  self-contradiction between the resequencing narrative and the phase bodies?
- Method: 3 auditors (dependency-graph / acceptance-criteria / claims-and-estimates)
  × adversarial verification of every candidate (17 raw, incl. 3 self-seeded),
  refute-by-default.

## Result: 17 raw → 2 confirmed (both LOW), both fixed

### [AC-1] LOW — B.8 phase tag stale after the reorder moved to Phase 3
Appendix B tagged all of `test_lts_checkpoint_v2` (B.8) as Phase 2, but its last
case — "dof_data on-disk canonical order under reorder" — depends on the reorder
+ permutation table, which rev 5 moved to Phase 3 step 0. Phase 2 is bulk-only
(no fault dof_data), so that case can't run there.
**Fixed:** split the row — B.8 (refusal + round-trip + GTS-V1 cases) stays Phase 2;
new **B.8b** (the reorder case) tagged Phase 3. The binding Phase-2 acceptance line
already said "refusal paths (B.8)", so no gate changed.

### [AC-4] LOW — canary gate still tagged Phase 1 in the test matrix
The "extended existing" row tagged the "TPV104-spatial stations byte-identical
under lts=rate2 + still-GTS" reorder canary as `(1)`, but rev 5 moved it to Phase 3
step 0 (and Phase-1 acceptance lists it as "→ Phase 3").
**Fixed:** re-tagged `(3, step 0)`, wording updated to "still-GTS reorder canary".

## Refuted (representative — why the resequencing holds)
- **Reorder cluster-id source swap (np=1 local → serial) forces re-touching the 8
  P-006 consumers** — REFUTED. The swap is upstream at the *clustering* stage;
  the reorder and its consumers depend only on the resulting ids/permutation, so
  Phase 4 re-points the id source without touching the consumers.
- **Phase 3/4 estimates not bumped ⇒ total undercounts** — REFUTED. Phase 1's
  estimate line already accounts for the full 1.5 wk and states where the
  undelivered parts went; bumping Phase 3/4 would *double-count*. Leaving them is
  the more correct choice.
- **"Phase 2 unblocked" is false** — REFUTED. Phase 2 is bulk-only np=1 and needs
  none of the trio; the checkpoint's "serial-mesh element order" uses a documented
  np=1 local-order stand-in.
- **Checkpoint hash unstable across phases** — REFUTED as a soundness hole: the V2
  hash refuses mismatches by design, and Phase-2-era LTS checkpoints are
  throwaway across the code-version change.

## Verdict
**PASS with two LOW doc fixes applied.** The resequenced dependency graph is
acyclic, no phase depends on a later phase's runtime output, and the "unblocked /
moved" classifications hold. The only defects were two stale phase tags in the
Appendix-B test matrix (rev 5 had declared Appendix B "unchanged"); both are now
corrected and the rev-5 note updated to say so.
