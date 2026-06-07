# Fix Report — Cross-rank Phase 3 (against REVIEW_cross_rank_phase3_code_2026-06-07.md)

**Date:** 2026-06-07. Review verdict MERGEABLE, **0 critical**, 3 LOW. Proceeding after fix.

## Changes made
- **P3-002 (LOW, the valuable one):** the gate previously asserted only the reclassify EFFECT
  (face out of central set) — the `shared = N` count was log-only. Added a member
  `n_reclass_shared_` + public accessor `GetNReclassShared()` (set at the end of
  `BuildPerFaceCentralFluxMatrices_` from the per-rank `n_reclass_shared`), and asserted it in the
  gate: STRONG ⇒ `GetNReclassShared() == 1`, WEAK ⇒ `== 0`. This also confirms BOTH ranks reclassify
  the seam face (each reports 1), corroborating the symmetric-guard finding.
- **P3-001 (LOW, accepted):** the "a 2-rank seam face counts twice" label is honest for any np (a
  shared face is always a 2-element cut); no change.
- **P3-003 (LOW, accepted):** histogram-excludes-shared log wording already clarified in Phase 2
  (P2-004); no further change.

## Verification (re-run after fix)
- Gate `seas_test_bimaterial_seam_material_np2` (np=2): **135/135** (incl. `GetNReclassShared()==1`
  strong / `==0` weak).
- Regressions: parity **9/9**, dispatch **10/10**, mixed-flux-shared (np=2) **16/16**.

## Ready for next phase: YES (0 critical; the count gate strengthens the reclassify proof; tests green).
