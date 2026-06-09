# Fix Report — Cross-rank Phase 4 (against REVIEW_cross_rank_phase4_code_2026-06-07.md)

**Date:** 2026-06-07. Review verdict: **0 critical**, 2 moderate (both in TESTS), 2 low. The
reviewer independently confirmed the production code is correct on every high-risk axis (R-204
index spaces, R-401 peer Jacobian, interior byte-exactness, collective safety, `.find()+continue`).
All findings addressed. No critical ⇒ proceeding after fix.

## Changes made (all in `tests/parallel/test_bimaterial_seam_fault_np2.cpp`, except P4-004)
- **P4-001 (MODERATE — the flagged one):** the contrast gate now does a PLUS/MINUS-SENSITIVE check —
  the parallel shared fault's `Zp_plus`/`Zp_minus` are compared to the SERIAL interior oracle's
  UNSORTED pair (same geometry ⇒ same FaultBasis +/- convention). A plus/minus transposition in the
  shared loop's `e1_plus` swap would make parallel `Zp_plus` match the serial `Zp_MINUS` ⇒ the gate
  fails. (The old sorted `{min,max}` compare could not catch a transposition.) Confirmed it passes at
  worst rel 0.0 — i.e. +y/strong → `Zp_plus` on both serial and parallel, no transposition. The
  shared loop is the single +/- code path for all modes, so the Coefficient contrast gate guards the
  GridFunction gate's swap too. (Did NOT add the suggested `TamperSharedFaultElem1OnPlus` negative
  control — it would require a test-only mutator on the production class; the unsorted oracle compare
  already fails on a real transposition.)
- **P4-002 (MODERATE):** the constant + depthprofile gates now wrap their per-rank assertions in
  `if (nshr > 0)` with an `n_checked > 0` guard, so a rank with no shared fault face cannot pass
  vacuously (the global `glob_shr >= 1` still ensures some rank checks). The contrast/gridfunction/
  index_mapping gates already had this guard.
- **P4-003 (LOW):** removed the dead `const real_t v3y = skew ? 1.0 : 1.0;`; v3's y is `1.0`
  directly, with a comment that `skew` shifts the apex in x/z (keeping the fault in the y=0 plane).
- **P4-004 (LOW):** corrected the plan-doc Phase-4 text from `shared_fault_dof_offset_.at(sf)` to
  `find(sf)->second` with `find()+continue` (matches the implemented code, which never throws).

## Verification (re-run after fixes)
- Gate `seas_test_bimaterial_seam_fault_np2` (np=2): **24/24** (incl. the new plus/minus-sensitive
  contrast comparison at 0.0; depthprofile symmetric on the skewed peer at 0.0 = R-401).
- Regressions: perside-material **11/11** (interior fault byte-exact), parity **9/9**, central **6/6**,
  seam_material np=2 **135/135**, mixed-flux-shared np=2 **16/16**.

## Ready for next phase: YES (0 critical; the +/- transposition coverage gap is closed; tests green).
