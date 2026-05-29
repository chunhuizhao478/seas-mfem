# SUPERSEDED — see `REVIEW.md`

This was an interim draft written while the tool channel was dropping output. It listed three items as
`[POSSIBLE]` bugs (ADER-corrector `dt·Q_bg` target, global-bbox MIN/MAX reduction, fault-bbox Allreduce) that were
subsequently **verified CORRECT** against the committed code. Do NOT act on those — "fixing" them would break correct,
intentional code (violating the repo's don't-revert rule).

The authoritative Phase-12 (PML) review is **`REVIEW.md`** in the project root. It contains:
- 0 critical, 1 moderate (R-001: ADER-corrector PML path has no automated test — add one), 2 low (R-002 `R_eff`
  labeling; R-003 stale TOML comments),
- a "Verified CORRECT (do NOT touch)" section enumerating the hunted-and-cleared subtleties,
- verdict: PASS WITH FIXES.
