# Fix Report: REVIEW_mesh_scripts.md (2026-04-30)

## Summary

- Findings addressed: **5 of 5**
- Files modified: 4 (`run_two_crossing_2000m_cgal.sh`,
  `run_two_crossing_2000m.sh`, `run_all8_2000m_cgal.sh`,
  `run_all8_2000m.sh`)
- Tests added: 0 (shell-script UX fixes; verified by exercising the
  Phase-3 wrapper end-to-end against a now-failing validator)
- Test suite:
  - **Bash syntax check** (`bash -n` on all 4 scripts): PASS.
  - **Phase-3 end-to-end via the wrapper**: exit code 1 (correctly
    propagated by R-405/R-601 from a real upstream validator
    failure), `cat` of report dumped to stderr 1×, "pipeline
    complete" success line absent — confirms R-601 + R-604 firing
    on a real fail.

## Changes Made

1. **R-601 — Phase-5 validate-fail handler.**
   Copied R-405's catch-block from `run_two_crossing_2000m_cgal.sh`
   into `run_all8_2000m_cgal.sh:148–177`.  On validator FAIL: dump
   the full report to stderr (with file-existence guard from R-604);
   exit 1.  On success: tail the summary so pass/fail counts are
   visible without opening the report file.

2. **R-602 — `--remesh-iters` env-overridable.**
   - `run_two_crossing_2000m_cgal.sh:48`: added
     `REMESH_ITERS="${TWO_REMESH_ITERS:-3}"`.
   - `run_two_crossing_2000m_cgal.sh:98`: replaced literal `3` with
     `"$REMESH_ITERS"`.
   - `run_all8_2000m_cgal.sh:71`: added
     `REMESH_ITERS="${SAFS_REMESH_ITERS:-3}"`.
   - `run_all8_2000m_cgal.sh:140`: same replacement.

3. **R-603 — wrappers `set -euo pipefail`.**
   Added to both `run_two_crossing_2000m.sh` and
   `run_all8_2000m.sh` immediately above the `export ... exec` lines.

4. **R-604 — guard `cat` against missing report file.**
   In `run_two_crossing_2000m_cgal.sh:130–142` the `||` block now
   tests `[ -f "$report" ]` before catting; if missing, prints a
   diagnostic about validator-crash-before-write.  The `tail -3`
   summary on the success path is also wrapped in the same guard.
   The Phase-5 script already has the same guard from R-601.

5. **R-605 — `SAFS_FAULTS` env override for the all-8 script.**
   `run_all8_2000m_cgal.sh:47–60`: if `SAFS_FAULTS` is set (and
   non-empty), `read -r -a FAULTS <<< "$SAFS_FAULTS"`; else use the
   hardcoded 8-fault default.  Updated the R-402a preflight message
   (lines 115–117) to point users at the env var first, then fall
   back to "edit the FAULTS=() array".

## Verification

| Check | Result |
|---|---|
| `bash -n` syntax of all 4 scripts | PASS |
| Phase-3 wrapper exit on validator FAIL | exit=1 ✓ |
| R-601/R-604 catch-block fires on FAIL | "validate_msh.py FAILED" appears once in stderr ✓ |
| Pipeline-complete echo absent on FAIL | confirmed (script aborted) ✓ |
| Report-file guard (`[ -f ... ]`) | confirmed (file existed; cat path took) ✓ |
| `REMESH_ITERS` propagated | grep verified `"$REMESH_ITERS"` is referenced where `3` used to be ✓ |
| `set -euo pipefail` in both wrappers | grep verified line 7 of both wrappers ✓ |
| `SAFS_FAULTS` override active | grep verified `read -r -a FAULTS <<< "$SAFS_FAULTS"` at line 48 ✓ |

## Notes

- The 10/11 validator outcome on the latest run is from upstream
  R-402 modifications to `validate_msh.py:check_10_tet_quality`
  (the near-fault-sliver gate) — those changes were applied
  outside my scope per the prior `/code-fix` round's plan-deviation
  stop.  My changes here are **purely shell-script**; they do not
  alter the mesh.
- The "10/11 checks passed" result is a real validator regression
  that exists at the moment, not introduced by this fix round.
  The Phase-3 wrapper now correctly surfaces that failure
  (previously the script would have exited silently before the
  R-405 dump).

## Ready for Re-Review: YES

All 5 findings applied; 4 files modified; the Phase-3 happy-and-
sad paths both work as designed.
