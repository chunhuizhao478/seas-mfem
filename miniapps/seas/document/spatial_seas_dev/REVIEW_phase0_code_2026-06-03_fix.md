# Fix Report: Phase 0 review findings (REVIEW_phase0_code_2026-06-03.md) — 2026-06-03

## Summary
- Findings addressed: **4 of 4** (R-001, R-002, R-003 fixed this round; R-004 was already fixed in the review round).
- Files modified: `drivers/spatial_seas_driver.cpp` (only). No shared code touched; `Makefile` unchanged this round.
- Tests added: 0 new unit-test targets (all findings LOW; no CRITICAL/MODERATE requiring a test). R-001's behavior is covered by the CLI smoke (exit-2 on unreadable `--config`); a formal config unit test (`tests/unit/test_spatial_seas_config.cpp`) is owned by Phase 1 per the plan.
- Build: **clean, no errors/warnings** (worktree invocation).
- Smoke suite: all driver paths re-verified PASS (table below).

## Changes Made
1. **R-001** (dead `try/catch`; config-load failure aborted with exit 1 instead of the documented exit 2) → added `#include <fstream>` and a QD-driver-local readability pre-check (`std::ifstream probe(config_path)`) immediately before `LoadSpatialFrictionConfig`. A missing/unreadable `--config` now prints `ERROR: --config '<path>' is not readable.` on rank 0 and returns **2**. The shared parser is untouched (a readable-but-malformed TOML still aborts inside it, matching `spatial_dyn_driver.cpp`). This realizes the review's suggested fix.
2. **R-002** (banner echoed shared default `paraview_fault="hdf5"`, contradicting the QD "VTU/PVD, never vtkhdf" directive) → added a `pv_mode_echo` lambda that renders any `"hdf5"` mode as `hdf5 (Phase 6 will force vtu)` in the banner, applied to **both** `paraview_volume` and `paraview_fault` (extended to volume for the same root cause — a TOML may set volume="hdf5" too). Non-hdf5 modes (`vtu`/`off`) print verbatim with no annotation.
3. **R-003** (non-dry-run returned 0 while doing nothing → could read as a successful run) → the "Phase 0 skeleton not yet implemented" notice now goes to **stderr** (was stdout); the exit code stays **0** by design so the skeleton never fails a future build/smoke. stdout now carries only the legitimate banner.
4. **R-004** (banner `dt_initial` showed `-1 (auto)`) → already fixed in the review round (`dt_initial: auto`); no action.

## Verification
- [x] R-001: missing-file `--config` → **exit 2** + "is not readable" (was exit 1 via MFEM_ABORT).
- [x] R-002: default-hdf5 config banner shows `hdf5 (Phase 6 will force vtu)` for volume+fault; `--paraview-fault vtu --paraview-volume off` shows plain `vtu` / `off`.
- [x] R-003: non-dry-run → notice on stderr, banner on stdout, **exit 0**.
- [x] R-004: `dt_initial: auto` (already fixed).
- [x] Regression: good-config `--dry-run` still exits 0; MPI np=2 unchanged; missing `--config` still exit 2.

| Smoke | Result |
|---|---|
| missing `--config` | exit 2 ✔ |
| unreadable `--config <path>` | exit 2 + clean message ✔ (R-001) |
| good dry-run | banner + exit 0 ✔ |
| banner hdf5 annotation | volume+fault annotated ✔ (R-002) |
| banner vtu/off override | plain, no annotation ✔ (R-002) |
| non-dry-run | notice→stderr, exit 0 ✔ (R-003) |

## Unresolved Findings
- None.

## Pre-existing failure (unchanged; out of scope, not fixed)
- `seas_test_spatial_friction_config` still aborts in its `[NUM-1]` sub-case (parser↔test schema drift on `[material] kind="depth_profile_1d"`, `spatial_friction.cpp:1622`). This is in the no-touch shared parser/test, unrelated to Phase 0, and was documented in the review. Not addressed here (correctly out of scope). Recommend a separate ticket.

## Note on test-suite scope
The full `make all && make test` was not run end-to-end (heavy; contains the pre-existing failure above). My changes are confined to the new `spatial_seas_driver.cpp` — no shared object, no existing target, and no `all`/`test` aggregate is touched — so the suite's behavior is unchanged by this changeset. The driver's own behavior is covered by the CLI smokes above.

## Ready for Re-Review: YES
