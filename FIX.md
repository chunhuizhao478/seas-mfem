# Fix Report — Round 8 (2026-05-17)

## Summary

- Findings addressed: **8 of 8** (4 MODERATE + 3 LOW + 1 POSSIBLE)
- Files modified: 6
- Tests added: 1 sub-test + 11 new grep/runtime assertions
- Test suite: **PASS** — 148/148 default; 150/150 with
  `SEAS_TEST_RUNTIME_QSIZE_CHECK=1`; 94/94 BP5 PETSc TS restart;
  TPV104 driver builds clean.

### Files modified
- `miniapps/seas/io/tpv104_checkpoint.hpp` (R-104)
- `miniapps/seas/drivers/tpv104_driver.cpp` (R-106)
- `miniapps/seas/jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch`
  (R-101 + R-103)
- `miniapps/seas/jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch`
  (R-103 BP5 sibling)
- `miniapps/seas/tests/unit/test_tpv104_checkpoint.cpp`
  (R-102 + R-105 + R-108 + R-101/R-104 grep coverage)
- `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`
  (R-107)

---

## Changes Made

### R-101 MODERATE — Stale "no-restart gap" sbatch comments

Updated 5 stale comment blocks in
`jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch`:

- Lines 127–131 (chained-restart header) — "Today, segment_002 is a
  SECOND FRESH RUN (Phase B's --restart is silently ignored)" →
  "segment_002 is RESTARTED from segment_001's end-of-run V1 checkpoint."
- Lines 287–290 (Validation #1) — "is what restart would need to
  preserve" → "matches what Phase B's --restart load expects."
- Lines 306–309 (Validation #2) — "Even though restart was a no-op
  ... TPV104 lacks safety check" → "The --output-dir separation is
  enforced by the canonical-path collision check in tpv104_driver.cpp."
- Lines 325–328 (Validation #3) — "no safety check in the TPV104
  driver" → "the --output-dir separation + driver safety check kept
  Phase B's writes out of Phase A's directory."  Also notes the BEFORE-B
  capture introduced by R-103.
- Lines 343–347 (Validation #4) — "CONFIRMS the no-restart gap" →
  "proves Phase B's V1 restart actually fired."

Also renamed the SBATCH identifiers (R-101 sub-fix):
- `#SBATCH -J tpv104_restart_scaff` → `tpv104_restart_v1`
- `#SBATCH -o tpv104_restart_scaffolding_%j.out` → `tpv104_restart_v1_%j.out`
- `#SBATCH -e tpv104_restart_scaffolding_%j.err` → `tpv104_restart_v1_%j.err`
- Updated `PHASE_LOG` at line 286 to match the new `-o` path so
  Validation #4's grep still finds the log.

Sub-test 6 extended with 5 R-101 grep assertions (sbatch must NOT
contain "silently ignored" / "no-restart gap" / "no safety check in
the TPV104 driver" / "restart was a no-op" / must contain
`FAULT_A_BEFORE_B_SIZE`).

### R-102 MODERATE — Test header wrong sbatch filename

`tests/unit/test_tpv104_checkpoint.cpp:7` —
`tpv104_restart_test_dev_2hr.sbatch` → `tpv104_restart_test_v1_dev_2hr.sbatch`.

Sub-test 6 extended with self-grep that confirms the test file
references the real sbatch and does NOT reference the placeholder.
The placeholder search-string is built from two concatenated halves so
the assertion itself does not self-match.

### R-103 MODERATE — Validation #3 broken-by-design (both sbatches)

**TPV104** (`jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch`):
- Inserted `FAULT_A_BEFORE_B_SIZE` capture right after Phase A
  finishes (line 217–225, parallel to the existing md5sum
  capture).
- Updated Validation #3 to compare `FAULT_A_SIZE_AFTER_B` against
  `FAULT_A_BEFORE_B_SIZE` (the BEFORE-B capture) instead of the
  AFTER-B reference that was guaranteed equal.

**BP5** (`jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch`):
- Same `FAULT_A_BEFORE_B_SIZE` capture inserted at line 264–270,
  between Phase A wall-time print and Phase B start.
- Updated Validation #8 (lines 462–478) to compare against the
  BEFORE-B size.  Comment now references R-103 and explains that
  without the BEFORE-B capture the check would be a no-op.

### R-104 MODERATE — Removed `expected_Q_size = -1` sentinel

`io/tpv104_checkpoint.hpp`:
- `ReadTpv104CheckpointImpl` now requires `expected_Q_size >= 0`
  unconditionally (no `if (expected_Q_size >= 0)` gate around the
  size compare).
- Added a separate `MFEM_VERIFY(expected_Q_size >= 0, ...)` that
  fires loudly if a caller passes a negative value.
- Updated both the `internal::ReadTpv104CheckpointImpl` doc block
  and the public-overload doc block to remove the "pass -1 to skip"
  language and explain the R-104 hardening.

No existing caller used `-1` (Sub-tests 1 and 9 already pass the
real Q size; the driver passes `Q.Size()`), so the sentinel removal
is a contract tightening, not a behavioural break.

Sub-test 6 R-104 grep extended:
- Header MUST NOT contain `if (expected_Q_size >= 0)` (the old
  sentinel branch).
- Header MUST contain `expected_Q_size >= 0,` (the new MFEM_VERIFY
  required-positive check).

### R-105 LOW — Mark test Sub-tests 2/3 as NOT IMPLEMENTED

`test_tpv104_checkpoint.cpp:11–14` header — Sub-tests 2 and 3 were
listed in the inventory but never implemented (pre-existing; R-008
inherited).  Marked both as `(NOT IMPLEMENTED)` with the same caveat
(MFEM_VERIFY not catchable; source-level check still in place at the
header's `MFEM_VERIFY(file_rank == rank)` and
`MFEM_VERIFY(file_num_ranks == size)` lines).

### R-106 LOW — Driver block comment `bulk.vtkhdf`

`tpv104_driver.cpp:498–503` — block comment listed `bulk.vtkhdf`
(which TPV104 does not write).  Replaced with the actual file
inventory matching the user-facing R-003 error message:
`fault.vtkhdf / volume.vtkhdf / ParaView_bulk/volume.vtkhdf /
*_station_*.dat / *_checkpoint_r*.txt`.

### R-107 LOW — Plan doc "6 sub-tests" inventory

`petsc_ts_restart_plan_2026-05-16.md:897` — expanded the test
inventory to enumerate all 8 sub-tests by ID + role + which round-7
findings they gate.  Also notes the opt-in Sub-test 8 runtime check
via `SEAS_TEST_RUNTIME_QSIZE_CHECK=1`.

### R-108 POSSIBLE → IMPLEMENTED — Opt-in Sub-test 8 runtime check

`test_tpv104_checkpoint.cpp`:
- Added `<cstdlib>`, `<sys/wait.h>`, `<unistd.h>` includes.
- Added `Subtest8_WrongQSizeRuntime()` mirroring the BP5 Sub-test 14
  opt-in subprocess pattern.  Forks; in the child, calls
  `ReadTpv104Checkpoint` with `expected_Q_size = 100` (wrong); the
  R-104 MFEM_VERIFY fires, child aborts via SIGABRT.  In the parent,
  `waitpid` reaps the child and asserts `WIFSIGNALED(status)` or
  non-zero `WEXITSTATUS`.
- Stderr in child is redirected to `/dev/null` (`freopen`) so the
  deliberate abort message doesn't pollute the test summary.
- Before `fork()`, both `std::cout`/`std::cerr` and the C `stdout`/`stderr`
  buffers are flushed so the child does not inherit + re-emit the
  parent's pre-fork output (without the flush, when stdout is a pipe
  the child's exit duplicates parent output).
- Opt-in via env var `SEAS_TEST_RUNTIME_QSIZE_CHECK=1` (default-on
  would emit a deliberate abort message even with /dev/null redirect
  on some libc builds).

Verified:
- Default: SKIP, 148/148 PASS.
- Opt-in: 2 new assertions, 150/150 PASS.

---

## Verification

### Findings checklist

- [x] **R-101** — All 5 stale comments updated; SBATCH job/log names
  renamed; PHASE_LOG updated to match.  Sub-test 6 R-101 grep PASSES.
- [x] **R-102** — Test header points at the real sbatch.  Sub-test
  6 R-102 self-grep PASSES (uses split-literal trick to avoid
  self-matching).
- [x] **R-103** — Both BP5 and TPV104 sbatches now capture
  `FAULT_A_BEFORE_B_SIZE` before Phase B starts and compare against
  it after Phase B.  Sub-test 6 grep PASSES that TPV104 sbatch
  contains `FAULT_A_BEFORE_B_SIZE`.
- [x] **R-104** — `if (expected_Q_size >= 0)` sentinel branch
  removed; `MFEM_VERIFY(expected_Q_size >= 0, ...)` added.  Sub-test
  6 R-104 grep PASSES (header has no sentinel branch + has the
  required-positive verify).  Sub-test 8 runtime PASSES with opt-in.
- [x] **R-105** — Sub-tests 2 and 3 marked NOT IMPLEMENTED in the
  test header inventory.
- [x] **R-106** — Driver block comment lists actual file paths,
  matching the R-003 user-facing error message.
- [x] **R-107** — Plan doc inventory updated to enumerate all 8
  sub-tests + the opt-in runtime check.
- [x] **R-108** — Sub-test 8 added; default SKIP, opt-in PASS via
  `SEAS_TEST_RUNTIME_QSIZE_CHECK=1`.

### Test runs

- `make seas_test_tpv104_checkpoint && ./seas_test_tpv104_checkpoint`
  → **148/148 PASS, 0 FAIL** (was 139/139 before R-101..R-104 grep
  extensions added 9 more assertions).
- `SEAS_TEST_RUNTIME_QSIZE_CHECK=1 ./seas_test_tpv104_checkpoint`
  → **150/150 PASS, 0 FAIL** (Sub-test 8 contributes 2 additional
  assertions: waitpid reaped + child aborted/non-zero).
- `make test-bp5-petsc-ts-restart` → **94/94 PASS, 0 FAIL**.  No
  regression from R-103 BP5 sbatch edit (the BP5 sbatch is not
  exercised by the unit test, but the test confirms no source-level
  regressions).
- `make seas_tpv104_driver -j8` → **clean build**.

### Pre-existing failure noted (NOT in this round's scope)

- `make test-checkpoint` (BP5 V1 checkpoint test) still fails on
  `AntiplaneDomainOperator::ComputeTractionDiagnostics` /
  `IsFirstStepDebugEnabled`.  Files unmodified by round 7 or round 8;
  unrelated to this fix work.

---

## Notes on suggested-fix deviations

- **R-102 grep self-match**: the round-8 review's suggested test
  case used the literal string `"(jobs/tpv104/tpv104_restart_test_dev_2hr.sbatch)"`
  directly in the `MFEM_VERIFY` argument.  Applied verbatim, the
  assertion self-matched (the literal appeared in its own assertion
  text, so the test file always "referenced" the bad path).  Built
  the search string from two concatenated halves so the literal
  does not appear contiguously anywhere in the source.
- **R-108 stdout buffer inheritance**: the initial implementation
  did not flush `stdout`/`stderr` before fork().  When the test was
  invoked through a pipe (e.g., `conda run … | tail`), the child
  inherited the parent's pre-fork stdout buffer and on
  child-process exit re-emitted everything, producing duplicate
  sub-test banners.  Fixed by flushing both C and C++ stdio before
  `fork()`.  The 150/150 result was correct in both cases; only the
  visible output was cosmetically duplicated.

---

## Unresolved Findings

None.  All 8 findings from REVIEW.md round 8 are addressed.

---

## Ready for Re-Review: YES

Suggested re-review focus on the next round:
1. Confirm the BP5 sbatch's R-103 fix is correctly placed (the
   `FAULT_A_BEFORE_B_SIZE` capture lives between Phase A finish and
   Phase B start — verify by simulating a clobber).
2. Whether the R-104 hardening introduces any subtle issues for
   future callers — the public API now requires every Read caller
   to know the live Q size.  The driver passes `Q.Size()`; the test
   passes the constexpr `Qsize`.  Both are correct.  A future driver
   that does not pre-size `Q` before calling Read would now
   MFEM_VERIFY-fail (which is the desired behavior — but the failure
   message should still be clear enough).
3. R-108's opt-in Sub-test 8 runs `fork()` from a C++ test binary
   that includes MFEM.  On macOS conda mfem-dev this works; on
   Frontera intel/19 the behavior should be similar but un-verified.
   Consider gating Sub-test 8 with a `__APPLE__ || __linux__` macro
   if it ever flakes on a different platform.
