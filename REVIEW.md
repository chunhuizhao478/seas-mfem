# Code Review: Round 8 — Post-Round-7 TPV104 Fix Audit (2026-05-17)

## Review Scope

- Plan: `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`
  (V2 PETSc TS restart + Phase 4 TPV104 V1 restart)
- Files reviewed (round 7 fix outputs):
  - `miniapps/seas/io/tpv104_checkpoint.hpp` (R-007 dedup + R-002 expected_Q_size)
  - `miniapps/seas/drivers/tpv104_driver.cpp` (R-001/R-002/R-003/R-005/R-006)
  - `miniapps/seas/jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch` (R-004/R-009)
  - `miniapps/seas/tests/unit/test_tpv104_checkpoint.cpp` (R-008 extensions)
  - `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md` (R-006 plan note)
  - `miniapps/seas/Makefile` (test target wiring)
  - `miniapps/seas/jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch` (BP5 parity check for R-103)
- Domain context consulted:
  - `CLAUDE.md` (root)
  - `miniapps/seas/CLAUDE.md` (Phase 6.4 R-310 paths, FaultBasis convention)
  - `FIX.md` (round-7 implementer report)
- Round-7 round-trip verification: `make seas_test_tpv104_checkpoint && ./seas_test_tpv104_checkpoint`
  → 139/139 PASS, 0 FAIL.  Cross-overload Sub-test 9 PASS confirms R-007 dedup forwarding works.

This round is a FRESH adversarial pass — I am hunting for new bugs the
round-7 fixes introduced or for round-7 misses.  I am NOT just checking
off R-001 through R-009.

## Findings

### [R-101] MODERATE [tpv104_restart_test_v1_dev_2hr.sbatch:127–131,287–290,309,326–328,343–347] — Stale "no-restart gap" / "restart silently ignored" comments throughout the sbatch

**Category:** QUALITY (documentation drift that actively misleads)

**Description:**
The round-7 R-009 fix updated SOME stale comments to reflect that V1
restart now works, but missed at least FIVE stale comment blocks that
still claim restart is a no-op or that the driver has no safety check:

- **Line 127–131** (chained-restart header):
  > "Today, segment_002 is a SECOND FRESH RUN (Phase B's --restart is
  > silently ignored — see header)."
  WRONG.  V1 restart is fully wired; Phase B actually restarts.

- **Line 287–290** (Validation #1 comment):
  > "the basic output layout is what restart would need to preserve."
  Implies restart doesn't yet preserve it.  WRONG.

- **Line 306–309** (Validation #2 comment):
  > "Even though restart was a no-op, the --output-dir separation is
  > the operator's primary protection against clobber.  Same gate as
  > BP5's safety check (which TPV104 lacks)."
  Both clauses are WRONG.  Restart is not a no-op; TPV104 has the
  safety check (added at `tpv104_driver.cpp:495-579`).

- **Line 325–328** (Validation #3 comment):
  > "This is the only thing protecting V1 results today (no safety
  > check in the TPV104 driver)."
  WRONG.  TPV104 has the safety check.

- **Line 343–347** (Validation #4 comment):
  > "Acceptance check #4 — CONFIRMS the no-restart gap."
  WRONG.  Validation #4 confirms restart FIRED, not that it didn't.

A future debugger reading the sbatch will be misled into thinking
restart is broken.  Combined with the (also-stale) `SBATCH -J
tpv104_restart_scaff` job name and `tpv104_restart_scaffolding_%j.out`
log paths (lines 2–4), the artefact still presents itself as
"scaffolding only" when it is now a real end-to-end test.

**Trigger:**
Read the sbatch.

**Actual behavior:**
Comments contradict the validation logic.  E.g., line 343 says #4
"CONFIRMS the no-restart gap" but the actual #4 block (lines 348–365)
checks for `"TPV104 restart loaded:"` log line — the opposite signal.

**Expected behavior:**
All comments reflect the post-R-001..R-009 reality: restart works,
safety check is in place, validations gate end-to-end correctness.

**Suggested fix:**
```diff
@@ jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch:127-131
-# Chained-restart naming convention (matches the BP5 sbatch +
-# bp5_verification_full.cpp safety-check hint): single base dir +
-# segment_NNN subdirs.  Today, segment_002 is a SECOND FRESH RUN
-# (Phase B's --restart is silently ignored — see header).
+# Chained-restart naming convention (matches the BP5 sbatch +
+# bp5_verification_full.cpp safety-check hint): single base dir +
+# segment_NNN subdirs.  segment_002 is RESTARTED from segment_001's
+# end-of-run V1 checkpoint via `--restart ${RESULT_DIR_A}/${OUTPUT_PREFIX_A}`.

@@ :287-290
-# Acceptance check #1 — Phase A produced fault.vtkhdf.  Sanity check
-# that the driver is functioning at all and the basic output layout
-# is what restart would need to preserve.
+# Acceptance check #1 — Phase A produced fault.vtkhdf.  Sanity check
+# that the driver is functioning at all and Phase A's output layout
+# matches what Phase B's --restart load expects.

@@ :306-309
-# Acceptance check #2 — Phase B produced its own fault.vtkhdf in a
-# DIFFERENT directory.  Even though restart was a no-op, the
-# --output-dir separation is the operator's primary protection against
-# clobber.  Same gate as BP5's safety check (which TPV104 lacks).
+# Acceptance check #2 — Phase B produced its own fault.vtkhdf in a
+# DIFFERENT directory.  The --output-dir separation is enforced by
+# the canonical-path collision check in tpv104_driver.cpp:495-579;
+# this validation confirms Phase B's outputs landed in segment_002.

@@ :325-328
-# Acceptance check #3 — Phase A's fault.vtkhdf untouched after Phase B
-# finishes.  Confirms the --output-dir separation protects Phase A's
-# output from any partial writes by Phase B.  This is the only thing
-# protecting V1 results today (no safety check in the TPV104 driver).
+# Acceptance check #3 — Phase A's fault.vtkhdf size unchanged across
+# Phase B's run.  Confirms the --output-dir separation + driver
+# safety check kept Phase B's writes out of Phase A's directory.
+# (See R-103 — current implementation is broken-by-design; both stat
+# calls happen AFTER Phase B finishes.)

@@ :343-347
-# Acceptance check #4 — CONFIRMS the no-restart gap.  Phase B's stdout
-# SHOULD contain the "TPV104 restart loaded:" log line printed by the
+# Acceptance check #4 — proves Phase B's V1 restart actually fired.
+# Phase B's stdout MUST contain the "TPV104 restart loaded:" log line
+# printed by the
```

Also rename the SBATCH job/log identifiers:
```diff
@@ :1-4
 #!/bin/bash
-#SBATCH -J tpv104_restart_scaff
-#SBATCH -o tpv104_restart_scaffolding_%j.out
-#SBATCH -e tpv104_restart_scaffolding_%j.err
+#SBATCH -J tpv104_restart_v1
+#SBATCH -o tpv104_restart_v1_%j.out
+#SBATCH -e tpv104_restart_v1_%j.err
```
**WARNING**: the rename will require updating `PHASE_LOG` at line 286
to match (`tpv104_restart_v1_${SLURM_JOB_ID}.out` instead of
`tpv104_restart_scaffolding_${SLURM_JOB_ID}.out`), or else
Validation #4 will FAIL — it greps `${PHASE_LOG}` for the restart-loaded
line.

**Test case:**
Add a grep assertion to `test_tpv104_checkpoint.cpp` Sub-test 6 that
the sbatch does NOT contain the misleading strings:
```cpp
std::ifstream sb("jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch");
std::stringstream sbuf; sbuf << sb.rdbuf();
const std::string sbsrc = sbuf.str();
MFEM_VERIFY(sbsrc.find("silently ignored") == std::string::npos,
            "R-101: sbatch must not claim restart is silently ignored");
MFEM_VERIFY(sbsrc.find("no-restart gap") == std::string::npos,
            "R-101: sbatch must not claim there is a 'no-restart gap'");
MFEM_VERIFY(sbsrc.find("no safety check in the TPV104 driver")
            == std::string::npos,
            "R-101: sbatch must not claim TPV104 driver has no safety check");
MFEM_VERIFY(sbsrc.find("restart was a no-op") == std::string::npos,
            "R-101: sbatch must not claim restart was a no-op");
```

---

### [R-102] MODERATE [test_tpv104_checkpoint.cpp:7] — Header comment references the non-existent `tpv104_restart_test_dev_2hr.sbatch`

**Category:** BUG (documentation — same as the round-7 R-005 in driver, not propagated to the test)

**Description:**
`test_tpv104_checkpoint.cpp:7` reads:

```cpp
// (jobs/tpv104/tpv104_restart_test_dev_2hr.sbatch).
```

This is the identical wrong filename that round 7's R-005 fixed in the
driver (`tpv104_driver.cpp:534`).  The test file's header comment
references the SAME non-existent file.  The actual file is
`tpv104_restart_test_v1_dev_2hr.sbatch` (with the `_v1_` infix).

Sub-test 6's own R-005 grep check on the driver does not cover the
test file's own header, so the regression test for R-005 silently
allowed the same bug to survive in a sibling file.

**Trigger:**
Read the test header.  Try to `ls` the sbatch it points at.

**Actual behavior:**
`ls jobs/tpv104/tpv104_restart_test_dev_2hr.sbatch` → "No such file".

**Expected behavior:**
Header references the real sbatch filename.

**Suggested fix:**
```diff
@@ tests/unit/test_tpv104_checkpoint.cpp:7
-// (jobs/tpv104/tpv104_restart_test_dev_2hr.sbatch).
+// (jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch).
```

**Test case:**
Extend Sub-test 6 to grep its OWN source file too:
```cpp
// R-102: test header must reference the actual sbatch filename.
std::ifstream self("tests/unit/test_tpv104_checkpoint.cpp");
std::stringstream selfbuf; selfbuf << self.rdbuf();
const std::string selfsrc = selfbuf.str();
MFEM_VERIFY(
   selfsrc.find("tpv104_restart_test_v1_dev_2hr.sbatch")
   != std::string::npos,
   "R-102: test header must reference the actual sbatch name");
MFEM_VERIFY(
   selfsrc.find("(jobs/tpv104/tpv104_restart_test_dev_2hr.sbatch)")
   == std::string::npos,
   "R-102: test header must not point at a non-existent sbatch");
```

---

### [R-103] MODERATE [tpv104_restart_test_v1_dev_2hr.sbatch:288–341, bp5_restart_test_v2_dev_2hr.sbatch:440–474] — Validation #3 fault.vtkhdf-unchanged check is a no-op (broken-by-design)

**Category:** BUG (logic error inherited from BP5 sbatch; round-7 R-004 fixed the analogous station-file gap correctly, but Validation #3 still has the same defect)

**Description:**
The intent of TPV104 Validation #3 (and BP5 Validation #8) is to
verify Phase B did NOT clobber Phase A's `fault.vtkhdf`.  The
implementation:

```bash
# Line 294-297 (Validation #1, runs AFTER Phase B finishes):
FAULT_A_SIZE=$(stat ... "${FAULT_A}" ...)
...
# Line 332-335 (Validation #3, also AFTER Phase B finishes):
FAULT_A_SIZE_AFTER_B=$(stat ... "${FAULT_A}" ...)
if [ "${FAULT_A_SIZE_AFTER_B}" = "${FAULT_A_SIZE}" ]; then
   echo "  PASS — Phase A fault.vtkhdf size unchanged"
```

Both `stat` calls happen AFTER Phase B has already cleanly exited
(line 257–262 verifies `PHASE_B_RC == 0`).  After Phase B exits, NO
process is writing to Phase A's file — so the two `stat` calls are
GUARANTEED to return the same byte count.  The `if` always evaluates
true; Validation #3 always PASSes.

The check cannot catch a real clobber.

The new Validation #7 (R-004) does this correctly: it captures
md5sums BEFORE Phase B starts (line 211–213, between Phase A and
Phase B) and re-captures after (line 439–441).  The same pattern
should be applied to fault.vtkhdf.

The same defect exists in BP5's
`jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch:440–474` (Validation #8),
because the TPV104 sbatch was modelled on it.

**Trigger:**
Phase B intentionally clobbers Phase A's fault.vtkhdf (e.g., remove
the safety check, point `--output-dir` at Phase A's dir).  Validation
#3 still PASSes.

**Actual behavior:**
`FAULT_A_SIZE` and `FAULT_A_SIZE_AFTER_B` are both captured after the
clobber finished, so they match the clobbered size.  PASS is reported.

**Expected behavior:**
Capture `FAULT_A_SIZE` (or md5sum) BEFORE Phase B runs.  Re-capture
AFTER Phase B finishes.  Compare — that's the only way to detect a
clobber.

**Suggested fix:**
Move FAULT_A_SIZE capture to BEFORE the Phase B ibrun (parallel to
the new md5sum capture at lines 210–215):

```diff
@@ jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch:206-216
 # Capture Phase A station-file fingerprints BEFORE Phase B runs.
 ...
 echo "Captured $(wc -l < "${PHASE_A_STATION_HASHES_BEFORE_B}") Phase A"
 echo "  station-file md5sums before Phase B starts (R-004 unchanged check)."

+# Capture Phase A fault.vtkhdf size BEFORE Phase B runs.  Used by
+# Validation #3 (R-103) to confirm Phase B did not mutate Phase A's
+# fault output.  Without this, Validation #3 is a no-op (both stats
+# happen after Phase B exits and are guaranteed equal).
+FAULT_A_BEFORE_B_SIZE=$(stat -c %s "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
+                     || stat -f %z "${RESULT_DIR_A}/fault.vtkhdf" 2>/dev/null \
+                     || echo 0)
+echo "Captured Phase A fault.vtkhdf size: ${FAULT_A_BEFORE_B_SIZE} bytes"
+echo "  (before Phase B starts; R-103 unchanged check)."
```

Then update Validation #3:

```diff
@@ jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch:332-335
 FAULT_A_SIZE_AFTER_B=$(stat -c %s "${FAULT_A}" 2>/dev/null \
                     || stat -f %z "${FAULT_A}" 2>/dev/null \
                     || echo 0)
-if [ "${FAULT_A_SIZE_AFTER_B}" = "${FAULT_A_SIZE}" ]; then
-   echo "  PASS — Phase A fault.vtkhdf size unchanged (${FAULT_A_SIZE})"
+if [ "${FAULT_A_SIZE_AFTER_B}" = "${FAULT_A_BEFORE_B_SIZE}" ]; then
+   echo "  PASS — Phase A fault.vtkhdf size unchanged"
+   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, still ${FAULT_A_SIZE_AFTER_B} after B)"
 else
-   echo "  FAIL — Phase A fault.vtkhdf was modified during Phase B"
-   echo "         (was ${FAULT_A_SIZE}, now ${FAULT_A_SIZE_AFTER_B})"
+   echo "  FAIL — Phase A fault.vtkhdf was modified during Phase B"
+   echo "         (was ${FAULT_A_BEFORE_B_SIZE} before B, now ${FAULT_A_SIZE_AFTER_B})"
    exit 5
 fi
```

Apply the equivalent fix to the BP5 sbatch's Validation #8 at
`jobs/bp5/bp5_restart_test_v2_dev_2hr.sbatch:440–474`.

**Test case:**
Local synthetic fixture:
```bash
mkdir -p /tmp/rfix && cd /tmp/rfix
echo "phase A original" > fault.vtkhdf
FAULT_BEFORE=$(stat -c %s fault.vtkhdf)
# Simulate Phase B clobbering Phase A:
echo "phase B clobber padded out to a different byte count xxx" > fault.vtkhdf
FAULT_AFTER=$(stat -c %s fault.vtkhdf)
[ "${FAULT_BEFORE}" = "${FAULT_AFTER}" ] && echo "OLD CHECK: PASS (wrong)" \
                                          || echo "NEW CHECK: FAIL (correct)"
```

---

### [R-104] MODERATE [tpv104_checkpoint.hpp:118–121,161–168,295–298] — `expected_Q_size = -1` sentinel allows silent bypass in production code

**Category:** ASSUMPTION (defensive-design hole introduced by round-7 R-002)

**Description:**
The round-7 R-002 fix correctly adds an `int expected_Q_size`
parameter to both Read overloads.  But it also adds a sentinel
escape: passing `-1` SKIPS the check.

```cpp
// io/tpv104_checkpoint.hpp:161-168
if (expected_Q_size >= 0)
{
   MFEM_VERIFY(Q_size == expected_Q_size, ...);
}
```

The header doc comment justifies the sentinel as "Unit-test callers
that round-trip a synthetic Q (no live mesh) pass -1 to skip the
check."  But the unit tests DO pass the actual size — Sub-test 1 at
line 157 passes `expected_Q_size=Qsize` (= 27); Sub-test 9 at line
484 passes `expected_Q_size=Qsize` (= 18).  No caller uses `-1`.  The
sentinel exists only as a footgun: a future driver developer who
forgets to pass `Q.Size()` (e.g., passes `0` or refactors to default
the arg to `-1`) silently re-enables the original R-002 bug.

The round-7 reviewer's R-002 suggested fix explicitly called for a
REQUIRED parameter: "Extend the signature of `ReadTpv104Checkpoint(...)`
... to accept an extra parameter `int expected_Q_size`."  Sentinel
defaults defeat that intent.

**Trigger:**
A future driver edit passes `-1`, or passes `expected_Q_size` from a
variable that ends up zero or negative (e.g., uninitialised `int`).

**Actual behavior:**
With `-1`, no MFEM_VERIFY.  Wrong-mesh restart silently accepted —
the exact bug R-002 was meant to close.

**Expected behavior:**
No sentinel.  Require every caller to pass a non-negative size.

**Suggested fix:**
Remove the sentinel branch; assert that the parameter is non-negative
and the file's size matches:

```diff
@@ io/tpv104_checkpoint.hpp:161-168
    int Q_size = 0;
    read_tag("Q_size"); in >> Q_size;
-   if (expected_Q_size >= 0)
-   {
-      MFEM_VERIFY(Q_size == expected_Q_size,
-                  "ReadTpv104Checkpoint: Q size mismatch: file has "
-                  << Q_size << " doubles but the current driver expects "
-                  << expected_Q_size << " (NUM_STATE * ndof_total). "
-                  "Different mesh, polynomial order, or partition?");
-   }
+   MFEM_VERIFY(expected_Q_size >= 0,
+               "ReadTpv104Checkpoint: caller passed expected_Q_size="
+               << expected_Q_size << " (must be >= 0; this parameter "
+               "is REQUIRED to gate wrong-mesh restart per R-002).");
+   MFEM_VERIFY(Q_size == expected_Q_size,
+               "ReadTpv104Checkpoint: Q size mismatch: file has "
+               << Q_size << " doubles but the current driver expects "
+               << expected_Q_size << " (NUM_STATE * ndof_total). "
+               "Different mesh, polynomial order, or partition?");
    Q.SetSize(Q_size);
```

Update the doc-comment to remove the sentinel rationale:

```diff
@@ io/tpv104_checkpoint.hpp:117-122
-/// Single body for ReadTpv104Checkpoint (raw rank/size).  Both public
-/// overloads forward here.  `expected_Q_size` enforces R-002: if it is
-/// >= 0, the file's Q_size must equal it; if < 0, the check is skipped
-/// (legacy unit-test callers that want to round-trip a synthetic Q
-/// pass -1).  All production callers pass the driver's
-/// NUM_STATE * ndof_total via Q.Size() before this function resizes Q.
+/// Single body for ReadTpv104Checkpoint (raw rank/size).  Both public
+/// overloads forward here.  `expected_Q_size` enforces R-002:
+/// MUST be >= 0; the file's Q_size MUST equal it or the function
+/// aborts via MFEM_VERIFY.  All callers pass the live wave-field
+/// size (NUM_STATE * ndof_total) so wrong-mesh restart fails loudly.
```

Same edit to the public-overload comment at lines 295-298.

Unit-test changes — no source changes needed; Sub-test 1 and Sub-test
9 already pass the right size.

**Test case:**
Sub-test 6 grep that no `>= 0` sentinel branch remains:
```cpp
// R-104: no sentinel branch — Q_size check must be unconditional.
MFEM_VERIFY(
   hsrc.find("if (expected_Q_size >= 0)") == std::string::npos,
   "R-104: ReadTpv104CheckpointImpl must not have a sentinel "
   "branch that skips the Q-size check");
```

---

### [R-105] LOW [test_tpv104_checkpoint.cpp:9–14] — Header lists Sub-tests 2 and 3 that are not implemented

**Category:** QUALITY (documentation drift, pre-existing but inherited by R-008 extension)

**Description:**
The header comment lists 8 sub-tests (1, 2, 3, 4, 5, 6, 7, 9 — note 8
is skipped, which is intentional).  But `Subtest2_*` and `Subtest3_*`
functions do not exist; `main()` calls only 1, 4, 5, 6, 7, 9.  The
header header lies about coverage.

This is pre-existing (sub-tests 2 and 3 were never implemented), but
the round-7 R-008 extension added 7 and 9 to the same lying list
without removing or marking the absent 2/3 entries.

**Trigger:**
Read the header.

**Suggested fix:**
Either implement 2 and 3 (wrong-rank and wrong-num-ranks guards —
same SKIP pattern as 4/5 since MFEM_VERIFY isn't catchable), or mark
them explicitly as NOT IMPLEMENTED:

```diff
@@ tests/unit/test_tpv104_checkpoint.cpp:11-14
-//   2. Wrong-rank guard: a checkpoint written for rank=R must NOT
-//      be readable as rank=R+1.
-//   3. Wrong-num-ranks guard: a checkpoint written under N ranks
-//      must NOT be readable under M != N ranks.
+//   2. (NOT IMPLEMENTED) Wrong-rank guard — same MFEM_VERIFY
+//      non-catchable caveat as Sub-tests 4/5; covered by source-level
+//      MFEM_VERIFY in ReadTpv104CheckpointImpl.
+//   3. (NOT IMPLEMENTED) Wrong-num-ranks guard — same caveat as #2.
```

---

### [R-106] LOW [tpv104_driver.cpp:500] — Block comment still references `bulk.vtkhdf`

**Category:** QUALITY (documentation drift — round-7 R-003 fixed the user-facing error message but missed the surrounding block comment)

**Description:**
The block comment above the safety-check block reads:

```cpp
// Without this gate a user who
// re-uses the same --output-dir across a restart would silently
// lose Phase A's fault.vtkhdf / volume.vtkhdf / bulk.vtkhdf /
// station files / checkpoint files.
```

The actual user-facing error message (lines 528–531, fixed by R-003)
correctly says `ParaView_bulk/volume.vtkhdf` and `*_station_*.dat`.
The block comment was not updated, so it still says `bulk.vtkhdf`
(which TPV104 does not write — per Phase 6.4 R-310 the secondary
collection writes `<output_dir>/ParaView_bulk/volume.vtkhdf`).

**Suggested fix:**
```diff
@@ drivers/tpv104_driver.cpp:498-501
    // Without this gate a user who
    // re-uses the same --output-dir across a restart would silently
-   // lose Phase A's fault.vtkhdf / volume.vtkhdf / bulk.vtkhdf /
-   // station files / checkpoint files.  TPV104's TWO ParaView
+   // lose Phase A's fault.vtkhdf / volume.vtkhdf /
+   // ParaView_bulk/volume.vtkhdf / *_station_*.dat /
+   // *_checkpoint_r*.txt files.  TPV104's TWO ParaView
    // collections (pv_out + pv_bulk_out) make the clobber risk
```

---

### [R-107] LOW [petsc_ts_restart_plan_2026-05-16.md:897] — "6 sub-tests" outdated count

**Category:** QUALITY (plan doc inventory drift)

**Description:**
The plan doc R-006 entry says:

> `tests/unit/test_tpv104_checkpoint.cpp` — 6 sub-tests covering
> round-trip, wrong-format guard, and driver-grep.

After round 7 the test has Sub-tests 1, 4, 5, 6, 7, 9 (6 implemented,
2 documented-but-not).  R-008 added Sub-test 9 (the cross-overload
round-trip) — a substantively new piece of coverage worth naming in
the inventory.

**Suggested fix:**
```diff
@@ debug_document/.../petsc_ts_restart_plan_2026-05-16.md:897
-  - `tests/unit/test_tpv104_checkpoint.cpp` — 6 sub-tests covering round-trip, wrong-format guard, and driver-grep.
+  - `tests/unit/test_tpv104_checkpoint.cpp` — 6 implemented sub-tests covering V1 round-trip (#1), wrong-format/dof-size/Q-size guards (#4/#5/#7 — SKIP, MFEM_VERIFY non-catchable), driver+header source-grep covering R-001..R-007 (#6), and cross-overload byte-identity round-trip (#9, gates R-007 dedup).
```

---

### [R-108] POSSIBLE MODERATE [tpv104_checkpoint.hpp:163-168 + test_tpv104_checkpoint.cpp:Sub-test 7] — R-002 runtime behavior not tested; only static grep coverage

**Category:** EDGE_CASE (test coverage gap)

**Description:**
Sub-test 6's R-002 grep verifies the SOURCE CODE contains the string
`"Q size mismatch"`.  But it does not verify the runtime check FIRES.
If a future refactor introduces an off-by-one in the comparison
(e.g., `Q_size == expected_Q_size + 1`), the grep still passes but
the check is broken.

Sub-test 7 is SKIPed for the same reason as Sub-tests 4/5 —
MFEM_VERIFY isn't catchable.  But BP5's Sub-test 14 demonstrates an
opt-in subprocess-based runtime test pattern (`SEAS_TEST_RUNTIME_SAFETY_CHECK=1`
env var) that could provide actual runtime coverage.

**Trigger:**
Refactor `if (expected_Q_size >= 0) MFEM_VERIFY(Q_size == expected_Q_size,...)`
to `MFEM_VERIFY(Q_size != expected_Q_size, ...)` (inverted operator).
All current tests still PASS — the grep finds the string, Sub-tests
1 and 9 pass because the size match would invert into a fail.  Wait
— Sub-tests 1 and 9 would actually catch this because they pass the
CORRECT size which would now MFEM_VERIFY-fail (`Q_size != expected_Q_size`
is false → abort).  OK so this specific refactor would be caught.

But a subtler refactor — `MFEM_VERIFY(Q_size >= 0, ...)` (drop the
== comparison) — would silently pass all current tests because the
grep finds "Q size mismatch", Sub-tests 1/9 trivially pass any
Q_size >= 0, and the runtime check would never fire on wrong sizes.

**Suggested fix:**
Add an opt-in subprocess test (mirror of BP5's Sub-test 14):

```cpp
static void Subtest8_WrongQSizeRuntime()
{
   std::cout << "\n--- Sub-test 8: wrong-Q-size runtime guard (opt-in) ---\n";
   if (!std::getenv("SEAS_TEST_RUNTIME_QSIZE_CHECK"))
   {
      std::cout << "  INFO: SKIP — opt-in via SEAS_TEST_RUNTIME_QSIZE_CHECK=1\n";
      return;
   }
   // Write a checkpoint with Q_size = 27, then fork+exec self with
   // a magic argv that re-runs READ via raw-MPI overload with
   // expected_Q_size = 100.  The child MUST exit non-zero (MFEM_VERIFY
   // abort).  Parent reaps the child and asserts exit != 0.
   // (Pattern lifted from test_bp5_petsc_ts_restart.cpp Sub-test 14.)
}
```

**POSSIBLE flag rationale:** I cannot prove the current grep+SKIP
combination is insufficient — Sub-tests 1 and 9 do exercise the
correct-size path.  But the gap exists for INVERTED OR ABSENT
runtime check refactors.  Downgrade to LOW if you accept the risk.

**Test case:**
See suggested fix.  Manually invoke with
`SEAS_TEST_RUNTIME_QSIZE_CHECK=1 ./seas_test_tpv104_checkpoint`
on a build with the subprocess pattern implemented.

---

## Summary

- Critical issues: **0**
- Moderate issues: **4** (R-101, R-102, R-103, R-104)
- Low issues: **3** (R-105, R-106, R-107)
- POSSIBLE issues: **1** (R-108)
- Plan compliance: **FULL** for R-001..R-009; **PARTIAL** for round-8
  audit (R-103 reveals a long-standing logic flaw in BOTH BP5 and
  TPV104 sbatches that pre-dates round 7; R-101 documentation drift
  partly caused by round 7's incomplete sweep).
- Verdict: **PASS WITH FIXES** — no Critical issues; round 7's
  fixes for R-001..R-009 are all correctly applied and tested.
  The MODERATE findings are documentation drift (R-101, R-102) and
  pre-existing logic flaws now exposed by the round-7 improvements
  (R-103, R-104).  None block cluster submission, but each one is a
  silent correctness gap that should be closed before the BP5/TPV104
  pair is treated as fully verified.

## Answer to the three review focus questions (still applicable from round 7)

1. **Does restart work in both BP5 and TPV104?**
   - BP5: YES (V2 PETSc TS restart, 94/94 unit-test PASS).
   - TPV104: YES (V1 checkpoint, 139/139 unit-test PASS, driver builds
     clean, cross-overload Sub-test 9 confirms R-007 dedup forwards
     correctly).  Open: R-104 sentinel hole; R-006 secondary-collection
     limitation documented + warned.

2. **Do all saved quantities continue to be saved in a separate file?**
   - LANDING is correctly separated (driver safety check enforces; both
     ParaView collections route by `output_dir`).
   - VERIFICATION of separation: station files now correctly checked
     (Validation #5/#6, #7 md5sum unchanged).  fault.vtkhdf check is
     INERT (R-103) — same flaw in BP5 sbatch.

3. **Is the sbatch job created for cluster verification with guards?**
   - TPV104: YES — `tpv104_restart_test_v1_dev_2hr.sbatch`, dev queue,
     2h, 8N×400r, 7 acceptance checks (exit codes 3–9).  Validations
     #5/#6/#7 now actually exercise station-file isolation post R-004.
     But: Validation #3 is inert (R-103); 5 stale comments mislead
     (R-101); job name + log path strings still say "scaffolding"
     (R-101).

## Unreviewed Areas

- The actual checkpoint round-trip on a real TPV104 mesh + 8N parallel
  decomposition (only unit-tested with synthetic 18- and 27-DOF
  vectors at rank=0/size=1).  Cluster-only verification.
- Whether `seas_driver.cpp` (BP5 production) needs the same restart
  hooks (per `seas/CLAUDE.md` Phase 4 deferred deviation, it doesn't
  use `ParaViewOutput`).  Same out-of-scope conclusion as round 7.
- BP5 sbatch's Validation #8 (the equivalent of TPV104's broken
  Validation #3) — flagged at R-103, fix should be applied to BOTH
  sbatches in the next fix round.
- `make test-checkpoint` (BP5 V1 checkpoint test) has a pre-existing
  compile error in `AntiplaneDomainOperator::ComputeTractionDiagnostics`
  / `IsFirstStepDebugEnabled`.  Files unmodified by round 7; out of
  scope here.
