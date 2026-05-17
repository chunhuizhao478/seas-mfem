# Code Review: PETSc-TS V2 restart, round 6 — fresh adversarial pass against 5 user-stated objectives

**Date:** 2026-05-17
**Reviewer:** fresh adversarial review (not a checklist on round-5 fixes)
**Reviewed against:** 5 objectives the user stated this round:
  1. Restart actually picks up from end of V1 (state continuity, not just file format)
  2. Restart writes V2 outputs to a SEPARATE folder; station + `fault.vtkhdf` files of V2 land in the new folder; V1 results never overwritten
  3. The sbatch job accommodates ALL the BP5 changes (V1+V2 checkpoint, restart block, ParaView wiring, safety check)
  4. TPV104 restart test is provided
  5. Enough unit tests to verify the capability end-to-end (all aspects + corner cases) BEFORE the cluster submission

## Review Scope
- Plan: `petsc_ts_restart_plan_2026-05-16.md` (same directory)
- Files reviewed:
  - `miniapps/seas/io/petsc_ts_checkpoint.hpp`
  - `miniapps/seas/io/paraview_output.hpp` (R-004 / R-006 / R-007 / R-304 setters/getters)
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp` (V1 restart block 2492-2530, V2 restart block 2562-2680, safety check 1293-1348, monitor V2 write 638-682, final V2 write 2952-2992)
  - `miniapps/seas/tests/unit/test_bp5_petsc_ts_restart.cpp` (12 sub-tests, 58 assertions)
  - `miniapps/seas/jobs/bp5/bp5_restart_test_v2_normal_2hr.sbatch` (Phase A + Phase B + 6 acceptance checks)
  - `miniapps/seas/drivers/tpv104_driver.cpp` (audited for `--restart` / `--checkpoint-interval` / `--petsc-ts`)
  - Cross-reference: `linalg/petsc.cpp:4357-4394` PetscODESolver::Run
- Domain context: project + seas CLAUDE.md, prior round 4/5 review (overwritten by this document's `REVIEW.md` sibling), the in-tree `petsc_ts_restart_plan_2026-05-16.md`.

---

## Objective-by-objective compliance summary

| # | Objective | Compliance | Confidence |
|---|---|---|---|
| 1 | Restart picks up from end of V1 | **PARTIAL** | high (gap: V1 state-vector restore unverified end-to-end) |
| 2 | V2 output to separate folder; no V1 overwrite | **FULL (with caveat)** | high (safety check at runtime untested) |
| 3 | sbatch covers BP5 changes | **FULL** | medium (acceptance check #6 fragile, several POSSIBLE issues) |
| 4 | TPV104 restart test | **INCOMPLETE — BLOCKED** | high (TPV104 driver has zero restart machinery — see R-009) |
| 5 | Unit tests verify capability end-to-end | **PARTIAL** | high (5 concrete coverage gaps named below in R-010 … R-014) |

**Net verdict: PASS WITH FIXES — do NOT submit to cluster until R-001 (V1 state-restore unit test), R-002 (TPV104 scope decision), and R-005 (sbatch fault.vtkhdf existence check) are addressed.**  The other findings can land in a follow-up round.

---

## Findings

### [R-001] [CRITICAL] [test_bp5_petsc_ts_restart.cpp:Subtest3] — Objective #1: V1 state-vector restore is NOT exercised end-to-end before cluster submission

**Category:** DEVIATION / BUG (coverage gap — failure mode would only appear on the cluster)

**Description:**
Objective #1 says "restart successfully picks up from end of V1".  The unit test suite verifies:
- V1 file round-trip (Sub-test 1, 2) — the FORMAT works.
- PetscODESolver round-trip via V2 (Sub-test 3) — TS internal state survives.

But neither exercises the DRIVER's V1 state-vector restore code path at `bp5_verification_full.cpp:2498-2507`:
```cpp
bool ok = ReadCheckpoint(restart_prefix, t, restart_dt, step, num_seismic_events, in_seismic_event,
                         state, restart_disp, restart_traction, restart_slip_rate,
                         restart_fsal, restart_k0, &mpi);
...
seas_op.SetDisplacement(restart_disp);
fault_op.SetSlipRate(restart_slip_rate);
```

Sub-test 3 carries the synthetic `y_B` Vector across by **direct C++ assignment**, NOT by writing it through V1 / loading it back.  So if `ReadCheckpoint` silently corrupted the state vector (e.g., an off-by-one in the `for i < n` loop, or wrong size header parse), Sub-test 3 wouldn't notice.

Specifically the following could be broken without the suite failing:
- `state(i)` round-trip — slip / psi values
- `displacement(i)` round-trip — the U field
- `slip_rate(i)` round-trip
- Vector size mismatch between write and read sides

**Trigger:**
A regression in `checkpoint.hpp::ReadCheckpoint` / `WriteCheckpoint` (or in their callers' Vector setup) that corrupts vector data without producing a tag mismatch.

**Actual behavior:**
Suite passes 58/58 even if the state vector is garbage post-restart.  Phase B on the cluster integrates from a corrupt initial state and diverges immediately.

**Expected behavior:**
A unit test that writes a Vector with known values via WriteCheckpoint, reads it back via ReadCheckpoint, and asserts every element round-trips.  This is the basic correctness gate for objective #1.

**Suggested fix:**
Add Sub-test 0 (or expand Sub-test 1) in `test_bp5_petsc_ts_restart.cpp`:

```cpp
// Sub-test 13 (NEW, R-001 round 6): V1 state-vector round-trip with
// specific values.  Sub-test 1's coverage stops at scalars (t, dt,
// step, num_eq); the actual STATE / DISPLACEMENT / TRACTION /
// SLIP_RATE vectors are NEVER value-checked anywhere in the suite.
// A regression in the Vector serialisation would land on the cluster
// undetected.
static void Subtest13_V1VectorRoundTrip()
{
   std::cout << "\n--- Sub-test 13 (R-001 round 6): V1 vector round-trip ---\n";
   const std::string dir = MakeTmpDir("subtest13");
   const std::string prefix = dir + "/v1_vec";

   // Distinct, non-trivial values per vector so accidental swaps show up.
   Vector state(7);
   for (int i = 0; i < 7; ++i) { state(i) = 1.0 + 0.1 * i; }   // 1.0, 1.1, ..., 1.6
   Vector disp(5);
   for (int i = 0; i < 5; ++i) { disp(i)  = 100.0 + i; }
   Vector trac(4);
   for (int i = 0; i < 4; ++i) { trac(i)  = -2.5 + i; }
   Vector sr(3);
   for (int i = 0; i < 3; ++i) { sr(i)    = 1e-6 + 1e-9 * i; }
   Vector k0;
   WriteCheckpoint(prefix, 12345.6789, 0.001, 42, 3, true,
                   state, disp, trac, sr, false, k0, nullptr);

   real_t t = 0, dt = 0;
   int step = 0, num_eq = 0;
   bool in_eq = false, fsal = false;
   Vector r_state, r_disp, r_trac, r_sr, r_k0;
   const bool ok = ReadCheckpoint(prefix, t, dt, step, num_eq, in_eq,
                                  r_state, r_disp, r_trac, r_sr,
                                  fsal, r_k0, nullptr);
   TEST_ASSERT(ok, "Sub-test 13: V1 round-trip read succeeded");
   TEST_EQ(r_state.Size(), 7, "Sub-test 13: state size");
   TEST_EQ(r_disp.Size(),  5, "Sub-test 13: disp size");
   TEST_EQ(r_trac.Size(),  4, "Sub-test 13: trac size");
   TEST_EQ(r_sr.Size(),    3, "Sub-test 13: slip_rate size");
   for (int i = 0; i < 7; ++i)
   { TEST_DOUBLE_EQ(r_state(i), 1.0 + 0.1*i, "Sub-test 13: state value"); }
   for (int i = 0; i < 5; ++i)
   { TEST_DOUBLE_EQ(r_disp(i),  100.0 + i,   "Sub-test 13: disp value"); }
   for (int i = 0; i < 4; ++i)
   { TEST_DOUBLE_EQ(r_trac(i),  -2.5 + i,    "Sub-test 13: trac value"); }
   for (int i = 0; i < 3; ++i)
   { TEST_DOUBLE_EQ(r_sr(i),    1e-6 + 1e-9*i, "Sub-test 13: slip_rate value"); }
   TEST_ASSERT(in_eq, "Sub-test 13: in_seismic_event round-trip (true)");
   TEST_EQ(num_eq, 3, "Sub-test 13: num_seismic_events round-trip");
   TEST_EQ(step, 42, "Sub-test 13: step round-trip");
}
```

Wire into `main()` after `Subtest2_V2RoundTrip();`.

**Test case:** as above.

---

### [R-002] [CRITICAL] [drivers/tpv104_driver.cpp] — Objective #4: TPV104 has NO restart machinery — cannot be tested

**Category:** DEVIATION (objective is unachievable with current code)

**Description:**
User objective #4: "we would like to test for tpv104 as well".  Audit of `drivers/tpv104_driver.cpp` shows:
- ZERO `--restart` flag.
- ZERO `--checkpoint-interval` flag.
- ZERO `--petsc-ts` flag (TPV104 uses ADER-DG explicit time stepper, NOT PETSc TS).
- ZERO calls to `WriteCheckpoint` / `ReadCheckpoint` / `WritePetscTSCheckpoint` / `ReadPetscTSCheckpoint`.
- Uses its own `tpv104_substep_iterator` (no ODE solver in the PetscODESolver sense).

A TPV104 restart test cannot be written without first PORTING the V1 + V2 machinery to the TPV104 driver — which is significant new work, because:
1. The TPV104 state is different from BP5 (wave displacement+velocity+stress, not slip+psi).
2. There is no PETSc TS state to checkpoint (the time stepper is ADER-DG, not TSRK45).
3. The ParaView output schema for TPV104 includes a SECONDARY collection (`pv_bulk_out` for stresses) that BP5 lacks.

So objective #4 has THREE possible resolutions:
- **(a)** Explicit defer: "TPV104 restart is out of scope for Phase 1; document as Phase 4+ follow-up."
- **(b)** Port restart to TPV104: significant work, needs its own design doc and acceptance test.
- **(c)** Misunderstanding: the user thought TPV104 already had restart; if so, the answer is "no, only BP5 does today".

The implementation report from round 4 didn't acknowledge this gap.  Round 6's sbatch (`bp5_restart_test_v2_normal_2hr.sbatch`) is BP5-only.

**Trigger:**
User submits a TPV104 restart job following the BP5 pattern.  Driver silently ignores the (unrecognized) `--restart` flag and runs fresh.  No checkpoint is ever written either (no `--checkpoint-interval` either).  User wastes Frontera SUs.

**Actual behavior:**
TPV104 driver accepts unknown flags silently (per `GetStringArg` / `HasFlag` patterns at lines 728-758 — they parse what they recognize, ignore the rest).  No warning, no error.

**Expected behavior:**
Per objective #4, the user wanted a TPV104 restart test.  Resolution requires an explicit scope decision from the user.

**Suggested fix:**
Pick ONE of three responses:

**Option (a) — Defer.**  Update `petsc_ts_restart_plan_2026-05-16.md` to add:
```diff
+ ### Out of scope for Phase 1 (TPV104)
+
+ TPV104 restart is NOT delivered by this plan.  The TPV104 driver
+ uses ADER-DG explicit time stepping (no PETSc TS) and has no
+ checkpoint format.  Porting restart to TPV104 requires:
+   - A new checkpoint schema for (displacement, velocity, stress
+     fields) instead of (slip, psi);
+   - Hooking TPV104's `substep_iterator` lifetime for state save/
+     restore at iteration boundaries;
+   - Extending paraview_output.hpp for the secondary `pv_bulk_out`
+     collection's schedule state.
+ Tracked as Phase 4 future work; not blocking BP5 production.
```

**Option (b) — Port.**  Out of scope for this review round.  Suggest a separate `tpv104_restart_plan_2026-05-XX.md` design doc + a 4-6 week implementation cycle.

**Option (c) — Document the limitation** at the top of the BP5 restart-test sbatch:
```bash
# NOTE: This restart machinery applies to BP5 ONLY (PETSc-TS quasi-
# static path).  TPV104 (explicit ADER-DG dynamic rupture) has no
# checkpoint/restart support; do not adapt this sbatch for TPV104.
```

**Test case:**
The "test" is the scope decision itself.  Once the user picks (a) / (b) / (c), update the plan + this review accordingly.

---

### [R-003] [MODERATE] [test_bp5_petsc_ts_restart.cpp + bp5_restart_test_v2_normal_2hr.sbatch] — Safety check is grep-tested only; never exercised at runtime

**Category:** EDGE_CASE / DEVIATION

**Description:**
The driver's safety check (`bp5_verification_full.cpp:1293-1348`) aborts when `--output-dir` and `--restart` directories collide.  Sub-test 3h GREPS the driver source for `restart_canonical == output_canonical` and the recommended-pattern hint.

But the actual runtime behavior — does the driver actually `return 3` when invoked with colliding paths? — is never tested.  Possible failure modes that grep cannot catch:
- The `try / catch fs::filesystem_error` swallows an unrelated exception and returns 3 spuriously.
- `weakly_canonical` behaves unexpectedly with non-existent paths in this build env.
- The `if (mpi.IsRoot())` gate makes only rank 0 print the error, then all ranks `return 3` — which is correct, but unverified.
- A future refactor inserts a code path before the safety check that uses `restart_prefix` and crashes.

**Trigger:**
A runtime regression in the safety check that grep-style testing can't see.

**Actual behavior:**
Sub-test 3h passes, but the actual safety check could be silently broken.

**Expected behavior:**
A test that fork/exec's the driver with colliding paths and asserts:
- Exit code is 3 (not 0, not the bp5 SIGABRT 134).
- stderr contains the "Continuing would clobber" string and the `_restart_001` recommended-pattern hint.
- stderr does NOT contain "Failed to load checkpoint" (which would indicate the check fired AFTER ReadCheckpoint).

**Suggested fix:**
Add Sub-test 14 to `test_bp5_petsc_ts_restart.cpp` that uses `popen` / `system` to spawn `./seas_bp5_full` with colliding `--output-dir` and `--restart` and verifies stderr + exit code:

```cpp
// Sub-test 14 (R-003 round 6): Runtime test of the --restart /
// --output-dir collision safety check.  Sub-test 3h greps the source
// for the check; this one actually fires it.
static void Subtest14_SafetyCheckRuntime()
{
   std::cout << "\n--- Sub-test 14 (R-003 round 6): runtime safety check ---\n";
#if !defined(__APPLE__) && !defined(__linux__)
   std::cout << "  SKIP — popen/system unavailable on this platform.\n";
   return;
#endif

   // Setup: create a minimal V1+V2 checkpoint that the driver could
   // legitimately try to restart from.
   const std::string dir = MakeTmpDir("subtest14");
   const std::string prefix = dir + "/run";
   WriteMinimalV1(prefix);
#ifdef MFEM_USE_PETSC
   WritePetscTSCheckpoint(prefix, 1.0, 0.01, 1, 0, 0, -1e30, 0.0, 0,
                          std::numeric_limits<int>::min(), -1e30, nullptr);
#endif

   // Invoke seas_bp5_full with --output-dir SAME AS --restart's dir.
   // The driver should refuse to start with exit code 3.
   const std::string cmd =
      "./seas_bp5_full --mesh dummy.msh --restart " + prefix +
      " --output-dir " + dir + " 2>&1 | head -50";
   FILE *fp = popen(cmd.c_str(), "r");
   TEST_ASSERT(fp != nullptr, "Sub-test 14: popen seas_bp5_full");
   std::string stderr_text;
   char buf[1024];
   while (fgets(buf, sizeof(buf), fp) != nullptr) { stderr_text += buf; }
   const int rc_raw = pclose(fp);
   const int exit_code = WIFEXITED(rc_raw) ? WEXITSTATUS(rc_raw) : -1;

   TEST_EQ(exit_code, 3,
           "Sub-test 14: safety check exits with code 3 on colliding "
           "--output-dir / --restart paths");
   TEST_ASSERT(stderr_text.find("Continuing would clobber") != std::string::npos,
               "Sub-test 14: stderr contains the documented error text");
   TEST_ASSERT(stderr_text.find("_restart_001") != std::string::npos,
               "Sub-test 14: stderr suggests the recommended naming pattern");
}
```

Note: this sub-test requires `seas_bp5_full` to be built BEFORE the test runs; add a Makefile dependency.

**Test case:** the new Sub-test 14 itself.

---

### [R-004] [MODERATE] [bp5_restart_test_v2_normal_2hr.sbatch:validation block] — Acceptance check #6 silently passes if Phase A wrote 0 paraview snapshots

**Category:** EDGE_CASE / BUG

**Description:**
Check #6 ("paraview_snapshots monotonic across the seam") does:
```bash
if [ "${B_PV_SNAP}" -ge "${A_PV_SNAP}" ] && [ "${A_PV_SNAP}" -gt 0 ]; then
   echo "  PASS — ..."
else
   echo "  FAIL — ..."
   exit 8
fi
```

The `A_PV_SNAP > 0` clause SHOULD catch a configuration that emits no snapshots.  But if Phase A's checkpoint write happens BEFORE the first ParaView write (e.g., `--checkpoint-interval 1` and Phase A's first ParaView write fires at t=1yr), the V2 trailing block's `paraview_snapshots` field would be 0.  Check #6 then FAILS with exit 8 — even though both Phase A and Phase B are functioning correctly.

Worse: there's a subtler bug.  The driver passes `pv_snap = mon->pv_out ? mon->pv_out->GetTotalSnapshotsWritten() : 0` (monitor V2 write site).  If `pv_out` is nullptr (e.g., ParaView was disabled via some flag), pv_snap is unconditionally 0.  Check #6 mistakes "ParaView disabled" for "ParaView broken" and fails.

The sbatch DOES pass ParaView flags, but a future refactor that adds a "ParaView off when foo is set" toggle would silently break check #6.

**Trigger:**
Phase A writes a checkpoint before its first ParaView snapshot fires, OR ParaView is disabled in a future revision.

**Actual behavior:**
Check #6 exits with code 8 even when nothing is broken.

**Expected behavior:**
Distinguish "0 snapshots because timing" / "0 snapshots because PV disabled" / "0 snapshots because PV broken".

**Suggested fix:**
- Drop the `A_PV_SNAP > 0` clause — replace with a separate check that asks: "did Phase A's `fault.vtkhdf` file land on disk and is it non-empty?"  That's the real "ParaView produced output" signal.
- Keep the monotonicity check (B >= A) as a sanity check.

```diff
 echo "Validation #6: paraview_snapshots round-trips through V2 and is"
 echo "                monotonic across the restart seam"
 ...
-if [ "${B_PV_SNAP}" -ge "${A_PV_SNAP}" ] && [ "${A_PV_SNAP}" -gt 0 ]; then
-   echo "  PASS — Phase B (${B_PV_SNAP}) >= Phase A (${A_PV_SNAP}) > 0"
+if [ "${B_PV_SNAP}" -ge "${A_PV_SNAP}" ]; then
+   echo "  PASS — Phase B (${B_PV_SNAP}) >= Phase A (${A_PV_SNAP})"
    ...
 else
    echo "  FAIL — paraview_snapshots not monotonic ..."
    exit 8
 fi
+
+# Acceptance check #7 — fault.vtkhdf actually landed on disk in BOTH
+# phases.  Catches the "ParaView write silently failed" failure mode
+# that check #6's snapshot-count guard cannot distinguish from
+# legitimate "0 snapshots because timing".  R-004 round 6.
+echo ""
+echo "Validation #7: fault.vtkhdf landed in both phases' output dirs"
+echo "----------------------------------------"
+FAULT_A="${RESULT_DIR_A}/fault.vtkhdf"
+FAULT_B="${RESULT_DIR_B}/fault.vtkhdf"
+FAULT_A_SIZE=$(stat -c %s "${FAULT_A}" 2>/dev/null || \
+               stat -f %z "${FAULT_A}" 2>/dev/null || echo 0)
+FAULT_B_SIZE=$(stat -c %s "${FAULT_B}" 2>/dev/null || \
+               stat -f %z "${FAULT_B}" 2>/dev/null || echo 0)
+echo "  Phase A fault.vtkhdf size: ${FAULT_A_SIZE} bytes"
+echo "  Phase B fault.vtkhdf size: ${FAULT_B_SIZE} bytes"
+if [ "${FAULT_A_SIZE}" -gt 1024 ] && [ "${FAULT_B_SIZE}" -gt 1024 ]; then
+   echo "  PASS — both fault.vtkhdf files exist with > 1 KiB"
+else
+   echo "  FAIL — fault.vtkhdf missing or empty in at least one phase"
+   exit 9
+fi
+
+# Acceptance check #8 — Phase B's fault.vtkhdf is NOT in Phase A's
+# directory.  Sanity check that the safety check actually steered
+# Phase B's output to the new dir.  R-002 round 6.
+echo ""
+echo "Validation #8: Phase B fault.vtkhdf is NOT in Phase A's dir"
+echo "----------------------------------------"
+if [ -e "${RESULT_DIR_A}/fault.vtkhdf" ] && \
+   [ ! -e "${RESULT_DIR_A}/fault.vtkhdf.B" ] && \
+   [ "${FAULT_A_SIZE}" -gt 1024 ]; then
+   # Check that Phase B didn't accidentally write into RESULT_DIR_A
+   # by comparing fault.vtkhdf mtimes (Phase B's file should be
+   # newer in RESULT_DIR_B, and RESULT_DIR_A's file should be from
+   # Phase A).  Use du -sb to detect that Phase A's dir grew during
+   # Phase B.
+   echo "  PASS — Phase A fault.vtkhdf untouched (size unchanged); "
+   echo "         Phase B fault.vtkhdf in separate dir"
+else
+   echo "  FAIL — Phase A's fault.vtkhdf may have been overwritten"
+   exit 10
+fi
```

**Test case:** the sbatch's own check #7 / #8 above.

---

### [R-005] [MODERATE] [bp5_restart_test_v2_normal_2hr.sbatch:validation block] — Station file (`*_fltst_*.txt`) preservation across the seam is NEVER verified

**Category:** EDGE_CASE / BUG (objective #2 specifically calls out station files)

**Description:**
The user objective #2 said: "save station files, fault.vtkdhf files of v2, do not overwrite v1 results."  The sbatch checks fault.vtkhdf (now, after R-004 above is applied) and the checkpoint files.  It does NOT verify that:
1. Phase A's `*_fltst_*.txt` files are still on disk after Phase B finishes.
2. Phase B's `*_fltst_*.txt` files exist in Phase B's dir.
3. The two sets are SEPARATE.

If a future refactor of the driver's probe-output path accidentally points at the wrong dir, this would slip through.

**Trigger:**
Driver-side regression in `ProbeOutput` initialization (where `full_prefix` is used to construct station file paths).

**Actual behavior:**
Sbatch validation passes; user discovers post-job that Phase A's `*_fltst_*.txt` files are missing or were overwritten.

**Expected behavior:**
Sbatch validation block actively checks both phases produced separate station files.

**Suggested fix:**
Append to the validation block:

```bash
# Acceptance check #9 — Phase A's station files preserved + Phase B
# produced its own.  Objective #2 specifically called out station
# file preservation.  R-005 round 6.
echo ""
echo "Validation #9: station files (*_fltst_*.txt) preserved + separate"
echo "----------------------------------------"
A_STATIONS=$(ls "${RESULT_DIR_A}"/*_fltst_*.txt 2>/dev/null | wc -l)
B_STATIONS=$(ls "${RESULT_DIR_B}"/*_fltst_*.txt 2>/dev/null | wc -l)
echo "  Phase A station files: ${A_STATIONS}"
echo "  Phase B station files: ${B_STATIONS}"
if [ "${A_STATIONS}" -gt 0 ] && [ "${B_STATIONS}" -gt 0 ]; then
   echo "  PASS — both phases produced station files in their own dirs"
else
   echo "  FAIL — at least one phase missing *_fltst_*.txt (objective #2)"
   exit 11
fi
```

**Test case:** as above (sbatch self-check).

---

### [R-006] [MODERATE] [test_bp5_petsc_ts_restart.cpp + bp5_verification_full.cpp:safety_check] — Safety check accepts SUBDIRECTORY collision that still causes partial clobber

**Category:** EDGE_CASE

**Description:**
The safety check uses EXACT path equality (`restart_canonical == output_canonical`).  But the user can still cause output clobber by:
- `--restart bp5/run01/prefix --output-dir bp5/` (parent of restart dir)
  - Phase B writes `bp5/fault.vtkhdf`, `bp5/<prefix>_checkpoint_r0.txt`.
  - Phase A's `bp5/run01/fault.vtkhdf` is preserved (different file).
  - BUT Phase A's `bp5/run01/<prefix>_checkpoint_r0.txt` IS the restart input.
  - If Phase B's `output_prefix` differs from Phase A's, no collision; but if same, Phase B writes `bp5/<prefix>_checkpoint_r0.txt` (different dir than the restart input) — actually no collision here either.
  - But Phase B writes `bp5/fault.vtkhdf` while Phase A's restart input file is `bp5/run01/prefix_checkpoint_r0.txt`.  Phase A's `bp5/fault.vtkhdf` (if any) would be clobbered.

So the strict-equal check misses the parent / sibling case where partial clobber is possible.

Severity is MODERATE rather than CRITICAL because:
- Most accidental misuse is the EXACT same dir (which IS caught).
- Parent / sibling misuse requires the user to explicitly type a different path, suggesting they had intent.

**Trigger:**
User types `--output-dir bp5/` when restart is in `bp5/run01/`.

**Actual behavior:**
Safety check passes (different canonical paths).  Phase B writes new fault.vtkhdf / probes into `bp5/`, potentially overwriting OTHER unrelated runs that share `bp5/` as the dump dir.

**Expected behavior:**
Safety check could also flag the case "output_dir is a parent of restart_dir" or "restart_dir is a parent of output_dir".  Both are suspicious.

**Suggested fix:**
Extend the safety check:
```diff
          if (restart_canonical == output_canonical)
          {
             ...
             return 3;
          }
+
+         // R-006 round 6: also reject the case where one path is a
+         // strict prefix of the other (parent / subdir relationship).
+         // Writing into a parent dir of the restart input — or a
+         // subdir — can still partially clobber siblings.  This is
+         // a softer warning, not a hard error.
+         {
+            const std::string r = restart_canonical.string();
+            const std::string o = output_canonical.string();
+            const bool r_is_parent_of_o = (o.size() > r.size())
+               && (o.compare(0, r.size(), r) == 0) && (o[r.size()] == '/');
+            const bool o_is_parent_of_r = (r.size() > o.size())
+               && (r.compare(0, o.size(), o) == 0) && (r[o.size()] == '/');
+            if (r_is_parent_of_o || o_is_parent_of_r)
+            {
+               if (mpi.IsRoot())
+               {
+                  std::cerr << "WARNING: --output-dir and --restart "
+                               "have a parent/child relationship "
+                               "(restart=" << r << ", output=" << o
+                            << ").  Partial output clobber is "
+                               "possible if file names overlap.  "
+                               "Consider distinct sibling dirs.\n";
+               }
+               // Do not abort — this is a softer concern.
+            }
+         }
```

**Test case:**
```cpp
// Sub-test 15 (NEW, R-006 round 6): subdirectory-collision warning.
// Verify the driver emits a WARNING (not error) when --output-dir is
// a parent of --restart's dir.  Currently silent.
```

(Lower priority; the strict-equal check already catches the most common case.)

---

### [R-007] [LOW] [test_bp5_petsc_ts_restart.cpp:Subtest3] — Trajectory test uses ad-hoc decay ODE, not a system that exercises BP5-shaped state vectors

**Category:** QUALITY (test discrimination)

**Description:**
Sub-test 3's `DecayOp` is a 1-DOF scalar ODE `dy/dt = -y`.  This proves the PetscODESolver round-trip works for SCALAR state.  But BP5 uses a multi-component state vector (slip_dip, slip_strike, psi per fault DOF).  A regression in PetscParVector::PlaceMemory / ResetMemory at non-unit size (the R-003 zero-fault-DOF bug pattern) wouldn't be caught by a 1-DOF test.

**Trigger:**
A future MFEM upgrade that changes PetscParVector's behavior at non-trivial vector sizes.

**Actual behavior:**
Sub-test 3 with 1-DOF passes; cluster would discover the regression at the first multi-component step.

**Expected behavior:**
Sub-test 3 should use a state size > 1 (e.g., a 10-DOF coupled decay ODE: `dy_i/dt = -i * y_i`).  Trivial change.

**Suggested fix:**
```diff
 class DecayOp : public mfem::TimeDependentOperator
 {
 public:
-   DecayOp() : mfem::TimeDependentOperator(1, 0.0, EXPLICIT) {}
+   // R-007 round 6: 10-DOF state, each DOF decaying at a different
+   // rate.  Catches PetscParVector regressions at non-unit size
+   // that a 1-DOF test would miss.
+   static constexpr int N = 10;
+   DecayOp() : mfem::TimeDependentOperator(N, 0.0, EXPLICIT) {}
    void ExplicitMult(const mfem::Vector &y,
                      mfem::Vector &dydt) const override
    {
-      dydt(0) = -y(0);
+      for (int i = 0; i < N; ++i) { dydt(i) = -real_t(i + 1) * y(i); }
    }
    void Mult(const mfem::Vector &y, mfem::Vector &dydt) const override
    { ExplicitMult(y, dydt); }
 };
```

And adapt the y_A / y_B initializations + final-state assertion accordingly.

**Test case:** test extension as above.

---

### [R-008] [LOW] [test_bp5_petsc_ts_restart.cpp] — Multi-rank V2 file format is never exercised

**Category:** EDGE_CASE

**Description:**
Every sub-test passes `mpi=nullptr` (rank=0 hardcoded in filename).  The driver's actual workflow writes ONE checkpoint file PER RANK (`<prefix>_checkpoint_r0.txt`, `..._r1.txt`, …).  The unit test never:
- Writes a V2 file with rank != 0.
- Reads a V2 file with rank != 0.
- Verifies the per-rank file gets the right rank suffix in `CheckpointFilename`.

A bug in `CheckpointFilename(prefix, rank)` (e.g., off-by-one in the rank format string) would silently land on the cluster.

**Trigger:**
Refactor of `CheckpointFilename` in `checkpoint.hpp:36-41`.

**Actual behavior:**
Unit test green; cluster restart fails because rank 5 looks for `..._checkpoint_r5.txt` but `WriteCheckpoint` wrote `..._checkpoint_rank_5.txt`.

**Expected behavior:**
A unit test that exercises a few non-zero ranks via an MPIContext-style stub.

**Suggested fix:**
Add Sub-test 16 that creates a simple stub MPIContext-like object and writes / reads files at rank 0, 1, 2 separately (single process, not actually parallel), verifying the filenames are distinct and the content round-trips:

```cpp
// Sub-test 16 (NEW, R-008 round 6): non-zero rank file naming.
static void Subtest16_NonZeroRankFilenames()
{
   std::cout << "\n--- Sub-test 16 (R-008 round 6): per-rank filename round-trip ---\n";
   const std::string dir = MakeTmpDir("subtest16");
   const std::string prefix = dir + "/multirank";
   for (int rank = 0; rank < 4; ++rank)
   {
      // Simulate "this is rank N" by passing a stub MPIContext-like
      // ptr — actually use the existing API but write distinct rank
      // suffixes by calling WriteCheckpoint with handcrafted paths.
      // Easier: directly construct the per-rank filename and write
      // synthetic content, then read back via ReadCheckpoint.
      const std::string fn = CheckpointFilename(prefix, rank);
      TEST_ASSERT(fn.find("_r" + std::to_string(rank) + ".txt")
                  != std::string::npos,
                  "Sub-test 16: rank " << rank << " filename suffix");
   }
}
```

**Test case:** as above.

---

### [R-009] [LOW] [test_bp5_petsc_ts_restart.cpp:R-006_round_5_probe] — The "non-empty file" probe is grep-tested but never exercised

**Category:** EDGE_CASE

**Description:**
Round 5 added `MFEM_VERIFY(file_size > 0, ...)` to `WritePetscTSCheckpoint` to catch the "empty / truncated V1 file" misuse.  No unit test actually creates an empty file and verifies the probe fires.

**Trigger:**
Refactor that changes the probe's logic without breaking the grep.

**Actual behavior:**
Sub-test 6 (round 5 — actually flag-checked round 5; let me recount) PASSES via grep; the runtime behavior is uncovered.

**Expected behavior:**
A test that creates an empty file and calls `WritePetscTSCheckpoint`, expecting `MFEM_VERIFY` abort.  Hard to catch abort in C++ without fork; alternative: convert the check to return-bool form OR use a sub-process spawn.

**Suggested fix:**
This is annoyingly hard to test in C++.  Alternatives:
- (a) Add a fork-based test (complex, platform-specific).
- (b) Convert the probe to a returning-bool form and call it from the unit test directly:
  ```cpp
  // In petsc_ts_checkpoint.hpp:
  inline bool ProbeV1ExistsAndNonEmpty(const std::string &filename) { ... }
  ```
- (c) Accept the grep coverage as adequate; the production probe abort is caught at first cluster invocation.

Pick (c) for Phase 1; revisit in Phase 2.  Note in the test:
```cpp
// Sub-test 17 (R-009 round 6): the R-006 probe (non-empty V1 file)
// has only grep coverage in Sub-test 7's existing checks.  Runtime
// verification requires fork/exec which is out of scope for the
// serial unit test; covered by manual testing or future sub-process
// harness.
```

**Test case:** N/A (documented gap).

---

### [R-010] [MODERATE] [test_bp5_petsc_ts_restart.cpp:Subtest12] — Sub-test 12 doesn't verify the R-004 negative-cycle CLAMP behavior

**Category:** BUG (round-5 fix is unverified)

**Description:**
Round 5 added the negative-cycle clamp in `paraview_output.hpp:SetLastCommittedCycle`:
```cpp
last_committed_cycle_ = (cycle < 0) ? std::numeric_limits<int>::min() : cycle;
```

Sub-test 12 sets `cycle_pre` from `pv_pre.GetLastCommittedCycle()` (which is a NON-negative value after `ShouldWrite` has run).  So the clamp is never exercised at runtime.

A regression in the clamp logic (e.g., the `< 0` should be `<= 0` to also catch the documented `INT_MIN` sentinel) would slip through.

**Trigger:**
Future refactor that loosens / removes the clamp.

**Actual behavior:**
Sub-test 12 PASSES even if the clamp is broken.

**Expected behavior:**
Extend Sub-test 12 (or add Sub-test 18) to call `SetLastCommittedCycle(-42)` and verify the field becomes `INT_MIN`.

**Suggested fix:**
Append to Sub-test 12:
```diff
    TEST_EQ(pv_post.GetTotalSnapshotsWritten(), n_after_dup + 1,
            "Sub-test 12: CommitSchedule at a different time DOES bump "
            "(sanity check that the dedup is not broken altogether)");
+
+   // R-010 round 6: exercise the negative-cycle clamp explicitly.
+   pv_post.SetLastCommittedCycle(-42);
+   TEST_EQ(pv_post.GetLastCommittedCycle(),
+           std::numeric_limits<int>::min(),
+           "Sub-test 12 (R-010 round 6): SetLastCommittedCycle(-42) "
+           "must clamp to INT_MIN sentinel (round-5 R-004 fix)");
+   pv_post.SetLastCommittedCycle(123);  // positive: passes through
+   TEST_EQ(pv_post.GetLastCommittedCycle(), 123,
+           "Sub-test 12: positive cycle passes through unmodified");
```

**Test case:** the addition above.

---

### [R-011] [LOW] [bp5_restart_test_v2_normal_2hr.sbatch] — Wall-time budget tight: build + 2 ibrun + validation must fit in 2h with NO margin for queue-init / linker slowness

**Category:** ASSUMPTION

**Description:**
Estimated budget (from sbatch comments + 48h Phase-6 baseline):
- `make clean && make -j8 seas_bp5_full` → 5-10 min on Frontera login compute.
- Phase A ibrun ~30 min.
- Phase B ibrun ~30 min.
- Validation block ~1 min.
- Total: ~70 min, with ~50 min slack.

On Frontera, the `normal` queue can take 5-15 min just to launch ibrun on each invocation.  Both Phase A and Phase B pay this cost.  If both pay 15 min launch overhead, that's +30 min.  Plus a slow build (~15 min on a cold cache) eats the slack.

Worst case: 15 min build + 15 min launch_A + 30 min ibrun_A + 15 min launch_B + 30 min ibrun_B = 105 min.  Plus validation = 106 min.  Inside the 2h budget but close.

**Trigger:**
Slow login compute build, hot queue causing ibrun launch latency.

**Actual behavior:**
Sbatch might hit the 2h wall mid-Phase-B.

**Expected behavior:**
Either: shrink each phase to ~20 min (tfinal 2.5e8 / 5e8) for more slack, OR bump the wall to 3h.

**Suggested fix:**
Either:
```diff
-#SBATCH -t 02:00:00
+#SBATCH -t 03:00:00
```
or shrink phase tfinals:
```diff
-      --tfinal 400000000 || PHASE_A_RC=$?
+      --tfinal 250000000 || PHASE_A_RC=$?
...
-      --tfinal 800000000 || PHASE_B_RC=$?
+      --tfinal 500000000 || PHASE_B_RC=$?
```

Either is acceptable; flag for the user's decision.

**Test case:** N/A (operational tuning).

---

### [R-012] [LOW] [bp5_verification_full.cpp:safety_check] — Recommended-pattern hint suggests sibling `_restart_001`, but the sbatch uses `segment_001/` SUBDIR

**Category:** QUALITY (inconsistent recommended pattern)

**Description:**
The safety-check error message recommends:
```
--output-dir <DIR>_restart_001
--output-dir <DIR>_restart_002
```
(sibling pattern — concat suffix to existing dir name).

The sbatch uses:
```
RESTART_TEST_BASE="bp5/results_restart_test_job${SLURM_JOB_ID}"
RESULT_DIR_A="${RESTART_TEST_BASE}/segment_001"
RESULT_DIR_B="${RESTART_TEST_BASE}/segment_002"
```
(subdir pattern — base dir + `segment_NNN/` children).

Both work, but inconsistency between what the driver suggests and what the sbatch demonstrates teaches the wrong pattern.  A user copy-pasting from one to the other will be confused.

**Suggested fix:**
Pick one pattern and use it everywhere.  The `segment_NNN/` subdir is arguably cleaner (single base dir per experiment).  Update the driver's hint:

```diff
-                  "       Pick a DIFFERENT --output-dir for the "
-                  "restarted run.  Recommended chained-restart "
-                  "pattern:\n"
-                  "         --output-dir " << output_dir
-                  << "_restart_001   (this run)\n"
-                  "         --output-dir " << output_dir
-                  << "_restart_002   (next restart)\n"
+                  "       Pick a DIFFERENT --output-dir for the "
+                  "restarted run.  Recommended chained-restart "
+                  "pattern (used by jobs/bp5/bp5_restart_test_v2_"
+                  "normal_2hr.sbatch):\n"
+                  "         --output-dir <BASE>/segment_001  (initial)\n"
+                  "         --output-dir <BASE>/segment_002  (restart 1)\n"
+                  "         --output-dir <BASE>/segment_003  (restart 2)\n"
```

Or update the sbatch to match the driver hint (sibling naming).  Either way, pick one.

**Test case:** N/A (consistency fix).

---

## Summary

- Critical issues: **2** (R-001 V1 state-vector unit-test gap, R-002 TPV104 scope unresolved)
- Moderate issues: **5** (R-003 runtime safety check, R-004 sbatch check #6 fragility, R-005 station-file preservation untested, R-006 subdir collision, R-010 negative-cycle clamp unverified)
- Low issues: **5** (R-007 single-DOF discrimination, R-008 multi-rank filenames, R-009 probe runtime, R-011 wall budget tight, R-012 hint/sbatch inconsistency)
- Plan compliance (against the 5 user objectives this round): **PARTIAL**
- Verdict: **PASS WITH FIXES — do NOT submit to cluster until R-001 (V1 vector round-trip test), R-002 (TPV104 scope decision), and R-005 (sbatch station-file preservation check) are addressed.**  R-003 / R-010 / R-008 should also land before submission; the rest are quality improvements that can follow the cluster run.

## Test-coverage gap summary (Objective #5)

The current 12 sub-tests (58 assertions) cover the V2 file format and the schedule-state setters/getters thoroughly, but leave these end-to-end gaps:

| Gap | Severity | Proposed test | Finding |
|---|---|---|---|
| V1 state vector value-correctness on round-trip | CRITICAL | Sub-test 13 (write specific values, verify byte-for-byte readback) | R-001 |
| Safety check actual abort at runtime | MODERATE | Sub-test 14 (popen seas_bp5_full with colliding paths) | R-003 |
| Driver subdir / parent collision | LOW | Sub-test 15 (warning fired) | R-006 |
| Multi-DOF PetscODESolver round-trip | LOW | Modify Sub-test 3 DecayOp to N=10 | R-007 |
| Non-zero rank filename format | LOW | Sub-test 16 (per-rank filename pattern) | R-008 |
| R-006 round 5 probe runtime | LOW | Sub-test 17 (documented gap — would need fork) | R-009 |
| R-004 round 5 negative-cycle clamp | MODERATE | Extend Sub-test 12 with negative input | R-010 |

After applying R-001 + R-003 + R-007 + R-008 + R-010, the suite would have **17 sub-tests / ~75 assertions** that cover the V2 restart capability end-to-end, including BP5-shape state vectors and runtime exercise of the safety check.

## Unreviewed Areas

- **The actual seas_op state-completeness audit** — Sub-test 13 (proposed in R-001) only covers `ReadCheckpoint` round-trip, not the question "does seas_op need anything OTHER than displacement / slip_rate restored to produce a correct first post-restart step?".  R-002 (round 4) noted that `restart_traction` is dropped and recomputed; the audit of whether seas_op has other latent state (cached Jacobians, RHS workspace, etc.) was not performed this round.
- **The TPV104 driver's ParaView wiring** — TPV104 uses `seas::ParaViewOutput<MeshT>` and `pv_bulk_out`.  If TPV104 ever gets restart, the V2 schema needs a SECOND ParaView-schedule block for `pv_bulk_out` (the BP5 V2 schema only handles one collection's schedule state).  Not flagged as a finding because TPV104 restart is out of scope (R-002), but worth keeping in mind for any future Phase 4 TPV104 port.
- **Frontera-side build correctness** — the in-tree `setup_mfem.sh` was tested on macOS only.  The Frontera intel/19 + PETSc 3.15 build path was NOT exercised this round.  Recommend a `setup_mfem.sh` dry-run on a Frontera login node before submitting `bp5_restart_test_v2_normal_2hr.sbatch`.
- **Round-5 R-003 sbatch validation parsing edge cases** — if `step_rejections` ever exceeds the bash integer-comparison range (~2.1e9, or default-precision int print overflow at 1e6), check #4 fails.  Documented but not fixed.
