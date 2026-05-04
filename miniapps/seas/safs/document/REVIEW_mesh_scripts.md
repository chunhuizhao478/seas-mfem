# Code Review: CGAL mesh-generation scripts (2026-04-30)

## Review Scope

- **Plan:** `PLAN_cgal_corefine.md` (Phases 3–5).
- **Files reviewed:**
  - `mesh/run_two_crossing_2000m_cgal.sh` (Phase-3 canonical)
  - `mesh/run_two_crossing_2000m.sh`       (Phase-4 wrapper)
  - `mesh/run_all8_2000m_cgal.sh`          (Phase-5 canonical)
  - `mesh/run_all8_2000m.sh`               (Phase-5 wrapper)
  - `mesh/generate_safs_mesh.py`           (STL-discovery path only)
- **Domain context:** prior reviews in `REVIEW_phase01..02_v2_fix.md`,
  `REVIEW_phase3.md` (R-405 wiring), `REVIEW_autorefine_mode.md`
  (the new code path's separate findings).
- **Live verification:**
  - I read the four shell scripts and the relevant
    `generate_safs_mesh.py` STL-discovery code.  I did NOT run the
    all-8 pipeline end-to-end (would hit the documented Phase-5
    self-intersection blocker on `safs_mult_banning`).

I assumed at least 3 bugs and found 5 (2 MODERATE, 3 LOW).  None
of the findings are CRITICAL — Phase-3 + the Phase-4 wrapper produce
a correct 11/11 mesh on Mill Creek × SBMT-SAF as verified earlier;
the issues are inconsistencies between the Phase-3 and Phase-5
canonical scripts, plus UX / hygiene issues in the wrappers.

---

## Findings

### [R-601] [MODERATE] [run_all8_2000m_cgal.sh:148–157] — validator-failure handler not propagated from Phase-3 (R-405) to Phase-5

**Category:** DEVIATION (between two near-identical scripts)

**Description:**
`run_two_crossing_2000m_cgal.sh:117–136` (Phase-3, R-405 applied)
catches a non-zero exit from `validate_msh.py`, dumps the full
report to stderr, and exits 1.  On success it tails the report
summary so the user sees pass/fail counts without opening the file.

`run_all8_2000m_cgal.sh:148–157` (Phase-5, R-405 not applied) just
calls `validate_msh.py` with the report path.  `set -euo pipefail`
on line 32 propagates the non-zero exit, so the script aborts —
but with no report dump on stderr.  The user sees:

```
==> 6/7 validate_msh
[corefine_faults]: ... (some opaque output) ...
$ # script exited 1, no clue which validator check failed
```

This is the exact silent-pass-through diagnostic gap R-405 was
written to close.  The all-8 build is the harder case (more
faults, more checks at risk); it should at LEAST match the
Phase-3 script's diagnostics.

**Trigger:** Any all-8 invocation where validate_msh.py reports
a check failure.

**Actual behavior:** Script aborts with no diagnostic on which
check failed.

**Expected behavior:** Same diagnostic surface as Phase-3: full
report dumped to stderr on FAIL, summary tail on success.

**Suggested fix:**

```diff
@@ run_all8_2000m_cgal.sh:148–157
 echo "==> 6/7 validate_msh"
+# R-405 (parity with Phase-3): surface validator FAIL by dumping
+# the full report to stderr.  Even on success, print the summary
+# tail so pass/fail counts are visible without opening the report.
 python validate_msh.py \
     --msh "$OUTDIR/output/safs_all8_${RES}m.msh" \
     --transform-json "$OUTDIR/transform.json" \
     --bbox-json "$OUTDIR/bbox.json" \
     --provenance-json "$OUTDIR/output/fault_provenance.json" \
     --domain-box-json "$OUTDIR/output/domain_box.json" \
     --sizing-json "$OUTDIR/output/sizing.json" \
     "${INCL[@]/--include-fault/--expected-faults}" \
-    --report "$OUTDIR/output/validation_report.txt"
+    --report "$OUTDIR/output/validation_report.txt" || {
+        echo "" >&2
+        echo "validate_msh.py FAILED — full report follows:" >&2
+        if [ -f "$OUTDIR/output/validation_report.txt" ]; then
+            cat "$OUTDIR/output/validation_report.txt" >&2
+        else
+            echo "(no report file written; validator likely crashed before report write)" >&2
+        fi
+        exit 1
+    }
+if [ -f "$OUTDIR/output/validation_report.txt" ]; then
+    tail -3 "$OUTDIR/output/validation_report.txt"
+fi
```

**Test case:**
```bash
# Build a deliberately-broken all-8 mesh (e.g., delete a fault STL
# between step 3 and step 4 so generate_safs_mesh produces output
# that fails validator check 8 provenance partition).  Run the
# all-8 script.  Assert: exit code 1; stderr contains "validation
# report follows:" and the offending check ID.
bash run_all8_2000m_cgal.sh 2>&1 | grep -q "validation_report.txt FAILED — full report follows:"
```

---

### [R-602] [MODERATE] [run_two_crossing_2000m_cgal.sh:97 + run_all8_2000m_cgal.sh:128] — `--remesh-iters 3` is hardcoded; no env-var override even though every other corefine_faults parameter is overridable

**Category:** ASSUMPTION

**Description:**
Both scripts have:

```bash
"$COREFINE_BIN" \
    --in-stl-dir "$OUTDIR/stl_raw" \
    --out-stl-dir "$OUTDIR/stl_conformal" "${INCL[@]}" \
    --clearance-m "$CLEARANCE" \
    --target-edge-m "$TARGET_EDGE" \
    --remesh-iters 3
```

`CLEARANCE`, `TARGET_EDGE`, and indeed every other tunable in the
script come from env vars with sensible defaults
(`TWO_CLEARANCE`, `TWO_TARGET_EDGE`, etc.).  But `--remesh-iters`
is the literal string `3` — the user cannot override it without
editing the script.

This is inconsistent with the rest of the script and surprising
when comparing different remesh-iter counts during validator
tuning (e.g., to investigate the R-401 + R-403 tube-uniformity
findings from `REVIEW_phase3.md`, which call out remesh quality
as a knob).

**Trigger:** A user trying to A/B compare different remesh-iter
counts.

**Actual behavior:** Stuck on 3.

**Expected behavior:** Configurable via env, like every other
parameter.

**Suggested fix:**

```diff
@@ run_two_crossing_2000m_cgal.sh:47
 TARGET_EDGE="${TWO_TARGET_EDGE:-1000}"
+REMESH_ITERS="${TWO_REMESH_ITERS:-3}"
@@ run_two_crossing_2000m_cgal.sh:91–97
 echo "==> 3/7 corefine_faults (CGAL; target_edge=${TARGET_EDGE} m)"
 "$COREFINE_BIN" \
     --in-stl-dir "$OUTDIR/stl_raw" \
     --out-stl-dir "$OUTDIR/stl_conformal" "${INCL[@]}" \
     --clearance-m "$CLEARANCE" \
     --target-edge-m "$TARGET_EDGE" \
-    --remesh-iters 3
+    --remesh-iters "$REMESH_ITERS"
```

Same diff structure for `run_all8_2000m_cgal.sh` with
`SAFS_REMESH_ITERS`.

**Test case:**
```bash
# Override the env var; verify it propagates to the binary.
TWO_REMESH_ITERS=0 bash mesh/run_two_crossing_2000m_cgal.sh 2>&1 \
    | grep "remesh-iters" | grep -q "0"
```

---

### [R-603] [LOW] [run_two_crossing_2000m.sh + run_all8_2000m.sh] — wrappers don't enable `set -euo pipefail`

**Category:** QUALITY

**Description:**
Both wrappers:

```bash
#!/usr/bin/env bash
# ... comment ...
export TWO_OUTPUT_SUFFIX="${TWO_OUTPUT_SUFFIX-}"
exec "$(dirname "$0")/run_two_crossing_2000m_cgal.sh" "$@"
```

No `set -euo pipefail`.  In practice the only failure mode is
`exec` failing if the inner script doesn't exist (unlikely; the
canonical script lives next to the wrapper).  But strict mode
is the standard convention for production shell scripts in this
repo (every other shell script in `mesh/` has it).

**Trigger:** N/A under normal operation.

**Actual behavior:** Wrappers silently inherit the inner script's
mode.

**Expected behavior:** Wrappers explicitly set strict mode.

**Suggested fix:**

```diff
 #!/usr/bin/env bash
 # ... comment ...
+set -euo pipefail
 export TWO_OUTPUT_SUFFIX="${TWO_OUTPUT_SUFFIX-}"
 exec "$(dirname "$0")/run_two_crossing_2000m_cgal.sh" "$@"
```

Same diff for `run_all8_2000m.sh`.

---

### [R-604] [LOW] [run_two_crossing_2000m_cgal.sh:130–135] — `cat <report>` on validator failure doesn't guard against the report file being missing

**Category:** EDGE_CASE

**Description:**

```bash
python validate_msh.py ... --report "$OUTDIR/output/validation_report.txt" || {
    echo "" >&2
    echo "validate_msh.py FAILED — full report follows:" >&2
    cat "$OUTDIR/output/validation_report.txt" >&2
    exit 1
}
```

If `validate_msh.py` crashes BEFORE writing the report (e.g.,
import error, argument-parsing exception, or kill signal), the
report file does not exist.  `cat <missing>` then emits
"cat: <path>: No such file or directory" to stderr — confusing
diagnostic.  Worse, on `set -e`, the `cat` failure inside the `||`
block is suppressed (we're in a `||` consequent), so the script
continues to `exit 1` — net behaviour is OK, but the user sees a
misleading "follows:" header followed by a missing-file message
instead of the actual error from the validator.

**Trigger:** Validator crashes before writing the report (rare).

**Suggested fix:**

```diff
 python validate_msh.py ... --report "$OUTDIR/output/validation_report.txt" || {
     echo "" >&2
     echo "validate_msh.py FAILED — full report follows:" >&2
-    cat "$OUTDIR/output/validation_report.txt" >&2
+    if [ -f "$OUTDIR/output/validation_report.txt" ]; then
+        cat "$OUTDIR/output/validation_report.txt" >&2
+    else
+        echo "  (no report file written; validate_msh likely "
+        echo "   crashed before report write — re-run with the "
+        echo "   stderr captured to see the real error)" >&2
+    fi
     exit 1
 }
```

**Test case:**
```bash
# Force validate_msh.py to crash before writing the report by
# passing a non-existent --msh path.
TWO_FAULTS_OVERRIDE=... bash mesh/run_two_crossing_2000m_cgal.sh 2>&1 \
    | grep -q "no report file written"
```

---

### [R-605] [LOW] [run_all8_2000m_cgal.sh:42–51] — `FAULTS` array hardcoded; user must edit script source to skip at-risk faults

**Category:** QUALITY (UX)

**Description:**
The all-8 script ships with all 8 faults pre-listed in `FAULTS=(...)`,
and the at-risk preflight (R-402a) tells the user to "edit the
FAULTS=(...) array near the top of this script" to skip them.

This forces the user to modify a checked-in script — the script
becomes a working file rather than a stable artifact.  Tracked
versions diverge between users.  The 2-fault script has the same
hard-coding (`FAULTS=(safs_sbmt_millcreek safs_sbmt_saf)`) but in
the 2-fault case both faults are documented-clean, so editing
isn't required.

For the all-8 case, an env-var override `SAFS_FAULTS` (space-
separated) would let the user run a 6-fault subset (excluding
`safs_mult_banning` + `safs_pmfz_pinto`) without script
modification.

**Trigger:** A user trying to skip the at-risk faults.

**Suggested fix:**

```diff
@@ run_all8_2000m_cgal.sh:42–51
-FAULTS=(
-    safs_coav_missioncreek
-    safs_mjvs_saf
-    safs_mult_banning
-    safs_mult_ssaf_banning
-    safs_pmfz_pinto
-    safs_sbmt_millcreek
-    safs_sbmt_missioncreek
-    safs_sbmt_saf
-)
+# Default to all 8 CFM faults at RES=2000.  Override with
+# `SAFS_FAULTS="<space-separated>"` to run a subset (e.g., to
+# skip the at-risk safs_mult_banning + safs_pmfz_pinto without
+# editing this script).
+if [ -n "${SAFS_FAULTS-}" ]; then
+    read -r -a FAULTS <<< "$SAFS_FAULTS"
+else
+    FAULTS=(
+        safs_coav_missioncreek
+        safs_mjvs_saf
+        safs_mult_banning
+        safs_mult_ssaf_banning
+        safs_pmfz_pinto
+        safs_sbmt_millcreek
+        safs_sbmt_missioncreek
+        safs_sbmt_saf
+    )
+fi
```

The R-402a preflight would also benefit from updating the
"To skip them, edit ..." message to mention the env var.

---

## Summary

- Critical issues: **0**
- Moderate issues: **2** (R-601 missing validate-fail handler in
  Phase-5; R-602 `--remesh-iters` hardcoded)
- Low issues: **3** (R-603 wrappers no strict-mode; R-604 cat-on-
  fail no file guard; R-605 FAULTS not env-overridable)
- Plan compliance: **PARTIAL**.
  - Phase 3 (`run_two_crossing_2000m_cgal.sh`): FULL compliance,
    11/11 PASS verified earlier.
  - Phase 4 wrapper: FULL compliance.
  - Phase 5 canonical: PARTIAL — the script exists and has the
    correct structure, but R-601 means the Phase-4 R-405 fix
    didn't propagate, and R-605 forces the user to edit the
    script to bypass the at-risk default.
- Verdict: **PASS WITH FIXES**.  The pipeline produces correct
  output on the Phase-3 happy path; findings are inconsistencies
  between the two canonical scripts plus wrapper hygiene.

## Why no CRITICAL findings?

The previous reviews (autorefine_mode, Phase-1/2 fix rounds) have
already cleaned up the C++ binary's correctness gaps.  The shell
scripts here are thin orchestrators that pass the right env-var
defaults and CLI flags to known-good Python + C++ binaries — most
of the reach into incorrect output requires those underlying
components to misbehave, which is covered by the other reviews.

The only ways the SHELL scripts themselves can produce a wrong
mesh are:
- silently changing the at-risk-fault default (R-605 makes that
  hard);
- silently skipping the validator (no script does this — both
  invoke validate_msh.py);
- silently using a stale STL from a prior run.  Verified
  `generate_safs_mesh.py` looks up STLs by EXACT short-name
  (lines 191, 298), not by glob — so stale files for unrelated
  faults are ignored even if they linger from a prior run.  Not a
  finding.

## Suggested fix order

1. **R-601** (Phase-5 validate-fail handler) — copy-paste from
   Phase-3, ~10 lines.
2. **R-602** (`--remesh-iters` env override) — three-line change
   in each script.
3. **R-603, R-604, R-605** — quality cleanups; can land
   together.

## Unreviewed Areas

- **`generate_safs_mesh.py`** is invoked from the scripts but its
  internal correctness is out of scope here (it's a Python script
  the plan says is "deliberately untouched" in Phases 0–3).  See
  `REVIEW_phase3.md` R-401, R-403, R-404, R-406 for the deferred
  size-field findings.
- **`validate_msh.py`** internal logic is out of scope; deferred
  R-401 + R-402 from `REVIEW_phase3.md` cover the gate-tightening.
- **Behaviour on macOS bash 3.2 vs Linux bash 5.x.**  The scripts
  use array syntax that works in both, but I did not verify the
  `set -u` interaction with empty array expansion under bash 3.2.
  Most users run via conda env which bundles bash 4.x or higher;
  not flagging as a finding.
- **Concurrency: two simultaneous invocations of the same script.**
  Both write to the same `$OUTDIR`; results are racey.  Out of
  scope for a single-user mesh-build script.
- **The `audit_ts_quality.py` step's CSV is never read by any
  later step.**  Informational only; not a finding.
