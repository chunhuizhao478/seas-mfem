# Code Review: Caliper perfgraph PLAN (fresh review, 2026-05-29)

> **Note on location/name:** the `/code-review` skill specifies `REVIEW.md` in
> the project root, but `miniapps/seas/REVIEW.md` already exists (phase 4–11
> reviews) and would be clobbered. This review is therefore the `*_check.md`
> companion to `caliper_perfgraph_plan.md`, matching the plan/check/fix
> convention used by the code-fix workflow. The `/code-fix` agent should read
> THIS file.

> **What is under review:** the PLAN document
> `document/caliper_perfgraph_dev/caliper_perfgraph_plan.md` — there is no
> implementation yet. Findings are defects *in the plan* (wrong/missing
> instructions, recipe gaps, completeness holes) that would mislead the
> implementer or yield a non-functional perfgraph.

## Review Scope
- Plan: `miniapps/seas/document/caliper_perfgraph_dev/caliper_perfgraph_plan.md`
- Reviewed against: `dynamic/wave_operator.{hpp,inl}`,
  `dynamic/bimaterial_wave_operator.{hpp,inl}`, `drivers/spatial_dyn_driver.cpp`,
  `general/annotation.hpp`, `config/defaults.mk`, `Makefile`,
  `dynamic/*_substep_iterator.*`.
- Domain context: `miniapps/seas/CLAUDE.md` (extreme-care files; no-local-mesh
  rule), memories `makefile-no-header-deps-stale-o`,
  `worktree-build-mfem-dir-override`.
- Verified-correct plan claims (NOT findings): (i) `mfem.hpp:29` →
  `annotation.hpp`, and all three target files include `mfem.hpp`, so the macros
  are in scope; (ii) `BimaterialWaveOperator` overrides only
  `ComputeMaxDt/SetMixedFluxMode/UsePrecomputedFaceFluxes/FluxForElem_/Interior
  FaceFlux_/SharedInteriorFaceFlux_/ApplyElementJacobian_` — so annotating the
  base `Mult/AdvanceADER/ComputeADERSubStepStates/Compute*FaceFluxRHS` DOES cover
  the matrix path; (iii) the seas Makefile links via `$(MFEM_LIBS)`, which
  carries `CALIPER_LIB` once MFEM is Caliper-built.

## Findings

### [R-001] [MODERATE] [caliper_perfgraph_plan.md §Phase 3 build recipe] — Caliper rebuild omits the seas-object `make clean`; perfgraph silently comes up empty

**Category:** BUG

**Description:**
Phase 3 step 3 is `cd miniapps/seas && make seas_spatial_dyn_driver` after
flipping `MFEM_USE_CALIPER=YES`. The plan itself documents (in the Local
Unit-Test Verification section) that the seas Makefile has **no header
dependencies**, so editing/headers-flipping does not force recompilation of
existing `.o`. The recipe fails to apply that same rule to itself: on Frontera a
non-Caliper `seas_spatial_dyn_driver` almost always already exists, so an
incremental `make` sees the unchanged `.cpp` files and **reuses the stale `.o`
that compiled every `MFEM_PERF_SCOPE` as a no-op**. The driver then links
`libcaliper` but emits **zero seas regions** — an empty/library-only perfgraph,
with no error to signal the problem.

**Trigger:**
A pre-existing (non-Caliper) seas build, then: reconfigure MFEM with Caliper →
`make seas_spatial_dyn_driver` (incremental, no clean).

**Actual behavior:**
Stale no-op seas objects are relinked; `CALI_CONFIG=runtime-report` shows only
MFEM-library regions (if any), none of the `seas::*` regions from Phases 1–2.

**Expected behavior:**
The seas translation units recompile against the Caliper-enabled `config.hpp`
so the `MFEM_PERF_SCOPE` regions are live.

**Suggested fix:** (edit the plan recipe)
```diff
- 3. **Rebuild the seas driver:** `cd miniapps/seas && make seas_spatial_dyn_driver`.
+ 3. **Rebuild the seas driver (force-clean — same no-header-deps hazard as the
+    local verification):** `cd miniapps/seas && make clean && make seas_spatial_dyn_driver`.
+    Skipping `make clean` here is the #1 cause of an EMPTY perfgraph: stale
+    objects compiled before the Caliper flip expanded every MFEM_PERF_SCOPE to
+    a no-op.
```

**Test case (verification check, not a Python unit test — this is a build/runtime concern):**
```bash
# After the Caliper build + a short run with CALI_CONFIG=runtime-report:
#   PASS iff at least one seas region appears.
grep -q 'seas::WaveOperator::Mult' runtime_report.txt \
  || { echo "FAIL R-001: no seas regions — stale no-op objects (missing make clean)"; exit 1; }
```

---

### [R-002] [MODERATE] [caliper_perfgraph_plan.md §Phase 1/Phase 2] — Friction solve (`iterator.Advance`) is not annotated; the perfgraph cannot separate friction cost from wave-propagation cost

**Category:** DEVIATION (from option-A intent "annotate the dynamic-rupture hot loops")

**Description:**
The coupled rate-and-state / LSW friction solve runs inside
`iterator.Advance(dof_data, fault_coords, ...)` at `spatial_dyn_driver.cpp:449`
(virtual dispatch into `tpv205/tpv104/tpv102_substep_iterator`). This is a
primary hot path of a dynamic-rupture solve (a per-QP nonlinear solve every
ADER sub-step). The plan annotates the wave-operator methods, the macro-step
wrapper (`AdvanceADERWithSubStep`), and the step loop — but **not** the friction
solve. Its time is therefore folded invisibly into
`seas::spatial_dyn::AdvanceADERWithSubStep` together with
`ComputeADERSubStepStates` + `AdvanceADER`. The single most common question a
perfgraph is built to answer for a coupled solver — *"how much time is in
friction vs. wave propagation?"* — is unanswerable as planned.

**Trigger:**
Any perfgraph run; the friction cost has no dedicated node.

**Actual behavior:**
No `seas::*friction*` region; friction time is hidden inside the macro-step
wrapper.

**Expected behavior:**
A distinct region around the friction sub-step solve.

**Suggested fix:** (add to plan Phase 2; the implementer then adds one region)
Primary (covers all iterator subclasses with a single edit at the call site):
```diff
  // drivers/spatial_dyn_driver.cpp, around line 449
+ MFEM_PERF_BEGIN("seas::spatial_dyn::friction_substep");
  iterator.Advance(dof_data, fault_coords,
                   ...);
+ MFEM_PERF_END("seas::spatial_dyn::friction_substep");
```
Note for the fix agent: `MFEM_PERF_BEGIN/END` are not RAII — if a future change
lets `iterator.Advance` throw between them the region would be left open. The
exception-safe alternative is `MFEM_PERF_SCOPE("seas::FrictionIterator::Advance")`
at the top of each concrete `Advance` body in
`dynamic/{tpv205,tpv104,tpv102}_substep_iterator.cpp`. Either is acceptable;
the call-site form is one edit and `Advance` failing is currently fatal anyway.

**Test case (verification check):**
```bash
grep -q 'friction_substep' runtime_report.txt \
  || { echo "FAIL R-002: friction solve not broken out in the perfgraph"; exit 1; }
```

---

### [R-003] [LOW] [caliper_perfgraph_plan.md §Local Unit-Test Verification] — the test-run loop masks a failing exit code

**Category:** QUALITY

**Description:**
The verification loop
`for T in ...; do make ... "$T" && ./"$T" || { echo "FAIL: $T"; break; } done`
prints `FAIL` and `break`s on a failed test, but the loop — and the surrounding
shell — still terminate with exit status 0. An automated/CI harness (or an
inattentive operator) would read the run as a PASS despite a failed unit test,
defeating the purpose of the gate.

**Trigger:**
Any one of the 6 unit tests fails (e.g., a real regression, or the stale-`.o`
SIGABRT the plan warns about).

**Actual behavior:**
Loop exits 0; "verification passed" is inferred incorrectly.

**Expected behavior:**
Non-zero exit on any test failure.

**Suggested fix:**
```diff
+ rc=0
  for T in \
      seas_test_bimaterial_wave_operator_parity \
      ... ; do
-   make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN "$T" \
-     && ./"$T" || { echo "FAIL: $T"; break; }
+   make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN "$T" \
+     && ./"$T" || { echo "FAIL: $T"; rc=1; break; }
  done
+ [ "$rc" -eq 0 ] || { echo "VERIFICATION FAILED"; exit 1; }
```

**Test case:** n/a (shell-recipe hygiene; the fix is self-demonstrating —
inject a failing test and confirm non-zero exit).

---

### [R-004] [LOW] [caliper_perfgraph_plan.md §Phase 2 / Testing Strategy] — Phase 2 driver annotations have NO behavioral test coverage; only "links + diff" guards them

**Category:** ASSUMPTION

**Description:**
The annotations in Phase 2 live in `drivers/spatial_dyn_driver.cpp`
(`AdvanceADERWithSubStep_Spatial` + the time loop), which is driver-local code
not exercised by any unit test — the 6 listed tests cover the WaveOperator
library only. The plan acknowledges "the driver isn't behaviorally tested," but
buries it; the *only* real guard that a Phase-2 edit is harmless is the
"annotations-only `git diff`" acceptance criterion. That criterion should be
called out as the mandatory Phase-2 gate, not a footnote, since a misplaced
macro (e.g., inserted outside the loop body, or swallowing a statement) would
still "link" and still "pass" the wave-operator tests.

**Trigger:**
A Phase-2 macro inserted at the wrong scope.

**Actual behavior:**
Plan implies the unit-test set covers Phase 2; it does not.

**Expected behavior:**
Plan states explicitly that Phase 2's sole local gate is the annotations-only
diff + clean compile/link, and instructs the reviewer to diff
`spatial_dyn_driver.cpp` line-by-line.

**Suggested fix:**
```diff
  ### Acceptance Criteria   (Phase 2)
+ - [ ] **Sole behavioral gate:** there is NO unit test over
+   `spatial_dyn_driver.cpp`; Phase 2 is verified ONLY by (a) clean compile +
+   link and (b) a line-by-line `git diff` confirming exactly two inserted
+   `MFEM_PERF_*` lines and no other change. The reviewer MUST inspect the diff.
```

**Test case:** n/a (process gate).

---

### [R-005] [LOW] [POSSIBLE] [caliper_perfgraph_plan.md §Phase 3 sbatch reconciliation] — merge target for reuniting the annotations with the `TPV31_PERFGRAPH` sbatch hook is unspecified

**Category:** ASSUMPTION

**Description:**
The `TPV31_PERFGRAPH=1`→`CALI_CONFIG` block lives in the
`system/spatial_dyn_driver` working tree (uncommitted), NOT in this worktree
branch. The plan says "reconcile at merge / ensure exactly one copy" but never
states the merge target. If the worktree branch were merged somewhere other than
`system/spatial_dyn_driver` (e.g., a fresh branch off origin/main), the
annotations and the runnable sbatch hook would live on different branches and
the feature would be un-demonstrable end-to-end. This is POSSIBLE rather than
definite because the intended merge target is probably `system/spatial_dyn_driver`
— but the plan should pin it.

**Suggested fix:**
```diff
  ### Files to Modify   (Phase 3 reconciliation)
+ - **Merge target:** this worktree branch (`worktree-feature-caliper-perfgraph`)
+   MUST be merged into `system/spatial_dyn_driver`, where the `TPV31_PERFGRAPH`
+   sbatch hook and the p1/p2 TPV31 jobs already live (currently uncommitted on
+   that branch's working tree — commit them first). Merging elsewhere separates
+   the annotations from their runner.
```

**Test case:** n/a (merge-process note).

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 empty-perfgraph build-recipe gap; R-002 friction
  solve not annotated)
- Low issues: 3 (R-003 masked exit code; R-004 Phase-2 coverage; R-005 merge
  target)
- Plan compliance: **PARTIAL** — the plan faithfully covers the wave-operator +
  driver-loop scope it was given and its codebase claims check out, but (a) the
  build recipe would yield an empty perfgraph on the common incremental-build
  path, and (b) it omits the friction solve, which undercuts the perfgraph's
  main purpose for a coupled solver.
- Verdict: **PASS WITH FIXES** — apply R-001 and R-002 before implementing;
  R-003–R-005 are cheap and should go in the same pass. No finding blocks the
  annotation approach itself (which is sound and zero-impact on the default
  build).

## Unreviewed Areas
- The actual Caliper-enabled build + perfgraph emission was NOT exercised (no
  Caliper in this checkout; no local full-mesh runs per project rule). R-001 is
  reasoned from the Makefile's no-header-deps behavior, not observed.
- `CALI_CONFIG` spec strings (`runtime-report(calc.inclusive=true)`,
  `hatchet-region-profile(output=...)`) were checked for plausibility against
  Caliper's ConfigManager conventions but not run; the implementer should
  smoke-test them once on a Caliper build.
- The concrete `Advance` bodies in `tpv205/tpv104/tpv102_substep_iterator.cpp`
  were located but not line-audited (only relevant to R-002's alternative fix).
