# Code Review: Phase 4 — D-rename (ComputeSAFSParams→ComputeParams) — 2026-05-27

> NOTE: The pre-existing `REVIEW.md` in this directory is an **unrelated** TPV205 review
> (round 3, 2026-04-27, incl. CRITICAL R-016) and has been **left intact**. This Phase-4
> review lives in its own file. Point `/code-fix` at THIS file.

## Review Scope
- **Plan:** `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`, §Phase 4 (lines 1387–1413)
- **Worktree under review:** `/Users/chunhuizhao/projects/seas-mfem-safs` (branch `safs` @ `9f7f683`) — NOTE: the rename was applied here, NOT in the primary working dir `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver`. The `/code-fix` agent must operate on the `seas-mfem-safs` worktree.
- **Files reviewed (18):** `config/seas_config.hpp`, `config/seas_config_bridge.hpp`, `drivers/spatial_dyn_driver.cpp`, `fault/fault_geometry.hpp`, `fault/fault_geometry_safs.inl`, `fault/fault_geometry_safs_templated.inl`, `fault/rate_state_fault.hpp`, `spatial/code/spatial_stress.cpp`, `spatial/code/spatial_stress.hpp`, `tests/unit/test_compute_safs_params.cpp`, `tests/unit/test_spatial_stress_bundle.cpp`, `Makefile`, 5 `safs/.../config/spatial_friction_slip_weakening_safs_projected_stress*.toml`, `spatial/code/scripts/verify_constant_tensor_projection.py`
- **Domain context:** root `CLAUDE.md`, `miniapps/seas/CLAUDE.md` ("Files Requiring Extreme Care" table; "don't revert a fix"; build/test conventions), implementation completion report (conversation).

## Independent verification performed (not merely trusting the report)
- Tree-wide grep with **no** extension filter: only **16 `.md`** files retain old symbols; **zero** non-`.md` files (rules out `.h/.cc/.cxx/.txx/.ipp/.tpp` that the implementer's `--include='*.hpp,*.cpp,*.inl'` enumeration could have missed).
- `class FaultGeometry` (`fault_geometry.hpp:59`) has **no base class** → `ComputeParams`/`HasParams` cannot hide/override a base virtual.
- Both `ComputeParams` overloads (non-template `StressField3D` + templated `StressSource`) and member `params_computed_` intact with default args (`fault_geometry.hpp:580–642`).
- Zero-collateral proof: every one of the 18 files == `git baseline blob + exact token substitution` (no incidental edits).
- Compile verified across **all** renamed headers via 6 targets, incl. the extreme-care `rate_state_fault.hpp` via `seas_test_safs_mode_wiring` (EXIT=0).
- `seas_test_compute_safs_params`: 13/13 pass post-rename incl. T-65-5 (templated == sidecar); T-65-4 numerics byte-identical. Driver binary symbol audit: 5 `ComputeParams`, **0** `SAFSParams`.

**No correctness/logic bugs were found in the rename itself.** All findings below are plan-compliance and verification-completeness items. (Adversarial passes 1–3 were executed; the rename is mechanically clean, so manufacturing additional "bugs" would violate review rules 4/6.)

## Findings

### [R-001] MODERATE [acceptance][drivers/spatial_dyn_driver.cpp + tree] — SAFS+RS smoke (plan-required acceptance) was not executed

**Category:** DEVIATION

**Description:**
Plan §Phase 4 "Acceptance Criteria" lists TWO gates: (a) `seas_test_compute_safs_params` (incl. T-65-5) passes, and (b) "the SAFS+RS smoke ... pass unchanged after the rename." Gate (a) was executed and passes. Gate (b) was **not run** (deferred as "heavyweight"). The rename is strongly evidenced behavior-neutral (byte-identical T-65-4 numerics, symbol-only diff, clean compile), but the plan-mandated integration check is unverified, and `rate_state_fault.hpp` is on the `miniapps/seas/CLAUDE.md` "Files Requiring Extreme Care" list whose policy is "Any change to these files requires running full verification tests."

**Trigger:** Running the renamed `seas_spatial_dyn_driver` end-to-end with an RS config; or `make test`.

**Actual behavior:** Acceptance gate (b) not exercised; "passes unchanged" asserted by argument, not by run.

**Expected behavior:** SAFS+RS smoke (and ideally `make test`) run and confirmed identical to the pre-rename baseline.

**Suggested fix:** Run the smoke + unit suite in the `safs` worktree; record results. No source change.
```diff
# (verification, not a code edit) — in seas-mfem-safs/miniapps/seas, `conda activate mfem-dev`:
+ make test                      # full unit suite (incl. rate_state_fault.hpp consumers)
+ mpirun -np 8 ./seas_spatial_dyn_driver \
+   --config safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml \
+   --print-derived            # confirm clean exit + aging-iterator dispatch
```

**Test case:**
```python
def test_R001_safs_rs_smoke_unchanged():
    # baseline: probe/checkpoint output from a pre-rename short SAFS-RS run (small tfinal)
    # action:   run the renamed driver with identical config + ranks + seed
    # assert:   per-DOF psi, slip1/2, V1/2, tau*_corr, sigma_n at final step are
    #           bit-identical (or within RK tol) to the baseline
    assert smoke_outputs_equal(baseline_run, renamed_run)
```

---

### [R-002] LOW [acceptance][*.md docs] — Literal acceptance grep does NOT "return nothing"

**Category:** DEVIATION

**Description:**
Plan §Phase 4 acceptance: `grep -rn 'ComputeSAFSParams\|HasSAFSParams' miniapps/seas` "returns nothing." Actual: **85 matches across 16 `.md` docs** (debug history, reviews, and the PLAN doc itself, which documents the old→new rename). This is a **user-approved** scope decision (rename code + active comments only; leave docs), so the deviation is justified — renaming the PLAN doc would make it self-contradictory and renaming debug records rewrites history. The plan's literal acceptance text is unmet and should be reconciled so it is not falsely reported as failing.

**Trigger:** `grep -rn 'ComputeSAFSParams\|HasSAFSParams' miniapps/seas`

**Actual behavior:** 85 hits in 16 `.md` files.

**Expected behavior (reconciled):** Acceptance scoped to code/active source returns nothing (verified: 0 in `*.hpp/*.cpp/*.inl/Makefile/*.toml/*.py`).

**Suggested fix:** Update the plan's acceptance wording to scope the grep to code/active source (preferred — preserves the docs), rather than renaming the 16 docs.
```diff
- - [ ] `grep -rn 'ComputeSAFSParams\|HasSAFSParams' miniapps/seas` returns nothing.
+ - [ ] `grep -rn 'ComputeSAFSParams\|HasSAFSParams' miniapps/seas --include='*.hpp' --include='*.cpp' --include='*.inl' --include='*.toml' --include='*.py' --include='Makefile'` returns nothing.
+ #   (Per user scope decision 2026-05-27: 16 historical .md docs — incl. this plan —
+ #    intentionally retain the old names; renaming the plan would make it self-referentially wrong.)
```

**Test case:**
```python
def test_R002_grep_zero_in_code_scope():
    import subprocess
    out = subprocess.run(
        ["grep","-rIn","ComputeSAFSParams\\|HasSAFSParams\\|safs_params_computed_",
         "--include=*.hpp","--include=*.cpp","--include=*.inl","--include=*.toml",
         "--include=*.py","--include=Makefile","."],
        cwd="miniapps/seas", capture_output=True, text=True)
    assert out.stdout.strip() == ""
```

---

### [R-003] LOW [process][seas-mfem-safs worktree] — Change is uncommitted, isolated to one worktree, and git is broken inside it

**Category:** QUALITY

**Description:**
The rename exists only as **uncommitted edits** in the `seas-mfem-safs` worktree (branch `safs`). The primary worktree `seas-mfem-spatial-dyn-driver` (branch `system/spatial_dyn_driver`, same commit `9f7f683`) still has the OLD names. Additionally, the `safs` worktree's `.git` file points to a stale gitdir (`/Users/chunhuizhao/projects/seas-mfem/.git/worktrees/seas-mfem-safs`, nonexistent), so `git status/diff/commit` **fail from inside the safs worktree** — the change cannot be committed there until repaired. Risk: rename lost on `make clean`/checkout, or the two branches silently diverge.

**Trigger:** `git -C /Users/chunhuizhao/projects/seas-mfem-safs status` → `fatal: not a git repository: .../seas-mfem/.git/worktrees/seas-mfem-safs`.

**Actual behavior:** Rename un-committable from the worktree where it lives; absent on the primary branch.

**Expected behavior:** Worktree git pointer repaired; rename committed to the intended branch.

**Suggested fix:** Non-destructive worktree repair from the real repo, then commit on the intended branch.
```diff
+ # the real repo is seas-mfem-spatial-dyn-driver/.git; its registry already points correctly
+ git -C /Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver worktree repair \
+     /Users/chunhuizhao/projects/seas-mfem-safs
+ # then, in the safs worktree, review `git diff` and commit the 18-file rename.
```

**Test case:** (process check)
```python
def test_R003_git_usable_in_safs_worktree():
    import subprocess
    r = subprocess.run(["git","status","--porcelain"],
                        cwd="/Users/chunhuizhao/projects/seas-mfem-safs",
                        capture_output=True, text=True)
    assert r.returncode == 0   # git works in the worktree (post-repair)
```

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — plan-required SAFS+RS smoke not executed)
- Low issues: 2 (R-002 literal acceptance-grep deviation [user-approved]; R-003 uncommitted/worktree git breakage)
- **Plan compliance:** PARTIAL — the *code* rename is FULL and correct (all code + active comments; both `ComputeParams` overloads + `params_computed_` intact; guardrail req 2 honored — no R-001 sign flip, no `kind`-string or projection-rule changes). The *acceptance* is partial: smoke deferred (R-001) and literal grep scoped-to-code by approved decision (R-002).
- **Verdict:** PASS WITH FIXES — no source changes required for correctness. To fully satisfy the plan: (1) run the SAFS+RS smoke + `make test` [R-001], (2) reconcile the acceptance-grep wording [R-002], (3) repair the worktree git pointer and commit [R-003].

## Unreviewed Areas
- **SAFS+RS smoke runtime behavior** — not executed (R-001); behavior-neutrality argued from byte-identical T-65-4 numerics, symbol-only diff, and clean compilation, but not run-verified.
- **Full `make test` suite** — not run; compile coverage of the renamed headers verified via 6 representative targets (incl. `seas_test_safs_mode_wiring` for extreme-care `rate_state_fault.hpp`). `seas_test_spatial_stress_bundle`'s runtime is a **pre-existing** BP5-fixture segfault (excluded from `make test`), unrelated to this rename.
- **16 `.md` docs** — intentionally not renamed per user scope decision; out of correctness scope.
