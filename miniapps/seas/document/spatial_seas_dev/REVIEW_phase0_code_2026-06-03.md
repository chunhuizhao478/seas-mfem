# Code Review: Phase 0 implementation — `spatial_seas_driver.cpp` + Makefile (2026-06-03)

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md`, **Phase 0** ("Build target, driver skeleton, dry-run").
- Files reviewed (this changeset only):
  - `drivers/spatial_seas_driver.cpp` — **new** (Phase 0 QD skeleton).
  - `Makefile` — **modified** (added `SPATIAL_SEAS_DRIVER_SRC/OBJ` vars, the object-compile rule, and the `seas_spatial_seas_driver` link target).
  - `extern/toml11` — symlinked to the main checkout (gitignored; build prerequisite per project memory).
- Domain context consulted: `CLAUDE.md` (both root + worktree), project memory (worktree build invocation + toml11 symlink; no-touch constraint on shared `bp5/domain/fault/solver` + `spatial_friction` parser semantics), reference driver `drivers/spatial_dyn_driver.cpp` (the CLI-helper + config-merge structure Phase 0 says to mirror).
- Reviewer stance: adversarial — three passes (plan-compliance, bug-hunt, quality). The build was exercised (clean, no warnings) and all documented CLI paths were run.

## Build & smoke evidence (worktree)
Built with the memory invocation: `make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN seas_spatial_seas_driver` → **clean, no errors/warnings**.

| Smoke | Command | Result |
|---|---|---|
| missing `--config` | `./seas_spatial_seas_driver` | `ERROR: --config ... required`, **exit 2** ✔ |
| good config dry-run | `--config <rs_safs>.toml --dry-run` | full banner, **exit 0** ✔ |
| CLI overrides | `--mesh --tfinal --output-dir --paraview-fault vtu --paraview-volume off --checkpoint-every 7` | every field overridden in banner ✔ |
| non-dry-run | `--config <rs>.toml` | "Phase 0 skeleton … not yet implemented", **exit 0** ✔ |
| MPI np=2 dry-run | `mpirun -np 2 … --dry-run` | exactly **one** rank-0 banner, exit 0 ✔ |
| missing file | `--config /no/such.toml` | parser `MFEM_ABORT`, **exit 1** (see R-001) |
| malformed TOML | `--config <bad>.toml` | parser `MFEM_ABORT`, **exit 1** (see R-001) |

Config used for smokes: `safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml` (a rate-state spatial-schema config — the BP5/`config/safs_qd/` QD configs are not authored until Phase 7, so an existing parseable RS config was used; the dry-run never loads the mesh, so this is safe locally).

## Findings

### [R-001] LOW [spatial_seas_driver.cpp:main] — graceful config-load failure path is unreachable (parser aborts instead of throwing)

**Category:** ASSUMPTION / DEVIATION (matches reference driver)

**Description:**
The driver wraps `LoadSpatialFrictionConfig` in `try { … } catch (...) { print "ERROR: failed to load TOML config" ; return 2; }`. But `LoadSpatialFrictionConfig` / `ParseSpatialFrictionConfigString` (`spatial_friction.cpp:1688`, `:1711`) report errors with `MFEM_ABORT` (missing file at `:1698`, toml parse error at `:1722`), **not** C++ exceptions. So neither a missing file nor malformed TOML is caught: the process aborts (exit 1) with a raw MFEM stack instead of the intended `return 2`. The `catch(...)` is effectively dead code.

**Trigger:** `--config /no/such/file.toml` or any syntactically-invalid TOML.

**Actual behavior:** raw `MFEM_ABORT` banner, exit 1.

**Expected behavior (per the catch block's intent):** `ERROR: failed to load TOML config '<path>'.` on rank 0, then exit 2.

**Assessment:** This is **identical to the established `spatial_dyn_driver.cpp:595-610` pattern** — the reference driver has the same dead `try/catch`. So it is consistent with the codebase and not a Phase-0 regression. It is flagged LOW because (a) the only reachable clean error return is the missing-`--config` case (which works, exit 2), and (b) the catch is harmless defensive code should the shared parser ever migrate to exceptions. **Do not** "fix" it by editing the shared parser to throw (no-touch constraint on shared `spatial_friction`).

**Suggested fix (optional, QD-driver-local only):** pre-check existence before the parse so the common case gets a clean message + exit 2, and keep the catch for future-proofing:
```diff
+   {
+      std::ifstream probe(config_path);
+      if (!probe.good())
+      {
+         if (rank == 0) { std::cerr << "ERROR: --config '" << config_path
+                                    << "' is not readable.\n"; }
+#ifdef MFEM_USE_MPI
+         MPI_Finalize();
+#endif
+         return 2;
+      }
+   }
    spatial::SpatialFrictionConfig cfg;
    try { cfg = spatial::LoadSpatialFrictionConfig(config_path); }
```
(Malformed-TOML would still abort via the shared parser — acceptable and matching the reference.)

**Test case:**
```cpp
// test_R001_qd_missing_config_file_clean_exit:
//   run `seas_spatial_seas_driver --config /no/such.toml --dry-run`
//   EXPECT exit code 2 (not 1) and stderr contains "not readable".
// (Currently would FAIL — documents the gap; only enable if the pre-check
//  fix is applied.)
```

---

### [R-002] LOW [spatial_seas_driver.cpp:banner] — banner echoes shared default `paraview_fault = "hdf5"`, contradicting the QD "VTU/PVD, never vtkhdf" directive

**Category:** DEVIATION (plan consistency) — out of Phase-0 enforcement scope

**Description:**
The banner prints `cfg.output.paraview_fault`, whose shared-struct default (`OutputSpec`, `spatial_friction.hpp:169`) is `"hdf5"`. The QD plan's I/O scope change (2026-06-03) mandates **VTU/PVD output, never vtkhdf** for the QD driver. At Phase 0 the driver only *echoes* the config (the ParaView wiring + the vtu enforcement is Phase 6), so this is not yet a behavioral bug — but an operator reading the Phase 0 banner sees `paraview fault: hdf5` and may believe the QD driver will write `.vtkhdf`.

**Trigger:** any dry-run with a config that omits `[output].paraview_fault` (i.e. inherits the shared default).

**Actual behavior:** `paraview fault:   hdf5`.

**Expected behavior:** the QD driver should ultimately force `vtu`. At Phase 0, at minimum the banner should not imply hdf5 is the QD default.

**Suggested fix (defer to Phase 6, or a one-line Phase-0 note):**
```diff
+      // QD writes VTU/PVD (plan I/O scope 2026-06-03); the shared OutputSpec
+      // default is "hdf5" — Phase 6 forces vtu.  Flag it so the banner is honest.
+      const std::string pv_fault_eff =
+         (cfg.output.paraview_fault == "hdf5") ? "hdf5 (Phase 6 will force vtu)"
+                                               : cfg.output.paraview_fault;
-                << "paraview fault:   " << cfg.output.paraview_fault << "\n"
+                << "paraview fault:   " << pv_fault_eff << "\n"
```

**Test case:** N/A at Phase 0 (echo-only); the binding test belongs in Phase 6 ("QD fault output mode resolves to Vtu regardless of TOML `paraview_fault`").

---

### [R-003] LOW [spatial_seas_driver.cpp:main] — non-dry-run invocation returns 0 while doing nothing

**Category:** QUALITY / EDGE_CASE

**Description:**
A non-dry-run run (no `--dry-run`) prints the "Phase 0 skeleton: not yet implemented" notice and returns **0**. A CI/automation wrapper that only inspects the exit code would read this as a successful simulation. This is an intentional skeleton choice (so any future smoke target that invokes the binary without `--dry-run` does not fail the suite), and the message is explicit, so it is LOW. No aggregate (`all`/`test`) references the target today, so nothing depends on either choice yet.

**Suggested fix (optional):** keep exit 0 but route the notice to `stderr`, or return a distinct sentinel (e.g. 64 `EX_USAGE`) until the time loop lands. Recommend leaving as-is and revisiting when Phase 2+ fills the body.

**Test case:** N/A (design choice, not a correctness bug).

---

### [R-004] LOW [spatial_seas_driver.cpp:banner] — `dt_initial` sentinel display — FIXED this round

**Category:** QUALITY

**Description:** the first cut printed `dt_initial:       -1  (auto)` (raw sentinel + label). Polished to print `dt_initial:       auto` when `dt_initial < 0`, else the value with units. Rebuilt + re-verified (`dt_initial:       auto`). Recorded for completeness; no action needed.

---

## Summary
- Critical issues: 0
- Moderate issues: 0
- Low issues: 4 (R-001 dead catch [matches reference]; R-002 banner hdf5 default [Phase-6 scope]; R-003 non-dry-run exit 0 [by design]; R-004 banner polish [fixed])
- Plan compliance (Phase 0): **FULL** —
  - ✔ New `drivers/spatial_seas_driver.cpp` skeleton: MPI init, CLI parse (`--config`, `--dry-run`, `--mesh`, `--tfinal`, `--restart`, `--checkpoint-every`, QD `--paraview-*` subset), TOML load via `LoadSpatialFrictionConfig`, rank-0 banner, `MPI_Finalize`. Dynamic-only flags (`--ader-order/--mixed-flux/--pml*/--time-integrator`) correctly **omitted**.
  - ✔ CLI sentinels match the spatial convention (empty string / `<0` ⇒ keep TOML; CLI wins).
  - ✔ `--dry-run` exits after the config echo + derived-parameter print, **before** any operator/mesh construction.
  - ✔ Build-from-worktree invocation documented in the header comment.
  - ✔ `Makefile`: `SPATIAL_SEAS_DRIVER_SRC/OBJ` mirror the dynamic block; object rule uses `$(TOML_FLAGS) -DSEAS_USE_MPI`; link set is the **minimal non-wave** set `$(SPATIAL_FRICTION_OBJ)` (the only out-of-TU symbol referenced; proven by `seas_test_spatial_friction_config`'s identical link set). **No wave / Godunov / PML / ADER objects.**
  - ✔ Acceptance #1 (`make seas_spatial_seas_driver` succeeds) and #2 (`--dry-run` prints config, exit 0) **met**.
- Verdict: **PASS** — Phase 0 is implementable and implemented; the 4 findings are LOW (one already fixed; two are deliberately deferred to their owning phase; one matches the reference driver).

## Non-interference with existing drivers/tests (Acceptance #3)
Acceptance criterion #3 is "`make all && make test` still pass (no regression)". The full `make all && make test` was **not** run end-to-end (it builds the entire wave stack + every driver and is heavy; see the pre-existing failure below). Instead, non-interference was verified structurally + by spot build:
- The new target `seas_spatial_seas_driver` is **not** referenced by any `all` / `test` aggregate (grep-confirmed: only the 2 new var defs + the obj-rule + the link target mention it). Existing aggregates are byte-unchanged.
- `make -n` parses cleanly for `seas_spatial_dyn_driver`, `seas_driver`, and `seas_test_spatial_friction_config` after the Makefile edit (no rule corruption).
- My edits add only new symbols; no existing variable, rule, or source file was modified (`git status`: only `M Makefile`, `?? drivers/spatial_seas_driver.cpp`).
- The shared `spatial/code/spatial_friction.{hpp,cpp}` are **untouched** (git-confirmed) — the QD driver links the same `spatial_friction.o` the dynamic driver and the config test already use.

## Pre-existing failure (NOT introduced by this changeset — do not fix here)
- `seas_test_spatial_friction_config` **aborts (exit 134)** in sub-case `T_23_numerics_selectors_parse` ("[NUM-1]"): it calls `ParseSpatialFrictionConfigString` **directly in the main process** (not via the `ParseAbortsInChild` fork) with a `[material] kind="depth_profile_1d"` fixture (`test_spatial_friction_config.cpp:1074-1076`) that provides only `profile_csv` and **no** `[[material_profile.layer]]` array, which the parser now requires (`spatial_friction.cpp:1622`, `MFEM_VERIFY(root.contains("material_profile"))`). This is a **parser↔test schema drift on this branch**, entirely within shared spatial-friction code. It is unrelated to Phase 0 (which neither touches that file nor that code path), and it sits in the no-touch shared parser — so it must NOT be fixed under Phase 0. Recommend a separate ticket: update the NUM-1 fixture to supply a `[[material_profile.layer]]` table (or wrap it in `ParseAbortsInChild` if the abort is now intended).

## Unreviewed Areas
- Full `make all && make test` not executed (heavy; non-interference verified structurally as above; one pre-existing unrelated failure documented).
- The QD `config/safs_qd/` configs do not exist yet (Phase 7) — Phase 0 was smoked against an existing RS spatial config.
- Phases 1-9 behavior (operators, RK45 loop, ParaView VTU wiring, checkpoint) — out of Phase 0 scope; the skeleton's non-dry-run path is a placeholder by design.
