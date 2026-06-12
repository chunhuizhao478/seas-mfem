# Code Review (Round 4 — post-fix, post-Phase-4): Free-Surface Slice — 2026-06-02

> Rounds 1–2 reviewed the **plan**; Round 3 reviewed the **written Phases 1–3** and flagged
> the entirely-missing Phase 4 (R-301) plus that the trace-equivalence and collective paths
> were *unverified by execution*. This round audits the **fix + Phase-4 implementation** and,
> crucially, the changes were **built and run** (np=1 and np=2) — which surfaced a CRITICAL
> latent bug that no amount of static review had caught. Finding IDs continue at R-401+.
> Rounds 1–3 are superseded (git history retains them).

## Review Scope
- Files reviewed (final state, all built):
  - `io/free_surface_output.hpp` (Phase 1 + the R-401 fix)
  - `drivers/spatial_dyn_driver.cpp` (R-302/R-303 fixes from Round 3)
  - `tests/unit/test_free_surface_slice.cpp` (new, Phase 4 — MPI class test)
  - `tests/unit/test_spatial_friction_config.cpp` (+5 Phase-4 parse tests)
  - `Makefile` (target wiring), `verify_spatial_dyn_smoke_safs.py` (Check F)
- Verification performed (worktree build: `MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN
  MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN`, `conda activate mfem-dev`):
  - `seas_test_free_surface_slice`: **PASS at np=1 AND np=2** (Dimension==2, GlobalNE==4,
    trace ≤ 1e-12, schedule, file writes).
  - `seas_test_spatial_friction_config` (5 new FS tests): **PASS** (defaults, off round-trip,
    illegal-value/dt≤0/negative-attr aborts).
  - `seas_spatial_dyn_driver`: **compiles** with all fixes.
  - `verify_spatial_dyn_smoke_safs.py`: `py_compile` OK.

## Findings

### [R-401] [CRITICAL — FOUND & FIXED] `io/free_surface_output.hpp` — sub-space `L2_FECollection` built with the parent dimension (3) instead of the submesh dimension (2) → SIGSEGV in every run

**Category:** BUG (plan defect, inherited by the implementation)

**Description:**
Phase 1 §3 step 4 (and the original code) built the **submesh** velocity/rank
`L2_FECollection` with `dim = 3` (the parent's spatial dimension). The submesh is codim-1
(2D). MFEM's boundary→submesh L2 transfer (`SubMeshUtils::BuildVdofToVdofMap`, reached from
the `ParTransferMap` ctor) requires the sub collection's dimension to be the **submesh**
dimension; a `dim=3` collection produces inconsistent trace bookkeeping and **dereferences a
null pointer**. Because the slice is **default-ON** and the driver always builds the velocity
`ParTransferMap` at construction, this would have **SIGSEGV'd at startup on every TPV/SAFS
run** — including the committed Frontera production jobs. Static review (Rounds 1–3) could not
catch it; only building + running did. MFEM's own boundary-transfer unit test
(`tests/unit/mesh/test_submesh.cpp` via `create_fec(..., submesh->Dimension())`) confirms the
correct dimension.

**Trigger:** constructing `FreeSurfaceOutput` on any 3D parent (i.e. every run).

**Actual behavior (pre-fix):** `Signal 11 (SIGSEGV)` in
`mfem::SubMeshUtils::BuildVdofToVdofMap` ← `ParTransferMap::ParTransferMap` ←
`FreeSurfaceOutput` ctor, reproduced at np=1 and np=2.

**Expected behavior:** the transfer map builds; the sliced velocity equals the analytic trace.

**Fix applied:**
```diff
-      sub_vel_fec_ = std::make_unique<L2_FECollection>(order, 3, BasisType::GaussLobatto);
-      sub_rank_fec_ = std::make_unique<L2_FECollection>(0, 3, BasisType::GaussLobatto);
+      const int sub_dim = submesh_.Dimension();   // codim-1: parent_dim - 1
+      sub_vel_fec_  = std::make_unique<L2_FECollection>(order, sub_dim, BasisType::GaussLobatto);
+      sub_rank_fec_ = std::make_unique<L2_FECollection>(0,     sub_dim, BasisType::GaussLobatto);
```
(vdim stays 3 — a 3-component velocity living on the 2D surface. The **parent** spaces keep
`dim=3`, which is correct.)

**Test that now covers it:** `test_free_surface_slice.cpp` — constructs the writer on a
2×2×2 box and asserts `Dimension()==2`, `GlobalNE()==4`, and trace ≤ 1e-12. It SIGSEGV'd
before the fix and passes after (np=1 and np=2). This is the regression guard.

---

### [R-402] [LOW — HARDENED] `io/free_surface_output.hpp` — `FreeSurfaceOutput` was implicitly movable while holding an interior pointer

**Category:** EDGE_CASE / ROBUSTNESS

**Description:**
The fix changed `submesh_` from `unique_ptr<ParSubMesh>` to a **value member** (so the
submesh is built by guaranteed copy elision, avoiding any reliance on MFEM's incomplete
`ParSubMesh` move — a risk the plan itself flagged). Side effect: `pv_dc_` and the sub FE
spaces store interior pointers into `submesh_`, so the class must not relocate — yet it was
left implicitly movable. The driver holds it via `unique_ptr<FreeSurfaceOutput>` and never
moves the object, so no live bug, but the latent footgun was closed.

**Fix applied:** deleted copy and move:
```cpp
   FreeSurfaceOutput(const FreeSurfaceOutput &)            = delete;
   FreeSurfaceOutput &operator=(const FreeSurfaceOutput &) = delete;
   FreeSurfaceOutput(FreeSurfaceOutput &&)                 = delete;
   FreeSurfaceOutput &operator=(FreeSurfaceOutput &&)      = delete;
```
Confirmed the driver (`std::make_unique<FreeSurfaceOutput>(...)`, `fs_out.reset()`) still
builds — make_unique constructs in place and the `unique_ptr` move moves the pointer, not the
object.

---

### [R-403] [LOW] [POSSIBLE] `test_free_surface_slice.cpp` — the np=2 case may not exercise a *truly* empty-local-submesh rank

**Category:** ASSUMPTION (test coverage)

**Description:**
The production motivation for the `-np 2` test is to exercise ranks that own **zero** local
free-surface faces (interior ranks, the common case at scale). On a 2×2×2 box partitioned in
two, METIS most likely splits by a plane and gives **both** ranks some top quads — so the
genuinely empty-local path (local `submesh.GetNE()==0`, then collective `Transfer`/`Save`) may
not actually be hit. The code is correct by MFEM's design (`ParSubMesh::Dimension()` is a
topological property = parent_dim−1, valid even with 0 local elements, so `sub_dim==2` on
every rank), and np=2 *did* exercise the collective transfer across a partition — but the
zero-local-faces rank is not provably covered.

**Trigger:** an interior rank with no `kTop` faces.

**Suggested follow-up (not blocking):** add an np=3 or np=4 case, or a taller box
(`MakeCartesian3D(1,1,4)`) tagged only on the very top, so at least one rank is guaranteed
top-face-free; assert it still returns `GlobalNE()==expected` and trace ≤ 1e-12. Frontera at
production np is the ultimate coverage.

**Test case (sketch):**
```cpp
// MakeCartesian3D(1,1,4) tagged kTop only on z==z_max under -np 4 ⇒ ≥1 rank owns no top
// face; assert GlobalNE()==1 on all ranks and the global trace error ≤ 1e-12.
```

---

### [R-404] [MODERATE] [PRE-EXISTING, UNRELATED] `test_spatial_friction_config.cpp` aborts at the `depth_profile_1d` material test — blocks the full parse suite on TOML builds

**Category:** BUG (pre-existing; NOT caused by this change)

**Description:**
With `SEAS_USE_TOML` enabled, `seas_test_spatial_friction_config` aborts (SIGABRT) in
`T_45_material_kinds` at `[material].kind="depth_profile_1d"` + `profile_csv=...`: the test
expects a direct parse to **succeed**, but the parser (`spatial_friction.cpp:1640`) now
**requires** a `[[material_profile.layer]]` array and aborts. This is a test/parser mismatch
that predates this branch's free-surface work (my diff adds only output-block keys and 5
end-of-suite FS tests — neither touches material parsing). It was masked in the worktree
because the vendored `extern/toml11` header is absent there, so the suite was compiled with
TOML **disabled** and `main()` early-returned ("SEAS_USE_TOML not defined — skipping").

**Impact:** the 5 new FS parse tests are correct and **pass** (verified by temporarily moving
them ahead of `T_45`), but the *full* `make test-spatial-friction-config` cannot run green on
a TOML build until `T_45` (or the parser) is reconciled.

**Recommendation (separate from free-surface; needs owner decision per CLAUDE.md "don't fix
unrelated"):** either update `T_45` to supply a `[[material_profile.layer]]` block (if the
parser requirement is intended) or relax the parser (if `profile_csv` alone should suffice).
Not fixed here.

**Environment note:** a fresh git worktree does not receive the gitignored `extern/toml11`
vendored header; it was symlinked from the main checkout to build/run the parse tests locally.
Document this in the worktree setup, or the config suite silently no-ops.

---

## Summary
- Critical issues: 1 — **R-401, found by execution and fixed** (sub-FEC dimension; would have
  SIGSEGV'd every default-ON run). This is the headline: it validates Round 3's insistence
  that Phase 4 (build+run) was the only way to catch trace/collective defects.
- Moderate issues: 1 — R-404, **pre-existing and unrelated** (material `depth_profile_1d`
  test/parser mismatch blocking the full TOML parse suite). Noted, not fixed.
- Low issues: 2 — R-402 (movability hardened), R-403 (empty-local-rank coverage gap; code
  correct by design, follow-up test suggested).
- Plan compliance: **FULL** for Phases 1–4 as built and tested; the plan's prescribed sub-FEC
  `dim=3` was a defect corrected here (R-401).
- Verdict: **PASS.** All free-surface code builds; the class test passes at np=1 and np=2
  (trace ≤ 1e-12), the 5 parse tests pass, the driver compiles, the smoke script is wired.
  Two non-blocking follow-ups remain: R-403 (stronger empty-local coverage) and R-404
  (pre-existing material-test reconciliation — owner decision).

## Unreviewed Areas
- **Live Frontera behavior** at production np (hundreds of ranks; genuinely interior,
  top-face-free ranks): out of scope per `feedback-no-local-mesh-runs`; the np=2 run + R-401
  fix substantially de-risk it, but R-403's stronger coverage or a Frontera smoke is the final
  check before scale-up.
- **HDF5 (`Mode::Hdf5`) back end** of the slice: not exercised (VTU is the default and the
  tested path); the branch mirrors the volume writer and is `MFEM_USE_HDF5`-guarded.
- The smoke-script Check F's `velocity`/`mpi_rank` array detection was validated by logic +
  `py_compile`, not against a real run's `.vtu` (no production run locally, per policy).
