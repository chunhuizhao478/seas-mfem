# Implementation Plan: Caliper "perfgraph" profiling for the TPV31 dynamic-rupture path

## Overview
Add MFEM Caliper region annotations (`MFEM_PERF_SCOPE`) to the dynamic-rupture
hot loops so that, on a `MFEM_USE_CALIPER=YES` build, `seas_spatial_dyn_driver`
emits a per-region timing tree (a "perfgraph") via the runtime `CALI_CONFIG`
env var. On the default build (`MFEM_USE_CALIPER=NO`) every annotation is a
compile-time no-op, so numerics and performance are bit-for-bit unchanged. The
plan also documents the Caliper-enabled build recipe and the `CALI_CONFIG`
runtime usage (the sbatch hook `TPV31_PERFGRAPH=1` is already half-wired).

Work happens in worktree `worktree-feature-caliper-perfgraph` (branched from
`system/spatial_dyn_driver @ e0f8f15`). **No source logic changes** — only
annotation macros + documentation.

## Background (verified facts)
- `general/annotation.hpp` defines `MFEM_PERF_FUNCTION`, `MFEM_PERF_BEGIN(s)`,
  `MFEM_PERF_END(s)`, `MFEM_PERF_SCOPE(name)`. Under `#ifdef MFEM_USE_CALIPER`
  they wrap Caliper (`CALI_CXX_MARK_*`); otherwise each is **defined empty**
  (`annotation.hpp:28-32`) → a bare `;` after macro expansion.
- `mfem.hpp:29` includes `general/annotation.hpp`. Both
  `dynamic/wave_operator.hpp:15` and `drivers/spatial_dyn_driver.cpp:47`
  include `mfem.hpp`, and `dynamic/bimaterial_wave_operator.hpp` includes
  `wave_operator.hpp`. **The macros are already in scope in all three target
  files — no new `#include` is required** (Phase implementers MUST confirm,
  not assume).
- Caliper is OFF in this checkout: `config/defaults.mk:181` `MFEM_USE_CALIPER =
  NO`; build vars `CALIPER_DIR/CALIPER_OPT/CALIPER_LIB` at
  `config/defaults.mk:544-555`.
- `BimaterialWaveOperator` (matrix path) inherits `Mult`, `AdvanceADER`,
  `ComputeADERSubStepStates`, `ComputeVolumeRHS`, `ComputeFaceFluxRHS`,
  `ComputeSharedFaceFluxRHS` from `WaveOperator` unchanged; it overrides only
  `InteriorFaceFlux_`, `SharedInteriorFaceFlux_`, `ApplyElementJacobian_`,
  `ComputeMaxDt`. So annotating the **base** methods profiles BOTH the scalar
  (p1/O2 job) and matrix (TPV31) paths.

## Constraints
- **Files requiring extreme care** (`miniapps/seas/CLAUDE.md`):
  `dynamic/wave_operator.{hpp,inl}` and `drivers/spatial_dyn_driver.cpp`. Any
  change requires the full unit-test verification below.
- **Zero numerical/perf impact on the default build.** The only permitted edits
  are inserting annotation macros: `MFEM_PERF_SCOPE("...")` at method/loop-body
  scope, and the one `MFEM_PERF_BEGIN("...")`/`MFEM_PERF_END("...")` pair around
  the friction call site (Phase 2 item 3). No reordering, no new variables in
  hot paths, no control-flow changes. On `MFEM_USE_CALIPER=NO` every such macro
  is empty, so the emitted object code MUST be identical.
- **Do NOT edit `general/annotation.hpp` or `config/defaults.mk`.** The Caliper
  switch is set at configure time (`make config MFEM_USE_CALIPER=YES`), never by
  committing `=YES` (that would force every build to need libcaliper).
- **No full-mesh runs locally** (memory `feedback-no-local-mesh-runs`). Local
  verification = compile + unit tests only. The Caliper-ON build and actual
  perfgraph emission are validated on Frontera (documented, not run here).
- **Granularity rule:** annotate only methods called O(steps)·O(substeps)
  times. Do **NOT** annotate per-face / per-element hooks
  (`InteriorFaceFlux_`, `SharedInteriorFaceFlux_`, `ApplyElementJacobian_`,
  `FluxForElem_`) — they are called inside the assembly face/element loops
  (`wave_operator.inl:2851, 3417, 4399, 5155`), i.e. O(faces) ≈ millions of
  times per run; a Caliper region there adds push/pop overhead per call and
  distorts the very profile we want. The per-step aggregate methods
  (`ComputeFaceFluxRHS`/`ComputeSharedFaceFluxRHS`) capture the flux-assembly
  cost without that overhead.
- **Region-name convention:** `"seas::<Class>::<Method>"` (string literal passed
  to `MFEM_PERF_SCOPE`), so the perfgraph node labels are stable and readable.

## Phase 1: Annotate WaveOperator / BimaterialWaveOperator hot methods

### Goal
Every coarse-grained dynamic-rupture compute method emits a named Caliper region
on a Caliper build; the default build is byte-unchanged.

### Files to Modify
- `dynamic/wave_operator.inl` — insert `MFEM_PERF_SCOPE(...)` as the **first
  statement** in each method body (before the existing `MFEM_VERIFY` guards;
  the RAII scope then spans the whole body).
- `dynamic/bimaterial_wave_operator.inl` — annotate the matrix-path
  `ComputeMaxDt` override only.

### Detailed Requirements
Insert exactly these, one per method, as the first line inside `{`:

| File | Line (current) | Method | Insert |
|------|----------------|--------|--------|
| wave_operator.inl | 716 | `Mult` | `MFEM_PERF_SCOPE("seas::WaveOperator::Mult");` |
| wave_operator.inl | 804 | `ComputeVolumeRHS` | `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeVolumeRHS");` |
| wave_operator.inl | 1198 | `ComputeADERSubStepStates` | `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeADERSubStepStates");` |
| wave_operator.inl | 2284 | `ComputeFaceFluxRHS` | `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeFaceFluxRHS");` |
| wave_operator.inl | 2934 | `ComputeSharedFaceFluxRHS` | `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeSharedFaceFluxRHS");` |
| wave_operator.inl | 5283 | `AdvanceADER` | `MFEM_PERF_SCOPE("seas::WaveOperator::AdvanceADER");` |
| wave_operator.inl | 5530 | `ComputeMaxDt` | `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeMaxDt");` |
| bimaterial_wave_operator.inl | 515 | `ComputeMaxDt` (override) | `MFEM_PERF_SCOPE("seas::BimaterialWaveOperator::ComputeMaxDt");` |

(Line numbers are advisory — match on the method signature, which is unique, not
the line.)

### Explicitly NOT annotated (and why)
- `WaveOperator::InteriorFaceFlux_` (1053), `SharedInteriorFaceFlux_` (1073),
  and the `BimaterialWaveOperator` overrides `InteriorFaceFlux_` (471),
  `SharedInteriorFaceFlux_` (489), `ApplyElementJacobian_` (501): per-face /
  per-element, called millions of times — see the Granularity constraint. Their
  cost is already captured by the enclosing `ComputeFaceFluxRHS` /
  `ComputeSharedFaceFluxRHS` regions.

### Edge Cases to Handle
- Methods are `const` and templated: `MFEM_PERF_SCOPE` expands to a local RAII
  guard (or to nothing) — valid in both contexts. No `mutable`/`this` issues.
- Early returns / `MFEM_ABORT` inside a method: the RAII region still closes
  correctly on scope exit (Caliper) or is absent (no-op build). No special
  handling.
- `MFEM_PERF_SCOPE("x");` on a non-Caliper build becomes `;` — confirm it
  compiles with `-Werror` (it does: empty statement).

### Acceptance Criteria
- [ ] Default build (`MFEM_USE_CALIPER=NO`) compiles with no new warnings.
- [ ] All 8 regions inserted exactly as listed; no other edits to the two files
      (verify with `git diff --stat` = 2 files; `git diff` shows only inserted
      `MFEM_PERF_SCOPE` lines).
- [ ] `git diff` contains no change to any `InteriorFaceFlux_` /
      `SharedInteriorFaceFlux_` / `ApplyElementJacobian_` / `FluxForElem_` body.
- [ ] Local unit tests (see "Local Unit-Test Verification") pass identically to
      the pre-change baseline.

### Dependencies
- Depends on: nothing.
- Required by: Phase 3 (build recipe references these region names).

## Phase 2: Annotate the spatial driver step loop + macro-step function

### Goal
The driver's per-macro-step work, the ADER sub-step orchestration, AND the
friction sub-step solve each emit a Caliper region, giving the perfgraph its
top-level structure (`step` → `AdvanceADERWithSubStep` →
{`friction_substep`, WaveOperator regions from Phase 1}). Separating the
friction solve from wave propagation is the primary reason to profile a coupled
dynamic-rupture solver, so it is a first-class region, not folded into the
macro-step wrapper.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`:
  1. **`AdvanceADERWithSubStep_Spatial`** (free function, opens at line 375;
     body `{` at line ~386): insert as the first statement in the body
     `MFEM_PERF_SCOPE("seas::spatial_dyn::AdvanceADERWithSubStep");`
  2. **Main time loop** (`for (int step = step0; step < nsteps; ++step)` at
     line 2727; body `{` at line 2728): insert as the **first statement inside
     the loop body** `MFEM_PERF_SCOPE("seas::spatial_dyn::step");`
     (a fresh RAII region per iteration → Caliper aggregates "step" across all
     iterations and nests the per-step callees under it).
  3. **Friction sub-step solve** — the per-QP LSW/RS solve inside
     `iterator.Advance(...)` (call site at line ~449), a primary coupled-system
     hot path currently invisible in the graph. Bracket the call:
     insert `MFEM_PERF_BEGIN("seas::spatial_dyn::friction_substep");`
     immediately BEFORE the `iterator.Advance(dof_data, fault_coords, ...);`
     statement and `MFEM_PERF_END("seas::spatial_dyn::friction_substep");`
     immediately AFTER it.
     - Use BEGIN/END (not `MFEM_PERF_SCOPE`) because this brackets a single call
       site, not a function body. `iterator.Advance` failing is currently fatal,
       so the non-RAII "region left open on throw" risk is moot.
     - Exception-safe alternative (more edits, equivalent graph): put
       `MFEM_PERF_SCOPE("seas::FrictionIterator::Advance");` at the top of each
       concrete `Advance` body in
       `dynamic/{tpv205,tpv104,tpv102}_substep_iterator.cpp`. TPV31 dispatches
       to the LSW (tpv205-style) iterator; annotate all three for completeness.
       Pick ONE approach, not both (double-counting).

### Edge Cases to Handle
- The loop body has an early `if (dt_step <= 0.0) { break; }` (line 2730): place
  the `MFEM_PERF_SCOPE` **before** it so the region opens at iteration start; the
  `break` exits the scope cleanly (region closes for that iteration).
- `AdvanceADERWithSubStep_Spatial` is only reached on the ADER path; the RK path
  (not used by TPV31) is separate and out of scope for this plan.

### Acceptance Criteria
- [ ] Default build compiles, no new warnings.
- [ ] **Sole behavioral gate (mandatory):** there is NO unit test over
      `spatial_dyn_driver.cpp`; the 6 tests below cover the WaveOperator library
      only. Phase 2 is therefore verified ONLY by (a) clean compile + link and
      (b) the reviewer inspecting `git diff drivers/spatial_dyn_driver.cpp`
      line-by-line. A misplaced macro would still compile/link/pass those tests,
      so the diff inspection is the real guard — do not skip it.
- [ ] `git diff drivers/spatial_dyn_driver.cpp` shows exactly the two inserted
      `MFEM_PERF_SCOPE` lines (items 1–2) PLUS the
      `MFEM_PERF_BEGIN`/`MFEM_PERF_END` pair around `iterator.Advance` (item 3),
      and nothing else. (If the exception-safe alternative for item 3 was
      chosen, the diff is instead in the `*_substep_iterator.cpp` files —
      one `MFEM_PERF_SCOPE` line each — and the driver shows only items 1–2.)
- [ ] `seas_spatial_dyn_driver` builds (default, Caliper OFF) from the worktree.
- [ ] `--dry-run --verify-dispatch` behavior is unchanged (not run on a full
      mesh; covered by the driver building + the existing dispatch unit logic).

### Dependencies
- Depends on: nothing (independent of Phase 1, but share the same verification).
- Required by: Phase 3.

## Phase 3: Caliper build recipe + runtime usage + sbatch reconciliation

### Goal
A reproducible recipe to (a) build MFEM+seas with Caliper, (b) emit a perfgraph
at runtime, documented in this `caliper_perfgraph_dev/` folder; and the sbatch
hook reconciled into the worktree.

### Files to Create
- `document/caliper_perfgraph_dev/README_caliper_build_and_run.md` — the recipe
  below, verbatim, as the user-facing runbook.

### Files to Modify
- (Reconciliation only) `jobs/tpv31_spatial/tpv31_p2_aderO3_noflux_50m_flex.sbatch`
  — the `TPV31_PERFGRAPH=1` → `CALI_CONFIG` block already exists in the
  `system/spatial_dyn_driver` working tree (uncommitted, created alongside this
  effort). It is **NOT** in this worktree's base commit. At merge time, ensure
  exactly one copy of that block lands. The implementer MUST NOT duplicate it;
  if the sbatch is absent in the worktree, leave a note in the fix report rather
  than recreating it (avoid divergent copies).

- **Merge target (REQUIRED):** this worktree branch
  (`worktree-feature-caliper-perfgraph`) MUST be merged into
  `system/spatial_dyn_driver`, where the `TPV31_PERFGRAPH` sbatch hook and the
  p1/p2 TPV31 jobs live (currently uncommitted on that branch's working tree —
  commit them there first). Merging anywhere else (e.g. a fresh branch off
  `origin/main`) separates the code annotations from their runnable sbatch hook
  and leaves the feature un-demonstrable end-to-end.

### Detailed Requirements — Caliper build recipe (documented, run on Frontera)
1. **Obtain Caliper** (Frontera): `spack install caliper +adiak` (or build from
   `github.com/LLNL/Caliper` with CMake). Record the install prefix as
   `CALIPER_DIR`. (Adiak optional; enables run metadata.)
2. **Reconfigure MFEM with Caliper** (in the MFEM build dir):
   `make config MFEM_USE_CALIPER=YES CALIPER_DIR=<prefix> [ADIAK_DIR=<prefix>]`
   then `make -j`. This sets `MFEM_USE_CALIPER` in the generated `config.hpp`,
   flipping every `MFEM_PERF_*` macro live and linking `-lcaliper`
   (`config/defaults.mk:544-555`).
3. **Rebuild the seas driver (force-clean — same no-header-deps hazard as the
   local verification):** `cd miniapps/seas && make clean && make seas_spatial_dyn_driver`.
   Skipping `make clean` here is the #1 cause of an EMPTY perfgraph: objects
   compiled before the Caliper flip expanded every `MFEM_PERF_*` to a no-op, so
   without a clean the driver relinks those stale objects — it links libcaliper
   but emits zero `seas::*` regions. (Same applies if rebuilding the p1/p2 jobs'
   driver.)
4. **Run with a Caliper config** (the sbatch sets these when `TPV31_PERFGRAPH=1`):
   - Text region tree to stderr: `CALI_CONFIG="runtime-report(calc.inclusive=true)"`
   - Call-graph JSON for Hatchet:
     `CALI_CONFIG="hatchet-region-profile(output=<OUT>/caliper_perfgraph.json)"`
     then render with the `hatchet` Python package
     (`GraphFrame.from_caliperreader(...)` → `.tree()` / `.to_dot()`).
   - Both can be combined comma-separated (as the sbatch block does).

### Acceptance Criteria
- [ ] `README_caliper_build_and_run.md` exists and the recipe references the
      exact region names from Phases 1–2.
- [ ] The recipe states the two prerequisites (Caliper build + the annotations
      from Phases 1–2) and that on a non-Caliper build the hook is inert.
- [ ] No duplicate `TPV31_PERFGRAPH` block (verified at merge).

### Dependencies
- Depends on: Phases 1 & 2 (region names).
- Required by: nothing.

## Local Unit-Test Verification (MANDATORY — runs after Phases 1 & 2)

This is the only local gate. It proves the annotated **default** build
(`MFEM_USE_CALIPER=NO`) is functionally identical to baseline.

### Build environment
- `conda activate mfem-dev`.
- We are in the worktree; `config/config.mk` and `libmfem.a` are **absent here**
  (build-generated, not committed). Per memory `worktree-build-mfem-dir-override`,
  build seas targets against the MAIN repo's MFEM artifacts:
  ```
  MAIN=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver
  cd /Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/feature-caliper-perfgraph/miniapps/seas
  ```
  All `make` invocations below take:
  `MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN`
  (keep `MFEM_DIR` RELATIVE `../..`; an absolute `MFEM_DIR` wrongly switches the
  source tree — memory `worktree-build-mfem-dir-override`).

### Force-clean rebuild (REQUIRED — template hazard)
`WaveOperator`/`BimaterialWaveOperator` are **templates** instantiated in every
test TU; their bodies live in the edited `.inl` files. The Makefile has **no
header dependencies** (memory `makefile-no-header-deps-stale-o`), so editing
`.inl`/`.cpp` will NOT recompile the dependent `.o`/test objects — a stale
instantiation links against the new ABI and SIGABRTs. Therefore **`make clean`
before building the tests** (do not rely on incremental make):
```
make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN clean
```

### Tests to build + run (each builds then `./run`s)
These directly exercise the annotated methods (Mult / AdvanceADER /
ComputeADERSubStepStates / flux RHS / ComputeMaxDt) on both scalar and matrix
paths:
```
rc=0
for T in \
    seas_test_bimaterial_wave_operator_parity \
    seas_test_wave_operator \
    seas_test_wave_operator_spatial_derivative \
    seas_test_ader_ck_predictor \
    seas_test_ader_linear_wave_equivalence \
    seas_test_ader_tpv102_smoke ; do
  make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN "$T" \
    && ./"$T" || { echo "FAIL: $T"; rc=1; break; }
done
[ "$rc" -eq 0 ] || { echo "VERIFICATION FAILED"; exit 1; }
```
And confirm the driver itself links (Phase 2 target):
```
make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN seas_spatial_dyn_driver
```

### Pass criteria
- [ ] `make clean` + all 6 unit tests build and exit 0 (same pass/fail set as
      the pre-change baseline on `e0f8f15`; the 3 known pre-existing failures in
      memory `preexisting-worktree-test-failures-2026-05` — if they appear — are
      NOT in this set, so this set should be fully green).
- [ ] `seas_spatial_dyn_driver` links.
- [ ] `git diff` is annotations-only (Phases 1–2 acceptance criteria).
- [ ] (Optional, strongest) byte-identical object proof: build
      `dynamic/wave_operator.o` before and after the edit on the default build
      and confirm only the annotation source lines differ — i.e. the macro is
      truly empty. (A clean rebuild + green tests is sufficient; this is extra.)

## Testing Strategy
- **Phase 1 & 2:** the unit-test set above is the regression gate — it covers
  the exact annotated methods on scalar + matrix + ADER predictor/corrector.
  Because the macros are no-ops on this build, any test delta would indicate an
  accidental logic edit, not the annotations.
- **Phase 3:** documentation review; the live Caliper build + perfgraph emission
  is a Frontera step (out of local scope) — the recipe is the deliverable, the
  run is operator-executed.
- **No new unit test is required** (annotations are no-ops locally; there is
  nothing testable about a no-op without a Caliper build). If a future
  Caliper-enabled CI lane is added, a smoke test asserting ≥1 region named
  `seas::WaveOperator::Mult` appears in the Caliper output would be the natural
  addition — noted as a follow-up, not in scope.

## Risk Assessment
- **Accidental logic edit in an "extreme care" file.** Mitigation: the only
  allowed change is inserting annotation macros (one `MFEM_PERF_SCOPE` per
  method/loop body; one `MFEM_PERF_BEGIN`/`MFEM_PERF_END` pair around the
  friction call site); `git diff` review + the unit-test set catch any
  deviation. Reviewer must reject any non-macro hunk.
- **Stale-`.o` ABI SIGABRT** (template + no header deps). Mitigation: mandatory
  `make clean` before the verification build (above). This is the most likely
  way to get a *false* test failure; do not skip it.
- **Per-face annotation creeping in** and distorting/inflating the profile.
  Mitigation: the explicit "NOT annotated" list + a `git diff` check that no
  per-face hook body changed.
- **Caliper region overhead even at coarse granularity** is negligible
  (O(steps·substeps) push/pops), but the recipe should note that
  `runtime-report` adds small overhead; for production timing runs use
  `hatchet-region-profile` (lower overhead) and compare against a
  no-`CALI_CONFIG` wall-clock baseline.
- **sbatch block duplication at merge.** Mitigation: Phase 3 reconciliation note
  — the `TPV31_PERFGRAPH` hook already exists in the main tree; do not recreate.
- **Caliper unavailable on Frontera.** Mitigation: the recipe gives both `spack`
  and from-source paths; if neither is feasible, the fallback (TACC REMORA, no
  code change) is recorded in the session history as option (B).
