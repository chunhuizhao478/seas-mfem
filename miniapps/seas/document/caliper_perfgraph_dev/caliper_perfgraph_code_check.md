# Code Review: Caliper perfgraph IMPLEMENTATION (round 2 — implemented code)

> **Location/name:** the `/code-review` skill specifies `REVIEW.md` in the
> project root, but `miniapps/seas/REVIEW.md` already exists and would be
> clobbered. This is the code-review companion (distinct from the plan-review
> `caliper_perfgraph_check.md`). The `/code-fix` agent should read THIS file.

## Review Scope
- Plan: `document/caliper_perfgraph_dev/caliper_perfgraph_plan.md`
- Files reviewed (the implemented diff — 12 insertions, 0 deletions):
  `dynamic/wave_operator.inl`, `dynamic/bimaterial_wave_operator.inl`,
  `drivers/spatial_dyn_driver.cpp`, `document/caliper_perfgraph_dev/README_caliper_build_and_run.md`
- Domain context: `miniapps/seas/CLAUDE.md`; the ADER vs RK call structure in
  `wave_operator.inl`.
- Verified-correct (NOT findings): diff is annotations-only (12 added / 0
  removed, all `MFEM_PERF_*`); all 8 region names match the README table;
  `friction_substep` BEGIN/END is balanced (both inside `if (n_total_fault_qps
  > 0)`, no return/break between them); 6 wave-op/ADER unit tests pass on the
  default (no-op) build.

## Findings

### [R-001] [MODERATE] [wave_operator.inl — ADER flux methods] — The annotated flux regions are RK-path-only; TPV31's actual (ADER) face-flux cost is unannotated, so the perfgraph for the primary use case has no flux breakdown

**Category:** DEVIATION (coverage gap vs the stated goal: a useful TPV31 perfgraph)

**Description:**
TPV31 runs the **ADER** path, not the RK/`Mult` path. The implemented
annotations cover `Mult` (`:716`), `ComputeFaceFluxRHS` (`:2284`), and
`ComputeSharedFaceFluxRHS` (`:2934`) — but those are reached **only** through
`WaveOperator::Mult`, which the driver calls **only on the RK path**
(`spatial_dyn_driver.cpp:2756` "RK4/RK45 coupled stepper: drives wave.Mult
directly"). The ADER macro-step (`AdvanceADER`, `:5283`) instead calls, at
`:5327-5332`:
- `ComputeADERVolumeUpdate(I, rhs)` (`:2260`) → which delegates to the annotated
  `ComputeVolumeRHS` (so VOLUME cost *is* captured), but
- `ComputeADERFaceFluxRHS(I, dt, rhs)` (`:3466`) — **NOT annotated**, and
- `ComputeADERSharedFaceFluxRHS(I, dt, rhs)` (`:4486`) — **NOT annotated**.

Net effect on a Caliper build of TPV31: three annotated regions
(`seas::WaveOperator::Mult`, `::ComputeFaceFluxRHS`, `::ComputeSharedFaceFluxRHS`)
**never fire**, and the dominant DG cost — interior + shared face-flux assembly
— is invisible (folded opaquely into `seas::WaveOperator::AdvanceADER`). The
perfgraph cannot answer "how much time in face-flux assembly?" for the exact
benchmark this feature targets.

**Trigger:**
Caliper build, TPV31 (ADER) run, inspect `runtime-report`.

**Actual behavior:**
`AdvanceADER` shows volume (`ComputeVolumeRHS`) but no face-flux child;
`Mult`/`ComputeFaceFluxRHS`/`ComputeSharedFaceFluxRHS` are absent (0 calls).

**Expected behavior:**
The ADER face-flux assembly is a named region under `AdvanceADER`.

**Suggested fix:** annotate the two ADER face methods (keep the existing
Mult-path annotations — they are correct for RK runs).
```diff
  // wave_operator.inl  ~:3469 (ComputeADERFaceFluxRHS body opening)
  void WaveOperator<MeshType>::ComputeADERFaceFluxRHS(const Vector &I,
                                                      real_t dt,
                                                      Vector &rhs) const
  {
+    MFEM_PERF_SCOPE("seas::WaveOperator::ComputeADERFaceFluxRHS");
     MFEM_VERIFY(dt > 0.0,
                 "ComputeADERFaceFluxRHS: dt must be > 0, got " << dt);
```
```diff
  // wave_operator.inl  ~:4489 (ComputeADERSharedFaceFluxRHS body opening)
  void WaveOperator<MeshType>::ComputeADERSharedFaceFluxRHS(const Vector &I,
                                                            real_t dt,
                                                            Vector &rhs) const
  {
+    MFEM_PERF_SCOPE("seas::WaveOperator::ComputeADERSharedFaceFluxRHS");
     if constexpr (!IsParallelMesh<MeshType>::value)
```
(Optional but recommended for a complete ADER picture: also annotate
`ComputeADERVolumeUpdate` (`:2260`) so the ADER volume frame is explicit rather
than appearing as a bare `ComputeVolumeRHS` under `AdvanceADER`.)

**Test case (verification check — Caliper build required):**
```bash
# TPV31 (ADER) runtime-report MUST contain the ADER face regions:
grep -q 'seas::WaveOperator::ComputeADERFaceFluxRHS'       runtime_report.txt || { echo "FAIL R-001: ADER interior face flux unprofiled"; exit 1; }
grep -q 'seas::WaveOperator::ComputeADERSharedFaceFluxRHS' runtime_report.txt || { echo "FAIL R-001: ADER shared face flux unprofiled"; exit 1; }
```

---

### [R-002] [LOW] [README_caliper_build_and_run.md — "What was instrumented"] — Nesting diagram is wrong for the ADER path (claims `Mult` and a `step`-level `friction_substep`)

**Category:** QUALITY (documentation inaccuracy that will mislead graph interpretation)

**Description:**
The README states the nesting is
`step → AdvanceADERWithSubStep → {ComputeADERSubStepStates, AdvanceADER (→ Mult → {ComputeVolumeRHS, ComputeFaceFluxRHS, ComputeSharedFaceFluxRHS})}, with friction_substep as a sibling under step`.
Two errors for the TPV31/ADER path: (1) `Mult` is **never** on the ADER path
(RK-only), and `ComputeFaceFluxRHS`/`ComputeSharedFaceFluxRHS` are reached only
via `Mult` — so the drawn chain `AdvanceADER → Mult → {...flux...}` does not
occur for TPV31; (2) `friction_substep` is **not** a sibling under `step` — the
`iterator.Advance` call it brackets is inside `AdvanceADERWithSubStep_Spatial`
(driver `:449`, within the `:375`-function), so it nests under
`AdvanceADERWithSubStep`, not `step`.

**Trigger:** reading the README to interpret a TPV31 perfgraph.

**Actual / Expected:** the diagram should reflect the ADER reality (and the new
regions from R-001):
```
step
└─ AdvanceADERWithSubStep
   ├─ ComputeADERSubStepStates
   ├─ friction_substep
   └─ AdvanceADER
      ├─ ComputeVolumeRHS            (via ComputeADERVolumeUpdate)
      ├─ ComputeADERFaceFluxRHS       (after R-001)
      └─ ComputeADERSharedFaceFluxRHS (after R-001)
ComputeMaxDt   (once, at setup — top level)
Mult → {ComputeVolumeRHS, ComputeFaceFluxRHS, ComputeSharedFaceFluxRHS}  (RK path only; absent for TPV31)
```

**Suggested fix:** replace the README "expected nesting" paragraph with the tree
above, and add a one-line note that `Mult`/`ComputeFaceFluxRHS`/
`ComputeSharedFaceFluxRHS` appear only on `--time-integrator rk4|rk45` runs.

**Test case:** n/a (documentation).

---

### [R-003] [LOW] [POSSIBLE] [wave_operator.inl — ComputeADERSubStepStates internals] — ADER predictor's spatial-derivative work (`ApplySpatialDerivative`) is unprofiled

**Category:** ASSUMPTION (depth-of-coverage)

**Description:**
`ComputeADERSubStepStates` (annotated, `:1198`) drives the ADER Cauchy–Kovalewski
predictor, whose per-direction inner kernel is `ApplySpatialDerivative`
(`:892`), called in the predictor loop. It is not separately annotated, so the
predictor shows as one opaque block. This is POSSIBLE/optional: the predictor's
total IS captured by `ComputeADERSubStepStates`, and `ApplySpatialDerivative`
is finer-grained (called O(order × dim) per sub-step) — annotating it is a
reasonable next level of detail but adds modest region overhead. Flagged so the
fix agent can decide; not required for a first-cut perfgraph.

**Suggested fix (optional):**
```diff
  // wave_operator.inl  ~:892 (ApplySpatialDerivative body opening)
  void WaveOperator<MeshType>::ApplySpatialDerivative(int dir, ...) const
  {
+    MFEM_PERF_SCOPE("seas::WaveOperator::ApplySpatialDerivative");
     ...
```
(Verify the call frequency first — if it is invoked per-element inside a loop
rather than once per (direction, sub-step), DO NOT annotate it, per the plan's
granularity rule.)

**Test case:** n/a (optional depth).

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — ADER face-flux cost unprofiled for the primary
  TPV31 use case)
- Low issues: 2 (R-002 README nesting wrong; R-003 optional predictor depth)
- Plan compliance: **FULL vs the plan as written** — every annotation the plan
  listed was implemented correctly and the diff is annotations-only. The gap is
  that the **plan's region list itself** targeted the RK/`Mult` flux methods,
  not the ADER flux methods TPV31 actually executes; the implementation
  faithfully inherited that gap. R-001 corrects it at the code level.
- Verdict: **PASS WITH FIXES** — apply R-001 (and the dependent R-002 doc
  update) so the TPV31 perfgraph actually shows flux cost. R-003 is optional.

## Unreviewed Areas
- Caliper-ON behavior was not executed (no Caliper in this checkout; no local
  full-mesh runs per project rule). R-001's "regions never fire" conclusion is
  derived from the static call graph (`AdvanceADER:5327-5332` calls the ADER
  flux methods; `Mult` is RK-only per driver `:2756`), not observed at runtime.
- The RK path itself (`Mult` → annotated flux methods) was not exercised; those
  annotations are presumed correct for `--time-integrator rk4|rk45` but
  untested here (TPV31 does not use them).
