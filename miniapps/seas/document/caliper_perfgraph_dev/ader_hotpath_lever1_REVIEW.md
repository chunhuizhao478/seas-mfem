# Code Review: ADER hot-path Lever 1 implementation (fresh audit) — 2026-05-31

Round-2 adversarial review of the Lever-1 code (DerivMode `Cached` derivative
cache), re-run from scratch on the on-disk files. Honest headline up front: the
math is **correct** — the equivalence test shows Cached == OnTheFly to ≤2.5e-15
and the default path is byte-exact (bimaterial parity 9/9). There is **no
critical correctness bug**. The findings are real coverage/wiring/robustness gaps
that must be closed before the optimization can be trusted in the production
(recursion + parallel) path or deliver any speedup on Frontera.

## Review Scope
- Plan: `ader_hotpath_optimization_plan_2026-05-31.md` (+ `..._REVIEW.md`, `..._fix.md`)
- Files reviewed: `dynamic/elem_derivative_cache.hpp`, `dynamic/wave_operator.{hpp,inl}`
  (DerivMode/Cached branch/SetDerivMode), `tests/unit/test_wave_operator_spatial_derivative.cpp`,
  `tests/unit/bench_ader_hotpath.cpp`, `Makefile`, `drivers/spatial_dyn_driver.cpp`
- Domain context: `miniapps/seas/CLAUDE.md` (extreme-care files, no-local-mesh),
  memories (no-local-full-mesh, makefile-no-header-deps)

## Findings

### [R-001] [MODERATE] [wave_operator.inl — ApplySpatialDerivative/AdvanceADER] — the Cached path is correctness-tested ONLY at the single-call level; the production CK recursion in Cached mode is never asserted, only timed.

**Category:** EDGE_CASE / ASSUMPTION (test coverage)

**Description:**
`TestCachedEquivalence` (test_wave_operator_spatial_derivative.cpp:158-164)
toggles DerivMode around a bare `ApplySpatialDerivative` only. The PRODUCTION
consumers are the CK recursion — `ComputeADERSubStepStates`,
`ComputeADERTimeIntegrated`, `AdvanceADER` — which call `ApplySpatialDerivative`
(order−1)×3 times each into member scratch buffers. NO test runs those in Cached
mode and compares to OnTheFly. The existing `seas_test_ader_ck_predictor` /
`seas_test_ader_linear_wave_equivalence` run the DEFAULT (OnTheFly) only. So a
Cached-mode bug that manifests only through the recursion (buffer sizing/aliasing,
a stale-cache read, wrong scratch reuse) would pass every current test — the
bench (bench_ader_hotpath.cpp) exercises the recursion but only **times** it and
checksums the LAST mode, never comparing the two.

**Trigger:** any defect in the cached kernel that surfaces only under repeated
in-recursion calls.

**Actual behavior:** the recursion's Cached correctness is unverified.

**Expected behavior:** an end-to-end gate `AdvanceADER`(Cached) ≈ (OnTheFly) and
`ComputeADERSubStepStates`(Cached) ≈ (OnTheFly) to ≤1e-12 relative.

**Suggested fix:** extend `TestCachedEquivalence` (after the per-`dir` block):
```diff
+      // R-001: end-to-end recursion equivalence (production path).
+      wave.SetAbsorbingBackground(/*Q_bg=*/std::array<real_t,NUM_STATE>{}.data());
+      const real_t dt = 1e-4;
+      Vector Q2(N); FillRandomQ(Q2, 4242u + order);
+      auto adv = [&](DerivMode mode){ wave.SetDerivMode(mode); Vector qn;
+                                      wave.AdvanceADER(Q2, dt, order, qn); return qn; };
+      Vector qn_otf = adv(DerivMode::OnTheFly), qn_cac = adv(DerivMode::Cached);
+      wave.SetDerivMode(DerivMode::OnTheFly);
+      double mr=0,md=0; for(int i=0;i<N;i++){mr=std::max(mr,std::abs(qn_otf[i]));
+                          md=std::max(md,std::abs(qn_otf[i]-qn_cac[i]));}
+      TEST_LE(md/(mr+1e-300), 1e-12, "AdvanceADER cached==onthefly (rel)");
```
(plus the analogous `ComputeADERSubStepStates` `Q_per_node` comparison).

**Test case:**
```python
def test_R001_advance_ader_cached_equiv():
    qn_otf = advance_ader(Q, dt, order, mode="onthefly")
    qn_cac = advance_ader(Q, dt, order, mode="cached")
    assert rel_err(qn_otf, qn_cac) <= 1e-12      # recursion-level, not single-call
```

---

### [R-002] [MODERATE] [drivers/spatial_dyn_driver.cpp:1107] — Lever 1 has NO production activation path; the driver never enables Cached, so the measured 4–5× is unrealized on Frontera.

**Category:** DEVIATION (incomplete wiring)

**Description:**
`grep SetDerivMode` matches only `wave_operator.{hpp,inl}` + the test + the bench.
The driver builds `wave_ptr` (spatial_dyn_driver.cpp:1072/1104) and aliases
`WaveOperator<ParMesh> &wave = *wave_ptr;` (1107) but never calls
`wave.SetDerivMode(DerivMode::Cached)`, and there is no `--deriv-cache` CLI / TOML
knob. With the default `OnTheFly`, every TPV31 production run is byte-identical to
before — i.e. Lever 1 currently delivers ZERO speedup where it matters. This is
"staged-off by design," but the activation step is unbuilt and untracked, so the
effort is inert until it lands.

**Trigger:** running `seas_spatial_dyn_driver` (any TPV31 job) — it uses OnTheFly.

**Actual behavior:** no way to turn Cached on outside a unit binary.

**Expected behavior:** an opt-in (default still OnTheFly) that calls
`SetDerivMode(Cached)` right after operator construction, e.g. `--deriv-cache`.

**Suggested fix (driver, after line 1107):**
```diff
+   // Lever 1 (opt-in; default OnTheFly preserves byte-exact behaviour).
+   if (HasFlag(argc, argv, "--deriv-cache"))
+   {
+      wave.SetDerivMode(seas::DerivMode::Cached);
+      if (Mpi::Root()) { std::cout << "[deriv] Cached D_d^e enabled\n"; }
+   }
```

**Test case:**
```python
def test_R002_driver_exposes_deriv_cache_flag():
    # dry-run/dispatch banner with --deriv-cache reports "Cached D_d^e enabled";
    # without it, the operator stays OnTheFly.
    assert "Cached" in run_driver("--deriv-cache --verify-dispatch")
    assert "Cached" not in run_driver("--verify-dispatch")
```

---

### [R-003] [MODERATE] [POSSIBLE] [test_wave_operator_spatial_derivative.cpp] — Cached path has no ParMesh (parallel) coverage; production is `WaveOperator<ParMesh>`.

**Category:** EDGE_CASE

**Description:**
`TestCachedEquivalence` constructs serial `WaveOperator<Mesh>` /
`BimaterialWaveOperator<Mesh>` only. Production runs `WaveOperator<ParMesh>`.
`BuildElementDerivativeOperators` takes `*fes_` as `const FiniteElementSpace&`
(ParFiniteElementSpace binds to the base) and the cached kernel is element-LOCAL
(no collective), so it *should* be correct in parallel — but it is unverified,
and because R-002 leaves Cached unwired, the first parallel Cached execution would
be on Frontera with no prior test. Combined with R-001, the parallel recursion is
doubly uncovered.

**Trigger:** `WaveOperator<ParMesh>` with Cached at np>1.

**Suggested fix:** add an `np=2` parallel equivalence test (mirror an existing
`tests/parallel/` harness): build a `ParMesh`, assert Cached==OnTheFly for
`ApplySpatialDerivative` AND `AdvanceADER` to ≤1e-12 on each rank.

**Test case:**
```python
def test_R003_parmesh_cached_equiv_np2():
    # mpirun -np 2: per-rank max|cached - onthefly| <= 1e-12 for ApplySpatialDerivative
    run_mpi("seas_test_wave_operator_cached_parallel", np=2)
```

---

### [R-004] [LOW] [wave_operator.inl — ApplySpatialDerivative Cached branch] — missing the homogeneous-order guard the OnTheFly path has; a heterogeneous mesh corrupts silently instead of aborting.

**Category:** ASSUMPTION

**Description:**
OnTheFly asserts `MFEM_VERIFY(ndof == ndof_per_el_, ...)` per element. The Cached
branch uses `const int ndof = D.Height();` with `dof_offset = e * ndof_per_el_`
and NO such guard. On a (contract-forbidden) variable-`GetDof()` mesh, the cached
path would read/write wrong global offsets silently while OnTheFly aborts loud.
The class is homogeneous-only, so this is latent, but the two kernels must be
equally defensive.

**Suggested fix:**
```diff
       const DenseMatrix &D = elem_deriv_op_[e][dir];
       const int ndof = D.Height();
+      MFEM_VERIFY(ndof == ndof_per_el_,
+                  "ApplySpatialDerivative(Cached) assumes homogeneous elements: "
+                  "elem=" << e << " ndof=" << ndof << " != " << ndof_per_el_);
       const int dof_offset = e * ndof_per_el_;
```

**Test case:** n/a (would require a heterogeneous-order mesh the operator rejects
at construction; covered by inspection + symmetry with the OnTheFly guard).

---

### [R-005] [LOW] [wave_operator.inl — SetDerivMode] — labelled "Idempotent" but rebuilds the whole cache on every call; a per-step caller pays O(ne·ndof³) each step.

**Category:** QUALITY (latent performance footgun)

**Description:**
`SetDerivMode(Cached)` unconditionally calls `BuildElementDerivativeOperators`
every invocation (the header docstring claims "Idempotent"). Result is idempotent
but the work is not — calling it inside the time loop (an easy mistake when wiring
R-002) would rebuild the cache every macro-step, dwarfing the speedup.

**Suggested fix (early-out):**
```diff
 void WaveOperator<MeshType>::SetDerivMode(DerivMode m)
 {
+   const bool cache_ready =
+      (static_cast<int>(elem_deriv_op_.size()) == ne_);
+   if (m == deriv_mode_ &&
+       (m == DerivMode::OnTheFly || cache_ready)) { return; }   // true no-op
    if (m == DerivMode::Cached)
    {
```

**Test case:**
```python
def test_R005_setderivmode_idempotent_no_rebuild():
    wave.SetDerivMode(Cached); addr1 = id(wave.elem_deriv_op_[0][0].data)
    wave.SetDerivMode(Cached); addr2 = id(wave.elem_deriv_op_[0][0].data)
    assert addr1 == addr2   # second call must not reallocate/rebuild
```

---

## Summary
- Critical issues: 0 (the cached kernel is numerically correct; default path byte-exact)
- Moderate issues: 3 (R-001 recursion-level test gap, R-002 no production wiring,
  R-003 no parallel coverage)
- Low issues: 2 (R-004 missing guard, R-005 rebuild-not-early-out)
- Plan compliance: PARTIAL — Lever 1 kernel + equivalence (single-call) + R-004
  budget DONE; the recursion/parallel correctness gate and the production
  activation path are missing.
- Verdict: **PASS WITH FIXES** — the implementation is correct and safe-by-default,
  but R-001 (recursion-level equivalence test) and R-002 (driver wiring) must land
  before Lever 1 can be trusted in/around the production path or show any Frontera
  speedup. R-003–R-005 close the remaining gaps.

## Unreviewed Areas
- The R-004 abort path (over-budget `SetDerivMode(Cached)`) is code-correct but
  not unit-tested (the harness has no death-test facility) — verified by inspection.
- Lever 2 (merged CK recursion) is not implemented yet; out of scope for this
  Lever-1 audit.
