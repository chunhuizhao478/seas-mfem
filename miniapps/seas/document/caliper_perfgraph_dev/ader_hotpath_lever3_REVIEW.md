# Code Review: ADER hot-path Lever 3 (cached volume-RHS operator) — 2026-05-31

Fresh adversarial audit of the Lever-3 code, re-run from scratch on the on-disk
files. Honest headline: the cached `ComputeVolumeRHS` is **numerically correct** —
the equivalence test shows cached == OnTheFly to ≤5.6e-16 (scalar + bimaterial),
the recursion test (`AdvanceADER`, which now calls the cached path) still passes
at 1.4e-15, Lever-2 bit-exact parity stays 0.000e+00, and the default path is
byte-exact (bimaterial parity 9/9). The `S_d^e` math was re-derived and matches;
`FluxForElem_(e)` is used (NOT the `(1,1,1)` placeholder — avoids the
`matrix-path-fluxforelem-placeholder-leak` class). **No critical correctness
bug.** The findings are coverage gaps — the top one matters because the
production benchmark (TPV31) exercises a path the tests do not.

## Review Scope
- Plan: `ader_hotpath_optimization_plan_2026-05-31.md` §4 Lever 3 (+ `..._lever3_impl.md`)
- Files reviewed: `dynamic/elem_derivative_cache.hpp` (`BuildElementVolumeOperators`),
  `dynamic/wave_operator.{hpp,inl}` (`elem_volume_op_`, Cached `ComputeVolumeRHS`,
  `SetDerivMode`), `tests/unit/test_wave_operator_spatial_derivative.cpp`,
  `tests/unit/bench_ader_hotpath.cpp`
- Domain context: `miniapps/seas/CLAUDE.md`, memories (placeholder-leak,
  no-local-full-mesh)

## Findings

### [R-001] [MODERATE] [POSSIBLE] [wave_operator.inl ComputeVolumeRHS Cached + test] — the cached path is validated ONLY with homogeneous material (constant `A_d^e`); the production TPV31 benchmark is heterogeneous (`depth_profile_1d`), where `A_d^e` varies per element — and that is exactly the material-dependent part of the new code.

**Category:** EDGE_CASE / ASSUMPTION (coverage)

**Description:**
`TestCachedEquivalence` builds the bimaterial operator with
`MaterialField::MakeCoefficient(ConstantCoefficient…)` (test line ~294-296), so
every element gets the SAME `A_d^e`. The cached `ComputeVolumeRHS` pairs each
element's `S_d^e` with `FluxForElem_(e).GetReferenceStarMatrix(d)`; with a
constant material, a defect that mis-pairs the per-element flux (e.g. wrong `e`,
or applying `A` before/after the wrong reduction) would be INVISIBLE because all
`A_d^e` are identical. The geometry `S_d^e` does vary across the test's tets, so
the geometry path is covered — but the **per-element material application is
not**. Production (`[material].kind="depth_profile_1d"`) has element-varying
`A_d^e`, so a heterogeneous-only bug would silently corrupt the actual run while
every current test stays green. The code *structurally* uses `FluxForElem_(e)`
correctly (low risk), but on a numerical-physics code targeting this exact
benchmark, the consequence (silent wrong physics) warrants the test.

**Trigger:** a Cached run on an element-varying material (TPV31 depth_profile_1d).

**Expected behavior:** cached `ComputeVolumeRHS` == OnTheFly to ≤1e-12 with
**varying** `A_d^e`.

**Suggested fix (add to `TestCachedEquivalence`, alongside the homogeneous case):**
```cpp
// R-001: heterogeneous material -> per-element-varying A_d^e.
FunctionCoefficient lam_v([](const Vector &x){ return 32.04e9*(1.0+0.5*x(0)); });
FunctionCoefficient mu_v ([](const Vector &x){ return 32.04e9*(1.0+0.3*x(1)); });
ConstantCoefficient  rho_v(2670.0);
BimaterialWaveOperator<Mesh> wh(
   mesh, order, MaterialField::MakeCoefficient(&lam_v, &mu_v, &rho_v), bc);
const int Nh = NUM_STATE * wh.GetFESpace().GetNDofs();
double vrelh = VolumeRHSCachedVsOnTheFlyMaxRelErr(wh, Nh, 8484u + order);
TEST_LE(vrelh, tol, "heterogeneous volume RHS (varying A_d): cached==onthefly");
```

**Test case:**
```python
def test_R001_volume_rhs_heterogeneous_material():
    wh = bimaterial_op(material=position_varying(lambda_, mu, rho))   # A_d^e varies per element
    assert rel_err(volrhs(wh, mode="cached"), volrhs(wh, mode="onthefly")) <= 1e-12
```

---

### [R-002] [MODERATE] [test coverage] — cached `ComputeVolumeRHS` has no ParMesh (parallel) coverage and no `Mult` (RK) path coverage; both are reachable in production with `--deriv-cache`.

**Category:** EDGE_CASE

**Description:**
(1) The equivalence test uses serial `Mesh`; production is `WaveOperator<ParMesh>`.
The cached volume kernel is element-local (no MPI), so it *should* be correct in
parallel, but it is unverified — same gap as Lever 1's deferred R-003.
(2) `ComputeVolumeRHS` is also called from `Mult` (`wave_operator.inl:725`, the
RK/instantaneous path), which the cached branch now alters when `--deriv-cache`
is set on an RK job. The equivalence test exercises `ComputeVolumeRHS_ForTest`
and `AdvanceADER` (ADER), but NOT `Mult`. So the Cached+RK volume path is untested.

**Trigger:** `WaveOperator<ParMesh>` with `--deriv-cache` at np>1; or an RK job
(`--time-integrator rk45`) with `--deriv-cache`.

**Suggested fix:** add an np=2 parallel equivalence test (mirror a
`tests/parallel/` harness) asserting cached==OnTheFly `ComputeVolumeRHS` per rank
to ≤1e-12; and either a `Mult`-path equivalence assertion or an explicit note
that `--deriv-cache` is validated for the ADER path only.

**Test case:**
```python
def test_R002_volume_rhs_parmesh_np2():
    # mpirun -np 2: per-rank max|cached - onthefly| <= 1e-12 for ComputeVolumeRHS
    run_mpi("seas_test_wave_operator_cached_parallel", np=2)
```

---

### [R-003] [LOW] [wave_operator.inl SetDerivMode] — the combined `2×` cache-budget guard is not independently tested; the R-004 byte test only checks the single-cache helper.

**Category:** QUALITY (test coverage)

**Description:** `SetDerivMode(Cached)` now checks
`bytes = 2 * ElementDerivativeCacheBytes(...)` (D_d^e + S_d^e). The existing
`R-004 cache bytes = 3*ne*ndof^2*8` assertion validates the SINGLE-cache helper
only, so the doubled budget is asserted only by inspection. (The abort path is
also not unit-tested — the harness has no death tests; noted as inherent.)
**Suggested fix:** assert the combined size, e.g.
`TEST_LE(|2*ElementDerivativeCacheBytes(ne,ndof) - expected_combined|, 0.0, …)`.

---

### [R-004] [LOW] [wave_operator.inl ComputeVolumeRHS Cached] — the cached branch lacks a homogeneous-order guard (`ndof == ndof_per_el_`); a heterogeneous-order mesh would use the wrong `dof_offset` silently.

**Category:** ASSUMPTION

**Description:** the cached branch uses `dof_offset = e * ndof_per_el_` with
`ndof = S.Height()` and no guard. The OnTheFly `ComputeVolumeRHS` *also* lacks
this guard (uses `fe->GetDof()` + the same offset), so this is **consistent**,
not a regression — but it is the same class-wide homogeneous-order assumption
flagged in Lever 1's deferred R-004. Contract-protected (the operator is
homogeneous-only). **Suggested fix:** add
`MFEM_VERIFY(ndof == ndof_per_el_, …)` for defense-in-depth symmetry with the
ApplySpatialDerivative Cached branch (if Lever-1 R-004 is fixed, fix here too).

---

### [R-005] [LOW] [bench_ader_hotpath.cpp] — `rhs_vol` accumulates unboundedly across reps AND across both DerivMode measures (never zeroed); harmless to timing but a latent footgun.

**Category:** QUALITY

**Description:** `rhs_vol` is declared once and `ComputeVolumeRHS_ForTest` does
`+=`, so over the auto-grown rep count (and both `measure(OnTheFly)` /
`measure(Cached)` calls) it grows without bound. No overflow at these scales and
timing is value-independent, so results are unaffected — but it is not zeroed per
measure and not included in the checksum. **Suggested fix:** zero `rhs_vol`
before each `measure` (outside the timed block) and add it to the final checksum,
or document that the accumulation is intentional and inert.

---

## Summary
- Critical issues: 0 (the cached volume kernel is numerically correct; default byte-exact)
- Moderate issues: 2 (R-001 heterogeneous-material coverage — the production path;
  R-002 parallel + Mult-path coverage)
- Low issues: 3 (R-003 combined-budget test, R-004 homogeneous guard, R-005 bench nit)
- Plan compliance: FULL — Lever 3 implements "precompute the volume operator …
  to a few mat-vecs"; PML loop correctly left uncached per the plan.
- Verdict: **PASS WITH FIXES** — the implementation is correct and default-safe;
  **R-001 (heterogeneous-material equivalence test) should land before trusting
  `--deriv-cache` on the heterogeneous TPV31 production run**. R-002 closes the
  parallel/RK gaps; R-003–R-005 are minor.

## Unreviewed Areas
- The Caliper-build perfgraph emission of the new path (Frontera-only).
- The `Mult` (RK) caller's `dQdt` pre-zeroing contract was assumed unchanged
  (pre-existing; the cached branch matches the OnTheFly `+=` semantics).
