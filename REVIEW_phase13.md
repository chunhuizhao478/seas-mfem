# Code Review: Phase 13 — Split scalar / matrix into separate WaveOperator classes (Round 2, 2026-05-29, fresh adversarial pass + full local test run)

## Review Scope
- Plan: `document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md` §Phase 13 (lines 2725–3030)
- Files reviewed (this round, re-read from scratch):
  - `dynamic/wave_operator.hpp`, `dynamic/wave_operator.inl` — the four virtual hooks, all dispatch/deposit sites, the scalar default bodies, the 8 R-001 fault imposed-state sites, `ComputeMaxDt`, `SetMixedFluxMode`, ctor.
  - `dynamic/bimaterial_wave_operator.hpp`, `dynamic/bimaterial_wave_operator.inl` — the new subclass (ctor delegation + `(1,1,1)` poison, the four overrides, the moved builders).
  - `drivers/spatial_dyn_driver.cpp` — the scalar/matrix `wave_ptr` branch (1062–1101), the R-107 reflection/PML guard, the `cp_seed`/`FaultFaceFlux` seed (1544).
  - `tests/unit/test_bimaterial_wave_operator_parity.cpp` (C-6, new), `tests/unit/test_phaseh_wave_operator_constant_parity.cpp` (C-5, repointed), `tests/unit/test_wave_operator.cpp` (repointed).
  - `Makefile` — build targets + the `make test` run chain.
- Domain context: `CLAUDE.md` (byte-exact contract; "no local full-mesh runs"), memory notes (GodunovFluxPool dedup floor; placeholder-leak class; Phase-13 subclass split), git `HEAD~1` (the pre-split single-class form, used as the faithfulness oracle).

## Method (what is new in this round)
1. **Re-executed all three review passes from scratch** on the changed files (skill rule #10).
2. **Verified the "moved verbatim" claim mechanically:** `diff -w` of the bimaterial bodies extracted from `HEAD~1:wave_operator.inl` against `bimaterial_wave_operator.inl` (class-qualifier normalised) → **identical** (only the added `template<>` lines differ). The old inline interior-face dispatch + deposit (`HEAD~1` ComputeFaceFluxRHS 3283–3360) matches the new scalar-default + bimaterial-override split byte-for-byte, including the per-side deposit (`-F_h_e1`→Elem1, `+F_h_e2`→Elem2) and the bimaterial flux args (`Q_self=Q_e1, Q_nbr=Q_e2`, not swapped).
3. **Built and RAN the full Phase-13 acceptance suite locally** (`mfem-dev`, MFEM lib from the main checkout; sources from this worktree). All unit-level acceptance criteria pass — see "Test evidence" below. (The byte-exact full-mesh TPV/BP5 gate remains out of local scope per CLAUDE.md.)

## Test evidence (run this round)
| Test | Result |
|---|---|
| `seas_test_bimaterial_wave_operator_parity` (C-6, fault-bearing) | **9/9 PASS** (Mult, AdvanceADER, ComputeMaxDt; finite, field-scale rel < 1e-9) |
| `seas_test_phaseh_wave_operator_constant_parity` (C-5) np=1 | **19/19 PASS** |
| `seas_test_phaseh_wave_operator_constant_parity` (C-5) np=4 | **52/52 PASS** (exercises the shared-face bimaterial interior path) |
| `seas_test_wave_operator` | **20/20 PASS** (incl. R-003 matrix×mixed-flux abort) |
| `seas_test_godunov_flux_bimaterial` | **21/21 PASS** |

## What was verified correct (the fix agent should NOT touch this)
- **Scalar byte-exactness is construction-guaranteed.** `ComputeVolumeRHS` reads `FluxForElem_(e).GetReferenceStarMatrix(d)`; on the scalar class `FluxForElem_(e)==flux_` and `GetReferenceStarMatrix(d)` returns `ref_star_[d]` which is built by the *same* `GodunovFlux::BuildJacobian(d,·)` that fills the cached `Ax_/Ay_/Az_` — bit-identical. The scalar hook bodies (`InteriorFaceFlux_`/`SharedInteriorFaceFlux_`/`ApplyElementJacobian_`) and the scalar `ComputeMaxDt` are verbatim the pre-split code (mixed-flux per-face choice via `mf_on_ && central_flux_face_set_.count(face)` preserved).
- **`(1,1,1)` poison is provably dead on the matrix path.** The inherited `flux_` is read in the base ctor *only* via `flux_(λ,μ,ρ)` + `flux_.BuildJacobian(0/1/2, Ax_/Ay_/Az_)`; those `Ax_/Ay_/Az_` are then dead (every consumer routes through `FluxForElem_`), and `h_min_` is geometric (from `mesh_`). No virtual hook is invoked during base construction. Confirmed by full grep: the only bare-`flux_` material uses left in `wave_operator.inl` are the overridden scalar hook bodies + the dead ctor `BuildJacobian` + the overridden scalar `ComputeMaxDt`.
- **All 8 R-001 fault imposed-state sites** route through `FluxForElem_(elem_plus/elem_minus/e1/qa.local_elem)` with the correct per-side selection (`elem1_on_plus ? e1 : e2`) and per-side deposit. Faithful to `HEAD~1`.
- **Driver matrix branch** constructs `std::make_unique<BimaterialWaveOperator<ParMesh>>(...)` held as `unique_ptr<WaveOperator<ParMesh>>`; base dtor is virtual (via `TimeDependentOperator`), so the subclass dtor frees `owned_flux_pool_` — no leak.
- **No compile stragglers:** the removed het ctor (`WaveOperator(.,.,MaterialField,.)`) and the moved accessors (`GetPerElementMaterial`, …) are referenced only on `BimaterialWaveOperator` objects (C-5/test_wave_operator, repointed). The scalar class is `MaterialField`-free.

## Findings

### [R-001] [MODERATE] [Makefile / `test:` target] — C-6 not run by `make test` — **RESOLVED (verified this round)**

**Status:** FIXED in the reviewed commit. `test-bimaterial-wave-operator-parity` is now in the `make test` prerequisite chain at **Makefile:3729** (immediately after `test-phaseh-wave-operator-constant-parity:3728`), and its run rule is at Makefile:4028. No action.

---

### [R-002] [MODERATE] [spatial_dyn_driver.cpp reflection/PML cp] — matrix path computed `cp = sqrt(0/0) = NaN` — **RESOLVED (verified this round)**

**Status:** FIXED. The R-107 reflection warning + PML `pml_cp` are now guarded by `if (material.mode == MaterialField::Mode::Constant)` at **spatial_dyn_driver.cpp:1132** (comment cites "Phase 13 (REVIEW R-002)"). The matrix path no longer emits the bogus `cp_max (0 s)` warning. No action. (See R-007 for a remaining sibling of this same NaN class that is harmless but unguarded.)

---

### [R-003] [MODERATE→RESOLVED] [test C-6 metric] — field-scale vs per-component "to 1e-9" — **RECONCILED IN PLAN**

**Status:** RECONCILED. The plan's C-6 acceptance line (PLAN §Phase 13, lines 2988–2999) now carries the field-scale-metric rationale (round_sig(·,6) dedup floor) as the documented acceptance, and the test header documents it. The deviation is recorded in both the plan and the test, which is exactly what the prior finding required. The metric was tripwire-validated and the test passes at 3e-10 (Mult). No action; do NOT revert to per-component (it would false-fail at ~1e10-scale moduli).

---

### [R-004] [LOW] [test_bimaterial_wave_operator_parity.cpp] — homogeneous C-6 cannot detect a per-side / true-bimaterial fault-flux error

**Category:** EDGE_CASE — **STILL OPEN (accepted)**

**Description:**
C-6 (and C-5) use a **homogeneous** material, so `FluxForElem_(elem_plus) == FluxForElem_(elem_minus)` and the per-side selection / per-side bi-material flux asymmetry (`side0` vs `side1`, `A_e1` vs `A_e2`) is never exercised. A genuine bimaterial fault error (or an `e1`/`e2` swap) would pass. Per-side correctness rests on faithfulness to `HEAD~1`/hrs-ref (verified by diff this round), not on a test. The C-6 header (lines 28–38) already documents this scope. Note: the underlying bimaterial flux composition was ported in Phase 9 and is out of Phase-13 scope; Phase 13 only relocated it (move verified byte-faithful).

**Suggested fix (follow-up, non-blocking):** a true-bimaterial regression fixture (different `(λ,μ,ρ)` across the fault, with a hand-checked or hrs-ref-anchored expected flux) would close the gap. Not required for Phase-13 sign-off.

**Test case:** N/A (coverage note; LOW).

---

### [R-005] [LOW] [coverage] — substep fault-flux sites (per-QP `elem_plus_qq` at `wave_operator.inl:4266`; shared `qa.local_elem` at `:5256`) untested on the matrix path

**Category:** EDGE_CASE — **STILL OPEN (accepted)**

**Description:**
C-6 exercises `Mult` (interior-fault site, :2712) and `AdvanceADER` inline-fault (:3930) on a serial mesh. The per-substep dispatch sites — `ComputeADERFaceFluxRHS` per-QP averaged branch (:4266, reached only when `substep_I_imp_*` is set) and `ComputeADERSharedFaceFluxRHS` (:5256, parallel-only) — are not exercised on a bimaterial operator by any unit test. Routing was verified correct by inspection (per-side `FluxForElem_` selection identical to the tested sites) and by faithfulness to `HEAD~1`. The C-6 header documents this.

**Suggested fix (follow-up, non-blocking):** an np>1 bimaterial substep regression case.

**Test case:** N/A (coverage note; LOW).

---

### [R-006] [LOW] [bimaterial_wave_operator.{hpp,inl} / wave_operator UsePrecomputedFaceFluxes] — `UsePrecomputedFaceFluxes` is a material-dependent site NOT closed by the Phase-13 "by construction" invariant

**Category:** ASSUMPTION / BUG (latent) — **NEW this round**

**Description:**
`WaveOperator::UsePrecomputedFaceFluxes(true)` passes the inherited `flux_` to `PrecomputedFaceFluxes::Init(... flux_ ...)` and `InitSharedFaces(... flux_ ...)` (`wave_operator.inl:683, :702`). It is **not** `virtual` and is **not** overridden on `BimaterialWaveOperator`. On a bimaterial object `flux_` is the `(1,1,1)` poison, so enabling precomputed face fluxes would build interior-face flux matrices from `(1,1,1)` — a silent placeholder leak of exactly the class Phase 13 claims to make "impossible by construction." It also bypasses `InteriorFaceFlux_` at runtime (the precomputed branch at `:2833`/`:3389` deposits directly), so the C-6 tripwire would NOT catch it. Precomputed face fluxes is mathematically a scalar-only (single-material) optimisation and is conceptually incompatible with the matrix path, exactly like mixed flux (which IS structurally blocked via the `SetMixedFluxMode` override, R-003).

**Severity rationale (LOW, not higher):** currently **unreachable** — the spatial driver never calls `UsePrecomputedFaceFluxes`, and every test caller (`test_arm1_*`, `test_arm2_*`, `test_adjacent_triangle_*`) uses a scalar `WaveOperator`. So there is no live wrong-physics path today. It is a defense-in-depth gap that contradicts the plan's stated invariant.

**Trigger:** `BimaterialWaveOperator<…> w(...); w.UsePrecomputedFaceFluxes(true);` then any `Mult`/`AdvanceADER`.

**Actual behavior:** silently builds/uses precomputed interior-face flux matrices from the `(1,1,1)` poison `flux_`; the per-face bimaterial matrices are ignored on the precomputed path.

**Expected behavior:** abort loudly (mirror `SetMixedFluxMode`), since the matrix path replaces the interior-face flux with the bi-material Riemann solve.

**Suggested fix:** make the base method virtual and override on the subclass to abort.
```diff
// dynamic/wave_operator.hpp
-   void UsePrecomputedFaceFluxes(bool enable);
+   virtual void UsePrecomputedFaceFluxes(bool enable);
```
```diff
// dynamic/bimaterial_wave_operator.hpp  (public section, next to SetMixedFluxMode)
+   void UsePrecomputedFaceFluxes(bool enable) override
+   {
+      MFEM_VERIFY(!enable,
+                  "BimaterialWaveOperator::UsePrecomputedFaceFluxes: "
+                  "precomputed face fluxes are scalar-only; ...");
+   }
```

**Test case:**
```cpp
// tests/unit/test_wave_operator.cpp (next to the R-003 sibling)
void TestR006MatrixPrecomputedFluxAborts()
{
   BimaterialWaveOperator wave_het(*mesh, order,
      MaterialField::MakeConstant(lambda,mu,rho), bc);
   // EXPECT: wave_het.UsePrecomputedFaceFluxes(true) aborts (RunAbortsInChild);
   //         wave_het.UsePrecomputedFaceFluxes(false) does NOT abort.
}
```

---

### [R-007] [LOW] [spatial_dyn_driver.cpp:1544 cp_seed/FaultFaceFlux] — matrix path builds `FaultFaceFlux(0, NaN, NaN)` (harmless, but the same NaN class R-002 was fixed for)

**Category:** QUALITY / robustness — **NEW this round**

**Description:**
On the matrix path `material` is `Mode::Coefficient`, so `material.{lambda,mu,rho}_const == 0` (`heterogeneous_material.hpp:77-79`, `MakeCoefficient` leaves them at the struct default). The seed impedances at `spatial_dyn_driver.cpp:1544-1548` therefore compute `cp_seed = sqrt((0+0)/0) = NaN`, `cs_seed = NaN`, and construct `FaultFaceFlux fault_flux(0.0, NaN, NaN)`. This is the **same NaN-from-zero-const class** that R-002 fixed for the reflection/PML cp — but here it was not guarded.

**Why it is LOW, not a live bug:** `FaultFaceFlux`'s ctor does not assert (`fault_face_flux.cpp:36-41`), so there is no abort; and its scalar members `rho_/cp_/cs_/Zp_/Zs_` are **write-only** — set in the ctor, declared `private` with no getter, and **never read** anywhere (the physics consumes the per-DOF `DOFData.{Zp_plus,eta_p,…}` set by `InitializeFaultDOFs_Spatial`). Verified by grep: the only references to `Zp_/Zs_/cp_/cs_/rho_` are the ctor assignments and the private declarations. So the NaN dead-ends. The seed is genuinely "overwritten by InitializeFaultDOFs_Spatial" as the comment claims.

**Trigger:** any `interior_flux="matrix"` run with a fault (e.g. TPV31).

**Actual behavior:** `FaultFaceFlux` holds NaN seed impedances (never consumed).

**Expected behavior:** for cleanliness / to avoid a future trap if any code starts reading the scalar seed, derive the seed from a representative element or guard on Constant.

**Suggested fix (defensive, non-blocking):** derive `lam_seed/mu_seed/rho_seed` from element 0 via `material.EvalAt(0, *T0, ip0, …)` on the non-Constant path (same pattern `BuildGodunovFluxPool_` uses); keep `material.*_const` on the Constant path (byte-identical).

**Test case:** pseudocode — `MakeCoefficient` leaves `*_const == 0`, so the naive `sqrt((0+0)/0)` is NaN; after the fix the seed is finite (`EvalAt(elem0)` → `cp_seed > 0 && isfinite`).

---

## Summary
- Critical issues: **0**
- Moderate issues: **0 open** (R-001, R-002 verified FIXED this round; R-003 reconciled in the plan)
- Low issues: **4** — R-004 (homogeneous-only fault coverage, accepted), R-005 (substep matrix sites untested, accepted), R-006 (NEW: `UsePrecomputedFaceFluxes` defense-in-depth gap, unreachable today), R-007 (NEW: matrix-path `FaultFaceFlux(0,NaN,NaN)` seed, harmless/dead members)
- Plan compliance: **FULL** — every Phase-13 requirement is implemented and faithful: scalar/matrix class split, four virtual hooks with verbatim scalar defaults, 8 R-001 fault sites routed per-element, `(1,1,1)` poison contained, mixed-flux structurally blocked on the matrix path, driver dispatch, C-6 added + repointed C-5/test_wave_operator, Makefile wired (build + `make test`). The "moved verbatim" claim was mechanically verified (diff vs `HEAD~1`). Scalar byte-exactness is construction-guaranteed but NOT full-mesh-run-verified locally (forbidden per CLAUDE.md) — confirm on Frontera against a pre-Phase-13 checkpoint before production sign-off.
- Verdict: **PASS** — the core split is correct, faithful, and fully green on the local acceptance suite (C-5 19/19 + 52/52 np=4, C-6 9/9, test_wave_operator 20/20, godunov_flux_bimaterial 21/21). The two prior MODERATE findings are fixed. R-006 and R-007 are optional LOW defensive hardening (neither is a live wrong-physics path today); R-004/R-005 are accepted coverage gaps for a future true-bimaterial / np>1 substep regression test.

## Unreviewed Areas
- **Scalar TPV205/102/104 + BP5 byte-identical station/benchmark traces** — requires full-mesh runs (forbidden locally per CLAUDE.md / memory). Byte-exactness is construction-guaranteed (FluxForElem_==flux_, GetReferenceStarMatrix==Ax_, verbatim hook bodies) but not run-verified; confirm on Frontera before production sign-off.
- **True heterogeneous (non-homogeneous-Coefficient) matrix correctness** (TPV31 / depth profile) — the bimaterial flux composition is Phase-9 physics (moved-not-changed by Phase 13; move verified byte-faithful); its end-to-end correctness is a Phase-10 deliverable, out of Phase-13 scope.
- **np>1 matrix path** — `ExchangeBiMaterialNeighbours_` is a local-side stub (moved-not-fixed per the plan; warns for laterally-heterogeneous Coefficient + shared faces). Lateral-heterogeneous parallel runs are explicitly deferred; TPV31 is depth-only/seam-continuous so the stub is correct there.

---

## Round-2 fix status (applied by /code-fix, 2026-05-29)
- **R-006 — FIXED** in this worktree: `WaveOperator::UsePrecomputedFaceFluxes` made `virtual` (`wave_operator.hpp`); `BimaterialWaveOperator` overrides it to `MFEM_VERIFY(!enable, …)` (`bimaterial_wave_operator.hpp`); new `TestR006MatrixPrecomputedFluxAborts` in `test_wave_operator.cpp` (matrix enable aborts; matrix disable is a no-op; scalar control omitted because the base precomputed-flux path is tet-only and the test mesh is hex). `seas_test_wave_operator` 23/23.
- **R-007 — FIXED** in this worktree: matrix-path seed now derived from element 0 via `material.EvalAt` (finite); scalar path keeps `*_const` byte-identically (`spatial_dyn_driver.cpp`). Driver object recompiles clean.
- **R-001/R-002/R-003** — already resolved/reconciled in the reviewed commit (no action).
- **R-004/R-005** — accepted coverage gaps; already documented in the C-6 test header (no code change).
- Regression after fixes: C-6 9/9, C-5 19/19 (np=1) + 52/52 (np=4), test_wave_operator 23/23, all green.
