# Implementation Plan: Non-fault interior face geometry/shape cache (FaceCache)

**Date:** 2026-06-24. **Author:** architect (Claude). **Status:** plan → review → implement.

## Overview

The optimized run (51463476) is 2.43× faster than baseline (830→342 s). The new
#1 cost is `WaveOperator::ComputeADERFaceFluxRHS` **compute** at **148.7 s = 43 %**
of runtime, and it is *balanced* (max/avg = 1.13), so reducing it cuts wall-clock
directly. Its per-QP, per-face inner loop recomputes **geometry-only** quantities
every macro-step — `GetFaceElementTransformations(f)`, `CalcOrtho` (face normal),
and `CalcShape` (both sides) — which are **constant across all 1445 steps** on the
affine straight-sided tet mesh. This is the exact anti-pattern `--deriv-cache`
(Lever 1/3) fixed for the volume kernels (ApplySpatialDerivative 222→17 s,
ComputeVolumeRHS 84→23 s). This plan adds the analogous cache for the **non-fault
interior** face path: precompute `{normal, weight, shape1, shape2}` once at setup
and reuse them, eliminating the per-step geometry recompute.

The production matrix path (`BimaterialWaveOperator`, `interior_flux="matrix"`)
**cannot** use the existing scalar `PrecomputedFaceFluxes` (it overrides
`UsePrecomputedFaceFluxes` to abort — it would leak the `(1,1,1)` `flux_`
sentinel), which is *why* it falls through to the per-QP path. This cache works
*with* the per-face bimaterial flux matrices: it caches only the geometry feeding
`InteriorFaceFlux_`, leaving the flux dispatch untouched.

## Constraints

- **Numerics-preserving, ≤ 1e-12 relative (NOT bit-exact).** Like the
  deriv-cache, precomputing the normal/shape re-associates round-off. The TPV*/BP5
  byte-exact regression contract is preserved by leaving the default OFF and the
  on-the-fly path byte-identical (see Phase 2).
- **Editable surface only:** `dynamic/` (`wave_operator.{hpp,inl}`, new
  `dynamic/face_geom_cache.hpp`, `bimaterial_wave_operator.*`) + the SAFS sbatch +
  a new test. **NO-TOUCH:** `bp5/ bp1/ bp2/ domain/ fault/ solver/
  friction/dieterich_ruina.hpp`.
- **Do NOT touch the fault-face branch** (`face_is_fault`, the avg-mode block at
  `wave_operator.inl:4400-4612`, the per-QP `T_can`/Riemann/friction). Caching its
  geometry interacts with the documented MPI-determinism hazard
  (OPT-FACEFLUX-FAULT-ROT-REUSE, rejected) and the σ_n fingerprint sensitivity.
  Fault faces are a tiny fraction of all faces; leave them on the per-QP path.
- **Do NOT touch the shared (cross-rank) face corrector**
  (`ComputeADERSharedFaceFluxRHS`) in this plan — that region is imbalance-bound,
  a *partitioning* problem, out of scope here.
- **Additive fast-path, existing slow-path unchanged.** The cache is a new branch
  that `continue`s; every non-cached face (fault, boundary, free-surface,
  absorbing, shared, scalar-operator) takes the existing path verbatim.
- **[REVIEW R-001, CRITICAL] Mutually exclusive with `UsePrecomputedFaceFluxes`.**
  The interior dispatch is an if/else ladder (`wave_operator.inl:4734` precomputed
  → `:4803` plain `InteriorFaceFlux_`). The cached helper replicates ONLY the
  `:4803` `InteriorFaceFlux_` branch, which is a DIFFERENT algorithm from
  `precomputed_face_fluxes_.AddInteriorFaceRhs`. The fast-path MUST be gated
  `use_face_cache_ && !use_precomputed_face_fluxes_`, and `SetUseFaceCache(true)`
  MUST `MFEM_VERIFY(!use_precomputed_face_fluxes_)` (and vice-versa in
  `UsePrecomputedFaceFluxes`), mirroring the R-1203 `SetMixedFluxMode`↔
  precomputed-flux cross-check. Latent today (SAFS never enables precomputed
  fluxes; the bimaterial class aborts on enable) but a hard correctness guard.
- **[REVIEW R-006, single-reader invariant]** `face_geom_cache_` is read in
  EXACTLY ONE call site (`ComputeADERFaceFluxRHS`). The RK path
  (`ComputeFaceFluxRHS`) and the shared corrector must NOT read it (central-flux
  faces, populated under RK `mf_on_`, would be mishandled — the cache stores no
  central/upwind distinction). Add a code comment stating this and a grep-based
  drift test asserting no second reader appears. (Under ADER, `mixed_flux` is
  forced `None` by the `ComputeMaxDt` guard, so `central_flux_face_set_` is empty
  anyway — belt-and-suspenders.)
- **Geometry-only cache shared by scalar + bimaterial** (mirrors
  `elem_derivative_cache.hpp`: "independent of NUM_STATE and of the
  scalar/bimaterial dispatch").
- **SAFS speedup + fingerprint are Expanse-only** (project memory: local =
  construction-order + small-fixture parity; fingerprint needs Frontera/Expanse).
  Local validation proves *correctness* (cached==onthefly ≤1e-12); the *speedup*
  and *production fingerprint* are confirmed on Expanse.

## Mesh/discretisation context (the math)

Per non-fault interior face `f` with adjacent elements `e1,e2`, the ADER
corrector adds the DG numerical-flux face term to the strong-form RHS:

    rhs_{e1} -= Σ_q w_q  φ^{e1}_i(x_q)  F^{+}(I_self(x_q), I_nbr(x_q), n(x_q))
    rhs_{e2} += Σ_q w_q  φ^{e2}_i(x_q)  F^{-}(...)

where `I` is the time-integrated state, `F^{±}` is the per-side numerical flux
(`InteriorFaceFlux_`), `n` is the **unit** outward normal of `e1`, `w_q =
ω_q·|J_face(x_q)|` is the physical quadrature weight, and `φ` are the L2 basis
functions. On an affine (straight-sided) tet mesh, `n(x_q)`, `|J_face(x_q)|`, and
`φ(x_q)` are **independent of the macro-step** — only `I_self/I_nbr` and the
resulting `F^{±}` change per step. The cache stores `{n_q, w_q, φ^{e1}(x_q),
φ^{e2}(x_q)}`; the per-step work reduces to the gather (`Σ_i φ_i I_i`), the flux
call, and the scatter — no `GetFaceElementTransformations`/`CalcOrtho`/`CalcShape`.

---

## Phase 1: FaceGeomCache infrastructure + isolated equivalence test

### Goal
A standalone, unit-tested builder produces per-(face,QP) cached geometry that
reproduces the on-the-fly `CalcOrtho`/`CalcShape`/weight values to ≤1e-12 — with
ZERO change to any runtime code path yet.

### Files to Create
- `dynamic/face_geom_cache.hpp` — the cache struct + free-function builder
  (header-only, mirrors `elem_derivative_cache.hpp`).
- `tests/unit/test_face_geom_cache.cpp` — isolated equivalence test.

### Detailed Requirements
1. Define the cache entry and container:
   ```cpp
   namespace mfem { namespace seas {
   struct FaceGeomEntry {
      int e1 = -1, e2 = -1;        // Elem1No, Elem2No (e2 >= 0: interior)
      int ndof1 = 0, ndof2 = 0;    // dofs per element (homogeneous => equal)
      int nqp = 0;
      std::vector<std::array<real_t,3>> nor;   // [q] UNIT normal of e1
      std::vector<real_t>               w;     // [q] ω_q·|J_face|
      std::vector<mfem::Vector>         shape1; // [q] φ^{e1}(x_q), size ndof1
      std::vector<mfem::Vector>         shape2; // [q] φ^{e2}(x_q), size ndof2
   };
   } }
   ```
2. Builder signature (free function, like `BuildElementDerivativeOperators`):
   ```cpp
   inline void BuildNonFaultInteriorFaceGeomCache(
      mfem::Mesh &mesh,
      const mfem::FiniteElementSpace &fes,
      const std::vector<int> &face_bdr_attr,        // WaveOperator::face_bdr_attr_
      int fault_attr,                                // bc_.fault_attr (0 => none)
      const std::set<int> &shared_mesh_face_set,     // skip cross-rank faces
      int quad_order,                                // 2*fe_order (NON-fault rule)
      std::unordered_map<int, FaceGeomEntry> &out);  // keyed by mesh face index
   ```
3. The builder MUST iterate faces with the **identical** predicate sequence as
   `ComputeADERFaceFluxRHS:3962-3990` so the cached set is exactly the faces the
   fast-path will own: skip `!ftr`; skip `e2<0 && shared_mesh_face_set.count(f)`;
   `is_boundary = (e2<0 && bdr_attr>0)`; skip `e2<0 && bdr_attr==0`;
   `face_is_fault = !is_boundary && bdr_attr==fault_attr && fault_attr>0`. **Only
   cache faces with `e2>=0 && !face_is_fault`** (interior, non-fault, non-shared).
   Boundary faces are EXCLUDED from this cache (Phase 1 scope).
4. For each cached face, use the **same** quadrature
   `IntRules.Get(ftr->GetGeometryType(), quad_order)` with `quad_order==2*order_`
   and the **same** `CalcOrtho`/normalise/`w=ip.weight*nor_len`/`Loc1.Transform`/
   `CalcShape` calls as the runtime loop (lines 3994-4008 and the shape2 site), so
   values match to round-off. Store the UNIT normal (post-normalisation) and
   `w = ip.weight * nor_len` (so the scatter `w*shape*F` is unchanged).
5. The builder is **geometry-only**: no `NUM_STATE`, no material, no flux. One
   cache serves `WaveOperator` and `BimaterialWaveOperator`.

### Edge Cases to Handle
- **Curved / non-affine mesh** (`mesh.GetNodes() != nullptr` with order ≥ 2): the
  cache is still *correct* (it stores the QP values), but the builder MUST NOT
  assume affine. It stores per-QP values regardless, so it is valid for curved
  meshes too — but document that the *speed benefit* assumes affine.
- **Serial vs ParMesh:** builder takes `Mesh&` (base) so it works for both;
  `shared_mesh_face_set` is empty in serial.
- **Empty cache** (no qualifying faces, e.g. single-element mesh): `out` is empty;
  the fast-path is simply never taken.
- **Heterogeneous FE order:** out of scope — `MFEM_VERIFY` all elements share
  `ndof_per_el` (the operator already assumes this; assert it in the builder).

### Acceptance Criteria
- [ ] `test_face_geom_cache`: on a 4×4×4 tet box, for every cached face and QP,
      the cached `{nor, w, shape1, shape2}` match a fresh on-the-fly
      `GetFaceElementTransformations`+`CalcOrtho`+`CalcShape` to ≤1e-12 relative.
- [ ] **[REVIEW R-004]** `nor` is compared **component-wise WITH SIGN** (not just
      magnitude): the SCALAR `WaveOperator::InteriorFaceFlux_` consumes `nor`
      directly (`wave_operator.inl:1247`), so a sign flip (invisible to the
      bimaterial per-face-matrix path) would corrupt the scalar path.
- [ ] The set of cached face indices equals the set the runtime predicate selects
      (assert by re-running the predicate in the test).
- [ ] **[REVIEW R-002]** On a 2-tet **fault** mesh (`fault_attr=3`, e.g.
      `test_bimaterial_wave_operator_parity.cpp::BuildTwoTetFaultMesh`), assert the
      fault face index is **absent** from the cache (fault-exclusion is exercised
      here, since the Phase-2 end-to-end fixtures may be fault-free).
- [ ] `make test` (existing unit suite) still passes (no runtime path touched).

### Dependencies
- Depends on: nothing. Required by: Phase 2.

---

## Phase 2: Wire the cached fast-path into ComputeADERFaceFluxRHS (flag-gated)

### Goal
With `--face-cache` on, `ComputeADERFaceFluxRHS` processes non-fault interior
faces via the cache (no per-step `GetFaceElementTransformations`/`CalcOrtho`/
`CalcShape`), producing results ≤1e-12 from the on-the-fly path; with the flag
off (default) the function is byte-identical to today.

### Files to Modify
- `dynamic/wave_operator.hpp` — add the cache member, the mode flag, the setter,
  and `#include "face_geom_cache.hpp"`.
- `dynamic/wave_operator.inl` — build the cache in the setter; add the fast-path
  branch at the top of the `ComputeADERFaceFluxRHS` face loop.
- `drivers/spatial_dyn_driver.cpp` — parse `--face-cache`, call the setter.
- `tests/parallel/test_bimaterial_deriv_cache_parity.cpp` — extend with a
  face-cache parity check.
- the SAFS caliper sbatch(es) — add `SAFS_FACE_CACHE` to `OPT_ARGS`.

### Detailed Requirements
1. `wave_operator.hpp` additions (private):
   ```cpp
   bool use_face_cache_ = false;                       // default OFF
   std::unordered_map<int, FaceGeomEntry> face_geom_cache_;
   ```
   Public setter (mirrors `SetDerivMode`):
   ```cpp
   /// Opt-in non-fault interior face geometry/shape cache (≤1e-12, NOT
   /// bit-exact; default OFF preserves the byte-exact on-the-fly path).
   /// Builds the cache immediately; idempotent.
   /// [REVIEW R-003] Lifecycle: the predicate reads bc_.fault_attr,
   /// face_bdr_attr_, and shared_mesh_face_set_ — ALL ctor-final (verified:
   /// SetMixedFluxMode does NOT mutate them).  So this may be called as EARLY as
   /// SetDerivMode (spatial_dyn_driver.cpp:1302); it does NOT have the late
   /// R-1408 lifecycle rule of SetMixedFluxMode.  Assert face_bdr_attr_ is
   /// populated.
   /// [REVIEW R-001] MFEM_VERIFY(!use_precomputed_face_fluxes_) — mutually
   /// exclusive with UsePrecomputedFaceFluxes (different interior algorithm).
   void SetUseFaceCache(bool enable);
   bool UsingFaceCache() const { return use_face_cache_; }
   ```
2. `SetUseFaceCache(true)` calls `BuildNonFaultInteriorFaceGeomCache(mesh_, *fes_,
   face_bdr_attr_, bc_.fault_attr, shared_mesh_face_set_, 2*order_,
   face_geom_cache_)` and sets `use_face_cache_=true`. `SetUseFaceCache(false)`
   clears the map and the flag.
3. In `ComputeADERFaceFluxRHS`, at the **top of the `for (f...)` loop body**
   (after `ftr` is fetched, before the existing per-QP work), insert:
   ```cpp
   if (use_face_cache_ && !use_precomputed_face_fluxes_) {   // [REVIEW R-001]
      auto it = face_geom_cache_.find(f);
      if (it != face_geom_cache_.end()) {
         ComputeADERFaceFluxRHS_CachedInterior_(it->second, I_data, rhs);
         continue;                  // skip the on-the-fly path for this face
      }
   }
   ```
   where the new private helper does ONLY what lines 4010-4823's non-fault
   interior `else` branch does, using cached geometry:
   ```cpp
   void ComputeADERFaceFluxRHS_CachedInterior_(
       const FaceGeomEntry &fc, const real_t *I_data, Vector &rhs) const;
   ```
   Body (byte-for-algorithm-identical to the slow path, cached geometry):
   ```
   const int o1 = fc.e1 * ndof_per_el_, o2 = fc.e2 * ndof_per_el_;
   for q in 0..fc.nqp:
      gather I_self[c]  = Σ_i fc.shape1[q](i) * I_data[c*ndof_total_ + o1 + i]
      gather I_nbr[c]   = Σ_i fc.shape2[q](i) * I_data[c*ndof_total_ + o2 + i]
      real_t F_h_e1[NUM_STATE], F_h_e2[NUM_STATE];
      InteriorFaceFlux_(facekey, I_self, I_nbr, fc.nor[q].data(), F_h_e1, F_h_e2);
      for c, i:  rhs[c*ndof_total_+o1+i] -= fc.w[q]*fc.shape1[q](i)*F_h_e1[c];
                 rhs[c*ndof_total_+o2+i] += fc.w[q]*fc.shape2[q](i)*F_h_e2[c];
   ```
   `facekey` = the mesh face index `f` (store it in `FaceGeomEntry` so the helper
   can pass it to `InteriorFaceFlux_`, which keys the per-face bimaterial matrix
   on it). **Add `int face_index` to `FaceGeomEntry`.**
4. The helper MUST replicate the EXACT loop order (`c` outer in gather, `i` inner;
   scatter `c`,`i`) of the slow path so cached==onthefly holds to ≤1e-12 and the
   `InteriorFaceFlux_` per-side semantics are preserved.
5. `--face-cache` parse in the driver mirrors `--deriv-cache`
   (`spatial_dyn_driver.cpp:1299-1308`); rank-0 log line; independent of the
   other flags.
6. sbatch **[REVIEW R-005]**: add `--face-cache` to the **upwind_ADER** caliper
   sbatch only — that is the run whose 148.7 s/43 % `ComputeADERFaceFluxRHS` this
   targets. On the **mixedflux_rk** caliper sbatches the flag is a measurement
   no-op (RK uses `ComputeFaceFluxRHS`, never reads the cache) AND wastes
   build/memory, so gate it: `if [[ "${SAFS_TIME_INTEGRATOR}" == "ader" &&
   "${SAFS_FACE_CACHE:-1}" == "1" ]]; then OPT_ARGS+=" --face-cache"; fi` — same
   integrator-aware pattern as `--shared-ck-recursion`.

### Interfaces
- New: `SetUseFaceCache(bool)`, `UsingFaceCache()`,
  `ComputeADERFaceFluxRHS_CachedInterior_(const FaceGeomEntry&, const real_t*,
  Vector&) const`.
- Unchanged: `InteriorFaceFlux_`, the on-the-fly loop, all fault/boundary/shared
  handling.

### Edge Cases to Handle
- **`SetAbsorbingBackground` not yet called / Q_bg semantics:** the cached
  interior helper does NOT touch `bulk_bg_`/absorbing (those are boundary/PML
  faces, excluded from the cache). No interaction with the R-1505 hoist.
- **Fault faces present (`bc_.fault_attr>0`):** they are NOT in the cache → take
  the unchanged on-the-fly fault path. Verify the cache excludes them.
- **Cache built before fault bookkeeping final:** guard with the same R-1408
  lifecycle note as `SetMixedFluxMode`; assert in the setter that
  `face_bdr_attr_` is populated (non-empty when `GetNumFaces()>0`).
- **PML enabled:** PML damping operates on volume/boundary, not interior non-fault
  faces — no interaction. (Confirm during implementation: the cached helper omits
  the PML branch, which only fires for boundary/absorbing faces.)
- **Mesh with no interior faces:** empty cache → fast-path never taken → identical
  to on-the-fly.

### Acceptance Criteria
- [ ] Extended `test_bimaterial_deriv_cache_parity` (np2 + np4): with
      `SetUseFaceCache(true)`, `AdvanceADER` (full predictor-corrector, multi-step)
      matches `SetUseFaceCache(false)` to ≤1e-12 relative, for a two-material
      fixture at FE order {1,2}, combined with both DerivModes.
- [ ] **[REVIEW R-004]** The SCALAR `WaveOperator` path is ALSO covered (it is the
      path that consumes `fc.nor`): add a scalar-operator face-cache parity check
      (the existing harness builds a bimaterial operator, whose `InteriorFaceFlux_`
      ignores `nor`). A serial scalar `WaveOperator` AdvanceADER cached==onthefly
      ≤1e-12 check suffices.
- [ ] **[REVIEW R-002, CRITICAL]** A FAULT-BEARING fixture (e.g. extend
      `test_bimaterial_wave_operator_parity` serial 2-tet fault mesh, or a new
      np2 seam-fault test): with a live fault, `SetUseFaceCache(true)` vs `(false)`
      `AdvanceADER` matches ≤1e-12 AND the fault face still takes the on-the-fly
      path (cross-check the result is unchanged from flag-off). The Phase-1
      builder test already asserts the fault face is absent from the cache.
- [ ] With the flag OFF, `ComputeADERFaceFluxRHS` output is **bit-identical** to
      pre-change (no cached path entered) — assert by an unchanged run of the
      existing `test_wave_operator_spatial_derivative` / parallel suites.
- [ ] `TestSharedCKParity` and the full `make test` + `make test-parallel` pass.
- [ ] (Expanse, follow-up) caliper `ComputeADERFaceFluxRHS` region drops; total
      wall-clock falls; SAFS fingerprint within tolerance before production use.

### Dependencies
- Depends on: Phase 1. Required by: (future) Phase 3 (boundary/free-surface faces;
  RK `ComputeFaceFluxRHS`; Expanse profiling).

---

## Testing Strategy
- **Phase 1:** isolated geometry equivalence (`test_face_geom_cache`) — pure
  builder vs on-the-fly, ≤1e-12.
- **Phase 2:** end-to-end parity (`test_bimaterial_deriv_cache_parity`, extended)
  — `AdvanceADER` cached-face==onthefly ≤1e-12, np2/np4, two-material, FE {1,2};
  plus a flag-OFF bit-identical guard. Reuse the existing harness
  (`RelDiff`, the two-material fixture, the FE-order sweep) — it already covers
  the deriv-cache; add a `SetUseFaceCache` axis.
- **Regression:** `make test` + `make test-parallel` (incl. the existing
  `test-bimaterial-deriv-cache-parity` and `TestSharedCKParity`).
- **Validation of the win:** Expanse caliper re-run (region delta on
  `ComputeADERFaceFluxRHS`) + SAFS fingerprint — explicitly out of local scope.

## Risk Assessment
- **Numerical drift > 1e-12:** if the builder's QP/normalisation order differs
  from the runtime loop, parity fails. Mitigation: copy the exact call sequence;
  the test catches it.
- **Cached set ≠ runtime-selected set:** a face cached but handled differently
  (or vice versa) double-counts or drops a flux. Mitigation: Phase-1 acceptance
  asserts set equality against the runtime predicate; the `continue` makes the
  fast-path and slow-path mutually exclusive per face.
- **Fault/boundary leakage:** accidentally caching a fault or boundary face would
  bypass the fault/BC physics. Mitigation: builder predicate excludes them;
  Phase-2 test uses a fault-bearing fixture variant to confirm fault faces still
  take the on-the-fly path (cross-check vs flag-off).
- **Hidden per-QP state in the slow path** (e.g. `SEAS_DIAG_FAULT_FLUX` diagnostic
  dumps at 4823+, the C-2B probe): the cached helper omits diagnostics — acceptable
  (diagnostics are a debug build option, not production), but document it.
- **Lifecycle (R-1408):** building the cache before fault/face bookkeeping is final
  yields a stale cache. Mitigation: setter asserts populated `face_bdr_attr_`;
  driver calls `SetUseFaceCache` after `SetAbsorbingBackground`/fault setup, same
  ordering as `--deriv-cache`.
- **Memory:** ~`n_faces · nqp · (3 + 1 + 2·ndof)` doubles/rank. For SAFS P1
  (~12k faces/rank, 3 QP, ndof=4): ~12k·3·12·8 B ≈ 3.5 MB/rank — negligible, well
  under the deriv-cache R-004 budget. No budget guard needed but log the size.
- **Uncertain payoff:** the win depends on the per-QP geometry fraction of the
  148 s. Strong prior: deriv-cache removed the identical anti-pattern from the
  volume kernels for a 13× / 3.6× reduction. If the Expanse re-profile shows
  < ~10 % wall reduction, the flag stays opt-in and Phase 3 is reconsidered.
