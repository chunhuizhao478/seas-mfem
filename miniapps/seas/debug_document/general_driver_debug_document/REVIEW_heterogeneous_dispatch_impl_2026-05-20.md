# Code Review: Phase H Stage 2 implementation (TPV31 heterogeneous dispatch) — 2026-05-20

## Review Scope
- Plan: `debug_document/general_driver_debug_document/PLAN_heterogeneous_volume_bc_fault_dispatch_2026-05-20.md`
- Files reviewed (vs `git diff HEAD`):
  - `dynamic/godunov_flux_pool.cpp` (build-from-exact prerequisite)
  - `dynamic/wave_operator.hpp` (`FluxForElem_`, `ApplyJacobianPerElementDOF_` decls)
  - `dynamic/wave_operator.inl` (Groups A/B/C dispatch + guard relaxation)
  - `tests/unit/test_phaseh_wave_operator_constant_parity.cpp` (C-1..C-6)
  - `tests/unit/test_phaseh_heterogeneous_correctness.cpp` (H-1/H-2/H-3, new)
  - `tests/unit/test_phaser_dispatch_smoke.cpp` (death-test flip)
- Domain context: `miniapps/seas/CLAUDE.md` (Phase R / Phase H Stage 2 sections,
  byte-parity regression contract, fluctuation-Q), the prior plan review
  `REVIEW_heterogeneous_dispatch_plan_2026-05-20.md`.
- Verification performed: line-by-line read of every changed dispatch site;
  built and ran all three tests (serial): heterogeneous_correctness 10/10,
  constant_parity 30/30 (incl. C-6 fault Mult/AdvanceADER ~2.8e-12 / 1.2e-12),
  smoke 7/7.

## Correctness of the dispatch logic (audited, no bug found)
The following were checked explicitly and are correct:
- **A1 volume** (`ComputeVolumeRHS`): per-element `FluxForElem_(e).GetReferenceStarMatrix(0/1/2)`; bit-identical to scalar `Ax_/Ay_/Az_` on Mode::Constant (C-5 reports bit-diffs=0).
- **A2/A3 ADER CK** (`ApplyJacobianPerElementDOF_`): DOF blocking is exact because `ndof_total_ = ne_ * ndof_per_el_` (`wave_operator.inl:51`); per-element loop covers every DOF once; `D_next` is zeroed before the d-loop. Byte-identical on Constant (C-5/C-6).
- **B1/B2 BC**: `FluxForElem_(e1)` at all four RK4 + four ADER sites; `e1` is the one-sided owning element.
- **C1/C2/C3/C4 fault** (incl. the per-QP-batched second ADER-local block at `:5043`): `elem_plus/elem_minus = elem1_on_plus ? e1/e2 : e2/e1`, and the assembly (`:3506-3531`, `:4737-4760`, `:5067-5092`) maps `F_h_plus`→plus-element / `F_h_minus`→minus-element consistently. Shared faults correctly use `FluxForElem_(e1)` only (never off-rank `e2`).
- **Build-from-exact** (`godunov_flux_pool.cpp:106`): cached flux from exact `(λ,μ,ρ)`, dedup key still rounded — invariant #2 holds (C-5 `>6-sig` case passes).
- **Guard relaxation**: bi-material-across-fault is guarded by **always-on** `MFEM_VERIFY` (`fault_face_flux.cpp:332/472/743/868`), not debug-only `MFEM_ASSERT` — robust in release builds.

## Findings

### [R-001] MODERATE [tests/unit/test_phaseh_heterogeneous_correctness.cpp::H_2/H_3] — BC and fault per-element ROUTING is never exercised through Mult/AdvanceADER; the H-2/H-3 gates only check ingredients

**Category:** DEVIATION (from plan §6 acceptance intent: gates must "prove the fix does something real")

**Description:**
The implementation is correct, but the acceptance gates for Group B (BC) and
Group C (fault) do **not** actually exercise the per-element routing through
the operator:
- The Mode::Constant parity cases (`C_5_nonfault_parity` mixed BC, `C_6_fault_parity`)
  cannot detect a B/C routing regression: on Mode::Constant
  `FluxForElem_(e) == flux_` by construction, so a site that wrongly used
  `flux_` instead of `FluxForElem_(e1)` would still produce byte-identical
  output and PASS C-5/C-6.
- `H_2_boundary_per_element_impedance` and `H_3_fault_bulk_side_depth_varies`
  only check **ingredients**: that `GodunovFlux::AbsorbingTotal/FreeSurfaceTotal/
  Interior` differ between two materials (lines 283-286, 361-365) and that the
  per-element material array spans both layers (lines 314-316, 430-432). They
  never call `wave.Mult(...)` / `AdvanceADER(...)` on a heterogeneous operator,
  so they do not verify that the BC/fault dispatch *uses* `FluxForElem_(e1)`.

Contrast with **A1 (volume)**, which IS properly gated: `H_1_volume_layered_slab`
runs `ComputeVolumeRHS_ForTest` on the heterogeneous operator and asserts each
element's RHS equals the homogeneous operator of *that element's* material
(`:217-245`). A revert of the volume site to `flux_` would make H-1 fail. No
such operator-level check exists for B or C.

Concretely: for the heterogeneous (`Mode::Coefficient`) ctor, `flux_` is the
`(1,1,1)` placeholder (cp ≈ 1.73 m/s). If any BC or fault site were reverted to
`flux_`, heterogeneous physics would be catastrophically wrong, yet **every
existing test would still pass**. The plan (§6) lists H-2/H-3 as the gates that
"prove the fix does something real" for B and C; as implemented they prove the
*ingredients* are material-sensitive, not that the operator routes through them.

**Trigger:** A future edit reverts any Group-B or Group-C site from
`FluxForElem_(e1/elem_plus/elem_minus)` back to `flux_` (or a refactor drops the
per-element fetch). All current tests pass; heterogeneous BC/fault physics is
silently wrong.

**Actual behavior:** No automated gate fails on a B/C routing regression.

**Expected behavior:** An operator-level gate analogous to H-1 that runs Mult on
a heterogeneous mesh and asserts a boundary-/fault-adjacent element matches the
homogeneous operator built from *that element's* material (and would diverge if
the placeholder `flux_` were used).

**Suggested fix:** Strengthen H-2 and H-3 to be operator-level, mirroring H-1.
For H-2, build a two-layer slab with absorbing BC and pick a boundary-adjacent
element whose face-neighbours are all in the SAME layer (so its `Mult` RHS is
element-/layer-local), then assert the heterogeneous op equals the homogeneous
op of that layer at that element's DOFs:
```diff
+   // (c) operator-level: hetero Mult at a boundary element whose neighbours
+   //     are all in the same layer must equal the homogeneous op of that
+   //     layer — fails if the BC site uses the (1,1,1) placeholder flux_.
+   WaveOperator<Mesh> wave_homoA(mesh, k_order, k_lamA, k_muA, k_rhoA, bc);
+   WaveOperator<Mesh> wave_homoB(mesh, k_order, k_lamB, k_muB, k_rhoB, bc);
+   real_t Qbg[NUM_STATE] = {0};
+   wave.SetAbsorbingBackground(Qbg);
+   wave_homoA.SetAbsorbingBackground(Qbg);
+   wave_homoB.SetAbsorbingBackground(Qbg);
+   Vector Q(NUM_STATE*wave.GetScalarNDof()); FillQDeterministic(Q);
+   Vector dh(Q.Size()), da(Q.Size()), db(Q.Size());
+   wave.Mult(Q,dh); wave_homoA.Mult(Q,da); wave_homoB.Mult(Q,db);
+   // a deep-layer-A boundary element e_A (centroid z<z_split, all nbrs in A):
+   TEST_ASSERT(ElemMaxRelDiff(dh,da,e_A,npe,nt) < 1e-9 &&
+               ElemMaxRelDiff(dh,db,e_A,npe,nt) > 1e-3,
+               "hetero BC Mult at a layer-A boundary elem == homoA, != homoB");
```
For H-3, extend the `C_6_fault_parity` scaffold to a `Mode::Coefficient` fault
mesh where a chosen fault face and its two elements lie entirely in one layer,
seed the DOFData identically, run `Mult`/`AdvanceADER`, and assert that face's
DOFs equal the homogeneous op of that layer (and differ from the other layer).
This catches a fault-site revert to the placeholder `flux_`.

**Test case:**
```python
def test_R001_bc_routing_through_mult_is_per_element():
    # Mode::Coefficient two-layer slab, absorbing BC.
    # e_A = a z<z_split boundary element with all face-neighbours in layer A.
    # hetero.Mult(Q)[e_A] must equal homoA.Mult(Q)[e_A] (<1e-9 rel) and
    # differ from homoB.Mult(Q)[e_A] (>1e-3). A revert of the BC site to the
    # (1,1,1) placeholder flux_ makes hetero != homoA -> test fails.
    assert relerr(hetero_mult[e_A], homoA_mult[e_A]) < 1e-9
    assert relerr(hetero_mult[e_A], homoB_mult[e_A]) > 1e-3

def test_R001_fault_routing_through_mult_is_per_element():
    # Coefficient fault mesh; fault face F and its 2 elements all in layer A.
    # hetero.Mult at F's DOFs == homoA.Mult; reverting C1/C3 to flux_ fails it.
    assert relerr(hetero_mult[F_dofs], homoA_mult[F_dofs]) < 1e-9
```

---

### [R-002] LOW [tests/unit/test_phaseh_heterogeneous_correctness.cpp::H_2/H_3 messages] — TEST_ASSERT messages overclaim what is verified

**Category:** QUALITY (misleading test descriptions mask the R-001 gap)

**Description:**
The H-2/H-3 assertion messages assert routing/operator behavior that the test
does not actually check — e.g. `"per-element BC impedance is depth-correct"`
(`:315-316`) and `"shallow vs deep fault QPs feed FluxForElem_ different
Jacobians"` (`:431-432`). The test only verifies that the per-element material
array spans two layers, not that `FluxForElem_` is consulted by the dispatch.
A reader (or a future auditor) trusting these messages would believe the
routing is gated when it is not — which is how R-001 went unnoticed.

**Suggested fix:** Either implement R-001 (making the messages true) or soften
the messages to state what is actually checked:
```diff
-   TEST_ASSERT(saw_A && saw_B,
-               "absorbing-slab op has both layer-A and layer-B per-element "
-               "material (per-element BC impedance is depth-correct)");
+   TEST_ASSERT(saw_A && saw_B,
+               "absorbing-slab op stores both layer-A and layer-B per-element "
+               "material (INGREDIENT only; operator-level BC routing is gated "
+               "by test_R001_bc_routing_through_mult_is_per_element)");
```

**Test case:** N/A (documentation/message correctness; covered by R-001's test).

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001)
- Low issues: 1 (R-002)
- Plan compliance: **FULL on code** (Groups A/B/C, guard relaxation,
  build-from-exact, parity extension, H-1/H-2/H-3 all present and the dispatch
  logic is correct); **PARTIAL on test strength** (B/C operator-level routing
  not gated — R-001).
- Verdict: **PASS WITH FIXES** — no correctness bug in the shipped dispatch
  (verified by reading + all tests pass). R-001 is a regression-prevention gap
  that should be closed before this is treated as a hardened gate, given
  `wave_operator.*` is a "files requiring extreme care" target per CLAUDE.md.

## Unreviewed Areas
- **np=4 parallel parity** (`C_2_mult_parity_parallel`, shared-fault C2/C4):
  run serially only here (np=1). The macOS conda MPI stack has a known
  `Bus error` at MPI launch for these tests (per `spatial_dyn_heterogeneous_
  riemann_fix_2026-05-19.md`), so the np=4 gate must be confirmed on a cluster.
  The shared-fault sites (C2 `:4128`, C4 `:5707`) were read and are correct, but
  not executed multi-rank in this review.
- **End-to-end TPV31 dry-run / short run** through `seas_spatial_dyn_driver`
  (plan §Phase 4 acceptance): not run here (needs the 50 m or coarse mesh + the
  driver). Recommended as the final integration gate.
- BP5 (`seas_driver.cpp`): constructs no `WaveOperator`; unaffected (confirmed).
