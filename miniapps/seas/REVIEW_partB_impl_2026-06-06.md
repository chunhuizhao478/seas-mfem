# Code Review: Part B implementation (bi-material FAULT per-side Riemann)

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md` (Part B, B0-B3)
- Files reviewed: `dynamic/fault_face_flux.{hpp,cpp}` (B2), `dynamic/wave_operator.hpp`
  (virtual SetFaultFlux + AssignFault no-op), `dynamic/bimaterial_wave_operator.{hpp,inl}`
  (B1 + SetFaultFlux override), `dynamic/heterogeneous_material.{hpp,cpp}` (B3 factory),
  `spatial/code/spatial_friction.{hpp,cpp}` (B3 config), `drivers/spatial_dyn_driver.cpp`
  (B1 call + B3 build), the two new tests, `Makefile`.
- Domain context: CLAUDE.md (sign conventions, no-local-full-mesh), the B0 verdict doc,
  the completion report.
- Verified CORRECT (not findings):
  - **B1 +/- convention matches the dispatch.** Mult block (`wave_operator.inl`) sets
    `Q_plus = elem1_on_plus ? Q_self(Elem1) : Q_nbr(Elem2)`; B1 assigns
    `Zp_plus from Elem1 when e1_plus` — consistent, no side swap.
  - **B0 math** (sigma_n_trial == welded star; R/T exact) — `seas_test_bimaterial_fault_riemann`.
  - **eps-offset = perpendicular** (J^-1 . n_inward): depth unchanged for depth-profile
    (byte-exact equal), step disambiguated for halfspace — `seas_test_bimaterial_fault_perside_material`.
  - **Lifecycle order:** driver calls AssignFault AFTER SetFaultDOFData; the
    `MFEM_VERIFY(ir.GetNPoints()==nbf)` guards mis-indexing.
  - **B2 guards** relax only with the flag; scalar path still aborts (exit 134).
  - **Sidecar material is Mode::Coefficient** (`MakeMaterialField`), so AssignFault's EvalAt
    works on the live SAFS matrix path (the GridFunction abort is NOT reachable — see R-002).

## Findings

### [R-001] MODERATE — bimaterial_wave_operator/SetFaultFlux: the operator→flag wiring (the whole B2 integration) is untested

**Category:** EDGE_CASE (test completeness for a load-bearing integration)

**Description:**
B2 relies on `BimaterialWaveOperator::SetFaultFlux` calling `ff->SetPerSideFluxApplied(true)`
(via the newly-`virtual` base method) so the matrix bi-material fault does NOT abort.  No unit
test exercises this: `seas_test_bimaterial_fault_riemann` sets the flag MANUALLY
(`ff.SetPerSideFluxApplied(true)`), not through the operator.  If a refactor dropped the
`virtual` keyword, or the override, every matrix bi-material-fault run would abort at runtime
and the unit suite would stay green.

**Trigger:** remove `virtual` from `WaveOperator::SetFaultFlux`, or remove the bimaterial
override — all current tests still pass.

**Actual behavior:** the operator→flag path is only validated by the driver building, not by a test.

**Expected behavior:** a test asserts `BimaterialWaveOperator::SetFaultFlux(&ff)` leaves
`ff.GetPerSideFluxApplied() == true`, and a scalar `WaveOperator::SetFaultFlux(&ff)` leaves it
`false`.

**Suggested fix (append a block to `tests/unit/test_bimaterial_fault_perside_material.cpp`,
which already builds both operator types):**
```cpp
   // B2 wiring: the matrix operator's SetFaultFlux affirms per-side-A; scalar does not.
   {
      FunctionCoefficient lc([](const Vector&){return 3.0e10;});
      FunctionCoefficient mc([](const Vector&){return 3.0e10;});
      ConstantCoefficient rc(2670.0);
      MaterialField mat = MaterialField::MakeCoefficient(&lc, &mc, &rc);
      BimaterialWaveOperator<Mesh> opb(mesh, order, mat, bc);
      FaultFaceFlux ffb(2670.0, 6000.0, 3464.0);
      opb.SetFaultFlux(&ffb);
      TEST_ASSERT(ffb.GetPerSideFluxApplied(),
                  "BimaterialWaveOperator::SetFaultFlux affirms per-side-A (B2 wiring)");
      WaveOperator<Mesh> ops(mesh, order, 3.0e10, 3.0e10, 2670.0, bc);
      FaultFaceFlux ffs(2670.0, 6000.0, 3464.0);
      ops.SetFaultFlux(&ffs);
      TEST_ASSERT(!ffs.GetPerSideFluxApplied(),
                  "scalar WaveOperator::SetFaultFlux leaves the flag false (guard intact)");
   }
```
(NOTE: confirm the scalar `WaveOperator<Mesh>` ctor signature `(mesh, order, lambda, mu, rho, bc)`
against `wave_operator.hpp` before using; adjust if it differs.)

**Test case:** the snippet above IS the test.

---

### [R-002] LOW — bimaterial_wave_operator.inl: AssignFaultSidePerMaterialImpedances aborts on a GridFunction-mode material

**Category:** ASSUMPTION (defensive robustness; not currently reachable)

**Description:**
`AssignFault` calls `material_->EvalAt(...)`, which `MFEM_ABORT`s in `Mode::GridFunction`.
No current driver path feeds a GridFunction material to the operator (`MakeGridFunction` has
zero callers; depth-profile/halfspace/sidecar are all `Mode::Coefficient`, constant is
`Mode::Constant`), so the abort is unreachable today.  But the driver calls AssignFault
UNCONDITIONALLY on the matrix path, and `material_` is a general `MaterialField*`; a future
GridFunction material would crash here instead of degrading gracefully.

**Trigger:** a `Mode::GridFunction` material on the matrix path (no such config today).

**Actual behavior:** would MFEM_ABORT inside the per-side assignment.

**Expected behavior:** skip the per-side assignment for GridFunction mode (leave the
driver-seeded values; a nodal GF material has no closed-form per-side eps-offset and is
typically continuous/fault-symmetric anyway).

**Suggested fix (in `AssignFaultSidePerMaterialImpedances`, after the null check):**
```diff
   if (material_ == nullptr) { return; }
+  // GridFunction-mode material has no QP-evaluable per-side eps-offset (EvalAt aborts);
+  // skip and keep the driver-seeded values.  Constant/Coefficient are handled below.
+  if (material_->mode == MaterialField::Mode::GridFunction) { return; }
```

**Test case:** not added (no GridFunction material fixture exists; LOW + unreachable, rule 4).
The guard is a no-op for the tested Coefficient/Constant paths.

---

### [R-003] LOW — fault_face_flux.cpp: the 5 bimaterial-guard abort messages are now misleading

**Category:** QUALITY

**Description:**
The guards now read `per_side_flux_applied_ || (homog_ok...)`, but the message still says
"Extend GodunovFlux to per-side A before running this configuration."  The per-side path now
EXISTS; the guard fires only when the caller did not affirm it (the scalar path).  The message
should point at the real cause.

**Trigger:** a bi-material fault on the scalar (non-matrix) path.

**Suggested fix (the message is identical in all 5 guards; update via replace_all):**
```diff
-               "before running this configuration.");
+               "before running this configuration, OR run on the matrix interior-flux "
+               "path (BimaterialWaveOperator), which affirms per-side A via "
+               "SetPerSideFluxApplied and is the supported bi-material-fault path.");
```
(Match the exact surrounding text of the shared message string.)

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — untested operator→flag wiring)
- Low issues: 2 (R-002 GridFunction-abort guard [unreachable]; R-003 stale guard message)
- Plan compliance: FULL for B0-B3; the shared-fault peer exchange (R-101) is a
  plan-sanctioned documented deferral.
- Verdict: PASS WITH FIXES — no critical bugs.  B1's +/- convention, the eps-offset, the B0
  math, and the guard relaxation are correct and tested; apply R-001 (closes the load-bearing
  integration gap), R-002 and R-003 (robustness/clarity).

## Unreviewed Areas
- Parallel (np>1) shared-fault per-side material — deferred (R-101); not runnable locally.
- TPV6 production behavior (the actual bi-material rupture) — Frontera, user-submitted; needs a
  Part C config (out of Part B scope).
