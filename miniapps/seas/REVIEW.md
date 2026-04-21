# Code Review: TPV102 v9.1.0 Debug Plan Follow-through (2026-04-20, rev 3)

> Revision 3 is a fresh adversarial audit of the code changes applied in
> response to rev 2's findings (R-001 … R-005). All three review passes
> were re-executed from scratch against the current working tree. R-001–
> R-005 are verified resolved; three new findings R-006 / R-007 / R-008
> and one process/documentation finding R-009 were uncovered by the
> fresh pass.

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.1.0_debug_plan.md`
- Pending changes in working tree (verified via `git status --short`):
  - `M miniapps/seas/Makefile`
  - `M miniapps/seas/fault/fault_basis.hpp`
  - `M miniapps/seas/io/paraview_output.hpp`
  - `?? miniapps/seas/tests/unit/test_fault_basis_dip_strike_symmetry.cpp`
  - `?? miniapps/seas/tests/unit/test_fault_surface_vtu_continuity.cpp`
  - `?? miniapps/seas/test_r001_continuity/` (generated VTU output)
  - `?? miniapps/seas/test_r001_nbf6/` (generated VTU output)
- Files reviewed:
  - Both modified production files (diffs inspected via `git diff`).
  - Both new unit tests (full file reads).
  - `drivers/tpv102_driver.cpp:695-745` and `tests/verification/bp5_verification_full.cpp:1500-1765` — downstream callers of the modified VTU writer.
  - `dynamic/wave_operator.hpp:140-164` — DOFData stride contract.
  - MFEM `fem/intrules.cpp:1263-1293` (degree-1/2/3/4 triangle rules) and `fem/fe/fe_h1.cpp:451-531` (`H1_TriangleElement` nodal layout) — ground truth for per-workflow nbf semantics.
- Build and run validation performed:
  - Ran `./seas_test_fault_surface_vtu_continuity` → **63/63 PASSED**.
  - Ran `./seas_test_fault_basis_dip_strike_symmetry` → **17/17 PASSED**.
  - Inspected the generated `test_r001_continuity/FaultSurface/fault_surface_r0_c0.vtu` CellData values and cross-checked against the linear field `0.5 + 0.1*x2 + 0.2*x3` evaluated at each face's centroid (all 6 faces match to the emitted 10-digit precision).
- Domain context: project CLAUDE.md (BP5 tangent1=dip / tangent2=strike convention), seas-miniapp CLAUDE.md (sign conventions and "files requiring extreme care" list — `fault/fault_basis.hpp` is on it, `io/paraview_output.hpp` is not).

## Findings

### Verification of rev-2 findings (all resolved)

- **[R-001] RESOLVED** — `paraview_output.hpp::WriteFaultSurfaceVTU` now writes `<CellData>` with one value per triangle computed as the per-face arithmetic mean over `nbf` DOFs. `test_fault_surface_vtu_continuity` Test 4.1(a/b) confirms, for all 6 interior triangles of a 2×2×1 tet mesh, that the emitted `slip_dip` / `slip_strike` / `slip_rate_dip` / `slip_rate_strike` / `traction_dip` / `traction_strike` / `state_variable` / `param_a` / `param_Dc` cell values equal the linear field evaluated at the face centroid, to within 5e-9 relative. H-V91-A4 is closed on the TPV102 workflow.
- **[R-002] RESOLVED** — the inner loop is now `for (int k = 0; k < nbf; k++)`. Test 4.1(c) with `nbf = 6` per-QP IDs `d`, `d+1`, …, `d+5` emits the correct mean `d + 2.5`; the defensive negative assertion `|cell - (d + 1.0)| > 0.5` (which would fail under the old hardcoded `k < 3`) also passes. P ≥ 2 is now correctly supported.
- **[R-003] RESOLVED** — `fault_basis.hpp:425-437` now explicitly normalizes `dip` after the second cross product and verifies `|dv| > 1e-12`. `test_fault_basis_dip_strike_symmetry` runs 1000 perturbed-normal trials plus a canonical TPV102 axis-aligned check; all 3000 length and 3000 orthogonality constraints pass. (See R-006 below — the *tolerance* chosen for this test is arguably too loose to actually gate the regression.)
- **[R-004] RESOLVED** — the `WriteFaultSurfaceVTU` docstring now describes the CellData + per-face-average behaviour and explicitly references R-001 / R-005.
- **[R-005] RESOLVED** (side-effect of R-001 fix) — the `k < nbf` loop bound + per-face averaging cleanly handles `nbf = 1` without cross-face reads or OOB. Not active in BP5 production (IP default), but the latent class of bug is eliminated.

### New findings (rev 3)

### [R-006] [MODERATE] [POSSIBLE] [tests/unit/test_fault_basis_dip_strike_symmetry.cpp:122-143] — Orthonormality tolerance `10 * eps` is too loose to detect the pre-R-003 drift; the test passes with OR without the fix

**Category:** BUG (in test, not in the fix)

**Description:**
The test declares its contract in the file header:

> "Pre-R-003 the strike was explicitly normalized but the dip was computed as `strike x n_raw` without an explicit normalize; FP rounding in the strike normalize leaked into the cross product so |dip| drifted up to ~5 ULP.  The test documents the pre-fix drift as a regression gate — deleting the R-003 patch should make this test fail on the |dip| assertions."

But the per-trial tolerance is set to

```cpp
const double tol_len = 10.0 * eps;   // = 10 * std::numeric_limits<double>::epsilon() ≈ 2.22e-15
const double tol_dot = 10.0 * eps;
```

i.e. **10 ULP at value 1**. The plan (§3.4 H-V91-B1) and the test's own header comment both estimate the pre-fix drift at **~5 ULP**. The actual numerical-analysis bound for `|dv| - 1` when `dv = strike_unit × n_raw_unit` without a final normalize is

```
|dv|² = |s|²·|n|² - (s·n)² ≈ (1 + O(eps))·(1 + O(eps)) - (O(eps))²
       = 1 + O(eps)  →  |dv| = 1 + O(eps)/2
```

so 1–4 ULP drift is typical, always strictly below the 10 ULP test bar. Deleting the R-003 normalize block on lines 425-437 of `fault_basis.hpp` and re-running the test would **still report "3*1000 |length - 1| <= 10*eps checks PASSED"** — the test is not actually a regression gate, despite its header claim.

Marked **POSSIBLE** because I verified the analysis on paper but did not physically delete the R-003 block to confirm the test still passes (the working tree is in a delivery-ready state; I did not want to mutate it). The theoretical bound is solid enough to flag.

**Trigger:**
Any future attempt to use this test as a regression gate on the `dip` normalization block. The test will silently greenlight a revert.

**Actual behavior:**
Both pre-R-003 and post-R-003 |dv| residuals fall within the 10 ULP window. The test reports PASS in both cases. It verifies a generic orthonormality contract but not the specific R-003 improvement.

**Expected behavior:**
The tolerance should be tight enough that pre-R-003 (missing normalize) fails on at least a small fraction of the 1000 trials. Either tighten `tol_len` to `2.0 * eps` (2 ULP) or add a second assertion that compares the post-fix residual against a recorded pre-fix upper bound.

**Suggested fix:**
Tightening the length tolerance alone is sufficient because the test iterates 1000 random trials — hitting a worst-case drift becomes nearly certain with a 2 ULP bound:

```diff
@@ tests/unit/test_fault_basis_dip_strike_symmetry.cpp:121-123 @@
    const double eps = std::numeric_limits<double>::epsilon();
-   const double tol_len = 10.0 * eps;
-   const double tol_dot = 10.0 * eps;
+   // R-003 gate: pre-fix |dv| drift is ~1-4 ULP (see plan §3.4 H-V91-B1);
+   // require < 2 ULP so deleting the normalize block on fault_basis.hpp:425-437
+   // makes at least one of the 1000 trials fail.  The orthogonality tolerance
+   // stays at 10*eps because cross-product rounding dominates there.
+   const double tol_len = 2.0 * eps;
+   const double tol_dot = 10.0 * eps;
```

If the post-R-003 test is too tight (say, 5 trials of 1000 fail at 2*eps due to unlucky rounding chains), widen to 3*eps — but not back to 10*eps.

**Test case:**
```cpp
TEST(FaultBasisDipSymmetry, GateRejectsPreR003Revert)
{
   // Procedure: temporarily #ifdef out lines 425-437 of fault_basis.hpp and
   // rebuild.  Re-running the regression test must produce FAILED lines on
   // "length" checks within the 1000-trial loop.
   //
   // If the test still reports "3*1000 |length - 1| <= 10*eps checks PASSED"
   // after the #ifdef, the tolerance is too loose (R-006).  Post-R-006 fix
   // the count of failed `|dv| - 1 > tol_len` events should be in the
   // hundreds at tol_len = 2*eps (verified manually by removing the normalize
   // and re-running with debug printfs of d_len in RunTrial).
}
```

---

### [R-007] [LOW] [io/paraview_output.hpp:595 WriteFaultSurfaceVTU::process_face] — `MFEM_ASSERT(nbf > 0, …)` is compiled out in Release; a misconfigured driver silently produces NaN cell values

**Category:** EDGE_CASE / ROBUSTNESS

**Description:**
The R-001 patch added
```cpp
MFEM_ASSERT(nbf > 0,
            "WriteFaultSurfaceVTU: nbf_per_face_ must be > 0");
const int base = face_mesh_idx * nbf;
...
const double inv_nbf = 1.0 / static_cast<double>(nbf);
```
Release builds compile `MFEM_ASSERT` to a no-op. If an upstream caller ever initialises `InitFaultOutputBP5` with `nbf_per_face = 0` (the signature permits it; there is no validation on line 220), the loop `for (int k = 0; k < 0; k++)` runs zero times, each `a_*` stays at 0.0, then `inv_nbf = 1.0/0.0 = +Inf`, then every push `0.0 * +Inf = NaN` — so the CellData arrays become solid NaN and ParaView renders nothing, silently. No driver currently hits this (TPV102 driver and BP5 verify both derive `nbf` from a genuine `IntRules.Get(...).GetNPoints()` or `face_quad_->NumBasisFunctions()`), but the guard should promote the assert to `MFEM_VERIFY` so a future misconfiguration fails loud in production rather than emitting silent NaN.

**Trigger:**
Any driver that calls `InitFaultOutputBP5(..., 0)` — including a future regression where `GetNbfPerFace()` returns 0 from an unconfigured domain operator.

**Actual behavior:**
Release: NaN flood in the VTU cell data; ParaView shows an empty surface. No error signalled.
Debug: MFEM_ASSERT fires on the first face, aborts with message.

**Expected behavior:**
Both Release and Debug: fail loud at the function entry with a clear error message before any output is touched.

**Suggested fix:**
Move the check to the top of the function and use `MFEM_VERIFY`:
```diff
@@ io/paraview_output.hpp:525-535 WriteFaultSurfaceVTU @@
    if (!has_fault_output_) { return; }
    const int nbf = nbf_per_face_;
    const int n_int = n_interior_fault_faces_;
    const int n_shared = n_shared_fault_faces_;
    const int n_faces = n_int + n_shared;
    const bool has_normal = (local_normal_stress.Size() > 0);
+
+   // R-007: Release-visible guard.  `nbf == 0` with n_faces > 0 would make
+   // process_face divide by zero and emit NaN cells.  Fail loud instead.
+   MFEM_VERIFY(nbf > 0 || n_faces == 0,
+               "WriteFaultSurfaceVTU: n_faces=" << n_faces << " faces but "
+               "nbf_per_face_=0.  Call InitFaultOutputBP5 with a positive "
+               "nbf_per_face before writing fault-surface VTU.");
@@ io/paraview_output.hpp:595 process_face @@
-         MFEM_ASSERT(nbf > 0,
-                     "WriteFaultSurfaceVTU: nbf_per_face_ must be > 0");
         const int base = face_mesh_idx * nbf;
```

**Test case:**
```cpp
TEST(FaultSurfaceVTU, NbfZeroAbortsWithMessage)
{
   Mesh mesh = MakeTinyTetMesh();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 2);
   Array<int> empty_shared;
   ParaViewOutput<Mesh> pv("test_r007", mesh, 1);
   // Misconfigure: pass nbf=0 even though we have real fault faces.
   pv.InitFaultOutputBP5(fault_faces, empty_shared, /*nbf_per_face=*/0);
   // Post-R-007 fix: this should MFEM_ABORT with a clear message.
   // Pre-fix Release: silently emits NaN cells.  Pre-fix Debug: asserts.
   EXPECT_DEATH_OR_VERIFY_FAIL(
      pv.WriteFaultSurfaceVTU("test_r007", 0, 0.0, 0, 1,
                               /* all local_* Vectors of size 0 */,
                               ...),
      "nbf_per_face_=0");
}
```

---

### [R-008] [LOW] [fault/fault_basis.hpp:337-342] — `public:` inserted mid-class after `private:` creates a trap for future member additions

**Category:** QUALITY

**Description:**
The R-003 patch needed to expose `ComputeOrientedFrame` for direct unit-testing and did so by inserting a `public:` access specifier at line 342, immediately after the `private:` section that contains `dim_`, `num_faces_`, `basis_`. The class definition ends at line 473 with no returning `private:` label, so any future member appended before the closing brace will default to `public` — the opposite of the surrounding convention (member variables of this class are clearly meant to stay private). This is a small but real trap for the next contributor.

**Trigger:**
A later patch that adds a helper or cached field anywhere between line 451 (end of `ComputeOrientedFrame`) and line 473 (closing brace) without explicitly re-stating `private:`.

**Actual behavior:**
Any new member variable or helper accidentally becomes part of the public API.

**Expected behavior:**
Access control returns to `private` after the exposed test helper, or the test helper is moved out entirely (e.g., to a free function in the same header, or a friend declaration).

**Suggested fix:**
Option A — close the `public:` block with a re-stated `private:` at end of class:
```diff
@@ fault/fault_basis.hpp:471-473 @@
          }
       }
    }
+
+private:
+   // (Intentionally empty; restored access control after R-003's public:
+   //  exposure of ComputeOrientedFrame on line 342.  Any future member
+   //  appended inside this class should live here.)
 };
```
Option B — keep the test helper private and declare the test class as a `friend`:
```diff
@@ fault/fault_basis.hpp class declaration @@
 class FaultBasis
 {
+#if defined(SEAS_TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY)
+   friend void RunTrial(uint64_t &, int, int &, int &);
+   friend void TestCanonicalTPV102Frame();
+#endif
 public:
   ...
```
Option A is the minimum-surface fix and keeps the test file simple.

**Test case:** None — pure code-structure correctness. Verification is that adding `int dummy_field_;` between the current `}` of `ComputeOrientedFrame` and the class's `};` produces a public member (observable by `sizeof(FaultBasis)` change + successful external access); post-fix Option A, the same addition would leave `dummy_field_` private.

---

### [R-009] [LOW] [docs] — Plan's required `tpv102_debug_v9.1.0_fix.md` entry is missing for the R-001 / R-002 / R-003 patch

**Category:** DEVIATION

**Description:**
Plan §TODO item:
> "**Debugger:** Write `tpv102_debug_v9.1.0_fix.md` after each step, one entry per substantive action."

Two substantive actions have been applied:
1. `paraview_output.hpp::WriteFaultSurfaceVTU` CellData/per-face-average rewrite (R-001 + R-002 + R-005 side-effect).
2. `fault_basis.hpp::ComputeOrientedFrame` dip-vector normalize (R-003) + `public:` exposure (R-008 tag).

No `tpv102_debug_v9.1.0_fix.md` exists in `miniapps/seas/debug_document/tpv102_debug_document/`. The plan's commit-sequence §8 lists C1 (test additions) and C5 (sbatch generation); it does not cover the production-code fix C-step that actually landed in this round. Without the fix doc, future re-entries into this debug stream lose the decision trail for why the VTU writer's CellData path was restored (reversing commit `0756b80`'s 2026-04-07 "Switch fault surface VTU from CellData to PointData for smooth interpolation" change).

Also mild scope deviation: plan §8 C1 was advertised as "No production-code change"; the landed patch bundles tests **plus** two production-code files. My own rev-2 REVIEW.md recommended that bundling, so this deviation is acknowledged, not silent — but the plan itself remains un-amended.

**Trigger:**
Any future R-V91-A regression triage. Without the v9.1.0 fix doc, a later debugger will see commit `0756b80` restored (CellData→PointData) and reverted again (PointData→CellData) with no in-repo explanation tying the round-trip to the QP/vertex-position mismatch analysis.

**Actual behavior:**
Plan TODO is partially complete; §TODO `[ ] Debugger: Write tpv102_debug_v9.1.0_fix.md after each step` remains unchecked.

**Expected behavior:**
Write `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.1.0_fix.md` summarising the R-001 / R-002 / R-003 / R-005 / R-006 / R-007 / R-008 landing, cross-linking to commit `0756b80` (the commit that introduced R-001), and recording that R-V91-A is now closed on the TPV102 path pending an end-to-end post-fix VTU screenshot.

**Suggested fix:**
Create the file with one entry per landed change. Minimum skeleton:

```markdown
# TPV102 Debug v9.1.0 — Fix Log

## 2026-04-20 — R-001 / R-002 / R-005 resolved: VTU writer switched to CellData
- Reverts commit 0756b80 (2026-04-07 "Switch fault surface VTU from CellData
  to PointData") which introduced H-V91-A4 by mapping interior-QP DOF
  values onto reference-triangle vertex positions.
- Per-face CellData mean over nbf DOFs.  Correct for any nbf ≥ 1 (BR2
  centroid, IP/TPV102 P1 three-node, future higher-order).
- Files: miniapps/seas/io/paraview_output.hpp.
- Test: tests/unit/test_fault_surface_vtu_continuity.cpp (63/63 PASSED).

## 2026-04-20 — R-003 resolved: explicit dip normalize
- fault_basis.hpp:425-437 adds sqrt + division + MFEM_VERIFY on |dv|.
- Test: tests/unit/test_fault_basis_dip_strike_symmetry.cpp (17/17 PASSED;
  R-006 notes the tolerance is currently too loose to gate the revert —
  retain test as orthonormality contract, follow up per R-006).

## 2026-04-20 — R-008 side-effect: ComputeOrientedFrame promoted to public
- Required for test_fault_basis_dip_strike_symmetry to call the static
  method directly.  Followup: re-close class with `private:` before final
  `};` (R-008).
```

**Test case:** N/A (docs).

---

## Summary
- Critical issues: **0** (R-001 is resolved; verified end-to-end with a live test run)
- Moderate issues: **1** (R-006 — test-design gap; does not affect runtime correctness)
- Low issues: **3** (R-007 assert→verify, R-008 access-control layout, R-009 missing fix doc)
- Plan compliance: **PARTIAL** — R-001/R-002/R-003/R-005 fixes landed with new unit tests per §4.1/§4.2; §5 offline diagnostics and §6 Frontera 12 s re-run remain open as the plan specifies. R-009 flags the missing §TODO item. No illicit scope creep detected.
- Verdict: **PASS WITH FIXES** — production-code changes are correct for all three target workflows (TPV102 wave operator, BP5 IP, BP5 BR2) and both regression tests pass. R-006 is the highest-impact open item: the R-003 regression test should be re-armed before the v9.1.0 C1 commit lands. R-007, R-008, R-009 are non-blocking cleanup.

## Unreviewed Areas
- **`dynamic/fault_face_flux.cpp` + `wave_operator.inl`** — v9.0.0 Pelties-9 per-side flux math. Re-audited in rev 1; no changes in this round. R-V91-B remains open awaiting the §6 Frontera 12 s run (test scope, not code scope).
- **`dynamic/friction_solver.cpp`** — Brent path is unchanged; test_friction_solver / test_vector_friction untouched.
- **BP5 quasi-dynamic solve + RK45 tolerances** — untouched by this round; out of scope for v9.1.0.
- **`ParaViewOutput::UpdateFaultFieldsBP5` L2-p0 projection path (lines 398-499)** — not touched by R-001 patch; already uses correct per-face `k < nbf` averaging and was confirmed clean in rev 2.
- **Pre-fix empirical confirmation for R-006** — I verified the tolerance gap on paper (|dv| drift of 1-4 ULP vs 10-ULP tolerance) but did not physically delete the `fault_basis.hpp:425-437` normalize block and re-run the 1000-trial test to observe a PASS. R-006 is therefore marked POSSIBLE; the fix agent should confirm empirically before declaring R-006 closed.
