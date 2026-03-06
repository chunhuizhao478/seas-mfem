# Phase 2a Implementation Review: FaultBasis and DomainOperator Interface

**Date**: 2026-03-03 (revision 3 — all code issues and test gaps resolved)
**Reviewing**: `fullelasticity_phase2a_plan_03022026.md` vs. current implementation
**Files checked**: `fault/fault_basis.hpp`, `domain/domain_operator.hpp`,
`tests/unit/test_fault_basis.cpp`, `tests/unit/test_domain_operator_interface.cpp`,
`tests/unit/test_macros.hpp`, `Makefile`
**References**: SCEC BP5-QD spec, Tandem `src/geometry/Curvilinear.cpp:facetBasis()`,
Tandem `app/localoperator/AdapterBase.cpp:prepare()`,
Tandem `app/localoperator/ElasticityAdapter.cpp`

---

## 1. Overall Assessment

Phase 2a is complete. The `FaultBasis` class and `DomainOperator` interface
extensions are implemented and tested. The FaultBasis algorithm correctly
follows Tandem's `facetBasis()`. The DomainOperator defaults are backward
compatible with existing 2D antiplane operators.

All code issues from revision 1 have been resolved. Test coverage has been
expanded from 8 to 13 FaultBasis test functions (187 assertions) and domain
interface tests migrated to shared `test_macros.hpp`.

**Changes since revision 1:**
- 2D collinearity check added (Issue 2.1 → COMPLETED)
- Parameter renamed `dof` → `fi` (Issue 2.2 → COMPLETED)
- MFEM_ASSERT bounds checking added to all methods (Issue 2.3 → COMPLETED)
- `test_domain_operator_interface.cpp` migrated to `test_macros.hpp` (Issue 4.1 → COMPLETED)
- EmbedSlip branch hoisted outside loop (Issue 4.2 → COMPLETED)
- Plan updated: NumComponents/Dimension marked as pre-existing (Issue 3.1 → COMPLETED)
- Plan updated: test counts corrected to "13+5 after review fixes" (Issue 3.2 → COMPLETED)
- 3 new FaultBasis tests added (rev 1→2): EmptyFaultFaces, MultiFace2D, OrientationFlip2D
- GetFaultBasis nullptr check added to domain interface tests 4 and 5

**Changes since revision 2:**
- 2 new FaultBasis tests added: TestNonUnitInputVectors (test 12), TestNearCollinearUpNormal (test 13)
- Plan updated to reflect 13 test functions with 187 assertions
- Non-unit ref_normal/up test gap → COMPLETED (test 12)
- Collinear up/n boundary test gap → COMPLETED (test 13)

**All Phase 2a items resolved.** 4 Phase 2b concerns moved to Phase 2b plan
as carry-forward CF-1 through CF-4 (see 7.1, summary table).

---

## 2. Bugs / Potential Bugs

### 2.1 ~~[LOW] 2D FaultBasis lacks collinearity check for `up ∥ n`~~ COMPLETED

**File**: `fault/fault_basis.hpp:147-148`

**Fixed.** The 2D branch now has a collinearity check matching the 3D branch:
```cpp
MFEM_VERIFY(std::abs(cross) > 1e-12,
            "2D: Up vector and normal are nearly collinear");
```

This matches Tandem's behavior of throwing for the collinear case in both
2D and 3D.

### 2.2 ~~[INFO] ProjectTraction/EmbedSlip/NormalStress parameter named `dof`~~ COMPLETED

**Fixed.** All four methods now use `fi` (face index) as the parameter name,
matching the plan interface:
```cpp
const FaultBasisData &GetBasis(int fi) const
void ProjectTraction(int fi, ...) const
void EmbedSlip(int fi, ...) const
real_t NormalStress(int fi, ...) const
```

### 2.3 ~~[INFO] No bounds checking on face index~~ COMPLETED

**Fixed.** All four methods now have MFEM_ASSERT bounds checking:
```cpp
MFEM_ASSERT(fi >= 0 && fi < num_faces_,
            "FaultBasis::GetBasis: face index " << fi
            << " out of range [0, " << num_faces_ << ")");
```
Same pattern in `ProjectTraction`, `EmbedSlip`, `NormalStress`.

---

## 3. Deviations from Plan

### 3.1 ~~NumComponents() and Dimension() are pure virtual, not defaulted~~ COMPLETED

**Fixed in plan.** The plan now correctly classifies `NumComponents()` and
`Dimension()` as "Pre-existing pure virtual methods (unchanged)" and clearly
separates them from the "NEW Phase 2a methods" section. Only
`NumSlipComponents()`, `GetFaultCoords2D()`, `GetFaultBasis()`, and
`GetOffFaultDisplacement()` are listed as Phase 2a additions.

### 3.2 ~~Test count: "113 tests" vs 8 test functions~~ COMPLETED

**Fixed in plan.** The plan status line now reads: "13 FaultBasis test functions
with 187 assertions + 5 DomainOperator interface test functions with 32
assertions". The test section lists all 13 FaultBasis tests and 5 DomainOperator
test functions.

The plan's "Files Created/Modified" table has also been updated to
"13 test functions, 187 assertions".

### 3.3 ~~Parameter naming: `fi` (plan) vs `dof` (code)~~ COMPLETED

**Fixed.** Code now uses `fi` matching the plan. See Issue 2.2 above.

### 3.4 ~~[INFO] Plan `GetBasis` interface uses `face_idx`, code uses `fi`~~ COMPLETED

**Fixed in plan.** The plan now shows `GetBasis(int fi)` (line 70), consistent
with `ProjectTraction(int fi, ...)`, `EmbedSlip(int fi, ...)`, and the code.

---

## 4. Code Quality Issues

### 4.1 ~~[LOW] test_domain_operator_interface.cpp duplicates test macros~~ COMPLETED

**Fixed.** The file now uses `#include "test_macros.hpp"` (line 22) and the
shared `TEST_ASSERT`, `TEST_NEAR`, `TEST_PRINT_RESULTS()` macros. Duplicate
macro definitions have been removed. Test functions changed from `bool` return
to `void` return (matching the other test files). The `RUN_TEST` wrapper macro
has been replaced by direct function calls in `main()`.

**Note**: The shared `test_macros.hpp` `TEST_ASSERT` does NOT return early on
failure — it increments `num_failed` and continues. This means if a mesh
creation fails (`TEST_ASSERT(mesh != nullptr, ...)`), subsequent code in the
same function will still execute and potentially crash on the null mesh.
In practice, `CreateTestMesh()` should never fail, so this is a theoretical
concern only. If early-exit-on-failure behavior is desired in the future,
consider an `ABORT_ON_FAILURE` variant macro.

### 4.2 ~~[INFO] EmbedSlip has dimension branch inside loop~~ COMPLETED

**Fixed.** `EmbedSlip` now has the dimension branch outside the loop,
matching `ProjectTraction` style:
```cpp
for (int d = 0; d < dim_; d++)
{
   delta_u_global[d] = slip_local[0] * b.tangent1[d];
}
if (dim_ == 3)
{
   for (int d = 0; d < dim_; d++)
   {
      delta_u_global[d] += slip_local[1] * b.tangent2[d];
   }
}
```

---

## 5. Correctness Verification Against Tandem

### 5.1 FaultBasis Algorithm: CORRECT

Matches Tandem's `Curvilinear::facetBasis()`:

| Step | Tandem | SEAS-MFEM | Match |
|------|--------|-----------|-------|
| Normal computation | `det(J) * J^{-T} * ref_normal` | `CalcOrtho(J, n_raw)` | ✓ (equivalent) |
| Orientation | `dot(ref_normal, normal) < 0 → flip` | `dot(n_raw, ref_normal) < 0 → Neg()` | ✓ |
| Normalize | `n / |n|` | `n_raw / n_len` | ✓ |
| Strike | `normalize(up × n)` | `normalize(up × n)` | ✓ |
| Dip | `s × n` | `s × n` | ✓ |
| Storage order | `[n, dip, strike]` columns | `normal, tangent1=dip, tangent2=strike` | ✓ |
| 2D tangent | `sgn(up×n) * perp(n)` | `(cross >= 0 ? 1 : -1) * perp(n)` | ✓* |
| Collinearity (3D) | `throw if |s| < 10000*eps` | `MFEM_VERIFY(s_len > 1e-12)` | ✓ |
| Collinearity (2D) | `throw if |sgn| < colinear_tol` | `MFEM_VERIFY(|cross| > 1e-12)` | ✓ (fixed, see 2.1) |

*The 2D sign convention produces identical results for all non-degenerate cases.
The collinear case is now caught by MFEM_VERIFY in both 2D and 3D.

### 5.2 Sign Flip Handling: Subtle Difference (Not a Bug)

Tandem (AdapterBase.cpp):
1. Reorients normal to match `ref_normal`
2. Computes basis from the reoriented normal
3. **Additionally flips the entire basis matrix** when `sign_flipped == true`

SEAS-MFEM:
1. Reorients normal to match `ref_normal`
2. Computes basis from the reoriented normal
3. **No additional flip**

Tandem's extra flip is for its DG formulation where the two sides of a fault
face see different orientations. SEAS-MFEM's approach is correct for
single-sided operations (one consistent normal direction). For Phase 2b
(`ElasticityDomainOperator`), the DG numerical flux will need to handle
side selection separately — the FaultBasis need not encode it.

### 5.3 Per-Face vs Per-Quadrature-Point: Known Simplification

Tandem stores one basis per (face, quadrature point), supporting curved faces
with variable geometry. SEAS-MFEM stores one basis per face (constant normal
per planar face). For the BP5 planar fault this is equivalent. Already
documented in Phase 1 review.

### 5.4 Collinearity Tolerance

- Tandem: `10000 * std::numeric_limits<double>::epsilon()` ≈ `2.2e-12`
- SEAS-MFEM: `1e-12` (both 3D and 2D, see Issue 2.1 fix)

Both are adequate. The values are close.

### 5.5 ProjectTraction / EmbedSlip: CORRECT

SEAS-MFEM uses direct dot products, which is equivalent to Tandem's
matrix-based projection for piecewise-constant bases on planar faces:
- `tau_local = basis^T * tau_global` ↔ `tau[i] = traction · tangent_i`
- `Δu = basis * (0, s1, s2)^T` ↔ `delta_u = s1*tangent1 + s2*tangent2`

Tandem additionally computes a fault mass matrix and its inverse for
integrating tractions from quadrature points to DOFs. This is needed for
higher-order elements with multiple DOFs per face. For piecewise-constant
DG (one DOF per face), the mass matrix is diagonal and the direct approach
is equivalent.

### 5.6 DomainOperator Extensions: CORRECT

- `NumSlipComponents()` defaults to 1 (backward compatible) ✓
- `GetFaultCoords2D()` default calls `GetFaultDepths()`, sets x2=0 ✓
- `GetFaultBasis()` defaults to `nullptr` ✓
- `GetOffFaultDisplacement()` default sets size to 0 ✓
- Forward declaration used for `FaultBasis` (no coupling) ✓

### 5.7 NormalStress Convention: CORRECT

`sigma_n = -traction · normal` → positive in compression.
This matches the convention expected by rate-and-state friction laws
where `sigma_n > 0` means compressive normal stress.

### 5.8 Makefile: CORRECT

- `fault/fault_basis.hpp` in `FAULT_HEADERS` ✓
- `test_fault_basis` object depends on `$(FAULT_HEADERS) $(DOMAIN_HEADERS)` ✓
- `test_domain_interface` object depends on `$(SEAS_HEADERS)` ✓
- Both targets in `SEQ_MINIAPPS` and `test:` target ✓

---

## 6. Test Coverage Assessment

### FaultBasis Tests (`test_fault_basis.cpp`) — 13 test functions, 187 assertions

| # | Test | What's Verified | Status |
|---|------|-----------------|--------|
| 1 | `TestAxisAligned3D` | n,t1,t2 vectors, orthonormality, SCEC convention | ✓ |
| 2 | `TestProjectTraction` | Global → local projection, pure-normal isolation | ✓ |
| 3 | `TestEmbedSlip` | Local → global embedding, dip/strike/mixed | ✓ |
| 4 | `TestRoundTrip` | Project → embed consistency, normal discarded | ✓ |
| 5 | `TestNormalStress` | Tensile, compressive, purely tangential | ✓ |
| 6 | `Test2DDegeneracy` | Single tangent, project, embed, normal stress | ✓ |
| 7 | `TestOrientationFlip` | Opposite normals, orthonormality preserved (3D) | ✓ |
| 8 | `TestMultiFace3D` | 3 faces: orthonormality, consistency, indexing | ✓ |
| 9 | `TestEmptyFaultFaces` | Compute with 0 faces, NumFaces()==0 | ✓ |
| 10 | `TestMultiFace2D` | 4-element quad mesh, 3 interior faces, 2D consistency | ✓ |
| 11 | `TestOrientationFlip2D` | 2D ref_normal flip, orthonormality preserved | ✓ |
| 12 | `TestNonUnitInputVectors` | Scaled ref_normal/up produce same basis (3D+2D) | ✓ NEW |
| 13 | `TestNearCollinearUpNormal` | Near-collinear up/n still orthonormal (3D+2D) | ✓ NEW |

### DomainOperator Interface Tests (`test_domain_operator_interface.cpp`) — 5 functions, 32 assertions

| # | Test | What's Verified | Status |
|---|------|-----------------|--------|
| 1 | `test_num_slip_components_antiplane` | `NumSlipComponents() == 1` | ✓ |
| 2 | `test_get_fault_coords_2d_default` | x2=0, x3=depths | ✓ |
| 3 | `test_off_fault_displacement_default` | Returns empty | ✓ |
| 4 | `test_polymorphic_dispatch` | Base pointer dispatch + GetFaultBasis nullptr | ✓ UPDATED |
| 5 | `test_bdrload_interface_methods` | BdrLoad inherits defaults + GetFaultBasis nullptr | ✓ UPDATED |

Note: The plan previously listed "6 tests" counting GetFaultBasis nullptr as
a separate case; the plan now correctly lists 5 test functions.

### Test Gaps (Updated)

- [x] **GetFaultBasis default returns nullptr**: Tested in both
      `test_polymorphic_dispatch` (line 171) and `test_bdrload_interface_methods`
      (line 217). COMPLETED.
- [x] **Empty fault_faces array**: `TestEmptyFaultFaces` added (test 9).
      Verifies `NumFaces()==0`, `Dimension()==3`, `NumTangentComponents()==2`.
      COMPLETED.
- [x] **2D multi-face mesh**: `TestMultiFace2D` added (test 10). 4-element
      quad mesh with 3 interior faces, verifies 2D branch consistency.
      COMPLETED.
- [x] **2D orientation flip**: `TestOrientationFlip2D` added (test 11).
      Exercises 2D `cross` sign logic under ref_normal flip. COMPLETED.
- [x] **Non-unit ref_normal and up vectors**: `TestNonUnitInputVectors`
      added (test 12). Scaled ref_normal (7.5x) and up (42x, 3.14x, 100x)
      produce identical basis to unit inputs, tested in both 3D and 2D.
      COMPLETED.
- [x] **Collinear up and normal**: `TestNearCollinearUpNormal` added (test 13).
      Tests near-collinear boundary (up just above `1e-12` tolerance) in both
      3D and 2D. Verifies orthonormality is preserved even with nearly-degenerate
      input. Uses relaxed `1e-10` tolerance for dot products. COMPLETED.
      Note: Exact-collinear case cannot be tested without MFEM_USE_EXCEPTIONS
      (MFEM_VERIFY calls `abort()`).
- ~~**FaultBasis rotated fault**~~: SKIP — moved to Phase 2b plan as CF-1.

---

## 7. Suggestions

### 7.1 ~~Phase 2b Readiness Concerns~~ SKIP — moved to Phase 2b plan

All three items below have been transferred to
`fullelasticity_phase2b_plan_03022026.md` as carry-forward items CF-2, CF-3, CF-4.

1. ~~**DG two-sided operations**~~: → Phase 2b CF-2
2. ~~**Mass matrix for traction integration**~~: → Phase 2b CF-3
3. ~~**Per-face granularity**~~: → Phase 2b CF-4

### 7.2 Compilation Verification

Run the full test suite to verify no regressions:
```bash
make seas_test_fault_basis seas_test_domain_interface && \
make test-fault-basis && make test-domain-interface && \
make test-friction && make test-psi-state
```

---

## 8. Summary Table

| Item | Status | Severity | Action |
|------|--------|----------|--------|
| FaultBasis algorithm | **Correct** | -- | Matches Tandem |
| ProjectTraction | **Correct** | -- | Direct dot product OK for DG0 |
| EmbedSlip | **Correct** | -- | Direct embedding OK for DG0 |
| NormalStress convention | **Correct** | -- | Positive compression ✓ |
| DomainOperator defaults | **Correct** | -- | Backward compatible |
| Forward declaration | **Correct** | -- | No domain→fault coupling |
| SCEC convention encoding | **Correct** | -- | Tested in TestAxisAligned3D |
| Makefile | **Correct** | -- | Dependencies wired correctly |
| 2D collinearity check | **COMPLETED** | ~~LOW~~ | MFEM_VERIFY added (see 2.1) |
| Parameter naming `fi` | **COMPLETED** | ~~INFO~~ | Renamed from `dof` (see 2.2) |
| Bounds checking on face index | **COMPLETED** | ~~INFO~~ | MFEM_ASSERT added (see 2.3) |
| test_domain_interface macros | **COMPLETED** | ~~LOW~~ | Migrated to test_macros.hpp (see 4.1) |
| EmbedSlip branch style | **COMPLETED** | ~~INFO~~ | Branch hoisted outside loop (see 4.2) |
| Plan: NumComponents/Dimension | **COMPLETED** | ~~INFO~~ | Plan corrected (see 3.1) |
| Plan: test counts | **COMPLETED** | ~~INFO~~ | Plan updated to 13+5 (see 3.2) |
| GetFaultBasis nullptr test | **COMPLETED** | ~~INFO~~ | Added in tests 4 and 5 |
| Empty fault_faces test | **COMPLETED** | ~~INFO~~ | TestEmptyFaultFaces added (test 9) |
| 2D multi-face test | **COMPLETED** | ~~INFO~~ | TestMultiFace2D added (test 10) |
| 2D orientation flip test | **COMPLETED** | ~~INFO~~ | TestOrientationFlip2D added (test 11) |
| Non-unit ref_normal/up test | **COMPLETED** | ~~INFO~~ | TestNonUnitInputVectors added (test 12) |
| Collinear up/n boundary test | **COMPLETED** | ~~INFO~~ | TestNearCollinearUpNormal added (test 13) |
| Plan: GetBasis `face_idx` vs `fi` | **COMPLETED** | ~~INFO~~ | Plan now uses `fi` (see 3.4) |
| Plan: "113 tests" in files table | **COMPLETED** | ~~INFO~~ | Plan updated to 13+187 (see 3.2) |
| Rotated fault test | **SKIP** | INFO | Moved to Phase 2b CF-1 |
| Sign flip for DG two-sided ops | **SKIP** | INFO | Moved to Phase 2b CF-2 |
| Mass matrix for higher-order DG | **SKIP** | INFO | Moved to Phase 2b CF-3 |
| Per-face granularity | **SKIP** | INFO | Moved to Phase 2b CF-4 |
