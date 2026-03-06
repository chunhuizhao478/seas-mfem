# Phase 2b Implementation Review: DG Elasticity BR2 Integrator + ElasticityDomainOperator

**Date**: 2026-03-03
**Last updated**: 2026-03-03 (revision 4 — all code quality issues resolved, all test suggestions addressed)
**Reviewing**: `fullelasticity_phase2b_plan_03022026.md` vs. current implementation
**Files checked**: `integrator/dg_elasticity_br2_integrator.hpp`,
`domain/elasticity_operator.hpp`, `tests/unit/test_elasticity_br2.cpp`,
`tests/unit/test_elasticity_operator.cpp`, `Makefile`
**References**: SCEC BP5-QD spec, Tandem `app/localoperator/Elasticity.h`,
Tandem `app/kernels/elasticity.py`, Tandem `app/localoperator/AdapterBase.cpp`

---

## 1. Overall Assessment

Phase 2b is complete — all code issues from previous revisions have been
resolved, and test coverage has been significantly expanded.

**Revision 4 status:**
- **0 open code issues** (all 8 from revisions 1-3 are COMPLETED, including §4.6 naming)
- **0 missing critical tests** (all 6 test gaps now covered)
- **1 suggested test improvement** remaining (full-system patch test)
- Implementation is ready for integration testing with the SEAS time-stepper

**Change summary since revision 3:**
- §4.6 (sigma\_ naming in operator): COMPLETED — member renamed to `epsilon_` (line 151)
- §6 traction non-zero assertion: COMPLETED — `TEST_ASSERT(trac_norm > 1e-6, ...)`
  now present in both TestNonZeroSlipIP (line 400) and TestNonZeroSlipBR2 (line 450)

---

## 2. Bugs / Potential Bugs — All Resolved

### 2.1 [COMPLETED] BR2 face matrix DOF ordering

Code uses correct `i * ndof + k` (component-major, byNODES) in all 5 blocks
(lines 489, 527, 565, 603, 791). Consistent with MFEM's `DGElasticityIntegrator`.
New TestDOFOrdering test verifies positive definiteness and symmetry at order 1.
New TestPatchTest verifies BR2 and IP agree for continuous functions. ✓

### 2.2 [COMPLETED] IsFaultFace3D bounds check

Full x1/x2/x3 bounds check present (lines 264-277). ✓

### 2.3 [COMPLETED] BR2 penalty geometry-dependent

Geometry-dependent penalty in all 4 locations: interior integrator (line 128),
boundary integrator (line 649), BR2 slip (line 613), BR2 Dirichlet (line 1042).
New TestTetElements validates tet assembly end-to-end. ✓

### 2.4 [INFO] ComputeTraction uses average flux only

ComputeTraction computes `{σ(u)}·n` (average of interior gradients).
Not a practical concern: both BP1 and BP5 use DG order 1 operationally,
where interior gradients are non-zero and the average flux correctly captures
traction. DG0 is only used in unit tests for matrix validation. ✓

---

## 3. Deviations from Plan

### 3.1 Plan status line needs update

**Plan**: "COMPLETE (14 BR2 tests + 32 elasticity operator tests passing)"

**Current actual**:
- `test_elasticity_br2.cpp`: **10 test functions** with **32 assertions** (780 lines)
- `test_elasticity_operator.cpp`: **12 test functions** with **~42 assertions** (651 lines)

**Recommendation**: Update plan status to:
"COMPLETE (10 BR2 test functions with 32 assertions + 12 elasticity operator
test functions with ~42 assertions)"

### 3.2 IsFaultFace3D signature differs from plan

**Plan**: `bool IsFaultFace3D(const Vector &center) const`
**Code**: `bool IsFaultFace3D(FaceElementTransformations *FTr) const`

Reasonable API change (encapsulates coordinate extraction). Full bounds
check is present. ✓

### 3.3 File line counts

| File | Plan says | Actual |
|------|-----------|--------|
| `dg_elasticity_br2_integrator.hpp` | ~850 lines | 819 lines ✓ |
| `elasticity_operator.hpp` | ~900 lines | 1336 lines ✗ |
| `test_elasticity_br2.cpp` | ~406 lines | 780 lines ✗ |
| `test_elasticity_operator.cpp` | ~375 lines | 651 lines ✗ |

The elasticity operator is larger than planned (~1336 vs ~900, extra ~436 lines
from BR2 Dirichlet loading and detailed traction computation). Test files
grew significantly from expanded test coverage.

---

## 4. Code Quality Issues — All Resolved

### 4.1 [COMPLETED] TestNormal duplication

Boundary class inherits from interior (line 91). No duplication. ✓

### 4.2 [COMPLETED] sigma\_ mutation in AssembleFaceMatrix

Both integrators use local `const real_t sigma` (lines 128, 649). Thread-safe. ✓

### 4.3 [COMPLETED] mu\_all\_ naming

Consistent naming `mu_all` (no underscore) in both integrators (lines 149, 665). ✓

### 4.4 [COMPLETED] SetupBoundaryMarkers dead code

Dead center computation removed. Method is now clean (lines 203-221):
```cpp
for (int be = 0; be < mesh_.GetNBE(); be++)
{
   int attr = mesh_.GetBdrAttribute(be);
   if (attr == 3 || attr == 4)
   {
      dirichlet_bdr_marker_[attr - 1] = 1;
   }
}
```
No unused variables or abandoned code. ✓

### 4.5 [COMPLETED] Dead sigma\_ member in BR2 integrator

The `sigma_` member has been removed from `DGElasticityBR2Integrator` entirely.
Class members (lines 67-70) are now only `lambda_`, `mu_`, `epsilon_`,
`elem_mass_inv_`, `dim_`. Constructor body is empty `{}`. ✓

### 4.6 [COMPLETED] sigma\_ naming in ElasticityDomainOperator

Member renamed from `sigma_` to `epsilon_` (line 151):
```cpp
real_t epsilon_;  // SIPG sign = -1
```
Initialized at line 86: `epsilon_ = -1.0;  // SIPG`. Passed as the `epsilon`
argument to integrator constructors. Naming now matches DG convention. ✓

---

## 5. Correctness Verification Against Tandem

### 5.1 test\_normal operator: CORRECT ✓
Formula matches Tandem `elasticity.py` identically.

### 5.2 BR2 penalty: CORRECT ✓
`2*dim` for hex, `dim+1` for tet. Matches Tandem's `NumFacets`.

### 5.3 Lifting computation: CORRECT ✓
Math follows scalar BR2 generalized for vector elasticity.

### 5.4 Sign conventions: CORRECT ✓
All 4 interior blocks and boundary block match plan.

### 5.5 Slip BC (IP): CORRECT ✓
Embedding, sign, symmetry, penalty, byNODES ordering all verified.

### 5.6 Slip BC (BR2): CORRECT ✓
Lifted slip with test\_normal coupling, geometry-dependent penalty.

### 5.7 Traction computation: CORRECT ✓
Average flux approach, strain/stress computation, FaultBasis projection.
Appropriate for DG1 (operational order).

### 5.8 Dirichlet loading: CORRECT ✓
Both IP and BR2 methods, ±x2 walls with u\_D = (0, ±Vp·t/2, 0).

### 5.9 Boundary condition table: CORRECT ✓
All 6 faces with correct attribute mapping.

### 5.10 FE space setup: CORRECT ✓
DG\_FECollection, vdim=3, Ordering::byNODES, scalar FES for BR2 Minv.

### 5.11 Fault detection: CORRECT ✓
Geometric detection with full x1/x2/x3 bounds.

### 5.12 Makefile: CORRECT ✓
All targets, dependencies, and header lists verified.

---

## 6. Test Coverage Assessment

### Covered (all 6 previous gaps now addressed):

**BR2 integrator tests** (10 functions, 32 assertions):
- [x] Face matrix symmetry (SIPG ε=-1) — TestSymmetry
- [x] Face matrix non-zero norm — TestNonZero
- [x] Boundary face symmetry — TestBoundarySymmetry
- [x] Full bilinear form assembly + symmetry — TestFullAssembly
- [x] Order 0 face matrix (6×6) — TestOrder0
- [x] BP5 material parameters — TestMaterialParams
- [x] **DOF ordering at order 1** (positive definiteness + symmetry) — TestDOFOrdering *(NEW)*
- [x] **BR2 vs IP at order 0** (same block structure) — TestBR2vsIPOrder0 *(NEW)*
- [x] **Patch test** (constant u → zero, linear u → BR2=IP) — TestPatchTest *(NEW)*
- [x] **Tet elements** (assembly, symmetry, correct size) — TestTetElements *(NEW)*

**Elasticity operator tests** (12 functions, ~42 assertions):
- [x] Construction: NumComponents, Dimension, NumSlipComponents — TestConstruction
- [x] Fault detection at x1=0 — TestFaultDetection
- [x] FaultBasis orthonormality through operator — TestFaultBasis
- [x] Zero slip + zero loading → zero displacement — TestZeroSlipEquilibrium
- [x] Zero slip → near-zero traction — TestTractionExtraction
- [x] Stiffness assembly succeeds (IP and BR2) — TestStiffnessAssembly
- [x] FaultGeometry 3D with BP5Params — TestFaultGeometry3D
- [x] **Non-zero slip → non-zero displacement (IP)** — TestNonZeroSlipIP *(NEW)*
- [x] **Non-zero slip → non-zero displacement (BR2)** — TestNonZeroSlipBR2 *(NEW)*
- [x] **Dirichlet loading at t>0 (IP and BR2)** — TestDirichletLoading *(NEW)*
- [x] **BR2 vs IP comparison** (same fault count, similar norms) — TestBR2vsIP *(NEW)*
- [x] **BR2 default method** — TestBR2Default *(NEW)*

### Minor improvements suggested:

1. **[COMPLETED] Traction non-zero assertion**: Both TestNonZeroSlipIP (line 400)
   and TestNonZeroSlipBR2 (line 450) now assert `trac_norm > 1e-6`:
   ```cpp
   TEST_ASSERT(trac_norm > 1e-6, "IP: non-zero slip produces non-zero traction");
   TEST_ASSERT(trac_norm > 1e-6, "BR2: non-zero slip produces non-zero traction");
   ```
   ✓

2. **[LOW] Full-system patch test**: TestPatchTest validates the face
   integrator matrix action (K\_faces \* u). A stronger variant would verify
   that the full stiffness K (volume + faces) applied to a linear displacement
   gives the correct boundary-only residual. Not critical — current coverage
   is already good.

---

## 7. Phase 2a Carry-Forward Items (from plan §CF-1 to CF-4)

### CF-1: Rotated fault FaultBasis test — OPEN
All FaultBasis tests use axis-aligned faults. TestPatchTest provides
good DOF ordering validation but doesn't exercise non-axis-aligned geometry.

### CF-2: DG two-sided sign handling — PARTIALLY ADDRESSED
TestNonZeroSlipIP and TestNonZeroSlipBR2 verify that non-zero slip produces
non-zero displacement, implying sign handling works. A more targeted test
with known analytical solution would be stronger.

### CF-3: Mass matrix for higher-order DG traction — OPEN
For DG order > 0, traction at multiple face quadrature points needs mass
matrix projection. Current centroid evaluation is adequate for DG1 with
one DOF per element.

### CF-4: Per-face vs per-quadrature-point basis — OPEN
FaultBasis stores one basis per face (planar faces). Curved faults would
need per-quadrature-point storage.

---

## 8. Design Observations

### 8.1 mutable pattern: Acceptable
MFEM's non-const Coefficient::Eval() requires `mutable`. Standard MFEM pattern.

### 8.2 Lazy stiffness assembly: Good
Assembled on first Solve(), not in constructor.

### 8.3 Two FE spaces: Correct
Scalar FES for BR2 mass matrix inverses, vector FES for the problem.

### 8.4 Solver setup: Adequate
Serial: CG + Gauss-Seidel. Parallel: CG + HypreILU (or MUMPS).
Tolerances reasonable (rtol=1e-12).

---

## 9. Summary Table

| Item | Status | Severity | Section |
|------|--------|----------|---------|
| BR2 face matrix DOF ordering | **COMPLETED** | ~~CRITICAL~~ | §2.1 |
| IsFaultFace3D bounds check | **COMPLETED** | ~~MEDIUM~~ | §2.2 |
| BR2 penalty geometry-dependent | **COMPLETED** | ~~MEDIUM~~ | §2.3 |
| ComputeTraction (DG0 theoretical) | Not practical | INFO | §2.4 |
| TestNormal code duplication | **COMPLETED** | ~~LOW~~ | §4.1 |
| sigma\_ mutated in AssembleFaceMatrix | **COMPLETED** | ~~LOW~~ | §4.2 |
| mu\_all\_ naming inconsistency | **COMPLETED** | ~~INFO~~ | §4.3 |
| SetupBoundaryMarkers dead code | **COMPLETED** | ~~INFO~~ | §4.4 |
| Dead sigma\_ member (BR2 integrator) | **COMPLETED** | ~~INFO~~ | §4.5 |
| sigma\_ naming in operator | **COMPLETED** | ~~INFO~~ | §4.6 |
| Plan status / line counts | Deviation | INFO | §3.1, §3.3 |
| test\_normal formula | Correct | — | §5.1 |
| BR2 penalty (per element type) | Correct | — | §5.2 |
| Lifting computation | Correct | — | §5.3 |
| Sign conventions | Correct | — | §5.4 |
| Slip BC (IP) | Correct | — | §5.5 |
| Slip BC (BR2) | Correct | — | §5.6 |
| Traction computation | Correct | — | §5.7 |
| Dirichlet loading | Correct | — | §5.8 |
| Boundary condition table | Correct | — | §5.9 |
| FE space setup | Correct | — | §5.10 |
| Makefile | Correct | — | §5.12 |
| Traction non-zero assertion | **COMPLETED** | ~~LOW~~ | §6 |
| Full-system patch test | Suggested | LOW | §6 |
| CF-1: Rotated fault test | Open | LOW | §7 |
| CF-2: Two-sided sign handling | Partially addressed | MEDIUM | §7 |
| CF-3: Higher-order traction | Open | LOW | §7 |
| CF-4: Per-quad-point basis | Open | LOW | §7 |

---

## Remaining items summary

| Category | Count | Details |
|----------|-------|---------|
| Open code bugs | **0** | All resolved |
| Open code quality issues | **0** | All resolved (§4.6 naming fixed) |
| Missing critical tests | **0** | All 6 gaps covered |
| Suggested test improvements | **1 LOW** (§6 full-system patch test) |
| Plan deviations | **2 INFO** (§3.1, §3.3) |
| Carry-forward items | **4** (CF-1 to CF-4) |

**Compared to revision 3**: 2 more items moved to COMPLETED:
- §4.6 (`sigma_` → `epsilon_` rename in operator)
- §6 traction non-zero assertion (now present in both TestNonZeroSlipIP and TestNonZeroSlipBR2)

Phase 2b implementation is fully complete with no open bugs, no open code
quality issues, and comprehensive test coverage (10 BR2 tests + 12 operator
tests). Ready for integration testing with the SEAS time-stepper.
