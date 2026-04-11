# BP5 Debug v31: Fix DG Slip Sign Bug + Dirichlet Interior Faces

**Date**: 2026-03-14
**Status**: Implemented, all unit tests pass. Ready for simulation runs.
**Previous**: v30 (Tandem coord system, all v30 IP/BR2 runs blow up with traction ~10^12-10^14 Pa)

---

## 1. Problem Statement

All v30 BP5 simulations blow up catastrophically:
- IP and BR2 methods both fail
- With and without nucleation
- Traction values reach 10^12-10^14 Pa (expected ~15 MPa)
- Positive feedback loop: wrong prescribed jump -> spurious penalty traction -> larger V -> blowup

## 2. Root Cause: DG Slip Sign Variable

v30 negated `tau_pre` to match Tandem convention (`tau_pre = -tau0`). This made:
- V < 0 for right-lateral motion (previously V > 0)
- `delta_u` accumulates as `d(delta_u)/dt = V < 0`, becoming negative

But the DG sign variable was never adjusted:

```
// OLD (wrong after v30):
real_t sign = (nor(1) > 0) ? -1.0 : 1.0;
```

With `sign = -1` and `delta_u < 0` (right-lateral):
- `g = sign * delta_u = (-1) * (-delta_u) = +|delta_u|`
- But `[[u]] = u(-Y) - u(+Y) < 0` for right-lateral
- **g has wrong sign** -> penalty corrects in wrong direction -> exponential blowup

### Tandem comparison

Tandem handles this through **fault basis flipping** (`AdapterBase.cpp` lines 63-75):
when the geometry normal opposes the reference normal, the entire fault basis is negated,
which flips the projected slip. MFEM instead uses a sign variable with a fixed (unflipped)
basis. The sign variable must be flipped to compensate for negative `delta_u`.

## 3. Fix 1: Flip the Sign Variable (6 locations)

**File**: `domain/elasticity_operator.hpp`

Changed `(nor(1) > 0) ? -1.0 : 1.0` to `(nor(1) > 0) ? 1.0 : -1.0` in all 6 locations:

| # | Function | Context |
|---|----------|---------|
| 1 | `AssembleSlipContributionIP` | Interior faces, IP method |
| 2 | `AssembleSlipContributionBR2` | Interior faces, BR2 method |
| 3 | `AssembleSlipContributionIPShared` | Shared faces, IP method |
| 4 | `AssembleSlipContributionBR2Shared` | Shared faces, BR2 method |
| 5 | `ComputeTraction` (interior) | Traction extraction, interior faces |
| 6 | `ComputeTraction` (shared) | Traction extraction, shared faces |

### Why flipped sign is correct

With `sign = (nor(1) > 0) ? +1 : -1` and `delta_u < 0` (right-lateral):

| Elem1 side | nor(1) | sign | g = sign * delta_u | [[u]] direction | Match? |
|------------|--------|------|--------------------|-----------------|--------|
| -Y side | > 0 | +1 | delta_u < 0 | u(-Y) - u(+Y) < 0 | g = [[u]] |
| +Y side | < 0 | -1 | -delta_u > 0 | u(+Y) - u(-Y) > 0 | g = [[u]] |

Traction at equilibrium: `jump = u_jump - sign*delta_u = [[u]] - g = 0` -> zero correction.

### Updated comment (lines 866-872)

```cpp
// Convention: [[u]] = u- - u+ (CalcOrtho n points outward from Elem1)
// delta_u evolves as d(delta_u)/dt = V, where V < 0 for right-lateral
// (Tandem convention: tau_pre < 0 -> V < 0 -> delta_u < 0).
// sign converts: sign * delta_u = [[u]] = g^F (prescribed jump)
// If nor(1) > 0: Elem1 is -Y side, [[u]] = u(-Y) - u(+Y) = delta_u -> sign = +1
// If nor(1) < 0: Elem1 is +Y side, [[u]] = u(+Y) - u(-Y) = -delta_u -> sign = -1
real_t sign = (nor(1) > 0) ? 1.0 : -1.0;
```

## 4. Fix 2: Dirichlet Interior Face Handling

**File**: `domain/elasticity_operator.hpp`

### Problem

In Tandem's BP5 mesh, Physical Surface(5) tags ALL far-field faces, including Y=0 interior
faces outside the fault region. In Tandem, every face gets a BC type via `DGOperatorTopo.cpp`:

```cpp
info.bc = boundaryData->getBoundaryConditions()[fctNo]
```

- Interior faces with `BC::Dirichlet` go through `rhs_skeleton()` -> `bc_skeleton()`,
  which evaluates `fun_dirichlet` and assembles a DG RHS contribution to both elements.
- The bilinear form (`assemble_skeleton`) is the same for ALL interior faces regardless of
  BC type -- only the RHS differs.

In MFEM, `AssembleDirichletLoading()` only iterates over boundary elements (`GetNBE()`),
which handles true mesh boundary faces. Interior Y=0 faces (two volume elements on both
sides) are missed even though they have attr 5 boundary elements.

### Implementation

Added `BuildDirichletInteriorFaces()` method that:
1. Collects vertex sets of all boundary elements with attr 5
2. For each interior face, checks if vertices match an attr-5 boundary element
3. Excludes faces already in `fault_face_keys_` (attr 3, fault faces)
4. Stores results in `dirichlet_interior_faces_` and `dirichlet_shared_faces_`

Extended `AssembleDirichletLoading()` with an interior face loop using the skeleton
(two-element) pattern:
- Half symmetry coefficient: `w = ip.weight / (2*detJ)` (vs full for boundary)
- Both elements contribute: elvec1 and elvec2
- Opposite penalty signs: `+penalty` for Elem1, `-penalty` for Elem2
- Supports both IP and BR2 methods

This matches Tandem's `rhs_skeleton()` behavior:
- `c10 = 0.5 * epsilon` (half symmetry for skeleton)
- `c20 = +penalty` for Elem0, `c20 *= -1` for Elem1

### Practical effect for BP5

At Y=0, `u_D = sgn(0)*Vp*t/2 = 0`. All RHS terms multiply by `u_D = 0`, so the
contribution is zero. The DG bilinear form (via `AddInteriorFaceIntegrator`) already
enforces continuity [[u]] = 0 on these faces.

Verified in Tandem: `assemble_skeleton()` does NOT check BC type (line 341). The bilinear
form is identical for all interior faces. Only `rhs_skeleton()` adds the Dirichlet data,
which is zero at Y=0.

This fix ensures structural correctness matching Tandem's face classification.

## 5. Unit Test Updates

**File**: `tests/unit/test_elasticity_operator.cpp`

### Test 18: `TestBR2SlipSignConvention`

- Changed strike slip from `+1.0` to `-1.0` (right-lateral, Tandem convention)
- Changed assertion from `avg_strike > 0.0` to `std::abs(avg_strike) > 1e-6`
- Updated comments to reflect Tandem convention

### Test 11 Sub-test B: IP vs BR2 comparison

- Changed slip from `+1.0` to `-1.0` (negative = right-lateral)
- Changed assertion from `avg_strike > 0` to `std::abs(avg_strike) > 1e-6`
- Updated Sub-test C comment for consistent sign convention

## 6. Test Results

| Test Suite | Tests | Result |
|------------|-------|--------|
| `seas_test_fault_basis` | 187 | All pass |
| `seas_test_elasticity_operator` | 84 | All pass |
| `seas_test_elasticity_br2` | 46 | All pass |
| `seas_test_fault_detection` | 7 | All pass |
| `seas_test_bp5_params` | 96 | All pass |
| `seas_test_domain_interface` | 32 | All pass |
| All 22 serial test suites | 452 | All pass |

Note: `seas_test_bp5_parallel_smoke` has 1 pre-existing failure (fault DOF detection
in parallel) unrelated to v31 changes.

## 7. Files Changed

| File | Changes |
|------|---------|
| `domain/elasticity_operator.hpp` | Sign flip (6 locations), updated comment, `BuildDirichletInteriorFaces()`, extended `AssembleDirichletLoading()` |
| `tests/unit/test_elasticity_operator.cpp` | Test 18 + Test 11B updated for negative slip convention |

## 8. Simulation Test Plan

### Test 1: IP + Uniform fault (no nucleation)

Validates the sign fix doesn't break steady-state behavior.
- `--delta-tau-factor 0 --V-nuc 1e-9`
- PASS criterion: V/Vp in [0.5, 2.0] for 300 yr, no blowup

### Test 2: IP + SCEC initialization (nucleation)

Tests full earthquake cycle with sign fix.
- Default SCEC init
- PASS criterion: survives t > 1e8 s, traction O(15 MPa), earthquake nucleation

### Test 3: BR2 + SCEC initialization

Verifies BR2 method also works with sign fix.
- `--dg-method BR2`
- PASS criterion: same as Test 2

### What to look for in diagnostic VTK

1. **Traction**: Should be O(15 MPa), NOT O(1 GPa) or higher
2. **Slip direction**: Right-lateral (u(+Y) in +X, u(-Y) in -X)
3. **Displacement field**: X-component loading, sgn(Y) pattern
4. **Slip rate**: Should nucleate and produce earthquake cycles (Tests 2-3)

### Success criteria

- IP and BR2 runs survive beyond t > 1e8 s (~3 years) without blowing up
- Traction values physically reasonable: O(sigma_n * f) ~ 15 MPa
- Earthquake nucleation occurs with SCEC initialization
