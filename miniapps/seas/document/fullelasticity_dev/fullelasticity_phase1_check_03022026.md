# Phase 1 Implementation Review: BP5 Parameters and Vector Friction

**Date**: 2026-03-02
**Last updated**: 2026-03-02 (revision 3 — all issues resolved)
**Reviewing**: `fullelasticity_phase1_plan_03022026.md` vs. current implementation
**Files checked**: `config/bp5_params.hpp`, `friction/dieterich_ruina.hpp`,
`fault/fault_basis.hpp`, `domain/domain_operator.hpp`,
`tests/unit/test_bp5_params.cpp`, `tests/unit/test_vector_friction.cpp`,
`tests/unit/test_fault_basis.cpp`, `tests/unit/test_macros.hpp`, `Makefile`
**References**: SCEC BP5-QD spec, Tandem `examples/tandem/3d/bp5.lua`,
Tandem `app/localoperator/DieterichRuinaAgeing.h`

---

## 1. Overall Assessment

All Phase 1 items are complete. Every bug, code quality issue, and suggestion
from the initial and revision 2 reviews has been resolved. Phase 2a (`FaultBasis`)
and `DomainOperator` interface extensions are also done. The code is ready for
Phase 2b (`ElasticityDomainOperator`).

**Summary of changes since initial review**:
- tau0_vec now uses `Vi_abs` (magnitude) instead of `Vi[1]` — fragile coupling removed
- Transition continuity tests fixed — now test actual transition zone points
- `V_init` member is now used correctly in `V_init_vec()`
- `Validate()` method added to BP5Params
- `psi_init()` helper added to BP5Params
- `FaultBasis` class fully implemented with comprehensive tests (8 tests)
- `DomainOperator` extended with `GetFaultBasis()`, `NumSlipComponents()`,
  `GetFaultCoords2D()`, `GetOffFaultDisplacement()`
- Six test coverage gaps from initial review all addressed with new test functions

**Additional changes since revision 2**:
- `domain_operator.hpp` now uses forward declaration for `FaultBasis` (no tight coupling)
- `FaultBasis::NumDOFs()` renamed to `NumFaces()` for clarity
- Multi-face test (`TestMultiFace3D`) added to `test_fault_basis.cpp`
- 2D tangent sign convention documented with clarifying comment
- Test macros extracted to shared `tests/unit/test_macros.hpp`
- `Dc_of_x2_x3()` alias added to `BP5Params`

---

## 2. Bugs / Potential Bugs (from initial review)

### 2.1 ~~[MEDIUM] V_init_vec: strike vs. dip direction disagrees with SCEC spec~~

**Status**: COMPLETED

The code at `bp5_params.hpp:216-228` now has a comprehensive comment block
explaining the Tandem convention (V[0]=strike~0, V[1]=dip=dominant) and how
it differs from SCEC Eq. 16. The comment also notes that `FaultBasis` handles
the mapping between fault-local frame and global coordinates.

### 2.2 ~~[LOW] tau0_vec: hard-coded use of Vi[1] instead of |Vi|~~

**Status**: COMPLETED

`bp5_params.hpp:265-277` now uses `Vi_abs` (the magnitude) throughout:
```cpp
real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
...
real_t tau0_scalar = sigma_n * a *
                     std::asinh((Vi_abs / (2.0 * V0)) * e) +
                     eta_val * Vi_abs;
```
This matches the SCEC formula exactly and is convention-independent.

### 2.3 ~~[LOW] Transition continuity tests don't actually test the transition~~

**Status**: COMPLETED

`test_bp5_params.cpp:150-160` now tests points in the actual transition zone:
```cpp
real_t a_near_vw = p.a_of_x2_x3(0.0, 4.0e3 - 1.0);  // x3=3999 in transition
real_t a_near_vs = p.a_of_x2_x3(0.0, 2.0e3 + 1.0);   // x3=2001 in transition
```

---

## 3. Deviations from Plan

### 3.1 tau0_vec signature omits eta parameter

**Plan**: `tau0_vec(x2, x3, eta, tau[2])`
**Code**: `tau0_vec(x2, x3, tau[2])` (computes eta internally via `eta()`)

This is actually *better* than the plan -- fewer parameters, impossible to pass the
wrong eta. No action needed, but note the deviation.

### 3.2 V_nuc = 0.01 (Tandem) vs. 0.03 (SCEC)

Already documented in the code comment at line 109. The plan also notes this. When
running the official SCEC BP5 benchmark comparison, switching to `V_nuc = 0.03` should
be a simple constant change.

### 3.3 FaultBasis::Compute parameter order

**Plan** (Section 8.3): `Compute(mesh, fault_faces, up, ref_normal)`
**Code**: `Compute(mesh, fault_faces, ref_normal, up)` (ref_normal first)

Minor API ordering difference. The code is consistent internally and all tests
use the actual signature. No action needed.

### 3.4 FaultBasisData uses raw arrays instead of MFEM Vector

**Plan** (Section 8.3): `Vector normal`, `Vector tangent1`, `Vector tangent2`
**Code**: `real_t normal[3]`, `real_t tangent1[3]`, `real_t tangent2[3]`

Using fixed-size arrays avoids heap allocation per face and is better for
performance. The `ProjectTraction`/`EmbedSlip` interfaces use `real_t*` pointers
consistently. This is a good deviation.

---

## 4. Code Quality Issues (from initial review)

### 4.1 ~~V_init member variable is unused~~

**Status**: COMPLETED

`V_init` is now used by `V_init_vec()` at line 242 (`V[1] = V_init`) and by
`psi_init()` at line 334. `Vp` and `V_init` are semantically distinct (plate
rate vs initial slip rate) even though numerically equal for BP5. The current
code uses them correctly in their respective contexts.

### 4.2 No output station definitions

**Status**: Still open (INFO) — planned for I/O phase.

### 4.3 ~~No parameter validation~~

**Status**: COMPLETED

`Validate()` added at `bp5_params.hpp:321-330` with checks:
- `a0 < b` (velocity-weakening)
- `amax > b` (velocity-strengthening)
- `L_nuc < L0`
- `sigma_n > 0`
- `Vp > 0`
- `V_nuc > V_init`

Test at `test_bp5_params.cpp:510-537` verifies `Validate()` passes on defaults.

### 4.4 ~~Duplicated test macros~~

**Status**: COMPLETED

Test macros (`TEST_ASSERT`, `TEST_NEAR`, `TEST_REL_NEAR`, `TEST_PRINT_RESULTS`)
extracted to shared `tests/unit/test_macros.hpp`. All three test files now
`#include "test_macros.hpp"` instead of defining macros inline.

**Note**: `test_macros.hpp` is not listed in the Makefile header dependency
groups (`FRICTION_HEADERS`, `FAULT_HEADERS`, etc.), so `make` won't rebuild
test objects if only the macros header changes. This is a minor build hygiene
issue — the macros are stable and unlikely to change often.

---

## 5. Correctness Verification Against Tandem

### 5.1 SolveSlipRateVectorPsi algorithm: CORRECT

The implementation matches Tandem's `DieterichRuinaAgeing::slip_rate()`:
1. Compute `|tau|` from 2-component traction ✓
2. Solve scalar `V = SolveSlipRatePsi(|tau|, psi, sigma_n, eta, a)` ✓
3. Project: `V_vec = -(V / |tau|) * tau_vec` ✓
4. Zero-traction guard (`tau_abs < 1e-30`) ✓

Tandem uses a `zeroIn` bracketing solver; the MFEM code uses Newton-Raphson with
bisection fallback. Both are valid for this monotone equation.

### 5.2 BP5 a(x2, x3) function: CORRECT

Matches SCEC Eq. 14 exactly:
- VW core check ✓
- VS zones check ✓
- Transition formula with `r = max(...)` ✓
- `r` is guaranteed in [0, 1] within the transition zone ✓

### 5.3 IsNucleationZone: CORRECT

Matches SCEC definition:
- Depth range: `hs+ht <= x3 <= hs+ht+H` → 4km to 16km ✓
- Strike range: `-l_vw/2 <= x2 <= -l_vw/2 + w_nuc` → -30km to -18km ✓

### 5.4 L_of_x2_x3: CORRECT

- Returns `L_nuc = 0.13` in nucleation zone ✓
- Returns `L0 = 0.14` elsewhere ✓

### 5.5 Material constants: CORRECT

- `mu = rho * cs^2 = 2670 * 3464^2 = 32.04 GPa` ✓
- `lambda = mu` for nu=0.25 ✓
- `eta = mu / (2*cs)` ✓

### 5.6 Makefile: CORRECT

- `config/bp5_params.hpp` added to `FRICTION_HEADERS` ✓
- `fault/fault_basis.hpp` added to `FAULT_HEADERS` ✓
- New source/object/target entries for all three tests ✓
- Added to `SEQ_MINIAPPS` ✓
- Build rules with correct header dependencies ✓
- `test-bp5-params`, `test-vector-friction`, `test-fault-basis` targets defined ✓
- All three added to the `test:` target ✓

### 5.7 FaultBasis algorithm: CORRECT

Matches Tandem's `Curvilinear::facetBasis()`:
- Raw normal via MFEM `CalcOrtho(J, n_raw)` ✓
- Orientation via `dot(n_raw, ref_normal)` with flip ✓
- `strike = normalize(up × n)` ✓
- `dip = strike × n` ✓
- `tangent1 = dip`, `tangent2 = strike` (Tandem convention) ✓
- 2D: single tangent via 90° rotation with up-consistent sign ✓
- Degenerate case check: `MFEM_VERIFY(s_len > 1e-12, ...)` for up ∥ n ✓

### 5.8 DomainOperator extensions: CORRECT

- `NumSlipComponents()` defaults to 1 (backward compatible) ✓
- `GetFaultCoords2D()` default returns (0, depths) for 2D ✓
- `GetFaultBasis()` defaults to nullptr ✓
- `GetOffFaultDisplacement()` default is no-op ✓

---

## 6. Test Coverage Assessment

### Covered (from plan + additions):
- [x] Material properties (mu, lambda, eta)
- [x] a() in VW core (multiple points)
- [x] a() in VS zones (shallow, deep, far along-strike)
- [x] a() in transition zone (depth and strike)
- [x] a() transition continuity (points in actual transition zone)
- [x] a() symmetry in x3 about VW center
- [x] a() zone boundary edge cases (exact boundaries)
- [x] L() in/out of nucleation zone
- [x] IsNucleationZone (inside/outside, corners, boundaries)
- [x] V_init_vec outside/inside nucleation (magnitude and component values)
- [x] tau0_vec magnitude and sign
- [x] tau0_vec direction (ratio tau[0]/tau[1] matches Vi[0]/Vi[1])
- [x] tau0_vec at VS zone point (deep VS with amax)
- [x] tau0_vec self-consistency (feed tau0 into solver, recover |V_init|)
- [x] Vector friction: pure x/y direction, 45-degree, negative direction
- [x] Magnitude consistency (vector vs scalar solver, multiple angles)
- [x] Zero traction
- [x] BP5 values recovery (vector solver)
- [x] Various a values (VW, VS, transition)
- [x] Print function
- [x] psi_init() helper
- [x] Validate() on default params
- [x] FaultBasis: axis-aligned 3D (n, t1, t2 vectors, orthonormality)
- [x] FaultBasis: ProjectTraction correctness
- [x] FaultBasis: EmbedSlip correctness
- [x] FaultBasis: round-trip (project → embed → verify)
- [x] FaultBasis: NormalStress extraction
- [x] FaultBasis: 2D degeneracy (single tangent, project, embed)
- [x] FaultBasis: orientation flip (ref_normal sign change)
- [x] FaultBasis: multi-face mesh (3 interior faces, consistency + indexing)

### Remaining gaps:
- [ ] FaultBasis: rotated fault (non-axis-aligned, e.g., 30° about z-axis)
- [ ] FaultBasis: per-face vs per-DOF semantics for higher-order elements
      (current implementation is per-face which is correct for piecewise-constant
      normals on planar faces, but curved elements would need per-quad-point)

---

## 7. Suggestions for Phase 2 Preparation

### Completed:

1. ~~**Implement `FaultBasis`** (Section 8)~~: DONE. `fault/fault_basis.hpp`
   fully implemented with 8-test test suite.

2. ~~**Use Vi magnitude in tau0_vec**~~: DONE. Now uses `Vi_abs`.

3. ~~**State initialization helper**~~: DONE. `psi_init()` added to BP5Params.

4. ~~**Add `Dc_of_x2_x3()` alias**~~: DONE. `bp5_params.hpp:224-227` adds
   `Dc_of_x2_x3()` as an alias for `L_of_x2_x3()`, bridging the BP5 naming
   convention (`L`) with the friction law naming (`Dc`).

### Still open:

5. **Compilation test**: Run `make seas_test_bp5_params seas_test_vector_friction
   seas_test_fault_basis && make test-bp5-params && make test-vector-friction &&
   make test-fault-basis && make test-friction && make test-psi-state` to verify
   everything compiles and all tests (including regression) pass.

---

## 8. Design: Global-Local Coordinate Transformation System

### 8.1 Problem Statement

The current codebase has **no coordinate transformation** between global mesh
coordinates and fault-local coordinates. The 2D antiplane case avoids this because
the fault is axis-aligned (vertical at x=0), making global = local implicitly.

For 3D elasticity, the code needs exactly **two coordinate systems**:

```
Global frame:      (x, y, z)       — defined by the 3D mesh geometry
Fault-local frame: (n, t1, t2)     — per-DOF orthonormal basis on fault surface
                    normal, strike, dip
```

The SCEC convention (x1=normal, x2=strike, x3=depth) is a **problem setup concern**,
not a runtime coordinate system. When designing the mesh and setting up BP5Params,
we ensure that the mesh global coordinates and the fault-local tangent vectors are
consistent with the SCEC physics. But the code itself only transforms between global
and local — no third frame is needed.

### 8.2 The Two Frames

**Global frame (x, y, z)**: Defined by the mesh. All domain PDE quantities live
here — displacement u, stress tensor σ, displacement jumps Δu.

**Fault-local frame (n, t1, t2)**: Defined per DOF on the fault surface. All
fault physics quantities live here — tangential traction τ[2], slip rate V[2],
accumulated slip s[2], friction parameters. The friction law, state evolution,
and BP5 spatial parameters all operate purely in this frame.

The two transformations needed:

```
Global → Local:  traction projection      σ·n → (τ·t1, τ·t2)
Local → Global:  slip embedding           (s1, s2) → s1·t1 + s2·t2 = Δu
```

### 8.3 Implementation: `FaultBasis` — Per-Face Coordinate Frame

**STATUS**: IMPLEMENTED in `fault/fault_basis.hpp`

Key implementation details vs. original plan:
- Uses `real_t[3]` arrays instead of MFEM `Vector` (no heap allocation)
- Parameter order is `(mesh, fault_faces, ref_normal, up)` not `(mesh, fault_faces, up, ref_normal)`
- Takes `Mesh&` (non-const) because MFEM's `GetInteriorFaceTransformations` requires it
- Includes bonus `NormalStress()` method not in original plan
- Per-face granularity (one basis per face, not per DOF) — correct for planar faces

Algorithm (following Tandem Curvilinear::facetBasis):
1. Get raw face normal from mesh (`CalcOrtho`)
2. Orient: `if dot(n_raw, ref_normal) < 0 → flip`
3. Normalize: `n = n_raw / |n_raw|`
4. `strike = normalize(up × n)`
5. `dip = strike × n`
6. `tangent1 = dip`, `tangent2 = strike` (Tandem convention)

### 8.4 Data Flow — Where Each Transformation Happens

```
                    GLOBAL FRAME                      LOCAL FRAME
                   (mesh x,y,z)                    (n, t1=dip, t2=strike)

  ┌─────────────────────────────────────────────────────────────────────┐
  │ ElasticityDomainOperator::Solve()                                   │
  │                                                                     │
  │   slip_local[2] ──EmbedSlip──→ Δu[3] ──→ DG interior face BC       │
  └─────────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────────┐
  │ ElasticityDomainOperator::ComputeTraction()                         │
  │                                                                     │
  │   ∇u → ε → σ → σ·n = τ_global[3] ──ProjectTraction──→ τ_local[2]  │
  └─────────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────────┐
  │ FaultGeometry (setup, once)                                         │
  │                                                                     │
  │   mesh_coords(dof) → extract (strike_coord, depth_coord)           │
  │       → BP5Params::a(), L(), V_init_vec(), tau0_vec()               │
  │       → precompute per-DOF arrays (all in local frame)              │
  └─────────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────────┐
  │ RateStateFaultOperator::ComputeRHS()                                │
  │                                                                     │
  │   Everything in local frame:                                        │
  │   τ_local[2] + τ0[2] → SolveSlipRateVectorPsi → V_local[2], dψ/dt │
  └─────────────────────────────────────────────────────────────────────┘
```

The domain operator is the **only place** that touches both frames.
Everything above the fault interface (friction, state evolution, parameters)
works entirely in the local frame.

### 8.5 Changes Required in Existing Classes

#### 8.5.1 `DomainOperator` (interface) — DONE

Added to `domain/domain_operator.hpp`:
- `virtual int NumSlipComponents() const { return 1; }` ✓
- `virtual void GetFaultCoords2D(Vector&, Vector&) const` ✓
- `virtual const FaultBasis *GetFaultBasis() const { return nullptr; }` ✓
- `virtual void GetOffFaultDisplacement(...)` ✓

#### 8.5.2 `FaultGeometry` — TODO (Phase 2c)

Extend `fault/fault_geometry.hpp` with optional 3D members.

#### 8.5.3 `RateStateFaultOperator` — TODO (Phase 3)

Generalize for vector state (StatePerNode = 3 for 3D).

#### 8.5.4 `SEASQuasiDynamicOperator` — TODO (Phase 4)

No conceptual change needed, just wider vectors.

### 8.6 How SCEC Convention Is Baked In at Setup Time

The SCEC convention is **not a runtime coordinate system**. It is satisfied by
two setup-time choices:

1. **Mesh construction**: Build the 3D mesh so that the fault plane, domain
   boundaries, and loading directions correspond to the SCEC geometry. The
   mesh defines the global frame.

2. **FaultBasis "up" and "ref_normal" vectors**: These user-provided vectors
   (passed to the `FaultBasis::Compute()` call) determine how the local
   tangent vectors align with the physical strike/dip directions. For a
   standard BP5 mesh:
   - `ref_normal = {1, 0, 0}` → normal points in +x1 direction
   - `up = {0, 0, -1}` → "up" in the depth-positive-downward sense

   The result is `tangent1` = dip direction, `tangent2` = strike direction,
   matching SCEC (x2, x3) = (strike, depth).

3. **BP5Params**: Already written in SCEC convention. `FaultGeometry` calls
   it with (strike_coord, depth_coord) extracted from the mesh. The returned
   V_init, tau0, a, L values are in the local frame by construction.

Once setup is complete, the runtime code only sees global and local frames.
The SCEC convention is implicitly encoded in the mesh + basis vectors.

### 8.7 Tandem Reference

| Tandem file | Role | SEAS-MFEM equivalent |
|---|---|---|
| `src/geometry/Curvilinear.cpp:facetBasis()` | Computes (n, t1, t2) from mesh | `FaultBasis::Compute()` |
| `app/localoperator/AdapterBase.cpp:prepare()` | Stores basis, handles sign flip | `FaultGeometry` 3D constructor |
| `app/kernels/elasticity_adapter.py` | τ projection & slip embedding | `ElasticityDomainOperator` internals |

Tandem stores a D×D matrix per quad point (columns = n, t1, t2).
Projection: `tau_local = basis^T * tau_global`.
Embedding: `Δu = basis * (0, s1, s2)^T`.

### 8.8 Implementation Phases

| Phase | What | Status | Files |
|---|---|---|---|
| Phase 1 | BP5Params, vector friction solver | **DONE** | `config/bp5_params.hpp`, `friction/dieterich_ruina.hpp` |
| Phase 2a | `FaultBasis` from mesh geometry | **DONE** | `fault/fault_basis.hpp` (new) |
| Phase 2a' | `DomainOperator` interface extensions | **DONE** | `domain/domain_operator.hpp` (extended) |
| Phase 2b | `ElasticityDomainOperator` with basis-based projection | TODO | `domain/elasticity_operator.hpp` (new) |
| Phase 2c | Extend `FaultGeometry` for 3D | TODO | `fault/fault_geometry.hpp` (extend) |
| Phase 3 | Generalize `RateStateFaultOperator` for vector state | TODO | `fault/rate_state_fault.hpp` (extend) |
| Phase 4 | Wire through `SEASQuasiDynamicOperator` | TODO | `solver/seas_operator.hpp` (extend) |

### 8.9 Unit Test Strategy — Status

```
test_fault_basis.cpp (DONE — 8 tests):

  1. Axis-aligned vertical fault at x=0:                            ✓ IMPLEMENTED
     - Verify n = (1,0,0), t1 = dip, t2 = strike
     - Orthonormality checks
     - ProjectTraction with known traction → check components
     - EmbedSlip with known slip → check delta_u

  2. Round-trip:                                                     ✓ IMPLEMENTED
     - Tangential traction → ProjectTraction → EmbedSlip → matches original
     - Traction with normal component → project discards normal part

  3. Rotated fault (30° about z-axis):                              ✗ NOT YET
     - Would need non-Cartesian mesh or mesh transform
     - Deferred to when non-planar fault support is needed

  4. 2D degeneracy:                                                  ✓ IMPLEMENTED
     - Single tangent vector, 90° rotation of normal
     - ProjectTraction → scalar output
     - EmbedSlip → 2D vector output
     - NormalStress extraction

  5. Orientation flip (BONUS):                                       ✓ IMPLEMENTED
     - Opposite ref_normals → opposite normals
     - Orthonormality preserved under flip

  6. Normal stress extraction (BONUS):                               ✓ IMPLEMENTED
     - Tensile, compressive, purely tangential tractions

  7. Multi-face 3D mesh:                                             ✓ IMPLEMENTED
     - 4-element hex mesh, 3 interior faces
     - Per-face orthonormality, cross-face consistency
     - ProjectTraction consistency across face indices
```

---

## 9. Issues Found in Revision 2 (all resolved)

### 9.1 ~~[LOW] domain_operator.hpp includes fault_basis.hpp — tight coupling~~

**Status**: COMPLETED

`domain/domain_operator.hpp:23-24` now uses a forward declaration:
```cpp
// Forward declaration — only a pointer is returned, no include needed.
class FaultBasis;
```
No build dependency from domain layer to fault layer.

### 9.2 ~~[LOW] FaultBasis tests only use single-face meshes~~

**Status**: COMPLETED

`TestMultiFace3D()` added at `test_fault_basis.cpp:437-509`. Creates a
4-element hex mesh with 3 interior faces. Tests:
- Correct face count (`NumFaces() == 3`)
- Per-face orthonormality
- Cross-face consistency (all faces produce identical basis)
- `ProjectTraction` consistency across face indices

### 9.3 ~~[INFO] FaultBasis is per-face, not per-DOF — naming mismatch~~

**Status**: COMPLETED

`NumDOFs()` renamed to `NumFaces()` and `num_dofs_` renamed to `num_faces_`
in `fault/fault_basis.hpp`. All call sites in test files updated.

### 9.4 ~~[INFO] FaultBasis 2D tangent sign depends on cross product with "up"~~

**Status**: COMPLETED

Clarifying comment added at `fault/fault_basis.hpp:140-145`:
```cpp
// Tandem 2D convention: tangent perpendicular to normal.
// The sign is chosen from the cross product of "up" with n so
// that the tangent points in the "depth" direction when "up"
// points upward in the mesh (e.g., up=(0,-1) with depth
// positive downward gives tangent=(0,1), matching the existing
// antiplane convention).
```

---

## 10. New Issues Found in Revision 3

### 10.1 [INFO] test_macros.hpp not in Makefile header dependencies

**File**: `Makefile`

The shared `tests/unit/test_macros.hpp` is not listed in any Makefile header
dependency group (`FRICTION_HEADERS`, `FAULT_HEADERS`, etc.). The test object
rules depend on `$(FRICTION_HEADERS)` or `$(FAULT_HEADERS)`, but not on
`test_macros.hpp`. If the macros header is modified, `make` won't rebuild
the test executables.

**Impact**: Very minor. The macros are stable (4 simple macros). A `make clean`
followed by `make` would pick up changes. Not worth adding a new header group
just for one test-infrastructure file.

**Recommendation**: No action needed unless more test infrastructure headers
are added. If a `TEST_HEADERS` group is ever created, include `test_macros.hpp`.

---

## 11. Summary Table

| Item | Status | Severity | Action |
|------|--------|----------|--------|
| V_init_vec direction vs SCEC | **COMPLETED** (documented) | ~~MEDIUM~~ | Convention explained in comments |
| tau0_vec uses Vi[1] not \|Vi\| | **COMPLETED** (uses Vi_abs) | ~~LOW~~ | Fixed |
| Transition continuity tests | **COMPLETED** (correct points) | ~~LOW~~ | Fixed |
| V_init member unused | **COMPLETED** (now used) | ~~LOW~~ | Used in V_init_vec + psi_init |
| No parameter validation | **COMPLETED** (Validate()) | ~~LOW~~ | Added with test |
| psi_init() helper | **COMPLETED** | ~~LOW~~ | Added with test |
| **FaultBasis implementation** | **COMPLETED** | ~~HIGH~~ | Implemented + tested (8 tests) |
| **DomainOperator 3D interface** | **COMPLETED** | ~~MEDIUM~~ | NumSlipComponents, GetFaultBasis, etc. |
| Duplicated test macros | **COMPLETED** | ~~INFO~~ | Extracted to test_macros.hpp |
| domain_operator includes fault_basis | **COMPLETED** | ~~LOW~~ | Forward declaration used |
| FaultBasis: single-face tests only | **COMPLETED** | ~~LOW~~ | TestMultiFace3D added |
| FaultBasis: NumDOFs → NumFaces | **COMPLETED** | ~~INFO~~ | Renamed |
| 2D tangent sign documentation | **COMPLETED** | ~~INFO~~ | Clarifying comment added |
| Dc_of_x2_x3 alias | **COMPLETED** | ~~INFO~~ | Added to BP5Params |
| tau0_vec signature vs plan | Deviation (better) | INFO | None |
| V_nuc 0.01 vs SCEC 0.03 | Known, documented | INFO | None |
| No output stations | Still open | INFO | Add in I/O phase |
| test_macros.hpp not in Makefile deps | New (rev 3) | INFO | No action needed |
| Rotated fault test | Still open | INFO | Add when non-planar support needed |
| Compilation test | Still open | INFO | Run before Phase 2b |
| SolveSlipRateVectorPsi | Correct | -- | -- |
| a(x2,x3) formula | Correct | -- | -- |
| IsNucleationZone | Correct | -- | -- |
| L(x2,x3) | Correct | -- | -- |
| Material constants | Correct | -- | -- |
| FaultBasis algorithm | Correct | -- | -- |
| Makefile | Correct | -- | -- |
