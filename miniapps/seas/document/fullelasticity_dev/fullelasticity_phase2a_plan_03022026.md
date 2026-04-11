# Phase 2a: FaultBasis and DomainOperator Interface Extensions

**Date**: 2026-03-02
**Status**: COMPLETE (13 FaultBasis test functions with 187 assertions + 5 DomainOperator interface test functions with 32 assertions)
**Prerequisite**: Phase 1 core complete (BP5Params, vector friction)

---

## Overview

Phase 2a implements the coordinate transformation system needed for 3D fault surfaces:

1. **`FaultBasis`** — Per-face orthonormal coordinate frames for projecting between global and fault-local frames
2. **`DomainOperator` interface extensions** — Virtual methods enabling 3D domain operators

These are prerequisite infrastructure for the `ElasticityDomainOperator` (Phase 2b).

---

## Part 1: FaultBasis

**File**: `fault/fault_basis.hpp` (NEW)

### Problem Statement

The 2D antiplane case avoids coordinate transforms because the fault is axis-aligned (vertical at x=0), making global = local implicitly. For 3D elasticity, we need two transforms:

```
Global → Local:  traction projection      σ·n → (τ·t1, τ·t2)
Local → Global:  slip embedding           (s1, s2) → s1·t1 + s2·t2 = Δu
```

### Algorithm (Following Tandem's `Curvilinear::facetBasis()`)

For each interior fault face:

1. Get raw face normal from mesh: `CalcOrtho(Jacobian, n_raw)`
2. Orient: `if dot(n_raw, ref_normal) < 0 → flip`
3. Normalize: `n = n_raw / |n_raw|`
4. Compute strike direction: `strike = normalize(up × n)`
5. Compute dip direction: `dip = strike × n`
6. Store: `tangent1 = dip`, `tangent2 = strike` (Tandem convention)

For 2D (single tangent): 90° rotation of normal with sign from `up`.

### Data Structure

```cpp
struct FaultBasisData {
   real_t normal[3];    // Unit normal (oriented by ref_normal)
   real_t tangent1[3];  // Dip direction
   real_t tangent2[3];  // Strike direction
};
```

Uses fixed-size `real_t[3]` arrays (no heap allocation per face).

### Class Interface

```cpp
class FaultBasis {
public:
   void Compute(Mesh &mesh, const Array<int> &fault_faces,
                const Vector &ref_normal, const Vector &up);

   int NumFaces() const;
   int Dimension() const;
   int NumTangentComponents() const;  // 1 for 2D, 2 for 3D

   const FaultBasisData &GetBasis(int fi) const;

   // Coordinate transforms
   void ProjectTraction(int fi, const real_t *T_global, real_t *tau_local) const;
   void EmbedSlip(int fi, const real_t *slip_local, real_t *delta_u) const;
   real_t NormalStress(int fi, const real_t *T_global) const;
};
```

### ProjectTraction

For 3D: `tau_local = [T_global · tangent1, T_global · tangent2]` (2 components)
For 2D: `tau_local = [T_global · tangent1]` (1 component)

### EmbedSlip

For 3D: `delta_u = slip_local[0] * tangent1 + slip_local[1] * tangent2` (no-opening constraint: normal component = 0)
For 2D: `delta_u = slip_local[0] * tangent1`

### Data Flow in the Full System

```
ElasticityDomainOperator::Solve()
  slip_local[2] ──EmbedSlip──→ Δu[3] ──→ DG interior face BC

ElasticityDomainOperator::ComputeTraction()
  ∇u → ε → σ → σ·n = T_global[3] ──ProjectTraction──→ τ_local[2]
```

The domain operator is the **only place** touching both frames. Everything above (friction, state evolution, parameters) works entirely in the local frame.

---

## Part 2: DomainOperator Interface Extensions

**File**: `domain/domain_operator.hpp` (EXTENDED)

### New Virtual Methods (Phase 2a Additions Only)

Note: `NumComponents()` and `Dimension()` are **pre-existing pure virtual** methods
(not Phase 2a additions). They remain `= 0` and are already overridden by
`AntiplaneDomainOperator`. Only the methods below were added in Phase 2a:

```cpp
template <typename MeshType = Mesh>
class DomainOperator {
public:
   // Pre-existing pure virtual methods (unchanged)
   virtual int NumComponents() const = 0;    // 1 antiplane, 3 for 3D
   virtual int Dimension() const = 0;         // 2 for 2D, 3 for 3D
   virtual void Solve(...) = 0;
   virtual void ComputeTraction(...) = 0;
   virtual int GetNumFaultDOFs() const = 0;
   virtual void GetFaultDepths(Vector &depths) const = 0;

   // NEW Phase 2a methods (with backward-compatible defaults)
   virtual int NumSlipComponents() const { return 1; }       // 2 for 3D
   virtual const FaultBasis *GetFaultBasis() const { return nullptr; }
   virtual void GetFaultCoords2D(Vector &x2, Vector &x3) const {
      // Default: x2 = 0, x3 = depths (backward compatible for 2D)
   }
   virtual void GetOffFaultDisplacement(const std::vector<Vector> &points,
                                         Vector &displacements) const {}
};
```

### Key Design

- **Forward declaration** of `FaultBasis` (not `#include`) to avoid coupling domain layer to fault layer
- **Default implementations** return backward-compatible values (1 slip component, nullptr basis, etc.)
- **No changes needed** in existing `AntiplaneDomainOperator` — inherits defaults

---

## Tests

### FaultBasis Tests (`tests/unit/test_fault_basis.cpp`) — 13 test functions, 187 assertions

1. **Axis-aligned 3D**: n=(1,0,0), t1=dip, t2=strike, orthonormality
2. **ProjectTraction**: Known global traction → correct local components
3. **EmbedSlip**: Known local slip → correct global displacement jump
4. **Round-trip**: project → embed → verify consistency
5. **NormalStress**: Tensile, compressive, purely tangential tractions
6. **2D degeneracy**: Single tangent, project/embed in 2D
7. **Orientation flip (3D)**: Opposite `ref_normal` → opposite normals, preserved orthonormality
8. **Multi-face mesh (3D)**: 4-element hex mesh, 3 interior faces, cross-face consistency
9. **Empty fault_faces**: Graceful handling of zero-face input
10. **Multi-face mesh (2D)**: 4-element quad mesh, 3 interior faces, 2D consistency
11. **Orientation flip (2D)**: Opposite `ref_normal` in 2D, orthonormality preserved
12. **Non-unit input vectors**: Scaled ref_normal and up produce same basis as unit inputs (3D + 2D)
13. **Near-collinear up/n**: Nearly-parallel up/normal still produces orthonormal basis (3D + 2D)

### DomainOperator Interface Tests (`tests/unit/test_domain_operator_interface.cpp`) — 5 test functions, 32 assertions

1. NumSlipComponents() returns 1 for antiplane
2. GetFaultCoords2D() returns (0, depths) by default
3. GetOffFaultDisplacement() returns empty by default
4. Polymorphic dispatch through base pointer (includes GetFaultBasis nullptr check)
5. BdrLoad operator inherits defaults (includes GetFaultBasis nullptr check)

---

## SCEC Convention Encoding

The SCEC convention (x1=normal, x2=strike, x3=depth) is **not** a runtime coordinate system. It is satisfied by two setup-time choices:

1. **Mesh construction**: Fault plane at x1=0, domain aligned with SCEC axes
2. **FaultBasis parameters**:
   - `ref_normal = {1, 0, 0}` → normal in +x1 direction
   - `up = {0, 0, -1}` → depth-positive-downward

Result: `tangent1`=dip, `tangent2`=strike, matching SCEC (x3, x2).

---

## Files Created/Modified

| File | Action |
|------|--------|
| `fault/fault_basis.hpp` | NEW — FaultBasis class |
| `domain/domain_operator.hpp` | EXTENDED — new virtual methods |
| `tests/unit/test_fault_basis.cpp` | NEW — 13 test functions, 187 assertions |
| `Makefile` | EXTENDED — new test target |

## Tandem Reference

| SEAS-MFEM | Tandem |
|-----------|--------|
| `FaultBasis::Compute()` | `Curvilinear::facetBasis()` |
| `ProjectTraction()` | `basis^T * tau_global` |
| `EmbedSlip()` | `basis * (0, s1, s2)^T` |
| Per-face storage | Per-quad-point D×D matrix |
