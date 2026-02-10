# Analysis of Tandem Boundary Condition Implementation for SEAS Benchmarks

## Executive Summary

This report analyzes how Tandem implements boundary conditions for the BP1 and BP2 SEAS benchmarks, with a focus on understanding where the plate loading (Vp) is applied.

**Key Finding:** Tandem does **NOT** apply the plate loading as a Dirichlet BC at the far-field x-boundaries. Instead, it uses a combination of:
1. Dirichlet BC at **top/bottom** boundaries (z = 0 and z = Lz) with u = Vp/2 * t
2. Natural BC (zero traction) at the **right** boundary (x = Lx)
3. Fault BC at the left boundary (x = 0) for z < Wf
4. Natural BC (zero traction) at the left boundary (x = 0) for z ≥ Wf (below fault)

This approach differs significantly from a naive implementation that applies Vp at x = ±∞.

---

## 1. Investigation Methodology

### Files Examined

**Configuration Files:**
- `/Users/chunhuizhao/projects/tandem/examples/tandem/2d/bp1.lua` - BP1 scenario definition
- `/Users/chunhuizhao/projects/tandem/examples/tandem/2d/bp1_sym.toml` - BP1 symmetric configuration
- `/Users/chunhuizhao/projects/tandem/examples/tandem/2d/bp1_sym.geo` - BP1 mesh geometry

**Source Code:**
- `/Users/chunhuizhao/projects/tandem/src/io/GlobalSimplexMeshBuilder.cpp` - BC type mapping
- `/Users/chunhuizhao/projects/tandem/app/localoperator/Poisson.cpp` - BC implementation
- `/Users/chunhuizhao/projects/tandem/app/form/SeasQDOperator.cpp` - Quasi-dynamic solver

### Note on BP2

Tandem does not have a specific BP2 example in its repository. However, BP2 and BP1 share the same fundamental domain structure (2D antiplane shear), differing only in:
- Critical slip distance (Dc = 0.004 m for BP2, 0.008 m for BP1)
- Some geometric parameters

The boundary condition approach would be identical for both benchmarks.

---

## 2. Tandem Mesh Geometry (bp1_sym.geo)

```
Domain: [0, 400km] × [-400km, 0]  (right half of full domain)

          z = 0 (free surface)
          ┌─────────────────────────────┐
          │                             │
          │                             │
     x=0  │      Physical Curve(1)      │  x = 400km
  (fault) │      (Dirichlet BC)         │  (right boundary)
          │                             │
          │                             │
          └─────────────────────────────┘
          z = -400km (bottom)

Physical Curve Mapping:
  - Physical Curve(1) = BC::Dirichlet (enum value 1)
  - Physical Curve(3) = BC::Fault (enum value 3)
  - Physical Curve(5) = BC::Natural (enum value 5)
```

### Mesh Point Layout

```cpp
// From bp1_sym.geo
Point(2) = {w, 0, 0, h};      // (400, 0) - top-right corner
Point(3) = {w, -d, 0, h};     // (400, -400) - bottom-right corner
Point(5) = {0, -d, 0, h};     // (0, -400) - bottom-left corner
Point(6) = {0, 0, 0, hf};     // (0, 0) - top-left corner (fault top)
Point(7) = {0, -d4, 0, hf};   // (0, -40) - fault bottom
Point(8) = {0, -d3, 0, hf};   // (0, -18) - fault zone point
Point(9) = {0, -d2, 0, hf};   // (0, -16) - fault zone point
Point(10) = {0, -d1, 0, hf};  // (0, -15) - VW/VS transition
```

### Boundary Assignments

| Physical Curve | Lines | Location | BC Type |
|----------------|-------|----------|---------|
| 1 | 2, 4 | Top (z=0) and Bottom (z=-400km) | Dirichlet |
| 3 | 8, 9, 10, 11 | Left boundary, z=0 to z=-40km | Fault |
| 5 | 3, 7 | Right boundary + Left below z=-40km | Natural |

---

## 3. Boundary Condition Implementation

### 3.1 Dirichlet BC (Physical Curve 1)

**Location:** Top and bottom boundaries (z = 0 and z = -Lz)

**Implementation in bp1.lua:**
```lua
function BP1:boundary(x, y, t)
    if x > 1.0 then
        return self.Vp/2.0 * t      -- Right side of domain
    elseif x < -1.0 then
        return -self.Vp/2.0 * t     -- Left side (unused in half-domain)
    else
        return self.Vp * t          -- Transition region near fault
    end
end

-- Symmetric version (used in bp1_sym):
function bp1_sym:boundary(x, y, t)
    return self.Vp/2.0 * t          -- Uniform displacement everywhere
end
```

**Physical Interpretation:**
- The displacement at top/bottom boundaries is set to u = Vp/2 * t
- This represents the cumulative slip from far-field plate motion
- Since this is a half-domain (x ≥ 0), the relative motion is Vp/2

### 3.2 Fault BC (Physical Curve 3)

**Location:** Left boundary, z = 0 to z = -40km (rate-state fault region)

**Implementation in Poisson.cpp:**
```cpp
bool Poisson::bc_boundary(std::size_t fctNo, BC bc, double f_q_raw[]) const {
    auto f_q = Matrix<double>(f_q_raw, 1, tensor::f_q::Shape[0]);
    if (bc == BC::Fault) {
        fun_slip(fctNo, f_q, true);
        for (std::size_t q = 0; q < tensor::f_q::Shape[0]; ++q) {
            f_q(0, q) *= 0.5;  // Half-slip for symmetric domain
        }
    } else if (bc == BC::Dirichlet) {
        fun_dirichlet(fctNo, f_q, true);
    } else {
        return false;
    }
    return true;
}
```

**Physical Interpretation:**
- The fault boundary uses slip as the boundary value
- The factor of 0.5 accounts for the symmetric half-domain formulation
- The slip is computed from the rate-and-state friction evolution

### 3.3 Natural BC (Physical Curve 5)

**Location:**
- Right boundary (x = Lx)
- Left boundary below z = -40km (below rate-state region)

**Implementation in Poisson.cpp:**
```cpp
bool Poisson::assemble_boundary(...) const {
    if (info.bc == BC::Natural) {
        return false;  // No contribution to stiffness matrix
    }
    // ... Dirichlet/Fault handling
}

bool Poisson::rhs_boundary(...) const {
    if (!bc_boundary(fctNo, info.bc, f_q_raw)) {
        return false;  // bc_boundary returns false for Natural BC
    }
    // ... RHS assembly
}
```

**Physical Interpretation:**
- Natural BC = zero traction (Neumann BC with ∂u/∂n = 0)
- At right boundary: represents the "far-field" where stress perturbations vanish
- Below fault: represents the creeping region where no friction law is applied

---

## 4. Physical Interpretation of Tandem's Approach

### 4.1 Why Top/Bottom Dirichlet Instead of Far-Field?

The Tandem approach is physically motivated by the half-space formulation:

```
                Far-field (x → ∞)
                      ↑
                      │ u → Vp*t/2
    ┌─────────────────┼─────────────────┐
    │                 │                 │
    │   Material      │   Material      │
    │   moves with    │   moves with    │
    │   fault         │   plate         │
    │                 │                 │
    ├─────────────────┤                 │
    │     FAULT       │                 │
    │   (slip = δ)    │                 │
    ├─────────────────┤                 │
    │   Creeping      │                 │
    │   region        │                 │
    │   (free slip)   │                 │
    └─────────────────┴─────────────────┘
```

**Key insight:** In a 2D antiplane problem with a vertical fault:
- The loading comes from relative motion between the two half-spaces
- This manifests as displacement at z = 0 and z = Lz, not at x = Lx
- The right boundary (x = Lx) is the "far-field" where traction → 0

### 4.2 Why Natural BC Below the Fault?

Below the rate-state fault region (z > Wf = 40km), the fault creeps freely:
- No friction law is applied (not velocity-weakening)
- The fault can slip freely to accommodate plate motion
- This is implemented as zero traction (Natural BC)

---

## 5. Comparison with Current MFEM Implementation

### Current Implementation (antiplane_operator.hpp)

```cpp
void ApplyEssentialBCs(real_t time, const Vector &slip_bc, GridFuncType &u)
{
    real_t u_farfield = 0.5 * Vp_ * time;

    // Far-field left BC: u = -Vp*t/2 at x = -Lx
    // Far-field right BC: u = +Vp*t/2 at x = +Lx  ← INCORRECT
    // Fault slip BC: u = slip/2 on fault
    // Bottom BC: u = 0  ← INCORRECT
}
```

### Issues with Current Approach

1. **Far-field BC at x = Lx:** Should be Natural BC (zero traction), not Dirichlet
2. **Bottom BC at z = Lz:** Should be Dirichlet with u = Vp/2 * t, not u = 0
3. **Below fault region:** Should be Natural BC (zero traction) for z > Wf

---

## 6. Recommended Corrections for MFEM Implementation

### 6.1 Updated Boundary Attribute Convention

```cpp
struct BP2BoundaryAttributes
{
    static constexpr int FAULT = 1;          // x = 0, z ∈ [0, Wf]
    static constexpr int CREEP = 2;          // x = 0, z ∈ [Wf, Lz] (below fault)
    static constexpr int FARFIELD_RIGHT = 3; // x = Lx (zero traction)
    static constexpr int FREE_SURFACE = 4;   // z = 0 (zero traction)
    static constexpr int BOTTOM = 5;         // z = Lz (plate loading)
};
```

### 6.2 Corrected Boundary Conditions

| Boundary | Old BC | New BC (Tandem-style) |
|----------|--------|----------------------|
| Fault (x=0, z<Wf) | Dirichlet: u = slip/2 | Dirichlet: u = slip/2 (same) |
| Below fault (x=0, z≥Wf) | N/A | Natural: ∂u/∂n = 0 |
| Right (x=Lx) | Dirichlet: u = Vp*t/2 | Natural: ∂u/∂n = 0 |
| Free surface (z=0) | Natural | Natural: ∂u/∂z = 0 (same) |
| Bottom (z=Lz) | Dirichlet: u = 0 | Dirichlet: u = Vp*t/2 |

### 6.3 Updated ApplyEssentialBCs

```cpp
void ApplyEssentialBCs(real_t time, const Vector &slip_bc, GridFuncType &u)
{
    // Only essential (Dirichlet) BCs - Natural BCs handled automatically

    // 1. Fault slip BC: u = slip/2 on fault (x = 0, z < Wf)
    for (int i = 0; i < fault_dofs_.Size(); i++) {
        u(fault_dofs_[i]) = 0.5 * slip_bc(i);
    }

    // 2. Bottom BC: u = Vp*t/2 at z = Lz (plate loading)
    real_t u_loading = 0.5 * Vp_ * time;
    for (int i = 0; i < bottom_dofs_.Size(); i++) {
        u(bottom_dofs_[i]) = u_loading;
    }

    // 3. Free surface (z=0): Natural BC - no action needed
    // 4. Right boundary (x=Lx): Natural BC - no action needed
    // 5. Below fault (x=0, z≥Wf): Natural BC - no action needed
}
```

---

## 7. Full Domain vs Half Domain Formulations

### 7.1 Half Domain (Symmetric) Formulation

Tandem's `bp1_sym` example uses a **half domain** with symmetry:
- Domain: x ∈ [0, Lx], z ∈ [0, Lz]
- Fault at x = 0 (left boundary)
- Exploits symmetry: only solve one side of the fault

**Half Domain Boundary Conditions:**
| Boundary | Condition |
|----------|-----------|
| Fault (x=0, z<Wf) | Dirichlet: u = slip/2 |
| Creep (x=0, z≥Wf) | Natural (zero traction) |
| Right (x=Lx) | Natural (zero traction) |
| Free surface (z=0) | Natural (zero traction) |
| Bottom (z=Lz) | Dirichlet: u = Vp·t/2 |

### 7.2 Full Domain Formulation (Recommended)

Tandem's `bp1` (non-symmetric) example uses the **full domain**:
- Domain: x ∈ [-Lx, +Lx], z ∈ [0, Lz]
- Fault at x = 0 (internal interface)
- Opposite displacements on each side of the fault

**Full Domain Boundary Conditions:**
| Boundary | Condition |
|----------|-----------|
| Fault left (x=0⁻, z<Wf) | Dirichlet: u = -slip/2 |
| Fault right (x=0⁺, z<Wf) | Dirichlet: u = +slip/2 |
| Creep left (x=0⁻, z≥Wf) | Natural (zero traction) |
| Creep right (x=0⁺, z≥Wf) | Natural (zero traction) |
| Far-field left (x=-Lx) | Natural (zero traction) |
| Far-field right (x=+Lx) | Natural (zero traction) |
| Free surface (z=0) | Natural (zero traction) |
| Bottom (z=Lz) | Dirichlet: u = sign(x)·Vp·t/2 |

**Tandem's boundary function for full domain (bp1.lua):**
```lua
function BP1:boundary(x, y, t)
    if x > 1.0 then
        return self.Vp/2.0 * t      -- Right side: +Vp/2 * t
    elseif x < -1.0 then
        return -self.Vp/2.0 * t     -- Left side: -Vp/2 * t
    else
        return self.Vp * t          -- Near fault: transition
    end
end
```

### 7.3 Why Full Domain?

The full domain formulation is preferred because:
1. **Physical completeness**: Models both sides of the fault explicitly
2. **No symmetry assumption**: Required for non-symmetric problems (e.g., dipping faults)
3. **Correct traction computation**: Traction τ = μ·∂u/∂x includes contributions from both sides
4. **Consistent with Tandem**: Matches Tandem's general approach

### 7.4 DG Implementation for Full Domain

In DG (Discontinuous Galerkin), the full domain is natural:
- The fault at x=0 can be an **interior interface** with a prescribed jump
- Or the mesh can be **split at x=0** with the fault as boundary faces

**Option A: Interior Fault (Tandem approach)**
- Single mesh spanning [-Lx, +Lx]
- Fault faces at x≈0 are interior faces
- Slip imposed as jump condition: [[u]] = u⁺ - u⁻ = slip

**Option B: Split Mesh (simpler for MFEM)**
- Mesh has a "seam" at x=0 with duplicate nodes
- Fault faces are boundary faces of each half
- Left side: u = -slip/2, Right side: u = +slip/2

---

## 8. MFEM Full Domain Implementation

### 8.1 Mesh Generation

Create mesh for x ∈ [-Lx, +Lx], z ∈ [0, Lz] with fault at x=0:

```cpp
// Create mesh centered at x=0
Mesh::MakeCartesian2D(2*nx, nz, Element::QUADRILATERAL, true, 2*Lx, Lz);
// Shift x-coordinates: x → x - Lx
for (int i = 0; i < mesh.GetNV(); i++) {
    real_t *v = mesh.GetVertex(i);
    v[0] -= Lx;  // Now x ∈ [-Lx, +Lx]
}
```

### 8.2 Boundary Attributes

```cpp
struct BP2BoundaryAttributes {
    static constexpr int FAULT_LEFT = 1;     // x=0⁻, z<Wf: Dirichlet u=-slip/2
    static constexpr int FAULT_RIGHT = 2;    // x=0⁺, z<Wf: Dirichlet u=+slip/2
    static constexpr int CREEP_LEFT = 3;     // x=0⁻, z≥Wf: Natural
    static constexpr int CREEP_RIGHT = 4;    // x=0⁺, z≥Wf: Natural
    static constexpr int FARFIELD_LEFT = 5;  // x=-Lx: Natural
    static constexpr int FARFIELD_RIGHT = 6; // x=+Lx: Natural
    static constexpr int FREE_SURFACE = 7;   // z=0: Natural
    static constexpr int BOTTOM = 8;         // z=Lz: Dirichlet u=sign(x)·Vp·t/2
};
```

### 8.3 Position-Dependent Bottom BC

At the bottom boundary (z=Lz), displacement varies with x:
```cpp
real_t GetBottomBC(real_t x, real_t time, real_t Vp) {
    if (x > 0) return +0.5 * Vp * time;
    else if (x < 0) return -0.5 * Vp * time;
    else return 0.0;  // At fault, transition
}
```

---

## 9. Summary

| Aspect | Half Domain (bp1_sym) | Full Domain (bp1) | MFEM Implementation |
|--------|----------------------|-------------------|---------------------|
| Domain | x ∈ [0, Lx] | x ∈ [-Lx, +Lx] | **Full domain** |
| Fault location | x=0 (boundary) | x=0 (internal) | x=0 (split boundary) |
| Fault BC | u = slip/2 | u = ±slip/2 | u = ±slip/2 |
| Bottom BC | u = Vp·t/2 | u = sign(x)·Vp·t/2 | u = sign(x)·Vp·t/2 |
| Far-field | Natural (x=Lx) | Natural (x=±Lx) | Natural (x=±Lx) |

**MFEM Implementation**: Uses full domain with split mesh at x=0, following Tandem's `bp1` approach.

---

## Appendix: Tandem Source Code References

### A.1 BC Enum Mapping (GlobalSimplexMeshBuilder.cpp:80-96)
```cpp
switch (tag) {
    case static_cast<long>(BC::None):     bc = BC::None; break;
    case static_cast<long>(BC::Dirichlet): bc = BC::Dirichlet; break;
    case static_cast<long>(BC::Fault):     bc = BC::Fault; break;
    case static_cast<long>(BC::Natural):   bc = BC::Natural; break;
    default: ++unknownBC; break;
}
```

### A.2 Boundary Function Application (SeasQDOperator.cpp:57-61)
```cpp
void SeasQDOperator::solve(double time, BlockView const& state_view) {
    dgop_->set_slip(adapter_->slip_bc(state_view));
    if (fun_boundary_) {
        dgop_->set_dirichlet((*fun_boundary_)(time));
    }
    linear_solver_.update_rhs(*dgop_);
    linear_solver_.solve();
}
```

### A.3 BP1 Parameters (bp1.lua)
```lua
BP1.Vp = 1e-9       -- Plate velocity: 1 nm/s
BP1.H = 15.0        -- VW zone depth: 15 km
BP1.h = 3.0         -- Transition width: 3 km
-- Wf = 40 km (implicit from mesh)
```
