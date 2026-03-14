# BP5 Debug v21: Boundary Loading Fix — Tandem Has Natural BC on Top/Bottom

**Date**: 2026-03-13
**Fix**: H25 — Correct boundary condition assignment to match Tandem

## Problem: MFEM Applies Dirichlet on ALL Faces, Tandem Does NOT

### Tandem's Actual BP5 BCs (from `examples/tandem/3d/bp5.geo`)

```
diri() = Surface{:};       // start with ALL surfaces
diri() -= top();           // REMOVE z=0  (Earth's surface)
diri() -= bottom();        // REMOVE z=-100 (deep boundary)
diri() -= fault();         // REMOVE fault

Physical Surface(1) = {bottom(),top()};  → BC::Natural = 1  → ZERO TRACTION
Physical Surface(3) = {fault()};         → BC::Fault = 3    → Rate-state friction
Physical Surface(5) = {diri()};          → BC::Dirichlet = 5 → plate loading
```

BC enum from `src/form/BC.h`:
```cpp
enum class BC : int { None = 0, Natural = 1, Fault = 3, Dirichlet = 5 };
```

The physical surface numbers ARE the BC enum values. So:
- **Surface 1 (top + bottom) → Natural (zero traction)**
- **Surface 3 (fault) → Fault (friction)**
- **Surface 5 (far-field vertical faces) → Dirichlet (plate loading)**

### Coordinate Mapping (Tandem → MFEM)

| Tandem | MFEM | Range |
|--------|------|-------|
| X (along-strike) | y | [-200, 200] km |
| Y (fault-normal) | x | [-100, 100] km |
| Z (depth, negative down) | -z (depth, positive down) | [-100, 0] → [0, 100] |

Both have z=0 as Earth's surface. The z-sign flip doesn't change surface identification.

### Tandem BC Mapping to MFEM Attributes

| MFEM attr | MFEM face | Tandem equivalent | Tandem BC |
|-----------|-----------|-------------------|-----------|
| 1 | x = -Lx (-100) | y = Y0 (-100) | Part of `diri()` → **Dirichlet** |
| 2 | x = +Lx (+100) | y = Y1 (+100) | Part of `diri()` → **Dirichlet** |
| 3 | y = +Ly (+200) | x = X1 (+200) | Part of `diri()` → **Dirichlet** |
| 4 | y = -Ly (-200) | x = X0 (-200) | Part of `diri()` → **Dirichlet** |
| 5 | z = 0 (surface) | z = 0 = `top()` | **Natural (zero traction)** |
| 6 | z = Lz (deep) | z = -100 = `bottom()` | **Natural (zero traction)** |

### The `boundary_linear` Misunderstanding

MFEM code comments said: "No free surface — consistent with Tandem's `boundary_linear=true`"

Tandem's docs (`docs/first-model/parameters.rst`) say:
> **boundary_linear**: Assert that boundary is a linear function of time (i.e. boundary(x, t) = f(x) t). Default = false.

This is just an **optimization flag** for time stepping — it has nothing to do with which faces are Dirichlet. The previous H10 fix incorrectly interpreted this as "apply Dirichlet everywhere."

### MFEM's Current Implementation (WRONG)

`SetupBoundaryMarkers()` marks ALL attrs 1-6 as Dirichlet, including z=0 (free surface) and z=Lz (deep boundary). This over-constrains the domain.

### Impact

Applying Dirichlet loading on top (free surface) and bottom:
- Over-constrains the elastic domain
- Prevents free relaxation at the surface
- Artificial stress accumulation near surface and at depth
- Could cause incorrect recurrence interval (observed 1.65 years vs expected 300-400 years)

## Fix H25: Configurable BC Mode

Add a `BCMode` enum to `ElasticityDomainOperator`:

| Mode | Dirichlet attrs | Natural attrs | Description |
|------|-----------------|---------------|-------------|
| `FarField` (default) | 1-4 (x=±Lx, y=±Ly) | 5-6 (z=0, z=Lz) | Dirichlet on far-field vertical faces |
| `XOnly` | 1-2 (x=±Lx) | 3-6 | Antiplane-style, load only fault-normal faces |
| `AllDirichlet` | 1-6 | none | Previous (wrong) implementation, kept for comparison |

### Code Changes

**File**: `domain/elasticity_operator.hpp`
- Add `BCMode` enum
- Add `bc_mode_` member and constructor parameter
- `SetupBoundaryMarkers()`: use `bc_mode_` to select which attrs are Dirichlet
- Update comments to reflect correct Tandem analysis

**File**: `tests/verification/bp5_verification_full.cpp`
- Add `--bc-mode` CLI option (tandem, x-only, all-dirichlet)
- Default to `Tandem` mode

**File**: `bp5/mesh/bp5.geo`
- Update comments to reflect correct BC assignments

## Tandem Boundary Function Reference

From `examples/tandem/3d/bp5.lua`:
```lua
function BP5:boundary(x, y, z, t)
    local Vh = self.Vp * t
    if y > 1 then       -- fault-normal positive side (MFEM x > 0)
        Vh = Vh / 2.0
    elseif y < -1 then  -- fault-normal negative side (MFEM x < 0)
        Vh = -Vh / 2.0
    end
    return Vh, 0, 0     -- along-strike displacement only
end
```

This function is only evaluated on Physical Surface 5 (Dirichlet faces = far-field vertical boundaries). It is NOT evaluated on Physical Surface 1 (top/bottom = Natural).
