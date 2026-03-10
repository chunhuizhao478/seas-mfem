# BP5 Debug v2: Fix Boundary Loading Direction and Up Vector

## Problem

BP5 simulation runs to completion but results mismatch Tandem benchmark data.
Shear stress variation is opposite at some locations.

## Root Cause 1 — Dirichlet BC applied at wrong boundary walls

**File**: `elasticity_operator.hpp`, `AssembleDirichletLoading()` and `SetupBoundaries()`

### The bug

Loading was applied at attr 3,4 (y=±Ly, along-strike walls) with `u_D[1] = ±Vp*t/2`.
This creates ∂u_y/∂y → ε_{yy} → σ_{yy} (normal stress in y-direction).

The fault traction T = σ·n with n=(1,0,0) picks up only σ_{x*} components.
σ_{yy} does NOT contribute to fault traction at all. The loading was mechanically
disconnected from the fault.

### Correct behavior

Loading should be at attr 1,2 (x=±Lx, fault-normal walls) with along-strike
displacement u_y = ±Vp*t/2. This creates:
- ∂u_y/∂x → ε_{xy} → σ_{xy}
- Fault traction: T_y = σ_{xy} · n_x = σ_{xy}

This directly drives shear traction on the fault.

### Coordinate convention mapping

| Direction      | Tandem | MFEM  |
|---------------|--------|-------|
| Along-strike  | x      | y     |
| Fault-normal  | y      | x     |
| Depth         | z      | z     |

Tandem applies `u_x = ±Vp*t/2` at ±y walls (fault-normal walls in Tandem coords).
Mapping to MFEM: Tandem ±y (fault-normal) → MFEM ±x (attr 1,2).

### Fix

In `SetupBoundaries()`:
- Changed `attr == 3 || attr == 4` → `attr == 1 || attr == 2`

In `AssembleDirichletLoading()`:
- Changed `attr != 3 && attr != 4` → `attr != 1 && attr != 2`
- Changed `attr == 3` → `attr == 2` (positive x side = positive loading)
- Changed `attr == 4` → `attr == 1` (negative x side = negative loading)

## Root Cause 2 — Up vector reversed from Tandem convention

**File**: `elasticity_operator.hpp`, `SetupFaultInfo()` line 258

### The bug

`up = (0, 0, -1)` gives:
- `strike = up × n = (0,0,-1) × (1,0,0) = (0,-1,0)`
- Positive slip in strike direction = left-lateral (wrong)

### Correct behavior

Tandem uses `up = (0, 0, 1)` (default). With MFEM's `n = (1,0,0)`:
- `strike = (0,0,1) × (1,0,0) = (0,1,0)`
- Positive slip = right-lateral, matching BP5 convention

### Fix

Changed `up(2) = -1.0` → `up(2) = 1.0`

### Downstream effects

This flips the strike direction, which affects:
- `ProjectTraction()`: tau_local[1] (strike component) sign flips
- `EmbedSlip()`: slip[1] (strike component) maps to correct y-direction
- `rate_state_fault.hpp`: below-fault prescribed rate `rate(1) = Vp_bp5_` — with
  corrected strike=(0,1,0), positive Vp is correct (right-lateral plate motion)

## Not Fixed — DG penalty term in ComputeTraction

MFEM computes `t = {σ}·n` (average stress only). Tandem includes a DG penalty
term: `t = {σ}·n + penalty*([[u]] - δ)`. However, the BP1 antiplane implementation
also lacks the penalty term and works correctly. This is not a critical cause of
mismatch and was not implemented.

## Verification

1. Build: `conda activate mfem-dev && make -j4`
2. Unit tests: `./seas_test_vector_friction` (should pass 45/45)
3. Run BP5 with dx=1000m mesh
4. Compare shear stress time histories against Tandem benchmark data
