# BP5 Debug v8: Slip Rate Sign Convention — Positive Feedback Instability

## Problem

After v7 fixes (robust bracket handling + smaller initial dt), the BP5 simulation runs without crashing but V_max grows unboundedly:

```
Step       Time [yr]        dt [s]     V_max [m/s]     EQs
----------------------------------------------------------------
     1            0.00     1.477e-02       3.749e-02       1
     2            0.00     1.257e-02       4.744e-02       1
     ...
    23            0.00     2.095e-03       2.459e+01       1
```

V_max grows from 0.03 → 25 m/s in 23 steps (~0.1 seconds). For quasi-dynamic BP5:
- V_max should be bounded by τ/η ≈ 15 MPa / 4.62 MPa·s/m ≈ 3.25 m/s
- V = 25 m/s requires τ_total ≈ 115 MPa — far exceeding the initial τ₀ ≈ 21 MPa

**Diagnosis**: The elastic traction **reinforces** τ_pre instead of opposing it (positive feedback).

## Root Cause

### Sign Convention Mismatch Between Slip Rate Direction and DG Slip BC

The coupling loop:

1. **Slip rate**: `V_vec = -(V_abs / |τ|) * τ` → V is **antiparallel** to τ
2. **tau_pre**: `tau = -τ₀ * V_init / |V_init|` → τ_pre is **antiparallel** to V_init
3. **Consistency**: V_init is parallel to V (both antiparallel to τ), so slip accumulates in +strike direction
4. **EmbedSlip**: Positive strike slip `S` → `delta_u = (0, S, 0)` in global frame
5. **DG slip BC**: `u_jump = sign * delta_u` where `sign = (nor(0) > 0) ? -1 : 1`
   - This maps positive strike slip to **right-lateral** motion
6. **Elastic traction**: Right-lateral slip produces **negative** shear traction (stress drop: σ_xy < 0)
   - `T_y = {σ}·n_x < 0` → `traction_strike < 0`
7. **Total stress**: `τ_total = τ_pre(-) + traction(-) = -(τ₀ + |traction|)`
   - **Both negative → reinforce → |τ_total| increases → V increases → UNSTABLE**

### Detailed 1D Traction Sign Derivation

For right-lateral slip (u_jump_y = -S where S > 0, via sign = -1):
- Element 1 (-x side): u_y(0⁻) = -S/2, u_y(-L) = 0 → du_y/dx = -S/(2L)
- Element 2 (+x side): u_y(0⁺) = +S/2, u_y(+L) = 0 → du_y/dx = -S/(2L)
- avg_grad(y,x) = -S/(2L)
- ε_xy = -S/(4L)
- σ_xy = 2μ · ε_xy = -μS/(2L) < 0
- T_y = σ_xy · basis.normal_x = -μS/(2L) < 0
- tau_strike = T_y · tangent2_y = -μS/(2L) < 0

### Why Tandem Works (Left-Lateral Convention)

Tandem uses the **same** V = -(V/|τ|)*τ convention, but without the MFEM `sign` factor in the slip BC:

1. **Tandem**: f_q = slip · fault_basis (no sign flip)
   - Positive strike slip → **left-lateral** motion (u₀_y > u₁_y where u₀ is -x side)
   - Left-lateral → **positive** elastic traction
   - `τ_total = τ_pre(-) + traction(+)` → **partial cancellation → STABLE**

2. **MFEM**: u_jump = sign · delta_u (with sign flip)
   - Positive strike slip → **right-lateral** motion
   - Right-lateral → **negative** elastic traction
   - `τ_total = τ_pre(-) + traction(-)` → **reinforcement → UNSTABLE**

The MFEM `sign` factor, which ensures consistent element ordering, flips the physical slip direction compared to Tandem.

## Fix

Make V **parallel** to τ (instead of antiparallel), and τ_pre **parallel** to V_init. This restores correct feedback:

- V parallel to τ (+strike) → positive slip → right-lateral (via `sign`)
- Right-lateral → negative elastic traction
- `τ_total = τ_pre(+) + traction(-) → decreasing → STABLE`

### Change 1: `miniapps/seas/friction/dieterich_ruina.hpp` (line 396-397)

```cpp
// BEFORE (antiparallel — causes positive feedback with MFEM sign convention):
V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];

// AFTER (parallel — correct stable feedback):
V_vec[0] = (V_abs / tau_abs) * tau_vec[0];
V_vec[1] = (V_abs / tau_abs) * tau_vec[1];
```

Update comment (line 371-372): "Slip velocity is parallel to traction (slip occurs in the direction of driving stress)."

### Change 2: `miniapps/seas/config/bp5_params.hpp` (line 305-306)

```cpp
// BEFORE (antiparallel to V_init):
tau[0] = -tau0_scalar * Vi[0] / Vi_abs;
tau[1] = -tau0_scalar * Vi[1] / Vi_abs;

// AFTER (parallel to V_init):
tau[0] = tau0_scalar * Vi[0] / Vi_abs;
tau[1] = tau0_scalar * Vi[1] / Vi_abs;
```

Update comment (line 274): "Pre-stress direction is parallel to initial velocity."

## Why Only These Two Changes Suffice

| Component | Impact | Change needed? |
|-----------|--------|----------------|
| Scalar friction solver (`SolveSlipRatePsi`) | Returns V_abs > 0 | No |
| Aging law (`Rate(V_abs, psi, Dc)`) | Takes V_abs | No |
| `InitialStatePsi(tau0, V_init, ...)` | Scalar, sign-independent | No |
| Stress equilibrium check | Uses |τ_total| vs σ_n·f + η·V | No |
| Below-fault rate `(0, Vp)` | Already positive strike | No |
| Slip assembly (`sign * delta_u`) | Sign convention unchanged | No |
| `ComputeTraction` (`sign * delta_u`) | Same sign convention | No |
| `EmbedSlip` / `ProjectTraction` | Direction-agnostic transforms | No |

### Verification of Initialization Consistency

With the fix, at t=0 (slip = 0, traction ≈ 0):
- `τ_pre = (0, +τ₀)` with τ₀ ≈ 21 MPa (positive strike)
- `τ_total = τ_pre ≈ (0, +τ₀)`, |τ_total| = τ₀
- `V_vec = (V_abs / τ₀) · (0, +τ₀) = (0, +V_abs)` → V positive, matches V_init ✓
- `psi₀ = InitialStatePsi(τ₀, V_nuc, σ_n, η, a)` → same scalar computation ✓
- Stress equilibrium: |τ_total| = σ_n·f(V_abs, ψ) + η·V_abs → unchanged ✓

## Expected Result After Fix

- V_max during earthquake: peaks at ~1–5 m/s (bounded by τ/η)
- V_max during interseismic: ~V_p = 1e-9 m/s
- Earthquake cycle: nucleation → rupture → arrest → healing → repeat

## Verification Steps

1. Build: `conda activate mfem-dev && make seas_bp5_full`
2. Unit tests: existing tests should pass (stress equilibrium, basis transforms)
3. Smoke test: small BP5 mesh — V_max stabilizes during earthquake
4. Frontera: full 1000m mesh job — complete earthquake cycle without blowup
