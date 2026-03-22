# BP5 Debug v50: Production Defaults — Automatic CFL dt + V-guard

**Date**: 2026-03-22
**Status**: IMPLEMENTED. CFL-aware dt_init and V-guard are now automatic defaults. No manual flags needed for production runs.
**Previous**: v49 (CFL stability analysis, Phase 1-3 testing across all configurations)
**Branch**: `feature/elasticity`

---

## 1. Background

v49 identified and resolved two stability issues:

- **Issue A** (1000m serial crash): IP penalty stiffness scales as 1/h, making the default dt_init=0.13s too large for h=1000m. The RK45 stages cascade exponentially, producing GPa-level traction and segfault.
- **Issue B** (2500m parallel crash): Same cascade mechanism, amplified by shared face DG coupling in parallel.

Both issues were fixed by `--dt-init 0.05 --v-guard 100` flags, proven across all configurations (serial, 8/48/400 ranks, 1000m/2500m meshes). v50 makes these fixes the automatic default.

---

## 2. Code Changes

### 2.1 File: `tests/verification/bp5_verification_full.cpp`

**Change 1: h_min computation** (line 579-581, after mesh loading)

```cpp
// v50: Compute minimum element size for CFL-aware dt
real_t h_min, h_max, kappa_min, kappa_max;
pmesh.GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
```

Uses MFEM's `Mesh::GetCharacteristics()` which iterates all elements and computes `h = pow(|detJ|, 1/dim)` per element, returning the global minimum. Printed to stdout:
```
h_min = XXX m, h_max = XXX m
```

**Change 2: CFL-aware dt_init** (lines 1017-1031)

Old default:
```cpp
dt_init = min(1e3, 0.01 * L_nuc / max(V_max_init, 1e-20))  // = 0.13s always
```

New default:
```cpp
real_t dt_V   = min(1e3, 0.01 * L_nuc / max(V_max_init, 1e-20));
real_t dt_CFL = 2.0 * params.eta() * h_min / (4.0 * params.mu());
real_t dt_init = min(dt_V, dt_CFL);
```

| Mesh | dt_V | dt_CFL | dt_init | Limiting factor |
|------|------|--------|---------|-----------------|
| 1000m | 0.13s | ~0.072s | ~0.072s | CFL |
| 2500m | 0.13s | ~0.180s | 0.13s | V (physics) |
| 5000m | 0.13s | ~0.361s | 0.13s | V (physics) |

The CFL formula: `dt_CFL = C * eta * h_min / (beta * mu)` where:
- C = 2.0 (safety factor below empirical z_crit of 2.5-3.0)
- eta = mu/(2*cs) = 4.62e6 Pa*s/m (radiation damping)
- beta = 4.0 (geometry/penalty factor for IP at p=2 with x3 correction)
- mu = 32.04e9 Pa (shear modulus)

The `--dt-init <value>` flag still works as a manual override for testing.

**Change 3: V-guard default ON** (lines 1053-1065)

Old behavior: V-guard OFF unless `--v-guard <factor>` flag provided.

New behavior: V-guard ON by default with factor=100. New flags:
- `--v-guard <factor>` — override the default factor
- `--no-v-guard` — disable V-guard entirely (for testing)

```cpp
if (!no_v_guard)
{
   real_t factor = (v_guard_factor > 0) ? v_guard_factor : 100.0;
   ode_solver.SetVGuard(factor);
}
```

### 2.2 No changes to other files

- `config/bp5_params.hpp`: `mu()` (line 54) and `eta()` (line 61) already exist.
- `solver/time_stepper.hpp`: V-guard mechanism unchanged (SetVGuard, check at each RK stage, halve dt on rejection).
- `domain/elasticity_operator.hpp`: No changes.

---

## 3. Expected Output

For a 1000m mesh production run (no manual flags):
```
ParMesh: 63451 global elements
h_min = XXX m, h_max = XXX m
...
dt_V = 0.13 s, dt_CFL = 0.072 s
Initial dt: 0.072 s (CFL-limited)
V-guard: ON (factor=100)
```

For a 2500m mesh:
```
h_min = XXX m, h_max = XXX m
...
dt_V = 0.13 s, dt_CFL = 0.18 s
Initial dt: 0.13 s (V-limited)
V-guard: ON (factor=100)
```

---

## 4. Test Plan

| Test | Config | Expected |
|------|--------|----------|
| `bp5_v50_1000m_production.sbatch` | 1000m, 400 ranks, NO manual flags | dt≈0.072s, V-guard ON, stable (matches v49m) |

---

## 5. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v49 | CFL analysis: identified cascade root cause, tested V-guard + dt fixes | Done |
| v49f-m | Confirmed fixes across all configurations (serial, parallel, 1000m, 2500m) | All STABLE |
| **v50** | **CFL dt + V-guard made automatic production defaults** | **Implemented** |
