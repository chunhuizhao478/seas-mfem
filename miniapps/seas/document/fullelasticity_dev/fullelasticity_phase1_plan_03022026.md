# Phase 1: BP5 Parameters, Vector Friction, and FaultBasis

**Date**: 2026-03-02
**Status**: COMPLETE (235+ tests passing)

## Context
Extending the SEAS miniapp from 2D antiplane (BP1/BP2) to 3D full elasticity (BP5-QD). Phase 1 creates the BP5 parameter definitions, adds vector slip rate solving to the friction law, implements `FaultBasis` for coordinate transforms, extends `DomainOperator` interface for 3D, and verifies everything with unit tests.

## Files to Create

### 1. `config/bp5_params.hpp` — BP5 parameter struct

Model after `config/bp2_params.hpp` (line-for-line same style). Contains `struct BP5Params`:

**Constants (SI units, from SCEC BP5 spec Table 1):**
- Material: `rho=2670`, `cs=3464`, `nu=0.25`, `mu()=rho*cs²`, `lambda()=2*nu*mu/(1-2*nu)`
- Friction: `V0=1e-6`, `f0=0.6`, `b=0.03`, `sigma_n=25e6`
- Slip distance: `L0=0.14`, `L_nuc=0.13`
- Rate-state a: `a0=0.004`, `amax=0.04`
- Geometry (m): `hs=2e3`, `ht=2e3`, `H=12e3`, `l_vw=60e3`, `w_nuc=12e3`, `Wf=40e3`, `lf=100e3`
- Loading: `Vp=1e-9`, `V_init=1e-9`, `V_zero=1e-20`, `V_nuc=0.01`
- Time: `t_final=1800*seconds_per_year`

**Methods:**

`a_of_x2_x3(x2, x3)` — BP5 spec Eq. 14 (3 zones):
```
VW core:   a0 when (hs+ht <= x3 <= hs+ht+H) AND (|x2| <= l_vw/2)
VS zones:  amax when x3 <= hs OR x3 >= hs+2*ht+H OR |x2| >= l_vw/2+ht
Transition: a0 + r*(amax-a0), r = max(|x3-hs-ht-H/2|-H/2, |x2|-l_vw/2) / ht
```
Note: x3 is depth-positive-downward (SCEC convention), range [0, Wf].

`L_of_x2_x3(x2, x3)` — Returns `L_nuc=0.13` in nucleation zone, `L0=0.14` elsewhere.

`IsNucleationZone(x2, x3)` — True when `hs+ht <= x3 <= hs+ht+H` AND `-l_vw/2 <= x2 <= -l_vw/2+w_nuc`.

`V_init_vec(x2, x3, V[2])` — Following Tandem's bp5.lua:
- Nucleation zone: `V = {V_zero, V_nuc}` (strike≈0, dip=0.01)
- Elsewhere: `V = {V_zero, Vp}` (strike≈0, dip=plate rate)

`tau0_vec(x2, x3, eta, tau[2])` — Pre-stress for initial equilibrium. Following Tandem's tau_pre():
```cpp
real_t Vi1, Vi2;  V_init_vec(x2, x3, {Vi1, Vi2});
real_t Vi = sqrt(Vi1*Vi1 + Vi2*Vi2);
real_t a = a_of_x2_x3(x2, x3);
real_t psi_ss = f0 + b*log(V0/Vp);  // steady-state psi at plate rate
real_t e = exp(psi_ss / a);
real_t tau0 = sigma_n * a * asinh((Vi2/(2*V0)) * e) + eta * Vi2;
tau[0] = -tau0 * Vi1 / Vi;
tau[1] = -tau0 * Vi2 / Vi;
```

`eta()` — Radiation damping: `mu() / (2*cs)`.

`Print(os)` — Formatted parameter output (same style as BP2Params::Print).

### 2. `friction/dieterich_ruina.hpp` — Add `SolveSlipRateVectorPsi()`

Add one new method to `DieterichRuinaFriction` (after `SolveSlipRatePsi`, ~line 383):

```cpp
/// Solve for 2-component slip rate given 2-component traction and scalar psi.
///
/// Algorithm (following Tandem's DieterichRuinaAgeing::slip_rate):
/// 1. tau_abs = ||tau_vec||
/// 2. V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a)
/// 3. V_vec = -(V_abs / tau_abs) * tau_vec
///
/// The negative sign means slip velocity is anti-parallel to traction
/// (friction opposes applied stress, slip occurs in stress direction).
void SolveSlipRateVectorPsi(const real_t tau_vec[2], real_t psi,
                             real_t sigma_n, real_t eta, real_t a,
                             real_t V_vec[2],
                             int *iterations = nullptr) const
{
    real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);
    if (tau_abs < 1e-30) {
        V_vec[0] = 0.0;
        V_vec[1] = 0.0;
        if (iterations) *iterations = 0;
        return;
    }
    real_t V_abs = SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a, iterations);
    V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
    V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];
}
```

This reuses the existing scalar Newton solver — no new numerical code needed.

### 3. `tests/unit/test_bp5_params.cpp` — BP5 parameter tests

Uses existing test macros (TEST_ASSERT, TEST_NEAR, TEST_REL_NEAR) from `test_friction_law.cpp`.

**Tests:**
1. **Material properties**: mu = rho*cs² = 32.04 GPa, lambda = mu (for nu=0.25), eta = mu/(2*cs)
2. **a() in VW core**: a(0, 10e3) = a0 = 0.004 (center of VW zone)
3. **a() in VS zones**: a(0, 0) = amax, a(0, 50e3) = amax, a(60e3, 10e3) = amax
4. **a() in transition**: a(0, 3e3) should be between a0 and amax (hs < x3 < hs+ht)
5. **a() transition continuity**: a at VW/transition boundary = a0, at transition/VS boundary = amax
6. **L() in nucleation zone**: L(-25e3, 10e3) = 0.13
7. **L() outside nucleation**: L(0, 10e3) = 0.14
8. **IsNucleationZone**: true for (-25e3, 10e3), false for (0, 10e3)
9. **V_init_vec outside nucleation**: V ≈ (V_zero, Vp)
10. **V_init_vec in nucleation**: V ≈ (V_zero, V_nuc)
11. **tau0_vec direction**: anti-parallel to V_init (negative sign)
12. **tau0_vec self-consistency**: Feed tau0_vec + sigma_n + psi_ss into SolveSlipRateVectorPsi, verify it returns V_init

### 4. `tests/unit/test_vector_friction.cpp` — Vector friction tests

**Tests:**
1. **Pure x-direction stress**: tau = (tau_x, 0) → V should be anti-parallel: V = (-V_abs, 0)
2. **Pure y-direction stress**: tau = (0, tau_y) → V = (0, -V_abs)
3. **45-degree stress**: tau = (t, t) → V direction at 225° (anti-parallel)
4. **Magnitude consistency**: ||V_vec|| from vector solver == V_scalar from SolveSlipRatePsi(||tau||)
5. **Zero traction**: tau = (0, 0) → V = (0, 0)
6. **Known BP5 values**: Use BP5Params, compute psi_ss, tau0, and verify V_init recovery

### 5. `Makefile` — Add new targets

Add to FRICTION_HEADERS: `config/bp5_params.hpp`

Add new source/object/target entries following the existing pattern:
```makefile
TEST_BP5_PARAMS_SRC = tests/unit/test_bp5_params.cpp
TEST_VECTOR_FRICTION_SRC = tests/unit/test_vector_friction.cpp
```

Add to SEQ_MINIAPPS: `seas_test_bp5_params seas_test_vector_friction`

Add build rules (depend on FRICTION_HEADERS), test targets:
```makefile
test-bp5-params: seas_test_bp5_params
	./seas_test_bp5_params

test-vector-friction: seas_test_vector_friction
	./seas_test_vector_friction
```

Add both to the `test:` target.

## Verification

1. `conda activate mfem-dev`
2. `cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas`
3. `make seas_test_bp5_params seas_test_vector_friction` — must compile cleanly
4. `make test-bp5-params` — all BP5 param tests pass
5. `make test-vector-friction` — all vector friction tests pass
6. `make test-friction` — existing friction tests still pass (regression check)
7. `make test-psi-state` — existing psi-state tests still pass

## Key Design Notes

- **No existing behavior changes**: `SolveSlipRateVectorPsi` is purely additive to `DieterichRuinaFriction`. All existing scalar methods untouched.
- **Tandem convention for V_init**: Following the plan document, nucleation V_nuc = 0.01 (Tandem's bp5.lua), not 0.03 (SCEC spec literal text). The SCEC spec says 0.03, but Tandem community codes use 0.01. This can be changed by adjusting the `V_nuc` constant.
- **SCEC coordinate convention**: x1=fault-normal, x2=along-strike, x3=depth (positive downward). All BP5Params methods use this convention.
- **Units**: All SI (meters, seconds, Pascals) — consistent with BP2Params.

---

## Implementation Status

| Component | File | Tests | Status |
|-----------|------|-------|--------|
| BP5 parameters | `config/bp5_params.hpp` | 77 tests | COMPLETE |
| Vector friction | `friction/dieterich_ruina.hpp` | 45 tests | COMPLETE |
| Test macros | `tests/unit/test_macros.hpp` | — | COMPLETE |
| Makefile updates | `Makefile` | — | COMPLETE |

All Phase 1 tests pass. All pre-existing tests (friction, psi-state, fault detection, domain interface) continue to pass — zero regressions.

See `fullelasticity_phase1_check_03022026.md` for detailed code review.
See `fullelasticity_phase2a_plan_03022026.md` for FaultBasis and DomainOperator interface extensions (Phase 2a).
