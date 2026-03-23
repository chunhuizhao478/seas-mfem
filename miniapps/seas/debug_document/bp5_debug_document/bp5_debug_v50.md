# BP5 Debug v50: Production Defaults + First Earthquake + p=4 Nucleation Failure

**Date**: 2026-03-22
**Status**: 1000m p=2 production run completed first earthquake successfully (V_peak=0.79 m/s, 91% of Tandem). 2500m p=4 FAILS — nucleation V decays instead of growing.
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

## 4. v50+ Code Change: Order-Aware CFL Beta

The original CFL formula used `beta=4.0` hardcoded for p=2. For higher orders, the IP penalty coefficient `c_N_1 = p*(p+dim-1)/dim` scales quadratically, requiring a larger beta:

```cpp
real_t c_N_1 = order * (order + dim - 1.0) / dim;
real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;  // c_N_1 at p=2 = 8/3
real_t beta = 4.0 * c_N_1 / c_N_1_ref;
real_t dt_CFL = 2.0 * params.eta() * h_min / (beta * params.mu());
```

| Order | c_N_1 | beta | dt_CFL (h=2000m) | dt_CFL (h=3200m) |
|-------|-------|------|------------------|------------------|
| p=1 | 1.0 | 1.5 | ~0.19s | ~0.31s |
| p=2 | 2.67 | 4.0 | ~0.072s | ~0.12s |
| p=4 | 8.0 | 12.0 | ~0.024s | ~0.039s |
| p=6 | 16.0 | 24.0 | ~0.012s | ~0.019s |

---

## 5. First Earthquake Results: 1000m p=2 (v50 Production)

### 5.1 Run Configuration

- Mesh: bp5_tandem.msh (1000m), 400 ranks, 8 nodes
- Order: p=2, IP method, MUMPS BLR solver
- h_min = 812.6 m, h_max = 53799 m
- dt_CFL = 0.0586s (CFL-limited), V-guard ON (factor=100)
- V-guard triggered 2 rejections at step 0 (initial dt too large for first stages)

### 5.2 Earthquake Timeline

The first earthquake nucleated at (x2=-24km, x3=10km) and propagated toward the center:

| Event | Time (s) | V_max location |
|-------|----------|----------------|
| Nucleation start | 0 | x2=-24, x3=10 (V=0.01) |
| V=0.032 at nucleation | 28.4 | x2=-24 |
| V=0.1 at nucleation | 40.4 | x2=-24 |
| **Peak V=0.33 at nucleation** | **47.0** | x2=-24 (post-peak) |
| **Peak V=0.41 at x2=-16** | **50.8** | x2=-16, x3=10 |
| Rupture reaches center | 54.6 | x2=0, x3=10 (V>1e-5) |
| **Peak V=0.79 at center** | **61.0** | x2=0, x3=10 |
| Post-peak healing | 62.3 | V declining to 0.60 |

### 5.3 Comparison with Tandem p6 Reference

**Event-aligned timing (at nucleation station x2=-24, x3=10):**

| V threshold | Tandem (s) | Our v50 (s) | Delay |
|-------------|-----------|-------------|-------|
| 0.032 m/s | 24.6 | 28.4 | +3.7s |
| 0.100 m/s | 35.2 | 40.4 | +5.2s |
| 0.316 m/s | 40.5 | 46.4 | +5.9s |

**Event-aligned timing (at propagation station x2=-16, x3=10):**

| V threshold | Tandem (s) | Our v50 (s) | Delay |
|-------------|-----------|-------------|-------|
| 0.032 m/s | 34.8 | 39.9 | +5.1s |
| 0.100 m/s | 38.7 | 44.4 | +5.7s |
| 0.316 m/s | 42.0 | 48.3 | +6.2s |

**Peak V comparison across stations:**

| Station | Our Peak V | Tandem Peak V | Ratio | Timing Delay |
|---------|-----------|---------------|-------|-------------|
| x2=-24, x3=10 (nucleation) | 0.33 m/s | 0.42 m/s | 80% | +5.1s |
| x2=-16, x3=10 (propagation) | 0.41 m/s | 0.59 m/s | 69% | +5.8s |
| **x2=0, x3=10 (center)** | **0.79 m/s** | **0.87 m/s** | **91%** | **+7.8s** |

**Stress drop comparison:**

| Station | Our delta-tau | Tandem delta-tau | Match |
|---------|-------------|-----------------|-------|
| x2=-24 (nucleation) | 13.0 MPa | 12.8 MPa | 101% |
| x2=-16 (propagation) | 11.6 MPa | 11.3 MPa | 103% |
| x2=0 (center) | 10.4 MPa | 11.3 MPa | 92% (still healing) |

### 5.4 Assessment

**Correct behavior:**
- Nucleation at the right location (x2=-24, x3=10) with correct initial V=0.01
- Rupture propagation direction: x2=-24 → x2=-16 → x2=0 (correct)
- Peak V at center: 0.79 m/s (91% of Tandem's 0.87 m/s)
- Stress drops match within 1-3% at nucleation and propagation stations
- No instability, no blowup — V-guard + CFL-dt stabilization works

**Discrepancies:**
- 4-8 second timing delay that grows with distance from nucleation
- Near-field peak V deficit (69-80% of Tandem at x2=-16 and x2=-24)
- Center match is best (91%) suggesting the propagating rupture converges

**Likely causes of timing delay:**
- Different mesh resolution (1000m p=2 vs 4000m p=6)
- V-guard rejections at step 0 slowing early nucleation
- Different effective DOF spacing (500m vs 667m)

---

## 6. p=4 Nucleation Failure: 2500m Mesh (v50a)

### 6.1 Run Configuration

- Mesh: bp5_tandem_2500m.msh, 200 ranks, 4 nodes
- Order: p=4, IP method
- nbf_per_face = 15 DOFs per fault face
- CFL: c_N_1=8.0, beta=12.0

### 6.2 Observation: V Decays Instead of Growing

At the nucleation station (x2=-24, x3=10):

| Time (s) | Our V (m/s) | Tandem p4 V (m/s) | Problem |
|----------|-------------|-------------------|---------|
| 0.000 | 0.01000 | 0.01000 | Match |
| 0.001 | 0.01005 ↑ | 0.01000 | Brief increase |
| 0.006 | **0.00974** ↓ | 0.01001 | **STARTS DROPPING** |
| 0.017 | **0.00856** ↓ | 0.01002 | Continuing to decay |
| 0.060 | **0.00603** ↓ | 0.01005 | 40% below initial |
| 2.05 | **0.00251** ↓ | ~0.010 | **75% below initial** |

**The nucleation FAILS at p=4.** V decays from 0.01 to 0.0025 while Tandem p4 shows steady growth at 0.01+.

### 6.3 What is Correct

- Initial conditions match exactly: tau=21.148 MPa, V=0.01 m/s
- The first time step (t=0→0.001s) shows V briefly increasing to 0.01005 (correct direction)
- All non-nucleation stations show V ≈ plate rate (1e-9) as expected

### 6.4 What is Wrong

After the first step, V reverses direction and decays monotonically. The traction field is no longer sustaining V=0.01 — it is actively decelerating the nucleation zone.

### 6.5 Candidate Root Causes

**Candidate 1: IP penalty over-stiffening at p=4**
- c_N_1(p=4) = 8.0 vs c_N_1(p=2) = 2.67 — penalty is 3× stronger
- Stronger penalty means the elastic solve enforces the slip constraint more rigidly
- If the prescribed slip is slightly inconsistent with the stress field, the penalty correction OPPOSES the motion
- At p=2, the penalty is weak enough that the physics (friction, radiation damping) dominates
- At p=4, the penalty may dominate the physics, creating artificial resistance

**Candidate 2: Multi-DOF slip distribution error at p=4**
- nbf=15 per face (vs nbf=6 at p=2)
- Each face DOF carries an independent slip value via `EmbedSlip`
- If the L2 projection or `InterpolateToQuadPoints` distributes slip incorrectly with 15 nodes, the per-DOF penalty creates internal inconsistency
- The penalty fights the mismatched DOFs, creating resistance

**Candidate 3: FaceQuadrature accuracy at p=4**
- The `FaceQuadrature` class stores quadrature points and weights
- At p=4, the face quadrature needs order 2*4+1=9, requiring more quadrature points
- If the quadrature is insufficient, the penalty/symmetry integrals are inaccurate

**Candidate 4: Psi initialization mismatch at p=4**
- With 15 DOFs per face, each DOF gets its own psi (state variable)
- If the psi initialization depends on the traction at each DOF, and the traction computation at p=4 has errors, the initial psi could be wrong
- Wrong psi means the friction law fights the initial V, pulling it toward a different value

**Candidate 5: Tandem uses a different DG formulation**
- Tandem may use NIPG (epsilon=+1) while we use SIPG (epsilon=-1)
- Tandem may use a different penalty scaling
- The nucleation sensitivity to these differences could be amplified at higher order

### 6.6 Diagnostic Plan

1. **Compare traction at t=0+ between p=2 and p=4**: On the same mesh (2500m), run both p=2 and p=4 with `--diag-first-traction` and compare the traction at the nucleation station. If p=4 traction is significantly different from p=2 at the same physical point, the IP formulation has an order-dependent error.

2. **Test BR2 at p=4**: BR2 doesn't have the IP penalty issue. If BR2 p=4 nucleates correctly, the IP penalty formulation is the problem.

3. **Test on 4000m mesh (same as Tandem)**: Our v50b test uses the 4000m mesh with p=6. If p=6 4000m also fails, it confirms a systematic high-order issue.

4. **Reduce penalty factor**: Test with a penalty safety factor < 1 (e.g., 0.5) to see if reducing the penalty restores nucleation. If it does, the penalty magnitude is the root cause.

---

## 7. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v49 | CFL analysis: identified cascade root cause, tested V-guard + dt fixes | Done |
| v49f-m | Confirmed fixes across all configurations (serial, parallel, 1000m, 2500m) | All STABLE |
| **v50** | **CFL dt + V-guard made automatic production defaults** | **Implemented** |
| **v50 prod** | **1000m p=2 first earthquake: V_peak=0.79 m/s (91% of Tandem)** | **SUCCESS** |
| v50+ | Order-aware CFL beta scaling: beta(p) = 4.0 * c_N_1(p) / c_N_1(2) | Implemented |
| **v50a** | **2500m p=4: nucleation FAILS — V decays 0.01→0.0025 in 2s** | **BUG** |
| v50b | 4000m p=6: pending results | Pending |
