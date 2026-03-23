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

Three parallel tests on the same 2500m mesh to isolate the root cause:

1. **v50c: BR2 at p=4** — if BR2 nucleates, IP penalty is the problem
2. **v50d: IP at p=2** — if p=2 nucleates on same mesh, confirms p=4-specific
3. **v50e: IP at p=4, penalty_factor=0.5** — if half-penalty nucleates, magnitude is the issue

### 6.7 Diagnostic Results

| Test | Method | Order | Penalty | V trajectory | Result |
|------|--------|-------|---------|-------------|--------|
| **v50c** | BR2 | p=4 | N/A | 0.010 → 0.0109 (peak step 6) → **0.0000224** (step 420) | **FAILS** |
| **v50d** | IP | p=2 | 1.0 | 0.010 → **0.320** (step 3049, growing fast) | **WORKS** |
| **v50e** | IP | p=4 | **0.5** | 0.010 → **0.0173** (step 420, growing) | **WORKS** |

**Detailed V evolution:**

v50c (BR2 p=4) — initial growth then collapse:
```
Step 1:  V=0.01009 ↑  dt=0.097s
Step 6:  V=0.01093 ↑  dt=0.290s  ← peak
Step 7:  V=0.01090 ↓  dt=0.350s  ← reversal begins
Step 10: V=0.01011 ↓  dt=0.345s
Step 15: V=0.00811 ↓  dt=0.451s
Step 420: V=0.0000224 ↓  dt=84s  ← collapsed, interseismic dt
```

v50d (p=2 IP) — healthy nucleation:
```
Step 1:  V=0.01012 ↑  dt=0.016s
Step 15: V=0.01082 ↑  dt=0.041s
Step 3049: V=0.320 ↑  dt=0.015s  ← earthquake in progress
```

v50e (p=4 IP half penalty) — slow but positive nucleation:
```
Step 1:  V=0.01009 ↑  dt=0.011s
Step 15: V=0.01050 ↑  dt=0.029s
Step 420: V=0.01732 ↑  dt=0.018s  ← still growing
```

### 6.8 Analysis of Diagnostic Results

**Key finding 1: BR2 also fails at p=4** — rules out IP penalty as SOLE cause. BR2 uses lifting operators (no explicit penalty), yet V also reverses and collapses.

**Key finding 2: p=2 works on same mesh** — confirms issue is p=4-specific, not mesh-specific.

**Key finding 3: Half-penalty IP p=4 works** — reducing penalty restores nucleation.

**Key finding 4: BR2 V reversal correlates with dt growth** — V grows for 6 steps (dt 0.10→0.29s), reverses at dt≈0.35s.

### 6.9 Deep Root Cause Analysis: Catastrophic Cancellation in Traction Recovery

After tracing the full computation pipeline (ODE → elastic solve → traction → friction), the root cause is identified as **catastrophic cancellation in the operator-split traction recovery**.

#### 6.9.1 The Pipeline

The SEAS ODE right-hand side (`SEASQuasiDynamicOperator::Mult`) does:
1. `GetSlip(state, slip)` — extract per-DOF slip from state vector
2. `Solve(t, slip, u)` — assemble K*u = f(slip) and solve
3. `ComputeTraction(u, slip, traction)` — recover per-DOF fault traction from displacement
4. `ComputeRHS(traction, state, rate)` — friction law gives dV/dt, dpsi/dt

The fault operator then combines per-DOF analytical pre-stress `tau_pre` with the elastic traction:
```
tau_total = tau_pre + traction_elastic
```
where `tau_pre` is ~20 MPa (purely from coordinates, order-independent) and `traction_elastic` is a small perturbation computed from the DG displacement field.

#### 6.9.2 The Cancellation Mechanism

The IP traction formula at each quadrature point is:
```
T_q = {sigma(u)·n} - penalty * sign * ((u1-u2) - sign*delta_u)
```

This is the difference of two terms:
- **Stress traction**: `{sigma(u)·n}` — average stress from both elements
- **Penalty correction**: `penalty * (jump - prescribed_slip)` — enforces the DG constraint

The penalty at p=4 is 3× larger than at p=2 (because c_N_1(4)/c_N_1(2) = 8.0/2.67 = 3.0). With stronger penalty:
- The solver enforces `[u] ≈ delta_u` MORE TIGHTLY (smaller residual)
- BOTH terms `{sigma·n}` and `penalty*(jump-delta_u)` are individually large
- The net traction T_q is their DIFFERENCE — a small residual of two large, nearly cancelling numbers

At p=2, the penalty is moderate, so the cancellation ratio is manageable. At p=4, the penalty is 3× stronger, and the cancellation becomes severe. Any perturbation in the displacement field (solver tolerance, element coupling, polynomial constraints) causes a disproportionately large **relative** error in the net traction.

For nucleation, the net traction perturbation at each time step must have the correct sign (reducing stress to allow friction weakening to drive V increase). At p=4, the cancellation error can flip the sign of individual DOF tractions, creating artificial resistive forces that decelerate V.

#### 6.9.3 Why This Explains All Observations

1. **IP p=4 full penalty FAILS**: penalty 3× → cancellation 3× → net traction sign errors → V decays
2. **IP p=4 half penalty WORKS**: penalty 1.5× → cancellation reduced → net traction remains accurate → V grows
3. **BR2 p=4 FAILS**: BR2 uses lifting operators for stabilization. At p=4, the lifted function is a 35-DOF polynomial per element, providing higher effective stiffness. While the mechanism differs from IP penalty, the effective coupling is still stronger at p=4, creating analogous cancellation in the BR2 traction formula. The BR2 failure is delayed (V grows for 6 steps before reversing) because BR2's effective stiffness is lower than IP at full penalty.
4. **p=2 on same mesh WORKS**: c_N_1 = 2.67 → moderate penalty → manageable cancellation → correct traction signs

#### 6.9.4 Why Tandem Does Not Have This Problem

Tandem uses the SAME DG formulation (SIPG) with the SAME penalty formula. The critical difference is the **solver architecture**:

**Tandem**: Monolithic formulation — the displacement, slip rate, and state variable are coupled in a single system. The traction is computed WITHIN the solver iteration, maintaining algebraic consistency between the K-matrix penalty and the traction penalty. The cancellation is avoided because both terms come from the same linear algebra.

**SEAS-MFEM**: Operator-split formulation:
1. Slip → elastic solve → displacement (K*u = f)
2. Displacement → post-process → traction (separate evaluation)
3. Traction → friction law → dV/dt

The traction is POST-PROCESSED from the displacement field in a separate evaluation. The penalty in K and the penalty in traction recovery are computed independently. Even though they use the same formula, floating-point evaluation on different inputs (matrix assembly vs quad-point evaluation) introduces small discrepancies that are amplified by the cancellation.

#### 6.9.5 Confirmed: What is NOT the Root Cause

The deep analysis confirmed these are NOT the cause:
- **FaceQuadrature L2 projection**: Uses M_ref_inv with 2p+1 quadrature. For flat faces on linear tets, both the mass matrix and projection integral are exact. The traction (degree ≤ p polynomial) is projected exactly.
- **tau0/V_init computation**: Purely analytical functions of coordinates, order-independent.
- **Psi initialization**: Computed from `InitialStatePsi(|tau_pre + 0|, V_init, sigma_n, eta, a)` where elastic traction is exactly zero at t=0 (zero slip → zero displacement). Order-independent.
- **Sign convention**: Consistent across all paths for flat faces.
- **nbf_per_face computation**: Correctly set to `(p+1)*(p+2)/2`.

### 6.10 Clarification: What is Monolithic and What is Not

The previous analysis (Section 6.9.4) incorrectly implied our code is fully operator-split.
In fact:

| Step | Code | Monolithic? |
|------|------|-------------|
| Assemble K (domain + face integrals + penalty) | MFEM DG integrators | ✓ |
| Assemble f (slip RHS with penalty) | Our `AssembleSlipContributionIP` | ✓ |
| Solve K*u = f | MUMPS direct solver | ✓ |
| **Recover traction T from u** | **Our `ComputeTraction`** | **✗ POST-PROCESSING** |

The solve IS monolithic — MUMPS gives u exactly (up to solver tolerance ~1e-14).
The problem is Step 4: `ComputeTraction` re-evaluates the DG traction formula from
scratch using u, rather than extracting it from the linear system.

The penalty in K (computed by `DGElasticityIPPenaltyIntegrator::AssembleFaceMatrix`)
and the penalty in `ComputeTraction` (our own quad-point evaluation) are computed by
**different code paths** on **different inputs**. Even though they use the same
mathematical formula, numerical evaluation differences are amplified by the cancellation.

**Why this was not caught in previous Tandem comparisons**: The comparison focused on
formulation matching (SIPG ✓, penalty formula ✓, sign convention ✓, time stepper ✓).
The traction recovery method — HOW traction is extracted after the solve — was never
compared. We assumed same-formula consistency; the catastrophic cancellation mechanism
was not considered.

### 6.11 Fix Strategies (Revised)

**Strategy 1: Stress-only traction (`--traction-stress-only`) — simplest diagnostic**

The IP traction is:
```
T_q = {σ(u)·n} - penalty * sign * ((u1-u2) - sign*delta_u)
      ^^^^^^^^    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      stress       penalty correction (should be ≈0 for perfect solve)
```

Since MUMPS solves K*u = f exactly, u already satisfies `[u] ≈ delta_u` on fault faces
(enforced by the penalty in K). The penalty correction `penalty*(jump - delta_u)` should
be ≈ 0. If we drop it and use only `{σ(u)·n}`, we get the physical traction without
any cancellation.

**Risk**: Dropping the correction removes DG discretization error information.
For a perfect direct solve this should be fine. For iterative solvers it would not be.

**Implementation**: Add flag `--traction-stress-only` that skips the penalty correction
term in `ComputeTraction` for both interior and shared fault faces.

**Strategy 2: Weak-form traction (`--traction-weak-form`) — recommended fix**

Instead of re-evaluating the traction formula from scratch, reuse the SAME DG face
integrators that K uses. For each fault face:

1. Call `DGElasticityIntegrator::AssembleFaceMatrix` (consistency + symmetry) on the
   fault face → get element matrix `K_face`
2. Call `DGElasticityIPPenaltyIntegrator::AssembleFaceMatrix` (penalty) → get `P_face`
3. Compute face operator application: `(K_face + P_face) * u_local`
4. Subtract the slip RHS: `- f_slip_local`
5. The result is the traction contribution, algebraically consistent with K and f

This avoids cancellation because `K_face` and `P_face` are the SAME matrices used in
the global K assembly. The penalty in the operator and the penalty in the traction
come from identical code paths and identical numerical evaluations.

**Implementation**: Add flag `--traction-weak-form` that replaces the pointwise
traction evaluation with face-integrator-based evaluation.

**Strategy 3: penalty_factor tuning (workaround)**

Use `--penalty-factor 0.5` for p≥4 runs. Reduces cancellation enough for nucleation.
Not a proper fix but allows immediate production runs.

### 6.12 Test Plan

Two parallel tests on 2500m mesh at p=4:

| Job | Strategy | Flag | What it tests |
|-----|----------|------|---------------|
| **v50f** | Stress-only | `--traction-stress-only` | Is penalty correction the culprit? |
| **v50g** | Weak-form | `--traction-weak-form` | Does algebraic consistency fix p=4? |

**Decision tree:**
- v50f works → penalty correction is the problem; stress-only may be sufficient
- v50g works → algebraic consistency is the proper fix
- Both work → v50g preferred (more principled, includes DG error info)
- Both fail → root cause is elsewhere (psi init, friction law, etc.)

Also test both at p=2 to verify they don't break the working case:
- v50f should reproduce v50d (p=2 same mesh) results
- v50g should reproduce v50d results

---

## 7. Code Changes

### 7.1 `--penalty-factor` flag (v50a+)

Added `--penalty-factor <value>` command-line flag that scales the IP penalty uniformly in:
- Bilinear form K matrix (`DGElasticityIPPenaltyIntegrator`)
- RHS slip assembly (`AssembleSlipContributionIP`, `AssembleSlipContributionIPShared`)
- Traction computation (`ComputeTraction` IP path, interior and shared faces)

Both K and f are scaled by the same factor, maintaining K-f consistency.

### 7.2 Order-aware CFL beta (v50+)

```cpp
real_t c_N_1 = order * (order + dim - 1.0) / dim;
real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;
real_t beta = 4.0 * c_N_1 / c_N_1_ref;
```

### 7.3 `--traction-stress-only` flag (v50f)

Disables penalty correction in `ComputeTraction` IP path for both interior and shared
fault faces. The traction becomes `T_q = {σ(u)·n}` only. The penalty still appears in
K (stiffness matrix) and f (RHS) — only the traction recovery changes.

### 7.4 `--traction-weak-form` flag (v50g)

Replaces the pointwise traction evaluation in `ComputeTraction` with face-integrator-based
evaluation. For each fault face:
1. Get the face element matrix from `DGElasticityIntegrator` + `DGElasticityIPPenaltyIntegrator`
2. Apply to local displacement DOFs: `T_face = (K_face + P_face) * u_local - f_slip_local`
3. L2-project the result to per-DOF traction

This ensures the penalty in traction recovery uses the EXACT SAME numerical values as
the penalty in the stiffness matrix assembly.

---

## 8. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v49 | CFL analysis: identified cascade root cause, tested V-guard + dt fixes | Done |
| v49f-m | Confirmed fixes across all configurations (serial, parallel, 1000m, 2500m) | All STABLE |
| **v50** | **CFL dt + V-guard made automatic production defaults** | **Implemented** |
| **v50 prod** | **1000m p=2 first earthquake: V_peak=0.79 m/s (91% of Tandem)** | **SUCCESS** |
| v50+ | Order-aware CFL beta scaling: beta(p) = 4.0 * c_N_1(p) / c_N_1(2) | Implemented |
| **v50a** | **2500m p=4 IP: nucleation FAILS — V decays 0.01→0.0025 in 2s** | **BUG** |
| v50b | 4000m p=6: pending results | Pending |
| **v50c** | **2500m p=4 BR2: ALSO FAILS — V decays 0.01→0.00002** | **BUG** |
| **v50d** | **2500m p=2 IP: WORKS — V grows to 0.32 (nucleation healthy)** | **WORKS** |
| **v50e** | **2500m p=4 IP half-penalty: WORKS — V grows to 0.017** | **WORKS** |
| v50+ | Root cause: cancellation in post-processed traction (solve is monolithic, recovery is not) | ANALYZED |
| v50f | Strategy 1: stress-only traction test at p=4 | Planned |
| v50g | Strategy 2: weak-form traction test at p=4 | Planned |
