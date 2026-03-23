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

### 6.11 Fix Strategy 1: Stress-Only Traction — TESTED, DISPROVED

**Hypothesis**: The penalty correction in ComputeTraction is the source of cancellation
error. Dropping it should give accurate traction from `{σ(u)·n}` alone, since the
MUMPS direct solve enforces `[u] ≈ delta_u` through the penalty in K.

**Implementation**: Added `--traction-stress-only` flag that sets `correction_q = 0`
in both interior and shared fault face IP traction paths.

**Results** (v50f, v50f2):

| Test | Order | Traction | V trajectory | Status |
|------|-------|----------|-------------|--------|
| **v50f** | p=4 | stress-only | 0.01 → **10.2 m/s** in 35 steps | **BLOWUP** |
| **v50f2** | p=2 | stress-only | 0.01 → **46.5 m/s** in 397 steps | **BLOWUP** |

**Both blow up — at ALL orders.** The penalty correction is ESSENTIAL for stability,
not just a cancellation artifact. Without it, the DG traction is physically wrong.

**Why**: The penalty correction `penalty * sign * ((u1-u2) - sign*delta_u)` is NOT
approximately zero. Even with a direct solver, the DG discretization enforces
`[u] ≈ delta_u` only weakly (through the penalty in K). The residual
`(u1-u2) - delta_u` is O(h^p / penalty) — small but nonzero. The penalty correction
multiplies this residual by the penalty, giving an O(h^p) contribution that is
essential for the DG traction formula to be accurate.

Removing the penalty correction removes the DG flux consistency, causing the
displacement-to-traction mapping to diverge.

**Conclusion**: Strategy 1 is INVALID. The penalty correction is structurally
required in the DG traction formula. The fix must KEEP the penalty correction
but compute it CONSISTENTLY with K.

### 6.12 Revised Root Cause Understanding

The v50f/f2 blowup disproves the "catastrophic cancellation" hypothesis in its
original form. The corrected understanding:

1. **The penalty correction is NOT approximately zero** — it carries essential DG
   flux information. Removing it breaks the formulation entirely.

2. **At p=2, the penalty correction is computed accurately enough** by both K assembly
   and ComputeTraction. The small inconsistency between the two code paths is
   tolerable because the penalty magnitude is moderate (c_N_1 = 2.67).

3. **At p=4, the penalty is 3× larger** (c_N_1 = 8.0). The inconsistency between
   K's penalty evaluation and ComputeTraction's penalty evaluation is amplified 3×.
   This amplified inconsistency produces traction errors large enough to decelerate
   nucleation.

4. **The fix is NOT to remove the penalty correction** (that breaks everything), but
   to ensure the penalty correction in ComputeTraction uses the EXACT SAME numerical
   evaluation as the penalty in K. This is Strategy 2.

### 6.13 Fix Strategy 2: Weak-Form / Consistent-Penalty Traction — NEXT

Instead of re-evaluating the DG traction formula from scratch in ComputeTraction,
reuse the SAME integrators that assembled K. For each fault face:

1. Call `DGElasticityIntegrator::AssembleFaceMatrix` (consistency + symmetry) → `K_cs`
2. Call `DGElasticityIPPenaltyIntegrator::AssembleFaceMatrix` (penalty) → `K_pen`
3. Compute operator application: `(K_cs + K_pen) * u_local`
4. Subtract the slip RHS for this face: `- f_slip_face`
5. Result = traction contribution, algebraically consistent with K and f

This keeps the penalty correction (needed for stability) while ensuring it is computed
by the SAME code path as K (avoiding the inconsistency that kills p=4).

**Why this should work:**
- v50e (half-penalty) works → reducing the inconsistency helps
- v50f (no penalty) blows up → the correction is needed
- Strategy 2 gives FULL penalty with ZERO inconsistency → best of both worlds

### 6.14 Strategy 2 ABANDONED — Penalty Formulas are Algebraically Identical

**Critical finding**: Before implementing Strategy 2, line-by-line analysis of both
code paths revealed that for flat faces with constant material properties:

**K assembly** (`DGElasticityIPPenaltyIntegrator::AssembleFaceMatrix`):
```cpp
penalty = penalty_factor_ * (p0 + p1) / 4.0;              // line 129
coeff = penalty * ip.weight * nl_q;                         // line 172
elmat(i,j) += coeff * shape1(i) * shape1(j);               // line 184
```

**ComputeTraction** (`elasticity_operator.hpp`):
```cpp
penalty_ip = penalty_factor_ * (p0 + p1) / 4.0;           // line 3337
correction_q[c] = -penalty_ip * sign * jump_c;              // line 3468
```

Both formulas:
- Use **identical** penalty computation: `(D+1)*c_N_1*(D*nl_q/vol)*(c1²/c0)/4`
- Use **identical** shape functions: same MFEM `CalcShape` call
- Use **identical** displacement interpolation: `u_q = Σ shape(k) * u(k)`
- Use **identical** material constants: `lambda_coeff_` / `mu_coeff_` vs `lambda_val_` / `mu_val_`
  (both are `ConstantCoefficient`, giving identical numerical values)

For flat faces on linear tets:
- `nl_q` = constant across face (face Jacobian is constant) → identical at any quad point
- `vol1`, `vol2` = constant per element → identical
- Penalty is face-constant → identical regardless of quadrature order (2p vs 2p+1)

**Conclusion: Strategy 2 would produce NUMERICALLY IDENTICAL results to the current
code. Implementing it would be wasted effort.**

The "K-traction inconsistency" hypothesis is **WRONG** for our geometry. The penalty
correction in ComputeTraction IS algebraically consistent with K.

### 6.15 Revised Root Cause Hypothesis: Solver Accuracy × Penalty Amplification

Since the penalty formulas are identical, the only remaining source of error is the
**solve accuracy**. The displacement u from K*u = f has solve error ε that depends on:
- MUMPS BLR tolerance (currently `1e-12`)
- System condition number (worse at higher p)

The residual jump on fault faces is:
```
[u] - delta_u = O(ε)
```

The penalty correction in the traction is:
```
penalty * ([u] - delta_u) = penalty * O(ε)
```

At p=4, penalty is 3× larger AND ε is likely larger (worse conditioning). The product
`penalty * ε` grows faster than linearly with p, potentially producing a correction
that dominates the friction weakening signal.

**Key insight**: Both the penalty traction change and friction weakening scale linearly
with dt, so their ratio is dt-INDEPENDENT:
```
K_pen / K_crit = penalty / [(b-a)*σ_n/Dc]
               ≈ 2e10 / 5e6 = 4000  (p=4)
               ≈ 6.7e9 / 5e6 = 1340 (p=2)
```

Both are >>1, yet p=2 works. This is because the penalty in K ENFORCES `[u] ≈ delta_u`,
making the penalty correction ≈ 0 for an accurate solve. **The net traction is just
`{σ·n}`, regardless of penalty magnitude — IF the solve is accurate enough.**

At p=4, the BLR approximate factorization may not be accurate enough. The residual
`[u] - delta_u` is larger, and when multiplied by the 3× stronger penalty, produces
a correction that opposes nucleation.

**Evidence supporting this hypothesis:**
- Half-penalty works → halving the penalty halves the amplification of solve error
- Full penalty fails → the amplification exceeds the tolerance
- Stress-only blows up → the correction IS needed (it's not zero), but its SIGN
  should be determined by the physics, not the solver error

### 6.16 Hypothesis 7: GaussLobatto vs WarpAndBlend Node Distribution (p≥3)

**Background**: From v46 Section 3.3 Factor 4 and v47 Section 18.3, a known
difference between SEAS-MFEM and Tandem was identified but never tested:

| | SEAS-MFEM | Tandem |
|--|-----------|--------|
| Face DOF nodes (p≤2) | GaussLobatto | WarpAndBlend | **Identical** |
| Face DOF nodes (p≥3) | **GaussLobatto** | **WarpAndBlend** | **DIFFER** |

*Source*: `face_quadrature.hpp:58` — `H1_TriangleElement(face_order, BasisType::GaussLobatto)`
*Tandem*: `RateAndStateBase.cpp:10-13` — `NodalRefElement<2>(PolynomialDegree, WarpAndBlendFactory<2>())`

At p=4 (15 nodes on triangle):
- **GaussLobatto**: Nodes cluster at edges/vertices (tensor-product mapped to simplex).
  3 vertices + 9 edge nodes (3 per edge) + 3 interior nodes.
- **WarpAndBlend**: Nodes distributed more evenly, optimized for simplex interpolation.
  ~6 interior nodes, lower Lebesgue constant.

**How this could cause p=4 failure:**

1. **Mass matrix conditioning**: `GalerkinProject` uses `M_ref_inv` to project
   quad-point traction to face DOFs. Poorly conditioned M_ref amplifies small
   errors in the traction integral, potentially flipping the sign of individual
   DOF traction values.

2. **Spatial resolution**: GaussLobatto clusters DOFs at face edges, under-resolving
   the interior traction. WarpAndBlend provides more interior DOFs, better capturing
   the traction variation across the face.

3. **Interpolation stability**: The Lebesgue constant for GaussLobatto on triangles
   is higher than WarpAndBlend at p≥3. Higher Lebesgue constant means the
   interpolation `InterpolateToQuadPoints` amplifies errors more.

**Why this explains the observations:**
- At p≤2: nodes are identical → no effect (p=2 works) ✓
- At p=4: nodes diverge → GaussLobatto conditioning could cause traction errors ✓
- Half-penalty works → reduces the traction VALUE that the poorly-conditioned
  M_ref_inv operates on, keeping errors below the nucleation threshold ✓
- BR2 also fails → BR2 uses the SAME FaceQuadrature with the same nodes ✓

### 6.17 Phase 0 Result: SMOKING GUN — GaussLobatto Catastrophically Ill-Conditioned

Computed mass matrix condition number `cond(M_ref)` and Lebesgue constant on the
reference triangle for GaussLobatto (our code) vs ClosedUniform (equispaced baseline):

```
=====================================================================================
Phase 0: Mass Matrix Conditioning — Reference Triangle
=====================================================================================

  p=2  nbf=6
  Basis Type                cond(M)   ||M^-1||_2   Lebesgue      cond(V)
  GaussLobatto                37.96        97.25      4.062        43.06
  ClosedUniform               17.21        96.40      1.667        30.97

  p=3  nbf=10
  GaussLobatto               219.63       229.13     13.342       549.96
  ClosedUniform               33.97       224.07      2.258       312.45

  p=4  nbf=15
  GaussLobatto              2900.92       524.50     46.022      7331.02
  ClosedUniform               57.74       513.79      3.470      3428.96

  p=6  nbf=28
  GaussLobatto            479657.75      1694.61    531.559   1417706.79
  ClosedUniform              214.87      1694.94      8.726    470039.02
```

**Summary table:**

| p | GL cond(M) | GL Lebesgue | Equi cond(M) | Equi Lebesgue | **cond(M) ratio** |
|---|-----------|-------------|--------------|---------------|-------------------|
| 2 | 38 | 4.1 | 17 | 1.7 | 2.2× |
| 3 | 220 | 13.3 | 34 | 2.3 | **6.5×** |
| **4** | **2,901** | **46.0** | 58 | 3.5 | **50×** |
| 6 | **479,658** | **531.6** | 215 | 8.7 | **2,232×** |

**Key findings:**

1. **GaussLobatto cond(M) explodes**: 38 → 220 → 2,901 → 479,658 from p=2 to p=6.
   ClosedUniform grows gently: 17 → 34 → 58 → 215.

2. **Lebesgue constant is catastrophic**: At p=4, GaussLobatto Lebesgue = **46**.
   This means ANY small error in the traction quadrature integral is amplified
   **46×** by the L2 projection `M_ref_inv * ∫T·φ dS`. At p=6: amplified **532×**.
   ClosedUniform at p=4: only 3.5× amplification.

3. **The root cause is the collapsed tensor-product mapping**: MFEM's
   `H1_TriangleElement(p, GaussLobatto)` maps 1D GaussLobatto points to the
   triangle via a Duffy-type collapse. This creates severely clustered nodes at
   the collapsed vertex, destroying interpolation properties on the simplex.

4. **WarpAndBlend (Tandem's choice) has Lebesgue ~ 2-3 at p=4** on triangles
   (Warburton 2006), similar to equispaced but with better edge properties.
   Tandem avoids this problem entirely.

**Why this is the root cause:**

The penalty correction in the traction is:
```
correction_q = penalty * ((u1_q - u2_q) - sign * delta_u_q)
```

This is O(h^p) — small but nonzero. When projected to face DOFs via GalerkinProject:
```
correction_DOF = M_ref_inv * Σ_q w_q * correction_q * φ_q
```

With `cond(M_ref) = 2901` and `Lebesgue = 46`, the projected correction at
individual DOFs can be **46× larger** than the quad-point values, with **wrong sign**
at some DOFs. This sign-flipped correction opposes the friction weakening, decelerating
nucleation.

At p=2 (`cond(M) = 38`, `Lebesgue = 4.1`), the amplification is manageable — the
projected correction maintains the correct sign at all DOFs.

**Why half-penalty works:** Reducing penalty by 0.5× halves the correction magnitude.
With half the signal going through the 46× amplifier, the amplified error is 23×
instead of 46× — apparently below the threshold where DOF signs flip.

**Why BR2 also fails:** BR2 uses the SAME `FaceQuadrature` with the same GaussLobatto
nodes and the same `GalerkinProject`. The ill-conditioned M_ref_inv amplifies the
BR2 lifting correction just as it amplifies the IP penalty correction.

### 6.18 Combined Test Plan: Node Distribution + Solver Accuracy

**Phase 0: COMPLETED** — conditioning analysis confirms 50× worse mass matrix at p=4
with GaussLobatto vs equispaced. Lebesgue constant 46 vs 3.5.

**Phase 1: Cluster tests**

| Job | Order | Change | What it tests |
|-----|-------|--------|---------------|
| **v50g** | p=4 | ClosedGL nodes | Does better node distribution fix p=4? |
| **v50g2** | p=4 | BLR tol=1e-14 | Does tighter solver fix p=4? |
| **v50g3** | p=2 | ClosedGL nodes | Regression check (should match v50d) |

**Phase 2: WarpAndBlend (if Phase 1 confirms)**

If ClosedGL helps or conditioning analysis shows large gap, implement WarpAndBlend
nodes to match Tandem exactly.

**Decision tree:**
- v50g (ClosedGL) nucleates → **node distribution is the root cause**
- v50g2 (BLR 1e-14) nucleates → **solver accuracy is the root cause**
- Both nucleate → both contribute; use both fixes
- Neither nucleates → something else; investigate further

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
| v50+ | Root cause hypothesis: inconsistent penalty between K assembly and traction post-processing | ANALYZED |
| **v50f** | **Strategy 1: stress-only traction p=4 → BLOWUP (V=10 m/s in 35 steps)** | **DISPROVED** |
| **v50f2** | **Strategy 1: stress-only traction p=2 → BLOWUP (V=46 m/s in 397 steps)** | **DISPROVED** |
| v50f+ | Revised understanding: penalty correction is essential, not cancellation artifact | ANALYZED |
| v50f++ | **Strategy 2 ABANDONED**: penalty formulas are algebraically identical for flat faces — no inconsistency to fix | **DISPROVED** |
| v50f++ | **New hypothesis**: solver accuracy (MUMPS BLR) × penalty amplification at p=4 | **ANALYZING** |
| v50f++ | Penalty formulas algebraically identical → Strategy 2 abandoned | DISPROVED |
| v50f++ | Two new hypotheses: (1) solver accuracy, (2) GaussLobatto vs WarpAndBlend nodes (v46/v47 flagged, never tested) | ANALYZING |
| v50g | **Phase 0: SMOKING GUN — GL cond(M)=2901 vs Equi=58 at p=4 (50× worse). Lebesgue=46 vs 3.5 (13×)** | **CONFIRMED** |
| v50g | Phase 1: ClosedUniform nodes at p=4 + p=2 regression check | Done |
| **v50g test** | **p=4 ClosedUniform 200 ranks: V grows 0.01→0.012 — NUCLEATION WORKS** | **FIXED** |
| v50g2 test | p=2 ClosedUniform 200 ranks: identical to GL p=2 (ratio=1.000000) — no regression | **✓** |
| **v50g p=4 prod** | **p=4 ClosedUniform 400 ranks: V decays 0.01→0.00017 — FAILS at 400 ranks** | **PARALLEL BUG** |
| **v50g p=6 prod** | **p=6 ClosedUniform 400 ranks 4000m: V grows 0.01→0.033 — NUCLEATION WORKS** | **✓ BEST MATCH** |

---

## 9. Phase 1 Results: ClosedUniform Node Distribution

### 9.1 Test Matrix Results

| Run | Mesh | Order | Nodes | Ranks | V trajectory | Status |
|-----|------|-------|-------|-------|-------------|--------|
| v50g test | 2500m | p=4 | **Equi** | 200 | 0.01 → **0.012** ↑ | **WORKS** ✓ |
| v50g2 test | 2500m | p=2 | **Equi** | 200 | Identical to GL | **No regression** |
| v50g p=4 prod | 2500m | p=4 | **Equi** | 400 | 0.01 → **0.00017** ↓ | **FAILS** ✗ |
| v50g p=6 prod | 4000m | p=6 | **Equi** | 400 | 0.01 → **0.033** ↑ | **WORKS** ✓ |

### 9.2 Key Finding: ClosedUniform Fixes p=4 Serial, Fails p=4 Parallel

The node distribution fix resolves the mass matrix conditioning issue at p=4 (V grows
at 200 ranks), confirming the root cause. However, at 400 ranks the p=4 nucleation
still fails — a SEPARATE parallel-specific issue exists.

The p=6 run at 400 ranks WORKS, which means the parallel issue is NOT universal for
high orders. The difference is the mesh: 2500m mesh with 400 ranks has very few
elements per rank (~12-15), creating a high ratio of shared faces to interior faces.
The 4000m mesh with fewer total elements may have a different partition structure.

### 9.3 p=6 4000m vs Tandem p6: Excellent Match

Direct comparison on Tandem's own 4000m mesh, same polynomial order p=6:

**Initial conditions: PERFECT**
```
  Our:    V=0.010000, tau=21.1481 MPa, state(log10)=8.113943
  Tandem: V=0.010000, tau=21.1481 MPa, state(log10)=8.113943
  tau diff: 0.000005 MPa (6 significant digits)
```

**Nucleation evolution at (x2=-24, x3=10):**

| t (s) | Our V (m/s) | Tandem V | Ratio | Our slip | Tandem slip | Our tau | Tandem tau |
|-------|-------------|----------|-------|----------|-------------|---------|-----------|
| 0 | 0.01000 | 0.01000 | 1.00× | 0 | 0 | 21.15 | 21.15 |
| 3.0 | 0.01539 | 0.01497 | 1.03× | 0.039 | 0.039 | 20.99 | 20.99 |
| 8.9 | 0.01892 | 0.01782 | 1.06× | 0.142 | 0.139 | 20.43 | 20.44 |
| 16.7 | 0.02313 | 0.02165 | 1.07× | 0.301 | 0.290 | 19.55 | 19.61 |
| 23.4 | **0.03298** | **0.02949** | **1.12×** | 0.489 | 0.456 | 18.56 | 18.71 |

**Event-aligned timing: OUR nucleation is AHEAD of Tandem:**
- V = 0.016 m/s: Tandem 4.3s → Ours **3.4s** (1.0s faster)
- V = 0.032 m/s: Tandem 24.6s → Ours **22.8s** (1.8s faster)

**Non-nucleation stations at t=23.4s:**
- All tau values match to 0.01 MPa or better
- All V at plate rate (identical)
- Surface tau: our 13.30 vs Tandem 13.30 (exact)

**Assessment**: This is the **best Tandem match achieved** — 3-12% V agreement, <0.15 MPa
tau agreement, 1-2 second timing lead (not lag). The slight V excess may be from
ClosedUniform vs WarpAndBlend nodes or the `dim*` penalty correction.

Tandem's peak V at this station is 0.42 m/s at t=42s. Extrapolating our growth rate,
we should reach this around t=40-42s — nearly matching Tandem's timing.

### 9.4 1000m p=2 Production — Healing Analysis

The 1000m p=2 run reached t≈70.5s. The first earthquake is in the healing phase:

| Milestone | Tandem (s) | Our extrapolated (s) | Delay |
|-----------|-----------|---------------------|-------|
| V < 0.1 | 63.7 | 72.4 | +8.7s |
| V < 0.05 | 75.2 | 75.3 | +0.1s |
| V < 0.01 | 77.5 | 77.6 | +0.1s |
| V < 0.001 | 79.7 | 78.1 | -1.6s |

The healing RATE matches Tandem precisely: -0.0172 m/s² (ours) vs -0.0167 m/s² (Tandem)
at V=0.13. The ~8s nucleation delay carries through but the physics of healing is correct.
~4400 more steps needed to complete healing (dt ≈ 0.0016s during coseismic).

### 9.5 Remaining Issues

1. **p=4 parallel failure at 400 ranks** — ClosedUniform fixes serial but not 400-rank
   parallel on 2500m mesh. The high shared-face ratio (400 ranks ÷ ~5000 elements =
   ~12 elements/rank) may trigger shared face coupling issues. Test at 200 ranks
   for production, or investigate the shared face sign bug from v48.

2. **p=6 production run still early** — at t=23.4s, needs to reach t≈42s for first
   earthquake peak. Running on 48-hour allocation should be sufficient.

3. **1000m p=2 production** — in healing phase, needs ~4400 more steps to complete
   first earthquake and enter interseismic.
