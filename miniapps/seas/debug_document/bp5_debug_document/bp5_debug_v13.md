# BP5 Debug v13: In-Depth Hidden Error Analysis

**Date**: 2026-03-13 (updated), 2026-03-12 (initial analysis)
**Status**: Fixes H1-H4 implemented & compiled; H7 added and implemented
**Previous**: v12 (BR2 face_int scaling fix with w_centroid)

---

## 1. Executive Summary

After reviewing 12 rounds of debugging, the current BP5 1000m BR2 implementation produces
reasonable earthquake cycle behavior (~8 events over ~750 years) but exhibits:

1. **Event timing mismatch** with Tandem benchmark data
2. **Interseismic output intervals** larger than SCEC specification (~0.3 yr vs spec's ~0.1 yr)
3. **IP method blowup** (traction runaway from step ~1, never exits nucleation)
4. **Dip-slip divergence** at depth stations (z=10-22km) growing over time
5. **Small-number digit cutoff** in log10(V_dip) and log10(state) columns

This document identifies **7 hidden errors** not caught in v1-v12, traces each to specific
code locations, and proposes targeted fixes. Fixes H1-H4 and H7 have been implemented.

---

## 2. Parameter Verification: MFEM vs SCEC Spec vs Tandem

### 2.1 Material Properties (All Match)

| Parameter | MFEM (`bp5_params.hpp`) | SCEC Spec | Tandem (`bp5.lua`) | Status |
|-----------|------------------------|-----------|--------------------|----|
| rho | 2670.0 kg/m^3 | 2670 kg/m^3 | 2.670 g/cm^3 | OK |
| cs | 3464.0 m/s | 3.464 km/s | 3.464 km/s | OK |
| nu | 0.25 | 0.25 | 0.25 | OK |
| mu | rho*cs^2 = 32.04 GPa | rho*cs^2 | cs^2*rho0 | OK |
| lambda | mu (for nu=0.25) | mu | 2*nu*mu/(1-2nu) = mu | OK |
| eta | mu/(2*cs) = 4.625 MPa*s/m | mu/(2*cs) | cs*rho0/2 | OK |

### 2.2 Friction Parameters (All Match)

| Parameter | MFEM | SCEC | Tandem | Status |
|-----------|------|------|--------|--------|
| V0 | 1e-6 m/s | 1e-6 m/s | 1e-6 m/s | OK |
| f0 | 0.6 | 0.6 | 0.6 | OK |
| b | 0.03 | 0.03 | 0.03 | OK |
| L0 (Dc) | 0.14 m | 0.14 m | 0.14 m | OK |
| L_nuc | 0.13 m | 0.13 m | 0.13 m | OK |
| a0 | 0.004 | 0.004 | 0.004 | OK |
| amax | 0.04 | 0.04 | 0.04 | OK |
| sigma_n | 25e6 Pa | 25 MPa | 25.0 MPa | OK |

### 2.3 Loading & Initial Conditions (DISCREPANCIES FOUND)

| Parameter | MFEM | SCEC Spec | Tandem | Status |
|-----------|------|-----------|--------|--------|
| Vp | 1e-9 m/s | 1e-9 m/s | 1e-9 m/s | OK |
| V_init | 1e-9 m/s | 1e-9 m/s | 1e-9 m/s | OK |
| V_zero | 1e-20 m/s | 1e-20 m/s | 1e-20 m/s | OK |
| **V_nuc** | **0.03 m/s** | **0.03 m/s** | **0.01 m/s** | **MISMATCH** |
| **delta_tau (QD)** | **eta*V_nuc** | **eta*V_i** | **None** | **MISMATCH** |

### 2.4 Geometric Parameters (All Match)

| Parameter | MFEM | SCEC | Tandem | Status |
|-----------|------|------|--------|--------|
| hs | 2 km | 2 km | 2 km | OK |
| ht | 2 km | 2 km | 2 km | OK |
| H | 12 km | 12 km | 12 km | OK |
| l_vw | 60 km | 60 km | 60 km | OK |
| w_nuc | 12 km | 12 km | 12 km | OK |
| Wf | 40 km | 40 km | 40 km | OK |
| lf | 100 km | 100 km | 100 km | OK |

### 2.5 Time Stepping Parameters

| Parameter | MFEM | Tandem | Status |
|-----------|------|--------|--------|
| Method | DOPRI5(4) | DOPRI5(4) via PETSc | OK |
| atol | 1e-7 | 1e-7 | OK |
| rtol | 1e-50 (disabled) | 1e-50 (disabled) | OK |
| Error norm | L-infinity | L-infinity | OK |
| dt_init | min(1e3, 0.01*L_nuc/V_init) | 0.0001 s | **Different** |
| **dt_max** | **0.1 year (was 0.5)** | **PETSc default** | **Fixed (H4)** |
| dt_min | 1e-6 s | PETSc default | OK |

---

## 3. Hidden Errors Identified

### H1: PsiToTheta Uses Global Dc Instead of Per-DOF Dc (Output Bug) -- FIXED

**Location**: `dieterich_ruina.hpp:261-264`

```cpp
real_t PsiToTheta(real_t psi) const
{
   return (cp_.Dc / cp_.V0) * std::exp((psi - cp_.f0) / cp_.b);
}
```

**Problem**: `cp_.Dc` is the global constant set at construction time:

```cpp
// bp5_verification_full.cpp:471
fc.Dc = params.L0;  // = 0.14
DieterichRuinaFriction friction(fc);
```

In BP5, Dc varies spatially: L0=0.14m outside nucleation, L_nuc=0.13m inside nucleation.
`PsiToTheta()` always uses `cp_.Dc = 0.14`, even for nucleation zone DOFs where Dc=0.13.

**Mathematical impact**:
- theta = (Dc/V0) * exp((psi - f0) / b)
- Using Dc=0.14 instead of 0.13 gives theta values that are 0.14/0.13 = **1.077x too large**
- log10(theta) error: log10(0.14/0.13) = **+0.033** (constant bias)
- This affects the "log10(state)" output column at nucleation zone stations

**Affected stations**: Those overlapping the nucleation zone:
- `fltst_strk-24dp+10` (x2=-24km, x3=10km) -- inside nucleation zone
- `fltst_strk-16dp+10` (x2=-16km, x3=10km) -- near nucleation zone boundary

**Does NOT affect dynamics**: PsiToTheta is only called for output (`GetTheta()`).
The time integration uses psi directly and calls `evolution_->Rate(V_abs, psi, Dc)` with
the correct per-DOF Dc from `Dc_values_(i)`.

**Fix**: Add a Dc parameter to PsiToTheta, or add a per-DOF version in `rate_state_fault.hpp`.

---

### H2: Output Interval Logic Causes Sparse Interseismic Output -- FIXED

**Location**: `bp5_benchmark_output.hpp:297-311`

```cpp
static real_t OutputInterval(real_t V_max)
{
   if (V_max > 1e-3)      return 0.001;                            // Coseismic
   else if (V_max > 1e-6) return 0.1;                              // Nucleation
   else                   return 0.01 * BP5Params::seconds_per_year; // Interseismic
}
```

**Problem**: The interseismic output interval is 0.01 year (315,576 s), but this is a
**minimum threshold**, not a guaranteed interval. Output is only written at accepted
ODE time steps. During deep interseismic, the adaptive RK45 takes steps up to
dt_max = 0.5 year, so the **effective** output interval is ~0.3-0.5 years.

**SCEC spec output requirements** (Section 4.1):
- Aseismic (max(V) < 1e-3 m/s): ~0.1 year
- Seismic (max(V) >= 1e-3 m/s): ~0.1 s

**Observed in MFEM results**: Interseismic output intervals alternate between ~0.085 yr
and ~0.202 yr (combined ~0.287 yr per pair). This is governed by the ODE step size, not
the 0.01-year threshold.

**Tandem comparison**: Tandem outputs every ~0.302 years during interseismic (controlled
by its adaptive output strategy with `t_max = 0.3 yr`).

**Impact**: The effective output frequency is comparable to Tandem's (~0.3 yr), but does
not match the SCEC spec's recommended ~0.1 yr. The "greater time interval" the user
observed is caused by dt_max being too large (0.5 yr), allowing the ODE solver to skip
past the output interval.

**Fix**: Reduce dt_max to 0.1 year (matching SCEC spec), or implement sub-stepping output
that writes at intermediate times within a large ODE step.

---

### H3: IP Method Penalty Lacks Elasticity Tensor Coupling (Blowup Root Cause) -- FIXED

**Location**: `elasticity_operator.hpp`, ComputeTraction() IP path

```cpp
if (method_ == DGMethod::IP)
{
   real_t kappa = (order_ + 1) * (order_ + 1);
   real_t detJ1 = FTr->Elem1->Weight();
   real_t detJ2 = FTr->Elem2->Weight();
   real_t nor_sq = nor * nor;
   real_t penalty = kappa * nor_sq * (1.0/(2.0*detJ1) + 1.0/(2.0*detJ2));

   for (int c = 0; c < dim; c++)
   {
      correction[c] = penalty * (u_jump[c] - sign * delta_u[c]);
   }
}
```

**Problem**: The IP traction correction applies a **scalar** penalty to each displacement
component independently: `correction[c] = penalty * jump[c]`. This is a "component-wise"
penalty that does NOT couple the displacement components through the elasticity tensor.

For 3D vector elasticity, the correct IP penalty should apply the 4th-order elasticity
tensor to the jump, giving a traction-like correction:

```
correction_i = penalty * C_{ij1k} * n_1 * n_j * jump_k / ||n||^2
```

The current scalar penalty effectively uses `correction_i = penalty * delta_ij * jump_j`,
which misses the lambda coupling between components and the mu shearing terms. This
creates an **inconsistency** between the bilinear form (which uses the full tensor in
`dg_elasticity_br2_integrator.hpp`) and the traction evaluation (which uses scalar penalty).

**Mathematical analysis**: For a fault face with normal n=(1,0,0) and Poisson ratio nu=0.25
(lambda=mu):
- Correct IP traction for jump_2 (strike): T_2 = mu * penalty * jump_2 (shear)
- Correct IP traction for jump_1 (normal): T_1 = (lambda+2mu) * penalty * jump_1
- Current code: T_c = penalty * jump_c (same scaling for ALL components)

The current penalty does not distinguish between shear and normal jumps, so it over-penalizes
shear jumps and under-penalizes normal jumps (for nu=0.25, it's ~2x wrong for normal).

**Why it causes blowup**: During the nucleation phase, the rapid slip creates large
displacement jumps. The scalar penalty fails to properly enforce the fault constraint,
causing traction oscillations that feed back into the friction law. The positive feedback
amplifies V_max from 0.039 to 1457 m/s in 608 steps without advancing past t=0.

**Contrast with BR2**: The BR2 method uses the full elasticity tensor coupling through
the `test_normal` operator:
```cpp
real_t tn = lambda_val_ * (u==s ? 1.0 : 0.0) * basis.normal[i]
   + mu_val_ * ((i==u ? 1.0 : 0.0) * basis.normal[s]
               + (i==s ? 1.0 : 0.0) * basis.normal[u]);
```
This is why BR2 works but IP does not.

**Fix**: Replace the scalar IP penalty with an elasticity-tensor-coupled penalty:
```cpp
for (int i = 0; i < dim; i++) {
   correction[i] = 0.0;
   for (int u = 0; u < dim; u++) {
      real_t tn = lambda_val_ * basis.normal[i] * basis.normal[u]
         + mu_val_ * ((i==u ? 1.0 : 0.0) * nor_sq_inv
                     + basis.normal[i] * basis.normal[u]);
      // Simplified: needs proper test_normal formulation
      correction[i] += penalty * tn * jump[u];
   }
}
```
Or more precisely, match the BR2 test_normal formulation but using the IP penalty
coefficient instead of the lifting operator.

---

### H4: dt_max = 0.5 Year Is Too Large -- FIXED

**Location**: `time_stepper.hpp:151`, `bp5_verification_full.cpp:568`

```cpp
dt_max_(0.5 * 3.15576e7),  // 0.5 year
// ...
ode_solver.SetDtMax(0.5 * BP5Params::seconds_per_year);
```

**Problem**: During deep interseismic (V_max ~ 1e-9 m/s), the adaptive RK45 increases dt
toward dt_max = 0.5 year (~1.58e7 s). This causes two issues:

1. **Output gaps**: Since output is only written at accepted time steps, the effective
   interseismic output interval becomes ~0.3-0.5 years, exceeding the SCEC spec's ~0.1 yr.

2. **Accumulated drift**: Large time steps during the interseismic phase may cause the
   error estimator to accumulate small drifts in the dip-slip component. The dip-slip
   dynamics are very weak (V_dip ~ 1e-20 m/s), so even small per-step errors can
   dominate over many cycles.

**Tandem comparison**: PETSc's TS adapter with DOPRI5(4) appears to take steps of ~0.3
years during interseismic (based on the ~0.302-year output intervals in the Tandem
benchmark data). This suggests PETSc's internal dt_max or error control limits steps
to ~0.3 years naturally.

**Fix**: Reduce dt_max to 0.1 year:
```cpp
ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);
```
This ensures output at ~0.1-year intervals (matching SCEC spec) and reduces per-step
error accumulation in weak dip-slip dynamics.

---

### H5: V_nuc Mismatch Explains Event Timing Discrepancy (Not a Bug)

**Location**: `bp5_params.hpp:112`

```cpp
real_t V_nuc = 0.03;  // SCEC spec
```

**vs Tandem** (`bp5.lua`):
```lua
return self.Vzero, 0.01   -- Tandem uses 0.01 m/s
```

**Impact on earthquake timing**: MFEM's nucleation perturbation is **3x stronger** than
Tandem's (V_nuc=0.03 vs 0.01). Additionally, MFEM adds `delta_tau = eta * V_nuc` in the
nucleation pre-stress (SCEC spec Eq. 23), while Tandem does NOT add this term.

The combined effect:
- delta_tau = eta * 0.03 = 4.625e6 * 0.03 = **138.75 kPa** overstress in MFEM
- Tandem has **zero** extra overstress (no delta_tau)
- MFEM's first earthquake nucleates **earlier** than Tandem's

This explains the event timing mismatch the user observed. It is NOT a bug -- MFEM
correctly follows the SCEC spec. The discrepancy is because the benchmark comparison
data comes from Tandem, which deviates from the spec.

**Recommendation**: For benchmark verification against the SCEC community, keep V_nuc=0.03
and delta_tau (SCEC spec). For comparison against Tandem specifically, a separate test run
with V_nuc=0.01 and no delta_tau would isolate other discrepancies from this parameter
difference. Could add a `--tandem-compat` CLI flag.

---

### H6: Pre-Stress Formula Differences (Not a Bug, Informational)

**MFEM** (`bp5_params.hpp:290-306`):
```cpp
real_t psi_ss = f0 + b * std::log(V0 / Vp);
real_t tau0_scalar = sigma_n * a * std::asinh((Vi_abs / (2.0 * V0)) * e)
                     + eta_val * Vi_abs;
if (IsNucleationZone(x2, x3))
   tau0_scalar += eta_val * Vi_abs;  // delta_tau = eta * V_i (QD)
```

**Tandem** (`bp5.lua`):
```lua
local tau0 = sn * ax * math.asinh((Vi2 / (2.0 * self.V0)) * e)
             + self:eta(x, y, z) * Vi2
-- No delta_tau added
```

Key differences:
1. MFEM uses `Vi_abs` (magnitude), Tandem uses `Vi2` (strike component only)
   - For non-nucleation: Vi_abs ≈ Vi2 (V_zero=1e-20 negligible) -- no practical difference
   - For nucleation: Vi_abs ≈ V_nuc ≈ Vi2 (V_zero negligible) -- no practical difference
2. MFEM adds delta_tau (per SCEC spec), Tandem does not
3. MFEM direction: positive (parallel to V), Tandem: negative (antiparallel to V)
   - Both are internally consistent within their own stress convention

**Conclusion**: The pre-stress formulas are equivalent modulo the delta_tau term and
sign convention (which is handled consistently). No fix needed.

---

### H7: Artificial 1e-30 Floor in SolveSlipRateVectorPsi -- FIXED

**Location**: `dieterich_ruina.hpp:394`

```cpp
if (tau_abs < 1e-30)
{
   V_vec[0] = 0.0;
   V_vec[1] = 0.0;
   ...
   return;
}
```

**Problem**: The friction solver has a hard cutoff at `tau_abs < 1e-30`, setting the
slip rate vector to exactly zero. This creates a discontinuity in the ODE right-hand
side: when `tau_abs` crosses the `1e-30` threshold, `V` jumps from exactly 0 to a small
positive value. While `1e-30` Pa is far below any physical stress, the discontinuity is
unnecessary and inconsistent with Tandem's approach.

**Tandem comparison**: Tandem's `DieterichRuinaAgeing.h:81-120` has **no floor at all**.
It uses the Brent bracket `[0, tauAbs/eta]` directly and computes the slip direction
as `-(V / tauAbs) * tauAbsVec` (line 119). The ratio `V/tauAbs` is well-defined for
any `tauAbs > 0` because:
- `V_abs` lies in `[0, tauAbs/eta]` from the Brent solver
- As `tauAbs → 0`, `V_abs → 0` and `V_abs/tauAbs → 1/(sigma_n * df/dV|_{V=0} + eta)`,
  which is finite: `df/dV|_{V=0} = a * exp(psi/a) / (2*V0)`.
- In practice, `tauAbs > 0` is guaranteed by the pre-stress `tau_pre`.

**Impact**: The `1e-30` floor does not visibly affect the current simulation results
(all physical stress values are far above `1e-30`), but it represents a code design
inconsistency with Tandem and an unnecessary discontinuity in the ODE RHS.

**Fix**: Change the guard from `tau_abs < 1e-30` to `tau_abs <= 0.0` (exact zero only):

```cpp
if (tau_abs <= 0.0)
{
   // Exactly zero traction -- no slip. Matches Tandem (no artificial floor).
   V_vec[0] = 0.0;
   V_vec[1] = 0.0;
   if (iterations) { *iterations = 0; }
   return;
}
```

**Note**: The separate `1e-30` floors in the output code (`bp5_benchmark_output.hpp:358-359`)
are retained -- these only prevent `log10(0) = -inf` in the output columns and do not
affect dynamics.

---

## 4. Root Cause Analysis: Dip-Slip Divergence at Depth

The comparison plots show that surface stations (z=0) match Tandem well, while depth
stations (z=10-22km) show progressive dip-slip divergence. Root causes:

1. **dt_max too large (H4)**: Small dip-slip errors (~1e-20 m/s * 0.5 year ≈ 1e-12 m)
   accumulate over 750 years with O(1000) time steps. Over 8 earthquake cycles, these
   compound to produce visible divergence in cumulative dip-slip.

2. **PsiToTheta Dc error (H1)**: Incorrect theta output at nucleation zone stations
   creates apparent divergence in the state variable comparison.

3. **V_nuc and delta_tau differences (H5, H6)**: Different nucleation strength causes
   different earthquake timing, which shifts the entire time series. Even if the dynamics
   are identical, a timing offset creates large apparent differences in interseismic
   quantities when plotted against absolute time.

4. **Possible DG traction noise at depth**: v10 identified mesh-transition DG traction
   noise in dip-slip. The smooth mesh grading (v10 fix) may not completely resolve this
   at 1000m resolution. Higher-order elements or finer resolution near the VW-VS
   transition could help.

---

## 5. Root Cause Analysis: IP Method Blowup

The IP method output (`bp5_1000m_ip_7597070.out`) shows:
- V_max grows monotonically from 0.039 to 1457 m/s in 608 steps
- Time never advances past t ≈ 0 years
- TRACTION BLOWUP begins at step 176, starting at Rank 50 DOF 15
- Traction magnitudes reach 6.7 GPa (vs normal stress of 25 MPa)

**Root cause**: H3 (scalar IP penalty for vector elasticity). The penalty correction
`correction[c] = penalty * jump[c]` does not couple displacement components through
the elasticity tensor. During nucleation:

1. Large strike-slip velocity → large displacement jump in x2
2. Scalar penalty applies same correction to all components
3. Missing lambda coupling means normal traction (x1 direction) is not properly
   balanced by the jump correction
4. Traction error at specific DOFs grows exponentially
5. Friction solver responds with unrealistic slip rates
6. Positive feedback loop: more slip → more traction error → more slip

The BR2 method avoids this because its lifting operator naturally includes the full
elasticity tensor coupling through the `test_normal` operator.

---

## 6. Output Interval Detailed Comparison

### MFEM Output Logic
```cpp
// Triggered when: time - last_write_time >= OutputInterval(V_max) * 0.99
// But output is only written at accepted ODE time steps.
// Effective interval = max(ODE_step_size, OutputInterval)
```

### Observed Intervals (from station strk+00dp+00)

**MFEM BR2 (successful run, ~750 years)**:
- Interseismic: alternating 0.085 yr / 0.202 yr (avg 0.143 yr/point)
- Total: 7,551 data points

**Tandem benchmark (1800 years)**:
- Interseismic: uniform ~0.302 yr
- Total: 36,057 data points

**SCEC spec recommendation**:
- Interseismic: ~0.1 yr
- Seismic: ~0.1 s
- Total: 10,000-100,000 time steps

The MFEM interseismic output is governed by the ODE step size (~0.3 yr at dt_max=0.5yr),
not the output interval threshold (0.01 yr). Reducing dt_max to 0.1 yr (fix H4) will
bring the output interval in line with the SCEC spec.

---

## 7. Implemented Fixes

### Fix H1: Per-DOF PsiToTheta -- IMPLEMENTED

**File**: `dieterich_ruina.hpp:265-269` -- Added overloaded PsiToTheta with Dc parameter:

```cpp
/// Convert psi to theta with per-DOF Dc: theta = (Dc/V0)*exp((psi - f0)/b).
real_t PsiToTheta(real_t psi, real_t Dc) const
{
   return (Dc / cp_.V0) * std::exp((psi - cp_.f0) / cp_.b);
}
```

**File**: `rate_state_fault.hpp:476-486` -- Use per-DOF Dc in GetTheta:

```cpp
void GetTheta(const Vector &state, Vector &theta) const
{
   // ...
   if (use_psi_)
   {
      for (int i = 0; i < num_nodes_; i++)
      {
         real_t psi = state(i * StatePerNode + PsiIndex);
         if constexpr (SlipComponents == 2)
         {
            // BP5: use per-DOF Dc for correct theta conversion
            real_t Dc = Dc_values_(i);
            theta(i) = dr_friction_->PsiToTheta(psi, Dc);
         }
         else
         {
            theta(i) = dr_friction_->PsiToTheta(psi);
         }
      }
   }
   // ...
}
```

### Fix H2 + H4: Reduce dt_max and Adjust Output Interval -- IMPLEMENTED

**File**: `bp5_verification_full.cpp:568`:
```cpp
ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);  // Was 0.5
```

**File**: `time_stepper.hpp:151` -- Also updated default dt_max from `0.5 * 3.15576e7` to `0.1 * 3.15576e7`.

**File**: `bp5_benchmark_output.hpp:297-311` -- Adjusted output intervals to SCEC spec:
```cpp
static real_t OutputInterval(real_t V_max)
{
   if (V_max > 1e-3)      return 0.1;                             // Coseismic: 0.1s (SCEC spec)
   else if (V_max > 1e-6) return 0.1;                             // Nucleation: 0.1s
   else                   return 0.1 * BP5Params::seconds_per_year; // Interseismic: 0.1 yr
}
```

Note: Coseismic changed from 0.001s to 0.1s per SCEC spec. The 0.001s interval generates
enormous output files (>100K points per earthquake) without adding useful information for
benchmark comparison. SCEC spec recommends ~0.1s.

### Fix H3: Elasticity-Tensor-Coupled IP Penalty -- IMPLEMENTED

**File**: `elasticity_operator.hpp`, ComputeTraction() IP path (both interior faces
and shared faces in parallel):

Replaced the scalar penalty with the full `test_normal` tensor-coupled penalty,
matching the BR2 formulation. The implemented code uses the exact bilinear form
assembly formula:

```cpp
if (method_ == DGMethod::IP)
{
   real_t kappa = (order_ + 1) * (order_ + 1);
   real_t detJ1 = FTr->Elem1->Weight();
   real_t detJ2 = FTr->Elem2->Weight();
   real_t nor_sq = nor * nor;
   real_t ip_coeff = kappa * nor_sq * (1.0/(2.0*detJ1) + 1.0/(2.0*detJ2));

   real_t jump[3];
   for (int c = 0; c < dim; c++)
      jump[c] = u_jump[c] - sign * delta_u[c];

   for (int i = 0; i < dim; i++)
   {
      correction[i] = 0.0;
      for (int u = 0; u < dim; u++)
      {
         for (int s = 0; s < dim; s++)
         {
            real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
               + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                           + (i == s ? 1.0 : 0.0) * basis.normal[u]);
            correction[i] += ip_coeff * tn * basis.normal[s] * jump[u];
         }
      }
   }
}
```

This reuses the same `test_normal` formulation as BR2, ensuring consistency between
the bilinear form and the traction evaluation. Applied to both the interior face path
(line ~1969) and the shared face path (line ~2220).

---

### Fix H7: Remove Artificial tau_abs Floor in SolveSlipRateVectorPsi -- IMPLEMENTED

**File**: `dieterich_ruina.hpp:394`

Changed `tau_abs < 1e-30` to `tau_abs <= 0.0` to match Tandem's approach (no artificial
floor). See H7 description above for rationale.

---

## 8. Implementation Status

| Fix | Issue | Status | Files Modified |
|-----|-------|--------|----------------|
| H1 | PsiToTheta per-DOF Dc | IMPLEMENTED | `dieterich_ruina.hpp`, `rate_state_fault.hpp` |
| H2 | Output intervals to SCEC spec | IMPLEMENTED | `bp5_benchmark_output.hpp` |
| H3 | IP tensor-coupled penalty | IMPLEMENTED | `elasticity_operator.hpp` (2 paths) |
| H4 | dt_max 0.5yr → 0.1yr | IMPLEMENTED | `bp5_verification_full.cpp`, `time_stepper.hpp` |
| H5 | V_nuc mismatch | No fix needed | Informational (SCEC vs Tandem) |
| H6 | Pre-stress formula | No fix needed | Informational |
| H7 | tau_abs 1e-30 floor | IMPLEMENTED | `dieterich_ruina.hpp` |

All code changes compile successfully (`make -j seas_bp5_full`).

H5 and H6 are informational; no code changes needed. For Tandem comparison,
consider adding a `--tandem-compat` flag that sets V_nuc=0.01 and disables delta_tau.

---

## 9. Summary of v12 + v13 Fix Status

The v12 fix (adding `w_centroid` to face_int and keeping `nor(s)` instead of
`basis.normal[s]`) has been applied to the current codebase. The blowup in the Mar 12
cluster runs was likely from jobs submitted with v11 code (before v12 was committed).
The current code (commit `2cffe7c`) should have the correct BR2 traction scaling.

v13 fixes (H1-H4, H7) have been implemented on top of v12:
- H1: `PsiToTheta(psi, Dc)` overload + per-DOF Dc in `GetTheta()` for BP5
- H2: Output intervals adjusted to SCEC spec (coseismic 0.1s, interseismic 0.1yr)
- H3: IP penalty now uses full elasticity tensor coupling (both interior + shared faces)
- H4: dt_max reduced from 0.5yr to 0.1yr (time_stepper default + bp5_verification_full)
- H7: `tau_abs < 1e-30` guard changed to `tau_abs <= 0.0` (matches Tandem, no artificial floor)

**Recommendation**: Submit a new BR2 run with all v13 fixes to verify improvements
in output accuracy, output intervals, and dip-slip drift. Then test IP method to
verify H3 fix resolves the traction blowup.

---

## 10. Verification Plan

After implementing fixes:

1. **Smoke test** (inline mesh, 1 year): Verify no NaN/blowup, check output format
2. **Short run** (inline mesh, 100 years): Verify earthquake nucleation and output intervals
3. **Full run** (1000m Gmsh mesh, 1800 years, BR2): Compare with Tandem benchmark data
4. **IP comparison** (1000m Gmsh mesh, 100 years, IP): Verify IP no longer blows up
5. **Tandem-compat run** (V_nuc=0.01, no delta_tau): Isolate parameter effects from code bugs
