# BP5 Debug v46: Multi-DOF Nucleation Zone — Tandem Comparison & Revert

## Problem

The v46b run (p=2, IP, multi-DOF) produces an immediate earthquake within ~15 seconds of t=0, while:
- p=1 with identical initialization does NOT produce an immediate earthquake
- Tandem reference shows first earthquake at ~250 years

The immediate earthquake creates ~5.7m of strike slip at the nucleation station, causing a permanent offset in strike-direction comparisons with Tandem.

## Root Cause Analysis

**Two interacting issues:**

### 1. SCEC psi initialization (genuine overstress)

The v46b run uses the default `psi_init_mode = "scec"` (not `--psi-init tandem`). With SCEC initialization:
- psi0 = f0 + b*ln(V0/V_init) = 0.8072 everywhere (same as psi_ss at plate rate)
- In the nucleation zone, tau_pre includes delta_tau = eta * V_nuc = 138.7 kPa
- The delta_tau is genuine overstress (NOT absorbed into psi)
- This produces V = 0.049 m/s at nucleation DOFs at t=0 (not V_nuc = 0.03)

### 2. Within-face DOF discontinuity at nucleation boundary (p≥2 only)

At p=2, each face has nbf=6 DOFs at GaussLobatto nodes. Faces straddling the nucleation zone boundary have:
- **Inside DOFs**: V_init = V_nuc = 0.03, tau_pre with delta_tau, Dc = 0.13
- **Outside DOFs**: V_init = 1e-9, tau_pre without delta_tau, Dc = 0.14

The V_init ratio across the boundary is 3×10^7. When this discontinuity is interpolated to quad points via the polynomial basis, the IP penalty (which scales as p²/h) amplifies it into stress concentrations that trigger runaway nucleation.

**At p=1** (nbf=1): only one DOF per face, no within-face variation. The nucleation zone boundary is handled face-by-face with smooth jumps. The 0.65% overstress decays via radiation damping without triggering rupture:
- Max V at nucleation station: 0.052 m/s (peaks at t=0.2s, then decays)

**At p=2** (nbf=6): within-face polynomial interpolation of the 7-order-of-magnitude V gradient creates stress concentrations that compound the 0.65% overstress:
- Max V at nucleation station: 0.514 m/s at t=14.8s (full coseismic earthquake)

## Evidence

Data comparison at nucleation station (x2=-24km, x3=10km):

| Quantity | p=1 (v44d) | p=2 (v46b) |
|----------|-----------|-----------|
| V at t=0 | 0.049 m/s | 0.049 m/s |
| Max V | 0.052 m/s (t=0.2s) | 0.514 m/s (t=14.8s) |
| Earthquake? | No | Yes (within 15s) |
| Final slip (1yr) | 0.87 m | 5.72 m |

Both runs use identical SCEC initialization. The difference is entirely due to multi-DOF.

---

## Attempted Fix: Face-Centroid Nucleation Zone Detection (REVERTED)

An initial fix was attempted: determine nucleation zone status at the face centroid, then apply consistent parameters to all DOFs on that face:

- **`a(x2,x3)`**: Per-DOF (smoothly varying with 2km transition) — unchanged
- **`IsNucleationZone`**: Per-face centroid — prevents within-face discontinuity
- **`V_init`, `delta_tau`, `Dc`**: Based on face centroid — uniform within face

### Why it was reverted

A detailed comparison with Tandem's implementation revealed that the face-centroid approach diverges from how Tandem handles this problem. **Tandem uses per-DOF parameter evaluation at all polynomial orders**, and avoids instability through a fundamentally different initialization approach (see below).

---

## Tandem's Approach: Per-DOF Evaluation + Equilibrium Initialization

### Parameter Evaluation: Per-DOF (NOT Per-Face-Centroid)

Tandem evaluates ALL fault parameters (a, L/Dc, V_init, tau_pre, eta, sigma_n) at each individual DOF's physical coordinates.

**Source**: `tandem/app/localoperator/RateAndState.h:36-43`
```cpp
void set_params(param_fun_t pfun) {
    auto num_nodes = fault_.storage().size();  // = numFaces * nbf
    law_.set_num_nodes(num_nodes);
    for (std::size_t index = 0; index < num_nodes; ++index) {
        auto params = pfun(fault_.storage()[index].template get<Coords>());
        //                 ↑ Each DOF's own physical (x,y,z) coordinates
        law_.set_params(index, params);
    }
}
```

This means on a face straddling the nucleation boundary, different DOFs on the same face CAN get different L, V_init, and tau_pre values. Tandem explicitly allows within-face parameter discontinuities.

### Nucleation Zone Boundary: Epsilon Tolerance

Tandem uses an epsilon-based tolerance (`eps=1e-3`, i.e., 1mm) to slightly expand the nucleation zone boundary.

**Source**: `tandem/examples/tandem/3d/bp5.lua:49-57, 105-107`
```lua
function BP5:in_nucleation(x, y, z)
    local d = -z
    local s = x
    local eps = self.eps  -- 1e-3 for bp5_outside
    if self.h_s + self.h_t <= d+eps and d-eps <= self.h_s + self.h_t + self.H
       and -self.l/2.0 <= s+eps and s-eps <= -self.l/2.0 + self.w then
        return true
    end
    return false
end

bp5_exact = BP5.new({eps=0.0})     -- exact boundary
bp5_outside = BP5.new({eps=1e-3})  -- expand by 1mm (DEFAULT)
bp5_inside = BP5.new({eps=-1e-3})  -- shrink by 1mm
```

Default config (`bp5.toml:4`): `scenario = "bp5_outside"` — the nucleation zone is expanded by 1mm in all directions so DOFs within 1mm outside the nominal boundary are treated as inside.

### Initialization: NO delta_tau, Equilibrium psi_init

The critical difference: Tandem's initialization creates exact stress equilibrium at every DOF.

| Aspect | Tandem | Our Code (SCEC default) |
|--------|--------|-------------------------|
| V_nuc | 0.01 m/s | 0.03 m/s |
| delta_tau | **Not defined** (bp5.lua has no delta_tau function) | eta * V_nuc = 138.7 kPa |
| tau_pre formula | `sn*a*asinh(Vi/(2V0)*exp(psi_ss/a)) + eta*Vi` | Same + delta_tau in nucleation zone |
| psi_init | Computed from total stress at t=0 → exact equilibrium | Fixed: f0 + b*ln(V0/V_init) → does NOT absorb delta_tau |
| Initial V at nucleation | Exactly V_nuc = 0.01 m/s | 0.049 m/s (overstressed) |
| Within-face V ratio | 1e7 (0.01 vs 1e-9) | 3e7 (0.03 vs 1e-9) |

**Tandem's tau_pre** (`bp5.lua:94-103`):
```lua
function BP5:tau_pre(x, y, z)
    local Vi1, Vi2 = self:Vinit(x, y, z)
    local Vi = math.sqrt(Vi1^2 + Vi2^2)
    local sn = self:sn_pre(x, y, z)
    local ax = self:a(x, y, z)
    local e = math.exp((self.f0 + self.b * math.log(self.V0 / self.Vp)) / ax)
    local tau0 = sn * ax * math.asinh((Vi2 / (2.0 * self.V0)) * e)
                 + self:eta(x, y, z) * Vi2
    -- NOTE: NO delta_tau added. eta*Vi2 is equilibrium radiation damping.
    return -tau0 * Vi1 / Vi, -tau0 * Vi2 / Vi
end
```

**Tandem's psi_init** (`DieterichRuinaAgeing.h:51-63`):
```cpp
double psi_init(std::size_t index, double sn,
                std::array<double, TangentialComponents> const& tau) const {
    double snAbs = -sn + p_[index].get<SnPre>();
    double tauAbs = norm(tau + p_[index].get<TauPre>());
    auto Vi = norm(p_[index].get<Vinit>());
    auto a = p_[index].get<A>();
    auto eta = p_[index].get<Eta>();
    double s = sinh((tauAbs - eta * Vi) / (a * snAbs));
    return a * log((2.0 * cp_.V0 / Vi) * s);
    // psi computed from total stress → exact equilibrium with Vi
}
```

**Mathematical proof of equilibrium**: With zero slip at t=0, elastic traction tau ≈ 0. Then:
```
tauAbs = |tau_pre| = sn*a*asinh(Vi/(2V0)*exp(psi_ss/a)) + eta*Vi
tauAbs - eta*Vi = sn*a*asinh(Vi/(2V0)*exp(psi_ss/a))
sinh((tauAbs - eta*Vi)/(a*sn)) = Vi/(2V0)*exp(psi_ss/a)
2V0/Vi * sinh(...) = exp(psi_ss/a)
psi_init = a * log(exp(psi_ss/a)) = psi_ss
```

Result: psi_init = psi_ss = f0 + b*ln(V0/Vp) everywhere, and V = Vi at every DOF. **No transient.**

### Why Tandem Can Allow Within-Face Discontinuities

1. **No delta_tau overstress**: No sharp 138.7 kPa stress jump at the nucleation boundary. The stress field varies smoothly because tau_pre is computed from Vi using the smooth asinh function.

2. **psi_init absorbs stress**: psi is computed from the actual total stress at each DOF. Even with different Vi on different DOFs of the same face, psi adjusts so V = Vi everywhere → exact equilibrium, no runaway.

3. **Smaller V_nuc**: 0.01 vs 0.03 → less severe V ratio at boundary (1e7 vs 3e7).

4. **SBP-SAT vs IP**: Tandem uses SBP-SAT (summation by parts — simultaneous approximation term), not IP (Interior Penalty). SBP-SAT does not have the p²/h penalty scaling that amplifies within-face discontinuities in IP.

---

## Decision: Revert Face-Centroid, Follow Tandem's Per-DOF Approach

### What was reverted:
- `config/bp5_params.hpp`: Removed `V_init_vec_nuc(is_nuc, V)`, `tau0_vec_nuc(x2, x3, is_nuc, tau)`, `Dc_nuc(is_nuc)` overloads. Restored original `V_init_vec()` and `tau0_vec()` implementations.
- `fault/fault_geometry.hpp`: Restored `ComputeBP5Params()` to per-DOF evaluation loop (matching Tandem).

### What was kept:
- `domain/domain_operator.hpp`: `GetNbfPerFace()` and `GetNumFaultFaces()` virtual methods with defaults — needed for multi-DOF plan Phase 2+.
- `domain/elasticity_operator.hpp`: `override` on `GetNumFaultFaces()` and `GetNbfPerFace()` — structural improvement.
- `fault/fault_geometry.hpp`: `nbf_per_face_` and `num_fault_faces_` member variables — structural, read from domain operator.

### Proper fix for p≥2 instability:

The p≥2 immediate earthquake at SCEC initialization is caused by the combination of:
1. delta_tau overstress (SCEC-specific, not present in Tandem)
2. Within-face V gradient amplified by IP penalty

The correct fix is NOT to smooth parameters across faces, but to use **`--psi-init tandem`** which:
- Removes delta_tau overstress (delta_tau_factor=0)
- Uses Tandem-style V_nuc=0.01 (not 0.03)
- Computes psi_init from actual stress for exact equilibrium
- This independently prevents the immediate earthquake without changing parameter evaluation

### Verification after revert:
- Build: ✅ (seas_bp5_full compiles)
- Tests: ✅ (46/46 pass)
- p=1 behavior: unchanged (per-DOF at nbf=1 is identical to face-centroid at nbf=1)

---

## Phase 2: Adopt Tandem Initialization Parameters

After reverting the face-centroid approach, we adopted Tandem's initialization parameters exactly. This addresses the COMMON mismatch between MFEM (both p=1 and p=2) and Tandem reference data.

### Changes Made

#### 1. `bp5_params.hpp` — Default parameter values

| Parameter | Old (SCEC) | New (Tandem) | Tandem source |
|-----------|-----------|-------------|---------------|
| `V_nuc` | 0.03 m/s | **0.01 m/s** | `bp5.lua:72` |
| `delta_tau_factor` | 1.0 | **0.0** | `bp5.lua` has no `delta_tau` function |

With these defaults, tau_pre is computed as:
```
tau_pre = sn * a * asinh(Vi/(2V0) * exp(psi_ss/a)) + eta * Vi
```
No delta_tau is added. This matches Tandem's `bp5.lua:94-103` exactly.

Old SCEC values can be recovered via command-line flags:
```
--V-nuc 0.03 --delta-tau-factor 1
```

#### 2. `rate_state_fault.hpp` — Default psi_init mode

| Setting | Old | New |
|---------|-----|-----|
| `scec_psi_init_` | `true` (SCEC fixed psi) | **`false`** (Tandem InitialStatePsi) |

With `scec_psi_init_ = false`, psi is computed from the actual stress at t=0 via `InitialStatePsi()`, which finds psi such that the system is in exact equilibrium at V = V_init.

Note: With delta_tau_factor = 0, both SCEC and Tandem psi modes give the same result:
- SCEC: psi = f0 + b*ln(V0/V_init) = f0 + b*ln(V0/Vp) = 0.8072
- Tandem: psi = a*ln(2*V0/Vi * sinh((|tau_pre|-eta*Vi)/(a*sn))) = 0.8072 (cancels to same value)

#### 3. `bp5_verification_full.cpp` — Default mode

| Setting | Old | New |
|---------|-----|-----|
| `psi_init_mode_str` | `"scec"` | **`"tandem"`** |

Override with `--psi-init-mode scec` for SCEC behavior.

### Mathematical Proof of Equilibrium

With the Tandem defaults, every DOF starts in exact stress equilibrium:

1. tau_pre = sn*a*asinh(Vi/(2V0)*exp(psi_ss/a)) + eta*Vi
2. At t=0 with zero slip, elastic traction ≈ 0
3. tau_total = tau_pre + 0 = tau_pre
4. psi_init = psi_ss = f0 + b*ln(V0/Vp) = 0.8072
5. Friction solve: tau_total = sn*a*asinh(V/(2V0)*exp(psi/a)) + eta*V → V = Vi

At nucleation DOFs (Vi = 0.01): V = 0.01 m/s exactly (no overstress).
At non-nucleation DOFs (Vi = 1e-9): V = 1e-9 m/s exactly.

Compare with old SCEC defaults:
- tau_pre included delta_tau = eta*V_nuc = 138.7 kPa (overstress)
- With psi = 0.8072 (not absorbing delta_tau), V = 0.049 m/s at nucleation (>> V_nuc = 0.03)

### Unit Tests Added

4 new test suites (23 new tests, total 119):

1. **TestTandemInitEquilibrium** (12 tests): Verifies V = V_init and psi = psi_ss at 6 spatial points (VW, nucleation, VS, transitions).

2. **TestTandemTauPreFormula** (4 tests): Verifies our tau_pre matches Tandem's bp5.lua formula at nucleation and non-nucleation DOFs.

3. **TestScecModeOverride** (4 tests): Verifies SCEC mode can be recovered with explicit overrides (V_nuc=0.03, delta_tau_factor=1.0) — confirms overstress and InitialStatePsi absorption.

4. **TestPsiInitNumerical** (2 tests): Verifies psi_init = 0.8072 matches Tandem's computed value.

### Verification

- Build: ✅ (seas_bp5_full, seas_test_bp5_params)
- Unit tests: ✅ (119/119 pass)
- Full test suite: ✅ (all tests pass)

---

## Phase 3: v46 Results Analysis — p=2 Divergence Investigation

**Date**: 2026-03-20
**Runs analyzed**: v46a (p=1, 1000m, IP, mumps-blr), v46b (p=2, 1000m, IP, mumps-blr)
**Status**: Both runs still executing; partial data analyzed (p1 at ~999yr / 4 events, p2 at ~77yr / 1 event). Full-run plots available.

---

### 3.1 Visual Summary from Plots

Both runs use the Tandem initialization defaults (V_nuc=0.01, delta_tau_factor=0, psi_init_mode=tandem). The initial earthquake bug from Phase 1 is eliminated — both start in equilibrium.

**Closeup (0–0.1 yr): Excellent p=2 match**

At all 10 stations, the closeup plots show p=2 tracking Tandem's p=4 reference closely, often overlapping with p=1:

| Station | Closeup quality (p=2 vs Tandem) | Notes |
|---------|-------------------------------|-------|
| strk+00dp+00 | ✅ Excellent | V_strike, tau, state all match |
| strk+00dp+10 | ✅ Excellent | Nearly identical to Tandem |
| strk+00dp+22 | ⚠️ Good with offset | V_dip shows oscillation at VW-VS boundary |
| strk+16dp+00 | ✅ Excellent | Both p1 and p2 close to Tandem |
| strk+16dp+10 | ✅ Excellent | Very close match |
| strk+36dp+00 | ⚠️ Moderate | tau_dip and V_dip show ~15% deviation from t=0 |
| strk-16dp+00 | ✅ Excellent | Good agreement |
| strk-16dp+10 | ✅ Excellent | Nearly identical |
| strk-24dp+10 | ✅ Excellent | Near nucleation center, good match |
| strk-36dp+00 | ⚠️ Moderate | V_dip elevated, tau_dip offset from t=0 |

**Full run (0–1800 yr): Progressive p=2 divergence**

Over multiple earthquake cycles, p=2 shows growing divergence from Tandem:

| Observable | p=1 vs Tandem | p=2 vs Tandem |
|------------|--------------|--------------|
| Event timing (1st) | ~250yr ≈ reference | ~250yr ≈ reference |
| Event timing (later) | Slight drift, stays close | **Progressive drift, events shift by ~10-20yr by event 5+** |
| Recurrence interval | ~250yr (stable) | **Shortening trend** visible |
| Strike slip | Good match all events | Diverges ~event 3-4 |
| Dip slip | Good match | **Larger divergence from event 2+** |
| Shear stress (dip) | Matches Tandem pattern | **Oscillatory, different amplitude** |
| State variable (psi) | Matches evolution | Good early, diverges late |

**Key observation**: The divergence is WORSE at off-center and deep stations:
- strk±36dp+00: Dip-slip component diverges significantly by ~750yr
- strk+00dp+22: VW-VS boundary shows dip-stress drift
- strk+16dp+00: Dip-slip offset grows progressively

---

### 3.2 Quantitative Time-Stepping Analysis

**Step count comparison (same simulation period):**

| Period | p=1 steps | p=2 steps | p2/p1 ratio |
|--------|-----------|-----------|-------------|
| 0–70 yr | 1,093 | 2,873 | **2.6×** |
| Full run to data end | 63,804 (999yr) | 16,592 (77yr) | N/A (different durations) |

**Time step statistics (interseismic, t=1–70yr):**

| Statistic | p=1 | p=2 | Ratio |
|-----------|-----|-----|-------|
| Average dt | 1.991e6 s (0.063 yr) | 7.574e5 s (0.024 yr) | **p2 is 2.6× smaller** |
| Max dt | 3.156e6 s (0.100 yr) | 1.968e6 s (0.062 yr) | **p2 capped 38% lower** |
| Min dt | 1.309e5 s | 3.229e4 s | p2 needs 4× smaller min |

**Implications**: p=2 is running 2.6× slower in wall-clock time per simulation year, and the smaller time steps are driven by the error controller finding larger local truncation errors.

---

### 3.3 Root Cause Analysis: Why p=2 Diverges More

Five contributing factors identified, ranked by likely impact:

#### Factor 1 (HIGH): 1/3 IP Penalty — More Damaging at p=2

**Finding**: MFEM's IP penalty is exactly 1/3 of Tandem's across all polynomial orders.

| Code | Geometric ratio in penalty | Physical meaning |
|------|---------------------------|-----------------|
| Tandem | `area / volume = A_phys / V_phys` | Physical face-to-volume ratio |
| MFEM | `nl_q / Weight() = 2A_phys / (6V_phys) = A_phys / (3V_phys)` | Reference element conventions |

The factor of 3 comes from MFEM's reference tetrahedron conventions:
- `CalcOrtho |nor| = 2 × A_phys` (reference triangle area = 1/2)
- `Weight() = det(J) = 6 × V_phys` (reference tet volume = 1/6)
- Ratio: `(2A)/(6V) = A/(3V)` instead of `A/V`

**Why it's worse at p=2**: The theoretical minimum penalty for coercivity scales as `c_N_1 = p(p+2)/3`:
- p=1: c_N_1 = 1.0, effective penalty = 1/3 of theoretical
- p=2: c_N_1 = 2.67, effective penalty = 2.67/3 = 0.89 of theoretical (CLOSER to threshold)
- p=4: c_N_1 = 8.0, effective penalty = 8/3 = 2.67 of theoretical

The penalty formula: `eta_F = (D+1) × c_N_1 × (A/(3V)) × c1²/c0`

With the 1/3 factor, the stability margin shrinks as p increases. At p=1, other stabilizing mechanisms (rate-state friction damping, radiation damping term eta*V) compensate. At p=2, the margin is thinner, and numerical errors in traction computation are less damped by the penalty.

**Status**: The v45 debug document explains why simply multiplying by 3 doesn't work (it disrupts the consistency-penalty balance since only penalty, not consistency/symmetry terms, is at 1/3 strength). The proper fix requires understanding how the DG bilinear form terms scale together.

*Source*: `elasticity_operator.hpp:1076-1079` (bilinear form), `3092-3100` (ComputeTraction)
*Tandem ref*: `Elasticity.cpp:278-291`, `InverseInequality.h:27-29`

#### Factor 2 (HIGH): L∞ Error Norm Over 6× More DOFs

**Finding**: The RK45 adaptive time stepper uses L∞ error norm: `err = max_i |err_i| / atol`.

At p=2, there are 6 DOFs per fault face (vs 1 at p=1), meaning the state vector is 6× longer. The L∞ norm picks the single WORST error across all DOFs. With 6× more DOFs:

1. More opportunities for a DOF at a parameter transition to have large local error
2. DOFs near the nucleation zone boundary (VW-VS transition) see sharp parameter gradients
3. A single "bad" DOF forces globally smaller time steps

**Quantitative evidence**: The max dt achieved by p=2 (0.062 yr) is 38% smaller than p=1 (0.100 yr = dt_max), suggesting p=2 frequently hits the error tolerance ceiling while p=1 coasts at dt_max.

The smaller time steps mean MORE time steps per earthquake cycle, and each step accumulates floating-point errors differently. Over 7+ events, these accumulate into visible timing drift.

**Possible fix**: Consider using L2 (RMS) error norm instead of L∞, which would average over DOFs rather than picking the worst. Tandem uses L∞ too (`-ts_adapt_wnormtype infinity`), but Tandem's p=2 default may work because of the full-strength penalty and SBP-SAT.

*Source*: `time_stepper.hpp:338-386` (error norm), `bp5_verification_full.cpp:971-974` (tolerances)
*Tandem ref*: `examples/options/rk45.cfg` (atol=1e-7, rtol=1e-50, wnormtype=infinity)

#### Factor 3 (MEDIUM): Stiffer ODE System

The coupled elasticity + rate-state ODE system becomes stiffer at p=2:

1. **Penalty stiffness**: The effective spring constant from the IP penalty scales as `p²/h`. At p=2 this is 2.67× stiffer, creating faster-decaying modes that the RK45 stages capture as larger error estimates.

2. **More spatial resolution**: p=2 resolves within-element stress gradients that p=1 cannot see. These sharper features require smaller time steps to integrate accurately.

3. **Friction nonlinearity**: The rate-state friction law has exponential dependence on stress: `V = 2V0*sinh(tau/(a*sn))*exp(-psi/a)`. At p=2, the traction at individual DOFs can be sharper (not face-averaged), leading to more extreme V values that dominate the error estimate.

#### Factor 4 (MEDIUM): Node Distribution — GaussLobatto vs WarpAndBlend

**Finding**: SEAS-MFEM uses GaussLobatto nodes on triangles; Tandem uses WarpAndBlend nodes.

At p=2, both give 6 nodes but at different positions:
- **GaussLobatto**: 3 vertices + 3 edge midpoints
- **WarpAndBlend**: Optimized for interpolation (lower Lebesgue constant)

The Lebesgue constant affects the stability of the L2 projection (smaller is better). WarpAndBlend nodes are specifically designed for simplex interpolation and may provide better conditioning for the mass matrix inverse used in `GalerkinProject`.

For BP5's flat fault faces, both should give spectrally convergent results, but the conditioning difference could explain why p=2 accumulates slightly different errors cycle-to-cycle.

*Source*: `face_quadrature.hpp:58` — `H1_TriangleElement(face_order, BasisType::GaussLobatto)`
*Tandem ref*: `RateAndStateBase.cpp:10-13` — `NodalRefElement<2>(PolynomialDegree, WarpAndBlendFactory<2>())`

#### Factor 5 (LOW): dt_max Difference

| | SEAS-MFEM | Tandem (QD mode) |
|--|-----------|-------------------|
| dt_max | 0.1 yr (hardcoded) | PETSC_MAX_REAL (unlimited) |

In Tandem's QD mode, `cfl_time_step()` returns `std::nullopt` → no dt_max is set. Tandem's RK45 can take arbitrarily large steps if the error is small enough.

Our p=1 hits dt_max (0.1yr) frequently during interseismic periods. Our p=2 never reaches it (max was 0.062yr). So dt_max is NOT the cause of p=2 divergence, but it means p=1 is artificially capped — possibly hiding the fact that it could take even larger steps (and potentially diverge more).

*Source*: `bp5_verification_full.cpp:974` — `rk45.SetMaxDt(0.1 * year_s)`
*Tandem ref*: `SEAS.cpp:198-201` — QD mode returns nullopt, no dt_max

---

### 3.4 Algorithmic Comparison Table

| Aspect | SEAS-MFEM | Tandem | Match? |
|--------|-----------|--------|--------|
| RK method | Dormand-Prince 5(4), 7-stage FSAL | Same (PETSc `5dp`) | ✅ |
| Error norm | L∞ weighted | L∞ weighted | ✅ |
| atol | 1e-7 | 1e-7 | ✅ |
| rtol | 1e-50 | 1e-50 | ✅ |
| safety factor | 0.9 | 0.9 | ✅ |
| reject_safety | 0.5 | 0.5 | ✅ |
| growth_max | 10.0 | 10.0 | ✅ |
| shrink_min | 0.1 | 0.1 | ✅ |
| dt_max (QD) | **0.1 yr** | **unlimited** | ❌ |
| CFL constraint | None | None | ✅ |
| DG method | IP (SIPG) | **SBP-SAT** | ❌ fundamental |
| Penalty c_N formula | p(p+D-1)/D | p(p+D-1)/D | ✅ |
| Penalty geometric factor | **A/(3V)** | **A/V** | ❌ factor of 3 |
| Interior penalty averaging | (p0+p1)/4 | (p(0)+p(1))/4 | ✅ |
| Traction projection | L2 (GalerkinProject) | L2 (minv * E^T * W * nl) | ✅ (flat faces) |
| Mass matrix | Consistent (not lumped) | Consistent (not lumped) | ✅ |
| Face node distribution | **GaussLobatto** | **WarpAndBlend** | ❌ |
| Fault basis rotation | Per-face (constant) | Per-quad-point | ✅ (flat fault) |
| Jacobian in projection | Reference (no nl_q) | Physical (with nl_q) | ✅ (flat faces) |
| Polynomial order | Runtime (`--order`) | Compile-time template | N/A |

---

### 3.5 Specific Numerical Evidence of Early Divergence

**At strk-36dp+00 (x2=-36km, x3=0km) at t≈0.04 yr:**

| Quantity | p=1 | p=2 | Deviation |
|----------|-----|-----|-----------|
| V_dip (log10) | -7.784 | -7.709 | p2 is 19% higher |
| tau_dip (MPa) | 1.279 | 1.235 | p2 is 3.4% lower |
| V_strike (log10) | -6.679 | -6.645 | p2 is 8.6% higher |

This is at t=0.04yr — well before any earthquake. The elastic solution quality already differs, with p=2 showing systematically higher slip rates and slightly different shear stress. This early-time deviation is NOT from time stepping but from the **spatial discretization** (traction computation, penalty, DG solution quality).

**Initial conditions are identical** at all stations (V_strike=1e-9, tau=13.273 MPa, psi=8.146128). The divergence emerges from the very first elastic solve.

---

### 3.6 Diagnosis: Most Likely Cause of Progressive Divergence

The evidence points to a combination of Factor 1 (1/3 penalty) and Factor 2 (L∞ error norm) as the primary drivers:

1. **The 1/3 penalty under-stabilization** produces a slightly different DG elasticity solution at p=2 vs Tandem's SBP-SAT, especially in the dip-direction stress (where the problem is geometrically asymmetric). This manifests as the systematic ~3-20% deviations seen at off-center stations from t=0.

2. **The smaller time steps** forced by L∞ over 6× more DOFs create a different error accumulation pattern. Over many earthquake cycles, the cumulative effect shifts event timing progressively.

3. The divergence is NOT from initialization (initial conditions match exactly), NOT from time stepping parameters (identical tolerances), and NOT from a code bug — it's from the **fundamental difference between IP with 1/3 penalty and SBP-SAT**.

---

### 3.7 Recommended Investigation Priorities

#### Priority 1: Understand the penalty-consistency balance
The v45 debug doc explains that simply multiplying penalty by 3 disrupts the consistency-penalty balance. Need to investigate whether MFEM's DG bilinear form (consistency + symmetry + penalty) is correctly balanced in reference-element terms, or whether there's a systematic scaling issue affecting ALL three terms.

**Key question**: Does MFEM's `DGElasticityIntegrator` produce correct physical integrals for all three terms, or is only the penalty affected by the 1/3 factor?

#### Priority 2: Test with dt_max = unlimited
Verify that removing the dt_max=0.1yr cap doesn't affect p=1 results (it should already be hitting dt_max regularly). If p=1 quality improves, it suggests the cap was hiding convergence issues.

#### Priority 3: Test with L2 (RMS) error norm
Switch the RK45 error norm from L∞ to RMS. This would average errors over all DOFs rather than picking the worst, potentially allowing p=2 to take larger time steps while maintaining accuracy on average.

#### Priority 4: WarpAndBlend nodes
Test switching from GaussLobatto to WarpAndBlend nodes for the fault face quadrature. This matches Tandem exactly and may improve the L2 projection conditioning.

#### Priority 5: Full penalty analysis at p=2
Compute the actual penalty values at specific faces and compare with Tandem's values. Verify whether the 1/3 factor brings the penalty below the coercivity threshold at p=2.

---

## Key Tandem Source Files Referenced

| File | What it shows |
|------|---------------|
| `tandem/examples/tandem/3d/bp5.lua` | All BP5 parameters defined as Lua functions of (x,y,z). V_nuc=0.01, no delta_tau, eps-based nucleation boundary. |
| `tandem/examples/tandem/3d/bp5.toml` | Default scenario = `bp5_outside` (eps=+1e-3). |
| `tandem/app/localoperator/RateAndStateBase.cpp:24-35` | Physical coordinates computed at each DOF position per face. |
| `tandem/app/localoperator/RateAndState.h:36-43` | `set_params()` iterates over ALL DOFs, calls param function with per-DOF coordinates. |
| `tandem/app/localoperator/RateAndState.h:122-146` | `init()`: computes psi_init per-DOF from actual traction + tau_pre. |
| `tandem/app/tandem/FrictionConfig.h:97-109` | `param_fun()` wraps Lua functions, evaluated at physical coordinate x. |
| `tandem/app/localoperator/DieterichRuinaAgeing.h:39-63` | `set_params()` stores per-DOF params; `psi_init()` computes equilibrium psi from total stress. |
