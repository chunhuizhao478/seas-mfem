# BP5 Debug v22: Deep VS Zone Lockup — Root Cause Analysis

**Date**: 2026-03-14
**Status**: Diagnosis complete, diagnostic tests proposed
**Previous**: v21 (FarField BC fix H25)
**Results**: `results_1000m_bp5_v21_med_mblr`
**Plots**: `plots_results_v21_med_mblr`

---

## 1. Problem Statement

MFEM produces **one earthquake at t≈0** and then no subsequent events over 300 years. Tandem produces 8-9 regular earthquake cycles with ~200 year recurrence over 1800 years.

After the single earthquake, the fault remains locked with negligible stress buildup:
- Stress loading rate at z=0 km: 0.0015 MPa/yr
- Stress loading rate at z=10 km: 0.015 MPa/yr
- Tandem's loading rate: ~0.07 MPa/yr

At these rates, the next event would take thousands of years.

---

## 2. Key Finding: Deep VS Zone Locks Up

The velocity-strengthening zone at z=22 km (a=0.04, b=0.03, a-b=0.01) stops creeping after the premature earthquake.

**Direct evidence from station strk+00dp+22 (x2=0, x3=22 km):**

| Time (yr) | V/Vp | θ/θ_ss | τ_strike (MPa) | Status |
|-----------|------|--------|-----------------|--------|
| 0.00 | 1.000 | 1.0 | 13.273 | Initial equilibrium |
| 0.01 | 810× | 0.001 | 14.748 | Post-seismic transient |
| 5 | 1.39 | 0.4 | 12.976 | Decelerating through Vp |
| 50 | **0.138** | **4.5** | 12.423 | **Locked** |
| 100 | **0.087** | **8.2** | 12.417 | **Locked** |
| 200 | **0.078** | **12.1** | 12.588 | **Locked** |
| 295 | **0.095** | **11.5** | 12.751 | **Locked** (barely recovering) |

**Expected behavior**: The VS zone should maintain V ≈ Vp (velocity-strengthening friction drives stable creep at plate rate).

**Actual behavior**: V drops to 8-14% of Vp. The state variable θ grows to 12× its steady-state value, making friction much stronger than in steady state.

### 2.1 Why the VS Zone Locks Up

The aging law creates a positive feedback loop when stress drops below steady state:

1. V drops below Vp
2. dθ/dt = 1 - V·θ/Dc > 0 → θ grows
3. Growing θ increases friction strength
4. Higher friction → lower V for the same stress
5. Lower V → θ grows even faster
6. The fault "freezes" at V ≈ 0.08×Vp

### 2.2 Why It Can't Recover

The velocity-strengthening is **weak**: (a-b)·σ_n·ln(10) = 0.01 × 25 MPa × 2.303 = **0.576 MPa** for a 10× change in V.

This means:
- Steady-state τ at V = Vp: 13.273 MPa
- Steady-state τ at V = 0.1×Vp: 12.697 MPa
- Steady-state τ at V = 0.01×Vp: 12.122 MPa
- Only 0.576 MPa separates V=Vp from V=0.1×Vp

The actual stress at t=50 yr is 12.42 MPa — 0.85 MPa below steady state at Vp. This is enough to keep V well below Vp.

Furthermore, the fault is NOT in steady state. At t=50 yr, θ = 4.5×θ_ss, so the friction strength at V=Vp with the inflated θ is **14.40 MPa** — a **2.0 MPa deficit** from the actual stress of 12.42 MPa. The boundary loading at ~0.005-0.015 MPa/yr would take 130-400 years to close this gap, but θ keeps growing during that time, making the gap worse.

---

## 3. The Premature Earthquake

### 3.1 What Happens

The v20 psi fix (H24) sets ψ(0) = f₀ + b·ln(V₀/V_init) everywhere (SCEC Eq. 18). Combined with the δτ overstress in the nucleation zone (SCEC Eq. 23), this produces V_max = 0.049 m/s at t=0, triggering an immediate earthquake.

The earthquake produces:
- ~5.2 m of coseismic slip at the hypocenter (z=0-10 km)
- ~1.2 m of coseismic slip at z=22 km (deep VS zone)
- Stress increase at z=22 km: 13.27 → 14.75 MPa

At t=0, the boundary displacement is **zero** (u_D = sgn(x)·Vp·0/2 = 0).

### 3.2 The Cascade

1. Post-seismic transient: V at z=22 km spikes to 810×Vp
2. Fast creep accumulates ~1.4 m additional slip in weeks
3. Total deep slip (~2.6 m) **far exceeds boundary displacement (~0)**
4. Stress drops from 14.75 → 12.4 MPa (below initial 13.27 MPa)
5. V drops below Vp at ~5 years
6. θ starts growing above θ_ss (positive feedback)
7. Deep VS zone locks up by ~50 years
8. Without deep creep, only weak boundary loading remains

### 3.3 Why Tandem Doesn't Have This Problem

In Tandem, `InitialStatePsi` absorbs δτ into ψ, so the system starts in near-equilibrium at V = V_init. The first event nucleates naturally at ~150 years after boundary + deep creep loading builds sufficient stress. At that point:
- Boundary has accumulated ~2.4 m of displacement per side
- The coseismic slip doesn't far exceed the accumulated boundary motion
- The post-earthquake stress deficit at depth is smaller
- The VS zone can recover to near plate-rate creep

---

## 4. Structural Differences: Tandem vs MFEM

### 4.1 Fault Depth Extent

| Aspect | Tandem | MFEM |
|--------|--------|------|
| Fault depth extent | **0-40 km** (stops at Wf) | **0-40 km** (stops at Wf) |
| Below fault | Domain is bonded (continuous) | Domain is bonded (continuous) |
| Deep VS zone (18-40 km) | Rate-state, a=0.04 | Rate-state, a=0.04 (same) |

Both codes limit the fault to 40 km depth. In MFEM's bp5.geo, the fault surface is
`Rectangle(fault) = {0, -l_f/2, 0, W_f, l_f}` with `W_f = 40` (km). The
`IsFaultFace3D` code also enforces `center(2) <= Wf_ + tol` (Wf_ = 40e3 m).
Below 40 km, the domain is simply connected in both codes.

**Note:** The `ComputeRHS` code in `rate_state_fault.hpp` has a branch
`if (depths(i) > Wf_bp5_ + 1.0)` that prescribes Vp creep. This is **dead code** —
no fault nodes exist below 40 km, so the condition never triggers.

### 4.2 Nucleation Parameters

| Parameter | Tandem | MFEM |
|-----------|--------|------|
| V_nuc (nucleation zone) | 0.01 m/s | 0.03 m/s |
| V_init (elsewhere) | 1e-9 m/s | 1e-9 m/s |
| V_zero (dip component) | 1e-20 m/s | 1e-20 m/s |
| δτ (QD) | η × 0.01 ≈ 0.046 MPa | η × 0.03 ≈ 0.139 MPa |
| ψ initialization | Absorbs δτ into ψ | SCEC-correct: ψ = f₀+b·ln(V₀/V_init) |
| First earthquake timing | ~150 years (natural) | t=0 (immediate from δτ) |

Tandem uses V_nuc = 0.01 m/s (from bp5.lua line 73). MFEM uses V_nuc = 0.03 m/s (bp5_params.hpp line 112). Both are within SCEC spec range but produce different δτ magnitudes.

### 4.3 Boundary Conditions

Both use the same approach:
- Dirichlet on far-field vertical faces: u = (0, sgn(x)·Vp·t/2, 0)
- Natural (zero traction) on top (z=0) and bottom (z=Lz)

MFEM attrs 1-4 = Dirichlet, attrs 5-6 = Natural (FarField BCMode, H25 fix).

### 4.4 Deep Fault Handling

Tandem does **NOT** have prescribed plate-rate creep. The entire fault (0-40 km) uses rate-state friction with depth-dependent `a` parameter. The VS zone (18-40 km) naturally creeps at ~Vp because a > b makes it velocity-strengthening.

MFEM has an explicit cutoff at `Wf_bp5_ + 1.0 = 41001 m`:
```cpp
if (depths(i) > Wf_bp5_ + 1.0)
{
    rate(i * StatePerNode + 0) = 0.0;       // dip rate = 0
    rate(i * StatePerNode + 1) = Vp_bp5_;   // strike rate = Vp
    rate(i * StatePerNode + PsiIndex) = 0.0; // no state evolution
    continue;
}
```

---

## 5. Boundary Loading Analysis

### 5.1 Analytical Estimate

For a 3D domain [-Lx,Lx] × [-Ly,Ly] × [0,Lz] with locked fault at x=0:

Simple 2D antiplane: dτ/dt = μ·Vp/(2·Lx) = 32.04e9 × 1e-9 / (2 × 100e3) = **0.005 MPa/yr**

3D correction (fault finite depth): lower by factor ~a/(π·L) ≈ 0.06 → **~0.0003-0.001 MPa/yr**

### 5.2 Observed Loading Rates (from MFEM data, post-earthquake)

| Station | Depth | Loading Rate | Notes |
|---------|-------|-------------|-------|
| strk+00dp+00 | 0 km | 0.0015 MPa/yr | Surface, complicated by dip-stress rotation |
| strk+00dp+10 | 10 km | 0.015 MPa/yr | VW core, fault locked |
| strk+00dp+22 | 22 km | ~0 (decreasing) | Deep VS, fault locked |

The observed rate at z=10 km (0.015 MPa/yr) is in the right ballpark for boundary-only loading plus some contribution from the prescribed deep zone (41-100 km).

### 5.3 Comparison with Tandem

Tandem's loading rate is ~0.07 MPa/yr. The difference from MFEM's boundary-only rate (~0.005-0.015 MPa/yr) comes from the **deep VS zone creep**. In Tandem, the VS zone (18-40 km) creeps at ~Vp, creating differential motion with the locked VW zone. This deep creep is the dominant loading mechanism, not the far-field boundary.

In MFEM after the premature earthquake, the VS zone locks up, removing this dominant loading contribution. Only the weak far-field boundary loading remains.

### 5.4 Surface Station Stress Rotation Anomaly

At strk+00dp+00 (z=0 km, surface), after the earthquake:
- τ_strike: 13.27 → -9.12 MPa (went negative!)
- τ_dip: 0 → 13.30 MPa (huge spurious dip stress!)

This stress "rotation" from strike to dip is not observed at deeper stations (z=10, z=22 km have τ_dip ≈ 0 throughout). It may be a free-surface effect or a traction computation issue specific to the z=0 surface station. This deserves investigation but is likely not the primary cause of the non-cycling behavior.

---

## 6. Friction Parameters Reference

### 6.1 SCEC BP5 Parameters (from bp5_params.hpp)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Shear modulus | μ | 32.04 GPa |
| Shear wave speed | c_s | 3464 m/s |
| Normal stress | σ_n | 25 MPa |
| Reference velocity | V₀ | 1e-6 m/s |
| Reference friction | f₀ | 0.6 |
| VW core a | a₀ | 0.004 |
| VS zone a | a_max | 0.04 |
| State evolution b | b | 0.03 |
| Characteristic slip | L₀ (Dc) | 0.14 m |
| Plate rate | V_p | 1e-9 m/s |
| Radiation damping | η = μ/(2c_s) | 4.626e6 Pa·s/m |
| Fault width | Wf | 40 km |
| Seismogenic depth hs | h_s | 2 km |
| Transition width | h_t | 2 km |
| VW core height | H | 12 km |

### 6.2 Depth Zones

| Depth Range | Zone | a | a-b | Behavior |
|-------------|------|---|-----|----------|
| 0-2 km | Shallow VS | 0.04 | +0.01 | Velocity-strengthening |
| 2-4 km | Shallow transition | 0.004-0.04 | varies | Transition |
| 4-16 km | VW core | 0.004 | -0.026 | Velocity-weakening |
| 16-18 km | Deep transition | 0.004-0.04 | varies | Transition |
| 18-40 km | Deep VS | 0.04 | +0.01 | Velocity-strengthening |
| >40 km | No fault | N/A | N/A | Domain bonded (no slip) |

### 6.3 Steady-State Stress at z=22 km (a=0.04)

| V/Vp | τ_ss (MPa) |
|------|-----------|
| 1.0 | 13.273 |
| 0.5 | 13.100 |
| 0.1 | 12.697 |
| 0.01 | 12.122 |

Only 0.576 MPa separates V=Vp from V=0.1×Vp in steady state.

---

## 7. Diagnostic Tests

### Test 1: Run with delta_tau = 0 (no nucleation perturbation)

**Purpose**: Determine if the deep VS zone maintains V ≈ Vp without earthquake disruption.

**Method**: Set `delta_tau_factor = 0` in BP5Params (or pass a CLI flag). Run for ~500 years.

**Expected outcomes**:
- If VS zone maintains V ≈ Vp → the premature earthquake is the sole cause of lockup
- If VS zone still drifts below Vp → there is a separate boundary loading problem

**What to check**:
- V at z=22 km over time (should stay near Vp)
- Stress at z=10 km over time (should build at ~0.07 MPa/yr from combined boundary + deep creep)
- Compare loading rate with Tandem's pre-first-event loading rate

### Test 2: Boundary loading traction verification

**Purpose**: Verify the DG Dirichlet boundary loading produces correct traction at the fault.

**Method**: Solve elasticity at t=1yr with zero slip, compute traction at fault center.

**Analytical prediction**: For locked fault at x=0, boundary u_y = sgn(x)·Vp·1yr/2:
- τ_strike ≈ μ·Vp·1yr/(2·Lx) = 32.04e9 × 1e-9 × 3.15e7 / (2×100e3) ≈ **5050 Pa ≈ 0.005 MPa**

**Implementation**: The `--diag-vtk` code path already does this solve (bp5_verification_full.cpp line 688). Add traction computation:
```cpp
domain.ComputeTraction(u_diag, zero_slip, diag_traction);
// Print traction at fault center and compare with analytical estimate
```

**What to check**:
- Is the traction magnitude ~0.005 MPa at z=10 km?
- Is it purely along-strike (τ_dip ≈ 0)?
- Does it match the analytical estimate within a factor of 2?

### ~~Test 3: Match Tandem fault depth (40 km only)~~ — NOT NEEDED

The fault already extends only to 40 km in both the mesh (`bp5.geo: W_f = 40`)
and the code (`IsFaultFace3D: z <= Wf_`). Below 40 km the domain is bonded,
matching Tandem. The prescribed Vp creep in `ComputeRHS` is dead code.

### Test 3 (renumbered): Run with absorbed δτ (Tandem-style initialization)

**Purpose**: Verify that matching Tandem's initialization produces matching behavior.

**Method**: In the BP5 vector path of `Init()`, use `InitialStatePsi(tau_abs, V_abs_init, sigma_n, eta, a)` instead of `bp5_params_.psi_init()`. This absorbs δτ into ψ.

**What to check**:
- First event timing (~150 years, matching Tandem?)
- Subsequent recurrence interval (~200 years?)
- Deep VS zone V ≈ Vp throughout interseismic?
- Loading rate at z=10 km matches Tandem?

---

## 8. Possible Root Causes (Ranked by Likelihood)

### 8.1 Primary: Premature earthquake disrupts deep VS equilibrium (HIGH confidence)

The v20 psi fix (H24) creates SCEC-correct initial conditions where δτ is genuine overstress. Combined with zero initial boundary displacement, this triggers an immediate earthquake. The earthquake creates a large stress deficit at depth that the boundary loading cannot recover, causing the deep VS zone to lock up through the aging law positive feedback.

**Evidence**: Direct observation of deep VS lockup in station data. The mechanism is physically consistent with the weak velocity-strengthening (a-b=0.01) and slow boundary loading rate.

### 8.2 Secondary: Boundary loading may be too weak (MEDIUM confidence)

The observed boundary-only loading rate (~0.005-0.015 MPa/yr) is consistent with analytical estimates for the 3D geometry. However, I have not verified the DG BR2 Dirichlet loading implementation with a precise numerical test. A factor-of-2 error in the boundary loading would be masked by the premature earthquake issue but could affect subsequent cycling even with correct initialization.

**Evidence**: Analytical estimate matches observed rate within ~3×, but a precise verification test is needed (Test 2).

### 8.3 ~~Tertiary: Fault depth mismatch with Tandem~~ — RULED OUT

~~MFEM's fault extends to 100 km with prescribed Vp creep below 41 km.~~

**Correction**: Both MFEM and Tandem limit the fault to 40 km depth. The mesh
(`bp5.geo: W_f = 40`) and the fault detection code (`IsFaultFace3D: z <= Wf_`)
both enforce this. Below 40 km the domain is bonded in both codes. The prescribed
Vp creep code in `ComputeRHS` is dead code — no fault nodes exist at z > 40 km.

### 8.4 Possible: V_nuc mismatch (LOW confidence)

MFEM uses V_nuc = 0.03 m/s, Tandem uses V_nuc = 0.01 m/s. This affects the δτ magnitude (0.139 vs 0.046 MPa) and the initial slip rate in the nucleation zone. However, this is a small quantitative difference and unlikely to be the primary cause.

---

## 9. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H12 | Friction, output, penalty, BC, nucleation fixes | Done |
| H13 | MUMPS print level 0→1 | Done (v17) |
| H14 | Post-solve residual check | Done (v17), fixed in v18 (H18) |
| H15 | AMG SetSystemsOptions(3) → SetElasticityOptions | Done (v17), upgraded (v18) |
| H16 | Runtime --solver flag (default: cg) | Done (v17) |
| H17 | Traction monitoring at fault stations | Done (v17) |
| H18 | Global residual norms via MPI_Allreduce | Done (v18) |
| H19 | CG RelTol 1e-12 → 1e-10 | Done (v18) |
| H20 | GMRES+BlockILU solver option | Done (v18) |
| H21 | MUMPS BLR, SuperLU, STRUMPACK solver options | Done (v19) |
| H22 | GMRES+BlockILU solver option, faster than CG+AMG | Done (v19) |
| H23 | MUMPS BLR set as default solver | Done (v19) |
| H24 | Fix InitialStatePsi absorbing delta_tau → use psi_init() | Done (v20) |
| H25 | FarField BC mode (attrs 1-4 Dirichlet, 5-6 Natural) | Done (v21) |

---

## 10. Files Referenced

| File | Relevance |
|------|-----------|
| `config/bp5_params.hpp` | BP5 parameters, a(x2,x3), tau0_vec, V_init_vec |
| `fault/rate_state_fault.hpp` | ComputeRHS (deep creep cutoff at Wf+1), Init (psi initialization) |
| `fault/fault_geometry.hpp` | ComputeBP5Params, tau_pre computation |
| `domain/elasticity_operator.hpp` | AssembleDirichletLoading (BR2 path), ComputeTraction, BCMode |
| `solver/seas_operator.hpp` | SEASQuasiDynamicOperator::Mult (time propagation to Solve) |
| `solver/time_stepper.hpp` | DormandPrinceRK45 (ODE integrator) |
| `integrator/dg_elasticity_br2_integrator.hpp` | DGElasticityBR2BoundaryIntegrator |
| `tests/verification/bp5_verification_full.cpp` | Main driver, CLI options, diagnostic VTK |
| `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/bp5.lua` | Tandem BP5 config |
| `/Users/chunhuizhao/projects/tandem/examples/tandem/3d/bp5.geo` | Tandem BP5 mesh (fault to 40 km) |
| `/Users/chunhuizhao/projects/tandem/app/localoperator/Elasticity.cpp` | Tandem DG boundary handling |
