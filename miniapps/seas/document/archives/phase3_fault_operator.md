# Phase 3: Fault Operator (Serial)

## Overview

This phase implements the fault operator that manages rate-and-state friction computation, fault geometry, and state variable evolution on the fault surface. **Following Tandem's implementation** as the reference.

## Reference Implementation

**Tandem source files** (in `/Users/chunhuizhao/projects/tandem/`):
- `app/localoperator/DieterichRuinaAgeing.h` - Ageing law friction implementation
- `app/localoperator/RateAndState.h` - Rate-and-state local operator
- `app/form/FrictionOperator.h` - Friction operator wrapper
- `src/util/Zero.h`, `Zero.cpp` - Brent's method root finder
- `examples/tandem/2d/bp1.lua` - BP1/BP2 configuration

**Theory document**: See `dg_antiplane_theory.md`:
- **Section 6**: Traction computation: τ = μ · {{∂u/∂x}}
- **Section 4.4.2 / 5.3.2**: Fault slip contribution to RHS

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 3.1 Fault geometry | `fault_geometry.hpp` | `test_fault_operator.cpp` |
| 3.2 Friction law | `friction_law.hpp`, `dieterich_ruina.hpp` | `test_friction.cpp` |
| 3.3 State evolution | `state_evolution.hpp` | `test_friction.cpp` |
| 3.4 Rate-state fault operator | `rate_state_fault.hpp` | `test_fault_operator.cpp` |
| 3.5 BP2 parameters | `bp2_params.hpp` | `test_fault_operator.cpp` |

## Verification Checkpoint

Initial state matches benchmark data (tau0 ~ 26.546 MPa, theta(0) computed from equilibrium).

---

## BP2 Fault Configuration

### Fault Parameters (from SCEC BP2 specification)

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Normal stress | σ_n | 50 | MPa |
| Rate-state a (VW zone) | a0 | 0.010 | - |
| Rate-state a (VS zone) | amax | 0.025 | - |
| Rate-state b | b | 0.015 | - |
| Characteristic length | Dc | 0.004 | m |
| Plate rate | Vp | 10^-9 | m/s |
| Initial slip rate | Vinit | 10^-9 | m/s |
| Reference slip rate | V0 | 10^-6 | m/s |
| Reference friction | f0 | 0.6 | - |
| VW zone depth | H | 15 | km |
| VW-VS transition width | h | 3 | km |
| Rate-state fault depth | Wf | 40 | km |

### Radiation Damping

In the implementation, radiation damping is computed as:
```
η = μ / (2·cs) = 32.04e9 / (2 · 3464) ≈ 4.63 MPa·s/m
```

This is a method on `BP2Params`:
```cpp
real_t eta() const { return mu() / (2.0 * cs); }
```

---

## Friction Law: Dieterich-Ruina Regularized (Following Tandem)

**Implementation Note:** The actual implementation class is `DieterichRuinaFriction` in `friction/dieterich_ruina.hpp`, using the θ (theta) state variable form.

### Friction Coefficient

From Tandem `DieterichRuinaAgeing.h` line 73-79:

**Regularized friction coefficient** (Lapusta et al. form):
```
f(V, ψ) = a · asinh[(V / 2V₀) · exp(ψ/a)]
```

where:
- ψ = state variable (natural log form)
- V = slip rate
- V₀ = reference slip rate
- a = direct effect parameter (depth-dependent)

**Shear strength**:
```
τ_f = σ_n · f(V, ψ) + η · V
```

### State Evolution (Ageing Law)

From Tandem `DieterichRuinaAgeing.h` lines 122-125:

```cpp
double state_rhs(std::size_t index, double V, double psi) const {
    double myL = p_[index].get<L>();
    return cp_.b * cp_.V0 / myL * (exp((cp_.f0 - psi) / cp_.b) - V / cp_.V0);
}
```

**Mathematical form**:
```
dψ/dt = (b·V₀/L) · [exp((f₀ - ψ)/b) - V/V₀]
```

This is the **Ageing (Dieterich-Ruina) evolution law** for state variable ψ.

**Alternative form** (in terms of θ = exp(ψ/b)):
```
dθ/dt = 1 - V·θ/L
```

### Depth-Dependent Parameter a(z)

Implemented as `BP2Params::a_of_z(z)` where z is the depth coordinate (z=0 at surface, z<0 at depth):

```cpp
real_t a_of_z(real_t z) const
{
   real_t depth = std::abs(z);  // Convert to positive depth
   if (depth < H)
      return a0;
   else if (depth < H + h)
      return a0 + (amax - a0) * (depth - H) / h;
   else
      return amax;
}
```

Where H = 15 km, h = 3 km, a0 = 0.010, amax = 0.025.

Also provides `IsVelocityWeakening(z)` which returns `a_of_z(z) < b`.

---

## Slip Rate Solving (Brent's Method)

### Stress Balance at Fault

The quasi-dynamic stress balance at each fault point is:

```
┌───────────────────────────────────────────────────────────────────┐
│  τ_total = τ_0 + τ_qs = σ_n · f(V, ψ) + η · V                    │
└───────────────────────────────────────────────────────────────────┘
```

where:
- τ_total = total shear stress
- τ_0 = pre-stress (uniform, computed at t=0 for equilibrium)
- τ_qs = quasi-static traction from domain solve (from `dg_antiplane_theory.md` Section 6)
- σ_n = normal stress (50 MPa for BP2)
- f(V, ψ) = friction coefficient
- η = radiation damping coefficient = μ/(2cs)
- V = slip rate (unknown to solve for)

### Total Stress Decomposition

```
τ_total = τ_0 + τ_qs

where τ_qs = μ · {{∂u/∂x}} is computed by DomainOperator::ComputeTraction()
```

### Force Balance

The slip rate V is found by solving the force balance equation:

```
τ = σ_n · f(V, ψ) + η · V
```

where τ = |τ_total| is the absolute shear traction.

**Residual function** (from Tandem `DieterichRuinaAgeing.h` lines 81-120):

```
R(V) = τ - σ_n · f(V, ψ) - η · V = 0
```

### Slip Rate Solving (Newton-Raphson)

The MFEM SEAS implementation uses **Newton-Raphson** with safeguarded bounds (implemented in `DieterichRuinaFriction::SolveSlipRate`):

- Initial guess: Previous V or V₀
- Lower bound: `V_min_ = 1e-30` (avoid division by zero)
- Upper bound: `V_max = tau / eta` (physical limit when friction is zero)
- Convergence: relative tolerance 1e-12, max 100 iterations
- Safeguard: if Newton step goes out of bounds, halve the step (bisection-style fallback)

The derivative `dR/dV` is available analytically:
```
R(V) = τ - σ_n · f(V, θ) - η · V = 0
dR/dV = -σ_n · df/dV - η
```

where `df/dV = a / √(V² + (2V₀·exp(-ψ/a))²)` and `ψ = f₀ + b·ln(V₀·θ/Dc)`.

Special case: when `σ_n ≤ 0`, the friction law degenerates to viscous sliding: `V = τ / η`.

---

## Initial State Computation

### Pre-Stress τ₀

Implemented as `BP2Params::tau0()`. The pre-stress is computed at steady state with `V = V_init` in the VS zone (where `a = amax`):

```cpp
real_t tau0() const
{
   // At steady state: theta_ss = Dc / V_init
   // f_ss = amax * asinh[(V_init / 2V0) * exp((f0 + b*ln(V0/V_init)) / amax)]
   real_t log_term = f0 + b * std::log(V0 / V_init);
   real_t exp_term = std::exp(log_term / amax);
   real_t arg = (V_init / (2.0 * V0)) * exp_term;
   real_t f_ss = amax * std::asinh(arg);
   return sigma_n * f_ss + eta() * V_init;
}
```

For BP2 parameters: **τ₀ ≈ 26.546 MPa**

### Initial θ (State Variable)

Implemented in `DieterichRuinaFriction::InitialState()`. Given τ₀ and V_init, the initial θ is computed from the stress balance equation:

```
τ₀ = σ_n · f(V_init, θ₀) + η · V_init
```

The implementation:
1. Computes `τ_eff = τ₀ - η · V_init`
2. Computes the equivalent ψ: `ψ = a · ln[2V₀/V_init · sinh(τ_eff / (a·σ_n))]`
3. Converts to θ: `θ = (Dc/V₀) · exp((ψ - f₀) / b)`
4. Includes sinh overflow protection (threshold: `f/a > 700`)
5. Clamps result: `θ ≥ θ_min = 1e-30`

---

## Detailed Component Design

### 3.1 Fault Geometry

```cpp
// fault/fault_geometry.hpp

template <typename MeshType = Mesh>
class FaultGeometry
{
public:
   /// Constructor takes domain operator (not mesh + face array)
   /// Extracts fault DOF info and precomputes depth-dependent params
   FaultGeometry(DomainOperator<MeshType> &domain_op, const BP2Params &params);

   int NumFaultDOFs() const { return num_fault_dofs_; }
   const Vector &GetDepths() const { return depths_; }
   const Vector &GetAValues() const { return a_values_; }
   const Vector &GetEtaValues() const { return eta_values_; }
   const BP2Params &GetParams() const { return params_; }

   /// Find DOF index closest to target depth
   int FindNearestDOF(real_t target_depth) const;

   /// Check if depth is velocity-weakening (a(z) < b)
   bool IsVelocityWeakening(real_t z) const;

   /// Get VW/VS transition depths
   real_t GetVWDepth() const { return -params_.H; }
   real_t GetVSDepth() const { return -(params_.H + params_.h); }

   /// Get DOF indices in VW/VS zones
   void GetVWDOFs(Array<int> &vw_dofs) const;
   void GetVSDOFs(Array<int> &vs_dofs) const;

   /// Print fault geometry info
   void Print(std::ostream &os = mfem::out) const;

private:
   BP2Params params_;
   int num_fault_dofs_;
   Vector depths_;       // z-coordinates from domain operator
   Vector a_values_;     // a(z) precomputed via params_.a_of_z()
   Vector eta_values_;   // η (constant for BP2, from params_.eta())

   void ComputeDepthDependentParams();
};
```

Key differences from original design:
- Constructor takes `DomainOperator<MeshType>&` instead of `Mesh& + Array<int>&`
- No `GetFaceAreas()` method (not needed for midpoint DOF design)
- No `IsRateStateRegion()` / `ComputeA()` (moved to `BP2Params`)
- Added `GetEtaValues()`, `FindNearestDOF()`, `GetVWDOFs()`, `GetVSDOFs()`
- Template parameter `MeshType` for serial/parallel support

### 3.2 Friction Law

The friction law uses the θ (theta) state variable form (equivalent to Tandem's ψ form). See Phase 1 for the full `FrictionLaw` interface and `DieterichRuinaFriction` implementation.

Key methods used by the fault operator:
- `FrictionCoefficient(V, theta, a)` - regularized friction f(V, θ)
- `FaultStrength(V, theta, sigma_n, eta, a)` - τ_f = σ_n·f + η·V
- `SolveSlipRate(tau, theta, sigma_n, eta, a)` - Newton-Raphson solve
- `InitialState(tau0, V_init, sigma_n, eta, a)` - compute θ₀ from equilibrium

The `DieterichRuinaFriction` constructor takes a `Constants` struct:
```cpp
DieterichRuinaFriction::Constants fc;
fc.V0 = params.V0;   // 1e-6 m/s
fc.f0 = params.f0;   // 0.6
fc.b = params.b;     // 0.015
fc.Dc = params.Dc;   // 0.004 m
DieterichRuinaFriction friction(fc);
```

### 3.3 State Evolution

State evolution is separated from friction into `AgingLaw` and `SlipLaw` classes (see Phase 1):

```cpp
AgingLaw aging;   // dθ/dt = 1 - V·θ/Dc
SlipLaw slip;     // dθ/dt = -V·θ/Dc · ln(V·θ/Dc)
```

Both implement `StateEvolution::Rate(V, theta, Dc)` used by the fault operator.

### 3.4 Rate-State Fault Operator

```cpp
// fault/rate_state_fault.hpp

template <typename MeshType = Mesh>
class RateStateFaultOperator
{
public:
   /// State layout constants
   static constexpr int StatePerNode = 2;  ///< slip + theta per node
   static constexpr int SlipIndex = 0;     ///< Index of slip in per-node state
   static constexpr int ThetaIndex = 1;    ///< Index of theta in per-node state

   /// Constructor with abstract interfaces
   RateStateFaultOperator(FaultGeometry<MeshType> *geom,
                          FrictionLaw *friction,
                          StateEvolution *evolution,
                          const BP2Params &params);

   // --- State size ---
   int StateSize() const { return num_nodes_ * StatePerNode; }
   int NumNodes() const { return num_nodes_; }

   // --- Initialization (Tandem interface) ---
   void PreInit(Vector &state);                        // Set slip=0, theta=Dc/V_init
   real_t Init(const Vector &traction, Vector &state); // Compute θ₀ from equilibrium

   // --- Time stepping ---
   real_t ComputeRHS(const Vector &traction, const Vector &state, Vector &rate);

   // --- State accessors ---
   void GetSlip(const Vector &state, Vector &slip) const;
   void GetTheta(const Vector &state, Vector &theta) const;
   void SetSlip(const Vector &slip, Vector &state);
   void SetTheta(const Vector &theta, Vector &state);

   // --- Monitoring ---
   const Vector &GetSlipRate() const { return slip_rate_; }
   real_t GetMaxSlipRate() const { return V_max_; }
   real_t GetTau0() const { return tau0_; }

   // --- Component accessors ---
   const FaultGeometry<MeshType> *GetGeometry() const;
   const FrictionLaw *GetFrictionLaw() const;
   const StateEvolution *GetEvolution() const;
   const BP2Params &GetParams() const;

   // --- Verification ---
   real_t VerifyStressEquilibrium(const Vector &traction, const Vector &state,
                                  real_t tol = 1e-10) const;
   void PrintState(const Vector &state, std::ostream &os = mfem::out) const;

private:
   FaultGeometry<MeshType> *geom_;
   FrictionLaw *friction_;
   StateEvolution *evolution_;
   BP2Params params_;
   int num_nodes_;
   real_t tau0_;       // Pre-stress [Pa]
   real_t V_max_;      // Max slip rate from last evaluation
   Vector slip_rate_;  // Cached slip rate [NumNodes()]
};
```

Key design points:
- Uses abstract `FrictionLaw*` and `StateEvolution*` (not concrete `DieterichRuinaAgeing`)
- Template `<typename MeshType>` for serial/parallel support
- `SlipIndex=0`, `ThetaIndex=1` named constants
- Per-node η from `geom_->GetEtaValues()` (not a single value)
- DOFs below `Wf` get prescribed plate rate `Vp` (no state evolution)
- `VerifyStressEquilibrium` checks `τ = σ_n · f(V, θ) + η · V` at all nodes
- `PrintState` outputs summary including slip/theta/V ranges

### 3.5 BP2 Parameters

```cpp
// config/bp2_params.hpp

namespace mfem {
namespace seas {

struct BP2Params {
    // Material properties
    real_t rho = 2670.0;           // Density [kg/m³]
    real_t cs = 3464.0;            // Shear wave speed [m/s]
    real_t mu() const { return rho * cs * cs; }  // Shear modulus [Pa] (~32.04 GPa)

    // Rate-and-state friction parameters
    real_t V0 = 1.0e-6;            // Reference slip rate [m/s]
    real_t f0 = 0.6;               // Reference friction coefficient
    real_t b = 0.015;              // State evolution effect parameter
    real_t Dc = 0.004;             // Critical slip distance [m]
    real_t a0 = 0.010;             // Rate-state a in VW zone
    real_t amax = 0.025;           // Rate-state a in VS zone

    // Stress parameters
    real_t sigma_n = 50.0e6;       // Normal stress [Pa]
    real_t eta() const { return mu() / (2.0 * cs); }  // Radiation damping [Pa·s/m]

    // Loading
    real_t Vp = 1.0e-9;            // Plate rate [m/s]
    real_t V_init = 1.0e-9;        // Initial slip rate [m/s]

    // Geometry
    real_t H = 15000.0;            // VW zone depth [m]
    real_t h = 3000.0;             // VW-VS transition width [m]
    real_t Wf = 40000.0;           // Rate-state fault depth [m]

    // Simulation
    static constexpr real_t seconds_per_year = 365.25 * 24.0 * 3600.0;
    real_t t_final = 1200.0 * seconds_per_year;  // Final time [s] (1200 years)

    // --- Methods ---

    /// Depth-dependent a(z): a0 in VW, linear transition, amax in VS
    real_t a_of_z(real_t z) const;

    /// Check if depth is velocity-weakening (a(z) < b)
    bool IsVelocityWeakening(real_t z) const;

    /// Pre-stress tau0 at steady state with V = V_init in VS zone
    real_t tau0() const;

    /// 12 output station depths for BP2 benchmark [m]
    static void GetOutputDepths(Array<real_t> &depths);

    /// 5 benchmark comparison depths [m]
    static void GetBenchmarkDepths(Array<real_t> &depths);

    /// Print all parameters
    void Print(std::ostream &os = mfem::out) const;
};

} // namespace seas
} // namespace mfem
```

Key differences from original design:
- `mu` is a computed method `mu()` (= `rho * cs * cs`), not a stored field
- `Dc` replaces `L` as the parameter name
- `eta()` = `mu() / (2 * cs)`, not `sqrt(mu * rho) / 2`
- `tau0()` replaces `ComputeTau0()`
- No `Lx` / `Lz` domain geometry fields (those are in `BP2MeshGenerator`)
- `t_final` = 1200 years (not 3000 years)
- `seconds_per_year` static constexpr for unit conversion
- Added `a_of_z()`, `IsVelocityWeakening()`, `GetOutputDepths()`, `GetBenchmarkDepths()`, `Print()` methods

---

## Traction Computation Interface

The fault operator receives traction from the domain operator. See `dg_antiplane_theory.md` Section 6.

### Traction Formula

```
τ = μ · {{∂u/∂x}} = μ · (∂u⁺/∂x + ∂u⁻/∂x) / 2
```

where:
- μ = shear modulus (physical units, Pa)
- ∂u⁺/∂x = gradient from right side of fault (x > 0)
- ∂u⁻/∂x = gradient from left side of fault (x < 0)
- {{·}} = average operator

### Traction Method Signature

The domain operator's `ComputeTraction` takes 3 arguments (displacement, slip BC interpolation, and output traction):

```cpp
void ComputeTraction(const SEASGridFunction<MeshType> &displacement,
                     const Vector &slip_bc,
                     Vector &traction);
```

The `slip_bc` vector contains the fault-boundary-interpolated slip values (computed by `InterpolateSlipBC(slip, slip_bc)`), which is used internally for consistent traction evaluation.

### Traction Data Flow

```
┌──────────────────┐     slip δ      ┌──────────────────┐
│   FaultOperator  │ ──────────────► │  DomainOperator  │
│  (rate-state)    │                 │  (DG Laplace)    │
│                  │ ◄────────────── │                  │
│  V, dθ/dt        │    traction τ   │  u solution      │
└──────────────────┘                 └──────────────────┘
```

---

## Quasi-Dynamic Coupling

### Sequence per Time Step

1. **Extract slip** from fault state
2. **Interpolate slip** to fault boundary DOFs (via `InterpolateSlipBC`)
3. **Solve elastic BVP** with current slip as interior condition
4. **Compute traction** on fault from displacement solution (3-arg with `slip_bc`)
5. **Compute state derivatives** via friction RHS (traction is τ_qs only; τ₀ is added inside `ComputeRHS`)

```cpp
/// Quasi-dynamic RHS evaluation (Phase 4)
void SEASOperator::RHS(real_t time, const Vector &state, Vector &rate)
{
    // 1. Extract slip from state
    fault_op_->GetSlip(state, slip_);

    // 2. Interpolate slip to boundary DOFs
    domain_op_->InterpolateSlipBC(slip_, slip_bc_);

    // 3. Solve domain problem with current slip
    domain_op_->Solve(time, slip_bc_, displacement_);

    // 4. Compute quasi-static traction on fault
    domain_op_->ComputeTraction(displacement_, slip_bc_, traction_);

    // 5. Compute fault state derivatives
    //    RateStateFaultOperator::ComputeRHS adds tau0_ internally:
    //    τ_total = τ₀ + τ_qs (traction from domain solve)
    fault_op_->ComputeRHS(traction_, state, rate);
}
```

---

## Output Stations

BP2 requires output at 12 stations at specific depths, implemented as `BP2Params::GetOutputDepths()`:

| Station | Depth (km) | z coordinate (m) |
|---------|------------|-------------------|
| 1 | 0.0 | 0.0 |
| 2 | 2.4 | -2400.0 |
| 3 | 4.8 | -4800.0 |
| 4 | 7.2 | -7200.0 |
| 5 | 9.6 | -9600.0 |
| 6 | 12.0 | -12000.0 |
| 7 | 14.4 | -14400.0 |
| 8 | 16.8 | -16800.0 |
| 9 | 19.2 | -19200.0 |
| 10 | 24.0 | -24000.0 |
| 11 | 28.8 | -28800.0 |
| 12 | 36.0 | -36000.0 |

A subset of 5 depths with available reference data is provided by `BP2Params::GetBenchmarkDepths()`: z = 0, -4800, -12000, -16800, -24000 m.

The nearest DOF to each station is found using `FaultGeometry::FindNearestDOF(target_depth)`.

---

## Unit Test Specifications

All tests are in `tests/unit/test_fault_operator.cpp`. Uses a custom test framework (not Google Test) with `TEST_ASSERT`, `TEST_NEAR`, `TEST_REL_NEAR` macros.

### Test Groups

**TestPreStress** (3 tests):
- Pre-stress τ₀ ≈ 26.546 MPa (via `params.tau0()`)
- Radiation damping η ≈ 4.63 MPa·s/m (via `params.eta()`)
- Shear modulus μ ≈ 32.04 GPa (via `params.mu()`)

**TestAofZ** (6 tests):
- `a(z=0)` = a0 (VW zone)
- `a(z=-H)` = a0 (top of transition)
- `a(z=-16.5km)` = (a0+amax)/2 (middle of transition)
- `a(z=-18km)` = amax (VS zone)
- `a(z=-30km)` = amax (deep VS)
- `IsVelocityWeakening()` at z=0, -10km (true), -30km (false)

**TestFaultGeometry** (5 tests):
- Positive number of fault DOFs
- Depths in expected range [-Lz, 0]
- All a values match `a_of_z()` profile
- All eta values constant
- `FindNearestDOF(-10km)` returns valid index

**TestRateStateFaultOperator** (6+ tests):
- `StateSize() = 2 * NumFaultDOFs()`
- `PreInit`: all slips zero, all theta positive
- `Init`: τ₀ ≈ 26.546 MPa, V_max > 0, theta in [100, 1e9]
- Stress equilibrium satisfied after Init (error < 1e-8)
- `ComputeRHS`: slip rates > 0, dθ/dt reasonable magnitude
- State access roundtrip (Set/Get slip, theta)

**TestInitialStateValues** (3 tests):
- Initial θ at a=a0: stress equilibrium within 1e-10
- Initial θ at a=amax: stress equilibrium within 1e-10
- Slip rate recovery: SolveSlipRate(τ₀, θ₀) matches V_init within 1e-6

**TestStateEvolutionConsistency** (3 tests):
- Steady state: θ = Dc/V → dθ/dt = 0
- Nearly locked: V ~ 0 → dθ/dt ≈ 1 (healing)
- Fast slip: V = 1 m/s → dθ/dt < 0 (weakening)

---

## File Organization

```
fault/
├── fault_geometry.hpp       # FaultGeometry<MeshType> template class
├── rate_state_fault.hpp     # RateStateFaultOperator<MeshType> template class

friction/
├── friction_law.hpp         # FrictionLaw abstract interface
├── dieterich_ruina.hpp      # DieterichRuinaFriction (regularized, Newton-Raphson solver)
├── state_evolution.hpp      # StateEvolution abstract + AgingLaw, SlipLaw

config/
├── bp2_params.hpp           # BP2Params struct with a_of_z, tau0, eta, etc.

tests/unit/
├── test_fault_operator.cpp  # All Phase 3 tests (~26 tests across 6 groups)
```

**Note:** The slip rate root finder (Newton-Raphson) is implemented inline in `dieterich_ruina.hpp::SolveSlipRate()`. There is no separate root finder file.

---

## Acceptance Criteria

| Component | Tests | Pass Criteria |
|-----------|-------|---------------|
| BP2 Parameters | 3 (TestPreStress) | τ₀ ≈ 26.546 MPa (±0.01), η ≈ 4.63 MPa·s/m, μ ≈ 32.04 GPa |
| Depth profile a(z) | 6 (TestAofZ) | Piecewise linear profile, VW/VS identification |
| Fault geometry | 5 (TestFaultGeometry) | DOFs > 0, depths in range, a/eta match params |
| Fault operator | 6+ (TestRateStateFaultOperator) | PreInit zeros, Init equilibrium, ComputeRHS positive V |
| Initial state | 3 (TestInitialStateValues) | Stress equilibrium < 1e-10, V recovery < 1e-6 |
| State evolution | 3 (TestStateEvolutionConsistency) | Steady state, healing, weakening behaviors |

**Total: ~26 tests across 6 test groups.**

### Expected Initial Values (BP2)

| Parameter | Expected Value |
|-----------|----------------|
| Pre-stress τ₀ | 26.546 MPa |
| Radiation damping η | 4.63 MPa·s/m |
| Shear modulus μ | 32.04 GPa |
| Initial θ (z=0, a=a0) | Computed from τ₀, V_init via `InitialState()` |
| Initial θ (VS zone, a=amax) | Steady-state: Dc/V_init = 4e6 s |
| Initial V | 10⁻⁹ m/s (everywhere) |
