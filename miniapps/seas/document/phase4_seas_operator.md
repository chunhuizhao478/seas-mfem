# Phase 4: SEAS Operator (Serial)

## Overview

This phase implements the SEAS (Sequences of Earthquakes and Aseismic Slip) time-dependent operator that couples the domain solver with the fault operator for quasi-dynamic earthquake cycle simulations.

## Dependencies

- **Phase 2** (`phase2_domain_operator.md`): Domain operator for solving elastic BVP
- **Phase 3** (`phase3_fault_operator.md`): Fault operator for rate-and-state friction
- **Theory**: `dg_antiplane_theory.md` for DG formulation

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 4.1 Quasi-dynamic operator | `solver/seas_operator.hpp` | `tests/unit/test_quasi_dynamic.cpp` |
| 4.2 Time stepper integration | `solver/time_stepper.hpp` | `tests/unit/test_quasi_dynamic.cpp` |
| 4.3 Short simulation test | - | `tests/unit/test_bp2_short.cpp` |

## Verification Checkpoint

1-year simulation runs without crash, output is qualitatively correct.

---

## Quasi-Dynamic Formulation

### Coupling Equations

The quasi-dynamic approach couples:
1. **Domain problem**: Solve Laplace equation with slip BC (Phase 2)
2. **Fault evolution**: Rate-and-state friction with radiation damping (Phase 3)

**Stress balance at fault** (see `phase3_fault_operator.md`):
```
tau = tau0 + tau_qs = sigma_n * f(V, theta) + eta * V
```

where:
- `tau0` = pre-stress (uniform, computed at t=0)
- `tau_qs` = quasi-static shear stress from domain solve (see `dg_antiplane_theory.md` Section 6)
- `eta = mu/(2*cs)` = radiation damping coefficient
- `f(V, theta)` = regularized friction coefficient (see `phase3_fault_operator.md`)

### Operator Coupling Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    SEAS Quasi-Dynamic Operator                          │
│                                                                         │
│  State Vector: [δ₀, θ₀, δ₁, θ₁, ..., δₙ, θₙ]                          │
│                  └─slip─┘  └─state─┘                                   │
│                                                                         │
│  ┌─────────────┐                           ┌─────────────────────────┐ │
│  │ Extract     │      slip δ               │ DomainOperator         │ │
│  │ slip from   ├──────────────────────────►│ (Phase 2)              │ │
│  │ state       │                           │                         │ │
│  └─────────────┘                           │ • Solve: ∇²u = 0       │ │
│                                            │ • Fault BC: [[u]] = δ  │ │
│                                            │ • Bottom: u = Vp·t/2   │ │
│                                            └───────────┬─────────────┘ │
│                                                        │               │
│                                            displacement u              │
│                                                        │               │
│                                                        ▼               │
│  ┌─────────────────────────┐              ┌─────────────────────────┐ │
│  │ FaultOperator           │◄─────────────┤ ComputeTraction()       │ │
│  │ (Phase 3)               │  traction τ  │ τ = μ · {{∂u/∂x}}       │ │
│  │                         │              └─────────────────────────┘ │
│  │ • Stress balance        │                                          │
│  │ • Newton root finder    │                                          │
│  │ • State evolution       │                                          │
│  │                         │                                          │
│  │ Output: V, dθ/dt        │                                          │
│  └───────────┬─────────────┘                                          │
│              │                                                         │
│              ▼                                                         │
│  ┌─────────────────────────┐                                          │
│  │ Rate Vector:            │                                          │
│  │ [V₀, dθ₀/dt, V₁, ...]   │                                          │
│  └─────────────────────────┘                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

**State evolution** (aging law):
```
dθ/dt = 1 - V·θ/Dc

where θ is the physical state variable [s], related to ψ by: ψ = f₀ + b·ln(V₀·θ/Dc)
```

**Slip evolution**:
```
dδ/dt = V
```

---

## ODE System

### State Vector Layout

The state vector contains slip and state variable for each fault DOF:

```
state = [δ₀, θ₀, δ₁, θ₁, δ₂, θ₂, ..., δₙ₋₁, θₙ₋₁]
         └── Fault DOF 0 ──┘  └── DOF 1 ──┘        └── DOF n-1 ──┘
```

Size: 2 × (number of fault DOFs)

### Rate Vector Layout

```
rate = [V₀, dθ₀/dt, V₁, dθ₁/dt, ..., Vₙ₋₁, dθₙ₋₁/dt]
```

where:
- V_i = slip rate at fault DOF i (computed from stress balance)
- dθ_i/dt = state rate at fault DOF i (from aging law: 1 - V·θ/Dc)

### MFEM ODE Solver Integration

The SEAS operator inherits from `TimeDependentOperator`:

```cpp
class SEASQuasiDynamicOperator : public TimeDependentOperator {
public:
    /// Constructor: size = 2 * num_fault_dofs
    SEASQuasiDynamicOperator(DomainOpType *domain,
                              RateStateFaultOperator<MeshType> *fault)
        : TimeDependentOperator(fault->StateSize()),
          domain_(domain), fault_(fault) {}

    /// Main ODE function: compute dy/dt given y
    /// Called by RK4Solver or other ODE solver
    void Mult(const Vector &state, Vector &rate) const override;
};
```

### Available ODE Solvers

| Solver | Order | Adaptive | Notes |
|--------|-------|----------|-------|
| `RK4Solver` (MFEM) | 4 | No | Classic Runge-Kutta, fixed dt |
| `DormandPrinceRK45` (custom) | 5(4) | Yes | DOPRI5 with error-based step control (**recommended**) |

**`RK4Solver`** is used for initial testing with external `AdaptiveTimeStepper` for V-based dt control.

**`DormandPrinceRK45`** is a custom embedded RK integrator in `solver/time_stepper.hpp`:
- 7-stage FSAL (First Same As Last) design
- Error estimated from 5th-order vs 4th-order difference
- Weighted L-infinity error norm with atol/rtol control
- PI controller for step size: `dt_new = safety * dt * err_norm^(-1/5)`
- Designed following Tandem's approach (PETSc TS `-ts_rk_type 5dp`)
- Returns `bool` from `Step()`: true = accepted, false = rejected (caller retries)

---

## Detailed Component Design

### 4.1 Quasi-Dynamic Operator

```cpp
// solver/seas_operator.hpp

namespace mfem {
namespace seas {

/// SEAS quasi-dynamic time integration operator
/// Inherits from TimeDependentOperator for use with MFEM ODE solvers
template <typename MeshType = Mesh>
class SEASQuasiDynamicOperator : public TimeDependentOperator {
public:
    using DomainOpType = AntiplaneDomainOperator<MeshType>;
    using GridFuncType = typename DomainOpType::GridFuncType;

    SEASQuasiDynamicOperator(DomainOpType *domain,
                              RateStateFaultOperator<MeshType> *fault);

    /// Initialize state at t=0
    void SetInitialCondition(Vector &state);

    /// Compute d(state)/dt = RHS(t, state)
    /// state = [delta_0, theta_0, delta_1, theta_1, ...]
    void Mult(const Vector &state, Vector &rate) const override;

    /// Get current displacement solution
    const GridFuncType &GetDisplacement() const { return *u_gf_; }

    /// Get traction at fault
    const Vector &GetTraction() const { return traction_; }

    /// Get maximum slip rate
    real_t GetMaxSlipRate() const { return fault_->GetMaxSlipRate(); }

    /// Get the domain operator
    const DomainOpType *GetDomain() const { return domain_; }

    /// Get the fault operator
    const RateStateFaultOperator<MeshType> *GetFault() const { return fault_; }

private:
    DomainOpType *domain_;
    RateStateFaultOperator<MeshType> *fault_;

    std::unique_ptr<GridFuncType> u_gf_;

    mutable Vector slip_;
    mutable Vector traction_;
};

} // namespace seas
} // namespace mfem
```

### Implementation Details

```cpp
template <typename MeshType>
void SEASQuasiDynamicOperator<MeshType>::Mult(
    const Vector &state, Vector &rate) const
{
    // 1. Extract slip from state vector
    fault_->GetSlip(state, slip_);

    // 2. Solve domain problem with slip BC
    //    ∇²u = 0 with [[u]] = slip on fault (all x=0 faces)
    domain_->Solve(t, slip_, *u_gf_);

    // 3. Compute traction at fault from displacement (3-arg: disp, slip, traction)
    //    τ_qs = μ * {{∂u/∂x}} + μ * κ * h⁻¹ * ([[u]] - δ)
    domain_->ComputeTraction(*u_gf_, slip_, traction_);

    // 4. Compute fault RHS (slip rate and state rate)
    //    Stress balance: τ₀ + τ_qs = σ_n * f(V, θ) + η * V
    //    State evolution: dθ/dt = G(V, θ)
    fault_->ComputeRHS(traction_, state, rate);
}

template <typename MeshType>
void SEASQuasiDynamicOperator<MeshType>::SetInitialCondition(Vector &state)
{
    // Guard: verify state vector size
    MFEM_VERIFY(state.Size() == fault_->StateSize(),
                "State vector size mismatch: got " << state.Size()
                << ", expected " << fault_->StateSize());

    // Phase 1: Pre-initialize (slip=0, theta=placeholder)
    fault_->PreInit(state);

    // Phase 2: Solve domain with zero slip to get initial traction
    fault_->GetSlip(state, slip_);
    domain_->Solve(0.0, slip_, *u_gf_);
    domain_->ComputeTraction(*u_gf_, slip_, traction_);

    // Phase 3: Initialize theta from stress equilibrium
    //   tau0 + traction = sigma_n * f(V_init, theta) + eta * V_init
    real_t V_max = fault_->Init(traction_, state);

    // Phase 4: Verify initial slip rate
    // Re-solve domain and recompute traction with updated state
    fault_->GetSlip(state, slip_);
    domain_->Solve(0.0, slip_, *u_gf_);
    domain_->ComputeTraction(*u_gf_, slip_, traction_);

    // Compute RHS to populate slip rates
    Vector rate_temp(fault_->StateSize());
    fault_->ComputeRHS(traction_, state, rate_temp);
    V_max = fault_->GetMaxSlipRate();

    // Verify stress equilibrium
    real_t eq_error = fault_->VerifyStressEquilibrium(traction_, state);
    MFEM_VERIFY(eq_error < 1e-6,
                "Initial stress equilibrium error too large: " << eq_error);

    // Verify initial slip rate is close to V_init (10% tolerance for coarse meshes)
    const BP2Params &params = fault_->GetParams();
    real_t V_rel_err = std::abs(V_max - params.V_init) /
                       std::max(params.V_init, 1e-30);
    MFEM_VERIFY(V_rel_err < 0.1,
                "Initial V_max differs from V_init by " << V_rel_err * 100 << "%");
}
```

### 4.2 Adaptive Time Stepper

```cpp
// solver/time_stepper.hpp

namespace mfem {
namespace seas {

/// Adaptive time stepper for SEAS simulations
/// Adjusts time step based on maximum slip rate
class AdaptiveTimeStepper {
public:
    AdaptiveTimeStepper();

    /// Configuration setters
    void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }
    void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }
    void SetTargetSlipRateMax(real_t V_target) { V_target_ = V_target; }
    void SetGrowthFactor(real_t factor) { growth_factor_ = factor; }
    void SetInitialDt(real_t dt0) { dt_ = dt0; }

    /// Compute new time step based on slip rate
    real_t ComputeNewDt(real_t V_max, real_t dt_current) const {
        real_t dt_new;

        if (V_max < 1e-12) {
            // Very slow slip: use maximum dt
            dt_new = dt_max_;
        } else if (V_max > V_target_) {
            // Rapid slip: reduce time step
            dt_new = dt_current * V_target_ / V_max;
        } else {
            // Slow slip: can increase time step
            dt_new = std::min(dt_current * growth_factor_, dt_max_);
        }

        // Clamp to bounds
        return std::max(dt_min_, std::min(dt_max_, dt_new));
    }

    /// Get/set current time step
    real_t GetDt() const { return dt_; }
    void SetDt(real_t dt) { dt_ = dt; }

    /// Configuration getters
    real_t GetDtMin() const { return dt_min_; }
    real_t GetDtMax() const { return dt_max_; }
    real_t GetTargetSlipRateMax() const { return V_target_; }

private:
    real_t dt_min_ = 1e-6;          // Minimum time step (seconds)
    real_t dt_max_ = 3.15e7;        // Maximum time step (1 year)
    real_t V_target_ = 1e-6;        // Target max slip rate for dt control
    real_t growth_factor_ = 1.5;    // Factor for increasing dt during slow slip
    real_t dt_ = 1e3;               // Current time step
};

} // namespace seas
} // namespace mfem
```

### 4.2b Dormand-Prince RK45 (Error-Based Adaptive)

```cpp
// solver/time_stepper.hpp

class DormandPrinceRK45 {
public:
    DormandPrinceRK45();

    // Configuration
    void SetAbsTol(real_t atol);    // Default: 1e-7
    void SetRelTol(real_t rtol);    // Default: 1e-7
    void SetSafety(real_t safety);  // Default: 0.9
    void SetGrowthMax(real_t gmax); // Default: 5.0
    void SetShrinkMin(real_t smin); // Default: 0.2
    void SetDtMin(real_t dt_min);   // Default: 1e-6
    void SetDtMax(real_t dt_max);   // Default: 0.1 year
    void SetDt(real_t dt);          // Default: 1e3

    /// Allocate stage vectors
    void Init(TimeDependentOperator &op);

    /// Take one adaptive step. Returns true if accepted, false if rejected.
    /// On rejection, dt_ is shrunk and caller should retry without advancing step counter.
    bool Step(TimeDependentOperator &op, Vector &state, real_t &t, real_t &dt);

    // Accessors
    real_t GetDt() const;
    int GetTotalRejections() const;

private:
    Vector k_[7];     // 7 stage vectors (FSAL)
    Vector y_tmp_;    // Temporary solution (5th-order)
    Vector err_;      // Error estimate (b - b*)
    bool initialized_; // Whether k_[0] is valid from previous step

    // DOPRI5(4) Butcher tableau coefficients (static constexpr)
    // ... (standard Dormand-Prince coefficients)
};
```

Key features:
- **FSAL**: `k_[6]` from accepted step becomes `k_[0]` for next step (saves 1 evaluation)
- **Error norm**: Weighted L-infinity: `max(|err_i| / (atol + rtol * |y_i|))`
- **Step control**: `dt_new = safety * dt * err_norm^(-1/5)`, clamped by growth/shrink factors
- **Diagnostics**: Logs when dt stuck near dt_min (warns about stiff DOFs)

### Time Step Control Strategy

**`AdaptiveTimeStepper`** (V-based, for use with `RK4Solver`):

| Slip Rate (m/s) | Time Step | Description |
|-----------------|-----------|-------------|
| V < 10^-12 | dt_max (1 year) | Very slow interseismic |
| 10^-12 < V < 10^-9 | Gradual increase | Slow interseismic |
| 10^-9 < V < 10^-6 | ~dt_current | Normal interseismic |
| 10^-6 < V < 10^-3 | Decrease | Nucleation |
| V > 10^-3 | dt_min | Coseismic |

**`DormandPrinceRK45`** (error-based): Step size is controlled automatically by the local truncation error estimate. No explicit V-based control needed.

---

## Integration Tests

All tests use a `TestComponents` helper struct that creates mesh (via `BP2MeshGenerator::Create`), domain operator, fault geometry, friction law, state evolution, fault operator, and SEAS operator. Uses custom test framework (`TEST_ASSERT`, `TEST_NEAR`, `TEST_REL_NEAR` macros).

### test_quasi_dynamic.cpp (19 tests)

**Operator Construction and Initialization:**
1. `TestOperatorConstruction` — State size = 2 * num_fault_dofs
2. `TestInitialCondition` — Slip=0, theta>0, V_max ~ V_init, stress equilibrium < 1e-6
3. `TestMultDeterministic` — Same input gives same output

**Time Stepping (RK4):**
4. `TestSteadyStateSlip` — V ~ V_init after 10 steps of 1000s (5% tolerance)
5. `TestStressBalanceMaintained` — Stress equilibrium < 1e-4 during 5 steps
6. `TestSlipConservation` — Slip within order-of-magnitude of V_init*t
7. `TestThetaPositive` — θ remains positive over 20 steps
8. `TestRK4WithAdaptiveDt` — Runs to t=1e5 s with AdaptiveTimeStepper

**Below-Wf Loading:**
9. `TestExtendedFaultFullDepth` — Fault extends below Wf, has DOFs in both zones
10. `TestPrescribedVpBelowWf` — dslip/dt = Vp, dθ/dt = 0 below Wf
11. `TestSlipAccumulatesAtVp` — Below-Wf slip = Vp*t (RK4 exact for linear)
12. `TestUniformSlipDisplacement` — Uniform slip → ±δ/2 displacement
13. `TestBelowWfSlipProducesTraction` — Below-Wf slip generates traction above Wf

**Adaptive Time Stepper:**
14. `TestAdaptiveTimeStepper` — V-based dt control: slow→dt_max, rapid→reduce, bounds

**Dormand-Prince RK45:**
15. `TestRK45ExponentialDecay` — dy/dt = -λy, accuracy < 1e-5
16. `TestRK45HarmonicOscillator` — Energy conservation < 1e-6 over 5 periods
17. `TestRK45StepRejection` — Logistic growth: confirms step rejections occur
18. `TestRK45FSAL` — FSAL property gives accurate result over 10 steps
19. `TestRK45ToleranceControl` — Tighter tolerance → smaller error, more steps

### test_bp2_short.cpp (6 tests)

**1-year BP2 simulation** with RK4 + AdaptiveTimeStepper on coarse mesh (4×8, order 1):
1. Initial V_max ~ V_init (10% tolerance)
2. Simulation completes within step limit
3. No coseismic events (V_max < 1e-3)
4. θ remains positive throughout
5. Stress in [15, 40] MPa range
6. Mean slip > 0 and within 10x of V_init*t

---

## Main Driver Structure (Planned)

**Note:** No `seas.cpp` main driver exists yet. The following shows the intended structure based on the test code patterns:

```cpp
// seas.cpp (planned - Phase 5)

#include "domain/antiplane_operator.hpp"
#include "fault/rate_state_fault.hpp"
#include "solver/seas_operator.hpp"
#include "solver/time_stepper.hpp"

int main(int argc, char *argv[]) {
    BP2Params params;

    // Create BP2 mesh
    BP2MeshGenerator::Parameters mesh_params;
    mesh_params.Lx = 50.0e3;
    mesh_params.Lz = 50.0e3;
    mesh_params.Wf = params.Wf;
    mesh_params.nx = 4;
    mesh_params.nz = 8;
    auto mesh = BP2MeshGenerator::Create(mesh_params);

    // Create domain operator
    AntiplaneDomainOperator<Mesh> domain(*mesh, 1, params.mu(), params.Vp, params.Wf);

    // Create fault components
    FaultGeometry<Mesh> fault_geom(domain, params);

    DieterichRuinaFriction::Constants fc;
    fc.V0 = params.V0;  fc.f0 = params.f0;
    fc.b  = params.b;   fc.Dc = params.Dc;
    DieterichRuinaFriction friction(fc);
    AgingLaw aging;

    RateStateFaultOperator<Mesh> fault_op(&fault_geom, &friction, &aging, params);

    // Create SEAS operator
    SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

    // Initialize state
    Vector state(fault_op.StateSize());
    seas_op.SetInitialCondition(state);

    // Option A: RK4 + AdaptiveTimeStepper (V-based dt control)
    AdaptiveTimeStepper stepper;
    RK4Solver ode_solver;
    ode_solver.Init(seas_op);

    real_t t = 0.0;
    while (t < params.t_final) {
        real_t dt = stepper.GetDt();
        if (t + dt > params.t_final) { dt = params.t_final - t; }
        ode_solver.Step(state, t, dt);
        stepper.SetDt(stepper.ComputeNewDt(seas_op.GetMaxSlipRate(), dt));
    }

    // Option B: DormandPrinceRK45 (error-based adaptive)
    // DormandPrinceRK45 rk45;
    // rk45.SetAbsTol(1e-7);  rk45.SetRelTol(1e-7);
    // rk45.Init(seas_op);
    // while (t < params.t_final) {
    //     real_t dt;
    //     if (rk45.Step(seas_op, state, t, dt)) { /* accepted */ }
    // }

    return 0;
}
```

---

## Expected Behavior

### Earthquake Cycle Phases

1. **Interseismic** (majority of time):
   - V ~ V_init = 10^-9 m/s
   - theta slowly increasing
   - tau slowly increasing due to loading

2. **Nucleation** (brief, before earthquake):
   - V accelerates (10^-9 to 10^-3 m/s)
   - theta drops rapidly
   - Localized acceleration on fault

3. **Coseismic** (very brief, ~seconds):
   - V ~ 1 m/s
   - theta very small
   - Rapid stress drop

4. **Postseismic** (days to years):
   - V decays from high values
   - theta recovers
   - Afterslip

### Expected Timing

| Event | Expected Time |
|-------|---------------|
| First earthquake | ~90 years |
| Recurrence interval | ~90 years |
| Coseismic duration | ~10 seconds |

---

## File Organization

```
solver/
├── seas_operator.hpp          # SEASQuasiDynamicOperator<MeshType> (TimeDependentOperator)
└── time_stepper.hpp           # AdaptiveTimeStepper + DormandPrinceRK45

tests/unit/
├── test_quasi_dynamic.cpp     # 19 tests: coupling, time stepping, below-Wf, RK45
└── test_bp2_short.cpp         # 6 tests: 1-year BP2 simulation
```

---

## Acceptance Criteria

| Component | Tests | Pass Criteria |
|-----------|-------|---------------|
| Operator construction | 2 (TestOperatorConstruction) | State size correct |
| Initial condition | 4 (TestInitialCondition) | Slip=0, θ>0, V~V_init, equilibrium<1e-6 |
| Steady-state | 2 (TestSteadyStateSlip) | V within 5% of V_init for 10,000s |
| Stress balance | 1 (TestStressBalanceMaintained) | <1e-4 relative during stepping |
| Conservation | 1 (TestSlipConservation) | Slip within order-of-magnitude of V_init*t |
| θ positive | 1 (TestThetaPositive) | θ>0 over 20 steps |
| Adaptive dt | 5 (TestAdaptiveTimeStepper) | V-based control, bounds respected |
| RK4+adaptive | 4 (TestRK4WithAdaptiveDt) | Completes, V/θ reasonable |
| Determinism | 1 (TestMultDeterministic) | Same input→same output |
| Below-Wf | 5 tests | Fault extends below Wf, prescribed V=Vp, traction transfer |
| Displacement | 4 (TestUniformSlipDisplacement) | ±δ/2 symmetry |
| RK45 accuracy | 5 tests (exponential, oscillator, logistic, FSAL, tolerance) | Error < tolerance |
| BP2 short | 6 (TestBP2ShortSimulation) | 1-year run, stress [15,40] MPa, V<1e-3 |

**Total: ~25 tests in `test_quasi_dynamic.cpp` + 6 tests in `test_bp2_short.cpp`.**

Unit tests use coarse meshes (4×8 or 4×10, order 1) with relaxed tolerances suitable for CI.
