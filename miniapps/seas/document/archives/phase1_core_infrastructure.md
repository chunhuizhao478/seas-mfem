# Phase 1: Core Infrastructure (Serial)

## Overview

This phase establishes the foundational components for rate-and-state friction computation. All components are implemented in serial mode first for concept verification.

## Dependencies

- None (first phase)

## Dependent Phases

- **Phase 3** (`phase3_fault_operator.md`): Uses friction law and state evolution from this phase
- **Phase 4** (`phase4_seas_operator.md`): Uses initial state computation

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 1.1 Friction law interface | `friction_law.hpp` | `test_friction_law.cpp` |
| 1.2 Dieterich-Ruina implementation | `dieterich_ruina.hpp` | Included above |
| 1.3 State evolution laws | `state_evolution.hpp` | `test_state_evolution.cpp` |
| 1.4 Slip rate solver | In `dieterich_ruina.hpp` | `test_slip_rate_solver.cpp` |
| 1.5 Initial state computation | In `dieterich_ruina.hpp` | `test_initial_state.cpp` |

## Verification Checkpoint

All unit tests pass for friction components.

---

## Detailed Component Design

### 1.1 Friction Law Interface

```cpp
// friction/friction_law.hpp

class FrictionLaw {
public:
    virtual ~FrictionLaw() = default;

    /// Compute friction coefficient f(V, theta)
    /// Uses scalar 'a' instead of Vector params for efficiency.
    /// For BP2, only 'a' varies spatially; other parameters
    /// (b, Dc, V0, f0) are stored in the Constants struct.
    virtual real_t FrictionCoefficient(real_t V, real_t theta, real_t a) const = 0;

    /// Compute fault strength F = sigma_n * f(V, theta)
    /// Default implementation provided (derived classes need not override)
    virtual real_t FaultStrength(real_t V, real_t theta, real_t sigma_n,
                                  real_t a) const
    {
        return sigma_n * FrictionCoefficient(V, theta, a);
    }

    /// Solve for slip rate V given stress tau and state theta
    /// Solves: tau = sigma_n * f(V, theta) + eta * V
    /// Optional iterations output for convergence diagnostics
    virtual real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                                  real_t eta, real_t a,
                                  int *iterations = nullptr) const = 0;

    /// Compute initial state theta(0) from equilibrium
    virtual real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                                 real_t eta, real_t a) const = 0;

    /// Compute df/dV for Newton solver and implicit time stepping
    virtual real_t FrictionDerivativeV(real_t V, real_t theta, real_t a) const = 0;

    /// Compute df/dtheta for implicit time stepping
    virtual real_t FrictionDerivativeTheta(real_t V, real_t theta,
                                           real_t a) const = 0;

    /// Parameter accessors
    virtual real_t GetV0() const = 0;
    virtual real_t GetF0() const = 0;
    virtual real_t GetB() const = 0;
    virtual real_t GetDc() const = 0;
};
```

### 1.2 Dieterich-Ruina Implementation

```cpp
// friction/dieterich_ruina.hpp

class DieterichRuinaFriction : public FrictionLaw {
public:
    // Global constants
    struct Constants {
        real_t V0 = 1.0e-6;    // Reference slip rate [m/s]
        real_t f0 = 0.6;       // Reference friction coefficient
        real_t b  = 0.015;     // State evolution parameter
        real_t Dc = 0.004;     // Critical slip distance [m] (BP2 value)
    };

    explicit DieterichRuinaFriction(const Constants &c) : cp_(c) {}

    /// Default constructor with BP2 default values
    DieterichRuinaFriction()
        : DieterichRuinaFriction(Constants{1.0e-6, 0.6, 0.015, 0.004}) {}

    /// f(V, theta) = a * asinh[(V / 2V0) * exp((f0 + b*ln(V0*theta/Dc)) / a)]
    /// V and theta are clamped to V_min_=1e-30 and theta_min_=1e-30
    real_t FrictionCoefficient(real_t V, real_t theta, real_t a) const override {
        V = std::max(V, V_min_);
        theta = std::max(theta, theta_min_);
        real_t log_arg = cp_.V0 * theta / cp_.Dc;
        real_t exp_arg = (cp_.f0 + cp_.b * std::log(log_arg)) / a;
        real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(exp_arg);
        return a * std::asinh(sinh_arg);
    }

    /// df/dV using chain rule on asinh
    real_t FrictionDerivativeV(real_t V, real_t theta, real_t a) const override;

    /// df/dtheta using chain rule on asinh
    real_t FrictionDerivativeTheta(real_t V, real_t theta, real_t a) const override;

    /// Solve tau = sigma_n * f(V, theta) + eta * V for V
    /// Uses Newton-Raphson with bisection-style bound safeguards.
    /// Bounds: V_lo = 1e-30, V_hi = min(tau/eta, 100 m/s).
    /// Convergence: rel_change < 1e-12 or |F| < tol*tau, max_iter = 100.
    /// Handles sigma_n <= 0: returns tau/eta (viscous sliding) or 0.
    real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                          real_t eta, real_t a,
                          int *iterations = nullptr) const override;

    /// Compute theta(0) such that tau0 = sigma_n * f(V_init, theta(0)) + eta * V_init
    /// Includes sinh overflow protection: uses exp(x)/2 approximation when f/a > 700
    real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                         real_t eta, real_t a) const override;

    /// Parameter accessors
    real_t GetV0() const override { return cp_.V0; }
    real_t GetF0() const override { return cp_.f0; }
    real_t GetB() const override { return cp_.b; }
    real_t GetDc() const override { return cp_.Dc; }

    /// Constants struct access
    const Constants &GetConstants() const { return cp_; }
    void SetConstants(const Constants &c) { cp_ = c; }

private:
    Constants cp_;

    static constexpr real_t V_min_ = 1.0e-30;      // Minimum slip rate
    static constexpr real_t theta_min_ = 1.0e-30;   // Minimum state variable
};
```

### 1.3 State Evolution Laws

```cpp
// friction/state_evolution.hpp

class StateEvolution {
public:
    virtual ~StateEvolution() = default;

    /// Compute d(theta)/dt given current V and theta
    virtual real_t Rate(real_t V, real_t theta, real_t Dc) const = 0;

    /// Compute steady-state theta for given V
    virtual real_t SteadyState(real_t V, real_t Dc) const = 0;

    /// Compute dG/dV for implicit time stepping
    virtual real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const = 0;

    /// Compute dG/dtheta for implicit time stepping
    virtual real_t RateDerivativeTheta(real_t V, real_t theta, real_t Dc) const = 0;

    /// Get the name of this evolution law for output/logging
    virtual const char *GetName() const = 0;
};

/// Aging law: d(theta)/dt = 1 - V*theta/Dc
class AgingLaw : public StateEvolution {
public:
    real_t Rate(real_t V, real_t theta, real_t Dc) const override {
        return 1.0 - V * theta / Dc;
    }

    real_t SteadyState(real_t V, real_t Dc) const override {
        MFEM_ASSERT(V > 0.0, "Slip rate must be positive for steady state");
        return Dc / V;
    }

    real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const override {
        return -theta / Dc;
    }

    real_t RateDerivativeTheta(real_t V, real_t theta, real_t Dc) const override {
        return -V / Dc;
    }

    const char *GetName() const override { return "AgingLaw"; }
};

/// Slip law: d(theta)/dt = -V*theta/Dc * ln(V*theta/Dc)
/// Includes log(0) guard: x = max(V*theta/Dc, 1e-50)
class SlipLaw : public StateEvolution {
public:
    real_t Rate(real_t V, real_t theta, real_t Dc) const override {
        real_t x = V * theta / Dc;
        if (x < 1e-50) { x = 1e-50; }
        return -x * std::log(x);
    }

    real_t SteadyState(real_t V, real_t Dc) const override {
        MFEM_ASSERT(V > 0.0, "Slip rate must be positive for steady state");
        return Dc / V;
    }

    real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const override {
        real_t x = V * theta / Dc;
        if (x < 1e-50) { x = 1e-50; }
        return -(theta / Dc) * (std::log(x) + 1.0);
    }

    real_t RateDerivativeTheta(real_t V, real_t theta, real_t Dc) const override {
        real_t x = V * theta / Dc;
        if (x < 1e-50) { x = 1e-50; }
        return -(V / Dc) * (std::log(x) + 1.0);
    }

    const char *GetName() const override { return "SlipLaw"; }
};
```

---

## Unit Test Specifications

### Friction Law Tests

```cpp
// tests/unit/test_friction_law.cpp

/// Test 1: Friction coefficient at reference state
/// At V = V0, theta = Dc/V0 (steady state), f should equal f0
TEST(FrictionLaw, ReferenceState) {
    DieterichRuinaFriction::Constants c;
    c.V0 = 1e-6; c.f0 = 0.6; c.b = 0.015; c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    real_t a = 0.015;  // a = b (neutral)
    real_t V = c.V0;
    real_t theta = c.Dc / c.V0;  // steady state

    real_t f = law.FrictionCoefficient(V, theta, a);

    // At reference: f = a * asinh[(V0/2V0) * exp((f0 + b*ln(1))/a)]
    //             = a * asinh[0.5 * exp(f0/a)]
    real_t expected = a * std::asinh(0.5 * std::exp(c.f0 / a));
    EXPECT_NEAR(f, expected, 1e-12);
}

/// Test 2: Friction increases with slip rate (direct effect)
TEST(FrictionLaw, DirectEffect) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    real_t a = 0.010;
    real_t theta = 1000.0;  // Fixed state

    real_t V1 = 1e-9, V2 = 1e-6;
    real_t f1 = law.FrictionCoefficient(V1, theta, a);
    real_t f2 = law.FrictionCoefficient(V2, theta, a);

    EXPECT_GT(f2, f1);  // Higher V -> higher f (at fixed theta)
}

/// Test 3: Friction decreases with state (evolution effect)
TEST(FrictionLaw, EvolutionEffect) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    real_t a = 0.010;
    real_t V = 1e-9;

    real_t theta1 = 100.0, theta2 = 10000.0;
    real_t f1 = law.FrictionCoefficient(V, theta1, params);
    real_t f2 = law.FrictionCoefficient(V, theta2, params);

    EXPECT_GT(f2, f1);  // Higher theta -> higher f (at fixed V)
}

/// Test 4: Velocity-weakening in nucleation zone (a < b)
TEST(FrictionLaw, VelocityWeakening) {
    DieterichRuinaFriction::Constants c;
    c.b = 0.015;
    DieterichRuinaFriction law(c);

    real_t a = 0.010;  // a < b

    real_t V1 = 1e-10, V2 = 1e-8;
    real_t theta1_ss = c.Dc / V1;
    real_t theta2_ss = c.Dc / V2;

    real_t f1 = law.FrictionCoefficient(V1, theta1_ss, a);
    real_t f2 = law.FrictionCoefficient(V2, theta2_ss, a);

    // For a < b, steady-state friction decreases with V
    EXPECT_GT(f1, f2);
}

/// Test 5: Velocity-strengthening in creeping zone (a > b)
TEST(FrictionLaw, VelocityStrengthening) {
    DieterichRuinaFriction::Constants c;
    c.b = 0.015;
    DieterichRuinaFriction law(c);

    real_t a = 0.025;  // a > b

    real_t V1 = 1e-10, V2 = 1e-8;
    real_t theta1_ss = c.Dc / V1;
    real_t theta2_ss = c.Dc / V2;

    real_t f1 = law.FrictionCoefficient(V1, theta1_ss, a);
    real_t f2 = law.FrictionCoefficient(V2, theta2_ss, a);

    // For a > b, steady-state friction increases with V
    EXPECT_LT(f1, f2);
}
```

### Slip Rate Solver Tests

```cpp
// tests/unit/test_slip_rate_solver.cpp

/// Test 1: Solver recovers known slip rate
TEST(SlipRateSolver, RecoverKnownRate) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    real_t a = 0.015;
    real_t sigma_n = 50e6;  // 50 MPa
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);  // mu/(2cs)

    // Pick a V and theta, compute tau, then recover V
    real_t V_true = 1e-8;
    real_t theta = 5000.0;

    real_t f = law.FrictionCoefficient(V_true, theta, a);
    real_t tau = sigma_n * f + eta * V_true;

    real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a);

    EXPECT_NEAR(V_solved, V_true, V_true * 1e-10);
}

/// Test 2: Solver handles zero radiation damping
TEST(SlipRateSolver, ZeroRadiationDamping) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    real_t a = 0.015;
    real_t sigma_n = 50e6;
    real_t eta = 0.0;  // No radiation damping

    real_t V_true = 1e-7;
    real_t theta = 3000.0;

    real_t f = law.FrictionCoefficient(V_true, theta, a);
    real_t tau = sigma_n * f;

    real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a);

    EXPECT_NEAR(V_solved, V_true, V_true * 1e-10);
}

/// Test 3: Solver handles extreme slip rates
TEST(SlipRateSolver, ExtremeRates) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    real_t a = 0.015;
    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);

    // Very slow (interseismic)
    real_t V_slow = 1e-12;
    real_t theta_slow = 1e8;
    real_t f_slow = law.FrictionCoefficient(V_slow, theta_slow, a);
    real_t tau_slow = sigma_n * f_slow + eta * V_slow;
    real_t V_solved_slow = law.SolveSlipRate(tau_slow, theta_slow, sigma_n, eta, a);
    EXPECT_NEAR(V_solved_slow, V_slow, V_slow * 1e-8);

    // Fast (coseismic)
    real_t V_fast = 1.0;  // 1 m/s
    real_t theta_fast = 0.01;
    real_t f_fast = law.FrictionCoefficient(V_fast, theta_fast, a);
    real_t tau_fast = sigma_n * f_fast + eta * V_fast;
    real_t V_solved_fast = law.SolveSlipRate(tau_fast, theta_fast, sigma_n, eta, a);
    EXPECT_NEAR(V_solved_fast, V_fast, V_fast * 1e-8);
}

/// Test 4: Solver convergence count
TEST(SlipRateSolver, ConvergenceEfficiency) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    // Should converge in < 50 iterations for any reasonable input
    // (Implementation should track iteration count)
}
```

### Initial State Tests

```cpp
// tests/unit/test_initial_state.cpp

/// Test 1: Initial state satisfies stress balance
TEST(InitialState, StressBalance) {
    DieterichRuinaFriction::Constants c;
    c.Dc = 0.004;  // BP2 value
    DieterichRuinaFriction law(c);

    real_t a = 0.015;  // a = b (neutral)
    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // Compute pre-stress tau0 (BP2 formula at amax)
    real_t amax = 0.025;
    real_t theta_ss = c.Dc / V_init;
    real_t f_init = law.FrictionCoefficient(V_init, theta_ss, amax);
    real_t tau0 = sigma_n * f_init + eta * V_init;

    // Compute initial state
    real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);

    // Verify: tau0 = sigma_n * f(V_init, theta0) + eta * V_init
    real_t f_check = law.FrictionCoefficient(V_init, theta0, a);
    real_t tau_check = sigma_n * f_check + eta * V_init;

    EXPECT_NEAR(tau_check, tau0, tau0 * 1e-10);
}

/// Test 2: Initial state varies with depth (a varies)
TEST(InitialState, DepthDependence) {
    DieterichRuinaFriction::Constants c;
    c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // Compute tau0 at amax
    real_t amax = 0.025;
    real_t theta_ss = c.Dc / V_init;
    real_t f_max = law.FrictionCoefficient(V_init, theta_ss, amax);
    real_t tau0 = sigma_n * f_max + eta * V_init;

    // Compute theta0 at different depths (different a values)
    std::vector<real_t> a_values = {0.010, 0.015, 0.020, 0.025};
    std::vector<real_t> theta0_values;

    for (real_t a : a_values) {
        real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);
        theta0_values.push_back(theta0);

        // Verify stress balance
        real_t f_check = law.FrictionCoefficient(V_init, theta0, a);
        real_t tau_check = sigma_n * f_check + eta * V_init;
        EXPECT_NEAR(tau_check, tau0, tau0 * 1e-10);
    }

    // theta0 should vary with a
    EXPECT_NE(theta0_values[0], theta0_values[3]);
}

/// Test 3: BP2-specific initial state values
TEST(InitialState, BP2Values) {
    // Compare with benchmark data initial values
    // At z=0: state_log10 = 3.602... -> theta ~ 4000 s

    DieterichRuinaFriction::Constants c;
    c.V0 = 1e-6; c.f0 = 0.6; c.b = 0.015; c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // tau0 from benchmark: shear_stress ~ 26.546 MPa
    real_t tau0_expected = 26.546e6;

    // Compute using our formulas
    real_t amax = 0.025;
    real_t theta_ss_max = c.Dc / V_init;
    real_t f_max = law.FrictionCoefficient(V_init, theta_ss_max, amax);
    real_t tau0_computed = sigma_n * f_max + eta * V_init;

    // tau0 should match benchmark
    EXPECT_NEAR(tau0_computed / 1e6, tau0_expected / 1e6, 0.01);

    // Compute theta0 at a = a0 = 0.010 (surface)
    real_t a0 = 0.010;
    real_t theta0_computed = law.InitialState(tau0_computed, V_init, sigma_n, eta, a0);

    // theta0 should be in reasonable range
    EXPECT_GT(theta0_computed, 1000.0);
    EXPECT_LT(theta0_computed, 10000.0);
}
```

### State Evolution Tests

```cpp
// tests/unit/test_state_evolution.cpp

/// Test 1: Aging law at steady state
TEST(StateEvolution, AgingLawSteadyState) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta_ss = aging.SteadyState(V, Dc);

    // At steady state: d(theta)/dt = 0
    real_t rate = aging.Rate(V, theta_ss, Dc);
    EXPECT_NEAR(rate, 0.0, 1e-20);
}

/// Test 2: Aging law increases theta when V*theta < Dc
TEST(StateEvolution, AgingLawIncreasing) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta = 100.0;  // V*theta = 1e-7 < Dc = 0.004

    real_t rate = aging.Rate(V, theta, Dc);
    EXPECT_GT(rate, 0.0);  // theta should increase
}

/// Test 3: Aging law decreases theta when V*theta > Dc
TEST(StateEvolution, AgingLawDecreasing) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1.0;        // 1 m/s (coseismic)
    real_t theta = 1000.0; // V*theta = 1000 >> Dc = 0.004

    real_t rate = aging.Rate(V, theta, Dc);
    EXPECT_LT(rate, 0.0);  // theta should decrease rapidly
}

/// Test 4: Slip law at steady state
TEST(StateEvolution, SlipLawSteadyState) {
    SlipLaw slip;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta_ss = slip.SteadyState(V, Dc);

    real_t rate = slip.Rate(V, theta_ss, Dc);
    EXPECT_NEAR(rate, 0.0, 1e-20);
}

/// Test 5: Evolution timescale
TEST(StateEvolution, Timescale) {
    AgingLaw aging;
    real_t Dc = 0.004;
    real_t V = 1e-9;

    // Characteristic timescale = Dc/V
    real_t timescale = Dc / V;
    EXPECT_NEAR(timescale, 4e6, 1.0);  // ~46 days
}
```

---

## Acceptance Criteria

| Component | Required Tests | Pass Criteria |
|-----------|----------------|---------------|
| Friction coefficient | 5+ tests | All pass with tolerance 1e-10 |
| Slip rate solver | 4+ tests | Convergence in <50 iterations |
| Initial state | 3+ tests | Stress balance within 1e-10 |
| State evolution | 5+ tests | All pass |

---

## BP2 Key Formulas Reference

**Pre-stress**:
```
tau0 = sigma_n * amax * asinh[(V_init/2V0) * exp((f0+b*ln(V0/V_init))/amax)] + eta*V_init
     ~ 26.546 MPa
```

**Initial state**:
```
theta(z,0) = (Dc/V0) * exp{(a(z)/b) * ln[(2V0/V_init) * sinh((tau0-eta*V_init)/(a(z)*sigma_n))] - f0/b}
```

**Radiation damping**:
```
eta = mu/(2cs) = 32.04e9 / (2 * 3464) ~ 4.625e6 Pa*s/m
```
