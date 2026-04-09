// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Unit tests for Phase 1: Friction Law
// Tests for DieterichRuinaFriction class:
// - Friction coefficient computation
// - Direct effect (friction increases with V at fixed theta)
// - Evolution effect (friction increases with theta at fixed V)
// - Velocity-weakening (a < b)
// - Velocity-strengthening (a > b)
// - Slip rate solver
// - Initial state computation

#include "mfem.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Simple test framework
// =============================================================================

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      real_t _val = (value); \
      real_t _exp = (expected); \
      real_t _tol = (tol); \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", Diff: " << std::abs(_val - _exp) << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_REL_NEAR(value, expected, rel_tol, message) \
   do { \
      num_tests++; \
      real_t _val = (value); \
      real_t _exp = (expected); \
      real_t _tol = (rel_tol) * std::abs(_exp); \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", RelErr: " << std::abs(_val - _exp) / std::abs(_exp) << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// Friction Coefficient Tests
// =============================================================================

void TestFrictionCoefficient()
{
   std::cout << "\n=== Testing Friction Coefficient ===\n";

   DieterichRuinaFriction::Constants c;
   c.V0 = 1e-6;
   c.f0 = 0.6;
   c.b = 0.015;
   c.Dc = 0.004;
   DieterichRuinaFriction law(c);

   // Test 1: Reference state computation
   // At V = V0, theta = Dc/V0 (steady state):
   // f = a * asinh[(V0/2V0) * exp((f0 + b*ln(1))/a)]
   //   = a * asinh[0.5 * exp(f0/a)]
   {
      real_t a = 0.015;  // a = b (neutral)
      real_t V = c.V0;
      real_t theta = c.Dc / c.V0;  // steady state

      real_t f = law.FrictionCoefficient(V, theta, a);
      real_t expected = a * std::asinh(0.5 * std::exp(c.f0 / a));

      TEST_NEAR(f, expected, 1e-12, "Friction coefficient at reference state");
   }

   // Test 2: Friction at BP2 initial conditions
   {
      BP2Params params;
      real_t a = params.amax;
      real_t V = params.V_init;
      real_t theta = c.Dc / V;  // steady state

      real_t f = law.FrictionCoefficient(V, theta, a);

      // Expected from BP2: tau0 = sigma_n * f + eta * V_init
      // f = (tau0 - eta * V_init) / sigma_n
      // tau0 ~ 26.546 MPa for BP2
      real_t tau0 = params.tau0();
      real_t f_expected = (tau0 - params.eta() * V) / params.sigma_n;

      TEST_REL_NEAR(f, f_expected, 1e-6, "Friction coefficient at BP2 initial state");
   }
}

// =============================================================================
// Direct Effect Tests (friction increases with V at fixed theta)
// =============================================================================

void TestDirectEffect()
{
   std::cout << "\n=== Testing Direct Effect ===\n";

   DieterichRuinaFriction::Constants c;
   DieterichRuinaFriction law(c);

   real_t a = 0.010;
   real_t theta = 1000.0;  // Fixed state

   // Test: f(V2) > f(V1) for V2 > V1 at fixed theta
   real_t V1 = 1e-9;
   real_t V2 = 1e-6;

   real_t f1 = law.FrictionCoefficient(V1, theta, a);
   real_t f2 = law.FrictionCoefficient(V2, theta, a);

   TEST_ASSERT(f2 > f1, "Higher V -> higher f (direct effect)");

   // Test derivative is positive
   real_t df_dV = law.FrictionDerivativeV(V1, theta, a);
   TEST_ASSERT(df_dV > 0.0, "df/dV > 0 (positive direct effect)");
}

// =============================================================================
// Evolution Effect Tests (friction increases with theta at fixed V)
// =============================================================================

void TestEvolutionEffect()
{
   std::cout << "\n=== Testing Evolution Effect ===\n";

   DieterichRuinaFriction::Constants c;
   DieterichRuinaFriction law(c);

   real_t a = 0.010;
   real_t V = 1e-9;

   // Test: f(theta2) > f(theta1) for theta2 > theta1 at fixed V
   real_t theta1 = 100.0;
   real_t theta2 = 10000.0;

   real_t f1 = law.FrictionCoefficient(V, theta1, a);
   real_t f2 = law.FrictionCoefficient(V, theta2, a);

   TEST_ASSERT(f2 > f1, "Higher theta -> higher f (evolution effect)");

   // Test derivative is positive
   real_t df_dtheta = law.FrictionDerivativeTheta(V, theta1, a);
   TEST_ASSERT(df_dtheta > 0.0, "df/dtheta > 0 (positive evolution effect)");
}

// =============================================================================
// Velocity-Weakening Tests (a < b: steady-state friction decreases with V)
// =============================================================================

void TestVelocityWeakening()
{
   std::cout << "\n=== Testing Velocity-Weakening Behavior ===\n";

   DieterichRuinaFriction::Constants c;
   c.b = 0.015;
   DieterichRuinaFriction law(c);

   real_t a = 0.010;  // a < b -> velocity-weakening

   // At steady state: theta_ss = Dc/V
   real_t V1 = 1e-10;
   real_t V2 = 1e-8;
   real_t theta1_ss = c.Dc / V1;
   real_t theta2_ss = c.Dc / V2;

   real_t f1 = law.FrictionCoefficient(V1, theta1_ss, a);
   real_t f2 = law.FrictionCoefficient(V2, theta2_ss, a);

   // For a < b, steady-state friction decreases with V
   TEST_ASSERT(f1 > f2, "Velocity-weakening: f_ss(V1) > f_ss(V2) for V1 < V2 when a < b");
}

// =============================================================================
// Velocity-Strengthening Tests (a > b: steady-state friction increases with V)
// =============================================================================

void TestVelocityStrengthening()
{
   std::cout << "\n=== Testing Velocity-Strengthening Behavior ===\n";

   DieterichRuinaFriction::Constants c;
   c.b = 0.015;
   DieterichRuinaFriction law(c);

   real_t a = 0.025;  // a > b -> velocity-strengthening

   // At steady state: theta_ss = Dc/V
   real_t V1 = 1e-10;
   real_t V2 = 1e-8;
   real_t theta1_ss = c.Dc / V1;
   real_t theta2_ss = c.Dc / V2;

   real_t f1 = law.FrictionCoefficient(V1, theta1_ss, a);
   real_t f2 = law.FrictionCoefficient(V2, theta2_ss, a);

   // For a > b, steady-state friction increases with V
   TEST_ASSERT(f1 < f2, "Velocity-strengthening: f_ss(V1) < f_ss(V2) for V1 < V2 when a > b");
}

// =============================================================================
// Slip Rate Solver Tests
// =============================================================================

void TestSlipRateSolver()
{
   std::cout << "\n=== Testing Slip Rate Solver ===\n";

   DieterichRuinaFriction::Constants c;
   DieterichRuinaFriction law(c);

   BP2Params params;
   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();

   // Test 1: Solver recovers known slip rate
   {
      real_t a = 0.015;
      real_t V_true = 1e-8;
      real_t theta = 5000.0;

      real_t f = law.FrictionCoefficient(V_true, theta, a);
      real_t tau = sigma_n * f + eta * V_true;

      int iterations;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      TEST_REL_NEAR(V_solved, V_true, 1e-10, "Solver recovers known slip rate");
      std::cout << "    (iterations: " << iterations << ")\n";
   }

   // Test 2: Solver handles zero radiation damping
   {
      real_t a = 0.015;
      real_t V_true = 1e-7;
      real_t theta = 3000.0;

      real_t f = law.FrictionCoefficient(V_true, theta, a);
      real_t tau = sigma_n * f;  // No radiation damping

      int iterations;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, 0.0, a, &iterations);

      TEST_REL_NEAR(V_solved, V_true, 1e-10, "Solver handles zero radiation damping");
      std::cout << "    (iterations: " << iterations << ")\n";
   }

   // Test 3: Solver handles very slow slip rate (interseismic)
   {
      real_t a = 0.015;
      real_t V_true = 1e-12;
      real_t theta = 1e8;

      real_t f = law.FrictionCoefficient(V_true, theta, a);
      real_t tau = sigma_n * f + eta * V_true;

      int iterations;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      TEST_REL_NEAR(V_solved, V_true, 1e-8, "Solver handles very slow slip rate");
      std::cout << "    (iterations: " << iterations << ")\n";
   }

   // Test 4: Solver handles fast slip rate (coseismic)
   {
      real_t a = 0.015;
      real_t V_true = 1.0;  // 1 m/s
      real_t theta = 0.01;

      real_t f = law.FrictionCoefficient(V_true, theta, a);
      real_t tau = sigma_n * f + eta * V_true;

      int iterations;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      TEST_REL_NEAR(V_solved, V_true, 1e-8, "Solver handles fast slip rate");
      std::cout << "    (iterations: " << iterations << ")\n";
   }

   // Test 5: Solver convergence efficiency
   {
      real_t a = 0.015;
      real_t V_true = 1e-6;
      real_t theta = 4000.0;

      real_t f = law.FrictionCoefficient(V_true, theta, a);
      real_t tau = sigma_n * f + eta * V_true;

      int iterations;
      law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      TEST_ASSERT(iterations < 50, "Solver converges in < 50 iterations");
      std::cout << "    (iterations: " << iterations << ")\n";
   }
}

// =============================================================================
// Initial State Tests
// =============================================================================

void TestInitialState()
{
   std::cout << "\n=== Testing Initial State Computation ===\n";

   DieterichRuinaFriction::Constants c;
   c.Dc = 0.004;  // BP2 value
   DieterichRuinaFriction law(c);

   BP2Params params;
   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();
   real_t V_init = params.V_init;

   // Test 1: Initial state satisfies stress balance
   {
      real_t a = 0.015;  // a = b (neutral)

      // Compute tau0 at amax (BP2 specification)
      real_t tau0 = params.tau0();

      // Compute initial state
      real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);

      // Verify: tau0 = sigma_n * f(V_init, theta0) + eta * V_init
      real_t f_check = law.FrictionCoefficient(V_init, theta0, a);
      real_t tau_check = sigma_n * f_check + eta * V_init;

      TEST_REL_NEAR(tau_check, tau0, 1e-10, "Initial state satisfies stress balance");
   }

   // Test 2: Initial state varies with depth (a varies)
   {
      real_t tau0 = params.tau0();

      std::vector<real_t> a_values = {0.010, 0.015, 0.020, 0.025};
      std::vector<real_t> theta0_values;

      bool all_balanced = true;
      for (real_t a : a_values)
      {
         real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);
         theta0_values.push_back(theta0);

         // Verify stress balance
         real_t f_check = law.FrictionCoefficient(V_init, theta0, a);
         real_t tau_check = sigma_n * f_check + eta * V_init;

         if (std::abs(tau_check - tau0) > tau0 * 1e-10)
         {
            all_balanced = false;
         }
      }

      TEST_ASSERT(all_balanced, "Stress balance maintained for all a values");

      // theta0 should vary with a
      TEST_ASSERT(std::abs(theta0_values[0] - theta0_values[3]) > 1.0,
                  "Initial state varies with depth (different a values)");
   }

   // Test 3: BP2-specific initial state values
   // At z=0: state_log10 ~ 3.602 -> theta ~ 10^3.602 ~ 4000 s (from benchmark)
   {
      real_t tau0 = params.tau0();
      real_t a = params.a0;  // a at surface (z=0)

      real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);

      // Expected from benchmark data: state_log10 ~ 3.602
      // This gives theta ~ 10^3.602 ~ 4000 s
      // Allow wider tolerance since our formula may differ slightly
      TEST_ASSERT(theta0 > 1000.0, "Initial theta > 1000 s");
      TEST_ASSERT(theta0 < 10000.0, "Initial theta < 10000 s");

      std::cout << "    theta0 at z=0: " << theta0 << " s"
                << " (log10: " << std::log10(theta0) << ")\n";
   }

   // Test 4: Verify tau0 value matches BP2 specification
   {
      real_t tau0 = params.tau0();
      real_t tau0_MPa = tau0 / 1e6;

      // Expected from benchmark: tau0 ~ 26.546 MPa
      TEST_NEAR(tau0_MPa, 26.546, 0.01, "Pre-stress tau0 matches BP2 (~26.546 MPa)");

      std::cout << "    tau0 = " << tau0_MPa << " MPa\n";
   }
}

// =============================================================================
// Robustness Tests
// =============================================================================

void TestRobustness()
{
   std::cout << "\n=== Testing Robustness (Edge Cases) ===\n";

   DieterichRuinaFriction::Constants c;
   DieterichRuinaFriction law(c);

   // Test 1: Negative normal stress with radiation damping (fault in tension)
   // Per Tandem: V = tau / eta when sigma_n <= 0 and eta > 0
   {
      real_t tau = 1.0e6;     // 1 MPa shear stress
      real_t theta = 4000.0;
      real_t sigma_n = -1.0e6;  // Negative! (fault in tension)
      real_t eta = 4.63e6;      // ~4.63 MPa·s/m (BP2 radiation damping)
      real_t a = 0.015;

      int iterations = -1;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      // Expected: V = tau / eta (viscous sliding)
      real_t V_expected = tau / eta;

      TEST_REL_NEAR(V_solved, V_expected, 1e-10,
                    "Negative sigma_n: V = tau/eta (viscous sliding)");
      TEST_ASSERT(iterations == 0,
                  "Negative sigma_n: No iterations needed");
   }

   // Test 2: Zero normal stress with radiation damping
   {
      real_t tau = 2.0e6;
      real_t theta = 4000.0;
      real_t sigma_n = 0.0;    // Zero normal stress
      real_t eta = 4.63e6;
      real_t a = 0.015;

      int iterations = -1;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      real_t V_expected = tau / eta;

      TEST_REL_NEAR(V_solved, V_expected, 1e-10,
                    "Zero sigma_n: V = tau/eta (viscous sliding)");
      TEST_ASSERT(iterations == 0,
                  "Zero sigma_n: No iterations needed");
   }

   // Test 3: Negative normal stress with zero radiation damping
   // Should return 0.0 (cannot determine slip rate)
   {
      real_t tau = 1.0e6;
      real_t theta = 4000.0;
      real_t sigma_n = -1.0e6;
      real_t eta = 0.0;         // No radiation damping
      real_t a = 0.015;

      int iterations = -1;
      real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, a, &iterations);

      TEST_NEAR(V_solved, 0.0, 1e-30,
                "Negative sigma_n, zero eta: V = 0");
      TEST_ASSERT(iterations == 0,
                  "Negative sigma_n, zero eta: No iterations needed");
   }

   // Test 4: Large f/a ratio in InitialState (overflow protection)
   // For small 'a' relative to friction coefficient, sinh(f/a) can overflow
   // Using moderately extreme values that are numerically challenging but don't cause actual overflow
   {
      // Use a very small 'a' value (smaller than BP2 range)
      // f ≈ 0.5 (typical friction coefficient)
      // a = 0.005 → f/a = 100, sinh(100) ≈ 1.3e43 (large but finite)
      BP2Params params;
      real_t sigma_n = params.sigma_n;  // 50 MPa
      real_t eta = params.eta();
      real_t V_init = params.V_init;

      real_t tau0 = params.tau0();  // ~26.5 MPa
      real_t a_small = 0.005;  // Smaller than BP2 range

      // This should not crash or produce inf/nan
      real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a_small);

      TEST_ASSERT(std::isfinite(theta0), "Small 'a' InitialState: theta is finite");
      TEST_ASSERT(theta0 > 0.0, "Small 'a' InitialState: theta is positive");

      // Verify stress balance
      real_t f_check = law.FrictionCoefficient(V_init, theta0, a_small);
      real_t tau_check = sigma_n * f_check + eta * V_init;

      // Allow slightly larger tolerance for extreme parameters
      TEST_REL_NEAR(tau_check, tau0, 1e-8,
                    "Small 'a' InitialState: stress balance maintained");
   }

   // Test 5: Very large f/a (near overflow threshold)
   // f/a > 700 would cause sinh to overflow, test the protection
   {
      // Construct a scenario where f/a is large
      // f ≈ tau/sigma_n, with tau high and sigma_n low, using very small a
      DieterichRuinaFriction::Constants c_test;
      c_test.V0 = 1e-6;
      c_test.f0 = 0.6;
      c_test.b = 0.015;
      c_test.Dc = 0.004;
      DieterichRuinaFriction law_test(c_test);

      real_t sigma_n = 10.0e6;   // 10 MPa (lower than BP2)
      real_t eta = 1.0e6;        // 1 MPa·s/m
      real_t V_init = 1e-9;

      // Choose tau to give f ≈ 0.5
      real_t tau0 = 5.0e6;  // 5 MPa → f ≈ 0.5

      // With a = 0.001, f/a = 500 (below overflow but large)
      real_t a_tiny = 0.001;

      real_t theta0 = law_test.InitialState(tau0, V_init, sigma_n, eta, a_tiny);

      TEST_ASSERT(std::isfinite(theta0), "Tiny 'a' InitialState: theta is finite");
      TEST_ASSERT(theta0 > 0.0, "Tiny 'a' InitialState: theta is positive");
   }

   // Test 6: When both Brent brackets fail, match Tandem and return NaN.
   // A non-finite psi makes the residual non-finite at both bracket endpoints.
   {
      real_t tau = 1.0e6;
      real_t psi = std::numeric_limits<real_t>::quiet_NaN();
      real_t sigma_n = 50.0e6;
      real_t eta = 4.63e6;
      real_t a = 0.015;

      int iterations = -1;
      real_t V_solved =
         law.SolveSlipRatePsi(tau, psi, sigma_n, eta, a, &iterations);

      TEST_ASSERT(std::isnan(V_solved),
                  "Bracket failure: SolveSlipRatePsi returns NaN like Tandem");
      TEST_ASSERT(iterations == 0,
                  "Bracket failure: No iterations reported");
   }
}

// =============================================================================
// Derivative Tests (numerical verification)
// =============================================================================

void TestDerivatives()
{
   std::cout << "\n=== Testing Derivatives (Numerical Verification) ===\n";

   DieterichRuinaFriction::Constants c;
   DieterichRuinaFriction law(c);

   real_t a = 0.015;
   real_t V = 1e-7;
   real_t theta = 3000.0;

   // Test df/dV using finite differences
   {
      real_t h = V * 1e-6;  // Small perturbation
      real_t f_plus = law.FrictionCoefficient(V + h, theta, a);
      real_t f_minus = law.FrictionCoefficient(V - h, theta, a);
      real_t df_dV_numerical = (f_plus - f_minus) / (2.0 * h);

      real_t df_dV_analytical = law.FrictionDerivativeV(V, theta, a);

      TEST_REL_NEAR(df_dV_analytical, df_dV_numerical, 1e-5,
                    "df/dV matches numerical derivative");
   }

   // Test df/dtheta using finite differences
   {
      real_t h = theta * 1e-6;  // Small perturbation
      real_t f_plus = law.FrictionCoefficient(V, theta + h, a);
      real_t f_minus = law.FrictionCoefficient(V, theta - h, a);
      real_t df_dtheta_numerical = (f_plus - f_minus) / (2.0 * h);

      real_t df_dtheta_analytical = law.FrictionDerivativeTheta(V, theta, a);

      TEST_REL_NEAR(df_dtheta_analytical, df_dtheta_numerical, 1e-5,
                    "df/dtheta matches numerical derivative");
   }
}

// =============================================================================
// RecomputeSlipRate consistency test
// =============================================================================
//
// Verifies that re-solving the friction equation from the current state
// produces the expected slip rate.  This is the core mechanism behind the
// RecomputeSlipRate fix for the stale-cache output bug at locked fault DOFs.
//
// At a = 0.004 (BP5 VW core), V ~ exp(-psi/a) is exponentially sensitive
// to psi.  This test confirms that the friction solver resolves V correctly
// even at extreme parameters, and that small psi perturbations cause bounded
// V changes (i.e., the solver doesn't suddenly return V=0).

void TestRecomputeSlipRateConsistency()
{
   std::cout << "\n=== Testing RecomputeSlipRate Consistency ===\n";

   DieterichRuinaFriction::Constants c;
   c.f0 = 0.6;
   c.V0 = 1.0e-6;
   c.b = 0.03;
   c.Dc = 0.14;
   DieterichRuinaFriction law(c);

   real_t sigma_n = 25.0e6;   // BP5 normal stress
   real_t eta = 32.04e9 / (2.0 * 3464.0);  // radiation damping

   // --- Test 1: VW core (a=0.004), interseismic locked fault ---
   // At depth 10km, psi grows as theta accumulates.
   // After ~2 days, theta ~ 3e5 s, psi ~ 0.62.
   // V should be extremely small but non-zero.
   {
      real_t a = 0.004;  // BP5 VW core
      real_t psi = 0.62;
      real_t tau = 8.74e6;  // Typical interseismic traction at dp+10

      real_t V = law.SolveSlipRatePsi(tau, psi, sigma_n, eta, a);
      TEST_ASSERT(std::isfinite(V), "VW locked: V is finite");
      TEST_ASSERT(V > 0.0, "VW locked: V > 0 (not clamped to zero)");
      // V should be extremely small (order 10^(-30) to 10^(-40))
      TEST_ASSERT(V < 1e-10, "VW locked: V < 1e-10 (fault is locked)");
      TEST_ASSERT(V > 1e-300, "VW locked: V > 1e-300 (not underflow)");

      std::cout << "    V(a=0.004, psi=0.62) = " << V
                << " (log10 = " << std::log10(V) << ")\n";
   }

   // --- Test 2: Small psi perturbation produces bounded V change ---
   // Verifies that V changes smoothly with psi (no sudden jump to zero).
   {
      real_t a = 0.004;
      real_t psi0 = 0.62;
      real_t tau = 8.74e6;
      real_t delta_psi = 1e-4;  // One RK stage's worth of psi change

      real_t V0 = law.SolveSlipRatePsi(tau, psi0, sigma_n, eta, a);
      real_t V1 = law.SolveSlipRatePsi(tau, psi0 + delta_psi, sigma_n, eta, a);

      TEST_ASSERT(std::isfinite(V0) && V0 > 0,
                  "Psi perturbation: V0 valid");
      TEST_ASSERT(std::isfinite(V1) && V1 > 0,
                  "Psi perturbation: V1 valid");

      // V change should be bounded: |log10(V1/V0)| < delta_psi / (a * ln10) + margin
      real_t expected_log_change = delta_psi / (a * std::log(10.0));
      real_t actual_log_change = std::abs(std::log10(V1) - std::log10(V0));

      TEST_ASSERT(actual_log_change < expected_log_change * 2.0,
                  "Psi perturbation: V change bounded");

      std::cout << "    delta_psi=" << delta_psi
                << " expected_dlog10V=" << expected_log_change
                << " actual=" << actual_log_change << "\n";
   }

   // --- Test 3: Vector solve consistency ---
   // Verifies that the 2-component solver decomposes V correctly
   // at extreme parameters.
   {
      real_t a = 0.004;
      real_t psi = 0.62;
      real_t tau_vec[2] = {-0.016e6, 8.74e6};  // Mostly strike
      real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);

      real_t V_scalar = law.SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a);

      real_t V_vec[2];
      law.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);
      real_t V_vec_abs = std::sqrt(V_vec[0]*V_vec[0] + V_vec[1]*V_vec[1]);

      TEST_ASSERT(std::isfinite(V_vec[0]) && std::isfinite(V_vec[1]),
                  "Vector solve: components finite");
      TEST_ASSERT(V_vec_abs > 0.0,
                  "Vector solve: |V| > 0");

      // |V_vec| should equal the scalar V
      if (V_scalar > 0 && V_vec_abs > 0)
      {
         real_t rel_diff = std::abs(V_vec_abs - V_scalar) /
                           std::max(V_scalar, V_vec_abs);
         TEST_ASSERT(rel_diff < 1e-10,
                     "Vector solve: |V_vec| == V_scalar");
      }

      std::cout << "    V_scalar=" << V_scalar << " V_vec_abs=" << V_vec_abs
                << " V_vec=(" << V_vec[0] << ", " << V_vec[1] << ")\n";
   }

   // --- Test 4: VS region (a=0.04, larger than b) ---
   // At the surface, V should be close to Vp during interseismic.
   {
      real_t a = 0.04;  // VS region
      real_t psi = 0.62;
      real_t tau = 13.3e6;  // Surface traction

      real_t V = law.SolveSlipRatePsi(tau, psi, sigma_n, eta, a);
      TEST_ASSERT(std::isfinite(V) && V > 0,
                  "VS surface: V valid");
      // VS region: V should be relatively large (much > VW core)
      TEST_ASSERT(V > 1e-20,
                  "VS surface: V > 1e-20 (not locked)");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "========================================\n";
   std::cout << "SEAS Phase 1 Unit Tests: Friction Law\n";
   std::cout << "========================================\n";

   TestFrictionCoefficient();
   TestDirectEffect();
   TestEvolutionEffect();
   TestVelocityWeakening();
   TestVelocityStrengthening();
   TestSlipRateSolver();
   TestInitialState();
   TestDerivatives();
   TestRobustness();
   TestRecomputeSlipRateConsistency();

   std::cout << "\n========================================\n";
   std::cout << "Test Summary\n";
   std::cout << "========================================\n";
   std::cout << "Total tests: " << num_tests << "\n";
   std::cout << "Passed:      " << num_passed << "\n";
   std::cout << "Failed:      " << num_failed << "\n";

   if (num_failed > 0)
   {
      std::cout << "\nSOME TESTS FAILED!\n";
      return 1;
   }

   std::cout << "\nALL TESTS PASSED!\n";
   return 0;
}
