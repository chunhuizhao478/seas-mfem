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

// Unit tests for psi-space state variable integration.
// Tests that psi-space (logarithmic) and theta-space (physical) formulations
// produce identical results, following Tandem's approach.

#include "mfem.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Simple test framework (matches existing test patterns)
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
// Test 1: Psi-Theta round-trip conversion
// =============================================================================

void TestPsiThetaConversion()
{
   std::cout << "\n=== Test 1: Psi-Theta Conversion ===\n";

   DieterichRuinaFriction::Constants c;
   c.V0 = 1e-6;
   c.f0 = 0.6;
   c.b = 0.015;
   c.Dc = 0.004;
   DieterichRuinaFriction law(c);

   // Test round-trip at multiple theta values spanning the physical range
   std::vector<real_t> theta_values = {0.004, 1.0, 1000.0, 1e6, 4e6};

   for (real_t theta_orig : theta_values)
   {
      real_t psi = law.ThetaToPsi(theta_orig);
      real_t theta_back = law.PsiToTheta(psi);
      TEST_REL_NEAR(theta_back, theta_orig, 1e-12,
                    "theta->psi->theta round-trip");
   }

   // Test round-trip starting from psi values
   std::vector<real_t> psi_values = {0.3, 0.5, 0.6, 0.7, 0.9};

   for (real_t psi_orig : psi_values)
   {
      real_t theta = law.PsiToTheta(psi_orig);
      real_t psi_back = law.ThetaToPsi(theta);
      TEST_REL_NEAR(psi_back, psi_orig, 1e-12,
                    "psi->theta->psi round-trip");
   }

   // Verify formula: psi = f0 + b*ln(V0*theta/Dc)
   {
      real_t theta = 4000.0;
      real_t psi = law.ThetaToPsi(theta);
      real_t psi_expected = c.f0 + c.b * std::log(c.V0 * theta / c.Dc);
      TEST_NEAR(psi, psi_expected, 1e-15, "ThetaToPsi matches formula");
   }

   // Verify formula: theta = (Dc/V0)*exp((psi-f0)/b)
   {
      real_t psi = 0.6;
      real_t theta = law.PsiToTheta(psi);
      real_t theta_expected = (c.Dc / c.V0) * std::exp((psi - c.f0) / c.b);
      TEST_NEAR(theta, theta_expected, 1e-15, "PsiToTheta matches formula");
   }
}

// =============================================================================
// Test 2: AgingLawPsi steady state
// =============================================================================

void TestAgingLawPsiSteadyState()
{
   std::cout << "\n=== Test 2: AgingLawPsi Steady State ===\n";

   real_t b = 0.015;
   real_t V0 = 1e-6;
   real_t f0 = 0.6;
   real_t Dc = 0.004;

   AgingLawPsi aging_psi(b, V0, f0);

   // At steady state, dpsi/dt = 0 for psi_ss = f0 + b*ln(V0/V)
   std::vector<real_t> V_values = {1e-12, 1e-9, 1e-6, 1e-3, 1.0};

   for (real_t V : V_values)
   {
      real_t psi_ss = aging_psi.SteadyState(V, Dc);
      real_t rate = aging_psi.Rate(V, psi_ss, Dc);
      // Scale tolerance by the magnitude of the largest term in the rate
      real_t scale = std::max(b * V0 / Dc, b * V / Dc);
      TEST_NEAR(rate, 0.0, 1e-12 * scale,
                "dpsi/dt = 0 at steady state");
   }

   // Verify psi_ss = f0 + b*ln(V0/V)
   {
      real_t V = 1e-9;
      real_t psi_ss = aging_psi.SteadyState(V, Dc);
      real_t psi_expected = f0 + b * std::log(V0 / V);
      TEST_NEAR(psi_ss, psi_expected, 1e-15, "psi_ss matches formula");
   }

   // Verify steady-state psi matches steady-state theta converted
   {
      DieterichRuinaFriction::Constants c{V0, f0, b, Dc};
      DieterichRuinaFriction law(c);
      AgingLaw aging_theta;

      real_t V = 1e-9;
      real_t theta_ss = aging_theta.SteadyState(V, Dc);
      real_t psi_from_theta = law.ThetaToPsi(theta_ss);
      real_t psi_ss = aging_psi.SteadyState(V, Dc);

      TEST_REL_NEAR(psi_ss, psi_from_theta, 1e-12,
                    "psi_ss = ThetaToPsi(theta_ss)");
   }
}

// =============================================================================
// Test 3: AgingLawPsi rate matches theta-space via chain rule
// =============================================================================

void TestAgingLawPsiRate()
{
   std::cout << "\n=== Test 3: AgingLawPsi Rate vs Theta-Space ===\n";

   real_t b = 0.015;
   real_t V0 = 1e-6;
   real_t f0 = 0.6;
   real_t Dc = 0.004;

   AgingLawPsi aging_psi(b, V0, f0);
   AgingLaw aging_theta;

   DieterichRuinaFriction::Constants c{V0, f0, b, Dc};
   DieterichRuinaFriction law(c);

   // dpsi/dt = (b/theta) * dtheta/dt
   // because psi = f0 + b*ln(V0*theta/Dc), so dpsi/dtheta = b/theta
   struct TestCase { real_t V; real_t theta; };
   std::vector<TestCase> cases = {
      {1e-12, 1e6},
      {1e-9, 4000.0},
      {1e-6, 4.0},
      {1e-3, 0.004},
      {1.0, 0.001}
   };

   for (auto &tc : cases)
   {
      real_t psi = law.ThetaToPsi(tc.theta);

      real_t dtheta_dt = aging_theta.Rate(tc.V, tc.theta, Dc);
      real_t dpsi_dt_from_theta = (b / tc.theta) * dtheta_dt;
      real_t dpsi_dt_direct = aging_psi.Rate(tc.V, psi, Dc);

      real_t scale = std::max(std::abs(dpsi_dt_from_theta), 1e-20);
      TEST_NEAR(dpsi_dt_direct, dpsi_dt_from_theta, 1e-10 * scale,
                "dpsi/dt = (b/theta)*dtheta/dt");
   }
}

// =============================================================================
// Test 4: Friction coefficient equivalence (psi vs theta)
// =============================================================================

void TestFrictionPsiEquivalence()
{
   std::cout << "\n=== Test 4: Friction Coefficient Psi Equivalence ===\n";

   DieterichRuinaFriction::Constants c;
   c.V0 = 1e-6;
   c.f0 = 0.6;
   c.b = 0.015;
   c.Dc = 0.004;
   DieterichRuinaFriction law(c);

   std::vector<real_t> a_values = {0.010, 0.015, 0.020, 0.025};
   std::vector<real_t> V_values = {1e-12, 1e-9, 1e-6, 1e-3};
   std::vector<real_t> theta_values = {0.01, 100.0, 4000.0, 1e6};

   for (real_t a : a_values)
   {
      for (real_t V : V_values)
      {
         for (real_t theta : theta_values)
         {
            real_t psi = law.ThetaToPsi(theta);

            real_t f_theta = law.FrictionCoefficient(V, theta, a);
            real_t f_psi = law.FrictionCoefficientPsi(V, psi, a);

            TEST_REL_NEAR(f_psi, f_theta, 1e-12,
                          "f(V,psi) == f(V,theta) after conversion");
         }
      }
   }
}

// =============================================================================
// Test 5: SolveSlipRatePsi matches SolveSlipRate
// =============================================================================

void TestSolveSlipRatePsi()
{
   std::cout << "\n=== Test 5: SolveSlipRatePsi Equivalence ===\n";

   DieterichRuinaFriction::Constants c;
   c.V0 = 1e-6;
   c.f0 = 0.6;
   c.b = 0.015;
   c.Dc = 0.004;
   DieterichRuinaFriction law(c);

   BP2Params params;
   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();

   struct TestCase { real_t V_true; real_t theta; real_t a; };
   std::vector<TestCase> cases = {
      {1e-12, 1e8, 0.015},     // interseismic
      {1e-9, 4000.0, 0.010},   // velocity-weakening
      {1e-6, 4.0, 0.015},      // reference
      {1e-3, 0.01, 0.020},     // coseismic
      {1.0, 0.001, 0.025},     // fast slip
   };

   for (auto &tc : cases)
   {
      real_t psi = law.ThetaToPsi(tc.theta);

      // Compute tau from known V, theta
      real_t f = law.FrictionCoefficient(tc.V_true, tc.theta, tc.a);
      real_t tau = sigma_n * f + eta * tc.V_true;

      // Solve in theta-space
      real_t V_theta = law.SolveSlipRate(tau, tc.theta, sigma_n, eta, tc.a);

      // Solve in psi-space
      real_t V_psi = law.SolveSlipRatePsi(tau, psi, sigma_n, eta, tc.a);

      TEST_REL_NEAR(V_psi, V_theta, 1e-8,
                    "SolveSlipRatePsi == SolveSlipRate");
   }
}

// =============================================================================
// Test 6: InitialStatePsi matches InitialState after conversion
// =============================================================================

void TestInitialStatePsi()
{
   std::cout << "\n=== Test 6: InitialStatePsi Equivalence ===\n";

   DieterichRuinaFriction::Constants c;
   c.V0 = 1e-6;
   c.f0 = 0.6;
   c.b = 0.015;
   c.Dc = 0.004;
   DieterichRuinaFriction law(c);

   BP2Params params;
   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();
   real_t V_init = params.V_init;
   real_t tau0 = params.tau0();

   std::vector<real_t> a_values = {0.010, 0.015, 0.020, 0.025};

   for (real_t a : a_values)
   {
      real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, a);
      real_t psi0 = law.InitialStatePsi(tau0, V_init, sigma_n, eta, a);

      // Convert psi0 to theta and compare
      real_t theta_from_psi = law.PsiToTheta(psi0);
      TEST_REL_NEAR(theta_from_psi, theta0, 1e-10,
                    "PsiToTheta(InitialStatePsi) == InitialState");

      // Verify stress balance with psi0
      real_t f_psi = law.FrictionCoefficientPsi(V_init, psi0, a);
      real_t tau_check = sigma_n * f_psi + eta * V_init;
      TEST_REL_NEAR(tau_check, tau0, 1e-10,
                    "InitialStatePsi satisfies stress balance");
   }
}

// =============================================================================
// Test 7: Forward Euler time integration in both spaces
// =============================================================================

void TestTimeIntegration()
{
   std::cout << "\n=== Test 7: Time Integration Equivalence ===\n";

   real_t b = 0.015;
   real_t V0 = 1e-6;
   real_t f0 = 0.6;
   real_t Dc = 0.004;

   DieterichRuinaFriction::Constants c{V0, f0, b, Dc};
   DieterichRuinaFriction law(c);

   AgingLaw aging_theta;
   AgingLawPsi aging_psi(b, V0, f0);

   BP2Params params;
   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();
   real_t a = 0.015;

   // Initial conditions: stress equilibrium
   real_t tau0 = params.tau0();
   real_t theta = law.InitialState(tau0, params.V_init, sigma_n, eta, a);
   real_t psi = law.ThetaToPsi(theta);
   real_t slip_theta = 0.0;
   real_t slip_psi = 0.0;

   // Apply a small perturbation to tau to drive evolution
   real_t tau = tau0 + 1.0e3;  // +1 kPa perturbation

   // Forward Euler with small time step
   int N = 1000;
   real_t dt = 1.0;  // 1 second steps

   for (int step = 0; step < N; step++)
   {
      // Theta-space
      real_t V_theta = law.SolveSlipRate(tau, theta, sigma_n, eta, a);
      real_t dtheta_dt = aging_theta.Rate(V_theta, theta, Dc);
      slip_theta += V_theta * dt;
      theta += dtheta_dt * dt;

      // Psi-space
      real_t V_psi = law.SolveSlipRatePsi(tau, psi, sigma_n, eta, a);
      real_t dpsi_dt = aging_psi.Rate(V_psi, psi, Dc);
      slip_psi += V_psi * dt;
      psi += dpsi_dt * dt;
   }

   // Compare final results
   real_t theta_from_psi = law.PsiToTheta(psi);

   // Forward Euler accumulates small differences between parameterizations
   TEST_REL_NEAR(slip_psi, slip_theta, 1e-6,
                 "Slip matches after time integration");
   TEST_REL_NEAR(theta_from_psi, theta, 1e-6,
                 "Theta matches after time integration");

   std::cout << "    Final slip (theta-space): " << slip_theta << " m\n";
   std::cout << "    Final slip (psi-space):   " << slip_psi << " m\n";
   std::cout << "    Final theta (theta-space): " << theta << " s\n";
   std::cout << "    Final theta (psi-space):   " << theta_from_psi << " s\n";
}

// =============================================================================
// Test 8: SlipLawPsi equivalence
// =============================================================================

void TestSlipLawPsi()
{
   std::cout << "\n=== Test 8: SlipLawPsi Equivalence ===\n";

   real_t b = 0.015;
   real_t V0 = 1e-6;
   real_t f0 = 0.6;
   real_t Dc = 0.004;

   SlipLaw slip_theta;
   SlipLawPsi slip_psi(b, V0, f0);

   DieterichRuinaFriction::Constants c{V0, f0, b, Dc};
   DieterichRuinaFriction law(c);

   // Test steady state matches
   {
      real_t V = 1e-9;
      real_t psi_ss = slip_psi.SteadyState(V, Dc);
      real_t theta_ss = slip_theta.SteadyState(V, Dc);
      real_t psi_from_theta = law.ThetaToPsi(theta_ss);

      TEST_REL_NEAR(psi_ss, psi_from_theta, 1e-12,
                    "SlipLawPsi steady state matches");
   }

   // Test rate via chain rule: dpsi/dt = (b/theta) * dtheta/dt
   struct TestCase { real_t V; real_t theta; };
   std::vector<TestCase> cases = {
      {1e-9, 4000.0},
      {1e-6, 4.0},
      {1e-3, 0.01}
   };

   for (auto &tc : cases)
   {
      real_t psi = law.ThetaToPsi(tc.theta);

      real_t dtheta_dt = slip_theta.Rate(tc.V, tc.theta, Dc);
      real_t dpsi_dt_from_theta = (b / tc.theta) * dtheta_dt;
      real_t dpsi_dt_direct = slip_psi.Rate(tc.V, psi, Dc);

      real_t scale = std::max(std::abs(dpsi_dt_from_theta), 1e-20);
      TEST_NEAR(dpsi_dt_direct, dpsi_dt_from_theta, 1e-10 * scale,
                "SlipLawPsi rate matches chain rule");
   }

   // Test dpsi/dt = 0 at steady state
   {
      real_t V = 1e-9;
      real_t psi_ss = slip_psi.SteadyState(V, Dc);
      real_t rate = slip_psi.Rate(V, psi_ss, Dc);
      TEST_NEAR(rate, 0.0, 1e-10 * b * V / Dc,
                "SlipLawPsi dpsi/dt = 0 at steady state");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "========================================\n";
   std::cout << "SEAS Unit Tests: Psi-Space State Variable\n";
   std::cout << "========================================\n";

   TestPsiThetaConversion();
   TestAgingLawPsiSteadyState();
   TestAgingLawPsiRate();
   TestFrictionPsiEquivalence();
   TestSolveSlipRatePsi();
   TestInitialStatePsi();
   TestTimeIntegration();
   TestSlipLawPsi();

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
