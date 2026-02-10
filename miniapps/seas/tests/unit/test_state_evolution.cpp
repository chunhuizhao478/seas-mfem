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

// Unit tests for Phase 1: State Evolution Laws
// Tests for AgingLaw and SlipLaw classes:
// - Steady state computation
// - Rate computation at steady state (should be zero)
// - Rate sign when above/below steady state
// - Derivatives
// - Timescale characteristics

#include "mfem.hpp"
#include "../../friction/state_evolution.hpp"
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
      if (_tol < 1e-30) _tol = 1e-30; \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", RelErr: " << std::abs(_val - _exp) / (std::abs(_exp) + 1e-30) << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// Aging Law Tests
// =============================================================================

void TestAgingLaw()
{
   std::cout << "\n=== Testing Aging Law ===\n";

   AgingLaw aging;
   BP2Params params;
   real_t Dc = params.Dc;

   // Test 1: Steady state computation
   {
      real_t V = 1e-9;
      real_t theta_ss = aging.SteadyState(V, Dc);
      real_t expected = Dc / V;

      TEST_NEAR(theta_ss, expected, 1e-12 * expected, "AgingLaw: SteadyState = Dc/V");
   }

   // Test 2: Rate at steady state should be zero
   {
      real_t V = 1e-9;
      real_t theta_ss = aging.SteadyState(V, Dc);
      real_t rate = aging.Rate(V, theta_ss, Dc);

      TEST_NEAR(rate, 0.0, 1e-20, "AgingLaw: Rate = 0 at steady state");
   }

   // Test 3: Rate is positive when V*theta < Dc (below steady state)
   // This means theta will increase (fault heals)
   {
      real_t V = 1e-9;
      real_t theta = 100.0;  // V*theta = 1e-7 << Dc = 0.004

      real_t rate = aging.Rate(V, theta, Dc);

      TEST_ASSERT(rate > 0.0, "AgingLaw: dtheta/dt > 0 when V*theta < Dc (healing)");
   }

   // Test 4: Rate is negative when V*theta > Dc (above steady state)
   // This means theta will decrease (contacts broken by slip)
   {
      real_t V = 1.0;  // 1 m/s (coseismic)
      real_t theta = 1000.0;  // V*theta = 1000 >> Dc = 0.004

      real_t rate = aging.Rate(V, theta, Dc);

      TEST_ASSERT(rate < 0.0, "AgingLaw: dtheta/dt < 0 when V*theta > Dc (weakening)");
   }

   // Test 5: Healing when fault is locked (V ~ 0)
   // dtheta/dt = 1 - V*theta/Dc ≈ 1 when V → 0
   {
      real_t V = 1e-20;  // Essentially locked
      real_t theta = 1000.0;

      real_t rate = aging.Rate(V, theta, Dc);

      TEST_NEAR(rate, 1.0, 1e-10, "AgingLaw: dtheta/dt ≈ 1 when fault is locked");
   }

   // Test 6: Derivative dG/dV = -theta/Dc
   {
      real_t V = 1e-6;
      real_t theta = 4000.0;

      real_t dG_dV = aging.RateDerivativeV(V, theta, Dc);
      real_t expected = -theta / Dc;

      TEST_NEAR(dG_dV, expected, 1e-12 * std::abs(expected),
                "AgingLaw: dG/dV = -theta/Dc");
   }

   // Test 7: Derivative dG/dtheta = -V/Dc
   {
      real_t V = 1e-6;
      real_t theta = 4000.0;

      real_t dG_dtheta = aging.RateDerivativeTheta(V, theta, Dc);
      real_t expected = -V / Dc;

      TEST_NEAR(dG_dtheta, expected, 1e-12 * std::abs(expected),
                "AgingLaw: dG/dtheta = -V/Dc");
   }

   // Test 8: Characteristic timescale
   // The characteristic evolution time is tau = Dc/V
   {
      real_t V = 1e-9;  // Plate rate
      real_t timescale = Dc / V;

      // For BP2: Dc = 0.004 m, V = 1e-9 m/s
      // timescale = 0.004 / 1e-9 = 4e6 s ≈ 46 days
      real_t expected_days = 46.3;  // 4e6 / (24*3600)
      real_t timescale_days = timescale / (24.0 * 3600.0);

      TEST_NEAR(timescale_days, expected_days, 1.0, "AgingLaw: Characteristic timescale ~46 days");
      std::cout << "    (timescale: " << timescale_days << " days)\n";
   }
}

// =============================================================================
// Slip Law Tests
// =============================================================================

void TestSlipLaw()
{
   std::cout << "\n=== Testing Slip Law ===\n";

   SlipLaw slip;
   BP2Params params;
   real_t Dc = params.Dc;

   // Test 1: Steady state computation (same as aging law)
   {
      real_t V = 1e-9;
      real_t theta_ss = slip.SteadyState(V, Dc);
      real_t expected = Dc / V;

      TEST_NEAR(theta_ss, expected, 1e-12 * expected, "SlipLaw: SteadyState = Dc/V");
   }

   // Test 2: Rate at steady state should be zero
   // At steady state: V*theta/Dc = 1, so ln(1) = 0, rate = 0
   {
      real_t V = 1e-9;
      real_t theta_ss = slip.SteadyState(V, Dc);
      real_t rate = slip.Rate(V, theta_ss, Dc);

      TEST_NEAR(rate, 0.0, 1e-15, "SlipLaw: Rate = 0 at steady state");
   }

   // Test 3: Rate is positive when V*theta < Dc (x < 1 means ln(x) < 0)
   // dtheta/dt = -x*ln(x) > 0 when x < 1
   {
      real_t V = 1e-9;
      real_t theta = 100.0;  // x = V*theta/Dc = 1e-9 * 100 / 0.004 = 2.5e-5 << 1

      real_t rate = slip.Rate(V, theta, Dc);

      TEST_ASSERT(rate > 0.0, "SlipLaw: dtheta/dt > 0 when V*theta < Dc");
   }

   // Test 4: Rate is negative when V*theta > Dc (x > 1 means ln(x) > 0)
   // dtheta/dt = -x*ln(x) < 0 when x > 1
   {
      real_t V = 1.0;  // 1 m/s
      real_t theta = 1.0;  // x = V*theta/Dc = 1.0 * 1.0 / 0.004 = 250 >> 1

      real_t rate = slip.Rate(V, theta, Dc);

      TEST_ASSERT(rate < 0.0, "SlipLaw: dtheta/dt < 0 when V*theta > Dc");
   }

   // Test 5: Slip law requires slip for evolution (no healing when locked)
   // When V → 0, dtheta/dt → 0 (unlike aging law)
   // For x → 0: -x*ln(x) → 0 (L'Hopital)
   {
      real_t V = 1e-15;  // Very slow
      real_t theta = 1000.0;

      real_t rate = slip.Rate(V, theta, Dc);

      // Should be small (approaches 0 as V → 0)
      TEST_ASSERT(std::abs(rate) < 1e-8, "SlipLaw: Rate → 0 when V → 0 (no healing without slip)");
   }

   // Test 6: Verify derivatives numerically
   {
      real_t V = 1e-7;
      real_t theta = 3000.0;
      real_t h = V * 1e-6;  // Use relative step size for better accuracy

      // Numerical derivative dG/dV
      real_t rate_plus = slip.Rate(V + h, theta, Dc);
      real_t rate_minus = slip.Rate(V - h, theta, Dc);
      real_t dG_dV_numerical = (rate_plus - rate_minus) / (2.0 * h);

      real_t dG_dV_analytical = slip.RateDerivativeV(V, theta, Dc);

      // Use 1% tolerance for numerical derivative comparison
      TEST_REL_NEAR(dG_dV_analytical, dG_dV_numerical, 1e-2,
                    "SlipLaw: dG/dV matches numerical derivative");
   }

   // Test 7: Verify derivatives numerically for dG/dtheta
   {
      real_t V = 1e-7;
      real_t theta = 3000.0;
      real_t h = 1.0;  // Small perturbation in theta

      real_t rate_plus = slip.Rate(V, theta + h, Dc);
      real_t rate_minus = slip.Rate(V, theta - h, Dc);
      real_t dG_dtheta_numerical = (rate_plus - rate_minus) / (2.0 * h);

      real_t dG_dtheta_analytical = slip.RateDerivativeTheta(V, theta, Dc);

      TEST_REL_NEAR(dG_dtheta_analytical, dG_dtheta_numerical, 1e-4,
                    "SlipLaw: dG/dtheta matches numerical derivative");
   }
}

// =============================================================================
// Comparison Tests (Aging vs Slip Law)
// =============================================================================

void TestLawComparison()
{
   std::cout << "\n=== Testing Aging vs Slip Law Comparison ===\n";

   AgingLaw aging;
   SlipLaw slip;
   BP2Params params;
   real_t Dc = params.Dc;

   // Test 1: Both laws have same steady state
   {
      real_t V = 1e-9;

      real_t theta_ss_aging = aging.SteadyState(V, Dc);
      real_t theta_ss_slip = slip.SteadyState(V, Dc);

      TEST_NEAR(theta_ss_aging, theta_ss_slip, 1e-12,
                "Both laws have same steady state");
   }

   // Test 2: Both laws give zero rate at steady state
   {
      real_t V = 1e-9;
      real_t theta_ss = Dc / V;

      real_t rate_aging = aging.Rate(V, theta_ss, Dc);
      real_t rate_slip = slip.Rate(V, theta_ss, Dc);

      TEST_NEAR(rate_aging, 0.0, 1e-15, "Aging law: zero rate at steady state");
      TEST_NEAR(rate_slip, 0.0, 1e-15, "Slip law: zero rate at steady state");
   }

   // Test 3: Aging law heals when locked, slip law does not
   {
      real_t V = 1e-20;  // Essentially locked
      real_t theta = 1000.0;

      real_t rate_aging = aging.Rate(V, theta, Dc);
      real_t rate_slip = slip.Rate(V, theta, Dc);

      // Aging law: dtheta/dt ≈ 1 (healing)
      TEST_NEAR(rate_aging, 1.0, 1e-10, "Aging law heals when locked (dtheta/dt ≈ 1)");

      // Slip law: dtheta/dt ≈ 0 (no healing without slip)
      TEST_ASSERT(std::abs(rate_slip) < 1e-10, "Slip law does not heal when locked");
   }

   // Test 4: Both laws weaken during rapid slip
   {
      real_t V = 1.0;  // Coseismic
      real_t theta = 1000.0;

      real_t rate_aging = aging.Rate(V, theta, Dc);
      real_t rate_slip = slip.Rate(V, theta, Dc);

      TEST_ASSERT(rate_aging < 0.0, "Aging law weakens during rapid slip");
      TEST_ASSERT(rate_slip < 0.0, "Slip law weakens during rapid slip");
   }
}

// =============================================================================
// Integration Test: State Evolution Over Time
// =============================================================================

void TestTimeIntegration()
{
   std::cout << "\n=== Testing State Evolution Over Time ===\n";

   AgingLaw aging;
   BP2Params params;
   real_t Dc = params.Dc;

   // Test: Forward Euler integration converges to steady state
   {
      real_t V = 1e-9;
      real_t theta = 100.0;  // Start below steady state
      real_t theta_ss = Dc / V;  // Steady state

      // Integrate for one characteristic time
      real_t dt = 1e4;  // 10000 seconds
      int n_steps = 1000;

      for (int i = 0; i < n_steps; ++i)
      {
         real_t rate = aging.Rate(V, theta, Dc);
         theta += dt * rate;
      }

      // Should be close to steady state after long time
      real_t total_time = dt * n_steps;  // 1e7 seconds
      real_t timescale = Dc / V;  // ~4e6 seconds

      std::cout << "    Initial theta: 100 s\n";
      std::cout << "    Final theta:   " << theta << " s\n";
      std::cout << "    Steady state:  " << theta_ss << " s\n";
      std::cout << "    Integration time: " << total_time / timescale << " timescales\n";

      // After 2.5 timescales, should be within 10% of steady state
      real_t rel_error = std::abs(theta - theta_ss) / theta_ss;
      TEST_ASSERT(rel_error < 0.1, "State converges to steady state (rel error < 10%)");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "========================================\n";
   std::cout << "SEAS Phase 1 Unit Tests: State Evolution\n";
   std::cout << "========================================\n";

   TestAgingLaw();
   TestSlipLaw();
   TestLawComparison();
   TestTimeIntegration();

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
