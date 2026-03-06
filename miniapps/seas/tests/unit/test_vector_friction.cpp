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

// Unit tests for vector slip rate solver (SolveSlipRateVectorPsi).
// Tests direction, magnitude consistency, and self-consistency with BP5 params.

#include "test_macros.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../config/bp5_params.hpp"

#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test: Direction of vector slip rate
// =============================================================================
void TestDirection()
{
   std::cout << "\n=== Vector Slip Rate Direction ===\n";

   DieterichRuinaFriction::Constants fc;
   fc.V0 = 1.0e-6;
   fc.f0 = 0.6;
   fc.b = 0.03;
   fc.Dc = 0.14;
   DieterichRuinaFriction friction(fc);

   real_t sigma_n = 25.0e6;
   real_t eta = 2670.0 * 3464.0 * 3464.0 / (2.0 * 3464.0); // mu/(2*cs)
   real_t a = 0.004;  // VW zone
   real_t psi = fc.f0 + fc.b * std::log(fc.V0 / 1e-9);

   // Test 1: Pure x-direction stress
   {
      real_t tau_vec[2] = {20.0e6, 0.0};
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);

      TEST_ASSERT(V_vec[0] < 0.0, "Pure x-stress: V[0] < 0 (anti-parallel)");
      TEST_NEAR(V_vec[1], 0.0, 1e-30, "Pure x-stress: V[1] = 0");
   }

   // Test 2: Pure y-direction stress
   {
      real_t tau_vec[2] = {0.0, 20.0e6};
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);

      TEST_NEAR(V_vec[0], 0.0, 1e-30, "Pure y-stress: V[0] = 0");
      TEST_ASSERT(V_vec[1] < 0.0, "Pure y-stress: V[1] < 0 (anti-parallel)");
   }

   // Test 3: 45-degree stress
   {
      real_t tau_val = 20.0e6;
      real_t tau_vec[2] = {tau_val, tau_val};
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);

      // Both components should be negative (anti-parallel to positive tau)
      TEST_ASSERT(V_vec[0] < 0.0, "45-deg: V[0] < 0");
      TEST_ASSERT(V_vec[1] < 0.0, "45-deg: V[1] < 0");

      // Components should be equal (same magnitude in both directions)
      TEST_REL_NEAR(V_vec[0], V_vec[1], 1e-12,
                     "45-deg: V[0] = V[1] (symmetric)");
   }

   // Test 4: Negative stress direction
   {
      real_t tau_vec[2] = {-15.0e6, 10.0e6};
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);

      // V should be anti-parallel to tau: V[0] > 0, V[1] < 0
      TEST_ASSERT(V_vec[0] > 0.0,
                  "(-,+) stress: V[0] > 0 (anti-parallel)");
      TEST_ASSERT(V_vec[1] < 0.0,
                  "(-,+) stress: V[1] < 0 (anti-parallel)");
   }
}

// =============================================================================
// Test: Magnitude consistency with scalar solver
// =============================================================================
void TestMagnitudeConsistency()
{
   std::cout << "\n=== Magnitude Consistency ===\n";

   DieterichRuinaFriction::Constants fc;
   fc.V0 = 1.0e-6;
   fc.f0 = 0.6;
   fc.b = 0.03;
   fc.Dc = 0.14;
   DieterichRuinaFriction friction(fc);

   real_t sigma_n = 25.0e6;
   real_t eta = 2670.0 * 3464.0 * 3464.0 / (2.0 * 3464.0);

   // Use a = 0.01 (moderate transition zone value) to ensure tau values
   // are above the minimum friction force. For a=0.004, the minimum
   // solvable tau is ~15 MPa; for a=0.01 it's ~6 MPa.
   real_t a = 0.01;
   real_t psi = fc.f0 + fc.b * std::log(fc.V0 / 1e-9);

   // Test various stress magnitudes and directions
   real_t tau_mags[] = {15.0e6, 18.0e6, 20.0e6, 25.0e6};
   real_t angles[] = {0.0, M_PI / 6.0, M_PI / 4.0, M_PI / 3.0, M_PI / 2.0,
                      M_PI, 3.0 * M_PI / 2.0};

   for (real_t tau_mag : tau_mags)
   {
      // Scalar solution
      real_t V_scalar = friction.SolveSlipRatePsi(tau_mag, psi, sigma_n,
                                                   eta, a);

      for (real_t angle : angles)
      {
         real_t tau_vec[2] = {tau_mag * std::cos(angle),
                              tau_mag * std::sin(angle)};
         real_t V_vec[2];
         friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);

         real_t V_vec_mag = std::sqrt(V_vec[0] * V_vec[0] +
                                      V_vec[1] * V_vec[1]);

         TEST_REL_NEAR(V_vec_mag, V_scalar, 1e-10,
                        "||V_vec|| = V_scalar for tau=" +
                        std::to_string(tau_mag / 1e6) + "MPa, angle=" +
                        std::to_string(angle));
      }
   }
}

// =============================================================================
// Test: Zero traction
// =============================================================================
void TestZeroTraction()
{
   std::cout << "\n=== Zero Traction ===\n";

   DieterichRuinaFriction friction;

   real_t tau_vec[2] = {0.0, 0.0};
   real_t V_vec[2];
   int iterations;

   friction.SolveSlipRateVectorPsi(tau_vec, 0.6, 25.0e6, 1e6, 0.01,
                                    V_vec, &iterations);

   TEST_NEAR(V_vec[0], 0.0, 1e-30, "Zero traction: V[0] = 0");
   TEST_NEAR(V_vec[1], 0.0, 1e-30, "Zero traction: V[1] = 0");
   TEST_NEAR((real_t)iterations, 0.0, 0.5, "Zero traction: 0 iterations");
}

// =============================================================================
// Test: BP5 values — verify V_init recovery with vector solver
// =============================================================================
void TestBP5Values()
{
   std::cout << "\n=== BP5 Values Recovery ===\n";

   BP5Params p;

   DieterichRuinaFriction::Constants fc;
   fc.V0 = p.V0;
   fc.f0 = p.f0;
   fc.b = p.b;
   fc.Dc = p.L0;
   DieterichRuinaFriction friction(fc);

   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

   // Outside nucleation: feed tau0_vec through vector solver
   {
      real_t x2 = 0.0, x3 = 10.0e3;
      real_t tau[2];
      p.tau0_vec(x2, x3, tau);

      // The vector solver expects the total traction as positive input.
      // tau0_vec returns negative values (anti-parallel to V).
      // The solver returns V anti-parallel to input tau.
      // So -tau fed to solver should give V in the direction of V_init.
      real_t neg_tau[2] = {-tau[0], -tau[1]};

      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(neg_tau, psi_ss, p.sigma_n, p.eta(),
                                       p.a_of_x2_x3(x2, x3), V_vec);

      // V_vec should be anti-parallel to neg_tau, i.e., parallel to tau,
      // which is anti-parallel to V_init. So V_vec should be anti-parallel
      // to V_init. Actually: -tau (positive) -> solver -> V = -V_init direction.
      // The magnitude should match.
      real_t V_mag = std::sqrt(V_vec[0] * V_vec[0] + V_vec[1] * V_vec[1]);
      real_t Vi[2];
      p.V_init_vec(x2, x3, Vi);
      real_t Vi_mag = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);

      TEST_REL_NEAR(V_mag, Vi_mag, 1e-6,
                     "BP5 outside nucleation: |V_solved| = |V_init|");
   }

   // In nucleation zone
   {
      real_t x2 = -25.0e3, x3 = 10.0e3;

      // Use L_nuc for nucleation
      DieterichRuinaFriction::Constants fc_nuc = fc;
      fc_nuc.Dc = p.L_nuc;
      DieterichRuinaFriction friction_nuc(fc_nuc);

      real_t tau[2];
      p.tau0_vec(x2, x3, tau);
      real_t neg_tau[2] = {-tau[0], -tau[1]};

      real_t V_vec[2];
      friction_nuc.SolveSlipRateVectorPsi(neg_tau, psi_ss, p.sigma_n, p.eta(),
                                           p.a_of_x2_x3(x2, x3), V_vec);

      real_t V_mag = std::sqrt(V_vec[0] * V_vec[0] + V_vec[1] * V_vec[1]);
      real_t Vi[2];
      p.V_init_vec(x2, x3, Vi);
      real_t Vi_mag = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);

      TEST_REL_NEAR(V_mag, Vi_mag, 1e-6,
                     "BP5 nucleation zone: |V_solved| = |V_init|");
   }
}

// =============================================================================
// Test: Various a values (VW and VS zones)
// =============================================================================
void TestVariousAValues()
{
   std::cout << "\n=== Various a Values ===\n";

   DieterichRuinaFriction::Constants fc;
   fc.V0 = 1.0e-6;
   fc.f0 = 0.6;
   fc.b = 0.03;
   fc.Dc = 0.14;
   DieterichRuinaFriction friction(fc);

   real_t sigma_n = 25.0e6;
   real_t eta = 2670.0 * 3464.0 * 3464.0 / (2.0 * 3464.0);
   real_t psi = fc.f0 + fc.b * std::log(fc.V0 / 1e-9);

   real_t tau_vec[2] = {15.0e6, 10.0e6};
   real_t tau_mag = std::sqrt(tau_vec[0] * tau_vec[0] +
                              tau_vec[1] * tau_vec[1]);

   // Test with a = a0 (VW zone)
   {
      real_t a = 0.004;
      real_t V_scalar = friction.SolveSlipRatePsi(tau_mag, psi, sigma_n,
                                                   eta, a);
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);
      real_t V_vec_mag = std::sqrt(V_vec[0] * V_vec[0] +
                                   V_vec[1] * V_vec[1]);
      TEST_REL_NEAR(V_vec_mag, V_scalar, 1e-10,
                     "a=0.004 (VW): vector and scalar magnitudes agree");
   }

   // Test with a = amax (VS zone)
   {
      real_t a = 0.04;
      real_t V_scalar = friction.SolveSlipRatePsi(tau_mag, psi, sigma_n,
                                                   eta, a);
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);
      real_t V_vec_mag = std::sqrt(V_vec[0] * V_vec[0] +
                                   V_vec[1] * V_vec[1]);
      TEST_REL_NEAR(V_vec_mag, V_scalar, 1e-10,
                     "a=0.04 (VS): vector and scalar magnitudes agree");
   }

   // Test with a in transition
   {
      real_t a = 0.015;
      real_t V_scalar = friction.SolveSlipRatePsi(tau_mag, psi, sigma_n,
                                                   eta, a);
      real_t V_vec[2];
      friction.SolveSlipRateVectorPsi(tau_vec, psi, sigma_n, eta, a, V_vec);
      real_t V_vec_mag = std::sqrt(V_vec[0] * V_vec[0] +
                                   V_vec[1] * V_vec[1]);
      TEST_REL_NEAR(V_vec_mag, V_scalar, 1e-10,
                     "a=0.015 (transition): vector and scalar magnitudes agree");
   }
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  Vector Friction Solver Unit Tests\n";
   std::cout << "========================================\n";

   TestDirection();
   TestMagnitudeConsistency();
   TestZeroTraction();
   TestBP5Values();
   TestVariousAValues();

   TEST_PRINT_RESULTS();
   return (num_failed > 0) ? 1 : 0;
}
