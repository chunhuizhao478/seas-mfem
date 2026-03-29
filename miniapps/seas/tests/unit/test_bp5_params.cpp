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

// Unit tests for BP5 parameters.
// Tests the spatial parameter functions a(x2,x3), L(x2,x3),
// IsNucleationZone(), V_init_vec(), and tau0_vec().

#include "test_macros.hpp"
#include "../../config/bp5_params.hpp"
#include "../../friction/dieterich_ruina.hpp"

#include <cstdlib>
#include <string>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test: Material properties
// =============================================================================
void TestMaterialProperties()
{
   std::cout << "\n=== Material Properties ===\n";

   BP5Params p;

   // mu = rho * cs^2 = 2670 * 3464^2 = 32.04 GPa
   real_t expected_mu = 2670.0 * 3464.0 * 3464.0;
   TEST_REL_NEAR(p.mu(), expected_mu, 1e-12, "mu = rho * cs^2");
   TEST_REL_NEAR(p.mu() / 1e9, 32.03862240, 1e-4, "mu ~ 32.04 GPa");

   // For nu = 0.25: lambda = 2*0.25*mu / (1 - 2*0.25) = mu
   TEST_REL_NEAR(p.lambda(), p.mu(), 1e-12, "lambda = mu for nu=0.25");

   // eta = mu / (2*cs)
   real_t expected_eta = expected_mu / (2.0 * 3464.0);
   TEST_REL_NEAR(p.eta(), expected_eta, 1e-12, "eta = mu / (2*cs)");
}

// =============================================================================
// Test: a(x2, x3) spatial function
// =============================================================================
void TestAFunction()
{
   std::cout << "\n=== a(x2, x3) Function ===\n";

   BP5Params p;

   // VW core center: x2=0, x3=10e3 (hs+ht=4e3 <= 10e3 <= 16e3=hs+ht+H)
   TEST_NEAR(p.a_of_x2_x3(0.0, 10.0e3), p.a0, 1e-15,
             "a(0, 10km) = a0 (VW core center)");

   // VW core edges
   TEST_NEAR(p.a_of_x2_x3(0.0, 4.0e3), p.a0, 1e-15,
             "a(0, 4km) = a0 (VW top edge, x3=hs+ht)");
   TEST_NEAR(p.a_of_x2_x3(0.0, 16.0e3), p.a0, 1e-15,
             "a(0, 16km) = a0 (VW bottom edge, x3=hs+ht+H)");
   TEST_NEAR(p.a_of_x2_x3(29.9e3, 10.0e3), p.a0, 1e-15,
             "a(29.9km, 10km) = a0 (VW right edge)");
   TEST_NEAR(p.a_of_x2_x3(-29.9e3, 10.0e3), p.a0, 1e-15,
             "a(-29.9km, 10km) = a0 (VW left edge)");

   // VS zones
   TEST_NEAR(p.a_of_x2_x3(0.0, 0.0), p.amax, 1e-15,
             "a(0, 0) = amax (shallow VS, x3 <= hs)");
   TEST_NEAR(p.a_of_x2_x3(0.0, 1.0e3), p.amax, 1e-15,
             "a(0, 1km) = amax (shallow VS)");
   TEST_NEAR(p.a_of_x2_x3(0.0, 2.0e3), p.amax, 1e-15,
             "a(0, 2km) = amax (x3 = hs boundary)");
   TEST_NEAR(p.a_of_x2_x3(0.0, 20.0e3), p.amax, 1e-15,
             "a(0, 20km) = amax (deep VS, x3 >= hs+2*ht+H)");
   TEST_NEAR(p.a_of_x2_x3(0.0, 35.0e3), p.amax, 1e-15,
             "a(0, 35km) = amax (deep VS)");
   TEST_NEAR(p.a_of_x2_x3(35.0e3, 10.0e3), p.amax, 1e-15,
             "a(35km, 10km) = amax (far along-strike VS)");

   // Transition zone: x3 between hs and hs+ht (depth transition)
   real_t a_mid = p.a_of_x2_x3(0.0, 3.0e3); // x3=3km, between hs=2km and hs+ht=4km
   TEST_ASSERT(a_mid > p.a0 && a_mid < p.amax,
               "a(0, 3km) in transition (a0 < a < amax)");

   // Transition zone: along-strike transition
   real_t a_strike = p.a_of_x2_x3(31.0e3, 10.0e3); // |x2|=31km > l/2=30km
   TEST_ASSERT(a_strike > p.a0 && a_strike < p.amax,
               "a(31km, 10km) in transition (along-strike)");

   // Transition continuity: just outside VW in transition (r→0)
   // x3=3999 < hs+ht=4000, so in transition with r close to 0 → a ~ a0
   real_t a_near_vw = p.a_of_x2_x3(0.0, 4.0e3 - 1.0);
   TEST_NEAR(a_near_vw, p.a0, 1e-2 * (p.amax - p.a0),
             "a near VW boundary ~ a0 (transition continuity)");

   // Transition continuity: just outside VS in transition (r→1)
   // x3=2001 > hs=2000, so in transition with r close to 1 → a ~ amax
   real_t a_near_vs = p.a_of_x2_x3(0.0, 2.0e3 + 1.0);
   TEST_NEAR(a_near_vs, p.amax, 1e-2 * (p.amax - p.a0),
             "a near VS boundary ~ amax (transition continuity)");

   // Symmetry: a(x2, x3) = a(-x2, x3)
   TEST_NEAR(p.a_of_x2_x3(20.0e3, 8.0e3), p.a_of_x2_x3(-20.0e3, 8.0e3),
             1e-15, "a is symmetric in x2");
}

// =============================================================================
// Test: L(x2, x3) and nucleation zone
// =============================================================================
void TestNucleationAndL()
{
   std::cout << "\n=== Nucleation Zone and L(x2, x3) ===\n";

   BP5Params p;

   // Nucleation zone: hs+ht <= x3 <= hs+ht+H AND -l/2 <= x2 <= -l/2+w
   // i.e., 4e3 <= x3 <= 16e3 AND -30e3 <= x2 <= -18e3

   // Inside nucleation zone
   TEST_ASSERT(p.IsNucleationZone(-25.0e3, 10.0e3),
               "(-25km, 10km) is in nucleation zone");
   TEST_ASSERT(p.IsNucleationZone(-30.0e3, 4.0e3),
               "(-30km, 4km) is in nucleation zone (corner)");
   TEST_ASSERT(p.IsNucleationZone(-18.0e3, 16.0e3),
               "(-18km, 16km) is in nucleation zone (corner)");

   // Outside nucleation zone
   TEST_ASSERT(!p.IsNucleationZone(0.0, 10.0e3),
               "(0, 10km) is NOT in nucleation zone");
   TEST_ASSERT(!p.IsNucleationZone(-25.0e3, 3.0e3),
               "(-25km, 3km) is NOT in nucleation zone (too shallow)");
   TEST_ASSERT(!p.IsNucleationZone(-25.0e3, 17.0e3),
               "(-25km, 17km) is NOT in nucleation zone (too deep)");
   TEST_ASSERT(!p.IsNucleationZone(-17.0e3, 10.0e3),
               "(-17km, 10km) is NOT in nucleation zone (x2 > -l/2+w)");
   TEST_ASSERT(!p.IsNucleationZone(-31.0e3, 10.0e3),
               "(-31km, 10km) is NOT in nucleation zone (x2 < -l/2)");
   TEST_ASSERT(!p.IsNucleationZone(25.0e3, 10.0e3),
               "(25km, 10km) is NOT in nucleation zone (positive x2)");

   // L function
   TEST_NEAR(p.L_of_x2_x3(-25.0e3, 10.0e3), p.L_nuc, 1e-15,
             "L(-25km, 10km) = L_nuc (in nucleation)");
   TEST_NEAR(p.L_of_x2_x3(0.0, 10.0e3), p.L0, 1e-15,
             "L(0, 10km) = L0 (outside nucleation)");
   TEST_NEAR(p.L_of_x2_x3(0.0, 0.0), p.L0, 1e-15,
             "L(0, 0) = L0 (surface)");

   // Tandem stock semantics use bp5_outside with eps=1e-3, so points lying
   // just outside the nominal rectangle by less than eps are still classified
   // as inside. Switching eps back to 0 reproduces bp5_exact behavior.
   TEST_ASSERT(p.IsNucleationZone(-17.9999995e3, 10.0e3),
               "bp5_outside: point 0.5mm outside strike boundary is still inside");
   TEST_ASSERT(p.IsNucleationZone(-25.0e3, 16.0000005e3),
               "bp5_outside: point 0.5mm below depth boundary is still inside");

   p.nucleation_eps = 0.0;
   TEST_ASSERT(!p.IsNucleationZone(-17.9999995e3, 10.0e3),
               "bp5_exact: point outside strike boundary is excluded");
   TEST_ASSERT(!p.IsNucleationZone(-25.0e3, 16.0000005e3),
               "bp5_exact: point outside depth boundary is excluded");
}

// =============================================================================
// Test: V_init_vec
// =============================================================================
void TestVinitVec()
{
   std::cout << "\n=== V_init_vec ===\n";

   BP5Params p;
   real_t V[2];

   // Outside nucleation zone
   p.V_init_vec(0.0, 10.0e3, V);
   TEST_NEAR(V[0], p.V_zero, 1e-30, "V[0] = V_zero outside nucleation");
   TEST_NEAR(V[1], p.V_init, 1e-30, "V[1] = V_init outside nucleation");

   // In nucleation zone
   p.V_init_vec(-25.0e3, 10.0e3, V);
   TEST_NEAR(V[0], p.V_zero, 1e-30, "V[0] = V_zero in nucleation");
   TEST_NEAR(V[1], p.V_nuc, 1e-15, "V[1] = V_nuc in nucleation");

   // Check magnitude: outside nucleation, |V| ~ Vp (V_zero negligible)
   p.V_init_vec(0.0, 10.0e3, V);
   real_t V_mag = std::sqrt(V[0] * V[0] + V[1] * V[1]);
   TEST_REL_NEAR(V_mag, p.V_init, 1e-10, "|V| ~ V_init outside nucleation");

   // In nucleation, |V| ~ V_nuc
   p.V_init_vec(-25.0e3, 10.0e3, V);
   V_mag = std::sqrt(V[0] * V[0] + V[1] * V[1]);
   TEST_REL_NEAR(V_mag, p.V_nuc, 1e-10, "|V| ~ V_nuc in nucleation");
}

// =============================================================================
// Test: tau0_vec pre-stress
// =============================================================================
void TestTau0Vec()
{
   std::cout << "\n=== tau0_vec Pre-stress ===\n";

   BP5Params p;
   real_t tau[2];

   // Pre-stress outside nucleation zone
   p.tau0_vec(0.0, 10.0e3, tau);

   // tau should be finite and non-zero
   real_t tau_mag = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);
   TEST_ASSERT(tau_mag > 0.0, "|tau0| > 0 outside nucleation");
   TEST_ASSERT(tau_mag < 1e8, "|tau0| < 100 MPa (physically reasonable)");

   // tau is negated (Tandem convention: traction on -Y face for right-lateral)
   // V_init = (V_zero, Vp) → tau = -tau0 * V/|V| → nearly pure negative strike
   TEST_ASSERT(tau[1] < 0.0, "tau[1] < 0 (negated, Tandem convention)");

   // Pre-stress in nucleation zone
   p.tau0_vec(-25.0e3, 10.0e3, tau);
   real_t tau_nuc_mag = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);
   TEST_ASSERT(tau_nuc_mag > 0.0, "|tau0| > 0 in nucleation");

   // Nucleation zone should have higher stress (higher V_nuc)
   TEST_ASSERT(tau_nuc_mag > tau_mag,
               "|tau0_nuc| > |tau0_outside| (higher stress in nucleation)");

   // Self-consistency: feed tau0 into friction solver, should recover V_init
   std::cout << "\n  --- tau0_vec self-consistency check ---\n";

   // Create friction law with BP5 constants
   DieterichRuinaFriction::Constants fc;
   fc.V0 = p.V0;
   fc.f0 = p.f0;
   fc.b = p.b;
   fc.Dc = p.L0;  // Use L0 for outside nucleation
   DieterichRuinaFriction friction(fc);

   // Outside nucleation zone
   {
      real_t x2 = 0.0, x3 = 10.0e3;
      real_t Vi[2];
      p.V_init_vec(x2, x3, Vi);
      real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);

      p.tau0_vec(x2, x3, tau);

      // tau0_vec gives the pre-stress. In the friction solver, the input traction
      // is -tau (since tau opposes slip, the traction driving slip is -tau).
      // Actually: tau0_vec is the pre-stress that, combined with friction,
      // produces the initial slip rate. The friction solver expects positive
      // traction and returns slip rate anti-parallel to it.
      // So we feed |tau| and should get |V_init| back.
      real_t tau_total_mag = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);

      // Compute psi at steady state with Vp
      real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

      // The scalar friction solver should return V_init_magnitude
      real_t V_solved = friction.SolveSlipRatePsi(
         tau_total_mag, psi_ss, p.sigma_n, p.eta(), p.a_of_x2_x3(x2, x3));

      TEST_REL_NEAR(V_solved, Vi_abs, 1e-6,
                     "Friction solver recovers |V_init| outside nucleation");
   }

   // In nucleation zone: with Tandem defaults (delta_tau_factor = 0),
   // tau_pre is computed for exact equilibrium with Vi and psi_ss.
   // The friction solver should recover V = V_nuc exactly.
   //
   // Tandem reference: bp5.lua lines 94-103 — tau_pre = sn*a*asinh(...) + eta*Vi
   // with NO delta_tau addition.
   {
      real_t x2 = -25.0e3, x3 = 10.0e3;
      real_t Vi[2];
      p.V_init_vec(x2, x3, Vi);
      real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);

      p.tau0_vec(x2, x3, tau);
      real_t tau_total_mag = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);

      // Use L_nuc for Dc in nucleation
      DieterichRuinaFriction::Constants fc_nuc;
      fc_nuc.V0 = p.V0;
      fc_nuc.f0 = p.f0;
      fc_nuc.b = p.b;
      fc_nuc.Dc = p.L_nuc;
      DieterichRuinaFriction friction_nuc(fc_nuc);

      real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

      real_t V_solved = friction_nuc.SolveSlipRatePsi(
         tau_total_mag, psi_ss, p.sigma_n, p.eta(), p.a_of_x2_x3(x2, x3));

      // With Tandem defaults (delta_tau_factor=0), V_solved == V_nuc (equilibrium)
      TEST_REL_NEAR(V_solved, Vi_abs, 1e-6,
                     "Friction solver recovers V_nuc in nucleation (Tandem equilibrium)");

      // Also verify via InitialStatePsi round-trip:
      // psi from InitialStatePsi should equal psi_ss
      real_t psi_from_stress = friction_nuc.InitialStatePsi(
         tau_total_mag, Vi_abs, p.sigma_n, p.eta(), p.a_of_x2_x3(x2, x3));
      TEST_REL_NEAR(psi_from_stress, psi_ss, 1e-10,
                     "InitialStatePsi recovers psi_ss at nucleation (Tandem equilibrium)");
   }
}

// =============================================================================
// Test: Quantitative pre-stress magnitude against SCEC analytical values
// =============================================================================
void TestTau0VecQuantitative()
{
   std::cout << "\n=== tau0_vec Quantitative (SCEC Spec) ===\n";

   BP5Params p;
   real_t tau[2];

   // Analytical computation:
   // psi_ss = f0 + b * ln(V0/Vp) = 0.6 + 0.03 * ln(1e-6/1e-9) = 0.807233
   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);
   real_t eta = p.eta();

   // ----------- VW zone (a=0.004), outside nucleation -----------
   // x2=0, x3=10e3 → a=0.004
   {
      real_t x2 = 0.0, x3 = 10.0e3;
      real_t a = p.a_of_x2_x3(x2, x3);
      TEST_NEAR(a, 0.004, 1e-15, "a(0, 10km) = 0.004 (VW zone)");

      real_t Vi_abs = std::sqrt(p.V_zero * p.V_zero + p.V_init * p.V_init);
      real_t e = std::exp(psi_ss / a);
      real_t expected_tau = p.sigma_n * a *
         std::asinh((Vi_abs / (2.0 * p.V0)) * e) + eta * Vi_abs;

      p.tau0_vec(x2, x3, tau);
      real_t tau_mag = std::sqrt(tau[0]*tau[0] + tau[1]*tau[1]);

      // Expected ~19.49 MPa
      TEST_REL_NEAR(tau_mag, expected_tau, 1e-10,
                     "tau0 magnitude matches analytical (VW, a=0.004)");
      TEST_REL_NEAR(tau_mag / 1e6, 19.49, 1e-2,
                     "tau0 ~ 19.49 MPa (VW zone)");
   }

   // ----------- VS zone (a=0.04), outside nucleation -----------
   // x2=0, x3=30e3 → a=0.04 (deep VS)
   {
      real_t x2 = 0.0, x3 = 30.0e3;
      real_t a = p.a_of_x2_x3(x2, x3);
      TEST_NEAR(a, 0.04, 1e-15, "a(0, 30km) = 0.04 (VS zone)");

      real_t Vi_abs = std::sqrt(p.V_zero * p.V_zero + p.V_init * p.V_init);
      real_t e = std::exp(psi_ss / a);
      real_t expected_tau = p.sigma_n * a *
         std::asinh((Vi_abs / (2.0 * p.V0)) * e) + eta * Vi_abs;

      p.tau0_vec(x2, x3, tau);
      real_t tau_mag = std::sqrt(tau[0]*tau[0] + tau[1]*tau[1]);

      // Expected ~13.27 MPa
      TEST_REL_NEAR(tau_mag, expected_tau, 1e-10,
                     "tau0 magnitude matches analytical (VS, a=0.04)");
      TEST_REL_NEAR(tau_mag / 1e6, 13.27, 1e-2,
                     "tau0 ~ 13.27 MPa (VS zone)");
   }

   // ----------- Nucleation zone (a=0.004), Tandem defaults (NO delta_tau) ------
   // x2=-25e3, x3=10e3 → nucleation zone, a=0.004
   // With delta_tau_factor=0 (Tandem default):
   //   tau = sn * a * asinh(Vi/(2V0) * exp(psi_ss/a)) + eta * Vi
   // No delta_tau added. Matches Tandem bp5.lua lines 94-102.
   {
      real_t x2 = -25.0e3, x3 = 10.0e3;
      real_t a = p.a_of_x2_x3(x2, x3);
      TEST_NEAR(a, 0.004, 1e-15, "a(-25km, 10km) = 0.004 (nucleation in VW)");
      TEST_ASSERT(p.IsNucleationZone(x2, x3), "(-25km, 10km) is nucleation");

      real_t Vi_abs = std::sqrt(p.V_zero * p.V_zero + p.V_nuc * p.V_nuc);
      real_t e = std::exp(psi_ss / a);

      // Equilibrium tau: sn * a * asinh(Vi/(2V0) * e) + eta * Vi
      // NO delta_tau (delta_tau_factor = 0)
      real_t expected_tau = p.sigma_n * a *
         std::asinh((Vi_abs / (2.0 * p.V0)) * e) + eta * Vi_abs;

      p.tau0_vec(x2, x3, tau);
      real_t tau_mag = std::sqrt(tau[0]*tau[0] + tau[1]*tau[1]);

      TEST_REL_NEAR(tau_mag, expected_tau, 1e-10,
                     "tau0 magnitude matches analytical (nucleation, no delta_tau)");

      // Verify delta_tau_factor is 0 (Tandem default)
      TEST_NEAR(p.delta_tau_factor, 0.0, 1e-15,
                "delta_tau_factor = 0 (Tandem default, no overstress)");

      // Verify nucleation tau > non-nucleation tau
      // (because Vi_nuc=0.01 >> Vi_init=1e-9, so friction is higher)
      real_t tau_non_nuc[2];
      p.tau0_vec(0.0, 10.0e3, tau_non_nuc);
      real_t tau_non_nuc_mag = std::sqrt(tau_non_nuc[0]*tau_non_nuc[0] +
                                          tau_non_nuc[1]*tau_non_nuc[1]);
      TEST_ASSERT(tau_mag > tau_non_nuc_mag,
                  "Nucleation tau > non-nucleation tau (higher V_nuc)");
   }
}

// =============================================================================
// Test: Pre-stress → friction solver → V_init round-trip at multiple a values
// =============================================================================
void TestPreStressFrictionRoundTrip()
{
   std::cout << "\n=== Pre-stress Friction Round-Trip ===\n";

   BP5Params p;
   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

   // Test at multiple (x2, x3) points covering VW, transition, and VS zones
   struct TestPoint { real_t x2, x3; const char *label; };
   TestPoint points[] = {
      {      0.0, 10.0e3, "VW center (a=0.004)"},
      {      0.0,  3.0e3, "Shallow transition"},
      {      0.0, 19.0e3, "Deep transition"},
      {      0.0, 30.0e3, "Deep VS (a=0.04)"},
      {  31.0e3, 10.0e3, "Strike transition"},
      { -10.0e3, 10.0e3, "VW off-center"},
   };

   for (auto &pt : points)
   {
      real_t a = p.a_of_x2_x3(pt.x2, pt.x3);
      real_t Vi_abs = std::sqrt(p.V_zero * p.V_zero + p.V_init * p.V_init);

      real_t tau[2];
      p.tau0_vec(pt.x2, pt.x3, tau);
      real_t tau_mag = std::sqrt(tau[0]*tau[0] + tau[1]*tau[1]);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = p.V0; fc.f0 = p.f0; fc.b = p.b;
      fc.Dc = p.L_of_x2_x3(pt.x2, pt.x3);
      DieterichRuinaFriction friction(fc);

      real_t V_solved = friction.SolveSlipRatePsi(
         tau_mag, psi_ss, p.sigma_n, p.eta(), a);

      std::string msg = "Round-trip V at " + std::string(pt.label) +
                        " (a=" + std::to_string(a) + ")";
      TEST_REL_NEAR(V_solved, Vi_abs, 1e-6, msg.c_str());
   }
}

// =============================================================================
// Test: Print function (just verify it doesn't crash)
// =============================================================================
void TestPrint()
{
   std::cout << "\n=== Print ===\n";

   BP5Params p;
   std::ostringstream oss;
   p.Print(oss);
   std::string output = oss.str();

   TEST_ASSERT(output.find("BP5") != std::string::npos,
               "Print output contains 'BP5'");
   TEST_ASSERT(output.find("32.0") != std::string::npos,
               "Print output contains mu ~ 32.0 GPa");
   TEST_ASSERT(output.find("25") != std::string::npos,
               "Print output contains sigma_n = 25 MPa");

   num_tests++;
   num_passed++;
   std::cout << "  PASSED: Print completes without error\n";
}

// =============================================================================
// Test: a() depth symmetry within VW zone
// =============================================================================
void TestADepthSymmetry()
{
   std::cout << "\n=== a() Depth Symmetry ===\n";

   BP5Params p;

   // VW core center depth: hs + ht + H/2 = 2e3 + 2e3 + 6e3 = 10e3
   real_t center = p.hs + p.ht + p.H / 2.0;

   // Points symmetric about center within VW zone
   for (real_t d : {1.0e3, 2.0e3, 4.0e3, 5.9e3})
   {
      real_t a_above = p.a_of_x2_x3(0.0, center - d);
      real_t a_below = p.a_of_x2_x3(0.0, center + d);
      std::string msg = "a symmetric about center for d=" +
                        std::to_string(d / 1e3) + "km";
      TEST_NEAR(a_above, a_below, 1e-15, msg.c_str());
   }
}

// =============================================================================
// Test: Zone boundary edge cases
// =============================================================================
void TestZoneBoundaryEdgeCases()
{
   std::cout << "\n=== Zone Boundary Edge Cases ===\n";

   BP5Params p;

   // Exact VW boundaries: x3=hs+ht=4e3, x3=hs+ht+H=16e3
   TEST_NEAR(p.a_of_x2_x3(0.0, p.hs + p.ht), p.a0, 1e-15,
             "a at x3=hs+ht (VW top edge) = a0");
   TEST_NEAR(p.a_of_x2_x3(0.0, p.hs + p.ht + p.H), p.a0, 1e-15,
             "a at x3=hs+ht+H (VW bottom edge) = a0");

   // Exact along-strike VW boundary: |x2|=l_vw/2
   TEST_NEAR(p.a_of_x2_x3(p.l_vw / 2.0, 10.0e3), p.a0, 1e-15,
             "a at |x2|=l_vw/2 (VW strike edge) = a0");
   TEST_NEAR(p.a_of_x2_x3(-p.l_vw / 2.0, 10.0e3), p.a0, 1e-15,
             "a at |x2|=-l_vw/2 (VW strike edge) = a0");

   // Exact VS boundaries: x3=hs=2e3, x3=hs+2*ht+H=18e3, |x2|=l_vw/2+ht
   TEST_NEAR(p.a_of_x2_x3(0.0, p.hs), p.amax, 1e-15,
             "a at x3=hs (shallow VS edge) = amax");
   TEST_NEAR(p.a_of_x2_x3(0.0, p.hs + 2.0 * p.ht + p.H), p.amax, 1e-15,
             "a at x3=hs+2ht+H (deep VS edge) = amax");
   TEST_NEAR(p.a_of_x2_x3(p.l_vw / 2.0 + p.ht, 10.0e3), p.amax, 1e-15,
             "a at |x2|=l_vw/2+ht (far strike VS edge) = amax");
}

// =============================================================================
// Test: tau0_vec at VS zone (self-consistency with amax)
// =============================================================================
void TestTau0VecVSZone()
{
   std::cout << "\n=== tau0_vec at VS Zone ===\n";

   BP5Params p;
   real_t tau[2];

   // Point in deep VS zone: (0, 30e3) where a=amax
   p.tau0_vec(0.0, 30.0e3, tau);
   real_t tau_mag = std::sqrt(tau[0] * tau[0] + tau[1] * tau[1]);

   TEST_ASSERT(tau_mag > 0.0, "|tau0| > 0 in VS zone");
   TEST_ASSERT(tau_mag < 1e8, "|tau0| < 100 MPa in VS zone");
   TEST_ASSERT(tau[1] < 0.0, "tau[1] < 0 in VS zone (negated, Tandem convention)");

   // Verify a=amax at this point
   TEST_NEAR(p.a_of_x2_x3(0.0, 30.0e3), p.amax, 1e-15,
             "a(0, 30km) = amax (VS zone confirmed)");

   // Self-consistency: friction solver should recover V_init
   DieterichRuinaFriction::Constants fc;
   fc.V0 = p.V0;
   fc.f0 = p.f0;
   fc.b = p.b;
   fc.Dc = p.L0;
   DieterichRuinaFriction friction(fc);

   real_t Vi[2];
   p.V_init_vec(0.0, 30.0e3, Vi);
   real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

   real_t V_solved = friction.SolveSlipRatePsi(
      tau_mag, psi_ss, p.sigma_n, p.eta(), p.amax);

   TEST_REL_NEAR(V_solved, Vi_abs, 1e-6,
                  "Friction recovers |V_init| in VS zone (amax)");
}

// =============================================================================
// Test: tau0_vec direction (tau[0]/tau[1] ~ Vi[0]/Vi[1])
// =============================================================================
void TestTau0VecDirection()
{
   std::cout << "\n=== tau0_vec Direction ===\n";

   BP5Params p;
   real_t tau[2], Vi[2];

   // Outside nucleation zone
   p.V_init_vec(0.0, 10.0e3, Vi);
   p.tau0_vec(0.0, 10.0e3, tau);

   // tau should be parallel to Vi: tau = |tau| * Vi/|Vi|
   // So tau[0]/tau[1] should equal Vi[0]/Vi[1]
   real_t Vi_ratio = Vi[0] / Vi[1];
   real_t tau_ratio = tau[0] / tau[1];
   TEST_REL_NEAR(tau_ratio, Vi_ratio, 1e-10,
                  "tau direction matches V_init direction (outside nuc)");

   // In nucleation zone
   p.V_init_vec(-25.0e3, 10.0e3, Vi);
   p.tau0_vec(-25.0e3, 10.0e3, tau);

   Vi_ratio = Vi[0] / Vi[1];
   tau_ratio = tau[0] / tau[1];
   TEST_REL_NEAR(tau_ratio, Vi_ratio, 1e-10,
                  "tau direction matches V_init direction (in nuc)");
}

// =============================================================================
// Test: psi_init() helper
// =============================================================================
void TestPsiInit()
{
   std::cout << "\n=== psi_init() ===\n";

   BP5Params p;

   real_t expected = p.f0 + p.b * std::log(p.V0 / p.V_init);
   TEST_NEAR(p.psi_init(), expected, 1e-15,
             "psi_init = f0 + b*ln(V0/V_init)");

   // Verify it's a reasonable value (should be positive)
   TEST_ASSERT(p.psi_init() > 0.0, "psi_init > 0");
   TEST_ASSERT(p.psi_init() > p.f0, "psi_init > f0 (since V0 > V_init)");
}

// =============================================================================
// Test: Validate() on default params
// =============================================================================
void TestValidate()
{
   std::cout << "\n=== Validate() ===\n";

   BP5Params p;

   // Default params should pass validation without throwing
   bool passed = true;
   try
   {
      p.Validate();
   }
   catch (...)
   {
      passed = false;
   }
   TEST_ASSERT(passed, "Default BP5 params pass Validate()");

   // Verify the individual conditions
   TEST_ASSERT(p.a0 < p.b, "a0 < b (velocity-weakening condition)");
   TEST_ASSERT(p.amax > p.b, "amax > b (velocity-strengthening condition)");
   TEST_ASSERT(p.L_nuc < p.L0, "L_nuc < L0");
   TEST_ASSERT(p.sigma_n > 0, "sigma_n > 0");
   TEST_ASSERT(p.Vp > 0, "Vp > 0");
   TEST_ASSERT(p.V_nuc >= p.V_init, "V_nuc >= V_init");
}

// =============================================================================
// Test: Tandem initialization produces exact equilibrium everywhere
// =============================================================================
void TestTandemInitEquilibrium()
{
   std::cout << "\n=== Tandem Init Equilibrium (V = V_init at all DOFs) ===\n";

   BP5Params p;
   // Verify defaults match Tandem
   TEST_NEAR(p.V_nuc, 0.01, 1e-15, "V_nuc = 0.01 (Tandem default)");
   TEST_NEAR(p.delta_tau_factor, 0.0, 1e-15,
             "delta_tau_factor = 0.0 (Tandem default)");

   DieterichRuinaFriction::Constants fc;
   fc.V0 = p.V0; fc.f0 = p.f0; fc.b = p.b; fc.Dc = p.L0;
   DieterichRuinaFriction friction(fc);

   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

   // Test at multiple points: VW center, nucleation zone, VS zone, transition
   struct TestPoint { real_t x2, x3; const char *label; bool is_nuc; };
   TestPoint points[] = {
      {      0.0, 10.0e3, "VW center", false},
      { -25.0e3, 10.0e3, "nucleation zone center", true},
      { -30.0e3,  4.0e3, "nucleation zone corner", true},
      {      0.0, 30.0e3, "deep VS", false},
      {      0.0,  3.0e3, "shallow transition", false},
      {  31.0e3, 10.0e3, "strike transition", false},
   };

   for (auto &pt : points)
   {
      real_t Vi[2], tau[2];
      p.V_init_vec(pt.x2, pt.x3, Vi);
      p.tau0_vec(pt.x2, pt.x3, tau);

      real_t Vi_abs = std::sqrt(Vi[0]*Vi[0] + Vi[1]*Vi[1]);
      real_t tau_mag = std::sqrt(tau[0]*tau[0] + tau[1]*tau[1]);
      real_t a = p.a_of_x2_x3(pt.x2, pt.x3);
      real_t Dc = p.L_of_x2_x3(pt.x2, pt.x3);

      // Set correct Dc for this point
      DieterichRuinaFriction::Constants fc_pt;
      fc_pt.V0 = p.V0; fc_pt.f0 = p.f0; fc_pt.b = p.b; fc_pt.Dc = Dc;
      DieterichRuinaFriction friction_pt(fc_pt);

      // With psi_ss, friction solver should return V = Vi (equilibrium)
      real_t V_solved = friction_pt.SolveSlipRatePsi(
         tau_mag, psi_ss, p.sigma_n, p.eta(), a);

      std::string msg = "V = V_init at " + std::string(pt.label);
      TEST_REL_NEAR(V_solved, Vi_abs, 1e-6, msg.c_str());

      // Also verify InitialStatePsi gives psi_ss
      real_t psi_computed = friction_pt.InitialStatePsi(
         tau_mag, Vi_abs, p.sigma_n, p.eta(), a);
      msg = "psi_init = psi_ss at " + std::string(pt.label);
      TEST_REL_NEAR(psi_computed, psi_ss, 1e-10, msg.c_str());
   }
}

// =============================================================================
// Test: Tandem tau_pre matches Tandem's bp5.lua formula exactly
// =============================================================================
void TestTandemTauPreFormula()
{
   std::cout << "\n=== Tandem tau_pre Formula Match ===\n";

   BP5Params p;
   real_t eta = p.eta();

   // Tandem's bp5.lua formula:
   //   e = exp((f0 + b * log(V0 / Vp)) / a)
   //   tau0 = sn * a * asinh((Vi2 / (2*V0)) * e) + eta * Vi2
   //   tau_pre = (-tau0 * Vi1/Vi, -tau0 * Vi2/Vi)

   // Test at nucleation zone center: Vi = (1e-20, 0.01)
   {
      real_t x2 = -25.0e3, x3 = 10.0e3;
      real_t a = p.a_of_x2_x3(x2, x3);

      real_t Vi1 = p.V_zero;  // 1e-20
      real_t Vi2 = p.V_nuc;   // 0.01
      real_t Vi = std::sqrt(Vi1*Vi1 + Vi2*Vi2);

      // Tandem's formula
      real_t e = std::exp((p.f0 + p.b * std::log(p.V0 / p.Vp)) / a);
      real_t tau0_tandem = p.sigma_n * a * std::asinh((Vi2 / (2.0 * p.V0)) * e)
                          + eta * Vi2;
      real_t tau_tandem[2] = {-tau0_tandem * Vi1 / Vi,
                               -tau0_tandem * Vi2 / Vi};

      // Our formula
      real_t tau_ours[2];
      p.tau0_vec(x2, x3, tau_ours);

      TEST_REL_NEAR(tau_ours[0], tau_tandem[0], 1e-10,
                     "tau[0] matches Tandem formula (nucleation)");
      TEST_REL_NEAR(tau_ours[1], tau_tandem[1], 1e-10,
                     "tau[1] matches Tandem formula (nucleation)");
   }

   // Test outside nucleation: Vi = (1e-20, 1e-9)
   {
      real_t x2 = 0.0, x3 = 10.0e3;
      real_t a = p.a_of_x2_x3(x2, x3);

      real_t Vi1 = p.V_zero;
      real_t Vi2 = p.V_init;
      real_t Vi = std::sqrt(Vi1*Vi1 + Vi2*Vi2);

      real_t e = std::exp((p.f0 + p.b * std::log(p.V0 / p.Vp)) / a);
      real_t tau0_tandem = p.sigma_n * a * std::asinh((Vi2 / (2.0 * p.V0)) * e)
                          + eta * Vi2;
      real_t tau_tandem[2] = {-tau0_tandem * Vi1 / Vi,
                               -tau0_tandem * Vi2 / Vi};

      real_t tau_ours[2];
      p.tau0_vec(x2, x3, tau_ours);

      TEST_REL_NEAR(tau_ours[0], tau_tandem[0], 1e-10,
                     "tau[0] matches Tandem formula (outside nuc)");
      TEST_REL_NEAR(tau_ours[1], tau_tandem[1], 1e-10,
                     "tau[1] matches Tandem formula (outside nuc)");
   }
}

// =============================================================================
// Test: SCEC mode can be recovered with explicit parameter overrides
// =============================================================================
void TestScecModeOverride()
{
   std::cout << "\n=== SCEC Mode Override ===\n";

   // Create params with SCEC overrides
   BP5Params p;
   p.V_nuc = 0.03;
   p.delta_tau_factor = 1.0;
   p.Validate();

   real_t tau_nuc[2];
   p.tau0_vec(-25.0e3, 10.0e3, tau_nuc);
   real_t tau_nuc_mag = std::sqrt(tau_nuc[0]*tau_nuc[0] + tau_nuc[1]*tau_nuc[1]);

   real_t psi_ss = p.f0 + p.b * std::log(p.V0 / p.Vp);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = p.V0; fc.f0 = p.f0; fc.b = p.b; fc.Dc = p.L_nuc;
   DieterichRuinaFriction friction(fc);

   // With SCEC psi (fixed = psi_ss) and delta_tau, V > V_nuc (overstressed)
   real_t V_solved = friction.SolveSlipRatePsi(
      tau_nuc_mag, psi_ss, p.sigma_n, p.eta(), p.a_of_x2_x3(-25.0e3, 10.0e3));

   real_t Vi_abs = std::sqrt(p.V_zero * p.V_zero + p.V_nuc * p.V_nuc);
   TEST_ASSERT(V_solved > Vi_abs,
               "SCEC mode: V > V_nuc (delta_tau overstress present)");

   // Verify the overstress is significant (~2x V_nuc)
   TEST_ASSERT(V_solved > 1.5 * Vi_abs,
               "SCEC mode: V > 1.5 * V_nuc (significant overstress)");

   // Verify InitialStatePsi absorbs the overstress (Tandem psi_init)
   real_t psi_absorbed = friction.InitialStatePsi(
      tau_nuc_mag, Vi_abs, p.sigma_n, p.eta(), p.a_of_x2_x3(-25.0e3, 10.0e3));
   TEST_ASSERT(psi_absorbed > psi_ss,
               "SCEC mode: InitialStatePsi > psi_ss (absorbs delta_tau)");

   // With the absorbed psi, V should recover V_nuc
   real_t V_absorbed = friction.SolveSlipRatePsi(
      tau_nuc_mag, psi_absorbed, p.sigma_n, p.eta(),
      p.a_of_x2_x3(-25.0e3, 10.0e3));
   TEST_REL_NEAR(V_absorbed, Vi_abs, 1e-6,
                  "SCEC mode: absorbed psi recovers V_nuc");
}

// =============================================================================
// Test: psi_init() numerical value matches Tandem
// =============================================================================
void TestPsiInitNumerical()
{
   std::cout << "\n=== psi_init Numerical Value ===\n";

   BP5Params p;

   // Both Tandem and SCEC compute the same value:
   // psi = f0 + b*ln(V0/Vp) = 0.6 + 0.03*ln(1e-6/1e-9) = 0.6 + 0.03*6.9078 = 0.8072
   real_t expected = 0.6 + 0.03 * std::log(1e-6 / 1e-9);
   TEST_REL_NEAR(p.psi_init(), expected, 1e-12,
                  "psi_init = f0 + b*ln(V0/V_init) = 0.8072");

   // Verify numerical value to 4 decimal places
   TEST_REL_NEAR(p.psi_init(), 0.80723, 1e-4,
                  "psi_init ~ 0.8072 (Tandem verified value)");
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  BP5 Parameters Unit Tests\n";
   std::cout << "========================================\n";

   TestMaterialProperties();
   TestAFunction();
   TestNucleationAndL();
   TestVinitVec();
   TestTau0Vec();
   TestTau0VecQuantitative();
   TestPreStressFrictionRoundTrip();
   TestPrint();
   TestADepthSymmetry();
   TestZoneBoundaryEdgeCases();
   TestTau0VecVSZone();
   TestTau0VecDirection();
   TestPsiInit();
   TestValidate();
   TestTandemInitEquilibrium();
   TestTandemTauPreFormula();
   TestScecModeOverride();
   TestPsiInitNumerical();

   TEST_PRINT_RESULTS();
   return (num_failed > 0) ? 1 : 0;
}
