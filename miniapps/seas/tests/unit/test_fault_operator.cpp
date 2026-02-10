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

// Unit tests for Phase 3: Fault Operator
// Tests for FaultGeometry and RateStateFaultOperator classes:
// - Fault geometry extraction and depth-dependent parameters
// - State initialization (PreInit, Init)
// - RHS computation
// - Stress equilibrium verification
// - Pre-stress tau0 verification (~26.546 MPa)

#include "mfem.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <memory>

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
// Helper: Create a simple test mesh with fault boundary
// =============================================================================

/// Create a simple 2D mesh for testing
/// Domain: x in [-Lx, Lx], z in [-Lz, 0]
/// Fault at x = 0 is an interior interface
std::unique_ptr<Mesh> CreateTestMesh(int nx = 4, int nz = 8)
{
   // Create a rectangular mesh
   // Full domain from -Lx to +Lx
   real_t Lx = 50.0e3;   // 50 km half-width
   real_t Lz = 50.0e3;   // 50 km depth

   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian2D(2 * nx, nz, Element::QUADRILATERAL,
                            true, 2.0 * Lx, Lz));

   // Shift mesh so that:
   // x: [-Lx, +Lx]
   // z: [-Lz, 0]
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *coords = mesh->GetVertex(i);
      coords[0] -= Lx;   // x: [0, 2Lx] -> [-Lx, +Lx]
      coords[1] -= Lz;   // z: [0, Lz] -> [-Lz, 0]
   }

   mesh->FinalizeTopology();
   mesh->Finalize();

   return mesh;
}

// =============================================================================
// BP2 Pre-stress Tests
// =============================================================================

void TestPreStress()
{
   std::cout << "\n=== Testing Pre-stress Computation ===\n";

   BP2Params params;

   // Test 1: Pre-stress tau0 matches BP2 specification (~26.546 MPa)
   {
      real_t tau0 = params.tau0();
      real_t tau0_MPa = tau0 / 1.0e6;

      // Expected: ~26.546 MPa (from BP2 benchmark)
      TEST_NEAR(tau0_MPa, 26.546, 0.01, "Pre-stress tau0 ~ 26.546 MPa");
      std::cout << "    tau0 = " << tau0_MPa << " MPa\n";
   }

   // Test 2: Radiation damping coefficient
   {
      real_t eta = params.eta();
      real_t eta_MPa_s_m = eta / 1.0e6;

      // Expected: eta = sqrt(mu * rho) / 2 ~ 4.63 MPa*s/m
      TEST_NEAR(eta_MPa_s_m, 4.63, 0.01, "Radiation damping eta ~ 4.63 MPa·s/m");
      std::cout << "    eta = " << eta_MPa_s_m << " MPa·s/m\n";
   }

   // Test 3: Shear modulus
   {
      real_t mu = params.mu();
      real_t mu_GPa = mu / 1.0e9;

      // Expected: mu = rho * cs^2 ~ 32.04 GPa
      TEST_NEAR(mu_GPa, 32.04, 0.01, "Shear modulus mu ~ 32.04 GPa");
      std::cout << "    mu = " << mu_GPa << " GPa\n";
   }
}

// =============================================================================
// a(z) Depth Profile Tests
// =============================================================================

void TestAofZ()
{
   std::cout << "\n=== Testing a(z) Depth Profile ===\n";

   BP2Params params;

   // Test 1: a at surface (z = 0, in VW zone)
   {
      real_t a = params.a_of_z(0.0);
      TEST_NEAR(a, params.a0, 1e-10, "a(z=0) = a0 = 0.010 (VW zone)");
   }

   // Test 2: a at depth = H (boundary of VW zone)
   {
      real_t z = -params.H;  // -15 km
      real_t a = params.a_of_z(z);
      TEST_NEAR(a, params.a0, 1e-10, "a(z=-H) = a0 = 0.010 (top of transition)");
   }

   // Test 3: a at depth = H + h/2 (middle of transition)
   {
      real_t z = -(params.H + params.h / 2);  // -16.5 km
      real_t a = params.a_of_z(z);
      real_t expected = params.a0 + (params.amax - params.a0) / 2.0;
      TEST_NEAR(a, expected, 1e-10, "a(z=-16.5km) = (a0+amax)/2 (middle of transition)");
   }

   // Test 4: a at depth = H + h (bottom of transition, fully VS)
   {
      real_t z = -(params.H + params.h);  // -18 km
      real_t a = params.a_of_z(z);
      TEST_NEAR(a, params.amax, 1e-10, "a(z=-18km) = amax = 0.025 (VS zone)");
   }

   // Test 5: a at depth = 30 km (deep VS)
   {
      real_t z = -30.0e3;
      real_t a = params.a_of_z(z);
      TEST_NEAR(a, params.amax, 1e-10, "a(z=-30km) = amax = 0.025 (deep VS)");
   }

   // Test 6: VW zone identification
   {
      TEST_ASSERT(params.IsVelocityWeakening(0.0),
                  "z=0 is velocity-weakening (a < b)");
      TEST_ASSERT(params.IsVelocityWeakening(-10.0e3),
                  "z=-10km is velocity-weakening (a < b)");
      TEST_ASSERT(!params.IsVelocityWeakening(-30.0e3),
                  "z=-30km is velocity-strengthening (a > b)");
   }
}

// =============================================================================
// Fault Geometry Tests (with domain operator)
// =============================================================================

void TestFaultGeometry()
{
   std::cout << "\n=== Testing Fault Geometry ===\n";

   BP2Params params;
   auto mesh = CreateTestMesh(4, 8);

   // Create domain operator
   AntiplaneDomainOperator<Mesh> domain(*mesh, 1, params.mu(),
                                         params.Vp, params.Wf);

   // Create fault geometry
   FaultGeometry<Mesh> fault_geom(domain, params);

   // Test 1: Number of fault DOFs
   {
      int num_dofs = fault_geom.NumFaultDOFs();
      TEST_ASSERT(num_dofs > 0, "Fault has positive number of DOFs");
      std::cout << "    num_fault_dofs = " << num_dofs << "\n";
   }

   // Test 2: Depths are in expected range
   {
      const Vector &depths = fault_geom.GetDepths();
      real_t z_min = depths.Min();
      real_t z_max = depths.Max();

      // Depths should be between -Lz and 0.
      // The fault interface at x=0 extends the full mesh depth.
      // DOFs below -Wf get prescribed plate rate Vp (not rate-state).
      real_t Lz = 50.0e3;  // matches CreateTestMesh
      TEST_ASSERT(z_max <= 0.0, "Max depth <= 0 (at or below surface)");
      TEST_ASSERT(z_min >= -Lz, "Min depth >= -Lz (within mesh domain)");

      std::cout << "    Depth range: [" << z_max/1000.0 << ", "
                << z_min/1000.0 << "] km\n";
   }

   // Test 3: a values match depth profile
   {
      const Vector &depths = fault_geom.GetDepths();
      const Vector &a_values = fault_geom.GetAValues();

      bool all_match = true;
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         real_t z = depths(i);
         real_t a_expected = params.a_of_z(z);
         if (std::abs(a_values(i) - a_expected) > 1e-12)
         {
            all_match = false;
            break;
         }
      }

      TEST_ASSERT(all_match, "All a values match a(z) profile");
   }

   // Test 4: eta values are constant (for BP2)
   {
      const Vector &eta_values = fault_geom.GetEtaValues();
      real_t eta_expected = params.eta();

      bool all_match = true;
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         if (std::abs(eta_values(i) - eta_expected) > 1e-12)
         {
            all_match = false;
            break;
         }
      }

      TEST_ASSERT(all_match, "All eta values are constant");
   }

   // Test 5: FindNearestDOF works
   {
      // Find DOF closest to z = -10 km
      int idx = fault_geom.FindNearestDOF(-10.0e3);
      TEST_ASSERT(idx >= 0 && idx < fault_geom.NumFaultDOFs(),
                  "FindNearestDOF returns valid index");
   }

   // Print geometry info
   std::cout << "\n    Fault Geometry Info:\n";
   fault_geom.Print(std::cout);
}

// =============================================================================
// Rate-State Fault Operator Tests
// =============================================================================

void TestRateStateFaultOperator()
{
   std::cout << "\n=== Testing Rate-State Fault Operator ===\n";

   BP2Params params;
   auto mesh = CreateTestMesh(4, 8);

   // Create domain operator
   AntiplaneDomainOperator<Mesh> domain(*mesh, 1, params.mu(),
                                         params.Vp, params.Wf);

   // Create fault geometry
   FaultGeometry<Mesh> fault_geom(domain, params);

   // Create friction law and state evolution
   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   // Create fault operator
   RateStateFaultOperator<Mesh> fault_op(&fault_geom, &friction, &aging, params);

   // Test 1: State size
   {
      int state_size = fault_op.StateSize();
      int expected = fault_geom.NumFaultDOFs() * 2;  // slip + theta per node
      TEST_ASSERT(state_size == expected, "State size = 2 * num_nodes");
   }

   // Test 2: PreInit
   {
      Vector state(fault_op.StateSize());
      fault_op.PreInit(state);

      // Check that all slips are zero
      Vector slip;
      fault_op.GetSlip(state, slip);
      TEST_NEAR(slip.Norml2(), 0.0, 1e-15, "PreInit: all slips are zero");

      // Check that theta values are positive (placeholder)
      Vector theta;
      fault_op.GetTheta(state, theta);
      TEST_ASSERT(theta.Min() > 0.0, "PreInit: all theta values positive");
   }

   // Test 3: Init with zero traction
   {
      Vector state(fault_op.StateSize());
      fault_op.PreInit(state);

      // Zero traction (no quasi-static contribution)
      Vector traction(fault_geom.NumFaultDOFs());
      traction = 0.0;

      real_t V_max = fault_op.Init(traction, state);

      // Check tau0
      TEST_NEAR(fault_op.GetTau0() / 1e6, 26.546, 0.01,
                "Init: tau0 ~ 26.546 MPa");

      // Check V_max (should be close to V_init at initialization)
      TEST_ASSERT(V_max > 0.0, "Init: V_max > 0");

      // Check theta values are reasonable
      Vector theta;
      fault_op.GetTheta(state, theta);
      TEST_ASSERT(theta.Min() > 100.0, "Init: theta values > 100 s");
      TEST_ASSERT(theta.Max() < 1e9, "Init: theta values < 1e9 s");

      std::cout << "    V_max = " << V_max << " m/s\n";
      std::cout << "    theta range: [" << theta.Min() << ", " << theta.Max() << "] s\n";
   }

   // Test 4: Stress equilibrium after Init
   {
      Vector state(fault_op.StateSize());
      fault_op.PreInit(state);

      Vector traction(fault_geom.NumFaultDOFs());
      traction = 0.0;

      fault_op.Init(traction, state);

      real_t max_error = fault_op.VerifyStressEquilibrium(traction, state);
      TEST_ASSERT(max_error < 1e-8, "Stress equilibrium satisfied after Init");
      std::cout << "    Max stress balance error: " << max_error << "\n";
   }

   // Test 5: ComputeRHS
   {
      Vector state(fault_op.StateSize());
      fault_op.PreInit(state);

      Vector traction(fault_geom.NumFaultDOFs());
      traction = 0.0;

      fault_op.Init(traction, state);

      Vector rate(fault_op.StateSize());
      real_t V_max = fault_op.ComputeRHS(traction, state, rate);

      // Check that slip rates are positive
      Vector slip_rate;
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         real_t V = rate(i * 2 + 0);  // SlipIndex = 0
         TEST_ASSERT(V > 0.0, "ComputeRHS: slip rate > 0");
         if (V <= 0.0) break;  // Don't flood with errors
      }

      // Check that state rates have correct signs
      // At initial state close to steady state, dtheta/dt ~ 0
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         real_t dtheta_dt = rate(i * 2 + 1);  // ThetaIndex = 1
         // Shouldn't be extremely large
         TEST_ASSERT(std::abs(dtheta_dt) < 1e10,
                     "ComputeRHS: dtheta/dt reasonable magnitude");
         if (std::abs(dtheta_dt) >= 1e10) break;
      }

      std::cout << "    V_max from RHS: " << V_max << " m/s\n";
   }

   // Test 6: State access methods roundtrip
   {
      Vector state(fault_op.StateSize());
      fault_op.PreInit(state);

      // Create some test data
      Vector slip_in(fault_geom.NumFaultDOFs());
      Vector theta_in(fault_geom.NumFaultDOFs());
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         slip_in(i) = i * 0.001;
         theta_in(i) = 1000.0 + i * 100.0;
      }

      // Set and get back
      fault_op.SetSlip(slip_in, state);
      fault_op.SetTheta(theta_in, state);

      Vector slip_out, theta_out;
      fault_op.GetSlip(state, slip_out);
      fault_op.GetTheta(state, theta_out);

      real_t slip_diff = 0.0, theta_diff = 0.0;
      for (int i = 0; i < fault_geom.NumFaultDOFs(); i++)
      {
         slip_diff = std::max(slip_diff, std::abs(slip_in(i) - slip_out(i)));
         theta_diff = std::max(theta_diff, std::abs(theta_in(i) - theta_out(i)));
      }

      TEST_NEAR(slip_diff, 0.0, 1e-15, "State access roundtrip: slip");
      TEST_NEAR(theta_diff, 0.0, 1e-15, "State access roundtrip: theta");
   }
}

// =============================================================================
// Initial State Verification
// =============================================================================

void TestInitialStateValues()
{
   std::cout << "\n=== Testing Initial State Values ===\n";

   BP2Params params;

   // Create friction law
   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);

   real_t sigma_n = params.sigma_n;
   real_t eta = params.eta();
   real_t V_init = params.V_init;
   real_t tau0 = params.tau0();

   // Test 1: Initial theta at surface (a = a0)
   {
      real_t a = params.a0;
      real_t theta0 = friction.InitialState(tau0, V_init, sigma_n, eta, a);

      // Verify stress balance
      real_t f = friction.FrictionCoefficient(V_init, theta0, a);
      real_t tau_check = sigma_n * f + eta * V_init;

      TEST_REL_NEAR(tau_check, tau0, 1e-10,
                    "Initial state (a=a0): stress equilibrium");
      std::cout << "    theta0(z=0) = " << theta0 << " s"
                << " (log10: " << std::log10(theta0) << ")\n";
   }

   // Test 2: Initial theta at VS zone (a = amax)
   {
      real_t a = params.amax;
      real_t theta0 = friction.InitialState(tau0, V_init, sigma_n, eta, a);

      // Verify stress balance
      real_t f = friction.FrictionCoefficient(V_init, theta0, a);
      real_t tau_check = sigma_n * f + eta * V_init;

      TEST_REL_NEAR(tau_check, tau0, 1e-10,
                    "Initial state (a=amax): stress equilibrium");
      std::cout << "    theta0(VS zone) = " << theta0 << " s"
                << " (log10: " << std::log10(theta0) << ")\n";
   }

   // Test 3: Slip rate recovery at initial state
   {
      real_t a = params.a0;
      real_t theta0 = friction.InitialState(tau0, V_init, sigma_n, eta, a);
      real_t V_solved = friction.SolveSlipRate(tau0, theta0, sigma_n, eta, a);

      TEST_REL_NEAR(V_solved, V_init, 1e-6,
                    "Initial state: solved V matches V_init");
   }
}

// =============================================================================
// State Evolution Consistency
// =============================================================================

void TestStateEvolutionConsistency()
{
   std::cout << "\n=== Testing State Evolution Consistency ===\n";

   BP2Params params;
   AgingLaw aging;

   // Test 1: At steady state (theta = Dc/V), dtheta/dt = 0
   {
      real_t V = 1e-9;
      real_t theta_ss = params.Dc / V;
      real_t rate = aging.Rate(V, theta_ss, params.Dc);

      TEST_NEAR(rate, 0.0, 1e-15, "At steady state: dtheta/dt = 0");
   }

   // Test 2: When locked (V ~ 0), dtheta/dt = 1 (healing)
   {
      real_t V = 1e-20;  // Nearly locked
      real_t theta = 1000.0;
      real_t rate = aging.Rate(V, theta, params.Dc);

      TEST_NEAR(rate, 1.0, 1e-10, "When locked: dtheta/dt ~ 1 (healing)");
   }

   // Test 3: During fast slip, theta decreases (weakening)
   {
      real_t V = 1.0;  // Fast slip
      real_t theta = 4000.0;
      real_t rate = aging.Rate(V, theta, params.Dc);

      TEST_ASSERT(rate < 0.0, "During fast slip: dtheta/dt < 0 (weakening)");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "========================================\n";
   std::cout << "SEAS Phase 3 Unit Tests: Fault Operator\n";
   std::cout << "========================================\n";

   TestPreStress();
   TestAofZ();
   TestFaultGeometry();
   TestRateStateFaultOperator();
   TestInitialStateValues();
   TestStateEvolutionConsistency();

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
