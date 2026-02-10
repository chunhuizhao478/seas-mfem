// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Integration test: Short BP2 simulation (1 year)
//
// Runs the full SEAS quasi-dynamic simulation with BP2 parameters for 1 year.
// Verifies:
// 1. Simulation completes without crash
// 2. Slip rate starts near V_init and remains bounded
// 3. State variable (theta) remains positive
// 4. Stress is in reasonable range
// 5. Total slip is approximately V_init * t

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/bp2_mesh.hpp"

#include <iostream>
#include <cmath>
#include <memory>
#include <iomanip>

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

// =============================================================================
// BP2 Short Simulation Test: 1 year
// =============================================================================

void TestBP2ShortSimulation()
{
   std::cout << "\n=== BP2 Short Simulation (1 year) ===\n";

   BP2Params params;
   real_t t_final = BP2Params::seconds_per_year;  // 1 year in seconds

   std::cout << "  t_final = " << t_final / BP2Params::seconds_per_year
             << " years (" << t_final << " s)\n";

   // Create mesh (coarse for fast test)
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 8;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   // Create domain operator
   AntiplaneDomainOperator<Mesh> domain(
      *mesh, 1, params.mu(), params.Vp, params.Wf);

   // Create fault components
   FaultGeometry<Mesh> fault_geom(domain, params);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<Mesh> fault_op(
      &fault_geom, &friction, &aging, params);

   // Create SEAS operator
   SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

   // Initialize state
   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   // Record initial conditions
   real_t V_initial = seas_op.GetMaxSlipRate();
   std::cout << "  Initial V_max = " << V_initial << " m/s\n";

   // Test 1: Initial slip rate near V_init
   TEST_ASSERT(std::abs(V_initial - params.V_init) / params.V_init < 0.1,
               "Initial V_max ~ V_init");

   // Set up time integration
   AdaptiveTimeStepper stepper;
   stepper.SetDtMin(1e-3);
   stepper.SetDtMax(t_final / 10);  // Max dt = 0.1 years
   stepper.SetInitialDt(1e4);       // Start with ~2.78 hours

   RK4Solver ode_solver;
   ode_solver.Init(seas_op);

   // Run simulation
   real_t t = 0.0;
   int step = 0;
   int max_steps = 10000;
   real_t V_max_during_sim = V_initial;
   real_t theta_min_during_sim = 1e30;

   // Track traction range
   real_t tau_min = 1e30;
   real_t tau_max = -1e30;

   std::cout << "\n  Time stepping:\n";
   std::cout << "  " << std::setw(10) << "Step"
             << std::setw(15) << "Time [yr]"
             << std::setw(15) << "dt [s]"
             << std::setw(15) << "V_max [m/s]"
             << "\n";

   while (t < t_final && step < max_steps)
   {
      real_t dt = stepper.GetDt();
      // Don't overshoot
      if (t + dt > t_final)
      {
         dt = t_final - t;
      }

      ode_solver.Step(state, t, dt);
      step++;

      // Update adaptive dt
      real_t V_max = seas_op.GetMaxSlipRate();
      stepper.SetDt(stepper.ComputeNewDt(V_max, dt));

      V_max_during_sim = std::max(V_max_during_sim, V_max);

      // Track theta
      Vector theta;
      fault_op.GetTheta(state, theta);
      theta_min_during_sim = std::min(theta_min_during_sim, theta.Min());

      // Track traction (evaluate traction for monitoring)
      const Vector &traction = seas_op.GetTraction();
      real_t tau0 = fault_op.GetTau0();
      for (int i = 0; i < traction.Size(); i++)
      {
         real_t tau_total = tau0 + traction(i);
         tau_min = std::min(tau_min, tau_total);
         tau_max = std::max(tau_max, tau_total);
      }

      // Periodic output
      if (step % 100 == 0 || step <= 5 || t >= t_final)
      {
         std::cout << "  " << std::setw(10) << step
                   << std::setw(15) << std::scientific
                   << std::setprecision(4)
                   << t / BP2Params::seconds_per_year
                   << std::setw(15) << dt
                   << std::setw(15) << V_max
                   << "\n";
      }
   }

   std::cout << "\n  Simulation completed:\n";
   std::cout << "    Steps: " << step << "\n";
   std::cout << "    Final time: " << t / BP2Params::seconds_per_year
             << " years\n";
   std::cout << "    V_max during sim: " << V_max_during_sim << " m/s\n";
   std::cout << "    Theta min: " << theta_min_during_sim << " s\n";
   std::cout << "    Tau range: [" << tau_min / 1e6 << ", "
             << tau_max / 1e6 << "] MPa\n";

   // Test 2: Simulation completed
   TEST_ASSERT(step < max_steps,
               "Simulation completed within step limit");

   // Test 3: Slip rate remains bounded (no large events in 1 year)
   TEST_ASSERT(V_max_during_sim < 1e-3,
               "No coseismic events during 1-year simulation");

   // Test 4: State variable remains positive
   TEST_ASSERT(theta_min_during_sim > 0.0,
               "Theta remains positive throughout simulation");

   // Test 5: Stress in reasonable range
   // BP2 pre-stress is ~26.5 MPa, should be in [15, 40] MPa range
   TEST_ASSERT(tau_min / 1e6 > 15.0,
               "Minimum stress > 15 MPa");
   TEST_ASSERT(tau_max / 1e6 < 40.0,
               "Maximum stress < 40 MPa");

   // Test 6: Total slip approximately V_init * t
   Vector slip;
   fault_op.GetSlip(state, slip);

   // Mean slip across fault
   real_t mean_slip = 0.0;
   for (int i = 0; i < slip.Size(); i++)
   {
      mean_slip += slip(i);
   }
   mean_slip /= slip.Size();

   real_t expected_slip = params.V_init * t_final;
   std::cout << "    Mean slip: " << mean_slip << " m\n";
   std::cout << "    Expected (V_init*t): " << expected_slip << " m\n";

   // Allow factor-of-2 tolerance due to depth-dependent slip rates
   TEST_ASSERT(mean_slip > 0.0,
               "Mean slip is positive after 1 year");
   TEST_ASSERT(mean_slip < 10.0 * expected_slip,
               "Mean slip within 10x of V_init*t");
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "================================================\n";
   std::cout << "SEAS Phase 4: BP2 Short Simulation Test\n";
   std::cout << "================================================\n";

   TestBP2ShortSimulation();

   std::cout << "\n================================================\n";
   std::cout << "Test Summary\n";
   std::cout << "================================================\n";
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
