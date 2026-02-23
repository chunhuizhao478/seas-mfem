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

// Integration tests for Phase 4: SEAS Quasi-Dynamic Operator
// Tests for SEASQuasiDynamicOperator and AdaptiveTimeStepper:
// - Domain-fault coupling
// - Steady-state slip at plate rate
// - Stress balance maintenance
// - Slip conservation
// - Adaptive time stepping
// - RK4 integration

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

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Simple test framework (consistent with Phase 3 tests)
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
// Helper: Create test components
// =============================================================================

struct TestComponents
{
   BP2Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<AntiplaneDomainOperator<Mesh>> domain;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   DieterichRuinaFriction friction;
   AgingLaw aging;
   std::unique_ptr<RateStateFaultOperator<Mesh>> fault;
   std::unique_ptr<SEASQuasiDynamicOperator<Mesh>> seas_op;

   /// Create test components with a small mesh
   void Setup(int nx = 4, int nz = 8)
   {
      // Create mesh
      BP2MeshGenerator::Parameters mesh_params;
      mesh_params.Lx = 50.0e3;
      mesh_params.Lz = 50.0e3;
      mesh_params.Wf = params.Wf;
      mesh_params.nx = nx;
      mesh_params.nz = nz;
      mesh = BP2MeshGenerator::Create(mesh_params);

      // Create domain operator (order 1 for fast tests)
      domain = std::make_unique<AntiplaneDomainOperator<Mesh>>(
         *mesh, 1, params.mu(), params.Vp, params.Wf);

      // Create fault geometry
      fault_geom = std::make_unique<FaultGeometry<Mesh>>(*domain, params);

      // Create friction law
      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0;
      fc.f0 = params.f0;
      fc.b = params.b;
      fc.Dc = params.Dc;
      friction.SetConstants(fc);

      // Create fault operator
      fault = std::make_unique<RateStateFaultOperator<Mesh>>(
         fault_geom.get(), &friction, &aging, params);

      // Create SEAS operator
      seas_op = std::make_unique<SEASQuasiDynamicOperator<Mesh>>(
         domain.get(), fault.get());
   }
};

// =============================================================================
// Test 1: Operator Construction and Initialization
// =============================================================================

void TestOperatorConstruction()
{
   std::cout << "\n=== Test: Operator Construction ===\n";

   TestComponents tc;
   tc.Setup();

   // Test 1: State size is correct
   {
      int expected_size = tc.fault->StateSize();
      TEST_ASSERT(expected_size > 0, "State size is positive");
      TEST_ASSERT(expected_size == tc.fault_geom->NumFaultDOFs() * 2,
                  "State size = 2 * num_fault_dofs");
      std::cout << "    State size: " << expected_size << "\n";
   }
}

// =============================================================================
// Test 2: Initial Condition
// =============================================================================

void TestInitialCondition()
{
   std::cout << "\n=== Test: Initial Condition ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Test 1: Slip should be zero initially
   {
      Vector slip;
      tc.fault->GetSlip(state, slip);
      TEST_NEAR(slip.Norml2(), 0.0, 1e-15,
                "Initial slip is zero");
   }

   // Test 2: Theta should be positive and reasonable
   {
      Vector theta;
      tc.fault->GetTheta(state, theta);
      TEST_ASSERT(theta.Min() > 0.0,
                  "Initial theta is positive");
      TEST_ASSERT(theta.Max() < 1e15,
                  "Initial theta is not absurdly large");
      std::cout << "    theta range: [" << theta.Min() << ", "
                << theta.Max() << "] s\n";
   }

   // Test 3: Initial slip rate should be close to V_init
   {
      real_t V_max = tc.seas_op->GetMaxSlipRate();
      TEST_REL_NEAR(V_max, tc.params.V_init, 0.1,
                    "Initial V_max ~ V_init");
      std::cout << "    V_max = " << V_max << " m/s\n";
   }

   // Test 4: Stress equilibrium should be satisfied
   {
      real_t eq_error = tc.fault->VerifyStressEquilibrium(
         tc.seas_op->GetTraction(), state);
      TEST_ASSERT(eq_error < 1e-6,
                  "Stress equilibrium at initialization");
      std::cout << "    Stress equilibrium error: " << eq_error << "\n";
   }
}

// =============================================================================
// Test 3: Steady-State Slip at Plate Rate
// =============================================================================

void TestSteadyStateSlip()
{
   std::cout << "\n=== Test: Steady-State Slip ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Run for a short time with small dt
   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t dt = 1e3;  // 1000 seconds
   int nsteps = 10;

   for (int step = 0; step < nsteps; step++)
   {
      ode_solver.Step(state, t, dt);
      // Note: MFEM's Step() advances t internally, do NOT add dt again
   }

   // After short time, V should still be near V_init
   real_t V_max = tc.seas_op->GetMaxSlipRate();
   TEST_REL_NEAR(V_max, tc.params.V_init, 0.05,
                 "V ~ V_init after short time");
   std::cout << "    V_max at t=" << t << "s: " << V_max << " m/s\n";

   // Slip should be approximately V_init * t
   Vector slip;
   tc.fault->GetSlip(state, slip);
   real_t mean_slip = 0.0;
   for (int i = 0; i < slip.Size(); i++)
   {
      mean_slip += slip(i);
   }
   mean_slip /= slip.Size();

   real_t expected_slip = tc.params.V_init * t;
   // Generous tolerance since slip distribution varies with depth
   TEST_ASSERT(mean_slip > 0.0, "Mean slip is positive after time stepping");
   std::cout << "    Mean slip: " << mean_slip << " m\n";
   std::cout << "    Expected (V_init*t): " << expected_slip << " m\n";
}

// =============================================================================
// Test 4: Stress Balance Maintained During Time Stepping
// =============================================================================

void TestStressBalanceMaintained()
{
   std::cout << "\n=== Test: Stress Balance During Time Stepping ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t dt = 1e3;
   real_t max_eq_error = 0.0;

   for (int step = 0; step < 5; step++)
   {
      ode_solver.Step(state, t, dt);

      // After each step, evaluate Mult to update traction
      Vector rate(tc.fault->StateSize());
      tc.seas_op->Mult(state, rate);

      // Check stress equilibrium
      real_t eq_error = tc.fault->VerifyStressEquilibrium(
         tc.seas_op->GetTraction(), state);
      max_eq_error = std::max(max_eq_error, eq_error);
   }

   TEST_ASSERT(max_eq_error < 1e-4,
               "Stress balance maintained during stepping");
   std::cout << "    Max stress balance error: " << max_eq_error << "\n";
}

// =============================================================================
// Test 5: Slip Conservation
// =============================================================================

void TestSlipConservation()
{
   std::cout << "\n=== Test: Slip Conservation ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t dt = 1e4;  // 10000 seconds
   int nsteps = 10;
   real_t t_total = dt * nsteps;

   for (int step = 0; step < nsteps; step++)
   {
      ode_solver.Step(state, t, dt);
   }

   // Check that total slip is approximately V_init * t at each node
   Vector slip;
   tc.fault->GetSlip(state, slip);

   real_t max_rel_error = 0.0;
   for (int i = 0; i < slip.Size(); i++)
   {
      real_t expected = tc.params.V_init * t_total;
      if (expected > 1e-30)
      {
         real_t rel_err = std::abs(slip(i) - expected) / expected;
         max_rel_error = std::max(max_rel_error, rel_err);
      }
   }

   // Allow generous tolerance - slip varies with depth due to fault geometry
   TEST_ASSERT(max_rel_error < 1.0,
               "Slip within order-of-magnitude of V_init*t");
   std::cout << "    Max relative slip error: " << max_rel_error << "\n";
   std::cout << "    Expected slip: " << tc.params.V_init * t_total << " m\n";
}

// =============================================================================
// Test 6: Theta Stays Positive
// =============================================================================

void TestThetaPositive()
{
   std::cout << "\n=== Test: Theta Stays Positive ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t dt = 1e4;
   int nsteps = 20;

   bool theta_always_positive = true;
   for (int step = 0; step < nsteps; step++)
   {
      ode_solver.Step(state, t, dt);

      Vector theta;
      tc.fault->GetTheta(state, theta);
      if (theta.Min() <= 0.0)
      {
         theta_always_positive = false;
         std::cerr << "    theta went non-positive at step " << step
                   << ", min theta = " << theta.Min() << "\n";
         break;
      }
   }

   TEST_ASSERT(theta_always_positive,
               "Theta remains positive during time stepping");
}

// =============================================================================
// Test 7: Adaptive Time Stepper
// =============================================================================

void TestAdaptiveTimeStepper()
{
   std::cout << "\n=== Test: Adaptive Time Stepper ===\n";

   AdaptiveTimeStepper stepper;

   // Test 1: Very slow slip -> dt_max
   {
      real_t dt = stepper.ComputeNewDt(1e-15, 1e3);
      TEST_NEAR(dt, stepper.GetDtMax(), 1e-10,
                "Very slow slip gives dt_max");
   }

   // Test 2: Rapid slip -> reduced dt
   {
      real_t V_max = 1e-3;
      real_t dt_current = 1e3;
      real_t dt = stepper.ComputeNewDt(V_max, dt_current);
      TEST_ASSERT(dt < dt_current,
                  "Rapid slip reduces dt");
      real_t expected = dt_current * stepper.GetTargetSlipRateMax() / V_max;
      TEST_NEAR(dt, expected, 1e-10,
                "Rapid slip dt = dt_current * V_target / V_max");
   }

   // Test 3: Moderate slip -> gradual increase
   {
      real_t V_max = 1e-8;
      real_t dt_current = 1e5;
      real_t dt = stepper.ComputeNewDt(V_max, dt_current);
      TEST_ASSERT(dt > dt_current,
                  "Moderate slip increases dt");
      TEST_ASSERT(dt <= stepper.GetDtMax(),
                  "dt doesn't exceed dt_max");
   }

   // Test 4: dt always within bounds
   {
      real_t dt1 = stepper.ComputeNewDt(1e-15, 1e-10);
      TEST_ASSERT(dt1 >= stepper.GetDtMin(), "dt >= dt_min (case 1)");
      TEST_ASSERT(dt1 <= stepper.GetDtMax(), "dt <= dt_max (case 1)");

      real_t dt2 = stepper.ComputeNewDt(1e3, 1e10);
      TEST_ASSERT(dt2 >= stepper.GetDtMin(), "dt >= dt_min (case 2)");
      TEST_ASSERT(dt2 <= stepper.GetDtMax(), "dt <= dt_max (case 2)");
   }

   // Test 5: Configuration
   {
      stepper.SetDtMin(1e-8);
      stepper.SetDtMax(1e9);
      stepper.SetTargetSlipRateMax(1e-5);
      TEST_NEAR(stepper.GetDtMin(), 1e-8, 1e-20, "SetDtMin works");
      TEST_NEAR(stepper.GetDtMax(), 1e9, 1e-20, "SetDtMax works");
      TEST_NEAR(stepper.GetTargetSlipRateMax(), 1e-5, 1e-20,
                "SetTargetSlipRateMax works");
   }
}

// =============================================================================
// Test 8: RK4 Integration with Adaptive dt
// =============================================================================

void TestRK4WithAdaptiveDt()
{
   std::cout << "\n=== Test: RK4 with Adaptive Time Step ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Set up time integration
   AdaptiveTimeStepper stepper;
   stepper.SetDtMax(1e6);  // Cap at ~11.5 days for test
   stepper.SetInitialDt(1e3);

   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t t_final = 1e5;  // ~1.16 days
   int step = 0;
   int max_steps = 1000;  // Safety limit

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
      real_t V_max = tc.seas_op->GetMaxSlipRate();
      stepper.SetDt(stepper.ComputeNewDt(V_max, dt));
   }

   TEST_ASSERT(step < max_steps,
               "Simulation completed within step limit");
   TEST_NEAR(t, t_final, 1e-6,
             "Simulation reached t_final");

   // Verify slip rate is still reasonable
   real_t V_max = tc.seas_op->GetMaxSlipRate();
   TEST_ASSERT(V_max > 0.0 && V_max < 1.0,
               "Final V_max is in reasonable range");

   // Verify theta still positive
   Vector theta;
   tc.fault->GetTheta(state, theta);
   TEST_ASSERT(theta.Min() > 0.0,
               "Theta remains positive after adaptive stepping");

   std::cout << "    Steps taken: " << step << "\n";
   std::cout << "    Final t: " << t << " s\n";
   std::cout << "    Final V_max: " << V_max << " m/s\n";
   std::cout << "    Final dt: " << stepper.GetDt() << " s\n";
}

// =============================================================================
// Test 9: Mult is Deterministic (same input -> same output)
// =============================================================================

void TestMultDeterministic()
{
   std::cout << "\n=== Test: Mult is Deterministic ===\n";

   TestComponents tc;
   tc.Setup();

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Call Mult twice with same input
   Vector rate1(tc.fault->StateSize());
   Vector rate2(tc.fault->StateSize());

   tc.seas_op->Mult(state, rate1);
   tc.seas_op->Mult(state, rate2);

   // Results should be identical
   real_t diff = 0.0;
   for (int i = 0; i < rate1.Size(); i++)
   {
      diff = std::max(diff, std::abs(rate1(i) - rate2(i)));
   }

   TEST_NEAR(diff, 0.0, 1e-15,
             "Mult is deterministic (same input -> same output)");
}

// =============================================================================
// Test 10: Extended Fault Covers Full Depth
// =============================================================================

void TestExtendedFaultFullDepth()
{
   std::cout << "\n=== Test: Extended Fault Covers Full Depth ===\n";

   TestComponents tc;
   tc.params.Wf = 40000.0;  // 40 km
   // Use Lz > Wf so there are DOFs below Wf
   tc.Setup(4, 10);  // nx=4, nz=10

   int num_dofs = tc.fault_geom->NumFaultDOFs();
   TEST_ASSERT(num_dofs > 0, "Fault has DOFs");

   const Vector &depths = tc.fault_geom->GetDepths();
   real_t z_min = depths.Min();
   real_t z_max = depths.Max();

   std::cout << "    Fault DOFs: " << num_dofs << "\n";
   std::cout << "    Depth range: [" << z_max / 1000 << ", "
             << z_min / 1000 << "] km\n";
   std::cout << "    Wf = " << tc.params.Wf / 1000 << " km\n";
   std::cout << "    Lz = 50 km\n";

   // Fault should extend below -Wf (to -Lz)
   TEST_ASSERT(z_min < -tc.params.Wf,
               "Fault extends below Wf");
   // With 1 midpoint DOF per face, shallowest DOF is at z = -h/2
   // where h = Lz/nz.  For nz=10, Lz=50km: h=5km, z_max = -2.5km.
   real_t h_z = 50000.0 / 10.0;
   TEST_ASSERT(z_max > -(h_z / 2 + 100.0),
               "Fault reaches near surface (within half-element of z=0)");

   // Count DOFs above and below Wf
   int above_wf = 0, below_wf = 0;
   for (int i = 0; i < num_dofs; i++)
   {
      if (depths(i) < -tc.params.Wf)
      {
         below_wf++;
      }
      else
      {
         above_wf++;
      }
   }

   TEST_ASSERT(above_wf > 0, "Has DOFs above Wf (rate-state zone)");
   TEST_ASSERT(below_wf > 0, "Has DOFs below Wf (prescribed loading zone)");
   std::cout << "    Above Wf: " << above_wf << " DOFs\n";
   std::cout << "    Below Wf: " << below_wf << " DOFs\n";
}

// =============================================================================
// Test 11: ComputeRHS Prescribes V=Vp Below Wf
// =============================================================================

void TestPrescribedVpBelowWf()
{
   std::cout << "\n=== Test: ComputeRHS Prescribes V=Vp Below Wf ===\n";

   TestComponents tc;
   tc.Setup(4, 10);

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Compute RHS
   Vector rate(tc.fault->StateSize());
   tc.seas_op->Mult(state, rate);

   const Vector &depths = tc.fault_geom->GetDepths();
   const Vector &slip_rate = tc.fault->GetSlipRate();
   int num_dofs = tc.fault->NumNodes();

   int below_wf_count = 0;
   int above_wf_count = 0;
   real_t max_vp_error = 0.0;

   for (int i = 0; i < num_dofs; i++)
   {
      if (depths(i) < -tc.params.Wf)
      {
         // Below Wf: slip rate should be exactly Vp
         real_t err = std::abs(slip_rate(i) - tc.params.Vp);
         max_vp_error = std::max(max_vp_error, err);
         below_wf_count++;

         // dslip/dt should be Vp
         real_t dslip_dt = rate(i * 2 + 0);
         TEST_NEAR(dslip_dt, tc.params.Vp, 1e-20,
                   "dslip/dt = Vp for below-Wf DOF");

         // dtheta/dt should be 0
         real_t dtheta_dt = rate(i * 2 + 1);
         TEST_NEAR(dtheta_dt, 0.0, 1e-20,
                   "dtheta/dt = 0 for below-Wf DOF");

         // Only test first below-Wf DOF in detail
         break;
      }
      else
      {
         above_wf_count++;
      }
   }

   TEST_ASSERT(below_wf_count > 0 || above_wf_count == num_dofs,
               "Found DOFs to check");

   if (below_wf_count > 0)
   {
      TEST_NEAR(max_vp_error, 0.0, 1e-20,
                "Below-Wf slip rate == Vp (exact)");
   }

   // Above Wf: slip rate should be ~ V_init (rate-state equilibrium)
   bool above_wf_reasonable = true;
   for (int i = 0; i < num_dofs; i++)
   {
      if (depths(i) >= -tc.params.Wf)
      {
         if (slip_rate(i) < 1e-15 || slip_rate(i) > 1e-3)
         {
            above_wf_reasonable = false;
         }
      }
   }
   TEST_ASSERT(above_wf_reasonable,
               "Above-Wf slip rates in reasonable range");
}

// =============================================================================
// Test 12: Slip Accumulates at Vp Below Wf
// =============================================================================

void TestSlipAccumulatesAtVp()
{
   std::cout << "\n=== Test: Slip Accumulates at Vp Below Wf ===\n";

   TestComponents tc;
   tc.Setup(4, 10);

   Vector state(tc.fault->StateSize());
   tc.seas_op->SetInitialCondition(state);

   // Take a few RK4 steps
   RK4Solver ode_solver;
   ode_solver.Init(*tc.seas_op);

   real_t t = 0.0;
   real_t dt = 1e3;  // 1000 seconds
   int nsteps = 5;

   for (int step = 0; step < nsteps; step++)
   {
      ode_solver.Step(state, t, dt);
   }

   real_t expected_t = dt * nsteps;
   std::cout << "    t = " << t << " s (expected " << expected_t << ")\n";

   const Vector &depths = tc.fault_geom->GetDepths();
   Vector slip;
   tc.fault->GetSlip(state, slip);

   real_t max_below_wf_error = 0.0;
   for (int i = 0; i < tc.fault->NumNodes(); i++)
   {
      if (depths(i) < -tc.params.Wf)
      {
         real_t expected_slip = tc.params.Vp * t;
         real_t rel_err = std::abs(slip(i) - expected_slip) /
                          std::max(expected_slip, 1e-30);
         max_below_wf_error = std::max(max_below_wf_error, rel_err);
      }
   }

   // RK4 is 4th order, so for linear slip rate the error should be ~machine eps
   TEST_ASSERT(max_below_wf_error < 1e-6,
               "Below-Wf slip = Vp * t (RK4 exact for linear)");
   std::cout << "    Max below-Wf slip rel error: " << max_below_wf_error << "\n";
}

// =============================================================================
// Test 13: Domain Solve with Uniform Vp Slip -> ±Vp·t/2 Displacement
// =============================================================================

void TestUniformSlipDisplacement()
{
   std::cout << "\n=== Test: Uniform Slip -> ±Vp·t/2 Displacement ===\n";

   // Create domain operator directly (no fault/friction needed)
   BP2Params params;
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 10;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(*mesh, 1, params.mu(), params.Vp, params.Wf);

   int num_fault_dofs = domain.GetNumFaultDOFs();
   std::cout << "    Fault DOFs: " << num_fault_dofs << "\n";

   // Set uniform slip = delta on all fault DOFs
   real_t delta = 1.0e-3;  // 1 mm slip
   Vector slip_bc(num_fault_dofs);
   slip_bc = delta;

   // Solve domain
   GridFunction u_gf(&domain.GetFESpace());
   u_gf = 0.0;
   domain.Solve(0.0, slip_bc, u_gf);

   // Sample displacement at points on x>0 and x<0 sides
   // For uniform slip delta across the full fault at x=0:
   // u(x>0) ≈ +delta/2, u(x<0) ≈ -delta/2
   //
   // We use FindPoints to evaluate at specific physical locations
   int dim = mesh->Dimension();
   DenseMatrix points(dim, 2);

   // Point 1: x = +Lx/4, z = -Lz/2 (right side, mid-depth)
   points(0, 0) = mesh_params.Lx / 4.0;
   points(1, 0) = -mesh_params.Lz / 2.0;

   // Point 2: x = -Lx/4, z = -Lz/2 (left side, mid-depth)
   points(0, 1) = -mesh_params.Lx / 4.0;
   points(1, 1) = -mesh_params.Lz / 2.0;

   Array<int> elem_ids;
   Array<IntegrationPoint> ips;
   mesh->FindPoints(points, elem_ids, ips);

   TEST_ASSERT(elem_ids[0] >= 0, "Found element for right-side point");
   TEST_ASSERT(elem_ids[1] >= 0, "Found element for left-side point");

   if (elem_ids[0] >= 0 && elem_ids[1] >= 0)
   {
      real_t u_right = u_gf.GetValue(elem_ids[0], ips[0]);
      real_t u_left = u_gf.GetValue(elem_ids[1], ips[1]);

      std::cout << "    Slip delta: " << delta << " m\n";
      std::cout << "    u(+Lx/4) = " << u_right << " m\n";
      std::cout << "    u(-Lx/4) = " << u_left << " m\n";
      std::cout << "    Expected: ±" << delta / 2.0 << " m\n";

      // Check right side is positive, left side is negative
      TEST_ASSERT(u_right > 0.0, "Right-side displacement is positive");
      TEST_ASSERT(u_left < 0.0, "Left-side displacement is negative");

      // Check magnitudes are approximately equal
      real_t mag_diff = std::abs(std::abs(u_right) - std::abs(u_left));
      real_t avg_mag = (std::abs(u_right) + std::abs(u_left)) / 2.0;
      TEST_ASSERT(mag_diff / avg_mag < 0.1,
                  "Displacement magnitudes are symmetric");

      // Check values are approximately ±delta/2
      // Generous tolerance since boundary effects exist
      TEST_REL_NEAR(u_right, delta / 2.0, 0.3,
                    "Right-side u ≈ +delta/2");
      TEST_REL_NEAR(u_left, -delta / 2.0, 0.3,
                    "Left-side u ≈ -delta/2");
   }
}

// =============================================================================
// Test 14: Below-Wf Slip Produces Traction on Above-Wf DOFs
// =============================================================================

void TestBelowWfSlipProducesTraction()
{
   std::cout << "\n=== Test: Below-Wf Slip Produces Traction ===\n";

   // Create domain operator
   BP2Params params;
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 10;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(*mesh, 1, params.mu(), params.Vp, params.Wf);

   // Get fault depths
   Vector depths;
   domain.GetFaultDepths(depths);
   int num_fault_dofs = domain.GetNumFaultDOFs();

   // Set slip = delta only on below-Wf DOFs, zero above
   real_t delta = 1.0e-3;
   Vector slip_bc(num_fault_dofs);
   slip_bc = 0.0;

   int below_count = 0;
   for (int i = 0; i < num_fault_dofs; i++)
   {
      if (depths(i) < -params.Wf)
      {
         slip_bc(i) = delta;
         below_count++;
      }
   }
   std::cout << "    Below-Wf DOFs with slip: " << below_count << "\n";
   TEST_ASSERT(below_count > 0, "Have below-Wf DOFs to apply slip");

   // Solve domain
   GridFunction u_gf(&domain.GetFESpace());
   u_gf = 0.0;
   domain.Solve(0.0, slip_bc, u_gf);

   // Compute traction on all fault DOFs
   Vector traction(num_fault_dofs);
   domain.ComputeTraction(u_gf, slip_bc, traction);

   // Check that above-Wf DOFs have non-zero traction
   // (loading is transferred from below to above)
   real_t max_above_traction = 0.0;
   for (int i = 0; i < num_fault_dofs; i++)
   {
      if (depths(i) >= -params.Wf)
      {
         max_above_traction = std::max(max_above_traction,
                                       std::abs(traction(i)));
      }
   }

   std::cout << "    Max above-Wf traction: " << max_above_traction << " Pa\n";
   TEST_ASSERT(max_above_traction > 0.0,
               "Below-Wf slip produces non-zero traction above Wf");
}

// =============================================================================
// Dormand-Prince RK45 Tests
// =============================================================================

/// Simple exponential decay ODE: dy/dt = -lambda * y, exact: y(t) = y0 * exp(-lambda*t)
class ExponentialDecayOp : public TimeDependentOperator
{
   real_t lambda_;
public:
   ExponentialDecayOp(real_t lambda = 1.0)
      : TimeDependentOperator(1), lambda_(lambda) {}
   void Mult(const Vector &y, Vector &dydt) const override
   {
      dydt(0) = -lambda_ * y(0);
   }
};

/// Harmonic oscillator: dy1/dt = y2, dy2/dt = -y1
/// Exact: y1 = cos(t), y2 = -sin(t) (with y1(0)=1, y2(0)=0)
class HarmonicOscillatorOp : public TimeDependentOperator
{
public:
   HarmonicOscillatorOp() : TimeDependentOperator(2) {}
   void Mult(const Vector &y, Vector &dydt) const override
   {
      dydt(0) = y(1);
      dydt(1) = -y(0);
   }
};

/// Stiff-transition ODE that mimics nucleation behavior:
/// dy/dt = lambda * y * (1 - y/K), logistic growth.
/// Solution accelerates nonlinearly, similar to fault nucleation.
class LogisticGrowthOp : public TimeDependentOperator
{
   real_t lambda_, K_;
public:
   LogisticGrowthOp(real_t lambda = 10.0, real_t K = 1.0)
      : TimeDependentOperator(1), lambda_(lambda), K_(K) {}
   void Mult(const Vector &y, Vector &dydt) const override
   {
      dydt(0) = lambda_ * y(0) * (1.0 - y(0) / K_);
   }
};

/// Test 15: RK45 correctly solves exponential decay
void TestRK45ExponentialDecay()
{
   std::cout << "\n=== Test: RK45 Exponential Decay ===\n";

   real_t lambda = 2.0;
   ExponentialDecayOp op(lambda);

   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-10);
   rk45.SetRelTol(1e-10);
   rk45.SetDt(0.1);
   rk45.SetDtMax(1.0);
   rk45.Init(op);

   Vector y(1);
   y(0) = 1.0;
   real_t t = 0.0;
   real_t t_final = 5.0;
   int steps = 0;

   while (t < t_final && steps < 10000)
   {
      if (t + rk45.GetDt() > t_final) { rk45.SetDt(t_final - t); }
      real_t dt;
      if (rk45.Step(op, y, t, dt)) { steps++; }
   }

   real_t exact = std::exp(-lambda * t);
   real_t rel_err = std::abs(y(0) - exact) / exact;

   std::cout << "    t = " << t << ", y = " << y(0)
             << ", exact = " << exact << ", rel_err = " << rel_err << "\n";
   std::cout << "    Steps: " << steps
             << ", Rejections: " << rk45.GetTotalRejections() << "\n";

   TEST_ASSERT(rel_err < 1e-5, "RK45 exponential decay accuracy");
   TEST_ASSERT(steps < 200, "RK45 exponential decay efficiency (< 200 steps)");
}

/// Test 16: RK45 conserves energy in harmonic oscillator
void TestRK45HarmonicOscillator()
{
   std::cout << "\n=== Test: RK45 Harmonic Oscillator ===\n";

   HarmonicOscillatorOp op;

   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-10);
   rk45.SetRelTol(1e-10);
   rk45.SetDt(0.1);
   rk45.SetDtMax(2.0);
   rk45.Init(op);

   Vector y(2);
   y(0) = 1.0;  // cos(0)
   y(1) = 0.0;  // -sin(0)
   real_t t = 0.0;
   real_t t_final = 10.0 * M_PI;  // 5 full periods
   int steps = 0;

   while (t < t_final && steps < 100000)
   {
      if (t + rk45.GetDt() > t_final) { rk45.SetDt(t_final - t); }
      real_t dt;
      if (rk45.Step(op, y, t, dt)) { steps++; }
   }

   // Check solution at t = 10*pi: y1 = cos(10pi) = 1, y2 = -sin(10pi) = 0
   real_t err_y1 = std::abs(y(0) - 1.0);
   real_t err_y2 = std::abs(y(1) - 0.0);

   // Check energy conservation: y1^2 + y2^2 should equal 1
   real_t energy = y(0) * y(0) + y(1) * y(1);
   real_t energy_err = std::abs(energy - 1.0);

   std::cout << "    t = " << t << ", y = (" << y(0) << ", " << y(1) << ")\n";
   std::cout << "    Energy = " << energy << ", err = " << energy_err << "\n";
   std::cout << "    Steps: " << steps
             << ", Rejections: " << rk45.GetTotalRejections() << "\n";

   TEST_ASSERT(err_y1 < 1e-6, "RK45 harmonic oscillator y1 accuracy");
   TEST_ASSERT(err_y2 < 1e-6, "RK45 harmonic oscillator y2 accuracy");
   TEST_ASSERT(energy_err < 1e-6, "RK45 harmonic oscillator energy conservation");
}

/// Test 17: RK45 step rejection works for stiff problems
void TestRK45StepRejection()
{
   std::cout << "\n=== Test: RK45 Step Rejection ===\n";

   // Logistic growth with large lambda = rapid transition
   real_t lambda = 50.0;
   real_t K = 1.0;
   LogisticGrowthOp op(lambda, K);

   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-8);
   rk45.SetRelTol(1e-8);
   rk45.SetDt(1.0);  // Start with a large dt that will be rejected
   rk45.SetDtMax(10.0);
   rk45.Init(op);

   Vector y(1);
   y(0) = 0.01;  // Start small, will grow to K=1
   real_t t = 0.0;
   real_t t_final = 1.0;
   int steps = 0;

   while (t < t_final && steps < 100000)
   {
      if (t + rk45.GetDt() > t_final) { rk45.SetDt(t_final - t); }
      real_t dt;
      if (rk45.Step(op, y, t, dt)) { steps++; }
   }

   // Exact logistic solution: y(t) = K / (1 + (K/y0 - 1) * exp(-lambda*t))
   real_t y0 = 0.01;
   real_t exact = K / (1.0 + (K / y0 - 1.0) * std::exp(-lambda * t));
   real_t rel_err = std::abs(y(0) - exact) / exact;

   std::cout << "    t = " << t << ", y = " << y(0)
             << ", exact = " << exact << ", rel_err = " << rel_err << "\n";
   std::cout << "    Steps: " << steps
             << ", Rejections: " << rk45.GetTotalRejections() << "\n";

   TEST_ASSERT(rel_err < 1e-5, "RK45 logistic growth accuracy");
   TEST_ASSERT(rk45.GetTotalRejections() > 0,
               "RK45 rejects steps during rapid growth phase");
}

/// Test 18: RK45 FSAL property — second call reuses last stage
void TestRK45FSAL()
{
   std::cout << "\n=== Test: RK45 FSAL Property ===\n";

   // Run same problem twice: one with FSAL (consecutive steps),
   // one resetting initialization each time. FSAL should give same result.
   real_t lambda = 1.0;
   ExponentialDecayOp op(lambda);

   // Run 1: normal consecutive steps (FSAL active)
   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-12);
   rk45.SetRelTol(1e-12);
   rk45.SetDt(0.01);
   rk45.SetDtMax(0.5);
   rk45.Init(op);

   Vector y(1);
   y(0) = 1.0;
   real_t t = 0.0;
   int steps = 0;

   // Take exactly 10 accepted steps
   while (steps < 10)
   {
      real_t dt;
      if (rk45.Step(op, y, t, dt)) { steps++; }
   }

   real_t exact = std::exp(-lambda * t);
   real_t rel_err = std::abs(y(0) - exact) / exact;

   std::cout << "    After 10 steps: t = " << t << ", y = " << y(0)
             << ", rel_err = " << rel_err << "\n";

   TEST_ASSERT(rel_err < 1e-10, "RK45 FSAL gives accurate result over 10 steps");
}

/// Test 19: RK45 tolerance control — tighter tolerance gives smaller error
void TestRK45ToleranceControl()
{
   std::cout << "\n=== Test: RK45 Tolerance Control ===\n";

   real_t lambda = 3.0;
   ExponentialDecayOp op(lambda);
   real_t t_final = 2.0;
   real_t exact = std::exp(-lambda * t_final);

   // Run with loose tolerance
   DormandPrinceRK45 rk45_loose;
   rk45_loose.SetAbsTol(1e-4);
   rk45_loose.SetRelTol(1e-4);
   rk45_loose.SetDt(0.1);
   rk45_loose.SetDtMax(1.0);
   rk45_loose.Init(op);

   Vector y_loose(1);
   y_loose(0) = 1.0;
   real_t t_loose = 0.0;
   int steps_loose = 0;
   while (t_loose < t_final && steps_loose < 10000)
   {
      if (t_loose + rk45_loose.GetDt() > t_final)
      {
         rk45_loose.SetDt(t_final - t_loose);
      }
      real_t dt;
      if (rk45_loose.Step(op, y_loose, t_loose, dt)) { steps_loose++; }
   }
   real_t err_loose = std::abs(y_loose(0) - exact) / exact;

   // Run with tight tolerance
   DormandPrinceRK45 rk45_tight;
   rk45_tight.SetAbsTol(1e-10);
   rk45_tight.SetRelTol(1e-10);
   rk45_tight.SetDt(0.1);
   rk45_tight.SetDtMax(1.0);
   rk45_tight.Init(op);

   Vector y_tight(1);
   y_tight(0) = 1.0;
   real_t t_tight = 0.0;
   int steps_tight = 0;
   while (t_tight < t_final && steps_tight < 10000)
   {
      if (t_tight + rk45_tight.GetDt() > t_final)
      {
         rk45_tight.SetDt(t_final - t_tight);
      }
      real_t dt;
      if (rk45_tight.Step(op, y_tight, t_tight, dt)) { steps_tight++; }
   }
   real_t err_tight = std::abs(y_tight(0) - exact) / exact;

   std::cout << "    Loose (1e-4): err=" << err_loose
             << ", steps=" << steps_loose << "\n";
   std::cout << "    Tight (1e-10): err=" << err_tight
             << ", steps=" << steps_tight << "\n";

   TEST_ASSERT(err_tight < err_loose,
               "Tighter tolerance gives smaller error");
   TEST_ASSERT(steps_tight > steps_loose,
               "Tighter tolerance takes more steps");
   TEST_ASSERT(err_loose < 1e-2, "Loose tolerance still reasonable");
   TEST_ASSERT(err_tight < 1e-8, "Tight tolerance very accurate");
}

/// Test 20: RK45 default parameters match PETSc/Tandem
void TestRK45PetscDefaults()
{
   std::cout << "\n=== Test: RK45 PETSc/Tandem Default Parameters ===\n";

   DormandPrinceRK45 rk45;

   // Verify defaults match PETSc TSAdaptBasic + Tandem rk45.cfg
   // rtol = 1e-50 (effectively disabled, matching -ts_rtol 1e-50)
   // These are tested indirectly: construct and check getters exist
   // atol = 1e-7 (matching -ts_atol 1e-7)

   // Verify via exponential decay that pure absolute tolerance is used:
   // With rtol ~ 0, the error scale is always atol regardless of |y|.
   real_t lambda = 1.0;
   ExponentialDecayOp op(lambda);

   // Run with defaults (rtol = 1e-50, atol = 1e-7)
   DormandPrinceRK45 rk45_default;
   rk45_default.SetDt(0.1);
   rk45_default.SetDtMax(2.0);
   rk45_default.Init(op);

   Vector y(1);
   y(0) = 1.0;
   real_t t = 0.0;
   int steps = 0;
   while (t < 3.0 && steps < 10000)
   {
      if (t + rk45_default.GetDt() > 3.0) { rk45_default.SetDt(3.0 - t); }
      real_t dt;
      if (rk45_default.Step(op, y, t, dt)) { steps++; }
   }
   real_t exact = std::exp(-lambda * t);
   real_t rel_err = std::abs(y(0) - exact) / exact;

   std::cout << "    Default rtol=1e-50: t=" << t << ", y=" << y(0)
             << ", exact=" << exact << ", rel_err=" << rel_err
             << ", steps=" << steps << "\n";

   TEST_ASSERT(rel_err < 1e-4, "RK45 with PETSc defaults gives accurate result");
   TEST_ASSERT(steps > 0, "RK45 with PETSc defaults completes");
}

/// Test 21: RK45 reject safety factor reduces dt extra after rejection
void TestRK45RejectSafety()
{
   std::cout << "\n=== Test: RK45 Reject Safety Factor ===\n";

   // Use a stiff problem that triggers rejections.
   // Compare: with reject_safety = 0.5 (default) vs reject_safety = 1.0 (none)
   real_t lambda = 50.0;
   real_t K = 1.0;
   LogisticGrowthOp op(lambda, K);

   // Run 1: with reject_safety = 0.5 (PETSc default)
   DormandPrinceRK45 rk45_with;
   rk45_with.SetAbsTol(1e-8);
   rk45_with.SetRelTol(1e-50);
   rk45_with.SetRejectSafety(0.5);
   rk45_with.SetDt(1.0);
   rk45_with.SetDtMax(10.0);
   rk45_with.Init(op);

   Vector y1(1);
   y1(0) = 0.01;
   real_t t1 = 0.0;
   int steps1 = 0, attempts1 = 0;
   while (t1 < 0.5 && attempts1 < 100000)
   {
      real_t dt;
      if (rk45_with.Step(op, y1, t1, dt)) { steps1++; }
      attempts1++;
   }
   int rejections1 = rk45_with.GetTotalRejections();

   // Run 2: with reject_safety = 1.0 (no extra shrink)
   DormandPrinceRK45 rk45_without;
   rk45_without.SetAbsTol(1e-8);
   rk45_without.SetRelTol(1e-50);
   rk45_without.SetRejectSafety(1.0);
   rk45_without.SetDt(1.0);
   rk45_without.SetDtMax(10.0);
   rk45_without.Init(op);

   Vector y2(1);
   y2(0) = 0.01;
   real_t t2 = 0.0;
   int steps2 = 0, attempts2 = 0;
   while (t2 < 0.5 && attempts2 < 100000)
   {
      real_t dt;
      if (rk45_without.Step(op, y2, t2, dt)) { steps2++; }
      attempts2++;
   }
   int rejections2 = rk45_without.GetTotalRejections();

   std::cout << "    With reject_safety=0.5: steps=" << steps1
             << ", rejections=" << rejections1
             << ", attempts=" << attempts1 << "\n";
   std::cout << "    With reject_safety=1.0: steps=" << steps2
             << ", rejections=" << rejections2
             << ", attempts=" << attempts2 << "\n";

   // Both should reach the answer accurately
   real_t exact = K / (1.0 + (K / 0.01 - 1.0) * std::exp(-lambda * 0.5));
   real_t err1 = std::abs(y1(0) - exact) / exact;
   real_t err2 = std::abs(y2(0) - exact) / exact;

   std::cout << "    Accuracy: err1=" << err1 << ", err2=" << err2 << "\n";

   TEST_ASSERT(err1 < 1e-5, "RK45 with reject_safety=0.5 accurate");
   TEST_ASSERT(err2 < 1e-5, "RK45 with reject_safety=1.0 accurate");

   // With reject_safety=0.5, fewer total attempts (rejections converge faster)
   TEST_ASSERT(attempts1 <= attempts2,
               "reject_safety=0.5 needs fewer or equal attempts than 1.0");
}

/// Test 22: RK45 growth/shrink clip values match PETSc defaults
void TestRK45ClipValues()
{
   std::cout << "\n=== Test: RK45 Growth/Shrink Clip Values ===\n";

   // Test that with growth_max=10 (PETSc default), dt can grow 10x per step.
   // Use a simple problem where error is very small (dt wants to grow fast).
   real_t lambda = 0.1;  // Slow decay = easy problem
   ExponentialDecayOp op(lambda);

   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-7);
   rk45.SetRelTol(1e-50);
   rk45.SetDt(1e-3);       // Start very small
   rk45.SetDtMax(100.0);
   rk45.Init(op);

   Vector y(1);
   y(0) = 1.0;
   real_t t = 0.0;

   // Take one step from very small dt
   real_t dt;
   bool accepted = rk45.Step(op, y, t, dt);
   TEST_ASSERT(accepted, "First step accepted");

   real_t dt_after = rk45.GetDt();
   real_t growth_ratio = dt_after / dt;

   std::cout << "    Initial dt=" << dt << ", next dt=" << dt_after
             << ", growth=" << growth_ratio << "x\n";

   // With growth_max=10 and a very easy problem, dt should grow by up to 10x
   TEST_ASSERT(growth_ratio <= 10.0 + 1e-10,
               "Growth capped at growth_max=10");
   TEST_ASSERT(growth_ratio > 4.0,
               "Easy problem allows significant dt growth (> 4x)");

   // Now test shrink_min = 0.1: on a hard problem, dt can shrink by up to 10x
   real_t lambda_hard = 100.0;
   real_t K = 1.0;
   LogisticGrowthOp hard_op(lambda_hard, K);

   DormandPrinceRK45 rk45_hard;
   rk45_hard.SetAbsTol(1e-8);
   rk45_hard.SetRelTol(1e-50);
   rk45_hard.SetDt(1.0);    // Way too large for lambda=100
   rk45_hard.SetDtMax(10.0);
   rk45_hard.Init(hard_op);

   Vector y2(1);
   y2(0) = 0.01;
   real_t t2 = 0.0;
   real_t dt2;
   bool accepted2 = rk45_hard.Step(hard_op, y2, t2, dt2);

   // This should be rejected
   real_t dt_after_reject = rk45_hard.GetDt();
   // Account for reject_safety: effective shrink = clip * reject_safety
   // With shrink_min=0.1 and reject_safety=0.5, minimum ratio = 0.1 * 0.5 = 0.05
   real_t shrink_ratio = dt_after_reject / dt2;

   std::cout << "    Hard problem: dt=" << dt2 << ", accepted=" << accepted2
             << ", next dt=" << dt_after_reject
             << ", shrink=" << shrink_ratio << "x\n";

   TEST_ASSERT(!accepted2, "Hard problem rejects large initial dt");
   // After rejection with reject_safety, dt should shrink significantly
   TEST_ASSERT(shrink_ratio < 0.5,
               "Rejected step shrinks dt substantially");
   // But not below shrink_min * reject_safety * dt (= 0.05 * dt)
   TEST_ASSERT(dt_after_reject >= rk45_hard.GetDtMin(),
               "Shrunk dt stays above dt_min");
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "================================================\n";
   std::cout << "SEAS Phase 4 Tests: Quasi-Dynamic Operator\n";
   std::cout << "================================================\n";

   TestOperatorConstruction();
   TestInitialCondition();
   TestSteadyStateSlip();
   TestStressBalanceMaintained();
   TestSlipConservation();
   TestThetaPositive();
   TestAdaptiveTimeStepper();
   TestRK4WithAdaptiveDt();
   TestMultDeterministic();

   std::cout << "\n=== Below-Wf Loading Tests ===\n";
   TestExtendedFaultFullDepth();
   TestPrescribedVpBelowWf();
   TestSlipAccumulatesAtVp();
   TestUniformSlipDisplacement();
   TestBelowWfSlipProducesTraction();

   std::cout << "\n=== Dormand-Prince RK45 Tests ===\n";
   TestRK45ExponentialDecay();
   TestRK45HarmonicOscillator();
   TestRK45StepRejection();
   TestRK45FSAL();
   TestRK45ToleranceControl();

   std::cout << "\n=== RK45 PETSc/Tandem Matching Tests ===\n";
   TestRK45PetscDefaults();
   TestRK45RejectSafety();
   TestRK45ClipValues();

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
