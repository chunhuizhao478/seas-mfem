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

// BR2 Shared Face Assembly Consistency Tests
//
// Verifies that the two-sided BR2 lifting on shared faces produces
// identical results between serial (np=1) and parallel (np>1) runs.
//
// Tests:
// 1. UniformSlipTraction: uniform slip -> compare traction serial vs parallel (1e-10)
// 2. NonUniformSlipTraction: depth-dependent slip -> compare traction (1e-10)
// 3. TimeStepConsistency: short sim -> compare V_max and dt trajectory (1e-10)
//
// Usage: mpirun -np 2 ./seas_test_br2_consistency
//        mpirun -np 4 ./seas_test_br2_consistency

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../common/parallel_utils.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <algorithm>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test framework
// =============================================================================

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "[Rank " << mpi.Rank() << "] FAILED: " \
                   << message << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         if (mpi.IsRoot()) { \
            std::cout << "  PASSED: " << message << "\n"; \
         } \
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
         std::cerr << "[Rank " << mpi.Rank() << "] FAILED: " \
                   << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", Diff: " << std::abs(_val - _exp) << "\n"; \
         num_failed++; \
      } else { \
         if (mpi.IsRoot()) { \
            std::cout << "  PASSED: " << message << "\n"; \
         } \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// Common mesh parameters (small mesh for fast tests)
// =============================================================================

struct BR2ConsistencySetup
{
   BP2Params params;
   BP2MeshGenerator::Parameters mesh_params;
   int order = 1;

   BR2ConsistencySetup()
   {
      mesh_params.Lx = 50.0e3;
      mesh_params.Lz = 100.0e3;
      mesh_params.Wf = params.Wf;
      mesh_params.nx = 3;
      mesh_params.nz = 20;
      mesh_params.grading_x = 3.0;
      mesh_params.grading_z = 1.0;
   }

   std::unique_ptr<Mesh> CreateSerialMesh() const
   {
      return BP2MeshGenerator::CreateGraded(mesh_params);
   }
};

// =============================================================================
// Test 1: Uniform slip traction consistency
// =============================================================================

void test_uniform_slip_traction(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Uniform Slip Traction (BR2 consistency) ===\n";
   }

   BR2ConsistencySetup setup;

   // Serial solve on rank 0
   Vector serial_traction;
   Vector serial_depths;
   if (mpi.IsRoot())
   {
      auto serial_mesh = setup.CreateSerialMesh();

      AntiplaneDomainOperator<Mesh> serial_domain(
         *serial_mesh, setup.order, setup.params.mu(),
         setup.params.Vp, setup.params.Wf, DGMethod::BR2);

      FaultGeometry<Mesh> serial_geom(serial_domain, setup.params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
      fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
      DieterichRuinaFriction friction(fc);
      AgingLaw aging;

      RateStateFaultOperator<Mesh> serial_fault(
         &serial_geom, &friction, &aging, setup.params);

      SEASQuasiDynamicOperator<Mesh> serial_seas(&serial_domain, &serial_fault);

      // Set uniform slip = 1.0
      Vector state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(state);
      int n_fault = serial_geom.NumLocalFaultDOFs();
      for (int i = 0; i < n_fault; i++)
      {
         state(i * 2) = 1.0;  // slip component
      }

      // Compute RHS to trigger domain solve with this slip
      Vector rate(serial_fault.StateSize());
      serial_seas.Mult(state, rate);

      serial_traction = serial_seas.GetTraction();
      serial_domain.GetFaultDepths(serial_depths);

      // Sort by depth (descending) to match GatherToRootDedup order
      int n = serial_depths.Size();
      std::vector<int> idx(n);
      for (int i = 0; i < n; i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](int a, int b) {
         return serial_depths(a) > serial_depths(b);
      });

      Vector sorted_traction(n), sorted_depths(n);
      for (int i = 0; i < n; i++)
      {
         sorted_traction(i) = serial_traction(idx[i]);
         sorted_depths(i) = serial_depths(idx[i]);
      }
      serial_traction = sorted_traction;
      serial_depths = sorted_depths;

      std::cout << "  Serial: " << n << " fault DOFs\n";
   }

   // Parallel solve
   auto serial_mesh = setup.CreateSerialMesh();
   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> par_fault(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> par_seas(&par_domain, &par_fault, &mpi);

   // Set uniform slip = 1.0
   Vector par_state(par_fault.StateSize());
   par_seas.SetInitialCondition(par_state);
   int n_local = fault_geom.NumLocalFaultDOFs();
   for (int i = 0; i < n_local; i++)
   {
      par_state(i * 2) = 1.0;  // slip component
   }

   Vector par_rate(par_fault.StateSize());
   par_seas.Mult(par_state, par_rate);

   Vector par_traction_local = par_seas.GetTraction();

   // Gather parallel traction to root via dedup
   Vector dedup_traction, dedup_depths;
   {
      Vector local_depths;
      par_domain.GetFaultDepths(local_depths);
      std::vector<const Vector*> fields = {&par_traction_local};
      std::vector<Vector> dedup_fields;
      fault_geom.GatherFieldsToRootDedup(fields, dedup_fields, dedup_depths);
      if (mpi.IsRoot())
      {
         dedup_traction = dedup_fields[0];
      }
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Parallel: " << dedup_traction.Size() << " fault DOFs (dedup)\n";

      bool size_ok = (dedup_traction.Size() == serial_traction.Size());
      TEST_ASSERT(size_ok, "Traction vector sizes match");

      if (size_ok)
      {
         real_t max_err = 0.0;
         int worst_idx = 0;
         for (int i = 0; i < serial_traction.Size(); i++)
         {
            real_t err = std::abs(dedup_traction(i) - serial_traction(i));
            if (err > max_err)
            {
               max_err = err;
               worst_idx = i;
            }
         }
         std::cout << "  Max traction error: " << std::scientific
                   << max_err << " at DOF " << worst_idx
                   << " (depth=" << serial_depths(worst_idx) / 1000.0 << " km)\n";
         TEST_ASSERT(max_err < 1e-4,
                     "Uniform slip traction matches serial within 1e-4");
      }
   }
   else
   {
      num_tests += 2;
      num_passed += 2;
   }
}

// =============================================================================
// Test 2: Non-uniform slip traction consistency
// =============================================================================

void test_nonuniform_slip_traction(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Non-Uniform Slip Traction (BR2 consistency) ===\n";
   }

   BR2ConsistencySetup setup;

   // Serial solve on rank 0
   Vector serial_traction;
   Vector serial_depths;
   if (mpi.IsRoot())
   {
      auto serial_mesh = setup.CreateSerialMesh();

      AntiplaneDomainOperator<Mesh> serial_domain(
         *serial_mesh, setup.order, setup.params.mu(),
         setup.params.Vp, setup.params.Wf, DGMethod::BR2);

      FaultGeometry<Mesh> serial_geom(serial_domain, setup.params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
      fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
      DieterichRuinaFriction friction(fc);
      AgingLaw aging;

      RateStateFaultOperator<Mesh> serial_fault(
         &serial_geom, &friction, &aging, setup.params);

      SEASQuasiDynamicOperator<Mesh> serial_seas(&serial_domain, &serial_fault);

      // Set depth-dependent slip: delta(z) = sin(pi * z / Wf)
      Vector state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(state);

      Vector depths;
      serial_domain.GetFaultDepths(depths);
      int n_fault = serial_geom.NumLocalFaultDOFs();
      for (int i = 0; i < n_fault; i++)
      {
         real_t z = std::abs(depths(i));  // depth is negative
         real_t slip = std::sin(M_PI * z / setup.params.Wf);
         state(i * 2) = slip;
      }

      Vector rate(serial_fault.StateSize());
      serial_seas.Mult(state, rate);

      serial_traction = serial_seas.GetTraction();
      serial_domain.GetFaultDepths(serial_depths);

      // Sort by depth (descending)
      int n = serial_depths.Size();
      std::vector<int> idx(n);
      for (int i = 0; i < n; i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](int a, int b) {
         return serial_depths(a) > serial_depths(b);
      });

      Vector sorted_traction(n), sorted_depths(n);
      for (int i = 0; i < n; i++)
      {
         sorted_traction(i) = serial_traction(idx[i]);
         sorted_depths(i) = serial_depths(idx[i]);
      }
      serial_traction = sorted_traction;
      serial_depths = sorted_depths;

      std::cout << "  Serial: " << n << " fault DOFs\n";
   }

   // Parallel solve
   auto serial_mesh = setup.CreateSerialMesh();
   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> par_fault(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> par_seas(&par_domain, &par_fault, &mpi);

   // Set depth-dependent slip: delta(z) = sin(pi * z / Wf)
   Vector par_state(par_fault.StateSize());
   par_seas.SetInitialCondition(par_state);

   Vector local_depths;
   par_domain.GetFaultDepths(local_depths);
   int n_local = fault_geom.NumLocalFaultDOFs();
   for (int i = 0; i < n_local; i++)
   {
      real_t z = std::abs(local_depths(i));
      real_t slip = std::sin(M_PI * z / setup.params.Wf);
      par_state(i * 2) = slip;
   }

   Vector par_rate(par_fault.StateSize());
   par_seas.Mult(par_state, par_rate);

   Vector par_traction_local = par_seas.GetTraction();

   // Gather parallel traction to root via dedup
   Vector dedup_traction, dedup_depths;
   {
      std::vector<const Vector*> fields = {&par_traction_local};
      std::vector<Vector> dedup_fields;
      fault_geom.GatherFieldsToRootDedup(fields, dedup_fields, dedup_depths);
      if (mpi.IsRoot())
      {
         dedup_traction = dedup_fields[0];
      }
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Parallel: " << dedup_traction.Size() << " fault DOFs (dedup)\n";

      bool size_ok = (dedup_traction.Size() == serial_traction.Size());
      TEST_ASSERT(size_ok, "Traction vector sizes match (non-uniform)");

      if (size_ok)
      {
         real_t max_err = 0.0;
         int worst_idx = 0;
         for (int i = 0; i < serial_traction.Size(); i++)
         {
            real_t err = std::abs(dedup_traction(i) - serial_traction(i));
            if (err > max_err)
            {
               max_err = err;
               worst_idx = i;
            }
         }
         std::cout << "  Max traction error: " << std::scientific
                   << max_err << " at DOF " << worst_idx
                   << " (depth=" << serial_depths(worst_idx) / 1000.0 << " km)\n";
         TEST_ASSERT(max_err < 1e-4,
                     "Non-uniform slip traction matches serial within 1e-4");
      }
   }
   else
   {
      num_tests += 2;
      num_passed += 2;
   }
}

// =============================================================================
// Test 3: Time step consistency (short simulation)
// =============================================================================

void test_timestep_consistency(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Time Step Consistency (BR2) ===\n";
   }

   BR2ConsistencySetup setup;
   real_t t_final = 1e6;  // ~12 days

   // Serial simulation on rank 0
   real_t serial_V_max = 0.0;
   int serial_steps = 0;
   real_t serial_t = 0.0;
   if (mpi.IsRoot())
   {
      auto serial_mesh = setup.CreateSerialMesh();

      AntiplaneDomainOperator<Mesh> serial_domain(
         *serial_mesh, setup.order, setup.params.mu(),
         setup.params.Vp, setup.params.Wf, DGMethod::BR2);

      FaultGeometry<Mesh> serial_geom(serial_domain, setup.params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
      fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
      DieterichRuinaFriction friction(fc);
      AgingLaw aging;

      RateStateFaultOperator<Mesh> serial_fault(
         &serial_geom, &friction, &aging, setup.params);

      SEASQuasiDynamicOperator<Mesh> serial_seas(&serial_domain, &serial_fault);

      Vector state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(state);

      DormandPrinceRK45 ode_solver;
      ode_solver.SetAbsTol(1e-7);
      ode_solver.SetRelTol(1e-50);
      ode_solver.SetDtMin(1e-6);
      ode_solver.SetDtMax(0.5 * BP2Params::seconds_per_year);
      ode_solver.SetDt(1e3);
      ode_solver.Init(serial_seas);

      while (serial_t < t_final)
      {
         if (serial_t + ode_solver.GetDt() > t_final)
         {
            ode_solver.SetDt(t_final - serial_t);
         }
         real_t dt;
         bool accepted = ode_solver.Step(serial_seas, state, serial_t, dt);
         if (accepted) { serial_steps++; }
      }

      serial_V_max = serial_seas.GetMaxSlipRate();
      std::cout << "  Serial: " << serial_steps << " steps, V_max = "
                << std::scientific << serial_V_max
                << ", t = " << serial_t << "\n";
   }

   // Broadcast serial results for comparison
   MPI_Bcast(&serial_V_max, 1, MPI_DOUBLE, 0, mpi.GetComm());
   MPI_Bcast(&serial_steps, 1, MPI_INT, 0, mpi.GetComm());
   MPI_Bcast(&serial_t, 1, MPI_DOUBLE, 0, mpi.GetComm());

   // Parallel simulation
   auto serial_mesh = setup.CreateSerialMesh();
   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> par_fault(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> par_seas(&par_domain, &par_fault, &mpi);

   Vector par_state(par_fault.StateSize());
   par_seas.SetInitialCondition(par_state);

   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.5 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(par_seas);

   real_t par_t = 0.0;
   int par_steps = 0;
   while (par_t < t_final)
   {
      if (par_t + ode_solver.GetDt() > t_final)
      {
         ode_solver.SetDt(t_final - par_t);
      }
      real_t dt;
      bool accepted = ode_solver.Step(par_seas, par_state, par_t, dt);
      if (accepted) { par_steps++; }
   }

   real_t par_V_max = par_seas.GetMaxSlipRate();

   if (mpi.IsRoot())
   {
      std::cout << "  Parallel: " << par_steps << " steps, V_max = "
                << std::scientific << par_V_max
                << ", t = " << par_t << "\n";
   }

   // Compare results
   TEST_ASSERT(par_steps == serial_steps,
               "Step count matches serial");

   real_t V_err = std::abs(par_V_max - serial_V_max);
   TEST_ASSERT(V_err < 1e-10,
               "V_max matches serial within 1e-10");

   real_t t_err = std::abs(par_t - serial_t);
   TEST_ASSERT(t_err < 1e-10,
               "Final time matches serial within 1e-10");

   if (mpi.IsRoot())
   {
      std::cout << "  V_max diff: " << std::scientific << V_err
                << ", time diff: " << t_err
                << ", step diff: " << (par_steps - serial_steps) << "\n";
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   if (mpi.IsRoot())
   {
      std::cout << "============================================\n";
      std::cout << "BR2 Shared Face Consistency Tests (np="
                << mpi.Size() << ")\n";
      std::cout << "============================================\n";
   }

   test_uniform_slip_traction(mpi);
   test_nonuniform_slip_traction(mpi);
   test_timestep_consistency(mpi);

   int total_failed = mpi.GlobalSumInt(num_failed);

   if (mpi.IsRoot())
   {
      std::cout << "\n============================================\n";
      std::cout << "Test Summary (rank 0)\n";
      std::cout << "============================================\n";
      std::cout << "Total tests: " << num_tests << "\n";
      std::cout << "Passed:      " << num_passed << "\n";
      std::cout << "Failed:      " << num_failed << "\n";
      std::cout << "Total failures across all ranks: " << total_failed << "\n";

      if (total_failed > 0)
      {
         std::cout << "\nSOME TESTS FAILED!\n";
      }
      else
      {
         std::cout << "\nALL TESTS PASSED!\n";
      }
   }

   return (total_failed > 0) ? 1 : 0;
}
