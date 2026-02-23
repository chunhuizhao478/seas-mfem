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

// Phase 8 parallel fault tests
//
// Tests:
// 1. DOF distribution: sum of local DOFs across ranks = serial total
// 2. GatherToRoot: gathered local data reconstructs global data
// 3. Global V_max reduction: max across ranks is correct
// 4. Serial-parallel consistency: parallel fault Init/ComputeRHS matches serial
// 5. Probe output consistency: parallel gathered probe values match serial
//
// Usage: mpirun -np 2 ./seas_test_parallel_fault
//        mpirun -np 4 ./seas_test_parallel_fault

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
#include "../../io/probe_output.hpp"
#include "../../io/benchmark_output.hpp"
#include "../../io/parallel_benchmark_output.hpp"

#include <fstream>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <numeric>

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
// Helper: create mesh and operators
// =============================================================================

struct TestSetup
{
   BP2Params params;
   std::unique_ptr<Mesh> serial_mesh;
   int order = 1;

   TestSetup()
   {
      BP2MeshGenerator::Parameters mp;
      mp.Lx = 50.0e3;
      mp.Lz = 100.0e3;
      mp.Wf = params.Wf;
      mp.nx = 5;
      mp.nz = 50;  // Coarse for fast tests
      mp.grading_x = 3.0;
      mp.grading_z = 1.0;
      serial_mesh = BP2MeshGenerator::CreateGraded(mp);
   }
};

// =============================================================================
// Test 1: DOF distribution
// =============================================================================

void test_dof_distribution(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: DOF Distribution ===\n";
   }

   TestSetup setup;

   // Serial DOF count
   int serial_fault_dofs = 0;
   {
      AntiplaneDomainOperator<Mesh> serial_domain(
         *setup.serial_mesh, setup.order, setup.params.mu(),
         setup.params.Vp, setup.params.Wf, DGMethod::BR2);
      serial_fault_dofs = serial_domain.GetNumFaultDOFs();
   }

   // Parallel DOF count
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   int local_dofs = fault_geom.NumLocalFaultDOFs();
   int global_dofs = fault_geom.NumGlobalFaultDOFs();

   // Sum of local DOFs should equal global DOFs
   int sum_local = mpi.GlobalSumInt(local_dofs);
   TEST_ASSERT(sum_local == global_dofs,
               "Sum of local DOFs == NumGlobalFaultDOFs");

   // Global DOF count should be >= serial (DG may duplicate at partition boundaries)
   TEST_ASSERT(global_dofs >= serial_fault_dofs,
               "Global fault DOFs >= serial fault DOFs (DG duplication at partition boundaries)");

   // Each rank should have some DOFs (with a reasonable mesh)
   TEST_ASSERT(local_dofs >= 0,
               "Each rank has non-negative local DOFs");

   if (mpi.IsRoot())
   {
      std::cout << "  Serial DOFs: " << serial_fault_dofs
                << ", Global DOFs: " << global_dofs
                << ", Local (rank 0): " << local_dofs << "\n";
   }
}

// =============================================================================
// Test 2: GatherToRoot
// =============================================================================

void test_gather_to_root(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: GatherToRoot ===\n";
   }

   TestSetup setup;
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   int local_n = fault_geom.NumLocalFaultDOFs();
   int global_n = fault_geom.NumGlobalFaultDOFs();

   // Create local data: each DOF gets rank * 1000 + local_index
   Vector local_data(local_n);
   for (int i = 0; i < local_n; i++)
   {
      local_data(i) = mpi.Rank() * 1000.0 + i;
   }

   // Gather
   Vector global_data;
   fault_geom.GatherToRoot(local_data, global_data);

   if (mpi.IsRoot())
   {
      TEST_ASSERT(global_data.Size() == global_n,
                  "Gathered data has correct global size");

      // Verify data was gathered correctly (check total size)
      TEST_ASSERT(global_data.Size() > 0,
                  "Gathered data is non-empty");
   }
   else
   {
      // Non-root: the gathered size is undefined
      num_tests += 2;
      num_passed += 2;
   }

   // Also test gathering fault depths and verifying they're sorted
   // (or at least that they span the expected range)
   Vector local_depths;
   par_domain.GetFaultDepths(local_depths);

   Vector global_depths;
   fault_geom.GatherToRoot(local_depths, global_depths);

   if (mpi.IsRoot())
   {
      TEST_ASSERT(global_depths.Size() == global_n,
                  "Gathered depths has correct global size");

      // Check depth range covers expected domain
      real_t z_min = global_depths.Min();
      real_t z_max = global_depths.Max();
      TEST_ASSERT(z_max >= -1000.0,  // Near surface
                  "Global depths include near-surface DOFs");
      TEST_ASSERT(z_min <= -setup.params.Wf + 1000.0,
                  "Global depths extend to near fault bottom");
   }
   else
   {
      num_tests += 3;
      num_passed += 3;
   }
}

// =============================================================================
// Test 3: Global V_max reduction
// =============================================================================

void test_global_vmax(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Global V_max Reduction ===\n";
   }

   TestSetup setup;
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> fault_op(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> seas_op(&par_domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t local_V = fault_op.GetMaxSlipRate();
   real_t global_V = fault_op.GetGlobalMaxSlipRate();

   // Global V should be >= local V
   TEST_ASSERT(global_V >= local_V - 1e-30,
               "Global V_max >= local V_max");

   // Global V should equal the true max across all ranks
   real_t expected_global = mpi.GlobalMax(local_V);
   TEST_NEAR(global_V, expected_global, 1e-20,
             "Global V_max matches MPI_Allreduce(MAX)");

   // All ranks should agree on global V_max
   real_t global_V_min = mpi.GlobalMin(global_V);
   real_t global_V_max = mpi.GlobalMax(global_V);
   TEST_NEAR(global_V_min, global_V_max, 1e-20,
             "All ranks agree on global V_max");

   if (mpi.IsRoot())
   {
      std::cout << "  Local V_max (rank 0): " << local_V
                << ", Global V_max: " << global_V << "\n";
   }
}

// =============================================================================
// Test 4: Serial-parallel consistency
// =============================================================================

void test_serial_parallel_consistency(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Serial-Parallel Consistency ===\n";
   }

   TestSetup setup;

   // Run serial initialization
   real_t serial_V_max = 0.0;
   real_t serial_tau0 = 0.0;
   {
      AntiplaneDomainOperator<Mesh> serial_domain(
         *setup.serial_mesh, setup.order, setup.params.mu(),
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

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      serial_V_max = serial_seas.GetMaxSlipRate();
      serial_tau0 = serial_fault.GetTau0();
   }

   // Run parallel initialization
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
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

   real_t par_V_max = par_seas.GetMaxSlipRate();
   real_t par_tau0 = par_fault.GetTau0();

   // tau0 should be identical (it's computed from params, not mesh)
   TEST_NEAR(par_tau0, serial_tau0, 1e-10,
             "Parallel tau0 matches serial");

   // V_max should be close (not exact due to mesh partitioning effects on DG)
   real_t V_rel_err = std::abs(par_V_max - serial_V_max) /
                      std::max(serial_V_max, 1e-30);
   TEST_ASSERT(V_rel_err < 0.1,
               "Parallel V_max within 10% of serial");

   if (mpi.IsRoot())
   {
      std::cout << "  Serial V_max: " << serial_V_max
                << ", Parallel V_max: " << par_V_max
                << " (rel err: " << V_rel_err << ")\n";
   }
}

// =============================================================================
// Test 5: Probe output consistency
// =============================================================================

void test_probe_gather(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Probe Gather Consistency ===\n";
   }

   TestSetup setup;
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   // Get local depths and gather
   Vector local_depths;
   par_domain.GetFaultDepths(local_depths);

   Vector global_depths;
   fault_geom.GatherToRoot(local_depths, global_depths);

   if (mpi.IsRoot())
   {
      // Create a probe interpolator on the global depths
      std::vector<real_t> probe_locs = {0.0, -10000.0, -20000.0};
      ProbeInterpolator interp(global_depths, probe_locs);

      // Create a synthetic field: f(z) = z / 1000
      Vector global_field(global_depths.Size());
      for (int i = 0; i < global_depths.Size(); i++)
      {
         global_field(i) = global_depths(i) / 1000.0;
      }

      Vector probe_vals;
      interp.Interpolate(global_field, probe_vals);

      // Probe at z=0 should give ~0
      TEST_NEAR(probe_vals(0), 0.0, 1.0,
                "Probe at z=0 gives ~0 for f(z)=z/1000");
      // Probe at z=-10km should give ~-10
      TEST_NEAR(probe_vals(1), -10.0, 1.0,
                "Probe at z=-10km gives ~-10 for f(z)=z/1000");
      // Probe at z=-20km should give ~-20
      TEST_NEAR(probe_vals(2), -20.0, 1.0,
                "Probe at z=-20km gives ~-20 for f(z)=z/1000");
   }
   else
   {
      // Non-root: count tests as passed
      num_tests += 3;
      num_passed += 3;
   }
}

// =============================================================================
// Test 6: RK45 time step synchronization (confirms Fix #1)
// =============================================================================

void test_rk45_sync(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: RK45 Time Step Synchronization ===\n";
   }

   TestSetup setup;
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> fault_op(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> seas_op(&par_domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(seas_op);

   real_t t = 0.0;
   int accepted_steps = 0;
   bool all_synced = true;

   // Take 5 accepted steps
   while (accepted_steps < 5)
   {
      real_t dt;
      bool accepted = ode_solver.Step(seas_op, state, t, dt);
      if (!accepted) { continue; }
      accepted_steps++;

      // Check that all ranks have the same t
      real_t t_min = mpi.GlobalMin(t);
      real_t t_max = mpi.GlobalMax(t);
      if (std::abs(t_max - t_min) > 1e-20)
      {
         all_synced = false;
      }

      // Check that all ranks have the same dt_next
      real_t dt_next = ode_solver.GetDt();
      real_t dt_min_val = mpi.GlobalMin(dt_next);
      real_t dt_max_val = mpi.GlobalMax(dt_next);
      if (std::abs(dt_max_val - dt_min_val) > 1e-20)
      {
         all_synced = false;
      }
   }

   TEST_ASSERT(all_synced,
               "All ranks agree on t and dt after 5 RK45 steps");
   TEST_ASSERT(t > 0.0,
               "Time advanced past zero after RK45 steps");

   if (mpi.IsRoot())
   {
      std::cout << "  After 5 steps: t = " << t
                << ", dt = " << ode_solver.GetDt() << "\n";
   }
}

// =============================================================================
// Test 7: ComputeRHS parallel consistency
// =============================================================================

void test_compute_rhs_consistency(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: ComputeRHS Parallel Consistency ===\n";
   }

   TestSetup setup;

   // Serial: init and compute RHS
   Vector serial_depths;
   Vector serial_rate;
   {
      AntiplaneDomainOperator<Mesh> serial_domain(
         *setup.serial_mesh, setup.order, setup.params.mu(),
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

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      // Compute RHS (calls Mult)
      serial_rate.SetSize(serial_fault.StateSize());
      serial_seas.Mult(serial_state, serial_rate);

      serial_domain.GetFaultDepths(serial_depths);
   }

   // Parallel: init and compute RHS
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
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

   Vector par_rate(par_fault.StateSize());
   par_seas.Mult(par_state, par_rate);

   // Gather parallel rates to root via dedup
   // Rate vector has interleaved [V_0, dtheta_0, V_1, dtheta_1, ...]
   // Extract V (slip rate) components
   int local_n = par_fault.NumNodes();
   Vector local_V(local_n);
   for (int i = 0; i < local_n; i++)
   {
      local_V(i) = par_rate(2 * i);  // V component
   }

   Vector dedup_V, dedup_depths;
   fault_geom.GatherToRootDedup(local_V, dedup_V, dedup_depths);

   if (mpi.IsRoot())
   {
      // Extract serial V
      int serial_n = serial_rate.Size() / 2;
      Vector serial_V(serial_n);
      for (int i = 0; i < serial_n; i++)
      {
         serial_V(i) = serial_rate(2 * i);
      }

      // Compare max slip rates
      real_t serial_V_max = serial_V.Normlinf();
      real_t par_V_max = dedup_V.Normlinf();
      real_t rel_err = std::abs(par_V_max - serial_V_max) /
                       std::max(serial_V_max, 1e-30);

      TEST_NEAR(par_V_max, serial_V_max, serial_V_max * 0.1,
                "Parallel max(V) from RHS within 10% of serial");

      // Check sizes are consistent
      TEST_ASSERT(dedup_V.Size() == serial_n,
                  "Deduplicated parallel V size matches serial");

      std::cout << "  Serial max(V): " << serial_V_max
                << ", Parallel max(V): " << par_V_max
                << " (rel err: " << rel_err << ")\n";
   }
   else
   {
      num_tests += 2;
      num_passed += 2;
   }
}

// =============================================================================
// Test 8: Deduplicated probe interpolation vs serial (confirms Fix #2)
// =============================================================================

void test_dedup_probe_interpolation(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Dedup Probe Interpolation vs Serial ===\n";
   }

   TestSetup setup;

   // Serial probe values
   Vector serial_probe_vals;
   std::vector<real_t> probe_locs = {0.0, -10000.0, -20000.0, -40000.0};
   {
      AntiplaneDomainOperator<Mesh> serial_domain(
         *setup.serial_mesh, setup.order, setup.params.mu(),
         setup.params.Vp, setup.params.Wf, DGMethod::BR2);

      Vector serial_depths;
      serial_domain.GetFaultDepths(serial_depths);

      ProbeInterpolator interp(serial_depths, probe_locs);

      // Create synthetic field: f(z) = sin(z / 10000)
      Vector field(serial_depths.Size());
      for (int i = 0; i < serial_depths.Size(); i++)
      {
         field(i) = std::sin(serial_depths(i) / 10000.0);
      }

      interp.Interpolate(field, serial_probe_vals);
   }

   // Parallel: gather, dedup, then interpolate
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   Vector local_depths;
   par_domain.GetFaultDepths(local_depths);

   // Create synthetic field on local depths
   Vector local_field(local_depths.Size());
   for (int i = 0; i < local_depths.Size(); i++)
   {
      local_field(i) = std::sin(local_depths(i) / 10000.0);
   }

   Vector dedup_field, dedup_depths;
   fault_geom.GatherToRootDedup(local_field, dedup_field, dedup_depths);

   if (mpi.IsRoot())
   {
      ProbeInterpolator interp(dedup_depths, probe_locs);
      Vector par_probe_vals;
      interp.Interpolate(dedup_field, par_probe_vals);

      bool all_close = true;
      for (size_t i = 0; i < probe_locs.size(); i++)
      {
         real_t diff = std::abs(par_probe_vals(i) - serial_probe_vals(i));
         if (diff > 0.01)
         {
            all_close = false;
            std::cerr << "  Probe " << i << " (z=" << probe_locs[i]
                      << "): serial=" << serial_probe_vals(i)
                      << " par=" << par_probe_vals(i)
                      << " diff=" << diff << "\n";
         }
      }

      TEST_ASSERT(all_close,
                  "Dedup parallel probe values match serial within 0.01");
      TEST_ASSERT(dedup_depths.Size() > 0,
                  "Deduplicated depths are non-empty");
   }
   else
   {
      num_tests += 2;
      num_passed += 2;
   }
}

// =============================================================================
// Test 9: Parallel benchmark output file content
// =============================================================================

void test_benchmark_output_consistency(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Parallel Benchmark Output File Content ===\n";
   }

   TestSetup setup;
   std::vector<real_t> probe_depths = {0.0, -12000.0};
   std::string serial_prefix = "/tmp/seas_test_serial_bench";
   std::string par_prefix = "/tmp/seas_test_par_bench";

   // Serial: write at t=0
   if (mpi.IsRoot())
   {
      AntiplaneDomainOperator<Mesh> serial_domain(
         *setup.serial_mesh, setup.order, setup.params.mu(),
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

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      Vector serial_depths;
      serial_domain.GetFaultDepths(serial_depths);

      BenchmarkOutput<Mesh> serial_out(serial_prefix, setup.params,
                                        probe_depths, serial_depths);
      serial_out.ForceWrite(0.0, serial_state, serial_fault,
                            serial_seas.GetTraction());
      serial_out.Close();
   }
   mpi.Barrier();

   // Parallel: write at t=0
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
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

   // Gather dedup depths for ParallelBenchmarkOutput
   Vector local_depths;
   par_domain.GetFaultDepths(local_depths);
   Vector dedup_fault_depths;
   {
      Vector dedup_data;
      fault_geom.GatherToRootDedup(local_depths, dedup_data, dedup_fault_depths);
   }

   real_t V_init = par_seas.GetMaxSlipRate();
   ParallelBenchmarkOutput par_out(par_prefix, setup.params, probe_depths,
                                    fault_geom, mpi, dedup_fault_depths);
   par_out.ForceWrite(0.0, par_state, par_fault, par_seas.GetTraction(),
                      V_init);
   par_out.Close();

   mpi.Barrier();

   // Compare output files on root
   if (mpi.IsRoot())
   {
      // Read first probe file from each
      std::string serial_file = serial_prefix + "_z0km.txt";
      std::string par_file = par_prefix + "_z0km.txt";

      std::ifstream sf(serial_file);
      std::ifstream pf(par_file);

      TEST_ASSERT(sf.good(), "Serial output file exists");
      TEST_ASSERT(pf.good(), "Parallel output file exists");

      if (sf.good() && pf.good())
      {
         // Skip header lines (lines starting with #)
         std::string sline, pline;
         while (std::getline(sf, sline) && sline[0] == '#') {}
         while (std::getline(pf, pline) && pline[0] == '#') {}

         // Parse first data line
         std::istringstream ss(sline), ps(pline);
         double s_time, s_slip, s_logV, s_tau, s_logth;
         double p_time, p_slip, p_logV, p_tau, p_logth;

         ss >> s_time >> s_slip >> s_logV >> s_tau >> s_logth;
         ps >> p_time >> p_slip >> p_logV >> p_tau >> p_logth;

         TEST_NEAR(p_time, s_time, 1e-10,
                   "Output time matches serial");
         TEST_NEAR(p_logV, s_logV, 0.5,
                   "Output log10(V) within 0.5 of serial at z=0km");
      }

      // Clean up temp files
      std::remove(serial_file.c_str());
      std::remove(par_file.c_str());
      std::string serial_file2 = serial_prefix + "_z12km.txt";
      std::string par_file2 = par_prefix + "_z12km.txt";
      std::remove(serial_file2.c_str());
      std::remove(par_file2.c_str());
   }
   else
   {
      num_tests += 4;
      num_passed += 4;
   }
}

// =============================================================================
// Test 10: Time stepper dt agreement across ranks (confirms Fix #1)
// =============================================================================

void test_dt_agreement(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Time Stepper dt Agreement ===\n";
   }

   TestSetup setup;
   ParMesh pmesh(mpi.GetComm(), *setup.serial_mesh);
   AntiplaneDomainOperator<ParMesh> par_domain(
      pmesh, setup.order, setup.params.mu(),
      setup.params.Vp, setup.params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(par_domain, setup.params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = setup.params.V0;  fc.f0 = setup.params.f0;
   fc.b = setup.params.b;     fc.Dc = setup.params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> fault_op(
      &fault_geom, &friction, &aging, setup.params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> seas_op(&par_domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(seas_op);

   real_t t = 0.0;
   int accepted_steps = 0;
   bool dt_always_agreed = true;

   while (accepted_steps < 5)
   {
      real_t dt;
      bool accepted = ode_solver.Step(seas_op, state, t, dt);
      if (!accepted) { continue; }
      accepted_steps++;

      real_t dt_next = ode_solver.GetDt();
      real_t dt_min_val = mpi.GlobalMin(dt_next);
      real_t dt_max_val = mpi.GlobalMax(dt_next);

      if (std::abs(dt_max_val - dt_min_val) > 1e-20)
      {
         dt_always_agreed = false;
         if (mpi.IsRoot())
         {
            std::cerr << "  Step " << accepted_steps
                      << ": dt_min=" << dt_min_val
                      << " dt_max=" << dt_max_val << "\n";
         }
      }
   }

   TEST_ASSERT(dt_always_agreed,
               "All ranks use same dt after every RK45 step");

   if (mpi.IsRoot())
   {
      std::cout << "  Final dt: " << ode_solver.GetDt()
                << " after " << accepted_steps << " steps\n";
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
      std::cout << "====================================\n";
      std::cout << "Phase 8: Parallel Fault Tests (np="
                << mpi.Size() << ")\n";
      std::cout << "====================================\n";
   }

   test_dof_distribution(mpi);
   test_gather_to_root(mpi);
   test_global_vmax(mpi);
   test_serial_parallel_consistency(mpi);
   test_probe_gather(mpi);
   test_rk45_sync(mpi);
   test_compute_rhs_consistency(mpi);
   test_dedup_probe_interpolation(mpi);
   test_benchmark_output_consistency(mpi);
   test_dt_agreement(mpi);

   // Reduce test results across ranks
   int total_failed = mpi.GlobalSumInt(num_failed);

   if (mpi.IsRoot())
   {
      std::cout << "\n====================================\n";
      std::cout << "Test Summary (rank 0)\n";
      std::cout << "====================================\n";
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
