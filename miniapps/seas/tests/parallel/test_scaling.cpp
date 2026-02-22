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

// Phase 9: Scaling Tests
//
// Tests (reporting only — always PASS):
// 1. StrongScaling: fixed mesh, measure wall time for 50 RK45 steps
// 2. WeakScaling: scale nz with sqrt(np), measure wall time
// 3. CommunicationOverhead: run 20 steps, report timing baseline
//
// Usage: mpirun -np 2 ./seas_test_scaling
//        mpirun -np 4 ./seas_test_scaling

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

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test framework
// =============================================================================

static int num_tests = 0;
static int num_passed = 0;

#define TEST_REPORT(message) \
   do { \
      num_tests++; \
      num_passed++; \
      if (mpi.IsRoot()) { \
         std::cout << "  PASSED: " << message << "\n"; \
      } \
   } while (0)

// =============================================================================
// Helper: run N accepted RK45 steps, return wall time
// =============================================================================

real_t RunSteps(MPIContext &mpi, int mesh_nx, int mesh_nz, int target_steps)
{
   BP2Params params;

   BP2MeshGenerator::Parameters mp;
   mp.Lx = 50.0e3;  mp.Lz = 100.0e3;  mp.Wf = params.Wf;
   mp.nx = mesh_nx;  mp.nz = mesh_nz;
   mp.grading_x = 3.0;  mp.grading_z = 1.0;

   auto serial_mesh = BP2MeshGenerator::CreateGraded(mp);
   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   int order = 1;
   AntiplaneDomainOperator<ParMesh> domain(
      pmesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;  fc.f0 = params.f0;
   fc.b = params.b;     fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<ParMesh> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   SEASQuasiDynamicOperator<ParMesh> seas_op(&domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-7);
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(seas_op);

   mpi.Barrier();
   double t_start = MPI_Wtime();

   real_t t = 0.0;
   int steps = 0;
   while (steps < target_steps)
   {
      real_t dt;
      bool accepted = ode_solver.Step(seas_op, state, t, dt);
      if (accepted) { steps++; }
   }

   mpi.Barrier();
   double t_end = MPI_Wtime();

   return t_end - t_start;
}

// =============================================================================
// Test 1: Strong scaling
// =============================================================================

void test_strong_scaling(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Strong Scaling ===\n";
   }

   int nx = 5, nz = 125;
   int target_steps = 50;

   real_t wall_time = RunSteps(mpi, nx, nz, target_steps);

   if (mpi.IsRoot())
   {
      std::cout << "  Mesh: " << 2*nx << "x" << nz
                << " (" << 2*nx*nz << " elements)\n";
      std::cout << "  Steps: " << target_steps << "\n";
      std::cout << "  np=" << mpi.Size()
                << ": wall time = " << std::fixed << std::setprecision(3)
                << wall_time << " s\n";
      std::cout << "  (Compare with different np values for speedup)\n";
   }

   TEST_REPORT("Strong scaling benchmark completed");
}

// =============================================================================
// Test 2: Weak scaling
// =============================================================================

void test_weak_scaling(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Weak Scaling ===\n";
   }

   int nx = 5;
   // Scale nz with sqrt(np) to keep work per rank roughly constant
   int base_nz = 50;
   int nz = static_cast<int>(base_nz * std::sqrt(static_cast<double>(mpi.Size())));
   int target_steps = 50;

   real_t wall_time = RunSteps(mpi, nx, nz, target_steps);

   if (mpi.IsRoot())
   {
      std::cout << "  Mesh: " << 2*nx << "x" << nz
                << " (" << 2*nx*nz << " elements, scaled with sqrt(np))\n";
      std::cout << "  Steps: " << target_steps << "\n";
      std::cout << "  np=" << mpi.Size()
                << ": wall time = " << std::fixed << std::setprecision(3)
                << wall_time << " s\n";
      std::cout << "  (Ideal: constant time as np increases)\n";
   }

   TEST_REPORT("Weak scaling benchmark completed");
}

// =============================================================================
// Test 3: Communication overhead baseline
// =============================================================================

void test_communication_overhead(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Communication Overhead Baseline ===\n";
   }

   int nx = 5, nz = 125;
   int target_steps = 20;

   real_t wall_time = RunSteps(mpi, nx, nz, target_steps);
   real_t time_per_step = wall_time / target_steps;

   if (mpi.IsRoot())
   {
      std::cout << "  Mesh: " << 2*nx << "x" << nz << "\n";
      std::cout << "  Steps: " << target_steps << "\n";
      std::cout << "  np=" << mpi.Size()
                << ": total = " << std::fixed << std::setprecision(3)
                << wall_time << " s"
                << ", per step = " << std::setprecision(4)
                << time_per_step << " s\n";
      std::cout << "  (Use MPI profiling tools for detailed breakdown)\n";
   }

   TEST_REPORT("Communication overhead baseline recorded");
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
      std::cout << "Phase 9: Scaling Tests (np=" << mpi.Size() << ")\n";
      std::cout << "============================================\n";
   }

   test_strong_scaling(mpi);
   test_weak_scaling(mpi);
   test_communication_overhead(mpi);

   if (mpi.IsRoot())
   {
      std::cout << "\n============================================\n";
      std::cout << "Scaling Summary\n";
      std::cout << "============================================\n";
      std::cout << "Total benchmarks: " << num_tests << "\n";
      std::cout << "All completed: " << num_passed << "\n";
      std::cout << "\nALL TESTS PASSED! (reporting tests)\n";
   }

   return 0;
}
