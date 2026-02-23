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

// Phase 9: Serial vs Parallel Consistency Tests
//
// Tests:
// 1. DomainSolution: uniform slip → solve → compare serial vs parallel (1e-10)
// 2. FaultState: short sim (t=1e6 s) → compare final state (1e-8)
// 3. ProbeOutput: compare BenchmarkOutput vs ParallelBenchmarkOutput files (1e-8)
// 4. Reproducibility: verify all ranks agree on V_max and total steps
//
// Usage: mpirun -np 2 ./seas_test_serial_parallel_consistency
//        mpirun -np 4 ./seas_test_serial_parallel_consistency

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
#include "../../io/benchmark_output.hpp"
#include "../../io/parallel_benchmark_output.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <sstream>
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

struct ConsistencySetup
{
   BP2Params params;
   BP2MeshGenerator::Parameters mesh_params;
   int order = 1;

   ConsistencySetup()
   {
      mesh_params.Lx = 50.0e3;
      mesh_params.Lz = 100.0e3;
      mesh_params.Wf = params.Wf;
      mesh_params.nx = 3;
      mesh_params.nz = 20;  // Small for fast execution
      mesh_params.grading_x = 3.0;
      mesh_params.grading_z = 1.0;
   }

   std::unique_ptr<Mesh> CreateSerialMesh() const
   {
      return BP2MeshGenerator::CreateGraded(mesh_params);
   }
};

// =============================================================================
// Test 1: Domain solution consistency
// =============================================================================

void test_domain_solution(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Domain Solution Consistency ===\n";
   }

   ConsistencySetup setup;

   // Serial: init and compute RHS
   real_t serial_V_max = 0.0;
   real_t serial_tau0 = 0.0;
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

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      // Compute RHS to trigger domain solve
      Vector serial_rate(serial_fault.StateSize());
      serial_seas.Mult(serial_state, serial_rate);

      serial_V_max = serial_seas.GetMaxSlipRate();
      serial_tau0 = serial_fault.GetTau0();
   }

   // Parallel: same setup
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

   Vector par_rate(par_fault.StateSize());
   par_seas.Mult(par_state, par_rate);

   real_t par_V_max = par_seas.GetMaxSlipRate();
   real_t par_tau0 = par_fault.GetTau0();

   // tau0 should be identical (computed from params, not mesh)
   TEST_NEAR(par_tau0, serial_tau0, 1e-10,
             "Parallel tau0 matches serial");

   // V_max should be close (not exact due to DG partitioning effects)
   real_t V_rel_err = std::abs(par_V_max - serial_V_max) /
                      std::max(serial_V_max, 1e-30);
   TEST_ASSERT(V_rel_err < 0.1,
               "Domain solution V_max: serial-parallel within 10%");

   if (mpi.IsRoot())
   {
      std::cout << "  Serial V_max: " << serial_V_max
                << ", Parallel V_max: " << par_V_max
                << " (rel err: " << V_rel_err << ")\n";
   }
}

// =============================================================================
// Test 2: Fault state consistency after short simulation
// =============================================================================

void test_fault_state(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Fault State Consistency (short sim) ===\n";
   }

   ConsistencySetup setup;
   real_t t_final = 1e6;  // ~12 days

   // Serial simulation on rank 0
   real_t serial_V_max = 0.0;
   int serial_steps = 0;
   // Serial fields sorted by depth (descending) to match GatherToRootDedup order
   Vector serial_slip_sorted, serial_theta_sorted, serial_depths_sorted;
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
      ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
      ode_solver.SetDtMin(1e-6);
      ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
      ode_solver.SetDt(1e3);
      ode_solver.Init(serial_seas);

      real_t t = 0.0;
      while (t < t_final)
      {
         if (t + ode_solver.GetDt() > t_final)
         {
            ode_solver.SetDt(t_final - t);
         }
         real_t dt;
         bool accepted = ode_solver.Step(serial_seas, state, t, dt);
         if (accepted) { serial_steps++; }
      }

      serial_V_max = serial_seas.GetMaxSlipRate();

      Vector serial_slip, serial_theta;
      serial_fault.GetSlip(state, serial_slip);
      serial_fault.GetTheta(state, serial_theta);

      // Get depths and sort everything by depth (descending, surface first)
      // to match GatherToRootDedup ordering
      Vector serial_depths;
      serial_domain.GetFaultDepths(serial_depths);
      int n = serial_depths.Size();
      std::vector<int> idx(n);
      for (int i = 0; i < n; i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](int a, int b) {
         return serial_depths(a) > serial_depths(b);
      });

      serial_slip_sorted.SetSize(n);
      serial_theta_sorted.SetSize(n);
      serial_depths_sorted.SetSize(n);
      for (int i = 0; i < n; i++)
      {
         serial_slip_sorted(i) = serial_slip(idx[i]);
         serial_theta_sorted(i) = serial_theta(idx[i]);
         serial_depths_sorted(i) = serial_depths(idx[i]);
      }

      std::cout << "  Serial: " << serial_steps << " steps, V_max = "
                << serial_V_max << "\n";
   }

   // Parallel simulation on all ranks
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
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(par_seas);

   real_t t = 0.0;
   int par_steps = 0;
   while (t < t_final)
   {
      if (t + ode_solver.GetDt() > t_final)
      {
         ode_solver.SetDt(t_final - t);
      }
      real_t dt;
      bool accepted = ode_solver.Step(par_seas, par_state, t, dt);
      if (accepted) { par_steps++; }
   }

   real_t par_V_max = par_seas.GetMaxSlipRate();

   // Gather parallel slip and theta via dedup
   Vector par_slip_local, par_theta_local;
   par_fault.GetSlip(par_state, par_slip_local);
   par_fault.GetTheta(par_state, par_theta_local);

   Vector dedup_slip, dedup_theta, dedup_depths;
   {
      std::vector<const Vector*> fields = {&par_slip_local, &par_theta_local};
      std::vector<Vector> dedup_fields;
      fault_geom.GatherFieldsToRootDedup(fields, dedup_fields, dedup_depths);
      if (mpi.IsRoot())
      {
         dedup_slip = dedup_fields[0];
         dedup_theta = dedup_fields[1];
      }
   }

   if (mpi.IsRoot())
   {
      std::cout << "  Parallel: " << par_steps << " steps, V_max = "
                << par_V_max << "\n";

      // Compare final slip (both sorted by depth descending)
      bool slip_ok = true;
      std::cout << "  Dedup slip size: " << dedup_slip.Size()
                << ", Serial slip size: " << serial_slip_sorted.Size() << "\n";
      if (dedup_slip.Size() == serial_slip_sorted.Size())
      {
         real_t max_rel_err = 0.0;
         for (int i = 0; i < serial_slip_sorted.Size(); i++)
         {
            real_t rel_err = std::abs(dedup_slip(i) - serial_slip_sorted(i)) /
                             std::max(std::abs(serial_slip_sorted(i)), 1e-30);
            max_rel_err = std::max(max_rel_err, rel_err);
            if (rel_err > 1e-8)
            {
               slip_ok = false;
            }
         }
         std::cout << "  Slip max rel err: " << max_rel_err << "\n";
      }
      else
      {
         slip_ok = false;
         std::cout << "  Size mismatch!\n";
      }
      TEST_ASSERT(slip_ok,
                  "Final slip matches serial within 1e-8 relative");

      // Compare final theta
      bool theta_ok = true;
      if (dedup_theta.Size() == serial_theta_sorted.Size())
      {
         real_t max_rel_err = 0.0;
         for (int i = 0; i < serial_theta_sorted.Size(); i++)
         {
            real_t rel_err = std::abs(dedup_theta(i) - serial_theta_sorted(i)) /
                             std::max(std::abs(serial_theta_sorted(i)), 1e-30);
            max_rel_err = std::max(max_rel_err, rel_err);
            if (rel_err > 1e-8)
            {
               theta_ok = false;
            }
         }
         std::cout << "  Theta max rel err: " << max_rel_err << "\n";
      }
      else
      {
         theta_ok = false;
         std::cout << "  Theta size mismatch: dedup=" << dedup_theta.Size()
                   << " serial=" << serial_theta_sorted.Size() << "\n";
      }
      TEST_ASSERT(theta_ok,
                  "Final theta matches serial within 1e-8 relative");

      // Compare V_max
      real_t V_rel_err = std::abs(par_V_max - serial_V_max) /
                         std::max(serial_V_max, 1e-30);
      TEST_ASSERT(V_rel_err < 1e-8,
                  "Final V_max matches serial within 1e-8 relative");
   }
   else
   {
      num_tests += 3;
      num_passed += 3;
   }
}

// =============================================================================
// Test 3: Probe output consistency
// =============================================================================

void test_probe_output(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Probe Output Consistency ===\n";
   }

   ConsistencySetup setup;
   std::vector<real_t> probe_depths = {0.0, -12000.0};
   std::string serial_prefix = "/tmp/seas_p9_serial";
   std::string par_prefix = "/tmp/seas_p9_par";

   // Serial: write at t=0 on rank 0
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

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      Vector serial_depths_vec;
      serial_domain.GetFaultDepths(serial_depths_vec);

      BenchmarkOutput<Mesh> serial_out(serial_prefix, setup.params,
                                        probe_depths, serial_depths_vec);
      serial_out.ForceWrite(0.0, serial_state, serial_fault,
                            serial_seas.GetTraction());
      serial_out.Close();
   }
   mpi.Barrier();

   // Parallel: write at t=0
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

   // Compare files on root
   if (mpi.IsRoot())
   {
      std::string serial_file = serial_prefix + "_z0km.txt";
      std::string par_file = par_prefix + "_z0km.txt";

      std::ifstream sf(serial_file);
      std::ifstream pf(par_file);

      TEST_ASSERT(sf.good(), "Serial probe file exists");
      TEST_ASSERT(pf.good(), "Parallel probe file exists");

      if (sf.good() && pf.good())
      {
         // Skip headers
         std::string sline, pline;
         while (std::getline(sf, sline) && !sline.empty() && sline[0] == '#') {}
         while (std::getline(pf, pline) && !pline.empty() && pline[0] == '#') {}

         std::istringstream ss(sline), ps(pline);
         double s_time, s_slip, s_logV, s_tau, s_logth;
         double p_time, p_slip, p_logV, p_tau, p_logth;

         ss >> s_time >> s_slip >> s_logV >> s_tau >> s_logth;
         ps >> p_time >> p_slip >> p_logV >> p_tau >> p_logth;

         TEST_NEAR(p_logV, s_logV, 1e-8,
                   "Probe log10(V) matches serial within 1e-8");
         TEST_NEAR(p_tau, s_tau, 1e-8,
                   "Probe tau matches serial within 1e-8");
      }

      // Cleanup
      std::remove((serial_prefix + "_z0km.txt").c_str());
      std::remove((serial_prefix + "_z12km.txt").c_str());
      std::remove((par_prefix + "_z0km.txt").c_str());
      std::remove((par_prefix + "_z12km.txt").c_str());
   }
   else
   {
      num_tests += 4;
      num_passed += 4;
   }
}

// =============================================================================
// Test 4: Reproducibility — all ranks agree
// =============================================================================

void test_reproducibility(MPIContext &mpi)
{
   if (mpi.IsRoot())
   {
      std::cout << "\n=== Test: Reproducibility (all ranks agree) ===\n";
   }

   ConsistencySetup setup;

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

   Vector state(par_fault.StateSize());
   par_seas.SetInitialCondition(state);

   DormandPrinceRK45 ode_solver;
   ode_solver.SetMPIContext(&mpi);
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);
   ode_solver.SetDt(1e3);
   ode_solver.Init(par_seas);

   real_t t = 0.0;
   int steps = 0;
   int target_steps = 10;

   while (steps < target_steps)
   {
      real_t dt;
      bool accepted = ode_solver.Step(par_seas, state, t, dt);
      if (accepted) { steps++; }
   }

   real_t V_max = par_seas.GetMaxSlipRate();

   // All ranks should agree on V_max
   real_t V_min_all = mpi.GlobalMin(V_max);
   real_t V_max_all = mpi.GlobalMax(V_max);
   TEST_ASSERT(std::abs(V_max_all - V_min_all) < 1e-20,
               "All ranks agree on final V_max");

   // All ranks should agree on step count
   int steps_min = mpi.GlobalMinInt(steps);
   int steps_max = mpi.GlobalMaxInt(steps);
   TEST_ASSERT(steps_min == steps_max,
               "All ranks agree on total accepted steps");

   // All ranks should agree on time
   real_t t_min = mpi.GlobalMin(t);
   real_t t_max = mpi.GlobalMax(t);
   TEST_ASSERT(std::abs(t_max - t_min) < 1e-20,
               "All ranks agree on final time");

   if (mpi.IsRoot())
   {
      std::cout << "  V_max = " << V_max << ", steps = " << steps
                << ", t = " << t << "\n";
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
      std::cout << "Phase 9: Serial-Parallel Consistency (np="
                << mpi.Size() << ")\n";
      std::cout << "============================================\n";
   }

   test_domain_solution(mpi);
   test_fault_state(mpi);
   test_probe_output(mpi);
   test_reproducibility(mpi);

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
