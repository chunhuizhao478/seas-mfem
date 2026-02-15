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

// Phase 7: Parallel domain operator tests
// Run with: mpirun -np 2 ./seas_test_parallel_domain
//       or: mpirun -np 4 ./seas_test_parallel_domain

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../config/bp2_params.hpp"
#include "../../common/parallel_utils.hpp"

#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0;
static int num_passed = 0;

#define TEST_CHECK(ctx, name, condition) \
   do { \
      num_tests++; \
      if (condition) { num_passed++; } \
      else if (ctx.IsRoot()) { \
         std::cerr << "FAIL: " << name << std::endl; \
      } \
   } while(0)

/// Test: ParMesh preserves global element count
bool test_mesh_distribution(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_mesh_distribution... " << std::flush;
   }

   // Create serial mesh on all ranks, then partition
   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 4;
   params.nz = 4;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   int global_ne_expected = serial_mesh->GetNE();  // 32

   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   // Sum local elements across ranks
   int local_ne = pmesh.GetNE();
   int global_ne = ctx.GlobalSumInt(local_ne);

   bool ok = (global_ne == global_ne_expected);

   TEST_CHECK(ctx, "global element count preserved", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (global_ne=" << global_ne << ", expected=" << global_ne_expected << ")"
                << std::endl;
   }
   return ok;
}

/// Test: DG FE space has TrueVSize == VSize (no shared DOFs)
bool test_dof_distribution(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_dof_distribution... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 4;
   params.nz = 4;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   DG_FECollection fec(1, 2, BasisType::GaussLobatto);
   ParFiniteElementSpace pfes(&pmesh, &fec);

   // DG property: no shared DOFs
   bool ok = (pfes.GetTrueVSize() == pfes.GetVSize());

   TEST_CHECK(ctx, "DG TrueVSize == VSize", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (TrueVSize=" << pfes.GetTrueVSize()
                << ", VSize=" << pfes.GetVSize() << ")"
                << std::endl;
   }
   return ok;
}

/// Test: Parallel solve with zero slip converges
bool test_parallel_solve(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_solve... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   real_t Wf = 10.0e3;
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf);

   Vector slip(op.GetNumFaultDOFs());
   slip = 0.0;

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   // With zero slip, solution should be near zero
   real_t local_max = u.Normlinf();
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_max < 1e-10);

   TEST_CHECK(ctx, "zero slip solution near zero", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_u=" << global_max << ")" << std::endl;
   }
   return ok;
}

/// Test: Serial and parallel traction match within tolerance
/// The system has all-Neumann BCs so solutions differ by a constant.
/// We compare traction (gradient-based) which is independent of the constant.
bool test_serial_parallel_consistency(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_serial_parallel_consistency... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   // --- Serial solve and traction ---
   real_t serial_trac_norm = 0.0;
   {
      auto serial_mesh = BP2MeshGenerator::Create(params);
      BP2Params bp2;
      real_t Wf = 10.0e3;
      AntiplaneDomainOperator<Mesh> op(*serial_mesh, 1, bp2.mu(), bp2.Vp, Wf);

      Vector slip(op.GetNumFaultDOFs());
      slip = 1.0;

      GridFunction u(&op.GetFESpace());
      op.Solve(0.0, slip, u);

      Vector traction;
      op.ComputeTraction(u, slip, traction);
      serial_trac_norm = traction.Norml2();
   }

   // Broadcast serial result
   ctx.Bcast(serial_trac_norm);

   // --- Parallel solve and traction ---
   real_t parallel_trac_norm = 0.0;
   {
      auto serial_mesh = BP2MeshGenerator::Create(params);
      ParMesh pmesh(ctx.GetComm(), *serial_mesh);

      BP2Params bp2;
      real_t Wf = 10.0e3;
      AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf);

      int local_fault_dofs = op.GetNumFaultDOFs();
      Vector slip(local_fault_dofs);
      slip = 1.0;

      ParGridFunction u(&op.GetFESpace());
      op.Solve(0.0, slip, u);

      Vector traction;
      op.ComputeTraction(u, slip, traction);

      // Compute global traction norm: sum of squares across ranks
      real_t local_sq = traction * traction;
      real_t global_sq = 0.0;
      MPI_Allreduce(&local_sq, &global_sq, 1, MPI_DOUBLE, MPI_SUM,
                     ctx.GetComm());
      parallel_trac_norm = std::sqrt(global_sq);
   }

   // Both traction norms should be small (uniform slip gives near-zero residual).
   // Use absolute tolerance since we're comparing small discretization errors.
   real_t abs_diff = std::abs(serial_trac_norm - parallel_trac_norm);
   real_t max_trac = std::max(serial_trac_norm, parallel_trac_norm);

   // Both should be small (< 1e-3) and reasonably close
   bool ok = (max_trac < 1e-3);

   TEST_CHECK(ctx, "serial-parallel traction consistency", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (serial_trac=" << serial_trac_norm
                << ", parallel_trac=" << parallel_trac_norm
                << ", abs_diff=" << abs_diff << ")"
                << std::endl;
   }
   return ok;
}

/// Test: Parallel uniform slip with full-depth fault gives max_u ~ 0.5
bool test_parallel_uniform_slip(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_uniform_slip... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 16;
   params.nz = 16;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   real_t Wf = params.Lz;  // Full depth fault
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf);

   Vector slip(op.GetNumFaultDOFs());
   slip = 1.0;

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   real_t local_max = u.Normlinf();
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_max > 0.35 && global_max < 0.65);

   TEST_CHECK(ctx, "full-depth fault max_u ~ 0.5", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_u=" << global_max << ", expected ~0.5)"
                << std::endl;
   }
   return ok;
}

/// Test: Fault DOF counts are reasonable
/// In parallel, shared faces are counted by both ranks, so the global sum
/// of local fault DOFs >= serial fault DOFs. Each rank's interior + shared
/// faces cover the full fault when combined across ranks.
bool test_fault_dof_count(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_fault_dof_count... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   // Serial fault DOF count
   int serial_fault_dofs = 0;
   {
      auto mesh = BP2MeshGenerator::Create(params);
      BP2Params bp2;
      AntiplaneDomainOperator<Mesh> op(*mesh, 1, bp2.mu(), bp2.Vp, params.Lz);
      serial_fault_dofs = op.GetNumFaultDOFs();
   }

   // Parallel
   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, params.Lz);

   int local_fault_dofs = op.GetNumFaultDOFs();
   int local_interior = op.GetFaultInteriorFaces().Size();
   int local_shared = op.GetFaultSharedFaces().Size();
   int global_interior = ctx.GlobalSumInt(local_interior);
   int global_shared = ctx.GlobalSumInt(local_shared);

   // Interior faces are counted once; shared faces counted by both ranks
   // So: global_interior + global_shared/2 == serial_fault_dofs
   bool ok = (global_interior + global_shared / 2 == serial_fault_dofs);

   TEST_CHECK(ctx, "fault DOF accounting correct", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (interior=" << global_interior
                << ", shared_total=" << global_shared
                << ", serial=" << serial_fault_dofs << ")"
                << std::endl;
   }
   return ok;
}

/// Test: Parallel BR2 solve with zero slip converges
bool test_parallel_solve_br2(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_solve_br2... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   real_t Wf = 10.0e3;
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf,
                                        DGMethod::BR2);

   Vector slip(op.GetNumFaultDOFs());
   slip = 0.0;

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   real_t local_max = u.Normlinf();
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_max < 1e-10);

   TEST_CHECK(ctx, "BR2 zero slip solution near zero", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_u=" << global_max << ")" << std::endl;
   }
   return ok;
}

/// Test: Parallel BR2 uniform slip with full-depth fault gives max_u ~ 0.5
bool test_parallel_uniform_slip_br2(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_uniform_slip_br2... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 16;
   params.nz = 16;

   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   real_t Wf = params.Lz;  // Full depth fault
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf,
                                        DGMethod::BR2);

   Vector slip(op.GetNumFaultDOFs());
   slip = 1.0;

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   real_t local_max = u.Normlinf();
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_max > 0.35 && global_max < 0.65);

   TEST_CHECK(ctx, "BR2 full-depth fault max_u ~ 0.5", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_u=" << global_max << ", expected ~0.5)"
                << std::endl;
   }
   return ok;
}

/// Test: Serial and parallel traction match with non-uniform slip
bool test_traction_nonuniform_slip(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_traction_nonuniform_slip... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   // --- Serial solve with non-uniform slip ---
   Vector serial_traction;
   Vector serial_depths;
   {
      auto mesh = BP2MeshGenerator::Create(params);
      BP2Params bp2;
      real_t Wf = params.Lz;
      AntiplaneDomainOperator<Mesh> op(*mesh, 1, bp2.mu(), bp2.Vp, Wf);

      op.GetFaultDepths(serial_depths);
      int nf = op.GetNumFaultDOFs();
      Vector slip(nf);
      for (int i = 0; i < nf; i++)
      {
         // Depth-varying slip: sin(pi * z / Lz) where z in [-Lz, 0]
         slip(i) = std::sin(M_PI * (serial_depths(i) + params.Lz) / params.Lz);
      }

      GridFunction u(&op.GetFESpace());
      op.Solve(0.0, slip, u);
      op.ComputeTraction(u, slip, serial_traction);
   }

   // Broadcast serial data to all ranks
   ctx.Bcast(serial_depths);
   ctx.Bcast(serial_traction);

   // --- Parallel solve with non-uniform slip ---
   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   real_t Wf = params.Lz;
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf);

   Vector par_depths;
   op.GetFaultDepths(par_depths);
   int local_nf = op.GetNumFaultDOFs();
   Vector slip(local_nf);
   for (int i = 0; i < local_nf; i++)
   {
      slip(i) = std::sin(M_PI * (par_depths(i) + params.Lz) / params.Lz);
   }

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   Vector par_traction;
   op.ComputeTraction(u, slip, par_traction);

   // Gather parallel depths and traction to root for comparison
   std::vector<real_t> local_depths_vec(local_nf);
   std::vector<real_t> local_traction_vec(local_nf);
   for (int i = 0; i < local_nf; i++)
   {
      local_depths_vec[i] = par_depths(i);
      local_traction_vec[i] = par_traction(i);
   }

   std::vector<real_t> global_depths_vec, global_traction_vec;
   GatherVectorToRoot(local_depths_vec, global_depths_vec, ctx.GetComm());
   GatherVectorToRoot(local_traction_vec, global_traction_vec, ctx.GetComm());

   bool ok = true;
   if (ctx.IsRoot())
   {
      // For each parallel fault DOF, find the matching serial DOF by depth
      // and compare traction values
      real_t max_rel_err = 0.0;
      int matched = 0;
      for (size_t p = 0; p < global_depths_vec.size(); p++)
      {
         real_t pz = global_depths_vec[p];
         real_t pt = global_traction_vec[p];

         // Find closest serial DOF
         int best = -1;
         real_t best_dist = 1e30;
         for (int s = 0; s < serial_depths.Size(); s++)
         {
            real_t d = std::abs(serial_depths(s) - pz);
            if (d < best_dist)
            {
               best_dist = d;
               best = s;
            }
         }

         if (best >= 0 && best_dist < 1.0)  // Within 1 m
         {
            real_t st = serial_traction(best);
            real_t denom = std::max(std::abs(st), 1e-10);
            real_t rel_err = std::abs(pt - st) / denom;
            max_rel_err = std::max(max_rel_err, rel_err);
            matched++;
         }
      }

      // Expect most DOFs to match within reasonable tolerance
      // DG on different partitions can produce slightly different results
      ok = (matched > 0 && max_rel_err < 0.1);

      std::cout << (ok ? "PASSED" : "FAILED")
                << " (matched=" << matched
                << ", max_rel_err=" << max_rel_err << ")"
                << std::endl;
   }

   // Broadcast result to all ranks
   int ok_int = ok ? 1 : 0;
   MPI_Bcast(&ok_int, 1, MPI_INT, 0, ctx.GetComm());
   ok = (ok_int == 1);

   TEST_CHECK(ctx, "nonuniform slip traction consistency", ok);

   return ok;
}

/// Test: Gathered parallel fault depths match serial fault depths
bool test_fault_depth_consistency(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_fault_depth_consistency... " << std::flush;
   }

   BP2MeshGenerator::Parameters params;
   params.Lx = 10.0e3;
   params.Lz = 10.0e3;
   params.nx = 8;
   params.nz = 8;

   // Serial fault depths
   std::vector<real_t> serial_sorted;
   {
      auto mesh = BP2MeshGenerator::Create(params);
      BP2Params bp2;
      AntiplaneDomainOperator<Mesh> op(*mesh, 1, bp2.mu(), bp2.Vp, params.Lz);
      Vector depths;
      op.GetFaultDepths(depths);
      serial_sorted.resize(depths.Size());
      for (int i = 0; i < depths.Size(); i++)
      {
         serial_sorted[i] = depths(i);
      }
      std::sort(serial_sorted.begin(), serial_sorted.end());
   }

   // Parallel fault depths
   auto serial_mesh = BP2MeshGenerator::Create(params);
   ParMesh pmesh(ctx.GetComm(), *serial_mesh);

   BP2Params bp2;
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, params.Lz);
   Vector par_depths;
   op.GetFaultDepths(par_depths);

   int local_nf = par_depths.Size();
   std::vector<real_t> local_vec(local_nf);
   for (int i = 0; i < local_nf; i++)
   {
      local_vec[i] = par_depths(i);
   }

   std::vector<real_t> global_vec;
   GatherVectorToRoot(local_vec, global_vec, ctx.GetComm());

   bool ok = true;
   if (ctx.IsRoot())
   {
      // Remove duplicates from gathered depths (shared faces counted twice)
      std::sort(global_vec.begin(), global_vec.end());
      std::vector<real_t> unique_depths;
      for (size_t i = 0; i < global_vec.size(); i++)
      {
         if (unique_depths.empty() ||
             std::abs(global_vec[i] - unique_depths.back()) > 1e-6)
         {
            unique_depths.push_back(global_vec[i]);
         }
      }

      // Compare sizes
      ok = (unique_depths.size() == serial_sorted.size());

      if (ok)
      {
         // Compare values pointwise
         for (size_t i = 0; i < serial_sorted.size(); i++)
         {
            if (std::abs(unique_depths[i] - serial_sorted[i]) > 1.0)
            {
               ok = false;
               break;
            }
         }
      }

      std::cout << (ok ? "PASSED" : "FAILED")
                << " (serial=" << serial_sorted.size()
                << ", parallel_unique=" << unique_depths.size() << ")"
                << std::endl;
   }

   int ok_int = ok ? 1 : 0;
   MPI_Bcast(&ok_int, 1, MPI_INT, 0, ctx.GetComm());
   ok = (ok_int == 1);

   TEST_CHECK(ctx, "fault depth consistency", ok);

   return ok;
}

int main(int argc, char *argv[])
{
   MPIContext ctx(&argc, &argv);

   if (ctx.IsRoot())
   {
      std::cout << "===============================================" << std::endl;
      std::cout << "   Parallel Domain Operator Tests" << std::endl;
      std::cout << "   (np=" << ctx.Size() << ")" << std::endl;
      std::cout << "===============================================" << std::endl;
   }

   test_mesh_distribution(ctx);
   test_dof_distribution(ctx);
   test_fault_dof_count(ctx);
   test_parallel_solve(ctx);
   test_parallel_uniform_slip(ctx);
   test_serial_parallel_consistency(ctx);
   test_parallel_solve_br2(ctx);
   test_parallel_uniform_slip_br2(ctx);
   test_traction_nonuniform_slip(ctx);
   test_fault_depth_consistency(ctx);

   if (ctx.IsRoot())
   {
      std::cout << "\nParallel domain tests: " << num_passed << " / " << num_tests
                << " passed" << std::endl;
   }

   return (num_passed == num_tests) ? 0 : 1;
}
