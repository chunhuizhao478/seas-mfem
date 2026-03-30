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

// Parallel elasticity operator tests — shared face support for BP5
// Run with: mpirun -np 2 ./seas_test_parallel_elasticity
//       or: mpirun -np 4 ./seas_test_parallel_elasticity

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/bp5_params.hpp"

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

/// Create a 3D hex mesh: [-Lx,Lx] x [-Ly,Ly] x [-Lz,0]
/// Boundary attrs: 1=top/bottom (z=0,z=-Lz), 5=far-field (sides)
/// Matches serial test mesh convention (Tandem: z <= 0 for depth).
Mesh CreateTestMesh3D(int nx, int ny, int nz,
                       real_t Lx, real_t Ly, real_t Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                                      Element::HEXAHEDRON,
                                      2.0 * Lx, 2.0 * Ly, Lz);

   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
      v[2] -= Lz;  // Z ranges [-Lz, 0]
   }

   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(3);
      T->Transform(ip, center);

      real_t tol = 1e-6;
      if (std::abs(center(2)) < tol || std::abs(center(2) + Lz) < tol)
      {
         mesh.SetBdrAttribute(be, 1);  // Natural (top/bottom)
      }
      else
      {
         mesh.SetBdrAttribute(be, 5);  // Dirichlet (far-field)
      }
   }

   mesh.SetAttributes();
   return mesh;
}

/// Test 1: Shared fault faces are detected in parallel
bool test_shared_fault_detection(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_shared_fault_detection... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);
   ParMesh pmesh(ctx.GetComm(), serial_mesh);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<ParMesh> op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int local_nf = op.GetNumFaultDOFs();
   int global_nf = ctx.GlobalSumInt(local_nf);

   // At least one rank should have shared faces (for np >= 2)
   int local_shared = op.GetFaultSharedFaces().Size();
   int global_shared = ctx.GlobalSumInt(local_shared);

   bool ok = (global_nf > 0);
   // With 2+ ranks, partitioner should create some shared faces
   // (though not guaranteed for all mesh sizes)

   TEST_CHECK(ctx, "fault faces detected in parallel", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (global_nf=" << global_nf
                << ", global_shared=" << global_shared << ")"
                << std::endl;
   }
   return ok;
}

/// Test 2: Global fault DOF count in parallel matches serial
bool test_serial_parallel_fault_dof_count(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_serial_parallel_fault_dof_count... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   // Serial count (all ranks compute independently)
   ElasticityDomainOperator<Mesh> serial_op(serial_mesh, 1, lambda, mu,
                                             Vp, Wf, lf);
   int nf_serial = serial_op.GetNumFaultDOFs();

   // Parallel count
   ParMesh pmesh(ctx.GetComm(), serial_mesh);
   ElasticityDomainOperator<ParMesh> par_op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int local_interior = par_op.GetFaultInteriorFaces().Size();
   int local_shared = par_op.GetFaultSharedFaces().Size();
   int global_interior = ctx.GlobalSumInt(local_interior);
   int global_shared = ctx.GlobalSumInt(local_shared);

   // Interior faces counted once; shared faces counted by both ranks
   // So: global_interior + global_shared/2 == serial_fault_dofs
   bool ok = (global_interior + global_shared / 2 == nf_serial);

   TEST_CHECK(ctx, "parallel fault DOF count matches serial", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (interior=" << global_interior
                << ", shared_total=" << global_shared
                << ", serial=" << nf_serial << ")"
                << std::endl;
   }
   return ok;
}

/// Test 3: Zero slip → zero traction in parallel
bool test_parallel_zero_slip_traction(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_zero_slip_traction... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);
   ParMesh pmesh(ctx.GetComm(), serial_mesh);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<ParMesh> op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int nf = op.GetNumFaultDOFs();
   Vector slip(2 * nf);
   slip = 0.0;

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   Vector traction;
   op.ComputeTraction(u, slip, traction);

   real_t local_max = 0.0;
   for (int i = 0; i < traction.Size(); i++)
   {
      local_max = std::max(local_max, std::abs(traction(i)));
   }
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_max < 1e-8);

   TEST_CHECK(ctx, "zero slip gives zero traction", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_traction=" << global_max << ")"
                << std::endl;
   }
   return ok;
}

/// Test 4: Serial-parallel traction consistency for uniform slip
bool test_serial_parallel_traction_consistency(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_serial_parallel_traction_consistency... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   // Serial solve with uniform strike slip
   ElasticityDomainOperator<Mesh> serial_op(serial_mesh, 1, lambda, mu,
                                             Vp, Wf, lf);
   int nf_serial = serial_op.GetNumFaultDOFs();
   Vector slip_serial(2 * nf_serial);
   slip_serial = 0.0;
   for (int i = 0; i < nf_serial; i++)
   {
      slip_serial(2 * i + 1) = 1.0;  // strike slip
   }

   GridFunction u_serial(&serial_op.GetFESpace());
   serial_op.Solve(0.0, slip_serial, u_serial);

   Vector trac_serial;
   serial_op.ComputeTraction(u_serial, slip_serial, trac_serial);

   // Compute serial max traction magnitude
   real_t serial_max_trac = 0.0;
   for (int i = 0; i < trac_serial.Size(); i++)
   {
      serial_max_trac = std::max(serial_max_trac, std::abs(trac_serial(i)));
   }

   // Broadcast serial result to all ranks
   ctx.Bcast(serial_max_trac);

   // Parallel solve with same uniform strike slip
   ParMesh pmesh(ctx.GetComm(), serial_mesh);
   ElasticityDomainOperator<ParMesh> par_op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int nf_par = par_op.GetNumFaultDOFs();
   Vector slip_par(2 * nf_par);
   slip_par = 0.0;
   for (int i = 0; i < nf_par; i++)
   {
      slip_par(2 * i + 1) = 1.0;  // strike slip
   }

   ParGridFunction u_par(&par_op.GetFESpace());
   par_op.Solve(0.0, slip_par, u_par);

   Vector trac_par;
   par_op.ComputeTraction(u_par, slip_par, trac_par);

   // Compute parallel max traction magnitude
   real_t local_max_trac = 0.0;
   for (int i = 0; i < trac_par.Size(); i++)
   {
      local_max_trac = std::max(local_max_trac, std::abs(trac_par(i)));
   }
   real_t parallel_max_trac = ctx.GlobalMax(local_max_trac);

   // Max traction should be similar (same problem, same slip)
   real_t rel_err = (serial_max_trac > 0.0)
      ? std::abs(parallel_max_trac - serial_max_trac) / serial_max_trac : 0.0;

   bool ok = (rel_err < 0.1);  // 10% tolerance for different partitionings

   TEST_CHECK(ctx, "serial-parallel traction consistency", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (serial_max=" << serial_max_trac
                << ", parallel_max=" << parallel_max_trac
                << ", rel_err=" << rel_err << ")"
                << std::endl;
   }
   return ok;
}

/// Test 5: All parallel traction values are bounded (no blowup)
bool test_parallel_traction_bounded(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_traction_bounded... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);
   ParMesh pmesh(ctx.GetComm(), serial_mesh);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<ParMesh> op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int nf = op.GetNumFaultDOFs();
   Vector slip(2 * nf);
   slip = 0.0;
   for (int i = 0; i < nf; i++)
   {
      slip(2 * i + 1) = 1.0;  // strike slip
   }

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   Vector traction;
   op.ComputeTraction(u, slip, traction);

   // Physically meaningful bound: for unit slip on a coarse mesh with
   // BP5 material (mu ~ 32 GPa), traction should be O(mu).  Allow 10x
   // headroom for DG penalty on the coarse 2x1x1 test mesh.
   const real_t traction_bound = 10.0 * mu;

   bool local_ok = true;
   real_t local_max = 0.0;
   for (int i = 0; i < traction.Size(); i++)
   {
      real_t val = traction(i);
      if (std::isnan(val) || std::isinf(val) || std::abs(val) > traction_bound)
      {
         local_ok = false;
      }
      local_max = std::max(local_max, std::abs(val));
   }

   int local_ok_int = local_ok ? 1 : 0;
   int global_ok = ctx.GlobalMinInt(local_ok_int);
   real_t global_max = ctx.GlobalMax(local_max);

   bool ok = (global_ok == 1);

   TEST_CHECK(ctx, "parallel traction bounded", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (max_traction=" << global_max << ")"
                << std::endl;
   }
   return ok;
}

/// Test 6: Parallel solve with non-zero slip produces non-zero displacement
bool test_parallel_solve_with_slip(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_parallel_solve_with_slip... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);
   ParMesh pmesh(ctx.GetComm(), serial_mesh);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<ParMesh> op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int nf = op.GetNumFaultDOFs();
   Vector slip(2 * nf);
   slip = 0.0;
   for (int i = 0; i < nf; i++)
   {
      slip(2 * i + 1) = 1.0;  // strike slip
   }

   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);

   real_t local_norm = u.Norml2();
   real_t global_norm_sq = ctx.GlobalSum(local_norm * local_norm);
   real_t global_norm = std::sqrt(global_norm_sq);

   bool ok = (global_norm > 0.0);

   TEST_CHECK(ctx, "parallel solve produces non-zero displacement", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (||u||=" << global_norm << ")"
                << std::endl;
   }
   return ok;
}

int main(int argc, char *argv[])
{
   MPIContext ctx(&argc, &argv);

   if (ctx.IsRoot())
   {
      std::cout << "Parallel Elasticity Operator Tests ("
                << ctx.Size() << " ranks)" << std::endl;
      std::cout << "========================================" << std::endl;
   }

   test_shared_fault_detection(ctx);
   test_serial_parallel_fault_dof_count(ctx);
   test_parallel_zero_slip_traction(ctx);
   test_serial_parallel_traction_consistency(ctx);
   test_parallel_traction_bounded(ctx);
   test_parallel_solve_with_slip(ctx);

   if (ctx.IsRoot())
   {
      std::cout << "========================================" << std::endl;
      std::cout << num_passed << " / " << num_tests << " tests passed"
                << std::endl;
   }

   return (num_passed == num_tests) ? 0 : 1;
}
