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
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>

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

/// Test 3: Owned fault DOFs match serial unique count and shared ghosts sync
bool test_owned_fault_layout(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_owned_fault_layout... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<Mesh> serial_op(serial_mesh, 1, lambda, mu,
                                             Vp, Wf, lf);
   const int nf_serial = serial_op.GetNumFaultDOFs();

   ParMesh pmesh(ctx.GetComm(), serial_mesh);
   ElasticityDomainOperator<ParMesh> par_op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   const int local_owned = par_op.GetNumOwnedFaultDOFs();
   const int global_owned = ctx.GlobalSumInt(local_owned);

   Vector owned_slip(2 * local_owned);
   for (int i = 0; i < local_owned; i++)
   {
      owned_slip(2 * i) = -1000.0 - 10.0 * ctx.Rank() - i;
      owned_slip(2 * i + 1) = 1000.0 + 10.0 * ctx.Rank() + i;
   }

   Vector local_slip;
   par_op.ExpandOwnedToLocalFault(owned_slip, local_slip, 2);

   Array<HYPRE_BigInt> global_face_ids;
   pmesh.GetGlobalFaceIndices(global_face_ids);

   std::vector<HYPRE_BigInt> local_gids;
   std::vector<real_t> local_vals;
   for (int i = 0; i < par_op.GetFaultSharedFaces().Size(); i++)
   {
      const int sf = par_op.GetFaultSharedFaces()[i];
      const int lf_idx = pmesh.GetSharedFace(sf);
      const int dof = par_op.GetFaultInteriorFaces().Size() + i;
      local_gids.push_back(global_face_ids[lf_idx]);
      local_vals.push_back(local_slip(2 * dof));
      local_vals.push_back(local_slip(2 * dof + 1));
   }

   int local_shared = static_cast<int>(local_gids.size());
   std::vector<int> recv_counts(ctx.Size(), 0);
#ifdef SEAS_USE_MPI
   MPI_Gather(&local_shared, 1, MPI_INT,
              ctx.IsRoot() ? recv_counts.data() : nullptr,
              1, MPI_INT, 0, ctx.GetComm());
#endif

   bool shared_ok = true;
   int total_shared = 0;
   std::vector<int> displs(ctx.Size(), 0);
   if (ctx.IsRoot())
   {
      for (int r = 0; r < ctx.Size(); r++)
      {
         displs[r] = total_shared;
         total_shared += recv_counts[r];
      }
   }

   std::vector<HYPRE_BigInt> global_gids(total_shared);
   std::vector<real_t> global_vals(2 * total_shared);
#ifdef SEAS_USE_MPI
   MPI_Gatherv(local_gids.data(), local_shared, HYPRE_MPI_BIG_INT,
               ctx.IsRoot() ? global_gids.data() : nullptr,
               ctx.IsRoot() ? recv_counts.data() : nullptr,
               ctx.IsRoot() ? displs.data() : nullptr,
               HYPRE_MPI_BIG_INT, 0, ctx.GetComm());

   std::vector<int> recv_counts_vals(ctx.Size(), 0), displs_vals(ctx.Size(), 0);
   if (ctx.IsRoot())
   {
      for (int r = 0; r < ctx.Size(); r++)
      {
         recv_counts_vals[r] = 2 * recv_counts[r];
         displs_vals[r] = 2 * displs[r];
      }
   }
   MPI_Gatherv(local_vals.data(), 2 * local_shared, MPI_DOUBLE,
               ctx.IsRoot() ? global_vals.data() : nullptr,
               ctx.IsRoot() ? recv_counts_vals.data() : nullptr,
               ctx.IsRoot() ? displs_vals.data() : nullptr,
               MPI_DOUBLE, 0, ctx.GetComm());
#endif

   if (ctx.IsRoot())
   {
      std::map<HYPRE_BigInt, std::pair<real_t, real_t>> first_seen;
      for (int i = 0; i < total_shared; i++)
      {
         const auto gid = global_gids[i];
         const std::pair<real_t, real_t> val =
            {global_vals[2 * i], global_vals[2 * i + 1]};
         auto it = first_seen.find(gid);
         if (it == first_seen.end())
         {
            first_seen.emplace(gid, val);
         }
         else if (std::abs(it->second.first - val.first) > 1e-12 ||
                  std::abs(it->second.second - val.second) > 1e-12)
         {
            shared_ok = false;
            break;
         }
      }
   }

   int shared_ok_int = shared_ok ? 1 : 0;
   shared_ok_int = ctx.GlobalMinInt(shared_ok_int);

   const bool ok = (global_owned == nf_serial) && (shared_ok_int == 1);

   TEST_CHECK(ctx, "owned BP5 fault layout follows serial unique faces", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (global_owned=" << global_owned
                << ", serial=" << nf_serial << ")"
                << std::endl;
   }
   return ok;
}

/// Test 4: Canonical DOF permutation — parallel owned coords match serial
///
/// The strongest invariant: serial fault DOF coordinates must exactly match
/// the parallel owned coordinates (which pass through the canonical
/// permutation in RestrictToOwnedFault).  Gather all owned (x2, x3) pairs
/// to root, sort both sets lexicographically, and compare.
bool test_canonical_dof_coords_match_serial(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_canonical_dof_coords_match_serial... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   // --- Serial reference coordinates ---
   ElasticityDomainOperator<Mesh> serial_op(serial_mesh, 1, lambda, mu,
                                             Vp, Wf, lf);
   Vector serial_x2, serial_x3;
   serial_op.GetFaultCoords2D(serial_x2, serial_x3);
   const int nf_serial = serial_op.GetNumFaultDOFs();

   // --- Parallel owned coordinates (through canonical permutation) ---
   ParMesh pmesh(ctx.GetComm(), serial_mesh);
   ElasticityDomainOperator<ParMesh> par_op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   Vector local_x2_full, local_x3_full;
   par_op.GetFaultCoords2D(local_x2_full, local_x3_full);

   // RestrictToOwnedFault applies the canonical permutation
   Vector owned_x2, owned_x3;
   par_op.RestrictToOwnedFault(local_x2_full, owned_x2);
   par_op.RestrictToOwnedFault(local_x3_full, owned_x3);
   const int local_owned = par_op.GetNumOwnedFaultDOFs();

   // --- Gather owned coordinates to root ---
   std::vector<int> recv_counts(ctx.Size(), 0);
   std::vector<int> displs(ctx.Size(), 0);
#ifdef SEAS_USE_MPI
   MPI_Gather(&local_owned, 1, MPI_INT,
              ctx.IsRoot() ? recv_counts.data() : nullptr,
              1, MPI_INT, 0, ctx.GetComm());
#endif

   int total_owned = 0;
   if (ctx.IsRoot())
   {
      for (int r = 0; r < ctx.Size(); r++)
      {
         displs[r] = total_owned;
         total_owned += recv_counts[r];
      }
   }

   std::vector<real_t> all_x2(total_owned), all_x3(total_owned);
#ifdef SEAS_USE_MPI
   MPI_Gatherv(owned_x2.GetData(), local_owned, MPI_DOUBLE,
               ctx.IsRoot() ? all_x2.data() : nullptr,
               ctx.IsRoot() ? recv_counts.data() : nullptr,
               ctx.IsRoot() ? displs.data() : nullptr,
               MPI_DOUBLE, 0, ctx.GetComm());
   MPI_Gatherv(owned_x3.GetData(), local_owned, MPI_DOUBLE,
               ctx.IsRoot() ? all_x3.data() : nullptr,
               ctx.IsRoot() ? recv_counts.data() : nullptr,
               ctx.IsRoot() ? displs.data() : nullptr,
               MPI_DOUBLE, 0, ctx.GetComm());
#endif

   // --- Compare on root: sort both sets by (x2, x3) and check ---
   bool coords_ok = true;
   if (ctx.IsRoot())
   {
      // Build and sort serial coordinate pairs
      std::vector<std::pair<real_t, real_t>> serial_pairs(nf_serial);
      for (int i = 0; i < nf_serial; i++)
      {
         serial_pairs[i] = {serial_x2(i), serial_x3(i)};
      }
      std::sort(serial_pairs.begin(), serial_pairs.end());

      // Build and sort parallel coordinate pairs
      std::vector<std::pair<real_t, real_t>> par_pairs(total_owned);
      for (int i = 0; i < total_owned; i++)
      {
         par_pairs[i] = {all_x2[i], all_x3[i]};
      }
      std::sort(par_pairs.begin(), par_pairs.end());

      if (nf_serial != total_owned)
      {
         coords_ok = false;
      }
      else
      {
         for (int i = 0; i < nf_serial; i++)
         {
            if (std::abs(serial_pairs[i].first - par_pairs[i].first) > 1e-10 ||
                std::abs(serial_pairs[i].second - par_pairs[i].second) > 1e-10)
            {
               coords_ok = false;
               if (i < 5) // print first few mismatches
               {
                  std::cerr << "  mismatch DOF " << i
                            << ": serial=(" << serial_pairs[i].first
                            << "," << serial_pairs[i].second
                            << ") par=(" << par_pairs[i].first
                            << "," << par_pairs[i].second << ")\n";
               }
            }
         }
      }
   }

   int ok_int = coords_ok ? 1 : 0;
   ok_int = ctx.GlobalMinInt(ok_int);
   const bool ok = (ok_int == 1);

   TEST_CHECK(ctx, "canonical DOF coords match serial", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (serial=" << nf_serial
                << ", par_owned=" << total_owned << ")"
                << std::endl;
   }
   return ok;
}

/// Test 5: Serial-parallel displacement match for non-uniform slip
///
/// Core proof-of-bug for shared-face parameterization: assemble K and b
/// in serial and parallel with the SAME non-uniform slip (step function
/// at the fault midpoint, mimicking nucleation boundary). Compare the
/// global ||u||_inf.  Any difference > O(ε) proves that shared face
/// parameterization corrupts the RHS.
bool test_serial_parallel_displacement_match(MPIContext &ctx)
{
   if (ctx.IsRoot())
   {
      std::cout << "  test_serial_parallel_displacement_match... " << std::flush;
   }

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   auto serial_mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;
   real_t lf = 2.0 * Ly;

   // --- Serial: solve with non-uniform strike slip ---
   ElasticityDomainOperator<Mesh> serial_op(serial_mesh, 1, lambda, mu,
                                             Vp, Wf, lf);
   int nf_serial = serial_op.GetNumFaultDOFs();
   Vector serial_x2, serial_x3;
   serial_op.GetFaultCoords2D(serial_x2, serial_x3);

   // Non-uniform slip: step function at x2=0 (mimics nucleation boundary)
   Vector slip_serial(2 * nf_serial);
   slip_serial = 0.0;
   for (int i = 0; i < nf_serial; i++)
   {
      real_t x2 = serial_x2(i);
      slip_serial(2 * i + 1) = (x2 > 0.0) ? 1.0 : 0.01;  // strike slip
   }

   GridFunction u_serial(&serial_op.GetFESpace());
   serial_op.Solve(0.0, slip_serial, u_serial);
   real_t serial_u_inf = u_serial.Normlinf();

   // --- Parallel: solve with same non-uniform slip ---
   ParMesh pmesh(ctx.GetComm(), serial_mesh);
   ElasticityDomainOperator<ParMesh> par_op(pmesh, 1, lambda, mu, Vp, Wf, lf);

   int nf_par = par_op.GetNumFaultDOFs();
   Vector par_x2, par_x3;
   par_op.GetFaultCoords2D(par_x2, par_x3);

   Vector slip_par(2 * nf_par);
   slip_par = 0.0;
   for (int i = 0; i < nf_par; i++)
   {
      real_t x2 = par_x2(i);
      slip_par(2 * i + 1) = (x2 > 0.0) ? 1.0 : 0.01;  // same step function
   }

   ParGridFunction u_par(&par_op.GetFESpace());
   par_op.Solve(0.0, slip_par, u_par);
   real_t local_u_inf = u_par.Normlinf();
   real_t parallel_u_inf = ctx.GlobalMax(local_u_inf);

   // Broadcast serial result
   ctx.Bcast(serial_u_inf);

   real_t rel_err = (serial_u_inf > 0.0)
      ? std::abs(parallel_u_inf - serial_u_inf) / serial_u_inf : 0.0;

   // Tight tolerance: should be < 1e-10 for partition-independent assembly.
   // If > 1e-6, shared face parameterization is corrupting the solve.
   bool ok = (rel_err < 1e-10);

   TEST_CHECK(ctx, "serial-parallel displacement match (non-uniform slip)", ok);

   if (ctx.IsRoot())
   {
      std::cout << (ok ? "PASSED" : "FAILED")
                << " (serial_u=" << std::scientific << std::setprecision(12)
                << serial_u_inf
                << ", parallel_u=" << parallel_u_inf
                << ", rel_err=" << rel_err << ")"
                << std::endl;
   }
   return ok;
}

/// Test 6: Zero slip → zero traction in parallel
/// (unchanged from before)
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

/// Test 6: Serial-parallel traction consistency for uniform slip
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

/// Test 7: All parallel traction values are bounded (no blowup)
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

/// Test 8: Parallel solve with non-zero slip produces non-zero displacement
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
   test_owned_fault_layout(ctx);
   test_canonical_dof_coords_match_serial(ctx);
   test_serial_parallel_displacement_match(ctx);
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
