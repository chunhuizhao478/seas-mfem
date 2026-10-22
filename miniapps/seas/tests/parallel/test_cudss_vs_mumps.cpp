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

// test_cudss_vs_mumps.cpp — GPU plan G3a acceptance: cuDSS against MUMPS on
// the BP5 stiffness K (SIPG, tets), both through solver/qd_linear_solver.hpp.
//
//   mpirun -np N ./seas_test_cudss_vs_mumps [--mesh FILE] [--mesh-scale S]
//          [--order P] [--nrhs 10] [--seed K] [--tol 1e-10] [--device cuda]
//          [--estimate-only]
//
// Default: each back end analyses and factors K once; then --nrhs stage
// right-hand sides b_k = B·δ_k + D·g(t_k) (random slip up to 10 m on every
// owned fault DOF, random t_k in [0, 1500 yr]) are solved by both, on
// device vectors:
//   1. ‖K u − b‖₂/‖b‖₂ < tol for every cuDSS solution (MUMPS's printed
//      alongside);
//   2. ‖u_cuDSS − u_MUMPS‖∞/‖u_MUMPS‖∞ < tol for every b_k.
// Printed: setup time of each, the device memory cuDSS took (drop in free
// memory over its setup), cuDSS's own analysis estimates (LU nonzeros,
// permanent and peak device and host memory) and the time per solve.
//
// --estimate-only: no factorization and no MUMPS; cuDSS's analysis phase
// alone on K and its estimates, for meshes whose factor may not fit
// (500 m, 250 m).
//
// A build without cuDSS (MFEM_USE_CUDSS=NO) runs the MUMPS half, says so,
// and exits 2.  Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../config/bp5_params.hpp"
#include "../../solver/qd_linear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

MPI_Comm comm = MPI_COMM_WORLD;
int myrank = 0;
int failures = 0;

void DeviceSync()
{
   if (Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK))
   {
      MFEM_DEVICE_SYNC;
   }
}

real_t GlobalNorm2(const Vector &v)
{
   return std::sqrt(InnerProduct(comm, v, v));
}

real_t GlobalMaxAbs(const Vector &v)
{
   real_t a = v.Normlinf();
   MPI_Allreduce(MPI_IN_PLACE, &a, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
   return a;
}

double MaxOverRanks(double x)
{
   MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MAX, comm);
   return x;
}

#ifdef MFEM_USE_CUDSS
#define CUDSS_TEST_CHECK(x)                                               \
   do {                                                                   \
      const cudssStatus_t st_ = (x);                                      \
      MFEM_VERIFY(st_ == CUDSS_STATUS_SUCCESS, "cuDSS call failed: " #x   \
                  " [code " << static_cast<int>(st_) << "]");             \
   } while (0)

/// cuDSS's analysis phase on K, set up as MFEM's CuDSSSolver does
/// (linalg/cudss.cpp: SPD, full view, each rank's rows of K merged from
/// hypre's diag and offd blocks, row distribution, the same communication
/// layer), and what it estimates for the factorization.  Collective.
struct CuDSSEstimate
{
   double seconds = 0.0;      ///< analysis (reordering + symbolic)
   int64_t lu_nnz = -1;       ///< nonzeros in the factor
   int64_t mem[16] = {};      ///< CUDSS_DATA_MEMORY_ESTIMATES
};

CuDSSEstimate EstimateCuDSS(const HypreParMatrix &K)
{
   CuDSSEstimate est;
   MPI_Comm kcomm = K.GetComm();
   hypre_ParCSRMatrix *parcsr = K;
   K.HypreRead();
   hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
   K.HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr);
#endif
   const int n_global = static_cast<int>(K.GetGlobalNumRows());
   const int row_start = static_cast<int>(parcsr->first_row_index);
   const int row_end = row_start + K.Height() - 1;
   const int nnz = static_cast<int>(hypre_CSRMatrixNumNonzeros(csr));

   cudssHandle_t handle;
   cudssConfig_t config;
   cudssData_t data;
   CUDSS_TEST_CHECK(cudssCreate(&handle));
   CUDSS_TEST_CHECK(cudssSetStream(handle, 0));
   const char *comm_lib = GetEnv("CUDSS_COMM_LIB");
#ifdef MFEM_CUDSS_COMM_LIB
   if (comm_lib == nullptr) { comm_lib = MFEM_CUDSS_COMM_LIB; }
#endif
   CUDSS_TEST_CHECK(cudssSetCommLayer(handle, comm_lib));
   CUDSS_TEST_CHECK(cudssConfigCreate(&config));
   CUDSS_TEST_CHECK(cudssDataCreate(handle, &data));
#if CUDSS_VERSION >= 800
   CUDSS_TEST_CHECK(cudssDataSet(handle, data, CUDSS_DATA_COMM_HOST, &kcomm,
                                 sizeof(MPI_Comm *)));
   const auto index_t = CUDSS_R_32I;
   const auto value_t = CUDSS_R_64F;
#else
   CUDSS_TEST_CHECK(cudssDataSet(handle, data, CUDSS_DATA_COMM, &kcomm,
                                 sizeof(MPI_Comm *)));
   const auto index_t = CUDA_R_32I;
   const auto value_t = CUDA_R_64F;
#endif
   static_assert(sizeof(real_t) == 8, "the estimate assumes double");

   cudssMatrix_t A, x, b;
#if CUDSS_VERSION >= 800
   CUDSS_TEST_CHECK(cudssMatrixCreateCsr(
                       &A, n_global, n_global, nnz, csr->i, NULL, csr->j,
                       csr->data, index_t, index_t, value_t, CUDSS_MTYPE_SPD,
                       CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO));
#else
   CUDSS_TEST_CHECK(cudssMatrixCreateCsr(
                       &A, n_global, n_global, nnz, csr->i, NULL, csr->j,
                       csr->data, index_t, value_t, CUDSS_MTYPE_SPD,
                       CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO));
#endif
   CUDSS_TEST_CHECK(cudssMatrixSetDistributionRow1d(A, row_start, row_end));
   CUDSS_TEST_CHECK(cudssMatrixCreateDn(&x, n_global, 1, n_global, NULL,
                                        value_t, CUDSS_LAYOUT_COL_MAJOR));
   CUDSS_TEST_CHECK(cudssMatrixCreateDn(&b, n_global, 1, n_global, NULL,
                                        value_t, CUDSS_LAYOUT_COL_MAJOR));
   CUDSS_TEST_CHECK(cudssMatrixSetDistributionRow1d(x, row_start, row_end));
   CUDSS_TEST_CHECK(cudssMatrixSetDistributionRow1d(b, row_start, row_end));

   MPI_Barrier(kcomm);
   const double t0 = MPI_Wtime();
   CUDSS_TEST_CHECK(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, config, data,
                                 A, x, b));
   DeviceSync();
   est.seconds = MaxOverRanks(MPI_Wtime() - t0);

   size_t written = 0;
   CUDSS_TEST_CHECK(cudssDataGet(handle, data, CUDSS_DATA_LU_NNZ, &est.lu_nnz,
                                 sizeof(est.lu_nnz), &written));
   CUDSS_TEST_CHECK(cudssDataGet(handle, data, CUDSS_DATA_MEMORY_ESTIMATES,
                                 est.mem, sizeof(est.mem), &written));

   CUDSS_TEST_CHECK(cudssMatrixDestroy(A));
   CUDSS_TEST_CHECK(cudssMatrixDestroy(x));
   CUDSS_TEST_CHECK(cudssMatrixDestroy(b));
   CUDSS_TEST_CHECK(cudssDataDestroy(handle, data));
   CUDSS_TEST_CHECK(cudssConfigDestroy(config));
   CUDSS_TEST_CHECK(cudssDestroy(handle));
   hypre_CSRMatrixDestroy(csr);
   return est;
}

/// Print the estimates of every rank (rank order) with the largest first.
void PrintEstimate(const CuDSSEstimate &est)
{
   int nranks;
   MPI_Comm_size(comm, &nranks);
   const int nv = 5;
   const double mine[nv] = {static_cast<double>(est.lu_nnz),
                            static_cast<double>(est.mem[0]),
                            static_cast<double>(est.mem[1]),
                            static_cast<double>(est.mem[2]),
                            static_cast<double>(est.mem[3])
                           };
   std::vector<double> all(nv * nranks);
   MPI_Gather(mine, nv, MPI_DOUBLE, all.data(), nv, MPI_DOUBLE, 0, comm);
   if (myrank != 0) { return; }
   const double gib = 1024.0 * 1024.0 * 1024.0;
   std::printf("  cuDSS analysis: %.3f s (slowest rank)\n", est.seconds);
   std::printf("  %-6s %16s %14s %14s %14s %14s\n", "rank", "LU nonzeros",
               "dev perm GiB", "dev peak GiB", "host perm GiB",
               "host peak GiB");
   for (int r = 0; r < nranks; r++)
   {
      const double *a = all.data() + nv * r;
      std::printf("  %-6d %16.0f %14.3f %14.3f %14.3f %14.3f\n", r, a[0],
                  a[1] / gib, a[2] / gib, a[3] / gib, a[4] / gib);
   }
   std::printf("  (CUDSS_DATA_LU_NNZ; CUDSS_DATA_MEMORY_ESTIMATES entries "
               "0-3: permanent and peak device, permanent and peak host "
               "memory in bytes)\n");
}
#endif // MFEM_USE_CUDSS

void Check(const std::string &name, real_t err, real_t bound)
{
   const bool ok = (err < bound);
   if (!ok) { failures++; }
   if (myrank == 0)
   {
      std::printf("  %-58s %10.3e  < %10.3e  %s\n", name.c_str(), err, bound,
                  ok ? "PASS" : "FAIL");
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   myrank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();

   std::string mesh_file = "bp5/mesh/reference/bp5_tandem_exact.msh";
   real_t mesh_scale = 1000.0;
   int order = 1;
   int nrhs = 10;
   int seed = 2026;
   real_t tol = 1e-10;
#ifdef MFEM_USE_CUDSS
   std::string device_str = "cuda";
#else
   std::string device_str = "cpu";
#endif
   bool estimate_only = false;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Gmsh BP5 mesh.");
   args.AddOption(&mesh_scale, "-s", "--mesh-scale", "Coordinate scale.");
   args.AddOption(&order, "-o", "--order", "Polynomial order.");
   args.AddOption(&nrhs, "-n", "--nrhs", "Right-hand sides per factor.");
   args.AddOption(&seed, "-r", "--seed", "Random seed.");
   args.AddOption(&tol, "-t", "--tol", "Residual and agreement tolerance.");
   args.AddOption(&device_str, "-d", "--device", "MFEM device string.");
   args.AddOption(&estimate_only, "-e", "--estimate-only", "-no-e",
                  "--no-estimate-only",
                  "cuDSS analysis estimates only (no factorization, no "
                  "MUMPS).");
   args.Parse();
   if (!args.Good())
   {
      if (myrank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }

   // One GPU per rank on a node, as seas_bp5_full chooses it.
   int device_id = 0;
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
   if (device_str.find("cuda") != std::string::npos ||
       device_str.find("hip") != std::string::npos)
   {
      const int num_gpu = Device::GetDeviceCount();
      if (num_gpu > 1)
      {
         MPI_Comm node_comm;
         MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                             &node_comm);
         int node_rank = 0;
         MPI_Comm_rank(node_comm, &node_rank);
         MPI_Comm_free(&node_comm);
         device_id = node_rank % num_gpu;
      }
   }
#endif
   Device device(device_str, device_id);
   if (myrank == 0) { device.Print(); }

   // ---- Mesh and K (as seas_bp5_full) ----
   Mesh smesh(mesh_file.c_str(), 1, 1);
   if (mesh_scale != 1.0)
   {
      for (int i = 0; i < smesh.GetNV(); i++)
      {
         real_t *v = smesh.GetVertex(i);
         for (int d = 0; d < 3; d++) { v[d] *= mesh_scale; }
      }
      smesh.SetAttributes();
   }
   ParMesh pmesh(comm, smesh);
   smesh.Clear();

   BP5Params params;
   LinearElastic model(params.lambda(), params.mu());
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(params.Vp);
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, model, params.Vp, params.Wf, params.lf, bc,
      DGMethod::IP, SolverType::MUMPS);
   double t0 = MPI_Wtime();
   domain.AssembleMatrixOnly();
   const double assembly_s = MaxOverRanks(MPI_Wtime() - t0);
   const HypreParMatrix &K = *domain.GetStiffnessMatrix();
   const long long global_ne = pmesh.GetGlobalNE();   // collective
   long long K_nnz = 0;
   {
      hypre_ParCSRMatrix *h = K;
      K_nnz = hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixDiag(h))
              + hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixOffd(h));
      MPI_Allreduce(MPI_IN_PLACE, &K_nnz, 1, MPI_LONG_LONG, MPI_SUM, comm);
   }
   if (myrank == 0)
   {
      std::cout << "G3a cuDSS vs MUMPS: " << mesh_file << ", p = " << order
                << ", " << nranks << " rank(s), " << global_ne
                << " elements, K " << K.GetGlobalNumRows() << " rows, "
                << K_nnz << " nonzeros (assembly " << assembly_s << " s)\n";
   }

#ifdef MFEM_USE_CUDSS
   {
      const CuDSSEstimate est = EstimateCuDSS(K);
      PrintEstimate(est);
   }
   if (estimate_only)
   {
      if (myrank == 0) { std::cout << "estimate only: done\n"; }
      return 0;
   }
#else
   if (estimate_only)
   {
      if (myrank == 0)
      {
         std::cout << "--estimate-only needs cuDSS (MFEM_USE_CUDSS=YES)\n";
      }
      return 2;
   }
#endif

   // ---- Stage right-hand sides b = B·δ + D·g(t) ----
   domain.BuildQDStageOperators(false);
   const int nown = domain.GetNumOwnedFaultDOFs();
   const int n = K.Height();
   const real_t yr = BP5Params::seconds_per_year;
   std::mt19937_64 rng_local(seed + 7919 * myrank);  // owned values
   std::mt19937_64 rng_global(seed);                  // same on every rank
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   std::vector<Vector> b(nrhs);
   Vector slip(2 * nown);
   for (int k = 0; k < nrhs; k++)
   {
      for (int i = 0; i < slip.Size(); i++) { slip(i) = 10.0 * U(rng_local); }
      const real_t t = 1500.0 * yr * 0.5 * (1.0 + U(rng_global));
      b[k].UseDevice(true);
      domain.QDStageRHS(t, slip, b[k]);
   }

   // ---- Back ends: one analysis and factorization each ----
   std::unique_ptr<Solver> mumps, cudss;
   QDLinearSolverSetup s_mumps, s_cudss;
   mumps = MakeQDLinearSolver(QDLinearSolverKind::MUMPS, K, 0, 0.0, &s_mumps);
   if (myrank == 0) { std::cout << QDLinearSolverSummary(s_mumps) << "\n"; }
   const bool have_cudss = QDLinearSolverAvailable(QDLinearSolverKind::CUDSS);
   if (have_cudss)
   {
      cudss = MakeQDLinearSolver(QDLinearSolverKind::CUDSS, K, 0, 0.0,
                                 &s_cudss);
      if (myrank == 0) { std::cout << QDLinearSolverSummary(s_cudss) << "\n"; }
   }
   else if (myrank == 0)
   {
      std::cout << "cuDSS: not in this build (MFEM_USE_CUDSS=NO); the MUMPS "
                << "half only\n";
   }

   // ---- Solves ----
   Vector um(n), uc(n), r(n);
   um.UseDevice(true);
   uc.UseDevice(true);
   r.UseDevice(true);
   std::vector<double> t_mumps(nrhs), t_cudss(nrhs, 0.0);
   real_t worst_res_c = 0.0, worst_res_m = 0.0, worst_diff = 0.0;
   if (myrank == 0)
   {
      std::printf("  %-4s %14s %14s %14s %12s %12s\n", "rhs", "res cuDSS",
                  "res MUMPS", "|uc-um|/|um|", "MUMPS ms", "cuDSS ms");
   }
   for (int k = 0; k < nrhs; k++)
   {
      const real_t bnorm = GlobalNorm2(b[k]);
      um = 0.0;
      DeviceSync();
      MPI_Barrier(comm);
      t0 = MPI_Wtime();
      mumps->Mult(b[k], um);
      DeviceSync();
      t_mumps[k] = MaxOverRanks(MPI_Wtime() - t0);
      K.Mult(um, r);
      r -= b[k];
      const real_t res_m = GlobalNorm2(r) / bnorm;
      worst_res_m = std::max(worst_res_m, res_m);

      real_t res_c = -1.0, diff = -1.0;
      if (have_cudss)
      {
         uc = 0.0;
         DeviceSync();
         MPI_Barrier(comm);
         t0 = MPI_Wtime();
         cudss->Mult(b[k], uc);
         DeviceSync();
         t_cudss[k] = MaxOverRanks(MPI_Wtime() - t0);
         K.Mult(uc, r);
         r -= b[k];
         res_c = GlobalNorm2(r) / bnorm;
         subtract(uc, um, r);
         diff = GlobalMaxAbs(r) / GlobalMaxAbs(um);
         worst_res_c = std::max(worst_res_c, res_c);
         worst_diff = std::max(worst_diff, diff);
      }
      if (myrank == 0)
      {
         std::printf("  %-4d %14.3e %14.3e %14.3e %12.3f %12.3f\n", k, res_c,
                     res_m, diff, 1e3 * t_mumps[k], 1e3 * t_cudss[k]);
      }
   }
   auto median = [](std::vector<double> v)
   {
      std::sort(v.begin(), v.end());
      return v.empty() ? 0.0 : v[v.size() / 2];
   };
   if (myrank == 0)
   {
      std::printf("  solve time per right-hand side, median (slowest rank): "
                  "MUMPS %.3f ms, cuDSS %.3f ms\n", 1e3 * median(t_mumps),
                  1e3 * median(t_cudss));
   }

   Check("MUMPS ||K u - b||_2/||b||_2, worst (reference)", worst_res_m, tol);
   if (!have_cudss)
   {
      if (myrank == 0) { std::cout << "cuDSS checks SKIPPED (no cuDSS)\n"; }
      return 2;
   }
   Check("cuDSS ||K u - b||_2/||b||_2, worst", worst_res_c, tol);
   Check("cuDSS vs MUMPS ||uc - um||_inf/||um||_inf, worst", worst_diff, tol);
   if (myrank == 0)
   {
      std::cout << (failures == 0 ? "ALL PASSED" : "FAILURES") << " ("
                << failures << " failed)\n";
   }
   return failures == 0 ? 0 : 1;
}
