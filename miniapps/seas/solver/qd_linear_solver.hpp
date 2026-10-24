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

// qd_linear_solver.hpp — the direct solve K u = b of a quasi-dynamic stage
// behind one interface (GPU plan G3).  K does not change during a run, so a
// back end analyses and factors it once, in MakeQDLinearSolver, and every
// stage is one Mult.  Back ends:
//
//   MUMPS, MUMPS_BLR   CPU (MUMPS through MFEM), the production solver.
//   CUDSS              GPU (MFEM's CuDSSSolver, MFEM_USE_CUDSS).  Under MPI
//                      it gives cuDSS each rank's rows of K (row
//                      distribution, cuDSS's multi-GPU mode, also on one
//                      rank) and loads the communication layer from
//                      CUDSS_COMM_LIB or the path MFEM was built with.
//                      Mult reads and writes device vectors.
//
// Both are given K as symmetric positive definite, the full matrix.

#ifndef MFEM_SEAS_QD_LINEAR_SOLVER_HPP
#define MFEM_SEAS_QD_LINEAR_SOLVER_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

enum class QDLinearSolverKind { MUMPS, MUMPS_BLR, CUDSS };

/// What setting up a back end cost, for the log.
struct QDLinearSolverSetup
{
   std::string name;
   double seconds = 0.0;      ///< analysis + factorization, slowest rank
   double device_bytes = -1;  ///< device memory taken, largest rank (-1: CPU)
};

/// True when this MFEM build has the back end @a kind.
inline bool QDLinearSolverAvailable(QDLinearSolverKind kind)
{
   switch (kind)
   {
#ifdef MFEM_USE_MUMPS
      case QDLinearSolverKind::MUMPS:
      case QDLinearSolverKind::MUMPS_BLR: return true;
#endif
#ifdef MFEM_USE_CUDSS
      case QDLinearSolverKind::CUDSS: return true;
#endif
      default: return false;
   }
}

namespace detail
{
/// Free device memory in bytes (0 without a GPU backend).
inline double FreeDeviceBytes()
{
#if defined(MFEM_USE_CUDA)
   if (Device::Allows(Backend::CUDA_MASK))
   {
      size_t free_b = 0, total_b = 0;
      MFEM_GPU_CHECK(cudaMemGetInfo(&free_b, &total_b));
      return static_cast<double>(free_b);
   }
#endif
   return 0.0;
}
} // namespace detail

/// Build the back end @a kind for @a K: analysis and factorization happen
/// here.  @a mumps_print_level and @a blr_tol apply to the MUMPS back ends.
/// Collective over the communicator of @a K.
inline std::unique_ptr<Solver> MakeQDLinearSolver(
   QDLinearSolverKind kind, const HypreParMatrix &K, int mumps_print_level,
   double blr_tol, QDLinearSolverSetup *setup = nullptr)
{
   MFEM_VERIFY(QDLinearSolverAvailable(kind),
               "QD linear solver: this MFEM build has no "
               << (kind == QDLinearSolverKind::CUDSS ? "cuDSS (MFEM_USE_CUDSS)"
                                                     : "MUMPS (MFEM_USE_MUMPS)"));
   MPI_Comm comm = K.GetComm();
   MPI_Barrier(comm);
   const double free0 = detail::FreeDeviceBytes();
   const double t0 = MPI_Wtime();
   std::unique_ptr<Solver> solver;
   std::string name;
   switch (kind)
   {
#ifdef MFEM_USE_MUMPS
      case QDLinearSolverKind::MUMPS:
      case QDLinearSolverKind::MUMPS_BLR:
      {
         auto *mumps = new MUMPSSolver(comm);
         mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
         mumps->SetPrintLevel(mumps_print_level);
         if (kind == QDLinearSolverKind::MUMPS_BLR) { mumps->SetBLRTol(blr_tol); }
         mumps->SetOperator(K);
         solver.reset(mumps);
         name = (kind == QDLinearSolverKind::MUMPS) ? "MUMPS" : "MUMPS BLR";
         break;
      }
#endif
#ifdef MFEM_USE_CUDSS
      case QDLinearSolverKind::CUDSS:
      {
         auto *cudss = new CuDSSSolver(comm);
         cudss->SetMatrixSymType(CuDSSSolver::SYMMETRIC_POSITIVE_DEFINITE);
         cudss->SetMatrixViewType(CuDSSSolver::FULL);
         cudss->SetOperator(K);
         solver.reset(cudss);
         int nranks;
         MPI_Comm_size(comm, &nranks);
         name = "cuDSS (" + std::to_string(nranks) + " rank"
                + (nranks > 1 ? "s" : "") + ", row distribution)";
         break;
      }
#endif
      default: break;
   }
   if (Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK))
   {
      MFEM_DEVICE_SYNC;
   }
   const double seconds = MPI_Wtime() - t0;
   const double used = (kind == QDLinearSolverKind::CUDSS)
                       ? free0 - detail::FreeDeviceBytes() : -1.0;
   if (setup)
   {
      setup->name = name;
      MPI_Allreduce(&seconds, &setup->seconds, 1, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(&used, &setup->device_bytes, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   return solver;
}

/// One line for the log, e.g. "[QD-SOLVER] cuDSS (1 rank, row distribution):
/// analysis + factorization 1.234 s (slowest rank), device memory 2.10 GiB
/// (largest rank, ...)".
inline std::string QDLinearSolverSummary(const QDLinearSolverSetup &s)
{
   std::ostringstream o;
   o << "[QD-SOLVER] " << s.name << ": analysis + factorization "
     << std::fixed << std::setprecision(3) << s.seconds << " s (slowest rank)";
   if (s.device_bytes >= 0.0)
   {
      o << ", device memory " << std::setprecision(2)
        << s.device_bytes / (1024.0 * 1024.0 * 1024.0)
        << " GiB (largest rank, drop in free memory over the setup)";
   }
   return o.str();
}

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_LINEAR_SOLVER_HPP
