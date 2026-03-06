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

#ifndef MFEM_SEAS_BP5_PARALLEL_OUTPUT_HPP
#define MFEM_SEAS_BP5_PARALLEL_OUTPUT_HPP

#include "mfem.hpp"
#include "bp5_benchmark_output.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../common/mpi_context.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Parallel wrapper for BP5BenchmarkOutput.
///
/// Gathers local fault data to root rank via FaultGeometry::GatherToRoot,
/// splits vector fields into per-component scalars, then root writes output
/// using the serial BP5BenchmarkOutput class. Only root rank opens files.
///
/// Follows ParallelBenchmarkOutput pattern exactly.
class ParallelBP5BenchmarkOutput
{
   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   /// @brief Construct parallel BP5 benchmark output.
   ///
   /// @param prefix Output file prefix
   /// @param params BP5 benchmark parameters
   /// @param stations Probe station list
   /// @param fault_geom Fault geometry with parallel gather support
   /// @param mpi_ctx MPI context
   /// @param global_x2 Globally gathered along-strike coords (root only)
   /// @param global_x3 Globally gathered depth coords (root only)
   /// @param global_tau_pre_dip Globally gathered tau_pre dip component (root)
   /// @param global_tau_pre_strike Globally gathered tau_pre strike component (root)
   ParallelBP5BenchmarkOutput(
      const std::string &prefix,
      const BP5Params &params,
      const std::vector<Probe2DInterpolator::Station> &stations,
      FaultGeometry<ParMesh> &fault_geom,
      MPIContext &mpi_ctx,
      const Vector &global_x2,
      const Vector &global_x3,
      const Vector &global_tau_pre_dip,
      const Vector &global_tau_pre_strike)
      : fault_geom_(fault_geom),
        mpi_ctx_(mpi_ctx),
        last_write_time_(-1e30)
   {
      if (mpi_ctx_.IsRoot())
      {
         bench_out_ = std::make_unique<BP5BenchmarkOutput<Mesh>>(
            prefix, params, stations, global_x2, global_x3);
         bench_out_->SetTauPre(global_tau_pre_dip, global_tau_pre_strike);
      }
   }

   /// @brief Write output if adaptive schedule requires it.
   ///
   /// Gathers local fault data to root, then root writes.
   /// All ranks must call this collectively.
   ///
   /// @param time Current simulation time [s]
   /// @param state Full state vector (local)
   /// @param fault Fault operator (BP5, SlipComponents=2)
   /// @param traction Local traction vector [2*N_local interleaved]
   /// @param global_V_max Global maximum slip rate (already reduced)
   /// @return true if data was written
   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<ParMesh, 2> &fault,
              const Vector &traction, real_t global_V_max)
   {
      // Synchronize write decision across all ranks
      real_t dt_out = BP5BenchmarkOutput<Mesh>::OutputInterval(global_V_max);
      int should_write =
         (time - last_write_time_ >= dt_out * kOutputTimeTolerance) ? 1 : 0;

#ifdef SEAS_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &should_write, 1, MPI_INT, MPI_MAX,
                    mpi_ctx_.GetComm());
#endif
      if (!should_write) { return false; }

      // Extract local fault quantities
      Vector local_slip, local_theta;
      fault.GetSlip(state, local_slip);       // [2*N_local interleaved]
      fault.GetTheta(state, local_theta);     // [N_local]
      const Vector &local_V = fault.GetSlipRate(); // [2*N_local interleaved]

      int N = fault.NumNodes();

      // Split interleaved vectors into per-component scalars
      Vector local_slip_dip(N), local_slip_strike(N);
      Vector local_V_dip(N), local_V_strike(N);
      Vector local_trac_dip(N), local_trac_strike(N);

      for (int i = 0; i < N; i++)
      {
         local_slip_dip(i)    = local_slip(2 * i + 0);
         local_slip_strike(i) = local_slip(2 * i + 1);
         local_V_dip(i)       = local_V(2 * i + 0);
         local_V_strike(i)    = local_V(2 * i + 1);
         local_trac_dip(i)    = traction(2 * i + 0);
         local_trac_strike(i) = traction(2 * i + 1);
      }

      // Gather 7 scalar fields to root
      Vector g_slip_dip, g_slip_strike, g_theta;
      Vector g_V_dip, g_V_strike;
      Vector g_trac_dip, g_trac_strike;

      fault_geom_.GatherToRoot(local_slip_dip, g_slip_dip);
      fault_geom_.GatherToRoot(local_slip_strike, g_slip_strike);
      fault_geom_.GatherToRoot(local_theta, g_theta);
      fault_geom_.GatherToRoot(local_V_dip, g_V_dip);
      fault_geom_.GatherToRoot(local_V_strike, g_V_strike);
      fault_geom_.GatherToRoot(local_trac_dip, g_trac_dip);
      fault_geom_.GatherToRoot(local_trac_strike, g_trac_strike);

      // Root writes using globally gathered data
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         bench_out_->WriteFromGlobalData(
            time, g_slip_dip, g_slip_strike, g_theta,
            g_V_dip, g_V_strike, g_trac_dip, g_trac_strike);
      }

      last_write_time_ = time;
      return true;
   }

   /// @brief Force a write at the current state (always flushes).
   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<ParMesh, 2> &fault,
                   const Vector &traction, real_t global_V_max)
   {
      last_write_time_ = -1e30;
      Write(time, state, fault, traction, global_V_max);
      Flush();
   }

   /// Flush all output files.
   void Flush()
   {
      if (mpi_ctx_.IsRoot() && bench_out_) { bench_out_->Flush(); }
   }

   /// Close all output files.
   void Close()
   {
      if (mpi_ctx_.IsRoot() && bench_out_) { bench_out_->Close(); }
   }

private:
   FaultGeometry<ParMesh> &fault_geom_;
   MPIContext &mpi_ctx_;
   std::unique_ptr<BP5BenchmarkOutput<Mesh>> bench_out_;
   real_t last_write_time_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_PARALLEL_OUTPUT_HPP
