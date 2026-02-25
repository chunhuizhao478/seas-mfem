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

#ifndef MFEM_SEAS_PARALLEL_BENCHMARK_OUTPUT_HPP
#define MFEM_SEAS_PARALLEL_BENCHMARK_OUTPUT_HPP

#include "mfem.hpp"
#include "benchmark_output.hpp"
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

/// @brief Parallel wrapper for BenchmarkOutput.
///
/// Gathers local fault data to root rank via FaultGeometry::GatherToRootDedup,
/// deduplicates partition-boundary duplicate DOFs, then root writes output
/// using the serial BenchmarkOutput class. Only root rank opens files.
class ParallelBenchmarkOutput
{
   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   /// @brief Construct parallel benchmark output.
   ///
   /// @param prefix Output file prefix
   /// @param params BP2 benchmark parameters
   /// @param probe_depths Probe depths in meters
   /// @param fault_geom Fault geometry with parallel gather support
   /// @param mpi_ctx MPI context
   /// @param dedup_fault_depths Deduplicated global fault depths (root only)
   ParallelBenchmarkOutput(const std::string &prefix,
                            const BP2Params &params,
                            const std::vector<real_t> &probe_depths,
                            FaultGeometry<ParMesh> &fault_geom,
                            MPIContext &mpi_ctx,
                            const Vector &dedup_fault_depths)
      : fault_geom_(fault_geom),
        mpi_ctx_(mpi_ctx),
        last_write_time_(-1e30)
   {
      if (mpi_ctx_.IsRoot())
      {
         bench_out_ = std::make_unique<BenchmarkOutput<Mesh>>(
            prefix, params, probe_depths, dedup_fault_depths);
      }
   }

   /// @brief Write output if the adaptive schedule requires it.
   ///
   /// Gathers local fault data, deduplicates on root, then writes.
   ///
   /// @param time Current simulation time [s]
   /// @param state Full state vector (local)
   /// @param fault Fault operator (provides slip rate, traction data)
   /// @param traction Local traction vector from domain solve
   /// @param global_V_max Global maximum slip rate (already reduced)
   /// @return true if data was written
   /// @pre global_V_max must be globally reduced (same value on all ranks).
   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<ParMesh> &fault,
              const Vector &traction, real_t global_V_max)
   {
      // Compute whether this rank thinks we should write
      real_t dt_out = BenchmarkOutput<Mesh>::OutputInterval(global_V_max);
      int should_write = (time - last_write_time_ >= dt_out * kOutputTimeTolerance) ? 1 : 0;

      // Synchronize decision across all ranks to prevent MPI_Gatherv deadlock.
      // If ANY rank decides to write, ALL ranks must participate in the gather.
#ifdef SEAS_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &should_write, 1, MPI_INT, MPI_MAX,
                    mpi_ctx_.GetComm());
#endif
      if (!should_write)
      {
         return false;
      }

      // Extract local fault quantities
      Vector local_slip, local_theta;
      fault.GetSlip(state, local_slip);
      fault.GetTheta(state, local_theta);
      const Vector &local_V = fault.GetSlipRate();

      // Gather and deduplicate all fields on root
      std::vector<const Vector*> local_fields = {
         &local_slip, &local_theta, &local_V, &traction
      };
      std::vector<Vector> dedup_fields;
      Vector dedup_depths;
      fault_geom_.GatherFieldsToRootDedup(local_fields, dedup_fields,
                                           dedup_depths);

      // Root writes using deduplicated data
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         bench_out_->WriteFromGlobalData(time,
                                          dedup_fields[0],  // slip
                                          dedup_fields[1],  // theta
                                          dedup_fields[2],  // V
                                          dedup_fields[3],  // traction
                                          fault.GetTau0(),
                                          fault.GetParams().eta());
      }

      last_write_time_ = time;
      return true;
   }

   /// @brief Force a write at the current state (always flushes to disk).
   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<ParMesh> &fault,
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
   std::unique_ptr<BenchmarkOutput<Mesh>> bench_out_;
   real_t last_write_time_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PARALLEL_BENCHMARK_OUTPUT_HPP
