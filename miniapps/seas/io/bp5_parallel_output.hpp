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
   /// @param nbf_per_face Number of fault DOFs per face in gathered ordering
   /// @param face_basis_type BasisType used for fault face nodes
   ParallelBP5BenchmarkOutput(
      const std::string &prefix,
      const BP5Params &params,
      const std::vector<Probe2DInterpolator::Station> &stations,
      FaultGeometry<ParMesh> &fault_geom,
      MPIContext &mpi_ctx,
      const Vector &global_x2,
      const Vector &global_x3,
      const Vector &global_tau_pre_dip,
      const Vector &global_tau_pre_strike,
      int nbf_per_face = 1,
      int face_basis_type = BasisType::GaussLobatto)
      : fault_geom_(fault_geom),
        mpi_ctx_(mpi_ctx),
        nbf_per_face_(nbf_per_face),
        last_write_time_(-1e30)
   {
      if (mpi_ctx_.IsRoot())
      {
         // Build face deduplication mapping to handle shared-face duplicates
         // at partition boundaries. Without this, the contiguous-per-face
         // assumption in TryBuildExactMatch would break.
         if (nbf_per_face > 1)
         {
            FaultGeometry<ParMesh>::BuildFaceDedupMap(
               global_x2, global_x3, nbf_per_face,
               dedup_face_indices_);

            int num_raw = global_x2.Size();
            int num_dedup = static_cast<int>(dedup_face_indices_.size())
                            * nbf_per_face;
            if (num_dedup < num_raw)
            {
               std::cout << "  Face dedup: " << num_raw / nbf_per_face
                         << " raw faces → "
                         << dedup_face_indices_.size()
                         << " unique faces ("
                         << (num_raw - num_dedup) / nbf_per_face
                         << " duplicates removed)\n";
            }

            // Deduplicate coordinates and tau_pre
            Vector dedup_x2, dedup_x3, dedup_tau_dip, dedup_tau_strike;
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_x2, dedup_face_indices_, nbf_per_face, dedup_x2);
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_x3, dedup_face_indices_, nbf_per_face, dedup_x3);
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_tau_pre_dip, dedup_face_indices_, nbf_per_face,
               dedup_tau_dip);
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_tau_pre_strike, dedup_face_indices_, nbf_per_face,
               dedup_tau_strike);

            bench_out_ = std::make_unique<BP5BenchmarkOutput<Mesh>>(
               prefix, params, stations, dedup_x2, dedup_x3,
               nbf_per_face, face_basis_type);
            bench_out_->SetTauPre(dedup_tau_dip, dedup_tau_strike);
         }
         else
         {
            // nbf=1: no face structure to deduplicate, use raw data directly
            bench_out_ = std::make_unique<BP5BenchmarkOutput<Mesh>>(
               prefix, params, stations, global_x2, global_x3,
               nbf_per_face, face_basis_type);
            bench_out_->SetTauPre(global_tau_pre_dip, global_tau_pre_strike);
         }
      }
   }

   void EnableTractionDecompositionOutput()
   {
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         bench_out_->EnableTractionDecompositionOutput();
      }
   }

   void EnableJumpResidualOutput()
   {
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         bench_out_->EnableJumpResidualOutput();
      }
   }

   /// Print station mapping diagnostics on root.
   void PrintDiagnostics(const Vector &global_x2, const Vector &global_x3,
                         std::ostream &os = std::cout) const
   {
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         if (!dedup_face_indices_.empty())
         {
            Vector d_x2, d_x3;
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_x2, dedup_face_indices_, nbf_per_face_, d_x2);
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               global_x3, dedup_face_indices_, nbf_per_face_, d_x3);
            bench_out_->PrintDiagnostics(d_x2, d_x3, os);
         }
         else
         {
            bench_out_->PrintDiagnostics(global_x2, global_x3, os);
         }
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

      // Root writes using globally gathered data (deduplicated if needed)
      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         if (!dedup_face_indices_.empty())
         {
            Vector d_sd, d_ss, d_th, d_vd, d_vs, d_td, d_ts;
            auto dedup = [&](const Vector &raw, Vector &out) {
               FaultGeometry<ParMesh>::ApplyFaceDedupMap(
                  raw, dedup_face_indices_, nbf_per_face_, out);
            };
            dedup(g_slip_dip, d_sd);
            dedup(g_slip_strike, d_ss);
            dedup(g_theta, d_th);
            dedup(g_V_dip, d_vd);
            dedup(g_V_strike, d_vs);
            dedup(g_trac_dip, d_td);
            dedup(g_trac_strike, d_ts);
            bench_out_->WriteFromGlobalData(
               time, d_sd, d_ss, d_th, d_vd, d_vs, d_td, d_ts);
         }
         else
         {
            bench_out_->WriteFromGlobalData(
               time, g_slip_dip, g_slip_strike, g_theta,
               g_V_dip, g_V_strike, g_trac_dip, g_trac_strike);
         }
      }

      last_write_time_ = time;
      return true;
   }

   /// Write station-level traction decomposition at the current time.
   ///
   /// Inputs are local interleaved traction components [2*N_local]:
   ///   traction_stress = physical stress contribution
   ///   traction_correction = correction contribution added to the stress part
   void WriteTractionDecomposition(real_t time,
                                   const Vector &traction_stress,
                                   const Vector &traction_correction)
   {
      if (!bench_out_ && !mpi_ctx_.IsRoot())
      {
         // Non-root still participates in gather below.
      }

      const int N = fault_geom_.NumLocalFaultDOFs();
      MFEM_ASSERT(traction_stress.Size() == 2 * N,
                  "traction_stress size mismatch");
      MFEM_ASSERT(traction_correction.Size() == 2 * N,
                  "traction_correction size mismatch");

      Vector local_stress_dip(N), local_stress_strike(N);
      Vector local_corr_dip(N), local_corr_strike(N);
      for (int i = 0; i < N; i++)
      {
         local_stress_dip(i) = traction_stress(2 * i + 0);
         local_stress_strike(i) = traction_stress(2 * i + 1);
         local_corr_dip(i) = traction_correction(2 * i + 0);
         local_corr_strike(i) = traction_correction(2 * i + 1);
      }

      Vector g_stress_dip, g_stress_strike, g_corr_dip, g_corr_strike;
      fault_geom_.GatherToRoot(local_stress_dip, g_stress_dip);
      fault_geom_.GatherToRoot(local_stress_strike, g_stress_strike);
      fault_geom_.GatherToRoot(local_corr_dip, g_corr_dip);
      fault_geom_.GatherToRoot(local_corr_strike, g_corr_strike);

      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         if (!dedup_face_indices_.empty())
         {
            Vector d_sd, d_ss, d_cd, d_cs;
            auto dedup = [&](const Vector &raw, Vector &out) {
               FaultGeometry<ParMesh>::ApplyFaceDedupMap(
                  raw, dedup_face_indices_, nbf_per_face_, out);
            };
            dedup(g_stress_dip, d_sd);
            dedup(g_stress_strike, d_ss);
            dedup(g_corr_dip, d_cd);
            dedup(g_corr_strike, d_cs);
            bench_out_->WriteTractionDecompositionFromGlobalData(
               time, d_sd, d_ss, d_cd, d_cs);
         }
         else
         {
            bench_out_->WriteTractionDecompositionFromGlobalData(
               time, g_stress_dip, g_stress_strike, g_corr_dip, g_corr_strike);
         }
      }
   }

   /// Write station-level jump residual at the current time.
   ///
   /// Input is a local interleaved residual field [2*N_local]:
   ///   jump_residual = [[u]] - delta projected into local [dip, strike]
   void WriteJumpResidual(real_t time, const Vector &jump_residual)
   {
      if (!bench_out_ && !mpi_ctx_.IsRoot())
      {
         // Non-root still participates in gather below.
      }

      const int N = fault_geom_.NumLocalFaultDOFs();
      MFEM_ASSERT(jump_residual.Size() == 2 * N,
                  "jump_residual size mismatch");

      Vector local_res_dip(N), local_res_strike(N);
      for (int i = 0; i < N; i++)
      {
         local_res_dip(i) = jump_residual(2 * i + 0);
         local_res_strike(i) = jump_residual(2 * i + 1);
      }

      Vector g_res_dip, g_res_strike;
      fault_geom_.GatherToRoot(local_res_dip, g_res_dip);
      fault_geom_.GatherToRoot(local_res_strike, g_res_strike);

      if (mpi_ctx_.IsRoot() && bench_out_)
      {
         if (!dedup_face_indices_.empty())
         {
            Vector d_rd, d_rs;
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               g_res_dip, dedup_face_indices_, nbf_per_face_, d_rd);
            FaultGeometry<ParMesh>::ApplyFaceDedupMap(
               g_res_strike, dedup_face_indices_, nbf_per_face_, d_rs);
            bench_out_->WriteJumpResidualFromGlobalData(
               time, d_rd, d_rs);
         }
         else
         {
            bench_out_->WriteJumpResidualFromGlobalData(
               time, g_res_dip, g_res_strike);
         }
      }
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
   int nbf_per_face_ = 1;
   std::vector<int> dedup_face_indices_;  // Face dedup mapping (root only)
   std::unique_ptr<BP5BenchmarkOutput<Mesh>> bench_out_;
   real_t last_write_time_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_PARALLEL_OUTPUT_HPP
