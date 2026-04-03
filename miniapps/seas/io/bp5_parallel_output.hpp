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
#include "probe_output.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../common/mpi_context.hpp"

#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

/// @brief Distributed parallel BP5 benchmark output (Tandem-style).
///
/// Each rank locates probe stations on its local fault DOFs and faces.
/// Ownership is assigned via MPI_Allreduce: the rank with minimum geometric
/// distance owns each station (ties broken by lower rank number). Each rank
/// writes only its owned stations from local data — no global gather needed
/// during time stepping.
///
/// This matches Tandem's BoundaryProbeWriter architecture:
///   - BoundaryPointLocator assigns each probe to a rank
///   - clean_duplicate_probes() uses MPI_Allreduce(MIN) for ownership
///   - Each rank writes its probes independently
class ParallelBP5BenchmarkOutput
{
   static constexpr real_t kOutputTimeTolerance = 0.99;

   using Station = Probe2DInterpolator::Station;

public:
   ParallelBP5BenchmarkOutput(
      const std::string &prefix,
      const BP5Params &params,
      const std::vector<Station> &stations,
      FaultGeometry<ParMesh> &fault_geom,
      MPIContext &mpi_ctx,
      const Vector &local_x2,
      const Vector &local_x3,
      const Vector &local_tau_pre_dip,
      const Vector &local_tau_pre_strike,
      int nbf_per_face = 1,
      int face_basis_type = BasisType::GaussLobatto)
      : fault_geom_(fault_geom),
        mpi_ctx_(mpi_ctx),
        prefix_(prefix),
        stations_(stations),
        eta_(params.eta()),
        last_write_time_(-1e30)
   {
      int num_stations = static_cast<int>(stations.size());

      // Step 1: Each rank creates a local interpolator using LOCAL coords.
      // This finds the nearest DOF and (if nbf >= 3) tries exact face match.
      Probe2DInterpolator local_interp(
         local_x2, local_x3, stations, nbf_per_face, face_basis_type, false);

      // Step 2: Tandem-style distributed ownership via MPI_Allreduce.
      // Each rank reports its distance to each station; the global minimum
      // determines ownership. Ties broken by lower rank number.
      owned_stations_.clear();
      std::vector<real_t> local_dist(num_stations);
      for (int s = 0; s < num_stations; s++)
      {
         local_dist[s] = local_interp.GetMatchDistance(s);
         // Ranks without any fault DOFs get infinite distance
         if (local_interp.GetNearestDOF(s) < 0 &&
             !local_interp.HasExactMatch(s))
         {
            local_dist[s] = std::numeric_limits<real_t>::max();
         }
      }

      std::vector<real_t> global_min_dist(num_stations);
#ifdef SEAS_USE_MPI
      MPI_Allreduce(local_dist.data(), global_min_dist.data(),
                    num_stations, MPI_DOUBLE, MPI_MIN, mpi_ctx.GetComm());
#else
      global_min_dist = local_dist;
#endif

      // Tie-break by rank number (matching Tandem's clean_duplicate_probes)
      int rank = 0;
#ifdef SEAS_USE_MPI
      MPI_Comm_rank(mpi_ctx.GetComm(), &rank);
#endif
      std::vector<int> local_owner_rank(num_stations);
      for (int s = 0; s < num_stations; s++)
      {
         local_owner_rank[s] =
            (local_dist[s] == global_min_dist[s])
               ? rank : std::numeric_limits<int>::max();
      }
      std::vector<int> global_owner_rank(num_stations);
#ifdef SEAS_USE_MPI
      MPI_Allreduce(local_owner_rank.data(), global_owner_rank.data(),
                    num_stations, MPI_INT, MPI_MIN, mpi_ctx.GetComm());
#else
      global_owner_rank = local_owner_rank;
#endif

      // Step 3: Store ownership and create probe files for owned stations.
      local_tau_pre_dip_ = local_tau_pre_dip;
      local_tau_pre_strike_ = local_tau_pre_strike;
      local_interp_ = std::make_unique<Probe2DInterpolator>(
         local_x2, local_x3, stations, nbf_per_face, face_basis_type, false);

      std::vector<std::string> columns = {
         "time(s)", "slip_strike(m)", "slip_dip(m)",
         "log10(V_strike)(m/s)", "log10(V_dip)(m/s)",
         "tau_strike(MPa)", "tau_dip(MPa)", "log10(state)(s)"
      };

      for (int s = 0; s < num_stations; s++)
      {
         if (global_owner_rank[s] == rank)
         {
            owned_stations_.push_back(s);

            real_t x2_km = stations[s].x2 / 1000.0;
            real_t x3_km = stations[s].x3 / 1000.0;
            std::ostringstream desc;
            desc << "BP5-QD time series at " << stations[s].name
                 << " (x2=" << x2_km << "km, x3=" << x3_km << "km)";
            std::string filename = prefix + "_" + stations[s].name + ".txt";
            auto probe = std::make_unique<ProbeOutput>(
               filename, columns, desc.str());
            probes_.push_back(std::move(probe));
         }
      }

      // Print ownership summary on root
      if (mpi_ctx.IsRoot())
      {
         std::cout << "  Distributed probe ownership:\n";
         for (int s = 0; s < num_stations; s++)
         {
            std::cout << "    " << stations[s].name
                      << " → rank " << global_owner_rank[s]
                      << " (dist=" << global_min_dist[s] << " m";
            if (global_owner_rank[s] == rank &&
                local_interp.HasExactMatch(s))
            {
               std::cout << ", exact";
            }
            std::cout << ")\n";
         }
      }
   }

   void EnableTractionDecompositionOutput()
   {
      if (owned_stations_.empty()) { return; }
      if (!decomp_probes_.empty()) { return; }

      std::vector<std::string> columns = {
         "time(s)",
         "tau_stress_strike(MPa)", "tau_stress_dip(MPa)",
         "tau_corr_strike(MPa)", "tau_corr_dip(MPa)",
         "tau_total_strike(MPa)", "tau_total_dip(MPa)"
      };

      for (int idx = 0; idx < static_cast<int>(owned_stations_.size()); idx++)
      {
         int s = owned_stations_[idx];
         std::string filename =
            prefix_ + "_tracdec_" + stations_[s].name + ".txt";
         std::ostringstream desc;
         desc << "BP5 traction decomposition at " << stations_[s].name;
         decomp_probes_.push_back(
            std::make_unique<ProbeOutput>(filename, columns, desc.str()));
      }
   }

   void EnableJumpResidualOutput()
   {
      if (owned_stations_.empty()) { return; }
      if (!jump_res_probes_.empty()) { return; }

      std::vector<std::string> columns = {
         "time(s)", "jump_res_strike", "jump_res_dip", "jump_res_mag"
      };

      for (int idx = 0; idx < static_cast<int>(owned_stations_.size()); idx++)
      {
         int s = owned_stations_[idx];
         std::string filename =
            prefix_ + "_jumpres_" + stations_[s].name + ".txt";
         std::ostringstream desc;
         desc << "BP5 jump residual at " << stations_[s].name;
         jump_res_probes_.push_back(
            std::make_unique<ProbeOutput>(filename, columns, desc.str()));
      }
   }

   void PrintDiagnostics(const Vector &local_x2, const Vector &local_x3,
                         std::ostream &os = std::cout) const
   {
      if (!owned_stations_.empty() && local_interp_)
      {
         local_interp_->PrintDiagnostics(stations_, local_x2, local_x3, os);
      }
   }

   /// @brief Write output if adaptive schedule requires it.
   ///
   /// Each rank evaluates owned stations from local data — no gather.
   /// All ranks must call collectively for the write-decision sync.
   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<ParMesh, 2> &fault,
              const Vector &traction, real_t global_V_max)
   {
      real_t dt_out = BP5BenchmarkOutput<Mesh>::OutputInterval(global_V_max);
      int should_write =
         (time - last_write_time_ >= dt_out * kOutputTimeTolerance) ? 1 : 0;
#ifdef SEAS_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &should_write, 1, MPI_INT, MPI_MAX,
                    mpi_ctx_.GetComm());
#endif
      if (!should_write) { return false; }

      if (!owned_stations_.empty())
      {
         // Extract local data
         Vector local_slip, local_theta;
         fault.GetSlip(state, local_slip);
         fault.GetTheta(state, local_theta);
         const Vector &local_V = fault.GetSlipRate();
         int N = fault.NumNodes();

         Vector slip_dip(N), slip_strike(N), V_dip(N), V_strike(N);
         Vector trac_dip(N), trac_strike(N);
         for (int i = 0; i < N; i++)
         {
            slip_dip(i)    = local_slip(2*i);
            slip_strike(i) = local_slip(2*i+1);
            V_dip(i)       = local_V(2*i);
            V_strike(i)    = local_V(2*i+1);
            trac_dip(i)    = traction(2*i);
            trac_strike(i) = traction(2*i+1);
         }

         // Write owned stations from local data
         for (int idx = 0; idx < static_cast<int>(owned_stations_.size()); idx++)
         {
            int s = owned_stations_[idx];

            real_t sd = -local_interp_->EvaluateScalar(slip_dip, s);
            real_t ss = -local_interp_->EvaluateScalar(slip_strike, s);
            real_t vd = std::abs(local_interp_->EvaluateScalar(V_dip, s));
            real_t vs = std::abs(local_interp_->EvaluateScalar(V_strike, s));

            real_t tau_d = -(local_interp_->EvaluateScalar(local_tau_pre_dip_, s)
                           + local_interp_->EvaluateScalar(trac_dip, s)
                           + eta_ * local_interp_->EvaluateScalar(V_dip, s)) / 1e6;
            real_t tau_s = -(local_interp_->EvaluateScalar(local_tau_pre_strike_, s)
                           + local_interp_->EvaluateScalar(trac_strike, s)
                           + eta_ * local_interp_->EvaluateScalar(V_strike, s)) / 1e6;

            real_t th = local_interp_->EvaluateScalar(local_theta, s);

            // Distinguish non-finite (NaN, inf) from genuine V=0.
            // Without this, both map to -300 via the "v > 0.0" check.
            auto safe_log10 = [](real_t v) -> real_t {
               if (!std::isfinite(v)) { return -999.0; }
               return v > 0.0 ? std::log10(v) : -300.0;
            };
            std::vector<real_t> row = {
               time, ss, sd,
               safe_log10(vs), safe_log10(vd),
               tau_s, tau_d,
               safe_log10(th)
            };
            probes_[idx]->WriteStep(row);
         }
      }

      last_write_time_ = time;
      return true;
   }

   void WriteTractionDecomposition(real_t time,
                                   const Vector &traction_stress,
                                   const Vector &traction_correction)
   {
      if (decomp_probes_.empty() || owned_stations_.empty()) { return; }

      const int N = fault_geom_.NumLocalFaultDOFs();
      Vector sd(N), ss(N), cd(N), cs(N);
      for (int i = 0; i < N; i++)
      {
         sd(i) = traction_stress(2*i);
         ss(i) = traction_stress(2*i+1);
         cd(i) = traction_correction(2*i);
         cs(i) = traction_correction(2*i+1);
      }

      for (int idx = 0; idx < static_cast<int>(owned_stations_.size()); idx++)
      {
         int s = owned_stations_[idx];
         real_t tsd = -local_interp_->EvaluateScalar(sd, s) / 1e6;
         real_t tss = -local_interp_->EvaluateScalar(ss, s) / 1e6;
         real_t tcd = -local_interp_->EvaluateScalar(cd, s) / 1e6;
         real_t tcs = -local_interp_->EvaluateScalar(cs, s) / 1e6;

         std::vector<real_t> row = {
            time, tss, tsd, tcs, tcd, tss + tcs, tsd + tcd
         };
         decomp_probes_[idx]->WriteStep(row);
      }
   }

   void WriteJumpResidual(real_t time, const Vector &jump_residual)
   {
      if (jump_res_probes_.empty() || owned_stations_.empty()) { return; }

      const int N = fault_geom_.NumLocalFaultDOFs();
      Vector rd(N), rs(N);
      for (int i = 0; i < N; i++)
      {
         rd(i) = jump_residual(2*i);
         rs(i) = jump_residual(2*i+1);
      }

      for (int idx = 0; idx < static_cast<int>(owned_stations_.size()); idx++)
      {
         int s = owned_stations_[idx];
         real_t res_d = local_interp_->EvaluateScalar(rd, s);
         real_t res_s = local_interp_->EvaluateScalar(rs, s);
         real_t res_mag = std::sqrt(res_d*res_d + res_s*res_s);

         std::vector<real_t> row = { time, res_s, res_d, res_mag };
         jump_res_probes_[idx]->WriteStep(row);
      }
   }

   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<ParMesh, 2> &fault,
                   const Vector &traction, real_t global_V_max)
   {
      last_write_time_ = -1e30;
      Write(time, state, fault, traction, global_V_max);
      Flush();
   }

   void Flush()
   {
      for (auto &p : probes_) { p->Flush(); }
      for (auto &p : decomp_probes_) { p->Flush(); }
      for (auto &p : jump_res_probes_) { p->Flush(); }
   }

   void Close()
   {
      for (auto &p : probes_) { p->Close(); }
      for (auto &p : decomp_probes_) { p->Close(); }
      for (auto &p : jump_res_probes_) { p->Close(); }
   }

private:
   FaultGeometry<ParMesh> &fault_geom_;
   MPIContext &mpi_ctx_;
   std::string prefix_;
   std::vector<Station> stations_;
   real_t eta_ = 0.0;

   // Distributed ownership: indices into stations_ that this rank owns
   std::vector<int> owned_stations_;

   // Local interpolator and tau_pre (used only for owned stations)
   std::unique_ptr<Probe2DInterpolator> local_interp_;
   Vector local_tau_pre_dip_;
   Vector local_tau_pre_strike_;

   // Per-owned-station output files
   std::vector<std::unique_ptr<ProbeOutput>> probes_;
   std::vector<std::unique_ptr<ProbeOutput>> decomp_probes_;
   std::vector<std::unique_ptr<ProbeOutput>> jump_res_probes_;

   real_t last_write_time_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_PARALLEL_OUTPUT_HPP
