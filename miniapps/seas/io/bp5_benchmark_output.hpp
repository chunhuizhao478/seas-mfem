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

#ifndef MFEM_SEAS_BP5_BENCHMARK_OUTPUT_HPP
#define MFEM_SEAS_BP5_BENCHMARK_OUTPUT_HPP

#include "mfem.hpp"
#include "probe_output.hpp"
#include "../config/bp5_params.hpp"
#include "../fault/rate_state_fault.hpp"

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>

namespace mfem
{
namespace seas
{

/// @brief 2D nearest-DOF probe matching for BP5 fault stations.
///
/// Unlike ProbeInterpolator (1D depth interpolation), this uses
/// simple nearest-neighbor matching in 2D (x2, x3) space.
/// No interpolation needed for DG nodes.
class Probe2DInterpolator
{
public:
   struct Station { std::string name; real_t x2; real_t x3; };

   /// @brief Construct from fault DOF coordinates and station list.
   ///
   /// @param fault_x2 Along-strike coordinates of fault DOFs [num_dofs]
   /// @param fault_x3 Depth coordinates of fault DOFs [num_dofs]
   /// @param stations List of probe stations with (name, x2, x3)
   Probe2DInterpolator(const Vector &fault_x2, const Vector &fault_x3,
                        const std::vector<Station> &stations)
      : num_stations_(static_cast<int>(stations.size())),
        num_dofs_(fault_x2.Size())
   {
      MFEM_VERIFY(fault_x2.Size() == fault_x3.Size(),
                  "Coordinate arrays must have the same size");

      nearest_dof_.resize(num_stations_);
      match_distance_.resize(num_stations_);

      for (int s = 0; s < num_stations_; s++)
      {
         FindNearest(fault_x2, fault_x3,
                     stations[s].x2, stations[s].x3, s);
      }
   }

   /// Get the DOF index closest to the given station.
   int GetNearestDOF(int station_idx) const
   {
      return nearest_dof_[station_idx];
   }

   /// Get the Euclidean distance from station to its matched DOF.
   real_t GetMatchDistance(int station_idx) const
   {
      return match_distance_[station_idx];
   }

   /// Number of stations.
   int NumStations() const { return num_stations_; }

private:
   int num_stations_;
   int num_dofs_;
   std::vector<int> nearest_dof_;
   std::vector<real_t> match_distance_;

   void FindNearest(const Vector &x2, const Vector &x3,
                    real_t target_x2, real_t target_x3, int s)
   {
      if (num_dofs_ == 0)
      {
         nearest_dof_[s] = -1;
         match_distance_[s] = std::numeric_limits<real_t>::max();
         return;
      }

      int best = 0;
      real_t best_dist = std::numeric_limits<real_t>::max();

      for (int i = 0; i < num_dofs_; i++)
      {
         real_t dx = x2(i) - target_x2;
         real_t dz = x3(i) - target_x3;
         real_t dist = std::sqrt(dx * dx + dz * dz);
         if (dist < best_dist)
         {
            best_dist = dist;
            best = i;
         }
      }

      nearest_dof_[s] = best;
      match_distance_[s] = best_dist;
   }
};

/// @brief BP5 SCEC benchmark output with 8-column vector format.
///
/// Writes time series at 2D fault stations in SCEC BP5 format:
///   t(s), slip_strike(m), slip_dip(m), log10(V_strike)(m/s),
///   log10(V_dip)(m/s), tau_strike(MPa), tau_dip(MPa), log10(state)(s)
///
/// Component ordering: Internal MFEM [dip=0, strike=1] -> Output [strike, dip]
///
/// Follows BenchmarkOutput pattern: same adaptive output, same
/// ProbeOutput usage, same ForceWrite/Flush/Close interface.
template <typename MeshType = Mesh>
class BP5BenchmarkOutput
{
   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   using Station = Probe2DInterpolator::Station;

   /// @brief Construct benchmark output manager for BP5.
   ///
   /// @param prefix Output file prefix (e.g., "output/bp5_full")
   /// @param params BP5 benchmark parameters
   /// @param stations Probe station list
   /// @param fault_x2 Along-strike coords of fault DOFs [num_dofs]
   /// @param fault_x3 Depth coords of fault DOFs [num_dofs]
   BP5BenchmarkOutput(const std::string &prefix,
                      const BP5Params &params,
                      const std::vector<Station> &stations,
                      const Vector &fault_x2, const Vector &fault_x3)
      : prefix_(prefix),
        stations_(stations),
        interpolator_(fault_x2, fault_x3, stations),
        last_write_time_(-1e30)
   {
      // Store per-DOF tau_pre components for WriteFromGlobalData
      int n = fault_x2.Size();
      tau_pre_dip_.SetSize(n);
      tau_pre_strike_.SetSize(n);
      tau_pre_dip_ = 0.0;
      tau_pre_strike_ = 0.0;

      std::vector<std::string> columns = {
         "time(s)", "slip_strike(m)", "slip_dip(m)",
         "log10(V_strike)(m/s)", "log10(V_dip)(m/s)",
         "tau_strike(MPa)", "tau_dip(MPa)", "log10(state)(s)"
      };

      for (size_t i = 0; i < stations.size(); i++)
      {
         real_t x2_km = stations[i].x2 / 1000.0;
         real_t x3_km = stations[i].x3 / 1000.0;

         std::ostringstream desc;
         desc << "BP5-QD time series at " << stations[i].name
              << " (x2=" << x2_km << "km, x3=" << x3_km << "km)";

         std::string filename = prefix + "_" + stations[i].name + ".txt";
         auto probe = std::make_unique<ProbeOutput>(
            filename, columns, desc.str());
         probes_.push_back(std::move(probe));
      }
   }

   /// @brief Set pre-stress components for WriteFromGlobalData path.
   ///
   /// Called once after construction with globally gathered tau_pre.
   /// @param tau_pre_dip Per-DOF dip pre-stress [num_dofs]
   /// @param tau_pre_strike Per-DOF strike pre-stress [num_dofs]
   void SetTauPre(const Vector &tau_pre_dip, const Vector &tau_pre_strike)
   {
      tau_pre_dip_ = tau_pre_dip;
      tau_pre_strike_ = tau_pre_strike;
   }

   /// @brief Write output if adaptive schedule requires it.
   ///
   /// @param time Current simulation time [s]
   /// @param state Full state vector
   /// @param fault Fault operator (BP5, SlipComponents=2)
   /// @param traction Traction from domain solve [2*N interleaved]
   /// @param V_max Maximum slip rate [m/s]
   /// @return true if data was written
   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<MeshType, 2> &fault,
              const Vector &traction, real_t V_max)
   {
      real_t dt_out = OutputInterval(V_max);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }

      Vector slip, theta;
      fault.GetSlip(state, slip);       // [2*N interleaved: dip, strike]
      fault.GetTheta(state, theta);     // [N]
      const Vector &slip_rate = fault.GetSlipRate(); // [2*N interleaved]

      const Vector &tau_pre = fault.GetGeometry()->GetTauPre(); // [2*N]

      WriteRow(time, slip, theta, slip_rate, traction, tau_pre);

      last_write_time_ = time;
      return true;
   }

   /// @brief Force a write at the current state (always flushes).
   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<MeshType, 2> &fault,
                   const Vector &traction, real_t V_max)
   {
      last_write_time_ = -1e30;
      Write(time, state, fault, traction, V_max);
      Flush();
   }

   /// @brief Write from pre-gathered global component data (parallel path).
   ///
   /// All fields are per-DOF scalars of size M (globally gathered):
   /// Uses stored tau_pre_dip_ and tau_pre_strike_ for stress computation.
   void WriteFromGlobalData(real_t time,
                            const Vector &global_slip_dip,
                            const Vector &global_slip_strike,
                            const Vector &global_theta,
                            const Vector &global_V_dip,
                            const Vector &global_V_strike,
                            const Vector &global_trac_dip,
                            const Vector &global_trac_strike)
   {
      for (int s = 0; s < interpolator_.NumStations(); s++)
      {
         int dof = interpolator_.GetNearestDOF(s);
         if (dof < 0) { continue; }

         // Negate slip for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t slip_dip    = -global_slip_dip(dof);
         real_t slip_strike = -global_slip_strike(dof);

         real_t V_dip    = std::abs(global_V_dip(dof));
         real_t V_strike = std::abs(global_V_strike(dof));

         // Negate for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t tau_dip    = -(tau_pre_dip_(dof) +
                               global_trac_dip(dof)) / 1e6;
         real_t tau_strike = -(tau_pre_strike_(dof) +
                               global_trac_strike(dof)) / 1e6;

         real_t th = global_theta(dof);

         // Output: strike first, dip second (SCEC convention)
         std::vector<real_t> row = {
            time,
            slip_strike,
            slip_dip,
            V_strike > 0.0 ? std::log10(V_strike) : -300.0,
            V_dip > 0.0 ? std::log10(V_dip) : -300.0,
            tau_strike,
            tau_dip,
            th > 0.0 ? std::log10(th) : -300.0
         };
         probes_[s]->WriteStep(row);
      }

      last_write_time_ = time;
   }

   /// Flush all output files.
   void Flush()
   {
      for (auto &p : probes_) { p->Flush(); }
   }

   /// Close all output files.
   void Close()
   {
      for (auto &p : probes_) { p->Close(); }
   }

   /// Number of probes.
   int NumProbes() const { return static_cast<int>(probes_.size()); }

   /// @brief Compute adaptive output interval based on slip rate.
   /// Same thresholds as BenchmarkOutput.
   static real_t OutputInterval(real_t V_max)
   {
      if (V_max > 1e-3)
      {
         return 0.1;    // Coseismic: every 0.1 s (SCEC spec)
      }
      else if (V_max > 1e-6)
      {
         return 0.1;    // Nucleation: every 0.1 s
      }
      else
      {
         return 0.1 * BP5Params::seconds_per_year;  // Interseismic: ~0.1 yr (SCEC spec)
      }
   }

   /// @brief Default 10 SCEC BP5 on-fault stations (Section 4.1 of spec).
   static std::vector<Station> DefaultStations()
   {
      return {
         {"fltst_strk-36dp+00", -36e3,  0.0},
         {"fltst_strk-16dp+00", -16e3,  0.0},
         {"fltst_strk+00dp+00",   0.0,  0.0},
         {"fltst_strk+16dp+00",  16e3,  0.0},
         {"fltst_strk+36dp+00",  36e3,  0.0},
         {"fltst_strk-24dp+10", -24e3, 10e3},
         {"fltst_strk-16dp+10", -16e3, 10e3},
         {"fltst_strk+00dp+10",   0.0, 10e3},
         {"fltst_strk+16dp+10",  16e3, 10e3},
         {"fltst_strk+00dp+22",   0.0, 22e3},
      };
   }

private:
   std::string prefix_;
   std::vector<Station> stations_;
   Probe2DInterpolator interpolator_;
   std::vector<std::unique_ptr<ProbeOutput>> probes_;
   real_t last_write_time_;

   // Pre-stress components for WriteFromGlobalData path
   Vector tau_pre_dip_;
   Vector tau_pre_strike_;

   /// Write one row per station from interleaved data (serial path).
   void WriteRow(real_t time,
                 const Vector &slip,       // [2*N interleaved]
                 const Vector &theta,      // [N]
                 const Vector &slip_rate,  // [2*N interleaved]
                 const Vector &traction,   // [2*N interleaved]
                 const Vector &tau_pre)    // [2*N interleaved]
   {
      for (int s = 0; s < interpolator_.NumStations(); s++)
      {
         int dof = interpolator_.GetNearestDOF(s);
         if (dof < 0) { continue; }

         // Internal: index 0 = dip, index 1 = strike
         // Negate slip for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t slip_dip    = -slip(2 * dof + 0);
         real_t slip_strike = -slip(2 * dof + 1);

         real_t V_dip    = std::abs(slip_rate(2 * dof + 0));
         real_t V_strike = std::abs(slip_rate(2 * dof + 1));

         // Total shear stress = tau_pre + elastic_traction
         // Negate for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t tau_dip = -(tau_pre(2 * dof + 0) +
                            traction(2 * dof + 0)) / 1e6;
         real_t tau_strike = -(tau_pre(2 * dof + 1) +
                               traction(2 * dof + 1)) / 1e6;

         real_t th = theta(dof);

         // Output: strike first, dip second (SCEC convention)
         std::vector<real_t> row = {
            time,
            slip_strike,
            slip_dip,
            V_strike > 0.0 ? std::log10(V_strike) : -300.0,
            V_dip > 0.0 ? std::log10(V_dip) : -300.0,
            tau_strike,
            tau_dip,
            th > 0.0 ? std::log10(th) : -300.0
         };
         probes_[s]->WriteStep(row);
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_BENCHMARK_OUTPUT_HPP
