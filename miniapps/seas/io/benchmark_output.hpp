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

#ifndef MFEM_SEAS_BENCHMARK_OUTPUT_HPP
#define MFEM_SEAS_BENCHMARK_OUTPUT_HPP

#include "mfem.hpp"
#include "probe_output.hpp"
#include "../config/bp2_params.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../solver/seas_operator.hpp"

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace mfem
{
namespace seas
{

/// @brief Handles SCEC benchmark format output for SEAS simulations.
///
/// Writes time series data at specified probe depths in the SCEC-standard
/// format:  time(s), slip(m), log10(slip_rate)(m/s), shear_stress(MPa),
///          log10(state)(s)
///
/// Output files follow the SCEC naming convention:
///   {prefix}_{benchmark}_z{depth}km.txt
///
/// Adaptive output frequency based on maximum slip rate:
///   V < 1e-6 m/s:  interseismic (every 0.01 years)
///   1e-6 < V < 1e-3: nucleation (every 0.1 seconds)
///   V > 1e-3 m/s:  coseismic (every 0.001 seconds)
template <typename MeshType = Mesh>
class BenchmarkOutput
{
   /// Tolerance factor to avoid writing twice at the same time due to
   /// floating-point rounding in adaptive time stepping.
   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   /// @brief Construct benchmark output manager.
   ///
   /// @param prefix Output file prefix (e.g., "output/mfem_bp2qd")
   /// @param params BP2 benchmark parameters
   /// @param probe_depths Probe depths in meters (z=0 at surface, z<0 at depth)
   /// @param fault_depths Fault DOF depths from the domain operator
   BenchmarkOutput(const std::string &prefix,
                   const BP2Params &params,
                   const std::vector<real_t> &probe_depths,
                   const Vector &fault_depths)
      : prefix_(prefix),
        probe_depths_(probe_depths),
        interpolator_(fault_depths, probe_depths),
        last_write_time_(-1e30)
   {
      // Create one output file per probe
      std::vector<std::string> columns = {
         "time(s)", "slip(m)", "log10(slip_rate)(m/s)",
         "shear_stress(MPa)", "log10(state)(s)"
      };

      for (size_t i = 0; i < probe_depths.size(); i++)
      {
         std::string filename = MakeFilename(prefix, probe_depths[i]);
         real_t depth_km = std::abs(probe_depths[i]) / 1000.0;

         std::ostringstream desc;
         desc << "BP2-QD time series at z = " << depth_km << " km";

         auto probe = std::make_unique<ProbeOutput>(
            filename, columns, desc.str());
         probes_.push_back(std::move(probe));
      }
   }

   /// @brief Write output if the adaptive schedule requires it.
   ///
   /// Checks whether enough time has elapsed since the last write
   /// based on the current phase (interseismic/nucleation/coseismic).
   ///
   /// @param time Current simulation time [s]
   /// @param state Full state vector [slip_0, theta_0, ...]
   /// @param fault Fault operator (provides slip rate, traction data)
   /// @param traction Traction vector from domain solve
   /// @return true if data was written
   bool Write(real_t time, const Vector &state,
              const RateStateFaultOperator<MeshType> &fault,
              const Vector &traction)
   {
      real_t V_max = fault.GetMaxSlipRate();
      real_t dt_out = OutputInterval(V_max);

      if (time - last_write_time_ < dt_out * kOutputTimeTolerance) { return false; }

      // Extract fault quantities
      Vector slip, theta;
      fault.GetSlip(state, slip);
      fault.GetTheta(state, theta);

      const Vector &slip_rate = fault.GetSlipRate();
      real_t tau0 = fault.GetTau0();

      // Interpolate to probe locations and write
      Vector probe_slip(interpolator_.NumProbes());
      Vector probe_theta(interpolator_.NumProbes());
      Vector probe_V(interpolator_.NumProbes());
      Vector probe_tau(interpolator_.NumProbes());

      interpolator_.Interpolate(slip, probe_slip);
      interpolator_.Interpolate(theta, probe_theta);
      interpolator_.Interpolate(slip_rate, probe_V);
      interpolator_.Interpolate(traction, probe_tau);

      for (int p = 0; p < interpolator_.NumProbes(); p++)
      {
         MFEM_ASSERT(probe_V(p) >= 0, "Negative slip rate at probe " << p << ": " << probe_V(p));
         MFEM_ASSERT(probe_theta(p) >= 0, "Negative state variable at probe " << p << ": " << probe_theta(p));
         real_t V = std::max(probe_V(p), 1e-30);
         real_t th = std::max(probe_theta(p), 1e-30);
         real_t tau_total = tau0 + probe_tau(p);

         std::vector<real_t> row = {
            time,
            probe_slip(p),
            std::log10(V),
            tau_total / 1e6,  // Convert Pa to MPa
            std::log10(th)
         };
         probes_[p]->WriteStep(row);
      }

      last_write_time_ = time;
      return true;
   }

   /// @brief Force a write at the current state (ignoring schedule, always flushes).
   void ForceWrite(real_t time, const Vector &state,
                   const RateStateFaultOperator<MeshType> &fault,
                   const Vector &traction)
   {
      last_write_time_ = -1e30;  // Reset to force write
      Write(time, state, fault, traction);
      Flush();
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

   /// @brief Compute adaptive output interval based on slip rate.
   ///
   /// @param V_max Maximum slip rate [m/s]
   /// @return Output interval [s]
   static real_t OutputInterval(real_t V_max)
   {
      if (V_max > 1e-3)
      {
         return 0.001;  // Coseismic: every 0.001 s
      }
      else if (V_max > 1e-6)
      {
         return 0.1;    // Nucleation: every 0.1 s
      }
      else
      {
         return 0.01 * BP2Params::seconds_per_year;  // Interseismic: every 0.01 yr
      }
   }

   /// @brief Write from pre-gathered global data (for parallel output).
   ///
   /// Allows parallel code to gather data and pass it directly.
   void WriteFromGlobalData(real_t time,
                            const Vector &global_slip,
                            const Vector &global_theta,
                            const Vector &global_V,
                            const Vector &global_traction,
                            real_t tau0)
   {
      Vector probe_slip(interpolator_.NumProbes());
      Vector probe_theta(interpolator_.NumProbes());
      Vector probe_V(interpolator_.NumProbes());
      Vector probe_tau(interpolator_.NumProbes());

      interpolator_.Interpolate(global_slip, probe_slip);
      interpolator_.Interpolate(global_theta, probe_theta);
      interpolator_.Interpolate(global_V, probe_V);
      interpolator_.Interpolate(global_traction, probe_tau);

      for (int p = 0; p < interpolator_.NumProbes(); p++)
      {
         real_t V = std::max(probe_V(p), 1e-30);
         real_t th = std::max(probe_theta(p), 1e-30);
         real_t tau_total = tau0 + probe_tau(p);

         std::vector<real_t> row = {
            time,
            probe_slip(p),
            std::log10(V),
            tau_total / 1e6,
            std::log10(th)
         };
         probes_[p]->WriteStep(row);
      }

      last_write_time_ = time;
   }

   /// Number of probes.
   int NumProbes() const { return static_cast<int>(probes_.size()); }

private:
   std::string prefix_;
   std::vector<real_t> probe_depths_;
   ProbeInterpolator interpolator_;
   std::vector<std::unique_ptr<ProbeOutput>> probes_;
   real_t last_write_time_;

   /// @brief Generate SCEC-convention filename.
   ///
   /// E.g. "output/mfem_bp2qd_z0km.txt", "output/mfem_bp2qd_z4.8km.txt"
   static std::string MakeFilename(const std::string &prefix, real_t depth_m)
   {
      real_t depth_km = std::abs(depth_m) / 1000.0;

      std::ostringstream oss;
      oss << prefix << "_z";

      // Format depth: integer if exact, else one decimal place
      if (std::abs(depth_km - std::round(depth_km)) < 1e-6)
      {
         oss << static_cast<int>(std::round(depth_km));
      }
      else
      {
         oss << std::fixed << std::setprecision(1) << depth_km;
      }
      oss << "km.txt";
      return oss.str();
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BENCHMARK_OUTPUT_HPP
