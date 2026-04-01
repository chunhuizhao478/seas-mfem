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
#include <utility>

namespace mfem
{
namespace seas
{

/// @brief 2D fault-station matching/interpolation for BP5 output.
///
/// The legacy path snaps each station to the nearest fault DOF. That is
/// exact only when the requested station coincides with a fault node.
/// Tandem instead locates the containing fault face and evaluates the
/// finite-element field at the requested point. This class now does the
/// same when the gathered fault data preserves per-face ordering
/// (`nbf_per_face >= 3` and contiguous face blocks). If exact face
/// interpolation cannot be constructed, it falls back to nearest-DOF.
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
                       const std::vector<Station> &stations,
                       int nbf_per_face = 1,
                       int face_basis_type = BasisType::GaussLobatto,
                       bool emit_warnings = true)
      : num_stations_(static_cast<int>(stations.size())),
        num_dofs_(fault_x2.Size()),
        nbf_per_face_(nbf_per_face),
        face_basis_type_(face_basis_type)
   {
      MFEM_VERIFY(fault_x2.Size() == fault_x3.Size(),
                  "Coordinate arrays must have the same size");

      nearest_dof_.resize(num_stations_);
      match_distance_.resize(num_stations_);
      exact_match_.resize(num_stations_, false);
      face_start_.resize(num_stations_, -1);
      exact_weights_.resize(num_stations_);

      for (int s = 0; s < num_stations_; s++)
      {
         FindNearest(fault_x2, fault_x3,
                     stations[s].x2, stations[s].x3, s);
         TryBuildExactMatch(fault_x2, fault_x3,
                            stations[s].x2, stations[s].x3, s);
      }

      // Warn about stations using nearest-DOF snap instead of exact
      // face interpolation — the output method is different.
      if (emit_warnings)
      {
         for (int s = 0; s < num_stations_; s++)
         {
            if (!exact_match_[s] && nbf_per_face_ >= 3)
            {
               std::cout << "  WARNING: Station " << stations[s].name
                         << " at (" << stations[s].x2 << ", "
                         << stations[s].x3 << ") using nearest-DOF snap "
                         << "(dist=" << match_distance_[s]
                         << " m) instead of exact face interpolation.\n";
            }
         }
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

   /// Whether the station uses exact face-local interpolation.
   bool HasExactMatch(int station_idx) const
   {
      return exact_match_[station_idx];
   }

   /// Get the face start index for a station (valid only if HasExactMatch).
   int GetFaceStart(int station_idx) const { return face_start_[station_idx]; }

   /// Get the exact interpolation weights for a station.
   const Vector &GetExactWeights(int station_idx) const
   {
      return exact_weights_[station_idx];
   }

   /// Get nbf_per_face used by this interpolator.
   int GetNbfPerFace() const { return nbf_per_face_; }

   /// Print diagnostic information for all stations.
   ///
   /// For each station, prints whether exact interpolation is active,
   /// the face_start index, weights, nearest DOF, match distance,
   /// and the coordinates of the contributing DOFs.
   void PrintDiagnostics(const std::vector<Station> &stations,
                         const Vector &x2, const Vector &x3,
                         std::ostream &os = std::cout) const
   {
      os << "\n=== Probe2DInterpolator Diagnostics ===\n";
      os << "  num_dofs = " << num_dofs_
         << ", nbf_per_face = " << nbf_per_face_
         << ", num_faces = "
         << (nbf_per_face_ > 0 ? num_dofs_ / nbf_per_face_ : 0)
         << "\n";
      if (nbf_per_face_ >= 3 && num_dofs_ % nbf_per_face_ != 0)
      {
         os << "  WARNING: num_dofs % nbf_per_face != 0 ("
            << num_dofs_ << " % " << nbf_per_face_
            << " = " << (num_dofs_ % nbf_per_face_)
            << ") — exact interpolation disabled\n";
      }

      for (int s = 0; s < num_stations_; s++)
      {
         os << "\n  Station " << s;
         if (s < static_cast<int>(stations.size()))
         {
            os << " [" << stations[s].name
               << "] target=(" << stations[s].x2
               << ", " << stations[s].x3 << ")";
         }
         os << "\n";

         os << "    nearest_dof = " << nearest_dof_[s]
            << ", match_distance = " << match_distance_[s] << " m\n";

         if (nearest_dof_[s] >= 0 && nearest_dof_[s] < num_dofs_)
         {
            os << "    nearest coords = ("
               << x2(nearest_dof_[s]) << ", "
               << x3(nearest_dof_[s]) << ")\n";
         }

         os << "    exact_match = "
            << (exact_match_[s] ? "true" : "false") << "\n";

         if (exact_match_[s])
         {
            int start = face_start_[s];
            os << "    face_start = " << start << "\n";
            os << "    weights = [";
            for (int k = 0; k < exact_weights_[s].Size(); k++)
            {
               if (k > 0) { os << ", "; }
               os << exact_weights_[s](k);
            }
            os << "]\n";
            os << "    face DOF coords:\n";
            int nbf = std::min(nbf_per_face_, num_dofs_ - start);
            for (int k = 0; k < nbf; k++)
            {
               os << "      DOF " << (start + k) << ": ("
                  << x2(start + k) << ", " << x3(start + k) << ")\n";
            }

            // Validate: check triangle non-degeneracy
            if (nbf >= 3)
            {
               real_t ax = x2(start + 1) - x2(start);
               real_t az = x3(start + 1) - x3(start);
               real_t bx = x2(start + 2) - x2(start);
               real_t bz = x3(start + 2) - x3(start);
               real_t det = ax * bz - az * bx;
               os << "    triangle det = " << det;
               if (std::abs(det) < 1e-14)
               {
                  os << " *** DEGENERATE ***";
               }
               os << "\n";
            }
         }
      }
      os << "=== End Diagnostics ===\n\n";
   }

   /// Detect duplicate faces in the global coordinate array.
   ///
   /// Returns the number of face pairs that share the same vertex
   /// coordinates (within tolerance). This indicates shared-face
   /// duplication from the parallel gather.
   static int CountDuplicateFaces(const Vector &x2, const Vector &x3,
                                  int nbf_per_face, real_t tol = 1.0,
                                  std::ostream *os = nullptr)
   {
      if (nbf_per_face < 3) { return 0; }
      int ndofs = x2.Size();
      if (ndofs % nbf_per_face != 0) { return 0; }
      int nfaces = ndofs / nbf_per_face;
      int duplicates = 0;

      for (int i = 0; i < nfaces; i++)
      {
         int si = i * nbf_per_face;
         for (int j = i + 1; j < nfaces; j++)
         {
            int sj = j * nbf_per_face;
            // Compare first vertex (sufficient for face identification)
            real_t dx = x2(si) - x2(sj);
            real_t dz = x3(si) - x3(sj);
            if (std::sqrt(dx * dx + dz * dz) < tol)
            {
               // Confirm with second vertex
               real_t dx1 = x2(si + 1) - x2(sj + 1);
               real_t dz1 = x3(si + 1) - x3(sj + 1);
               if (std::sqrt(dx1 * dx1 + dz1 * dz1) < tol)
               {
                  duplicates++;
                  if (os)
                  {
                     *os << "  Duplicate faces: " << i << " and " << j
                         << " at (" << x2(si) << ", " << x3(si) << ")\n";
                  }
               }
            }
         }
      }
      return duplicates;
   }

   /// Evaluate a scalar field [num_dofs] at the given station.
   real_t EvaluateScalar(const Vector &field, int station_idx) const
   {
      MFEM_ASSERT(field.Size() == num_dofs_,
                  "EvaluateScalar size mismatch: got " << field.Size()
                  << ", expected " << num_dofs_);

      if (exact_match_[station_idx])
      {
         const int start = face_start_[station_idx];
         const Vector &w = exact_weights_[station_idx];
         real_t value = 0.0;
         for (int k = 0; k < w.Size(); k++)
         {
            value += w(k) * field(start + k);
         }
         return value;
      }

      const int dof = nearest_dof_[station_idx];
      return (dof >= 0) ? field(dof) : 0.0;
   }

   /// Evaluate one component of an interleaved field [ncomp*num_dofs].
   real_t EvaluateInterleaved(const Vector &field, int station_idx,
                              int comp, int ncomp = 2) const
   {
      MFEM_ASSERT(field.Size() == ncomp * num_dofs_,
                  "EvaluateInterleaved size mismatch: got " << field.Size()
                  << ", expected " << ncomp * num_dofs_);
      MFEM_ASSERT(comp >= 0 && comp < ncomp,
                  "Component index out of range");

      if (exact_match_[station_idx])
      {
         const int start = face_start_[station_idx];
         const Vector &w = exact_weights_[station_idx];
         real_t value = 0.0;
         for (int k = 0; k < w.Size(); k++)
         {
            value += w(k) * field(ncomp * (start + k) + comp);
         }
         return value;
      }

      const int dof = nearest_dof_[station_idx];
      return (dof >= 0) ? field(ncomp * dof + comp) : 0.0;
   }

private:
   int num_stations_;
   int num_dofs_;
   int nbf_per_face_;
   int face_basis_type_;
   std::vector<int> nearest_dof_;
   std::vector<real_t> match_distance_;
   std::vector<bool> exact_match_;
   std::vector<int> face_start_;
   std::vector<Vector> exact_weights_;

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

   static int OrderFromNbf(int nbf)
   {
      if (nbf == 1) { return 0; }

      const real_t disc = std::sqrt(1.0 + 8.0 * nbf);
      const int p = static_cast<int>(std::llround((disc - 3.0) / 2.0));
      return ((p + 1) * (p + 2) / 2 == nbf) ? p : -1;
   }

   static bool ComputeReferenceIP(real_t x0, real_t z0,
                                  real_t x1, real_t z1,
                                  real_t x2, real_t z2,
                                  real_t xp, real_t zp,
                                  IntegrationPoint &ip)
   {
      const real_t ax = x1 - x0;
      const real_t az = z1 - z0;
      const real_t bx = x2 - x0;
      const real_t bz = z2 - z0;
      const real_t px = xp - x0;
      const real_t pz = zp - z0;
      const real_t det = ax * bz - az * bx;
      if (std::abs(det) < 1e-14) { return false; }

      const real_t r = (px * bz - pz * bx) / det;
      const real_t s = (ax * pz - az * px) / det;
      const real_t l0 = 1.0 - r - s;
      const real_t l1 = r;
      const real_t l2 = s;
      const real_t tol = 1e-10;
      if (l0 < -tol || l1 < -tol || l2 < -tol) { return false; }

      ip.x = r;
      ip.y = s;
      ip.weight = 0.0;
      return true;
   }

   void TryBuildExactMatch(const Vector &x2, const Vector &x3,
                           real_t target_x2, real_t target_x3, int s)
   {
      if (nbf_per_face_ < 3 || num_dofs_ % nbf_per_face_ != 0)
      {
         return;
      }

      const int order = OrderFromNbf(nbf_per_face_);
      if (order < 1) { return; }

      H1_TriangleElement face_fe(order, face_basis_type_);
      Vector shape(nbf_per_face_);
      const int num_faces = num_dofs_ / nbf_per_face_;

      for (int f = 0; f < num_faces; f++)
      {
         const int start = f * nbf_per_face_;
         IntegrationPoint ip;
         if (!ComputeReferenceIP(x2(start + 0), x3(start + 0),
                                 x2(start + 1), x3(start + 1),
                                 x2(start + 2), x3(start + 2),
                                 target_x2, target_x3, ip))
         {
            continue;
         }

         face_fe.CalcShape(ip, shape);
         exact_match_[s] = true;
         face_start_[s] = start;
         exact_weights_[s] = shape;
         match_distance_[s] = 0.0;
         return;
      }
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
                      const Vector &fault_x2, const Vector &fault_x3,
                      int nbf_per_face = 1,
                      int face_basis_type = BasisType::GaussLobatto)
      : prefix_(prefix),
        stations_(stations),
        interpolator_(fault_x2, fault_x3, stations,
                      nbf_per_face, face_basis_type),
        last_write_time_(-1e30),
        eta_(params.eta())
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

   /// Enable writing separate station files for traction decomposition.
   ///
   /// Files are written as:
   ///   <prefix>_tracdec_<station>.txt
   /// with columns:
   ///   t, tau_stress_strike, tau_stress_dip, tau_corr_strike,
   ///   tau_corr_dip, tau_total_strike, tau_total_dip
   void EnableTractionDecompositionOutput()
   {
      if (!decomp_probes_.empty()) { return; }

      std::vector<std::string> columns = {
         "time(s)",
         "tau_stress_strike(MPa)", "tau_stress_dip(MPa)",
         "tau_corr_strike(MPa)", "tau_corr_dip(MPa)",
         "tau_total_strike(MPa)", "tau_total_dip(MPa)"
      };

      for (size_t i = 0; i < stations_.size(); i++)
      {
         real_t x2_km = stations_[i].x2 / 1000.0;
         real_t x3_km = stations_[i].x3 / 1000.0;

         std::ostringstream desc;
         desc << "BP5 traction decomposition at " << stations_[i].name
              << " (x2=" << x2_km << "km, x3=" << x3_km << "km)";

         std::string filename =
            prefix_ + "_tracdec_" + stations_[i].name + ".txt";
         auto probe = std::make_unique<ProbeOutput>(
            filename, columns, desc.str());
         decomp_probes_.push_back(std::move(probe));
      }
   }

   /// Enable writing station files for fault jump residual diagnostics.
   ///
   /// Files are written as:
   ///   <prefix>_jumpres_<station>.txt
   /// with columns:
   ///   t, jump_res_strike, jump_res_dip, jump_res_mag
   void EnableJumpResidualOutput()
   {
      if (!jump_res_probes_.empty()) { return; }

      std::vector<std::string> columns = {
         "time(s)",
         "jump_res_strike(m)", "jump_res_dip(m)",
         "jump_res_mag(m)"
      };

      for (size_t i = 0; i < stations_.size(); i++)
      {
         real_t x2_km = stations_[i].x2 / 1000.0;
         real_t x3_km = stations_[i].x3 / 1000.0;

         std::ostringstream desc;
         desc << "BP5 jump residual at " << stations_[i].name
              << " (x2=" << x2_km << "km, x3=" << x3_km << "km)";

         std::string filename =
            prefix_ + "_jumpres_" + stations_[i].name + ".txt";
         auto probe = std::make_unique<ProbeOutput>(
            filename, columns, desc.str());
         jump_res_probes_.push_back(std::move(probe));
      }
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
         if (dof < 0 && !interpolator_.HasExactMatch(s)) { continue; }

         // Negate slip for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t slip_dip =
            -interpolator_.EvaluateScalar(global_slip_dip, s);
         real_t slip_strike =
            -interpolator_.EvaluateScalar(global_slip_strike, s);

         real_t V_dip =
            std::abs(interpolator_.EvaluateScalar(global_V_dip, s));
         real_t V_strike =
            std::abs(interpolator_.EvaluateScalar(global_V_strike, s));

         // v55: tau_hat = tau_pre + elastic_traction + eta*V (Tandem convention)
         // V is signed from GetSlipRate; eta*V uses the raw signed values.
         real_t tau_dip =
            -(interpolator_.EvaluateScalar(tau_pre_dip_, s) +
              interpolator_.EvaluateScalar(global_trac_dip, s) +
              eta_ * interpolator_.EvaluateScalar(global_V_dip, s)) / 1e6;
         real_t tau_strike =
            -(interpolator_.EvaluateScalar(tau_pre_strike_, s) +
              interpolator_.EvaluateScalar(global_trac_strike, s) +
              eta_ * interpolator_.EvaluateScalar(global_V_strike, s)) / 1e6;

         real_t th = interpolator_.EvaluateScalar(global_theta, s);

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

   /// Write station-level traction decomposition from globally gathered data.
   void WriteTractionDecompositionFromGlobalData(
      real_t time,
      const Vector &global_trac_stress_dip,
      const Vector &global_trac_stress_strike,
      const Vector &global_trac_corr_dip,
      const Vector &global_trac_corr_strike)
   {
      if (decomp_probes_.empty()) { return; }

      for (int s = 0; s < interpolator_.NumStations(); s++)
      {
         int dof = interpolator_.GetNearestDOF(s);
         if (dof < 0 && !interpolator_.HasExactMatch(s)) { continue; }

         real_t tau_stress_dip =
            -interpolator_.EvaluateScalar(global_trac_stress_dip, s) / 1e6;
         real_t tau_stress_strike =
            -interpolator_.EvaluateScalar(global_trac_stress_strike, s) / 1e6;
         real_t tau_corr_dip =
            -interpolator_.EvaluateScalar(global_trac_corr_dip, s) / 1e6;
         real_t tau_corr_strike =
            -interpolator_.EvaluateScalar(global_trac_corr_strike, s) / 1e6;

         std::vector<real_t> row = {
            time,
            tau_stress_strike,
            tau_stress_dip,
            tau_corr_strike,
            tau_corr_dip,
            tau_stress_strike + tau_corr_strike,
            tau_stress_dip + tau_corr_dip
         };
         decomp_probes_[s]->WriteStep(row);
      }
   }

   /// Write station-level fault jump residual from globally gathered data.
   void WriteJumpResidualFromGlobalData(
      real_t time,
      const Vector &global_jump_res_dip,
      const Vector &global_jump_res_strike)
   {
      if (jump_res_probes_.empty()) { return; }

      for (int s = 0; s < interpolator_.NumStations(); s++)
      {
         int dof = interpolator_.GetNearestDOF(s);
         if (dof < 0 && !interpolator_.HasExactMatch(s)) { continue; }

         const real_t res_dip =
            interpolator_.EvaluateScalar(global_jump_res_dip, s);
         const real_t res_strike =
            interpolator_.EvaluateScalar(global_jump_res_strike, s);
         const real_t res_mag =
            std::sqrt(res_dip * res_dip + res_strike * res_strike);

         std::vector<real_t> row = {
            time,
            res_strike,
            res_dip,
            res_mag
         };
         jump_res_probes_[s]->WriteStep(row);
      }
   }

   /// Flush all output files.
   void Flush()
   {
      for (auto &p : probes_) { p->Flush(); }
      for (auto &p : decomp_probes_) { p->Flush(); }
      for (auto &p : jump_res_probes_) { p->Flush(); }
   }

   /// Close all output files.
   void Close()
   {
      for (auto &p : probes_) { p->Close(); }
      for (auto &p : decomp_probes_) { p->Close(); }
      for (auto &p : jump_res_probes_) { p->Close(); }
   }

   /// Number of probes.
   int NumProbes() const { return static_cast<int>(probes_.size()); }

   /// Access the interpolator (for diagnostics).
   const Probe2DInterpolator &GetInterpolator() const { return interpolator_; }

   /// Print station mapping diagnostics (delegates to interpolator).
   ///
   /// @param fault_x2 The same x2 coordinate vector used to construct the output
   /// @param fault_x3 The same x3 coordinate vector used to construct the output
   void PrintDiagnostics(const Vector &fault_x2, const Vector &fault_x3,
                         std::ostream &os = std::cout) const
   {
      interpolator_.PrintDiagnostics(stations_, fault_x2, fault_x3, os);

      // Check for duplicate faces
      int nbf = interpolator_.GetNbfPerFace();
      if (nbf >= 3)
      {
         os << "Checking for duplicate faces in gathered coordinates...\n";
         int dups = Probe2DInterpolator::CountDuplicateFaces(
            fault_x2, fault_x3, nbf, 1.0, &os);
         os << "  Total duplicate face pairs: " << dups << "\n";
         if (dups > 0)
         {
            os << "  WARNING: " << dups << " shared faces are duplicated "
               << "in the global gather. Exact interpolation may use stale "
               << "data from the wrong rank's copy.\n";
         }
      }
   }

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
   std::vector<std::unique_ptr<ProbeOutput>> decomp_probes_;
   std::vector<std::unique_ptr<ProbeOutput>> jump_res_probes_;
   real_t last_write_time_;

   // Pre-stress components for WriteFromGlobalData path
   Vector tau_pre_dip_;
   Vector tau_pre_strike_;

   // v55: Radiation damping coefficient for tau_hat output (matching Tandem)
   real_t eta_ = 0.0;

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
         if (dof < 0 && !interpolator_.HasExactMatch(s)) { continue; }

         // Internal: index 0 = dip, index 1 = strike
         // Negate slip for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         real_t slip_dip =
            -interpolator_.EvaluateInterleaved(slip, s, 0);
         real_t slip_strike =
            -interpolator_.EvaluateInterleaved(slip, s, 1);

         real_t V_dip =
            std::abs(interpolator_.EvaluateInterleaved(slip_rate, s, 0));
         real_t V_strike =
            std::abs(interpolator_.EvaluateInterleaved(slip_rate, s, 1));

         // v55: Total shear stress = tau_hat = tau_pre + elastic_traction + eta*V
         // Matches Tandem's DieterichRuinaBase::tau_hat (line 76):
         //   tau_hat = tau + TauPre + eta * V
         // Negate for SCEC output: internal convention uses negative
         // for right-lateral, SCEC expects positive.
         // V is signed (from GetSlipRate), eta*V uses absolute values.
         real_t tau_dip =
            -(interpolator_.EvaluateInterleaved(tau_pre, s, 0) +
              interpolator_.EvaluateInterleaved(traction, s, 0) +
              eta_ * interpolator_.EvaluateInterleaved(slip_rate, s, 0)) / 1e6;
         real_t tau_strike =
            -(interpolator_.EvaluateInterleaved(tau_pre, s, 1) +
              interpolator_.EvaluateInterleaved(traction, s, 1) +
              eta_ * interpolator_.EvaluateInterleaved(slip_rate, s, 1)) / 1e6;

         real_t th = interpolator_.EvaluateScalar(theta, s);

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
