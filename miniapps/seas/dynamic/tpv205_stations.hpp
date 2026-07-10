// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/tpv205_stations.hpp — SCEC TPV205 (LSW) on-fault station-trace
// writer.  Extracted verbatim from `dynamic/tpv205_setup.hpp` so the SAME
// writer can be shared by the native `tpv205_driver` (which re-includes
// this header from tpv205_setup.hpp) AND by `seas_spatial_dyn_driver`
// (which includes this lean header directly, without dragging the
// native-driver fault-init / surface-station code into its translation
// unit).  Mirrors the lean `dynamic/tpv31_stations.hpp`.
//
// One file per station; 9-column SCEC TPV5 trace layout:
//
//   t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate
//      v-shear-stress  n-stress  mu_eff
//
// where (canonical SEAS / BP5 frame, shared with TPV31/102/104):
//   h = horizontal = along-strike = component 2  (slip2 / V2 / tau2_corr)
//   v = vertical   = along-dip    = component 1  (slip1 / V1 / tau1_corr)
//   n-stress       = sigma_n_corr  (Pa, compression POSITIVE — internal
//                    convention; NOT the SCEC extension-positive sign)
//   mu_eff         = LSW μ(δ) at δ = sqrt(slip1² + slip2²)
//
// TPV205 is right-lateral strike-slip, so the v-* (dip) channels stay
// near zero.

#ifndef MFEM_SEAS_TPV205_STATIONS_HPP
#define MFEM_SEAS_TPV205_STATIONS_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"            // DOFData
#include "station_nearest_tiebreak.hpp"   // FindNearestFaultDOFLex / ResolveStationOwnerLex
#include "tpv205_friction.hpp"            // LSWFrictionCoefficient_TPV205
#include "../config/tpv205_params.hpp"    // StationTPV205, kStationsTPV205

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Station definition for TPV205 fault output (mirrors TPV104Station).
struct TPV205Station
{
   real_t along_strike;  ///< Along-strike coordinate x2 [m]
   real_t down_dip;      ///< Down-dip coordinate x3 [m] (positive)
   std::string name;     ///< SCEC trace-filename suffix
};

/// @brief Return the 16 SCEC TPV205 on-fault stations.
inline std::vector<TPV205Station> DefaultStations_TPV205()
{
   std::vector<TPV205Station> out;
   out.reserve(16);
   for (const auto &s : kStationsTPV205)
   {
      out.push_back({s.x2, s.x3, std::string(s.label)});
   }
   return out;
}

/// @brief Find the nearest fault DOF index to a given station location.
/// Delegates to `FindNearestFaultDOFLex` (station_nearest_tiebreak.hpp):
/// same distance metric as the historical scan, plus the deterministic
/// lexicographic (x, z, y) tie-break for stations equidistant from several
/// QPs (np4_attractor_root_cause_2026-07-10.md).  Non-tied stations
/// resolve to the identical QP as before.
inline int FindNearestDOF_TPV205(const TPV205Station &station,
                                 const std::vector<Vector> &fault_coords,
                                 int ndof)
{
   return FindNearestFaultDOFLex(fault_coords, ndof,
                                 station.along_strike, station.down_dip);
}

/// @brief Station output writer for TPV205 fault-trace data.
///
/// Writes one file per station.  Columns match the SCEC benchmark-trace
/// layout used by TPV102/TPV104 (BP5-canonical-frame remap):
///   t, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate,
///   v-shear-stress, n-stress, μ_eff
/// where h = strike = component 2, v = dip = component 1, and the
/// trailing column is the LSW effective friction coefficient
///   μ_eff(δ) = μ_s − (μ_s − μ_d) · min(δ/d_c, 1)
/// at the per-QP slip magnitude δ = sqrt(slip1² + slip2²).  μ_eff
/// replaces the rate-and-state ψ column used in TPV104 traces.
class TPV205StationWriter
{
public:
   /// Open station files (serial — all DOFs are local).
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV205Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof)
   {
      stations_ = stations;
      station_dof_.resize(stations.size());
      station_owns_.assign(stations.size(), true);
      files_.resize(stations.size());

      for (size_t s = 0; s < stations.size(); ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV205(stations[s],
                                                 fault_coords, ndof);
         OpenStationFile(s, output_dir, prefix);
      }
   }

#ifdef MFEM_USE_MPI
   /// Open station files with MPI ownership resolution.
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV205Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof, MPI_Comm comm)
   {
      stations_ = stations;
      const int nstations = static_cast<int>(stations.size());
      station_dof_.resize(nstations);
      station_owns_.assign(nstations, false);
      files_.resize(nstations);

      std::vector<real_t> local_dist(nstations,
                                     std::numeric_limits<real_t>::max());
      for (int s = 0; s < nstations; ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV205(stations[s],
                                                 fault_coords, ndof);
         if (station_dof_[s] >= 0 && station_dof_[s] < ndof)
         {
            const real_t dx = fault_coords[station_dof_[s]](0)
                              - stations[s].along_strike;
            const real_t dz = std::abs(fault_coords[station_dof_[s]](2))
                              - stations[s].down_dip;
            local_dist[s] = std::sqrt(dx * dx + dz * dz);
         }
      }

      std::vector<real_t> global_min_dist(nstations);
      MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);

      int my_rank, nprocs_loc;
      MPI_Comm_rank(comm, &my_rank);
      MPI_Comm_size(comm, &nprocs_loc);
      for (int s = 0; s < nstations; ++s)
      {
         const real_t tie_tol =
            std::max<real_t>(static_cast<real_t>(1e-10),
                             static_cast<real_t>(1e-9)
                               * std::max<real_t>(global_min_dist[s], 1.0));
         // R-201 anchored re-scan (REVIEW_station_tiebreak_2026-07-10.md):
         // candidates ≡ ranks owning a QP within tie_tol of the GLOBAL
         // minimum distance; the submitted pick is the local lexicographic
         // minimum over exactly that set — partition-invariant even in the
         // osculating band at the window edge.  Owner = rank holding the
         // lexicographically smallest candidate COORDINATE
         // (np4_attractor_root_cause_2026-07-10.md; the fail-loud
         // winning-rank VERIFY lives inside ResolveStationOwnerLex).
         const int lex_dof = FindNearestFaultDOFLex(
            fault_coords, ndof, stations[s].along_strike,
            stations[s].down_dip, global_min_dist[s], tie_tol);
         if (lex_dof >= 0) { station_dof_[s] = lex_dof; }
         const bool is_candidate = lex_dof >= 0;
         const real_t RMAX = std::numeric_limits<real_t>::max();
         const real_t bx = is_candidate ? fault_coords[station_dof_[s]](0)
                                        : RMAX;
         const real_t bz = is_candidate ? fault_coords[station_dof_[s]](2)
                                        : RMAX;
         const real_t by = is_candidate ? fault_coords[station_dof_[s]](1)
                                        : RMAX;
         if (ResolveStationOwnerLex(comm, is_candidate, bx, bz, by,
                                    "TPV205StationWriter",
                                    stations[s].name))
         {
            station_owns_[s] = true;
            OpenStationFile(s, output_dir, prefix);
         }
      }
   }
#endif

   /// Write a snapshot at time `t`.  SCEC trace column order, with the
   /// trailing channel being μ_eff(δ) for LSW (instead of ψ for R&S).
   void WriteStep(real_t t, const std::vector<DOFData> &dof_data)
   {
      for (size_t s = 0; s < stations_.size(); ++s)
      {
         if (!station_owns_[s]) { continue; }
         const int idx = station_dof_[s];
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         if (!files_[s].is_open()) { continue; }

         const DOFData &d = dof_data[idx];

         // R-017 / R-016: route through the canonical LSW helper so the
         // strength-barrier short-circuit (R-002) applies uniformly
         // across station traces, ParaView, and the friction solve.
         // Reads the LSW-native fields directly — no repurposing.
         const real_t delta = std::sqrt(d.slip1 * d.slip1
                                        + d.slip2 * d.slip2);
         const real_t mu_eff =
            LSWFrictionCoefficient_TPV205(delta,
                                          d.lsw_mu_s, d.lsw_mu_d,
                                          d.lsw_d_c);

         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip2 << " "
                   << d.V2 << " "
                   << d.tau2_corr << " "
                   << d.slip1 << " "
                   << d.V1 << " "
                   << d.tau1_corr << " "
                   << d.sigma_n_corr << " "
                   << mu_eff << "\n";
         files_[s].flush();
      }
   }

   void Flush()
   {
      for (auto &f : files_) { if (f.is_open()) { f.flush(); } }
   }

   void Close()
   {
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   void OpenStationFile(size_t s, const std::string &output_dir,
                        const std::string &prefix)
   {
      const std::string fname = output_dir + "/" + prefix + "_station_"
                                + stations_[s].name + ".dat";
      files_[s].open(fname);
      if (files_[s].is_open())
      {
         files_[s] << "# TPV205 station: " << stations_[s].name
                   << " (along_strike=" << stations_[s].along_strike
                   << ", down_dip=" << stations_[s].down_dip << ")\n";
         files_[s] << "# Columns (SCEC TPV5 trace layout, LSW): "
                      "t  h-slip  h-slip-rate  h-shear-stress  "
                      "v-slip  v-slip-rate  v-shear-stress  "
                      "n-stress  mu_eff\n";
         files_[s] << "# BP5 convention: h = strike = component 2, "
                      "v = dip = component 1.  TPV205 is right-lateral "
                      "strike-slip; v-* channels remain near 0.  "
                      "mu_eff is the LSW μ(δ) at the current slip "
                      "magnitude δ = sqrt(slip1² + slip2²).\n";
      }
   }

   std::vector<TPV205Station> stations_;
   std::vector<int> station_dof_;
   std::vector<bool> station_owns_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV205_STATIONS_HPP
