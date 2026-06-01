// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/tpv104_stations.hpp — SCEC TPV104 (rate-and-state, slip law +
// strong-rate weakening) on-fault station-trace writer.  Extracted verbatim
// from `dynamic/tpv104_setup.hpp` so the SAME writer can be shared by the
// native `tpv104_driver` (which re-includes this header from
// tpv104_setup.hpp) AND by `seas_spatial_dyn_driver` (which includes this
// lean header directly, without dragging the native-driver fault-init /
// V_w side-channel / surface-station code into its translation unit).
// Mirrors the lean `dynamic/tpv31_stations.hpp`.
//
// One file per station; SCEC TPV104 trace layout (BP5 / Tandem canonical
// frame — h = strike = component 2, v = dip = component 1; TPV104 is pure
// strike-slip):
//
//   t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate
//      v-shear-stress  n-stress  psi

#ifndef MFEM_SEAS_TPV104_STATIONS_HPP
#define MFEM_SEAS_TPV104_STATIONS_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"            // DOFData
#include "../config/tpv104_params.hpp"    // kStationsTPV104

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

/// Station definition for TPV104 fault output.
/// Field names mirror the TPV102 pattern for drop-in consistency.
struct TPV104Station
{
   real_t along_strike;  ///< Along-strike coordinate x2 [m]
   real_t down_dip;      ///< Down-dip coordinate x3 [m] (positive distance)
   std::string name;     ///< SCEC trace-filename suffix (e.g. "x2_0_x3_7.5")
};

/// @brief Return the nine SCEC-required TPV104 fault stations.
///
/// Coordinates and SCEC trace-filename suffixes come from
/// `config/tpv104_params.hpp::kStationsTPV104`.  The nine (x2, x3) pairs
/// are the canonical TPV104 benchmark-trace stations used by Phase-3
/// probe-diff tooling (Step 13).
inline std::vector<TPV104Station> DefaultStations_TPV104()
{
   std::vector<TPV104Station> out;
   out.reserve(9);
   for (const auto &s : kStationsTPV104)
   {
      out.push_back({s.x2, s.x3, std::string(s.label)});
   }
   return out;
}

/// @brief Find the nearest fault DOF index to a given station location.
///
/// Linear scan over `fault_coords`; matches the TPV102 pattern.
///
/// @param[in] station       Station coordinates (along_strike, down_dip).
/// @param[in] fault_coords  Per-QP physical coordinates.
/// @param[in] ndof          Number of fault DOFs.
/// @return Index of the nearest DOF, or -1 if `ndof <= 0`.
inline int FindNearestDOF_TPV104(const TPV104Station &station,
                                 const std::vector<Vector> &fault_coords,
                                 int ndof)
{
   if (ndof <= 0) { return -1; }

   int best = 0;
   real_t best_dist = std::numeric_limits<real_t>::max();
   for (int i = 0; i < ndof; ++i)
   {
      const real_t dx = fault_coords[i](0) - station.along_strike;
      const real_t dz = std::abs(fault_coords[i](2)) - station.down_dip;
      const real_t dist2 = dx * dx + dz * dz;
      if (dist2 < best_dist)
      {
         best_dist = dist2;
         best = i;
      }
   }
   return best;
}

/// @brief Station output writer for TPV104 fault-trace data.
///
/// Writes one file per station.  Columns match the SCEC benchmark-trace
/// layout (§4.4 item 5):
///   t, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate,
///   v-shear-stress, n-stress, psi
/// Under BP5 / Tandem convention `horizontal = strike = component 2`
/// and `vertical = dip = component 1`, so the columns map to
///   t, slip2, V2, tau2_corr, slip1, V1, tau1_corr, sigma_n_corr, psi.
class TPV104StationWriter
{
public:
   /// Open station files (serial version — all DOFs are local).
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV104Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof)
   {
      stations_ = stations;
      station_dof_.resize(stations.size());
      station_owns_.assign(stations.size(), true);
      files_.resize(stations.size());

      for (size_t s = 0; s < stations.size(); ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV104(stations[s],
                                                 fault_coords, ndof);
         OpenStationFile(s, output_dir, prefix);
      }
   }

#ifdef MFEM_USE_MPI
   /// Open station files with MPI ownership resolution.  Only the rank
   /// with the globally nearest DOF opens the file (same pattern as
   /// `TPV102StationWriter::Open(... comm)`).
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV104Station> &stations,
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
         station_dof_[s] = FindNearestDOF_TPV104(stations[s],
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
         // R4-007 (review round 4): scale the tie-break tolerance by the
         // global-min distance so stations that happen to land far from
         // any QP (e.g. on a coarser mesh) still resolve a unique owner.
         // The previous absolute 1e-10 tolerance silently disqualified
         // every rank for a station whose |local - global_min| was
         // driven by mesh resolution, not by numerical round-off.
         const real_t tie_tol =
            std::max<real_t>(static_cast<real_t>(1e-10),
                             static_cast<real_t>(1e-9)
                               * std::max<real_t>(global_min_dist[s], 1.0));
         const bool is_candidate = std::abs(local_dist[s]
                                            - global_min_dist[s]) < tie_tol;
         const int candidate_rank = is_candidate ? my_rank : nprocs_loc;
         int winning_rank;
         MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT,
                       MPI_MIN, comm);
         // R4-007: MFEM_VERIFY that SOME rank won — otherwise the
         // station is silently dropped and the trace file never opens,
         // which only shows up at Phase-3 diff time.
         MFEM_VERIFY(winning_rank < nprocs_loc,
                     "TPV104StationWriter::Open: station "
                     << stations[s].name
                     << " has no winning rank (tie-break tolerance "
                     << tie_tol << " dropped every candidate). "
                     "Check global_min_dist = " << global_min_dist[s]);
         if (is_candidate && my_rank == winning_rank)
         {
            station_owns_[s] = true;
            OpenStationFile(s, output_dir, prefix);
         }
      }
   }
#endif

   /// Write a snapshot at time `t`.  SCEC trace column order.
   void WriteStep(real_t t, const std::vector<DOFData> &dof_data)
   {
      for (size_t s = 0; s < stations_.size(); ++s)
      {
         if (!station_owns_[s]) { continue; }
         const int idx = station_dof_[s];
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         if (!files_[s].is_open()) { continue; }

         const DOFData &d = dof_data[idx];

         // SCEC column order:  t, h-slip, h-slip-rate, h-shear-stress,
         //                     v-slip, v-slip-rate, v-shear-stress,
         //                     n-stress, psi.
         // h = strike = component 2, v = dip = component 1.
         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip2 << " "
                   << d.V2 << " "
                   << d.tau2_corr << " "
                   << d.slip1 << " "
                   << d.V1 << " "
                   << d.tau1_corr << " "
                   << d.sigma_n_corr << " "
                   << d.psi << "\n";
         // Flush each write: a Slurm SIGKILL at wall-time otherwise
         // discards the buffered probe data.
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
         files_[s] << "# TPV104 station: " << stations_[s].name
                   << " (along_strike=" << stations_[s].along_strike
                   << ", down_dip=" << stations_[s].down_dip << ")\n";
         files_[s] << "# Columns (SCEC TPV104 trace layout): "
                      "t  h-slip  h-slip-rate  h-shear-stress  "
                      "v-slip  v-slip-rate  v-shear-stress  "
                      "n-stress  psi\n";
         files_[s] << "# BP5 convention: h = strike = component 2, "
                      "v = dip = component 1.  TPV104 is pure strike-slip "
                      "so v-slip / v-slip-rate / v-shear-stress ≈ 0.\n";
      }
   }

   std::vector<TPV104Station> stations_;
   std::vector<int> station_dof_;
   std::vector<bool> station_owns_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_STATIONS_HPP
