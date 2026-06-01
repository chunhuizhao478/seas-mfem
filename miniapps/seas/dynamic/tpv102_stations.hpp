// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/tpv102_stations.hpp — SCEC TPV102 (rate-and-state) on-fault
// station-trace writer.  Extracted verbatim from `dynamic/tpv102_setup.hpp`
// so the SAME writer can be shared by the native `tpv102_driver` (which
// re-includes this header from tpv102_setup.hpp) AND by
// `seas_spatial_dyn_driver` (which includes this lean header directly,
// without dragging the native-driver fault-init / nucleation / surface-
// station code into its translation unit).  Mirrors the lean
// `dynamic/tpv31_stations.hpp`.
//
// One file per station; column layout (BP5 / Tandem canonical frame —
// 1 = dip = tangent1, 2 = strike = tangent2; TPV102 is pure strike-slip):
//
//   t  slip1  slip2  V1  V2  tau1  tau2  sigma_n  log10_theta
//
// The function names DefaultStations() / FindNearestDOF() are GENERIC (no
// _TPV102 suffix) because the native driver and tests call them by that
// name; preserved verbatim to keep native behaviour byte-identical.

#ifndef MFEM_SEAS_TPV102_STATIONS_HPP
#define MFEM_SEAS_TPV102_STATIONS_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"            // DOFData
#include "../config/tpv102_params.hpp"    // TPV102Params (psi -> theta)

#include <algorithm>
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

/// Station definition for fault output.
struct TPV102Station
{
   real_t along_strike;  ///< Along-strike coordinate [m]
   real_t down_dip;      ///< Down-dip coordinate [m]
   std::string name;     ///< Station identifier
};

/// @brief Return the 9 SCEC-required fault stations for TPV102.
///
/// SCEC spec requires time histories at these (along_strike, down_dip) locations.
inline std::vector<TPV102Station> DefaultStations()
{
   return {
      // 9 SCEC fault stations (x = along_strike [m], y = down_dip [m])
      {0.0,   3e3,   "flt_0_3"},
      {0.0,   7.5e3, "flt_0_7.5"},
      {0.0,   12e3,  "flt_0_12"},
      {9e3,   7.5e3, "flt_9_7.5"},
      {12e3,  3e3,   "flt_12_3"},
      {12e3,  12e3,  "flt_12_12"},
      {-9e3,  7.5e3, "flt_n9_7.5"},
      {-12e3, 3e3,   "flt_n12_3"},
      {-12e3, 12e3,  "flt_n12_12"},
   };
}

/// @brief Find the nearest fault DOF index to a given station location.
///
/// @param[in] station  Station coordinates (along_strike, down_dip).
/// @param[in] fault_coords  Fault DOF coordinates.
/// @param[in] ndof  Number of fault DOFs.
/// @return Index of the nearest DOF, or -1 if no DOFs.
inline int FindNearestDOF(const TPV102Station &station,
                          const std::vector<Vector> &fault_coords,
                          int ndof)
{
   if (ndof <= 0) { return -1; }

   int best = 0;
   real_t best_dist = std::numeric_limits<real_t>::max();

   for (int i = 0; i < ndof; i++)
   {
      real_t dx = fault_coords[i](0) - station.along_strike;
      real_t dz = std::abs(fault_coords[i](2)) - station.down_dip;
      real_t dist2 = dx*dx + dz*dz;
      if (dist2 < best_dist)
      {
         best_dist = dist2;
         best = i;
      }
   }
   return best;
}

/// @brief Station output writer for TPV102 fault data.
///
/// Writes time-series of (time, slip1, slip2, slip_rate, tau1, tau2, log10_theta)
/// for each station. One file per station.
class TPV102StationWriter
{
public:
   /// Open station output files.
   /// @param[in] output_dir  Directory for output files.
   /// @param[in] prefix  File name prefix.
   /// @param[in] stations  Station definitions.
   /// @param[in] fault_coords  Fault DOF coordinates.
   /// @param[in] ndof  Number of fault DOFs.
   /// Open station files (serial version — all DOFs are local).
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV102Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof)
   {
      stations_ = stations;
      station_dof_.resize(stations.size());
      station_owns_.assign(stations.size(), true);
      files_.resize(stations.size());

      for (size_t s = 0; s < stations.size(); s++)
      {
         station_dof_[s] = FindNearestDOF(stations[s], fault_coords, ndof);
         OpenStationFile(s, output_dir, prefix);
      }
   }

#ifdef MFEM_USE_MPI
   /// Open station files with MPI ownership resolution (R-004 fix).
   /// Only the rank with the globally nearest DOF opens the file.
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV102Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof, MPI_Comm comm)
   {
      stations_ = stations;
      int nstations = static_cast<int>(stations.size());
      station_dof_.resize(nstations);
      station_owns_.assign(nstations, false);
      files_.resize(nstations);

      // Compute local distance to each station
      std::vector<real_t> local_dist(nstations,
                                     std::numeric_limits<real_t>::max());
      for (int s = 0; s < nstations; s++)
      {
         station_dof_[s] = FindNearestDOF(stations[s], fault_coords, ndof);
         if (station_dof_[s] >= 0 && station_dof_[s] < ndof)
         {
            real_t dx = fault_coords[station_dof_[s]](0) - stations[s].along_strike;
            real_t dz = std::abs(fault_coords[station_dof_[s]](2))
                      - stations[s].down_dip;
            local_dist[s] = std::sqrt(dx*dx + dz*dz);
         }
      }

      // Global min distance across all ranks.  R-401 fix: match real_t at
      // compile time via MPITypeMap (same pattern used at all R-303 sites).
      // Hardcoded MPI_DOUBLE silently corrupts local_dist/global_min_dist
      // on MFEM_USE_SINGLE builds (real_t = float, stride 4 bytes).
      std::vector<real_t> global_min_dist(nstations);
      MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);

      // R-005 fix: Tiebreaker — among equidistant ranks, lowest rank ID wins.
      // Prevents multiple ranks writing the same station file.
      int my_rank, nprocs_loc;
      MPI_Comm_rank(comm, &my_rank);
      MPI_Comm_size(comm, &nprocs_loc);
      for (int s = 0; s < nstations; s++)
      {
         bool is_candidate = std::abs(local_dist[s] - global_min_dist[s]) < 1e-10;
         int candidate_rank = is_candidate ? my_rank : nprocs_loc;
         int winning_rank;
         MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT, MPI_MIN, comm);
         if (is_candidate && my_rank == winning_rank)
         {
            station_owns_[s] = true;
            OpenStationFile(s, output_dir, prefix);
         }
      }
   }
#endif

   /// Write a snapshot at the current time.
   void WriteStep(real_t t, const std::vector<DOFData> &dof_data)
   {
      for (size_t s = 0; s < stations_.size(); s++)
      {
         if (!station_owns_[s]) { continue; }
         int idx = station_dof_[s];
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         if (!files_[s].is_open()) { continue; }

         const DOFData &d = dof_data[idx];

         // Convert psi to theta: theta = (Dc/V0) * exp((psi - f0) / b)
         real_t theta = (d.Dc / TPV102Params::V0)
                      * std::exp((d.psi - TPV102Params::f0) / TPV102Params::b);
         real_t log10_theta = std::log10(std::max(theta, 1e-300));

         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip1 << " "
                   << d.slip2 << " "
                   << d.V1 << " "
                   << d.V2 << " "
                   << d.tau1_corr << " "
                   << d.tau2_corr << " "
                   << d.sigma_n_corr << " "
                   << log10_theta << "\n";
         // Flush on every write: output_interval makes this ~30/rank/hour,
         // and a Slurm SIGKILL at the wall-time limit otherwise discards
         // all buffered probe data (v1 debug: all .dat files ended up 0 B).
         files_[s].flush();
      }
   }

   /// Flush all open file handles.
   void Flush()
   {
      for (auto &f : files_) { if (f.is_open()) { f.flush(); } }
   }

   /// Close all file handles.
   void Close()
   {
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   void OpenStationFile(size_t s, const std::string &output_dir,
                        const std::string &prefix)
   {
      std::string fname = output_dir + "/" + prefix + "_station_"
                        + stations_[s].name + ".dat";
      files_[s].open(fname);
      if (files_[s].is_open())
      {
         files_[s] << "# TPV102 station: " << stations_[s].name
                   << " (along_strike=" << stations_[s].along_strike
                   << ", down_dip=" << stations_[s].down_dip << ")\n";
         files_[s] << "# Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta\n";
         files_[s] << "# BP5 convention: 1 = dip (tangent1), 2 = strike (tangent2). TPV102 is pure strike-slip so V1/slip1/tau1 ~ 0 "
                      "and the interesting physics is in column 2.\n";
      }
   }

   std::vector<TPV102Station> stations_;
   std::vector<int> station_dof_;
   std::vector<bool> station_owns_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_STATIONS_HPP
