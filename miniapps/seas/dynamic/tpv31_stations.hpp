// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/tpv31_stations.hpp — SCEC TPV31 on-fault station-trace writer
// for `seas_spatial_dyn_driver`.
//
// Mirrors `TPV205StationWriter` (dynamic/tpv205_setup.hpp): one file per
// station, same 9-column internal layout consumed by
// `tpv31/visualize_results.py` (MFEM vs SCEC EQdyna / SeisSol overlay):
//
//   t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate
//      v-shear-stress  n-stress  mu_eff
//
// where (canonical SEAS frame, shared with TPV205/102/104):
//   h = horizontal = along-strike = component 2  (slip2 / V2 / tau2_corr)
//   v = vertical   = along-dip    = component 1  (slip1 / V1 / tau1_corr)
//   n-stress       = sigma_n_corr  (Pa, compression POSITIVE — internal
//                    convention; NOT the SCEC extension-positive sign)
//   mu_eff         = LSW μ(δ) at δ = sqrt(slip1² + slip2²)
//
// Stresses are written in raw Pa (not MPa) to match the TPV205 writer
// and the comparison harness (which compares columns directly, no unit
// conversion).  TPV31 is pure right-lateral strike-slip, so the v-*
// (dip) channels stay near zero.
//
// The 30 SCEC on-fault stations (TPV31_32_Description_v03.pdf Part 5)
// are the 3 × 10 grid of along-strike {0, 6, 12} km × down-dip
// {0, 0.2, 0.5, 1.0, 2.4, 3.0, 5.0, 7.5, 10.0, 12.0} km positions.
// The TPV31 TOML is rotated into the canonical SEAS frame, so a station
// at spec (along-strike S, depth D, 0) lives at mesh (S, 0, -D); the
// nearest-DOF search uses |z| as the down-dip depth (matching the
// TPV205 writer).

#ifndef MFEM_SEAS_TPV31_STATIONS_HPP
#define MFEM_SEAS_TPV31_STATIONS_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "station_nearest_tiebreak.hpp"   // FindNearestFaultDOFLex / ResolveStationOwnerLex   // DOFData

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Station definition for TPV31 fault output (mirrors TPV205Station).
struct TPV31Station
{
   real_t along_strike;  ///< Along-strike coordinate x [m]
   real_t down_dip;      ///< Down-dip coordinate (depth) [m] (positive)
   std::string name;     ///< SCEC trace-filename suffix (faultst###dp###)
};

/// @brief Effective LSW friction coefficient at slip magnitude `delta`.
///
/// μ(δ) = μ_s − (μ_s − μ_d) · min(δ/d_c, 1).  Generic (no TPV205
/// strength-barrier sentinel — TPV31 has no μ_s barrier; slip is pinned
/// to zero at the fault border by the slip boundary condition instead).
inline real_t LSWFrictionCoefficient_TPV31(real_t delta,
                                           real_t mu_s, real_t mu_d,
                                           real_t d_c)
{
   if (d_c <= 0.0)   { return mu_s; }
   if (delta <= 0.0) { return mu_s; }
   if (delta >= d_c) { return mu_d; }
   return mu_s - (mu_s - mu_d) * (delta / d_c);
}

/// @brief Return the 30 SCEC TPV31 on-fault stations (spec Part 5).
inline std::vector<TPV31Station> DefaultStations_TPV31()
{
   const real_t strikes_m[3]  = { 0.0, 6000.0, 12000.0 };
   const real_t depths_m[10]  = { 0.0, 200.0, 500.0, 1000.0, 2400.0,
                                  3000.0, 5000.0, 7500.0, 10000.0, 12000.0 };

   std::vector<TPV31Station> out;
   out.reserve(30);
   char buf[32];
   for (real_t st : strikes_m)
   {
      for (real_t dp : depths_m)
      {
         // SCEC naming: 3-digit zero-padded hundreds-of-metres code.
         const int st_code = static_cast<int>(std::lround(st / 100.0));
         const int dp_code = static_cast<int>(std::lround(dp / 100.0));
         std::snprintf(buf, sizeof(buf), "faultst%03ddp%03d",
                       st_code, dp_code);
         out.push_back({ st, dp, std::string(buf) });
      }
   }
   return out;
}

/// @brief Find the nearest fault DOF index to a given station location.
///
/// Matches the TPV205 search: dx = x − along_strike, dz = |z| − down_dip
/// (depth is |z| in the canonical mesh frame, z < 0 below surface).
/// Delegates to `FindNearestFaultDOFLex` (station_nearest_tiebreak.hpp):
/// same distance metric as the historical scan, plus the deterministic
/// lexicographic (x, z, y) tie-break for stations equidistant from several
/// QPs (np4_attractor_root_cause_2026-07-10.md).  Non-tied stations
/// resolve to the identical QP as before.
inline int FindNearestDOF_TPV31(const TPV31Station &station,
                                const std::vector<Vector> &fault_coords,
                                int ndof)
{
   return FindNearestFaultDOFLex(fault_coords, ndof,
                                 station.along_strike, station.down_dip);
}

/// @brief Station output writer for TPV31 fault-trace data.
///
/// One file per station, same column layout as TPV205StationWriter (the
/// trailing channel is the LSW μ_eff(δ)).  TPV31 is right-lateral
/// strike-slip; the v-* (dip) channels remain near 0.
class TPV31StationWriter
{
public:
   /// Open station files (serial — all DOFs are local).
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV31Station> &stations,
             const std::vector<Vector> &fault_coords,
             int ndof)
   {
      stations_ = stations;
      station_dof_.resize(stations.size());
      station_owns_.assign(stations.size(), true);
      files_.resize(stations.size());

      for (size_t s = 0; s < stations.size(); ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV31(stations[s],
                                                fault_coords, ndof);
         OpenStationFile(s, output_dir, prefix);
      }
   }

#ifdef MFEM_USE_MPI
   /// Open station files with MPI ownership resolution (the owning rank
   /// is the one whose nearest DOF is globally closest; ties broken by
   /// lowest rank).  Mirrors TPV205StationWriter::Open.
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV31Station> &stations,
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
         station_dof_[s] = FindNearestDOF_TPV31(stations[s],
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
                                    "TPV31StationWriter",
                                    stations[s].name))
         {
            station_owns_[s] = true;
            OpenStationFile(s, output_dir, prefix);
         }
      }
   }
#endif

   /// Write a snapshot at time `t`.  Column order matches the TPV205
   /// writer, with the trailing channel being μ_eff(δ) for LSW.
   void WriteStep(real_t t, const std::vector<DOFData> &dof_data)
   {
      for (size_t s = 0; s < stations_.size(); ++s)
      {
         if (!station_owns_[s]) { continue; }
         const int idx = station_dof_[s];
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         if (!files_[s].is_open()) { continue; }

         const DOFData &d = dof_data[idx];

         const real_t delta = std::sqrt(d.slip1 * d.slip1
                                        + d.slip2 * d.slip2);
         const real_t mu_eff =
            LSWFrictionCoefficient_TPV31(delta,
                                         d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);

         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip2 << " "        // h-slip       (strike)
                   << d.V2 << " "           // h-slip-rate
                   << d.tau2_corr << " "    // h-shear-stress
                   << d.slip1 << " "        // v-slip       (dip)
                   << d.V1 << " "           // v-slip-rate
                   << d.tau1_corr << " "    // v-shear-stress
                   << d.sigma_n_corr << " " // n-stress (Pa, compression+)
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
         files_[s] << "# TPV31 station: " << stations_[s].name
                   << " (along_strike=" << stations_[s].along_strike
                   << ", down_dip=" << stations_[s].down_dip << ")\n";
         files_[s] << "# Columns (SCEC TPV31 trace layout, LSW): "
                      "t  h-slip  h-slip-rate  h-shear-stress  "
                      "v-slip  v-slip-rate  v-shear-stress  "
                      "n-stress  mu_eff\n";
         files_[s] << "# Canonical frame: h = strike = component 2, "
                      "v = dip = component 1.  TPV31 is right-lateral "
                      "strike-slip; v-* channels remain near 0.  "
                      "Stresses in Pa, n-stress compression-positive; "
                      "mu_eff is the LSW μ(δ) at δ = sqrt(slip1² + slip2²).\n";
      }
   }

   std::vector<TPV31Station> stations_;
   std::vector<int> station_dof_;
   std::vector<bool> station_owns_;
   std::vector<std::ofstream> files_;
};

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_TPV31_STATIONS_HPP
