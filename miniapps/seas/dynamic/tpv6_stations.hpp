// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// SCEC TPV6/TPV7 per-side ON-FAULT station writer (Part C / C2).
//
// Per the TPV6/7 spec (Harris/Day et al. 2007, Part I) the on-fault split-node
// stations report, on EACH side of the bi-material fault, the absolute particle
// DISPLACEMENT and VELOCITY (NOT slip / slip-rate), plus the fault traction.
// One pair of files per station — nearside (STRONG/fast material) and farside
// (WEAK/slow material) — matching the DRDG3D reference layout in
// tpv6/benchmark_data/scec_drdg3d so tpv6/visualize_results.py overlays directly.
//
// Columns (MKS; one row per output step):
//   t  h-disp  h-vel  h-stress  v-disp  v-vel  v-stress  n-disp  n-vel  n-stress
//   h = along-strike, v = along-dip (down-dip), n = fault-normal.
//   Velocities m/s, displacements m, STRESSES IN MPa.
//
// Sign / side conventions (see tpv6/benchmark_data/README.md):
//   * nearside = the STRONG side (larger Zp); farside = WEAK side.  Resolved per
//     DOF from DOFData::Zp_{plus,minus} so it is correct regardless of the
//     config's n_y (which physical y-side the strong material is on).
//   * n-stress is COMPRESSION-POSITIVE (this repo's convention; the DRDG3D
//     reference is compression-negative — visualize_results.py flips the reference
//     so both align).
//   * Per-side velocity = the IMPOSED (split-node) Godunov velocity captured into
//     DOFData::v_imp_{plus,minus} by FaultFaceFlux::Evaluate*; displacement is the
//     trapezoidal time integral d += 0.5*(v_prev + v)*dt accumulated here.
//
// Restart NOTE: like the other tpv*_stations writers, Open() truncates the trace
// files, so the displacement accumulator restarts from 0 on a resumed run (the
// trace covers only the resumed segment).  TPV6 is short (single-segment), so this
// is acceptable; full restart-continuity (append + last-row read) is a follow-up.

#ifndef MFEM_SEAS_TPV6_STATIONS_HPP
#define MFEM_SEAS_TPV6_STATIONS_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"   // DOFData

#include <array>
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

/// On-fault station for TPV6/7 (split-node, per-side output).
struct TPV6Station
{
   real_t along_strike;  ///< along-strike x [m]
   real_t down_dip;      ///< down-dip depth [m] (positive)
   std::string name;     ///< filename id, e.g. "x2_0_x3_0" (matches drdg3d)
};

/// @brief The 5 SCEC TPV6 on-fault stations (drdg3d reference set):
///   (strike, depth) km = (0,0), (-12,0), (12,0), (-12,7.5), (12,7.5).
inline std::vector<TPV6Station> DefaultStations_TPV6()
{
   return {
      { 0.0,       0.0,    "x2_0_x3_0" },
      { -12000.0,  0.0,    "x2_-12_x3_0" },
      {  12000.0,  0.0,    "x2_12_x3_0" },
      { -12000.0,  7500.0, "x2_-12_x3_7.5" },
      {  12000.0,  7500.0, "x2_12_x3_7.5" },
   };
}

/// @brief Nearest fault DOF to a station (dx = x - strike, dz = |z| - depth).
inline int FindNearestDOF_TPV6(const TPV6Station &station,
                               const std::vector<Vector> &fault_coords, int ndof)
{
   if (ndof <= 0) { return -1; }
   int best = 0;
   real_t best_dist = std::numeric_limits<real_t>::max();
   for (int i = 0; i < ndof; ++i)
   {
      const real_t dx = fault_coords[i](0) - station.along_strike;
      const real_t dz = std::abs(fault_coords[i](2)) - station.down_dip;
      const real_t dist2 = dx * dx + dz * dz;
      if (dist2 < best_dist) { best_dist = dist2; best = i; }
   }
   return best;
}

/// @brief Per-side on-fault station writer (two files per station).
class TPV6StationWriter
{
public:
   void Open(const std::string &output_dir, const std::string &prefix,
             const std::vector<TPV6Station> &stations,
             const std::vector<Vector> &fault_coords, int ndof)
   {
      Init_(stations);
      for (size_t s = 0; s < stations.size(); ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV6(stations[s], fault_coords, ndof);
         station_owns_[s] = true;
         OpenStationFiles_(s, output_dir, prefix);
      }
   }

#ifdef MFEM_USE_MPI
   void Open(const std::string &output_dir, const std::string &prefix,
             const std::vector<TPV6Station> &stations,
             const std::vector<Vector> &fault_coords, int ndof, MPI_Comm comm)
   {
      Init_(stations);
      const int n = static_cast<int>(stations.size());
      std::vector<real_t> local_dist(n, std::numeric_limits<real_t>::max());
      for (int s = 0; s < n; ++s)
      {
         station_dof_[s] = FindNearestDOF_TPV6(stations[s], fault_coords, ndof);
         if (station_dof_[s] >= 0 && station_dof_[s] < ndof)
         {
            const real_t dx = fault_coords[station_dof_[s]](0) - stations[s].along_strike;
            const real_t dz = std::abs(fault_coords[station_dof_[s]](2)) - stations[s].down_dip;
            local_dist[s] = std::sqrt(dx * dx + dz * dz);
         }
      }
      std::vector<real_t> gmin(n);
      MPI_Allreduce(local_dist.data(), gmin.data(), n,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
      int my_rank, nproc;
      MPI_Comm_rank(comm, &my_rank); MPI_Comm_size(comm, &nproc);
      for (int s = 0; s < n; ++s)
      {
         const real_t tol = std::max<real_t>(real_t(1e-10),
                              real_t(1e-9) * std::max<real_t>(gmin[s], 1.0));
         const bool cand = std::abs(local_dist[s] - gmin[s]) < tol;
         const int cand_rank = cand ? my_rank : nproc;
         int win; MPI_Allreduce(&cand_rank, &win, 1, MPI_INT, MPI_MIN, comm);
         MFEM_VERIFY(win < nproc, "TPV6StationWriter::Open: station "
                     << stations[s].name << " has no owning rank.");
         if (cand && my_rank == win)
         {
            station_owns_[s] = true;
            OpenStationFiles_(s, output_dir, prefix);
         }
      }
   }
#endif

   /// Write one snapshot.  Per-side velocity from DOFData::v_imp_{plus,minus};
   /// displacement integrated trapezoidally between calls.
   void WriteStep(real_t t, const std::vector<DOFData> &dof_data)
   {
      const real_t dt = have_prev_ ? (t - t_prev_) : 0.0;
      for (size_t s = 0; s < stations_.size(); ++s)
      {
         if (!station_owns_[s]) { continue; }
         const int idx = station_dof_[s];
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         const DOFData &d = dof_data[idx];

         // near = STRONG (larger Zp); far = WEAK.  Robust to the n_y choice.
         const bool minus_is_near = (d.Zp_minus >= d.Zp_plus);
         const real_t *vn = minus_is_near ? d.v_imp_minus : d.v_imp_plus;
         const real_t *vf = minus_is_near ? d.v_imp_plus  : d.v_imp_minus;

         // Trapezoidal displacement integral, per side, per component
         // (component order [0]=normal, [1]=dip, [2]=strike).
         for (int k = 0; k < 3; ++k)
         {
            near_disp_[s][k] += 0.5 * (near_vprev_[s][k] + vn[k]) * dt;
            far_disp_[s][k]  += 0.5 * (far_vprev_[s][k]  + vf[k]) * dt;
            near_vprev_[s][k] = vn[k];
            far_vprev_[s][k]  = vf[k];
         }

         // Single-valued fault traction (Pa -> MPa), compression-positive sigma_n.
         const real_t h_stress = d.tau2_corr   / 1.0e6;   // strike
         const real_t v_stress = d.tau1_corr   / 1.0e6;   // dip
         const real_t n_stress = d.sigma_n_corr / 1.0e6;  // normal (compression +)

         WriteRow_(files_near_[s], t, near_disp_[s], vn, h_stress, v_stress, n_stress);
         WriteRow_(files_far_[s],  t, far_disp_[s],  vf, h_stress, v_stress, n_stress);
      }
      t_prev_ = t;
      have_prev_ = true;
   }

   void Flush()
   {
      for (auto &f : files_near_) { if (f.is_open()) { f.flush(); } }
      for (auto &f : files_far_)  { if (f.is_open()) { f.flush(); } }
   }
   void Close()
   {
      for (auto &f : files_near_) { if (f.is_open()) { f.close(); } }
      for (auto &f : files_far_)  { if (f.is_open()) { f.close(); } }
   }

private:
   void Init_(const std::vector<TPV6Station> &stations)
   {
      stations_ = stations;
      const size_t n = stations.size();
      station_dof_.assign(n, -1);
      station_owns_.assign(n, false);
      files_near_.resize(n);
      files_far_.resize(n);
      near_disp_.assign(n, {0.0, 0.0, 0.0});
      far_disp_.assign(n, {0.0, 0.0, 0.0});
      near_vprev_.assign(n, {0.0, 0.0, 0.0});
      far_vprev_.assign(n, {0.0, 0.0, 0.0});
      t_prev_ = 0.0;
      have_prev_ = false;
   }

   // disp/v component order: [0]=normal, [1]=dip, [2]=strike.  Output column
   // order: t, h(strike), v(dip), n(normal) triples.
   void WriteRow_(std::ofstream &f, real_t t, const std::array<real_t, 3> &disp,
                  const real_t *vel, real_t hs, real_t vs, real_t ns)
   {
      if (!f.is_open()) { return; }
      f << std::scientific << std::setprecision(10)
        << t << " "
        << disp[2] << " " << vel[2] << " " << hs << " "   // h = strike
        << disp[1] << " " << vel[1] << " " << vs << " "   // v = dip
        << disp[0] << " " << vel[0] << " " << ns << "\n"; // n = normal
      f.flush();
   }

   void OpenStationFiles_(size_t s, const std::string &output_dir,
                          const std::string &prefix)
   {
      OpenOne_(files_near_[s], output_dir, prefix, "nearside", s);
      OpenOne_(files_far_[s],  output_dir, prefix, "farside",  s);
   }
   void OpenOne_(std::ofstream &f, const std::string &output_dir,
                 const std::string &prefix, const char *side, size_t s)
   {
      const std::string fname = output_dir + "/" + prefix + "_" + side + "_"
                                + stations_[s].name + ".dat";
      f.open(fname);
      if (f.is_open())
      {
         f << "# TPV6/7 station " << stations_[s].name << " " << side
           << " (along_strike=" << stations_[s].along_strike
           << ", down_dip=" << stations_[s].down_dip << ")\n";
         f << "# nearside = STRONG (larger Zp) side; farside = WEAK side.\n";
         f << "# Columns: t h-disp h-vel h-stress v-disp v-vel v-stress "
              "n-disp n-vel n-stress\n";
         f << "# h=strike, v=dip, n=normal; vel m/s, disp m, stress MPa, "
              "n-stress COMPRESSION-POSITIVE (drdg3d ref is negative — the viz "
              "flips it).  Velocity = imposed split-node Godunov velocity; "
              "disp = trapezoidal time-integral.\n";
      }
   }

   std::vector<TPV6Station> stations_;
   std::vector<int> station_dof_;
   std::vector<bool> station_owns_;
   std::vector<std::ofstream> files_near_, files_far_;
   std::vector<std::array<real_t, 3>> near_disp_, far_disp_;
   std::vector<std::array<real_t, 3>> near_vprev_, far_vprev_;
   real_t t_prev_ = 0.0;
   bool have_prev_ = false;
};

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_TPV6_STATIONS_HPP
