// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV5 / TPV205 benchmark setup: fault-DOF init, patch-dependent
// pre-stress, station output writers.
//
// Mirrors `dynamic/tpv104_setup.hpp` with these LSW-specific changes:
//   - Per-QP `μ_s(x, z)` via `ComputeMuS_TPV205` (0.677 inside rupture
//     area, 10000 outside as strength barrier — SCEC TPV5 §11).
//   - Per-QP `μ_d` and `d_c` are uniform (0.525 / 0.40 m).
//   - Per-QP `τ_strike(x, z)` via `ComputeTau2_0_TPV205` — three 3 km
//     square patches override the 70 MPa background:
//       nucleation  (0,    7.5 km)  →  81.6 MPa   (§7)
//       right       (+7.5, 7.5 km)  →  62.0 MPa   (§8)
//       left        (-7.5, 7.5 km)  →  78.0 MPa   (§9)
//   - Nucleation channels (`tau1_nuc`, `tau2_nuc`, `sigma_n_nuc`) are
//     ZERO at init AND are NEVER updated — TPV5 nucleation is via the
//     static patch pre-stress, not a time-varying perturbation.
//
// DOFData field repurposing for LSW (no header edits required):
//   `data.a`   ← μ_s (static friction; per-QP — barrier flag)
//   `data.psi` ← μ_d (dynamic friction; uniform)
//   `data.Dc`  ← d_c (slip-weakening critical distance; uniform)
//   `data.tau2_0` ← patch-dependent strike pre-stress
//   `data.slip1, data.slip2` ← accumulated slip; |δ| = sqrt(slip1²+slip2²)
//
// Coordinate frame ([C1]): pre-stress fields are written in the BP5 /
// Tandem canonical frame (`tangent1 = dip, tangent2 = strike`) so the
// strike-slip pre-stress goes into `tau2_0`, mirroring TPV102 / TPV104.

#ifndef MFEM_SEAS_TPV205_SETUP_HPP
#define MFEM_SEAS_TPV205_SETUP_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "fault_face_flux.hpp"
#include "tpv205_friction.hpp"
#include "../config/tpv205_params.hpp"

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

/// @brief Initialize DOFData array for TPV205 (linear slip-weakening).
///
/// Populates impedances, patch-dependent background pre-stress, per-QP
/// LSW friction parameters (μ_s with strength-barrier mask, μ_d, d_c),
/// and zeros the nucleation channels (TPV5 has no time-varying
/// nucleation perturbation).
///
/// Pre-stress convention (BP5 / Tandem canonical frame):
///   - `sigma_n0 = 120 MPa` (positive = compression)
///   - `tau1_0   = 0`       (no dip pre-stress — pure strike-slip at t=0)
///   - `tau2_0   = patch-dependent` (70 / 81.6 / 78.0 / 62.0 MPa)
///
/// LSW field repurposing:
///   - `d.a   = μ_s(x, z)`  (0.677 inside 30 km × 15 km rupture area,
///                           10000 outside — strength barrier per §11)
///   - `d.psi = μ_d`        (0.525, uniform)
///   - `d.Dc  = d_c`        (0.40 m, uniform)
///
/// @param[out] dof_data      Resized to `ndof`; each entry populated.
/// @param[in]  ndof          Number of fault DOFs (QPs).
/// @param[in]  fault_coords  Per-QP physical coordinates.
inline void InitializeFaultDOFs_TPV205(std::vector<DOFData> &dof_data,
                                       int ndof,
                                       const std::vector<Vector> &fault_coords)
{
   MFEM_VERIFY(ndof >= 0,
               "InitializeFaultDOFs_TPV205: ndof must be >= 0; got "
               << ndof);
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= ndof,
               "InitializeFaultDOFs_TPV205: fault_coords size "
               << fault_coords.size() << " < ndof " << ndof);

   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData &d = dof_data[i];

      // Impedances (homogeneous medium).
      d.Zp_plus  = TPV205Params::Zp;
      d.Zp_minus = TPV205Params::Zp;
      d.Zs_plus  = TPV205Params::Zs;
      d.Zs_minus = TPV205Params::Zs;
      d.eta_p    = TPV205Params::eta_p;
      d.eta_s    = TPV205Params::eta_s;

      // Fault coordinates: x = along-strike, |z| = down-dip depth.
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));

      // Patch-dependent background pre-stress (BP5 frame: strike → comp 2).
      d.sigma_n0 = TPV205Params::sigma_n;
      d.tau1_0   = 0.0;
      d.tau2_0   = ComputeTau2_0_TPV205(along_strike, down_dip);

      // Nucleation channels stay at 0 throughout the run — TPV5 has no
      // time-varying perturbation.  Zero them explicitly to defend
      // against a reused vector carrying stale bits.
      d.sigma_n_nuc = 0.0;
      d.tau1_nuc    = 0.0;
      d.tau2_nuc    = 0.0;

      // R-016: LSW-native fields.  Populates the dedicated slots so every
      // LSW reader (iterator, station writer, ParaView, the new
      // EvaluateADER_LSW) consumes them directly — no field repurposing,
      // no rate-and-state code path can accidentally read these as
      // a / psi / Dc.
      d.lsw_mu_s = ComputeMuS_TPV205(along_strike, down_dip);  // μ_s
      d.lsw_mu_d = ComputeMuD_TPV205(along_strike, down_dip);  // μ_d
      d.lsw_d_c  = ComputeDc_TPV205(along_strike, down_dip);   // d_c

      // Defensive: zero the rate-and-state slots.  TPV205 must NEVER
      // run rate-and-state.  If a future regression sends TPV205 data
      // into Evaluate / EvaluateTotal, the FrictionSolver guard sees
      // a == 0 and aborts loudly rather than computing nonsense
      // (REVIEW R-016 was caused by Brent silently running on the
      // repurposed values).
      d.a   = 0.0;
      d.psi = 0.0;
      d.Dc  = 0.0;

      // Initial slip rate — at rest.  TPV5 nucleation is spontaneous:
      // the static initial shear stress in the nucleation patch
      // (81.6 MPa) exceeds the static yield (0.677 × 120 = 81.24 MPa),
      // so the fault begins to slip immediately at t = 0+ inside that
      // patch.  No initial-velocity offset is needed.
      d.slip_rate = TPV205Params::V_ini;
      d.V1        = 0.0;
      d.V2        = TPV205Params::V_ini;
      d.slip1     = 0.0;
      d.slip2     = 0.0;

      // Initial corrected traction = background (no Riemann perturbation
      // yet; the wave-operator dispatch overwrites these at the first
      // ADER step).
      d.tau1_corr     = 0.0;
      d.tau2_corr     = d.tau2_0;
      d.sigma_n_corr  = TPV205Params::sigma_n;
   }
}

/// @brief Initialize the Q state vector for TPV205 (fluctuation-Q only).
///
/// Q = 0: the wave field represents perturbations from the background
/// equilibrium.  Background stress lives in DOFData side-channel
/// (`sigma_n0`, `tau1_0`, `tau2_0`); bulk Q carries fluctuation only.
inline void InitializeState_TPV205(Vector &Q, int ndof_total)
{
   MFEM_VERIFY(ndof_total >= 0,
               "InitializeState_TPV205: ndof_total must be >= 0; got "
               << ndof_total);
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
}

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
inline int FindNearestDOF_TPV205(const TPV205Station &station,
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
         const bool is_candidate = std::abs(local_dist[s]
                                            - global_min_dist[s]) < tie_tol;
         const int candidate_rank = is_candidate ? my_rank : nprocs_loc;
         int winning_rank;
         MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT,
                       MPI_MIN, comm);
         MFEM_VERIFY(winning_rank < nprocs_loc,
                     "TPV205StationWriter::Open: station "
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

/// @brief Surface station definition for TPV205 free-surface output.
///
/// R-010: the second coordinate is the FAULT-PERPENDICULAR distance
/// (mesh y axis), NOT depth.  In the project's coordinate system z is
/// depth and the y-axis is fault-normal; off-fault surface stations
/// live on z = 0 at non-zero y.  Field is named `fault_perp` to make
/// this explicit and avoid the foot-gun where a contributor might
/// supply a depth value.
struct TPV205SurfaceStation
{
   real_t x;            ///< Along-strike coordinate [m]
   real_t fault_perp;   ///< Fault-perpendicular coordinate (mesh y axis) [m]
   std::string name;
};

/// @brief Default TPV205 surface stations — INTENTIONALLY EMPTY.
///
/// R-006: SCEC TPV5 §III.2 references off-fault receivers in
/// `TPV5_forwebsite.pdf` Figure 2 that have not been transcribed yet.
/// Inheriting the TPV102 / TPV104 layout (six stations on the free
/// surface at ±9 / ±6 km fault-perpendicular distance) is unsafe — those
/// coordinates are not the SCEC TPV5 list and would silently produce
/// wrong-station traces under the default invocation.  Until the
/// canonical TPV5 surface-station list is added here, return an empty
/// vector so the SurfaceStationWriter writes nothing rather than
/// fabricating receivers at the wrong positions.  The driver emits a
/// banner warning when the list is empty.
inline std::vector<TPV205SurfaceStation> DefaultSurfaceStations_TPV205()
{
   return {};
}

/// @brief Surface station writer (mirrors TPV104SurfaceStationWriter).
class TPV205SurfaceStationWriter
{
public:
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV205SurfaceStation> &stations,
             Mesh &mesh, FiniteElementSpace &fes)
   {
      stations_ = stations;
      fes_ = &fes;
      ndof_per_el_ = fes.GetFE(0)->GetDof();
      ndof_total_ = mesh.GetNE() * ndof_per_el_;
      const int npts = static_cast<int>(stations.size());

      DenseMatrix point_mat(3, npts);
      for (int s = 0; s < npts; ++s)
      {
         point_mat(0, s) = stations[s].x;
         point_mat(1, s) = stations[s].fault_perp;
         point_mat(2, s) = 0.0;
      }

      elem_ids_.SetSize(npts);
      ip_refs_.SetSize(npts);
      mesh.FindPoints(point_mat, elem_ids_, ip_refs_, false);

      files_.resize(npts);
      for (int s = 0; s < npts; ++s)
      {
         if (elem_ids_[s] < 0) { continue; }
         const std::string fname = output_dir + "/" + prefix
                                   + "_station_" + stations[s].name
                                   + ".dat";
         files_[s].open(fname);
         if (files_[s].is_open())
         {
            files_[s] << "# TPV205 surface station: " << stations[s].name
                      << " (x=" << stations[s].x
                      << ", fault_perp=" << stations[s].fault_perp << ")\n";
            files_[s] << "# Columns: time vx vy vz\n";
         }
      }
   }

   void WriteStep(real_t t, const Vector &Q)
   {
      if (!fes_) { return; }

      for (int s = 0; s < static_cast<int>(stations_.size()); ++s)
      {
         const int e = elem_ids_[s];
         if (e < 0 || !files_[s].is_open()) { continue; }

         const FiniteElement *fe = fes_->GetFE(e);
         const int ndof = fe->GetDof();
         Vector shape(ndof);
         fe->CalcShape(ip_refs_[s], shape);

         const int dof_offset = e * ndof_per_el_;
         real_t vx = 0, vy = 0, vz = 0;
         for (int i = 0; i < ndof; ++i)
         {
            vx += shape(i) * Q[VX * ndof_total_ + dof_offset + i];
            vy += shape(i) * Q[VY * ndof_total_ + dof_offset + i];
            vz += shape(i) * Q[VZ * ndof_total_ + dof_offset + i];
         }

         files_[s] << std::scientific << std::setprecision(10)
                   << t << " " << vx << " " << vy << " " << vz << "\n";
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
   std::vector<TPV205SurfaceStation> stations_;
   FiniteElementSpace *fes_ = nullptr;
   int ndof_per_el_ = 0;
   int ndof_total_ = 0;
   Array<int> elem_ids_;
   Array<IntegrationPoint> ip_refs_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV205_SETUP_HPP
