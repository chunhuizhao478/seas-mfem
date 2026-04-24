// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 benchmark setup: fault-DOF init, pre-stress, nucleation-
// prestress channels, V_w side-channel, station output writers.
//
// R-006 (2026-04-24) Step 3 of TPV104 implementation.  Duplicates the
// `dynamic/tpv102_setup.hpp` layout (per §2.5 +
// feedback_tpv102_bp5_no_shared_edit: no shared helper between TPV102 and
// TPV104 — structurally independent for future quasi-dynamic + dynamic
// merges) with these TPV104-specific changes:
//
//   - Per-QP `a(x, z)` via `ComputeA_TPV104` (§4.1 TPV104 boxcar).
//   - Per-QP ψ_init via `ComputeInitialPsiTPV104(a_i)`  (§3.13 anchor
//     5.6359184e-01 for a_in = 0.01).
//   - `Dc = TPV104Params::L = 0.4 m`  (SCEC L, stored in DOFData::Dc).
//   - `V_ini = 1e-16 m/s`  (much lower than TPV102's 1e-12; SCEC TPV104).
//   - Pre-stress: `tau1_0 = 0`, `tau2_0 = 40 MPa`, `sigma_n0 = 120 MPa`
//     (§3.11 positive compression; §4.1 SCEC values).
//   - Nucleation channels zeroed at init (Step 6's accumulator writes
//     `tau2_nuc` incrementally per ADER sub-step).
//
// Fluctuation-Q-only (§3.10 directive): this header ships ONLY the
// mode-1 fluctuation-Q initialiser.  There is NO `InitializeStateTotal_TPV104`
// or `ApplyNucleationTotalPrestress_TPV104` — total-Q mode is dropped
// for TPV104 so bulk Q carries fluctuation only and pre-stress lives
// exclusively in the fault-CS side-channel (DOFData.sigma_n0, tau*_0).
// T_TPV104_SETUP_5 enforces this at link/grep time.
//
// Coordinate frame ([C1]): pre-stress fields are written in the BP5 /
// Tandem canonical frame (`tangent1 = dip, tangent2 = strike`), mirroring
// `dynamic/tpv102_setup.hpp:62-78` under R-801 Option A.  TPV104 is pure
// strike-slip, so `tau2_0 = tau_ini` and `V2 = V_ini`; `tau1_0` and `V1`
// remain zero.  T_TPV104_SETUP_4 pins this layout; T_TPV104_SETUP_3
// independently verifies the `FaultBasis::ComputeOrientedFrame` output
// for the TPV104 reference normal + up vectors so a future BP5 edit
// that flips the convention would trip this gate immediately ([C2]).

#ifndef MFEM_SEAS_TPV104_SETUP_HPP
#define MFEM_SEAS_TPV104_SETUP_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "fault_face_flux.hpp"
#include "../config/tpv104_params.hpp"

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

/// @brief Initialize DOFData array for TPV104.
///
/// Populates impedances, background stress, per-QP direct-effect
/// parameter, critical slip distance, initial state variable, and
/// initial corrected traction for every fault DOF.  Nucleation channels
/// (`sigma_n_nuc`, `tau1_nuc`, `tau2_nuc`) are zeroed — the Step-6
/// accumulator writes them incrementally at each ADER sub-step.
///
/// Pre-stress convention (BP5 / Tandem canonical frame, per R-801
/// Option A and §3.11 positive-compression decision):
///   - `sigma_n0 = 120 MPa` (positive compression)
///   - `tau1_0   = 0`       (no dip pre-stress — pure strike-slip)
///   - `tau2_0   = 40 MPa`  (along-strike pre-stress into tangent2)
///
/// @param[out] dof_data      Resized to `ndof`; each entry populated.
/// @param[in]  ndof          Number of fault DOFs (QPs).
/// @param[in]  fault_coords  Per-QP physical coordinates (x, y, z).
///                           `x = along-strike`, `|z| = down-dip`.
inline void InitializeFaultDOFs_TPV104(std::vector<DOFData> &dof_data,
                                       int ndof,
                                       const std::vector<Vector> &fault_coords)
{
   MFEM_VERIFY(ndof >= 0,
               "InitializeFaultDOFs_TPV104: ndof must be >= 0; got "
               << ndof);
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= ndof,
               "InitializeFaultDOFs_TPV104: fault_coords size "
               << fault_coords.size() << " < ndof " << ndof);

   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData &d = dof_data[i];

      // Impedances (homogeneous half-space — ρ, cp, cs from TPV104Params).
      d.Zp_plus  = TPV104Params::Zp;
      d.Zp_minus = TPV104Params::Zp;
      d.Zs_plus  = TPV104Params::Zs;
      d.Zs_minus = TPV104Params::Zs;
      d.eta_p    = TPV104Params::eta_p;
      d.eta_s    = TPV104Params::eta_s;

      // Background pre-stress in the BP5 / Tandem canonical frame
      // (tangent1 = dip, tangent2 = strike).  TPV104 is pure strike-slip
      // so `tau1_0 = 0` (no dip) and `tau2_0 = tau_ini` (along-strike).
      // `sigma_n0 > 0` = compression (geology convention — §3.11).
      d.sigma_n0 = TPV104Params::sigma_n;
      d.tau1_0   = 0.0;
      d.tau2_0   = TPV104Params::tau_ini;

      // Nucleation channels zeroed at init — the Step-6
      // `ApplyNucleationIncremental_TPV104` accumulator writes `tau2_nuc`
      // per sub-step.  `tau1_nuc` and `sigma_n_nuc` remain 0 throughout
      // (TPV104 is pure strike-slip; no normal-stress nucleation —
      // see R3-001 sentinel test in test_tpv104_nucleation.cpp).
      d.sigma_n_nuc = 0.0;
      d.tau1_nuc    = 0.0;
      d.tau2_nuc    = 0.0;

      // Fault coordinates: x = along-strike, |z| = down-dip depth.
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));

      // Per-QP direct-effect parameter (SCEC TPV104 boxcar a(x, z)).
      d.a  = ComputeA_TPV104(along_strike, down_dip);
      // SCEC TPV104 uses L = 0.4 m; DOFData stores it as `Dc`.
      d.Dc = TPV104Params::L;

      // Initial state variable from equilibrium inversion (§3.13 anchor).
      d.psi = ComputeInitialPsiTPV104(d.a);

      // Initial slip rate — pure strike-slip, into tangent2.
      d.slip_rate = TPV104Params::V_ini;
      d.V1        = 0.0;
      d.V2        = TPV104Params::V_ini;
      d.slip1     = 0.0;
      d.slip2     = 0.0;

      // Initial corrected traction = background (no Riemann perturbation
      // yet).  Same tangent convention as pre-stress.
      d.tau1_corr     = 0.0;
      d.tau2_corr     = TPV104Params::tau_ini;
      d.sigma_n_corr  = TPV104Params::sigma_n;
   }
}

/// @brief Initialize the Q state vector for TPV104 (fluctuation-Q only).
///
/// Q = 0: the wave field represents perturbations from the background
/// equilibrium.  Background stress (`tau_ini`, `sigma_n`) lives only in
/// the DOFData side-channel fields (`sigma_n0`, `tau1_0`, `tau2_0`) —
/// bulk Q carries fluctuation, matching the reference FVW runtime's
/// `initialStressInFaultCS` storage.  This prevents double-counting when
/// `sigma_n_total = sigma_n0 + sigma_n_trial` is evaluated (the trial
/// term comes from Q).
///
/// §3.10 directive: total-Q mode is NOT supported for TPV104.  There is
/// no `InitializeStateTotal_TPV104` companion.
///
/// @param[out] Q           Resized to `NUM_STATE * ndof_total`; set to 0.
/// @param[in]  ndof_total  Total DOFs per scalar component.
inline void InitializeState_TPV104(Vector &Q, int ndof_total)
{
   MFEM_VERIFY(ndof_total >= 0,
               "InitializeState_TPV104: ndof_total must be >= 0; got "
               << ndof_total);
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
}

/// @brief Populate the driver-owned per-QP V_w side-channel.
///
/// §4.5 directive: `V_w[i]` lives outside DOFData (the Extreme-Care
/// file `fault_face_flux.hpp` MUST NOT be edited).  The driver owns a
/// `std::vector<real_t> V_w` sized `num_fault_total` and populated once
/// at init.  Every `SlipLawSRWPsi::Rate_SRW` /
/// `UpdateStateAnalyticSlipLawSRW` call in the sub-step iterator reads
/// `V_w[i]` alongside `dof_data[i]`.
///
/// @param[out] V_w           Resized to `fault_coords.size()`.
/// @param[in]  fault_coords  Per-QP physical coordinates.
inline void PopulateVwSideChannel_TPV104(std::vector<real_t> &V_w,
                                         const std::vector<Vector> &fault_coords)
{
   const int n = static_cast<int>(fault_coords.size());
   V_w.resize(n);
   for (int i = 0; i < n; ++i)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      V_w[i] = ComputeVw_TPV104(along_strike, down_dip);
   }
}

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

/// Surface station definition for TPV104 free-surface output.
struct TPV104SurfaceStation
{
   real_t x;           ///< Along-strike coordinate [m]
   real_t z;           ///< Fault-normal coordinate [m]
   std::string name;
};

/// @brief Return the six SCEC TPV104 surface stations (same layout as
/// TPV102 — SCEC shares the same fault-normal layout across both).
inline std::vector<TPV104SurfaceStation> DefaultSurfaceStations_TPV104()
{
   return {
      { 0.0,    9.0e3, "surf_0_9"    },
      { 0.0,   -9.0e3, "surf_0_n9"   },
      { 12.0e3,  6.0e3, "surf_12_6"   },
      { 12.0e3, -6.0e3, "surf_12_n6"  },
      {-12.0e3,  6.0e3, "surf_n12_6"  },
      {-12.0e3, -6.0e3, "surf_n12_n6" },
   };
}

/// @brief Surface station writer for TPV104 free-surface data.
/// Uses `Mesh::FindPoints()` for element containment; works on both
/// `Mesh` and `ParMesh`.  Writes (t, vx, vy, vz) per station.
class TPV104SurfaceStationWriter
{
public:
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV104SurfaceStation> &stations,
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
         // Surface at Z = 0; X = along-strike, Y = fault-normal.
         point_mat(0, s) = stations[s].x;
         point_mat(1, s) = stations[s].z;
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
            files_[s] << "# TPV104 surface station: " << stations[s].name
                      << " (x=" << stations[s].x
                      << ", z=" << stations[s].z << ")\n";
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

   /// R8-003: explicit Close() symmetric to TPV104StationWriter::Close.
   /// The driver must call this before MPI_Finalize so that file
   /// handles are released BEFORE MPI shuts down (guarding the code
   /// path if surface output is ever routed through MPI-IO).  Without
   /// an explicit Close(), the destructor would run AFTER MPI_Finalize
   /// which is UB for MPI-IO handles.
   void Close()
   {
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   std::vector<TPV104SurfaceStation> stations_;
   FiniteElementSpace *fes_ = nullptr;
   int ndof_per_el_ = 0;
   int ndof_total_ = 0;
   Array<int> elem_ids_;
   Array<IntegrationPoint> ip_refs_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_SETUP_HPP
