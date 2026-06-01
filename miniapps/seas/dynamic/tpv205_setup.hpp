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
//
// NOTE: the on-fault station writer (TPV205Station / DefaultStations_TPV205
// / FindNearestDOF_TPV205 / TPV205StationWriter) was extracted to the lean
// `dynamic/tpv205_stations.hpp` (re-included below) so that header can be
// shared with `seas_spatial_dyn_driver` without dragging this file's
// fault-init / surface-station code into the spatial driver's translation
// unit.  Native-driver behaviour is unchanged — the symbols are still
// available here via the re-include.

#ifndef MFEM_SEAS_TPV205_SETUP_HPP
#define MFEM_SEAS_TPV205_SETUP_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "fault_face_flux.hpp"
#include "tpv205_friction.hpp"
#include "../config/tpv205_params.hpp"
#include "tpv205_stations.hpp"   // TPV205Station / DefaultStations_TPV205 / TPV205StationWriter (extracted)

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

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

// -------------------------------------------------------------------------
// On-fault station output (TPV205Station, DefaultStations_TPV205,
// FindNearestDOF_TPV205, TPV205StationWriter) lives in the lean
// `tpv205_stations.hpp`, re-included at the top of this file so the native
// driver continues to see these symbols unchanged.
// -------------------------------------------------------------------------

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
