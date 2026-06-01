// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 benchmark setup: spatial parameters, nucleation, station output.

#ifndef MFEM_SEAS_TPV102_SETUP_HPP
#define MFEM_SEAS_TPV102_SETUP_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "fault_face_flux.hpp"
#include "../config/tpv102_params.hpp"
#include "tpv102_stations.hpp"   // TPV102Station / DefaultStations / FindNearestDOF / TPV102StationWriter (extracted)

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Initialize DOFData array for TPV102 benchmark.
///
/// Sets impedances, background stress, friction parameters (spatially varying a),
/// and initial state variable psi from equilibrium for each fault DOF.
///
/// @param[out] dof_data  Array of DOFData to initialize (one per fault DOF).
/// @param[in] ndof  Number of fault DOFs.
/// @param[in] fault_coords  Fault DOF coordinates: fault_coords[i] = (x, y, z).
///   x = along-strike, z = depth (z <= 0, surface at z = 0).
///   down_dip = |z| (positive distance from surface).
inline void InitializeFaultDOFs(std::vector<DOFData> &dof_data, int ndof,
                                const std::vector<Vector> &fault_coords)
{
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= ndof,
               "fault_coords size " << fault_coords.size() << " < ndof " << ndof);

   dof_data.resize(ndof);

   for (int i = 0; i < ndof; i++)
   {
      DOFData &d = dof_data[i];

      // Impedances (homogeneous)
      d.Zp_plus  = TPV102Params::Zp;
      d.Zp_minus = TPV102Params::Zp;
      d.Zs_plus  = TPV102Params::Zs;
      d.Zs_minus = TPV102Params::Zs;
      d.eta_p    = TPV102Params::Zp / 2.0;
      d.eta_s    = TPV102Params::eta_s;

      // Background stress (uniform, in the BP5 canonical fault-local frame).
      // Under R-801 Option A the whole TPV102 pipeline uses BP5's
      // FaultBasis convention — tangent1 = dip, tangent2 = strike
      // (Tandem convention, see fault/fault_basis.hpp:54-56).  For TPV102's
      // vertical planar fault at y=0 with ref_normal=(0,-1,0) and up=(0,0,1):
      //   can_t1 = dip    = (0, 0, -1)   (-z = down into earth)
      //   can_t2 = strike = (+1, 0, 0)   (+x = along strike)
      // TPV102 is pure strike-slip, so the along-strike pre-stress and
      // initial slip rate live in COMPONENT 2 (tangent2), not component 1.
      // Pre-R-801, this file wrote strike into component 1 under the
      // GodunovFlux::BuildFrame convention (t1=x=strike); that made the
      // interior-fault path self-consistent but collided with the shared-
      // fault path (which always used BP5's canonical frame via R-701),
      // producing mixed semantics for DOFData.V1/V2/tau1_corr/tau2_corr
      // across QPs on the same fault.  Option A puts every QP on BP5's
      // convention, restoring a single source of truth.
      d.sigma_n0 = TPV102Params::sigma_n;
      d.tau1_0   = 0.0;                      // no dip pre-stress
      d.tau2_0   = TPV102Params::tau_ini;    // along-strike pre-stress

      // v9.4.0 Commit 1 / R-003: persistent nucleation channel.
      // Evaluate() now reads data.tau*_nuc / sigma_n_nuc into the
      // friction input.  DOFData default-ctor already zeroes these
      // fields, but zeroing explicitly here defends against a reused
      // vector carrying stale bits from a prior run.  ApplyNucleation
      // Prestress overwrites tau2_nuc per call; tau1_nuc / sigma_n_nuc
      // stay zero for TPV102 (pure strike-slip, no normal-stress nuc).
      d.sigma_n_nuc = 0.0;
      d.tau1_nuc    = 0.0;
      d.tau2_nuc    = 0.0;

      // Fault coordinates: x = along-strike, z = depth
      real_t along_strike = fault_coords[i](0);
      real_t down_dip = std::abs(fault_coords[i](2));  // depth as positive distance

      // Spatially varying direct effect parameter
      d.a  = ComputeA(along_strike, down_dip);
      // Phase 11a: the substep iterator now reads the state-evolution
      // parameter per-DOF from d.b (was the scalar state_evo_.GetB()).  Set it
      // to TPV102's scalar b so native TPV102 is byte-identical to the
      // pre-Phase-11 path (AgingLawPsi in tpv102_driver.cpp uses TPV102Params::b).
      d.b  = TPV102Params::b;
      d.Dc = TPV102Params::Dc;

      // Initial state from equilibrium
      d.psi = ComputeInitialPsi(d.a);

      // Initially locked — TPV102 initial slip rate is purely along-strike
      // (mode II).  Under BP5 convention strike lives in component 2.
      d.slip_rate = TPV102Params::V_ini;
      d.V1 = 0.0;                            // no dip slip rate
      d.V2 = TPV102Params::V_ini;            // along-strike initial slip rate
      d.slip1 = 0.0;
      d.slip2 = 0.0;

      // Initial corrected traction = background (no perturbation).  Same
      // component assignment as the pre-stress: strike in component 2.
      d.tau1_corr = 0.0;
      d.tau2_corr = TPV102Params::tau_ini;
      d.sigma_n_corr = TPV102Params::sigma_n;
   }
}

/// @brief Initialize the Q state vector for TPV102 equilibrium.
///
/// Q = 0: the wave field represents perturbations from the background
/// equilibrium state. Background stress (tau_ini, sigma_n) lives only
/// in DOFData pre-stress fields (sigma_n0, tau1_0, tau2_0), not in Q.
/// This prevents double-counting when the fault coupling computes
/// sigma_n_total = sigma_n0 + sigma_n_trial (where sigma_n_trial comes from Q).
///
/// @param[out] Q  State vector (9 * ndof_total), set to zero.
/// @param[in] ndof_total  Total DOFs per scalar component.
inline void InitializeState(Vector &Q, int ndof_total)
{
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
}

/// @brief Apply nucleation perturbation to fault DOF shear stress.
///
/// Adds delta_tau(x, z, t) to the background along-strike shear pre-stress
/// for each DOF.  Under R-801 Option A the BP5 canonical frame is in force
/// everywhere, so the along-strike pre-stress lives in `tau2_0` (tangent2 =
/// strike), NOT `tau1_0` (tangent1 = dip).  This modifies the pre-stress
/// term that enters the trial traction computation.
///
/// @param[in,out] dof_data  Fault DOF data array.
/// @param[in] ndof  Number of fault DOFs.
/// @param[in] fault_coords  Fault DOF physical coordinates.
/// @param[in] t  Current simulation time [s].
inline void ApplyNucleation(std::vector<DOFData> &dof_data, int ndof,
                            const std::vector<Vector> &fault_coords,
                            real_t t)
{
   for (int i = 0; i < ndof; i++)
   {
      real_t along_strike = fault_coords[i](0);
      real_t down_dip = std::abs(fault_coords[i](2));

      real_t dtau = NucleationPerturbation(along_strike, down_dip, t);
      // BP5 convention (tangent2 = strike): along-strike pre-stress = tau2_0.
      dof_data[i].tau2_0 = TPV102Params::tau_ini + dtau;
   }
}

// -------------------------------------------------------------------------
// On-fault station output (TPV102Station, DefaultStations, FindNearestDOF,
// TPV102StationWriter) lives in the lean `tpv102_stations.hpp`, re-included
// at the top of this file so the native driver continues to see these
// symbols unchanged.
// -------------------------------------------------------------------------

/// Surface station definition for TPV102 free-surface output.
struct TPV102SurfaceStation
{
   real_t x;       ///< Along-strike coordinate [m]
   real_t z;       ///< Fault-normal coordinate [m] (z<0 = −side, z>0 = +side)
   std::string name;
};

/// @brief Return the 6 SCEC-required surface stations for TPV102.
///
/// Located on the free surface (y=0) at specified (x, z) coordinates.
/// These require interpolating the wave state Q, not fault DOF lookup.
inline std::vector<TPV102SurfaceStation> DefaultSurfaceStations()
{
   return {
      // (x, z) in km → meters. SCEC spec: z=y (fault-normal at surface)
      {0.0,    9e3,  "surf_0_9"},
      {0.0,   -9e3,  "surf_0_n9"},
      {12e3,   6e3,  "surf_12_6"},
      {12e3,  -6e3,  "surf_12_n6"},
      {-12e3,  6e3,  "surf_n12_6"},
      {-12e3, -6e3,  "surf_n12_n6"},
   };
}

/// @brief Surface station writer for TPV102 free-surface data.
///
/// Uses Mesh::FindPoints() for proper element containment and reference
/// coordinate computation. Works with both Mesh and ParMesh (polymorphic).
/// Writes (time, vx, vy, vz) per station.
class TPV102SurfaceStationWriter
{
public:
   /// Open surface station output files and locate points in the mesh.
   void Open(const std::string &output_dir,
             const std::string &prefix,
             const std::vector<TPV102SurfaceStation> &stations,
             Mesh &mesh, FiniteElementSpace &fes)
   {
      stations_ = stations;
      fes_ = &fes;
      ndof_per_el_ = fes.GetFE(0)->GetDof();
      ndof_total_ = mesh.GetNE() * ndof_per_el_;
      int npts = static_cast<int>(stations.size());

      // Build point matrix (dim x npts) for FindPoints
      DenseMatrix point_mat(3, npts);
      for (int s = 0; s < npts; s++)
      {
         // Surface at Z=0; X=along-strike, Y=fault-normal
         point_mat(0, s) = stations[s].x;
         point_mat(1, s) = stations[s].z;
         point_mat(2, s) = 0.0;
      }

      // FindPoints: returns elem_ids and reference IntegrationPoints.
      // For ParMesh, automatically handles MPI ownership.
      // elem_ids[s] = -1 if point not found on this rank.
      elem_ids_.SetSize(npts);
      ip_refs_.SetSize(npts);
      mesh.FindPoints(point_mat, elem_ids_, ip_refs_, false);

      // Open files only for stations found on this rank
      files_.resize(npts);
      for (int s = 0; s < npts; s++)
      {
         if (elem_ids_[s] < 0) { continue; }
         std::string fname = output_dir + "/" + prefix + "_station_"
                           + stations[s].name + ".dat";
         files_[s].open(fname);
         if (files_[s].is_open())
         {
            files_[s] << "# TPV102 surface station: " << stations[s].name
                      << " (x=" << stations[s].x
                      << ", z=" << stations[s].z << ")\n";
            files_[s] << "# Columns: time vx vy vz\n";
         }
      }
   }

   /// Write velocity at surface stations from the wave state Q.
   void WriteStep(real_t t, const Vector &Q)
   {
      if (!fes_) { return; }

      for (int s = 0; s < static_cast<int>(stations_.size()); s++)
      {
         int e = elem_ids_[s];
         if (e < 0 || !files_[s].is_open()) { continue; }

         const FiniteElement *fe = fes_->GetFE(e);
         int ndof = fe->GetDof();
         Vector shape(ndof);
         fe->CalcShape(ip_refs_[s], shape);

         int dof_offset = e * ndof_per_el_;
         real_t vx = 0, vy = 0, vz = 0;
         for (int i = 0; i < ndof; i++)
         {
            vx += shape(i) * Q[VX * ndof_total_ + dof_offset + i];
            vy += shape(i) * Q[VY * ndof_total_ + dof_offset + i];
            vz += shape(i) * Q[VZ * ndof_total_ + dof_offset + i];
         }

         files_[s] << std::scientific << std::setprecision(10)
                   << t << " " << vx << " " << vy << " " << vz << "\n";
         // F1 parity: flush every write so a Slurm SIGKILL at the wall
         // limit does not discard buffered probe data.
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
   std::vector<TPV102SurfaceStation> stations_;
   FiniteElementSpace *fes_ = nullptr;
   int ndof_per_el_ = 0;
   int ndof_total_ = 0;
   Array<int> elem_ids_;
   Array<IntegrationPoint> ip_refs_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_SETUP_HPP
