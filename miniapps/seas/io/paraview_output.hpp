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

#ifndef MFEM_SEAS_PARAVIEW_OUTPUT_HPP
#define MFEM_SEAS_PARAVIEW_OUTPUT_HPP

#include "mfem.hpp"
#include "../config/bp5_params.hpp"

#include <string>
#include <cmath>
#include <memory>
#include <sys/stat.h>
#include <vector>
#include <type_traits>

namespace mfem
{
namespace seas
{

// Helper to select GridFunction vs ParGridFunction based on MeshType.
template <typename MeshType> struct GFType
{ using type = GridFunction; using FESType = FiniteElementSpace; };

#ifdef MFEM_USE_MPI
template <> struct GFType<ParMesh>
{ using type = ParGridFunction; using FESType = ParFiniteElementSpace; };
#endif

/// @brief Manages ParaView-compatible visualization output for SEAS simulations.
///
/// Outputs domain fields (displacement) directly and fault fields (slip,
/// slip rate, traction, state variable) as piecewise-constant (L2 p=0)
/// GridFunctions on the domain mesh.  Fault values are mapped from owned
/// fault DOF vectors to the two volume elements adjacent to each fault face.
///
/// Output frequency is adaptive, keyed to V_max (slip rate).  Defaults:
///   - coseismic  (V > 1e-3 m/s)  : every 0.01 s
///   - nucleation (V > 1e-6)       : every 1.0 s
///   - interseismic                : every 1 year
/// Thresholds and intervals are runtime-configurable via GetSchedule()
/// (see AdaptiveSchedule).  A hysteresis_factor > 1 suppresses regime
/// chattering near the thresholds.
///
/// @tparam MeshType  Mesh (serial) or ParMesh (parallel)
template <typename MeshType = Mesh>
class ParaViewOutput
{
   using GF  = typename GFType<MeshType>::type;
   using FES = typename GFType<MeshType>::FESType;

   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   /// Step-based output interval (write every N steps).  Set to 0 to
   /// use the adaptive V_max or fixed-dt schedule.
   int output_every_n_steps = 0;

   /// Fixed time interval between writes (seconds).  Set to 0 to use
   /// the adaptive V_max schedule.  Takes precedence over V_max schedule
   /// but not over step-based interval.
   real_t fixed_dt = 0.0;

   /// Runtime-configurable adaptive-schedule parameters.  Defaults match
   /// the legacy hardcoded OutputInterval thresholds/intervals and have
   /// hysteresis_factor = 1.0 (no hysteresis), so existing callers see
   /// bit-identical behavior.
   ///
   /// Regime indices: 0 = interseismic, 1 = nucleation, 2 = coseismic.
   struct AdaptiveSchedule
   {
      // V thresholds (m/s).  Coseismic enter when V exceeds v_coseismic;
      // nucleation enter when V exceeds v_nucleation.
      real_t v_coseismic       = 1e-3;
      real_t v_nucleation      = 1e-6;
      // Hysteresis factor (>= 1).  Leaving a regime requires V to drop
      // below (threshold / hysteresis_factor).  1.0 disables hysteresis.
      real_t hysteresis_factor = 1.0;
      // Output intervals (seconds) used in each regime.
      real_t dt_coseismic      = 0.01;
      real_t dt_nucleation     = 1.0;
      real_t dt_interseismic   = 1.0 * BP5Params::seconds_per_year;

      /// Return the output interval (seconds) for a given regime.
      /// V is accepted for signature symmetry with NextRegime but is
      /// unused — the cadence is stable within a regime.
      real_t Interval(real_t /*V*/, int regime) const
      {
         switch (regime)
         {
            case 2:  return dt_coseismic;
            case 1:  return dt_nucleation;
            default: return dt_interseismic;
         }
      }

      /// State machine mapping (V, prev_regime) -> new regime.
      ///
      /// Entry comparisons are STRICT > (not >=) to preserve byte
      /// compatibility with the legacy OutputInterval which used
      /// V_max > 1e-3 / V_max > 1e-6.  Exit comparisons are strict <,
      /// so an exact-threshold V stays in its current regime.
      ///
      /// Non-finite V (NaN/Inf) is treated as coseismic, so we keep
      /// dense output during a blowup.
      int NextRegime(real_t V, int prev) const
      {
         const real_t V_co_enter = v_coseismic;
         const real_t V_co_exit  = v_coseismic  / hysteresis_factor;
         const real_t V_nu_enter = v_nucleation;
         const real_t V_nu_exit  = v_nucleation / hysteresis_factor;

         if (!std::isfinite(V)) { return 2; }

         switch (prev)
         {
            case 0: // interseismic
               if (V > V_co_enter) { return 2; }
               if (V > V_nu_enter) { return 1; }
               return 0;
            case 1: // nucleation
               if (V > V_co_enter) { return 2; }
               if (V < V_nu_exit)  { return 0; }
               return 1;
            case 2: // coseismic
               if (V < V_co_exit)
               {
                  // Drop-past-threshold semantics (intentional):
                  // the target uses the STRICT ENTRY threshold, so a V
                  // landing in [V_nu_exit, V_nu_enter] goes straight to
                  // regime 0 (interseismic) — we do NOT transfer any
                  // accumulated nucleation-hysteresis credit from the
                  // outbound path.  Rationale: after a completed
                  // rupture, V is usually already deep in interseismic
                  // territory; returning to nucleation only when V is
                  // unambiguously above the entry threshold avoids
                  // over-sampling post-rupture ring-down.
                  if (V > V_nu_enter) { return 1; }
                  return 0;
               }
               return 2;
            default: return 0;
         }
      }

      /// One-shot invariant check.  Call once after applying any CLI
      /// overrides (NOT inside NextRegime / Interval — those are
      /// hot-path).  Aborts with a readable message on bad inputs.
      void Validate() const
      {
         MFEM_VERIFY(hysteresis_factor >= 1.0,
                     "AdaptiveSchedule: hysteresis_factor must be >= 1.0, got "
                     << hysteresis_factor);
         MFEM_VERIFY(v_coseismic > v_nucleation && v_nucleation > 0.0,
                     "AdaptiveSchedule: require v_coseismic (" << v_coseismic
                     << ") > v_nucleation (" << v_nucleation << ") > 0");
         MFEM_VERIFY(dt_coseismic > 0.0 && dt_nucleation > 0.0 &&
                     dt_interseismic > 0.0,
                     "AdaptiveSchedule: all dt_* must be positive");
      }
   };

   /// Mutable accessor so drivers can override thresholds/intervals at
   /// configuration time.  Call Validate() after tweaking.
   AdaptiveSchedule &GetSchedule() { return adaptive_; }
   const AdaptiveSchedule &GetSchedule() const { return adaptive_; }

   /// @brief Construct the ParaView output manager.
   ParaViewOutput(const std::string &prefix,
                  MeshType &mesh,
                  int order)
      : mesh_(mesh),
        order_(order),
        last_write_time_(-1e30),
        pv_(prefix, &mesh)
   {
      pv_.SetPrefixPath(prefix);
      pv_.SetDataFormat(VTKFormat::BINARY);
      pv_.SetHighOrderOutput(true);
      pv_.SetLevelsOfDetail(order);
   }

   // ---------------------------------------------------------------
   //  Domain field registration
   // ---------------------------------------------------------------

   /// Register a domain-mesh GridFunction for output (non-owning).
   void RegisterDomainField(const std::string &name, GF *gf)
   {
      pv_.RegisterField(name, gf);
   }

   // ---------------------------------------------------------------
   //  BP5 fault field output (2-component tangential + 1 state)
   // ---------------------------------------------------------------

   /// @brief Initialize L2-p0 fault fields on the domain mesh for BP5.
   ///
   /// Creates 7 scalar GridFunctions (slip_dip, slip_strike, slip_rate_dip,
   /// slip_rate_strike, traction_dip, traction_strike, state_variable) as
   /// piecewise-constant on the volume mesh.  Only elements adjacent to
   /// fault faces carry nonzero values.
   ///
   /// @param fault_interior_faces  Interior fault face indices (local faces)
   /// @param fault_shared_faces    Shared fault face indices (for ParMesh)
   /// @param nbf_per_face          Basis functions per fault face (e.g. 3 for p=1 tri)
   void InitFaultOutputBP5(const Array<int> &fault_interior_faces,
                           const Array<int> &fault_shared_faces,
                           int nbf_per_face)
   {
      nbf_per_face_ = nbf_per_face;
      int dim = mesh_.Dimension();

      // L2 p=0 on domain mesh — one DOF per element
      fault_fec_ = std::make_unique<L2_FECollection>(0, dim);
      fault_fes_ = MakeFES(mesh_, fault_fec_.get());

      // Allocate 7 scalar fields
      auto make_gf = [&]() {
         auto gf = std::make_unique<GF>(fault_fes_.get());
         *gf = 0.0;
         return gf;
      };
      fault_slip_dip_        = make_gf();
      fault_slip_strike_     = make_gf();
      fault_slip_rate_dip_   = make_gf();
      fault_slip_rate_strike_= make_gf();
      fault_trac_dip_        = make_gf();
      fault_trac_strike_     = make_gf();
      fault_state_           = make_gf();
      fault_normal_stress_   = make_gf();

      // Build fault face → (elem1, elem2) mapping for interior faces.
      // The driver list may include orphan / shared-as-bdr faces (kept to
      // preserve fi*nbf index alignment with fault_coords); for those the
      // interior-transformation lookup returns nullptr and we record -1 so
      // Set/UpdateFaultFieldsBP5 can skip the L2-p0 scatter.
      int n_int = fault_interior_faces.Size();
      fault_face_elem1_.resize(n_int);
      fault_face_elem2_.resize(n_int);
      for (int i = 0; i < n_int; i++)
      {
         int face = fault_interior_faces[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr)
         {
            fault_face_elem1_[i] = FTr->Elem1No;
            fault_face_elem2_[i] = FTr->Elem2No;
         }
         else
         {
            fault_face_elem1_[i] = -1;
            fault_face_elem2_[i] = -1;
         }
      }
      n_interior_fault_faces_ = n_int;

      // Shared faces: only Elem1 is local (Elem2 is on neighbor rank).
      // Map to local Elem1 only.
      n_shared_fault_faces_ = 0;
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         n_shared_fault_faces_ = fault_shared_faces.Size();
         fault_shared_elem1_.resize(n_shared_fault_faces_);
         for (int i = 0; i < n_shared_fault_faces_; i++)
         {
            int sf = fault_shared_faces[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            fault_shared_elem1_[i] = FTr->Elem1No;
         }
#endif
      }

      // Friction parameter / coordinate fields (static)
      fault_param_a_   = make_gf();
      fault_param_Dc_  = make_gf();
      fault_coord_x2_  = make_gf();
      fault_coord_x3_  = make_gf();

      // Register fields
      pv_.RegisterField("slip_dip",         fault_slip_dip_.get());
      pv_.RegisterField("slip_strike",      fault_slip_strike_.get());
      pv_.RegisterField("slip_rate_dip",    fault_slip_rate_dip_.get());
      pv_.RegisterField("slip_rate_strike", fault_slip_rate_strike_.get());
      pv_.RegisterField("traction_dip",     fault_trac_dip_.get());
      pv_.RegisterField("traction_strike",  fault_trac_strike_.get());
      pv_.RegisterField("state_variable",   fault_state_.get());
      pv_.RegisterField("normal_stress",    fault_normal_stress_.get());
      pv_.RegisterField("param_a",          fault_param_a_.get());
      pv_.RegisterField("param_Dc",         fault_param_Dc_.get());
      pv_.RegisterField("fault_x2",         fault_coord_x2_.get());
      pv_.RegisterField("fault_x3",         fault_coord_x3_.get());

      fault_interior_faces_ = fault_interior_faces;
      fault_shared_faces_ = fault_shared_faces;

      has_fault_output_ = true;
   }

   /// @brief Set static friction parameters and coordinates on fault elements.
   ///
   /// Call once after InitFaultOutputBP5. Vectors are in local (all faces)
   /// layout with 1 component per DOF.
   void SetFaultParamsBP5(const Vector &local_a,
                          const Vector &local_Dc,
                          const Vector &local_x2,
                          const Vector &local_x3)
   {
      if (!has_fault_output_) { return; }

      *fault_param_a_  = 0.0;
      *fault_param_Dc_ = 0.0;
      *fault_coord_x2_ = 0.0;
      *fault_coord_x3_ = 0.0;

      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;

      for (int fi = 0; fi < n_int; fi++)
      {
         int base = fi * nbf;
         real_t avg_a = 0, avg_Dc = 0, avg_x2 = 0, avg_x3 = 0;
         for (int k = 0; k < nbf; k++)
         {
            avg_a  += local_a(base + k);
            avg_Dc += local_Dc(base + k);
            avg_x2 += local_x2(base + k);
            avg_x3 += local_x3(base + k);
         }
         real_t inv = 1.0 / nbf;
         avg_a *= inv; avg_Dc *= inv; avg_x2 *= inv; avg_x3 *= inv;

         int e1 = fault_face_elem1_[fi];
         int e2 = fault_face_elem2_[fi];
         // Skip orphan / shared-as-bdr entries kept only for index alignment.
         if (e1 < 0) { continue; }
         (*fault_param_a_)(e1)  = avg_a;
         (*fault_param_Dc_)(e1) = avg_Dc;
         (*fault_coord_x2_)(e1) = avg_x2;
         (*fault_coord_x3_)(e1) = avg_x3;
         if (e2 >= 0)
         {
            (*fault_param_a_)(e2)  = avg_a;
            (*fault_param_Dc_)(e2) = avg_Dc;
            (*fault_coord_x2_)(e2) = avg_x2;
            (*fault_coord_x3_)(e2) = avg_x3;
         }
      }

      for (int si = 0; si < n_shared_fault_faces_; si++)
      {
         int fi = n_int + si;
         int base = fi * nbf;
         real_t avg_a = 0, avg_Dc = 0, avg_x2 = 0, avg_x3 = 0;
         for (int k = 0; k < nbf; k++)
         {
            avg_a  += local_a(base + k);
            avg_Dc += local_Dc(base + k);
            avg_x2 += local_x2(base + k);
            avg_x3 += local_x3(base + k);
         }
         real_t inv = 1.0 / nbf;
         avg_a *= inv; avg_Dc *= inv; avg_x2 *= inv; avg_x3 *= inv;

         int e1 = fault_shared_elem1_[si];
         (*fault_param_a_)(e1)  = avg_a;
         (*fault_param_Dc_)(e1) = avg_Dc;
         (*fault_coord_x2_)(e1) = avg_x2;
         (*fault_coord_x3_)(e1) = avg_x3;
      }
   }

   /// @brief Update all fault L2-p0 fields from owned-DOF vectors.
   ///
   /// Maps fault DOF values to the adjacent volume elements.  Each fault
   /// face has `nbf_per_face` DOFs; we average them for the element value.
   /// Vectors are "local" (all fault faces on this rank, interior+shared),
   /// with 2-component layout: [comp0_dof0, comp0_dof1, ..., comp1_dof0, ...]
   ///
   /// @param local_slip       Slip vector (2 * num_local_fault_dofs)
   /// @param local_slip_rate  Slip rate vector (2 * num_local_fault_dofs)
   /// @param local_traction   Traction vector (2 * num_local_fault_dofs)
   /// @param local_state      State (psi) vector (1 * num_local_fault_dofs)
   /// @param local_normal_stress  Normal stress vector (1 * num_local_fault_dofs),
   ///                             may be empty if elastic sigma_n is disabled
   void UpdateFaultFieldsBP5(const Vector &local_slip,
                             const Vector &local_slip_rate,
                             const Vector &local_traction,
                             const Vector &local_state,
                             const Vector &local_normal_stress = Vector())
   {
      if (!has_fault_output_) { return; }

      // Zero all fields (most elements won't have fault data)
      *fault_slip_dip_ = 0.0;        *fault_slip_strike_ = 0.0;
      *fault_slip_rate_dip_ = 0.0;   *fault_slip_rate_strike_ = 0.0;
      *fault_trac_dip_ = 0.0;        *fault_trac_strike_ = 0.0;
      *fault_state_ = 0.0;
      *fault_normal_stress_ = 0.0;
      const bool has_normal = (local_normal_stress.Size() > 0);

      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;

      // Interior faces: assign to both elem1 and elem2
      for (int fi = 0; fi < n_int; fi++)
      {
         int base = fi * nbf;  // local fault DOF offset
         real_t s_d = 0, s_s = 0, sr_d = 0, sr_s = 0;
         real_t tr_d = 0, tr_s = 0, psi = 0, sn = 0;

         for (int k = 0; k < nbf; k++)
         {
            int dof = base + k;
            s_d  += local_slip(2 * dof);
            s_s  += local_slip(2 * dof + 1);
            sr_d += local_slip_rate(2 * dof);
            sr_s += local_slip_rate(2 * dof + 1);
            tr_d += local_traction(2 * dof);
            tr_s += local_traction(2 * dof + 1);
            psi  += local_state(dof);
            if (has_normal) { sn += local_normal_stress(dof); }
         }
         real_t inv = 1.0 / nbf;
         s_d *= inv; s_s *= inv; sr_d *= inv; sr_s *= inv;
         tr_d *= inv; tr_s *= inv; psi *= inv; sn *= inv;

         int e1 = fault_face_elem1_[fi];
         int e2 = fault_face_elem2_[fi];
         // Skip orphan / shared-as-bdr entries kept only for index alignment.
         if (e1 < 0) { continue; }
         (*fault_slip_dip_)(e1)         = s_d;
         (*fault_slip_strike_)(e1)      = s_s;
         (*fault_slip_rate_dip_)(e1)    = sr_d;
         (*fault_slip_rate_strike_)(e1) = sr_s;
         (*fault_trac_dip_)(e1)         = tr_d;
         (*fault_trac_strike_)(e1)      = tr_s;
         (*fault_state_)(e1)            = psi;
         if (has_normal) { (*fault_normal_stress_)(e1) = sn; }
         if (e2 >= 0)
         {
            (*fault_slip_dip_)(e2)         = s_d;
            (*fault_slip_strike_)(e2)      = s_s;
            (*fault_slip_rate_dip_)(e2)    = sr_d;
            (*fault_slip_rate_strike_)(e2) = sr_s;
            (*fault_trac_dip_)(e2)         = tr_d;
            (*fault_trac_strike_)(e2)      = tr_s;
            (*fault_state_)(e2)            = psi;
            if (has_normal) { (*fault_normal_stress_)(e2) = sn; }
         }
      }

      // Shared faces: only elem1 is local
      for (int si = 0; si < n_shared_fault_faces_; si++)
      {
         int fi = n_int + si;  // local fault face index (after interior)
         int base = fi * nbf;
         real_t s_d = 0, s_s = 0, sr_d = 0, sr_s = 0;
         real_t tr_d = 0, tr_s = 0, psi = 0, sn = 0;

         for (int k = 0; k < nbf; k++)
         {
            int dof = base + k;
            s_d  += local_slip(2 * dof);
            s_s  += local_slip(2 * dof + 1);
            sr_d += local_slip_rate(2 * dof);
            sr_s += local_slip_rate(2 * dof + 1);
            tr_d += local_traction(2 * dof);
            tr_s += local_traction(2 * dof + 1);
            psi  += local_state(dof);
            if (has_normal) { sn += local_normal_stress(dof); }
         }
         real_t inv = 1.0 / nbf;
         s_d *= inv; s_s *= inv; sr_d *= inv; sr_s *= inv;
         tr_d *= inv; tr_s *= inv; psi *= inv; sn *= inv;

         int e1 = fault_shared_elem1_[si];
         (*fault_slip_dip_)(e1)         = s_d;
         (*fault_slip_strike_)(e1)      = s_s;
         (*fault_slip_rate_dip_)(e1)    = sr_d;
         (*fault_slip_rate_strike_)(e1) = sr_s;
         (*fault_trac_dip_)(e1)         = tr_d;
         (*fault_trac_strike_)(e1)      = tr_s;
         (*fault_state_)(e1)            = psi;
         if (has_normal) { (*fault_normal_stress_)(e1) = sn; }
      }
   }

   // ---------------------------------------------------------------
   //  Fault surface VTU output (proper 2D face geometry)
   // ---------------------------------------------------------------

   /// @brief Write fault surface as a triangle VTU, one averaged value per cell.
   ///
   /// Each rank writes fault_surface_r{rank}_c{cycle}.vtu containing the
   /// fault-face triangles and one CellData value per triangle.  Each cell
   /// value is the arithmetic mean of the DG field evaluated at the face's
   /// `nbf_per_face_` quadrature/nodal points.  Averaging (rather than
   /// per-vertex interpolation from QP-at-vertex coordinates) avoids the
   /// R-001 / H-V91-A4 speckle artefact that arose when QP values were
   /// written to reference-triangle vertex positions directly (plan
   /// v9.1.0 §2.3).  It also side-steps the BR2 latent bug R-005 where
   /// the old hardcoded `k<3` loop crossed face boundaries for `nbf=1`.
   /// Rank 0 also writes a .pvtu index.
   ///
   /// @param prefix   Output directory path
   /// @param cycle    Time step number
   /// @param time     Simulation time
   /// @param rank     MPI rank
   /// @param nranks   Total number of MPI ranks
   void WriteFaultSurfaceVTU(
      const std::string &prefix,
      int cycle, real_t time, int rank, int nranks,
      const Vector &local_slip,
      const Vector &local_slip_rate,
      const Vector &local_traction,
      const Vector &local_state,
      const Vector &local_normal_stress,
      const Vector &local_a,
      const Vector &local_Dc,
      const Vector &local_x2,
      const Vector &local_x3)
   {
      if (!has_fault_output_) { return; }
      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;
      const int n_shared = n_shared_fault_faces_;
      const int n_faces = n_int + n_shared;
      const bool has_normal = (local_normal_stress.Size() > 0);

      // R-007 (v9.1.0 rev 3): hard-fail on nbf<=0 in both Debug and Release.
      // The per-face average below divides by nbf; nbf=0 would produce
      // inv_nbf = 1/0 = +Inf and 0*Inf = NaN in every CellData entry.
      // InitFaultOutputBP5 assigns nbf_per_face_ = caller-supplied value
      // with no sanity check, so guard here at the call boundary.
      MFEM_VERIFY(nbf > 0,
                  "WriteFaultSurfaceVTU: nbf_per_face_ must be > 0, got "
                  << nbf);

      // Collect face vertices (geometry only) and per-CELL (one entry per
      // output triangle) averaged field values.  Pre-v9.1.0 this emitted
      // per-VERTEX values drawn from per-QP DOFs at reference-triangle
      // corner positions, which is incorrect when the QPs live at interior
      // barycentric points (TPV102 wave operator, order=1 triangle rule =
      // `(1/6,1/6), (2/3,1/6), (1/6,2/3)`) — see plan v9.1.0 §2.3 R-001
      // / H-V91-A4.  The per-cell average is exact for any barycentric
      // quadrature of degree ≤ linear and side-steps the BR2 latent R-005
      // cross-face read (the loop bound is now `k<nbf`, not hardcoded 3).
      std::vector<std::array<double,3>> vertices;
      std::vector<std::array<int,3>> triangles;
      // Per-cell field values (one entry per output triangle)
      std::vector<double> c_sd, c_ss, c_srd, c_srs, c_td, c_ts, c_psi, c_sn;
      std::vector<double> c_a, c_Dc, c_x2, c_x3;

      auto process_face = [&](int fi, int face_mesh_idx, bool is_shared)
      {
         FaceElementTransformations *FTr = nullptr;
         if (!is_shared)
         {
            int face = fault_interior_faces_[fi];
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }
         else
         {
            if constexpr (std::is_same_v<MeshType, ParMesh>)
            {
#ifdef MFEM_USE_MPI
               int sf = fi;
               FTr = mesh_.GetSharedFaceTransformations(sf);
#endif
            }
         }
         if (!FTr) { return; }

         // Reference triangle vertices: (0,0), (1,0), (0,1).  Kept for
         // geometry; the per-cell data array below holds the single
         // averaged value for this triangle.
         int base_vert = static_cast<int>(vertices.size());
         const double ref_tri[3][2] = {{0,0}, {1,0}, {0,1}};
         for (int v = 0; v < 3; v++)
         {
            IntegrationPoint ip;
            ip.Set2(ref_tri[v][0], ref_tri[v][1]);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);
            vertices.push_back({coords(0), coords(1), coords(2)});
         }
         triangles.push_back({base_vert, base_vert+1, base_vert+2});

         // Face-averaged cell-data values (over this face's nbf DOFs).
         // Loop bound `k<nbf` handles every DG method:
         //   BR2 (nbf=1): single centroid value, no mixing (R-005 fix).
         //   IP  (nbf=3): mean of three vertex values.
         //   TPV102 wave operator (nbf=3): mean of three interior QPs.
         //   Higher order (nbf=6, ...): mean of all DOFs (R-002 fix).
         // (nbf>0 enforced at function entry via MFEM_VERIFY — R-007.)
         const int base = face_mesh_idx * nbf;
         double a_sd = 0.0, a_ss = 0.0, a_srd = 0.0, a_srs = 0.0;
         double a_td = 0.0, a_ts = 0.0, a_psi = 0.0, a_sn = 0.0;
         double a_a  = 0.0, a_Dc = 0.0, a_x2  = 0.0, a_x3 = 0.0;
         for (int k = 0; k < nbf; k++)
         {
            const int d = base + k;
            a_sd  += local_slip(2*d);
            a_ss  += local_slip(2*d+1);
            a_srd += local_slip_rate(2*d);
            a_srs += local_slip_rate(2*d+1);
            a_td  += local_traction(2*d);
            a_ts  += local_traction(2*d+1);
            a_psi += local_state(d);
            a_sn  += has_normal ? local_normal_stress(d) : 0.0;
            a_a   += local_a.Size()  > 0 ? local_a(d)  : 0.0;
            a_Dc  += local_Dc.Size() > 0 ? local_Dc(d) : 0.0;
            a_x2  += local_x2.Size() > 0 ? local_x2(d) : 0.0;
            a_x3  += local_x3.Size() > 0 ? local_x3(d) : 0.0;
         }
         const double inv_nbf = 1.0 / static_cast<double>(nbf);
         c_sd.push_back(a_sd  * inv_nbf);
         c_ss.push_back(a_ss  * inv_nbf);
         c_srd.push_back(a_srd * inv_nbf);
         c_srs.push_back(a_srs * inv_nbf);
         c_td.push_back(a_td  * inv_nbf);
         c_ts.push_back(a_ts  * inv_nbf);
         c_psi.push_back(a_psi * inv_nbf);
         c_sn.push_back(a_sn  * inv_nbf);
         c_a.push_back(a_a   * inv_nbf);
         c_Dc.push_back(a_Dc  * inv_nbf);
         c_x2.push_back(a_x2  * inv_nbf);
         c_x3.push_back(a_x3  * inv_nbf);
      };

      // Interior faces
      for (int i = 0; i < n_int; i++)
      {
         process_face(i, i, false);
      }
      // Shared faces
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < n_shared; i++)
         {
            int sf_idx = fault_shared_faces_[i];
            process_face(sf_idx, n_int + i, true);
         }
#endif
      }

      // Write VTU (create directory if needed)
      std::string fault_dir = prefix + "/FaultSurface";
      if (rank == 0)
      {
         ::mkdir(fault_dir.c_str(), 0755);  // ignore error if exists
      }
#ifdef MFEM_USE_MPI
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
         MPI_Barrier(mesh_.GetComm());
      }
#endif
      std::string vtu_name = fault_dir + "/fault_surface_r"
                           + std::to_string(rank)
                           + "_c" + std::to_string(cycle) + ".vtu";
      std::ofstream vtu(vtu_name);
      vtu << std::setprecision(10);
      int npts = static_cast<int>(vertices.size());
      int ncells = static_cast<int>(triangles.size());

      vtu << "<?xml version=\"1.0\"?>\n";
      vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      vtu << "<UnstructuredGrid>\n";
      vtu << "<Piece NumberOfPoints=\"" << npts
          << "\" NumberOfCells=\"" << ncells << "\">\n";

      // Points
      vtu << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             "format=\"ascii\">\n";
      for (auto &v : vertices)
      {
         vtu << v[0] << " " << v[1] << " " << v[2] << "\n";
      }
      vtu << "</DataArray></Points>\n";

      // Cells
      vtu << "<Cells>\n";
      vtu << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
      for (auto &t : triangles)
      {
         vtu << t[0] << " " << t[1] << " " << t[2] << "\n";
      }
      vtu << "</DataArray>\n";
      vtu << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; i++) { vtu << (i+1)*3 << "\n"; }
      vtu << "</DataArray>\n";
      vtu << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; i++) { vtu << 5 << "\n"; } // VTK_TRIANGLE
      vtu << "</DataArray>\n";
      vtu << "</Cells>\n";

      // Cell data — one value per output triangle (speckle-free per
      // R-001 / v9.1.0 §2.3 at the cost of sub-triangle shading, which
      // the DG face-local basis doesn't carry anyway).
      vtu << "<CellData>\n";
      auto write_field = [&](const char* name, const std::vector<double>& vals) {
         vtu << "<DataArray type=\"Float64\" Name=\"" << name
             << "\" format=\"ascii\">\n";
         for (double v : vals) { vtu << v << "\n"; }
         vtu << "</DataArray>\n";
      };
      write_field("slip_dip",         c_sd);
      write_field("slip_strike",      c_ss);
      write_field("slip_rate_dip",    c_srd);
      write_field("slip_rate_strike", c_srs);
      write_field("traction_dip",     c_td);
      write_field("traction_strike",  c_ts);
      write_field("state_variable",   c_psi);
      write_field("normal_stress",    c_sn);
      write_field("param_a",          c_a);
      write_field("param_Dc",         c_Dc);
      write_field("fault_x2",         c_x2);
      write_field("fault_x3",         c_x3);
      vtu << "</CellData>\n";

      vtu << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
      vtu.close();

      // Rank 0 writes PVTU index + appends to PVD time series
      if (rank == 0)
      {
         std::string pvtu_rel = "fault_surface_c"
                              + std::to_string(cycle) + ".pvtu";
         std::string pvtu_name = fault_dir + "/" + pvtu_rel;
         std::ofstream pvtu(pvtu_name);
         pvtu << "<?xml version=\"1.0\"?>\n";
         pvtu << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n";
         pvtu << "<PUnstructuredGrid GhostLevel=\"0\">\n";
         pvtu << "<PPoints><PDataArray type=\"Float64\" "
                 "NumberOfComponents=\"3\"/></PPoints>\n";
         pvtu << "<PCellData>\n";
         const char* fields[] = {
            "slip_dip","slip_strike","slip_rate_dip","slip_rate_strike",
            "traction_dip","traction_strike","state_variable","normal_stress",
            "param_a","param_Dc","fault_x2","fault_x3"
         };
         for (auto f : fields) {
            pvtu << "<PDataArray type=\"Float64\" Name=\"" << f << "\"/>\n";
         }
         pvtu << "</PCellData>\n";
         for (int r = 0; r < nranks; r++) {
            pvtu << "<Piece Source=\"fault_surface_r" << r
                 << "_c" << cycle << ".vtu\"/>\n";
         }
         pvtu << "</PUnstructuredGrid>\n</VTKFile>\n";
         pvtu.close();

         // Append to PVD time-series file
         fault_pvd_entries_.push_back({time, pvtu_rel});
         WriteFaultPVD(fault_dir);
      }
   }

   // Stored shared face indices for surface VTU output
   Array<int> fault_shared_faces_;
   Array<int> fault_interior_faces_;

   // ---------------------------------------------------------------
   //  Save methods
   // ---------------------------------------------------------------

   /// Save based on adaptive V_max schedule.  All ranks must call
   /// collectively (ParaViewDataCollection::Save is MPI-collective).
   /// @param V_max  Global maximum slip rate (must be globally reduced
   ///               BEFORE calling so all ranks see the same value).
   bool Save(int cycle, real_t time, real_t V_max)
   {
      if (output_every_n_steps > 0)
      {
         // Step-based: write only at multiples of the interval
         if (cycle % output_every_n_steps == 0)
         {
            last_v_max_ = V_max;
            return ForceSaveImpl(cycle, time);
         }
         return false;
      }
      // Time-based: fixed dt or adaptive V_max schedule (with hysteresis).
      const int new_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : adaptive_.Interval(V_max, new_regime);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }
      current_regime_ = new_regime;
      last_v_max_     = V_max;
      return ForceSaveImpl(cycle, time);
   }

   /// Schedule check without writing the volume PVD.  Returns true on the
   /// cycles/times when Save() would write, and advances last_write_time_
   /// (and current_regime_, last_v_max_) so subsequent scheduling stays
   /// consistent.  Use this when only the fault-surface VTU is wanted and
   /// the volume mesh+fields are suppressed to save disk space.
   bool ShouldWrite(int cycle, real_t time, real_t V_max)
   {
      if (output_every_n_steps > 0)
      {
         if (cycle % output_every_n_steps == 0)
         {
            last_write_time_ = time;
            last_v_max_      = V_max;
            return true;
         }
         return false;
      }
      const int new_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : adaptive_.Interval(V_max, new_regime);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }
      current_regime_  = new_regime;
      last_write_time_ = time;
      last_v_max_      = V_max;
      return true;
   }

   /// Read-only schedule check: returns true iff Save() (or ShouldWrite())
   /// would write at (cycle, time, V_max), without mutating last_write_time_,
   /// current_regime_, or last_v_max_.  Use this to gate expensive per-step
   /// packing work: then call the two-arg CommitSchedule(time, V_max) to
   /// advance state after the writes are committed.  The single-arg shim
   /// CommitSchedule(time) uses the last V_max recorded by Save/ShouldWrite
   /// (NOT by Peek — Peek is strictly read-only) and is correct for any
   /// call site that uses the default hysteresis_factor = 1.0 (where the
   /// regime is a stateless function of V).  Fault-only callers that set
   /// hysteresis_factor > 1 MUST use the two-arg CommitSchedule.
   bool PeekShouldWrite(int cycle, real_t time, real_t V_max) const
   {
      if (output_every_n_steps > 0)
      {
         return (cycle % output_every_n_steps == 0);
      }
      // Compute the prospective regime locally WITHOUT mutating state.
      const int prospective_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : adaptive_.Interval(V_max, prospective_regime);
      return (time - last_write_time_ >= dt_out * kOutputTimeTolerance);
   }

   /// Advance `last_write_time_` AND `current_regime_` AND `last_v_max_`
   /// in one atomic step after PeekShouldWrite has confirmed a write.
   /// Use this two-arg form from callers that skip Save (e.g. the
   /// fault-only path that only writes the fault-surface VTU) so the
   /// regime state machine keeps advancing — otherwise
   /// `current_regime_` would be pinned at 0 and hysteresis would
   /// silently no-op.
   void CommitSchedule(real_t time, real_t V_max)
   {
      current_regime_  = adaptive_.NextRegime(V_max, current_regime_);
      last_write_time_ = time;
      last_v_max_      = V_max;
   }

   /// Back-compat single-arg shim preserved for existing call sites
   /// (TPV102 driver and legacy BP5 paths).  Delegates to the two-arg
   /// form using the last V_max recorded by Save/ShouldWrite; this is
   /// correct for callers that either (a) use the default
   /// hysteresis_factor = 1.0 (regime is stateless in V), or (b) call
   /// Save/ShouldWrite just before this.  Callers that want hysteresis
   /// to advance through a PeekShouldWrite-gated path MUST call the
   /// two-arg overload instead.
   void CommitSchedule(real_t time) { CommitSchedule(time, last_v_max_); }

   /// Force a save at the current state.
   void ForceSave(int cycle, real_t time)
   {
      ForceSaveImpl(cycle, time);
   }

   /// Adaptive output interval based on maximum slip rate (legacy API,
   /// still used by test_io.cpp and any caller that wants the default
   /// thresholds/intervals without constructing a ParaViewOutput).  For
   /// runtime-configurable thresholds, mutate GetSchedule() instead.
   static real_t OutputInterval(real_t V_max)
   {
      AdaptiveSchedule s;  // defaults = legacy behavior
      return s.Interval(V_max, s.NextRegime(V_max, 0));
   }

   void SetDataFormat(VTKFormat fmt) { pv_.SetDataFormat(fmt); }
   void SetHighOrderOutput(bool enable) { pv_.SetHighOrderOutput(enable); }
   void SetLevelsOfDetail(int lod) { pv_.SetLevelsOfDetail(lod); }
   bool HasFaultOutput() const { return has_fault_output_; }

private:
   MeshType &mesh_;
   int order_;
   real_t last_write_time_;
   ParaViewDataCollection pv_;

   // Adaptive-schedule state.  adaptive_ holds the tunable thresholds
   // and intervals (GetSchedule() exposes it for CLI overrides);
   // current_regime_ threads hysteresis across writes; last_v_max_
   // backs the single-arg CommitSchedule(time) shim.
   AdaptiveSchedule adaptive_;
   int    current_regime_ = 0;   // 0=interseismic, 1=nucleation, 2=coseismic
   real_t last_v_max_     = 0.0;

   // Fault L2-p0 output
   bool has_fault_output_ = false;
   int  nbf_per_face_ = 1;
   int  n_interior_fault_faces_ = 0;
   int  n_shared_fault_faces_ = 0;

   std::unique_ptr<FiniteElementCollection> fault_fec_;
   std::unique_ptr<FES> fault_fes_;

   std::unique_ptr<GF> fault_slip_dip_;
   std::unique_ptr<GF> fault_slip_strike_;
   std::unique_ptr<GF> fault_slip_rate_dip_;
   std::unique_ptr<GF> fault_slip_rate_strike_;
   std::unique_ptr<GF> fault_trac_dip_;
   std::unique_ptr<GF> fault_trac_strike_;
   std::unique_ptr<GF> fault_state_;
   std::unique_ptr<GF> fault_normal_stress_;

   // Friction parameter and coordinate fields (static, written once)
   std::unique_ptr<GF> fault_param_a_;
   std::unique_ptr<GF> fault_param_Dc_;
   std::unique_ptr<GF> fault_coord_x2_;
   std::unique_ptr<GF> fault_coord_x3_;

   // Interior fault face → element mapping
   std::vector<int> fault_face_elem1_;
   std::vector<int> fault_face_elem2_;
   // Shared fault face → local element mapping (Elem1 only)
   std::vector<int> fault_shared_elem1_;

   // Fault surface PVD time-series entries (rank 0 only)
   struct PVDEntry { real_t time; std::string pvtu_file; };
   std::vector<PVDEntry> fault_pvd_entries_;

   bool ForceSaveImpl(int cycle, real_t time)
   {
      pv_.SetCycle(cycle);
      pv_.SetTime(time);
      pv_.Save();
      last_write_time_ = time;
      return true;
   }

   /// Write (or overwrite) the fault surface PVD file with all entries so far.
   void WriteFaultPVD(const std::string &fault_dir)
   {
      std::string pvd_name = fault_dir + "/fault_surface.pvd";
      std::ofstream pvd(pvd_name, std::ios::trunc);
      pvd << std::setprecision(17);
      pvd << "<?xml version=\"1.0\"?>\n";
      pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
      pvd << "<Collection>\n";
      for (const auto &e : fault_pvd_entries_)
      {
         pvd << "<DataSet timestep=\"" << e.time
             << "\" file=\"" << e.pvtu_file << "\"/>\n";
      }
      pvd << "</Collection>\n</VTKFile>\n";
      pvd.close();
   }

   // Factory for FES: serial vs parallel
   static std::unique_ptr<FES> MakeFES(MeshType &mesh,
                                       FiniteElementCollection *fec)
   {
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         return std::make_unique<ParFiniteElementSpace>(&mesh, fec);
#else
         MFEM_ABORT("ParMesh requires MFEM_USE_MPI");
         return nullptr;
#endif
      }
      else
      {
         return std::make_unique<FiniteElementSpace>(&mesh, fec);
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PARAVIEW_OUTPUT_HPP
