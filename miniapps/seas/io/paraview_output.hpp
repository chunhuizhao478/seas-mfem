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
/// Output frequency is adaptive, keyed to V_max (slip rate):
///   - coseismic  (V > 1e-3 m/s)  : every 0.01 s
///   - nucleation (V > 1e-6)       : every 1.0 s
///   - interseismic                : every 1 year
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

      // Build fault face → (elem1, elem2) mapping for interior faces
      int n_int = fault_interior_faces.Size();
      fault_face_elem1_.resize(n_int);
      fault_face_elem2_.resize(n_int);
      for (int i = 0; i < n_int; i++)
      {
         int face = fault_interior_faces[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         fault_face_elem1_[i] = FTr->Elem1No;
         fault_face_elem2_[i] = FTr->Elem2No;
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
         (*fault_param_a_)(e1)  = avg_a;  (*fault_param_a_)(e2)  = avg_a;
         (*fault_param_Dc_)(e1) = avg_Dc; (*fault_param_Dc_)(e2) = avg_Dc;
         (*fault_coord_x2_)(e1) = avg_x2; (*fault_coord_x2_)(e2) = avg_x2;
         (*fault_coord_x3_)(e1) = avg_x3; (*fault_coord_x3_)(e2) = avg_x3;
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
         (*fault_slip_dip_)(e1)         = s_d;
         (*fault_slip_dip_)(e2)         = s_d;
         (*fault_slip_strike_)(e1)      = s_s;
         (*fault_slip_strike_)(e2)      = s_s;
         (*fault_slip_rate_dip_)(e1)    = sr_d;
         (*fault_slip_rate_dip_)(e2)    = sr_d;
         (*fault_slip_rate_strike_)(e1) = sr_s;
         (*fault_slip_rate_strike_)(e2) = sr_s;
         (*fault_trac_dip_)(e1)         = tr_d;
         (*fault_trac_dip_)(e2)         = tr_d;
         (*fault_trac_strike_)(e1)      = tr_s;
         (*fault_trac_strike_)(e2)      = tr_s;
         (*fault_state_)(e1)            = psi;
         (*fault_state_)(e2)            = psi;
         if (has_normal)
         {
            (*fault_normal_stress_)(e1) = sn;
            (*fault_normal_stress_)(e2) = sn;
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
            return ForceSaveImpl(cycle, time);
         }
         return false;
      }
      // Time-based: fixed dt or adaptive V_max schedule
      real_t dt_out = (fixed_dt > 0.0) ? fixed_dt : OutputInterval(V_max);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }
      return ForceSaveImpl(cycle, time);
   }

   /// Force a save at the current state.
   void ForceSave(int cycle, real_t time)
   {
      ForceSaveImpl(cycle, time);
   }

   /// Adaptive output interval based on maximum slip rate.
   static real_t OutputInterval(real_t V_max)
   {
      if (V_max > 1e-3)
      {
         return 0.01;  // Coseismic: every 0.01 s
      }
      else if (V_max > 1e-6)
      {
         return 1.0;   // Nucleation: every 1.0 s
      }
      else
      {
         return 1.0 * BP5Params::seconds_per_year;  // Interseismic: every 1 year
      }
   }

   void SetDataFormat(VTKFormat fmt) { pv_.SetDataFormat(fmt); }
   void SetHighOrderOutput(bool enable) { pv_.SetHighOrderOutput(enable); }
   bool HasFaultOutput() const { return has_fault_output_; }

private:
   MeshType &mesh_;
   int order_;
   real_t last_write_time_;
   ParaViewDataCollection pv_;

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

   bool ForceSaveImpl(int cycle, real_t time)
   {
      pv_.SetCycle(cycle);
      pv_.SetTime(time);
      pv_.Save();
      last_write_time_ = time;
      return true;
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
