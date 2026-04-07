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
   //  Fault surface VTU output (proper 2D face geometry)
   // ---------------------------------------------------------------

   /// @brief Write fault surface as a triangle VTU with per-cell field data.
   ///
   /// Each rank writes fault_surface_r{rank}_c{cycle}.vtu containing the
   /// actual fault face triangles and per-face averaged field values.
   /// Rank 0 also writes a .pvtu index.  This replaces the L2-p0 volume
   /// projection which creates scattered-dot artifacts in ParaView.
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

      // Collect face vertices and per-face averaged fields
      std::vector<std::array<double,3>> vertices;
      std::vector<std::array<int,3>> triangles;
      // Per-face field values
      std::vector<double> f_sd, f_ss, f_srd, f_srs, f_td, f_ts, f_psi, f_sn;
      std::vector<double> f_a, f_Dc, f_x2, f_x3;

      auto process_face = [&](int fi, int face_mesh_idx, bool is_shared)
      {
         // Get face vertex coordinates
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
               int sf = fi;  // shared face index
               FTr = mesh_.GetSharedFaceTransformations(sf);
#endif
            }
         }
         if (!FTr) { return; }

         // Get face vertex positions (3 vertices for a triangle)
         int base_vert = static_cast<int>(vertices.size());
         const IntegrationRule &nir =
            IntRules.Get(FTr->FaceGeom, 1);  // linear nodes
         for (int v = 0; v < nbf; v++)
         {
            const IntegrationPoint &ip = nir.IntPoint(v);
            FTr->Face->SetIntPoint(&ip);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);
            vertices.push_back({coords(0), coords(1), coords(2)});
         }
         // Triangle connectivity
         if (nbf == 3)
         {
            triangles.push_back({base_vert, base_vert+1, base_vert+2});
         }
         else
         {
            // Fallback for non-triangle faces (shouldn't happen for tet mesh)
            for (int v = 1; v < nbf - 1; v++)
            {
               triangles.push_back({base_vert, base_vert+v, base_vert+v+1});
            }
         }

         // Average DOF values for this face
         int base = face_mesh_idx * nbf;
         double sd=0,ss=0,srd=0,srs=0,td=0,ts=0,psi=0,sn=0;
         double a=0,Dc=0,x2=0,x3=0;
         for (int k = 0; k < nbf; k++)
         {
            int d = base + k;
            sd  += local_slip(2*d);
            ss  += local_slip(2*d+1);
            srd += local_slip_rate(2*d);
            srs += local_slip_rate(2*d+1);
            td  += local_traction(2*d);
            ts  += local_traction(2*d+1);
            psi += local_state(d);
            if (has_normal) { sn += local_normal_stress(d); }
            if (local_a.Size() > 0) { a += local_a(d); }
            if (local_Dc.Size() > 0) { Dc += local_Dc(d); }
            if (local_x2.Size() > 0) { x2 += local_x2(d); }
            if (local_x3.Size() > 0) { x3 += local_x3(d); }
         }
         double inv = 1.0 / nbf;
         f_sd.push_back(sd*inv);   f_ss.push_back(ss*inv);
         f_srd.push_back(srd*inv); f_srs.push_back(srs*inv);
         f_td.push_back(td*inv);   f_ts.push_back(ts*inv);
         f_psi.push_back(psi*inv); f_sn.push_back(sn*inv);
         f_a.push_back(a*inv);     f_Dc.push_back(Dc*inv);
         f_x2.push_back(x2*inv);   f_x3.push_back(x3*inv);
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
         // Simple mkdir -p equivalent (safe to call redundantly)
         std::string cmd = "mkdir -p " + fault_dir;
         (void)system(cmd.c_str());
      }
#ifdef MFEM_USE_MPI
      MPI_Barrier(mesh_.GetComm());
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

      // Cell data
      vtu << "<CellData>\n";
      auto write_field = [&](const char* name, const std::vector<double>& vals) {
         vtu << "<DataArray type=\"Float64\" Name=\"" << name
             << "\" format=\"ascii\">\n";
         for (double v : vals) { vtu << v << "\n"; }
         vtu << "</DataArray>\n";
      };
      write_field("slip_dip", f_sd);
      write_field("slip_strike", f_ss);
      write_field("slip_rate_dip", f_srd);
      write_field("slip_rate_strike", f_srs);
      write_field("traction_dip", f_td);
      write_field("traction_strike", f_ts);
      write_field("state_variable", f_psi);
      write_field("normal_stress", f_sn);
      write_field("param_a", f_a);
      write_field("param_Dc", f_Dc);
      write_field("fault_x2", f_x2);
      write_field("fault_x3", f_x3);
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
