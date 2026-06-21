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

#ifndef MFEM_SEAS_FREE_SURFACE_OUTPUT_HPP
#define MFEM_SEAS_FREE_SURFACE_OUTPUT_HPP

#include "mfem.hpp"
#include "free_surface_vtu_binary.hpp"  // fsvtu:: single-file VTU writer (Vtu mode)
#include <string>
#include <memory>
#include <filesystem>
#include <cstring>   // std::memcpy (UpdateVelocity); not relied on transitively
#include <cstdio>    // std::rename (atomic PVD rewrite, Vtu mode)
#include <sys/stat.h> // ::mkdir (rank-0 Vtu output dir)
#include <vector>
#include <utility>   // std::pair (PVD entries)
#include <fstream>   // std::ofstream (PVD)
#include <iomanip>   // std::setprecision (PVD)

#ifdef MFEM_USE_MPI

namespace mfem
{
namespace seas
{

/// @brief Writes the free-surface (z = z_max) boundary of a 3D SAFS/TPV domain
///        as a codim-1 ParaView slice carrying the volume `velocity` field.
///
/// The slice is built once with `mfem::ParSubMesh::CreateFromBoundary` on the
/// free-surface boundary attribute(s) and refreshed each output step through a
/// cached `mfem::ParTransferMap`.  It writes on its own fixed-dt schedule,
/// independent of the volume/bulk/fault collections and of the
/// `paraview_enabled` master gate (the driver wires this up in Phase 3).
///
/// Velocity transfer reuses the same byNODES, vdim=3 GaussLobatto L2 layout the
/// volume path uses (spatial_dyn_driver.cpp:2310), so the velocity DOF block can
/// be `memcpy`-ed straight out of the state vector `Q` at offset
/// `VX * ndof_total` (length `3 * ndof_total`).  The parent L2 spaces are forced
/// to `BasisType::GaussLobatto`: the boundary->submesh L2 transfer asserts the
/// *parent* DG basis is GaussLobatto (`submesh_utils.cpp:117`).
///
/// All MFEM-collective calls (`CreateFromBoundary`, `ParTransferMap::Transfer`,
/// `DataCollection::Save`) must be reached by EVERY rank — including ranks with
/// zero local free-surface faces (an empty local submesh is valid).
class FreeSurfaceOutput
{
public:
   /// Output back end: per-rank VTU/PVTU/PVD or a single `.vtkhdf` file.
   enum class Mode { Vtu, Hdf5 };

   /// @param prefix           Output directory prefix
   ///                         (e.g. "<out>/ParaView_free_surface").
   /// @param parent           Parent volume ParMesh (must outlive this object).
   /// @param free_surface_attrs Boundary attribute(s) of the free surface;
   ///                         non-empty (validated; caller resolves + guards).
   /// @param order            Polynomial order of the parent velocity space
   ///                         (matches the wave-operator scalar order).
   /// @param rank             MPI rank, written as the static `mpi_rank` field.
   /// @param fixed_dt         Fixed output cadence (seconds, > 0).
   /// @param collection_name  ParaView collection basename.
   /// @param mode             `Mode::Vtu` (default) or `Mode::Hdf5`.
   FreeSurfaceOutput(const std::string &prefix,
                     ParMesh           &parent,
                     const Array<int>  &free_surface_attrs,
                     int                order,
                     int                rank,
                     real_t             fixed_dt,
                     const std::string &collection_name = "free_surface",
                     Mode               mode = Mode::Vtu)
      // 1. Boundary submesh (collective).  Constructed in the init list via the
      //    helper so C++17 guaranteed copy elision builds it IN PLACE, never
      //    through ParSubMesh's (implicit) move ctor.  std::make_unique<
      //    ParSubMesh>(CreateFromBoundary(...)) would force that move, and
      //    MFEM's Mesh/ParMesh move semantics are incomplete (the plan's Risk
      //    Assessment flags this); copy elision sidesteps the concern entirely.
      //    Ranks with no local free-surface faces get an empty local submesh
      //    and must still construct (collective).
      : submesh_(MakeBoundarySubmesh(parent, free_surface_attrs))
   {
      // GlobalNE via Allreduce(SUM) — NOT Reduce — so it is valid on EVERY rank
      // (the -np 2 Phase 4 test asserts it on all ranks, not just rank 0).
      long long local_ne = submesh_.GetNE();
      MPI_Allreduce(&local_ne, &global_ne_, 1, MPI_LONG_LONG, MPI_SUM,
                    parent.GetComm());

      // 2. Parent velocity space — MUST match the Q velocity block layout:
      //    L2 GaussLobatto, vdim=3, byNODES (mirrors spatial_dyn_driver.cpp:2310).
      parent_vel_fec_ = std::make_unique<L2_FECollection>(
                           order, 3, BasisType::GaussLobatto);
      parent_vel_fes_ = std::make_unique<ParFiniteElementSpace>(
                           &parent, parent_vel_fec_.get(), 3, Ordering::byNODES);
      parent_vel_gf_  = std::make_unique<ParGridFunction>(parent_vel_fes_.get());
      *parent_vel_gf_ = 0.0;

      // 3. Parent rank space — L2 order-0, GaussLobatto (NOT the default
      //    GaussLegendre at fe_coll.hpp:364): the boundary->submesh L2 transfer
      //    asserts the *parent* DG basis is GaussLobatto (submesh_utils.cpp:117).
      //    This intentionally differs from the volume-path pv_rank_fec
      //    (driver L2318, GaussLegendre) because that path never transfers
      //    through a submesh.
      parent_rank_fec_ = std::make_unique<L2_FECollection>(
                            0, 3, BasisType::GaussLobatto);
      parent_rank_fes_ = std::make_unique<ParFiniteElementSpace>(
                            &parent, parent_rank_fec_.get());
      parent_rank_gf_  = std::make_unique<ParGridFunction>(parent_rank_fes_.get());
      *parent_rank_gf_ = static_cast<real_t>(rank);

      // 4. Submesh spaces.  The submesh is codim-1, so its L2 collection's
      //    DIMENSION argument must be the SUBMESH dimension (parent_dim - 1),
      //    NOT the parent's 3 — this mirrors MFEM's own boundary-transfer test
      //    (tests/unit/mesh: sub fec uses `submesh->Dimension()`).  Passing 3
      //    here builds a collection whose trace bookkeeping is inconsistent with
      //    the 2D submesh and SIGSEGVs inside SubMeshUtils::BuildVdofToVdofMap.
      //    The vdim stays 3 (a 3-component velocity living on the 2D surface).
      const int sub_dim = submesh_.Dimension();
      sub_vel_fec_ = std::make_unique<L2_FECollection>(
                        order, sub_dim, BasisType::GaussLobatto);
      sub_vel_fes_ = std::make_unique<ParFiniteElementSpace>(
                        &submesh_, sub_vel_fec_.get(), 3, Ordering::byNODES);
      sub_vel_gf_  = std::make_unique<ParGridFunction>(sub_vel_fes_.get());
      *sub_vel_gf_ = 0.0;

      sub_rank_fec_ = std::make_unique<L2_FECollection>(
                         0, sub_dim, BasisType::GaussLobatto);
      sub_rank_fes_ = std::make_unique<ParFiniteElementSpace>(
                         &submesh_, sub_rank_fec_.get());
      sub_rank_gf_  = std::make_unique<ParGridFunction>(sub_rank_fes_.get());
      *sub_rank_gf_ = 0.0;

      // 5. Transfer maps — built once, reused every output step.
      vel_map_  = std::make_unique<ParTransferMap>(
                     ParSubMesh::CreateTransferMap(*parent_vel_gf_, *sub_vel_gf_));
      rank_map_ = std::make_unique<ParTransferMap>(
                     ParSubMesh::CreateTransferMap(*parent_rank_gf_,
                                                   *sub_rank_gf_));

      // 6. The rank field is static — transfer it once now.
      rank_map_->Transfer(*parent_rank_gf_, *sub_rank_gf_);

      // 7. Output back end.
      //    Mode::Hdf5 — MFEM ParaViewHDFDataCollection (single .vtkhdf/run).
      //    Mode::Vtu  — our consolidating single-file writer: ONE
      //                 `<collection>_c<cycle>.vtu` + a cumulative
      //                 `<collection>.pvd` per cycle (rank-0 gather), exactly
      //                 like the fault `fault_surface_c<cycle>.vtu` output —
      //                 NO per-rank `proc*.vtu` explosion.  pv_dc_ stays null in
      //                 Vtu mode; Save() drives WriteVtuCycle() instead.
      prefix_          = prefix;
      collection_name_ = collection_name;
      rank_            = rank;
      mode_            = mode;
#ifdef MFEM_USE_MPI
      comm_ = parent.GetComm();
      MPI_Comm_size(comm_, &nranks_);
#endif
      if (mode == Mode::Hdf5)
      {
#ifdef MFEM_USE_HDF5
         pv_dc_ = std::make_unique<ParaViewHDFDataCollection>(
                     collection_name, &submesh_);
         pv_dc_->SetPrefixPath(prefix);
         pv_dc_->SetDataFormat(VTKFormat::BINARY);
         pv_dc_->SetHighOrderOutput(true);
         pv_dc_->SetLevelsOfDetail(order);
         pv_dc_->RegisterField("velocity", sub_vel_gf_.get());
         pv_dc_->RegisterField("mpi_rank", sub_rank_gf_.get());
#else
         MFEM_ABORT("FreeSurfaceOutput Mode::Hdf5 needs MFEM_USE_HDF5=YES");
#endif
      }
      // Mode::Vtu: no DataCollection — WriteVtuCycle() handles it in Save().

      // 8. Schedule state.
      fixed_dt_        = fixed_dt;
      last_write_time_ = -1e30;
   }

   // Non-copyable AND non-movable: pv_dc_ and the submesh FE spaces store
   // interior pointers into submesh_ (a value member), so the object must never
   // relocate.  The driver holds it via unique_ptr<FreeSurfaceOutput>, so this
   // does not impede make_unique / reset / the unique_ptr's own move.
   FreeSurfaceOutput(const FreeSurfaceOutput &)            = delete;
   FreeSurfaceOutput &operator=(const FreeSurfaceOutput &) = delete;
   FreeSurfaceOutput(FreeSurfaceOutput &&)                 = delete;
   FreeSurfaceOutput &operator=(FreeSurfaceOutput &&)      = delete;

   /// Fixed-dt schedule check (read-only).  Mirrors the volume writer's 0.99
   /// tolerance so a frame due at `n*dt` is not suppressed by rounding.
   bool ShouldWrite(real_t time) const
   {
      return (time - last_write_time_) >= fixed_dt_ * kTol;
   }

   /// Refresh the sliced velocity from a contiguous [VX|VY|VZ] block.
   ///
   /// `vxvyvz_block` is the byNODES, vdim=3 velocity block of the state vector
   /// `Q` (Q.GetData() + VX*ndof_total); `ndof_per_component` is `ndof_total`.
   /// The byNODES layout makes this block bit-identical to the parent velocity
   /// GridFunction data, so a single memcpy refreshes it before transfer.
   void UpdateVelocity(const real_t *vxvyvz_block, int ndof_per_component)
   {
      MFEM_VERIFY(parent_vel_fes_->GetNDofs() == ndof_per_component,
                  "FreeSurfaceOutput: parent velocity NDofs mismatch");
      std::memcpy(parent_vel_gf_->GetData(), vxvyvz_block,
                  3 * static_cast<size_t>(ndof_per_component) * sizeof(real_t));
      vel_map_->Transfer(*parent_vel_gf_, *sub_vel_gf_);
   }

   /// Write the current slice at (cycle, time).  Collective (the Vtu-mode gather
   /// and the Hdf5-mode Save are both MPI-collective — every rank must call).
   void Save(int cycle, real_t time)
   {
      if (mode_ == Mode::Vtu)
      {
         WriteVtuCycle(cycle, time);     // single <coll>_c<cycle>.vtu + PVD
      }
      else
      {
         MFEM_VERIFY(pv_dc_, "FreeSurfaceOutput::Save: null collection");
         pv_dc_->SetCycle(cycle);
         pv_dc_->SetTime(time);
         pv_dc_->Save();
      }
      last_write_time_ = time;
   }

   /// Total surface element count across all ranks (valid on every rank).
   long long GlobalNE() const { return global_ne_; }
   /// Submesh topological dimension (2 for a surface extracted from a volume).
   int Dimension() const { return submesh_.Dimension(); }
   /// Read-only handle on the underlying ParaView collection (tests/banner).
   const ParaViewDataCollectionBase *GetDataCollection() const
   { return pv_dc_.get(); }
   /// Mutable submesh handle (tests).
   ParSubMesh &SubMesh() { return submesh_; }
   /// Read-only sliced velocity GridFunction (Phase 4 trace test).
   const ParGridFunction &SubVelocity() const { return *sub_vel_gf_; }

private:
   static constexpr real_t kTol = 0.99;

   // Build the boundary submesh, enforcing the non-empty contract.  Returns by
   // value so the ctor init-list use is a guaranteed-copy-elision in-place
   // construction — never a ParSubMesh move (see the ctor note).
   static ParSubMesh MakeBoundarySubmesh(ParMesh          &parent,
                                         const Array<int> &free_surface_attrs)
   {
      MFEM_VERIFY(free_surface_attrs.Size() > 0,
                  "FreeSurfaceOutput: free_surface_attrs must be non-empty "
                  "(caller must resolve + guard before construction)");
      return ParSubMesh::CreateFromBoundary(parent, free_surface_attrs);
   }

   // ---- Vtu mode (single consolidated VTU + PVD per cycle) -----------------

   // Extract this rank's free-surface triangles + per-vertex velocity from the
   // submesh and the sliced velocity GF.  Quad faces are split into two
   // triangles; geometry comes from the element transformation and velocity
   // from `sub_vel_gf_` sampled at each triangle's corner reference points
   // (exact for the order-1 fields these runs use; a linear sampling of a
   // higher-order field otherwise).
   fsvtu::LocalSurfacePack BuildLocalPack()
   {
      fsvtu::LocalSurfacePack p;
      // Corner reference coords, pre-split into triangles per face geometry.
      static const double kTri [1][3][2] = { { {0,0}, {1,0}, {0,1} } };
      static const double kQuad[2][3][2] = { { {0,0}, {1,0}, {1,1} },
                                             { {0,0}, {1,1}, {0,1} } };
      const int ne = submesh_.GetNE();
      for (int e = 0; e < ne; ++e)
      {
         const Geometry::Type geom = submesh_.GetElementBaseGeometry(e);
         int n_sub = 0;
         const double (*tris)[3][2] = nullptr;
         if      (geom == Geometry::TRIANGLE) { n_sub = 1; tris = kTri;  }
         else if (geom == Geometry::SQUARE)   { n_sub = 2; tris = kQuad; }
         else
         {
            MFEM_ABORT("FreeSurfaceOutput Vtu writer: unsupported free-surface "
                       "face geometry " << geom << " (expected TRIANGLE/SQUARE)");
         }
         ElementTransformation *T = submesh_.GetElementTransformation(e);
         for (int t = 0; t < n_sub; ++t)
         {
            for (int j = 0; j < 3; ++j)
            {
               IntegrationPoint ip;
               ip.Set2(tris[t][j][0], tris[t][j][1]);
               T->SetIntPoint(&ip);
               Vector x(3); T->Transform(ip, x);
               Vector v;    sub_vel_gf_->GetVectorValue(e, ip, v);
               p.vertices.push_back({x(0), x(1), x(2)});
               p.velocity.push_back({v(0), v(1), v(2)});
            }
            p.rank_cell.push_back(static_cast<double>(rank_));
         }
      }
      return p;
   }

   // Collective: gather the local pack to rank 0, which writes ONE VTU for this
   // cycle + rewrites the cumulative PVD.  Only rank 0 does I/O after the gather
   // (no trailing barrier — matches the fault binary path).
   void WriteVtuCycle(int cycle, real_t time)
   {
      fsvtu::LocalSurfacePack pack = BuildLocalPack();
#ifdef MFEM_USE_MPI
      fsvtu::GatheredSurfacePack g =
         fsvtu::GatherSurfacePackToRoot(pack, rank_, nranks_, comm_);
#else
      fsvtu::GatheredSurfacePack g =
         fsvtu::GatherSurfacePackToRoot(pack, rank_, 1);
#endif
      if (rank_ == 0)
      {
         ::mkdir(prefix_.c_str(), 0755);   // ignore EEXIST
         const std::string vtu_rel =
            collection_name_ + "_c" + std::to_string(cycle) + ".vtu";
         fsvtu::WriteSurfacePackVTU(prefix_ + "/" + vtu_rel, g,
                                    VTKFormat::BINARY);
         pvd_entries_.push_back({time, vtu_rel});
         WritePVD();
      }
   }

   // Rewrite the cumulative collection PVD (rank 0).  Atomic via .partial +
   // rename so a job killed mid-write never leaves a truncated PVD.
   void WritePVD() const
   {
      const std::string pvd_name = prefix_ + "/" + collection_name_ + ".pvd";
      const std::string tmp_name = pvd_name + ".partial";
      std::ofstream pvd(tmp_name, std::ios::trunc);
      pvd << std::setprecision(17);
      pvd << "<?xml version=\"1.0\"?>\n";
      pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
      for (const auto &e : pvd_entries_)
      {
         pvd << "<DataSet timestep=\"" << e.first
             << "\" file=\"" << e.second << "\"/>\n";
      }
      pvd << "</Collection>\n</VTKFile>\n";
      pvd.close();
      std::rename(tmp_name.c_str(), pvd_name.c_str());   // atomic on POSIX
   }

   // Declaration order == construction order; destruction is reverse order, so
   // pv_dc_ (which stores a Mesh* into submesh_) is destroyed BEFORE submesh_.
   // submesh_ is a value member (NOT a unique_ptr): it must be built by copy
   // elision (see the ctor note), and as the first-declared member it is also
   // destroyed last — after the pv_dc_/spaces/maps that point into it.
   ParSubMesh submesh_;

   std::unique_ptr<L2_FECollection>       parent_vel_fec_;
   std::unique_ptr<ParFiniteElementSpace> parent_vel_fes_;
   std::unique_ptr<ParGridFunction>       parent_vel_gf_;
   std::unique_ptr<L2_FECollection>       parent_rank_fec_;
   std::unique_ptr<ParFiniteElementSpace> parent_rank_fes_;
   std::unique_ptr<ParGridFunction>       parent_rank_gf_;

   std::unique_ptr<L2_FECollection>       sub_vel_fec_;
   std::unique_ptr<ParFiniteElementSpace> sub_vel_fes_;
   std::unique_ptr<ParGridFunction>       sub_vel_gf_;
   std::unique_ptr<L2_FECollection>       sub_rank_fec_;
   std::unique_ptr<ParFiniteElementSpace> sub_rank_fes_;
   std::unique_ptr<ParGridFunction>       sub_rank_gf_;

   std::unique_ptr<ParTransferMap> vel_map_;
   std::unique_ptr<ParTransferMap> rank_map_;

   std::unique_ptr<ParaViewDataCollectionBase> pv_dc_;

   real_t    fixed_dt_        = 0.0;
   real_t    last_write_time_ = -1e30;
   long long global_ne_       = 0;

   // Vtu-mode single-file writer state (unused in Hdf5 mode).
   Mode        mode_ = Mode::Vtu;
   std::string prefix_;
   std::string collection_name_;
   int         rank_   = 0;
#ifdef MFEM_USE_MPI
   MPI_Comm    comm_   = MPI_COMM_NULL;
   int         nranks_ = 1;
#endif
   std::vector<std::pair<real_t, std::string>> pvd_entries_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_FREE_SURFACE_OUTPUT_HPP
