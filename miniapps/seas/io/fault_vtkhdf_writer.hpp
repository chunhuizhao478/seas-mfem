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

#ifndef MFEM_SEAS_FAULT_VTKHDF_WRITER_HPP
#define MFEM_SEAS_FAULT_VTKHDF_WRITER_HPP

// -----------------------------------------------------------------------------
// Phase 2b of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md
//
// Single-file VTKHDF (HDF5) writer for the fault surface.  Replaces
// the per-cycle VTU stream (Phase 1) with one .vtkhdf per run, opened
// in ParaView 5.11+ with native time-series scrubbing.
//
// IMPLEMENTATION DEVIATION FROM PLAN §Phase 2b:
//   The plan specified a ParSubMesh-based per-rank fault submesh +
//   collective parallel HDF5 write.  This file instead reuses the
//   Phase 1 GatherFaultPackToRoot path and emits the .vtkhdf from
//   rank 0 only via `ParaViewHDFDataCollection` constructed with a
//   *serial* `mfem::Mesh` (not `ParMesh`).  Inside MFEM,
//   `ParaViewHDFDataCollection::EnsureVTKHDF` dynamic_casts the mesh:
//   a serial Mesh produces a serial VTKHDF object that performs no
//   MPI calls (`fem/datacollection.cpp:1374-1402`).  This is what we
//   want — only rank 0 ever calls `dc->Save()` and it must not enter
//   any collective HDF5 path.  Justification:
//     - ParSubMesh::CreateFromBoundary requires boundary attributes,
//       but fault faces are interior; promoting them to boundary
//       attributes would invalidate the parent mesh's attribute
//       contract used by the elasticity operator.
//     - The fault surface is a few thousand triangles even for the
//       largest BP5 runs; the parallel-HDF5 write throughput
//       advantage is negligible vs the much simpler reuse of the
//       already-tested Phase 1 gather.
//     - The user-visible result is identical: one .vtkhdf file per
//       run, ParaView-readable, ranks-in-the-same-file-via-gather.
//
// All ranks call collectively (the gather is collective).  Only rank
// 0 owns the file handle / dc lifecycle; non-root ranks return
// immediately after the gather.
// -----------------------------------------------------------------------------

#include "mfem.hpp"
#include "fault_vtu_binary.hpp"   // vtu::LocalFaultPack, GatheredFaultPack,
                                   // GatherFaultPackToRoot

#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include <sys/stat.h>

#ifndef MFEM_USE_HDF5
#error "fault_vtkhdf_writer.hpp requires MFEM_USE_HDF5=YES.  Either enable HDF5 in config/config.mk and rebuild MFEM, or use Phase 1's binary VTU path (the default in this build)."
#endif

#include "fem/datacollection.hpp"  // ParaViewHDFDataCollection
#include "mesh/vtkhdf.hpp"         // VTKHDF (FaultHDFStream)

namespace mfem
{
namespace seas
{
namespace vtkhdf
{

/// The rank-0 serial mesh of the gathered fault triangles (2-D elements in
/// 3-D space), one element per gathered triangle in gathered order.
inline std::unique_ptr<Mesh> BuildFaultTriangleMesh(const vtu::GatheredFaultPack &g)
{
   const int ncells = static_cast<int>(g.triangles.size());
   const int npts   = static_cast<int>(g.vertices.size());

   // Build a serial 2D-in-3D Mesh of triangles.  spaceDim=3 because
   // the fault sits in a 3D embedding space.
   auto mesh = std::make_unique<Mesh>(/*Dim=*/2, /*NVert=*/npts,
                                      /*NElem=*/ncells, /*NBdrElem=*/0,
                                      /*spaceDim=*/3);
   for (int i = 0; i < npts; ++i)
   {
      mesh->AddVertex(g.vertices[i][0], g.vertices[i][1],
                      g.vertices[i][2]);
   }
   const int triangle_attr = 1;
   for (int c = 0; c < ncells; ++c)
   {
      const int t[3] = { g.triangles[c][0], g.triangles[c][1],
                         g.triangles[c][2] };
      mesh->AddTriangle(t, triangle_attr);
   }
   // R-007 (REVIEW.md 2026-04-28): fix_orientation=true so any
   // ParaView consumer that requests vtkPolyDataNormals or surface
   // lighting gets consistent normals.  Costs nothing for cell-data
   // bit-exactness (the data layout is unchanged).
   //
   // R-111: fix_orientation MAY swap the order of vertices within a
   // triangle to enforce a positive Jacobian.  This is invariant for
   // SCALAR cell-data (one value per triangle, the only data this
   // writer emits).  If you ever add per-vertex or per-edge data, OR
   // vector cell-data with an implicit per-vertex direction, this
   // call must be revisited because the gathered (3*c, 3*c+1, 3*c+2)
   // ordering will no longer match the rendered triangle's vertex
   // ordering.
   mesh->FinalizeTriMesh(/*generate_edges=*/0, /*refine=*/0,
                         /*fix_orientation=*/true);
   return mesh;
}

/// @brief Persistent state for the single-file fault VTKHDF writer.
///
/// One instance per `seas::ParaViewOutput`.  All members live only on
/// rank 0 (non-root ranks have empty/null members).  The mesh and
/// per-field GridFunctions are constructed once on the first call to
/// `WriteFaultPackHdf` and reused across subsequent saves; only the
/// per-cell scalar values change between cycles, so a single Mesh and
/// fixed L2-p0 fields are enough.
struct FaultHDFState
{
   // Rank-0 mesh of triangles (one element per fault face, geometry
   // copied from the first gathered pack).
   std::unique_ptr<Mesh>                          mesh;
   std::unique_ptr<L2_FECollection>               fec;
   std::unique_ptr<FiniteElementSpace>            fes;

   // GridFunctions indexed parallel to `field_names`.
   std::vector<std::unique_ptr<GridFunction>>     gfs;
   std::vector<std::string>                       field_names;

   // The data collection — owns the .vtkhdf file handle.
   std::unique_ptr<ParaViewHDFDataCollection>     dc;

   // Phase 2d.3: chunk-filter selector applied to `dc` at Init time.
   // Default Deflate-level-6 reproduces Phase 2b lossless behaviour;
   // ZfpAccuracy switches to LLNL ZFP filter id 32013 (requires
   // MFEM_USE_H5Z_ZFP at build time and HDF5_PLUGIN_PATH at run time).
   ParaViewHDFDataCollection::HDFCompression compression_alg =
       ParaViewHDFDataCollection::HDFCompression::Deflate;
   double compression_param = 6.0;

   bool initialised = false;

   /// Construct the rank-0 mesh + L2-p0 fields from the first gathered
   /// pack.  Subsequent saves reuse this state.
   ///
   /// @param prefix  Output directory prefix (the .vtkhdf will live at
   ///                <prefix>/fault.vtkhdf).  The chunk filter
   ///                (Deflate level OR ZFP tolerance) is read from
   ///                `compression_alg` / `compression_param` on `*this`,
   ///                which the seas wrapper sets via
   ///                `SetFaultHDFCompression` BEFORE the first save.
   void Init(const std::string &prefix,
             const vtu::GatheredFaultPack &g)
   {
      // R-005 (REVIEW.md 2026-04-28): defensive checks at Init entry.
      // GatheredFaultPack invariants are normally enforced inside the
      // gather, but a future caller that constructs the pack by hand
      // could violate them.  Reading field_names[k] / vertices[3i+j]
      // OOB below would otherwise be undefined behaviour.
      MFEM_VERIFY(g.field_arrays.size() == g.field_names.size(),
                  "FaultHDFState::Init: field_arrays.size()="
                  << g.field_arrays.size() << " != field_names.size()="
                  << g.field_names.size());
      MFEM_VERIFY(g.vertices.size() == 3 * g.triangles.size(),
                  "FaultHDFState::Init: vertices.size()="
                  << g.vertices.size()
                  << " expected 3*triangles.size()=3*"
                  << g.triangles.size());
      for (std::size_t k = 0; k < g.field_arrays.size(); ++k)
      {
         MFEM_VERIFY(g.field_arrays[k].size() == g.triangles.size(),
                     "FaultHDFState::Init: field '" << g.field_names[k]
                     << "' has " << g.field_arrays[k].size()
                     << " entries; expected " << g.triangles.size());
      }
      mesh = BuildFaultTriangleMesh(g);

      // L2 p=0 finite element space — one DOF per element.
      fec = std::make_unique<L2_FECollection>(0, mesh->Dimension());
      fes = std::make_unique<FiniteElementSpace>(mesh.get(), fec.get());

      // One GridFunction per gathered field.
      gfs.clear();
      field_names = g.field_names;
      gfs.reserve(g.field_arrays.size());
      for (std::size_t k = 0; k < g.field_arrays.size(); ++k)
      {
         auto gf = std::make_unique<GridFunction>(fes.get());
         *gf = 0.0;
         gfs.push_back(std::move(gf));
      }

      // Data collection.  COMM_SELF so dc.Save() is a serial HDF5 op;
      // rank 0 owns the entire .vtkhdf file.
      // Renamed from "fault_surface" to "fault" per
      // PLAN_split_bulk_solutions_2026-05-12.  On-disk filename becomes
      // <prefix>/fault.vtkhdf; tests that grep for the literal name
      // were updated in the same commit.
      dc = std::make_unique<ParaViewHDFDataCollection>("fault",
                                                        mesh.get());
      dc->SetPrefixPath(prefix);
      dc->SetDataFormat(VTKFormat::BINARY);
      // R-107: route the user-selected chunk filter through ONE call.
      // SetHDFCompression(Deflate, level) syncs `compression` and
      // `compression_level` internally (R-102 fix), so a separate
      // SetCompression / SetCompressionLevel pair is no longer needed
      // and would just be over-written.  ZfpAccuracy mode applies ZFP
      // to FP datasets and uses `compression_level` (kept >= 0 by the
      // R-102 sync for Deflate, or by the EnableCompression fallback
      // for ZfpAccuracy) for integer-dataset deflate.
      dc->SetHDFCompression(compression_alg, compression_param);
      dc->SetHighOrderOutput(false);
      for (std::size_t k = 0; k < gfs.size(); ++k)
      {
         dc->RegisterField(field_names[k], gfs[k].get());
      }
      initialised = true;
   }
};

/// @brief Append one timestep to the rank-0 .vtkhdf file.
///
/// Collective on the parent ParaViewOutput's communicator: all ranks
/// call `GatherFaultPackToRoot` collectively; rank 0 then builds the
/// state on first use and writes the timestep.  Non-root ranks return
/// after the gather.
///
/// The chunk filter (Deflate level OR ZFP tolerance) is read from
/// `state.compression_alg` / `state.compression_param`, which the seas
/// wrapper configures via `SetFaultHDFCompression` BEFORE the first
/// save.  R-107: there is no per-call filter argument because the
/// underlying ParaViewHDFDataCollection is built lazily on first use
/// and locks in its filter at construction; a per-call argument would
/// silently no-op after the first call.
///
/// @param state   Persistent rank-0 state (constructed in-place on
///                first call).  Must outlive the program / data
///                collection lifetime.
/// @param prefix  Output directory.  The .vtkhdf is written to
///                `<prefix>/fault.vtkhdf`.
/// @param cycle   Time step index (passed to dc.SetCycle).
/// @param time    Simulation time (passed to dc.SetTime).
/// @param local   This rank's `LocalFaultPack`; gathered to rank 0.
/// @param rank    MPI rank.
/// @param nranks  Total ranks.
/// @param comm    Communicator (MPI_COMM_NULL on serial builds).
inline void WriteFaultPackHdf(FaultHDFState &state,
                              const std::string &prefix,
                              int cycle, real_t time,
                              const vtu::LocalFaultPack &local,
                              int rank, int nranks
#ifdef MFEM_USE_MPI
                              , MPI_Comm comm = MPI_COMM_NULL
#endif
                             )
{
   // 1. Collective gather to rank 0 (reuse Phase 1 helper).
   vtu::GatheredFaultPack g =
      vtu::GatherFaultPackToRoot(local, rank, nranks
#ifdef MFEM_USE_MPI
                                 , comm
#endif
                                );

   // 2. Non-root ranks finish here.  The .vtkhdf file is rank-0-only.
   if (rank != 0) { return; }

   // 3. First call on rank 0: construct mesh + GFs + dc.
   if (!state.initialised)
   {
      state.Init(prefix, g);
   }
   else
   {
      // Geometry is fixed across cycles; only check that the gathered
      // shape is unchanged (same number of cells / fields).  If a
      // future caller resizes the fault, that's a contract violation.
      const int n_cells_now =
         static_cast<int>(state.mesh->GetNE());
      MFEM_VERIFY(static_cast<int>(g.triangles.size()) == n_cells_now,
                  "WriteFaultPackHdf: gathered cell count changed from "
                  << n_cells_now << " to " << g.triangles.size()
                  << " — fault geometry is expected to be static.");
      MFEM_VERIFY(g.field_names == state.field_names,
                  "WriteFaultPackHdf: gathered field-name list changed "
                  "between cycles.");
   }

   // 4. Update the per-element GF values from the gathered cell data.
   for (std::size_t k = 0; k < state.gfs.size(); ++k)
   {
      MFEM_VERIFY(g.field_arrays[k].size() ==
                  static_cast<std::size_t>(state.fes->GetNE()),
                  "WriteFaultPackHdf: field '" << g.field_names[k]
                  << "' size mismatch with mesh element count.");
      GridFunction &gf = *state.gfs[k];
      for (int e = 0; e < state.fes->GetNE(); ++e)
      {
         // L2-p0: one DOF per element, vdof index = element id.
         gf(e) = g.field_arrays[k][e];
      }
   }

   // 5. Append timestep.
   state.dc->SetCycle(cycle);
   state.dc->SetTime(time);
   state.dc->Save();
}

/// @brief Faster writer of the same <prefix>/fault.vtkhdf as
/// WriteFaultPackHdf, for a fault whose geometry and field list never
/// change.  Plan item 1 of the Nitsche profiling plan (nitsche-profiling,
/// Chunhui Zhao 2026-10-17).  Per frame it saves:
///   - the gather: the geometry crosses to rank 0 once; later frames send
///     only the field values, all fields in one MPI_Gatherv;
///   - the per-cell GridFunction evaluation of SaveGridFunction: values go
///     straight to VTKHDF::SavePointData (each cell value at the cell's 3
///     points, exactly what SaveGridFunction writes for an L2 p = 0 field);
///   - the re-compression of every dataset's last, partly filled 0.5 MB
///     chunk: new datasets are chunked to their first append
///     (VTKHDF::SetChunkToFirstAppend);
///   - the constant fields (constant[k] != 0): written at the first frame
///     only; later steps keep their Steps/PointDataOffsets entry at the
///     first frame's data, which VTKHDF readers follow.
/// The values, point order, field names and step layout are otherwise those
/// of WriteFaultPackHdf with the same filter.  Collective; the file lives on
/// rank 0.
class FaultHDFStream
{
public:
   /// Chunk filter, as FaultHDFState::compression_alg / compression_param.
   /// Call before the first Write.
   void SetCompression(ParaViewHDFDataCollection::HDFCompression alg, double param)
   {
      alg_ = alg;
      param_ = param;
   }

   void Write(const std::string &prefix, real_t time, const vtu::LocalFaultPack &local,
              const std::vector<char> &constant, int rank, int nranks
#ifdef MFEM_USE_MPI
              , MPI_Comm comm = MPI_COMM_NULL
#endif
             )
   {
      const double t0 = Now_();
#ifdef MFEM_USE_MPI
      const bool par = nranks > 1 && comm != MPI_COMM_NULL;
#else
      const bool par = false;
#endif
      const std::size_t nf = local.field_arrays.size();
      if (frames_ == 0)
      {
         MFEM_VERIFY(constant.size() == nf, "FaultHDFStream: " << constant.size()
                     << " constant flags for " << nf << " fields");
         constant_ = constant;
         names_ = local.field_names;
         ncell_ = static_cast<int>(local.triangles.size());
         for (std::size_t k = 0; k < nf; k++) { if (!constant_[k]) { var_.push_back(static_cast<int>(k)); } }
         vtu::GatheredFaultPack g = vtu::GatherFaultPackToRoot(local, rank, nranks
#ifdef MFEM_USE_MPI
                                                                , comm
#endif
                                                               );
#ifdef MFEM_USE_MPI
         if (par)
         {
            counts_.assign(rank == 0 ? nranks : 0, 0);
            MPI_Gather(&ncell_, 1, MPI_INT, rank == 0 ? counts_.data() : nullptr, 1, MPI_INT, 0, comm);
         }
#endif
         gather_sec_ += Now_() - t0;
         if (rank == 0) { Open_(prefix, g); }
         if (rank == 0) { Save_(time, g.field_arrays); }
      }
      else
      {
         MFEM_VERIFY(static_cast<int>(local.triangles.size()) == ncell_ && local.field_names == names_,
                     "FaultHDFStream: the fault geometry or field list changed after the first frame");
         if (!par)
         {
            gather_sec_ += Now_() - t0;
            if (rank == 0) { Save_(time, local.field_arrays); }
         }
#ifdef MFEM_USE_MPI
         else
         {
            const int nv = static_cast<int>(var_.size());
            send_.resize(static_cast<std::size_t>(nv) * ncell_);
            for (int v = 0; v < nv; v++)
            {
               std::copy(local.field_arrays[var_[v]].begin(), local.field_arrays[var_[v]].end(),
                         send_.begin() + static_cast<std::size_t>(v) * ncell_);
            }
            std::vector<int> cnt, dsp;
            int total = 0;
            if (rank == 0)
            {
               cnt.resize(nranks);
               dsp.resize(nranks);
               for (int r = 0; r < nranks; r++)
               {
                  cnt[r] = nv * counts_[r];
                  dsp[r] = total;
                  total += cnt[r];
               }
               recv_.resize(total);
            }
            MPI_Gatherv(send_.data(), static_cast<int>(send_.size()), MPI_DOUBLE,
                        rank == 0 ? recv_.data() : nullptr, rank == 0 ? cnt.data() : nullptr,
                        rank == 0 ? dsp.data() : nullptr, MPI_DOUBLE, 0, comm);
            if (rank == 0)
            {
               // rank r's block: nv fields x counts_[r] cells, field-major
               arrays_.assign(nf, std::vector<double>());
               for (int v = 0; v < nv; v++) { arrays_[var_[v]].resize(ncell_total_); }
               int cell0 = 0;
               for (int r = 0; r < nranks; r++)
               {
                  for (int v = 0; v < nv; v++)
                  {
                     const double *b = recv_.data() + dsp[r] + static_cast<std::size_t>(v) * counts_[r];
                     std::copy(b, b + counts_[r], arrays_[var_[v]].begin() + cell0);
                  }
                  cell0 += counts_[r];
               }
            }
            gather_sec_ += Now_() - t0;
            if (rank == 0) { Save_(time, arrays_); }
         }
#endif
      }
      frames_++;
   }

   int Frames() const { return frames_; }
   /// Seconds in the gather (all ranks) and in the HDF5 write (rank 0).
   double GatherSeconds() const { return gather_sec_; }
   double WriteSeconds() const { return write_sec_; }

private:
   static double Now_()
   {
      return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
   }

   /// Rank 0: the mesh and the file, with the filter mapping of
   /// ParaViewHDFDataCollection::TSave.
   void Open_(const std::string &prefix, const vtu::GatheredFaultPack &g)
   {
      ncell_total_ = static_cast<int>(g.triangles.size());
      mesh_ = BuildFaultTriangleMesh(g);
      std::string dir = prefix;
      if (!dir.empty() && dir.back() != '/') { dir += '/'; }
      if (!dir.empty()) { ::mkdir(dir.c_str(), 0775); }
      file_ = std::make_unique<VTKHDF>(dir + "fault.vtkhdf");
      const bool compression = (alg_ == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
                               || param_ >= 0.0;
      const int level = (alg_ == ParaViewHDFDataCollection::HDFCompression::Deflate && param_ >= 0.0)
                        ? static_cast<int>(param_) : -1;
      if (alg_ == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
      {
         file_->SetZfpAccuracy(param_);
      }
      else { file_->SetCompressionAlgorithm(VTKHDF::CompressionAlgorithm::Deflate); }
      if (compression) { file_->EnableCompression(level >= 0 ? level : 6); }
      else { file_->DisableCompression(); }
      file_->SetChunkToFirstAppend(true);
      // the data collection writes its fields in name order (std::map)
      order_.resize(names_.size());
      for (std::size_t k = 0; k < order_.size(); k++) { order_[k] = static_cast<int>(k); }
      std::sort(order_.begin(), order_.end(), [&](int a, int b) { return names_[a] < names_[b]; });
      // points per cell of SaveMesh(mesh, high_order = false, ref = 1)
      ppc_ = GlobGeometryRefiner.Refine(Geometry::TRIANGLE, 1, 1)->RefPts.GetNPoints();
   }

   /// Rank 0: one step.  arrays[k] = per-cell values of field k (constant
   /// fields are read at the first frame only).
   void Save_(real_t time, const std::vector<std::vector<double>> &arrays)
   {
      const double t0 = Now_();
      file_->SaveMesh<double>(*mesh_, /*high_order=*/false, /*ref=*/1);
      for (int k : order_)
      {
         if (frames_ > 0 && constant_[k]) { continue; }
         const std::vector<double> &a = arrays[k];
         MFEM_VERIFY(static_cast<int>(a.size()) == ncell_total_, "FaultHDFStream: field '"
                     << names_[k] << "' has " << a.size() << " cells, the mesh " << ncell_total_);
         pts_.resize(static_cast<std::size_t>(ppc_) * a.size());
         for (std::size_t c = 0; c < a.size(); c++)
         {
            for (int i = 0; i < ppc_; i++) { pts_[ppc_ * c + i] = a[c]; }
         }
         file_->SavePointData<double>(names_[k], pts_, 1);
      }
      file_->UpdateSteps(time);
      file_->Flush();
      write_sec_ += Now_() - t0;
   }

   ParaViewHDFDataCollection::HDFCompression alg_ =
      ParaViewHDFDataCollection::HDFCompression::Deflate;
   double param_ = 6.0;
   int frames_ = 0, ncell_ = 0, ncell_total_ = 0, ppc_ = 0;
   std::vector<std::string> names_;
   std::vector<char> constant_;
   std::vector<int> var_, order_, counts_;
   std::vector<double> send_, recv_, pts_;
   std::vector<std::vector<double>> arrays_;
   std::unique_ptr<Mesh> mesh_;
   std::unique_ptr<VTKHDF> file_;
   double gather_sec_ = 0.0, write_sec_ = 0.0;
};

} // namespace vtkhdf
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_VTKHDF_WRITER_HPP
