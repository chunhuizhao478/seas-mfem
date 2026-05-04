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

#include <memory>
#include <string>
#include <vector>

#ifndef MFEM_USE_HDF5
#error "fault_vtkhdf_writer.hpp requires MFEM_USE_HDF5=YES.  Either enable HDF5 in config/config.mk and rebuild MFEM, or use Phase 1's binary VTU path (the default in this build)."
#endif

#include "fem/datacollection.hpp"  // ParaViewHDFDataCollection

namespace mfem
{
namespace seas
{
namespace vtkhdf
{

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

   bool initialised = false;

   /// Construct the rank-0 mesh + L2-p0 fields from the first gathered
   /// pack.  Subsequent saves reuse this state.
   ///
   /// @param prefix       Output directory prefix (the .vtkhdf will
   ///                     live at <prefix>/fault_surface.vtkhdf).
   /// @param compression  zlib compression level (0..9); 0 disables.
   void Init(const std::string &prefix,
             const vtu::GatheredFaultPack &g,
             int compression_level)
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
      const int ncells = static_cast<int>(g.triangles.size());
      const int npts   = static_cast<int>(g.vertices.size());

      // Build a serial 2D-in-3D Mesh of triangles.  spaceDim=3 because
      // the fault sits in a 3D embedding space.
      mesh = std::make_unique<Mesh>(/*Dim=*/2, /*NVert=*/npts,
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
      mesh->FinalizeTriMesh(/*generate_edges=*/0, /*refine=*/0,
                            /*fix_orientation=*/true);

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
      dc = std::make_unique<ParaViewHDFDataCollection>("fault_surface",
                                                        mesh.get());
      dc->SetPrefixPath(prefix);
      dc->SetDataFormat(VTKFormat::BINARY);
      // R-001 (REVIEW.md 2026-04-28): plan §Phase 2 §4 mandates HDF5
      // compression always ON.  HDF5's internal compression is
      // independent of MFEM_USE_ZLIB (per fem/datacollection.hpp:673).
      // A `compression_level` of 0 or negative falls back to the plan-
      // mandated level 3 ("fast"); a positive value is honoured.
      dc->SetCompression(true);
      dc->SetCompressionLevel(compression_level > 0 ? compression_level : 3);
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
/// @param state       Persistent rank-0 state (constructed in-place on
///                    first call).  Must outlive the program / data
///                    collection lifetime.
/// @param prefix      Output directory.  The .vtkhdf is written to
///                    `<prefix>/fault_surface.vtkhdf`.
/// @param cycle       Time step index (passed to dc.SetCycle).
/// @param time        Simulation time (passed to dc.SetTime).
/// @param compression zlib compression level (0..9); ignored after
///                    first call (locked in at Init time).
/// @param local       This rank's `LocalFaultPack`; gathered to rank 0.
/// @param rank        MPI rank.
/// @param nranks      Total ranks.
/// @param comm        Communicator (MPI_COMM_NULL on serial builds).
inline void WriteFaultPackHdf(FaultHDFState &state,
                              const std::string &prefix,
                              int cycle, real_t time,
                              int compression_level,
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
      state.Init(prefix, g, compression_level);
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

} // namespace vtkhdf
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_VTKHDF_WRITER_HPP
