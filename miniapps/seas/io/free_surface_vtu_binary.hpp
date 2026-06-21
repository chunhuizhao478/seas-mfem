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

#ifndef MFEM_SEAS_FREE_SURFACE_VTU_BINARY_HPP
#define MFEM_SEAS_FREE_SURFACE_VTU_BINARY_HPP

// -----------------------------------------------------------------------------
// Single-file free-surface VTU writer.
//
// Mirrors io/fault_vtu_binary.hpp (the fault's `fault_surface_c<cycle>.vtu`
// writer): each output cycle, every rank packs its local free-surface geometry,
// a variable-length MPI gather concatenates to rank 0, and rank 0 emits ONE
// binary VTU (raw appended, UInt64 headers).  This replaces MFEM's per-rank
// `ParaViewDataCollection` proc*.vtu explosion (N_ranks files per cycle) with a
// single consolidated file, exactly like the fault output.
//
// DIFFERENCE FROM THE FAULT WRITER: the fault carries per-triangle CellData
// scalars (slip/traction, averaged).  The free surface instead carries the
// VELOCITY as a 3-component POINT-data vector (the ground-motion field), so
// ParaView shows it as a true vector glyph.  Geometry is triangle-only (each
// triangle owns 3 unique, unshared vertices: vertices.size()==3*ntri) — the
// caller triangulates quad faces into two triangles before packing, so a hex
// free surface and a tet free surface share this writer.  The unshared-vertex
// layout represents the DG-discontinuous L2 velocity exactly (no inter-element
// averaging) and makes the gather index remap a constant offset add.
// -----------------------------------------------------------------------------

#include "mfem.hpp"
#include "mesh/vtk.hpp"          // VTKFormat, VTKByteOrder, VTKGeometry::TRIANGLE

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mfem
{
namespace seas
{
namespace fsvtu
{

/// One rank's free-surface contribution for a single write cycle.
///
/// Layout invariant: `vertices.size() == velocity.size() == 3 * ntri()` and
/// triangle `i` owns vertices {3*i, 3*i+1, 3*i+2} (no inter-cell sharing) —
/// identical to the fault `LocalFaultPack` so the gather remap is a constant
/// offset add.  `velocity[k]` is the 3-component velocity at vertex `k`
/// (POINT data).  `rank_cell[i]` is the owning MPI rank of triangle `i`
/// (CellData; one value per triangle).
struct LocalSurfacePack
{
   std::vector<std::array<double, 3>> vertices;   // 3*ntri
   std::vector<std::array<double, 3>> velocity;   // 3*ntri (per-vertex)
   std::vector<double>                rank_cell;  // ntri
   int ntri() const { return static_cast<int>(rank_cell.size()); }
};

/// Result of GatherSurfacePackToRoot.  Populated on rank 0; empty elsewhere.
struct GatheredSurfacePack
{
   std::vector<std::array<double, 3>> vertices;
   std::vector<std::array<double, 3>> velocity;
   std::vector<double>                rank_cell;
   std::vector<std::array<int, 3>>    triangles;
   int ntri() const { return static_cast<int>(rank_cell.size()); }
};

// =============================================================================
// MPI gather (collective on `comm`)
// =============================================================================

/// Concatenate per-rank `LocalSurfacePack`s onto rank 0.  All ranks call
/// collectively.  Non-root ranks get an empty struct.  Triangle indices are
/// rebuilt on root as {3c, 3c+1, 3c+2} over the gathered global vertex array.
inline GatheredSurfacePack
GatherSurfacePackToRoot(const LocalSurfacePack &local, int rank, int nranks
#ifdef MFEM_USE_MPI
                        , MPI_Comm comm = MPI_COMM_NULL
#endif
                       )
{
   GatheredSurfacePack out;
   MFEM_VERIFY(local.vertices.size() == 3u * local.rank_cell.size() &&
               local.velocity.size() == 3u * local.rank_cell.size(),
               "GatherSurfacePackToRoot: pack invariant violated"
               " (vertices/velocity must be 3*ntri)");

#ifdef MFEM_USE_MPI
   if (nranks <= 1 || comm == MPI_COMM_NULL)
#else
   (void)rank; (void)nranks;
#endif
   {
      out.vertices  = local.vertices;
      out.velocity  = local.velocity;
      out.rank_cell = local.rank_cell;
      const int nt = out.ntri();
      out.triangles.resize(nt);
      for (int c = 0; c < nt; ++c) { out.triangles[c] = {3*c, 3*c+1, 3*c+2}; }
      return out;
   }

#ifdef MFEM_USE_MPI
   const int local_nt = local.ntri();

   std::vector<int> counts(nranks, 0);
   MPI_Allgather(&local_nt, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

   int total = 0;
   std::vector<int> cnt_v9(nranks, 0), displ_v9(nranks, 0),
                    cnt_cell(nranks, 0), displ_cell(nranks, 0);
   if (rank == 0)
   {
      int accv = 0, accc = 0;
      for (int r = 0; r < nranks; ++r)
      {
         cnt_v9[r]     = counts[r] * 9;   displ_v9[r]   = accv; accv += cnt_v9[r];
         cnt_cell[r]   = counts[r];       displ_cell[r] = accc; accc += counts[r];
      }
      total = accc;
   }

   // Flatten geometry + velocity, 9 doubles per triangle (3 verts * 3 comps).
   std::vector<double> lg(static_cast<std::size_t>(local_nt) * 9);
   std::vector<double> lv(static_cast<std::size_t>(local_nt) * 9);
   for (int i = 0; i < local_nt; ++i)
   {
      for (int j = 0; j < 3; ++j)
      {
         for (int d = 0; d < 3; ++d)
         {
            lg[9*i + 3*j + d] = local.vertices[3*i + j][d];
            lv[9*i + 3*j + d] = local.velocity[3*i + j][d];
         }
      }
   }

   std::vector<double> rg, rv, rrank;
   if (rank == 0)
   {
      rg.resize(static_cast<std::size_t>(total) * 9);
      rv.resize(static_cast<std::size_t>(total) * 9);
      rrank.resize(total);
   }
   MPI_Gatherv(lg.data(), local_nt * 9, MPI_DOUBLE,
               rank == 0 ? rg.data() : nullptr,
               rank == 0 ? cnt_v9.data() : nullptr,
               rank == 0 ? displ_v9.data() : nullptr, MPI_DOUBLE, 0, comm);
   MPI_Gatherv(lv.data(), local_nt * 9, MPI_DOUBLE,
               rank == 0 ? rv.data() : nullptr,
               rank == 0 ? cnt_v9.data() : nullptr,
               rank == 0 ? displ_v9.data() : nullptr, MPI_DOUBLE, 0, comm);
   MPI_Gatherv(const_cast<double *>(local.rank_cell.data()), local_nt, MPI_DOUBLE,
               rank == 0 ? rrank.data() : nullptr,
               rank == 0 ? cnt_cell.data() : nullptr,
               rank == 0 ? displ_cell.data() : nullptr, MPI_DOUBLE, 0, comm);

   if (rank == 0)
   {
      out.vertices.resize(static_cast<std::size_t>(total) * 3);
      out.velocity.resize(static_cast<std::size_t>(total) * 3);
      out.rank_cell = std::move(rrank);
      out.triangles.resize(total);
      for (int c = 0; c < total; ++c)
      {
         for (int j = 0; j < 3; ++j)
         {
            out.vertices[3*c + j] = { rg[9*c + 3*j + 0], rg[9*c + 3*j + 1],
                                      rg[9*c + 3*j + 2] };
            out.velocity[3*c + j] = { rv[9*c + 3*j + 0], rv[9*c + 3*j + 1],
                                      rv[9*c + 3*j + 2] };
         }
         out.triangles[c] = {3*c, 3*c + 1, 3*c + 2};
      }
   }
   return out;
#endif // MFEM_USE_MPI
}

// =============================================================================
// Pack -> VTU emit (rank 0 only)
// =============================================================================

/// Emit a single VTU at @a path: triangle UnstructuredGrid with a 3-component
/// POINT-data `velocity` vector and a CellData `mpi_rank` scalar.  BINARY uses
/// the raw-appended UInt64 idiom (byte-identical encoding to WriteFaultPackVTU);
/// ASCII is a regex-parseable fallback for unit tests.  A zero-triangle pack
/// still emits a valid (empty) VTU so the PVD entry stays well-formed.
inline void WriteSurfacePackVTU(const std::string &path,
                                const GatheredSurfacePack &g,
                                VTKFormat format = VTKFormat::BINARY)
{
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
   static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
                 "WriteSurfacePackVTU requires a little-endian host");
#endif
   const int npts   = static_cast<int>(g.vertices.size());
   const int ncells = static_cast<int>(g.triangles.size());
   MFEM_VERIFY(static_cast<int>(g.velocity.size()) == npts,
               "WriteSurfacePackVTU: velocity.size() != points");
   MFEM_VERIFY(static_cast<int>(g.rank_cell.size()) == ncells,
               "WriteSurfacePackVTU: rank_cell.size() != cells");

   std::ofstream vtu(path, std::ios::binary);
   MFEM_VERIFY(vtu.is_open(), "WriteSurfacePackVTU: failed to open " << path);

   if (format == VTKFormat::ASCII)
   {
      vtu << std::setprecision(10);
      vtu << "<?xml version=\"1.0\"?>\n";
      vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      vtu << "<UnstructuredGrid>\n";
      vtu << "<Piece NumberOfPoints=\"" << npts
          << "\" NumberOfCells=\"" << ncells << "\">\n";
      vtu << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\""
             " format=\"ascii\">\n";
      for (const auto &v : g.vertices)
      { vtu << v[0] << " " << v[1] << " " << v[2] << "\n"; }
      vtu << "</DataArray></Points>\n<Cells>\n"
             "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
      for (const auto &t : g.triangles)
      { vtu << t[0] << " " << t[1] << " " << t[2] << "\n"; }
      vtu << "</DataArray>\n"
             "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; ++i) { vtu << (i + 1) * 3 << "\n"; }
      vtu << "</DataArray>\n"
             "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; ++i)
      { vtu << static_cast<int>(VTKGeometry::TRIANGLE) << "\n"; }
      vtu << "</DataArray>\n</Cells>\n";
      vtu << "<PointData Vectors=\"velocity\"><DataArray type=\"Float64\""
             " Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      for (const auto &v : g.velocity)
      { vtu << v[0] << " " << v[1] << " " << v[2] << "\n"; }
      vtu << "</DataArray></PointData>\n";
      vtu << "<CellData><DataArray type=\"Float64\" Name=\"mpi_rank\""
             " format=\"ascii\">\n";
      for (double r : g.rank_cell) { vtu << r << "\n"; }
      vtu << "</DataArray></CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
      return;
   }

   // Binary raw-appended (matches WriteFaultPackVTU).  Block order:
   //   0 Points(Float64,3c) 1 connectivity(Int32) 2 offsets(Int32)
   //   3 types(UInt8) 4 velocity(Float64,3c POINT) 5 mpi_rank(Float64 CELL)
   struct AppendedBlock { std::vector<unsigned char> bytes; };
   std::vector<AppendedBlock> blocks;
   auto add_block = [&](const void *src, std::size_t n)
   {
      AppendedBlock b; b.bytes.resize(n);
      if (n > 0) { std::memcpy(b.bytes.data(), src, n); }
      blocks.push_back(std::move(b));
   };

   {
      std::vector<double> flat(static_cast<std::size_t>(npts) * 3);
      for (int i = 0; i < npts; ++i)
      { flat[3*i] = g.vertices[i][0]; flat[3*i+1] = g.vertices[i][1];
        flat[3*i+2] = g.vertices[i][2]; }
      add_block(flat.data(), flat.size() * sizeof(double));
   }
   {
      std::vector<int32_t> flat(static_cast<std::size_t>(ncells) * 3);
      for (int c = 0; c < ncells; ++c)
      { flat[3*c] = g.triangles[c][0]; flat[3*c+1] = g.triangles[c][1];
        flat[3*c+2] = g.triangles[c][2]; }
      add_block(flat.data(), flat.size() * sizeof(int32_t));
   }
   {
      std::vector<int32_t> flat(ncells);
      for (int c = 0; c < ncells; ++c) { flat[c] = (c + 1) * 3; }
      add_block(flat.data(), flat.size() * sizeof(int32_t));
   }
   {
      std::vector<uint8_t> flat(ncells,
                                static_cast<uint8_t>(VTKGeometry::TRIANGLE));
      add_block(flat.data(), flat.size() * sizeof(uint8_t));
   }
   {
      std::vector<double> flat(static_cast<std::size_t>(npts) * 3);
      for (int i = 0; i < npts; ++i)
      { flat[3*i] = g.velocity[i][0]; flat[3*i+1] = g.velocity[i][1];
        flat[3*i+2] = g.velocity[i][2]; }
      add_block(flat.data(), flat.size() * sizeof(double));
   }
   add_block(g.rank_cell.data(), g.rank_cell.size() * sizeof(double));

   std::vector<uint64_t> off(blocks.size(), 0);
   {
      uint64_t pos = 0;
      for (std::size_t k = 0; k < blocks.size(); ++k)
      { off[k] = pos; pos += sizeof(uint64_t) + blocks[k].bytes.size(); }
   }

   vtu << "<?xml version=\"1.0\"?>\n";
   vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\""
       << " byte_order=\"" << VTKByteOrder() << "\""
       << " header_type=\"UInt64\">\n";
   vtu << "<UnstructuredGrid>\n";
   vtu << "<Piece NumberOfPoints=\"" << npts
       << "\" NumberOfCells=\"" << ncells << "\">\n";
   vtu << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\""
       << " format=\"appended\" offset=\"" << off[0] << "\"/></Points>\n";
   vtu << "<Cells>\n"
       << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\""
          " offset=\"" << off[1] << "\"/>\n"
       << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\""
          " offset=\"" << off[2] << "\"/>\n"
       << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\""
          " offset=\"" << off[3] << "\"/>\n</Cells>\n";
   vtu << "<PointData Vectors=\"velocity\">"
       << "<DataArray type=\"Float64\" Name=\"velocity\""
          " NumberOfComponents=\"3\" format=\"appended\""
          " offset=\"" << off[4] << "\"/></PointData>\n";
   vtu << "<CellData>"
       << "<DataArray type=\"Float64\" Name=\"mpi_rank\" format=\"appended\""
          " offset=\"" << off[5] << "\"/></CellData>\n";
   vtu << "</Piece>\n</UnstructuredGrid>\n";

   vtu << "<AppendedData encoding=\"raw\">\n_";
   for (const auto &b : blocks)
   {
      const uint64_t len = b.bytes.size();
      vtu.write(reinterpret_cast<const char *>(&len), sizeof(len));
      if (len > 0)
      {
         vtu.write(reinterpret_cast<const char *>(b.bytes.data()),
                   static_cast<std::streamsize>(len));
      }
   }
   vtu << "\n</AppendedData>\n</VTKFile>\n";
}

} // namespace fsvtu
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FREE_SURFACE_VTU_BINARY_HPP
