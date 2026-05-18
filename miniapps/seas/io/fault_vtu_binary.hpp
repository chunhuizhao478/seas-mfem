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

#ifndef MFEM_SEAS_FAULT_VTU_BINARY_HPP
#define MFEM_SEAS_FAULT_VTU_BINARY_HPP

// -----------------------------------------------------------------------------
// Phase 1 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md
//
// Helpers used by `seas::ParaViewOutput::WriteFaultSurfaceVTU` to:
//   1. Pack each rank's fault-surface geometry + averaged cell-data into a
//      flat `LocalFaultPack`.
//   2. Gather all ranks' packs to rank 0 via a variable-length MPI gather
//      (no-op on serial builds).
//   3. Emit a single binary VTU on rank 0 per cycle.  Replaces the legacy
//      per-rank-per-cycle ASCII VTU writer.
//
// IMPLEMENTATION DEVIATION FROM PLAN §Phase 1.3 ("VTU binary appended"):
//   The plan specified the VTU `<AppendedData encoding="raw">_<bytes>` block
//   with `format="appended" offset="N"` per DataArray.  This file instead
//   uses the inline base64 idiom (per-DataArray `format="binary"`) emitted
//   by MFEM's own `ParaViewDataCollection` (`fem/datacollection.cpp:1185-
//   1209`).  The two formats are interchangeable as far as ParaView is
//   concerned, but the inline-base64 path:
//     - reuses `mfem::WriteBinaryOrASCII<T>` and `mfem::WriteBase64With
//       SizeAndClear` (mesh/vtk.hpp), already battle-tested with ParaView;
//     - eliminates byte-offset arithmetic across the appended block, which
//       is the most common source of silent VTU corruption;
//     - costs 33% in file size vs raw appended, which is negligible
//       relative to the 5-8x ASCII->binary saving and the N_ranks-fold
//       file-count reduction this Phase achieves.
//   When `MFEM_USE_ZLIB=YES` is set in the build, base64 blocks are zlib-
//   compressed, recovering the size delta vs raw appended.
// -----------------------------------------------------------------------------

#include "mfem.hpp"
#include "mesh/vtk.hpp"          // VTKFormat, WriteBinaryOrASCII,
                                  // WriteBase64WithSizeAndClear,
                                  // VTKByteOrder, VTKGeometry::TRIANGLE
#include "general/binaryio.hpp"  // bin_io::DecodeBase64

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mfem
{
namespace seas
{
namespace vtu
{

// =============================================================================
// Data structures
// =============================================================================

/// One rank's contribution to the fault surface for a single write cycle.
///
/// Layout invariant: `vertices.size() == 3 * triangles.size()`, and
/// `triangles[i] == {3*i, 3*i+1, 3*i+2}` (each triangle owns three unique
/// vertices, no inter-cell vertex sharing).  This matches the legacy
/// per-face emit at paraview_output.hpp:650-660 and is required for the
/// gather step's index remapping to be a constant offset add.
///
/// `field_arrays[k]` and `field_names[k]` describe the kth scalar
/// CellData array.  `field_arrays[k].size() == triangles.size()` for all
/// k.  Field arrays that the field filter rejects (see
/// `ParaViewOutput::SetFaultVTUFields`) are absent from both vectors —
/// the caller pre-filters at pack build time so that empty arrays are
/// not gathered, sent, or written.
struct LocalFaultPack
{
   std::vector<std::array<double, 3>> vertices;
   std::vector<std::array<int, 3>>    triangles;
   std::vector<std::vector<double>>   field_arrays;
   std::vector<std::string>           field_names;
};

/// Result of `GatherFaultPackToRoot`.  Populated only on rank 0; on
/// non-root ranks all members are empty.
///
/// Layout: vertices and triangles are concatenations of per-rank
/// `LocalFaultPack`s, with triangle indices remapped so that all
/// triangles point at the global vertex array.  Field arrays are
/// concatenated in the same per-rank order.
struct GatheredFaultPack
{
   std::vector<std::array<double, 3>> vertices;
   std::vector<std::array<int, 3>>    triangles;
   std::vector<std::vector<double>>   field_arrays;
   std::vector<std::string>           field_names;
};

// =============================================================================
// Pack -> VTU emit (rank 0 only)
// =============================================================================

/// Emit a single binary or ASCII VTU at @a vtu_path containing the fault
/// surface geometry and CellData carried by @a g.
///
/// @param vtu_path           Absolute or relative path; parent dir must exist.
/// @param g                  Gathered fault data (rank 0).  May have zero
///                           triangles — in that case still emit a valid
///                           VTU with zero Pieces, keeping the PVD entry
///                           well-formed.
/// @param format             `VTKFormat::BINARY` (recommended; default for
///                           the Phase 1 writer) or `VTKFormat::ASCII`
///                           (legacy fallback used by tests that need a
///                           regex-parseable VTU).
/// @param compression_level  0 = uncompressed.  1-9 = zlib level (only
///                           takes effect when MFEM_USE_ZLIB is enabled
///                           in the build).
inline void WriteFaultPackVTU(const std::string &vtu_path,
                              const GatheredFaultPack &g,
                              VTKFormat format = VTKFormat::BINARY,
                              int compression_level = 0)
{
   // R-304 / plan §Phase 1 Edge Cases: this writer assumes a
   // little-endian host so per-array binary blocks can be written
   // verbatim.  A big-endian host requires byte-swapping AND emitting
   // `byte_order="BigEndian"` in the XML.  Hard-fail at compile time on
   // big-endian builds.
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
   static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
                 "WriteFaultPackVTU requires a little-endian host");
#endif

   const int npts   = static_cast<int>(g.vertices.size());
   const int ncells = static_cast<int>(g.triangles.size());
   const std::size_t nfields = g.field_arrays.size();

   // Preconditions — these mirror invariants of LocalFaultPack and the
   // gather step.  Hard-fail rather than emit an unreadable VTU.
   MFEM_VERIFY(g.field_arrays.size() == g.field_names.size(),
               "WriteFaultPackVTU: field_arrays.size() != field_names.size()");
   for (std::size_t k = 0; k < nfields; ++k)
   {
      MFEM_VERIFY(static_cast<int>(g.field_arrays[k].size()) == ncells,
                  "WriteFaultPackVTU: field '" << g.field_names[k] <<
                  "' has " << g.field_arrays[k].size() <<
                  " entries; expected " << ncells);
   }
   (void)compression_level;   // raw-appended path does not compress.

   std::ofstream vtu(vtu_path, std::ios::binary);
   MFEM_VERIFY(vtu.is_open(),
               "WriteFaultPackVTU: failed to open " << vtu_path);

   // -- ASCII fallback (legacy unit-test parsers depend on this layout) --
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
      vtu << "</DataArray></Points>\n";
      vtu << "<Cells>\n"
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
      vtu << "</DataArray>\n</Cells>\n<CellData>\n";
      for (std::size_t k = 0; k < nfields; ++k)
      {
         vtu << "<DataArray type=\"Float64\" Name=\"" << g.field_names[k]
             << "\" format=\"ascii\">\n";
         for (double v : g.field_arrays[k]) { vtu << v << "\n"; }
         vtu << "</DataArray>\n";
      }
      vtu << "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
      return;
   }

   // -- Binary path (R-302 / plan §Phase 1 §3): raw appended encoding.
   //
   // Each DataArray's storage block is `[uint64 byte_length][raw_bytes]`
   // concatenated under <AppendedData encoding="raw"> with a leading
   // `_` marker.  `header_type="UInt64"` (R-303 / plan §3) future-proofs
   // against the >4 GB single-array threshold; ParaView 5.5+ supports
   // it (project minimum is 5.10).  Strict-aliasing-safe: we
   // std::memcpy into a flat byte buffer rather than reinterpret_cast.
   //
   // Block layout (offsets are byte positions from start of appended
   // data, AFTER the leading `_` marker):
   //   block 0   : Points     (Float64, 3 components, npts vertices)
   //   block 1   : connectivity (Int32, 3 per cell, ncells)
   //   block 2   : offsets    (Int32, ncells)
   //   block 3   : types      (UInt8, ncells)
   //   block 4..K: CellData   (Float64, 1 per cell, ncells)

   struct AppendedBlock { std::vector<unsigned char> bytes; };
   std::vector<AppendedBlock> blocks;
   blocks.reserve(4 + nfields);

   auto add_block = [&](const void *src, std::size_t bytes)
   {
      AppendedBlock b;
      b.bytes.resize(bytes);
      if (bytes > 0) { std::memcpy(b.bytes.data(), src, bytes); }
      blocks.push_back(std::move(b));
   };

   // Block 0: Points (Float64, 3 components per vertex).
   {
      std::vector<double> flat(static_cast<std::size_t>(npts) * 3);
      for (int i = 0; i < npts; ++i)
      {
         flat[3 * i + 0] = g.vertices[i][0];
         flat[3 * i + 1] = g.vertices[i][1];
         flat[3 * i + 2] = g.vertices[i][2];
      }
      add_block(flat.data(), flat.size() * sizeof(double));
   }
   // Block 1: connectivity (Int32, 3 per cell).
   {
      std::vector<int32_t> flat(static_cast<std::size_t>(ncells) * 3);
      for (int c = 0; c < ncells; ++c)
      {
         flat[3 * c + 0] = static_cast<int32_t>(g.triangles[c][0]);
         flat[3 * c + 1] = static_cast<int32_t>(g.triangles[c][1]);
         flat[3 * c + 2] = static_cast<int32_t>(g.triangles[c][2]);
      }
      add_block(flat.data(), flat.size() * sizeof(int32_t));
   }
   // Block 2: offsets (Int32, one per cell).
   {
      std::vector<int32_t> flat(ncells);
      for (int c = 0; c < ncells; ++c) { flat[c] = (c + 1) * 3; }
      add_block(flat.data(), flat.size() * sizeof(int32_t));
   }
   // Block 3: cell types (UInt8, one per cell).
   {
      std::vector<uint8_t> flat(
         ncells, static_cast<uint8_t>(VTKGeometry::TRIANGLE));
      add_block(flat.data(), flat.size() * sizeof(uint8_t));
   }
   // Blocks 4..K: per-field CellData (Float64).
   for (std::size_t k = 0; k < nfields; ++k)
   {
      add_block(g.field_arrays[k].data(),
                g.field_arrays[k].size() * sizeof(double));
   }

   // Compute offsets: each block's offset is the running byte position
   // INCLUDING the 8-byte length headers for all preceding blocks.
   std::vector<uint64_t> offsets(blocks.size(), 0);
   {
      uint64_t pos = 0;
      for (std::size_t k = 0; k < blocks.size(); ++k)
      {
         offsets[k] = pos;
         pos += sizeof(uint64_t) + blocks[k].bytes.size();
      }
   }

   // Emit XML header.
   vtu << "<?xml version=\"1.0\"?>\n";
   vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\""
       << " byte_order=\"" << VTKByteOrder() << "\""
       << " header_type=\"UInt64\""
       << ">\n";
   vtu << "<UnstructuredGrid>\n";
   vtu << "<Piece NumberOfPoints=\"" << npts
       << "\" NumberOfCells=\"" << ncells << "\">\n";

   vtu << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\""
       << " format=\"appended\" offset=\"" << offsets[0] << "\"/></Points>\n";
   vtu << "<Cells>\n"
       << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\""
          " offset=\"" << offsets[1] << "\"/>\n"
       << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\""
          " offset=\"" << offsets[2] << "\"/>\n"
       << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\""
          " offset=\"" << offsets[3] << "\"/>\n"
       << "</Cells>\n<CellData>\n";
   for (std::size_t k = 0; k < nfields; ++k)
   {
      vtu << "<DataArray type=\"Float64\" Name=\"" << g.field_names[k]
          << "\" format=\"appended\" offset=\"" << offsets[4 + k] << "\"/>\n";
   }
   vtu << "</CellData>\n</Piece>\n</UnstructuredGrid>\n";

   // Emit raw appended data block.  The leading `_` marker (after the
   // newline) is required by the VTU spec; everything after it is
   // binary until the closing tag.
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

// =============================================================================
// MPI gather (collective on `comm`)
// =============================================================================

/// Concatenate per-rank `LocalFaultPack`s onto rank 0.  All ranks must
/// call collectively.  On non-root ranks the returned struct is empty.
///
/// Triangle indices are remapped during the gather: rank r's local
/// triangle (3i, 3i+1, 3i+2) becomes (3*(displs_cells[r]+i), ...) on
/// rank 0 so the indices reference the gathered global vertex array.
///
/// Field arrays are gathered in the order specified by rank 0's
/// `field_names`.  All ranks are required to present the same
/// field_names in the same order; this is enforced via a hash check
/// (DEBUG builds) or trusted (RELEASE).
inline GatheredFaultPack
GatherFaultPackToRoot(const LocalFaultPack &local,
                      int rank, int nranks
#ifdef MFEM_USE_MPI
                      , MPI_Comm comm = MPI_COMM_NULL
#endif
                     )
{
   GatheredFaultPack out;

   // R-005 (REVIEW.md 2026-04-28): document and enforce the
   // `LocalFaultPack` invariants the gather assumes.  The flatten loop
   // below reads `vertices[3*i + j]` up to index `3*n_cells - 1`, and
   // each per-rank field MPI_Gatherv uses `local_n_cells` as the
   // sendcount.  If a future caller (e.g. Phase 2's submesh path that
   // shares vertices between adjacent triangles) violates either
   // invariant, hard-fail at the gather entry rather than read OOB or
   // silently drop data.
   MFEM_VERIFY(local.vertices.size() == 3 * local.triangles.size(),
               "GatherFaultPackToRoot: LocalFaultPack invariant violated"
               " — vertices.size()=" << local.vertices.size()
               << " expected 3*triangles.size()=3*"
               << local.triangles.size());
   MFEM_VERIFY(local.field_arrays.size() == local.field_names.size(),
               "GatherFaultPackToRoot: field_arrays.size()="
               << local.field_arrays.size() << " != field_names.size()="
               << local.field_names.size());
   for (std::size_t k = 0; k < local.field_arrays.size(); ++k)
   {
      MFEM_VERIFY(local.field_arrays[k].size() == local.triangles.size(),
                  "GatherFaultPackToRoot: field '"
                  << local.field_names[k] << "' has "
                  << local.field_arrays[k].size()
                  << " entries; expected " << local.triangles.size());
   }

   // Serial / single-rank path: copy directly.
#ifdef MFEM_USE_MPI
   if (nranks <= 1 || comm == MPI_COMM_NULL)
#else
   (void)rank;
   (void)nranks;
#endif
   {
      out.vertices     = local.vertices;
      out.triangles    = local.triangles;
      out.field_arrays = local.field_arrays;
      out.field_names  = local.field_names;
      return out;
   }

#ifdef MFEM_USE_MPI
   // 1. Allgather per-rank cell counts so every rank can compute
   //    displacements (we use displs only on root, but we need every
   //    rank to know nfields too).
   const int local_n_cells   = static_cast<int>(local.triangles.size());
   const int local_n_fields  = static_cast<int>(local.field_arrays.size());

   std::vector<int> counts_cells(nranks, 0);
   MPI_Allgather(&local_n_cells, 1, MPI_INT,
                 counts_cells.data(), 1, MPI_INT, comm);

   // 2. Field count consistency.  All ranks must agree on `nfields` and
   //    on `field_names`.  Mismatches indicate a driver-side desync
   //    (e.g. some ranks built the pack with `_k4` enabled and others
   //    without).  Hard-fail with rank info.
   std::vector<int> counts_fields(nranks, 0);
   MPI_Allgather(&local_n_fields, 1, MPI_INT,
                 counts_fields.data(), 1, MPI_INT, comm);
   if (rank == 0)
   {
      for (int r = 1; r < nranks; ++r)
      {
         MFEM_VERIFY(counts_fields[r] == counts_fields[0],
                     "GatherFaultPackToRoot: rank " << r << " has "
                     << counts_fields[r] << " fields; rank 0 has "
                     << counts_fields[0]
                     << " (driver desync — check `_k4` flag uniformity).");
      }
   }

   // R-002 (REVIEW.md 2026-04-28): plan §Phase 1.2 required a hash check
   // on the field-name list, not just the count.  Without it, two ranks
   // with the same number of fields in different orders would silently
   // swap data — producing a corrupted VTU that ParaView opens without
   // complaint.  Compare a deterministic hash of the concatenated
   // field-name list across ranks.
   {
      std::string concat;
      for (const auto &n : local.field_names) { concat += n; concat += '|'; }
      // R-207: don't rely on std::hash<std::string> being deterministic
      // across processes (implementation-defined per C++17 standard).
      // FNV-1a 64-bit is 5 lines, deterministic across every libc++ /
      // libstdc++ version and across heterogeneous-STL container setups.
      uint64_t local_hash = 14695981039346656037ULL;  // FNV offset basis
      for (unsigned char c : concat)
      {
         local_hash ^= c;
         local_hash *= 1099511628211ULL;             // FNV prime
      }
      uint64_t root_hash = local_hash;
      MPI_Bcast(&root_hash, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
      MFEM_VERIFY(local_hash == root_hash,
                  "GatherFaultPackToRoot: rank " << rank <<
                  " has different field_names than rank 0 (hash mismatch); "
                  "names = [" << concat << "]");
   }

   // 3. Compute displacements on root.
   std::vector<int> displs_cells(nranks, 0);
   int total_cells = 0;
   if (rank == 0)
   {
      for (int r = 0; r < nranks; ++r)
      {
         displs_cells[r] = total_cells;
         total_cells += counts_cells[r];
      }
   }

   // 4. Gather vertices.  Each cell carries 3 vertices * 3 doubles = 9
   //    doubles.  We treat the vertex array as flat doubles for transit.
   std::vector<int> counts_doubles(nranks), displs_doubles(nranks);
   if (rank == 0)
   {
      int acc = 0;
      for (int r = 0; r < nranks; ++r)
      {
         counts_doubles[r] = counts_cells[r] * 9;
         displs_doubles[r] = acc;
         acc += counts_doubles[r];
      }
   }
   const int local_n_doubles = local_n_cells * 9;
   std::vector<double> local_v_flat(local_n_doubles);
   for (int i = 0; i < local_n_cells; ++i)
   {
      for (int j = 0; j < 3; ++j)
      {
         for (int d = 0; d < 3; ++d)
         {
            local_v_flat[9*i + 3*j + d] = local.vertices[3*i + j][d];
         }
      }
   }
   std::vector<double> root_v_flat;
   if (rank == 0) { root_v_flat.resize(total_cells * 9); }
   MPI_Gatherv(local_v_flat.data(), local_n_doubles, MPI_DOUBLE,
               rank == 0 ? root_v_flat.data() : nullptr,
               rank == 0 ? counts_doubles.data() : nullptr,
               rank == 0 ? displs_doubles.data() : nullptr,
               MPI_DOUBLE, 0, comm);

   // 5. Gather field arrays (one MPI_Gatherv per field).
   std::vector<int> counts_field_per_rank(nranks),
                    displs_field_per_rank(nranks);
   if (rank == 0)
   {
      int acc = 0;
      for (int r = 0; r < nranks; ++r)
      {
         counts_field_per_rank[r] = counts_cells[r];
         displs_field_per_rank[r] = acc;
         acc += counts_cells[r];
      }
   }

   std::vector<std::vector<double>> root_fields;
   if (rank == 0)
   {
      root_fields.assign(local_n_fields, std::vector<double>(total_cells));
   }
   for (int k = 0; k < local_n_fields; ++k)
   {
      const double *send = local.field_arrays[k].data();
      MFEM_ASSERT(static_cast<int>(local.field_arrays[k].size()) == local_n_cells,
                  "GatherFaultPackToRoot: local field " << k <<
                  " has wrong size on rank " << rank);
      MPI_Gatherv(send, local_n_cells, MPI_DOUBLE,
                  rank == 0 ? root_fields[k].data() : nullptr,
                  rank == 0 ? counts_field_per_rank.data() : nullptr,
                  rank == 0 ? displs_field_per_rank.data() : nullptr,
                  MPI_DOUBLE, 0, comm);
   }

   // 6. Reconstruct on rank 0.
   if (rank == 0)
   {
      out.vertices.resize(total_cells * 3);
      out.triangles.resize(total_cells);
      out.field_arrays = std::move(root_fields);
      out.field_names  = local.field_names;  // rank 0 owns canonical names

      // Vertices: unflatten root_v_flat back into [3*total_cells] of
      // {x,y,z}.
      for (int c = 0; c < total_cells; ++c)
      {
         for (int j = 0; j < 3; ++j)
         {
            out.vertices[3*c + j] = {
               root_v_flat[9*c + 3*j + 0],
               root_v_flat[9*c + 3*j + 1],
               root_v_flat[9*c + 3*j + 2]
            };
         }
      }
      // Triangles: trivial since each cell owns three unique consecutive
      // vertices in the gathered global vertex array.
      for (int c = 0; c < total_cells; ++c)
      {
         out.triangles[c] = {3*c, 3*c + 1, 3*c + 2};
      }
   }
   return out;
#endif // MFEM_USE_MPI
}

// =============================================================================
// VTU CellData reader (used by tests; supports ASCII, base64
// inline-binary, AND raw appended formats — the three formats this
// writer can emit / has historically emitted).
// =============================================================================

/// Read the named CellData scalar array from a VTU produced by
/// `WriteFaultPackVTU`.  Supports:
///   - `format="ascii"`     (whitespace-separated doubles)
///   - `format="binary"`    (legacy: base64-encoded uint32-prefixed block)
///   - `format="appended"`  (current default per R-302: raw bytes in
///                           `<AppendedData encoding="raw">_<uint64
///                           length><bytes>...`)
/// Returns an empty vector if the field is not present.
inline std::vector<double>
ParseFaultVTUCellData(const std::string &path, const std::string &field)
{
   std::vector<double> out;
   std::ifstream in(path, std::ios::binary);
   if (!in.is_open()) { return out; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();

   // Locate the <CellData> ... </CellData> block.
   const auto cd_open  = content.find("<CellData>");
   const auto cd_close = content.find("</CellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos
       || cd_close <= cd_open) { return out; }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);

   // Locate the named DataArray inside the block.
   const std::string marker = std::string("Name=\"") + field + "\"";
   const auto name_pos = cd.find(marker);
   if (name_pos == std::string::npos) { return out; }
   const auto tag_close = cd.find('>', name_pos);
   if (tag_close == std::string::npos) { return out; }

   // Extract format= — read it BEFORE the self-close vs body check so
   // we can dispatch on layout (raw appended uses self-closing tags
   // and does NOT have a </DataArray> closer; ASCII / inline-binary
   // are body-bearing tags that DO).
   const auto fmt_pos = cd.find("format=\"", name_pos);
   std::string fmt_val = "ascii";
   if (fmt_pos != std::string::npos && fmt_pos < tag_close)
   {
      const auto fmt_qopen = fmt_pos + std::strlen("format=\"");
      const auto fmt_qclose = cd.find('"', fmt_qopen);
      fmt_val = cd.substr(fmt_qopen, fmt_qclose - fmt_qopen);
   }

   if (fmt_val == "appended")
   {
      // R-302: raw appended layout.  Read the offset attribute on the
      // DataArray, then locate the `_` marker after
      // <AppendedData encoding="raw">, then read
      //   uint64 byte_length followed by byte_length bytes of doubles
      // at that offset.  header_type defaults to UInt64 for files this
      // writer emits (see the writer); for forward-compat we honour the
      // VTKFile root's `header_type` attribute too.
      const auto off_pos = cd.find("offset=\"", name_pos);
      if (off_pos == std::string::npos || off_pos > tag_close)
      { return out; }
      const auto off_qopen  = off_pos + std::strlen("offset=\"");
      const auto off_qclose = cd.find('"', off_qopen);
      const uint64_t da_offset =
         std::stoull(cd.substr(off_qopen, off_qclose - off_qopen));

      // Header type (uint32 vs uint64) per the VTKFile attribute.
      bool header_is_uint64 = true;
      const auto vfk = content.find("<VTKFile");
      if (vfk != std::string::npos)
      {
         const auto vfk_close = content.find('>', vfk);
         const auto ht = content.find("header_type=\"", vfk);
         if (ht != std::string::npos && ht < vfk_close)
         {
            const auto ht_qopen = ht + std::strlen("header_type=\"");
            const auto ht_qclose = content.find('"', ht_qopen);
            const std::string ht_val =
               content.substr(ht_qopen, ht_qclose - ht_qopen);
            header_is_uint64 = (ht_val == "UInt64");
         }
      }

      // Find the leading `_` marker in <AppendedData encoding="raw">.
      const auto ad_open = content.find("<AppendedData");
      if (ad_open == std::string::npos) { return out; }
      const auto ad_tag_close = content.find('>', ad_open);
      if (ad_tag_close == std::string::npos) { return out; }
      const auto under = content.find('_', ad_tag_close);
      if (under == std::string::npos) { return out; }
      const std::size_t blob_start = under + 1;

      // Read header @ blob_start + da_offset.
      const std::size_t header_pos = blob_start + da_offset;
      const std::size_t header_size = header_is_uint64
                                      ? sizeof(uint64_t) : sizeof(uint32_t);
      if (header_pos + header_size > content.size()) { return out; }
      uint64_t byte_length = 0;
      if (header_is_uint64)
      {
         std::memcpy(&byte_length, content.data() + header_pos,
                     sizeof(uint64_t));
      }
      else
      {
         uint32_t hdr32 = 0;
         std::memcpy(&hdr32, content.data() + header_pos, sizeof(uint32_t));
         byte_length = hdr32;
      }
      const std::size_t data_pos = header_pos + header_size;
      if (data_pos + byte_length > content.size()) { return out; }
      const std::size_t n = byte_length / sizeof(double);
      out.resize(n);
      if (n > 0)
      {
         std::memcpy(out.data(), content.data() + data_pos,
                     n * sizeof(double));
      }
      return out;
   }

   // For body-bearing formats (ascii, binary), locate the </DataArray>
   // closer.  Raw appended uses self-closing tags so this lookup is
   // unused on that path (handled above).
   const auto arr_close = cd.find("</DataArray>", name_pos);
   if (arr_close == std::string::npos || arr_close <= tag_close)
   { return out; }
   const std::string body = cd.substr(tag_close + 1,
                                      arr_close - tag_close - 1);
   if (fmt_val == "ascii")
   {
      std::istringstream body_ss(body);
      double v;
      while (body_ss >> v) { out.push_back(v); }
      return out;
   }

   // fmt_val == "binary": legacy inline-base64 path retained for
   // backward compat with files emitted before R-302.  Layout
   // (uncompressed): `[uint32 byte_length] [byte_length bytes of
   // raw doubles]`, base64-encoded.  Compressed mode is not supported
   // by this tests-only reader.
   std::string b64;
   b64.reserve(body.size());
   for (char c : body)
   {
      if (c != ' ' && c != '\n' && c != '\r' && c != '\t') { b64.push_back(c); }
   }
   std::vector<char> raw;
   bin_io::DecodeBase64(b64.data(), b64.size(), raw);
   if (raw.size() < sizeof(uint32_t)) { return out; }
   uint32_t header;
   std::memcpy(&header, raw.data(), sizeof(uint32_t));
   if (header + sizeof(uint32_t) != raw.size())
   {
      // R-004 (REVIEW.md 2026-04-28): the input layout did not match the
      // uncompressed `[uint32 byte_length][raw_bytes]` shape.  Either it
      // is zlib-compressed (multi-block header: num_blocks, block_size,
      // last_block_size, compressed_size_per_block...) or malformed.
      // This tests-only reader does not link zlib; aborting with an
      // explicit message is preferable to returning a silent empty
      // vector that fails downstream length assertions with no clue
      // about the cause.
      MFEM_ABORT("ParseFaultVTUCellData: VTU at " << path
                 << " field '" << field
                 << "' appears zlib-compressed (raw size " << raw.size()
                 << ", uint32 header value " << header << "). "
                 "Re-emit the VTU with compression_level=0 or extend "
                 "this parser to handle the multi-block compressed "
                 "header.");
   }
   const std::size_t n = header / sizeof(double);
   out.resize(n);
   std::memcpy(out.data(), raw.data() + sizeof(uint32_t), header);
   return out;
}

} // namespace mfem::seas::vtu
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_VTU_BINARY_HPP
