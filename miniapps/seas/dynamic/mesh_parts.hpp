// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/mesh_parts.hpp — the parallel mesh of the Nitsche dynamic-rupture
// driver with the serial mesh on rank 0 only (Chunhui Zhao 2026-10-22:
// partition-memory report rev1, option D; speed-next S1).  Rank 0 holds the
// serial mesh and its partition, cuts the mesh into one part per rank with
// MFEM's MeshPartitioner and sends each rank its part, as text in MFEM's
// parallel mesh format ("MFEM mesh v1.2", vertex coordinates to 17 digits,
// so they read back bit for bit), together with one integer per element of
// the part (the LTS clusters).  Every rank builds its ParMesh from its part
// with refine = false and fix_orientation = false, so the element and
// boundary vertex order is the one the serial mesh chose, as in
// ParMesh(comm, serial mesh, partition).
//
// The local elements of part r are the serial elements e with part[e] == r
// in increasing e (MeshPartitioner's order, the order of ParMesh(comm, serial
// mesh, partition)); rank 0 checks that for every element of every part
// (attribute and the coordinates of each vertex, exactly) before it sends
// the part.  The parallel mesh then has the same local vertices, elements,
// boundary elements, faces and DOF numbers as ParMesh(comm, serial mesh,
// partition); only its shared faces are listed group by group instead of in
// serial face order (tests/verification/test_mesh_parts.cpp).

#ifndef MFEM_SEAS_MESH_PARTS_HPP
#define MFEM_SEAS_MESH_PARTS_HPP

#include "mfem.hpp"

#include <climits>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// The serial elements of each part, in increasing element number:
/// elements of part r = ids[off[r]] .. ids[off[r+1]-1].
inline void MeshPartElements(int ne, const int *part, int nparts, std::vector<int> &off,
                             std::vector<int> &ids)
{
   off.assign(nparts + 1, 0);
   for (int e = 0; e < ne; e++)
   {
      MFEM_VERIFY(part[e] >= 0 && part[e] < nparts, "mesh parts: element " << e
                  << " is in part " << part[e] << " of " << nparts);
      off[part[e] + 1]++;
   }
   for (int r = 0; r < nparts; r++) { off[r + 1] += off[r]; }
   ids.resize(ne);
   std::vector<int> at(off.begin(), off.end() - 1);
   for (int e = 0; e < ne; e++) { ids[at[part[e]]++] = e; }
}

/// Check that the elements of `mp` (part r) are the serial elements ids[0..n-1]
/// of `smesh`, in that order: same attribute and, vertex by vertex, the same
/// coordinates bit for bit.  Tetrahedral meshes, or triangle meshes in the
/// (x, y) plane (2D in-plane design rev 3; Chunhui Zhao, 2026-10-27).
inline void CheckMeshPart(const Mesh &smesh, const MeshPart &mp, const int *ids, int n, int r)
{
   MFEM_VERIFY(mp.num_elements == n, "mesh parts: part " << r << " has " << mp.num_elements
               << " elements; the partition gives it " << n);
   const bool tri = smesh.Dimension() == 2;
   const Geometry::Type geom = tri ? Geometry::TRIANGLE : Geometry::TETRAHEDRON;
   const int nv = tri ? 3 : 4;   // vertices per element
   const Array<int> &tv = mp.entity_to_vertex[geom];
   MFEM_VERIFY(tv.Size() == nv * n, "mesh parts: part " << r << " is not all "
               << (tri ? "triangles" : "tetrahedra"));
   const int sdim = mp.space_dimension;
   MFEM_VERIFY(mp.vertex_coordinates.Size() == sdim * mp.num_vertices,
               "mesh parts: part " << r << " has no vertex coordinates");
   const bool map = mp.element_map.Size() == n;
   Array<int> sv;
   for (int k = 0; k < n; k++)
   {
      const int e = ids[k];
      const int bt = map ? mp.element_map[k] : k;
      smesh.GetElementVertices(e, sv);
      bool same = mp.attributes[k] == smesh.GetAttribute(e) && sv.Size() == nv;
      for (int j = 0; same && j < nv; j++)
      {
         const real_t *a = &mp.vertex_coordinates[sdim * tv[nv * bt + j]];
         const real_t *b = smesh.GetVertex(sv[j]);
         for (int d = 0; d < sdim; d++) { same = same && a[d] == b[d]; }
      }
      MFEM_VERIFY(same, "mesh parts: element " << k << " of part " << r
                  << " is not serial element " << e << " (attribute or vertices differ)");
   }
}

/// Collective over `comm`.  On rank 0: the serial mesh `smesh`, its partition
/// `part` into the ranks of `comm` and, when `with_values`, `values` (one
/// integer per serial element); ignored on the other ranks.  Returns this
/// rank's ParMesh; with `with_values`, `local_values` gets the values of its
/// local elements, in local element order.  `with_values` must be the same
/// on every rank.
inline std::unique_ptr<ParMesh> DistributeMeshParts(MPI_Comm comm, Mesh *smesh,
                                                   const int *part, bool with_values,
                                                   const int *values,
                                                   std::vector<int> *local_values)
{
   int rank, np;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &np);
   const int tag_text = 7101, tag_values = 7102;
   std::string mine;
   std::vector<int> mine_values;
   if (rank == 0)
   {
      MFEM_VERIFY(smesh && part, "mesh parts: rank 0 needs the serial mesh and its partition");
      MFEM_VERIFY(!with_values || values, "mesh parts: rank 0 needs the element values");
      std::vector<int> off, ids;
      MeshPartElements(smesh->GetNE(), part, np, off, ids);
      MeshPartitioner partitioner(*smesh, np, part);
      MeshPart mp;
      std::vector<int> v;
      // rank 0's own part last, so the other ranks start building theirs
      for (int i = 1; i <= np; i++)
      {
         const int r = i % np;
         const int n = off[r + 1] - off[r];
         partitioner.ExtractPart(r, mp);
         CheckMeshPart(*smesh, mp, ids.data() + off[r], n, r);
         std::ostringstream os;
         os.precision(17);
         mp.Print(os);
         std::string text = os.str();
         MFEM_VERIFY(text.size() < static_cast<size_t>(INT_MAX), "mesh parts: part " << r
                     << " is " << text.size() << " bytes, more than one MPI message holds; "
                     "use more ranks");
         if (with_values)
         {
            v.resize(n);
            for (int k = 0; k < n; k++) { v[k] = values[ids[off[r] + k]]; }
         }
         if (r == 0)
         {
            mine.swap(text);
            mine_values.swap(v);
         }
         else
         {
            MPI_Send(text.data(), static_cast<int>(text.size()), MPI_CHAR, r, tag_text, comm);
            if (with_values) { MPI_Send(v.data(), n, MPI_INT, r, tag_values, comm); }
         }
      }
   }
   else
   {
      MPI_Status st;
      int n = 0;
      MPI_Probe(0, tag_text, comm, &st);
      MPI_Get_count(&st, MPI_CHAR, &n);
      mine.resize(n);
      MPI_Recv(&mine[0], n, MPI_CHAR, 0, tag_text, comm, MPI_STATUS_IGNORE);
      if (with_values)
      {
         MPI_Probe(0, tag_values, comm, &st);
         MPI_Get_count(&st, MPI_INT, &n);
         mine_values.resize(n);
         MPI_Recv(mine_values.data(), n, MPI_INT, 0, tag_values, comm, MPI_STATUS_IGNORE);
      }
   }
   std::unique_ptr<ParMesh> pmesh;
   {
      std::istringstream is(mine);
      mine.clear();
      mine.shrink_to_fit();
      pmesh = std::make_unique<ParMesh>(comm, is, /*refine=*/false, /*generate_edges=*/1,
                                        /*fix_orientation=*/false);
   }
   if (with_values)
   {
      MFEM_VERIFY(static_cast<int>(mine_values.size()) == pmesh->GetNE(), "mesh parts: rank "
                  << rank << " has " << pmesh->GetNE() << " elements but "
                  << mine_values.size() << " element values");
      local_values->swap(mine_values);
   }
   return pmesh;
}

// ---------------------------------------------------------------- part files
// (speed-next S1 step 2, Chunhui Zhao 2026-10-22).  A directory of part
// files holds, for parts r = 0 .. N-1, mesh.<r> (MakeParFilename: six
// digits; the text DistributeMeshParts sends), values.<r> (the element
// values, one per line, when there are any) and mesh_parts.txt, one
// "key value" line per key, which the caller fills (what the parts were cut
// from and with which options) and checks when it reads them.

/// Rank 0 alone, no communication: cut `smesh` into the `nparts` parts of
/// `part` and write them to `dir`, with the element values when
/// `with_values`; every part is checked as in DistributeMeshParts.
inline void WriteMeshParts(const std::string &dir, Mesh &smesh, const int *part, int nparts,
                           bool with_values, const int *values)
{
   std::vector<int> off, ids;
   MeshPartElements(smesh.GetNE(), part, nparts, off, ids);
   MeshPartitioner partitioner(smesh, nparts, part);
   MeshPart mp;
   for (int r = 0; r < nparts; r++)
   {
      const int n = off[r + 1] - off[r];
      partitioner.ExtractPart(r, mp);
      CheckMeshPart(smesh, mp, ids.data() + off[r], n, r);
      const std::string f = MakeParFilename(dir + "/mesh.", r);
      std::ofstream os(f);
      os.precision(17);
      mp.Print(os);
      os.close();
      MFEM_VERIFY(os, "mesh parts: cannot write " << f);
      if (with_values)
      {
         const std::string fv = MakeParFilename(dir + "/values.", r);
         std::ofstream ov(fv);
         ov << n << "\n";
         for (int k = 0; k < n; k++) { ov << values[ids[off[r] + k]] << "\n"; }
         ov.close();
         MFEM_VERIFY(ov, "mesh parts: cannot write " << fv);
      }
   }
}

/// Collective over `comm`: every rank builds its ParMesh from its own part
/// file in `dir` (refine = false, fix_orientation = false, as
/// DistributeMeshParts) and, with `with_values`, reads its element values.
inline std::unique_ptr<ParMesh> ReadMeshPart(MPI_Comm comm, const std::string &dir,
                                            bool with_values, std::vector<int> *local_values)
{
   int rank;
   MPI_Comm_rank(comm, &rank);
   const std::string f = MakeParFilename(dir + "/mesh.", rank);
   std::ifstream is(f);
   MFEM_VERIFY(is, "mesh parts: rank " << rank << " cannot open " << f);
   auto pmesh = std::make_unique<ParMesh>(comm, is, /*refine=*/false, /*generate_edges=*/1,
                                          /*fix_orientation=*/false);
   if (with_values)
   {
      const std::string fv = MakeParFilename(dir + "/values.", rank);
      std::ifstream iv(fv);
      int n = -1;
      iv >> n;
      MFEM_VERIFY(iv && n == pmesh->GetNE(), "mesh parts: " << fv << " holds " << n
                  << " element values; rank " << rank << " has " << pmesh->GetNE()
                  << " elements");
      local_values->resize(n);
      for (int k = 0; k < n; k++) { iv >> (*local_values)[k]; }
      MFEM_VERIFY(iv, "mesh parts: cannot read the " << n << " values of " << fv);
   }
   return pmesh;
}

/// mesh_parts.txt in `dir`: one "key value" line per entry, in order.
inline void WriteMeshPartsInfo(const std::string &dir,
                               const std::vector<std::pair<std::string, std::string>> &kv)
{
   const std::string f = dir + "/mesh_parts.txt";
   std::ofstream os(f);
   for (const auto &e : kv) { os << e.first << " " << e.second << "\n"; }
   os.close();
   MFEM_VERIFY(os, "mesh parts: cannot write " << f);
}

/// Collective over `comm`: rank 0 reads mesh_parts.txt in `dir` and every
/// rank gets its entries (key, rest of the line).
inline std::map<std::string, std::string> ReadMeshPartsInfo(MPI_Comm comm,
                                                            const std::string &dir)
{
   int rank;
   MPI_Comm_rank(comm, &rank);
   std::string text;
   int n = 0;
   if (rank == 0)
   {
      const std::string f = dir + "/mesh_parts.txt";
      std::ifstream is(f);
      MFEM_VERIFY(is, "mesh parts: cannot open " << f);
      std::ostringstream ss;
      ss << is.rdbuf();
      text = ss.str();
      n = static_cast<int>(text.size());
   }
   MPI_Bcast(&n, 1, MPI_INT, 0, comm);
   text.resize(n);
   MPI_Bcast(&text[0], n, MPI_CHAR, 0, comm);
   std::map<std::string, std::string> kv;
   std::istringstream is(text);
   for (std::string line; std::getline(is, line);)
   {
      const size_t sp = line.find(' ');
      if (line.empty() || sp == std::string::npos) { continue; }
      kv[line.substr(0, sp)] = line.substr(sp + 1);
   }
   return kv;
}

/// 64-bit FNV-1a of a file's bytes, as 16 hex digits: names the mesh the
/// parts were cut from.
inline std::string FileFnv1a64(const std::string &path)
{
   std::ifstream is(path, std::ios::binary);
   MFEM_VERIFY(is, "mesh parts: cannot open " << path);
   unsigned long long h = 1469598103934665603ULL;
   std::vector<char> buf(1 << 20);
   while (is)
   {
      is.read(buf.data(), static_cast<std::streamsize>(buf.size()));
      const std::streamsize got = is.gcount();
      for (std::streamsize i = 0; i < got; i++)
      {
         h ^= static_cast<unsigned char>(buf[i]);
         h *= 1099511628211ULL;
      }
   }
   std::ostringstream os;
   os << std::hex << std::setw(16) << std::setfill('0') << h;
   return os.str();
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MESH_PARTS_HPP
