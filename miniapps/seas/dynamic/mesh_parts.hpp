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
#include <memory>
#include <sstream>
#include <string>
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
/// coordinates bit for bit.  Tetrahedral meshes.
inline void CheckMeshPart(const Mesh &smesh, const MeshPart &mp, const int *ids, int n, int r)
{
   MFEM_VERIFY(mp.num_elements == n, "mesh parts: part " << r << " has " << mp.num_elements
               << " elements; the partition gives it " << n);
   const Array<int> &tv = mp.entity_to_vertex[Geometry::TETRAHEDRON];
   MFEM_VERIFY(tv.Size() == 4 * n, "mesh parts: part " << r << " is not all tetrahedra");
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
      bool same = mp.attributes[k] == smesh.GetAttribute(e) && sv.Size() == 4;
      for (int j = 0; same && j < 4; j++)
      {
         const real_t *a = &mp.vertex_coordinates[sdim * tv[4 * bt + j]];
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

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MESH_PARTS_HPP
