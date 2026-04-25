// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SPDX-License-Identifier: BSD-3-Clause
//
// Fault-locality partitioning helper (G1) — TPV104 dynamic only.
//
// Goal
// ----
// Eliminate the dip-direction pollution observed at np ≥ 4 on the
// y-mirror-symmetric mesh (debug doc 2026-04-25_pm Section 13).  Root
// cause: ParMETIS subdivides along Y, cutting the y=0 fault plane;
// ghost-cell exchanges between +y and −y ranks introduce a topology-
// induced FP non-associativity that the rupture amplifies linearly
// to mm-scale slip_dip.
//
// Fix (G1, generalized)
// ----------------------
// For every fault face F with adjacent elements E_+(F), E_−(F),
// constrain the partition so both elements are co-resident on the
// same rank.  This is the "fault-locality" invariant.  Generalizes
// trivially to non-flat / curved / dipping / branched faults — the
// constraint is per-face, not geometry-dependent.
//
// Implementation: union-find pre-merge of fault-face element pairs
// into super-elements.  Call MFEM's standard partitioner on the
// original elements; for each super-element, all its members are
// reassigned to the rank of its "representative" (= element with
// smallest index in the union-find component).  ParMETIS still gets
// to load-balance; we just lock fault adjacencies post-hoc.
//
// Scope (per user directive 2026-04-25):
//   - Limited to dynamic/ TPV104 path only.
//   - BP5 driver and bp5/ / domain/ / fault/ / solver/ are NOT touched.
//   - Header-only; included from drivers/tpv104_driver.cpp.

#ifndef MFEM_SEAS_FAULT_LOCALITY_PARTITION_HPP
#define MFEM_SEAS_FAULT_LOCALITY_PARTITION_HPP

#include "mfem.hpp"
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Minimal Union-Find with path compression + union by rank.
/// Internal helper for BuildFaultLocalityPartitioning.
class FaultLocalityUF
{
public:
   explicit FaultLocalityUF(int n) : parent_(n), rank_(n, 0)
   {
      for (int i = 0; i < n; i++) { parent_[i] = i; }
   }

   int Find(int x)
   {
      while (parent_[x] != x)
      {
         parent_[x] = parent_[parent_[x]];   // path compression
         x = parent_[x];
      }
      return x;
   }

   void Union(int x, int y)
   {
      int rx = Find(x), ry = Find(y);
      if (rx == ry) { return; }
      if (rank_[rx] < rank_[ry]) { std::swap(rx, ry); }
      parent_[ry] = rx;
      if (rank_[rx] == rank_[ry]) { rank_[rx]++; }
   }

private:
   std::vector<int> parent_;
   std::vector<int> rank_;
};

/// Build a partitioning array of length `mesh.GetNE()` such that, for
/// every face index in `fault_face_indices`, its two adjacent elements
/// are assigned to the same rank.
///
/// Algorithm:
///   1. Run MFEM's standard partitioner (METIS via Mesh::GeneratePartitioning).
///   2. Build union-find over elements; merge each fault face's pair.
///   3. For each union-find component, assign all members to the rank of
///      the representative (lowest-index element).
///
/// On TPV104's flat fault, each fault face merges exactly 2 elements
/// (each tet has at most 1 fault face).  The reassignment moves
/// at most one element of each fault pair; load balance is preserved
/// to within a few elements.
///
/// For more complex fault topologies (branches, intersections, curves),
/// the union-find handles transitive merges automatically.  In the
/// limit where every element is fault-adjacent (impossible in practice),
/// all elements would land on one rank.
///
/// @param[in]  mesh                 Serial mesh.
/// @param[in]  fault_face_indices   Mesh face indices that are faults
///                                  (typically derived from the mesh
///                                  boundary attribute, e.g. attr=3).
/// @param[in]  np                   Target number of ranks.
/// @param[out] partitioning_out     Length-mesh.GetNE() array; caller
///                                  owns.  Indices in [0, np).
/// @param[out] n_pairs_relocated    Diagnostic: how many fault pairs
///                                  required a rank move.  Optional;
///                                  pass nullptr if not needed.
inline void BuildFaultLocalityPartitioning(
   Mesh &mesh,
   const Array<int> &fault_face_indices,
   int np,
   Array<int> &partitioning_out,
   int *n_pairs_relocated = nullptr)
{
   const int ne = mesh.GetNE();
   MFEM_VERIFY(np >= 1, "BuildFaultLocalityPartitioning: np must be >= 1");

   // Trivial case: serial.  All on rank 0.
   if (np == 1)
   {
      partitioning_out.SetSize(ne);
      partitioning_out = 0;
      if (n_pairs_relocated) { *n_pairs_relocated = 0; }
      return;
   }

   // Step 1: standard METIS partitioning as starting point.
   // Mesh::GeneratePartitioning returns a heap-allocated int* of length ne.
   int *base_part = mesh.GeneratePartitioning(np, /*part_method=*/1);
   MFEM_VERIFY(base_part != nullptr,
               "BuildFaultLocalityPartitioning: GeneratePartitioning failed");

   // Step 2: union-find over elements; merge each fault pair.
   FaultLocalityUF uf(ne);
   int n_uf_unions = 0;
   for (int i = 0; i < fault_face_indices.Size(); i++)
   {
      const int f = fault_face_indices[i];
      // GetFaceElementTransformations returns interior face info if
      // both sides are present.  For boundary or shared faces (parallel),
      // only one side is local; skip — the caller must rebuild the
      // serial-mesh face list before partitioning.
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr == nullptr) { continue; }
      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      if (e1 < 0 || e2 < 0) { continue; }
      if (uf.Find(e1) != uf.Find(e2))
      {
         uf.Union(e1, e2);
         n_uf_unions++;
      }
   }

   // Step 3: for each UF component, lock all members to the rep's rank.
   // Representative = element with smallest index in the component.
   // We track the rep's rank per component using a map from UF root.
   std::map<int, int> root_to_rank;
   int relocated = 0;
   for (int e = 0; e < ne; e++)
   {
      const int r = uf.Find(e);
      auto it = root_to_rank.find(r);
      if (it == root_to_rank.end())
      {
         // First time we see this component → its rank is base_part[e].
         // Since we iterate e in order, this is the lowest-index member.
         root_to_rank[r] = base_part[e];
      }
   }

   partitioning_out.SetSize(ne);
   for (int e = 0; e < ne; e++)
   {
      const int r = uf.Find(e);
      const int new_rank = root_to_rank[r];
      if (new_rank != base_part[e]) { relocated++; }
      partitioning_out[e] = new_rank;
   }

   delete[] base_part;
   if (n_pairs_relocated) { *n_pairs_relocated = relocated; }
}

/// Verify the fault-locality invariant: for every fault face F, the two
/// adjacent elements have the same rank.  Returns the number of
/// violating face pairs (0 = success).  Useful for tests / asserts.
inline int VerifyFaultLocality(Mesh &mesh,
                               const Array<int> &fault_face_indices,
                               const Array<int> &partitioning)
{
   int violations = 0;
   for (int i = 0; i < fault_face_indices.Size(); i++)
   {
      const int f = fault_face_indices[i];
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr == nullptr) { continue; }
      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      if (e1 < 0 || e2 < 0) { continue; }
      if (partitioning[e1] != partitioning[e2]) { violations++; }
   }
   return violations;
}

/// Load a partitioning array from a sidecar file written by
/// `build_symmirror_mesh.py --np N --emit-partition`.  File format:
///   Line 1:  "np <N>" — target rank count.
///   Line 2:  "ne <M>" — number of elements (must match mesh.GetNE()).
///   Lines 3..M+2:  one integer per line, partition[e] in [0, N).
///
/// On mismatch (np ≠ comm-size, ne ≠ mesh.GetNE()) returns false and
/// leaves partitioning_out empty.  Otherwise returns true with the
/// loaded partition.
inline bool LoadPartitioningFromFile(const std::string &path, int np_expected,
                                      int ne_expected,
                                      Array<int> &partitioning_out)
{
   std::ifstream in(path);
   if (!in.is_open()) { return false; }
   std::string tag; int np_file = -1, ne_file = -1;
   in >> tag >> np_file;
   if (tag != "np" || np_file != np_expected) { return false; }
   in >> tag >> ne_file;
   if (tag != "ne" || ne_file != ne_expected) { return false; }
   partitioning_out.SetSize(ne_file);
   for (int e = 0; e < ne_file; e++)
   {
      int p = -1;
      if (!(in >> p)) { partitioning_out.SetSize(0); return false; }
      if (p < 0 || p >= np_expected)
      {
         partitioning_out.SetSize(0);
         return false;
      }
      partitioning_out[e] = p;
   }
   return true;
}

/// Discover fault-face indices in a serial mesh by scanning interior
/// faces whose `bdr_attr` matches `fault_attr` (TPV104 convention:
/// fault_attr = 3).  Caller may use this output as the
/// `fault_face_indices` input to BuildFaultLocalityPartitioning.
///
/// MFEM stores boundary attributes on EXTERIOR faces by default.  For
/// the embedded fault face (Y=0 plane interior to a 2-sided element
/// pair), the .msh format places a triangle element with bdr_attr=3 at
/// the fault face position, and MFEM picks this up via
/// GetFaceInformation / GetBdrAttribute or via the mesh element table.
/// We replicate the same scan that wave_operator.inl does to populate
/// fault_interior_faces_.
inline void FindFaultFaceIndices(Mesh &mesh, int fault_attr,
                                 Array<int> &fault_face_indices_out)
{
   fault_face_indices_out.SetSize(0);
   // Iterate over boundary elements (which include embedded fault triangles
   // tagged with attr=fault_attr).  Get each bdr element's underlying face
   // index, check attribute, then verify the face is interior (2-sided).
   // Same pattern as wave_operator.inl:120-124.
   const int nbe = mesh.GetNBE();
   for (int b = 0; b < nbe; b++)
   {
      const int attr = mesh.GetBdrAttribute(b);
      if (attr != fault_attr) { continue; }
      const int f = mesh.GetBdrElementFaceIndex(b);
      if (f < 0) { continue; }
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr == nullptr) { continue; }
      if (ftr->Elem2No < 0) { continue; }   // not interior
      fault_face_indices_out.Append(f);
   }
}

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_FAULT_LOCALITY_PARTITION_HPP
