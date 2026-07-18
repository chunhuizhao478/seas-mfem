// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_layout.cpp — implementation of Phase 1 run-side layout (Appendix A.2 of
// PLAN_clustered_lts_ader_2026-07-18.md).

#include "lts_layout.hpp"

#include "mfem.hpp"   // MFEM_VERIFY

#include <algorithm>

namespace mfem
{
namespace seas
{

long long LtsLayout::num_owned_faces() const
{
   long long n = 0;
   for (const LtsCluster& c : clusters)
   {
      n += static_cast<long long>(c.face_ids.size());
   }
   return n;
}


LtsLayout BuildLtsLayout(const std::vector<int>&         cluster,
                         int                             num_clusters,
                         const std::vector<LtsFaceSpec>& faces)
{
   const int ne = static_cast<int>(cluster.size());
   MFEM_VERIFY(ne > 0, "BuildLtsLayout: cluster vector is empty.");
   MFEM_VERIFY(num_clusters >= 1,
               "BuildLtsLayout: num_clusters must be >= 1 (got "
               << num_clusters << ").");

   int max_c = 0;
   for (int i = 0; i < ne; ++i)
   {
      MFEM_VERIFY(cluster[i] >= 0 && cluster[i] < num_clusters,
                  "BuildLtsLayout: cluster[" << i << "] = " << cluster[i]
                  << " out of range [0," << num_clusters << ").");
      max_c = std::max(max_c, cluster[i]);
   }
   MFEM_VERIFY(max_c == num_clusters - 1,
               "BuildLtsLayout: num_clusters (" << num_clusters
               << ") != max cluster id + 1 (" << (max_c + 1)
               << ") — ids must be contiguous.");

   LtsLayout L;
   L.num_clusters = num_clusters;
   L.clusters.assign(static_cast<std::size_t>(num_clusters), LtsCluster());
   L.local_elems_per_cluster.assign(num_clusters, 0);
   L.local_fault_faces_per_cluster.assign(num_clusters, 0);

   // Element lists per cluster (ascending element id — deterministic).
   for (int e = 0; e < ne; ++e)
   {
      L.clusters[static_cast<std::size_t>(cluster[e])].elems.push_back(e);
      ++L.local_elems_per_cluster[cluster[e]];
   }

   // Mark coarse-side elements of diff-1 faces (the provider / buffer-owner set)
   // while assigning per-side face roles.
   std::vector<char> is_coarse_side(static_cast<std::size_t>(ne), 0);

   auto add_face = [&](int cl_id, int face_id, FaceRole role, int nbr)
   {
      LtsCluster& c = L.clusters[static_cast<std::size_t>(cl_id)];
      c.face_ids.push_back(face_id);
      c.faces.push_back(role);
      c.face_nbr.push_back(nbr);
   };

   for (const LtsFaceSpec& f : faces)
   {
      const int e1 = f.elem1;
      MFEM_VERIFY(e1 >= 0 && e1 < ne,
                  "BuildLtsLayout: face " << f.face_id << " elem1 " << e1
                  << " out of range [0," << ne << ").");
      const int c1 = cluster[e1];

      if (f.elem2 < 0)
      {
         // Boundary / rank seam: owned by e1's cluster, no neighbour.
         add_face(c1, f.face_id, FaceRole::Boundary, -1);
         continue;
      }
      const int e2 = f.elem2;
      MFEM_VERIFY(e2 >= 0 && e2 < ne,
                  "BuildLtsLayout: face " << f.face_id << " elem2 " << e2
                  << " out of range [0," << ne << ").");
      const int c2 = cluster[e2];

      if (f.is_fault)
      {
         MFEM_VERIFY(c1 == c2,
                     "BuildLtsLayout: fault face " << f.face_id
                     << " joins clusters " << c1 << " and " << c2
                     << " (fault faces must be diff 0).");
         // Fault faces are tracked separately in BOTH incident clusters' lists
         // (here c1 == c2, so once).  They are handled by the fault machinery,
         // not the bulk face sweep.
         LtsCluster& c = L.clusters[static_cast<std::size_t>(c1)];
         c.fault_faces.push_back(f.face_id);
         ++L.local_fault_faces_per_cluster[c1];
         continue;
      }

      const int d = c1 - c2;
      MFEM_VERIFY(d >= -1 && d <= 1,
                  "BuildLtsLayout: non-fault face " << f.face_id
                  << " joins clusters " << c1 << " and " << c2
                  << " differing by " << (d < 0 ? -d : d)
                  << " (maxdiff <= 1 violated).");

      if (d == 0)
      {
         // GTS face: owned once by the shared cluster (both sides scatter).
         add_face(c1, f.face_id, FaceRole::IntraClusterGTS, e2);
      }
      else if (d == 1)
      {
         // e1 coarser (larger cluster id => larger dt), e2 finer.
         add_face(c1, f.face_id, FaceRole::ProviderCoarseSkip, e2);
         add_face(c2, f.face_id, FaceRole::ConsumerFine,       e1);
         is_coarse_side[static_cast<std::size_t>(e1)] = 1;
      }
      else // d == -1: e2 coarser, e1 finer.
      {
         add_face(c2, f.face_id, FaceRole::ProviderCoarseSkip, e1);
         add_face(c1, f.face_id, FaceRole::ConsumerFine,       e2);
         is_coarse_side[static_cast<std::size_t>(e2)] = 1;
      }
   }

   // Derive the coarse-side element sets (ascending id) + dense slot maps.
   L.provider_slot_of_elem.assign(static_cast<std::size_t>(ne), -1);
   L.buffer_slot_of_elem.assign(static_cast<std::size_t>(ne), -1);
   for (int e = 0; e < ne; ++e)
   {
      if (!is_coarse_side[static_cast<std::size_t>(e)]) { continue; }
      L.provider_slot_of_elem[static_cast<std::size_t>(e)] =
         static_cast<int>(L.provider_elems.size());
      L.provider_elems.push_back(e);
      L.buffer_slot_of_elem[static_cast<std::size_t>(e)] =
         static_cast<int>(L.consumer_owner_elems.size());
      L.consumer_owner_elems.push_back(e);
   }

   return L;
}


LtsGlobalMeta ReduceGlobalMeta(
   const LtsLayout&                                    layout,
   const std::function<void(std::vector<long long>&)>& sum_across_ranks)
{
   const int nc = layout.num_clusters;
   MFEM_VERIFY(static_cast<int>(layout.local_elems_per_cluster.size()) == nc
               && static_cast<int>(layout.local_fault_faces_per_cluster.size())
                     == nc,
               "ReduceGlobalMeta: layout metadata size mismatch.");

   // Pack [elems(0..Nc), fault_faces(0..Nc)] and reduce in place.
   std::vector<long long> buf(static_cast<std::size_t>(2 * nc), 0);
   for (int c = 0; c < nc; ++c)
   {
      buf[static_cast<std::size_t>(c)]      = layout.local_elems_per_cluster[c];
      buf[static_cast<std::size_t>(nc + c)] = layout.local_fault_faces_per_cluster[c];
   }
   if (sum_across_ranks) { sum_across_ranks(buf); }
   MFEM_VERIFY(static_cast<int>(buf.size()) == 2 * nc,
               "ReduceGlobalMeta: sum_across_ranks resized the buffer.");

   LtsGlobalMeta meta;
   meta.global_elems.assign(static_cast<std::size_t>(nc), 0);
   meta.global_fault_faces.assign(static_cast<std::size_t>(nc), 0);
   for (int c = 0; c < nc; ++c)
   {
      meta.global_elems[static_cast<std::size_t>(c)]       = buf[static_cast<std::size_t>(c)];
      meta.global_fault_faces[static_cast<std::size_t>(c)] = buf[static_cast<std::size_t>(nc + c)];
   }
   return meta;
}

} // namespace seas
} // namespace mfem
