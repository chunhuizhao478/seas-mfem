// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/lts_partition.cpp — the METIS binding for the LTS-aware partition
// (Phase 4 Step 0, P-009).  This translation unit includes <metis.h> and
// NOTHING from mfem, so metis.h's `real_t` (float) does not collide with
// `mfem::real_t` (double).  See lts_partition.hpp for the contract.

#include "lts_partition.hpp"

#include <metis.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace mfem
{
namespace seas
{

std::vector<int> BuildLtsVertexWeights(const std::vector<int>& cluster,
                                       int num_clusters,
                                       const std::vector<double>* cell_cost,
                                       bool scalar_weights)
{
   const int nv = static_cast<int>(cluster.size());
   // Per-element integer work weight w_v >= 1 (METIS needs positive weights),
   // clamped to a ceiling so the scalar `2^shift` scaling stays in int range.
   auto weight_of = [&](int v) -> long long
   {
      long long w = 1;
      if (cell_cost && (*cell_cost)[v] > 0.0)
      {
         // REVIEW B-7: clamp the DOUBLE to the ceiling before llround — llround
         // on a value > LLONG_MAX is implementation-defined (may return a
         // negative sentinel, defeating the intended clamp).
         const double c = std::min((*cell_cost)[v],
                                   static_cast<double>(1LL << 20));
         w = std::llround(c);
      }
      if (w < 1) { w = 1; }
      if (w > (1LL << 20)) { w = (1LL << 20); }
      return w;
   };

   if (scalar_weights || num_clusters <= 1)
   {
      std::vector<int> vwgt(static_cast<std::size_t>(nv));
      const int maxC = (num_clusters > 0) ? num_clusters - 1 : 0;
      for (int v = 0; v < nv; ++v)
      {
         int shift = maxC - cluster[v];
         if (shift < 0)  { shift = 0; }
         if (shift > 30) { shift = 30; }
         long long w = weight_of(v) << shift;
         if (w > 2000000000LL) { w = 2000000000LL; }   // < INT_MAX
         vwgt[static_cast<std::size_t>(v)] = static_cast<int>(w);
      }
      return vwgt;
   }

   // Multi-constraint: each element contributes its weight to ONE constraint
   // (its cluster) and 0 to the others, so METIS balances every cluster's work
   // independently across ranks.
   std::vector<int> vwgt(static_cast<std::size_t>(nv) * num_clusters, 0);
   for (int v = 0; v < nv; ++v)
   {
      const int c = cluster[v];
      vwgt[static_cast<std::size_t>(v) * num_clusters + c] =
         static_cast<int>(weight_of(v));
   }
   return vwgt;
}

namespace
{
// Union-find over elements (path-halving + union-by-rank).
struct UnionFind
{
   std::vector<int> parent, rank;
   explicit UnionFind(int n) : parent(n), rank(n, 0)
   { for (int i = 0; i < n; ++i) { parent[i] = i; } }
   int Find(int x)
   { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
   void Union(int a, int b)
   {
      a = Find(a); b = Find(b);
      if (a == b) { return; }
      if (rank[a] < rank[b]) { std::swap(a, b); }
      parent[b] = a;
      if (rank[a] == rank[b]) { ++rank[a]; }
   }
};
}  // namespace

int VerifyLtsPartitionFaultLocality(
   const std::vector<std::pair<int, int>>& fault_pairs,
   const std::vector<int>& part)
{
   const int np = static_cast<int>(part.size());
   int violations = 0;
   for (const auto& pr : fault_pairs)
   {
      if (pr.first >= 0 && pr.second >= 0 && pr.first < np && pr.second < np
          && part[pr.first] != part[pr.second])
      {
         ++violations;
      }
   }
   return violations;
}

bool BuildLtsAwarePartition(
   int nvtxs,
   const std::vector<int>& xadj,
   const std::vector<int>& adjncy,
   const std::vector<int>& cluster,
   int num_clusters,
   int nparts,
   const std::vector<double>* cell_cost,
   const std::vector<std::pair<int, int>>& fault_pairs,
   const LtsPartitionOptions& opt,
   std::vector<int>& part_out)
{
   if (nvtxs <= 0) { part_out.clear(); return true; }
   part_out.assign(static_cast<std::size_t>(nvtxs), 0);
   if (nparts <= 1) { return true; }   // trivial: everything on rank 0

   if (static_cast<int>(cluster.size()) != nvtxs)
   { throw std::runtime_error("BuildLtsAwarePartition: cluster.size() != nvtxs"); }
   if (static_cast<int>(xadj.size()) != nvtxs + 1)
   { throw std::runtime_error("BuildLtsAwarePartition: xadj.size() != nvtxs+1"); }
   if (num_clusters < 1)
   { throw std::runtime_error("BuildLtsAwarePartition: num_clusters < 1"); }
   if (cell_cost && static_cast<int>(cell_cost->size()) != nvtxs)
   { throw std::runtime_error("BuildLtsAwarePartition: cell_cost.size() != nvtxs"); }

   // REVIEW B-2: every cluster id must be in [0, num_clusters); for the
   // multi-constraint partition each cluster must be NON-EMPTY (a zero-total
   // constraint is version-dependent in METIS balance refinement).
   // BuildLtsClustering guarantees contiguous 0-based ids, but validate here
   // since this is a public / independently-tested entry point.
   {
      std::vector<char> used(static_cast<std::size_t>(num_clusters), 0);
      for (int v = 0; v < nvtxs; ++v)
      {
         const int c = cluster[v];
         if (c < 0 || c >= num_clusters)
         {
            throw std::runtime_error(
               "BuildLtsAwarePartition: cluster id " + std::to_string(c)
               + " out of [0, " + std::to_string(num_clusters) + ")");
         }
         used[static_cast<std::size_t>(c)] = 1;
      }
      if (!opt.scalar_weights && num_clusters > 1)
      {
         for (int c = 0; c < num_clusters; ++c)
         {
            if (!used[static_cast<std::size_t>(c)])
            {
               throw std::runtime_error(
                  "BuildLtsAwarePartition: cluster " + std::to_string(c)
                  + " is empty (multi-constraint requires every cluster "
                    "non-empty; ids must be contiguous 0-based)");
            }
         }
      }
   }

   const bool multi = (!opt.scalar_weights) && (num_clusters > 1);
   idx_t ncon = multi ? static_cast<idx_t>(num_clusters) : 1;

   const std::vector<int> vwgt_i =
      BuildLtsVertexWeights(cluster, num_clusters, cell_cost, !multi);

   // Copy inputs into idx_t / real_t arrays (idx_t == int32 here; METIS wants
   // mutable pointers).  real_t is metis's float in this TU — no mfem collision.
   std::vector<idx_t> m_xadj(xadj.begin(), xadj.end());
   std::vector<idx_t> m_adj(adjncy.begin(), adjncy.end());
   std::vector<idx_t> m_vwgt(vwgt_i.begin(), vwgt_i.end());
   std::vector<idx_t> m_part(static_cast<std::size_t>(nvtxs), 0);
   idx_t m_nv = nvtxs, m_np = nparts, objval = 0;
   std::vector<real_t> ubvec(static_cast<std::size_t>(ncon),
                             static_cast<real_t>(opt.imbalance));

   idx_t options[METIS_NOPTIONS];
   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_SEED] = 0;   // deterministic: same mesh => same cut

   const int status = METIS_PartGraphKway(
      &m_nv, &ncon, m_xadj.data(), m_adj.data(), m_vwgt.data(),
      /*vsize=*/nullptr, /*adjwgt=*/nullptr, &m_np, /*tpwgts=*/nullptr,
      ubvec.data(), options, &objval, m_part.data());
   if (status != METIS_OK)
   {
      throw std::runtime_error(
         "BuildLtsAwarePartition: METIS_PartGraphKway failed (status "
         + std::to_string(status) + ")");
   }

   for (int v = 0; v < nvtxs; ++v)
   { part_out[static_cast<std::size_t>(v)] = static_cast<int>(m_part[v]); }

   // Fault locality (D-2): union-find each fault pair, then lock every component
   // to the rank of its lowest-index member (mirrors
   // BuildFaultLocalityPartitioning; no fault face is cut across ranks).
   if (!fault_pairs.empty())
   {
      UnionFind uf(nvtxs);
      for (const auto& pr : fault_pairs)
      {
         if (pr.first >= 0 && pr.second >= 0 && pr.first < nvtxs
             && pr.second < nvtxs)
         { uf.Union(pr.first, pr.second); }
      }
      std::vector<int> root_rank(static_cast<std::size_t>(nvtxs), -1);
      for (int e = 0; e < nvtxs; ++e)
      {
         const int r = uf.Find(e);
         if (root_rank[r] < 0) { root_rank[r] = part_out[e]; }  // lowest-index
      }
      for (int e = 0; e < nvtxs; ++e)
      { part_out[e] = root_rank[uf.Find(e)]; }
   }
   return true;
}

} // namespace seas
} // namespace mfem
