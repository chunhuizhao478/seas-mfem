// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_clustering.cpp — implementation of Phase 0 clustering (Appendix A.1 of
// PLAN_clustered_lts_ader_2026-07-18.md).

#include "lts_clustering.hpp"

#include "mfem.hpp"   // MFEM_ABORT, MFEM_VERIFY, mfem::Table

#include <algorithm>
#include <cmath>       // std::ldexp
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

namespace
{

// dt of cluster c for a given base step, = dt_base * 2^c, EXACT (ldexp on a
// power of two is exact).  Used identically for binning-edge comparison, cost,
// and the GAP-A1 assert so the three never disagree.
inline double cluster_dt(double dt_base, int c)
{
   return std::ldexp(dt_base, c);
}

// Rate-2 binning of one element against lambda-scaled edges.  INTEGER loop —
// never floor(log2(.)) — so the id is bit-reproducible (determinism
// requirement).  Element joins the largest c with dt_base * 2^c <= dt_e.
inline int bin_one(double dt_e_i, double dt_base, int max_clusters)
{
   int c = 0;
   double edge = dt_base;                 // = dt_base * 2^0
   while (2.0 * edge <= dt_e_i && c < max_clusters - 1)
   {
      edge *= 2.0;
      ++c;
   }
   return c;
}

// In-place SeisSol neighbour-difference fixpoint: adjacent clusters differ by
// <= 1; the two elements of a fault face differ by 0.  Monotone-decreasing
// (ids only lowered) clamp iterated to a fixed point.  Fault edges are handled
// in their own pass (diff = 0, forcing equality); the same edge appearing in
// the adjacency graph is skipped there.
void maxdiff_fixpoint(std::vector<int>&                       cluster,
                      const mfem::Table&                      elem_to_elem,
                      const std::set<std::pair<int, int>>&    fault_edges)
{
   const int ne = static_cast<int>(cluster.size());
   bool changed = true;
   while (changed)
   {
      changed = false;

      // Adjacency pass: neighbours differ by at most 1 (skip fault edges).
      for (int i = 0; i < ne; ++i)
      {
         const int  *nbrs = elem_to_elem.GetRow(i);
         const int   nnbr = elem_to_elem.RowSize(i);
         for (int k = 0; k < nnbr; ++k)
         {
            const int j = nbrs[k];
            if (j < 0 || j >= ne || j == i) { continue; }
            const auto key = (i < j) ? std::make_pair(i, j)
                                     : std::make_pair(j, i);
            if (fault_edges.count(key)) { continue; }   // handled below at diff 0
            // cluster[i] <= cluster[j] + 1  and symmetric.
            if (cluster[i] > cluster[j] + 1) { cluster[i] = cluster[j] + 1; changed = true; }
            if (cluster[j] > cluster[i] + 1) { cluster[j] = cluster[i] + 1; changed = true; }
         }
      }

      // Fault pass: both sides of a fault face share one cluster (diff = 0).
      for (const auto& e : fault_edges)
      {
         const int i = e.first, j = e.second;
         const int m = std::min(cluster[i], cluster[j]);
         if (cluster[i] != m) { cluster[i] = m; changed = true; }
         if (cluster[j] != m) { cluster[j] = m; changed = true; }
      }
   }
}

// modeled_cost = sum_e cellCost_e / (2^{c_e} * dt_base)  (Appendix A.1).
// Lower cost == fewer element-updates per unit simulated time.
double compute_cost(const std::vector<int>&    cluster,
                    double                     dt_base,
                    const std::vector<double>* cell_cost)
{
   double cost = 0.0;
   const double inv_base = 1.0 / dt_base;
   for (std::size_t i = 0; i < cluster.size(); ++i)
   {
      const double w = cell_cost ? (*cell_cost)[i] : 1.0;
      // 1 / (2^c * dt_base) = ldexp(inv_base, -c), exact.
      cost += w * std::ldexp(inv_base, -cluster[i]);
   }
   return cost;
}

// num_clusters = max id + 1.  maxdiff <= 1 with a non-empty cluster 0
// guarantees the used ids are contiguous 0..max (an empty middle cluster would
// require a jump of 2 across some adjacency edge, which the fixpoint forbids).
int count_clusters(const std::vector<int>& cluster)
{
   int mx = 0;
   for (int c : cluster) { mx = std::max(mx, c); }
   return mx + 1;
}

} // anonymous namespace


LtsClustering BuildLtsClustering(
   const std::vector<double>&              dt_e,
   const mfem::Table&                      elem_to_elem,
   const std::vector<std::pair<int, int>>& fault_face_elem_pairs,
   const LtsClusteringOptions&             opt)
{
   const int ne = static_cast<int>(dt_e.size());
   MFEM_VERIFY(ne > 0, "BuildLtsClustering: dt_e is empty.");
   MFEM_VERIFY(opt.rate == 2,
               "BuildLtsClustering: only rate == 2 is supported in v1 (got "
               << opt.rate << ").");
   MFEM_VERIFY(opt.max_clusters >= 1,
               "BuildLtsClustering: max_clusters must be >= 1.");
   MFEM_VERIFY(elem_to_elem.Size() == ne,
               "BuildLtsClustering: elem_to_elem has " << elem_to_elem.Size()
               << " rows but dt_e has " << ne << " elements.");
   if (opt.cell_cost)
   {
      MFEM_VERIFY(static_cast<int>(opt.cell_cost->size()) == ne,
                  "BuildLtsClustering: cell_cost length "
                  << opt.cell_cost->size() << " != dt_e length " << ne << ".");
   }

   double dt_min = std::numeric_limits<double>::infinity();
   for (int i = 0; i < ne; ++i)
   {
      MFEM_VERIFY(dt_e[i] > 0.0 && std::isfinite(dt_e[i]),
                  "BuildLtsClustering: dt_e[" << i << "] = " << dt_e[i]
                  << " is not a positive finite number.");
      dt_min = std::min(dt_min, dt_e[i]);
   }

   // Fault edges as an unordered-pair set (i<j) for O(log) lookup/enumeration.
   std::set<std::pair<int, int>> fault_edges;
   for (const auto& p : fault_face_elem_pairs)
   {
      const int i = p.first, j = p.second;
      MFEM_VERIFY(i >= 0 && i < ne && j >= 0 && j < ne,
                  "BuildLtsClustering: fault pair (" << i << "," << j
                  << ") out of range [0," << ne << ").");
      if (i != j) { fault_edges.insert(i < j ? std::make_pair(i, j)
                                             : std::make_pair(j, i)); }
   }

   // ---- lambda selection: fixed, or a (0.5, 1] grid scan minimising cost. ----
   auto build_at_lambda = [&](double lambda, std::vector<int>& cluster) -> double
   {
      const double dt_base = lambda * dt_min;
      cluster.assign(ne, 0);
      for (int i = 0; i < ne; ++i)
      {
         cluster[i] = bin_one(dt_e[i], dt_base, opt.max_clusters);
      }
      maxdiff_fixpoint(cluster, elem_to_elem, fault_edges);
      return compute_cost(cluster, dt_base, opt.cell_cost);
   };

   double best_lambda = opt.wiggle_scan ? 1.0 : opt.lambda_fixed;
   std::vector<int> cluster;
   if (opt.wiggle_scan)
   {
      // Scan lambda from 1.00 down to just above 0.50 in 0.01 steps; keep the
      // STRICT minimum cost, so on a tie the largest lambda (encountered first)
      // wins — deterministic.
      double best_cost = std::numeric_limits<double>::infinity();
      std::vector<int> tmp;
      for (int step = 0; step <= 49; ++step)      // lambda = 1.00, 0.99, ..., 0.51
      {
         const double lambda = 1.0 - 0.01 * step;
         const double cost = build_at_lambda(lambda, tmp);
         if (cost < best_cost)
         {
            best_cost   = cost;
            best_lambda = lambda;
            cluster     = tmp;
         }
      }
   }
   else
   {
      MFEM_VERIFY(opt.lambda_fixed > 0.5 && opt.lambda_fixed <= 1.0,
                  "BuildLtsClustering: lambda_fixed must be in (0.5, 1]; got "
                  << opt.lambda_fixed);
      build_at_lambda(opt.lambda_fixed, cluster);
   }

   const double dt_base = best_lambda * dt_min;
   double cost = compute_cost(cluster, dt_base, opt.cell_cost);

   // ---- auto-merge: collapse the top level toward nc_cap, cost-gated. ----
   // baseline = the uncapped (chosen-lambda) cost; each merge must keep the
   // cost within (1 + merge_loss_tol) * baseline.  Merging reassigns cluster
   // (Nc-1) -> (Nc-2): CFL-safe (ids only decrease) and maxdiff-preserving
   // (the whole top level collapses onto the level below).
   if (opt.nc_cap > 0)
   {
      const double baseline = cost;
      int nc = count_clusters(cluster);
      while (nc > opt.nc_cap && nc > 1)
      {
         std::vector<int> merged = cluster;
         const int top = nc - 1;
         for (int i = 0; i < ne; ++i)
         {
            if (merged[i] == top) { merged[i] = top - 1; }
         }
         const double merged_cost = compute_cost(merged, dt_base, opt.cell_cost);
         if (merged_cost <= (1.0 + opt.merge_loss_tol) * baseline)
         {
            cluster = std::move(merged);
            cost    = merged_cost;
            --nc;
         }
         else
         {
            break;   // next merge would blow the cost budget; stop above nc_cap
         }
      }
   }

   LtsClustering out;
   out.cluster      = std::move(cluster);
   out.num_clusters = count_clusters(out.cluster);
   out.dt_base      = dt_base;
   out.lambda       = best_lambda;
   out.modeled_cost = cost;

   // GAP-A1 (production-path assert): every element steps no faster than its
   // own CFL limit, using the identical cluster_dt expression as binning.
   for (int i = 0; i < ne; ++i)
   {
      MFEM_VERIFY(cluster_dt(out.dt_base, out.cluster[i]) <= dt_e[i],
                  "BuildLtsClustering: element " << i << " assigned cluster "
                  << out.cluster[i] << " (dt = "
                  << cluster_dt(out.dt_base, out.cluster[i])
                  << ") exceeds its CFL limit dt_e = " << dt_e[i]
                  << " — clustering is CFL-unstable (should be impossible).");
   }
   // Contiguity sanity (maxdiff <= 1 with non-empty cluster 0 guarantees it).
   {
      std::vector<char> seen(out.num_clusters, 0);
      for (int c : out.cluster) { seen[c] = 1; }
      for (int c = 0; c < out.num_clusters; ++c)
      {
         MFEM_VERIFY(seen[c],
                     "BuildLtsClustering: cluster id " << c
                     << " is unused — ids are not contiguous (maxdiff bug).");
      }
   }
   return out;
}


std::vector<long long> LtsClustering::cells_per_cluster() const
{
   std::vector<long long> h(num_clusters, 0);
   for (int c : cluster) { ++h[c]; }
   return h;
}


double ArithmeticMeanSpeedup(const LtsClustering& cl)
{
   // (1/N) * sum_c cells_c * 2^c  = SeisSol's printed "clustered LTS" statistic.
   const auto h = cl.cells_per_cluster();
   long long N = 0;
   double num = 0.0;
   for (int c = 0; c < cl.num_clusters; ++c)
   {
      N   += h[c];
      num += static_cast<double>(h[c]) * std::ldexp(1.0, c);   // h_c * 2^c
   }
   return (N > 0) ? num / static_cast<double>(N) : 0.0;
}


double HarmonicUpdateSpeedup(const LtsClustering& cl)
{
   // N / sum_c cells_c / 2^c  = the honest GTS-vs-LTS element-update ratio.
   const auto h = cl.cells_per_cluster();
   long long N = 0;
   double den = 0.0;
   for (int c = 0; c < cl.num_clusters; ++c)
   {
      N   += h[c];
      den += static_cast<double>(h[c]) * std::ldexp(1.0, -c);  // h_c / 2^c
   }
   return (den > 0.0) ? static_cast<double>(N) / den : 0.0;
}

} // namespace seas
} // namespace mfem
