// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for BuildLtsClustering — Phase 0 / Appendix B.1 of
// PLAN_clustered_lts_ader_2026-07-18.md.  Pure, single-process, no HDF5.
//
// Cases: integer-loop bin edges (dt at a lower edge joins THAT cluster,
// dt_cluster == dt_e; one ulp below joins c-1); single-element mesh;
// all-equal-dt -> 1 cluster; maxdiff fixpoint on a chain; fault diff == 0;
// lambda-scan (a fat bin just above an edge makes the scan pick lambda < 1
// with strictly lower cost; CFL holds; RAW mode == the lambda=1 baseline);
// auto-merge (num_clusters <= nc_cap, cost within tolerance, ids contiguous,
// maxdiff preserved, nc_cap <= 0 reproduces the uncapped layout bit-for-bit).

#include "mfem.hpp"

#include "../../dynamic/lts_clustering.hpp"

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define CHECK(cond, msg) \
   do { ++num_tests; \
      if (!(cond)) { std::cerr << "FAILED: " << msg << " (line " << __LINE__ \
                               << ")\n"; ++num_failed; } \
      else { ++num_passed; } } while (0)

namespace
{

// Build an mfem::Table (symmetric element adjacency) from an adjacency list.
mfem::Table make_table(const std::vector<std::vector<int>>& adj)
{
   const int n = static_cast<int>(adj.size());
   mfem::Table t;
   t.MakeI(n);
   for (int i = 0; i < n; ++i) { for (int j : adj[i]) { (void)j; t.AddAColumnInRow(i); } }
   t.MakeJ();
   for (int i = 0; i < n; ++i) { for (int j : adj[i]) { t.AddConnection(i, j); } }
   t.ShiftUpI();
   return t;
}

// Undirected chain 0-1-2-...-(n-1).
std::vector<std::vector<int>> chain_adj(int n)
{
   std::vector<std::vector<int>> adj(n);
   for (int i = 0; i + 1 < n; ++i) { adj[i].push_back(i + 1); adj[i + 1].push_back(i); }
   return adj;
}

// Max neighbour-cluster difference over the graph, ignoring fault edges.
int max_neighbor_diff(const std::vector<int>& cl,
                      const std::vector<std::vector<int>>& adj)
{
   int mx = 0;
   for (std::size_t i = 0; i < adj.size(); ++i)
      for (int j : adj[i]) { mx = std::max(mx, std::abs(cl[i] - cl[j])); }
   return mx;
}

double cluster_dt(double dt_base, int c) { return std::ldexp(dt_base, c); }

} // namespace


int main()
{
   // ---- T1: integer-loop bin edges (lower edge inclusive; ulp below -> c-1) --
   {
      std::vector<double> dt = {1.0, 2.0, 4.0, 8.0};      // dt_min = 1 => dt_base = 1
      auto adj = chain_adj(4);
      LtsClusteringOptions opt; opt.wiggle_scan = false; opt.lambda_fixed = 1.0; opt.nc_cap = 0;
      auto r = BuildLtsClustering(dt, make_table(adj), {}, opt);
      CHECK(r.dt_base == 1.0, "T1 dt_base == dt_min");
      CHECK((r.cluster == std::vector<int>{0, 1, 2, 3}), "T1 bins 0,1,2,3");
      for (int i = 0; i < 4; ++i)
      { CHECK(cluster_dt(r.dt_base, r.cluster[i]) == dt[i], "T1 dt_cluster == dt_e at exact edge"); }

      // one ulp below 8.0 -> cluster 2 (not 3)
      std::vector<double> dt2 = {1.0, 2.0, 4.0, std::nextafter(8.0, 0.0)};
      auto r2 = BuildLtsClustering(dt2, make_table(adj), {}, opt);
      CHECK(r2.cluster[3] == 2, "T1 ulp below the edge joins the lower cluster");
   }

   // ---- T2: single element -> exactly one cluster --------------------------
   {
      std::vector<double> dt = {5.0};
      auto r = BuildLtsClustering(dt, make_table({{}}), {});
      CHECK(r.num_clusters == 1 && r.cluster == std::vector<int>{0}, "T2 single element");
      CHECK(r.dt_base == 5.0, "T2 dt_base == the lone dt");
   }

   // ---- T3: all-equal dt -> one cluster ------------------------------------
   {
      std::vector<double> dt = {3.0, 3.0, 3.0, 3.0};
      LtsClusteringOptions opt; opt.wiggle_scan = false; opt.lambda_fixed = 1.0; opt.nc_cap = 0;
      auto r = BuildLtsClustering(dt, make_table(chain_adj(4)), {}, opt);
      CHECK(r.num_clusters == 1, "T3 all-equal -> 1 cluster");
   }

   // ---- T4: maxdiff fixpoint climbs by 1 per hop from a fine element -------
   {
      std::vector<double> dt = {1.0, 100.0, 100.0, 100.0};   // raw bins 0,6,6,6
      auto adj = chain_adj(4);
      LtsClusteringOptions opt; opt.wiggle_scan = false; opt.lambda_fixed = 1.0; opt.nc_cap = 0;
      auto r = BuildLtsClustering(dt, make_table(adj), {}, opt);
      CHECK((r.cluster == std::vector<int>{0, 1, 2, 3}), "T4 fixpoint 0,1,2,3");
      CHECK(max_neighbor_diff(r.cluster, adj) <= 1, "T4 maxdiff <= 1");
      for (int i = 0; i < 4; ++i)
      { CHECK(cluster_dt(r.dt_base, r.cluster[i]) <= dt[i], "T4 CFL-safe"); }
   }

   // ---- T5: a fault face forces its two elements into one cluster ----------
   {
      std::vector<double> dt = {1.0, 64.0};                  // raw bins 0, 6
      auto adj = chain_adj(2);
      LtsClusteringOptions opt; opt.wiggle_scan = false; opt.lambda_fixed = 1.0; opt.nc_cap = 0;
      auto with_fault = BuildLtsClustering(dt, make_table(adj), {{0, 1}}, opt);
      CHECK((with_fault.cluster == std::vector<int>{0, 0}), "T5 fault pair -> same cluster");
      auto no_fault = BuildLtsClustering(dt, make_table(adj), {}, opt);
      CHECK((no_fault.cluster == std::vector<int>{0, 1}), "T5 no fault -> maxdiff-1 gradient");
   }

   // ---- T6: lambda-scan picks lambda < 1 when a fat bin sits just above an edge
   {
      // Chain: [dt=1] then 1000 elements at dt=3.9 (just below 4 = 2^2).
      const int nfat = 1000;
      std::vector<double> dt; dt.push_back(1.0);
      for (int i = 0; i < nfat; ++i) { dt.push_back(3.9); }
      auto adj = chain_adj(static_cast<int>(dt.size()));
      auto tab = make_table(adj);

      LtsClusteringOptions raw; raw.wiggle_scan = false; raw.lambda_fixed = 1.0; raw.nc_cap = 0;
      auto r_raw = BuildLtsClustering(dt, tab, {}, raw);       // lambda = 1 baseline

      LtsClusteringOptions scan; scan.wiggle_scan = true; scan.nc_cap = 0;
      auto r_scan = BuildLtsClustering(dt, tab, {}, scan);

      CHECK(r_scan.lambda < 1.0, "T6 scan picks lambda < 1");
      CHECK(r_scan.modeled_cost < r_raw.modeled_cost, "T6 scan cost strictly lower than lambda=1");
      CHECK(max_neighbor_diff(r_scan.cluster, adj) <= 1, "T6 maxdiff preserved after scan");
      for (std::size_t i = 0; i < dt.size(); ++i)
      { CHECK(cluster_dt(r_scan.dt_base, r_scan.cluster[i]) <= dt[i], "T6 CFL-safe at chosen lambda"); }

      // RAW mode (wiggle off, no merge) == the lambda=1 plain algorithm.
      CHECK(r_raw.lambda == 1.0, "T6 RAW mode lambda == 1");

      // CFL holds for several explicit lambda in the grid (each build asserts internally).
      for (double lam : {1.0, 0.97, 0.51})
      {
         LtsClusteringOptions o; o.wiggle_scan = false; o.lambda_fixed = lam; o.nc_cap = 0;
         auto rr = BuildLtsClustering(dt, tab, {}, o);
         for (std::size_t i = 0; i < dt.size(); ++i)
         { CHECK(cluster_dt(rr.dt_base, rr.cluster[i]) <= dt[i], "T6 CFL-safe per lambda"); }
      }
   }

   // ---- T7: auto-merge collapses the top toward nc_cap, cost-gated ---------
   {
      // Monotone chain: blocks at dt = 1,2,4,8 (1000 each) then 16,32,64,128 (10 each).
      std::vector<double> dt;
      auto add_block = [&](double v, int n) { for (int i = 0; i < n; ++i) { dt.push_back(v); } };
      add_block(1, 1000); add_block(2, 1000); add_block(4, 1000); add_block(8, 1000);
      add_block(16, 10); add_block(32, 10); add_block(64, 10); add_block(128, 10);
      auto adj = chain_adj(static_cast<int>(dt.size()));
      auto tab = make_table(adj);

      LtsClusteringOptions uncapped; uncapped.wiggle_scan = false; uncapped.lambda_fixed = 1.0; uncapped.nc_cap = 0;
      auto r_unc = BuildLtsClustering(dt, tab, {}, uncapped);
      CHECK(r_unc.num_clusters == 8, "T7 uncapped -> 8 clusters");
      const double baseline = r_unc.modeled_cost;

      LtsClusteringOptions capped; capped.wiggle_scan = false; capped.lambda_fixed = 1.0;
      capped.nc_cap = 4; capped.merge_loss_tol = 0.05;
      auto r_cap = BuildLtsClustering(dt, tab, {}, capped);
      CHECK(r_cap.num_clusters <= 4, "T7 merged num_clusters <= nc_cap");
      CHECK(r_cap.num_clusters == 4, "T7 tolerance permits reaching nc_cap");
      CHECK(max_neighbor_diff(r_cap.cluster, adj) <= 1, "T7 maxdiff preserved post-merge");
      CHECK(r_cap.modeled_cost <= (1.0 + 0.05) * baseline, "T7 cost within (1+tol)*baseline");
      // ids contiguous
      {
         std::vector<char> seen(r_cap.num_clusters, 0);
         for (int c : r_cap.cluster) { seen[c] = 1; }
         bool contig = true;
         for (char s : seen) { contig = contig && s; }
         CHECK(contig, "T7 ids contiguous");
      }
      // nc_cap <= 0 reproduces the uncapped layout bit-for-bit.
      CHECK(r_cap.dt_base == r_unc.dt_base, "T7 merge does not change dt_base");
      LtsClusteringOptions raw2 = uncapped; raw2.nc_cap = -1;
      auto r_raw2 = BuildLtsClustering(dt, tab, {}, raw2);
      CHECK(r_raw2.cluster == r_unc.cluster, "T7 nc_cap<=0 reproduces uncapped bit-for-bit");
   }

   // ---- T8 (REVIEW C-2): the COST GATE binds before nc_cap, and the gate is
   //      cumulative-vs-baseline, not per-step.  A 4-cluster chain with a fine
   //      cluster 0 dominating the cost: 200 cells dt=1 (c0), 8 dt=2 (c1),
   //      8 dt=4 (c2), 8 dt=8 (c3).  dt_base=1, so
   //        baseline B = 200 + 8/2 + 8/4 + 8/8 = 207.
   //      Merging the top level doubles those cells' cost contribution:
   //        merge c3->c2 : +8/8 = +1  -> cost 208
   //        merge c2->c1 : +16/4 = +4 -> cost 212
   //      With tol = 0.01 (budget = 1.01*207 = 209.07): merge 1 is accepted
   //      (208 <= 209.07) but merge 2 is rejected against the FIXED baseline
   //      (212 > 209.07) even though nc_cap = 1 asks to collapse further.  So the
   //      cost gate stops at Nc = 3.  A PER-STEP gate (compare to the previous
   //      cost 208) would give 212/208 = 1.019 > 1.01 and ALSO reject here — but
   //      the invariant that DISTINGUISHES them is that the FINAL cost must stay
   //      <= (1+tol)*BASELINE, which a per-step gate can violate on a longer
   //      escalating chain; we assert that invariant plus the early stop.
   {
      std::vector<double> dt;
      std::vector<std::vector<int>> adj;
      auto push_block = [&](int count, double dtv)
      { for (int k = 0; k < count; ++k) { dt.push_back(dtv); } };
      push_block(200, 1.0); push_block(8, 2.0); push_block(8, 4.0); push_block(8, 8.0);
      const int n = static_cast<int>(dt.size());   // 224
      adj = chain_adj(n);                            // connected, consecutive clusters
      mfem::Table tab = make_table(adj);

      LtsClusteringOptions unc; unc.wiggle_scan = false; unc.lambda_fixed = 1.0; unc.nc_cap = 0;
      auto r_unc = BuildLtsClustering(dt, tab, {}, unc);
      CHECK(r_unc.num_clusters == 4, "T8 uncapped -> 4 clusters");
      CHECK(std::abs(r_unc.modeled_cost - 207.0) < 1e-9, "T8 baseline cost == 207");

      LtsClusteringOptions capped = unc;
      capped.nc_cap = 1; capped.merge_loss_tol = 0.01;   // ask to collapse to 1
      auto r_cap = BuildLtsClustering(dt, tab, {}, capped);
      // Cost gate stops the merge at Nc=3, NOT the nc_cap=1 target.
      CHECK(r_cap.num_clusters == 3, "T8 cost gate binds before nc_cap (Nc=3, not 1)");
      CHECK(r_cap.num_clusters > capped.nc_cap, "T8 final Nc exceeds nc_cap (cost-limited)");
      // The load-bearing invariant: final cost within (1+tol)*BASELINE (a per-step
      // gate can silently overshoot this on an escalating chain).
      CHECK(r_cap.modeled_cost <= (1.0 + capped.merge_loss_tol) * r_unc.modeled_cost,
            "T8 final cost within (1+tol)*baseline");
      CHECK(std::abs(r_cap.modeled_cost - 208.0) < 1e-9, "T8 exactly one merge accepted (cost 208)");
      CHECK(max_neighbor_diff(r_cap.cluster, adj) <= 1, "T8 maxdiff preserved");
   }

   // ---- speedup statistics on the SeisSol benchmark histogram -------------
   {
      // Reconstruct a layout with the known SeisSol cell counts and check the
      // two speedup helpers reproduce the published/verified numbers.
      const long long cells[11] =
      { 123, 822, 3346, 13917, 171472, 487273, 192506, 87246, 312888, 49565, 136 };
      LtsClustering cl; cl.num_clusters = 11; cl.dt_base = 1.0;
      for (int c = 0; c < 11; ++c) { for (long long k = 0; k < cells[c]; ++k) { cl.cluster.push_back(c); } }
      const double harm = HarmonicUpdateSpeedup(cl);
      const double arith = ArithmeticMeanSpeedup(cl);
      CHECK(std::abs(harm - 38.73) < 0.05, "harmonic speedup ~= 38.73 (SeisSol mesh)");
      CHECK(std::abs(arith - 111.85) < 0.05, "arithmetic-mean speedup ~= 111.85 (SeisSol log)");
   }

   std::cout << "\ntest_lts_clustering: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed.\n";
   return num_failed == 0 ? 0 : 1;
}
