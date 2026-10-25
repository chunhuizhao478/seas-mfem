// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_partition.cpp — LTS Phase 4 Step 0 (P-009): the cluster-weighted
// METIS partition companion (dynamic/lts_partition.{hpp,cpp}).
//
// Unit-testable slice (the balance QUALITY at scale is Frontera-only):
//   (1) BuildLtsVertexWeights — multi-constraint one-hot-by-cluster weights and
//       the scalar 2^(maxC-cluster) fallback, incl. clamping;
//   (2) BuildLtsAwarePartition — on a small connected grid: partition in range
//       [0, nparts), the fault-locality invariant (every fault pair co-located),
//       the nparts<=1 trivial case, and the scalar-weight fallback;
//   (3) L2 (Chunhui Zhao 2026-10-25): BuildLtsBlockWeights against hand counts,
//       its rank sums against the patch elements + 2 x interior + shared face
//       blocks counted from a cut, ComputeLtsBlockBalance by hand, and the
//       block-weight cut (in range, fault pairs together, a better block
//       balance than the tet cut on a graded grid).

#include "../../dynamic/lts_partition.hpp"

#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace mfem::seas;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

namespace
{
// Build a P×Q structured grid's CSR element adjacency (4-neighbour, symmetric,
// self-loop-free) + fault pairs along the vertical seam at column `seam_col`.
void BuildGrid(int P, int Q, std::vector<int>& xadj, std::vector<int>& adjncy,
               std::vector<int>& cluster, std::vector<std::pair<int,int>>& fault,
               int seam_col)
{
   const int N = P * Q;
   auto id = [Q](int i, int j) { return i * Q + j; };
   xadj.assign(N + 1, 0);
   adjncy.clear();
   cluster.assign(N, 0);
   fault.clear();
   for (int i = 0; i < P; ++i)
   {
      for (int j = 0; j < Q; ++j)
      {
         const int v = id(i, j);
         // cluster 0 left of the seam, cluster 1 at/after it.
         cluster[v] = (j < seam_col) ? 0 : 1;
         if (i > 0)     { adjncy.push_back(id(i - 1, j)); }
         if (i < P - 1) { adjncy.push_back(id(i + 1, j)); }
         if (j > 0)     { adjncy.push_back(id(i, j - 1)); }
         if (j < Q - 1) { adjncy.push_back(id(i, j + 1)); }
         xadj[v + 1] = static_cast<int>(adjncy.size());
         // fault pair across the seam edge (horizontal neighbours j=seam_col-1,seam_col).
         if (j == seam_col && seam_col > 0)
         { fault.emplace_back(id(i, seam_col - 1), id(i, seam_col)); }
      }
   }
}
// A P×Q grid graded in columns: cluster c holds `width[c]` columns from the
// left (the finest at the fault seam on the left edge, fault pairs between
// columns 0 and 1), the coarsest the rest.
void BuildGradedGrid(int P, int Q, const std::vector<int>& width, std::vector<int>& xadj,
                     std::vector<int>& adjncy, std::vector<int>& cluster,
                     std::vector<std::pair<int,int>>& fault)
{
   BuildGrid(P, Q, xadj, adjncy, cluster, fault, 1);
   const int nc = static_cast<int>(width.size());
   for (int i = 0; i < P; ++i)
   {
      for (int j = 0, c = 0, edge = width[0]; j < Q; ++j)
      {
         while (j >= edge && c < nc - 1) { ++c; edge += width[c]; }
         cluster[i * Q + j] = c;
      }
   }
}

// P_c(r) + 2 I_c(r) + S_c(r) counted from the cut, as NitscheLtsBlocks lays
// out the blocks: flat [c * nparts + r].
std::vector<long long> CountBlocks(const std::vector<int>& xadj, const std::vector<int>& adjncy,
                                   const std::vector<int>& cluster, int nc,
                                   const std::vector<int>& part, int nparts)
{
   std::vector<long long> B(static_cast<std::size_t>(nc) * nparts, 0);
   const int nv = static_cast<int>(cluster.size());
   for (int v = 0; v < nv; ++v)
   {
      std::set<int> patch = {cluster[v]};
      for (int k = xadj[v]; k < xadj[v + 1]; ++k) { patch.insert(cluster[adjncy[k]]); }
      for (int c : patch) { B[static_cast<std::size_t>(c) * nparts + part[v]] += 1; }
      for (int k = xadj[v]; k < xadj[v + 1]; ++k)
      {
         const int f = adjncy[k];
         std::set<int> sides = {cluster[v], cluster[f]};
         for (int c : sides)
         {
            if (part[f] != part[v])
            { B[static_cast<std::size_t>(c) * nparts + part[v]] += 1; }   // shared
            else if (v < f)
            { B[static_cast<std::size_t>(c) * nparts + part[v]] += 2; }   // interior
         }
      }
   }
   return B;
}
}  // namespace

int main()
{
   // ---- T1: multi-constraint vertex weights (one-hot by cluster) ----
   {
      const std::vector<int> cluster = {0, 1, 2, 1, 0};
      const int nc = 3;
      const std::vector<int> w = BuildLtsVertexWeights(cluster, nc, nullptr, false);
      CHECK(w.size() == cluster.size() * nc, "multi vwgt size");
      bool ok = true;
      for (std::size_t v = 0; v < cluster.size(); ++v)
      {
         for (int c = 0; c < nc; ++c)
         {
            const int expect = (cluster[v] == c) ? 1 : 0;
            if (w[v * nc + c] != expect) { ok = false; }
         }
      }
      CHECK(ok, "multi vwgt is one-hot (w=1) at the element's cluster, 0 else");
   }

   // ---- T2: scalar fallback weights = 2^(maxC - cluster) ----
   {
      const std::vector<int> cluster = {0, 1, 2, 3};
      const int nc = 4;   // maxC = 3
      const std::vector<int> w = BuildLtsVertexWeights(cluster, nc, nullptr, true);
      CHECK(w.size() == cluster.size(), "scalar vwgt size");
      CHECK(w[0] == 8 && w[1] == 4 && w[2] == 2 && w[3] == 1,
            "scalar vwgt = 2^(maxC - cluster)");
   }

   // ---- T3: weight from cell_cost, with clamping ----
   {
      const std::vector<int> cluster = {0, 0, 0};
      const std::vector<double> cost = {3.6, -5.0, 1e30};   // round / floor / clamp
      const std::vector<int> w = BuildLtsVertexWeights(cluster, 1, &cost, true);
      CHECK(w[0] == 4, "cell_cost rounds (3.6 -> 4)");
      CHECK(w[1] == 1, "non-positive cell_cost clamps to 1");
      CHECK(w[2] >= 1 && w[2] <= 2000000000, "huge cell_cost clamps into int range");
   }

   // ---- T4: partition a 4×6 grid into 3 ranks; fault seam at col 3 ----
   {
      std::vector<int> xadj, adjncy, cluster;
      std::vector<std::pair<int,int>> fault;
      BuildGrid(/*P=*/4, /*Q=*/6, xadj, adjncy, cluster, fault, /*seam_col=*/3);
      const int N = 24;
      std::vector<int> part;
      LtsPartitionOptions opt;
      bool ok = BuildLtsAwarePartition(N, xadj, adjncy, cluster, /*nc=*/2,
                                       /*nparts=*/3, nullptr, fault, opt, part);
      CHECK(ok, "multi-constraint partition succeeded");
      CHECK(static_cast<int>(part.size()) == N, "partition length == N");
      bool in_range = true;
      for (int p : part) { if (p < 0 || p >= 3) { in_range = false; } }
      CHECK(in_range, "every partition id in [0, nparts)");
      CHECK(VerifyLtsPartitionFaultLocality(fault, part) == 0,
            "fault pairs are co-located (fault-locality invariant)");
      CHECK(!fault.empty(), "fixture actually has fault pairs (non-vacuous)");

      // B-6: the invariant above is vacuous if the partition collapsed to one
      // rank.  Assert it is non-trivial.
      std::set<int> ranks(part.begin(), part.end());
      CHECK(ranks.size() >= 2, "partition uses >1 rank (non-trivial)");
   }

   // ---- T4b: the fault-lock is DETERMINISTICALLY exercised.  Two triangles
   // (cliques {0,1,2} and {3,4,5}) bridged by the single edge (2,3), which is
   // ALSO the fault pair; the balanced 2-way min-cut MUST sever the bridge, so
   // the raw METIS partition cuts the fault pair and the lock must repair it. ----
   {
      const std::vector<int> xadj   = {0, 2, 4, 7, 10, 12, 14};
      const std::vector<int> adjncy = {1,2, 0,2, 0,1,3, 2,4,5, 3,5, 3,4};
      const std::vector<int> cluster(6, 0);
      const std::vector<std::pair<int,int>> fault = {{2, 3}};
      const std::vector<std::pair<int,int>> nofault;
      LtsPartitionOptions opt;

      std::vector<int> raw, locked;
      BuildLtsAwarePartition(6, xadj, adjncy, cluster, 1, 2, nullptr, nofault,
                             opt, raw);
      CHECK(VerifyLtsPartitionFaultLocality(fault, raw) > 0,
            "raw METIS severs the bridge fault pair (lock has work to do)");
      BuildLtsAwarePartition(6, xadj, adjncy, cluster, 1, 2, nullptr, fault,
                             opt, locked);
      CHECK(VerifyLtsPartitionFaultLocality(fault, locked) == 0,
            "the fault-lock repairs the cut pair (lock is exercised)");
   }

   // ---- T5: scalar-weight fallback also partitions + preserves fault locality --
   {
      std::vector<int> xadj, adjncy, cluster;
      std::vector<std::pair<int,int>> fault;
      BuildGrid(4, 6, xadj, adjncy, cluster, fault, 3);
      const int N = 24;
      std::vector<int> part;
      LtsPartitionOptions opt; opt.scalar_weights = true;
      bool ok = BuildLtsAwarePartition(N, xadj, adjncy, cluster, 2, 3, nullptr,
                                       fault, opt, part);
      CHECK(ok, "scalar-fallback partition succeeded");
      CHECK(VerifyLtsPartitionFaultLocality(fault, part) == 0,
            "scalar-fallback preserves fault locality");
   }

   // ---- T6: nparts <= 1 is the trivial all-rank-0 partition ----
   {
      std::vector<int> xadj = {0}, adjncy, cluster;
      std::vector<std::pair<int,int>> fault;
      std::vector<int> part;
      LtsPartitionOptions opt;
      // nvtxs=0 edge case
      CHECK(BuildLtsAwarePartition(0, xadj, adjncy, cluster, 1, 4, nullptr, fault,
                                   opt, part) && part.empty(),
            "nvtxs=0 => empty partition");
      // nparts=1
      std::vector<int> xa = {0, 0, 0, 0}, ad, cl = {0, 0, 0};
      std::vector<int> p1;
      CHECK(BuildLtsAwarePartition(3, xa, ad, cl, 1, 1, nullptr, fault, opt, p1)
            && p1.size() == 3 && p1[0] == 0 && p1[1] == 0 && p1[2] == 0,
            "nparts=1 => all rank 0");
   }

   // ---- T7: ComputeLtsPartitionImbalance (Phase-5 realization diagnostic) ----
   {
      // (a) perfectly balanced: each cluster split evenly across 2 ranks.
      {
         std::vector<int> cluster = {0, 0, 1, 1};
         std::vector<int> part    = {0, 1, 0, 1};
         std::vector<std::pair<int,int>> nofault;
         auto pib = ComputeLtsPartitionImbalance(cluster, 2, part, 2, nullptr,
                                                 nofault);
         CHECK(pib.per_cluster.size() == 2,
               "imb: per_cluster length == num_clusters");
         CHECK(std::abs(pib.overall - 1.0) < 1e-12,
               "imb: balanced => overall == 1");
         CHECK(std::abs(pib.per_cluster[0] - 1.0) < 1e-12,
               "imb: balanced => cluster0 == 1");
         CHECK(std::abs(pib.per_cluster[1] - 1.0) < 1e-12,
               "imb: balanced => cluster1 == 1");
         CHECK(std::abs(pib.fault - 1.0) < 1e-12,
               "imb: no fault pairs => fault == 1");
      }
      // (b) cluster 0 piled onto rank 0 => cluster-0 imbalance = 2, worst = c0.
      {
         std::vector<int> cluster = {0, 0, 1, 1};
         std::vector<int> part    = {0, 0, 0, 1};
         std::vector<std::pair<int,int>> nofault;
         auto pib = ComputeLtsPartitionImbalance(cluster, 2, part, 2, nullptr,
                                                 nofault);
         CHECK(std::abs(pib.per_cluster[0] - 2.0) < 1e-12,
               "imb: piled cluster0 imbalance == 2");
         CHECK(std::abs(pib.per_cluster[1] - 1.0) < 1e-12,
               "imb: cluster1 still balanced");
         CHECK(pib.worst_cluster == 0, "imb: worst cluster is c0");
         CHECK(std::abs(pib.worst_cluster_imbalance - 2.0) < 1e-12,
               "imb: worst imbalance == 2");
         CHECK(pib.overall > 1.0, "imb: overall > 1 when piled");
      }
      // (c) balanced BULK but fault work concentrated on one rank (SAFS gap).
      {
         std::vector<int> cluster = {0, 0, 0, 0};
         std::vector<int> part    = {0, 0, 1, 1};
         std::vector<std::pair<int,int>> fault = {{0, 1}};   // both on rank 0
         auto pib = ComputeLtsPartitionImbalance(cluster, 1, part, 2, nullptr,
                                                 fault);
         CHECK(std::abs(pib.per_cluster[0] - 1.0) < 1e-12,
               "imb: single-cluster bulk balanced");
         CHECK(std::abs(pib.fault - 2.0) < 1e-12,
               "imb: fault work all on one rank => fault == 2");
      }
      // (d) defensive: size mismatch / empty inputs => neutral 1.0.
      {
         std::vector<int> cluster = {0, 1};
         std::vector<int> part_bad = {0};        // wrong length
         std::vector<std::pair<int,int>> nofault;
         auto pib = ComputeLtsPartitionImbalance(cluster, 2, part_bad, 2, nullptr,
                                                 nofault);
         CHECK(std::abs(pib.overall - 1.0) < 1e-12,
               "imb: size mismatch => neutral overall 1");
         std::vector<int> empty;
         auto pib2 = ComputeLtsPartitionImbalance(empty, 0, empty, 2, nullptr,
                                                  nofault);
         CHECK(std::abs(pib2.overall - 1.0) < 1e-12,
               "imb: empty => neutral overall 1");
      }
   }

   // ---- T8 (L2): block weights against hand counts ----
   {
      // element 0 in cluster 1 with neighbours 1, 2 (cluster 1) and 3
      // (cluster 2), its fourth face on the boundary (the note's example);
      // element 4 in cluster 3 with four neighbours 5..8 in cluster 3
      const std::vector<int> xadj   = {0, 3, 4, 5, 6, 10, 11, 12, 13, 14};
      const std::vector<int> adjncy = {1, 2, 3, 0, 0, 0, 5, 6, 7, 8, 4, 4, 4, 4};
      const std::vector<int> cluster = {1, 1, 1, 2, 3, 3, 3, 3, 3};
      const int nc = 4;
      const std::vector<int> w = BuildLtsBlockWeights(xadj, adjncy, cluster, nc, nullptr);
      CHECK(w.size() == cluster.size() * nc, "blocks: size nvtxs * num_clusters");
      auto row = [&](int v) { return std::vector<int>(w.begin() + v * nc, w.begin() + (v + 1) * nc); };
      CHECK(row(0) == std::vector<int>({0, 4, 2, 0}), "blocks: seam tet w_1 = 1 + 3, w_2 = 1 + 1");
      CHECK(row(1) == std::vector<int>({0, 2, 0, 0}), "blocks: one face inside cluster 1 -> 2");
      CHECK(row(3) == std::vector<int>({0, 2, 2, 0}), "blocks: cluster-2 tet facing cluster 1");
      CHECK(row(4) == std::vector<int>({0, 0, 0, 5}), "blocks: deep tet = 1 + 4");
      // cell_cost: x 100, plus round(100 x) on the own cluster; <= 0 adds nothing
      std::vector<double> extra(cluster.size(), 0.0);
      extra[0] = 0.5;
      extra[3] = -2.0;
      const std::vector<int> wc = BuildLtsBlockWeights(xadj, adjncy, cluster, nc, &extra);
      CHECK(std::vector<int>(wc.begin(), wc.begin() + nc) == std::vector<int>({0, 450, 200, 0}),
            "blocks: cell_cost 0.5 -> 100 w + 50 on the own cluster");
      CHECK(std::vector<int>(wc.begin() + 3 * nc, wc.begin() + 4 * nc)
            == std::vector<int>({0, 200, 200, 0}), "blocks: negative cell_cost adds nothing");
      bool threw = false;
      try
      {
         const std::vector<int> bad = {1, 2, 3, 9, 0, 0, 5, 6, 7, 8, 4, 4, 4, 4};
         BuildLtsBlockWeights(xadj, bad, cluster, nc, nullptr);
      }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw, "blocks: a neighbour out of range throws");
   }

   // ---- T9 (L2): rank sums of the block weights = P + 2I + S of the cut ----
   {
      std::vector<int> xadj, adjncy, cluster;
      std::vector<std::pair<int,int>> fault;
      BuildGradedGrid(12, 16, {2, 3, 4, 7}, xadj, adjncy, cluster, fault);
      const int nc = 4, nparts = 5, nv = 12 * 16;
      const std::vector<int> w = BuildLtsBlockWeights(xadj, adjncy, cluster, nc, nullptr);
      bool ok = true;
      for (int cut = 0; cut < 2; ++cut)
      {
         std::vector<int> part(nv);
         for (int v = 0; v < nv; ++v)
         { part[v] = cut == 0 ? (v * 7 + v / 5) % nparts : ((v / 16) / 3 + (v % 16) / 6) % nparts; }
         const std::vector<long long> B = CountBlocks(xadj, adjncy, cluster, nc, part, nparts);
         for (int c = 0; c < nc; ++c)
         {
            for (int r = 0; r < nparts; ++r)
            {
               long long sum = 0;
               for (int v = 0; v < nv; ++v) { if (part[v] == r) { sum += w[v * nc + c]; } }
               if (sum != B[static_cast<std::size_t>(c) * nparts + r]) { ok = false; }
            }
         }
      }
      CHECK(ok, "blocks: per rank and cluster, sum w_c = patch + 2 interior + shared blocks");
   }

   // ---- T10 (L2): ComputeLtsBlockBalance by hand ----
   {
      // rank 0: W_0 = 2, W_1 = 1; rank 1: W_0 = 0, W_1 = 1.  M_0 = 2, M_1 = 3,
      // interval = M_0 + M_1 = 5; at the mean 1 + 2 = 3.
      const std::vector<int> w = {2, 0, 0, 1, 0, 1};
      const std::vector<int> part = {0, 1, 0};
      const LtsBlockBalance b = ComputeLtsBlockBalance(w, 2, part, 2);
      CHECK(b.per_cluster.size() == 2 && std::abs(b.per_cluster[0] - 2.0) < 1e-12
            && std::abs(b.per_cluster[1] - 1.0) < 1e-12, "balance: per cluster 2 and 1");
      CHECK(std::abs(b.tick_model - 5.0 / 3.0) < 1e-12, "balance: tick model 5/3");
      // one element per cluster on different ranks: M_0 = 1, M_1 = 1, at the
      // mean 0.5 + 1
      const LtsBlockBalance e = ComputeLtsBlockBalance({1, 0, 0, 1}, 2, {0, 1}, 2);
      CHECK(std::abs(e.per_cluster[0] - 2.0) < 1e-12 && std::abs(e.tick_model - 4.0 / 3.0) < 1e-12,
            "balance: one element per cluster on different ranks");
      bool threw = false;
      try { ComputeLtsBlockBalance(w, 2, {0, 2, 0}, 2); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw, "balance: a rank out of range throws");
      threw = false;
      try { ComputeLtsBlockBalance(w, 2, {0, 1}, 2); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw, "balance: weights of the wrong size throw");
   }

   // ---- T11 (L2): the block-weight cut ----
   {
      std::vector<int> xadj, adjncy, cluster;
      std::vector<std::pair<int,int>> fault;
      BuildGradedGrid(40, 40, {2, 4, 8, 26}, xadj, adjncy, cluster, fault);
      const int nc = 4, nparts = 8, nv = 1600;
      LtsPartitionOptions tets, blocks;
      blocks.block_weights = true;
      std::vector<int> pt, pb;
      BuildLtsAwarePartition(nv, xadj, adjncy, cluster, nc, nparts, nullptr, fault, tets, pt);
      CHECK(BuildLtsAwarePartition(nv, xadj, adjncy, cluster, nc, nparts, nullptr, fault,
                                   blocks, pb), "block cut succeeded");
      bool in_range = static_cast<int>(pb.size()) == nv;
      for (int r : pb) { if (r < 0 || r >= nparts) { in_range = false; } }
      CHECK(in_range, "block cut: every part id in [0, nparts)");
      CHECK(VerifyLtsPartitionFaultLocality(fault, pb) == 0, "block cut: fault pairs together");
      CHECK(std::set<int>(pb.begin(), pb.end()).size() == static_cast<std::size_t>(nparts),
            "block cut: every rank gets elements");
      CHECK(pb != pt, "block cut differs from the tet cut");
      const std::vector<int> w = BuildLtsBlockWeights(xadj, adjncy, cluster, nc, nullptr);
      const LtsBlockBalance bt = ComputeLtsBlockBalance(w, nc, pt, nparts);
      const LtsBlockBalance bb = ComputeLtsBlockBalance(w, nc, pb, nparts);
      double mt = 0.0, mb = 0.0;
      for (int c = 0; c < nc; ++c)
      {
         mt = std::max(mt, bt.per_cluster[c]);
         mb = std::max(mb, bb.per_cluster[c]);
      }
      std::printf("  T11 graded 40x40 grid, 8 parts: block work max/mean over clusters "
                  "tets %.3f, blocks %.3f; tick model tets %.3f, blocks %.3f\n",
                  mt, mb, bt.tick_model, bb.tick_model);
      CHECK(bb.tick_model < bt.tick_model, "block cut: lower tick model on block work");
      CHECK(mb < mt, "block cut: lower worst block imbalance per cluster");
      // with a cell_cost (the driver's --fault-weight term) the cut still runs
      std::vector<double> extra(nv, 0.0);
      for (const auto& pr : fault) { extra[pr.first] += 1.0; extra[pr.second] += 1.0; }
      std::vector<int> pf;
      CHECK(BuildLtsAwarePartition(nv, xadj, adjncy, cluster, nc, nparts, &extra, fault,
                                   blocks, pf)
            && VerifyLtsPartitionFaultLocality(fault, pf) == 0,
            "block cut with cell_cost: succeeded, fault pairs together");
      LtsPartitionOptions both = blocks;
      both.scalar_weights = true;
      bool threw = false;
      try { BuildLtsAwarePartition(nv, xadj, adjncy, cluster, nc, nparts, nullptr, fault, both, pf); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw, "block_weights with scalar_weights throws");
   }

   std::printf("test_lts_partition: %d checks, %d failures\n", g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}
