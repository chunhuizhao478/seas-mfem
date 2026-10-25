// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/lts_partition.hpp — LTS-aware cluster-weighted graph partition
// (Phase 4 Step 0, P-009 of PLAN_clustered_lts_ader_2026-07-18.md).
//
// Deliberately free of <metis.h> AND <mfem.hpp>: metis.h typedefs `real_t` as
// `float`, which collides with `mfem::real_t` (double).  The actual METIS
// binding lives in lts_partition.cpp, which includes <metis.h> and NOTHING from
// mfem.  Callers (the driver) build the CSR element-adjacency from the mesh on
// the mfem side and pass raw int arrays through this interface.
//
// Goal: when a fault-embedded regional mesh is partitioned for MPI, the fine
// clusters (near the fault, tiny elements) do FAR more work per rank than the
// coarse far-field clusters.  A single-constraint partition that balances total
// element COUNT leaves the fine-cluster work wildly imbalanced under LTS.  This
// partition balances EACH cluster's work independently (multi-constraint METIS,
// ncon = num_clusters), then locks both elements of every fault face onto one
// rank (fault locality, D-2) so no fault face is cut across ranks.

#ifndef MFEM_SEAS_LTS_PARTITION_HPP
#define MFEM_SEAS_LTS_PARTITION_HPP

#include <cstddef>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

struct LtsPartitionOptions
{
   /// METIS `ubvec` — the allowed load imbalance per constraint (1.05 = 5%).
   double imbalance = 1.05;
   /// When true (or on METIS multi-constraint failure), fall back to a SINGLE
   /// constraint with geometric per-cluster weights `2^(maxC - cluster[v])`,
   /// so a finer element (smaller cluster id ⇒ more steps) carries more weight.
   bool scalar_weights = false;
   /// L2 (Chunhui Zhao 2026-10-25): when true, the constraint weights are the
   /// block work of BuildLtsBlockWeights (what the LTS sweep of each cluster
   /// pays) instead of the tets of BuildLtsVertexWeights, and `cell_cost` is
   /// read as BuildLtsBlockWeights reads it.  Excludes `scalar_weights`.
   bool block_weights = false;
};

/// Build the METIS vertex weights (exposed for unit testing).
///   multi-constraint (`scalar_weights == false`, the default):
///     `vwgt[v*num_clusters + c] = w_v` if `cluster[v] == c` else `0`,
///   single-constraint (`scalar_weights == true`):
///     `vwgt[v] = w_v * 2^(maxC - cluster[v])`,
/// where `w_v = max(1, round(cell_cost[v]))` (or `1` when `cell_cost == nullptr`)
/// and `maxC = num_clusters - 1`.  Weights are clamped to `>= 1` (METIS requires
/// positive integer weights) and to a safe ceiling so the `2^(...)` scale cannot
/// overflow int.  Returns a flat vector of length `nvtxs*num_clusters`
/// (multi-constraint) or `nvtxs` (scalar).
std::vector<int> BuildLtsVertexWeights(
   const std::vector<int>& cluster,
   int num_clusters,
   const std::vector<double>* cell_cost,
   bool scalar_weights);

/// Block-work weights (L2): constraint c of element e gets the dense blocks of
/// cluster c's LTS sweep (solver/nitsche_lts_run.hpp, NitscheLtsBlocks) that e
/// carries,
///   w_c(e) = [e in the patch of c] + #{face neighbours f of e :
///                                       cluster[e] == c or cluster[f] == c},
/// the patch of c being its elements and their face neighbours.  The first
/// term is e's diagonal block, the second gives every face with a side in c
/// one unit on each of its two elements, so a face inside a rank counts 2
/// (its block and its transpose) and a cut face 1 on each rank: the rank sums
/// of w_c are the patch elements + 2 x interior face blocks + shared face
/// blocks of cluster c for ANY cut.  Values 0..5 on tets (faces without a
/// neighbour carry no block).
///   xadj / adjncy : CSR face adjacency, as for BuildLtsAwarePartition.
///   cell_cost     : nullptr => the weights above.  Otherwise the extra work
///                   of each element in units of one block (e.g. the fault-
///                   face term), added on its own cluster, and everything is
///                   scaled by 100 so a fractional extra survives rounding:
///                   100 w_c(e) + [c == cluster[e]] round(100 cell_cost[e])
///                   (cell_cost <= 0 adds nothing; capped at 2^20).
/// Returns a flat vector of length `nvtxs*num_clusters`.  Throws
/// std::runtime_error on inconsistent sizes or ids.
std::vector<int> BuildLtsBlockWeights(
   const std::vector<int>& xadj,
   const std::vector<int>& adjncy,
   const std::vector<int>& cluster,
   int num_clusters,
   const std::vector<double>* cell_cost);

/// Partition `nvtxs` elements into `nparts` ranks with cluster-balanced weights,
/// then lock every fault pair onto one rank (fault locality).
///
///   xadj / adjncy : CSR element adjacency — `xadj` length `nvtxs+1`, `adjncy`
///                   length `xadj[nvtxs]` (each element's face-neighbours; the
///                   graph must be symmetric and self-loop-free, as METIS
///                   requires).
///   cluster       : per-element cluster id in `[0, num_clusters)`.
///   cell_cost     : per-element work weight (`nullptr` ⇒ uniform 1); with
///                   `opt.block_weights` the extra work per element of
///                   BuildLtsBlockWeights.
///   fault_pairs   : `(e1, e2)` element pairs across each fault face; both are
///                   forced onto the rep's rank AFTER the METIS cut (mirrors
///                   BuildFaultLocalityPartitioning).
///   part_out      : length-`nvtxs` partition in `[0, nparts)`.
///
/// Returns true on success (METIS + fault-lock).  `nparts <= 1` is the trivial
/// all-rank-0 partition (returns true).  On a METIS error the function throws
/// `std::runtime_error` (the caller may catch and fall back to a scalar-weight
/// or MFEM partition).
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
   std::vector<int>& part_out);

/// Verify the fault-locality invariant on a partition: every fault pair's two
/// elements share a rank.  Returns the number of violating pairs (0 = success).
int VerifyLtsPartitionFaultLocality(
   const std::vector<std::pair<int, int>>& fault_pairs,
   const std::vector<int>& part);

/// Partition-quality diagnostic (Phase 5 realization-fraction attribution).
///
/// For a computed `part` over `nparts` ranks, reports the load imbalance
/// `max_rank_work / mean_rank_work` (1.0 = perfect; the amount above 1.0 is the
/// fraction of the busiest rank the others sit idle) at three granularities:
///
///   overall     : LTS-weighted per-sync work `w_v * 2^(maxC - cluster_v)`
///                 (cluster c steps 2^(maxC-c) times per sync).  The aggregate
///                 idle-time predictor.
///   per_cluster : cluster-c work in isolation — exactly what the multi-
///                 constraint METIS targets, and the DIRECT per-substep-balance
///                 predictor.  The fine clusters (c=0,1) dominate wall-clock, so
///                 their imbalance is the number that decides the realization
///                 fraction.
///   fault       : count of fault-incident elements per rank (a friction-work
///                 proxy).  The bulk multi-constraint does NOT model friction
///                 work, so a high `fault` with low `per_cluster` flags the case
///                 where a fault-owning rank is friction-bound while its bulk
///                 clusters look balanced (the SAFS ncon=2 bulk+friction gap).
struct LtsPartitionImbalance
{
   double overall = 1.0;
   std::vector<double> per_cluster;    ///< length num_clusters
   double fault = 1.0;                 ///< 1.0 when there are no fault pairs
   int    worst_cluster = -1;          ///< argmax over per_cluster (-1 if none)
   double worst_cluster_imbalance = 1.0;
};

/// Compute the imbalance diagnostic above.  Pure reduction over `part` — no
/// METIS, no MPI.  `cell_cost == nullptr` ⇒ uniform per-element weight 1.
/// Ranks/clusters out of range are skipped defensively; a zero-total quantity
/// reports imbalance 1.0.
LtsPartitionImbalance ComputeLtsPartitionImbalance(
   const std::vector<int>& cluster,
   int num_clusters,
   const std::vector<int>& part,
   int nparts,
   const std::vector<double>* cell_cost,
   const std::vector<std::pair<int, int>>& fault_pairs);

/// Block-work balance of a cut (L2 diagnostic): `weights` is a flat
/// `nvtxs*num_clusters` vector (BuildLtsBlockWeights without cell_cost, or
/// any per-cluster work), `part` the cut.
///   per_cluster : max/mean over ranks of each cluster's work (1.0 when a
///                 cluster carries none).
///   tick_model  : the cost of one sync interval of 2^L ticks (L =
///                 num_clusters - 1; tick type k ends clusters 0..k and waits
///                 for the busiest rank), sum_{k<L} 2^(L-k-1) M_k + M_L with
///                 M_k = max_r sum_{c<=k} W_c(r) (the interval of the LTS
///                 report's tick model, solver/nitsche_lts_setup.hpp), over
///                 the same sum with every rank at the mean (1.0 = perfect).
/// Pure reduction; throws std::runtime_error on inconsistent sizes or ranks.
struct LtsBlockBalance
{
   std::vector<double> per_cluster;   ///< length num_clusters
   double tick_model = 1.0;
};

LtsBlockBalance ComputeLtsBlockBalance(
   const std::vector<int>& weights,
   int num_clusters,
   const std::vector<int>& part,
   int nparts);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_PARTITION_HPP
