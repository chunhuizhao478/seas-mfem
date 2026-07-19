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

/// Partition `nvtxs` elements into `nparts` ranks with cluster-balanced weights,
/// then lock every fault pair onto one rank (fault locality).
///
///   xadj / adjncy : CSR element adjacency — `xadj` length `nvtxs+1`, `adjncy`
///                   length `xadj[nvtxs]` (each element's face-neighbours; the
///                   graph must be symmetric and self-loop-free, as METIS
///                   requires).
///   cluster       : per-element cluster id in `[0, num_clusters)`.
///   cell_cost     : per-element work weight (`nullptr` ⇒ uniform 1).
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

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_PARTITION_HPP
