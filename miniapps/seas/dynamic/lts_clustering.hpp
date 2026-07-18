// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_clustering.hpp — Phase 0 of the clustered-LTS plan
//   (document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md, Appendix A.1).
//
// PURE FUNCTIONS, NO MPI.  Bins elements into rate-2 clusters by their allowed
// CFL time step, enforces the SeisSol neighbour-difference rule (adjacent
// clusters differ by <= 1; the two elements of a fault face differ by 0),
// optionally rescales the bin edges by a wiggle factor lambda and merges the
// top clusters under a cost model, and reports the modelled cost.
//
// Determinism (a hard requirement — the run-side layout that consumes cluster
// ids must be a pure function of the serial mesh + material + config, or the
// cross-rank fault tripwire aborts): binning uses an INTEGER comparison loop,
// never floor(log2(.)), so cluster ids are bit-reproducible across platforms.
//
// This header has ZERO dependence on the wave operator, the driver, or MPI;
// it takes the per-element allowed dt, the element face-adjacency graph, and
// the fault-face element pairs, and returns the cluster assignment.

#ifndef MFEM_SEAS_LTS_CLUSTERING_HPP
#define MFEM_SEAS_LTS_CLUSTERING_HPP

#include "mfem.hpp"   // mfem::Table

#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// Tunables for BuildLtsClustering (Appendix A.1).
struct LtsClusteringOptions
{
   int    rate         = 2;      ///< cluster c steps at dt_base * rate^c (only 2 supported in v1)
   int    max_clusters = 32;     ///< hard ceiling on raw bins (pre-merge)
   bool   wiggle_scan  = true;   ///< scan lambda in (0.5, 1] step 0.01; false => lambda_fixed
   double lambda_fixed = 1.0;    ///< used when wiggle_scan == false
   int    nc_cap       = 6;      ///< auto-merge target; <= 0 disables merge (RAW mode)
   double merge_loss_tol = 0.05; ///< accept a merge while cost <= (1+tol)*cost(uncapped)
   /// Per-element partition/update cost weight (defaults to uniform 1.0 when null).
   /// Must have the same length as dt_e when non-null.
   const std::vector<double>* cell_cost = nullptr;
};

/// Result of BuildLtsClustering (Appendix A.1).
struct LtsClustering
{
   std::vector<int> cluster;   ///< per element; post-clamp, post-merge; 0-based contiguous
   int    num_clusters = 0;    ///< Nc after cap/auto-merge
   double dt_base      = 0.0;  ///< = lambda * min(dt_e); cluster c steps at dt_base * 2^c
   double lambda       = 1.0;  ///< chosen wiggle factor (1.0 in RAW mode)
   double modeled_cost = 0.0;  ///< sum_e cellCost_e / (2^{c_e} * dt_base), post-clamp/merge

   /// Histogram: cells[c] = number of elements in cluster c (size num_clusters).
   std::vector<long long> cells_per_cluster() const;
};

/// Build the rate-2 cluster assignment from per-element allowed dt.
///
/// Inputs:
///   dt_e                 - allowed CFL time step per element (dt_e[i] > 0).
///   elem_to_elem         - element face-adjacency graph (row i lists element
///                          i's face neighbours; symmetric — as produced by
///                          mfem::Mesh::ElementToElementTable()).  Boundary
///                          faces / rank seams simply have fewer neighbours.
///   fault_face_elem_pairs- the (elem1, elem2) pairs joined by a fault face;
///                          these are forced into the SAME cluster (diff == 0).
///   opt                  - tunables (Appendix A.1).
///
/// Aborts (MFEM_ABORT) on: empty dt_e, non-positive dt_e[i], a fault pair with
/// an out-of-range element id, or cell_cost length mismatch.
///
/// Guarantees (asserted internally, GAP-A1): for every element i,
///   dt_base * 2^{cluster[i]} <= dt_e[i]   (each element steps no faster than
/// its own CFL limit), using the identical integer-edge comparison as binning.
LtsClustering BuildLtsClustering(
   const std::vector<double>&               dt_e,
   const mfem::Table&                       elem_to_elem,
   const std::vector<std::pair<int, int>>&  fault_face_elem_pairs,
   const LtsClusteringOptions&              opt = LtsClusteringOptions());

/// SeisSol's "theoretical speedup" statistic: the ARITHMETIC mean of
/// dt_cluster/dt_min over all cells = (1/N) * sum_c cells_c * 2^c.  This is the
/// number SeisSol prints ("Theoretical speedup ... clustered LTS"); it is NOT
/// the workload ratio (it overweights coarse cells).  Provided only for the
/// Phase-0 cross-check against SeisSol's log.
double ArithmeticMeanSpeedup(const LtsClustering& cl);

/// The HONEST workload ratio: element-updates under GTS-at-dt_base divided by
/// element-updates under the clustering = N / sum_c cells_c / 2^c.  This is the
/// go/no-go number.  (Uses uniform per-cell weight; cost-weighted variants use
/// modeled_cost directly.)
double HarmonicUpdateSpeedup(const LtsClustering& cl);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_CLUSTERING_HPP
