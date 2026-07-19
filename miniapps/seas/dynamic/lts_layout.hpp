// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_layout.hpp — run-side clustered-LTS layout (Appendix A.2 of
//   document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md).
//
// PURE FUNCTION, NO MPI: given the per-(local-)element cluster ids and the
// rank-local face connectivity, produce the per-cluster element lists, the
// per-cluster face lists tagged with roles (the role-driven-sweep tables of
// the Buffer design), the derived coarse-side element sets (`provider_elems`
// retain D(k); `consumer_owner_elems` own accumulate buffers) with dense slot
// maps, and the LOCAL per-cluster element / fault-face counts.  The GLOBAL
// (broadcast) counts that gate collectives are produced by a separate driver-
// side reduction into `LtsGlobalMeta` (see `ReduceGlobalMeta`).
//
// FaceRole semantics (maxdiff <= 1, so a non-fault non-boundary face joins two
// clusters differing by 0 or 1; a fault face joins two elements of the SAME
// cluster, diff 0):
//   IntraClusterGTS   : neighbour in the SAME cluster — visited once, both
//                       sides scatter as in GTS.
//   ConsumerFine      : from the FINE element's side of a diff-1 face — the
//                       fine side does the two-sided work (integrates the coarse
//                       neighbour's Taylor D(k) over its sub-interval) AND fills
//                       the coarse element's accumulate buffer.
//   ProviderCoarseSkip: from the COARSE element's side of a diff-1 face — the
//                       coarse side SKIPS the face in its own sweep (it consumes
//                       the buffer at its correct instead) and RETAINS D(k) for
//                       the finer neighbour.
//   Boundary          : no local neighbour (domain boundary or, in Phase 1, a
//                       rank seam — proper seam handling is Phase 4).
//   Fault             : a fault face (both sides same cluster); tracked in the
//                       cluster's `fault_faces` and not swept as a bulk face.

#ifndef MFEM_SEAS_LTS_LAYOUT_HPP
#define MFEM_SEAS_LTS_LAYOUT_HPP

#include "lts_stepper.hpp"   // LtsGlobalMeta (global metadata target)

#include <cstdint>
#include <cstring>
#include <functional>
#include <vector>

namespace mfem
{
namespace seas
{

/// Role of one (element, face) incidence in a cluster's owned-face sweep.
enum class FaceRole
{
   IntraClusterGTS,
   ConsumerFine,
   ProviderCoarseSkip,
   Boundary,
   Fault
};

/// One rank-local face and the two elements it joins.  `elem2 < 0` marks a face
/// with no local neighbour (domain boundary or, in Phase 1, a rank seam).  Fault
/// faces set `is_fault = true` (both elements share a cluster).
struct LtsFaceSpec
{
   int  face_id = -1;
   int  elem1   = -1;
   int  elem2   = -1;
   bool is_fault = false;
};

/// Per-cluster, per-rank layout (Appendix A.2).  `faces` is the role vector
/// (verbatim A.2); `face_ids` / `face_nbr` are the parallel face id and
/// neighbour-element arrays the Phase-2 role-driven sweep needs.
struct LtsCluster
{
   std::vector<int>      elems;       ///< local element ids in this cluster
   std::vector<int>      face_ids;    ///< local face ids owned by this cluster's sweep
   std::vector<FaceRole> faces;       ///< role per owned face (parallel to face_ids)
   std::vector<int>      face_nbr;    ///< neighbour local elem id per owned face (-1 = none)
   std::vector<int>      fault_faces; ///< local fault face ids in this cluster
};

/// Full rank-local layout (Appendix A.2).
struct LtsLayout
{
   int num_clusters = 0;
   std::vector<LtsCluster> clusters;             ///< size num_clusters

   // Coarse-side element sets (a coarse element that neighbours >=1 finer
   // element is BOTH a D(k) provider and a buffer owner; for rate 2 the two sets
   // coincide, but they carry independent dense slot maps because Phase 2 sizes
   // the D(k) store and the accumulate buffers differently).
   std::vector<int> provider_elems;        ///< sorted, unique coarse-side elems
   std::vector<int> consumer_owner_elems;  ///< sorted, unique coarse-side elems

   // Dense element -> slot maps, sized [n_local_elems], -1 where absent.
   std::vector<int> provider_slot_of_elem; ///< index into provider_elems, else -1
   std::vector<int> buffer_slot_of_elem;   ///< index into consumer_owner_elems, else -1

   // LOCAL metadata (the global broadcast counts are produced by ReduceGlobalMeta).
   std::vector<long long> local_elems_per_cluster;       ///< size num_clusters
   std::vector<long long> local_fault_faces_per_cluster; ///< size num_clusters

   /// Convenience: total owned faces across clusters (for tests / logging).
   long long num_owned_faces() const;
};

/// Build the rank-local layout from cluster ids + the rank-local face list.
///
///   cluster       - per local element cluster id in [0, num_clusters).
///   num_clusters  - Nc (== max(cluster)+1; validated).
///   faces         - every rank-local face exactly once (any order); boundary
///                   faces have elem2 < 0; fault faces set is_fault.
///
/// Deterministic: cluster element lists, provider/consumer sets, and slot maps
/// are all in ascending element-id order; each cluster's face list is in the
/// input face order.  Aborts on: empty cluster vector, out-of-range cluster id,
/// a non-fault diff>1 face (maxdiff violated), or a fault face whose two
/// elements are in different clusters (fault diff!=0 violated).
LtsLayout BuildLtsLayout(const std::vector<int>&          cluster,
                         int                              num_clusters,
                         const std::vector<LtsFaceSpec>&  faces);

/// Reduce the layout's LOCAL per-cluster counts into a broadcast-consistent
/// LtsGlobalMeta (element + fault-face counts per cluster).  `sum_across_ranks`
/// is the caller's MPI_Allreduce(SUM) hook: it receives a length-2*num_clusters
/// buffer (elems[0..Nc), then fault_faces[0..Nc)) and must reduce it in place
/// across all ranks.  Pass nullptr for the serial (single-rank) case.  The other
/// LtsGlobalMeta fields (num_state, drop_fault_predictor_exchange) keep their
/// defaults and are set by the caller if needed.
LtsGlobalMeta ReduceGlobalMeta(
   const LtsLayout&                              layout,
   const std::function<void(std::vector<long long>&)>& sum_across_ranks);

/// (Checkpoint V2) Deterministic 64-bit FNV-1a hash of a clustering, for the
/// restart layout-match check (Phase 2, `io` checkpoint V2 schema).  Hashed, in
/// order, as little-endian bytes: `rate` (int32), `num_clusters` (int32), each
/// `cluster_ids[i]` (int32, in SERIAL-mesh element order — the caller supplies
/// that order), then the raw IEEE-754 bit patterns of `dt_base` and `lambda`.
/// Cross-platform reproducible (integer loop + explicit LE byte order; no float
/// arithmetic on the ids), so a restart with a changed layout is refused rather
/// than silently continued.
inline std::uint64_t LtsLayoutHash(int rate, int num_clusters,
                                   const std::vector<int>& cluster_ids,
                                   double dt_base, double lambda)
{
   std::uint64_t h = 14695981039346656037ULL;         // FNV offset basis
   auto mix_byte = [&](std::uint8_t b)
   { h ^= static_cast<std::uint64_t>(b); h *= 1099511628211ULL; };  // FNV prime
   auto mix_i32 = [&](std::int32_t v)
   {
      const std::uint32_t u = static_cast<std::uint32_t>(v);
      for (int k = 0; k < 4; ++k) { mix_byte(static_cast<std::uint8_t>((u >> (8 * k)) & 0xFF)); }
   };
   auto mix_f64 = [&](double d)
   {
      std::uint64_t u; std::memcpy(&u, &d, sizeof(u));
      for (int k = 0; k < 8; ++k) { mix_byte(static_cast<std::uint8_t>((u >> (8 * k)) & 0xFF)); }
   };
   mix_i32(rate);
   mix_i32(num_clusters);
   for (int c : cluster_ids) { mix_i32(c); }
   mix_f64(dt_base);
   mix_f64(lambda);
   return h;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_LAYOUT_HPP
