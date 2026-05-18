// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/godunov_flux_pool.hpp — Phase H.1 of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// `GodunovFluxPool`: a per-element cache of `GodunovFlux` instances,
// deduplicated by the rounded `(λ, μ, ρ)` triple at the element
// centroid.  The heterogeneous-material `WaveOperator(MaterialField,
// BoundaryConfig)` constructor (Phase H.2) builds one of these from
// per-element material evaluations so the ADER hot loop can call
// `flux_pool_->At(elem).Interior(...)` instead of using a single
// scalar-material `flux_` member.
//
// Existing `GodunovFlux` (one per WaveOperator, scalar material) is
// preserved verbatim — TPV/BP5 runs that don't set a `flux_pool_`
// continue to use it.  This file adds an ADDITIVE class only; nothing
// else in `dynamic/godunov_flux.hpp` is touched.

#ifndef MFEM_SEAS_DYNAMIC_GODUNOV_FLUX_POOL_HPP
#define MFEM_SEAS_DYNAMIC_GODUNOV_FLUX_POOL_HPP

#include "mfem.hpp"

#include "godunov_flux.hpp"

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

namespace mfem
{
namespace seas
{

/// Per-element flux cache deduplicated by rounded (lambda, mu, rho).
///
/// Usage (Phase H.2 driver / wave_operator):
///   GodunovFluxPool pool;
///   std::vector<std::array<real_t, 3>> per_elem_lmr(ne);
///   for (int e = 0; e < ne; ++e) { ... evaluate (lam, mu, rho) per element ...
///                                  per_elem_lmr[e] = {lam, mu, rho}; }
///   pool.Build(ne, per_elem_lmr, /*dedup_sig_figs=*/6);
///
///   ... in the hot loop ...
///   const GodunovFlux& flux_e = pool.At(elem);
///   flux_e.Interior(nor, Q_self, Q_nbr, F_h);
///
/// On a typical SAFS 1000 m mesh (~4M tets, smooth basin material),
/// `dedup_sig_figs = 6` produces ≤ ~10000 unique triples → ≤ ~10000
/// `GodunovFlux` instances → memory budget << 2 GB (plan §Risk
/// Assessment).
class GodunovFluxPool
{
public:
   GodunovFluxPool() = default;

   /// Build the pool from per-element (lambda, mu, rho).
   /// Deduplicates triples to `dedup_sig_figs` significant figures
   /// (default 6).  Aborts if any element has rho <= 0 or
   /// (lambda + 2*mu) <= 0.
   ///
   /// @param ne                Total number of elements (size of per_elem_lmr).
   /// @param per_elem_lmr      [ne] (lambda, mu, rho) at the element centroid.
   /// @param dedup_sig_figs    Significant figures used to round before
   ///                          deduplication.  Default 6.
   void Build(int ne,
              const std::vector<std::array<real_t, 3>>& per_elem_lmr,
              int dedup_sig_figs = 6);

   /// Read the cached GodunovFlux for element `e`.
   /// Aborts if e is out of range or the pool is not yet built.
   const GodunovFlux& At(int e) const
   {
      MFEM_ASSERT(!unique_fluxes_.empty(),
                  "GodunovFluxPool::At: pool not built (call Build() first)");
      MFEM_ASSERT(e >= 0 && e < static_cast<int>(elem_to_flux_idx_.size()),
                  "GodunovFluxPool::At(" << e << "): element index out of "
                  "range [0, " << elem_to_flux_idx_.size() << ")");
      return *unique_fluxes_[elem_to_flux_idx_[e]];
   }

   int NumElements() const
   { return static_cast<int>(elem_to_flux_idx_.size()); }

   /// Number of unique (rounded) (λ, μ, ρ) triples — the size of the
   /// underlying flux cache.  Test-only accessor.
   int NumUniqueTriples() const
   { return static_cast<int>(unique_fluxes_.size()); }

   /// Test-only: flux index for element `e`.
   int FluxIndexFor(int e) const
   {
      MFEM_ASSERT(e >= 0 && e < static_cast<int>(elem_to_flux_idx_.size()),
                  "GodunovFluxPool::FluxIndexFor(" << e << "): out of range");
      return elem_to_flux_idx_[e];
   }

private:
   std::vector<int>                          elem_to_flux_idx_;
   std::vector<std::unique_ptr<GodunovFlux>> unique_fluxes_;
};

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_DYNAMIC_GODUNOV_FLUX_POOL_HPP
