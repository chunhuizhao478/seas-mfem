// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/godunov_flux_pool.cpp — Phase H.1 implementation.

#include "godunov_flux_pool.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>

namespace mfem
{
namespace seas
{

namespace
{

/// Round a real_t to `sig` significant figures, preserving sign and
/// near-zero values.  Used as the dedup key for (λ, μ, ρ).
real_t round_sig(real_t v, int sig)
{
   if (v == 0.0 || !std::isfinite(v)) { return v; }
   const real_t mag = std::abs(v);
   const real_t log_mag = std::log10(mag);
   // R-312 / R-406: range-guard the input.  For |log_mag| > ~290 the
   // exponent in std::pow can overflow to +inf, producing NaN in the
   // rounded result and corrupting the dedup key.  Silently returning
   // the unrounded value (round-3 R-312) would also break dedup
   // determinism — two inputs that should land in the same bucket
   // would land in different ones.  Material values in SAFS are
   // O(1) – O(1e10), far inside the safe window; anything outside is
   // almost certainly corrupted input and we MFEM_VERIFY-abort.
   MFEM_VERIFY(std::isfinite(log_mag) && std::abs(log_mag) <= 290.0,
               "GodunovFluxPool::round_sig: |log10(|v|)| = "
               << std::abs(log_mag)
               << " is outside the safe rounding range; got v = " << v
               << ".  Expected O(1)..O(1e10) for material values.");
   const real_t scale = std::pow(static_cast<real_t>(10.0),
                                 sig - 1 - std::floor(log_mag));
   const real_t r = std::round(v * scale) / scale;
   return r;
}

/// Encode (lambda, mu, rho) rounded to `sig` figures into a string
/// key — using a string avoids hash-collision pitfalls with a custom
/// composite real_t key.
std::string make_key(const std::array<real_t, 3>& lmr, int sig)
{
   std::ostringstream oss;
   oss.precision(sig + 2);
   oss << round_sig(lmr[0], sig) << "|"
       << round_sig(lmr[1], sig) << "|"
       << round_sig(lmr[2], sig);
   return oss.str();
}

}  // namespace

void GodunovFluxPool::Build(
   int ne,
   const std::vector<std::array<real_t, 3>>& per_elem_lmr,
   int dedup_sig_figs)
{
   MFEM_VERIFY(ne >= 0,
               "GodunovFluxPool::Build: ne must be >= 0; got " << ne);
   MFEM_VERIFY(static_cast<int>(per_elem_lmr.size()) == ne,
               "GodunovFluxPool::Build: per_elem_lmr.size() ("
               << per_elem_lmr.size() << ") must equal ne (" << ne << ")");
   MFEM_VERIFY(dedup_sig_figs >= 1 && dedup_sig_figs <= 15,
               "GodunovFluxPool::Build: dedup_sig_figs ("
               << dedup_sig_figs
               << ") must be in [1, 15]");

   elem_to_flux_idx_.assign(ne, -1);
   unique_fluxes_.clear();

   std::unordered_map<std::string, int> key_to_idx;
   key_to_idx.reserve(static_cast<size_t>(ne) / 4 + 16);

   for (int e = 0; e < ne; ++e)
   {
      const auto& lmr = per_elem_lmr[e];
      const real_t lam = lmr[0];
      const real_t mu  = lmr[1];
      const real_t rho = lmr[2];
      MFEM_VERIFY(rho > 0.0,
                  "GodunovFluxPool::Build: rho (" << rho
                  << ") must be > 0 at element " << e);
      MFEM_VERIFY(lam + 2.0 * mu > 0.0,
                  "GodunovFluxPool::Build: (lambda + 2*mu) ("
                  << (lam + 2.0 * mu) << ") must be > 0 at element " << e);

      const std::string key = make_key(lmr, dedup_sig_figs);
      auto it = key_to_idx.find(key);
      if (it == key_to_idx.end())
      {
         const int new_idx = static_cast<int>(unique_fluxes_.size());
         // Use the rounded values for the cached flux so per-element
         // queries are exactly consistent with the dedup decision.
         const real_t rl = round_sig(lam, dedup_sig_figs);
         const real_t rm = round_sig(mu,  dedup_sig_figs);
         const real_t rr = round_sig(rho, dedup_sig_figs);
         unique_fluxes_.emplace_back(std::make_unique<GodunovFlux>(rl, rm, rr));
         key_to_idx.emplace(key, new_idx);
         elem_to_flux_idx_[e] = new_idx;
      }
      else
      {
         elem_to_flux_idx_[e] = it->second;
      }
   }
}

}  // namespace seas
}  // namespace mfem
