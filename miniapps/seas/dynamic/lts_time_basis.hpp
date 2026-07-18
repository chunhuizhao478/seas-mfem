// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_time_basis.hpp — Taylor time-integration for clustered LTS (Appendix A.4
//   of PLAN_clustered_lts_ader_2026-07-18.md).
//
// A provider (coarse) element's ADER predictor produces a RAW/unscaled Taylor
// derivative stack D(0..order-1): the element's state at a time tau after the
// provider's expansion point is
//
//     Q(tau) = sum_{k=0}^{order-1}  D(k) * tau^k / k!
//
// A finer neighbour that consumes this state over a sub-interval [a, b]
// (a, b measured from the SAME expansion point) needs the TIME INTEGRAL
//
//     integral_a^b Q(tau) dtau
//       = sum_k D(k) / k! * (b^{k+1} - a^{k+1}) / (k+1)
//       = sum_k D(k) * (b^{k+1} - a^{k+1}) / (k+1)!
//       = sum_k coeff[k] * D(k),   coeff[k] = (b^{k+1} - a^{k+1}) / (k+1)!
//
// This is the ONLY place the raw stacks are folded against a sub-interval; the
// coefficients depend only on (a, b, order), never on the field data, so this
// function is pure and fully unit-testable (Appendix B.3).

#ifndef MFEM_SEAS_LTS_TIME_BASIS_HPP
#define MFEM_SEAS_LTS_TIME_BASIS_HPP

#include "mfem.hpp"   // mfem::real_t, MFEM_ASSERT

namespace mfem
{
namespace seas
{

/// Compute the Taylor sub-interval coefficients coeff[k] = (b^{k+1}-a^{k+1})/(k+1)!
/// for k = 0..order-1 into `coeff` (length >= order).  Exposed for tests /
/// callers that want the pure coefficients without touching field data.
inline void TaylorIntegralCoeffs(mfem::real_t a, mfem::real_t b, int order,
                                 mfem::real_t* coeff)
{
   MFEM_ASSERT(order >= 1, "TaylorIntegralCoeffs: order must be >= 1.");
   // Build a^{k+1}, b^{k+1} incrementally and divide by (k+1)! incrementally.
   // factorial_{k+1} = (k+1)!  ->  running product.
   mfem::real_t apow = a;                 // a^{k+1}, starts at a^1
   mfem::real_t bpow = b;                 // b^{k+1}, starts at b^1
   mfem::real_t fact = mfem::real_t(1);   // (k+1)! , starts at 1! = 1
   for (int k = 0; k < order; ++k)
   {
      // here apow = a^{k+1}, bpow = b^{k+1}, fact = (k+1)!
      coeff[k] = (bpow - apow) / fact;
      apow *= a;
      bpow *= b;
      fact *= static_cast<mfem::real_t>(k + 2);   // -> (k+2)! for the next step
   }
}

/// Fold a RAW Taylor derivative stack `dk` against the sub-interval [a, b].
///
///   dk   : length order*n, laid out [k][j] = dk[k*n + j]  (k the Taylor index,
///          j the flattened (component, node) index; in production n = NUM_STATE
///          * ndof_per_block).
///   out  : length n; OVERWRITTEN with the time integral over [a, b].
///
/// a, b are relative to the provider's expansion point (its last predict time);
/// b >= a is expected but not required (b < a just yields the negated integral).
inline void IntegrateTaylor(mfem::real_t a, mfem::real_t b,
                            const mfem::real_t* dk,
                            int order, int n,
                            mfem::real_t* out)
{
   MFEM_ASSERT(order >= 1, "IntegrateTaylor: order must be >= 1.");
   MFEM_ASSERT(n >= 0, "IntegrateTaylor: n must be >= 0.");

   // k = 0 term (overwrite), then accumulate k = 1..order-1.
   mfem::real_t coeff[64];
   MFEM_ASSERT(order <= 64, "IntegrateTaylor: order exceeds coeff buffer (64).");
   TaylorIntegralCoeffs(a, b, order, coeff);

   const mfem::real_t c0 = coeff[0];
   for (int j = 0; j < n; ++j) { out[j] = c0 * dk[j]; }
   for (int k = 1; k < order; ++k)
   {
      const mfem::real_t ck = coeff[k];
      const mfem::real_t* dkk = dk + static_cast<std::size_t>(k) * n;
      for (int j = 0; j < n; ++j) { out[j] += ck * dkk[j]; }
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_TIME_BASIS_HPP
