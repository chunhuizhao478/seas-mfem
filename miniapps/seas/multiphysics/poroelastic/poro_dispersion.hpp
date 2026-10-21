// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_dispersion.hpp — closed-form Biot wave speeds (poroelastic plan rev 3,
// step P1.2; poro-derivation-rev1.2 §4.5, eq. (4.6)).
//
// Plane waves without drag (the high-frequency limit, which sets the time
// step).  The P speeds c solve
//
//    det( [[lambda_u + 2 mu, alpha M], [alpha M, M]] - c^2 [[rho, rho_f], [rho_f, m]] ) = 0,
//
// a quadratic in s = c^2:  a s^2 - q s + r = 0 with
//    a = rho m - rho_f^2,   q = (lambda_u + 2 mu) m + M rho - 2 alpha M rho_f,
//    r = (lambda_u + 2 mu) M - (alpha M)^2 = (lambda + 2 mu) M.
// The larger root (q + sqrt(q^2 - 4 a r)) / (2 a) is the fast P speed; the
// smaller is taken as r / (a s_fast), the product of the roots, which avoids
// the cancellation of the minus sign.  S: c_s^2 = mu / (rho - rho_f^2 / m) at
// high frequency, mu / rho at low frequency.  The undrained and drained P
// speeds are sqrt((lambda_u + 2 mu) / rho) and sqrt((lambda + 2 mu) / rho).

#ifndef MFEM_SEAS_PORO_DISPERSION_HPP
#define MFEM_SEAS_PORO_DISPERSION_HPP

#include "poro_material.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{
namespace poro
{

struct BiotSpeeds
{
   real_t c_fast = 0.0;          ///< fast P, high frequency (sets the time step)
   real_t c_slow = 0.0;          ///< slow P, high frequency
   real_t c_s_hf = 0.0;          ///< S, high frequency
   real_t c_s_lf = 0.0;          ///< S, low frequency
   real_t c_p_undrained = 0.0;   ///< P with lambda_u (seismic frequencies)
   real_t c_p_drained = 0.0;     ///< P with lambda
};

/// Eq. (4.6) and the limits above, from one element's constants.
inline BiotSpeeds ComputeBiotSpeeds(const BiotConstants &c)
{
   const real_t K11 = c.lambda_u + 2.0 * c.mu;
   const real_t K12 = c.alpha * c.M;
   const real_t a = c.rho * c.m - c.rho_f * c.rho_f;
   const real_t q = K11 * c.m + c.M * c.rho - 2.0 * K12 * c.rho_f;
   const real_t r = (c.lambda + 2.0 * c.mu) * c.M;
   const real_t s_fast = (q + std::sqrt(q * q - 4.0 * a * r)) / (2.0 * a);
   const real_t s_slow = r / (a * s_fast);
   BiotSpeeds s;
   s.c_fast = std::sqrt(s_fast);
   s.c_slow = std::sqrt(s_slow);
   s.c_s_hf = std::sqrt(c.mu / (c.rho - c.rho_f * c.rho_f / c.m));
   s.c_s_lf = std::sqrt(c.mu / c.rho);
   s.c_p_undrained = std::sqrt(K11 / c.rho);
   s.c_p_drained = std::sqrt((c.lambda + 2.0 * c.mu) / c.rho);
   return s;
}

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_DISPERSION_HPP
