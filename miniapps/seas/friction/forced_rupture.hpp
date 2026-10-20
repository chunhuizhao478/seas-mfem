// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// forced_rupture.hpp — smoothed forced rupture for linear slip weakening
// (SCEC TPV26/27 description v13, Parts 4 and 5; rheology plan rev 1,
// section 8).
//
// Each fault point gets a time of forced rupture T(r), r its distance to the
// hypocenter:
//
//    T(r) = r / (0.7 Vs) + 0.081 r_c / (0.7 Vs) (1 / (1 - (r/r_c)^2) - 1),  r < r_c
//    T(r) = 1.0e9,                                                          r >= r_c
//
// and from T on the friction coefficient is forced down over t0:
//
//    mu = mu_s + (mu_d - mu_s) max(f1, f2),   f1 = min(D/Dc, 1),
//    f2 = 0 (t < T),  (t - T)/t0 (T <= t < T + t0),  1 (t >= T + t0).
//
// The forced front starts at 0.7 Vs at the hypocenter and stops at r_c.  r
// is the straight-line distance in 3-D, which on a planar fault is the
// in-plane distance of the description.  The arithmetic of
// ForcedRupture::FrictionCoefficient is that of the SAFS helper
// spatial::LSWFrictionCoefficient_ForcedRupture (spatial_friction.hpp), in
// the same order (test U11), without its strength-barrier sentinel: the
// Nitsche fault is only the fault faces, so it needs none.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_FORCED_RUPTURE_HPP
#define MFEM_SEAS_FORCED_RUPTURE_HPP

#include "mfem.hpp"

#include <cmath>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

struct ForcedRuptureSpec
{
   real_t hypocenter[3] = {0.0, 0.0, 0.0};   ///< global frame [m]
   real_t r_crit = 0.0;   ///< radius of the forced zone r_c [m]
   real_t vs = 0.0;       ///< shear-wave speed Vs of the formula [m/s]
   real_t t0 = 0.0;       ///< forced-rupture decay time [s]
};

class ForcedRupture
{
public:
   /// T of a point that is never forced (the description's value).
   static constexpr real_t kNeverForced = 1.0e9;
   /// The constants of the description's T(r): the front speed at the
   /// hypocenter as a fraction of Vs, and the smoothing near r_c.
   static constexpr real_t kSpeedFraction = 0.7;
   static constexpr real_t kSmoothing = 0.081;

   /// "" when the spec is usable: r_crit, vs, t0 > 0 and finite.
   static std::string Validate(const ForcedRuptureSpec &s)
   {
      std::ostringstream e;
      for (int i = 0; i < 3; i++)
      {
         if (!std::isfinite(s.hypocenter[i])) { e << "hypocenter must be finite"; return e.str(); }
      }
      if (!(s.r_crit > 0.0) || !std::isfinite(s.r_crit)) { e << "r_crit_m must be > 0"; }
      else if (!(s.vs > 0.0) || !std::isfinite(s.vs)) { e << "vs_m_s must be > 0"; }
      else if (!(s.t0 > 0.0) || !std::isfinite(s.t0)) { e << "t0_s must be > 0"; }
      return e.str();
   }

   /// T(r) at the point x (global frame).
   static real_t Time(const ForcedRuptureSpec &s, const real_t x[3])
   {
      real_t r2 = 0.0;
      for (int i = 0; i < 3; i++)
      {
         const real_t d = x[i] - s.hypocenter[i];
         r2 += d * d;
      }
      return TimeAtDistance(s, std::sqrt(r2));
   }

   static real_t TimeAtDistance(const ForcedRuptureSpec &s, real_t r)
   {
      if (r >= s.r_crit) { return kNeverForced; }
      const real_t v = kSpeedFraction * s.vs;
      const real_t q = r / s.r_crit;
      return r / v + kSmoothing * s.r_crit / v * (1.0 / (1.0 - q * q) - 1.0);
   }

   /// mu with forced rupture at slip measure D and time t.
   static real_t FrictionCoefficient(real_t D, real_t mu_s, real_t mu_d, real_t Dc,
                                     real_t t, real_t T, real_t t0)
   {
      real_t f1;
      if (D <= 0.0) { f1 = 0.0; }
      else if (D >= Dc) { f1 = 1.0; }
      else { f1 = D / Dc; }
      real_t f2;
      if (t < T) { f2 = 0.0; }
      else if (t < T + t0) { f2 = (t - T) / t0; }
      else { f2 = 1.0; }
      const real_t f = (f1 > f2) ? f1 : f2;
      return mu_s + (mu_d - mu_s) * f;
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FORCED_RUPTURE_HPP
