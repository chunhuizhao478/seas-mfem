// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_forced_rupture.hpp — option O4 of the complex-geometries plan
// (rev 2, Sections 5.5 and 10.3; Chunhui Zhao 2026-10-20): forced-rupture
// nucleation, for TPV22/23 (the same form as TPV26/27).
//
// Deck (Nitsche driver only; absent = every bit as before):
//
//   [forced_rupture]
//   hypocenter_m      = [-10000.0, 0.0, -10000.0]   # global X, Y, Z
//   r_crit_m          = 3000.0
//   rupture_speed_m_s = 2424.8       # e.g. 0.7 Vs
//   smoothing         = 0.081
//   t0_s              = 0.5
//
// Each fault DOF gets a forcing time from its 3-D distance r to the
// hypocentre,
//    T(r) = r / v + smoothing r_c / v (1 / (1 - (r / r_c)^2) - 1),  r < r_c,
//    T    = +inf (never forced),                                    r >= r_c,
// and f2(t) = clamp((t - T) / t0, 0, 1).  The friction coefficient is
//    mu = mu_s + (mu_d - mu_s) max(f1, f2),   f1 = min(D / Dc, 1),
// which NitscheWaveOperator::SetForcedRupture obtains without touching the
// closure: since f2 <= 1, max(min(D/Dc, 1), f2) = min(max(D/Dc, f2), 1),
// so the law is handed D_eff = max(D, f2 Dc) (D = the slip path with O3).
// f2 uses the time the fault advance receives (the fault cluster's time
// under local time stepping).  Slip weakening only.

#ifndef MFEM_SEAS_NITSCHE_OPT_FORCED_RUPTURE_HPP
#define MFEM_SEAS_NITSCHE_OPT_FORCED_RUPTURE_HPP

#include "mfem.hpp"
#include "nitsche_opt_common.hpp"

#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

struct ForcedRuptureConfig
{
   static constexpr const char *kTable = "forced_rupture";
   real_t hyp[3] = {0.0, 0.0, 0.0};
   real_t r_crit = 0.0, v_r = 0.0, smoothing = 0.0, t0 = 0.0;

   /// Forcing time at distance r from the hypocentre (+inf for r >= r_crit).
   real_t Time(real_t r) const
   {
      if (!(r < r_crit)) { return std::numeric_limits<real_t>::infinity(); }
      const real_t q = r / r_crit;
      return r / v_r + smoothing * r_crit / v_r * (1.0 / (1.0 - q * q) - 1.0);
   }
   /// Forcing time at position x.
   real_t TimeAt(const real_t x[3]) const
   {
      return Time(std::sqrt((x[0] - hyp[0]) * (x[0] - hyp[0]) + (x[1] - hyp[1]) * (x[1] - hyp[1])
                            + (x[2] - hyp[2]) * (x[2] - hyp[2])));
   }
};

inline ForcedRuptureConfig ParseForcedRupture(const toml::value &t)
{
   using T = ForcedRuptureConfig;
   OptRefuseUnknownKeys(t, T::kTable, {"hypocenter_m", "r_crit_m", "rupture_speed_m_s",
                                       "smoothing", "t0_s"});
   T c;
   const std::vector<real_t> h = OptRealArray(t, T::kTable, "hypocenter_m");
   MFEM_VERIFY(h.size() == 3, "[" << T::kTable << "] hypocenter_m needs 3 values (X, Y, Z)");
   for (int i = 0; i < 3; i++) { c.hyp[i] = h[i]; }
   c.r_crit = OptReal(t, T::kTable, "r_crit_m");
   c.v_r = OptReal(t, T::kTable, "rupture_speed_m_s");
   c.smoothing = OptReal(t, T::kTable, "smoothing");
   c.t0 = OptReal(t, T::kTable, "t0_s");
   MFEM_VERIFY(c.r_crit > 0.0 && c.v_r > 0.0 && c.t0 > 0.0 && c.smoothing >= 0.0,
               "[" << T::kTable << "]: r_crit_m, rupture_speed_m_s and t0_s must be > 0 and "
               "smoothing >= 0 (got " << c.r_crit << ", " << c.v_r << ", " << c.t0 << ", "
               << c.smoothing << ")");
   return c;
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_FORCED_RUPTURE_HPP
