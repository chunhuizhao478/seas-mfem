// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pulse_load.hpp — borehole pressure of the pulse-power experiments.
//
// Exact port of ElkPulseLoadExperiment::value (farms branch cdms,
// src/functions/pulsepower/ElkPulseLoadExperiment.C), used by both models:
//
//   D(x)  distance in mm from x to the electrode-gap segment on the axis
//         through discharge_center, the gap of length `gap` centred at z_c
//         (inside the gap: radial distance only)
//   Pp    = 0.1e6 * B * (eta EM)^kappa / D^m      [Pa]   (EM in kJ)
//         or the fixed `peak_pressure` when it is > 0
//   p(t)  = Pp (e^{-a tau} - e^{-b tau}) / (e^{-a tr} - e^{-b tr}),
//         tau = fmod(t, Tw), for t <= Np Tw, else 0;
//         optionally max(p, p_min)
//
// The sign convention is MOOSE's Pressure BC: the traction is -p n.

#ifndef MFEM_PFF_PULSE_LOAD_HPP
#define MFEM_PFF_PULSE_LOAD_HPP

#include "mfem.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace pff
{

struct PulseLoadParams
{
   real_t shape_alpha = 0.0;        ///< a [1/s]
   real_t shape_beta = 0.0;         ///< b [1/s]
   real_t rise_time = 0.0;          ///< t_r [s]
   real_t pulse_duration = 0.0;     ///< window T_w [s]
   real_t EM = 0.0;                 ///< stored energy per pulse [kJ]
   real_t gap = 0.0;                ///< electrode gap [m]
   real_t efficiency = 1.0;         ///< eta
   real_t kappa = 0.35;             ///< fitting_param_alpha
   real_t m_exp = 1.0;              ///< fitting_param_exponent (MOOSE default 1)
   real_t base_factor = 9000.0;     ///< B (MOOSE default 9000; inputs use 8000)
   real_t center[3] = {0.0, 0.0, 0.0};
   int num_pulses = 1;
   real_t peak_pressure = 0.0;      ///< > 0 overrides the distance law
   bool use_min_pressure = false;
   real_t min_pressure = 0.0;
};

class PulseLoad
{
public:
   explicit PulseLoad(const PulseLoadParams &p) : P_(p) {}

   const PulseLoadParams &Params() const { return P_; }

   /// Distance to the gap segment in mm (ElkPulseLoadExperiment's Dsensor).
   real_t DistanceMM(const real_t *x) const
   {
      const real_t dx = x[0] - P_.center[0], dy = x[1] - P_.center[1];
      const real_t z = x[2];
      const real_t zlo = P_.center[2] - P_.gap / 2.0, zhi = P_.center[2] + P_.gap / 2.0;
      real_t dz = 0.0;
      if (z < zlo) { dz = z - zlo; }
      else if (z > zhi) { dz = z - zhi; }
      return 1000.0 * std::sqrt(dx*dx + dy*dy + dz*dz);
   }

   /// Peak pressure at x [Pa].
   real_t Peak(const real_t *x) const
   {
      if (P_.peak_pressure > 0.0) { return P_.peak_pressure; }
      const real_t D = DistanceMM(x);
      return (0.1 * 1e6) * P_.base_factor / std::pow(D, P_.m_exp) *
             std::pow(P_.efficiency * P_.EM, P_.kappa);
   }

   /// Normalised time history in [0, 1] (peak 1 at tau = t_r).
   real_t Shape(real_t t) const
   {
      if (t > P_.num_pulses * P_.pulse_duration) { return 0.0; }
      const real_t tau = std::fmod(t, P_.pulse_duration);
      const real_t a = P_.shape_alpha, b = P_.shape_beta, tr = P_.rise_time;
      return (std::exp(-a * tau) - std::exp(-b * tau)) /
             (std::exp(-a * tr) - std::exp(-b * tr));
   }

   /// p(t, x) [Pa]; x has 3 components (pad z = 0 in 2D, as MOOSE does).
   real_t Value(real_t t, const real_t *x) const
   {
      real_t p = Peak(x) * Shape(t);
      if (P_.use_min_pressure) { p = std::max(p, P_.min_pressure); }
      return p;
   }

private:
   PulseLoadParams P_;
};

/// The rates (a, b) of a double exponential with its peak at t_r and 10 % of the
/// peak at t_d (spe_hftc paper eq. 14), by Newton on b with
/// a = a(b) from ln(b/a)/(b-a) = t_r.  Used to derive inputs, not at run time.
inline bool ShapeRatesFromRiseAndDuration(real_t tr, real_t td, real_t &a, real_t &b)
{
   // a(b): solve ln(b/a) = tr (b - a) for a < b by bisection.
   auto a_of_b = [tr](real_t bb)
   {
      real_t lo = 1e-12 * bb, hi = bb * (1.0 - 1e-12);
      for (int it = 0; it < 200; it++)
      {
         const real_t mid = 0.5 * (lo + hi);
         const real_t f = std::log(bb / mid) - tr * (bb - mid);
         if (f > 0.0) { lo = mid; } else { hi = mid; }
      }
      return 0.5 * (lo + hi);
   };
   auto ratio = [&](real_t bb)
   {
      const real_t aa = a_of_b(bb);
      return (std::exp(-aa * td) - std::exp(-bb * td)) /
             (std::exp(-aa * tr) - std::exp(-bb * tr)) - 0.1;
   };
   // bracket b in (1/tr, 1e4/tr): larger b -> smaller a -> longer tail
   real_t lo = 1.0001 / tr, hi = 1e4 / tr;
   real_t flo = ratio(lo), fhi = ratio(hi);
   if (flo * fhi > 0.0) { return false; }
   for (int it = 0; it < 200; it++)
   {
      const real_t mid = std::sqrt(lo * hi);
      const real_t fm = ratio(mid);
      if ((fm > 0.0) == (flo > 0.0)) { lo = mid; flo = fm; } else { hi = mid; }
   }
   b = std::sqrt(lo * hi);
   a = a_of_b(b);
   return true;
}

} // namespace pff
} // namespace mfem

#endif
