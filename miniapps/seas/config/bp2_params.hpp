// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_BP2_PARAMS_HPP
#define MFEM_SEAS_BP2_PARAMS_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// Parameters for the SCEC SEAS BP2-QD benchmark problem.
///
/// BP2 is a 2D antiplane shear problem with rate-and-state friction.
/// Reference: https://strike.scec.org/cvws/seas/download/SEAS_BP2_QD.pdf
///
/// All values are in SI units:
/// - Length: meters [m]
/// - Time: seconds [s]
/// - Stress: Pascals [Pa]
/// - Velocity: meters/second [m/s]
struct BP2Params
{
   // =========================================================================
   // Material properties
   // =========================================================================

   /// Rock density [kg/m^3]
   real_t rho = 2670.0;

   /// Shear wave speed [m/s]
   real_t cs = 3464.0;

   /// Shear modulus: mu = rho * cs^2 [Pa]
   real_t mu() const { return rho * cs * cs; }

   // =========================================================================
   // Rate-and-state friction parameters
   // =========================================================================

   /// Reference slip rate [m/s]
   real_t V0 = 1.0e-6;

   /// Reference friction coefficient [-]
   real_t f0 = 0.6;

   /// State evolution effect parameter b [-]
   real_t b = 0.015;

   /// Critical slip distance [m] (key BP2 parameter, differs from BP1)
   real_t Dc = 0.004;

   /// Direct effect parameter in velocity-weakening zone (z < H) [-]
   real_t a0 = 0.010;

   /// Direct effect parameter in velocity-strengthening zone (z >= H+h) [-]
   real_t amax = 0.025;

   // =========================================================================
   // Stress parameters
   // =========================================================================

   /// Normal stress [Pa] (50 MPa)
   real_t sigma_n = 50.0e6;

   /// Radiation damping coefficient: eta = mu / (2 * cs) [Pa*s/m]
   real_t eta() const { return mu() / (2.0 * cs); }

   // =========================================================================
   // Loading parameters
   // =========================================================================

   /// Plate rate (far-field loading velocity) [m/s]
   real_t Vp = 1.0e-9;

   /// Initial slip rate [m/s] (same as plate rate for this benchmark)
   real_t V_init = 1.0e-9;

   // =========================================================================
   // Geometric parameters
   // =========================================================================

   /// Depth to bottom of velocity-weakening zone [m] (15 km)
   real_t H = 15000.0;

   /// Width of VW-to-VS transition zone [m] (3 km)
   real_t h = 3000.0;

   /// Depth to bottom of rate-and-state fault [m] (40 km)
   real_t Wf = 40000.0;

   // =========================================================================
   // Simulation parameters
   // =========================================================================

   /// Seconds per year [s]
   static constexpr real_t seconds_per_year = 365.25 * 24.0 * 3600.0;

   /// Final simulation time [s] (1200 years)
   real_t t_final = 1200.0 * seconds_per_year;

   // =========================================================================
   // Methods
   // =========================================================================

   /// Compute the depth-dependent parameter a(z).
   ///
   /// Z-coordinate convention: z=0 at surface, z<0 at depth
   /// The depth magnitude is |z|:
   ///
   /// a(z) = a0                                      for 0 <= |z| < H
   ///      = a0 + (amax - a0) * (|z| - H) / h        for H <= |z| < H + h
   ///      = amax                                    for H + h <= |z| < Wf
   ///
   /// @param[in] z Depth coordinate [m] (z=0 at surface, z<0 at depth)
   /// @return Rate-and-state parameter a [-]
   real_t a_of_z(real_t z) const
   {
      // Use absolute value since z is negative for depth
      real_t depth = std::abs(z);

      if (depth < H)
      {
         return a0;
      }
      else if (depth < H + h)
      {
         return a0 + (amax - a0) * (depth - H) / h;
      }
      else
      {
         return amax;
      }
   }

   /// Check if a depth is in the velocity-weakening zone (a < b).
   ///
   /// @param[in] z Depth coordinate [m] (z=0 at surface, z<0 at depth)
   /// @return True if a(z) < b (velocity-weakening)
   bool IsVelocityWeakening(real_t z) const
   {
      return a_of_z(z) < b;
   }

   /// Compute the pre-stress tau0 (uniform across the fault).
   ///
   /// Pre-stress is computed at steady state with V = V_init at the
   /// velocity-strengthening zone (z >= H + h, where a = amax).
   ///
   /// tau0 = sigma_n * amax * asinh[(V_init / 2V0) * exp((f0 + b*ln(V0/V_init)) / amax)]
   ///        + eta * V_init
   ///
   /// @return Pre-stress [Pa]
   real_t tau0() const
   {
      // At steady state: theta_ss = Dc / V_init
      // The friction coefficient at steady state is:
      // f_ss = amax * asinh[(V_init / 2V0) * exp((f0 + b * ln(V0 * theta_ss / Dc)) / amax)]
      //      = amax * asinh[(V_init / 2V0) * exp((f0 + b * ln(V0 / V_init)) / amax)]

      real_t log_term = f0 + b * std::log(V0 / V_init);
      real_t exp_term = std::exp(log_term / amax);
      real_t arg = (V_init / (2.0 * V0)) * exp_term;
      real_t f_ss = amax * std::asinh(arg);

      return sigma_n * f_ss + eta() * V_init;
   }

   /// Output stations (depths) for BP2 benchmark [m].
   /// 12 stations at z = 0, -2.4, -4.8, -7.2, -9.6, -12, -14.4, -16.8, -19.2, -24, -28.8, -36 km
   /// (z=0 at surface, z<0 at depth)
   static void GetOutputDepths(Array<real_t> &depths)
   {
      depths.SetSize(12);
      depths[0] = 0.0;
      depths[1] = -2400.0;
      depths[2] = -4800.0;
      depths[3] = -7200.0;
      depths[4] = -9600.0;
      depths[5] = -12000.0;
      depths[6] = -14400.0;
      depths[7] = -16800.0;
      depths[8] = -19200.0;
      depths[9] = -24000.0;
      depths[10] = -28800.0;
      depths[11] = -36000.0;
   }

   /// Benchmark data depths (subset with available reference data) [m].
   /// (z=0 at surface, z<0 at depth)
   static void GetBenchmarkDepths(Array<real_t> &depths)
   {
      depths.SetSize(5);
      depths[0] = 0.0;
      depths[1] = -4800.0;
      depths[2] = -12000.0;
      depths[3] = -16800.0;
      depths[4] = -24000.0;
   }

   /// Print parameters to output stream.
   void Print(std::ostream &os = mfem::out) const
   {
      os << "BP2 Parameters:\n";
      os << "  Material:\n";
      os << "    rho = " << rho << " kg/m^3\n";
      os << "    cs  = " << cs << " m/s\n";
      os << "    mu  = " << mu() / 1e9 << " GPa\n";
      os << "  Friction:\n";
      os << "    V0  = " << V0 << " m/s\n";
      os << "    f0  = " << f0 << "\n";
      os << "    a0  = " << a0 << "\n";
      os << "    amax = " << amax << "\n";
      os << "    b   = " << b << "\n";
      os << "    Dc  = " << Dc << " m\n";
      os << "  Stress:\n";
      os << "    sigma_n = " << sigma_n / 1e6 << " MPa\n";
      os << "    eta     = " << eta() << " Pa*s/m\n";
      os << "    tau0    = " << tau0() / 1e6 << " MPa\n";
      os << "  Loading:\n";
      os << "    Vp     = " << Vp << " m/s\n";
      os << "    V_init = " << V_init << " m/s\n";
      os << "  Geometry:\n";
      os << "    H  = " << H / 1000.0 << " km\n";
      os << "    h  = " << h / 1000.0 << " km\n";
      os << "    Wf = " << Wf / 1000.0 << " km\n";
      os << "  Simulation:\n";
      os << "    t_final = " << t_final / seconds_per_year << " years\n";
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP2_PARAMS_HPP
