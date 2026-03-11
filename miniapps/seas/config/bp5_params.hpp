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

#ifndef MFEM_SEAS_BP5_PARAMS_HPP
#define MFEM_SEAS_BP5_PARAMS_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// Parameters for the SCEC SEAS BP5-QD benchmark problem.
///
/// BP5 is a 3D full elasticity problem with vector slip on a 2D fault plane.
/// Reference: https://strike.scec.org/cvws/seas/download/SEAS_BP5_QD.pdf
///
/// Coordinate convention (SCEC):
///   x1 = fault-normal
///   x2 = along-strike
///   x3 = depth (positive downward)
///
/// All values are in SI units:
/// - Length: meters [m]
/// - Time: seconds [s]
/// - Stress: Pascals [Pa]
/// - Velocity: meters/second [m/s]
struct BP5Params
{
   // =========================================================================
   // Material properties
   // =========================================================================

   /// Rock density [kg/m^3]
   real_t rho = 2670.0;

   /// Shear wave speed [m/s]
   real_t cs = 3464.0;

   /// Poisson's ratio [-]
   real_t nu = 0.25;

   /// Shear modulus: mu = rho * cs^2 [Pa]
   real_t mu() const { return rho * cs * cs; }

   /// First Lame parameter: lambda = 2*nu*mu/(1-2*nu) [Pa]
   /// For nu=0.25: lambda = mu
   real_t lambda() const { return 2.0 * nu * mu() / (1.0 - 2.0 * nu); }

   /// Radiation damping coefficient: eta = mu / (2 * cs) [Pa*s/m]
   real_t eta() const { return mu() / (2.0 * cs); }

   // =========================================================================
   // Rate-and-state friction parameters
   // =========================================================================

   /// Reference slip rate [m/s]
   real_t V0 = 1.0e-6;

   /// Reference friction coefficient [-]
   real_t f0 = 0.6;

   /// State evolution effect parameter b [-]
   real_t b = 0.03;

   /// Critical slip distance outside nucleation zone [m]
   real_t L0 = 0.14;

   /// Critical slip distance in nucleation zone [m]
   real_t L_nuc = 0.13;

   /// Direct effect parameter in velocity-weakening zone [-]
   real_t a0 = 0.004;

   /// Direct effect parameter in velocity-strengthening zone [-]
   real_t amax = 0.04;

   // =========================================================================
   // Stress parameters
   // =========================================================================

   /// Effective normal stress [Pa] (25 MPa)
   real_t sigma_n = 25.0e6;

   // =========================================================================
   // Loading parameters
   // =========================================================================

   /// Plate rate (far-field loading velocity) [m/s]
   real_t Vp = 1.0e-9;

   /// Initial slip rate [m/s] (same as plate rate)
   real_t V_init = 1.0e-9;

   /// Minimum velocity to avoid log(0) [m/s]
   real_t V_zero = 1.0e-20;

   /// Nucleation zone initial slip rate [m/s]
   /// SCEC BP5 spec Section 3: V_i = 0.03 m/s in the favorable nucleation zone.
   /// Same for both BP5-QD and BP5-FD; only δτ in pre-stress differs.
   /// (Tandem's bp5.lua uses 0.01, which deviates from the SCEC standard.)
   real_t V_nuc = 0.03;

   // =========================================================================
   // Geometric parameters (all in meters)
   // =========================================================================

   /// Width of shallow velocity-strengthening zone [m] (2 km)
   real_t hs = 2.0e3;

   /// Width of VW-VS transition zone [m] (2 km)
   real_t ht = 2.0e3;

   /// Width of uniform velocity-weakening region [m] (12 km)
   real_t H = 12.0e3;

   /// Length of uniform velocity-weakening region [m] (60 km)
   real_t l_vw = 60.0e3;

   /// Width of favorable nucleation zone [m] (12 km)
   real_t w_nuc = 12.0e3;

   /// Depth to bottom of rate-and-state fault [m] (40 km)
   real_t Wf = 40.0e3;

   /// Length of rate-and-state fault [m] (100 km)
   real_t lf = 100.0e3;

   // =========================================================================
   // Simulation parameters
   // =========================================================================

   /// Seconds per year [s]
   static constexpr real_t seconds_per_year = 365.25 * 24.0 * 3600.0;

   /// Final simulation time [s] (1800 years)
   real_t t_final = 1800.0 * seconds_per_year;

   // =========================================================================
   // Spatial parameter functions
   // =========================================================================

   /// Compute the 2D rate-and-state parameter a(x2, x3).
   ///
   /// BP5 spec Eq. 14:
   ///   a = a0      in VW core:  (hs+ht <= x3 <= hs+ht+H) AND (|x2| <= l/2)
   ///   a = amax    in VS zones: x3 <= hs OR x3 >= hs+2*ht+H OR |x2| >= l/2+ht
   ///   a = a0 + r*(amax-a0) in transition, where
   ///       r = max(|x3 - hs - ht - H/2| - H/2, |x2| - l/2) / ht
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m] (positive downward, range [0, Wf])
   /// @return Rate-and-state parameter a [-]
   real_t a_of_x2_x3(real_t x2, real_t x3) const
   {
      real_t abs_x2 = std::abs(x2);
      real_t half_l = l_vw / 2.0;

      // VW core: (hs+ht <= x3 <= hs+ht+H) AND (|x2| <= l/2)
      if (x3 >= hs + ht && x3 <= hs + ht + H && abs_x2 <= half_l)
      {
         return a0;
      }

      // VS zones: x3 <= hs OR x3 >= hs+2*ht+H OR |x2| >= l/2+ht
      if (x3 <= hs || x3 >= hs + 2.0 * ht + H || abs_x2 >= half_l + ht)
      {
         return amax;
      }

      // Transition zone
      real_t r = std::max(std::abs(x3 - hs - ht - H / 2.0) - H / 2.0,
                          abs_x2 - half_l) / ht;
      return a0 + r * (amax - a0);
   }

   /// Check if a point is in the favorable nucleation zone.
   ///
   /// Nucleation zone: (hs+ht <= x3 <= hs+ht+H) AND (-l/2 <= x2 <= -l/2+w)
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m] (positive downward)
   /// @return True if in nucleation zone
   bool IsNucleationZone(real_t x2, real_t x3) const
   {
      real_t half_l = l_vw / 2.0;
      return (x3 >= hs + ht && x3 <= hs + ht + H &&
              x2 >= -half_l && x2 <= -half_l + w_nuc);
   }

   /// Compute spatially varying critical slip distance L(x2, x3).
   ///
   /// L = L_nuc (0.13 m) in nucleation zone, L0 (0.14 m) elsewhere.
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m]
   /// @return Critical slip distance [m]
   real_t L_of_x2_x3(real_t x2, real_t x3) const
   {
      if (IsNucleationZone(x2, x3))
      {
         return L_nuc;
      }
      return L0;
   }

   /// Alias for L_of_x2_x3() using friction law naming convention.
   ///
   /// The DieterichRuinaFriction::Constants struct uses `Dc` (critical slip
   /// distance). BP5 calls the same quantity `L`. This alias reduces confusion
   /// when initializing per-DOF friction constants from BP5 parameters.
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m]
   /// @return Critical slip distance Dc [m]
   real_t Dc_of_x2_x3(real_t x2, real_t x3) const
   {
      return L_of_x2_x3(x2, x3);
   }

   /// Compute initial velocity vector V_init(x2, x3).
   ///
   /// Following Tandem's bp5.lua convention:
   /// - Nucleation zone: V = (V_zero, V_nuc)
   /// - Elsewhere: V = (V_zero, V_init)
   ///
   /// Component ordering convention (fault-local frame):
   ///   V[0] = along-dip ≈ 0 (V_zero placeholder)
   ///   V[1] = along-strike = dominant slip component
   ///
   /// This matches Tandem's bp5.lua where tangent1=dip, tangent2=strike.
   /// Note: SCEC Eq. 16 writes V=(V_strike, V_dip)=(dominant, ≈0), which
   /// is the opposite ordering. The FaultBasis class handles the mapping
   /// between this fault-local frame and global coordinates.
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m]
   /// @param[out] V Initial velocity vector [m/s] (2 components)
   void V_init_vec(real_t x2, real_t x3, real_t V[2]) const
   {
      V[0] = V_zero;
      if (IsNucleationZone(x2, x3))
      {
         V[1] = V_nuc;
      }
      else
      {
         V[1] = V_init;
      }
   }

   /// Compute pre-stress vector tau0(x2, x3) for initial conditions.
   ///
   /// SCEC BP5 spec Eq. 20 (outside nucleation zone):
   ///   tau0 = sigma_n * a * asinh(V_init/(2V0) * exp(psi_ss/a)) + eta*V_init
   ///
   /// SCEC BP5 spec Eq. 23 (inside nucleation zone):
   ///   tau0_i = sigma_n * a * asinh(V_i/(2V0) * exp(psi_ss/a)) + eta*V_i + delta_tau
   ///   where delta_tau = eta*V_i for BP5-QD, delta_tau = 0 for BP5-FD
   ///
   /// The extra delta_tau in QD means the nucleation zone is initially
   /// overstressed, driving V above V_i and accelerating nucleation.
   ///
   /// Pre-stress direction is parallel to initial velocity.
   ///
   /// @param[in] x2 Along-strike coordinate [m]
   /// @param[in] x3 Depth coordinate [m]
   /// @param[out] tau Pre-stress vector [Pa] (2 components)
   void tau0_vec(real_t x2, real_t x3, real_t tau[2]) const
   {
      real_t Vi[2];
      V_init_vec(x2, x3, Vi);

      real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);

      real_t a = a_of_x2_x3(x2, x3);
      real_t eta_val = eta();

      // Steady-state psi at plate rate
      real_t psi_ss = f0 + b * std::log(V0 / Vp);

      // Scalar pre-stress magnitude (using velocity magnitude)
      real_t e = std::exp(psi_ss / a);
      real_t tau0_scalar = sigma_n * a *
                           std::asinh((Vi_abs / (2.0 * V0)) * e) +
                           eta_val * Vi_abs;

      // BP5-QD: add delta_tau = eta * V_i in nucleation zone (Eq. 23)
      if (IsNucleationZone(x2, x3))
      {
         tau0_scalar += eta_val * Vi_abs;
      }

      // Direction: parallel to initial velocity
      tau[0] = tau0_scalar * Vi[0] / Vi_abs;
      tau[1] = tau0_scalar * Vi[1] / Vi_abs;
   }

   /// Print parameters to output stream.
   void Print(std::ostream &os = mfem::out) const
   {
      os << "BP5 Parameters:\n";
      os << "  Material:\n";
      os << "    rho    = " << rho << " kg/m^3\n";
      os << "    cs     = " << cs << " m/s\n";
      os << "    nu     = " << nu << "\n";
      os << "    mu     = " << mu() / 1e9 << " GPa\n";
      os << "    lambda = " << lambda() / 1e9 << " GPa\n";
      os << "  Friction:\n";
      os << "    V0     = " << V0 << " m/s\n";
      os << "    f0     = " << f0 << "\n";
      os << "    a0     = " << a0 << "\n";
      os << "    amax   = " << amax << "\n";
      os << "    b      = " << b << "\n";
      os << "    L0     = " << L0 << " m\n";
      os << "    L_nuc  = " << L_nuc << " m\n";
      os << "  Stress:\n";
      os << "    sigma_n = " << sigma_n / 1e6 << " MPa\n";
      os << "    eta     = " << eta() << " Pa*s/m\n";
      os << "  Loading:\n";
      os << "    Vp     = " << Vp << " m/s\n";
      os << "    V_init = " << V_init << " m/s\n";
      os << "    V_nuc  = " << V_nuc << " m/s\n";
      os << "  Geometry:\n";
      os << "    hs     = " << hs / 1000.0 << " km\n";
      os << "    ht     = " << ht / 1000.0 << " km\n";
      os << "    H      = " << H / 1000.0 << " km\n";
      os << "    l_vw   = " << l_vw / 1000.0 << " km\n";
      os << "    w_nuc  = " << w_nuc / 1000.0 << " km\n";
      os << "    Wf     = " << Wf / 1000.0 << " km\n";
      os << "    lf     = " << lf / 1000.0 << " km\n";
      os << "  Simulation:\n";
      os << "    t_final = " << t_final / seconds_per_year << " years\n";
   }

   /// Validate physical consistency of parameters.
   void Validate() const
   {
      MFEM_VERIFY(a0 < b, "BP5: a0 must be < b for velocity-weakening");
      MFEM_VERIFY(amax > b, "BP5: amax must be > b for velocity-strengthening");
      MFEM_VERIFY(L_nuc < L0, "BP5: L_nuc should be < L0");
      MFEM_VERIFY(sigma_n > 0, "BP5: sigma_n must be positive");
      MFEM_VERIFY(Vp > 0, "BP5: Vp must be positive");
      MFEM_VERIFY(V_nuc > V_init, "BP5: V_nuc should exceed V_init");
   }

   /// Steady-state psi at plate rate (initial state variable, same everywhere).
   /// SCEC: theta(0) = L/V_init → psi = f0 + b*ln(V0/V_init)
   real_t psi_init() const { return f0 + b * std::log(V0 / V_init); }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_PARAMS_HPP
