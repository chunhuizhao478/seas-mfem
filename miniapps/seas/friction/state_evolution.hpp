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

#ifndef MFEM_SEAS_STATE_EVOLUTION_HPP
#define MFEM_SEAS_STATE_EVOLUTION_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// Abstract base class for state evolution laws.
///
/// State evolution laws describe how the state variable theta evolves
/// with slip rate V and time. The general form is:
///   dtheta/dt = G(V, theta, Dc)
///
/// where Dc is the critical slip distance.
class StateEvolution
{
public:
   virtual ~StateEvolution() = default;

   /// Compute the state evolution rate dtheta/dt.
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] Dc Critical slip distance [m]
   /// @return Rate of change of state variable [dimensionless]
   virtual real_t Rate(real_t V, real_t theta, real_t Dc) const = 0;

   /// Compute the steady-state value of theta for a given slip rate.
   ///
   /// At steady state, dtheta/dt = 0, so theta_ss = f(V, Dc).
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] Dc Critical slip distance [m]
   /// @return Steady-state value of theta [s]
   virtual real_t SteadyState(real_t V, real_t Dc) const = 0;

   /// Compute the derivative dG/dV for implicit time stepping.
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] Dc Critical slip distance [m]
   /// @return Derivative of Rate with respect to V
   virtual real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const = 0;

   /// Compute the derivative dG/dtheta for implicit time stepping.
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] Dc Critical slip distance [m]
   /// @return Derivative of Rate with respect to theta
   virtual real_t RateDerivativeTheta(real_t V, real_t theta,
                                       real_t Dc) const = 0;

   /// Get the name of this evolution law for output/logging.
   virtual const char *GetName() const = 0;
};

/// Aging law: dtheta/dt = 1 - V*theta/Dc
///
/// The aging law describes time-dependent healing of the fault surface.
/// When the fault is locked (V=0), theta increases linearly with time.
/// At steady state (V*theta = Dc), dtheta/dt = 0.
///
/// Physical interpretation:
/// - theta represents the average age of contacts on the fault
/// - Contacts heal (theta increases) when the fault is stationary
/// - Contacts are renewed (theta decreases) when slip occurs
class AgingLaw : public StateEvolution
{
public:
   /// Compute dtheta/dt = 1 - V*theta/Dc.
   real_t Rate(real_t V, real_t theta, real_t Dc) const override
   {
      return 1.0 - V * theta / Dc;
   }

   /// Compute steady-state theta = Dc/V.
   real_t SteadyState(real_t V, real_t Dc) const override
   {
      MFEM_ASSERT(V > 0.0, "Slip rate must be positive for steady state");
      return Dc / V;
   }

   /// Compute dG/dV = -theta/Dc.
   real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const override
   {
      return -theta / Dc;
   }

   /// Compute dG/dtheta = -V/Dc.
   real_t RateDerivativeTheta(real_t V, real_t theta, real_t Dc) const override
   {
      return -V / Dc;
   }

   const char *GetName() const override { return "AgingLaw"; }
};

/// Slip law: dtheta/dt = -V*theta/Dc * ln(V*theta/Dc)
///
/// The slip law describes slip-dependent evolution of the state variable.
/// Evolution only occurs during slip (no healing when V=0).
///
/// Physical interpretation:
/// - State evolves proportionally to slip rate
/// - No time-dependent healing when the fault is locked
class SlipLaw : public StateEvolution
{
public:
   /// Compute dtheta/dt = -V*theta/Dc * ln(V*theta/Dc).
   real_t Rate(real_t V, real_t theta, real_t Dc) const override
   {
      real_t x = V * theta / Dc;
      // Avoid log(0) - use a small positive value
      if (x < 1e-50)
      {
         x = 1e-50;
      }
      return -x * std::log(x);
   }

   /// Compute steady-state theta = Dc/V.
   /// Note: Same as aging law at steady state.
   real_t SteadyState(real_t V, real_t Dc) const override
   {
      MFEM_ASSERT(V > 0.0, "Slip rate must be positive for steady state");
      return Dc / V;
   }

   /// Compute dG/dV = -(theta/Dc) * (ln(V*theta/Dc) + 1).
   real_t RateDerivativeV(real_t V, real_t theta, real_t Dc) const override
   {
      real_t x = V * theta / Dc;
      if (x < 1e-50)
      {
         x = 1e-50;
      }
      return -(theta / Dc) * (std::log(x) + 1.0);
   }

   /// Compute dG/dtheta = -(V/Dc) * (ln(V*theta/Dc) + 1).
   real_t RateDerivativeTheta(real_t V, real_t theta, real_t Dc) const override
   {
      real_t x = V * theta / Dc;
      if (x < 1e-50)
      {
         x = 1e-50;
      }
      return -(V / Dc) * (std::log(x) + 1.0);
   }

   const char *GetName() const override { return "SlipLaw"; }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STATE_EVOLUTION_HPP
