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

#ifndef MFEM_SEAS_FRICTION_LAW_HPP
#define MFEM_SEAS_FRICTION_LAW_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// Abstract base class for rate-and-state friction laws.
///
/// This class defines the interface for computing friction coefficients,
/// fault strength, and solving for slip rate given stress conditions.
/// Implementations should provide specific friction formulations such as
/// Dieterich-Ruina (regularized) or slip-weakening laws.
///
/// The friction law relates shear stress (tau) to slip rate (V) and state (theta):
///   tau = sigma_n * f(V, theta) + eta * V
///
/// where:
///   - sigma_n: normal stress
///   - f(V, theta): friction coefficient
///   - eta: radiation damping coefficient (mu / 2*cs for quasi-dynamic)
///   - V: slip rate
///   - theta: state variable
class FrictionLaw
{
public:
   virtual ~FrictionLaw() = default;

   /// Compute the friction coefficient f(V, theta).
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] a Rate-and-state direct effect parameter (depth-dependent)
   ///
   /// @note Implementation uses scalar 'a' instead of Vector params for
   ///       efficiency. For BP2, only 'a' varies spatially; other parameters
   ///       (b, Dc, V0, f0) are stored in the Constants struct.
   ///
   /// @return Friction coefficient (dimensionless)
   virtual real_t FrictionCoefficient(real_t V, real_t theta, real_t a) const = 0;

   /// Compute the fault strength F = sigma_n * f(V, theta).
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Fault strength [Pa]
   virtual real_t FaultStrength(real_t V, real_t theta, real_t sigma_n,
                                real_t a) const
   {
      return sigma_n * FrictionCoefficient(V, theta, a);
   }

   /// Solve for slip rate V given shear stress tau and state theta.
   ///
   /// Solves the equation: tau = sigma_n * f(V, theta) + eta * V
   ///
   /// @param[in] tau Total shear stress [Pa]
   /// @param[in] theta State variable [s]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] eta Radiation damping coefficient [Pa*s/m]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @param[out] iterations Optional: number of iterations used (nullptr to ignore)
   /// @return Slip rate V [m/s]
   virtual real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                                real_t eta, real_t a,
                                int *iterations = nullptr) const = 0;

   /// Compute initial state theta(0) for stress equilibrium.
   ///
   /// Given pre-stress tau0 and initial slip rate V_init, compute the
   /// state variable theta such that:
   ///   tau0 = sigma_n * f(V_init, theta) + eta * V_init
   ///
   /// @param[in] tau0 Pre-stress / initial shear stress [Pa]
   /// @param[in] V_init Initial slip rate [m/s]
   /// @param[in] sigma_n Normal stress [Pa]
   /// @param[in] eta Radiation damping coefficient [Pa*s/m]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Initial state variable theta [s]
   virtual real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                               real_t eta, real_t a) const = 0;

   /// Compute the derivative df/dV for Newton solver.
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Derivative of friction coefficient with respect to V [s/m]
   virtual real_t FrictionDerivativeV(real_t V, real_t theta, real_t a) const = 0;

   /// Compute the derivative df/dtheta for implicit time stepping.
   ///
   /// @param[in] V Slip rate [m/s]
   /// @param[in] theta State variable [s]
   /// @param[in] a Rate-and-state direct effect parameter
   /// @return Derivative of friction coefficient with respect to theta [1/s]
   virtual real_t FrictionDerivativeTheta(real_t V, real_t theta,
                                          real_t a) const = 0;

   /// Get the reference slip rate V0.
   virtual real_t GetV0() const = 0;

   /// Get the reference friction coefficient f0.
   virtual real_t GetF0() const = 0;

   /// Get the state evolution parameter b.
   virtual real_t GetB() const = 0;

   /// Get the critical slip distance Dc.
   virtual real_t GetDc() const = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_LAW_HPP
