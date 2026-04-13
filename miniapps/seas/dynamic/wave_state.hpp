// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_WAVE_STATE_HPP
#define MFEM_SEAS_WAVE_STATE_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// Number of state components in the 3D velocity-stress system.
static constexpr int NUM_STATE = 9;

/// State vector index convention for Q = [sigma, v].
/// Matches SeisSol ordering (ElasticSetup.h) and the plan Section 2.1.1.
enum QIndex : int
{
   SXX = 0,  ///< Normal stress sigma_xx
   SYY = 1,  ///< Normal stress sigma_yy
   SZZ = 2,  ///< Normal stress sigma_zz
   SXY = 3,  ///< Shear stress sigma_xy
   SYZ = 4,  ///< Shear stress sigma_yz
   SXZ = 5,  ///< Shear stress sigma_xz
   VX  = 6,  ///< Particle velocity v_x
   VY  = 7,  ///< Particle velocity v_y
   VZ  = 8   ///< Particle velocity v_z
};

/// Compute elastic energy density: E = (1/2) sigma : epsilon + (1/2) rho |v|^2.
///
/// For isotropic linear elasticity with Lame parameters lambda, mu:
///   E_strain = (1/(4*mu)) * (sigma:sigma - lambda/(3*lambda+2*mu) * (tr sigma)^2)
///   E_kinetic = (1/2) * rho * |v|^2
///
/// @param[in] Q  State vector (9 components at a point)
/// @param[in] lambda  First Lame parameter [Pa]
/// @param[in] mu  Shear modulus [Pa]
/// @param[in] rho  Density [kg/m^3]
/// @return Total energy density [J/m^3]
inline real_t EnergyDensity(const real_t *Q,
                            real_t lambda, real_t mu, real_t rho)
{
   // Strain energy: E_strain = sigma:epsilon / 2
   //   epsilon = (1/(2*mu)) * (sigma - (lambda/(3*lambda+2*mu)) * tr(sigma) * I)
   // So E_strain = (1/(4*mu)) * (|sigma|^2 - lambda/(3*lambda+2*mu) * tr(sigma)^2)
   real_t tr_sig = Q[SXX] + Q[SYY] + Q[SZZ];
   real_t sig_sq = Q[SXX]*Q[SXX] + Q[SYY]*Q[SYY] + Q[SZZ]*Q[SZZ]
                 + 2.0*(Q[SXY]*Q[SXY] + Q[SYZ]*Q[SYZ] + Q[SXZ]*Q[SXZ]);
   real_t inv4mu = 0.25 / mu;
   real_t E_strain = inv4mu * (sig_sq - lambda / (3.0*lambda + 2.0*mu) * tr_sig*tr_sig);

   // Kinetic energy
   real_t E_kinetic = 0.5 * rho * (Q[VX]*Q[VX] + Q[VY]*Q[VY] + Q[VZ]*Q[VZ]);

   return E_strain + E_kinetic;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_WAVE_STATE_HPP
