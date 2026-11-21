// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cfl.hpp — stable step of the central difference scheme with stiffness-
// proportional (Kelvin-Voigt) damping zeta:
//
//   dt_crit = (2 / omega_max) (sqrt(1 + xi^2) - xi),   xi = zeta omega_max / 2
//
// (Belytschko, Liu & Moran, Nonlinear Finite Elements, eq. 6.3.36, with the
// damping ratio of the highest mode xi = zeta omega / 2 for C = zeta K).
// omega_max^2 is the largest eigenvalue of M^{-1} K. Damage only lowers it.

#ifndef MFEM_PFF_CFL_HPP
#define MFEM_PFF_CFL_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace pff
{

inline real_t CriticalStep(real_t omega_max, real_t zeta)
{
   const real_t xi = 0.5 * zeta * omega_max;
   return 2.0 / omega_max * (std::sqrt(1.0 + xi * xi) - xi);
}

} // namespace pff
} // namespace mfem

#endif
