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

#ifndef MFEM_SEAS_PML_LAYER_HPP
#define MFEM_SEAS_PML_LAYER_HPP

#include "mfem.hpp"
#include "wave_state.hpp"

namespace mfem
{
namespace seas
{

/// @brief Convolutional Perfectly Matched Layer (CPML) for absorbing boundaries.
///
/// Implements Eq. (16)-(17) from the plan: adds a damping term -d(x) D Q
/// to the wave equation RHS. The damping coefficient d(x) varies cubically
/// from zero at the inner PML boundary to d_max at the outer domain boundary.
///
/// The damping is directional: each PML face absorbs waves in its normal
/// direction only. At corners where PML regions overlap, the per-component
/// damping is d_c = d_x*D_x[c] + d_y*D_y[c] + d_z*D_z[c] (R-004 fix).
///
/// Reference: Komatitsch & Martin (2007), unsplit convolutional PML.
class PMLLayer
{
public:
   /// @brief Construct a PML layer.
   ///
   /// @param[in] x_min  Lower corner of the computational domain.
   /// @param[in] x_max  Upper corner of the computational domain.
   /// @param[in] thickness  PML layer thickness [m] (5-10 km for BP5).
   /// @param[in] cp  P-wave speed [m/s] (for computing d_max).
   /// @param[in] target_R  Target reflection coefficient (default 1e-3).
   /// @param[in] dirs  Which directions to apply PML (bitmask: 1=x, 2=y, 4=z).
   ///                  Default 7 = all directions.
   PMLLayer(const Vector &x_min, const Vector &x_max,
            real_t thickness, real_t cp,
            real_t target_R = 1e-3, int dirs = 7);

   /// @brief Compute the PML damping at a physical point.
   ///
   /// Returns three directional damping values (d_x, d_y, d_z).
   /// At corners, multiple directions may have nonzero damping.
   ///
   /// @param[in] x, y, z  Physical coordinates.
   /// @param[out] dx  Damping coefficient in x-direction.
   /// @param[out] dy  Damping coefficient in y-direction.
   /// @param[out] dz  Damping coefficient in z-direction.
   void ComputeDamping(real_t x, real_t y, real_t z,
                       real_t &dx, real_t &dy, real_t &dz) const;

   /// @name Accessors
   ///@{
   real_t GetThickness() const { return L_pml_; }
   real_t GetDmax() const { return d_max_; }
   ///@}

   /// Damping direction matrices D_x, D_y, D_z.
   /// D_x[c] = 1 if component c has an x-index, 0 otherwise.
   /// E.g., D_x = {1, 0, 0, 1, 0, 1, 1, 0, 0} for
   ///       {SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ}.
   static const int Dx[NUM_STATE];
   static const int Dy[NUM_STATE];
   static const int Dz[NUM_STATE];

private:
   Vector x_min_, x_max_;
   real_t L_pml_;
   real_t d_max_;
   int dirs_;

   /// Cubic damping profile: d(s) = d_max * (s / L_pml)^3
   /// where s = distance from inner PML boundary.
   real_t DampingProfile(real_t dist) const;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PML_LAYER_HPP
