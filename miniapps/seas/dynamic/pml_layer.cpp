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

#include "pml_layer.hpp"
#include <cmath>
#include <algorithm>

namespace mfem
{
namespace seas
{

// Damping direction matrices (plan Section 3.2.5, Eq. 16).
// D_x[c] = 1 if component c involves x-direction, 0 otherwise.
//                        SXX SYY SZZ SXY SYZ SXZ VX VY VZ
const int PMLLayer::Dx[NUM_STATE] = {1,  0,  0,  1,  0,  1,  1, 0, 0};
const int PMLLayer::Dy[NUM_STATE] = {0,  1,  0,  1,  1,  0,  0, 1, 0};
const int PMLLayer::Dz[NUM_STATE] = {0,  0,  1,  0,  1,  1,  0, 0, 1};

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
PMLLayer::PMLLayer(const Vector &x_min, const Vector &x_max,
                   real_t thickness, real_t cp,
                   real_t target_R, int dirs)
   : x_min_(x_min), x_max_(x_max), L_pml_(thickness), dirs_(dirs)
{
   MFEM_VERIFY(thickness > 0, "PML thickness must be positive, got " << thickness);
   MFEM_VERIFY(cp > 0, "P-wave speed must be positive, got " << cp);
   MFEM_VERIFY(target_R > 0 && target_R < 1,
               "Target reflection must be in (0,1), got " << target_R);
   MFEM_VERIFY(x_min.Size() == 3 && x_max.Size() == 3,
               "Domain bounds must be 3D vectors");

   // d_max = (3 * cp) / (2 * L_pml) * ln(1/R_0)   (Eq. 17)
   d_max_ = (3.0 * cp) / (2.0 * L_pml_) * std::log(1.0 / target_R);
}

// ---------------------------------------------------------------------------
// Cubic damping profile
// ---------------------------------------------------------------------------
real_t PMLLayer::DampingProfile(real_t dist) const
{
   if (dist <= 0.0) { return 0.0; }
   real_t ratio = dist / L_pml_;
   if (ratio > 1.0) { ratio = 1.0; }
   return d_max_ * ratio * ratio * ratio;  // cubic profile
}

// ---------------------------------------------------------------------------
// Compute directional damping at a physical point
// ---------------------------------------------------------------------------
void PMLLayer::ComputeDamping(real_t x, real_t y, real_t z,
                              real_t &dx, real_t &dy, real_t &dz) const
{
   dx = dy = dz = 0.0;

   // x-direction PML (near x_min and x_max boundaries)
   if (dirs_ & 1)
   {
      real_t dist_lo = (x_min_(0) + L_pml_) - x;  // distance INTO PML from inner boundary
      real_t dist_hi = x - (x_max_(0) - L_pml_);
      if (dist_lo > 0) { dx = DampingProfile(dist_lo); }
      if (dist_hi > 0) { dx = std::max(dx, DampingProfile(dist_hi)); }
   }

   // y-direction PML
   if (dirs_ & 2)
   {
      real_t dist_lo = (x_min_(1) + L_pml_) - y;
      real_t dist_hi = y - (x_max_(1) - L_pml_);
      if (dist_lo > 0) { dy = DampingProfile(dist_lo); }
      if (dist_hi > 0) { dy = std::max(dy, DampingProfile(dist_hi)); }
   }

   // z-direction PML
   if (dirs_ & 4)
   {
      real_t dist_lo = (x_min_(2) + L_pml_) - z;
      real_t dist_hi = z - (x_max_(2) - L_pml_);
      if (dist_lo > 0) { dz = DampingProfile(dist_lo); }
      if (dist_hi > 0) { dz = std::max(dz, DampingProfile(dist_hi)); }
   }
}

} // namespace seas
} // namespace mfem
