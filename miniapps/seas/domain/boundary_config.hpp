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

#ifndef MFEM_SEAS_BOUNDARY_CONFIG_HPP
#define MFEM_SEAS_BOUNDARY_CONFIG_HPP

#include "mfem.hpp"
#include <functional>
#include <map>
#include <set>

namespace mfem
{
namespace seas
{

/// Dirichlet boundary function: (x, time) → displacement vector.
///
/// DirichletFunc captures parameters BY VALUE at construction.
/// Changing Vp later requires creating a new lambda.
///
/// The function receives the spatial coordinate and current time,
/// and writes the prescribed displacement into the output vector.
/// The output vector is pre-sized to dim (3 for 3D).
using DirichletFunc = std::function<void(const Vector &x, real_t t, Vector &u)>;

/// @brief Boundary configuration for the elasticity domain operator.
///
/// Specifies which mesh boundary attributes correspond to Dirichlet,
/// Natural, and fault boundaries, plus the Dirichlet loading functions.
///
/// No hardcoded presets — the user always provides attribute numbers
/// explicitly (via TOML config or constructor arguments).
struct BoundaryConfig
{
   std::set<int> dirichlet_attrs;  ///< Boundary attrs with Dirichlet BC
   std::set<int> natural_attrs;    ///< Boundary attrs with Natural (zero-traction) BC
   std::set<int> absorbing_attrs;  ///< Boundary attrs with first-order absorbing BC
   int fault_attr = 3;             ///< Boundary attr for fault faces

   /// Per-attribute Dirichlet functions. Key = attr number.
   /// Attrs not in this map use default_dirichlet_func.
   std::map<int, DirichletFunc> dirichlet_funcs;

   /// Default Dirichlet function applied to all Dirichlet attrs
   /// not listed in dirichlet_funcs.
   DirichletFunc default_dirichlet_func;
};

/// Create the standard BP5 Dirichlet loading function.
///
/// Tandem bp5.lua boundary(x,y,z,t):
///   y > 1000 m:   u_D = ( Vp*t/2, 0, 0)
///   y < -1000 m:  u_D = (-Vp*t/2, 0, 0)
///   |y| <= 1000:  u_D = ( Vp*t,   0, 0)
///
/// Vp is captured by value.
inline DirichletFunc MakeBP5DirichletFunc(real_t Vp)
{
   return [Vp](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(3);
      u = 0.0;
      real_t y = x(1);
      real_t Vh = Vp * t;
      if (y > 1000.0)       { Vh *= 0.5; }
      else if (y < -1000.0) { Vh *= -0.5; }
      u(0) = Vh;
   };
}

/// Create the SAF far-field plate-loading function for a fault box whose
/// y-coordinates are NOT centred on 0 (e.g. a UTM mesh, y ~ 3.7e6 m).
///
/// The BP5 function above keys the +/-Vp/2 antisymmetry off the ABSOLUTE y with
/// a +/-1000 m threshold (BP5's fault sits at y=0).  On a non-origin-centred
/// mesh every y exceeds +1000, so BP5 degenerates to UNIFORM +Vp/2 -> no
/// differential plate motion (REVIEW R-001).  This function instead measures y
/// relative to the box mid-plane y0 and imposes an antisymmetric strike-parallel
/// (x-component) displacement:
///   yr = y - y0 >= 0 ("front") : u_D = ( Vp*t/2, 0, 0)
///   yr = y - y0 <  0 ("back")  : u_D = (-Vp*t/2, 0, 0)
/// i.e. per-side +/-Vp/2, relative plate rate Vp across the box.
///
/// Planar approximation: this loads u_X, which equals fault-strike-parallel
/// motion only when the box x-axis is (approximately) the fault strike — the
/// assumption the user accepted for this SAF setup.  Vp and y0 are captured by
/// value (y0 = box mid-plane, computed from the mesh bounding box by the caller).
inline DirichletFunc MakeSAFDirichletFunc(real_t Vp, real_t y0)
{
   return [Vp, y0](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(3);
      u = 0.0;
      const real_t yr = x(1) - y0;
      real_t Vh = Vp * t;
      if (yr >= 0.0) { Vh *= 0.5; }
      else           { Vh *= -0.5; }
      u(0) = Vh;
   };
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BOUNDARY_CONFIG_HPP
