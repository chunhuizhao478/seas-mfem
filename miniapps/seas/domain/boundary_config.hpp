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
/// The function receives the spatial coordinate and current time,
/// and writes the prescribed displacement into the output vector.
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
   int fault_attr = 3;             ///< Boundary attr for fault faces

   /// Per-attribute Dirichlet functions. Key = attr number.
   /// Attrs not in this map use default_dirichlet_func.
   std::map<int, DirichletFunc> dirichlet_funcs;

   /// Default Dirichlet function applied to all Dirichlet attrs
   /// not listed in dirichlet_funcs.
   DirichletFunc default_dirichlet_func;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BOUNDARY_CONFIG_HPP
