// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/nucleation_factory.hpp — Phase 7 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// `MakeNucleation` maps `cfg.nucleation.kind` (and the absent-[nucleation]
// case) to the correct `INucleationMethod` concrete, resolving its per-DOF
// targets from the per-DOF coords + (dip, strike) basis.

#ifndef MFEM_SEAS_NUCLEATION_FACTORY_HPP
#define MFEM_SEAS_NUCLEATION_FACTORY_HPP

#include "mfem.hpp"

#include "nucleation_method.hpp"
#include "../spatial/code/spatial_friction.hpp"   // SpatialFrictionConfig

#include <memory>

namespace mfem
{
namespace seas
{

/// @brief Build the nucleation method selected by `cfg.nucleation`.
///
/// `cfg.nucleation.enabled == false` (absent `[nucleation]` block, e.g.
/// TPV205) → `StaticOverstress`.  Otherwise dispatch on `cfg.nucleation.kind`:
///   GradualOverstress                → GaussianGradualOverstress
///   GradualOverstressCompactCircular → CompactCircularGradualOverstress
///   InstantaneousOverstressCircular  → InstantaneousOverstressCircular
///
/// DEVIATION (plan §5.2): the plan's signature takes
/// `const std::vector<Vector>&` / `const std::vector<FaultBasisRow>&`, but
/// `FaultBasisRow` does not exist in this codebase — the driver and the
/// existing resolvers use a flat `Vector dof_coords_3d` (3N) + a `(9, N)`
/// column-major `DenseMatrix dof_basis` (rows 0..2 normal, 3..5 dip, 6..8
/// strike).  We use those concrete types so MakeNucleation forwards directly
/// to the existing `Resolve*` overloads with no conversion.
std::unique_ptr<INucleationMethod> MakeNucleation(
   const spatial::SpatialFrictionConfig& cfg,
   const Vector&                         dof_coords_3d,
   const DenseMatrix&                    dof_basis);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NUCLEATION_FACTORY_HPP
