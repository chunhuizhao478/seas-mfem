// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Implementation of MakeNucleation (Phase 7).

#include "nucleation_factory.hpp"

namespace mfem
{
namespace seas
{

std::unique_ptr<INucleationMethod> MakeNucleation(
   const spatial::SpatialFrictionConfig& cfg,
   const Vector&                         dof_coords_3d,
   const DenseMatrix&                    dof_basis)
{
   // Absent [nucleation] (enabled == false) → static (TPV205 path): the
   // pre-stress / patches already live in tau_pre_; nothing to perturb.
   if (!cfg.nucleation.enabled)
   {
      return std::make_unique<StaticOverstress>();
   }

   switch (cfg.nucleation.kind)
   {
      case spatial::NucleationKind::GradualOverstress:
         // SAFS + RS Gaussian path — byte-identical to the driver's pre-Phase-7
         // inline ResolveGradualOverstress + ApplyGradualOverstressIncrement.
         return std::make_unique<GaussianGradualOverstress>(
            spatial::ResolveGradualOverstress(
               cfg.nucleation.gradual_overstress, /*enabled=*/true,
               dof_coords_3d, dof_basis),
            cfg.nucleation.gradual_overstress.T_nuc_s);

      case spatial::NucleationKind::GradualOverstressCompactCircular:
         return std::make_unique<CompactCircularGradualOverstress>(
            spatial::ResolveGradualOverstressCompactCircular(
               cfg.nucleation.compact_circular, /*enabled=*/true,
               dof_coords_3d, dof_basis),
            cfg.nucleation.compact_circular.T_nuc_s);

      case spatial::NucleationKind::InstantaneousOverstressCircular:
         return std::make_unique<InstantaneousOverstressCircular>(
            spatial::ResolveInstantaneousOverstressCircular(
               cfg.nucleation.instantaneous_circular, /*enabled=*/true,
               dof_coords_3d, dof_basis));
   }

   MFEM_ABORT("MakeNucleation: unhandled cfg.nucleation.kind = "
              << static_cast<int>(cfg.nucleation.kind));
   return nullptr;
}

} // namespace seas
} // namespace mfem
