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
   const DenseMatrix&                    dof_basis,
   const std::function<real_t(real_t, real_t, real_t)>& mu_at_xyz)
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
         // Phase 10 (TPV31): forward the per-point mu lookup so the spec-p.7
         // mu(depth)/mu_ref amplitude scaling is applied (no-op when the
         // callback is empty or mu_ref_pa <= 0).
         return std::make_unique<InstantaneousOverstressCircular>(
            spatial::ResolveInstantaneousOverstressCircular(
               cfg.nucleation.instantaneous_circular, /*enabled=*/true,
               dof_coords_3d, dof_basis, mu_at_xyz));
   }

   MFEM_ABORT("MakeNucleation: unhandled cfg.nucleation.kind = "
              << static_cast<int>(cfg.nucleation.kind));
   return nullptr;
}

} // namespace seas
} // namespace mfem
