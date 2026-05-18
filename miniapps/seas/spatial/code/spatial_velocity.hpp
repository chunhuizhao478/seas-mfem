// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_velocity.hpp — Phase 2 of spatial_dynamic_rupture_plan.md rev-3.
//
// Provides a thin SpatialVelocityBundle that owns three DataField3D
// readers (Vp, Vs, density), three {Lambda,Mu,Rho}FromSidecar
// Coefficients, and exposes a MaterialField::MakeCoefficient(...)
// factory.
//
// Plan reference: §Phase 2 of spatial_dynamic_rupture_plan.md.
// Schema reference for [velocity] block:
//   safs/project_7.0_alternative/document/spatial_friction_config_schema.md

#ifndef MFEM_SEAS_SPATIAL_VELOCITY_HPP
#define MFEM_SEAS_SPATIAL_VELOCITY_HPP

#include "mfem.hpp"

#include "../../dynamic/heterogeneous_material.hpp"
#include "../../io/data_field_3d.hpp"
#include "../../io/material_coefficients.hpp"

#include "spatial_friction.hpp"   // VelocitySpec + VelocityModel

#include <memory>
#include <string>

namespace mfem
{
namespace seas
{
namespace spatial
{

/// Owns three DataField3D readers (Vp, Vs, density) plus three
/// {Lambda,Mu,Rho}FromSidecar Coefficients and lets the driver build a
/// MaterialField::MakeCoefficient(...) handle from them.  The driver
/// keeps the bundle alive for the lifetime of the MaterialField (the
/// Coefficient pointers stored in MaterialField are non-owning).
struct SpatialVelocityBundle
{
   std::unique_ptr<DataField3D>      vp_field;
   std::unique_ptr<DataField3D>      vs_field;
   std::unique_ptr<DataField3D>      rho_field;
   std::unique_ptr<LambdaFromSidecar> lambda_coef;
   std::unique_ptr<MuFromSidecar>     mu_coef;
   std::unique_ptr<RhoFromSidecar>    rho_coef;

   MaterialField MakeMaterialField() const
   {
      return MaterialField::MakeCoefficient(lambda_coef.get(),
                                            mu_coef.get(),
                                            rho_coef.get());
   }
};

/// Resolve the on-disk sidecar path for a (model, dataset_root) pair.
/// Rev-3 R-003 correction: filenames have NO mesh_tag suffix; one
/// sidecar per CVM model:
///
///   CVMH                 -> <root>/velocity/results/cvmh/velocity_safs.h5
///   CVMS_4_26_M01        -> <root>/velocity/results/cvm_s4.26.m01/velocity_safs.h5
///   MultiscaleStatewise  -> <root>/velocity/results/multiscale_statewise_cvm/velocity_safs.h5
///
/// If spec.override_path is non-empty, the function returns it verbatim
/// (bypasses model resolution).
std::string ResolveSpatialVelocitySidecarPath(const VelocitySpec& spec);

/// Load the bundle.  Performs file-existence + ContainsBBox check vs
/// pmesh on all three fields via FieldProjector::AbortContainmentFailure.
/// Aborts with a clean "file not found at <path>" message if the
/// resolved sidecar is missing.
SpatialVelocityBundle LoadSpatialVelocityBundle(const VelocitySpec& spec,
                                                ParMesh&            pmesh);

/// Serial-mesh overload for unit tests (the bbox check is done against
/// the supplied Mesh; the rest of the bundle construction is identical).
SpatialVelocityBundle LoadSpatialVelocityBundle(const VelocitySpec& spec,
                                                mfem::Mesh&         mesh);

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif // MFEM_SEAS_SPATIAL_VELOCITY_HPP
