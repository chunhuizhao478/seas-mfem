// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_stress.cpp — Phase 3 + Phase 3b implementation.

#include "spatial_stress.hpp"

#include "../../fault/fault_geometry_safs.inl"
// The templated ComputeParams<StressSource> overload (Phase 3b) is
// defined in this companion inline file:
#include "../../fault/fault_geometry_safs_templated.inl"

#include <limits>

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
//  Phase 3: ApplyCsmStressSidecar
// =====================================================================

namespace
{

template <typename MeshT>
void apply_csm_impl(const StressSpec& spec, FaultGeometry<MeshT>& geom
#ifdef MFEM_USE_MPI
                    , MPI_Comm node_comm = MPI_COMM_NULL
#endif
                    )
{
   MFEM_VERIFY(spec.kind == StressSourceKind::SidecarHDF5,
               "ApplyCsmStressSidecar: StressSpec.kind must be SidecarHDF5");
   MFEM_VERIFY(!spec.sidecar_path.empty(),
               "ApplyCsmStressSidecar: spec.sidecar_path is empty");

   // Reject BP2 / antiplane FaultGeometry up front (before any accessor
   // that could read uninitialised per-DOF state).  The 3-D / BP5 ctor
   // is the only one that populates per-DOF coords, basis, and zero-
   // normal counts; ComputeParams<StressField3D> later re-checks
   // is_bp5_ but failing there mid-projection produces a less precise
   // diagnostic.
   MFEM_VERIFY(geom.IsBP5(),
               "ApplyCsmStressSidecar: requires the 3-D / BP5 "
               "FaultGeometry ctor.  The BP2 / antiplane ctor cannot "
               "be used here.");

   // R-001 / SAFS-mode pre-flight: any zero-normal fallback would
   // silently project to zero.  Mirror the operator-side guard message.
   MFEM_VERIFY(geom.NumZeroNormalFallbacks() == 0,
               "ApplyCsmStressSidecar: FaultGeometry reports "
               << geom.NumZeroNormalFallbacks()
               << " zero-normal fallback DOF(s); the SAFS-mode "
               "projection would silently produce sigma_n = 0 "
               "at those DOFs.  Fix the fault-DOF basis (see "
               "ElasticityDomainOperator::GetFaultDOFBasis) before "
               "applying the CSM sidecar.");

   // The StressField3D is TRANSIENT: it destructs (and, in shared-
   // memory mode, frees its six node-shared windows) when this function
   // returns — the per-DOF projection below is its only consumer.
#ifdef MFEM_USE_MPI
   StressField3D field(spec.sidecar_path, OOBPolicy::Abort, node_comm);
#else
   StressField3D field(spec.sidecar_path);
#endif
   geom.ComputeParams(field,
                          spec.pore_pressure.P_p_pa,
                          spec.pore_pressure.P_p_grad_pa_per_m,
                          spec.pore_pressure.min_sigma_n_pa);
   MFEM_VERIFY(geom.HasParams(),
               "ApplyCsmStressSidecar: ComputeParams did not "
               "set HasParams() = true");
}

}  // namespace

void ApplyCsmStressSidecar(const StressSpec&       spec,
                           FaultGeometry<ParMesh>& geom)
{
   apply_csm_impl<ParMesh>(spec, geom);
}

#ifdef MFEM_USE_MPI
void ApplyCsmStressSidecar(const StressSpec&       spec,
                           FaultGeometry<ParMesh>& geom,
                           MPI_Comm                node_comm)
{
   apply_csm_impl<ParMesh>(spec, geom, node_comm);
}
#endif

void ApplyCsmStressSidecar(const StressSpec&         spec,
                           FaultGeometry<mfem::Mesh>& geom)
{
   apply_csm_impl<mfem::Mesh>(spec, geom);
}

// =====================================================================
//  Phase 3b: ConstantTensorStressSource
// =====================================================================

ConstantTensorStressSource::ConstantTensorStressSource(
   real_t sxx, real_t syy, real_t szz,
   real_t sxy, real_t syz, real_t sxz)
{
   sigma_.SetSize(3, 3);
   sigma_(0, 0) = sxx;  sigma_(1, 1) = syy;  sigma_(2, 2) = szz;
   sigma_(0, 1) = sxy;  sigma_(1, 0) = sxy;
   sigma_(1, 2) = syz;  sigma_(2, 1) = syz;
   sigma_(0, 2) = sxz;  sigma_(2, 0) = sxz;

   constexpr real_t inf = std::numeric_limits<real_t>::infinity();
   bbox_ = { -inf, inf, -inf, inf, -inf, inf };
}

mfem::DenseMatrix ConstantTensorStressSource::Evaluate(real_t /*x*/,
                                                       real_t /*y*/,
                                                       real_t /*z*/) const
{
   return sigma_;
}

// =====================================================================
//  Phase 10 (TPV31): DepthProportionalToShearModulusStressSource
// =====================================================================

DepthProportionalToShearModulusStressSource::
DepthProportionalToShearModulusStressSource(
   real_t sxx_per_mu, real_t syy_per_mu, real_t szz_per_mu,
   real_t sxy_per_mu, real_t syz_per_mu, real_t sxz_per_mu,
   real_t mu_ref_pa,
   MuAtFn mu_at_xyz)
   : sxx_per_mu_(sxx_per_mu), syy_per_mu_(syy_per_mu),
     szz_per_mu_(szz_per_mu),
     sxy_per_mu_(sxy_per_mu), syz_per_mu_(syz_per_mu),
     sxz_per_mu_(sxz_per_mu),
     mu_ref_pa_(mu_ref_pa),
     mu_at_xyz_(std::move(mu_at_xyz))
{
   MFEM_VERIFY(mu_ref_pa_ > 0.0,
               "DepthProportionalToShearModulusStressSource: mu_ref_pa "
               "must be > 0; got " << mu_ref_pa_);
   MFEM_VERIFY(static_cast<bool>(mu_at_xyz_),
               "DepthProportionalToShearModulusStressSource: mu_at_xyz "
               "callback must be non-null.");
   constexpr real_t inf = std::numeric_limits<real_t>::infinity();
   bbox_ = { -inf, inf, -inf, inf, -inf, inf };
}

mfem::DenseMatrix DepthProportionalToShearModulusStressSource::Evaluate(
   real_t x, real_t y, real_t z) const
{
   const real_t mu_local = mu_at_xyz_(x, y, z);
   MFEM_VERIFY(mu_local > 0.0,
               "DepthProportionalToShearModulusStressSource: mu_at_xyz("
               << x << ", " << y << ", " << z << ") returned "
               << mu_local << " (must be > 0).");
   const real_t scale = mu_local / mu_ref_pa_;
   mfem::DenseMatrix sigma(3, 3);
   sigma(0, 0) = sxx_per_mu_ * scale;
   sigma(1, 1) = syy_per_mu_ * scale;
   sigma(2, 2) = szz_per_mu_ * scale;
   sigma(0, 1) = sxy_per_mu_ * scale;  sigma(1, 0) = sigma(0, 1);
   sigma(1, 2) = syz_per_mu_ * scale;  sigma(2, 1) = sigma(1, 2);
   sigma(0, 2) = sxz_per_mu_ * scale;  sigma(2, 0) = sigma(0, 2);
   return sigma;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
