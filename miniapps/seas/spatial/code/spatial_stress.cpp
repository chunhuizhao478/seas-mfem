// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_stress.cpp — Phase 3 + Phase 3b implementation.

#include "spatial_stress.hpp"

#include "../../fault/fault_geometry_safs.inl"
// The templated ComputeSAFSParams<StressSource> overload (Phase 3b) is
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
void apply_csm_impl(const StressSpec& spec, FaultGeometry<MeshT>& geom)
{
   MFEM_VERIFY(spec.kind == StressSourceKind::SidecarHDF5,
               "ApplyCsmStressSidecar: StressSpec.kind must be SidecarHDF5");
   MFEM_VERIFY(!spec.sidecar_path.empty(),
               "ApplyCsmStressSidecar: spec.sidecar_path is empty");

   // Reject BP2 / antiplane FaultGeometry up front (before any accessor
   // that could read uninitialised per-DOF state).  The 3-D / BP5 ctor
   // is the only one that populates per-DOF coords, basis, and zero-
   // normal counts; ComputeSAFSParams<StressField3D> later re-checks
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

   StressField3D field(spec.sidecar_path);
   geom.ComputeSAFSParams(field,
                          spec.pore_pressure.P_p_pa,
                          spec.pore_pressure.P_p_grad_pa_per_m,
                          spec.pore_pressure.min_sigma_n_pa);
   MFEM_VERIFY(geom.HasSAFSParams(),
               "ApplyCsmStressSidecar: ComputeSAFSParams did not "
               "set HasSAFSParams() = true");
}

}  // namespace

void ApplyCsmStressSidecar(const StressSpec&       spec,
                           FaultGeometry<ParMesh>& geom)
{
   apply_csm_impl<ParMesh>(spec, geom);
}

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

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
