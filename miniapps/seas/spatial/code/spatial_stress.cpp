// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_stress.cpp — Phase 3 + Phase 3b implementation.

#include "spatial_stress.hpp"

#include "../../fault/fault_geometry_safs.inl"
// The templated ComputeParams<StressSource> overload (Phase 3b) is
// defined in this companion inline file:
#include "../../fault/fault_geometry_safs_templated.inl"

#include <algorithm>
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

   StressField3D field(spec.sidecar_path);
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

// =====================================================================
//  Phase 1 (TPV26/27): Tpv2627DepthStressSource
// =====================================================================

Tpv2627DepthStressSource::Tpv2627DepthStressSource(
   real_t rho, real_t g, real_t water_density,
   real_t b11, real_t b33, real_t b13,
   real_t omega_top_m, real_t omega_bot_m)
   : rho_(rho), g_(g), wd_(water_density),
     b11_(b11), b33_(b33), b13_(b13),
     omega_top_(omega_top_m), omega_bot_(omega_bot_m)
{
   MFEM_VERIFY(rho_ > 0.0,
               "Tpv2627DepthStressSource: rho must be > 0; got " << rho_);
   MFEM_VERIFY(g_ > 0.0,
               "Tpv2627DepthStressSource: g must be > 0; got " << g_);
   // R-007: a negative water_density gives Pf < 0, which flips the sign of the
   // effective vertical stress (sigma22 + Pf) and hence of sigma13 — silently
   // producing a LEFT-lateral background.  0 is the sanctioned "no pore
   // pressure" setting (total-stress profile).
   MFEM_VERIFY(wd_ >= 0.0,
               "Tpv2627DepthStressSource: water_density must be >= 0 "
               "(0 disables pore pressure); got " << wd_);
   MFEM_VERIFY(omega_bot_ > omega_top_,
               "Tpv2627DepthStressSource: omega_bot_m ("
               << omega_bot_ << ") must be > omega_top_m ("
               << omega_top_ << ").");
   constexpr real_t inf = std::numeric_limits<real_t>::infinity();
   bbox_ = { -inf, inf, -inf, inf, -inf, inf };
}

mfem::DenseMatrix Tpv2627DepthStressSource::Evaluate(real_t /*x*/,
                                                     real_t /*y*/,
                                                     real_t z) const
{
   // depth = max(0, -z); z<0 is down in the code frame (see PLAN §2.4).
   const real_t depth = std::max(static_cast<real_t>(0.0), -z);

   const real_t Pf      = wd_ * g_ * depth;          // hydrostatic fluid pressure
   const real_t sigma22 = -rho_ * g_ * depth;        // vertical/lithostatic (tension +)

   // Omega taper (spec §1.2).
   real_t Omega;
   if (depth <= omega_top_)      { Omega = 1.0; }
   else if (depth >= omega_bot_) { Omega = 0.0; }
   else { Omega = (omega_bot_ - depth) / (omega_bot_ - omega_top_); }

   const real_t sp = sigma22 + Pf;   // effective vertical (spec "sigma22 + Pf")
   const real_t sigma11 = Omega * (b11_ * sp - Pf) + (1.0 - Omega) * sigma22; // strike
   const real_t sigma33 = Omega * (b33_ * sp - Pf) + (1.0 - Omega) * sigma22; // fault-normal
   const real_t sigma13 = Omega * (b13_ * sp);                                // on-fault shear

   // Frame remap (spec -> code): xx=sigma11, yy=sigma33, zz=sigma22,
   // xy=sigma13, yz=xz=0.  The spec tensor is tension-positive; the
   // StressSource3D concept requires compression-POSITIVE, so negate all
   // components (this also bakes in the right-lateral-positive tau_strike
   // sign — see the header note and the golden unit test).
   mfem::DenseMatrix sigma(3, 3);
   sigma = 0.0;
   sigma(0, 0) = -sigma11;                        // xx (strike)
   sigma(1, 1) = -sigma33;                        // yy (fault-normal)
   sigma(2, 2) = -sigma22;                        // zz (vertical)
   sigma(0, 1) = -sigma13;  sigma(1, 0) = -sigma13;  // xy (on-fault shear)
   // sigma(1,2)=sigma(2,1)=sigma(0,2)=sigma(2,0)=0 (no yz/xz shear).
   return sigma;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
