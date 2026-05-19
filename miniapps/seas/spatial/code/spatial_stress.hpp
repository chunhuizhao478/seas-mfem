// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_stress.hpp — Phase 3 + Phase 3b of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// Consolidated header (rev-3 §"What is genuinely new") that provides:
//   - the StressSource3D *concept* documentation (no virtual base)
//   - ApplyCsmStressSidecar (Phase 3, CSM HDF5 path)
//   - ConstantTensorStressSource (Phase 3b, D-1 testing primitive)
//
// StressSource3D concept (duck-typed):
//
//   Any class S satisfying
//       mfem::DenseMatrix S::Evaluate(real_t x, real_t y, real_t z) const
//       const std::array<real_t, 6>& S::BBox() const
//       bool S::ContainsBBox(real_t xmin, real_t xmax,
//                            real_t ymin, real_t ymax,
//                            real_t zmin, real_t zmax,
//                            real_t eps = 0.0) const
//   satisfies StressSource3D.  Implementations:
//      - mfem::seas::StressField3D                (Phase 3, CSM HDF5)
//      - mfem::seas::spatial::ConstantTensorStressSource (Phase 3b, D-1)
//
//   The templated FaultGeometry::ComputeParams<StressSource>(...)
//   overload (declared in fault/fault_geometry.hpp, defined in
//   fault/fault_geometry_safs_templated.inl) consumes any S satisfying
//   the concept.  Overload resolution always prefers the existing
//   non-templated ComputeParams(const StressField3D&, ...) overload
//   for StressField3D arguments, preserving the BP5 path byte-exact.

#ifndef MFEM_SEAS_SPATIAL_STRESS_HPP
#define MFEM_SEAS_SPATIAL_STRESS_HPP

#include "mfem.hpp"

#include "../../fault/fault_geometry.hpp"
#include "../../dynamic/heterogeneous_material.hpp"   // Phase R.5: MaterialField
#include "../../io/stress_field_3d.hpp"

#include "spatial_friction.hpp"   // StressSpec, PorePressureSpec

#include <array>
#include <functional>
#include <limits>

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
//  Phase 3: CSM HDF5 sidecar wrapper
// =====================================================================

/// Loads StressField3D from the sidecar and calls
/// geom.ComputeParams(field, P_p_pa, ...).  After return,
/// geom.HasParams() == true; geom.sigma_n_per_dof() and
/// geom.GetTauPre() are populated.
///
/// Pre-flight check: aborts with the same precise message as the
/// RateStateFaultOperator::SetSAFSMode R-001 guard if
/// geom.NumZeroNormalFallbacks() != 0.
///
/// Aborts on sidecar load failure or projection failure (no silent
/// fallback).
void ApplyCsmStressSidecar(const StressSpec&       spec,
                           FaultGeometry<ParMesh>& geom);

/// Serial-mesh overload (BP2 fault path; not used by the SAFS dynamic
/// driver but exposed for symmetry / completeness).
void ApplyCsmStressSidecar(const StressSpec&     spec,
                           FaultGeometry<mfem::Mesh>& geom);

// =====================================================================
//  Phase 3b: ConstantTensorStressSource (D-1 testing primitive)
// =====================================================================

/// Constant Cauchy stress tensor in EAST-NORTH-UP frame, compression
/// POSITIVE (SEAS convention; matches
/// stress/code/hickman_and_zoback_*).  Satisfies the StressSource3D
/// concept.
///
/// The rotation onto the per-DOF fault normal happens inside
/// FaultGeometry::ComputeParams<StressSource>(...), NOT in this
/// class.  Evaluate(x, y, z) returns the same DenseMatrix for every
/// query — the tensor is constant in space.
///
/// BBox is ±infinity in every direction; ContainsBBox is always true.
class ConstantTensorStressSource
{
public:
   /// All six independent components of the symmetric Cauchy tensor.
   /// Symmetry (sigma_xy == sigma_yx, etc.) is enforced by construction.
   ConstantTensorStressSource(real_t sxx, real_t syy, real_t szz,
                              real_t sxy, real_t syz, real_t sxz);

   /// Returns the same constant tensor for every (x, y, z).  Never
   /// aborts.
   mfem::DenseMatrix Evaluate(real_t x, real_t y, real_t z) const;

   /// BBox is the canonical "infinite" box {-inf, +inf, ...}.
   const std::array<real_t, 6>& BBox() const { return bbox_; }

   /// Always true — a constant-tensor source is defined everywhere.
   bool ContainsBBox(real_t /*xmin*/, real_t /*xmax*/,
                     real_t /*ymin*/, real_t /*ymax*/,
                     real_t /*zmin*/, real_t /*zmax*/,
                     real_t /*eps*/ = 0.0) const { return true; }

   /// Direct access to the 3x3 tensor (test convenience).
   const mfem::DenseMatrix& Tensor() const { return sigma_; }

private:
   mfem::DenseMatrix sigma_;   // 3x3, populated in ctor.
   std::array<real_t, 6> bbox_;
};

// =====================================================================
//  ConstantTensorWithPatchesStressSource (TPV205-style square patches)
// =====================================================================

/// Background constant Cauchy tensor (EAST-NORTH-UP, compression
/// POSITIVE) with zero or more rectangular patches that override
/// individual stress components at t = 0.  Satisfies the
/// StressSource3D concept.
///
/// Per-DOF semantics: `Evaluate(x, y, z)` returns the background tensor;
/// then for each patch (in document order) whose `contains(x,y,z)` is
/// true, any non-NaN `sigma_*_pa` override replaces that component in
/// the returned tensor.  Last-match wins on a per-component basis,
/// matching the TPV205 spec where the three patches are spatially
/// disjoint.
///
/// SCEC TPV205 encoding (3 patches at depth = 7.5 km, half-side 1.5 km,
/// fault y = 0):
///   * Nucleation (0, 0, 7500),    sigma_xy_pa = 81.6e6  (yield-exceeding)
///   * Left       (-7500, 0, 7500), sigma_xy_pa = 78.0e6
///   * Right      (+7500, 0, 7500), sigma_xy_pa = 62.0e6
///
/// BBox is ±infinity; ContainsBBox is always true (the background tensor
/// is defined everywhere, patches only modify a measure-zero region).
class ConstantTensorWithPatchesStressSource
{
public:
   ConstantTensorWithPatchesStressSource(real_t bg_sxx, real_t bg_syy,
                                         real_t bg_szz, real_t bg_sxy,
                                         real_t bg_syz, real_t bg_sxz,
                                         std::vector<StressPatch> patches);

   mfem::DenseMatrix Evaluate(real_t x, real_t y, real_t z) const;

   const std::array<real_t, 6>& BBox() const { return bbox_; }
   bool ContainsBBox(real_t /*xmin*/, real_t /*xmax*/,
                     real_t /*ymin*/, real_t /*ymax*/,
                     real_t /*zmin*/, real_t /*zmax*/,
                     real_t /*eps*/ = 0.0) const { return true; }

   const std::vector<StressPatch>& Patches() const { return patches_; }
   const mfem::DenseMatrix& BackgroundTensor() const { return bg_sigma_; }

private:
   mfem::DenseMatrix       bg_sigma_;
   std::vector<StressPatch> patches_;
   std::array<real_t, 6>   bbox_;
};

// =====================================================================
//  Phase R.5: DepthProportionalToShearModulusStressSource (TPV31)
// =====================================================================

/// TPV31-style pre-stress: a constant Cauchy tensor (specified at a
/// reference shear modulus `mu_ref_pa`) scaled at each spatial point
/// by `mu(point) / mu_ref`, where `mu(point)` is evaluated from a
/// MaterialField at the point's location.
///
/// Per SCEC TPV31 spec p. 6: the stress is depth-dependent and
/// proportional to the local shear modulus.  Components are stored
/// in MPa in the TOML config and converted to Pa during parsing
/// (see DepthProportionalStressSpec).
///
/// Satisfies the StressSource3D concept: `Evaluate(x, y, z)` returns
/// the tensor at that point; `BBox` is ±infinity; `ContainsBBox` is
/// always true.
///
/// `mu_at_xyz` is a callback that returns the shear modulus at a
/// physical point (x, y, z).  The driver constructs this from a
/// `DepthProfile1DMaterial` wrapper or any other coordinate-aware
/// source — keeps this class decoupled from MaterialField's coupling
/// to MFEM ElementTransformation (which would otherwise force the
/// caller to manufacture a synthetic transformation per stress query).
class DepthProportionalToShearModulusStressSource
{
public:
   using MuAtFn = std::function<real_t(real_t, real_t, real_t)>;

   /// @param sxx_per_mu...sxz_per_mu  Six scaled components in Pa
   ///        (already multiplied by 1e6 from the TOML MPa input).
   /// @param mu_ref_pa                Reference shear modulus (Pa).
   /// @param mu_at_xyz                Callback: returns mu(x,y,z) in Pa.
   DepthProportionalToShearModulusStressSource(
      real_t sxx_per_mu, real_t syy_per_mu, real_t szz_per_mu,
      real_t sxy_per_mu, real_t syz_per_mu, real_t sxz_per_mu,
      real_t mu_ref_pa,
      MuAtFn mu_at_xyz);

   /// Returns the per-point Cauchy tensor (3x3) in EAST-NORTH-UP
   /// frame, compression POSITIVE.
   mfem::DenseMatrix Evaluate(real_t x, real_t y, real_t z) const;

   const std::array<real_t, 6>& BBox() const { return bbox_; }

   bool ContainsBBox(real_t /*xmin*/, real_t /*xmax*/,
                     real_t /*ymin*/, real_t /*ymax*/,
                     real_t /*zmin*/, real_t /*zmax*/,
                     real_t /*eps*/ = 0.0) const { return true; }

   real_t MuRefPa() const { return mu_ref_pa_; }

private:
   real_t sxx_per_mu_, syy_per_mu_, szz_per_mu_;
   real_t sxy_per_mu_, syz_per_mu_, sxz_per_mu_;
   real_t mu_ref_pa_;
   MuAtFn mu_at_xyz_;
   std::array<real_t, 6> bbox_;
};

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif // MFEM_SEAS_SPATIAL_STRESS_HPP
