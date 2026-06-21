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
#include "../../io/stress_field_3d.hpp"
#include "../../config/bp5_params.hpp"   // BP5Params (Bp5AnalyticStressSource)

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
//  Phase 3c (QD): Bp5AnalyticStressSource (BP5-parity prestress)
// =====================================================================

/// Per-DOF stress source that reproduces BP5's analytic steady-state
/// prestress (config/bp5_params.hpp `tau0_vec` + effective `sigma_n`), so
/// the quasi-dynamic driver reaches the SAME `tau_pre_` as the analytic
/// `ComputeBP5Params` path through the generic
/// `FaultGeometry::ComputeParams<StressSource>` Cauchy projection.  A
/// uniform `FaultLocalPrestress` cannot match it because BP5's `tau0` is
/// `a(x2,x3)`-heterogeneous and carries a nucleation `delta_tau` (R-007).
/// Verification-only (BP5 / Phase 7 parity); NOT a SAF production prestress
/// (see PLAN_spatial_seas_quasidynamic_driver §Phase 3c, R-007 / OBS-1).
///
/// BP5's prestress is defined in the FAULT-LOCAL frame (dip, strike,
/// normal), whereas the `StressSource3D` concept hands `ComputeParams` a
/// Cauchy tensor that it projects onto the per-DOF fault basis.  For the
/// planar BP5 fault the basis `(n, t1=dip, t2=strike)` is CONSTANT, so we
/// take it from `FaultGeometry::fault_dof_basis()` (any column — assert
/// constness at the call site) and embed the fault-local prestress as
///   S = sigma_n (n⊗n) + tau_dip (t1⊗n + n⊗t1) + tau_strike (t2⊗n + n⊗t2),
/// which the projection inside `ComputeParams<StressSource>`
///   sigma_n = n·S·n,  tau_dip = t1·S·n,  tau_strike = t2·S·n
/// recovers EXACTLY on an orthonormal basis, reproducing
/// `ComputeBP5Params`' `tau_pre_` to round-off.  Taking the basis from the
/// geometry (rather than hard-coding sign conventions) keeps the embedded
/// prestress in the SAME frame as the elastic traction — the CLAUDE.md
/// sign-convention guard against a mixed-frame `tau_pre + traction`.
///
/// Coordinate mapping: `ComputeParams` calls `Evaluate` with the global
/// fault-DOF coords (`dof_coords_3d_`).  BP5's 2-D fault coords are
///   x2 = along-strike = X = x,   x3 = depth = -Z = -z
/// (the `GetFaultCoords2D` depth flip; elasticity_operator_traction.inl:105).
/// `Evaluate` is constant per fault DOF in space (depends only on x2,x3);
/// `BBox`/`ContainsBBox` are the canonical infinite box.
class Bp5AnalyticStressSource
{
public:
   /// @param bp5  BP5 spatial parameters (same struct `ComputeBP5Params` uses).
   /// @param n    fault normal   (`dof_basis_` rows 0-2; unit).
   /// @param t1   dip    tangent (`dof_basis_` rows 3-5; unit).
   /// @param t2   strike tangent (`dof_basis_` rows 6-8; unit).
   Bp5AnalyticStressSource(const BP5Params& bp5,
                           const real_t n[3], const real_t t1[3],
                           const real_t t2[3])
      : bp5_(bp5)
   {
      for (int d = 0; d < 3; ++d) { n_[d] = n[d]; t1_[d] = t1[d]; t2_[d] = t2[d]; }
      const real_t inf = std::numeric_limits<real_t>::infinity();
      bbox_ = { -inf, inf, -inf, inf, -inf, inf };
   }

   /// Embed the BP5 fault-local prestress at (x,y,z) into a symmetric
   /// Cauchy tensor (compression POSITIVE) that projects back to
   /// (sigma_n, tau_dip, tau_strike) on the stored fault basis.
   mfem::DenseMatrix Evaluate(real_t x, real_t /*y*/, real_t z) const
   {
      const real_t x2 = x;     // along-strike (X)
      const real_t x3 = -z;    // depth (= -Z; GetFaultCoords2D depth flip)
      real_t tau[2];
      bp5_.tau0_vec(x2, x3, tau);          // tau[0]=dip, tau[1]=strike
      const real_t tdip = tau[0], tstr = tau[1];
      const real_t sn = bp5_.sigma_n;      // effective sigma_n, compression +

      mfem::DenseMatrix S(3, 3);
      for (int r = 0; r < 3; ++r)
      {
         for (int c = 0; c < 3; ++c)
         {
            S(r, c) = sn   * n_[r]  * n_[c]
                    + tdip * (t1_[r] * n_[c] + n_[r] * t1_[c])
                    + tstr * (t2_[r] * n_[c] + n_[r] * t2_[c]);
         }
      }
      return S;
   }

   /// Defined everywhere — the canonical infinite box.
   const std::array<real_t, 6>& BBox() const { return bbox_; }

   /// Always true — a BP5-analytic source is defined on the whole fault.
   bool ContainsBBox(real_t /*xmin*/, real_t /*xmax*/,
                     real_t /*ymin*/, real_t /*ymax*/,
                     real_t /*zmin*/, real_t /*zmax*/,
                     real_t /*eps*/ = 0.0) const { return true; }

private:
   BP5Params bp5_;
   real_t n_[3], t1_[3], t2_[3];
   std::array<real_t, 6> bbox_;
};

// =====================================================================
//  Phase 10 (TPV31): DepthProportionalToShearModulusStressSource
// =====================================================================

/// TPV31-style pre-stress: a constant Cauchy tensor (specified at a
/// reference shear modulus `mu_ref_pa`) scaled at each spatial point by
/// `mu(point) / mu_ref`, where `mu(point)` is evaluated by the
/// `mu_at_xyz` callback (the driver builds it from a
/// `DepthProfile1DMaterial::eval_at_xyz`).  Per SCEC TPV31 spec p. 6.
///
/// Components are stored in Pa (already converted from the TOML MPa
/// input by the parser, see DepthProportionalStressSpec).  Satisfies the
/// StressSource3D concept: Evaluate / BBox (±inf) / ContainsBBox (true).
class DepthProportionalToShearModulusStressSource
{
public:
   using MuAtFn = std::function<real_t(real_t, real_t, real_t)>;

   /// @param sxx_per_mu...sxz_per_mu  Six scaled components in Pa.
   /// @param mu_ref_pa                Reference shear modulus (Pa).
   /// @param mu_at_xyz                Callback: returns mu(x,y,z) in Pa.
   DepthProportionalToShearModulusStressSource(
      real_t sxx_per_mu, real_t syy_per_mu, real_t szz_per_mu,
      real_t sxy_per_mu, real_t syz_per_mu, real_t sxz_per_mu,
      real_t mu_ref_pa,
      MuAtFn mu_at_xyz);

   /// Per-point Cauchy tensor (3x3), EAST-NORTH-UP, compression POSITIVE.
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
