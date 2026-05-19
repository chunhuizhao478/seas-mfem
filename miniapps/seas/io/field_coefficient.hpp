// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// field_coefficient.hpp — bind a DataField3D to MFEM via the
// `Coefficient` interface, plus a `FieldProjector` driver that
// projects sidecar fields onto a `(Par)FiniteElementSpace` with the
// schema-v1 invariants (interpolation-only mesh-bbox guard,
// post-projection sanity bounds, one-time-load call counter).
//
// See miniapps/seas/document/features_dev/data_projection_feature_plan_v2.md
// (Phase 4) for the contract.

#ifndef MFEM_SEAS_FIELD_COEFFICIENT_HPP
#define MFEM_SEAS_FIELD_COEFFICIENT_HPP

#include "mfem.hpp"

#include "data_field_3d.hpp"
#include "stress_field_3d.hpp"

#include <atomic>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

// Forward declarations — keep field_coefficient.hpp's dep graph free of
// fault/ to avoid an include cycle (fault_geometry.hpp can include this
// header transitively via the SAFS-mode plumbing).
template <typename MeshType> class FaultGeometry;

/// Scalar `mfem::Coefficient` backed by a `DataField3D`.
///
/// Example:
///   DataField3D vp("/path/to/sidecar.h5", "Vp");
///   FieldCoefficient vp_coef(vp);
///   gf.ProjectCoefficient(vp_coef);
///
/// `scale` and `offset` apply at evaluation: returned value =
/// scale * field(x, y, z) + offset.  No silent flooring is provided;
/// schema-v1 forbids any clamp policy at the runtime layer.
class FieldCoefficient : public mfem::Coefficient
{
public:
   FieldCoefficient(const DataField3D& field,
                    real_t scale = 1.0,
                    real_t offset = 0.0)
      : field_(field)
      , scale_(scale)
      , offset_(offset)
   { }

   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      mfem::Vector x(3);
      T.Transform(ip, x);
      return scale_ * field_.Evaluate(x[0], x[1], x[2]) + offset_;
   }

   const DataField3D& Field() const { return field_; }
   real_t Scale()  const { return scale_; }
   real_t Offset() const { return offset_; }

private:
   const DataField3D& field_;
   real_t scale_;
   real_t offset_;
};


/// Driver for one-shot field projection onto MFEM finite-element
/// spaces.  Performs:
///   1. mesh-vs-data bbox containment check (Phase 4 §3 in plan v2)
///   2. ProjectCoefficient
///   3. post-projection min/max bounds check against the field's
///      declared `[min_value, max_value]`.
/// Aborts loudly on violation of any of (1) or (3).
class FieldProjector
{
public:
   struct Result
   {
      /// Base-class pointer that holds either a `mfem::GridFunction`
      /// (serial path) or a `mfem::ParGridFunction` (parallel path).
      /// Callers on the parallel path may `dynamic_cast` to recover
      /// `ParGridFunction*`.
      std::shared_ptr<mfem::GridFunction> gf;
      real_t min_value;   ///< observed in projected GF
      real_t max_value;
   };

   /// Project a single named field. Calls `Field`'s `ContainsBBox`
   /// upfront against the mesh's global bbox; aborts if the mesh
   /// is not fully contained.
   ///
   /// `scale` and `offset` are forwarded to `FieldCoefficient` so
   /// callers can compose unit conversions / additive shifts at
   /// projection time.
   ///
   /// Increments the static call counter (`CallCount`).
   static Result Project(const DataField3D& field,
                         mfem::ParFiniteElementSpace& target_fes,
                         real_t scale = 1.0,
                         real_t offset = 0.0);

   /// Like Project, but for a serial `mfem::FiniteElementSpace`. Used
   /// in the unit tests to keep them MPI-free.
   static Result ProjectSerial(const DataField3D& field,
                               mfem::FiniteElementSpace& target_fes,
                               real_t scale = 1.0,
                               real_t offset = 0.0);

   /// Convenience for the velocity field set: returns three GFs in
   /// the order {ρ, λ, μ} given Vp, Vs, ρ in the sidecar.
   ///   μ = ρ V_s²,    λ = ρ V_p² − 2 μ.
   /// Per-field source bounds enforced at sidecar load + pre-flight;
   /// derived (λ, μ) bounds enforced at the supplied limits.
   struct VelocityFields
   {
      std::shared_ptr<mfem::ParGridFunction> rho;
      std::shared_ptr<mfem::ParGridFunction> lambda;
      std::shared_ptr<mfem::ParGridFunction> mu;
      std::shared_ptr<mfem::ParGridFunction> vp;
      std::shared_ptr<mfem::ParGridFunction> vs;
      real_t min_vp, max_vp, min_vs, max_vs;
      real_t min_rho, max_rho;
      real_t min_lambda, max_lambda, min_mu, max_mu;
   };

   /// Defaults track the relaxed source-field bounds set by the
   /// CVM-H writer (Vs_min=100, Vp_min=1000, rho_min=1000; see
   /// build_velocity_cvmh.DEFAULT_*_MIN_*).  Worst-case derived
   /// mu = rho_min * Vs_min^2 = 1.0e7 Pa.  lambda can be ~0 in the
   /// limit Vp -> sqrt(2)*Vs; the lower bound is therefore 0.
   /// Pass tighter values explicitly when projecting a sidecar that
   /// was built with stricter source bounds.
   ///
   /// ``interp`` selects the in-sidecar interpolation scheme used at
   /// every projector evaluation point (mesh DOF coord).  Default is
   /// ``InterpMode::Trilinear`` (matches the pre-existing schema-v1
   /// behaviour exactly).  ``InterpMode::CatmullRom`` uses a
   /// tensor-product 64-voxel cardinal cubic — see
   /// ``data_field_3d.hpp`` for the contract.
   static VelocityFields ProjectVelocity(
      const std::string&           sidecar_path,
      mfem::ParFiniteElementSpace& target_fes,
      real_t lambda_min_pa = 0.0,
      real_t lambda_max_pa = 1.0e12,
      real_t mu_min_pa     = 1.0e7,
      real_t mu_max_pa     = 1.0e11,
      InterpMode interp    = InterpMode::Trilinear);

   /// Six-component bulk Cauchy stress projection (Phase 6 §3 of
   /// PLAN_onfaultstress.md).
   ///
   /// Loads the six schema-v1 sigma_* fields from `sidecar_path`
   /// via a single `StressField3D` ctor (one HDF5 open per
   /// component), then projects each onto `target_fes` independently
   /// via the existing per-field `Project` path.  Six independent
   /// `Project` calls means the static `CallCount` increments by
   /// six per `ProjectStress` invocation.
   ///
   /// Sign convention: pure pass-through, compression POSITIVE
   /// (SEAS internal).  The bulk-path source-site flip lives in
   /// Phase 3's `bulk_stress_tensor_field`; this projector does
   /// no sign manipulation (R-501 / R-502).
   struct StressFields
   {
      // Components in schema-v1 canonical order: xx, yy, zz, xy, yz, xz.
      std::shared_ptr<mfem::ParGridFunction> sigma_xx;
      std::shared_ptr<mfem::ParGridFunction> sigma_yy;
      std::shared_ptr<mfem::ParGridFunction> sigma_zz;
      std::shared_ptr<mfem::ParGridFunction> sigma_xy;
      std::shared_ptr<mfem::ParGridFunction> sigma_yz;
      std::shared_ptr<mfem::ParGridFunction> sigma_xz;
      real_t min_sigma_xx, max_sigma_xx;
      real_t min_sigma_yy, max_sigma_yy;
      real_t min_sigma_zz, max_sigma_zz;
      real_t min_sigma_xy, max_sigma_xy;
      real_t min_sigma_yz, max_sigma_yz;
      real_t min_sigma_xz, max_sigma_xz;
   };

   /// Project all six bulk-stress components onto `target_fes`.
   /// See `StressFields` docstring for the sign-convention contract.
   ///
   /// `interp` selects the in-sidecar interpolation scheme used by
   /// each component reader, propagated to all six in lock-step via
   /// `StressField3D::SetInterpMode`.
   ///
   /// `scale` and `offset` are forwarded uniformly to all six inner
   /// `Project(...)` calls (mirrors the sibling `Project` /
   /// `ProjectVelocity` API).  Default values (1.0, 0.0) preserve
   /// the raw pass-through behaviour; non-default values are
   /// applied uniformly to every component (R-904).
   static StressFields ProjectStress(
      const std::string&           sidecar_path,
      mfem::ParFiniteElementSpace& target_fes,
      InterpMode interp            = InterpMode::Trilinear,
      real_t scale                 = 1.0,
      real_t offset                = 0.0);

   /// Rotate a sidecar bulk stress tensor onto the per-DOF fault basis
   /// to produce SEAS-internal pre-stress slots (Phase 6 §4 of
   /// PLAN_onfaultstress.md).
   ///
   /// Inputs from the sidecar are already in SEAS internal sign
   /// convention (compression POSITIVE, Pa).  The rotation projects
   /// the Cauchy tensor onto the fault-local frame:
   ///
   ///   sigma_n_per_dof(i)     = n_i^T σ(x_i) n_i
   ///                            − (P_p_pa + P_p_grad_pa_per_m
   ///                                          * max(0, -z_i))
   ///   tau_pre_per_dof(2*i)   = − t1_i^T σ(x_i) n_i   (dip,    BP5 t1)
   ///   tau_pre_per_dof(2*i+1) = − t2_i^T σ(x_i) n_i   (strike, BP5 t2)
   ///
   /// where (t1, t2, n) is the canonical SEAS fault-local frame
   /// (CLAUDE.md "fault-local tangent frame" rule: t1 = dip, t2 = strike).
   /// The leading minus sign on tau{1,2} matches the native
   /// TPV102/104/205 driver convention `DOFData::tau{1,2}_0 = +tau_ini`
   /// for σ_xy > 0 right-lateral driving stress (the raw T = σ·n
   /// projection returns the +y-side traction on the −y side, which is
   /// the Newton's-3rd-law mirror of the driving stress).  See R-001
   /// in debug_document/general_driver_debug_document/tpv102_tpv104_review.md
   /// and the implementation comment in field_coefficient.cpp.
   /// `sigma_n_total = n·S·n` is sign-invariant under n → −n.
   ///
   /// Pore pressure is subtracted from the normal stress: this
   /// implements the standard effective normal stress σ_n − P_p with
   /// the depth-dependent gradient `P_p_grad_pa_per_m * max(0, -z_i)`
   /// (clamped at the free surface, z_i > 0).
   ///
   /// Negative or below-threshold effective normal stresses are NOT
   /// clamped by default (`min_sigma_n_pa = 0.0`).  Callers can opt
   /// into a Pa-valued floor that warns once per rank and clamps
   /// sigma_n_per_dof(i) ← max(sigma_n_per_dof(i), min_sigma_n_pa)
   /// only when the value is below the floor (plan §1944-1951).
   ///
   /// @param field            Loaded six-component sidecar (Phase 6 §1).
   /// @param fault_geom       Source of the per-DOF coords and basis
   ///                         (Phase 6.A, populated only on the
   ///                         3-D / BP5 ctor path).
   /// @param[out] sigma_n_per_dof   Sized to NumFaultDOFs.
   /// @param[out] tau_pre_per_dof   Sized to 2 * NumFaultDOFs.
   /// @param P_p_pa                  Constant pore pressure term [Pa].
   /// @param P_p_grad_pa_per_m       Depth gradient of pore pressure
   ///                                 [Pa / m]; effective P_p at z is
   ///                                 P_p_pa + grad * max(0, -z).
   /// @param min_sigma_n_pa          Optional floor [Pa]; default 0 → no clamp.
   template <typename MeshType>
   static void ProjectFaultPreStress(
      const StressField3D&        field,
      const FaultGeometry<MeshType>& fault_geom,
      mfem::Vector&               sigma_n_per_dof,
      mfem::Vector&               tau_pre_per_dof,
      real_t                      P_p_pa = 0.0,
      real_t                      P_p_grad_pa_per_m = 0.0,
      real_t                      min_sigma_n_pa = 0.0);

   /// Test instrumentation: returns the number of times Project /
   /// ProjectSerial / ProjectVelocity / ProjectStress has been
   /// called since the process started OR since the last
   /// ResetCallCount().  Used by the one-time-load contract test
   /// (T-4-7); each ProjectStress call advances CallCount by 6
   /// (one per component).
   static int  CallCount()      { return call_count_.load(); }
   static void ResetCallCount() { call_count_.store(0); }

private:
   /// Compute mesh bbox over a serial mesh.
   static void ComputeMeshBBoxSerial(const mfem::Mesh& mesh,
                                     real_t& xmin, real_t& xmax,
                                     real_t& ymin, real_t& ymax,
                                     real_t& zmin, real_t& zmax);

   /// Compute mesh bbox over a parallel mesh (MPI_Allreduce).
   static void ComputeMeshBBoxParallel(const mfem::ParMesh& mesh,
                                       real_t& xmin, real_t& xmax,
                                       real_t& ymin, real_t& ymax,
                                       real_t& zmin, real_t& zmax);

   /// Compute (min, max) of a `GridFunction` over a serial mesh.
   static void ComputeGridFunctionMinMaxSerial(
      const mfem::GridFunction& gf, real_t& lo, real_t& hi);

   /// Compute (min, max) of a `ParGridFunction` (MPI_Allreduce).
   static void ComputeGridFunctionMinMaxParallel(
      const mfem::ParGridFunction& gf, real_t& lo, real_t& hi);

   /// Throw the "post-projection out of declared bounds" abort.
   static void AbortRangeFailure(
      const std::string& field_name,
      real_t observed_lo, real_t observed_hi,
      real_t declared_lo, real_t declared_hi);

   static std::atomic<int> call_count_;

public:
   /// Throw the "mesh not contained" abort with a formatted bbox table.
   ///
   /// Exposed publicly (was private) so external consumers can reuse the
   /// same containment-failure abort message instead of duplicating its
   /// formatting.  Consumers: spatial/code/spatial_velocity.cpp (Spatial
   /// Phase 2 of spatial_dynamic_rupture_plan.md rev-3 — sidecar
   /// containment pre-flight) and spatial/code/spatial_stress.cpp
   /// (Spatial Phase 3 — CSM stress sidecar wrapper).
   static void AbortContainmentFailure(
      const DataField3D& field,
      real_t mxmin, real_t mxmax,
      real_t mymin, real_t mymax,
      real_t mzmin, real_t mzmax);
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FIELD_COEFFICIENT_HPP
