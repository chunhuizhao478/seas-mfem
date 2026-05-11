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

#include <atomic>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

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

   /// Test instrumentation: returns the number of times Project /
   /// ProjectSerial / ProjectVelocity has been called since the
   /// process started OR since the last ResetCallCount().  Used by
   /// the one-time-load contract test (T-4-7).
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

   /// Throw the "mesh not contained" abort with a formatted bbox table.
   static void AbortContainmentFailure(
      const DataField3D& field,
      real_t mxmin, real_t mxmax,
      real_t mymin, real_t mymax,
      real_t mzmin, real_t mzmax);

   /// Throw the "post-projection out of declared bounds" abort.
   static void AbortRangeFailure(
      const std::string& field_name,
      real_t observed_lo, real_t observed_hi,
      real_t declared_lo, real_t declared_hi);

   static std::atomic<int> call_count_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FIELD_COEFFICIENT_HPP
