// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// heterogeneous_material.hpp — `MaterialField` adapter for the data-
// projection feature.
//
// `MaterialField` is the runtime container that the WaveOperator
// (and consumer adapters) reads to look up `(λ, μ, ρ)` at every
// (element, dof) — or at a quadrature point in Coefficient mode.
// Three modes are supported:
//
//   - `Mode::Constant`     — homogeneous case, scalar values.  This
//                             is the back-compat case for every
//                             existing benchmark driver.
//   - `Mode::GridFunction` — heterogeneous case, three
//                             `ParGridFunction`s (one per parameter)
//                             that share a `ParFiniteElementSpace`
//                             (the same FES used by the wave field).
//   - `Mode::Coefficient`  — heterogeneous case, three
//                             `mfem::Coefficient*` pointers (typically
//                             the `LambdaFromSidecar` / `MuFromSidecar`
//                             / `RhoFromSidecar` classes defined in
//                             `io/material_coefficients.hpp`, but any
//                             `mfem::Coefficient` works).  Evaluated
//                             directly at each quadrature point with
//                             no FE-interpolation step — this is the
//                             recommended mode for new SAFS / CVM-H
//                             drivers because it preserves the
//                             sidecar's trilinear accuracy at every
//                             qpoint (whereas `Mode::GridFunction`
//                             smears it through the P1 basis).
//
// Accessor contracts:
//   * `At(elem, dof, ...)`   — Constant + GridFunction.  Aborts in
//                              Mode::Coefficient (which has no nodal
//                              representation; use `EvalAt` instead).
//   * `EvalAt(elem, T, ip,
//             ..., ...)`      — Hot-path qpoint accessor.  Required
//                              for Mode::Coefficient; available in
//                              Mode::Constant (returns the constant
//                              triple).  Aborts in Mode::GridFunction.
//   * `MaxCpInElement(elem,
//                     T*)`   — CFL helper.  In Mode::Coefficient
//                              walks the element's quadrature points
//                              via `T` (must be non-null).  In the
//                              other modes the `T*` argument is
//                              ignored.

#ifndef MFEM_SEAS_HETEROGENEOUS_MATERIAL_HPP
#define MFEM_SEAS_HETEROGENEOUS_MATERIAL_HPP

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>

namespace mfem
{
namespace seas
{

struct MaterialField
{
   enum class Mode : int
   {
      Constant     = 0,
      GridFunction = 1,
      Coefficient  = 2
   };

   Mode mode = Mode::Constant;

   // Scalar fallback (mode == Constant).
   real_t lambda_const = 0.0;
   real_t mu_const     = 0.0;
   real_t rho_const    = 0.0;

   // Heterogeneous (mode == GridFunction).  All three GFs share the
   // SAME `ParFiniteElementSpace` (the same one used by the wave
   // operator's state).  Per-element DOF index 0..(ndof_per_elem-1)
   // is the index used by `At`.
   std::shared_ptr<mfem::ParGridFunction> rho_gf;
   std::shared_ptr<mfem::ParGridFunction> lambda_gf;
   std::shared_ptr<mfem::ParGridFunction> mu_gf;

   // Heterogeneous (mode == Coefficient).  Non-owning pointers; the
   // caller (driver) owns the underlying objects and must keep them
   // alive at least until this MaterialField is destroyed.
   mfem::Coefficient* lambda_coef = nullptr;
   mfem::Coefficient* mu_coef     = nullptr;
   mfem::Coefficient* rho_coef    = nullptr;

   /// Convenience: build a constant-mode MaterialField.
   static MaterialField MakeConstant(real_t lambda, real_t mu, real_t rho)
   {
      MaterialField m;
      m.mode = Mode::Constant;
      m.lambda_const = lambda;
      m.mu_const     = mu;
      m.rho_const    = rho;
      return m;
   }

   /// Convenience: build a GridFunction-mode MaterialField.  Caller
   /// retains shared ownership; the struct stores shared_ptr copies.
   /// All three pointers MUST be non-null (R-004 round-2 tightening).
   static MaterialField MakeGridFunction(
      std::shared_ptr<mfem::ParGridFunction> rho_gf,
      std::shared_ptr<mfem::ParGridFunction> lambda_gf,
      std::shared_ptr<mfem::ParGridFunction> mu_gf)
   {
      MFEM_VERIFY(rho_gf && lambda_gf && mu_gf,
                  "MakeGridFunction: all three ParGridFunction "
                  "shared_ptrs must be non-null.");
      MaterialField m;
      m.mode      = Mode::GridFunction;
      m.rho_gf    = std::move(rho_gf);
      m.lambda_gf = std::move(lambda_gf);
      m.mu_gf     = std::move(mu_gf);
      return m;
   }

   /// Convenience: build a Coefficient-mode MaterialField from three
   /// non-owning `mfem::Coefficient` pointers.  All three MUST be
   /// non-null.  Defined in `heterogeneous_material.cpp`.
   static MaterialField MakeCoefficient(mfem::Coefficient* lambda,
                                        mfem::Coefficient* mu,
                                        mfem::Coefficient* rho);

   /// Element-and-DOF-local accessor used inside ADER assembly.
   /// In Constant mode this is branch-free and returns the scalar
   /// triple; in GridFunction mode it indexes the three GFs by
   /// element-local DOF.  In Mode::Coefficient it ABORTS (no nodal
   /// representation; callers must use `EvalAt` with a qpoint).
   inline void At(int elem, int dof,
                  real_t& lambda_out, real_t& mu_out, real_t& rho_out) const
   {
      if (mode == Mode::Constant)
      {
         lambda_out = lambda_const;
         mu_out     = mu_const;
         rho_out    = rho_const;
         return;
      }
      if (mode == Mode::Coefficient)
      {
         MFEM_ABORT("MaterialField::At called in Mode::Coefficient — "
                    "use EvalAt(elem, T, ip, ...) instead (Coefficient "
                    "mode has no nodal representation).");
      }
      MFEM_ASSERT(rho_gf && lambda_gf && mu_gf,
                  "MaterialField::At: GridFunction mode requires all "
                  "three GFs to be set");
      mfem::ParFiniteElementSpace* fes = rho_gf->ParFESpace();
      MFEM_ASSERT(fes != nullptr,
                  "MaterialField::At: rho_gf has null ParFESpace");
      mfem::Array<int> vdofs;
      fes->GetElementDofs(elem, vdofs);
      MFEM_ASSERT(dof >= 0 && dof < vdofs.Size(),
                  "MaterialField::At: dof out of range for element");
      const int idx = vdofs[dof];
      lambda_out = (*lambda_gf)(idx);
      mu_out     = (*mu_gf)(idx);
      rho_out    = (*rho_gf)(idx);
   }

   /// Quadrature-point accessor used by Mode::Coefficient consumers
   /// (the dynamic-rupture ADER hot path; the quasi-dynamic DG
   /// integrators access the Coefficient directly via the operator's
   /// view pointers).  `elem` is the local element index; `T` is the
   /// `ElementTransformation` for `elem`; `ip` is a reference-element
   /// `IntegrationPoint`.
   ///
   ///   - Mode::Constant     → returns the constant triple (T/ip unused).
   ///   - Mode::Coefficient  → calls Eval(T, ip) on each of the three
   ///                          non-null Coefficient*.
   ///   - Mode::GridFunction → MFEM_ABORT (callers must use At() in
   ///                          that mode; Coefficient is the only mode
   ///                          that meaningfully consumes a qpoint).
   inline void EvalAt(int /*elem*/,
                      mfem::ElementTransformation& T,
                      const mfem::IntegrationPoint& ip,
                      real_t& lambda_out, real_t& mu_out,
                      real_t& rho_out) const
   {
      if (mode == Mode::Constant)
      {
         lambda_out = lambda_const;
         mu_out     = mu_const;
         rho_out    = rho_const;
         return;
      }
      if (mode == Mode::Coefficient)
      {
         MFEM_ASSERT(lambda_coef && mu_coef && rho_coef,
                     "MaterialField::EvalAt: Coefficient mode requires "
                     "all three Coefficient* pointers to be non-null");
         lambda_out = lambda_coef->Eval(T, ip);
         mu_out     = mu_coef->Eval(T, ip);
         rho_out    = rho_coef->Eval(T, ip);
         return;
      }
      // Mode::GridFunction
      MFEM_ABORT("MaterialField::EvalAt called in Mode::GridFunction — "
                 "use At(elem, dof, ...) instead.");
   }

   /// Maximum c_p over the DOFs (or quadrature points, in Coefficient
   /// mode) in the given element.  `T` is required in Mode::Coefficient
   /// and ignored otherwise.  Definition lives in
   /// `heterogeneous_material.cpp`.
   real_t MaxCpInElement(int elem,
                         mfem::ElementTransformation* T = nullptr) const;
};

// =========================================================================
// Phase R.3 step 6 (PLAN_phase_R_exact_bimaterial_riemann_rev3.md):
// piecewise-1D depth-varying material profile (TPV31).
// =========================================================================

/// One layer of a piecewise-1D depth-varying material profile.
/// Layers MUST be ordered shallow-to-deep (smaller `depth_top_m` first);
/// contiguous (layer[i].depth_bot_m == layer[i+1].depth_top_m, possibly
/// with a repeat at a discontinuity).
///
/// `interp = "linear"` linearly interpolates (vp, vs, rho) between this
/// layer's TOP values and the NEXT layer's top values; "constant" holds
/// the layer's values across the whole [depth_top_m, depth_bot_m] range.
struct DepthProfileLayer
{
   real_t      depth_top_m  = 0.0;
   real_t      depth_bot_m  = 0.0;
   real_t      vp_ms        = 0.0;
   real_t      vs_ms        = 0.0;
   real_t      rho_kgm3     = 0.0;
   std::string interp;          // "constant" | "linear"
};

/// Wrapper that owns three `mfem::FunctionCoefficient` instances and
/// bundles them into a `MaterialField` (Mode::Coefficient).
///
/// Use the wrapper's `field` member at all callsites; do NOT extract
/// the raw coefficient pointers separately — their lifetimes are tied
/// to the wrapper's destruction.
///
/// Jump tie-break (REVIEW R-007): at exactly the boundary `y = y_jump`
/// between two layers, the DEEPER layer (larger depth_top_m) claims the
/// point.  Enforced in the lookup via `if (depth < layer.depth_top_m)
/// continue;`.  This rule MUST be matched by the analogous cohesion
/// `C_0(y)` evaluation so a single y_jump DOF receives consistent
/// material + friction values.
struct DepthProfile1DMaterial
{
   std::unique_ptr<mfem::FunctionCoefficient> lambda;
   std::unique_ptr<mfem::FunctionCoefficient> mu;
   std::unique_ptr<mfem::FunctionCoefficient> rho;
   MaterialField                              field;

   /// Coordinate-only evaluator: (x, y, z) → (lambda, mu, rho).  Reuses
   /// the same layer walk + axis-to-depth conversion that drives the
   /// three `FunctionCoefficient` instances above.  Provided for
   /// callsites that need per-point material without an
   /// `ElementTransformation` (e.g., the TPV31
   /// `DepthProportionalToShearModulusStressSource`'s `mu_at_xyz`
   /// callback).  Populated by `MakeDepthProfile1DMaterial`; safe to
   /// call from any thread.
   std::function<void(real_t x, real_t y, real_t z,
                      real_t& lambda_out,
                      real_t& mu_out,
                      real_t& rho_out)>     eval_at_xyz;

   // Non-copyable, non-movable (the MaterialField stores raw Coefficient*
   // pointers to the unique_ptr-owned members; moving would invalidate
   // them).
   DepthProfile1DMaterial() = default;
   DepthProfile1DMaterial(const DepthProfile1DMaterial&) = delete;
   DepthProfile1DMaterial& operator=(const DepthProfile1DMaterial&) = delete;
   DepthProfile1DMaterial(DepthProfile1DMaterial&&) = delete;
   DepthProfile1DMaterial& operator=(DepthProfile1DMaterial&&) = delete;
};

/// Build a piecewise-1D depth-varying material profile.
///
/// @param layers       Shallow-to-deep ordered list (non-empty).
/// @param depth_axis   One of 'x', 'y', 'z' — selects which physical
///                     coordinate is the depth coordinate.  TPV31 uses
///                     'y' (depth is the y-axis, increasing downward
///                     per SCEC TPV31 spec p. 3).
///
/// Aborts on invalid input: empty `layers`, non-monotonic ordering,
/// non-positive (vp, vs, rho), invalid `interp`, invalid `depth_axis`.
std::unique_ptr<DepthProfile1DMaterial> MakeDepthProfile1DMaterial(
   const std::vector<DepthProfileLayer>& layers,
   char depth_axis);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_HETEROGENEOUS_MATERIAL_HPP
