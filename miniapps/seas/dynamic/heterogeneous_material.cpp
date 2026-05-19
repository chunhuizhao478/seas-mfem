// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// heterogeneous_material.cpp — out-of-line definitions for
// `MaterialField` members that cannot be header-inline:
//   * `MakeCoefficient` — trivial factory but kept out-of-line to
//     match `MaxCpInElement`'s linkage.
//   * `MaxCpInElement` — Coefficient mode walks the element's
//     quadrature points via `IntRules`, which is too heavy to
//     inline at every call site.

#include "heterogeneous_material.hpp"

#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

MaterialField MaterialField::MakeCoefficient(mfem::Coefficient* lambda,
                                             mfem::Coefficient* mu,
                                             mfem::Coefficient* rho)
{
   MFEM_VERIFY(lambda && mu && rho,
               "MakeCoefficient: all three Coefficient pointers must be "
               "non-null.");
   MaterialField m;
   m.mode        = Mode::Coefficient;
   m.lambda_coef = lambda;
   m.mu_coef     = mu;
   m.rho_coef    = rho;
   return m;
}


real_t MaterialField::MaxCpInElement(
   int elem, mfem::ElementTransformation* T) const
{
   if (mode == Mode::Constant)
   {
      return std::sqrt((lambda_const + 2.0 * mu_const) / rho_const);
   }
   if (mode == Mode::GridFunction)
   {
      MFEM_ASSERT(rho_gf && lambda_gf && mu_gf,
                  "MaterialField::MaxCpInElement: GridFunction mode "
                  "requires all three GFs to be set");
      mfem::ParFiniteElementSpace* fes = rho_gf->ParFESpace();
      MFEM_ASSERT(fes != nullptr,
                  "MaterialField::MaxCpInElement: rho_gf has null "
                  "ParFESpace");
      mfem::Array<int> vdofs;
      fes->GetElementDofs(elem, vdofs);
      // R-003 round-4: sentinel = -inf (matches the Coefficient path
      // below).  A 0.0 sentinel would silently produce dt = h / 0 = inf
      // downstream if all samples return NaN (corrupted sidecar with
      // rho == 0 etc.).
      real_t cp_max = -std::numeric_limits<real_t>::infinity();
      for (int d = 0; d < vdofs.Size(); ++d)
      {
         const int idx = vdofs[d];
         const real_t la  = (*lambda_gf)(idx);
         const real_t mu_ = (*mu_gf)(idx);
         const real_t rho = (*rho_gf)(idx);
         const real_t cp  = std::sqrt((la + 2.0 * mu_) / rho);
         if (cp > cp_max) { cp_max = cp; }
      }
      return cp_max;
   }
   // Mode::Coefficient — sample at the element's quadrature points.
   MFEM_VERIFY(T != nullptr,
               "MaterialField::MaxCpInElement(elem, T) called in "
               "Mode::Coefficient with T == nullptr; an "
               "ElementTransformation is required to evaluate Coefficient "
               "at qpoints.");
   MFEM_VERIFY(lambda_coef && mu_coef && rho_coef,
               "MaterialField::MaxCpInElement in Mode::Coefficient "
               "requires all three Coefficient* pointers to be non-null.");

   // R-005 round-2: use a fixed `default_qorder = 2` quadrature
   // rule.  c_p varies slowly in space; an order-4 rule (11 points
   // for tetrahedra, 8 for hexahedra) is more than enough to bound
   // the per-element maximum.  Future callers needing higher
   // accuracy can extend MaxCpInElement with an explicit qorder
   // argument.
   const int default_qorder = 2;
   const mfem::IntegrationRule& ir = mfem::IntRules.Get(
      T->GetGeometryType(), 2 * default_qorder);

   // R-002 round-4: save the caller's IntegrationPoint pointer so we
   // can restore it before returning.  Coefficient subclasses that
   // rely on T.GetIntPoint() (rather than the explicit ip argument)
   // would otherwise observe our last internal qpoint, silently
   // corrupting any subsequent T-aware evaluation by the caller.
   //
   // Note on GetIntPoint(): MFEM's GetIntPoint() returns *IntPoint,
   // which is UB if the caller never called SetIntPoint.  The
   // &-of-dereference is treated as identity by every mainstream
   // compiler (no actual load), giving us back the internal pointer
   // — including nullptr if it was never set.  The conditional
   // restore at the bottom of this function then no-ops in that
   // case.
   const mfem::IntegrationPoint* const saved_ip = &T->GetIntPoint();

   real_t cp_max = -std::numeric_limits<real_t>::infinity();
   for (int q = 0; q < ir.GetNPoints(); ++q)
   {
      const mfem::IntegrationPoint& ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      const real_t la  = lambda_coef->Eval(*T, ip);
      const real_t mu_ = mu_coef    ->Eval(*T, ip);
      const real_t rho = rho_coef   ->Eval(*T, ip);
      const real_t cp  = std::sqrt((la + 2.0 * mu_) / rho);
      if (cp > cp_max) { cp_max = cp; }
   }
   if (saved_ip) { T->SetIntPoint(saved_ip); }
   return cp_max;
}

// =========================================================================
// Phase R.3 step 6 — DepthProfile1DMaterial constructor.
// =========================================================================
std::unique_ptr<DepthProfile1DMaterial> MakeDepthProfile1DMaterial(
   const std::vector<DepthProfileLayer>& layers,
   char depth_axis)
{
   MFEM_VERIFY(!layers.empty(),
               "MakeDepthProfile1DMaterial: layers list must be non-empty.");
   MFEM_VERIFY(depth_axis == 'x' || depth_axis == 'y' || depth_axis == 'z',
               "MakeDepthProfile1DMaterial: depth_axis must be 'x', 'y', "
               "or 'z'; got '" << depth_axis << "'.");

   // Validate monotonic ordering, positive material constants, valid interp.
   for (std::size_t i = 0; i < layers.size(); ++i)
   {
      const auto& L = layers[i];
      MFEM_VERIFY(L.vp_ms > 0.0,
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has vp_ms = " << L.vp_ms << " (must be > 0).");
      MFEM_VERIFY(L.vs_ms > 0.0,
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has vs_ms = " << L.vs_ms << " (must be > 0).");
      // REVIEW R-005 (round-4): Lamé conversion is
      //   lambda = rho * (vp^2 - 2*vs^2)
      // which is non-negative iff vp >= sqrt(2)*vs (Poisson ratio >= 0).
      // Real earth materials always satisfy this; catching the common
      // vp/vs swap at parse time gives a clear diagnostic instead of a
      // far-removed `BimaterialFlux::BuildGodunovStateFaceLocal`
      // "lam+2*mu must be positive" abort.
      MFEM_VERIFY(L.vp_ms >= std::sqrt(2.0) * L.vs_ms,
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has vp_ms=" << L.vp_ms
                  << " < sqrt(2)*vs_ms=" << (std::sqrt(2.0) * L.vs_ms)
                  << " (Poisson ratio < 0; almost certainly a vp/vs "
                  "swap in the TOML config — real earth materials have "
                  "vp > sqrt(2)*vs).");
      MFEM_VERIFY(L.rho_kgm3 > 0.0,
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has rho_kgm3 = " << L.rho_kgm3 << " (must be > 0).");
      MFEM_VERIFY(L.depth_bot_m >= L.depth_top_m,
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has depth_bot_m (" << L.depth_bot_m
                  << ") < depth_top_m (" << L.depth_top_m << ").");
      MFEM_VERIFY(L.interp == "constant" || L.interp == "linear",
                  "MakeDepthProfile1DMaterial: layer " << i
                  << " has interp = '" << L.interp
                  << "' (must be 'constant' or 'linear').");
      if (i > 0)
      {
         MFEM_VERIFY(L.depth_top_m >= layers[i-1].depth_top_m,
                     "MakeDepthProfile1DMaterial: layers must be ordered "
                     "shallow-to-deep (layer " << i << " depth_top_m="
                     << L.depth_top_m << " < layer " << (i-1)
                     << " depth_top_m=" << layers[i-1].depth_top_m << ").");
         // REVIEW R-004 (round-4): require contiguity — depths in a
         // GAP between consecutive layers would silently take the
         // shallower layer's material (the lookup walks shallow-to-deep
         // and the chosen-index stays at the shallower layer for any
         // depth in the gap).  Tolerance 1e-9 lets users supply equal
         // top/bot values for a true discontinuity at a single point.
         MFEM_VERIFY(L.depth_top_m <= layers[i-1].depth_bot_m + 1e-9,
                     "MakeDepthProfile1DMaterial: layers must be "
                     "contiguous (layer " << (i-1) << " depth_bot_m="
                     << layers[i-1].depth_bot_m << " < layer " << i
                     << " depth_top_m=" << L.depth_top_m
                     << ").  Gap depths would silently take the shallower "
                     "layer's properties; use repeated top values for "
                     "true discontinuities instead.");
      }
   }

   const int axis_idx = depth_axis - 'x';   // 'x'=0, 'y'=1, 'z'=2

   auto wrapper = std::unique_ptr<DepthProfile1DMaterial>(
      new DepthProfile1DMaterial());

   // Copy the layers into a vector captured by the lambdas.  Use
   // shared_ptr<const vector> so all three coefficients see the same
   // table; the wrapper owns the coefficients (via unique_ptr) and they
   // in turn own a copy of the shared_ptr.
   auto layers_shared = std::make_shared<const std::vector<DepthProfileLayer>>(
      layers);

   // Per-coordinate evaluator that walks the layers and computes the
   // appropriate (lam, mu, rho) component via Lamé conversion:
   //   mu     = rho * vs^2
   //   lambda = rho * (vp^2 - 2 * vs^2)
   //
   // Component selector: 0 = lambda, 1 = mu, 2 = rho.  We can't capture
   // an enum by value through three lambdas + a templated functor without
   // adding ceremony, so each component lambda inlines the formula.
   auto eval_at_depth = [layers_shared](real_t depth,
                                        int component) -> real_t
   {
      // Tie-break R-007: at the jump y = y_jump, the DEEPER layer
      // wins.  Walk shallow-to-deep, skip if depth strictly less than
      // top; the FIRST layer that depth strictly exceeds top wins (and
      // because of the skip, the deeper layer is preferred on equality).
      // The reversed scan order is required: we want the LATEST layer
      // whose top <= depth.
      const auto& layers_v = *layers_shared;
      int chosen = -1;
      for (std::size_t i = 0; i < layers_v.size(); ++i)
      {
         // REVIEW R-007 enforcement: strict `<` so the boundary point
         // y = y_jump is NOT claimed by this layer (the next, deeper,
         // layer claims it instead).
         if (depth < layers_v[i].depth_top_m) { continue; }
         chosen = static_cast<int>(i);
      }
      if (chosen < 0)
      {
         // Above the shallowest layer's top — clamp to the shallowest
         // layer's values (this happens e.g. when y < 0 in TPV31's
         // half-space spec; the spec defines depth >= 0 so this clamp
         // is a defensive no-op for valid inputs).
         chosen = 0;
      }

      const auto& L = layers_v[chosen];
      real_t vp = L.vp_ms, vs = L.vs_ms, rho = L.rho_kgm3;

      if (L.interp == "linear" &&
          static_cast<std::size_t>(chosen) + 1 < layers_v.size())
      {
         const auto& Lnext = layers_v[chosen + 1];
         const real_t span = Lnext.depth_top_m - L.depth_top_m;
         if (span > 0.0)
         {
            real_t t = (depth - L.depth_top_m) / span;
            if (t < 0.0) { t = 0.0; }
            if (t > 1.0) { t = 1.0; }
            vp  = (1.0 - t) * L.vp_ms    + t * Lnext.vp_ms;
            vs  = (1.0 - t) * L.vs_ms    + t * Lnext.vs_ms;
            rho = (1.0 - t) * L.rho_kgm3 + t * Lnext.rho_kgm3;
         }
      }

      const real_t mu     = rho * vs * vs;
      const real_t lambda = rho * (vp * vp - 2.0 * vs * vs);

      switch (component)
      {
         case 0: return lambda;
         case 1: return mu;
         case 2: return rho;
         default:
            MFEM_ABORT("DepthProfile1DMaterial: invalid component "
                       << component << " (expected 0=lambda, 1=mu, 2=rho).");
      }
      return 0.0;  // unreachable
   };

   // Axis-to-depth conversion.  Depth is always POSITIVE in the layer
   // table (`depth_top_m` / `depth_bot_m` are spec-positive depths).
   // For the canonical SEAS mesh frame (CLAUDE.md "Canonical Coordinate
   // System": z = 0 at surface, z < 0 below) `axis = 'z'` requires
   // negating to recover positive depth.  For axes 'x' and 'y' the raw
   // value is used directly (TPV31 originally encoded depth on the
   // y-axis with positive y = depth; that legacy semantic stays
   // available so existing tests and pre-rotation configs keep working).
   auto axis_value_to_depth =
      [](real_t raw, int axis) -> real_t
   {
      if (axis == 2)
      {
         // Canonical: depth = max(0, -z).  Clamp so DOFs above the free
         // surface (z > 0, hypothetical) reuse the shallowest layer.
         const real_t d = -raw;
         return (d > 0.0) ? d : 0.0;
      }
      return raw;
   };

   // FunctionCoefficient signature is `real_t(*)(const Vector&)`.  We
   // pass an axis-aware lambda for each component.
   wrapper->lambda = std::make_unique<mfem::FunctionCoefficient>(
      [eval_at_depth, axis_idx, axis_value_to_depth]
      (const mfem::Vector& x) -> real_t
      {
         return eval_at_depth(axis_value_to_depth(x(axis_idx), axis_idx), 0);
      });
   wrapper->mu = std::make_unique<mfem::FunctionCoefficient>(
      [eval_at_depth, axis_idx, axis_value_to_depth]
      (const mfem::Vector& x) -> real_t
      {
         return eval_at_depth(axis_value_to_depth(x(axis_idx), axis_idx), 1);
      });
   wrapper->rho = std::make_unique<mfem::FunctionCoefficient>(
      [eval_at_depth, axis_idx, axis_value_to_depth]
      (const mfem::Vector& x) -> real_t
      {
         return eval_at_depth(axis_value_to_depth(x(axis_idx), axis_idx), 2);
      });

   // R-003: coordinate-only evaluator for callsites that need per-point
   // material without an ElementTransformation (e.g., the TPV31
   // DepthProportionalToShearModulusStressSource's mu_at_xyz callback).
   // Reuses the same eval_at_depth + axis_value_to_depth lambdas as the
   // three FunctionCoefficient instances above, so values are identical
   // by construction.
   wrapper->eval_at_xyz =
      [eval_at_depth, axis_idx, axis_value_to_depth]
      (real_t x, real_t y, real_t z,
       real_t& lambda_out, real_t& mu_out, real_t& rho_out) -> void
      {
         const real_t coords[3] = { x, y, z };
         const real_t depth =
            axis_value_to_depth(coords[axis_idx], axis_idx);
         lambda_out = eval_at_depth(depth, 0);
         mu_out     = eval_at_depth(depth, 1);
         rho_out    = eval_at_depth(depth, 2);
      };

   wrapper->field = MaterialField::MakeCoefficient(
      wrapper->lambda.get(), wrapper->mu.get(), wrapper->rho.get());

   return wrapper;
}

} // namespace seas
} // namespace mfem
