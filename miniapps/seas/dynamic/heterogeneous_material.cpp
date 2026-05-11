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

} // namespace seas
} // namespace mfem
