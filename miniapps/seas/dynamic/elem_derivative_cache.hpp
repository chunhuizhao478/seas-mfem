// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// elem_derivative_cache.hpp — Lever 1 of the ADER hot-path optimization
// (document/caliper_perfgraph_dev/ader_hotpath_optimization_plan_*.md).
//
// Builds, ONCE at setup, the per-element fused derivative operators
//
//     D_d^e = M_e^{-1} K_d^e ,   K_d^e[i,j] = Σ_q w_q φ_i(x_q) ∂_{x_d} φ_j(x_q)
//
// so that the element-local L2-projected spatial derivative
//   dQ_c^e = D_d^e · Q_c^e
// reproduces WaveOperator::ApplySpatialDerivative's on-the-fly quadrature kernel
// to round-off, with NO per-call CalcShape / CalcPhysDShape / Jacobian inversion.
//
// This is the ISOLATED new logic for the "flagged in-place" landing strategy —
// the operator stores the result and reads it from the `DerivMode::Cached`
// branch; the on-the-fly path is untouched.  Geometry-only: independent of
// NUM_STATE and of the (scalar/bimaterial) material dispatch, so a single cache
// serves both WaveOperator and BimaterialWaveOperator.
//
// REVIEW R-002: the cached apply is NOT bit-identical to the on-the-fly kernel
// (the setup `M^{-1}·K_d` matrix product re-associates the round-off); the
// equivalence gate is ≤ 1e-12 relative, not bit-equality.

#ifndef MFEM_SEAS_ELEM_DERIVATIVE_CACHE_HPP
#define MFEM_SEAS_ELEM_DERIVATIVE_CACHE_HPP

#include "mfem.hpp"

#include <array>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Build per-element D_d^e = M_e^{-1} K_d^e for d = 0,1,2.
///
/// Uses the SAME integration rule `IntRules.Get(geom, quad_order)` and the SAME
/// `CalcShape`/`CalcPhysDShape` that WaveOperator::ApplySpatialDerivative uses
/// on the fly (`quad_order == 2*order_`), so the cached matvec matches the
/// quadrature kernel up to floating-point re-association.
///
/// @param fes            Homogeneous-order L2 FE space (ParFiniteElementSpace
///                       passes as its FiniteElementSpace base).
/// @param quad_order     Integration-rule order; pass `2 * fe_order` to match
///                       the on-the-fly kernel exactly.
/// @param elem_mass_inv  Per-element inverse mass matrices M_e^{-1} (already
///                       assembled by WaveOperator::AssembleElementMassInverse).
/// @param[out] out       out[e][d] = the ndof×ndof matrix D_d^e.
inline void BuildElementDerivativeOperators(
   const FiniteElementSpace &fes,
   int quad_order,
   const std::vector<DenseMatrix> &elem_mass_inv,
   std::vector<std::array<DenseMatrix, 3>> &out)
{
   const int ne = fes.GetNE();
   MFEM_VERIFY(static_cast<int>(elem_mass_inv.size()) == ne,
               "BuildElementDerivativeOperators: elem_mass_inv size "
               << elem_mass_inv.size() << " != ne = " << ne);

   out.resize(ne);

   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const int ndof = fe->GetDof();

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), quad_order);
      const int nqp = ir.GetNPoints();

      // K_d^e[i,j] = Σ_q w_q φ_i(x_q) ∂_{x_d} φ_j(x_q), one matrix per direction.
      DenseMatrix Kd[3];
      for (int d = 0; d < 3; d++) { Kd[d].SetSize(ndof); Kd[d] = 0.0; }

      Vector shape(ndof);
      DenseMatrix dshape(ndof, 3);

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         const real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);
         fe->CalcPhysDShape(*Tr, dshape);

         for (int i = 0; i < ndof; i++)
         {
            const real_t w_shape_i = w * shape(i);
            for (int d = 0; d < 3; d++)
            {
               for (int j = 0; j < ndof; j++)
               {
                  Kd[d](i, j) += w_shape_i * dshape(j, d);
               }
            }
         }
      }

      // D_d^e = M_e^{-1} · K_d^e.
      const DenseMatrix &Minv = elem_mass_inv[e];
      MFEM_VERIFY(Minv.Height() == ndof && Minv.Width() == ndof,
                  "BuildElementDerivativeOperators: elem " << e
                  << " mass-inverse is " << Minv.Height() << "x"
                  << Minv.Width() << ", expected " << ndof << "x" << ndof);

      for (int d = 0; d < 3; d++)
      {
         out[e][d].SetSize(ndof);
         Mult(Minv, Kd[d], out[e][d]);   // out = Minv * Kd[d]
      }
   }
}

/// @brief Bytes one rank's D_d^e cache occupies: 3 · ne · ndof² · sizeof(real_t).
/// Used by the R-004 budget guard before building the cache.
inline std::size_t ElementDerivativeCacheBytes(int ne, int ndof_per_el)
{
   return static_cast<std::size_t>(3)
        * static_cast<std::size_t>(ne)
        * static_cast<std::size_t>(ndof_per_el)
        * static_cast<std::size_t>(ndof_per_el)
        * sizeof(real_t);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ELEM_DERIVATIVE_CACHE_HPP
