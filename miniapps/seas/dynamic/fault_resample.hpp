// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2 of the fault-dealiasing plan
// (document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md):
// the fault-face RESAMPLE operator R (mirrors SeisSol's init::resample).
//
// R is the L2 projector onto the degree-N polynomial space, expressed in
// per-QP value space (plan §4.2):
//
//     R = V (Vᵀ W V)⁻¹ Vᵀ W ,   V[q,i] = φ_i(x_q),   W = diag(w_q),
//
// where {φ_i} is any basis of the degree-N (= `order`) polynomial space on the
// reference fault face and {x_q, w_q} are the points/weights of the fault
// quadrature `ir` (the over-integrated rule from Phase 1).  R maps per-QP
// values → per-QP values:
//   - reproduces every mode 0..N EXACTLY (including the top mode N — it does
//     NOT annihilate genuine degree-N content), and
//   - discards the > N content that the over-integrated GPs expose but the
//     degree-N DOF basis cannot represent.
//
// Properties (verified in tests/unit/test_fault_resample.cpp):
//   * idempotent  R² = R,
//   * rank (N+1)(N+2)/2 on a triangle  (= #DOF),
//   * W-self-adjoint:  W R = (W R)ᵀ  (R is the W-orthogonal projector),
//   * R = I when #QP = #DOF (the minimal mass-matrix rule) ⇒ a no-op without
//     over-integration — over-integration is the prerequisite (plan §3.3/§4.2).
//
// R is basis-independent (any basis spanning the degree-N space gives the same
// R), and for affine faces the |J_F| scale cancels in the projector, so ONE
// reference R is shared by every fault face of a given (geometry, order,
// rule).  This header builds R; Phase 3 applies it to the per-step accumulated
// state increment (Δψ for rate-state; the slip-rate magnitude for LSW).

#ifndef MFEM_SEAS_FAULT_RESAMPLE_HPP
#define MFEM_SEAS_FAULT_RESAMPLE_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// Build the degree-`order` fault-face resample (L2 projection) matrix R for
/// the reference face geometry `geom` at the quadrature points of `ir`.
///
/// @param[in]  geom   reference face geometry (e.g. Geometry::TRIANGLE for tets)
/// @param[in]  order  DG polynomial order N (the degree the slip/state lives in)
/// @param[in]  ir     fault-face quadrature rule (the over-integrated rule;
///                    its #points must be >= the degree-N DOF count)
/// @param[out] R      `nq x nq` resample matrix (nq = ir.GetNPoints()).
///                    Apply per face as `out = R * in` over the face's QP block.
inline void BuildFaultResampleMatrix(Geometry::Type geom, int order,
                                     const IntegrationRule &ir, DenseMatrix &R)
{
   MFEM_VERIFY(order >= 1,
               "BuildFaultResampleMatrix: order must be >= 1, got " << order);

   // Degree-`order` full polynomial space P_order on the 2D reference face.
   // (L2 on a triangle is P_order with (order+1)(order+2)/2 dofs; the basis
   // choice is immaterial — R is basis-independent.)
   L2_FECollection fec(order, /*dim=*/2);
   const FiniteElement *fe = fec.FiniteElementForGeometry(geom);
   MFEM_VERIFY(fe != nullptr,
               "BuildFaultResampleMatrix: no FE for the requested geometry");

   const int ndof = fe->GetDof();      // (N+1)(N+2)/2 on a triangle
   const int nq   = ir.GetNPoints();
   MFEM_VERIFY(nq >= ndof,
               "BuildFaultResampleMatrix: rule has " << nq << " QPs < "
               << ndof << " DOFs — the degree-N space is not represented "
               "(rule too coarse for this order).");

   // V[q,i] = φ_i(x_q);  WV[q,i] = w_q φ_i(x_q).
   DenseMatrix V(nq, ndof), WV(nq, ndof);
   Vector shape(ndof);
   for (int q = 0; q < nq; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      fe->CalcShape(ip, shape);
      const real_t w = ip.weight;
      for (int i = 0; i < ndof; i++)
      {
         V(q, i)  = shape(i);
         WV(q, i) = w * shape(i);
      }
   }

   // M = Vᵀ W V  (ndof x ndof, SPD), then M ← M⁻¹.
   DenseMatrix M(ndof, ndof);
   MultAtB(V, WV, M);
   M.Invert();

   // R = V M⁻¹ (W V)ᵀ = V (Vᵀ W V)⁻¹ Vᵀ W   (nq x nq).
   DenseMatrix A(nq, ndof);
   Mult(V, M, A);          // A = V M⁻¹
   R.SetSize(nq, nq);
   MultABt(A, WV, R);      // R = A (WV)ᵀ
}

/// Per-face apply: `out[k] = Σ_j R(k,j) in[j]` for a single fault face's QP
/// block (length nq = R.Height()).  `in` and `out` may NOT alias.
inline void ApplyFaultResample(const DenseMatrix &R,
                               const real_t *in, real_t *out)
{
   const int nq = R.Height();
   for (int k = 0; k < nq; k++)
   {
      real_t s = 0.0;
      for (int j = 0; j < nq; j++) { s += R(k, j) * in[j]; }
      out[k] = s;
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_RESAMPLE_HPP
