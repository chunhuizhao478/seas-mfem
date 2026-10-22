// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dg_poro_div_integrator.hpp — poroelastic option, step P2.1 (plan rev 3).
//
// DGPoroDivIntegrator assembles the divergence form B_c of derivation
// rev 1.2, eq. (6.3), for one coefficient field c, with trial field a and
// test field b in the vector L2 space of u:
//
//   B_c(a, b) = sum_K int_K c (div a)(div b)
//             - sum_F int_F ( {c div a} [b.n] + {c div b} [a.n] )
//             + sum_F int_F sigma_w [a.n] [b.n]        (fluid block only)
//
// It is used three times (derivation section 6.4): c = alpha M for Z_wu
// (trial u, test q; Z_uw = Z_wu^T) with no penalty, and c = M with the
// penalty sigma_w of eq. (6.5) for Z_ww.  The face sums run over interior
// faces, which include the permeable fault faces: the slip d is tangential,
// so B_c carries no slip load (derivation section 7.1).
//
// The face block follows the pseudocode of derivation section 2.3, with side
// signs s_1 = +1 (MFEM's element 1, whose outward normal is n) and s_2 = -1:
//
//   Ntest_i  = s_i b_i.n          Dtest_i  = c_i div b_i
//   Ntrial_j = s_j a_j.n          Dtrial_j = c_j div a_j
//   block(i,j) += wq |J_F| [ -1/2 Ntest_i Dtrial_j      consistency
//                            -1/2 Dtest_i Ntrial_j      symmetry
//                            + sw Ntest_i Ntrial_j ]    penalty
//
// which is the (c0, c1, c2) = (-s_i/2, -s_j/2, s_i s_j sw) pattern of the
// elastic AssembleFaceBlock (integrator/dg_elasticity_ip_combined_integrator.hpp)
// with the signs folded into N.  Each entry is formed from products that are
// the same numbers in the transposed entry, so an interior face matrix (and a
// boundary one with both normal terms) is symmetric bit for bit, without the
// elastic integrator's averaging with the transpose.
//
// Boundary faces (derivation section 6.6) keep the one-sided terms the caller
// names at construction (BoundaryTerms below); MFEM's AddBdrFaceIntegrator
// with a marker picks the faces.  Free-surface and absorbing faces get no
// integrator.
//
// AssembleSourceRHS assembles the source loads of eq. (6.6), which have the
// pattern of the consistency half of B_c with a known scalar s in place of
// div a:
//
//   F(b) = sum_K int_K c s div b - sum_F int_F {c s} [b.n]
//
// F_V takes s = V with c = alpha M on the solid rows and c = M on the fluid
// rows; F_init takes s = p_init with c = 1 on the fluid rows.  On a boundary
// face it is the one-sided -int_F c s (b.n), and the caller passes only the
// faces whose load section 6.6 lists: Dirichlet-u faces for the solid rows of
// F_V, no-flow faces for the fluid rows of F_V and F_init.
// DGPoroDivSourceIntegrator wraps it for ParLinearForm.
//
// AssembleDirichletRHS (added in step P2.2) assembles the load of the
// kTrialNormal rows on a Dirichlet-u face, -int_F c (div b)(g.n): the
// kTrialNormal term with the prescribed displacement g in place of the trial
// field.  It is the g part of the coupling term -int alpha M (div q)
// ((u - g).n) of section 6.6, so with u = g on the face the two cancel and the
// face stays consistent once w != 0.  DGPoroDivDirichletIntegrator wraps it
// for ParLinearForm.
//
// Quadrature: faces use 2p + 1 points, as the elastic face integrator; the
// volume uses 2 OrderGrad, as MFEM's ElasticityIntegrator.  Both are exact
// for a coefficient that is constant on each element (PoroMaterialSample) on
// affine tetrahedra.

#ifndef MFEM_SEAS_DG_PORO_DIV_INTEGRATOR_HPP
#define MFEM_SEAS_DG_PORO_DIV_INTEGRATOR_HPP

#include "mfem.hpp"

#include <algorithm>

namespace mfem
{
namespace seas
{

class DGPoroDivIntegrator : public BilinearFormIntegrator
{
public:
   /// One-sided terms kept on a boundary face (derivation section 6.6), as
   /// bit flags.  An interior face always has all of its terms.
   enum BoundaryTerms : int
   {
      kInteriorOnly = 0, ///< a boundary face aborts
      kTestNormal   = 1, ///< -int_F c (div a)(b.n): no-flow rows of Z_wu
      kTrialNormal  = 2, ///< -int_F c (div b)(a.n): Dirichlet-u rows of Z_wu
      kBothNormal   = 3  ///< both: no-flow rows of Z_ww, with the penalty
   };

   /// @param c          the coefficient field (alpha M or M); kept, not owned.
   /// @param kappa_w    the fluid penalty factor of eq. (6.5); 0 = no penalty
   ///                   (the coupling blocks).
   /// @param bdr_terms  the boundary terms, for an integrator added with
   ///                   AddBdrFaceIntegrator.
   DGPoroDivIntegrator(Coefficient &c, real_t kappa_w = 0.0,
                       int bdr_terms = kInteriorOnly)
      : c_(&c), kappa_w_(kappa_w), bdr_terms_(bdr_terms)
   {
      MFEM_VERIFY(kappa_w >= 0.0,
                  "DGPoroDivIntegrator: kappa_w (" << kappa_w << ") must be >= 0.");
      MFEM_VERIFY(bdr_terms >= kInteriorOnly && bdr_terms <= kBothNormal,
                  "DGPoroDivIntegrator: unknown boundary terms " << bdr_terms);
      MFEM_VERIFY(kappa_w == 0.0 || bdr_terms == kInteriorOnly
                  || bdr_terms == kBothNormal,
                  "DGPoroDivIntegrator: the fluid penalty on a boundary face "
                  "belongs to the no-flow rows, which keep both normal terms "
                  "(derivation section 6.6).");
   }

   using BilinearFormIntegrator::AssembleElementMatrix;
   using BilinearFormIntegrator::AssembleFaceMatrix;

   /// Volume part of B_c on one element: c B_e^T W B_e (derivation
   /// section 6.4), in MFEM's vector layout (row p*ndof + k is component p
   /// of shape function k).
   void AssembleElementMatrix(const FiniteElement &el,
                              ElementTransformation &Tr,
                              DenseMatrix &elmat) override
   {
      const int dim = el.GetDim();
      const int nd = el.GetDof();
      const int n = dim * nd;
      elmat.SetSize(n);
      elmat = 0.0;

      const IntegrationRule &ir =
         IntRules.Get(el.GetGeomType(), 2 * Tr.OrderGrad(&el));

      DenseMatrix dshape(nd, dim), dshape_adj(nd, dim), adjJ(dim);
      Vector D(n);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcDShape(ip, dshape);
         CalcAdjugate(Tr.Jacobian(), adjJ);
         Mult(dshape, adjJ, dshape_adj);        // grad phi * detJ

         // div(phi_k e_p) = d phi_k / d x_p = dshape_adj(k, p) / detJ, and
         // dx = detJ dxi, so the weight is c / detJ.
         const real_t W = ip.weight * c_->Eval(Tr, ip) / Tr.Weight();
         for (int p = 0; p < dim; p++)
         {
            for (int k = 0; k < nd; k++) { D(p * nd + k) = dshape_adj(k, p); }
         }
         for (int s = 0; s < n; s++)
         {
            for (int r = 0; r < n; r++) { elmat(r, s) += W * (D(r) * D(s)); }
         }
      }
   }

   /// Face part of B_c (derivation section 2.3 pseudocode), plus the fluid
   /// penalty when kappa_w > 0.  Interior face: rows and columns are side 1's
   /// vector DOFs, then side 2's.  Boundary face: the one-sided terms of
   /// bdr_terms (c0 = -1 for kTestNormal, c1 = -1 for kTrialNormal) and the
   /// one-sided penalty.
   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override
   {
      const int dim = el1.GetDim();
      const int nd1 = el1.GetDof();
      const int nd2 = (Trans.Elem2No >= 0) ? el2.GetDof() : 0;
      const bool is_interior = (nd2 > 0);
      MFEM_VERIFY(is_interior || bdr_terms_ != kInteriorOnly,
                  "DGPoroDivIntegrator: a boundary face reached an integrator "
                  "built without boundary terms.");
      const int n1 = dim * nd1;
      const int n = dim * (nd1 + nd2);
      elmat.SetSize(n);
      elmat = 0.0;

      const int p = std::max(el1.GetOrder(), is_interior ? el2.GetOrder() : 0);
      const IntegrationRule &ir = IntRules.Get(Trans.GetGeometryType(), 2 * p + 1);

      // Consistency and symmetry coefficients (the sign s_i is in N).
      const real_t t0 = is_interior ? -0.5 : ((bdr_terms_ & kTestNormal) ? -1.0 : 0.0);
      const real_t t1 = is_interior ? -0.5 : ((bdr_terms_ & kTrialNormal) ? -1.0 : 0.0);

      Vector N(n), D(n), nor(dim), nhat(dim);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);
         CalcOrtho(Trans.Jacobian(), nor);      // |nor| = |J_F|, out of side 1
         const real_t nl_q = nor.Norml2();
         nhat.Set(1.0 / nl_q, nor);

         real_t c1, detJ1, c2 = 0.0, detJ2 = 0.0;
         SideRows_(el1, *Trans.Elem1, Trans.GetElement1IntPoint(), nhat, 1.0,
                   N, D, 0, c1, detJ1);
         if (is_interior)
         {
            SideRows_(el2, *Trans.Elem2, Trans.GetElement2IntPoint(), nhat, -1.0,
                      N, D, n1, c2, detJ2);
         }
         const real_t sw = (kappa_w_ > 0.0)
                           ? Penalty(dim, p, nl_q, detJ1, c1, detJ2, c2, is_interior)
                           : 0.0;

         const real_t W = ip.weight * nl_q;    // wq |J_F|
         for (int s = 0; s < n; s++)
         {
            for (int r = 0; r < n; r++)
            {
               real_t v = t0 * (N(r) * D(s)) + t1 * (D(r) * N(s));
               if (sw != 0.0) { v += sw * (N(r) * N(s)); }
               elmat(r, s) += W * v;
            }
         }
      }
   }

   /// sigma_w of eq. (6.5) at one face point, with c in place of M (the
   /// fluid block has c = M):
   ///   interior: kappa_w/4 sum_i (D+1) p(p+D-1)/D * D|J_F|/|J_K_i| * c_i,
   ///   boundary: kappa_w (D+1) p(p+D-1)/D * D|J_F|/|J_K_1| * c_1.
   /// with D = dim (3 on tetrahedra), p the larger order of the two sides,
   /// nl_q = |J_F| and detJ_i = |J_K_i|.  The arithmetic is that of the
   /// elastic ComputePenalty with (D lambda + 2 mu)^2 / (2 mu) replaced by c.
   real_t Penalty(int dim, int p, real_t nl_q, real_t detJ1, real_t c1,
                  real_t detJ2, real_t c2, bool is_interior) const
   {
      const real_t c_N_1 = p * (p + dim - 1.0) / dim;
      const real_t p0 = (dim + 1) * c_N_1 * (real_t(dim) * nl_q / detJ1) * c1;
      if (!is_interior) { return kappa_w_ * p0; }
      const real_t p1 = (dim + 1) * c_N_1 * (real_t(dim) * nl_q / detJ2) * c2;
      return kappa_w_ * (p0 + p1) / 4.0;
   }

   /// Element part of a source load of eq. (6.6): int_K c s div b.
   void AssembleSourceRHS(const FiniteElement &el, ElementTransformation &Tr,
                          Coefficient &s, Vector &elvect) const
   {
      const int dim = el.GetDim();
      const int nd = el.GetDof();
      elvect.SetSize(dim * nd);
      elvect = 0.0;

      const IntegrationRule &ir =
         IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 1);

      DenseMatrix dshape(nd, dim), dshape_adj(nd, dim), adjJ(dim);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcDShape(ip, dshape);
         CalcAdjugate(Tr.Jacobian(), adjJ);
         Mult(dshape, adjJ, dshape_adj);        // grad phi * detJ = grad phi |dx/dxi|
         const real_t W = ip.weight * c_->Eval(Tr, ip) * s.Eval(Tr, ip);
         for (int p = 0; p < dim; p++)
         {
            for (int k = 0; k < nd; k++) { elvect(p * nd + k) += W * dshape_adj(k, p); }
         }
      }
   }

   /// Face part of a source load of eq. (6.6).  Interior face:
   /// -int_F {c s} [b.n], side 1's rows then side 2's (on a shared face
   /// ParLinearForm keeps side 1's).  Boundary face: -int_F c s (b.n).
   void AssembleSourceRHS(const FiniteElement &el1, const FiniteElement &el2,
                          FaceElementTransformations &Trans, Coefficient &s,
                          Vector &elvect) const
   {
      const int dim = el1.GetDim();
      const int nd1 = el1.GetDof();
      const int nd2 = (Trans.Elem2No >= 0) ? el2.GetDof() : 0;
      const bool is_interior = (nd2 > 0);
      const int n1 = dim * nd1;
      const int n = dim * (nd1 + nd2);
      elvect.SetSize(n);
      elvect = 0.0;

      const int p = std::max(el1.GetOrder(), is_interior ? el2.GetOrder() : 0);
      const IntegrationRule &ir = IntRules.Get(Trans.GetGeometryType(), 2 * p + 1);

      Vector N(n), D(n), nor(dim), nhat(dim);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);
         CalcOrtho(Trans.Jacobian(), nor);
         const real_t nl_q = nor.Norml2();
         nhat.Set(1.0 / nl_q, nor);

         real_t c1, detJ1, c2 = 0.0, detJ2 = 0.0;
         const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
         SideRows_(el1, *Trans.Elem1, eip1, nhat, 1.0, N, D, 0, c1, detJ1);
         // {c s} on an interior face, c s on a boundary face.
         real_t cs = c1 * s.Eval(*Trans.Elem1, eip1);
         if (is_interior)
         {
            const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();
            SideRows_(el2, *Trans.Elem2, eip2, nhat, -1.0, N, D, n1, c2, detJ2);
            cs = 0.5 * (cs + c2 * s.Eval(*Trans.Elem2, eip2));
         }
         const real_t W = -ip.weight * nl_q * cs;
         for (int r = 0; r < n; r++) { elvect(r) += W * N(r); }
      }
   }

   /// Dirichlet load of the kTrialNormal rows on a boundary face:
   /// -int_F c (div b)(g.n), with g evaluated on side 1 at the face points,
   /// as the elastic Dirichlet load evaluates it (2p + 1 points).
   void AssembleDirichletRHS(const FiniteElement &el,
                             FaceElementTransformations &Trans,
                             VectorCoefficient &g, Vector &elvect) const
   {
      MFEM_VERIFY(Trans.Elem2No < 0, "DGPoroDivIntegrator::AssembleDirichletRHS: "
                  "an interior face; the Dirichlet load belongs to boundary faces.");
      const int dim = el.GetDim();
      const int n = dim * el.GetDof();
      elvect.SetSize(n);
      elvect = 0.0;

      const IntegrationRule &ir =
         IntRules.Get(Trans.GetGeometryType(), 2 * el.GetOrder() + 1);

      Vector N(n), D(n), nor(dim), nhat(dim), gq(dim);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);
         CalcOrtho(Trans.Jacobian(), nor);
         const real_t nl_q = nor.Norml2();
         nhat.Set(1.0 / nl_q, nor);

         real_t c1, detJ1;
         const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
         SideRows_(el, *Trans.Elem1, eip1, nhat, 1.0, N, D, 0, c1, detJ1);
         g.Eval(gq, *Trans.Elem1, eip1);
         const real_t W = -ip.weight * nl_q * (gq * nhat);
         for (int r = 0; r < n; r++) { elvect(r) += W * D(r); }
      }
   }

private:
   /// One side's rows at a face point: N = s_i phi_k nhat_p and
   /// D = c_i div(phi_k e_p), at offset off; also returns c_i and detJ_i.
   void SideRows_(const FiniteElement &el, ElementTransformation &T,
                  const IntegrationPoint &eip, const Vector &nhat, real_t sign,
                  Vector &N, Vector &D, int off, real_t &c, real_t &detJ) const
   {
      const int dim = el.GetDim();
      const int nd = el.GetDof();
      Vector shape(nd);
      DenseMatrix dshape(nd, dim), dshape_adj(nd, dim), adjJ(dim);
      el.CalcShape(eip, shape);
      el.CalcDShape(eip, dshape);
      CalcAdjugate(T.Jacobian(), adjJ);
      Mult(dshape, adjJ, dshape_adj);
      detJ = T.Weight();
      c = c_->Eval(T, eip);
      const real_t cd = c / detJ;
      for (int p = 0; p < dim; p++)
      {
         for (int k = 0; k < nd; k++)
         {
            N(off + p * nd + k) = sign * shape(k) * nhat(p);
            D(off + p * nd + k) = cd * dshape_adj(k, p);
         }
      }
   }

   Coefficient *c_;
   real_t kappa_w_;
   int bdr_terms_;
};

/// A source load of eq. (6.6) as a LinearFormIntegrator, so that
/// ParLinearForm assembles it with MFEM's element, interior-face (local and
/// shared) and marked boundary-face loops.  The form owns and deletes each
/// integrator it is given, so each Add...Integrator call needs its own
/// instance.  B and s are kept, not owned.
class DGPoroDivSourceIntegrator : public LinearFormIntegrator
{
public:
   DGPoroDivSourceIntegrator(const DGPoroDivIntegrator &B, Coefficient &s)
      : B_(B), s_(s) {}

   void AssembleRHSElementVect(const FiniteElement &el,
                               ElementTransformation &Tr,
                               Vector &elvect) override
   { B_.AssembleSourceRHS(el, Tr, s_, elvect); }

   void AssembleRHSElementVect(const FiniteElement &el,
                               FaceElementTransformations &Tr,
                               Vector &elvect) override
   { B_.AssembleSourceRHS(el, el, Tr, s_, elvect); }

   void AssembleRHSElementVect(const FiniteElement &el1,
                               const FiniteElement &el2,
                               FaceElementTransformations &Tr,
                               Vector &elvect) override
   { B_.AssembleSourceRHS(el1, el2, Tr, s_, elvect); }

private:
   const DGPoroDivIntegrator &B_;
   Coefficient &s_;
};

/// The Dirichlet load AssembleDirichletRHS as a LinearFormIntegrator, for
/// ParLinearForm::AddBdrFaceIntegrator with the Dirichlet faces' marker (one
/// instance per call, as above).  B and g are kept, not owned.
class DGPoroDivDirichletIntegrator : public LinearFormIntegrator
{
public:
   DGPoroDivDirichletIntegrator(const DGPoroDivIntegrator &B, VectorCoefficient &g)
      : B_(B), g_(g) {}

   void AssembleRHSElementVect(const FiniteElement &, ElementTransformation &,
                               Vector &) override
   {
      MFEM_ABORT("DGPoroDivDirichletIntegrator: a boundary-face load; add it "
                 "with AddBdrFaceIntegrator.");
   }

   void AssembleRHSElementVect(const FiniteElement &el,
                               FaceElementTransformations &Tr,
                               Vector &elvect) override
   { B_.AssembleDirichletRHS(el, Tr, g_, elvect); }

private:
   const DGPoroDivIntegrator &B_;
   VectorCoefficient &g_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_PORO_DIV_INTEGRATOR_HPP
