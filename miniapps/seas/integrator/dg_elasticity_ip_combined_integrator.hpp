// Combined IP DG elasticity face integrator matching Tandem's structure.
//
// Computes the full skeleton/boundary face matrix (consistency + symmetry +
// penalty) in one pass. Also provides AssembleSlipFaceRHS() for the fault
// slip RHS contribution using the SAME traction operator formula, guaranteeing
// K-b consistency by construction.
//
// Reference: Tandem app/kernels/elasticity.py
//   - assembleSurface (bilinear form)
//   - rhsFacet (slip RHS)
//   - tractionTest / test_normal (shared helpers)
//   - lift_ip (IP penalty lifting)

#ifndef MFEM_SEAS_DG_ELASTICITY_IP_COMBINED_INTEGRATOR_HPP
#define MFEM_SEAS_DG_ELASTICITY_IP_COMBINED_INTEGRATOR_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

class DGElasticityIPCombinedIntegrator : public BilinearFormIntegrator
{
public:
   /// @param lambda First Lame parameter
   /// @param mu Shear modulus
   /// @param dim Spatial dimension (3 for BP5)
   /// @param epsilon SIPG sign (-1 for SIPG)
   /// @param penalty_factor Multiplier on penalty (default 1.0)
   DGElasticityIPCombinedIntegrator(Coefficient &lambda, Coefficient &mu,
                                     int dim, real_t epsilon,
                                     real_t penalty_factor = 1.0)
      : lambda_(lambda), mu_(mu), dim_(dim),
        epsilon_(epsilon), penalty_factor_(penalty_factor) {}

   using BilinearFormIntegrator::AssembleFaceMatrix;

   /// Assemble the full face matrix: consistency + symmetry + penalty.
   /// Mirrors Tandem's assembleSurface kernel.
   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override
   {
      const int ndof1 = el1.GetDof();
      const int ndof2 = (Trans.Elem2No >= 0) ? el2.GetDof() : 0;
      const int ndofs = ndof1 + ndof2;
      const int nvdofs = dim_ * ndofs;
      const bool is_interior = (ndof2 > 0);

      elmat.SetSize(nvdofs);
      elmat = 0.0;

      // Quadrature: use 2p+1 matching Tandem's MinQuadOrder for all terms
      const int order = 2 * std::max(el1.GetOrder(),
                                      ndof2 ? el2.GetOrder() : 0) + 1;
      const IntegrationRule &ir = IntRules.Get(Trans.GetGeometryType(), order);

      Vector shape1(ndof1), shape2(std::max(ndof2, 1));
      DenseMatrix dshape1_ref(ndof1, dim_), dshape2_ref(std::max(ndof2,1), dim_);
      DenseMatrix dshape1_adj(ndof1, dim_), dshape2_adj(std::max(ndof2,1), dim_);
      DenseMatrix adjJ(dim_);
      Vector nor(dim_);

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);
         const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();

         el1.CalcShape(eip1, shape1);
         el1.CalcDShape(eip1, dshape1_ref);
         CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);
         Mult(dshape1_ref, adjJ, dshape1_adj);

         CalcOrtho(Trans.Jacobian(), nor);
         real_t nl_q = nor.Norml2();

         real_t detJ1 = Trans.Elem1->Weight();
         real_t lam1 = lambda_.Eval(*Trans.Elem1, eip1);
         real_t mu1 = mu_.Eval(*Trans.Elem1, eip1);

         real_t detJ2 = detJ1;
         real_t lam2 = lam1, mu2 = mu1;

         if (is_interior)
         {
            const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();
            el2.CalcShape(eip2, shape2);
            el2.CalcDShape(eip2, dshape2_ref);
            CalcAdjugate(Trans.Elem2->Jacobian(), adjJ);
            Mult(dshape2_ref, adjJ, dshape2_adj);

            detJ2 = Trans.Elem2->Weight();
            lam2 = lambda_.Eval(*Trans.Elem2, eip2);
            mu2 = mu_.Eval(*Trans.Elem2, eip2);
         }

         // Penalty (Tandem formula)
         real_t penalty = ComputePenalty(el1, el2, detJ1, detJ2,
                                         lam1, mu1, nl_q, is_interior);

         // Coefficients matching Tandem's assembleSurface
         // Skeleton: c0=-0.5, c1=eps/2, c2=+penalty (same side)
         //           c0=+0.5, c1=-eps/2, c2=-penalty (cross side)
         // Boundary: c0=-1, c1=eps, c2=+penalty
         real_t w_q = ip.weight;

         // Block (0,0): test=elem1, trial=elem1
         AssembleFaceBlock(
            ndof1, ndof1, 0, 0,
            shape1, shape1,
            dshape1_adj, dshape1_adj,
            lam1, mu1, detJ1,
            lam1, mu1, detJ1,
            nor, nl_q, w_q,
            is_interior ? -0.5 : -1.0,  // c0
            is_interior ? 0.5*epsilon_ : epsilon_,  // c1
            penalty,  // c2 (same side)
            elmat);

         if (!is_interior) { continue; }

         // Block (0,1): test=elem1, trial=elem2
         AssembleFaceBlock(
            ndof1, ndof2, 0, dim_*ndof1,
            shape1, shape2,
            dshape1_adj, dshape2_adj,
            lam2, mu2, detJ2,  // trial = elem2 for consistency
            lam1, mu1, detJ1,  // test = elem1 for symmetry (but c1 uses trial side)
            nor, nl_q, w_q,
            -0.5,               // c0 (same as block 00 — test is elem1)
            -0.5*epsilon_,       // c1 (negated — trial is elem2)
            -penalty,            // c2 (cross side)
            elmat);

         // Block (1,0): test=elem2, trial=elem1
         AssembleFaceBlock(
            ndof2, ndof1, dim_*ndof1, 0,
            shape2, shape1,
            dshape2_adj, dshape1_adj,
            lam1, mu1, detJ1,  // trial = elem1
            lam2, mu2, detJ2,  // test = elem2
            nor, nl_q, w_q,
            0.5,                // c0 (negated — test is elem2)
            0.5*epsilon_,        // c1 (trial is elem1)
            -penalty,            // c2 (cross side)
            elmat);

         // Block (1,1): test=elem2, trial=elem2
         AssembleFaceBlock(
            ndof2, ndof2, dim_*ndof1, dim_*ndof1,
            shape2, shape2,
            dshape2_adj, dshape2_adj,
            lam2, mu2, detJ2,
            lam2, mu2, detJ2,
            nor, nl_q, w_q,
            0.5,                // c0
            -0.5*epsilon_,       // c1
            penalty,             // c2 (same side)
            elmat);
      }
   }

   /// Assemble the slip RHS contribution for a single fault face.
   ///
   /// Computes: b_x = c1 * tractionTest(x, f_q) * w + c2 * φ_x * f_q * nl * w
   /// for both elements sharing the face.
   ///
   /// This uses the EXACT same traction operator as AssembleFaceMatrix,
   /// guaranteeing K-b consistency.
   ///
   /// @param fe1, fe2 Finite elements on each side
   /// @param Trans Face transformation
   /// @param slip_3d Prescribed 3D displacement jump at quad points [dim*nq]
   ///               (already sign-corrected: sign * EmbedSlip(slip))
   /// @param elvec1 Output RHS for element 1 [ndof1*dim]
   /// @param elvec2 Output RHS for element 2 [ndof2*dim]
   void AssembleSlipFaceRHS(const FiniteElement &fe1,
                             const FiniteElement &fe2,
                             FaceElementTransformations &Trans,
                             const Vector &slip_3d,
                             Vector &elvec1, Vector &elvec2) const
   {
      const int ndof1 = fe1.GetDof();
      const int ndof2 = fe2.GetDof();
      elvec1.SetSize(ndof1 * dim_); elvec1 = 0.0;
      elvec2.SetSize(ndof2 * dim_); elvec2 = 0.0;

      const int order = 2 * std::max(fe1.GetOrder(), fe2.GetOrder()) + 1;
      const IntegrationRule &ir = IntRules.Get(Trans.GetGeometryType(), order);
      int nq = ir.GetNPoints();

      Vector shape1(ndof1), shape2(ndof2);
      DenseMatrix dshape1_ref(ndof1, dim_), dshape2_ref(ndof2, dim_);
      DenseMatrix dshape1_adj(ndof1, dim_), dshape2_adj(ndof2, dim_);
      DenseMatrix adjJ(dim_);
      Vector nor(dim_);

      for (int q = 0; q < nq; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Trans.SetAllIntPoints(&ip);
         const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
         const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

         fe1.CalcShape(eip1, shape1);
         fe1.CalcDShape(eip1, dshape1_ref);
         CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);
         Mult(dshape1_ref, adjJ, dshape1_adj);

         fe2.CalcShape(eip2, shape2);
         fe2.CalcDShape(eip2, dshape2_ref);
         CalcAdjugate(Trans.Elem2->Jacobian(), adjJ);
         Mult(dshape2_ref, adjJ, dshape2_adj);

         CalcOrtho(Trans.Jacobian(), nor);
         real_t nl_q = nor.Norml2();
         real_t w_q = ip.weight;

         real_t detJ1 = Trans.Elem1->Weight();
         real_t detJ2 = Trans.Elem2->Weight();
         real_t lam1 = lambda_.Eval(*Trans.Elem1, eip1);
         real_t mu1 = mu_.Eval(*Trans.Elem1, eip1);
         real_t lam2 = lambda_.Eval(*Trans.Elem2, eip2);
         real_t mu2 = mu_.Eval(*Trans.Elem2, eip2);

         real_t penalty = ComputePenalty(fe1, fe2, detJ1, detJ2,
                                          lam1, mu1, nl_q, true);

         // Slip at this quad point
         real_t f_q[3];
         for (int c = 0; c < dim_; c++)
         {
            f_q[c] = slip_3d(c * nq + q);
         }

         // f_q dot nor (for traction operator contraction)
         real_t f_dot_n = 0.0;
         for (int d = 0; d < dim_; d++) { f_dot_n += f_q[d] * nor(d); }

         // Tandem rhsFacet coefficients:
         //   Side 0: c1 = eps/2,  c2 = +penalty
         //   Side 1: c1 = eps/2,  c2 = -penalty
         // (c1 is the SAME for both sides — Tandem rhs_skeleton lines 618,632-633)
         real_t c1 = 0.5 * epsilon_;

         // --- Element 1 RHS (side 0) ---
         // b1[k,p] += c1 * tractionTest(0, f_q) * w + penalty * φ0[k] * f_q[p] * nl * w
         for (int k = 0; k < ndof1; k++)
         {
            // tractionTest(0, f_q) contracted to [k,p]:
            // = λ1 * Dx1[k,p] * (f·n) + μ1 * (Dx1[k,:]·n * f[p] + Dx1[k,:]·f * n[p])
            real_t grad_dot_n = 0.0;
            real_t grad_dot_f = 0.0;
            for (int d = 0; d < dim_; d++)
            {
               grad_dot_n += dshape1_adj(k, d) * nor(d);
               grad_dot_f += dshape1_adj(k, d) * f_q[d];
            }

            for (int p = 0; p < dim_; p++)
            {
               real_t trac_test = lam1 * dshape1_adj(k, p) * f_dot_n
                  + mu1 * (grad_dot_n * f_q[p] + grad_dot_f * nor(p));

               int idx = p * ndof1 + k;
               elvec1(idx) += c1 * trac_test * w_q / detJ1;
               elvec1(idx) += penalty * w_q * nl_q * shape1(k) * f_q[p];
            }
         }

         // --- Element 2 RHS (side 1) ---
         // b2[k,p] += c1 * tractionTest(1, f_q) * w - penalty * φ1[k] * f_q[p] * nl * w
         for (int k = 0; k < ndof2; k++)
         {
            real_t grad_dot_n = 0.0;
            real_t grad_dot_f = 0.0;
            for (int d = 0; d < dim_; d++)
            {
               grad_dot_n += dshape2_adj(k, d) * nor(d);
               grad_dot_f += dshape2_adj(k, d) * f_q[d];
            }

            for (int p = 0; p < dim_; p++)
            {
               real_t trac_test = lam2 * dshape2_adj(k, p) * f_dot_n
                  + mu2 * (grad_dot_n * f_q[p] + grad_dot_f * nor(p));

               int idx = p * ndof2 + k;
               elvec2(idx) += c1 * trac_test * w_q / detJ2;
               elvec2(idx) -= penalty * w_q * nl_q * shape2(k) * f_q[p];
            }
         }
      }
   }

private:
   Coefficient &lambda_;
   Coefficient &mu_;
   int dim_;
   real_t epsilon_;
   real_t penalty_factor_;

   /// Compute IP penalty matching Tandem's formula.
   /// penalty = (p0+p1)/4 for interior, p0 for boundary.
   /// p_i = (D+1) * c_N_1 * (D*nl_q/detJ_i) * (c1^2/c0)
   real_t ComputePenalty(const FiniteElement &fe1,
                          const FiniteElement &fe2,
                          real_t detJ1, real_t detJ2,
                          real_t lam, real_t mu_val,
                          real_t nl_q, bool is_interior) const
   {
      int p = std::max(fe1.GetOrder(),
                       is_interior ? fe2.GetOrder() : fe1.GetOrder());
      real_t c_N_1 = p * (p + dim_ - 1.0) / dim_;
      real_t c0 = 2.0 * mu_val;
      real_t c1 = dim_ * lam + 2.0 * mu_val;
      real_t ratio = c1 * c1 / c0;

      real_t p0 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / detJ1) * ratio;
      if (!is_interior)
      {
         return penalty_factor_ * p0;
      }
      real_t p1 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / detJ2) * ratio;
      return penalty_factor_ * (p0 + p1) / 4.0;
   }

   /// Assemble one (test, trial) block of the face matrix.
   ///
   /// Adds to elmat:
   ///   c0 * φ_test[k] * w * T_op_trial[l,u,p] / detJ_trial
   /// + c1 * φ_trial[l] * w * T_op_test[k,p,u] / detJ_test
   /// + c2 * φ_test[k] * φ_trial[l] * δ[p,u] * nl * w
   void AssembleFaceBlock(
      int ndof_test, int ndof_trial,
      int offset_test, int offset_trial,
      const Vector &shape_test, const Vector &shape_trial,
      const DenseMatrix &dshape_adj_test,
      const DenseMatrix &dshape_adj_trial,
      real_t lam_trial, real_t mu_trial, real_t detJ_trial,
      real_t lam_test, real_t mu_test, real_t detJ_test,
      const Vector &nor, real_t nl_q, real_t w_q,
      real_t c0, real_t c1, real_t c2,
      DenseMatrix &elmat) const
   {
      // Precompute dshape_adj · nor for test and trial
      Vector grad_dot_n_test(ndof_test), grad_dot_n_trial(ndof_trial);
      dshape_adj_test.Mult(nor, grad_dot_n_test);
      dshape_adj_trial.Mult(nor, grad_dot_n_trial);

      for (int l = 0; l < ndof_trial; l++)
      {
         for (int u = 0; u < dim_; u++)
         {
            // T_op_trial[l,u,p]: traction operator for trial DOF l,
            // trial component u, evaluated against test component p
            for (int k = 0; k < ndof_test; k++)
            {
               for (int p = 0; p < dim_; p++)
               {
                  // Consistency: c0 * φ_test[k] * T_op_trial[l,u,p] * w / detJ_trial
                  real_t T_trial_lup =
                     lam_trial * dshape_adj_trial(l, u) * nor(p)
                     + mu_trial * ((u == p ? 1.0 : 0.0) * grad_dot_n_trial(l)
                                   + dshape_adj_trial(l, p) * nor(u));
                  real_t consistency = c0 * shape_test(k) * T_trial_lup
                                      * w_q / detJ_trial;

                  // Symmetry: c1 * φ_trial[l] * T_op_test[k,p,u] * w / detJ_test
                  real_t T_test_kpu =
                     lam_test * dshape_adj_test(k, p) * nor(u)
                     + mu_test * ((p == u ? 1.0 : 0.0) * grad_dot_n_test(k)
                                  + dshape_adj_test(k, u) * nor(p));
                  real_t symmetry = c1 * shape_trial(l) * T_test_kpu
                                    * w_q / detJ_test;

                  // Penalty: c2 * φ_test[k] * φ_trial[l] * δ[p,u] * nl * w
                  real_t pen = 0.0;
                  if (p == u)
                  {
                     pen = c2 * shape_test(k) * shape_trial(l)
                           * nl_q * w_q;
                  }

                  int row = offset_test + p * ndof_test + k;
                  int col = offset_trial + u * ndof_trial + l;
                  elmat(row, col) += consistency + symmetry + pen;
               }
            }
         }
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_ELASTICITY_IP_COMBINED_INTEGRATOR_HPP
