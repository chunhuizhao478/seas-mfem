// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// split_elasticity_integrator.hpp — small-strain elasticity with a degraded
// tension/compression split, as a NonlinearFormIntegrator for H1 vector
// spaces (byNODES element layout).
//
// MOOSE/RACCOON: ADStressDivergenceTensors with ADComputeSmallStrain,
// SmallDeformationIsotropicElasticity (decomposition NONE / SPECTRAL /
// VOLDEV) and ComputeSmallDeformationStress.  The pointwise stress and
// energies are ComputeSplit (constitutive/strain_split.hpp); 2D is plane
// strain with lambda = K - 2G/3 (LIBMESH_DIM = 3).  g = g(d) with d
// interpolated at the quadrature point from the nodal damage field.
//
// Residual  R_a,i = sum_qp w |J| sigma_ij(eps(u), d) dN_a/dx_j
// Tangent   d R / d u: analytic for NONE (g C) and VOLDEV
//           (g C + (1 - g) K H(-tr eps) I x I); SPECTRAL through
//           SplitTangentApply (derivative of the eigenprojection), checked
//           against finite differences in tests/unit/test_strain_split.cpp.

#ifndef MFEM_PFF_SPLIT_ELASTICITY_INTEGRATOR_HPP
#define MFEM_PFF_SPLIT_ELASTICITY_INTEGRATOR_HPP

#include "mfem.hpp"
#include "../constitutive/strain_split.hpp"
#include "../phasefield/phase_field_functions.hpp"

namespace mfem
{
namespace pff
{

class SplitElasticityIntegrator : public NonlinearFormIntegrator
{
public:
   SplitElasticityIntegrator(real_t K, real_t G, Decomposition dec,
                             const PowerDegradation &deg, const GridFunction &d)
      : K_(K), G_(G), dec_(dec), deg_(deg), d_(d) {}

   /// Strain (3x3, plane strain padded) at the current integration point of
   /// Tr from the element dofs; dshape is the physical shape gradient.
   static void Strain(int dim, int dof, const DenseMatrix &dshape, const Vector &elfun,
                      real_t *eps3)
   {
      real_t gu[9] = {0.0};
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++)
         {
            real_t s = 0.0;
            for (int a = 0; a < dof; a++) { s += elfun(a + i * dof) * dshape(a, j); }
            gu[i + 3 * j] = s;
         }
      for (int i = 0; i < 9; i++) { eps3[i] = 0.0; }
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++)
         {
            eps3[i + 3 * j] = 0.5 * (gu[i + 3 * j] + gu[j + 3 * i]);
         }
   }

   real_t Degradation(ElementTransformation &Tr, const IntegrationPoint &ip) const
   {
      return deg_.g(d_.GetValue(Tr, ip));
   }

   void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr,
                              const Vector &elfun, Vector &elvect) override
   {
      const int dof = el.GetDof(), dim = el.GetDim();
      dshape_.SetSize(dof, dim);
      elvect.SetSize(dof * dim);
      elvect = 0.0;
      const IntegrationRule &ir = Rule(el, Tr);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysDShape(Tr, dshape_);
         real_t eps[9];
         Strain(dim, dof, dshape_, elfun, eps);
         SplitResult r;
         ComputeSplit(dec_, eps, K_, G_, Degradation(Tr, ip), r);
         const real_t w = ip.weight * Tr.Weight();
         for (int a = 0; a < dof; a++)
            for (int i = 0; i < dim; i++)
            {
               real_t s = 0.0;
               for (int j = 0; j < dim; j++) { s += r.sigma[i + 3 * j] * dshape_(a, j); }
               elvect(a + i * dof) += w * s;
            }
      }
   }

   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr,
                            const Vector &elfun, DenseMatrix &elmat) override
   {
      if (dec_ == Decomposition::Spectral)
      {
         AssembleSpectralGrad(el, Tr, elfun, elmat);
         return;
      }
      const int dof = el.GetDof(), dim = el.GetDim();
      dshape_.SetSize(dof, dim);
      elmat.SetSize(dof * dim);
      elmat = 0.0;
      const real_t lambda = K_ - 2.0 * G_ / 3.0;
      const IntegrationRule &ir = Rule(el, Tr);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysDShape(Tr, dshape_);
         const real_t g = Degradation(Tr, ip);
         // C_ijkl = g (lambda d_ij d_kl + G (d_ik d_jl + d_il d_jk)) + c_vol d_ij d_kl
         real_t c_vol = 0.0;
         if (dec_ == Decomposition::VolDev)
         {
            real_t eps[9];
            Strain(dim, dof, dshape_, elfun, eps);
            const real_t tr = eps[0] + eps[4] + eps[8];
            if (tr < 0.0) { c_vol = (1.0 - g) * K_; }
         }
         const real_t lam = g * lambda + c_vol, mu = g * G_;
         const real_t w = ip.weight * Tr.Weight();
         for (int a = 0; a < dof; a++)
            for (int b = 0; b < dof; b++)
            {
               real_t gg = 0.0;
               for (int j = 0; j < dim; j++) { gg += dshape_(a, j) * dshape_(b, j); }
               for (int i = 0; i < dim; i++)
                  for (int k = 0; k < dim; k++)
                  {
                     // sum_jl dN_a/dx_j C_ijkl dN_b/dx_l
                     real_t v = lam * dshape_(a, i) * dshape_(b, k) +
                                mu * dshape_(a, k) * dshape_(b, i);
                     if (i == k) { v += mu * gg; }
                     elmat(a + i * dof, b + k * dof) += w * v;
                  }
            }
      }
   }

   /// Element averages of psi_e^active (the CONSTANT MONOMIAL aux variable
   /// RACCOON transfers to the fracture app): sum w|J| psi+ / sum w|J|.
   void ElementAverageActiveEnergy(const FiniteElement &el, ElementTransformation &Tr,
                                   const Vector &elfun, real_t &avg)
   {
      const int dof = el.GetDof(), dim = el.GetDim();
      dshape_.SetSize(dof, dim);
      const IntegrationRule &ir = Rule(el, Tr);
      real_t num = 0.0, den = 0.0;
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysDShape(Tr, dshape_);
         real_t eps[9];
         Strain(dim, dof, dshape_, elfun, eps);
         SplitResult r;
         ComputeSplit(dec_, eps, K_, G_, 1.0, r);
         const real_t w = ip.weight * Tr.Weight();
         num += w * r.psi_active;
         den += w;
      }
      avg = num / den;
   }

private:
   void AssembleSpectralGrad(const FiniteElement &el, ElementTransformation &Tr,
                             const Vector &elfun, DenseMatrix &elmat)
   {
      const int dof = el.GetDof(), dim = el.GetDim();
      dshape_.SetSize(dof, dim);
      elmat.SetSize(dof * dim);
      elmat = 0.0;
      const IntegrationRule &ir = Rule(el, Tr);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysDShape(Tr, dshape_);
         real_t eps[9];
         Strain(dim, dof, dshape_, elfun, eps);
         const real_t g = Degradation(Tr, ip);
         const real_t w = ip.weight * Tr.Weight();
         for (int b = 0; b < dof; b++)
            for (int k = 0; k < dim; k++)
            {
               // deps = sym(e_k (x) grad N_b)
               real_t de[9] = {0.0};
               for (int l = 0; l < dim; l++)
               {
                  de[k + 3*l] += 0.5 * dshape_(b, l);
                  de[l + 3*k] += 0.5 * dshape_(b, l);
               }
               real_t ds[9];
               SplitTangentApply(dec_, eps, K_, G_, g, de, ds);
               for (int a = 0; a < dof; a++)
                  for (int i = 0; i < dim; i++)
                  {
                     real_t v = 0.0;
                     for (int j = 0; j < dim; j++) { v += ds[i + 3*j] * dshape_(a, j); }
                     elmat(a + i * dof, b + k * dof) += w * v;
                  }
            }
      }
   }

   const IntegrationRule &Rule(const FiniteElement &el, ElementTransformation &Tr) const
   {
      if (IntRule) { return *IntRule; }
      return IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());
   }

   real_t K_, G_;
   Decomposition dec_;
   PowerDegradation deg_;
   const GridFunction &d_;
   DenseMatrix dshape_;
};

} // namespace pff
} // namespace mfem

#endif
