// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_PFF_DAMAGE_DIFFUSION_INTEGRATOR_HPP
#define MFEM_PFF_DAMAGE_DIFFUSION_INTEGRATOR_HPP

#include "mfem.hpp"

namespace mfem
{

namespace pff
{

/** @brief Bilinear integrator for the damage diffusion term.
 *
 * Implements the bilinear form:
 *   ∫_Ω (Gc·l/c₀) ∇d · ∇δd dΩ
 *
 * This is equivalent to a standard diffusion integrator with coefficient
 * k = Gc * l / c0.
 *
 * For the AT2 model: c0 = 2.
 */
class DamageDiffusionIntegrator : public BilinearFormIntegrator
{
public:
   /** @brief Construct damage diffusion integrator.
    *
    * @param[in] Gc Fracture toughness [J/m^2]
    * @param[in] l Regularization length [m]
    * @param[in] c0 Normalization constant (default: 2.0 for AT2)
    */
   DamageDiffusionIntegrator(real_t Gc, real_t l, real_t c0 = 2.0)
      : coeff_(Gc * l / c0) {}

   /// Get the diffusion coefficient
   real_t GetCoefficient() const { return coeff_; }

   /** @brief Assemble element matrix for the diffusion term.
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[out] elmat Element matrix
    */
   void AssembleElementMatrix(const FiniteElement &el,
                              ElementTransformation &Tr,
                              DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      int dim = el.GetDim();

      elmat.SetSize(nd, nd);
      elmat = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderGrad(&el);
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      DenseMatrix dshape(nd, dim);
      DenseMatrix dshapedxt(nd, dim);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcDShape(ip, dshape);

         // Transform gradient to physical coordinates
         Mult(dshape, Tr.InverseJacobian(), dshapedxt);

         real_t w = ip.weight * Tr.Weight() * coeff_;

         // Add ∫ k * ∇φ_i · ∇φ_j dΩ
         AddMult_a_AAt(w, dshapedxt, elmat);
      }
   }

private:
   real_t coeff_;  ///< Diffusion coefficient: Gc * l / c0
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_DAMAGE_DIFFUSION_INTEGRATOR_HPP
