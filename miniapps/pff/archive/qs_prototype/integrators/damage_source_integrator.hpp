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

#ifndef MFEM_PFF_DAMAGE_SOURCE_INTEGRATOR_HPP
#define MFEM_PFF_DAMAGE_SOURCE_INTEGRATOR_HPP

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../materials/degradation_function.hpp"

namespace mfem
{

namespace pff
{

/** @brief Nonlinear integrator for the damage source term.
 *
 * Implements the weak form:
 *   ∫_Ω [2·Gc/(c₀·l)·d + g'(d)·H] δd dΩ
 *
 * where:
 *   - d is the damage variable
 *   - H is the strain energy history (driving force)
 *   - g(d) = (1-d)^p * (1-η) + η is the degradation function
 *   - g'(d) = -p * (1-d)^(p-1) * (1-η)
 *
 * The Jacobian (element gradient) is:
 *   ∫_Ω [2·Gc/(c₀·l) + g''(d)·H] φ_i φ_j dΩ
 */
class DamageSourceIntegrator : public NonlinearFormIntegrator
{
public:
   /** @brief Construct damage source integrator.
    *
    * @param[in] mat Material parameters
    * @param[in] H Strain energy history field
    */
   DamageSourceIntegrator(const PFFMaterialParameters &mat,
                          const GridFunction &H)
      : mat_(mat), H_(H), g_(mat.eta, mat.p),
        reaction_coeff_(mat.GetDamageReactionCoeff()) {}

   /** @brief Assemble element residual vector.
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for d
    * @param[out] elvect Element residual vector
    */
   void AssembleElementVector(const FiniteElement &el,
                              ElementTransformation &Tr,
                              const Vector &elfun,
                              Vector &elvect) override
   {
      int nd = el.GetDof();
      elvect.SetSize(nd);
      elvect = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderW();
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      Vector shape(nd);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcShape(ip, shape);

         // Evaluate d at quadrature point
         real_t d_val = 0.0;
         for (int j = 0; j < nd; j++)
         {
            d_val += shape(j) * elfun(j);
         }

         // Clamp d to [0, 1]
         d_val = std::max(0.0, std::min(1.0, d_val));

         // Evaluate H at quadrature point
         real_t H_val = H_.GetValue(Tr, ip);

         // Compute source term: 2*Gc/(c0*l)*d + g'(d)*H
         real_t g_prime = g_.EvalDerivative(d_val);
         real_t source = reaction_coeff_ * d_val + g_prime * H_val;

         real_t w = ip.weight * Tr.Weight();

         // Add contribution to element vector
         for (int j = 0; j < nd; j++)
         {
            elvect(j) += w * source * shape(j);
         }
      }
   }

   /** @brief Assemble element Jacobian matrix.
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for d
    * @param[out] elmat Element Jacobian matrix
    */
   void AssembleElementGrad(const FiniteElement &el,
                            ElementTransformation &Tr,
                            const Vector &elfun,
                            DenseMatrix &elmat) override
   {
      int nd = el.GetDof();
      elmat.SetSize(nd, nd);
      elmat = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderW();
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      Vector shape(nd);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcShape(ip, shape);

         // Evaluate d at quadrature point
         real_t d_val = 0.0;
         for (int j = 0; j < nd; j++)
         {
            d_val += shape(j) * elfun(j);
         }

         // Clamp d to [0, 1]
         d_val = std::max(0.0, std::min(1.0, d_val));

         // Evaluate H at quadrature point
         real_t H_val = H_.GetValue(Tr, ip);

         // Compute Jacobian coefficient: 2*Gc/(c0*l) + g''(d)*H
         real_t g_double_prime = g_.EvalSecondDerivative(d_val);
         real_t jac_coeff = reaction_coeff_ + g_double_prime * H_val;

         real_t w = ip.weight * Tr.Weight() * jac_coeff;

         // Add mass-like term: ∫ coeff * φ_i * φ_j dΩ
         AddMult_a_VVt(w, shape, elmat);
      }
   }

   /** @brief Get element energy for this integrator.
    *
    * Energy contribution: ∫_Ω [Gc/(c0*l) * d² + g(d) * H] dΩ
    *
    * @param[in] el Finite element
    * @param[in] Tr Element transformation
    * @param[in] elfun Element DOF values for d
    * @return Element energy
    */
   real_t GetElementEnergy(const FiniteElement &el,
                           ElementTransformation &Tr,
                           const Vector &elfun) override
   {
      int nd = el.GetDof();
      real_t energy = 0.0;

      const IntegrationRule *ir = IntRule;
      if (ir == nullptr)
      {
         int order = 2 * el.GetOrder() + Tr.OrderW();
         ir = &IntRules.Get(el.GetGeomType(), order);
      }

      Vector shape(nd);

      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Tr.SetIntPoint(&ip);

         el.CalcShape(ip, shape);

         // Evaluate d at quadrature point
         real_t d_val = 0.0;
         for (int j = 0; j < nd; j++)
         {
            d_val += shape(j) * elfun(j);
         }

         // Clamp d to [0, 1]
         d_val = std::max(0.0, std::min(1.0, d_val));

         // Evaluate H at quadrature point
         real_t H_val = H_.GetValue(Tr, ip);

         // Energy: Gc/(c0*l) * d^2 + g(d) * H
         real_t g_val = g_.Eval(d_val);
         real_t local_energy = 0.5 * reaction_coeff_ * d_val * d_val + g_val * H_val;

         real_t w = ip.weight * Tr.Weight();
         energy += w * local_energy;
      }

      return energy;
   }

private:
   const PFFMaterialParameters &mat_;
   const GridFunction &H_;
   DegradationFunction g_;
   real_t reaction_coeff_;  ///< 2 * Gc / (c0 * l)
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_DAMAGE_SOURCE_INTEGRATOR_HPP
