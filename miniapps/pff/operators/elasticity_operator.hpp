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

#ifndef MFEM_PFF_ELASTICITY_OPERATOR_HPP
#define MFEM_PFF_ELASTICITY_OPERATOR_HPP

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../integrators/degraded_elasticity_integrator.hpp"

namespace mfem
{

namespace pff
{

/** @brief Operator wrapper for the degraded elasticity residual and Jacobian.
 *
 * R(u) = ∫_Ω σ(u,d) : ε(v) dΩ,  with σ = g(d) σ⁺ + σ⁻
 *
 * Boundary conditions are handled by enforcing:
 *   - residual on essential true dofs: r_i = u_i - u_prescribed_i
 *   - gradient rows/cols eliminated with identity on essential true dofs
 */
class ElasticityOperator : public Operator
{
public:
   /** @brief Construct elasticity operator.
    *
    * @param[in] fes Finite element space for displacement (vector H1)
    * @param[in] mat Material parameters
    * @param[in] ess_bdr Array marking essential boundary attributes
    */
   ElasticityOperator(ParFiniteElementSpace &fes,
                      const PFFMaterialParameters &mat,
                      Array<int> &ess_bdr)
      : Operator(fes.TrueVSize()), fes_(fes), mat_(mat),
        damage_gf_(nullptr), nlf_(nullptr), u_prescribed_(height)
   {
      // Build nonlinear form with degraded elasticity integrator
      nlf_ = new ParNonlinearForm(&fes_);

      // Compute essential true dofs from boundary markers
      ess_tdof_list_.SetSize(0);
      fes_.GetEssentialTrueDofs(ess_bdr, ess_tdof_list_);
      nlf_->SetEssentialTrueDofs(ess_tdof_list_);

      // Default prescribed displacement is zero
      u_prescribed_ = 0.0;
   }

   ~ElasticityOperator() override { delete nlf_; }

   /** @brief Set the current damage field.
    *
    * @param[in] d Damage GridFunction (defined on a compatible FE space)
    */
   void SetDamageField(const ParGridFunction &d)
   {
      // Recreate the nonlinear form so the integrator holds a reference to 'd'
      delete nlf_;
      nlf_ = new ParNonlinearForm(&fes_);
      nlf_->AddDomainIntegrator(new DegradedElasticityIntegrator(mat_, d));
      nlf_->SetEssentialTrueDofs(ess_tdof_list_);
      damage_gf_ = &d;
   }

   /** @brief Compute residual R(u) = internal forces with BC enforcement.
    *
    * On essential true dofs: r_i = u_i - u_prescribed_i
    *
    * @param[in] u_true Displacement true DOF vector
    * @param[out] r Residual vector
    */
   void Mult(const Vector &u_true, Vector &r) const override
   {
      MFEM_VERIFY(nlf_ != nullptr, "Nonlinear form not initialized");
      r.SetSize(height);
      nlf_->Mult(u_true, r);

      // Enforce Dirichlet on essential true dofs: r_i = u_i - u_prescribed_i
      if (ess_tdof_list_.Size() > 0)
      {
         for (int i = 0; i < ess_tdof_list_.Size(); i++)
         {
            const int tdof = ess_tdof_list_[i];
            r(tdof) = u_true(tdof) - u_prescribed_(tdof);
         }
      }
   }

   /** @brief Get gradient K = dR/du with elimination on essential rows/cols.
    *
    * @param[in] u_true Displacement true DOF vector
    * @return Reference to Jacobian operator
    */
   Operator &GetGradient(const Vector &u_true) const override
   {
      MFEM_VERIFY(nlf_ != nullptr, "Nonlinear form not initialized");
      Operator &J = nlf_->GetGradient(u_true);

      // Eliminate essential rows/cols (identity on diagonal)
      if (ess_tdof_list_.Size() > 0)
      {
         HypreParMatrix *J_hyp = dynamic_cast<HypreParMatrix*>(&J);
         if (J_hyp)
         {
            HypreParMatrix *Je = J_hyp->EliminateRowsCols(ess_tdof_list_);
            delete Je; // not needed further
         }
         else
         {
            // Serial fallback (mainly for unit tests)
            SparseMatrix *J_sp = dynamic_cast<SparseMatrix*>(&J);
            if (J_sp)
            {
               J_sp->EliminateBC(ess_tdof_list_, Operator::DiagonalPolicy::DIAG_ONE);
            }
         }
      }
      return J;
   }

   /** @brief Compute positive strain energy density and project to psi_active.
    *
    * This computes ψ⁺ (tensile strain energy density) at quadrature points
    * and averages it over each element.
    *
    * @param[in] u_true Displacement true DOF vector
    * @param[out] psi_active GridFunction to store ψ⁺ values
    */
   void ComputeStrainEnergyDensity(const Vector &u_true,
                                   ParGridFunction &psi_active) const
   {
      // We evaluate element-wise ψ⁺ using the integrator and set a constant
      // value per element in psi_active's FESpace.
      MFEM_VERIFY(damage_gf_ != nullptr, "Damage field must be set before use");
      DegradedElasticityIntegrator integ(mat_, *damage_gf_);

      ParFiniteElementSpace *fes = dynamic_cast<ParFiniteElementSpace*>(
         psi_active.FESpace());
      MFEM_VERIFY(fes != nullptr, "psi_active must be a ParGridFunction");

      Mesh *mesh = fes_.GetParMesh();
      MFEM_VERIFY(mesh != nullptr, "ParMesh required");

      // Build a GridFunction for u on the FE space
      ParGridFunction u_gf(&fes_);
      u_gf.SetFromTrueDofs(u_true);

      psi_active = 0.0;
      Vector elfun, qp_vals;
      Array<int> vdofs;

      const int NE = fes_.GetNE();
      for (int el = 0; el < NE; el++)
      {
         const FiniteElement &fel = *fes_.GetFE(el);
         ElementTransformation &Tr = *fes_.GetElementTransformation(el);
         fes_.GetElementVDofs(el, vdofs);
         u_gf.GetSubVector(vdofs, elfun);

         // Compute ψ⁺ at quadrature points and average over element
         integ.ComputePositiveStrainEnergy(fel, Tr, elfun, qp_vals);

         // Weighted average with integration weights
         const IntegrationRule *ir = &IntRules.Get(fel.GetGeomType(),
                                                   2 * fel.GetOrder() + 3);
         MFEM_ASSERT(ir->GetNPoints() == qp_vals.Size(),
                     "Quadrature size mismatch");

         real_t int_val = 0.0, meas = 0.0;
         for (int i = 0; i < ir->GetNPoints(); i++)
         {
            const IntegrationPoint &ip = ir->IntPoint(i);
            Tr.SetIntPoint(&ip);
            const real_t w = ip.weight * Tr.Weight();
            int_val += w * qp_vals(i);
            meas += w;
         }
         const real_t avg = (meas > 0.0) ? (int_val / meas) : 0.0;

         // Set element values in psi_active (assume constant per element)
         Array<int> ddofs;
         psi_active.FESpace()->GetElementDofs(el, ddofs);
         for (int i = 0; i < ddofs.Size(); i++)
         {
            psi_active(ddofs[i]) = avg;
         }
      }
      // For L2/DG spaces this is element-local; for H1, MFEM handles assembly.
   }

   /** @brief Set essential true dofs list explicitly.
    *
    * @param[in] ess_tdof_list Array of essential true DOF indices
    */
   void SetEssentialTrueDofs(const Array<int> &ess_tdof_list)
   {
      ess_tdof_list_.SetSize(ess_tdof_list.Size());
      ess_tdof_list_.Assign(ess_tdof_list);
      nlf_->SetEssentialTrueDofs(ess_tdof_list_);
   }

   /** @brief Set prescribed displacement on essential true dofs.
    *
    * @param[in] u_prescribed Vector of prescribed values (size = TrueVSize)
    */
   void SetPrescribedDisplacement(const Vector &u_prescribed)
   {
      MFEM_VERIFY(u_prescribed.Size() == height,
                  "Prescribed displacement size mismatch");
      u_prescribed_ = u_prescribed;
   }

   /// Get essential true dofs list
   const Array<int> &GetEssentialTrueDofs() const { return ess_tdof_list_; }

private:
   ParFiniteElementSpace &fes_;
   const PFFMaterialParameters &mat_;
   mutable const ParGridFunction *damage_gf_; // pointer to external damage field
   mutable ParNonlinearForm *nlf_;            // internal nonlinear form
   Array<int> ess_tdof_list_;                 // essential true dof list
   mutable Vector u_prescribed_;              // prescribed displacement (true dofs)
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_ELASTICITY_OPERATOR_HPP
