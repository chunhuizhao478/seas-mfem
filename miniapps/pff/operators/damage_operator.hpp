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

#ifndef MFEM_PFF_DAMAGE_OPERATOR_HPP
#define MFEM_PFF_DAMAGE_OPERATOR_HPP

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../integrators/damage_diffusion_integrator.hpp"
#include "../integrators/damage_source_integrator.hpp"

namespace mfem
{

namespace pff
{

/** @brief Operator for the damage (phase field) equation.
 *
 * Residual: R(d) = K_diff * d + R_src(d; H)
 *   where K_diff corresponds to ∫ (Gc·l/c0) ∇d·∇δd dΩ and
 *         R_src  corresponds to ∫ [2Gc/(c0 l) d + g'(d) H] δd dΩ.
 *
 * Gradient: K = K_diff + dR_src/dd
 */
class DamageOperator : public Operator
{
public:
   /** @brief Construct damage operator.
    *
    * @param[in] fes Finite element space for damage (scalar H1)
    * @param[in] mat Material parameters
    */
   DamageOperator(ParFiniteElementSpace &fes,
                  const PFFMaterialParameters &mat)
      : Operator(fes.TrueVSize()), fes_(fes), mat_(mat),
        H_gf_(&fes), nlf_(nullptr)
   {
      // Assemble the diffusion (linear) part once
      ParBilinearForm k_form(&fes_);
      k_form.AddDomainIntegrator(
         new DamageDiffusionIntegrator(mat_.Gc, mat_.l, mat_.c0));
      k_form.Assemble();
      k_form.Finalize();
      diffusion_op_.Reset(k_form.ParallelAssemble());

      // Nonlinear source term, depends on H (history)
      nlf_ = new ParNonlinearForm(&fes_);
      H_gf_ = 0.0;
      nlf_->AddDomainIntegrator(new DamageSourceIntegrator(mat_, H_gf_));
   }

   ~DamageOperator() override { delete nlf_; }

   /** @brief Update the strain energy history as the running maximum.
    *
    * H = max(H, psi_active)
    *
    * @param[in] psi_active Current positive strain energy density
    */
   void UpdateStrainEnergyHistory(const ParGridFunction &psi_active)
   {
      // Project psi_active to H_gf_ space if needed
      ParGridFunction tmp(&fes_);
      tmp.ProjectGridFunction(psi_active);

      // H = max(H, tmp)
      const int n = H_gf_.Size();
      for (int i = 0; i < n; i++)
      {
         H_gf_(i) = std::max(H_gf_(i), tmp(i));
      }
   }

   /** @brief Compute residual R(d) = K_diff * d + R_src(d).
    *
    * @param[in] d_true Damage true DOF vector
    * @param[out] r Residual vector
    */
   void Mult(const Vector &d_true, Vector &r) const override
   {
      r.SetSize(height);
      r = 0.0;

      // Linear diffusion part
      if (auto *A = diffusion_op_.Is<HypreParMatrix>())
      {
         A->Mult(d_true, r);
      }
      else if (auto *S = diffusion_op_.Is<SparseMatrix>())
      {
         S->Mult(d_true, r);
      }

      // Nonlinear source part
      Vector r_src(height);
      nlf_->Mult(d_true, r_src);
      r += r_src;
   }

   /** @brief Get gradient K = K_diff + dR_src/dd.
    *
    * @param[in] d_true Damage true DOF vector
    * @return Reference to Jacobian operator
    */
   Operator &GetGradient(const Vector &d_true) const override
   {
      Operator &Jsrc = nlf_->GetGradient(d_true);

      // Combine with diffusion matrix
      if (auto *A = diffusion_op_.Is<HypreParMatrix>())
      {
         if (auto *Js = dynamic_cast<HypreParMatrix*>(&Jsrc))
         {
            HypreParMatrix *Sum = Add(1.0, *A, 1.0, *Js);
            K_.Reset(Sum); // K_ owns Sum
            return *K_.Ptr();
         }
      }
      else if (auto *S = diffusion_op_.Is<SparseMatrix>())
      {
         if (auto *Js = dynamic_cast<SparseMatrix*>(&Jsrc))
         {
            SparseMatrix *Sum = Add(1.0, *S, 1.0, *Js);
            K_.Reset(Sum);
            return *K_.Ptr();
         }
      }

      // Fallback: return source Jacobian only
      return Jsrc;
   }

   /** @brief Project d to satisfy irreversibility: d >= d_old, and clamp to [0,1].
    *
    * @param[in,out] d Damage vector to modify
    * @param[in] d_old Previous damage vector
    */
   void ApplyIrreversibility(Vector &d, const Vector &d_old) const
   {
      MFEM_VERIFY(d.Size() == d_old.Size(), "size mismatch");
      for (int i = 0; i < d.Size(); i++)
      {
         real_t v = std::max(d(i), d_old(i));
         if (v < 0.0) { v = 0.0; }
         if (v > 1.0) { v = 1.0; }
         d(i) = v;
      }
   }

   /// Get access to the strain energy history field
   ParGridFunction &GetStrainEnergyHistory() { return H_gf_; }
   const ParGridFunction &GetStrainEnergyHistory() const { return H_gf_; }

private:
   ParFiniteElementSpace &fes_;
   const PFFMaterialParameters &mat_;
   mutable ParGridFunction H_gf_;      // history field in H1 for simplicity
   mutable ParNonlinearForm *nlf_;
   mutable OperatorPtr diffusion_op_;  // assembled diffusion operator
   mutable OperatorPtr K_;             // assembled total gradient when requested
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_DAMAGE_OPERATOR_HPP
