// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_wave_operator.hpp — Nitsche dynamic rupture, plan item N-01.
//
// P3 content: the fault-free explicit elastodynamic operator
//
//    M a = F(t) - K u,    a = d2u/dt2,
//
// P4: optional first-order absorbing (dashpot) boundary on
// BoundaryConfig::absorbing_attrs, adding a velocity-only damping C:
//
//    M a + C v + K u = F(t).
//
// With central difference (gamma = 1/2) the step solves
// (M + dt/2 C) A^{n+1} = F(t_{n+1}) - K U^{n+1} - C U*  (U* = predicted
// velocity), so that M A^{n+1} + C Udot^{n+1} = F - K U^{n+1} holds exactly
// with Udot^{n+1} = U* + dt/2 A^{n+1}.  Without absorbing faces nothing
// changes (bit-identical to P3).
//
// as an MFEM SecondOrderTimeDependentOperator, to be advanced by MFEM's
// CentralDifferenceSolver (Newmark beta = 0, gamma = 1/2; see
// solver/nitsche_cfl.hpp).  K is the IP stiffness of ElasticityDomainOperator
// built through the Coefficient constructor (C-09) from the per-element
// NitscheMaterialSample; M^{-1} is NitscheMass (C-02) from the same sample.
// The fault (SetFaultState, slip, ...) is not part of P3.
//
// NitscheMaterialSample lives in nitsche_material_sample.hpp (moved there in
// P3 to break an include cycle with nitsche_mass.hpp); it is re-exported here.
//
// Scope (Chunhui Zhao, 2026-10-13): the Nitsche path uses tetrahedral meshes only.

#ifndef MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP
#define MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP

#include "mfem.hpp"
#include "nitsche_material_sample.hpp"
#include "nitsche_mass.hpp"
#include "elasticity_operator.hpp"
#include "../integrator/dg_absorbing_bc_integrator.hpp"

#include <cmath>
#include <limits>
#include <memory>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// Fault-free explicit elastodynamics: d2u/dt2 = M^{-1} (F(t) - K u).
class NitscheWaveOperator : public SecondOrderTimeDependentOperator
{
public:
   /// `material` is sampled once per element; K and M share that sample.
   /// `bc` must carry no fault faces in P3.  Dirichlet data: F(t) comes from
   /// ElasticityDomainOperator::AssembleDirichletOnlyRHS with Vp = 0, i.e.
   /// F = 0 (its BP5 shortcut returns early when time * Vp == 0).
   NitscheWaveOperator(ParMesh &mesh, int order, const MaterialField &material,
                       const BoundaryConfig &bc,
                       const DomainConfig &config = {})
      : SecondOrderTimeDependentOperator(0, 0.0),
        mesh_(mesh), order_(order)
   {
      sample_ = std::make_unique<NitscheMaterialSample>(mesh, material);

      // Scalar fallbacks for the operator's scalar-only accessors: global
      // element means of the sample (not used by the IP assembly).
      real_t sums[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < sample_->LambdaGF().Size(); ++i)
      {
         sums[0] += sample_->LambdaGF()(i);
         sums[1] += sample_->MuGF()(i);
      }
      sums[2] = mesh.GetNE();
      MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPITypeMap<real_t>::mpi_type,
                    MPI_SUM, mesh.GetComm());
      const real_t lam_fb = sums[0] / sums[2];
      const real_t mu_fb  = sums[1] / sums[2];

      elast_ = std::make_unique<ElasticityDomainOperator<ParMesh>>(
         mesh, order, sample_->Lambda(), sample_->Mu(), lam_fb, mu_fb,
         /*Vp=*/0.0, /*Wf=*/0.0, /*lf=*/0.0, bc, DGMethod::IP,
         SolverType::CG_AMG, config);
      MFEM_VERIFY(elast_->GetNumFaultDOFs() == 0,
                  "NitscheWaveOperator (P3): the mesh must carry no fault "
                  "faces; found " << elast_->GetNumFaultDOFs()
                  << " fault DOFs.");

      RhoPerElement(*sample_, rho_e_);
      mass_ = std::make_unique<NitscheMass>(mesh, order, rho_e_);
      BuildDampingBlocks_();

      const int n = elast_->GetFESpace().GetVSize();
      MFEM_VERIFY(mass_->VSize() == n,
                  "NitscheWaveOperator: mass size " << mass_->VSize()
                  << " != stiffness vector size " << n);
      height = width = n;
      Kx_.SetSize(n);
      F_.SetSize(n);
   }

   /// Time step used to check the Newmark coefficients in ImplicitSolve.
   /// With absorbing faces it also builds (M + dt/2 C)^{-1}.
   void SetTimeStep(real_t dt)
   {
      dt_ = dt;
      if (!damping_blocks_.empty())
      {
         mass_damped_ = std::make_unique<NitscheMass>(
            mesh_, order_, rho_e_, &damping_blocks_, 0.5 * dt);
      }
   }
   real_t GetTimeStep() const { return dt_; }

   /// First-order form (Operator::Mult) is not used: second-order solvers
   /// call the three-argument Mult and ImplicitSolve below.
   void Mult(const Vector &, Vector &) const override
   {
      MFEM_ABORT("NitscheWaveOperator is second order; use "
                 "Mult(x, dxdt, y) / ImplicitSolve with a "
                 "SecondOrderODESolver.");
   }

   /// y = M^{-1}(F(t) - K x), t = GetTime().  Called once by NewmarkSolver
   /// for the initial acceleration A^0; advances nothing.
   void Mult(const Vector &x, const Vector &dxdt, Vector &y) const override
   {
      Acceleration_(x, &dxdt, *mass_, y);
   }

   /// Central difference (beta = 0, gamma = 1/2): fac0 = beta dt^2 = 0 and
   /// fac1 = gamma dt = dt/2, so the "implicit" solve is explicit:
   /// k = M^{-1}(F(t) - K x) with t = GetTime() = t_{n+1}.  One call per step
   /// (the time stamp below; P6 hangs the fault advance on it).
   void ImplicitSolve(const real_t fac0, const real_t fac1,
                      const Vector &x, const Vector &dxdt, Vector &k) override
   {
      MFEM_VERIFY(fac0 == 0.0,
                  "NitscheWaveOperator::ImplicitSolve: fac0 = " << fac0
                  << " != 0; only the central-difference scheme "
                  "(Newmark beta = 0) is supported.");
      MFEM_VERIFY(fac1 == 0.5 * dt_,
                  "NitscheWaveOperator::ImplicitSolve: fac1 = " << fac1
                  << " != dt/2 = " << 0.5 * dt_
                  << "; call SetTimeStep(dt) with the solver's dt "
                  "(Newmark gamma = 1/2).");
      const real_t t = GetTime();
      MFEM_VERIFY(t != last_solve_time_,
                  "NitscheWaveOperator::ImplicitSolve called twice at t = "
                  << t << "; exactly one call per time step is allowed.");
      last_solve_time_ = t;
      if (damping_blocks_.empty())
      {
         Acceleration_(x, nullptr, *mass_, k);
      }
      else
      {
         MFEM_VERIFY(mass_damped_ && mass_damped_->HalfDt() == fac1,
                     "NitscheWaveOperator::ImplicitSolve: the damped mass "
                     "(M + dt/2 C) was built for dt/2 = "
                     << (mass_damped_ ? mass_damped_->HalfDt() : 0.0)
                     << " but fac1 = " << fac1 << "; call SetTimeStep(dt).");
         Acceleration_(x, &dxdt, *mass_damped_, k);
      }
   }

   const NitscheMaterialSample &Sample() const { return *sample_; }
   ParMesh &GetMesh() const { return mesh_; }
   int GetOrder() const { return order_; }
   ElasticityDomainOperator<ParMesh> &Elasticity() { return *elast_; }
   const NitscheMass &Mass() const { return *mass_; }
   /// (M + dt/2 C)^{-1}; null without absorbing faces or before SetTimeStep.
   const NitscheMass *DampedMass() const { return mass_damped_.get(); }
   /// Per-element dashpot blocks C_e (element-local byNODES, 3 nd x 3 nd).
   const std::map<int, DenseMatrix> &DampingBlocks() const { return damping_blocks_; }
   bool HasDamping() const { return !damping_blocks_.empty(); }

   /// Cv = C v (zero without absorbing faces).
   void ApplyDamping(const Vector &v, Vector &Cv) const
   {
      Cv.SetSize(v.Size());
      Cv = 0.0;
      if (damping_blocks_.empty()) { return; }
      const int n = mass_->NumScalarDofs();
      Vector x, y;
      for (const auto &kv : damping_blocks_)
      {
         const Array<int> &dofs = mass_->ElementScalarDofs(kv.first);
         const int nd = dofs.Size();
         x.SetSize(3 * nd);
         y.SetSize(3 * nd);
         for (int c = 0; c < 3; ++c)
            for (int k = 0; k < nd; ++k) { x(c * nd + k) = v(c * n + dofs[k]); }
         kv.second.Mult(x, y);
         for (int c = 0; c < 3; ++c)
            for (int k = 0; k < nd; ++k) { Cv(c * n + dofs[k]) += y(c * nd + k); }
      }
   }

private:
   /// Shared by Mult and ImplicitSolve so both give identical bits:
   /// a = Minv (F(t) - K x - C v) with t = GetTime().  Without absorbing
   /// faces the C v term is skipped entirely (P3 bits).
   void Acceleration_(const Vector &x, const Vector *v,
                      const NitscheMass &Minv, Vector &a) const
   {
      elast_->ApplyStiffness(x, Kx_);
      elast_->AssembleDirichletOnlyRHS(F_, GetTime());
      F_ -= Kx_;
      if (!damping_blocks_.empty() && v)
      {
         ApplyDamping(*v, Cv_);
         F_ -= Cv_;
      }
      Minv.Apply(F_, a);
   }

   /// C_e for every element with an absorbing face (N-06), material from the
   /// sample (not the operator's fallback scalars).
   void BuildDampingBlocks_()
   {
      damping_blocks_.clear();
      const auto &fes = elast_->GetFESpace();
      const ParFiniteElementSpace *mfes = sample_->RhoGF().ParFESpace();
      Array<int> mdofs;
      for (int be : elast_->GetAbsorbingBdrElements())
      {
         FaceElementTransformations *Tr = mesh_.GetBdrFaceTransformations(be);
         MFEM_VERIFY(Tr != nullptr, "NitscheWaveOperator: no face "
                     "transformation for absorbing boundary element " << be);
         const int e = Tr->Elem1No;
         const FiniteElement &fe = *fes.GetFE(e);
         mfes->GetElementDofs(e, mdofs);
         const real_t lam = sample_->LambdaGF()(mdofs[0]);
         const real_t mu  = sample_->MuGF()(mdofs[0]);
         const real_t rho = sample_->RhoGF()(mdofs[0]);
         const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
         const real_t cs = std::sqrt(mu / rho);
         const int nd = fe.GetDof();
         auto it = damping_blocks_.find(e);
         if (it == damping_blocks_.end())
         {
            it = damping_blocks_.emplace(e, DenseMatrix(3 * nd)).first;
            it->second = 0.0;
         }
         DGAbsorbingBCIntegrator::AddFaceBlock(fe, *Tr, rho, cp, cs, it->second);
      }
   }

   ParMesh &mesh_;
   int order_;
   std::unique_ptr<NitscheMaterialSample> sample_;
   std::unique_ptr<ElasticityDomainOperator<ParMesh>> elast_;
   std::unique_ptr<NitscheMass> mass_;
   std::unique_ptr<NitscheMass> mass_damped_;
   std::map<int, DenseMatrix> damping_blocks_;
   Vector rho_e_;
   real_t dt_ = 0.0;
   real_t last_solve_time_ = std::numeric_limits<real_t>::quiet_NaN();
   mutable Vector Kx_, F_, Cv_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP
