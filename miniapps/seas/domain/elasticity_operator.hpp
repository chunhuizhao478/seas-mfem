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

#ifndef MFEM_SEAS_ELASTICITY_OPERATOR_HPP
#define MFEM_SEAS_ELASTICITY_OPERATOR_HPP

#include "mfem.hpp"
#include "domain_operator.hpp"
#include "antiplane_operator.hpp"  // For DGMethod enum
#include "../common/seas_types.hpp"
#include "../fault/fault_basis.hpp"
#include "../integrator/dg_elasticity_br2_integrator.hpp"

#include <memory>
#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

/// @brief DG Elasticity domain operator for 3D vector elasticity (BP5)
///
/// Solves the 3D linear elasticity problem:
///   ∇·σ(u) = 0, σ = λ·tr(ε)·I + 2μ·ε
///
/// using DG method with fault slip as interior jump BC.
///
/// Coordinate convention (SCEC BP5):
///   x1 = fault-normal
///   x2 = along-strike
///   x3 = depth (positive downward in SCEC, but mapped to mesh z-axis)
///
/// - Fault at x1=0 is an interior interface
/// - Slip imposed as jump [[u]] on fault interior faces
/// - Dirichlet loading on ±x2 walls: u₂ = ±Vp·t/2
/// - Free surface / natural BC on other boundaries
///
/// Supports both BR2 (default, matching Tandem) and IP DG methods.
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class ElasticityDomainOperator : public DomainOperator<MeshType>
{
public:
   using Base = DomainOperator<MeshType>;
   using FESpaceType = typename Base::FESpaceType;
   using GridFuncType = typename Base::GridFuncType;
   using BilinFormType = typename Base::BilinFormType;
   using LinFormType = typename Base::LinFormType;

   /// @brief Construct the 3D DG elasticity operator
   ///
   /// @param mesh The 3D mesh (hex or tet)
   /// @param order Polynomial order for DG
   /// @param lambda First Lame parameter [Pa]
   /// @param mu Shear modulus [Pa]
   /// @param Vp Plate rate [m/s]
   /// @param Wf Fault depth [m] (rate-state zone)
   /// @param lf Fault length [m] (along-strike extent)
   /// @param method DG method (BR2 default, IP alternative)
   ElasticityDomainOperator(MeshType &mesh, int order,
                             real_t lambda, real_t mu,
                             real_t Vp, real_t Wf, real_t lf,
                             DGMethod method = DGMethod::BR2)
      : mesh_(mesh), order_(order),
        lambda_val_(lambda), mu_val_(mu),
        Vp_(Vp), Wf_(Wf), lf_(lf),
        method_(method),
        lambda_coeff_(lambda), mu_coeff_(mu),
        mass_inv_computed_(false),
        fault_depths_computed_(false),
        fault_coords_computed_(false)
   {
      MFEM_VERIFY(mesh_.Dimension() == 3, "ElasticityDomainOperator requires 3D mesh");

      epsilon_ = -1.0;  // SIPG

      SetupFESpace();
      SetupBoundaryMarkers();
      SetupFaultInfo();
      SetupSolver();

      // Precompute fault depths/coordinates eagerly
      {
         Vector tmp;
         GetFaultDepths(tmp);
      }

      if (method_ == DGMethod::BR2)
      {
         PrecomputeMassInverse();
      }
   }

   ~ElasticityDomainOperator() override = default;

   int NumComponents() const override { return 3; }
   int Dimension() const override { return 3; }
   int NumSlipComponents() const override { return 2; }

   void Solve(real_t time, const Vector &slip_bc,
              GridFuncType &displacement) override;

   void ComputeTraction(const GridFuncType &displacement,
                        const Vector &slip_bc,
                        Vector &traction) override;

   FESpaceType &GetFESpace() override { return *fes_; }
   const FESpaceType &GetFESpace() const override { return *fes_; }

   MeshType &GetMesh() override { return mesh_; }
   const MeshType &GetMesh() const override { return mesh_; }

   real_t GetShearModulus() const override { return mu_val_; }

   int GetNumFaultDOFs() const override { return num_fault_dofs_; }

   void GetFaultDepths(Vector &depths) const override;

   void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const override;

   const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

   const FaultBasis *GetFaultBasis() const override { return &fault_basis_; }

   const Array<int> &GetFaultInteriorFaces() const { return fault_interior_faces_; }
   const Array<int> &GetFaultSharedFaces() const { return fault_shared_faces_; }

   DGMethod GetMethod() const { return method_; }
   int GetOrder() const { return order_; }
   real_t GetPlateRate() const { return Vp_; }
   real_t GetFaultDepthLimit() const { return Wf_; }
   real_t GetFaultLength() const { return lf_; }
   real_t GetLambda() const { return lambda_val_; }

private:
   MeshType &mesh_;
   int order_;
   real_t lambda_val_, mu_val_;
   real_t Vp_, Wf_, lf_;
   DGMethod method_;
   real_t epsilon_;  // SIPG sign = -1

   // Coefficients (mutable: used in const assembly methods, MFEM Coefficient::Eval is non-const)
   mutable ConstantCoefficient lambda_coeff_, mu_coeff_;

   // FE spaces
   std::unique_ptr<DG_FECollection> fec_;
   std::unique_ptr<FESpaceType> fes_;            // vector DG (vdim=3)
   std::unique_ptr<FiniteElementSpace> scalar_fes_;  // scalar DG (for BR2 Minv)

   // BR2 mass matrix inverses
   mutable std::vector<DenseMatrix> elem_mass_inv_;
   mutable bool mass_inv_computed_;

   // Stiffness matrix
   mutable bool stiffness_assembled_ = false;
   mutable std::unique_ptr<BilinFormType> cached_a_;
   mutable OperatorHandle cached_Ah_;
   mutable std::unique_ptr<Solver> cached_prec_;

   // Solver
   mutable std::unique_ptr<Solver> solver_;
   mutable Vector X_, B_;
   mutable std::unique_ptr<GSSmoother> serial_prec_;

   // Fault data
   Array<int> fault_interior_faces_;
   Array<int> fault_shared_faces_;
   Array<int> fault_dofs_;
   int num_fault_dofs_ = 0;
   FaultBasis fault_basis_;

   mutable Vector fault_depths_;
   mutable bool fault_depths_computed_;
   mutable Vector fault_x2_, fault_x3_;
   mutable bool fault_coords_computed_;

   // Boundary markers
   mutable Array<int> dirichlet_bdr_marker_;

   // ========================================================================
   // Setup methods
   // ========================================================================

   void SetupFESpace()
   {
      fec_ = std::make_unique<DG_FECollection>(order_, 3, BasisType::GaussLobatto);
      // Vector DG space with 3 components, byNODES ordering
      fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get(), 3, Ordering::byNODES);
      // Scalar DG space for BR2 mass matrix inverses
      scalar_fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
   }

   void SetupBoundaryMarkers()
   {
      // Identify Dirichlet boundaries (loading walls at ±x1, fault-normal)
      // We mark boundaries that will receive Dirichlet loading
      int num_bdr = mesh_.bdr_attributes.Size() > 0 ? mesh_.bdr_attributes.Max() : 0;
      dirichlet_bdr_marker_.SetSize(num_bdr);
      dirichlet_bdr_marker_ = 0;

      // Mark Dirichlet boundaries by attribute
      // Convention: attributes 1,2 are ±x1 walls (fault-normal, Dirichlet loading)
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr == 1 || attr == 2)
         {
            dirichlet_bdr_marker_[attr - 1] = 1;
         }
      }
   }

   void SetupFaultInfo()
   {
      fault_interior_faces_.SetSize(0);

      // Detect fault faces at x1 ≈ 0
      int num_faces = mesh_.GetNumFaces();
      for (int f = 0; f < num_faces; f++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (FTr == nullptr) { continue; }

         if (IsFaultFace3D(FTr))
         {
            fault_interior_faces_.Append(f);
         }
      }

      // Shared faces (parallel): coordinate-based detection.
      fault_shared_faces_.SetSize(0);
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.ExchangeFaceNbrData();
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (pfes) { pfes->ExchangeFaceNbrData(); }

         for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
         {
            if (IsFaultFace3DShared(sf))
            {
               fault_shared_faces_.Append(sf);
            }
         }
#endif
      }

      num_fault_dofs_ = fault_interior_faces_.Size() + fault_shared_faces_.Size();

      fault_dofs_.SetSize(num_fault_dofs_);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         fault_dofs_[i] = i;
      }

      // Compute FaultBasis for coordinate transforms
      if (num_fault_dofs_ > 0)
      {
         Vector ref_normal(3);
         ref_normal = 0.0;
         ref_normal(0) = 1.0;  // x1 is fault-normal

         Vector up(3);
         up = 0.0;
         up(2) = 1.0;  // Up direction (matches Tandem convention)

         fault_basis_.Compute(mesh_, fault_interior_faces_, ref_normal, up);

         // Append shared faces to FaultBasis
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            fault_basis_.AppendSharedFaces(mesh_, fault_shared_faces_,
                                           ref_normal, up);
#endif
         }
      }
   }

   bool IsFaultFace3D(FaceElementTransformations *FTr) const
   {
      const IntegrationPoint &ip = Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector center(3);
      FTr->Face->Transform(ip, center);

      // Fault is at x1 = 0, within the fault zone bounds
      const real_t tol = 1e-10 * std::max(Wf_, 1.0);
      return std::abs(center(0)) < tol
          && std::abs(center(1)) <= lf_ / 2.0 + tol
          && center(2) >= -tol
          && center(2) <= Wf_ + tol;
   }

   bool IsFaultFace3DShared(int shared_face) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(shared_face);
         if (FTr == nullptr) { return false; }

         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->Face->SetIntPoint(&ip);
         Vector center(3);
         FTr->Face->Transform(ip, center);

         const real_t tol = 1e-10 * std::max(Wf_, 1.0);
         return std::abs(center(0)) < tol
             && std::abs(center(1)) <= lf_ / 2.0 + tol
             && center(2) >= -tol
             && center(2) <= Wf_ + tol;
#endif
      }
      return false;
   }

   void SetupSolver()
   {
      X_.SetSize(fes_->GetTrueVSize());
      B_.SetSize(fes_->GetTrueVSize());

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         // Parallel solver created in AssembleStiffness
      }
      else
      {
         auto *cg = new CGSolver();
         cg->SetRelTol(1e-12);
         cg->SetAbsTol(0.0);
         cg->SetMaxIter(10000);
         cg->SetPrintLevel(0);
         solver_.reset(cg);
      }
   }

   // ========================================================================
   // Mass inverse precomputation (for BR2)
   // ========================================================================

   void PrecomputeMassInverse() const
   {
      if (mass_inv_computed_) { return; }

      int ne = mesh_.GetNE();
      elem_mass_inv_.resize(ne);

      for (int e = 0; e < ne; e++)
      {
         // Use scalar FE space for mass matrix (same for all components)
         const FiniteElement *fe = scalar_fes_->GetFE(e);
         ElementTransformation *T = mesh_.GetElementTransformation(e);

         int ndof = fe->GetDof();
         DenseMatrix M(ndof);

         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                   2 * fe->GetOrder());
         Vector shape(ndof);
         M = 0.0;

         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            real_t w = ip.weight * T->Weight();
            for (int k = 0; k < ndof; k++)
            {
               for (int l = 0; l < ndof; l++)
               {
                  M(k, l) += w * shape(k) * shape(l);
               }
            }
         }

         elem_mass_inv_[e].SetSize(ndof);
         DenseMatrixInverse M_inv(M);
         M_inv.GetInverseMatrix(elem_mass_inv_[e]);
      }

#ifdef MFEM_USE_MPI
      // In parallel, also compute mass inverses for face-neighbor elements.
      // These are needed by DGElasticityBR2Integrator::AssembleFaceMatrix when
      // ParBilinearForm::AssembleSharedFaces passes a shared face with
      // Trans.Elem2No >= mesh.GetNE().
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         MFEM_VERIFY(pfes, "Expected ParFiniteElementSpace in parallel");
         pfes->ExchangeFaceNbrData();

         ParMesh *pmesh = pfes->GetParMesh();
         int nel_nbr = pmesh->GetNFaceNeighborElements();
         elem_mass_inv_.resize(ne + nel_nbr);

         for (int i = 0; i < nel_nbr; i++)
         {
            const FiniteElement *fe = pfes->GetFaceNbrFE(i);
            ElementTransformation *T =
               pmesh->GetFaceNbrElementTransformation(i);

            int ndof = fe->GetDof();
            DenseMatrix M(ndof);

            const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                      2 * fe->GetOrder());
            Vector shape(ndof);
            M = 0.0;

            for (int j = 0; j < ir.GetNPoints(); j++)
            {
               const IntegrationPoint &ip = ir.IntPoint(j);
               T->SetIntPoint(&ip);
               fe->CalcShape(ip, shape);

               real_t w = ip.weight * T->Weight();
               for (int k = 0; k < ndof; k++)
               {
                  for (int l = 0; l < ndof; l++)
                  {
                     M(k, l) += w * shape(k) * shape(l);
                  }
               }
            }

            elem_mass_inv_[ne + i].SetSize(ndof);
            DenseMatrixInverse M_inv(M);
            M_inv.GetInverseMatrix(elem_mass_inv_[ne + i]);
         }
      }
#endif

      mass_inv_computed_ = true;
   }

   // ========================================================================
   // Stiffness assembly
   // ========================================================================

   void AssembleStiffness() const
   {
      // In parallel, exchange face neighbor data before assembly.
      // ParBilinearForm::AssembleSharedFaces needs this to access
      // face neighbor elements and DOFs.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.ExchangeFaceNbrData();
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (pfes) { pfes->ExchangeFaceNbrData(); }
#endif
      }

      cached_a_ = std::make_unique<BilinFormType>(fes_.get());

      // Volume term: ∫_Ω σ(u):ε(v) dV
      cached_a_->AddDomainIntegrator(
         new ElasticityIntegrator(lambda_coeff_, mu_coeff_));

      if (method_ == DGMethod::BR2)
      {
         if (!mass_inv_computed_)
         {
            PrecomputeMassInverse();
         }
         // BR2 interior face integrator with elasticity tensor coupling
         cached_a_->AddInteriorFaceIntegrator(
            new DGElasticityBR2Integrator(
               lambda_coeff_, mu_coeff_, epsilon_, elem_mass_inv_, 3));

         // BR2 boundary face integrator for Dirichlet walls
         if (dirichlet_bdr_marker_.Size() > 0)
         {
            cached_a_->AddBdrFaceIntegrator(
               new DGElasticityBR2BoundaryIntegrator(
                  lambda_coeff_, mu_coeff_, epsilon_, elem_mass_inv_, 3),
               dirichlet_bdr_marker_);
         }
      }
      else  // IP
      {
         real_t kappa = (order_ + 1) * (order_ + 1);
         cached_a_->AddInteriorFaceIntegrator(
            new DGElasticityIntegrator(lambda_coeff_, mu_coeff_, epsilon_, kappa));

         if (dirichlet_bdr_marker_.Size() > 0)
         {
            cached_a_->AddBdrFaceIntegrator(
               new DGElasticityIntegrator(lambda_coeff_, mu_coeff_, epsilon_, kappa),
               dirichlet_bdr_marker_);
         }
      }

      cached_a_->Assemble();
      cached_a_->Finalize();

      // Set up solver operator
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         cached_Ah_.SetType(Operator::Hypre_ParCSR);
         cached_a_->ParallelAssemble(cached_Ah_);

#ifdef MFEM_USE_MUMPS
         auto *mumps = new MUMPSSolver(mesh_.GetComm());
         mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
         mumps->SetPrintLevel(0);
         mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
         solver_.reset(mumps);
#else
         auto *cg = new CGSolver(mesh_.GetComm());
         cg->SetRelTol(1e-12);
         cg->SetAbsTol(0.0);
         cg->SetMaxIter(10000);
         cg->SetPrintLevel(0);

         auto *ilu = new HypreILU();
         ilu->SetLevelOfFill(1);
         ilu->SetPrintLevel(0);
         ilu->SetOperator(*cached_Ah_.As<HypreParMatrix>());
         cached_prec_.reset(ilu);

         cg->SetPreconditioner(*cached_prec_);
         cg->SetOperator(*cached_Ah_.As<HypreParMatrix>());
         solver_.reset(cg);
#endif
#endif
      }
      else
      {
         SparseMatrix &A = cached_a_->SpMat();
         serial_prec_ = std::make_unique<GSSmoother>(A);
         static_cast<CGSolver*>(solver_.get())->SetPreconditioner(*serial_prec_);
         solver_->SetOperator(A);
      }

      stiffness_assembled_ = true;
   }

   // ========================================================================
   // Slip BC assembly — IP method
   // ========================================================================

   void AssembleSlipContributionIP(Vector &rhs, const Vector &slip_bc) const
   {
      int dim = 3;
      real_t kappa = (order_ + 1) * (order_ + 1);

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         // Embed local slip to global frame
         real_t slip_local[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(fi, slip_local, delta_u);

         real_t max_delta = std::max({std::abs(delta_u[0]),
                                       std::abs(delta_u[1]),
                                       std::abs(delta_u[2])});
         if (max_delta < 1e-15) { continue; }

         // Get DOFs (vector space)
         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         int ndof1 = fe1->GetDof();
         int ndof2 = fe2->GetDof();

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2 * face_order + 1);

         Vector elvec1(vdofs1.Size()), elvec2(vdofs2.Size());
         elvec1 = 0.0;
         elvec2 = 0.0;

         for (int p = 0; p < ir.GetNPoints(); p++)
         {
            const IntegrationPoint &ip = ir.IntPoint(p);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            // Face normal
            Vector nor(dim);
            CalcOrtho(FTr->Jacobian(), nor);

            // Determine slip sign based on normal orientation
            // Normal points from elem1 to elem2; slip jump convention:
            // [[u]] = u+ - u- = delta_u (from side with n pointing outward)
            real_t sign = (nor(0) > 0) ? -1.0 : 1.0;

            // Shapes
            Vector shape1(ndof1), shape2(ndof2);
            fe1->CalcShape(eip1, shape1);
            fe2->CalcShape(eip2, shape2);

            // Physical gradients
            DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
            fe1->CalcDShape(eip1, dshape1_ref);
            fe2->CalcDShape(eip2, dshape2_ref);

            DenseMatrix adjJ1(dim), adjJ2(dim);
            CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);
            CalcAdjugate(FTr->Elem2->Jacobian(), adjJ2);

            DenseMatrix dshape1_adj(ndof1, dim), dshape2_adj(ndof2, dim);
            Mult(dshape1_ref, adjJ1, dshape1_adj);
            Mult(dshape2_ref, adjJ2, dshape2_adj);

            real_t detJ1 = FTr->Elem1->Weight();
            real_t detJ2 = FTr->Elem2->Weight();
            real_t w1 = ip.weight / (2.0 * detJ1);
            real_t w2 = ip.weight / (2.0 * detJ2);
            real_t nor_sq = nor * nor;
            real_t wq_penalty = kappa * nor_sq * (w1 + w2);

            // For each element, each DOF k, each component i:
            // Symmetry term: σ * [σ(φ_k e_i)·n]_u * delta_u[u]
            // Penalty term: kappa * φ_k * delta_u[i]
            for (int k = 0; k < ndof1; k++)
            {
               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  // [σ(φ_k e_i)·n]_u = λ*(∂φ_k/∂x_i)*n_u + μ*(δ_iu*∇φ_k·n + ∂φ_k/∂x_u*n_i)
                  Vector dn1(1);
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape1_adj(k, d) * nor(d);
                  }
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape1_adj(k, i) * nor(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape1_adj(k, u) * nor(i));
                     sym_val += trac_iu * sign * delta_u[u];
                  }

                  // byNODES ordering: DOF k, component i
                  int idx = i * ndof1 + k;
                  elvec1(idx) += epsilon_ * sym_val * w1;
                  elvec1(idx) += wq_penalty * sign * delta_u[i] * shape1(k);
               }
            }

            for (int k = 0; k < ndof2; k++)
            {
               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape2_adj(k, d) * nor(d);
                  }
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape2_adj(k, i) * nor(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape2_adj(k, u) * nor(i));
                     sym_val += trac_iu * sign * delta_u[u];
                  }

                  int idx = i * ndof2 + k;
                  elvec2(idx) += epsilon_ * sym_val * w2;
                  elvec2(idx) -= wq_penalty * sign * delta_u[i] * shape2(k);
               }
            }
         }

         // Add to global RHS
         for (int j = 0; j < vdofs1.Size(); j++)
         {
            int gj = vdofs1[j];
            if (gj >= 0) { rhs(gj) += elvec1(j); }
            else { rhs(-1 - gj) -= elvec1(j); }
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            int gj = vdofs2[j];
            if (gj >= 0) { rhs(gj) += elvec2(j); }
            else { rhs(-1 - gj) -= elvec2(j); }
         }
      }
   }

   // ========================================================================
   // Slip BC assembly — BR2 method
   // ========================================================================

   void AssembleSlipContributionBR2(Vector &rhs, const Vector &slip_bc) const
   {
      int dim = 3;
      if (!mass_inv_computed_) { PrecomputeMassInverse(); }

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         // BR2 penalty: num faces per element (geometry-dependent)
         Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
         real_t penalty = (geom == Geometry::TETRAHEDRON)
                              ? real_t(dim + 1) : real_t(2 * dim);

         // Embed local slip to global frame
         real_t slip_local[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(fi, slip_local, delta_u);

         real_t max_delta = std::max({std::abs(delta_u[0]),
                                       std::abs(delta_u[1]),
                                       std::abs(delta_u[2])});
         if (max_delta < 1e-15) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         int ndof1 = fe1->GetDof();
         int ndof2 = fe2->GetDof();

         const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
         const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2 * face_order + 1);
         int nqp = ir.GetNPoints();

         // Precompute shapes and normals
         DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
         Vector nor_all(dim * nqp);
         Vector w_all(nqp);
         Vector slip_sign_all(nqp);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            Vector s1q(shapes1.GetColumn(q), ndof1);
            Vector s2q(shapes2.GetColumn(q), ndof2);
            fe1->CalcShape(eip1, s1q);
            fe2->CalcShape(eip2, s2q);

            Vector nor_q(&nor_all[q * dim], dim);
            CalcOrtho(FTr->Jacobian(), nor_q);

            w_all[q] = ip.weight;
            slip_sign_all[q] = (nor_q(0) > 0) ? -1.0 : 1.0;
         }

         // ===================================================================
         // Compute BR2 lifted slip with elasticity tensor coupling
         //
         // Step 1: face_int[elem][m,s] = sum_q E[elem][m,q] * delta_u_signed[s,q] * w[q]
         //   where delta_u_signed[s,q] = sum_c slip_sign * delta_u[c] * ... no, simpler:
         //   Actually for the RHS, the "source" is the imposed slip delta_u.
         //   The lifting formula is:
         //     f_lifted[elem][u,m,s] = 0.5 * Minv[m,o] * sum_q E[o,q] * delta_u_signed[u] * n[s,q] * w[q]
         //   Then f_lifted_q[i,q] = 0.5 * sum_{elem,s,m} test_normal_{iu,sq} * E[elem][m,q] * f_lifted[elem][u,m,s]
         //
         // This follows exactly the scalar BR2 pattern from antiplane but with
         // delta_u being 3-component and test_normal coupling.
         // ===================================================================

         // Compute face integrals for each element
         // face_int[elem][u*dim+s, m] = sum_q E[m,q] * sign_q * delta_u[u] * n[s,q] * w[q]
         DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
         face_int1 = 0.0;
         face_int2 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            real_t wq = w_all[q];
            real_t sign_q = slip_sign_all[q];

            for (int u = 0; u < dim; u++)
            {
               real_t delta_u_signed = sign_q * delta_u[u];
               for (int s = 0; s < dim; s++)
               {
                  real_t n_s = nor_all[q * dim + s];
                  real_t factor = delta_u_signed * n_s * wq;

                  for (int m = 0; m < ndof1; m++)
                  {
                     face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                  }
                  for (int m = 0; m < ndof2; m++)
                  {
                     face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                  }
               }
            }
         }

         // f_lifted[elem][u,s,m] = 0.5 * Minv[m,o] * face_int[elem][u*dim+s, o]
         DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
         MultABt(face_int1, Minv1, f_lifted1);
         f_lifted1 *= 0.5;
         MultABt(face_int2, Minv2, f_lifted2);
         f_lifted2 *= 0.5;

         // f_lifted_q[i,q] = 0.5 * sum_{elem,u,s,m} test_normal_{iu,s} * E[elem][m,q] * f_lifted[elem][u,s,m]
         DenseMatrix f_lifted_q(dim, nqp);
         f_lifted_q = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            Vector n_q(dim);
            for (int d = 0; d < dim; d++) { n_q(d) = nor_all[q * dim + d]; }

            for (int i = 0; i < dim; i++)
            {
               real_t sum = 0.0;
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                     + (i == s ? 1.0 : 0.0) * n_q(u));

                     // Eval lifted: sum_m E[elem][m,q] * f_lifted[elem][u*dim+s, m]
                     real_t eval1 = 0.0, eval2 = 0.0;
                     for (int m = 0; m < ndof1; m++)
                     {
                        eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                     }
                     for (int m = 0; m < ndof2; m++)
                     {
                        eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                     }
                     sum += tn * (eval1 + eval2);
                  }
               }
               f_lifted_q(i, q) = 0.5 * sum;
            }
         }

         // ===================================================================
         // Assemble RHS:
         // b_k^i += c1 * sum_u [σ_i(φ_k e_i)·n]_u * delta_u_signed[u] * w_elem
         //        ± penalty * φ_k * f_lifted_q[i,q] * w
         // ===================================================================

         Vector elvec1(vdofs1.Size()), elvec2(vdofs2.Size());
         elvec1 = 0.0;
         elvec2 = 0.0;

         real_t c1 = epsilon_ * 0.5;  // -0.5 for SIPG

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            real_t wq = w_all[q];
            real_t sign_q = slip_sign_all[q];

            Vector nor_q(dim);
            for (int d = 0; d < dim; d++) { nor_q[d] = nor_all[q * dim + d]; }

            // Gradients
            DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
            fe1->CalcDShape(eip1, dshape1_ref);
            fe2->CalcDShape(eip2, dshape2_ref);

            DenseMatrix adjJ1(dim), adjJ2(dim);
            CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);
            CalcAdjugate(FTr->Elem2->Jacobian(), adjJ2);

            DenseMatrix dshape1_adj(ndof1, dim), dshape2_adj(ndof2, dim);
            Mult(dshape1_ref, adjJ1, dshape1_adj);
            Mult(dshape2_ref, adjJ2, dshape2_adj);

            real_t detJ1 = FTr->Elem1->Weight();
            real_t detJ2 = FTr->Elem2->Weight();

            // Element 1 consistency + BR2
            for (int k = 0; k < ndof1; k++)
            {
               real_t grad_dot_n = 0.0;
               for (int d = 0; d < dim; d++)
               {
                  grad_dot_n += dshape1_adj(k, d) * nor_q(d);
               }

               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape1_adj(k, i) * nor_q(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape1_adj(k, u) * nor_q(i));
                     sym_val += trac_iu * sign_q * delta_u[u];
                  }

                  int idx = i * ndof1 + k;  // byNODES ordering
                  elvec1(idx) += c1 * wq * sym_val / detJ1;
                  elvec1(idx) += penalty * wq * shapes1(k, q) * f_lifted_q(i, q);
               }
            }

            // Element 2 consistency + BR2
            for (int k = 0; k < ndof2; k++)
            {
               real_t grad_dot_n = 0.0;
               for (int d = 0; d < dim; d++)
               {
                  grad_dot_n += dshape2_adj(k, d) * nor_q(d);
               }

               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape2_adj(k, i) * nor_q(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape2_adj(k, u) * nor_q(i));
                     sym_val += trac_iu * sign_q * delta_u[u];
                  }

                  int idx = i * ndof2 + k;
                  elvec2(idx) += c1 * wq * sym_val / detJ2;
                  elvec2(idx) += (-penalty) * wq * shapes2(k, q) * f_lifted_q(i, q);
               }
            }
         }

         // Add to global RHS
         for (int j = 0; j < vdofs1.Size(); j++)
         {
            int gj = vdofs1[j];
            if (gj >= 0) { rhs(gj) += elvec1(j); }
            else { rhs(-1 - gj) -= elvec1(j); }
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            int gj = vdofs2[j];
            if (gj >= 0) { rhs(gj) += elvec2(j); }
            else { rhs(-1 - gj) -= elvec2(j); }
         }
      }
   }

   // ========================================================================
   // Shared Face Slip Assembly — IP method (Parallel)
   // ========================================================================

   void AssembleSlipContributionIPShared(Vector &rhs, const Vector &slip_bc,
                                         int interior_face_count) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int dim = 3;
         real_t kappa = (order_ + 1) * (order_ + 1);

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            // Slip index: interior faces first, then shared faces
            int slip_idx = interior_face_count + i;
            real_t slip_local[2] = {slip_bc(2 * slip_idx),
                                    slip_bc(2 * slip_idx + 1)};
            real_t delta_u[3];
            fault_basis_.EmbedSlip(slip_idx, slip_local, delta_u);

            real_t max_delta = std::max({std::abs(delta_u[0]),
                                          std::abs(delta_u[1]),
                                          std::abs(delta_u[2])});
            if (max_delta < 1e-15) { continue; }

            // Only Elem1 is local
            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            int ndof1 = fe1->GetDof();

            // Get Elem2 (face-neighbor) FE for detJ2
            auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

            int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom,
                                                     2 * face_order + 1);

            Vector elvec1(vdofs1.Size());
            elvec1 = 0.0;

            for (int p = 0; p < ir.GetNPoints(); p++)
            {
               const IntegrationPoint &ip = ir.IntPoint(p);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

               Vector nor(dim);
               CalcOrtho(FTr->Jacobian(), nor);

               real_t sign = (nor(0) > 0) ? -1.0 : 1.0;

               Vector shape1(ndof1);
               fe1->CalcShape(eip1, shape1);

               DenseMatrix dshape1_ref(ndof1, dim);
               fe1->CalcDShape(eip1, dshape1_ref);

               DenseMatrix adjJ1(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);

               DenseMatrix dshape1_adj(ndof1, dim);
               Mult(dshape1_ref, adjJ1, dshape1_adj);

               real_t detJ1 = FTr->Elem1->Weight();
               real_t detJ2 = FTr->Elem2->Weight();
               real_t w1 = ip.weight / (2.0 * detJ1);
               real_t w2 = ip.weight / (2.0 * detJ2);
               real_t nor_sq = nor * nor;
               real_t wq_penalty = kappa * nor_sq * (w1 + w2);

               for (int k = 0; k < ndof1; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape1_adj(k, d) * nor(d);
                  }

                  for (int ci = 0; ci < dim; ci++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac_iu = lambda_val_ * dshape1_adj(k, ci) * nor(u)
                           + mu_val_ * ((ci == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape1_adj(k, u) * nor(ci));
                        sym_val += trac_iu * sign * delta_u[u];
                     }

                     int idx = ci * ndof1 + k;
                     elvec1(idx) += epsilon_ * sym_val * w1;
                     elvec1(idx) += wq_penalty * sign * delta_u[ci] * shape1(k);
                  }
               }
            }

            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
#endif
      }
   }

   // ========================================================================
   // Shared Face Slip Assembly — BR2 method (Parallel)
   // ========================================================================

   void AssembleSlipContributionBR2Shared(Vector &rhs, const Vector &slip_bc,
                                           int interior_face_count) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int dim = 3;
         if (!mass_inv_computed_) { PrecomputeMassInverse(); }

         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

         for (int fi = 0; fi < fault_shared_faces_.Size(); fi++)
         {
            int sf = fault_shared_faces_[fi];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int slip_idx = interior_face_count + fi;
            real_t slip_local[2] = {slip_bc(2 * slip_idx),
                                    slip_bc(2 * slip_idx + 1)};
            real_t delta_u[3];
            fault_basis_.EmbedSlip(slip_idx, slip_local, delta_u);

            real_t max_delta = std::max({std::abs(delta_u[0]),
                                          std::abs(delta_u[1]),
                                          std::abs(delta_u[2])});
            if (max_delta < 1e-15) { continue; }

            // BR2 penalty
            Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t penalty = (geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            int ndof1 = fe1->GetDof();

            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            int ndof2 = fe2->GetDof();

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

            int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom,
                                                     2 * face_order + 1);
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_all(dim * nqp);
            Vector w_all(nqp);
            Vector slip_sign_all(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

               Vector s1q(shapes1.GetColumn(q), ndof1);
               Vector s2q(shapes2.GetColumn(q), ndof2);
               fe1->CalcShape(eip1, s1q);
               fe2->CalcShape(eip2, s2q);

               Vector nor_q(&nor_all[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nor_q);

               w_all[q] = ip.weight;
               slip_sign_all[q] = (nor_q(0) > 0) ? -1.0 : 1.0;
            }

            // Two-sided BR2 lifting (both elements contribute to lifting)
            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_all[q];
               real_t sign_q = slip_sign_all[q];

               for (int u = 0; u < dim; u++)
               {
                  real_t delta_u_signed = sign_q * delta_u[u];
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_all[q * dim + s];
                     real_t factor = delta_u_signed * n_s * wq;

                     for (int m = 0; m < ndof1; m++)
                     {
                        face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                     }
                     for (int m = 0; m < ndof2; m++)
                     {
                        face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                     }
                  }
               }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            f_lifted1 *= 0.5;
            MultABt(face_int2, Minv2, f_lifted2);
            f_lifted2 *= 0.5;

            // f_lifted_q[i,q] with elasticity tensor coupling
            DenseMatrix f_lifted_q(dim, nqp);
            f_lifted_q = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               Vector n_q(dim);
               for (int d = 0; d < dim; d++) { n_q(d) = nor_all[q * dim + d]; }

               for (int ci = 0; ci < dim; ci++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(ci)
                           + mu_val_ * ((ci == u ? 1.0 : 0.0) * n_q(s)
                                        + (ci == s ? 1.0 : 0.0) * n_q(u));

                        real_t eval1 = 0.0, eval2 = 0.0;
                        for (int m = 0; m < ndof1; m++)
                        {
                           eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                        }
                        for (int m = 0; m < ndof2; m++)
                        {
                           eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                        }
                        sum += tn * (eval1 + eval2);
                     }
                  }
                  f_lifted_q(ci, q) = 0.5 * sum;
               }
            }

            // Assemble RHS (Elem1 only — local element)
            Vector elvec1(vdofs1.Size());
            elvec1 = 0.0;

            real_t c1 = epsilon_ * 0.5;

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

               real_t wq = w_all[q];
               real_t sign_q = slip_sign_all[q];

               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_all[q * dim + d]; }

               DenseMatrix dshape1_ref(ndof1, dim);
               fe1->CalcDShape(eip1, dshape1_ref);

               DenseMatrix adjJ1(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);

               DenseMatrix dshape1_adj(ndof1, dim);
               Mult(dshape1_ref, adjJ1, dshape1_adj);

               real_t detJ1 = FTr->Elem1->Weight();

               for (int k = 0; k < ndof1; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape1_adj(k, d) * nor_q(d);
                  }

                  for (int ci = 0; ci < dim; ci++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac_iu = lambda_val_ * dshape1_adj(k, ci) * nor_q(u)
                           + mu_val_ * ((ci == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape1_adj(k, u) * nor_q(ci));
                        sym_val += trac_iu * sign_q * delta_u[u];
                     }

                     int idx = ci * ndof1 + k;
                     elvec1(idx) += c1 * wq * sym_val / detJ1;
                     elvec1(idx) += penalty * wq * shapes1(k, q) * f_lifted_q(ci, q);
                  }
               }
            }

            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
#endif
      }
   }

   // ========================================================================
   // Dirichlet loading assembly
   // ========================================================================

   void AssembleDirichletLoading(Vector &rhs, real_t time) const
   {
      // Plate loading on ±x1 walls (fault-normal): u₂ = ±Vp·t/2
      // We need to assemble the DG Dirichlet BC contribution.
      // For each boundary face on the ±x1 walls, add:
      //   b[k,i] += c0 * [σ(φ_k e_i)·n]_u * u_D_u * (1/detJ)
      //           + penalty * φ_k * u_D_i
      //
      // Attr 1 = -x wall (negative fault-normal side) → u_y = -Vp*t/2
      // Attr 2 = +x wall (positive fault-normal side) → u_y = +Vp*t/2
      // This creates ∂u_y/∂x → ε_xy → σ_xy which drives shear on the fault.

      if (std::abs(time * Vp_) < 1e-30) { return; }

      int dim = 3;

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr != 1 && attr != 2) { continue; }  // Only ±x1 walls (fault-normal)

         // Dirichlet value: u = (0, ±Vp*t/2, 0)
         real_t u_D[3] = {0.0, 0.0, 0.0};
         if (attr == 2) { u_D[1] = Vp_ * time / 2.0; }   // +x1
         else           { u_D[1] = -Vp_ * time / 2.0; }   // -x1

         // Get face transformation
         int face_idx, face_info;
         mesh_.GetBdrElementFace(be, &face_idx, &face_info);
         FaceElementTransformations *FTr =
            mesh_.GetFaceElementTransformations(face_idx);
         if (FTr == nullptr) { continue; }

         Array<int> vdofs;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs);
         const FiniteElement *fe = scalar_fes_->GetFE(FTr->Elem1No);
         int ndof = fe->GetDof();

         int face_order = fe->GetOrder();
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2 * face_order + 1);

         Vector elvec(vdofs.Size());
         elvec = 0.0;

         if (method_ == DGMethod::IP)
         {
            real_t kappa = (order_ + 1) * (order_ + 1);
            for (int p = 0; p < ir.GetNPoints(); p++)
            {
               const IntegrationPoint &ip = ir.IntPoint(p);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip = FTr->GetElement1IntPoint();

               Vector nor(dim);
               CalcOrtho(FTr->Jacobian(), nor);

               Vector shape(ndof);
               fe->CalcShape(eip, shape);

               DenseMatrix dshape_ref(ndof, dim);
               fe->CalcDShape(eip, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof, dim);
               Mult(dshape_ref, adjJ, dshape_adj);

               real_t detJ = FTr->Elem1->Weight();
               real_t w = ip.weight / detJ;
               real_t nor_sq = nor * nor;
               real_t wq_penalty = kappa * nor_sq * w;

               for (int k = 0; k < ndof; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape_adj(k, d) * nor(d);
                  }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac = lambda_val_ * dshape_adj(k, i) * nor(u)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape_adj(k, u) * nor(i));
                        sym_val += trac * u_D[u];
                     }

                     int idx = i * ndof + k;
                     elvec(idx) += epsilon_ * sym_val * w;
                     elvec(idx) += wq_penalty * u_D[i] * shape(k);
                  }
               }
            }
         }
         else  // BR2
         {
            const DenseMatrix &Minv = elem_mass_inv_[FTr->Elem1No];
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes(ndof, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip = FTr->GetElement1IntPoint();

               Vector sq(shapes.GetColumn(q), ndof);
               fe->CalcShape(eip, sq);

               Vector nq(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nq);
               w_arr[q] = ip.weight;
            }

            // Compute lifted Dirichlet: similar to slip but for boundary
            DenseMatrix face_int(dim * dim, ndof);
            face_int = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_arr[q];
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_arr[q * dim + s];
                     real_t factor = u_D[u] * n_s * wq;
                     for (int m = 0; m < ndof; m++)
                     {
                        face_int(u * dim + s, m) += shapes(m, q) * factor;
                     }
                  }
               }
            }

            DenseMatrix f_lifted(dim * dim, ndof);
            MultABt(face_int, Minv, f_lifted);
            // No 0.5 for boundary (full factor)

            // f_lifted_q[i,q]
            DenseMatrix fl_q(dim, nqp);
            fl_q = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               Vector n_q(dim);
               for (int d = 0; d < dim; d++) { n_q(d) = nor_arr[q * dim + d]; }

               for (int i = 0; i < dim; i++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                        + (i == s ? 1.0 : 0.0) * n_q(u));
                        real_t ev = 0.0;
                        for (int m = 0; m < ndof; m++)
                        {
                           ev += shapes(m, q) * f_lifted(u * dim + s, m);
                        }
                        sum += tn * ev;
                     }
                  }
                  fl_q(i, q) = sum;
               }
            }

            // BR2 penalty: num faces per element (geometry-dependent)
            Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip = FTr->GetElement1IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof, dim);
               fe->CalcDShape(eip, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem1->Weight();

               for (int k = 0; k < ndof; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape_adj(k, d) * nor_q(d);
                  }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac = lambda_val_ * dshape_adj(k, i) * nor_q(u)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape_adj(k, u) * nor_q(i));
                        sym_val += trac * u_D[u];
                     }

                     int idx = i * ndof + k;
                     elvec(idx) += epsilon_ * wq * sym_val / detJ;
                     elvec(idx) += br2_pen * wq * shapes(k, q) * fl_q(i, q);
                  }
               }
            }
         }

         for (int j = 0; j < vdofs.Size(); j++)
         {
            int gj = vdofs[j];
            if (gj >= 0) { rhs(gj) += elvec(j); }
            else { rhs(-1 - gj) -= elvec(j); }
         }
      }
   }
};

// ============================================================================
// Out-of-class implementation of remaining methods
// ============================================================================

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDepths(Vector &depths) const
{
   if (!fault_depths_computed_)
   {
      fault_depths_.SetSize(num_fault_dofs_);

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->Face->SetIntPoint(&ip);
         Vector coords(3);
         FTr->Face->Transform(ip, coords);

         // x3 = depth (z-coordinate)
         fault_depths_(i) = coords(2);
      }

      // Shared fault faces (parallel only)
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            const IntegrationPoint &ip =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);

            int idx = fault_interior_faces_.Size() + i;
            fault_depths_(idx) = coords(2);
         }
#endif
      }

      fault_depths_computed_ = true;
   }

   depths = fault_depths_;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultCoords2D(
   Vector &coords_x2, Vector &coords_x3) const
{
   if (!fault_coords_computed_)
   {
      fault_x2_.SetSize(num_fault_dofs_);
      fault_x3_.SetSize(num_fault_dofs_);

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->Face->SetIntPoint(&ip);
         Vector coords(3);
         FTr->Face->Transform(ip, coords);

         // x2 = along-strike (y-coordinate), x3 = depth (z-coordinate)
         fault_x2_(i) = coords(1);
         fault_x3_(i) = coords(2);
      }

      // Shared fault faces (parallel only)
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            const IntegrationPoint &ip =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);

            int idx = fault_interior_faces_.Size() + i;
            fault_x2_(idx) = coords(1);
            fault_x3_(idx) = coords(2);
         }
#endif
      }

      fault_coords_computed_ = true;
   }

   coords_x2 = fault_x2_;
   coords_x3 = fault_x3_;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::Solve(
   real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
   if (!stiffness_assembled_)
   {
      AssembleStiffness();
   }

   // Build RHS
   LinFormType b(fes_.get());
   b.Assemble();
   Vector &rhs = b;

   // Add slip contributions (interior + shared faces)
   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);
      AssembleSlipContributionIPShared(rhs, slip_bc,
                                       fault_interior_faces_.Size());
   }
   else
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
      AssembleSlipContributionBR2Shared(rhs, slip_bc,
                                        fault_interior_faces_.Size());
   }

   // Add Dirichlet loading
   AssembleDirichletLoading(rhs, time);

   X_ = 0.0;
   B_ = rhs;

   solver_->Mult(B_, X_);

   // Check convergence
   {
      auto *cg = dynamic_cast<CGSolver*>(solver_.get());
      if (cg && !cg->GetConverged())
      {
         mfem::err << "WARNING: CG did not converge after "
                   << cg->GetNumIterations() << " iterations, final norm = "
                   << cg->GetFinalNorm() << "\n";
      }
   }

   displacement = X_;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction)
{
   int dim = 3;
   traction.SetSize(2 * num_fault_dofs_);
   traction = 0.0;

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      int face = fault_interior_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      // Evaluate face centroid
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());
      FTr->SetAllIntPoints(&ip);
      const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
      const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

      // Get DOFs
      Array<int> vdofs1, vdofs2;
      fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
      fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

      const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
      const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
      int ndof1 = fe1->GetDof();
      int ndof2 = fe2->GetDof();

      // Get displacement values
      Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
      displacement.GetSubVector(vdofs1, u1_all);
      displacement.GetSubVector(vdofs2, u2_all);

      // Compute gradient from each side: grad_u[d1, d2] = du_{d1}/dx_{d2}
      DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
      fe1->CalcDShape(eip1, dshape1_ref);
      fe2->CalcDShape(eip2, dshape2_ref);

      DenseMatrix Jinv1(dim), Jinv2(dim);
      CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
      CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

      DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
      Mult(dshape1_ref, Jinv1, dshape1_phys);
      Mult(dshape2_ref, Jinv2, dshape2_phys);

      // Compute ∇u on each side (byNODES: u1_all has ndof1*3 entries)
      // u_component_c at DOF k: u1_all(c * ndof1 + k)
      DenseMatrix grad1(dim, dim), grad2(dim, dim);
      grad1 = 0.0;
      grad2 = 0.0;

      for (int c = 0; c < dim; c++)
      {
         for (int d = 0; d < dim; d++)
         {
            real_t val1 = 0.0, val2 = 0.0;
            for (int k = 0; k < ndof1; k++)
            {
               val1 += dshape1_phys(k, d) * u1_all(c * ndof1 + k);
            }
            for (int k = 0; k < ndof2; k++)
            {
               val2 += dshape2_phys(k, d) * u2_all(c * ndof2 + k);
            }
            grad1(c, d) = val1;
            grad2(c, d) = val2;
         }
      }

      // Average gradient
      DenseMatrix avg_grad(dim, dim);
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++)
            avg_grad(i, j) = 0.5 * (grad1(i, j) + grad2(i, j));

      // Strain tensor: ε_{ij} = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
      DenseMatrix strain(dim, dim);
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++)
            strain(i, j) = 0.5 * (avg_grad(i, j) + avg_grad(j, i));

      // Stress tensor: σ_{ij} = λ * tr(ε) * δ_{ij} + 2μ * ε_{ij}
      real_t tr_eps = strain(0, 0) + strain(1, 1) + strain(2, 2);

      DenseMatrix stress(dim, dim);
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++)
            stress(i, j) = lambda_val_ * tr_eps * (i == j ? 1.0 : 0.0)
                           + 2.0 * mu_val_ * strain(i, j);

      // Traction: T = σ · n (using fault basis normal)
      const auto &basis = fault_basis_.GetBasis(fi);
      real_t T_global[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < dim; i++)
      {
         for (int j = 0; j < dim; j++)
         {
            T_global[i] += stress(i, j) * basis.normal[j];
         }
      }

      // === DG penalty correction: T -= penalty * ([[u]] - δ) ===
      // This matches Tandem's traction formula:
      //   t = {σ}·n + c0 * (E_q[0]*u[0] - E_q[1]*u[1] - f_q)
      // where c0 = -penalty, so t = {σ}·n - penalty * ([[u]] - δ)

      // 1. Compute displacement values at face centroid
      Vector shape1(ndof1), shape2(ndof2);
      fe1->CalcShape(eip1, shape1);
      fe2->CalcShape(eip2, shape2);

      real_t u1_val[3] = {0.0, 0.0, 0.0};
      real_t u2_val[3] = {0.0, 0.0, 0.0};
      for (int c = 0; c < dim; c++)
      {
         for (int k = 0; k < ndof1; k++)
         {
            u1_val[c] += shape1(k) * u1_all(c * ndof1 + k);
         }
         for (int k = 0; k < ndof2; k++)
         {
            u2_val[c] += shape2(k) * u2_all(c * ndof2 + k);
         }
      }

      // 2. Displacement jump [[u]] = u1 - u2
      real_t u_jump[3];
      for (int c = 0; c < dim; c++)
      {
         u_jump[c] = u1_val[c] - u2_val[c];
      }

      // 3. Prescribed slip in global frame
      real_t slip_local[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
      real_t delta_u[3];
      fault_basis_.EmbedSlip(fi, slip_local, delta_u);

      // 4. Sign correction (same convention as slip assembly)
      Vector nor(dim);
      CalcOrtho(FTr->Jacobian(), nor);
      real_t sign = (nor(0) > 0) ? -1.0 : 1.0;

      // 5. Compute penalty correction based on DG method
      real_t correction[3] = {0.0, 0.0, 0.0};

      if (method_ == DGMethod::IP)
      {
         // IP penalty: kappa * |nor|^2 * (1/(2*detJ1) + 1/(2*detJ2))
         real_t kappa = (order_ + 1) * (order_ + 1);
         real_t detJ1 = FTr->Elem1->Weight();
         real_t detJ2 = FTr->Elem2->Weight();
         real_t nor_sq = nor * nor;
         real_t penalty = kappa * nor_sq * (1.0 / (2.0 * detJ1) + 1.0 / (2.0 * detJ2));

         for (int c = 0; c < dim; c++)
         {
            correction[c] = penalty * (u_jump[c] - sign * delta_u[c]);
         }
      }
      else  // BR2
      {
         if (!mass_inv_computed_) { PrecomputeMassInverse(); }

         Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
         real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                  ? real_t(dim + 1) : real_t(2 * dim);

         const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
         const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

         // Compute the jump to penalize: [[u]] - sign * δ
         real_t jump[3];
         for (int c = 0; c < dim; c++)
         {
            jump[c] = u_jump[c] - sign * delta_u[c];
         }

         // BR2 lifting at face centroid
         // face_int[u*dim+s, m] = shape[m] * jump[u] * nor[s]
         DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
         face_int1 = 0.0;
         face_int2 = 0.0;
         for (int u = 0; u < dim; u++)
         {
            for (int s = 0; s < dim; s++)
            {
               for (int m = 0; m < ndof1; m++)
               {
                  face_int1(u * dim + s, m) = shape1(m) * jump[u] * nor(s);
               }
               for (int m = 0; m < ndof2; m++)
               {
                  face_int2(u * dim + s, m) = shape2(m) * jump[u] * nor(s);
               }
            }
         }

         // f_lifted = 0.5 * face_int * Minv^T
         DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
         MultABt(face_int1, Minv1, f_lifted1);
         f_lifted1 *= 0.5;
         MultABt(face_int2, Minv2, f_lifted2);
         f_lifted2 *= 0.5;

         // Evaluate f_lifted_q at centroid with elasticity tensor coupling
         for (int i = 0; i < dim; i++)
         {
            real_t sum = 0.0;
            for (int u = 0; u < dim; u++)
            {
               for (int s = 0; s < dim; s++)
               {
                  real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
                     + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                                + (i == s ? 1.0 : 0.0) * basis.normal[u]);
                  real_t eval1 = 0.0, eval2 = 0.0;
                  for (int m = 0; m < ndof1; m++)
                  {
                     eval1 += shape1(m) * f_lifted1(u * dim + s, m);
                  }
                  for (int m = 0; m < ndof2; m++)
                  {
                     eval2 += shape2(m) * f_lifted2(u * dim + s, m);
                  }
                  sum += tn * (eval1 + eval2);
               }
            }
            correction[i] = br2_penalty * 0.5 * sum;
         }
      }

      // 6. Apply correction: T -= penalty * ([[u]] - δ)
      for (int c = 0; c < dim; c++)
      {
         T_global[c] -= correction[c];
      }

      // Project to local frame: (tau_dip, tau_strike)
      real_t tau_local[2];
      fault_basis_.ProjectTraction(fi, T_global, tau_local);
      traction(2 * fi)     = tau_local[0];
      traction(2 * fi + 1) = tau_local[1];
   }

   // Shared fault faces (parallel only)
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      ParGridFunction par_u(pfes);
      par_u = displacement;
      pfes->ExchangeFaceNbrData();
      par_u.ExchangeFaceNbrData();

      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }

         int trac_idx = fault_interior_faces_.Size() + i;

         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->SetAllIntPoints(&ip);
         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         // Elem1 (local)
         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         int ndof1 = fe1->GetDof();

         Array<int> vdofs1;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         Vector u1_all(vdofs1.Size());
         par_u.GetSubVector(vdofs1, u1_all);

         // Elem2 (face-neighbor)
         int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
         int ndof2 = fe2->GetDof();

         // Get face-neighbor DOF values (vector space)
         Array<int> vdofs2;
         pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
         const Vector &nbr_data = par_u.FaceNbrData();
         Vector u2_all(vdofs2.Size());
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            u2_all(j) = nbr_data(vdofs2[j]);
         }

         // Compute gradients
         DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
         fe1->CalcDShape(eip1, dshape1_ref);
         fe2->CalcDShape(eip2, dshape2_ref);

         DenseMatrix Jinv1(dim), Jinv2(dim);
         CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
         CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

         DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
         Mult(dshape1_ref, Jinv1, dshape1_phys);
         Mult(dshape2_ref, Jinv2, dshape2_phys);

         // Compute ∇u on each side (byNODES ordering)
         DenseMatrix grad1(dim, dim), grad2(dim, dim);
         grad1 = 0.0;
         grad2 = 0.0;

         for (int c = 0; c < dim; c++)
         {
            for (int d = 0; d < dim; d++)
            {
               real_t val1 = 0.0, val2 = 0.0;
               for (int k = 0; k < ndof1; k++)
               {
                  val1 += dshape1_phys(k, d) * u1_all(c * ndof1 + k);
               }
               for (int k = 0; k < ndof2; k++)
               {
                  val2 += dshape2_phys(k, d) * u2_all(c * ndof2 + k);
               }
               grad1(c, d) = val1;
               grad2(c, d) = val2;
            }
         }

         // Average gradient → strain → stress
         DenseMatrix avg_grad(dim, dim);
         for (int ci = 0; ci < dim; ci++)
            for (int cj = 0; cj < dim; cj++)
               avg_grad(ci, cj) = 0.5 * (grad1(ci, cj) + grad2(ci, cj));

         DenseMatrix strain(dim, dim);
         for (int ci = 0; ci < dim; ci++)
            for (int cj = 0; cj < dim; cj++)
               strain(ci, cj) = 0.5 * (avg_grad(ci, cj) + avg_grad(cj, ci));

         real_t tr_eps = strain(0, 0) + strain(1, 1) + strain(2, 2);

         DenseMatrix stress(dim, dim);
         for (int ci = 0; ci < dim; ci++)
            for (int cj = 0; cj < dim; cj++)
               stress(ci, cj) = lambda_val_ * tr_eps * (ci == cj ? 1.0 : 0.0)
                                 + 2.0 * mu_val_ * strain(ci, cj);

         // Traction: T = σ · n
         const auto &basis = fault_basis_.GetBasis(trac_idx);
         real_t T_global[3] = {0.0, 0.0, 0.0};
         for (int ci = 0; ci < dim; ci++)
            for (int cj = 0; cj < dim; cj++)
               T_global[ci] += stress(ci, cj) * basis.normal[cj];

         // DG penalty correction: T -= penalty * ([[u]] - δ)
         Vector shape1(ndof1), shape2(ndof2);
         fe1->CalcShape(eip1, shape1);
         fe2->CalcShape(eip2, shape2);

         real_t u1_val[3] = {0.0, 0.0, 0.0};
         real_t u2_val[3] = {0.0, 0.0, 0.0};
         for (int c = 0; c < dim; c++)
         {
            for (int k = 0; k < ndof1; k++)
               u1_val[c] += shape1(k) * u1_all(c * ndof1 + k);
            for (int k = 0; k < ndof2; k++)
               u2_val[c] += shape2(k) * u2_all(c * ndof2 + k);
         }

         real_t u_jump[3];
         for (int c = 0; c < dim; c++)
            u_jump[c] = u1_val[c] - u2_val[c];

         real_t slip_local[2] = {slip_bc(2 * trac_idx),
                                 slip_bc(2 * trac_idx + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(trac_idx, slip_local, delta_u);

         Vector nor(dim);
         CalcOrtho(FTr->Jacobian(), nor);
         real_t sign = (nor(0) > 0) ? -1.0 : 1.0;

         real_t correction[3] = {0.0, 0.0, 0.0};

         if (method_ == DGMethod::IP)
         {
            real_t kappa = (order_ + 1) * (order_ + 1);
            real_t detJ1 = FTr->Elem1->Weight();
            real_t detJ2 = FTr->Elem2->Weight();
            real_t nor_sq = nor * nor;
            real_t pen = kappa * nor_sq * (1.0 / (2.0 * detJ1) + 1.0 / (2.0 * detJ2));

            for (int c = 0; c < dim; c++)
               correction[c] = pen * (u_jump[c] - sign * delta_u[c]);
         }
         else  // BR2
         {
            if (!mass_inv_computed_) { PrecomputeMassInverse(); }

            Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                     ? real_t(dim + 1) : real_t(2 * dim);

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

            real_t jump[3];
            for (int c = 0; c < dim; c++)
               jump[c] = u_jump[c] - sign * delta_u[c];

            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;
            for (int u = 0; u < dim; u++)
            {
               for (int s = 0; s < dim; s++)
               {
                  for (int m = 0; m < ndof1; m++)
                     face_int1(u * dim + s, m) = shape1(m) * jump[u] * nor(s);
                  for (int m = 0; m < ndof2; m++)
                     face_int2(u * dim + s, m) = shape2(m) * jump[u] * nor(s);
               }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            f_lifted1 *= 0.5;
            MultABt(face_int2, Minv2, f_lifted2);
            f_lifted2 *= 0.5;

            for (int ci = 0; ci < dim; ci++)
            {
               real_t sum = 0.0;
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[ci]
                        + mu_val_ * ((ci == u ? 1.0 : 0.0) * basis.normal[s]
                                     + (ci == s ? 1.0 : 0.0) * basis.normal[u]);
                     real_t eval1 = 0.0, eval2 = 0.0;
                     for (int m = 0; m < ndof1; m++)
                        eval1 += shape1(m) * f_lifted1(u * dim + s, m);
                     for (int m = 0; m < ndof2; m++)
                        eval2 += shape2(m) * f_lifted2(u * dim + s, m);
                     sum += tn * (eval1 + eval2);
                  }
               }
               correction[ci] = br2_penalty * 0.5 * sum;
            }
         }

         for (int c = 0; c < dim; c++)
            T_global[c] -= correction[c];

         real_t tau_local[2];
         fault_basis_.ProjectTraction(trac_idx, T_global, tau_local);
         traction(2 * trac_idx)     = tau_local[0];
         traction(2 * trac_idx + 1) = tau_local[1];
      }
#endif
   }
}

// Convenience type alias
using SerialElasticityDomainOperator = ElasticityDomainOperator<Mesh>;

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ELASTICITY_OPERATOR_HPP
