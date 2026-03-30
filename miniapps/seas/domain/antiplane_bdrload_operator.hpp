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

#ifndef MFEM_SEAS_ANTIPLANE_BDRLOAD_OPERATOR_HPP
#define MFEM_SEAS_ANTIPLANE_BDRLOAD_OPERATOR_HPP

/// @file antiplane_bdrload_operator.hpp
/// @brief DG antiplane operator with far-field Dirichlet loading for BP1.
///
/// This operator is an alternative to AntiplaneDomainOperator that drives
/// tectonic loading via Dirichlet BCs on the far-field left/right boundaries:
///   u(x=-Lx, z, t) = -Vp/2 * t   (left far-field)
///   u(x=+Lx, z, t) = +Vp/2 * t   (right far-field)
///
/// This matches Tandem's BP1 boundary condition setup where:
///   Physical Curve(1) = Natural  (top + bottom)
///   Physical Curve(3) = Fault    (x=0, z=0 to z=-Wf)
///   Physical Curve(5) = Dirichlet (far-field + deep extension below Wf)
///
/// The stiffness matrix includes DG boundary-face penalty on far-field boundaries,
/// and the RHS includes DGDirichletLFIntegrator contributions from the
/// time-dependent far-field BC at each Solve() call.

#include "mfem.hpp"
#include "domain_operator.hpp"
#include "bp2_mesh.hpp"
#include "seas_boundary_tags.hpp"  // For FaultBoundaryData
#include "../common/seas_types.hpp"
#include "../integrator/dg_br2_integrator.hpp"
#include "antiplane_operator.hpp"  // For DGMethod enum

#include <memory>

namespace mfem
{
namespace seas
{

/// @brief Time-dependent coefficient for far-field Dirichlet loading.
///
/// Returns u(x, z, t) = sign(x) * Vp/2 * t
/// Applied on left (x=-Lx) and right (x=+Lx) far-field boundaries.
class FarFieldLoadingCoefficient : public Coefficient
{
public:
   FarFieldLoadingCoefficient(real_t Vp) : Vp_(Vp), time_(0.0) {}

   void SetTime(real_t t) override { time_ = t; }

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      real_t sign_x = (x(0) >= 0.0) ? 1.0 : -1.0;
      return sign_x * 0.5 * Vp_ * time_;
   }

private:
   real_t Vp_;
   real_t time_;
};

/// @brief DG Antiplane domain operator with far-field Dirichlet loading.
///
/// Identical to AntiplaneDomainOperator except:
/// 1. AssembleStiffness() adds DG boundary-face penalty on far-field boundaries
/// 2. Solve() adds DGDirichletLFIntegrator RHS from the time-dependent far-field BC
///
/// Boundary conditions (matching Tandem's BP1 setup):
///   - Far-field left (x=-Lx): Dirichlet: u = -Vp/2 * t
///   - Far-field right (x=+Lx): Dirichlet: u = +Vp/2 * t
///   - Free surface (z=0): Natural BC (zero traction)
///   - Bottom (z=-Lz): Natural BC (zero traction)
///   - Fault (x=0, interior): Slip jump [[u]] = delta
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class AntiplaneBdrLoadOperator : public DomainOperator<MeshType>
{
public:
   using Base = DomainOperator<MeshType>;
   using FESpaceType = typename Base::FESpaceType;
   using GridFuncType = typename Base::GridFuncType;
   using BilinFormType = typename Base::BilinFormType;
   using LinFormType = typename Base::LinFormType;

   /// @brief Construct the DG antiplane operator with bottom Dirichlet loading.
   ///
   /// @param mesh The computational mesh (full domain, owned externally)
   /// @param order Polynomial order for DG finite elements
   /// @param mu Shear modulus [Pa]
   /// @param Vp Plate rate [m/s]
   /// @param Wf Fault depth [m] - used for fault geometry, not for loading cutoff
   /// @param method DG method (IP or BR2), default is IP
   /// @param fault_tag Gmsh Physical Curve tag for fault detection.
   ///                  If >= 1, uses tag-based detection from boundary elements.
   ///                  If -1 (default), uses legacy coordinate-based detection.
   AntiplaneBdrLoadOperator(MeshType &mesh, int order, real_t mu, real_t Vp,
                            real_t Wf = 40.0e3, DGMethod method = DGMethod::IP,
                            int fault_tag = -1);

   ~AntiplaneBdrLoadOperator() override = default;

   int NumComponents() const override { return 1; }
   int Dimension() const override { return mesh_.Dimension(); }

   void Solve(real_t time, const Vector &slip_bc,
              GridFuncType &displacement) override;

   void ComputeTraction(const GridFuncType &displacement,
                        const Vector &slip_bc,
                        Vector &traction,
                        Vector *normal_traction = nullptr) override;

   FESpaceType &GetFESpace() override { return *fes_; }
   const FESpaceType &GetFESpace() const override { return *fes_; }

   MeshType &GetMesh() override { return mesh_; }
   const MeshType &GetMesh() const override { return mesh_; }

   real_t GetShearModulus() const override { return mu_; }
   int GetNumFaultDOFs() const override { return num_fault_dofs_; }
   void GetFaultDepths(Vector &depths) const override;
   const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

   const Array<int> &GetFaultInteriorFaces() const { return fault_interior_faces_; }
   const Array<int> &GetFaultSharedFaces() const { return fault_shared_faces_; }

   real_t GetPlateRate() const { return Vp_; }
   real_t GetFaultDepth() const { return Wf_; }
   int GetOrder() const { return order_; }
   DGMethod GetMethod() const { return method_; }

   /// @brief Get the fault tag used for detection (-1 = coordinate-based)
   int GetFaultTag() const { return fault_tag_; }

   real_t InterpolateSlipBC(real_t z, const Vector &slip_bc) const;

private:
   MeshType &mesh_;
   int order_;
   real_t mu_;
   real_t Vp_;
   real_t Wf_;

   // Fault detection: tag-based (>= 1) or coordinate-based (-1)
   int fault_tag_;

   DGMethod method_;
   real_t sigma_;       // SIPG sign (-1)
   real_t br2_penalty_; // BR2 penalty (D+1 = 3)

   // Precomputed element mass matrix inverses for BR2
   mutable std::vector<DenseMatrix> elem_mass_inv_;
   mutable bool mass_inv_computed_;

   std::unique_ptr<FiniteElementCollection> fec_;
   std::unique_ptr<FESpaceType> fes_;

   // Far-field boundary marker for Dirichlet BC (LEFT + RIGHT)
   mutable Array<int> farfield_bdr_;

   // Interior fault faces (at x=0)
   Array<int> fault_interior_faces_;
   Array<int> fault_shared_faces_;

   // Fault DOF information
   Array<int> fault_dofs_;
   int num_fault_dofs_;

   mutable Vector fault_depths_;
   mutable bool fault_depths_computed_;

   // Solver components
   mutable std::unique_ptr<Solver> solver_;
   std::unique_ptr<Solver> prec_;
   mutable std::unique_ptr<GSSmoother> serial_prec_;
   mutable Vector X_, B_;

   // Cached stiffness matrix
   mutable bool stiffness_assembled_ = false;
   mutable std::unique_ptr<BilinFormType> cached_a_;
   mutable OperatorHandle cached_Ah_;
   mutable std::unique_ptr<Solver> cached_prec_;

   void AssembleStiffness() const;

   void SetupBoundaryMarkers();
   void SetupFESpace();
   void SetupFaultInfo();
   void SetupSolver();

   bool IsFaultFace(int face) const;
   bool IsFaultFaceShared(int shared_face) const;
   void GetFaceCenter(int face, real_t &x, real_t &z) const;
   void GetFaceCenterShared(int shared_face, real_t &x, real_t &z) const;

   void PrecomputeMassMatrixInverses() const;
   const DenseMatrix& GetMassMatrixInverse(int elem) const;
   void ComputeLiftingOperator(int elem, FaceElementTransformations &FTr,
                               real_t jump, Vector &lift) const;

   void AssembleSlipContributionIP(Vector &rhs, const Vector &slip_bc) const;
   void AssembleSlipContributionBR2(Vector &rhs, const Vector &slip_bc) const;
   void AssembleSlipContributionIPShared(Vector &rhs, const Vector &slip_bc,
                                         int interior_face_count) const;
   void AssembleSlipContributionBR2Shared(Vector &rhs, const Vector &slip_bc,
                                          int interior_face_count) const;
};

// ============================================================================
// Implementation
// ============================================================================

template <typename MeshType>
AntiplaneBdrLoadOperator<MeshType>::AntiplaneBdrLoadOperator(
   MeshType &mesh, int order, real_t mu, real_t Vp, real_t Wf, DGMethod method,
   int fault_tag)
   : mesh_(mesh), order_(order), mu_(mu), Vp_(Vp), Wf_(Wf),
     fault_tag_(fault_tag),
     method_(method),
     fault_depths_computed_(false),
     mass_inv_computed_(false)
{
   sigma_ = -1.0;  // SIPG
   br2_penalty_ = 3.0;  // D+1 for 2D

   SetupBoundaryMarkers();
   SetupFESpace();
   SetupFaultInfo();
   SetupSolver();

   // Precompute fault depths eagerly (same reason as base operator)
   {
      Vector depths_tmp;
      GetFaultDepths(depths_tmp);
   }

   if (method_ == DGMethod::BR2)
   {
      PrecomputeMassMatrixInverses();
   }
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::SetupBoundaryMarkers()
{
   // Set up far-field boundary markers for Dirichlet BC (left + right)
   // Matches Tandem's BC: Physical Curve(5) = Dirichlet on far-field
   int num_bdr = mesh_.bdr_attributes.Max();
   farfield_bdr_.SetSize(num_bdr);
   farfield_bdr_ = 0;
   farfield_bdr_[BP2BoundaryAttributes::FARFIELD_LEFT - 1] = 1;
   farfield_bdr_[BP2BoundaryAttributes::FARFIELD_RIGHT - 1] = 1;
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::SetupFESpace()
{
   int dim = mesh_.Dimension();
   fec_ = std::make_unique<DG_FECollection>(order_, dim, BasisType::GaussLobatto);
   fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get());
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::GetFaceCenter(
   int face, real_t &x, real_t &z) const
{
   FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
   if (FTr == nullptr) { return; }

   IntegrationPoint ip;
   ip.x = 0.5;

   FTr->SetAllIntPoints(&ip);
   Vector coords(mesh_.Dimension());
   FTr->Face->Transform(ip, coords);

   x = coords(0);
   z = coords(1);
}

template <typename MeshType>
bool AntiplaneBdrLoadOperator<MeshType>::IsFaultFace(int face) const
{
   real_t x, z;
   GetFaceCenter(face, x, z);

   // All interior faces at x ~ 0 are fault faces (no Wf cutoff)
   const real_t tol = 1e-10 * std::max(Wf_, 1.0);
   return std::abs(x) < tol;
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::GetFaceCenterShared(
   int shared_face, real_t &x, real_t &z) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      FaceElementTransformations *FTr =
         mesh_.GetSharedFaceTransformations(shared_face);
      if (FTr == nullptr) { return; }

      IntegrationPoint ip;
      ip.x = 0.5;

      FTr->SetAllIntPoints(&ip);
      Vector coords(mesh_.Dimension());
      FTr->Face->Transform(ip, coords);

      x = coords(0);
      z = coords(1);
#endif
   }
}

template <typename MeshType>
bool AntiplaneBdrLoadOperator<MeshType>::IsFaultFaceShared(int shared_face) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      real_t x, z;
      GetFaceCenterShared(shared_face, x, z);
      const real_t tol = 1e-10 * std::max(Wf_, 1.0);
      return std::abs(x) < tol;
#endif
   }
   return false;
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::SetupFaultInfo()
{
   fault_interior_faces_.SetSize(0);

   if (fault_tag_ >= 1)
   {
      // Tag-based detection: scan boundary elements for the fault tag
      fault_interior_faces_ =
         FaultBoundaryData::FindFaultInteriorFaces(mesh_, fault_tag_);
   }
   else
   {
      // Legacy coordinate-based detection
      int num_faces = mesh_.GetNumFaces();
      for (int f = 0; f < num_faces; f++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (FTr == nullptr) { continue; }

         if (IsFaultFace(f))
         {
            fault_interior_faces_.Append(f);
         }
      }
   }

   // Shared faces: coordinate-based for both paths (see antiplane_operator.hpp)
   fault_shared_faces_.SetSize(0);
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      mesh_.ExchangeFaceNbrData();
      fes_->ExchangeFaceNbrData();

      for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
      {
         if (IsFaultFaceShared(sf))
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
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::SetupSolver()
{
   X_.SetSize(fes_->GetTrueVSize());
   B_.SetSize(fes_->GetTrueVSize());

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // Parallel solver created in AssembleStiffness()
#endif
   }
   else
   {
      auto *cg = new CGSolver();
      cg->SetRelTol(1e-12);
      cg->SetAbsTol(0.0);
      cg->SetMaxIter(5000);
      cg->SetPrintLevel(0);
      solver_.reset(cg);
   }
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::GetFaultDepths(Vector &depths) const
{
   if (!fault_depths_computed_)
   {
      fault_depths_.SetSize(num_fault_dofs_);

      int idx = 0;
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         IntegrationPoint ip;
         ip.x = 0.5;
         FTr->SetAllIntPoints(&ip);
         Vector coords(mesh_.Dimension());
         FTr->Face->Transform(ip, coords);
         fault_depths_(idx++) = coords(1);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            IntegrationPoint ip;
            ip.x = 0.5;
            FTr->SetAllIntPoints(&ip);
            Vector coords(mesh_.Dimension());
            FTr->Face->Transform(ip, coords);
            fault_depths_(idx++) = coords(1);
         }
#endif
      }

      fault_depths_computed_ = true;
   }

   depths = fault_depths_;
}

// ============================================================================
// Stiffness Assembly — adds far-field boundary DG face penalty
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::AssembleStiffness() const
{
   cached_a_ = std::make_unique<BilinFormType>(fes_.get());

   ConstantCoefficient one(1.0);
   cached_a_->AddDomainIntegrator(new DiffusionIntegrator(one));

   if (method_ == DGMethod::IP)
   {
      real_t rep_kappa = (order_ + 1) * (order_ + 1);
      cached_a_->AddInteriorFaceIntegrator(
         new DGDiffusionIntegrator(one, sigma_, rep_kappa));
   }
   else  // BR2
   {
      if (!mass_inv_computed_)
      {
         PrecomputeMassMatrixInverses();
      }
      cached_a_->AddInteriorFaceIntegrator(
         new BR2InteriorFaceIntegrator(one, sigma_, elem_mass_inv_,
                                       mesh_.Dimension()));
   }

   // *** KEY DIFFERENCE: Add far-field boundary face penalty for Dirichlet BC ***
   // Uses IP-style penalty on boundary faces (same as SolveMMS pattern)
   real_t bdr_kappa = (order_ + 1) * (order_ + 1);
   cached_a_->AddBdrFaceIntegrator(
      new DGDiffusionIntegrator(one, sigma_, bdr_kappa), farfield_bdr_);

   cached_a_->Assemble();
   cached_a_->Finalize();

   // Set up solver
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      cached_Ah_.SetType(Operator::Hypre_ParCSR);
      cached_a_->ParallelAssemble(cached_Ah_);

#ifdef MFEM_USE_MUMPS
      auto *mumps = new MUMPSSolver(mesh_.GetComm());
      mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
      mumps->SetPrintLevel(0);
      // Note: For large-scale SEAS problems (>1M DOFs), MUMPS may need
      // ICNTL(14) > 20% to avoid INFOG(1)=-9 workspace errors. The
      // MUMPSSolver API does not currently expose per-instance ICNTL
      // overrides; if hits occur, increase the default in linalg/mumps.cpp.
      mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
      solver_.reset(mumps);
#else
      auto *cg = new CGSolver(mesh_.GetComm());
      cg->SetRelTol(1e-12);
      cg->SetAbsTol(0.0);
      cg->SetMaxIter(5000);
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

// ============================================================================
// Solve — adds far-field Dirichlet RHS contribution
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::Solve(
   real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
   if (!stiffness_assembled_)
   {
      AssembleStiffness();
   }

   // Build RHS
   LinFormType b(fes_.get());

   // *** KEY DIFFERENCE: Add far-field Dirichlet BC RHS contribution ***
   FarFieldLoadingCoefficient bc_coeff(Vp_);
   bc_coeff.SetTime(time);
   ConstantCoefficient mu_coeff(mu_);
   // Note: DGDirichletLFIntegrator uses mu_coeff for the diffusion coefficient.
   // Since the stiffness matrix uses coefficient=1.0 (dimensionless), the
   // Dirichlet LF integrator also uses coefficient=1.0 for consistency.
   ConstantCoefficient one(1.0);
   real_t bdr_kappa = (order_ + 1) * (order_ + 1);
   b.AddBdrFaceIntegrator(
      new DGDirichletLFIntegrator(bc_coeff, one, sigma_, bdr_kappa),
      farfield_bdr_);

   b.Assemble();

   Vector &rhs = b;

   // Add slip contributions (same as base operator)
   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);
   }
   else
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      int interior_count = fault_interior_faces_.Size();
      if (method_ == DGMethod::IP)
      {
         AssembleSlipContributionIPShared(rhs, slip_bc, interior_count);
      }
      else
      {
         AssembleSlipContributionBR2Shared(rhs, slip_bc, interior_count);
      }
#endif
   }

   X_ = 0.0;
   B_ = rhs;

   if (!std::isfinite(B_.Norml2()))
   {
      mfem::err << "ERROR: RHS has NaN/Inf before solve (bdrload operator).\n";
   }

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

   if (!std::isfinite(X_.Norml2()))
   {
      mfem::err << "WARNING: Solver produced NaN/Inf solution.\n";
      X_ = 0.0;
   }

   displacement = X_;
}

// ============================================================================
// InterpolateSlipBC — same as base operator
// ============================================================================

template <typename MeshType>
real_t AntiplaneBdrLoadOperator<MeshType>::InterpolateSlipBC(
   real_t z, const Vector &slip_bc) const
{
   if (num_fault_dofs_ == 0) { return 0.0; }
   if (slip_bc.Size() == 0) { return 0.0; }

   Vector depths;
   GetFaultDepths(depths);

   int n = depths.Size();

   int idx_lo = -1, idx_hi = -1;
   real_t z_lo = -std::numeric_limits<real_t>::max();
   real_t z_hi = std::numeric_limits<real_t>::max();

   int closest_idx = 0;
   real_t min_dist = std::abs(depths(0) - z);

   for (int i = 0; i < n; i++)
   {
      real_t di = depths(i);
      real_t dist = std::abs(di - z);
      if (dist < min_dist)
      {
         min_dist = dist;
         closest_idx = i;
      }
      if (di <= z && di > z_lo)
      {
         z_lo = di;
         idx_lo = i;
      }
      if (di >= z && di < z_hi)
      {
         z_hi = di;
         idx_hi = i;
      }
   }

   if (idx_lo >= 0 && idx_hi >= 0 && idx_lo != idx_hi &&
       idx_lo < slip_bc.Size() && idx_hi < slip_bc.Size())
   {
      real_t dz = z_hi - z_lo;
      if (std::abs(dz) > 1e-30)
      {
         real_t t = (z - z_lo) / dz;
         return (1.0 - t) * slip_bc(idx_lo) + t * slip_bc(idx_hi);
      }
   }

   if (closest_idx < slip_bc.Size())
   {
      return slip_bc(closest_idx);
   }
   return 0.0;
}

// ============================================================================
// ComputeTraction — identical to base operator
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement, const Vector &slip_bc,
   Vector &traction, Vector *normal_traction)
{
   // Antiplane has no normal traction component; ignore the parameter.
   (void)normal_traction;
   traction.SetSize(num_fault_dofs_);
   traction = 0.0;

   GradientGridFunctionCoefficient grad_u(&displacement);

   int idx = 0;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      IntegrationPoint ip;
      ip.x = 0.5;

      FTr->SetAllIntPoints(&ip);
      const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
      const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

      Vector grad1(mesh_.Dimension()), grad2(mesh_.Dimension());
      grad_u.Eval(grad1, *FTr->Elem1, eip1);
      grad_u.Eval(grad2, *FTr->Elem2, eip2);

      real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));
      traction(idx++) = mu_ * avg_dudx;
   }

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

         IntegrationPoint ip;
         ip.x = 0.5;

         FTr->SetAllIntPoints(&ip);
         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         // Elem1 gradient (local)
         const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
         Vector grad1(mesh_.Dimension());
         {
            DenseMatrix dshape(fe1->GetDof(), mesh_.Dimension());
            fe1->CalcDShape(eip1, dshape);

            DenseMatrix Jinv(mesh_.Dimension());
            CalcInverse(FTr->Elem1->Jacobian(), Jinv);

            DenseMatrix dshape_phys(fe1->GetDof(), mesh_.Dimension());
            Mult(dshape, Jinv, dshape_phys);

            Array<int> dofs1;
            fes_->GetElementDofs(FTr->Elem1No, dofs1);
            Vector u1(fe1->GetDof());
            par_u.GetSubVector(dofs1, u1);

            grad1 = 0.0;
            for (int k = 0; k < fe1->GetDof(); k++)
            {
               for (int d = 0; d < mesh_.Dimension(); d++)
               {
                  grad1(d) += dshape_phys(k, d) * u1(k);
               }
            }
         }

         // Elem2 gradient (face-neighbor)
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(
            FTr->Elem2No - mesh_.GetNE());
         Vector grad2(mesh_.Dimension());
         {
            DenseMatrix dshape(fe2->GetDof(), mesh_.Dimension());
            fe2->CalcDShape(eip2, dshape);

            DenseMatrix Jinv(mesh_.Dimension());
            CalcInverse(FTr->Elem2->Jacobian(), Jinv);

            DenseMatrix dshape_phys(fe2->GetDof(), mesh_.Dimension());
            Mult(dshape, Jinv, dshape_phys);

            Vector u2(fe2->GetDof());
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            Array<int> dofs2;
            pfes->GetFaceNbrElementVDofs(nbr_idx, dofs2);
            const Vector &nbr_data = par_u.FaceNbrData();
            for (int k = 0; k < fe2->GetDof(); k++)
            {
               u2(k) = nbr_data(dofs2[k]);
            }

            grad2 = 0.0;
            for (int k = 0; k < fe2->GetDof(); k++)
            {
               for (int d = 0; d < mesh_.Dimension(); d++)
               {
                  grad2(d) += dshape_phys(k, d) * u2(k);
               }
            }
         }

         real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));
         traction(idx++) = mu_ * avg_dudx;
      }
#endif
   }
}

// ============================================================================
// BR2 Lifting Implementation — same as base operator
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::PrecomputeMassMatrixInverses() const
{
   if (mass_inv_computed_) { return; }

   int ne = mesh_.GetNE();
   elem_mass_inv_.resize(ne);

   for (int i = 0; i < ne; i++)
   {
      const FiniteElement *fe = fes_->GetFE(i);
      ElementTransformation *T = mesh_.GetElementTransformation(i);

      int ndof = fe->GetDof();
      DenseMatrix M(ndof);

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());
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

      elem_mass_inv_[i].SetSize(ndof);
      DenseMatrixInverse M_inv(M);
      M_inv.GetInverseMatrix(elem_mass_inv_[i]);
   }

#ifdef MFEM_USE_MPI
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
         ElementTransformation *T = pmesh->GetFaceNbrElementTransformation(i);

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

template <typename MeshType>
const DenseMatrix& AntiplaneBdrLoadOperator<MeshType>::GetMassMatrixInverse(
   int elem) const
{
   if (!mass_inv_computed_)
   {
      PrecomputeMassMatrixInverses();
   }
   return elem_mass_inv_[elem];
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::ComputeLiftingOperator(
   int elem, FaceElementTransformations &FTr, real_t jump, Vector &lift) const
{
   const FiniteElement *fe = fes_->GetFE(elem);
   int ndof = fe->GetDof();
   lift.SetSize(ndof);
   lift = 0.0;

   bool is_elem1 = (elem == FTr.Elem1No);

   int face_order = fe->GetOrder();
   const IntegrationRule &ir = IntRules.Get(FTr.FaceGeom, 2 * face_order + 1);

   Vector shape(ndof);
   Vector int_result(ndof);
   int_result = 0.0;

   for (int p = 0; p < ir.GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir.IntPoint(p);
      FTr.SetAllIntPoints(&ip);

      const IntegrationPoint &eip = is_elem1 ?
                                    FTr.GetElement1IntPoint() : FTr.GetElement2IntPoint();

      fe->CalcShape(eip, shape);

      real_t face_weight = ip.weight * FTr.Face->Weight();
      real_t factor = (FTr.Elem2No >= 0) ? 0.5 : 1.0;
      for (int k = 0; k < ndof; k++)
      {
         int_result(k) += factor * shape(k) * jump * face_weight;
      }
   }

   const DenseMatrix &Minv = GetMassMatrixInverse(elem);
   Minv.Mult(int_result, lift);
}

// ============================================================================
// Slip Assembly — identical to base operator
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::AssembleSlipContributionIP(
   Vector &rhs, const Vector &slip_bc) const
{
   real_t kappa = (order_ + 1) * (order_ + 1);

   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      Array<int> dofs1, dofs2;
      fes_->GetElementDofs(FTr->Elem1No, dofs1);
      fes_->GetElementDofs(FTr->Elem2No, dofs2);

      const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
      const FiniteElement *fe2 = fes_->GetFE(FTr->Elem2No);

      int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
      const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);

      Vector elvec1(dofs1.Size()), elvec2(dofs2.Size());
      elvec1 = 0.0;
      elvec2 = 0.0;

      real_t slip_phys = slip_bc(i);
      if (std::abs(slip_phys) < 1e-15) { continue; }

      for (int p = 0; p < ir.GetNPoints(); p++)
      {
         const IntegrationPoint &ip = ir.IntPoint(p);
         FTr->SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         Vector nor(mesh_.Dimension());
         CalcOrtho(FTr->Jacobian(), nor);

         real_t slip_imposed = (nor(0) > 0) ? -slip_phys : slip_phys;

         Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
         fe1->CalcShape(eip1, shape1);
         fe2->CalcShape(eip2, shape2);

         DenseMatrix dshape1(fe1->GetDof(), mesh_.Dimension());
         DenseMatrix dshape2(fe2->GetDof(), mesh_.Dimension());
         fe1->CalcDShape(eip1, dshape1);
         fe2->CalcDShape(eip2, dshape2);

         DenseMatrix adjJ1(mesh_.Dimension()), adjJ2(mesh_.Dimension());
         CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);
         CalcAdjugate(FTr->Elem2->Jacobian(), adjJ2);

         DenseMatrix dshape1_adj(fe1->GetDof(), mesh_.Dimension());
         DenseMatrix dshape2_adj(fe2->GetDof(), mesh_.Dimension());
         Mult(dshape1, adjJ1, dshape1_adj);
         Mult(dshape2, adjJ2, dshape2_adj);

         Vector dn1(fe1->GetDof()), dn2(fe2->GetDof());
         dshape1_adj.Mult(nor, dn1);
         dshape2_adj.Mult(nor, dn2);

         real_t detJ1 = FTr->Elem1->Weight();
         real_t detJ2 = FTr->Elem2->Weight();
         real_t w1 = ip.weight / (2.0 * detJ1);
         real_t w2 = ip.weight / (2.0 * detJ2);
         real_t nor_sq = nor * nor;
         real_t wq_penalty = kappa * nor_sq * (w1 + w2);

         for (int k = 0; k < dofs1.Size(); k++)
         {
            elvec1(k) += sigma_ * dn1(k) * slip_imposed * w1;
            elvec1(k) += wq_penalty * slip_imposed * shape1(k);
         }

         for (int k = 0; k < dofs2.Size(); k++)
         {
            elvec2(k) += sigma_ * dn2(k) * slip_imposed * w2;
            elvec2(k) -= wq_penalty * slip_imposed * shape2(k);
         }
      }

      for (int k = 0; k < dofs1.Size(); k++)
      {
         rhs(dofs1[k]) += elvec1(k);
      }
      for (int k = 0; k < dofs2.Size(); k++)
      {
         rhs(dofs2[k]) += elvec2(k);
      }
   }
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::AssembleSlipContributionBR2(
   Vector &rhs, const Vector &slip_bc) const
{
   int dim = mesh_.Dimension();

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      int face = fault_interior_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      Array<int> dofs1, dofs2;
      fes_->GetElementDofs(FTr->Elem1No, dofs1);
      fes_->GetElementDofs(FTr->Elem2No, dofs2);

      const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
      const FiniteElement *fe2 = fes_->GetFE(FTr->Elem2No);

      int ndof1 = fe1->GetDof();
      int ndof2 = fe2->GetDof();

      int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
      const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);
      int nqp = ir.GetNPoints();

      const DenseMatrix &Minv1 = GetMassMatrixInverse(FTr->Elem1No);
      const DenseMatrix &Minv2 = GetMassMatrixInverse(FTr->Elem2No);

      real_t slip_phys = slip_bc(fi);
      if (std::abs(slip_phys) < 1e-15) { continue; }

      DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
      Vector nor_all(dim * nqp);
      Vector w_all(nqp);
      Vector slip_all(nqp);

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         FTr->SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         Vector shape1_q(shapes1.GetColumn(q), ndof1);
         Vector shape2_q(shapes2.GetColumn(q), ndof2);
         fe1->CalcShape(eip1, shape1_q);
         fe2->CalcShape(eip2, shape2_q);

         Vector nor_q(&nor_all[q * dim], dim);
         CalcOrtho(FTr->Jacobian(), nor_q);

         w_all[q] = ip.weight;
         slip_all[q] = (nor_q(0) > 0) ? -slip_phys : slip_phys;
      }

      real_t max_slip = 0.0;
      for (int q = 0; q < nqp; q++)
      {
         max_slip = std::max(max_slip, std::abs(slip_all[q]));
      }
      if (max_slip < 1e-15) { continue; }

      // Compute BR2 lifted slip
      DenseMatrix f_lifted1(ndof1, dim), f_lifted2(ndof2, dim);
      f_lifted1 = 0.0;
      f_lifted2 = 0.0;

      DenseMatrix face_int1(ndof1, dim), face_int2(ndof2, dim);
      face_int1 = 0.0;
      face_int2 = 0.0;

      for (int q = 0; q < nqp; q++)
      {
         real_t wq = w_all[q];
         real_t slip_q = slip_all[q];

         for (int j = 0; j < dim; j++)
         {
            real_t n_j = nor_all[q * dim + j];
            real_t factor = slip_q * n_j * wq;

            for (int m = 0; m < ndof1; m++)
            {
               face_int1(m, j) += shapes1(m, q) * factor;
            }
            for (int m = 0; m < ndof2; m++)
            {
               face_int2(m, j) += shapes2(m, q) * factor;
            }
         }
      }

      for (int j = 0; j < dim; j++)
      {
         Vector fi1_j(ndof1), fi2_j(ndof2);
         Vector fl1_j(ndof1), fl2_j(ndof2);

         face_int1.GetColumn(j, fi1_j);
         face_int2.GetColumn(j, fi2_j);

         Minv1.Mult(fi1_j, fl1_j);
         Minv2.Mult(fi2_j, fl2_j);

         fl1_j *= 0.5;
         fl2_j *= 0.5;

         f_lifted1.SetCol(j, fl1_j);
         f_lifted2.SetCol(j, fl2_j);
      }

      Vector f_lifted_q(nqp);
      for (int q = 0; q < nqp; q++)
      {
         real_t sum = 0.0;
         for (int j = 0; j < dim; j++)
         {
            real_t n_j = nor_all[q * dim + j];
            real_t contrib1 = 0.0, contrib2 = 0.0;

            for (int l = 0; l < ndof1; l++)
            {
               contrib1 += shapes1(l, q) * f_lifted1(l, j);
            }
            for (int l = 0; l < ndof2; l++)
            {
               contrib2 += shapes2(l, q) * f_lifted2(l, j);
            }

            sum += n_j * (contrib1 + contrib2);
         }
         f_lifted_q[q] = 0.5 * sum;
      }

      real_t penalty = br2_penalty_;

      Vector elvec1(ndof1), elvec2(ndof2);
      elvec1 = 0.0;
      elvec2 = 0.0;

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         FTr->SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         real_t wq = w_all[q];
         real_t slip_q = slip_all[q];

         if (std::abs(slip_q) < 1e-15) { continue; }

         Vector nor_q(dim);
         for (int j = 0; j < dim; j++)
         {
            nor_q[j] = nor_all[q * dim + j];
         }

         DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
         fe1->CalcDShape(eip1, dshape1_ref);
         fe2->CalcDShape(eip2, dshape2_ref);

         DenseMatrix adjJ1(dim), adjJ2(dim);
         CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);
         CalcAdjugate(FTr->Elem2->Jacobian(), adjJ2);

         DenseMatrix dshape1_adj(ndof1, dim), dshape2_adj(ndof2, dim);
         Mult(dshape1_ref, adjJ1, dshape1_adj);
         Mult(dshape2_ref, adjJ2, dshape2_adj);

         Vector dn1(ndof1), dn2(ndof2);
         dshape1_adj.Mult(nor_q, dn1);
         dshape2_adj.Mult(nor_q, dn2);

         real_t detJ1 = FTr->Elem1->Weight();
         real_t detJ2 = FTr->Elem2->Weight();

         real_t c1 = sigma_ * 0.5;

         for (int k = 0; k < ndof1; k++)
         {
            elvec1(k) += c1 * wq * dn1(k) * slip_q / detJ1;
         }
         for (int k = 0; k < ndof2; k++)
         {
            elvec2(k) += c1 * wq * dn2(k) * slip_q / detJ2;
         }

         for (int k = 0; k < ndof1; k++)
         {
            elvec1(k) += penalty * wq * shapes1(k, q) * f_lifted_q[q];
         }
         for (int k = 0; k < ndof2; k++)
         {
            elvec2(k) += (-penalty) * wq * shapes2(k, q) * f_lifted_q[q];
         }
      }

      for (int k = 0; k < ndof1; k++)
      {
         rhs(dofs1[k]) += elvec1(k);
      }
      for (int k = 0; k < ndof2; k++)
      {
         rhs(dofs2[k]) += elvec2(k);
      }
   }
}

// ============================================================================
// Shared Face Slip Assembly — identical to base operator
// ============================================================================

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::AssembleSlipContributionIPShared(
   Vector &rhs, const Vector &slip_bc, int interior_face_count) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      real_t kappa = (order_ + 1) * (order_ + 1);

      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }

         int slip_idx = interior_face_count + i;
         real_t slip_phys = slip_bc(slip_idx);
         if (std::abs(slip_phys) < 1e-15) { continue; }

         Array<int> dofs1;
         fes_->GetElementDofs(FTr->Elem1No, dofs1);

         const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
         int face_order = fe1->GetOrder();
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);

         Vector elvec1(dofs1.Size());
         elvec1 = 0.0;

         for (int p = 0; p < ir.GetNPoints(); p++)
         {
            const IntegrationPoint &ip = ir.IntPoint(p);
            FTr->SetAllIntPoints(&ip);

            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

            Vector nor(mesh_.Dimension());
            CalcOrtho(FTr->Jacobian(), nor);

            real_t slip_imposed = (nor(0) > 0) ? -slip_phys : slip_phys;

            Vector shape1(fe1->GetDof());
            fe1->CalcShape(eip1, shape1);

            DenseMatrix dshape1(fe1->GetDof(), mesh_.Dimension());
            fe1->CalcDShape(eip1, dshape1);

            DenseMatrix adjJ1(mesh_.Dimension());
            CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);

            DenseMatrix dshape1_adj(fe1->GetDof(), mesh_.Dimension());
            Mult(dshape1, adjJ1, dshape1_adj);

            Vector dn1(fe1->GetDof());
            dshape1_adj.Mult(nor, dn1);

            real_t detJ1 = FTr->Elem1->Weight();
            real_t detJ2 = FTr->Elem2->Weight();
            real_t w1 = ip.weight / (2.0 * detJ1);
            real_t nor_sq = nor * nor;
            real_t w2 = ip.weight / (2.0 * detJ2);
            real_t wq_penalty = kappa * nor_sq * (w1 + w2);

            for (int k = 0; k < dofs1.Size(); k++)
            {
               elvec1(k) += sigma_ * dn1(k) * slip_imposed * w1;
               elvec1(k) += wq_penalty * slip_imposed * shape1(k);
            }
         }

         for (int k = 0; k < dofs1.Size(); k++)
         {
            rhs(dofs1[k]) += elvec1(k);
         }
      }
#endif
   }
}

template <typename MeshType>
void AntiplaneBdrLoadOperator<MeshType>::AssembleSlipContributionBR2Shared(
   Vector &rhs, const Vector &slip_bc, int interior_face_count) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      int dim = mesh_.Dimension();
      real_t penalty = br2_penalty_;

      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

      for (int fi = 0; fi < fault_shared_faces_.Size(); fi++)
      {
         int sf = fault_shared_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }

         int slip_idx = interior_face_count + fi;
         real_t slip_phys = slip_bc(slip_idx);
         if (std::abs(slip_phys) < 1e-15) { continue; }

         Array<int> dofs1;
         fes_->GetElementDofs(FTr->Elem1No, dofs1);

         const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
         int ndof1 = fe1->GetDof();

         int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
         int ndof2 = fe2->GetDof();

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);
         int nqp = ir.GetNPoints();

         const DenseMatrix &Minv1 = GetMassMatrixInverse(FTr->Elem1No);
         const DenseMatrix &Minv2 = GetMassMatrixInverse(FTr->Elem2No);

         DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
         Vector nor_all(dim * nqp);
         Vector w_all(nqp);
         Vector slip_all(nqp);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            Vector shape1_q(shapes1.GetColumn(q), ndof1);
            fe1->CalcShape(eip1, shape1_q);

            Vector shape2_q(shapes2.GetColumn(q), ndof2);
            fe2->CalcShape(eip2, shape2_q);

            Vector nor_q(&nor_all[q * dim], dim);
            CalcOrtho(FTr->Jacobian(), nor_q);

            w_all[q] = ip.weight;
            slip_all[q] = (nor_q(0) > 0) ? -slip_phys : slip_phys;
         }

         DenseMatrix f_lifted1(ndof1, dim), f_lifted2(ndof2, dim);
         f_lifted1 = 0.0;
         f_lifted2 = 0.0;

         DenseMatrix face_int1(ndof1, dim), face_int2(ndof2, dim);
         face_int1 = 0.0;
         face_int2 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            for (int j = 0; j < dim; j++)
            {
               real_t n_j = nor_all[q * dim + j];
               real_t factor = slip_all[q] * n_j * w_all[q];
               for (int m = 0; m < ndof1; m++)
               {
                  face_int1(m, j) += shapes1(m, q) * factor;
               }
               for (int m = 0; m < ndof2; m++)
               {
                  face_int2(m, j) += shapes2(m, q) * factor;
               }
            }
         }

         for (int j = 0; j < dim; j++)
         {
            Vector fi1_j(ndof1), fi2_j(ndof2);
            Vector fl1_j(ndof1), fl2_j(ndof2);

            face_int1.GetColumn(j, fi1_j);
            face_int2.GetColumn(j, fi2_j);

            Minv1.Mult(fi1_j, fl1_j);
            Minv2.Mult(fi2_j, fl2_j);

            fl1_j *= 0.5;
            fl2_j *= 0.5;

            f_lifted1.SetCol(j, fl1_j);
            f_lifted2.SetCol(j, fl2_j);
         }

         Vector f_lifted_q(nqp);
         for (int q = 0; q < nqp; q++)
         {
            real_t sum = 0.0;
            for (int j = 0; j < dim; j++)
            {
               real_t n_j = nor_all[q * dim + j];
               real_t contrib1 = 0.0, contrib2 = 0.0;
               for (int l = 0; l < ndof1; l++)
               {
                  contrib1 += shapes1(l, q) * f_lifted1(l, j);
               }
               for (int l = 0; l < ndof2; l++)
               {
                  contrib2 += shapes2(l, q) * f_lifted2(l, j);
               }
               sum += n_j * (contrib1 + contrib2);
            }
            f_lifted_q[q] = 0.5 * sum;
         }

         Vector elvec1(ndof1);
         elvec1 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

            real_t wq = w_all[q];
            real_t slip_q = slip_all[q];
            if (std::abs(slip_q) < 1e-15) { continue; }

            Vector nor_q(dim);
            for (int j = 0; j < dim; j++)
            {
               nor_q[j] = nor_all[q * dim + j];
            }

            DenseMatrix dshape1_ref(ndof1, dim);
            fe1->CalcDShape(eip1, dshape1_ref);

            DenseMatrix adjJ1(dim);
            CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);

            DenseMatrix dshape1_adj(ndof1, dim);
            Mult(dshape1_ref, adjJ1, dshape1_adj);

            Vector dn1(ndof1);
            dshape1_adj.Mult(nor_q, dn1);

            real_t detJ1 = FTr->Elem1->Weight();
            real_t c1 = sigma_ * 0.5;

            for (int k = 0; k < ndof1; k++)
            {
               elvec1(k) += c1 * wq * dn1(k) * slip_q / detJ1;
               elvec1(k) += penalty * wq * shapes1(k, q) * f_lifted_q[q];
            }
         }

         for (int k = 0; k < ndof1; k++)
         {
            rhs(dofs1[k]) += elvec1(k);
         }
      }
#endif
   }
}

// Convenience type aliases
using SerialBdrLoadOperator = AntiplaneBdrLoadOperator<Mesh>;
#ifdef MFEM_USE_MPI
using ParallelBdrLoadOperator = AntiplaneBdrLoadOperator<ParMesh>;
#endif

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ANTIPLANE_BDRLOAD_OPERATOR_HPP
