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

#ifndef MFEM_SEAS_ANTIPLANE_OPERATOR_HPP
#define MFEM_SEAS_ANTIPLANE_OPERATOR_HPP

#include "mfem.hpp"
#include "domain_operator.hpp"
#include "bp2_mesh.hpp"  // For BP2BoundaryAttributes
#include "../common/seas_types.hpp"
#include "../integrator/dg_br2_integrator.hpp"

#include <memory>

namespace mfem
{
namespace seas
{

/// @brief DG method selection for interior penalty formulations
///
/// Following Tandem implementation (see /Users/chunhuizhao/projects/tandem/src/form/DGCurvilinearCommon.h)
enum class DGMethod
{
   IP,   ///< Interior Penalty (SIPG) with per-face penalty following Tandem
   BR2   ///< Bassi-Rebay 2 with lifting operator
};

/// @brief DG Antiplane domain operator for 2D scalar Laplace equation (FULL DOMAIN)
///
/// Solves the antiplane shear problem using Discontinuous Galerkin (DG) method:
///   ∇²u = ∂²u/∂x² + ∂²u/∂z² = 0
///
/// This follows Tandem's full domain approach:
/// - Full domain: x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
/// - z = 0 is the free surface, z < 0 represents greater depth
/// - Fault at x = 0 is an INTERIOR interface (not a boundary)
/// - Slip is imposed as a jump condition on interior faces: [[u]] = slip
///
/// Boundary conditions:
///   - Far-field left (x=-Lx): Natural BC (zero traction)
///   - Far-field right (x=+Lx): Natural BC (zero traction)
///   - Free surface (z=0): Natural BC (zero traction)
///   - Bottom (z=-Lz): Natural BC (zero traction)
///
/// Tectonic loading: Prescribed slip rate Vp on fault (x=0) below Wf,
/// following SCEC SEAS benchmark specification.
///
/// In DG formulation, Dirichlet BCs are imposed weakly through face integrators:
/// - BilinearForm: AddBdrFaceIntegrator(DGDiffusionIntegrator)
/// - LinearForm: AddBdrFaceIntegrator(DGDirichletLFIntegrator)
///
/// The traction on the fault is computed from the interior faces at x=0:
///   τ = μ * {{∂u/∂x}} (average of x-gradient from both sides)
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class AntiplaneDomainOperator : public DomainOperator<MeshType>
{
public:
   // Inherit type aliases from base
   using Base = DomainOperator<MeshType>;
   using FESpaceType = typename Base::FESpaceType;
   using GridFuncType = typename Base::GridFuncType;
   using BilinFormType = typename Base::BilinFormType;
   using LinFormType = typename Base::LinFormType;

   /// @brief Construct the DG antiplane operator for full domain
   ///
   /// @param mesh The computational mesh (full domain, owned externally)
   /// @param order Polynomial order for DG finite elements
   /// @param mu Shear modulus [Pa]
   /// @param Vp Plate rate [m/s] for far-field boundary condition
   /// @param Wf Fault depth [m] - slip applied from z=0 to z=-Wf (positive value)
   /// @param method DG method (IP or BR2), default is IP
   AntiplaneDomainOperator(MeshType &mesh, int order, real_t mu, real_t Vp,
                           real_t Wf = 40.0e3, DGMethod method = DGMethod::IP);

   /// Destructor
   ~AntiplaneDomainOperator() override = default;

   // Interface implementation
   int NumComponents() const override { return 1; }
   int Dimension() const override { return mesh_.Dimension(); }

   void Solve(real_t time, const Vector &slip_bc,
              GridFuncType &displacement) override;

   /// @brief Solve with custom Dirichlet BC (for MMS testing)
   ///
   /// Solves ∇²u = 0 with given Dirichlet BC on the bottom boundary
   /// and zero slip. This allows MMS tests to verify the solver
   /// convergence with an exact solution.
   ///
   /// @param dirichlet_bc Coefficient defining the Dirichlet BC at the bottom
   /// @param displacement Output displacement field
   void SolveMMS(Coefficient &dirichlet_bc, GridFuncType &displacement);

   void ComputeTraction(const GridFuncType &displacement,
                        const Vector &slip_bc,
                        Vector &traction) override;

   FESpaceType &GetFESpace() override { return *fes_; }
   const FESpaceType &GetFESpace() const override { return *fes_; }

   MeshType &GetMesh() override { return mesh_; }
   const MeshType &GetMesh() const override { return mesh_; }

   real_t GetShearModulus() const override { return mu_; }

   int GetNumFaultDOFs() const override { return num_fault_dofs_; }

   void GetFaultDepths(Vector &depths) const override;

   const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

   /// @brief Get the array of fault interior face indices
   const Array<int> &GetFaultInteriorFaces() const { return fault_interior_faces_; }

   /// @brief Get the array of fault shared face indices (parallel only)
   const Array<int> &GetFaultSharedFaces() const { return fault_shared_faces_; }

   /// @brief Get plate rate
   real_t GetPlateRate() const { return Vp_; }

   /// @brief Get fault depth
   real_t GetFaultDepth() const { return Wf_; }

   /// @brief Get the polynomial order
   int GetOrder() const { return order_; }

   /// @brief Get the DG method (IP or BR2)
   DGMethod GetMethod() const { return method_; }

   /// @brief Interpolate slip BC value at a given depth
   ///
   /// This is used by the DG slip jump coefficient to get the slip
   /// value at any point on the fault interface.
   real_t InterpolateSlipBC(real_t z, const Vector &slip_bc) const;

private:
   MeshType &mesh_;
   int order_;
   real_t mu_;   // Shear modulus
   real_t Vp_;   // Plate rate
   real_t Wf_;   // Fault depth

   // DG method selection
   DGMethod method_;

   // DG parameters
   real_t sigma_;  // SIPG sign parameter (-1 for SIPG, +1 for NIPG)

   // For IP method: per-face penalty computed from Tandem's formula
   // For BR2 method: constant penalty sigma = D+1 = 3
   real_t br2_penalty_;  // BR2 penalty (constant = D+1)

   // Precomputed element mass matrix inverses for BR2
   mutable std::vector<DenseMatrix> elem_mass_inv_;
   mutable bool mass_inv_computed_;

   std::unique_ptr<FiniteElementCollection> fec_;
   std::unique_ptr<FESpaceType> fes_;

   // No Dirichlet boundary markers needed - all boundaries have Natural BC

   // Interior fault faces (at x=0)
   Array<int> fault_interior_faces_;

   // Shared fault faces (at x=0, parallel only)
   Array<int> fault_shared_faces_;

   // Fault DOF information (for traction computation and output)
   Array<int> fault_dofs_;
   int num_fault_dofs_;

   // Fault DOF depths (z-coordinates)
   mutable Vector fault_depths_;
   mutable bool fault_depths_computed_;

   // Solver components
   std::unique_ptr<Solver> solver_;
   std::unique_ptr<Solver> prec_;
   mutable std::unique_ptr<GSSmoother> serial_prec_;  // For serial builds - must outlive solve
   mutable Vector X_, B_;

   // Cached stiffness matrix (assembled once, reused across Solve calls)
   mutable bool stiffness_assembled_ = false;
   mutable std::unique_ptr<BilinFormType> cached_a_;
   // Parallel: cached HypreParMatrix and preconditioner
   mutable OperatorHandle cached_Ah_;
   mutable std::unique_ptr<Solver> cached_prec_;

   /// @brief Assemble stiffness matrix (once, on first Solve call)
   void AssembleStiffness() const;

   // Setup methods
   void SetupBoundaryMarkers();
   void SetupFESpace();
   void SetupFaultInfo();
   void SetupSolver();

   /// @brief Check if an interior face is on the fault (x ≈ 0 and z > -Wf)
   ///
   /// Fault extends from z=0 (surface) to z=-Wf (depth)
   bool IsFaultFace(int face) const;

   /// @brief Check if a shared face is on the fault (parallel only)
   bool IsFaultFaceShared(int shared_face) const;

   /// @brief Get the center coordinates of an interior face
   void GetFaceCenter(int face, real_t &x, real_t &z) const;

   /// @brief Get the center coordinates of a shared face (parallel only)
   void GetFaceCenterShared(int shared_face, real_t &x, real_t &z) const;

   // ========================================================================
   // BR2 Lifting Operator
   // ========================================================================

   /// @brief Precompute element mass matrix inverses for BR2
   void PrecomputeMassMatrixInverses() const;

   /// @brief Get element mass matrix inverse (precomputed)
   const DenseMatrix& GetMassMatrixInverse(int elem) const;

   /// @brief Compute BR2 lifting operator for a face
   ///
   /// r_e|_K = M_K^{-1} * ∫_e φ_K ⊗ n [[v]] ds
   /// For interior faces, use factor 0.5 on each side
   /// For boundary faces, use full factor
   ///
   /// @param elem Element index
   /// @param face_idx Local face index in element
   /// @param jump The jump value [[v]] at the face
   /// @param lift Output: lifting operator contribution to element DOFs
   void ComputeLiftingOperator(int elem, FaceElementTransformations &FTr,
                               real_t jump, Vector &lift) const;

   // ========================================================================
   // Method-specific assembly
   // ========================================================================

   /// @brief Assemble slip contribution using IP method
   void AssembleSlipContributionIP(Vector &rhs, const Vector &slip_bc) const;

   /// @brief Assemble slip contribution using BR2 method
   void AssembleSlipContributionBR2(Vector &rhs, const Vector &slip_bc) const;

   /// @brief Assemble slip contribution on shared faces using IP method (parallel only)
   void AssembleSlipContributionIPShared(Vector &rhs, const Vector &slip_bc,
                                         int interior_face_count) const;

   /// @brief Assemble slip contribution on shared faces using BR2 method (parallel only)
   void AssembleSlipContributionBR2Shared(Vector &rhs, const Vector &slip_bc,
                                          int interior_face_count) const;
};

// ============================================================================
// Implementation
// ============================================================================

template <typename MeshType>
AntiplaneDomainOperator<MeshType>::AntiplaneDomainOperator(
   MeshType &mesh, int order, real_t mu, real_t Vp, real_t Wf, DGMethod method)
   : mesh_(mesh), order_(order), mu_(mu), Vp_(Vp), Wf_(Wf),
     method_(method),
     fault_depths_computed_(false),
     mass_inv_computed_(false)
{
   // SIPG sign parameter (always -1 for symmetric interior penalty)
   sigma_ = -1.0;

   // BR2 penalty: σ = D + 1 = 3 for 2D
   // See dg_antiplane_theory.md Section 4.2
   br2_penalty_ = 3.0;

   SetupBoundaryMarkers();
   SetupFESpace();
   SetupFaultInfo();
   SetupSolver();

   // Precompute fault depths eagerly to avoid invalidating
   // FaceElementTransformations pointers during slip assembly.
   // GetFaultDepths calls GetInteriorFaceTransformations which returns a
   // pointer to a shared internal object (Mesh::FaceElemTr). If depths are
   // computed lazily inside AssembleSlipContribution*, the call chain
   // InterpolateSlipBC -> GetFaultDepths -> GetInteriorFaceTransformations
   // clobbers the caller's FTr pointer, corrupting the first fault face.
   {
      Vector depths_tmp;
      GetFaultDepths(depths_tmp);
   }

   // Precompute mass matrix inverses for BR2
   if (method_ == DGMethod::BR2)
   {
      PrecomputeMassMatrixInverses();
   }
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SetupBoundaryMarkers()
{
   // All boundaries have Natural BC (zero traction):
   // - FARFIELD_LEFT (x=-Lx): zero traction at far-field
   // - FARFIELD_RIGHT (x=+Lx): zero traction at far-field
   // - FREE_SURFACE (z=0): zero traction at free surface (Earth's surface)
   // - BOTTOM (z=-Lz): zero traction (loading applied via prescribed Vp on fault)
   //
   // Tectonic loading is applied by prescribing slip rate = Vp on the fault
   // interface (x=0) below the rate-state zone (z < -Wf), following the
   // SCEC SEAS benchmark specification.
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SetupFESpace()
{
   int dim = mesh_.Dimension();

   // Use DG (L2) finite element collection for discontinuous Galerkin
   // This is required by the documentation: "Scalar L2 DG space"
   // see also https://github.com/mfem/mfem/blob/master/examples/ex17.cpp
   fec_ = std::make_unique<DG_FECollection>(order_, dim, BasisType::GaussLobatto);

   // Create finite element space (scalar, vdim=1)
   // Works for both serial (Mesh/FiniteElementSpace) and parallel
   // (ParMesh/ParFiniteElementSpace) via the FESpaceType alias
   fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get());

   // In DG, there are no essential DOFs - all BCs are imposed weakly
   // through face integrators
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::GetFaceCenter(
   int face, real_t &x, real_t &z) const
{
   // Get face transformation
   FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
   if (FTr == nullptr) { return; }

   // Use face centroid
   IntegrationPoint ip;
   ip.x = 0.5;

   FTr->SetAllIntPoints(&ip);
   Vector coords(mesh_.Dimension());
   FTr->Face->Transform(ip, coords);

   x = coords(0);
   z = coords(1);
}

template <typename MeshType>
bool AntiplaneDomainOperator<MeshType>::IsFaultFace(int face) const
{
   real_t x, z;
   GetFaceCenter(face, x, z);

   // Check if face is at x ≈ 0 (full fault interface from z=0 to z=-Lz)
   // Rate-state friction applies for z > -Wf; prescribed Vp loading for z < -Wf
   const real_t tol = 1e-10 * std::max(Wf_, 1.0);
   return std::abs(x) < tol;
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::GetFaceCenterShared(
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
bool AntiplaneDomainOperator<MeshType>::IsFaultFaceShared(int shared_face) const
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
void AntiplaneDomainOperator<MeshType>::SetupFaultInfo()
{
   // In full domain, the fault is an interior interface at x = 0.
   // We identify interior faces that lie on x = 0 within the fault zone (z < Wf).
   fault_interior_faces_.SetSize(0);

   // Iterate over all interior faces
   int num_faces = mesh_.GetNumFaces();
   for (int f = 0; f < num_faces; f++)
   {
      // Check if this is an interior face
      FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }  // Boundary face

      // Check if face is on the fault
      if (IsFaultFace(f))
      {
         fault_interior_faces_.Append(f);
      }
   }

   // Also iterate shared faces in parallel
   fault_shared_faces_.SetSize(0);
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // Must exchange face neighbor data before accessing shared faces
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

   // One fault DOF per face, evaluated at the face midpoint.
   // Includes both interior and shared fault faces.
   num_fault_dofs_ = fault_interior_faces_.Size() + fault_shared_faces_.Size();

   // Create fault DOF array (these are just indices for the fault state vector)
   fault_dofs_.SetSize(num_fault_dofs_);
   for (int i = 0; i < num_fault_dofs_; i++)
   {
      fault_dofs_[i] = i;
   }
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SetupSolver()
{
   // Allocate work vectors
   X_.SetSize(fes_->GetTrueVSize());
   B_.SetSize(fes_->GetTrueVSize());

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // Parallel solver: CGSolver with BoomerAMG preconditioner.
      // The DG penalty/stabilization terms make the system SPD even with
      // all-Neumann BCs, so AMG is safe and much faster than plain smoothing.
      // CG tolerance must be much tighter than ODE atol (1e-7) to avoid
      // solver noise contaminating the RK45 error estimate. Tandem uses 1e-12.
      auto *cg = new CGSolver(mesh_.GetComm());
      cg->SetRelTol(1e-12);
      cg->SetAbsTol(0.0);
      cg->SetMaxIter(2000);
      cg->SetPrintLevel(-1);

      solver_.reset(cg);
#endif
   }
   else
   {
      // Serial solver setup
      auto *cg = new CGSolver();
      cg->SetRelTol(1e-12);
      cg->SetAbsTol(0.0);
      cg->SetMaxIter(2000);
      cg->SetPrintLevel(-1);

      solver_.reset(cg);
   }
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::GetFaultDepths(Vector &depths) const
{
   if (!fault_depths_computed_)
   {
      fault_depths_.SetSize(num_fault_dofs_);

      int idx = 0;
      // Interior fault faces
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

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::AssembleStiffness() const
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

   cached_a_->Assemble();
   cached_a_->Finalize();

   // Set up solver operator and preconditioner
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      cached_Ah_.SetType(Operator::Hypre_ParCSR);
      cached_a_->ParallelAssemble(cached_Ah_);

      auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
      amg->SetPrintLevel(0);
      cached_prec_.reset(amg);

      auto *cg = static_cast<CGSolver*>(solver_.get());
      cg->SetPreconditioner(*cached_prec_);
      cg->SetOperator(*cached_Ah_.As<HypreParMatrix>());
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

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::Solve(
   real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
   // Assemble stiffness matrix once (it doesn't change between time steps)
   if (!stiffness_assembled_)
   {
      AssembleStiffness();
   }

   // Build RHS (changes every call due to different slip_bc)
   LinFormType b(fes_.get());
   b.Assemble();

   Vector &rhs = b;

   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);
   }
   else  // BR2
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
   }

   // Add shared face slip contributions (parallel only)
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

   // Solve using cached operator
   solver_->Mult(B_, X_);

   // Copy solution to GridFunction
   displacement = X_;
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SolveMMS(
   Coefficient &dirichlet_bc, GridFuncType &displacement)
{
   // For MMS, apply Dirichlet BCs on ALL boundaries so the exact solution
   // is fully determined. The physical operator only has Dirichlet on the
   // bottom, but MMS needs it everywhere for a well-posed verification.

   // Mark all boundaries for Dirichlet BC
   int num_bdr = mesh_.bdr_attributes.Max();
   Array<int> all_bdr_marker(num_bdr);
   all_bdr_marker = 1;

   // Create bilinear form
   BilinFormType a(fes_.get());
   ConstantCoefficient one(1.0);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));

   // For MMS, use IP penalty on boundary faces for both methods to ensure
   // consistency with DGDirichletLFIntegrator on the RHS. The method-specific
   // stabilization (IP vs BR2) is applied on interior faces only.
   real_t bdr_kappa = (order_ + 1) * (order_ + 1);

   if (method_ == DGMethod::IP)
   {
      a.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(one, sigma_, bdr_kappa));
   }
   else  // BR2
   {
      if (!mass_inv_computed_)
      {
         PrecomputeMassMatrixInverses();
      }
      a.AddInteriorFaceIntegrator(
         new BR2InteriorFaceIntegrator(one, sigma_, elem_mass_inv_, mesh_.Dimension()));
   }

   // Boundary face integrator: always use standard DGDiffusionIntegrator
   // to be consistent with DGDirichletLFIntegrator for the Dirichlet RHS
   a.AddBdrFaceIntegrator(
      new DGDiffusionIntegrator(one, sigma_, bdr_kappa), all_bdr_marker);

   a.Assemble();
   a.Finalize();

   // RHS with exact solution as Dirichlet BC on all boundaries, no slip
   LinFormType b(fes_.get());
   b.AddBdrFaceIntegrator(
      new DGDirichletLFIntegrator(dirichlet_bc, one, sigma_, bdr_kappa),
      all_bdr_marker);
   b.Assemble();

   X_ = 0.0;

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // Parallel: assemble to HypreParMatrix
      OperatorHandle Ah;
      Ah.SetType(Operator::Hypre_ParCSR);
      a.ParallelAssemble(Ah);

      B_ = b;

      // Use HypreSmoother as preconditioner for CGSolver
      auto smoother = std::make_unique<HypreSmoother>(
         *Ah.As<HypreParMatrix>());

      auto *cg = static_cast<CGSolver*>(solver_.get());
      cg->SetPreconditioner(*smoother);
      cg->SetOperator(*Ah.As<HypreParMatrix>());
      cg->Mult(B_, X_);
#endif
   }
   else
   {
      SparseMatrix &A = a.SpMat();
      B_ = b;

      serial_prec_ = std::make_unique<GSSmoother>(A);
      static_cast<CGSolver*>(solver_.get())->SetPreconditioner(*serial_prec_);
      solver_->SetOperator(A);
      solver_->Mult(B_, X_);
   }

   displacement = X_;
}

template <typename MeshType>
real_t AntiplaneDomainOperator<MeshType>::InterpolateSlipBC(
   real_t z, const Vector &slip_bc) const
{
   // Linear interpolation of slip values at fault DOF depths.
   // This provides O(h²) accuracy for smooth slip profiles,
   // compared to O(h) for nearest-neighbor.

   if (num_fault_dofs_ == 0) { return 0.0; }
   if (slip_bc.Size() == 0) { return 0.0; }

   // Get fault depths (precomputed)
   Vector depths;
   GetFaultDepths(depths);

   int n = depths.Size();

   // Find the bracketing interval for linear interpolation.
   // Depths may not be sorted, so find the two closest points.
   int idx_lo = -1, idx_hi = -1;
   real_t z_lo = -std::numeric_limits<real_t>::max();
   real_t z_hi = std::numeric_limits<real_t>::max();

   // Also track closest point for fallback (extrapolation)
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

   // If both brackets found and they are distinct, interpolate linearly
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

   // Exact match or single bracket: use the closest value
   if (closest_idx < slip_bc.Size())
   {
      return slip_bc(closest_idx);
   }
   return 0.0;
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement, const Vector &slip_bc,
   Vector &traction)
{
   // DG numerical flux traction on fault faces.
   //
   // Following Tandem's grad_u kernel (poisson.py:147-149):
   //   grad_u = 0.5*(K∇u₁ + K∇u₂) + c0 * ([[u]] - δ) * n̂
   //   traction = μ * grad_u · n̂
   //
   // For IP method: gradient only (penalty residual is too large on coarse mesh)
   // For BR2 method: include penalty with c0 = -penalty = -(D+1) = -3
   //   The BR2 lifting produces a much smaller penalty residual than IP,
   //   making the penalty correction safe even on coarse meshes.

   traction.SetSize(num_fault_dofs_);
   traction = 0.0;

   // Create a gradient coefficient
   GradientGridFunctionCoefficient grad_u(&displacement);

   int idx = 0;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);

      if (FTr == nullptr) { continue; }

      // Single midpoint evaluation, consistent with GetFaultDepths().
      IntegrationPoint ip;
      ip.x = 0.5;

      FTr->SetAllIntPoints(&ip);
      const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
      const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

      // Evaluate gradient in both adjacent elements
      Vector grad1(mesh_.Dimension()), grad2(mesh_.Dimension());
      grad_u.Eval(grad1, *FTr->Elem1, eip1);
      grad_u.Eval(grad2, *FTr->Elem2, eip2);

      // Average gradient (for symmetric flux)
      real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));

      real_t tau_face = mu_ * avg_dudx;

      traction(idx++) = tau_face;
   }

   // Shared fault faces (parallel only)
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // For shared faces, we need face-neighbor data for Elem2 gradient.
      // Exchange face-neighbor data first.
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
         Vector grad1(mesh_.Dimension());
         {
            const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
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
         Vector grad2(mesh_.Dimension());
         {
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(
               FTr->Elem2No - mesh_.GetNE());
            DenseMatrix dshape(fe2->GetDof(), mesh_.Dimension());
            fe2->CalcDShape(eip2, dshape);

            DenseMatrix Jinv(mesh_.Dimension());
            CalcInverse(FTr->Elem2->Jacobian(), Jinv);

            DenseMatrix dshape_phys(fe2->GetDof(), mesh_.Dimension());
            Mult(dshape, Jinv, dshape_phys);

            // Get face-neighbor DOF values
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
// BR2 Lifting Implementation
// ============================================================================

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::PrecomputeMassMatrixInverses() const
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

      // Compute element mass matrix
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

      // Invert mass matrix
      elem_mass_inv_[i].SetSize(ndof);
      DenseMatrixInverse M_inv(M);
      M_inv.GetInverseMatrix(elem_mass_inv_[i]);
   }

#ifdef MFEM_USE_MPI
   // In parallel, also compute mass inverses for face-neighbor elements.
   // These are needed by BR2InteriorFaceIntegrator::AssembleFaceMatrix when
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
const DenseMatrix& AntiplaneDomainOperator<MeshType>::GetMassMatrixInverse(
   int elem) const
{
   if (!mass_inv_computed_)
   {
      PrecomputeMassMatrixInverses();
   }
   return elem_mass_inv_[elem];
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::ComputeLiftingOperator(
   int elem, FaceElementTransformations &FTr, real_t jump, Vector &lift) const
{
   const FiniteElement *fe = fes_->GetFE(elem);
   int ndof = fe->GetDof();
   lift.SetSize(ndof);
   lift = 0.0;

   // Get the integration point mapping
   bool is_elem1 = (elem == FTr.Elem1No);

   // Face integration
   int face_order = fe->GetOrder();
   const IntegrationRule &ir = IntRules.Get(FTr.FaceGeom, 2 * face_order + 1);

   Vector shape(ndof);
   Vector int_result(ndof);
   int_result = 0.0;

   for (int p = 0; p < ir.GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir.IntPoint(p);
      FTr.SetAllIntPoints(&ip);

      // Get integration point in element reference coordinates
      const IntegrationPoint &eip = is_elem1 ?
                                    FTr.GetElement1IntPoint() : FTr.GetElement2IntPoint();

      // Evaluate shape functions at this point
      fe->CalcShape(eip, shape);

      // Face Jacobian gives the face length scaling
      real_t face_weight = ip.weight * FTr.Face->Weight();

      // Accumulate: ∫_e φ * jump * ds
      // For interior faces with two sides, use factor 0.5
      real_t factor = (FTr.Elem2No >= 0) ? 0.5 : 1.0;
      for (int k = 0; k < ndof; k++)
      {
         int_result(k) += factor * shape(k) * jump * face_weight;
      }
   }

   // Apply mass matrix inverse: lift = M^{-1} * int_result
   const DenseMatrix &Minv = GetMassMatrixInverse(elem);
   Minv.Mult(int_result, lift);
}

// ============================================================================
// Method-Specific Slip Assembly
// ============================================================================

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionIP(
   Vector &rhs, const Vector &slip_bc) const
{
   // IP method slip contribution for SIPG formulation.
   // RHS for imposed jump [[u]] = δ: σ⟨δ, {∇v·n}⟩_F + κ⟨δ, [[v]]⟩_F
   //
   // The symmetry (gradient) term uses element-specific weights:
   //   w_k = ip.weight / (2 * detJ_k)
   // The penalty (jump) term uses the combined weight from both elements:
   //   wq = kappa * |nor|² * (w1 + w2)
   // matching MFEM's DGDiffusionIntegrator which applies the same combined
   // penalty weight to all four blocks of the face matrix.
   real_t kappa = (order_ + 1) * (order_ + 1);

   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);

      if (FTr == nullptr) { continue; }

      // Get DOFs for elements on both sides of the face
      Array<int> dofs1, dofs2;
      fes_->GetElementDofs(FTr->Elem1No, dofs1);
      fes_->GetElementDofs(FTr->Elem2No, dofs2);

      // Get the finite elements
      const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
      const FiniteElement *fe2 = fes_->GetFE(FTr->Elem2No);

      // Use face integration rule
      int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
      const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);

      // Element-local vectors for accumulation
      Vector elvec1(dofs1.Size()), elvec2(dofs2.Size());
      elvec1 = 0.0;
      elvec2 = 0.0;

      // Each fault face has exactly one DOF (the midpoint).
      // Use this face's slip value directly — no depth interpolation.
      real_t slip_phys = slip_bc(i);
      if (std::abs(slip_phys) < 1e-15) { continue; }

      for (int p = 0; p < ir.GetNPoints(); p++)
      {
         const IntegrationPoint &ip = ir.IntPoint(p);
         FTr->SetAllIntPoints(&ip);

         const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

         // Face normal (points from elem1 to elem2)
         // CalcOrtho gives nor with |nor| = det(J_face) = face length in 2D
         Vector nor(mesh_.Dimension());
         CalcOrtho(FTr->Jacobian(), nor);

         // Determine imposed jump based on normal direction
         real_t slip_imposed = (nor(0) > 0) ? -slip_phys : slip_phys;

         // Shape functions
         Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
         fe1->CalcShape(eip1, shape1);
         fe2->CalcShape(eip2, shape2);

         // Gradient of shape functions (for symmetry term)
         DenseMatrix dshape1(fe1->GetDof(), mesh_.Dimension());
         DenseMatrix dshape2(fe2->GetDof(), mesh_.Dimension());
         fe1->CalcDShape(eip1, dshape1);
         fe2->CalcDShape(eip2, dshape2);

         // Use adjugate transformation for gradient
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

         // Weights following MFEM's DGDiffusionIntegrator (bilininteg.cpp:3829-3895)
         //
         // The symmetry (gradient) term uses element-specific weights:
         //   w_k = ip.weight / (2 * detJ_k) for each element k
         // This matches the bilinear form's gradient averaging.
         //
         // The penalty (jump) term must use the COMBINED weight from both
         // elements, because MFEM's DGDiffusionIntegrator accumulates
         //   wq = kappa * |nor|² * (w1 + w2)
         // and applies this same wq to ALL four blocks of jmat.
         // Using separate weights (w1 for elem1, w2 for elem2) gives only
         // half the penalty, causing the slip BC to be under-constrained.
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

      // Add element contributions to global RHS
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
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionBR2(
   Vector &rhs, const Vector &slip_bc) const
{
   // BR2 method slip contribution following Tandem's rhs_lift_skeleton and rhsFacet
   //
   // The RHS contribution for slip δ on interior face has two parts:
   // 1. Consistency: ε * 0.5 * ∫_F K∇φ·n * δ ds (same as IP)
   // 2. BR2 lifting: penalty * ∫_F φ * f_lifted_q ds
   //
   // where f_lifted_q is computed using BR2 lifting:
   //   f_lifted[i] = 0.5 * Minv[i] * ∫_F E[i] * δ * n * w ds
   //   f_lifted_q = 0.5 * n · (K * E[0] * f_lifted[0] + K * E[1] * f_lifted[1])
   //
   // For constant K=1, the penalty is the same mesh-dependent penalty used in
   // the bilinear form to ensure consistency.

   int dim = mesh_.Dimension();

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      int face = fault_interior_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);

      if (FTr == nullptr) { continue; }

      // Get DOFs and finite elements
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

      // Get mass matrix inverses
      const DenseMatrix &Minv1 = GetMassMatrixInverse(FTr->Elem1No);
      const DenseMatrix &Minv2 = GetMassMatrixInverse(FTr->Elem2No);

      // Each fault face has exactly one DOF (the midpoint).
      // Use this face's slip value directly — no depth interpolation.
      real_t slip_phys = slip_bc(fi);
      if (std::abs(slip_phys) < 1e-15) { continue; }

      // Precompute shapes, normals, and slip at quadrature points
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

         // Shape functions
         Vector shape1_q(shapes1.GetColumn(q), ndof1);
         Vector shape2_q(shapes2.GetColumn(q), ndof2);
         fe1->CalcShape(eip1, shape1_q);
         fe2->CalcShape(eip2, shape2_q);

         // Normal
         Vector nor_q(&nor_all[q * dim], dim);
         CalcOrtho(FTr->Jacobian(), nor_q);

         // Quadrature weight
         w_all[q] = ip.weight;

         // Sign based on normal direction (normal points from elem1 to elem2)
         slip_all[q] = (nor_q(0) > 0) ? -slip_phys : slip_phys;
      }

      // Check if there's any slip on this face
      real_t max_slip = 0.0;
      for (int q = 0; q < nqp; q++)
      {
         max_slip = std::max(max_slip, std::abs(slip_all[q]));
      }
      if (max_slip < 1e-15) { continue; }

      // ===================================================================
      // Compute BR2 lifted slip: f_lifted_q following Tandem's rhs_lift_skeleton
      // f_lifted[i][l,j] = 0.5 * Minv[i][l,m] * Σ_q E[i][m,q] * δ[q] * w[q] * n[j,q]
      // f_lifted_q[q] = 0.5 * n[j,q] * (E[0][l,q] * f_lifted[0][l,j] +
      //                                  E[1][l,q] * f_lifted[1][l,j])
      // ===================================================================

      // Compute f_lifted for each element
      // f_lifted[i] has shape (ndof_i, dim)
      DenseMatrix f_lifted1(ndof1, dim), f_lifted2(ndof2, dim);
      f_lifted1 = 0.0;
      f_lifted2 = 0.0;

      // First compute the face integral: ∫_F E * δ * n * w ds
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

      // Apply mass inverse: f_lifted = 0.5 * Minv * face_int
      // f_lifted[l,j] = 0.5 * Minv[l,m] * face_int[m,j]
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

      // Compute f_lifted_q at each quadrature point
      // f_lifted_q[q] = 0.5 * n[j,q] * (E[0][l,q] * f_lifted[0][l,j] +
      //                                  E[1][l,q] * f_lifted[1][l,j])
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

      // ===================================================================
      // Assemble RHS following Tandem's rhsFacet:
      // b[k] += c1 * ∫_F K∇φ_k·n * δ ds + c2 * ∫_F φ_k * f_lifted_q ds
      //
      // For element 1 (side 0):
      //   c1 = 0.5 * epsilon = -0.5 (SIPG)
      //   c2 = +penalty
      // For element 2 (side 1):
      //   c1 = 0.5 * epsilon = -0.5
      //   c2 = -penalty (opposite sign for element 2)
      // ===================================================================

      // Use the same BR2 penalty (sigma = D+1 = 3) as the bilinear form integrator
      // This ensures consistency between the bilinear form and RHS assembly.
      // The BR2 lifting already provides the correct h-scaling through M^{-1}.
      real_t penalty = br2_penalty_;  // D+1 = 3

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

         // Normal
         Vector nor_q(dim);
         for (int j = 0; j < dim; j++)
         {
            nor_q[j] = nor_all[q * dim + j];
         }

         // Gradient of shape functions (physical)
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

         // Consistency term: c1 * K∇φ·n * δ (c1 = sigma_ * 0.5 for averaging)
         real_t c1 = sigma_ * 0.5;  // = -0.5 for SIPG

         // Element 1 consistency
         for (int k = 0; k < ndof1; k++)
         {
            elvec1(k) += c1 * wq * dn1(k) * slip_q / detJ1;
         }

         // Element 2 consistency
         for (int k = 0; k < ndof2; k++)
         {
            elvec2(k) += c1 * wq * dn2(k) * slip_q / detJ2;
         }

         // BR2 lifting term: c2 * φ * f_lifted_q
         // Element 1: c2 = +penalty
         // Element 2: c2 = -penalty (opposite sign)
         for (int k = 0; k < ndof1; k++)
         {
            elvec1(k) += penalty * wq * shapes1(k, q) * f_lifted_q[q];
         }
         for (int k = 0; k < ndof2; k++)
         {
            elvec2(k) += (-penalty) * wq * shapes2(k, q) * f_lifted_q[q];
         }
      }

      // Add element contributions to global RHS
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
// Shared Face Slip Assembly (Parallel)
// ============================================================================

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionIPShared(
   Vector &rhs, const Vector &slip_bc, int interior_face_count) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // For shared faces, only Elem1 is local. We assemble one-sided
      // contributions — the other rank independently handles its Elem1.
      real_t kappa = (order_ + 1) * (order_ + 1);

      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);

         if (FTr == nullptr) { continue; }

         // Slip index: interior faces first, then shared faces
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
            // Same weights as serial interior face assembly
            real_t w1 = ip.weight / (2.0 * detJ1);
            real_t nor_sq = nor * nor;
            // Full penalty with both element weights (Elem2 info available via FTr)
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
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionBR2Shared(
   Vector &rhs, const Vector &slip_bc, int interior_face_count) const
{
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // BR2 shared face assembly: two-sided (Elem1 + Elem2)
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

         // Get Elem2 (face-neighbor) data
         int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
         int ndof2 = fe2->GetDof();

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);
         int nqp = ir.GetNPoints();

         const DenseMatrix &Minv1 = GetMassMatrixInverse(FTr->Elem1No);
         const DenseMatrix &Minv2 = GetMassMatrixInverse(FTr->Elem2No);

         // Precompute shapes, normals
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

         // Two-sided BR2 lifting (Elem1 + Elem2)
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

         // Apply mass inverse with factor 0.5 (same as interior)
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

         // Compute f_lifted_q (two-sided: Elem1 + Elem2 contributions)
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

         // Assemble RHS (Elem1 only)
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
using SerialAntiplaneOperator = AntiplaneDomainOperator<Mesh>;
#ifdef MFEM_USE_MPI
using ParallelAntiplaneOperator = AntiplaneDomainOperator<ParMesh>;
#endif

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ANTIPLANE_OPERATOR_HPP
