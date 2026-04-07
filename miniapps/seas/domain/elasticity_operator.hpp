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
#include "../integrator/dg_elasticity_ip_penalty_integrator.hpp"
#include "../integrator/dg_elasticity_ip_combined_integrator.hpp"
#include "../fault/face_quadrature.hpp"

#include <memory>
#include <cmath>
#include <limits>
#include <set>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace mfem
{
namespace seas
{

/// Linear solver type for the elasticity domain operator
enum class SolverType { CG_AMG, MUMPS, MUMPS_BLR, GMRES_BlockILU, GMRES_AMG, SUPERLU, STRUMPACK };

/// Boundary condition mode for the elasticity domain operator
///
/// Controls which boundary faces receive Dirichlet plate-rate loading
/// vs Natural (zero-traction) BC.
///
/// Tandem's BP5 (bp5.geo): top/bottom are Natural, far-field vertical
/// faces are Dirichlet. The physical surface numbers in Tandem's mesh
/// ARE the BC enum values (Natural=1, Fault=3, Dirichlet=5).
enum class BCMode
{
   /// Far-field: Dirichlet on attrs 1-4 (x=+-Lx, y=+-Ly),
   /// Natural on attrs 5-6 (z=0 free surface, z=Lz deep boundary)
   FarField,
   /// Antiplane-style: Dirichlet on attrs 1-2 (x=+-Lx) only,
   /// Natural on attrs 3-6
   XOnly,
   /// Legacy: Dirichlet on all attrs 1-6 (previous wrong implementation)
   AllDirichlet
};

/// @brief DG Elasticity domain operator for 3D vector elasticity (BP5)
///
/// Solves the 3D linear elasticity problem:
///   ∇·σ(u) = 0, σ = λ·tr(ε)·I + 2μ·ε
///
/// using DG method with fault slip as interior jump BC.
///
/// Coordinate convention (matches Tandem / SCEC BP5):
///   X = along-strike (SCEC x2)
///   Y = fault-normal (SCEC x1)
///   Z = depth, negative downward (Z=0 at surface)
///
/// - Fault at Y=0 is an interior interface
/// - Slip imposed as jump [[u]] on fault interior faces
/// - Boundary loading: u_X = sgn(Y) * Vp * t / 2 on Dirichlet faces
/// - Boundary conditions (Tandem Physical Surface tags):
///   - Tag 1 = Natural (top Z=0 + bottom Z=Z0)
///   - Tag 3 = Fault (Y=0 interior)
///   - Tag 5 = Dirichlet (far-field: Y=±Y1, X=±X1)
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
   /// @param solver_type Linear solver: CG_AMG, MUMPS, or GMRES_BlockILU
   /// @param bc_mode Boundary condition mode (FarField default)
   ElasticityDomainOperator(MeshType &mesh, int order,
                             real_t lambda, real_t mu,
                             real_t Vp, real_t Wf, real_t lf,
                             DGMethod method = DGMethod::BR2,
                             SolverType solver_type = SolverType::MUMPS_BLR,
                             BCMode bc_mode = BCMode::FarField,
                             int face_basis_type = BasisType::GaussLobatto)
      : mesh_(mesh), order_(order),
        lambda_val_(lambda), mu_val_(mu),
        Vp_(Vp), Wf_(Wf), lf_(lf),
        method_(method), solver_type_(solver_type),
        bc_mode_(bc_mode),
        check_residual_(false),
        lambda_coeff_(lambda), mu_coeff_(mu),
        mass_inv_computed_(false),
        fault_depths_computed_(false),
        fault_coords_computed_(false),
        face_basis_type_(face_basis_type)
   {
      MFEM_VERIFY(mesh_.Dimension() == 3, "ElasticityDomainOperator requires 3D mesh");

      epsilon_ = -1.0;  // SIPG

      SetupFESpace();
      SetupBoundaryMarkers();
      SetupFaultInfo();
      RunStartupFaceAudit();
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
                        Vector &traction,
                        Vector *normal_traction = nullptr) override;

   /// Compute total traction together with its stress and correction parts.
   /// The decomposition satisfies traction = traction_stress + traction_correction
   /// in the local fault basis (dip, strike) at each fault DOF.
   void ComputeTractionComponents(const GridFuncType &displacement,
                                 const Vector &slip_bc,
                                 Vector &traction,
                                 Vector &traction_stress,
                                 Vector &traction_correction,
                                 Vector *normal_traction = nullptr);

   /// Compute traction together with stress part, correction part, and
   /// fault jump residual in the local fault basis.
   /// jump_residual is stored as [dip, strike] interleaved and has units of m.
   void ComputeTractionDiagnostics(const GridFuncType &displacement,
                                   const Vector &slip_bc,
                                   Vector &traction,
                                   Vector &traction_stress,
                                   Vector &traction_correction,
                                   Vector &jump_residual,
                                   Vector *normal_traction = nullptr,
                                   Vector *normal_stress = nullptr,
                                   Vector *normal_correction = nullptr);

   /// Assemble only the fault-slip RHS contribution into the DG displacement
   /// space. This is a test/debug utility for checking K*u against b(slip)
   /// without Dirichlet loading.
   void AssembleSlipOnlyRHS(Vector &rhs, const Vector &slip_bc) const;

   /// Assemble only the Dirichlet loading RHS contribution into the DG
   /// displacement space. Debug utility for first-step RHS comparison.
   void AssembleDirichletOnlyRHS(Vector &rhs, real_t time) const;

   FESpaceType &GetFESpace() override { return *fes_; }
   const FESpaceType &GetFESpace() const override { return *fes_; }

   MeshType &GetMesh() override { return mesh_; }
   const MeshType &GetMesh() const override { return mesh_; }

   real_t GetShearModulus() const override { return mu_val_; }

   int GetNumFaultDOFs() const override { return num_fault_dofs_; }
   int GetNumOwnedFaultDOFs() const override { return num_owned_fault_dofs_; }

   void GetFaultDepths(Vector &depths) const override;

   void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const override;

   void RestrictToOwnedFault(const Vector &local_data,
                             Vector &owned_data,
                             int comps_per_dof = 1) const override;

   void ExpandOwnedToLocalFault(const Vector &owned_data,
                                Vector &local_data,
                                int comps_per_dof = 1) const override;

   const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

   const FaultBasis *GetFaultBasis() const override { return &fault_basis_; }

   const Array<int> &GetFaultInteriorFaces() const { return fault_interior_faces_; }
   const Array<int> &GetFaultSharedFaces() const { return fault_shared_faces_; }
   const Array<int> &GetOwnedFaultFaceMap() const
   { return owned_fault_face_to_local_face_; }
   int GetNumOwnedFaultFaces() const { return num_owned_fault_faces_; }

   DGMethod GetMethod() const { return method_; }
   int GetOrder() const { return order_; }
   real_t GetPlateRate() const { return Vp_; }
   real_t GetFaultDepthLimit() const { return Wf_; }
   real_t GetFaultLength() const { return lf_; }
   real_t GetLambda() const { return lambda_val_; }
   BCMode GetBCMode() const { return bc_mode_; }

   int GetNumFaultFaces() const override { return num_fault_faces_; }
   int GetNbfPerFace() const override { return nbf_per_face_; }
   const FaceQuadrature *GetFaceQuadrature() const { return face_quad_.get(); }

   /// Enable/disable post-solve residual check (||K*x - b|| / ||b||)
   void SetCheckResidual(bool check) { check_residual_ = check; }

   /// Set MUMPS-BLR tolerance (default 1e-10). Lower = more accurate, more memory.
   /// Only affects MUMPS_BLR solver type. Must be called BEFORE first Solve().
   void SetBLRTol(real_t tol) { blr_tol_ = tol; }

   /// v49: Use 2p quadrature instead of 2p+1 for fault face integration.
   /// Tests whether the extra quad point causes instability at p>=2.
   void SetMatchQuadOrder(bool v) { match_quad_order_ = v; }

   /// v50a: Scale IP penalty by this factor (1.0 = default, <1.0 = reduced).
   /// Used to diagnose whether over-stiff penalty causes nucleation failure at high p.
   void SetPenaltyFactor(real_t f) { penalty_factor_ = f; }

   /// v50g: Set face DOF node type for FaceQuadrature.
   /// Must be called BEFORE Init() (which creates FaceQuadrature).
   /// BasisType::GaussLobatto (default), BasisType::ClosedUniform, etc.
   void SetFaceBasisType(int bt) { face_basis_type_ = bt; }

   struct FirstStepDebugConfig
   {
      bool enabled = false;
      int target_rank = -1;
      std::string output_dir = ".";
      std::set<int> target_fault_faces;
   };

   void SetFirstStepDebugConfig(const FirstStepDebugConfig &cfg)
   {
      first_step_debug_ = cfg;
      debug_target_elems_cached_ = false;
   }

   bool IsFirstStepDebugEnabled() const { return first_step_debug_.enabled; }

private:
   MeshType &mesh_;
   int order_;
   real_t lambda_val_, mu_val_;
   real_t Vp_, Wf_, lf_;
   DGMethod method_;
   SolverType solver_type_;
   BCMode bc_mode_;
   bool check_residual_;  // Post-solve residual check
   real_t blr_tol_ = 1e-12;  // MUMPS-BLR factorization tolerance (v48: tightened from 1e-10)

   // Reference normal for skeleton Dirichlet orientation sign.
   // Matches Tandem's ref_normal from bp5.toml (default (0,-1,0) for BP5).
   // Used by ComputeSkeletonDirichletSign() to determine f_q sign on
   // interior/shared Dirichlet faces, same role as DGCurvilinearCommon.h:97.
   Vector ref_normal_;

   bool match_quad_order_ = false;       // Use 2p instead of 2p+1 quadrature
   real_t penalty_factor_ = 1.0;  // v50a: scale IP penalty (1.0=default)
   int face_basis_type_ = BasisType::GaussLobatto;  // v50g: face DOF node type

   enum class DebugAssemblePhase { None, Slip, Dirichlet };
   mutable FirstStepDebugConfig first_step_debug_;
   mutable bool first_step_debug_done_ = false;
   mutable DebugAssemblePhase debug_phase_ = DebugAssemblePhase::None;
   mutable real_t debug_time_ = 0.0;
   mutable bool debug_face_header_written_ = false;
   mutable bool debug_elem_header_written_ = false;
   mutable bool debug_jump_header_written_ = false;
   mutable bool debug_trac_header_written_ = false;
   mutable bool debug_target_elems_cached_ = false;
   mutable std::set<int> debug_target_elems_;

   void ComputeTractionImpl(const GridFuncType &displacement,
                            const Vector &slip_bc,
                            Vector &traction,
                            Vector *normal_traction,
                            Vector *traction_stress_out,
                            Vector *traction_correction_out,
                            Vector *jump_residual_out,
                            Vector *normal_stress_out = nullptr,
                            Vector *normal_correction_out = nullptr);

   int DebugRank() const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         return rank;
#endif
      }
      return 0;
   }

   bool DebugEnabledForTime(real_t time) const
   {
      if (!first_step_debug_.enabled || first_step_debug_done_) { return false; }
      // Target: dump at first accepted step (last RK45 stage).
      // With dt_init=0.02, the last stage is at t=0.02.
      // Threshold skips intermediate RK45 stages (t < 0.019).
      if (time < 0.019) { return false; }
      return (first_step_debug_.target_rank < 0 ||
              DebugRank() == first_step_debug_.target_rank);
   }

   const char *DebugPhaseName() const
   {
      switch (debug_phase_)
      {
         case DebugAssemblePhase::Slip: return "slip";
         case DebugAssemblePhase::Dirichlet: return "dirichlet";
         default: return "none";
      }
   }

   std::string DebugFilePath(const std::string &stem) const
   {
      std::ostringstream oss;
      oss << first_step_debug_.output_dir;
      if (!first_step_debug_.output_dir.empty() &&
          first_step_debug_.output_dir.back() != '/')
      {
         oss << "/";
      }
      oss << stem << "_r" << DebugRank() << ".csv";
      return oss.str();
   }

   const std::set<int> &DebugTargetLocalElements() const
   {
      if (debug_target_elems_cached_) { return debug_target_elems_; }
      debug_target_elems_.clear();
      const int num_local = mesh_.GetNE();
      const int interior_count = fault_interior_faces_.Size();
      for (int fault_idx : first_step_debug_.target_fault_faces)
      {
         if (fault_idx < 0) { continue; }
         FaceElementTransformations *FTr = nullptr;
         if (fault_idx < interior_count)
         {
            int face = fault_interior_faces_[fault_idx];
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            int sh_idx = fault_idx - interior_count;
            if (sh_idx >= 0 && sh_idx < fault_shared_faces_.Size())
            {
               FTr = mesh_.GetSharedFaceTransformations(fault_shared_faces_[sh_idx]);
            }
         }
         if (!FTr) { continue; }
         if (FTr->Elem1No >= 0 && FTr->Elem1No < num_local) { debug_target_elems_.insert(FTr->Elem1No); }
         if (FTr->Elem2No >= 0 && FTr->Elem2No < num_local) { debug_target_elems_.insert(FTr->Elem2No); }
      }
      debug_target_elems_cached_ = true;
      return debug_target_elems_;
   }

   bool DebugShouldDumpFace(int logical_fault_idx,
                            FaceElementTransformations *FTr) const
   {
      if (!DebugEnabledForTime(debug_time_) || FTr == nullptr) { return false; }
      if (logical_fault_idx >= 0 &&
          first_step_debug_.target_fault_faces.count(logical_fault_idx) > 0)
      {
         return true;
      }

      const auto &target_elems = DebugTargetLocalElements();
      return target_elems.count(FTr->Elem1No) > 0 ||
             target_elems.count(FTr->Elem2No) > 0;
   }

   void DebugDumpFaceData(int logical_fault_idx,
                          int mesh_face_idx,
                          const char *face_kind,
                          FaceElementTransformations *FTr,
                          const Vector *phys_y_qp,
                          const Vector *input_qp,
                          const Vector &elvec1,
                          const Vector *elvec2) const
   {
      if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

      // Truncate on first write to discard stale data from previous runs
      std::ofstream out(DebugFilePath("first_step_face_rhs"),
                        debug_face_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_face_header_written_)
      {
         out << "time,phase,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_face_header_written_ = true;
      }

      const int elem1 = FTr ? FTr->Elem1No : -1;
      const int elem2 = FTr ? FTr->Elem2No : -1;
      auto write_row = [&](const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << DebugPhaseName() << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << i0 << ","
             << i1 << ","
             << value << "\n";
      };

      if (phys_y_qp)
      {
         for (int q = 0; q < phys_y_qp->Size(); q++)
         {
            write_row("phys_y", q, -1, (*phys_y_qp)(q));
         }
      }
      if (input_qp)
      {
         const int nq = input_qp->Size() / 3;
         for (int q = 0; q < nq; q++)
         {
            for (int c = 0; c < 3; c++)
            {
               write_row("input_qp", q, c, (*input_qp)(c * nq + q));
            }
         }
      }
      for (int j = 0; j < elvec1.Size(); j++)
      {
         write_row("elvec1", j, -1, elvec1(j));
      }
      if (elvec2)
      {
         for (int j = 0; j < elvec2->Size(); j++)
         {
            write_row("elvec2", j, -1, (*elvec2)(j));
         }
      }
   }

   void DebugDumpElementVector(const char *quantity, const Vector &vec) const
   {
      if (!DebugEnabledForTime(debug_time_)) { return; }

      std::ofstream out(DebugFilePath("first_step_elem_data"),
                        debug_elem_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_elem_header_written_)
      {
         out << "time,quantity,elem,component,local_dof,vdof,value\n";
         debug_elem_header_written_ = true;
      }

      const auto &target_elems = DebugTargetLocalElements();
      for (int elem : target_elems)
      {
         Array<int> vdofs;
         fes_->GetElementVDofs(elem, vdofs);
         const int ndof = scalar_fes_->GetFE(elem)->GetDof();
         for (int j = 0; j < vdofs.Size(); j++)
         {
            const int vdof = vdofs[j];
            const int lid = (vdof >= 0) ? vdof : (-1 - vdof);
            const real_t value = (vdof >= 0) ? vec(lid) : -vec(lid);
            out << std::setprecision(17) << debug_time_ << ","
                << quantity << ","
                << elem << ","
                << (j / ndof) << ","
                << (j % ndof) << ","
                << vdof << ","
                << value << "\n";
         }
      }
   }

   /// Dump per-element and per-face K matrix contributions for target elements.
   /// Each contribution is dumped separately so we can compare against Tandem's
   /// assemble_volume / assemble_skeleton / assemble_boundary piece by piece.
   /// MPI-safe: only local elements and local interior faces.
   void DebugDumpKContributions() const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      const int dim = 3;
      const auto &target_elems = DebugTargetLocalElements();
      if (target_elems.empty()) { return; }

      auto write_mat = [](std::ofstream &os, const char *source, int id,
                          int e1, int e2, const char *block,
                          const DenseMatrix &M)
      {
         for (int i = 0; i < M.Height(); i++)
         {
            for (int j = 0; j < M.Width(); j++)
            {
               real_t v = M(i, j);
               if (std::abs(v) > 1e-30)
               {
                  os << source << "," << id << "," << e1 << "," << e2
                     << "," << block << "," << i << "," << j << ","
                     << std::setprecision(17) << v << "\n";
               }
            }
         }
      };

      // 1. Volume contributions: ElasticityIntegrator on target elements
      {
         std::ofstream out(DebugFilePath("first_step_K_volume"), std::ios::trunc);
         if (out)
         {
            out << "source,elem,elem1,elem2,block,row,col,value\n";
            ElasticityIntegrator vol_integ(lambda_coeff_, mu_coeff_);
            for (int elem : target_elems)
            {
               const FiniteElement *fe = scalar_fes_->GetFE(elem);
               ElementTransformation *eltrans =
                  mesh_.GetElementTransformation(elem);
               DenseMatrix K_vol;
               vol_integ.AssembleElementMatrix(*fe, *eltrans, K_vol);
               write_mat(out, "volume", elem, elem, -1, "A00", K_vol);
            }
         }
      }

      // 2. Interior face contributions: DGElasticityIPCombinedIntegrator
      {
         std::ofstream out(DebugFilePath("first_step_K_skeleton"),
                           std::ios::trunc);
         if (out)
         {
            out << "source,face,elem1,elem2,block,row,col,value\n";
            DGElasticityIPCombinedIntegrator face_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int nfaces = mesh_.GetNumFaces();
            for (int f = 0; f < nfaces; f++)
            {
               FaceElementTransformations *FTr =
                  mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }

               bool e1_target = target_elems.count(FTr->Elem1No) > 0;
               bool e2_target = (FTr->Elem2No >= 0 &&
                                 target_elems.count(FTr->Elem2No) > 0);
               if (!e1_target && !e2_target) { continue; }

               const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
               const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
               DenseMatrix K_face;
               face_integ.AssembleFaceMatrix(*fe1, *fe2, *FTr, K_face);

               int n1 = fe1->GetDof() * dim;
               int n2 = fe2->GetDof() * dim;
               DenseMatrix A00(n1, n1), A01(n1, n2), A10(n2, n1), A11(n2, n2);
               K_face.GetSubMatrix(0, n1, 0, n1, A00);
               K_face.GetSubMatrix(0, n1, n1, n1 + n2, A01);
               K_face.GetSubMatrix(n1, n1 + n2, 0, n1, A10);
               K_face.GetSubMatrix(n1, n1 + n2, n1, n1 + n2, A11);

               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A00", A00);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A01", A01);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A10", A10);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A11", A11);
            }

            // Boundary face contributions (if any target element touches one)
            for (int be = 0; be < mesh_.GetNBE(); be++)
            {
               int attr = mesh_.GetBdrAttribute(be);
               if (dirichlet_bdr_marker_.Size() > 0 &&
                   dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

               int face_idx, face_info_val;
               mesh_.GetBdrElementFace(be, &face_idx, &face_info_val);
               FaceElementTransformations *FTr =
                  mesh_.GetFaceElementTransformations(face_idx);
               if (!FTr || target_elems.count(FTr->Elem1No) == 0)
               { continue; }

               const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
               DenseMatrix K_bdr;
               face_integ.AssembleFaceMatrix(*fe1, *fe1, *FTr, K_bdr);
               write_mat(out, "boundary_face", face_idx, FTr->Elem1No, -1,
                         "A00", K_bdr);
            }
         }
      }
   }

   /// Local-only variant: dumps interior fault faces only (no MPI exchange).
   /// Safe to call on a single rank without deadlocking.
   void DebugDumpFaultJumpsLocal(const GridFuncType &displacement,
                                 const Vector &slip_bc) const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_jump"),
                        debug_jump_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_jump_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,q,component,value\n";
         debug_jump_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int q, int c, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << "," << logical_fault_idx << ","
             << mesh_face_idx << "," << elem1 << "," << elem2 << ","
             << record << "," << q << "," << c << "," << value << "\n";
      };

      const int dim = 3;

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         const int ndof1 = fe1->GetDof();
         const int ndof2 = fe2->GetDof();
         const int nq = face_quad_->NumQuadPoints();

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         const IntegrationRule &ir = face_quad_->GetQuadRule();
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &fip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&fip);

            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));

            Vector s1(ndof1), s2(ndof2);
            fe1->CalcShape(FTr->GetElement1IntPoint(), s1);
            fe2->CalcShape(FTr->GetElement2IntPoint(), s2);
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++) { u1q += s1(k) * u1_all(c * ndof1 + k); }
               for (int k = 0; k < ndof2; k++) { u2q += s2(k) * u2_all(c * ndof2 + k); }
               const real_t slipq = delta_u_quad(c * nq + q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "u1_q", q, c, u1q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "u2_q", q, c, u2q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "slip_q", q, c, slipq);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "jump_minus_slip", q, c,
                         (u1q - u2q) - slipq);
            }
         }
      }
   }

   /// Local-only variant: dumps interior fault face traction (no MPI exchange).
   /// Local-only variant: computes traction per-face inline for interior
   /// faces only. No ComputeTraction call (which uses MPI for shared faces).
   void DebugDumpFaultTractionLocal(const GridFuncType &displacement,
                                     const Vector &slip_bc)
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_trac"),
                        debug_trac_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_trac_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_trac_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << "," << logical_fault_idx << ","
             << mesh_face_idx << "," << elem1 << "," << elem2 << ","
             << record << "," << i0 << "," << i1 << "," << value << "\n";
      };

      const int dim = 3;
      const int nbf = nbf_per_face_;

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad, nl_q;
         trac_integ.ComputeTractionAtQuadPoints(
            *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad, T_quad, nullptr, &nl_q);

         const int nq = T_quad.Size() / dim;
         const IntegrationRule &ir = IntRules.Get(
            FTr->GetGeometryType(), 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1);
         for (int q = 0; q < nq; q++)
         {
            FTr->SetAllIntPoints(&ir.IntPoint(q));
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "nl_q", q, -1, nl_q(q));
            for (int c = 0; c < dim; c++)
            {
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "traction_q", q, c,
                         T_quad(c * nq + q));
            }
         }

         // Note: projected DOF-level traction (traction_dip/strike/normal)
         // is omitted here because it requires ComputeTraction which uses MPI.
         // The per-QP traction_q values above are sufficient for cross-code
         // comparison — project offline if needed.
      }
   }

   void DebugDumpFaultJumps(const GridFuncType &displacement,
                            const Vector &slip_bc) const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_jump"),
                        debug_jump_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_jump_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,q,component,value\n";
         debug_jump_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int q, int c, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << q << ","
             << c << ","
             << value << "\n";
      };

      auto dump_face = [&](int logical_fault_idx, int mesh_face_idx,
                           const char *face_kind, FaceElementTransformations *FTr,
                           const FiniteElement *fe1, const Vector &u1_all,
                           const FiniteElement *fe2, const Vector &u2_all)
      {
         if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

         const int dim = 3;
         const int ndof1 = fe1->GetDof();
         const int ndof2 = fe2->GetDof();
         const int nq = face_quad_->NumQuadPoints();

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(logical_fault_idx, slip_bc, delta_u_quad);

         const IntegrationRule &ir = face_quad_->GetQuadRule();
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &fip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&fip);

            // Physical coordinates at this QP (for cross-code matching)
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));

            Vector s1(ndof1), s2(ndof2);
            fe1->CalcShape(FTr->GetElement1IntPoint(), s1);
            fe2->CalcShape(FTr->GetElement2IntPoint(), s2);
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++) { u1q += s1(k) * u1_all(c * ndof1 + k); }
               for (int k = 0; k < ndof2; k++) { u2q += s2(k) * u2_all(c * ndof2 + k); }
               const real_t slipq = delta_u_quad(c * nq + q);
               const real_t jump_minus_slip = (u1q - u2q) - slipq;
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "u1_q", q, c, u1q);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "u2_q", q, c, u2q);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "slip_q", q, c, slipq);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No,
                         "jump_minus_slip", q, c, jump_minus_slip);
            }
         }
      };

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);
         dump_face(fi, face, "fault_interior", FTr,
                   scalar_fes_->GetFE(FTr->Elem1No), u1_all,
                   scalar_fes_->GetFE(FTr->Elem2No), u2_all);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (!pfes) { return; }

         ParGridFunction par_u(pfes);
         par_u = displacement;
         pfes->ExchangeFaceNbrData();
         par_u.ExchangeFaceNbrData();

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int logical_fault_idx = fault_interior_faces_.Size() + i;
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr = mesh_.GetSharedFaceTransformations(sf);
            if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { continue; }

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            Vector u1_all(vdofs1.Size());
            par_u.GetSubVector(vdofs1, u1_all);

            const int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            Array<int> vdofs2;
            pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
            const Vector &nbr_data = par_u.FaceNbrData();
            Vector u2_all(vdofs2.Size());
            for (int j = 0; j < vdofs2.Size(); j++)
            {
               u2_all(j) = nbr_data(vdofs2[j]);
            }

            dump_face(logical_fault_idx, mesh_.GetSharedFace(sf), "fault_shared", FTr,
                      scalar_fes_->GetFE(FTr->Elem1No), u1_all, fe2, u2_all);
         }
#endif
      }
   }

   void DebugDumpFaultTraction(const GridFuncType &displacement,
                               const Vector &slip_bc)
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      Vector traction, normal_traction;
      ComputeTraction(displacement, slip_bc, traction, &normal_traction);

      std::ofstream out(DebugFilePath("first_step_face_trac"),
                        debug_trac_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_trac_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_trac_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << i0 << ","
             << i1 << ","
             << value << "\n";
      };

      auto dump_face = [&](int logical_fault_idx, int mesh_face_idx,
                           const char *face_kind, FaceElementTransformations *FTr,
                           const FiniteElement *fe1, const Vector &u1_all,
                           const FiniteElement *fe2, const Vector &u2_all)
      {
         if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

         const int dim = 3;
         const int nbf = nbf_per_face_;
         const int nq = face_quad_->NumQuadPoints();

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(logical_fault_idx, slip_bc, delta_u_quad);

         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad, nl_q;
         trac_integ.ComputeTractionAtQuadPoints(
            *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad, T_quad, nullptr, &nl_q);

         // Compute physical coordinates at each QP for cross-code matching
         const IntegrationRule &ir_trac = IntRules.Get(
            FTr->GetGeometryType(), 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1);
         for (int q = 0; q < nq; q++)
         {
            FTr->SetAllIntPoints(&ir_trac.IntPoint(q));
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "nl_q", q, -1, nl_q(q));
            for (int c = 0; c < dim; c++)
            {
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "traction_q", q, c,
                         T_quad(c * nq + q));
            }
         }

         for (int kk = 0; kk < nbf; kk++)
         {
            const int dof_idx = logical_fault_idx * nbf + kk;
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "traction_dip", kk, -1,
                      traction(2 * dof_idx));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "traction_strike", kk, -1,
                      traction(2 * dof_idx + 1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "normal_traction", kk, -1,
                      normal_traction(dof_idx));
         }
      };

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);
         dump_face(fi, face, "fault_interior", FTr,
                   scalar_fes_->GetFE(FTr->Elem1No), u1_all,
                   scalar_fes_->GetFE(FTr->Elem2No), u2_all);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (!pfes) { return; }

         ParGridFunction par_u(pfes);
         par_u = displacement;
         pfes->ExchangeFaceNbrData();
         par_u.ExchangeFaceNbrData();

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int logical_fault_idx = fault_interior_faces_.Size() + i;
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr = mesh_.GetSharedFaceTransformations(sf);
            if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { continue; }

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            Vector u1_all(vdofs1.Size());
            par_u.GetSubVector(vdofs1, u1_all);

            const int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            Array<int> vdofs2;
            pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
            const Vector &nbr_data = par_u.FaceNbrData();
            Vector u2_all(vdofs2.Size());
            for (int j = 0; j < vdofs2.Size(); j++)
            {
               u2_all(j) = nbr_data(vdofs2[j]);
            }

            dump_face(logical_fault_idx, mesh_.GetSharedFace(sf), "fault_shared", FTr,
                      scalar_fes_->GetFE(FTr->Elem1No), u1_all, fe2, u2_all);
         }
#endif
      }
   }

   /// Build 3D slip at quadrature points for a single fault face.
   ///
   /// Shared by production assembly (AssembleSlipContributionIP) and debug
   /// diagnostics (DebugDumpFaultJumps, DebugDumpFaultTraction).  Single
   /// source of truth for the interpolation + embedding + sign chain.
   ///
   /// @param logical_fault_idx  Face index in the combined interior+shared list
   /// @param slip_bc            Full local slip vector [2 * num_fault_dofs]
   /// @param[out] delta_u_quad  3D slip at QPs [dim * nq], sign-corrected
   void BuildSlipAtQuadPoints(int logical_fault_idx, const Vector &slip_bc,
                              Vector &delta_u_quad) const
   {
      const int dim = 3;
      const int nbf = nbf_per_face_;
      const auto &basis = fault_basis_.GetBasis(logical_fault_idx);
      // Tandem convention: sign is baked into the basis vectors.
      // No separate sign factor needed.

      if (!basis.qp_data.empty())
      {
         Vector slip_tang(2 * nbf);
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = logical_fault_idx * nbf + kk;
            slip_tang(0 * nbf + kk) = slip_bc(2 * dof_idx);
            slip_tang(1 * nbf + kk) = slip_bc(2 * dof_idx + 1);
         }
         Vector slip_tang_q;
         face_quad_->InterpolateToQuadPoints(2, slip_tang, slip_tang_q);
         int nqp = slip_tang_q.Size() / 2;
         delta_u_quad.SetSize(dim * nqp);
         for (int q = 0; q < nqp; q++)
         {
            real_t sl_q[2] = {slip_tang_q(q), slip_tang_q(nqp + q)};
            real_t du[3];
            fault_basis_.EmbedSlipQP(logical_fault_idx, q, sl_q, du);
            for (int c = 0; c < dim; c++)
               delta_u_quad(c * nqp + q) = du[c];
         }
      }
      else
      {
         Vector delta_u_nodal(dim * nbf);
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = logical_fault_idx * nbf + kk;
            real_t slip_local[2] = {slip_bc(2 * dof_idx),
                                    slip_bc(2 * dof_idx + 1)};
            real_t du[3];
            fault_basis_.EmbedSlip(logical_fault_idx, slip_local, du);
            for (int c = 0; c < dim; c++)
               delta_u_nodal(c * nbf + kk) = du[c];
         }
         face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal, delta_u_quad);
      }
   }

   // ---- Facet BC classification (single source of truth) ----
   // Mirrors Tandem's per-facet BC enum (DGOperatorTopo FacetInfo.bc).
   // Built once at startup from mesh tags, propagated exactly via MPI.
   enum class FacetBC : int8_t { None = 0, Fault = 1, Dirichlet = 2 };

   std::vector<FacetBC> face_bc_;         // [mesh_.GetNumFaces()] interior faces
   std::vector<FacetBC> shared_face_bc_;  // [mesh_.GetNSharedFaces()] shared faces

   // Legacy face arrays derived from the BC tables (kept for assembly loops)
   Array<int> fault_interior_faces_;
   Array<int> fault_shared_faces_;
   Array<int> dirichlet_interior_faces_;
   Array<int> dirichlet_shared_faces_;

   // Element-pair keys for fault faces (populated during BC table build)
   std::set<long> fault_face_keys_;

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
#ifdef MFEM_USE_SUPERLU
   mutable std::unique_ptr<SuperLURowLocMatrix> cached_superlu_mat_;
#endif
#ifdef MFEM_USE_STRUMPACK
   mutable std::unique_ptr<STRUMPACKRowLocMatrix> cached_strumpack_mat_;
#endif

   // Fault data
   Array<int> fault_dofs_;
   int num_fault_dofs_ = 0;
   int num_owned_fault_dofs_ = 0;
   FaultBasis fault_basis_;

   // Multi-DOF fault discretization (v45, Phase 2)
   std::unique_ptr<FaceQuadrature> face_quad_;
   int nbf_per_face_ = 1;   // BR2: 1, IP: (p+1)(p+2)/2 on triangle faces
   int num_fault_faces_ = 0; // number of fault faces (interior + shared)
   int num_owned_fault_faces_ = 0;

   Array<int> owned_fault_face_to_local_face_;
   Array<int> full_fault_face_to_owned_face_;
   Array<int> owned_fault_dof_to_local_dof_;

   /// Canonical face key: sorted global vertex IDs for triangular faces.
   /// BP5/tet-specific (3 vertices per face). Would need extension for
   /// quad faces (hex meshes) if used in that context.
   struct FaceVertexKey
   {
      HYPRE_BigInt v[3] = {-1, -1, -1};
      bool operator<(const FaceVertexKey &o) const
      {
         if (v[0] != o.v[0]) return v[0] < o.v[0];
         if (v[1] != o.v[1]) return v[1] < o.v[1];
         return v[2] < o.v[2];
      }
      bool operator==(const FaceVertexKey &o) const
      {
         return v[0] == o.v[0] && v[1] == o.v[1] && v[2] == o.v[2];
      }
   };

   /// Build a canonical face key from a local face index (sorted global vertex IDs).
   FaceVertexKey MakeFaceKey(int local_face,
                              const Array<HYPRE_BigInt> &gvert) const
   {
      Array<int> verts;
      mesh_.GetFaceVertices(local_face, verts);
      FaceVertexKey key;
      for (int j = 0; j < 3 && j < verts.Size(); j++)
      {
         key.v[j] = gvert[verts[j]];
      }
      // Sort to canonical order
      if (key.v[0] > key.v[1]) { std::swap(key.v[0], key.v[1]); }
      if (key.v[1] > key.v[2]) { std::swap(key.v[1], key.v[2]); }
      if (key.v[0] > key.v[1]) { std::swap(key.v[0], key.v[1]); }
      return key;
   }

   struct SharedFaultFaceBlock
   {
      FaceVertexKey key;
      int face_idx = -1;
   };

   struct SharedFaultCommBlock
   {
      int neighbor_rank = -1;
      std::vector<int> send_owned_faces;
      std::vector<int> recv_local_faces;
   };

   std::vector<SharedFaultCommBlock> shared_fault_comm_blocks_;

   /// Per-face canonical DOF permutation (Tandem sorted-simplex convention).
   /// Flat array: canonical_to_local_perm_[face_idx * nbf_per_face_ + canonical_k]
   ///           = mfem_local_k
   /// where canonical order = sorted by ascending global vertex ID.
   /// For p=1 triangles: DOF k = vertex k, so vertex perm = DOF perm.
   Array<int> canonical_to_local_perm_;

   mutable Vector fault_depths_;
   mutable bool fault_depths_computed_;
   mutable Vector fault_x2_, fault_x3_;
   mutable bool fault_coords_computed_;

   // Boundary markers
   mutable Array<int> dirichlet_bdr_marker_;

   // ========================================================================
   // Helpers
   // ========================================================================

   /// Compute the skeleton Dirichlet orientation sign at the current face integration point.
   ///
   /// Matches Tandem DGCurvilinearCommon.h:95-99: on skeleton (interior)
   /// faces, negate f_q when the face normal opposes ref_normal_.
   /// This is the same role as sign_flipped for fault slip, generalized
   /// to use the full dot product with the stored ref_normal_ vector.
   ///
   /// @param FTr Face transformation (SetAllIntPoints must have been called)
   /// @return +1.0 if face normal aligns with ref_normal_, -1.0 otherwise
   real_t ComputeSkeletonDirichletSign(
      FaceElementTransformations *FTr) const
   {
      Vector nor(3);
      CalcOrtho(FTr->Jacobian(), nor);
      real_t dot_ref = 0.0;
      for (int d = 0; d < ref_normal_.Size(); d++)
      {
         dot_ref += nor(d) * ref_normal_(d);
      }
      return (dot_ref < 0.0) ? -1.0 : 1.0;
   }

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
      // Tandem Physical Surface tags:
      //   Tag 1 = Natural (top Z=0 + bottom Z=Z0) → zero traction
      //   Tag 3 = Fault (Y=0 interior) → handled separately
      //   Tag 5 = Dirichlet (far-field: Y=±Y1, X=±X1) → plate loading

      int num_bdr = mesh_.bdr_attributes.Size() > 0 ? mesh_.bdr_attributes.Max() : 0;
      dirichlet_bdr_marker_.SetSize(num_bdr);
      dirichlet_bdr_marker_ = 0;

      if (bc_mode_ == BCMode::AllDirichlet)
      {
         mfem::out << "\n  *** WARNING: AllDirichlet BC mode is legacy and known "
                   << "to be incorrect for BP5. ***\n"
                   << "  *** Use BCMode::FarField (default) for production runs. "
                   << "***\n\n";
         for (int i = 0; i < num_bdr; i++)
         {
            dirichlet_bdr_marker_[i] = 1;
         }
      }
      else if (bc_mode_ == BCMode::FarField || bc_mode_ == BCMode::XOnly)
      {
         // Mark only attr 5 as Dirichlet (Tandem tag for far-field faces)
         for (int i = 0; i < num_bdr; i++)
         {
            int attr = i + 1;
            if (attr == 5) { dirichlet_bdr_marker_[i] = 1; }
         }
      }
   }

   /// Build the facet BC tables from mesh tags. Single source of truth.
   ///
   /// Mirrors Tandem's facet-BC model: import tags → propagate exactly →
   /// store per-facet class → consume per-facet class.
   ///
   /// 1. Scan boundary elements with attr 3 (fault) and attr 5 (Dirichlet)
   /// 2. Mark interior faces directly in face_bc_
   /// 3. Collect shared-face canonical vertex keys, Allgatherv across ranks
   /// 4. Mark shared faces in shared_face_bc_
   /// 5. Validate: attr 3 and attr 5 must exist; fault ∩ Dirichlet = ∅
   /// 6. Derive legacy face arrays for downstream assembly
   void BuildFacetBCTables()
   {
      const int num_faces = mesh_.GetNumFaces();
      face_bc_.assign(num_faces, FacetBC::None);

      int num_shared = 0;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         num_shared = mesh_.GetNSharedFaces();
#endif
      }
      shared_face_bc_.assign(num_shared, FacetBC::None);

      // ---- Check global attr existence ----
      int local_has3 = 0, local_has5 = 0;
      // Scan actual boundary elements instead of mesh_.bdr_attributes. In
      // ParMesh, internal tagged faces can appear in GetNBE() even when the
      // summary attribute list does not include their tag.
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr == 3) { local_has3 = 1; }
         if (attr == 5) { local_has5 = 1; }
      }
      int global_has3 = local_has3, global_has5 = local_has5;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         MPI_Allreduce(MPI_IN_PLACE, &global_has3, 1, MPI_INT,
                       MPI_MAX, mesh_.GetComm());
         MPI_Allreduce(MPI_IN_PLACE, &global_has5, 1, MPI_INT,
                       MPI_MAX, mesh_.GetComm());
#endif
      }
      // If the mesh has no fault tags (attr 3), this is a non-BP5 problem
      // (unit test, pure elasticity, etc.). Skip fault classification but
      // still build Dirichlet faces if attr 5 is present.
      // If the mesh HAS fault tags, enforce BP5 requirement: attr 5 must
      // also be present.
      bool has_fault = (global_has3 > 0);
      bool has_dirichlet = (global_has5 > 0);
      if (has_fault)
      {
         MFEM_VERIFY(has_dirichlet,
            "ERROR: Mesh has fault faces (attr 3) but no Dirichlet faces "
            "(attr 5). BP5 requires both Physical Surface tags.");
      }
      if (!has_fault && !has_dirichlet) { return; }

      // ---- Build reverse map: local face → shared face index ----
      std::unordered_map<int, int> lface_to_sface;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int sf = 0; sf < num_shared; sf++)
         {
            lface_to_sface[mesh_.GetSharedFace(sf)] = sf;
         }
#endif
      }

      // ---- Phase 1: Scan boundary elements, mark local faces ----
      // Collect shared-face tags locally (only 1 rank has the bdr element)
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr != 3 && attr != 5) { continue; }
         FacetBC bc = (attr == 3) ? FacetBC::Fault : FacetBC::Dirichlet;
         int face_idx = mesh_.GetBdrElementFaceIndex(be);

         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face_idx);
         if (FTr != nullptr)
         {
            // Interior face: mark in face_bc_ table
            MFEM_VERIFY(face_bc_[face_idx] == FacetBC::None ||
                        face_bc_[face_idx] == bc,
               "ERROR: Interior face " << face_idx << " tagged as both "
               "attr 3 and attr 5.");
            face_bc_[face_idx] = bc;

            if (bc == FacetBC::Fault)
            {
               int e1 = FTr->Elem1No;
               int e2 = FTr->Elem2No;
               long key = (long)std::min(e1, e2) * mesh_.GetNE()
                          + std::max(e1, e2);
               fault_face_keys_.insert(key);
            }
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            auto it = lface_to_sface.find(face_idx);
            if (it != lface_to_sface.end())
            {
               int sf = it->second;
               MFEM_VERIFY(shared_face_bc_[sf] == FacetBC::None ||
                           shared_face_bc_[sf] == bc,
                  "ERROR: Shared face " << sf << " tagged as both "
                  "attr 3 and attr 5.");
               shared_face_bc_[sf] = bc;
            }
         }
      }

      // ---- Phase 2: Propagate shared-face tags across MPI ranks ----
      // A shared face's boundary element exists on only ONE rank. The
      // other rank must discover it via canonical vertex key exchange.
      // We propagate fault and Dirichlet tags in a single Allgatherv,
      // packing (key, bc_class) per face.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         // Pack: 4 values per locally-tagged shared face (3 key + bc_class)
         std::vector<HYPRE_BigInt> local_flat;
         for (int sf = 0; sf < num_shared; sf++)
         {
            if (shared_face_bc_[sf] == FacetBC::None) { continue; }
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey fk = MakeFaceKey(lf, gvert);
            local_flat.push_back(fk.v[0]);
            local_flat.push_back(fk.v[1]);
            local_flat.push_back(fk.v[2]);
            local_flat.push_back(static_cast<HYPRE_BigInt>(shared_face_bc_[sf]));
         }

         int lc = static_cast<int>(local_flat.size());
         int nranks = 1;
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         std::vector<int> rc(nranks), dp(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, rc.data(), 1, MPI_INT,
                       mesh_.GetComm());
         int tot = 0;
         for (int r = 0; r < nranks; r++) { dp[r] = tot; tot += rc[r]; }
         std::vector<HYPRE_BigInt> all(tot);
         MPI_Allgatherv(local_flat.data(), lc, HYPRE_MPI_BIG_INT,
                        all.data(), rc.data(), dp.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Build global key → bc map
         std::map<FaceVertexKey, FacetBC> global_shared_bc;
         for (int i = 0; i < tot; i += 4)
         {
            FaceVertexKey k;
            k.v[0] = all[i]; k.v[1] = all[i+1]; k.v[2] = all[i+2];
            FacetBC bc = static_cast<FacetBC>(all[i+3]);
            auto [it, inserted] = global_shared_bc.emplace(k, bc);
            MFEM_VERIFY(inserted || it->second == bc,
               "ERROR: Shared face key (" << k.v[0] << "," << k.v[1]
               << "," << k.v[2] << ") has conflicting tags from "
               "different ranks.");
         }

         // Fill untagged local shared faces from the global map
         for (int sf = 0; sf < num_shared; sf++)
         {
            if (shared_face_bc_[sf] != FacetBC::None) { continue; }
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey key = MakeFaceKey(lf, gvert);
            auto it = global_shared_bc.find(key);
            if (it != global_shared_bc.end())
            {
               shared_face_bc_[sf] = it->second;
            }
         }
#endif
      }

      // ---- Phase 3: Derive legacy face arrays from BC tables ----
      fault_interior_faces_.SetSize(0);
      fault_shared_faces_.SetSize(0);
      dirichlet_interior_faces_.SetSize(0);
      dirichlet_shared_faces_.SetSize(0);

      for (int f = 0; f < num_faces; f++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (!FTr) { continue; }
         if (face_bc_[f] == FacetBC::Fault)
         {
            fault_interior_faces_.Append(f);
         }
         else if (face_bc_[f] == FacetBC::Dirichlet)
         {
            dirichlet_interior_faces_.Append(f);
         }
      }
      for (int sf = 0; sf < num_shared; sf++)
      {
         if (shared_face_bc_[sf] == FacetBC::Fault)
         {
            fault_shared_faces_.Append(sf);
         }
         else if (shared_face_bc_[sf] == FacetBC::Dirichlet)
         {
            dirichlet_shared_faces_.Append(sf);
         }
      }

      // ---- Phase 4: Assert no None-classified face in any assembly array ----
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         MFEM_VERIFY(face_bc_[fault_interior_faces_[i]] == FacetBC::Fault,
            "ERROR: fault_interior_faces_[" << i << "] = "
            << fault_interior_faces_[i] << " has BC="
            << static_cast<int>(face_bc_[fault_interior_faces_[i]])
            << ", expected Fault.");
      }
      for (int i = 0; i < dirichlet_interior_faces_.Size(); i++)
      {
         MFEM_VERIFY(face_bc_[dirichlet_interior_faces_[i]] == FacetBC::Dirichlet,
            "ERROR: dirichlet_interior_faces_[" << i << "] = "
            << dirichlet_interior_faces_[i] << " has BC="
            << static_cast<int>(face_bc_[dirichlet_interior_faces_[i]])
            << ", expected Dirichlet.");
      }
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         MFEM_VERIFY(shared_face_bc_[fault_shared_faces_[i]] == FacetBC::Fault,
            "ERROR: fault_shared_faces_[" << i << "] = "
            << fault_shared_faces_[i] << " has BC="
            << static_cast<int>(shared_face_bc_[fault_shared_faces_[i]])
            << ", expected Fault.");
      }
      for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
      {
         MFEM_VERIFY(shared_face_bc_[dirichlet_shared_faces_[i]] == FacetBC::Dirichlet,
            "ERROR: dirichlet_shared_faces_[" << i << "] = "
            << dirichlet_shared_faces_[i] << " has BC="
            << static_cast<int>(shared_face_bc_[dirichlet_shared_faces_[i]])
            << ", expected Dirichlet.");
      }

      // ---- Phase 4b: Audit y=0 face classification ----
      // Detect interior faces on y=0 that are FacetBC::None (unclassified).
      // These faces get DG penalty enforcing zero jump, conflicting with
      // adjacent fault slip or Dirichlet displacement.
      {
         int rank = 0;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Comm_rank(mesh_.GetComm(), &rank);
#endif
         }

         int local_y0_none = 0, local_y0_fault = 0, local_y0_dir = 0;
         int local_y0_none_shared = 0;
         std::vector<std::array<double, 4>> none_faces;  // x, y, z, face_idx

         // Interior faces
         for (int f = 0; f < num_faces; f++)
         {
            auto *FTr = mesh_.GetInteriorFaceTransformations(f);
            if (!FTr) { continue; }
            const IntegrationPoint &ip =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ip);
            Vector fc(3);
            FTr->Face->Transform(ip, fc);
            if (std::abs(fc(1)) > 1.0) { continue; }  // not on y=0

            if (face_bc_[f] == FacetBC::None)
            {
               local_y0_none++;
               none_faces.push_back({fc(0), fc(1), fc(2),
                                     static_cast<double>(f)});
            }
            else if (face_bc_[f] == FacetBC::Fault) { local_y0_fault++; }
            else if (face_bc_[f] == FacetBC::Dirichlet) { local_y0_dir++; }
         }

         // Shared faces — count by classification only (avoid
         // GetSharedFaceTransformations which can fail for boundary faces).
         // Use GetSharedFace to get the local face index, then compute
         // centroid via the local face geometry.
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            for (int sf = 0; sf < num_shared; sf++)
            {
               if (shared_face_bc_[sf] != FacetBC::None) { continue; }

               // Get centroid of shared face via local face index
               int lf = mesh_.GetSharedFace(sf);
               auto *face_tr = mesh_.GetFaceTransformation(lf);
               if (!face_tr) { local_y0_none_shared++; continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(face_tr->GetGeometryType());
               face_tr->SetIntPoint(&ip);
               Vector fc(3);
               face_tr->Transform(ip, fc);
               if (std::abs(fc(1)) > 1.0) { continue; }

               local_y0_none_shared++;
               none_faces.push_back({fc(0), fc(1), fc(2),
                                     static_cast<double>(lf)});
            }
#endif
         }

         // Reduce counts
         int global_counts[4] = {local_y0_fault, local_y0_dir,
                                  local_y0_none, local_y0_none_shared};
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, global_counts, 4, MPI_INT,
                          MPI_SUM, mesh_.GetComm());
#endif
         }

         if (rank == 0)
         {
            std::cout << "\n=== Y=0 Face Classification Audit ===\n"
                      << "  Fault (attr=3):      " << global_counts[0] << "\n"
                      << "  Dirichlet (attr=5):  " << global_counts[1] << "\n"
                      << "  NONE (interior):     " << global_counts[2] << "\n"
                      << "  NONE (shared):       " << global_counts[3] << "\n";
            if (global_counts[2] + global_counts[3] > 0)
            {
               std::cout << "  *** WARNING: " << global_counts[2] + global_counts[3]
                         << " y=0 faces have NO BC! These enforce zero jump "
                         << "via DG penalty, conflicting with fault slip.\n";
            }
            else
            {
               std::cout << "  OK: all y=0 interior faces classified.\n";
            }
            std::cout << "=====================================\n\n";
         }

         // Dump unclassified face locations to CSV
         if (!none_faces.empty())
         {
            std::string fname = "face_audit_r" + std::to_string(rank) + ".csv";
            std::ofstream ofs(fname);
            ofs << "x,y,z,face_idx,bc\n";
            for (auto &nf : none_faces)
            {
               ofs << nf[0] << "," << nf[1] << ","
                   << nf[2] << "," << static_cast<int>(nf[3])
                   << ",None\n";
            }
            // Also dump fault faces near the boundary for context
            for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
            {
               int f = fault_interior_faces_[fi];
               auto *FTr = mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->Transform(ip, fc);
               ofs << fc(0) << "," << fc(1) << ","
                   << fc(2) << "," << f << ",Fault\n";
            }
            for (int di = 0; di < dirichlet_interior_faces_.Size(); di++)
            {
               int f = dirichlet_interior_faces_[di];
               auto *FTr = mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->Transform(ip, fc);
               ofs << fc(0) << "," << fc(1) << ","
                   << fc(2) << "," << f << ",Dirichlet\n";
            }
         }
      }

      // ---- Phase 5: Exact canonical-key validation ----
      ValidateFacetBCTables();
   }

   /// Allgather a local set of FaceVertexKeys and return the global union.
   std::set<FaceVertexKey> AllgatherKeys(
      const std::set<FaceVertexKey> &local_keys) const
   {
      std::set<FaceVertexKey> global_keys = local_keys;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int lc = static_cast<int>(local_keys.size());
         int nranks = 1;
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         std::vector<int> rc(nranks), dp(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, rc.data(), 1, MPI_INT,
                       mesh_.GetComm());
         int tot = 0;
         for (int r = 0; r < nranks; r++) { dp[r] = tot; tot += rc[r]; }
         std::vector<HYPRE_BigInt> lf(3 * lc);
         int idx = 0;
         for (auto &k : local_keys)
         {
            lf[3*idx] = k.v[0]; lf[3*idx+1] = k.v[1]; lf[3*idx+2] = k.v[2];
            idx++;
         }
         std::vector<int> rc3(nranks), dp3(nranks);
         for (int r = 0; r < nranks; r++)
         { rc3[r] = 3*rc[r]; dp3[r] = 3*dp[r]; }
         std::vector<HYPRE_BigInt> af(3 * tot);
         MPI_Allgatherv(lf.data(), 3*lc, HYPRE_MPI_BIG_INT,
                        af.data(), rc3.data(), dp3.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());
         for (int i = 0; i < tot; i++)
         {
            FaceVertexKey k;
            k.v[0] = af[3*i]; k.v[1] = af[3*i+1]; k.v[2] = af[3*i+2];
            global_keys.insert(k);
         }
#endif
      }
      return global_keys;
   }

   /// Exact canonical-key validation of recovered face sets.
   void ValidateFacetBCTables() const
   {
      Array<HYPRE_BigInt> gvert;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.GetGlobalVertexIndices(gvert);
#endif
      }
      else
      {
         gvert.SetSize(mesh_.GetNV());
         for (int i = 0; i < mesh_.GetNV(); i++) { gvert[i] = i; }
      }

      // Build tagged key sets from boundary elements (ground truth).
      // Separate tagged attr-5 faces into recoverable (interior/shared)
      // and exterior (handled by AddBdrFaceIntegrator, not recovered here).
      std::set<FaceVertexKey> local_tagged_fault;
      std::set<FaceVertexKey> local_tagged_dir_recoverable;  // interior/shared only
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr != 3 && attr != 5) { continue; }
         int face_idx = mesh_.GetBdrElementFaceIndex(be);
         FaceVertexKey key = MakeFaceKey(face_idx, gvert);
         if (attr == 3)
         {
            local_tagged_fault.insert(key);
         }
         else
         {
            // Is this face interior or shared? If so, it's recoverable.
            bool is_interior = (mesh_.GetInteriorFaceTransformations(face_idx)
                                != nullptr);
            bool is_shared = false;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               if (!is_interior)
               {
                  // Check if face_idx corresponds to a shared face
                  for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
                  {
                     if (mesh_.GetSharedFace(sf) == face_idx)
                     {
                        is_shared = true;
                        break;
                     }
                  }
               }
#endif
            }
            if (is_interior || is_shared)
            {
               local_tagged_dir_recoverable.insert(key);
            }
         }
      }
      std::set<FaceVertexKey> global_tagged_fault =
         AllgatherKeys(local_tagged_fault);
      std::set<FaceVertexKey> global_tagged_dir_recoverable =
         AllgatherKeys(local_tagged_dir_recoverable);

      // Build recovered key sets from the derived face arrays
      std::set<FaceVertexKey> local_recovered_fault, local_recovered_dir;
      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
         local_recovered_fault.insert(MakeFaceKey(fault_interior_faces_[fi], gvert));
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
         local_recovered_dir.insert(MakeFaceKey(dirichlet_interior_faces_[fi], gvert));
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
            local_recovered_fault.insert(
               MakeFaceKey(mesh_.GetSharedFace(fault_shared_faces_[i]), gvert));
         for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
            local_recovered_dir.insert(
               MakeFaceKey(mesh_.GetSharedFace(dirichlet_shared_faces_[i]), gvert));
#endif
      }
      std::set<FaceVertexKey> global_recovered_fault =
         AllgatherKeys(local_recovered_fault);
      std::set<FaceVertexKey> global_recovered_dir =
         AllgatherKeys(local_recovered_dir);

      // Exact fault key equality: tagged == recovered (bidirectional)
      for (auto &k : global_tagged_fault)
      {
         MFEM_VERIFY(global_recovered_fault.count(k) > 0,
            "ERROR: Tagged attr-3 face key (" << k.v[0] << "," << k.v[1]
            << "," << k.v[2] << ") was not recovered as a fault face.");
      }
      for (auto &k : global_recovered_fault)
      {
         MFEM_VERIFY(global_tagged_fault.count(k) > 0,
            "ERROR: Recovered fault face key (" << k.v[0] << "," << k.v[1]
            << "," << k.v[2] << ") has no matching attr-3 boundary element.");
      }

      // Exact Dirichlet key equality for recoverable faces (bidirectional)
      for (auto &k : global_tagged_dir_recoverable)
      {
         MFEM_VERIFY(global_recovered_dir.count(k) > 0,
            "ERROR: Tagged attr-5 interior/shared face key (" << k.v[0]
            << "," << k.v[1] << "," << k.v[2]
            << ") was not recovered as a Dirichlet face.");
      }
      for (auto &k : global_recovered_dir)
      {
         MFEM_VERIFY(global_tagged_dir_recoverable.count(k) > 0,
            "ERROR: Recovered Dirichlet face key (" << k.v[0] << ","
            << k.v[1] << "," << k.v[2]
            << ") has no matching recoverable attr-5 boundary element.");
      }

      // Fault ∩ Dirichlet = ∅
      for (auto &k : global_recovered_fault)
      {
         MFEM_VERIFY(!global_recovered_dir.count(k),
            "ERROR: Face key (" << k.v[0] << "," << k.v[1] << ","
            << k.v[2] << ") classified as both fault AND Dirichlet.");
      }

      // Summary
      bool is_root = true;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank; MPI_Comm_rank(mesh_.GetComm(), &rank);
         is_root = (rank == 0);
#endif
      }
      if (is_root)
      {
         mfem::out << "  Tag validation: "
                   << global_tagged_fault.size() << " fault keys (exact), "
                   << global_recovered_dir.size() << " Dirichlet keys (exact, "
                   << global_tagged_dir_recoverable.size() << " recoverable)"
                   << " (OK)\n";
         mfem::out << "  Fault faces: "
                   << fault_interior_faces_.Size() << " interior + "
                   << fault_shared_faces_.Size() << " shared (this rank)\n";
         mfem::out << "  Dirichlet faces: "
                   << dirichlet_interior_faces_.Size() << " interior + "
                   << dirichlet_shared_faces_.Size() << " shared (this rank)\n";
      }
   }

   void SetupFaultInfo()
   {
      // Build facet BC tables (single source of truth) and derive
      // legacy face arrays. Includes exact canonical-key validation.
      BuildFacetBCTables();

      num_fault_faces_ = fault_interior_faces_.Size() + fault_shared_faces_.Size();

      // Exchange face-neighbor data for shared face assembly
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.ExchangeFaceNbrData();
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (pfes) { pfes->ExchangeFaceNbrData(); }
#endif
      }

      // Multi-DOF fault quadrature
      // IP: use the same nodal triangle order as the volume space, matching
      // Tandem's fault discretization even at p=1.
      // BR2: keep the legacy face-averaged path (nbf=1).
      int face_fe_order = (method_ == DGMethod::IP) ? order_ : 0;
      face_quad_ = std::make_unique<FaceQuadrature>(face_fe_order,
                                                     std::max(order_, 1),
                                                     Geometry::TRIANGLE,
                                                     face_basis_type_);
      nbf_per_face_ = face_quad_->NumBasisFunctions();
      num_fault_dofs_ = num_fault_faces_ * nbf_per_face_;
      BuildOwnedFaultLayout();

      // Diagnostic: report owned vs local fault DOFs across all ranks
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0, nranks = 1;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         int local_shared = fault_shared_faces_.Size();
         int global_shared = 0;
         MPI_Reduce(&local_shared, &global_shared, 1, MPI_INT, MPI_SUM, 0,
                    mesh_.GetComm());
         int global_owned = 0;
         MPI_Reduce(&num_owned_fault_dofs_, &global_owned, 1, MPI_INT, MPI_SUM,
                    0, mesh_.GetComm());
         int global_local = 0;
         MPI_Reduce(&num_fault_dofs_, &global_local, 1, MPI_INT, MPI_SUM,
                    0, mesh_.GetComm());
         // Check how many faces have non-identity permutation
         int local_permuted = 0;
         for (int fi = 0; fi < num_fault_faces_; fi++)
         {
            for (int k = 0; k < nbf_per_face_; k++)
            {
               if (canonical_to_local_perm_[fi * nbf_per_face_ + k] != k)
               {
                  local_permuted++;
                  break;
               }
            }
         }
         int global_permuted = 0;
         MPI_Reduce(&local_permuted, &global_permuted, 1, MPI_INT, MPI_SUM, 0,
                    mesh_.GetComm());
         if (rank == 0)
         {
            mfem::out << "  Owned fault layout: "
                      << "global_shared_faces=" << global_shared
                      << ", global_owned_dofs=" << global_owned
                      << ", global_local_dofs=" << global_local
                      << ", faces_with_nontrivial_perm=" << global_permuted
                      << "\n";
         }
#endif
      }

      fault_dofs_.SetSize(num_fault_dofs_);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         fault_dofs_[i] = i;
      }

      // Compute FaultBasis for coordinate transforms
      if (num_fault_dofs_ > 0)
      {
         // Tandem coordinate system:
         //   X = along-strike, Y = fault-normal, Z = depth (negative down)
         //   Fault at Y = 0
         //   ref_normal = (0, -1, 0) matches Tandem's convention
         //   Up = (0, 0, 1)
         Vector ref_normal(3);
         ref_normal = 0.0;
         ref_normal(1) = -1.0;  // Y = fault-normal, pointing -Y

         // Store for skeleton Dirichlet sign computation
         ref_normal_ = ref_normal;

         Vector up(3);
         up = 0.0;
         up(2) = 1.0;  // +Z = upward

         fault_basis_.Compute(mesh_, fault_interior_faces_, ref_normal, up);

         // Append shared faces to FaultBasis
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            fault_basis_.AppendSharedFaces(mesh_, fault_shared_faces_,
                                           ref_normal, up);
#endif
         }

         // Per-quad-point basis (matching Tandem AdapterBase::prepare).
         // Uses the same 2p+1 face quadrature order as the combined integrator.
         // Detect face geometry from the first available fault face (interior
         // or shared), so ranks with only shared faces get the right rule.
         const int face_quad_order = 2 * order_ + 1;
         Geometry::Type face_geom = Geometry::TRIANGLE;  // default for tets
         if (fault_interior_faces_.Size() > 0)
         {
            auto *ftr0 = mesh_.GetInteriorFaceTransformations(
               fault_interior_faces_[0]);
            if (ftr0) { face_geom = ftr0->GetGeometryType(); }
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            if (fault_shared_faces_.Size() > 0)
            {
               auto *ftr0 = mesh_.GetSharedFaceTransformations(
                  fault_shared_faces_[0]);
               if (ftr0) { face_geom = ftr0->GetGeometryType(); }
            }
#endif
         }
         const IntegrationRule &face_ir =
            IntRules.Get(face_geom, face_quad_order);

         fault_basis_.ComputeQPBasis(mesh_, fault_interior_faces_,
                                      ref_normal, up, face_ir);

         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            fault_basis_.ComputeQPBasisShared(mesh_, fault_shared_faces_,
                                               ref_normal, up, face_ir,
                                               fault_interior_faces_.Size());
#endif
         }
      }
   }

   /// Startup face audit (validation-only). Verifies shared-face classification
   /// agreement between neighboring MPI ranks. No coordinate-based logic —
   /// only validates that the tag-recovered sets are self-consistent.
   void RunStartupFaceAudit()
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0, nranks = 1;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         MPI_Comm_size(mesh_.GetComm(), &nranks);

         // Build classification lookup for shared faces:
         // 0=unclassified, 1=fault, 2=dirichlet
         std::map<int, int> shared_cls;
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
            shared_cls[fault_shared_faces_[i]] = 1;
         for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
            shared_cls[dirichlet_shared_faces_[i]] = 2;

         // Pack: 4 values per classified shared face (3 key + cls)
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         std::vector<HYPRE_BigInt> local_flat;
         for (auto &[sf, cls] : shared_cls)
         {
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey fk = MakeFaceKey(lf, gvert);
            local_flat.push_back(fk.v[0]);
            local_flat.push_back(fk.v[1]);
            local_flat.push_back(fk.v[2]);
            local_flat.push_back(cls);
         }

         int lc = static_cast<int>(local_flat.size());
         std::vector<int> counts(nranks), displs(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, counts.data(), 1, MPI_INT,
                        mesh_.GetComm());
         int total = 0;
         for (int r = 0; r < nranks; r++)
         {
            displs[r] = total;
            total += counts[r];
         }
         std::vector<HYPRE_BigInt> all(total);
         MPI_Allgatherv(local_flat.data(), lc, HYPRE_MPI_BIG_INT,
                         all.data(), counts.data(), displs.data(),
                         HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Check: each key must appear exactly twice (once per rank sharing
         // the face), and both ranks must agree on the classification.
         std::map<FaceVertexKey, std::vector<int>> kmap;
         for (int i = 0; i < total; i += 4)
         {
            FaceVertexKey k;
            k.v[0] = all[i]; k.v[1] = all[i+1]; k.v[2] = all[i+2];
            kmap[k].push_back(static_cast<int>(all[i+3]));
         }
         int mismatch = 0, single_rank = 0;
         for (auto &[k, entries] : kmap)
         {
            if (entries.size() < 2)
            {
               // Face seen from only one rank — valid in some MFEM
               // partitioning layouts where only the owning rank
               // contributes the shared face to the allgather.
               single_rank++;
               continue;
            }
            // Check that all ranks contributing this face agree
            for (size_t e = 1; e < entries.size(); e++)
            {
               if (entries[e] != entries[0])
               {
                  mismatch++;
                  if (rank == 0 && mismatch <= 10)
                  {
                     mfem::out << "  AUDIT MISMATCH: key=("
                               << k.v[0] << "," << k.v[1] << "," << k.v[2]
                               << ") cls=" << entries[0] << " vs cls="
                               << entries[e] << "\n";
                  }
                  break;
               }
            }
         }
         MFEM_VERIFY(mismatch == 0,
            "ERROR: " << mismatch << " shared faces have classification "
            "mismatch between neighboring ranks.");
         if (rank == 0)
         {
            int paired = static_cast<int>(kmap.size()) - single_rank;
            mfem::out << "  AUDIT: shared-face agreement OK ("
                      << paired << " paired, "
                      << single_rank << " single-rank)\n";
         }
#endif
      }
   }

   void BuildOwnedFaultLayout()
   {
      owned_fault_face_to_local_face_.SetSize(0);
      full_fault_face_to_owned_face_.SetSize(num_fault_faces_);
      full_fault_face_to_owned_face_ = -1;
      owned_fault_dof_to_local_dof_.SetSize(0);
      shared_fault_comm_blocks_.clear();

      const int num_interior = fault_interior_faces_.Size();

      for (int face_idx = 0; face_idx < num_interior; face_idx++)
      {
         full_fault_face_to_owned_face_[face_idx] =
            owned_fault_face_to_local_face_.Size();
         owned_fault_face_to_local_face_.Append(face_idx);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);

         // Build canonical vertex-key for each local shared fault face
         std::vector<FaceVertexKey> local_shared_keys(fault_shared_faces_.Size());
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int sf = fault_shared_faces_[i];
            const int local_face = mesh_.GetSharedFace(sf);
            local_shared_keys[i] = MakeFaceKey(local_face, gvert);
         }

         // Allgather: collect all shared fault face keys + originating rank
         int comm_size = 1;
         MPI_Comm_size(mesh_.GetComm(), &comm_size);
         const int local_shared_count = fault_shared_faces_.Size();
         std::vector<int> recv_counts(comm_size, 0), recv_displs(comm_size, 0);
         MPI_Allgather(&local_shared_count, 1, MPI_INT,
                       recv_counts.data(), 1, MPI_INT, mesh_.GetComm());

         int total_shared = 0;
         for (int r = 0; r < comm_size; r++)
         {
            recv_displs[r] = total_shared;
            total_shared += recv_counts[r];
         }

         // Pack as flat array: 3 HYPRE_BigInt per face (sorted vertex IDs)
         std::vector<HYPRE_BigInt> local_flat(3 * local_shared_count);
         for (int i = 0; i < local_shared_count; i++)
         {
            local_flat[3*i]   = local_shared_keys[i].v[0];
            local_flat[3*i+1] = local_shared_keys[i].v[1];
            local_flat[3*i+2] = local_shared_keys[i].v[2];
         }
         std::vector<int> recv3(comm_size), disp3(comm_size);
         for (int r = 0; r < comm_size; r++)
         {
            recv3[r] = 3 * recv_counts[r];
            disp3[r] = 3 * recv_displs[r];
         }
         std::vector<HYPRE_BigInt> all_flat(3 * total_shared);
         MPI_Allgatherv(local_flat.data(), 3 * local_shared_count,
                        HYPRE_MPI_BIG_INT,
                        all_flat.data(), recv3.data(), disp3.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Build map: face key → list of ranks that have it
         std::map<FaceVertexKey, std::vector<int>> key_ranks;
         for (int r = 0; r < comm_size; r++)
         {
            for (int j = 0; j < recv_counts[r]; j++)
            {
               int idx = recv_displs[r] + j;
               FaceVertexKey k;
               k.v[0] = all_flat[3*idx];
               k.v[1] = all_flat[3*idx+1];
               k.v[2] = all_flat[3*idx+2];
               key_ranks[k].push_back(r);
            }
         }
         for (auto &kv : key_ranks)
         {
            auto &ranks = kv.second;
            std::sort(ranks.begin(), ranks.end());
            ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
         }

         std::unordered_map<int, std::vector<SharedFaultFaceBlock>> send_by_rank;
         std::unordered_map<int, std::vector<SharedFaultFaceBlock>> recv_by_rank;

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int sf = fault_shared_faces_[i];
            const int local_face = mesh_.GetSharedFace(sf);
            const int face_idx = num_interior + i;
            const FaceVertexKey key = local_shared_keys[i];
            auto key_it = key_ranks.find(key);
            MFEM_VERIFY(key_it != key_ranks.end() && !key_it->second.empty(),
                        "Missing shared fault ownership info for face key");
            const auto &sharing_ranks = key_it->second;
            const int owner_rank = sharing_ranks.front();

            if (rank == owner_rank)
            {
               const int owned_face =
                  owned_fault_face_to_local_face_.Size();
               full_fault_face_to_owned_face_[face_idx] = owned_face;
               owned_fault_face_to_local_face_.Append(face_idx);
               for (int other_rank : sharing_ranks)
               {
                  if (other_rank == rank) { continue; }
                  send_by_rank[other_rank].push_back({key, owned_face});
               }
            }
            else
            {
               recv_by_rank[owner_rank].push_back({key, face_idx});
            }
         }

         std::set<int> neighbors;
         for (const auto &kv : send_by_rank) { neighbors.insert(kv.first); }
         for (const auto &kv : recv_by_rank) { neighbors.insert(kv.first); }

         for (int neighbor : neighbors)
         {
            auto &send_faces = send_by_rank[neighbor];
            auto &recv_faces = recv_by_rank[neighbor];

            std::sort(send_faces.begin(), send_faces.end(),
                      [](const SharedFaultFaceBlock &a,
                         const SharedFaultFaceBlock &b)
                      {
                         return a.key < b.key;
                      });
            std::sort(recv_faces.begin(), recv_faces.end(),
                      [](const SharedFaultFaceBlock &a,
                         const SharedFaultFaceBlock &b)
                      {
                         return a.key < b.key;
                      });

            SharedFaultCommBlock block;
            block.neighbor_rank = neighbor;
            block.send_owned_faces.reserve(send_faces.size());
            block.recv_local_faces.reserve(recv_faces.size());
            for (const auto &entry : send_faces)
            {
               block.send_owned_faces.push_back(entry.face_idx);
            }
            for (const auto &entry : recv_faces)
            {
               block.recv_local_faces.push_back(entry.face_idx);
            }
            shared_fault_comm_blocks_.push_back(std::move(block));
         }
#endif
      }

      num_owned_fault_faces_ = owned_fault_face_to_local_face_.Size();
      num_owned_fault_dofs_ = num_owned_fault_faces_ * nbf_per_face_;

      // ------------------------------------------------------------------
      // Build per-face canonical DOF permutation (Tandem sorted-simplex).
      // For each fault face, sort the face vertices by global vertex ID.
      // canonical_to_local_perm_[face][k] = MFEM local DOF index for
      // the k-th vertex in sorted-global-ID order.
      // For p=1 triangles: DOF k = vertex k, so vertex perm = DOF perm.
      // ------------------------------------------------------------------
      canonical_to_local_perm_.SetSize(num_fault_faces_ * nbf_per_face_);

      // Get global vertex IDs (parallel) or use local indices (serial).
      // Use int64_t for sorting to avoid HYPRE dependency in serial builds.
      std::vector<int64_t> global_vert_ids;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvi;
         mesh_.GetGlobalVertexIndices(gvi);
         global_vert_ids.resize(gvi.Size());
         for (int i = 0; i < gvi.Size(); i++)
         {
            global_vert_ids[i] = static_cast<int64_t>(gvi[i]);
         }
#endif
      }
      else
      {
         global_vert_ids.resize(mesh_.GetNV());
         for (int i = 0; i < mesh_.GetNV(); i++)
         {
            global_vert_ids[i] = static_cast<int64_t>(i);
         }
      }

      for (int fi = 0; fi < num_fault_faces_; fi++)
      {
         // Get local face index in the mesh
         int local_face_idx;
         if (fi < num_interior)
         {
            local_face_idx = fault_interior_faces_[fi];
         }
         else
         {
            const int si = fi - num_interior;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
               local_face_idx = mesh_.GetSharedFace(fault_shared_faces_[si]);
            }
            else
            {
               MFEM_ABORT("Shared fault face in serial mesh");
               local_face_idx = -1;
            }
         }

         // Get face vertices (local indices)
         Array<int> vert;
         mesh_.GetFaceVertices(local_face_idx, vert);
         const int nv = vert.Size();

         // Build (global_id, mfem_local_index) pairs and sort by global ID
         // (Tandem sorted-simplex convention)
         std::vector<std::pair<int64_t, int>> gid_idx(nv);
         for (int k = 0; k < nv; k++)
         {
            gid_idx[k] = {global_vert_ids[vert[k]], k};
         }
         std::sort(gid_idx.begin(), gid_idx.end());

         // canonical_to_local_perm_[fi * nbf + canonical_k] = mfem_local_k
         if (nbf_per_face_ == 1)
         {
            // BR2: single centroid DOF per face, no vertex permutation
            canonical_to_local_perm_[fi] = 0;
         }
         else
         {
            // IP (nbf == nv for p=1): reorder by sorted global vertex IDs
            for (int k = 0; k < nv && k < nbf_per_face_; k++)
            {
               canonical_to_local_perm_[fi * nbf_per_face_ + k] =
                  gid_idx[k].second;
            }
            // For higher-order DOFs beyond vertices (p>=2): identity
            for (int k = nv; k < nbf_per_face_; k++)
            {
               canonical_to_local_perm_[fi * nbf_per_face_ + k] = k;
            }
         }
      }

      // Build owned DOF → local DOF map (with permutation applied)
      owned_fault_dof_to_local_dof_.SetSize(num_owned_fault_dofs_);
      for (int owned_face = 0; owned_face < num_owned_fault_faces_; owned_face++)
      {
         const int local_face = owned_fault_face_to_local_face_[owned_face];
         for (int kk = 0; kk < nbf_per_face_; kk++)
         {
            const int mfem_kk =
               canonical_to_local_perm_[local_face * nbf_per_face_ + kk];
            owned_fault_dof_to_local_dof_[owned_face * nbf_per_face_ + kk] =
               local_face * nbf_per_face_ + mfem_kk;
         }
      }
   }

   /// Facet BC table lookup for interior faces.
   FacetBC GetFaceBC(int face_idx) const
   {
      return (face_idx >= 0 && face_idx < static_cast<int>(face_bc_.size()))
             ? face_bc_[face_idx] : FacetBC::None;
   }

   /// Facet BC table lookup for shared faces.
   FacetBC GetSharedFaceBC(int shared_face) const
   {
      return (shared_face >= 0 &&
              shared_face < static_cast<int>(shared_face_bc_.size()))
             ? shared_face_bc_[shared_face] : FacetBC::None;
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

      // v57: Enable full neighbor-block assembly for DG shared faces.
      // With keep_nbr_block=false (MFEM default), shared face contributions
      // only assemble Elem1 rows into a rectangular local matrix. This causes
      // the global HypreParMatrix to have partition-dependent NNZ and slight
      // numerical differences in the stiffness matrix, producing ~6%
      // displacement errors across different MPI rank counts.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         cached_a_->KeepNbrBlock(true);
      }

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
         // v55: Combined integrator matching Tandem's assembleSurface kernel.
         // All three DG terms (consistency + symmetry + penalty) computed in
         // one pass with the same traction operator, same quadrature (2p+1),
         // and same penalty formula. K, fault slip b, Dirichlet b, and traction
         // recovery all use this same integrator.
         cached_a_->AddInteriorFaceIntegrator(
            new DGElasticityIPCombinedIntegrator(
               lambda_coeff_, mu_coeff_, 3, epsilon_, penalty_factor_));

         if (dirichlet_bdr_marker_.Size() > 0)
         {
            cached_a_->AddBdrFaceIntegrator(
               new DGElasticityIPCombinedIntegrator(
                  lambda_coeff_, mu_coeff_, 3, epsilon_, penalty_factor_),
               dirichlet_bdr_marker_);
         }
      }

      // v57: Assemble with skip_zeros=0 to ensure partition-independent NNZ.
      // With skip_zeros=1 (default), shared vs interior faces produce
      // slightly different face matrices (due to elem1/elem2 swap changing
      // floating-point order), causing O(ε) entries to be zero-skipped
      // on one partition but not another → different K matrix.
      cached_a_->Assemble(0);
      cached_a_->Finalize();

      // K assembly diagnostics: count faces/elements processed.
      // All ranks participate in MPI_Allreduce — MPI-safe.
      {
         long long local_volume = mesh_.GetNE();

         long long local_interior_face = 0;
         for (int f = 0; f < mesh_.GetNumFaces(); f++)
         {
            if (mesh_.GetInteriorFaceTransformations(f)) { local_interior_face++; }
         }

         // Count TRUE one-sided boundary faces (not interior) that K processes.
         // A boundary element is a true boundary face only if its underlying
         // mesh face does NOT have two local parent elements.
         // Use FaceIsInterior (checks Elem2No >= 0) instead of
         // GetBdrFaceTransformations (which can deadlock on ParMesh).
         long long local_bdr_face_true = 0;   // attr-5 true boundary faces
         long long local_bdr_all_true = 0;    // all-attr true boundary faces
         long long local_bdr_elem_total = mesh_.GetNBE();  // all boundary elements
         for (int be = 0; be < mesh_.GetNBE(); be++)
         {
            int face_idx, face_info_val;
            mesh_.GetBdrElementFace(be, &face_idx, &face_info_val);
            // Skip if the underlying face is interior (two local parents)
            // or shared (FaceIsTrueInterior includes shared faces)
            if (mesh_.FaceIsInterior(face_idx)) { continue; }
            local_bdr_all_true++;
            int attr = mesh_.GetBdrAttribute(be);
            if (dirichlet_bdr_marker_.Size() > 0 &&
                dirichlet_bdr_marker_[attr - 1] == 1)
            {
               local_bdr_face_true++;
            }
         }

         long long local_shared = 0;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            local_shared = mesh_.GetNSharedFaces();
#endif
         }

         long long counts[5] = {local_volume, local_interior_face, local_shared,
                                local_bdr_face_true, local_bdr_all_true};
         long long global[5] = {0};

         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(counts, global, 5, MPI_LONG_LONG, MPI_SUM,
                          mesh_.GetComm());
#endif
         }
         else
         {
            for (int i = 0; i < 5; i++) { global[i] = counts[i]; }
         }

         // Count true one-sided boundary faces per attribute.
         // Uses FaceIsTrueInterior to skip interior and shared faces
         // (same filter as the main boundary count above).
         int local_max_attr = mesh_.GetNBE() > 0 ? mesh_.bdr_attributes.Max() : 0;
         int max_attr = local_max_attr;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(&local_max_attr, &max_attr, 1, MPI_INT, MPI_MAX,
                          mesh_.GetComm());
#endif
         }
         std::vector<long long> local_bdr_by_attr(max_attr + 1, 0);
         for (int be = 0; be < mesh_.GetNBE(); be++)
         {
            int face_idx2, face_info_val2;
            mesh_.GetBdrElementFace(be, &face_idx2, &face_info_val2);
            if (mesh_.FaceIsInterior(face_idx2)) { continue; }
            int attr = mesh_.GetBdrAttribute(be);
            if (attr >= 1 && attr <= max_attr) { local_bdr_by_attr[attr]++; }
         }
         std::vector<long long> global_bdr_by_attr(max_attr + 1, 0);
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(local_bdr_by_attr.data(), global_bdr_by_attr.data(),
                          max_attr + 1, MPI_LONG_LONG, MPI_SUM, mesh_.GetComm());
#endif
         }
         else
         {
            global_bdr_by_attr = local_bdr_by_attr;
         }

         if (DebugRank() == 0)
         {
            mfem::out << "  [K-DIAG] Volume elements:            " << global[0] << "\n";
            mfem::out << "  [K-DIAG] Interior faces (K):         " << global[1] << "\n";
            mfem::out << "  [K-DIAG] Shared faces (K):           " << global[2]
                      << " (each face counted by both ranks)\n";
            mfem::out << "  [K-DIAG] Bdr faces w/ K (Dirichlet): " << global[3]
                      << " (true one-sided, attr 5)\n";
            mfem::out << "  [K-DIAG] Bdr faces (all attrs):      " << global[4]
                      << " (true one-sided, any attr)\n";
            mfem::out << "  [K-DIAG] Bdr elements total:         " << local_bdr_elem_total
                      << " (local rank, includes interior)\n";
            for (int a = 1; a <= max_attr; a++)
            {
               if (global_bdr_by_attr[a] > 0)
               {
                  mfem::out << "  [K-DIAG]   attr " << a << ": "
                            << global_bdr_by_attr[a] << " bdr faces"
                            << (dirichlet_bdr_marker_.Size() >= a &&
                                dirichlet_bdr_marker_[a-1] == 1
                                ? " (Dirichlet K)" : " (NO K integrator)")
                            << "\n";
               }
            }
            mfem::out << "  [K-DIAG] Total K contributions:  "
                      << global[0] + global[1] + global[2] + global[3]
                      << " (vol + interior + shared + bdr_dir)\n";
         }
      }

      // Set up solver operator
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         cached_Ah_.SetType(Operator::Hypre_ParCSR);
         cached_a_->ParallelAssemble(cached_Ah_);

         // K·v diagnostic: apply global K to all-ones vector,
         // dump result at target elements for cross-code comparison.
         if (first_step_debug_.enabled)
         {
            auto *K = cached_Ah_.As<HypreParMatrix>();
            int local_size = K->Height();
            HypreParVector ones(K->GetComm(), K->GetGlobalNumRows(),
                                K->GetRowStarts());
            ones = 1.0;
            HypreParVector Kv(K->GetComm(), K->GetGlobalNumRows(),
                              K->GetRowStarts());
            K->Mult(ones, Kv);

            // Dump K·1 for target elements
            if (DebugRank() == first_step_debug_.target_rank ||
                first_step_debug_.target_rank < 0)
            {
               std::string kv_path = DebugFilePath("first_step_Kv");
               std::ofstream out(kv_path, std::ios::trunc);
               out << "elem,component,local_dof,vdof,Kv_value\n";
               out << std::setprecision(17);

               const auto &target_elems = DebugTargetLocalElements();
               for (int e : target_elems)
               {
                  Array<int> vdofs;
                  fes_->GetElementVDofs(e, vdofs);
                  for (int j = 0; j < vdofs.Size(); j++)
                  {
                     int gj = vdofs[j];
                     int sign = 1;
                     if (gj < 0) { gj = -1 - gj; sign = -1; }
                     int comp = j / (vdofs.Size() / 3);
                     int ldof = j % (vdofs.Size() / 3);
                     double val = (gj < local_size) ? Kv(gj) * sign : 0.0;
                     out << e << "," << comp << "," << ldof << ","
                         << vdofs[j] << "," << val << "\n";
                  }
               }
               mfem::out << "  [K-DIAG] K·1 dumped to " << kv_path
                         << " (" << target_elems.size() << " elements)\n";
            }

            // Print global norms (collective, rank 0 prints)
            double local_norm1 = 0.0, local_norm2sq = 0.0, local_norminf = 0.0;
            for (int i = 0; i < local_size; i++)
            {
               double v = std::abs(Kv(i));
               local_norm1 += v;
               local_norm2sq += v * v;
               if (v > local_norminf) { local_norminf = v; }
            }
            double global_norm1, global_norm2sq, global_norminf;
            MPI_Reduce(&local_norm1, &global_norm1, 1, MPI_DOUBLE, MPI_SUM,
                        0, mesh_.GetComm());
            MPI_Reduce(&local_norm2sq, &global_norm2sq, 1, MPI_DOUBLE, MPI_SUM,
                        0, mesh_.GetComm());
            MPI_Reduce(&local_norminf, &global_norminf, 1, MPI_DOUBLE, MPI_MAX,
                        0, mesh_.GetComm());
            int myrank;
            MPI_Comm_rank(mesh_.GetComm(), &myrank);
            if (myrank == 0)
            {
               mfem::out << std::setprecision(15);
               mfem::out << "  [K-DIAG] ||K·1||_1   = " << global_norm1 << "\n";
               mfem::out << "  [K-DIAG] ||K·1||_2   = "
                         << std::sqrt(global_norm2sq) << "\n";
               mfem::out << "  [K-DIAG] ||K·1||_inf = " << global_norminf << "\n";
            }
         }

#ifdef MFEM_USE_MUMPS
         if (solver_type_ == SolverType::MUMPS)
         {
            auto *mumps = new MUMPSSolver(mesh_.GetComm());
            mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
            mumps->SetPrintLevel(1);
            mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(mumps);
         }
         else if (solver_type_ == SolverType::MUMPS_BLR)
         {
            auto *mumps = new MUMPSSolver(mesh_.GetComm());
            mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
            mumps->SetPrintLevel(1);
            mumps->SetBLRTol(blr_tol_);
            mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(mumps);
         }
         else
#endif
#ifdef MFEM_USE_SUPERLU
         if (solver_type_ == SolverType::SUPERLU)
         {
            cached_superlu_mat_.reset(
               new SuperLURowLocMatrix(*cached_Ah_.As<HypreParMatrix>()));
            auto *superlu = new SuperLUSolver(mesh_.GetComm());
            superlu->SetOperator(*cached_superlu_mat_);
            superlu->SetColumnPermutation(superlu::PARMETIS);
            superlu->SetIterativeRefine(superlu::SLU_DOUBLE);
            superlu->SetSymmetricPattern(true);
            superlu->SetPrintStatistics(true);
            solver_.reset(superlu);
         }
         else
#endif
#ifdef MFEM_USE_STRUMPACK
         if (solver_type_ == SolverType::STRUMPACK)
         {
            cached_strumpack_mat_.reset(
               new STRUMPACKRowLocMatrix(*cached_Ah_.As<HypreParMatrix>()));
            auto *strumpack = new STRUMPACKSolver(mesh_.GetComm());
            strumpack->SetOperator(*cached_strumpack_mat_);
            strumpack->SetCompression(strumpack::CompressionType::BLR);
            strumpack->SetCompressionRelTol(1e-10);
            strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
            strumpack->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);
            strumpack->SetPrintFactorStatistics(true);
            strumpack->SetPrintSolveStatistics(true);
            solver_.reset(strumpack);
         }
         else
#endif
         if (solver_type_ == SolverType::GMRES_BlockILU)
         {
            auto *gmres = new GMRESSolver(mesh_.GetComm());
            gmres->SetRelTol(1e-10);
            gmres->SetAbsTol(0.0);
            gmres->SetMaxIter(1000);
            gmres->SetKDim(50);
            gmres->SetPrintLevel(0);

            // BlockILU with element-sized blocks — natural for DG
            int block_size = fes_->GetTypicalFE()->GetDof() * 3;  // vdim=3
            cached_prec_.reset(new BlockILU(block_size,
               BlockILU::Reordering::MINIMUM_DISCARDED_FILL));

            gmres->SetPreconditioner(*cached_prec_);
            gmres->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(gmres);
         }
         else if (solver_type_ == SolverType::GMRES_AMG)
         {
            // GMRES + BoomerAMG: for DG at high p where:
            //   - MUMPS runs out of memory (dense DG factorization)
            //   - GMRES+BlockILU doesn't converge (local preconditioner too weak)
            //   - CG+AMG fails because AMG is non-SPD for DG matrices
            // GMRES tolerates non-SPD preconditioners while still benefiting
            // from AMG's multilevel global coarse-grid correction.
            auto *gmres = new GMRESSolver(mesh_.GetComm());
            gmres->SetRelTol(1e-10);
            gmres->SetAbsTol(0.0);
            gmres->SetMaxIter(2000);
            gmres->SetKDim(100);
            gmres->SetPrintLevel(1);

            auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
            amg->SetPrintLevel(0);
            // Skip SetElasticityOptions for DG — the near-null-space
            // setup assumes continuous FEM DOF connectivity and can
            // produce incorrect coarsening for DG sparsity patterns.
            // Plain AMG still provides effective multilevel preconditioning.
            cached_prec_.reset(amg);

            gmres->SetPreconditioner(*cached_prec_);
            gmres->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(gmres);
         }
         else
         {
            auto *cg = new CGSolver(mesh_.GetComm());
            cg->SetRelTol(1e-10);
            cg->SetAbsTol(0.0);
            cg->SetMaxIter(10000);
            cg->SetPrintLevel(0);

            auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
            amg->SetElasticityOptions(dynamic_cast<ParFiniteElementSpace*>(fes_.get()));
            amg->SetPrintLevel(0);
            cached_prec_.reset(amg);

            cg->SetPreconditioner(*cached_prec_);
            cg->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(cg);
         }
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
      int nbf = nbf_per_face_;

      // v55: Use the combined integrator's AssembleSlipFaceRHS to guarantee
      // K-b consistency. Same traction operator, quadrature, and penalty as K.
      DGElasticityIPCombinedIntegrator slip_integrator(
         lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         // Check if any slip is non-zero
         bool all_zero = true;
         for (int kk = 0; kk < nbf && all_zero; kk++)
         {
            int dof_idx = fi * nbf + kk;
            if (std::abs(slip_bc(2 * dof_idx)) >= 1e-15 ||
                std::abs(slip_bc(2 * dof_idx + 1)) >= 1e-15)
            { all_zero = false; }
         }
         if (all_zero) { continue; }

         // Build slip at quad points (Tandem evaluate_slip)
         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);

         // Assemble slip RHS using the combined integrator
         Vector elvec1, elvec2;
         slip_integrator.AssembleSlipFaceRHS(
            *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);
         DebugDumpFaceData(fi, face, "fault_interior", FTr, nullptr,
                           &delta_u_quad, elvec1, &elvec2);

         // Scatter into global RHS
         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

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
         int quad_order_br2 = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_br2);
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
            slip_sign_all[q] = (nor_q(1) > 0) ? 1.0 : -1.0;
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

         real_t c1 = epsilon_ * 0.5;  // symmetry data: −∫ g^F · {{σ(v)·n}} ds

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

            // Element 1 symmetry data + BR2 lifting data
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

            // Element 2 symmetry data + BR2 lifting data
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
         int nbf = nbf_per_face_;

         // v55: Use combined integrator for K-b consistency (shared faces)
         DGElasticityIPCombinedIntegrator slip_integrator(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int slip_idx = interior_face_count + i;
            bool all_zero = true;
            for (int kk = 0; kk < nbf && all_zero; kk++)
            {
               int dof_idx = slip_idx * nbf + kk;
               if (std::abs(slip_bc(2 * dof_idx)) >= 1e-15 ||
                   std::abs(slip_bc(2 * dof_idx + 1)) >= 1e-15)
               { all_zero = false; }
            }
            if (all_zero) { continue; }

            Vector delta_u_quad;
            BuildSlipAtQuadPoints(slip_idx, slip_bc, delta_u_quad);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

            // Assemble using combined integrator — only use elvec1 (local elem)
            Vector elvec1, elvec2;
            slip_integrator.AssembleSlipFaceRHS(
               *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);
            DebugDumpFaceData(slip_idx, mesh_.GetSharedFace(sf), "fault_shared",
                              FTr, nullptr, &delta_u_quad, elvec1, &elvec2);

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
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
            int quad_order_br2_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_br2_sh);
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
               slip_sign_all[q] = (nor_q(1) > 0) ? 1.0 : -1.0;
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

            real_t c1 = epsilon_ * 0.5;  // symmetry data: −∫ g^F · {{σ(v)·n}} ds

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
      // BP5 Dirichlet loading: u_D = (sgn(y)*Vp*t/2, 0, 0) on far-field,
      //                        u_D = (Vp*t, 0, 0) on y=0 skeleton faces.
      // Applied on attr-5 faces only (controlled by BCMode).
      //
      // DG Dirichlet BC contribution:
      //   b[k,i] += c0 * [σ(φ_k e_i)·n]_u * u_D_u * (1/detJ)
      //           + penalty * φ_k * u_D_i

      if (std::abs(time * Vp_) < 1e-30) { return; }

      int dim = 3;

      // Build set of face indices already handled as interior/shared Dirichlet.
      // These must NOT also be processed as boundary faces (would double-load
      // with incompatible DG formulas: boundary uses c1=epsilon, skeleton uses
      // c1=0.5*epsilon).
      std::set<int> dir_interior_set;
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
      {
         dir_interior_set.insert(dirichlet_interior_faces_[fi]);
      }
      std::set<int> dir_shared_set;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int fi = 0; fi < dirichlet_shared_faces_.Size(); fi++)
         {
            dir_shared_set.insert(mesh_.GetSharedFace(dirichlet_shared_faces_[fi]));
         }
#endif
      }

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

         // Get face transformation
         int face_idx, face_info;
         mesh_.GetBdrElementFace(be, &face_idx, &face_info);

         // Skip faces already handled as interior or shared Dirichlet
         if (dir_interior_set.count(face_idx) > 0) { continue; }
         if (dir_shared_set.count(face_idx) > 0) { continue; }

         FaceElementTransformations *FTr =
            mesh_.GetFaceElementTransformations(face_idx);
         if (FTr == nullptr) { continue; }

         Array<int> vdofs;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs);
         const FiniteElement *fe = scalar_fes_->GetFE(FTr->Elem1No);
         int ndof = fe->GetDof();

         int face_order = fe->GetOrder();
         int quad_order_bdr = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_bdr);

         Vector elvec(vdofs.Size());
         elvec = 0.0;

         if (method_ == DGMethod::IP)
         {
            // v55: Use combined integrator's AssembleBoundaryFaceRHS
            // (same formula as K's boundary face, Tandem's rhs_boundary)
            DGElasticityIPCombinedIntegrator dir_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int quad_order_dir = 2 * fe->GetOrder() + 1;
            const IntegrationRule &ir_dir = IntRules.Get(
               FTr->FaceGeom, quad_order_dir);
            int nq_dir = ir_dir.GetNPoints();

            // v55: Evaluate Tandem's boundary(x,y,z,t) at each quad point
            Vector u_D_3d(dim * nq_dir);
            u_D_3d = 0.0;
            Vector phys_y_qp(nq_dir);
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
               real_t y = phys(1);
               phys_y_qp(q) = y;
               // Tandem bp5.lua boundary(x,y,z,t):
               real_t Vh = Vp_ * time;
               if (y > 1000.0) { Vh *= 0.5; }
               else if (y < -1000.0) { Vh *= -0.5; }
               u_D_3d(0 * nq_dir + q) = Vh;
            }

            Vector elvec_dir;
            dir_integ.AssembleBoundaryFaceRHS(*fe, *FTr, u_D_3d, elvec_dir);
            DebugDumpFaceData(-1, face_idx, "dirichlet_boundary", FTr,
                              &phys_y_qp, &u_D_3d, elvec_dir, nullptr);

            for (int j = 0; j < elvec_dir.Size(); j++)
               elvec(j) += elvec_dir(j);
         }
         else  // BR2
         {
            // BR2: compute u_D from face centroid (legacy)
            const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_c);
            Vector fc_br2(dim);
            FTr->Face->Transform(ip_c, fc_br2);
            real_t Vh_br2 = Vp_ * time;
            if (fc_br2(1) > 1000.0) { Vh_br2 *= 0.5; }
            else if (fc_br2(1) < -1000.0) { Vh_br2 *= -0.5; }
            real_t u_D[3] = {Vh_br2, 0.0, 0.0};

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

      // ---------------------------------------------------------------
      // Interior faces with Dirichlet BC (skeleton pattern).
      // In Tandem's BP5 mesh, Physical Surface(5) includes Y=0 faces
      // outside the fault. These are interior faces in MFEM and need
      // the DG skeleton RHS contribution from both elements.
      // Tandem's bp5.lua boundary(): at Y=0, neither y>1 nor y<-1
      // triggers, so u_D = (Vp*t, 0, 0) — the full plate velocity.
      // This represents a locked fault (no relative slip).
      // ---------------------------------------------------------------
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
      {
         int f = dirichlet_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (FTr == nullptr) { continue; }

         // Evaluate bp5.lua boundary(x,y,z,t) at each quad point:
         //   y > 1000 m (1 km):  u_D = (Vp*t/2, 0, 0)
         //   y < -1000 m:        u_D = (-Vp*t/2, 0, 0)
         //   |y| <= 1000 m:      u_D = (Vp*t, 0, 0)
         //
         // Orientation sign applied per QP via ComputeSkeletonDirichletSign,
         // matching Tandem DGCurvilinearCommon.h:92-99.

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         int ndof1 = fe1->GetDof();
         int ndof2 = fe2->GetDof();

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         int quad_order_dir = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_dir);

         Vector elvec1(vdofs1.Size()), elvec2(vdofs2.Size());
         elvec1 = 0.0;
         elvec2 = 0.0;

         if (method_ == DGMethod::IP)
         {
            DGElasticityIPCombinedIntegrator dir_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_dir = IntRules.Get(
               FTr->FaceGeom, quad_order_dir);
            int nq_dir = ir_dir.GetNPoints();

            // Evaluate boundary(x,y,z,t) at each QP with per-QP orientation
            // sign, matching Tandem DGCurvilinearCommon.h:92-99.
            Vector u_D_3d(dim * nq_dir);
            u_D_3d = 0.0;
            Vector phys_y_qp(nq_dir);
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
               real_t y = phys(1);
               phys_y_qp(q) = y;
               real_t Vh = Vp_ * time;
               if (y > 1000.0) { Vh *= 0.5; }
               else if (y < -1000.0) { Vh *= -0.5; }
               // Per-QP orientation sign (Tandem DGCurvilinearCommon.h:97-98)
               real_t dir_sign = ComputeSkeletonDirichletSign(FTr);
               u_D_3d(0 * nq_dir + q) = dir_sign * Vh;
            }

            Vector ev1, ev2;
            dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);
            DebugDumpFaceData(-1, f, "dirichlet_interior", FTr,
                              &phys_y_qp, &u_D_3d, ev1, &ev2);

            for (int j = 0; j < ev1.Size(); j++) { elvec1(j) += ev1(j); }
            for (int j = 0; j < ev2.Size(); j++) { elvec2(j) += ev2(j); }
         }
         else  // BR2
         {
            // BR2 path: compute u_D_int from face centroid (legacy)
            const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_c);
            Vector fc_br2(dim);
            FTr->Face->Transform(ip_c, fc_br2);
            real_t Vh_br2 = Vp_ * time;
            if (fc_br2(1) > 1000.0) { Vh_br2 *= 0.5; }
            else if (fc_br2(1) < -1000.0) { Vh_br2 *= -0.5; }
            // Orientation sign for BR2 (face-constant, affine faces)
            const IntegrationPoint &ip_s = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ip_s);
            real_t dir_sign_br2 = ComputeSkeletonDirichletSign(FTr);
            real_t u_D_int[3] = {dir_sign_br2 * Vh_br2, 0.0, 0.0};

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);

               Vector sq1(shapes1.GetColumn(q), ndof1);
               fe1->CalcShape(FTr->GetElement1IntPoint(), sq1);

               Vector sq2(shapes2.GetColumn(q), ndof2);
               fe2->CalcShape(FTr->GetElement2IntPoint(), sq2);

               Vector nq(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nq);
               w_arr[q] = ip.weight;
            }

            // Compute lifted Dirichlet for skeleton:
            // face_int_k[u*dim+s, m] = sum_q shape_k_m(q) * u_D[u] * n[s] * w[q] * 0.5
            //
            // CRITICAL: Both face_int1 and face_int2 use the SAME sign (+=).
            // This matches the bilinear form (DGElasticityBR2Integrator) and the
            // slip assembly (AssembleSlipContributionBR2), which both use:
            //   IntFace_11 += shapes1 * factor
            //   IntFace_21 += shapes2 * factor  (same sign)
            //
            // The u_D_int is a prescribed JUMP (u1-u2), not a per-element value.
            // The BR2 lifting of the jump into both elements uses the same
            // reference normal. The -sign for elem2's penalty is handled later
            // via "elvec2 -= penalty * ..." (not here in face_int).
            //
            // BUG FIX (v37): Previously face_int2 used -= (opposite sign), which:
            //   - In v34 (no cross-element): pushed elem2 in WRONG direction
            //     (resisting correct loading, causing delayed recurrence)
            //   - In v36 (cross-element + 0.5): caused eval1+eval2 cancellation
            //     → zero penalty → fault lockup
            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_arr[q];
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_arr[q * dim + s];
                     real_t factor = u_D_int[u] * n_s * wq * 0.5;
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
            MultABt(face_int2, Minv2, f_lifted2);

            // BR2 penalty
            Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            // Combined BR2 lifted function evaluation at each quadrature point.
            // Uses cross-element lifting: f_lifted_q = C:n:(r_1 + r_2)
            // where r_k = Minv_k * face_int_k is the lift into element k.
            // The 0.5 factor matches the bilinear form's L_q = 0.5 * sum
            // structure (face_int already contains the {ψ} average factor 0.5).
            DenseMatrix f_lifted_q(dim, nqp);
            f_lifted_q = 0.0;
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

                        // Cross-element evaluation: sum over BOTH elements
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

            // Elem1 contribution
            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof1, dim);
               fe1->CalcDShape(eip1, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof1, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem1->Weight();

               for (int k = 0; k < ndof1; k++)
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
                        sym_val += trac * u_D_int[u];
                     }

                     int idx = i * ndof1 + k;
                     elvec1(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                     elvec1(idx) += br2_pen * wq * shapes1(k, q) * f_lifted_q(i, q);
                  }
               }
            }

            // Elem2 contribution (opposite penalty sign)
            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof2, dim);
               fe2->CalcDShape(eip2, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem2->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof2, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem2->Weight();

               for (int k = 0; k < ndof2; k++)
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
                        sym_val += trac * u_D_int[u];
                     }

                     int idx = i * ndof2 + k;
                     elvec2(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                     elvec2(idx) -= br2_pen * wq * shapes2(k, q) * f_lifted_q(i, q);
                  }
               }
            }
         }

         // Scatter to global RHS
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

      // ---------------------------------------------------------------
      // Shared faces with Dirichlet BC (parallel only).
      // Same skeleton pattern as interior Dirichlet faces, but only
      // elem1 is local — elem2 is on a neighboring rank.
      // Each rank contributes its local elem1 side; the neighboring rank
      // handles that same face's elem2 side as ITS elem1.
      // ---------------------------------------------------------------
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

         for (int fi = 0; fi < dirichlet_shared_faces_.Size(); fi++)
         {
            int sf = dirichlet_shared_faces_[fi];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            // v55: Per-quad-point boundary evaluation (done below in IP block)

            // Only elem1 is local
            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            int ndof1 = fe1->GetDof();

            // Elem2 (face-neighbor)
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            int ndof2 = fe2->GetDof();

            int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
            int quad_order_dir_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_dir_sh);

            Vector elvec1(vdofs1.Size());
            elvec1 = 0.0;

            if (method_ == DGMethod::IP)
            {
               DGElasticityIPCombinedIntegrator dir_integ(
                  lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

               int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
               const IntegrationRule &ir_dir = IntRules.Get(
                  FTr->FaceGeom, quad_order_dir);
               int nq_dir = ir_dir.GetNPoints();

               // Per-QP boundary evaluation with orientation sign
               Vector u_D_3d(dim * nq_dir);
               u_D_3d = 0.0;
               Vector phys_y_qp(nq_dir);
               for (int q = 0; q < nq_dir; q++)
               {
                  const IntegrationPoint &ipq = ir_dir.IntPoint(q);
                  FTr->SetAllIntPoints(&ipq);
                  Vector phys(dim);
                  FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
                  real_t y = phys(1);
                  phys_y_qp(q) = y;
                  real_t Vh = Vp_ * time;
                  if (y > 1000.0) { Vh *= 0.5; }
                  else if (y < -1000.0) { Vh *= -0.5; }
                  real_t dir_sign = ComputeSkeletonDirichletSign(FTr);
                  u_D_3d(0 * nq_dir + q) = dir_sign * Vh;
               }

               Vector ev1, ev2;
               dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);
               DebugDumpFaceData(-1, mesh_.GetSharedFace(sf), "dirichlet_shared",
                                 FTr, &phys_y_qp, &u_D_3d, ev1, &ev2);
               for (int j = 0; j < ev1.Size(); j++) { elvec1(j) += ev1(j); }
               // ev2 goes to neighbor rank, not used here
            }
            else  // BR2
            {
               if (!mass_inv_computed_) { PrecomputeMassInverse(); }

               // BR2: compute u_D_int from face centroid (legacy)
               const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
               FTr->Face->SetIntPoint(&ip_c);
               Vector fc_br2(dim);
               FTr->Face->Transform(ip_c, fc_br2);
               real_t Vh_br2 = Vp_ * time;
               if (fc_br2(1) > 1000.0) { Vh_br2 *= 0.5; }
               else if (fc_br2(1) < -1000.0) { Vh_br2 *= -0.5; }
               // Orientation sign for BR2 (face-constant, affine faces)
               FTr->SetAllIntPoints(&ip_c);
               real_t dir_sign_br2 = ComputeSkeletonDirichletSign(FTr);
               real_t u_D_int[3] = {dir_sign_br2 * Vh_br2, 0.0, 0.0};

               const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
               const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];
               int nqp = ir.GetNPoints();

               // Precompute shapes and normals
               DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
               Vector nor_arr(dim * nqp), w_arr(nqp);

               for (int q = 0; q < nqp; q++)
               {
                  const IntegrationPoint &ip = ir.IntPoint(q);
                  FTr->SetAllIntPoints(&ip);

                  Vector sq1(shapes1.GetColumn(q), ndof1);
                  fe1->CalcShape(FTr->GetElement1IntPoint(), sq1);

                  Vector sq2(shapes2.GetColumn(q), ndof2);
                  fe2->CalcShape(FTr->GetElement2IntPoint(), sq2);

                  Vector nq(&nor_arr[q * dim], dim);
                  CalcOrtho(FTr->Jacobian(), nq);
                  w_arr[q] = ip.weight;
               }

               // Compute lifted Dirichlet for skeleton (both sides)
               DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
               face_int1 = 0.0;
               face_int2 = 0.0;
               for (int q = 0; q < nqp; q++)
               {
                  real_t wq = w_arr[q];
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t n_s = nor_arr[q * dim + s];
                        real_t factor = u_D_int[u] * n_s * wq * 0.5;
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
               MultABt(face_int2, Minv2, f_lifted2);

               // BR2 penalty
               Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
               real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                    ? real_t(dim + 1) : real_t(2 * dim);

               // Cross-element lifted function evaluation
               DenseMatrix f_lifted_q(dim, nqp);
               f_lifted_q = 0.0;
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

               // Elem1 contribution only
               for (int q = 0; q < nqp; q++)
               {
                  const IntegrationPoint &ip = ir.IntPoint(q);
                  FTr->SetAllIntPoints(&ip);
                  const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

                  real_t wq = w_arr[q];
                  Vector nor_q(dim);
                  for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

                  DenseMatrix dshape_ref(ndof1, dim);
                  fe1->CalcDShape(eip1, dshape_ref);
                  DenseMatrix adjJ(dim);
                  CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
                  DenseMatrix dshape_adj(ndof1, dim);
                  Mult(dshape_ref, adjJ, dshape_adj);
                  real_t detJ = FTr->Elem1->Weight();

                  for (int k = 0; k < ndof1; k++)
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
                           sym_val += trac * u_D_int[u];
                        }

                        int idx = i * ndof1 + k;
                        elvec1(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                        elvec1(idx) += br2_pen * wq * shapes1(k, q) * f_lifted_q(i, q);
                     }
                  }
               }
            }

            // Scatter elem1 to global RHS
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
      fault_depths_ = 0.0;

      // Nodal rule for per-DOF coordinate evaluation.
      // At nbf=1 (BR2 / order-0 face space): single centroid point.
      // At nbf>1 (IP nodal fault space): face nodes on the reference triangle.
      const IntegrationRule &nir = face_quad_->GetNodalRule();
      int nbf = nbf_per_face_;

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         for (int kk = 0; kk < nbf; kk++)
         {
            const IntegrationPoint &nip = nir.IntPoint(kk);
            FTr->Face->SetIntPoint(&nip);
            Vector coords(3);
            FTr->Face->Transform(nip, coords);

            // Depth: -Z (Z is negative downward in Tandem, depth is positive)
            fault_depths_(i * nbf + kk) = -coords(2);
         }
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

            int face_idx = fault_interior_faces_.Size() + i;
            for (int kk = 0; kk < nbf; kk++)
            {
               const IntegrationPoint &nip = nir.IntPoint(kk);
               FTr->Face->SetIntPoint(&nip);
               Vector coords(3);
               FTr->Face->Transform(nip, coords);

               fault_depths_(face_idx * nbf + kk) = -coords(2);
            }
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
      fault_x2_ = 0.0;
      fault_x3_ = 0.0;

      // Nodal rule for per-DOF coordinate evaluation
      const IntegrationRule &nir = face_quad_->GetNodalRule();
      int nbf = nbf_per_face_;

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         for (int kk = 0; kk < nbf; kk++)
         {
            const IntegrationPoint &nip = nir.IntPoint(kk);

            // Map face integration point to physical coordinates via Face
            FTr->Face->SetIntPoint(&nip);
            Vector coords(3);
            FTr->Face->Transform(nip, coords);

            // Tandem: X=along-strike=coords(0), depth=-Z=-coords(2)
            fault_x2_(i * nbf + kk) = coords(0);
            fault_x3_(i * nbf + kk) = -coords(2);
         }
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

            int face_idx = fault_interior_faces_.Size() + i;
            for (int kk = 0; kk < nbf; kk++)
            {
               const IntegrationPoint &nip = nir.IntPoint(kk);
               FTr->Face->SetIntPoint(&nip);
               Vector coords(3);
               FTr->Face->Transform(nip, coords);

               fault_x2_(face_idx * nbf + kk) = coords(0);
               fault_x3_(face_idx * nbf + kk) = -coords(2);
            }
         }
#endif
      }

      fault_coords_computed_ = true;
   }

   coords_x2 = fault_x2_;
   coords_x3 = fault_x3_;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::RestrictToOwnedFault(
   const Vector &local_data, Vector &owned_data, int comps_per_dof) const
{
   MFEM_VERIFY(comps_per_dof > 0, "comps_per_dof must be positive");
   MFEM_VERIFY(local_data.Size() == comps_per_dof * num_fault_dofs_,
               "Local fault vector size mismatch: got " << local_data.Size()
               << ", expected " << comps_per_dof * num_fault_dofs_);

   owned_data.SetSize(comps_per_dof * num_owned_fault_dofs_);
   for (int owned_dof = 0; owned_dof < num_owned_fault_dofs_; owned_dof++)
   {
      const int local_dof = owned_fault_dof_to_local_dof_[owned_dof];
      for (int c = 0; c < comps_per_dof; c++)
      {
         owned_data(comps_per_dof * owned_dof + c) =
            local_data(comps_per_dof * local_dof + c);
      }
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ExpandOwnedToLocalFault(
   const Vector &owned_data, Vector &local_data, int comps_per_dof) const
{
   MFEM_VERIFY(comps_per_dof > 0, "comps_per_dof must be positive");
   MFEM_VERIFY(owned_data.Size() == comps_per_dof * num_owned_fault_dofs_,
               "Owned fault vector size mismatch: got " << owned_data.Size()
               << ", expected " << comps_per_dof * num_owned_fault_dofs_);

   local_data.SetSize(comps_per_dof * num_fault_dofs_);
   local_data = 0.0;

   for (int owned_dof = 0; owned_dof < num_owned_fault_dofs_; owned_dof++)
   {
      const int local_dof = owned_fault_dof_to_local_dof_[owned_dof];
      for (int c = 0; c < comps_per_dof; c++)
      {
         local_data(comps_per_dof * local_dof + c) =
            owned_data(comps_per_dof * owned_dof + c);
      }
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      if (shared_fault_comm_blocks_.empty()) { return; }

      const int block_size = comps_per_dof * nbf_per_face_;
      std::vector<Vector> send_buffers(shared_fault_comm_blocks_.size());
      std::vector<Vector> recv_buffers(shared_fault_comm_blocks_.size());
      std::vector<MPI_Request> requests;
      requests.reserve(2 * shared_fault_comm_blocks_.size());

      for (int bi = 0; bi < static_cast<int>(shared_fault_comm_blocks_.size()); bi++)
      {
         const auto &block = shared_fault_comm_blocks_[bi];

         if (!block.recv_local_faces.empty())
         {
            recv_buffers[bi].SetSize(block_size * block.recv_local_faces.size());
            MPI_Request req;
            MPI_Irecv(recv_buffers[bi].GetData(), recv_buffers[bi].Size(),
                      MPI_DOUBLE, block.neighbor_rank, 27183,
                      mesh_.GetComm(), &req);
            requests.push_back(req);
         }

         if (!block.send_owned_faces.empty())
         {
            send_buffers[bi].SetSize(block_size * block.send_owned_faces.size());
            for (int j = 0; j < static_cast<int>(block.send_owned_faces.size()); j++)
            {
               const int owned_face = block.send_owned_faces[j];
               for (int kk = 0; kk < nbf_per_face_; kk++)
               {
                  const int owned_dof = owned_face * nbf_per_face_ + kk;
                  for (int c = 0; c < comps_per_dof; c++)
                  {
                     send_buffers[bi](j * block_size + kk * comps_per_dof + c) =
                        owned_data(comps_per_dof * owned_dof + c);
                  }
               }
            }

            MPI_Request req;
            MPI_Isend(send_buffers[bi].GetData(), send_buffers[bi].Size(),
                      MPI_DOUBLE, block.neighbor_rank, 27183,
                      mesh_.GetComm(), &req);
            requests.push_back(req);
         }
      }

      if (!requests.empty())
      {
         MPI_Waitall(static_cast<int>(requests.size()), requests.data(),
                     MPI_STATUSES_IGNORE);
      }

      for (int bi = 0; bi < static_cast<int>(shared_fault_comm_blocks_.size()); bi++)
      {
         const auto &block = shared_fault_comm_blocks_[bi];
         const Vector &recv = recv_buffers[bi];
         if (recv.Size() == 0) { continue; }

         for (int j = 0; j < static_cast<int>(block.recv_local_faces.size()); j++)
         {
            const int local_face = block.recv_local_faces[j];
            // Received data is in canonical order (sorted by global vertex ID).
            // Apply this rank's canonical→local permutation so the values
            // land at the correct MFEM-local DOF positions.
            for (int canonical_kk = 0; canonical_kk < nbf_per_face_; canonical_kk++)
            {
               const int mfem_kk =
                  canonical_to_local_perm_[local_face * nbf_per_face_ + canonical_kk];
               const int local_dof = local_face * nbf_per_face_ + mfem_kk;
               for (int c = 0; c < comps_per_dof; c++)
               {
                  local_data(comps_per_dof * local_dof + c) =
                     recv(j * block_size + canonical_kk * comps_per_dof + c);
               }
            }
         }
      }
#endif
   }
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

   const bool debug_first_step = (method_ == DGMethod::IP) &&
                                 DebugEnabledForTime(time);
   debug_time_ = time;

   // Production assembly path: same on ALL ranks (MPI-safe).
   // Debug dumps happen AFTER the solve, gated by rank — no separate
   // assembly path that could cause MPI divergence.
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

   // Snapshot slip-only RHS before adding Dirichlet.
   // Taken on ALL ranks when debug is enabled (needed for MPI-collective norms).
   Vector rhs_slip_snapshot;
   const bool debug_norm = first_step_debug_.enabled && (time > 0.0);
   if (debug_first_step || debug_norm)
   {
      rhs_slip_snapshot = rhs;  // copy before Dirichlet is added
   }

   // Add Dirichlet loading
   AssembleDirichletLoading(rhs, time);

   X_ = 0.0;
   B_ = rhs;

   solver_->Mult(B_, X_);

   // Post-solve residual check: ||K*x - b|| / ||b|| (global norms)
   if (check_residual_)
   {
      Vector R_(B_.Size());
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         cached_Ah_.As<HypreParMatrix>()->Mult(X_, R_);
         R_ -= B_;
         // Use global norms (MPI AllReduce) for meaningful relative residual
         real_t local_res2 = R_ * R_;
         real_t local_rhs2 = B_ * B_;
         real_t global_res2 = 0.0, global_rhs2 = 0.0;
         MPI_Allreduce(&local_res2, &global_res2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&local_rhs2, &global_rhs2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         real_t global_res = std::sqrt(global_res2);
         real_t global_rhs = std::sqrt(global_rhs2);
         real_t rel_res = (global_rhs > 0.0) ? global_res / global_rhs : global_res;
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         if (rel_res > 1e-8 && rank == 0)
         {
            mfem::err << "RESIDUAL WARNING: global ||K*x-b||/||b|| = " << rel_res
                      << " (threshold 1e-8)\n";
         }
#endif
      }
      else
      {
         cached_a_->SpMat().Mult(X_, R_);
         R_ -= B_;
         real_t res_norm = R_.Norml2();
         real_t rhs_norm = B_.Norml2();
         real_t rel_res = (rhs_norm > 0.0) ? res_norm / rhs_norm : res_norm;
         if (rel_res > 1e-8)
         {
            mfem::err << "RESIDUAL WARNING: ||K*x-b||/||b|| = " << rel_res
                      << " (threshold 1e-8)\n";
         }
      }
   }

   // Check for NaN/Inf in solution
   {
      real_t u_max = X_.Normlinf();
      MFEM_VERIFY(std::isfinite(u_max),
         "Domain solve produced NaN/Inf in displacement (||u||_inf = "
         << u_max << ")");
   }

   // Check convergence and log solver info when RHS is large
   {
      auto *iter_solver = dynamic_cast<IterativeSolver*>(solver_.get());
      if (iter_solver)
      {
         if (!iter_solver->GetConverged())
         {
            mfem::err << "WARNING: iterative solver did not converge after "
                      << iter_solver->GetNumIterations() << " iterations, final norm = "
                      << iter_solver->GetFinalNorm() << "\n";
         }
      }
   }

   displacement = X_;

   // Debug dumps: rank-local file I/O only, no MPI collectives.
   // Placed AFTER displacement = X_ so the GridFunction has the current solution.
   if (debug_first_step)
   {
      // Compute Dirichlet-only RHS by subtraction (no extra assembly)
      Vector rhs_dir(rhs.Size());
      subtract(rhs, rhs_slip_snapshot, rhs_dir);

      DebugDumpElementVector("rhs_slip", rhs_slip_snapshot);
      DebugDumpElementVector("rhs_dirichlet", rhs_dir);
      DebugDumpElementVector("rhs_total", rhs);
      DebugDumpElementVector("u", X_);

      // Dump fault jumps and traction using the solved displacement.
      // Interior faces only — no ExchangeFaceNbrData (MPI-safe).
      DebugDumpFaultJumpsLocal(displacement, slip_bc);
      DebugDumpFaultTractionLocal(displacement, slip_bc);
      DebugDumpKContributions();
      first_step_debug_done_ = true;
   }

   // Global norm comparison — OUTSIDE rank-gated block.
   // All ranks must participate in MPI_Reduce (no deadlock).
   // Condition uses only rank-independent flags (enabled + time).
   // Prints at EVERY Mult() call during first step to trace RK45 stages.
   {
      static int norm_count = 0;
      if (norm_count < 10 && first_step_debug_.enabled && time > 0.0)
      {
         norm_count++;
         int myrank;
         MPI_Comm_rank(mesh_.GetComm(), &myrank);
         if (myrank == 0)
         {
            mfem::out << "  [NORM] Stage " << norm_count
                      << " at time = " << std::setprecision(17) << time << "\n";
         }

         auto print_global_norm = [&](const char *label, const Vector &v) {
            double local_n1 = 0.0, local_n2sq = 0.0, local_ninf = 0.0;
            for (int i = 0; i < v.Size(); i++)
            {
               double a = std::abs(v(i));
               local_n1 += a;
               local_n2sq += a * a;
               if (a > local_ninf) { local_ninf = a; }
            }
            double g1, g2sq, ginf;
            MPI_Reduce(&local_n1, &g1, 1, MPI_DOUBLE, MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_n2sq, &g2sq, 1, MPI_DOUBLE, MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_ninf, &ginf, 1, MPI_DOUBLE, MPI_MAX, 0, mesh_.GetComm());
            if (myrank == 0)
            {
               mfem::out << std::setprecision(15);
               mfem::out << "  [NORM] ||" << label << "||_1   = " << g1 << "\n";
               mfem::out << "  [NORM] ||" << label << "||_2   = "
                         << std::sqrt(g2sq) << "\n";
               mfem::out << "  [NORM] ||" << label << "||_inf = " << ginf << "\n";
            }
         };

         print_global_norm("b_slip", rhs_slip_snapshot);
         print_global_norm("b_total", rhs);
         print_global_norm("u", X_);
      }
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector *normal_traction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       nullptr, nullptr, nullptr);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::AssembleSlipOnlyRHS(
   Vector &rhs, const Vector &slip_bc) const
{
   rhs.SetSize(fes_->GetVSize());
   rhs = 0.0;

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
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::AssembleDirichletOnlyRHS(
   Vector &rhs, real_t time) const
{
   rhs.SetSize(fes_->GetVSize());
   rhs = 0.0;
   AssembleDirichletLoading(rhs, time);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionComponents(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector &traction_stress,
   Vector &traction_correction,
   Vector *normal_traction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       &traction_stress, &traction_correction, nullptr);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionDiagnostics(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector &traction_stress,
   Vector &traction_correction,
   Vector &jump_residual,
   Vector *normal_traction,
   Vector *normal_stress,
   Vector *normal_correction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       &traction_stress, &traction_correction, &jump_residual,
                       normal_stress, normal_correction);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionImpl(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector *normal_traction,
   Vector *traction_stress_out,
   Vector *traction_correction_out,
   Vector *jump_residual_out,
   Vector *normal_stress_out,
   Vector *normal_correction_out)
{
   int dim = 3;
   traction.SetSize(2 * num_fault_dofs_);
   traction = 0.0;

    if (traction_stress_out)
    {
       traction_stress_out->SetSize(2 * num_fault_dofs_);
       *traction_stress_out = 0.0;
    }
    if (traction_correction_out)
    {
       traction_correction_out->SetSize(2 * num_fault_dofs_);
       *traction_correction_out = 0.0;
    }
    if (jump_residual_out)
    {
       jump_residual_out->SetSize(2 * num_fault_dofs_);
       *jump_residual_out = 0.0;
    }
   if (normal_stress_out)
   {
      normal_stress_out->SetSize(num_fault_dofs_);
      *normal_stress_out = 0.0;
   }
   if (normal_correction_out)
   {
      normal_correction_out->SetSize(num_fault_dofs_);
      *normal_correction_out = 0.0;
   }

   // v51: Elastic normal traction for sigma_n feedback
   if (normal_traction)
   {
      normal_traction->SetSize(num_fault_dofs_);
      *normal_traction = 0.0;
   }

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      int face = fault_interior_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      // ============================================================
      // v35: Per-quadrature-point traction evaluation with face averaging
      //
      // Previous code evaluated {σ·n̂} at face centroid only and used
      // face-averaged shapes for BR2 correction. This is exact for p=1
      // but wrong for p≥2 (stress varies across face, lifted function
      // has spatial variation). Now we evaluate everything at each face
      // quadrature point and face-average the result.
      // ============================================================

      // Face centroid IP (used for diagnostic coordinate output only)
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());

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

      // Fault basis
      const auto &basis = fault_basis_.GetBasis(fi);

      // Element Jacobian inverses (constant for linear tets/hexes)
      // Must set an integration point first so Jacobian() is valid.
      FTr->SetAllIntPoints(&ip);
      DenseMatrix Jinv1(dim), Jinv2(dim);
      CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
      CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

      // Face quadrature rule
      int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
      int quad_order_trac = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
      const IntegrationRule &ir_trac = IntRules.Get(
         FTr->GetGeometryType(), quad_order_trac);
      int nqp = ir_trac.GetNPoints();

      // Face-averaged accumulators (used for BR2 path and diagnostics)
      real_t T_global[3] = {0.0, 0.0, 0.0};
      real_t T_stress[3] = {0.0, 0.0, 0.0};
      real_t correction[3] = {0.0, 0.0, 0.0};
      real_t residual_avg[3] = {0.0, 0.0, 0.0};
      real_t sum_wq = 0.0;

      if (method_ == DGMethod::IP)
      {
         // Tandem-style traction recovery using combined integrator.
         // ComputeTractionAtQuadPoints (per-qp geometry, same penalty as K)
         // + ProjectTractionToFaultDOFs (nl_q-weighted L2 projection)
         // When decomposition diagnostics are requested, uses the decomposed
         // variant to also output stress, correction, and jump residual.
         int nbf = nbf_per_face_;
         bool need_decomp = traction_stress_out || traction_correction_out ||
                            jump_residual_out ||
                            normal_stress_out || normal_correction_out;

         // Build sign-corrected slip at quad points (Tandem evaluate_slip).
         // 1. Collect tangential slip components (dip, strike) per DOF
         // 2. Interpolate to quad points using face basis functions e_q
         // 3. Embed at each QP using per-QP tangent frame
         // This matches Tandem's tensor contraction:
         //   slip_q[p,q] = e_q[l,q] * fault_basis_q[p,o,q] * slip[l,n] * copy_slip[n,o]
         int qo_slip = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
         const IntegrationRule &ir_slip = IntRules.Get(
            FTr->GetGeometryType(), qo_slip);
         int nqp_slip = ir_slip.GetNPoints();
         Vector delta_u_quad_t;

         if (!basis.qp_data.empty())
         {
            // Per-QP tangent embedding (matches Tandem for any nbf)
            // Step 1: tangential slip components per DOF
            Vector slip_tang(2 * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               slip_tang(0 * nbf + kk) = slip_bc(2 * dof_idx);      // dip
               slip_tang(1 * nbf + kk) = slip_bc(2 * dof_idx + 1);  // strike
            }
            // Step 2: interpolate tangential components to QPs
            Vector slip_tang_q;
            face_quad_->InterpolateToQuadPoints(2, slip_tang, slip_tang_q);
            // Step 3: embed at each QP using per-QP tangent frame
            delta_u_quad_t.SetSize(dim * nqp_slip);
            for (int q = 0; q < nqp_slip; q++)
            {
               real_t sl_q[2] = {slip_tang_q(q), slip_tang_q(nqp_slip + q)};
               real_t du[3];
               fault_basis_.EmbedSlipQP(fi, q, sl_q, du);
               for (int c = 0; c < dim; c++)
                  delta_u_quad_t(c * nqp_slip + q) = du[c];
            }
         }
         else
         {
            // Fallback: centroid tangents (exact for flat faces, used by BR2)
            Vector delta_u_nodal_t(dim * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               real_t sl[2] = {slip_bc(2 * dof_idx), slip_bc(2 * dof_idx + 1)};
               real_t du[3];
               fault_basis_.EmbedSlip(fi, sl, du);
               for (int c = 0; c < dim; c++)
                  delta_u_nodal_t(c * nbf + kk) = du[c];
            }
            face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal_t, delta_u_quad_t);
         }

         // Step 1: Traction at quad points (Tandem compute_traction)
         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad_new, nl_q_vec;
         Vector T_stress_quad_dec, T_corr_quad_dec, R_quad_dec;
         if (need_decomp)
         {
            trac_integ.ComputeTractionAtQuadPointsDecomposed(
               *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_t,
               T_quad_new,
               traction_stress_out ? &T_stress_quad_dec : nullptr,
               traction_correction_out ? &T_corr_quad_dec : nullptr,
               jump_residual_out ? &R_quad_dec : nullptr,
               nullptr, &nl_q_vec);
         }
         else
         {
            trac_integ.ComputeTractionAtQuadPoints(
               *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_t,
               T_quad_new, nullptr, &nl_q_vec);
         }
         int nqp_new = T_quad_new.Size() / dim;

         // Step 2: Project to fault DOFs (Tandem evaluate_traction)
         int quad_order_new = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
         const IntegrationRule &ir_new = IntRules.Get(
            FTr->GetGeometryType(), quad_order_new);
         const DenseMatrix &e_q = face_quad_->BasisAtQuadPoints();

         const auto *qpd = basis.qp_data.empty() ? nullptr : &basis.qp_data;

         // Unified projection: project all local components through one
         // L2 solve (matching Tandem's single-kernel approach).
         // When normal_traction is requested, project 3 components
         // (normal, dip, strike); otherwise project 2 (dip, strike).
         int ncomp_proj = normal_traction ? 3 : 2;
         real_t basis_vecs[3][3];
         if (normal_traction)
         {
            for (int d = 0; d < 3; d++)
            {
               basis_vecs[0][d] = basis.normal[d];
               basis_vecs[1][d] = basis.tangent1[d];
               basis_vecs[2][d] = basis.tangent2[d];
            }
         }
         else
         {
            for (int d = 0; d < 3; d++)
            {
               basis_vecs[0][d] = basis.tangent1[d];
               basis_vecs[1][d] = basis.tangent2[d];
               basis_vecs[2][d] = 0.0;
            }
         }

         Vector trac_local_new;
         DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
            dim, ncomp_proj, T_quad_new, nl_q_vec, ir_new, nbf, e_q,
            basis_vecs, trac_local_new, qpd);

         // Store per-DOF shear traction
         int tang_offset = normal_traction ? 1 : 0;
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = fi * nbf_per_face_ + kk;
            traction(2 * dof_idx)     = trac_local_new((tang_offset + 0) * nbf + kk);
            traction(2 * dof_idx + 1) = trac_local_new((tang_offset + 1) * nbf + kk);
         }

         // Normal stress: extracted from the unified 3-component projection
         // above (Tandem-style: normal + dip + strike through one L2 solve).
         // Sign convention: positive in compression = -(T · n̂).
         if (normal_traction)
         {
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_traction)(dof_idx) = -trac_local_new(0 * nbf + kk);
            }
         }

         // Project decomposition components to fault DOFs (if requested)
         if (traction_stress_out)
         {
            Vector stress_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, T_stress_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, stress_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*traction_stress_out)(2 * dof_idx)     = stress_local(0 * nbf + kk);
               (*traction_stress_out)(2 * dof_idx + 1) = stress_local(1 * nbf + kk);
            }
         }
         if (traction_correction_out)
         {
            Vector corr_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, T_corr_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, corr_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*traction_correction_out)(2 * dof_idx)     = corr_local(0 * nbf + kk);
               (*traction_correction_out)(2 * dof_idx + 1) = corr_local(1 * nbf + kk);
            }
         }
         if (jump_residual_out)
         {
            Vector res_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, R_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, res_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*jump_residual_out)(2 * dof_idx)     = res_local(0 * nbf + kk);
               (*jump_residual_out)(2 * dof_idx + 1) = res_local(1 * nbf + kk);
            }
         }

         // Normal decomposition: project stress and correction onto normal
         // Same sign convention as total normal traction: -(T · n̂)
         if (normal_stress_out && normal_traction)
         {
            Vector ns_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 1, T_stress_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs, ns_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_stress_out)(dof_idx) = -ns_local(kk);
            }
         }
         if (normal_correction_out && normal_traction)
         {
            Vector nc_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 1, T_corr_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs, nc_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_correction_out)(dof_idx) = -nc_local(kk);
            }
         }

      }
      else  // BR2
      {
         // BR2 always uses nbf=1, so face-index == DOF-index
         real_t slip_local_br2[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(fi, slip_local_br2, delta_u);

         if (!mass_inv_computed_) {
            PrecomputeMassInverse();
            // Re-fetch FTr: PrecomputeMassInverse calls GetElementTransformation
            // for every element, which may invalidate MFEM's cached face transforms.
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }

         Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
         real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                  ? real_t(dim + 1) : real_t(2 * dim);

         const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
         const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

         // Precompute shapes, normals, and face integrals for BR2 lifting
         DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
         Vector nor_arr(dim * nqp), w_arr(nqp);

         DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
         face_int1 = 0.0;
         face_int2 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &fip = ir_trac.IntPoint(q);
            FTr->SetAllIntPoints(&fip);
            const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();

            Vector s1q(shapes1.GetColumn(q), ndof1);
            Vector s2q(shapes2.GetColumn(q), ndof2);
            fe1->CalcShape(eip1_q, s1q);
            fe2->CalcShape(eip2_q, s2q);

            Vector nor_q(&nor_arr[q * dim], dim);
            CalcOrtho(FTr->Jacobian(), nor_q);
            w_arr[q] = fip.weight;

            // Per-quadrature-point displacement jump for lifting
            real_t jump_q[3] = {0.0, 0.0, 0.0};
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++)
                  u1q += s1q(k) * u1_all(c * ndof1 + k);
               for (int k = 0; k < ndof2; k++)
                  u2q += s2q(k) * u2_all(c * ndof2 + k);
               jump_q[c] = (u1q - u2q) - delta_u[c];
            }

            for (int u = 0; u < dim; u++)
               for (int s = 0; s < dim; s++)
               {
                  for (int m = 0; m < ndof1; m++)
                     face_int1(u * dim + s, m) +=
                        w_arr[q] * s1q(m) * jump_q[u] * nor_q(s);
                  for (int m = 0; m < ndof2; m++)
                     face_int2(u * dim + s, m) +=
                        w_arr[q] * s2q(m) * jump_q[u] * nor_q(s);
               }
         }

         // Lift: f_lifted = 0.5 * Minv * face_int
         DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
         MultABt(face_int1, Minv1, f_lifted1);
         f_lifted1 *= 0.5;
         MultABt(face_int2, Minv2, f_lifted2);
         f_lifted2 *= 0.5;

         // Pass 2: evaluate traction at each quadrature point
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &fip = ir_trac.IntPoint(q);
            FTr->SetAllIntPoints(&fip);
            const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();
            real_t wq = w_arr[q];

            // {σ·n̂} at quadrature point q
            DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
            fe1->CalcDShape(eip1_q, dshape1_ref);
            fe2->CalcDShape(eip2_q, dshape2_ref);
            DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
            Mult(dshape1_ref, Jinv1, dshape1_phys);
            Mult(dshape2_ref, Jinv2, dshape2_phys);

            DenseMatrix grad1(dim, dim), grad2(dim, dim);
            grad1 = 0.0; grad2 = 0.0;
            for (int c = 0; c < dim; c++)
               for (int d = 0; d < dim; d++)
               {
                  for (int k = 0; k < ndof1; k++)
                     grad1(c, d) += dshape1_phys(k, d)
                                    * u1_all(c * ndof1 + k);
                  for (int k = 0; k < ndof2; k++)
                     grad2(c, d) += dshape2_phys(k, d)
                                    * u2_all(c * ndof2 + k);
               }

            real_t T_stress_q[3] = {0.0, 0.0, 0.0};
            for (int ci = 0; ci < dim; ci++)
            {
               for (int cj = 0; cj < dim; cj++)
               {
                  real_t ag = 0.5 * (grad1(ci, cj) + grad2(ci, cj));
                  real_t ag_t = 0.5 * (grad1(cj, ci) + grad2(cj, ci));
                  real_t eps_ij = 0.5 * (ag + ag_t);
                  real_t tr_contrib = (ci == cj)
                     ? lambda_val_ * (0.5*((grad1(0,0)+grad2(0,0))
                        + (grad1(1,1)+grad2(1,1))
                        + (grad1(2,2)+grad2(2,2)))) : 0.0;
                  real_t stress_ij = tr_contrib + 2.0 * mu_val_ * eps_ij;
                  T_stress_q[ci] += stress_ij * basis.normal[cj];
               }
            }

            // BR2 correction at quadrature point q:
            // Evaluate lifted function using per-point shapes (not avg_shapes)
            real_t correction_q[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < dim; i++)
            {
               real_t sum = 0.0;
               for (int u = 0; u < dim; u++)
                  for (int s = 0; s < dim; s++)
                  {
                     real_t tn = lambda_val_
                           * (u == s ? 1.0 : 0.0) * basis.normal[i]
                        + mu_val_
                           * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                              + (i == s ? 1.0 : 0.0) * basis.normal[u]);
                     real_t eval1 = 0.0, eval2 = 0.0;
                     for (int m = 0; m < ndof1; m++)
                        eval1 += shapes1(m, q)
                                 * f_lifted1(u * dim + s, m);
                     for (int m = 0; m < ndof2; m++)
                        eval2 += shapes2(m, q)
                                 * f_lifted2(u * dim + s, m);
                     sum += tn * (eval1 + eval2);
                  }
               correction_q[i] = br2_penalty * 0.5 * sum;
            }

            // Accumulate face average
            for (int c = 0; c < dim; c++)
            {
               T_stress[c] += wq * T_stress_q[c];
               correction[c] += wq * correction_q[c];
            }
            sum_wq += wq;
         }
      }

      // BR2 path: normalize face average and store (BR2 uses nbf=1 always)
      if (method_ != DGMethod::IP)
      {
         if (sum_wq > 0.0)
         {
            for (int c = 0; c < dim; c++)
            {
               T_stress[c] /= sum_wq;
               correction[c] /= sum_wq;
               T_global[c] = T_stress[c] - correction[c];
            }
         }

         // Project to local frame: (tau_dip, tau_strike)
         // BR2 always has nbf_per_face_=1, so dof_idx = fi
         real_t tau_local[2];
         fault_basis_.ProjectTraction(fi, T_global, tau_local);
         traction(2 * fi)     = tau_local[0];
         traction(2 * fi + 1) = tau_local[1];
         if (traction_stress_out || traction_correction_out)
         {
            real_t tau_stress_local[2], tau_corr_local[2];
            fault_basis_.ProjectTraction(fi, T_stress, tau_stress_local);
            real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
            fault_basis_.ProjectTraction(fi, corr_neg, tau_corr_local);
            if (traction_stress_out)
            {
               (*traction_stress_out)(2 * fi) = tau_stress_local[0];
               (*traction_stress_out)(2 * fi + 1) = tau_stress_local[1];
            }
            if (traction_correction_out)
            {
               (*traction_correction_out)(2 * fi) = tau_corr_local[0];
               (*traction_correction_out)(2 * fi + 1) = tau_corr_local[1];
            }
         }
         if (jump_residual_out)
         {
            real_t res_local[2];
            fault_basis_.ProjectTraction(fi, residual_avg, res_local);
            (*jump_residual_out)(2 * fi) = res_local[0];
            (*jump_residual_out)(2 * fi + 1) = res_local[1];
         }
         // v51: elastic normal traction for sigma_n feedback
         if (normal_traction)
         {
            (*normal_traction)(fi) = fault_basis_.NormalStress(fi, T_global);
         }
      }
   }

   // Shared fault faces (parallel only)
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      int rank;
      MPI_Comm_rank(mesh_.GetComm(), &rank);
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

         // v35: Per-quadrature-point traction (same approach as interior faces)
         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());

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

         Array<int> vdofs2;
         pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
         const Vector &nbr_data = par_u.FaceNbrData();
         Vector u2_all(vdofs2.Size());
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            u2_all(j) = nbr_data(vdofs2[j]);
         }

         // Fault basis
         const auto &basis = fault_basis_.GetBasis(trac_idx);

         // Element Jacobian inverses (constant for linear tets)
         DenseMatrix Jinv1(dim), Jinv2(dim);
         CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
         CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

         // Face quadrature rule
         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         int quad_order_trac_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir_trac = IntRules.Get(
            FTr->GetGeometryType(), quad_order_trac_sh);
         int nqp = ir_trac.GetNPoints();

         // Face-averaged accumulators (used for BR2 path and diagnostics)
         real_t T_global[3] = {0.0, 0.0, 0.0};
         real_t T_stress[3] = {0.0, 0.0, 0.0};
         real_t correction[3] = {0.0, 0.0, 0.0};
         real_t residual_avg[3] = {0.0, 0.0, 0.0};
         real_t sum_wq = 0.0;


         if (method_ == DGMethod::IP)
         {
            // Tandem-style traction for shared faces (same as interior).
            // When decomposition is requested, uses decomposed variant.
            int nbf_sh = nbf_per_face_;
            bool need_decomp_sh = traction_stress_out ||
                                  traction_correction_out || jump_residual_out ||
                                  normal_stress_out || normal_correction_out;
            const auto &basis_sh = fault_basis_.GetBasis(trac_idx);

            // Build sign-corrected slip at quad points (shared faces).
            // Same Tandem-style per-QP tangent embedding as interior faces.
            int quad_order_sh_slip = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_sh_slip = IntRules.Get(
               FTr->GetGeometryType(), quad_order_sh_slip);
            int nqp_sh_slip = ir_sh_slip.GetNPoints();
            Vector delta_u_quad_sh;

            if (!basis_sh.qp_data.empty())
            {
               // Per-QP tangent embedding (Tandem evaluate_slip)
               Vector slip_tang_sh(2 * nbf_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_sh + kk;
                  slip_tang_sh(0 * nbf_sh + kk) = slip_bc(2 * dof_idx);      // dip
                  slip_tang_sh(1 * nbf_sh + kk) = slip_bc(2 * dof_idx + 1);  // strike
               }
               Vector slip_tang_q_sh;
               face_quad_->InterpolateToQuadPoints(2, slip_tang_sh, slip_tang_q_sh);
               delta_u_quad_sh.SetSize(dim * nqp_sh_slip);
               for (int q = 0; q < nqp_sh_slip; q++)
               {
                  real_t sl_q[2] = {slip_tang_q_sh(q),
                                    slip_tang_q_sh(nqp_sh_slip + q)};
                  real_t du[3];
                  fault_basis_.EmbedSlipQP(trac_idx, q, sl_q, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_quad_sh(c * nqp_sh_slip + q) = du[c];
               }
            }
            else
            {
               // Fallback: centroid tangents (BR2 path)
               Vector delta_u_nodal_sh(dim * nbf_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_sh + kk;
                  real_t sl[2] = {slip_bc(2 * dof_idx), slip_bc(2 * dof_idx + 1)};
                  real_t du[3];
                  fault_basis_.EmbedSlip(trac_idx, sl, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_nodal_sh(c * nbf_sh + kk) = du[c];
               }
               face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal_sh, delta_u_quad_sh);
            }

            DGElasticityIPCombinedIntegrator trac_integ_sh(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
            Vector T_quad_sh, nl_q_sh;
            Vector T_stress_quad_sh, T_corr_quad_sh, R_quad_sh;
            if (need_decomp_sh)
            {
               trac_integ_sh.ComputeTractionAtQuadPointsDecomposed(
                  *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_sh,
                  T_quad_sh,
                  traction_stress_out ? &T_stress_quad_sh : nullptr,
                  traction_correction_out ? &T_corr_quad_sh : nullptr,
                  jump_residual_out ? &R_quad_sh : nullptr,
                  nullptr, &nl_q_sh);
            }
            else
            {
               trac_integ_sh.ComputeTractionAtQuadPoints(
                  *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_sh,
                  T_quad_sh, nullptr, &nl_q_sh);
            }
            int nqp_sh = T_quad_sh.Size() / dim;

            int quad_order_sh2 = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_sh2 = IntRules.Get(
               FTr->GetGeometryType(), quad_order_sh2);
            const DenseMatrix &e_q_sh = face_quad_->BasisAtQuadPoints();

            const auto *qpd_sh = basis_sh.qp_data.empty() ? nullptr : &basis_sh.qp_data;

            // Unified projection for shared faces (same as interior)
            int ncomp_sh = normal_traction ? 3 : 2;
            real_t basis_vecs_sh[3][3];
            if (normal_traction)
            {
               for (int d = 0; d < 3; d++)
               {
                  basis_vecs_sh[0][d] = basis_sh.normal[d];
                  basis_vecs_sh[1][d] = basis_sh.tangent1[d];
                  basis_vecs_sh[2][d] = basis_sh.tangent2[d];
               }
            }
            else
            {
               for (int d = 0; d < 3; d++)
               {
                  basis_vecs_sh[0][d] = basis_sh.tangent1[d];
                  basis_vecs_sh[1][d] = basis_sh.tangent2[d];
                  basis_vecs_sh[2][d] = 0.0;
               }
            }

            Vector trac_local_sh;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, ncomp_sh, T_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
               basis_vecs_sh, trac_local_sh, qpd_sh);

            int tang_off_sh = normal_traction ? 1 : 0;
            for (int kk = 0; kk < nbf_sh; kk++)
            {
               int dof_idx = trac_idx * nbf_per_face_ + kk;
               traction(2 * dof_idx)     = trac_local_sh((tang_off_sh + 0) * nbf_sh + kk);
               traction(2 * dof_idx + 1) = trac_local_sh((tang_off_sh + 1) * nbf_sh + kk);
            }
            if (normal_traction)
            {
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_traction)(dof_idx) = -trac_local_sh(0 * nbf_sh + kk);
               }
            }

            // Project decomposition to fault DOFs (shared faces)
            if (traction_stress_out)
            {
               Vector stress_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, T_stress_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, stress_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*traction_stress_out)(2 * dof_idx)     = stress_local_sh(0 * nbf_sh + kk);
                  (*traction_stress_out)(2 * dof_idx + 1) = stress_local_sh(1 * nbf_sh + kk);
               }
            }
            if (traction_correction_out)
            {
               Vector corr_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, T_corr_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, corr_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*traction_correction_out)(2 * dof_idx)     = corr_local_sh(0 * nbf_sh + kk);
                  (*traction_correction_out)(2 * dof_idx + 1) = corr_local_sh(1 * nbf_sh + kk);
               }
            }
            if (jump_residual_out)
            {
               Vector res_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, R_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, res_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*jump_residual_out)(2 * dof_idx)     = res_local_sh(0 * nbf_sh + kk);
                  (*jump_residual_out)(2 * dof_idx + 1) = res_local_sh(1 * nbf_sh + kk);
               }
            }

            // Normal decomposition for shared faces
            if (normal_stress_out && normal_traction)
            {
               Vector ns_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 1, T_stress_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh, ns_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_stress_out)(dof_idx) = -ns_local_sh(kk);
               }
            }
            if (normal_correction_out && normal_traction)
            {
               Vector nc_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 1, T_corr_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh, nc_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_correction_out)(dof_idx) = -nc_local_sh(kk);
               }
            }

         }
         else  // BR2
         {
            // BR2 always uses nbf=1, so face-index == DOF-index
            real_t slip_local_br2[2] = {slip_bc(2 * trac_idx),
                                        slip_bc(2 * trac_idx + 1)};
            real_t delta_u[3];
            fault_basis_.EmbedSlip(trac_idx, slip_local_br2, delta_u);

            if (!mass_inv_computed_) {
               PrecomputeMassInverse();
               // Re-fetch FTr: PrecomputeMassInverse invalidates cached transforms.
               FTr = mesh_.GetSharedFaceTransformations(sf);
            }

            Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                     ? real_t(dim + 1) : real_t(2 * dim);

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &fip = ir_trac.IntPoint(q);
               FTr->SetAllIntPoints(&fip);
               const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();

               Vector s1q(shapes1.GetColumn(q), ndof1);
               Vector s2q(shapes2.GetColumn(q), ndof2);
               fe1->CalcShape(eip1_q, s1q);
               fe2->CalcShape(eip2_q, s2q);

               Vector nor_q(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nor_q);
               w_arr[q] = fip.weight;

               real_t jump_q[3] = {0.0, 0.0, 0.0};
               for (int c = 0; c < dim; c++)
               {
                  real_t u1q = 0.0, u2q = 0.0;
                  for (int k = 0; k < ndof1; k++)
                     u1q += s1q(k) * u1_all(c * ndof1 + k);
                  for (int k = 0; k < ndof2; k++)
                     u2q += s2q(k) * u2_all(c * ndof2 + k);
                  jump_q[c] = (u1q - u2q) - delta_u[c];
               }

               for (int u = 0; u < dim; u++)
                  for (int s = 0; s < dim; s++)
                  {
                     for (int m = 0; m < ndof1; m++)
                        face_int1(u * dim + s, m) +=
                           w_arr[q] * s1q(m) * jump_q[u] * nor_q(s);
                     for (int m = 0; m < ndof2; m++)
                        face_int2(u * dim + s, m) +=
                           w_arr[q] * s2q(m) * jump_q[u] * nor_q(s);
                  }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            f_lifted1 *= 0.5;
            MultABt(face_int2, Minv2, f_lifted2);
            f_lifted2 *= 0.5;

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &fip = ir_trac.IntPoint(q);
               FTr->SetAllIntPoints(&fip);
               const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();
               real_t wq = w_arr[q];

               DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
               fe1->CalcDShape(eip1_q, dshape1_ref);
               fe2->CalcDShape(eip2_q, dshape2_ref);
               DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
               Mult(dshape1_ref, Jinv1, dshape1_phys);
               Mult(dshape2_ref, Jinv2, dshape2_phys);

               DenseMatrix grad1(dim, dim), grad2(dim, dim);
               grad1 = 0.0; grad2 = 0.0;
               for (int c = 0; c < dim; c++)
                  for (int d = 0; d < dim; d++)
                  {
                     for (int k = 0; k < ndof1; k++)
                        grad1(c, d) += dshape1_phys(k, d)
                                       * u1_all(c * ndof1 + k);
                     for (int k = 0; k < ndof2; k++)
                        grad2(c, d) += dshape2_phys(k, d)
                                       * u2_all(c * ndof2 + k);
                  }

               real_t T_stress_q[3] = {0.0, 0.0, 0.0};
               for (int ci = 0; ci < dim; ci++)
                  for (int cj = 0; cj < dim; cj++)
                  {
                     real_t ag = 0.5 * (grad1(ci, cj) + grad2(ci, cj));
                     real_t ag_t = 0.5 * (grad1(cj, ci) + grad2(cj, ci));
                     real_t eps_ij = 0.5 * (ag + ag_t);
                     real_t tr_contrib = (ci == cj)
                        ? lambda_val_ * (0.5*((grad1(0,0)+grad2(0,0))
                           + (grad1(1,1)+grad2(1,1))
                           + (grad1(2,2)+grad2(2,2)))) : 0.0;
                     real_t stress_ij = tr_contrib + 2.0 * mu_val_ * eps_ij;
                     T_stress_q[ci] += stress_ij * basis.normal[cj];
                  }

               real_t correction_q[3] = {0.0, 0.0, 0.0};
               for (int ci = 0; ci < dim; ci++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_
                              * (u == s ? 1.0 : 0.0) * basis.normal[ci]
                           + mu_val_
                              * ((ci == u ? 1.0 : 0.0) * basis.normal[s]
                                 + (ci == s ? 1.0 : 0.0) * basis.normal[u]);
                        real_t eval1 = 0.0, eval2 = 0.0;
                        for (int m = 0; m < ndof1; m++)
                           eval1 += shapes1(m, q)
                                    * f_lifted1(u * dim + s, m);
                        for (int m = 0; m < ndof2; m++)
                           eval2 += shapes2(m, q)
                                    * f_lifted2(u * dim + s, m);
                        sum += tn * (eval1 + eval2);
                     }
                  correction_q[ci] = br2_penalty * 0.5 * sum;
               }

               for (int c = 0; c < dim; c++)
               {
                  T_stress[c] += wq * T_stress_q[c];
                  correction[c] += wq * correction_q[c];
               }
               sum_wq += wq;
            }
         }

         // BR2 path: normalize face average and store (BR2 uses nbf=1 always)
         if (method_ != DGMethod::IP)
         {
            if (sum_wq > 0.0)
            {
               for (int c = 0; c < dim; c++)
               {
                  T_stress[c] /= sum_wq;
                  correction[c] /= sum_wq;
                  residual_avg[c] /= sum_wq;
                  T_global[c] = T_stress[c] - correction[c];
               }
            }

            // BR2 always has nbf_per_face_=1, so dof_idx = trac_idx
            real_t tau_local[2];
            fault_basis_.ProjectTraction(trac_idx, T_global, tau_local);
            traction(2 * trac_idx)     = tau_local[0];
            traction(2 * trac_idx + 1) = tau_local[1];
            if (traction_stress_out || traction_correction_out)
            {
               real_t tau_stress_local[2], tau_corr_local[2];
               fault_basis_.ProjectTraction(trac_idx, T_stress, tau_stress_local);
               real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
               fault_basis_.ProjectTraction(trac_idx, corr_neg, tau_corr_local);
               if (traction_stress_out)
               {
                  (*traction_stress_out)(2 * trac_idx) = tau_stress_local[0];
                  (*traction_stress_out)(2 * trac_idx + 1) = tau_stress_local[1];
               }
               if (traction_correction_out)
               {
                  (*traction_correction_out)(2 * trac_idx) = tau_corr_local[0];
                  (*traction_correction_out)(2 * trac_idx + 1) = tau_corr_local[1];
               }
            }
            if (jump_residual_out)
            {
               real_t res_local[2];
               fault_basis_.ProjectTraction(trac_idx, residual_avg, res_local);
               (*jump_residual_out)(2 * trac_idx) = res_local[0];
               (*jump_residual_out)(2 * trac_idx + 1) = res_local[1];
            }
            // v51: elastic normal traction for sigma_n feedback
            if (normal_traction)
            {
               (*normal_traction)(trac_idx) =
                  fault_basis_.NormalStress(trac_idx, T_global);
            }
         }
      }
#endif
   }

   // Diagnostic: check for traction blowup
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      int rank;
      MPI_Comm_rank(mesh_.GetComm(), &rank);
      const int num_interior_faces = fault_interior_faces_.Size();
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         real_t tau_mag = std::sqrt(traction(2*i)*traction(2*i) +
                                    traction(2*i+1)*traction(2*i+1));
         if (tau_mag > 1e9 || std::isnan(tau_mag))
         {
            // Convert DOF index to face index (nbf_per_face_ DOFs per face)
            int face_i = i / nbf_per_face_;
            // Get face coordinates for diagnostics
            Vector face_center(3);
            face_center = 0.0;
            if (face_i < num_interior_faces)
            {
               int face_idx = fault_interior_faces_[face_i];
               FaceElementTransformations *FTr =
                  mesh_.GetInteriorFaceTransformations(face_idx);
               if (FTr)
               {
                  const IntegrationPoint &ip =
                     Geometries.GetCenter(FTr->GetGeometryType());
                  FTr->Face->SetIntPoint(&ip);
                  FTr->Face->Transform(ip, face_center);
               }
            }
            else
            {
               int shared_idx = face_i - num_interior_faces;
               int sf = fault_shared_faces_[shared_idx];
               FaceElementTransformations *FTr =
                  mesh_.GetSharedFaceTransformations(sf);
               if (FTr)
               {
                  const IntegrationPoint &ip =
                     Geometries.GetCenter(FTr->GetGeometryType());
                  FTr->Face->SetIntPoint(&ip);
                  FTr->Face->Transform(ip, face_center);
               }
            }
            mfem::out << "[Rank " << rank << "] TRACTION BLOWUP: DOF " << i
                      << " (face " << face_i << "/"
                      << (face_i < num_interior_faces ? "interior" : "shared")
                      << ")"
                      << " tau_mag=" << tau_mag
                      << " tau=(" << traction(2*i) << ","
                      << traction(2*i+1) << ")"
                      << " at x=(" << face_center(0) << ","
                      << face_center(1) << "," << face_center(2) << ")"
                      << " slip=(" << slip_bc(2*i) << ","
                      << slip_bc(2*i+1) << ")\n";
         }
      }
   }
#endif

}

// Convenience type alias
using SerialElasticityDomainOperator = ElasticityDomainOperator<Mesh>;

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ELASTICITY_OPERATOR_HPP
