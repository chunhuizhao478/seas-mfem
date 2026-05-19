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
#include "../common/fault_scatter.hpp"
#include "../constitutive/constitutive_model.hpp"
#include "../constitutive/linear_elastic.hpp"
#include "boundary_config.hpp"
#include "domain_config.hpp"
#include "../fault/fault_basis.hpp"
#include "../integrator/dg_elasticity_br2_integrator.hpp"
#include "../integrator/dg_elasticity_ip_penalty_integrator.hpp"
#include "../integrator/dg_elasticity_ip_combined_integrator.hpp"
#include "../fault/face_quadrature.hpp"
#include "matrix_stats.hpp"

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

   /// @brief Construct with explicit ConstitutiveModel + BoundaryConfig (Phase 4+ path).
   ///
   /// The model must outlive this operator (non-owning pointer stored).
   /// BoundaryConfig is required — no default. User specifies boundary
   /// attributes explicitly from their mesh file.
   ElasticityDomainOperator(MeshType &mesh, int order,
                             const ConstitutiveModel &model,
                             real_t Vp, real_t Wf, real_t lf,
                             const BoundaryConfig &bdr_config,
                             DGMethod method = DGMethod::BR2,
                             SolverType solver_type = SolverType::MUMPS_BLR,
                             const DomainConfig &config = {})
      : mesh_(mesh), order_(order),
        model_(&model),
        bdr_config_(bdr_config),
        Vp_(Vp), Wf_(Wf), lf_(lf),
        method_(method), solver_type_(solver_type),
        bc_mode_(BCMode::FarField),
        check_residual_(config.check_residual),
        mass_inv_computed_(false),
        fault_depths_computed_(false),
        fault_coords_computed_(false),
        face_basis_type_(config.face_basis_type),
        penalty_factor_(config.penalty_factor),
        blr_tol_(config.blr_tol),
        match_quad_order_(config.match_quad_order)
   {
      const auto *le = dynamic_cast<const LinearElastic *>(model_);
      MFEM_VERIFY(le, "ElasticityDomainOperator currently requires LinearElastic");
      lambda_val_ = le->GetLambda();
      mu_val_ = le->GetMu();
      lambda_coeff_.constant = lambda_val_;
      mu_coeff_.constant = mu_val_;

      InitOperator();
   }

   /// @brief Legacy constructor: scalar lambda/mu (deprecated, delegates).
   ///
   /// Creates an internal LinearElastic model and builds BoundaryConfig
   /// internally from BCMode with legacy BP5 attribute numbers.
   /// Existing tests and drivers continue to work unchanged.
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
      owned_model_ = std::make_unique<LinearElastic>(lambda, mu);
      model_ = owned_model_.get();

      // Build BoundaryConfig from legacy BCMode numbers.
      bdr_config_.fault_attr = 3;
      bdr_config_.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
      if (bc_mode == BCMode::AllDirichlet)
      {
         mfem::out << "\n  *** WARNING: AllDirichlet BC mode is legacy and known "
                   << "to be incorrect for BP5. Use BCMode::FarField. ***\n\n";
         MFEM_VERIFY(mesh_.bdr_attributes.Size() > 0,
                     "BCMode::AllDirichlet requires boundary elements, "
                     "but mesh has none");
         int max_attr = mesh_.bdr_attributes.Max();
         for (int a = 1; a <= max_attr; a++)
         {
            bdr_config_.dirichlet_attrs.insert(a);
         }
      }
      else if (bc_mode == BCMode::XOnly)
      {
         bdr_config_.dirichlet_attrs = {1, 2};
         bdr_config_.natural_attrs = {3, 4, 5, 6};
      }
      else // FarField (default, correct for BP5)
      {
         bdr_config_.dirichlet_attrs = {5};
         bdr_config_.natural_attrs = {1};
      }

      InitOperator();
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
                                   Vector *normal_correction = nullptr) override;

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

   void GetFaultDOFCoords3D(Vector &dof_coords_3d) const override;

   void GetFaultDOFBasis(DenseMatrix &dof_basis) const override;

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

   /// Access the constitutive model.
   const ConstitutiveModel &GetModel() const { return *model_; }

   // Legacy setters (deprecated — use DomainConfig in the new constructor).
   // Kept for backward compatibility with existing drivers/tests.
   void SetCheckResidual(bool check) { check_residual_ = check; }
   void SetBLRTol(real_t tol) { blr_tol_ = tol; }
   void SetMatchQuadOrder(bool v) { match_quad_order_ = v; }
   void SetPenaltyFactor(real_t f) { penalty_factor_ = f; }
   void SetFaceBasisType(int bt) { face_basis_type_ = bt; }

   /// Enable one-shot matrix analysis right after ParallelAssemble.
   /// Computes nnz, symmetry, diagonal sanity, etc., and writes a JSON
   /// summary to @a json_path (rank 0).  Implies MUMPS print-level 2
   /// when a MUMPS solver is selected (so fill-in stats are logged too).
   void SetMatrixStatsConfig(bool enable, const std::string &json_path)
   {
      matrix_stats_enabled_ = enable;
      matrix_stats_path_    = json_path;
   }
   bool MatrixStatsEnabled() const { return matrix_stats_enabled_; }

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

   bool IsFirstStepDebugEnabled() const override { return first_step_debug_.enabled; }


#include "elasticity_operator_verify.inl"


private:
   /// Common initialization shared by both constructors.
   /// Evaluate the Dirichlet function for a given boundary attribute.
   /// Uses per-attr function if available, otherwise default.
   /// Aborts if no function is set (no silent fallback).
   void EvalDirichletFunc(int attr, const Vector &x, real_t t,
                          Vector &u_D) const
   {
      auto it = bdr_config_.dirichlet_funcs.find(attr);
      if (it != bdr_config_.dirichlet_funcs.end())
      {
         it->second(x, t, u_D);
         return;
      }
      if (bdr_config_.default_dirichlet_func)
      {
         bdr_config_.default_dirichlet_func(x, t, u_D);
         return;
      }
      MFEM_ABORT("No DirichletFunc set for boundary attr " << attr
                 << ". Set default_dirichlet_func or per-attr func "
                    "in BoundaryConfig.");
   }

   /// Common initialization shared by both constructors.
   void InitOperator()
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

   MeshType &mesh_;
   int order_;

   // Constitutive model: model_ is always valid (non-owning pointer).
   // owned_model_ is set only when the legacy constructor creates it.
   const ConstitutiveModel *model_ = nullptr;
   std::unique_ptr<LinearElastic> owned_model_;

   // Cached scalar values from the model (for legacy code paths)
   real_t lambda_val_ = 0.0, mu_val_ = 0.0;

   // Boundary configuration (Phase 4+)
   BoundaryConfig bdr_config_;
   real_t Vp_, Wf_, lf_;
   DGMethod method_;
   SolverType solver_type_;
   BCMode bc_mode_;  // Legacy — kept for backward compatibility
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

   // One-shot matrix analysis (enabled via SetMatrixStatsConfig).
   bool matrix_stats_enabled_ = false;
   std::string matrix_stats_path_;

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


#include "elasticity_operator_debug.inl"


   // ---- Facet BC classification (single source of truth) ----
   // Mirrors Tandem's per-facet BC enum (DGOperatorTopo FacetInfo.bc).
   // Built once at startup from mesh tags, propagated exactly via MPI.
   enum class FacetBC : int8_t { None = 0, Fault = 1, Dirichlet = 2 };

   std::vector<FacetBC> face_bc_;         // [mesh_.GetNumFaces()] interior faces
   std::vector<int> face_bc_attr_;       // [mesh_.GetNumFaces()] boundary attr per face
   std::vector<FacetBC> shared_face_bc_;  // [mesh_.GetNSharedFaces()] shared faces
   std::vector<int> shared_face_bc_attr_; // [mesh_.GetNSharedFaces()] boundary attr

   // Legacy face arrays derived from the BC tables (kept for assembly loops)
   Array<int> fault_interior_faces_;
   Array<int> fault_shared_faces_;
   Array<int> dirichlet_interior_faces_;
   Array<int> dirichlet_shared_faces_;
   // Per-face boundary attr (parallel to dirichlet_interior/shared_faces_)
   Array<int> dirichlet_interior_attrs_;
   Array<int> dirichlet_shared_attrs_;

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

   // SharedFaultCommBlock defined in common/fault_scatter.hpp
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


#include "elasticity_operator_setup.inl"


#include "elasticity_operator_assembly.inl"

};


#include "elasticity_operator_traction.inl"


// Convenience type alias
using SerialElasticityDomainOperator = ElasticityDomainOperator<Mesh>;

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ELASTICITY_OPERATOR_HPP
