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
#include <iomanip>

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
                                   Vector *normal_traction = nullptr);

   /// Assemble only the fault-slip RHS contribution into the DG displacement
   /// space. This is a test/debug utility for checking K*u against b(slip)
   /// without Dirichlet loading.
   void AssembleSlipOnlyRHS(Vector &rhs, const Vector &slip_bc) const;

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

   /// Enable/disable traction decomposition diagnostic.
   /// When enabled, ComputeTraction prints T_stress, T_penalty, jump, and penalty
   /// for each fault DOF (helps identify stress vs penalty instability sources).
   void SetDiagTractionDecomp(bool enable) { diag_traction_decomp_ = enable; }

   /// Set MUMPS-BLR tolerance (default 1e-10). Lower = more accurate, more memory.
   /// Only affects MUMPS_BLR solver type. Must be called BEFORE first Solve().
   void SetBLRTol(real_t tol) { blr_tol_ = tol; }

   /// v49: Use 2p quadrature instead of 2p+1 for fault face integration.
   /// Tests whether the extra quad point causes instability at p>=2.
   void SetMatchQuadOrder(bool v) { match_quad_order_ = v; }

   /// v49: Diagnostic - compare CalcOrtho normals with FaultBasis normals.
   void SetDiagNormals(bool v) { diag_normals_ = v; }

   /// v49: Diagnostic - dump traction values at first evaluation (zero slip).
   void SetDiagFirstTraction(bool v) { diag_first_traction_ = v; }

   /// v50a: Scale IP penalty by this factor (1.0 = default, <1.0 = reduced).
   /// Used to diagnose whether over-stiff penalty causes nucleation failure at high p.
   void SetPenaltyFactor(real_t f) { penalty_factor_ = f; }

   /// v51: Dump per-component traction (global x,y,z + local dip,strike) at first
   /// non-zero-slip ComputeTraction call. Isolates cross-component coupling source.
   void SetDiagDipTraction(bool v) { diag_dip_traction_ = v; }

   /// v51: Dump z-component of displacement at fault face quad points after solve.
   /// Tests K-f consistency: for pure strike-slip, u_z should be exactly 0.
   void SetDiagUzFault(bool v) { diag_uz_fault_ = v; }

   /// v52: Traction coherence diagnostic.
   /// Tests whether ComputeTraction is consistent with the solve by measuring:
   ///   1. Penalty correction magnitude (η*[[u]]-slip) — should be ~0 if solve is exact
   ///   2. Dip contamination from stress vs penalty separately
   ///   3. Solver fault residual |[[u]] - slip| per quad point
   /// Triggers once when slip is non-trivial, prints summary, then done.
   void SetDiagTractionCoherence(bool v) { diag_traction_coherence_ = v; }

   /// v52: Enable RHS z-component diagnostic (fires once after first non-trivial slip)
   void SetDiagRhsZ(bool v) { diag_rhs_z_ = v; }

   /// Diagnostic: output per-QP traction data at tip faces, matching Tandem [TND-TQ] format.
   /// Fires once at first ComputeTraction call. Each rank prints its own tip faces.
   void SetDiagTndTQ(bool v) { diag_tnd_tq_ = v; }

   /// v50g: Set face DOF node type for FaceQuadrature.
   /// Must be called BEFORE Init() (which creates FaceQuadrature).
   /// BasisType::GaussLobatto (default), BasisType::ClosedUniform, etc.
   void SetFaceBasisType(int bt) { face_basis_type_ = bt; }

private:
   MeshType &mesh_;
   int order_;
   real_t lambda_val_, mu_val_;
   real_t Vp_, Wf_, lf_;
   DGMethod method_;
   SolverType solver_type_;
   BCMode bc_mode_;
   bool check_residual_;  // Post-solve residual check
   bool diag_traction_decomp_ = false;  // Print traction decomposition (stress vs penalty)
   real_t blr_tol_ = 1e-12;  // MUMPS-BLR factorization tolerance (v48: tightened from 1e-10)
   mutable int diag_face_call_ = 0;  // Face consistency diagnostic: trigger on call #2 (non-zero slip)

   // v49 Phase 1 diagnostic flags
   bool match_quad_order_ = false;       // Use 2p instead of 2p+1 quadrature
   bool diag_normals_ = false;           // Compare CalcOrtho vs FaultBasis normals
   bool diag_first_traction_ = false;    // Dump traction at first zero-slip evaluation
   mutable bool diag_normals_done_ = false;
   mutable bool diag_first_traction_done_ = true;  // v58: disabled by default, was flooding output
   real_t penalty_factor_ = 1.0;  // v50a: scale IP penalty (1.0=default)
   mutable bool used_coord_fallback_ = false;  // Set when coordinate-based fault detection is used
   bool has_fault_attr_ = false;         // True when the mesh globally contains fault attr 3
   bool diag_dip_traction_ = false;     // v51: dump per-component traction (global xyz)
   mutable bool diag_dip_traction_done_ = false;
   bool diag_uz_fault_ = false;         // v51: dump u_z at fault faces after solve
   mutable bool diag_uz_fault_done_ = false;
   bool diag_traction_coherence_ = false;  // v52: traction coherence diagnostic
   mutable bool diag_traction_coherence_done_ = false;
   bool diag_rhs_z_ = false;               // v52: dump f_z components of RHS
   mutable bool diag_rhs_z_done_ = false;
   mutable bool diag_matrix_norm_done_ = false;  // v57 MPI diagnostic
   mutable bool diag_slip_embed_done_ = true;    // v58 FaultBasis diagnostic (disabled by default)
   bool diag_tnd_tq_ = false;                     // Per-QP tip traction (Tandem [TND-TQ] comparison)
   mutable bool diag_tnd_tq_done_ = false;
   int face_basis_type_ = BasisType::GaussLobatto;  // v50g: face DOF node type

   void ComputeTractionImpl(const GridFuncType &displacement,
                            const Vector &slip_bc,
                            Vector &traction,
                            Vector *normal_traction,
                            Vector *traction_stress_out,
                            Vector *traction_correction_out,
                            Vector *jump_residual_out);

   // Tag-based fault face detection (matches Tandem's Physical Surface approach)
   Array<int> fault_tagged_faces_;      // Interior face indices from mesh tags
   std::set<long> fault_face_keys_;     // Element-pair keys for fast lookup
   std::set<int> fault_shared_tagged_;  // Shared face indices from mesh tags

   // Dirichlet interior faces: interior faces with attr 5 that are NOT fault faces.
   // In Tandem's BP5 mesh, Physical Surface(5) tags all far-field boundaries including
   // Y=0 faces outside the fault region. These are interior faces in MFEM and need
   // Dirichlet RHS contributions via the skeleton (two-element) pattern.
   Array<int> dirichlet_interior_faces_;
   Array<int> dirichlet_shared_faces_;

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
   Array<int> fault_interior_faces_;
   Array<int> fault_shared_faces_;
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

   /// Build tag-based fault face lookup from mesh boundary element attributes.
   ///
   /// Uses MFEM's direct GetBdrElementFaceIndex(be) to map each attr-3
   /// boundary element to its face index in O(1), avoiding fragile vertex-set
   /// matching. This is the closest MFEM equivalent to Tandem's direct
   /// facet-tag import (GlobalSimplexMeshBuilder.cpp:72).
   void BuildFaultTaggedFaces()
   {
      fault_tagged_faces_.SetSize(0);
      fault_face_keys_.clear();
      fault_shared_tagged_.clear();

      int fault_attr = 3;
      bool has_fault_attr_local = false;
      for (int i = 0; i < mesh_.bdr_attributes.Size(); i++)
      {
         if (mesh_.bdr_attributes[i] == fault_attr)
         {
            has_fault_attr_local = true;
            break;
         }
      }

      has_fault_attr_ = has_fault_attr_local;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int local = has_fault_attr_local ? 1 : 0;
         int global = 0;
         MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, mesh_.GetComm());
         has_fault_attr_ = (global != 0);
#endif
      }
      if (!has_fault_attr_) { return; }

      // Build reverse map: local face index → shared face index.
      // Needed to tag shared faces directly from boundary elements.
      std::unordered_map<int, int> lface_to_sface;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
         {
            lface_to_sface[mesh_.GetSharedFace(sf)] = sf;
         }
#endif
      }

      // Collect fault face indices directly from boundary elements using
      // GetBdrElementFaceIndex — O(1) per boundary element, no vertex
      // matching. Each attr-3 boundary element maps to exactly one face.
      std::set<int> fault_face_set;
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         if (mesh_.GetBdrAttribute(be) != fault_attr) { continue; }
         int face_idx = mesh_.GetBdrElementFaceIndex(be);
         fault_face_set.insert(face_idx);
      }

      // Classify each tagged face as interior or shared
      for (int face_idx : fault_face_set)
      {
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face_idx);
         if (FTr != nullptr)
         {
            // Interior face: both adjacent elements are local
            fault_tagged_faces_.Append(face_idx);
            int e1 = FTr->Elem1No;
            int e2 = FTr->Elem2No;
            long key = (long)std::min(e1, e2) * mesh_.GetNE()
                       + std::max(e1, e2);
            fault_face_keys_.insert(key);
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            // Check if this face is a shared face at a partition boundary
            auto it = lface_to_sface.find(face_idx);
            if (it != lface_to_sface.end())
            {
               fault_shared_tagged_.insert(it->second);
            }
         }
      }

      // v58 fix: Exchange shared fault face tags between neighboring ranks.
      //
      // Internal boundary elements (attr=3) only exist on ONE rank per
      // shared face.  The rank without the boundary element never detects
      // the face as a fault face, causing its elem2 RHS contribution to
      // be silently dropped in AssembleSlipContributionIPShared.
      //
      // Fix: identify each face by its sorted global vertex IDs (canonical
      // triplet).  Each rank broadcasts the vertex signatures of its
      // detected shared fault faces.  Other ranks match their shared faces
      // against this global set.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvert_tag;
         mesh_.GetGlobalVertexIndices(gvert_tag);

         // Build canonical vertex keys for locally-detected shared fault faces
         std::vector<FaceVertexKey> local_fault_keys;
         for (int sf : fault_shared_tagged_)
         {
            int lf = mesh_.GetSharedFace(sf);
            local_fault_keys.push_back(MakeFaceKey(lf, gvert_tag));
         }

         // Allgather: collect all fault face keys from all ranks
         int local_count = static_cast<int>(local_fault_keys.size());
         int nranks = 1;
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         std::vector<int> recv_counts(nranks), displs(nranks);
         MPI_Allgather(&local_count, 1, MPI_INT,
                       recv_counts.data(), 1, MPI_INT, mesh_.GetComm());
         int total = 0;
         for (int r = 0; r < nranks; r++)
         {
            displs[r] = total;
            total += recv_counts[r];
         }

         // Pack as flat array (3 HYPRE_BigInt per face)
         std::vector<HYPRE_BigInt> local_flat(3 * local_count);
         for (int i = 0; i < local_count; i++)
         {
            local_flat[3*i]   = local_fault_keys[i].v[0];
            local_flat[3*i+1] = local_fault_keys[i].v[1];
            local_flat[3*i+2] = local_fault_keys[i].v[2];
         }
         std::vector<int> recv3(nranks), disp3(nranks);
         for (int r = 0; r < nranks; r++)
         {
            recv3[r] = 3 * recv_counts[r];
            disp3[r] = 3 * displs[r];
         }
         std::vector<HYPRE_BigInt> all_flat(3 * total);
         MPI_Allgatherv(local_flat.data(), 3 * local_count, HYPRE_MPI_BIG_INT,
                        all_flat.data(), recv3.data(), disp3.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Build set of all known fault face keys
         std::set<FaceVertexKey> global_fault_keys;
         for (int i = 0; i < total; i++)
         {
            FaceVertexKey k;
            k.v[0] = all_flat[3*i];
            k.v[1] = all_flat[3*i+1];
            k.v[2] = all_flat[3*i+2];
            global_fault_keys.insert(k);
         }

         // Check each of MY shared faces against the global set
         for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
         {
            if (fault_shared_tagged_.count(sf) > 0) { continue; }
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey key = MakeFaceKey(lf, gvert_tag);
            if (global_fault_keys.count(key) > 0)
            {
               fault_shared_tagged_.insert(sf);
            }
         }
#endif
      }
   }

   /// Build Dirichlet interior face list from mesh boundary element attributes.
   ///
   /// Uses the same direct GetBdrElementFaceIndex approach as
   /// BuildFaultTaggedFaces. For attr 5 (far-field Dirichlet), excluding
   /// faces already tagged as fault.
   void BuildDirichletInteriorFaces()
   {
      dirichlet_interior_faces_.SetSize(0);
      dirichlet_shared_faces_.SetSize(0);

      int dirichlet_attr = 5;
      bool has_dirichlet_attr = false;
      for (int i = 0; i < mesh_.bdr_attributes.Size(); i++)
      {
         if (mesh_.bdr_attributes[i] == dirichlet_attr)
         {
            has_dirichlet_attr = true;
            break;
         }
      }
      if (!has_dirichlet_attr) { return; }

      // Build reverse map: local face → shared face index
      std::unordered_map<int, int> lface_to_sface;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
         {
            lface_to_sface[mesh_.GetSharedFace(sf)] = sf;
         }
#endif
      }

      // Direct face lookup from boundary elements
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         if (mesh_.GetBdrAttribute(be) != dirichlet_attr) { continue; }
         int face_idx = mesh_.GetBdrElementFaceIndex(be);

         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face_idx);
         if (FTr != nullptr)
         {
            // Exclude fault faces by element-pair key
            int e1 = FTr->Elem1No;
            int e2 = FTr->Elem2No;
            long key = (long)std::min(e1, e2) * mesh_.GetNE()
                       + std::max(e1, e2);
            if (fault_face_keys_.count(key) > 0) { continue; }

            dirichlet_interior_faces_.Append(face_idx);
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            auto it = lface_to_sface.find(face_idx);
            if (it != lface_to_sface.end())
            {
               if (fault_shared_tagged_.count(it->second) > 0) { continue; }
               dirichlet_shared_faces_.Append(it->second);
            }
         }
      }
   }

   void SetupFaultInfo()
   {
      fault_interior_faces_.SetSize(0);

      // Build tag-based fault face lookup from mesh Physical Surface attributes.
      // This matches Tandem's BC-based face classification.
      BuildFaultTaggedFaces();

      // Build Dirichlet interior face list (must come after BuildFaultTaggedFaces
      // since it uses fault_face_keys_ to exclude fault faces)
      BuildDirichletInteriorFaces();

      // Detect fault faces using tags (or coordinate fallback)
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

      num_fault_faces_ = fault_interior_faces_.Size() + fault_shared_faces_.Size();

      if (has_fault_attr_)
      {
         // Validate: the number of recovered tagged fault faces should match
         // the number of boundary elements with attr 3 (globally).
         int local_tagged_bdr = 0;
         for (int be = 0; be < mesh_.GetNBE(); be++)
         {
            if (mesh_.GetBdrAttribute(be) == 3) { local_tagged_bdr++; }
         }
         int global_tagged_bdr = local_tagged_bdr;
         // Count BOTH interior and shared tagged faces
         int local_tagged_total = fault_tagged_faces_.Size()
            + static_cast<int>(fault_shared_tagged_.size());
         int global_tagged_total = local_tagged_total;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, &global_tagged_bdr, 1, MPI_INT,
                          MPI_SUM, mesh_.GetComm());
            MPI_Allreduce(MPI_IN_PLACE, &global_tagged_total, 1, MPI_INT,
                          MPI_SUM, mesh_.GetComm());
#endif
         }
         bool is_root = true;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            int rank;
            MPI_Comm_rank(mesh_.GetComm(), &rank);
            is_root = (rank == 0);
#endif
         }
         // In MPI, shared faces are counted on both ranks, so
         // global_tagged_total may exceed global_tagged_bdr.
         // A genuine problem is when tagged_total < tagged_bdr.
         if (is_root && global_tagged_bdr > 0 &&
             global_tagged_total < global_tagged_bdr)
         {
            mfem::out << "  WARNING: Tag-based fault recovery: "
                      << global_tagged_total
                      << " tagged faces (interior+shared) vs "
                      << global_tagged_bdr << " boundary elements with attr 3. "
                      << "Some fault faces may be missing.\n";
         }
         if (is_root && global_tagged_bdr > 0)
         {
            mfem::out << "  Tag-based fault faces: "
                      << fault_tagged_faces_.Size() << " interior + "
                      << fault_shared_tagged_.size() << " shared (this rank), "
                      << global_tagged_bdr << " boundary elements globally\n";
         }
      }

      // Emit coordinate-fallback warning once, on root only
      if (used_coord_fallback_ && !has_fault_attr_)
      {
         bool is_root = true;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            int rank;
            MPI_Comm_rank(mesh_.GetComm(), &rank);
            is_root = (rank == 0);
#endif
         }
         if (is_root)
         {
            mfem::out << "  WARNING: Using coordinate-based fault detection "
                      << "(mesh has no fault attr 3). This does not follow "
                      << "Tandem's tag-based classification and can hide mesh "
                      << "tagging errors.\n";
         }
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
         // For p=1: nbf_per_face_ == nv (3 vertices = 3 DOFs)
         for (int k = 0; k < nv && k < nbf_per_face_; k++)
         {
            canonical_to_local_perm_[fi * nbf_per_face_ + k] = gid_idx[k].second;
         }
         // For higher-order DOFs beyond vertices (p>=2): identity for now
         for (int k = nv; k < nbf_per_face_; k++)
         {
            canonical_to_local_perm_[fi * nbf_per_face_ + k] = k;
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

   bool IsFaultFace3D(FaceElementTransformations *FTr) const
   {
      // Follow Tandem: if the mesh carries the fault tag, trust tag-based
      // classification and do not silently fall back to coordinates.
      if (has_fault_attr_)
      {
         int e1 = FTr->Elem1No;
         int e2 = FTr->Elem2No;
         long key = (long)std::min(e1, e2) * mesh_.GetNE() + std::max(e1, e2);
         return fault_face_keys_.count(key) > 0;
      }

      // Fallback: coordinate-based detection (Tandem convention only).
      // This can hide bad mesh tagging — prefer tag-based detection.
      used_coord_fallback_ = true;

      // Fault at Y=0, X in [-lf/2, lf/2], Z in [-Wf, 0]
      const IntegrationPoint &ip = Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector center(3);
      FTr->Face->Transform(ip, center);

      const real_t tol = 1e-10 * std::max(Wf_, 1.0);

      return std::abs(center(1)) < tol
          && std::abs(center(0)) <= lf_ / 2.0 + tol
          && center(2) >= -Wf_ - tol
          && center(2) <= tol;
   }

   bool IsFaultFace3DShared(int shared_face) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         // Follow Tandem: if the mesh carries the fault tag, trust tag-based
         // classification and do not silently fall back to coordinates.
         if (has_fault_attr_)
         {
            return fault_shared_tagged_.count(shared_face) > 0;
         }

         // Fallback: coordinate-based (Tandem convention only)
         // This can hide bad mesh tagging — see warning in IsFaultFace3D.
         // Fault at Y=0, X in [-lf/2,lf/2], Z in [-Wf,0]
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(shared_face);
         if (FTr == nullptr) { return false; }

         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->Face->SetIntPoint(&ip);
         Vector center(3);
         FTr->Face->Transform(ip, center);

         const real_t tol = 1e-10 * std::max(Wf_, 1.0);

         return std::abs(center(1)) < tol
             && std::abs(center(0)) <= lf_ / 2.0 + tol
             && center(2) >= -Wf_ - tol && center(2) <= tol;
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

      // Set up solver operator
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         // v57 face-count diagnostic: how many faces does each assembly path see?
         {
            int local_interior = mesh_.GetNumFaces();
            int interior_with_nbr = 0;
            for (int f = 0; f < local_interior; f++)
            {
               auto *FTr = mesh_.GetInteriorFaceTransformations(f);
               if (FTr) { interior_with_nbr++; }
            }
            int local_shared = mesh_.GetNSharedFaces();
            int local_bdr = mesh_.GetNBE();
            int total_local_faces = mesh_.GetNumFaces();

            int global_interior = 0, global_shared = 0;
            int global_bdr = 0, global_total = 0;
            MPI_Reduce(&interior_with_nbr, &global_interior, 1, MPI_INT,
                       MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_shared, &global_shared, 1, MPI_INT,
                       MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_bdr, &global_bdr, 1, MPI_INT,
                       MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&total_local_faces, &global_total, 1, MPI_INT,
                       MPI_SUM, 0, mesh_.GetComm());

            // Local sparse matrix NNZ (before ParallelAssemble)
            int local_mat_nnz = cached_a_->SpMat().NumNonZeroElems();
            long long global_mat_nnz = 0;
            long long ll_nnz = local_mat_nnz;
            MPI_Reduce(&ll_nnz, &global_mat_nnz, 1, MPI_LONG_LONG,
                       MPI_SUM, 0, mesh_.GetComm());

            int rank = 0;
            MPI_Comm_rank(mesh_.GetComm(), &rank);
            // [MPI-DIAG] face assembly counts — disabled by default (v58)
            (void)global_interior; (void)global_shared;
            (void)global_bdr; (void)global_total; (void)global_mat_nnz;
         }

         cached_Ah_.SetType(Operator::Hypre_ParCSR);
         cached_a_->ParallelAssemble(cached_Ah_);

         // v57 diagnostic: P matrix and HypreParMatrix global dimensions
         {
            auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
            auto *Kh = cached_Ah_.As<HypreParMatrix>();
            HYPRE_BigInt glob_rows = Kh->GetGlobalNumRows();
            HYPRE_BigInt glob_cols = Kh->GetGlobalNumCols();
            HYPRE_BigInt glob_nnz = Kh->NNZ();

            // Check P matrix: for DG, should be identity
            int local_dofs = pfes->GetVSize();
            int true_dofs = pfes->GetTrueVSize();

            int rank = 0;
            MPI_Comm_rank(mesh_.GetComm(), &rank);
            // Compute sum(|K_ij|) (L1 norm of all entries) — partition-independent test
            hypre_ParCSRMatrix *hA = (hypre_ParCSRMatrix*)(*Kh);
            hypre_CSRMatrix *d = hypre_ParCSRMatrixDiag(hA);
            hypre_CSRMatrix *o = hypre_ParCSRMatrixOffd(hA);
            double *dd = hypre_CSRMatrixData(d);
            double *od = hypre_CSRMatrixData(o);
            int dnnz = hypre_CSRMatrixNumNonzeros(d);
            int onnz = hypre_CSRMatrixNumNonzeros(o);
            real_t local_abs_sum = 0.0;
            for (int ii = 0; ii < dnnz; ii++) local_abs_sum += std::abs(dd[ii]);
            for (int ii = 0; ii < onnz; ii++) local_abs_sum += std::abs(od[ii]);
            real_t global_abs_sum = 0.0;
            MPI_Reduce(&local_abs_sum, &global_abs_sum, 1, MPI_DOUBLE,
                       MPI_SUM, 0, mesh_.GetComm());

            // [MPI-DIAG] HypreParMatrix dims — disabled by default (v58)
            (void)glob_rows; (void)glob_cols; (void)glob_nnz;
            (void)global_abs_sum; (void)local_dofs; (void)true_dofs;
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
      diag_face_call_++;

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

         const auto &basis_slip = fault_basis_.GetBasis(fi);
         real_t sign = basis_slip.sign_flipped ? 1.0 : -1.0;

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
         if (!basis_slip.qp_data.empty())
         {
            // Per-QP tangent embedding: interpolate tangential components,
            // then embed using per-QP tangent frame
            Vector slip_tang(2 * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               slip_tang(0 * nbf + kk) = slip_bc(2 * dof_idx);      // dip
               slip_tang(1 * nbf + kk) = slip_bc(2 * dof_idx + 1);  // strike
            }
            Vector slip_tang_q;
            face_quad_->InterpolateToQuadPoints(2, slip_tang, slip_tang_q);
            int nqp = slip_tang_q.Size() / 2;
            delta_u_quad.SetSize(dim * nqp);
            for (int q = 0; q < nqp; q++)
            {
               real_t sl_q[2] = {slip_tang_q(q), slip_tang_q(nqp + q)};
               real_t du[3];
               fault_basis_.EmbedSlipQP(fi, q, sl_q, du);
               for (int c = 0; c < dim; c++)
                  delta_u_quad(c * nqp + q) = sign * du[c];
            }
         }
         else
         {
            // Fallback: centroid tangents (BR2 path)
            Vector delta_u_nodal(dim * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               real_t slip_local[2] = {slip_bc(2 * dof_idx),
                                       slip_bc(2 * dof_idx + 1)};
               real_t du[3];
               fault_basis_.EmbedSlip(fi, slip_local, du);
               for (int c = 0; c < dim; c++)
                  delta_u_nodal(c * nbf + kk) = du[c];
            }
            face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal, delta_u_quad);
            for (int i = 0; i < delta_u_quad.Size(); i++)
               delta_u_quad(i) *= sign;
         }

         // v58 diagnostic: print FaultBasis + delta_u_quad for first 3 interior faces
         if (!diag_slip_embed_done_ && fi < 3)
         {
            const IntegrationPoint &ip_c =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_c);
            Vector fc(3);
            FTr->Face->Transform(ip_c, fc);
            int rank = 0;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               MPI_Comm_rank(mesh_.GetComm(), &rank);
#endif
            }
            int nqp_diag = delta_u_quad.Size() / dim;
            mfem::out << "[SLIP-EMBED] INTERIOR fi=" << fi
                      << " rank=" << rank
                      << " centroid=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                      << " sign_flipped=" << basis_slip.sign_flipped
                      << " n=(" << basis_slip.normal[0] << ","
                      << basis_slip.normal[1] << "," << basis_slip.normal[2] << ")"
                      << " t1=(" << basis_slip.tangent1[0] << ","
                      << basis_slip.tangent1[1] << "," << basis_slip.tangent1[2] << ")"
                      << " t2=(" << basis_slip.tangent2[0] << ","
                      << basis_slip.tangent2[1] << "," << basis_slip.tangent2[2] << ")"
                      << " du_q0=(" << delta_u_quad(0) << ","
                      << delta_u_quad(nqp_diag) << ","
                      << delta_u_quad(2*nqp_diag) << ")"
                      << "\n";
         }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);

         // Assemble slip RHS using the combined integrator
         Vector elvec1, elvec2;
         slip_integrator.AssembleSlipFaceRHS(
            *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);

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
            const auto &basis_slip = fault_basis_.GetBasis(slip_idx);
            real_t sign = basis_slip.sign_flipped ? 1.0 : -1.0;

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
            if (!basis_slip.qp_data.empty())
            {
               Vector slip_tang(2 * nbf);
               for (int kk = 0; kk < nbf; kk++)
               {
                  int dof_idx = slip_idx * nbf + kk;
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
                  fault_basis_.EmbedSlipQP(slip_idx, q, sl_q, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_quad(c * nqp + q) = sign * du[c];
               }
            }
            else
            {
               Vector delta_u_nodal(dim * nbf);
               for (int kk = 0; kk < nbf; kk++)
               {
                  int dof_idx = slip_idx * nbf + kk;
                  real_t slip_local[2] = {slip_bc(2 * dof_idx),
                                          slip_bc(2 * dof_idx + 1)};
                  real_t du[3];
                  fault_basis_.EmbedSlip(slip_idx, slip_local, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_nodal(c * nbf + kk) = du[c];
               }
               face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal, delta_u_quad);
               for (int j = 0; j < delta_u_quad.Size(); j++)
                  delta_u_quad(j) *= sign;
            }

            // v58 diagnostic: print FaultBasis + delta_u_quad for ALL shared fault faces
            if (!diag_slip_embed_done_)
            {
               const IntegrationPoint &ip_c =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->Face->SetIntPoint(&ip_c);
               Vector fc(3);
               FTr->Face->Transform(ip_c, fc);
               int rank = 0;
               MPI_Comm_rank(mesh_.GetComm(), &rank);
               int nqp_diag = delta_u_quad.Size() / dim;
               mfem::out << "[SLIP-EMBED] SHARED i=" << i
                         << " rank=" << rank
                         << " centroid=(" << fc(0) << "," << fc(1)
                         << "," << fc(2) << ")"
                         << " sign_flipped=" << basis_slip.sign_flipped
                         << " n=(" << basis_slip.normal[0] << ","
                         << basis_slip.normal[1] << ","
                         << basis_slip.normal[2] << ")"
                         << " t1=(" << basis_slip.tangent1[0] << ","
                         << basis_slip.tangent1[1] << ","
                         << basis_slip.tangent1[2] << ")"
                         << " t2=(" << basis_slip.tangent2[0] << ","
                         << basis_slip.tangent2[1] << ","
                         << basis_slip.tangent2[2] << ")"
                         << " du_q0=(" << delta_u_quad(0) << ","
                         << delta_u_quad(nqp_diag) << ","
                         << delta_u_quad(2*nqp_diag) << ")"
                         << "\n";
            }

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

            // Assemble using combined integrator — only use elvec1 (local elem)
            Vector elvec1, elvec2;
            slip_integrator.AssembleSlipFaceRHS(
               *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
         // Only mark done if we actually printed (had non-zero slip)
         {
            real_t slip_max_local = slip_bc.Normlinf();
            if (!diag_slip_embed_done_ && slip_max_local > 1e-20)
            {
               diag_slip_embed_done_ = true;
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
      // BP5 Dirichlet loading: u = (0, sgn(x)·Vp·t/2, 0)
      // Applied on Dirichlet-marked boundaries only (controlled by BCMode).
      // sgn(x) determined from face centroid x-coordinate.
      //
      // DG Dirichlet BC contribution:
      //   b[k,i] += c0 * [σ(φ_k e_i)·n]_u * u_D_u * (1/detJ)
      //           + penalty * φ_k * u_D_i

      if (std::abs(time * Vp_) < 1e-30) { return; }

      int dim = 3;

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

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
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Face->SetIntPoint(&ipq);
               FTr->Face->Transform(ipq, phys);
               real_t y = phys(1);
               // Tandem bp5.lua boundary(x,y,z,t):
               real_t Vh = Vp_ * time;
               if (y > 1.0) { Vh *= 0.5; }
               else if (y < -1.0) { Vh *= -0.5; }
               u_D_3d(0 * nq_dir + q) = Vh;
            }

            Vector elvec_dir;
            dir_integ.AssembleBoundaryFaceRHS(*fe, *FTr, u_D_3d, elvec_dir);

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
            if (fc_br2(1) > 1.0) { Vh_br2 *= 0.5; }
            else if (fc_br2(1) < -1.0) { Vh_br2 *= -0.5; }
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

         // v55: Evaluate Tandem's boundary function at each quad point,
         // matching bp5.lua boundary(x,y,z,t) exactly:
         //   y > 1:  u_D = (Vp*t/2, 0, 0)
         //   y < -1: u_D = (-Vp*t/2, 0, 0)
         //   else:   u_D = (Vp*t, 0, 0)
         //
         // The DG skeleton RHS uses this as f_q (prescribed data at the face).
         // Previously we computed u_D from element centroid Y-signs, which is
         // mathematically equivalent but evaluates at different points.

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
            // v55: Use combined integrator's AssembleSlipFaceRHS for
            // Dirichlet interior faces (skeleton pattern, same as fault slip)
            DGElasticityIPCombinedIntegrator dir_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_dir = IntRules.Get(
               FTr->FaceGeom, quad_order_dir);
            int nq_dir = ir_dir.GetNPoints();

            // v55: Evaluate Tandem's boundary(x,y,z,t) at each quad point
            Vector u_D_3d(dim * nq_dir);
            u_D_3d = 0.0;
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Face->SetIntPoint(&ipq);
               FTr->Face->Transform(ipq, phys);
               real_t y = phys(1);
               // Tandem bp5.lua boundary(x,y,z,t):
               real_t Vh = Vp_ * time;
               if (y > 1.0) { Vh *= 0.5; }
               else if (y < -1.0) { Vh *= -0.5; }
               // else: Vh = Vp*t (full rate for |y| <= 1)
               u_D_3d(0 * nq_dir + q) = Vh;  // X component
               // Y and Z components = 0
            }

            Vector ev1, ev2;
            dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);

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
            if (fc_br2(1) > 1.0) { Vh_br2 *= 0.5; }
            else if (fc_br2(1) < -1.0) { Vh_br2 *= -0.5; }
            real_t u_D_int[3] = {Vh_br2, 0.0, 0.0};

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
               // v55: Use combined integrator (skeleton pattern, only elem1)
               DGElasticityIPCombinedIntegrator dir_integ(
                  lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

               int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
               const IntegrationRule &ir_dir = IntRules.Get(
                  FTr->FaceGeom, quad_order_dir);
               int nq_dir = ir_dir.GetNPoints();

               // v55: Evaluate Tandem's boundary(x,y,z,t) at each quad point
               Vector u_D_3d(dim * nq_dir);
               u_D_3d = 0.0;
               for (int q = 0; q < nq_dir; q++)
               {
                  const IntegrationPoint &ipq = ir_dir.IntPoint(q);
                  FTr->SetAllIntPoints(&ipq);
                  Vector phys(dim);
                  FTr->Face->SetIntPoint(&ipq);
                  FTr->Face->Transform(ipq, phys);
                  real_t y = phys(1);
                  real_t Vh = Vp_ * time;
                  if (y > 1.0) { Vh *= 0.5; }
                  else if (y < -1.0) { Vh *= -0.5; }
                  u_D_3d(0 * nq_dir + q) = Vh;
               }

               Vector ev1, ev2;
               dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);
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
               if (fc_br2(1) > 1.0) { Vh_br2 *= 0.5; }
               else if (fc_br2(1) < -1.0) { Vh_br2 *= -0.5; }
               real_t u_D_int[3] = {Vh_br2, 0.0, 0.0};

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

   // Add slip contributions (interior + shared faces)
   real_t rhs_before_slip = rhs.Normlinf();

   // v57 RHS decomposition diagnostic: capture interior-only and shared-only norms
   Vector rhs_after_interior;
   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);

      // Snapshot after interior-only (before shared)
      if (!diag_matrix_norm_done_)
      {
         rhs_after_interior.SetSize(rhs.Size());
         rhs_after_interior = rhs;
      }

      AssembleSlipContributionIPShared(rhs, slip_bc,
                                       fault_interior_faces_.Size());
   }
   else
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
      if (!diag_matrix_norm_done_)
      {
         rhs_after_interior.SetSize(rhs.Size());
         rhs_after_interior = rhs;
      }
      AssembleSlipContributionBR2Shared(rhs, slip_bc,
                                        fault_interior_faces_.Size());
   }
   real_t rhs_after_slip = rhs.Normlinf();

   // Snapshot total slip RHS (interior + shared, before Dirichlet)
   Vector rhs_total_slip;
   if (!diag_matrix_norm_done_ && rhs_after_slip > 1e-30)
   {
      rhs_total_slip.SetSize(rhs.Size());
      rhs_total_slip = rhs;
   }

   // v52: RHS z-component diagnostic — capture f_z after slip, before Dirichlet
   // byNODES ordering: rhs[0..N-1]=x, rhs[N..2N-1]=y, rhs[2N..3N-1]=z
   int N_scalar = fes_->GetNDofs();
   Vector rhs_slip_snapshot;

   // Synchronize trigger across all ranks to avoid MPI deadlock
   // (ranks without fault DOFs may have rhs_after_slip == 0)
   bool diag_rhs_z_trigger = (diag_rhs_z_ && !diag_rhs_z_done_ && rhs_after_slip > 0.0);
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      int local_trigger = diag_rhs_z_trigger ? 1 : 0;
      int global_trigger = 0;
      MPI_Allreduce(&local_trigger, &global_trigger, 1, MPI_INT, MPI_MAX,
                    mesh_.GetComm());
      diag_rhs_z_trigger = (global_trigger > 0);
#endif
   }

   if (diag_rhs_z_trigger)
   {
      rhs_slip_snapshot.SetSize(rhs.Size());
      rhs_slip_snapshot = rhs;
   }

   // Add Dirichlet loading
   AssembleDirichletLoading(rhs, time);

   // v52: RHS z-component diagnostic — fire once after first non-trivial RHS
   if (diag_rhs_z_trigger)
   {
      diag_rhs_z_done_ = true;

      // Helper: compute per-component L2 norms (local)
      auto component_norms = [&](const Vector &v, real_t &nx, real_t &ny, real_t &nz)
      {
         nx = ny = nz = 0.0;
         for (int i = 0; i < N_scalar; i++)
         {
            nx += v(i) * v(i);
            ny += v(N_scalar + i) * v(N_scalar + i);
            nz += v(2 * N_scalar + i) * v(2 * N_scalar + i);
         }
      };

      // Slip-only contribution
      real_t slip_nx, slip_ny, slip_nz;
      component_norms(rhs_slip_snapshot, slip_nx, slip_ny, slip_nz);

      // Dirichlet-only contribution (total minus slip)
      Vector rhs_diri(rhs.Size());
      subtract(rhs, rhs_slip_snapshot, rhs_diri);
      real_t diri_nx, diri_ny, diri_nz;
      component_norms(rhs_diri, diri_nx, diri_ny, diri_nz);

      // Total RHS
      real_t tot_nx, tot_ny, tot_nz;
      component_norms(rhs, tot_nx, tot_ny, tot_nz);

      // Also compute max absolute z-values
      real_t slip_zmax = 0.0, diri_zmax = 0.0, tot_zmax = 0.0;
      for (int i = 0; i < N_scalar; i++)
      {
         slip_zmax = std::max(slip_zmax, std::abs(rhs_slip_snapshot(2 * N_scalar + i)));
         diri_zmax = std::max(diri_zmax, std::abs(rhs_diri(2 * N_scalar + i)));
         tot_zmax  = std::max(tot_zmax,  std::abs(rhs(2 * N_scalar + i)));
      }

      // MPI reduce for global norms
      bool is_root = true;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         is_root = (rank == 0);
         // Sum of squares for L2 norms
         real_t local_vals[9] = {slip_nx, slip_ny, slip_nz,
                                  diri_nx, diri_ny, diri_nz,
                                  tot_nx,  tot_ny,  tot_nz};
         real_t global_vals[9];
         MPI_Allreduce(local_vals, global_vals, 9, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         slip_nx = global_vals[0]; slip_ny = global_vals[1]; slip_nz = global_vals[2];
         diri_nx = global_vals[3]; diri_ny = global_vals[4]; diri_nz = global_vals[5];
         tot_nx  = global_vals[6]; tot_ny  = global_vals[7]; tot_nz  = global_vals[8];
         // Max for max-abs
         real_t local_max[3] = {slip_zmax, diri_zmax, tot_zmax};
         real_t global_max[3];
         MPI_Allreduce(local_max, global_max, 3, MPI_DOUBLE, MPI_MAX,
                       mesh_.GetComm());
         slip_zmax = global_max[0]; diri_zmax = global_max[1]; tot_zmax = global_max[2];
#endif
      }

      // Take sqrt for L2 norms
      slip_nx = std::sqrt(slip_nx); slip_ny = std::sqrt(slip_ny); slip_nz = std::sqrt(slip_nz);
      diri_nx = std::sqrt(diri_nx); diri_ny = std::sqrt(diri_ny); diri_nz = std::sqrt(diri_nz);
      tot_nx  = std::sqrt(tot_nx);  tot_ny  = std::sqrt(tot_ny);  tot_nz  = std::sqrt(tot_nz);

      if (is_root)
      {
         mfem::out << "\n[DIAG-RHS-Z] RHS z-component diagnostic (t=" << time << " s)\n"
            << "  DOF ordering: byNODES, N_scalar=" << N_scalar << "\n"
            << "  Slip contribution (f_slip):\n"
            << "    ||f_x|| = " << slip_nx << "\n"
            << "    ||f_y|| = " << slip_ny << "\n"
            << "    ||f_z|| = " << slip_nz << "\n"
            << "    ||f_z||/||f_x|| = " << (slip_nx > 0 ? slip_nz / slip_nx : 0.0) << "\n"
            << "    max|f_z| = " << slip_zmax << "\n"
            << "  Dirichlet contribution (f_diri):\n"
            << "    ||f_x|| = " << diri_nx << "\n"
            << "    ||f_y|| = " << diri_ny << "\n"
            << "    ||f_z|| = " << diri_nz << "\n"
            << "    ||f_z||/||f_x|| = " << (diri_nx > 0 ? diri_nz / diri_nx : 0.0) << "\n"
            << "    max|f_z| = " << diri_zmax << "\n"
            << "  Total RHS (f_slip + f_diri):\n"
            << "    ||f_x|| = " << tot_nx << "\n"
            << "    ||f_y|| = " << tot_ny << "\n"
            << "    ||f_z|| = " << tot_nz << "\n"
            << "    ||f_z||/||f_x|| = " << (tot_nx > 0 ? tot_nz / tot_nx : 0.0) << "\n"
            << "    max|f_z| = " << tot_zmax << "\n"
            << "  Interpretation:\n"
            << "    If ||f_z||/||f_x|| >> 0: RHS has z-forcing (assembly bug, Hypothesis A)\n"
            << "    If ||f_z||/||f_x|| ~ 0: f_z=0, dip comes from solver (Hypothesis B)\n"
            << std::endl;
      }
   }

   // Diagnostic: check for RHS blowup
   real_t rhs_final = rhs.Normlinf();
   real_t slip_max = slip_bc.Normlinf();
   if (std::isnan(rhs_final))
   {
      int rank = 0;
#ifdef MFEM_USE_MPI
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         MPI_Comm_rank(mesh_.GetComm(), &rank);
      }
#endif
      mfem::out << "[Rank " << rank << "] RHS NaN DETECTED: before_slip="
                << rhs_before_slip << " after_slip=" << rhs_after_slip
                << " final=" << rhs_final
                << " slip_max=" << slip_max
                << " time=" << time << "\n";
   }

   X_ = 0.0;
   B_ = rhs;

   // v57 MPI diagnostic: print ||K|| and ||b|| once on first non-trivial solve
   // Use global trigger to avoid MPI deadlock (ranks without fault DOFs
   // may have zero local RHS).
   {
      bool local_trigger = !diag_matrix_norm_done_ && rhs.Normlinf() > 1e-30;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int lt = local_trigger ? 1 : 0, gt = 0;
         MPI_Allreduce(&lt, &gt, 1, MPI_INT, MPI_MAX, mesh_.GetComm());
         local_trigger = (gt > 0);
#endif
      }
      if (local_trigger && !diag_matrix_norm_done_)
      {
      diag_matrix_norm_done_ = true;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *Kh = cached_Ah_.As<HypreParMatrix>();
         // Global ||b||
         real_t local_b2 = B_ * B_;
         real_t global_b2 = 0.0;
         MPI_Allreduce(&local_b2, &global_b2, 1, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         // Global NNZ and Frobenius norm via hypre
         hypre_ParCSRMatrix *hA = (hypre_ParCSRMatrix*)(*Kh);
         // Diagonal block
         hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(hA);
         hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(hA);
         real_t local_K2 = 0.0;
         int diag_nnz = hypre_CSRMatrixNumNonzeros(diag);
         int offd_nnz = hypre_CSRMatrixNumNonzeros(offd);
         double *diag_data = hypre_CSRMatrixData(diag);
         double *offd_data = hypre_CSRMatrixData(offd);
         for (int i = 0; i < diag_nnz; i++)
            local_K2 += diag_data[i] * diag_data[i];
         for (int i = 0; i < offd_nnz; i++)
            local_K2 += offd_data[i] * offd_data[i];
         real_t global_K2 = 0.0;
         MPI_Allreduce(&local_K2, &global_K2, 1, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         int local_nnz = diag_nnz + offd_nnz;
         long long global_nnz = 0;
         long long ll_nnz = local_nnz;
         MPI_Allreduce(&ll_nnz, &global_nnz, 1, MPI_LONG_LONG, MPI_SUM,
                       mesh_.GetComm());
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         // RHS decomposition: interior-only vs shared vs total (slip only, no Dirichlet)
         real_t local_b_int2 = 0.0, local_b_shared2 = 0.0;
         if (rhs_after_interior.Size() > 0 && rhs_total_slip.Size() > 0)
         {
            local_b_int2 = rhs_after_interior * rhs_after_interior;
            // shared contribution = total_slip - interior
            Vector rhs_shared(rhs_total_slip.Size());
            subtract(rhs_total_slip, rhs_after_interior, rhs_shared);
            local_b_shared2 = rhs_shared * rhs_shared;
         }
         real_t global_b_int2 = 0.0, global_b_shared2 = 0.0;
         MPI_Allreduce(&local_b_int2, &global_b_int2, 1, MPI_DOUBLE,
                       MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&local_b_shared2, &global_b_shared2, 1, MPI_DOUBLE,
                       MPI_SUM, mesh_.GetComm());

         // Local sparse NNZ (before ParallelAssemble)
         int local_sparse_nnz = cached_a_->SpMat().NumNonZeroElems();
         long long global_sparse_nnz = 0;
         long long ll_sparse = local_sparse_nnz;
         MPI_Allreduce(&ll_sparse, &global_sparse_nnz, 1, MPI_LONG_LONG,
                       MPI_SUM, mesh_.GetComm());

         // [MPI-DIAG] Stiffness matrix/RHS norms — disabled by default (v58)
         (void)global_K2; (void)global_b2; (void)global_b_int2;
         (void)global_b_shared2; (void)global_nnz; (void)global_sparse_nnz;
#endif
      }
      }  // if (local_trigger)
   }  // scope for MPI diagnostic

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

   // v52: Solution u_z diagnostic — companion to RHS diagnostic
   if (diag_rhs_z_ && diag_rhs_z_done_)
   {
      // Only fire once: reset the flag to prevent repeat (done_ was set above)
      // Use a static to fire exactly once
      static bool u_z_diag_fired = false;
      if (!u_z_diag_fired)
      {
         u_z_diag_fired = true;
         real_t ux2 = 0.0, uy2 = 0.0, uz2 = 0.0;
         real_t ux_max = 0.0, uz_max = 0.0;
         for (int i = 0; i < N_scalar; i++)
         {
            ux2 += X_(i) * X_(i);
            uy2 += X_(N_scalar + i) * X_(N_scalar + i);
            uz2 += X_(2 * N_scalar + i) * X_(2 * N_scalar + i);
            ux_max = std::max(ux_max, std::abs(X_(i)));
            uz_max = std::max(uz_max, std::abs(X_(2 * N_scalar + i)));
         }

         bool is_root = true;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            int rank;
            MPI_Comm_rank(mesh_.GetComm(), &rank);
            is_root = (rank == 0);
            real_t local_sum[3] = {ux2, uy2, uz2};
            real_t global_sum[3];
            MPI_Allreduce(local_sum, global_sum, 3, MPI_DOUBLE, MPI_SUM,
                          mesh_.GetComm());
            ux2 = global_sum[0]; uy2 = global_sum[1]; uz2 = global_sum[2];
            real_t local_mx[2] = {ux_max, uz_max};
            real_t global_mx[2];
            MPI_Allreduce(local_mx, global_mx, 2, MPI_DOUBLE, MPI_MAX,
                          mesh_.GetComm());
            ux_max = global_mx[0]; uz_max = global_mx[1];
#endif
         }

         real_t ux_norm = std::sqrt(ux2);
         real_t uy_norm = std::sqrt(uy2);
         real_t uz_norm = std::sqrt(uz2);

         if (is_root)
         {
            mfem::out << "[DIAG-RHS-Z] Solution u decomposition (same time step):\n"
               << "  ||u_x|| = " << ux_norm << "  (along-strike)\n"
               << "  ||u_y|| = " << uy_norm << "  (fault-normal)\n"
               << "  ||u_z|| = " << uz_norm << "  (dip/vertical)\n"
               << "  ||u_z||/||u_x|| = " << (ux_norm > 0 ? uz_norm / ux_norm : 0.0) << "\n"
               << "  max|u_x| = " << ux_max << "  max|u_z| = " << uz_max << "\n"
               << "  max|u_z|/max|u_x| = " << (ux_max > 0 ? uz_max / ux_max : 0.0) << "\n"
               << "  Interpretation:\n"
               << "    If f_z=0 but u_z>>0: solver introduces dip (Hypothesis B)\n"
               << "    If f_z>>0 and u_z>>0: assembly bug forces dip (Hypothesis A)\n"
               << std::endl;
         }
      }
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
         if (rhs_final > 1e15)
         {
            int rank = 0;
#ifdef MFEM_USE_MPI
            if constexpr (IsParallelMesh<MeshType>::value)
            {
               MPI_Comm_rank(mesh_.GetComm(), &rank);
            }
#endif
            mfem::out << "[Rank " << rank << "] SOLVER: iters="
                      << iter_solver->GetNumIterations()
                      << " converged=" << iter_solver->GetConverged()
                      << " final_norm=" << iter_solver->GetFinalNorm()
                      << " ||u||_inf=" << X_.Normlinf() << "\n";
         }
      }
   }

   // Diagnostic: check for displacement blowup
   {
      real_t u_max = X_.Normlinf();
      if (u_max > 1e6 || std::isnan(u_max))
      {
         int rank = 0;
#ifdef MFEM_USE_MPI
         if constexpr (IsParallelMesh<MeshType>::value)
         {
            MPI_Comm_rank(mesh_.GetComm(), &rank);
         }
#endif
         mfem::out << "[Rank " << rank << "] DISPLACEMENT BLOWUP: ||u||_inf="
                   << u_max << "\n";
         mfem::out << "[Rank " << rank << "] ||RHS||_inf=" << B_.Normlinf()
                   << "\n";
      }
   }

   displacement = X_;
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
   Vector *normal_traction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       &traction_stress, &traction_correction, &jump_residual);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionImpl(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector *normal_traction,
   Vector *traction_stress_out,
   Vector *traction_correction_out,
   Vector *jump_residual_out)
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

   // v51: Elastic normal traction for sigma_n feedback
   if (normal_traction)
   {
      normal_traction->SetSize(num_fault_dofs_);
      *normal_traction = 0.0;
   }

   // v49: Diagnostic - compare CalcOrtho normals with FaultBasis normals
   if (diag_normals_ && !diag_normals_done_)
   {
      int rank = 0;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         MPI_Comm_rank(mesh_.GetComm(), &rank);
      }

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         auto *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!FTr) { continue; }

         const IntegrationPoint &ip_diag =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->SetAllIntPoints(&ip_diag);
         Vector nor(3);
         CalcOrtho(FTr->Jacobian(), nor);
         real_t nl = nor.Norml2();

         const auto &basis = fault_basis_.GetBasis(fi);
         real_t dot = 0;
         for (int d = 0; d < 3; d++)
         {
            dot += (nor(d) / nl) * basis.normal[d];
         }

         // Print if normals are not aligned (|dot| != 1) or for first few faces
         if (std::abs(std::abs(dot) - 1.0) > 1e-10 || fi < 5)
         {
            mfem::out << "[NOR-DIAG] rank=" << rank << " fi=" << fi
               << " CalcOrtho=(" << nor(0)/nl << "," << nor(1)/nl
               << "," << nor(2)/nl << ")"
               << " basis.n=(" << basis.normal[0] << ","
               << basis.normal[1] << "," << basis.normal[2] << ")"
               << " dot=" << dot
               << (std::abs(std::abs(dot) - 1.0) > 1e-10
                   ? " *** MISMATCH ***" : "")
               << "\n";
         }
      }
      // Also check shared faces if parallel
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         int base_idx = fault_interior_faces_.Size();
         for (int si = 0; si < fault_shared_faces_.Size(); si++)
         {
            int sf = fault_shared_faces_[si];
            auto *FTr = mesh_.GetSharedFaceTransformations(sf);
            if (!FTr) { continue; }

            const IntegrationPoint &ip_diag =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ip_diag);
            Vector nor(3);
            CalcOrtho(FTr->Jacobian(), nor);
            real_t nl = nor.Norml2();

            int trac_idx = base_idx + si;
            const auto &basis = fault_basis_.GetBasis(trac_idx);
            real_t dot = 0;
            for (int d = 0; d < 3; d++)
            {
               dot += (nor(d) / nl) * basis.normal[d];
            }

            if (std::abs(std::abs(dot) - 1.0) > 1e-10 || si < 5)
            {
               mfem::out << "[NOR-DIAG] rank=" << rank
                  << " shared si=" << si
                  << " CalcOrtho=(" << nor(0)/nl << "," << nor(1)/nl
                  << "," << nor(2)/nl << ")"
                  << " basis.n=(" << basis.normal[0] << ","
                  << basis.normal[1] << "," << basis.normal[2] << ")"
                  << " dot=" << dot
                  << (std::abs(std::abs(dot) - 1.0) > 1e-10
                      ? " *** MISMATCH ***" : "")
                  << "\n";
            }
         }
      }
      diag_normals_done_ = true;
   }

   // v52: Traction coherence diagnostic accumulators
   bool coherence_active = diag_traction_coherence_ &&
                            !diag_traction_coherence_done_;
   // Check if slip is non-trivial (skip zero-slip evaluations).
   // IMPORTANT: must be globally consistent (all ranks agree) because the
   // summary section uses MPI collectives. Ranks without fault DOFs have
   // empty slip_bc, so local-only check would deadlock.
   if (coherence_active)
   {
      int local_has_slip = 0;
      for (int i = 0; i < slip_bc.Size(); i++)
      {
         if (std::abs(slip_bc(i)) > 1e-20) { local_has_slip = 1; break; }
      }
      int global_has_slip = local_has_slip;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&local_has_slip, &global_has_slip, 1,
                        MPI_INT, MPI_MAX, mesh_.GetComm());
#endif
      }
      if (!global_has_slip) { coherence_active = false; }
   }
   // Accumulate across all fault faces (interior + shared)
   real_t coh_sum_stress_dip2 = 0.0, coh_sum_stress_strike2 = 0.0;
   real_t coh_sum_corr_dip2 = 0.0, coh_sum_corr_strike2 = 0.0;
   real_t coh_sum_res2 = 0.0;  // ||[[u]] - slip||^2 at quad points
   real_t coh_max_res = 0.0;   // max |[[u]] - slip| component
   real_t coh_max_corr_dip = 0.0;
   real_t coh_max_stress_dip = 0.0;
   int coh_n_qp = 0;

   // Per-face centroid data for station-level reporting
   struct CohFaceData {
      real_t x2, x3;               // face centroid in fault coords
      real_t stress_dip, stress_strike;  // face-averaged stress-only
      real_t corr_dip, corr_strike;      // face-averaged penalty correction
      real_t max_res;              // max |[[u]]-slip| on this face
   };
   std::vector<CohFaceData> coh_face_data;

   // [MFEM-TQ] Pre-compute global vertex indices once (may be collective)
   Array<HYPRE_BigInt> gvert_tq_;
   if (diag_tnd_tq_ && !diag_tnd_tq_done_)
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         mesh_.GetGlobalVertexIndices(gvert_tq_);
      }
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

      // v51: Dump u_z at fault face centroid (one-time diagnostic)
      if (diag_uz_fault_ && !diag_uz_fault_done_)
      {
         // Check if any slip is non-zero
         bool has_slip = false;
         int nbf_check = nbf_per_face_;
         for (int kk = 0; kk < nbf_check; kk++)
         {
            int di = fi * nbf_check + kk;
            if (di * 2 + 1 < slip_bc.Size() &&
                (std::abs(slip_bc(2*di)) > 1e-20 ||
                 std::abs(slip_bc(2*di+1)) > 1e-20))
            { has_slip = true; break; }
         }
         if (has_slip || fi == 0)
         {
            // Evaluate u at face centroid
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            Vector shape1(ndof1), shape2(ndof2);
            fe1->CalcShape(eip1, shape1);
            fe2->CalcShape(eip2, shape2);

            real_t u1[3] = {0,0,0}, u2[3] = {0,0,0};
            for (int c = 0; c < dim; c++)
               for (int k = 0; k < ndof1; k++)
                  u1[c] += shape1(k) * u1_all(c * ndof1 + k);
            for (int c = 0; c < dim; c++)
               for (int k = 0; k < ndof2; k++)
                  u2[c] += shape2(k) * u2_all(c * ndof2 + k);

            Vector fc(3);
            FTr->Face->SetIntPoint(&ip);
            FTr->Face->Transform(ip, fc);

            // u_z is component 2 (dip direction for BP5)
            // For pure strike-slip, u_z should be exactly 0
            real_t uz_avg = 0.5 * (u1[2] + u2[2]);
            real_t uz_jump = u1[2] - u2[2];
            real_t ux_avg = 0.5 * (u1[0] + u2[0]);

            if (std::abs(uz_avg) > 1e-20 || fi < 5)
            {
               mfem::out << "[UZ-FAULT] fi=" << fi
                  << " loc=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                  << " u1=(" << u1[0] << "," << u1[1] << "," << u1[2] << ")"
                  << " u2=(" << u2[0] << "," << u2[1] << "," << u2[2] << ")"
                  << " uz_avg=" << uz_avg
                  << " uz_jump=" << uz_jump
                  << " ux_avg=" << ux_avg
                  << " |uz/ux|=" << (std::abs(ux_avg) > 1e-30 ?
                     std::abs(uz_avg/ux_avg) : 0.0)
                  << std::endl;
            }
         }
      }

      // Fault basis
      const auto &basis = fault_basis_.GetBasis(fi);

      // v55: Sign from FaultBasis sign_flipped (general, not BP5-specific)
      real_t sign = basis.sign_flipped ? 1.0 : -1.0;

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
         // + ProjectTractionToFaultDOFs (nl_q-weighted, sign_flipped)
         // When decomposition diagnostics are requested, uses the decomposed
         // variant to also output stress, correction, and jump residual.
         int nbf = nbf_per_face_;
         bool need_decomp = traction_stress_out || traction_correction_out ||
                            jump_residual_out || coherence_active ||
                            diag_traction_decomp_ ||
                            !diag_first_traction_done_;  // force decomp for tip diagnostic

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
                  delta_u_quad_t(c * nqp_slip + q) = sign * du[c];
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
            for (int j = 0; j < delta_u_quad_t.Size(); j++)
               delta_u_quad_t(j) *= sign;
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
               (traction_stress_out || coherence_active || diag_traction_decomp_)
                  ? &T_stress_quad_dec : nullptr,
               (traction_correction_out || coherence_active || diag_traction_decomp_)
                  ? &T_corr_quad_dec : nullptr,
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

         // [MFEM-TQ] per-QP diagnostic at tip faces (interior), matching Tandem [TND-TQ]
         if (diag_tnd_tq_ && !diag_tnd_tq_done_)
         {
            // Compute face centroid to check if this is a tip face
            const IntegrationPoint &ip_tq =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_tq);
            Vector fc_tq(3);
            FTr->Face->Transform(ip_tq, fc_tq);
            real_t cx_m = fc_tq(0);  // meters

            if (std::abs(cx_m) > 49000.0)  // tip: |x| > 49 km
            {
               int rank_tq = 0;
               if constexpr (IsParallelMesh<MeshType>::value)
               {
#ifdef MFEM_USE_MPI
                  MPI_Comm_rank(mesh_.GetComm(), &rank_tq);
#endif
               }

               FaceVertexKey fkey = MakeFaceKey(face, gvert_tq_);

               // Quad rule (same as used for traction)
               int qo_tq = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
               const IntegrationRule &ir_tq = IntRules.Get(
                  FTr->GetGeometryType(), qo_tq);

               // Element volumes
               real_t vol0_tq = FTr->Elem1->Weight();
               real_t vol1_tq = FTr->Elem2->Weight();

               // Face area
               real_t area_tq = 0.0;
               for (int q = 0; q < nqp_new; q++)
                  area_tq += ir_tq.IntPoint(q).weight * nl_q_vec(q);

               // Per-QP output
               Vector sh1_tq(ndof1), sh2_tq(ndof2);
               for (int q = 0; q < nqp_new; q++)
               {
                  const IntegrationPoint &fip = ir_tq.IntPoint(q);
                  FTr->SetAllIntPoints(&fip);
                  fe1->CalcShape(FTr->GetElement1IntPoint(), sh1_tq);
                  fe2->CalcShape(FTr->GetElement2IntPoint(), sh2_tq);

                  // QP physical coordinates
                  Vector cq(3);
                  FTr->Face->SetIntPoint(&fip);
                  FTr->Face->Transform(fip, cq);

                  // Normal at this QP
                  Vector nor_tq(dim);
                  CalcOrtho(FTr->Jacobian(), nor_tq);
                  real_t nl_tq = nor_tq.Norml2();
                  real_t ny_tq = nor_tq(1) / nl_tq;

                  // Penalty at this QP
                  real_t detJ1_tq = FTr->Elem1->Weight();
                  real_t detJ2_tq = FTr->Elem2->Weight();
                  real_t pen_tq = trac_integ.GetPenalty(
                     *fe1, *fe2, detJ1_tq, detJ2_tq,
                     lambda_val_, mu_val_, nl_tq);

                  // u0_y, u1_y (y = component 1)
                  real_t u0y = 0.0, u1y = 0.0;
                  for (int k = 0; k < ndof1; k++)
                     u0y += sh1_tq(k) * u1_all(1 * ndof1 + k);
                  for (int k = 0; k < ndof2; k++)
                     u1y += sh2_tq(k) * u2_all(1 * ndof2 + k);

                  real_t slip_y = delta_u_quad_t(1 * nqp_new + q);
                  real_t jump_y = u0y - u1y;
                  real_t Ty = T_quad_new(1 * nqp_new + q);

                  mfem::out << std::scientific << std::setprecision(10)
                     << "[MFEM-TQ] r=" << rank_tq
                     << " fct=" << face << " q=" << q
                     << " key=(" << fkey.v[0] << "," << fkey.v[1]
                     << "," << fkey.v[2] << ")"
                     << " cx=" << cq(0) << " cz=" << -cq(2)
                     << " u0_y=" << u0y << " u1_y=" << u1y
                     << " slip_y=" << slip_y
                     << " jump_y=" << jump_y
                     << " Ty=" << Ty
                     << " pen=" << pen_tq
                     << " ny=" << ny_tq
                     << " area=" << area_tq
                     << " vol0=" << vol0_tq << " vol1=" << vol1_tq
                     << " interior"
                     << "\n";
               }
            }
         }

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
            basis_vecs, basis.sign_flipped, trac_local_new, qpd);

         // Store per-DOF shear traction
         int tang_offset = normal_traction ? 1 : 0;
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = fi * nbf_per_face_ + kk;
            traction(2 * dof_idx)     = trac_local_new((tang_offset + 0) * nbf + kk);
            traction(2 * dof_idx + 1) = trac_local_new((tang_offset + 1) * nbf + kk);
         }

         // v58 face-level tip diagnostic: compare mirror left/right tip faces
         if (!diag_first_traction_done_)
         {
            // Compute face centroid
            const IntegrationPoint &ip_diag =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_diag);
            Vector fc_diag(3);
            FTr->Face->Transform(ip_diag, fc_diag);
            real_t cx = fc_diag(0), cz = -fc_diag(2);

            // Only dump for tip faces: |x| > 49km AND depth < 2.5km
            if (std::abs(cx) > 49000.0 && cz < 2500.0)
            {
               // Element volumes
               real_t detJ1_diag = FTr->Elem1->Weight();
               real_t detJ2_diag = FTr->Elem2->Weight();

               // Face area (sum of nl_q * w_q)
               real_t face_area = 0.0;
               for (int q = 0; q < nqp_new; q++)
                  face_area += ir_new.IntPoint(q).weight * nl_q_vec(q);

               // DG penalty (same formula as K assembly)
               real_t nl0 = nl_q_vec(0);
               real_t pen = trac_integ.GetPenalty(
                  *fe1, *fe2, detJ1_diag, detJ2_diag,
                  lambda_val_, mu_val_, nl0);

               // Traction decomposition at quad point 0
               real_t T_stress_q0[3] = {0,0,0}, T_corr_q0[3] = {0,0,0};
               if (T_stress_quad_dec.Size() > 0)
               {
                  for (int c = 0; c < dim; c++)
                  {
                     T_stress_q0[c] = T_stress_quad_dec(c * nqp_new);
                     T_corr_q0[c] = T_corr_quad_dec(c * nqp_new);
                  }
               }

               // Jump at quad point 0: u1-u2-slip
               DenseMatrix shape1_d(fe1->GetDof()), shape2_d(fe2->GetDof());
               {
                  const IntegrationPoint &ip0 = ir_new.IntPoint(0);
                  FTr->SetAllIntPoints(&ip0);
                  Vector s1(fe1->GetDof()), s2(fe2->GetDof());
                  fe1->CalcShape(FTr->GetElement1IntPoint(), s1);
                  fe2->CalcShape(FTr->GetElement2IntPoint(), s2);
                  real_t jump_q0[3] = {0,0,0};
                  for (int c = 0; c < dim; c++)
                  {
                     real_t u1q = 0, u2q = 0;
                     for (int k = 0; k < fe1->GetDof(); k++)
                        u1q += s1(k) * u1_all(c * fe1->GetDof() + k);
                     for (int k = 0; k < fe2->GetDof(); k++)
                        u2q += s2(k) * u2_all(c * fe2->GetDof() + k);
                     real_t fq = delta_u_quad_t(c * nqp_new);
                     jump_q0[c] = u1q - u2q - fq;
                  }

                  int rank = 0;
                  if constexpr (IsParallelMesh<MeshType>::value)
                  {
#ifdef MFEM_USE_MPI
                     MPI_Comm_rank(mesh_.GetComm(), &rank);
#endif
                  }
                  mfem::out << std::scientific << std::setprecision(6)
                     << "[TIP-FACE] rank=" << rank
                     << " fi=" << fi << " interior"
                     << " centroid=(" << cx << "," << fc_diag(1)
                     << "," << fc_diag(2) << ")"
                     << " sign_flipped=" << basis.sign_flipped
                     << " detJ1=" << detJ1_diag
                     << " detJ2=" << detJ2_diag
                     << " face_area=" << face_area
                     << " penalty=" << pen
                     << " nl_q0=" << nl0
                     << "\n    T_stress_q0=(" << T_stress_q0[0]
                     << "," << T_stress_q0[1] << "," << T_stress_q0[2] << ")"
                     << " T_corr_q0=(" << T_corr_q0[0]
                     << "," << T_corr_q0[1] << "," << T_corr_q0[2] << ")"
                     << "\n    jump_q0=(" << jump_q0[0]
                     << "," << jump_q0[1] << "," << jump_q0[2] << ")"
                     << " slip_q0=(" << delta_u_quad_t(0)
                     << "," << delta_u_quad_t(nqp_new)
                     << "," << delta_u_quad_t(2*nqp_new) << ")"
                     << "\n    trac_local=(";
                  for (int kk = 0; kk < nbf; kk++)
                  {
                     int di = fi * nbf_per_face_ + kk;
                     mfem::out << "(" << traction(2*di) << ","
                               << traction(2*di+1) << ")";
                     if (kk < nbf-1) mfem::out << " ";
                  }
                  mfem::out << ")\n";
               }
            }
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
               basis_vecs + tang_offset, basis.sign_flipped, stress_local, qpd);
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
               basis_vecs + tang_offset, basis.sign_flipped, corr_local, qpd);
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
               basis_vecs + tang_offset, basis.sign_flipped, res_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*jump_residual_out)(2 * dof_idx)     = res_local(0 * nbf + kk);
               (*jump_residual_out)(2 * dof_idx + 1) = res_local(1 * nbf + kk);
            }
         }

         // v52: Coherence diagnostic from decomposed quad-point data
         if (coherence_active && T_stress_quad_dec.Size() > 0)
         {
            for (int q = 0; q < nqp_new; q++)
            {
               real_t T_s_q[3], T_c_q[3];
               for (int c = 0; c < dim; c++)
               {
                  T_s_q[c] = T_stress_quad_dec(c * nqp_new + q);
                  T_c_q[c] = T_corr_quad_dec(c * nqp_new + q);
               }
               real_t tau_s[2], tau_c[2];
               fault_basis_.ProjectTraction(fi, T_s_q, tau_s);
               fault_basis_.ProjectTraction(fi, T_c_q, tau_c);

               coh_sum_stress_dip2 += tau_s[0] * tau_s[0];
               coh_sum_stress_strike2 += tau_s[1] * tau_s[1];
               coh_sum_corr_dip2 += tau_c[0] * tau_c[0];
               coh_sum_corr_strike2 += tau_c[1] * tau_c[1];
               coh_max_stress_dip = std::max(coh_max_stress_dip,
                                              std::abs(tau_s[0]));
               coh_max_corr_dip = std::max(coh_max_corr_dip,
                                            std::abs(tau_c[0]));
               coh_n_qp++;
            }

            FTr->SetAllIntPoints(&ip);
            Vector fc(3);
            FTr->Face->SetIntPoint(&ip);
            FTr->Face->Transform(ip, fc);

            // Face-averaged stress/correction for per-face recording
            real_t T_s_avg[3] = {0,0,0}, T_c_avg[3] = {0,0,0};
            real_t sw = 0.0;
            for (int q = 0; q < nqp_new; q++)
            {
               real_t wq = ir_new.IntPoint(q).weight;
               for (int c = 0; c < dim; c++)
               {
                  T_s_avg[c] += wq * T_stress_quad_dec(c * nqp_new + q);
                  T_c_avg[c] += wq * T_corr_quad_dec(c * nqp_new + q);
               }
               sw += wq;
            }
            if (sw > 0.0)
            {
               for (int c = 0; c < dim; c++)
               {
                  T_s_avg[c] /= sw;
                  T_c_avg[c] /= sw;
               }
            }
            real_t tau_s_face[2], tau_c_face[2];
            fault_basis_.ProjectTraction(fi, T_s_avg, tau_s_face);
            fault_basis_.ProjectTraction(fi, T_c_avg, tau_c_face);

            coh_face_data.push_back({fc(0), std::abs(fc(2)),
                                      tau_s_face[0], tau_s_face[1],
                                      tau_c_face[0], tau_c_face[1],
                                      0.0});
         }

         // Traction decomposition diagnostic (face-averaged)
         if (diag_traction_decomp_ && T_stress_quad_dec.Size() > 0)
         {
            real_t T_s_avg[3] = {0,0,0}, T_c_avg[3] = {0,0,0};
            real_t sw = 0.0;
            for (int q = 0; q < nqp_new; q++)
            {
               real_t wq = ir_new.IntPoint(q).weight;
               for (int c = 0; c < dim; c++)
               {
                  T_s_avg[c] += wq * T_stress_quad_dec(c * nqp_new + q);
                  T_c_avg[c] += wq * T_corr_quad_dec(c * nqp_new + q);
               }
               sw += wq;
            }
            if (sw > 0.0)
            {
               for (int c = 0; c < dim; c++) { T_s_avg[c] /= sw; T_c_avg[c] /= sw; }
            }
            real_t tau_stress_local[2], tau_corr_local[2];
            fault_basis_.ProjectTraction(fi, T_s_avg, tau_stress_local);
            fault_basis_.ProjectTraction(fi, T_c_avg, tau_corr_local);

            FTr->SetAllIntPoints(&ip);
            Vector fc(3);
            FTr->Face->SetIntPoint(&ip);
            FTr->Face->Transform(ip, fc);

            mfem::out << "  TRAC_DECOMP interior DOF=" << fi
                      << " x=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                      << " stress_dip=" << tau_stress_local[0]
                      << " stress_strike=" << tau_stress_local[1]
                      << " corr_dip=" << tau_corr_local[0]
                      << " corr_strike=" << tau_corr_local[1]
                      << "\n";
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
               jump_q[c] = (u1q - u2q) - sign * delta_u[c];
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

         // Traction decomposition diagnostic
         if (diag_traction_decomp_)
         {
            real_t tau_stress_local[2], tau_corr_local[2];
            fault_basis_.ProjectTraction(fi, T_stress, tau_stress_local);
            real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
            fault_basis_.ProjectTraction(fi, corr_neg, tau_corr_local);

            FTr->SetAllIntPoints(&ip);
            Vector fc(3);
            FTr->Face->SetIntPoint(&ip);
            FTr->Face->Transform(ip, fc);

            mfem::out << "  TRAC_DECOMP interior DOF=" << fi
                      << " x=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                      << " stress_dip=" << tau_stress_local[0]
                      << " stress_strike=" << tau_stress_local[1]
                      << " corr_dip=" << tau_corr_local[0]
                      << " corr_strike=" << tau_corr_local[1];
            mfem::out << "\n";
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

         // v55: Sign from FaultBasis sign_flipped (general)
         real_t sign = basis.sign_flipped ? 1.0 : -1.0;

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
                                  coherence_active || diag_traction_decomp_;
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
                     delta_u_quad_sh(c * nqp_sh_slip + q) = sign * du[c];
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
               for (int j = 0; j < delta_u_quad_sh.Size(); j++)
                  delta_u_quad_sh(j) *= sign;
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
                  (traction_stress_out || coherence_active || diag_traction_decomp_)
                     ? &T_stress_quad_sh : nullptr,
                  (traction_correction_out || coherence_active || diag_traction_decomp_)
                     ? &T_corr_quad_sh : nullptr,
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

            // [MFEM-TQ] per-QP diagnostic at tip faces (shared)
            if (diag_tnd_tq_ && !diag_tnd_tq_done_)
            {
               const IntegrationPoint &ip_tq =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->Face->SetIntPoint(&ip_tq);
               Vector fc_tq(3);
               FTr->Face->Transform(ip_tq, fc_tq);
               real_t cx_m = fc_tq(0);

               if (std::abs(cx_m) > 49000.0)
               {
                  int lf_sh = mesh_.GetSharedFace(sf);
                  FaceVertexKey fkey = MakeFaceKey(lf_sh, gvert_tq_);

                  int qo_sh = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
                  const IntegrationRule &ir_tq_sh = IntRules.Get(
                     FTr->GetGeometryType(), qo_sh);

                  real_t vol0_sh = FTr->Elem1->Weight();
                  real_t vol1_sh = FTr->Elem2->Weight();

                  real_t area_sh = 0.0;
                  for (int q = 0; q < nqp_sh; q++)
                     area_sh += ir_tq_sh.IntPoint(q).weight * nl_q_sh(q);

                  Vector sh1_tq(ndof1), sh2_tq(ndof2);
                  for (int q = 0; q < nqp_sh; q++)
                  {
                     const IntegrationPoint &fip = ir_tq_sh.IntPoint(q);
                     FTr->SetAllIntPoints(&fip);
                     fe1->CalcShape(FTr->GetElement1IntPoint(), sh1_tq);
                     fe2->CalcShape(FTr->GetElement2IntPoint(), sh2_tq);

                     Vector cq(3);
                     FTr->Face->SetIntPoint(&fip);
                     FTr->Face->Transform(fip, cq);

                     Vector nor_tq(dim);
                     CalcOrtho(FTr->Jacobian(), nor_tq);
                     real_t nl_tq = nor_tq.Norml2();
                     real_t ny_tq = nor_tq(1) / nl_tq;

                     real_t detJ1_sh = FTr->Elem1->Weight();
                     real_t detJ2_sh = FTr->Elem2->Weight();
                     real_t pen_tq = trac_integ_sh.GetPenalty(
                        *fe1, *fe2, detJ1_sh, detJ2_sh,
                        lambda_val_, mu_val_, nl_tq);

                     real_t u0y = 0.0, u1y = 0.0;
                     for (int k = 0; k < ndof1; k++)
                        u0y += sh1_tq(k) * u1_all(1 * ndof1 + k);
                     for (int k = 0; k < ndof2; k++)
                        u1y += sh2_tq(k) * u2_all(1 * ndof2 + k);

                     real_t slip_y = delta_u_quad_sh(1 * nqp_sh + q);
                     real_t jump_y = u0y - u1y;
                     real_t Ty = T_quad_sh(1 * nqp_sh + q);

                     mfem::out << std::scientific << std::setprecision(10)
                        << "[MFEM-TQ] r=" << rank
                        << " fct=" << sf << " q=" << q
                        << " key=(" << fkey.v[0] << "," << fkey.v[1]
                        << "," << fkey.v[2] << ")"
                        << " cx=" << cq(0) << " cz=" << -cq(2)
                        << " u0_y=" << u0y << " u1_y=" << u1y
                        << " slip_y=" << slip_y
                        << " jump_y=" << jump_y
                        << " Ty=" << Ty
                        << " pen=" << pen_tq
                        << " ny=" << ny_tq
                        << " area=" << area_sh
                        << " vol0=" << vol0_sh << " vol1=" << vol1_sh
                        << " shared"
                        << "\n";
                  }
               }
            }

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
               basis_vecs_sh, basis_sh.sign_flipped, trac_local_sh, qpd_sh);

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
                  basis_vecs_sh + tang_off_sh, basis_sh.sign_flipped, stress_local_sh, qpd_sh);
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
                  basis_vecs_sh + tang_off_sh, basis_sh.sign_flipped, corr_local_sh, qpd_sh);
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
                  basis_vecs_sh + tang_off_sh, basis_sh.sign_flipped, res_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*jump_residual_out)(2 * dof_idx)     = res_local_sh(0 * nbf_sh + kk);
                  (*jump_residual_out)(2 * dof_idx + 1) = res_local_sh(1 * nbf_sh + kk);
               }
            }

            // Traction decomposition diagnostic (shared faces)
            if (diag_traction_decomp_ && T_stress_quad_sh.Size() > 0)
            {
               real_t T_s_avg[3] = {0,0,0}, T_c_avg[3] = {0,0,0};
               real_t sw = 0.0;
               for (int q = 0; q < nqp_sh; q++)
               {
                  real_t wq = ir_sh2.IntPoint(q).weight;
                  for (int c = 0; c < dim; c++)
                  {
                     T_s_avg[c] += wq * T_stress_quad_sh(c * nqp_sh + q);
                     T_c_avg[c] += wq * T_corr_quad_sh(c * nqp_sh + q);
                  }
                  sw += wq;
               }
               if (sw > 0.0)
               {
                  for (int c = 0; c < dim; c++) { T_s_avg[c] /= sw; T_c_avg[c] /= sw; }
               }
               real_t tau_stress_local[2], tau_corr_local[2];
               fault_basis_.ProjectTraction(trac_idx, T_s_avg, tau_stress_local);
               fault_basis_.ProjectTraction(trac_idx, T_c_avg, tau_corr_local);

               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->SetIntPoint(&ip);
               FTr->Face->Transform(ip, fc);

               mfem::out << "  TRAC_DECOMP shared DOF=" << trac_idx
                         << " rank=" << rank
                         << " x=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                         << " stress_dip=" << tau_stress_local[0]
                         << " stress_strike=" << tau_stress_local[1]
                         << " corr_dip=" << tau_corr_local[0]
                         << " corr_strike=" << tau_corr_local[1]
                         << "\n";
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
                  jump_q[c] = (u1q - u2q) - sign * delta_u[c];
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

            // Traction decomposition diagnostic (shared faces)
            if (diag_traction_decomp_)
            {
               real_t tau_stress_local[2], tau_corr_local[2];
               fault_basis_.ProjectTraction(trac_idx, T_stress, tau_stress_local);
               real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
               fault_basis_.ProjectTraction(trac_idx, corr_neg, tau_corr_local);

               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->SetIntPoint(&ip);
               FTr->Face->Transform(ip, fc);

               mfem::out << "  TRAC_DECOMP shared DOF=" << trac_idx
                         << " rank=" << rank
                         << " x=(" << fc(0) << "," << fc(1) << "," << fc(2) << ")"
                         << " stress_dip=" << tau_stress_local[0]
                         << " stress_strike=" << tau_stress_local[1]
                         << " corr_dip=" << tau_corr_local[0]
                         << " corr_strike=" << tau_corr_local[1];
               mfem::out << "\n";
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

   // Mark [MFEM-TQ] diagnostic as done after processing all faces
   if (diag_tnd_tq_ && !diag_tnd_tq_done_) { diag_tnd_tq_done_ = true; }

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

   // v49: Diagnostic - dump traction at first evaluation with zero slip
   if (diag_first_traction_ && !diag_first_traction_done_)
   {
      int rank = 0;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         MPI_Comm_rank(mesh_.GetComm(), &rank);
      }

      // Check if this is first stage (all slip ~ 0)
      real_t max_slip = 0;
      for (int i = 0; i < slip_bc.Size(); i++)
      {
         max_slip = std::max(max_slip, std::abs(slip_bc(i)));
      }

      if (max_slip < 1e-20)
      {
         // First evaluation with zero slip - dump traction
         int count = 0;
         for (int i = 0; i < num_fault_dofs_; i++)
         {
            real_t tau_dip = traction(2*i);
            real_t tau_strike = traction(2*i+1);
            real_t tau_mag = std::sqrt(tau_dip*tau_dip + tau_strike*tau_strike);

            // Print DOFs with non-negligible traction (> 1 Pa)
            // or first few DOFs for reference
            if (tau_mag > 1.0 || i < 3)
            {
               mfem::out << "[SEED-TRAC] rank=" << rank << " DOF=" << i
                  << " tau=(" << tau_dip << "," << tau_strike << ")"
                  << " mag=" << tau_mag << "\n";
               count++;
            }
         }
         mfem::out << "[SEED-TRAC] rank=" << rank
            << " total_DOFs_with_tau>1Pa: " << count
            << " / " << num_fault_dofs_ << "\n";
         diag_first_traction_done_ = true;
      }

      // v51: Mark u_z diagnostic as done
      if (diag_uz_fault_ && !diag_uz_fault_done_)
      {
         bool any_slip = false;
         for (int i = 0; i < slip_bc.Size(); i++)
         {
            if (std::abs(slip_bc(i)) > 1e-20) { any_slip = true; break; }
         }
         if (any_slip) { diag_uz_fault_done_ = true; }
      }

      // v51: Mark dip traction diagnostic as done after first complete pass
      if (diag_dip_traction_ && !diag_dip_traction_done_)
      {
         // Check if any face had non-zero slip
         bool any_slip = false;
         for (int i = 0; i < slip_bc.Size(); i++)
         {
            if (std::abs(slip_bc(i)) > 1e-20) { any_slip = true; break; }
         }
         if (any_slip)
         {
            diag_dip_traction_done_ = true;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               int rank;
               MPI_Comm_rank(mesh_.GetComm(), &rank);
               mfem::out << "[DIP-TRAC] rank=" << rank << " diagnostic complete\n";
#endif
            }
         }
      }
   }

   // v52: Print traction coherence summary (MPI-reduced)
   // NOTE: ALL ranks must enter this block when coherence_active is true
   // (coherence_active is globally consistent via MPI_Allreduce above).
   // Do NOT gate on local coh_n_qp — ranks without fault faces have 0
   // quad points but must still participate in MPI collectives.
   if (coherence_active)
   {
      real_t g_stress_dip2 = coh_sum_stress_dip2;
      real_t g_stress_strike2 = coh_sum_stress_strike2;
      real_t g_corr_dip2 = coh_sum_corr_dip2;
      real_t g_corr_strike2 = coh_sum_corr_strike2;
      real_t g_res2 = coh_sum_res2;
      real_t g_max_res = coh_max_res;
      real_t g_max_corr_dip = coh_max_corr_dip;
      real_t g_max_stress_dip = coh_max_stress_dip;
      int g_n_qp = coh_n_qp;
      bool is_root = true;

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&coh_sum_stress_dip2, &g_stress_dip2, 1,
                        MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&coh_sum_stress_strike2, &g_stress_strike2, 1,
                        MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&coh_sum_corr_dip2, &g_corr_dip2, 1,
                        MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&coh_sum_corr_strike2, &g_corr_strike2, 1,
                        MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&coh_sum_res2, &g_res2, 1,
                        MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&coh_max_res, &g_max_res, 1,
                        MPI_DOUBLE, MPI_MAX, mesh_.GetComm());
         MPI_Allreduce(&coh_max_corr_dip, &g_max_corr_dip, 1,
                        MPI_DOUBLE, MPI_MAX, mesh_.GetComm());
         MPI_Allreduce(&coh_max_stress_dip, &g_max_stress_dip, 1,
                        MPI_DOUBLE, MPI_MAX, mesh_.GetComm());
         int local_nqp = coh_n_qp;
         MPI_Allreduce(&local_nqp, &g_n_qp, 1,
                        MPI_INT, MPI_SUM, mesh_.GetComm());
         int rank;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         is_root = (rank == 0);
#endif
      }

      if (is_root && g_n_qp > 0)
      {
         real_t rms_stress_dip = std::sqrt(g_stress_dip2 / g_n_qp);
         real_t rms_stress_strike = std::sqrt(g_stress_strike2 / g_n_qp);
         real_t rms_corr_dip = std::sqrt(g_corr_dip2 / g_n_qp);
         real_t rms_corr_strike = std::sqrt(g_corr_strike2 / g_n_qp);
         real_t rms_res = std::sqrt(g_res2 / g_n_qp);

         mfem::out << "\n[TRACTION-COHERENCE] Solve vs Traction Extraction "
                   << "Diagnostic (v52)\n";
         mfem::out << "  Quad points sampled: " << g_n_qp << "\n";
         mfem::out << "\n  Stress-only traction {sigma.n} (should carry "
                   << "physical signal):\n";
         mfem::out << "    RMS tau_strike(stress) = "
                   << rms_stress_strike << " Pa\n";
         mfem::out << "    RMS tau_dip(stress)    = "
                   << rms_stress_dip << " Pa\n";
         mfem::out << "    max |tau_dip(stress)|  = "
                   << g_max_stress_dip << " Pa\n";
         mfem::out << "    dip/strike ratio (RMS) = "
                   << (rms_stress_strike > 1e-30 ?
                       rms_stress_dip / rms_stress_strike * 100.0 : 0.0)
                   << " %\n";
         mfem::out << "\n  Penalty correction eta*(jump-slip) (should be "
                   << "~0 if solve is exact):\n";
         mfem::out << "    RMS tau_strike(corr)   = "
                   << rms_corr_strike << " Pa\n";
         mfem::out << "    RMS tau_dip(corr)      = "
                   << rms_corr_dip << " Pa\n";
         mfem::out << "    max |tau_dip(corr)|    = "
                   << g_max_corr_dip << " Pa\n";
         mfem::out << "    corr/stress ratio (strike) = "
                   << (rms_stress_strike > 1e-30 ?
                       rms_corr_strike / rms_stress_strike * 100.0 : 0.0)
                   << " %\n";
         mfem::out << "    corr/stress ratio (dip)    = "
                   << (rms_stress_dip > 1e-30 ?
                       rms_corr_dip / rms_stress_dip * 100.0 : 0.0)
                   << " %\n";
         mfem::out << "\n  Solver fault residual |[[u]] - slip|:\n";
         mfem::out << "    RMS residual = " << rms_res << " m\n";
         mfem::out << "    max residual = " << g_max_res << " m\n";
         mfem::out << "\n  Interpretation:\n";
         mfem::out << "    If dip/strike(stress) >> 0: DG solution "
                   << "itself has cross-component coupling\n";
         mfem::out << "    If corr/stress >> 0: penalty amplifies "
                   << "solver residual into spurious traction\n";
         mfem::out << "    If both are small but total dip is large: "
                   << "coherence problem in stress evaluation\n";
      }

      // Per-station nearest-face report
      // BP5 benchmark stations: (x2 [m], x3 [m])
      struct StationDef { const char *name; real_t x2; real_t x3; };
      StationDef stations[] = {
         {"strk-36dp+00", -36e3,  0.0},
         {"strk-16dp+00", -16e3,  0.0},
         {"strk+00dp+00",   0.0,  0.0},
         {"strk+16dp+00",  16e3,  0.0},
         {"strk+36dp+00",  36e3,  0.0},
         {"strk-24dp+10", -24e3, 10e3},
         {"strk-16dp+10", -16e3, 10e3},
         {"strk+00dp+10",   0.0, 10e3},
         {"strk+16dp+10",  16e3, 10e3},
         {"strk+00dp+22",   0.0, 22e3},
      };
      int n_stations = 10;

      // Each rank finds its nearest face to each station
      // Pack: [dist, stress_dip, stress_strike, corr_dip, corr_strike, max_res]
      const int n_fields = 6;
      std::vector<double> local_best(n_stations * n_fields, 1e30);
      for (int s = 0; s < n_stations; s++)
      {
         local_best[s * n_fields + 0] = 1e30;  // distance (sentinel)
      }

      for (size_t f = 0; f < coh_face_data.size(); f++)
      {
         const auto &fd = coh_face_data[f];
         for (int s = 0; s < n_stations; s++)
         {
            real_t dx = fd.x2 - stations[s].x2;
            real_t dz = fd.x3 - stations[s].x3;
            real_t dist = std::sqrt(dx*dx + dz*dz);
            if (dist < local_best[s * n_fields + 0])
            {
               local_best[s * n_fields + 0] = dist;
               local_best[s * n_fields + 1] = fd.stress_dip;
               local_best[s * n_fields + 2] = fd.stress_strike;
               local_best[s * n_fields + 3] = fd.corr_dip;
               local_best[s * n_fields + 4] = fd.corr_strike;
               local_best[s * n_fields + 5] = fd.max_res;
            }
         }
      }

      // MPI reduce: pick the rank with smallest distance for each station
      std::vector<double> global_best = local_best;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         // Use MPI_MINLOC to find rank with smallest distance per station
         // Pack distance + rank into MPI_DOUBLE_INT struct
         struct { double val; int rank; } local_dr[10], global_dr[10];
         int my_rank;
         MPI_Comm_rank(mesh_.GetComm(), &my_rank);
         for (int s = 0; s < n_stations; s++)
         {
            local_dr[s].val = local_best[s * n_fields + 0];
            local_dr[s].rank = my_rank;
         }
         MPI_Allreduce(local_dr, global_dr, n_stations,
                        MPI_DOUBLE_INT, MPI_MINLOC, mesh_.GetComm());

         // Broadcast each winning rank's data
         for (int s = 0; s < n_stations; s++)
         {
            MPI_Bcast(&local_best[s * n_fields], n_fields,
                       MPI_DOUBLE, global_dr[s].rank, mesh_.GetComm());
         }
         global_best = local_best;
#endif
      }

      if (is_root)
      {
         mfem::out << "\n  Per-station decomposition (nearest fault face):\n";
         mfem::out << "  " << std::setw(16) << "Station"
                   << std::setw(10) << "dist(m)"
                   << std::setw(14) << "stress_dip"
                   << std::setw(14) << "stress_strk"
                   << std::setw(12) << "dip/strk%"
                   << std::setw(14) << "corr_dip"
                   << std::setw(14) << "corr_strk"
                   << std::setw(12) << "max_res"
                   << "\n";
         for (int s = 0; s < n_stations; s++)
         {
            double dist = global_best[s * n_fields + 0];
            double sd   = global_best[s * n_fields + 1];
            double ss   = global_best[s * n_fields + 2];
            double cd   = global_best[s * n_fields + 3];
            double cs   = global_best[s * n_fields + 4];
            double mr   = global_best[s * n_fields + 5];
            double ratio = (std::abs(ss) > 1e-30)
                           ? std::abs(sd) / std::abs(ss) * 100.0 : 0.0;
            mfem::out << "  " << std::setw(16) << stations[s].name
                      << std::setw(10) << std::fixed << std::setprecision(0)
                      << dist
                      << std::setw(14) << std::scientific << std::setprecision(3)
                      << sd
                      << std::setw(14) << ss
                      << std::setw(12) << std::fixed << std::setprecision(1)
                      << ratio
                      << std::setw(14) << std::scientific << std::setprecision(3)
                      << cd
                      << std::setw(14) << cs
                      << std::setw(12) << mr
                      << "\n";
         }
         mfem::out << std::defaultfloat;
         mfem::out << "\n[TRACTION-COHERENCE] Done.\n\n";
      }
      diag_traction_coherence_done_ = true;
   }
}

// Convenience type alias
using SerialElasticityDomainOperator = ElasticityDomainOperator<Mesh>;

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ELASTICITY_OPERATOR_HPP
