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

#ifndef MFEM_SEAS_WAVE_OPERATOR_HPP
#define MFEM_SEAS_WAVE_OPERATOR_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "godunov_flux.hpp"
#include "pml_layer.hpp"
#include "fault_face_flux.hpp"
#include "../domain/boundary_config.hpp"
#include "../fault/fault_basis.hpp"
#include "../common/seas_types.hpp"

#include <memory>
#include <vector>
#include <map>
#include <set>

namespace mfem
{
namespace seas
{

/// @brief DG wave operator for the 3D velocity-stress elastic wave equation.
///
/// Solves dQ/dt + A dQ/dx + B dQ/dy + C dQ/dz = 0 using DG with
/// Godunov upwind flux. Inherits only TimeDependentOperator (R-001 fix),
/// NOT DomainOperator.
///
/// Templated on MeshType (Mesh or ParMesh) for serial/parallel support.
/// Uses FESpaceForMesh trait for automatic FE space type resolution.
///
/// The state vector Q has 9 components per DOF:
///   [sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz, v_x, v_y, v_z]
/// stored as Q[c * ndof_total + local_dof] (component-major).
///
/// Reference: Dumbser & Kaser (2006), de la Puente et al. (2009).
template <typename MeshType = Mesh>
class WaveOperator : public TimeDependentOperator
{
   using FESpaceType = FESpaceForMesh<MeshType>;

public:
   /// @brief Construct a WaveOperator on the given mesh.
   WaveOperator(MeshType &mesh, int order,
                real_t lambda, real_t mu, real_t rho,
                const BoundaryConfig &bc);

   ~WaveOperator() override;

   /// @brief Compute dQ/dt = (M^{-1}) * (-Face + Vol).
   void Mult(const Vector &Q, Vector &dQdt) const override;

   /// @name Accessors
   ///@{
   int GetOrder() const { return order_; }
   int GetNDof() const { return ndof_per_el_; }
   int NumElements() const { return ne_; }

   const GodunovFlux &GetFlux() const { return flux_; }
   const FESpaceType &GetFESpace() const { return *fes_; }

   int GetScalarNDof() const { return ndof_total_; }
   real_t ComputeMaxDt(real_t cfl) const;

   const FaultBasis *GetFaultBasis() const { return fault_basis_.get(); }
   int GetNumFaultDOFs() const { return num_fault_dofs_; }

   const DenseMatrix &GetElementMassInverse(int e) const { return elem_mass_inv_[e]; }

   void SetPML(PMLLayer *pml) { pml_layer_ = pml; }
   const PMLLayer *GetPML() const { return pml_layer_; }

   /// Set/get FaultFaceFlux for fault face dispatch (R-002 fix).
   void SetFaultFlux(FaultFaceFlux *ff) { fault_flux_ = ff; }
   FaultFaceFlux *GetFaultFlux() { return fault_flux_; }

   /// Get shared face boundary attributes (for driver's shared fault DOFData collection).
   const std::vector<int> &GetSharedFaceBdrAttr() const { return shared_face_bdr_attr_; }

   /// @name Fault-face geometry (single source of truth, matches BP5 pattern)
   ///
   /// The wave operator owns the canonical list of fault faces.  Both the
   /// physics (flux dispatch, state evolution) and downstream consumers
   /// (driver-level fault_coords / DOFData initialization, ParaView output)
   /// read these getters.  Building the list twice with different filters
   /// is what produced the v1 visualization/coord-mismatch bug — do not
   /// reconstruct in callers.
   ///
   /// Layout: local fault DOFs occupy indices [0, nbf_per_face_ *
   /// fault_interior_faces_.Size()) followed by shared-fault DOFs in
   /// [nbf_per_face_ * fault_interior_faces_.Size(), ...).
   ///@{
   const Array<int> &GetFaultInteriorFaces() const { return fault_interior_faces_; }
   const Array<int> &GetFaultSharedFaces() const { return fault_shared_faces_; }
   int GetNbfPerFace() const { return nbf_per_face_; }
   int GetNumLocalFaultQPs() const
   { return fault_interior_faces_.Size() * nbf_per_face_; }
   int GetNumSharedFaultQPs() const
   { return fault_shared_faces_.Size() * nbf_per_face_; }
   int GetNumTotalFaultQPs() const
   { return GetNumLocalFaultQPs() + GetNumSharedFaultQPs(); }
   ///@}

   /// Set fault DOF data and build face→DOFData index mapping.
   ///
   /// Uses the wave operator's owned fault-face lists.  The data array's
   /// layout must match: local-interior QPs first (in fault_interior_faces_
   /// order), then shared-fault QPs (in fault_shared_faces_ order).
   ///
   /// @param[in] data           Pointer to DOFData array sized to
   ///                           GetNumTotalFaultQPs().
   /// @param[in] nqp_per_face   Number of QPs per fault face.
   void SetFaultDOFData(std::vector<DOFData> *data, int nqp_per_face = 0)
   {
      fault_dof_data_ = data;
      fault_face_dof_offset_.clear();
      shared_fault_dof_offset_.clear();
      if (!data || nqp_per_face <= 0) { return; }

      MFEM_VERIFY(nqp_per_face == nbf_per_face_ || nbf_per_face_ == 0,
                  "nqp_per_face inconsistent with fault-face list setup");
      nbf_per_face_ = nqp_per_face;

      // Local interior fault faces: index i → offset i*nqp_per_face
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         fault_face_dof_offset_[fault_interior_faces_[i]] = i * nqp_per_face;
      }

      // Shared fault faces: placed after all local-interior QPs
      const int shared_base = fault_interior_faces_.Size() * nqp_per_face;
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         shared_fault_dof_offset_[fault_shared_faces_[i]] =
            shared_base + i * nqp_per_face;
      }
   }

   MeshType &GetMesh() { return mesh_; }
   const MeshType &GetMesh() const { return mesh_; }

   /// @brief R-101 fix: verify shared-fault DOFData consistency across ranks.
   ///
   /// For every shared fault QP, both ranks that own it carry an independent
   /// DOFData entry.  R-001's (+,-) canonicalisation is supposed to keep
   /// those two entries bit-identical, but the claim rests on an MFEM
   /// invariant (identical face normals on both ranks) that the v1 fix
   /// report flagged as unverified.  This method gathers all shared-fault
   /// DOFData across ranks, matches pairs by face centroid, and asserts
   /// bit-equality (within `tol`).  On mismatch, it calls `MFEM_ABORT`.
   ///
   /// Intended usage: call once from the driver after the first RK4 step.
   /// Cost is O(n_shared_fault_global) communication + memory, one-shot.
   ///
   /// No-op on serial builds.
   ///
   /// @param[in] tol  Absolute tolerance for per-field equality; default
   ///                 1e-10 matches the double-precision noise floor of
   ///                 one Evaluate call.
   void VerifySharedFaultDOFDataConsistency(real_t tol = 1e-10) const;
   ///@}

private:
   MeshType &mesh_;
   int order_;
   int ndof_per_el_;
   int ne_;
   int ndof_total_;

   std::unique_ptr<L2_FECollection> fec_;
   std::unique_ptr<FESpaceType> fes_;

   GodunovFlux flux_;
   BoundaryConfig bc_;

   DenseMatrix Ax_, Ay_, Az_;
   std::vector<DenseMatrix> elem_mass_inv_;

   std::unique_ptr<FaultBasis> fault_basis_;
   int num_fault_dofs_ = 0;

   real_t h_min_;
   std::vector<int> face_bdr_attr_;

   PMLLayer *pml_layer_ = nullptr;
   FaultFaceFlux *fault_flux_ = nullptr;
   std::vector<DOFData> *fault_dof_data_ = nullptr;
   std::map<int, int> fault_face_dof_offset_;  ///< face_index → DOFData start index
   std::map<int, int> shared_fault_dof_offset_;  ///< shared_face_index → DOFData start index
   std::vector<int> shared_face_bdr_attr_;  ///< shared face boundary attr (0 = regular interior)
   std::set<int> shared_mesh_face_set_;  ///< mesh face indices that are shared (ParMesh only)
   /// Per-shared-face peer-rank record.
   ///
   /// R-001 fix: the rank with the lower ID is the canonical "+" owner of a
   /// shared fault face, so both ranks feed the same `(Q+, Q-)` into
   /// `Evaluate` and update `DOFData` identically.
   ///
   /// R-107 fix: split `resolved` from `peer_rank` so a ctor failure
   /// (GetSharedFaceTransformations returned null, or `fn` lookup couldn't
   /// find the face neighbor) is distinguishable from a successful lookup
   /// at rank 0.  `ComputeSharedFaceFluxRHS` (R-102) aborts if a shared
   /// fault face has `resolved == false`.
   struct SharedFacePeer
   {
      bool resolved = false;
      int  peer_rank = -1;
   };
   std::vector<SharedFacePeer> shared_face_peer_;

   // Canonical fault-face geometry lists (built in constructor).
   // See GetFaultInteriorFaces / GetFaultSharedFaces for layout contract.
   Array<int> fault_interior_faces_;     ///< mesh face indices, 2-sided & non-shared
   Array<int> fault_shared_faces_;       ///< shared-face indices (sf), ParMesh only
   int nbf_per_face_ = 0;                ///< QPs per fault face (set by SetFaultDOFData)

   /// Persistent ghost exchange state (R-001/R-003 fix).
   bool ghost_initialized_ = false;
   /// Cached MPI rank (R-109 fix; 0 on serial builds).  Used by R-001's
   /// shared-fault (+,-) canonicalisation inside ComputeSharedFaceFluxRHS;
   /// caching avoids re-querying pmesh.GetMyRank() on every RK4 stage.
   int my_rank_ = 0;
#ifdef MFEM_USE_MPI
   /// Reusable ParGridFunction for ghost exchange (mutable: used in const Mult).
   mutable std::unique_ptr<ParGridFunction> ghost_gf_;
#endif

   void ComputeVolumeRHS(const Vector &Q, Vector &rhs) const;
   void ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const;
   void ComputeSharedFaceFluxRHS(const Vector &Q, Vector &rhs) const;
   void ApplyMassInverse(Vector &dQdt) const;
   void AssembleElementMassInverse();
   void ApplyPMLDamping(const Vector &Q, Vector &rhs) const;

   enum class FaceBC { Interior, Absorbing, FreeSurface, Fault };
   FaceBC ClassifyBoundaryFace(int bdr_attr) const;
};

// Implementation in wave_operator.inl
#include "wave_operator.inl"

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_WAVE_OPERATOR_HPP
