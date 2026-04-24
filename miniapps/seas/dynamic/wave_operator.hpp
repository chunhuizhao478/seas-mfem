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
#include "precomputed_face_fluxes.hpp"
#include "shared_fault_key.hpp"
#include "../domain/boundary_config.hpp"
#include "../fault/fault_basis.hpp"
#include "../common/seas_types.hpp"
#include "seas_diag_rank.hpp"

#include <algorithm>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <sstream>

namespace mfem
{
namespace seas
{

/// Free-surface boundary-condition flux variant (I-04).
/// - Gamma:   historical gamma-mirror ghost state (default).
/// - Godunov: characteristic Godunov projection (SeisSol parity).
/// Controlled at runtime via WaveOperator::SetFreeSurfaceBCMode.  Default
/// is Gamma so drivers that don't set it reproduce pre-v9.3.0 behaviour.
enum class FreeSurfaceBCMode : int { Gamma = 0, Godunov = 1 };

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

   /// I-04: select which free-surface BC flux variant is dispatched by
   /// ComputeFaceFluxRHS at FaceBC::FreeSurface.  Default Gamma keeps
   /// pre-v9.3.0 byte-identical output; set to Godunov to use the
   /// characteristic projection.
   void SetFreeSurfaceBCMode(FreeSurfaceBCMode m) { free_surface_bc_mode_ = m; }
   FreeSurfaceBCMode GetFreeSurfaceBCMode() const { return free_surface_bc_mode_; }

   /// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
   /// Phase 2a (§6.3, R4-001 + R5-002 FIXES): opt-in dispatch switch for
   /// local non-fault face flux.  Default OFF — drivers that don't set it
   /// reproduce the pre-patch runtime Godunov path bit-identically.
   ///
   /// On the first call with `enable == true`, populates `fault_face_set_`
   /// from `fault_interior_faces_` / `fault_shared_faces_` (R5-002 late
   /// population per §6.3a) and calls `precomputed_face_fluxes_.Init(...)`.
   /// Subsequent calls just flip the flag; `Init` is idempotent via
   /// `IsInitialized()` guard.
   ///
   /// `UsingPrecomputedFaceFluxes()` exposes the current flag value for
   /// test introspection (e.g. `test_precomputed_fluxes_phase2a_switch`).
   ///
   /// `GetPrecomputedFaceFluxes()` is a test-only accessor used by
   /// `test_R5_002_fault_face_set_populated_at_init_time` to verify that
   /// `fault_face_set_` faces are absent from `face_elem_to_entry_`.
   void UsePrecomputedFaceFluxes(bool enable);
   bool UsingPrecomputedFaceFluxes() const { return use_precomputed_face_fluxes_; }
   const PrecomputedFaceFluxes &GetPrecomputedFaceFluxes() const
   { return precomputed_face_fluxes_; }
   const std::set<int> &GetFaultFaceSet() const { return fault_face_set_; }

   /// R-003 (v9.5.0): opt in to Arm 3d topology-based shape-table
   /// precomputation.  Default OFF.  Must be set BEFORE the first call
   /// to `UsePrecomputedFaceFluxes(true)` — the flag is read at Init
   /// time and ignored afterwards (Init is idempotent).  Production
   /// drivers leave it off; the Arm 3d unit test opts in so its
   /// AddInteriorFaceRhsFull / AddBoundaryFaceRhsFull probes have the
   /// shape tables they need.
   void EnableArm3dShapeTables(bool enable) { arm3d_tables_enabled_ = enable; }
   bool Arm3dShapeTablesEnabled() const { return arm3d_tables_enabled_; }

   /// I-06 migration: supply a bulk background state used by all total-Q
   /// BC variants (R-I06-001 absorbing, R-I06-005 free-surface, R-I06-007
   /// PML).  When set:
   ///   - FaceBC::Absorbing      → GodunovFlux::AbsorbingTotal(Q_self, Q_bg)
   ///   - FaceBC::FreeSurface γ  → GodunovFlux::FreeSurfaceTotal(Q_self, Q_bg)
   ///   - FaceBC::FreeSurface G  → GodunovFlux::FreeSurfaceGodunovTotal(Q_self, Q_bg)
   ///   - ApplyPMLDamping        → damps (Q - Q_bg) toward zero
   /// When not set (default), all paths use the pre-migration
   /// fluctuation-Q semantics (Q_ghost = 0, Q damped toward 0).
   ///
   /// R-I06-round3 R-004: the background buffer is OWNED by the wave
   /// operator (NUM_STATE doubles stored inline).  The caller passes a
   /// 9-component array and the values are copied in, so no external
   /// buffer lifetime management is required.  Setting to nullptr
   /// clears the background (reverts to fluctuation semantics).
   /// For TPV102, Q_bg is the uniform pre-stress tensor in global
   /// coordinates.
   void SetAbsorbingBackground(const real_t *Q_bg)
   {
      if (Q_bg == nullptr)
      {
         has_bulk_bg_ = false;
         std::fill(bulk_bg_, bulk_bg_ + NUM_STATE, real_t(0));
      }
      else
      {
         has_bulk_bg_ = true;
         for (int c = 0; c < NUM_STATE; c++) { bulk_bg_[c] = Q_bg[c]; }
      }
   }
   const real_t *GetAbsorbingBackground() const
   { return has_bulk_bg_ ? bulk_bg_ : nullptr; }

   /// Get shared face boundary attributes (for driver's shared fault DOFData collection).
   const std::vector<int> &GetSharedFaceBdrAttr() const { return shared_face_bdr_attr_; }

   /// ADER Phase 1: L2-projected element-local spatial derivative.
   ///
   /// Computes `dQ_dxdir[c,i] = (M_e^{-1} · K_d^e · Q_c)[i]` on every element,
   /// where `K_d^e[i,j] = ∫_e φ_i · ∂_{x_dir} φ_j dV` and `M_e^{-1}` is the
   /// cached element mass inverse.  No inter-element flux coupling — this is
   /// the pure element-local building block for the ADER Cauchy-Kovalevskaya
   /// recursion (Phase 3).
   ///
   /// Quadrature uses the same `IntRules.Get(geom, 2*order_)` rule as
   /// `ComputeVolumeRHS`, so the projection is consistent with the volume
   /// integral.
   ///
   /// @param[in]  dir      Spatial direction in {0, 1, 2}.
   /// @param[in]  Q        State vector of size NUM_STATE * ndof_total_.
   /// @param[out] dQ_dxdir Output spatial derivative, same size as Q.
   ///                       Resized if empty on entry.
   void ApplySpatialDerivative(int dir,
                               const Vector &Q,
                               Vector &dQ_dxdir) const;

   /// ADER Phase 3: time-integrated state via Cauchy-Kovalevskaya recursion.
   ///
   /// Computes `I = ∫_0^{dt} Q(x, t+τ) dτ` element-locally (no inter-element
   /// flux coupling) by Taylor-expanding `Q(t+τ)` and replacing time
   /// derivatives with spatial derivatives via `∂_t Q = -Σ_d A_d ∂_{x_d} Q`:
   ///   D(0)     = Q
   ///   D(k+1)   = L(D(k)) = -Σ_d A_d · ∂_{x_d} D(k)   (element-local)
   ///   I        = Σ_{k=0}^{O-1} (dt^{k+1} / (k+1)!) · D(k)
   ///
   /// Storage cost: three Q-sized buffers (ping-pong D_curr/D_next + dQ_dxd).
   /// PML damping is NOT included in the predictor (the plan defers PML to
   /// the corrector step in a follow-up phase).
   ///
   /// @param[in]  Q      State at t, size NUM_STATE * ndof_total_.
   /// @param[in]  dt     Time step.  dt ≤ 0 yields `I = 0`.
   /// @param[in]  order  ADER order O in {2, 3, 4}.
   /// @param[out] I      Time-integrated state, same size as Q.  Units:
   ///                     [Q] × time.  Resized if empty on entry.
   void ComputeADERTimeIntegrated(const Vector &Q,
                                  real_t dt,
                                  int order,
                                  Vector &I) const;

   /// ADER Phase 4: volume-integral contribution of the ADER corrector.
   ///
   /// Adds `Σ_d ∫_e (∇φ_i)_d · (A_d · I)[c] dV` to `rhs` on every element
   /// (same strong-form DG volume integrand as `ComputeVolumeRHS`, applied
   /// to the time-integrated state `I` rather than instantaneous `Q`).
   ///
   /// Additive contract: `rhs` is NOT zeroed; the caller is responsible.
   /// For `I = Q · dt` (explicit-Euler predictor), this routine is
   /// bit-identical to `dt * ComputeVolumeRHS(Q, rhs)` — the two share a
   /// single integration kernel.
   ///
   /// @param[in]  I    Time-integrated state of size NUM_STATE * ndof_total_.
   /// @param[in,out] rhs  RHS vector of size NUM_STATE * ndof_total_.  Resized if empty.
   void ComputeADERVolumeUpdate(const Vector &I, Vector &rhs) const;

   /// ADER Phase 6: one-step predictor-corrector update of Q.
   ///
   /// Computes `Q_new = Q + M^{-1} [ Vol(I) - FaceFlux(I) ]` where
   ///   I   = ∫_0^{dt} Q(t+τ) dτ  (Phase 3 Cauchy-Kovalevskaya predictor)
   ///   Vol(I)       = ADER volume contribution (Phase 4)
   ///   FaceFlux(I)  = ADER face fluxes — interior, boundary, and fault
   ///                  (Phases 5-6; fault faces use `EvaluateADER[Total]`)
   ///
   /// Replaces a 4-stage RK4 step with a single predictor-corrector sweep.
   /// On a linear (no-friction) problem, `AdvanceADER(Q, dt, 2, Q_new)`
   /// matches the 4-stage RK4 output to O(dt²); O(dt^min(O,4)) for higher
   /// ADER order O.  RK4 path is unchanged (this is a parallel entry point).
   ///
   /// @param[in]  Q      State at t.  Size NUM_STATE * ndof_total_.
   /// @param[in]  dt     Time step.  Must be > 0.
   /// @param[in]  order  ADER order in {2, 3, 4}.
   /// @param[out] Q_new  State at t + dt.  Resized if empty.  Must be
   ///                    DISTINCT from Q (same aliasing constraint as
   ///                    `ComputeADERTimeIntegrated`).
   void AdvanceADER(const Vector &Q, real_t dt, int order,
                    Vector &Q_new) const;

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

   /// R-801 fix: look up the FaultBasis index for an interior fault face by
   /// its mesh face index.  Returns -1 if the face is not an interior fault
   /// face (i.e., not in fault_interior_faces_).
   ///
   /// FaultBasis stores basis data for [interior fault faces, then shared
   /// fault faces]; interior faces occupy positions [0, nfi), so the
   /// returned index is directly usable as `fault_basis_->GetBasis(idx)`.
   int LookupInteriorFaultBasisIndex(int mesh_face_idx) const
   {
      auto it = fault_interior_face_to_basis_idx_.find(mesh_face_idx);
      return (it == fault_interior_face_to_basis_idx_.end()) ? -1 : it->second;
   }

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
   /// DOFData entry.  Under R-701's canonical-fault-frame path (both ranks
   /// run Evaluate in the same pre-sign-flip frame derived from
   /// `FaultBasis`'s `sign_flipped` bit), the two entries are bit-identical
   /// by construction.  This method gathers all shared-fault DOFData across
   /// ranks, matches pairs by face centroid, and asserts bit-equality
   /// (within the relative `tol`).  On mismatch, it calls `MFEM_ABORT`.
   ///
   /// Intended usage: call once from the driver after the first RK4 step
   /// as regression insurance.  Under R-701 a mismatch indicates a new bug
   /// in the canonical-frame reconstruction, the sign_flipped swap, or the
   /// driver's RK4 averaging — it is NO LONGER the v1–v5 symptom.
   /// Cost is O(n_shared_fault_global) communication + memory, one-shot.
   ///
   /// No-op on serial builds.
   ///
   /// @param[in] tol  RELATIVE tolerance for per-field equality (R-502):
   ///                 `|a - b| / max(|a|, |b|, 1.0) > tol` triggers abort.
   ///                 Default 1e-10 tolerates the double-precision noise
   ///                 floor on Pa-scale fields (~1e+8 × ULP ≈ 2e-8).
   ///                 Under R-701's canonical frame this is expected to
   ///                 report `max_rel_diff=0` always; it remains as
   ///                 regression insurance.
   void VerifySharedFaultDOFDataConsistency(real_t tol = 1e-10) const;
   ///@}

#ifdef SEAS_TEST_INTERNAL
   /// Test-only accessors gated by -DSEAS_TEST_INTERNAL.  NEVER defined
   /// in production drivers; used only by unit tests in tests/unit and
   /// tests/parallel to probe internal bookkeeping.  See pepper-bug
   /// unit-test plan 2026-04-22 Phases 1, 2, 4.
   const std::map<int, int> &
   GetFaultInteriorFaceToBasisIdx() const
   { return fault_interior_face_to_basis_idx_; }

   const std::vector<bool> &
   GetSharedFaultElem1OnPlus() const
   { return shared_fault_elem1_on_plus_; }

   const std::map<int, int> &
   GetFaultFaceDofOffset() const
   { return fault_face_dof_offset_; }

   const std::map<int, int> &
   GetSharedFaultDofOffset() const
   { return shared_fault_dof_offset_; }

   /// Test-only mutator: tamper with `shared_fault_elem1_on_plus_[idx]`
   /// for tampering-demonstration assertions (plan §1 acceptance).
   void TamperSharedFaultElem1OnPlus(int idx)
   {
      MFEM_VERIFY(idx >= 0 &&
                  idx < static_cast<int>(shared_fault_elem1_on_plus_.size()),
                  "TamperSharedFaultElem1OnPlus: idx out of range");
      shared_fault_elem1_on_plus_[idx] = !shared_fault_elem1_on_plus_[idx];
   }
#endif

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

   /// I-04: free-surface BC flux dispatch mode.  Defaults to Gamma so
   /// setup-free drivers keep pre-v9.3.0 output.
   FreeSurfaceBCMode free_surface_bc_mode_ = FreeSurfaceBCMode::Gamma;

   /// TPV102 Phase 2a (§6.3): opt-in flag + cached precomputed flux tables.
   /// `precomputed_face_fluxes_` is mutable because the const dispatch
   /// paths (`Mult`, `AdvanceADER` via `ComputeADERFaceFluxRHS`) apply its
   /// matrices to accumulate into `rhs`; `Init` is only called from the
   /// non-const `UsePrecomputedFaceFluxes(true)` public setter.
   /// R4-002 / R5-002 FIX: `fault_face_set_` is populated inside
   /// `UsePrecomputedFaceFluxes(true)` (see §6.3a), NOT in the ctor;
   /// it holds the union of `fault_interior_faces_` and
   /// `fault_shared_faces_` so `PrecomputedFaceFluxes::Init` can skip
   /// fault faces (they go through `FaultFaceFlux`).
   bool use_precomputed_face_fluxes_ = false;
   mutable PrecomputedFaceFluxes precomputed_face_fluxes_;
   std::set<int> fault_face_set_;

   /// R-003 (v9.5.0): opt-in flag for the Arm 3d topology-based
   /// shape-table precomputation.  When true, the first
   /// `UsePrecomputedFaceFluxes(true)` call asks `Init` to populate
   /// `FaceEntry::shape_self` / `shape_nbr` / `w_qp`; when false the
   /// tables are left empty (production default — the Full path is
   /// not wired into dispatch).  Read only at Init time.
   bool arm3d_tables_enabled_ = false;

   /// I-06: bulk background state for total-Q BC dispatch
   /// (R-I06-001 absorbing / R-I06-005 free-surface / R-I06-007 PML).
   /// R-I06-round3 R-004: owned inline buffer (NUM_STATE doubles);
   /// `has_bulk_bg_` gates the total-Q code paths.  When false, all BC
   /// dispatch falls through to the pre-migration fluctuation semantics
   /// (Q_ghost = 0, Q damped toward 0).  `SetAbsorbingBackground(ptr)`
   /// copies values in, removing any external-buffer lifetime concern.
   bool has_bulk_bg_ = false;
   real_t bulk_bg_[NUM_STATE] = {0};
   std::vector<DOFData> *fault_dof_data_ = nullptr;
   std::map<int, int> fault_face_dof_offset_;  ///< face_index → DOFData start index
   std::map<int, int> shared_fault_dof_offset_;  ///< shared_face_index → DOFData start index
   std::vector<int> shared_face_bdr_attr_;  ///< shared face boundary attr (0 = regular interior)
   std::set<int> shared_mesh_face_set_;  ///< mesh face indices that are shared (ParMesh only)
   // R-705 fix: SharedFacePeer struct and shared_face_peer_ member removed.
   // Under R-701 the canonical fault frame (reconstructed from FaultBasis)
   // plus a geometry-based `elem1_on_plus_side` swap flag make Evaluate
   // inputs bit-identical on both ranks without any owner/non-owner MPI
   // exchange.  No peer-rank resolution needed.
   /// For each shared-fault face (indexed by position in fault_shared_faces_),
   /// true iff this rank's local Elem1 sits on the canonical "+" side of
   /// the fault (opposite to where ref_normal points).  The flag is derived
   /// purely from local geometry (element centroid vs face centroid vs
   /// ref_normal), so it is robust to whatever CalcOrtho orientation
   /// convention MFEM uses for shared faces.  Populated in the ctor.
   std::vector<bool> shared_fault_elem1_on_plus_;

   // Canonical fault-face geometry lists (built in constructor).
   // See GetFaultInteriorFaces / GetFaultSharedFaces for layout contract.
   Array<int> fault_interior_faces_;     ///< mesh face indices, 2-sided & non-shared
   Array<int> fault_shared_faces_;       ///< shared-face indices (sf), ParMesh only
   int nbf_per_face_ = 0;                ///< QPs per fault face (set by SetFaultDOFData)

   /// R-801 fix: mesh face index → position in fault_interior_faces_ (=
   /// FaultBasis interior basis index).  Built in the ctor; lets
   /// ComputeFaceFluxRHS reconstruct the BP5 canonical frame for each
   /// interior fault QP (same frame convention as the shared-fault path),
   /// so DOFData.V1/V2/tau1_corr/tau2_corr/slip1/slip2 carry the SAME
   /// physical meaning (component 1 = dip, component 2 = strike) on
   /// every fault QP — interior or shared.  Previously the interior
   /// branch used `GodunovFlux::BuildFrame` (t1 = strike, t2 = up),
   /// creating a convention split between interior and shared fault QPs.
   std::map<int, int> fault_interior_face_to_basis_idx_;

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

   // Test-visibility accessors — exposed for the Arm 1 localization probes
   // (Phase 3 STOP investigation).  Production code uses `Mult` /
   // `AdvanceADER`; these wrappers forward to the private impls so unit
   // tests can probe volume-only / face-only / mass-inverse-only paths
   // without reassembling the call graph.
public:
   void ComputeVolumeRHS_ForTest(const Vector &Q, Vector &rhs) const
   { ComputeVolumeRHS(Q, rhs); }
   void ComputeFaceFluxRHS_ForTest(const Vector &Q, Vector &rhs) const
   { ComputeFaceFluxRHS(Q, rhs); }
   void ApplyMassInverse_ForTest(Vector &rhs) const
   { ApplyMassInverse(rhs); }

private:
   void ComputeVolumeRHS(const Vector &Q, Vector &rhs) const;
   void ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const;
   void ComputeSharedFaceFluxRHS(const Vector &Q, Vector &rhs) const;
   void ApplyMassInverse(Vector &dQdt) const;
   void AssembleElementMassInverse();
   void ApplyPMLDamping(const Vector &Q, Vector &rhs) const;

   /// ADER Phase 6 internal helpers — face flux contribution on a
   /// time-integrated state `I` with step size `dt`.  The fault-face
   /// branches dispatch to `FaultFaceFlux::EvaluateADER[Total]`; all
   /// non-fault branches apply the standard Godunov flux treating `I`
   /// as a time-integrated state, with the background state `bulk_bg_`
   /// (when set) also scaled by `dt` to get the correctly time-
   /// integrated boundary flux.  Additive into `rhs`.
   void ComputeADERFaceFluxRHS(const Vector &I, real_t dt,
                               Vector &rhs) const;
   void ComputeADERSharedFaceFluxRHS(const Vector &I, real_t dt,
                                     Vector &rhs) const;

   enum class FaceBC { Interior, Absorbing, FreeSurface, Fault };
   FaceBC ClassifyBoundaryFace(int bdr_attr) const;
};

// Implementation in wave_operator.inl
#include "wave_operator.inl"

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_WAVE_OPERATOR_HPP
