// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_FAULT_GEOMETRY_HPP
#define MFEM_SEAS_FAULT_GEOMETRY_HPP

#include "mfem.hpp"
#include "../config/bp2_params.hpp"
#include "../config/bp5_params.hpp"
#include "../domain/domain_operator.hpp"
#include "../common/seas_types.hpp"
#include "../common/mpi_context.hpp"
#include <iomanip>

// Forward declarations to avoid pulling FieldProjector / StressField3D
// into every translation unit that uses FaultGeometry.  Phase 6 §5's
// ComputeParams calls FieldProjector::ProjectFaultPreStress which
// is defined in io/field_coefficient.cpp with explicit instantiations
// for Mesh and ParMesh; callers that exercise SAFS-mode must link
// field_coefficient.o.

namespace mfem
{
namespace seas
{
class StressField3D;
class FieldProjector;
} // seas
} // mfem

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Manages fault geometry and depth-dependent parameters.
///
/// This class extracts fault DOF information from the domain operator
/// and precomputes depth-dependent parameters (a(z), eta) for efficient
/// access during time stepping.
///
/// Following Tandem's approach for spatial parameter handling.
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class FaultGeometry
{
public:
   /// @brief Construct fault geometry from domain operator.
   ///
   /// @param domain_op Domain operator providing fault DOF information
   /// @param params BP2 benchmark parameters
   FaultGeometry(DomainOperator<MeshType> &domain_op, const BP2Params &params,
                  MPIContext *mpi_ctx = nullptr)
      : params_(params), mpi_ctx_(mpi_ctx)
   {
      // Use the owned fault view for the friction/state ODE.
      num_fault_dofs_ = domain_op.GetNumOwnedFaultDOFs();

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         num_local_fault_dofs_ = num_fault_dofs_;
         // Compute global DOF count
         if (mpi_ctx_)
         {
            num_global_fault_dofs_ = mpi_ctx_->GlobalSumInt(num_local_fault_dofs_);
            ComputeGatherInfo();
         }
         else
         {
            num_global_fault_dofs_ = num_local_fault_dofs_;
         }
      }
      else
      {
         num_local_fault_dofs_ = num_fault_dofs_;
         num_global_fault_dofs_ = num_fault_dofs_;
      }

      if (num_global_fault_dofs_ == 0)
      {
         MFEM_WARNING("FaultGeometry: No fault DOFs found on any rank");
         return;
      }

      if (num_fault_dofs_ == 0)
      {
         // This rank has no fault DOFs — expected with many ranks.
         return;
      }

      // Get depths from domain operator and restrict to the owned view.
      Vector local_depths;
      domain_op.GetFaultDepths(local_depths);
      domain_op.RestrictToOwnedFault(local_depths, depths_);

      // Compute depth-dependent parameters
      ComputeDepthDependentParams();
   }

   /// @brief Construct fault geometry for 3D (BP5) with spatially varying params.
   ///
   /// @param domain_op Domain operator providing 2D fault coordinates
   /// @param params BP5 benchmark parameters (2D spatially varying a, L, etc.)
   FaultGeometry(DomainOperator<MeshType> &domain_op, const BP5Params &params,
                  MPIContext *mpi_ctx = nullptr)
      : bp5_params_(params), mpi_ctx_(mpi_ctx), is_bp5_(true)
   {
      num_fault_dofs_ = domain_op.GetNumOwnedFaultDOFs();
      nbf_per_face_ = domain_op.GetNbfPerFace();
      num_fault_faces_ = (nbf_per_face_ > 0) ? num_fault_dofs_ / nbf_per_face_ : 0;
      num_local_fault_dofs_ = num_fault_dofs_;
      num_global_fault_dofs_ = num_fault_dofs_;

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         if (mpi_ctx_)
         {
            num_global_fault_dofs_ = mpi_ctx_->GlobalSumInt(num_local_fault_dofs_);
            ComputeGatherInfo();
         }
      }

      if (num_fault_dofs_ == 0) { return; }

      // Get 2D fault coordinates on the owned fault view.
      Vector local_x2, local_x3;
      domain_op.GetFaultCoords2D(local_x2, local_x3);
      domain_op.RestrictToOwnedFault(local_x2, coords_x2_);
      domain_op.RestrictToOwnedFault(local_x3, coords_x3_);

      // Also store depths for compatibility
      depths_.SetSize(num_fault_dofs_);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         depths_(i) = coords_x3_(i);
      }

      // Phase 6.A — per-DOF global 3-D coordinates and (n, t1, t2) basis.
      // Write-once at init; read only by SAFS-mode consumers (Phase 6 §4-§7).
      // BP5/BP2/TPV102 code paths do not read these members, so the BP5
      // bit-exact contract is preserved.
      ComputePerDOFCoordsAndBasis_(domain_op);

      // Precompute per-DOF parameters using BP5 2D functions
      ComputeBP5Params();
   }

   /// @brief Phase 5a of spatial_dynamic_rupture_plan.md (rev-3): NEW BP5
   /// ctor overload that accepts pre-built per-DOF arrays directly,
   /// avoiding the need to instantiate a throw-away
   /// `ElasticityDomainOperator` (R-110 Option A).
   ///
   /// Lets the SAFS dynamic-rupture driver build a `FaultGeometry` from
   /// the data it already walked out of `WaveOperator` (interior +
   /// shared fault face arrays) without paying for a MUMPS_BLR
   /// factorisation just to read 2D fault coordinates.
   ///
   /// - `dof_coords_3d` is the OWNED-fault view, interleaved
   ///   `[x_0, y_0, z_0, x_1, y_1, z_1, ...]` of length `3 * N`.
   /// - `dof_basis` is a `9 x N` matrix; columns `3*i .. 3*i+2` hold
   ///   the (n, t1, t2) basis at fault DOF i (BP5 / Tandem convention:
   ///   t1 = dip, t2 = strike).
   /// - `dof_to_elem(i)` is the bulk-element index owning fault DOF i;
   ///   consumed by `dynamic/spatial_setup.hpp::InitializeFaultDOFs_-
   ///   Spatial` to query the MaterialField at the right element.
   /// - `dof_ips` is OPTIONAL.  When non-empty its size MUST equal `N`
   ///   and each entry is the reference `IntegrationPoint` at fault
   ///   DOF i; populates `fault_dof_ip_` so IP-aware downstream
   ///   helpers can evaluate Coefficients at the actual fault QP.
   ///   When empty, `fault_dof_ip_` stays empty and IP-aware readers
   ///   abort with a clear message.
   ///
   /// `coords_x2_` / `coords_x3_` are NaN-filled (R-310): the legacy
   /// BP5 (along-strike, along-dip) convention does NOT apply to SAFS
   /// curvilinear faults, so feeding world (x, z) into
   /// `bp5_params_.a_of_x2_x3` would silently produce garbage.  All
   /// BP5 analytic per-DOF arrays (`a_values_`, `dc_values_`,
   /// `V_init_vec_`, BP5-side `tau_pre_`, `eta_values_`) are also
   /// NaN-filled so any premature read from a non-SAFS path fails
   /// loud rather than silently consuming out-of-range values.
   ///
   /// `depths_(i)` is populated from world z (the value the absorbing
   /// / free-surface BC dispatch expects, R-404).  For a planar BP5
   /// fault aligned with the z-axis the legacy ctor's
   /// `depths_(i) = coords_x3_(i)` (along-dip x3) happens to equal
   /// world z; for SAFS curvilinear faults this IS world depth and
   /// callers that previously assumed BP5-x3 semantics see different
   /// values.
   ///
   /// **The legacy `(DomainOperator&, BP5Params&, MPIContext*)` ctor
   /// above is preserved verbatim** (§TPV Benchmark Isolation gate
   /// `T-FAULTGEO-LEGACY-CTOR-BYTE-IDENTICAL`).
   FaultGeometry(const BP5Params &params,
                 const Vector &dof_coords_3d,
                 const DenseMatrix &dof_basis,
                 const Array<int> &dof_to_elem,
                 int nbf_per_face,
                 MPIContext *mpi_ctx = nullptr,
                 const std::vector<mfem::IntegrationPoint> &dof_ips =
                    std::vector<mfem::IntegrationPoint>())
      : bp5_params_(params), mpi_ctx_(mpi_ctx), is_bp5_(true)
   {
      MFEM_VERIFY(nbf_per_face >= 1,
                  "FaultGeometry(new BP5 ctor): nbf_per_face must be "
                  ">= 1; got " << nbf_per_face);
      const int N = dof_coords_3d.Size() / 3;
      MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
                  "FaultGeometry(new BP5 ctor): dof_coords_3d.Size() ("
                  << dof_coords_3d.Size() << ") must be 3 * N");
      MFEM_VERIFY(dof_basis.Height() == 9 && dof_basis.Width() == N,
                  "FaultGeometry(new BP5 ctor): dof_basis shape ("
                  << dof_basis.Height() << ", " << dof_basis.Width()
                  << ") must be (9, " << N << ")");
      MFEM_VERIFY(dof_to_elem.Size() == N,
                  "FaultGeometry(new BP5 ctor): dof_to_elem.Size() ("
                  << dof_to_elem.Size() << ") must equal N (" << N << ")");
      MFEM_VERIFY(dof_ips.empty()
                  || static_cast<int>(dof_ips.size()) == N,
                  "FaultGeometry(new BP5 ctor): dof_ips.size() ("
                  << dof_ips.size() << ") must be either 0 or N (" << N
                  << ")");

      num_fault_dofs_         = N;
      nbf_per_face_           = nbf_per_face;
      num_fault_faces_        = (nbf_per_face_ > 0) ? N / nbf_per_face_ : 0;
      num_local_fault_dofs_   = N;
      num_global_fault_dofs_  = N;

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         if (mpi_ctx_)
         {
            num_global_fault_dofs_ =
               mpi_ctx_->GlobalSumInt(num_local_fault_dofs_);
            ComputeGatherInfo();
         }
      }

      dof_to_elem_new_ctor_ = dof_to_elem;
      if (!dof_ips.empty()) { fault_dof_ip_ = dof_ips; }

      if (N == 0) { return; }

      dof_coords_3d_ = dof_coords_3d;
      dof_basis_     = dof_basis;

      constexpr real_t k_nan = std::numeric_limits<real_t>::quiet_NaN();
      coords_x2_.SetSize(N);
      coords_x3_.SetSize(N);
      depths_.SetSize(N);
      for (int i = 0; i < N; ++i)
      {
         coords_x2_(i) = k_nan;
         coords_x3_(i) = k_nan;
         depths_(i)    = dof_coords_3d(3 * i + 2);
      }

      // Eagerly size + NaN-fill the BP5 analytic per-DOF arrays
      // (R-402): any premature reader (Print, IsVelocityWeakening,
      // GetAValues, GetEtaValues, GetTauPre, GetVInit) sees a
      // consistent NaN-loud state instead of an out-of-bounds crash.
      a_values_.SetSize(N);
      eta_values_.SetSize(N);
      dc_values_.SetSize(N);
      tau_pre_.SetSize(2 * N);
      V_init_vec_.SetSize(2 * N);
      for (int i = 0; i < N; ++i)
      {
         a_values_(i)            = k_nan;
         eta_values_(i)          = k_nan;
         dc_values_(i)           = k_nan;
         tau_pre_(2 * i + 0)     = k_nan;
         tau_pre_(2 * i + 1)     = k_nan;
         V_init_vec_(2 * i + 0)  = k_nan;
         V_init_vec_(2 * i + 1)  = k_nan;
      }

      // R-304: count degenerate normal / tangent1 columns so the
      // SAFS R-001 pre-flight guard in ApplyCsmStressSidecar can
      // gate on `NumZeroNormalFallbacks() == 0`.
      num_dof_basis_fallbacks_   = 0;
      num_zero_normal_fallbacks_ = 0;
      num_t1_fallbacks_          = 0;
      constexpr real_t k_norm_tol = static_cast<real_t>(1e-10);
      for (int i = 0; i < N; ++i)
      {
         const real_t nx  = dof_basis(0, i);
         const real_t ny  = dof_basis(1, i);
         const real_t nz  = dof_basis(2, i);
         const real_t t1x = dof_basis(3, i);
         const real_t t1y = dof_basis(4, i);
         const real_t t1z = dof_basis(5, i);
         const real_t n_norm  = std::sqrt(nx*nx + ny*ny + nz*nz);
         const real_t t1_norm = std::sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
         if (n_norm < k_norm_tol)
         {
            ++num_zero_normal_fallbacks_;
            ++num_dof_basis_fallbacks_;
         }
         if (t1_norm < k_norm_tol) { ++num_t1_fallbacks_; }
      }
   }

   /// @brief Number of fault DOFs (local in parallel, total in serial).
   int NumFaultDOFs() const { return num_fault_dofs_; }

   /// @brief Number of local fault DOFs (same as NumFaultDOFs).
   int NumLocalFaultDOFs() const { return num_local_fault_dofs_; }

   /// @brief Number of global fault DOFs (sum across all ranks).
   int NumGlobalFaultDOFs() const { return num_global_fault_dofs_; }

   /// @brief Gather local fault data to root rank (rank 0).
   ///
   /// In serial mode, simply copies local_data to global_data.
   /// In parallel mode, uses MPI_Gatherv to collect data on root.
   ///
   /// @param[in] local_data Local fault data [NumLocalFaultDOFs()]
   /// @param[out] global_data Global fault data [NumGlobalFaultDOFs()] (valid on root only)
   void GatherToRoot(const Vector &local_data, Vector &global_data) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef SEAS_USE_MPI
         if (!mpi_ctx_)
         {
            global_data = local_data;
            return;
         }
         if (mpi_ctx_->IsRoot())
         {
            global_data.SetSize(num_global_fault_dofs_);
         }
         MPI_Gatherv(local_data.GetData(), num_local_fault_dofs_, MPI_DOUBLE,
                      mpi_ctx_->IsRoot() ? global_data.GetData() : nullptr,
                      recv_counts_.data(), recv_displs_.data(), MPI_DOUBLE,
                      0, mpi_ctx_->GetComm());
#endif
      }
      else
      {
         global_data = local_data;
      }
   }

   /// @brief Gather and deduplicate fault data on root.
   ///
   /// In parallel DG, shared fault faces at partition boundaries produce
   /// duplicate DOFs. This method gathers data and depths, then on root
   /// sorts by depth, merges duplicates (averaging their field values),
   /// and returns deduplicated results.
   ///
   /// @param[in] local_data Local fault data [NumLocalFaultDOFs()]
   /// @param[out] dedup_data Deduplicated data sorted by depth (root only)
   /// @param[out] dedup_depths Deduplicated depths sorted (root only)
   void GatherToRootDedup(const Vector &local_data,
                           Vector &dedup_data,
                           Vector &dedup_depths) const
   {
      Vector raw_data, raw_depths;
      GatherToRoot(local_data, raw_data);
      GatherToRoot(depths_, raw_depths);

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         if (mpi_ctx_ && mpi_ctx_->IsRoot())
         {
            DeduplicateByDepth(raw_depths, raw_data, dedup_depths, dedup_data);
         }
      }
      else
      {
         dedup_data = raw_data;
         dedup_depths = raw_depths;
      }
   }

   /// @brief Gather multiple fields and depths, deduplicate on root.
   ///
   /// Convenience method for gathering several fault fields at once.
   /// All fields are deduplicated using the same depth-based merging.
   void GatherFieldsToRootDedup(
      const std::vector<const Vector*> &local_fields,
      std::vector<Vector> &dedup_fields,
      Vector &dedup_depths) const
   {
      // Gather depths
      Vector raw_depths;
      GatherToRoot(depths_, raw_depths);

      // Gather all fields
      int nf = static_cast<int>(local_fields.size());
      std::vector<Vector> raw_fields(nf);
      for (int f = 0; f < nf; f++)
      {
         GatherToRoot(*local_fields[f], raw_fields[f]);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         if (mpi_ctx_ && mpi_ctx_->IsRoot())
         {
            DeduplicateMultipleByDepth(raw_depths, raw_fields,
                                       dedup_depths, dedup_fields);
         }
         else
         {
            dedup_fields = raw_fields;
            dedup_depths = raw_depths;
         }
      }
      else
      {
         dedup_fields = raw_fields;
         dedup_depths = raw_depths;
      }
   }

   /// @brief Get the MPI context (may be nullptr in serial).
   MPIContext *GetMPIContext() const { return mpi_ctx_; }

   /// @brief Get depths at each fault DOF (z coordinate, negative below surface).
   const Vector &GetDepths() const { return depths_; }

   /// @brief Get rate-state parameter a at each fault DOF (precomputed from depth).
   const Vector &GetAValues() const { return a_values_; }

   /// @brief Get radiation damping η at each fault DOF.
   const Vector &GetEtaValues() const { return eta_values_; }

   /// @brief Get the BP2 parameters.
   const BP2Params &GetParams() const { return params_; }

   /// @brief Get the BP5 parameters (only valid if constructed with BP5Params).
   const BP5Params &GetBP5Params() const { return bp5_params_; }

   /// @brief Whether this was constructed for BP5 (3D, spatially varying).
   bool IsBP5() const { return is_bp5_; }

   /// @brief Get critical slip distance (Dc/L) at each fault DOF.
   const Vector &GetDcValues() const { return dc_values_; }

   /// @brief Get pre-stress vector at fault DOFs [2*NumFaultDOFs].
   /// Layout: [tau_dip_0, tau_strike_0, tau_dip_1, tau_strike_1, ...]
   const Vector &GetTauPre() const { return tau_pre_; }

   /// @brief Get initial velocity at fault DOFs [2*NumFaultDOFs].
   /// Layout: [V_dip_0, V_strike_0, V_dip_1, V_strike_1, ...]
   const Vector &GetVInit() const { return V_init_vec_; }

   /// @brief Get 2D fault coordinates.
   const Vector &GetCoordsX2() const { return coords_x2_; }
   const Vector &GetCoordsX3() const { return coords_x3_; }

   /// @brief Get per-DOF global 3-D (x, y, z) coordinates [3 * NumFaultDOFs].
   ///
   /// Interleaved layout: entry `3*i + d` is the d-th coordinate of fault
   /// DOF `i`. Populated by the BP5/3-D constructor only; for BP2/antiplane
   /// the returned vector is empty.
   ///
   /// Phase 6.A — used by FieldProjector::ProjectFaultPreStress and
   /// FaultGeometry::ComputeParams for sidecar lookups.
   const Vector &fault_dof_coords_3d() const { return dof_coords_3d_; }

   /// @brief Get per-DOF orthonormal fault basis [9 x NumFaultDOFs].
   ///
   /// Column `i` is `[n_i; t1_i; t2_i]` (BP5 / FaultBasis Tandem
   /// convention: t1 = dip, t2 = strike). The sign-flip is baked into
   /// each vector — callers use them directly. Re-orthonormalised via
   /// Gram-Schmidt at init (Phase 2 `basis_to_node` convention).
   ///
   /// Empty in BP2/antiplane paths.
   const DenseMatrix &fault_dof_basis() const { return dof_basis_; }

   /// @brief Number of fault DOFs whose basis was replaced by a fallback —
   ///        sum of zero-normal and t1-degenerate cases.
   ///
   /// Reported once at init from the BP5 ctor. Always 0 for planar faults
   /// with healthy fault basis on every owned DOF; positive when either:
   ///   (a) the operator's GetFaultDOFBasis returned a zero-length normal
   ///       for some DOF (e.g. boundary or non-fault DOFs in some
   ///       fixtures) — counted in `NumZeroNormalFallbacks()`; the basis
   ///       slot is zeroed and SAFS projection at that DOF would yield
   ///       zero pre-stress and zero sigma_n, so SAFS mode must refuse
   ///       to enable in that case (see `RateStateFaultOperator::SetSAFSMode`).
   ///   (b) the projected dip t1 dropped below 1e-12 and a fallback was
   ///       derived from the reference up-vector — counted in
   ///       `NumT1Fallbacks()`; the slot has a valid (but FaultBasis-
   ///       sign-flip-undefined) basis.  SAFS-mode users may proceed
   ///       with a runtime warning.
   int NumDOFBasisFallbacks() const { return num_dof_basis_fallbacks_; }

   /// @brief Of `NumDOFBasisFallbacks()`, how many were zero-normal
   ///        (basis slot zeroed; SAFS would project to 0).
   int NumZeroNormalFallbacks() const { return num_zero_normal_fallbacks_; }

   /// @brief Of `NumDOFBasisFallbacks()`, how many were t1-degeneracy
   ///        (basis slot valid but sign-flip-undefined).
   int NumT1Fallbacks() const { return num_t1_fallbacks_; }

   /// @brief Phase 5a (rev-3) — size of the per-DOF reference
   /// `IntegrationPoint` cache populated by the new BP5 ctor.
   ///
   /// Returns 0 when this `FaultGeometry` was built by the LEGACY BP5
   /// ctor `(DomainOperator&, BP5Params&, MPIContext*)` (which never
   /// touches `fault_dof_ip_`).  Returns `N` (= NumFaultDOFs) when
   /// built by the new BP5 ctor with a non-empty `dof_ips` argument.
   ///
   /// Use this to gate IP-aware downstream readers (e.g. the IP-aware
   /// overload of `InitializeFaultDOFs_Spatial` in
   /// `dynamic/spatial_setup.hpp`); when 0, fall back to the
   /// centroid-based overload.
   int fault_dof_ip_size() const
   { return static_cast<int>(fault_dof_ip_.size()); }

   /// @brief Phase 5a — accessor for the per-DOF reference
   /// `IntegrationPoint` cache.  Aborts if the cache is empty
   /// (legacy ctor was used or new ctor with `dof_ips.empty()`).
   const mfem::IntegrationPoint &fault_dof_ip(int i) const
   {
      MFEM_ASSERT(!fault_dof_ip_.empty(),
                  "FaultGeometry::fault_dof_ip(): cache is empty.  "
                  "This FaultGeometry was built by the LEGACY BP5 "
                  "ctor (or the NEW ctor with an empty dof_ips "
                  "argument); IP-aware readers must gate on "
                  "fault_dof_ip_size() > 0 before reading.");
      MFEM_ASSERT(i >= 0 && i < static_cast<int>(fault_dof_ip_.size()),
                  "FaultGeometry::fault_dof_ip(" << i
                  << "): out of range.");
      return fault_dof_ip_[i];
   }

   /// @brief Phase 5a — accessor for the per-DOF bulk-element ownership
   /// cache populated by the new BP5 ctor.  Empty when built by the
   /// legacy ctor.  Consumed by
   /// `dynamic/spatial_setup.hpp::InitializeFaultDOFs_Spatial` to
   /// look up `MaterialField::EvalAt` at the right element per DOF.
   const Array<int> &dof_to_elem_new_ctor() const
   { return dof_to_elem_new_ctor_; }

   /// @brief Phase 5a — size of the per-DOF bulk-element ownership cache
   /// populated by the new BP5 ctor.  Returns 0 when built by the
   /// legacy ctor.
   int dof_to_elem_new_ctor_size() const
   { return dof_to_elem_new_ctor_.Size(); }

   /// @brief Phase 5a — element index for fault DOF `i` (new BP5 ctor
   /// path).  Aborts if the cache is empty (legacy ctor was used).
   int dof_to_elem_new_ctor(int i) const
   {
      MFEM_ASSERT(dof_to_elem_new_ctor_.Size() > 0,
                  "FaultGeometry::dof_to_elem_new_ctor(i): cache is "
                  "empty.  This FaultGeometry was built by the LEGACY "
                  "BP5 ctor; readers must gate on "
                  "dof_to_elem_new_ctor_size() > 0 before reading.");
      MFEM_ASSERT(i >= 0 && i < dof_to_elem_new_ctor_.Size(),
                  "FaultGeometry::dof_to_elem_new_ctor(" << i
                  << "): out of range.");
      return dof_to_elem_new_ctor_[i];
   }

   /// @brief Get per-DOF normal stress [NumFaultDOFs].
   ///
   /// Populated by `ComputeParams` only — empty in the standard BP5
   /// path (which uses the scalar `bp5_params_.sigma_n`).
   const Vector &sigma_n_per_dof() const { return sigma_n_per_dof_; }

   /// @brief Whether ComputeParams has been invoked successfully.
   bool HasParams() const { return params_computed_; }

   /// @brief Phase 6 §5 — pre-stress initialisation from a six-component
   /// stress sidecar.  Used project-wide (SAFS, TPV102/104/205, TPV31)
   /// via the templated overload below; this overload handles the
   /// `StressField3D` sidecar input.  Renamed from `ComputeSAFSParams`
   /// to `ComputeParams` for general use (per `tpv102_tpv104_review.md`
   /// R-001 directive: "no special case for the TPV sign convention —
   /// it IS the project-wide convention").
   ///
   /// Parallel slot to ComputeBP5Params: keeps the analytic spatial
   /// `a(x2, x3)`, `Dc(x2, x3)`, `V_init(x2, x3)` and `eta` from
   /// `bp5_params_`, but replaces the analytic `tau0_vec(x2, x3)` and
   /// scalar `bp5_params_.sigma_n` with sidecar-sourced per-DOF values.
   ///
   /// On return, `tau_pre_` and `sigma_n_per_dof_` are populated; the
   /// remaining BP5-state arrays are unchanged from `ComputeBP5Params`.
   /// `params_computed_` is set to `true` so consumers can branch on
   /// it via `HasParams()`.
   ///
   /// Sign convention (CLAUDE.md project-wide, matches native
   /// TPV102/104/205 drivers): positive `tau_pre(2*i+1)` represents
   /// right-lateral driving stress in the +strike direction.  See
   /// `FieldProjector::ProjectFaultPreStress` for the projection
   /// implementation.
   ///
   /// @param field             Six-component sidecar reader.
   /// @param P_p_pa            Constant pore-pressure offset [Pa].
   /// @param P_p_grad_pa_per_m Depth gradient of pore pressure
   ///                           [Pa / m]; effective P_p at z is
   ///                           P_p_pa + grad * max(0, -z).
   /// @param min_sigma_n_pa    Optional Pa-valued floor on the
   ///                           effective normal stress; default 0
   ///                           means no clamp.
   void ComputeParams(const StressField3D& field,
                      real_t P_p_pa = 0.0,
                      real_t P_p_grad_pa_per_m = 0.0,
                      real_t min_sigma_n_pa = 0.0);

   /// @brief Templated overload for any StressSource3D-conformant type
   /// (Phase 3b of spatial_dynamic_rupture_plan.md rev-3).
   ///
   /// Any class `S` exposing
   ///   `mfem::DenseMatrix S::Evaluate(real_t x, real_t y, real_t z) const`
   /// satisfies the concept and is accepted by this template.  Overload
   /// resolution always selects the non-templated
   /// `ComputeParams(const StressField3D&, ...)` overload above for
   /// `StressField3D` arguments (a non-template wins by C++
   /// overload-ranking rules), so the BP5 byte-exact contract is
   /// preserved — this template fires for *other* sources, e.g.
   /// `mfem::seas::spatial::ConstantTensorStressSource` (TPV102/104,
   /// SAFS) and `ConstantTensorWithPatchesStressSource` (TPV205).
   ///
   /// The body lives in `fault/fault_geometry_safs_templated.inl`.
   /// Callers that need the templated overload must include both
   /// `fault_geometry.hpp` AND `fault_geometry_safs_templated.inl`
   /// in the same translation unit so the compiler can instantiate.
   /// BP5 / TPV callers that only consume the non-templated path
   /// continue to include `fault_geometry.hpp` (+ the existing
   /// `fault_geometry_safs.inl` non-template body) and never
   /// instantiate this template.
   template <typename StressSource>
   void ComputeParams(const StressSource& source,
                      real_t P_p_pa = 0.0,
                      real_t P_p_grad_pa_per_m = 0.0,
                      real_t min_sigma_n_pa = 0.0);

   /// @brief Find the DOF index closest to a target depth.
   ///
   /// @param target_depth Target depth (z coordinate, negative for below surface)
   /// @return Index of the closest DOF
   int FindNearestDOF(real_t target_depth) const
   {
      if (num_fault_dofs_ == 0)
      {
         return -1;
      }

      int closest_idx = 0;
      real_t min_dist = std::abs(depths_(0) - target_depth);

      for (int i = 1; i < num_fault_dofs_; i++)
      {
         real_t dist = std::abs(depths_(i) - target_depth);
         if (dist < min_dist)
         {
            min_dist = dist;
            closest_idx = i;
         }
      }

      return closest_idx;
   }

   /// @brief Check if a DOF is in the velocity-weakening zone.
   ///
   /// Uses precomputed a_values_ which are correct for both BP2 and BP5.
   ///
   /// @param dof_idx DOF index
   /// @return True if a(dof_idx) < b (velocity-weakening)
   bool IsVelocityWeakening(int dof_idx) const
   {
      real_t b_val = is_bp5_ ? bp5_params_.b : params_.b;
      return a_values_(dof_idx) < b_val;
   }

   /// @brief Get the VW/VS transition depth (top of transition zone).
   ///
   /// Returns the depth H where the transition from VW to VS begins.
   /// Only valid for BP2 (1D depth profile). For BP5, the VW zone is 2D.
   real_t GetVWDepth() const
   {
      MFEM_VERIFY(!is_bp5_,
                   "GetVWDepth() not applicable for BP5 (2D VW zone)");
      return -params_.H;
   }

   /// @brief Get the full VS depth (bottom of transition zone).
   ///
   /// Returns the depth H+h where fully VS behavior begins.
   /// Only valid for BP2 (1D depth profile). For BP5, the VS zone is 2D.
   real_t GetVSDepth() const
   {
      MFEM_VERIFY(!is_bp5_,
                   "GetVSDepth() not applicable for BP5 (2D VW zone)");
      return -(params_.H + params_.h);
   }

   /// @brief Get indices of DOFs in the velocity-weakening zone.
   void GetVWDOFs(Array<int> &vw_dofs) const
   {
      vw_dofs.SetSize(0);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         if (IsVelocityWeakening(i))
         {
            vw_dofs.Append(i);
         }
      }
   }

   /// @brief Get indices of DOFs in the velocity-strengthening zone.
   void GetVSDOFs(Array<int> &vs_dofs) const
   {
      vs_dofs.SetSize(0);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         if (!IsVelocityWeakening(i))
         {
            vs_dofs.Append(i);
         }
      }
   }

   /// @brief Print fault geometry information.
   void Print(std::ostream &os = mfem::out) const
   {
      os << "Fault Geometry:\n";
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         os << "  Local fault DOFs: " << num_local_fault_dofs_ << "\n";
         os << "  Global fault DOFs: " << num_global_fault_dofs_ << "\n";
      }
      else
      {
         os << "  Number of DOFs: " << num_fault_dofs_ << "\n";
      }

      if (num_fault_dofs_ > 0)
      {
         // Local statistics (rank 0 only — global reduction would deadlock
         // since Print() is called only on root)
         real_t z_min = depths_.Min();
         real_t z_max = depths_.Max();

         real_t b_val = is_bp5_ ? bp5_params_.b : params_.b;
         int vw_count = 0;
         for (int i = 0; i < num_fault_dofs_; i++)
         {
            if (a_values_(i) < b_val) { vw_count++; }
         }

         real_t a_min = a_values_.Min();
         real_t a_max = a_values_.Max();

         os << "  Depth range: [" << z_max / 1000.0 << ", "
            << z_min / 1000.0 << "] km (local rank)\n";
         os << "  VW DOFs: " << vw_count << " (local rank)\n";
         os << "  VS DOFs: " << num_fault_dofs_ - vw_count << " (local rank)\n";
         os << "  a range: [" << a_min << ", " << a_max << "] (local rank)\n";
         os << "  eta: " << eta_values_(0) / 1e6 << " MPa·s/m\n";
      }
   }

private:
   BP2Params params_;
   BP5Params bp5_params_;
   MPIContext *mpi_ctx_ = nullptr;
   bool is_bp5_ = false;
   int num_fault_dofs_;
   int nbf_per_face_ = 1;         // basis functions per face
   int num_fault_faces_ = 0;      // number of fault faces
   int num_local_fault_dofs_ = 0;
   int num_global_fault_dofs_ = 0;
   Vector depths_;      // z-coordinates of fault DOFs
   Vector a_values_;    // a for each DOF (from depth in BP2, from (x2,x3) in BP5)
   Vector eta_values_;  // η for each DOF
   Vector dc_values_;   // Dc/L for each DOF (BP5: spatially varying)
   Vector tau_pre_;     // Pre-stress [2*N for BP5, N for BP2]
   Vector V_init_vec_;  // Initial velocity [2*N for BP5]
   Vector coords_x2_;   // Along-strike coordinate
   Vector coords_x3_;   // Depth coordinate

   // Phase 6.A: per-DOF 3-D coordinates and basis (BP5/3-D ctor only)
   Vector       dof_coords_3d_;   // [3 * num_fault_dofs_]
   DenseMatrix  dof_basis_;       // [9 x num_fault_dofs_], col i = [n_i; t1_i; t2_i]
   int          num_dof_basis_fallbacks_ = 0;     // sum of both kinds
   int          num_zero_normal_fallbacks_ = 0;   // basis slot zeroed
   int          num_t1_fallbacks_ = 0;            // up-vector fallback

   // Phase 6 §5: per-DOF normal stress (populated by ComputeParams)
   Vector sigma_n_per_dof_;       // [num_fault_dofs_]
   bool   params_computed_ = false;

   // Phase 5a of spatial_dynamic_rupture_plan.md (rev-3): per-DOF
   // reference IntegrationPoint cache + per-DOF bulk-element ownership
   // cache, BOTH populated only by the NEW BP5 ctor `(BP5Params&,
   // Vector& coords, DenseMatrix& basis, Array<int>& dof_to_elem, int,
   // MPIContext*, std::vector<IntegrationPoint>&)`.
   //
   // The LEGACY BP5 ctor `(DomainOperator&, BP5Params&, MPIContext*)`
   // does NOT touch these — they stay default-constructed (empty).
   // Downstream IP-aware readers (e.g. the IP-aware overload of
   // `dynamic/spatial_setup.hpp::InitializeFaultDOFs_Spatial`) MUST
   // gate on `fault_dof_ip_size() > 0` and call back to the centroid-
   // based path when zero, so the legacy BP5 / TPV byte-exact
   // contract is unaffected.
   std::vector<mfem::IntegrationPoint> fault_dof_ip_;
   Array<int>                          dof_to_elem_new_ctor_;

   // MPI gather info (parallel only)
   std::vector<int> recv_counts_;
   std::vector<int> recv_displs_;

   /// @brief Compute recv_counts and recv_displs for MPI_Gatherv.
   void ComputeGatherInfo()
   {
#ifdef SEAS_USE_MPI
      if (!mpi_ctx_) { return; }
      int size = mpi_ctx_->Size();
      recv_counts_.resize(size);
      recv_displs_.resize(size);

      MPI_Allgather(&num_local_fault_dofs_, 1, MPI_INT,
                     recv_counts_.data(), 1, MPI_INT,
                     mpi_ctx_->GetComm());

      recv_displs_[0] = 0;
      for (int i = 1; i < size; i++)
      {
         recv_displs_[i] = recv_displs_[i-1] + recv_counts_[i-1];
      }
#endif
   }

   /// @brief Deduplicate gathered data by depth.
   ///
   /// Sorts indices by depth, merges entries within tolerance (averaging
   /// their values). This handles DG partition-boundary duplicates.
   static void DeduplicateByDepth(const Vector &raw_depths,
                                   const Vector &raw_data,
                                   Vector &dedup_depths,
                                   Vector &dedup_data,
                                   real_t depth_tol = 1.0)
   {
      int n = raw_depths.Size();
      if (n == 0) { dedup_depths.SetSize(0); dedup_data.SetSize(0); return; }

      // Sort indices by depth (descending, surface first)
      std::vector<int> idx(n);
      for (int i = 0; i < n; i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](int a, int b) {
         return raw_depths(a) > raw_depths(b);
      });

      // Merge duplicates within tolerance
      std::vector<real_t> d_depths, d_data;
      d_depths.reserve(n);
      d_data.reserve(n);

      int i = 0;
      while (i < n)
      {
         real_t z_sum = raw_depths(idx[i]);
         real_t v_sum = raw_data(idx[i]);
         int count = 1;
         int j = i + 1;
         while (j < n && std::abs(raw_depths(idx[j]) - raw_depths(idx[i])) < depth_tol)
         {
            z_sum += raw_depths(idx[j]);
            v_sum += raw_data(idx[j]);
            count++;
            j++;
         }
         d_depths.push_back(z_sum / count);
         d_data.push_back(v_sum / count);
         i = j;
      }

      int m = static_cast<int>(d_depths.size());
      dedup_depths.SetSize(m);
      dedup_data.SetSize(m);
      for (int k = 0; k < m; k++)
      {
         dedup_depths(k) = d_depths[k];
         dedup_data(k) = d_data[k];
      }
   }

   /// @brief Deduplicate multiple fields at once using the same depth merging.
   static void DeduplicateMultipleByDepth(
      const Vector &raw_depths,
      const std::vector<Vector> &raw_fields,
      Vector &dedup_depths,
      std::vector<Vector> &dedup_fields,
      real_t depth_tol = 1.0)
   {
      int n = raw_depths.Size();
      int nf = static_cast<int>(raw_fields.size());
      if (n == 0)
      {
         dedup_depths.SetSize(0);
         dedup_fields.resize(nf);
         for (int f = 0; f < nf; f++) { dedup_fields[f].SetSize(0); }
         return;
      }

      // Sort indices by depth (descending)
      std::vector<int> idx(n);
      for (int i = 0; i < n; i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](int a, int b) {
         return raw_depths(a) > raw_depths(b);
      });

      // Merge duplicates
      std::vector<real_t> d_depths;
      std::vector<std::vector<real_t>> d_fields(nf);
      d_depths.reserve(n);
      for (int f = 0; f < nf; f++) { d_fields[f].reserve(n); }

      int i = 0;
      while (i < n)
      {
         real_t z_sum = raw_depths(idx[i]);
         std::vector<real_t> v_sums(nf, 0.0);
         for (int f = 0; f < nf; f++)
         {
            v_sums[f] = raw_fields[f](idx[i]);
         }
         int count = 1;
         int j = i + 1;
         while (j < n && std::abs(raw_depths(idx[j]) - raw_depths(idx[i])) < depth_tol)
         {
            z_sum += raw_depths(idx[j]);
            for (int f = 0; f < nf; f++)
            {
               v_sums[f] += raw_fields[f](idx[j]);
            }
            count++;
            j++;
         }
         d_depths.push_back(z_sum / count);
         for (int f = 0; f < nf; f++)
         {
            d_fields[f].push_back(v_sums[f] / count);
         }
         i = j;
      }

      int m = static_cast<int>(d_depths.size());
      dedup_depths.SetSize(m);
      for (int k = 0; k < m; k++) { dedup_depths(k) = d_depths[k]; }

      dedup_fields.resize(nf);
      for (int f = 0; f < nf; f++)
      {
         dedup_fields[f].SetSize(m);
         for (int k = 0; k < m; k++) { dedup_fields[f](k) = d_fields[f][k]; }
      }
   }

public:
   /// @brief Build face deduplication mapping from gathered 2D coordinates.
   ///
   /// In parallel DG, shared fault faces at partition boundaries produce
   /// duplicate face blocks in the gathered DOF array. This method identifies
   /// duplicate faces by comparing face centroids (average of DOF coordinates)
   /// within a tolerance, and returns the list of unique face start indices.
   ///
   /// Usage:
   ///   1. Call once with gathered x2, x3 coordinates to get the mapping
   ///   2. Use ApplyFaceDedupMap() to deduplicate any gathered field
   ///
   /// @param raw_x2 Gathered along-strike coordinates [M]
   /// @param raw_x3 Gathered depth coordinates [M]
   /// @param nbf_per_face Number of DOFs per face (must divide M evenly)
   /// @param[out] unique_face_indices Index of first DOF for each unique face
   /// @param tol Coordinate tolerance for centroid matching [m]
   static void BuildFaceDedupMap(
      const Vector &raw_x2, const Vector &raw_x3,
      int nbf_per_face,
      std::vector<int> &unique_face_indices,
      real_t tol = 1.0)
   {
      int M = raw_x2.Size();
      unique_face_indices.clear();
      if (M == 0 || nbf_per_face < 1) { return; }
      if (M % nbf_per_face != 0)
      {
         // Cannot form complete faces; return all DOFs as individual "faces"
         for (int i = 0; i < M; i++) { unique_face_indices.push_back(i); }
         return;
      }

      int num_faces = M / nbf_per_face;

      // Compute face centroids
      std::vector<real_t> cx2(num_faces), cx3(num_faces);
      for (int f = 0; f < num_faces; f++)
      {
         real_t s2 = 0.0, s3 = 0.0;
         for (int k = 0; k < nbf_per_face; k++)
         {
            s2 += raw_x2(f * nbf_per_face + k);
            s3 += raw_x3(f * nbf_per_face + k);
         }
         cx2[f] = s2 / nbf_per_face;
         cx3[f] = s3 / nbf_per_face;
      }

      // Mark unique faces: a face is duplicate if its centroid matches
      // a previously seen face within tolerance.
      std::vector<bool> is_dup(num_faces, false);
      for (int f = 0; f < num_faces; f++)
      {
         if (is_dup[f]) { continue; }
         unique_face_indices.push_back(f * nbf_per_face);
         // Mark later faces with matching centroid as duplicates
         for (int g = f + 1; g < num_faces; g++)
         {
            if (is_dup[g]) { continue; }
            if (std::abs(cx2[f] - cx2[g]) < tol &&
                std::abs(cx3[f] - cx3[g]) < tol)
            {
               is_dup[g] = true;
            }
         }
      }
   }

   /// @brief Apply a face dedup mapping to a gathered scalar field.
   ///
   /// Extracts the unique face DOF blocks from a raw gathered vector.
   ///
   /// @param raw Raw gathered field [M]
   /// @param unique_face_indices From BuildFaceDedupMap (first DOF of each unique face)
   /// @param nbf_per_face DOFs per face
   /// @param[out] dedup Deduplicated field [num_unique * nbf_per_face]
   static void ApplyFaceDedupMap(
      const Vector &raw,
      const std::vector<int> &unique_face_indices,
      int nbf_per_face,
      Vector &dedup)
   {
      int num_unique = static_cast<int>(unique_face_indices.size());
      dedup.SetSize(num_unique * nbf_per_face);
      for (int u = 0; u < num_unique; u++)
      {
         int start = unique_face_indices[u];
         for (int k = 0; k < nbf_per_face; k++)
         {
            dedup(u * nbf_per_face + k) = raw(start + k);
         }
      }
   }

private:
   /// @brief Phase 6.A: populate per-DOF global 3-D coords and (n, t1, t2)
   /// basis from the domain operator, restrict to owned DOFs, and
   /// re-orthonormalise via Gram-Schmidt.
   ///
   /// Gram-Schmidt re-orthonormalises ALL THREE input vectors so the
   /// resulting (n, t1, t2) matches the elasticity operator's per-face
   /// FaultBasis sign convention element-wise (R-001 contract: the
   /// SAFS pre-stress must live in the same frame as the elastic
   /// traction produced by `FaultBasis::ProjectTraction`):
   ///
   ///   n_i  ← n_i / |n_i|
   ///   t1_i ← (t1_i − (t1_i·n_i) n_i),  then  t1_i ← t1_i / |t1_i|
   ///   t2_i ← (t2_i − (t2_i·n_i) n_i − (t2_i·t1_i) t1_i),
   ///                                          then  t2_i ← t2_i / |t2_i|
   ///
   /// Note: the cross-product form `t2 = n × t1` *loses* the
   /// FaultBasis sign-flip when `sign_flipped == true` (see
   /// fault_basis.hpp:463-471); projecting the input t2 preserves it.
   ///
   /// Degenerate-case fallback (plan §1768-1773): if the projected t1
   /// drops below 1e-12 in magnitude (sub-vertical or sub-horizontal
   /// face where the up-vector lies in the face plane), derive a fresh
   /// t1 from the reference up-vector via `t1 = n × (up × n)`, then
   /// re-derive `t2 = n × t1_fallback` because the original input t2
   /// is unreliable when its companion t1 was degenerate (R-103). The
   /// number of fallbacks is counted and reported.
   void ComputePerDOFCoordsAndBasis_(DomainOperator<MeshType> &domain_op)
   {
      // R-205: reset fallback counters BEFORE the early-return branch
      // so the antiplane-fallback path (full_coords/full_basis empty)
      // is observably "no fallbacks" rather than carrying over whatever
      // a previous invocation left.  Kept here at function entry so any
      // future re-use of this function does not leak stale counts into
      // the SAFS guard at rate_state_fault.hpp:214.
      num_dof_basis_fallbacks_   = 0;
      num_zero_normal_fallbacks_ = 0;
      num_t1_fallbacks_          = 0;

      // Gather full local (interior + shared) coords + basis from the operator
      Vector full_coords;
      DenseMatrix full_basis;
      domain_op.GetFaultDOFCoords3D(full_coords);
      domain_op.GetFaultDOFBasis(full_basis);

      // If the operator did not populate either (e.g. antiplane fallback),
      // leave both members empty.
      if (full_coords.Size() == 0 || full_basis.Height() == 0)
      {
         dof_coords_3d_.SetSize(0);
         dof_basis_.SetSize(0, 0);
         return;
      }

      const int num_full = full_basis.Width();
      MFEM_VERIFY(full_coords.Size() == 3 * num_full,
                  "ComputePerDOFCoordsAndBasis_: coords size "
                  << full_coords.Size() << " inconsistent with basis cols "
                  << num_full);

      // Restrict 3-D coords to the owned fault DOF view (matches the size
      // of coords_x2_ / coords_x3_).
      domain_op.RestrictToOwnedFault(full_coords, dof_coords_3d_, 3);

      // Restrict the 9-row dense basis matrix to the owned view by
      // packing into a Vector with comps_per_dof = 9 and unpacking.
      Vector full_basis_flat(9 * num_full);
      for (int j = 0; j < num_full; j++)
      {
         for (int r = 0; r < 9; r++)
         {
            full_basis_flat(9 * j + r) = full_basis(r, j);
         }
      }
      Vector owned_basis_flat;
      domain_op.RestrictToOwnedFault(full_basis_flat, owned_basis_flat, 9);

      const int num_owned = num_fault_dofs_;
      MFEM_VERIFY(owned_basis_flat.Size() == 9 * num_owned,
                  "ComputePerDOFCoordsAndBasis_: owned basis flat size "
                  << owned_basis_flat.Size()
                  << " != 9 * num_owned " << 9 * num_owned);

      dof_basis_.SetSize(9, num_owned);

      // Reference up vector — matches the elasticity_operator setup.
      const real_t up_ref[3] = {0.0, 0.0, 1.0};

      // R-006 Gram-Schmidt re-orthonormalisation, mirroring Phase 2
      // `basis_to_node` (project_to_fault_stress.py).
      // R-205: counter reset moved to function preamble so the early-
      // return branch also observes a clean count.
      for (int i = 0; i < num_owned; i++)
      {
         real_t n[3] = { owned_basis_flat(9 * i + 0),
                         owned_basis_flat(9 * i + 1),
                         owned_basis_flat(9 * i + 2) };
         real_t t1[3] = { owned_basis_flat(9 * i + 3),
                          owned_basis_flat(9 * i + 4),
                          owned_basis_flat(9 * i + 5) };
         real_t t2[3] = { owned_basis_flat(9 * i + 6),
                          owned_basis_flat(9 * i + 7),
                          owned_basis_flat(9 * i + 8) };

         // n_i ← n_i / |n_i|
         const real_t n_len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
         if (n_len <= 1e-12)
         {
            // Zero-length normal: this DOF was not visited by any owned
            // fault face during basis population (e.g. boundary / unowned
            // DOFs in non-SAFS fixtures like elasticity_operator_tests).
            // BP5/TPV/antiplane paths do not read these per-DOF arrays.
            // Mark the slot with a sentinel zero basis and continue rather
            // than aborting; SAFS mode is required to refuse if any such
            // fallback occurred (see RateStateFaultOperator::SetSAFSMode).
            for (int d = 0; d < 9; d++) { dof_basis_(d, i) = 0.0; }
            num_zero_normal_fallbacks_++;
            num_dof_basis_fallbacks_++;
            continue;
         }
         const real_t inv_n = 1.0 / n_len;
         n[0] *= inv_n; n[1] *= inv_n; n[2] *= inv_n;

         // t1 ← t1 − (t1·n) n
         const real_t t1_dot_n = t1[0]*n[0] + t1[1]*n[1] + t1[2]*n[2];
         t1[0] -= t1_dot_n * n[0];
         t1[1] -= t1_dot_n * n[1];
         t1[2] -= t1_dot_n * n[2];
         real_t t1_len = std::sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);

         bool used_t1_fallback = false;
         if (t1_len < 1e-12)
         {
            // Degenerate fallback (plan §1768-1773): derive t1 from up.
            // t1 = n × (up × n) = up − (up·n) n  (since up is a unit vector
            // here, but to be safe normalise after).
            const real_t up_dot_n = up_ref[0]*n[0] + up_ref[1]*n[1]
                                    + up_ref[2]*n[2];
            t1[0] = up_ref[0] - up_dot_n * n[0];
            t1[1] = up_ref[1] - up_dot_n * n[1];
            t1[2] = up_ref[2] - up_dot_n * n[2];
            t1_len = std::sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
            MFEM_VERIFY(t1_len > 1e-12,
                        "ComputePerDOFCoordsAndBasis_: degenerate fallback at "
                        "DOF " << i << " — up vector nearly parallel to normal");
            num_t1_fallbacks_++;
            num_dof_basis_fallbacks_++;
            used_t1_fallback = true;
         }
         const real_t inv_t1 = 1.0 / t1_len;
         t1[0] *= inv_t1; t1[1] *= inv_t1; t1[2] *= inv_t1;

         // R-103: when the fallback re-derived t1 from `up`, the
         // original FaultBasis t2 is also unreliable (it was derived
         // from the same near-zero `up × n` that made t1 degenerate),
         // so the Gram-Schmidt projection below would collapse t2 to
         // near zero and trip the t2_len > 1e-12 abort.  In that case
         // re-derive t2 from the cross product `n × t1_fallback` —
         // the FaultBasis sign convention was undefined for this DOF
         // anyway, so sign loss is acceptable.
         if (used_t1_fallback)
         {
            t2[0] = n[1]*t1[2] - n[2]*t1[1];
            t2[1] = n[2]*t1[0] - n[0]*t1[2];
            t2[2] = n[0]*t1[1] - n[1]*t1[0];
         }

         // R-001 (post-fix): Gram-Schmidt re-orthonormalise the input
         // t2 against (n, t1) instead of redefining it as n × t1.
         //
         // The cross-product form `t2 = n × t1` *loses* the FaultBasis
         // sign-flip convention (Tandem: `t2_face = -strike_can` when
         // the CalcOrtho normal is anti-aligned with `ref_normal_` —
         // see fault_basis.hpp:463-471).  Dropping the sign-flip would
         // place the SAFS pre-stress in a different frame from the
         // elastic traction that `FaultBasis::ProjectTraction`
         // produces, leading to a mixed-frame `tau_vec = tau_pre +
         // traction` inside `RateStateFaultOperator` and reversed
         // slip direction on every sign-flipped DOF.
         //
         // Projecting the input `t2` preserves the sign convention
         // baked into FaultBasis output while still enforcing
         // orthonormality (defensive against FP drift from upstream).
         const real_t t2_dot_n  = t2[0]*n[0]  + t2[1]*n[1]  + t2[2]*n[2];
         const real_t t2_dot_t1 = t2[0]*t1[0] + t2[1]*t1[1] + t2[2]*t1[2];
         t2[0] -= t2_dot_n * n[0] + t2_dot_t1 * t1[0];
         t2[1] -= t2_dot_n * n[1] + t2_dot_t1 * t1[1];
         t2[2] -= t2_dot_n * n[2] + t2_dot_t1 * t1[2];
         const real_t t2_len = std::sqrt(t2[0]*t2[0] + t2[1]*t2[1]
                                         + t2[2]*t2[2]);
         MFEM_VERIFY(t2_len > 1e-12,
                     "ComputePerDOFCoordsAndBasis_: degenerate t2 at DOF "
                     << i);
         const real_t inv_t2 = 1.0 / t2_len;
         t2[0] *= inv_t2; t2[1] *= inv_t2; t2[2] *= inv_t2;

         for (int d = 0; d < 3; d++)
         {
            dof_basis_(d,     i) = n[d];
            dof_basis_(3 + d, i) = t1[d];
            dof_basis_(6 + d, i) = t2[d];
         }
      }

      if (num_dof_basis_fallbacks_ > 0)
      {
         mfem::err << "FaultGeometry WARNING: "
                   << num_zero_normal_fallbacks_ << " zero-normal + "
                   << num_t1_fallbacks_ << " t1-degenerate (of "
                   << num_owned << " owned fault DOFs) hit the basis "
                   "fallback path.  Zero-normal slots have a zeroed "
                   "basis and would silently project to zero pre-stress "
                   "and zero sigma_n in SAFS mode.  t1-degenerate slots "
                   "have a sign-flip-undefined basis (see lines 961-976) "
                   "that can place the projected sidecar stress "
                   "ANTIPARALLEL to the elastic traction at that DOF, "
                   "creating the positive-feedback failure documented in "
                   "CLAUDE.md (debug v8).  SetSAFSMode refuses to enable "
                   "in either case (see rate_state_fault.hpp:214).\n";
      }
   }

   /// @brief Compute 2D spatially varying parameters for BP5.
   ///
   /// Following Tandem's approach: ALL parameters are evaluated at each
   /// individual DOF's physical coordinates (per-DOF evaluation).
   ///
   /// Tandem reference: RateAndState.h:36-43 — set_params() iterates over
   /// all DOFs (numFaultFaces * nbf) and calls the Lua parameter function
   /// with each DOF's physical (x,y,z) coordinates.
   ///
   /// At p>=2 (multi-DOF), DOFs on the same face straddling the nucleation
   /// zone boundary will get different Dc, V_init, tau_pre values. This
   /// within-face discontinuity is handled correctly when combined with
   /// Tandem-style equilibrium initialization (--psi-init tandem), which
   /// absorbs all stress into psi so no overstress exists.
   ///
   /// Note: With SCEC initialization (delta_tau overstress), this per-DOF
   /// discontinuity combined with IP penalty can trigger instability at p>=2.
   /// The proper fix is to use --psi-init tandem, not to smooth the
   /// parameters (see bp5_debug_v46.md for analysis).
   void ComputeBP5Params()
   {
      a_values_.SetSize(num_fault_dofs_);
      eta_values_.SetSize(num_fault_dofs_);
      dc_values_.SetSize(num_fault_dofs_);
      tau_pre_.SetSize(2 * num_fault_dofs_);
      V_init_vec_.SetSize(2 * num_fault_dofs_);

      real_t eta = bp5_params_.eta();

      for (int i = 0; i < num_fault_dofs_; i++)
      {
         real_t x2 = coords_x2_(i);
         real_t x3 = coords_x3_(i);

         a_values_(i) = bp5_params_.a_of_x2_x3(x2, x3);
         eta_values_(i) = eta;
         dc_values_(i) = bp5_params_.Dc_of_x2_x3(x2, x3);

         real_t tau[2];
         bp5_params_.tau0_vec(x2, x3, tau);
         tau_pre_(2 * i)     = tau[0];
         tau_pre_(2 * i + 1) = tau[1];

         real_t Vi[2];
         bp5_params_.V_init_vec(x2, x3, Vi);
         V_init_vec_(2 * i)     = Vi[0];
         V_init_vec_(2 * i + 1) = Vi[1];

         // v58 diagnostic: dump full state for DOFs near fault tip
         // Disabled by default — enable via code flag if needed.
         if (false && std::abs(x2) > 45000.0 && x3 < 3000.0)
         {
            int rank = mpi_ctx_ ? mpi_ctx_->Rank() : 0;
            real_t psi_ss = bp5_params_.f0
                          + bp5_params_.b * std::log(bp5_params_.V0 / bp5_params_.Vp);
            mfem::out << std::scientific << std::setprecision(10)
                      << "[TIP-INIT] rank=" << rank
                      << " dof=" << i
                      << " x2=" << x2
                      << " x3=" << x3
                      << " a=" << a_values_(i)
                      << " b=" << bp5_params_.b
                      << " Dc=" << dc_values_(i)
                      << " f0=" << bp5_params_.f0
                      << " sigma_n=" << bp5_params_.sigma_n
                      << " eta=" << eta
                      << " V_init=(" << Vi[0] << "," << Vi[1] << ")"
                      << " tau_pre=(" << tau[0] << "," << tau[1] << ")"
                      << " psi_ss=" << psi_ss
                      << " psi_ss/a=" << psi_ss / a_values_(i)
                      << "\n";
         }
      }
   }

   /// @brief Compute depth-dependent parameters.
   void ComputeDepthDependentParams()
   {
      a_values_.SetSize(num_fault_dofs_);
      eta_values_.SetSize(num_fault_dofs_);

      // Radiation damping is constant for BP2
      real_t eta = params_.eta();

      for (int i = 0; i < num_fault_dofs_; i++)
      {
         // Compute a(z) from depth profile
         a_values_(i) = params_.a_of_z(depths_(i));

         // Constant radiation damping
         eta_values_(i) = eta;
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_GEOMETRY_HPP
