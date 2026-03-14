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
      // Get fault DOF count and depths from domain operator
      num_fault_dofs_ = domain_op.GetNumFaultDOFs();

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

      // Get depths from domain operator
      domain_op.GetFaultDepths(depths_);

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
      num_fault_dofs_ = domain_op.GetNumFaultDOFs();
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

      // Get 2D fault coordinates
      domain_op.GetFaultCoords2D(coords_x2_, coords_x3_);

      // Also store depths for compatibility
      depths_.SetSize(num_fault_dofs_);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         depths_(i) = coords_x3_(i);
      }

      // Precompute per-DOF parameters using BP5 2D functions
      ComputeBP5Params();
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

   /// @brief Compute 2D spatially varying parameters for BP5.
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
