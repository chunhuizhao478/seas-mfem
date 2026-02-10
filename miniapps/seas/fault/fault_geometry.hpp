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
#include "../domain/domain_operator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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
   FaultGeometry(DomainOperator<MeshType> &domain_op, const BP2Params &params)
      : params_(params)
   {
      // Get fault DOF count and depths from domain operator
      num_fault_dofs_ = domain_op.GetNumFaultDOFs();

      if (num_fault_dofs_ == 0)
      {
         MFEM_WARNING("FaultGeometry: No fault DOFs found in domain operator");
         return;
      }

      // Get depths from domain operator
      domain_op.GetFaultDepths(depths_);

      // Compute depth-dependent parameters
      ComputeDepthDependentParams();
   }

   /// @brief Number of fault DOFs (for state vector sizing).
   int NumFaultDOFs() const { return num_fault_dofs_; }

   /// @brief Get depths at each fault DOF (z coordinate, negative below surface).
   const Vector &GetDepths() const { return depths_; }

   /// @brief Get rate-state parameter a at each fault DOF (precomputed from depth).
   const Vector &GetAValues() const { return a_values_; }

   /// @brief Get radiation damping η at each fault DOF.
   const Vector &GetEtaValues() const { return eta_values_; }

   /// @brief Get the BP2 parameters.
   const BP2Params &GetParams() const { return params_; }

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

   /// @brief Check if a depth is in the velocity-weakening zone.
   ///
   /// @param z Depth coordinate (negative below surface)
   /// @return True if a(z) < b (velocity-weakening)
   bool IsVelocityWeakening(real_t z) const
   {
      return params_.a_of_z(z) < params_.b;
   }

   /// @brief Get the VW/VS transition depth (top of transition zone).
   ///
   /// Returns the depth H where the transition from VW to VS begins.
   real_t GetVWDepth() const { return -params_.H; }

   /// @brief Get the full VS depth (bottom of transition zone).
   ///
   /// Returns the depth H+h where fully VS behavior begins.
   real_t GetVSDepth() const { return -(params_.H + params_.h); }

   /// @brief Get indices of DOFs in the velocity-weakening zone.
   void GetVWDOFs(Array<int> &vw_dofs) const
   {
      vw_dofs.SetSize(0);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         if (IsVelocityWeakening(depths_(i)))
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
         if (!IsVelocityWeakening(depths_(i)))
         {
            vs_dofs.Append(i);
         }
      }
   }

   /// @brief Print fault geometry information.
   void Print(std::ostream &os = mfem::out) const
   {
      os << "Fault Geometry:\n";
      os << "  Number of DOFs: " << num_fault_dofs_ << "\n";

      if (num_fault_dofs_ > 0)
      {
         // Find depth range
         real_t z_min = depths_.Min();
         real_t z_max = depths_.Max();
         os << "  Depth range: [" << z_max / 1000.0 << ", "
            << z_min / 1000.0 << "] km\n";

         // Count VW and VS DOFs
         int vw_count = 0;
         for (int i = 0; i < num_fault_dofs_; i++)
         {
            if (a_values_(i) < params_.b) { vw_count++; }
         }
         os << "  VW DOFs: " << vw_count << "\n";
         os << "  VS DOFs: " << num_fault_dofs_ - vw_count << "\n";

         // Show a(z) range
         real_t a_min = a_values_.Min();
         real_t a_max = a_values_.Max();
         os << "  a range: [" << a_min << ", " << a_max << "]\n";

         // Show eta (should be constant for BP2)
         os << "  eta: " << eta_values_(0) / 1e6 << " MPa·s/m\n";
      }
   }

private:
   BP2Params params_;
   int num_fault_dofs_;
   Vector depths_;      // z-coordinates of fault DOFs
   Vector a_values_;    // a(z) for each DOF
   Vector eta_values_;  // η for each DOF

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
