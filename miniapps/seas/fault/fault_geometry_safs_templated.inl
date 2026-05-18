// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_geometry_safs_templated.inl — out-of-class definition of the
// TEMPLATED `FaultGeometry<MeshType>::ComputeSAFSParams<StressSource>`
// overload (Phase 3b of spatial_dynamic_rupture_plan.md rev-3).
//
// This file is split from `fault_geometry_safs.inl` so that:
//   1. The existing non-templated overload's body (defined in
//      `fault_geometry_safs.inl`) is preserved VERBATIM — no edits.
//   2. The templated body is a *separate* function with its own copy
//      of the per-DOF projection loop.  It does NOT call into the
//      existing `FieldProjector::ProjectFaultPreStress` helper (which
//      is hard-coded to `StressField3D` and stays unchanged for the
//      BP5 byte-exact contract).
//
// Callers that want the templated overload include both this file AND
// `fault_geometry.hpp`.  BP5 / TPV callers that only want the
// non-templated path include `fault_geometry.hpp` + `fault_geometry_safs.inl`
// (the existing pattern) and never instantiate this template.

#ifndef MFEM_SEAS_FAULT_GEOMETRY_SAFS_TEMPLATED_INL
#define MFEM_SEAS_FAULT_GEOMETRY_SAFS_TEMPLATED_INL

#include "fault_geometry.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{

template <typename MeshType>
template <typename StressSource>
void FaultGeometry<MeshType>::ComputeSAFSParams(
   const StressSource& source,
   real_t P_p_pa,
   real_t P_p_grad_pa_per_m,
   real_t min_sigma_n_pa)
{
   MFEM_VERIFY(is_bp5_,
               "FaultGeometry::ComputeSAFSParams<StressSource>: only the "
               "3-D / BP5 constructor populates per-DOF coords / basis; "
               "the BP2 ctor cannot be used in SAFS mode.");

   if (num_fault_dofs_ == 0)
   {
      // Match the non-templated overload's early-out exactly: do NOT
      // touch tau_pre_ here.  (See fault_geometry_safs.inl §37–42.)
      sigma_n_per_dof_.SetSize(0);
      safs_params_computed_ = true;
      return;
   }

   // BP5 analytic per-DOF arrays must be populated before SAFS-mode
   // overwrites tau_pre_ (matches the non-templated overload pattern;
   // §fault_geometry_safs.inl).
   if (a_values_.Size() != num_fault_dofs_)
   {
      ComputeBP5Params();
   }

   const Vector&      coords = dof_coords_3d_;
   const DenseMatrix& basis  = dof_basis_;

   MFEM_VERIFY(coords.Size() == 3 * num_fault_dofs_,
               "ComputeSAFSParams<StressSource>: dof_coords_3d_ size "
               << coords.Size() << " != 3 * num_fault_dofs " <<
               (3 * num_fault_dofs_));
   MFEM_VERIFY(basis.Height() == 9 && basis.Width() == num_fault_dofs_,
               "ComputeSAFSParams<StressSource>: dof_basis_ shape ("
               << basis.Height() << ", " << basis.Width()
               << ") != (9, " << num_fault_dofs_ << ")");

   sigma_n_per_dof_.SetSize(num_fault_dofs_);
   tau_pre_.SetSize(2 * num_fault_dofs_);
   sigma_n_per_dof_ = 0.0;
   tau_pre_ = 0.0;

   int clamp_count = 0;
   for (int i = 0; i < num_fault_dofs_; ++i)
   {
      const real_t x = coords(3 * i + 0);
      const real_t y = coords(3 * i + 1);
      const real_t z = coords(3 * i + 2);

      const real_t n[3]  = { basis(0, i), basis(1, i), basis(2, i) };
      const real_t t1[3] = { basis(3, i), basis(4, i), basis(5, i) };
      const real_t t2[3] = { basis(6, i), basis(7, i), basis(8, i) };

      // StressSource3D concept: Evaluate(x, y, z) returns a 3x3
      // symmetric Cauchy tensor, compression POSITIVE.
      const DenseMatrix S = source.Evaluate(x, y, z);

      real_t Sn[3];
      for (int r = 0; r < 3; ++r)
      {
         Sn[r] = S(r, 0) * n[0] + S(r, 1) * n[1] + S(r, 2) * n[2];
      }
      const real_t sigma_n_total = n[0]*Sn[0]  + n[1]*Sn[1]  + n[2]*Sn[2];
      const real_t tau1          = t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2];
      const real_t tau2          = t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2];

      const real_t depth_below = std::max(static_cast<real_t>(0.0), -z);
      const real_t P_p         = P_p_pa + P_p_grad_pa_per_m * depth_below;
      real_t sigma_n_eff       = sigma_n_total - P_p;

      if (min_sigma_n_pa > 0.0 && sigma_n_eff < min_sigma_n_pa)
      {
         sigma_n_eff = min_sigma_n_pa;
         clamp_count++;
      }

      sigma_n_per_dof_(i)        = sigma_n_eff;
      tau_pre_(2 * i)            = tau1;    // dip
      tau_pre_(2 * i + 1)        = tau2;    // strike
   }

   if (clamp_count > 0)
   {
      mfem::out << "FaultGeometry::ComputeSAFSParams<StressSource>: clamped "
                << clamp_count << " / " << num_fault_dofs_
                << " fault DOFs to min_sigma_n_pa = "
                << min_sigma_n_pa << " Pa\n";
   }

   safs_params_computed_ = true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_GEOMETRY_SAFS_TEMPLATED_INL
