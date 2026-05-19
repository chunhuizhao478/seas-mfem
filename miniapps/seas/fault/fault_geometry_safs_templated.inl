// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_geometry_safs_templated.inl — out-of-class definition of the
// TEMPLATED `FaultGeometry<MeshType>::ComputeParams<StressSource>`
// overload (Phase 3b of spatial_dynamic_rupture_plan.md rev-3;
// renamed from `ComputeSAFSParams` to `ComputeParams` for general use
// across TPV102/104/205, TPV31, and SAFS).
//
// This file is split from `fault_geometry_safs.inl` so that:
//   1. The non-templated overload's body lives next to it in
//      `fault_geometry_safs.inl`.
//   2. The templated body is a *separate* function with its own copy
//      of the per-DOF projection loop.  It does NOT call into the
//      `FieldProjector::ProjectFaultPreStress` helper (which is
//      hard-coded to `StressField3D`).  The sign convention applied
//      to `tau1` / `tau2` MUST match the non-templated overload —
//      see below.
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
void FaultGeometry<MeshType>::ComputeParams(
   const StressSource& source,
   real_t P_p_pa,
   real_t P_p_grad_pa_per_m,
   real_t min_sigma_n_pa)
{
   MFEM_VERIFY(is_bp5_,
               "FaultGeometry::ComputeParams<StressSource>: only the "
               "3-D / BP5 constructor populates per-DOF coords / basis; "
               "the BP2 ctor cannot be used here.");

   if (num_fault_dofs_ == 0)
   {
      // Match the non-templated overload's early-out exactly: do NOT
      // touch tau_pre_ here.  (See fault_geometry_safs.inl §37–42.)
      sigma_n_per_dof_.SetSize(0);
      params_computed_ = true;
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
               "ComputeParams<StressSource>: dof_coords_3d_ size "
               << coords.Size() << " != 3 * num_fault_dofs " <<
               (3 * num_fault_dofs_));
   MFEM_VERIFY(basis.Height() == 9 && basis.Width() == num_fault_dofs_,
               "ComputeParams<StressSource>: dof_basis_ shape ("
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
      // Sign convention (CLAUDE.md project-wide, set by the proven native
      // TPV102/104/205 drivers): positive `tau2_0` represents the
      // right-lateral driving stress in the +strike direction, e.g.,
      //   drivers/tpv102_driver.cpp + dynamic/tpv102_setup.hpp:
      //      d.tau2_0 = TPV102Params::tau_ini   ( = +75e6 for σ_xy = +75e6 )
      //   drivers/tpv205_driver.cpp:
      //      d.tau2_0 = ComputeTau2_0_TPV205(...)   ( returns +tau_back etc. )
      //
      // The raw Cauchy projection T = σ·n with n = (0,-1,0) returns the
      // traction the +y side exerts on the −y side (which is in the
      // OPPOSITE sense to the driving-stress convention used by the
      // native drivers — Newton's 3rd-law mirror).  Negate `tau1` /
      // `tau2` so the projection result feeds `DOFData::tau{1,2}_0`
      // in the same sign convention as the native TPV setup.
      //
      // `sigma_n_total = n·S·n` is unchanged (it is sign-invariant
      // under n → −n).
      //
      // See debug_document/general_driver_debug_document/tpv102_tpv104_review.md R-001.
      const real_t tau1          = -(t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2]);
      const real_t tau2          = -(t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2]);

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
      mfem::out << "FaultGeometry::ComputeParams<StressSource>: clamped "
                << clamp_count << " / " << num_fault_dofs_
                << " fault DOFs to min_sigma_n_pa = "
                << min_sigma_n_pa << " Pa\n";
   }

   params_computed_ = true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_GEOMETRY_SAFS_TEMPLATED_INL
