// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_geometry_safs.inl — out-of-class definition of
// FaultGeometry<MeshType>::ComputeParams (Phase 6 §5 of
// PLAN_onfaultstress.md; renamed from `ComputeSAFSParams` to
// `ComputeParams` for general use across TPV102/104/205, TPV31, and
// SAFS).
//
// This file is split from `fault_geometry.hpp` so that the include of
// `io/field_coefficient.hpp` (which transitively pulls HDF5 / sidecar
// readers) only happens in the TUs that actually use the sidecar
// stress source.  Callers must include this header in addition to
// `fault_geometry.hpp` before invoking the non-templated overload.

#ifndef MFEM_SEAS_FAULT_GEOMETRY_SAFS_INL
#define MFEM_SEAS_FAULT_GEOMETRY_SAFS_INL

#include "fault_geometry.hpp"
#include "../io/field_coefficient.hpp"
#include "../io/stress_field_3d.hpp"

namespace mfem
{
namespace seas
{

template <typename MeshType>
void FaultGeometry<MeshType>::ComputeParams(
   const StressField3D& field,
   real_t P_p_pa,
   real_t P_p_grad_pa_per_m,
   real_t min_sigma_n_pa)
{
   MFEM_VERIFY(is_bp5_,
               "FaultGeometry::ComputeParams: only the 3-D / BP5 "
               "constructor populates per-DOF coords / basis; the BP2 "
               "ctor cannot be used here.");

   if (num_fault_dofs_ == 0)
   {
      sigma_n_per_dof_.SetSize(0);
      params_computed_ = true;
      return;
   }

   // Make sure the BP5 analytic per-DOF arrays (a, eta, Dc, V_init,
   // and the BP5 analytic tau_pre_) are populated.  ComputeParams
   // overwrites tau_pre_ but keeps the BP5 analytic forms for the
   // remaining state vectors (plan §1869).
   if (a_values_.Size() != num_fault_dofs_)
   {
      ComputeBP5Params();
   }

   // Rotate the sidecar stress tensor onto the per-DOF fault basis.
   // Sign convention: positive tau{1,2}_pre represents driving stress
   // in the +dip / +strike direction, matching the native
   // TPV102/104/205 drivers' DOFData::tau{1,2}_0 = +tau_ini convention.
   // The sign flip is applied INSIDE FieldProjector::ProjectFaultPreStress
   // — see io/field_coefficient.cpp.
   FieldProjector::ProjectFaultPreStress<MeshType>(
      field, *this,
      sigma_n_per_dof_,
      tau_pre_,
      P_p_pa,
      P_p_grad_pa_per_m,
      min_sigma_n_pa);

   params_computed_ = true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_GEOMETRY_SAFS_INL
