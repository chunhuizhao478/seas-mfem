// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_cfl.hpp — the explicit time-step bound of the poroelastic option
// (poroelastic plan rev 3, step P1.2; poro-derivation-rev1.2 §9.4, eq. (9.6)):
//
//    dt = cfl * min_e h_e / (p^2 c_fast,e),
//
// ComputeMaxDt (solver/nitsche_cfl.hpp) with the fast P speed of eq. (4.6) in
// place of sqrt((lambda + 2 mu) / rho): the same element length h_e (6 V / A,
// NitscheElementLength), the same expression and a global MIN over ranks.
// c_fast,e comes from the sample's constants of element e (ComputeBiotSpeeds).
// The fitted gamma_w > 1/2 lowers the undamped Newmark bound only where the
// drag is stiff and the w modes are overdamped (derivation §9.4); test U2-18
// measures the stable constant.

#ifndef MFEM_SEAS_PORO_CFL_HPP
#define MFEM_SEAS_PORO_CFL_HPP

#include "mfem.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "poro_dispersion.hpp"
#include "poro_material.hpp"

#include <algorithm>
#include <limits>

namespace mfem
{
namespace seas
{
namespace poro
{

#ifdef MFEM_USE_MPI

/// dt = cfl * min_e h_e / (p^2 c_fast,e), global MIN over ranks.
inline real_t PoroMaxDt(const PoroMaterialSample &sample, ParMesh &mesh,
                        int order, real_t cfl)
{
   MFEM_VERIFY(order >= 1, "PoroMaxDt: order must be >= 1");
   const real_t p2 = real_t(order) * real_t(order);
   real_t local_dt_min = std::numeric_limits<real_t>::infinity();
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      const real_t c_e = ComputeBiotSpeeds(sample.Constants(e)).c_fast;
      const real_t h_e = NitscheElementLength(mesh, e);
      const real_t dt_e = cfl * h_e / (p2 * c_e);
      local_dt_min = std::min(local_dt_min, dt_e);
   }
   real_t dt_global = local_dt_min;
   MPI_Allreduce(&local_dt_min, &dt_global, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, mesh.GetComm());
   return dt_global;
}

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_CFL_HPP
