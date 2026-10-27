// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_cfl.hpp — Nitsche dynamic rupture, plan item N-02.
//
// (1) ComputeMaxDt: explicit time-step bound for the central-difference
//     (Newmark beta = 0, gamma = 1/2) elastodynamic step,
//
//        dt = cfl * min_e  h_e / (p^2 c_p,e),   c_p,e = sqrt((lam_e+2 mu_e)/rho_e),
//
//     with (lam_e, mu_e, rho_e) from the per-element NitscheMaterialSample and
//     h_e the per-element length of BimaterialWaveOperator
//     (dynamic/bimaterial_wave_operator.inl:250-276: inscribed diameter
//     6 V / A_total for tets).  Same structure as
//     BimaterialWaveOperator::ComputeMaxDt with p^2 in place of its factor;
//     global MIN over ranks.  The measured stable constant is cfl = C_CFL.
//
// (2) RunCentralDifference: n steps of MFEM's own CentralDifferenceSolver
//     (linalg/ode.hpp: NewmarkSolver(0.0, 0.5)).  No stepping arithmetic of
//     our own.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_NITSCHE_CFL_HPP
#define MFEM_SEAS_NITSCHE_CFL_HPP

#include "mfem.hpp"
#include "../domain/nitsche_wave_operator.hpp"

#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// Per-element CFL length: inscribed diameter 6 V / A_total for a tet
/// (vol^{1/3} if the face area vanishes), as BimaterialWaveOperator.
inline real_t NitscheElementLength(Mesh &mesh, int e)
{
   MFEM_VERIFY(mesh.GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
               "NitscheElementLength: element " << e << " is a "
               << Geometry::Name[mesh.GetElementBaseGeometry(e)]
               << "; the Nitsche path supports tetrahedra only.");
   const real_t vol = mesh.GetElementVolume(e);
   Array<int> vert;
   mesh.GetElementVertices(e, vert);
   Vector v0(mesh.GetVertex(vert[0]), 3);
   Vector v1(mesh.GetVertex(vert[1]), 3);
   Vector v2(mesh.GetVertex(vert[2]), 3);
   Vector v3(mesh.GetVertex(vert[3]), 3);
   auto tri_area = [](const Vector &a, const Vector &b,
                      const Vector &c) -> real_t {
      real_t e1[3] = {b(0)-a(0), b(1)-a(1), b(2)-a(2)};
      real_t e2[3] = {c(0)-a(0), c(1)-a(1), c(2)-a(2)};
      real_t cx = e1[1]*e2[2] - e1[2]*e2[1];
      real_t cy = e1[2]*e2[0] - e1[0]*e2[2];
      real_t cz = e1[0]*e2[1] - e1[1]*e2[0];
      return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
   };
   const real_t A_total = tri_area(v0, v1, v2) + tri_area(v0, v2, v3)
                        + tri_area(v0, v3, v1) + tri_area(v1, v2, v3);
   return (A_total > 0) ? 6.0 * vol / A_total : std::pow(vol, 1.0 / 3.0);
}

/// dt = cfl * min_e h_e / (p^2 c_p,e), global MIN over ranks.
inline real_t ComputeMaxDt(const NitscheMaterialSample &sample, ParMesh &mesh,
                           int order, real_t cfl)
{
   MFEM_VERIFY(order >= 1, "ComputeMaxDt: order must be >= 1");
   const ParFiniteElementSpace *fes = sample.RhoGF().ParFESpace();
   const real_t p2 = real_t(order) * real_t(order);
   real_t local_dt_min = std::numeric_limits<real_t>::infinity();
   Array<int> dofs;
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      fes->GetElementDofs(e, dofs);
      const real_t lam = sample.LambdaGF()(dofs[0]);
      const real_t mu  = sample.MuGF()(dofs[0]);
      const real_t rho = sample.RhoGF()(dofs[0]);
      const real_t cp_e = std::sqrt((lam + 2.0 * mu) / rho);
      const real_t h_e  = NitscheElementLength(mesh, e);
      const real_t dt_e = cfl * h_e / (p2 * cp_e);
      local_dt_min = std::min(local_dt_min, dt_e);
   }
   real_t dt_global = local_dt_min;
   MPI_Allreduce(&local_dt_min, &dt_global, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, mesh.GetComm());
   return dt_global;
}

/// n steps of MFEM's CentralDifferenceSolver from (U, Udot) at time t.
/// A fresh solver per call: its first Step evaluates A^0 = op.Mult once.
inline void RunCentralDifference(SecondOrderTimeDependentOperator &op,
                                 Vector &U, Vector &Udot, real_t &t,
                                 real_t dt, int nsteps)
{
   CentralDifferenceSolver ode;
   ode.Init(op);
   for (int n = 0; n < nsteps; ++n)
   {
      real_t dt_n = dt;
      ode.Step(U, Udot, t, dt_n);
   }
}

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_CFL_HPP
