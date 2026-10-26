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
//     6 V / A_total for tets; 4 A / P for triangles, 2D in-plane design
//     rev 3, derivation 10.3: twice the inradius in both dimensions, so the
//     deck's cfl means the same in a 2D and a 3D run).  Same structure as
//     BimaterialWaveOperator::ComputeMaxDt with p^2 in place of its factor;
//     global MIN over ranks.  The measured stable constant is cfl = C_CFL.
//
// (2) RunCentralDifference: n steps of MFEM's own CentralDifferenceSolver
//     (linalg/ode.hpp: NewmarkSolver(0.0, 0.5)).  No stepping arithmetic of
//     our own.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only; triangles in 2D
// (2D in-plane design rev 3, phase C).

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
/// (vol^{1/3} if the face area vanishes), as BimaterialWaveOperator; for a
/// triangle 4 A / P (twice the inradius; sqrt(A) if the perimeter vanishes).
inline real_t NitscheElementLength(Mesh &mesh, int e)
{
   const Geometry::Type g = mesh.GetElementBaseGeometry(e);
   MFEM_VERIFY(g == Geometry::TETRAHEDRON || g == Geometry::TRIANGLE,
               "NitscheElementLength: element " << e << " is a "
               << Geometry::Name[g] << "; the Nitsche path supports "
               "tetrahedra (3D) and triangles (2D) only.");
   const real_t vol = mesh.GetElementVolume(e);
   Array<int> vert;
   mesh.GetElementVertices(e, vert);
   if (g == Geometry::TRIANGLE)
   {
      const int sd = mesh.SpaceDimension();
      auto dist = [&](int i, int j) {
         const real_t *xi = mesh.GetVertex(vert[i]), *xj = mesh.GetVertex(vert[j]);
         real_t d2 = 0.0;
         for (int c = 0; c < sd; c++) { d2 += (xi[c] - xj[c]) * (xi[c] - xj[c]); }
         return std::sqrt(d2);
      };
      const real_t P = dist(0, 1) + dist(1, 2) + dist(2, 0);
      return (P > 0) ? 4.0 * vol / P : std::sqrt(vol);
   }
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

/// The stable constant C(p) of dt = cfl C(p) h_min / (p^2 c_p), p = 1..4:
/// 2 / (omega_max dt_unit) on a structured simplex mesh, measured by power
/// iteration and by bisection on blow-up and rounded down.  3D: P3 on 24
/// tets per hex (the drivers' kCeig).  2D: 4 triangles per square
/// (tests/2d_inplane/test_2d_nitsche.cpp N2; 2D in-plane design rev 3
/// section 5.5: the 3D column does not carry over, triangles allow
/// 1.6 to 1.9 times more with h = 4 A / P).
inline real_t NitscheStableConstant(int dim, int p)
{
   static const real_t c3[5] = {0.0, 0.113, 0.212, 0.282, 0.332};
   static const real_t c2[5] = {0.0, 0.215, 0.365, 0.461, 0.526};
   MFEM_VERIFY((dim == 2 || dim == 3) && p >= 1 && p <= 4,
               "NitscheStableConstant: dim " << dim << ", p " << p
               << "; measured for dim 2 or 3 and p = 1..4 only");
   return (dim == 3) ? c3[p] : c2[p];
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
