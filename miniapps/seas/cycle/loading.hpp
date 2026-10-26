// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// loading.hpp — 2D in-plane design rev 3, phase A (A4): the Dirichlet
// loading of the fully dynamic cycle, one function for both regimes.
//
// [cycle.loading] (cycle_config.hpp, LoadingConfig):
//   kind = "bp5" (default)  MakeBP5DirichletFunc(Vp), unchanged: the fault
//                           plane y = 0, slip along x, band 1 km
//   kind = "plate"          the plate motion of deck vectors, with
//                           d = (x - p) . m:
//       d >  b:  u = +(Vp t / 2) s
//       d < -b:  u = -(Vp t / 2) s
//       else:    u = sgn(-m . ref) Vp t s
// p plane_point, m plane_normal (unit), s slip_dir (unit, in the plane), b
// band, ref the [fault_geometry] reference normal.  Beyond the band the
// side m points to moves by +(Vp t / 2) s and the other side by the
// opposite.  Inside the band the value is the jump across the
// prescribed-slip (creep) faces: the operator multiplies it by
// ComputeSkeletonDirichletSign, +1 where the face's Elem1-outward normal
// points along the oriented normal n (n . ref > 0 on a planar fault), so
// the imposed jump is u(side opposite to n) - u(side n points to).  That
// jump equals the far-field difference u(side of m) - u(other side) =
// Vp t s when m points away from n, and its negative when m points along
// n: hence the factor sgn(-m . ref), which LoadingConfig::band_sign holds.
// The value is unchanged by (m, s) -> (-m, -s), the same plate motion.
//
// BP5 is p = 0, m = (0, 1, 0), s = (1, 0, 0), b = 1000 m with ref
// (0, -1, 0): band_sign = +1 and the "plate" values equal the BP5 function's
// bit for bit (tests/2d_inplane/test_plate_loading.cpp).  Components of s
// that are exactly 0 are written as +0.0, as the BP5 function writes them,
// so no -0.0 appears.  u(x, t) = t u(x, 1), as QDPhase::CheckLinearLoading_
// requires.  x and u have the mesh's dimension: on a 2-D mesh x = (x, y)
// and the deck's vectors have z = 0 (phase E: x(2) was read past the end).
//
// Scope (Chunhui Zhao, 2026-10-25).

#ifndef MFEM_SEAS_CYCLE_LOADING_HPP
#define MFEM_SEAS_CYCLE_LOADING_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../domain/boundary_config.hpp"

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The [cycle.loading] function with plate rate Vp (see the file comment).
inline DirichletFunc MakePlateLoading(const LoadingConfig &L, real_t Vp)
{
   if (!L.Plate()) { return MakeBP5DirichletFunc(Vp); }
   const std::array<real_t, 3> p = L.plane_point, m = L.plane_normal, s = L.slip_dir;
   const real_t b = L.band, band_sign = L.band_sign;
   return [Vp, p, m, s, b, band_sign](const Vector &x, real_t t, Vector &u)
   {
      const int dim = x.Size();
      u.SetSize(dim);
      u = 0.0;
      const real_t d = (dim == 3)
                       ? (x(0) - p[0]) * m[0] + (x(1) - p[1]) * m[1] + (x(2) - p[2]) * m[2]
                       : (x(0) - p[0]) * m[0] + (x(1) - p[1]) * m[1];
      real_t Vh = Vp * t;
      if (d > b)       { Vh *= 0.5; }
      else if (d < -b) { Vh *= -0.5; }
      else             { Vh *= band_sign; }
      for (int i = 0; i < dim; i++)
      {
         if (s[i] != 0.0) { u(i) = Vh * s[i]; }
      }
   };
}

/// The cycle's loading: [cycle.loading] at [cycle].plate_rate.
inline DirichletFunc MakeCycleLoading(const CycleConfig &cfg)
{
   return MakePlateLoading(cfg.loading, cfg.plate_rate);
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_LOADING_HPP
