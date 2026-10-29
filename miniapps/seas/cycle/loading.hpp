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
// CheckCreepInBand (Chunhui Zhao, 2026-10-29): QDPhase refuses a mesh whose
// prescribed-slip faces (interior Dirichlet faces) leave the band, where the
// function gives a block value instead of the jump.
//
// Scope (Chunhui Zhao, 2026-10-25).

#ifndef MFEM_SEAS_CYCLE_LOADING_HPP
#define MFEM_SEAS_CYCLE_LOADING_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>

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

#ifdef MFEM_USE_MPI
/// The prescribed-slip faces against the band (Chunhui Zhao, 2026-10-29,
/// creep extensions in the UCERF 2-D model): the loading gives a face the
/// jump band_sign Vp t s only within |d| <= b of the plane; a prescribed-slip
/// face beyond the band would take a block value +-(Vp t / 2) s as its jump,
/// with no error.  The census counts, over all ranks, the interior boundary
/// elements with an attribute in `dirichlet` (the prescribed-slip faces: the
/// walls are exterior) with a vertex beyond the band, and the largest |d| of
/// those vertices.  The "bp5" kind is the plane y = 0 with b = 1 km, the
/// LoadingConfig defaults.  Collective.
struct CreepBandCensus
{
   long long faces = 0;    ///< prescribed-slip faces (a shared one may count twice)
   long long beyond = 0;   ///< of those, faces with a vertex at |d| > b
   real_t max_d = 0.0;     ///< the largest |d| of a vertex beyond the band [m]
};

inline CreepBandCensus CountCreepBeyondBand(ParMesh &pm, const std::set<int> &dirichlet,
                                            const LoadingConfig &L)
{
   const int dim = pm.Dimension();
   const std::array<real_t, 3> p = L.plane_point, m = L.plane_normal;
   CreepBandCensus c;
   Array<int> v;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (dirichlet.count(pm.GetBdrAttribute(be)) == 0) { continue; }
      if (pm.GetFaceInformation(pm.GetBdrElementFaceIndex(be)).IsBoundary()) { continue; }
      c.faces++;
      pm.GetBdrElementVertices(be, v);
      bool out = false;
      for (int k = 0; k < v.Size(); k++)
      {
         const real_t *x = pm.GetVertex(v[k]);
         real_t d = 0.0;
         for (int i = 0; i < dim; i++) { d += (x[i] - p[i]) * m[i]; }
         if (std::abs(d) > L.band)
         {
            out = true;
            c.max_d = std::max(c.max_d, std::abs(d));
         }
      }
      if (out) { c.beyond++; }
   }
   MPI_Comm comm = pm.GetComm();
   MPI_Allreduce(MPI_IN_PLACE, &c.faces, 1, MPI_LONG_LONG, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &c.beyond, 1, MPI_LONG_LONG, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &c.max_d, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   return c;
}

/// Refuses a mesh whose prescribed-slip faces leave the band (see above).
inline void CheckCreepInBand(ParMesh &pm, const std::set<int> &dirichlet, const LoadingConfig &L)
{
   const CreepBandCensus c = CountCreepBeyondBand(pm, dirichlet, L);
   MFEM_VERIFY(c.beyond == 0, "[cycle.loading]: " << c.beyond << " of " << c.faces
               << " prescribed-slip faces (interior Dirichlet faces) have a vertex beyond the "
               "band |d| <= " << L.band << " m of the loading plane (largest |d| " << c.max_d
               << " m); there the loading is a block value +-(Vp t / 2) s, not the jump "
               "Vp t s: put the prescribed-slip faces on the plane or widen band");
}
#endif

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_LOADING_HPP
