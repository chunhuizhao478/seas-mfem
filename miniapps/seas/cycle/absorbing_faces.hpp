// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// absorbing_faces.hpp — fully dynamic BP5 cycle, P5 (plan rev 3, Sections
// 6.4 and 11.2): the outer boundary faces that absorb in the dynamic regime
// ([cycle.absorbing]; Chunhui Zhao, 2026-10-26: dashpots on all outer walls,
// since the walls and the bottom send the event's waves back to the fault).
//
// The dynamic regime builds W on its own copy of the QD regime's ParMesh: the
// same nodes, elements and partition, so the same DOFs, with new boundary
// attributes on the outer faces.  The QD regime keeps its mesh and its
// boundaries, so nothing in QD changes.  RetagOuterFaces, on that copy:
//
//   for every boundary element on an exterior face (one element: the fault
//   and the creep plane are interior faces and are never touched):
//     attribute in `natural`, outward unit normal n, n . up < min_up
//         -> natural_to: absorbing, no load (the bottom; the free surface,
//            n . up >= min_up, keeps its attribute and stays traction-free)
//     attribute in `dirichlet`, walls on
//         -> dirichlet_to: absorbing, with the frozen wall load that
//            DynPhase::Enter forms (hook H3, plan 11.2)
//   SetAttributes (collective); the counts are global
//
// The outward normal is the face's Jacobian normal, flipped if it points into
// its element (checked against the element centroid).  No solver header
// (rule M4).

#ifndef MFEM_SEAS_CYCLE_ABSORBING_FACES_HPP
#define MFEM_SEAS_CYCLE_ABSORBING_FACES_HPP

#include "mfem.hpp"

#include <array>
#include <cmath>
#include <set>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The new attributes and the global face counts of RetagOuterFaces.
struct OuterFaceTags
{
   int natural_to = 0;          ///< attribute of the absorbing traction-free faces
   int dirichlet_to = 0;        ///< attribute of the absorbing walls (0: walls off)
   long long natural = 0;       ///< exterior faces moved to natural_to
   long long dirichlet = 0;     ///< exterior faces moved to dirichlet_to
   long long free_surface = 0;  ///< exterior `natural` faces kept (n . up >= min_up)
};

#ifdef MFEM_USE_MPI
/// Collective.  `m` is the dynamic regime's copy of the mesh; see the header.
inline OuterFaceTags RetagOuterFaces(ParMesh &m, const std::set<int> &natural,
                                     const std::set<int> &dirichlet, bool walls,
                                     const std::array<real_t, 3> &up, real_t min_up)
{
   const int dim = m.Dimension();
   MFEM_VERIFY(dim == 3 && m.SpaceDimension() == 3,
               "RetagOuterFaces: 3D meshes only (dimension " << dim << ")");
   const real_t un = std::sqrt(up[0] * up[0] + up[1] * up[1] + up[2] * up[2]);
   MFEM_VERIFY(un > 0.0, "RetagOuterFaces: the up vector is zero");
   MFEM_VERIFY(min_up > -1.0 && min_up <= 1.0,
               "RetagOuterFaces: min_up = " << min_up << " must be in (-1, 1]");

   OuterFaceTags r;
   const int amax = m.bdr_attributes.Size() ? m.bdr_attributes.Max() : 0;   // global
   r.natural_to = amax + 1;
   r.dirichlet_to = walls ? amax + 2 : 0;

   long long cnt[3] = {0, 0, 0};
   Vector n(3), xf(3), xe(3);
   const IntegrationPoint &cf = Geometries.GetCenter(Geometry::TRIANGLE);
   for (int be = 0; be < m.GetNBE(); be++)
   {
      const int a = m.GetBdrAttribute(be);
      const bool nat = natural.count(a) > 0;
      const bool dir = walls && dirichlet.count(a) > 0;
      if (!nat && !dir) { continue; }
      if (!m.GetFaceInformation(m.GetBdrElementFaceIndex(be)).IsBoundary()) { continue; }
      if (dir)
      {
         m.SetBdrAttribute(be, r.dirichlet_to);
         cnt[1]++;
         continue;
      }
      MFEM_VERIFY(m.GetBdrElementGeometry(be) == Geometry::TRIANGLE,
                  "RetagOuterFaces: boundary element " << be << " is not a triangle");
      FaceElementTransformations *T = m.GetBdrFaceTransformations(be);
      MFEM_VERIFY(T, "RetagOuterFaces: no face transformation for boundary element " << be);
      T->SetAllIntPoints(&cf);
      CalcOrtho(T->Jacobian(), n);
      T->Transform(cf, xf);
      ElementTransformation &E = *T->Elem1;
      E.Transform(Geometries.GetCenter(E.GetGeometryType()), xe);
      real_t out = 0.0, nn = 0.0, nu = 0.0;
      for (int c = 0; c < 3; c++)
      {
         out += (xf(c) - xe(c)) * n(c);
         nn += n(c) * n(c);
         nu += n(c) * up[c];
      }
      nu /= std::sqrt(nn) * un;
      if (out < 0.0) { nu = -nu; }   // the Jacobian normal pointed inward
      if (nu >= min_up)
      {
         cnt[2]++;
         continue;
      }
      m.SetBdrAttribute(be, r.natural_to);
      cnt[0]++;
   }
   MPI_Allreduce(MPI_IN_PLACE, cnt, 3, MPI_LONG_LONG, MPI_SUM, m.GetComm());
   r.natural = cnt[0];
   r.dirichlet = cnt[1];
   r.free_surface = cnt[2];
   m.SetAttributes();
   return r;
}
#endif // MFEM_USE_MPI

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_ABSORBING_FACES_HPP
