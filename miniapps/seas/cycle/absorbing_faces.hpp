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
//
// 2-D in-plane (2D in-plane design rev 3, section 6.4; Chunhui Zhao,
// 2026-10-27): on a triangle mesh in the (x, y) plane the faces are segments
// and n and up are taken in the plane, so `up` must have an in-plane part
// (BP3: up = [0, 1], the surface y = 0 faces up; the bottom y = -400 km has
// n . up = -1).  An up of +-e_z (the TPV205-2D frame) has none and is
// refused.
//
// Every outer face (faces = "all"; Chunhui Zhao, 2026-10-27: a box around
// the UCERF 2-D fault traces, loaded through its walls in QD, absorbs on
// every wall in the dynamic regime, with no free surface): with
// keep_free_surface = false no exterior `natural` face is kept, so `up` and
// `min_up` are not read (a map-view mesh's up, +-e_z, is not refused), and
// with walls on every exterior `dirichlet` face is a wall as above.

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
/// keep_free_surface = false: every exterior `natural` face absorbs ("all");
/// `up` and `min_up` are then not read.
inline OuterFaceTags RetagOuterFaces(ParMesh &m, const std::set<int> &natural,
                                     const std::set<int> &dirichlet, bool walls,
                                     const std::array<real_t, 3> &up, real_t min_up,
                                     bool keep_free_surface = true)
{
   const int dim = m.Dimension();
   MFEM_VERIFY((dim == 3 || dim == 2) && m.SpaceDimension() == dim,
               "RetagOuterFaces: 3-D meshes, or 2-D meshes in the (x, y) plane (dimension "
               << dim << ", space dimension " << m.SpaceDimension() << ")");
   real_t un = 0.0;   // |up| in the mesh's space
   if (keep_free_surface)
   {
      for (int c = 0; c < dim; c++) { un += up[c] * up[c]; }
      un = std::sqrt(un);
      MFEM_VERIFY(un > 0.0, "RetagOuterFaces: the up vector (" << up[0] << ", " << up[1] << ", "
                  << up[2] << ") has no component in the mesh's " << dim << "-D space"
                  << (dim == 2 ? " (on a 2-D mesh it must lie in the (x, y) plane)" : ""));
      MFEM_VERIFY(min_up > -1.0 && min_up <= 1.0,
                  "RetagOuterFaces: min_up = " << min_up << " must be in (-1, 1]");
   }

   OuterFaceTags r;
   const int amax = m.bdr_attributes.Size() ? m.bdr_attributes.Max() : 0;   // global
   r.natural_to = amax + 1;
   r.dirichlet_to = walls ? amax + 2 : 0;

   long long cnt[3] = {0, 0, 0};
   Vector n(dim), xf(dim), xe(dim);
   const Geometry::Type fgeom = (dim == 3) ? Geometry::TRIANGLE : Geometry::SEGMENT;
   const IntegrationPoint &cf = Geometries.GetCenter(fgeom);
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
      if (!keep_free_surface)
      {
         m.SetBdrAttribute(be, r.natural_to);
         cnt[0]++;
         continue;
      }
      MFEM_VERIFY(m.GetBdrElementGeometry(be) == fgeom, "RetagOuterFaces: boundary element "
                  << be << " is not a " << (dim == 3 ? "triangle" : "segment"));
      FaceElementTransformations *T = m.GetBdrFaceTransformations(be);
      MFEM_VERIFY(T, "RetagOuterFaces: no face transformation for boundary element " << be);
      T->SetAllIntPoints(&cf);
      CalcOrtho(T->Jacobian(), n);
      T->Transform(cf, xf);
      ElementTransformation &E = *T->Elem1;
      E.Transform(Geometries.GetCenter(E.GetGeometryType()), xe);
      real_t out = 0.0, nn = 0.0, nu = 0.0;
      for (int c = 0; c < dim; c++)
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
