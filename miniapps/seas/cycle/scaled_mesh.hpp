// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// scaled_mesh.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.1):
// LoadScaledMesh and the mesh report the P0 driver stub prints.
//
// LoadScaledMesh reads the mesh and scales it exactly as seas_bp5_full does
// (tests/verification/bp5_verification_full.cpp:1441-1451): the vertices are
// multiplied by the scale before the ParMesh is built, so the QD and the
// dynamic operators see the same coordinates.  Mesh-only code: no solver
// header (rule M4).

#ifndef MFEM_SEAS_CYCLE_SCALED_MESH_HPP
#define MFEM_SEAS_CYCLE_SCALED_MESH_HPP

#include "mfem.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// Reads a 3-D tetrahedral mesh and multiplies every vertex by `scale`
/// (1000 turns the km mesh into m).  Refuses other element types and meshes
/// with a nodal grid function, whose nodes the vertex scaling would miss.
inline std::unique_ptr<Mesh> LoadScaledMesh(const std::string &path, real_t scale)
{
   MFEM_VERIFY(scale > 0.0, "LoadScaledMesh: scale = " << scale << " must be > 0");
   auto mesh = std::make_unique<Mesh>(path.c_str(), 1, 1);
   MFEM_VERIFY(mesh->Dimension() == 3, "LoadScaledMesh: " << path << " is not 3-D");
   for (int e = 0; e < mesh->GetNE(); e++)
   {
      MFEM_VERIFY(mesh->GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                  "LoadScaledMesh: tetrahedral meshes only; element " << e << " is "
                  << Geometry::Name[mesh->GetElementBaseGeometry(e)]);
   }
   MFEM_VERIFY(mesh->GetNodes() == nullptr,
               "LoadScaledMesh: " << path << " has a nodal grid function; "
               "scaling its vertices would leave the nodes unscaled");
   if (scale != 1.0)
   {
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         real_t *v = mesh->GetVertex(i);
         v[0] *= scale;
         v[1] *= scale;
         v[2] *= scale;
      }
      mesh->SetAttributes();   // as seas_bp5_full
   }
   return mesh;
}

/// Boundary faces of one attribute: how many, how many are exterior (one
/// adjacent element) or interior (two), and their bounding box.
struct BdrAttrCount
{
   int attr = 0;
   long long faces = 0;
   long long exterior = 0;
   long long interior = 0;
   long long triangles = 0;
   real_t lo[3] = {std::numeric_limits<real_t>::max(),
                   std::numeric_limits<real_t>::max(),
                   std::numeric_limits<real_t>::max()};
   real_t hi[3] = {std::numeric_limits<real_t>::lowest(),
                   std::numeric_limits<real_t>::lowest(),
                   std::numeric_limits<real_t>::lowest()};
};

/// One BdrAttrCount per boundary attribute of the serial mesh, by attribute.
inline std::map<int, BdrAttrCount> CountBoundary(const Mesh &mesh)
{
   std::map<int, BdrAttrCount> out;
   Array<int> verts;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      const int a = mesh.GetBdrAttribute(be);
      BdrAttrCount &c = out[a];
      c.attr = a;
      c.faces++;
      if (mesh.FaceIsInterior(mesh.GetBdrElementFaceIndex(be))) { c.interior++; }
      else { c.exterior++; }
      if (mesh.GetBdrElementGeometry(be) == Geometry::TRIANGLE) { c.triangles++; }
      mesh.GetBdrElementVertices(be, verts);
      for (int k = 0; k < verts.Size(); k++)
      {
         const real_t *x = mesh.GetVertex(verts[k]);
         for (int d = 0; d < 3; d++)
         {
            c.lo[d] = std::min(c.lo[d], x[d]);
            c.hi[d] = std::max(c.hi[d], x[d]);
         }
      }
   }
   return out;
}

/// Fault DOFs of the IP layout on one triangular face: the (p+1)(p+2)/2 nodal
/// points of the order-p triangle (domain/elasticity_operator_setup.inl:706-745).
inline int FaultDOFsPerFace(int order) { return (order + 1) * (order + 2) / 2; }

/// Fault faces held by this rank after the partition (each boundary element
/// goes to exactly one rank, mesh/pmesh.cpp BuildLocalBoundary) and how many
/// of them lie on a face shared with another rank.
struct LocalFaultFaces
{
   long long faces = 0;
   long long shared = 0;
};

inline LocalFaultFaces CountLocalFaultFaces(const ParMesh &pmesh, int fault_attr)
{
   std::set<int> shared_faces;
   for (int i = 0; i < pmesh.GetNSharedFaces(); i++)
   {
      shared_faces.insert(pmesh.GetSharedFace(i));
   }
   LocalFaultFaces c;
   for (int be = 0; be < pmesh.GetNBE(); be++)
   {
      if (pmesh.GetBdrAttribute(be) != fault_attr) { continue; }
      c.faces++;
      if (shared_faces.count(pmesh.GetBdrElementFaceIndex(be))) { c.shared++; }
   }
   return c;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_SCALED_MESH_HPP
