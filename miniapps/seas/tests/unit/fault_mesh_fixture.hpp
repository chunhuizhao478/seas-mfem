// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_mesh_fixture.hpp — shared small-mesh fixture for the spatial-seas
// (quasi-dynamic) unit tests.  Builds a structured hex mesh with a y=0 fault
// (Tandem convention: attr 3 = fault, 1 = Natural top/bottom, 5 = Dirichlet
// far-field).  Extracted from the copies in test_elasticity_heterogeneous_*.cpp
// (R-503).  Header-inline; each test executable links its own copy.

#ifndef MFEM_SEAS_TEST_FAULT_MESH_FIXTURE_HPP
#define MFEM_SEAS_TEST_FAULT_MESH_FIXTURE_HPP

#include "mfem.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{
namespace test
{

// Add attr-3 internal boundary elements at the y=0 interior faces (fault).
inline void AddFaultBoundaryElements(Mesh &mesh, real_t tol = 1e-6)
{
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *FTr = mesh.GetInteriorFaceTransformations(f);
      if (!FTr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector center(3);
      FTr->Face->Transform(ip, center);
      if (std::abs(center(1)) > tol) { continue; }

      Array<int> verts;
      mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 4)
      { mesh.AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3); }
      else if (verts.Size() == 3)
      { mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
}

// Domain [-Lx,Lx] x [-Ly,Ly] x [-Lz,0]; attr 1=Natural (z top/bottom),
// 5=Dirichlet (far-field walls), 3=fault (y=0).
inline Mesh CreateTestMesh3D(int nx, int ny, int nz,
                             real_t Lx, real_t Ly, real_t Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz, Element::HEXAHEDRON,
                                     2.0 * Lx, 2.0 * Ly, Lz);
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] -= Lx; v[1] -= Ly; v[2] -= Lz;
   }
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(3);
      T->Transform(ip, center);
      const real_t tol = 1e-6;
      if (std::abs(center(2)) < tol || std::abs(center(2) + Lz) < tol)
      { mesh.SetBdrAttribute(be, 1); }   // Natural
      else
      { mesh.SetBdrAttribute(be, 5); }   // Dirichlet
   }
   AddFaultBoundaryElements(mesh);
   return mesh;
}

} // namespace test
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TEST_FAULT_MESH_FIXTURE_HPP
