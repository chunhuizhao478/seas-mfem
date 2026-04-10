// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_BP5_MESH_UTILS_HPP
#define MFEM_SEAS_BP5_MESH_UTILS_HPP

#include "mfem.hpp"
#include <memory>
#include <cmath>
#include <algorithm>

namespace mfem
{
namespace seas
{

/// Create a 3D hex mesh for BP5 with Tandem boundary attributes.
/// Domain: [-Lx,Lx] x [-Ly,Ly] x [-Lz,0]
/// Boundary attributes (Tandem tags):
///   1 = Natural (z=0 top, z=-Lz bottom)
///   5 = Dirichlet (x=±Lx, y=±Ly far-field)
///
/// @param nx  Elements per half-domain in x
/// @param ny  Elements per half-domain in y
/// @param nz  Elements in z (full domain)
/// @param Lx  Half-domain size in x [m]
/// @param Ly  Half-domain size in y [m]
/// @param Lz  Domain depth [m] (positive value; mesh spans [-Lz, 0])
inline std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz)
{
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                            Element::HEXAHEDRON,
                            2.0 * Lx, 2.0 * Ly, Lz));

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
      v[2] -= Lz;  // Z ranges [-Lz, 0]
   }

   const real_t tol = 1e-6 * std::max({Lx, Ly, Lz});

   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      Array<int> vertices;
      mesh->GetBdrElementVertices(i, vertices);

      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int j = 0; j < vertices.Size(); j++)
      {
         const real_t *v = mesh->GetVertex(vertices[j]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= vertices.Size();
      cy /= vertices.Size();
      cz /= vertices.Size();

      int attr;
      if (std::abs(cz) < tol || std::abs(cz + Lz) < tol)
      {
         attr = 1;  // Natural (top z=0, bottom z=-Lz)
      }
      else
      {
         attr = 5;  // Dirichlet (far-field x=±Lx, y=±Ly)
      }

      mesh->SetBdrAttribute(i, attr);
   }

   // Mark interior faces at y=0 as fault boundary elements (attr 3).
   // Without this, the mesh has no fault DOFs and the simulation produces
   // trivially zero output.
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      auto *FTr = mesh->GetInteriorFaceTransformations(f);
      if (!FTr) { continue; }
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector center(3);
      FTr->Face->Transform(ip, center);
      if (std::abs(center(1)) > tol) { continue; }

      Array<int> verts;
      mesh->GetFaceVertices(f, verts);
      if (verts.Size() == 4)
      {
         mesh->AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3);
      }
      else if (verts.Size() == 3)
      {
         mesh->AddBdrTriangle(verts[0], verts[1], verts[2], 3);
      }
   }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP5_MESH_UTILS_HPP
