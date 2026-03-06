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

// Unit tests for BP5 inline mesh creation:
//   - Element count
//   - Boundary attributes
//   - Domain extent
//   - Fault interface detection

#include "mfem.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <set>
#include <limits>
#include <memory>

using namespace mfem;

// ============================================================================
// CreateBP5InlineMesh: inline 3D hex mesh for BP5 smoke tests
// ============================================================================

/// Create a 3D hex mesh for BP5 with proper boundary attributes.
///
/// Domain: [-Lx,Lx] x [-Ly,Ly] x [0,Lz]
/// Fault plane at x=0 (internal interface between adjacent elements).
///
/// Boundary attributes:
///   1 = x = -Lx (far-field, natural BC)
///   2 = x = +Lx (far-field, natural BC)
///   3 = y = +Ly (Dirichlet: plate loading +)
///   4 = y = -Ly (Dirichlet: plate loading -)
///   5 = z = 0   (free surface, natural BC)
///   6 = z = Lz  (deep boundary, natural BC)
///
/// @param nx Number of elements in x per half-domain (total 2*nx in x)
/// @param ny Number of elements in y per half-domain (total 2*ny in y)
/// @param nz Number of elements in z
/// @param Lx Half-domain size in x [m]
/// @param Ly Half-domain size in y [m]
/// @param Lz Domain depth [m]
std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz)
{
   // Create mesh: MakeCartesian3D makes [0, 2*Lx] x [0, 2*Ly] x [0, Lz]
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                            Element::HEXAHEDRON,
                            2.0 * Lx, 2.0 * Ly, Lz));

   // Shift to center: x -> x - Lx, y -> y - Ly (z stays [0, Lz])
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
   }

   // Assign boundary attributes based on face location
   const real_t tol = 1e-6 * std::max({Lx, Ly, Lz});

   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      // Get face center
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
      if (std::abs(cx - (-Lx)) < tol)
      {
         attr = 1;  // x = -Lx
      }
      else if (std::abs(cx - Lx) < tol)
      {
         attr = 2;  // x = +Lx
      }
      else if (std::abs(cy - Ly) < tol)
      {
         attr = 3;  // y = +Ly (Dirichlet +)
      }
      else if (std::abs(cy - (-Ly)) < tol)
      {
         attr = 4;  // y = -Ly (Dirichlet -)
      }
      else if (std::abs(cz) < tol)
      {
         attr = 5;  // z = 0 (free surface)
      }
      else if (std::abs(cz - Lz) < tol)
      {
         attr = 6;  // z = Lz (deep boundary)
      }
      else
      {
         attr = 1;  // fallback
      }

      mesh->SetBdrAttribute(i, attr);
   }

   // Ensure bdr_attributes array is up-to-date
   mesh->SetAttributes();

   return mesh;
}

// ============================================================================
// Test: Element count
// ============================================================================

void TestCreateBP5InlineMesh_ElementCount()
{
   std::cout << "\n=== Test: CreateBP5InlineMesh_ElementCount ===\n";

   int nx = 2, ny = 2, nz = 1;
   auto mesh = CreateBP5InlineMesh(nx, ny, nz, 200e3, 100e3, 100e3);

   int expected_ne = (2 * nx) * (2 * ny) * nz;  // 4*4*1 = 16
   TEST_ASSERT(mesh->GetNE() == expected_ne,
               "Element count = " + std::to_string(expected_ne));

   // Verify hex mesh
   TEST_ASSERT(mesh->Dimension() == 3, "Mesh dimension = 3");
   TEST_ASSERT(mesh->GetElement(0)->GetType() == Element::HEXAHEDRON,
               "Element type = HEXAHEDRON");

   // Also test a different configuration
   auto mesh2 = CreateBP5InlineMesh(4, 4, 2, 200e3, 100e3, 100e3);
   int expected2 = 8 * 8 * 2;  // 128
   TEST_ASSERT(mesh2->GetNE() == expected2,
               "Config B: element count = " + std::to_string(expected2));
}

// ============================================================================
// Test: Boundary attributes
// ============================================================================

void TestCreateBP5InlineMesh_BoundaryAttributes()
{
   std::cout << "\n=== Test: CreateBP5InlineMesh_BoundaryAttributes ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   // Collect all boundary attributes
   std::set<int> attrs;
   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      attrs.insert(mesh->GetBdrAttribute(i));
   }

   TEST_ASSERT(attrs.size() == 6, "6 distinct boundary attributes");
   TEST_ASSERT(attrs.count(1) == 1, "Attribute 1 present (x=-Lx)");
   TEST_ASSERT(attrs.count(2) == 1, "Attribute 2 present (x=+Lx)");
   TEST_ASSERT(attrs.count(3) == 1, "Attribute 3 present (y=+Ly)");
   TEST_ASSERT(attrs.count(4) == 1, "Attribute 4 present (y=-Ly)");
   TEST_ASSERT(attrs.count(5) == 1, "Attribute 5 present (z=0)");
   TEST_ASSERT(attrs.count(6) == 1, "Attribute 6 present (z=Lz)");

   // Count faces per attribute
   std::map<int, int> attr_count;
   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      attr_count[mesh->GetBdrAttribute(i)]++;
   }

   // For 4x4x1 mesh:
   // x faces: 4*1 = 4 per side → attr 1 and 2 each have 4
   // y faces: 4*1 = 4 per side → attr 3 and 4 each have 4
   // z faces: 4*4 = 16 per side → attr 5 and 6 each have 16
   TEST_ASSERT(attr_count[1] == 4, "Attr 1 (x=-Lx): 4 faces");
   TEST_ASSERT(attr_count[2] == 4, "Attr 2 (x=+Lx): 4 faces");
   TEST_ASSERT(attr_count[3] == 4, "Attr 3 (y=+Ly): 4 faces");
   TEST_ASSERT(attr_count[4] == 4, "Attr 4 (y=-Ly): 4 faces");
   TEST_ASSERT(attr_count[5] == 16, "Attr 5 (z=0): 16 faces");
   TEST_ASSERT(attr_count[6] == 16, "Attr 6 (z=Lz): 16 faces");
}

// ============================================================================
// Test: Domain extent
// ============================================================================

void TestCreateBP5InlineMesh_DomainExtent()
{
   std::cout << "\n=== Test: CreateBP5InlineMesh_DomainExtent ===\n";

   real_t Lx = 200e3, Ly = 100e3, Lz = 100e3;
   auto mesh = CreateBP5InlineMesh(2, 2, 1, Lx, Ly, Lz);

   real_t x_min = std::numeric_limits<real_t>::max();
   real_t x_max = -std::numeric_limits<real_t>::max();
   real_t y_min = x_min, y_max = x_max;
   real_t z_min = x_min, z_max = x_max;

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      const real_t *v = mesh->GetVertex(i);
      x_min = std::min(x_min, v[0]);
      x_max = std::max(x_max, v[0]);
      y_min = std::min(y_min, v[1]);
      y_max = std::max(y_max, v[1]);
      z_min = std::min(z_min, v[2]);
      z_max = std::max(z_max, v[2]);
   }

   TEST_NEAR(x_min, -Lx, 1.0, "x_min ≈ -200 km");
   TEST_NEAR(x_max, Lx, 1.0, "x_max ≈ +200 km");
   TEST_NEAR(y_min, -Ly, 1.0, "y_min ≈ -100 km");
   TEST_NEAR(y_max, Ly, 1.0, "y_max ≈ +100 km");
   TEST_NEAR(z_min, 0.0, 1.0, "z_min ≈ 0");
   TEST_NEAR(z_max, Lz, 1.0, "z_max ≈ 100 km");
}

// ============================================================================
// Test: Fault interface at x=0
// ============================================================================

void TestCreateBP5InlineMesh_FaultInterface()
{
   std::cout << "\n=== Test: CreateBP5InlineMesh_FaultInterface ===\n";

   // Use a mesh with even number of elements in x (so x=0 falls on element boundary)
   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   // Count internal faces near x=0
   int num_fault_faces = 0;
   int num_interior_faces = mesh->GetNumFaces() - mesh->GetNBE();
   const real_t tol = 1e-6 * 200e3;

   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      // Check if this is an interior face
      FaceElementTransformations *FTr = mesh->GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      // Get face center
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector center(3);
      FTr->Face->Transform(ip, center);

      if (std::abs(center(0)) < tol)
      {
         num_fault_faces++;
      }
   }

   // For a 4x4x1 mesh: x=0 splits at the middle, giving 4*1=4 fault faces
   TEST_ASSERT(num_fault_faces > 0,
               "Found interior faces at x≈0 (fault plane)");
   TEST_ASSERT(num_fault_faces == 4,
               "Expected 4 fault faces at x=0 for 4x4x1 mesh");
}

// ============================================================================
// Main
// ============================================================================

int main()
{
   std::cout << "BP5 Inline Mesh Unit Tests\n";
   std::cout << "=========================\n";

   TestCreateBP5InlineMesh_ElementCount();
   TestCreateBP5InlineMesh_BoundaryAttributes();
   TestCreateBP5InlineMesh_DomainExtent();
   TestCreateBP5InlineMesh_FaultInterface();

   TEST_PRINT_RESULTS();
   return num_failed;
}
