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

/// @file test_domain_operator_interface.cpp
/// @brief Unit tests for DomainOperator interface methods
///
/// Tests Phase 2a additions to the DomainOperator base class:
/// 1. NumSlipComponents() returns 1 for antiplane
/// 2. GetFaultCoords2D() returns (0, depths) by default
/// 3. GetOffFaultDisplacement() returns empty by default
/// 4. Polymorphic dispatch through base pointer
/// 5. GetFaultBasis() returns nullptr by default

#include "test_macros.hpp"
#include "../../domain/domain_operator.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/antiplane_bdrload_operator.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../config/bp2_params.hpp"

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test: NumSlipComponents() returns 1 for antiplane
// =============================================================================

void test_num_slip_components_antiplane()
{
   std::cout << "\n=== NumSlipComponents Antiplane ===\n";

   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   TEST_ASSERT(domain.NumSlipComponents() == 1,
               "Antiplane should have 1 slip component");
}

// =============================================================================
// Test: GetFaultCoords2D() default returns (0, depths)
// =============================================================================

void test_get_fault_coords_2d_default()
{
   std::cout << "\n=== GetFaultCoords2D Default ===\n";

   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   Vector coords_x2, coords_x3;
   domain.GetFaultCoords2D(coords_x2, coords_x3);

   int n = domain.GetNumFaultDOFs();
   TEST_ASSERT(coords_x2.Size() == n,
               "coords_x2 size should match NumFaultDOFs");
   TEST_ASSERT(coords_x3.Size() == n,
               "coords_x3 size should match NumFaultDOFs");

   // x2 (along-strike) should be all zeros for 2D antiplane
   for (int i = 0; i < n; i++)
   {
      TEST_NEAR(coords_x2(i), 0.0, 1e-14,
                "coords_x2 should be 0 for antiplane");
   }

   // x3 (depth) should match GetFaultDepths
   Vector depths;
   domain.GetFaultDepths(depths);
   for (int i = 0; i < n; i++)
   {
      TEST_NEAR(coords_x3(i), depths(i), 1e-14,
                "coords_x3 should match depth");
   }
}

// =============================================================================
// Test: GetOffFaultDisplacement() default returns empty
// =============================================================================

void test_off_fault_displacement_default()
{
   std::cout << "\n=== GetOffFaultDisplacement Default ===\n";

   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   std::vector<Vector> points(3);
   for (auto &p : points)
   {
      p.SetSize(2);
      p = 0.0;
   }

   Vector displacements;
   domain.GetOffFaultDisplacement(points, displacements);

   // Default implementation returns empty
   TEST_ASSERT(displacements.Size() == 0,
               "Default GetOffFaultDisplacement should return empty");
}

// =============================================================================
// Test: Polymorphic dispatch through base pointer
// =============================================================================

void test_polymorphic_dispatch()
{
   std::cout << "\n=== Polymorphic Dispatch ===\n";

   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> concrete(*mesh, order, mu, Vp, Wf);

   // Access through base pointer
   DomainOperator<Mesh> *base = &concrete;

   // All interface methods should work through base pointer
   TEST_ASSERT(base->NumComponents() == 1,
               "NumComponents through base should be 1");
   TEST_ASSERT(base->Dimension() == 2,
               "Dimension through base should be 2");
   TEST_ASSERT(base->NumSlipComponents() == 1,
               "NumSlipComponents through base should be 1");
   TEST_ASSERT(base->GetShearModulus() == mu,
               "GetShearModulus through base should match");
   TEST_ASSERT(base->GetNumFaultDOFs() > 0,
               "GetNumFaultDOFs through base should be > 0");

   // GetFaultCoords2D through base
   Vector x2, x3;
   base->GetFaultCoords2D(x2, x3);
   TEST_ASSERT(x2.Size() == base->GetNumFaultDOFs(),
               "GetFaultCoords2D through base: x2 size matches");

   // GetFaultBasis through base
   TEST_ASSERT(base->GetFaultBasis() == nullptr,
               "GetFaultBasis default should return nullptr");

   // GetOffFaultDisplacement through base
   std::vector<Vector> pts;
   Vector disp;
   base->GetOffFaultDisplacement(pts, disp);
   TEST_ASSERT(disp.Size() == 0,
               "GetOffFaultDisplacement through base should return empty");
}

// =============================================================================
// Test: BdrLoad operator also supports interface methods
// =============================================================================

void test_bdrload_interface_methods()
{
   std::cout << "\n=== BdrLoad Interface Methods ===\n";

   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneBdrLoadOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   // These inherited defaults should work on BdrLoad operator too
   TEST_ASSERT(domain.NumSlipComponents() == 1,
               "BdrLoad NumSlipComponents should be 1");

   Vector x2, x3;
   domain.GetFaultCoords2D(x2, x3);
   TEST_ASSERT(x2.Size() == domain.GetNumFaultDOFs(),
               "BdrLoad GetFaultCoords2D x2 size matches");

   // All x2 should be zero
   for (int i = 0; i < x2.Size(); i++)
   {
      TEST_NEAR(x2(i), 0.0, 1e-14,
                "BdrLoad coords_x2 should be 0");
   }

   // GetFaultBasis default
   TEST_ASSERT(domain.GetFaultBasis() == nullptr,
               "BdrLoad GetFaultBasis default should return nullptr");
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "========================================\n";
   std::cout << "  DomainOperator Interface Tests\n";
   std::cout << "========================================\n";

   test_num_slip_components_antiplane();
   test_get_fault_coords_2d_default();
   test_off_fault_displacement_default();
   test_polymorphic_dispatch();
   test_bdrload_interface_methods();

   TEST_PRINT_RESULTS();
   return (num_failed > 0) ? 1 : 0;
}
