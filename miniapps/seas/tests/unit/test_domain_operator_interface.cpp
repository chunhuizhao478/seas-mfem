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
/// Tests Phase 1b additions to the DomainOperator base class:
/// 1. NumSlipComponents() returns 1 for antiplane
/// 2. GetFaultCoords2D() returns (0, depths) by default
/// 3. GetOffFaultDisplacement() returns empty by default
/// 4. Polymorphic dispatch through base pointer

#include "mfem.hpp"
#include "../../domain/domain_operator.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/antiplane_bdrload_operator.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test infrastructure
// =============================================================================

static int num_tests_passed = 0;
static int num_tests_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << std::endl; \
         std::cerr << "  at " << __FILE__ << ":" << __LINE__ << std::endl; \
         num_tests_failed++; \
         return false; \
      } \
   } while (0)

#define TEST_ASSERT_NEAR(a, b, tol, message) \
   do { \
      real_t diff = std::abs((a) - (b)); \
      if (diff > tol) { \
         std::cerr << "FAILED: " << message << std::endl; \
         std::cerr << "  Expected: " << (b) << ", Got: " << (a) \
                   << ", Diff: " << diff << ", Tol: " << tol << std::endl; \
         std::cerr << "  at " << __FILE__ << ":" << __LINE__ << std::endl; \
         num_tests_failed++; \
         return false; \
      } \
   } while (0)

#define RUN_TEST(test_func) \
   do { \
      std::cout << "Running " << #test_func << "... "; \
      std::cout.flush(); \
      if (test_func()) { \
         std::cout << "PASSED" << std::endl; \
         num_tests_passed++; \
      } else { \
         std::cout << "FAILED" << std::endl; \
      } \
   } while (0)

// =============================================================================
// Test: NumSlipComponents() returns 1 for antiplane
// =============================================================================

bool test_num_slip_components_antiplane()
{
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   TEST_ASSERT(domain.NumSlipComponents() == 1,
               "Antiplane should have 1 slip component, got "
               << domain.NumSlipComponents());

   return true;
}

// =============================================================================
// Test: GetFaultCoords2D() default returns (0, depths)
// =============================================================================

bool test_get_fault_coords_2d_default()
{
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
               "coords_x2 size should be " << n << ", got " << coords_x2.Size());
   TEST_ASSERT(coords_x3.Size() == n,
               "coords_x3 size should be " << n << ", got " << coords_x3.Size());

   // x2 (along-strike) should be all zeros for 2D antiplane
   for (int i = 0; i < n; i++)
   {
      TEST_ASSERT_NEAR(coords_x2(i), 0.0, 1e-14,
                        "coords_x2[" << i << "] should be 0");
   }

   // x3 (depth) should match GetFaultDepths
   Vector depths;
   domain.GetFaultDepths(depths);
   for (int i = 0; i < n; i++)
   {
      TEST_ASSERT_NEAR(coords_x3(i), depths(i), 1e-14,
                        "coords_x3[" << i << "] should match depth");
   }

   return true;
}

// =============================================================================
// Test: GetOffFaultDisplacement() default returns empty
// =============================================================================

bool test_off_fault_displacement_default()
{
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
               "Default GetOffFaultDisplacement should return empty, got size "
               << displacements.Size());

   return true;
}

// =============================================================================
// Test: Polymorphic dispatch through base pointer
// =============================================================================

bool test_polymorphic_dispatch()
{
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
               "GetFaultCoords2D through base: x2 size mismatch");

   // GetOffFaultDisplacement through base
   std::vector<Vector> pts;
   Vector disp;
   base->GetOffFaultDisplacement(pts, disp);
   TEST_ASSERT(disp.Size() == 0,
               "GetOffFaultDisplacement through base should return empty");

   return true;
}

// =============================================================================
// Test: BdrLoad operator also supports interface methods
// =============================================================================

bool test_bdrload_interface_methods()
{
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
               "BdrLoad GetFaultCoords2D x2 size mismatch");

   // All x2 should be zero
   for (int i = 0; i < x2.Size(); i++)
   {
      TEST_ASSERT_NEAR(x2(i), 0.0, 1e-14,
                        "BdrLoad coords_x2[" << i << "] should be 0");
   }

   return true;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "=== SEAS DomainOperator Interface Tests ===" << std::endl;
   std::cout << std::endl;

   RUN_TEST(test_num_slip_components_antiplane);
   RUN_TEST(test_get_fault_coords_2d_default);
   RUN_TEST(test_off_fault_displacement_default);
   RUN_TEST(test_polymorphic_dispatch);
   RUN_TEST(test_bdrload_interface_methods);

   std::cout << std::endl;
   std::cout << "=== Results: " << num_tests_passed << " passed, "
             << num_tests_failed << " failed ===" << std::endl;

   return num_tests_failed > 0 ? 1 : 0;
}
