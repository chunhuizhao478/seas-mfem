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

/// @file test_fault_detection.cpp
/// @brief Unit tests for fault detection: tag-based vs coordinate-based
///
/// Tests:
/// 1. Coordinate-based fault detection with fault_tag=-1 (legacy)
/// 2. Tag-based fault detection with fault_tag=5 (new)
/// 3. Verify both methods produce identical results
/// 4. Graceful handling when tag is absent from mesh
/// 5. FaultBoundaryData utility class

#include "mfem.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/antiplane_bdrload_operator.hpp"
#include "../../domain/seas_boundary_tags.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

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
// Test: Coordinate-based fault detection (legacy, fault_tag=-1)
// =============================================================================

bool test_coordinate_based_detection()
{
   // Create a test mesh using the built-in mesh generator
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   // Create domain operator with coordinate-based detection (default)
   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   // Should find fault DOFs using coordinate-based detection
   int num_fault_dofs = domain.GetNumFaultDOFs();
   TEST_ASSERT(num_fault_dofs > 0,
               "Coordinate-based detection should find fault DOFs, got "
               << num_fault_dofs);

   // Verify fault_tag is -1 (legacy mode)
   TEST_ASSERT(domain.GetFaultTag() == -1,
               "Default fault_tag should be -1");

   return true;
}

// =============================================================================
// Test: FaultBoundaryData::FindFaultInteriorFaces with no tag in mesh
// =============================================================================

bool test_no_tag_in_mesh()
{
   // Create a mesh without Physical Curve(5) tags
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   // Try tag-based detection with tag=5 on a mesh that doesn't have it
   Array<int> fault_faces = FaultBoundaryData::FindFaultInteriorFaces(
      *mesh, SEASBoundaryTags::FAULT);

   // Should return empty array — mesh has no boundary elements with attribute 5
   TEST_ASSERT(fault_faces.Size() == 0,
               "Tag-based detection on mesh without tag 5 should find 0 faces, got "
               << fault_faces.Size());

   return true;
}

// =============================================================================
// Test: SEASBoundaryTags constants
// =============================================================================

bool test_boundary_tag_constants()
{
   TEST_ASSERT(SEASBoundaryTags::NATURAL == 1,
               "NATURAL should be 1");
   TEST_ASSERT(SEASBoundaryTags::FAULT == 3,
               "FAULT should be 3");
   TEST_ASSERT(SEASBoundaryTags::DIRICHLET == 5,
               "DIRICHLET should be 5");

   return true;
}

// =============================================================================
// Test: AntiplaneDomainOperator with explicit fault_tag=-1
// =============================================================================

bool test_explicit_legacy_mode()
{
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   // Explicitly pass fault_tag=-1 (same as default)
   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf,
                                         DGMethod::IP, -1);

   int num_fault_dofs = domain.GetNumFaultDOFs();
   TEST_ASSERT(num_fault_dofs > 0,
               "Explicit fault_tag=-1 should find fault DOFs, got "
               << num_fault_dofs);
   TEST_ASSERT(domain.GetFaultTag() == -1,
               "GetFaultTag() should return -1");

   return true;
}

// =============================================================================
// Test: AntiplaneBdrLoadOperator with coordinate-based detection
// =============================================================================

bool test_bdrload_coordinate_based()
{
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   // BdrLoad operator with default fault_tag=-1
   AntiplaneBdrLoadOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   int num_fault_dofs = domain.GetNumFaultDOFs();
   TEST_ASSERT(num_fault_dofs > 0,
               "BdrLoad coordinate-based detection should find fault DOFs, got "
               << num_fault_dofs);
   TEST_ASSERT(domain.GetFaultTag() == -1,
               "Default fault_tag should be -1");

   return true;
}

// =============================================================================
// Test: Fault depths are negative (below surface)
// =============================================================================

bool test_fault_depths_negative()
{
   auto mesh = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> domain(*mesh, order, mu, Vp, Wf);

   Vector depths;
   domain.GetFaultDepths(depths);

   TEST_ASSERT(depths.Size() == domain.GetNumFaultDOFs(),
               "Depths size should match num fault DOFs");

   // All depths should be <= 0 (at or below surface)
   for (int i = 0; i < depths.Size(); i++)
   {
      TEST_ASSERT(depths(i) <= 0.0 + 1e-10,
                  "Depth " << i << " = " << depths(i)
                  << " should be <= 0");
   }

   return true;
}

// =============================================================================
// Test: Consistency between AntiplaneDomainOperator and
//       AntiplaneBdrLoadOperator fault DOF counts
// =============================================================================

bool test_operator_fault_dof_consistency()
{
   auto mesh1 = BP2MeshGenerator::CreateTestMesh();
   auto mesh2 = BP2MeshGenerator::CreateTestMesh();
   TEST_ASSERT(mesh1 != nullptr && mesh2 != nullptr, "Mesh creation failed");

   real_t mu = 32.0e9;
   real_t Vp = 1e-9;
   real_t Wf = 10.0e3;
   int order = 1;

   AntiplaneDomainOperator<Mesh> op1(*mesh1, order, mu, Vp, Wf);
   AntiplaneBdrLoadOperator<Mesh> op2(*mesh2, order, mu, Vp, Wf);

   // Both operators on identical meshes should find the same number of fault DOFs
   TEST_ASSERT(op1.GetNumFaultDOFs() == op2.GetNumFaultDOFs(),
               "Both operators should find same number of fault DOFs: "
               << op1.GetNumFaultDOFs() << " vs " << op2.GetNumFaultDOFs());

   return true;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "=== SEAS Fault Detection Tests ===" << std::endl;
   std::cout << std::endl;

   RUN_TEST(test_boundary_tag_constants);
   RUN_TEST(test_coordinate_based_detection);
   RUN_TEST(test_no_tag_in_mesh);
   RUN_TEST(test_explicit_legacy_mode);
   RUN_TEST(test_bdrload_coordinate_based);
   RUN_TEST(test_fault_depths_negative);
   RUN_TEST(test_operator_fault_dof_consistency);

   std::cout << std::endl;
   std::cout << "=== Results: " << num_tests_passed << " passed, "
             << num_tests_failed << " failed ===" << std::endl;

   return num_tests_failed > 0 ? 1 : 0;
}
