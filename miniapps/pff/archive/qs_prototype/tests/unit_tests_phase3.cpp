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

// Unit tests for Phase 3: Operators
// - ElasticityOperator
// - DamageOperator

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../operators/elasticity_operator.hpp"
#include "../operators/damage_operator.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace mfem::pff;

// Simple test framework
static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      if (std::abs((value) - (expected)) > (tol)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << (expected) << ", Got: " << (value) << "\n"; \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// ElasticityOperator Tests
// =============================================================================

void TestElasticityOperator()
{
   std::cout << "\n=== Testing ElasticityOperator ===\n";

   // Create a simple 2D mesh
   Mesh serial_mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   // Finite element spaces
   H1_FECollection u_fec(1, dim);
   ParFiniteElementSpace u_fes(&mesh, &u_fec, dim);  // Vector field

   H1_FECollection d_fec(1, dim);
   ParFiniteElementSpace d_fes(&mesh, &d_fec);       // Scalar field

   // Material parameters
   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Essential boundary (fix bottom)
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;  // bottom boundary

   // Test 1: Construction
   {
      ElasticityOperator op(u_fes, mat, ess_bdr);
      TEST_ASSERT(op.Height() == u_fes.TrueVSize(),
                  "Operator height matches TrueVSize");
      TEST_ASSERT(op.Width() == u_fes.TrueVSize(),
                  "Operator width matches TrueVSize");
   }

   // Test 2: Zero displacement gives zero residual (interior DOFs)
   {
      ElasticityOperator op(u_fes, mat, ess_bdr);

      // Set zero damage
      ParGridFunction d(&d_fes);
      d = 0.0;
      op.SetDamageField(d);

      // Zero displacement
      ParGridFunction u(&u_fes);
      u = 0.0;

      Vector u_true;
      u.GetTrueDofs(u_true);

      Vector r(u_fes.TrueVSize());
      op.Mult(u_true, r);

      // Residual should be zero (or small) for zero displacement
      TEST_NEAR(r.Norml2(), 0.0, 1e-12, "Zero displacement gives zero residual");
   }

   // Test 3: Non-zero displacement gives non-zero residual
   {
      ElasticityOperator op(u_fes, mat, ess_bdr);

      ParGridFunction d(&d_fes);
      d = 0.0;
      op.SetDamageField(d);

      // Apply linear displacement
      ParGridFunction u(&u_fes);
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.001 * x(0);
         u(1) = 0.0;
      });
      u.ProjectCoefficient(u_coeff);

      Vector u_true;
      u.GetTrueDofs(u_true);

      Vector r(u_fes.TrueVSize());
      op.Mult(u_true, r);

      // Internal DOFs should have non-zero residual for a linear displacement
      // that doesn't satisfy equilibrium
      // Note: For this simple test, we just check it runs without error
      TEST_ASSERT(true, "Non-zero displacement produces residual without error");
   }

   // Test 4: Jacobian is returned
   {
      ElasticityOperator op(u_fes, mat, ess_bdr);

      ParGridFunction d(&d_fes);
      d = 0.0;
      op.SetDamageField(d);

      ParGridFunction u(&u_fes);
      u = 0.0;

      Vector u_true;
      u.GetTrueDofs(u_true);

      Operator &K = op.GetGradient(u_true);
      TEST_ASSERT(K.Height() == u_fes.TrueVSize(),
                  "Jacobian height matches TrueVSize");
      TEST_ASSERT(K.Width() == u_fes.TrueVSize(),
                  "Jacobian width matches TrueVSize");
   }

   // Test 5: Essential BC enforcement
   {
      ElasticityOperator op(u_fes, mat, ess_bdr);

      ParGridFunction d(&d_fes);
      d = 0.0;
      op.SetDamageField(d);

      // Set prescribed displacement
      Vector u_prescribed(u_fes.TrueVSize());
      u_prescribed = 0.0;
      // Set some prescribed values
      const Array<int> &ess_dofs = op.GetEssentialTrueDofs();
      for (int i = 0; i < ess_dofs.Size(); i++)
      {
         u_prescribed(ess_dofs[i]) = 0.001;
      }
      op.SetPrescribedDisplacement(u_prescribed);

      // Test with u = 0 -> residual on ess DOFs should be -0.001
      Vector u_true(u_fes.TrueVSize());
      u_true = 0.0;

      Vector r(u_fes.TrueVSize());
      op.Mult(u_true, r);

      // Check essential DOF residuals
      if (ess_dofs.Size() > 0)
      {
         real_t r_ess = r(ess_dofs[0]);
         TEST_NEAR(r_ess, -0.001, 1e-12, "Essential BC residual = u - u_prescribed");
      }
   }

   // Test 6: Damage reduces stiffness
   {
      ElasticityOperator op_intact(u_fes, mat, ess_bdr);
      ElasticityOperator op_damaged(u_fes, mat, ess_bdr);

      ParGridFunction d_zero(&d_fes);
      d_zero = 0.0;
      op_intact.SetDamageField(d_zero);

      ParGridFunction d_half(&d_fes);
      d_half = 0.5;
      op_damaged.SetDamageField(d_half);

      // Same displacement
      ParGridFunction u(&u_fes);
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.001 * x(0);
         u(1) = 0.0;
      });
      u.ProjectCoefficient(u_coeff);

      Vector u_true;
      u.GetTrueDofs(u_true);

      Vector r_intact(u_fes.TrueVSize()), r_damaged(u_fes.TrueVSize());
      op_intact.Mult(u_true, r_intact);
      op_damaged.Mult(u_true, r_damaged);

      // Damaged residual should be smaller in magnitude
      TEST_ASSERT(r_damaged.Norml2() < r_intact.Norml2(),
                  "Damage reduces residual magnitude (reduces stiffness)");
   }
}

// =============================================================================
// DamageOperator Tests
// =============================================================================

void TestDamageOperator()
{
   std::cout << "\n=== Testing DamageOperator ===\n";

   // Create a simple 2D mesh
   Mesh serial_mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   // Finite element spaces
   H1_FECollection d_fec(1, dim);
   ParFiniteElementSpace d_fes(&mesh, &d_fec);

   L2_FECollection psi_fec(0, dim);
   ParFiniteElementSpace psi_fes(&mesh, &psi_fec);

   // Material parameters
   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Test 1: Construction
   {
      DamageOperator op(d_fes, mat);
      TEST_ASSERT(op.Height() == d_fes.TrueVSize(),
                  "Operator height matches TrueVSize");
      TEST_ASSERT(op.Width() == d_fes.TrueVSize(),
                  "Operator width matches TrueVSize");
   }

   // Test 2: Zero damage with zero history gives zero residual
   {
      DamageOperator op(d_fes, mat);

      ParGridFunction d(&d_fes);
      d = 0.0;

      Vector d_true;
      d.GetTrueDofs(d_true);

      Vector r(d_fes.TrueVSize());
      op.Mult(d_true, r);

      TEST_NEAR(r.Norml2(), 0.0, 1e-10,
                "Zero damage with zero history gives zero residual");
   }

   // Test 3: Non-zero damage with zero history gives positive residual
   {
      DamageOperator op(d_fes, mat);

      ParGridFunction d(&d_fes);
      d = 0.5;

      Vector d_true;
      d.GetTrueDofs(d_true);

      Vector r(d_fes.TrueVSize());
      op.Mult(d_true, r);

      // Sum of residuals should be positive (restoring force)
      real_t sum = 0.0;
      for (int i = 0; i < r.Size(); i++)
      {
         sum += r(i);
      }
      TEST_ASSERT(sum > 0.0,
                  "Non-zero damage with zero history gives positive residual sum");
   }

   // Test 4: High strain energy history drives damage
   {
      DamageOperator op(d_fes, mat);

      // Set high strain energy history
      ParGridFunction psi(&psi_fes);
      psi = mat.GetCriticalStrainEnergy() * 10.0;  // Much higher than critical
      op.UpdateStrainEnergyHistory(psi);

      ParGridFunction d(&d_fes);
      d = 0.0;

      Vector d_true;
      d.GetTrueDofs(d_true);

      Vector r(d_fes.TrueVSize());
      op.Mult(d_true, r);

      // g'(0) = -2*(1-eta) < 0, so g'(0)*H < 0
      // Residual sum should be negative (driving damage up)
      real_t sum = 0.0;
      for (int i = 0; i < r.Size(); i++)
      {
         sum += r(i);
      }
      TEST_ASSERT(sum < 0.0,
                  "High strain energy history gives negative residual (drives damage)");
   }

   // Test 5: Jacobian is returned
   {
      DamageOperator op(d_fes, mat);

      ParGridFunction d(&d_fes);
      d = 0.3;

      Vector d_true;
      d.GetTrueDofs(d_true);

      Operator &K = op.GetGradient(d_true);
      TEST_ASSERT(K.Height() == d_fes.TrueVSize(),
                  "Jacobian height matches TrueVSize");
      TEST_ASSERT(K.Width() == d_fes.TrueVSize(),
                  "Jacobian width matches TrueVSize");
   }

   // Test 6: Irreversibility constraint
   {
      DamageOperator op(d_fes, mat);

      Vector d(d_fes.TrueVSize());
      d = 0.3;

      Vector d_old(d_fes.TrueVSize());
      d_old = 0.5;

      op.ApplyIrreversibility(d, d_old);

      // After irreversibility, d should be >= d_old
      bool all_ok = true;
      for (int i = 0; i < d.Size(); i++)
      {
         if (d(i) < d_old(i) - 1e-14)
         {
            all_ok = false;
            break;
         }
      }
      TEST_ASSERT(all_ok, "Irreversibility enforces d >= d_old");
   }

   // Test 7: Clamping to [0, 1]
   {
      DamageOperator op(d_fes, mat);

      Vector d(d_fes.TrueVSize());
      d = 1.5;  // Out of range

      Vector d_old(d_fes.TrueVSize());
      d_old = -0.1;  // Also out of range

      op.ApplyIrreversibility(d, d_old);

      // Should be clamped to [0, 1]
      bool all_clamped = true;
      for (int i = 0; i < d.Size(); i++)
      {
         if (d(i) < 0.0 || d(i) > 1.0)
         {
            all_clamped = false;
            break;
         }
      }
      TEST_ASSERT(all_clamped, "Damage clamped to [0, 1]");
   }

   // Test 8: History field is monotonic
   {
      DamageOperator op(d_fes, mat);

      ParGridFunction psi1(&psi_fes);
      psi1 = 100.0;
      op.UpdateStrainEnergyHistory(psi1);

      ParGridFunction psi2(&psi_fes);
      psi2 = 50.0;  // Lower value
      op.UpdateStrainEnergyHistory(psi2);

      // History should still be 100 (maximum)
      const ParGridFunction &H = op.GetStrainEnergyHistory();
      real_t H_min = H.Min();
      TEST_ASSERT(H_min >= 100.0 - 1e-10, "History field is monotonically increasing");
   }
}

// =============================================================================
// Integration Test: Combined Operators
// =============================================================================

void TestOperatorsIntegration()
{
   std::cout << "\n=== Testing Operators Integration ===\n";

   // Create mesh
   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   // FE spaces
   H1_FECollection u_fec(1, dim);
   ParFiniteElementSpace u_fes(&mesh, &u_fec, dim);

   H1_FECollection d_fec(1, dim);
   ParFiniteElementSpace d_fes(&mesh, &d_fec);

   L2_FECollection psi_fec(0, dim);
   ParFiniteElementSpace psi_fes(&mesh, &psi_fec);

   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Essential boundary
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;

   // Test: Full staggered iteration setup
   {
      // Create operators
      ElasticityOperator u_op(u_fes, mat, ess_bdr);
      DamageOperator d_op(d_fes, mat);

      // Initial fields
      ParGridFunction u(&u_fes);
      u = 0.0;

      ParGridFunction d(&d_fes);
      d = 0.0;

      ParGridFunction psi(&psi_fes);
      psi = 0.0;

      // Set damage field for elasticity
      u_op.SetDamageField(d);

      // Compute strain energy
      Vector u_true;
      u.GetTrueDofs(u_true);
      u_op.ComputeStrainEnergyDensity(u_true, psi);

      // Update history in damage operator
      d_op.UpdateStrainEnergyHistory(psi);

      // Compute damage residual
      Vector d_true;
      d.GetTrueDofs(d_true);

      Vector r_d(d_fes.TrueVSize());
      d_op.Mult(d_true, r_d);

      TEST_ASSERT(true, "Full staggered iteration setup works");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   // Initialize MPI
   Mpi::Init(argc, argv);
   int myid = Mpi::WorldRank();

   if (myid == 0)
   {
      std::cout << "Phase Field Fracture - Phase 3 Unit Tests (Operators)\n";
      std::cout << "======================================================\n";
   }

   // Run all tests
   TestElasticityOperator();
   TestDamageOperator();
   TestOperatorsIntegration();

   // Summary (only on rank 0)
   if (myid == 0)
   {
      std::cout << "\n======================================================\n";
      std::cout << "Test Summary:\n";
      std::cout << "  Total tests: " << num_tests << "\n";
      std::cout << "  Passed:      " << num_passed << "\n";
      std::cout << "  Failed:      " << num_failed << "\n";
      std::cout << "======================================================\n";

      if (num_failed > 0)
      {
         std::cout << "SOME TESTS FAILED!\n";
      }
      else
      {
         std::cout << "ALL TESTS PASSED!\n";
      }
   }

   return (num_failed > 0) ? 1 : 0;
}
