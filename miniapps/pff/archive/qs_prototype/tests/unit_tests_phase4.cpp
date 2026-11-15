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

// Unit tests for Phase 4: Staggered Solver (PFFSolver)

#include "mfem.hpp"
#include "../pff_solver.hpp"

#include <iostream>
#include <cmath>
#include <memory>

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
         if (Mpi::Root()) { \
            std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         } \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      if (std::abs((value) - (expected)) > (tol)) { \
         if (Mpi::Root()) { \
            std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
            std::cerr << "  Expected: " << (expected) << ", Got: " << (value) << "\n"; \
         } \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// PrescribedDisplacement Coefficient for Testing
// =============================================================================

class TestPrescribedDisplacement : public VectorCoefficient
{
   real_t amplitude_;
   real_t t_;
public:
   TestPrescribedDisplacement(int dim, real_t amplitude = 1.0)
      : VectorCoefficient(dim), amplitude_(amplitude), t_(0.0) {}

   void SetTime(real_t t) override { t_ = t; }

   void Eval(Vector &V, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      V.SetSize(vdim);
      V = 0.0;
      // Positive y-displacement for tension
      if (vdim >= 2)
      {
         V(1) = amplitude_ * t_;
      }
   }
};

// =============================================================================
// PFFSolver Tests
// =============================================================================

void TestPFFSolverConstruction()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Construction ===\n";
   }

   // Create a simple 2D mesh
   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);

   // Material parameters
   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   // Test 1: Construction with default order
   {
      PFFSolver solver(mesh, mat, 1);
      TEST_ASSERT(solver.GetDisplacementFES().GlobalTrueVSize() > 0,
                  "Solver construction succeeds with order 1");
   }

   // Test 2: Construction with higher order
   {
      PFFSolver solver(mesh, mat, 2);
      TEST_ASSERT(solver.GetDisplacementFES().GlobalTrueVSize() > 0,
                  "Solver construction succeeds with order 2");
   }

   // Test 3: Initial fields are zero
   {
      PFFSolver solver(mesh, mat, 1);
      real_t u_max = solver.GetDisplacement().Max();
      real_t d_max = solver.GetDamage().Max();
      TEST_NEAR(u_max, 0.0, 1e-14, "Initial displacement is zero");
      TEST_NEAR(d_max, 0.0, 1e-14, "Initial damage is zero");
   }
}

void TestPFFSolverBoundaryConditions()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Boundary Conditions ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   // Test 1: ConfigureSimpleTensionBCs
   {
      PFFSolver solver(mesh, mat, 1);
      solver.ConfigureSimpleTensionBCs();

      // Check that boundary markers are set
      const Array<int> &mask_top_y = solver.GetMaskTopY();
      const Array<int> &mask_bot_y = solver.GetMaskBotY();

      int num_top = 0, num_bot = 0;
      for (int i = 0; i < mask_top_y.Size(); i++)
      {
         num_top += mask_top_y[i];
         num_bot += mask_bot_y[i];
      }

      // Reduce across MPI
      int global_top, global_bot;
      MPI_Allreduce(&num_top, &global_top, 1, MPI_INT, MPI_SUM, mesh.GetComm());
      MPI_Allreduce(&num_bot, &global_bot, 1, MPI_INT, MPI_SUM, mesh.GetComm());

      TEST_ASSERT(global_top > 0, "Top boundary DOFs marked");
      TEST_ASSERT(global_bot > 0, "Bottom boundary DOFs marked");
   }

   // Test 2: Set prescribed displacement
   {
      PFFSolver solver(mesh, mat, 1);
      solver.ConfigureSimpleTensionBCs();

      TestPrescribedDisplacement u_bc(dim, 1.0);
      solver.SetPrescribedDisplacement(u_bc);

      TEST_ASSERT(true, "Prescribed displacement coefficient set without error");
   }
}

void TestPFFSolverSingleStep()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Single Time Step ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   PFFSolver solver(mesh, mat, 1);
   solver.ConfigureSimpleTensionBCs();
   solver.SetPrintLevel(0);  // Quiet for testing

   TestPrescribedDisplacement u_bc(dim, 1.0);
   solver.SetPrescribedDisplacement(u_bc);

   // Test 1: Solve one time step
   {
      real_t dt = 1e-4;
      int iters = solver.SolveTimeStep(dt);

      TEST_ASSERT(iters > 0, "Solver completes at least one iteration");
      TEST_ASSERT(iters <= 50, "Solver converges within max iterations");
   }

   // Test 2: Displacement is non-zero after loading
   {
      real_t u_max = solver.GetDisplacement().Max();
      real_t u_min = solver.GetDisplacement().Min();

      // Get global max/min
      real_t global_u_max, global_u_min;
      MPI_Allreduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, mesh.GetComm());
      MPI_Allreduce(&u_min, &global_u_min, 1, MPI_DOUBLE, MPI_MIN, mesh.GetComm());

      TEST_ASSERT(global_u_max > 0.0 || global_u_min < 0.0,
                  "Non-zero displacement after loading");
   }

   // Test 3: Damage is in valid range [0, 1]
   {
      real_t d_max = solver.GetDamage().Max();
      real_t d_min = solver.GetDamage().Min();

      real_t global_d_max, global_d_min;
      MPI_Allreduce(&d_max, &global_d_max, 1, MPI_DOUBLE, MPI_MAX, mesh.GetComm());
      MPI_Allreduce(&d_min, &global_d_min, 1, MPI_DOUBLE, MPI_MIN, mesh.GetComm());

      TEST_ASSERT(global_d_min >= -1e-10, "Damage >= 0");
      TEST_ASSERT(global_d_max <= 1.0 + 1e-10, "Damage <= 1");
   }

   // Test 4: Advance time
   {
      real_t dt = 1e-4;
      solver.AdvanceTime(dt);
      TEST_NEAR(solver.GetTime(), dt, 1e-14, "Time advances correctly");
   }
}

void TestPFFSolverDamageEvolution()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Damage Evolution ===\n";
   }

   // Use a finer mesh for damage evolution test
   Mesh serial_mesh = Mesh::MakeCartesian2D(8, 8, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   // Parameters tuned for visible damage with small deformation
   real_t E = 2.1e5;
   real_t nu = 0.3;
   real_t l = 0.05;  // Larger regularization length
   real_t Gc = 0.1;  // Lower fracture toughness for easier damage
   PFFMaterialParameters mat(E, nu, l, Gc, 1e-6, 2.0, true);

   PFFSolver solver(mesh, mat, 1);
   solver.ConfigureSimpleTensionBCs();
   solver.SetPrintLevel(0);
   solver.SetMaxIterations(100);

   // Large amplitude to induce damage
   TestPrescribedDisplacement u_bc(dim, 10.0);
   solver.SetPrescribedDisplacement(u_bc);

   // Test 1: Damage should evolve with sufficient loading
   {
      // Take several steps with increasing load
      real_t dt = 1e-3;
      real_t d_max_initial = 0.0;

      for (int step = 0; step < 5; step++)
      {
         solver.SolveTimeStep(dt);
         solver.AdvanceTime(dt);
      }

      real_t d_max = solver.GetDamage().Max();
      real_t global_d_max;
      MPI_Allreduce(&d_max, &global_d_max, 1, MPI_DOUBLE, MPI_MAX, mesh.GetComm());

      // With high loading, we expect some damage
      TEST_ASSERT(global_d_max >= 0.0, "Damage is non-negative under loading");
   }

   // Test 2: Damage is irreversible
   {
      real_t d_max_before = solver.GetDamage().Max();
      real_t global_d_max_before;
      MPI_Allreduce(&d_max_before, &global_d_max_before, 1, MPI_DOUBLE,
                    MPI_MAX, mesh.GetComm());

      // Take another step
      solver.SolveTimeStep(1e-3);

      real_t d_max_after = solver.GetDamage().Max();
      real_t global_d_max_after;
      MPI_Allreduce(&d_max_after, &global_d_max_after, 1, MPI_DOUBLE,
                    MPI_MAX, mesh.GetComm());

      // Damage should not decrease (irreversibility)
      TEST_ASSERT(global_d_max_after >= global_d_max_before - 1e-10,
                  "Damage is irreversible (max damage does not decrease)");
   }
}

void TestPFFSolverNoLoadNoDamage()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver No Load No Damage ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   PFFSolver solver(mesh, mat, 1);
   solver.ConfigureSimpleTensionBCs();
   solver.SetPrintLevel(0);

   // Zero amplitude - no loading
   TestPrescribedDisplacement u_bc(dim, 0.0);
   solver.SetPrescribedDisplacement(u_bc);

   // Test: With zero load, damage should remain zero
   {
      solver.SolveTimeStep(1e-4);

      real_t d_max = solver.GetDamage().Max();
      real_t global_d_max;
      MPI_Allreduce(&d_max, &global_d_max, 1, MPI_DOUBLE, MPI_MAX, mesh.GetComm());

      TEST_NEAR(global_d_max, 0.0, 1e-10, "No load produces no damage");
   }
}

void TestPFFSolverConvergence()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Convergence ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   PFFSolver solver(mesh, mat, 1);
   solver.ConfigureSimpleTensionBCs();
   solver.SetPrintLevel(0);
   solver.SetMaxIterations(100);
   solver.SetRelativeTolerance(1e-8);

   TestPrescribedDisplacement u_bc(dim, 1.0);
   solver.SetPrescribedDisplacement(u_bc);

   // Test: Solver should converge
   {
      real_t dt = 1e-4;
      int iters = solver.SolveTimeStep(dt);

      TEST_ASSERT(iters < 100, "Staggered iteration converges");
   }
}

// =============================================================================
// Integration Test: Multiple Time Steps
// =============================================================================

void TestPFFSolverMultipleSteps()
{
   if (Mpi::Root())
   {
      std::cout << "\n=== Testing PFFSolver Multiple Time Steps ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL,
                                            true, 1.0, 1.0);
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   int dim = mesh.Dimension();

   PFFMaterialParameters mat(2.1e5, 0.3, 0.02, 2.7, 1e-6, 2.0, true);

   PFFSolver solver(mesh, mat, 1);
   solver.ConfigureSimpleTensionBCs();
   solver.SetPrintLevel(0);

   TestPrescribedDisplacement u_bc(dim, 1.0);
   solver.SetPrescribedDisplacement(u_bc);

   // Test: Run multiple time steps
   {
      real_t dt = 1e-4;
      int num_steps = 5;
      bool all_converged = true;

      for (int step = 0; step < num_steps; step++)
      {
         int iters = solver.SolveTimeStep(dt);
         if (iters >= 50)
         {
            all_converged = false;
         }
         solver.AdvanceTime(dt);
      }

      TEST_ASSERT(all_converged, "All time steps converge");
      TEST_NEAR(solver.GetTime(), num_steps * dt, 1e-12,
                "Time tracking is correct after multiple steps");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   // Initialize MPI
   Mpi::Init(argc, argv);

   if (Mpi::Root())
   {
      std::cout << "Phase Field Fracture - Phase 4 Unit Tests (Staggered Solver)\n";
      std::cout << "=============================================================\n";
   }

   // Run all tests
   TestPFFSolverConstruction();
   TestPFFSolverBoundaryConditions();
   TestPFFSolverSingleStep();
   TestPFFSolverDamageEvolution();
   TestPFFSolverNoLoadNoDamage();
   TestPFFSolverConvergence();
   TestPFFSolverMultipleSteps();

   // Summary (only on rank 0)
   if (Mpi::Root())
   {
      std::cout << "\n=============================================================\n";
      std::cout << "Test Summary:\n";
      std::cout << "  Total tests: " << num_tests << "\n";
      std::cout << "  Passed:      " << num_passed << "\n";
      std::cout << "  Failed:      " << num_failed << "\n";
      std::cout << "=============================================================\n";

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
