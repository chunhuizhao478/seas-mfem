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

// Unit tests for Phase 2: Integrators
// - DegradedElasticityIntegrator
// - DamageDiffusionIntegrator
// - DamageSourceIntegrator

#include "mfem.hpp"
#include "../materials/pff_material.hpp"
#include "../materials/degradation_function.hpp"
#include "../materials/spectral_decomposition.hpp"
#include "../integrators/degraded_elasticity_integrator.hpp"
#include "../integrators/damage_diffusion_integrator.hpp"
#include "../integrators/damage_source_integrator.hpp"

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
// DegradedElasticityIntegrator Tests
// =============================================================================

void TestDegradedElasticityIntegrator()
{
   std::cout << "\n=== Testing DegradedElasticityIntegrator ===\n";

   // Create a simple 2D mesh
   Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true, 1.0, 1.0);
   int dim = mesh.Dimension();

   // Finite element spaces
   H1_FECollection u_fec(1, dim);
   FiniteElementSpace u_fes(&mesh, &u_fec, dim);  // Vector field for displacement

   L2_FECollection d_fec(0, dim);
   FiniteElementSpace d_fes(&mesh, &d_fec);  // Scalar field for damage

   // Material parameters
   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Test 1: Zero displacement should give zero residual
   {
      GridFunction damage(&d_fes);
      damage = 0.0;

      GridFunction u(&u_fes);
      u = 0.0;

      NonlinearForm nlf(&u_fes);
      nlf.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, damage));

      Vector residual(u_fes.GetTrueVSize());
      nlf.Mult(u, residual);

      TEST_NEAR(residual.Norml2(), 0.0, 1e-12,
                "Zero displacement gives zero residual");
   }

   // Test 2: Zero damage should give standard linear elasticity behavior
   {
      GridFunction damage(&d_fes);
      damage = 0.0;

      // Create uniform strain state
      GridFunction u(&u_fes);
      u = 0.0;

      // Apply simple displacement
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.001 * x(0);  // Small strain in x
         u(1) = 0.0;
      });
      u.ProjectCoefficient(u_coeff);

      // Energy should be positive
      NonlinearForm nlf(&u_fes);
      nlf.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, damage));

      real_t energy = nlf.GetGridFunctionEnergy(u);
      TEST_ASSERT(energy > 0.0, "Positive strain energy for non-zero displacement");
   }

   // Test 3: Full damage (d=1) should reduce stiffness to eta
   {
      GridFunction damage(&d_fes);
      damage = 1.0;  // Fully damaged

      GridFunction u(&u_fes);
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.001 * x(0);
         u(1) = 0.0;
      });
      u.ProjectCoefficient(u_coeff);

      NonlinearForm nlf_damaged(&u_fes);
      nlf_damaged.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, damage));
      real_t energy_damaged = nlf_damaged.GetGridFunctionEnergy(u);

      // Compare with undamaged
      GridFunction damage_zero(&d_fes);
      damage_zero = 0.0;

      NonlinearForm nlf_intact(&u_fes);
      nlf_intact.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, damage_zero));
      real_t energy_intact = nlf_intact.GetGridFunctionEnergy(u);

      // Damaged energy should be much smaller (scaled by eta for tensile part)
      TEST_ASSERT(energy_damaged < energy_intact,
                  "Full damage reduces strain energy");
      TEST_ASSERT(energy_damaged > 0.0,
                  "Residual stiffness gives positive energy");
   }

   // Test 4: Jacobian should be symmetric
   {
      GridFunction damage(&d_fes);
      damage = 0.3;

      GridFunction u(&u_fes);
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.001 * x(0);
         u(1) = -0.0003 * x(1);
      });
      u.ProjectCoefficient(u_coeff);

      NonlinearForm nlf(&u_fes);
      nlf.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, damage));

      // Get Jacobian
      Operator &K = nlf.GetGradient(u);
      SparseMatrix *K_sp = dynamic_cast<SparseMatrix*>(&K);

      if (K_sp)
      {
         // Check symmetry by comparing K*x with K^T*x for random vectors
         bool symmetric = true;
         real_t sym_tol = 1e-8;
         srand(42);

         for (int test = 0; test < 5 && symmetric; test++)
         {
            Vector x(K_sp->Height()), Kx(K_sp->Height()), Ktx(K_sp->Height());
            for (int i = 0; i < x.Size(); i++)
            {
               x(i) = rand() / (real_t)RAND_MAX - 0.5;
            }

            K_sp->Mult(x, Kx);
            K_sp->MultTranspose(x, Ktx);

            // For symmetric matrix: K*x should equal K^T*x
            real_t diff = 0.0;
            real_t norm = 0.0;
            for (int i = 0; i < Kx.Size(); i++)
            {
               diff += (Kx(i) - Ktx(i)) * (Kx(i) - Ktx(i));
               norm += Kx(i) * Kx(i) + Ktx(i) * Ktx(i);
            }
            if (std::sqrt(diff) > sym_tol * std::sqrt(norm) + 1e-14)
            {
               symmetric = false;
            }
         }
         TEST_ASSERT(symmetric, "Jacobian matrix is symmetric");
      }
   }

   // Test 5: Positive strain energy at quadrature points
   {
      GridFunction damage(&d_fes);
      damage = 0.0;

      GridFunction u(&u_fes);
      VectorFunctionCoefficient u_coeff(dim, [](const Vector &x, Vector &u) {
         u(0) = 0.002 * x(0);
         u(1) = 0.001 * x(1);
      });
      u.ProjectCoefficient(u_coeff);

      DegradedElasticityIntegrator integ(mat, damage);

      // Check energy on first element
      const FiniteElement &el = *u_fes.GetFE(0);
      ElementTransformation *Tr = mesh.GetElementTransformation(0);

      Array<int> vdofs;
      u_fes.GetElementVDofs(0, vdofs);
      Vector elfun(vdofs.Size());
      u.GetSubVector(vdofs, elfun);

      Vector psi_active;
      integ.ComputePositiveStrainEnergy(el, *Tr, elfun, psi_active);

      bool all_positive = true;
      for (int i = 0; i < psi_active.Size(); i++)
      {
         if (psi_active(i) < -1e-14)
         {
            all_positive = false;
            break;
         }
      }
      TEST_ASSERT(all_positive, "Positive strain energy at all quadrature points");
   }
}

// =============================================================================
// DamageDiffusionIntegrator Tests
// =============================================================================

void TestDamageDiffusionIntegrator()
{
   std::cout << "\n=== Testing DamageDiffusionIntegrator ===\n";

   // Create mesh and FE space
   Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true, 1.0, 1.0);
   int dim = mesh.Dimension();

   H1_FECollection fec(1, dim);
   FiniteElementSpace fes(&mesh, &fec);

   // Test 1: Coefficient computation
   {
      real_t Gc = 100.0, l = 0.01, c0 = 2.0;
      DamageDiffusionIntegrator integ(Gc, l, c0);
      real_t expected = Gc * l / c0;
      TEST_NEAR(integ.GetCoefficient(), expected, 1e-12, "Coefficient = Gc*l/c0");
   }

   // Test 2: Compare with standard DiffusionIntegrator
   {
      real_t Gc = 100.0, l = 0.01, c0 = 2.0;
      real_t coeff = Gc * l / c0;

      DamageDiffusionIntegrator pff_integ(Gc, l, c0);
      ConstantCoefficient const_coeff(coeff);
      DiffusionIntegrator std_integ(const_coeff);

      // Compare element matrices on first element
      const FiniteElement &el = *fes.GetFE(0);
      ElementTransformation *Tr = mesh.GetElementTransformation(0);

      DenseMatrix elmat_pff, elmat_std;
      pff_integ.AssembleElementMatrix(el, *Tr, elmat_pff);
      std_integ.AssembleElementMatrix(el, *Tr, elmat_std);

      real_t diff = 0.0;
      for (int i = 0; i < elmat_pff.Height(); i++)
      {
         for (int j = 0; j < elmat_pff.Width(); j++)
         {
            diff += std::abs(elmat_pff(i, j) - elmat_std(i, j));
         }
      }
      TEST_NEAR(diff, 0.0, 1e-10, "Matches standard DiffusionIntegrator");
   }

   // Test 3: Matrix is symmetric
   {
      DamageDiffusionIntegrator integ(100.0, 0.01);

      const FiniteElement &el = *fes.GetFE(0);
      ElementTransformation *Tr = mesh.GetElementTransformation(0);

      DenseMatrix elmat;
      integ.AssembleElementMatrix(el, *Tr, elmat);

      bool symmetric = true;
      for (int i = 0; i < elmat.Height(); i++)
      {
         for (int j = i + 1; j < elmat.Width(); j++)
         {
            if (std::abs(elmat(i, j) - elmat(j, i)) > 1e-12)
            {
               symmetric = false;
               break;
            }
         }
         if (!symmetric) break;
      }
      TEST_ASSERT(symmetric, "Element matrix is symmetric");
   }

   // Test 4: Matrix is positive semi-definite
   {
      DamageDiffusionIntegrator integ(100.0, 0.01);

      const FiniteElement &el = *fes.GetFE(0);
      ElementTransformation *Tr = mesh.GetElementTransformation(0);

      DenseMatrix elmat;
      integ.AssembleElementMatrix(el, *Tr, elmat);

      // Check that x^T K x >= 0 for random vectors
      bool psd = true;
      srand(42);
      for (int k = 0; k < 10; k++)
      {
         Vector x(elmat.Height());
         for (int i = 0; i < x.Size(); i++)
         {
            x(i) = rand() / (real_t)RAND_MAX - 0.5;
         }
         Vector Kx(elmat.Height());
         elmat.Mult(x, Kx);
         real_t xKx = x * Kx;
         if (xKx < -1e-12)
         {
            psd = false;
            break;
         }
      }
      TEST_ASSERT(psd, "Element matrix is positive semi-definite");
   }
}

// =============================================================================
// DamageSourceIntegrator Tests
// =============================================================================

void TestDamageSourceIntegrator()
{
   std::cout << "\n=== Testing DamageSourceIntegrator ===\n";

   // Create mesh and FE spaces
   Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL, true, 1.0, 1.0);
   int dim = mesh.Dimension();

   H1_FECollection fec(1, dim);
   FiniteElementSpace fes(&mesh, &fec);

   L2_FECollection l2_fec(0, dim);
   FiniteElementSpace l2_fes(&mesh, &l2_fec);

   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Test 1: Zero damage and zero history gives zero residual
   {
      GridFunction H(&l2_fes);
      H = 0.0;

      GridFunction d(&fes);
      d = 0.0;

      NonlinearForm nlf(&fes);
      nlf.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      Vector residual(fes.GetTrueVSize());
      nlf.Mult(d, residual);

      TEST_NEAR(residual.Norml2(), 0.0, 1e-12,
                "Zero damage and history gives zero residual");
   }

   // Test 2: Non-zero damage with zero history should give positive residual
   // (because of the 2Gc/(c0l)*d term)
   {
      GridFunction H(&l2_fes);
      H = 0.0;

      GridFunction d(&fes);
      d = 0.5;

      NonlinearForm nlf(&fes);
      nlf.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      Vector residual(fes.GetTrueVSize());
      nlf.Mult(d, residual);

      // Residual should be positive (driving damage back to zero)
      real_t sum = 0.0;
      for (int i = 0; i < residual.Size(); i++)
      {
         sum += residual(i);
      }
      TEST_ASSERT(sum > 0.0, "Positive damage with zero H gives positive residual");
   }

   // Test 3: High strain energy history should drive damage
   {
      GridFunction H(&l2_fes);
      H = mat.GetCriticalStrainEnergy() * 10.0;  // Much higher than critical

      GridFunction d(&fes);
      d = 0.0;

      NonlinearForm nlf(&fes);
      nlf.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      Vector residual(fes.GetTrueVSize());
      nlf.Mult(d, residual);

      // g'(0) = -2*(1-eta) < 0, so g'(0)*H < 0
      // Residual should be negative (driving damage up)
      real_t sum = 0.0;
      for (int i = 0; i < residual.Size(); i++)
      {
         sum += residual(i);
      }
      TEST_ASSERT(sum < 0.0, "High H drives damage evolution (negative residual)");
   }

   // Test 4: Jacobian correctness (finite difference check using directional derivative)
   {
      GridFunction H(&l2_fes);
      H = 1000.0;

      GridFunction d(&fes);
      d = 0.3;

      NonlinearForm nlf(&fes);
      nlf.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      // Get analytical Jacobian
      Operator &K = nlf.GetGradient(d);

      // Compute finite difference vs analytical directional derivative
      // For direction v: (R(d + eps*v) - R(d)) / eps ≈ K * v
      real_t eps = 1e-7;
      Vector d_vec;
      d.GetTrueDofs(d_vec);

      Vector r0(fes.GetTrueVSize());
      nlf.Mult(d, r0);

      bool fd_match = true;
      srand(123);

      for (int test = 0; test < 3 && fd_match; test++)
      {
         // Random direction vector
         Vector v(d_vec.Size());
         for (int i = 0; i < v.Size(); i++)
         {
            v(i) = rand() / (real_t)RAND_MAX - 0.5;
         }

         // Finite difference: (R(d + eps*v) - R(d)) / eps
         Vector d_pert(d_vec);
         d_pert.Add(eps, v);
         GridFunction d_gf(&fes);
         d_gf.SetFromTrueDofs(d_pert);

         Vector r1(fes.GetTrueVSize());
         nlf.Mult(d_gf, r1);

         Vector fd_deriv(r0.Size());
         for (int i = 0; i < r0.Size(); i++)
         {
            fd_deriv(i) = (r1(i) - r0(i)) / eps;
         }

         // Analytical: K * v
         Vector Kv(r0.Size());
         K.Mult(v, Kv);

         // Compare
         real_t diff = 0.0, norm = 0.0;
         for (int i = 0; i < r0.Size(); i++)
         {
            diff += (fd_deriv(i) - Kv(i)) * (fd_deriv(i) - Kv(i));
            norm += Kv(i) * Kv(i) + fd_deriv(i) * fd_deriv(i);
         }
         if (std::sqrt(diff) > 1e-4 * std::sqrt(norm) + 1e-10)
         {
            fd_match = false;
         }
      }
      TEST_ASSERT(fd_match, "Jacobian matches finite difference");
   }

   // Test 5: Energy is consistent
   {
      GridFunction H(&l2_fes);
      H = 500.0;

      GridFunction d(&fes);
      d = 0.4;

      NonlinearForm nlf(&fes);
      nlf.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      real_t energy = nlf.GetGridFunctionEnergy(d);
      TEST_ASSERT(energy > 0.0, "Energy is positive");
   }
}

// =============================================================================
// Integration Test: Combined Integrators
// =============================================================================

void TestIntegratorsIntegration()
{
   std::cout << "\n=== Testing Integrators Integration ===\n";

   // Create mesh
   Mesh mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL, true, 1.0, 1.0);
   int dim = mesh.Dimension();

   // FE spaces
   H1_FECollection u_fec(1, dim);
   FiniteElementSpace u_fes(&mesh, &u_fec, dim);

   H1_FECollection d_fec(1, dim);
   FiniteElementSpace d_fes(&mesh, &d_fec);

   L2_FECollection H_fec(0, dim);
   FiniteElementSpace H_fes(&mesh, &H_fec);

   PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

   // Test: Full system setup doesn't crash
   {
      GridFunction d(&d_fes);
      d = 0.0;

      GridFunction H(&H_fes);
      H = 0.0;

      // Elasticity form
      NonlinearForm elasticity(&u_fes);
      elasticity.AddDomainIntegrator(new DegradedElasticityIntegrator(mat, d));

      // Damage form (diffusion + source)
      BilinearForm damage_diff(&d_fes);
      damage_diff.AddDomainIntegrator(new DamageDiffusionIntegrator(
         mat.Gc, mat.l, mat.c0));

      NonlinearForm damage_src(&d_fes);
      damage_src.AddDomainIntegrator(new DamageSourceIntegrator(mat, H));

      // Assemble
      damage_diff.Assemble();

      GridFunction u(&u_fes);
      u = 0.0;

      Vector r_u(u_fes.GetTrueVSize());
      elasticity.Mult(u, r_u);

      Vector r_d(d_fes.GetTrueVSize());
      damage_src.Mult(d, r_d);

      TEST_ASSERT(true, "Full system setup and assembly works");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "Phase Field Fracture - Phase 2 Unit Tests (Integrators)\n";
   std::cout << "========================================================\n";

   // Run all tests
   TestDegradedElasticityIntegrator();
   TestDamageDiffusionIntegrator();
   TestDamageSourceIntegrator();
   TestIntegratorsIntegration();

   // Summary
   std::cout << "\n========================================================\n";
   std::cout << "Test Summary:\n";
   std::cout << "  Total tests: " << num_tests << "\n";
   std::cout << "  Passed:      " << num_passed << "\n";
   std::cout << "  Failed:      " << num_failed << "\n";
   std::cout << "========================================================\n";

   if (num_failed > 0)
   {
      std::cout << "SOME TESTS FAILED!\n";
      return 1;
   }
   else
   {
      std::cout << "ALL TESTS PASSED!\n";
      return 0;
   }
}
