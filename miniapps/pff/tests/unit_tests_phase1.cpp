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

// Unit tests for Phase 1: Material Models and Utilities
// - DegradationFunction
// - SpectralDecomposition
// - PFFMaterialParameters

#include "mfem.hpp"
#include "../materials/degradation_function.hpp"
#include "../materials/spectral_decomposition.hpp"
#include "../materials/pff_material.hpp"

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
// DegradationFunction Tests
// =============================================================================

void TestDegradationFunction()
{
   std::cout << "\n=== Testing DegradationFunction ===\n";

   // Test 1: Default constructor
   {
      DegradationFunction g;
      TEST_NEAR(g.GetEta(), 1e-6, 1e-12, "Default eta");
      TEST_NEAR(g.GetP(), 2.0, 1e-12, "Default p");
   }

   // Test 2: g(0) = 1
   {
      DegradationFunction g;
      real_t val = g.Eval(0.0);
      TEST_NEAR(val, 1.0, 1e-10, "g(0) should be 1");
   }

   // Test 3: g(1) = eta
   {
      real_t eta = 1e-6;
      DegradationFunction g(eta);
      real_t val = g.Eval(1.0);
      TEST_NEAR(val, eta, 1e-12, "g(1) should be eta");
   }

   // Test 4: g'(0) = -2*(1-eta) for p=2
   {
      real_t eta = 1e-6;
      DegradationFunction g(eta, 2.0);
      real_t dg = g.EvalDerivative(0.0);
      real_t expected = -2.0 * (1.0 - eta);
      TEST_NEAR(dg, expected, 1e-10, "g'(0) for p=2");
   }

   // Test 5: g is monotonically decreasing
   {
      DegradationFunction g;
      bool monotonic = true;
      for (real_t d = 0.0; d < 1.0; d += 0.1)
      {
         if (g.EvalDerivative(d) >= 0.0)
         {
            monotonic = false;
            break;
         }
      }
      TEST_ASSERT(monotonic, "g should be monotonically decreasing");
   }

   // Test 6: g''(d) > 0 for d in (0, 1) (convexity)
   {
      DegradationFunction g;
      bool convex = true;
      for (real_t d = 0.01; d < 0.99; d += 0.1)
      {
         if (g.EvalSecondDerivative(d) <= 0.0)
         {
            convex = false;
            break;
         }
      }
      TEST_ASSERT(convex, "g should be convex (g'' > 0)");
   }

   // Test 7: Different exponents
   {
      DegradationFunction g3(1e-6, 3.0);
      real_t val = g3.Eval(0.5);
      real_t expected = std::pow(0.5, 3.0) * (1.0 - 1e-6) + 1e-6;
      TEST_NEAR(val, expected, 1e-10, "g(0.5) for p=3");
   }

   // Test 8: Verify numerical derivative
   {
      DegradationFunction g;
      real_t d = 0.3;
      real_t h = 1e-8;
      real_t numerical_deriv = (g.Eval(d + h) - g.Eval(d - h)) / (2.0 * h);
      real_t analytical_deriv = g.EvalDerivative(d);
      TEST_NEAR(numerical_deriv, analytical_deriv, 1e-6,
                "Numerical vs analytical derivative");
   }
}

// =============================================================================
// CrackGeometricFunction Tests
// =============================================================================

void TestCrackGeometricFunction()
{
   std::cout << "\n=== Testing CrackGeometricFunction ===\n";

   // Test 1: alpha(0) = 0
   {
      real_t val = CrackGeometricFunction::Eval(0.0);
      TEST_NEAR(val, 0.0, 1e-12, "alpha(0) should be 0");
   }

   // Test 2: alpha(1) = 1
   {
      real_t val = CrackGeometricFunction::Eval(1.0);
      TEST_NEAR(val, 1.0, 1e-12, "alpha(1) should be 1");
   }

   // Test 3: alpha'(d) = 2*d
   {
      real_t d = 0.5;
      real_t dval = CrackGeometricFunction::EvalDerivative(d);
      TEST_NEAR(dval, 2.0 * d, 1e-12, "alpha'(0.5) should be 1");
   }

   // Test 4: c0 = 2 for AT2
   {
      real_t c0 = CrackGeometricFunction::GetNormalizationConstant();
      TEST_NEAR(c0, 2.0, 1e-12, "c0 should be 2 for AT2");
   }
}

// =============================================================================
// SpectralDecomposition Tests
// =============================================================================

void TestSpectralDecomposition2D()
{
   std::cout << "\n=== Testing SpectralDecomposition (2D) ===\n";

   // Test 1: Pure tension
   {
      DenseMatrix strain(2, 2);
      strain = 0.0;
      strain(0, 0) = 0.01;  // Tensile in x
      strain(1, 1) = -0.0025;  // Poisson contraction

      DenseMatrix strain_pos(2, 2), strain_neg(2, 2);
      SpectralDecomposition<2>::Decompose(strain, strain_pos, strain_neg);

      // Positive strain should capture the tensile eigenvalue
      real_t psi_pos = SpectralDecomposition<2>::PositiveEnergy(strain, 1.0, 1.0);
      real_t psi_neg = SpectralDecomposition<2>::NegativeEnergy(strain, 1.0, 1.0);

      TEST_ASSERT(psi_pos > 0.0, "Positive energy should be > 0 for tension");
      TEST_ASSERT(psi_neg >= 0.0, "Negative energy should be >= 0");
   }

   // Test 2: Pure compression
   {
      DenseMatrix strain(2, 2);
      strain = 0.0;
      strain(0, 0) = -0.01;  // Compressive
      strain(1, 1) = 0.0025;

      real_t psi_pos = SpectralDecomposition<2>::PositiveEnergy(strain, 1.0, 1.0);
      real_t psi_neg = SpectralDecomposition<2>::NegativeEnergy(strain, 1.0, 1.0);

      // For compression, negative energy should dominate or be comparable
      TEST_ASSERT(psi_neg >= 0.0, "Negative energy should be >= 0 for compression");
   }

   // Test 3: Additivity: epsilon_+ + epsilon_- = epsilon
   {
      DenseMatrix strain(2, 2);
      strain(0, 0) = 0.005;
      strain(0, 1) = 0.002;
      strain(1, 0) = 0.002;
      strain(1, 1) = -0.003;

      DenseMatrix strain_pos(2, 2), strain_neg(2, 2);
      SpectralDecomposition<2>::Decompose(strain, strain_pos, strain_neg);

      for (int i = 0; i < 2; i++)
      {
         for (int j = 0; j < 2; j++)
         {
            real_t sum = strain_pos(i, j) + strain_neg(i, j);
            TEST_NEAR(sum, strain(i, j), 1e-10, "Additivity of decomposition");
         }
      }
   }

   // Test 4: Symmetry preservation
   {
      DenseMatrix strain(2, 2);
      strain(0, 0) = 0.003;
      strain(0, 1) = 0.001;
      strain(1, 0) = 0.001;
      strain(1, 1) = -0.002;

      DenseMatrix strain_pos(2, 2), strain_neg(2, 2);
      SpectralDecomposition<2>::Decompose(strain, strain_pos, strain_neg);

      TEST_NEAR(strain_pos(0, 1), strain_pos(1, 0), 1e-12,
                "strain_pos should be symmetric");
      TEST_NEAR(strain_neg(0, 1), strain_neg(1, 0), 1e-12,
                "strain_neg should be symmetric");
   }

   // Test 5: Energy non-negativity
   {
      // Generate some random symmetric strains
      srand(42);
      for (int k = 0; k < 10; k++)
      {
         DenseMatrix strain(2, 2);
         strain(0, 0) = (rand() / (real_t)RAND_MAX - 0.5) * 0.01;
         strain(1, 1) = (rand() / (real_t)RAND_MAX - 0.5) * 0.01;
         strain(0, 1) = (rand() / (real_t)RAND_MAX - 0.5) * 0.005;
         strain(1, 0) = strain(0, 1);

         real_t psi_pos = SpectralDecomposition<2>::PositiveEnergy(strain, 1.0, 1.0);
         real_t psi_neg = SpectralDecomposition<2>::NegativeEnergy(strain, 1.0, 1.0);

         TEST_ASSERT(psi_pos >= -1e-14, "Positive energy >= 0");
         TEST_ASSERT(psi_neg >= -1e-14, "Negative energy >= 0");
      }
   }
}

void TestSpectralDecomposition3D()
{
   std::cout << "\n=== Testing SpectralDecomposition (3D) ===\n";

   // Test 1: Uniaxial tension
   {
      real_t eps_xx = 0.01;
      real_t nu = 0.25;
      DenseMatrix strain(3, 3);
      strain = 0.0;
      strain(0, 0) = eps_xx;
      strain(1, 1) = -nu * eps_xx;
      strain(2, 2) = -nu * eps_xx;

      DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
      SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

      // The x-direction should be in the positive part
      TEST_ASSERT(strain_pos(0, 0) > 0.0, "Tensile strain in positive part");
   }

   // Test 2: Hydrostatic tension
   {
      real_t eps = 0.005;
      DenseMatrix strain(3, 3);
      strain = 0.0;
      strain(0, 0) = eps;
      strain(1, 1) = eps;
      strain(2, 2) = eps;

      real_t psi_pos = SpectralDecomposition<3>::PositiveEnergy(strain, 1.0, 1.0);
      real_t psi_neg = SpectralDecomposition<3>::NegativeEnergy(strain, 1.0, 1.0);

      TEST_ASSERT(psi_pos > 0.0, "Positive energy > 0 for hydrostatic tension");
      TEST_NEAR(psi_neg, 0.0, 1e-14, "Negative energy = 0 for hydrostatic tension");
   }

   // Test 3: Hydrostatic compression
   {
      real_t eps = -0.005;
      DenseMatrix strain(3, 3);
      strain = 0.0;
      strain(0, 0) = eps;
      strain(1, 1) = eps;
      strain(2, 2) = eps;

      real_t psi_pos = SpectralDecomposition<3>::PositiveEnergy(strain, 1.0, 1.0);
      real_t psi_neg = SpectralDecomposition<3>::NegativeEnergy(strain, 1.0, 1.0);

      TEST_NEAR(psi_pos, 0.0, 1e-14, "Positive energy = 0 for hydrostatic compression");
      TEST_ASSERT(psi_neg > 0.0, "Negative energy > 0 for hydrostatic compression");
   }

   // Test 4: Additivity
   {
      DenseMatrix strain(3, 3);
      strain(0, 0) = 0.003;
      strain(1, 1) = -0.002;
      strain(2, 2) = 0.001;
      strain(0, 1) = strain(1, 0) = 0.001;
      strain(0, 2) = strain(2, 0) = -0.0005;
      strain(1, 2) = strain(2, 1) = 0.0002;

      DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
      SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            real_t sum = strain_pos(i, j) + strain_neg(i, j);
            TEST_NEAR(sum, strain(i, j), 1e-10, "3D additivity");
         }
      }
   }

   // Test 5: Stress computation consistency
   {
      real_t lambda = 1.0, mu = 1.0;
      DenseMatrix strain(3, 3);
      strain = 0.0;
      strain(0, 0) = 0.01;
      strain(1, 1) = -0.003;
      strain(2, 2) = 0.002;

      DenseMatrix stress_pos(3, 3), stress_neg(3, 3);
      SpectralDecomposition<3>::PositiveStress(strain, lambda, mu, stress_pos);
      SpectralDecomposition<3>::NegativeStress(strain, lambda, mu, stress_neg);

      // Stress should also be symmetric
      for (int i = 0; i < 3; i++)
      {
         for (int j = i + 1; j < 3; j++)
         {
            TEST_NEAR(stress_pos(i, j), stress_pos(j, i), 1e-12,
                      "stress_pos symmetry");
            TEST_NEAR(stress_neg(i, j), stress_neg(j, i), 1e-12,
                      "stress_neg symmetry");
         }
      }
   }
}

// =============================================================================
// PFFMaterialParameters Tests
// =============================================================================

void TestPFFMaterialParameters()
{
   std::cout << "\n=== Testing PFFMaterialParameters ===\n";

   // Test 1: Default constructor
   {
      PFFMaterialParameters mat;
      TEST_ASSERT(mat.E > 0.0, "E should be positive");
      TEST_ASSERT(mat.nu > -1.0 && mat.nu < 0.5, "nu should be in valid range");
      TEST_ASSERT(mat.Gc > 0.0, "Gc should be positive");
   }

   // Test 2: Lame parameters from E and nu
   {
      real_t E = 40e6, nu = 0.25;
      PFFMaterialParameters mat(E, nu, 5e-5, 6.43e6);

      real_t expected_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      real_t expected_mu = E / (2.0 * (1.0 + nu));
      real_t expected_K = E / (3.0 * (1.0 - 2.0 * nu));

      TEST_NEAR(mat.lambda, expected_lambda, 1e-6, "Lambda from E, nu");
      TEST_NEAR(mat.mu, expected_mu, 1e-6, "Mu from E, nu");
      TEST_NEAR(mat.K, expected_K, 1e-6, "K from E, nu");
   }

   // Test 3: Fracture toughness from strength
   {
      real_t E = 40e6, l = 5e-5, sigma_t = 6.43e6;
      PFFMaterialParameters mat(E, 0.25, l, sigma_t);

      real_t expected_Gc = 8.0 * l * sigma_t * sigma_t / (3.0 * E);
      TEST_NEAR(mat.Gc, expected_Gc, 1e-10, "Gc from strength formula");
   }

   // Test 4: Explicit Gc constructor
   {
      real_t E = 40e6, nu = 0.25, l = 5e-5, Gc = 100.0;
      PFFMaterialParameters mat(E, nu, l, Gc, 1e-6, 2.0, true);

      TEST_NEAR(mat.Gc, Gc, 1e-10, "Explicit Gc");
   }

   // Test 5: Derived quantities
   {
      PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

      real_t diff_coeff = mat.GetDamageDiffusionCoeff();
      real_t expected_diff = mat.Gc * mat.l / mat.c0;
      TEST_NEAR(diff_coeff, expected_diff, 1e-12, "Diffusion coefficient");

      real_t react_coeff = mat.GetDamageReactionCoeff();
      real_t expected_react = 2.0 * mat.Gc / (mat.c0 * mat.l);
      TEST_NEAR(react_coeff, expected_react, 1e-12, "Reaction coefficient");

      real_t psi_c = mat.GetCriticalStrainEnergy();
      real_t expected_psi_c = mat.Gc / (mat.c0 * mat.l);
      TEST_NEAR(psi_c, expected_psi_c, 1e-12, "Critical strain energy");
   }

   // Test 6: DegradationFunction creation
   {
      PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6, 1e-5, 3.0);
      DegradationFunction g = mat.GetDegradationFunction();

      TEST_NEAR(g.GetEta(), 1e-5, 1e-12, "Degradation function eta");
      TEST_NEAR(g.GetP(), 3.0, 1e-12, "Degradation function p");
   }

   // Test 7: Validation
   {
      PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);
      TEST_ASSERT(mat.Validate(), "Valid parameters should pass validation");
   }

   // Test 8: Relationship between Gc and sigma_t
   {
      // If we set Gc explicitly and then compute sigma_t, we should get
      // consistency
      real_t E = 40e6, l = 5e-5;
      real_t sigma_t_orig = 6.43e6;

      PFFMaterialParameters mat(E, 0.25, l, sigma_t_orig);
      real_t Gc = mat.Gc;

      // Now create with explicit Gc and verify sigma_t is recovered
      PFFMaterialParameters mat2(E, 0.25, l, Gc, 1e-6, 2.0, true);
      TEST_NEAR(mat2.sigma_t, sigma_t_orig, sigma_t_orig * 1e-6,
                "sigma_t recovery from Gc");
   }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "Phase Field Fracture - Phase 1 Unit Tests\n";
   std::cout << "==========================================\n";

   // Run all tests
   TestDegradationFunction();
   TestCrackGeometricFunction();
   TestSpectralDecomposition2D();
   TestSpectralDecomposition3D();
   TestPFFMaterialParameters();

   // Summary
   std::cout << "\n==========================================\n";
   std::cout << "Test Summary:\n";
   std::cout << "  Total tests: " << num_tests << "\n";
   std::cout << "  Passed:      " << num_passed << "\n";
   std::cout << "  Failed:      " << num_failed << "\n";
   std::cout << "==========================================\n";

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
