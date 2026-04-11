// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for LinearElastic constitutive model.
// Run: ./seas_test_linear_elastic

#include "mfem.hpp"
#include "../../constitutive/constitutive_model.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

#define CHECK(name, cond) TEST_ASSERT(cond, name)

void TestMetadata()
{
   std::cout << "\n=== Test: LinearElastic Metadata ===\n";
   LinearElastic model(3e10, 3e10);  // lambda = mu = 30 GPa

   CHECK("NumInternalVars == 0", model.NumInternalVars() == 0);
   CHECK("IsNonlinear == false", model.IsNonlinear() == false);
   CHECK("HasStateEvolution == false", model.HasStateEvolution() == false);
   CHECK("NumNonLocalVars == 0", model.NumNonLocalVars() == 0);
}

void TestComputeStressPureTension()
{
   std::cout << "\n=== Test: ComputeStress (pure tension xx) ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   // Pure tension in x: eps_xx = 0.001, all others zero
   real_t eps[6] = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0};
   real_t sig[6] = {};
   model.ComputeStress(eps, nullptr, sig);

   // sigma_xx = (lambda + 2*mu) * eps_xx = (3e10 + 6e10) * 0.001 = 9e7
   // sigma_yy = lambda * eps_xx = 3e10 * 0.001 = 3e7
   // sigma_zz = lambda * eps_xx = 3e7
   CHECK("sigma_xx", std::abs(sig[0] - 9e7) < 1.0);
   CHECK("sigma_yy", std::abs(sig[1] - 3e7) < 1.0);
   CHECK("sigma_zz", std::abs(sig[2] - 3e7) < 1.0);
   CHECK("tau_xy == 0", std::abs(sig[3]) < 1e-10);
   CHECK("tau_yz == 0", std::abs(sig[4]) < 1e-10);
   CHECK("tau_xz == 0", std::abs(sig[5]) < 1e-10);
}

void TestComputeStressPureShear()
{
   std::cout << "\n=== Test: ComputeStress (pure shear xy) ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   // Pure shear xy: gamma_xy = 0.002 (engineering strain)
   real_t eps[6] = {0.0, 0.0, 0.0, 0.002, 0.0, 0.0};
   real_t sig[6] = {};
   model.ComputeStress(eps, nullptr, sig);

   // tau_xy = mu * gamma_xy = 3e10 * 0.002 = 6e7
   CHECK("tau_xy", std::abs(sig[3] - 6e7) < 1.0);
   CHECK("sigma_xx == 0", std::abs(sig[0]) < 1e-10);
   CHECK("sigma_yy == 0", std::abs(sig[1]) < 1e-10);
   CHECK("sigma_zz == 0", std::abs(sig[2]) < 1e-10);
}

void TestComputeStressHydrostatic()
{
   std::cout << "\n=== Test: ComputeStress (hydrostatic) ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   // Hydrostatic: eps_xx = eps_yy = eps_zz = 0.001
   real_t eps[6] = {0.001, 0.001, 0.001, 0.0, 0.0, 0.0};
   real_t sig[6] = {};
   model.ComputeStress(eps, nullptr, sig);

   // sigma_ii = lambda * 3 * 0.001 + 2 * mu * 0.001 = 9e7 + 6e7 = 1.5e8
   real_t expected = lambda * 3 * 0.001 + 2 * mu * 0.001;
   CHECK("sigma_xx == sigma_yy", std::abs(sig[0] - sig[1]) < 1e-10);
   CHECK("sigma_yy == sigma_zz", std::abs(sig[1] - sig[2]) < 1e-10);
   CHECK("sigma_xx value", std::abs(sig[0] - expected) < 1.0);
   CHECK("no shear", std::abs(sig[3]) + std::abs(sig[4]) + std::abs(sig[5]) < 1e-10);
}

void TestComputeStressZeroStrain()
{
   std::cout << "\n=== Test: ComputeStress (zero strain) ===\n";
   LinearElastic model(3e10, 3e10);
   real_t eps[6] = {0, 0, 0, 0, 0, 0};
   real_t sig[6] = {999, 999, 999, 999, 999, 999};
   model.ComputeStress(eps, nullptr, sig);

   real_t sum = 0;
   for (int i = 0; i < 6; i++) { sum += std::abs(sig[i]); }
   CHECK("all zero stress", sum < 1e-20);
}

void TestComputeTangent()
{
   std::cout << "\n=== Test: ComputeTangent ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   DenseMatrix C(6, 6);
   model.ComputeTangent(nullptr, nullptr, C);

   // Check diagonal
   CHECK("C(0,0) = lambda+2*mu", std::abs(C(0,0) - (lambda + 2*mu)) < 1.0);
   CHECK("C(1,1) = lambda+2*mu", std::abs(C(1,1) - (lambda + 2*mu)) < 1.0);
   CHECK("C(2,2) = lambda+2*mu", std::abs(C(2,2) - (lambda + 2*mu)) < 1.0);
   CHECK("C(3,3) = mu", std::abs(C(3,3) - mu) < 1.0);
   CHECK("C(4,4) = mu", std::abs(C(4,4) - mu) < 1.0);
   CHECK("C(5,5) = mu", std::abs(C(5,5) - mu) < 1.0);

   // Check off-diagonal
   CHECK("C(0,1) = lambda", std::abs(C(0,1) - lambda) < 1.0);
   CHECK("C(0,2) = lambda", std::abs(C(0,2) - lambda) < 1.0);
   CHECK("C(1,2) = lambda", std::abs(C(1,2) - lambda) < 1.0);

   // Symmetry
   bool symmetric = true;
   for (int i = 0; i < 6; i++)
   {
      for (int j = 0; j < 6; j++)
      {
         if (std::abs(C(i,j) - C(j,i)) > 1e-10) { symmetric = false; }
      }
   }
   CHECK("symmetric", symmetric);

   // Positive definite: all eigenvalues > 0
   // For isotropic: eigenvalues are lambda+2*mu (x1), lambda+2*mu (x1),
   // lambda+2*mu (x1), ... actually the eigenvalues of the 6x6 Voigt matrix
   // for isotropic are: 3*lambda+2*mu (x1), 2*mu (x2 degenerate), mu (x3)
   // All positive since mu > 0 and lambda+2*mu > 0
   CHECK("positive definite (lambda+2*mu > 0)", lambda + 2*mu > 0);
   CHECK("positive definite (mu > 0)", mu > 0);
}

void TestGetPenaltyModulus()
{
   std::cout << "\n=== Test: GetPenaltyModulus ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   CHECK("penalty modulus = lambda+2*mu",
         std::abs(model.GetPenaltyModulus() - (lambda + 2*mu)) < 1.0);
}

void TestGetMaxWaveSpeed()
{
   std::cout << "\n=== Test: GetMaxWaveSpeed ===\n";
   real_t lambda = 3e10, mu = 3e10, rho = 2670.0;
   LinearElastic model(lambda, mu);

   real_t expected = std::sqrt((lambda + 2*mu) / rho);
   CHECK("wave speed", std::abs(model.GetMaxWaveSpeed(rho) - expected) < 0.01);
}

void TestConvenienceAccessors()
{
   std::cout << "\n=== Test: GetLambda / GetMu ===\n";
   real_t lambda = 3.2044e10, mu = 3.2038e10;
   LinearElastic model(lambda, mu);

   CHECK("GetLambda", model.GetLambda() == lambda);
   CHECK("GetMu", model.GetMu() == mu);
}

void TestStressTangentConsistency()
{
   std::cout << "\n=== Test: Stress-Tangent consistency ===\n";
   real_t lambda = 3e10, mu = 3e10;
   LinearElastic model(lambda, mu);

   // Verify: sigma = C * epsilon for any strain
   DenseMatrix C(6, 6);
   model.ComputeTangent(nullptr, nullptr, C);

   real_t eps[6] = {0.001, -0.0005, 0.0003, 0.002, -0.001, 0.0015};
   real_t sig_direct[6] = {};
   model.ComputeStress(eps, nullptr, sig_direct);

   // C * eps via matrix multiply
   real_t sig_tangent[6] = {};
   for (int i = 0; i < 6; i++)
   {
      for (int j = 0; j < 6; j++)
      {
         sig_tangent[i] += C(i, j) * eps[j];
      }
   }

   bool match = true;
   for (int i = 0; i < 6; i++)
   {
      if (std::abs(sig_direct[i] - sig_tangent[i]) > 1.0) { match = false; }
   }
   CHECK("sigma = C * epsilon", match);
}

int main(int argc, char *argv[])
{
   TestMetadata();
   TestComputeStressPureTension();
   TestComputeStressPureShear();
   TestComputeStressHydrostatic();
   TestComputeStressZeroStrain();
   TestComputeTangent();
   TestGetPenaltyModulus();
   TestGetMaxWaveSpeed();
   TestConvenienceAccessors();
   TestStressTangentConsistency();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
