// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests of the tension/compression split (constitutive/strain_split.hpp):
//  - psi+ + psi- = psi0 and sigma(g = 1) = sigma0 for every decomposition
//  - sigma = d psi / d eps (central differences at fixed g) on random strains
//  - limiting states: pure tension, pure compression, repeated eigenvalues
//  - plane-strain padding: lambda = K - 2G/3 as with LIBMESH_DIM = 3 in MOOSE

#include "mfem.hpp"
#include "../../constitutive/strain_split.hpp"
#include "../test_framework.hpp"

#include <random>

using namespace mfem;
using namespace mfem::pff;

static const real_t E = 2.80e9, nu = 0.35;
static const real_t K = E / 3.0 / (1.0 - 2.0 * nu), G = E / 2.0 / (1.0 + nu);

static void RandomSymmetric(std::mt19937 &rng, real_t scale, real_t *e)
{
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   for (int j = 0; j < 3; j++)
      for (int i = 0; i <= j; i++) { e[i + 3*j] = e[j + 3*i] = scale * U(rng); }
}

static real_t Psi(Decomposition dec, const real_t *e, real_t g)
{
   SplitResult r; ComputeSplit(dec, e, K, G, g, r); return r.psi;
}

static const char *Name(Decomposition d)
{
   return d == Decomposition::None ? "none" : d == Decomposition::Spectral ? "spectral" : "voldev";
}

static void TestConsistency()
{
   std::mt19937 rng(1234);
   for (Decomposition dec : {Decomposition::None, Decomposition::Spectral, Decomposition::VolDev})
   {
      std::cout << "\n[" << Name(dec) << ": energy split and stress = d psi/d eps]\n";
      real_t max_split_err = 0.0, max_g1_err = 0.0, max_fd_err = 0.0;
      for (int trial = 0; trial < 1000; trial++)
      {
         real_t e[9]; RandomSymmetric(rng, 1e-3, e);
         const real_t tr = e[0] + e[4] + e[8];
         real_t ee = 0.0; for (int i = 0; i < 9; i++) { ee += e[i]*e[i]; }
         const real_t lambda = K - 2.0*G/3.0;
         const real_t psi0 = 0.5*lambda*tr*tr + G*ee;

         SplitResult r1; ComputeSplit(dec, e, K, G, 1.0, r1);
         const real_t scale = std::abs(psi0) + 1.0;
         if (dec != Decomposition::None)
         {
            max_split_err = std::max(max_split_err,
                                     std::abs(r1.psi_active + r1.psi_inactive - psi0) / scale);
         }
         else
         {
            max_split_err = std::max(max_split_err, std::abs(r1.psi_active - psi0) / scale);
         }
         for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
            {
               const real_t s0 = lambda*tr*(i == j) + 2.0*G*e[i + 3*j];
               max_g1_err = std::max(max_g1_err, std::abs(r1.sigma[i + 3*j] - s0) / (G * 1e-3));
            }

         const real_t g = 0.37;
         SplitResult r; ComputeSplit(dec, e, K, G, g, r);
         const real_t h = 1e-9;
         for (int j = 0; j < 3; j++)
            for (int i = 0; i <= j; i++)
            {
               // symmetric perturbation of (i,j) and (j,i): d psi = sigma_ij de_ij (+ ji)
               real_t ep[9], em[9];
               for (int k = 0; k < 9; k++) { ep[k] = em[k] = e[k]; }
               ep[i + 3*j] += h; em[i + 3*j] -= h;
               if (i != j) { ep[j + 3*i] += h; em[j + 3*i] -= h; }
               const real_t fd = (Psi(dec, ep, g) - Psi(dec, em, g)) / (2*h);
               const real_t an = (i == j) ? r.sigma[i + 3*j] : 2.0 * r.sigma[i + 3*j];
               max_fd_err = std::max(max_fd_err, std::abs(fd - an) / (G * 1e-3));
            }
      }
      TEST_ASSERT(max_split_err < 1e-12, "psi+ + psi- = psi0 on 1000 random strains");
      TEST_ASSERT(max_g1_err < 1e-12, "g = 1 gives the undamaged stress");
      TEST_ASSERT(max_fd_err < 1e-5, "sigma = d psi/d eps at g = 0.37 (central differences)");
      std::cout << "    max errors: split " << max_split_err << ", g=1 " << max_g1_err
                << ", fd " << max_fd_err << "\n";
   }
}

static void TestLimits()
{
   std::cout << "\n[spectral: limiting states]\n";
   const real_t g = 1e-3;
   // all principal strains positive: everything active, sigma = g sigma0
   real_t et[9] = {2e-4, 0, 0, 0, 1e-4, 0, 0, 0, 3e-5};
   SplitResult r; ComputeSplit(Decomposition::Spectral, et, K, G, g, r);
   TEST_NEAR(r.psi_inactive, 0.0, 1e-9, "pure tension: psi- = 0");
   TEST_REL(r.sigma[0], g * ((K - 2*G/3)*3.3e-4 + 2*G*2e-4), 1e-12, "pure tension: sigma = g sigma0");
   // all negative: nothing degraded
   real_t ec[9] = {-2e-4, 0, 0, 0, -1e-4, 0, 0, 0, -3e-5};
   ComputeSplit(Decomposition::Spectral, ec, K, G, g, r);
   TEST_NEAR(r.psi_active, 0.0, 1e-12, "pure compression: psi+ = 0");
   TEST_REL(r.sigma[0], (K - 2*G/3)*(-3.3e-4) + 2*G*(-2e-4), 1e-12, "pure compression: sigma = sigma0");
   // repeated eigenvalues (hydrostatic, and a double root) stay finite and continuous
   real_t eh[9] = {1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4};
   ComputeSplit(Decomposition::Spectral, eh, K, G, 1.0, r);
   TEST_REL(r.psi_active, 0.5*K*9e-8 + 0.0 + (0.5*(K - 2*G/3)*9e-8 + G*3e-8 - 0.5*K*9e-8), 1e-12,
            "hydrostatic tension: psi+ = psi0");
   real_t e2[9] = {1e-4, 0, 0, 0, 1e-4, 0, 0, 0, -5e-5};
   real_t e2p[9] = {1e-4, 1e-12, 0, 1e-12, 1e-4, 0, 0, 0, -5e-5};
   SplitResult ra, rb;
   ComputeSplit(Decomposition::Spectral, e2, K, G, 0.2, ra);
   ComputeSplit(Decomposition::Spectral, e2p, K, G, 0.2, rb);
   real_t dmax = 0.0; for (int i = 0; i < 9; i++) { dmax = std::max(dmax, std::abs(ra.sigma[i] - rb.sigma[i])); }
   TEST_ASSERT(dmax < 1e-6 * G * 1e-4 * 1e3, "double root: stress continuous under a 1e-12 perturbation");
}

static void TestPlaneStrainPadding()
{
   std::cout << "\n[plane strain: 2x2 strain padded to 3x3]\n";
   const real_t e2[4] = {1e-4, 2e-5, 2e-5, -3e-5};
   real_t e3[9]; PadStrain3(2, e2, e3);
   TEST_NEAR(e3[8], 0.0, 0.0, "eps_zz = 0");
   SplitResult r; ComputeSplit(Decomposition::Spectral, e3, K, G, 1.0, r);
   const real_t lambda = K - 2*G/3, tr = 1e-4 - 3e-5;
   TEST_REL(r.sigma[8], lambda * tr, 1e-12, "sigma_zz = lambda tr(eps) (plane strain, lambda = K - 2G/3)");
}

int main()
{
   TestConsistency();
   TestLimits();
   TestPlaneStrainPadding();
   return TestSummary("test_strain_split");
}
