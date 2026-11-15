// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests: g(d), alpha(d), c0, damage source and the AT1 strength relation
// against the values of the MOOSE inputs (farms branch cdms).

#include "mfem.hpp"
#include "../../phasefield/phase_field_functions.hpp"
#include "../test_framework.hpp"

using namespace mfem;
using namespace mfem::pff;

static void TestDegradation()
{
   std::cout << "\n[degradation g(d) = (1-d)^2 (1-eta) + eta]\n";
   PowerDegradation g;   // p = 2, eta = 1e-6 (all inputs)
   TEST_NEAR(g.g(0.0), 1.0, 1e-15, "g(0) = 1");
   TEST_NEAR(g.g(1.0), 1e-6, 1e-18, "g(1) = eta");
   TEST_NEAR(g.dg(0.3), -2.0 * 0.7 * (1.0 - 1e-6), 1e-14, "g'(d) = -2(1-d)(1-eta)");
   TEST_NEAR(g.d2g(0.3), 2.0 * (1.0 - 1e-6), 1e-14, "g''(d) = 2(1-eta)");
   const real_t h = 1e-6, d = 0.42;
   TEST_NEAR(g.dg(d), (g.g(d + h) - g.g(d - h)) / (2*h), 1e-8, "g' matches central difference");
}

static void TestCrackGeometric()
{
   std::cout << "\n[crack geometric function and c0 = 4 int sqrt(alpha)]\n";
   for (CrackModel m : {CrackModel::AT1, CrackModel::AT2})
   {
      CrackGeometric a; a.model = m;
      // midpoint rule on 1e5 cells
      const int n = 100000; real_t s = 0.0;
      for (int i = 0; i < n; i++) { s += std::sqrt(a.alpha((i + 0.5) / n)); }
      s *= 4.0 / n;
      TEST_NEAR(a.c0(), s, 1e-6, m == CrackModel::AT1 ? "AT1 c0 = 8/3" : "AT2 c0 = 2");
   }
}

static void TestAT1Strength()
{
   std::cout << "\n[AT1 strength relation, spe_hftc values]\n";
   const real_t E = 2.80e9, ft = 69e6, l = 6e-4;
   const real_t Gc = PhaseFieldModel::AT1ToughnessFromStrength(E, ft, l);
   // elasticity_absorbing.i: Gc_const = 8 l ft^2/(3E); fracture.i default 2720.5714285714284
   TEST_NEAR(Gc, 2720.5714285714284, 1e-9, "Gc = 8 l ft^2 / (3E) = 2720.57 J/m^2");
   TEST_NEAR(PhaseFieldModel::AT1Strength(E, Gc, l), ft, 1e-6, "sigma_c = ft");

   PhaseFieldModel pf; pf.Gc = Gc; pf.l = l;
   pf.crack.model = CrackModel::AT1; pf.driving = DrivingEnergy::Active;
   // homogeneous uniaxial stress: psi+ = sigma^2/(2E); the source vanishes at
   // d = 0 when sigma = ft / sqrt(1 - eta)
   const real_t Dc = pf.OnsetDrivingEnergy();
   TEST_REL(Dc, ft * ft / (2.0 * E) / (1.0 - 1e-6), 1e-12, "onset energy = ft^2/(2E(1-eta))");
   TEST_NEAR(pf.Source(0.0, Dc), 0.0, 1e-9 * Dc, "source(0, D_c) = 0");
   TEST_ASSERT(pf.Source(0.0, 0.99 * Dc) > 0.0, "below onset the source is positive (d stays at its bound)");
   TEST_ASSERT(pf.Source(0.0, 1.01 * Dc) < 0.0, "above onset the source is negative (d grows)");
   TEST_NEAR(pf.DiffusionCoefficient(), 2.0 * Gc * l / (8.0/3.0), 1e-12, "diffusion 2 Gc l / c0 (ADPFFDiffusion)");
}

static void TestDrivingConventions()
{
   std::cout << "\n[driving energy conventions]\n";
   PhaseFieldModel spe; spe.Gc = 2720.57; spe.l = 6e-4;
   spe.driving = DrivingEnergy::Active;
   PhaseFieldModel cm = spe; cm.driving = DrivingEnergy::HalfHistoryEnhanced;
   // cmame: H = 2 psi+ (dry, p = 0) and psi uses 0.5 g H, so D equals psi+
   const real_t psi_plus = 3.0e5;
   TEST_NEAR(cm.DrivingFrom(2.0 * psi_plus), spe.DrivingFrom(psi_plus), 1e-9,
             "0.5 * (2 psi+) = psi+: the conventions agree for a dry history");
}

int main()
{
   TestDegradation();
   TestCrackGeometric();
   TestAT1Strength();
   TestDrivingConventions();
   return TestSummary("test_phase_field_functions");
}
