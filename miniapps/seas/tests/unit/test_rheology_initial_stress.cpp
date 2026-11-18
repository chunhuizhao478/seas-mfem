// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Rheology plan rev 1, unit tests U8-U9 (no mesh): the effective initial
// stress of rheology/initial_stress.hpp for TPV26/27.
//
//   U8  sigma_0'(x) against the description's formulas written in the SCEC
//       frame (v13, p. 6-7), the taper, sigma_0' = 0 at the surface, the
//       DRDG3D t = 0 fault stresses, and a general up vector by rotation
//   U9  the initial state is below yield everywhere; sqrt(J2)/Y = 0.930 at
//       10 km (the description's 93 %) and 0.943 at 15 km
//
// Run: ./seas_test_rheology_initial_stress

#include "mfem.hpp"
#include "../../rheology/rheology_model.hpp"
#include "../../rheology/initial_stress.hpp"
#include "../../rheology/drucker_prager.hpp"
#include "test_macros.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace
{

// TPV26/27 constants (description v13, Parts 2 and 4)
const real_t kRho = 2670.0, kRhoF = 1000.0, kG = 9.8;
const real_t kB11 = 0.926793, kB33 = 1.073206, kB13 = -0.169029;

InitialStressSpec Tpv26Spec()
{
   InitialStressSpec s;
   s.kind = InitialStressSpec::Kind::LithostaticRatio;
   s.rho_rock = kRho;
   s.rho_fluid = kRhoF;
   s.gravity = kG;
   s.surface_elevation = 0.0;
   s.taper_top_depth = 15000.0;
   s.taper_bottom_depth = 20000.0;
   // code frame: X strike (SCEC x), Y fault normal (SCEC z), Z up (SCEC -y)
   const real_t B[9] = {kB11, kB13, 0.0, kB13, kB33, 0.0, 0.0, 0.0, 0.0};
   for (int i = 0; i < 9; i++) { s.B[i] = B[i]; }
   return s;
}

/// The description's formulas in the SCEC frame at depth y (total stress).
void ScecFormulas(real_t y, real_t &s11, real_t &s22, real_t &s33, real_t &s13, real_t &pf)
{
   real_t om;
   if (y <= 15000.0) { om = 1.0; }
   else if (y <= 20000.0) { om = (20000.0 - y) / 5000.0; }
   else { om = 0.0; }
   pf = kRhoF * kG * y;
   s22 = -kRho * kG * y;
   s11 = om * (kB11 * (s22 + pf) - pf) + (1.0 - om) * s22;
   s33 = om * (kB33 * (s22 + pf) - pf) + (1.0 - om) * s22;
   s13 = om * (kB13 * (s22 + pf));
}

real_t MaxAbs9(const real_t a[9])
{
   real_t m = 0.0;
   for (int i = 0; i < 9; i++) { m = std::max(m, std::abs(a[i])); }
   return m;
}

void Rotate9(const real_t R[9], const real_t A[9], real_t out[9])
{
   real_t T[9];
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         T[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { T[3 * i + j] += R[3 * i + k] * A[3 * k + j]; }
      }
   }
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         out[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { out[3 * i + j] += T[3 * i + k] * R[3 * j + k]; }
      }
   }
}

} // namespace

// ---------------------------------------------------------------- U8
void TestU8Formulas()
{
   std::cout << "\n=== U8: sigma_0' against the description's formulas ===\n";
   InitialStress is(Tpv26Spec());
   const real_t depths[] = {0.0, 1000.0, 5000.0, 10000.0, 14999.0, 15000.0, 17500.0,
                            20000.0, 25000.0
                           };
   const real_t xy[][2] = {{-5000.0, 0.0}, {12000.0, -350.0}};
   real_t worst = 0.0;
   for (real_t d : depths)
   {
      for (const auto &p : xy)
      {
         const real_t x[3] = {p[0], p[1], -d};
         real_t s[9];
         is.EvalEffective(x, s);
         real_t s11, s22, s33, s13, pf;
         ScecFormulas(d, s11, s22, s33, s13, pf);
         // code frame: XX = 11, YY = 33, ZZ = 22, XY = 13; effective adds pf
         const real_t ref[9] = {s11 + pf, s13, 0.0, s13, s33 + pf, 0.0, 0.0, 0.0, s22 + pf};
         real_t diff[9];
         for (int i = 0; i < 9; i++) { diff[i] = s[i] - ref[i]; }
         worst = std::max(worst, MaxAbs9(diff) / std::max(MaxAbs9(ref), real_t(1.0)));
      }
   }
   TEST_ASSERT(worst <= 1e-12, "U8 lithostatic_ratio = description formulas (1e-12)");

   TEST_ASSERT(is.Taper(15000.0) == 1.0 && is.Taper(17500.0) == 0.5 && is.Taper(20000.0) == 0.0,
               "U8 taper 1, 0.5, 0 at 15, 17.5, 20 km");
   const real_t x0[3] = {3000.0, 100.0, 0.0};
   real_t s0[9];
   is.EvalEffective(x0, s0);
   TEST_ASSERT(MaxAbs9(s0) == 0.0, "U8 sigma_0' = 0 at the surface");
   const real_t xa[3] = {3000.0, 100.0, 1e-9};   // mesh round-off above the surface
   is.EvalEffective(xa, s0);
   TEST_ASSERT(MaxAbs9(s0) == 0.0, "U8 a point 1e-9 m above the surface is at depth 0");
}

void TestU8Drdg3d()
{
   std::cout << "\n=== U8: DRDG3D t = 0 fault stresses ===\n";
   // tpv26/benchmark_data/scec_drdg3d/tpv26_drdg3d_x2_*_x3_*.txt, first row,
   // h-shear-stress and n-stress (MPa; effective, compression negative).
   // 15 km is left out: it is the kink of the taper, where DRDG3D's
   // element interpolation reads 0.4 % low (benchmark_data README).
   struct Row { real_t depth, h_shear, n_stress; };
   const Row rows[] = {{5000.0, 13.83164, -87.82045}, {10000.0, 27.66329, -175.6409}};
   InitialStress is(Tpv26Spec());
   for (const Row &r : rows)
   {
      const real_t x[3] = {-5000.0, 0.0, -r.depth};
      real_t s[9];
      is.EvalEffective(x, s);
      // fault normal nhat = (0, -1, 0): along-strike shear sigma_XY, normal sigma_YY
      TEST_ASSERT(std::abs(s[1] / 1e6 - r.h_shear) <= 1e-6 * r.h_shear,
                  "U8 shear at " << r.depth / 1000 << " km = DRDG3D (7 digits)");
      TEST_ASSERT(std::abs(s[4] / 1e6 - r.n_stress) <= 1e-6 * std::abs(r.n_stress),
                  "U8 normal at " << r.depth / 1000 << " km = DRDG3D (7 digits)");
   }
}

void TestU8GeneralUp()
{
   std::cout << "\n=== U8: general up vector by rotation ===\n";
   // rotation R (row-major), orthonormal, det +1
   const real_t a = 0.3, b = -0.7, c = 1.1;
   const real_t Rx[9] = {1, 0, 0, 0, std::cos(a), -std::sin(a), 0, std::sin(a), std::cos(a)};
   const real_t Ry[9] = {std::cos(b), 0, std::sin(b), 0, 1, 0, -std::sin(b), 0, std::cos(b)};
   const real_t Rz[9] = {std::cos(c), -std::sin(c), 0, std::sin(c), std::cos(c), 0, 0, 0, 1};
   real_t T[9], R[9];
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         T[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { T[3 * i + j] += Ry[3 * i + k] * Rz[3 * k + j]; }
      }
   }
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         R[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { R[3 * i + j] += Rx[3 * i + k] * T[3 * k + j]; }
      }
   }
   InitialStressSpec s = Tpv26Spec();
   InitialStressSpec sr = s;
   Rotate9(R, s.B, sr.B);
   for (int i = 0; i < 3; i++) { sr.up[i] = R[3 * i + 2]; }   // R e_Z
   sr.surface_elevation = 0.0;
   InitialStress is(s), isr(sr);
   real_t worst = 0.0;
   for (real_t d : {0.0, 4000.0, 12000.0, 16500.0, 30000.0})
   {
      const real_t x[3] = {2500.0, -800.0, -d};
      real_t xr[3];
      for (int i = 0; i < 3; i++) { xr[i] = R[3 * i] * x[0] + R[3 * i + 1] * x[1] + R[3 * i + 2] * x[2]; }
      real_t sg[9], sgr[9], ref[9], diff[9];
      is.EvalEffective(x, sg);
      isr.EvalEffective(xr, sgr);
      Rotate9(R, sg, ref);
      for (int i = 0; i < 9; i++) { diff[i] = sgr[i] - ref[i]; }
      worst = std::max(worst, MaxAbs9(diff) / std::max(MaxAbs9(ref), real_t(1.0)));
   }
   TEST_ASSERT(worst <= 1e-12, "U8 sigma_0'(R x) = R sigma_0'(x) R^T for a rotated up vector");

   InitialStressSpec bad = Tpv26Spec();
   bad.B[8] = 0.1;   // vertical component
   TEST_ASSERT(InitialStress::Validate(bad).find("must be horizontal") != std::string::npos,
               "U8 a non-horizontal B is refused");
}

// ---------------------------------------------------------------- U9
void TestU9BelowYield()
{
   std::cout << "\n=== U9: initial state below yield (TPV27) ===\n";
   InitialStress is(Tpv26Spec());
   DruckerPragerParams p;
   p.cohesion = 1.36e6;
   p.bulk_friction = 0.1934;
   p.relaxation_time = 0.03;
   DruckerPrager dp(p);
   auto ratio = [&](real_t d)
   {
      const real_t x[3] = {0.0, 0.0, -d};
      real_t s6[6], sm, dev[6], tau;
      is.EvalEffectiveSym6(x, s6);
      DruckerPrager::Invariants(s6, sm, dev, tau);
      return tau / dp.YieldStress(sm);
   };
   real_t worst = 0.0, worst_d = 0.0;
   for (int k = 0; k <= 3000; k++)
   {
      const real_t d = 10.0 * k;
      const real_t r = ratio(d);
      if (r > worst) { worst = r; worst_d = d; }
   }
   std::cout << "  sqrt(J2)/Y: 10 km " << ratio(10000.0) << ", 15 km " << ratio(15000.0)
             << ", largest " << worst << " at " << worst_d << " m\n";
   TEST_ASSERT(worst < 1.0, "U9 below yield at every depth 0-30 km");
   TEST_ASSERT(std::abs(ratio(10000.0) - 0.930) < 5e-4, "U9 0.930 at 10 km (description: 93 %)");
   TEST_ASSERT(std::abs(ratio(15000.0) - 0.943) < 5e-4, "U9 0.943 at 15 km");
}

int main()
{
   std::cout << "Rheology initial-stress tests (plan rev 1, U8-U9)\n";
   TestU8Formulas();
   TestU8Drdg3d();
   TestU8GeneralUp();
   TestU9BelowYield();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
