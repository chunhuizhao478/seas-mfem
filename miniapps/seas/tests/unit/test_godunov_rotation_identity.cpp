// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.1 (regression gate): rotation identity on the
// Voigt-9 stress/velocity state vector.
//
// Asserts `max|T · Tinv − I| ≤ 16 ULP` on the full 9x9 rotation block
// produced by GodunovFlux::BuildRotation/BuildRotationInverse for:
//   1. The BP5 canonical frame  ((0,-1,0), (0,0,-1), (1,0,0))
//        — expected bit-exact (signed permutation Q; integer arithmetic)
//   2. Three non-axis-aligned frames obtained by small Euler-angle
//      perturbations from BP5 (15 deg, 30 deg, 45 deg total rotation).
//      — expected <= 16 ULP per REVIEW.md R-004 (9-term FMA sum).
//
// Rationale (plan §4.1 + §A): The v9.2.0 rev-1 draft listed
// H-V92-R (rotation stress-block asymmetry) as a rank-1 candidate.
// REVIEW.md R-001 proved mathematically that `(T·Tinv)[ij,i'j'] = δ`
// exactly for any orthonormal Q.  This test encodes that proof as a
// regression gate.  If the rotation block is ever regressed (wrong
// Voigt symmetrization sign, wrong factor of 2, etc.) this test
// fails first.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// ULP helper — epsilon near 1.0 is ~2.22e-16 in double precision.
// ---------------------------------------------------------------------------
static constexpr double kEps = 2.220446049250313e-16;
static constexpr int    kUlpBudget = 16;  // REVIEW R-004
static constexpr double kTolUlp    = kUlpBudget * kEps;

// ---------------------------------------------------------------------------
// Measure max|T·Tinv − I| across all 81 entries.
// ---------------------------------------------------------------------------
static double MaxDeviationFromIdentity(const real_t n[3],
                                        const real_t t1[3],
                                        const real_t t2[3])
{
   DenseMatrix T, Tinv;
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   GodunovFlux::BuildRotation        (n, t1, t2, T);

   DenseMatrix M(NUM_STATE, NUM_STATE);
   Mult(T, Tinv, M);

   double max_dev = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      for (int j = 0; j < NUM_STATE; j++)
      {
         double expected = (i == j) ? 1.0 : 0.0;
         double dev = std::abs(M(i, j) - expected);
         if (dev > max_dev) { max_dev = dev; }
      }
   }
   return max_dev;
}

// ---------------------------------------------------------------------------
// Rotate the basis (n, t1, t2) about a given axis by an angle.
// Rodrigues' formula, applied to each basis vector independently.
// ---------------------------------------------------------------------------
static void RotateBasis(const real_t axis[3], double theta,
                        const real_t n_in[3], const real_t t1_in[3],
                        const real_t t2_in[3],
                        real_t n_out[3], real_t t1_out[3], real_t t2_out[3])
{
   const double c = std::cos(theta);
   const double s = std::sin(theta);
   const double one_minus_c = 1.0 - c;
   const double ax = axis[0], ay = axis[1], az = axis[2];

   auto rot = [&](const real_t v[3], real_t w[3])
   {
      // Rodrigues: w = v*c + (axis x v)*s + axis*(axis.v)*(1-c)
      double cross[3] = {
         ay * v[2] - az * v[1],
         az * v[0] - ax * v[2],
         ax * v[1] - ay * v[0]
      };
      double dot = ax * v[0] + ay * v[1] + az * v[2];
      w[0] = v[0] * c + cross[0] * s + ax * dot * one_minus_c;
      w[1] = v[1] * c + cross[1] * s + ay * dot * one_minus_c;
      w[2] = v[2] * c + cross[2] * s + az * dot * one_minus_c;
   };

   rot(n_in,  n_out);
   rot(t1_in, t1_out);
   rot(t2_in, t2_out);
}

// ---------------------------------------------------------------------------
// Check one frame: print max deviation and assert ≤ 16 ULP.
// axis_aligned=true uses stricter 0-ULP assertion (expected bit-exact).
// ---------------------------------------------------------------------------
static void CheckFrame(const real_t n[3], const real_t t1[3],
                       const real_t t2[3], const std::string &label,
                       bool axis_aligned)
{
   double dev = MaxDeviationFromIdentity(n, t1, t2);
   std::cout << "  frame [" << label << "]\n";
   std::cout << "    n  = (" << n[0]  << ", " << n[1]  << ", " << n[2]  << ")\n";
   std::cout << "    t1 = (" << t1[0] << ", " << t1[1] << ", " << t1[2] << ")\n";
   std::cout << "    t2 = (" << t2[0] << ", " << t2[1] << ", " << t2[2] << ")\n";
   std::cout << "    max|T·Tinv - I| = " << std::scientific
             << std::setprecision(3) << dev
             << " (" << dev / kEps << " ulp)\n";

   if (axis_aligned)
   {
      TEST_ASSERT(dev == 0.0,
                  std::string("BIT-EXACT identity on axis-aligned frame [")
                  + label + "]");
   }
   else
   {
      TEST_ASSERT(dev <= kTolUlp,
                  std::string("max|T·Tinv - I| <= 16 ULP on frame [")
                  + label + "]");
   }
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 §4.1 Regression: Rotation Identity ===\n";

   // -------- Frame 1: BP5 canonical (axis-aligned signed permutation) ------
   // ref_normal = (0,-1,0), up = (0,0,1):
   //   n       = (0, -1, 0)
   //   strike  = up x n = (0,0,1) x (0,-1,0) = (1, 0, 0) = tangent2
   //   dip     = strike x n = (1,0,0) x (0,-1,0) = (0, 0, -1) = tangent1
   {
      real_t n[3]  = {0.0, -1.0, 0.0};
      real_t t1[3] = {0.0,  0.0, -1.0};
      real_t t2[3] = {1.0,  0.0,  0.0};
      std::cout << "\n-- BP5 canonical frame (axis-aligned, expect bit-exact) --\n";
      CheckFrame(n, t1, t2, "BP5_canonical", /*axis_aligned=*/true);
   }

   // -------- Frames 2-4: non-axis-aligned perturbations of BP5 -------------
   const real_t n0[3]  = {0.0, -1.0, 0.0};
   const real_t t10[3] = {0.0,  0.0, -1.0};
   const real_t t20[3] = {1.0,  0.0,  0.0};

   struct PerturbSpec { real_t axis[3]; double theta_deg; const char *label; };
   PerturbSpec specs[] = {
      {{1.0, 0.0, 0.0}, 15.0, "BP5+15deg_about_x"},
      {{0.0, 1.0, 0.0}, 30.0, "BP5+30deg_about_y"},
      {{std::sqrt(1.0/3.0), std::sqrt(1.0/3.0), std::sqrt(1.0/3.0)},
        45.0, "BP5+45deg_about_111"},
   };
   std::cout << "\n-- Non-axis-aligned perturbations (expect <= 16 ULP) --\n";
   for (auto &sp : specs)
   {
      real_t n[3], t1[3], t2[3];
      double theta = sp.theta_deg * M_PI / 180.0;
      RotateBasis(sp.axis, theta, n0, t10, t20, n, t1, t2);
      CheckFrame(n, t1, t2, sp.label, /*axis_aligned=*/false);
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
