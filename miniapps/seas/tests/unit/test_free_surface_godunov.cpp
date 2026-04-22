// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface.
//
// Validates that GodunovFlux::FreeSurfaceGodunov (Godunov characteristic
// projection) and GodunovFlux::FreeSurface (gamma-mirror) produce
// algebraically equivalent flux on a flat free surface.
//
// The equivalence is algebraically exact at ANY tilt angle: let
// delta := Q_ghost(gamma) - Q_god(projection).  Then A^- * delta = 0
// by construction of Q_god, so ApplySplitFlux gives identical output
// to the last FP bit.  Any tilt-dependent divergence indicates a bug
// in BuildFrame's Gram-Schmidt orthonormalisation or in the Q_god
// projection.  Hence the tolerance is uniform 10 ULP per component
// across all tilt angles (R-008 — replaces the looser multi-tier
// tolerance of the original v9.3.0 Phase 1 AC).
//
// Tests:
//   1. BIT-EXACT check on nor = (0, 0, 1) with axis-unit Q_self
//      (R-004 direct check).
//   2. 10 ULP per component on nor = (0, 0, 1) for 50 random Q_self.
//   3. 10 ULP per component at tilts 10, 30, 45, 75 degrees from
//      vertical.  Algebraic equivalence is frame-invariant (R-008).
//   4. Regression: sizeof(GodunovFlux) must not exceed the ctor-only
//      layout (no new cached compliance member per R-004).

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
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
// ULP helper — eps near 1.0 ~ 2.22e-16 in double precision.
//
// "10 ULP per component" interpretation (R-008): each component of the
// flux vector is compared at a tolerance scaled by the MAX magnitude
// ACROSS the vector.  Using a per-component scale is defeated by
// components that are identically zero in exact arithmetic (a rotation +
// split-flux chain leaves ~O(eps * ||F||) noise in such components,
// which is meaningless "infinite" ULP relative to 0).  The max-norm
// scale is the standard practice for "the two computations agree to
// 10 ULP in the flux vector" and matches the existing v9.2 regression
// tests (e.g. tests/unit/test_godunov_rotation_identity.cpp uses the
// same per-matrix-entry reasoning with an absolute tolerance budget).
// ---------------------------------------------------------------------------
static constexpr double kEps = 2.220446049250313e-16;
static constexpr int    kUlpBudget = 10;  // per R-008

// ---------------------------------------------------------------------------
// Compare F_gamma vs F_godunov in 10 ULP per component (max-norm scaled).
// ---------------------------------------------------------------------------
static bool FluxesAgreeUlp(const real_t F_gamma[NUM_STATE],
                           const real_t F_god[NUM_STATE],
                           double &max_ulp_out,
                           int &worst_comp_out)
{
   // Scale = max absolute value across both flux vectors.
   double scale = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
   {
      scale = std::max({scale, std::abs(F_gamma[c]), std::abs(F_god[c])});
   }
   if (scale == 0.0) { scale = 1.0; }  // both fluxes exactly zero

   max_ulp_out = 0.0;
   worst_comp_out = -1;
   bool ok = true;
   for (int c = 0; c < NUM_STATE; c++)
   {
      double err = std::abs(F_gamma[c] - F_god[c]);
      double ulp = err / (scale * kEps);
      if (ulp > max_ulp_out) { max_ulp_out = ulp; worst_comp_out = c; }
      if (err > kUlpBudget * kEps * scale) { ok = false; }
   }
   return ok;
}

// ---------------------------------------------------------------------------
// Rodrigues rotation of a vector about an axis by an angle.
// ---------------------------------------------------------------------------
static void RotateVec(const real_t axis[3], double theta,
                      const real_t v[3], real_t w[3])
{
   const double c = std::cos(theta);
   const double s = std::sin(theta);
   const double one_minus_c = 1.0 - c;
   const double ax = axis[0], ay = axis[1], az = axis[2];
   double cross[3] = {
      ay * v[2] - az * v[1],
      az * v[0] - ax * v[2],
      ax * v[1] - ay * v[0]
   };
   double dot = ax * v[0] + ay * v[1] + az * v[2];
   w[0] = v[0] * c + cross[0] * s + ax * dot * one_minus_c;
   w[1] = v[1] * c + cross[1] * s + ay * dot * one_minus_c;
   w[2] = v[2] * c + cross[2] * s + az * dot * one_minus_c;
}

// ---------------------------------------------------------------------------
// Run the equivalence check at a given normal for a population of states.
// ---------------------------------------------------------------------------
static bool CheckAtNormal(const GodunovFlux &flux, const real_t nor[3],
                          const std::vector<std::array<real_t, NUM_STATE>> &Qs,
                          const std::string &label)
{
   double worst_ulp = 0.0;
   int worst_q = -1, worst_c = -1;
   bool all_ok = true;
   for (size_t k = 0; k < Qs.size(); k++)
   {
      const real_t *Q = Qs[k].data();
      real_t F_gamma[NUM_STATE];
      real_t F_god  [NUM_STATE];
      flux.FreeSurface       (nor, Q, F_gamma);
      flux.FreeSurfaceGodunov(nor, Q, F_god);
      double max_ulp; int bad_c;
      bool ok = FluxesAgreeUlp(F_gamma, F_god, max_ulp, bad_c);
      if (max_ulp > worst_ulp) { worst_ulp = max_ulp; worst_q = (int)k; worst_c = bad_c; }
      if (!ok) { all_ok = false; }
   }
   std::cout << "  [" << label << "] worst ULP = " << std::scientific
             << std::setprecision(3) << worst_ulp
             << " at Q-index " << worst_q << ", component " << worst_c
             << "  (budget = " << kUlpBudget << ")\n";
   TEST_ASSERT(all_ok,
               std::string("10 ULP per component [") + label + "]");
   return all_ok;
}

// ---------------------------------------------------------------------------
// Test 1: R-004 direct bit-check at nor = (0, 0, 1), Q_self with SXX=1.
// ---------------------------------------------------------------------------
void TestR004DirectCheck(const GodunovFlux &flux)
{
   std::cout << "\n-- Test 1: R-004 direct check --\n";
   real_t nor[3] = {0.0, 0.0, 1.0};
   real_t Q[NUM_STATE] = {0.0};
   Q[SXX] = 1.0;

   real_t F_gamma[NUM_STATE];
   real_t F_god  [NUM_STATE];
   flux.FreeSurface       (nor, Q, F_gamma);
   flux.FreeSurfaceGodunov(nor, Q, F_god);

   double worst_err = 0.0;
   int    worst_c   = -1;
   double scale     = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
   {
      double err = std::abs(F_gamma[c] - F_god[c]);
      if (err > worst_err) { worst_err = err; worst_c = c; }
      scale = std::max({scale, std::abs(F_gamma[c]), std::abs(F_god[c])});
   }
   if (scale == 0.0) { scale = 1.0; }
   double ulp = worst_err / (scale * kEps);
   std::cout << "  worst |F_gamma - F_god| = " << std::scientific
             << worst_err << " (" << ulp << " ULP of max-norm) at component "
             << worst_c << "\n";

   // R-004 AC (review-incorporation plan line 232-238): within 1 ULP per
   // component, with a max-norm scale across all 9 flux components.
   // Fix for REVIEW.md R-I04-001 (which flagged the previous single-
   // component-scale, 16-ULP-tolerance assertion as contradicting its
   // own "1 ULP" message).
   TEST_ASSERT(worst_err <= 1 * kEps * scale,
               "R-004: FreeSurfaceGodunov == FreeSurface on Q[SXX]=1 "
               "nor=(0,0,1) within 1 ULP of the max-norm flux scale");
}

// ---------------------------------------------------------------------------
// Test 2: 50 random Q_self at nor = (0, 0, 1).
// ---------------------------------------------------------------------------
void TestAxisAlignedRandom(const GodunovFlux &flux)
{
   std::cout << "\n-- Test 2: 50 random Q_self at nor = (0, 0, 1) --\n";

   std::mt19937_64 rng(42);
   std::uniform_real_distribution<double> stress_d(-1.0e8, 1.0e8);
   std::uniform_real_distribution<double> vel_d   (-1.0e2, 1.0e2);

   std::vector<std::array<real_t, NUM_STATE>> Qs;
   Qs.reserve(50);
   for (int k = 0; k < 50; k++)
   {
      std::array<real_t, NUM_STATE> Q{};
      for (int c = 0; c < VX; c++)  { Q[c] = stress_d(rng); }
      for (int c = VX; c < NUM_STATE; c++) { Q[c] = vel_d(rng); }
      Qs.push_back(Q);
   }

   real_t nor[3] = {0.0, 0.0, 1.0};
   CheckAtNormal(flux, nor, Qs, "nor=(0,0,1)");
}

// ---------------------------------------------------------------------------
// Test 3: tilts 10, 30, 45, 75 degrees from vertical (R-008 uniform 10 ULP).
// ---------------------------------------------------------------------------
void TestTiltedNormals(const GodunovFlux &flux)
{
   std::cout << "\n-- Test 3: tilted normals (10, 30, 45, 75 deg) --\n";

   std::mt19937_64 rng(123);
   std::uniform_real_distribution<double> stress_d(-1.0e8, 1.0e8);
   std::uniform_real_distribution<double> vel_d   (-1.0e2, 1.0e2);

   std::vector<std::array<real_t, NUM_STATE>> Qs;
   Qs.reserve(50);
   for (int k = 0; k < 50; k++)
   {
      std::array<real_t, NUM_STATE> Q{};
      for (int c = 0; c < VX; c++)  { Q[c] = stress_d(rng); }
      for (int c = VX; c < NUM_STATE; c++) { Q[c] = vel_d(rng); }
      Qs.push_back(Q);
   }

   const double tilt_angles_deg[] = {10.0, 30.0, 45.0, 75.0};
   real_t z_axis[3] = {0.0, 0.0, 1.0};
   real_t x_axis[3] = {1.0, 0.0, 0.0};

   for (double tilt_deg : tilt_angles_deg)
   {
      real_t nor[3];
      double theta = tilt_deg * M_PI / 180.0;
      RotateVec(x_axis, theta, z_axis, nor);
      // Normalize to kill rounding error.
      double nlen = std::sqrt(nor[0]*nor[0] + nor[1]*nor[1] + nor[2]*nor[2]);
      nor[0] /= nlen; nor[1] /= nlen; nor[2] /= nlen;

      std::string label = "tilt=" + std::to_string((int)tilt_deg) + "deg";
      CheckAtNormal(flux, nor, Qs, label);
   }
}

// ---------------------------------------------------------------------------
// Test 4: R-004 regression — the FreeSurfaceGodunov compliance block
// must remain INLINED (invZp/invZs), not cached as a new DenseMatrix.
//
// Original design: the class layout is 7 real_t (lambda, mu, rho, cp, cs,
// Zp, Zs) + 3 DenseMatrix (Ax, Ax_plus, Ax_minus).  R-I04-004 required
// that `FreeSurfaceGodunov` did NOT add a cached 9×9 compliance block.
//
// Updated (ADER I-05 Phase 2): the class now also carries 3 reference
// star matrices (A_x, A_y, A_z — `ref_star_`) for the CK recursion.
// That's +3 DenseMatrix entries, legitimately added for a DIFFERENT
// purpose.  The R-004 guard still fires if a FOURTH extra DenseMatrix
// appears (e.g. a re-added compliance-block cache).
// ---------------------------------------------------------------------------
void TestR004NoNewMember()
{
   std::cout << "\n-- Test 4: R-004 regression — no new cached member --\n";

   constexpr std::size_t field_floor =
      7 * sizeof(real_t) + 6 * sizeof(DenseMatrix);  // +3 for ADER ref_star_
   constexpr std::size_t field_ceiling_hint =
      field_floor + sizeof(DenseMatrix);  // would indicate a re-added cache

   std::cout << "  sizeof(GodunovFlux)  = " << sizeof(GodunovFlux) << " bytes\n";
   std::cout << "  field floor (7 real + 6 DenseMatrix) ~= "
             << field_floor << " bytes\n";
   std::cout << "  ceiling hint  (+DenseMatrix)          ~= "
             << field_ceiling_hint << " bytes\n";

   // Don't hard-assert an exact byte size — alignment/padding vary.
   // Assert the class is smaller than the ceiling hint to catch the
   // specific regression R-004 guards against (re-added compliance cache).
   TEST_ASSERT(sizeof(GodunovFlux) < field_ceiling_hint,
               "R-004: GodunovFlux did not acquire a new "
               "DenseMatrix-sized cached member beyond ref_star_");
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 1 (I-04): "
             << "FreeSurfaceGodunov vs FreeSurface equivalence ===\n";

   // TPV102 homogeneous material (matches config/tpv102_params.hpp).
   const real_t lambda = 32.04e9;
   const real_t mu     = 32.04e9;
   const real_t rho    = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   std::cout << "  Material: lambda=" << lambda
             << ", mu=" << mu << ", rho=" << rho << "\n";
   std::cout << "  Zp = " << flux.GetZp()
             << ",  Zs = " << flux.GetZs() << "\n";

   TestR004DirectCheck (flux);
   TestAxisAlignedRandom(flux);
   TestTiltedNormals   (flux);
   TestR004NoNewMember ();

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
