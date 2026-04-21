// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.1.0 §4.2 — post-R-003 FaultBasis orthonormality contract.
//
// Verifies that after the R-003 fix (explicit `dip = normalize(dip)` in
// `FaultBasis::ComputeOrientedFrame`, plan v9.1.0 §3.4 H-V91-B1),
// the dip and strike tangents satisfy an orthonormal basis with the
// face normal to ~10 ULP across 1000 perturbed-normal samples.
//
// Pre-R-003 the strike was explicitly normalized but the dip was
// computed as `strike x n_raw` without an explicit normalize; FP
// rounding in the strike normalize leaked into the cross product so
// |dip| drifted up to ~5 ULP.  The test documents the pre-fix drift
// as a regression gate — deleting the R-003 patch should make this
// test fail on the |dip| assertions.

#include "mfem.hpp"
#include "../../fault/fault_basis.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Test framework
// ---------------------------------------------------------------------------

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) <= t_) { num_passed++; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ \
      << ", |diff|=" << std::abs(v_-e_) << ", tol=" << t_ << ")\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// Lightweight deterministic PRNG (xorshift64*) — no <random> dependency
// ---------------------------------------------------------------------------

static inline uint64_t xorshift64(uint64_t &state)
{
   state ^= state << 13;
   state ^= state >> 7;
   state ^= state << 17;
   return state;
}

static inline double unit_rand(uint64_t &state)
{
   // Uniform in [-1, 1] from 53 bits of entropy.
   const uint64_t u = xorshift64(state) >> 11;
   const double x = static_cast<double>(u) / (1ULL << 53);
   return 2.0 * x - 1.0;
}

// ---------------------------------------------------------------------------
// 3-vector helpers
// ---------------------------------------------------------------------------

static inline double norm3(const real_t v[3])
{
   return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double dot3(const real_t a[3], const real_t b[3])
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ---------------------------------------------------------------------------
// One trial: perturb the ref_normal by a bounded random 3-vector,
// normalize, feed to ComputeOrientedFrame, assert orthonormality.
// ---------------------------------------------------------------------------

static void RunTrial(uint64_t &prng_state, int trial_id,
                     int &fail_acc_len, int &fail_acc_ortho)
{
   // TPV102 canonical choice: ref_normal = (0, -1, 0), up = (0, 0, 1).
   const int dim = 3;
   Vector ref_normal(dim);
   ref_normal(0) = 0.0; ref_normal(1) = -1.0; ref_normal(2) = 0.0;
   Vector up(dim);
   up(0) = 0.0; up(1) = 0.0; up(2) = 1.0;

   // Build n_raw = ref_normal + delta, ||delta|| <= 0.1, then normalize.
   // Draw delta components iid uniform in [-0.1/sqrt(3), +0.1/sqrt(3)]
   // so ||delta||_2 <= 0.1 with high probability.
   const double delta_bound = 0.1 / std::sqrt(3.0);
   Vector n_raw(dim);
   double dnorm2 = 0.0;
   for (int d = 0; d < dim; d++)
   {
      const double delta = delta_bound * unit_rand(prng_state);
      n_raw(d) = ref_normal(d) + delta;
      dnorm2 += (n_raw(d) - ref_normal(d)) *
                (n_raw(d) - ref_normal(d));
   }
   // Deliberately leave n_raw unnormalized — ComputeOrientedFrame
   // normalizes internally (line `n_raw /= nl`).  But scale it to a
   // physically-plausible raw-length magnitude to exercise the nl
   // branch (not strictly required; matches CalcOrtho output scale).
   n_raw *= 2.7;  // arbitrary >0 scaling

   real_t normal[3], tangent1[3], tangent2[3];
   bool sign_flipped = false;
   real_t nl = 0.0;

   FaultBasis::ComputeOrientedFrame(n_raw, dim, ref_normal, up,
                                    normal, tangent1, tangent2,
                                    sign_flipped, nl);

   const double eps = std::numeric_limits<double>::epsilon();
   // R-006 (v9.1.0 rev 3): tol_len tightened from 10*eps to 2*eps.
   // Post-R-003 the dip is explicitly normalized (`dv *= 1.0/|dv|`), so
   // |dv| == 1 to ~1 ULP; 2*eps catches the pre-R-003 drift (1-4 ULP per
   // numerical analysis) that 10*eps would silently let through.  A
   // revert of the R-003 normalize would make this test fail on the
   // length checks.  tol_dot stays at 10*eps: orthogonality is a
   // different invariant with a broader ULP band.
   const double tol_len = 2.0 * eps;
   const double tol_dot = 10.0 * eps;

   const double n_len = norm3(normal);
   const double t1_len = norm3(tangent1);   // dip
   const double t2_len = norm3(tangent2);   // strike
   const double dot_n_t1 = dot3(normal, tangent1);
   const double dot_n_t2 = dot3(normal, tangent2);
   const double dot_t1_t2 = dot3(tangent1, tangent2);

   // The two lengths and three pairwise orthogonalities are the
   // orthonormality contract.  We aggregate pass/fail counts rather than
   // emitting 3000 TEST_NEAR lines — the harness PASS count remains the
   // aggregate (n_trials lengths + n_trials orthos per channel).

   if (std::abs(n_len - 1.0) > tol_len) { fail_acc_len++; }
   if (std::abs(t1_len - 1.0) > tol_len) { fail_acc_len++; }
   if (std::abs(t2_len - 1.0) > tol_len) { fail_acc_len++; }

   if (std::abs(dot_n_t1) > tol_dot)  { fail_acc_ortho++; }
   if (std::abs(dot_n_t2) > tol_dot)  { fail_acc_ortho++; }
   if (std::abs(dot_t1_t2) > tol_dot) { fail_acc_ortho++; }

   // On the first trial, emit concrete TEST_NEAR lines so the log
   // records actual FP residuals.
   if (trial_id == 0)
   {
      TEST_NEAR(n_len,  1.0, tol_len, "trial 0 |normal|  == 1");
      TEST_NEAR(t1_len, 1.0, tol_len, "trial 0 |tangent1 = dip|    == 1 (R-003)");
      TEST_NEAR(t2_len, 1.0, tol_len, "trial 0 |tangent2 = strike| == 1");
      TEST_NEAR(dot_n_t1,  0.0, tol_dot, "trial 0 dot(normal, dip)    == 0");
      TEST_NEAR(dot_n_t2,  0.0, tol_dot, "trial 0 dot(normal, strike) == 0");
      TEST_NEAR(dot_t1_t2, 0.0, tol_dot, "trial 0 dot(dip, strike)    == 0");
   }
}

// ---------------------------------------------------------------------------
// Additional canonical check: exact TPV102 axis-aligned frame.
// ref_normal = (0, -1, 0), up = (0, 0, 1) => strike = (1, 0, 0),
//                                            dip    = (0, 0, -1).
// ---------------------------------------------------------------------------

static void TestCanonicalTPV102Frame()
{
   std::cout << "\n=== Canonical TPV102 axis-aligned frame ===\n";
   const int dim = 3;
   Vector ref_normal(dim);
   ref_normal(0) = 0.0; ref_normal(1) = -1.0; ref_normal(2) = 0.0;
   Vector up(dim);
   up(0) = 0.0; up(1) = 0.0; up(2) = 1.0;
   Vector n_raw(ref_normal);

   real_t normal[3], tangent1[3], tangent2[3];
   bool sign_flipped = false;
   real_t nl = 0.0;
   FaultBasis::ComputeOrientedFrame(n_raw, dim, ref_normal, up,
                                    normal, tangent1, tangent2,
                                    sign_flipped, nl);

   // Exact equalities in IEEE-754 for an axis-aligned construction.
   const double tol = 1.0e-15;
   TEST_NEAR(normal[0],    0.0, tol, "canonical normal[0]");
   TEST_NEAR(normal[1],   -1.0, tol, "canonical normal[1]");
   TEST_NEAR(normal[2],    0.0, tol, "canonical normal[2]");
   TEST_NEAR(tangent1[0],  0.0, tol, "canonical tangent1(dip)[0]");
   TEST_NEAR(tangent1[1],  0.0, tol, "canonical tangent1(dip)[1]");
   TEST_NEAR(tangent1[2], -1.0, tol, "canonical tangent1(dip)[2]");
   TEST_NEAR(tangent2[0],  1.0, tol, "canonical tangent2(strike)[0]");
   TEST_NEAR(tangent2[1],  0.0, tol, "canonical tangent2(strike)[1]");
   TEST_NEAR(tangent2[2],  0.0, tol, "canonical tangent2(strike)[2]");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.1.0 §4.2 — FaultBasis orthonormality (R-003)\n";
   std::cout << "========================================\n";

   TestCanonicalTPV102Frame();

   std::cout << "\n=== 1000 perturbed-normal trials ===\n";
   uint64_t state = 0x9E3779B97F4A7C15ULL;  // arbitrary non-zero seed
   int fail_len = 0, fail_ortho = 0;
   const int n_trials = 1000;
   for (int t = 0; t < n_trials; t++)
   {
      RunTrial(state, t, fail_len, fail_ortho);
   }

   // Aggregate PASS if zero failures across 3*n_trials length checks
   // and 3*n_trials orthogonality checks.
   num_tests++;
   if (fail_len == 0)
   {
      num_passed++;
      std::cout << "  PASSED: 3*" << n_trials
                << " |length - 1| <= 2*eps checks (R-003, R-006)\n";
   }
   else
   {
      num_failed++;
      std::cout << "  FAILED: " << fail_len << " of " << (3*n_trials)
                << " length checks exceeded 2*eps (R-003, R-006)\n";
   }
   num_tests++;
   if (fail_ortho == 0)
   {
      num_passed++;
      std::cout << "  PASSED: 3*" << n_trials
                << " |dot_pair| <= 10*eps checks\n";
   }
   else
   {
      num_failed++;
      std::cout << "  FAILED: " << fail_ortho << " of " << (3*n_trials)
                << " orthogonality checks exceeded 10*eps\n";
   }

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}
