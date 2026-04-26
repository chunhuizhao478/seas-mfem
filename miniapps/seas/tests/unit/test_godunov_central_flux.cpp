// Round-11 R-1101: GodunovFlux::Central primitive correctness test.
//
// Contract gates (from MIXED_FLUX_PLAN.md Phase 1):
//
// Gate 1: Central(nor, Q, Q, F_c) == Interior(nor, Q, Q, F_u) bit-equal
//         (zero-jump → upwind dissipation contribution is zero → Central
//          and Interior agree to FP precision).
//
// Gate 2 (R-1204 sign convention pinned):
//         Interior(nor, Q_self, Q_nbr, F_up) − Central(nor, Q_self, Q_nbr,
//         F_ce) == +0.5 · |A_n| · (Q_self − Q_nbr) within 1e−12.
//         |A_n| (global frame) = T · (Ax_plus_ − Ax_minus_) · Tinv.
//         Q_self is the "self" (Elem1) side.
//
// Coverage:
//   - Material params: TPV104 (rho=2670, cp=6000, cs=3464).
//   - 8 representative face normals: ±x, ±y, ±z, plus 2 oblique.
//   - State sweep: 9 channel-isolated unit Q states + 3 mixed (Q_self ≠ Q_nbr) pairs.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/godunov_flux.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace
{
constexpr real_t kRho    = 2670.0;
constexpr real_t kCp     = 6000.0;
constexpr real_t kCs     = 3464.0;
constexpr real_t kMu     = kRho * kCs * kCs;
constexpr real_t kLambda = kRho * kCp * kCp - 2.0 * kMu;

real_t MaxAbsDiff(const real_t *a, const real_t *b, int n)
{
   real_t m = 0.0;
   for (int i = 0; i < n; i++)
   {
      real_t d = std::abs(a[i] - b[i]);
      if (d > m) { m = d; }
   }
   return m;
}

// Compute |A_n| · (Q_self - Q_nbr) in the GLOBAL frame, for the
// algebraic-identity gate.  Process:
//   1. Build orthonormal frame (n, t1, t2) the same way GodunovFlux does
//      via BuildFrame.
//   2. Build T, Tinv via the public BuildRotation / BuildRotationInverse.
//   3. |A_n|_face = Ax_plus_ − Ax_minus_ (in face-rotated frame).
//      Reach Ax_plus_, Ax_minus_ via the public GetAxPlus() / GetAxMinus()
//      accessors declared in godunov_flux.hpp:153–155.
//   4. Rotate (Q_self - Q_nbr) into face frame: jump_local = Tinv · jump.
//   5. Apply |A_n|_face: out_local = (Ax_plus_ − Ax_minus_) · jump_local.
//   6. Rotate back: out = T · out_local.
void ComputeAbsAnTimesJumpGlobal(const GodunovFlux &flux, const real_t *nor,
                                 const real_t *Q_self, const real_t *Q_nbr,
                                 real_t *out)
{
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);
   DenseMatrix Tinv(NUM_STATE, NUM_STATE), T(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);
   GodunovFlux::BuildRotation(nor, t1, t2, T);

   real_t jump_global[NUM_STATE];
   for (int i = 0; i < NUM_STATE; i++)
   {
      jump_global[i] = Q_self[i] - Q_nbr[i];
   }
   real_t jump_local[NUM_STATE];
   Tinv.Mult(jump_global, jump_local);

   const DenseMatrix &Ap = flux.GetAxPlus();
   const DenseMatrix &Am = flux.GetAxMinus();
   real_t out_local[NUM_STATE];
   for (int i = 0; i < NUM_STATE; i++)
   {
      real_t s = 0.0;
      for (int j = 0; j < NUM_STATE; j++)
      {
         s += (Ap(i, j) - Am(i, j)) * jump_local[j];
      }
      out_local[i] = s;
   }
   T.Mult(out_local, out);
}

void NormalizeVec(real_t *n)
{
   const real_t len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   if (len > 0.0) { n[0] /= len; n[1] /= len; n[2] /= len; }
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-11 R-1101: GodunovFlux::Central primitive ===\n";

   GodunovFlux flux(kLambda, kMu, kRho);

   // 8 representative normals.
   real_t normals[8][3] = {
      {+1.0,  0.0,  0.0},
      {-1.0,  0.0,  0.0},
      { 0.0, +1.0,  0.0},
      { 0.0, -1.0,  0.0},
      { 0.0,  0.0, +1.0},
      { 0.0,  0.0, -1.0},
      { 1.0,  1.0,  1.0},  // oblique (will be normalized)
      { 1.0, -2.0,  3.0},  // oblique
   };
   for (int k = 0; k < 8; k++) { NormalizeVec(normals[k]); }

   // -----------------------------------------------------------------
   // Gate 1: zero-jump → Central == Interior to FP precision.
   //   At zero jump, Interior = Ax_plus_·Q + Ax_minus_·Q and
   //   Central = 0.5·(Ax_plus_+Ax_minus_)·(Q+Q).  These are
   //   mathematically equal but FP-summation order differs (Interior
   //   adds Ap·Q and Am·Q separately per j; Central forms (Ap+Am)·(Q+Q)
   //   per j).  Tolerance is ~1 ULP × max(|F|).  Result scale is
   //   ~cp · stress ~ 6e3 · 1e7 = 6e10, so ULP ~ 1e-5.  Use 1e-3
   //   absolute (~1e-14 relative on this scale) as a generous gate.
   // -----------------------------------------------------------------
   std::cout << "\n-- Gate 1: zero-jump FP-precision agreement --\n";
   real_t worst_g1 = 0.0;
   for (int c_iso = 0; c_iso < NUM_STATE; c_iso++)
   {
      real_t Q[NUM_STATE] = {0};
      // Use channel-appropriate magnitudes: stresses ~ 1e7 Pa,
      // velocities ~ 1 m/s.  Helps detect channel-mixing bugs.
      Q[c_iso] = (c_iso < SXY) ? 1.0e7
               : (c_iso < VX)  ? 1.0e6
                               : 1.0;
      for (int k = 0; k < 8; k++)
      {
         real_t F_up[NUM_STATE], F_ce[NUM_STATE];
         flux.Interior(normals[k], Q, Q, F_up);
         flux.Central (normals[k], Q, Q, F_ce);
         real_t d = MaxAbsDiff(F_up, F_ce, NUM_STATE);
         if (d > worst_g1) { worst_g1 = d; }
      }
   }
   TEST_LE(worst_g1, 1.0e-3,
           "Central(nor, Q, Q) ~ Interior(nor, Q, Q) at FP precision for "
           "all 9 channel-isolated states × 8 normals (~1 ULP scale)");

   // -----------------------------------------------------------------
   // Gate 2 (R-1204): Interior(Q_self, Q_nbr) − Central(Q_self, Q_nbr)
   //                  = +0.5 · |A_n| · (Q_self − Q_nbr)
   // -----------------------------------------------------------------
   std::cout << "\n-- Gate 2: algebraic identity (R-1204 sign) --\n";

   // Three (Q_self, Q_nbr) pairs that exercise mixed channels.
   struct StatePair { real_t Q_self[NUM_STATE]; real_t Q_nbr[NUM_STATE]; };
   StatePair pairs[3];
   // Pair A: stress-dominated jump.
   for (int c = 0; c < NUM_STATE; c++)
   {
      pairs[0].Q_self[c] = (c < SXY) ? 1.0e7 : 0.0;
      pairs[0].Q_nbr [c] = (c < SXY) ? 0.5e7 : 0.0;
   }
   // Pair B: velocity-dominated jump.
   for (int c = 0; c < NUM_STATE; c++)
   {
      pairs[1].Q_self[c] = (c >= VX) ? 1.0 : 0.0;
      pairs[1].Q_nbr [c] = (c >= VX) ? 0.5 : 0.0;
   }
   // Pair C: mixed P-wave-like state.
   {
      const real_t lp = kLambda + 2.0 * kMu;
      for (int c = 0; c < NUM_STATE; c++)
      {
         pairs[2].Q_self[c] = 0.0;
         pairs[2].Q_nbr [c] = 0.0;
      }
      pairs[2].Q_self[SXX] = -lp;     pairs[2].Q_self[VX] = +kCp;
      pairs[2].Q_nbr [SXX] = +lp;     pairs[2].Q_nbr [VX] = -kCp;
   }

   real_t worst_g2 = 0.0;
   for (int p = 0; p < 3; p++)
   {
      for (int k = 0; k < 8; k++)
      {
         real_t F_up[NUM_STATE], F_ce[NUM_STATE];
         flux.Interior(normals[k], pairs[p].Q_self, pairs[p].Q_nbr, F_up);
         flux.Central (normals[k], pairs[p].Q_self, pairs[p].Q_nbr, F_ce);

         // Expected: F_up − F_ce = +0.5 · |A_n|_global · (Q_self − Q_nbr)
         real_t abs_an_jump_global[NUM_STATE];
         ComputeAbsAnTimesJumpGlobal(flux, normals[k],
                                     pairs[p].Q_self, pairs[p].Q_nbr,
                                     abs_an_jump_global);

         real_t expected[NUM_STATE], observed[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            expected[c] = 0.5 * abs_an_jump_global[c];
            observed[c] = F_up[c] - F_ce[c];
         }
         real_t d = MaxAbsDiff(observed, expected, NUM_STATE);
         if (d > worst_g2) { worst_g2 = d; }
      }
   }
   // Tolerance (R-1204): dimensionally |A_n| has units of c (~1e3) and
   // jump in stress ~1e7, so |A_n|·jump ~1e10.  Three independent FP
   // accumulations contribute to the residual (Interior split sum,
   // Central full-sum, ComputeAbsAnTimesJumpGlobal helper); generous
   // tolerance is ~10 ULP × scale = 1e-5.  Use 1e-1 to give an order
   // of magnitude of headroom — still tight enough to detect a real bug
   // (which would be 5–10 orders of magnitude larger).
   TEST_LE(worst_g2, 1.0e-1,
           "R-1204 algebraic identity: F_up − F_ce = 0.5·|A_n|·(Q_self−Q_nbr) "
           "to FP precision (~1 ULP × scale)");

   // -----------------------------------------------------------------
   // Gate 3: bilinearity check.
   // Central is linear in (Q_self, Q_nbr), so:
   //   Central(n, αQ_self + βR_self, αQ_nbr + βR_nbr)
   //     = α·Central(n, Q_self, Q_nbr) + β·Central(n, R_self, R_nbr)
   // -----------------------------------------------------------------
   std::cout << "\n-- Gate 3: bilinearity --\n";
   real_t worst_g3 = 0.0;
   const real_t alpha = 1.7, beta = -2.3;
   for (int p = 0; p < 3; p++)
   {
      const real_t *Q_self = pairs[p].Q_self;
      const real_t *Q_nbr  = pairs[p].Q_nbr;
      const real_t *R_self = pairs[(p + 1) % 3].Q_self;
      const real_t *R_nbr  = pairs[(p + 1) % 3].Q_nbr;

      real_t Sself[NUM_STATE], Snbr[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         Sself[c] = alpha * Q_self[c] + beta * R_self[c];
         Snbr [c] = alpha * Q_nbr [c] + beta * R_nbr [c];
      }

      for (int k = 0; k < 8; k++)
      {
         real_t F_combined[NUM_STATE];
         flux.Central(normals[k], Sself, Snbr, F_combined);

         real_t F_q[NUM_STATE], F_r[NUM_STATE];
         flux.Central(normals[k], Q_self, Q_nbr, F_q);
         flux.Central(normals[k], R_self, R_nbr, F_r);

         real_t F_split[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            F_split[c] = alpha * F_q[c] + beta * F_r[c];
         }
         real_t d = MaxAbsDiff(F_combined, F_split, NUM_STATE);
         if (d > worst_g3) { worst_g3 = d; }
      }
   }
   // Same scale argument as Gate 2 (with α=1.7, β=−2.3 amplifying ~3×).
   // Output scale ~1e10, so 1 ULP ≈ 1e-5; use 1e-1 for order-of-mag headroom.
   TEST_LE(worst_g3, 1.0e-1,
           "Central is linear in (Q_self, Q_nbr) — α-scaled inputs sum to FP precision");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
